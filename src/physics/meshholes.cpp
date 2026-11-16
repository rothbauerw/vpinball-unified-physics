// license:GPLv3+

#include "core/stdafx.h"
#include "physics/meshholes.h"
#include "physics/contactsurfacehost.h"
#include "physics/contactsystem.h"
#include "core/player.h"
#include "ui/live/LiveUI.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <tuple>

// The cut (paper section 3 "Build"): the host's top faces within the outline's box plus 54 mm are split by a 54 mm plan grid; a grid cell the outline overlaps is split by the
// rays from the outline's centre through its vertices, and each wedge keeps what lies outside its outline edge (the cell minus the outline, as convex pieces whose shared
// edges match). Every new point takes its face's own plane height. No T-junction: every face (inside the region or not) with a new point on one of its edges is split at it
// (a fan from its centroid). A through hole gets one wall quad per rim edge, down by the host's thickness_mm. The host's face list is replaced, welded and linked again.
namespace
{
   constexpr double U = ContactSystem::kMetresPerUnit;
   constexpr double kPi = 3.14159265358979323846;
   constexpr double kGridM = 54e-3;     // the plan grid round a hole and the margin past its outline's box: 4 ball radii (paper section 3, INFERRED choice)
   constexpr double kExcessM = 0.2e-3;  // a round outline's vertex stands at most this far outside its circle (paper section 4, Q1 decided)
   constexpr double kSnapM = 1e-9;      // two points of the cut closer than this are one (the clip's double rounding; far below the stored float32's 0.03 um near the table)
   constexpr double kOnEdgeM = 1e-7;    // a new point this close to a face's edge lies on it (no T-junction)
   constexpr double kTopTolM = 1e-3;    // the host's top: the nearest face plane at or below the part's height plus this (stress v2's playfield_mesh is 0.05 VPU above the parts)
   constexpr double kUpNz = 0.5;        // a face facing up: |n z| over this (within 60 deg of horizontal; the stored winding is either way, rubbercontact.h SurfaceTri)

   struct P2 { double x, y; };
   struct P3 { double x, y, z; };
   struct T3 { P3 a, b, c; };
   struct Outline { std::string name; std::vector<P2> v; bool through; double z; };

   bool g_on = false;
   std::vector<Outline> g_outlines;

   void Notice(const std::string &msg)   // rule 8: on screen and in the log
   {
      PLOGW << msg;
      if (g_pplayer != nullptr && g_pplayer->m_liveUI != nullptr)
         g_pplayer->m_liveUI->PushNotification(msg, 8000);
   }

   double Cross(const P2 &a, const P2 &b, const P2 &c) { return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x); }
   double Area(const std::vector<P2> &p)   // signed, counter-clockwise positive
   {
      double s = 0.0;
      for (size_t i = 0; i < p.size(); i++)
         s += p[i].x * p[(i + 1) % p.size()].y - p[(i + 1) % p.size()].x * p[i].y;
      return 0.5 * s;
   }
   std::vector<P2> ClipLeft(const std::vector<P2> &p, const P2 &a, const P2 &b)   // what lies on or left of the directed line a -> b (convex p)
   {
      std::vector<P2> o;
      const size_t n = p.size();
      for (size_t i = 0; i < n; i++)
      {
         const P2 &s = p[i], &e = p[(i + 1) % n];
         const double fs = Cross(a, b, s), fe = Cross(a, b, e);
         if (fs >= 0.0)
            o.push_back(s);
         if ((fs >= 0.0) != (fe >= 0.0))
         {
            const double t = fs / (fs - fe);
            o.push_back({ s.x + t * (e.x - s.x), s.y + t * (e.y - s.y) });
         }
      }
      return o;
   }
   std::vector<P2> ClipConvex(std::vector<P2> p, const std::vector<P2> &h)   // p inside the counter-clockwise convex h
   {
      for (size_t i = 0; i < h.size() && p.size() >= 3; i++)
         p = ClipLeft(p, h[i], h[(i + 1) % h.size()]);
      return p;
   }
   std::vector<P2> Box(std::vector<P2> p, double x0, double x1, double y0, double y1)
   {
      p = ClipLeft(p, { x0, 1.0 }, { x0, 0.0 });   // x >= x0
      p = ClipLeft(p, { x1, 0.0 }, { x1, 1.0 });   // x <= x1
      p = ClipLeft(p, { 0.0, y0 }, { 1.0, y0 });   // y >= y0
      p = ClipLeft(p, { 1.0, y1 }, { 0.0, y1 });   // y <= y1
      return p;
   }
   P3 Normal(const T3 &t)
   {
      const double ux = t.b.x - t.a.x, uy = t.b.y - t.a.y, uz = t.b.z - t.a.z, vx = t.c.x - t.a.x, vy = t.c.y - t.a.y, vz = t.c.z - t.a.z;
      return { uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx };
   }
   bool Up(const T3 &t)
   {
      const P3 n = Normal(t);
      const double l = std::sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
      return l > 0.0 && std::abs(n.z) > kUpNz * l;
   }
   double PlaneZ(const T3 &t, double x, double y)   // (an up face only)
   {
      const P3 n = Normal(t);
      return t.a.z - (n.x * (x - t.a.x) + n.y * (y - t.a.y)) / n.z;
   }
   std::vector<P2> Plan(const T3 &t) { return { { t.a.x, t.a.y }, { t.b.x, t.b.y }, { t.c.x, t.c.y } }; }
   bool InTri(const std::vector<P2> &p, const P2 &q)   // on or inside, either winding
   {
      const double d0 = Cross(p[0], p[1], q), d1 = Cross(p[1], p[2], q), d2 = Cross(p[2], p[0], q);
      return (d0 >= 0.0 && d1 >= 0.0 && d2 >= 0.0) || (d0 <= 0.0 && d1 <= 0.0 && d2 <= 0.0);
   }
   bool SegCrosses(const P2 &s, const P2 &e, const std::vector<P2> &h)   // the segment s e runs through the convex h's interior
   {
      double t0 = 0.0, t1 = 1.0;
      for (size_t i = 0; i < h.size(); i++)
      {
         const double fs = Cross(h[i], h[(i + 1) % h.size()], s), fe = Cross(h[i], h[(i + 1) % h.size()], e);
         if (fs <= 0.0 && fe <= 0.0)
            return false;
         if (fs < 0.0 || fe < 0.0)
         {
            const double t = fs / (fs - fe);
            if (fs < 0.0)
               t0 = std::max(t0, t);
            else
               t1 = std::min(t1, t);
         }
      }
      return t1 - t0 > 1e-9;
   }
   bool Same(const P3 &a, const P3 &b) { return a.x == b.x && a.y == b.y && a.z == b.z; }

   struct Cut
   {
      std::vector<T3> f;      // the host's faces, m
      std::vector<P3> pool;   // the points of the current cut (snapped)
      P3 Snap(const P3 &p)
      {
         for (const P3 &q : pool)
            if (std::abs(q.x - p.x) <= kSnapM && std::abs(q.y - p.y) <= kSnapM && std::abs(q.z - p.z) <= kSnapM)
               return q;
         pool.push_back(p);
         return p;
      }
      void Emit(std::vector<T3> &out, const std::vector<P2> &poly, const T3 &plane, bool flip)   // a convex counter-clockwise plan polygon on the face's plane, as a fan
      {
         std::vector<P3> q;
         for (const P2 &p : poly)
         {
            const P3 s = Snap({ p.x, p.y, PlaneZ(plane, p.x, p.y) });
            if (q.empty() || !Same(q.back(), s))
               q.push_back(s);
         }
         while (q.size() > 1 && Same(q.front(), q.back()))
            q.pop_back();
         for (bool again = true; again && q.size() >= 3;)   // points on a straight side go (a neighbour that has them splits this face's edge: T-pass)
         {
            again = false;
            for (size_t i = 0; i < q.size(); i++)
            {
               const P3 &a = q[(i + q.size() - 1) % q.size()], &b = q[i], &c = q[(i + 1) % q.size()];
               const double l1 = std::hypot(b.x - a.x, b.y - a.y), l2 = std::hypot(c.x - b.x, c.y - b.y);
               if (std::abs(Cross({ a.x, a.y }, { b.x, b.y }, { c.x, c.y })) <= 1e-10 * l1 * l2)
               {
                  q.erase(q.begin() + i);
                  again = true;
                  break;
               }
            }
         }
         if (q.size() < 3)
            return;
         for (size_t k = 1; k + 1 < q.size(); k++)
            out.push_back(flip ? T3 { q[0], q[k + 1], q[k] } : T3 { q[0], q[k], q[k + 1] });
      }
   };

   // T-pass: every face with a new point inside one of its edges is split there (a fan from its centroid: no zero-area face, the face's own plane)
   void NoTJunction(std::vector<T3> &f, const std::vector<P3> &fresh)
   {
      std::vector<T3> out;
      out.reserve(f.size() + 64);
      for (const T3 &t : f)
      {
         const P3 v[3] = { t.a, t.b, t.c };
         double lo[3], hi[3];
         for (int k = 0; k < 3; k++)
         {
            lo[k] = std::min({ (&t.a.x)[k], (&t.b.x)[k], (&t.c.x)[k] }) - kOnEdgeM;
            hi[k] = std::max({ (&t.a.x)[k], (&t.b.x)[k], (&t.c.x)[k] }) + kOnEdgeM;
         }
         std::vector<P3> poly;
         bool split = false;
         for (int e = 0; e < 3; e++)
         {
            const P3 &p = v[e], &q = v[(e + 1) % 3];
            poly.push_back(p);
            const double dx = q.x - p.x, dy = q.y - p.y, dz = q.z - p.z, l2 = dx * dx + dy * dy + dz * dz;
            if (l2 <= 0.0)
               continue;
            std::vector<std::pair<double, P3>> on;
            for (const P3 &s : fresh)
            {
               if (s.x < lo[0] || s.x > hi[0] || s.y < lo[1] || s.y > hi[1] || s.z < lo[2] || s.z > hi[2] || Same(s, p) || Same(s, q))
                  continue;
               const double u = ((s.x - p.x) * dx + (s.y - p.y) * dy + (s.z - p.z) * dz) / l2;
               if (u <= 0.0 || u >= 1.0)
                  continue;
               const double ex = p.x + u * dx - s.x, ey = p.y + u * dy - s.y, ez = p.z + u * dz - s.z;
               if (ex * ex + ey * ey + ez * ez <= kOnEdgeM * kOnEdgeM)
                  on.push_back({ u, s });
            }
            std::sort(on.begin(), on.end(), [](const auto &a, const auto &b) { return a.first < b.first; });
            for (const auto &o : on)
               poly.push_back(o.second);
            split |= !on.empty();
         }
         if (!split)
         {
            out.push_back(t);
            continue;
         }
         const P3 g { (t.a.x + t.b.x + t.c.x) / 3.0, (t.a.y + t.b.y + t.c.y) / 3.0, (t.a.z + t.b.z + t.c.z) / 3.0 };
         for (size_t k = 0; k < poly.size(); k++)
            out.push_back({ g, poly[k], poly[(k + 1) % poly.size()] });
      }
      f.swap(out);
   }

   bool IsMeshKind(const std::string &k) { return k == "playfield" || k == "wall" || k == "ramp" || k == "primitive"; }

   // one outline into the host's faces; "" or the reason it is refused (nothing changed then)
   std::string CutOne(Cut &c, const Outline &o, double topZ, double thickM, int &nWalls, double &holeArea, std::vector<std::pair<P3, P3>> &bottoms)
   {
      const std::vector<P2> &H = o.v;
      const size_t n = H.size();
      P2 ctr { 0.0, 0.0 };
      double hx0 = H[0].x, hx1 = H[0].x, hy0 = H[0].y, hy1 = H[0].y;
      for (const P2 &p : H)
      {
         ctr.x += p.x / (double)n;
         ctr.y += p.y / (double)n;
         hx0 = std::min(hx0, p.x); hx1 = std::max(hx1, p.x); hy0 = std::min(hy0, p.y); hy1 = std::max(hy1, p.y);
      }
      const double rx0 = hx0 - kGridM, rx1 = hx1 + kGridM, ry0 = hy0 - kGridM, ry1 = hy1 + kGridM;
      const auto isTop = [&](const T3 &t) { return Up(t) && std::abs(PlaneZ(t, ctr.x, ctr.y) - topZ) <= kTopTolM; };
      const auto overlaps = [](const T3 &t, double x0, double x1, double y0, double y1) {
         return std::max({ t.a.x, t.b.x, t.c.x }) >= x0 && std::min({ t.a.x, t.b.x, t.c.x }) <= x1 && std::max({ t.a.y, t.b.y, t.c.y }) >= y0 && std::min({ t.a.y, t.b.y, t.c.y }) <= y1;
      };
      // refusals (paper section 6): the outline must lie wholly on the host's top, with no other face of the host in the wall's depth crossing it
      const double aH = Area(H);
      double cover = 0.0;
      for (const T3 &t : c.f)
      {
         if (!overlaps(t, hx0, hx1, hy0, hy1))
            continue;
         std::vector<P2> p = Plan(t);
         if (isTop(t))
         {
            if (Area(p) < 0.0)
               std::reverse(p.begin(), p.end());
            cover += std::abs(Area(ClipConvex(p, H)));
            continue;
         }
         const double zlo = std::min({ t.a.z, t.b.z, t.c.z }), zhi = std::max({ t.a.z, t.b.z, t.c.z });
         if (zhi < topZ - thickM + 0.5e-3 || zlo > topZ + kTopTolM)
            continue;
         bool inside = InTri(p, ctr);
         for (int k = 0; k < 3 && !inside; k++)
         {
            bool in = true;
            for (size_t i = 0; i < n; i++)
               in &= Cross(H[i], H[(i + 1) % n], p[k]) > 1e-9 * 1e-3;
            inside = in || SegCrosses(p[k], p[(k + 1) % 3], H);
         }
         if (inside)
            return "it crosses a face of the host that is not its top (a wall, a step, or a face within the wall's depth), or the mesh already has a hole there";
      }
      if (std::abs(cover - aH) > 1e-6 * aH)
         return "it runs past the edge of the host's top (top under it " + std::to_string(cover * 1e6) + " of " + std::to_string(aH * 1e6) + " mm2)";
      // the cut
      std::vector<double> xs, ys;
      for (double x = rx0; x < rx1 - 1e-9; x += kGridM)
         xs.push_back(x);
      xs.push_back(rx1);
      for (double y = ry0; y < ry1 - 1e-9; y += kGridM)
         ys.push_back(y);
      ys.push_back(ry1);
      std::vector<T3> out;
      std::vector<P3> rimZ(n);
      std::vector<bool> rimSet(n, false);
      for (const T3 &t : c.f)
      {
         if (!isTop(t) || !overlaps(t, rx0, rx1, ry0, ry1))
         {
            out.push_back(t);
            continue;
         }
         std::vector<P2> p = Plan(t);
         const bool flip = Area(p) < 0.0;
         if (flip)
            std::reverse(p.begin(), p.end());
         for (size_t i = 0; i < n; i++)   // each rim point's height from the top face it lies on
            if (!rimSet[i] && InTri(p, H[i]))
            {
               rimZ[i] = { H[i].x, H[i].y, PlaneZ(t, H[i].x, H[i].y) };
               rimSet[i] = true;
            }
         // outside the region: four convex pieces
         c.Emit(out, ClipLeft(p, { rx0, 0.0 }, { rx0, 1.0 }), t, flip);   // x <= rx0
         c.Emit(out, ClipLeft(p, { rx1, 1.0 }, { rx1, 0.0 }), t, flip);   // x >= rx1
         const std::vector<P2> mid = ClipLeft(ClipLeft(p, { rx0, 1.0 }, { rx0, 0.0 }), { rx1, 0.0 }, { rx1, 1.0 });
         c.Emit(out, ClipLeft(mid, { 1.0, ry0 }, { 0.0, ry0 }), t, flip);   // y <= ry0
         c.Emit(out, ClipLeft(mid, { 0.0, ry1 }, { 1.0, ry1 }), t, flip);   // y >= ry1
         // inside: the grid's cells, minus the outline
         for (size_t i = 0; i + 1 < xs.size(); i++)
            for (size_t j = 0; j + 1 < ys.size(); j++)
            {
               const std::vector<P2> cell = Box(p, xs[i], xs[i + 1], ys[j], ys[j + 1]);
               if (cell.size() < 3)
                  continue;
               if (std::abs(Area(ClipConvex(cell, H))) <= 1e-14)
               {
                  c.Emit(out, cell, t, flip);
                  continue;
               }
               for (size_t k = 0; k < n; k++)   // the wedge between the rays through vertices k and k+1, outside edge k
               {
                  const P2 &a = H[k], &b = H[(k + 1) % n];
                  std::vector<P2> w = ClipLeft(cell, ctr, a);
                  w = ClipLeft(w, b, ctr);
                  w = ClipLeft(w, b, a);
                  if (w.size() >= 3)
                     c.Emit(out, w, t, flip);
               }
            }
      }
      for (size_t i = 0; i < n; i++)
         if (!rimSet[i])
            return "a rim point found no top face under it";
      nWalls = 0;
      if (o.through)
         for (size_t i = 0; i < n; i++)   // the wall: one quad per rim edge, down by the host's thickness (normals into the hole for the counter-clockwise outline)
         {
            const P3 A = c.Snap(rimZ[i]), B = c.Snap(rimZ[(i + 1) % n]);
            const P3 Ab = c.Snap({ A.x, A.y, A.z - thickM }), Bb = c.Snap({ B.x, B.y, B.z - thickM });
            out.push_back({ A, B, Ab });
            out.push_back({ B, Bb, Ab });
            bottoms.push_back({ Ab, Bb });   // the wall's bottom edge: open by design (load check)
            nWalls++;
         }
      // no T-junction: the points this cut made, in every face of the host (the original vertices are not new)
      std::set<std::tuple<double, double, double>> orig;
      for (const T3 &t : c.f)
         for (const P3 *v : { &t.a, &t.b, &t.c })
            orig.insert({ v->x, v->y, v->z });
      std::vector<P3> fresh;
      for (const P3 &q : c.pool)
         if (!orig.count({ q.x, q.y, q.z }))
            fresh.push_back(q);
      NoTJunction(out, fresh);
      c.f.swap(out);
      c.pool.clear();
      std::vector<P2> hr(n);   // the hole as stored (float32 table units)
      for (size_t i = 0; i < n; i++)
         hr[i] = { (double)(float)(rimZ[i].x / U) * U, (double)(float)(rimZ[i].y / U) * U };
      holeArea = Area(hr);
      return std::string();
   }

   struct HostCheck { int open = 0; double top = 0.0; std::vector<std::pair<P3, P3>> edges; };   // edges: the open ones
   P3 Vert(const RubberContact::SurfaceTri &t, int k) { const RubberContact::V3 &v = k == 0 ? t.a : k == 1 ? t.b : t.c; return { v.x, v.y, v.z }; }
   bool OnSeg(const P3 &p, const P3 &a, const P3 &b)   // p on the segment a b (within the stored float32's rounding near the table)
   {
      constexpr double tol = 1e-6;
      const double dx = b.x - a.x, dy = b.y - a.y, dz = b.z - a.z, l2 = dx * dx + dy * dy + dz * dz;
      const double u = l2 > 0.0 ? std::clamp(((p.x - a.x) * dx + (p.y - a.y) * dy + (p.z - a.z) * dz) / l2, 0.0, 1.0) : 0.0;
      const double ex = a.x + u * dx - p.x, ey = a.y + u * dy - p.y, ez = a.z + u * dz - p.z;
      return ex * ex + ey * ey + ez * ez <= tol * tol;
   }
   HostCheck Measure(const RubberSurface &s)
   {
      HostCheck h;
      for (size_t t = 0; t < s.nb.size(); t++)
         for (int k = 0; k < 3; k++)
            if (s.nb[t][k] < 0)
            {
               h.open++;
               h.edges.push_back({ Vert(s.tris[t], k), Vert(s.tris[t], (k + 1) % 3) });
            }
      for (const RubberContact::SurfaceTri &t : s.tris)
      {
         const T3 q { { t.a.x, t.a.y, t.a.z }, { t.b.x, t.b.y, t.b.z }, { t.c.x, t.c.y, t.c.z } };
         if (Up(q))
            h.top += std::abs(Area(Plan(q)));
      }
      return h;
   }
}

void MeshHoles::Clear(bool on)
{
   g_on = on;
   g_outlines.clear();   // (the kickers add theirs in PhysicSetup, the turntables in Turntables::Setup, both after this)
}

bool MeshHoles::On()
{
   return g_on;
}

void MeshHoles::Add(const std::string &name, const std::vector<double> &xy, bool through, double partZ)
{
   Outline o { name, {}, through, partZ };
   for (size_t i = 0; i + 1 < xy.size(); i += 2)
      o.v.push_back({ xy[i], xy[i + 1] });
   if (o.v.size() >= 3 && Area(o.v) < 0.0)
      std::reverse(o.v.begin(), o.v.end());
   g_outlines.push_back(std::move(o));
}

std::vector<double> MeshHoles::Round(double x, double y, double r)
{
   int n = 24;
   while (r * (1.0 / std::cos(kPi / n) - 1.0) > kExcessM)
      n++;
   const double R = r / std::cos(kPi / n);   // vertices outside the circle: every edge touches it at its midpoint (paper section 4)
   std::vector<double> xy;
   for (int k = 0; k < n; k++)
   {
      xy.push_back(x + R * std::cos(2.0 * kPi * k / n));
      xy.push_back(y + R * std::sin(2.0 * kPi * k / n));
   }
   return xy;
}

void MeshHoles::Build()
{
   if (!g_on || g_outlines.empty())
      return;
   // the host under each outline: the static mesh host whose up face under the outline's centre is the nearest at or below the part's height
   std::map<ContactSurfaceHost *, std::vector<std::pair<const Outline *, double>>> byHost;
   for (const Outline &o : g_outlines)
   {
      if (o.v.size() < 3)
         continue;
      P2 ctr { 0.0, 0.0 };
      for (const P2 &p : o.v)
         ctr = { ctr.x + p.x / (double)o.v.size(), ctr.y + p.y / (double)o.v.size() };
      ContactSurfaceHost *best = nullptr;
      double bestZ = 0.0;
      for (ContactHost *const h : ContactHost::AllHosts())
      {
         ContactSurfaceHost *const s = dynamic_cast<ContactSurfaceHost *>(h);
         if (s == nullptr || s->MeshSurface() == nullptr || !IsMeshKind(s->MeshKind()))
            continue;
         for (const RubberContact::SurfaceTri &t : s->MeshSurface()->tris)
         {
            const T3 q { { t.a.x, t.a.y, t.a.z }, { t.b.x, t.b.y, t.b.z }, { t.c.x, t.c.y, t.c.z } };
            if (!Up(q) || !InTri(Plan(q), ctr))
               continue;
            const double z = PlaneZ(q, ctr.x, ctr.y);
            if (z <= o.z + kTopTolM && (best == nullptr || z > bestZ))
            {
               best = s;
               bestZ = z;
            }
         }
      }
      if (best == nullptr)
      {
         Notice("Mesh holes: " + o.name + ": no static mesh top found under it at or below its height: hole not cut");
         continue;
      }
      byHost[best].push_back({ &o, bestZ });
   }
   for (auto &[host, list] : byHost)
   {
      RubberSurface &s = *host->MeshSurface();
      const HostCheck before = Measure(s);
      const size_t facesBefore = s.tris.size();
      Cut c;
      for (const RubberContact::SurfaceTri &t : s.tris)
         c.f.push_back({ { t.a.x, t.a.y, t.a.z }, { t.b.x, t.b.y, t.b.z }, { t.c.x, t.c.y, t.c.z } });
      const double thickM = host->ThicknessMm() * 1e-3;
      int cut = 0;
      double expectTop = before.top;
      std::vector<std::pair<P3, P3>> bottoms;   // the walls' bottom edges (open by design)
      std::ostringstream per;
      std::vector<const Outline *> done;
      for (const auto &[o, topZ] : list)
      {
         std::string why;
         for (const Outline *d : done)   // overlapping outlines are joined in a later piece: refused here
            if (std::abs(Area(ClipConvex(o->v, d->v))) > 0.0)
               why = "it overlaps " + d->name + " (joining outlines is a later piece)";
         int nWalls = 0;
         double aHole = 0.0;
         if (why.empty())
            why = CutOne(c, *o, topZ, thickM, nWalls, aHole, bottoms);
         if (!why.empty())
         {
            Notice("Mesh holes: " + o->name + " not cut in the " + host->MeshKind() + " host: " + why);
            continue;
         }
         done.push_back(o);
         cut++;
         expectTop -= aHole;
         per << " " << o->name << " N " << o->v.size() << " (top z " << topZ * 1e3 << " mm);";
      }
      if (cut == 0)
         continue;
      RubberSurface ns;   // the cut faces, welded by their float bits as every host's (contactsurfacehost.cpp RubberSurface::Add); the uncut list is not kept
      for (const T3 &t : c.f)
         ns.Add(Vertex3Ds((float)(t.a.x / U), (float)(t.a.y / U), (float)(t.a.z / U)), Vertex3Ds((float)(t.b.x / U), (float)(t.b.y / U), (float)(t.b.z / U)),
                Vertex3Ds((float)(t.c.x / U), (float)(t.c.y / U), (float)(t.c.z / U)));
      s.tris = std::move(ns.tris);
      std::copy(ns.lo, ns.lo + 3, s.lo);
      std::copy(ns.hi, ns.hi + 3, s.hi);
      s.vid = std::move(ns.vid);
      s.weld = std::move(ns.weld);
      s.nb.clear();
      s.fold.clear();
      host->MeshRelink();   // shared edges again (contactsurfacehost.cpp Link) and the host's box
      const HostCheck after = Measure(s);
      // every open edge after the cut is a wall's bottom or a piece of an open edge the mesh had (its border, split where a grid line crosses it); any other is a gap the cut made
      int stray = 0;
      for (const auto &[a, b] : after.edges)
      {
         bool ok = false;
         for (const auto &[p, q] : bottoms)
            ok = ok || (OnSeg(a, p, q) && OnSeg(b, p, q));
         for (const auto &[p, q] : before.edges)
            ok = ok || (OnSeg(a, p, q) && OnSeg(b, p, q));
         stray += !ok;
      }
      PLOGI << "Mesh holes: " << host->MeshKind() << " host: faces " << facesBefore << " -> " << s.tris.size() << ", holes " << cut << ":" << per.str() << " open edges " << before.open << " -> " << after.open
            << " (" << bottoms.size() << " wall bottoms, " << stray << " stray), top area " << before.top << " -> " << after.top << " m2 (expected " << expectTop << ", diff " << (after.top - expectTop) << ")";
      if (stray > 0 || std::abs(after.top - expectTop) > 1e-9 * before.top)
         Notice("Mesh holes: load check failed on the " + host->MeshKind() + " host: " + std::to_string(stray) + " open edges that are neither a wall's bottom nor the mesh's own border, top area off by "
                + std::to_string((after.top - expectTop) * 1e6) + " mm2");
   }
}
