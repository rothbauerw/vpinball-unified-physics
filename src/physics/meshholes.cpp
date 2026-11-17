// license:GPLv3+

#include "core/stdafx.h"
#include "physics/meshholes.h"
#include "physics/contactsurfacehost.h"
#include "physics/contactsystem.h"
#include "physics/physicsini.h"
#include "core/player.h"
#include "ui/live/LiveUI.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <tuple>

// The cut (paper section 3 "Build"): outlines on one host whose plans overlap are one group, cut together (paper Q5: joined, no warning); an outline wholly inside another's floor
// (a slot in the groove) is cut after it, into that floor's faces. For a group, the host's top faces within the group's box plus 54 mm are split by a 54 mm plan grid; each grid cell
// minus the group's outlines is found exactly for any simple outline (every cell and outline edge split where they cross, a piece kept where it bounds what is left, the pieces chained
// into loops) and each loop, with any outline wholly inside it joined by a bridge, is ear-clipped. Every new point takes its face's own plane height. No T-junction: every face with a
// new point on one of its edges is split at it (a fan from its centroid). Walls (one quad per rim piece, down by the host's thickness_mm) only where the ball can pass through: the
// joined outline's largest inscribed circle is at least the ball's radius (paper section 3 step 4). Floors (the slot's ramps, the groove's strips) are faces under the top, given as
// plan triangles with a depth below the surface. The host's face list is replaced, welded and linked again.
namespace
{
   constexpr double U = ContactSystem::kMetresPerUnit;
   constexpr double kPi = 3.14159265358979323846;
   constexpr double kGridM = 54e-3;     // the plan grid round a hole and the margin past its outline's box: 4 ball radii (paper section 3, INFERRED choice)
   constexpr double kGridClearM = 1e-6; // no grid line passes this near an outline vertex, and no outline edge this near a grid node (the grid's origin is moved till none does; lead 2026-11-16)
   double g_excessM = 0.2e-3;           // a round outline's vertex stands at most this far outside its circle ([Features] hole_round_mm, physicsini.h; set by Clear each load)
   // the true curves of the outlines' round parts, for the drawing only (brief 6, user 2026-11-17: the drawn rim on the true curve, not the physics' polygon): each circle Round or a
   // rounded corner made (R > r: vertices at R, every edge tangent to the circle of radius r) and each slot end (R = r: vertices on the circle, chords); cleared by Clear
   struct Arc { double cx, cy, r, R; };
   std::vector<Arc> g_arcs;
   constexpr double kDrawTolM = 0.02e-3;   // the drawn curve's chord stands at most this far inside the true curve (playfieldcutdraw's kDrawTolM)
   constexpr double kSnapM = 0.1e-6;    // two points of the cut closer than this (3D) are one: the clip's double rounding and crossings a few float32 steps apart (the stored float32 is 0.03 um near
                                        // the table: a 68 nm wall face was degenerate in it, p3 run T2); an existing vertex of the region wins (Cut::pool is seeded with them)
   constexpr double kOnEdgeM = 1e-7;    // a new point this close to a face's edge lies on it (no T-junction)
   constexpr double kTopTolM = 1e-3;    // the host's top: the nearest face plane at or below the part's height plus this (stress v2's playfield_mesh is 0.05 VPU above the parts)
   constexpr double kIsTopTolM = 0.25e-3;   // a face of the cut's top: an up face whose plane at the outline's centre is within this of the top (under a 1 VPU slab's 0.54 mm; lead 2026-11-16)
   constexpr double kUpNz = 0.5;        // a face facing up: |n z| over this (within 60 deg of horizontal; the stored winding is either way, rubbercontact.h SurfaceTri)
   constexpr double kBallRM = 25.0 * U; // the ball the walls rule is judged for: the standard ball (the table carries no ball size)
   constexpr double kTurnRad = kPi / 12.0;   // a rounded corner's tangent polygon turns at most 15 deg per vertex (paper section 4: at least 6 per quarter, more by hole_round_mm)
   // [PlayfieldCut:<name>] shapes (slot-hole design 1, 4, 9; paper section 4 for the groove)
   constexpr double kRampM = 7.5e-3;          // the slot's drilled ramp: the outline narrows to a point this far from the end of the straight sides (user 2026-11-09; fixed)
   constexpr double kRampDepthM = 1.75e-3;    // the ramp floor's depth on the axis at the round end's far edge (u = a), rising linearly to 0 at the tip (slot-hole design 2: drill cross-section)
   constexpr double kWidthMaxMm = 5.03;       // slot width_mm default and largest: 0.198 in (PartData 7)
   constexpr double kWidthMinMm = 3.175;      // slot width_mm smallest: 1/8 in
   constexpr double kGrooveRM = 0.625 * 25.4e-3;      // the groove's half-round radius (PartData 7, 1.25 in diameter; MEASURED part data)
   constexpr double kGrooveDepthM = 0.188 * 25.4e-3;  // its depth (PartData 7, 0.188 in)
   constexpr int kGrooveStrips = 16;                  // straight strips tangent to the half-round (paper section 4: N 16, the smallest crease at the cost of 12)
   constexpr double kGrooveLenMinMm = 12.0 * 25.4;    // groove length_mm smallest and largest: 12 - 20 in (PartData 7)
   constexpr double kGrooveLenMaxMm = 20.0 * 25.4;

   struct P2 { double x, y; };
   struct P3 { double x, y, z; };
   struct T3 { P3 a, b, c; };
   struct FTri { P2 p[3]; double d[3]; };   // a face in plan with its depth below the surface at each corner (m); a vertical face has two corners at one plan point
   struct Outline
   {
      std::string name;
      std::vector<P2> v;                 // counter-clockwise
      double z = 0.0;                    // the part's height: the host is the nearest top at or below it
      std::string host;                  // the host named in the ini (empty: found at the rim points)
      MeshHoles::HoleKeys keys;          // draw, wall_depth_mm (cut_hole 0: never added)
      double passR = -1.0;               // >= 0: the largest circle through the hole's through part when smaller than the outline's (a chamfer's inner hole); the walls rule takes it
      std::vector<char> floored;         // per edge (v[i] -> v[i + 1]): a floor face meets the top along it (no wall, not an open rim)
      std::vector<FTri> floor, side;     // floor faces (up) and vertical faces (the groove's square end), depth below the surface
      std::vector<std::array<std::pair<P2, double>, 2>> openFloor;   // floor edges open by design (the ramp floor's edge over the slot's through part)
   };

   std::vector<Outline> g_outlines;
   std::vector<std::function<void()>> g_place;   // the parts' deferred placing steps (MeshHoles::Defer; run and emptied by PlaceParts, emptied by Clear)
   std::vector<MeshHoles::DrawnHole> g_drawn;   // every joined hole as Build cut it (the drawing's record)

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
   bool Eq(const P2 &a, const P2 &b) { return a.x == b.x && a.y == b.y; }
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
         if ((fs > 0.0 && fe < 0.0) || (fs < 0.0 && fe > 0.0))   // a crossing only where the sides strictly differ (Sutherland-Hodgman): a point on the line is kept once, as itself, never again as a t = 0 or t = 1 crossing
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
   std::vector<P2> PlanCCW(const T3 &t)
   {
      std::vector<P2> p = Plan(t);
      if (Area(p) < 0.0)
         std::reverse(p.begin(), p.end());
      return p;
   }
   bool InTri(const std::vector<P2> &p, const P2 &q)   // on or inside, either winding
   {
      const double d0 = Cross(p[0], p[1], q), d1 = Cross(p[1], p[2], q), d2 = Cross(p[2], p[0], q);
      return (d0 >= 0.0 && d1 >= 0.0 && d2 >= 0.0) || (d0 <= 0.0 && d1 <= 0.0 && d2 <= 0.0);
   }
   bool InPoly(const std::vector<P2> &h, const P2 &q)   // even-odd (any simple polygon)
   {
      bool in = false;
      for (size_t i = 0, j = h.size() - 1; i < h.size(); j = i++)
         if ((h[i].y > q.y) != (h[j].y > q.y) && q.x < (h[j].x - h[i].x) * (q.y - h[i].y) / (h[j].y - h[i].y) + h[i].x)
            in = !in;
      return in;
   }
   bool Same(const P3 &a, const P3 &b) { return a.x == b.x && a.y == b.y && a.z == b.z; }
   // the segments a1 a2 and b1 b2 cross at one point inside both (t, u in (0, 1)); the point from a's parameter, so a call with the same arguments in the same order gives the same bits
   bool Isect(const P2 &a1, const P2 &a2, const P2 &b1, const P2 &b2, P2 &x)
   {
      const double rx = a2.x - a1.x, ry = a2.y - a1.y, sx = b2.x - b1.x, sy = b2.y - b1.y, d = rx * sy - ry * sx;
      if (d == 0.0)
         return false;
      const double t = ((b1.x - a1.x) * sy - (b1.y - a1.y) * sx) / d, u = ((b1.x - a1.x) * ry - (b1.y - a1.y) * rx) / d;
      if (t <= 1e-12 || t >= 1.0 - 1e-12 || u <= 1e-12 || u >= 1.0 - 1e-12)
         return false;
      x = { a1.x + t * rx, a1.y + t * ry };
      return true;
   }
   bool SegX(const P2 &a, const P2 &b, const P2 &c, const P2 &d)   // proper crossing (no shared end)
   {
      const double d1 = Cross(c, d, a), d2 = Cross(c, d, b), d3 = Cross(a, b, c), d4 = Cross(a, b, d);
      return ((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0)) && ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0));
   }

   bool OnEdge(const P2 &a, const P2 &b, const P2 &q)   // q inside the segment a b (within 0.1 nm; not at its ends)
   {
      const double l = std::hypot(b.x - a.x, b.y - a.y), u = ((q.x - a.x) * (b.x - a.x) + (q.y - a.y) * (b.y - a.y)) / (l * l);
      return l > 0.0 && u > 1e-12 && u < 1.0 - 1e-12 && std::abs(Cross(a, b, q)) / l <= 1e-10;
   }

   using Group = std::vector<const Outline *>;
   bool InAny(const Group &g, const P2 &q, int skip = -1)
   {
      for (size_t k = 0; k < g.size(); k++)
         if ((int)k != skip && InPoly(g[k]->v, q))
            return true;
      return false;
   }
   double Param(const P2 &a, const P2 &b, const P2 &q) { return ((q.x - a.x) * (b.x - a.x) + (q.y - a.y) * (b.y - a.y)) / ((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y)); }
   // member m's edge e split where the other members' edges cross it (each crossing computed from the lower-numbered member's edge: the same bits on both), and where extra edges cross it
   std::vector<P2> SplitEdge(const Group &g, size_t m, size_t e, const std::vector<std::pair<P2, P2>> *extra)
   {
      const std::vector<P2> &h = g[m]->v;
      const P2 a = h[e], b = h[(e + 1) % h.size()];
      std::vector<std::pair<double, P2>> pts { { 0.0, a }, { 1.0, b } };
      P2 x;
      for (size_t k = 0; k < g.size(); k++)
      {
         if (k == m)
            continue;
         const std::vector<P2> &o = g[k]->v;
         for (size_t f = 0; f < o.size(); f++)
         {
            const P2 &c = o[f], &d = o[(f + 1) % o.size()];
            if (k > m ? Isect(a, b, c, d, x) : Isect(c, d, a, b, x))
               pts.push_back({ Param(a, b, x), x });
         }
      }
      if (extra != nullptr)
         for (const auto &[c, d] : *extra)
         {
            if (Isect(a, b, c, d, x))   // (outline edge first: the extra edge's split uses the same call)
               pts.push_back({ Param(a, b, x), x });
            if (OnEdge(a, b, c))   // a cell corner on the edge (a face's corner: the grid's own nodes are kept clear of the outlines)
               pts.push_back({ Param(a, b, c), c });
         }
      std::sort(pts.begin(), pts.end(), [](const auto &p, const auto &q) { return p.first < q.first; });
      std::vector<P2> out;
      for (const auto &p : pts)
         if (out.empty() || !Eq(out.back(), p.second))
            out.push_back(p.second);
      return out;
   }
   struct Piece { P2 a, b; size_t m, e; };   // a piece of member m's edge e, a -> b counter-clockwise along the outline
   // the joined outline's boundary (inside = false: the pieces outside every other member) or the pieces inside another member (inside = true; in: which)
   std::vector<Piece> Pieces(const Group &g, bool inside, std::vector<size_t> *in = nullptr)
   {
      std::vector<Piece> out;
      for (size_t m = 0; m < g.size(); m++)
         for (size_t e = 0; e < g[m]->v.size(); e++)
         {
            const std::vector<P2> s = SplitEdge(g, m, e, nullptr);
            for (size_t k = 0; k + 1 < s.size(); k++)
            {
               const P2 mid { 0.5 * (s[k].x + s[k + 1].x), 0.5 * (s[k].y + s[k + 1].y) };
               int which = -1;
               for (size_t j = 0; j < g.size() && which < 0; j++)
                  if (j != m && InPoly(g[j]->v, mid))
                     which = (int)j;
               if ((which >= 0) == inside)
               {
                  out.push_back({ s[k], s[k + 1], m, e });
                  if (in != nullptr)
                     in->push_back((size_t)std::max(which, 0));
               }
            }
         }
      return out;
   }

   // ear clipping of one simple counter-clockwise polygon (a bridged one has repeated points); false: no ear found
   bool EarClip(std::vector<P2> p, std::vector<std::array<P2, 3>> &tris)
   {
      while (p.size() > 3)
      {
         const size_t n = p.size();
         bool cut = false;
         for (size_t i = 0; i < n && !cut; i++)   // a vertex on a straight run: no triangle, it goes (a neighbour that has it splits this face's edge: T-pass)
         {
            const P2 &a = p[(i + n - 1) % n], &b = p[i], &c = p[(i + 1) % n];
            const double l1 = std::hypot(b.x - a.x, b.y - a.y), l2 = std::hypot(c.x - b.x, c.y - b.y);
            if (std::abs(Cross(a, b, c)) <= 1e-10 * l1 * l2 && (a.x - b.x) * (c.x - b.x) + (a.y - b.y) * (c.y - b.y) < 0.0)
            {
               p.erase(p.begin() + i);
               cut = true;
            }
         }
         if (cut)
            continue;
         for (size_t i = 0; i < n && !cut; i++)
         {
            const P2 &a = p[(i + n - 1) % n], &b = p[i], &c = p[(i + 1) % n];
            if (Cross(a, b, c) <= 0.0)
               continue;
            const std::vector<P2> t { a, b, c };
            bool empty = true;
            for (size_t j = 0; j < n && empty; j++)
               if (!Eq(p[j], a) && !Eq(p[j], b) && !Eq(p[j], c) && InTri(t, p[j]))
                  empty = false;
            if (!empty)
               continue;
            tris.push_back({ a, b, c });
            p.erase(p.begin() + i);
            cut = true;
         }
         if (!cut)
            return false;
      }
      if (p.size() == 3 && Cross(p[0], p[1], p[2]) > 0.0)
         tris.push_back({ p[0], p[1], p[2] });
      return true;
   }
   // the holes (clockwise loops) joined into the outer loop by bridges, rightmost hole first, each from its rightmost point to the nearest outer point it sees
   bool Bridge(std::vector<P2> &poly, std::vector<std::vector<P2>> holes)
   {
      std::sort(holes.begin(), holes.end(), [](const auto &a, const auto &b) {
         return std::max_element(a.begin(), a.end(), [](const P2 &p, const P2 &q) { return p.x < q.x; })->x > std::max_element(b.begin(), b.end(), [](const P2 &p, const P2 &q) { return p.x < q.x; })->x;
      });
      for (size_t hi = 0; hi < holes.size(); hi++)
      {
         const std::vector<P2> &h = holes[hi];
         const size_t k = (size_t)(std::max_element(h.begin(), h.end(), [](const P2 &p, const P2 &q) { return p.x < q.x; }) - h.begin());
         const P2 hp = h[k];
         int best = -1;
         double bestD = 1e30;
         for (size_t i = 0; i < poly.size(); i++)
         {
            const P2 v = poly[i];
            const double d = std::hypot(v.x - hp.x, v.y - hp.y);
            if (d >= bestD || d == 0.0)
               continue;
            bool ok = true;
            const auto test = [&](const std::vector<P2> &l) {
               for (size_t j = 0; j < l.size() && ok; j++)
                  if (SegX(hp, v, l[j], l[(j + 1) % l.size()]))
                     ok = false;
            };
            test(poly);
            for (size_t j = hi; j < holes.size(); j++)
               test(holes[j]);
            const P2 mid { 0.5 * (hp.x + v.x), 0.5 * (hp.y + v.y) };
            if (ok && InPoly(poly, mid))
               for (size_t j = hi; j < holes.size() && ok; j++)
                  ok = !InPoly(holes[j], mid);
            else
               ok = false;
            if (ok)
            {
               best = (int)i;
               bestD = d;
            }
         }
         if (best < 0)
            return false;
         std::vector<P2> ins;
         for (size_t j = 0; j <= h.size(); j++)
            ins.push_back(h[(k + j) % h.size()]);
         ins.push_back(poly[best]);
         poly.insert(poly.begin() + best + 1, ins.begin(), ins.end());
      }
      return true;
   }
   // the convex counter-clockwise polygon P less the group's outlines, as triangles; false: the pieces did not chain into loops or a loop could not be clipped
   bool Difference(const std::vector<P2> &P, const Group &g, std::vector<std::array<P2, 3>> &tris)
   {
      std::vector<std::pair<P2, P2>> cellEdges;
      for (size_t i = 0; i < P.size(); i++)
         cellEdges.push_back({ P[i], P[(i + 1) % P.size()] });
      std::vector<std::pair<P2, P2>> pieces;   // directed: the region on the left
      P2 x;
      for (const auto &[a, b] : cellEdges)   // the cell's edges, split where the outlines cross them (the crossing from the outline edge's call, as in SplitEdge)
      {
         std::vector<std::pair<double, P2>> pts { { 0.0, a }, { 1.0, b } };
         for (const Outline *o : g)
            for (size_t f = 0; f < o->v.size(); f++)
            {
               if (Isect(o->v[f], o->v[(f + 1) % o->v.size()], a, b, x))
                  pts.push_back({ Param(a, b, x), x });
               if (OnEdge(a, b, o->v[f]))   // an outline vertex on the cell's edge (a face's edge: a groove crease under a slot's tip)
                  pts.push_back({ Param(a, b, o->v[f]), o->v[f] });
            }
         std::sort(pts.begin(), pts.end(), [](const auto &p, const auto &q) { return p.first < q.first; });
         for (size_t k = 0; k + 1 < pts.size(); k++)
         {
            if (Eq(pts[k].second, pts[k + 1].second))
               continue;
            const P2 mid { 0.5 * (pts[k].second.x + pts[k + 1].second.x), 0.5 * (pts[k].second.y + pts[k + 1].second.y) };
            if (!InAny(g, mid))
               pieces.push_back({ pts[k].second, pts[k + 1].second });
         }
      }
      for (size_t m = 0; m < g.size(); m++)   // the outlines' edges inside the cell and outside the other outlines, reversed (the hole on the right)
         for (size_t e = 0; e < g[m]->v.size(); e++)
         {
            const std::vector<P2> s = SplitEdge(g, m, e, &cellEdges);
            for (size_t k = 0; k + 1 < s.size(); k++)
            {
               const P2 mid { 0.5 * (s[k].x + s[k + 1].x), 0.5 * (s[k].y + s[k + 1].y) };
               bool inCell = true;
               for (size_t i = 0; i < P.size() && inCell; i++)
               {
                  const P2 &pa = P[i], &pb = P[(i + 1) % P.size()];
                  if (std::hypot(pb.x - pa.x, pb.y - pa.y) >= kOnEdgeM)   // a cell edge shorter than kOnEdgeM has no direction to judge a side by: skip it
                     inCell = Cross(pa, pb, mid) > 0.0;
               }
               if (inCell && !InAny(g, mid, (int)m))
                  pieces.push_back({ s[k + 1], s[k] });
            }
         }
      std::multimap<std::pair<double, double>, size_t> from;
      for (size_t i = 0; i < pieces.size(); i++)
         from.insert({ { pieces[i].first.x, pieces[i].first.y }, i });
      std::vector<char> used(pieces.size(), 0);
      std::vector<std::vector<P2>> outer, holes;
      for (size_t i0 = 0; i0 < pieces.size(); i0++)
      {
         if (used[i0])
            continue;
         std::vector<P2> loop;
         size_t i = i0;
         while (true)
         {
            used[i] = 1;
            loop.push_back(pieces[i].first);
            const P2 end = pieces[i].second;
            if (Eq(end, pieces[i0].first))
               break;
            auto [lo, hi] = from.equal_range({ end.x, end.y });
            size_t next = pieces.size();
            for (auto it = lo; it != hi && next == pieces.size(); ++it)
               if (!used[it->second])
                  next = it->second;
            if (next == pieces.size())
               return false;
            i = next;
         }
         if (loop.size() < 3)
            continue;
         (Area(loop) > 0.0 ? outer : holes).push_back(loop);
      }
      for (std::vector<P2> &o : outer)
      {
         std::vector<std::vector<P2>> mine;
         for (const std::vector<P2> &h : holes)
            if (InPoly(o, h[0]))
               mine.push_back(h);
         if (!Bridge(o, mine) || !EarClip(o, tris))
            return false;
      }
      return true;
   }
   // the joined outline's largest inscribed circle radius (m): the deepest point inside it from its boundary pieces, on a 80 x 80 grid over its box, refined twice round the best
   double Inscribed(const Group &g, const std::vector<Piece> &bd)
   {
      double x0 = 1e30, x1 = -1e30, y0 = 1e30, y1 = -1e30;
      for (const Outline *o : g)
         for (const P2 &p : o->v)
         {
            x0 = std::min(x0, p.x); x1 = std::max(x1, p.x); y0 = std::min(y0, p.y); y1 = std::max(y1, p.y);
         }
      const auto depth = [&](const P2 &q) {
         if (!InAny(g, q))
            return -1.0;
         double d = 1e30;
         for (const Piece &s : bd)
         {
            const double u = std::clamp(Param(s.a, s.b, q), 0.0, 1.0);
            d = std::min(d, std::hypot(s.a.x + u * (s.b.x - s.a.x) - q.x, s.a.y + u * (s.b.y - s.a.y) - q.y));
         }
         return d;
      };
      double step = std::max(x1 - x0, y1 - y0) / 80.0, best = 0.0;
      P2 at { 0.5 * (x0 + x1), 0.5 * (y0 + y1) };
      for (double x = x0 + 0.5 * step; x < x1; x += step)
         for (double y = y0 + 0.5 * step; y < y1; y += step)
            if (const double d = depth({ x, y }); d > best)
            {
               best = d;
               at = { x, y };
            }
      for (int r = 0; r < 2; r++)
      {
         const P2 c = at;
         for (int i = -4; i <= 4; i++)
            for (int j = -4; j <= 4; j++)
               if (const double d = depth({ c.x + i * step / 4.0, c.y + j * step / 4.0 }); d > best)
               {
                  best = d;
                  at = { c.x + i * step / 4.0, c.y + j * step / 4.0 };
               }
         step /= 4.0;
      }
      return best;
   }

   struct Cut
   {
      std::vector<T3> f;      // the host's faces, m
      std::vector<P3> pool;   // the points of the current cut (snapped)
      P3 Snap(const P3 &p)
      {
         for (const P3 &q : pool)
            if ((q.x - p.x) * (q.x - p.x) + (q.y - p.y) * (q.y - p.y) + (q.z - p.z) * (q.z - p.z) < kSnapM * kSnapM)
               return q;
         pool.push_back(p);
         return p;
      }
      // a convex counter-clockwise plan polygon at the heights z(x, y), as a fan
      template <class Z> void EmitZ(std::vector<T3> &out, const std::vector<P2> &poly, const Z &z, bool flip)
      {
         std::vector<P3> q;
         for (const P2 &p : poly)
         {
            const P3 s = Snap({ p.x, p.y, z(p) });
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
               const double lac = std::hypot(c.x - a.x, c.y - a.y);   // b within the T-pass tolerance of the line a c: it goes (no face thinner than kOnEdgeM, lead 2026-11-16)
               if (lac <= 0.0 || std::abs(Cross({ a.x, a.y }, { b.x, b.y }, { c.x, c.y })) <= kOnEdgeM * lac)
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
      void Emit(std::vector<T3> &out, const std::vector<P2> &poly, const T3 &plane, bool flip)   // on the face's plane
      {
         EmitZ(out, poly, [&](const P2 &p) { return PlaneZ(plane, p.x, p.y); }, flip);
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
         std::vector<std::pair<double, P3>> on[3];   // each fresh point on one edge only: the nearest (lead 2026-11-16: one point on two edges of a sliver made the fan overlap)
         for (const P3 &s : fresh)
         {
            if (s.x < lo[0] || s.x > hi[0] || s.y < lo[1] || s.y > hi[1] || s.z < lo[2] || s.z > hi[2])
               continue;
            int best = -1;
            double bestD = kOnEdgeM * kOnEdgeM, bestU = 0.0;
            for (int e = 0; e < 3; e++)
            {
               const P3 &p = v[e], &q = v[(e + 1) % 3];
               const double dx = q.x - p.x, dy = q.y - p.y, dz = q.z - p.z, l2 = dx * dx + dy * dy + dz * dz;
               if (l2 <= 0.0 || Same(s, p) || Same(s, q))
                  continue;
               const double u = ((s.x - p.x) * dx + (s.y - p.y) * dy + (s.z - p.z) * dz) / l2;
               if (u <= 0.0 || u >= 1.0)
                  continue;
               const double ex = p.x + u * dx - s.x, ey = p.y + u * dy - s.y, ez = p.z + u * dz - s.z, d2 = ex * ex + ey * ey + ez * ez;
               if (d2 <= bestD)
               {
                  best = e;
                  bestD = d2;
                  bestU = u;
               }
            }
            if (best >= 0 && !(Same(s, v[0]) || Same(s, v[1]) || Same(s, v[2])))
               on[best].push_back({ bestU, s });
         }
         std::vector<P3> poly;
         bool split = false;
         for (int e = 0; e < 3; e++)
         {
            poly.push_back(v[e]);
            std::sort(on[e].begin(), on[e].end(), [](const auto &a, const auto &b) { return a.first < b.first; });
            for (const auto &o : on[e])
               poly.push_back(o.second);
            split |= !on[e].empty();
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

   // the depth of a floor face at q (affine over its plan triangle)
   double DepthIn(const FTri &t, const P2 &q)
   {
      const double A = Cross(t.p[0], t.p[1], t.p[2]);
      const double w0 = Cross(t.p[1], t.p[2], q) / A, w1 = Cross(t.p[2], t.p[0], q) / A, w2 = 1.0 - w0 - w1;
      return w0 * t.d[0] + w1 * t.d[1] + w2 * t.d[2];
   }
   // a point on an edge of the floor (within kOnEdgeM: a crossing on a slanted outline edge lies a rounding step off it, 2.5e-17 m in the switch test's groove + Round, offline brief 4)
   // takes the face nearest it
   bool FloorDepth(const Outline &o, const P2 &q, double &d)
   {
      const FTri *best = nullptr;
      double bestOut = kOnEdgeM;
      for (const FTri &t : o.floor)
      {
         if (InTri({ t.p[0], t.p[1], t.p[2] }, q))
         {
            d = DepthIn(t, q);
            return true;
         }
         const double sg = Cross(t.p[0], t.p[1], t.p[2]) > 0.0 ? 1.0 : -1.0;
         double out = 0.0;   // the point's largest distance outside one of the face's edges
         for (int k = 0; k < 3; k++)
         {
            const P2 &a = t.p[k], &b = t.p[(k + 1) % 3];
            if (const double l = std::hypot(b.x - a.x, b.y - a.y); l > 0.0)
               out = std::max(out, -sg * Cross(a, b, q) / l);
         }
         if (out <= bestOut)
         {
            bestOut = out;
            best = &t;
         }
      }
      if (best == nullptr)
         return false;
      d = DepthIn(*best, q);
      return true;
   }

   struct GroupResult { int nWalls = 0, nDrop = 0, nDeep = 0, nFloor = 0, nUnder = 0; double underCut = 0.0, holeArea = 0.0, holeAreaD = 0.0, dropUp = 0.0, missing = 0.0, floorUp = 0.0, inscribed = 0.0, gridOff = 0.0; };

   // one group of joined outlines into the host's faces; "" or the reason it is refused (nothing changed then)
   // bottoms: the walls' bottom edges; rims: open edges by design (the rim of a hole without walls, a floor's edge over a through part); both for the load check
   // dh: the hole as cut, for the drawing (meshholes.h DrawnHole: its rim loops, walls, wall and floor triangles)
   std::string CutGroup(Cut &c, const Group &g, double topZ, double thickM, GroupResult &r, std::vector<std::pair<P3, P3>> &bottoms, std::vector<std::pair<P3, P3>> &rims,
                        MeshHoles::DrawnHole &dh)
   {
      r = GroupResult();
      P2 ctr { 0.0, 0.0 };
      double hx0 = 1e30, hx1 = -1e30, hy0 = 1e30, hy1 = -1e30;
      size_t nv = 0;
      for (const Outline *o : g)
         for (const P2 &p : o->v)
         {
            ctr.x += p.x;
            ctr.y += p.y;
            nv++;
            hx0 = std::min(hx0, p.x); hx1 = std::max(hx1, p.x); hy0 = std::min(hy0, p.y); hy1 = std::max(hy1, p.y);
         }
      ctr = { ctr.x / (double)nv, ctr.y / (double)nv };
      const double rx0 = hx0 - kGridM, rx1 = hx1 + kGridM, ry0 = hy0 - kGridM, ry1 = hy1 + kGridM;
      const auto isTop = [&](const T3 &t) { return Up(t) && std::abs(PlaneZ(t, ctr.x, ctr.y) - topZ) <= kIsTopTolM; };
      const auto overlaps = [](const T3 &t, double x0, double x1, double y0, double y1) {
         return std::max({ t.a.x, t.b.x, t.c.x }) >= x0 && std::min({ t.a.x, t.b.x, t.c.x }) <= x1 && std::max({ t.a.y, t.b.y, t.c.y }) >= y0 && std::min({ t.a.y, t.b.y, t.c.y }) <= y1;
      };
      // an up face below the top within the wall's depth under the outline's box (a thin slab's underside): never refused. Cut by the joined outline at its own height (its rim open by
      // design) where the ball can reach below it: under a hole it can drop through (walls), or one whose floor goes deeper than the face (the groove, a slot's ramp floor: the whole
      // outline is cut, not only the deeper part); else left whole: never reached, and its faces would sit near the ball (user and lead 2026-11-16; not in the paper: the switch test
      // table's playfield_mesh is a 1 VPU slab)
      const auto isUnder = [&](const T3 &t) {
         return Up(t) && !isTop(t) && overlaps(t, hx0, hx1, hy0, hy1) && std::max({ t.a.z, t.b.z, t.c.z }) < topZ - kIsTopTolM && std::min({ t.a.z, t.b.z, t.c.z }) >= topZ - thickM;
      };
      const std::vector<Piece> bd = Pieces(g, false);   // the joined outline's boundary
      std::vector<size_t> inWhich;
      const std::vector<Piece> inner = Pieces(g, true, &inWhich);   // the members' edges inside another member
      for (const Piece &s : bd)
      {
         r.holeAreaD += 0.5 * (s.a.x * s.b.y - s.b.x * s.a.y);   // (r.holeArea: from the stored rim points, after the cut)
      }
      r.inscribed = Inscribed(g, bd);
      for (const Outline *o : g)   // a chamfered circle: the ball passes through its inner hole, not its outer outline
         if (o->passR >= 0.0)
            r.inscribed = std::min(r.inscribed, o->passR);
      const bool walls = r.inscribed >= kBallRM;   // rule 3 (lead 2026-11-16): walls only where the ball can pass through
      dh.walls = walls;
      // the drawing's rim: the boundary pieces chained into loops (each piece's start, in order; a piece's end meets the next one's start within kSnapM: the crossings of two members are
      // found once per member)
      {
         std::vector<char> used(bd.size(), 0);
         for (size_t s0 = 0; s0 < bd.size(); s0++)
         {
            if (used[s0])
               continue;
            std::vector<double> loop;
            std::vector<size_t> ps;   // the loop's pieces in order
            for (size_t s = s0; s < bd.size() && !used[s];)
            {
               used[s] = 1;
               ps.push_back(s);
               loop.push_back(bd[s].a.x / U);
               loop.push_back(bd[s].a.y / U);
               size_t nx = bd.size();
               double best = kSnapM;
               for (size_t k = 0; k < bd.size(); k++)
                  if (const double dd = std::hypot(bd[k].a.x - bd[s].b.x, bd[k].a.y - bd[s].b.y); (!used[k] || k == s0) && dd <= best)
                  {
                     best = dd;
                     nx = k;
                  }
               s = nx;
            }
            dh.rim.push_back(std::move(loop));
            // the drawn rim (drawing only): each corner of a tangent polygon (both its edges tangent to one circle of g_arcs, the corner at most R from its centre) replaced by the
            // true arc from the in edge's touching point to the out edge's (a clipped edge whose touching point is off it: its end, moved onto the circle); an edge with both
            // ends on a slot end's circle (R = r) drawn as the arc between them; every other point as cut. Arc points at most kDrawTolM inside the circle's chord
            constexpr double tol = 1e-8;   // a point within this of the circle (m): the clip's doubles
            std::vector<double> drawn;
            std::vector<char> fl;   // per drawn point: the edge from it is floored (no rim wall)
            const auto put = [&](const P2 &p, const char f) {
               const size_t k = drawn.size();
               if (k >= 2 && std::hypot(p.x - drawn[k - 2] * U, p.y - drawn[k - 1] * U) <= kSnapM)
                  return;
               drawn.push_back(p.x / U);
               drawn.push_back(p.y / U);
               fl.push_back(f);
            };
            const auto arc = [&](const Arc &C, const P2 &f0, const P2 &f1, const bool ends, const char f) {   // f0 to f1 counter-clockwise about C (ends: f0 and f1 too)
               const double t0 = std::atan2(f0.y - C.cy, f0.x - C.cx);
               double d = std::atan2(f1.y - C.cy, f1.x - C.cx) - t0;
               while (d < 0.0)
                  d += 2.0 * kPi;
               int m = std::max(1, (int)std::ceil(d / (2.0 * std::acos(std::max(-1.0, 1.0 - kDrawTolM / C.r)))));
               m += ends ? m % 2 : 0;   // a corner's arc (symmetric about its vertex): even, so one point lies straight in from the vertex (the distance to it is the vertex excess, run 6)
               for (int j = ends ? 0 : 1; j <= (ends ? m : m - 1); j++)
                  put({ C.cx + C.r * std::cos(t0 + d * j / m), C.cy + C.r * std::sin(t0 + d * j / m) }, f);
            };
            const auto touch = [&](const P2 &p, const P2 &q, const Arc &C, P2 &f) {   // the edge p q tangent to C (C on its left); f its touching point, or its end moved onto C
               const double lx = q.x - p.x, ly = q.y - p.y, l2 = lx * lx + ly * ly;
               if (l2 <= 0.0 || std::abs((lx * (C.cy - p.y) - ly * (C.cx - p.x)) / std::sqrt(l2) - C.r) > tol)
                  return false;
               const double t = ((C.cx - p.x) * lx + (C.cy - p.y) * ly) / l2;
               const P2 e = t < 0.0 ? p : q;
               const double de = std::hypot(e.x - C.cx, e.y - C.cy);
               f = t >= 0.0 && t <= 1.0 ? P2 { p.x + t * lx, p.y + t * ly } : P2 { C.cx + C.r * (e.x - C.cx) / de, C.cy + C.r * (e.y - C.cy) / de };
               return true;
            };
            const size_t n = ps.size();
            for (size_t i = 0; i < n; i++)
            {
               const P2 &a = bd[ps[(i + n - 1) % n]].a, &v = bd[ps[i]].a, &b = bd[ps[i]].b;
               const char f = g[bd[ps[i]].m]->floored[bd[ps[i]].e];
               bool done = false;
               for (const Arc &C : g_arcs)
                  if (P2 f0, f1; C.R > C.r && std::hypot(v.x - C.cx, v.y - C.cy) <= C.R + tol && touch(a, v, C, f0) && touch(v, b, C, f1))
                  {
                     arc(C, f0, f1, true, f);
                     done = true;
                     break;
                  }
               if (!done)
                  put(v, f);
               for (const Arc &C : g_arcs)
                  if (C.R == C.r && std::abs(std::hypot(v.x - C.cx, v.y - C.cy) - C.r) <= tol && std::abs(std::hypot(b.x - C.cx, b.y - C.cy) - C.r) <= tol)
                  {
                     arc(C, v, b, false, f);
                     break;
                  }
            }
            if (drawn.size() >= 4 && std::hypot(drawn[0] - drawn[drawn.size() - 2], drawn[1] - drawn.back()) * U <= kSnapM)
            {
               drawn.resize(drawn.size() - 2);
               fl.pop_back();
            }
            dh.drawRim.push_back(std::move(drawn));
            dh.drawFloored.push_back(std::move(fl));
         }
      }
      double along = 0.0;   // the drawing's wood: the rim length so far (the grain along the walls)
      // one triangle as emitted: each corner's place along and depth below the top; kind (meshholes.h DrawnHole::wood) 0 a floor, 1 a rim wall, 2 a wall under a floor's open edge, 3 a side
      const auto woodTri = [&](const T3 &t, const double *u, const double *depth, const int kind) {
         const P3 *p[3] = { &t.a, &t.b, &t.c };
         for (int k = 0; k < 3; k++)
            dh.wood.insert(dh.wood.end(), { (float)(p[k]->x / U), (float)(p[k]->y / U), (float)(p[k]->z / U), (float)(u[k] / U), (float)(depth[k] / U), (float)kind });
      };
      double floorMax = 0.0;   // the group's deepest floor point below the top
      for (const Outline *o : g)
         for (const FTri &t : o->floor)
            floorMax = std::max({ floorMax, t.d[0], t.d[1], t.d[2] });
      const auto cutUnder = [&](const T3 &t) { return isUnder(t) && (walls || topZ - std::max({ t.a.z, t.b.z, t.c.z }) < floorMax); };
      // the top faces under the region (before this cut), and the height of the top at a point
      std::vector<T3> base;
      for (const T3 &t : c.f)
         if (isTop(t) && overlaps(t, rx0, rx1, ry0, ry1))
            base.push_back(t);
      const auto baseZ = [&](const P2 &q, double &z) {
         for (const T3 &t : base)
            if (InTri(Plan(t), q))
            {
               z = PlaneZ(t, q.x, q.y);
               return true;
            }
         return false;
      };
      // refusals (paper section 6): the joined rim must lie on the host's top, and no face of the host that is not its top may cross it within the wall's depth. A smaller existing hole
      // wholly inside is cut over (piece 3a): the host's non-top faces wholly inside in plan go with the cut
      for (const Piece &s : bd)
      {
         const double l = std::hypot(s.b.x - s.a.x, s.b.y - s.a.y);
         if (l <= 0.0)
            continue;
         const P2 in { -(s.b.y - s.a.y) / l * 0.1e-3, (s.b.x - s.a.x) / l * 0.1e-3 };   // left of a counter-clockwise edge: inside
         const int m = std::max(1, (int)std::ceil(l / 0.5e-3));
         for (int k = 0; k <= m; k++)
         {
            const P2 q { s.a.x + (s.b.x - s.a.x) * k / m + in.x, s.a.y + (s.b.y - s.a.y) * k / m + in.y };
            double z;
            if (!baseZ(q, z))
               return "its rim (" + g[s.m]->name + " edge " + std::to_string(s.e) + ") runs off the host's top at (" + std::to_string(q.x * 1e3) + ", " + std::to_string(q.y * 1e3)
                  + ") mm (an existing hole or the mesh's edge crosses the rim)";
         }
      }
      double cover = 0.0;   // the top under the joined outline (an existing hole inside it has none): for the load check's top area
      std::vector<bool> drop(c.f.size(), false);
      for (size_t fi = 0; fi < c.f.size(); fi++)
      {
         const T3 &t = c.f[fi];
         if (!overlaps(t, hx0, hx1, hy0, hy1))
            continue;
         std::vector<P2> p = PlanCCW(t);
         if (isTop(t))
         {
            std::vector<std::array<P2, 3>> left;
            if (!Difference(p, g, left))
               return "the top face at (" + std::to_string(p[0].x * 1e3) + ", " + std::to_string(p[0].y * 1e3) + ") mm could not be cut (its pieces do not close)";
            double al = 0.0;
            for (const auto &q : left)
               al += 0.5 * Cross(q[0], q[1], q[2]);
            cover += std::abs(Area(p)) - al;
            continue;
         }
         if (isUnder(t))
         {
            if (!cutUnder(t))
               continue;
            std::vector<std::array<P2, 3>> left;
            if (!Difference(p, g, left))
               return "the underside face at (" + std::to_string(p[0].x * 1e3) + ", " + std::to_string(p[0].y * 1e3) + ") mm could not be cut (its pieces do not close)";
            double al = 0.0;
            for (const auto &q : left)
               al += 0.5 * Cross(q[0], q[1], q[2]);
            r.underCut += std::abs(Area(p)) - al;
            r.nUnder += std::abs(Area(p)) - al > 1e-14;
            continue;
         }
         const double zlo = std::min({ t.a.z, t.b.z, t.c.z }), zhi = std::max({ t.a.z, t.b.z, t.c.z });
         bool all = true, inside = false;
         for (int k = 0; k < 3; k++)
         {
            const bool in = InAny(g, p[k]);
            all &= in;
            inside |= in;
            for (const Piece &s : bd)
               inside |= SegX(p[k], p[(k + 1) % 3], s.a, s.b);
         }
         if (zlo > topZ + kTopTolM)
            continue;
         if (zhi < topZ - thickM + 0.5e-3)
         {
            r.nDeep += all;   // wholly inside but below the wall's depth: kept (logged)
            continue;
         }
         if (all)
         {
            drop[fi] = true;
            r.nDrop++;
            r.dropUp += Up(t) ? std::abs(Area(p)) : 0.0;   // (the load check's top area counts every up face)
            continue;
         }
         if (inside)
            return "it crosses a face of the host that is not its top (a wall, a step, or a face within the wall's depth) at its rim";
      }
      r.missing = r.holeAreaD - cover;
      const bool underWhole = r.underCut > 0.0 && std::abs(r.underCut - r.holeAreaD) <= 1e-9 * r.holeAreaD;   // the underside under the whole outline: its rim as stored, as the top's (below)
      // the grid: lines 54 mm apart, the origin moved by half a step and then by 1 mm at a time till no line passes within kGridClearM of an outline vertex and no outline edge within
      // kGridClearM of a grid node (the clip's crossings then never fall on a vertex: lead 2026-11-16)
      std::vector<double> xs, ys;
      for (int tries = 0; tries < 54; tries++)
      {
         const double off = 0.5 * kGridM + tries * 1e-3;
         xs = { rx0 };
         ys = { ry0 };
         for (double x = rx0 + off; x < rx1 - 1e-9; x += kGridM)
            xs.push_back(x);
         for (double y = ry0 + off; y < ry1 - 1e-9; y += kGridM)
            ys.push_back(y);
         xs.push_back(rx1);
         ys.push_back(ry1);
         bool clear = true;
         for (const Outline *o : g)
            for (size_t i = 0; i < o->v.size() && clear; i++)
            {
               const P2 &a = o->v[i], &b = o->v[(i + 1) % o->v.size()];
               for (const double x : xs)
                  clear &= std::abs(a.x - x) > kGridClearM;
               for (const double y : ys)
                  clear &= std::abs(a.y - y) > kGridClearM;
               const double l = std::hypot(b.x - a.x, b.y - a.y);
               for (const double x : xs)
                  for (const double y : ys)
                     if (const double u = Param(a, b, { x, y }); clear && l > 0.0 && u > 0.0 && u < 1.0)
                        clear = std::abs(Cross(a, b, { x, y })) / l > kGridClearM;
            }
         r.gridOff = off;
         if (clear)
            break;
      }
      // the cut. The region's existing vertices go into the snap pool first: a new point near one takes it (the faces outside the region keep it: no stray edge, p3 run T2)
      for (const T3 &t : c.f)
         if (overlaps(t, rx0 - kGridM, rx1 + kGridM, ry0 - kGridM, ry1 + kGridM))
            for (const P3 *v : { &t.a, &t.b, &t.c })
               c.Snap(*v);
      std::vector<T3> out;
      for (size_t fi = 0; fi < c.f.size(); fi++)
      {
         const T3 &t = c.f[fi];
         if (drop[fi])
            continue;
         if (!(isTop(t) || cutUnder(t)) || !overlaps(t, rx0, rx1, ry0, ry1))
         {
            out.push_back(t);
            continue;
         }
         std::vector<P2> p = Plan(t);
         const bool flip = Area(p) < 0.0;
         if (flip)
            std::reverse(p.begin(), p.end());
         // outside the region: four convex pieces
         c.Emit(out, ClipLeft(p, { rx0, 0.0 }, { rx0, 1.0 }), t, flip);   // x <= rx0
         c.Emit(out, ClipLeft(p, { rx1, 1.0 }, { rx1, 0.0 }), t, flip);   // x >= rx1
         const std::vector<P2> mid = ClipLeft(ClipLeft(p, { rx0, 1.0 }, { rx0, 0.0 }), { rx1, 0.0 }, { rx1, 1.0 });
         c.Emit(out, ClipLeft(mid, { 1.0, ry0 }, { 0.0, ry0 }), t, flip);   // y <= ry0
         c.Emit(out, ClipLeft(mid, { 0.0, ry1 }, { 1.0, ry1 }), t, flip);   // y >= ry1
         // inside: the grid's cells, minus the joined outline
         for (size_t i = 0; i + 1 < xs.size(); i++)
            for (size_t j = 0; j + 1 < ys.size(); j++)
            {
               const std::vector<P2> cell = Box(p, xs[i], xs[i + 1], ys[j], ys[j + 1]);
               if (cell.size() < 3 || std::abs(Area(cell)) <= 1e-14)
                  continue;
               std::vector<std::array<P2, 3>> tris;
               if (!Difference(cell, g, tris))
                  return "a grid cell at (" + std::to_string(xs[i] * 1e3) + ", " + std::to_string(ys[j] * 1e3) + ") mm could not be cut (its pieces do not close or clip)";
               for (const auto &q : tris)
                  c.Emit(out, { q[0], q[1], q[2] }, t, flip);
            }
      }
      // a rim piece's line on the top: split where it crosses a top face's edge (and a floor face's, for a wall under a floor), each point at its face's height less the floor's depth
      const auto rimLine = [&](const Piece &s, const Outline *under, std::vector<P3> &line) {
         std::vector<std::pair<double, P2>> pts { { 0.0, s.a }, { 1.0, s.b } };
         P2 x;
         const auto cross = [&](const P2 &c0, const P2 &c1) {
            if (Isect(s.a, s.b, c0, c1, x))
               pts.push_back({ Param(s.a, s.b, x), x });
         };
         for (const T3 &t : base)
         {
            const std::vector<P2> p = Plan(t);
            for (int k = 0; k < 3; k++)
               cross(p[k], p[(k + 1) % 3]);
         }
         if (under != nullptr)
            for (const FTri &t : under->floor)
               for (int k = 0; k < 3; k++)
                  cross(t.p[k], t.p[(k + 1) % 3]);
         std::sort(pts.begin(), pts.end(), [](const auto &p, const auto &q) { return p.first < q.first; });
         line.clear();
         for (const auto &[u, q] : pts)
         {
            double z, d = 0.0;
            if (!baseZ(q, z) || (under != nullptr && !FloorDepth(*under, q, d)))
               return false;
            const P3 s3 = c.Snap({ q.x, q.y, z - d });
            if (line.empty() || !Same(line.back(), s3))
               line.push_back(s3);
         }
         return true;
      };
      const auto wall = [&](const std::vector<P3> &line, const P2 &a, const P2 &b) {   // one quad per piece of the line, down by the host's thickness under the top
         for (size_t k = 0; k + 1 < line.size(); k++)
         {
            const P3 &A = line[k], &B = line[k + 1];
            double za, zb;
            baseZ({ A.x, A.y }, za);
            baseZ({ B.x, B.y }, zb);
            const P3 Ab = c.Snap({ A.x, A.y, za - thickM }), Bb = c.Snap({ B.x, B.y, zb - thickM });
            out.push_back({ A, B, Ab });
            out.push_back({ B, Bb, Ab });
            const double a0 = along, a1 = along + std::hypot(B.x - A.x, B.y - A.y);
            along = a1;
            const double u0[3] = { a0, a1, a0 }, d0[3] = { za - A.z, zb - B.z, za - Ab.z }, u1[3] = { a1, a1, a0 }, d1[3] = { zb - B.z, zb - Bb.z, za - Ab.z };
            woodTri(out[out.size() - 2], u0, d0, 1);
            woodTri(out.back(), u1, d1, 1);
            bottoms.push_back({ Ab, Bb });   // the wall's bottom edge: open by design (load check)
            r.nWalls++;
         }
         (void)a; (void)b;
      };
      std::vector<T3> under;   // the underside faces cut (before this cut): the hole's rim in them is open by design
      for (const T3 &t : c.f)
         if (cutUnder(t))
            under.push_back(t);
      for (const Piece &s : bd)
         for (const T3 &t : under)
            if (InTri(Plan(t), { 0.5 * (s.a.x + s.b.x), 0.5 * (s.a.y + s.b.y) }))
            {
               rims.push_back({ c.Snap({ s.a.x, s.a.y, PlaneZ(t, s.a.x, s.a.y) }), c.Snap({ s.b.x, s.b.y, PlaneZ(t, s.b.x, s.b.y) }) });
               break;
            }
      std::vector<P3> line;
      for (const Piece &s : bd)
      {
         if (g[s.m]->floored[s.e])
            continue;   // a floor meets the top here: closed
         if (!rimLine(s, nullptr, line))
            return "a rim point of " + g[s.m]->name + " found no top face under it";
         if (walls)
            wall(line, s.a, s.b);
         else
            for (size_t k = 0; k + 1 < line.size(); k++)
               rims.push_back({ line[k], line[k + 1] });   // no walls: the rim is the top faces' open edge (load check)
      }
      for (size_t i = 0; i < inner.size(); i++)   // a member's edge inside another member that has a floor: a wall from that floor down (walls on), else the floor's open edge
      {
         const Piece &s = inner[i];
         const Outline *under = g[inWhich[i]];
         if (under->floor.empty() || g[s.m]->floored[s.e])
            continue;
         const P2 mid { 0.5 * (s.a.x + s.b.x), 0.5 * (s.a.y + s.b.y) };
         bool third = false;   // (inside a third member too: open there)
         for (size_t k = 0; k < g.size(); k++)
            third |= k != inWhich[i] && k != s.m && InPoly(g[k]->v, mid);
         double dMid;
         if (third || !FloorDepth(*under, mid, dMid))
            continue;   // (over the other member's through part: open, no wall)
         if (!rimLine(s, under, line))
            return "a point of " + g[s.m]->name + "'s edge inside " + under->name + " found no floor under it";
         if (walls)
            wall(line, s.a, s.b);
         else
            for (size_t k = 0; k + 1 < line.size(); k++)
               rims.push_back({ line[k], line[k + 1] });
      }
      // floors: each member's floor faces less the other members, on the top faces they lie over (each piece flat: the top's plane less an affine depth)
      for (size_t m = 0; m < g.size(); m++)
      {
         const Outline &o = *g[m];
         Group others;
         for (size_t k = 0; k < g.size(); k++)
            if (k != m)
               others.push_back(g[k]);
         for (FTri t : o.floor)
         {
            if (Cross(t.p[0], t.p[1], t.p[2]) < 0.0)
            {
               std::swap(t.p[1], t.p[2]);
               std::swap(t.d[1], t.d[2]);
            }
            std::vector<std::array<P2, 3>> pieces;
            if (others.empty())
               pieces.push_back({ t.p[0], t.p[1], t.p[2] });
            else if (!Difference({ t.p[0], t.p[1], t.p[2] }, others, pieces))
               return "a floor face of " + o.name + " could not be cut by the outlines joined to it";
            for (const auto &pc : pieces)
               for (const T3 &b : base)
               {
                  const std::vector<P2> q = ClipConvex({ pc[0], pc[1], pc[2] }, PlanCCW(b));
                  if (q.size() < 3 || std::abs(Area(q)) <= 1e-14)
                     continue;
                  const size_t before = out.size();
                  c.EmitZ(out, q, [&](const P2 &p) { return PlaneZ(b, p.x, p.y) - DepthIn(t, p); }, false);
                  for (size_t k = before; k < out.size(); k++)
                  {
                     const T3 &f = out[k];
                     const double u[3] = { f.a.x + f.a.y, f.b.x + f.b.y, f.c.x + f.c.y };
                     const double dep[3] = { PlaneZ(b, f.a.x, f.a.y) - f.a.z, PlaneZ(b, f.b.x, f.b.y) - f.b.z, PlaneZ(b, f.c.x, f.c.y) - f.c.z };
                     woodTri(f, u, dep, 0);
                     if (Up(out[k]))   // the floor face's plan area as stored (float32 table units), as the measured top area counts it (load check)
                     {
                        std::vector<P2> s = Plan(out[k]);
                        for (P2 &v : s)
                           v = { (double)(float)(v.x / U) * U, (double)(float)(v.y / U) * U };
                        r.floorUp += std::abs(Area(s));
                     }
                     r.nFloor++;
                  }
               }
         }
         for (const FTri &t : o.side)   // vertical faces: kept where they bound wood (not inside another member)
         {
            const P2 g3 { (t.p[0].x + t.p[1].x + t.p[2].x) / 3.0, (t.p[0].y + t.p[1].y + t.p[2].y) / 3.0 };
            if (InAny(others, g3))
               continue;
            P3 v[3];
            double u[3], dep[3];
            for (int k = 0; k < 3; k++)
            {
               double z;
               if (!baseZ(t.p[k], z))
                  return "a side face of " + o.name + " found no top face over it";
               v[k] = c.Snap({ t.p[k].x, t.p[k].y, z - t.d[k] });
               u[k] = v[k].x + v[k].y;
               dep[k] = z - v[k].z;
            }
            out.push_back({ v[0], v[1], v[2] });
            woodTri(out.back(), u, dep, 3);
            r.nFloor++;
         }
         for (const auto &e : o.openFloor)   // a floor edge over a through part: open by design
         {
            if (InAny(others, e[0].first) || InAny(others, e[1].first))
               continue;
            // the edge's line as emitted: split where it crosses a top face's edge (the floor's pieces follow the top's planes), the depth linear along it
            const P2 &ea = e[0].first, &eb = e[1].first;
            std::vector<std::pair<double, P2>> pts { { 0.0, ea }, { 1.0, eb } };
            P2 x;
            for (const T3 &t : base)
            {
               const std::vector<P2> pl = Plan(t);
               for (int k = 0; k < 3; k++)
                  if (Isect(ea, eb, pl[k], pl[(k + 1) % 3], x))
                     pts.push_back({ Param(ea, eb, x), x });
            }
            std::sort(pts.begin(), pts.end(), [](const auto &p, const auto &q) { return p.first < q.first; });
            std::vector<P3> top3, bot3;
            std::vector<double> zs;
            for (const auto &[u, q] : pts)
            {
               double z;
               if (!baseZ(q, z))
                  return "an open floor edge of " + o.name + " found no top face over it";
               const P3 t3 = c.Snap({ q.x, q.y, z - (e[0].second + u * (e[1].second - e[0].second)) });
               if (!top3.empty() && Same(top3.back(), t3))
                  continue;   // a crossing at a point the line already has (snapped to it): no piece between them (a zero-length piece gave two zero-area faces, an open stray edge: brief 5)
               top3.push_back(t3);
               bot3.push_back(c.Snap({ q.x, q.y, z - thickM }));
               zs.push_back(z);
            }
            for (size_t k = 0; k + 1 < top3.size(); k++)
               if (walls)   // the walled through part closed under the floor's edge: an end face down to the walls' bottom (lead 2026-11-16)
               {
                  out.push_back({ top3[k], top3[k + 1], bot3[k] });
                  out.push_back({ top3[k + 1], bot3[k + 1], bot3[k] });
                  const double a0 = along, a1 = along + std::hypot(top3[k + 1].x - top3[k].x, top3[k + 1].y - top3[k].y);
                  along = a1;
                  const double u0[3] = { a0, a1, a0 }, d0[3] = { zs[k] - top3[k].z, zs[k + 1] - top3[k + 1].z, zs[k] - bot3[k].z };
                  const double u1[3] = { a1, a1, a0 }, d1[3] = { zs[k + 1] - top3[k + 1].z, zs[k + 1] - bot3[k + 1].z, zs[k] - bot3[k].z };
                  woodTri(out[out.size() - 2], u0, d0, 2);
                  woodTri(out.back(), u1, d1, 2);
                  bottoms.push_back({ bot3[k], bot3[k + 1] });
                  r.nWalls++;
               }
               else
                  rims.push_back({ top3[k], top3[k + 1] });
         }
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
      // the load check's hole area from the rim points as stored: every vertex the cut left on the joined outline's boundary (grid and face crossings included), in order along each
      // boundary piece, rounded to float32 table units as the faces are stored (the outline's own vertices alone miss the rounding of the points between them: p3 fix round)
      const auto rimArea = [&](const bool top) {
         double area = 0.0;
         for (const Piece &s : bd)
         {
            const double dx = s.b.x - s.a.x, dy = s.b.y - s.a.y, l2 = dx * dx + dy * dy;
            std::vector<std::pair<double, P2>> pts { { 0.0, s.a }, { 1.0, s.b } };
            for (const T3 &t : c.f)
            {
               if (!Up(t) || !overlaps(t, hx0, hx1, hy0, hy1))
                  continue;
               const double zHi = std::max({ t.a.z, t.b.z, t.c.z }), zLo = std::min({ t.a.z, t.b.z, t.c.z });
               if (top ? !isTop(t) : (isTop(t) || zHi >= topZ - kIsTopTolM || zLo < topZ - thickM - kIsTopTolM))
                  continue;
               for (const P3 *v : { &t.a, &t.b, &t.c })
               {
                  const double u = l2 > 0.0 ? ((v->x - s.a.x) * dx + (v->y - s.a.y) * dy) / l2 : 0.0;
                  if (u > 0.0 && u < 1.0 && std::abs(Cross(s.a, s.b, { v->x, v->y })) <= kOnEdgeM * std::sqrt(l2))
                     pts.push_back({ u, { v->x, v->y } });
               }
            }
            std::sort(pts.begin(), pts.end(), [](const auto &p, const auto &q) { return p.first < q.first; });
            std::vector<P2> f32;
            for (const auto &pu : pts)
            {
               const P2 q { (double)(float)(pu.second.x / U) * U, (double)(float)(pu.second.y / U) * U };
               if (f32.empty() || !Eq(f32.back(), q))
                  f32.push_back(q);
            }
            for (size_t k = 0; k + 1 < f32.size(); k++)
               area += 0.5 * (f32[k].x * f32[k + 1].y - f32[k + 1].x * f32[k].y);
         }
         return area;
      };
      r.holeArea = rimArea(true);
      if (underWhole)
         r.underCut = rimArea(false);
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

   void AddOutline(Outline o)   // counter-clockwise (the edge flags follow the edges)
   {
      const size_t n = o.v.size();
      if (o.floored.size() != n)
         o.floored.assign(n, 0);
      if (n >= 3 && Area(o.v) < 0.0)
      {
         std::vector<P2> v(n);
         std::vector<char> f(n);
         for (size_t k = 0; k < n; k++)
         {
            v[k] = o.v[n - 1 - k];
            f[k] = o.floored[(2 * n - 2 - k) % n];   // the new edge k is the old edge n - 2 - k, reversed
         }
         o.v.swap(v);
         o.floored.swap(f);
      }
      g_outlines.push_back(std::move(o));
   }

   // ---- [PlayfieldCut:<name>] (paper section 5) ----
   struct CutKeys
   {
      std::string name, shape;
      double x = 0.0, y = 0.0, angle = 0.0, lengthMm = -1.0, widthMm = -1.0, cornerMm = 0.0, radiusMm = -1.0, baseMm = -1.0, sideMm = -1.0, apexRMm = 0.0, baseRMm = 0.0, chamferDeg = 45.0, outerDMm = -1.0;
      bool xSet = false, ySet = false;
      MeshHoles::HoleKeys keys;   // cut_hole, draw, wall_depth_mm, host (MeshHoles::HoleKey)
   };
   std::vector<CutKeys> Parse(const std::filesystem::path &file)
   {
      std::vector<CutKeys> out;
      std::ifstream in(file);
      std::string line;
      CutKeys *cur = nullptr;
      while (std::getline(in, line))
      {
         const size_t c = line.find_first_of("#;");
         if (c != std::string::npos)
            line.resize(c);
         line = PhysicsIni::Trim(line);
         if (line.empty())
            continue;
         if (line.front() == '[' && line.back() == ']')
         {
            cur = nullptr;
            const std::string head = PhysicsIni::Trim(line.substr(1, line.size() - 2));
            if (PhysicsIni::Lower(head).rfind("playfieldcut:", 0) == 0)
            {
               out.emplace_back();
               cur = &out.back();
               cur->name = PhysicsIni::Trim(head.substr(13));
            }
            continue;
         }
         const size_t eq = line.find('=');
         if (cur == nullptr || eq == std::string::npos)
            continue;
         const std::string key = PhysicsIni::Lower(PhysicsIni::Trim(line.substr(0, eq)));
         const std::string raw = PhysicsIni::Trim(line.substr(eq + 1)), sval = PhysicsIni::Lower(raw);
         const double v = atof(sval.c_str());
         if (key == "shape") cur->shape = sval;
         else if (key == "x") { cur->x = v; cur->xSet = true; }
         else if (key == "y") { cur->y = v; cur->ySet = true; }
         else if (key == "angle") cur->angle = v;
         else if (key == "length_mm") cur->lengthMm = v;
         else if (key == "width_mm") cur->widthMm = v;
         else if (key == "corner_mm") cur->cornerMm = v;
         else if (key == "radius_mm") cur->radiusMm = v;            // circle: the hole's radius (with a chamfer: the inner hole's)
         else if (key == "outer_diameter_mm") cur->outerDMm = v;    // circle: the chamfer's diameter at the top, mm (absent = no chamfer); depth (outer / 2 - radius) tan(chamfer_deg)
         else if (key == "chamfer_deg") cur->chamferDeg = v;        // circle: the chamfer's angle from the top, deg (default 45; above 0, below 90); refused when as deep as the wall depth
         else if (key == "base_mm") cur->baseMm = v;                // triangle: the third side
         else if (key == "side_mm") cur->sideMm = v;                // triangle: the two equal sides
         else if (key == "apex_radius_mm") cur->apexRMm = v;        // triangle: the corner between the equal sides
         else if (key == "base_radius_mm") cur->baseRMm = v;        // triangle: the two corners on the base
         else if (MeshHoles::HoleKey(cur->keys, key, v, raw, "Playfield cut " + cur->name)) {}   // cut_hole, draw, wall_depth_mm, host: the static mesh (part name) the hole is cut in, its top at x, y the hole's height (paper R1; meshholes.h HoleKeys)
      }
      return out;
   }
   // a convex counter-clockwise polygon's corners rounded: each corner of radius r > 0 becomes the tangent polygon of its arc, 15 deg per vertex, the sides tangent where they meet it
   // (vertices at the arc's angles theta + phi (k + 1/2) / n, radius r / cos(phi / 2n): the piece 1 rule, 6 per quarter); "" or why not (the arcs do not fit the sides)
   std::string Rounded(const std::vector<P2> &c, const std::vector<double> &rad, std::vector<P2> &out, std::ostringstream &log)   // log: each rounded corner's edges and excess
   {
      const size_t n = c.size();
      for (size_t i = 0; i < n; i++)   // the arcs' tangent lengths must fit each side
      {
         const auto tl = [&](size_t k) {
            const P2 &a = c[(k + n - 1) % n], &b = c[k], &d = c[(k + 1) % n];
            const double ang = std::acos(std::clamp(((a.x - b.x) * (d.x - b.x) + (a.y - b.y) * (d.y - b.y)) / (std::hypot(a.x - b.x, a.y - b.y) * std::hypot(d.x - b.x, d.y - b.y)), -1.0, 1.0));
            return rad[k] / std::tan(0.5 * ang);
         };
         if (tl(i) + tl((i + 1) % n) > std::hypot(c[(i + 1) % n].x - c[i].x, c[(i + 1) % n].y - c[i].y) + 1e-12)
            return "its corner radii do not fit side " + std::to_string(i);
      }
      out.clear();
      for (size_t i = 0; i < n; i++)
      {
         const P2 &a = c[(i + n - 1) % n], &b = c[i], &d = c[(i + 1) % n];
         if (rad[i] <= 0.0)
         {
            out.push_back(b);
            continue;
         }
         const double l1 = std::hypot(b.x - a.x, b.y - a.y), l2 = std::hypot(d.x - b.x, d.y - b.y);
         const P2 n1 { (b.y - a.y) / l1, -(b.x - a.x) / l1 }, n2 { (d.y - b.y) / l2, -(d.x - b.x) / l2 };   // outward normals
         const double k = -rad[i] / (1.0 + n1.x * n2.x + n1.y * n2.y);
         const P2 ctr { b.x + k * (n1.x + n2.x), b.y + k * (n1.y + n2.y) };   // the arc's centre, r inside both sides
         const double th = std::atan2(n1.y, n1.x), phi = std::atan2(n1.x * n2.y - n1.y * n2.x, n1.x * n2.x + n1.y * n2.y);
         // paper section 4 'quarter arcs at the same rule' as a circle (Round): at least 6 per quarter, more till the vertex stands at most hole_round_mm outside the arc
         // (lead 2026-11-17: the fixed 15 deg step gave a 25 mm corner 0.216 mm, above the 0.2 mm default)
         int m = std::max(1, (int)std::lround(phi / kTurnRad));
         while (rad[i] * (1.0 / std::cos(0.5 * phi / m) - 1.0) > g_excessM)
            m++;
         const double R = rad[i] / std::cos(0.5 * phi / m);
         log << (log.tellp() > 0 ? ", " : "") << "corner " << i << ": radius " << rad[i] * 1e3 << " mm, " << m << " edges, excess " << (R - rad[i]) * 1e3 << " mm";
         g_arcs.push_back({ ctr.x, ctr.y, rad[i], R });   // (the drawing's true curve)
         for (int j = 0; j < m; j++)
            out.push_back({ ctr.x + R * std::cos(th + phi * (j + 0.5) / m), ctr.y + R * std::sin(th + phi * (j + 0.5) / m) });
      }
      return std::string();
   }
   // the host's face list replaced by f (m), welded by their float bits as every host's (contactsurfacehost.cpp RubberSurface::Add), linked again and its box recomputed; returns the
   // edges on more than two faces (Link)
   int Store(ContactSurfaceHost *const host, const std::vector<T3> &f)
   {
      RubberSurface &s = *host->MeshSurface();
      RubberSurface ns;
      for (const T3 &t : f)
         ns.Add(Vertex3Ds((float)(t.a.x / U), (float)(t.a.y / U), (float)(t.a.z / U)), Vertex3Ds((float)(t.b.x / U), (float)(t.b.y / U), (float)(t.b.z / U)),
                Vertex3Ds((float)(t.c.x / U), (float)(t.c.y / U), (float)(t.c.z / U)));
      s.tris = std::move(ns.tris);
      std::copy(ns.lo, ns.lo + 3, s.lo);
      std::copy(ns.hi, ns.hi + 3, s.hi);
      s.vid = std::move(ns.vid);
      s.weld = std::move(ns.weld);
      s.nb.clear();
      s.fold.clear();
      return host->MeshRelink();   // shared edges again (contactsurfacehost.cpp Link) and the host's box
   }

   // ---- overlap repair ([Features] mesh_overlap_repair; to-do N3 piece 5) ----
   constexpr double kPlaneM = 1e-6;       // two faces on one plane: every corner of each within this of the other's plane (user and lead 2026-11-16)
   constexpr double kOverlapM2 = 1e-12;   // two faces overlap when their common area exceeds this (1 um2; faces that only share an edge have none)
   P3 Sub(const P3 &a, const P3 &b) { return { a.x - b.x, a.y - b.y, a.z - b.z }; }
   double Dot(const P3 &a, const P3 &b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
   P3 Unit(const P3 &a) { const double l = std::sqrt(Dot(a, a)); return l > 0.0 ? P3 { a.x / l, a.y / l, a.z / l } : P3 { 0.0, 0.0, 0.0 }; }
   P3 CrossP(const P3 &a, const P3 &b) { return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x }; }
   // a face's own plane frame: origin its first corner, u v in the plane with u x v = its normal
   struct Frame
   {
      P3 o, u, v, n;
      explicit Frame(const T3 &t)
      {
         o = t.a;
         n = Unit(Normal(t));
         u = Unit(Sub(t.b, t.a));
         v = CrossP(n, u);
      }
      P2 To(const P3 &p) const { const P3 d = Sub(p, o); return { Dot(d, u), Dot(d, v) }; }
      P3 From(const P2 &q) const { return { o.x + q.x * u.x + q.y * v.x, o.y + q.x * u.y + q.y * v.y, o.z + q.x * u.z + q.y * v.z }; }
   };
   bool OnePlane(const T3 &a, const T3 &b)
   {
      const P3 na = Unit(Normal(a)), nb = Unit(Normal(b));
      if (std::abs(Dot(na, nb)) < 0.999)
         return false;
      for (const P3 *p : { &b.a, &b.b, &b.c })
         if (std::abs(Dot(Sub(*p, a.a), na)) > kPlaneM)
            return false;
      for (const P3 *p : { &a.a, &a.b, &a.c })
         if (std::abs(Dot(Sub(*p, b.a), nb)) > kPlaneM)
            return false;
      return true;
   }
   std::vector<P2> InFrame(const Frame &fr, const T3 &t)   // counter-clockwise in the frame
   {
      std::vector<P2> p { fr.To(t.a), fr.To(t.b), fr.To(t.c) };
      if (Area(p) < 0.0)
         std::reverse(p.begin(), p.end());
      return p;
   }
   double Area3(const T3 &t) { const P3 n = Normal(t); return 0.5 * std::sqrt(Dot(n, n)); }
}

void MeshHoles::RepairOverlaps(const bool on)
{
   if (!on)
      return;
   for (ContactHost *const h : ContactHost::AllHosts())
   {
      ContactSurfaceHost *const host = dynamic_cast<ContactSurfaceHost *>(h);
      if (host == nullptr || host->MeshSurface() == nullptr || !IsMeshKind(host->MeshKind()))
         continue;
      const auto t0 = std::chrono::steady_clock::now();
      const RubberSurface &s = *host->MeshSurface();
      std::vector<T3> f;
      f.reserve(s.tris.size());
      for (const RubberContact::SurfaceTri &t : s.tris)
         f.push_back({ { t.a.x, t.a.y, t.a.z }, { t.b.x, t.b.y, t.b.z }, { t.c.x, t.c.y, t.c.z } });
      const size_t n = f.size();
      if (n < 2)
         continue;
      // the pairs on one plane that overlap, found through a grid of boxes over the host (a face goes into every cell its box meets; at most 64 cells per axis, at least 2 mm)
      double lo[3] = { 1e30, 1e30, 1e30 }, hi[3] = { -1e30, -1e30, -1e30 };
      for (const T3 &t : f)
         for (const P3 *p : { &t.a, &t.b, &t.c })
            for (int k = 0; k < 3; k++)
            {
               lo[k] = std::min(lo[k], (&p->x)[k]);
               hi[k] = std::max(hi[k], (&p->x)[k]);
            }
      const double cell = std::max({ (hi[0] - lo[0]) / 64.0, (hi[1] - lo[1]) / 64.0, (hi[2] - lo[2]) / 64.0, 2e-3 });
      int nc[3];
      for (int k = 0; k < 3; k++)
         nc[k] = std::max(1, std::min(64, (int)std::ceil((hi[k] - lo[k]) / cell)));
      const auto at = [&](const double x, const int k) { return std::max(0, std::min(nc[k] - 1, (int)std::floor((x - lo[k]) / cell))); };
      std::vector<std::vector<int>> grid((size_t)nc[0] * nc[1] * nc[2]);
      std::vector<std::array<int, 6>> span(n);
      for (size_t i = 0; i < n; i++)
      {
         const T3 &t = f[i];
         for (int k = 0; k < 3; k++)
         {
            span[i][k] = at(std::min({ (&t.a.x)[k], (&t.b.x)[k], (&t.c.x)[k] }) - kPlaneM, k);
            span[i][k + 3] = at(std::max({ (&t.a.x)[k], (&t.b.x)[k], (&t.c.x)[k] }) + kPlaneM, k);
         }
         for (int x = span[i][0]; x <= span[i][3]; x++)
            for (int y = span[i][1]; y <= span[i][4]; y++)
               for (int z = span[i][2]; z <= span[i][5]; z++)
                  grid[((size_t)x * nc[1] + y) * nc[2] + z].push_back((int)i);
      }
      std::vector<std::vector<int>> nbr(n);   // the faces each overlaps on its plane
      std::vector<int> stamp(n, -1);
      std::vector<double> area(n);
      for (size_t i = 0; i < n; i++)
         area[i] = Area3(f[i]);
      size_t pairs = 0;
      for (size_t i = 0; i < n; i++)
      {
         if (area[i] <= 0.0)
            continue;
         const Frame fr(f[i]);
         const std::vector<P2> pi = InFrame(fr, f[i]);
         for (int x = span[i][0]; x <= span[i][3]; x++)
            for (int y = span[i][1]; y <= span[i][4]; y++)
               for (int z = span[i][2]; z <= span[i][5]; z++)
                  for (const int j : grid[((size_t)x * nc[1] + y) * nc[2] + z])
                  {
                     if (j <= (int)i || stamp[j] == (int)i || area[j] <= 0.0)
                        continue;
                     stamp[j] = (int)i;
                     if (!OnePlane(f[i], f[j]))
                        continue;
                     const std::vector<P2> common = ClipConvex(InFrame(fr, f[j]), pi);
                     if (common.size() < 3 || std::abs(Area(common)) <= kOverlapM2)
                        continue;
                     nbr[i].push_back(j);
                     nbr[(size_t)j].push_back((int)i);
                     pairs++;
                  }
      }
      const double ms0 = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
      if (pairs == 0)
      {
         PLOGI << "Mesh overlap repair: " << host->MeshKind() << " host " << host->Name() << ": faces " << n << ", no overlapping faces, " << ms0 << " ms";
         continue;
      }
      // the faces in an overlap, smallest first (ties: load order; lead 2026-11-16): each keeps what the faces before it that it overlaps do not cover (a face wholly covered is dropped).
      // Those faces' union is the union of what they kept, so each subtracts their whole faces: F less a triangle T = the parts of F outside T's edge k and inside its edges before k
      // (at most 3, disjoint; exact where an edge of T lies on one of F)
      std::vector<int> order;
      for (size_t i = 0; i < n; i++)
         if (!nbr[i].empty())
            order.push_back((int)i);
      std::sort(order.begin(), order.end(), [&](const int a, const int b) { return area[a] != area[b] ? area[a] < area[b] : a < b; });
      std::vector<char> done(n, 0), changed(n, 0);
      std::vector<std::vector<T3>> kept(n);
      std::vector<P3> pool;   // the points of the repair: a new point near one of these takes it
      for (const int i : order)
         for (const P3 *p : { &f[i].a, &f[i].b, &f[i].c })
            pool.push_back(*p);
      const size_t nOrig = pool.size();
      double areaBefore = 0.0, areaAfter = 0.0, removed = 0.0;
      for (size_t i = 0; i < n; i++)
         areaBefore += area[i];
      std::ostringstream where;
      int nDropped = 0, nClipped = 0;
      for (const int i : order)
      {
         done[i] = 1;
         std::vector<int> under;   // the faces before it that it overlaps
         for (const int j : nbr[i])
            if (done[j] && j != i)
               under.push_back(j);
         if (under.empty())
            continue;   // (kept whole)
         const Frame fr(f[i]);
         std::vector<std::vector<P2>> pieces { InFrame(fr, f[i]) };
         for (const int j : under)
         {
            const std::vector<P2> T = InFrame(fr, f[j]);
            std::vector<std::vector<P2>> next;
            for (const std::vector<P2> &P : pieces)
            {
               std::vector<P2> inside = P;
               for (int k = 0; k < 3 && inside.size() >= 3; k++)
               {
                  const std::vector<P2> out = ClipLeft(inside, T[(k + 1) % 3], T[k]);   // outside edge k (on its right)
                  if (out.size() >= 3 && std::abs(Area(out)) > kOverlapM2)
                     next.push_back(out);
                  inside = ClipLeft(inside, T[k], T[(k + 1) % 3]);
               }
            }
            pieces.swap(next);
         }
         changed[i] = 1;
         // a corner at a point of the faces it was clipped by (or of itself) takes that point; one on such a face's edge takes that edge's own point there (no T-junction across the
         // small step between the two planes); else it lies on its own plane. Corners within the T-pass tolerance of a straight run go (as EmitZ)
         const auto place = [&](const P2 &q) {
            for (const P3 &w3 : pool)
            {
               const P2 w = fr.To(w3);
               if ((w.x - q.x) * (w.x - q.x) + (w.y - q.y) * (w.y - q.y) < kSnapM * kSnapM && std::abs(Dot(Sub(w3, fr.o), fr.n)) <= kPlaneM)
                  return w3;
            }
            P3 p = fr.From(q);
            bool onEdge = false;
            for (size_t m = 0; m < under.size() && !onEdge; m++)
            {
               const T3 &t = f[under[m]];
               const P3 c3[3] = { t.a, t.b, t.c };
               for (int e = 0; e < 3 && !onEdge; e++)
               {
                  const P2 a = fr.To(c3[e]), b = fr.To(c3[(e + 1) % 3]);
                  const double l = std::hypot(b.x - a.x, b.y - a.y), u = l > 0.0 ? Param(a, b, q) : -1.0;
                  if (u > 0.0 && u < 1.0 && std::abs(Cross(a, b, q)) <= kOnEdgeM * l)
                  {
                     const P3 &A = c3[e], &B = c3[(e + 1) % 3];
                     p = { A.x + u * (B.x - A.x), A.y + u * (B.y - A.y), A.z + u * (B.z - A.z) };
                     onEdge = true;
                  }
               }
            }
            pool.push_back(p);
            return p;
         };
         double left = 0.0;
         for (const std::vector<P2> &P : pieces)
         {
            std::vector<P2> q2;
            for (const P2 &q : P)
               if (q2.empty() || !Eq(q2.back(), q))
                  q2.push_back(q);
            while (q2.size() > 1 && Eq(q2.front(), q2.back()))
               q2.pop_back();
            for (bool again = true; again && q2.size() >= 3;)
            {
               again = false;
               for (size_t k = 0; k < q2.size(); k++)
               {
                  const P2 &a = q2[(k + q2.size() - 1) % q2.size()], &b = q2[k], &c = q2[(k + 1) % q2.size()];
                  const double lac = std::hypot(c.x - a.x, c.y - a.y);
                  if (lac <= 0.0 || std::abs(Cross(a, b, c)) <= kOnEdgeM * lac)
                  {
                     q2.erase(q2.begin() + (std::ptrdiff_t)k);
                     again = true;
                     break;
                  }
               }
            }
            if (q2.size() < 3)
               continue;
            std::vector<P3> q3;
            for (const P2 &q : q2)
               q3.push_back(place(q));
            const bool flip = Dot(Normal({ q3[0], q3[1], q3[2] }), Normal(f[i])) < 0.0;   // (the face's own winding kept)
            for (size_t k = 1; k + 1 < q3.size(); k++)
            {
               const T3 t = flip ? T3 { q3[0], q3[k + 1], q3[k] } : T3 { q3[0], q3[k], q3[k + 1] };
               kept[i].push_back(t);
               left += Area3(t);
            }
         }
         const P3 g { (f[i].a.x + f[i].b.x + f[i].c.x) / 3.0, (f[i].a.y + f[i].b.y + f[i].c.y) / 3.0, (f[i].a.z + f[i].b.z + f[i].c.z) / 3.0 };
         (kept[i].empty() ? nDropped : nClipped)++;
         removed += area[i] - left;
         where << " face " << i << " centre (" << g.x * 1e3 << ", " << g.y * 1e3 << ", " << g.z * 1e3 << ") mm [" << g.x / U << ", " << g.y / U << " VPU] "
               << (kept[i].empty() ? std::string("dropped") : "clipped to " + std::to_string(kept[i].size()) + " faces") << ", under " << under.size() << " smaller faces, double cover removed "
               << (area[i] - left) * 1e6 << " mm2;";
      }
      std::vector<T3> out;
      for (size_t i = 0; i < n; i++)
         if (!changed[i])
            out.push_back(f[i]);
         else
            out.insert(out.end(), kept[i].begin(), kept[i].end());
      const std::vector<P3> fresh(pool.begin() + (std::ptrdiff_t)nOrig, pool.end());
      NoTJunction(out, fresh);
      for (const T3 &t : out)
         areaAfter += Area3(t);
      const int over = Store(host, out);
      int open = 0;
      for (const auto &nb : host->MeshSurface()->nb)
         for (int k = 0; k < 3; k++)
            open += nb[k] < 0;
      PLOGI << "Mesh overlap repair: " << host->MeshKind() << " host " << host->Name() << ": faces " << n << " -> " << out.size() << ", overlapping pairs " << pairs << ", faces dropped " << nDropped << ", clipped "
            << nClipped << ", face area " << areaBefore * 1e6 << " -> " << areaAfter * 1e6 << " mm2 (double cover removed " << removed * 1e6 << " mm2), open edges after " << open
            << ", edges on more than two faces after " << over << ", " << std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count() << " ms;" << where.str();
   }
}

void MeshHoles::Clear(const double roundMm)
{
   double mm = roundMm;
   if (!(mm >= 0.01 && mm <= 1.0))
   {
      mm = mm > 1.0 ? 1.0 : 0.01;
      std::ostringstream n;
      n << "Mesh holes: hole_round_mm " << roundMm << " is outside 0.01 .. 1: using " << mm << " (physics ini [Features])";
      Notice(n.str());
   }
   g_excessM = mm * 1e-3;
   PLOGI << "Mesh holes: hole_round_mm " << mm << " (a round hole's vertex at most this far outside its circle)";
   g_drawn.clear();
   g_arcs.clear();
   g_outlines.clear();   // (the Stars, Buttons and physical kickers add theirs in PlaceParts, the turntables in Turntables::Setup, the [PlayfieldCut:<name>] sections in AddIniCuts, all after this)
   g_place.clear();
}

void MeshHoles::Defer(std::function<void()> place)
{
   g_place.push_back(std::move(place));
}

void MeshHoles::PlaceParts()
{
   std::vector<std::function<void()>> place;
   place.swap(g_place);   // (emptied before the calls: no part pointer kept past the load)
   for (const std::function<void()> &f : place)
      f();
}

MeshHoles::Top MeshHoles::HostTop(const std::string &host, const double xU, const double yU, double &zU, double &thickM)
{
   // a named ramp host: not cut (Build says why for a ramp found at the rim; here before its top is looked up, so the refusal names the reason, not a missing top)
   if (!host.empty() && std::any_of(ContactHost::AllHosts().begin(), ContactHost::AllHosts().end(), [&](ContactHost *const h) {
          const ContactSurfaceHost *const s = dynamic_cast<ContactSurfaceHost *>(h);
          return s != nullptr && s->MeshKind() == "ramp" && PhysicsIni::Lower(s->Name()) == PhysicsIni::Lower(host);
       }))
      return Top::Ramp;
   bool found = false, named = false;
   for (ContactHost *const h : ContactHost::AllHosts())
   {
      const ContactSurfaceHost *const s = dynamic_cast<ContactSurfaceHost *>(h);
      const bool want = host.empty() ? h->IsPlayfieldKind()
                                     : s != nullptr && s->MeshSurface() != nullptr && IsMeshKind(s->MeshKind()) && PhysicsIni::Lower(s->Name()) == PhysicsIni::Lower(host);
      named |= want;
      if (want)
         if (double z; h->SurfaceHeightAt(xU, yU, z) && (!found || z > zU))
         {
            zU = z;
            found = true;
            thickM = (s != nullptr ? s->ThicknessMm() : RubberObjectSettings().thicknessMm) * 1e-3;
         }
   }
   return found ? Top::Found : named ? Top::NoTop : Top::NoMesh;
}

bool MeshHoles::PartOnHost(const std::string &part, const std::string &host, const double xU, const double yU, double &zU)
{
   double z = 0.0, thickM = 0.0;
   const Top t = HostTop(host, xU, yU, z, thickM);
   if (t == Top::Found)
   {
      zU = z;
      PLOGI << "Mesh holes: " << part << ": on host " << host << ", its top at x, y z " << z << " VPU";
      return true;
   }
   Notice(part + (t == Top::Ramp ? ": host ramp " + host + ": holes in ramps are not supported yet: hole not cut; the part stays at its own height. Name a playfield, wall or primitive in host (physics ini)"
                  : t == Top::NoMesh ? ": no static mesh named '" + host + "' (physics ini host): hole not cut; the part stays at its own height. Set host to a playfield, wall or primitive with collision"
                                     : ": host '" + host + "' has no top under its x, y: hole not cut; the part stays at its own height. Fix host or the part's place"));
   return false;
}

bool MeshHoles::HoleKey(HoleKeys &h, const std::string &lowerKey, const double v, const std::string &text, const std::string &where)
{
   if (lowerKey == "host")
      h.host = text;
   else if (lowerKey == "cut_hole")
      h.cut = v != 0.0;
   else if (lowerKey == "draw")
      h.draw = v != 0.0;
   else if (lowerKey == "wall_depth_mm")
   {
      h.wallDepthMm = v;
      if (!(v > 0.0))
      {
         h.wallDepthMm = -1.0;
         Notice(where + ": wall_depth_mm " + std::to_string(v) + " is not above 0: the host's thickness is used (physics ini)");
      }
   }
   else
      return false;
   return true;
}

void MeshHoles::Add(const std::string &name, const std::vector<double> &xy, double partZ, const HoleKeys &keys)
{
   if (!keys.cut)   // cut_hole = 0: no outline (every part's hole goes through here: one place for the key)
   {
      PLOGI << "Mesh holes: " << name << ": cut_hole = 0, no hole";
      return;
   }
   Outline o;
   o.name = name;
   o.z = partZ;
   o.host = keys.host;   // (named: cut in that mesh, Build's named-host branch)
   o.keys = keys;
   PLOGI << "Mesh holes: " << name << ": cut_hole 1, draw " << (keys.draw ? 1 : 0) << ", wall_depth_mm " << (keys.wallDepthMm > 0.0 ? std::to_string(keys.wallDepthMm) : std::string("host thickness"))
         << (keys.host.empty() ? std::string() : ", host " + keys.host);
   for (size_t i = 0; i + 1 < xy.size(); i += 2)
      o.v.push_back({ xy[i], xy[i + 1] });
   AddOutline(std::move(o));
}

std::vector<double> MeshHoles::Round(double x, double y, double r, int n)
{
   if (n <= 0)
   {
      n = 24;
      while (r * (1.0 / std::cos(kPi / n) - 1.0) > g_excessM)
         n++;
   }
   const double R = r / std::cos(kPi / n);   // vertices outside the circle: every edge touches it at its midpoint (paper section 4)
   g_arcs.push_back({ x, y, r, R });   // (the drawing's true curve)
   std::vector<double> xy;
   for (int k = 0; k < n; k++)
   {
      xy.push_back(x + R * std::cos(2.0 * kPi * k / n));
      xy.push_back(y + R * std::sin(2.0 * kPi * k / n));
   }
   return xy;
}

std::vector<double> MeshHoles::Offset(const std::vector<double> &xy, double d)
{
   std::vector<P2> p;
   for (size_t i = 0; i + 1 < xy.size(); i += 2)
      p.push_back({ xy[i], xy[i + 1] });
   if (p.size() >= 3 && Area(p) < 0.0)
      std::reverse(p.begin(), p.end());
   const size_t n = p.size();
   std::vector<double> o;
   for (size_t i = 0; i < n; i++)   // vertex i: where the offset lines of edges i-1 and i cross
   {
      const P2 &a = p[(i + n - 1) % n], &b = p[i], &c = p[(i + 1) % n];
      const double l1 = std::hypot(b.x - a.x, b.y - a.y), l2 = std::hypot(c.x - b.x, c.y - b.y);
      const P2 n1 { (b.y - a.y) / l1, -(b.x - a.x) / l1 }, n2 { (c.y - b.y) / l2, -(c.x - b.x) / l2 };   // outward normals (right of a counter-clockwise edge)
      const double k = d / (1.0 + n1.x * n2.x + n1.y * n2.y);   // b + d (n1 + n2) / (1 + n1 . n2) lies at distance d from both lines
      o.push_back(b.x + k * (n1.x + n2.x));
      o.push_back(b.y + k * (n1.y + n2.y));
   }
   return o;
}

void MeshHoles::AddIniCuts(const std::filesystem::path &file)
{
   if (file.empty() || !ContactSystem::Unified())   // (stock physics: no static mesh host to cut)
      return;
   std::vector<CutKeys> cuts = Parse(file);
   if (cuts.empty())
      return;
   for (const CutKeys &k : cuts)
   {
      const std::string where = "Playfield cut " + k.name + ": ";
      if (!k.xSet || !k.ySet)
      {
         Notice(where + "needs x and y: the cut is ignored (physics ini)");
         continue;
      }
      if (!k.keys.cut)
      {
         PLOGI << "Mesh holes: Cut " << k.name << ": cut_hole = 0, no hole";
         continue;
      }
      const std::string &host = k.keys.host;
      // the hole's height: with host named, that host's top at x, y (its highest up face there); else the highest playfield surface at x, y (HostTop, the parts' host key's lookup)
      double zU = 0.0, hostThickM = 0.0;   // (hostThickM: that host's thickness_mm, a chamfer's depth limit)
      const Top top = HostTop(host, k.x, k.y, zU, hostThickM);
      if (top == Top::Ramp)
      {
         Notice(where + "host ramp " + host + ": holes in ramps are not supported yet: hole not cut");
         continue;
      }
      if (top != Top::Found)
      {
         Notice(where + (host.empty() ? std::string("no playfield surface at its centre") : "the host '" + host + "' has no top under x, y") + ": the cut is ignored (physics ini x, y, host)");
         continue;
      }
      Outline o;
      o.name = "Cut " + k.name;
      o.z = zU * U;
      o.host = host;
      o.keys = k.keys;
      PLOGI << "Mesh holes: " << o.name << ": cut_hole 1, draw " << (k.keys.draw ? 1 : 0) << ", wall_depth_mm " << (k.keys.wallDepthMm > 0.0 ? std::to_string(k.keys.wallDepthMm) : std::string("host thickness"))
            << ", height " << zU << " VPU (" << (host.empty() ? "highest playfield surface" : "host " + host + "'s top") << " at x, y)";
      const double rot = k.angle * (kPi / 180.0), cx = k.x * U, cy = k.y * U;
      const P2 l { -std::sin(rot), std::cos(rot) }, p { -l.y, l.x };   // along (the trigger's rotation convention, wireshape.h Place) and across
      const auto W = [&](double s, double t) { return P2 { cx + l.x * s + p.x * t, cy + l.y * s + p.y * t }; };
      std::ostringstream log;
      if (k.shape == "circle")
      {
         if (k.radiusMm <= 0.0)
         {
            Notice(where + "a circle needs radius_mm above 0: the cut is ignored (physics ini)");
            continue;
         }
         if (k.outerDMm < 0.0)
         {
            const std::vector<double> xy = MeshHoles::Round(cx, cy, k.radiusMm * 1e-3);
            for (size_t i = 0; i + 1 < xy.size(); i += 2)
               o.v.push_back({ xy[i], xy[i + 1] });
            log << "circle radius " << k.radiusMm << " mm";
         }
         else
         {
            // a chamfered circle: a cone from the outer circle (outer_diameter_mm / 2, at the top) down to the inner hole (radius_mm) at depth (outer / 2 - r) tan(chamfer_deg), then
            // the inner hole's walls (when the ball can pass, the walls rule on the outer outline) down to the wall depth. Both rings tangent polygons with the outer ring's N (each
            // edge touches its circle at its midpoint); the cone: one face pair per edge, its inner edge an open floor edge over the through hole (as the slot's ramp floor's)
            const double ro = 0.5 * k.outerDMm * 1e-3, ri = k.radiusMm * 1e-3, depth = (ro - ri) * std::tan(k.chamferDeg * (kPi / 180.0));
            const double limM = k.keys.wallDepthMm > 0.0 ? k.keys.wallDepthMm * 1e-3 : hostThickM;
            if (!(ro > ri) || !(k.chamferDeg > 0.0 && k.chamferDeg < 90.0))
            {
               Notice(where + "a chamfer needs outer_diameter_mm above 2 radius_mm and chamfer_deg between 0 and 90: the cut is ignored (physics ini)");
               continue;
            }
            if (depth >= limM)
            {
               std::ostringstream n;
               n << where << "the chamfer is " << depth * 1e3 << " mm deep, not less than the " << (k.keys.wallDepthMm > 0.0 ? "wall depth " : "host's thickness ") << limM * 1e3 << " mm: the cut is ignored (physics ini)";
               Notice(n.str());
               continue;
            }
            const std::vector<double> xo = MeshHoles::Round(cx, cy, ro);
            const int nr = (int)(xo.size() / 2);
            const std::vector<double> xi = MeshHoles::Round(cx, cy, ri, nr);
            o.passR = ri;
            for (int i = 0; i < nr; i++)
            {
               const int j = (i + 1) % nr;
               const P2 a { xo[2 * i], xo[2 * i + 1] }, b { xo[2 * j], xo[2 * j + 1] }, ai { xi[2 * i], xi[2 * i + 1] }, bi { xi[2 * j], xi[2 * j + 1] };
               o.v.push_back(a);
               o.floored.push_back(1);
               o.floor.push_back({ { a, b, bi }, { 0.0, 0.0, depth } });
               o.floor.push_back({ { a, bi, ai }, { 0.0, depth, depth } });
               o.openFloor.push_back({ std::make_pair(ai, depth), std::make_pair(bi, depth) });
            }
            log << "circle radius " << k.radiusMm << " mm, chamfer " << k.chamferDeg << " deg to outer diameter " << k.outerDMm << " mm (" << depth * 1e3 << " mm deep, " << nr << " edges on both rings)";
         }
      }
      else if (k.shape == "rectangle" || k.shape == "triangle")
      {
         std::vector<P2> c;
         std::vector<double> rad;
         if (k.shape == "rectangle")
         {
            if (k.lengthMm <= 0.0 || k.widthMm <= 0.0)
            {
               Notice(where + "a rectangle needs length_mm and width_mm above 0: the cut is ignored (physics ini)");
               continue;
            }
            const double hl = 0.5 * k.lengthMm * 1e-3, hw = 0.5 * k.widthMm * 1e-3;
            double cr = k.cornerMm * 1e-3;
            if (cr < 0.0 || cr > std::min(hl, hw))
            {
               const double lim = cr < 0.0 ? 0.0 : std::min(hl, hw);
               std::ostringstream n;
               n << where << "corner_mm " << k.cornerMm << " is outside 0 .. " << lim * 1e3 << ": using " << lim * 1e3 << " (physics ini)";
               Notice(n.str());
               cr = lim;
            }
            c = { W(-hl, -hw), W(hl, -hw), W(hl, hw), W(-hl, hw) };
            rad.assign(4, cr);
            log << "rectangle " << k.lengthMm << " x " << k.widthMm << " mm, corner " << cr * 1e3 << " mm";
         }
         else
         {
            const double b = k.baseMm * 1e-3, sd = k.sideMm * 1e-3;
            if (b <= 0.0 || sd <= 0.5 * b)
            {
               Notice(where + "a triangle needs base_mm above 0 and side_mm above half of it: the cut is ignored (physics ini)");
               continue;
            }
            const double h = std::sqrt(sd * sd - 0.25 * b * b);   // x, y: the centroid; the apex along the angle's direction
            c = { W(-h / 3.0, -0.5 * b), W(2.0 * h / 3.0, 0.0), W(-h / 3.0, 0.5 * b) };
            rad = { std::max(0.0, k.baseRMm) * 1e-3, std::max(0.0, k.apexRMm) * 1e-3, std::max(0.0, k.baseRMm) * 1e-3 };
            log << "triangle base " << k.baseMm << " mm, sides " << k.sideMm << " mm, apex radius " << rad[1] * 1e3 << " mm, base radius " << rad[0] * 1e3 << " mm";
         }
         if (Area(c) < 0.0)
         {
            std::reverse(c.begin(), c.end());
            std::reverse(rad.begin(), rad.end());
         }
         std::ostringstream corners;
         if (const std::string why = Rounded(c, rad, o.v, corners); !why.empty())
         {
            Notice(where + why + ": the cut is ignored (physics ini)");
            continue;
         }
         if (!corners.str().empty())
            log << " (" << corners.str() << ")";
      }
      else if (k.shape == "slot")
      {
         if (k.lengthMm < 0.0)
         {
            Notice(where + "a slot needs length_mm (0 or more): the cut is ignored (physics ini)");
            continue;
         }
         double widthMm = k.widthMm < 0.0 ? kWidthMaxMm : k.widthMm;
         if (widthMm > kWidthMaxMm || widthMm < kWidthMinMm)
         {
            const double lim = widthMm > kWidthMaxMm ? kWidthMaxMm : kWidthMinMm;
            std::ostringstream n;
            n << where << "width_mm " << widthMm << " is outside " << kWidthMinMm << " .. " << kWidthMaxMm << ": using " << lim << " (physics ini)";
            Notice(n.str());
            widthMm = lim;
         }
         // the math slot's outline (slot-hole design 1): straight sides t = +-a; round end w^2 = a^2 - u^2 to U1 = L - sqrt(L^2 - a^2); the ramp w^2 = 2 U1 (L - u) to the
         // tip at u = L (tangent at U1). Ramp vertices on the curve at its direction's even steps (at most 15 deg) from u = a, where the floor starts, to the tip. The through part is the
         // straight part, the round ends and the ramp to u = a; the ramp floor from u = a (depth 1.75 mm on the axis) to the tip (0), each cross-section a V through the axis point
         // (the arc's two chords: the fewest faces; the ball meets the floor only near the tip, where both reach the top)
         const double a = 0.5 * widthMm * 1e-3, hs = 0.5 * k.lengthMm * 1e-3, L = kRampM, U1 = L - std::sqrt(L * L - a * a), cc = std::sqrt(2.0 * U1);
         const auto Wu = [&](double u) { return u <= U1 ? std::sqrt(std::max(0.0, a * a - u * u)) : cc * std::sqrt(std::max(0.0, L - u)); };
         std::vector<double> us { 0.0, U1, a };   // the end's vertices by u past the straight part's end (then the ramp's, then the tip)
         const double psiA = std::atan2(cc, 2.0 * std::sqrt(L - a));   // the ramp's direction from the axis at u = a (dw/du = -cc / (2 tau))
         const int nr = std::max(1, (int)std::ceil((0.5 * kPi - psiA) / kTurnRad));
         for (int j = 1; j < nr; j++)
         {
            const double psi = psiA + (0.5 * kPi - psiA) * j / nr, tau = cc / (2.0 * std::tan(psi));
            us.push_back(L - tau * tau);
         }
         std::vector<std::pair<P2, char>> st;   // (s, t) and whether the edge from it is floored, counter-clockwise from (-hs, -a)
         for (const double e : { 1.0, -1.0 })
         {
            for (size_t j = 0; j < us.size(); j++)   // from the side t = -e a round the end at s = e hs
               st.push_back({ { e * (hs + us[j]), -e * Wu(us[j]) }, (char)(us[j] >= a) });
            st.push_back({ { e * (hs + L), 0.0 }, 1 });   // the tip
            for (size_t j = us.size(); j-- > 0;)
               st.push_back({ { e * (hs + us[j]), e * Wu(us[j]) }, (char)(j >= 2 && us[j - 1] >= a) });
         }
         for (const auto &[q, f] : st)
         {
            o.v.push_back(W(q.x, q.y));
            o.floored.push_back(f);
         }
         for (const double e : { 1.0, -1.0 })   // the round ends' circles (vertices on them: the drawing's true curve; the ramp stays as cut)
         {
            const P2 ce = W(e * hs, 0.0);
            g_arcs.push_back({ ce.x, ce.y, a, a });
         }
         for (const double e : { 1.0, -1.0 })   // the ramp floors
         {
            std::vector<double> fu(us.begin() + 2, us.end());
            const auto D = [&](double u) { return kRampDepthM * (L - u) / (L - a); };
            for (size_t j = 0; j + 1 < fu.size(); j++)
            {
               const double u0 = fu[j], u1 = fu[j + 1];
               for (const double side : { 1.0, -1.0 })
               {
                  const P2 o0 = W(e * (hs + u0), side * Wu(u0)), o1 = W(e * (hs + u1), side * Wu(u1)), c0 = W(e * (hs + u0), 0.0), c1 = W(e * (hs + u1), 0.0);
                  o.floor.push_back({ { o0, o1, c1 }, { 0.0, 0.0, D(u1) } });
                  o.floor.push_back({ { o0, c1, c0 }, { 0.0, D(u1), D(u0) } });
               }
            }
            const double ul = fu.back();
            const P2 tip = W(e * (hs + L), 0.0), cl = W(e * (hs + ul), 0.0);
            for (const double side : { 1.0, -1.0 })
               o.floor.push_back({ { W(e * (hs + ul), side * Wu(ul)), tip, cl }, { 0.0, 0.0, D(ul) } });
            const P2 c0 = W(e * (hs + a), 0.0);   // the floor's edge over the through part, at u = a
            for (const double side : { 1.0, -1.0 })
               o.openFloor.push_back({ std::make_pair(W(e * (hs + a), side * Wu(a)), 0.0), std::make_pair(c0, kRampDepthM) });
         }
         log << "slot length " << k.lengthMm << " mm, width " << widthMm << " mm (round end to u " << U1 * 1e3 << " mm, ramp to " << L * 1e3 << " mm, floor from u " << a * 1e3 << " mm, " << nr << " ramp edges past it)";
      }
      else if (k.shape == "groove")
      {
         double lenMm = k.lengthMm;
         if (lenMm < kGrooveLenMinMm || lenMm > kGrooveLenMaxMm)
         {
            const double lim = lenMm < kGrooveLenMinMm ? kGrooveLenMinMm : kGrooveLenMaxMm;
            std::ostringstream n;
            n << where << "length_mm " << lenMm << " is outside " << kGrooveLenMinMm << " .. " << kGrooveLenMaxMm << ": using " << lim << " (physics ini)";
            Notice(n.str());
            lenMm = lim;
         }
         // x, y the plunger end on the axis, s from it up the table: angle 0 runs toward -y (the groove's own frame, Wg: the shapes' frame turned half a turn, the same handedness).
         // Cross-section: kGrooveStrips straight strips tangent to the half-round (paper section 4), creases j = 1 .. N-1 at the arc's angles -phi0 + j delta (radius Rg / cos(delta / 2)),
         // the rims j = 0, N where the outer strips meet the top. Full depth at s = 0 (the square end, a vertical face); the cross-section rises uniformly from there, its depth falling
         // linearly to 0 at s = length (user 2026-11-17: no run-out; every strip flat, crease j reaching the top at s = length dj / dmax, so the outline narrows to a point at the top end)
         const auto Wg = [&](double s, double t) { return W(-s, -t); };
         const double Rg = kGrooveRM, zc = Rg - kGrooveDepthM, phi0 = std::acos(zc / Rg), dl = 2.0 * phi0 / kGrooveStrips, Rr = Rg / std::cos(0.5 * dl);
         const int N = kGrooveStrips;
         std::vector<double> tj(N + 1), dj(N + 1);
         for (int j = 1; j < N; j++)
         {
            const double th = -phi0 + j * dl;
            tj[j] = Rr * std::sin(th);
            dj[j] = Rr * std::cos(th) - zc;
         }
         const double th0 = -phi0 + 0.5 * dl;
         tj[0] = (Rg - zc * std::cos(th0)) / std::sin(th0);
         tj[N] = -tj[0];
         dj[0] = dj[N] = 0.0;
         const double Lg = lenMm * 1e-3, dMax = dj[N / 2];
         std::vector<double> sj(N + 1);
         for (int j = 0; j <= N; j++)
            sj[j] = Lg * dj[j] / dMax;   // where crease j reaches the top (the rims j = 0, N: at s = 0)
         std::vector<std::pair<P2, char>> st { { { 0.0, tj[0] }, 1 } };
         for (int j = 1; j < N; j++)
            st.push_back({ { sj[j], tj[j] }, 1 });
         st.push_back({ { 0.0, tj[N] }, 1 });   // (the end edge back to the start: closed by the end face)
         for (const auto &[q, f] : st)
         {
            o.v.push_back(Wg(q.x, q.y));
            o.floored.push_back(f);
         }
         for (int j = 0; j < N; j++)   // strip j: A at the square end on crease j, B where crease j meets the top, C where crease j + 1 does, D at the square end on crease j + 1 (one plane)
         {
            const P2 A = Wg(0.0, tj[j]), B = Wg(sj[j], tj[j]), C = Wg(sj[j + 1], tj[j + 1]), D = Wg(0.0, tj[j + 1]);
            if (j > 0)   // (B = A at the first strip: no first face)
               o.floor.push_back({ { A, B, C }, { dj[j], 0.0, 0.0 } });
            if (j + 1 < N)   // (C = D at the last strip: no second face)
               o.floor.push_back({ { A, C, D }, { dj[j], 0.0, dj[j + 1] } });
         }
         for (int j = 1; j < N; j++)   // the square end: a vertical fan from the rim at j = 0
            o.side.push_back({ { Wg(0.0, tj[0]), Wg(0.0, tj[j]), Wg(0.0, tj[j + 1]) }, { 0.0, dj[j], dj[j + 1] } });
         log << "groove length " << lenMm << " mm up the table from x, y (angle 0: -y), depth tapering linearly to 0 at the top end, " << N << " strips, width " << 2.0 * tj[N] * 1e3 << " mm at the top, depth "
             << dMax * 1e3 << " mm at the crease";
      }
      else
      {
         Notice(where + (k.shape.empty() ? std::string("no shape") : "shape '" + k.shape + "' is not known") + ": the cut is ignored. The shapes are slot, groove, rectangle, circle and triangle (physics ini)");
         continue;
      }
      PLOGI << "Mesh holes: " << o.name << ": " << log.str() << ", " << o.v.size() << " outline points, " << o.floor.size() << " floor faces" << (o.host.empty() ? "" : ", host " + o.host);
      AddOutline(std::move(o));
   }
}

const std::vector<MeshHoles::DrawnHole> &MeshHoles::Drawn()
{
   return g_drawn;
}

void MeshHoles::Build()
{
   if (g_outlines.empty())
      return;
   // the host under each outline: the one the ini names, else found at its rim points (paper section 9 P2): each outline vertex takes the static mesh host whose up face under it is the
   // highest at or below the part's height plus kTopTolM; all vertices must give the same host. A different host inside the outline (a wall standing in the hole) is neither cut nor a reason to refuse.
   std::map<ContactSurfaceHost *, std::vector<const Outline *>> byHost;
   std::vector<ContactSurfaceHost *> order;
   for (const Outline &o : g_outlines)
   {
      if (o.v.size() < 3)
         continue;
      ContactSurfaceHost *host = nullptr;
      if (!o.host.empty())
      {
         for (ContactHost *const h : ContactHost::AllHosts())
            if (ContactSurfaceHost *const s = dynamic_cast<ContactSurfaceHost *>(h); s != nullptr && s->MeshSurface() != nullptr && IsMeshKind(s->MeshKind()) && PhysicsIni::Lower(s->Name()) == PhysicsIni::Lower(o.host))
               host = s;
         if (host == nullptr)
         {
            Notice("Mesh holes: " + o.name + ": no static mesh named '" + o.host + "' (physics ini host): hole not cut");
            continue;
         }
      }
      else
      {
         std::vector<ContactSurfaceHost *> hostAt(o.v.size(), nullptr);
         size_t missing = 0;
         for (size_t i = 0; i < o.v.size(); i++)
         {
            double bestZ = 0.0;
            for (ContactHost *const h : ContactHost::AllHosts())
            {
               ContactSurfaceHost *const s = dynamic_cast<ContactSurfaceHost *>(h);
               if (s == nullptr || s->MeshSurface() == nullptr || !IsMeshKind(s->MeshKind()))
                  continue;
               for (const RubberContact::SurfaceTri &t : s->MeshSurface()->tris)
               {
                  const T3 q { { t.a.x, t.a.y, t.a.z }, { t.b.x, t.b.y, t.b.z }, { t.c.x, t.c.y, t.c.z } };
                  if (!Up(q) || !InTri(Plan(q), o.v[i]))
                     continue;
                  const double z = PlaneZ(q, o.v[i].x, o.v[i].y);
                  if (z <= o.z + kTopTolM && (hostAt[i] == nullptr || z > bestZ))
                  {
                     hostAt[i] = s;
                     bestZ = z;
                  }
               }
            }
            missing += hostAt[i] == nullptr;
         }
         if (missing > 0)
         {
            Notice("Mesh holes: " + o.name + ": " + std::to_string(missing) + " of its " + std::to_string(o.v.size()) + " rim points have no static mesh top under them at or below its height: hole not cut");
            continue;
         }
         std::vector<ContactSurfaceHost *> hosts;
         for (ContactSurfaceHost *const s : hostAt)
            if (std::find(hosts.begin(), hosts.end(), s) == hosts.end())
               hosts.push_back(s);
         if (hosts.size() > 1)
         {
            std::string names;
            for (const ContactSurfaceHost *const s : hosts)
               names += (names.empty() ? "" : ", ") + s->MeshKind() + " " + s->Name();
            Notice("Mesh holes: " + o.name + " not cut: its rim points lie on different hosts (" + names + ")");
            continue;
         }
         host = hosts[0];
      }
      // a ramp host (named or found at the rim): not cut, in the physics or the drawing (user 2026-11-17: not allowed for now). A ramp's floor is sloped and, curved, not one plane:
      // the drawing's re-cut (playfieldcutdraw.cpp Fill) takes one flat plane, and the ramp's mesh is rebuilt while playing (ramp.cpp GenerateVertexBuffer)
      if (host->MeshKind() == "ramp")
      {
         Notice("Mesh holes: " + o.name + ": host ramp " + host->Name() + ": holes in ramps are not supported yet: hole not cut");
         continue;
      }
      if (byHost.find(host) == byHost.end())
         order.push_back(host);
      byHost[host].push_back(&o);
   }
   for (ContactSurfaceHost *const host : order)
   {
      const std::vector<const Outline *> &list = byHost[host];
      RubberSurface &s = *host->MeshSurface();
      const HostCheck before = Measure(s);
      const size_t facesBefore = s.tris.size();
      Cut c;
      for (const RubberContact::SurfaceTri &t : s.tris)
         c.f.push_back({ { t.a.x, t.a.y, t.a.z }, { t.b.x, t.b.y, t.b.z }, { t.c.x, t.c.y, t.c.z } });
      const double hostThickM = host->ThicknessMm() * 1e-3;   // the walls' depth unless an outline's wall_depth_mm says otherwise
      // pass 2: an outline wholly inside another's floor (a slot in the groove) is cut after it, into that floor's faces (paper Q5, lead 2026-11-16). The others: joined where their
      // plans overlap (paper Q5: one cutout, no warning)
      const size_t n = list.size();
      std::vector<char> second(n, 0);
      for (size_t i = 0; i < n; i++)
         for (size_t j = 0; j < n && !second[i]; j++)
            if (j != i && !list[j]->floor.empty())
            {
               bool all = true;
               for (const P2 &p : list[i]->v)
                  all &= InPoly(list[j]->v, p);
               second[i] = all;
            }
      const auto overlap = [](const Outline &a, const Outline &b) {
         for (const P2 &p : a.v)
            if (InPoly(b.v, p))
               return true;
         for (const P2 &p : b.v)
            if (InPoly(a.v, p))
               return true;
         for (size_t i = 0; i < a.v.size(); i++)
            for (size_t j = 0; j < b.v.size(); j++)
               if (SegX(a.v[i], a.v[(i + 1) % a.v.size()], b.v[j], b.v[(j + 1) % b.v.size()]))
                  return true;
         return false;
      };
      std::vector<size_t> grp(n);
      for (size_t i = 0; i < n; i++)
         grp[i] = i;
      const auto root = [&](size_t i) { while (grp[i] != i) i = grp[i]; return i; };
      for (size_t i = 0; i < n; i++)
         for (size_t j = i + 1; j < n; j++)
            if (!second[i] && !second[j] && overlap(*list[i], *list[j]))
               grp[std::max(root(i), root(j))] = std::min(root(i), root(j));
      std::vector<Group> groups;
      for (int pass = 0; pass < 2; pass++)
         for (size_t i = 0; i < n; i++)
            if (second[i] == pass && (pass == 1 || root(i) == i))
            {
               Group g;
               for (size_t j = 0; j < n; j++)
                  if (pass == 1 ? j == i : (!second[j] && root(j) == i))
                     g.push_back(list[j]);
               groups.push_back(g);
            }
      int cut = 0, nUnder = 0;
      double expectTop = before.top;
      std::vector<std::pair<P3, P3>> bottoms, rims;   // the walls' bottom edges and the open rims (open by design)
      std::ostringstream per;
      std::vector<std::pair<size_t, const Group *>> made;   // each drawn hole of this host (its g_drawn index) and its group: the floors as finally cut (below)
      for (const Group &g : groups)
      {
         std::string names;
         for (const Outline *o : g)
            names += (names.empty() ? "" : " + ") + o->name;
         // the top's height (p2's rule): the plane of the highest up face of the host under the first rim point at or below the part's height (after the cuts before it: a floor for
         // pass 2), taken at the group's centre (CutGroup's top test compares each face's plane there)
         P2 gc { 0.0, 0.0 };
         size_t gn = 0;
         for (const Outline *o : g)
            for (const P2 &q : o->v)
            {
               gc = { gc.x + q.x, gc.y + q.y };
               gn++;
            }
         gc = { gc.x / (double)gn, gc.y / (double)gn };
         double topZ = 0.0, best = 0.0;
         bool found = false;
         for (const T3 &t : c.f)
            if (Up(t) && InTri(Plan(t), g[0]->v[0]))
               if (const double z = PlaneZ(t, g[0]->v[0].x, g[0]->v[0].y); z <= g[0]->z + kTopTolM && (!found || z > best))
               {
                  best = z;
                  topZ = PlaneZ(t, gc.x, gc.y);
                  found = true;
               }
         double thickM = 0.0;   // the group's wall depth: the deepest of its outlines' (wall_depth_mm, else the host's thickness_mm)
         bool draw = false;     // drawn unless every joined outline says draw = 0
         for (const Outline *o : g)
         {
            thickM = std::max(thickM, o->keys.wallDepthMm > 0.0 ? o->keys.wallDepthMm * 1e-3 : hostThickM);
            draw = draw || o->keys.draw;
         }
         GroupResult r;
         MeshHoles::DrawnHole dh;
         dh.draw = draw;
         const size_t facesAt = c.f.size();
         const std::string why = found ? CutGroup(c, g, topZ, thickM, r, bottoms, rims, dh) : std::string("no top face under its first rim point");
         if (!why.empty())
         {
            Notice("Mesh holes: " + names + " not cut in the " + host->MeshKind() + " host: " + why);   // (not drawn either: drawn hole = physics hole)
            continue;
         }
         dh.name = names;
         dh.unit = U;
         dh.hostKind = host->MeshKind();
         dh.hostName = host->Name();
         dh.topZ = topZ / U;
         dh.thick = thickM / U;
         dh.physWallTris = 2 * r.nWalls;
         dh.onTop = !(g.size() == 1 && second[std::find(list.begin(), list.end(), g[0]) - list.begin()]);
         made.push_back({ g_drawn.size(), &g });
         g_drawn.push_back(std::move(dh));
         cut += (int)g.size();
         expectTop += -(r.holeArea - r.missing) - r.dropUp + r.floorUp - r.underCut;
         nUnder += r.nUnder;   // the top inside the joined outline (less an existing hole's opening), the deleted faces that face up, the floors
         per << " " << names << " (" << g.size() << " joined, wall depth " << thickM * 1e3 << " mm, draw " << (draw ? 1 : 0) << ", N";
         for (const Outline *o : g)
            per << " " << o->v.size();
         per << ", inscribed r " << r.inscribed * 1e3 << " mm -> " << (r.nWalls > 0 ? "walls " + std::to_string(r.nWalls) : std::string("no walls")) << ", floor faces " << r.nFloor << ", faces " << facesAt << " -> " << c.f.size()
             << ", grid origin +" << r.gridOff * 1e3 << " mm, top z " << topZ * 1e3 << " mm, underside faces cut " << r.nUnder << " (" << r.underCut * 1e6 << " mm2)";
         if (r.nDrop > 0 || r.nDeep > 0 || r.missing > 1e-6 * r.holeArea)
            per << ", over an existing hole: " << r.missing * 1e6 << " mm2 with no top, " << r.nDrop << " faces inside deleted (" << r.dropUp * 1e6 << " mm2 facing up), " << r.nDeep << " inside below the wall's depth kept";
         per << ");";
      }
      // the drawing's floors as Build finally leaves them (brief 6 item 2): a hole whose floor a later outline was cut into (a slot in the groove) takes, instead of the floor faces
      // its own cut emitted, the host's up faces now below its top inside its outline and no deeper than its deepest floor point, less those inside the inner outline (that hole's own
      // floors, drawn with it): the inner hole's opening is in the drawn floor
      for (const auto &[di, gi] : made)
      {
         MeshHoles::DrawnHole &dh = g_drawn[di];
         if (!dh.onTop)
            continue;
         std::vector<const Outline *> inner;
         for (const auto &[dj, gj] : made)
         {
            const Outline *const in = (*gj)[0];
            if (!g_drawn[dj].onTop && std::any_of(gi->begin(), gi->end(), [&](const Outline *o) { return InPoly(o->v, in->v[0]); }))
               inner.push_back(in);
         }
         if (inner.empty())
            continue;
         double floorMax = 0.0;
         for (const Outline *o : *gi)
            for (const FTri &t : o->floor)
               floorMax = std::max({ floorMax, t.d[0], t.d[1], t.d[2] });
         std::vector<float> wood;
         for (size_t i = 0; i + 18 <= dh.wood.size(); i += 18)
            if (dh.wood[i + 5] != 0.f)
               wood.insert(wood.end(), dh.wood.begin() + i, dh.wood.begin() + i + 18);
         const double top = dh.topZ * U;
         int nf = 0;
         for (const T3 &t : c.f)
         {
            const P2 m { (t.a.x + t.b.x + t.c.x) / 3.0, (t.a.y + t.b.y + t.c.y) / 3.0 };
            const double dm = top - (t.a.z + t.b.z + t.c.z) / 3.0;
            if (!Up(t) || dm <= 1e-6 || dm > floorMax + kSnapM || !std::any_of(gi->begin(), gi->end(), [&](const Outline *o) { return InPoly(o->v, m); })
                || std::any_of(inner.begin(), inner.end(), [&](const Outline *o) { return InPoly(o->v, m); }))
               continue;
            for (const P3 *p : { &t.a, &t.b, &t.c })
               wood.insert(wood.end(), { (float)(p->x / U), (float)(p->y / U), (float)(p->z / U), (float)((p->x + p->y) / U), (float)((top - p->z) / U), 0.f });
            nf++;
         }
         PLOGI << "Mesh holes: " << dh.name << ": drawn floor as finally cut (an inner hole cut into it): " << nf << " faces";
         dh.wood.swap(wood);
      }
      if (cut == 0)
         continue;
      const int over = Store(host, c.f);   // the cut faces (the uncut list is not kept), welded and linked again; edges on more than two faces: strays (overlapping faces)
      const HostCheck after = Measure(s);
      // every open edge after the cut is a wall's bottom, an open rim by design, or a piece of an open edge the mesh had (its border, split where a grid line crosses it); any other is a gap
      int stray = 0;
      for (const auto &[a, b] : after.edges)
      {
         bool ok = false;
         for (const auto &[p, q] : bottoms)
            ok = ok || (OnSeg(a, p, q) && OnSeg(b, p, q));
         for (const auto &[p, q] : rims)
            ok = ok || (OnSeg(a, p, q) && OnSeg(b, p, q));
         for (const auto &[p, q] : before.edges)
            ok = ok || (OnSeg(a, p, q) && OnSeg(b, p, q));
         stray += !ok;
         if (!ok && stray <= 20)   // where they are (brief 4: run a2's 13 strays were not found offline on a clean box)
            PLOGW << "Mesh holes: stray open edge " << stray << " on the " << host->MeshKind() << " host: (" << a.x / U << ", " << a.y / U << ", " << a.z / U << ") -> (" << b.x / U << ", "
                  << b.y / U << ", " << b.z / U << ") VPU";
      }
      PLOGI << "Mesh holes: " << host->MeshKind() << " host " << host->Name() << ": faces " << facesBefore << " -> " << s.tris.size() << ", holes " << cut << ":" << per.str() << " open edges " << before.open << " -> " << after.open
            << " (" << bottoms.size() << " wall bottoms, " << rims.size() << " open rim edges, " << stray << " stray, " << over << " edges on more than two faces), underside faces cut " << nUnder << ", top area " << before.top << " -> " << after.top << " m2 (expected " << expectTop << ", diff " << (after.top - expectTop) << ")";
      // The notice is for the edge checks only: a stray open edge or an edge on more than two faces is a real gap or overlap in the surface the ball rolls on. The top area (measured
      // from the faces as stored, expected from the cut) is logged in the line above, not a notice: it differs by float32 rounding, not by a gap (brief 5, turn table test 4: the cut's
      // faces in doubles give each level's area less the outline's exactly; stored as float32 they give +0.0212 mm2 per level, and the expected's float32 rim points 8.3e-4 mm2 less)
      if (stray > 0 || over > 0)
         Notice("Mesh holes: load check failed on the " + host->MeshKind() + " host: " + std::to_string(stray) + " open edges that are neither a wall's bottom, an open rim, nor the mesh's own border, "
                + std::to_string(over) + " edges on more than two faces");
   }
}
