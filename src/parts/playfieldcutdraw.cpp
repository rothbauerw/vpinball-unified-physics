// license:GPLv3+

// The playfield's cuts as drawn (design runs/2026-11-09-slot-hole/design.md section 8). The shape is the physics' (physics/playfieldcut.h, PlayfieldCuts::Drawn, the same numbers as the rim):
// the sharp outline (the rim's 0.25 mm rounding is not drawn), the through-hole's walls straight down the playfield's thickness, the ramp floors. All made once with the playfield's render
// buffers (Primitive::RenderSetup): the hole is cut in the playfield's own mesh, so every pass that draws the playfield (static prepass, dynamic pass, VR) draws it; the walls and floors are
// one mesh drawn right after the playfield (Primitive::Render) with its material and the wood texture assets/PlayfieldCutWood.png. Per frame: that one draw call.

#include "core/stdafx.h"
#include "parts/playfieldcutdraw.h"
#include "parts/primitive.h"
#include "physics/playfieldcut.h"
#include "core/VPApp.h"
#include "renderer/IndexBuffer.h"
#include "renderer/MeshBuffer.h"
#include "renderer/RenderDevice.h"
#include "renderer/Renderer.h"
#include "renderer/Shader.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <sstream>

namespace
{
   using PlayfieldCuts::DrawnSlot;

   // Curves are drawn as chords within kDrawTolM inside the true curve. kDrawTolM is a tenth of a pixel where the playfield is seen finest: in VR, the eye 30 cm from the playfield at about
   // 25 pixels per degree, a pixel is 0.21 mm; on a 4K cabinet screen (2160 pixels across a 20.25 in playfield) 0.24 mm. A tenth of a pixel is finer than antialiasing resolves, so the
   // chords look like the curve. On the 5.03 mm slot: 1 chord per round end of the outline (its 9.9 deg arc of 2.515 mm radius strays 0.009 mm), 19 per ramp (round the tip, whose
   // 0.43 mm radius is the outline's tightest curve).
   constexpr double kDrawTolM = 0.02e-3;
   constexpr double kTopTolM = 1e-6;   // a face is at the playfield's top over a cut when its corners are within this of it (the physics' flatness test, PlayfieldCuts::Setup)

   struct P2
   {
      double x, y;
   };
   double Cross(const P2 &o, const P2 &a, const P2 &b) { return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x); }

   // chords: a circle's over the angle dphi strays r (1 - cos(dphi / 2)) from it; the ramp's over dtau at most dtau^2 |P''| / 8 = dtau^2 / 4 (P'' = (-2, 0) in the slot frame)
   double ChordAngle(const double r) { return 2.0 * std::acos(1.0 - kDrawTolM / r); }
   int RampChords(const DrawnSlot &d) { return std::max(2, (int)std::ceil(d.tauJ / std::sqrt(kDrawTolM))); }

   // a point of cut d's frame (s along, t across, metres) in the world (VPX units)
   P2 World(const DrawnSlot &d, const double s, const double t) { return { (d.cx + d.lx * s + d.px * t) / d.unit, (d.cy + d.ly * s + d.py * t) / d.unit }; }

   // cut d's sharp outline in the world (VPX units), counter-clockwise in x, y: the straight sides, the round ends (radius a, to U1), the ramps (u = L - tau^2, w = c tau)
   std::vector<P2> Outline(const DrawnSlot &d)
   {
      const int nRound = std::max(1, (int)std::ceil(d.th1 / ChordAngle(d.a)));
      const int nRamp = RampChords(d);
      std::vector<P2> st;   // slot frame (s, t), metres
      for (const double e : { 1.0, -1.0 })   // the side t = -e a along +e s, then the end at +e: from the side to U1, round the ramp's tip, back to the other side
      {
         st.push_back({ -e * d.hs, -e * d.a });
         for (int k = 0; k <= nRound; k++)
         {
            const double ph = d.th1 * k / nRound;
            st.push_back({ e * (d.hs + d.a * std::sin(ph)), -e * d.a * std::cos(ph) });
         }
         for (int k = 1; k < nRamp; k++)
         {
            const double tau = d.tauJ * (2.0 * k / nRamp - 1.0);
            st.push_back({ e * (d.hs + d.L - tau * tau), e * d.c * tau });
         }
         for (int k = nRound; k >= 0; k--)
         {
            const double ph = d.th1 * k / nRound;
            st.push_back({ e * (d.hs + d.a * std::sin(ph)), e * d.a * std::cos(ph) });
         }
      }
      std::vector<P2> out;   // (the points shared by two parts once)
      for (const P2 &q : st)
      {
         const P2 w = World(d, q.x, q.y);
         if (out.empty() || std::hypot(w.x - out.back().x, w.y - out.back().y) > 1e-9)
            out.push_back(w);
      }
      while (out.size() > 1 && std::hypot(out[0].x - out.back().x, out[0].y - out.back().y) <= 1e-9)
         out.pop_back();
      return out;
   }

   // two convex polygons, counter-clockwise, overlap: no side of either has the other wholly on its outer side (touching is not overlapping)
   bool Overlap(const P2 *a, const size_t na, const P2 *b, const size_t nb)
   {
      const auto apart = [](const P2 *p, const size_t np, const P2 *q, const size_t nq) {
         for (size_t i = 0; i < np; i++)
            if (std::all_of(q, q + nq, [&](const P2 &v) { return Cross(p[i], p[(i + 1) % np], v) <= 0.0; }))
               return true;
         return false;
      };
      return !apart(a, na, b, nb) && !apart(b, nb, a, na);
   }

   // the segments p1 p2 and q1 q2 cross, or a point of one lies on the other (the caller leaves out the edges sharing an end with p1 p2)
   bool Meet(const P2 &p1, const P2 &p2, const P2 &q1, const P2 &q2)
   {
      const double d1 = Cross(q1, q2, p1), d2 = Cross(q1, q2, p2), d3 = Cross(p1, p2, q1), d4 = Cross(p1, p2, q2);
      if (((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0)) && ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0)))
         return true;
      const auto on = [](const P2 &a, const P2 &b, const P2 &v, const double d) {
         return d == 0.0 && std::min(a.x, b.x) <= v.x && v.x <= std::max(a.x, b.x) && std::min(a.y, b.y) <= v.y && v.y <= std::max(a.y, b.y);
      };
      return on(q1, q2, p1, d1) || on(q1, q2, p2, d2) || on(p1, p2, q1, d3) || on(p1, p2, q2, d4);
   }

   double Area(const std::vector<P2> &pt, const std::vector<int> &L)   // twice the signed area, counter-clockwise positive
   {
      double a = 0.0;
      for (size_t i = 0; i < L.size(); i++)
         a += Cross({ 0.0, 0.0 }, pt[L[i]], pt[L[(i + 1) % L.size()]]);
      return a;
   }

   // The playfield's faces sel (its top over the cuts) as one region less the outlines holes (world, counter-clockwise, convex, apart, inside the region). The region's boundary is kept (its
   // edges are the neighbour faces' own: no cracks); its inside is filled again by VPX's ear clipping (PolygonToTriangles) of one simple polygon per boundary loop, each hole joined to it
   // first by the shortest bridge that meets no edge. The other faces are copied to idx, the new triangles appended, wound as the faces were. A hole's point is a new vertex in verts with
   // the values (position, normal, texture coordinates) of the faces' plane: they must be one plane (an affine map of x, y), as on the implicit playfield quad. False (why): not covered.
   bool Fill(const Mesh &mesh, const std::vector<P2> &wp, const std::vector<char> &sel, const std::vector<std::vector<P2>> &holes, std::vector<Vertex3D_NoTex2> &verts,
             std::vector<unsigned int> &idx, std::string &why)
   {
      const size_t nTri = mesh.m_indices.size() / 3;
      std::map<std::pair<double, double>, int> placeOf;   // the region's corners by place (a mesh may have several vertices at one place)
      std::vector<P2> pt;                                 // the corners' places, then the holes' points
      std::vector<unsigned int> vtx;                      // their vertex in verts
      std::set<std::pair<int, int>> edges;                // the faces' edges, counter-clockwise in x, y
      int wind = 0, big = -1;                             // the faces' winding (+1 counter-clockwise in x, y); the largest face (the plane's values from its corners)
      double bigArea = 0.0;
      for (size_t t = 0; t < nTri; t++)
      {
         const unsigned int *I = &mesh.m_indices[3 * t];
         if (!sel[t])
         {
            idx.insert(idx.end(), I, I + 3);
            continue;
         }
         const double ar = Cross(wp[I[0]], wp[I[1]], wp[I[2]]);
         if (wind != 0 && (ar > 0.0 ? 1 : -1) != wind)
         {
            why = "the faces over the cuts are wound both ways";
            return false;
         }
         wind = ar > 0.0 ? 1 : -1;
         if (std::abs(ar) > bigArea)
         {
            bigArea = std::abs(ar);
            big = (int)t;
         }
         int p[3];
         for (int k = 0; k < 3; k++)
         {
            const auto [it, fresh] = placeOf.try_emplace({ wp[I[k]].x, wp[I[k]].y }, (int)pt.size());
            if (fresh)
            {
               pt.push_back(wp[I[k]]);
               vtx.push_back(I[k]);
            }
            p[k] = it->second;
         }
         if (wind < 0)
            std::swap(p[1], p[2]);
         for (int k = 0; k < 3; k++)
            edges.insert({ p[k], p[(k + 1) % 3] });
      }
      if (big < 0)
      {
         why = "no face at the top over the cuts";
         return false;
      }
      // the plane's values at a place: from the largest face's corners by their weights; every corner of the region must be on it
      const unsigned int *IB = &mesh.m_indices[3 * big];
      const Vertex3D_NoTex2 V[3] = { mesh.m_vertices[IB[0]], mesh.m_vertices[IB[1]], mesh.m_vertices[IB[2]] };
      const P2 K[3] = { wp[IB[0]], wp[IB[1]], wp[IB[2]] };
      const auto plane = [&](const P2 &q) {
         const double A = Cross(K[0], K[1], K[2]), w0 = Cross(q, K[1], K[2]) / A, w1 = Cross(K[0], q, K[2]) / A, w[3] = { w0, w1, 1.0 - w0 - w1 };
         const auto mix = [&](float Vertex3D_NoTex2::*f) { return (float)(w[0] * (double)(V[0].*f) + w[1] * (double)(V[1].*f) + w[2] * (double)(V[2].*f)); };
         Vertex3D_NoTex2 nv;
         nv.x = mix(&Vertex3D_NoTex2::x);
         nv.y = mix(&Vertex3D_NoTex2::y);
         nv.z = mix(&Vertex3D_NoTex2::z);
         nv.nx = mix(&Vertex3D_NoTex2::nx);
         nv.ny = mix(&Vertex3D_NoTex2::ny);
         nv.nz = mix(&Vertex3D_NoTex2::nz);
         nv.tu = mix(&Vertex3D_NoTex2::tu);
         nv.tv = mix(&Vertex3D_NoTex2::tv);
         return nv;
      };
      for (size_t i = 0; i < pt.size(); i++)
      {
         const Vertex3D_NoTex2 &a = verts[vtx[i]], b = plane(pt[i]);
         if (std::abs(a.tu - b.tu) > 1e-4f || std::abs(a.tv - b.tv) > 1e-4f || std::abs(a.nx - b.nx) > 1e-3f || std::abs(a.ny - b.ny) > 1e-3f || std::abs(a.nz - b.nz) > 1e-3f)
         {
            why = "the faces over the cuts have no one plane of texture coordinates and normals";
            return false;
         }
      }
      // the region's boundary: the edges whose reverse is no edge, and its loops (counter-clockwise: an outer edge; clockwise: a gap in the faces)
      std::map<int, int> next;
      for (const std::pair<int, int> &e : edges)
         if (edges.count({ e.second, e.first }) == 0 && !next.emplace(e.first, e.second).second)
         {
            why = "the faces over the cuts touch at a corner only";
            return false;
         }
      std::vector<std::vector<int>> loops;
      std::set<int> seen;
      for (const auto &[first, unused] : next)
      {
         if (seen.count(first))
            continue;
         std::vector<int> L;
         int s = first;
         do
         {
            if (!seen.insert(s).second || next.count(s) == 0)
            {
               why = "the faces over the cuts have a boundary that does not close";
               return false;
            }
            L.push_back(s);
            s = next[s];
         } while (s != first);
         if (Area(pt, L) <= 0.0)
         {
            why = "the faces over the cuts leave a gap";
            return false;
         }
         loops.push_back(std::move(L));
      }
      // each hole in the loop holding all its points (crossing number); its points clockwise
      std::vector<std::vector<std::vector<int>>> holesOf(loops.size());
      double holeArea = 0.0;
      for (const std::vector<P2> &O : holes)
      {
         const auto inLoop = [&](const std::vector<int> &L, const P2 &q) {
            bool in = false;
            for (size_t k = 0, l = L.size() - 1; k < L.size(); l = k++)
            {
               const P2 &a = pt[L[k]], &b = pt[L[l]];
               if ((a.y > q.y) != (b.y > q.y) && q.x < a.x + (q.y - a.y) * (b.x - a.x) / (b.y - a.y))
                  in = !in;
            }
            return in;
         };
         const auto owner = std::find_if(loops.begin(), loops.end(), [&](const std::vector<int> &L) { return std::all_of(O.begin(), O.end(), [&](const P2 &q) { return inLoop(L, q); }); });
         if (owner == loops.end())
         {
            why = "a cut runs past the edge of the faces round it";
            return false;
         }
         std::vector<int> h;
         for (size_t k = O.size(); k-- > 0;)
         {
            h.push_back((int)pt.size());
            pt.push_back(O[k]);
            Vertex3D_NoTex2 nv = plane(O[k]);
            const float nl = std::sqrt(nv.nx * nv.nx + nv.ny * nv.ny + nv.nz * nv.nz);
            if (nl > 0.f)
            {
               nv.nx /= nl;
               nv.ny /= nl;
               nv.nz /= nl;
            }
            vtx.push_back((unsigned int)verts.size());
            verts.push_back(nv);
         }
         holeArea -= Area(pt, h);
         holesOf[owner - loops.begin()].push_back(std::move(h));
      }
      // VPX's ear clipping works in floats: the points relative to the first hole's first point (0.1 um at a playfield's size)
      std::vector<RenderVertex> rgv(pt.size());
      for (size_t i = 0; i < pt.size(); i++)
      {
         rgv[i].x = (float)(pt[i].x - holes[0][0].x);
         rgv[i].y = (float)(pt[i].y - holes[0][0].y);
      }
      std::vector<char> twice(pt.size(), 0);   // a bridge's ends (met twice in the polygon): not a bridge's end again
      double loopArea = 0.0, triArea = 0.0;
      for (size_t li = 0; li < loops.size(); li++)
      {
         std::vector<int> P = loops[li];
         std::vector<std::vector<int>> &H = holesOf[li];
         loopArea += Area(pt, P);
         while (!H.empty())
         {
            int bh = -1, bi = -1, bj = -1;
            double best = 1e300;
            const auto clear = [&](const int a, const int b, const std::vector<int> &E) {
               for (size_t k = 0; k < E.size(); k++)
               {
                  const int q1 = E[k], q2 = E[(k + 1) % E.size()];
                  if (q1 != a && q1 != b && q2 != a && q2 != b && Meet(pt[a], pt[b], pt[q1], pt[q2]))
                     return false;
               }
               return true;
            };
            for (size_t h = 0; h < H.size(); h++)
               for (size_t i = 0; i < H[h].size(); i++)
                  for (size_t j = 0; j < P.size(); j++)
                  {
                     const int a = H[h][i], b = P[j];
                     const double dd = (pt[a].x - pt[b].x) * (pt[a].x - pt[b].x) + (pt[a].y - pt[b].y) * (pt[a].y - pt[b].y);
                     if (dd >= best || twice[b] || !clear(a, b, P) || !std::all_of(H.begin(), H.end(), [&](const std::vector<int> &E) { return clear(a, b, E); }))
                        continue;
                     best = dd;
                     bh = (int)h;
                     bi = (int)i;
                     bj = (int)j;
                  }
            if (bh < 0)
            {
               why = "no bridge from a cut to the faces' boundary";
               return false;
            }
            const std::vector<int> &Hb = H[bh];
            twice[P[bj]] = twice[Hb[bi]] = 1;
            std::vector<int> merged(P.begin(), P.begin() + bj + 1);   // ..., the polygon's point, round the hole from its point back to it, the polygon's point again, ...
            for (size_t k = 0; k <= Hb.size(); k++)
               merged.push_back(Hb[(bi + k) % Hb.size()]);
            merged.insert(merged.end(), P.begin() + bj, P.end());
            P.swap(merged);
            H.erase(H.begin() + bh);
         }
         // VPX's ear clipping takes the polygon clockwise in x, y (as VPX's counter-clockwise outlines with y down) and gives each triangle (a, c, b): counter-clockwise
         std::vector<unsigned int> poly(P.rbegin(), P.rend());
         std::vector<unsigned int> tri;
         PolygonToTriangles(rgv, poly, tri, false);
         if (tri.size() != 3 * (P.size() - 2))
         {
            why = "the ear clipping left a part unfilled";
            return false;
         }
         for (size_t k = 0; k < tri.size(); k += 3)
         {
            const double a = Cross(pt[tri[k]], pt[tri[k + 1]], pt[tri[k + 2]]);
            triArea += a;
            if (a < 0.0)
               triArea = -1e300;   // turned over: fails the area check
            idx.push_back(vtx[tri[k]]);
            idx.push_back(vtx[tri[k + (wind > 0 ? 1 : 2)]]);
            idx.push_back(vtx[tri[k + (wind > 0 ? 2 : 1)]]);
         }
      }
      if (std::abs(triArea - (loopArea - holeArea)) > 1e-6 * loopArea)
      {
         why = "the new triangles do not cover the faces less the cuts";
         return false;
      }
      return true;
   }

   // Cut d's walls and ramp floors in the wood (world, VPX units; normals into the hole), appended to v and idx. Texture: the grain along the wall at the texture's own scale, its height the
   // wall's depth: u = length along the wall (on the floor: along the slot) / rep, v = depth / thick (the floor shows the layer of wood at its depth, as the wall beside it).
   void Wood(const DrawnSlot &d, const double rep, std::vector<Vertex3D_NoTex2> &v, std::vector<unsigned int> &idx, int &nWall, int &nFloor, int &nSec, int &nAcross)
   {
      const size_t first = v.size();
      std::vector<std::array<double, 6>> at;   // each vertex in the slot frame (s, t, z, normal), from v[first] on: the winding is decided there, in doubles
      const auto vert = [&](const double s, const double t, const double z, const double ns, const double nt, const double nz, const double tu) {
         const P2 w = World(d, s, t);
         const double nx = d.lx * ns + d.px * nt, ny = d.ly * ns + d.py * nt, l = std::sqrt(nx * nx + ny * ny + nz * nz);
         Vertex3D_NoTex2 q;
         q.x = (float)w.x;
         q.y = (float)w.y;
         q.z = (float)((d.z0 + z) / d.unit);
         q.nx = (float)(nx / l);
         q.ny = (float)(ny / l);
         q.nz = (float)(nz / l);
         q.tu = (float)tu;
         q.tv = (float)(-z / d.thick);
         v.push_back(q);
         at.push_back({ s, t, z, ns / l, nt / l, nz / l });
         return (unsigned int)v.size() - 1;
      };
      // a triangle wound to be seen from the side its corners' normals face (the playfield's top is counter-clockwise seen from above); one of no area is left out
      const auto tri = [&](const unsigned int a, const unsigned int b, const unsigned int c) {
         const std::array<double, 6> &A = at[a - first], &B = at[b - first], &C = at[c - first];
         const double e1[3] = { B[0] - A[0], B[1] - A[1], B[2] - A[2] }, e2[3] = { C[0] - A[0], C[1] - A[1], C[2] - A[2] };
         const double n[3] = { e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0] };
         if (n[0] == 0.0 && n[1] == 0.0 && n[2] == 0.0)
            return 0;
         const double f = n[0] * (A[3] + B[3] + C[3]) + n[1] * (A[4] + B[4] + C[4]) + n[2] * (A[5] + B[5] + C[5]);
         idx.push_back(a);
         idx.push_back(f >= 0.0 ? b : c);
         idx.push_back(f >= 0.0 ? c : b);
         return 1;
      };
      // the ramp floor's cross-sections, u past the end of the straight part: its outline's chord ends (the hole's points), the through-hole's round end's chords where it cuts the floor
      // (U1 < u < a), and the tip (an arc of no width; its radius stays finite)
      const int nRamp = RampChords(d);
      const double dPh = ChordAngle(d.a);
      const int nRound = std::max(1, (int)std::ceil(d.th1 / dPh)), nCut = std::max(1, (int)std::ceil((0.5 * M_PI - d.th1) / dPh));
      std::vector<double> us { d.U1, d.a, d.L * (1.0 - 1e-7) };
      for (int k = 1; k < nRamp; k++)
      {
         const double tau = d.tauJ * (2.0 * k / nRamp - 1.0);
         if (tau > 0.0)
            us.push_back(d.L - tau * tau);
      }
      for (int k = 1; k < nCut; k++)
         us.push_back(d.a * std::sin(d.th1 + (0.5 * M_PI - d.th1) * k / nCut));
      std::sort(us.begin(), us.end());
      us.erase(std::unique(us.begin(), us.end(), [](const double x, const double y) { return std::abs(x - y) < 1e-12; }), us.end());
      struct Sec
      {
         double u, rho, zc, drho, dzc, p0, p1;   // the floor's arc and its change along u; from the angle p0 (the through-hole's edge, lower half, or the axis) to p1 (the outline)
      };
      std::vector<Sec> sec;
      nAcross = 1;
      for (const double u : us)
      {
         Sec q;
         q.u = u;
         d.FloorArc(u, q.rho, q.zc);
         const double h = std::min(1e-7, 0.5 * (d.L - u));
         double r1, z1;
         d.FloorArc(u + h, r1, z1);
         q.drho = (r1 - q.rho) / h;
         q.dzc = (z1 - q.zc) / h;
         q.p0 = -std::acos(std::min(1.0, d.C(u) / q.rho));
         q.p1 = std::atan2(-q.zc, d.W(u));
         nAcross = std::max(nAcross, (int)std::ceil((q.p1 - q.p0) / ChordAngle(q.rho)));
         sec.push_back(q);
      }
      nSec = (int)sec.size();
      // the floor point j of section q (u, t >= 0, z) and its normal into the cut: the arc's surface t^2 + (z - zc(u))^2 = rho(u)^2, its gradient turned inward
      const auto floorPt = [&](const Sec &q, const int j, double &t, double &z, double &nu, double &nt, double &nz) {
         if (j == 0)   // on the through-hole's wall (its top edge) or the axis: the same numbers as the wall's
         {
            t = d.C(q.u);
            z = q.zc - std::sqrt(std::max(0.0, q.rho * q.rho - t * t));
         }
         else if (j == nAcross)   // on the outline, at the top
         {
            t = d.W(q.u);
            z = 0.0;
         }
         else
         {
            const double ps = q.p0 + (q.p1 - q.p0) * j / nAcross;
            t = q.rho * std::cos(ps);
            z = q.zc + q.rho * std::sin(ps);
         }
         nu = (z - q.zc) * q.dzc + q.rho * q.drho;
         nt = -t;
         nz = q.zc - z;
      };
      // the walls: the through-hole's outline (the straight sides, the round ends of radius a), counter-clockwise; its top edge on the top up to U1, then on the floor
      struct WP
      {
         double s, t, z, ns, nt;
      };
      std::vector<WP> loop;
      for (const double e : { 1.0, -1.0 })
      {
         loop.push_back({ -e * d.hs, -e * d.a, 0.0, 0.0, e });   // the side t = -e a, from s = -e hs
         std::vector<WP> half;   // the round end at +e from that side (phi = 0) to the axis (phi = pi / 2)
         for (int k = 0; k <= nRound; k++)
         {
            const double ph = d.th1 * k / nRound;
            half.push_back({ d.a * std::sin(ph), d.a * std::cos(ph), 0.0, -std::sin(ph), std::cos(ph) });
         }
         for (const Sec &q : sec)
            if (q.u > d.U1 && q.u <= d.a)
            {
               double t, z, nu, nt, nz;
               floorPt(q, 0, t, z, nu, nt, nz);
               const double ph = std::asin(std::min(1.0, q.u / d.a));
               half.push_back({ q.u, t, z, -std::sin(ph), std::cos(ph) });
            }
         for (size_t k = 0; k < half.size(); k++)
            loop.push_back({ e * (d.hs + half[k].s), -e * half[k].t, half[k].z, e * half[k].ns, e * half[k].nt });
         for (size_t k = half.size() - 1; k-- > 0;)   // on round the other side (the axis once)
            loop.push_back({ e * (d.hs + half[k].s), e * half[k].t, half[k].z, e * half[k].ns, -e * half[k].nt });
      }
      const size_t base = v.size();   // (the loop ends where it began, its texture's seam)
      double len = 0.0;
      for (size_t k = 0; k < loop.size(); k++)
      {
         if (k > 0)
            len += std::hypot(loop[k].s - loop[k - 1].s, loop[k].t - loop[k - 1].t);
         vert(loop[k].s, loop[k].t, loop[k].z, loop[k].ns, loop[k].nt, 0.0, len / rep);
         vert(loop[k].s, loop[k].t, -d.thick, loop[k].ns, loop[k].nt, 0.0, len / rep);
      }
      nWall = 0;
      for (size_t k = 0; k + 1 < loop.size(); k++)
      {
         const unsigned int t0 = (unsigned int)(base + 2 * k), b0 = t0 + 1, t1 = t0 + 2, b1 = t0 + 3;
         nWall += tri(t0, b0, b1);
         nWall += tri(t0, b1, t1);
      }
      // the floors: at each end e, each side of the axis g, a grid of sections by points across
      nFloor = 0;
      for (const double e : { 1.0, -1.0 })
         for (const double g : { 1.0, -1.0 })
         {
            const size_t b = v.size();
            for (const Sec &q : sec)
               for (int j = 0; j <= nAcross; j++)
               {
                  double t, z, nu, nt, nz;
                  floorPt(q, j, t, z, nu, nt, nz);
                  const double s = e * (d.hs + q.u);
                  vert(s, g * t, z, e * nu, g * nt, nz, s / rep);
               }
            for (size_t i = 0; i + 1 < sec.size(); i++)
               for (int j = 0; j < nAcross; j++)
               {
                  const unsigned int p = (unsigned int)(b + i * (nAcross + 1) + j), q = p + (unsigned int)(nAcross + 1);
                  nFloor += tri(p, q, q + 1);
                  nFloor += tri(p, q + 1, p + 1);
               }
         }
   }

   struct Drawing
   {
      const Primitive *prim = nullptr;
      std::vector<Vertex3D_NoTex2> vertices;   // the table's own mesh while the cut one is drawn
      std::vector<unsigned int> indices;
      std::shared_ptr<MeshBuffer> wood;        // the walls and floors (world, VPX units)
      unsigned int woodIndices = 0;
      std::unique_ptr<Texture> tex;
   };
   Drawing g_drawing;
}

namespace PlayfieldCutDraw
{
   void Setup(Primitive *const pf, Renderer *const renderer)
   {
      const std::vector<DrawnSlot> cuts = PlayfieldCuts::Drawn(pf->GetPTable());
      if (cuts.empty() || g_drawing.prim != nullptr)
         return;
      Mesh &mesh = pf->m_mesh;
      const std::string where = "Playfield cuts, drawn playfield " + pf->GetName() + ": ";
      if (!mesh.m_animationFrames.empty())
      {
         PLOGW << where << "an animated mesh (its frames rewrite the vertices): drawn uncut";
         return;
      }
      // the corners in the world (VPX units)
      const Matrix3D &M = pf->RecalculateMatrices();
      const size_t nTri = mesh.m_indices.size() / 3;
      std::vector<P2> wp(mesh.m_vertices.size());
      std::vector<double> wz(mesh.m_vertices.size());
      for (size_t i = 0; i < wp.size(); i++)
      {
         const Vertex3D_NoTex2 &v = mesh.m_vertices[i];
         const Vertex3Ds w = M.MultiplyVectorNoPerspective(Vertex3Ds(v.x, v.y, v.z));
         wp[i] = { w.x, w.y };
         wz[i] = w.z;
      }
      // the outlines (a cut overlapping an earlier one is not drawn), and the faces over each: at the top (re-cut), under it (left: the playfield's underside), or neither (the drawn
      // playfield is not flat at the cut: the cut is not drawn, design section 6)
      std::vector<std::vector<P2>> outl;
      std::vector<const DrawnSlot *> drawn;
      std::vector<char> sel(nTri, 0);
      std::ostringstream names;
      int under = 0, underUp = 0;
      for (const DrawnSlot &d : cuts)
      {
         std::vector<P2> o = Outline(d);
         if (std::any_of(outl.begin(), outl.end(), [&](const std::vector<P2> &q) { return Overlap(o.data(), o.size(), q.data(), q.size()); }))
         {
            PLOGW << where << "cut " << d.name << " overlaps another cut: not drawn";
            continue;
         }
         P2 lo { 1e300, 1e300 }, hi { -1e300, -1e300 };
         for (const P2 &q : o)
         {
            lo = { std::min(lo.x, q.x), std::min(lo.y, q.y) };
            hi = { std::max(hi.x, q.x), std::max(hi.y, q.y) };
         }
         const double z0 = d.z0 / d.unit, tol = kTopTolM / d.unit;
         std::vector<int> top;
         int nUnder = 0, nUnderUp = 0;
         bool flat = true;
         for (size_t t = 0; t < nTri && flat; t++)
         {
            const unsigned int *I = &mesh.m_indices[3 * t];
            const P2 &A = wp[I[0]], &B = wp[I[1]], &C = wp[I[2]];
            if (std::max({ A.x, B.x, C.x }) <= lo.x || std::min({ A.x, B.x, C.x }) >= hi.x || std::max({ A.y, B.y, C.y }) <= lo.y || std::min({ A.y, B.y, C.y }) >= hi.y)
               continue;
            const double ar = Cross(A, B, C);
            if (ar == 0.0)   // seen edge-on from above: covers nothing
               continue;
            const P2 K[3] = { A, ar > 0.0 ? B : C, ar > 0.0 ? C : B };
            if (!Overlap(K, 3, o.data(), o.size()))
               continue;
            const double zHi = std::max({ wz[I[0]], wz[I[1]], wz[I[2]] }), zLo = std::min({ wz[I[0]], wz[I[1]], wz[I[2]] });
            if (zHi < z0 - tol)
            {
               nUnder++;
               nUnderUp += ar > 0.0;   // (wound as the playfield's top: drawn when seen from above)
            }
            else if (zLo >= z0 - tol && zHi <= z0 + tol)
               top.push_back((int)t);
            else
               flat = false;
         }
         if (!flat)
         {
            PLOGW << where << "cut " << d.name << ": the drawn playfield is not flat at its top over the cut: not drawn (design section 6)";
            continue;
         }
         for (const int t : top)
            sel[t] = 1;
         under += nUnder;
         underUp += nUnderUp;
         names << (outl.empty() ? "" : ", ") << d.name << " (" << o.size() << " outline points, over " << top.size() << " face(s))";
         outl.push_back(std::move(o));
         drawn.push_back(&d);
      }
      if (outl.empty())
         return;
      std::vector<Vertex3D_NoTex2> verts = mesh.m_vertices;
      std::vector<unsigned int> idx;
      std::string why;
      if (!Fill(mesh, wp, sel, outl, verts, idx, why))
      {
         PLOGW << where << why << ": the playfield is drawn uncut";
         return;
      }
      g_drawing.prim = pf;
      g_drawing.vertices = mesh.m_vertices;
      g_drawing.indices = mesh.m_indices;
      const size_t nV = mesh.m_vertices.size();
      mesh.m_vertices.swap(verts);
      mesh.m_indices.swap(idx);
      PLOGI << where << "re-cut with the outlines left out (chords within " << kDrawTolM * 1e3 << " mm inside the physics outline; the faces over them filled again by PolygonToTriangles): "
            << nTri << " triangles -> " << mesh.m_indices.size() / 3 << ", " << nV << " vertices -> " << mesh.m_vertices.size() << "; cuts " << names.str()
            << "; faces under the cuts left as they are: " << under << " (" << underUp << " wound to be seen from above)";
      // the walls and floors
      const std::filesystem::path file = g_app->m_fileLocator.GetAppPath(FileLocator::AppSubFolder::Assets, "PlayfieldCutWood.png");
      g_drawing.tex.reset(Texture::CreateFromFile(file));
      if (g_drawing.tex == nullptr || g_drawing.tex->m_width == 0 || g_drawing.tex->m_height == 0)
      {
         g_drawing.tex.reset();
         PLOGW << where << "no wood texture " << file.string() << ": the cuts' walls and floors are not drawn";
         return;
      }
      const double rep = drawn[0]->thick * (double)g_drawing.tex->m_width / (double)g_drawing.tex->m_height;   // the texture's height is the wall's depth; its width along the grain
      std::vector<Vertex3D_NoTex2> v;
      std::vector<unsigned int> wi;
      std::ostringstream o;
      for (const DrawnSlot *const d : drawn)
      {
         int nWall, nFloor, nSec, nAcross;
         Wood(*d, rep, v, wi, nWall, nFloor, nSec, nAcross);
         o << (d == drawn[0] ? "" : ", ") << d->name << " " << nWall << " wall + " << nFloor << " floor triangles (" << nSec << " sections x " << nAcross << " across per floor quarter)";
      }
      g_drawing.woodIndices = (unsigned int)wi.size();
      g_drawing.wood = std::make_shared<MeshBuffer>("PlayfieldCuts.Wood"s, std::make_shared<VertexBuffer>(renderer->m_renderDevice, (unsigned int)v.size(), (const float *)v.data()),
         std::make_shared<IndexBuffer>(renderer->m_renderDevice, wi), true);
      PLOGI << where << "walls and ramp floors: " << o.str() << "; " << v.size() << " vertices; walls " << drawn[0]->thick * 1e3 << " mm deep (the physics' through-hole); wood "
            << file.filename().string() << " " << g_drawing.tex->m_width << " x " << g_drawing.tex->m_height << " px: its height is the wall's depth, its width " << rep * 1e3
            << " mm along the grain";
   }

   void Draw(const Primitive *const pf, Renderer *const renderer, const Material *const mat)
   {
      if (g_drawing.prim != pf || g_drawing.wood == nullptr)
         return;
      RenderDevice *const rd = renderer->m_renderDevice;
      renderer->UpdateBasicShaderMatrix();   // (made in the world)
      rd->EnableAlphaBlend(false);
      rd->m_basicShader->SetBasic(mat, g_drawing.tex.get());
      rd->m_basicShader->SetTexture(ShaderUniform::tex_base_color, g_drawing.tex.get(), false, SamplerFilter::SF_ANISOTROPIC, SamplerAddressMode::SA_REPEAT, SamplerAddressMode::SA_REPEAT);
      rd->DrawMesh(rd->m_basicShader, false, pf->m_d.m_vPosition, pf->m_d.m_depthBias, g_drawing.wood, RenderDevice::TRIANGLELIST, 0, g_drawing.woodIndices);
   }

   void Release(Primitive *const pf)
   {
      if (g_drawing.prim != pf)
         return;
      pf->m_mesh.m_vertices.swap(g_drawing.vertices);
      pf->m_mesh.m_indices.swap(g_drawing.indices);
      g_drawing = Drawing();
   }
}
