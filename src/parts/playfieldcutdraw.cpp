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

   // cut d's sharp outline in the world (VPX units), counter-clockwise in x, y from the side t = -a at s = -hs: the end at +s (the round end of radius a to U1, the ramp u = L - tau^2,
   // w = c tau round its tip, the round end back to the side t = a), then the same turned half round (the end at -s); the points shared by two parts once
   std::vector<P2> Outline(const DrawnSlot &d)
   {
      const int nRound = std::max(1, (int)std::ceil(d.th1 / ChordAngle(d.a))), nRamp = RampChords(d);
      std::vector<P2> st;   // the end at +s in the slot frame (s, t), metres
      for (int k = 0; k <= nRound; k++)
         st.push_back({ d.hs + d.a * std::sin(d.th1 * k / nRound), -d.a * std::cos(d.th1 * k / nRound) });
      for (int k = 1; k < nRamp; k++)
      {
         const double tau = d.tauJ * (2.0 * k / nRamp - 1.0);
         st.push_back({ d.hs + d.L - tau * tau, d.c * tau });
      }
      for (int k = nRound; k >= 0; k--)
         st.push_back({ st[k].x, -st[k].y });
      std::vector<P2> out { World(d, -st.back().x, -st.back().y) };
      for (const double e : { 1.0, -1.0 })
         for (const P2 &q : st)
            if (const P2 w = World(d, e * q.x, e * q.y); std::hypot(w.x - out.back().x, w.y - out.back().y) > 1e-9)
               out.push_back(w);
      if (std::hypot(out[0].x - out.back().x, out[0].y - out.back().y) <= 1e-9)
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

   double Area(const std::vector<P2> &pt, const std::vector<int> &L)   // twice the signed area, counter-clockwise positive
   {
      double a = 0.0;
      for (size_t i = 0; i < L.size(); i++)
         a += Cross({ 0.0, 0.0 }, pt[L[i]], pt[L[(i + 1) % L.size()]]);
      return a;
   }

   // The playfield's faces sel (the corners wp, wz in the world) as one region less the outlines holes (counter-clockwise, convex, apart). The region's boundary is kept (its edges are the
   // neighbour faces' own: no cracks); its inside is filled again by VPX's ear clipping (PolygonToTriangles) of one simple polygon per boundary loop, each hole inside the loop joined to it
   // first by the shortest bridge that meets no edge. The other faces are copied to idx, the new triangles appended, wound as the faces are. A hole's point is a new vertex in verts with the
   // values (position, normal, texture coordinates) of the faces' plane: they must be one plane (an affine map of x, y) at the top z0 (within tol), as on the implicit playfield quad.
   // Returns why not, or nullptr. Every other failure (no face, a gap or a corner-only touch in them, a cut past their edge, no bridge) fails the one check at the end: n - 2 triangles per
   // polygon, none turned over, covering the region less the holes.
   const char *Fill(const Mesh &mesh, const std::vector<P2> &wp, const std::vector<double> &wz, const double z0, const double tol, const std::vector<char> &sel,
                    const std::vector<std::vector<P2>> &holes, std::vector<Vertex3D_NoTex2> &verts, std::vector<unsigned int> &idx)
   {
      static_assert(sizeof(Vertex3D_NoTex2) == 8 * sizeof(float));   // x, y, z, nx, ny, nz, tu, tv: read as 8 floats from &x
      const char *const noFill = "the faces over the cuts cannot be filled again round the cuts (none, a gap or a corner-only touch in them, a cut past their edge, or the ear clipping fails)";
      std::map<std::pair<double, double>, int> placeOf;   // the region's corners by place (a mesh may have several vertices at one place)
      std::vector<P2> pt;                                 // the corners' places, then the holes' points
      std::vector<unsigned int> vtx;                      // their vertex in verts
      std::set<std::pair<int, int>> edges;                // the faces' edges, counter-clockwise in x, y
      int big = -1;                                       // the largest face: the plane's values and the winding from its corners
      double bigArea = 0.0;
      for (size_t t = 0; t < mesh.m_indices.size() / 3; t++)
      {
         const unsigned int *I = &mesh.m_indices[3 * t];
         if (!sel[t])
         {
            idx.insert(idx.end(), I, I + 3);
            continue;
         }
         const double ar = Cross(wp[I[0]], wp[I[1]], wp[I[2]]);
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
         for (int k = 0; k < 3; k++)
            edges.insert(ar > 0.0 ? std::make_pair(p[k], p[(k + 1) % 3]) : std::make_pair(p[(k + 1) % 3], p[k]));
      }
      if (big < 0)
         return noFill;
      // the plane's values at a place: the largest face's corners' by their weights; every corner of the region must be on it
      const unsigned int *IB = &mesh.m_indices[3 * big];
      const P2 K[3] = { wp[IB[0]], wp[IB[1]], wp[IB[2]] };
      const double A = Cross(K[0], K[1], K[2]);
      const auto plane = [&](const P2 &q) {
         const double w0 = Cross(q, K[1], K[2]) / A, w1 = Cross(K[0], q, K[2]) / A, w[3] = { w0, w1, 1.0 - w0 - w1 };
         Vertex3D_NoTex2 nv;
         for (int f = 0; f < 8; f++)
            (&nv.x)[f] = (float)(w[0] * (&mesh.m_vertices[IB[0]].x)[f] + w[1] * (&mesh.m_vertices[IB[1]].x)[f] + w[2] * (&mesh.m_vertices[IB[2]].x)[f]);
         return nv;
      };
      for (size_t i = 0; i < pt.size(); i++)
      {
         const Vertex3D_NoTex2 a = verts[vtx[i]], b = plane(pt[i]);
         for (int f = 2; f < 8; f++)   // the height at the top within tol, the normal within 1e-3, the texture coordinates within 1e-4
            if (f == 2 ? std::abs(wz[vtx[i]] - z0) > tol : std::abs((&a.x)[f] - (&b.x)[f]) > (f < 6 ? 1e-3f : 1e-4f))
               return "the faces over the cuts are not one flat plane at the playfield's top (raised or curved there, design section 6; or a seam in normals or texture coordinates)";
      }
      // the region's boundary: the edges whose reverse is no edge, in loops, each counter-clockwise (an outer edge; a clockwise one is a gap in the faces)
      std::map<int, int> next;
      for (const auto &[a, b] : edges)
         if (edges.count({ b, a }) == 0)
            next.emplace(a, b);   // (a corner with two: one is lost, and a loop does not close)
      std::vector<std::vector<int>> loops;
      std::set<int> seen;
      for (const auto &[first, unused] : next)
      {
         if (seen.count(first))
            continue;
         std::vector<int> L;
         for (int s = first; L.empty() || s != first; s = next[s])
         {
            if (!seen.insert(s).second || next.count(s) == 0)
               return noFill;
            L.push_back(s);
         }
         if (Area(pt, L) <= 0.0)
            return noFill;
         loops.push_back(std::move(L));
      }
      // the holes' points, clockwise, each a new vertex on the plane
      std::vector<std::vector<int>> H;
      double holeArea = 0.0;
      for (const std::vector<P2> &O : holes)
      {
         H.emplace_back();
         for (size_t k = O.size(); k-- > 0;)
         {
            H.back().push_back((int)pt.size());
            pt.push_back(O[k]);
            vtx.push_back((unsigned int)verts.size());
            verts.push_back(plane(O[k]));
         }
         holeArea -= Area(pt, H.back());
      }
      // VPX's ear clipping works in floats: the points relative to the first hole's first point (0.1 um at a playfield's size)
      std::vector<RenderVertex> rgv(pt.size());
      for (size_t i = 0; i < pt.size(); i++)
      {
         rgv[i].x = (float)(pt[i].x - holes[0][0].x);
         rgv[i].y = (float)(pt[i].y - holes[0][0].y);
      }
      const auto inside = [&](const std::vector<int> &L, const P2 &q) {   // (crossing number)
         bool in = false;
         for (size_t k = 0, l = L.size() - 1; k < L.size(); l = k++)
         {
            const P2 &a = pt[L[k]], &b = pt[L[l]];
            if ((a.y > q.y) != (b.y > q.y) && q.x < a.x + (q.y - a.y) * (b.x - a.x) / (b.y - a.y))
               in = !in;
         }
         return in;
      };
      const auto clear = [&](const int a, const int b, const std::vector<int> &E) {   // the segment a b meets (crosses or touches: VPX's test) no edge of E but those ending at a or b
         for (size_t k = 0; k < E.size(); k++)
         {
            const int q1 = E[k], q2 = E[(k + 1) % E.size()];
            if (q1 != a && q1 != b && q2 != a && q2 != b && FLinesIntersect(&rgv[a], &rgv[b], &rgv[q1], &rgv[q2]))
               return false;
         }
         return true;
      };
      std::vector<char> twice(pt.size(), 0);   // a bridge's ends (met twice in the polygon): not a bridge's end again
      double loopArea = 0.0, triArea = 0.0;
      for (std::vector<int> &P : loops)
      {
         loopArea += Area(pt, P);
         std::vector<std::vector<int>> in;   // the holes with all their points inside the loop
         for (const std::vector<int> &h : H)
            if (std::all_of(h.begin(), h.end(), [&](const int q) { return inside(P, pt[q]); }))
               in.push_back(h);
         while (!in.empty())
         {
            int bh = -1, bi = -1, bj = -1;
            double best = 1e300;
            for (size_t h = 0; h < in.size(); h++)
               for (size_t i = 0; i < in[h].size(); i++)
                  for (size_t j = 0; j < P.size(); j++)
                  {
                     const int a = in[h][i], b = P[j];
                     const double dd = (pt[a].x - pt[b].x) * (pt[a].x - pt[b].x) + (pt[a].y - pt[b].y) * (pt[a].y - pt[b].y);
                     if (dd >= best || twice[b] || !clear(a, b, P) || !std::all_of(in.begin(), in.end(), [&](const std::vector<int> &E) { return clear(a, b, E); }))
                        continue;
                     best = dd;
                     bh = (int)h;
                     bi = (int)i;
                     bj = (int)j;
                  }
            if (bh < 0)
               return noFill;
            const std::vector<int> &Hb = in[bh];
            twice[P[bj]] = twice[Hb[bi]] = 1;
            std::vector<int> merged(P.begin(), P.begin() + bj + 1);   // ..., the polygon's point, round the hole from its point back to it, the polygon's point again, ...
            for (size_t k = 0; k <= Hb.size(); k++)
               merged.push_back(Hb[(bi + k) % Hb.size()]);
            merged.insert(merged.end(), P.begin() + bj, P.end());
            P.swap(merged);
            in.erase(in.begin() + bh);
         }
         // VPX's ear clipping takes the polygon clockwise in x, y (as VPX's counter-clockwise outlines with y down) and gives each triangle (a, c, b): counter-clockwise
         std::vector<unsigned int> poly(P.rbegin(), P.rend()), tri;
         PolygonToTriangles(rgv, poly, tri, false);
         if (tri.size() != 3 * (P.size() - 2))
            return noFill;
         for (size_t k = 0; k < tri.size(); k += 3)
         {
            const double a = Cross(pt[tri[k]], pt[tri[k + 1]], pt[tri[k + 2]]);
            triArea += a < 0.0 ? -1e300 : a;   // (turned over: fails the check)
            idx.insert(idx.end(), { vtx[tri[k]], vtx[tri[k + (A > 0.0 ? 1 : 2)]], vtx[tri[k + (A > 0.0 ? 2 : 1)]] });
         }
      }
      return std::abs(triArea - (loopArea - holeArea)) > 1e-6 * loopArea ? noFill : nullptr;
   }

   // Cut d's walls and ramp floors in the wood (world, VPX units; normals into the hole), appended to v and idx. The slot's frame (s, t, z up) turns as the world's (p is l turned by +90 deg,
   // PlayfieldCuts::Setup), so each quad is wound from its sides' directions to be seen from its normals' side, as the playfield's top (counter-clockwise seen from above). Texture: the
   // grain along the wall at the texture's own scale, its height the wall's depth: u = length along the wall (on the floor: along the slot) / rep, v = depth / thick (the floor shows the
   // layer of wood at its depth, as the wall beside it).
   void Wood(const DrawnSlot &d, const double rep, std::vector<Vertex3D_NoTex2> &v, std::vector<unsigned int> &idx)
   {
      const auto vert = [&](const double s, const double t, const double z, const double ns, const double nt, const double nz, const double tu) {
         const P2 w = World(d, s, t);
         const double nx = d.lx * ns + d.px * nt, ny = d.ly * ns + d.py * nt, l = std::sqrt(nx * nx + ny * ny + nz * nz);
         v.push_back({ (float)w.x, (float)w.y, (float)((d.z0 + z) / d.unit), (float)(nx / l), (float)(ny / l), (float)(nz / l), (float)tu, (float)(-z / d.thick) });
      };
      // the quads of rows x cols vertices from f0 (row i, column j at f0 + i cols + j): corners (i, j) (i + 1, j) (i + 1, j + 1) (i, j + 1) in turn, or the other way (flip)
      const auto grid = [&](const unsigned int f0, const unsigned int rows, const unsigned int cols, const bool flip) {
         for (unsigned int p = f0; p + cols + 1 < f0 + rows * cols; p++)
            if ((p - f0) % cols + 1 < cols)
               idx.insert(idx.end(), { p, flip ? p + 1 : p + cols, p + cols + 1, p, p + cols + 1, flip ? p + cols : p + 1 });
      };
      // the ramp floor's cross-sections, u past the end of the straight part: its outline's chord ends (the hole's points), the through-hole's round end's chords where it cuts the floor
      // (U1 < u < a), and the tip (an arc of no width; its radius stays finite)
      const int nRamp = RampChords(d);
      const double dPh = ChordAngle(d.a);
      const int nRound = std::max(1, (int)std::ceil(d.th1 / dPh)), nCut = std::max(1, (int)std::ceil((0.5 * M_PI - d.th1) / dPh));
      std::vector<double> us { d.U1, d.a, d.L * (1.0 - 1e-7) };
      for (int k = 1; k < nRamp; k++)
         if (const double tau = d.tauJ * (2.0 * k / nRamp - 1.0); tau > 0.0)
            us.push_back(d.L - tau * tau);
      for (int k = 1; k < nCut; k++)
         us.push_back(d.a * std::sin(d.th1 + (0.5 * M_PI - d.th1) * k / nCut));
      std::sort(us.begin(), us.end());
      us.erase(std::unique(us.begin(), us.end(), [](const double x, const double y) { return std::abs(x - y) < 1e-12; }), us.end());
      struct Sec
      {
         double u, rho, zc, drho, dzc, p0, p1;   // the floor's arc t^2 + (z - zc)^2 = rho^2 and its change along u; from the angle p0 (on the through-hole's wall, or the axis) to p1 (the outline)
      };
      std::vector<Sec> sec;
      int nA = 1;   // chords across each section (the most any needs)
      for (const double u : us)
      {
         Sec q { u };
         d.FloorArc(u, q.rho, q.zc);
         const double h = std::min(1e-7, 0.5 * (d.L - u));
         double r1, z1;
         d.FloorArc(u + h, r1, z1);
         q.drho = (r1 - q.rho) / h;
         q.dzc = (z1 - q.zc) / h;
         q.p0 = -std::acos(std::min(1.0, d.C(u) / q.rho));
         q.p1 = std::atan2(-q.zc, d.W(u));
         nA = std::max(nA, (int)std::ceil((q.p1 - q.p0) / ChordAngle(q.rho)));
         sec.push_back(q);
      }
      // the walls: the through-hole's outline (the straight sides, the round ends of radius a) from the side t = -a at s = -hs, counter-clockwise: the end at +s (E) from that side round to
      // the other, its top edge on the top up to U1, then on the floor (the sections U1 < u <= a, at the floor's angle p0); the end at -s, E turned half round; back to the start (the seam)
      struct WP
      {
         double s, t, z, ns, nt;
      };
      std::vector<WP> E;
      for (int k = 0; k <= nRound; k++)
      {
         const double ph = d.th1 * k / nRound;
         E.push_back({ d.hs + d.a * std::sin(ph), -d.a * std::cos(ph), 0.0, -std::sin(ph), std::cos(ph) });
      }
      for (const Sec &q : sec)
         if (q.u > d.U1 && q.u <= d.a)
            E.push_back({ d.hs + q.u, -q.rho * std::cos(q.p0), q.zc + q.rho * std::sin(q.p0), -q.u / d.a, q.rho * std::cos(q.p0) / d.a });
      for (size_t k = E.size() - 1; k-- > 0;)   // round the other side (the axis once)
         E.push_back({ E[k].s, -E[k].t, E[k].z, E[k].ns, -E[k].nt });
      const auto turn = [](const WP &p) { return WP { -p.s, -p.t, p.z, -p.ns, -p.nt }; };
      std::vector<WP> loop { turn(E.back()) };
      loop.insert(loop.end(), E.begin(), E.end());
      for (const WP &p : E)
         loop.push_back(turn(p));
      const unsigned int w0 = (unsigned int)v.size();
      double len = 0.0;
      for (size_t k = 0; k < loop.size(); k++)
      {
         len += k > 0 ? std::hypot(loop[k].s - loop[k - 1].s, loop[k].t - loop[k - 1].t) : 0.0;
         vert(loop[k].s, loop[k].t, loop[k].z, loop[k].ns, loop[k].nt, 0.0, len / rep);
         vert(loop[k].s, loop[k].t, -d.thick, loop[k].ns, loop[k].nt, 0.0, len / rep);
      }
      grid(w0, (unsigned int)loop.size(), 2, false);   // (top edge forward, then down: the hole on the left, seen from above)
      // the floors: at each end e, each side of the axis g, the sections by nA chords across; normal: the arc surface's gradient turned inward. Mirrored once (e g < 0): wound the other way
      for (const double e : { 1.0, -1.0 })
         for (const double g : { 1.0, -1.0 })
         {
            const unsigned int f0 = (unsigned int)v.size();
            for (const Sec &q : sec)
               for (int j = 0; j <= nA; j++)
               {
                  const double ps = q.p0 + (q.p1 - q.p0) * j / nA, s = e * (d.hs + q.u);
                  vert(s, g * q.rho * std::cos(ps), q.zc + q.rho * std::sin(ps), e * (std::sin(ps) * q.dzc + q.drho), -g * std::cos(ps), -std::sin(ps), s / rep);
               }
            grid(f0, (unsigned int)sec.size(), nA + 1, e * g < 0.0);
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
         const Vertex3Ds w = M.MultiplyVectorNoPerspective(Vertex3Ds(mesh.m_vertices[i].x, mesh.m_vertices[i].y, mesh.m_vertices[i].z));
         wp[i] = { w.x, w.y };
         wz[i] = w.z;
      }
      // each cut's outline (one overlapping an earlier cut is not drawn) and the faces it overlaps seen from above, but those wholly under the top (the playfield's underside, left as they
      // are): re-cut (Fill: they must be flat at the top)
      std::vector<std::vector<P2>> outl;
      std::vector<const DrawnSlot *> drawn;
      std::vector<char> sel(nTri, 0);
      std::ostringstream names;
      const double z0 = cuts[0].z0 / cuts[0].unit, tol = kTopTolM / cuts[0].unit;
      for (const DrawnSlot &d : cuts)
      {
         std::vector<P2> o = Outline(d);
         if (std::any_of(outl.begin(), outl.end(), [&](const std::vector<P2> &q) { return Overlap(o.data(), o.size(), q.data(), q.size()); }))
         {
            PLOGW << where << "cut " << d.name << " overlaps another cut: not drawn";
            continue;
         }
         for (size_t t = 0; t < nTri; t++)
         {
            const unsigned int *I = &mesh.m_indices[3 * t];
            const double ar = Cross(wp[I[0]], wp[I[1]], wp[I[2]]);   // (0: seen edge-on from above, covers nothing)
            const P2 K[3] = { wp[I[0]], wp[I[ar > 0.0 ? 1 : 2]], wp[I[ar > 0.0 ? 2 : 1]] };
            if (ar != 0.0 && std::max({ wz[I[0]], wz[I[1]], wz[I[2]] }) >= z0 - tol && Overlap(K, 3, o.data(), o.size()))
               sel[t] = 1;
         }
         names << (outl.empty() ? "" : ", ") << d.name << " (" << o.size() << " outline points)";
         outl.push_back(std::move(o));
         drawn.push_back(&d);
      }
      std::vector<Vertex3D_NoTex2> verts = mesh.m_vertices;
      std::vector<unsigned int> idx;
      if (const char *const why = Fill(mesh, wp, wz, z0, tol, sel, outl, verts, idx))
      {
         PLOGW << where << why << ": the playfield is drawn uncut";
         return;
      }
      g_drawing.prim = pf;
      g_drawing.vertices = std::exchange(mesh.m_vertices, std::move(verts));
      g_drawing.indices = std::exchange(mesh.m_indices, std::move(idx));
      PLOGI << where << "re-cut with the outlines left out (chords within " << kDrawTolM * 1e3 << " mm inside the physics outline; the faces over them filled again by PolygonToTriangles): "
            << nTri << " triangles -> " << mesh.m_indices.size() / 3 << ", " << g_drawing.vertices.size() << " vertices -> " << mesh.m_vertices.size() << "; cuts " << names.str();
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
      for (const DrawnSlot *const d : drawn)
         Wood(*d, rep, v, wi);
      g_drawing.woodIndices = (unsigned int)wi.size();
      g_drawing.wood = std::make_shared<MeshBuffer>("PlayfieldCuts.Wood"s, std::make_shared<VertexBuffer>(renderer->m_renderDevice, (unsigned int)v.size(), (const float *)v.data()),
         std::make_shared<IndexBuffer>(renderer->m_renderDevice, wi), true);
      PLOGI << where << "walls and ramp floors: " << wi.size() / 3 << " triangles, " << v.size() << " vertices; walls " << drawn[0]->thick * 1e3 << " mm deep (the physics' through-hole); wood "
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
