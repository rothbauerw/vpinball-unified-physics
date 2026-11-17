// license:GPLv3+

// The mesh holes as drawn: each hole exactly as the physics cut it (physics/meshholes.h MeshHoles::Drawn, recorded by MeshHoles::Build): the joined rim loops (the same points), and the
// physics' own wall, floor and side triangles. Made once with the playfield's render buffers (Primitive::RenderSetup): the rims are cut in the playfield's own mesh, so every pass that draws
// the playfield (static prepass, dynamic pass, VR) draws them; the walls and floors are one mesh drawn right after the playfield (Primitive::Render) with its material and the wood texture
// assets/PlayfieldCutWood.png. Per frame: that one draw call.

#include "core/stdafx.h"
#include "parts/playfieldcutdraw.h"
#include "parts/primitive.h"
#include "physics/meshholes.h"
#include "physics/physicsini.h"
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
   using MeshHoles::DrawnHole;

   constexpr double kTopTolM = 1e-6;   // a face is at the playfield's top over a cut when its corners are within this of it (the drawing's own flatness limit)

   struct P2
   {
      double x, y;
   };
   double Cross(const P2 &o, const P2 &a, const P2 &b) { return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x); }

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

   // The playfield's faces sel (the corners wp, wz in the world) as one region less the outlines holes (counter-clockwise, simple, apart). The region's boundary is kept (its edges are the
   // neighbour faces' own: no cracks); its inside is filled again by VPX's ear clipping (PolygonToTriangles) of one simple polygon per boundary loop, each hole inside the loop joined to it
   // first by the shortest bridge that meets no edge. The other faces are copied to idx, the new triangles appended, wound as the faces are. A hole's point is a new vertex in verts with the
   // values (position, normal, texture coordinates) of the faces' plane: they must be one plane (an affine map of x, y) at the top z0 (within tol), as on the implicit playfield quad.
   // Returns why not, or nullptr. Every other failure (no face, a gap in them, a cut past their edge, no bridge) fails the one check at the end: n - 2 triangles per
   // polygon, none turned over, covering the region less the holes.
   const char *Fill(const Mesh &mesh, const std::vector<P2> &wp, const std::vector<double> &wz, const double z0, const double tol, const std::vector<char> &sel,
                    const std::vector<std::vector<P2>> &holes, std::vector<Vertex3D_NoTex2> &verts, std::vector<unsigned int> &idx)
   {
      static_assert(sizeof(Vertex3D_NoTex2) == 8 * sizeof(float));   // x, y, z, nx, ny, nz, tu, tv: read as 8 floats from &x
      // which check refused (brief 4: the chamfered hole's 'drawn uncut' in run a2 named none of them)
      const auto noFill = [](const int k) -> const char * {
         static const char *const why[5] = { "the faces over the cuts cannot be filled again round the cuts: no face of the drawn top under them",
                                             "the faces over the cuts cannot be filled again round the cuts: their boundary does not close in loops (a cut past their edge)",
                                             "the faces over the cuts cannot be filled again round the cuts: a boundary loop runs clockwise (a gap in the faces)",
                                             "the faces over the cuts cannot be filled again round the cuts: no bridge from a hole to its loop clear of every edge",
                                             "the faces over the cuts cannot be filled again round the cuts: the ear clipping gave too few triangles" };
         return why[k];
      };
      static std::string areaWhy;   // the last check's numbers
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
         return noFill(0);
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
      // A corner where faces touch only at that corner (brief 5: a fan's centre, two of its faces under cuts and the two between them not) has two boundary edges out: the loop
      // through it goes on along the one turning most to the left (the same face's side: its inside is on the left), so each loop keeps to one side and a loop closes
      std::multimap<int, int> next;
      for (const auto &[a, b] : edges)
         if (edges.count({ b, a }) == 0)
            next.emplace(a, b);
      const auto after = [&](const int p, const int s) {   // the boundary edge out of s after the edge p s, or -1
         int best = -1;
         double bestTurn = -10.0;
         const auto [lo, hi] = next.equal_range(s);
         for (auto it = lo; it != hi; ++it)
         {
            const int w = it->second;
            const double dx = pt[s].x - pt[p].x, dy = pt[s].y - pt[p].y, ex = pt[w].x - pt[s].x, ey = pt[w].y - pt[s].y;
            const double turn = std::atan2(dx * ey - dy * ex, dx * ex + dy * ey);   // + left
            if (turn > bestTurn)
            {
               bestTurn = turn;
               best = w;
            }
         }
         return best;
      };
      std::vector<std::vector<int>> loops;
      std::set<std::pair<int, int>> seen;
      for (const auto &[a0, b0] : next)
      {
         if (seen.count({ a0, b0 }))
            continue;
         std::vector<int> L;
         for (std::pair<int, int> e { a0, b0 }; L.empty() || e != std::make_pair(a0, b0);)
         {
            if (!seen.insert(e).second)
               return noFill(1);
            L.push_back(e.first);
            const int w = after(e.first, e.second);
            if (w < 0)
               return noFill(1);
            e = { e.second, w };
         }
         if (Area(pt, L) <= 0.0)
            return noFill(2);
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
               return noFill(3);
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
            return noFill(4);
         for (size_t k = 0; k < tri.size(); k += 3)
         {
            const double a = Cross(pt[tri[k]], pt[tri[k + 1]], pt[tri[k + 2]]);
            triArea += a < 0.0 ? -1e300 : a;   // (turned over: fails the check)
            idx.insert(idx.end(), { vtx[tri[k]], vtx[tri[k + (A > 0.0 ? 1 : 2)]], vtx[tri[k + (A > 0.0 ? 2 : 1)]] });
         }
      }
      if (std::abs(triArea - (loopArea - holeArea)) <= 1e-6 * loopArea)
         return nullptr;
      std::ostringstream w;   // (a hole not wholly inside one loop is never bridged: its area is missing from the triangles; a turned-over triangle: -1e300)
      w << "the faces over the cuts cannot be filled again round the cuts: the triangles' area " << triArea << " is not the loops' " << loopArea << " less the holes' " << holeArea << " ("
        << loops.size() << " loops, " << H.size() << " holes)";
      areaWhy = w.str();
      return areaWhy.c_str();
   }

   // the convex hull of q, counter-clockwise (Andrew's monotone chain): the faces a rim loop overlaps are found by it (Overlap takes convex polygons; a face over the hull but not the
   // loop is filled again as it was)
   std::vector<P2> Hull(std::vector<P2> q)
   {
      std::sort(q.begin(), q.end(), [](const P2 &a, const P2 &b) { return a.x < b.x || (a.x == b.x && a.y < b.y); });
      std::vector<P2> h(2 * q.size());
      size_t k = 0;
      for (size_t i = 0; i < q.size(); i++)
      {
         while (k >= 2 && Cross(h[k - 2], h[k - 1], q[i]) <= 0.0)
            k--;
         h[k++] = q[i];
      }
      for (size_t i = q.size() - 1, t = k + 1; i-- > 0;)
      {
         while (k >= t && Cross(h[k - 2], h[k - 1], q[i]) <= 0.0)
            k--;
         h[k++] = q[i];
      }
      h.resize(k - 1);
      return h;
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
      const std::vector<DrawnHole> &holes = MeshHoles::Drawn();   // every joined hole as the physics cut it (MeshHoles::Build: in the physics engine's setup, before the render setup, every load)
      if (holes.empty() || g_drawing.prim != nullptr)
         return;
      Mesh &mesh = pf->m_mesh;
      const std::string where = "Mesh holes drawn, playfield " + pf->GetName() + ": ";
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
      // each hole's rim loops as the physics cut them, and the faces of the drawn top they overlap seen from above (their convex hull's), but those wholly under the top (the playfield's
      // underside, left as they are): re-cut (Fill: they must be flat at the top). The drawn top under a hole: the highest of the drawn faces' planes at its first rim point
      const double tol = kTopTolM / holes[0].unit;
      std::vector<std::vector<P2>> outl;
      std::vector<const DrawnHole *> owner;
      std::vector<char> sel(nTri, 0);
      double z0 = 0.0;
      for (const DrawnHole &h : holes)
      {
         std::ostringstream line;
         const bool same = PhysicsIni::Lower(h.hostName) == PhysicsIni::Lower(pf->GetName());
         line << where << h.name << ": cut by the physics in the " << h.hostKind << " host " << h.hostName << (same ? " (the drawn primitive)" : " (not the drawn primitive: drawn on it)")
              << ", walls " << (h.walls ? "yes, " : "no, ") << h.thick << " VPU deep, physics top z " << h.topZ << " VPU";
         if (!h.draw)
         {
            PLOGI << line.str() << "; draw = 0: cut in the physics, not drawn (the drawn top stays whole there)";
            continue;
         }
         if (!h.onTop || h.rim.empty() || h.rim[0].size() < 2)
         {
            PLOGI << line.str() << (h.onTop ? "; no rim recorded: not cut in the drawn top" : "; cut into another hole's floor: its rim is not cut in the drawn top, and that floor is drawn whole over it");
            continue;
         }
         const P2 q0 { h.rim[0][0], h.rim[0][1] };
         double zt = 0.0;
         bool found = false;
         for (size_t t = 0; t < nTri; t++)
         {
            const unsigned int *I = &mesh.m_indices[3 * t];
            const double A = Cross(wp[I[0]], wp[I[1]], wp[I[2]]);
            if (A == 0.0)
               continue;
            const double w0 = Cross(q0, wp[I[1]], wp[I[2]]) / A, w1 = Cross(wp[I[0]], q0, wp[I[2]]) / A, w2 = 1.0 - w0 - w1;
            if (w0 < 0.0 || w1 < 0.0 || w2 < 0.0)
               continue;
            const double z = w0 * wz[I[0]] + w1 * wz[I[1]] + w2 * wz[I[2]];
            if (!found || z > zt)
               zt = z;
            found = true;
         }
         if (!found)
         {
            PLOGI << line.str() << "; no face of the drawn primitive under its rim: not cut in the drawn top";
            continue;
         }
         if (outl.empty())
            z0 = zt;
         line << ", drawn top z " << zt << " VPU (drawn - physics " << zt - h.topZ << " VPU)";
         for (const std::vector<double> &L : h.rim)
         {
            std::vector<P2> o;
            for (size_t i = 0; i + 1 < L.size(); i += 2)
               o.push_back({ L[i], L[i + 1] });
            double ar = 0.0;
            for (size_t i = 0; i < o.size(); i++)
               ar += Cross({ 0.0, 0.0 }, o[i], o[(i + 1) % o.size()]);
            if (o.size() < 3 || ar <= 0.0)
            {
               line << "; a rim loop of " << o.size() << " points not counter-clockwise: not cut";
               continue;
            }
            const std::vector<P2> hull = Hull(o);
            for (size_t t = 0; t < nTri; t++)
            {
               const unsigned int *I = &mesh.m_indices[3 * t];
               const double a = Cross(wp[I[0]], wp[I[1]], wp[I[2]]);   // (0: seen edge-on from above, covers nothing)
               const P2 K[3] = { wp[I[0]], wp[I[a > 0.0 ? 1 : 2]], wp[I[a > 0.0 ? 2 : 1]] };
               if (a != 0.0 && std::max({ wz[I[0]], wz[I[1]], wz[I[2]] }) >= zt - tol && Overlap(K, 3, hull.data(), hull.size()))
                  sel[t] = 1;
            }
            outl.push_back(std::move(o));
            owner.push_back(&h);
         }
         PLOGI << line.str();
      }
      if (!outl.empty())
      {
         std::vector<Vertex3D_NoTex2> verts = mesh.m_vertices;
         std::vector<unsigned int> idx;
         const size_t base = verts.size();
         if (const char *const why = Fill(mesh, wp, wz, z0, tol, sel, outl, verts, idx))
         {
            PLOGW << where << why << ": the playfield is drawn uncut";
            return;
         }
         // each rim loop as drawn (Fill's new vertices, float, in the order Fill added them: each loop's points reversed) against the physics' rim points
         size_t at = base;
         for (size_t i = 0; i < outl.size(); i++)
         {
            double dmax = 0.0;
            for (size_t k = 0; k < outl[i].size(); k++)
            {
               const Vertex3D_NoTex2 &v = verts[at + outl[i].size() - 1 - k];   // (the primitive's own frame: to the world as its corners)
               const Vertex3Ds w = M.MultiplyVectorNoPerspective(Vertex3Ds(v.x, v.y, v.z));
               dmax = std::max(dmax, std::hypot(w.x - outl[i][k].x, w.y - outl[i][k].y));
            }
            at += outl[i].size();
            PLOGI << where << owner[i]->name << ": rim points drawn " << outl[i].size() << " vs physics " << outl[i].size()
                  << ", largest plan distance between matching vertices " << dmax << " VPU (the drawn vertex's float rounding)";
         }
         g_drawing.vertices = std::exchange(mesh.m_vertices, std::move(verts));
         g_drawing.indices = std::exchange(mesh.m_indices, std::move(idx));
         PLOGI << where << "re-cut with the physics' rim loops left out (the faces over them filled again by PolygonToTriangles): " << nTri << " triangles -> " << mesh.m_indices.size() / 3
               << ", " << g_drawing.vertices.size() << " vertices -> " << mesh.m_vertices.size() << ", " << outl.size() << " rim loops";
      }
      g_drawing.prim = pf;
      // the walls, floors and side faces: the physics' own triangles, each with its flat normal (a rim wall faces into the hole, a floor up, as the physics wound them; a wall under a
      // floor's open edge and a side face turned to face the hole's first rim loop's centre). Texture: the grain along the rim on walls, its height the wall's depth
      const std::filesystem::path file = g_app->m_fileLocator.GetAppPath(FileLocator::AppSubFolder::Assets, "PlayfieldCutWood.png");
      g_drawing.tex.reset(Texture::CreateFromFile(file));
      if (g_drawing.tex == nullptr || g_drawing.tex->m_width == 0 || g_drawing.tex->m_height == 0)
      {
         g_drawing.tex.reset();
         PLOGW << where << "no wood texture " << file.string() << ": the holes' walls and floors are not drawn";
         return;
      }
      const double rep = holes[0].thick * (double)g_drawing.tex->m_width / (double)g_drawing.tex->m_height;   // the texture's height is the wall's depth; its width along the grain
      std::vector<Vertex3D_NoTex2> v;
      std::vector<unsigned int> wi;
      for (const DrawnHole &h : holes)
      {
         if (!h.draw)   // draw = 0: no walls or floors drawn either
            continue;
         P2 ctr { 0.0, 0.0 };
         const size_t nr = h.rim.empty() ? 0 : h.rim[0].size() / 2;
         for (size_t i = 0; i < nr; i++)
            ctr = { ctr.x + h.rim[0][2 * i] / (double)nr, ctr.y + h.rim[0][2 * i + 1] / (double)nr };
         int walls = 0, others = 0;
         for (size_t i = 0; i + 18 <= h.wood.size(); i += 18)
         {
            const float *f = &h.wood[i];
            const int kind = (int)f[5];
            const double ux = f[6] - f[0], uy = f[7] - f[1], uz = f[8] - f[2], vx = f[12] - f[0], vy = f[13] - f[1], vz = f[14] - f[2];
            double nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
            const double l = std::sqrt(nx * nx + ny * ny + nz * nz);
            if (l <= 0.0)
               continue;
            nx /= l;
            ny /= l;
            nz /= l;
            const bool flip = kind >= 2 && nr > 0 && nx * (ctr.x - (f[0] + f[6] + f[12]) / 3.0) + ny * (ctr.y - (f[1] + f[7] + f[13]) / 3.0) < 0.0;
            for (const int k : { 0, flip ? 2 : 1, flip ? 1 : 2 })
            {
               const float *q = f + 6 * k;
               wi.push_back((unsigned int)v.size());
               v.push_back({ q[0], q[1], q[2], (float)(flip ? -nx : nx), (float)(flip ? -ny : ny), (float)(flip ? -nz : nz), (float)(q[3] / rep), (float)(q[4] / h.thick) });
            }
            (kind == 1 || kind == 2 ? walls : others)++;
         }
         // a hole the physics cut without walls (narrower than the ball: Button, Star, slot, a small cut): drawn walls all the same (user 2026-11-17), one quad per rim edge from the
         // top down to the wall depth (h.thick) at the physics top's height (h.topZ), facing into the hole (the hole lies left of every rim piece). Under a floor (the slot's ramps,
         // the groove's strips) the floor hides the part below it. Looks only: no physics face
         int drawnOnly = 0;
         if (!h.walls)
            for (const std::vector<double> &loop : h.rim)
            {
               const size_t n = loop.size() / 2;
               double along = 0.0;
               for (size_t i = 0; i < n; i++)
               {
                  const size_t j = (i + 1) % n;
                  const double ax = loop[2 * i], ay = loop[2 * i + 1], bx = loop[2 * j], by = loop[2 * j + 1], len = std::hypot(bx - ax, by - ay);
                  if (len <= 0.0)
                     continue;
                  const float nx = (float)(-(by - ay) / len), ny = (float)((bx - ax) / len);   // the left normal: into the hole
                  const float zt = (float)h.topZ, zb = (float)(h.topZ - h.thick);
                  const Vertex3D_NoTex2 At { (float)ax, (float)ay, zt, nx, ny, 0.f, (float)(along / rep), 0.f }, Bt { (float)bx, (float)by, zt, nx, ny, 0.f, (float)((along + len) / rep), 0.f };
                  const Vertex3D_NoTex2 Ab { (float)ax, (float)ay, zb, nx, ny, 0.f, (float)(along / rep), 1.f }, Bb { (float)bx, (float)by, zb, nx, ny, 0.f, (float)((along + len) / rep), 1.f };
                  along += len;
                  for (const Vertex3D_NoTex2 &q : { At, Bt, Bb, At, Bb, Ab })
                  {
                     wi.push_back((unsigned int)v.size());
                     v.push_back(q);
                  }
                  drawnOnly += 2;
               }
            }
         PLOGI << where << h.name << ": walls drawn " << walls << " triangles vs physics " << h.physWallTris << ", floors and sides drawn " << others << " triangles"
               << (h.walls ? std::string() : ", walls drawn without physics walls " + std::to_string(drawnOnly) + " triangles (" + std::to_string(h.thick) + " VPU deep)");
      }
      if (wi.empty())
         return;
      g_drawing.woodIndices = (unsigned int)wi.size();
      g_drawing.wood = std::make_shared<MeshBuffer>("PlayfieldCuts.Wood"s, std::make_shared<VertexBuffer>(renderer->m_renderDevice, (unsigned int)v.size(), (const float *)v.data()),
         std::make_shared<IndexBuffer>(renderer->m_renderDevice, wi), true);
      PLOGI << where << "walls, floors and sides: " << wi.size() / 3 << " triangles, " << v.size() << " vertices; wood " << file.filename().string() << " " << g_drawing.tex->m_width << " x "
            << g_drawing.tex->m_height << " px: its height is the wall's depth, its width " << rep << " VPU along the grain";
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
