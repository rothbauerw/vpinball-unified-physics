// license:GPLv3+

#pragma once

#include <filesystem>
#include <string>
#include <vector>

// Holes cut in a static mesh host's upward faces at load, in place, before the finder's tree (runs/2026-11-16-meshholes/paper.md sections 3, 4, 6).
// Each part that cuts one declares its closed outline during its setup (Add); Build cuts every outline into the static mesh host under it (the faces under the outline are
// replaced by faces around it plus a wall per rim edge), after every host exists and before ContactSurfaceHost::FinderBuild. The contact then reads the holes from the mesh
// itself (BarsTriangles, the shared static mesh path): no hole list at contact time. Pure SI (metres); the list is cleared when the physics engine is built.
namespace MeshHoles
{
   // A slot as the drawing needs it (parts/playfieldcutdraw.cpp, slot-hole design section 8). Frame: centre (cx, cy) at the top z0, along (lx, ly), across (px, py), metres (unit: metres
   // per VPX unit). Shape: a the half-width, hs the straight part's half-length, L the tip, U1 where the round end meets the ramp, wJ the half-width there; the ramp's outline
   // u = L - tau^2, w = c tau (tau from 0 at the tip to tauJ at U1); th1 the round end's angle at U1; thick the through-hole's depth
   struct DrawnSlot
   {
      std::string name;
      double unit = 1.0, cx = 0.0, cy = 0.0, z0 = 0.0, lx = 0.0, ly = 1.0, px = -1.0, py = 0.0;
      double a = 0.0, hs = 0.0, L = 0.0, U1 = 0.0, wJ = 0.0, c = 0.0, tauJ = 0.0, th1 = 0.0, thick = 0.0;
      std::vector<double> outlineU;   // the physics outline (meshholes.cpp, VPX units x0 y0 x1 y1 ..., counter-clockwise): the drawn hole (paper R5)
      double W(double u) const;   // the outline's half-width at u past the straight part's end (0 past the tip)
      double C(double u) const;   // the through-hole's half-width (0 past the round end, u >= a)
      // the ramp floor's cross-section at u (U1 < u < L): the circular arc t^2 + (z - zc)^2 = rho^2 through (+-W, 0) and (0, -depth), z up from the top
      void FloorArc(double u, double &rho, double &zc) const;
   };
   // a new load: the outline list and the drawn slots cleared
   void Clear();
   // [Features] mesh_overlap_repair = 1: in every static collidable mesh host (playfield, wall, ramp, primitive), the faces of one host that lie on one plane and overlap in it are made to
   // cover the surface once (a face wholly covered by smaller ones dropped, one partly covered clipped by them), then welded and linked again; log only. Before Build.
   void RepairOverlaps(bool on);
   // a closed outline x0 y0 x1 y1 ... (m); partZ: the part's height (m), the host is the nearest top below it. Walls down by the host's thickness only where the ball can pass through
   // (the joined outline's largest inscribed circle at least the ball's radius, paper section 3 step 4)
   void Add(const std::string &name, const std::vector<double> &xy, double partZ);
   // the table's [PlayfieldCut:<name>] sections (shapes slot, groove, rectangle, circle, triangle; keys in meshholes.cpp Parse) as outlines; file: the table's physics ini
   void AddIniCuts(const std::filesystem::path &file);
   // the [PlayfieldCut:<name>] shapes the drawing draws (the slots today, parts/playfieldcutdraw.cpp), each with its mesh outline (DrawnSlot::outlineU); cleared with the outlines
   std::vector<DrawnSlot> Drawn();
   // a round outline of radius r (m) about x y with edges tangent to the circle: N = max(24, smallest N with r (1 / cos(pi / N) - 1) <= 0.2 mm) (paper section 4)
   std::vector<double> Round(double x, double y, double r);
   // the closed outline x0 y0 x1 y1 ... (m) grown by d (m): each edge moved d outward, parallel, each vertex where its two edges' moved lines cross
   std::vector<double> Offset(const std::vector<double> &xy, double d);
   void Build();
}
