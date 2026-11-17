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
   // A joined hole as Build cut it, for the drawing (parts/playfieldcutdraw.cpp: drawn hole = physics hole). World, VPX units (metres / ContactSystem::kMetresPerUnit).
   struct DrawnHole
   {
      std::string name;                    // the joined outlines' names (a + b)
      double unit = 1.0;                   // metres per VPX unit
      std::string hostKind, hostName;      // the static mesh host the physics cut
      std::vector<std::vector<double>> rim;   // the joined outline's boundary loops as Build cut them (x0 y0 x1 y1 ..., counter-clockwise; Pieces' points in order, the same doubles / unit)
      double topZ = 0.0, thick = 0.0;      // the host's top at the hole and the walls' depth (VPX units)
      bool walls = false;                  // walls cut (the joined outline's inscribed circle at least the ball's radius)
      bool draw = true;                    // drawn (false: every joined outline says draw = 0)
      int physWallTris = 0;                // the wall triangles Build emitted (2 per wall quad, end faces included)
      bool onTop = true;                   // the rim lies on the host's top (false: an outline cut after the others into another's floor, a slot in the groove)
      std::vector<float> wood;             // the wall, floor and side triangles Build emitted (before the T-junction pass), 3 vertices each of 6 floats: x y z (VPX units), along (VPX
                                           // units: the rim length on walls, x + y elsewhere), depth below the top (VPX units), kind: 0 a floor (wound to face up), 1 a rim wall (wound to
                                           // face into the hole), 2 a wall under a floor's open edge, 3 a side face (2 and 3: wound either way)
   };
   // a new load: the outline list and the drawn holes cleared; roundMm: [Features] hole_round_mm, Round's vertex excess (physicsini.h)
   void Clear(double roundMm);
   // [Features] mesh_overlap_repair = 1: in every static collidable mesh host (playfield, wall, ramp, primitive), the faces of one host that lie on one plane and overlap in it are made to
   // cover the surface once (a face wholly covered by smaller ones dropped, one partly covered clipped by them), then welded and linked again; log only. Before Build.
   void RepairOverlaps(bool on);
   // a part's hole keys (to-do N11): read from [Triggers] / [Trigger:<name>] (Button, Star), [Kicker:<name>] (type physical), [Turntable:<name>] and [PlayfieldCut:<name>]
   struct HoleKeys
   {
      bool cut = true;              // cut_hole = 1 | 0: 1 (default) the part's hole is cut in the mesh under it; 0 no hole (no Add: the part sits on the whole surface)
      bool draw = true;             // draw = 1 | 0: 1 (default) the hole is drawn as the physics cut it; 0 cut in the physics but not drawn (the drawn top stays whole there)
      double wallDepthMm = -1.0;    // wall_depth_mm: the hole's walls' depth below the top, mm (above 0); absent = the cut host's thickness_mm (contacthost.h). Joined holes: the deepest
   };
   // true when lowerKey is one of HoleKeys' keys (its value read into h); a bad wall_depth_mm (0 or less) is kept as absent, with a notice
   bool HoleKey(HoleKeys &h, const std::string &lowerKey, double v, const std::string &where);
   // a closed outline x0 y0 x1 y1 ... (m); partZ: the part's height (m), the host is the nearest top below it. Walls down by the wall depth only where the ball can pass through
   // (the joined outline's largest inscribed circle at least the ball's radius, paper section 3 step 4)
   void Add(const std::string &name, const std::vector<double> &xy, double partZ, const HoleKeys &keys);
   // the table's [PlayfieldCut:<name>] sections (shapes slot, groove, rectangle, circle, triangle; keys in meshholes.cpp Parse) as outlines; file: the table's physics ini
   void AddIniCuts(const std::filesystem::path &file);
   // every joined hole Build cut on this load (the drawing's record; cleared with the outlines)
   const std::vector<DrawnHole> &Drawn();
   // a round outline of radius r (m) about x y with edges tangent to the circle: N = max(24, smallest N with r (1 / cos(pi / N) - 1) <= hole_round_mm, default 0.2 mm) (paper section 4)
   // n > 0: that many edges instead (a chamfer's inner ring takes its outer ring's N: one face pair per edge between them)
   std::vector<double> Round(double x, double y, double r, int n = 0);
   // the closed outline x0 y0 x1 y1 ... (m) grown by d (m): each edge moved d outward, parallel, each vertex where its two edges' moved lines cross
   std::vector<double> Offset(const std::vector<double> &xy, double d);
   void Build();
}
