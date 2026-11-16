// license:GPLv3+

#pragma once

#include <string>
#include <vector>

// Holes cut in a static mesh host's upward faces at load, in place, before the finder's tree ([Features] mesh_holes; runs/2026-11-16-meshholes/paper.md sections 3, 4, 6).
// Each part that cuts one declares its closed outline during its setup (Add); Build cuts every outline into the static mesh host under it (the faces under the outline are
// replaced by faces around it plus a wall per rim edge), after every host exists and before ContactSurfaceHost::FinderBuild. The contact then reads the holes from the mesh
// itself (BarsTriangles, the shared static mesh path): no hole list at contact time. Pure SI (metres); the list is cleared when the physics engine is built.
namespace MeshHoles
{
   // on: [Features] mesh_holes = 1 and the table has no running [PlayfieldCut:<name>] section (until piece 3 the math cuts and the mesh holes do not mix)
   void Clear(bool on);
   bool On();
   // a closed counter-clockwise outline x0 y0 x1 y1 ... (m); through: walls down by the host's thickness; partZ: the part's height (m), the host is the nearest top below it
   void Add(const std::string &name, const std::vector<double> &xy, bool through, double partZ);
   // a round outline of radius r (m) about x y with edges tangent to the circle: N = max(24, smallest N with r (1 / cos(pi / N) - 1) <= 0.2 mm) (paper section 4)
   std::vector<double> Round(double x, double y, double r);
   // the closed outline x0 y0 x1 y1 ... (m) grown by d (m): each edge moved d outward, parallel, each vertex where its two edges' moved lines cross
   std::vector<double> Offset(const std::vector<double> &xy, double d);
   void Build();
}
