#pragma once
// The mesh holes as drawn (playfieldcutdraw.cpp; design runs/2026-11-09-slot-hole/design.md section 8): each hole in the drawn top of the host the physics cut it in (the playfield, any
// primitive, a wall's top), the walls and the floors in the wood. Ramps are never cut (meshholes.cpp Build).
// Only for the holes the physics cut for this table (physics/meshholes.h): no hole in a host, nothing drawn or changed there.

#include <vector>

class Primitive;
class Surface;
class Renderer;
class Material;

namespace PlayfieldCutDraw
{
   // When a primitive's render buffers are made (Primitive::RenderSetup, after the physics cut the holes): its mesh re-cut with the rim loops of each hole cut in it left out (the
   // playfield host: the playfield primitive; a primitive host: the primitive of that name), the walls and floors made
   void Setup(Primitive *prim, Renderer *renderer);
   // When a wall's render buffers are made (Surface::RenderSetup): its top (top: its vertices, world, VPX units; idx: its triangles) re-cut with the rim loops of each hole cut in the
   // wall left out, the walls and floors made. true: top and idx replaced by the cut top
   bool SetupWallTop(const Surface *wall, Renderer *renderer, std::vector<Vertex3D_NoTex2> &top, std::vector<unsigned int> &idx);
   // Right after the part's top is drawn (Primitive::Render, not in its reflections-only pass; Surface::RenderWallsAtHeight, not dropped, not in the reflection pass): the walls and
   // floors, with the part's material and the wood texture
   void Draw(const void *part, Renderer *renderer, const Material *mat, const Vertex3Ds &center, float depthBias);
   // Primitive::RenderRelease, Surface::RenderRelease: a primitive's own mesh back (the table, and a save from the editor, never see the cut mesh), the walls and floors freed
   void Release(Primitive *prim);
   void Release(const Surface *wall);
}
