#pragma once
// The playfield's cuts as drawn (playfieldcutdraw.cpp; design runs/2026-11-09-slot-hole/design.md section 8): the hole in the drawn playfield, the walls and the ramp floors in the wood.
// Only for the cuts the physics read for this table (physics/playfieldcut.h, [Features] playfield_cuts): no cut active, nothing drawn or changed.

class Primitive;
class Renderer;
class Material;

namespace PlayfieldCutDraw
{
   // When the playfield's render buffers are made (Primitive::RenderSetup, after the physics read the cuts): its mesh re-cut with each cut's outline left out, the walls and floors made
   void Setup(Primitive *playfield, Renderer *renderer);
   // Right after the playfield is drawn (Primitive::Render; not in its reflections-only pass): the walls and floors, with the playfield's material and the wood texture
   void Draw(const Primitive *playfield, Renderer *renderer, const Material *mat);
   // Primitive::RenderRelease: the table's own mesh back (the table, and a save from the editor, never see the cut mesh), the walls and floors freed
   void Release(Primitive *playfield);
}
