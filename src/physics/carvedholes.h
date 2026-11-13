// license:GPLv3+

#pragma once

#include "physics/rubbercontact.h"

// The playfield's carved holes (KickerDesign.txt 1c; the turntable's CarvedHole made one list for every part that cuts one): round holes in the playfield where its faces do not press
// a ball; the nearest point of what is left is on the hole's rim (RubberContact::ProbeTriangles). Filled by Turntables::Setup (a disc and its gap) and Kicker::PhysicSetup (a physical
// kicker's cutout), cleared when the physics engine is built (ContactSurfaceHost::ClearPlayfieldHoles calls Clear), read by the playfield host's probes (contactsurfacehost.cpp). Metres. The list: carvedholes.cpp.
namespace CarvedHoles
{
   void Add(const RubberContact::CarvedHole &h);
   void Get(const RubberContact::CarvedHole *&holes, int &n);
   void Clear();
}
