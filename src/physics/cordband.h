// The cord in the engine: one band's model (bandcord.h) and what the contact system keeps for it between cycles (docs/GroupStepDesign.txt 8, 11)
#pragma once

#include "physics/bandcord.h"
#include <vector>

class HitBall;

struct CordBand
{
   BandCord::Part part;
   bool attached = false;          // the part is in a group's solve (AttachGroup done)
   std::vector<HitBall *> slots;   // the group's balls in the part's order, last cycle
   bool warnedFloor = false, warnedThrough = false;
};
