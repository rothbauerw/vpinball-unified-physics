// The cord in the engine: one band's model (bandcord.h) and what the contact system keeps for it between cycles (docs/GroupStepDesign.txt 8, 11)
#pragma once

#include "physics/bandcord.h"
#include <utility>
#include <vector>

class HitBall;

struct CordBand
{
   BandCord::Part part;
   bool attached = false;          // the part is in a group's solve (AttachGroup done)
   std::vector<HitBall *> slots;   // the group's balls in the part's order, last cycle
   bool warnedFloor = false, warnedHeld = false;
   int crossingsShown = 0;         // the part's crossings already told on screen and in the log
   int heldCycles = 0;             // cycles in a row whose step could not be solved: the balls are held where they are for the rest of each
   std::vector<HitBall *> onCord;  // the balls on the cord at the end of the last cycle
   std::vector<std::pair<HitBall *, float>> hits;   // a ball that came onto the cord this cycle, and its closing speed (VPX units): the script's Hit event, after the step
};
