// The cord in the engine: one band's model (cordfixed.h) and what the contact system keeps for it between cycles (docs/unified-physics/GroupStepDesign.txt 8, 11)
#pragma once

#include "physics/cordfixed.h"
#include <memory>
#include <utility>
#include <vector>

class HitBall;

struct CordBand
{
   std::unique_ptr<CordFixed::Part> fixed;   // the fixed-height cord (cordfixed.h)
   CordFixed::PlateTube fixedTube;           // what a ball meets on its plate
   static constexpr int kMaxBalls = 8;       // balls in one island with a cord; the others are left out that cycle
   bool PlateQuiet() const { return fixed->PlateQuiet(); }
   bool PlateAtRest() const { return fixed->PlateAtRest(); }   // at its rest angle and still (CordFixed::Part::PlateAtRest)
   bool OnCord(const int i) const { return fixed->OnCord(i); }
   double Gap(const int slot) const { return fixed->CordGap(slot); }
   double Closing(const int i, const GroupStep::V3 &v) const { return fixed->ClosingOnCord(i, v); }
   bool attached = false;          // the part is in a group's solve (AttachGroup done)
   std::vector<HitBall *> slots;   // the group's balls in the part's order, last cycle
   bool warnedHeld = false;
   int heldCycles = 0;             // cycles in a row whose step could not be solved: the balls are held where they are for the rest of each
   std::vector<HitBall *> onCord;  // the balls on the cord at the end of the last cycle
   std::vector<std::pair<HitBall *, float>> hits;   // a ball that came onto the cord this cycle, and its closing speed (VPX units): the script's Hit event, after the step
};
