#pragma once
// Test driver for the harness (docs/unified-physics/TestMatrixDesign.txt section 2). Off unless VPX_TEST_INPUT names an input file (read once at the first cycle of a table).
// The input file lists actions by physics millisecond; they are applied at the start of that cycle, before the physics step, so the same build and the same
// file give the same result whatever the frame rate or the machine load. The driver records the balls every cycle and checks the fixed limits it can check
// (TestMatrixDesign.txt section 1), writing balls.csv, report.txt and done.txt into VPX_TEST_OUT. It also times every cycle and frame (section 3 points 9 and 10):
// cycles.csv, frames.csv and the cost part of report.txt.

#include <cstdint>
#include <intrin.h>
#include <string>

namespace TestDriver
{
   void LoopStart();                // top of a physics loop pass, before the script timers
   void BeginCycle(uint32_t tMs);   // start of a 1 ms physics cycle, before the movers and the physics step
   void EndCycle(uint32_t tMs);     // end of the cycle, after PhysicsSimulateCycle
   void LoopEnd();                  // bottom of the physics loop pass
   void FrameStart();               // a new frame has begun (after the frame profiler's NewFrame): records the frame that ended
   bool Ended();                    // the input file's end was reached and the report is written

   // The timed parts of a cycle (section 3 point 10 a, e, f). Each part's time leaves out the script run inside it.
   enum Part
   {
      kStep, kMagnets, kTurntable,                    // the physics step (PhysicsSimulateCycle), the magnets' and the turntables' cycle steps (turntables: 0 since 2026-10-28, the disc steps as a moving part)
      kContact, kStockSearch, kMoveParts,             // inside the step: the contact system, the stock hit search, the movers' displacement
      kInput, kTimers, kDriver, kNudge, kVelocities, kWatchdog, kGuides, kTrails,   // outside the step
      kParts
   };
   struct Timed   // times one part of the cycle while in scope: { TestDriver::Timed t(TestDriver::kMagnets); Magnets::Step(); }
   {
      explicit Timed(Part p);
      ~Timed();
      Part part;
      int64_t t0;
      uint32_t script0;
   };

   // What the contact system did in this cycle (section 3 point 10 b, c, d), filled by the contact system while the driver times a cycle in detail.
   struct ContactCost
   {
      enum Kind { kSurfaces, kPairs, kBats, kBands, kMovingParts, kKinds };   // the rest of the contact time is "shared" (islands, step choice, integration, write-back)
      uint64_t kindTicks[kKinds] = {};
   };
   ContactCost *Contact();   // null unless this cycle is timed in detail
   const char *BallName(const void *hitBall);   // the driver's name of a ball (its HitBall), null for a ball the driver did not make
   // test tooling (VPX_TEST_IMPULSES=1): an impulse a source gave a ball, SI (J in N s, L in N m s about its centre), summed into impulses.csv (one row per ball and source each recorded cycle)
   void AddImpulse(const void *hitBall, const std::string &source, const double J[3], const double L[3]);
   // the same for a change made outside the group (stock collide and contact, the driver's sets): the ball's momentum and angular momentum now less before
   // (vel0: m_d.m_vel, L0: m_angularmomentum, VPX units, taken before the change), in SI
   void AddBallChange(const void *hitBall, const std::string &source, const float vel0[3], const float L0[3]);
   inline uint64_t Ticks() { return __rdtsc(); }
   struct KindTimer
   {
      KindTimer(ContactCost *cc, ContactCost::Kind kind) : c(cc), k(kind), t0(cc != nullptr ? Ticks() : 0) {}
      ~KindTimer() { if (c != nullptr) c->kindTicks[k] += Ticks() - t0; }
      ContactCost *c;
      ContactCost::Kind k;
      uint64_t t0;
   };
}
