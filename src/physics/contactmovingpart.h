// license:GPLv3+
//
// A moving part of the contact system: the engine-facing interface of anything that has its own bodies and surfaces, is stepped inside the islands of the balls that touch it, and
// raises its own events. Gates and spinners (HingeModel::Engine, contacthingepart.h) and the jet bumper (JetBumper::EnginePart, contactbumperpart.h) implement it; the contact system
// (contactsystem.cpp) and the host (ContactPartHost, contactparthost.h) know only this interface. The turntable disc (DiscModel::Engine, contactdiscpart.h) implements it too. A new part (a kicker door, a drop target, a plunger) implements it and nothing else changes.
//
// What the engine does with a part, in the order of a physics cycle:
//   SetGravity       the gravity the part feels (gravity minus the cabinet acceleration, so a nudge acts on it), once per cycle
//   Probe            where is the nearest surface of the part to a ball (distance, normal, curvature): engagement, release, the arrival-time rule
//   AdvanceAlone     a cycle with no ball near: the part moves by itself (a part at rest sleeps)
//   events           what happened (Hit, switch closed / opened, stops, spins), in order; the host fires them on the table part after the step
// Pure SI maths with no VPX types: the same code runs in the offline tests.

#pragma once

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "physics/rubbercontact.h"
#include "physics/contactsurface.h"
#include "physics/groupstep.h"
#include "physics/partrecords.h"

// What the engine keeps of every moving part between its calls (PartRecords: part of its state)
struct MovingPartState
{
   double steppedS = 0.0;           // how long the contact system stepped the part with its balls in this sub-interval (s); the host advances it alone for the rest
   bool steppedThisCycle = false;   // the contact system stepped the part with the balls of an island this cycle (the host's own advance skips it once)
   double eOwnJ = -1.0;             // energy skip test: what it can hold this cycle before a ball touches it (its kinetic energy and drive work, EnergyMove with no feed, plus squashJ), J; -1 = not set (ContactSystem ReachEnergy)
   double squashJ = 0.0;            // the energy reach: what its contacts hold at the cycle's end, J: each ball-part contact's stored energy by the contact law and its compliant stops' squash (GroupStep::Part::StopSquashJ); fed to the balls that can reach it (ContactSystem ReachEnergy)
   bool restTouched = false;        // a ball asleep in a pile touches it (the engine's pile sleep): the part does not sleep while one does (GroupStep::RestWatch)
};

class ContactMovingPart : public PartRecords, public MovingPartState
{
public:
   using V3 = RubberContact::V3;

   struct Event
   {
      enum Kind { Hit, SwitchClosed, SwitchOpened, LimitEOS, LimitBOS, Spin, BandHit, SlingFire, Dropped, Raised } kind;   // (Dropped, Raised: a drop target, droptargetsimple.h)
      int ball = -1;        // Hit, BandHit: the ball's index as the owner numbers them (a drop target's Hit: the ball's own id, Ball::m_id); a band's SwitchClosed / SwitchOpened: the blade's index
      double value = 0.0;   // Hit: normal approach speed (m/s); Limit: angular speed (deg/s); BandHit: the cord's largest deflection under the ball (m); SlingFire: 1 = the script's trigger
   };

   // A saved copy of a part's dynamic state, for rewinding an island (docs/unified-physics/archive/UnifiedContactDesign.txt section 11, point 2): plain bytes, no allocation. Its state records (PartRecords)
   // and its group part's, whole: the values that change while it runs, not its settings or shape.
   struct State
   {
      alignas(8) unsigned char bytes[2048];
      size_t used = 0;
   };

   std::vector<Event> events;

   ContactMovingPart()
   {
      KeepState<MovingPartState>(this);
      KeepQueue(&events);
   }
   virtual ~ContactMovingPart() = default;

   virtual void SetBall(const RubberContact::BallSI &b) = 0;   // the ball of the present call (balls may differ in size and mass)
   virtual void SetGravity(const V3 &g) = 0;
   virtual V3 Origin() const = 0;                              // the part's reference point (world, m): the host's box and the log
   virtual double Reach() const = 0;                           // the swept reach about it (m)
   virtual double MotionRange() const { return 2.0 * Reach(); }   // the farthest any point of its surface can move from where it is, over its whole range of motion (m); every moving point is within Reach of Origin
   virtual RubberContact::TubeProbe Probe(const V3 &c) const = 0;
   virtual void SurfaceVelocity(const V3 &point, V3 &v, V3 &w) const = 0;
   virtual int AdvanceAlone(double dt) = 0;
   // Save the dynamic state / put it back exactly: its records and its group part's, copied whole (one shared code for every part). Restoring cuts each event queue back to what it held at the save
   // (the replay makes the later ones again) and rebuilds what is derived from the state.
   void SaveState(State &out)
   {
      GroupStep::Part *const g = GroupPart();
      const size_t a = StateBytes(), b = g != nullptr ? g->StateBytes() : 0;
      assert(a + b <= sizeof out.bytes);
      SaveRecords(out.bytes);
      if (g != nullptr)
         g->SaveRecords(out.bytes + a);
      out.used = a + b;
   }
   void RestoreState(const State &in)
   {
      GroupStep::Part *const g = GroupPart();
      RestoreRecords(in.bytes);
      if (g != nullptr)
      {
         g->RestoreRecords(in.bytes + StateBytes());
         g->AfterRestore();
      }
   }
   // The fastest any point of the part's surfaces can move ALONG ITS OWN NORMAL (toward a ball) within the next cycle of cycleS seconds (m/s): its present speed plus what its own drive and springs
   // can add in a cycle. Balls are not counted (their contact is the contact). The engine adds it to the ball's closing speed when it decides to engage a ball and when it decides to let one go
   // (UnifiedContactDesign.txt 11 point 1, 12.1 d). A part may answer its whole point speed (an upper bound); a disc turning about its own axis answers 0 (its surfaces slide along themselves).
   virtual double SurfaceSpeedBound(double cycleS) const = 0;
   // the energy reach (design_reach.md sections 1, 4b; rubbercontact.h EnergyReach): the most any point of the part's surface can move relative to a point c over a cycle of T
   // seconds (m), its turn or travel bounded by energy: its own motion, its drives (gravity, springs, blades, cams, coil) and E, the kinetic energy (J) of the balls that can reach it
   // this cycle. work (if not null): what the part can hand to the balls in the cycle beyond E, its kinetic energy now plus its drives' work over that move (J). Each part its own law
   virtual double EnergyMove(const V3 &c, double E, double T, double *work) const = 0;
   // the one-pass energy bound (runs/2026-11-20-energy-bound/design.md, AGREED 2026-11-20 items 1-6, 4b, 5b): the most the part can hand to the balls over
   // a cycle of T seconds whatever they feed it (J), from its state now and its built data, no move: its energy in the group's book now (GroupStep::Part::Energy: its kinetic energy,
   // its springs, blades, cams and weight) less the least that book's stored energy takes over the poses it can reach between its stops (item 4b: one source, no per-spring formula),
   // + its drives over T (a coil's field energy now + V^2 T / R, a motor's largest power x T). Past its stops under the force that presses it there (weight, preload): a linear stop
   // adds aLin = F sqrt(2 / k) (J^0.5, item 5: E = E0 + aLin sqrt(E)), a Hertz stop (k d^1.5) aHertz = F (5 / (2 k))^0.4 (J^0.6, item 5b: E = E0 + aHertz E^0.4); the caller closes
   // them in one step (ContactSystem OnePassE). The squash its contacts and compliant stops hold now is squashJ, not here
   virtual double EnergyCap(double T, double &aLin, double &aHertz) const = 0;
   // between cycles (the host's advance, after the cycle): what starts at the next cycle's start (a jet's pending pulse: a jet's pulse at a cycle's start); cycleS: that cycle's length
   virtual void NextCycleStart(double /*cycleS*/) {}
   // the energy reach: a lower bound on the distance from c to the part at any pose it can take in the cycle of T seconds with E (J) fed to it (m; its region, design_reach.md
   // section 3 in its bounding form; a part whose region does not depend on its feed ignores E); -INFINITY = none (the move above alone). Grouping takes the larger of this less R
   // and the gap now less EnergyMove
   virtual double CycleRegionDist(const V3 &, double /*T*/, double /*E*/) const { return -INFINITY; }
   // the part cull (runs/2026-11-20-part-cull/design.md 1c, 2R, 2J): CullCycle once per cycle, before the grouping, with E the table's whole energy this cycle (J:
   // every ball's, every stored and every part's work, an upper bound on any group's feed) and the cycle T (s); then CullDist(c) is a lower bound on the distance (m) from c to any
   // touchable point of the part at any pose it can take in the cycle (the ball's radius not included). -INFINITY = no region: a candidate for every ball. Default: the
   // grouping's own region at the table's energy (CycleRegionDist: the bat, the hinge parts, the disc; the jet and the target precompute the same region once per cycle)
   virtual void CullCycle(const double E, const double T) { cullE = E; cullT = T; }
   virtual double CullDist(const V3 &c) const { return CycleRegionDist(c, cullT, cullE); }
   virtual std::string RegionNote(double /*E*/, double /*T*/) const { return std::string(); }   // (the log, once per table start: how its region is built; empty = nothing to tell)
   double cullT = 0.0, cullE = 0.0;   // (the cycle and the table's energy CullCycle was set for, s and J; not part of the state: set every cycle)
   // the part's state as a short text for the log (an engagement of a ball with a part that is moving is the case to look at)
   virtual std::string StateNote() const { return std::string(); }
   // A part on the group step (docs/unified-physics/GroupStepDesign.txt 18): the engine moves it with the balls in one solve and holds its shapes. Null: the part is on the old loop until it moves over. No part implements these yet.
   virtual GroupStep::Part *GroupPart() { return nullptr; }
   virtual int GroupSurfaces() const { return 1; }                                                   // its surfaces a ball can meet: one shape per ball and surface, each with its own contact and state
   virtual void GroupProbe(int, const V3 &, double, const double *, GroupStep::Hit &) const {}       // a surface for a ball centre and radius, its coordinates moved by trial (null = where they are); fills the hit's columns
   virtual bool GroupCycleNear(int, const V3 &, double, double, std::vector<int> &, int &) const { return false; }              // the part finder: the cycle's list (GroupStep::Shape::CycleNear)
   virtual int GroupCycleMisses(int, const V3 &, double, double, const std::vector<int> &) const { return 0; }                // (GroupStep::Shape::CycleMisses, check mode)
   // the part finder, the overlap check of free balls (ContactSystem CheckOverlaps): a lower bound on the distance from c to the part where it is now; each part its bounding
   // radius from its built pieces (group A: not Reach(), which carries margins). Here: none, always probed
   virtual double RegionDist(const V3 &) const { return -INFINITY; }
   virtual const RubberContact::RubberParams *GroupLaw(int, double &) const { return nullptr; }      // the surface's material and layer; null = the host's
   // the path test (GroupStep::Shape::PathReach): a surface of the part that never moves (a bumper's body, a disc turning about its own axis), the ball's centre moving straight
   // from s to s + d: 1 = it can touch it on the way, 0 = no (exact), -1 = no exact test (a surface that moves with the part: build 3, piece 6), FarStillFar probes as before
   virtual int GroupPathReach(int, const V3 &, const V3 &, double) const { return -1; }
   // the part finder, a surface that moves with the part (group A; GroupStep::Shape::PathMotion): the most any point of it moves relative to P over the step, the ball's centre
   // moving dLen and the part by trial (the solve's end move; null = no move); its probe (GroupProbe) must give the exact distance at any range. -1 = none: the two probes
   virtual double GroupPathMotion(int, const V3 &, double, const double *) const { return -1.0; }
   virtual double GroupPathSlope(int) const { return 1.0; }   // its probe's largest slope (GroupStep::Shape::PathSlope), built at load
   virtual double GroupPathCoverSag(int, int &) const { return -1.0; }   // the cover's largest sagitta and kind (the self-check, VPX_FINDER_CHECK; GroupStep::Shape::PathCoverSag); -1 = none
   virtual void CycleBegin() {}                                                                      // a group holding it is being built this cycle
   virtual void TouchBegan(int, int, int, double, const V3 &) {}                                     // this surface, this ball of the group (its index; its own id, Ball::m_id), closing at vn (m/s), the touch's normal (toward the ball): a touch began. A part's Hit is what it decides
   virtual bool Sleeping() const { return false; }                                                   // no ball near and nothing to step: no group this cycle
   virtual bool RestNumbers(double &, double &) const { return false; }                             // pile sleep's log: what AtRest is decided on, the move over the last cycle in a group and the surface speed (m, m/s); false = not given
   virtual bool AtRest() const { return false; }                                                     // nothing driven and nothing moved over the last whole cycle: it stays out of a ball's solve while no ball can reach it (a part that does not say so is always in)
   virtual GroupStep::RestWatch *Rest() { return nullptr; }                                          // its rest watch (GroupStep::RestWatch): the group restarts it where it is built, with RestPose
   virtual GroupStep::RestWatch::Pose RestPose() const { return {}; }                                // its one rest pose, alone and in a group: its coordinates and internal states
};
