// license:GPLv3+
//
// A moving part of the contact system: the engine-facing interface of anything that has its own bodies and surfaces, is stepped inside the islands of the balls that touch it, and
// raises its own events. Gates and spinners (HingeModel::Engine, contacthingepart.h) and the jet bumper (JetBumper::EnginePart, contactbumperpart.h) implement it; the contact system
// (contactsystem.cpp) and the host (ContactPartHost, contactparthost.h) know only this interface. The turntable disc (DiscModel::Engine, contactdiscpart.h) implements it too. A new part (a kicker door, a drop target, a plunger) implements it and nothing else changes.
//
// What the engine does with a part, in the order of a physics cycle:
//   SetGravity       the gravity the part feels (gravity minus the cabinet acceleration, so a nudge acts on it), once per cycle
//   Probe            where is the nearest surface of the part to a ball (distance, normal, curvature): engagement, release, the arrival-time rule
//   BallStep         the longest substep this ball's contact with the part allows (reduced mass of ball and part, contact stiffness, arrival time)
//   ContactItem      one substep of one ball's contact: the force and torque on the ball; the reaction is kept by the part and acts in StepInIsland
//   StepInIsland     one island substep: the part moves with the reactions of all the balls summed, in its own finer steps if it needs them
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

// Scale on the impact rule of the contact step (omega * h <= 0.2 in the contact system, 0.1 in a part's BallStep): VPX_OMEGAH_SCALE, default 1. A diagnostic knob for the step-size study of 2026-10-27.
inline double ContactOmegaHScale()
{
   static const double s = []() { const char *const e = getenv("VPX_OMEGAH_SCALE"); return e != nullptr && *e ? atof(e) : 1.0; }();
   return s;
}

namespace GroupStep { class Part; struct Hit; }

class ContactMovingPart
{
public:
   using V3 = RubberContact::V3;

   struct Event
   {
      enum Kind { Hit, SwitchClosed, SwitchOpened, LimitEOS, LimitBOS, Spin, BandHit, SlingFire } kind;
      int ball = -1;        // Hit, BandHit: the ball's index as the owner numbers them; a band's SwitchClosed / SwitchOpened: the blade's index
      double value = 0.0;   // Hit: normal approach speed (m/s); Limit: angular speed (deg/s); BandHit: the cord's largest deflection under the ball (m); SlingFire: 1 = the script's trigger
   };

   // A saved copy of a part's dynamic state, for rewinding an island (docs/unified-physics/archive/UnifiedContactDesign.txt section 11, point 2): plain bytes, no allocation. A part packs the values that change while it
   // runs (position, speed, switch and drive state, stop contacts, held forces), not its settings or shape.
   struct State
   {
      alignas(8) unsigned char bytes[2048];
      size_t used = 0;
   };
   template <class T> static void Put(State &s, const T &v)
   {
      assert(s.used + sizeof(T) <= sizeof(s.bytes));
      std::memcpy(s.bytes + s.used, &v, sizeof(T));
      s.used += sizeof(T);
   }
   struct StateReader
   {
      const State &s;
      size_t pos = 0;
      explicit StateReader(const State &st) : s(st) {}
      template <class T> void Get(T &v)
      {
         assert(pos + sizeof(T) <= s.used);
         std::memcpy(&v, s.bytes + pos, sizeof(T));
         pos += sizeof(T);
      }
   };

   std::vector<Event> events;
   double steppedS = 0.0;           // how long the contact system stepped the part with its balls in this sub-interval (s); the host advances it alone for the rest
   bool steppedThisCycle = false;   // the contact system stepped the part with the balls of an island this cycle (the host's own advance skips it once)

   virtual ~ContactMovingPart() = default;

   virtual void SetBall(const RubberContact::BallSI &b) = 0;   // the ball of the present call (balls may differ in size and mass)
   virtual void SetGravity(const V3 &g) = 0;
   virtual V3 Origin() const = 0;                              // the part's reference point (world, m): the host's box and the log
   virtual double Reach() const = 0;                           // the swept reach about it (m)
   virtual double MotionRange() const { return 2.0 * Reach(); }   // the farthest any point of its surface can move from where it is, over its whole range of motion (m); every moving point is within Reach of Origin
   virtual RubberContact::TubeProbe Probe(const V3 &c) const = 0;
   virtual RubberContact::ContactResult ContactItem(int ballIndex, SurfaceContact &it, const V3 &c, const V3 &v, const V3 &w, double h) = 0;
   virtual void SurfaceVelocity(const V3 &point, V3 &v, V3 &w) const = 0;
   virtual double BallStep(const V3 &c, const V3 &v, double hMax) const = 0;
   virtual void StepInIsland(double h) = 0;
   virtual int AdvanceAlone(double dt) = 0;
   // Save the dynamic state / put it back exactly. Restoring discards the events the part queued since (the replay makes them again) and rebuilds what is derived from the state.
   virtual void SaveState(State &out) const = 0;
   virtual void RestoreState(const State &in) = 0;
   // The fastest any point of the part's surfaces can move ALONG ITS OWN NORMAL (toward a ball) within the next cycle of cycleS seconds (m/s): its present speed plus what its own drive and springs
   // can add in a cycle. Balls are not counted (their contact is the contact). The engine adds it to the ball's closing speed when it decides to engage a ball and when it decides to let one go
   // (UnifiedContactDesign.txt 11 point 1, 12.1 d). A part may answer its whole point speed (an upper bound); a disc turning about its own axis answers 0 (its surfaces slide along themselves).
   virtual double SurfaceSpeedBound(double cycleS) const = 0;
   // the part's state as a short text for the log (an engagement of a ball with a part that is moving is the case to look at)
   virtual std::string StateNote() const { return std::string(); }
   // A part on the group step (docs/unified-physics/GroupStepDesign.txt 18): the engine moves it with the balls in one solve and holds its shapes. Null: the part is on the old loop until it moves over. No part implements these yet.
   virtual GroupStep::Part *GroupPart() { return nullptr; }
   virtual int GroupSurfaces() const { return 1; }                                                   // its surfaces a ball can meet: one shape per ball and surface, each with its own contact and state
   virtual void GroupProbe(int, const V3 &, double, const double *, GroupStep::Hit &) const {}       // a surface for a ball centre and radius, its coordinates moved by trial (null = where they are); fills the hit's columns
   virtual bool GroupNearPieces(int, const V3 &, double, double, double, double, std::vector<int> &) const { return false; }   // a surface of pieces: the step's list (GroupStep::Shape::NearPieces)
   virtual bool GroupNearHolds(int, const V3 &, double, double, const std::vector<int> &) const { return true; }                // (GroupStep::Shape::NearHolds)
   virtual const RubberContact::RubberParams *GroupLaw(int, double &) const { return nullptr; }      // the surface's material and layer; null = the host's
   virtual void CycleBegin() {}                                                                      // a group holding it is being built this cycle
   virtual void TouchBegan(int, int, double) {}                                                      // this surface, this ball of the group, closing at vn (m/s): a touch began. A part's Hit is what it decides
   virtual bool Sleeping() const { return false; }                                                   // no ball near and nothing to step: no group this cycle
   virtual bool RestNumbers(double &, double &) const { return false; }                             // pile sleep's log: what AtRest is decided on, the move over the last cycle in a group and the surface speed (m, m/s); false = not given
   virtual bool AtRest() const { return false; }                                                     // nothing driven and nothing moved over the last whole cycle: it stays out of a ball's solve while no ball can reach it (a part that does not say so is always in)
};
