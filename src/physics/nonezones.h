// license:GPLv3+

#pragma once

#include <cstdint>
#include <utility>
#include <vector>

class HitBall;
class Trigger;
class Kicker;
class IFireEvents;

namespace ContactSystem
{
   // None triggers (to-do 11.2, runs/2026-11-10-none-trigger/design.md): an invisible sensing zone, the trigger's polygon as a prism from its surface height to height + HitHeight, for the ball centre.
   // No force on the ball. Hit when a centre enters, Unhit when it leaves, tested on every ball's straight chord of each sub-interval; the events fire at the cycle's end in crossing order.
   // The zone list belongs to the physics engine (PhysicsEngine::m_noneZones): built in Trigger::PhysicSetup, cleared in Trigger::PhysicRelease (design 4m; code in parts/trigger.cpp).
   class NoneZones final
   {
   public:
      void SubStart(double startS);           // after SaveStart: every ball's centre at the start of the sub-interval (startS: its start in the cycle, s); a placed ball or a rechecked zone is judged here
      void SubTest(double startS, double hS); // after CheckOverlaps: each ball's chord from its start to its end of the sub-interval against every enabled zone
      void CycleEnd(double cycleS);           // after the last sub-interval: the cycle's events, in crossing order
      void Placed(const HitBall *ball);       // NotePlaced: a ball created or moved by the script or the driver
      void BallGone(const HitBall *ball);     // Player::DestroyBall: no Unhit, a zone it was in loses its state link when no other ball is listed in it, its queued events go
      void Add(Trigger *trigger, std::vector<Vertex2D> &&poly, float zlow); // Trigger::PhysicSetup: a None trigger's zone (polygon in VPX units, surface height)
      void Remove(const Trigger *trigger);    // Trigger::PhysicRelease
      void Recheck(const Trigger *trigger);   // Enabled went true or HitHeight changed: every ball judged where it is for this zone (the placement rule, design 4d, 4e)
      // A physical kicker's switch (KickerDesign.txt 1f, 7e; Requirements F 6): a switch with travel. It closes (Hit) when the ball's centre is inside the hole's circle and its bottom reaches
      // zClose; once closed it opens (Unhit) only when the bottom rises above zOpen (zOpen > zClose) or the centre leaves the circle. The same chord test, event queue and firing as a None
      // zone, the kicker as owner (BallCntOver reads the ball's m_vpVolObjs, the kicker's held ball follows it). Always on (a physical kicker has no Enabled, F 9). VPX units.
      // A virtual kicker's sensing zone (KickerDesign.txt 1g-1i; Requirements F 10; virt): the ball's centre inside the circle (cx, cy, r) and its height between zClose (the surface) and zOpen
      // (the surface + HitHeight), tested on the chord as a None zone. Entering, the kicker parks the ball if it is empty and enabled (Kicker::Park; else the ball is not listed, no Hit); the
      // parked ball is not judged out (a held ball stays at its place); a kicked ball gets its Unhit when its centre leaves.
      void AddKicker(Kicker *kicker, double cx, double cy, double r, double zClose, double zOpen, bool virt = false);   // Kicker::PhysicSetup
      void Remove(const Kicker *kicker);      // Kicker::PhysicRelease

   private:
      struct Zone
      {
         Trigger *trigger;               // a None trigger's zone (null for a kicker's)
         std::vector<Vertex2D> poly;     // the polygon, VPX units
         float x0, x1, y0, y1;           // its box
         float zlow;                     // the surface height under the trigger
         bool recheck;                   // every ball is judged where it is (the placement rule)
         Kicker *kicker = nullptr;       // a physical kicker's switch (trigger null): the circle (cx, cy, r); the height of the ball's bottom that closes it (zClose) and above which it opens (zOpen), VPX units
         double cx = 0.0, cy = 0.0, r = 0.0, zClose = 0.0, zOpen = 0.0;
         bool virt = false;              // a virtual kicker's sensing zone: zClose .. zOpen is the ball centre's height range
      };
      struct Event
      {
         double t;                       // time in the cycle, s
         Trigger *trigger;               // (null for a kicker's event)
         HitBall *ball;                  // null once the ball is destroyed
         bool hit;
         Kicker *kicker;
      };
      static bool Inside(const Zone &z, double x, double y, double zc, double R, bool held);   // the ball centre in the zone (R: the ball's radius; held: the ball is listed, for a kicker's switch travel)
      static bool Enabled(const Zone &z);
      static IFireEvents *Owner(const Zone &z);   // the zone's object as the ball's m_vpVolObjs lists it (IFireEvents *)
      void Judge(const Zone &z, HitBall *b, bool in, double t, const Vertex3Ds &at);   // at: the ball's centre where the change is judged (a virtual kicker parks it there)
      void ApplyPlacements(double t);
      void Chord(const Zone &z, HitBall *b, const Vertex3Ds &s, const Vertex3Ds &e, double startS, double hS);

      std::vector<Zone> m_zones;
      std::vector<std::pair<HitBall *, Vertex3Ds>> m_start;   // every ball's centre at the start of the sub-interval
      std::vector<const HitBall *> m_placed;                  // the balls the script or the driver put where they are since the last test
      std::vector<Event> m_events, m_firing;                  // this cycle's crossings; the ones being fired (a script may destroy a ball meanwhile)
      std::vector<double> m_cut;                              // a chord's crossing points (kept: no allocation per test)
      uint64_t m_ns = 0, m_cycles = 0, m_ballCycles = 0;      // VPX_CONTACT_TIMING: the zones' cost, logged when the table ends
   };
}
