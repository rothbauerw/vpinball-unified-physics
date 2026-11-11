// license:GPLv3+

#pragma once

#include <cstdint>
#include <utility>
#include <vector>

class HitBall;
class Trigger;

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

   private:
      struct Zone
      {
         Trigger *trigger;
         std::vector<Vertex2D> poly;     // the polygon, VPX units
         float x0, x1, y0, y1;           // its box
         float zlow;                     // the surface height under the trigger
         bool recheck;                   // every ball is judged where it is (the placement rule)
      };
      struct Event
      {
         double t;                       // time in the cycle, s
         Trigger *trigger;
         HitBall *ball;                  // null once the ball is destroyed
         bool hit;
      };
      static bool Inside(const Zone &z, double x, double y, double zc);   // the ball centre in the zone
      void Judge(const Zone &z, HitBall *b, bool in, double t);
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
