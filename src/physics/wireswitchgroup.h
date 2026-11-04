// license:GPLv3+
//
// A rollover wire switch as a part of the group step (runs/2026-11-03-switch-requirements/design-wire.md, topic 1): a rigid steel wire on the lever of a micro switch, one angle
// about the switch's hinge below the playfield (hingegroup.h, bodypart.h): a gate turned on its side.
//   spring   the micro switch's lever spring, a preloaded torsion spring tau = -k (theta - a), a < 0 its free angle beyond the rest stop (the hinge part's one-sided blade);
//            k and a set by the vertical force at the peak: spring_n at first touch and spring_full_n at full travel, the wire's own weight included
//   stops    rest (wire up, theta 0) and full travel (the peak stop_below_mm below the playfield): BodyStop, the Hertz steel stop with viscoelastic restitution, as gates
//   pivot    Coulomb pin friction, gravity with the cabinet's acceleration, as gates
//   switch   closes when the peak has travelled close_mm; opens when it rises back release_mm past that point (the micro switch's differential travel). Hit on close, Unhit on open
// Frame (hingegroup.h): the axis runs along the span (the lane), N points sideways from the wire to the hinge, so a positive theta presses the wire down from either end and a push along the span has no torque. Geometry INFERRED from a part photo (lead 2026-11-04): the 3/4 in leg goes down to the lever, which runs hinge_offset_mm sideways to the pivot, hinge_depth_mm below the playfield. Shape (this piece): the span the ball rolls over is one capsule; the legs
// carry mass only. The mesh centreline (wireshape.h) replaces the capsule when the part is wired. Pure SI.

#pragma once

#include <algorithm>
#include <cmath>

#include "physics/hingegroup.h"

namespace WireSwitchGroup
{
   using GroupStep::V3;
   using RubberContact::BallSI;

   struct Settings
   {
      double peakMm = 6.0;            // wire top above the playfield at rest
      double stopBelowMm = 1.0;       // full travel: the peak this far below the playfield
      double closeMm = 3.0;           // the switch closes at this travel of the peak
      double releaseMm = 0.5;         // and opens when the peak rises back this far past the close point
      double springN = 0.02;          // net vertical force at the peak at first touch, the wire's weight included
      double springFullN = 0.03;      // at full travel
      double wireDMm = 1.4;           // 0: the mesh's tube + WireThickness (when wired)
      double spanMm = 25.4, upLegMm = 19.05, downLegMm = 12.7;   // the pictured wire: 1 in span, 3/4 in leg up from the lever, 1/2 in leg down
      double hingeOffsetMm = 15.0;    // UNVERIFIED (no machine): the lever, sideways from the wire's plane to the pivot
      double hingeDepthMm = 13.0;     // UNVERIFIED: hinge below the playfield (the 3/4 in leg minus the 6 mm above)
      double stopArmMm = 5.0, stopE = 0.5, stopEdgeMm = 0.4;   // the lever on the switch body's stops, as the gate's
      double muPivot = 0.3;
      double density = 7850.0;
      double steelEGpa = 200.0;
      double elasticity = 0.3, friction = 0.02;   // ball contact: the wire's restitution at 1 m/s, friction
   };

   class Part : public HingeGroup::Part
   {
   public:
      Settings set;

      // surfacePoint: on the playfield under the span's middle; toHinge: horizontal, across the span, the way from the wire to the hinge (the span runs along E = toHinge x up)
      void Build(const Settings &s, const V3 &surfacePoint, const V3 &toHinge, const BallSI &b, const V3 &gravity)
      {
         set = s;
         const double mm = 1e-3, deg = GroupStep::kPi / 180.0;
         const V3 dir = toHinge / toHinge.Length();
         Begin(surfacePoint + dir * (s.hingeOffsetMm * mm) - V3(0.0, 0.0, s.hingeDepthMm * mm), dir, gravity);
         muPivot = s.muPivot;
         rho = 0.5 * s.wireDMm * mm;
         rPin = rho;
         const double lam = s.density * GroupStep::kPi * rho * rho;
         // body coordinates at theta 0: x along the span, y = -height above the hinge, z = distance along N from the hinge (the wire is beside it, z < 0)
         const double zc = -s.hingeOffsetMm * mm, yTop = -(s.hingeDepthMm + s.peakMm) * mm, yc = yTop + rho, hs = 0.5 * s.spanMm * mm;
         const V3 pNear(hs, yc, zc), pFar(-hs, yc, zc);
         AddWire(pNear, pFar, lam, true);
         AddWire(pNear + V3(0.0, s.upLegMm * mm, 0.0), pNear, lam, false);
         AddWire(pFar, pFar + V3(0.0, s.downLegMm * mm, 0.0), lam, false);
         peakY = yTop;
         peakZ = zc;
         const double reachM = std::sqrt(hs * hs + zc * zc + std::max(yTop * yTop, (yc + s.upLegMm * mm) * (yc + s.upLegMm * mm))) + rho;
         FinishShape(reachM);

         // the travel of the peak: rest stop at theta 0, full travel where it has dropped peak + stop_below, the switch's make and break
         q = rate = 0.0;
         thFull = AngleAtDrop((s.peakMm + s.stopBelowMm) * mm);
         switchMake = AngleAtDrop(s.closeMm * mm);
         switchBreak = AngleAtDrop(std::max(s.closeMm - s.releaseMm, 0.0) * mm);
         const auto stop = [&](const double at, const int side) {
            GroupStep::BodyStop st;
            st.at = at;
            st.side = side;
            st.lever = s.stopArmMm * mm;
            st.Setup(mp.I, std::sqrt(rho * s.stopEdgeMm * mm), s.stopE, 0.1);
            return st;
         };
         stopRest = PlaceStop(stop(0.0, -1));
         stopFull = PlaceStop(stop(thFull, 1));

         // the spring from the net vertical force at the peak: k (theta - a) - tau_g(theta) = F(theta) * arm(theta) at theta 0 and at full travel (tau_g: gravity's torque toward +theta)
         kBlade = aBlade = 0.0;
         const double t0 = s.springN * PeakArm(0.0) + StoredForce(0.0), tf = s.springFullN * PeakArm(thFull) + StoredForce(thFull);
         if (tf > t0 && t0 > 0.0)
         {
            kBlade = (tf - t0) / thFull;
            aBlade = -t0 / kBlade;
         }
         SetContactLaw(s.elasticity, s.friction, b);
         CycleBegin();
      }
      double PeakDrop(const double th) const { return PeakHeight(0.0) - PeakHeight(th); }   // how far the peak has travelled down, m
      double FullAngle() const { return thFull; }
      double MakeAngle() const { return switchMake; }
      double BreakAngle() const { return switchBreak; }

      // The trigger's Hit is the switch closing (swLog), not a touch
      bool NearStart(double) const override { return false; }
      void UpdateSwitches() override
      {
         if (!switchClosed && q >= switchMake)
         {
            switchClosed = true;
            swLog.push_back(1);
         }
         else if (switchClosed && q < switchBreak)
         {
            switchClosed = false;
            swLog.push_back(0);
         }
      }
      // bend out of the wire's plane at distance y from the axis: one wire, k = 3 E I / y^3, I = pi d^4 / 64
      double BendKb(const double y) const override
      {
         const double d = 2.0 * rho, I = GroupStep::kPi * d * d * d * d / 64.0;
         return 3.0 * set.steelEGpa * 1e9 * I / (y * y * y);
      }

   private:
      // the peak's height above the hinge at theta (Down(th) . Z = -cos, Front(th) . Z = sin) and its vertical lever, |d height / d theta|
      double PeakHeight(const double th) const { return -peakY * std::cos(th) + peakZ * std::sin(th); }
      double PeakArm(const double th) const { return std::fabs(peakY * std::sin(th) + peakZ * std::cos(th)); }
      // the angle where the peak has dropped d: -y cos(th) + z sin(th) = h0 - d, the root nearest 0 on the pressing side
      double AngleAtDrop(const double d) const
      {
         const double r = std::sqrt(peakY * peakY + peakZ * peakZ), phi = std::atan2(peakZ, -peakY), c = (PeakHeight(0.0) - d) / r;   // height = r cos(th - phi)
         return phi + std::acos(std::max(-1.0, std::min(1.0, c)));
      }
      double peakY = 0.0, peakZ = 0.0, thFull = 0.0, switchMake = 0.0, switchBreak = 0.0;
      int stopRest = -1, stopFull = -1;
   };
}
