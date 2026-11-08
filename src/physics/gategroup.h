// license:GPLv3+
//
// A gate as a part of the group step: the four real gate types (docs/unified-physics/GateSpinnerDesign.txt section 4; runs/2026-10-31-gate-spinner-disc) on one axis (hingegroup.h, bodypart.h).
//   wire W          pins at both ends; legs from the pin ends down to two bottoms at +-width/4, a middle V that rises to apex * depth. Each side of the V is its own contact
//   wire rectangle  legs from the pin ends down to a bottom bar
//   plate           trapezoid, top width 27 mm, bottom 18 mm, 20 mm long, hanging from a pin wire across its top
//   long plate      trapezoid 27 / 13.5 mm, 27 mm long, plus a lip at the bottom bent back by lip_bend
// Stops (steel on steel): one-way: a reverse stop at the rest angle and a forward stop at the open angle; two-way: stops at +- the open angle, the gate hangs along gravity. Hold: a driven stop that pushes it to the hold angle and goes back.
// Switch: closes when the gate is open far enough for the ball to pass, lead_deg before that; opens hyst_deg lower. Optional blade spring. Pure SI.

#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "physics/hingegroup.h"

namespace GateGroup
{
   using GroupStep::V3;
   using RubberContact::BallSI;

   enum class Type { WireW, WireRectangle, Plate, LongPlate };

   struct Settings
   {
      Type type = Type::WireW;
      double widthMm = 38.0;          // wire gates: the part's length
      double pivotHeightMm = 31.75;   // axis above the surface (1.25 in)
      double hitOffsetMm = 0.0;       // wire gates: where the lowest wire meets a rolling ball, above its centre
      double depthMm = 0.0;           // wire gates: > 0 overrides the depth from the pivot height
      double wireDMm = 1.2;
      double wApex = 0.45;            // wire W: middle apex rise as a fraction of the depth
      double bottomFrac = 1.0;        // wire rectangle: bottom width / top width
      double plateTopMm = 27.0, plateBottomMm = 18.0, plateLenMm = 20.0;
      double longTopMm = 27.0, longBottomMm = 13.5, longLenMm = 27.0;
      double plateThickMm = 0.65;
      double lipMm = 2.75, lipBendDeg = 15.0;
      bool oneWay = true;
      double restAngleDeg = 15.0;     // one-way rest (on the reverse stop)
      double openStopDeg = 90.0;
      double stopArmMm = 5.0, stopE = 0.5, bracketEdgeMm = 0.4;   // stop lever, restitution at 0.1 m/s, radius of the bracket edge it touches
      double muPivot = 0.3;
      double density = 7850.0;
      bool hasSwitch = true;          // ([Gates] switch = none)
      double switchLeadDeg = 5.0, switchHystDeg = 3.0, switchForceN = 0.0;
      double holdAngleDeg = 70.0, holdSpeedDps = 900.0;
      double steelEGpa = 200.0;       // wire and plate Young's modulus for the bending
      double hitWindowDeg = 5.0;      // a Hit only fires while the gate is within this of its start angle (the strike that opens it, user 2026-10-15)
      double elasticity = 0.3, friction = 0.02;   // ball contact: the part's restitution at 1 m/s, friction
   };

   class Part : public HingeGroup::Part
   {
   public:
      Settings set;

      // surfacePoint: on the surface under the axis; front: the way a ball passes from the front (horizontal)
      void Build(const Settings &s, const V3 &surfacePoint, const V3 &front, const BallSI &b, const V3 &gravity)
      {
         set = s;
         Begin(surfacePoint + V3(0.0, 0.0, s.pivotHeightMm * 1e-3), front, gravity);
         muPivot = s.muPivot;
         rPin = 0.5 * s.wireDMm * 1e-3;
         const double mm = 1e-3, lam = s.density * GroupStep::kPi * 0.25 * s.wireDMm * s.wireDMm * mm * mm;
         rho = 0.5 * s.wireDMm * mm;
         halfT = 0.5 * s.plateThickMm * mm;
         const double arm = s.stopArmMm * mm, sigma = s.density * s.plateThickMm * mm;
         const auto wire = [&](const V3 &p, const V3 &u, const bool touchable) { AddWire(p, u, lam, touchable); };
         double halfTop, depthLow;   // half width at the axis, depth of the lowest point below the axis
         if (s.type == Type::WireW || s.type == Type::WireRectangle)
         {
            const double hw = 0.5 * s.widthMm * mm, D = s.depthMm > 0.0 ? s.depthMm * mm : s.pivotHeightMm * mm - b.R - s.hitOffsetMm * mm;
            halfTop = hw;
            depthLow = D;
            if (s.type == Type::WireW)
            {
               const V3 p0(-hw, 0, 0), p1(-0.5 * hw, D, 0), p2(0, D * (1.0 - s.wApex), 0), p3(0.5 * hw, D, 0), p4(hw, 0, 0);
               wire(p0, p1, true); wire(p1, p2, true);
               surfaces = 2;   // a ball in the V touches both sides
               surfaceNow = 1;
               wire(p2, p3, true); wire(p3, p4, true);
               surfaceNow = 0;
            }
            else
            {
               const double hb = hw * s.bottomFrac;
               const V3 p0(-hw, 0, 0), p1(-hb, D, 0), p2(hb, D, 0), p3(hw, 0, 0);
               wire(p0, p1, true); wire(p1, p2, true); wire(p2, p3, true);
            }
            // pins into the bracket holes, the stop crank and the switch crank (outside the bracket: mass only)
            wire(V3(-hw - 4 * mm, 0, 0), V3(-hw, 0, 0), false);
            wire(V3(hw, 0, 0), V3(hw + 4 * mm, 0, 0), false);
            wire(V3(-hw - 4 * mm, 0, 0), V3(-hw - 4 * mm, arm, 0), false);
            wire(V3(hw + 4 * mm, 0, 0), V3(hw + 4 * mm, arm, 0), false);
         }
         else
         {
            const bool lng = s.type == Type::LongPlate;
            const double top = (lng ? s.longTopMm : s.plateTopMm) * mm, bot = (lng ? s.longBottomMm : s.plateBottomMm) * mm, L = (lng ? s.longLenMm : s.plateLenMm) * mm;
            halfTop = 0.5 * top;
            plateTop = top;
            plateBot = bot;
            plateLen = L;
            AddPlate(V3(-0.5 * top, 0, 0), V3(0.5 * top, 0, 0), V3(0.5 * bot, L, 0), V3(-0.5 * bot, L, 0), sigma);
            depthLow = L;
            if (lng)
            {
               const double a = s.lipBendDeg * GroupStep::kPi / 180.0, l = s.lipMm * mm;
               const V3 dir(0, std::cos(a), -std::sin(a));
               AddPlate(V3(-0.5 * bot, L, 0), V3(0.5 * bot, L, 0), V3(0.5 * bot, L, 0) + dir * l, V3(-0.5 * bot, L, 0) + dir * l, sigma);
               depthLow = L + l * std::cos(a);
            }
            wire(V3(-0.5 * top - 4 * mm, 0, 0), V3(0.5 * top + 4 * mm, 0, 0), true);   // the pin wire across the top (the ball can reach it on a bounce) and its ends in the bracket
         }
         FinishShape(std::sqrt(halfTop * halfTop + depthLow * depthLow) + 8e-3);

         // stops: the wire crank (radius wire/2) or the plate edge (t/2) on the bracket edge
         const double rPart = s.type == Type::Plate || s.type == Type::LongPlate ? halfT : rho, deg = GroupStep::kPi / 180.0;
         const auto stop = [&](const double angDeg, const int side) {
            GroupStep::BodyStopPlace st;
            st.at = angDeg * deg;
            st.side = side;
            st.lever = arm;
            st.Setup(mp.I, std::sqrt(rPart * s.bracketEdgeMm * mm), s.stopE, 0.1);
            return st;
         };
         q = s.oneWay ? s.restAngleDeg * deg : 0.0;
         stopRest = stopBack = -1;
         if (s.oneWay)
         {
            stopRest = PlaceStop(stop(s.restAngleDeg, -1));
            stopOpen = PlaceStop(stop(s.openStopDeg, 1));
         }
         else
         {
            stopOpen = PlaceStop(stop(s.openStopDeg, 1));
            stopBack = PlaceStop(stop(-s.openStopDeg, -1));
         }
         // the hold actuator: a driven stop on the closing side, parked at the rest angle (one-way) or behind the open stop's mirror (two-way)
         GroupStep::BodyStopPlace hold = stop(s.oneWay ? s.restAngleDeg : -s.openStopDeg, -1);
         hold.target = hold.at;
         holdPark = hold.at;
         stopHold = PlaceStop(hold);

         // the pass angle: the lowest point clears a ball on the surface; the switch closes lead_deg before it
         const double cp = (s.pivotHeightMm * mm - 2.0 * b.R) / depthLow;
         passDeg = cp >= 1.0 ? 0.0 : cp <= -1.0 ? 180.0 : std::acos(cp) / deg;
         switchMake = (passDeg - s.switchLeadDeg) * deg;
         switchBreak = switchMake - s.switchHystDeg * deg;
         if (s.hasSwitch && s.switchForceN > 0.0)   // the blade is pressed over the last 5 deg before the make angle; switchForceN at the make angle, at the crank's lever
         {
            aBlade = switchMake - 5.0 * deg;
            kBlade = s.switchForceN * arm / (5.0 * deg);
         }
         SetContactLaw(s.elasticity, s.friction, b);
         CycleBegin();
      }
      // The script's Open / Close
      void Hold(const bool open) { DriveStop(stopHold, open ? set.holdAngleDeg * GroupStep::kPi / 180.0 : holdPark, set.holdSpeedDps * GroupStep::kPi / 180.0); }
      bool CanPass() const { return passDeg < set.openStopDeg; }
      double PassAngleDeg() const { return passDeg; }

      // A Hit only fires near the start angle, read at the step's start: a light gate turns several degrees in the step a touch starts in
      bool NearStart(const double th) const override
      {
         const double start = set.oneWay ? set.restAngleDeg * GroupStep::kPi / 180.0 : HangAngle();
         return std::fabs(th - start) * 180.0 / GroupStep::kPi <= set.hitWindowDeg;
      }
      int LimitKind(const int stop) const override { return stop == stopOpen || stop == stopBack ? 1 : stop == stopRest || (stop == stopHold && set.oneWay) ? -1 : 0; }   // (a one-way gate's parked hold stop is at the rest angle)
      void UpdateSwitches() override
      {
         if (!set.hasSwitch)
            return;
         const double th = set.oneWay ? q : std::fabs(q);
         if (!switchClosed && th >= switchMake)
         {
            switchClosed = true;
            swLog.push_back(1);
         }
         else if (switchClosed && th < switchBreak)
         {
            switchClosed = false;
            swLog.push_back(0);
         }
      }
      // bend stiffness out of the gate's plane at distance y below the axis: wires k = 2 x 3 E I / y^3 (two legs share, I = pi d^4 / 64); plates k = 3 E (w t^3 / 12) / y^3 (w = mean width between the axis and y)
      double BendKb(const double y) const override
      {
         const double Ey = set.steelEGpa * 1e9;
         if (set.type == Type::WireW || set.type == Type::WireRectangle)
         {
            const double d = set.wireDMm * 1e-3, I2 = GroupStep::kPi * d * d * d * d / 64.0;
            return 2.0 * 3.0 * Ey * I2 / (y * y * y);
         }
         const double t = set.plateThickMm * 1e-3, f = std::min(1.0, y / plateLen), wMean = plateTop + (plateBot - plateTop) * 0.5 * f;
         return 3.0 * Ey * (wMean * t * t * t / 12.0) / (y * y * y);
      }

   private:
      double plateTop = 0.0, plateBot = 0.0, plateLen = 1.0, holdPark = 0.0, passDeg = 0.0, switchMake = 0.0, switchBreak = 0.0;
      int stopHold = -1, stopRest = -1, stopOpen = -1, stopBack = -1;
   };
}
