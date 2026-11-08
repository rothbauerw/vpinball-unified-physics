// license:GPLv3+
//
// A spinner as a part of the group step (docs/unified-physics/GateSpinnerDesign.txt section 5; runs/2026-10-31-gate-spinner-disc): a flat steel blade on an offset shaft, on one axis (hingegroup.h, bodypart.h).
//   blade   a rectangle width x height x thickness (default 39 x 31.75 x 0.8 mm), steel. The shaft runs across it with Z bends, so the axis lies IN THE PLANE of the blade, offset above the blade's centre (default 6 mm): the blade reaches height/2 - offset above the axis and height/2 + offset below it.
//           Gravity on the offset centre of mass is the pendulum that brings the blade back toward vertical after every spin.
//   shaft   pins beyond the blade into the bracket (mass only).
// Bend: the blade bends out of its plane like a cantilever from the shaft line, k = 3 E (width t^3 / 12) / y^3 at the distance y from the axis.
// Damping: the stock Damping (the fraction of the angular speed kept per 10 ms of free spin) as a viscous loss c = -I ln(Damping) / 0.01 s, plus the Coulomb friction of the pin wires in their holes.
// Switch: none (Spin when the angle passes switch_angle), opto (a flag interrupts a beam while the angle is within half_width of switch_angle: no force), leaf (a crank on the shaft presses a long thin metal blade over the top half of the turn, switch_angle +- contact_half_deg; the contacts are
// closed within switch_angle +- closed_half_deg; force leaf_force_n at full bend, small loss) or micro (the same crank on a lever with a preload micro_force_n). They take energy while pressed and give it back as they are released, so the spinner speeds up and slows down over its switch.
// Stops: only when angle_min / angle_max limit the turn (a spinner used as a flap); a normal spinner turns freely. Pure SI.

#pragma once

#include <algorithm>
#include <cmath>

#include "physics/hingegroup.h"

namespace SpinnerGroup
{
   using GroupStep::V3;
   using RubberContact::BallSI;

   enum class Switch { None, Opto, Leaf, Micro };

   struct Settings
   {
      double widthMm = 39.0;          // the part's length
      double heightMm = 31.75;        // 1.25 in
      double thickMm = 0.8;
      double offsetMm = 6.0;          // axis above the blade's centre, in the blade's plane
      double pivotHeightMm = 31.75;
      double wireDMm = 1.2;
      double density = 7850.0;
      double muPivot = 0.1;           // pin wire in its hole; 0.1 by feel in play (user 2026-10-15: 0.3 spun 1/2 - 3/4 as much as a real spinner)
      double muPivotSpring = 0.3;     // the pin wire in its hole, steel on steel (mu_pivot sets it too); no stock damping, the blade's air drag instead
      Switch sw = Switch::Leaf;
      double switchAngleDeg = 180.0;  // where the crank / the flag acts (body angle; 180 = blade upside down, so a blade at rest is free)
      double optoHalfWidthDeg = 10.0;
      // Leaf and micro switch (user 2026-10-15): the crank presses the long thin switch blade over the top half of a turn (switch_angle +- contact_half_deg, 90: the blade rests on the crank's circle at 90 / 270 deg from rest, so the bottom half never touches it);
      // the contacts are closed within switch_angle +- closed_half_deg (30: 150 - 210 deg from rest), and part again hyst_deg further out (the contact gap's hysteresis, no chatter).
      double camRMm = 3.0;            // crank radius
      double contactHalfDeg = 90.0, closedHalfDeg = 30.0, hystDeg = 1.0;
      double leafForceN = 0.02;       // at full bend (the crank pointing at the blade); default to be set by feel (user)
      double microForceN = 0.02;      // micro switch preload; to be set by feel
      double microKNPerM = 5.0;
      double switchLossNsPerM = 0.01; // the blade's loss
      double angleMinDeg = 0.0, angleMaxDeg = 0.0;   // equal: free turning
      double stopArmMm = 5.0, stopE = 0.5, bracketEdgeMm = 0.4;
      double steelEGpa = 200.0;
      double hitWindowDeg = 5.0;
      double elasticity = 0.3, friction = 0.02;
   };

   class Part : public HingeGroup::Part
   {
   public:
      Settings set;

      void Build(const Settings &s, const V3 &surfacePoint, const V3 &front, const BallSI &b, const V3 &gravity)
      {
         set = s;
         Begin(surfacePoint + V3(0.0, 0.0, s.pivotHeightMm * 1e-3), front, gravity);
         muPivot = s.muPivot;
         rPin = 0.5 * s.wireDMm * 1e-3;
         const double mm = 1e-3, lam = s.density * GroupStep::kPi * 0.25 * s.wireDMm * s.wireDMm * mm * mm, deg = GroupStep::kPi / 180.0;
         rho = 0.5 * s.wireDMm * mm;
         halfT = 0.5 * s.thickMm * mm;
         const double hw = 0.5 * s.widthMm * mm, H = s.heightMm * mm, e = s.offsetMm * mm, top = -(0.5 * H - e), bot = 0.5 * H + e;
         AddPlate(V3(-hw, top, 0), V3(hw, top, 0), V3(hw, bot, 0), V3(-hw, bot, 0), s.density * s.thickMm * mm);
         AddWire(V3(-hw - 5 * mm, 0, 0), V3(-hw, 0, 0), lam, false);   // pins into the bracket
         AddWire(V3(hw, 0, 0), V3(hw + 5 * mm, 0, 0), lam, false);
         FinishShape(std::sqrt(hw * hw + bot * bot) + 8e-3);
         width = 2.0 * hw;
         blowest = s.pivotHeightMm * mm - bot;
         // optional stops (a spinner used as a flap)
         stopMax = stopMin = -1;
         if (s.angleMaxDeg > s.angleMinDeg)
         {
            const auto stop = [&](const double angDeg, const int side) {
               GroupStep::BodyStopPlace st;
               st.at = angDeg * deg;
               st.side = side;
               st.lever = s.stopArmMm * mm;
               st.Setup(mp.I, std::sqrt(halfT * s.bracketEdgeMm * mm), s.stopE, 0.1);
               return st;
            };
            stopMax = PlaceStop(stop(s.angleMaxDeg, 1));
            stopMin = PlaceStop(stop(s.angleMinDeg, -1));
         }
         // the switch
         switchAngle = s.switchAngleDeg * deg;
         if (s.sw == Switch::Leaf || s.sw == Switch::Micro)
         {
            HingeBody::Cam cm;
            cm.phase = switchAngle;
            cm.r = s.camRMm * mm;
            cm.gap = cm.r * std::cos(s.contactHalfDeg * deg);   // the blade's distance from the shaft: the crank reaches it at +- contact_half_deg
            const double travel = std::max(cm.r - cm.gap, 1e-5);
            closeD = cm.r * std::cos(s.closedHalfDeg * deg) - cm.gap;                  // the bend at which the contacts touch
            openD = cm.r * std::cos((s.closedHalfDeg + s.hystDeg) * deg) - cm.gap;     // and part
            cm.k = s.sw == Switch::Leaf ? s.leafForceN / travel : s.microKNPerM;
            cm.F0 = s.sw == Switch::Micro ? s.microForceN : 0.0;
            cm.c = s.switchLossNsPerM;
            cams.push_back(cm);
         }
         bladeA = 0.5 * H - e;
         bladeB = bot;
         ApplyModel();
         prevWrapped = Wrap(q - switchAngle);
         SetContactLaw(s.elasticity, s.friction, b);
         CycleBegin();
      }
      double BladeLowestM() const { return blowest; }   // above the surface, hanging straight
      // the leaf (micro) switch rests on the shaft's actuator at every angle and rubs, its friction in the pin's limit (the same stick spring)
      double FrictionLimit(const double th, const double om, const V3 &ballForce) const override
      {
         const double L = HingeGroup::Part::FrictionLimit(th, om, ballForce);
         return L + rubNm;
      }
      // the blade's air drag, a flat plate turning about an axis in its plane: 0.5 rho Cd W w|w| (a^4 + b^4) / 4 at the step's mid rate (a loss, >= 0)
      double StoredChord(const double th, const double dq, const double h, double &slope, double &lossF, double &lossSlope) const override
      {
         const double chord = HingeGroup::Part::StoredChord(th, dq, h, slope, lossF, lossSlope);
         if (dragC > 0.0)
         {
            const double r = dq / h;
            lossF += dragC * r * std::fabs(r);
            lossSlope += 2.0 * dragC * std::fabs(r) / h;
         }
         return chord;
      }
      // bend stiffness of the blade at the distance y from the axis
      double BendKb(const double y) const override
      {
         const double Ey = set.steelEGpa * 1e9, t = set.thickMm * 1e-3;
         return 3.0 * Ey * (width * t * t * t / 12.0) / (y * y * y);
      }
      // a Hit only fires on a blade near rest (hanging within hit_window of vertical): a ball touching a spinning blade is not a new hit
      bool NearStart(const double th) const override { return std::fabs(Wrap(th - HangAngle())) <= set.hitWindowDeg * GroupStep::kPi / 180.0; }
      int LimitKind(const int stop) const override { return stop == stopMax ? 1 : stop == stopMin ? -1 : 0; }
      void UpdateSwitches() override
      {
         const double wr = Wrap(q - switchAngle);
         switch (set.sw)
         {
         case Switch::None:   // no switch modelled: Spin as the angle passes the switch angle (either way)
            if ((prevWrapped < 0.0) != (wr < 0.0) && std::fabs(wr - prevWrapped) < GroupStep::kPi)
               swLog.push_back(2);
            break;
         case Switch::Opto:
            SetSwitch(std::fabs(wr) <= set.optoHalfWidthDeg * GroupStep::kPi / 180.0);
            break;
         case Switch::Leaf:
         case Switch::Micro:
         {
            const double d = cams[0].Delta(q);
            if (!switchClosed && d >= closeD)
               SetSwitch(true);
            else if (switchClosed && d < openD)
               SetSwitch(false);
            break;
         }
         }
         prevWrapped = wr;
      }
      static double Wrap(double a)
      {
         a = std::fmod(a + GroupStep::kPi, 2.0 * GroupStep::kPi);
         if (a < 0.0)
            a += 2.0 * GroupStep::kPi;
         return a - GroupStep::kPi;
      }

   private:
      // the losses of the pin's law: the pin at mu_pivot_spring, the switch's rub and the blade's air drag (air 1.2 kg/m^3, flat plate Cd 1.2)
      void ApplyModel()
      {
         cViscous = 0.0;
         muPivot = set.muPivotSpring;
         for (HingeBody::Cam &cm : cams)
            cm.c = 0.0;
         const double leafF = set.sw == Switch::Leaf ? set.leafForceN : set.sw == Switch::Micro ? set.microForceN : 0.0;
         rubNm = leafF * muPivot * set.camRMm * 1e-3;
         dragC = 0.5 * 1.2 * 1.2 * width * (std::pow(bladeA, 4) + std::pow(bladeB, 4)) / 4.0;
      }
      void SetSwitch(const bool closed)
      {
         if (closed == switchClosed)
            return;
         switchClosed = closed;
         swLog.push_back(closed ? 1 : 0);
         if (closed)
            swLog.push_back(2);
      }
      double width = 0.0, blowest = 0.0, switchAngle = 0.0, closeD = 0.0, openD = 0.0, bladeA = 0.0, bladeB = 0.0, rubNm = 0.0, dragC = 0.0;
      int stopMax = -1, stopMin = -1;
   };
}
