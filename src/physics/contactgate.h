// license:GPLv3+
//
// A gate of the unified contact model (docs/GateSpinnerDesign.txt, section 4): the four real gate types as rigid bodies on one axis (contacthingebody.h).
//
// Pure SI maths with no VPX types (offline tests: runs/2026-10-14-hinge). The owner places the gate (axis point, front direction, the surface height), feeds the
// gravity it feels, steps it in the island of the balls touching it (or alone), and reads the angle and the events.
//
// Shape (body coordinates: x along the axis E, y down at theta = 0, z toward the front N at theta = 0; metres):
//   wire W          pins at both ends; legs from the pin ends down to two bottoms at +-width/4, a middle V that rises to apex * depth
//   wire rectangle  legs from the pin ends down to a bottom bar (bottom width = width * bottom_frac; < 1 is the deep trapezoid variant)
//   plate           trapezoid, top width 27 mm, bottom 18 mm, 20 mm long, hanging from a pin wire across its top
//   long plate      trapezoid 27 / 13.5 mm, 27 mm long, plus a lip at the bottom bent back by lip_bend so it stands vertical at the one-way rest angle
// A wire is a chain of straight capsules (tube radius wire_d / 2); a plate is a sheet: its mid-surface triangles with half the thickness as the probe offset.
// Wire depth: the lowest wire meets a rolling ball hit_offset above its centre: depth = pivot height - R - hit_offset (or depth_mm).
//
// Stops (steel on steel, contacthingebody.h): one-way: a reverse stop at the rest angle (the gate leans on it) and a forward stop at the open angle;
// two-way: stops at +- the open angle, the gate hangs at 0. Hold: a driven stop that pushes the gate to the hold angle (Open) and goes back (Close).
// Switch: closes when the gate is open far enough for the ball to pass, lead_deg before that; opens again hyst_deg lower. Optional blade spring.
// Ball contact, bending of the wire and plates, the tilt force and the step rule: contacthingepart.h (shared with the spinner).
// Bend stiffness out of the gate plane at depth y below the axis: wires k = 2 x 3 E I / y^3 (two legs share, I = pi d^4 / 64); plates k = 3 E (w t^3 / 12) / y^3
// (w = mean width between the axis and y).
// Events: Hit when a ball starts touching the gate (after it was more than release_gap away) while the gate is within hit_window of its start angle,
// SwitchClosed / SwitchOpened, LimitEOS / LimitBOS from the stops.

#pragma once

#include <cmath>
#include <algorithm>
#include <vector>

#include "physics/rubbercontact.h"
#include "physics/contacthingebody.h"
#include "physics/contacthingepart.h"

namespace GateModel
{
   using RubberContact::V3;
   using RubberContact::BallSI;
   using RubberContact::TubeProbe;
   using RubberContact::RubberParams;
   using RubberContact::RubberState;
   using RubberContact::ContactResult;
   using RubberContact::SurfaceTri;
   using HingePart::Event;
   using HingePart::BallLink;

   enum class Type { WireW, WireRectangle, Plate, LongPlate };
   enum class Switch { None, Leaf };

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
      double plateTopMm = 27.0, plateBottomMm = 18.0, plateLenMm = 20.0;            // plate
      double longTopMm = 27.0, longBottomMm = 13.5, longLenMm = 27.0;               // long plate
      double plateThickMm = 0.65;
      double lipMm = 2.75, lipBendDeg = 15.0;
      bool oneWay = true;
      double restAngleDeg = 15.0;     // one-way rest (on the reverse stop)
      double openStopDeg = 90.0;
      double stopArmMm = 5.0, stopE = 0.5, bracketEdgeMm = 0.4;   // stop lever, restitution at 0.1 m/s, radius of the bracket edge it touches
      double muPivot = 0.3;
      double density = 7850.0;
      Switch sw = Switch::Leaf;
      double switchLeadDeg = 5.0, switchHystDeg = 3.0, switchForceN = 0.0;
      double holdAngleDeg = 70.0, holdSpeedDps = 900.0;
      double steelEGpa = 200.0;       // wire and plate Young's modulus for the bending
      double releaseGapMm = 0.2;      // a new Hit needs the ball to have been this far from the gate
      double hitWindowDeg = 5.0;      // Hit only fires while the gate is within this of its start angle (the strike that opens it, user 2026-10-15)
      // ball contact material: steel on steel, restitution and friction from the part
      double elasticity = 0.3, friction = 0.02;
   };

   class Gate : public HingePart::Part
   {
   public:
      Settings set;

      // ---- build ----
      // axisPoint: the axis centre (world, m; z = surface + pivot height); front: the direction a ball passes from the front (horizontal)
      void Build(const Settings &s, const V3 &surfacePoint, const V3 &front, const BallSI &b)
      {
         set = s;
         ball = b;
         body = HingeBody::Body();
         body.SetFrame(surfacePoint + V3(0.0, 0.0, s.pivotHeightMm * 1e-3), front);
         body.muPivot = s.muPivot;
         body.rPin = 0.5 * s.wireDMm * 1e-3;
         ClearShape();
         const double mm = 1e-3, lam = s.density * M_PI * 0.25 * s.wireDMm * s.wireDMm * mm * mm;
         const double R = b.R;
         m_rho = 0.5 * s.wireDMm * mm;
         m_halfT = 0.5 * s.plateThickMm * mm;
         m_releaseGap = s.releaseGapMm * mm;
         auto wire = [&](const V3 &p, const V3 &q, bool touchable) { AddWire(p, q, lam, touchable); };
         auto plate = [&](const V3 &a, const V3 &bb, const V3 &c, const V3 &d) { AddPlate(a, bb, c, d, s.density * s.plateThickMm * mm); };
         const double arm = s.stopArmMm * mm;
         double halfTop, depthLow;   // half width at the axis, depth of the lowest point below the axis
         if (s.type == Type::WireW || s.type == Type::WireRectangle)
         {
            const double W = s.widthMm * mm, hw = 0.5 * W;
            const double D = s.depthMm > 0.0 ? s.depthMm * mm : s.pivotHeightMm * mm - R - s.hitOffsetMm * mm;
            halfTop = hw;
            depthLow = D;
            if (s.type == Type::WireW)
            {
               const V3 p0(-hw, 0, 0), p1(-0.5 * hw, D, 0), p2(0, D * (1.0 - s.wApex), 0), p3(0.5 * hw, D, 0), p4(hw, 0, 0);
               wire(p0, p1, true); wire(p1, p2, true); wire(p2, p3, true); wire(p3, p4, true);
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
            m_plateTop = top;
            m_plateBot = bot;
            m_plateLen = L;
            plate(V3(-0.5 * top, 0, 0), V3(0.5 * top, 0, 0), V3(0.5 * bot, L, 0), V3(-0.5 * bot, L, 0));
            depthLow = L;
            if (lng)
            {
               // the lip continues the plate past its bottom edge, bent back (toward -front) by lipBend: vertical when the plate hangs at that angle
               const double a = s.lipBendDeg * M_PI / 180.0, l = s.lipMm * mm;
               const V3 dir(0, std::cos(a), -std::sin(a));
               plate(V3(-0.5 * bot, L, 0), V3(0.5 * bot, L, 0), V3(0.5 * bot, L, 0) + dir * l, V3(-0.5 * bot, L, 0) + dir * l);
               depthLow = L + l * std::cos(a);
            }
            // the pin wire across the top (the ball can reach it on a bounce) and its ends in the bracket
            wire(V3(-0.5 * top - 4 * mm, 0, 0), V3(0.5 * top + 4 * mm, 0, 0), true);
         }
         m_depthLow = depthLow;
         m_halfTop = halfTop;
         m_reach = std::sqrt(halfTop * halfTop + depthLow * depthLow) + 8e-3;

         // stops: the touching part is the wire crank (radius wire/2) or the plate edge (t/2) on the bracket edge
         const double rPart = (s.type == Type::Plate || s.type == Type::LongPlate) ? m_halfT : m_rho;
         auto stop = [&](double angDeg, int side, double e) {
            HingeBody::Stop st;
            st.angle = angDeg * M_PI / 180.0;
            st.side = side;
            st.lever = arm;
            st.reff = std::sqrt(rPart * s.bracketEdgeMm * mm);
            HingeBody::SetupStop(st, body.mp.I, e, 0.1);
            return st;
         };
         body.theta = s.oneWay ? s.restAngleDeg * M_PI / 180.0 : 0.0;
         m_stopRest = m_stopOpen = m_stopBack = m_stopHold = -1;
         if (s.oneWay)
         {
            m_stopRest = body.PlaceStop(stop(s.restAngleDeg, -1, s.stopE));
            m_stopOpen = body.PlaceStop(stop(s.openStopDeg, 1, s.stopE));
         }
         else
         {
            m_stopOpen = body.PlaceStop(stop(s.openStopDeg, 1, s.stopE));
            m_stopBack = body.PlaceStop(stop(-s.openStopDeg, -1, s.stopE));
         }
         // the hold actuator: a driven stop on the closing side, parked at the rest angle (one-way) or behind the open stop's mirror (two-way)
         HingeBody::Stop hold = stop(s.oneWay ? s.restAngleDeg : -s.openStopDeg, -1, s.stopE);
         hold.target = hold.angle;
         m_holdPark = hold.angle;
         m_stopHold = body.PlaceStop(hold);

         // the pass angle: the lowest point clears a ball on the surface (cos = (pivot - 2R) / depth), the switch closes lead_deg before it
         const double cp = (s.pivotHeightMm * mm - 2.0 * R) / depthLow;
         m_passDeg = cp >= 1.0 ? 0.0 : (cp <= -1.0 ? 180.0 : std::acos(cp) * 180.0 / M_PI);
         m_canPass = m_passDeg < s.openStopDeg;
         m_switchMake = (m_passDeg - s.switchLeadDeg) * M_PI / 180.0;
         m_switchBreak = m_switchMake - s.switchHystDeg * M_PI / 180.0;
         m_switchClosed = false;
         if (s.sw == Switch::Leaf && s.switchForceN > 0.0)
         {
            // the blade is pressed over the last 5 deg before the make angle; switchForceN at the make angle, at the crank's lever
            HingeBody::Spring sp;
            sp.angle = m_switchMake - 5.0 * M_PI / 180.0;
            sp.side = 1;
            sp.k = s.switchForceN * arm / (5.0 * M_PI / 180.0);
            body.springs.push_back(sp);
         }

         SetupContactLaw(s.elasticity, s.friction);
         UpdateWorld();
      }

      // bend stiffness out of the gate's plane at distance y below the axis
      double BendKb(double y) const override
      {
         const double E = set.steelEGpa * 1e9;
         if (set.type == Type::WireW || set.type == Type::WireRectangle)
         {
            const double d = set.wireDMm * 1e-3, I = M_PI * d * d * d * d / 64.0;
            return 2.0 * 3.0 * E * I / (y * y * y);   // the two legs share the load
         }
         const double t = set.plateThickMm * 1e-3;
         const double f = std::min(1.0, y / m_plateLen);
         const double wMean = m_plateTop + (m_plateBot - m_plateTop) * 0.5 * f;   // mean width between the axis and the contact
         return 3.0 * E * (wMean * t * t * t / 12.0) / (y * y * y);
      }
      // a Hit only fires on a gate within hit_window of its start angle (the strike that opens it; user 2026-10-15)
      bool NearStart() const override
      {
         const double start = set.oneWay ? set.restAngleDeg : 0.0;
         return std::abs(AngleDeg() - start) <= set.hitWindowDeg;
      }

      bool CanPass() const { return m_canPass; }
      double PassAngleDeg() const { return m_passDeg; }

      // the switch and the stop events after a substep
      void UpdateSwitches() override
      {
         const double th = set.oneWay ? body.theta : std::abs(body.theta);
         if (set.sw != Switch::None)
         {
            if (!m_switchClosed && th >= m_switchMake)
            {
               m_switchClosed = true;
               events.push_back({ Event::SwitchClosed, -1, 0.0 });
            }
            else if (m_switchClosed && th < m_switchBreak)
            {
               m_switchClosed = false;
               events.push_back({ Event::SwitchOpened, -1, 0.0 });
            }
         }
         for (const auto &e : body.events)
         {
            if (e.stop == m_stopOpen || (e.stop == m_stopBack && !set.oneWay))
               events.push_back({ Event::LimitEOS, -1, e.speedDegS });
            else if (e.stop == m_stopRest)
               events.push_back({ Event::LimitBOS, -1, e.speedDegS });
         }
         body.events.clear();
      }

      // hold actuator (the script's Open / Close)
      void Hold(bool open)
      {
         if (m_stopHold < 0)
            return;
         body.DriveStop(m_stopHold, open ? set.holdAngleDeg * M_PI / 180.0 : m_holdPark, set.holdSpeedDps * M_PI / 180.0);
      }

      bool SwitchClosed() const { return m_switchClosed; }

   private:
      double m_plateTop = 0.0, m_plateBot = 0.0, m_plateLen = 1.0;
      double m_depthLow = 0.0, m_halfTop = 0.0;
      int m_stopRest = -1, m_stopOpen = -1, m_stopBack = -1, m_stopHold = -1;
      double m_holdPark = 0.0;
      double m_passDeg = 0.0, m_switchMake = 0.0, m_switchBreak = 0.0;
      bool m_canPass = true, m_switchClosed = false;
   };
}
