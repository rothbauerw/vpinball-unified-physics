// license:GPLv3+
//
// A spinner of the unified contact model (docs/GateSpinnerDesign.txt, section 5): a flat steel blade on an offset shaft, a rigid body on one axis
// (contacthingebody.h) with the shared ball contact and bending (contacthingepart.h).
//
// Pure SI maths with no VPX types (offline tests: runs/2026-10-14-hinge).
//
// Shape (body coordinates: x along the axis E, y down at theta = 0, z toward the front N at theta = 0; metres):
//   blade   a rectangle width x height x thickness (default 39 x 31.75 x 0.8 mm), steel. The shaft runs across it with Z bends, so the axis lies IN THE PLANE of
//           the blade, offset above the blade's centre (default 6 mm): the blade reaches height/2 - offset above the axis and height/2 + offset below it.
//           Gravity on the offset centre of mass is the pendulum that brings the blade back toward vertical after every spin.
//   shaft   pins beyond the blade into the bracket (mass only). The bracket is not collidable (a table that wants it uses a primitive).
// Bend: the blade bends out of its plane like a cantilever from the shaft line, k = 3 E (width t^3 / 12) / y^3 at the distance y from the axis (above or below).
// Damping: the stock Damping (the fraction of the angular speed kept per 10 ms of free spin) as a viscous loss c = -I ln(Damping) / 0.01 s, plus the Coulomb
// friction of the pin wires in their holes (contacthingebody.h).
// Switch: none (Spin when the angle passes switch_angle), opto (a flag interrupts a beam while the angle is within half_width of switch_angle: no force), leaf
// (a crank on the shaft presses a long thin metal blade over the top half of the turn, switch_angle +- contact_half_deg; the contacts are closed within
// switch_angle +- closed_half_deg; force leaf_force_n at full bend, small loss) or micro (the same crank on a lever with a preload micro_force_n). Leaf and
// micro are one spring-switch model (user 2026-10-14 / 15): they take energy while pressed and give it back as they are released,
// so the spinner speeds up and slows down over its switch; near rest the switch torque can hold the blade a little off vertical.
// Stops: only when angle_min / angle_max limit the turn (a spinner used as a flap); a normal spinner turns freely.
// Events: Hit when a ball strikes the blade near rest, Spin and SwitchClosed when the switch closes, SwitchOpened when it opens, LimitEOS / LimitBOS at stops.

#pragma once

#include <cmath>
#include <algorithm>
#include <vector>

#include "physics/rubbercontact.h"
#include "physics/contacthingebody.h"
#include "physics/contacthingepart.h"

namespace SpinnerModel
{
   using RubberContact::V3;
   using RubberContact::BallSI;
   using HingePart::Event;
   using HingePart::BallLink;

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
      double damping = 0.9879;        // stock Damping: angular speed kept per 10 ms of free spin
      Switch sw = Switch::Leaf;
      double switchAngleDeg = 180.0;  // where the crank / the flag acts (body angle; 180 = blade upside down, so a blade at rest is free)
      double optoHalfWidthDeg = 10.0;
      // Leaf and micro switch (user 2026-10-15): the crank presses the long thin switch blade over the top half of a turn (switch_angle +- contact_half_deg, 90:
      // the blade rests on the crank's circle at 90 / 270 deg from rest, so the bottom half never touches it); the contacts are closed within switch_angle +-
      // closed_half_deg (30: 150 - 210 deg from rest), and part again hyst_deg further out (the contact gap's hysteresis, no chatter).
      double camRMm = 3.0;            // crank radius
      double contactHalfDeg = 90.0, closedHalfDeg = 30.0, hystDeg = 1.0;
      double leafForceN = 0.02;       // at full bend (the crank pointing at the blade); default to be set by feel (user)
      double microForceN = 0.02;      // micro switch preload; to be set by feel
      double microKNPerM = 5.0;
      double switchLossNsPerM = 0.01; // the blade's loss
      double angleMinDeg = 0.0, angleMaxDeg = 0.0;   // equal: free turning
      double stopArmMm = 5.0, stopE = 0.5, bracketEdgeMm = 0.4;
      double steelEGpa = 200.0;
      double releaseGapMm = 0.2;
      double hitWindowDeg = 5.0;
      double elasticity = 0.3, friction = 0.02;
   };

   class Spinner : public HingePart::Part
   {
   public:
      Settings set;

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
         m_rho = 0.5 * s.wireDMm * mm;
         m_halfT = 0.5 * s.thickMm * mm;
         m_releaseGap = s.releaseGapMm * mm;
         const double hw = 0.5 * s.widthMm * mm, H = s.heightMm * mm, e = s.offsetMm * mm;
         const double top = -(0.5 * H - e), bot = 0.5 * H + e;
         AddPlate(V3(-hw, top, 0), V3(hw, top, 0), V3(hw, bot, 0), V3(-hw, bot, 0), s.density * s.thickMm * mm);
         AddWire(V3(-hw - 5 * mm, 0, 0), V3(-hw, 0, 0), lam, false);   // pins into the bracket
         AddWire(V3(hw, 0, 0), V3(hw + 5 * mm, 0, 0), lam, false);
         m_reach = std::sqrt(hw * hw + bot * bot) + 8e-3;
         m_width = 2.0 * hw;
         // viscous loss from the stock Damping (exact for a free exponential spin-down)
         body.c = (s.damping > 0.0 && s.damping < 1.0) ? -body.mp.I * std::log(s.damping) / 0.01 : 0.0;
         // optional stops (a spinner used as a flap)
         if (s.angleMaxDeg > s.angleMinDeg)
         {
            auto stop = [&](double angDeg, int side) {
               HingeBody::Stop st;
               st.angle = angDeg * M_PI / 180.0;
               st.side = side;
               st.lever = s.stopArmMm * mm;
               st.reff = std::sqrt(m_halfT * s.bracketEdgeMm * mm);
               HingeBody::SetupStop(st, body.mp.I, s.stopE, 0.1);
               return st;
            };
            m_stopMax = body.PlaceStop(stop(s.angleMaxDeg, 1));
            m_stopMin = body.PlaceStop(stop(s.angleMinDeg, -1));
         }
         // the switch
         m_switchAngle = s.switchAngleDeg * M_PI / 180.0;
         if (s.sw == Switch::Leaf || s.sw == Switch::Micro)
         {
            HingeBody::Cam cm;
            cm.phase = m_switchAngle;
            const double deg = M_PI / 180.0;
            cm.r = s.camRMm * mm;
            cm.gap = cm.r * std::cos(s.contactHalfDeg * deg);   // the blade's distance from the shaft: the crank reaches it at +- contact_half_deg
            const double travel = std::max(cm.r - cm.gap, 1e-5);
            m_closeD = cm.r * std::cos(s.closedHalfDeg * deg) - cm.gap;                  // the bend at which the contacts touch
            m_openD = cm.r * std::cos((s.closedHalfDeg + s.hystDeg) * deg) - cm.gap;     // and part
            cm.k = s.sw == Switch::Leaf ? s.leafForceN / travel : s.microKNPerM;
            cm.F0 = s.sw == Switch::Micro ? s.microForceN : 0.0;
            cm.c = s.switchLossNsPerM;
            m_cam = (int)body.cams.size();
            body.cams.push_back(cm);
         }
         m_switchClosed = false;
         m_prevWrapped = Wrap(body.theta - m_switchAngle);
         SetupContactLaw(s.elasticity, s.friction);
         UpdateWorld();
      }

      double BendKb(double y) const override
      {
         const double E = set.steelEGpa * 1e9, t = set.thickMm * 1e-3;
         return 3.0 * E * (m_width * t * t * t / 12.0) / (y * y * y);
      }
      // a Hit only fires on a blade near rest (hanging within hit_window of vertical): a ball touching a spinning blade is not a new hit
      bool NearStart() const override { return std::abs(Wrap(body.theta)) <= set.hitWindowDeg * M_PI / 180.0; }

      // the switch and the stop events after a substep
      void UpdateSwitches() override
      {
         const double wr = Wrap(body.theta - m_switchAngle);
         switch (set.sw)
         {
         case Switch::None:
            // no switch modelled: Spin as the angle passes the switch angle (either way)
            if ((m_prevWrapped < 0.0) != (wr < 0.0) && std::abs(wr - m_prevWrapped) < M_PI)
               events.push_back({ Event::Spin, -1, 0.0 });
            break;
         case Switch::Opto:
         {
            const bool in = std::abs(wr) <= set.optoHalfWidthDeg * M_PI / 180.0;
            SetSwitch(in);
            break;
         }
         case Switch::Leaf:
         case Switch::Micro:
         {
            const double d = body.cams[m_cam].Delta(body.theta);
            if (!m_switchClosed && d >= m_closeD)
               SetSwitch(true);
            else if (m_switchClosed && d < m_openD)
               SetSwitch(false);
            break;
         }
         }
         m_prevWrapped = wr;
         for (const auto &e : body.events)
         {
            if (e.stop == m_stopMax)
               events.push_back({ Event::LimitEOS, -1, e.speedDegS });
            else if (e.stop == m_stopMin)
               events.push_back({ Event::LimitBOS, -1, e.speedDegS });
         }
         body.events.clear();
      }

      void SaveExtra(State &o) const override
      {
         Put(o, m_switchClosed);
         Put(o, m_prevWrapped);
      }
      void LoadExtra(StateReader &r) override
      {
         r.Get(m_switchClosed);
         r.Get(m_prevWrapped);
      }

      bool SwitchClosed() const { return m_switchClosed; }
      static double Wrap(double a)
      {
         a = std::fmod(a + M_PI, 2.0 * M_PI);
         if (a < 0.0)
            a += 2.0 * M_PI;
         return a - M_PI;
      }

   private:
      double m_width = 0.0, m_switchAngle = 0.0, m_prevWrapped = 0.0, m_closeD = 0.0, m_openD = 0.0;
      int m_cam = -1, m_stopMax = -1, m_stopMin = -1;
      bool m_switchClosed = false;
      void SetSwitch(bool closed)
      {
         if (closed == m_switchClosed)
            return;
         m_switchClosed = closed;
         if (closed)
         {
            events.push_back({ Event::SwitchClosed, -1, 0.0 });
            events.push_back({ Event::Spin, -1, 0.0 });
         }
         else
            events.push_back({ Event::SwitchOpened, -1, 0.0 });
      }
   };
}
