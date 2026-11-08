// license:GPLv3+
//
// A flipper bat as a moving part of the contact system, on the group step (docs/unified-physics/GroupStepDesign.txt; runs/2026-10-31-flipper-engine/notes.md). The bat, its coil (plunger, link, power and hold
// windings), return spring and stops are BatGroup::Part (batgroup.h): the engine moves them with the balls in one solve, so an island with a bat always goes to the group. There is no other path: no
// old-loop contact, no impulse, no stock.
//
// What this class adds to the group part: what the engine needs of any moving part (the probe for engagement and the overlap replay, the speed bound, the saved state), the bat's own steps when no ball is
// near, and the pivot's nudge.
// Laws:
//  - ball contact: the shared rubber law (GroupStep::Item) on the bat outline at its present angle; its reaction on the bat is in the same solve.
//  - closing speed: the tip's present speed plus what the coil and spring can add in a cycle, plus the pivot's speed (SurfaceSpeedBound).
//  - the pivot's nudge: a given move of a coordinate with no unknown, a drive: the group carries it inside the step and counts its work on the balls.
// Energy: only the coil and the pivot's nudge add energy; the contact, the link's free play and the stops remove it.
// Pure SI maths with no VPX types: positions are world metres; the pivot is at z = 0 (the rubber heights are from the playfield surface). Angles: the table's frame (start angle, 0 along -y), swung by dirSign.

#pragma once

#include <algorithm>
#include <cmath>

#include "physics/contactmovingpart.h"
#include "physics/contactbat.h"
#include "physics/batgroup.h"

namespace BatPart
{
   using V3 = RubberContact::V3;

   class Part : public ContactMovingPart
   {
   public:
      ContactBat bat;                    // the bat's settings file (contactbat.h)
      BatGroup::Part gp;                 // the bat, its coil and its stops, as a part of the group
      V3 pivot;                          // pivot at rest, m (z = 0)
      RubberContact::BallSI ball;
      double alphaMax = 0.0;             // the largest angular acceleration the coil and the spring can give the bat, rad/s^2 (Configure)
      double leverMax = 0.0, springMax = 0.0;   // the crank's longest lever (m) and the return spring's largest torque (N m) over the stroke (Configure)
      long failedSteps = 0, aloneSteps = 0, energySteps = 0;     // a step of the bat alone could not be solved (the detector: rule 8); steps taken alone (cost); of them taken on energy alone (rule 8)
      int bookMs = 0;                    // the table's time, for the book of a cycle alone

      // After a settings change: the group part, the rubber, the release gap and the acceleration bound. The state is kept.
      void Configure()
      {
         gp.Configure(bat.set, FlipperModel::MakeRubberParams(bat.set));
         const FlipperSettings &s = bat.set;
         const FlipperDerived &p = bat.der;
         leverMax = 0.0;
         springMax = 0.0;
         for (int i = 0; i <= 100; i++)
         {
            const double theta = p.theta0 + s.stopDeg * FlipperModel::DEG * i / 100.0;
            const FlipperModel::Kin k = FlipperModel::Kinematics(s, p, theta);
            leverMax = std::max(leverMax, std::fabs(k.lever) / 1000.0);
            springMax = std::max(springMax, std::fabs(FlipperModel::SpringTorqueK(s, p, k)));
         }
         alphaMax = ((p.strengthP + p.strengthH) * leverMax + springMax) / std::max(p.iFlipper, 1e-9);
         gp.geo.pivotRest = pivot;
         gp.MovePivot();
      }

      // the bat's angle in the outline's frame (rad from +x, counter-clockwise) and its angular velocity about +z
      double PhiMath() const { return gp.geo.angleStart + gp.geo.dirSign * gp.phi - M_PI / 2.0; }
      double OmegaZ() const { return gp.geo.dirSign * gp.omega; }
      double TipReach() const { return gp.geo.ol.Lc + gp.geo.ol.re; }
      void SetCoil(const bool on)
      {
         gp.SetCoil(on);
         idle = false;
      }
      void NudgeStart(const double dispM, const double riseS, const double fallS)
      {
         gp.pulse.Start(dispM, riseS, fallS);   // (one already running is left to finish)
         idle = false;
      }

      // ---- the moving-part interface (contactmovingpart.h) ----
      void SetBall(const RubberContact::BallSI &b) override { ball = b; }
      void SetGravity(const V3 &) override {}   // the bat turns in the playfield's plane about an axis normal to it: gravity does not turn it
      V3 Origin() const override { return pivot; }
      double Reach() const override { return TipReach(); }

      RubberContact::TubeProbe Probe(const V3 &c) const override
      {
         const V3 &pv = gp.geo.pivot;
         const V3 cRel = c - pv;
         const RubberContact::Probe p = RubberContact::ProbeSurface(gp.rp, ball, gp.geo.ol, PhiMath(), cRel);
         RubberContact::TubeProbe pr;
         pr.point = p.point + pv;
         pr.dist = p.dist;
         pr.normal = p.dist > 1e-12 ? (cRel - p.point) / p.dist : V3(0.0, 0.0, 1.0);
         pr.delta = p.delta;
         pr.overlap = p.overlap;
         pr.reff = p.reff;
         pr.rho = 0.0;
         pr.seg = 0;
         return pr;
      }

      void SurfaceVelocity(const V3 &point, V3 &v, V3 &w) const override
      {
         const V3 &pv = gp.geo.pivot;
         const double wz = OmegaZ();
         v = V3(-wz * (point.y - pv.y) + gp.geo.vPivotX, wz * (point.x - pv.x), 0.0);
         w = V3(0.0, 0.0, wz);
      }

      // The old loop's calls: an island with a bat never goes there (GroupPart)
      RubberContact::ContactResult ContactItem(int, SurfaceContact &, const V3 &, const V3 &, const V3 &, double) override { assert(false); return RubberContact::ContactResult(); }
      double BallStep(const V3 &, const V3 &, double hMax) const override { assert(false); return hMax; }
      void StepInIsland(double) override { assert(false); }

      // ---- the group's calls ----
      GroupStep::Part *GroupPart() override { return &gp; }
      void GroupProbe(int, const V3 &c, const double R, const double *const trial, GroupStep::Hit &h) const override { gp.ProbeRubber(c, R, trial, h); }
      const RubberContact::RubberParams *GroupLaw(int, double &layerM) const override
      {
         layerM = gp.layerM;
         return &gp.rp;
      }
      void CycleBegin() override { idle = false; }
      void TouchBegan(int, const int ballIndex, const double vn) override { events.push_back({ Event::Hit, ballIndex, std::max(0.0, vn) }); }   // (the engine's touch: the release gap is its hysteresis)

      // A bat no ball is near: the same group with no balls, asleep while nothing changes (a part at rest sleeps, ModelStandards 3c: the coil unchanged, the state the same over a whole cycle)
      int AdvanceAlone(double dt) override
      {
         if (dt <= 0.0)
            return 0;
         if (idle && gp.solenoid == idleSol && !gp.pulse.active)
            return 0;
         const double phi0 = gp.phi, xp0 = gp.xp, iHigh0 = gp.iHigh, iLow0 = gp.iLow;
         const int steps = GroupStep::AloneCycle(gp, dt, bookMs, failedSteps, aloneSteps, energySteps);
         // a whole cycle in which nothing moved that any ball or the table could see (angle, speed, plunger, both coil currents): the bat is at rest
         const double r = TipReach();
         idle = dt >= 0.9e-3 && !gp.hitStopLatch && !gp.pulse.active && std::fabs(gp.phi - phi0) * r < 1e-9 && std::fabs(gp.omega) * r < 1e-6 && std::fabs(gp.xp - xp0) < 1e-9 && std::fabs(gp.vp) < 1e-6
                && std::fabs(gp.iHigh - iHigh0) < 1e-9 && std::fabs(gp.iLow - iLow0) < 1e-9;
         idleSol = gp.solenoid;
         return steps;
      }

      // The dynamic state, plain bytes (the settings, the outline and the coil input are not state)
      struct Dyn
      {
         double xp, vp, phi, omega, iHigh, iLow, aPhi, hitL, pendingStopJ;
         FlipperModel::PivotPulse pulse;
         double steppedS;
         bool rigid, wasPreEos, coilOffPrev, plungerPinned, phiPinned, hitStop, hitStopLatch, idle, idleSol, steppedThisCycle;
      };
      void SaveState(State &out) const override
      {
         Dyn d;
         std::memset(&d, 0, sizeof d);   // (the padding too: saved bytes must be the same for the same state)
         const BatGroup::Part &b = gp;
         d.xp = b.xp; d.vp = b.vp; d.phi = b.phi; d.omega = b.omega; d.iHigh = b.iHigh; d.iLow = b.iLow; d.aPhi = b.aPhi; d.hitL = b.hitL; d.pendingStopJ = b.pendingStopJ;
         d.pulse = b.pulse; d.steppedS = steppedS;
         d.rigid = b.rigid; d.wasPreEos = b.wasPreEos; d.coilOffPrev = b.coilOffPrev; d.plungerPinned = b.plungerPinned; d.phiPinned = b.phiPinned; d.hitStop = b.hitStop; d.hitStopLatch = b.hitStopLatch;
         d.idle = idle; d.idleSol = idleSol; d.steppedThisCycle = steppedThisCycle;
         out.used = 0;
         Put(out, d);
      }
      void RestoreState(const State &in) override
      {
         StateReader r(in);
         Dyn d;
         r.Get(d);
         BatGroup::Part &b = gp;
         b.xp = d.xp; b.vp = d.vp; b.phi = d.phi; b.omega = d.omega; b.iHigh = d.iHigh; b.iLow = d.iLow; b.aPhi = d.aPhi; b.hitL = d.hitL; b.pendingStopJ = d.pendingStopJ;
         b.pulse = d.pulse; steppedS = d.steppedS;
         b.rigid = d.rigid; b.wasPreEos = d.wasPreEos; b.coilOffPrev = d.coilOffPrev; b.plungerPinned = d.plungerPinned; b.phiPinned = d.phiPinned; b.hitStop = d.hitStop; b.hitStopLatch = d.hitStopLatch;
         idle = d.idle; idleSol = d.idleSol; steppedThisCycle = d.steppedThisCycle;
         b.MovePivot();
         events.clear();
      }

      // The largest angular acceleration the bat can have in the coming cycle, from its present state (UnifiedContactDesign 12.2; measured 2026-10-28: the full-coil bound of
      // about 209,000 rad/s^2 made a held bat count as able to reach 14 m/s at its tip within a cycle):
      //  - coil on with the plunger at its stop: the stop takes the coil's push; the bat can only move within the link's play, against the return spring
      //  - coil off and already off at the last step: the windings only decay, so the coil's force is at most its present force (no flux jump: that happens at the switch-off)
      //  - anything else (a flip under way, the coil just switched): the full coil at the longest lever
      // The coil changes only between cycles (the table's solenoid), so the state read here holds for the whole cycle.
      double AlphaBound() const
      {
         const FlipperSettings &s = bat.set;
         const FlipperDerived &p = bat.der;
         const double inertia = std::max(p.iFlipper, 1e-9);
         const double xpMm = gp.xp * 1e3;
         if (gp.solenoid && xpMm <= p.xMin)
            return springMax / inertia;
         if (!gp.solenoid && gp.coilOffPrev)
         {
            const double frac = (xpMm - p.xMin) / s.bobbinMm;
            const double gHigh = (s.x0P / (frac + s.x0P)) * (s.x0P / (frac + s.x0P)), gLow = (s.x0H / (frac + s.x0H)) * (s.x0H / (frac + s.x0H));
            const double fCoil = p.strengthP * gHigh * gp.iHigh * gp.iHigh + p.strengthH * gLow * gp.iLow * gp.iLow;
            return (fCoil * leverMax + springMax) / inertia;
         }
         return alphaMax;
      }

      // The fastest a point of the rubber can move toward a ball within cycleS: the tip's present speed plus what the coil and the spring can add (AlphaBound), plus the pivot's speed
      double SurfaceSpeedBound(double cycleS) const override
      {
         return (std::fabs(gp.omega) + AlphaBound() * cycleS) * TipReach() + std::fabs(gp.geo.vPivotX);
      }

      std::string StateNote() const override { return "phi " + std::to_string(gp.phi) + " rad, omega " + std::to_string(gp.omega) + " rad/s, coil " + (gp.solenoid ? "on" : "off"); }

   private:
      bool idle = false, idleSol = false;
   };
}
