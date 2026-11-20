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

   // what the engine side keeps (PartRecords: part of its state)
   struct EngineState
   {
      bool idle = false;
      GroupStep::RestWatch rest;            // the rest rule's watch, alone and over its last cycle in a group (RestPose)
   };

   class Part : public ContactMovingPart, public EngineState
   {
   public:
      Part() { KeepState<EngineState>(this); }
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

      // ---- the group's calls ----
      GroupStep::Part *GroupPart() override { return &gp; }
      void GroupProbe(int, const V3 &c, const double R, const double *const trial, GroupStep::Hit &h) const override { gp.ProbeRubber(c, R, trial, h); }
      double GroupPathMotion(int, const V3 &P, const double dLen, const double *const trial) const override { return gp.PathMotion(P, dLen, trial); }
      double RegionDist(const V3 &c) const override { return gp.RegionDist(c); }   // (group A: its bounding radius about the pivot)
      const RubberContact::RubberParams *GroupLaw(int, double &layerM) const override
      {
         layerM = gp.layerM;
         return &gp.rp;
      }
      void CycleBegin() override { idle = false; }
      // the one rest rule (GroupStep::RestWatch), as every part: driven while its coil or the pivot's pulse is on; the bat (at its tip) and the plunger still, and both windings' currents
      // died away (fractions of the held coil's current: full scale 1)
      bool Sleeping() const override { return idle; }
      GroupStep::RestWatch *Rest() override { return &rest; }
      GroupStep::RestWatch::Pose RestPose() const override { return GroupStep::RestWatch::Pose().Add(gp.phi, TipReach()).Add(gp.xp, 1.0).State(gp.iHigh, 1.0).State(gp.iLow, 1.0); }
      bool AtRest() const override
      {
         if (gp.solenoid || gp.pulse.active)
            return false;
         return idle || rest.Rests(RestPose(), RestSpeed());
      }
      void TouchBegan(int, const int ballIndex, const double vn) override { events.push_back({ Event::Hit, ballIndex, std::max(0.0, vn) }); }   // (the engine's touch: the release gap is its hysteresis)

      // A bat no ball is near: the same group with no balls, asleep while nothing changes (a part at rest sleeps, ModelStandards 3c: the coil unchanged, the state the same over a whole cycle)
      int AdvanceAlone(double dt) override
      {
         if (dt <= 0.0)
            return 0;
         if (idle && !gp.solenoid && !gp.pulse.active)   // (asleep only with the coil off: its coil or a pulse wakes it)
            return 0;
         const bool was = idle;
         const int steps = GroupStep::AloneCycle(gp, dt, bookMs, failedSteps, aloneSteps, energySteps);
         idle = rest.Step(RestPose(), RestSpeed(), gp.solenoid || gp.pulse.active, restTouched, dt);
         GroupStep::RestWatch::Book(3, was, idle);
         return steps;
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

      // the energy reach (design_reach.md section 1): the bat turns about the vertical axis through its pivot by at most EnergyReach(|omega|, E, I, AlphaBound) (the coil and the
      // spring as the drive, AlphaBound's present state holds for the cycle), at most pi; a point c moves relative to it by the chord at c's distance from that axis, plus the pivot's own
      // slide |vPivotX| T (given, as SurfaceSpeedBound reads it). The turn is bounded by the energy (EnergyMove), not clamped at the stops; work: as EnergyMove
      double CycleRegionDist(const V3 &c, const double T, double) const override { return gp.RegionDist(c) - std::fabs(gp.geo.vPivotX) * T; }   // (the energy reach: its bounding radius about the pivot at any angle, the pivot's slide off)
      // the energy reach, the coil law (runs/2026-11-14-reach-p1b/brief_coil.md item 2): a winding's current cannot jump: over T, i <= min(i0 + (V / L) T, V / R), as a
      // fraction of V / R: i0 + (R / L) T, L its smallest (the plunger farthest out, at rest). The solenoid changes only between cycles: on, both windings rise from their present
      // currents; off, they only decay from their present currents (just switched off, the hold winding carries the series loop's current, which it already holds: no swap).
      // The force at those currents with the plunger closed (the largest gap factor, 1), at the longest lever; with the coil on and the plunger at its stop, the stop takes the push
      bool Parked() const { return gp.solenoid && gp.xp <= bat.der.xMin * 1e-3; }   // coil on, the plunger at its stop (BatGroup pins xp at xMin in metres: the test holds at equality)
      double CoilAlphaBound(const double T) const
      {
         const FlipperSettings &s = bat.set;
         const FlipperDerived &p = bat.der;
         const double inertia = std::max(p.iFlipper, 1e-9);
         if (Parked())
            return springMax / inertia;
         double iP = gp.iHigh, iH = gp.iLow;
         if (gp.solenoid)
         {
            const double fracOut = std::max((p.xRest - p.xMin) / s.bobbinMm, 0.0);
            const double lP = (p.aP + p.bP / (fracOut + s.x0P)) * 1e-3, lH = (p.aH + p.bH / (fracOut + s.x0H)) * 1e-3;
            iP = std::min(iP + s.rP / std::max(lP, 1e-9) * T, std::max(iP, 1.0));
            iH = s.single ? 0.0 : std::min(iH + s.rH / std::max(lH, 1e-9) * T, std::max(iH, 1.0));
         }
         const double fCoil = p.strengthP * iP * iP + p.strengthH * iH * iH;
         return (fCoil * leverMax + springMax) / inertia;
      }
      // the energy its drives can supply over T (the law: a part's work over a cycle is at most its kinetic energy plus this), J: the coil's field energy now
      // 1/2 (L_P i_P^2 + L_H i_H^2) (amperes, L at its largest: the plunger closed) plus, the coil on, the supply's V^2 (1/R_P + 1/R_H) T, plus the return spring's stored energy
      // (its largest torque x the larger of the stroke and the turn from rest). Parked (coil on at its stop): the stop takes the coil, the spring only
      double DriveEnergy(const double T) const
      {
         const FlipperSettings &s = bat.set;
         const FlipperDerived &p = bat.der;
         const double spring = springMax * std::max(gp.PhiStop(), std::fabs(gp.phi));
         if (Parked())
            return spring;
         const double lP = (p.aP + p.bP / s.x0P) * 1e-3, lH = s.single ? 0.0 : (p.aH + p.bH / s.x0H) * 1e-3;
         const double aP = gp.iHigh * s.vSupply / s.rP, aH = s.single ? 0.0 : gp.iLow * s.vSupply / s.rH;
         const double supply = gp.solenoid ? s.vSupply * s.vSupply * (1.0 / s.rP + (s.single ? 0.0 : 1.0 / s.rH)) * T : 0.0;
         return 0.5 * (lP * aP * aP + lH * aH * aH) + supply + spring;
      }
      // the energy reach (design_reach.md section 1): the turn EnergyReach(|omega|, E, I, CoilAlphaBound), and at most the turn a bat with all its energy (its kinetic energy, the feed and
      // DriveEnergy) as speed makes in T: sqrt(w0^2 + 2 (E + D) / I) T (not clamped at the stops: past them the stops' own law only shortens it). work: 1/2 I w0^2 + the drives' work
      // over the turn, at most DriveEnergy
      double EnergyMove(const V3 &c, const double E, const double T, double *const work) const override
      {
         const double inertia = std::max(bat.der.iFlipper, 1e-9), alpha = CoilAlphaBound(T), w0 = std::fabs(gp.omega), D = DriveEnergy(T);
         const double dphi = DriveCappedReach(w0, E, inertia, alpha, D, T);
         if (work != nullptr)
            *work = 0.5 * inertia * w0 * w0 + std::min(inertia * alpha * dphi, D);
         const V3 &pv = gp.geo.pivot;
         const double rho = std::hypot(c.x - pv.x, c.y - pv.y);
         return 2.0 * rho * std::sin(0.5 * std::min(dphi, GroupStep::kPi)) + std::fabs(gp.geo.vPivotX) * T;
      }

      std::string StateNote() const override { return "phi " + std::to_string(gp.phi) + " rad, omega " + std::to_string(gp.omega) + " rad/s, coil " + (gp.solenoid ? "on" : "off"); }

   private:
      double RestSpeed() const { return std::fabs(gp.omega) * TipReach() + std::fabs(gp.vp) + std::fabs(gp.geo.vPivotX); }   // the bat's and plunger's surface speed bound, m/s
   };
}
