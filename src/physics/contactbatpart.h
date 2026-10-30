// license:GPLv3+
//
// A flipper bat as a moving part of the contact system (docs/UnifiedContactDesign.txt 12.2). The bat, its coil (plunger, link, power and hold windings), return spring and stops
// are ContactBat (contactbat.h, rothbauerw's real-data model); this part makes it one more ContactMovingPart, stepped inside the islands of the balls that touch it like a gate,
// a bumper or the turntable disc, and saved with the island's snapshot.
//
// Laws:
//  - ball contact: the bat rubber's law (RubberContact::Evaluate: Hertz on the rubber layer + damping + Mindlin + twist) on the bat outline at its present angle; the bat is a
//    moving surface (rotation about the pivot, plus the pivot's own motion during a nudge pulse).
//  - Newton's third law at every island substep: the torque of every ball on the bat is summed (tauBalls) and acts on the bat over the whole substep, in which the bat takes its
//    own steps (its link and stops are far stiffer than the ball contact: at most maxSubstepS).
//  - closing speed: the tip's present speed plus what the coil and spring can add in a cycle, plus the pivot's speed (SurfaceSpeedBound).
//  - the step a ball's contact needs: the reduced mass of the ball and the bat's effective mass at the contact (its inertia over the lever arm squared) with the rubber's
//    Hertz stiffness at the depth the relative closing speed gives (BallStep).
// Energy: only the coil (its work on the plunger: coilWorkJ), the spring's preload and the pivot nudge (a prescribed motion) add energy; the contact, the link's free play
// and the stops remove it.
// Pure SI maths with no VPX types: the same code runs in the offline tests (runs/2026-10-28-bat). Positions are world metres; the pivot is at z = 0 (the rubber heights
// rubber_z_lo_mm, rubber_height_mm are from the playfield surface, as the bat contact always had). Angles: the table's frame (start angle, 0 along -y), swung by dirSign.

#pragma once

#include <algorithm>
#include <cmath>

#include "physics/contactmovingpart.h"
#include "physics/contactbat.h"

namespace BatPart
{
   using V3 = RubberContact::V3;

   class Part : public ContactMovingPart
   {
   public:
      ContactBat bat;                    // settings and state of the bat, its coil and its stops (contactbat.h)
      V3 pivot;                          // pivot at rest, m (z = 0)
      RubberContact::Outline ol;         // bat outline: base radius, end radius, centre to centre, m
      double angleStart = 0.0;           // rad: the table's start angle (0 points along -y)
      double dirSign = 1.0;              // +1 / -1: the direction the bat swings toward its end stop
      RubberContact::RubberParams rp;    // the bat rubber (FlipperModel::MakeRubberParams of bat.set); Configure after the settings change
      double maxSubstepS = 20e-6;        // the bat's own step (coil, link, stops), at most
      double releaseGapM = 0.2e-3;       // a ball this far from the rubber counts as still touching for the Hit event (rubber_release_gap_mm)
      bool solenoid = false;             // input: the coil is energised (set by the table between cycles; not part of the saved state)
      bool hitStopLatch = false;         // the plunger reached the coil stop since the table last read it (EOS, the nudge of the other bat)
      double tauBalls = 0.0;             // torque of the balls in the present island substep, closing direction, N m
      RubberContact::BallSI ball;
      double alphaMax = 0.0;             // the largest angular acceleration the coil and the spring can give the bat, rad/s^2 (Configure)
      double leverMax = 0.0, springMax = 0.0;   // the crank's longest lever (m) and the return spring's largest torque (N m) over the stroke (Configure)

      // After a settings change: the rubber, the release gap and the acceleration bound
      void Configure()
      {
         rp = FlipperModel::MakeRubberParams(bat.set);
         releaseGapM = bat.set.rubberReleaseGapMm * 1e-3;
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
      }

      // the bat's angle in the outline's frame (rad from +x, counter-clockwise) and its angular velocity about +z
      double PhiMath() const { return angleStart + dirSign * bat.phi - M_PI / 2.0; }
      double OmegaZ() const { return dirSign * bat.omega; }
      // the pivot now (the nudge pulse moves it along x) and its velocity
      V3 PivotNow(double &vx) const
      {
         double d = 0.0, v = 0.0;
         bat.NudgeEval(bat.nzT, d, v);
         vx = v;
         return V3(pivot.x + d, pivot.y, pivot.z);
      }
      double TipReach() const { return ol.Lc + ol.re; }

      // ---- the moving-part interface (contactmovingpart.h) ----
      void SetBall(const RubberContact::BallSI &b) override { ball = b; }
      void SetGravity(const V3 &) override {}   // the bat turns in the playfield's plane about an axis normal to it: gravity does not turn it
      V3 Origin() const override { return pivot; }
      double Reach() const override { return TipReach(); }

      RubberContact::TubeProbe Probe(const V3 &c) const override
      {
         double vx;
         const V3 pv = PivotNow(vx);
         const V3 cRel = c - pv;
         const RubberContact::Probe p = RubberContact::ProbeSurface(rp, ball, ol, PhiMath(), cRel);
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
         double vx;
         const V3 pv = PivotNow(vx);
         const double wz = OmegaZ();
         v = V3(-wz * (point.y - pv.y) + vx, wz * (point.x - pv.x), 0.0);
         w = V3(0.0, 0.0, wz);
      }

      // One substep of one ball's contact: the rubber law in the pivot's frame (ball velocity relative to the moving pivot); the torque on the bat is kept for StepInIsland.
      // A Hit event when a contact starts (it.partTouching: the release gap is the hysteresis), with the closing speed of ball and surface.
      RubberContact::ContactResult ContactItem(int ballIndex, SurfaceContact &it, const V3 &c, const V3 &v, const V3 &w, double h) override
      {
         double vx;
         const V3 pv = PivotNow(vx);
         const V3 cRel = c - pv, vRel = v - V3(vx, 0.0, 0.0);
         const RubberContact::ContactResult r = RubberContact::Evaluate(rp, ball, ol, PhiMath(), OmegaZ(), cRel, vRel, w, h, it.st);
         if (r.inContact)
         {
            tauBalls += dirSign * r.torqueOnBat;   // (Evaluate gives it about +z; the bat's own frame closes along dirSign)
            if (!it.partTouching)
            {
               it.partTouching = true;
               V3 vs, ws;
               SurfaceVelocity(r.point + pv, vs, ws);
               events.push_back({ Event::Hit, ballIndex, std::max(0.0, -(v - vs).Dot(r.normal)) });
            }
         }
         else if (it.partTouching && r.dist - ball.R > releaseGapM)
            it.partTouching = false;
         return r;
      }

      // The longest substep this ball's contact allows: the stiffness the rubber's law has (Hertz with the thin bonded layer's stiffening, ContactCore: Fe = 4/3 E* a^3 / reff
      // (1 + beta (a/h)^3), a = sqrt(reff delta), so dFe/ddelta = 2 E* sqrt(reff delta) (1 + 2 beta (a/h)^3)) at the depth the relative closing speed reaches, on the reduced
      // mass of the ball and the bat's effective mass at the contact; omega h <= 0.05.
      // Why 0.05 (measured offline, runs/2026-10-28-bat): the bat takes the balls' torque of the substep's start state over the whole substep (the coupling of the island is first
      // order: its error grows with omega h, not its square). A ball the coil drives the bat back into (case 2, 1 m/s) needed omega h <= 0.05: at 0.1 the exit was 17 % off the
      // 1 us reference, at 0.05 0.35 % (report_scale05.txt); a predictor-corrector coupling at 0.1 still left 2.7 % at twice the force evaluations (report_pc.txt). Plain Hertz
      // with omega h <= 0.2 had let a bat contact step at 100 us, 9.5 % off (diag_report.txt).
      double BallStep(const V3 &c, const V3 &v, double hMax) const override
      {
         const RubberContact::TubeProbe pr = Probe(c);
         if (!pr.overlap && pr.dist - ball.R > 0.002)
            return hMax;   // not near: the arrival rule of the contact system bounds the step
         V3 vs, ws;
         SurfaceVelocity(pr.point, vs, ws);
         const double closing = std::max(std::fabs((v - vs).Dot(pr.normal)), 0.02);
         double vx;
         const V3 pv = PivotNow(vx);
         const V3 arm = pr.point - pv;
         const double lever = std::max(std::fabs(arm.x * pr.normal.y - arm.y * pr.normal.x), 1e-3);   // moment arm of the normal force about the pivot
         // the bat's mass at the contact in the contact's fast motion: the bat's own inertia; the plunger only when it is slaved to the bat (rigid: coil off, link taut). While
         // the coil drives, the plunger sits behind the link's spring and does not follow the bat within the contact's period.
         const double mBat = (bat.rigid ? bat.InertiaSI() : bat.der.iFlipper) / (lever * lever);
         const double mRed = ball.m * mBat / (ball.m + mBat);
         const double Estar = rp.eMpa * 1e6 / (1.0 - rp.nu * rp.nu);
         const double reff = std::max(pr.reff, 1e-4);
         const double dImpact = std::pow(15.0 * mRed * closing * closing / (16.0 * Estar * std::sqrt(reff)), 0.4);
         const double delta = std::max({ pr.overlap ? pr.delta : 0.0, dImpact, 1e-6 });
         const double ah = std::sqrt(reff * delta) / std::max(rp.thicknessMm * 1e-3, 1e-6);
         const double k = 2.0 * Estar * std::sqrt(reff * delta) * (1.0 + 2.0 * rp.betaLayer * ah * ah * ah);
         const double omega = std::sqrt(k / mRed);
         return std::min(hMax, 0.05 * ContactOmegaHScale() / omega);
      }

      // the bat's own steps within h, with the balls' torque held over all of them
      void Advance(double h, double tau)
      {
         const int n = std::max(1, (int)std::ceil(h / maxSubstepS - 1e-9));
         const double hs = h / n;
         for (int i = 0; i < n; i++)
         {
            bat.Step(hs, tau, solenoid);
            hitStopLatch = hitStopLatch || bat.plungerHitStop;
            AdvanceNudge(hs);
         }
      }
      void AdvanceNudge(double h)
      {
         if (!bat.nzActive)
            return;
         bat.nzT += h;
         if (bat.nzT >= bat.nzRiseS + bat.nzFallS)
            bat.nzActive = false;
      }

      void StepInIsland(double h) override
      {
         bat.idle = false;
         Advance(h, tauBalls);
         tauBalls = 0.0;
         steppedThisCycle = true;
         steppedS += h;
      }

      // A bat no ball is near: its own steps (as the flipper's mover did), asleep while nothing changes (dormancy: the coil unchanged, the state the same over a whole cycle)
      int AdvanceAlone(double dt) override
      {
         tauBalls = 0.0;
         if (dt <= 0.0)
            return 0;
         if (bat.IdleMatches(solenoid) && !bat.nzActive)
            return 0;
         const double phi0 = bat.phi, omega0 = bat.omega, xp0 = bat.xp, vp0 = bat.vp, iHigh0 = bat.iHigh, iLow0 = bat.iLow;
         const bool latch0 = hitStopLatch;
         hitStopLatch = false;
         const int n = std::max(1, (int)std::ceil(dt / maxSubstepS - 1e-9));
         Advance(dt, 0.0);
         const bool hitStop = hitStopLatch;
         hitStopLatch = hitStopLatch || latch0;
         // a whole cycle in which nothing moved (angle, speed, plunger, both coil currents) and the coil did not change: the bat is at rest
         bat.idle = dt >= 0.9e-3 && !hitStop && std::fabs(bat.phi - phi0) < 1e-12 && std::fabs(bat.omega - omega0) < 1e-9 && std::fabs(bat.xp - xp0) < 1e-12 && std::fabs(bat.vp - vp0) < 1e-9
                    && std::fabs(bat.iHigh - iHigh0) < 1e-12 && std::fabs(bat.iLow - iLow0) < 1e-12;
         if (bat.idle)
            bat.IdleCapture(solenoid);
         return n;
      }

      // The dynamic state, plain bytes (the settings, the outline and the solenoid input are not state)
      struct Dyn
      {
         double xp, vp, phi, omega, iHigh, iLow, lastAlpha, hitL, coilWorkJ;
         double idlePhi, idleOmega, idleXp, idleVp, idleIHigh, idleILow;
         double nzT, nzA, nzRiseS, nzFallS, nzBaseX;
         double tauBalls, steppedS;
         float nzVelU;
         bool rigid, wasPreEos, coilOffPrev, init, plungerHitStop, idle, idleSol, idleRigid, nzActive, hitStopLatch, steppedThisCycle;
      };
      void SaveState(State &out) const override
      {
         Dyn d;
         std::memset(&d, 0, sizeof d);   // (the padding too: saved bytes must be the same for the same state)
         const ContactBat &b = bat;
         d.xp = b.xp; d.vp = b.vp; d.phi = b.phi; d.omega = b.omega; d.iHigh = b.iHigh; d.iLow = b.iLow; d.lastAlpha = b.lastAlpha; d.hitL = b.hitL; d.coilWorkJ = b.coilWorkJ;
         d.idlePhi = b.idlePhi; d.idleOmega = b.idleOmega; d.idleXp = b.idleXp; d.idleVp = b.idleVp; d.idleIHigh = b.idleIHigh; d.idleILow = b.idleILow;
         d.nzT = b.nzT; d.nzA = b.nzA; d.nzRiseS = b.nzRiseS; d.nzFallS = b.nzFallS; d.nzBaseX = b.nzBaseX; d.nzVelU = b.nzVelU;
         d.tauBalls = tauBalls; d.steppedS = steppedS;
         d.rigid = b.rigid; d.wasPreEos = b.wasPreEos; d.coilOffPrev = b.coilOffPrev; d.init = b.init; d.plungerHitStop = b.plungerHitStop; d.idle = b.idle; d.idleSol = b.idleSol;
         d.idleRigid = b.idleRigid; d.nzActive = b.nzActive; d.hitStopLatch = hitStopLatch; d.steppedThisCycle = steppedThisCycle;
         out.used = 0;
         Put(out, d);
      }
      void RestoreState(const State &in) override
      {
         StateReader r(in);
         Dyn d;
         r.Get(d);
         ContactBat &b = bat;
         b.xp = d.xp; b.vp = d.vp; b.phi = d.phi; b.omega = d.omega; b.iHigh = d.iHigh; b.iLow = d.iLow; b.lastAlpha = d.lastAlpha; b.hitL = d.hitL; b.coilWorkJ = d.coilWorkJ;
         b.idlePhi = d.idlePhi; b.idleOmega = d.idleOmega; b.idleXp = d.idleXp; b.idleVp = d.idleVp; b.idleIHigh = d.idleIHigh; b.idleILow = d.idleILow;
         b.nzT = d.nzT; b.nzA = d.nzA; b.nzRiseS = d.nzRiseS; b.nzFallS = d.nzFallS; b.nzBaseX = d.nzBaseX; b.nzVelU = d.nzVelU;
         tauBalls = d.tauBalls; steppedS = d.steppedS;
         b.rigid = d.rigid; b.wasPreEos = d.wasPreEos; b.coilOffPrev = d.coilOffPrev; b.init = d.init; b.plungerHitStop = d.plungerHitStop; b.idle = d.idle; b.idleSol = d.idleSol;
         b.idleRigid = d.idleRigid; b.nzActive = d.nzActive; hitStopLatch = d.hitStopLatch; steppedThisCycle = d.steppedThisCycle;
         events.clear();
      }

      // The largest angular acceleration the bat can have in the coming cycle, from its present state (UnifiedContactDesign 12.2; measured 2026-10-28: the full-coil bound of
      // about 209,000 rad/s^2 made a held bat count as able to reach 14 m/s at its tip within a cycle, and the balls piled on it took 130-146 arrival substeps a cycle):
      //  - coil on with the plunger at its stop: the stop takes the coil's push; the bat can only move within the link's play, against the return spring
      //  - coil off and already off at the last step: the windings only decay, so the coil's force is at most its present force (no flux jump: that happens at the switch-off)
      //  - anything else (a flip under way, the coil just switched): the full coil at the longest lever
      // The coil changes only between cycles (the table's solenoid), so the state read here holds for the whole cycle.
      double AlphaBound() const
      {
         const FlipperSettings &s = bat.set;
         const FlipperDerived &p = bat.der;
         const double inertia = std::max(p.iFlipper, 1e-9);
         if (solenoid && bat.xp <= p.xMin)
            return springMax / inertia;
         if (!solenoid && bat.coilOffPrev)
         {
            const double frac = (bat.xp - p.xMin) / s.bobbinMm;
            const double gHigh = (s.x0P / (frac + s.x0P)) * (s.x0P / (frac + s.x0P)), gLow = (s.x0H / (frac + s.x0H)) * (s.x0H / (frac + s.x0H));
            const double fCoil = p.strengthP * gHigh * bat.iHigh * bat.iHigh + p.strengthH * gLow * bat.iLow * bat.iLow;
            return (fCoil * leverMax + springMax) / inertia;
         }
         return alphaMax;
      }

      // The fastest a point of the rubber can move toward a ball within cycleS: the tip's present speed plus what the coil and the spring can add (AlphaBound), plus the
      // pivot's speed
      double SurfaceSpeedBound(double cycleS) const override
      {
         double vx;
         PivotNow(vx);
         return (std::fabs(bat.omega) + AlphaBound() * cycleS) * TipReach() + std::fabs(vx);
      }

      std::string StateNote() const override { return "phi " + std::to_string(bat.phi) + " rad, omega " + std::to_string(bat.omega) + " rad/s, coil " + (solenoid ? "on" : "off"); }
   };
}
