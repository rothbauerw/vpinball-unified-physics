// license:GPLv3+
//
// A turntable disc as a moving part of the contact system (docs/UnifiedContactDesign.txt section 12.1): a flat disc set into the playfield, turned about its vertical axis by
// a DC motor, held back by its bearing, and turned by the balls on it. Islands that hold the disc use the shared arrival-time loop (SubstepsArrival) like a gate or a bumper.
//
// Laws:
//  - motor and bearing (unchanged from turntable.cpp): powered, I dw/dt = Ts V / Vrated - kTot w (the motor's torque-speed line, the viscous losses included in kTot);
//    unpowered, I dw/dt = -b w (bearing drag); both solved exactly over a step (exponential). A constant (Coulomb) friction torque frictionNm slows a turning disc and holds
//    a still one until the push of the step exceeds it.
//  - balls: the contact law of every surface (EvaluateTube: Hertz + damping + Mindlin + twist + rolling resistance) with the surface moving at w x r. The equal and opposite
//    torque about the axis acts on the disc in the same substep (Newton's third law at every substep, StepInIsland), not once per cycle.
//  - closing speed: a disc turning about its own axis maps its surfaces (top face, side) onto themselves; no point of them moves along its own normal, so the disc never
//    closes on a ball by its motion (SurfaceSpeedBound = 0).
// Energy: only the motor adds energy (its torque times the disc speed); the bearing and the contact remove it.
// Pure SI maths with no VPX types: the same code runs in the offline tests (runs/2026-10-28-disc). Angles and speeds are clockwise seen from above positive, the fixed frame of
// Turntables (TurntableSpeed, TurntableAngle); the sign of the voltage's effect is dir.

#pragma once

#include <algorithm>
#include <cmath>

#include "physics/contactmovingpart.h"

namespace TurntableDisc
{
   using V3 = RubberContact::V3;

   // The motor and bearing constants worked out from the ini keys (turntable.cpp Parse)
   struct Motor
   {
      // inputs
      double I = 0.0;            // kg m^2
      double freeRpm = 250.0;    // speed at ratedVolts with no load
      double ratedVolts = 12.0;
      double spinupS = 1.0;      // rest to 95 % of freeRpm at ratedVolts, no load
      double rampDownS = 4.0;    // freeRpm to 5 % with the motor off, no load
      double stallTorque = 0.0;  // N m at ratedVolts; 0 = from spinupS
      double frictionNm = 0.02;  // Coulomb torque of the bearing and gear train
      double dir = 1.0;          // +1: a positive voltage turns the disc clockwise seen from above
      // derived
      double Ts = 0.0, kTot = 0.0, b = 0.0, wRun = 0.0;
      double coastFrictionS = 0.0;   // how long the friction alone takes from freeRpm to 5 %, s (for the log)
      bool frictionFaster = false;   // the friction alone stops the disc sooner than rampDownS
      bool dragClamped = false;      // rampDownS asked for more drag than the motor's own slope: b = kTot

      void Derive()
      {
         constexpr double PI = 3.14159265358979323846;
         wRun = freeRpm * 2.0 * PI / 60.0;
         if (stallTorque > 0)
         {
            Ts = stallTorque;
            kTot = std::max(1e-9, (Ts - frictionNm) / wRun);
         }
         else
         {
            kTot = I / (spinupS / 3.0);   // 95 % = 3 time constants
            Ts = kTot * wRun + frictionNm;   // the running speed at ratedVolts is freeRpm with the friction included
         }
         // viscous drag b such that the unpowered disc, with the friction too, goes from freeRpm to 5 % in rampDownS:  I dw/dt = -b w - Tf
         const double w0 = wRun, Tf = frictionNm;
         const auto coast = [&](double bb) {
            if (Tf <= 0.0)
               return I / bb * std::log(20.0);
            if (bb < 1e-12)
               return 0.95 * w0 * I / Tf;
            return I / bb * std::log((w0 + Tf / bb) / (0.05 * w0 + Tf / bb));
         };
         double lo = 1e-9, hi = 1e3 * I;
         coastFrictionS = coast(1e-9);
         frictionFaster = coastFrictionS < rampDownS;
         if (frictionFaster)
            hi = lo;   // the friction alone stops the disc faster than rampDownS
         for (int i = 0; i < 80 && hi > lo; i++)
         {
            const double mid = std::sqrt(lo * hi);
            if (coast(mid) > rampDownS)
               lo = mid;
            else
               hi = mid;
         }
         b = 0.5 * (lo + hi);
         dragClamped = b > kTot;
         if (dragClamped)
            b = kTot;
      }
   };

   // Where the disc is and how a ball touches it
   struct Shape
   {
      V3 centre;                // centre of the top face, m
      double R = 0.0762;        // disc radius, m
      double zBottom = 0.0;     // bottom of the side when the disc stands proud of the playfield (= centre.z when it is flush: a disc with no side), m
      double edgeRhoM = 1e-3;   // the rim edge counts as rounded to this radius (Hertz effective radius, as ProbeTriangles' edgeRho)
      RubberContact::RubberParams rp;   // the textured surface against the steel ball
      double hLayer = 0.0054;   // layer thickness of the contact law, m
   };

   class Part : public ContactMovingPart
   {
   public:
      Motor motor;
      Shape shape;
      // state
      double volts = 0.0;       // the drive the script set (Turntable / TurntableDuty); an input, not changed by the step
      double omega = 0.0;       // rad/s, clockwise seen from above positive
      double theta = 0.0;       // rad, 0 .. 2 pi
      double tauBalls = 0.0;    // torque of the balls about the axis in the present substep (clockwise positive), N m; taken by StepInIsland
      RubberContact::BallSI ball;

      // One step of dt seconds with the balls' torque tau held over it: the motor law exactly, then the balls' angular impulse, then the Coulomb friction (as turntable.cpp did once per cycle)
      void Advance(const double dt, const double tau)
      {
         constexpr double PI = 3.14159265358979323846;
         const double vEff = volts * motor.dir;
         double wInf = 0.0, k = motor.b;
         if (vEff != 0.0)
         {
            wInf = motor.Ts * (vEff / motor.ratedVolts) / motor.kTot;
            k = motor.kTot;
         }
         const double wOld = omega;
         omega = wInf + (omega - wInf) * std::exp(-k * dt / motor.I) + tau * dt / motor.I;
         const double dwf = motor.frictionNm * dt / motor.I;   // the friction holds a still disc until the push of the step exceeds it
         omega = omega > 0.0 ? std::max(0.0, omega - dwf) : std::min(0.0, omega + dwf);
         theta += 0.5 * (wOld + omega) * dt;
         theta = std::fmod(theta, 2.0 * PI);
         if (theta < 0.0)
            theta += 2.0 * PI;
      }
      bool Asleep() const { return volts == 0.0 && omega == 0.0; }

      // ---- the moving-part interface (contactmovingpart.h) ----
      void SetBall(const RubberContact::BallSI &b) override { ball = b; }
      void SetGravity(const V3 &) override {}   // the disc turns about its own axis: its weight acts on the bearing, not about the axis
      V3 Origin() const override { return shape.centre; }
      double Reach() const override { return shape.R; }

      // The nearest point of the disc (a solid cylinder from zBottom up to the top face; a flush disc is the top face alone) to a ball centred at c. On the rim edge the
      // contact counts as rounded to edgeRhoM as on any face edge (the effective radius falls with the angle from the nearer face, fully rounded at 15 deg).
      RubberContact::TubeProbe Probe(const V3 &c) const override
      {
         RubberContact::TubeProbe pr;
         const double zTop = shape.centre.z, zBot = std::min(shape.zBottom, zTop);
         const double dx = c.x - shape.centre.x, dy = c.y - shape.centre.y;
         const double r = std::sqrt(dx * dx + dy * dy);
         const V3 radial = r > 1e-12 ? V3(dx / r, dy / r, 0.0) : V3(1.0, 0.0, 0.0);
         const V3 up(0.0, 0.0, 1.0);
         V3 p, nFace;
         bool edge = false;
         if (r <= shape.R && c.z <= zTop && c.z >= zBot && zTop > zBot)
         {
            // inside the solid: out through the nearer of the top face and the side
            const double toTop = zTop - c.z, toSide = shape.R - r;
            if (toTop <= toSide)
            {
               p = V3(c.x, c.y, zTop);
               nFace = up;
               pr.dist = -toTop;
            }
            else
            {
               p = V3(shape.centre.x + radial.x * shape.R, shape.centre.y + radial.y * shape.R, c.z);
               nFace = radial;
               pr.dist = -toSide;
            }
            pr.normal = nFace;
         }
         else
         {
            const double rc = std::min(r, shape.R), zc = std::min(std::max(c.z, zBot), zTop);
            p = V3(shape.centre.x + radial.x * rc, shape.centre.y + radial.y * rc, zc);
            const V3 d = c - p;
            pr.dist = d.Length();
            if (r <= shape.R)
               nFace = c.z >= zTop ? up : up * -1.0;   // over (or, for a flush disc, under) the face: a thin sheet, as the faces it replaces
            else if (c.z <= zTop && c.z >= zBot && zTop > zBot)
               nFace = radial;   // beside the side of a proud disc
            else
            {
               edge = true;   // past the rim and above (or below) the face: the edge
               nFace = std::abs(radial.Dot(d)) > std::abs(up.Dot(d)) && zTop > zBot ? radial : (c.z >= zTop ? up : up * -1.0);
            }
            pr.normal = (edge && pr.dist > 1e-12) ? d / pr.dist : nFace;
         }
         pr.point = p;
         pr.seg = 0;
         pr.rho = 0.0;
         pr.delta = ball.R - pr.dist;
         pr.overlap = pr.delta > 0.0;
         double edgeW = 0.0;   // 0 flat .. 1 fully rounded
         if (edge)
         {
            const double ang = std::acos(std::min(1.0, std::max(-1.0, pr.normal.Dot(nFace))));
            edgeW = std::min(1.0, std::max(0.0, (ang - 5.0 * M_PI / 180.0) / (10.0 * M_PI / 180.0)));
         }
         const double R1 = 1.0 / ((1.0 - edgeW) / ball.R + edgeW * (1.0 / ball.R + 1.0 / shape.edgeRhoM));
         pr.reff = std::sqrt(R1 * ball.R);
         return pr;
      }

      void SurfaceVelocity(const V3 &point, V3 &v, V3 &w) const override
      {
         v = V3(-omega * (point.y - shape.centre.y), omega * (point.x - shape.centre.x), 0.0);
         w = V3(0.0, 0.0, omega);
      }

      // One substep of one ball's contact: the shared law with the disc's surface velocity; the torque of the reaction about the axis is kept for StepInIsland
      RubberContact::ContactResult ContactItem(int, SurfaceContact &it, const V3 &c, const V3 &v, const V3 &w, double h) override
      {
         const RubberContact::TubeProbe pr = Probe(c);
         V3 vs, ws;
         if (pr.overlap)
            SurfaceVelocity(pr.point, vs, ws);
         const RubberContact::ContactResult r = RubberContact::EvaluateTube(shape.rp, ball, pr, shape.hLayer, v, w, h, it.st, vs, ws);
         if (r.inContact)
            tauBalls += ReactionTorque(pr, r);
         return r;
      }

      // The torque about the axis (clockwise positive) of the equal and opposite of a contact: the force at the contact point and the couples (twist, rolling resistance)
      double ReactionTorque(const RubberContact::TubeProbe &pr, const RubberContact::ContactResult &r) const
      {
         const V3 rB = r.normal * (-(ball.R - r.delta / 2.0));
         const V3 Mc = r.torqueOnBall - rB.Cross(r.force);
         return -((pr.point.x - shape.centre.x) * r.force.y - (pr.point.y - shape.centre.y) * r.force.x) - Mc.z;
      }

      // The ball's normal approach: the disc is held along its axis by its bearing and does not move along any normal, so the ball's contact needs nothing beyond the shared rule
      double BallStep(const V3 &, const V3 &, double hMax) const override { return hMax; }

      void StepInIsland(double h) override
      {
         Advance(h, tauBalls);
         tauBalls = 0.0;
         steppedThisCycle = true;
         steppedS += h;
      }

      int AdvanceAlone(double dt) override
      {
         tauBalls = 0.0;
         if (Asleep() || dt <= 0.0)
            return 0;
         Advance(dt, 0.0);
         return 1;
      }

      void SaveState(State &out) const override
      {
         out.used = 0;
         Put(out, steppedThisCycle);
         Put(out, steppedS);
         Put(out, omega);
         Put(out, theta);
         Put(out, tauBalls);
      }
      void RestoreState(const State &in) override
      {
         StateReader rd(in);
         rd.Get(steppedThisCycle);
         rd.Get(steppedS);
         rd.Get(omega);
         rd.Get(theta);
         rd.Get(tauBalls);
         events.clear();
      }

      double SurfaceSpeedBound(double) const override { return 0.0; }   // no point of the disc moves along its own normal (see the header)

      std::string StateNote() const override { return "omega " + std::to_string(omega) + " rad/s, volts " + std::to_string(volts); }
   };
}
