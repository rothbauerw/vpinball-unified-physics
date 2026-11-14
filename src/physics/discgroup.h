// license:GPLv3+
//
// A turntable disc as a part of the group step (docs/unified-physics/GroupStepDesign.txt 21; runs/2026-11-02-disc): a flat disc set into the playfield, turned about its vertical axis by a DC motor, held back by its bearing, turned by the balls on it. One unknown: the turn over a step (bodypart.h).
//   motor    powered, I dw/dt = Ts V / Vrated - kTot w - Tf sign(w): the torque-speed line as a held torque whose work is counted (DriveForce) and its slope as a viscous loss (Viscous). Energy enters only here.
//   bearing  unpowered, I dw/dt = -b w - Tf sign(w). Tf: the Coulomb friction, BodyPart's stick rule: a still disc stays until the other forces exceed it.
//   balls    the contact law of every surface with the surface moving at w x r: a carried surface of the group (the Hit's column), its reaction on the disc in the same solve.
// No stored energy, no stops. A disc turning about its own axis maps its surfaces onto themselves: no depth changes with the turn (slope 0), so the disc never closes on a ball by its motion, and its turn is on no shape.
// Pure SI, no VPX types. Angles and speeds are clockwise seen from above positive, the frame of Turntables (TurntableSpeed, TurntableAngle); the sign of the voltage's effect is dir.

#pragma once

#include <algorithm>
#include <cmath>

#include "physics/bodypart.h"

namespace DiscGroup
{
   using GroupStep::Hit;
   using GroupStep::V3;

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
      // the inertia is not a key (user 2026-11-04): the friction alone gives the spin-down, I = Tslide rampDownS / (0.95 freeRpm), at least the disc's own (a steel plate of radius radiusM)
      double dynRatio = 0.85;        // [Contact] mu_dynamic_ratio: the disc slides at dynRatio frictionNm (bodypart.h PinStep)
      double radiusM = 0.0762;
      static constexpr double kPlateThickM = 3.175e-3, kSteelRho = 7850.0;   // the disc's own inertia, the floor: 1/8 in steel
      double tSlide = 0.0;           // the bearing's sliding torque, N m
      bool floorApplied = false;     // the friction alone stops the plate's own inertia sooner than rampDownS: the drag makes up the rest

      // the viscous drag b such that the unpowered disc, with the Coulomb torque Tf too, goes from freeRpm to 5 % in rampDownS:  I dw/dt = -b w - Tf
      void SolveDrag(const double Tf)
      {
         const double w0 = wRun;
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

      void Derive()
      {
         // the motor's line holds the disc at freeRpm against the sliding torque, the speed it runs at; the friction alone coasts it to 5 % in rampDownS
         wRun = freeRpm * 2.0 * GroupStep::kPi / 60.0;
         tSlide = dynRatio * frictionNm;
         const double plate = 0.5 * kSteelRho * GroupStep::kPi * radiusM * radiusM * kPlateThickM * radiusM * radiusM;
         I = tSlide * rampDownS / (0.95 * wRun);
         floorApplied = I < plate;
         if (floorApplied)
            I = plate;
         if (stallTorque > 0)
         {
            Ts = stallTorque;
            kTot = std::max(1e-9, (Ts - tSlide) / wRun);
         }
         else
         {
            kTot = I / (spinupS / 3.0);   // 95 % = 3 time constants
            Ts = kTot * wRun + tSlide;
         }
         if (floorApplied)
            SolveDrag(tSlide);
         else
         {
            b = 0.0;   // exact: the friction alone takes rampDownS
            coastFrictionS = rampDownS;
            frictionFaster = dragClamped = false;
         }
      }
   };

   // its own state beside the body's (PartRecords)
   struct DiscState
   {
      double volts = 0.0;       // the drive the script set; an input, not changed by the step (kept so a replay runs on the drive it ran on)
   };

   class Part : public GroupStep::BodyPart, public DiscState
   {
   public:
      Part() { KeepState<DiscState>(this); }
      Motor motor;
      V3 centre;                // centre of the top face, m
      double R = 0.0762;        // disc radius, m
      double zBottom = 0.0;     // bottom of the side when the disc stands proud of the playfield (= centre.z when it is flush: a disc with no side), m
      double edgeRhoM = 1e-3;   // the rim edge counts as rounded to this radius (Hertz effective radius, as ProbeTriangles' edgeRho)

      void Build(const Motor &m, const V3 &faceCentre, const double radius, const double bottom, const double edgeRho)
      {
         motor = m;
         centre = faceCentre;
         R = radius;
         zBottom = bottom;
         edgeRhoM = edgeRho;
         mass = motor.I;
         reach = R;
         // [Features] path_test: the disc's surface as ProbeSurface measures it (a solid cylinder: bottom, side, top; a flush disc: the top face), revolved facets, sweep radius 0. Its
         // geometry never moves (turning about its own axis moves no point of the surface), so the cover is built once here
         const double zTop = centre.z, zBot = std::min(zBottom, zTop);
         if (zTop > zBot)
         {
            const double r[4] = { 0.0, R, R, 0.0 }, z[4] = { zBot, zBot, zTop, zTop }, rr[4] = {};
            pathCover.Build(V3(centre.x, centre.y, 0.0), r, z, rr, 4);
         }
         else
         {
            const double r[2] = { R, 0.0 }, z[2] = { zTop, zTop }, rr[2] = {};
            pathCover.Build(V3(centre.x, centre.y, 0.0), r, z, rr, 2);
         }
      }
      RubberContact::PathRevolve pathCover;   // [Features] path_test (Build)
      // [Features] path_test: the ball's centre moving straight from s to s + d: 1 = it can touch the disc on the way (or an end lies inside the solid: ProbeSurface's depth there is above 0),
      // 0 = no (exact: the facets' distance minus their sagitta)
      int PathReach(const V3 &s, const V3 &d, const double Rb, const bool prof, bool &band) const
      {
         band = false;
         const double zTop = centre.z, zBot = std::min(zBottom, zTop);
         for (const V3 &p : { s, s + d })
         {
            const double dx = p.x - centre.x, dy = p.y - centre.y;
            if (zTop > zBot && dx * dx + dy * dy <= R * R && p.z <= zTop && p.z >= zBot)
               return 1;
         }
         return pathCover.Reach(s, d, Rb, prof, band);
      }
      double Angle() const   // 0 .. 2 pi (q runs on: a disc has no angle in its forces)
      {
         const double a = std::fmod(q, 2.0 * GroupStep::kPi);
         return a < 0.0 ? a + 2.0 * GroupStep::kPi : a;
      }

      // ---- the body's terms ----
      double StoredEnergy(double) const override { return 0.0; }
      double StoredChord(double, double, double, double &slope, double &lossF, double &lossSlope) const override
      {
         slope = lossF = lossSlope = 0.0;
         return 0.0;
      }
      double StoredForce(double) const override { return 0.0; }
      double Viscous() const override { return volts != 0.0 ? motor.kTot : motor.b; }
      double DriveForce() const override { return motor.Ts * volts * motor.dir / motor.ratedVolts; }
      double FrictionLimit(double, double, const V3 &) const override { return motor.frictionNm; }
      // the stick spring that carries the bearing and gear train's friction is the shaft's torsion, G J / L (user 2026-11-04: a 1/4 in steel shaft, 30 mm; no key)
      static constexpr double kShaftG = 79e9, kShaftR = 3.175e-3, kShaftL = 0.030;
      double PinStiffness(double, double, const V3 &) const override { return kShaftG * 0.5 * GroupStep::kPi * kShaftR * kShaftR * kShaftR * kShaftR / kShaftL; }

      // The nearest point of the disc (a solid cylinder from zBottom up to the top face; a flush disc is the top face alone) to a ball of radius Rb centred at c. On the rim edge the
      // contact counts as rounded to edgeRhoM as on any face edge (the effective radius falls with the angle from the nearer face, fully rounded at 15 deg). The column: the point's speed and the surface's spin per unit rate; no slope.
      void ProbeSurface(const V3 &c, const double Rb, Hit &h) const
      {
         const double zTop = centre.z, zBot = std::min(zBottom, zTop);
         const double dx = c.x - centre.x, dy = c.y - centre.y;
         const double r = std::sqrt(dx * dx + dy * dy);
         const V3 radial = r > 1e-12 ? V3(dx / r, dy / r, 0.0) : V3(1.0, 0.0, 0.0);
         const V3 up(0.0, 0.0, 1.0);
         V3 p, nFace, n;
         double dist;
         bool edge = false;
         if (r <= R && c.z <= zTop && c.z >= zBot && zTop > zBot)
         {
            // inside the solid: out through the nearer of the top face and the side
            const double toTop = zTop - c.z, toSide = R - r;
            if (toTop <= toSide)
            {
               p = V3(c.x, c.y, zTop);
               nFace = up;
               dist = -toTop;
            }
            else
            {
               p = V3(centre.x + radial.x * R, centre.y + radial.y * R, c.z);
               nFace = radial;
               dist = -toSide;
            }
            n = nFace;
         }
         else
         {
            const double rc = std::min(r, R), zc = std::min(std::max(c.z, zBot), zTop);
            p = V3(centre.x + radial.x * rc, centre.y + radial.y * rc, zc);
            const V3 d = c - p;
            dist = d.Length();
            if (r <= R)
               nFace = c.z >= zTop ? up : up * -1.0;   // over (or, for a flush disc, under) the face: a thin sheet, as the faces it replaces
            else if (c.z <= zTop && c.z >= zBot && zTop > zBot)
               nFace = radial;   // beside the side of a proud disc
            else
            {
               edge = true;   // past the rim and above (or below) the face: the edge
               nFace = std::abs(radial.Dot(d)) > std::abs(up.Dot(d)) && zTop > zBot ? radial : (c.z >= zTop ? up : up * -1.0);
            }
            n = (edge && dist > 1e-12) ? d / dist : nFace;
         }
         double edgeW = 0.0;   // 0 flat .. 1 fully rounded
         if (edge)
         {
            const double ang = std::acos(std::min(1.0, std::max(-1.0, n.Dot(nFace))));
            edgeW = std::min(1.0, std::max(0.0, (ang - 5.0 * GroupStep::kPi / 180.0) / (10.0 * GroupStep::kPi / 180.0)));
         }
         const double R1 = 1.0 / ((1.0 - edgeW) / Rb + edgeW * (1.0 / Rb + 1.0 / edgeRhoM));
         h.delta = Rb - dist;
         h.n = n;
         h.reff = std::sqrt(R1 * Rb);
         h.curvDist = edge ? dist : 0.0;   // rolling over the edge the ball turns round it
         h.grad = 1.0;
         h.nc = 1;
         h.c[0].Set(0, 0.0, V3(-(p.y - centre.y), p.x - centre.x, 0.0), up);
      }
   };

   // The disc's surface, carried by its part
   class Shape : public GroupStep::Shape
   {
   public:
      const Part *owner = nullptr;
      void Probe(const V3 &c, const double R, Hit &h) const override { owner->ProbeSurface(c, R, h); }
   };
}
