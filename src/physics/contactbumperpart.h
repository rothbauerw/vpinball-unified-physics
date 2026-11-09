// license:GPLv3+
//
// The jet bumper as a moving part of the contact system, on the group step (runs/2026-11-02-bumper-group/notes.md). The skirt, the spoon as a spring on the follower distance, the ring and its coil are
// BumperGroup::Part (bumpergroup.h): the engine moves them with the balls in one solve, so an island with a bumper always goes to the group. There is no other path: no old-loop contact, no impulse, no stock.
//
// What this class adds to the group part: what the engine needs of any moving part (the probe for engagement and the overlap replay, the speed bound, the saved state), its three surfaces, its own steps when no
// ball is near, the script's commands, and its switch changes as events. Laws: ball contact is the shared law (GroupStep::Item) on each surface at its present pose; the reaction is in the same solve.
// Energy: only the coil adds energy; the contacts, the seats, the spoon's damper and the ring's stops remove it.
//
// World frame in metres, z up; the part is built at the foot of its axis (x, y, playfield height).

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>

#include "physics/contactmovingpart.h"
#include "physics/bumpergroup.h"

namespace JetBumper
{
   // what the engine side keeps (PartRecords: part of its state)
   struct EngineState
   {
      bool collidable = true;          // false: the bumper is not there at all (no contact, no switch, no coil, no events, parts at rest)
      bool idle = false;
      GroupStep::RestWatch rest;       // the rest rule's watch, alone and over its last cycle in a group (RestPose)
   };

   class EnginePart : public ContactMovingPart, public EngineState
   {
   public:
      EnginePart() { KeepState<EngineState>(this); }
      BumperGroup::Part gp;            // the skirt, the ring and the coil, as a part of the group
      GroupStep::Settings partSet;     // the gravity the part feels (SetGravity)
      BumperGroup::RevolveShape ringShape, bodyShape;
      V3 origin;                       // the foot of the bumper's axis on the playfield, world (m)
      RubberContact::BallSI ball;
      long failedSteps = 0, aloneSteps = 0, energySteps = 0;   // a step of the part alone could not be solved (the detector: rule 8); steps taken alone (cost); of them taken on energy alone (rule 8)
      int bookMs = 0;                  // the table's time, for the book of a cycle alone

      void Build(const Params &prm, const V3 &foot)
      {
         gp.Configure(prm, partSet);
         origin = foot;
         m_reach = 0.5 * prm.skirtOdMm * MM + 2.0e-3;
         const double edgeMin = prm.edgeMinMm * MM;
         ringShape.prof = &gp.cfg.ring0;
         ringShape.edgeMin = edgeMin;
         ringShape.owner = &gp;
         bodyShape.prof = &gp.cfg.body0;
         bodyShape.edgeMin = edgeMin;
      }

      // ---- the script ----
      // Collidable = False takes the bumper out of the world: parts return to rest, coil off, balls pass through, Solenoid ignored. True puts it back at rest.
      void SetCollidable(const bool c)
      {
         if (c == collidable)
            return;
         gp.Reset();
         collidable = c;
         idle = false;
      }
      // Bumper.Solenoid enabled[, volts[, duty]]: external mode sets the coil on/off; pulse mode fires one pulse on True; held mode ignores it (the caller warns once)
      bool Solenoid(const bool on, const double volts = -1.0, const double duty = 1.0)
      {
         if (!collidable || gp.cfg.p.mode == DriveMode::Held)
            return false;
         idle = false;
         if (gp.cfg.p.mode == DriveMode::External)
         {
            gp.extOn = on;
            gp.extVolts = on ? volts : -1.0;
            gp.extDuty = duty;
         }
         else if (on && !gp.pulseActive)
            gp.FirePulse();
         return true;
      }
      // The switch's changes so far as the part's events, in order: at a closure Hit first (no speed threshold), then the switch. The adapter calls it before the host takes the events.
      void Collect()
      {
         for (const char closed : gp.swLog)
         {
            if (closed)
            {
               events.push_back({ ContactMovingPart::Event::Hit, -1, 0.0 });
               events.push_back({ ContactMovingPart::Event::SwitchClosed, -1, 0.0 });
            }
            else
               events.push_back({ ContactMovingPart::Event::SwitchOpened, -1, 0.0 });
         }
         gp.swLog.clear();
      }

      // ---- the moving-part interface (contactmovingpart.h) ----
      void SetBall(const RubberContact::BallSI &b) override { ball = b; }
      void SetGravity(const V3 &g) override
      {
         const bool changed = (g - partSet.g).Length() > 1e-6;
         partSet.g = g;
         if (!gravityKnown)   // at rest on the table's own gravity from the start, not settling into it
         {
            gravityKnown = true;
            gp.Reset();
         }
         if (changed)
            idle = false;   // (a nudge wakes it)
      }
      V3 Origin() const override { return origin; }
      double Reach() const override { return m_reach; }
      // the ring: its speed now plus what the coil at full current in the closed position and the return spring can add in a cycle; the skirt: its speed (the rim moves faster than the centre)
      // plus what its springs can add. The bumper can start at any time (Solenoid), so even asleep the ring's acceleration counts.
      double SurfaceSpeedBound(double cycleS) const override
      {
         const Params &p = gp.cfg.p;
         const double ringAccel = (gp.cfg.coilD.strength + p.returnPreloadN + p.returnKNPerM * p.ringStrokeMm * MM) / p.ringMass;
         const double skirtAccel = (p.skirtPreloadN + p.skirtKNPerM * p.skirtClearanceMm * MM) / p.skirtMass;
         const double ring = std::abs(gp.vr) + ringAccel * cycleS;
         const double skirt = std::abs(gp.vz) + std::sqrt(gp.wx * gp.wx + gp.wy * gp.wy) * 0.5 * p.skirtOdMm * MM + skirtAccel * cycleS;
         return std::max(ring, skirt);
      }

      // the nearest surface as the engine's probe: dist = the ball centre's distance from the surface, rho 0 (the surface is a sheet, its thickness is in the profile)
      RubberContact::TubeProbe Probe(const V3 &c) const override
      {
         RubberContact::TubeProbe pr;
         const V3 rel = c - origin;
         if (!collidable || std::sqrt(rel.x * rel.x + rel.y * rel.y) > m_reach + ball.R + 0.05 || rel.z > 0.06 + ball.R)
         {
            pr.dist = INFINITY;
            pr.delta = -INFINITY;
            return pr;
         }
         GroupStep::Hit best;
         int surf = -1;
         for (int sf = 0; sf < kSurfaces; sf++)
         {
            GroupStep::Hit h;
            GroupProbe(sf, c, ball.R, nullptr, h);
            if (surf < 0 || h.delta > best.delta)
            {
               best = h;
               surf = sf;
            }
         }
         pr.seg = surf;
         pr.normal = best.n;
         pr.dist = ball.R - best.delta;
         pr.point = c - best.n * pr.dist;
         pr.rho = 0.0;
         pr.delta = best.delta;
         pr.overlap = best.delta > 0.0;
         pr.reff = best.reff;
         return pr;
      }

      // the surface nearest the point, and how it moves there
      void SurfaceVelocity(const V3 &point, V3 &v, V3 &w) const override
      {
         GroupStep::Hit best;
         double dBest = INFINITY;
         for (int sf = 0; sf < kSurfaces; sf++)
         {
            GroupStep::Hit h;
            GroupProbe(sf, point, 0.0, nullptr, h);
            if (std::fabs(h.delta) < dBest)
            {
               dBest = std::fabs(h.delta);
               best = h;
            }
         }
         v = V3();
         w = V3();
         for (int j = 0; j < best.nc; j++)
         {
            const double rate = gp.Coord(best.c[j].coord).rate;
            v = v + best.c[j].V() * rate;
            w = w + best.c[j].W() * rate;
         }
      }

      // The old loop's calls: an island with a bumper never goes there (GroupPart)
      RubberContact::ContactResult ContactItem(int, SurfaceContact &, const V3 &, const V3 &, const V3 &, double) override { assert(false); return RubberContact::ContactResult(); }
      double BallStep(const V3 &, const V3 &, double hMax) const override { assert(false); return hMax; }
      void StepInIsland(double) override { assert(false); }

      // ---- the group's calls ----
      GroupStep::Part *GroupPart() override { return &gp; }
      int GroupSurfaces() const override { return kSurfaces; }
      // 0 the skirt, 1 the ring, 2 the centre body (fixed: no columns). Not there: no depth, the same columns
      void GroupProbe(const int surface, const V3 &c, const double R, const double *const trial, GroupStep::Hit &h) const override
      {
         const V3 rel = c - origin;
         h.nc = 0;
         if (surface == 0)
            gp.ProbeSkirt(rel, R, trial, h);
         else if (surface == 1)
         {
            ringShape.trial = trial;
            ringShape.Probe(rel, R, h);
         }
         else
            bodyShape.Probe(rel, R, h);
         if (!collidable)
            h.delta = -1.0;
      }
      const RubberContact::RubberParams *GroupLaw(const int surface, double &layerM) const override
      {
         layerM = 0.0;
         return surface == 0 ? &gp.cfg.skirtMat : surface == 1 ? &gp.cfg.ringMat : &gp.cfg.bodyMat;
      }
      void CycleBegin() override { idle = false; }
      void TouchBegan(int, int, double) override {}   // (the bumper's Hit is its switch closing, not a touch)
      bool Sleeping() const override { return idle; }
      // The one rest rule (GroupStep::RestWatch) over its last cycle in a group: not driven (the coil's supply or a pulse), the skirt's and ring's move since the cycle's start and their speed now
      // within 1 nm, and the coil's current died away. A part out of the solve still moves alone, so a part caught at the top of a bounce only steps alone
      bool AtRest() const override
      {
         if (gp.CoilOn() || gp.pulseActive)
            return false;
         if (idle || !collidable)
            return true;
         return rest.Rests(RestPose(), RestSpeed());
      }
      GroupStep::RestWatch *Rest() override { return &rest; }
      // its one rest pose: the skirt (z, both tilts at the reach), the ring, and the coil's current (full scale: its steady current with the coil held on)
      GroupStep::RestWatch::Pose RestPose() const override
      {
         const double r = gp.reach;
         return GroupStep::RestWatch::Pose().Add(gp.z, 1.0).Add(gp.thx, r).Add(gp.thy, r).Add(gp.dr, 1.0).State(gp.iCoil, gp.cfg.coilD.iss);
      }

      // A bumper no ball is near: the same group with no balls, asleep while nothing is driven or moving, its state kept as it is (a part at rest sleeps, ModelStandards 3c)
      int AdvanceAlone(double dt) override
      {
         if (dt <= 0.0 || !collidable)
            return 0;
         if (idle && !gp.CoilOn())
            return 0;
         const bool was = idle;
         const int steps = GroupStep::AloneCycle(gp, dt, bookMs, failedSteps, aloneSteps, energySteps);
         // the one rest rule (GroupStep::RestWatch): driven while the coil's supply or a pulse is on; the skirt and the ring still for a whole cycle, the coil's current died away (RestPose)
         idle = rest.Step(RestPose(), RestSpeed(), gp.CoilOn() || gp.pulseActive, restTouched, dt);
         GroupStep::RestWatch::Book(0, was, idle);
         return steps;
      }

      std::string StateNote() const override
      {
         char buf[300];
         std::snprintf(buf, sizeof buf, "bumper: ring %.2f mm down at %.2f m/s, follower %.2f mm, skirt %.2f mm high tilt %.2f/%.2f deg, switch %d, coil %.2f A, asleep %d", gp.dr * 1e3, gp.vr, gp.FollowerD() * 1e3, gp.z * 1e3,
            gp.thx / DEG, gp.thy / DEG, (int)gp.closed, gp.iCoil, (int)idle);
         return buf;
      }

   private:
      static constexpr int kSurfaces = 3;
      double m_reach = 0.0;
      bool gravityKnown = false;
      double RestSpeed() const { return std::fabs(gp.vz) + (std::fabs(gp.wx) + std::fabs(gp.wy)) * gp.reach + std::fabs(gp.vr); }   // the skirt's and ring's surface speed bound, m/s
   };
}
