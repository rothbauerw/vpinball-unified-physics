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
         slopeSkirt = gp.cfg.skirt0.Slope();   // ([Features] part_finder: the profiles' largest slopes, for the path test)
         slopeRing = gp.cfg.ring0.Slope();
         origin = foot;
         m_reach = 0.5 * prm.skirtOdMm * MM + 2.0e-3;
         const double edgeMin = prm.edgeMinMm * MM;
         ringShape.prof = &gp.cfg.ring0;
         ringShape.edgeMin = edgeMin;
         ringShape.owner = &gp;
         bodyShape.prof = &gp.cfg.body0;
         bodyShape.edgeMin = edgeMin;
         // [Features] path_test: the fixed centre body's profile revolved into facets (world frame, once at load); the ring and skirt move with the part (build 3)
         std::vector<double> r, z, rr;
         for (const Revolve::Vtx &v : gp.cfg.body0.v)
         {
            r.push_back(v.r);
            z.push_back(v.z);
            rr.push_back(v.rr);
         }
         bodyCover.Build(origin, r.data(), z.data(), rr.data(), (int)r.size());
      }
      RubberContact::PathRevolve bodyCover;   // [Features] path_test (Build)
      double GroupPathCoverSag(const int surface, int &kind) const override   // (the self-check, VPX_FINDER_CHECK)
      {
         kind = 1;
         return surface == 2 && collidable ? bodyCover.MaxSag() : -1.0;
      }
      // [Features] path_test: surface 2, the centre body, never moves: its facets, reach R + the segment's larger sweep radius (RevolveShape's probe, owner null); not collidable: never touched
      int GroupPathReach(const int surface, const V3 &s, const V3 &d, const double R) const override
      {
         if (surface != 2)
            return -1;
         if (!collidable)
            return 0;
         const bool prof = ContactProf::CheckOn();
         bool band = false;
         const int r = bodyCover.Reach(s, d, R, prof, band);
         if (prof)
         {
            RubberContact::PathStats &ps = RubberContact::pathStats;
            ps.kTested[1]++;
            ps.kReached[1] += r > 0;
            ps.kBand[1] += band;
         }
         return r;
      }

      // [Features] part_finder, the path test's motion term (group A; GroupStep::Shape::PathMotion), P in the table's frame: 0 the skirt (BumperGroup::Part::SkirtMotion), 1 the ring:
      // it moves along the axis only, |d dr|; 2 the body is fixed (its own test, GroupPathReach)
      double slopeSkirt = 1.0, slopeRing = 1.0;   // Revolve::Slope of the skirt's and the ring's profiles (Build)
      double GroupPathSlope(const int surface) const override { return surface == 0 ? slopeSkirt : surface == 1 ? slopeRing : 1.0; }
      double GroupPathMotion(const int surface, const V3 &P, const double dLen, const double *const trial) const override
      {
         if (surface == 0)
            return gp.SkirtMotion(P - origin, dLen, trial);
         if (surface == 1)
            return trial != nullptr ? std::fabs(trial[BumperGroup::Part::kR]) : 0.0;
         return -1.0;
      }
      // the overlap check of free balls (part_finder, group A): |c - origin| less the bounding radius from its built profiles: the body; the ring plus its stroke (or its drop now, if more: the stop is a spring); the skirt plus its
      // pivot's height now (it turns about (0, 0, z))
      double RegionDist(const V3 &c) const override
      {
         const JetBumper::Bumper &b = gp.cfg;
         const double fr = std::max(BumperGroup::FarR(b.body0), std::max(BumperGroup::FarR(b.ring0) + std::max(b.p.ringStrokeMm * JetBumper::MM, std::fabs(gp.dr)), BumperGroup::FarR(b.skirt0) + std::fabs(gp.z)));
         return (c - origin).Length() - fr;
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

      // the energy reach (design_reach.md section 1, translating parts): the ring slides along the axis by at most EnergyReach(|vr|, E, ring mass, the coil at full current
      // plus the return spring plus gravity); the skirt's pivot slides by EnergyReach(|vz|, E, skirt mass, its springs plus gravity) and the skirt tilts by EnergyReach(|w|, E, I,
      // those springs' force at its farthest lever / I), a point c moving by its distance from the pivot times the tilt (at most pi) plus the slide (as SkirtMotion). The larger of the
      // two; each move at most what its energy allows (EnergyMove), not clamped at the stroke (the stops are springs). work: each body's kinetic energy now plus its drives' work over its move,
      // at most the energy they can supply
      void NextCycleStart(const double cycleS) override
      {
         if (GroupStep::jetPulseAtCycle)
            gp.CycleStart(cycleS);
      }
      double EnergyMove(const V3 &c, const double E, const double T, double *const work) const override
      {
         const Params &p = gp.cfg.p;
         const double g = partSet.g.Length();
         // the coil law (runs/2026-11-14-reach-p1b/brief_coil.md item 2): the current cannot jump, i <= min(i0 + V / L T, V / R), L its smallest (the ring up, the plunger farthest out);
         // the force at that current with the plunger closed (the largest gap factor). The coil's state for the cycle: [Features] jet_pulse_at_cycle on, a pulse starts only at a cycle's
         // start (BumperGroup::Part::CycleStart) and the script's command lands between cycles, so an off coil stays off this cycle (its current only decays: V 0); held mode follows the
         // switch inside the cycle, and with the key off a pulse fires inside the step whose switch closes: there an off coil is taken as switched on at the cycle's start
         const bool offStays = GroupStep::jetPulseAtCycle && p.mode != DriveMode::Held;
         const double V = gp.CoilOn() ? gp.CoilVolts() : offStays ? 0.0 : p.coil.vSupply;
         const double lMin = (gp.cfg.coilD.a + gp.cfg.coilD.b / (p.ringStrokeMm / p.coil.bobbinMm + p.plungerClosedFrac + p.coil.x0)) * 1e-3;
         const double iMax = std::min(std::fabs(gp.iCoil) + V / std::max(lMin, 1e-9) * T, std::max(V / p.coil.rOhm, std::fabs(gp.iCoil)));
         const double fCoil = CoilModel::Force(p.coil, gp.cfg.coilD, p.plungerClosedFrac, iMax);
         const double ringAccel = (fCoil + p.returnPreloadN + p.returnKNPerM * p.ringStrokeMm * MM) / p.ringMass + g;
         const double skirtF = p.skirtPreloadN + p.skirtKNPerM * p.skirtClearanceMm * MM, skirtAccel = skirtF / p.skirtMass + g;
         const double tiltAccel = (skirtF + p.skirtMass * g) * gp.reach / gp.TiltInertia(), wt = std::sqrt(gp.wx * gp.wx + gp.wy * gp.wy);
         // the law (a body's work over T is at most its kinetic energy plus what its drives can supply): the ring's drives: the coil's field energy now 1/2 L i^2 (L at its largest, the
         // plunger closed) plus the supply's V^2 / R T (coil on), the return spring's energy over the stroke (its largest force x the stroke) and gravity's m g over the move. The skirt's:
         // its springs' force at most skirtF over the clearance, gravity over its move. Each move also at most what a body with all that energy as speed makes in T, plus gravity's
         // 1/2 g T^2 (not clamped at the stroke or the clearance: past them the stops' own law only shortens it)
         const double Ep = std::max(E, 0.0), stroke = p.ringStrokeMm * MM, clear = p.skirtClearanceMm * MM, I = gp.TiltInertia();
         const double lMax = (gp.cfg.coilD.a + gp.cfg.coilD.b / (p.plungerClosedFrac + p.coil.x0)) * 1e-3;
         const double dRing = 0.5 * lMax * gp.iCoil * gp.iCoil + V * V / p.coil.rOhm * T + (p.returnPreloadN + p.returnKNPerM * stroke) * stroke;
         const double dSkirt = skirtF * clear;
         const double dr = DriveCappedReach(std::fabs(gp.vr), E, p.ringMass, ringAccel, dRing, T, 0.5 * g * T * T);
         const double dz = DriveCappedReach(std::fabs(gp.vz), E, p.skirtMass, skirtAccel, dSkirt, T, 0.5 * g * T * T);
         const double tiltG = p.skirtMass * g * gp.reach / I;
         const double dt = DriveCappedReach(wt, E, I, tiltAccel, dSkirt, T, 0.5 * tiltG * T * T);
         if (work != nullptr)
            *work = 0.5 * p.ringMass * gp.vr * gp.vr + std::min(p.ringMass * ringAccel * dr, dRing + p.ringMass * g * dr) + 0.5 * p.skirtMass * gp.vz * gp.vz + 0.5 * I * wt * wt
                  + std::min(p.skirtMass * skirtAccel * dz + I * tiltAccel * dt, dSkirt + p.skirtMass * g * dz + I * tiltG * dt);
         const V3 P = c - origin;
         const double rho = std::sqrt(P.x * P.x + P.y * P.y + (P.z - gp.z) * (P.z - gp.z));
         return std::max(dr, rho * std::min(dt, GroupStep::kPi) + dz);
      }

      // the nearest surface as the engine's probe: dist = the ball centre's distance from the surface, rho 0 (the surface is a sheet, its thickness is in the profile)
      RubberContact::TubeProbe Probe(const V3 &c) const override
      {
         RubberContact::TubeProbe pr;
         const V3 rel = c - origin;
         if (!collidable)
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
