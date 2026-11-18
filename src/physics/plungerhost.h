// license:GPLv3+
//
// The mechanical plunger's rod (plungerrod.h) as a moving part of the contact system, and its settings from the unified vpx_physics.ini (docs/unified-physics/PlungerDesign.txt;
// engine step 1, runs/2026-11-17-plunger/e1-keys/brief.md): the probe for engagement, the speed bound and energy reach, its own steps when no ball is near, the hand's calls between
// cycles (keys and script) and its LimitBOS / LimitEOS events. Modelled on the Button cap's (buttoncap.h Cap). No VPX types.

#pragma once

#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

#include "physics/plungerrod.h"
#include "physics/contactmovingpart.h"
#include "physics/contactpartsini.h"

namespace PlungerHost
{
   using V3 = RubberContact::V3;

   // [Plungers] for all, then [Plunger:<name>]: every PlungerRod::Settings field by its name in snake case (the comment beside each field in plungerrod.h says what it is and its source).
   // A value that cannot build a working rod is refused: the key named in `bad` (empty = all good). rubber = 0 (the opt-out) is read by ContactHost::PlungerWanted
   inline void Load(PlungerRod::Settings &s, const std::filesystem::path &file, const std::string &name, std::vector<std::string> &bad)
   {
      PartsIni::ReadSections(file, "plungers", "plunger:" + name, [&](const std::string &k, double v, const std::string &) {
         if (k == "main_wire_mm") s.mainWireMm = v;                      // main spring wire diameter, mm
         else if (k == "main_od_mm") s.mainOdMm = v;                     // main spring outside diameter, mm
         else if (k == "main_active_coils") s.mainActiveCoils = v;       // main spring active coils
         else if (k == "main_closed_per_end") s.mainClosedPerEnd = v;    // main spring closed coils at each end
         else if (k == "main_free_mm") s.mainFreeMm = v;                 // main spring free length, mm
         else if (k == "gap_mm") s.gapMm = v;                            // the main spring's length at rest (housing to tip washer), mm
         else if (k == "shear_modulus_gpa") s.shearModulusGPa = v;       // spring wire shear modulus, GPa
         else if (k == "steel_density") s.steelDensity = v;              // spring wire density, kg/m3
         else if (k == "outer_rate_nmm") s.outerRateNmm = v;             // outer (barrel) spring rate, N/mm
         else if (k == "outer_free_mm") s.outerFreeMm = v;               // outer spring free length, mm
         else if (k == "outer_solid_mm") s.outerSolidMm = v;             // outer spring solid length, mm
         else if (k == "outer_wire_mm") s.outerWireMm = v;               // outer spring wire diameter, mm
         else if (k == "rod_g") s.rodG = v;                              // rod mass, g
         else if (k == "knob_g") s.knobG = v;                            // knob mass, g
         else if (k == "tip_g") s.tipG = v;                              // tip mass, g
         else if (k == "washer_clip_g") s.washerClipG = v;               // washer and clip mass, g
         else if (k == "stop_e") s.stopE = v;                            // coil-bind stops' restitution at stop_v_ref
         else if (k == "stop_v_ref") s.stopVRef = v;                     // the speed of stop_e, m/s
         else if (k == "drag_ns_per_m") s.dragNsPerM = v;                // viscous drag on the rod, N s/m
         else if (k == "bore_mu") s.boreMu = v;                          // Coulomb friction of the rod in its bore
         else if (k == "rod_radius_mm") s.rodRadiusMm = v;               // the rod's radius, mm
         else if (k == "tip_radius_mm") s.tipRadiusMm = v;               // the tip's hemisphere radius, mm
         else if (k == "tip_cyl_mm") s.tipCylMm = v;                     // the tip's cylinder behind the hemisphere, mm
         else if (k == "tip_rubber_mm") s.tipRubberMm = v;               // the tip's rubber thickness for the layer term, mm (0 = thick)
         else if (k == "tip_e_mpa") s.tipEMpa = v;                       // the tip's rubber modulus, MPa
         else if (k == "tip_alpha") s.tipAlpha = v;                      // the tip's rubber damping, s/m
         else if (k == "tip_beta_layer") s.tipBetaLayer = v;             // the tip's layer term factor
         else if (k == "tip_mu") s.tipMu = v;                            // the tip's friction on the ball
         else if (k == "tip_shear_scale") s.tipShearScale = v;           // the tip's sideways stiffness factor
         else if (k == "tip_tan_delta") s.tipTanDelta = v;               // the tip's loss tangent
         else if (k == "pull_mm_s") s.pullMmS = v;                       // the hand's pull speed on a key, mm/s
         else if (ContactHost::IsSurfaceKey(k)) {}                       // (rubber = 0: the opt-out, ContactHost::PlungerWanted)
      });
      const auto need = [&](const bool ok, const char *key) {
         if (!ok)
            bad.emplace_back(key);
      };
      need(s.mainWireMm > 0.0, "main_wire_mm");
      need(s.mainOdMm > s.mainWireMm, "main_od_mm");
      need(s.mainActiveCoils > 0.0, "main_active_coils");
      need(s.mainClosedPerEnd >= 0.0, "main_closed_per_end");
      need(s.mainFreeMm > s.gapMm, "main_free_mm");   // a preload: the spring is compressed at rest
      need(s.gapMm > (s.mainActiveCoils + 2.0 * s.mainClosedPerEnd) * s.mainWireMm, "gap_mm");   // a travel: the gap exceeds the spring's solid length
      need(s.shearModulusGPa > 0.0, "shear_modulus_gpa");
      need(s.steelDensity > 0.0, "steel_density");
      need(s.outerRateNmm > 0.0, "outer_rate_nmm");
      need(s.outerSolidMm > 0.0, "outer_solid_mm");
      need(s.outerWireMm > 0.0, "outer_wire_mm");
      need(s.rodG >= 0.0 && s.knobG >= 0.0 && s.tipG >= 0.0 && s.washerClipG >= 0.0 && s.rodG + s.knobG + s.tipG + s.washerClipG > 0.0, "rod_g");
      need(s.stopE > 0.0 && s.stopE <= 1.0, "stop_e");
      need(s.stopVRef > 0.0, "stop_v_ref");
      need(s.dragNsPerM >= 0.0, "drag_ns_per_m");
      need(s.boreMu >= 0.0, "bore_mu");
      need(s.rodRadiusMm > 0.0, "rod_radius_mm");
      need(s.tipRadiusMm > 0.0, "tip_radius_mm");
      need(s.tipCylMm > 0.0, "tip_cyl_mm");
      need(s.tipRubberMm >= 0.0, "tip_rubber_mm");
      need(s.tipEMpa > 0.0, "tip_e_mpa");
      need(s.tipAlpha >= 0.0, "tip_alpha");
      need(s.tipMu >= 0.0, "tip_mu");
      need(s.pullMmS > 0.0, "pull_mm_s");
      if (bad.empty())   // the outer spring must reach its solid length forward of rest: its free length less solid exceeds its rest squeeze c0 = preload / k_o
      {
         const double d = s.mainWireMm, D = s.mainOdMm - d, km = s.shearModulusGPa * 1e3 * d * d * d * d / (8.0 * D * D * D * s.mainActiveCoils);   // N/mm
         need(s.outerFreeMm - s.outerSolidMm > km * (s.mainFreeMm - s.gapMm) / s.outerRateNmm, "outer_free_mm");
      }
   }

   // the engine side (PartRecords: part of its state)
   struct EngineState
   {
      bool idle = false;
      GroupStep::RestWatch rest;   // the rest rule's watch, alone and over its last cycle in a group (RestPose)
   };

   class Rod : public ContactMovingPart, public EngineState
   {
   public:
      Rod() { KeepState<EngineState>(this); }
      PlungerRod::Part part;
      RubberContact::BallSI ball;
      long failedSteps = 0, aloneSteps = 0, energySteps = 0;   // the part alone: a step not solved (rule 8), steps taken (cost), of them taken on energy alone (rule 8)
      int bookMs = 0;
      double slope = 1.0;
      V3 mid;                // the middle of the tip's travel (the host's box, Origin)
      bool retract = false;  // PullBackandRetract: ease back to rest once at full pull (design 2g)
      bool limitArmed = false;
      double lastQ = 0.0;

      void Build(const PlungerRod::Settings &s, const V3 &axis, const V3 &gravity, const V3 &apex, const RubberContact::BallSI &b)
      {
         ball = b;
         part.Build(s, axis, gravity, apex);
         slope = part.prof.Slope();
         mid = apex + axis * (0.5 * (part.travel + part.outerSolidS));
         idle = true;   // at rest on its springs' balance: nothing to step until a ball comes near, the hand moves or the gravity changes
      }
      // the largest force its springs and weight make anywhere in its travel (N)
      double MaxStoredForce() const { return part.km * (part.x0 + part.travel) + part.ko * (part.c0 - part.outerSolidS) + part.mass * part.g.Length(); }

      // ---- the hand (between cycles; keys and script, design 2b, 2g) ----
      void KeyDown(const bool andRetract)
      {
         part.KeyDown(part.pullRate);
         retract = andRetract;
         hand = Hand::Script;    // 6f: the script's (key's) hand holds; the sensor is ignored while it does
         sensorLatched = true;   // ... and after it, until the rod rests and the sensor reads rest (the 2d re-arm)
         Wake();
      }
      void KeyUp()
      {
         part.KeyUp();
         retract = false;
         if (hand == Hand::Script)
            hand = Hand::None;   // (the latch stays: re-armed by SensorStep)
         Wake();
      }

      // ---- the cabinet sensor's hand (engine step 2A; design 2c, 2d, 6f; runs/2026-11-17-plunger/e2a-sensor/brief.md) ----
      enum class Hand { None, Sensor, Script };
      Hand hand = Hand::None;       // who holds the hold stop: nobody, the sensor's hand, or the script's PullBack (6f)
      bool sensorLatched = false;   // 2d: after a release (or a script pull) the sensor is ignored until the rod rests and the sensor reads rest
      double sensorPrev = 0.0;      // the last reading's s, m (the first reading is taken from rest)
      double sensorSince = 0.0;     // s since the reading last changed (a USB sensor's sample is held about 8 ms between reports)
      static constexpr double kHomeM = 1e-5;         // m (0.01 mm, lead 2026-11-18, fix F2): the rod and the hold stop are "at rest position 0" within this when the hand lets go
      static constexpr double kReleaseSpeed = 0.5;   // m/s, design 2d (user 2026-11-17): above a hand-guided return (about 0.3), below a released knob (about 0.7 from 5 mm on the weakest spring)
      enum SensorEvent { kSensorNone, kSensorTake, kSensorReleaseSpeed, kSensorReleaseAhead, kSensorLetGoRest, kSensorRearm, kSensorLost };
      // Once per physics millisecond, after the cycle: VPX's filtered sensor state. present: a cabinet plunger (MechPlunger) with a sensor (H7); posPu: p.u. (0 rest, 1 full pull,
      // < 0 forward); velPuS: p.u./s (a velocity sensor's own when mapped); restBandPu: the rest band of the re-arm, p.u.; dt: the millisecond, s. Returns what the hand did (the log's)
      int SensorStep(const bool present, const double posPu, const double velPuS, const double restBandPu, const double dt)
      {
         PlungerRod::Part &p = part;
         const double s = std::max(posPu, 0.0) * p.travel;   // 2c: calibrated rest -> 0, full pull -> travel; a reading forward of rest is rest in this piece (palm smash: piece 2B)
         const double prev = sensorPrev, interval = sensorSince + dt;   // the previous reading, and the time since it first read so (one sensor sample)
         const bool changed = s != prev;
         sensorPrev = s;
         sensorSince = changed ? 0.0 : interval;
         if (!present)   // no sensor (or the override turned off): a hand it held lets go, the next reading starts from rest
         {
            sensorPrev = sensorSince = 0.0;
            if (hand != Hand::Sensor)
               return kSensorNone;
            p.KeyUp();
            hand = Hand::None;
            Wake();
            return kSensorLost;
         }
         if (hand == Hand::Script)   // 6f
            return kSensorNone;
         if (sensorLatched)   // 2d re-arm: the rod at rest (the engine's rest test, the one that lets it sleep) and the sensor reading rest within the band
         {
            if (!AtRest() || std::fabs(posPu) > restBandPu)
               return kSensorNone;
            sensorLatched = false;
            return kSensorRearm;
         }
         if (hand == Hand::Sensor)
         {
            // 2d release, only from a pulled reading (the hand held the knob behind rest: the previous reading > 0): the knob toward rest faster than kReleaseSpeed (the
            // sensor's speed x travel), or the reading ahead of the model rod (the rod has left the hand): the gap between the rod and the new reading, over the time since
            // the previous reading came, faster than the same speed (a hand-guided return keeps the rod within one sample's move of the reading)
            if (prev > 0.0)
            {
               const bool fast = -velPuS * p.travel > kReleaseSpeed;
               const bool ahead = changed && (p.q - s) / interval > kReleaseSpeed;
               if (fast || ahead)
               {
                  p.KeyUp();   // the hold stop is removed: the spring alone drives the rod
                  hand = Hand::None;
                  sensorLatched = true;
                  Wake();
                  return fast ? kSensorReleaseSpeed : kSensorReleaseAhead;
               }
            }
            const GroupStep::BodyStop &h = p.stops[p.stopHold];
            if (changed)   // 2c (fix F1): the hold stop toward the new reading at the speed between the last two distinct samples, the step over the time between them (a held
            {              // 125 Hz sample moves it over its 8 ms, not in a 1 ms burst; through VPX's filter each 1 ms reading is a sample, so the same rule serves both)
               p.DriveStop(p.stopHold, s, std::fabs(s - prev) / interval);
               Wake();
            }
            // fix F2: the hand guides the knob home: a reading at rest drives the hold stop to 0 (above), and the hand lets go only when the rod rests on it there: rod and stop at 0
            // within kHomeM, the rod's speed within the rest test's (RestWatch: kMoveM per cycle). Not when the reading first reads rest: the rod would fly on at the hand's speed
            else if (s == 0.0 && h.target == 0.0 && std::fabs(h.at) <= kHomeM && std::fabs(p.q) <= kHomeM
                     && std::fabs(p.rate) <= GroupStep::RestWatch::kMoveM / GroupStep::RestWatch::kCycleS)
            {
               p.KeyUp();
               hand = Hand::None;
               Wake();
               return kSensorLetGoRest;
            }
            return kSensorNone;
         }
         // the hand off: a pulled reading behind the rod takes the knob (the hold stop placed at the rod, then driven to the reading); a reading at rest leaves it off
         if (s > 0.0 && s > p.q)
         {
            p.SetStop(p.stopHold, p.q, true);
            p.DriveStop(p.stopHold, s, std::max(std::fabs(s - prev), s - p.q) / interval);   // (F1: the speed between the last two distinct samples; at least the rod's gap over that time)
            hand = Hand::Sensor;
            Wake();
            return kSensorTake;
         }
         return kSensorNone;
      }
      // the hand changed a stop: step it from the next cycle. Its CycleBegin (the hold stop's touch, an impact booked as the hand's drive) is the engine's, inside that cycle's book
      void Wake() { idle = false; }

      // After the cycle (the adapter's Publish): the stroke events, as stock arms them (a 1/50 travel band at each end, re-armed past twice that), and PullBackandRetract's return
      void Collect()
      {
         part.stopEvents.clear();
         const double q = part.q, band = part.travel / 50.0;
         if (limitArmed && q < band && lastQ >= band)   // reached rest going forward
         {
            events.push_back({ Event::LimitBOS, -1, std::fabs(part.rate) });
            limitArmed = false;
         }
         else if (limitArmed && q > part.travel - band)   // reached the full-pull stop
         {
            events.push_back({ Event::LimitEOS, -1, std::fabs(part.rate) });
            limitArmed = false;
         }
         else if (q > 2.0 * band && q < part.travel - 2.0 * band)
            limitArmed = true;
         lastQ = q;
         if (retract && q >= part.travel - 1e-4)   // at full pull (0.1 mm): the hand eases it back to rest at the pull speed
         {
            retract = false;
            part.Retract(part.pullRate);
            Wake();
         }
      }

      // ---- the moving-part interface (contactmovingpart.h) ----
      void SetBall(const RubberContact::BallSI &b) override { ball = b; }
      void SetGravity(const V3 &g) override
      {
         if ((g - part.g).Length() > 1e-6)
            idle = false;   // (a nudge wakes it)
         part.g = g;
      }
      V3 Origin() const override { return mid; }
      double Reach() const override { return BumperGroup::FarR(part.prof) + part.travel - part.outerSolidS; }   // every point of the tip at any pull
      double MotionRange() const override { return part.travel - part.outerSolidS; }
      double SurfaceSpeedBound(const double cycleS) const override { return std::fabs(part.rate) + part.pullRate + MaxStoredForce() / part.mass * cycleS; }
      // the energy reach (translating parts, as the Button cap): its springs and weight as the drive over its whole travel, the hand's slow pull inside the speed bound
      double EnergyMove(const V3 &, const double E, const double T, double *const work) const override
      {
         const double F = MaxStoredForce(), D = F * (part.travel - part.outerSolidS), v0 = std::fabs(part.rate) + part.pullRate;
         const double x = DriveCappedReach(v0, E, part.mass, F / part.mass, D, T);
         if (work != nullptr)
            *work = 0.5 * part.mass * v0 * v0 + std::min(F * x, D);
         return x;
      }
      RubberContact::TubeProbe Probe(const V3 &c) const override
      {
         GroupStep::Hit h;
         part.ProbeSurface(c, ball.R, 0.0, h);
         RubberContact::TubeProbe pr;
         pr.seg = 0;
         pr.normal = h.n;
         pr.dist = ball.R - h.delta;
         pr.point = c - h.n * pr.dist;
         pr.rho = 0.0;
         pr.delta = h.delta;
         pr.overlap = h.delta > 0.0;
         pr.reff = h.reff;
         return pr;
      }
      void SurfaceVelocity(const V3 &, V3 &v, V3 &w) const override
      {
         v = part.axis * part.rate;
         w = V3();
      }

      // ---- the group's calls ----
      GroupStep::Part *GroupPart() override { return &part; }
      int GroupSurfaces() const override { return 1; }
      void GroupProbe(int, const V3 &c, const double R, const double *const trial, GroupStep::Hit &h) const override { part.ProbeSurface(c, R, trial != nullptr ? trial[0] : 0.0, h); }
      double GroupPathMotion(int, const V3 &, double, const double *const trial) const override { return trial != nullptr ? std::fabs(trial[0]) : 0.0; }   // a slide: |dq|
      double RegionDist(const V3 &c) const override { return part.RegionDist(c); }
      double GroupPathSlope(int) const override { return slope; }
      const RubberContact::RubberParams *GroupLaw(int, double &layerM) const override
      {
         layerM = part.layerM;
         return &part.rp;
      }
      void CycleBegin() override
      {
         idle = false;
         part.CycleBegin();
      }
      bool Sleeping() const override { return idle; }
      bool RestNumbers(double &moveM, double &speedMS) const override
      {
         moveM = rest.CycleMoved(RestPose());
         speedMS = std::fabs(part.rate);
         return true;
      }
      // the one rest rule; never at rest while the hand moves a stop toward its target
      bool AtRest() const override { return !part.AnyStopMoving() && (idle || rest.Rests(RestPose(), std::fabs(part.rate))); }
      GroupStep::RestWatch *Rest() override { return &rest; }
      GroupStep::RestWatch::Pose RestPose() const override { return GroupStep::RestWatch::Pose().Add(part.q, 1.0); }
      // No ball near: the same group with no balls, asleep while nothing moves and the hand is still
      int AdvanceAlone(const double dt) override
      {
         if (dt <= 0.0 || idle)
            return 0;
         const bool was = idle;
         part.ClearLoads();
         part.CycleBegin();
         const int steps = GroupStep::AloneCycle(part, dt, bookMs, failedSteps, aloneSteps, energySteps);
         idle = rest.Step(RestPose(), std::fabs(part.rate), part.AnyStopMoving(), restTouched, dt);
         GroupStep::RestWatch::Book(1, was, idle);
         return steps;
      }

      std::string StateNote() const override
      {
         char buf[200];
         std::snprintf(buf, sizeof buf, "plunger: pull %.4f mm at %.4f m/s, on stop %d, asleep %d", part.q * 1e3, part.rate, part.PinStop(), (int)idle);
         return buf;
      }
   };
}
