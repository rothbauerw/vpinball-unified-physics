// license:GPLv3+
//
// The mechanical plunger's rod (plungerrod.h) as a moving part of the contact system, and its settings from the unified vpx_physics.ini (docs/unified-physics/PlungerDesign.txt;
// engine step 1, runs/2026-11-17-plunger/e1-keys/brief.md): the probe for engagement, the speed bound and energy reach, its own steps when no ball is near, the hand's calls between
// cycles (keys and script) and its LimitBOS / LimitEOS events. Modelled on the Button cap's (buttoncap.h Cap). No VPX types.

#pragma once

#include <algorithm>
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
         else if (k == "rub_n") s.rubN = v;                              // the rod and main spring rubbing in the sleeve: a constant Coulomb force, N
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
         else if (k == "release_accel") s.releaseAccel = v;              // the cabinet sensor: release when the fitted acceleration toward rest exceeds this, m/s^2 (2k b)
         else if (k == "release_accel_min_speed") s.releaseAccelMinSpeed = v;   // ... with the fitted speed toward rest over this, m/s (2k b)
         else if (k == "release_window_c") s.releaseWindowC = v;         // the fit's window constant C: W = sqrt(2 C q / release_window_noise) (2k a)
         else if (k == "release_window_noise") s.releaseWindowNoise = v; // the fit window's noise allowance, m/s^2, independent of release_accel (2k a, 14.3)
         else if (k == "release_window_min_ms") s.releaseWindowMinMs = v; // the fit window's floor, ms: a hand's own short-term acceleration (2k a, 14.3)
         else if (k == "sensor_filtered") s.sensorFiltered = v != 0.0;   // the cabinet sensor: 0 = its own reports (default), 1 = VPX's filtered position and speed (2i a)
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
      need(s.rubN >= 0.0, "rub_n");
      need(s.rodRadiusMm > 0.0, "rod_radius_mm");
      need(s.tipRadiusMm > 0.0, "tip_radius_mm");
      need(s.tipCylMm > 0.0, "tip_cyl_mm");
      need(s.tipRubberMm >= 0.0, "tip_rubber_mm");
      need(s.tipEMpa > 0.0, "tip_e_mpa");
      need(s.tipAlpha >= 0.0, "tip_alpha");
      need(s.tipMu >= 0.0, "tip_mu");
      need(s.pullMmS > 0.0, "pull_mm_s");
      need(s.releaseAccel > 0.0, "release_accel");
      need(s.releaseAccelMinSpeed >= 0.0, "release_accel_min_speed");
      need(s.releaseWindowC > 0.0, "release_window_c");
      need(s.releaseWindowNoise > 0.0, "release_window_noise");
      need(s.releaseWindowMinMs > 0.0, "release_window_min_ms");
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
         releaseAccel = s.releaseAccel;
         releaseAccelMinSpeed = s.releaseAccelMinSpeed;
         releaseWindowC = s.releaseWindowC;
         releaseWindowNoise = s.releaseWindowNoise;
         releaseWindowMin = s.releaseWindowMinMs * 1e-3;
         sensorFiltered = s.sensorFiltered;
         slope = part.prof.Slope();
         mid = apex + axis * (0.5 * (part.travel + part.outerSolidS));
         idle = true;   // at rest on its springs' balance: nothing to step until a ball comes near, the hand moves or the gravity changes
      }
      // the largest force its springs and weight make anywhere in its travel (N)
      double MaxStoredForce() const { return part.km * (part.x0 + part.travel) + part.ko * (part.c0 - part.outerSolidS) + part.mass * part.g.Length(); }

      // ---- the hand (between cycles; keys and script, design 2b, 2g) ----
      void KeyDown(const bool andRetract)
      {
         if (hand == Hand::Sensor)   // 6f: the script has priority: the sensor's tie lets go before the hold stop takes the knob
            part.Untie();
         part.KeyDown(part.pullRate);
         retract = andRetract;
         hand = Hand::Script;    // 6f: the script's (key's) hand holds; the sensor is ignored while it does
         resync = false;
         sensorLatched = true;   // ... and after it, until the rod rests (U2: then the tie resyncs to the reading)
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

      // ---- the cabinet sensor's hand (engine step 2A; design 2c, 2d, 6f; user U1, U2 2026-11-18, runs/2026-11-17-plunger/e2c-tie/brief.md) ----
      enum class Hand { None, Sensor, Script };
      Hand hand = Hand::None;       // who holds the knob: nobody, the sensor's tie (hold and push stops together at the reading, U1), or the script's PullBack (6f)
      bool sensorLatched = false;   // 2d: after a release (or a script pull) the sensor is ignored until the rod rests (U2)
      bool resync = false;          // U2: the tie is moving the rod from where it rested to the reading at the pull speed (until the stops reach the reading)
      double sensorPrev = 0.0;      // the s of the last distinct raw reading, m (the first reading is taken from rest)
      double sensorRawPrev = 0.0;   // the last distinct raw reading, p.u. (G2: a sample is a change of the raw reading, not of s: noise within the dead zone is samples at s = 0)
      double sensorSince = 0.0;     // s since the raw reading last changed (a USB sensor's sample is held about 8 ms between reports)
      double releaseAccel = 15.0, releaseAccelMinSpeed = 0.15, releaseWindowC = 10.6, releaseWindowNoise = 5.0, releaseWindowMin = 0.022;   // 2k: the Settings (min in s)
      // 2k a: the pulled reports (s > 0) since the tie took the knob, (time s, s m), oldest first, at most kFitMax; the board's step q (m, the smallest change seen between
      // two pulled reports, 0 until seen); the fit's result at the last report (speed and acceleration toward rest, the reports and span used; fitN 0 = no fit)
      static constexpr int kFitMax = 32;
      double fitT[kFitMax] = {}, fitS[kFitMax] = {};
      int fitCount = 0;
      double sensorClock = 0.0, boardStep = 0.0;
      double fitSpeed = 0.0, fitAccel = 0.0, fitSpan = 0.0;
      int fitN = 0;
      // e2m a: the board's report period T_b (s, 0 until measured): the median of the last kGapRing gaps between distinct raw reports (gaps over kSampleMaxS left out),
      // once kGapMin are seen. 15 gaps: about 45 ms of a 3 ms board, enough that one jittered gap does not move the median
      static constexpr int kGapRing = 15, kGapMin = 3;
      double gapS[kGapRing] = {};
      int gapCount = 0, gapNext = 0;
      double reportPeriod = 0.0;
      bool sensorFiltered = false;  // 2i a: Settings::sensorFiltered (the adapter reads VPX's filtered state instead of the sensor's own reports)
      static constexpr double kRestBandPu = 0.01;   // p.u. (D1, 2i c): the hand's rest band on the mapped reading, a noise floor only (stock's own, PlungerHandler.cpp:276);
                                                     // VPX's mapping has already removed its dead zone (SensorMapping.h:188-195), so it is not applied again here
      static constexpr double kSampleMaxS = 0.010;   // s (G2, lead 2026-11-18): VPX's own sensor timeout (PlungerHandler.cpp:290 repeats a sample older than 10 ms): a reading
                                                      // held longer is a hand at rest, not a slow move, so the time between two samples is at most this
      static constexpr double kReleaseSpeed = 0.5;   // m/s, design 2d (user 2026-11-17): above a hand-guided return (about 0.3), below a released knob (about 0.7 from 5 mm on the weakest spring)
      enum SensorEvent { kSensorNone, kSensorTake, kSensorReleaseSpeed, kSensorLetGoRest, kSensorRearm, kSensorLost, kSensorResync, kSensorReleaseAccel };
      // Once per physics millisecond, after the cycle: the sensor's state (2i a: its own reports, or VPX's filtered state with sensor_filtered). present: a cabinet plunger
      // (MechPlunger) with a sensor (H7); posPu: p.u. (0 rest, 1 full pull, < 0 forward); velPuS: p.u./s, used only when velOwn (a velocity sensor's own when mapped, or the
      // filtered speed); else the speed for the release test is the step between distinct reports over the time between them (2i a); dt: the millisecond, s. Returns what
      // the hand did (the log's)
      int SensorStep(const bool present, const double posPu, const double velPuS, const bool velOwn, const double dt)
      {
         PlungerRod::Part &p = part;
         // 2c: calibrated rest -> 0, full pull -> travel. D1 (6g): a reading within the rest band (kRestBandPu, a noise floor on the mapped reading: VPX's mapping already
         // removed its dead zone) is rest, s = 0 on both sides: a noisy sensor at rest neither ties the knob nor moves it. G1: outside the band the reading is measured from
         // the band's edge and rescaled, s = sign(p) (|p| - band) / (1 - band) x travel, so s is continuous at the edge and a full pull 1.0 still maps to the full travel.
         // D2 (2e): a reading forward of rest is the palm pushing the knob in, s < 0 on the same scale, no further than the outer spring's solid point (the filter reads to
         // -1/3 p.u. fully extended, past it)
         const double over = std::max(std::fabs(posPu) - kRestBandPu, 0.0) / (1.0 - kRestBandPu) * p.travel;
         const double s = std::max(posPu < 0.0 ? -over : over, p.outerSolidS);
         // G2: a sample is a distinct raw reading. The speed between samples is the step in s since the previous distinct raw reading over the time since it, that time at
         // most kSampleMaxS (VPX's sensor timeout: a reading held longer is a hand at rest, not a slow move)
         // e2m a (design 2c, the tie at the knob's speed between reports): after a still spell (no raw change for over 1.5 report periods: at least one report missed,
         // beyond the board's jitter) a step of more than one board step was made within the board's last report period, so it is timed over that period T_b
         // (reportPeriod, measured below); a one-step change (a slow move of one count per spell) keeps the time since the last change. Until T_b and q are measured: today's timing
         const bool rawChanged = posPu != sensorRawPrev;
         const double still = sensorSince + dt;   // s since the previous distinct raw reading
         const bool afterStill = rawChanged && reportPeriod > 0.0 && still > 1.5 * reportPeriod;
         const double prev = sensorPrev;
         const double interval = afterStill && boardStep > 0.0 && std::fabs(s - prev) > 1.5 * boardStep ? reportPeriod : std::min(still, kSampleMaxS);
         const bool changed = s != prev;   // the stops move only on a step in s (a raw change within the dead zone is a sample at rest, s stays 0)
         sensorClock += dt;
         fitN = 0;
         if (rawChanged)
         {
            sensorRawPrev = posPu;
            sensorPrev = s;
            sensorSince = 0.0;
            // e2m a: the board's report period T_b, the median of the last kGapRing gaps between distinct reports no longer than kSampleMaxS (a longer gap is a still
            // spell, not a report period): the typical gap, not the shortest (a 2 ms report among 3 ms ones would time steps 1.5 x fast). 0 until kGapMin gaps are seen
            if (still <= kSampleMaxS)
            {
               gapS[gapNext] = still;
               gapNext = (gapNext + 1) % kGapRing;
               gapCount = std::min(gapCount + 1, kGapRing);
               if (gapCount >= kGapMin)
               {
                  double g[kGapRing];
                  std::copy(gapS, gapS + gapCount, g);
                  std::nth_element(g, g + gapCount / 2, g + gapCount);
                  reportPeriod = g[gapCount / 2];
               }
            }
            // e2m c: a still spell longer than kSampleMaxS (VPX's sensor timeout: a hand at rest) ends the fit's motion: the fit restarts from this report (the reports
            // before the spell are a held knob, not this move: e2l, a 300 ms hold fitted with the first move's reports)
            if (still > kSampleMaxS)
               fitCount = 0;
            // 2k a: a pulled report joins the fit's list (a reading at or forward of rest empties it); the board's step q is the smallest change between two pulled reports
            if (s > 0.0)
            {
               if (fitCount > 0 && s != fitS[fitCount - 1])
               {
                  const double d = std::fabs(s - fitS[fitCount - 1]);
                  boardStep = boardStep > 0.0 ? std::min(boardStep, d) : d;
               }
               if (fitCount == kFitMax)
               {
                  for (int i = 1; i < kFitMax; i++)
                  {
                     fitT[i - 1] = fitT[i];
                     fitS[i - 1] = fitS[i];
                  }
                  fitCount--;
               }
               fitT[fitCount] = sensorClock;
               fitS[fitCount++] = s;
            }
            else
               fitCount = 0;
         }
         else
            sensorSince += dt;
         // the knob's speed toward rest for the release test, m/s: the sensor's own (velOwn), else the step in s between distinct reports over the time between them (0 between
         // reports: a release is seen on the report that shows it)
         const double inSpeed = velOwn ? -velPuS * p.travel : (changed ? (prev - s) / interval : 0.0);
         if (!present)   // no sensor (or the override turned off): a tie lets go, the next reading starts from rest
         {
            sensorPrev = sensorRawPrev = sensorSince = 0.0;
            resync = false;
            fitCount = 0;
            if (hand != Hand::Sensor)
               return kSensorNone;
            p.Untie();
            hand = Hand::None;
            Wake();
            return kSensorLost;
         }
         if (hand == Hand::Script)   // 6f (a tie let go in KeyDown)
            return kSensorNone;
         // U2: after a release (or a script pull) the sensor is ignored, a smash or pull included, until the model rod rests (the engine's rest test, the one that lets it
         // sleep), wherever the sensor reads. Then a reading outside the rest band ties the rod where it rests and moves it to the reading at the key's pull speed (no jump);
         // a reading at rest leaves it to the springs
         if (sensorLatched)
         {
            if (!AtRest())
               return kSensorNone;
            sensorLatched = false;
            if (s == 0.0)
               return kSensorRearm;
            FitRestart(s);
            p.Tie(p.q);
            p.DriveTie(s, p.pullRate);
            hand = Hand::Sensor;
            resync = true;
            Wake();
            return kSensorResync;
         }
         if (hand == Hand::Sensor)
         {
            // 2k MERGED (TODO 14.3, user 2026-11-18), only from a pulled reading (the previous reading > 0): one rule, the knob speeding up toward rest (a free knob: its
            // spring keeps it speeding up into rest; a hand guiding it back slows it), fitted acceleration toward rest over release_accel at a fitted speed over
            // release_accel_min_speed, or the knob toward rest faster than kReleaseSpeed. A quadratic fit of the pulled reports gives the speed and acceleration toward
            // rest; without a fit (under 4 reports, a reading at or forward of rest, or velOwn) the speed test takes the step between reports (the sensor's own speed with
            // velOwn): a board too coarse to fit keeps the speed test only (2k a)
            if (prev > 0.0)
            {
               const bool fit = changed && s > 0.0 && !velOwn && Fit();
               const bool fast = (fit ? fitSpeed : inSpeed) > kReleaseSpeed;
               const bool accel = fit && fitAccel > releaseAccel && fitSpeed > releaseAccelMinSpeed;
               if (fast || accel)
               {
                  p.Untie();   // the hand lets go: the spring alone drives the rod
                  hand = Hand::None;
                  resync = false;
                  sensorLatched = true;
                  Wake();
                  return fast ? kSensorReleaseSpeed : kSensorReleaseAccel;
               }
            }
            const GroupStep::BodyStop &h = p.stops[p.stopHold];
            if (!p.AnyStopMoving())   // U2: the stops have reached the reading: the tie follows normally from here
               resync = false;
            if (changed)   // U1, 2c: the tie toward the new reading (either side of rest, rest included) at the speed between samples (resync: the pull speed). The tie never jumps
            {
               p.DriveTie(s, resync ? p.pullRate : TieRate(s, prev, h.at, interval));
               Wake();
            }
            // U1: a reading at rest with the stops at 0 and the rod caught on them (its speed within the rest test's, RestWatch: kMoveM per cycle; the stops arrive at the
            // end of a cycle with the rod still at their speed, and the catch, restitution 0, stops it in the next): the tie lets go, the springs hold the rod at rest
            else if (s == 0.0 && h.target == 0.0 && h.at == 0.0 && std::fabs(p.rate) <= GroupStep::RestWatch::kMoveM / GroupStep::RestWatch::kCycleS)
            {
               p.Untie();
               hand = Hand::None;
               resync = false;
               Wake();
               return kSensorLetGoRest;
            }
            return kSensorNone;
         }
         // the hand off: a reading leaving the rest band (either side) ties the knob where the rod is, then drives it to the reading at 2c's rate from the first report
         // (2i c: no leaving-rest rule; a board's step at rest is removed by VPX's dead zone)
         if (s != 0.0)
         {
            FitRestart(s);
            p.Tie(p.q);
            p.DriveTie(s, TieRate(s, prev, p.q, interval));
            hand = Hand::Sensor;
            Wake();
            return kSensorTake;
         }
         return kSensorNone;
      }
      // 2k a: the tie takes the knob: the fit starts from this report (the pulls before it are another hand's)
      void FitRestart(const double s)
      {
         fitCount = 0;
         if (s > 0.0)
         {
            fitT[0] = sensorClock;
            fitS[fitCount++] = s;
         }
      }
      // 2k a: a least-squares quadratic s(t) over the pulled reports of the last W = sqrt(2 C q / release_window_noise) before the newest (at least 4), at the newest:
      // fitSpeed = -ds/dt, fitAccel = -d2s/dt2 (toward rest positive). False (no fit) under 4 reports or before the board's step q is known
      bool Fit()
      {
         if (fitCount < 4 || boardStep <= 0.0)
            return false;
         const double w = std::max(releaseWindowMin, std::sqrt(2.0 * releaseWindowC * boardStep / releaseWindowNoise)), tn = fitT[fitCount - 1];
         int i0 = 0;
         while (i0 < fitCount - 4 && fitT[i0] < tn - w - 1e-9)
            i0++;
         double S1 = 0.0, S2 = 0.0, S3 = 0.0, S4 = 0.0, Y0 = 0.0, Y1 = 0.0, Y2 = 0.0;   // tau in ms from the newest report (conditioning)
         const double S0 = fitCount - i0;
         for (int i = i0; i < fitCount; i++)
         {
            const double u = (fitT[i] - tn) * 1e3, u2 = u * u, y = fitS[i];
            S1 += u; S2 += u2; S3 += u2 * u; S4 += u2 * u2;
            Y0 += y; Y1 += u * y; Y2 += u2 * y;
         }
         // Cramer on the normal equations [S0 S1 S2; S1 S2 S3; S2 S3 S4] (b0 b1 b2) = (Y0 Y1 Y2)
         const double det = S0 * (S2 * S4 - S3 * S3) - S1 * (S1 * S4 - S3 * S2) + S2 * (S1 * S3 - S2 * S2);
         if (!(std::fabs(det) > 1e-12))
            return false;
         const double d1 = S0 * (Y1 * S4 - S3 * Y2) - Y0 * (S1 * S4 - S3 * S2) + S2 * (S1 * Y2 - Y1 * S2);
         const double d2 = S0 * (S2 * Y2 - Y1 * S3) - S1 * (S1 * Y2 - Y1 * S2) + Y0 * (S1 * S3 - S2 * S2);
         fitSpeed = -d1 / det * 1e3;           // m/ms -> m/s
         fitAccel = -2.0 * d2 / det * 1e6;     // m/ms^2 -> m/s^2
         fitN = fitCount - i0;
         fitSpan = tn - fitT[i0];
         return true;
      }
      // 2c's rule (F1): the step between the last two distinct samples over the time between them; at least the stops' gap to the new reading over that time (from the rod
      // when the tie takes it), so the tie reaches each sample by the next. restitution 0 (the hand's catch law, plungerrod.h): the rod takes the stops' speed.
      // e2m b: the catch-up is at most the key's pull speed over the knob's own (the agreed glide): the tie never moves faster than the knob plus the pull speed
      // (e2l: a lag closed in one 3 ms report drove the rod at 1.8 m/s behind a 1 m/s palm)
      double TieRate(const double s, const double prev, const double from, const double interval) const
      {
         const double step = std::fabs(s - prev);
         return std::min(std::max(step, std::fabs(s - from)), step + part.pullRate * interval) / interval;
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
      // 2i d: the hand's speed: the key's pull speed, or a hand stop's drive rate while it moves to its target when faster (the sensor's tie, up to about 1.3 m/s)
      double HandRate() const
      {
         double v = part.pullRate;
         for (const int i : { part.stopHold, part.stopPush })
            if (i >= 0 && part.stops[i].enabled && part.stops[i].at != part.stops[i].target)
               v = std::max(v, part.stops[i].rate);
         return v;
      }
      double SurfaceSpeedBound(const double cycleS) const override { return std::fabs(part.rate) + HandRate() + MaxStoredForce() / part.mass * cycleS; }
      // the energy reach (translating parts, as the Button cap): its springs and weight as the drive over its whole travel, the hand's move inside the speed bound
      double EnergyMove(const V3 &, const double E, const double T, double *const work) const override
      {
         const double F = MaxStoredForce(), D = F * (part.travel - part.outerSolidS), v0 = std::fabs(part.rate) + HandRate();
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
      // the part cull (design part-cull 2R): the tip's apex over the whole stroke s in [outerSolidS, travel] (pushed in to pulled out), widened to its place now plus or minus its move
      // this cycle (EnergyMove with the table's energy E) should it pass a stop (the stops are springs); every point of rod and tip within FarR(prof) of the apex: a capsule
      double cullLo = 0.0, cullHi = 0.0, cullFar = INFINITY;
      void CullCycle(const double E, const double T) override
      {
         const double m = EnergyMove(part.tipApex, E, T, nullptr);
         cullLo = std::min(part.outerSolidS, part.q - m);
         cullHi = std::max(part.travel, part.q + m);
         cullFar = BumperGroup::FarR(part.prof);
      }
      double CullDist(const V3 &c) const override
      {
         const double s = std::min(std::max((c - part.tipApex).Dot(part.axis), cullLo), cullHi);
         return (c - (part.tipApex + part.axis * s)).Length() - cullFar;
      }
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
