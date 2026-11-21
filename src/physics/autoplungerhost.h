// license:GPLv3+
//
// The auto plunger (autoplunger.h) as a moving part of the contact system, and its settings from the unified vpx_physics.ini (docs/unified-physics/PlungerDesign.txt
// section 4; engine step 3, runs/2026-11-17-plunger/a3-vpx/brief.md). Modelled on the mechanical rod's host (plungerhost.h Rod): the probe, the speed bound and energy
// reach, its own steps when no ball is near, the coil's Fire between cycles. No VPX types.

#pragma once

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

#include "physics/autoplunger.h"
#include "physics/plungerhost.h"

namespace AutoPlungerHost
{
   using V3 = RubberContact::V3;

   // [Plungers] for all, then [Plunger:<name>]: the auto plunger's keys (each one's source: the comment beside its field in autoplunger.h). The tip head takes the
   // mechanical rod's tip_* keys (plungerhost.h Load). A value that cannot build a working part is refused: the key named in `bad` (empty = all good)
   inline void Load(AutoPlunger::Settings &s, const std::filesystem::path &file, const std::string &name, std::vector<std::string> &bad)
   {
      PartsIni::ReadSections(file, "plungers", "plunger:" + name, [&](const std::string &k, double v, const std::string &lv) {
         if (k == "auto_energy_j") s.autoEnergyJ = v;                          // the coil's work over the stroke, J (ESTIMATED 1.0; per table)
         else if (k == "auto_head") s.head = lv == "tip" ? AutoPlunger::Head::Tip : AutoPlunger::Head::Fork;   // fork (default) or tip (design 4h)
         else if (k == "auto_mass_kg") s.massG = v * 1e3;                      // the moving mass (carriage, head, coil plunger, linkage), kg (ESTIMATED 0.080)
         else if (k == "auto_stroke_mm") s.strokeMm = v;                       // the carriage's travel, rest stop to end stop, mm (ESTIMATED 27)
         else if (k == "auto_return_preload_n") s.returnPreloadN = v;          // the return spring's force at rest, N (ESTIMATED 1.0)
         else if (k == "auto_return_rate_nmm") s.returnRateNmm = v;            // the return spring's rate, N/mm (ESTIMATED 0.1)
         else if (k == "rubber") { if (v == 0.0) ContactSystem::NoticeRetired("rubber = 0 ([Plungers] / [Plunger:<name>]) is ignored: no stock path on the unified engine"); }   // (once the opt-out to stock)
      });
      {
         std::vector<std::string> rodBad;   // the tip head's own keys (the rod's other keys are not used here)
         PlungerHost::Load(s.tip, file, name, rodBad);
         for (const std::string &k : rodBad)
            if (k.rfind("tip_", 0) == 0)
               bad.push_back(k);
      }
      const auto need = [&](const bool ok, const char *key) {
         if (!ok)
            bad.emplace_back(key);
      };
      need(s.autoEnergyJ > 0.0, "auto_energy_j");
      need(s.massG > 0.0, "auto_mass_kg");
      need(s.strokeMm > 0.0, "auto_stroke_mm");
      need(s.returnPreloadN >= 0.0, "auto_return_preload_n");
      need(s.returnRateNmm >= 0.0, "auto_return_rate_nmm");
   }

   // the engine side (PartRecords: part of its state)
   struct EngineState
   {
      bool idle = false;
      GroupStep::RestWatch rest;
   };

   class Auto : public ContactMovingPart, public EngineState
   {
   public:
      Auto() { KeepState<EngineState>(this); }
      AutoPlunger::Part part;
      RubberContact::BallSI ball;
      long failedSteps = 0, aloneSteps = 0, energySteps = 0;   // the part alone: a step not solved (rule 8), steps taken (cost), of them taken on energy alone (rule 8)
      int bookMs = 0;
      double tipSlope = 1.0;
      V3 mid;            // the head's reference point at the middle of the stroke (the host's box, Origin)
      double headR = 0.0;   // every point of the head within this of its reference point

      // axis: +s (toward the ball); up: the playfield's normal; cradle: the centre of a ball that touches the head at s = 0
      void Build(const AutoPlunger::Settings &s, const V3 &axis, const V3 &up, const V3 &gravity, const V3 &cradle, const RubberContact::BallSI &b)
      {
         ball = b;
         part.Build(s, axis, up, gravity, cradle, b.R, b.m);
         if (part.head == AutoPlunger::Head::Fork)
         {
            const V3 c = (part.tine[0] + part.tine[1]) * 0.5;
            headR = (part.tine[0] - c).Length() + part.rTine;
            mid = c + axis * (0.5 * part.stroke);
         }
         else
         {
            tipSlope = part.prof.Slope();
            headR = BumperGroup::FarR(part.prof);
            mid = part.tipApex + axis * (0.5 * part.stroke);
         }
         idle = true;   // on its rest stop: nothing to step until a ball comes near, the coil fires or the gravity changes
      }
      // the head's reference point now (the fork's tines' middle or the tip's apex)
      V3 HeadNow() const { return mid + part.axis * (part.q - 0.5 * part.stroke); }
      // the largest force on the carriage anywhere in its stroke (N): the drive, the return spring at the end, its weight
      double MaxForce() const { return part.force + part.pRet + part.kRet * part.stroke + part.mass * part.g.Length(); }

      // the coil (between cycles): the drive pushes the carriage to the end stop
      void Fire()
      {
         part.Fire();
         idle = false;
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
      double Reach() const override { return headR + part.stroke; }   // every point of the head at any travel
      double MotionRange() const override { return part.stroke; }
      double SurfaceSpeedBound(const double cycleS) const override { return std::fabs(part.rate) + MaxForce() / part.mass * cycleS; }
      // the energy reach (translating parts, as the rod): the largest force as the drive over the whole stroke
      double EnergyMove(const V3 &, const double E, const double T, double *const work) const override
      {
         const double F = MaxForce(), D = F * part.stroke, v0 = std::fabs(part.rate);
         const double x = DriveCappedReach(v0, E, part.mass, F / part.mass, D, T);
         if (work != nullptr)
            *work = 0.5 * part.mass * v0 * v0 + std::min(F * x, D);
         return x;
      }
      // the one-pass bound (EnergyCap, AGREED items 1, 4b, 5b): its energy in the book now (AutoPlunger: Energy = kinetic + StoredEnergy (return spring and weight) + pin) less
      // StoredEnergy's least over the stroke [0, stroke] (a convex quadratic: an end or where StoredForce is 0), + the coil: this model's coil is the stand-in drive F over the stroke
      // (autoplunger.h 4d, no current or resistance), so what it can still give is F (stroke - q) while it fires, 0 when off. Its stops are Hertz: aHertz (BodyPart::HertzSquashA)
      double EnergyCap(double, double &, double &aHertz) const override
      {
         const double hi = part.stroke;
         const auto at = [&](const double s) { return part.StoredEnergy(std::min(std::max(s, 0.0), hi)); };
         double least = std::min(at(0.0), at(hi));
         if (part.kRet > 0.0)
            least = std::min(least, at((part.GravityForce() - part.pRet) / part.kRet));
         aHertz = part.HertzSquashA();
         return part.Energy(nullptr, 0) - least + (part.firing ? part.force * std::max(hi - part.q, 0.0) : 0.0);
      }
      void HeadProbe(const int i, const V3 &c, const double R, const double dq, GroupStep::Hit &h) const
      {
         if (part.head == AutoPlunger::Head::Fork)
            part.ProbeTine(i, c, R, dq, h);
         else
            part.ProbeTip(c, R, dq, h);
      }
      RubberContact::TubeProbe Probe(const V3 &c) const override
      {
         GroupStep::Hit h;
         HeadProbe(0, c, ball.R, 0.0, h);
         if (part.head == AutoPlunger::Head::Fork)   // the nearer tine
         {
            GroupStep::Hit h1;
            HeadProbe(1, c, ball.R, 0.0, h1);
            if (h1.delta > h.delta)
               h = h1;
         }
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
      int GroupSurfaces() const override { return part.head == AutoPlunger::Head::Fork ? 2 : 1; }   // the fork: each tine its own contact
      void GroupProbe(const int i, const V3 &c, const double R, const double *const trial, GroupStep::Hit &h) const override { HeadProbe(i, c, R, trial != nullptr ? trial[0] : 0.0, h); }
      double GroupPathMotion(int, const V3 &, double, const double *const trial) const override { return trial != nullptr ? std::fabs(trial[0]) : 0.0; }   // a slide: |dq|
      double RegionDist(const V3 &c) const override { return (c - HeadNow()).Length() - (headR + part.stroke); }
      // the part cull (design part-cull 2R, as the rod): the head's reference point over the whole stroke s in [0, stroke], widened to its place now plus or minus its move this cycle
      // (EnergyMove with the table's energy E) should it pass a stop; every point of the head within headR of it: a capsule
      double cullLo = 0.0, cullHi = 0.0;
      void CullCycle(const double E, const double T) override
      {
         const double m = EnergyMove(mid, E, T, nullptr);
         cullLo = std::min(0.0, part.q - m);
         cullHi = std::max(part.stroke, part.q + m);
      }
      double CullDist(const V3 &c) const override
      {
         const V3 base = mid - part.axis * (0.5 * part.stroke);   // (s = 0)
         const double s = std::min(std::max((c - base).Dot(part.axis), cullLo), cullHi);
         return (c - (base + part.axis * s)).Length() - headR;
      }
      double GroupPathSlope(int) const override { return part.head == AutoPlunger::Head::Fork ? 1.0 : tipSlope; }   // a sphere's distance moves at most 1:1 with the slide
      const RubberContact::RubberParams *GroupLaw(int, double &layerM) const override
      {
         if (part.head == AutoPlunger::Head::Fork)
         {
            layerM = 0.0;
            return &part.forkRp;
         }
         layerM = part.tipLayerM;
         return &part.tipRp;
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
      // the one rest rule; never at rest while the coil drives
      bool AtRest() const override { return !part.firing && (idle || rest.Rests(RestPose(), std::fabs(part.rate))); }
      GroupStep::RestWatch *Rest() override { return &rest; }
      GroupStep::RestWatch::Pose RestPose() const override { return GroupStep::RestWatch::Pose().Add(part.q, 1.0); }
      // No ball near: the same group with no balls, asleep while nothing moves and the coil is off
      int AdvanceAlone(const double dt) override
      {
         if (dt <= 0.0 || idle)
            return 0;
         const bool was = idle;
         part.ClearLoads();
         part.CycleBegin();
         const int steps = GroupStep::AloneCycle(part, dt, bookMs, failedSteps, aloneSteps, energySteps);
         idle = rest.Step(RestPose(), std::fabs(part.rate), part.firing, restTouched, dt);
         GroupStep::RestWatch::Book(1, was, idle);
         return steps;
      }

      std::string StateNote() const override
      {
         char buf[200];
         std::snprintf(buf, sizeof buf, "auto plunger: travel %.4f mm at %.4f m/s, firing %d, asleep %d", part.q * 1e3, part.rate, (int)part.firing, (int)idle);
         return buf;
      }
   };
}
