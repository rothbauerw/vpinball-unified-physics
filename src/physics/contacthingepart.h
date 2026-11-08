// license:GPLv3+
//
// A part on one horizontal axis as a moving part of the contact system, on the group step (a gate, a spinner; runs/2026-11-02-gate/notes.md). Its shape, stops, friction and switch are a HingeGroup::Part (hingegroup.h, bodypart.h): the engine moves it with the balls in one solve,
// so an island with one always goes to the group. There is no other path: no old-loop contact, no impulse, no stock.
//
// What this class adds to the group part: what the engine needs of any moving part (the probe for engagement and the overlap replay, the speed bound, the saved state), its surfaces, its own steps when no ball is near, and its Hit, switch and stop events.
// Hit: a ball starts touching the part (the contact system keeps the release gap's hysteresis) while the part is near its start angle, read at the step's start: a light part turns several degrees in the step a touch starts in.
// World frame in metres, z up.

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>

#include "physics/contactmovingpart.h"
#include "physics/hingegroup.h"

namespace HingeModel
{
   class Engine : public ContactMovingPart
   {
   public:
      explicit Engine(HingeGroup::Part &part) : gp(part) {}
      HingeGroup::Part &gp;
      RubberContact::BallSI ball;
      bool collidable = true;          // false: the balls pass; it keeps moving under gravity and its drive
      long failedSteps = 0, aloneSteps = 0, energySteps = 0;   // the part alone: a step not solved (the detector: rule 8), steps taken (cost), of them taken on energy alone (rule 8)
      int bookMs = 0;                  // the table's time, for the book of a cycle alone

      void SetCollidable(const bool c) { collidable = c; }
      void Wake() { idle = false; }
      // The switch's changes, the spins and the stops' touches as the part's events. The adapter calls it before the host takes the events.
      void Collect()
      {
         for (const char e : gp.swLog)
            events.push_back({ e == 2 ? Event::Spin : e == 1 ? Event::SwitchClosed : Event::SwitchOpened, -1, 0.0 });
         gp.swLog.clear();
         for (const GroupStep::BodyStopEvent &e : gp.stopEvents)
         {
            const int kind = gp.LimitKind(e.stop);
            if (kind != 0)
               events.push_back({ kind > 0 ? Event::LimitEOS : Event::LimitBOS, -1, e.speed * 180.0 / GroupStep::kPi });
         }
         gp.stopEvents.clear();
      }

      // ---- the moving-part interface (contactmovingpart.h) ----
      void SetBall(const RubberContact::BallSI &b) override { ball = b; }
      void SetGravity(const V3 &g) override
      {
         if ((g - gp.g).Length() > 1e-6)
            idle = false;   // (a nudge wakes it)
         gp.g = g;
      }
      V3 Origin() const override { return gp.A; }
      double Reach() const override { return gp.reach; }
      // From its angle to the farther of its stops (a driven stop counted at the far end of its move), as a chord at its reach; a side with no stop: any angle
      double MotionRange() const override
      {
         double lo = INFINITY, hi = -INFINITY;
         for (const GroupStep::BodyStop &s : gp.stops)
            if (s.enabled)
            {
               if (s.side > 0)
                  hi = std::max(hi, s.rate != 0.0 ? std::max(s.at, s.target) : s.at);
               else
                  lo = std::min(lo, s.rate != 0.0 ? std::min(s.at, s.target) : s.at);
            }
         if (!std::isfinite(lo) || !std::isfinite(hi))
            return 2.0 * gp.reach;   // a side with no stop: any angle
         const double sweep = std::max(gp.q - lo, hi - gp.q);
         return sweep >= GroupStep::kPi ? 2.0 * gp.reach : 2.0 * gp.reach * std::sin(0.5 * std::max(sweep, 0.0));
      }
      // Its speed now plus what gravity, the blades and a driven stop can add in a cycle
      double SurfaceSpeedBound(const double cycleS) const override
      {
         double drive = 0.0;
         for (const GroupStep::BodyStop &s : gp.stops)
            drive = std::max(drive, std::fabs(s.SpeedNow()));
         return (std::fabs(gp.rate) + gp.MaxStoredTorque() / gp.mass * cycleS + drive) * gp.reach;
      }
      RubberContact::TubeProbe Probe(const V3 &c) const override
      {
         RubberContact::TubeProbe pr;
         if (!collidable)
         {
            pr.dist = INFINITY;
            pr.delta = -INFINITY;
            return pr;
         }
         GroupStep::Hit best;
         int surf = -1;
         for (int sf = 0; sf < gp.Surfaces(); sf++)
         {
            GroupStep::Hit h;
            gp.ProbeSurface(sf, c, ball.R, 0.0, h);
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
      void SurfaceVelocity(const V3 &point, V3 &v, V3 &w) const override
      {
         w = gp.E * gp.rate;
         v = w.Cross(point - gp.A);
      }
      // The old loop's calls: an island with this part never goes there (GroupPart)
      RubberContact::ContactResult ContactItem(int, SurfaceContact &, const V3 &, const V3 &, const V3 &, double) override { assert(false); return RubberContact::ContactResult(); }
      double BallStep(const V3 &, const V3 &, double hMax) const override { assert(false); return hMax; }
      void StepInIsland(double) override { assert(false); }

      // ---- the group's calls ----
      GroupStep::Part *GroupPart() override { return &gp; }
      int GroupSurfaces() const override { return gp.Surfaces(); }
      void GroupProbe(const int surface, const V3 &c, const double R, const double *const trial, GroupStep::Hit &h) const override
      {
         gp.ProbeSurface(surface, c, R, trial != nullptr ? trial[0] : 0.0, h);
         if (!collidable)
            h.delta = -1.0;
      }
      bool GroupNearPieces(const int surface, const V3 &c, const double R, const double reach, const double swing, const double minReach, std::vector<int> &list) const override
      {
         if (!gp.pieceSum)
            return false;
         gp.NearPieces(surface, c, R, reach, swing, minReach, list);
         return true;
      }
      bool GroupNearHolds(const int surface, const V3 &c, const double R, const double dq, const std::vector<int> &list) const override { return gp.NearHolds(surface, c, R, dq, list); }
      const RubberContact::RubberParams *GroupLaw(int, double &layerM) const override
      {
         layerM = 0.0;
         return &gp.rp;
      }
      void CycleBegin() override
      {
         idle = false;
         markQ = gp.q;
         marked = true;
         gp.CycleBegin();
      }
      void TouchBegan(int, const int ballIndex, const double vn) override
      {
         if (gp.NearStart(gp.qAtStepStart))
            events.push_back({ Event::Hit, ballIndex, vn });
      }
      bool Sleeping() const override { return idle; }
      // Nothing driven and nothing moved over the last cycle in a group: it stays out of a ball's solve while no ball can reach it. A part out of the solve still moves alone
      bool RestNumbers(double &moveM, double &speedMS) const override
      {
         moveM = marked ? std::fabs(gp.q - markQ) * gp.reach : -1.0;
         speedMS = std::fabs(gp.rate) * gp.reach;
         return true;
      }
      bool AtRest() const override
      {
         if (gp.AnyStopMoving())
            return false;
         return idle || (marked && std::fabs(gp.q - markQ) * gp.reach < 1e-9 && std::fabs(gp.rate) * gp.reach < 1e-6);
      }
      // No ball near: the same group with no balls, asleep while nothing is driven or moving, its state kept as it is
      int AdvanceAlone(const double dt) override
      {
         if (dt <= 0.0 || (idle && !gp.AnyStopMoving()))
            return 0;
         gp.ClearLoads();
         gp.CycleBegin();
         const int steps = GroupStep::AloneCycle(gp, dt, bookMs, failedSteps, aloneSteps, energySteps);
         idle = gp.NumDof() == 0 && gp.rate == 0.0 && !gp.AnyStopMoving();
         return steps;
      }

      // The dynamic state, plain bytes (the settings and the shape are not state). Restoring drops the events since the save (the replay makes them again) and keeps those before it
      struct Dyn
      {
         GroupStep::BodyPart::Dyn b;
         V3 g;
         double steppedS, markQ, prevWrapped;
         bool idle, steppedThisCycle, collidable, marked, switchClosed;
         uint32_t nSwLog, nStopEvents, nEvents;   // the events made before the save: a restore keeps them (only those since are made again by the replay)
      };
      void SaveState(State &out) const override
      {
         Dyn d;
         std::memset(&d, 0, sizeof d);   // (the padding too: saved bytes must be the same for the same state)
         gp.Fill(d.b);
         d.g = gp.g;
         d.steppedS = steppedS;
         d.markQ = markQ;
         d.prevWrapped = gp.prevWrapped;
         d.idle = idle;
         d.steppedThisCycle = steppedThisCycle;
         d.collidable = collidable;
         d.marked = marked;
         d.switchClosed = gp.switchClosed;
         d.nSwLog = (uint32_t)gp.swLog.size();
         d.nStopEvents = (uint32_t)gp.stopEvents.size();
         d.nEvents = (uint32_t)events.size();
         out.used = 0;
         Put(out, d);
      }
      void RestoreState(const State &in) override
      {
         StateReader r(in);
         Dyn d;
         r.Get(d);
         gp.Load(d.b);
         gp.g = d.g;
         gp.switchClosed = d.switchClosed;
         gp.prevWrapped = d.prevWrapped;
         gp.swLog.resize(std::min<size_t>(gp.swLog.size(), d.nSwLog));
         gp.stopEvents.resize(std::min<size_t>(gp.stopEvents.size(), d.nStopEvents));
         steppedS = d.steppedS;
         markQ = d.markQ;
         idle = d.idle;
         steppedThisCycle = d.steppedThisCycle;
         collidable = d.collidable;
         marked = d.marked;
         events.resize(std::min<size_t>(events.size(), d.nEvents));
      }
      std::string StateNote() const override
      {
         char buf[200];
         std::snprintf(buf, sizeof buf, "angle %.2f deg at %.2f rad/s, on stop %d, stuck %d, switch %d, asleep %d", gp.q * 180.0 / GroupStep::kPi, gp.rate, gp.PinStop(), (int)gp.Stuck(), (int)gp.switchClosed, (int)idle);
         return buf;
      }

   private:
      bool idle = false, marked = false;
      double markQ = 0.0;   // the angle at the start of its last cycle in a group (AtRest)
   };
}
