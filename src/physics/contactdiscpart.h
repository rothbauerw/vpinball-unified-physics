// license:GPLv3+
//
// A turntable disc as a moving part of the contact system, on the group step. Its motor, bearing and shape are a DiscGroup::Part (discgroup.h, bodypart.h): the engine moves it with the balls in one solve,
// so an island with one always goes to the group. There is no other path: no old-loop contact, no impulse, no stock.
//
// What this class adds to the group part: what the engine needs of any moving part (the probe for engagement, the saved state), its contact law, and its own steps when no ball is near.
// It raises no events: a disc comes from the ini and has no table object to raise them on. World frame in metres, z up.

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>

#include "physics/contactmovingpart.h"
#include "physics/discgroup.h"

namespace DiscModel
{
   // what the engine side keeps (PartRecords: part of its state)
   struct EngineState
   {
      bool idle = false;
      GroupStep::RestWatch rest;   // the rest rule's watch, alone and over its last cycle in a group (RestPose)
   };

   class Engine : public ContactMovingPart, public EngineState
   {
   public:
      Engine() { KeepState<EngineState>(this); }
      DiscGroup::Part gp;
      RubberContact::RubberParams rp;  // the textured surface against the steel ball
      double layerM = 0.0;
      RubberContact::BallSI ball;
      long failedSteps = 0, aloneSteps = 0, energySteps = 0;   // the disc alone: a step not solved (the detector: rule 8), steps taken (cost), of them taken on energy alone (rule 8)
      int bookMs = 0;                  // the table's time, for the book of a cycle alone

      // The script's drive, between cycles: a disc asleep wakes, and who has an unknown is decided again at the next CycleBegin
      void SetVolts(const double v)
      {
         gp.volts = v;
         idle = false;
      }

      // ---- the moving-part interface (contactmovingpart.h) ----
      void SetBall(const RubberContact::BallSI &b) override { ball = b; }
      void SetGravity(const V3 &) override {}   // the disc turns about its own axis: its weight acts on the bearing, not about the axis
      V3 Origin() const override { return gp.centre; }
      double Reach() const override { return gp.R; }
      double SurfaceSpeedBound(double) const override { return 0.0; }   // no point of the disc moves along its own normal (discgroup.h)
      double MotionRange() const override { return 0.0; }               // (nor over any time: turning moves no gap)
      RubberContact::TubeProbe Probe(const V3 &c) const override
      {
         GroupStep::Hit h;
         gp.ProbeSurface(c, ball.R, h);
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
      void SurfaceVelocity(const V3 &point, V3 &v, V3 &w) const override
      {
         v = V3(-gp.rate * (point.y - gp.centre.y), gp.rate * (point.x - gp.centre.x), 0.0);
         w = V3(0.0, 0.0, gp.rate);
      }
      // The old loop's calls: an island with this part never goes there (GroupPart)

      // ---- the group's calls ----
      GroupStep::Part *GroupPart() override { return &gp; }
      void GroupProbe(int, const V3 &c, const double R, const double *, GroupStep::Hit &h) const override { gp.ProbeSurface(c, R, h); }
      double GroupPathCoverSag(int, int &kind) const override   // (the self-check, VPX_FINDER_CHECK)
      {
         kind = 3;
         return gp.pathCover.MaxSag();
      }
      int GroupPathReach(int, const V3 &s, const V3 &d, const double R) const override   // [Features] path_test: the surface never moves (DiscGroup::Part::PathReach)
      {
         const bool prof = ContactProf::CheckOn();
         bool band = false;
         const int r = gp.PathReach(s, d, R, prof, band);
         if (prof)
         {
            RubberContact::PathStats &ps = RubberContact::pathStats;
            ps.kTested[3]++;
            ps.kReached[3] += r > 0;
            ps.kBand[3] += band;
         }
         return r;
      }
      const RubberContact::RubberParams *GroupLaw(int, double &layer) const override
      {
         layer = layerM;
         return &rp;
      }
      void CycleBegin() override
      {
         idle = false;
         gp.CycleBegin();
      }
      bool Sleeping() const override { return idle; }
      // Nothing moved over the last cycle in a group: it stays out of a ball's solve while no ball can reach it. A part out of the solve still moves alone
      bool RestNumbers(double &moveM, double &speedMS) const override
      {
         moveM = rest.CycleMoved(RestPose());
         speedMS = std::fabs(gp.rate) * gp.reach;
         return true;
      }
      bool AtRest() const override
      {
         return gp.volts == 0.0 && (idle || rest.Rests(RestPose(), std::fabs(gp.rate) * gp.reach));   // the one rest rule: the motor drives it; its move over its last cycle in a group, and its speed now
      }
      GroupStep::RestWatch *Rest() override { return &rest; }
      GroupStep::RestWatch::Pose RestPose() const override { return GroupStep::RestWatch::Pose().Add(gp.q, gp.reach); }   // its one rest pose: the angle at the reach (the motor's current is not modelled: no internal state)
      // No ball near: the same group with no balls, asleep while the friction holds it still; the script's drive wakes it
      int AdvanceAlone(const double dt) override
      {
         if (dt <= 0.0 || idle)
            return 0;
         const bool was = idle;
         gp.ClearLoads();
         gp.CycleBegin();
         const int steps = GroupStep::AloneCycle(gp, dt, bookMs, failedSteps, aloneSteps, energySteps);
         idle = rest.Step(RestPose(), std::fabs(gp.rate) * gp.reach, gp.volts != 0.0, restTouched, dt);   // the one rest rule: driven while the motor has volts
         GroupStep::RestWatch::Book(2, was, idle);
         return steps;
      }

      std::string StateNote() const override
      {
         char buf[160];
         std::snprintf(buf, sizeof buf, "%.2f rad/s, volts %.2f, stuck %d, asleep %d", gp.rate, gp.volts, (int)gp.Stuck(), (int)idle);
         return buf;
      }
   };
}
