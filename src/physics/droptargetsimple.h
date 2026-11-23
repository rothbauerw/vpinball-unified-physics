// license:GPLv3+
//
// A drop target, the SIMPLE INTERIM model (STOCK OFF piece 1; runs/2026-11-21-simple-drop/design.md sections 1-4, agreed with the user 2026-11-21). The enhanced model
// (docs/unified-physics/DropTargetDesign.txt, branch drop-axis) replaces it; this file is deleted then.
//
// The target's drawn mesh (hat included) as a hard surface on the shared contact path, shifted down by z (m; 0 = up, -zDown = down, 52 VPX units). z is ONE coordinate with NO
// unknown: its move over a step is given (GroupStep::CoordInfo::given, the bat pivot's pattern, GroupStepDesign 19.2): the drop speed x h, clipped to the travel left, so a move
// ends exactly at 0 or -zDown. No mass, no stored energy, no stops, no friction of its own; the given move's work on the balls is booked by the group as a drive (the fall: the
// drop spring's, the rise: the reset coil's). Contact energy is the shared law's.
// Phases (state): up (z = 0, still: asleep by the rest rule) / falling / down (z = -zDown, not solid) / waiting (the raise delay, physics clock) / rising.
// Departures from Requirements I (design section 1): no speed law (VPX's Threshold, in the adapter), no spring-driven fall, no rebound from the parts, no slot hole.
// World frame in metres, z up. Pure SI, no VPX types.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "physics/contactmovingpart.h"

namespace DropSimple
{
   using GroupStep::Hit;
   using GroupStep::V3;

   enum Phase { kUp = 0, kFalling = 1, kDown = 2, kWaiting = 3, kRising = 4 };

   // its state (PartRecords)
   struct DropState
   {
      double z = 0.0;       // the mesh's shift along world z (m): 0 = up, -zDown = down
      double rate = 0.0;    // its mean rate over the last step or share alone taken (m/s)
      double given = 0.0;   // this attempt's given move over the step (m; BeginStep)
      double waitS = 0.0;   // waiting: the raise delay left (s, physics clock)
      int phase = kUp;
   };

   class Part : public GroupStep::Part, public DropState
   {
   public:
      Part() { KeepState<DropState>(this); }

      // ---- shape and settings (not state) ----
      std::vector<RubberContact::SurfaceTri> tris;   // the drawn mesh at the up pose (world, m)
      V3 origin;                                     // the mesh's centre at the up pose (m)
      double edgeRho = 1e-3;                         // the mesh edges' radius for the contact (the host's hard-surface edge_radius_mm), m
      double meshR = 0.0;                            // the farthest mesh point from origin (m)
      V3 boxLo, boxHi;                               // the mesh's box at the up pose (world, m; the part cull)
      double zDown = 0.0;                            // the full travel down (m): stock's DROP_TARGET_LIMIT, 52 VPX units
      double speed = 0.0;                            // the drop speed (m/s): VPX DropSpeed, read live by the adapter
      std::vector<ContactMovingPart::Event> *events = nullptr;   // the engine's queue (Dropped, Raised)
      SlideMesh mesh;                                // its probe along world z (contactmovingpart.h; [Contact] target_patches)

      void Build(const std::vector<RubberContact::SurfaceTri> &faces, const std::vector<std::array<int, 3>> &vid, const double zDownM, const double edgeRhoM, const bool startDown)
      {
         tris = faces;
         mesh.Build(faces, vid, V3(0.0, 0.0, 1.0), V3(0.0, 0.0, 1.0), edgeRhoM);
         zDown = zDownM;
         edgeRho = edgeRhoM;
         V3 lo(INFINITY, INFINITY, INFINITY), hi(-INFINITY, -INFINITY, -INFINITY);
         for (const RubberContact::SurfaceTri &t : tris)
            for (const V3 &p : { t.a, t.b, t.c })
            {
               lo = V3(std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z));
               hi = V3(std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z));
            }
         origin = (lo + hi) * 0.5;
         boxLo = lo;
         boxHi = hi;
         meshR = 0.0;
         for (const RubberContact::SurfaceTri &t : tris)
            for (const V3 &p : { t.a, t.b, t.c })
               meshR = std::max(meshR, (p - origin).Length());
         z = startDown ? -zDown : 0.0;   // (a target saved dropped starts down, not solid: design 4h)
         phase = startDown ? kDown : kUp;
         rate = given = waitS = 0.0;
      }

      bool Solid() const { return phase != kDown && phase != kWaiting; }                       // down and waiting: the balls pass
      bool Busy() const { return phase == kFalling || phase == kWaiting || phase == kRising; }   // a drive acts (the fall, the raise delay's clock, the rise)

      // The move over a time t from the state now (m), and whether it reaches the end of the travel: the raise delay counts down first, then the rise. Clipped to the travel left
      double Move(double t, bool &ends) const
      {
         ends = false;
         int ph = phase;
         if (ph == kWaiting)
         {
            if (waitS >= t)
               return 0.0;
            t -= waitS;
            ph = kRising;
         }
         if (ph != kFalling && ph != kRising)
            return 0.0;
         const double left = ph == kFalling ? z + zDown : -z;
         double d = speed * t;
         if (d >= left)
         {
            d = left;
            ends = true;
         }
         return ph == kFalling ? -d : d;
      }
      // Take the move over t (a step accepted, or a share alone): at the end of the travel z is exactly the end and the phase changes; Dropped / Raised queued (fired after the step)
      void Advance(const double t)
      {
         bool ends;
         const double d = Move(t, ends);
         if (phase == kWaiting)
         {
            if (waitS >= t)
            {
               waitS -= t;
               rate = 0.0;
               return;
            }
            waitS = 0.0;
            phase = kRising;
         }
         rate = t > 0.0 ? d / t : 0.0;
         if (!ends)
         {
            z += d;
            return;
         }
         const bool fell = phase == kFalling;
         z = fell ? -zDown : 0.0;
         phase = fell ? kDown : kUp;
         if (events != nullptr)
            events->push_back({ fell ? ContactMovingPart::Event::Dropped : ContactMovingPart::Event::Raised, -1, 0.0 });
      }

      // The probe at shift z + dz: the ball's centre moved by -(z + dz) along world z into the up frame, the fixed mesh probed there; a pure translation leaves the normal as
      // it is. One column: the depth's slope along z is n_contact . ez, the point's speed ez per unit rate, no spin (as StandupFace::Part::ProbeAt along its travel)
      void ProbeAt(const V3 &c, const double R, const double dz, Hit &h) const
      {
         mesh.Probe(c, R, z + dz, h);
      }

      // ---- the group's interface: one coordinate, no unknown, a given move ----
      int NumDof() const override { return 0; }
      int NumCoords() const override { return 1; }
      int DofOf(int) const override { return -1; }
      GroupStep::CoordInfo Coord(int) const override
      {
         GroupStep::CoordInfo c;
         c.rate = rate;
         c.given = given;    // its move over the step (BeginStep): a drive
         c.invM = 0.0;       // no mass
         c.invMFree = 0.0;   // no body a ball can push off
         c.reach = 1.0;      // a translation: m per m
         return c;
      }
      double Energy(const GroupStep::Ball *, int) const override { return 0.0; }   // nothing stored: no mass, no spring
      bool BeginStep(const GroupStep::Ball *, int, const double h, double *, int) override
      {
         bool ends;
         given = Move(h, ends);
         return true;
      }
      bool Eval(const double *, int, double, V3 *, V3 *, double *, GroupStep::Jac) override { return true; }   // (no rows of its own)
      // the step ends where the travel or the raise delay does (the shared step rule's time left)
      void StepFacts(const GroupStep::Ball *, int, const GroupStep::Settings &, GroupStep::StepFactList &f) override
      {
         if (phase == kWaiting && waitS > 0.0)
            f.TimeLeft(waitS, 1);
         else if ((phase == kFalling || phase == kRising) && speed > 0.0)
         {
            const double left = phase == kFalling ? z + zDown : -z;
            if (left > 0.0)
               f.TimeLeft(left / speed, 2);
         }
      }
      void Accept(const double *, int, const double h, bool, double &driveJ, double &lossJ, bool &through) override
      {
         driveJ = lossJ = 0.0;   // (the group books the given move's work on the balls as a drive)
         through = false;
         Advance(h);
      }
      void ScaleSpeeds(const double s) override { rate *= s; }
      std::string StartNote() const override
      {
         char s[200];
         snprintf(s, sizeof s, " drop z %.17g rate %.17g given %.17g wait %.17g phase %d", z, rate, given, waitS, phase);
         return s;
      }
      std::string HistoryNote() const override
      {
         char s[120];
         snprintf(s, sizeof s, " drop z %.6g mm phase %d", z * 1e3, phase);
         return s;
      }
   };

   // what the engine side keeps (PartRecords: part of its state)
   struct EngineState
   {
      bool collidable = true;          // false: the balls pass (the script's Collidable)
      bool idle = false;
      GroupStep::RestWatch rest;       // the rest rule's watch
   };

   // The drop target as a moving part of the contact system (contactmovingpart.h): the probe, the speed and reach bounds, its Hit, Dropped and Raised events, its rest
   class Engine : public ContactMovingPart, public EngineState
   {
   public:
      Engine()
      {
         KeepState<EngineState>(this);
         gp.events = &events;
      }
      Part gp;
      RubberContact::BallSI ball;
      V3 back;   // the target's back (unit, horizontal): away from the ball's side; a touch whose normal points against it is on the front

      void Build(const std::vector<RubberContact::SurfaceTri> &faces, const std::vector<std::array<int, 3>> &vid, const V3 &backDir, const double zDownM, const double edgeRhoM,
                 const bool startDown, const RubberContact::BallSI &b)
      {
         gp.Build(faces, vid, zDownM, edgeRhoM, startDown);
         back = backDir / backDir.Length();
         ball = b;
      }
      // the table side (between cycles: the adapter's Publish and Fire, the script)
      void SetCollidable(const bool c) { collidable = c; }
      void SetSpeed(const double v) { gp.speed = std::max(v, 0.0); }
      double Z() const { return gp.z; }
      bool Down() const { return gp.phase == kDown; }
      // the fall starts from where it is (a front Hit at or over the threshold while up or rising; the script's IsDropped = true at any phase but down)
      void Drop()
      {
         if (gp.phase == kDown || gp.phase == kFalling)
            return;
         gp.phase = kFalling;
         Wake();
      }
      // the script's IsDropped = false: the raise delay on the physics clock, then the rise from where it is
      void Raise(const double delayS)
      {
         gp.phase = kWaiting;
         gp.waitS = std::max(delayS, 0.0);
         Wake();
      }

      // ---- the moving-part interface ----
      void SetBall(const RubberContact::BallSI &b) override { ball = b; }
      void SetGravity(const V3 &) override {}   // (no mass: gravity does not load it)
      V3 Origin() const override { return gp.origin; }
      double Reach() const override { return gp.meshR + gp.zDown; }
      RubberContact::TubeProbe Probe(const V3 &c) const override
      {
         RubberContact::TubeProbe pr;
         if (!collidable || !gp.Solid())
         {
            pr.dist = INFINITY;
            pr.delta = -INFINITY;
            return pr;
         }
         GroupStep::Hit h;
         gp.ProbeAt(c, ball.R, 0.0, h);
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
         v = V3(0.0, 0.0, gp.rate);
         w = V3();
      }
      // the drop speed while a drive acts (falling, waiting, rising); 0 up or down
      double SurfaceSpeedBound(double) const override { return gp.Busy() ? gp.speed : 0.0; }
      // the energy reach: its given travel over T. work: 0, as the bat pivot's given slide (contactbatpart.h CycleRegionDist): the given move is not bounded by energy
      double EnergyMove(const V3 &, double, const double T, double *const work) const override
      {
         if (work != nullptr)
            *work = 0.0;
         bool ends;
         return std::fabs(gp.Move(T, ends));
      }
      double EnergyCap(double, double &, double &) const override { return 0.0; }   // nothing in its book (Part::Energy 0)
      // the region over a cycle of T: the up mesh's box shifted along z over [z, z + its given move over T], grown by the edges' radius; none while it cannot be solid in the cycle
      bool SweptBox(const double T, V3 &lo, V3 &hi) const
      {
         if (!collidable || !(gp.Solid() || (gp.phase == kWaiting && gp.waitS < T)))
            return false;
         bool ends;
         const double d = gp.Move(T, ends), z0 = gp.z + std::min(d, 0.0), z1 = gp.z + std::max(d, 0.0);
         const V3 g(gp.edgeRho, gp.edgeRho, gp.edgeRho);
         lo = V3(gp.boxLo.x, gp.boxLo.y, gp.boxLo.z + z0) - g;
         hi = V3(gp.boxHi.x, gp.boxHi.y, gp.boxHi.z + z1) + g;
         return true;
      }
      static double BoxDist(const V3 &lo, const V3 &hi, const V3 &c)   // the point-to-box distance (0 inside)
      {
         const double dx = std::max({ lo.x - c.x, 0.0, c.x - hi.x }), dy = std::max({ lo.y - c.y, 0.0, c.y - hi.y }), dz = std::max({ lo.z - c.z, 0.0, c.z - hi.z });
         return std::sqrt(dx * dx + dy * dy + dz * dz);
      }
      double CycleRegionDist(const V3 &c, const double T, double) const override
      {
         V3 lo, hi;
         return SweptBox(T, lo, hi) ? BoxDist(lo, hi, c) : INFINITY;
      }
      V3 cullLo, cullHi;
      bool cullOn = false;
      void CullCycle(const double E, const double T) override
      {
         cullE = E;
         cullT = T;
         cullOn = SweptBox(T, cullLo, cullHi);
      }
      double CullDist(const V3 &c) const override { return cullOn ? BoxDist(cullLo, cullHi, c) : INFINITY; }
      std::string StateNote() const override
      {
         char buf[160];
         std::snprintf(buf, sizeof buf, "drop z %.4f mm at %.4f m/s, phase %d, wait %.4f s, asleep %d", gp.z * 1e3, gp.rate, gp.phase, gp.waitS, (int)idle);
         return buf;
      }

      // ---- the group's calls ----
      GroupStep::Part *GroupPart() override { return &gp; }
      int GroupSurfaces() const override { return 1; }
      void GroupProbe(int, const V3 &c, const double R, const double *const trial, GroupStep::Hit &h) const override
      {
         gp.ProbeAt(c, R, trial != nullptr ? trial[0] : 0.0, h);
         if (!collidable || !gp.Solid())
         {
            h.delta = -1.0;
            h.npat = 0;   // (no depth and no patches: no friction group either)
            h.npc = 0;
         }
      }
      bool GroupRestCentre(int, const V3 &c, const double *const move, V3 &out) const override   // (SlideMesh: the probe's own centre, c - (z + dz) a)
      {
         out = c - gp.mesh.a * (gp.z + (move != nullptr ? move[0] : 0.0));
         return true;
      }
      void GroupLinkWindows(int, const V3 &c0r, const V3 &c1r, const double R, std::vector<RubberContact::LinkWindow> &out) const override { gp.mesh.LinkWindows(c0r, c1r, R, out); }
      bool GroupPatchZero(int) const override { return gp.mesh.PatchZeroOn(); }   // (patch_zero: per-patch zeros keyed by face id, P3b piece 4)
      void CycleBegin() override { idle = false; }
      // Hit: a touch begins on the front (the touch's normal, toward the ball, against the back: drop-axis's test) while it is solid; the event carries the ball's own id.
      // The threshold and the drop are the adapter's, after the step
      void TouchBegan(int, int, const int ballId, const double vn, const V3 &n) override
      {
         if (gp.Solid() && -n.Dot(back) > 0.0)
            events.push_back({ Event::Hit, ballId, vn });
      }
      bool Sleeping() const override { return idle; }
      bool RestNumbers(double &moveM, double &speedMS) const override
      {
         moveM = rest.CycleMoved(RestPose());
         speedMS = std::fabs(gp.rate);
         return true;
      }
      bool AtRest() const override { return idle || (!gp.Busy() && rest.Rests(RestPose(), std::fabs(gp.rate))); }
      GroupStep::RestWatch *Rest() override { return &rest; }
      GroupStep::RestWatch::Pose RestPose() const override { return GroupStep::RestWatch::Pose().Add(gp.z, 1.0); }
      // No ball in its group: the given move alone (the fall, the raise delay, the rise); still and undriven, it sleeps
      int AdvanceAlone(const double dt) override
      {
         if (dt <= 0.0 || idle)
            return 0;
         gp.Advance(dt);
         idle = rest.Step(RestPose(), std::fabs(gp.rate), gp.Busy(), restTouched, dt);
         return 0;
      }

   private:
      void Wake()
      {
         idle = false;
         rest.Reset();
      }
   };
}
