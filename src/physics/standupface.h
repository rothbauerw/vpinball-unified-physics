// license:GPLv3+
//
// A stand-up target's face on its blade (piece 2; runs/2026-11-19-standup-p2/design.md sections 1 and 4): ONE massless coordinate x, the face's travel along its normal away from
// the ball (m), against the blade's spring k x, k = 3 E I / L^3 of the blade's sizes. No losses of its own, no tilt, twist, foam, plate or switch (pieces 3-6). The ball-face contact is
// the shared law (GroupStep::Item) on the target's own mesh, moved by x n.
//
// A massless coordinate (CoordInfo invM 0, with an unknown): its row is the force balance at the step's end, h (k x1 - the contacts' push at the end), written by the group
// (GroupStep::Group::Eval, the massless row) beside the part's own h k x1; the ball pays the blade's energy change (the group's massless energy correction). With no ball in its
// group the face rests at x = 0 and is not stepped (a massless spring with no load sits at its free length).
// World frame in metres, z up. Pure SI, no VPX types.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "physics/contactmovingpart.h"

namespace StandupFace
{
   using GroupStep::Hit;
   using GroupStep::V3;

   // its settings: [Targets] for every stand-up target, then [Target:<name>] for one (PartsIni::LoadTarget)
   struct Settings
   {
      bool faceTravel = true;          // face_travel: 1 = the face travels on its blade (this piece); 0 = piece 1's fixed face (the mesh as a fixed hard surface)
      double bladeThicknessIn = 0.020; // blade_thickness_in: the blade's thickness, inches (PartData T8c: 0.020; 0.025 is the stiffer blade). Measured part size
      double bladeLeverIn = 1.625;     // blade_lever_in: from the blade's clamp to the face centre, inches (PartData T8c). Measured part size
   };

   // its state (PartRecords)
   struct FaceState
   {
      double x = 0.0;      // the face's travel along n, away from the ball (m); 0 = at rest on the blade's free length
      double rate = 0.0;   // its mean rate over the last step taken (m/s); 0 alone
      bool unk = false;    // it has an unknown: a ball is in its group (CycleBegin, BeginStep); none alone (AdvanceAlone)
   };

   class Part : public GroupStep::Part, public FaceState
   {
   public:
      Part() { KeepState<FaceState>(this); }

      // ---- shape and settings (not state) ----
      std::vector<RubberContact::SurfaceTri> tris;   // the target's own mesh at rest (world, m): piece 1's faces
      V3 n;                                          // the travel direction, unit: the face's normal pointing away from the ball's side
      V3 origin;                                     // the mesh's centre at rest (m)
      double k = 888.0;                              // the blade's spring, N/m (3 E I / L^3)
      double edgeRho = 1e-3;                         // the mesh edges' radius for the contact (the host's hard-surface edge_radius_mm), m
      double meshR = 0.0;                            // the farthest mesh point from origin (m)
      V3 boxLo, boxHi;                               // the mesh's box at rest (world, m; the part cull)

      // the blade's spring from its sizes (PartData T8c, design section 1): a steel cantilever E 200 GPa (ASSUMED), 3/8 in wide, thick and long as set, loaded at the face centre
      static double BladeK(const double thicknessIn, const double leverIn)
      {
         constexpr double E = 200e9, in = 0.0254, w = 0.375 * in;
         const double t = thicknessIn * in, L = leverIn * in, I = w * t * t * t / 12.0;
         return 3.0 * E * I / (L * L * L);
      }
      SlideMesh mesh;                                // its probe along n (contactmovingpart.h: the per-patch path)

      void Build(const std::vector<RubberContact::SurfaceTri> &faces, const std::vector<std::array<int, 3>> &vid, const V3 &travel, const double kBlade, const double edgeRhoM)
      {
         tris = faces;
         mesh.Build(faces, vid, travel, travel / travel.Length() * -1.0, edgeRhoM);
         n = travel / travel.Length();
         k = kBlade;
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
         x = rate = 0.0;
         unk = false;
      }

      // The probe at travel x + dx (design section 4, as the gate's ProbeSurface): the ball's centre moved by -(x + dx) n into the rest frame, the fixed mesh probed there; a pure
      // translation leaves the normal as it is. One column: the depth's slope along x is n_contact . n (d/dx of R - dist(c - x n)), the point's speed n per unit rate, no spin
      void ProbeAt(const V3 &c, const double R, const double dx, Hit &h) const
      {
         mesh.Probe(c, R, x + dx, h);
      }

      // ---- the group's interface ----
      int NumDof() const override { return unk ? 1 : 0; }
      int MaxDof() const override { return 1; }
      int NumCoords() const override { return 1; }
      int DofOf(int) const override { return unk ? 0 : -1; }
      GroupStep::CoordInfo Coord(int) const override
      {
         GroupStep::CoordInfo c;
         c.rate = rate;
         c.given = 0.0;      // with no unknown: at rest
         c.invM = 0.0;       // massless (the group's massless row)
         c.invMFree = 0.0;
         c.reach = 1.0;      // a translation: m per m
         return c;
      }
      double Energy(const GroupStep::Ball *, int) const override { return 0.5 * k * x * x; }
      // the blade's energy change over the move x[off] (the group's massless energy correction): 1/2 k (x1^2 - x0^2) in the form with no cancellation; its slope k x1 into slope[off]
      double MasslessStoredChange(const double *xv, const int off, double *slope) const override
      {
         if (!unk)
            return 0.0;
         const double dx = xv[off];
         if (slope != nullptr)
            slope[off] += k * (x + dx);
         return k * dx * (x + 0.5 * dx);
      }
      // its unknown while a ball is in its group; the first guess is its row's (MasslessGuess, which the group calls after every part's BeginStep)
      bool BeginStep(const GroupStep::Ball *, const int nb, const double, double *xv, const int off) override
      {
         unk = nb > 0;
         if (unk)
            xv[off] = 0.0;
         return true;
      }
      // the first guess (design section 1 a: the massless face follows the ball): its row r(dx) = k (x + dx) - P(dx) = 0 solved exactly, P the contacts' end push by the full
      // law with the balls' guessed moves held (GroupStep::MasslessRow). Each contact's push falls as the face moves on (dP/ddx <= 0), so r rises (k > 0): one root on
      // dx >= -x, or none there and the face rests at -x (a contact cannot pull: k x1 = N1 >= 0). Safeguarded Newton from the move that keeps the touching contacts' depths
      // (the face's rest with none), inside a bracket: lo, where r < 0, and hi, where r >= 0; a point with r < 0 gives hi = -x + P / k (P falls beyond it, so r(hi) >= 0).
      // A Newton point outside the bracket: the rest when it has not been tried, else the bracket's middle. Done when the move changes by under the solver's tolerance on
      // depth (row.tol) or the bracket is that narrow
      void MasslessGuess(const int, const GroupStep::MasslessRow &row, double *xv, const int off) override
      {
         if (!unk)
            return;
         const double rest = -x;
         double lo = rest, hi = INFINITY, dx = std::max(rest, row.start);
         bool loNeg = false;   // r < 0 at lo (else lo is the rest, not yet tried)
         for (int it = 0; it < 200; it++)   // (bisection alone halves the bracket: 200 is past a double's precision; a guess, so no fault past it)
         {
            double dP;
            const double P = row.Push(dx, dP), r = k * (x + dx) - P;
            if (r == 0.0)
               break;
            if (r > 0.0 && !loNeg && dx <= rest)
            {
               dx = rest;
               break;
            }
            if (r < 0.0)
            {
               lo = dx;
               loNeg = true;
               hi = std::min(hi, rest + P / k);
            }
            else
               hi = dx;
            double nx = dx - r / (k - dP);
            if (!(nx > lo && nx < hi))
               nx = loNeg ? 0.5 * (lo + hi) : rest;
            const bool done = std::fabs(nx - dx) <= row.tol || (loNeg && hi - lo <= row.tol);
            dx = nx;
            if (done)
               break;
         }
         xv[off] = dx;
      }
      // its own share of its row: h k x1 (the group adds -h times the contacts' push at the step's end: the massless row, rule C)
      bool Eval(const double *xv, const int off, const double h, V3 *, V3 *, double *Res, GroupStep::Jac J) override
      {
         if (!unk)
            return true;
         Res[off] = h * k * (x + xv[off]);
         if (J)
            J[off][off] += h * k;
         return true;
      }
      // the ball on the blade: a spring of k on the ball's mass (the face has none)
      void StepFacts(const GroupStep::Ball *b, const int nb, const GroupStep::Settings &, GroupStep::StepFactList &f) override
      {
         double m = INFINITY;
         for (int i = 0; i < nb; i++)
            if (!b[i].held)
               m = std::min(m, b[i].m);
         if (std::isfinite(m) && m > 0.0)
            f.Stiffness(std::sqrt(k / m), false, 1);
      }
      void Accept(const double *xv, const int off, const double h, bool, double &driveJ, double &lossJ, bool &through) override
      {
         driveJ = lossJ = 0.0;
         through = false;
         if (!unk)
            return;
         const double dx = xv[off];
         x += dx;
         rate = h > 0.0 ? dx / h : 0.0;
      }
      void ScaleSpeeds(const double s) override { rate *= s; }
      std::string StartNote() const override
      {
         char s[160];
         snprintf(s, sizeof s, " face x %.17g rate %.17g unknown %d", x, rate, (int)unk);
         return s;
      }
      std::string HistoryNote() const override
      {
         char s[120];
         snprintf(s, sizeof s, " face x %.6g mm blade %.9g uJ", x * 1e3, 0.5 * k * x * x * 1e6);
         return s;
      }
   };

   // what the engine side keeps (PartRecords: part of its state)
   struct EngineState
   {
      bool collidable = true;          // false: the balls pass
      bool idle = false;
      GroupStep::RestWatch rest;       // the rest rule's watch over its last cycle in a group (RestPose)
   };

   // The face as a moving part of the contact system (contactmovingpart.h): the probe, the speed and energy bounds, its Hit event, its rest
   class Engine : public ContactMovingPart, public EngineState
   {
   public:
      Engine() { KeepState<EngineState>(this); }
      Part gp;
      RubberContact::BallSI ball;
      double travelReach = 0.0;        // the farthest the face can travel, for its box and reach: a ball at kMaxSpeed on the blade, v sqrt(m / k) (m)
      static constexpr double kMaxSpeed = 8.0;   // m/s: the fast-ball limit (VPX_SPEED_WATCH's and the driver's default)

      void Build(const std::vector<RubberContact::SurfaceTri> &faces, const std::vector<std::array<int, 3>> &vid, const V3 &travel, const double kBlade, const double edgeRhoM,
                 const RubberContact::BallSI &b)
      {
         gp.Build(faces, vid, travel, kBlade, edgeRhoM);
         ball = b;
         travelReach = kMaxSpeed * std::sqrt(b.m / kBlade);
      }
      void SetCollidable(const bool c) { collidable = c; }
      double Travel() const { return gp.x; }

      // ---- the moving-part interface ----
      void SetBall(const RubberContact::BallSI &b) override { ball = b; }
      void SetGravity(const V3 &) override {}   // (massless: gravity does not load it)
      V3 Origin() const override { return gp.origin; }
      double Reach() const override { return gp.meshR + travelReach; }
      RubberContact::TubeProbe Probe(const V3 &c) const override
      {
         RubberContact::TubeProbe pr;
         if (!collidable)
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
         v = gp.n * gp.rate;
         w = V3();
      }
      // toward a ball: its rate now, plus its whole travel back within the cycle (a massless face unloads at once; it cannot pass x = 0 toward the ball)
      double SurfaceSpeedBound(const double cycleS) const override { return std::fabs(gp.rate) + (cycleS > 0.0 ? gp.x / cycleS : 0.0); }
      // the energy reach: the face travels on [0, sqrt(x^2 + 2 E / k)] (energy enters only through the balls; the blade holds 1/2 k x^2 now, which it can hand back)
      double EnergyMove(const V3 &, const double E, double, double *const work) const override
      {
         const double xMax = std::sqrt(gp.x * gp.x + 2.0 * std::max(E, 0.0) / gp.k);
         if (work != nullptr)
            *work = 0.5 * gp.k * gp.x * gp.x;
         return std::max(gp.x, xMax - gp.x);
      }
      // the one-pass bound (EnergyCap, AGREED item 4b): its energy in the book now (Face::Energy: the blade's 1/2 k x^2; a massless face: no kinetic energy, no weight, no drive), least 0
      double EnergyCap(double, double &, double &) const override { return gp.Energy(nullptr, 0); }
      // the region over a cycle (design part-cull 2R, 2J; piece 2's massless face): the face moves only along n, x over the cycle in [0, xMax], xMax = sqrt(x^2 + 2 E / k)
      // (EnergyMove's bound on the blade spring; no plate in this piece, so no clamp): the mesh's box at rest swept along n over [0, xMax], grown by the edges' radius edgeRho
      void SweptBox(const double E, V3 &lo, V3 &hi) const
      {
         const double xMax = std::sqrt(gp.x * gp.x + 2.0 * std::max(E, 0.0) / gp.k);
         const V3 d = gp.n * xMax, g(gp.edgeRho, gp.edgeRho, gp.edgeRho);
         lo = V3(std::min(gp.boxLo.x, gp.boxLo.x + d.x), std::min(gp.boxLo.y, gp.boxLo.y + d.y), std::min(gp.boxLo.z, gp.boxLo.z + d.z)) - g;
         hi = V3(std::max(gp.boxHi.x, gp.boxHi.x + d.x), std::max(gp.boxHi.y, gp.boxHi.y + d.y), std::max(gp.boxHi.z, gp.boxHi.z + d.z)) + g;
      }
      static double BoxDist(const V3 &lo, const V3 &hi, const V3 &c)   // the point-to-box distance (0 inside)
      {
         const double dx = std::max({ lo.x - c.x, 0.0, c.x - hi.x }), dy = std::max({ lo.y - c.y, 0.0, c.y - hi.y }), dz = std::max({ lo.z - c.z, 0.0, c.z - hi.z });
         return std::sqrt(dx * dx + dy * dy + dz * dz);
      }
      // the grouping (E: the group's feed; it takes the larger of this less R and the gap less EnergyMove, so the swept box rules wherever it is the tighter) and the cull (E: the
      // table's, precomputed once per cycle) read the same region
      double CycleRegionDist(const V3 &c, double, const double E) const override
      {
         V3 lo, hi;
         SweptBox(E, lo, hi);
         return BoxDist(lo, hi, c);
      }
      V3 cullLo, cullHi;
      void CullCycle(const double E, const double T) override
      {
         cullE = E;
         cullT = T;
         SweptBox(E, cullLo, cullHi);
      }
      double CullDist(const V3 &c) const override { return BoxDist(cullLo, cullHi, c); }
      std::string StateNote() const override
      {
         char buf[160];
         std::snprintf(buf, sizeof buf, "face travel %.4f mm at %.4f m/s, blade %.0f N/m, asleep %d", gp.x * 1e3, gp.rate, gp.k, (int)idle);
         return buf;
      }

      // ---- the group's calls ----
      GroupStep::Part *GroupPart() override { return &gp; }
      int GroupSurfaces() const override { return 1; }
      void GroupProbe(int, const V3 &c, const double R, const double *const trial, GroupStep::Hit &h) const override
      {
         gp.ProbeAt(c, R, trial != nullptr ? trial[0] : 0.0, h);
         if (!collidable)
         {
            h.delta = -1.0;
            h.npat = 0;   // (no depth and no patches: no friction group either)
            h.npc = 0;
         }
      }
      bool GroupRestCentre(int, const V3 &c, const double *const move, V3 &out) const override   // (SlideMesh: the probe's own centre, c - (x + dx) n)
      {
         out = c - gp.mesh.a * (gp.x + (move != nullptr ? move[0] : 0.0));
         return true;
      }
      void GroupLinkWindows(int, const V3 &c0r, const V3 &c1r, const double R, std::vector<RubberContact::LinkWindow> &out) const override { gp.mesh.LinkWindows(c0r, c1r, R, out); }
      bool GroupPatchZero(int) const override { return gp.mesh.PatchZeroOn(); }   // (the per-patch zero: keyed by face id, P3b piece 4)
      void CycleBegin() override
      {
         idle = false;
         gp.unk = true;
      }
      void TouchBegan(int, const int ballIndex, int, const double vn, const V3 &) override { events.push_back({ Event::Hit, ballIndex, vn }); }
      bool Sleeping() const override { return idle; }
      bool RestNumbers(double &moveM, double &speedMS) const override
      {
         moveM = rest.CycleMoved(RestPose());
         speedMS = std::fabs(gp.rate);
         return true;
      }
      bool AtRest() const override { return idle || rest.Rests(RestPose(), std::fabs(gp.rate)); }
      GroupStep::RestWatch *Rest() override { return &rest; }
      GroupStep::RestWatch::Pose RestPose() const override { return GroupStep::RestWatch::Pose().Add(gp.x, 1.0); }
      // No ball in its group: the massless face sits at its free length with no unknown (design section 1); nothing to step
      int AdvanceAlone(const double dt) override
      {
         gp.unk = false;
         if (dt <= 0.0 || idle)
            return 0;
         gp.x = 0.0;
         gp.rate = 0.0;
         idle = rest.Step(RestPose(), 0.0, false, restTouched, dt);
         return 0;
      }
   };
}
