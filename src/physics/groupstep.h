// license:GPLv3+
//
// The group step (docs/GroupStepDesign.txt; ProjectAssessment-2026-10-27.txt item 0 b): the balls of a touching group and their contacts moved together in one solve,
// so that energy is kept by construction. Balls, fixed surfaces, ball pairs and moving parts.
//
// Laws (the shared contact law restated as stored energy plus losses, design pieces 4 and 5):
//  - normal: F = kH [delta^1.5 + beta reff^1.5 delta^3 / (h - delta)^3] (Hertz plus the thin bonded layer at the thickness left; no layer for steel), stored energy
//    V = kH [0.4 delta^2.5 + beta reff^1.5 W(delta)]. Its losses act on the mid-step compression rate: Hunt-Crossley 1.5 alpha F with F over the step (or Kuwabara-Kono gamma a), plus the
//    flat loss tanDelta sqrt(1.5 F m / delta) where normalLoss is on.
//  - shear: Mindlin's spring kt = 8 G* a shearScale as stored energy 1/2 kt |xi|^2 (xi: the spring's stretch in the contact plane, carried from step to step and turned
//    with the contact, its length kept); a slip past mu N is a loss; the spring's damper shearTan sqrt(kt mt) acts on the mid-step slip rate, spring and damper together
//    within mu N. Mindlin-Deresiewicz: a patch that grows holds the force, one that shrinks holds the stretch (rubbercontact.h ContactCore); what the spring holds only
//    falls by that, a loss, and no force acts along the compression for it.
//  - twist and rolling resistance: today's smooth laws (rubbercontact.h ContactCore, EvaluateTube) on the mid-step spin: losses.
// Step: each ball's move and spin change over the step are the unknowns (6 per ball: the shear spring ties travel and spin). Each stored energy's force is its discrete
// gradient over the move (Gonzalez: exact differences), every loss acts on mid-step rates, and the shear spring's force is k1 xibar over the step, at the stiffness of the
// step's end. Once the step's equations are solved: energy end = start + drive work - losses, at any step size. The step size is set by accuracy only (ChooseStep). A contact's effective radius, layer and material are held for the life of a touch (taken when it starts),
// so its energy is one function throughout.
// Moving parts (step 1b, GroupStepDesign.txt 10, 18): a Part describes its own bodies' coordinates and unknowns (after the balls'), its stored energy, losses, drives, grip and the facts of its step, and the
// group solves it with the balls (the band cord, cordfixed.h, is the first). A surface it carries says how a point on it moves with each coordinate (Hit's columns); the group does the rest.
// Pure SI maths, no VPX types, no allocation in the step: the offline tests runs/2026-10-29-group-step and runs/2026-10-29-group-band build it.

#pragma once

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "physics/rubbercontact.h"
#include "physics/pairform.h"

namespace GroupStep
{
   using V3 = RubberContact::V3;
   using RubberContact::RubberParams;
   constexpr double kPi = 3.14159265358979323846;

   // Storage of a group: grown between cycles (Group::Size, AddShape, AddPart, Configure), never inside a step; an index or a plain pointer to it, as an array
   template <class T> class Buf
   {
   public:
      void Grow(const size_t n)
      {
         if (v.size() < n)
            v.resize(n);
      }
      operator T *() { return v.data(); }
      operator const T *() const { return v.data(); }
      template <class I> T &operator[](const I i) { assert((size_t)i < v.size()); return v.data()[i]; }
      template <class I> const T &operator[](const I i) const { assert((size_t)i < v.size()); return v.data()[i]; }

   private:
      std::vector<T> v;
   };
   // A list only the group adds to (AddShape, AddPart: they size it); to anyone else it reads
   class Group;
   template <class T> class AddedList
   {
   public:
      template <class I> const T &operator[](const I i) const { return b[i]; }

   private:
      friend class Group;
      Buf<T> b;
   };
   // A square matrix by rows, J[r][c]; empty = none wanted
   class Jac
   {
   public:
      Jac() = default;
      Jac(std::nullptr_t) {}
      Jac(double *const data, const int stride) : d(data), n(stride) {}
      double *operator[](const int r) const { return d + r * n; }
      explicit operator bool() const { return d != nullptr; }

   private:
      double *d = nullptr;
      int n = 0;
   };

   // ---- 3 x 3 helpers ----
   struct M3
   {
      double a[3][3] = {};
      static M3 Outer(const V3 &u, const V3 &v)
      {
         M3 m;
         const double uu[3] = { u.x, u.y, u.z }, vv[3] = { v.x, v.y, v.z };
         for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
               m.a[i][j] = uu[i] * vv[j];
         return m;
      }
      static M3 Skew(const V3 &r)   // Skew(r) v = r x v
      {
         M3 m;
         m.a[0][1] = -r.z; m.a[0][2] = r.y;
         m.a[1][0] = r.z;  m.a[1][2] = -r.x;
         m.a[2][0] = -r.y; m.a[2][1] = r.x;
         return m;
      }
      static M3 Tangent(const V3 &n)   // I - n n^T
      {
         M3 m = Outer(n, n);
         for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
               m.a[i][j] = (i == j ? 1.0 : 0.0) - m.a[i][j];
         return m;
      }
      M3 operator*(const M3 &o) const
      {
         M3 m;
         for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
               m.a[i][j] = a[i][0] * o.a[0][j] + a[i][1] * o.a[1][j] + a[i][2] * o.a[2][j];
         return m;
      }
      M3 T() const
      {
         M3 m;
         for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
               m.a[i][j] = a[j][i];
         return m;
      }
   };

   // A position energy on the balls (a magnet's pull): what a ball's place stores, per unit of its Ball::fieldScale. Fixed for the cycle; it owns no state of a step. grad, hess: its slopes in the centre, when wanted
   class BallField
   {
   public:
      virtual ~BallField() = default;
      virtual double Eval(const V3 &c, double R, V3 *grad, M3 *hess) const = 0;
   };

   // ---- surfaces ----
   // What a surface reports for a ball centre: the signed depth (> 0 = squashed; the ball's radius plus a tube's radius minus the distance), the normal (surface -> ball
   // centre) and the Hertz effective radius. The depth must be continuous in the ball's position (the step takes exact differences of it).
   // How a surface point moves with one coordinate of the part that carries it. Plain numbers with no starting values, so that a fixed surface's hit costs nothing.
   struct Column
   {
      int coord;           // the part's coordinate number, the same whether or not it has an unknown now
      double slope;        // d(depth) / d(q); the group multiplies it by Hit::grad
      double v[3], w[3];   // the point's speed and the surface's spin per unit rate of q
      void Set(const int c, const double s, const V3 &vq, const V3 &wq)
      {
         coord = c;
         slope = s;
         v[0] = vq.x; v[1] = vq.y; v[2] = vq.z;
         w[0] = wq.x; w[1] = wq.y; w[2] = wq.z;
      }
      V3 V() const { return V3(v[0], v[1], v[2]); }
      V3 W() const { return V3(w[0], w[1], w[2]); }
   };
   struct Hit
   {
      double delta = -1.0;
      V3 n;
      double reff = 0.0;
      double curvDist = 0.0;   // > 0: the ball turns round a centre this far away while it rolls along the surface (a post, a tube): the step's turning rule
      double grad = 1.0;       // the depth's gradient is n x grad (below 1 where the contact patch spans a fold of the surface)
      double kb = 0.0;         // > 0: the surface bends in series with the contact here (N/m along the normal): a massless spring, taken with the law when a touch starts and held while it lasts
      int nc = 0;              // the columns in use; 0 = a fixed surface
      Column c[7];             // one rigid body's 6 coordinates and 1 given motion of its frame: a physical size, not a table limit
      Hit() {}
      Hit(const Hit &o) { *this = o; }
      Hit &operator=(const Hit &o)
      {
         delta = o.delta;
         n = o.n;
         reff = o.reff;
         curvDist = o.curvDist;
         grad = o.grad;
         kb = o.kb;
         nc = o.nc;
         for (int j = 0; j < nc; j++)
            c[j] = o.c[j];
         return *this;
      }
   };
   class Shape
   {
   public:
      RubberParams rp;         // the material (rubbercontact.h)
      double layerM = 0.0;     // the rubber's thickness for the layer term, m (0 = no layer term)
      int onlyBall = -1;       // >= 0: only this ball meets it (a surface's list of near faces is one ball's)
      int part = -1;           // the part that carries it, its place among the group's parts; -1 = fixed. Its Probe fills the same columns, in the same order, every time
      mutable const double *trial = nullptr;   // a carried surface: each of its part's coordinates moved by this for the probe; null = where they are now
      virtual ~Shape() = default;
      virtual void Probe(const V3 &c, double R, Hit &h) const = 0;
   };
   // A flat face through p with unit normal n (the ball's side); unbounded
   class PlaneShape : public Shape
   {
   public:
      V3 p, n;
      void Probe(const V3 &c, double R, Hit &h) const override
      {
         h.delta = R - (c - p).Dot(n);
         h.n = n;
         h.reff = R;
         h.curvDist = 0.0;
      }
   };
   // A tube of radius rho round a polyline centreline (a post's rubber ring, a wall rubber): the engine's own probe (rubbercontact.h ProbeTube), with the centreline's
   // precomputed bounds when given (built once by the owner: BuildTubePre over the nSeg + 1 points)
   class TubeShape : public Shape
   {
   public:
      static constexpr int kMaxPts = 129;
      V3 pts[kMaxPts];
      int nSeg = 0;                                   // segments: nSeg + 1 points
      double rho = 0.0;
      const RubberContact::TubePre *pre = nullptr;
      void Probe(const V3 &c, double R, Hit &h) const override
      {
         RubberContact::BallSI b;
         b.R = R;
         const RubberContact::TubeProbe pr = RubberContact::ProbeTube(b, pts, nSeg + 1, rho, c, pre);
         h.delta = pr.delta;
         h.n = pr.normal;
         h.reff = pr.reff;
         h.curvDist = pr.dist;
      }
   };

   struct Ball
   {
      V3 c, v, w;            // centre, velocity, spin (m, m/s, rad/s)
      double m = 0.080, R = 25.0 * 0.53975e-3;
      V3 fExt;               // a held outside force (a magnet not on a field, the magnets' ball-pair force, a part not yet moved over): a drive, its work counted
      double fieldScale = 0.0;   // how much of the group's fields it feels (Group::field); 0 = none
      double brakeRate = 0.0;    // a viscous brake on its spin across an axis (the eddy currents in a coil's field), 1/s; 0 = none. The caller's, held for the cycle
      V3 brakeAxis;              // the axis (unit): a spin about it is not braked
      double I() const { return 0.4 * m * R * R; }
   };

   struct Settings
   {
      V3 g = V3(0.0, 0.0, -9.81);
      RubberParams pair;     // ball on ball (steel)
      bool ballPairs = true;
      // the step (ChooseStep): accuracy only
      double omegaH = 0.6;                     // h <= omegaH / omega of the stiffest contact (normal, and shear when shearInStep)
      bool shearInStep = true;                 // the shear spring's stiffness counts in the step choice (false: the normal stiffness only; the solve is stable either way)
      bool slideFree = true;                   // a slide at the friction limit counts as steady (docs/GroupStepDesign.txt 16: its force is mu N, the spring's stiffness is not in the motion)
      // The step rule by motion (docs/GroupStepDesign.txt 9, 16, 17): a contact that is steady only needs its normal stiffness resolved to omegaHRest; an impact, a slip the spring can still catch or a fast
      // arrival keeps the strict rule above. Steady is motion, not load: its closing rate is within vSteady, the net force on its bodies cannot change that rate by more in the contact's own time (at a hit's
      // turnaround the rate is zero and the net force is the whole contact force; under a pile the net force is near zero whatever the load), and its slip is within vSteady or at the friction limit.
      // vSteady 0 = every contact strict. 0.02 m/s and 2.4 from 750 real hits against their 1 us runs (runs/2026-10-31-static-group).
      double vSteady = 0.02;                   // m/s
      double omegaHRest = 2.4;
      double turnMaxRad = 2.0 * kPi / 180.0;   // a ball rolling round a post or another ball turns at most this per step
      double tauShare = 0.25, fieldShare = 0.25;   // a part's own time constant is resolved to this share of it; a body moves at most this share of a distance over which its force changes (StepFactList; no part uses them yet; the bat's values)
      double vFloor = 0.01;                    // m/s: the closing speed below which a contact's expected depth is taken at this speed
      double hMax = 1e-3, hMin = 1e-9;
      // the solver
      double rollSmooth = 0.05;                // the rolling torque's smoothing through zero spin, rad/s ([Contact] roll_smooth_rad_s)
      int rollModel = 1;                       // 1 (default) = rolling and twist stick springs at the contact instead of the smoothed torques ([Contact] roll_model)
      double tolV = 1e-11;                     // solved when a round changes no ball's speed (or its rim speed R w) by more than this, m/s
      double tolNoise = 1e-9;                  // or when the changes stay below this for stallRounds rounds without halving: the forces' own noise (m/s)
      double tolEndRes = 0.0;                  // > 0: a step counts as solved only with its residual within this (m/s), and a noise stop keeps the point it read
      double tolRes = 1e-10, tolCheckMove = 1e-6;   // the residual in speed units that counts as solved; the move after which a round looks at the residual first, with no factorization (m/s)
      int stallRounds = 6;
      int maxRounds = 60;
      int noProgressRounds = 12;               // an attempt gives up (the step is halved, as at maxRounds) when this many rounds have not halved its best move
      double maxMoveM = 2e-3, maxSpinStep = 500.0;   // a round moves a ball at most this far, and changes its spin by at most this (rad/s)
      // diagnostics: 1 = also form the Jacobian by central differences of the residual each round and keep both (jacAna, jacFd of the last round); 2 = solve with
      // the difference Jacobian. 0 = off (normal).
      int jacMode = 0;
   };

   // Energy books (losses positive) and counters
   struct Books
   {
      double driveJ = 0.0, normalLossJ = 0.0, slipJ = 0.0, shearDampJ = 0.0, twistJ = 0.0, rollJ = 0.0;
      double partLossJ = 0.0;   // the parts' losses (their own books split them)
      double spinBrakeJ = 0.0;  // the balls' spin brake (Ball::brakeRate)
      double coilJ = 0.0;       // a coil's work on what the balls store: a current that changed since the last cycle, at the balls' places. It falls before this cycle's energy starts, so it is apart from driveJ
      long steps = 0, rounds = 0, halvings = 0, stopSteps = 0, noiseStops = 0;   // noiseStops: steps solved to the forces' noise (Settings::tolNoise)
      long energySteps = 0, energyBounces = 0, energyFree = 0;   // steps no solve could do: taken on energy alone, of them turned back at a step in the energy too high to pay, moved with no path (Group::EnergyStep)
      int maxRounds = 0;
      double hMinUsed = 1.0, maxDepthM = 0.0;
      bool failed = false, through = false;
      double Losses() const { return normalLossJ + slipJ + shearDampJ + twistJ + rollJ + partLossJ + spinBrakeJ; }
   };

   // VPX_GROUP_HFIX=<us> (a diagnostic, read once): every step of every group, balls alone, with parts or a part alone, is this long, not the step rule's (the energy's step sweep, and the
   // fixed-step reference of a behaviour check, e.g. 1). Its cycles may take more than the step guard's 3000 steps.
   inline double GroupFixedStepS()
   {
      static const double s = []() { const char *const e = getenv("VPX_GROUP_HFIX"); return e != nullptr && *e ? atof(e) * 1e-6 : 0.0; }();
      return s;
   }

   class Group;

   // What a coordinate of a part is, for the group. A coordinate may have an unknown or not (a body on a stop, stuck or held has none); a surface's column names it either way.
   struct CoordInfo
   {
      double rate = 0.0;       // its rate now
      double given = 0.0;      // with no unknown: its move over the step, from BeginStep on (0 = at rest; a driven stop's move). The work this does on the balls is booked as a drive
      double invM = 0.0;       // 1 / its mass or inertia while it has an unknown
      double invMFree = 0.0;   // the same while it has none: what a ball meets when it pushes the body off (0 = no such body)
      double reach = 0.0;      // the distance from which rate x reach is a speed, m: the solver's units, the far test
      double maxMove = 0.0;    // a Newton round moves it at most this; 0 = no limit
      double settle = 0.0;     // 0: its Newton move counts with the balls', as a speed at its reach. > 0: it stays out of theirs, and is solved when its own rate changes by no more than this (the cord's plate)
      int blockedSide = 0;     // with no unknown: it rests on a stop blocking this side (+1: larger q); 0 = held by friction only
   };

   // What a part says about its step: facts in its own numbers. The factors are the group's (Settings), in this one place. Kept as the shortest so far, so nothing is stored; why: the part's own number, for the trace.
   class StepFactList
   {
   public:
      explicit StepFactList(const Settings &s) : set(s), h(s.hMax) {}
      void Stiffness(const double omega, const bool steady, const int why)   // a stiffness in the motion: its natural frequency; steady: the motion it carries is at rest
      {
         if (omega > 0.0)
            Lim((steady ? set.omegaHRest : set.omegaH) / omega, why);
      }
      void Approach(const double gap, const double rate, const double omega, const int why) { Arrival(gap, rate, set.omegaH / omega, why); }   // closing on a stiffness of this frequency: reached in about two steps, not overshot
      void Arrival(const double gap, const double rate, const double shortest, const int why) { Lim(std::max(0.5 * gap / rate, shortest), why); }   // closing on something at this rate; no shorter than this
      void Turn(const double radius, const double speed, const int why) { Lim(set.turnMaxRad * radius / speed, why); }   // turning round a centre this far away at this speed
      void TimeLeft(const double t, const int why) { Lim(t, why); }   // something happens in this time
      void TimeConstant(const double tau, const int why) { Lim(set.tauShare * tau, why); }   // a state of its own relaxes with this time constant
      void Field(const double distance, const double speed, const double accel, const int why)   // a force that changes over this distance, on a body at this speed and acceleration
      {
         const double d = set.fieldShare * distance;
         if (accel > 1e-9)
            Lim((std::sqrt(speed * speed + 2.0 * accel * d) - speed) / accel, why);
         else if (speed > 1e-9)
            Lim(d / speed, why);
      }
      double H() const { return h; }
      int Rule() const { return rule; }

   private:
      void Lim(const double c, const int why)
      {
         if (c < h)
         {
            h = c;
            rule = why;
         }
      }
      const Settings &set;
      double h;
      int rule = 0;
   };

   // A moving part of a group (GroupStepDesign.txt piece 1; 10.1; 18): it describes its bodies, stored energy, losses, drives and grip, and the group moves them with
   // the balls in one solve. The unknowns x: per ball the move (3) and the spin change (3), then each part's own (NumDof of them, from off). The balls' residual
   // rows are 2 m dc / h - 2 m v - h F and I dw - h M; a part adds its forces on the balls (averaged over the step) to F and M and writes its own rows.
   // Every call of a step reads the balls' state at the step's start (BeginStep's b) and the move x.
   class Part
   {
   public:
      virtual ~Part() = default;
      // what it is
      virtual int NumDof() const = 0;
      virtual int MaxDof() const { return NumDof(); }   // the most its NumDof can be inside a cycle (the group's storage is sized between cycles)
      virtual int NumCoords() const = 0;                // constant between cycles
      virtual int DofOf(int coord) const = 0;           // the coordinate's place among its unknowns now; -1 = it has none now
      virtual CoordInfo Coord(int coord) const = 0;
      // the part's stored energy and its own bodies' kinetic energy at the balls' positions b (the balls' own energy is the group's)
      virtual double Energy(const Ball *b, int nb) const = 0;
      // a step of h from the state b: its start values, and the first guess of its own unknowns into x[off..]
      virtual bool BeginStep(const Ball *b, int nb, double h, double *x, int off) = 0;
      // the step's forces at x (see above); with J its entries of the residual's Jacobian (the balls' rows and columns in the group's layout). The group adds the balls' reaction on its own rows
      virtual bool Eval(const double *x, int off, double h, V3 *F, V3 *M, double *Res, Jac J) = 0;
      // between the solve's rounds: its grip impulses on the balls from the rates at x; returns their largest change (m/s)
      virtual double UpdateGrip(const double *, double) { return 0.0; }
      virtual bool InRange(const double *, int) const { return true; }   // the move x may be tried (no layer squashed to nothing)
      // the step's length: facts, from the state now (accuracy only; GroupStepDesign 10.3)
      virtual void StepFacts(const Ball *b, int nb, const Settings &s, StepFactList &f) = 0;
      // hard things: the events of a step, run by Group::Step
      virtual double StopExcess(const double *, int) const { return -INFINITY; }   // > 0: the move x carries a body past a hard stop by this much (the step is shortened to end on it; piece 1 item 7)
      virtual bool ArrivingAtStop() const { return true; }               // the move last tested by StopExcess is closing on what it came near
      virtual bool PinOnStop(const double *, int) { return false; }      // on its hard stop already and the move x just solved would take it on past: it stays there. False: nothing changed
      virtual bool FreedByStep() { return false; }                       // a body with no unknown: the forces of the step just solved take it off, so it has one and the step is solved again
      // the balls' reaction on it in this evaluation, before its Eval: per coordinate (one with no unknown too), and the force on the whole part (a pivot's friction goes with its load). Only a part carrying a surface hears it
      virtual void Loads(const double *, const V3 &) {}
      // the step taken (after an Eval at the solution, before the group moves its balls); onStop: it ended on the hard stop. Its drive work and losses of the
      // step, and whether a ball went through it
      virtual void Accept(const double *x, int off, double h, bool onStop, double &driveJ, double &lossJ, bool &through) = 0;
      virtual void AfterStep(const Group &) {}   // the step taken and the balls moved: what the part keeps of its own shapes' contacts
      // an energy step (Group::EnergyStep): the state kept and put back; its bodies' kinetic energy and all their speeds times s; what a step's forces left for the next one is put back
      virtual void SaveState() {}
      virtual void RestoreState() {}
      virtual double KineticJ() const { return 0.0; }
      virtual void ScaleSpeeds(double) {}
      virtual void AfterEnergyStep() {}
      virtual void HoldArcs() {}   // no solve in this attempt: what the part solves within itself stays as it is
      // diagnostics (rule 8): an attempt of a step failed and was halved (slot: its place in Group::halveDiag); the part keeps what its last round evaluated, Res
      // the residual of that round (its rows from off). Told only when an energy step is reported
      virtual void KeepFailNote(int, const double *, int) {}
      virtual std::string FailNote(int) const { return std::string(); }
      virtual void BallsMoved(const Ball *, int) {}              // the balls were moved with no step of its own
      virtual void ExplicitMove(const double *, int, double, double *) const {}   // its own unknowns in an energy step's explicit round
      virtual double EndResidual() const { return 0.0; }         // > 0: a step of a group holding this part counts as solved only with its residual within this (m/s; Settings::tolEndRes)
   };

   // One contact: ball a on shape s, or ball a on ball b. Its law is taken when a touch starts and held while it lasts.
   struct Item
   {
      int a = -1, b = -1, s = -1;
      bool touching = false;
      V3 xi;                          // the shear spring's stretch (in the contact plane)
      RubberContact::RollAngle ra;    // the rolling and twist springs' angles (roll_model 1)
      double dZero = 0.0;             // the depth its squash counts from: a ball found inside has stored nothing (no energy without the work that stored it; only ever down)
      // the held law
      double kH = 0.0, cL = 0.0, hL = 0.0, reff = 0.0, ktPerA = 0.0, E = 0.0, mRed = 0.0, mt = 0.0;
      double alpha = 0.0, gammaKK = 0.0, tanDelta = 0.0, shearTan = 0.0, mu = 0.0, rollingMu = 0.0;
      bool normalLoss = false, twist = false;
      double kb = 0.0;                // the surface's bend in series with the contact (0 = none): V(d) = min over u of Hertz(u) + 1/2 kb (d - u)^2, kH u^1.5 = kb (d - u)
      // the contact's share of the depth
      double Uc(const double d) const
      {
         if (kb <= 0.0 || d <= 0.0)
            return d;
         double u = std::min(d, std::pow(kb * d / kH, 2.0 / 3.0));   // above the root of a convex function: Newton comes down to it
         for (int i = 0; i < 30; i++)
         {
            const double su = std::sqrt(std::max(u, 0.0));
            const double du = (kH * u * su + kb * (u - d)) / (1.5 * kH * su + kb);
            u = std::max(0.0, u - du);
            if (std::fabs(du) <= 1e-14 * d)
               break;
         }
         return u;
      }
      double W(double d) const   // integral of x^3 / (h - x)^3 from 0 to d; for a shallow squash its series (the closed form cancels to rounding there)
      {
         const double s = d / hL;
         if (s < 0.01)
         {
            const double s4 = s * s * s * s;
            return hL * s4 * (0.25 + s * (0.6 + s * (1.0 + s * (10.0 / 7.0 + s * (15.0 / 8.0)))));
         }
         const double u = hL - d;
         return 1.5 * hL + 3.0 * hL * std::log(hL / u) + hL * hL * hL / (2.0 * u * u) - 3.0 * hL * hL / u + u;
      }
      double V(double d) const
      {
         if (kb > 0.0 && d > 0.0)
         {
            const double u = Uc(d);
            return 0.4 * kH * u * u * std::sqrt(u) + 0.5 * kb * (d - u) * (d - u);
         }
         return d > 0.0 ? 0.4 * kH * d * d * std::sqrt(d) + (cL > 0.0 ? cL * W(d) : 0.0) : 0.0;
      }
      double F(double d) const
      {
         if (d <= 0.0)
            return 0.0;
         if (kb > 0.0)
         {
            const double u = Uc(d);
            return kH * u * std::sqrt(u);
         }
         double f = kH * d * std::sqrt(d);
         if (cL > 0.0)
         {
            const double x = d / (hL - d);
            f += cL * x * x * x;
         }
         return f;
      }
      double K(double d) const
      {
         if (d <= 0.0)
            return 0.0;
         if (kb > 0.0)
         {
            const double kc = 1.5 * kH * std::sqrt(Uc(d));
            return kc * kb / (kc + kb);
         }
         double k = 1.5 * kH * std::sqrt(d);
         if (cL > 0.0)
         {
            const double u = hL - d, x = d / u;
            k += cL * 3.0 * x * x * (1.0 / u + d / (u * u));
         }
         return k;
      }
      // (V(d1) - V(d0)) / (d1 - d0), without the cancellation of the Hertz part; the midpoint force when the change is mostly rounding
      double Chord(double d0, double d1) const
      {
         if (d0 <= 0.0 && d1 <= 0.0)
            return 0.0;
         if (std::fabs(d1 - d0) < 1e-13)
            return F(0.5 * (d0 + d1));
         if (kb > 0.0)
         {
            if (d0 > 0.0 && d1 > 0.0 && std::fabs(d1 - d0) < 0.01 * std::max(d0, d1))   // a short move: the integral by quadrature, no cancellation (as the layer's, below)
            {
               const double dm = 0.5 * (d0 + d1), g = 0.3872983346207417 * (d1 - d0);
               return (5.0 * F(dm - g) + 8.0 * F(dm) + 5.0 * F(dm + g)) / 18.0;
            }
            return (V(d1) - V(d0)) / (d1 - d0);
         }
         double hz;
         if (d0 > 0.0 && d1 > 0.0)
         {
            const double s0 = std::sqrt(d0), s1 = std::sqrt(d1), s0q = d0, s1q = d1;
            hz = 0.4 * kH * (s1q * s1q + s1q * s1 * s0 + s1q * s0q + s1 * s0q * s0 + s0q * s0q) / (s1 + s0);
         }
         else
         {
            const double a = std::max(d0, 0.0), b = std::max(d1, 0.0);
            hz = 0.4 * kH * (b * b * std::sqrt(b) - a * a * std::sqrt(a)) / (d1 - d0);
         }
         double lw = 0.0;
         if (cL > 0.0)
         {
            if (d0 > 0.0 && d1 > 0.0 && std::fabs(d1 - d0) < 0.01 * (hL - std::max(d0, d1)))
            {
               const double dm = 0.5 * (d0 + d1), g = 0.3872983346207417 * (d1 - d0);
               const double q0 = (dm - g) / (hL - dm + g), q1 = dm / (hL - dm), q2 = (dm + g) / (hL - dm - g);
               lw = cL * (5.0 * q0 * q0 * q0 + 8.0 * q1 * q1 * q1 + 5.0 * q2 * q2 * q2) / 18.0;
            }
            else
               lw = cL * ((d1 > 0.0 ? W(d1) : 0.0) - (d0 > 0.0 ? W(d0) : 0.0)) / (d1 - d0);
         }
         return hz + lw;
      }
      double A(double d) const { return d > 0.0 ? std::sqrt(reff * Uc(d)) : 0.0; }   // contact patch radius
      double Kt(double d) const { return ktPerA * A(d); }
      // the rolling and twist springs of the patch (roll_model 1): tilting (4/3) E* a^3 and torsion (16/3) G* a^3 of a Hertz contact, G* = E* (1 - nu) / 2 with nu 0.3
      // (E* = 3 kH / (4 sqrt(reff)): the normal law's own contact modulus, a pair's compliances added)
      double Kr(double d) const { const double a = A(d); return kH / std::sqrt(reff) * a * a * a; }
      double Ktw(double d) const { const double a = A(d); return 16.0 / 3.0 * 0.35 * 0.75 * kH / std::sqrt(reff) * a * a * a; }
      double KtSlope(double d) const
      {
         if (d <= 0.0)
            return 0.0;
         if (kb > 0.0)
         {
            const double u = Uc(d), kc = 1.5 * kH * std::sqrt(u);
            return u > 0.0 ? 0.5 * ktPerA * std::sqrt(reff / u) * kb / (kc + kb) : 0.0;
         }
         return 0.5 * ktPerA * std::sqrt(reff / d);
      }
      double KKIntegral(const double d) const   // 2/3 sqrt(reff) d^1.5, 0 for no touch
      {
         if (d <= 0.0)
            return 0.0;
         const double u = Uc(d);
         return 2.0 / 3.0 * std::sqrt(reff) * u * std::sqrt(u);
      }
      // Viscoelastic damping over the step: the law's own integral at a steady rate, as Hunt-Crossley already is (runs/2026-10-31-gate-spinner-disc/patch_gs_rotor.py, without its bend in series). F = gamma a(d) d_dot, so the mean
      // coefficient over a move is gamma (G(d1) - G(d0)) / (d1 - d0). At the mid-step depth alone it rises as sqrt(d) from a touch's start: a slope with no bound, where Newton stalls. jacTerm: its slope for the Jacobian (d(c x (d1 - d0)) / d d1)
      double KKStep(const double d0, const double d1, const double dm, double &jacTerm) const
      {
         const double top = std::max(std::fabs(d0), std::fabs(d1));
         jacTerm = gammaKK * (d1 > 0.0 ? A(d1) * (kb > 0.0 ? kb / (kb + 1.5 * kH * std::sqrt(Uc(d1))) : 1.0) : 0.0);
         if (std::fabs(d1 - d0) > 1e-3 * top && top > 0.0)
            return gammaKK * (KKIntegral(d1) - KKIntegral(d0)) / (d1 - d0);
         double c = gammaKK * A(dm);
         if (kb > 0.0 && dm > 0.0)
            c *= kb / (kb + 1.5 * kH * std::sqrt(Uc(dm)));
         return c;
      }
      // the normal damping coefficient over a step, mid-step depth d and mean spring force fMean: force = c x compression rate. Hunt-Crossley's 1.5 alpha F is taken over the step, so its loss is
      // 1.5 alpha x the change in stored energy x the rate, the law's own integral at a steady rate; at the mid-step force alone the loss falls short of it by (omega h)^2 and the rebound comes out high
      double DampC(const double d, const double fMean) const
      {
         double c = gammaKK > 0.0 ? (d > 0.0 ? gammaKK * A(d) : 0.0) : 1.5 * alpha * fMean;
         if (kb > 0.0 && d > 0.0)   // the contact takes its share of the approach rate
            c *= kb / (kb + 1.5 * kH * std::sqrt(Uc(d)));
         if (d > 0.0 && normalLoss && tanDelta > 0.0)
            c += tanDelta * std::sqrt(1.5 * F(d) * mRed / d);
         return c;
      }
   };

   class Group
   {
   public:
      Settings set;
      Books books;
      int nb = 0, ns = 0, ni = 0, np = 0, nf = 0;
      Buf<Ball> ball;
      AddedList<const BallField *> field;   // the position energies on the balls: one energy, each ball's difference over its move exact
      // the magnets' ball-pair form (pairform.h, 2026-11-03): a stored energy E = -sum g q of the balls' places, fixed for the cycle; a and b index the group's balls
      // (b < 0: the other ball held at the term's o, outside this group). Its work over a move is the exact difference (discrete gradient), as the fields'
      Buf<PairForm::Term> form;
      int nform = 0;
      AddedList<const Shape *> shape;
      Buf<Item> item;
      AddedList<Part *> part;   // moving parts (the band); they read the balls through the step calls

      Group() { Size(8); }
      // Room for this many balls (before the balls are written, between cycles)
      void Size(const int balls)
      {
         const size_t n = (size_t)balls;
         ball.Grow(n);
         b1r.Grow(n);
         bmr.Grow(n);
         aLast.Grow(n);
         alphaLast.Grow(n);
         fNetS.Grow(n);
         evF.Grow(n);
         evM.Grow(n);
         evDc.Grow(n);
         evDw.Grow(n);
         evWm.Grow(n);
         evB1.Grow(n);
         evBm.Grow(n);
         farB1.Grow(n);
         farBm.Grow(n);
         esB0.Grow(n);
         esA0.Grow(n);
         esAl0.Grow(n);
         fieldE0.Grow(n);
         warmA.Grow(n);
         warmAl.Grow(n);
      }
      void AddField(const BallField *const f)
      {
         field.b.Grow((size_t)nf + 1);
         field.b[nf++] = f;
      }
      void AddForm(const PairForm::Term &t)
      {
         form.Grow((size_t)nform + 1);
         form[nform++] = t;
      }
      // What the balls' places store in the pair form; each ball's force from it (-slope) added to F when wanted
      double FormEnergy(const Ball *const b, V3 *const F = nullptr) const
      {
         double e = 0.0;
         for (int t = 0; t < nform; t++)
         {
            const PairForm::Term &T = form[t];
            const double ca[3] = { b[T.a].c.x, b[T.a].c.y, b[T.a].c.z };
            double cb[3], ga[3], gb[3];
            if (T.b >= 0)
               cb[0] = b[T.b].c.x, cb[1] = b[T.b].c.y, cb[2] = b[T.b].c.z;
            e -= T.g * PairForm::Q(T, ca, T.b >= 0 ? cb : nullptr, ga, gb);
            if (F != nullptr)
            {
               F[T.a] = F[T.a] + V3(ga[0], ga[1], ga[2]) * T.g;
               if (T.b >= 0)
                  F[T.b] = F[T.b] + V3(gb[0], gb[1], gb[2]) * T.g;
            }
         }
         return e;
      }
      void AddShape(const Shape *const s)
      {
         shape.b.Grow((size_t)ns + 1);
         shape.b[ns++] = s;
      }
      void AddPart(Part *const p)
      {
         part.b.Grow((size_t)np + 1);
         partOff.Grow((size_t)np + 1);
         part.b[np++] = p;
      }

      // After the balls and shapes are in: the contact list (every ball on every shape, every ball pair), with the touches already present
      void Configure()
      {
         {
            size_t n = set.ballPairs ? (size_t)nb * (nb - 1) / 2 : 0;
            for (int s = 0; s < ns; s++)
               n += shape[s]->onlyBall >= 0 ? 1 : nb;
            item.Grow(n);
            law.Grow(n);
            hNow.Grow(n);
            h1c.Grow(n);
            hmc.Grow(n);
            d0s.Grow(n);
            ie.Grow(n);
            farNow.Grow(n);
            prevSlipS.Grow(n);
            prevOnS.Grow(n);
            esI0.Grow(n);
            size_t nU = 6 * (size_t)nb;
            for (int p = 0; p < np; p++)
               nU += part[p]->MaxDof();
            coordOff.Grow(np);
            nCoord = 0;
            for (int p = 0; p < np; p++)
            {
               coordOff[p] = nCoord;
               nCoord += part[p]->NumCoords();
            }
            ci.Grow(nCoord);
            cUnk.Grow(nCoord);
            qEnd.Grow(nCoord);
            qMid.Grow(nCoord);
            qF.Grow(nCoord);
            loadF.Grow(np);
            carried = false;
            for (int s = 0; s < ns; s++)
               carried = carried || shape[s]->part >= 0;
            x.Grow(nU);
            resS.Grow(nU);
            stpS.Grow(nU);
            xPrevS.Grow(nU);
            stpPrevS.Grow(nU);
            trialS.Grow(nU);
            zeroS.Grow(nU);
            esRes.Grow(nU);
            idxS.Grow(nU);
            if ((int)nU > jStride)
            {
               jStride = (int)nU;
               jStore = Buf<double>();
               jStore.Grow(nU * nU);
               jacAnaS = Buf<double>();
               jacFdS = Buf<double>();
            }
            if (set.jacMode != 0)
            {
               const size_t nJ = (size_t)jStride;
               jacAnaS.Grow(nJ * nJ);
               jacFdS.Grow(nJ * nJ);
               xpS.Grow(nJ);
               rpS.Grow(nJ);
               rmS.Grow(nJ);
            }
         }
         set.tolEndRes = 0.0;
         for (int p = 0; p < np; p++)
            set.tolEndRes = std::max(set.tolEndRes, part[p]->EndResidual());
         ni = 0;
         for (int a = 0; a < nb; a++)
            for (int s = 0; s < ns; s++)
            {
               if (shape[s]->onlyBall >= 0 && shape[s]->onlyBall != a)
                  continue;
               Item &it = item[ni++];
               it = Item();
               it.a = a;
               it.s = s;
            }
         if (set.ballPairs)
            for (int a = 0; a < nb; a++)
               for (int b = a + 1; b < nb; b++)
               {
                  Item &it = item[ni++];
                  it = Item();
                  it.a = a;
                  it.b = b;
               }
         for (int k = 0; k < ni; k++)
         {
            Hit hp;
            Probe(item[k], ball, hp);
            if (hp.delta > 0.0)
            {
               TakeLaw(item[k], hp);
               item[k].touching = true;
            }
         }
         books = Books();
         stamp++;
      }

      // VPX_GROUP_BOOK_STEPS=2 (a diagnostic): every stored term of the group, one string: per ball place, kinetic, height; per item its depth, zero, squash, shear, roll; per part its energy
      std::string Terms() const
      {
         std::string o;
         char t[400];
         for (int i = 0; i < nb; i++)
         {
            const Ball &b = ball[i];
            snprintf(t, sizeof(t), "|b%d c=%.9g,%.9g,%.9g ke=%.9g pe=%.9g", i, b.c.x, b.c.y, b.c.z, (0.5 * b.m * b.v.Dot(b.v) + 0.5 * b.I() * b.w.Dot(b.w)) * 1e6, -b.m * set.g.Dot(b.c) * 1e6);
            o += t;
         }
         for (int k = 0; k < ni; k++)
         {
            const Item &it = item[k];
            Hit hp;
            Probe(it, ball, hp);
            const double vn = it.touching ? it.V(hp.delta) : 0.0, vs = it.touching ? 0.5 * it.Kt(hp.delta) * it.xi.Dot(it.xi) : 0.0;
            const double vr = it.touching && set.rollModel == 1 ? 0.5 * it.Kr(hp.delta) * it.ra.th.Dot(it.ra.th) + 0.5 * it.Ktw(hp.delta) * it.ra.tw * it.ra.tw : 0.0;
            snprintf(t, sizeof(t), "|i%d a%d b%d s%d t%d d=%.6g z=%.6g Vn=%.6g Vs=%.6g Vr=%.6g", k, it.a, it.b, it.s, it.touching ? 1 : 0, hp.delta * 1e6, it.dZero * 1e6, vn * 1e6, vs * 1e6, vr * 1e6);
            o += t;
         }
         for (int q = 0; q < np; q++)
         {
            snprintf(t, sizeof(t), "|p%d E=%.9g", q, part[q]->Energy(ball, nb) * 1e6);
            o += t;
         }
         snprintf(t, sizeof(t), "|form=%.9g", FormEnergy(ball) * 1e6);
         return o + t;
      }
      double Energy() const
      {
         double e = 0.0;
         for (int i = 0; i < nb; i++)
         {
            const Ball &b = ball[i];
            e += 0.5 * b.m * b.v.Dot(b.v) + 0.5 * b.I() * b.w.Dot(b.w) - b.m * set.g.Dot(b.c);
            if (nf > 0 && b.fieldScale != 0.0)
               e += FieldEnergy(b);
         }
         for (int k = 0; k < ni; k++)
            if (item[k].touching)
            {
               Hit hp;
               Probe(item[k], ball, hp);
               e += item[k].V(hp.delta) + 0.5 * item[k].Kt(hp.delta) * item[k].xi.Dot(item[k].xi);
               if (set.rollModel == 1)
                  e += 0.5 * item[k].Kr(hp.delta) * item[k].ra.th.Dot(item[k].ra.th) + 0.5 * item[k].Ktw(hp.delta) * item[k].ra.tw * item[k].ra.tw;
            }
         for (int p = 0; p < np; p++)
            e += part[p]->Energy(ball, nb);
         e += FormEnergy(ball);
         return e;
      }
      // What the ball's place stores in the group's fields; its slopes when wanted
      double FieldEnergy(const Ball &b, V3 *const grad = nullptr, M3 *const hess = nullptr) const
      {
         double e = 0.0;
         for (int f = 0; f < nf; f++)
         {
            V3 g;
            M3 hs;
            e += b.fieldScale * field[f]->Eval(b.c, b.R, grad != nullptr ? &g : nullptr, hess != nullptr ? &hs : nullptr);
            if (grad != nullptr)
               *grad = *grad + g * b.fieldScale;
            if (hess != nullptr)
               for (int i = 0; i < 3; i++)
                  for (int j = 0; j < 3; j++)
                     hess->a[i][j] += b.fieldScale * hs.a[i][j];
         }
         return e;
      }
      bool Touching(int k) const { return item[k].touching; }
      void ProbeItem(const int k, Hit &hp) const { Probe(item[k], ball, hp); }
      // Steps until dt is covered. capStep: the caller's own bound on a step; afterStep: after each. False: a step failed, and done says how far it got (time is lost, no energy can appear)
      template <class CapStep, class AfterStep>
      bool StepsOver(const double dt, double &done, const CapStep &capStep, const AfterStep &afterStep)
      {
         done = 0.0;
         int guard = 0;
         while (done < dt - 1e-12)
         {
            const double h = capStep(std::min(ChooseStep(), dt - done));
            const double taken = Step(h);
            if (books.failed || taken <= 0.0 || (++guard > 3000 && GroupFixedStepS() <= 0.0))
               return false;
            done += taken;
            afterStep();
         }
         return true;
      }
      // A carried surface's point speed now, from the hit's columns
      V3 SurfaceSpeed(const Shape &sh, const Hit &hp) const
      {
         V3 vs;
         for (int j = 0; j < hp.nc; j++)
            vs = vs + hp.c[j].V() * part[sh.part]->Coord(hp.c[j].coord).rate;
         return vs;
      }
      // The first guess starts afresh: a step taken again from a snapshot then gives the same result to the last digit
      void ClearGuess()
      {
         for (int i = 0; i < nb; i++)
            aLast[i] = alphaLast[i] = V3();
      }
      // After Configure: a part's own contact goes on as it was kept (touching, stretch, zero); one not touching that is already inside has stored nothing, its depth is its zero
      void SetShapeItem(const int k, const bool wasTouching, const double zero, const V3 &xi, const RubberContact::RollAngle &ra = RubberContact::RollAngle())
      {
         Item &it = item[k];
         Hit hp;
         it.dZero = 0.0;
         if (!wasTouching && zero <= 0.0)
         {
            Probe(it, ball, hp);
            it.dZero = std::max(hp.delta, 0.0);
         }
         else
            it.dZero = zero;
         it.xi = xi;
         it.ra = ra;
         Probe(it, ball, hp);
         it.touching = hp.delta > 0.0;
         if (it.touching)
            TakeLaw(it, hp);
         else
         {
            it.xi = V3();
            it.ra = RubberContact::RollAngle();
         }
      }

      // After Configure: a ball pair found inside with no work done counts from its zero, as a part's own contact does
      void SetPairItem(const int k, const double zero)
      {
         Item &it = item[k];
         Hit hp;
         it.dZero = zero;
         Probe(it, ball, hp);
         it.touching = hp.delta > 0.0;
         if (it.touching)
            TakeLaw(it, hp);
         else
         {
            it.xi = V3();
            it.ra = RubberContact::RollAngle();
         }
      }

      // What set the last ChooseStep (diagnostics for the offline trace; always kept): the rule (0 hMax, 1 normal stiffness, 2 shear stiffness, 3 arrival,
      // 4 turning, 5 floor 0.1 us, 7 a steady contact's normal stiffness at omegaHRest, 8 a part: partRule is its own), the contact (-1 = none) and that
      // contact's values at the choice
      struct StepDiag
      {
         int rule = 0, item = -1, partRule = -1;
         double delta = 0.0, dDot = 0.0, omN = 0.0, omT = 0.0;
      };
      mutable StepDiag diag;

      // The step choice from the state now (accuracy only)
      double ChooseStep() const
      {
         double h = set.hMax;
         StepDiag dg;
         V3 *const fNet = fNetS;   // the net force on each ball now: gravity, each touching contact's spring and shear stretch
         for (int i = 0; i < nb; i++)
         {
            fNet[i] = (set.g * ball[i].m) + ball[i].fExt;
            if (nf > 0 && ball[i].fieldScale != 0.0)
            {
               V3 gf;
               FieldEnergy(ball[i], &gf);
               fNet[i] = fNet[i] - gf;
            }
         }
         FormEnergy(ball, fNet);
         for (int k = 0; k < ni; k++)
         {
            const Item &it = item[k];
            Hit &hp = hNow[k];   // kept: the step's start probes the same positions
            Probe(it, ball, hp);
            if (it.touching && hp.delta > 0.0)
            {
               const V3 f = hp.n * it.F(hp.delta) - it.xi * it.Kt(hp.delta);
               fNet[it.a] = fNet[it.a] + f;
               if (it.b >= 0)
                  fNet[it.b] = fNet[it.b] - f;
            }
         }
         for (int k = 0; k < ni; k++)
         {
            const Item &it = item[k];
            Hit &hp = hNow[k];
            const Ball &A = ball[it.a];
            // a carried surface: its point's speed now, and what the ball meets through it: the part's body at the contact
            const int nq = it.b < 0 ? hp.nc : 0;
            V3 vs;
            double invQ = 0.0;
            for (int j = 0; j < nq; j++)
            {
               const Part &pq = *part[shape[it.s]->part];
               const CoordInfo q = pq.Coord(hp.c[j].coord);
               const double sl = hp.c[j].slope * hp.grad;
               vs = vs + hp.c[j].V() * q.rate;
               invQ += sl * sl * (pq.DofOf(hp.c[j].coord) >= 0 ? q.invM : q.blockedSide == 0 || q.blockedSide * sl > 0.0 ? q.invMFree : 0.0);
            }
            const V3 vr = it.b >= 0 ? A.v - ball[it.b].v : A.v - vs;
            const double dDot = -vr.Dot(hp.n);   // closing speed
            const double vc = std::max(dDot, set.vFloor);
            if (hp.delta <= 0.0 && -hp.delta >= vc * set.hMax + 1e-9)
               continue;   // cannot be reached within the longest step
            Item L = it;
            if (!it.touching)
               TakeLaw(L, hp);
            if (invQ > 0.0)
               L.mRed = 1.0 / (1.0 / A.m + invQ);
            // the depth this approach reaches (all the closing energy into the squash): the stiffness while it lasts
            const double eC = 0.5 * L.mRed * vc * vc;
            double lo = std::pow(2.5 * eC / L.kH, 0.4);   // the Hertz spring alone; a layer only makes it shallower
            if (L.cL > 0.0)
            {
               double hi = std::min(lo, L.hL);
               lo = 0.0;
               for (int r = 0; r < 50 && hi - lo > 1e-4 * hi; r++)
               {
                  const double mid = 0.5 * (lo + hi);
                  (L.V(mid) < eC ? lo : hi) = mid;
               }
            }
            if (L.kb > 0.0)
            {
               double uLo = 0.0, uHi = lo;
               for (int r = 0; r < 30 && uHi - uLo > 1e-3 * uHi; r++)
               {
                  const double u = 0.5 * (uLo + uHi), f = L.kH * u * std::sqrt(u);
                  (0.4 * f * u + 0.5 * f * f / L.kb < eC ? uLo : uHi) = u;
               }
               const double fE = L.kH * uHi * std::sqrt(uHi);
               lo = uHi + fE / L.kb;
            }
            const double dE = std::max(hp.delta, lo);
            const double ktE = L.Kt(dE);
            const double omN = std::sqrt(L.K(dE) / L.mRed), omT = std::sqrt(ktE / L.mt);
            // steady (Settings): by the motion, at the contact's own time scale. At the friction limit (Coulomb): a spring within mu N holds at most (mu N)^2 / (2 kt) of energy, so a slip with more than that
            // in 1/2 mt u^2 cannot be stopped by it: the contact slides, its force is mu N whatever the spring's stiffness, and the spring's time scale is not in the motion.
            bool steady = false;
            if (set.vSteady > 0.0)
            {
               if (hp.delta > 0.0)
               {
                  V3 u = A.v + A.w.Cross(hp.n * (-(A.R - 0.5 * hp.delta)));
                  if (it.b >= 0)
                     u = u - (ball[it.b].v + ball[it.b].w.Cross(hp.n * (ball[it.b].R - 0.5 * hp.delta)));
                  else if (nq > 0)
                     u = u - vs;
                  const V3 ut = u - hp.n * u.Dot(hp.n);
                  const double slipNow = ut.Length();
                  V3 aRel = fNet[it.a] * (1.0 / A.m);
                  if (it.b >= 0)
                     aRel = aRel - fNet[it.b] * (1.0 / ball[it.b].m);
                  const double omNow = std::sqrt(L.K(hp.delta) / L.mRed);
                  steady = std::fabs(dDot) <= set.vSteady && std::fabs(aRel.Dot(hp.n)) <= set.vSteady * omNow && (slipNow <= set.vSteady || (set.slideFree && slipNow * std::sqrt(ktE * L.mt) > L.mu * L.F(dE)));
               }
               else
                  steady = dDot <= set.vSteady;
            }
            const double om = steady ? omN : (set.shearInStep ? std::max(omN, omT) : omN);
            const double hs = om > 0.0 ? (steady ? set.omegaHRest : set.omegaH) / om : set.hMax;
            const auto lim = [&](const double c, const int rule) {
               if (c < h)
               {
                  h = c;
                  dg.rule = rule;
                  dg.item = k;
                  dg.delta = hp.delta;
                  dg.dDot = dDot;
                  dg.omN = omN;
                  dg.omT = omT;
               }
            };
            if (hp.delta < 0.0 && dDot > 0.0)   // reached in about two steps, not overshot
            {
               const double ha = 0.5 * -hp.delta / dDot;
               lim(std::max(ha, hs), ha > hs ? 3 : (set.shearInStep && omT > omN ? 2 : 1));
            }
            else
               lim(hs, steady ? 7 : (set.shearInStep && omT > omN ? 2 : 1));
            if (hp.curvDist > 0.0 && hp.delta > -1e-4)   // rolling round a post or a ball
            {
               const V3 vt = vr - hp.n * vr.Dot(hp.n);
               const double sp = vt.Length();
               if (sp > 1e-9)
                  lim(set.turnMaxRad * hp.curvDist / sp, 4);
            }
         }
         for (int p = 0; p < np; p++)
         {
            StepFactList facts(set);
            part[p]->StepFacts(ball, nb, set, facts);
            const double hp = std::max(facts.H(), 1e-7);
            if (hp < h)
            {
               h = hp;
               dg = StepDiag();
               dg.rule = 8;
               dg.partRule = facts.Rule();
            }
         }
         if (h < 1e-7)
            dg.rule = 5;
         diag = dg;
         nowStamp = stamp;
         return std::max(h, 1e-7);
      }

      static constexpr double kOnStop = 1e-10;   // (rad) a part this close to its hard stop is on it
      // One step of at most hTry (halved while the solve fails; shortened to end on a part's hard stop). Returns the time taken, 0 = failed.
      // The events of a step, in order, one loop for every part: a body with no unknown freed by the step's forces, or one already on the stop it would pass, changes the unknowns, and the step is solved again;
      // otherwise a stop passed shortens the step to end on it.
      double Step(double hTry)
      {
         double h = hTry;
         nHalveDiag = 0;
         for (;;)
         {
            if (Solve(h))
            {
               {
                  bool freed = false;
                  for (int p = 0; p < np; p++)
                     freed = part[p]->FreedByStep() || freed;
                  if (freed)
                     continue;
               }
               if (StopExcess() <= 0.0)
               {
                  Accept(h, false);
                  return h;
               }
               // the step that ends on the stop (Illinois false position on the step length), then solved again there (its end state is the one taken)
               const double *const zero = zeroS;
               {
                  bool pinnedOne = false;
                  for (int p = 0; p < np; p++)
                  {
                     const bool onIt = part[p]->StopExcess(zero, partOff[p]) > -kOnStop;
                     if (part[p]->StopExcess(x, partOff[p]) > 0.0 && onIt)   // (the solved move tested last: ArrivingAtStop and PinOnStop may use what it found)
                        pinnedOne = part[p]->PinOnStop(x, partOff[p]) || pinnedOne;
                  }
                  if (pinnedOne)
                     continue;
               }
               double lo = 0.0, flo = std::min(StopExcess(zero), -1e-30), hi = h, fhi = StopExcess(), hBelow = 0.0;
               int side = 0;
               bool solvedBelow = false;   // x holds the solve at hBelow (the search ended on it): no need to solve it again
               WarmFromSolve(h);
               for (int it = 0; it < 80; it++)
               {
                  double hm = (lo * fhi - hi * flo) / (fhi - flo);
                  if (!(hm > lo && hm < hi))
                     hm = 0.5 * (lo + hi);
                  if (!Solve(hm))
                  {
                     warmOn = false;
                     return EnergyStep(hTry);
                  }
                  WarmFromSolve(hm);
                  const double f = StopExcess();
                  if (f <= 0.0)
                  {
                     hBelow = hm;
                     bool arriving = true;
                     for (int p = 0; p < np; p++)
                        arriving = arriving && part[p]->ArrivingAtStop();
                     if (f > -kOnStop && arriving)
                     {
                        solvedBelow = true;
                        break;
                     }
                     lo = hm;
                     flo = f;
                     if (side == -1)
                        fhi *= 0.5;
                     side = -1;
                  }
                  else
                  {
                     hi = hm;
                     fhi = f;
                     if (side == 1)
                        flo *= 0.5;
                     side = 1;
                  }
               }
               const bool endOk = hBelow > 0.0 && (solvedBelow || Solve(hBelow));
               warmOn = false;
               if (!endOk)
                  return EnergyStep(hTry);
               Accept(hBelow, true);
               books.stopSteps++;
               return hBelow;
            }
            if (nHalveDiag < kHalveDiag)   // why this attempt failed (diagnostics)
            {
               HalveDiag &hd = halveDiag[nHalveDiag++];
               hd.why = failWhy;
               hd.h = h;
               for (int p = 0; p < np; p++)
                  part[p]->KeepFailNote(nHalveDiag - 1, resS, partOff[p]);
               hd.rounds = lastRounds;
               const RoundDiag &rl = roundDiag[std::max(0, std::min(lastRounds, kRounds) - 1)];
               hd.moveC = rl.moveC;
               hd.moveW = rl.moveW;
               for (int r = 0; r < std::min(lastRounds, kRounds); r++)   // over the failed attempt's rounds
               {
                  hd.slipFlips += roundDiag[r].nSlipFlip;
                  hd.touchFlips += roundDiag[r].nOnFlip;
                  hd.maxSliding = std::max(hd.maxSliding, roundDiag[r].nSlip);
               }
               for (int k = 0; k < 6; k++)   // the last six rounds' rim-speed changes
               {
                  const int r = std::min(lastRounds, kRounds) - 6 + k;
                  hd.lastW[k] = r >= 0 ? roundDiag[r].moveW : 0.0;
               }
            }
            h *= 0.5;
            books.halvings++;
            if (h < set.hMin)
               return EnergyStep(hTry);
         }
      }

      // diagnostics of the last solve: its rounds, and per round (the first kRounds) the largest Newton change of a ball's speed (moveC, m/s) and rim speed
      // (moveW = R |dw|, m/s), the share of the Newton move taken (scale), the contacts in play, those sliding, and how many changed between sticking and sliding
      // or between touching and not since the round before
      int lastRounds = 0;
      // why the last Step's attempts failed and were halved: P a part could not begin the step, B a part's forces could not be formed (the band's path would not
      // settle), G a first guess the layer check refused 60 times, N a residual not finite, L the linear solve
      // singular, R a Newton move the layer check refused 60 times, C not solved within maxRounds, E the end evaluation not finite
      struct HalveDiag { char why = '?'; double h = 0.0; int rounds = 0; double moveC = 0.0, moveW = 0.0; int slipFlips = 0, touchFlips = 0, maxSliding = 0; double lastW[6] = {}; };
      static constexpr int kHalveDiag = 8;
      HalveDiag halveDiag[kHalveDiag];
      int nHalveDiag = 0;
      char failWhy = '?';
      // back: a line-search round (the last move retaken shorter); grip: the parts' grip change after the move (m/s)
      struct RoundDiag { double moveC = 0.0, moveW = 0.0, scale = 1.0, grip = 0.0; int nActive = 0, nSlip = 0, nSlipFlip = 0, nOnFlip = 0, back = 0; };
      static constexpr int kRounds = 64;
      RoundDiag roundDiag[kRounds];
      // Settings::jacMode: the solve's own Jacobian and the central-difference one at the last round of the last solve (rows: per ball the momentum residual 3,
      // then the spin residual 3; columns: per ball the move 3, then the spin change 3)
      Jac JacAna() { return Jac(jacAnaS, jStride); }
      Jac JacFd() { return Jac(jacFdS, jStride); }
      int jacRound = 0;

   private:
      // scratch of the step being solved (no allocation)
      struct ItemEnd
      {
         bool active = false, slipping = false;
         double d0 = 0.0, d1 = 0.0;
         V3 xi1;
         RubberContact::RollAngle ra1;   // the rolling and twist springs' end angles (roll_model 1)
         double normalLoss = 0.0, slip = 0.0, shearDamp = 0.0, twist = 0.0, roll = 0.0;
      };
      Buf<Item> law;             // the law each contact uses in this step (held, or taken at the step's start)
      // Probes kept instead of repeated (a tube probe is most of a step's cost): the positions now (ChooseStep; valid while nowStamp == stamp), and the end and
      // mid-step positions of the move the range check accepted (the next Eval looks at the same positions)
      mutable Buf<Hit> hNow;
      mutable long nowStamp = -1;
      long stamp = 0;            // counts changes of the positions (Configure, Accept)
      Buf<Hit> h1c, hmc;
      bool trialCached = false;
      Buf<double> d0s;
      Buf<Ball> b1r, bmr;        // InRange's moved balls
      Buf<ItemEnd> ie;
      Buf<double> x;             // the unknowns: per ball the move (3) and the spin change (3), then the parts' own
      Buf<V3> aLast, alphaLast;  // each ball's acceleration and spin's over the last step (the next step's first guess)
      // what each function of a step works in
      mutable Buf<V3> fNetS;
      Buf<V3> evF, evM, evDc, evDw, evWm, esA0, esAl0;
      // the stop search's first guess: the accelerations of this step's last solve (same start state), so a re-solve at a nearby length starts near its answer
      Buf<V3> warmA, warmAl;
      bool warmOn = false;
      void WarmFromSolve(const double h)
      {
         for (int i = 0; i < nb; i++)
         {
            const V3 dc(x[6 * i], x[6 * i + 1], x[6 * i + 2]), dw(x[6 * i + 3], x[6 * i + 4], x[6 * i + 5]);
            warmA[i] = (dc - ball[i].v * h) * (2.0 / (h * h));   // dc = v h + a h^2 / 2
            warmAl[i] = dw * (1.0 / h);
         }
         warmOn = true;
      }
      Buf<Ball> evB1, evBm, farB1, farBm, esB0;
      Buf<Item> esI0;
      Buf<double> resS, stpS, xPrevS, stpPrevS, trialS, zeroS, esRes, xpS, rpS, rmS, jStore, jacAnaS, jacFdS;
      Buf<char> prevSlipS, prevOnS;
      Buf<int> idxS;
      int jStride = 0;
      double driveW = 0.0;
      Buf<int> partOff;
      // the parts' coordinates, one list (coordOff: each part's first): as they are this attempt, the unknown of each (-1 = none), the moves a probe is at, the balls' reaction on each and on each part in this evaluation
      Buf<int> coordOff, cUnk;
      int nCoord = 0;
      Buf<CoordInfo> ci;
      Buf<double> qEnd, qMid, qF;
      Buf<V3> loadF;
      bool carried = false;      // a shape of this group is carried by a part
      static constexpr int kGradMax = 3 + 7;   // a discrete gradient's entries: a ball's centre and a hit's columns (two balls: 6)
      int nUnk = 0;              // unknowns of the step being solved
      bool partEvalFailed = false;
      Buf<char> farNow;          // not touching and out of reach in this attempt
      Buf<double> fieldE0;       // what each ball's place stores at the start of the attempt
      // out of reach: the gap is more than both bodies' speeds can close in twice the step, with a margin. A body that picks up speed inside the step is what the check after the solve is for
      bool OutOfReach(const Item &it, const double d0, const double h, const Hit &h0) const
      {
         if (it.touching || d0 > 0.0)
            return false;
         double sp = ball[it.a].v.Length();
         if (it.b >= 0)
            sp += ball[it.b].v.Length();
         else
            for (int j = 0; j < h0.nc; j++)
            {
               const CoordInfo q = part[shape[it.s]->part]->Coord(h0.c[j].coord);
               sp += std::fabs(q.rate) * q.reach;
            }
         return -d0 > 2.0 * sp * h + 1e-4;
      }
      // where the coordinates are at the end and the middle of the move xv, for the probes of a carried surface
      void SetTrials(const double *xv)
      {
         for (int c = 0; c < nCoord; c++)
         {
            const double q = cUnk[c] >= 0 ? xv[cUnk[c]] : ci[c].given;
            qEnd[c] = q;
            qMid[c] = 0.5 * q;
         }
      }
      // the probes of a move: at its end and its middle (a carried surface moves with its part's coordinates)
      void ProbeMove(const Item &it, const Ball *b1, const Ball *bm, Hit &h1, Hit &hm) const
      {
         const Shape *const sh = it.b < 0 ? (const Shape *)shape[it.s] : nullptr;
         if (sh != nullptr && sh->part >= 0)
         {
            sh->trial = (const double *)qEnd + coordOff[sh->part];
            Probe(it, b1, h1);
            sh->trial = (const double *)qMid + coordOff[sh->part];
            Probe(it, bm, hm);
            sh->trial = nullptr;
            return;
         }
         Probe(it, b1, h1);
         Probe(it, bm, hm);
      }
      // the coordinates as this attempt has them (after the parts' BeginStep)
      void CoordsNow()
      {
         for (int p = 0; p < np; p++)
         {
            const int n = part[p]->NumCoords();
            for (int c = 0; c < n; c++)
            {
               const int i = coordOff[p] + c, d = part[p]->DofOf(c);
               ci[i] = part[p]->Coord(c);
               cUnk[i] = d >= 0 ? partOff[p] + d : -1;
            }
         }
      }
      // the parts' rows in the line search's speed units: over the mass, times the reach
      double CoordRes2(const double *Res) const
      {
         double r = 0.0;
         for (int c = 0; c < nCoord; c++)
            if (cUnk[c] >= 0)
            {
               const double s = Res[cUnk[c]] * ci[c].invM * ci[c].reach;
               r += s * s;
            }
         return r;
      }
      // After the solve: did a skipped contact get reached after all? Then the attempt runs again with none skipped.
      bool FarStillFar(const double *xv)
      {
         if (carried)
            SetTrials(xv);
         Ball *const b1 = farB1, *const bm = farBm;
         for (int i = 0; i < nb; i++)
         {
            const V3 dc(xv[6 * i], xv[6 * i + 1], xv[6 * i + 2]);
            b1[i] = bm[i] = ball[i];
            b1[i].c = ball[i].c + dc;
            bm[i].c = ball[i].c + dc * 0.5;
         }
         for (int k = 0; k < ni; k++)
            if (farNow[k])
            {
               Hit h1, hm;
               ProbeMove(item[k], b1, bm, h1, hm);
               if (h1.delta > 0.0 || hm.delta > 0.0)
                  return false;
            }
         return true;
      }

      // the number of unknowns (the parts' own placed after the balls')
      int NumUnk()
      {
         int n = 6 * nb;
         for (int p = 0; p < np; p++)
         {
            partOff[p] = n;
            n += part[p]->NumDof();
         }
         return n;
      }
      double StopExcess(const double *xv = nullptr) const
      {
         double e = -INFINITY;
         for (int p = 0; p < np; p++)
            e = std::max(e, part[p]->StopExcess(xv ? xv : x, partOff[p]));
         return e;
      }

      void Probe(const Item &it, const Ball *bs, Hit &hp) const
      {
         const Ball &A = bs[it.a];
         hp.kb = 0.0;
         hp.nc = 0;   // (hits are kept and reused: a fixed surface writes no columns, and must not inherit a carried one's)
         if (it.b < 0)
         {
            shape[it.s]->Probe(A.c, A.R, hp);
            hp.delta -= it.dZero;
            return;
         }
         const Ball &B = bs[it.b];
         const V3 d = A.c - B.c;
         const double D = d.Length();
         hp.delta = A.R + B.R - D - it.dZero;
         hp.grad = 1.0;
         hp.n = D > 0.0 ? d * (1.0 / D) : V3(0.0, 0.0, 1.0);
         hp.reff = A.R * B.R / (A.R + B.R);
         hp.curvDist = D;
      }
      // the law of a touch that starts: material, effective radius and layer from the surface (ball pairs: the pair's material, two equal bodies)
      void TakeLaw(Item &it, const Hit &hp) const
      {
         const Ball &A = ball[it.a];
         const RubberParams &rp = it.b >= 0 ? set.pair : shape[it.s]->rp;
         double E = rp.eMpa * 1e6;
         double Estar = E / (1.0 - rp.nu * rp.nu);
         const double G = E / (2.0 * (1.0 + rp.nu)), Gstar = G / (2.0 - rp.nu);
         double Gs = Gstar;
         if (it.b >= 0)   // two bodies of the same material: the compliances add
         {
            Estar *= 0.5;
            Gs *= 0.5;
         }
         it.reff = hp.reff;
         it.kb = it.b >= 0 ? 0.0 : hp.kb;
         it.kH = 4.0 / 3.0 * Estar * std::sqrt(it.reff);
         it.hL = it.b >= 0 ? 0.0 : shape[it.s]->layerM;
         it.cL = it.hL > 0.0 ? it.kH * rp.betaLayer * it.reff * std::sqrt(it.reff) : 0.0;
         it.ktPerA = 8.0 * Gs * rp.shearScale;
         it.E = E;
         double invMt = 1.0 / A.m + A.R * A.R / A.I();
         it.mRed = A.m;
         if (it.b >= 0)
         {
            const Ball &B = ball[it.b];
            it.mRed = A.m * B.m / (A.m + B.m);
            invMt += 1.0 / B.m + B.R * B.R / B.I();
         }
         it.mt = 1.0 / invMt;
         it.alpha = rp.alphaSPerM;
         it.gammaKK = rp.gammaKK;
         it.tanDelta = rp.tanDelta;
         it.normalLoss = rp.normalLoss;
         it.shearTan = rp.shearTan;
         it.mu = rp.mu;
         it.rollingMu = rp.rollingMu;
         it.twist = rp.twist;
      }
      // the depth a move may not reach: a layer squashed to nothing
      bool InRange(const double *xv)
      {
         if (carried)
            SetTrials(xv);
         trialCached = false;
         bool moved = false;
         for (int k = 0; k < ni; k++)
            if (!farNow[k] && law[k].cL > 0.0)
            {
               if (!moved)
               {
                  moved = true;
                  for (int i = 0; i < nb; i++)
                  {
                     b1r[i] = bmr[i] = ball[i];
                     const V3 dc(xv[6 * i], xv[6 * i + 1], xv[6 * i + 2]);
                     b1r[i].c = ball[i].c + dc;
                     bmr[i].c = ball[i].c + dc * 0.5;
                  }
               }
               const Ball *const b1 = b1r, *const bm = bmr;
               Hit &h1 = h1c[k], &hm = hmc[k];
               ProbeMove(item[k], b1, bm, h1, hm);
               const double lim = law[k].hL * (1.0 - 1e-9);
               if (h1.delta >= lim || hm.delta >= lim)
                  return false;
            }
         for (int p = 0; p < np; p++)
            if (!part[p]->InRange(xv, partOff[p]))
               return false;
         trialCached = true;   // the layered contacts' probes at this move are in h1c / hmc
         return true;
      }

      // Gonzalez: gm + (df - gm . dq) / (dq^T W dq) W dq, so that the result . dq = df exactly; a mismatch within the rounding of f is left out
      static void Gonzalez(const double *gm, double df, double dfNoise, const double *dq, const double *W, int n, double *out)
      {
         double gdq = 0.0, qq = 0.0;
         for (int a = 0; a < n; a++)
         {
            gdq += gm[a] * dq[a];
            qq += dq[a] * W[a] * dq[a];
         }
         const double cf = qq > 1e-300 && std::fabs(df - gdq) > dfNoise ? (df - gdq) / qq : 0.0;
         for (int a = 0; a < n; a++)
            out[a] = gm[a] + cf * W[a] * dq[a];
      }
      // The same elimination over the nonzero entries of each pivot row: a group's bodies touch few others
      static bool SolveLinear(const Jac A, double *b, int n, int *const idx)
      {
         for (int col = 0; col < n; col++)
         {
            int piv = col;
            double big = std::fabs(A[col][col]);
            for (int rw = col + 1; rw < n; rw++)
            {
               const double v = std::fabs(A[rw][col]);
               if (v > big)
               {
                  big = v;
                  piv = rw;
               }
            }
            if (big < 1e-300)
               return false;
            if (piv != col)
            {
               for (int k = col; k < n; k++)
                  std::swap(A[col][k], A[piv][k]);
               std::swap(b[col], b[piv]);
            }
            int nz = 0;
            const double *const ap = A[col];
            for (int k = col + 1; k < n; k++)
               if (ap[k] != 0.0)
                  idx[nz++] = k;
            const double ip = 1.0 / ap[col], bc = b[col];
            for (int rw = col + 1; rw < n; rw++)
            {
               double *const ar = A[rw];
               if (ar[col] == 0.0)
                  continue;
               const double f = ar[col] * ip;
               for (int q = 0; q < nz; q++)
                  ar[idx[q]] -= f * ap[idx[q]];
               ar[col] = 0.0;
               b[rw] -= f * bc;
            }
         }
         for (int rw = n - 1; rw >= 0; rw--)
         {
            double v = b[rw];
            const double *const ar = A[rw];
            for (int k = rw + 1; k < n; k++)
               if (ar[k] != 0.0)
                  v -= ar[k] * b[k];
            b[rw] = v / ar[rw];
         }
         return true;
      }
      static V3 MulV(const M3 &m, const V3 &v)
      {
         return V3(m.a[0][0] * v.x + m.a[0][1] * v.y + m.a[0][2] * v.z, m.a[1][0] * v.x + m.a[1][1] * v.y + m.a[1][2] * v.z, m.a[2][0] * v.x + m.a[2][1] * v.y + m.a[2][2] * v.z);
      }
      static void AddBlock(const Jac J, int r0, int c0, const M3 &m, double s)
      {
         for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
               J[r0 + i][c0 + j] += s * m.a[i][j];
      }

      // Forces of the step at the unknowns xv: the residual R (and the approximate Jacobian J when given), the contacts' end values into ie[]
      bool Eval(double h, const double *xv, double *Res, const Jac J, bool useCache)
      {
         const int nU = nUnk;
         V3 *const F = evF, *const M = evM, *const dc = evDc, *const dw = evDw, *const wm = evWm;
         Ball *const b1 = evB1, *const bm = evBm;
         double posScale = 0.1;
         for (int i = 0; i < nb; i++)
         {
            dc[i] = V3(xv[6 * i], xv[6 * i + 1], xv[6 * i + 2]);
            dw[i] = V3(xv[6 * i + 3], xv[6 * i + 4], xv[6 * i + 5]);
            wm[i] = ball[i].w + dw[i] * 0.5;
            b1[i] = bm[i] = ball[i];
            b1[i].c = ball[i].c + dc[i];
            bm[i].c = ball[i].c + dc[i] * 0.5;
            F[i] = set.g * ball[i].m + ball[i].fExt;
            M[i] = V3();
            posScale = std::max(posScale, std::fabs(ball[i].c.x) + std::fabs(ball[i].c.y) + std::fabs(ball[i].c.z));
         }
         if (J)
            for (int r = 0; r < nU; r++)
               for (int c = 0; c < nU; c++)
                  J[r][c] = 0.0;
         const double noiseD = 16.0 * 2.220446049250313e-16 * posScale;
         if (carried)
         {
            SetTrials(xv);
            for (int c = 0; c < nCoord; c++)
               qF[c] = 0.0;
            for (int p = 0; p < np; p++)
               loadF[p] = V3();
         }
         driveW = 0.0;
         for (int i = 0; i < nb; i++)
            driveW += ball[i].fExt.Dot(dc[i]);
         // the spin brake: on the mid-step spin across its axis, so it can only take energy, and that energy is what Accept books
         for (int i = 0; i < nb; i++)
            if (ball[i].brakeRate > 0.0)
            {
               const double c = ball[i].brakeRate * ball[i].I();
               M[i] = M[i] - (wm[i] - ball[i].brakeAxis * wm[i].Dot(ball[i].brakeAxis)) * c;
               if (J)
                  AddBlock(J, 6 * i + 3, 6 * i + 3, M3::Tangent(ball[i].brakeAxis), 0.5 * h * c);   // (the mid-step spin moves half as far as the end)
            }
         // the fields' stored energy: its discrete gradient over the ball's move (Gonzalez), so the work it does is the exact difference at any step length
         for (int i = 0; nf > 0 && i < nb; i++)
            if (ball[i].fieldScale != 0.0)
            {
               V3 gf;
               M3 hf;
               const double e1 = FieldEnergy(b1[i]);
               FieldEnergy(bm[i], &gf, J ? &hf : nullptr);
               const double gm[3] = { gf.x, gf.y, gf.z }, dq[3] = { dc[i].x, dc[i].y, dc[i].z }, W[3] = { 1.0, 1.0, 1.0 };
               double gd[3];
               Gonzalez(gm, e1 - fieldE0[i], 16.0 * 2.220446049250313e-16 * std::fabs(e1) + gf.Length() * noiseD, dq, W, 3, gd);
               F[i] = F[i] - V3(gd[0], gd[1], gd[2]);
               if (J)
                  AddBlock(J, 6 * i, 6 * i, hf, 0.5 * h);   // (the midpoint moves half as far as the end)
            }
         // the pair form: each term's discrete gradient over the move (Gonzalez), so its work is the exact change of -g q. Its curvature (about g / l, 0.4 N/mm against
         // the contacts' 3.5e3 N/mm) is left out of the Jacobian: it changes the rounds, not the answer
         for (int t = 0; t < nform; t++)
         {
            const PairForm::Term &T = form[t];
            const bool two = T.b >= 0;
            const double c0a[3] = { ball[T.a].c.x, ball[T.a].c.y, ball[T.a].c.z }, c1a[3] = { b1[T.a].c.x, b1[T.a].c.y, b1[T.a].c.z }, cma[3] = { bm[T.a].c.x, bm[T.a].c.y, bm[T.a].c.z };
            double c0b[3] = { 0, 0, 0 }, c1b[3] = { 0, 0, 0 }, cmb[3] = { 0, 0, 0 };
            if (two)
            {
               c0b[0] = ball[T.b].c.x, c0b[1] = ball[T.b].c.y, c0b[2] = ball[T.b].c.z;
               c1b[0] = b1[T.b].c.x, c1b[1] = b1[T.b].c.y, c1b[2] = b1[T.b].c.z;
               cmb[0] = bm[T.b].c.x, cmb[1] = bm[T.b].c.y, cmb[2] = bm[T.b].c.z;
            }
            double gm[6] = { 0, 0, 0, 0, 0, 0 };
            const double q0 = PairForm::Q(T, c0a, two ? c0b : nullptr), q1 = PairForm::Q(T, c1a, two ? c1b : nullptr);
            PairForm::Q(T, cma, two ? cmb : nullptr, gm, gm + 3);
            const double dq[6] = { dc[T.a].x, dc[T.a].y, dc[T.a].z, two ? dc[T.b].x : 0.0, two ? dc[T.b].y : 0.0, two ? dc[T.b].z : 0.0 };
            const double W[6] = { 1.0 / ball[T.a].m, 1.0 / ball[T.a].m, 1.0 / ball[T.a].m, two ? 1.0 / ball[T.b].m : 0.0, two ? 1.0 / ball[T.b].m : 0.0, two ? 1.0 / ball[T.b].m : 0.0 };
            double gd[6];
            Gonzalez(gm, q1 - q0, 16.0 * 2.220446049250313e-16 * std::fabs(q1) + noiseD, dq, W, two ? 6 : 3, gd);
            F[T.a] = F[T.a] + V3(gd[0], gd[1], gd[2]) * T.g;
            if (two)
               F[T.b] = F[T.b] + V3(gd[3], gd[4], gd[5]) * T.g;
         }

         for (int k = 0; k < ni; k++)
         {
            const Item &it = item[k];
            const Item &L = law[k];
            ItemEnd &e = ie[k];
            if (farNow[k])
            {
               e = ItemEnd();
               e.active = false;
               e.d0 = d0s[k];
               e.d1 = -1.0;
               continue;
            }
            e = ItemEnd();
            Hit h1, hm;
            if (useCache && L.cL > 0.0)
            {
               h1 = h1c[k];
               hm = hmc[k];
            }
            else
               ProbeMove(it, b1, bm, h1, hm);
            const double d0 = d0s[k], d1 = h1.delta, dm = hm.delta;
            e.d0 = d0;
            e.d1 = d1;
            if (d0 <= 0.0 && d1 <= 0.0 && dm <= 0.0)
               continue;
            e.active = true;
            const bool pr = it.b >= 0;
            const int na = pr ? 6 : 3;
            const int nq = pr ? 0 : hm.nc;   // a carried surface's columns
            const int qc0 = nq > 0 ? coordOff[shape[it.s]->part] : 0;
            const V3 &nm = hm.n;
            // the compression's discrete gradient over the move (ball a's centre, then ball b's, or the columns' coordinates)
            double gm[kGradMax] = { -nm.x * hm.grad, -nm.y * hm.grad, -nm.z * hm.grad, nm.x, nm.y, nm.z }, dq[kGradMax], W[kGradMax], gd[kGradMax];
            const Ball &A = ball[it.a];
            dq[0] = dc[it.a].x; dq[1] = dc[it.a].y; dq[2] = dc[it.a].z;
            W[0] = W[1] = W[2] = 1.0 / A.m;
            if (pr)
            {
               dq[3] = dc[it.b].x; dq[4] = dc[it.b].y; dq[5] = dc[it.b].z;
               W[3] = W[4] = W[5] = 1.0 / ball[it.b].m;
            }
            for (int j = 0; j < nq; j++)
            {
               const int cj = qc0 + hm.c[j].coord;
               gm[3 + j] = hm.c[j].slope * hm.grad;
               dq[3 + j] = qEnd[cj];
               W[3 + j] = cUnk[cj] >= 0 ? ci[cj].invM : 0.0;   // (a move given to it takes none of the correction: only the unknowns do)
            }
            Gonzalez(gm, d1 - d0, noiseD, dq, W, na + nq, gd);
            const double fe = L.Chord(d0, d1);
            double cD, kkJac = 0.0;
            if (L.gammaKK > 0.0)
            {
               cD = L.KKStep(d0, d1, dm, kkJac);
               if (dm > 0.0 && L.normalLoss && L.tanDelta > 0.0)
                  cD += L.tanDelta * std::sqrt(1.5 * L.F(dm) * L.mRed / dm);
            }
            else
               cD = L.DampC(dm, fe);
            const double fd = cD * (d1 - d0) / h;
            const double N = std::max(0.0, fe + fd);   // a contact pushes, it cannot pull (Hertz, no adhesion)
            e.normalLoss = (N - fe) * (d1 - d0);       // (cD (d1 - d0)^2 / h while it pushes; all the spring's release once damping would make it pull)

            // the contact point and the mid-step slip rate
            const double dmp = std::max(dm, 0.0);
            const V3 rA = nm * (-(A.R - 0.5 * dmp));
            const V3 rB = pr ? nm * (ball[it.b].R - 0.5 * dmp) : V3();
            V3 u = dc[it.a] * (1.0 / h) + wm[it.a].Cross(rA);
            if (pr)
               u = u - (dc[it.b] * (1.0 / h) + wm[it.b].Cross(rB));
            for (int j = 0; j < nq; j++)
               u = u - hm.c[j].V() * (dq[3 + j] / h);
            const V3 ut = u - nm * u.Dot(nm);
            // the shear spring: its stretch turned into the contact plane (length kept), the midpoint force, the Coulomb limit
            V3 xi0 = it.touching ? it.xi : V3();
            {
               const double l0 = xi0.Length();
               V3 xp = xi0 - nm * xi0.Dot(nm);
               const double lp = xp.Length();
               xi0 = lp > 1e-300 ? xp * (l0 / lp) : V3();
            }
            // the patch rule (rubbercontact.h ContactCore): a patch that grows holds the force, one that shrinks holds the stretch; the spring's fall is a loss
            const double k0 = L.Kt(d0), k1 = L.Kt(d1);
            const double held = 0.5 * k0 * xi0.Dot(xi0);
            if (k1 > k0)
               xi0 = xi0 * (k0 / k1);
            e.slip = held - 0.5 * k1 * xi0.Dot(xi0);
            // the patch: the spring k1 and the damper ct side by side, in series with the Coulomb slip. The trial force (no slip) is the spring's midpoint force and the
            // damper's; past mu N the patch slips along the trial and the spring's end stretch is the one that carries the limit force with the damper. The force is
            // continuous across stick and slip.
            V3 Ft, xi1, FsTr;
            double ct = 0.0, limS = 0.0, flTr = 0.0;
            bool slipping = false;
            if (k1 > 0.0)
            {
               const V3 xtr = xi0 + ut * h;
               ct = L.shearTan * std::sqrt(k1 * L.mt);
               Ft = (xi0 + xtr) * (0.5 * k1) + ut * ct;
               FsTr = Ft;
               const double lim = L.mu * N, fl = Ft.Length();
               limS = lim;
               flTr = fl;
               xi1 = xtr;
               if (fl > lim)
               {
                  Ft = fl > 0.0 ? Ft * (lim / fl) : V3();
                  // Ft = k1 (xi0 + xi1) / 2 + ct (xi1 - xi0) / h, so the slip xtr - xi1 = (trial - Ft) / (k1 / 2 + ct / h) lies along Ft: its loss Ft . slip >= 0
                  const double cs = 0.5 * k1 + ct / h;
                  xi1 = (Ft - xi0 * (0.5 * k1 - ct / h)) * (1.0 / cs);
                  slipping = true;
                  e.slipping = true;
                  e.slip += Ft.Dot(xtr - xi1);
               }
            }
            e.xi1 = xi1;
            {
               const V3 dxi = xi1 - xi0;   // the damper's loss over its own stretch (ct ut^2 h while it sticks)
               e.shearDamp = slipping ? ct * dxi.Dot(dxi) / h : ct * ut.Dot(ut) * h;
            }
            F[it.a] = F[it.a] - V3(gd[0], gd[1], gd[2]) * N - Ft;
            M[it.a] = M[it.a] + rA.Cross(Ft * -1.0);
            if (pr)
            {
               F[it.b] = F[it.b] - V3(gd[3], gd[4], gd[5]) * N + Ft;
               M[it.b] = M[it.b] + rB.Cross(Ft);
            }
            // on a carried surface's coordinates: the same work with the other sign
            double qf[7];
            for (int j = 0; j < nq; j++)
               qf[j] = -gd[3 + j] * N + Ft.Dot(hm.c[j].V());
            if (nq > 0)
               loadF[shape[it.s]->part] = loadF[shape[it.s]->part] + V3(gd[0], gd[1], gd[2]) * N + Ft;
            // twist and rolling resistance on the mid-step spin (relative to the other ball, or the surface)
            V3 wr = pr ? wm[it.a] - wm[it.b] : wm[it.a];
            for (int j = 0; j < nq; j++)
               wr = wr - hm.c[j].W() * (dq[3 + j] / h);
            const double am = L.A(dm);
            double kTw = 0.0, sRo = 0.0, rhoR = 0.0, wlR = 0.0, kkTw = 0.0;
            V3 wtR;
            M3 RmS;
            bool rollS = false;
            e.ra1 = RubberContact::RollAngle();
            if (set.rollModel == 1)
            {
               // the rolling and twist stick springs (roll_model 1): the patch's tilting and torsion springs with a damper side by side, in series with rolling
               // (twisting) at the limit, taken as the shear spring is (patch rule, midpoint torque, stick or roll on the total trial torque)
               RubberContact::RollAngle r0 = it.touching ? it.ra : RubberContact::RollAngle();
               {
                  const double l0 = r0.th.Length();
                  const V3 tp = r0.th - nm * r0.th.Dot(nm);
                  const double lp = tp.Length();
                  r0.th = lp > 1e-300 ? tp * (l0 / lp) : V3();
               }
               const double kr0 = L.Kr(d0), kr1 = L.Kr(d1), kw0 = L.Ktw(d0), kw1 = L.Ktw(d1);
               const double held = 0.5 * kr0 * r0.th.Dot(r0.th) + 0.5 * kw0 * r0.tw * r0.tw;
               if (kr1 > kr0)
                  r0.th = r0.th * (kr0 / kr1);
               if (kw1 > kw0)
                  r0.tw *= kw0 / kw1;
               e.roll = held - 0.5 * kr1 * r0.th.Dot(r0.th) - 0.5 * kw1 * r0.tw * r0.tw;   // the patch's change: a loss >= 0
               const double Ir = L.mRed * A.R * A.R;   // with 1.4: the ball's inertia about its contact point; with 0.4: about its centre (the contact's mass)
               if (L.twist && kw1 > 0.0 && N > 0.0 && am > 0.0)
               {
                  const double wn = wr.Dot(nm);
                  const double cw = L.shearTan * std::sqrt(kw1 * 0.4 * Ir);
                  const double twTr = r0.tw + wn * h;
                  double Mw = 0.5 * kw1 * (r0.tw + twTr) + cw * wn;
                  const double limW = 3.0 * kPi / 16.0 * L.mu * am * N;
                  double tw1 = twTr;
                  kkTw = 0.5 * kw1 * h + cw;   // the resisting torque's slope in wn while it sticks; 0 at the limit
                  if (std::fabs(Mw) > limW)
                  {
                     Mw = Mw > 0.0 ? limW : -limW;
                     const double cs = 0.5 * kw1 + cw / h;
                     tw1 = (Mw - r0.tw * (0.5 * kw1 - cw / h)) / cs;
                     e.twist = Mw * (twTr - tw1) + cw * (tw1 - r0.tw) * (tw1 - r0.tw) / h;
                     kkTw = 0.0;
                  }
                  else
                     e.twist = cw * wn * wn * h;
                  e.ra1.tw = tw1;
                  kTw = kkTw > 0.0 ? 1.0 : 0.0;   // (the Jacobian's gate)
                  M[it.a] = M[it.a] - nm * Mw;
                  if (pr)
                     M[it.b] = M[it.b] + nm * Mw;
                  for (int j = 0; j < nq; j++)
                     qf[j] += Mw * nm.Dot(hm.c[j].W());
               }
               if ((L.tanDelta > 0.0 || L.rollingMu > 0.0) && N > 0.0 && am > 0.0 && kr1 > 0.0)
               {
                  const V3 wt = wr - nm * wr.Dot(nm);
                  const double pa = N / (kPi * am * am);
                  const double muRr = std::max(std::min(2.34 * (pa / L.E) * L.tanDelta, RubberContact::kRollingMuMax), L.rollingMu);
                  const double lim = muRr * A.R * N;
                  const double cr = L.shearTan * std::sqrt(kr1 * 1.4 * Ir);
                  const V3 thTr = r0.th + wt * h;
                  V3 Mr = (r0.th + thTr) * (0.5 * kr1) + wt * cr;
                  V3 th1 = thTr;
                  const double ml = Mr.Length(), sl = 0.5 * kr1 * h + cr;
                  const M3 Pt = M3::Tangent(nm);
                  if (ml > lim)
                  {
                     // rolls: the torque at the limit along the trial; the roll thTr - th1 lies along it, its loss Mr . roll >= 0
                     Mr = Mr * (lim / ml);
                     const double cs = 0.5 * kr1 + cr / h;
                     th1 = (Mr - r0.th * (0.5 * kr1 - cr / h)) * (1.0 / cs);
                     const V3 dth = th1 - r0.th;
                     e.roll += Mr.Dot(thTr - th1) + cr * dth.Dot(dth) / h;
                     const V3 mh = Mr * (1.0 / lim);
                     const M3 Q = M3::Outer(mh, mh);
                     for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                           RmS.a[i][j] = lim / ml * sl * (Pt.a[i][j] - Q.a[i][j]);
                  }
                  else
                  {
                     e.roll += cr * wt.Dot(wt) * h;
                     for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                           RmS.a[i][j] = sl * Pt.a[i][j];
                  }
                  e.ra1.th = th1;
                  rollS = true;
                  M[it.a] = M[it.a] - Mr;
                  if (pr)
                     M[it.b] = M[it.b] + Mr;
                  for (int j = 0; j < nq; j++)
                     qf[j] += Mr.Dot(hm.c[j].W());
               }
            }
            else if (L.twist && N > 0.0 && am > 0.0)
            {
               const double wn = wr.Dot(nm);
               kTw = 3.0 * kPi / 16.0 * L.mu * am * N;
               const double Mtw = -kTw * std::tanh(wn / 2.0);
               M[it.a] = M[it.a] + nm * Mtw;
               if (pr)
                  M[it.b] = M[it.b] - nm * Mtw;
               for (int j = 0; j < nq; j++)
                  qf[j] -= Mtw * nm.Dot(hm.c[j].W());
               e.twist = -Mtw * wn * h;
            }
            if (set.rollModel != 1 && (L.tanDelta > 0.0 || L.rollingMu > 0.0) && N > 0.0 && am > 0.0)
            {
               const V3 wt = wr - nm * wr.Dot(nm);
               const double wl = wt.Length();
               const double pa = N / (kPi * am * am);
               const double muRr = std::max(std::min(2.34 * (pa / L.E) * L.tanDelta, RubberContact::kRollingMuMax), L.rollingMu);
               const double rho = muRr * A.R * N;
               sRo = wl > 1e-9 ? rho * std::tanh(wl / set.rollSmooth) / wl : rho / set.rollSmooth;
               rhoR = rho;
               wlR = wl;
               wtR = wt;
               const V3 Mr = wt * -sRo;
               M[it.a] = M[it.a] + Mr;
               if (pr)
                  M[it.b] = M[it.b] - Mr;
               for (int j = 0; j < nq; j++)
                  qf[j] -= Mr.Dot(hm.c[j].W());
               e.roll = sRo * wt.Dot(wt) * h;
            }
            for (int j = 0; j < nq; j++)
               qF[qc0 + hm.c[j].coord] += qf[j];

            if (J)
            {
               // normal: h dfN/dd1 g g^T on the centres; the chord's exact slope (F(d1) - chord) / (d1 - d0) (K / 2 for a vanishing change), the damping's
               // c / h and its change with the end depth (Hunt-Crossley 1.5 alpha x the chord's slope; Kuwabara-Kono gamma a' at dm, which moves half as far as d1)
               double kc;
               if (std::fabs(d1 - d0) > 1e-13)
                  kc = (L.F(d1) - fe) / (d1 - d0);
               else
                  kc = 0.5 * L.K(std::max(std::max(dm, d1), 0.0));
               const double cDp = L.gammaKK > 0.0 ? (dm > 0.0 ? 0.5 * L.gammaKK * 0.5 * std::sqrt(L.reff / dm) : 0.0) : 1.5 * L.alpha * kc;
               const double kn = fe + fd <= 0.0 ? 0.0 : L.gammaKK > 0.0 ? h * kc + kkJac : h * kc + cD + cDp * (d1 - d0);
               int ix[kGradMax];   // each gradient entry's unknown; -1 = a given move
               for (int p = 0; p < na; p++)
                  ix[p] = (p < 3 ? 6 * it.a : 6 * it.b) + p % 3;
               for (int j = 0; j < nq; j++)
                  ix[3 + j] = cUnk[qc0 + hm.c[j].coord];
               for (int p = 0; p < na + nq; p++)
                  if (ix[p] >= 0)
                     for (int q = 0; q < na + nq; q++)
                        if (ix[q] >= 0)
                           J[ix[p]][ix[q]] += kn * gm[p] * gm[q];
               // shear: the tangential force Ft changes with the slip rate ut (dFt/dut = Cm, in the contact plane) and with the depth d1 (dFt/dd1 = fd1).
               // ut = sp (dc_p / h + 1/2 Skew(r_p)^T dw_p) projected, d1 changes by -sq n . dc_q; on ball p the force is -sp Ft at r_p. So, for the residual of p:
               //    d/d dc_q: sp sq [Cm; Skew(r_p) Cm] - h sp sq [fd1; Skew(r_p) fd1] n^T,   d/d dw_q: 1/2 h sp sq [Cm; Skew(r_p) Cm] Skew(r_q)^T
               // The spring's force is k1 xi0 + k1 ut h / 2, and a growing patch holds k1 xi0 (= k0 xi0 before the hold): g = k1'(d1) ((k1 > k0 ? 0 : xi0) + ut h / 2) is its change with d1.
               // Sticking: Cm = (k1 h / 2 + ct) Pt, fd1 = g + k1' shearTan sqrt(mt / k1) ut / 2 (ct grows with k1 too).
               // Sliding (Ft = lim f^, f^ along the total trial force T = k1 xi0 + (k1 h / 2 + ct) ut): Cm = (lim / |T|) (k1 h / 2 + ct) (Pt - f^ f^^T),
               // fd1 = mu N'(d1) f^ + (lim / |T|) (Pt - f^ f^^T) (g + gc), N' = kn / h, the normal row's own slope.
               const M3 Pt = M3::Tangent(nm);
               const int nbod = pr ? 2 : 1;
               const int bi[2] = { it.a, it.b };
               const V3 rr[2] = { rA, rB };
               const double sg[2] = { 1.0, -1.0 };
               if (k1 > 0.0)
               {
                  M3 Cm;
                  V3 fd1;
                  const double kp = L.KtSlope(d1);
                  const V3 g = ((k1 > k0 ? V3() : xi0) + ut * (0.5 * h)) * kp;
                  const V3 gc = ut * (0.5 * L.shearTan * std::sqrt(L.mt / k1) * kp);
                  const double C = 0.5 * k1 * h + ct;
                  if (!slipping)
                  {
                     for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                           Cm.a[i][j] = C * Pt.a[i][j];
                     fd1 = g + gc;
                  }
                  else if (flTr > 0.0)
                  {
                     const V3 fh = FsTr * (1.0 / flTr);
                     const M3 Q = M3::Outer(fh, fh);
                     const double r = limS / flTr;
                     M3 Pf;
                     for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                        {
                           Pf.a[i][j] = Pt.a[i][j] - Q.a[i][j];
                           Cm.a[i][j] = r * C * Pf.a[i][j];
                        }
                     const double dN = kn / h;
                     const V3 gt = g + gc;
                     const V3 Pg(Pf.a[0][0] * gt.x + Pf.a[0][1] * gt.y + Pf.a[0][2] * gt.z, Pf.a[1][0] * gt.x + Pf.a[1][1] * gt.y + Pf.a[1][2] * gt.z,
                                 Pf.a[2][0] * gt.x + Pf.a[2][1] * gt.y + Pf.a[2][2] * gt.z);
                     fd1 = fh * (L.mu * dN) + Pg * r;
                  }
                  const M3 Dn = M3::Outer(fd1, nm);
                  for (int p = 0; p < nbod; p++)
                     for (int q = 0; q < nbod; q++)
                     {
                        const double ss = sg[p] * sg[q];
                        const M3 Sp = M3::Skew(rr[p]), SqT = M3::Skew(rr[q]).T();
                        AddBlock(J, 6 * bi[p], 6 * bi[q], Cm, ss);
                        AddBlock(J, 6 * bi[p], 6 * bi[q], Dn, -h * ss);
                        AddBlock(J, 6 * bi[p], 6 * bi[q] + 3, Cm * SqT, 0.5 * h * ss);
                        AddBlock(J, 6 * bi[p] + 3, 6 * bi[q], Sp * Cm, ss);
                        AddBlock(J, 6 * bi[p] + 3, 6 * bi[q], Sp * Dn, -h * ss);
                        AddBlock(J, 6 * bi[p] + 3, 6 * bi[q] + 3, Sp * Cm * SqT, 0.5 * h * ss);
                     }
                  // a carried surface's coordinates: ut has -v dq / h, d1 has +slope dq; the force on q is +Ft . v
                  for (int i = 0; i < nq; i++)
                  {
                     const int qi = cUnk[qc0 + hm.c[i].coord];
                     const V3 vi = hm.c[i].V();
                     if (qi >= 0)   // its row against the ball: -vi^T Cm + h (vi . fd1) n^T; against the spin: -1/2 h vi^T Cm Skew(r)^T
                     {
                        const V3 cTv = MulV(Cm.T(), vi);
                        const V3 rc = cTv * -1.0 + nm * (h * vi.Dot(fd1));
                        const V3 rw = MulV(M3::Skew(rA), cTv) * (-0.5 * h);
                        const int ia = 6 * it.a;
                        J[qi][ia] += rc.x; J[qi][ia + 1] += rc.y; J[qi][ia + 2] += rc.z;
                        J[qi][ia + 3] += rw.x; J[qi][ia + 4] += rw.y; J[qi][ia + 5] += rw.z;
                     }
                     const double li = hm.c[i].slope * hm.grad;
                     if (qi >= 0)   // the ball's rows against it: h dFt / dq = -Cm v + h slope fd1
                     {
                        const V3 dFq = MulV(Cm, vi) * -1.0 + fd1 * (h * li);
                        const V3 dMq = MulV(M3::Skew(rA), dFq);
                        const int ia = 6 * it.a;
                        J[ia][qi] += dFq.x; J[ia + 1][qi] += dFq.y; J[ia + 2][qi] += dFq.z;
                        J[ia + 3][qi] += dMq.x; J[ia + 4][qi] += dMq.y; J[ia + 5][qi] += dMq.z;
                     }
                     for (int j = 0; j < nq && qi >= 0; j++)
                     {
                        const int qj = cUnk[qc0 + hm.c[j].coord];
                        if (qj >= 0)
                           J[qi][qj] += vi.Dot(MulV(Cm, hm.c[j].V())) - h * vi.Dot(fd1) * hm.c[j].slope * hm.grad;
                     }
                  }
               }
               // twist and rolling on the spins
               if (kTw > 0.0)
               {
                  const double wn = wr.Dot(nm), sh = 1.0 / std::cosh(wn / 2.0);
                  const M3 nn = M3::Outer(nm, nm);
                  for (int p = 0; p < nbod; p++)
                     for (int q = 0; q < nbod; q++)
                        AddBlock(J, 6 * bi[p] + 3, 6 * bi[q] + 3, nn, set.rollModel == 1 ? h * 0.5 * kkTw * sg[p] * sg[q] : h * kTw * 0.25 * sh * sh * sg[p] * sg[q]);
                  const double kk = set.rollModel == 1 ? kkTw : kTw * 0.5 * sh * sh;   // dMtw / d(wn) = -kk; d(wn) / d(dq) = -nw / h
                  const int iw = 6 * it.a + 3;
                  for (int i = 0; i < nq; i++)
                  {
                     const int qi = cUnk[qc0 + hm.c[i].coord];
                     if (qi < 0)
                        continue;
                     const double nwi = nm.Dot(hm.c[i].W());
                     J[iw][qi] += -kk * nwi * nm.x; J[iw + 1][qi] += -kk * nwi * nm.y; J[iw + 2][qi] += -kk * nwi * nm.z;
                     J[qi][iw] += -0.5 * h * kk * nwi * nm.x; J[qi][iw + 1] += -0.5 * h * kk * nwi * nm.y; J[qi][iw + 2] += -0.5 * h * kk * nwi * nm.z;
                     for (int j = 0; j < nq; j++)
                     {
                        const int qj = cUnk[qc0 + hm.c[j].coord];
                        if (qj >= 0)
                           J[qi][qj] += kk * nwi * nm.Dot(hm.c[j].W());
                     }
                  }
               }
               if (sRo > 0.0 || rollS)
               {
                  // the rolling torque -rho tanh(|wt| / s) wt^ : its derivative along wt^ is rho sech^2 / s, across it rho tanh / |wt| (sRo); roll_model 1: the spring's (RmS)
                  M3 Rm;
                  if (rollS)
                     Rm = RmS;
                  else if (wlR > 1e-9)
                  {
                     const V3 wh = wtR * (1.0 / wlR);
                     const M3 Q = M3::Outer(wh, wh);
                     const double sh = 1.0 / std::cosh(wlR / set.rollSmooth);
                     const double along = rhoR / set.rollSmooth * sh * sh;
                     for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                           Rm.a[i][j] = along * Q.a[i][j] + sRo * (Pt.a[i][j] - Q.a[i][j]);
                  }
                  else
                     for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                           Rm.a[i][j] = sRo * Pt.a[i][j];
                  for (int p = 0; p < nbod; p++)
                     for (int q = 0; q < nbod; q++)
                        AddBlock(J, 6 * bi[p] + 3, 6 * bi[q] + 3, Rm, h * 0.5 * sg[p] * sg[q]);
                  const int iw = 6 * it.a + 3;
                  for (int i = 0; i < nq; i++)
                  {
                     const int qi = cUnk[qc0 + hm.c[i].coord];
                     if (qi < 0)
                        continue;
                     const V3 Re = MulV(Rm, hm.c[i].W());
                     J[iw][qi] += -Re.x; J[iw + 1][qi] += -Re.y; J[iw + 2][qi] += -Re.z;
                     J[qi][iw] += -0.5 * h * Re.x; J[qi][iw + 1] += -0.5 * h * Re.y; J[qi][iw + 2] += -0.5 * h * Re.z;
                     for (int j = 0; j < nq; j++)
                     {
                        const int qj = cUnk[qc0 + hm.c[j].coord];
                        if (qj >= 0)
                           J[qi][qj] += Re.Dot(hm.c[j].W());
                     }
                  }
               }
            }
         }
         if (carried)
            for (int p = 0; p < np; p++)
               if (part[p]->NumCoords() > 0)
                  part[p]->Loads((const double *)qF + coordOff[p], loadF[p]);
         for (int p = 0; p < np; p++)
            if (!part[p]->Eval(xv, partOff[p], h, F, M, Res, J))
            {
               partEvalFailed = true;
               return false;
            }
         if (carried)
            for (int c = 0; c < nCoord; c++)
            {
               if (cUnk[c] >= 0)
                  Res[cUnk[c]] -= h * qF[c];
               else
                  driveW -= qF[c] * qEnd[c];   // a given move does this work on the balls: a drive
            }
         for (int i = 0; i < nb; i++)
         {
            const Ball &b = ball[i];
            const V3 rc = dc[i] * (2.0 * b.m / h) - b.v * (2.0 * b.m) - F[i] * h;
            const V3 rw = dw[i] * b.I() - M[i] * h;
            Res[6 * i] = rc.x; Res[6 * i + 1] = rc.y; Res[6 * i + 2] = rc.z;
            Res[6 * i + 3] = rw.x; Res[6 * i + 4] = rw.y; Res[6 * i + 5] = rw.z;
            if (J)
               for (int a = 0; a < 3; a++)
               {
                  J[6 * i + a][6 * i + a] += 2.0 * b.m / h;
                  J[6 * i + 3 + a][6 * i + 3 + a] += b.I();
               }
         }
         for (int a = 0; a < nU; a++)
            if (!std::isfinite(Res[a]))
               return false;
         return true;
      }

      // Settings::jacMode: the residual's Jacobian at x by central differences, into jacFd (diagnostic; overwrites ie[], so an Eval at x must follow)
      void NumericJacobian(double h)
      {
         const int nU = nUnk;
         double *const xp = xpS, *const rp = rpS, *const rm = rmS;
         const Jac jacFd = JacFd();
         for (int c = 0; c < nU; c++)
         {
            for (int a = 0; a < nU; a++)
               xp[a] = x[a];
            const bool spin = c < 6 * nb && c % 6 >= 3;
            const double e = spin ? 1e-6 * (std::fabs(x[c]) + 1.0) : 1e-7 * (std::fabs(x[c]) + 1e-6);
            xp[c] = x[c] + e;
            const bool okP = Eval(h, xp, rp, nullptr, false);
            xp[c] = x[c] - e;
            const bool okM = Eval(h, xp, rm, nullptr, false);
            for (int r = 0; r < nU; r++)
               jacFd[r][c] = okP && okM ? (rp[r] - rm[r]) / (2.0 * e) : 0.0;
         }
      }

      // A step no solve can do: the implicit equation has no root (the energy has a step in it within the move, so the force that closes the energy grows as e / |dq| however short the step), or a
      // path that will not form. The same state would fail again for ever, so time goes on with the only law that needs no solve, the energy: every body takes the move the solve would have started
      // from (its speed and the force at the solve's first guess, one explicit round, stable within the step rule's omega h) with the path's elements kept, and whatever that costs in stored energy and height
      // is paid by all the speeds together (the discrete gradient's correction: the change acts along the motion). Energy that falls is booked as a loss, so nothing can be gained. What the speeds cannot pay
      // is taken shorter; speeds too small to pay for any move turn back: no move, no change in energy. The next step starts from a new state and goes back to the solve.
      double EnergyStep(const double hTry)
      {
         const double E0 = Energy();
         Ball *const b0 = esB0;
         Item *const i0 = esI0;
         V3 *const a0 = esA0, *const al0 = esAl0;
         const Books bk0 = books;
         for (int i = 0; i < nb; i++)
         {
            b0[i] = ball[i];
            a0[i] = aLast[i];
            al0[i] = alphaLast[i];
         }
         for (int k = 0; k < ni; k++)
            i0[k] = item[k];
         for (int p = 0; p < np; p++)
            part[p]->SaveState();
         const auto Restore = [&]() {
            for (int i = 0; i < nb; i++)
            {
               ball[i] = b0[i];
               aLast[i] = a0[i];
               alphaLast[i] = al0[i];
            }
            for (int k = 0; k < ni; k++)
               item[k] = i0[k];
            for (int p = 0; p < np; p++)
               part[p]->RestoreState();
            books = bk0;
            stamp++;
         };
         const auto Done = [&]() {
            for (int p = 0; p < np; p++)
               part[p]->AfterEnergyStep();
         };
         double *const Res = esRes;
         double h = hTry;
         for (int tries = 0; tries < 12; tries++, h *= 0.5)
         {
            partEvalFailed = false;
            const bool begun = BeginAttempt(h);
            for (int p = 0; p < np && begun; p++)
               part[p]->HoldArcs();
            bool formed = begun && Eval(h, x, Res, nullptr, trialCached);
            if (formed)   // one explicit round: the force at the guess, the bodies' own masses in place of Newton's matrix
            {
               for (int i = 0; i < nb; i++)
                  for (int a = 0; a < 3; a++)
                  {
                     x[6 * i + a] -= Res[6 * i + a] * h / (2.0 * ball[i].m);
                     x[6 * i + 3 + a] -= Res[6 * i + 3 + a] / ball[i].I();
                  }
               for (int p = 0; p < np; p++)
                  part[p]->ExplicitMove(Res, partOff[p], h, x);
               for (int i = 0; i < nb && formed; i++)
                  formed = x[6 * i] * x[6 * i] + x[6 * i + 1] * x[6 * i + 1] + x[6 * i + 2] * x[6 * i + 2] <= set.maxMoveM * set.maxMoveM;
               formed = formed && InRange(x) && StopExcess() <= 0.0 && Eval(h, x, Res, nullptr, trialCached);
            }
            if (formed)
            {
               Accept(h, false);
               const double ex = Energy() - E0;
               double K = 0.0;
               for (int i = 0; i < nb; i++)
                  K += 0.5 * ball[i].m * ball[i].v.Dot(ball[i].v) + 0.5 * ball[i].I() * ball[i].w.Dot(ball[i].w);
               for (int p = 0; p < np; p++)
                  K += part[p]->KineticJ();
               if (!(ex < K) && tries < 11)   // it cannot be paid: shorter (the path keeps its elements, so the energy has no step in it now)
               {
                  Restore();
                  continue;
               }
               if (!(ex < K))   // the speeds are too small to pay for any move
               {
                  Restore();
                  for (int i = 0; i < nb; i++)
                     ball[i].v = ball[i].v * -1.0;
                  for (int p = 0; p < np; p++)
                     part[p]->ScaleSpeeds(-1.0);
                  Done();
                  books.steps++;
                  books.energySteps++;
                  books.energyBounces++;
                  return hTry;
               }
               if (ex > 0.0)
               {
                  const double s = std::sqrt(1.0 - ex / K);
                  for (int i = 0; i < nb; i++)
                  {
                     ball[i].v = ball[i].v * s;
                     ball[i].w = ball[i].w * s;
                  }
                  for (int p = 0; p < np; p++)
                     part[p]->ScaleSpeeds(s);
               }
               if (Energy() - E0 > 1e-9)   // (paid only to rounding: a change too large for the speeds to close)
               {
                  Restore();
                  if (tries < 11)
                     continue;
                  for (int i = 0; i < nb; i++)
                     ball[i].v = ball[i].v * -1.0;
                  for (int p = 0; p < np; p++)
                     part[p]->ScaleSpeeds(-1.0);
                  Done();
                  books.steps++;
                  books.energySteps++;
                  books.energyBounces++;
                  return hTry;
               }
               const Books after = books;
               books = bk0;
               books.steps = after.steps;
               books.hMinUsed = after.hMinUsed;
               books.maxDepthM = after.maxDepthM;
               books.through = after.through;
               books.energySteps++;
               books.partLossJ += E0 - Energy();   // (only energy that has gone: what the speeds paid is closed to rounding)
               Done();
               stamp++;
               return h;
            }
            Restore();
         }
         // no move can be formed at any length (a force that cannot be formed at all): the balls go on by their speeds alone, the path follows them, and the same closure applies; what cannot be paid turns back
         Restore();
         for (int i = 0; i < nb; i++)
            ball[i].c = ball[i].c + ball[i].v * hTry;
         for (int p = 0; p < np; p++)
            part[p]->BallsMoved(ball, nb);
         stamp++;
         {
            const double ex = Energy() - E0;
            double K = 0.0;
            for (int i = 0; i < nb; i++)
               K += 0.5 * ball[i].m * ball[i].v.Dot(ball[i].v) + 0.5 * ball[i].I() * ball[i].w.Dot(ball[i].w);
            for (int p = 0; p < np; p++)
               K += part[p]->KineticJ();
            if (!(ex < K))
            {
               Restore();
               for (int i = 0; i < nb; i++)
                  ball[i].v = ball[i].v * -1.0;
               for (int p = 0; p < np; p++)
                  part[p]->ScaleSpeeds(-1.0);
               books.energyBounces++;
            }
            else if (ex > 0.0)
            {
               const double sf = std::sqrt(1.0 - ex / K);
               for (int i = 0; i < nb; i++)
               {
                  ball[i].v = ball[i].v * sf;
                  ball[i].w = ball[i].w * sf;
               }
               for (int p = 0; p < np; p++)
                  part[p]->ScaleSpeeds(sf);
            }
            books.partLossJ += E0 - Energy();
         }
         Done();
         books.steps++;
         books.energySteps++;
         books.energyFree++;
         stamp++;
         return hTry;
      }

      // The start of an attempt: the law of each touch, the first guess
      bool BeginAttempt(double h)
      {
         const int nU = nUnk = NumUnk();
         for (int k = 0; k < ni; k++)
         {
            Hit h0;
            if (nowStamp == stamp)
               h0 = hNow[k];
            else
               Probe(item[k], ball, h0);
            d0s[k] = h0.delta;
            farNow[k] = !farOff && OutOfReach(item[k], h0.delta, h, h0);
            if (farNow[k])
               continue;
            law[k] = item[k];
            if (!item[k].touching)
               TakeLaw(law[k], h0);
         }
         for (int i = 0; nf > 0 && i < nb; i++)
            fieldE0[i] = ball[i].fieldScale != 0.0 ? FieldEnergy(ball[i]) : 0.0;
         for (int i = 0; i < nb; i++)
         {
            // the first guess: the last step's accelerations carried on, or in the stop search this step's own from its last solve (a guess only: the solve's
            // answer does not depend on it)
            const V3 &ag = warmOn ? warmA[i] : aLast[i], &alg = warmOn ? warmAl[i] : alphaLast[i];
            const V3 g0 = (ball[i].v + ag * (0.5 * h)) * h, w0 = alg * h;
            x[6 * i] = g0.x; x[6 * i + 1] = g0.y; x[6 * i + 2] = g0.z;
            x[6 * i + 3] = w0.x; x[6 * i + 4] = w0.y; x[6 * i + 5] = w0.z;
         }
         for (int p = 0; p < np; p++)
            if (!part[p]->BeginStep(ball, nb, h, x, partOff[p]))
            {
               failWhy = 'P';
               return false;
            }
         CoordsNow();
         for (int it = 0; !InRange(x); it++)   // a first guess past a layer's range is shortened
         {
            if (it == 60)
            {
               failWhy = 'G';
               return false;
            }
            for (int a = 0; a < nU; a++)
               x[a] *= 0.5;
         }
         return true;
      }

      // Solve one step of h from the state now (the end values stay in x and ie[])
      bool Solve(double h)
      {
         farOff = false;
         if (!SolveOnce(h))
            return false;
         if (!FarStillFar(x))
         {
            farOff = true;
            const bool ok = SolveOnce(h);
            farOff = false;
            return ok;
         }
         return true;
      }
      bool farOff = false;
      bool SolveOnce(double h)
      {
         if (!BeginAttempt(h))
            return false;
         const int nU = nUnk;
         const Jac J(jStore, jStride);
         double *const Res = resS, *const stp = stpS;
         bool converged = false;
         int round = 1;
         char *const prevSlip = prevSlipS, *const prevOn = prevOnS;
         double *const xPrev = xPrevS, *const stpPrev = stpPrevS;
         double scalePrev = 0.0, rnPrev = INFINITY;
         int back = 0;
         bool gripMoved = false, endValuesHere = false;
         double lastCh = INFINITY;
         int stallCount = 0, sinceBest = 0;
         double stallMin = 0.0, bestMove = INFINITY;
         for (; round <= set.maxRounds; round++)
         {
            if (set.jacMode != 0)
               NumericJacobian(h);   // into jacFd (the Eval calls it makes are undone by the Eval below)
            partEvalFailed = false;
            if (set.jacMode == 0 && round > 1 && lastCh < set.tolCheckMove && !gripMoved && Eval(h, x, Res, nullptr, trialCached))
            {
               double rc = 0.0;
               for (int i = 0; i < nb; i++)
               {
                  const double ic = 1.0 / ball[i].m, iw = ball[i].R / ball[i].I();
                  for (int a = 0; a < 3; a++)
                     rc += Res[6 * i + a] * Res[6 * i + a] * ic * ic + Res[6 * i + 3 + a] * Res[6 * i + 3 + a] * iw * iw;
               }
               rc += CoordRes2(Res);
               if (rc <= set.tolRes * set.tolRes)
               {
                  converged = true;
                  endValuesHere = true;
                  break;
               }
            }
            if (!Eval(h, x, Res, J, trialCached))
            {
               failWhy = partEvalFailed ? 'B' : 'N';
               lastRounds = round;
               return false;
            }
            if (set.jacMode != 0)
            {
               for (int r = 0; r < nU; r++)
                  for (int c = 0; c < nU; c++)
                  {
                     JacAna()[r][c] = J[r][c];
                     if (set.jacMode == 2)
                        J[r][c] = JacFd()[r][c];
                  }
               jacRound = round;
            }
            // Line search: the residual in speed units (momentum / m, spin / (I / R)); a Newton move that did not reduce it is taken back and retaken half as
            // long (a solve can otherwise hop between two points where the law has a kink, e.g. the normal force's clamp at zero as a ball leaves).
            double rn = 0.0;
            for (int i = 0; i < nb; i++)
            {
               const double ic = 1.0 / ball[i].m, iw = ball[i].R / ball[i].I();
               for (int a = 0; a < 3; a++)
                  rn += Res[6 * i + a] * Res[6 * i + a] * ic * ic + Res[6 * i + 3 + a] * Res[6 * i + 3 + a] * iw * iw;
            }
            rn += CoordRes2(Res);
            rn = std::sqrt(rn);
            // (a part's grip changed since the last round: the residual moved for that too, so no comparison this round)
            if (round > 1 && !gripMoved && rn > rnPrev * (1.0 - 1e-4) && rn > set.tolV && back < 30)
            {
               bool ok = false;
               while (!ok && scalePrev > 1e-9)
               {
                  scalePrev *= 0.5;
                  for (int a = 0; a < nU; a++)
                     x[a] = xPrev[a] + stpPrev[a] * scalePrev;
                  ok = InRange(x);
               }
               back++;
               if (ok)
               {
                  if (round <= kRounds)
                  {
                     RoundDiag rb;
                     rb.scale = scalePrev;
                     rb.back = 1;
                     roundDiag[round - 1] = rb;
                  }
                  continue;
               }
               for (int a = 0; a < nU; a++)   // nothing shorter is in range: go on from where the move started
                  x[a] = xPrev[a];
               InRange(x);
            }
            back = 0;
            rnPrev = rn;
            RoundDiag rd;
            for (int k = 0; k < ni; k++)
            {
               const bool on = ie[k].d1 > 0.0;
               rd.nActive += ie[k].active ? 1 : 0;
               rd.nSlip += ie[k].slipping ? 1 : 0;
               if (round > 1)
               {
                  rd.nSlipFlip += ie[k].slipping != prevSlip[k] ? 1 : 0;
                  rd.nOnFlip += on != prevOn[k] ? 1 : 0;
               }
               prevSlip[k] = ie[k].slipping;
               prevOn[k] = on;
            }
            for (int a = 0; a < nU; a++)
               stp[a] = -Res[a];
            if (nU > 0 && !SolveLinear(J, stp, nU, idxS))
            {
               failWhy = 'L';
               lastRounds = round;
               return false;
            }
            double scale = 1.0;
            for (int i = 0; i < nb; i++)
            {
               const double sc = std::sqrt(stp[6 * i] * stp[6 * i] + stp[6 * i + 1] * stp[6 * i + 1] + stp[6 * i + 2] * stp[6 * i + 2]);
               const double sw = std::sqrt(stp[6 * i + 3] * stp[6 * i + 3] + stp[6 * i + 4] * stp[6 * i + 4] + stp[6 * i + 5] * stp[6 * i + 5]);
               if (sc > set.maxMoveM)
                  scale = std::min(scale, set.maxMoveM / sc);
               if (sw > set.maxSpinStep)
                  scale = std::min(scale, set.maxSpinStep / sw);
            }
            for (int c = 0; c < nCoord; c++)
               if (cUnk[c] >= 0 && ci[c].maxMove > 0.0 && std::fabs(stp[cUnk[c]]) > ci[c].maxMove)
                  scale = std::min(scale, ci[c].maxMove / std::fabs(stp[cUnk[c]]));
            double *const trial = trialS;
            int tries = 0;
            for (; tries < 60; tries++, scale *= 0.5)
            {
               for (int a = 0; a < nU; a++)
                  trial[a] = x[a] + stp[a] * scale;
               if (InRange(trial))
                  break;
            }
            if (tries == 60)
            {
               failWhy = 'R';
               lastRounds = round;
               return false;
            }
            bool moved = true;   // (a Newton move within the tolerance is not taken: the evaluation just made stays the step's end values, no further evaluation)
            double ch = 0.0;
            bool settled = true;
            for (int c = 0; c < nCoord; c++)
               if (cUnk[c] >= 0)
               {
                  const double rate = std::fabs(stp[cUnk[c]]) * scale / h;
                  if (ci[c].settle > 0.0)
                     settled = settled && rate <= ci[c].settle;
                  else
                     ch = std::max(ch, 2.0 * rate * ci[c].reach);
               }
            for (int i = 0; i < nb; i++)
            {
               const double sc = std::sqrt(stp[6 * i] * stp[6 * i] + stp[6 * i + 1] * stp[6 * i + 1] + stp[6 * i + 2] * stp[6 * i + 2]) * scale;
               const double sw = std::sqrt(stp[6 * i + 3] * stp[6 * i + 3] + stp[6 * i + 4] * stp[6 * i + 4] + stp[6 * i + 5] * stp[6 * i + 5]) * scale;
               ch = std::max(ch, std::max(2.0 * sc / h, sw * ball[i].R));
               rd.moveC = std::max(rd.moveC, 2.0 * sc / h);
               rd.moveW = std::max(rd.moveW, sw * ball[i].R);
            }
            moved = !(ch <= set.tolV && settled);
            lastCh = ch;
            for (int a = 0; a < nU; a++)
            {
               xPrev[a] = x[a];
               stpPrev[a] = stp[a];
               if (moved)
                  x[a] = trial[a];
            }
            scalePrev = scale;
            // the parts' grip on the new move (held through the next round's Newton move)
            double gch = 0.0;
            for (int p = 0; p < np; p++)
               gch = std::max(gch, part[p]->UpdateGrip(x, h));
            gripMoved = gch > set.tolV;
            rd.scale = scale;
            rd.grip = gch;
            if (round <= kRounds)
               roundDiag[round - 1] = rd;
            if (!moved && !gripMoved && set.tolEndRes > 0.0 && rn > set.tolEndRes)
            {
               failWhy = 'S';
               lastRounds = round;
               return false;
            }
            if (!moved && !gripMoved)
            {
               converged = true;
               endValuesHere = true;
               break;
            }
            // Solved to the forces' own noise: the moves (and grip changes) have stayed below tolNoise for stallRounds rounds without halving in that time. A
            // force cannot be resolved below its rounding; the energy left over is about m v x the move (6e-12 J at 2.5 m/s and 3e-11 m/s), far under the books' 1 uJ.
            {
               const double mv = std::max(ch, gch);
               // Not converging: the same answer as maxRounds would give (a force with a kink the rounds cannot settle on: a ball joining the cord at the very end of a step).
               if (mv < 0.5 * bestMove)
               {
                  bestMove = mv;
                  sinceBest = 0;
               }
               else if (++sinceBest >= set.noProgressRounds && mv > set.tolNoise)
                  break;
               if (settled && mv <= set.tolNoise)
               {
                  if (stallCount == 0 || mv < 0.5 * stallMin)
                  {
                     stallMin = mv;
                     stallCount = 1;
                  }
                  else if (++stallCount >= set.stallRounds)
                  {
                     if (set.tolEndRes > 0.0 && rn > set.tolEndRes)
                     {
                        failWhy = 'S';
                        lastRounds = round;
                        return false;
                     }
                     if (set.tolEndRes > 0.0)
                        for (int a = 0; a < nU; a++)
                           x[a] = xPrev[a];
                     converged = true;
                     books.noiseStops++;
                     break;
                  }
               }
               else
                  stallCount = 0;
            }
         }
         lastRounds = std::min(round, set.maxRounds);
         if (!converged)
         {
            failWhy = 'C';
            return false;
         }
         failWhy = 'E';
         if (endValuesHere)
            return true;
         return Eval(h, x, Res, nullptr, trialCached);   // the end values at the solution
      }

      void Accept(double h, bool onStop)
      {
         for (int p = 0; p < np; p++)   // before the balls move: a part reads their state at the step's start
         {
            double dJ = 0.0, lJ = 0.0;
            bool thr = false;
            part[p]->Accept(x, partOff[p], h, onStop, dJ, lJ, thr);
            books.driveJ += dJ;
            books.partLossJ += lJ;
            books.through = books.through || thr;
         }
         for (int i = 0; i < nb; i++)
         {
            Ball &b = ball[i];
            const V3 dc(x[6 * i], x[6 * i + 1], x[6 * i + 2]), dw(x[6 * i + 3], x[6 * i + 4], x[6 * i + 5]);
            const V3 v0 = b.v;
            if (b.brakeRate > 0.0)
            {
               const V3 wmid = b.w + dw * 0.5, wa = wmid - b.brakeAxis * wmid.Dot(b.brakeAxis);
               books.spinBrakeJ += b.brakeRate * b.I() * wa.Dot(wa) * h;
            }
            b.c = b.c + dc;
            b.v = dc * (2.0 / h) - b.v;
            b.w = b.w + dw;
            aLast[i] = (b.v - v0) * (1.0 / h);
            alphaLast[i] = dw * (1.0 / h);
         }
         for (int k = 0; k < ni; k++)
         {
            Item &it = item[k];
            const ItemEnd &e = ie[k];
            if (farNow[k])
               continue;
            if (e.d1 > 0.0)
            {
               if (!it.touching)
               {
                  const V3 keep = it.xi;
                  it = law[k];
                  it.xi = keep;
               }
               it.touching = true;
               it.xi = e.xi1;
               it.ra = e.ra1;
            }
            else
            {
               it.touching = false;
               it.xi = V3();
               it.ra = RubberContact::RollAngle();
            }
            if (it.dZero > 0.0)   // backing out, with nothing stored: the zero follows
               it.dZero = std::max(0.0, std::min(it.dZero, it.dZero + e.d1));
            books.normalLossJ += e.normalLoss;
            books.slipJ += e.slip;
            books.shearDampJ += e.shearDamp;
            books.twistJ += e.twist;
            books.rollJ += e.roll;
            books.maxDepthM = std::max(books.maxDepthM, e.d1);
            const double rMin = it.b >= 0 ? std::min(ball[it.a].R, ball[it.b].R) : ball[it.a].R;
            if (e.d1 >= rMin)
               books.through = true;
         }
         stamp++;
         books.driveJ += driveW;
         books.steps++;
         books.rounds += lastRounds;
         books.maxRounds = std::max(books.maxRounds, lastRounds);
         books.hMinUsed = std::min(books.hMinUsed, h);
         for (int p = 0; p < np; p++)
            part[p]->AfterStep(*this);
      }
   };

   // The one energy book of every group, with or without balls (a diagnostic, read once). VPX_GROUP_BOOK=<file>: one row per cycle: its books, and how much the energy rose over them.
   // VPX_GROUP_BOOK_STEPS=1: the rises per step too (an energy sum per step). A part's own split stays the part's.
   inline FILE *GroupBookFile()
   {
      static FILE *const f = []() -> FILE * {
         const char *const e = getenv("VPX_GROUP_BOOK");
         FILE *const o = e != nullptr && *e ? fopen(e, "w") : nullptr;
         if (o != nullptr)
            fprintf(o, "t_ms,alone,dt_us,nb,np,touch,steps,rounds,maxRounds,halvings,stopSteps,energySteps,failed,E0_mJ,E1_mJ,drive_uJ,loss_uJ,book_uJ,rises_uJ,maxRise_uJ,us,maxDepth_mm,coil_uJ,brake_uJ\n");
         return o;
      }();
      return f;
   }
   class CycleBook
   {
   public:
      explicit CycleBook(const Group &group) : g(group), f(GroupBookFile())
      {
         if (f == nullptr)
            return;
         static const bool on = []() { const char *const e = getenv("VPX_GROUP_BOOK_STEPS"); return e != nullptr && atoi(e) != 0; }();
         perStep = on;
         static const bool terms = []() { const char *const e = getenv("VPX_GROUP_BOOK_STEPS"); return e != nullptr && atoi(e) >= 2; }();
         if (terms)
            startTerms = g.Terms();
         e0 = g.Energy();
         t0 = std::chrono::steady_clock::now();
      }
      void AfterStep()
      {
         if (f == nullptr || !perStep)
            return;
         const double now = g.Energy() - e0 + g.books.Losses() - g.books.driveJ, rise = now - last;
         rises += std::max(rise, 0.0);
         maxRise = std::max(maxRise, rise);
         last = now;
      }
      void End(const int tMs, const bool alone, const double dt, const bool failed) const
      {
         if (f == nullptr)
            return;
         const double us = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count();
         const double e1 = g.Energy(), l = g.books.Losses(), d = g.books.driveJ;
         int touch = 0;
         for (int k = 0; k < g.ni; k++)
            touch += g.Touching(k) ? 1 : 0;
         fprintf(f, "%d,%d,%.3f,%d,%d,%d,%ld,%ld,%d,%ld,%ld,%ld,%d,%.9f,%.9f,%.6f,%.6f,%.6f,%.6f,%.6f,%.1f,%.5f,%.6f,%.6f\n", tMs, alone ? 1 : 0, dt * 1e6, g.nb, g.np, touch, g.books.steps, g.books.rounds, g.books.maxRounds, g.books.halvings, g.books.stopSteps, g.books.energySteps, failed ? 1 : 0,
                 e0 * 1e3, e1 * 1e3, d * 1e6, l * 1e6, (e1 - e0 + l - d) * 1e6, rises * 1e6, maxRise * 1e6, us, g.books.maxDepthM * 1e3, g.books.coilJ * 1e6, g.books.spinBrakeJ * 1e6);
         if (!startTerms.empty())   // (uJ, um, m)
            fprintf(f, "T,%d,%.3f,start%s\nT,%d,%.3f,end%s\n", tMs, dt * 1e6, startTerms.c_str(), tMs, dt * 1e6, g.Terms().c_str());
         fflush(f);
      }

   private:
      const Group &g;
      FILE *const f;
      bool perStep = false;
      double e0 = 0.0, last = 0.0, rises = 0.0, maxRise = 0.0;
      std::string startTerms;
      std::chrono::steady_clock::time_point t0;
   };

   // A part with no ball this cycle: the same solve, no balls; one function for every part. False: a step could not be solved (time is lost, no energy can appear: the next cycle goes on from here).
   // out: the group's books of this cycle, when wanted.
   inline bool AdvanceAlone(Part &p, const double dt, const double hMin, const int tMs, Books *const out = nullptr)
   {
      static Group g;
      g.nb = g.ns = g.np = 0;
      g.set.hMin = hMin;
      g.AddPart(&p);
      g.Configure();
      CycleBook book(g);
      double done = 0.0;
      const double hFix = GroupFixedStepS();
      const bool ok = g.StepsOver(dt, done, [&](const double h) { return hFix > 0.0 ? std::min(hFix, dt - done) : h; }, [&book]() { book.AfterStep(); });
      book.End(tMs, true, dt, !ok);
      if (out != nullptr)
         *out = g.books;
      return ok;
   }
}
