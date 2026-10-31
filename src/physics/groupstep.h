// license:GPLv3+
//
// The group step (docs/GroupStepDesign.txt; ProjectAssessment-2026-10-27.txt item 0 b): the balls of a touching group and their contacts moved together in one solve,
// so that energy is kept by construction. Build step 1a: balls, fixed surfaces and ball pairs (the band joins in step 1b).
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
// Moving parts (step 1b, GroupStepDesign.txt 10): a Part describes its own bodies' unknowns (after the balls'), its stored energy, losses, drives and grip, and the
// group solves it with the balls (the band cord, bandcord.h, is the first).
// Pure SI maths, no VPX types, no allocation in the step: the offline tests runs/2026-10-29-group-step and runs/2026-10-29-group-band build it.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstring>

#include "physics/rubbercontact.h"

namespace GroupStep
{
   using V3 = RubberContact::V3;
   using RubberContact::RubberParams;
   constexpr double kPi = 3.14159265358979323846;
   constexpr int kMaxBalls = 8, kMaxShapes = 32, kMaxParts = 4, kMaxPartDof = 2;
   constexpr int kMaxItems = kMaxBalls * kMaxShapes + kMaxBalls * (kMaxBalls - 1) / 2;
   constexpr int kMaxUnk = 6 * kMaxBalls + kMaxParts * kMaxPartDof;

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

   // ---- surfaces ----
   // What a surface reports for a ball centre: the signed depth (> 0 = squashed; the ball's radius plus a tube's radius minus the distance), the normal (surface -> ball
   // centre) and the Hertz effective radius. The depth must be continuous in the ball's position (the step takes exact differences of it).
   struct Hit
   {
      double delta = -1.0;
      V3 n;
      double reff = 0.0;
      double curvDist = 0.0;   // > 0: the ball turns round a centre this far away while it rolls along the surface (a post, a tube): the step's turning rule
   };
   class Shape
   {
   public:
      RubberParams rp;         // the material (rubbercontact.h)
      double layerM = 0.0;     // the rubber's thickness for the layer term, m (0 = no layer term)
      int onlyBall = -1;       // >= 0: only this ball meets it (a surface's list of near faces is one ball's)
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
      V3 fExt;               // a held outside force (a magnet, a part not yet moved over): a drive, its work counted
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
      // The step rule by motion (docs/GroupStepDesign.txt 9, 16): a contact that is touching and steady (its compression rate within vSteady, its slip rate within vSteady or at the friction limit) only
      // needs its normal stiffness resolved to omegaHRest; an impact, a slip the spring can still catch or an arrival keeps the strict rule above. vSteady 0 = every contact strict.
      // 0.05 m/s and 2.4 from the 2026-10-29 sweep: every exit as accurate as with the strict rule, resting and pinched balls 3-5x cheaper.
      double vSteady = 0.05;                   // m/s
      double omegaHRest = 2.4;
      double turnMaxRad = 2.0 * kPi / 180.0;   // a ball rolling round a post or another ball turns at most this per step
      double vFloor = 0.01;                    // m/s: the closing speed below which a contact's expected depth is taken at this speed
      double hMax = 1e-3, hMin = 1e-9;
      // the solver
      double tolV = 1e-11;                     // solved when a round changes no ball's speed (or its rim speed R w) by more than this, m/s
      double tolNoise = 1e-9;                  // or when the changes stay below this for stallRounds rounds without halving: the forces' own noise (m/s)
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
      long steps = 0, rounds = 0, halvings = 0, stopSteps = 0, noiseStops = 0;   // noiseStops: steps solved to the forces' noise (Settings::tolNoise)
      long energySteps = 0, energyBounces = 0, energyFree = 0;   // steps no solve could do: taken on energy alone, of them turned back at a step in the energy too high to pay, moved with no path (Group::EnergyStep)
      int maxRounds = 0;
      double hMinUsed = 1.0, maxDepthM = 0.0;
      bool failed = false, through = false;
      double Losses() const { return normalLossJ + slipJ + shearDampJ + twistJ + rollJ + partLossJ; }
   };

   class Group;

   // A moving part of a group (GroupStepDesign.txt piece 1; 10.1): it describes its bodies, stored energy, losses, drives and grip, and the group moves them with
   // the balls in one solve. The unknowns x: per ball the move (3) and the spin change (3), then each part's own (NumDof of them, from off). The balls' residual
   // rows are 2 m dc / h - 2 m v - h F and I dw - h M; a part adds its forces on the balls (averaged over the step) to F and M and writes its own rows.
   // Every call of a step reads the balls' state at the step's start (BeginStep's b) and the move x.
   class Part
   {
   public:
      virtual ~Part() = default;
      virtual int NumDof() const = 0;
      // the part's stored energy and its own bodies' kinetic energy at the balls' positions b (the balls' own energy is the group's)
      virtual double Energy(const Ball *b, int nb) const = 0;
      // its step choice from the state now (accuracy only; GroupStepDesign 10.3); rule: what set it, for the trace
      virtual double ChooseStep(const Ball *b, int nb, const Settings &s, int &rule) = 0;
      // a step of h from the state b: its start values, and the first guess of its own unknowns into x[off..]
      virtual bool BeginStep(const Ball *b, int nb, double h, double *x, int off) = 0;
      virtual bool InRange(const double *x, int off) const = 0;   // the move x may be tried (no layer squashed to nothing)
      // the step's forces at x (see above); with J its entries of the residual's Jacobian (the balls' rows and columns in the group's layout)
      virtual bool Eval(const double *x, int off, double h, V3 *F, V3 *M, double *Res, double (*J)[kMaxUnk]) = 0;
      // between the solve's rounds: its grip impulses on the balls from the rates at x; returns their largest change (m/s)
      virtual double UpdateGrip(const double *x, double h) = 0;
      // its share of a Newton move stp: the largest share (<= 1) that may be taken; whether its own unknowns have settled; its rows in the line search's speed units
      virtual double MoveScale(const double *stp, int off) const = 0;
      virtual bool Settled(const double *stp, int off, double scale, double h) const = 0;
      virtual double ResNorm2(const double *Res, int off) const = 0;
      // > 0: the move x carries a body past a hard stop by this much (the step is shortened to end on it; piece 1 item 7)
      virtual double StopExcess(const double *x, int off) const = 0;
      virtual bool PinOnStop(const double *, int) { return false; }   // on its hard stop already and the step would take it on past: it stays there. False: nothing changed
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
      virtual void BallsMoved(const Ball *, int) {}              // the balls were moved with no step of its own
      virtual void ExplicitMove(const double *, int, double, double *) const {}   // its own unknowns in an energy step's explicit round
   };

   // One contact: ball a on shape s, or ball a on ball b. Its law is taken when a touch starts and held while it lasts.
   struct Item
   {
      int a = -1, b = -1, s = -1;
      bool touching = false;
      V3 xi;                          // the shear spring's stretch (in the contact plane)
      double dZero = 0.0;             // the depth its squash counts from: a ball found inside has stored nothing (no energy without the work that stored it; only ever down)
      // the held law
      double kH = 0.0, cL = 0.0, hL = 0.0, reff = 0.0, ktPerA = 0.0, E = 0.0, mRed = 0.0, mt = 0.0;
      double alpha = 0.0, gammaKK = 0.0, tanDelta = 0.0, shearTan = 0.0, mu = 0.0, rollingMu = 0.0;
      bool normalLoss = false, twist = false;
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
      double V(double d) const { return d > 0.0 ? 0.4 * kH * d * d * std::sqrt(d) + (cL > 0.0 ? cL * W(d) : 0.0) : 0.0; }
      double F(double d) const
      {
         if (d <= 0.0)
            return 0.0;
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
         const double lw = cL > 0.0 ? cL * ((d1 > 0.0 ? W(d1) : 0.0) - (d0 > 0.0 ? W(d0) : 0.0)) / (d1 - d0) : 0.0;
         return hz + lw;
      }
      double A(double d) const { return d > 0.0 ? std::sqrt(reff * d) : 0.0; }   // contact patch radius
      double Kt(double d) const { return ktPerA * A(d); }
      double KtSlope(double d) const { return d > 0.0 ? 0.5 * ktPerA * std::sqrt(reff / d) : 0.0; }
      // the normal damping coefficient over a step, mid-step depth d and mean spring force fMean: force = c x compression rate. Hunt-Crossley's 1.5 alpha F is taken over the step, so its loss is
      // 1.5 alpha x the change in stored energy x the rate, the law's own integral at a steady rate; at the mid-step force alone the loss falls short of it by (omega h)^2 and the rebound comes out high
      double DampC(const double d, const double fMean) const
      {
         double c = gammaKK > 0.0 ? (d > 0.0 ? gammaKK * A(d) : 0.0) : 1.5 * alpha * fMean;
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
      int nb = 0, ns = 0, ni = 0, np = 0;
      Ball ball[kMaxBalls];
      const Shape *shape[kMaxShapes] = {};
      Item item[kMaxItems];
      Part *part[kMaxParts] = {};   // moving parts (the band); they read the balls through the step calls

      // After the balls and shapes are in: the contact list (every ball on every shape, every ball pair), with the touches already present
      void Configure()
      {
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

      double Energy() const
      {
         double e = 0.0;
         for (int i = 0; i < nb; i++)
         {
            const Ball &b = ball[i];
            e += 0.5 * b.m * b.v.Dot(b.v) + 0.5 * b.I() * b.w.Dot(b.w) - b.m * set.g.Dot(b.c);
         }
         for (int k = 0; k < ni; k++)
            if (item[k].touching)
            {
               Hit hp;
               Probe(item[k], ball, hp);
               e += item[k].V(hp.delta) + 0.5 * item[k].Kt(hp.delta) * item[k].xi.Dot(item[k].xi);
            }
         for (int p = 0; p < np; p++)
            e += part[p]->Energy(ball, nb);
         return e;
      }
      bool Touching(int k) const { return item[k].touching; }
      // The first guess starts afresh: a step taken again from a snapshot then gives the same result to the last digit
      void ClearGuess()
      {
         for (int i = 0; i < kMaxBalls; i++)
            aLast[i] = alphaLast[i] = V3();
      }
      // After Configure: a part's own contact goes on as it was kept (touching, stretch, zero); one not touching that is already inside has stored nothing, its depth is its zero
      void SetShapeItem(const int k, const bool wasTouching, const double zero, const V3 &xi)
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
         Probe(it, ball, hp);
         it.touching = hp.delta > 0.0;
         if (it.touching)
            TakeLaw(it, hp);
         else
            it.xi = V3();
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
            it.xi = V3();
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
         for (int k = 0; k < ni; k++)
         {
            const Item &it = item[k];
            Hit &hp = hNow[k];   // kept: the step's start probes the same positions
            Probe(it, ball, hp);
            const Ball &A = ball[it.a];
            const V3 vr = it.b >= 0 ? A.v - ball[it.b].v : A.v;
            const double dDot = -vr.Dot(hp.n);   // closing speed
            const double vc = std::max(dDot, set.vFloor);
            if (hp.delta <= 0.0 && -hp.delta >= vc * set.hMax + 1e-9)
               continue;   // cannot be reached within the longest step
            Item L = it;
            if (!it.touching)
               TakeLaw(L, hp);
            // the depth this approach reaches (all the closing energy into the squash): the stiffness while it lasts
            double lo = std::pow(1.25 * L.mRed * vc * vc / L.kH, 0.4);   // the Hertz spring alone; a layer only makes it shallower
            if (L.cL > 0.0)
            {
               double hi = std::min(lo, L.hL);
               lo = 0.0;
               for (int r = 0; r < 50 && hi - lo > 1e-4 * hi; r++)
               {
                  const double mid = 0.5 * (lo + hi);
                  (L.V(mid) < 0.5 * L.mRed * vc * vc ? lo : hi) = mid;
               }
            }
            const double dE = std::max(hp.delta, lo);
            const double ktE = L.Kt(dE);
            const double omN = std::sqrt(L.K(dE) / L.mRed), omT = std::sqrt(ktE / L.mt);
            // steady: touching, with its compression rate within vSteady, and the slip rate of its contact point within vSteady or at the friction limit
            bool steady = false;
            if (set.vSteady > 0.0 && hp.delta > 0.0 && std::fabs(dDot) <= set.vSteady)
            {
               const double dmp = hp.delta;
               V3 u = A.v + A.w.Cross(hp.n * (-(A.R - 0.5 * dmp)));
               if (it.b >= 0)
                  u = u - (ball[it.b].v + ball[it.b].w.Cross(hp.n * (ball[it.b].R - 0.5 * dmp)));
               const V3 ut = u - hp.n * u.Dot(hp.n);
               const double slipNow = ut.Length();
               // At the friction limit (Coulomb): a spring within mu N holds at most (mu N)^2 / (2 kt) of energy, so a slip with more than that in 1/2 mt u^2 cannot be stopped by it: the contact slides,
               // its force is mu N whatever the spring's stiffness, and the spring's time scale is not in the motion. Below that speed the spring can catch it, and the rule above resolves it again.
               steady = slipNow <= set.vSteady || (set.slideFree && slipNow * std::sqrt(ktE * L.mt) > L.mu * L.F(dE));
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
            int pr = 0;
            const double hp = part[p]->ChooseStep(ball, nb, set, pr);
            if (hp < h)
            {
               h = hp;
               dg = StepDiag();
               dg.rule = 8;
               dg.partRule = pr;
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
      double Step(double hTry)
      {
         double h = hTry;
         nHalveDiag = 0;
         for (;;)
         {
            if (Solve(h))
            {
               if (StopExcess() <= 0.0)
               {
                  Accept(h, false);
                  return h;
               }
               // the step that ends on the stop (Illinois false position on the step length), then solved again there (its end state is the one taken)
               double zero[kMaxUnk] = {};
               {
                  bool pinnedOne = false;
                  for (int p = 0; p < np; p++)
                     if (part[p]->StopExcess(x, partOff[p]) > 0.0 && part[p]->StopExcess(zero, partOff[p]) > -kOnStop)
                        pinnedOne = part[p]->PinOnStop(zero, partOff[p]) || pinnedOne;
                  if (pinnedOne)
                     continue;
               }
               double lo = 0.0, flo = std::min(StopExcess(zero), -1e-30), hi = h, fhi = StopExcess(), hBelow = 0.0;
               int side = 0;
               for (int it = 0; it < 80; it++)
               {
                  double hm = (lo * fhi - hi * flo) / (fhi - flo);
                  if (!(hm > lo && hm < hi))
                     hm = 0.5 * (lo + hi);
                  if (!Solve(hm))
                     return EnergyStep(hTry);
                  const double f = StopExcess();
                  if (f <= 0.0)
                  {
                     hBelow = hm;
                     if (f > -kOnStop)
                        break;
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
               if (hBelow <= 0.0 || !Solve(hBelow))
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
      double jacAna[kMaxUnk][kMaxUnk] = {}, jacFd[kMaxUnk][kMaxUnk] = {};
      int jacRound = 0;

   private:
      // scratch of the step being solved (no allocation)
      struct ItemEnd
      {
         bool active = false, slipping = false;
         double d0 = 0.0, d1 = 0.0;
         V3 xi1;
         double normalLoss = 0.0, slip = 0.0, shearDamp = 0.0, twist = 0.0, roll = 0.0;
      };
      Item law[kMaxItems];       // the law each contact uses in this step (held, or taken at the step's start)
      // Probes kept instead of repeated (a tube probe is most of a step's cost): the positions now (ChooseStep; valid while nowStamp == stamp), and the end and
      // mid-step positions of the move the range check accepted (the next Eval looks at the same positions)
      mutable Hit hNow[kMaxItems];
      mutable long nowStamp = -1;
      long stamp = 0;            // counts changes of the positions (Configure, Accept)
      Hit h1c[kMaxItems], hmc[kMaxItems];
      bool trialCached = false;
      double d0s[kMaxItems];
      Ball b1r[kMaxBalls], bmr[kMaxBalls];   // InRange's moved balls
      ItemEnd ie[kMaxItems];
      double x[kMaxUnk];         // the unknowns: per ball the move (3) and the spin change (3), then the parts' own
      V3 aLast[kMaxBalls], alphaLast[kMaxBalls];   // each ball's acceleration and spin's over the last step (the next step's first guess)
      double driveW = 0.0;
      int partOff[kMaxParts] = {};
      int nUnk = 0;              // unknowns of the step being solved
      bool partEvalFailed = false;

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
         trialCached = false;
         bool moved = false;
         for (int k = 0; k < ni; k++)
            if (law[k].cL > 0.0)
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
               Probe(item[k], b1, h1);
               Probe(item[k], bm, hm);
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
      static bool SolveLinear(double A[kMaxUnk][kMaxUnk], double *b, int n)
      {
         for (int col = 0; col < n; col++)
         {
            int piv = col;
            for (int rw = col + 1; rw < n; rw++)
               if (std::fabs(A[rw][col]) > std::fabs(A[piv][col]))
                  piv = rw;
            if (std::fabs(A[piv][col]) < 1e-300)
               return false;
            if (piv != col)
            {
               for (int k = 0; k < n; k++)
                  std::swap(A[col][k], A[piv][k]);
               std::swap(b[col], b[piv]);
            }
            for (int rw = col + 1; rw < n; rw++)
            {
               const double f = A[rw][col] / A[col][col];
               if (f == 0.0)
                  continue;
               for (int k = col; k < n; k++)
                  A[rw][k] -= f * A[col][k];
               b[rw] -= f * b[col];
            }
         }
         for (int rw = n - 1; rw >= 0; rw--)
         {
            double s = b[rw];
            for (int k = rw + 1; k < n; k++)
               s -= A[rw][k] * b[k];
            b[rw] = s / A[rw][rw];
         }
         return true;
      }
      static void AddBlock(double J[kMaxUnk][kMaxUnk], int r0, int c0, const M3 &m, double s)
      {
         for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
               J[r0 + i][c0 + j] += s * m.a[i][j];
      }

      // Forces of the step at the unknowns xv: the residual R (and the approximate Jacobian J when given), the contacts' end values into ie[]
      bool Eval(double h, const double *xv, double *Res, double (*J)[kMaxUnk], bool useCache)
      {
         const int nU = nUnk;
         V3 F[kMaxBalls], M[kMaxBalls], dc[kMaxBalls], dw[kMaxBalls], wm[kMaxBalls];
         Ball b1[kMaxBalls], bm[kMaxBalls];
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
         driveW = 0.0;
         for (int i = 0; i < nb; i++)
            driveW += ball[i].fExt.Dot(dc[i]);

         for (int k = 0; k < ni; k++)
         {
            const Item &it = item[k];
            const Item &L = law[k];
            ItemEnd &e = ie[k];
            e = ItemEnd();
            Hit h1, hm;
            if (useCache && L.cL > 0.0)
            {
               h1 = h1c[k];
               hm = hmc[k];
            }
            else
            {
               Probe(it, b1, h1);
               Probe(it, bm, hm);
            }
            const double d0 = d0s[k], d1 = h1.delta, dm = hm.delta;
            e.d0 = d0;
            e.d1 = d1;
            if (d0 <= 0.0 && d1 <= 0.0 && dm <= 0.0)
               continue;
            e.active = true;
            const bool pr = it.b >= 0;
            const int na = pr ? 6 : 3;
            const V3 &nm = hm.n;
            // the compression's discrete gradient over the move (ball a's centre, then ball b's)
            double gm[6] = { -nm.x, -nm.y, -nm.z, nm.x, nm.y, nm.z }, dq[6], W[6], gd[6];
            const Ball &A = ball[it.a];
            dq[0] = dc[it.a].x; dq[1] = dc[it.a].y; dq[2] = dc[it.a].z;
            W[0] = W[1] = W[2] = 1.0 / A.m;
            if (pr)
            {
               dq[3] = dc[it.b].x; dq[4] = dc[it.b].y; dq[5] = dc[it.b].z;
               W[3] = W[4] = W[5] = 1.0 / ball[it.b].m;
            }
            Gonzalez(gm, d1 - d0, noiseD, dq, W, na, gd);
            const double fe = L.Chord(d0, d1);
            const double cD = L.DampC(dm, fe);
            const double fd = cD * (d1 - d0) / h;
            const double N = std::max(0.0, fe + fd);
            e.normalLoss = cD * (d1 - d0) * (d1 - d0) / h;

            // the contact point and the mid-step slip rate
            const double dmp = std::max(dm, 0.0);
            const V3 rA = nm * (-(A.R - 0.5 * dmp));
            const V3 rB = pr ? nm * (ball[it.b].R - 0.5 * dmp) : V3();
            V3 u = dc[it.a] * (1.0 / h) + wm[it.a].Cross(rA);
            if (pr)
               u = u - (dc[it.b] * (1.0 / h) + wm[it.b].Cross(rB));
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
            V3 Fs, xi1, FsTr;
            double ct = 0.0, limS = 0.0, flTr = 0.0, dampShare = 1.0;
            bool slipping = false;
            if (k1 > 0.0)
            {
               const V3 xtr = xi0 + ut * h;
               Fs = (xi0 + xtr) * (0.5 * k1);
               FsTr = Fs;
               const double lim = L.mu * N, fl = Fs.Length();
               limS = lim;
               flTr = fl;
               if (fl > lim)
               {
                  Fs = fl > 0.0 ? Fs * (lim / fl) : V3();
                  xi1 = Fs * (2.0 / k1) - xi0;
                  slipping = true;
                  e.slipping = true;
                  e.slip += Fs.Dot(xtr - xi1);
               }
               else
               {
                  xi1 = xtr;
                  ct = L.shearTan * std::sqrt(k1 * L.mt);
                  dampShare = RubberContact::DamperShare(Fs, ut * ct, lim);
               }
            }
            e.xi1 = xi1;
            const V3 Ft = Fs + ut * (ct * dampShare);
            e.shearDamp = ct * dampShare * ut.Dot(ut) * h;
            const double fN = fe + fd;
            F[it.a] = F[it.a] - V3(gd[0], gd[1], gd[2]) * fN - Ft;
            M[it.a] = M[it.a] + rA.Cross(Ft * -1.0);
            if (pr)
            {
               F[it.b] = F[it.b] - V3(gd[3], gd[4], gd[5]) * fN + Ft;
               M[it.b] = M[it.b] + rB.Cross(Ft);
            }
            // twist and rolling resistance on the mid-step spin (relative to the other ball)
            const V3 wr = pr ? wm[it.a] - wm[it.b] : wm[it.a];
            const double am = L.A(dm);
            double kTw = 0.0, sRo = 0.0, rhoR = 0.0, wlR = 0.0;
            V3 wtR;
            if (L.twist && N > 0.0 && am > 0.0)
            {
               const double wn = wr.Dot(nm);
               kTw = 3.0 * kPi / 16.0 * L.mu * am * N;
               const double Mtw = -kTw * std::tanh(wn / 2.0);
               M[it.a] = M[it.a] + nm * Mtw;
               if (pr)
                  M[it.b] = M[it.b] - nm * Mtw;
               e.twist = -Mtw * wn * h;
            }
            if ((L.tanDelta > 0.0 || L.rollingMu > 0.0) && N > 0.0 && am > 0.0)
            {
               const V3 wt = wr - nm * wr.Dot(nm);
               const double wl = wt.Length();
               const double pa = N / (kPi * am * am);
               const double muRr = std::max(std::min(2.34 * (pa / L.E) * L.tanDelta, RubberContact::kRollingMuMax), L.rollingMu);
               const double rho = muRr * A.R * N;
               sRo = wl > 1e-9 ? rho * std::tanh(wl / 2.0) / wl : 0.5 * rho;
               rhoR = rho;
               wlR = wl;
               wtR = wt;
               const V3 Mr = wt * -sRo;
               M[it.a] = M[it.a] + Mr;
               if (pr)
                  M[it.b] = M[it.b] - Mr;
               e.roll = sRo * wt.Dot(wt) * h;
            }

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
               const double kn = h * kc + cD + cDp * (d1 - d0);
               for (int p = 0; p < na; p++)
                  for (int q = 0; q < na; q++)
                  {
                     const int rp = (p < 3 ? 6 * it.a : 6 * it.b) + p % 3, cq = (q < 3 ? 6 * it.a : 6 * it.b) + q % 3;
                     J[rp][cq] += kn * gm[p] * gm[q];
                  }
               // shear: the tangential force Ft changes with the slip rate ut (dFt/dut = Cm, in the contact plane) and with the depth d1 (dFt/dd1 = fd1).
               // ut = sp (dc_p / h + 1/2 Skew(r_p)^T dw_p) projected, d1 changes by -sq n . dc_q; on ball p the force is -sp Ft at r_p. So, for the residual of p:
               //    d/d dc_q: sp sq [Cm; Skew(r_p) Cm] - h sp sq [fd1; Skew(r_p) fd1] n^T,   d/d dw_q: 1/2 h sp sq [Cm; Skew(r_p) Cm] Skew(r_q)^T
               // The spring's force is k1 xi0 + k1 ut h / 2, and a growing patch holds k1 xi0 (= k0 xi0 before the hold): g = k1'(d1) ((k1 > k0 ? 0 : xi0) + ut h / 2) is its change with d1.
               // Sticking: Cm = (k1 h / 2 + ct) Pt, fd1 = g + k1' shearTan sqrt(mt / k1) ut / 2 (ct grows with k1 too).
               // Damper at the limit (Ft = Fs + s D, D = ct ut, |Ft| = lim): with s from that, dFt = P (dFs + s dD) + D lim dlim / (Ft . D), P = I - D Ft^T / (Ft . D).
               // Sliding (Ft = lim f^, f^ along the spring's trial force): Cm = (lim / |trial|) (k1 h / 2) (Pt - f^ f^^T), fd1 = mu N'(d1) f^ +
               // (lim / |trial|) (Pt - f^ f^^T) g, N' = K / 2 + c / h.
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
                  const double FD = Ft.Dot(ut) * ct;
                  if (!slipping && dampShare < 1.0 && FD > 1e-30)
                  {
                     const V3 D = ut * ct;
                     const V3 dd = g + gc * dampShare;
                     const double C = 0.5 * k1 * h + ct * dampShare;
                     const M3 Q = M3::Outer(D, Ft);
                     M3 P;
                     for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                           P.a[i][j] = (i == j ? 1.0 : 0.0) - Q.a[i][j] / FD;
                     const M3 PPt = P * Pt;
                     for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                           Cm.a[i][j] = C * PPt.a[i][j];
                     const double dN = 0.5 * L.K(std::max(std::max(dm, d1), 0.0)) + cD / h;
                     fd1 = dd + D * ((limS * L.mu * dN - Ft.Dot(dd)) / FD);
                  }
                  else if (!slipping)
                  {
                     const double C = 0.5 * k1 * h + ct;
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
                           Cm.a[i][j] = r * 0.5 * k1 * h * Pf.a[i][j];
                        }
                     const double dN = 0.5 * L.K(std::max(std::max(dm, d1), 0.0)) + cD / h;
                     const V3 Pg(Pf.a[0][0] * g.x + Pf.a[0][1] * g.y + Pf.a[0][2] * g.z, Pf.a[1][0] * g.x + Pf.a[1][1] * g.y + Pf.a[1][2] * g.z,
                                 Pf.a[2][0] * g.x + Pf.a[2][1] * g.y + Pf.a[2][2] * g.z);
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
               }
               // twist and rolling on the spins
               if (kTw > 0.0)
               {
                  const double wn = wr.Dot(nm), sh = 1.0 / std::cosh(wn / 2.0);
                  const M3 nn = M3::Outer(nm, nm);
                  for (int p = 0; p < nbod; p++)
                     for (int q = 0; q < nbod; q++)
                        AddBlock(J, 6 * bi[p] + 3, 6 * bi[q] + 3, nn, h * kTw * 0.25 * sh * sh * sg[p] * sg[q]);
               }
               if (sRo > 0.0)
               {
                  // the rolling torque -rho tanh(|wt| / 2) wt^ : its derivative along wt^ is rho sech^2 / 2, across it rho tanh / |wt| (sRo)
                  M3 Rm;
                  if (wlR > 1e-9)
                  {
                     const V3 wh = wtR * (1.0 / wlR);
                     const M3 Q = M3::Outer(wh, wh);
                     const double sh = 1.0 / std::cosh(wlR / 2.0), along = 0.5 * rhoR * sh * sh;
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
               }
            }
         }
         for (int p = 0; p < np; p++)
            if (!part[p]->Eval(xv, partOff[p], h, F, M, Res, J))
            {
               partEvalFailed = true;
               return false;
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
         double xp[kMaxUnk], rp[kMaxUnk], rm[kMaxUnk];
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
         static Ball b0[kMaxBalls];
         static Item i0[kMaxItems];
         V3 a0[kMaxBalls], al0[kMaxBalls];
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
         double Res[kMaxUnk];
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
               formed = InRange(x) && StopExcess() <= 0.0 && Eval(h, x, Res, nullptr, trialCached);
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
            law[k] = item[k];
            if (!item[k].touching)
               TakeLaw(law[k], h0);
         }
         for (int i = 0; i < nb; i++)
         {
            // the first guess: the last step's accelerations carried on (a guess only: the solve's answer does not depend on it)
            const V3 g0 = (ball[i].v + aLast[i] * (0.5 * h)) * h, w0 = alphaLast[i] * h;
            x[6 * i] = g0.x; x[6 * i + 1] = g0.y; x[6 * i + 2] = g0.z;
            x[6 * i + 3] = w0.x; x[6 * i + 4] = w0.y; x[6 * i + 5] = w0.z;
         }
         for (int p = 0; p < np; p++)
            if (!part[p]->BeginStep(ball, nb, h, x, partOff[p]))
            {
               failWhy = 'P';
               return false;
            }
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
         if (!BeginAttempt(h))
            return false;
         const int nU = nUnk;
         static double J[kMaxUnk][kMaxUnk];   // (scratch; the offline test and the engine step one group at a time)
         double Res[kMaxUnk], stp[kMaxUnk];
         bool converged = false;
         int round = 1;
         bool prevSlip[kMaxItems], prevOn[kMaxItems];
         double xPrev[kMaxUnk], stpPrev[kMaxUnk], scalePrev = 0.0, rnPrev = INFINITY;
         int back = 0;
         bool gripMoved = false, endValuesHere = false;
         int stallCount = 0, sinceBest = 0;
         double stallMin = 0.0, bestMove = INFINITY;
         for (; round <= set.maxRounds; round++)
         {
            if (set.jacMode != 0)
               NumericJacobian(h);   // into jacFd (the Eval calls it makes are undone by the Eval below)
            partEvalFailed = false;
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
                     jacAna[r][c] = J[r][c];
                     if (set.jacMode == 2)
                        J[r][c] = jacFd[r][c];
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
            for (int p = 0; p < np; p++)
               rn += part[p]->ResNorm2(Res, partOff[p]);
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
            if (nU > 0 && !SolveLinear(J, stp, nU))
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
            for (int p = 0; p < np; p++)
               scale = std::min(scale, part[p]->MoveScale(stp, partOff[p]));
            double trial[kMaxUnk];
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
            for (int i = 0; i < nb; i++)
            {
               const double sc = std::sqrt(stp[6 * i] * stp[6 * i] + stp[6 * i + 1] * stp[6 * i + 1] + stp[6 * i + 2] * stp[6 * i + 2]) * scale;
               const double sw = std::sqrt(stp[6 * i + 3] * stp[6 * i + 3] + stp[6 * i + 4] * stp[6 * i + 4] + stp[6 * i + 5] * stp[6 * i + 5]) * scale;
               ch = std::max(ch, std::max(2.0 * sc / h, sw * ball[i].R));
               rd.moveC = std::max(rd.moveC, 2.0 * sc / h);
               rd.moveW = std::max(rd.moveW, sw * ball[i].R);
            }
            bool settled = true;
            for (int p = 0; p < np; p++)
               settled = settled && part[p]->Settled(stp, partOff[p], scale, h);
            moved = !(ch <= set.tolV && settled);
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
            }
            else
            {
               it.touching = false;
               it.xi = V3();
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
}
