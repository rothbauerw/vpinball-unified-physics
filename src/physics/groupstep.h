// license:GPLv3+
//
// The group step (docs/unified-physics/GroupStepDesign.txt; ProjectAssessment-2026-10-27.txt item 0 b): the balls of a touching group and their contacts moved together in one solve,
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
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "physics/rubbercontact.h"
#include "physics/pairform.h"
#include "physics/partrecords.h"

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
         if (n > used)
            used = n;
         if (v.size() < n)
            v.resize(n);
      }
      void Wipe()   // the part this use grew back to its default (WipeRecord): the size stays, so no allocation
      {
         std::fill(v.begin(), v.begin() + (std::ptrdiff_t)std::min(used, v.size()), T());
         used = 0;
      }
      operator T *() { return v.data(); }
      operator const T *() const { return v.data(); }
      template <class I> T &operator[](const I i) { assert((size_t)i < v.size()); return v.data()[i]; }
      template <class I> const T &operator[](const I i) const { assert((size_t)i < v.size()); return v.data()[i]; }

   private:
      std::vector<T> v;
      size_t used = 0;   // the most this use has grown it to
   };
   template <class T> void WipeRecord(Buf<T> &b) { b.Wipe(); }
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
   // VPX_STEP_RULE counters (diagnostic only, read and reset by contactsystem.cpp per step; nothing reads them for the physics): massless-part solves and their cost
   struct EnvStats
   {
      long solves = 0, hits = 0, misses = 0, iters = 0, capsules = 0, terms = 0;   // (terms: hinge_piece_sum's union terms probed)
      double us = 0.0;   // (no longer filled: the massless part was removed 2026-11-05; solves..us and missBy stay 0, the STEP_RULE columns keep their place)
      long missBy[7] = {};   // the solves (misses) by caller (envWho): 0 a contact's probe, 1 residual end, 2 residual middle, 3 attempt start, 4 the step rule's rate, 5 the angle after a step, 6 energy
   };
   inline EnvStats envStats;

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
   // [Contact] hinge_piece_sum: one term of a contact's union (a wire piece +1, a joint shared by m pieces -(m - 1), a part's faces +1): its own depth, its
   // normal (world) and its column's slope; the group steps each term by its own exact difference (Eval)
   struct HitTerm
   {
      double d, slope, s;
      V3 n;
   };
   // a contact's terms: kInline held in the Hit (no allocation for a usual list), the rest on the heap, grown as needed (storage, not a limit). A Hit copies only its nt terms
   struct HitTerms
   {
      static constexpr int kInline = 12;
      HitTerm &operator[](const int i) { return i < kInline ? in[i] : ex[(size_t)(i - kInline)]; }
      const HitTerm &operator[](const int i) const { return i < kInline ? in[i] : ex[(size_t)(i - kInline)]; }
      void Fit(const int n) { if (n > kInline && ex.size() < (size_t)(n - kInline)) ex.resize((size_t)(n - kInline)); }
      HitTerm in[kInline];
      std::vector<HitTerm> ex;
   };
   struct Hit
   {
      double delta = -1.0;
      V3 n;
      double reff = 0.0;
      double curvDist = 0.0;   // > 0: the ball turns round a centre this far away while it rolls along the surface (a post, a tube): the step's turning rule
      double grad = 1.0;       // the depth's gradient is n x grad (below 1 where the contact patch spans a fold of the surface)
      double kb = 0.0;         // > 0: the surface bends in series with the contact here (N/m along the normal): a massless spring, taken with the law when a touch starts and held while it lasts
      double ktb = 0.0;        // > 0: the surface bends across the normal too, in series with the shear spring (N/m, the stiffer tangential direction; a wire switch's bent wire): held as kb
      int nc = 0;              // the columns in use; 0 = a fixed surface
      Column c[7];             // one rigid body's 6 coordinates and 1 given motion of its frame: a physical size, not a table limit
      int nt = 0;              // > 0: the union's terms (hinge_piece_sum), the same list in the same order at every place of the ball; delta, n, c are the union's
      HitTerms t;
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
         ktb = o.ktb;
         nc = o.nc;
         for (int j = 0; j < nc; j++)
            c[j] = o.c[j];
         nt = o.nt;
         t.Fit(nt);
         for (int j = 0; j < nt; j++)
            t[j] = o.t[j];
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
      virtual std::string Name() const { return std::string(); }   // diagnostics (VPX_GROUP_CLAMP_TRACE): its owner's name
      // [Contact] hinge_piece_sum, a surface made of pieces: the pieces a ball at c can reach while it moves at most reach and its part turns at most swing (rad), and
      // every piece within minReach of it (the far test's own bound), in the pieces' own order. False: not a surface of pieces (its probe takes no list)
      virtual bool NearPieces(const V3 &, double, double, double, double, std::vector<int> &) const { return false; }
      // the pieces off the list are off the ball at c with its part moved by dq (trial): a distance bound, no probe
      virtual bool NearHolds(const V3 &, double, double, const std::vector<int> &) const { return true; }
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
      double omegaH = 0.39;                    // h <= omegaH / omega of the stiffest unsteady contact's normal spring (the bounce; one clock: runs/2026-11-10-one-clock/design.md.
                                               // 0.39 = the old 0.6 on the shear clock for a ball on rubber, omT = 1.53 omN; a light part's faster rubbing no longer sets the step)
      bool gonzalezMass = true;                // each stored energy's correction (Gonzalez) split by mass, the kinetic-energy metric (ContactSettings gonzalez_mass)
      bool slideFree = true;                   // a slide at the friction limit counts as steady (docs/unified-physics/GroupStepDesign.txt 16: its force is mu N, the spring's stiffness is not in the motion)
      // The step rule by motion (docs/unified-physics/GroupStepDesign.txt 9, 16, 17): a contact that is steady only needs its normal stiffness resolved to omegaHRest; an impact, a slip the spring can still catch or a fast
      // arrival keeps the strict rule above. Steady is motion, not load: its closing rate is within vSteady, the net force on its bodies cannot change that rate by more in the contact's own time (at a hit's
      // turnaround the rate is zero and the net force is the whole contact force; under a pile the net force is near zero whatever the load), and its slip is within vSteady or at the friction limit.
      // vSteady 0 = every contact strict. 0.02 m/s and 2.4 from 750 real hits against their 1 us runs (runs/2026-10-31-static-group).
      double vSteady = 0.02;                   // m/s
      double omegaHRest = 2.4;
      // stepModel 1 ([Contact] step_model): a steady touching contact is stepped by the motion in its mode, not its frequency. The midpoint step keeps the mode's energy at any h; its error is the
      // phase, (omega h)^2 / 12 of the mode's amplitude A = sqrt((closing rate / omega)^2 + (net normal push / (m omega^2))^2) per radian, so h <= sqrt(12 epsRest / A) / omega, with omega at the real
      // depth. Never shorter than rule 0 (omegaHRest / omega); arrivals, impacts and slips keep their rules. 0 (default) = rule 0 alone.
      int stepModel = 1;
      double epsRest = 1e-7;                   // m: the phase error a steady contact may carry ([Contact] eps_rest_um, default 0.1 um)
      // swingCap ([Contact] step_swing_cap, with stepModel 1): the phase error is at most the mode's whole swing 2 A (the midpoint step keeps its energy, so its
      // amplitude), so a mode whose 2 A is within epsRest sets no limit; the other rules still bound the step. false (default) = the phase law alone.
      bool swingCap = true;
      // termArrival ([Contact] term_arrival): a contact made of terms (Hit::t, hinge_piece_sum) gives each term not yet pressing (d_t < 0) the arrival rule of any contact:
      // its closing rate dDot_t = -v . n_t + slope_t * rate (the term's exact gradient), and when closing the step is at most max(0.5 * -d_t / dDot_t, hs_t) with hs_t the
      // clock an untouched approaching contact gets there (rule 3). The union's depth sees only the deepest piece, so a second piece arriving inside a step is otherwise invisible.
      // true (default) = on; false = the union's rule alone.
      bool termArrival = true;
      // dampDirSlope ([Contact] damp_dir_slope): a contact made of terms damps along G = Fel / |Fel| (the pressing terms' summed elastic force); G turns as the terms'
      // forces change, at h dG/du = (I - G g^T) K / fe (K = sum_t h tk_t tg_t tg_t^T, g = G's ball part). Newton's slope then takes fd h dG/du and dd = G . dq's change
      // through G, (kd / h) G dq^T dG/du; N's slope (kv) the latter. The force is unchanged. true (default) = on; false = G held fixed in the slope.
      bool dampDirSlope = true;
      // leverRollSlope ([Contact] lever_roll_slope; runs/2026-11-08-lead/slope_terms_design.md): Newton's slope takes two terms of laws the residual already computes.
      // The contact point rA = -n (R - dm / 2) moves with the mid depth (d dm / du = gm / 2): the slip changes by (gm / 4) sum_p wm_p x n and the torque arm by (gm / 4) s_p n, on the
      // balls' force and spin rows. The rolling limit lim = muRr R N (roll_model 1, rolling) changes with N: muRr R (h dN/du) along the torque's direction, on the spin rows.
      // The force is unchanged. true (default) = on; false = both held fixed in the slope.
      bool leverRollSlope = true;
      // twistSlope ([Contact] twist_slope; runs/2026-11-08-lead/twist_slope_design.md): Newton's slope takes the twist torque Mw's change with the depth (roll_model 1).
      // Sticking, Mw = kw1 tw0' + (kw1 h / 2 + cw) wn: dMw/dd1 = kw1' (G tw0 + h wn / 2) + cw' wn, kw1' = 3 kw1 a'(d1) / a(d1), cw' = cw kw1' / (2 kw1), G = 0 while the patch
      // grows (it holds kw0 tw0), 1 otherwise; d d1/du = gm. At the limit, Mw = s 3 pi/16 mu a(dm) N: h dMw/du = s 3 pi/16 mu (a'(dm) N h gm / 2 + a(dm) h dN/du), h dN/du the
      // normal row's own slope (kn gm; a terms contact: kv). On the balls' spin rows (h s_p n) and a carried part's rows (-(n . W)). The torque is unchanged.
      // true (default) = on; false = only dMw/dwn (today's slope).
      bool twistSlope = true;
      // endGradSlope ([Contact] end_grad_slope; runs/2026-11-08-lead/end_gradient_design.md): a terms contact's (hinge_piece_sum) pressing term t gives h tk_t tg_t (d a1_t/du)^T,
      // the force along its own direction tg_t (the mid gradient: the row side) times its END depth's change over the unknowns. d a1_t/du is the end place's own gradient
      // (h1.t[t]: -n, slope), not the mid one: a ball crossing a wire's crest turns the normal within the step (the mid y slope 6 % off at 3497 ms on the switch table).
      // Taken in the terms' slope, G's turning (dampDirSlope) and N's slope kv (friction, lever-roll and twist limits). J is then not symmetric (it is the true derivative;
      // SolveLinear pivots on any matrix). The force is unchanged. true (default) = on; false = the mid gradient on both sides (today's slope).
      bool endGradSlope = true;
      // termForceRows ([Contact] term_force_rows; runs/2026-11-08-lead/term_rows_design.md): a terms contact's pressing term t pushes along gdT_t, its mid gradient tg_t
      // Gonzalez-corrected so that gdT_t . dq = a1_t - a0_t (the force sum_t s fT_t gdT_t, G = Fel / |Fel| and dd = G . dq are built on it). Newton's slope takes that same
      // direction as the row side of each term's slope: h tk_t gdT_t tgE_t^T in the terms' slope, in G's turning K (so gK and gW: dd's change through G is
      // (1 / fe) sum_t h tk_t [gdT_t . dq - dd (g . gdT_t)] tgE_t, 0 for one term) and in N's slope kv, kvE. The exact derivative of the residual with each term's
      // direction held fixed. The force is unchanged. true (default) = on; false = the mid gradient tg_t as the row (today's slope).
      bool termForceRows = true;
      // gripInStep ([Contact] grip_in_step; runs/2026-11-08-lead/cord_grip_design.md): a part's grip limit (the cord's mu Tm Gfree) is worked out inside each Eval from that
      // Eval's own mid tension and mid geometry, as a contact's mu N is, so the solve's residual is one fixed function of the unknowns; Newton's slope takes the sliding grip's
      // change through the limit (mu Gfree dTm/du). true (default) = on; false = the limit refreshed between the solve's rounds from the last round (Part::UpdateGrip), which
      // turns off the line search's backtrack, the cached-residual exit and the end checks in any round it moved.
      bool gripInStep = true;
      // dampingIntegral ([Contact] damping_integral; runs/2026-11-08-lead/damping_integral_design.md): every depth-dependent damping in its integrated form, the law's own work
      // over the step at a steady rate: c = (W(d1) - W(d0)) / (d1 - d0), W(d) = integral of c(s) from 0 (0 for no touch), short moves by quadrature; slope d(c (d1 - d0))/d d1 = c(d1).
      // Kuwabara-Kono (W = gamma 2/3 sqrt(reff) Uc^1.5), the loss tangent's dashpot (W = 0.8 tanDelta sqrt(1.5 kH m) Uc^1.25) and Hunt-Crossley (with a bend: 1.5 alpha 0.4 kH Uc^2.5,
      // the contact's own stored energy; without one, today's chord). The twist limit 3 pi/16 mu a N takes the patch a(fe) of the step's mean elastic force fe (the law's own
      // inverse: F(dE) = fe), not the mid depth's: no cusp where a contact arrives or leaves within a step; the twist is on while kw1 > 0 or a(fe) > 0, with no gate on dm.
      // Newton's slope takes their exact derivatives (through fe). true (default) = on; false = the committed forms (point values at dm, the KK chord/midpoint switch, a(dm)).
      bool dampingIntegral = true;
      // rollingPatch ([Contact] rolling_patch; runs/2026-11-08-lead/rolling_patch_design.md): rolling resistance's patch is the patch of the step's mean elastic force fe,
      // aE = A(dE) with F(dE) = fe (PatchOfForce, as the twist limit's under damping_integral), not the mid depth's a(dm): its pressure pa = N / (pi aE^2) and its gate
      // (rolling acts while N > 0 and aE > 0, no switch on dm). Newton's slope (lever_roll_slope) takes the limit lim = muRr R N's full change: R (muRr + N dmuRr/dN) h dN/du
      // + R N dmuRr/daE daE/dfe h dfe/du (muRr = 2.34 (pa / E) tanDelta between its floor rolling_resistance and its cap kRollingMuMax; at either, muRr's slope is 0).
      // roll_model 0 takes the same patch (no slope there, as today). true (default) = on; false = the mid depth's a(dm), the limit's slope muRr R h dN/du with muRr held.
      bool rollingPatch = true;
      // stickZoneTrial ([Contact] stick_zone_trial; runs/2026-11-08-lead/slip_reversal_design.md section 2b): Newton's move and slopes unchanged; before the line search's
      // full move (then 1/2, 1/4 ...), it tries, for each contact whose shear patch slides at the current point with N > 0, the fraction s of the move at the middle of the
      // part of (0, 1] where that contact's trial force lies inside its Coulomb limit: |a + s C dut| <= mu (N0 + s dN), a quadratic in s, with a = the trial force at the
      // current point, C = kS h / 2 + ct, kS, ct, the stretch and the contact's geometry (normal, lever) held at the current point, dut the slip rate's change over the
      // move (ut's own linear form in the move: centre / h, half the spin change x lever, carried columns / h), dN the normal row's own slope over the move. A sliding
      // patch's slope has no slope along the slip (Ft = mu N f^), so the full move overshoots a stick zone ~15x narrower (the Drain kicker 2-cycle). Accepted by the exact
      // residual as every move; candidates ascending, then 1, 1/2 ... A contact whose interval contains s = 1 (the full move already ends inside its zone) gives no
      // candidate; only one whose interval ends before s = 1 does. true (default) = on; false = the line search tries 1, 1/2 ... only.
      bool stickZoneTrial = true;
      // patchPath ([Contact] patch_path; runs/2026-11-08-lead/patch_state_design.md): every patch-sized coefficient of the step, the shear spring Kt = ktPerA a, its damper
      // ct = shearTan sqrt(Kt mt), the twist and rolling springs Ktw, Kr ~ a^3 and their dampers ~ sqrt(Ktw), sqrt(Kr), is its mean over the depth the step travels,
      // k_bar = (P(d1) - P(d0)) / (d1 - d0), P(d) = integral of k from 0 (0 for no touch; closed forms in u = Uc(d), dd = (1 + 1.5 kH sqrt(u) / kb) du; short moves by quadrature:
      // Item::MeanOver). The stretch carries over in two patch changes, each by the patch rule (a patch that grows holds the force, one that shrinks holds the stretch; the
      // spring's fall a loss): from Kt(d0) to k_bar before the step's move, from k_bar to Kt(d1) after it; the end stretch is stored on Kt(d1), the next step's start (no jump
      // across steps). Newton's slope takes the means' exact slopes d k_bar / d d1 = (k(d1) - k_bar) / (d1 - d0) (short moves: the quadrature's). true (default) = on;
      // false = the end depth's coefficients Kt(d1), Ktw(d1), Kr(d1) and their dampers for the whole step, with their slopes at d1.
      bool patchPath = true;
      // partGuessAccel ([Contact] part_guess_accel): a part's first guess of its coordinate's move carries its last step's acceleration, (rate + a h / 2) h, a the force of the
      // move taken over its mass, as a ball's guess does. Newton's start only: the law and the solved step are unchanged. false (default) = rate h.
      bool partGuessAccel = false;
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
      // the difference Jacobian. 0 = off (normal). Its default from VPX_GROUP_JAC_MODE (diagnostics, read once; unset = 0); with VPX_GROUP_ROUND_TRACE=1 each round prints the parts' rows and columns both ways
      int jacMode = []() { static const int m = []() { const char *const e = getenv("VPX_GROUP_JAC_MODE"); return e != nullptr ? std::max(0, std::min(2, atoi(e))) : 0; }(); return m; }();
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

   // The engine's step settings ([Contact] keys and the step floor, ContactSettings::ToStep, set when the ini is read): one copy, taken whole by every group the engine steps,
   // an island's (with its own pair and gravity on top) and a part's alone (AdvanceAlone). Its defaults are Settings' own (an offline program that never sets it steps on those)
   inline Settings &EngineSettings()
   {
      static Settings s;
      return s;
   }

   // One rest rule for every moving part (ModelStandards 3c: a part with nothing driving it and nothing touching it sleeps), alone and in a group alike:
   // nothing drives it (each part names its drives: a coil's supply, a pulse, a moving stop, a pending trigger, a motor), no ball touches it, over one whole cycle no point of its surface
   // moved more than kMoveM, nor could its surface speed move it more in the next cycle, and every internal state has died away (a coil's current, each winding: |i| <= kStateFrac of its
   // full scale, a level, not a change: a sleeping part is not stepped, so a current still decaying would freeze with it). Alone, the time watched adds up over the part's steps (a cycle can
   // come in shares); a move past the limit, a live state, a drive or a touch starts the watch again from where the part is. In a group the window is the group's cycle: the group restarts
   // the watch of every part a ball touched where it builds the group (Restart, the part's pose then is the one mark), with the same limits.
   // (GroupStepDesign beside the rest rule; runs/2026-11-09-sleep-fix/design.md)
   struct RestWatch
   {
      static constexpr double kMoveM = 1e-9;      // m, at the part's reach: the one limit
      static constexpr double kStateFrac = 1e-9;  // an internal state at rest: at most this fraction of its full scale (a coil's current with the coil held on). With the supply off a current only
                                                  // decays, so its level at a share's end is enough: at a 100 ms time constant a change test would pass it still at 1e-7 of full
      static constexpr double kCycleS = 1e-3;     // s: one physics cycle
      static constexpr int kMax = 4, kMaxState = 2;
      struct Pose   // the part's coordinates and the lever that turns each into a surface move (m per unit); its internal states and the full scale of each
      {
         int n = 0, ns = 0;
         double q[kMax] = {}, reach[kMax] = {};
         double s[kMaxState] = {}, full[kMaxState] = {};
         Pose &Add(const double v, const double r)
         {
            q[n] = v;
            reach[n] = r;
            n++;
            return *this;
         }
         Pose &State(const double v, const double fullScale)
         {
            s[ns] = v;
            full[ns] = fullScale;
            ns++;
            return *this;
         }
         bool Live() const   // an internal state not yet died away
         {
            for (int i = 0; i < ns; i++)
               if (std::fabs(s[i]) > kStateFrac * full[i])
                  return true;
            return false;
         }
      };
      Pose mark;
      double t = 0.0;
      bool on = false;
      bool cycled = false;   // the group restarted it (Restart): mark is the pose at the start of the part's last cycle in a group, the window AtRest judges
      Pose cycleMark;
      void Reset()
      {
         on = false;
         t = 0.0;
      }
      static double Moved(const Pose &p, const Pose &from)
      {
         double m = 0.0;
         for (int i = 0; i < p.n; i++)
            m += std::fabs(p.q[i] - from.q[i]) * p.reach[i];
         return m;
      }
      double Moved(const Pose &p) const { return Moved(p, mark); }
      static bool Still(const double movedM, const double speedMS) { return movedM <= kMoveM && speedMS * kCycleS <= kMoveM; }
      // a group holding the part is built this cycle (a ball touches it): the watch alone starts again, and the pose now is the mark of the group's cycle
      void Restart(const Pose &p)
      {
         Reset();
         cycleMark = p;
         cycled = true;
      }
      // the group's test (AtRest), with nothing driving the part: over its last cycle in a group it moved no more than the limit, its speed now could not, and its states have died away
      bool Rests(const Pose &p, const double speedMS) const { return cycled && !p.Live() && Still(Moved(p, cycleMark), speedMS); }
      // the move over its last cycle in a group (pile sleep's log); -1 = never in a group
      double CycleMoved(const Pose &p) const { return cycled ? Moved(p, cycleMark) : -1.0; }
      // a share dt of the part alone ended at pose p, surface speed speedMS (m/s): true = it rests. A live state restarts the watch as a drive does
      bool Step(const Pose &p, const double speedMS, const bool driven, const bool touched, const double dt)
      {
         if (driven || touched || p.Live())
         {
            Reset();
            return false;
         }
         if (!on || Moved(p) > kMoveM)
         {
            mark = p;
            t = 0.0;
            on = true;
            return false;
         }
         t += dt;
         return t >= kCycleS * (1.0 - 1e-9) && Still(0.0, speedMS);
      }
      // diagnostics (VPX_GROUP_STATS): per part kind (0 jet bumper, 1 gate / spinner / wire switch, 2 disc, 3 flipper, 4 sling plate) the steps alone taken, and how often one ended asleep after being awake
      static long *Counts()
      {
         static long c[10] = {};
         return c;
      }
      static void Book(const int kind, const bool was, const bool now)
      {
         Counts()[kind]++;
         if (!was && now)
            Counts()[5 + kind]++;
      }
   };


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

   // The step's inertia, one law for every body (balls and parts): over the move dq the momentum goes from p0 to p1 with dq = h Mb^-1 (p0 + p1) / 2,
   // Mb^-1 the mean of the start and end inverse inertia. The body's row is p1 - p0 - h Q (Q its forces over the move), its end rate p1 / M1.
   // A fixed mass (a ball, a gate's plate, the jet bumper's skirt and ring, the flipper's plunger) is the same law with Mb = m: the row 2 m dq / h - 2 m v, the end rate 2 dq / h - v
   template <class T> inline T InertiaRow(const T &dq, const T &v, const double m, const double h) { return dq * (2.0 * m / h) - v * (2.0 * m); }
   inline double InertiaSlope(const double m, const double h) { return 2.0 * m / h; }   // the row's slope along dq
   template <class T> inline T EndRate(const T &dq, const T &v, const double h) { return dq * (2.0 / h) - v; }
   // An inertia that varies with the pose (the flipper's bat, the sling plate): I0 at the step's start, I1 at its end, p0 = I0 w0
   struct MeanInertia
   {
      double I1, invMean, p0, p1;
      MeanInertia(const double I0, const double I1_, const double p0_, const double dq, const double h) : I1(I1_), invMean(0.5 * (1.0 / I0 + 1.0 / I1_)), p0(p0_), p1(2.0 * dq / (h * invMean) - p0_) {}
      double Row() const { return p1 - p0; }
      double Slope(const double h) const { return 2.0 / (h * invMean); }
      double EndRate() const { return p1 / I1; }
      // dT/dq at fixed p, T = p^2 / (2 I) (the part's force takes it off): the step's mean of p^2 / 2 times the part's chord (or slope, for a very small move) of 1/I over the move
      double KineticSlope(const double dInvI) const { return dInvI * 0.25 * (p0 * p0 + p1 * p1); }
   };

   // A moving part of a group (GroupStepDesign.txt piece 1; 10.1; 18): it describes its bodies, stored energy, losses, drives and grip, and the group moves them with
   // the balls in one solve. The unknowns x: per ball the move (3) and the spin change (3), then each part's own (NumDof of them, from off). The balls' residual
   // rows are 2 m dc / h - 2 m v - h F and I dw - h M; a part adds its forces on the balls (averaged over the step) to F and M and writes its own rows.
   // Every call of a step reads the balls' state at the step's start (BeginStep's b) and the move x.
   // Its values by lifetime (PartRecords): a part declares its state records, step work and attempt work; the group alone saves, restores and wipes them.
   class Part : public PartRecords
   {
   public:
      virtual ~Part() = default;
      bool guessAccel = false;   // [Contact] part_guess_accel (Settings::partGuessAccel, set by the group before BeginStep)
      bool gripInStep = true;    // [Contact] grip_in_step (Settings::gripInStep, set by the group before BeginStep): its grip limit from each Eval's own unknowns; false = UpdateGrip between rounds
      mutable const std::vector<int> *nearNow = nullptr;   // [Contact] hinge_piece_sum: the probe of a surface it carries takes only these pieces (this attempt's list); null = all
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
      // > 0: the move x carries a body past a hard stop by this much (the step is shortened to end on it; piece 1 item 7). mask: the stops counted (bit k = the part's stop k)
      virtual double StopExcess(const double *, int, unsigned = ~0u) const { return -INFINITY; }
      // the stops the move x carries it past, each on its own (bit k): only these shorten a step, so a stop a body is on and leaving never does
      virtual unsigned CrossingStops(const double *x, int off) const { return StopExcess(x, off) > 0.0 ? ~0u : 0u; }
      virtual bool ArrivingAtStop() const { return true; }               // the move last tested by StopExcess is closing on what it came near
      virtual bool PinOnStop(const double *, int) { return false; }      // on its hard stop already and the move x just solved would take it on past: it stays there. False: nothing changed
      virtual bool FreedByStep() { return false; }                       // a body with no unknown: the forces of the step just solved take it off, so it has one and the step is solved again
      // the balls' reaction on it in this evaluation, before its Eval: per coordinate (one with no unknown too), and the force on the whole part (a pivot's friction goes with its load). Only a part carrying a surface hears it
      virtual void Loads(const double *, const V3 &) {}
      // the step taken (after an Eval at the solution, before the group moves its balls); onStop: it ended on the hard stop. Its drive work and losses of the
      // step, and whether a ball went through it
      virtual void Accept(const double *x, int off, double h, bool onStop, double &driveJ, double &lossJ, bool &through) = 0;
      virtual void AfterStep(const Group &) {}   // the step taken and the balls moved: what the part keeps of its own shapes' contacts
      // an energy step (Group::EnergyStep): the state kept and put back (its state records, the same the engine's snapshot copies); its bodies' kinetic energy and all their speeds times s; what a step's forces left for the next one is put back
      void SaveState()
      {
         sv.resize(StateBytes());
         SaveRecords(sv.data());
      }
      void RestoreState()
      {
         RestoreRecords(sv.data());
         AfterRestore();
      }
      virtual void AfterRestore() {}   // its state put back: what is derived from it is made again (the bat's pivot)
      const unsigned char *SavedState() const { return sv.data(); }   // the net's saved records, in the order the part declared them (its notes)
      virtual double KineticJ() const { return 0.0; }
      virtual void ScaleSpeeds(double) {}
      virtual void AfterEnergyStep() {}
      virtual void HoldArcs() {}   // no solve in this attempt: what the part solves within itself stays as it is
      // diagnostics (rule 8): an attempt of a step failed and was halved (slot: its place in Group::halveDiag); the part keeps what its last round evaluated, Res
      // the residual of that round (its rows from off). Told only when an energy step is reported
      virtual void KeepFailNote(int, const double *, int) {}
      virtual std::string FailNote(int) const { return std::string(); }
      virtual std::string StartNote() const { return std::string(); }
      virtual std::string StopNote() const { return std::string(); }    // diagnostics (VPX_GROUP_STOP_TRACE): its state at a stop check   // after an energy step: the state it saved at the step's start, full precision (for an exact replay)
      virtual std::string SearchNote(const double *, int, double) const { return std::string(); }   // diagnostics (VPX_GROUP_SEARCH_TRACE): the part's state for the move x over h
      virtual std::string RowNote() const { return std::string(); }   // diagnostics (VPX_GROUP_ROUND_TRACE): the terms of its rows at its last evaluation
      virtual std::string LastProbeNote() const { return std::string(); }   // diagnostics (VPX_GROUP_ROUND_TRACE): what its last surface probe found (a wire's pieces)
      virtual std::string ContactNote(const double *, int, double, const double *) const { return std::string(); }
      virtual std::string HistoryNote() const { return std::string(); }   // diagnostics (VPX_GROUP_HISTORY_TRACE): its state and energies after an accepted step (read only)   // diagnostics (VPX_GROUP_CONTACT_TRACE): its state and row terms for the move x over h, its residual Res (read only)
      virtual void BallsMoved(const Ball *, int) {}              // the balls were moved with no step of its own
      virtual void ExplicitMove(const double *, int, double, double *) const {}   // its own unknowns in an energy step's explicit round
      virtual double EndResidual() const { return 0.0; }         // > 0: a step of a group holding this part counts as solved only with its residual within this (m/s; Settings::tolEndRes)

   private:
      std::vector<unsigned char> sv;   // the net's copy of its state (SaveState)
   };

   // One contact: ball a on shape s, or ball a on ball b. Its law is taken when a touch starts and held while it lasts.
   struct Item
   {
      int a = -1, b = -1, s = -1;
      int nearK = -1;                 // >= 0: inside an attempt its probes take the group's list nearL[nearK] (a surface of pieces, hinge_piece_sum)
      bool touching = false;
      V3 xi;                          // the shear spring's stretch (in the contact plane)
      RubberContact::RollAngle ra;    // the rolling and twist springs' angles (roll_model 1)
      double dZero = 0.0;             // the depth its squash counts from: a ball found inside has stored nothing (no energy without the work that stored it; only ever down)
      // the held law
      double kH = 0.0, cL = 0.0, hL = 0.0, reff = 0.0, ktPerA = 0.0, E = 0.0, mRed = 0.0, mt = 0.0;
      double alpha = 0.0, gammaKK = 0.0, tanDelta = 0.0, shearTan = 0.0, mu = 0.0, rollingMu = 0.0;
      bool normalLoss = false, twist = false;
      double kb = 0.0;                // the surface's bend in series with the contact (0 = none): V(d) = min over u of Hertz(u) + 1/2 kb (d - u)^2, kH u^1.5 = kb (d - u)
      double ktb = 0.0;               // the surface's bend across the normal in series with the shear spring (0 = none): Kt = Kc kt_b / (Kc + kt_b), Kc the patch's ktPerA a
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
      double ASlope(double d) const   // da/dd = 1/2 sqrt(reff / Uc) Uc'(d), Uc' = kb / (kc + kb) with the bend in series (kc = 1.5 kH sqrt(Uc)), else 1
      {
         if (d <= 0.0)
            return 0.0;
         const double u = Uc(d);
         if (u <= 0.0)
            return 0.0;
         return kb > 0.0 ? 0.5 * std::sqrt(reff / u) * kb / (1.5 * kH * std::sqrt(u) + kb) : 0.5 * std::sqrt(reff / u);
      }
      double Kt(double d) const
      {
         const double k = ktPerA * A(d);
         return ktb > 0.0 && k > 0.0 ? k * ktb / (k + ktb) : k;
      }
      // with the bend across (ktb > 0): the integral of Kt over the depth from 0, for its mean over a step (patch_path; Kt is no power of u then). With x = Kc / kt_b = ktU s / kt_b,
      // s = sqrt(u), ktU = ktPerA sqrt(reff), dd = (1 + b s) 2 s ds, b = 1.5 kH / kb (0 with no bend): 2 kt_b^3 / ktU^2 (F2(x) + beta F3(x)), beta = b kt_b / ktU,
      // F2 = x^2/2 - x + ln(1 + x), F3 = x^3/3 - x^2/2 + x - ln(1 + x); below x 0.05 by their series (the closed forms cancel there)
      double KtIntegral(const double d) const
      {
         if (d <= 0.0 || ktPerA <= 0.0)
            return 0.0;
         const double ktU = ktPerA * std::sqrt(reff), x = ktU * std::sqrt(Uc(d)) / ktb, beta = (kb > 0.0 ? 1.5 * kH / kb : 0.0) * ktb / ktU;
         double f2, f3;
         if (x < 0.05)
         {
            f2 = f3 = 0.0;
            double xp = x * x * x, sg = 1.0;   // F2 = sum (-1)^k x^(k+3) / (k+3), F3 = sum (-1)^k x^(k+4) / (k+4)
            for (int k = 0; k < 12; k++, sg = -sg)
            {
               f2 += sg * xp / (k + 3);
               f3 += sg * xp * x / (k + 4);
               xp *= x;
            }
         }
         else
         {
            const double l = std::log1p(x);
            f2 = x * x / 2.0 - x + l;
            f3 = x * x * x / 3.0 - x * x / 2.0 + x - l;
         }
         return 2.0 * ktb * ktb * ktb / (ktU * ktU) * (f2 + beta * f3);
      }
      // the rolling and twist springs of the patch (roll_model 1): tilting (4/3) E* a^3 and torsion (16/3) G* a^3 of a Hertz contact, G* = E* (1 - nu) / 2 with nu 0.3
      // (E* = 3 kH / (4 sqrt(reff)): the normal law's own contact modulus, a pair's compliances added)
      double Kr(double d) const { const double a = A(d); return kH / std::sqrt(reff) * a * a * a; }
      double Ktw(double d) const { const double a = A(d); return 16.0 / 3.0 * 0.35 * 0.75 * kH / std::sqrt(reff) * a * a * a; }
      double KtSlope(double d) const
      {
         if (d <= 0.0)
            return 0.0;
         double sl;
         if (kb > 0.0)
         {
            const double u = Uc(d), kc = 1.5 * kH * std::sqrt(u);
            sl = u > 0.0 ? 0.5 * ktPerA * std::sqrt(reff / u) * kb / (kc + kb) : 0.0;
         }
         else
            sl = 0.5 * ktPerA * std::sqrt(reff / d);
         if (ktb > 0.0)   // the bend across in series: dKt/dKc = (kt_b / (Kc + kt_b))^2
         {
            const double r = ktb / (ktPerA * A(d) + ktb);
            sl *= r * r;
         }
         return sl;
      }
      double KKIntegral(const double d) const   // 2/3 sqrt(reff) d^1.5, 0 for no touch
      {
         if (d <= 0.0)
            return 0.0;
         const double u = Uc(d);
         return 2.0 / 3.0 * std::sqrt(reff) * u * std::sqrt(u);
      }
      // Viscoelastic damping over the step: the law's own integral at a steady rate, as Hunt-Crossley already is (runs/2026-10-31-gate-spinner-disc/patch_gs_rotor.py, without its bend in series). F = gamma a(d) d_dot, so the mean
      // coefficient over a move is gamma (G(d1) - G(d0)) / (d1 - d0), on the end depths as the spring's chord (never the mid-step positions' depth, which a curved surface puts deeper:
      // runs/2026-11-08-lead/kk_damping_check.md). jacTerm: its exact slope for the Jacobian (d(c x (d1 - d0)) / d d1 = gamma a(d1) Uc'(d1))
      double KKRate(const double d) const { return A(d) * (kb > 0.0 ? kb / (kb + 1.5 * kH * std::sqrt(Uc(d))) : 1.0); }   // a(d) Uc'(d) = dG/dd: the coefficient per gamma at depth d, the bend in series taking its share
      double KKStep(const double d0, const double d1, double &jacTerm) const
      {
         jacTerm = gammaKK * (d1 > 0.0 ? KKRate(d1) : 0.0);
         return gammaKK * MeanOver(d0, d1, [this](const double d) { return KKRate(d); }, [this](const double d) { return KKIntegral(d); });
      }
      // [Contact] damping_integral 0: the committed form (the chord for a long move, the coefficient at the mid depth dm for a short one)
      double KKStepMid(const double d0, const double d1, const double dm, double &jacTerm) const
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
      // [Contact] damping_integral, patch_path: the mean of a depth-dependent coefficient c(d) over a move d0 -> d1 at a steady rate, (W(d1) - W(d0)) / (d1 - d0), W its integral
      // from 0 (rate: c(d), integral: W(d), both 0 for no touch); a short move (both touching, the change under 1 % of the depth) by 3-point Gauss quadrature (no cancellation).
      // W rises with d, so the mean is >= 0. The rate may return several coefficients at once (Pows). MeanSlope: the mean's exact change with d1, (c(d1) - mean) / (d1 - d0),
      // a short move's from the quadrature's own nodes (rateSlope: dc/dd)
      static bool ShortMove(const double d0, const double d1) { return d0 > 0.0 && d1 > 0.0 && std::fabs(d1 - d0) < 0.01 * std::max(d0, d1); }
      template <class Rate, class Integral> static auto MeanOver(const double d0, const double d1, const Rate &rate, const Integral &integral) -> decltype(rate(d0))
      {
         if (d0 <= 0.0 && d1 <= 0.0)
            return decltype(rate(d0))();
         if (ShortMove(d0, d1))
         {
            const double dm = 0.5 * (d0 + d1), g = 0.3872983346207417 * (d1 - d0);
            return (5.0 * rate(dm - g) + 8.0 * rate(dm) + 5.0 * rate(dm + g)) / 18.0;
         }
         return (integral(d1) - integral(d0)) / (d1 - d0);
      }
      template <class Rate, class RateSlope, class T> static T MeanSlope(const double d0, const double d1, const Rate &rate, const RateSlope &rateSlope, const T &mean)
      {
         if (d0 <= 0.0 && d1 <= 0.0)
            return T();
         if (ShortMove(d0, d1))   // d/dd1 of the nodes dm -+ g: 1/2 -+ 0.3873
         {
            const double dm = 0.5 * (d0 + d1), g = 0.3872983346207417 * (d1 - d0);
            return (5.0 * (0.5 - 0.3872983346207417) * rateSlope(dm - g) + 4.0 * rateSlope(dm) + 5.0 * (0.5 + 0.3872983346207417) * rateSlope(dm + g)) / 18.0;
         }
         return (rate(d1) - mean) / (d1 - d0);
      }
      // [Contact] patch_path: the patch-sized coefficients go as powers of the contact's compression u = Uc(d) (a = sqrt(reff u)): Kt = ktPerA sqrt(reff) u^0.5, its damper
      // shearTan sqrt(Kt mt) ~ u^0.25, Ktw and Kr ~ a^3 ~ u^1.5, their dampers ~ u^0.75. Pows holds the four (p05, p025, p15, p075)
      struct Pows
      {
         double p05 = 0.0, p025 = 0.0, p15 = 0.0, p075 = 0.0;
         friend Pows operator*(const double s, const Pows &a) { Pows r; r.p05 = s * a.p05; r.p025 = s * a.p025; r.p15 = s * a.p15; r.p075 = s * a.p075; return r; }
         friend Pows operator+(const Pows &a, const Pows &b) { Pows r; r.p05 = a.p05 + b.p05; r.p025 = a.p025 + b.p025; r.p15 = a.p15 + b.p15; r.p075 = a.p075 + b.p075; return r; }
         friend Pows operator-(const Pows &a, const Pows &b) { Pows r; r.p05 = a.p05 - b.p05; r.p025 = a.p025 - b.p025; r.p15 = a.p15 - b.p15; r.p075 = a.p075 - b.p075; return r; }
         friend Pows operator/(const Pows &a, const double s) { Pows r; r.p05 = a.p05 / s; r.p025 = a.p025 / s; r.p15 = a.p15 / s; r.p075 = a.p075 / s; return r; }
      };
      Pows PowRate(const double d) const   // u^p at depth d (0 for no touch)
      {
         Pows r;
         if (d <= 0.0)
            return r;
         const double u = Uc(d), s = std::sqrt(u), q = std::sqrt(s);
         r.p05 = s;
         r.p025 = q;
         r.p15 = u * s;
         r.p075 = s * q;
         return r;
      }
      Pows PowIntegral(const double d) const   // integral of u^p over the depth from 0: dd = (1 + b sqrt(u)) du, b = 1.5 kH / kb (0 with no bend), so u^(p+1)/(p+1) + b u^(p+1.5)/(p+1.5)
      {
         Pows r;
         if (d <= 0.0)
            return r;
         const double u = Uc(d), s = std::sqrt(u), q = std::sqrt(s), b = kb > 0.0 ? 1.5 * kH / kb : 0.0;
         r.p05 = u * s / 1.5 + b * u * u / 2.0;
         r.p025 = u * q / 1.25 + b * u * s * q / 1.75;
         r.p15 = u * u * s / 2.5 + b * u * u * u / 3.0;
         r.p075 = u * s * q / 1.75 + b * u * u * q / 2.25;
         return r;
      }
      Pows PowSlope(const double d) const   // d(u^p)/dd = p u^(p-1) Uc'(d) (a short move only: d > 0)
      {
         Pows r;
         if (d <= 0.0)
            return r;
         const double u = Uc(d), s = std::sqrt(u), q = std::sqrt(s), us = kb > 0.0 ? kb / (kb + 1.5 * kH * s) : 1.0;
         if (u <= 0.0)
            return r;
         r.p05 = 0.5 / s * us;
         r.p025 = 0.25 / (s * q) * us;
         r.p15 = 1.5 * s * us;
         r.p075 = 0.75 / q * us;
         return r;
      }
      // the means of u^p over the step's depth path d0 -> d1 and their slopes in d1 (patch_path)
      void PatchMeans(const double d0, const double d1, Pows &mean, Pows &slope) const
      {
         mean = MeanOver(d0, d1, [this](const double d) { return PowRate(d); }, [this](const double d) { return PowIntegral(d); });
         slope = MeanSlope(d0, d1, [this](const double d) { return PowRate(d); }, [this](const double d) { return PowSlope(d); }, mean);
      }
      double UcSlope(const double d) const { return d > 0.0 && kb > 0.0 ? kb / (kb + 1.5 * kH * std::sqrt(Uc(d))) : 1.0; }   // dUc/dd: the contact's share of a depth change
      // the loss tangent's dashpot on the contact's own compression u = Uc(d): c(u) = tanDelta sqrt(1.5 kH m) u^(1/4) (= tanDelta sqrt(1.5 F m / u), Hertz), per unit of d c(Uc) Uc';
      // W = 0.8 tanDelta sqrt(1.5 kH m) Uc^(5/4). LossStep: the mean over a move, jacTerm its slope c(Uc(d1)) Uc'(d1)
      double LossRate(const double d) const { return d > 0.0 ? tanDelta * std::sqrt(1.5 * kH * mRed) * std::sqrt(std::sqrt(Uc(d))) * UcSlope(d) : 0.0; }
      double LossIntegral(const double d) const { if (d <= 0.0) return 0.0; const double u = Uc(d); return 0.8 * tanDelta * std::sqrt(1.5 * kH * mRed) * u * std::sqrt(std::sqrt(u)); }
      double LossStep(const double d0, const double d1, double &jacTerm) const
      {
         jacTerm = LossRate(d1);
         return MeanOver(d0, d1, [this](const double d) { return LossRate(d); }, [this](const double d) { return LossIntegral(d); });
      }
      // Hunt-Crossley 1.5 alpha F on the contact's own compression through a bend in series: W = 1.5 alpha 0.4 kH Uc^2.5 (the contact's stored energy, not the bend's), per unit of d
      // 1.5 alpha F_c(Uc) Uc'. (Without a bend the caller keeps today's chord 1.5 alpha fe, the same law.) HCBendStep: the mean over a move, jacTerm its slope
      double HCBendRate(const double d) const { if (d <= 0.0) return 0.0; const double u = Uc(d); return 1.5 * alpha * kH * u * std::sqrt(u) * UcSlope(d); }
      double HCBendIntegral(const double d) const { if (d <= 0.0) return 0.0; const double u = Uc(d); return 1.5 * alpha * 0.4 * kH * u * u * std::sqrt(u); }
      double HCBendStep(const double d0, const double d1, double &jacTerm) const
      {
         jacTerm = HCBendRate(d1);
         return MeanOver(d0, d1, [this](const double d) { return HCBendRate(d); }, [this](const double d) { return HCBendIntegral(d); });
      }
      // [Contact] damping_integral: the patch of the depth dE whose elastic force is f, F(dE) = f (the law's own inverse), a = A(dE), and its slope da/df = A'(dE) / K(dE).
      // Hertz with or without a bend: Uc = (f / kH)^(2/3), a = sqrt(reff Uc), da/df = a / (3 f). With the layer (cL > 0): dE by Newton on the convex F from the Hertz depth above it
      void PatchOfForce(const double f, double &a, double &dadf) const
      {
         a = dadf = 0.0;
         if (f <= 0.0)
            return;
         if (cL <= 0.0 || kb > 0.0)   // (the bend's law has no layer: F = kH Uc^1.5)
         {
            a = std::sqrt(reff) * std::cbrt(f / kH);
            dadf = a / (3.0 * f);
            return;
         }
         double d = std::min(std::pow(f / kH, 2.0 / 3.0), hL * (1.0 - 1e-9));
         for (int i = 0; i < 60; i++)
         {
            const double k = K(d);
            if (k <= 0.0)
               break;
            const double dd = (F(d) - f) / k;
            d = std::max(d - dd, 0.5 * d);
            if (std::fabs(dd) <= 1e-14 * d)
               break;
         }
         a = A(d);
         const double k = K(d);
         dadf = k > 0.0 ? ASlope(d) / k : 0.0;
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

   // The group step's own values of one use (PartRecords: use work, with its buffers): set while it steps an island, never read by the next use
   struct GroupUseWork
   {
      int lastRounds = 0;
      int nHalveDiag = 0;
      char failWhy = '?';
      bool trialCached = false;
      bool warmOn = false;
      double driveW = 0.0;
      int nCoord = 0;
      bool carried = false;      // a shape of this group is carried by a part
      int nUnk = 0;              // unknowns of the step being solved
      bool partEvalFailed = false;
      bool nearOn = false;
      double nearGrow = 1.0;     // the reach's factor: x4 each time a solve's end left its list (the lists then grow to all the pieces)
      bool farOff = false;
   };

   class Group : public PartRecords, public GroupUseWork
   {
   public:
      // Its work by lifetime (PartRecords): every buffer and value a use writes is wiped at the start of the next (BeginUse). Kept: the settings; the stamps (a counter:
      // Configure moves it on, so no kept probe is taken as now); the books and diagnostics (the owner reads them after the use; reset where they start); the step history
      // (a log); the Jacobian's store (zeroed before each forming, Eval); the balls, shapes, parts and pair forms (the use's input, added after BeginUse)
      Group()
      {
         WipeEachUse<GroupUseWork>(this);
         WipeEachUse(&item);
         for (Buf<Item> *b : { &law, &esI0 })
            WipeEachUse(b);
         for (Buf<Hit> *b : { &hNow, &h1c, &hmc, &h0s })
            WipeEachUse(b);
         for (Buf<double> *b : { &d0s, &x, &resS, &stpS, &xPrevS, &stpPrevS, &trialS, &zeroS, &esRes, &xpS, &rpS, &rmS, &qEnd, &qMid, &qF, &fieldE0 })
            WipeEachUse(b);
         for (Buf<Ball> *b : { &b1r, &bmr, &evB1, &evBm, &farB1, &farBm, &esB0 })
            WipeEachUse(b);
         for (Buf<V3> *b : { &aLast, &alphaLast, &fNetS, &evF, &evM, &evDc, &evDw, &evWm, &esA0, &esAl0, &warmA, &warmAl, &loadF })
            WipeEachUse(b);
         for (Buf<char> *b : { &prevSlipS, &prevOnS, &farNow })
            WipeEachUse(b);
         for (Buf<int> *b : { &idxS, &partOff, &coordOff, &cUnk })
            WipeEachUse(b);
         WipeEachUse(&ie);
         WipeEachUse(&ci);
         WipeEachUse(&nearL);
         WipeEachUse(&tgS);
         WipeEachUse(&tkS);
         WipeEachUse(&tOnS);
         WipeEachUse(&tgES);
         WipeEachUse(&tgFS);
         WipeEachUse(&crossS);
         WipeEachUse(&szRec);
         WipeEachUse(&szCand);
         Size(8);
      }
      // The start of a use (an island's steps, a part's steps alone): no ball, shape or part yet, and nothing an earlier use wrote can be read
      void BeginUse()
      {
         nb = ns = ni = np = nf = nform = 0;
         WipeUseWork();
      }
      // VPX_GROUP_PHASE_TIMING=<file> (a diagnostic, read once; unset = off: no clock is read, nothing changes): one row per step, the wall time in each phase,
      // each phase's own (a phase entered inside another pauses it: a probe inside Eval counts as a probe, not as Eval). steady_clock (QueryPerformanceCounter on MSVC).
      enum Phase { phChoose, phBegin, phProbeEval, phProbeOther, phEvalRest, phLinear, phPartAfter, phAfter, phDiag, phN };
      struct PhaseClock
      {
         double us[phN] = {};
         long n[phN] = {};   // (the times each phase was entered)
         int cur = -1;
         std::chrono::steady_clock::time_point t0;
      };
      static FILE *PhaseFile()
      {
         static FILE *const f = []() -> FILE * {
            const char *const e = getenv("VPX_GROUP_PHASE_TIMING");
            FILE *const o = e != nullptr && *e ? fopen(e, "w") : nullptr;
            if (o != nullptr)
               fprintf(o, "t_ms,step,h_us,rounds,step_us,choose_us,begin_us,eval_probe_us,inrange_probe_us,eval_rest_us,linear_us,part_after_us,step_other_us,after_us,diag_us,attempts,evals,lin_solves,halvings,capsules,terms\n");
            return o;
         }();
         return f;
      }
      static PhaseClock &Clock() { static PhaseClock c; return c; }
      struct PhaseScope
      {
         int ph = -1, prev = -1;
         explicit PhaseScope(const int p)
         {
            if (PhaseFile() == nullptr)
               return;
            PhaseClock &c = Clock();
            const auto now = std::chrono::steady_clock::now();
            if (c.cur >= 0)
               c.us[c.cur] += std::chrono::duration<double, std::micro>(now - c.t0).count();
            prev = c.cur;
            ph = p;
            c.cur = p;
            c.n[p]++;
            c.t0 = now;
         }
         ~PhaseScope()
         {
            if (ph < 0)
               return;
            PhaseClock &c = Clock();
            const auto now = std::chrono::steady_clock::now();
            c.us[ph] += std::chrono::duration<double, std::micro>(now - c.t0).count();
            c.cur = prev;
            c.t0 = now;
         }
      };
      uint32_t phaseMs = 0;   // (the owner's game time, for the phase rows)
      long phaseStep = 0;

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
            h0s.Grow(n);
            ie.Grow(n);
            farNow.Grow(n);
            prevSlipS.Grow(n);
            prevOnS.Grow(n);
            szRec.Grow(n);
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
         overWhy = '?';
         stopBuf.clear();
         stopLines = 0;
         FILE *const pf = PhaseFile();
         while (done < dt - 1e-12)
         {
            const double hChoice = std::min(ChooseStep(), dt - done);
            double h;
            {
               PhaseScope ps(phDiag);   // (the owner's cap and its STEP_RULE row)
               h = capStep(hChoice);
            }
            double usIn[phN];
            const long rounds0 = books.rounds, halv0 = books.halvings, caps0 = envStats.capsules, terms0 = envStats.terms;
            std::chrono::steady_clock::time_point tStep;
            if (pf != nullptr)
            {
               for (int i = 0; i < phN; i++)
                  usIn[i] = Clock().us[i];
               tStep = std::chrono::steady_clock::now();
            }
            const double taken = Step(h);
            if (ClampFile() != nullptr)
               ClampFlush();
            double stepUs = 0.0, inStep = 0.0;
            if (pf != nullptr)
            {
               stepUs = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - tStep).count();
               for (int i = 0; i < phN; i++)
                  inStep += Clock().us[i] - usIn[i];
            }
            if (books.failed || taken <= 0.0 || (++guard > 3000 && GroupFixedStepS() <= 0.0))
            {
               overWhy = books.failed ? 'F' : taken <= 0.0 ? 'Z' : 'G';   // why the cycle was lost: F a failed step, Z a step of no length, G 3000 steps
               return false;
            }
            done += taken;
            {
               PhaseScope ps(phAfter);
               afterStep();
            }
            if (pf != nullptr)
            {
               PhaseClock &c = Clock();
               // (the counts: attempts = BeginAttempt calls, every solve of the step incl. failed ones, the halvings' and the stop search's; evals = Eval calls;
               // lin_solves = SolveLinear calls; capsules and terms: the hinge probes' capsule checks and union terms over the step)
               fprintf(pf, "%u,%ld,%.4f,%ld,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%ld,%ld,%ld,%ld,%ld,%ld\n", phaseMs, ++phaseStep, h * 1e6, books.rounds - rounds0, stepUs, c.us[phChoose], c.us[phBegin],
                       c.us[phProbeEval], c.us[phProbeOther], c.us[phEvalRest], c.us[phLinear], c.us[phPartAfter], stepUs - inStep, c.us[phAfter], c.us[phDiag],
                       c.n[phBegin], c.n[phEvalRest], c.n[phLinear], books.halvings - halv0, envStats.capsules - caps0, envStats.terms - terms0);
               for (int i = 0; i < phN; i++)
               {
                  c.us[i] = 0.0;
                  c.n[i] = 0;
               }
            }
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

      // What set the last ChooseStep (diagnostics for the offline trace; always kept): the rule (0 hMax, 1 normal stiffness, 3 arrival,
      // 4 turning, 5 floor 0.1 us, 7 a steady contact's normal stiffness at omegaHRest, 8 a part: partRule is its own), the contact (-1 = none) and that
      // contact's values at the choice
      struct StepDiag
      {
         int rule = 0, item = -1, partRule = -1;
         double delta = 0.0, dDot = 0.0, omN = 0.0;
         int unsteady = 0;                   // the steady test the contact failed, the first in order: 0 none (steady, or no test), 1 closing speed, 2 acceleration, 3 slip
         double unsteadyVal = 0.0, unsteadyLim = 0.0;   // that test's value and limit (m/s, m/s^2, m/s)
         int term = -1;                      // >= 0: the contact's term (Hit::t) whose arrival set the step (termArrival); delta and dDot are then that term's
      };
      mutable StepDiag diag;

      // The step choice from the state now (accuracy only)
      double ChooseStep() const
      {
         PhaseScope ps(phChoose);
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
            // termArrival: each term's closing rate by its own gradient, -v . n_t + slope_t * rate (the part's one coordinate, column 0)
            const bool terms = set.termArrival && hp.nt > 0 && nq > 0;
            const double qRate = terms ? part[shape[it.s]->part]->Coord(hp.c[0].coord).rate : 0.0;
            const auto termDot = [&](const HitTerm &ht) { return -A.v.Dot(ht.n) + ht.slope * qRate; };
            bool termNear = false;   // a term the longest step can reach though the union's depth cannot be
            for (int t = 0; terms && !termNear && t < hp.nt; t++)
               termNear = hp.t[t].d > 0.0 || -hp.t[t].d < std::max(termDot(hp.t[t]), set.vFloor) * set.hMax + 1e-9;
            if (hp.delta <= 0.0 && -hp.delta >= vc * set.hMax + 1e-9 && !termNear)
               continue;   // cannot be reached within the longest step
            Item L = it;
            if (!it.touching)
               TakeLaw(L, hp);
            if (invQ > 0.0)
               L.mRed = 1.0 / (1.0 / A.m + invQ);
            // the depth an approach at closing speed v reaches (all the closing energy into the squash): the stiffness while it lasts
            const auto reachDepth = [&](const double v) {
               const double eC = 0.5 * L.mRed * v * v;
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
               if (L.kb > 0.0)   // springs in series (the bend): the closing energy shared by them, the squash the Hertz part of it
               {
                  const double cS = 1.0 / L.kb;
                  double uLo = 0.0, uHi = lo;
                  for (int r = 0; r < 30 && uHi - uLo > 1e-3 * uHi; r++)
                  {
                     const double u = 0.5 * (uLo + uHi), f = L.kH * u * std::sqrt(u);
                     (0.4 * f * u + 0.5 * f * f * cS < eC ? uLo : uHi) = u;
                  }
                  const double fE = L.kH * uHi * std::sqrt(uHi);
                  lo = uHi + (L.kb > 0.0 ? fE / L.kb : 0.0);   // (the bend is in the depth; the part's give is in its solved angle, not in the depth)
               }
               return lo;
            };
            const double lo = reachDepth(vc);
            // the bounce's fastest rate at depth dE: m x'' + c x' + k x = 0 (c the law's own damping) decays at most at (c + sqrt(c^2 - 4 m k)) / 2m when overdamped,
            // else swings at omN. The mid-step damping takes a rate v to v (1 - z/2) / (1 + z/2), z = c h / m: on the spring's clock alone a damper past 2.56 times
            // critical reverses the rate in one step, a rebound the law has not (runs/2026-11-10-damped-clock/design.md)
            const auto fastRate = [&](const double dE, const double omN) {
               const double m = L.mRed, k = L.K(dE), c = L.DampC(dE, L.F(dE));
               return c * c > 4.0 * m * k ? (c + std::sqrt(c * c - 4.0 * m * k)) / (2.0 * m) : omN;
            };
            // an unsteady contact's clock at the expected depth dE: the bounce's fastest rate
            const auto unsteadyClock = [&](const double dE, double &omN) {
               omN = std::sqrt(L.K(dE) / L.mRed);
               const double r = fastRate(dE, omN);
               return r > 0.0 ? set.omegaH / r : set.hMax;
            };
            const double dE = std::max(hp.delta, lo);
            const double ktE = L.Kt(dE);
            const double omN = std::sqrt(L.K(dE) / L.mRed);
            // steady (Settings): by the motion, at the contact's own time scale. At the friction limit (Coulomb): a spring within mu N holds at most (mu N)^2 / (2 kt) of energy, so a slip with more than that
            // in 1/2 mt u^2 cannot be stopped by it: the contact slides, its force is mu N whatever the spring's stiffness, and the spring's time scale is not in the motion.
            bool steady = false;
            int why = 0;
            double whyVal = 0.0, whyLim = 0.0;
            if (set.vSteady > 0.0)
            {
               if (hp.delta > 0.0)
               {
                  V3 u = A.v + A.w.Cross(hp.n * (-(A.R - 0.5 * hp.delta)));
                  if (it.b >= 0)
                     u = u - (ball[it.b].v + ball[it.b].w.Cross(hp.n * (ball[it.b].R - 0.5 * hp.delta)));
                  else
                     u = u - vs;   // (a fixed surface: 0)
                  const V3 ut = u - hp.n * u.Dot(hp.n);
                  const double slipNow = ut.Length();
                  V3 aRel = fNet[it.a] * (1.0 / A.m);
                  if (it.b >= 0)
                     aRel = aRel - fNet[it.b] * (1.0 / ball[it.b].m);
                  const double omNow = std::sqrt(L.K(hp.delta) / L.mRed);
                  steady = std::fabs(dDot) <= set.vSteady && std::fabs(aRel.Dot(hp.n)) <= set.vSteady * omNow && (slipNow <= set.vSteady || (set.slideFree && slipNow * std::sqrt(ktE * L.mt) > L.mu * L.F(dE)));
                  if (!steady)
                  {
                     if (std::fabs(dDot) > set.vSteady)
                        why = 1, whyVal = dDot, whyLim = set.vSteady;
                     else if (std::fabs(aRel.Dot(hp.n)) > set.vSteady * omNow)
                        why = 2, whyVal = aRel.Dot(hp.n), whyLim = set.vSteady * omNow;
                     else
                        why = 3, whyVal = slipNow, whyLim = set.vSteady;
                  }
               }
               else
               {
                  steady = dDot <= set.vSteady;
                  if (!steady)
                     why = 1, whyVal = dDot, whyLim = set.vSteady;
               }
            }
            const double om = steady ? omN : fastRate(dE, omN);   // (a resting contact's rate is near 0 and the damping factor's size stays under 1: its spring's clock)
            double hs = om > 0.0 ? (steady ? set.omegaHRest : set.omegaH) / om : set.hMax;
            if (set.stepModel == 1 && steady && hp.delta > 0.0)   // by the motion in its mode (Settings::stepModel)
            {
               const double omR = std::sqrt(L.K(hp.delta) / L.mRed);
               V3 aRel = fNet[it.a] * (1.0 / A.m);
               if (it.b >= 0)
                  aRel = aRel - fNet[it.b] * (1.0 / ball[it.b].m);
               const double x = dDot / omR, y = aRel.Dot(hp.n) / (omR * omR), amp = std::sqrt(x * x + y * y);
               const double hm = amp > 0.0 ? std::sqrt(12.0 * set.epsRest / amp) / omR : set.hMax;
               hs = std::max(hs, std::min(hm, set.hMax));
               if (set.swingCap && 2.0 * amp <= set.epsRest)
                  hs = set.hMax;
            }
            const auto lim = [&](const double c, const int rule) {
               if (c < h)
               {
                  h = c;
                  dg.rule = rule;
                  dg.item = k;
                  dg.delta = hp.delta;
                  dg.dDot = dDot;
                  dg.omN = omN;
                  dg.unsteady = why;
                  dg.unsteadyVal = whyVal;
                  dg.unsteadyLim = whyLim;
                  dg.term = -1;
               }
            };
            if (hp.delta < 0.0 && dDot > 0.0)   // reached in about two steps, not overshot
            {
               const double ha = 0.5 * -hp.delta / dDot;
               lim(std::max(ha, hs), ha > hs ? 3 : 1);
            }
            else
               lim(hs, steady ? 7 : 1);
            // termArrival: each term not yet pressing is an approaching contact of its own: the arrival rule with the clock an untouched contact gets there
            // (unsteady, at the depth its own closing speed reaches), not the contact's (a resting contact's clock can be hMax)
            for (int t = 0; terms && t < hp.nt; t++)
            {
               const HitTerm &ht = hp.t[t];
               const double dDotT = termDot(ht);
               if (ht.d >= 0.0 || dDotT <= 0.0)
                  continue;
               const double ha = 0.5 * -ht.d / dDotT;
               if (ha >= h)
                  continue;   // (max(ha, hsT) >= h: no limit; the clock is not needed)
               double omNT = 0.0;
               const double hsT = unsteadyClock(reachDepth(std::max(dDotT, set.vFloor)), omNT);
               if (std::max(ha, hsT) < h)
               {
                  lim(std::max(ha, hsT), ha > hsT ? 3 : 1);
                  dg.term = t;
                  dg.delta = ht.d;
                  dg.dDot = dDotT;
                  dg.omN = omNT;
                  dg.unsteady = 1, dg.unsteadyVal = dDotT, dg.unsteadyLim = set.vSteady;
               }
            }
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
         const NearOff nearOff(nearOn);
         for (int p = 0; p < np; p++)   // what a part keeps over one step's attempts starts afresh (PartRecords)
            part[p]->WipeStepWork();
         double h = hTry;
         nHalveDiag = 0;
         failContactBuf.clear();
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
               // each stop on its own (one rule for every part): only the stops this step's move carries a body past count. A stop the body is on and leaving is none of
               // them: taken with the others it held the search's start at 0 and every step ended after ~1e-29 s (bench F4: 3000 steps, the cycle lost)
               crossS.assign(np, 0u);
               for (int p = 0; p < np; p++)
                  crossS[p] = part[p]->CrossingStops(x, partOff[p]);
               {
                  bool pinnedOne = false;
                  for (int p = 0; p < np; p++)
                  {
                     if (crossS[p] == 0u)
                        continue;
                     const bool onIt = part[p]->StopExcess(zero, partOff[p], crossS[p]) > -kOnStop;   // on a stop the move would carry it past
                     if (part[p]->StopExcess(x, partOff[p], crossS[p]) > 0.0 && onIt)   // (the solved move tested last: ArrivingAtStop and PinOnStop may use what it found)
                        pinnedOne = part[p]->PinOnStop(x, partOff[p]) || pinnedOne;
                  }
                  if (StopTraceOn() && stopLines < 12)   // diagnostics, off by default: each stop check that found a stop passed
                  {
                     char t[200];
                     snprintf(t, sizeof t, "  stop check h %.6g us pinned %d:", h * 1e6, pinnedOne ? 1 : 0);
                     stopBuf += t;
                     for (int p = 0; p < np; p++)
                     {
                        snprintf(t, sizeof t, " [p%d ex0 %.6g ex %.6g %s]", p, part[p]->StopExcess(zero, partOff[p]), part[p]->StopExcess(x, partOff[p]), part[p]->StopNote().c_str());
                        stopBuf += t;
                     }
                     stopBuf += "\n";
                     stopLines++;
                  }
                  if (pinnedOne)
                     continue;
               }
               double lo = 0.0, flo = std::min(StopExcessOf(zero), -1e-30), hi = h, fhi = StopExcessOf(x), hBelow = 0.0;
               int side = 0;
               bool solvedBelow = false;   // x holds the solve at hBelow (the search ended on it): no need to solve it again
               WarmFromSolve(h);
               for (int it = 0; it < 80; it++)
               {
                  double hm = (lo * fhi - hi * flo) / (fhi - flo);
                  if (!(hm > lo && hm < hi))
                     hm = 0.5 * (lo + hi);
                  searchTrace = SearchTraceOn();
                  searchBuf.clear();
                  const bool solvedHm = Solve(hm);
                  searchTrace = false;
                  if (!solvedHm)
                  {
                     if (SearchTraceOn())   // diagnostics, off by default: the failing stop-search solve, round by round
                        fprintf(stderr, "SEARCH FAIL why %c it %d hm %.6g us (lo %.6g hi %.6g h %.6g) stop excess lo %.4g hi %.4g\n%s", failWhy, it, hm * 1e6, lo * 1e6, hi * 1e6, h * 1e6, flo, fhi, searchBuf.c_str());
                     warmOn = false;
                     return EnergyStep(hTry);
                  }
                  WarmFromSolve(hm);
                  const double f = StopExcessOf(x);
                  if (f <= 0.0)
                  {
                     hBelow = hm;
                     bool arriving = true;
                     for (int p = 0; p < np; p++)   // (each searched part's last StopExcess was at this move)
                        arriving = arriving && (crossS[p] == 0u || part[p]->ArrivingAtStop());
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
               if (StopTraceOn() && stopLines <= 12)
               {
                  char t[120];
                  snprintf(t, sizeof t, "    -> hBelow %.6g us solvedBelow %d endOk %d\n", hBelow * 1e6, solvedBelow ? 1 : 0, endOk ? 1 : 0);
                  stopBuf += t;
               }
               if (!endOk)
                  return EnergyStep(hTry);
               Accept(hBelow, true);
               books.stopSteps++;
               return hBelow;
            }
            if (nHalveDiag < kHalveDiag)   // why this attempt failed (diagnostics)
            {
               if (RoundTraceOn() && nHalveDiag == 0)
                  failRoundBuf = roundBuf;
               if (ContactTraceNow() && nHalveDiag == 0)
                  failContactBuf = contactBuf;
               HalveDiag &hd = halveDiag[nHalveDiag++];
               hd = HalveDiag();   // (slots are reused: counts are this attempt's)
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
      // why the last Step's attempts failed and were halved: P a part could not begin the step, B a part's forces could not be formed (the band's path would not
      // settle), G a first guess the layer check refused 60 times, N a residual not finite, L the linear solve
      // singular, R a Newton move the layer check refused 60 times, C not solved within maxRounds, E the end evaluation not finite
      struct HalveDiag { char why = '?'; double h = 0.0; int rounds = 0; double moveC = 0.0, moveW = 0.0; int slipFlips = 0, touchFlips = 0, maxSliding = 0; double lastW[6] = {}; };
      static constexpr int kHalveDiag = 8;
      HalveDiag halveDiag[kHalveDiag];
      std::string roundBuf, failRoundBuf;   // VPX_GROUP_ROUND_TRACE: this attempt's rounds; the first failed attempt of the step
      mutable std::string probeTrace;       // VPX_GROUP_ROUND_TRACE: the carried surfaces' probes of the last Eval (end and mid, their pieces)
      // VPX_GROUP_JAC_WINDOW=<t0>-<t1> (ms; with VPX_GROUP_JAC_MODE=1; a diagnostic, read once; unset = off): for solves that converge in that window, each round's
      // residual and evaluations, and every Jacobian entry (balls and parts) where the solve's own (a) and central differences (f) differ by more than 10 % and 1e-6
      // (jacWinBuf, written by the owner; 600 lines in all)
      static bool JacWindow(const uint32_t ms)
      {
         static const double *const w = []() -> const double * {
            static double v[2];
            const char *const e = getenv("VPX_GROUP_JAC_WINDOW");
            return e != nullptr && sscanf(e, "%lf-%lf", &v[0], &v[1]) == 2 ? v : nullptr;
         }();
         return w != nullptr && ms >= w[0] && ms <= w[1];
      }
      std::string jacWinBuf, jacWinTry;
      long jacWinLines = 0, evalCount = 0, evalCountRound = 0, jacWinTryEntries = 0;   // (jacWinLines: Jacobian entry lines written, 2000 in all)
      bool jacWinDump = false;   // the main Eval of a round in the window: each contact's state into jacWinCon (the last round's is kept with the converged solve)
      std::string jacWinCon;
      // VPX_GROUP_CLAMP_TRACE=<file> (diagnostics, read once; unset = off: nothing is formed or written): each Newton round of a step that has a carried surface's
      // contact (a wire), the round's point (move or retake) and verdict (kept, rejected = taken back shorter, done), and each such contact's switches; a step is
      // written only if its attempts took more than 10 rounds in all or one failed
      static FILE *ClampFile()
      {
         static FILE *const f = []() -> FILE * {
            const char *const e = getenv("VPX_GROUP_CLAMP_TRACE");
            return e != nullptr && *e ? fopen(e, "w") : nullptr;
         }();
         return f;
      }
      bool clampDump = false;            // the main Eval of a round: each carried surface's contact into clampCon
      std::string clampCon, clampStep;   // this round's contact lines; this step's rounds
      int clampAttempt = 0, clampRounds = 0, clampFails = 0;
      bool clampRetake = false;          // the next round's point is a line-search retake
      std::string UnkName(const int u) const
      {
         if (u < 6 * nb)
         {
            static const char *const ax[6] = { "x", "y", "z", "wx", "wy", "wz" };
            return "b" + std::to_string(u / 6) + "." + ax[u % 6];
         }
         for (int c = 0; c < nCoord; c++)
            if (cUnk[c] == u)
               return "q" + std::to_string(c);
         return "u" + std::to_string(u);
      }
      // diagnostics, off by default (VPX_GROUP_CONTACT_TRACE=1): per round of an attempt, each contact's depth, normal, push, friction against mu N, stick or slide, and each part's
      // state and row terms (Part::ContactNote); the first failed attempt of a step kept in failContactBuf, printed for the first kContactTraceSteps energy steps of a run (contactsystem.cpp)
      static bool ContactTraceOn() { static const bool on = []() { const char *const e = getenv("VPX_GROUP_CONTACT_TRACE"); return e != nullptr && atoi(e) != 0; }(); return on; }
      static constexpr int kContactTraceSteps = 3;
      static int &ContactTraceLeft() { static int n = kContactTraceSteps; return n; }   // steps still to print
      static bool ContactTraceNow() { return ContactTraceOn() && ContactTraceLeft() > 0; }
      struct ContactDiag { V3 n; double dm = 0.0, N = 0.0, fe = 0.0, fd = 0.0, ft = 0.0, ftTr = 0.0, lim = 0.0, ut = 0.0, qN = 0.0, qF = 0.0; int nq = 0; };   // an Eval's values, copied (never read back by the step)
      std::vector<ContactDiag> cDiag;
      std::string contactBuf, failContactBuf;
      // diagnostics, off by default (VPX_GROUP_HISTORY_TRACE=1): a group with a part keeps its last kHistory accepted steps (Accept), one line each: the energy over the step against its
      // books' change, each part, each ball, each contact touching or in play; printed with the run's first energy-step warning (contactsystem.cpp). Read only
      static bool HistoryTraceOn() { static const bool on = []() { const char *const e = getenv("VPX_GROUP_HISTORY_TRACE"); return e != nullptr && atoi(e) != 0; }(); return on; }
      static constexpr int kHistory = 50;
      std::string history[kHistory];
      long historyN = 0;          // steps kept so far (the newest at (historyN - 1) % kHistory)
      unsigned histMs = 0;        // the caller's table time at the cycle's start (ms) and how far the cycle is (s): the steps' times
      double histClock = 0.0;
      bool histInEnergy = false;  // an Accept inside EnergyStep (its move may be restored)
      long histCycle = 0;         // the caller's cycle count (a new cycle: the group was gathered again from the table)
      // the group's energy as Energy() sums it, by where it sits (J): the balls' motion, their place (gravity and fields), the contacts' springs, each part, the pair form
      static constexpr int kHistParts = 8;
      struct HistE { double tot = 0.0, ke = 0.0, pe = 0.0, con = 0.0, form = 0.0, parts[kHistParts] = {}; int nb = 0, np = 0, ni = 0, touching = 0; long cycle = -1, eSteps = -1; size_t sig = 0; bool valid = false; };
      HistE HistEnergy() const
      {
         HistE s;
         for (int i = 0; i < nb; i++)
         {
            const Ball &b = ball[i];
            s.ke += 0.5 * b.m * b.v.Dot(b.v) + 0.5 * b.I() * b.w.Dot(b.w);
            s.pe += -b.m * set.g.Dot(b.c) + (nf > 0 && b.fieldScale != 0.0 ? FieldEnergy(b) : 0.0);
         }
         for (int k = 0; k < ni; k++)
            if (item[k].touching)
            {
               Hit hp;
               Probe(item[k], ball, hp);
               s.con += item[k].V(hp.delta) + 0.5 * item[k].Kt(hp.delta) * item[k].xi.Dot(item[k].xi);
               if (set.rollModel == 1)
                  s.con += 0.5 * item[k].Kr(hp.delta) * item[k].ra.th.Dot(item[k].ra.th) + 0.5 * item[k].Ktw(hp.delta) * item[k].ra.tw * item[k].ra.tw;
               s.touching++;
            }
         double ep = 0.0;
         for (int p = 0; p < np; p++)
         {
            const double e = part[p]->Energy(ball, nb);
            ep += e;
            if (p < kHistParts)
               s.parts[p] = e;
            s.sig = s.sig * 1000003u + (size_t)(const void *)part[p];   // which parts: the group's name across cycles
         }
         s.form = FormEnergy(ball);
         s.tot = s.ke + s.pe + s.con + ep + s.form;
         s.nb = nb;
         s.np = np;
         s.ni = ni;
         s.cycle = histCycle;
         s.eSteps = books.energySteps;
         s.valid = true;
         return s;
      }
      HistE histLast[kHistParts];   // the last accepted step's end of each group seen (by sig), oldest replaced
      int histLastNext = 0;
      std::string HistoryDump() const
      {
         std::string o;
         for (long i = std::max(0L, historyN - kHistory); i < historyN; i++)
            o += history[i % kHistory] + "\n";
         return o;
      }
      void HistoryKeep(const double h, const HistE &s0, const Books &bk0)
      {
         char t[700];
         const HistE s1 = HistEnergy();
         const double dDrive = books.driveJ - bk0.driveJ, dLoss = books.Losses() - bk0.Losses(), E0 = s0.tot, E1 = s1.tot;
         // the gap: this step's start against the same group's last accepted end (energy that came in between steps or cycles, outside any step's books), by where it sits
         HistE *prev = nullptr;
         for (HistE &l : histLast)
            if (l.valid && l.sig == s0.sig)
               prev = &l;
         std::string gap;
         if (prev == nullptr)
            gap = " | gap: first step of this group seen";
         else
         {
            const HistE &a = *prev;
            int g = snprintf(t, sizeof t, " | gap %.5g uJ (ballKE %.4g ballPE %.4g contacts %.4g form %.4g", (s0.tot - a.tot) * 1e6, (s0.ke - a.ke) * 1e6, (s0.pe - a.pe) * 1e6, (s0.con - a.con) * 1e6, (s0.form - a.form) * 1e6);
            for (int p = 0; p < std::min(s0.np, kHistParts) && g < 600; p++)
               g += snprintf(t + g, sizeof t - g, " p%d %.4g", p, (s0.parts[p] - a.parts[p]) * 1e6);
            snprintf(t + g, sizeof t - g, ")%s%s%s nb %d->%d ni %d->%d touching %d->%d", s0.cycle != a.cycle ? " NEW-CYCLE" : " same-cycle", s0.eSteps != a.eSteps ? " ENERGY-STEP-BETWEEN" : "",
                     a.np != s0.np ? " PARTS-CHANGED" : "", a.nb, s0.nb, a.ni, s0.ni, a.touching, s0.touching);
            gap = t;
         }
         if (prev == nullptr)
         {
            prev = &histLast[histLastNext];
            histLastNext = (histLastNext + 1) % kHistParts;
         }
         *prev = s1;
         int n = snprintf(t, sizeof t, "  #%ld grp %04zx t %u ms + %.3f us h %.4g us%s | E0 %.9g E1 %.9g uJ dE %.5g drive %.5g loss %.5g rise %.5g uJ (normal %.4g slip %.4g shear %.4g twist %.4g roll %.4g part %.4g brake %.4g)",
                          historyN, s0.sig & 0xffff, histMs, histClock * 1e6, h * 1e6, histInEnergy ? " ENERGY-STEP" : "", E0 * 1e6, E1 * 1e6, (E1 - E0) * 1e6, dDrive * 1e6, dLoss * 1e6, (E1 - E0 - dDrive + dLoss) * 1e6,
                          (books.normalLossJ - bk0.normalLossJ) * 1e6, (books.slipJ - bk0.slipJ) * 1e6, (books.shearDampJ - bk0.shearDampJ) * 1e6, (books.twistJ - bk0.twistJ) * 1e6,
                          (books.rollJ - bk0.rollJ) * 1e6, (books.partLossJ - bk0.partLossJ) * 1e6, (books.spinBrakeJ - bk0.spinBrakeJ) * 1e6);
         std::string line(t, std::max(0, std::min(n, (int)sizeof t - 1)));
         line += gap;
         for (int p = 0; p < np; p++)
            line += " | p" + std::to_string(p) + " " + part[p]->HistoryNote();
         for (int i = 0; i < nb; i++)
         {
            const Ball &b = ball[i];
            snprintf(t, sizeof t, " | b%d c %.4f %.4f %.4f mm speed %.5g m/s KEt %.5g KEw %.5g uJ", i, b.c.x * 1e3, b.c.y * 1e3, b.c.z * 1e3, b.v.Length(), 0.5 * b.m * b.v.Dot(b.v) * 1e6, 0.5 * b.I() * b.w.Dot(b.w) * 1e6);
            line += t;
         }
         for (int k = 0; k < ni; k++)
         {
            const Item &it = item[k];
            const ItemEnd &e = ie[k];
            if (!it.touching && !e.active)
               continue;
            Hit hp;
            Probe(it, ball, hp);
            const double vn = it.touching ? it.V(hp.delta) : 0.0, vs = it.touching ? 0.5 * it.Kt(hp.delta) * it.xi.Dot(it.xi) : 0.0;
            const double vr = it.touching && set.rollModel == 1 ? 0.5 * it.Kr(hp.delta) * it.ra.th.Dot(it.ra.th) + 0.5 * it.Ktw(hp.delta) * it.ra.tw * it.ra.tw : 0.0;
            const bool cd = e.active && (int)cDiag.size() >= ni;
            snprintf(t, sizeof t, " | c%d b%d %s%d t%d d %.5g um (solve d1 %.5g zero %.4g) n_z %.4f N %.5g Ft %.4g muN %.4g %s V %.5g Vs %.4g Vr %.4g uJ", k, it.a, it.b >= 0 ? "ball b" : "s", it.b >= 0 ? it.b : it.s,
                     it.touching ? 1 : 0, hp.delta * 1e6, e.d1 * 1e6, it.dZero * 1e6, hp.n.z, cd ? cDiag[k].N : 0.0, cd ? cDiag[k].ft : 0.0, cd ? cDiag[k].lim : 0.0, e.slipping ? "SLIDE" : "stick", vn * 1e6, vs * 1e6, vr * 1e6);
            line += t;
         }
         history[historyN % kHistory] = line;
         historyN++;
      }
      char overWhy = '?';                   // why StepsOver lost its cycle (F, Z, G)
      std::string stopBuf;                  // VPX_GROUP_STOP_TRACE=1: this cycle's stop checks (the first 12)
      int stopLines = 0;
      static bool StopTraceOn() { static const bool on = []() { const char *const e = getenv("VPX_GROUP_STOP_TRACE"); return e != nullptr && atoi(e) != 0; }(); return on; }
      // after an energy step: the balls and contacts it saved at the step's start, and the warm start then in hand, full precision (for an exact replay; m, m/s, rad/s)
      std::string StartNote() const
      {
         std::string msg;
         char z[700];
         for (int i = 0; i < nb; i++)
         {
            const Ball &b = esB0[i];
            snprintf(z, sizeof z, " [ball %d c %.17g %.17g %.17g v %.17g %.17g %.17g w %.17g %.17g %.17g m %.17g R %.17g fExt %.17g %.17g %.17g a %.17g %.17g %.17g alpha %.17g %.17g %.17g warmA %.17g %.17g %.17g warmAl %.17g %.17g %.17g warmOn %d]",
               i, b.c.x, b.c.y, b.c.z, b.v.x, b.v.y, b.v.z, b.w.x, b.w.y, b.w.z, b.m, b.R, b.fExt.x, b.fExt.y, b.fExt.z, esA0[i].x, esA0[i].y, esA0[i].z,
               esAl0[i].x, esAl0[i].y, esAl0[i].z, warmA[i].x, warmA[i].y, warmA[i].z, warmAl[i].x, warmAl[i].y, warmAl[i].z, (int)warmOn);
            msg += z;
         }
         for (int k = 0; k < ni; k++)
         {
            const Item &it = esI0[k];
            if (!it.touching && it.xi.x == 0.0 && it.xi.y == 0.0 && it.xi.z == 0.0)
               continue;
            snprintf(z, sizeof z, " [item %d a %d b %d s %d touching %d xi %.17g %.17g %.17g th %.17g %.17g %.17g tw %.17g dZero %.17g]",
               k, it.a, it.b, it.s, (int)it.touching, it.xi.x, it.xi.y, it.xi.z, it.ra.th.x, it.ra.th.y, it.ra.th.z, it.ra.tw, it.dZero);
            msg += z;
         }
         return msg;
      }
      // back: a line-search round (the last move retaken shorter); grip: the parts' grip change after the move (m/s)
      struct RoundDiag { double moveC = 0.0, moveW = 0.0, scale = 1.0, grip = 0.0; int nActive = 0, nSlip = 0, nSlipFlip = 0, nOnFlip = 0, back = 0; };
      static constexpr int kRounds = 64;
      RoundDiag roundDiag[kRounds];
      // Settings::jacMode: the solve's own Jacobian and the central-difference one at the last round of the last solve (rows: per ball the momentum residual 3,
      // then the spin residual 3; columns: per ball the move 3, then the spin change 3)
      Jac JacAna() { return Jac(jacAnaS, jStride); }
      Jac JacFd() { return Jac(jacFdS, jStride); }
      int jacRound = 0;

      // Test tooling, off by default (VPX_TEST_IMPULSES=1; read only: no value of a step depends on it): the impulse each force of the group gives each ball over the
      // accepted steps, h x the force and h x the torque about the ball's centre of the Eval the step was accepted on (implicit midpoint: m (v1 - v0) = h x their sum, to the
      // solve's tolerance). Keyed by group ball and source: a contact's shape Name() ("shape<s>" when it has none), "#<j>" for the pair with group ball j, "gravity",
      // "fext" (the held outside force), "field", "pairform", "brake" (the spin brake's torque), "part<p>" (a part's own force on a ball, the cord's), "energystep" (an
      // energy step's whole momentum change less gravity's: it has no force law). Summed over a use; the owner clears it (ImpClear) and keeps it only when it keeps the use's result.
      static bool ImpOn() { static const bool on = []() { const char *const e = getenv("VPX_TEST_IMPULSES"); return e != nullptr && atoi(e) != 0; }(); return on; }
      struct ImpSum { V3 J, L; };
      std::map<std::pair<int, std::string>, ImpSum> imp;
      void ImpClear() { imp.clear(); }
      void ImpAdd(const int b, const std::string &src, const V3 &J, const V3 &L)
      {
         if (J.x == 0.0 && J.y == 0.0 && J.z == 0.0 && L.x == 0.0 && L.y == 0.0 && L.z == 0.0)
            return;
         ImpSum &s = imp[std::make_pair(b, src)];
         s.J = s.J + J;
         s.L = s.L + L;
      }

   private:
      // ImpOn: the last Eval's forces by source. Per ball: gravity, fExt, the fields, the pair form (forces), the brake (torque); per item k, impIt[4k ..]: its force and
      // torque on ball a, then on ball b (a pair); per part and ball (p * nb + i): the part's own force and torque
      std::vector<V3> impG, impX, impFld, impForm, impBrake, impIt, impPF, impPM;
      void ImpAccept(const double h)
      {
         if (impG.size() < (size_t)nb || impIt.size() < 4 * (size_t)ni || impPF.size() < (size_t)np * nb)   // (no Eval of this shape before: a gap the check shows)
            return;
         for (int i = 0; i < nb; i++)
         {
            ImpAdd(i, "gravity", impG[i] * h, V3());
            ImpAdd(i, "fext", impX[i] * h, V3());
            ImpAdd(i, "field", impFld[i] * h, V3());
            ImpAdd(i, "pairform", impForm[i] * h, V3());
            ImpAdd(i, "brake", V3(), impBrake[i] * h);
         }
         for (int k = 0; k < ni; k++)
         {
            if (farNow[k])
               continue;
            const Item &it = item[k];
            const V3 *const o = &impIt[4 * (size_t)k];
            std::string name;
            if (it.b >= 0)
               name = "#" + std::to_string(it.b);
            else
            {
               if (shape[it.s] != nullptr)
                  name = shape[it.s]->Name();
               if (name.empty())
                  name = "shape" + std::to_string(it.s);
            }
            ImpAdd(it.a, name, o[0] * h, o[1] * h);
            if (it.b >= 0)
               ImpAdd(it.b, "#" + std::to_string(it.a), o[2] * h, o[3] * h);
         }
         for (int p = 0; p < np; p++)
            for (int i = 0; i < nb; i++)
               ImpAdd(i, "part" + std::to_string(p), impPF[(size_t)p * nb + i] * h, impPM[(size_t)p * nb + i] * h);
      }
      // an energy step's momentum change (EnergyStepCore): its own entry, the impulses of the Accepts inside it dropped (some are restored)
      double EnergyStep(const double hTry)
      {
         if (!ImpOn())
            return EnergyStepCore(hTry);
         const std::map<std::pair<int, std::string>, ImpSum> imp0 = imp;
         std::vector<V3> v0((size_t)nb), w0((size_t)nb);
         for (int i = 0; i < nb; i++)
         {
            v0[i] = ball[i].v;
            w0[i] = ball[i].w;
         }
         const double hs = EnergyStepCore(hTry);
         imp = imp0;
         for (int i = 0; i < nb; i++)
         {
            const V3 G = set.g * ball[i].m * hs;
            ImpAdd(i, "gravity", G, V3());
            ImpAdd(i, "energystep", (ball[i].v - v0[i]) * ball[i].m - G, (ball[i].w - w0[i]) * ball[i].I());
         }
         return hs;
      }

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
      Buf<double> d0s;
      Buf<Hit> h0s;   // the start's terms (hinge_piece_sum; nt 0 otherwise)
      Buf<Ball> b1r, bmr;        // InRange's moved balls
      Buf<ItemEnd> ie;
      Buf<double> x;             // the unknowns: per ball the move (3) and the spin change (3), then the parts' own
      Buf<V3> aLast, alphaLast;  // each ball's acceleration and spin's over the last step (the next step's first guess)
      // what each function of a step works in
      mutable Buf<V3> fNetS;
      Buf<V3> evF, evM, evDc, evDw, evWm, esA0, esAl0;
      // the stop search's first guess: the accelerations of this step's last solve (same start state), so a re-solve at a nearby length starts near its answer
      Buf<V3> warmA, warmAl;
      void WarmFromSolve(double)
      {
         for (int i = 0; i < nb; i++)
         {
            // the force and torque the solve used (Eval's at its move), over the mass and inertia: the mean accelerations of that move (m dv = F h at a solution),
            // bounded by the forces however short h is (from the move itself, (dc - v h) 2 / h^2 is the solve's tolerance over h^2 when h is near zero)
            warmA[i] = evF[i] * (1.0 / ball[i].m);
            warmAl[i] = evM[i] * (1.0 / ball[i].I());
         }
         warmOn = true;
      }
      Buf<Ball> evB1, evBm, farB1, farBm, esB0;
      Buf<Item> esI0;
      Buf<double> resS, stpS, xPrevS, stpPrevS, trialS, zeroS, esRes, xpS, rpS, rmS, jStore, jacAnaS, jacFdS;
      Buf<char> prevSlipS, prevOnS;
      Buf<int> idxS;
      int jStride = 0;
      Buf<int> partOff;
      // the parts' coordinates, one list (coordOff: each part's first): as they are this attempt, the unknown of each (-1 = none), the moves a probe is at, the balls' reaction on each and on each part in this evaluation
      Buf<int> coordOff, cUnk;
      Buf<CoordInfo> ci;
      Buf<double> qEnd, qMid, qF;
      Buf<V3> loadF;
      static constexpr int kGradMax = 3 + 7;   // a discrete gradient's entries: a ball's centre and a hit's columns (two balls: 6)
      // stick_zone_trial: each contact's linear model of its trial force along a move, written by a round's Eval (with the Jacobian) at its point
      struct StickRec
      {
         bool on = false;            // the shear patch slides at this point with N > 0
         V3 a, rA, rB, nm;           // the trial force kS xi0 + C ut; the levers; the normal
         double C = 0.0, mu = 0.0, N0 = 0.0;
         int A = -1, B = -1, nq = 0, nn = 0;
         V3 cV[kGradMax];            // a carried surface's column's contact-point velocity per unit rate
         int qU[kGradMax];           // its unknown (-1: a given move)
         double nG[kGradMax];        // N's slope on the unknown nIx (dN = sum nG d[nIx])
         int nIx[kGradMax];
      };
      Buf<StickRec> szRec;
      std::vector<double> szCand;
      Buf<char> farNow;          // not touching and out of reach in this attempt
      // [Contact] hinge_piece_sum: per contact, the pieces of its surface the ball can reach in this attempt (Shape::NearPieces), built at its start; every probe
      // of the attempt (start, middle, end, every Newton try) takes the same list. A piece off it has depth <= 0 at every place of the step (checked after the
      // solve: NearStillNear), so the union's sum, forces and energy are those of all its pieces. Outside an attempt (the step's choice, events) a probe takes all
      std::vector<std::vector<int>> nearL;
      std::vector<std::array<double, 4>> tgS;   // (Eval, hinge_piece_sum: each term's gradient, its chord slope and whether it presses; kept to save the allocation)
      std::vector<double> tkS;
      std::vector<char> tOnS;
      std::vector<std::array<double, 4>> tgES;   // (... each term's END gradient, the column side of its slope: endGradSlope)
      std::vector<std::array<double, 4>> tgFS;   // (... each term's force direction gdT, the row side of its slope: termForceRows)
      struct NearOff { bool &f; explicit NearOff(bool &b) : f(b) { f = false; } ~NearOff() { f = false; } };   // a step's lists end with it (Step, EnergyStep)
      static bool NearAllOn() { static const bool on = []() { const char *const e = getenv("VPX_HINGE_NEAR_ALL"); return e != nullptr && atoi(e) != 0; }(); return on; }   // diagnostics: every list = all pieces (A/B)
      // the lists of this attempt (after the parts' BeginStep: their rates now). The reach: the ball's travel 2 |v| h and the part's turn 2 |rate| h (times nearGrow),
      // and never less than the far test's bound (OutOfReach), so a contact that is not far has its deepest piece listed and its gap, normal and law are those of all
      // its pieces. The far test's second solve (farOff) and the diagnostic VPX_HINGE_NEAR_ALL take every piece
      void BuildNear(const double h)
      {
         nearOn = false;
         if (nearL.size() < (size_t)ni)
            nearL.resize((size_t)ni);
         const bool all = farOff || NearAllOn() || nearGrow > 1e6;
         for (int k = 0; k < ni; k++)
         {
            Item &it = item[k];
            it.nearK = -1;
            if (it.b >= 0 || shape[it.s]->part < 0)
               continue;
            const Shape *const sh = shape[it.s];
            const Part *const pp = part[sh->part];
            const Ball &A = ball[it.a];
            double rate = 0.0, sp = A.v.Length();
            for (int c = 0; c < pp->NumCoords(); c++)
            {
               const CoordInfo q = pp->Coord(c);
               rate += std::fabs(q.rate);
               sp += std::fabs(q.rate) * q.reach;
            }
            const double reach = all ? INFINITY : 2.0 * A.v.Length() * h * nearGrow + 1e-4, swing = all ? INFINITY : 2.0 * rate * h * nearGrow;
            if (sh->NearPieces(A.c, A.R, reach, swing, 2.0 * sp * h + 1e-4, nearL[(size_t)k]))
               it.nearK = k;
         }
         nearOn = true;
      }
      // after the solve: every piece off a list is off its ball at the move's middle and end (the start: by its building)
      bool NearStillNear(const double *xv)
      {
         if (carried)
            SetTrials(xv);
         for (int k = 0; k < ni; k++)
         {
            const Item &it = item[k];
            if (it.nearK < 0 || farNow[k])
               continue;
            const Shape *const sh = shape[it.s];
            const Ball &A = ball[it.a];
            const V3 dc(xv[6 * it.a], xv[6 * it.a + 1], xv[6 * it.a + 2]);
            const double *const qe = (const double *)qEnd + coordOff[sh->part], *const qm = (const double *)qMid + coordOff[sh->part];
            if (!sh->NearHolds(A.c + dc, A.R, qe[0], nearL[(size_t)k]) || !sh->NearHolds(A.c + dc * 0.5, A.R, qm[0], nearL[(size_t)k]))
               return false;
         }
         return true;
      }
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
         PhaseScope ps(Clock().cur == phEvalRest ? phProbeEval : phProbeOther);
         const Shape *const sh = it.b < 0 ? (const Shape *)shape[it.s] : nullptr;
         if (sh != nullptr && sh->part >= 0)
         {
            sh->trial = (const double *)qEnd + coordOff[sh->part];
            Probe(it, b1, h1);
            if (RoundTraceOn())   // diagnostics: each probe's pieces, for the round trace (Eval clears it)
               probeTrace += "    probe p" + std::to_string(sh->part) + " ball " + std::to_string(it.a) + " s" + std::to_string(it.s) + " end: " + part[sh->part]->LastProbeNote() + "\n";
            sh->trial = (const double *)qMid + coordOff[sh->part];
            Probe(it, bm, hm);
            if (RoundTraceOn())
               probeTrace += "    probe p" + std::to_string(sh->part) + " ball " + std::to_string(it.a) + " s" + std::to_string(it.s) + " mid: " + part[sh->part]->LastProbeNote() + "\n";
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
      std::vector<unsigned> crossS;   // per part, the stops the step's move carries it past (the stop search follows these only)
      double StopExcessOf(const double *xv) const
      {
         double e = -INFINITY;
         for (int p = 0; p < np; p++)
            if (crossS[p] != 0u)
               e = std::max(e, part[p]->StopExcess(xv, partOff[p], crossS[p]));
         return e;
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
         hp.ktb = 0.0;
         hp.nc = 0;   // (hits are kept and reused: a fixed surface writes no columns, and must not inherit a carried one's)
         if (it.b < 0)
         {
            const Shape *const sh = shape[it.s];
            const bool listed = nearOn && it.nearK >= 0 && sh->part >= 0;
            if (listed)
               part[sh->part]->nearNow = &nearL[(size_t)it.nearK];
            sh->Probe(A.c, A.R, hp);
            if (listed)
               part[sh->part]->nearNow = nullptr;
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
         it.ktb = it.b >= 0 ? 0.0 : hp.ktb;
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

      // Gonzalez: gm + (df - gm . dq) / (dq^T W dq) W dq, so that the result . dq = df exactly; a mismatch within the rounding of f is left out. iM: each coordinate's
      // inverse mass or inertia (0: it takes no correction). byMass: W = M / M_0, the kinetic-energy metric (one unit for metres and radians, so each body's speed changes
      // by its own move; M_0 the first coordinate's, so the weights are mass ratios: runs/2026-11-10-gonzalez-metric/design.md). Else W = iM: a light part took a heavy
      // body's mismatch as a false force (the 15441 ms halvings on a wire switch)
      static void Gonzalez(const double *gm, double df, double dfNoise, const double *dq, const double *iM, int n, bool byMass, double *out, bool *cfZero = nullptr)
      {
         const auto W = [&](int a) { return !byMass ? iM[a] : iM[a] > 0.0 ? iM[0] / iM[a] : 0.0; };
         double gdq = 0.0, qq = 0.0;
         for (int a = 0; a < n; a++)
         {
            gdq += gm[a] * dq[a];
            qq += dq[a] * W(a) * dq[a];
         }
         const double cf = qq > 1e-300 && std::fabs(df - gdq) > dfNoise ? (df - gdq) / qq : 0.0;
         for (int a = 0; a < n; a++)
            out[a] = gm[a] + cf * W(a) * dq[a];
         if (cfZero != nullptr)   // (VPX_GROUP_CLAMP_TRACE: the correction left out)
            *cfZero = cf == 0.0;
      }
      // The same elimination over the nonzero entries of each pivot row: a group's bodies touch few others
      static bool SolveLinear(const Jac A, double *b, int n, int *const idx)
      {
         PhaseScope ps(phLinear);
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
         PhaseScope ps(phEvalRest);
         evalCount++;
         if (RoundTraceOn())
            probeTrace.clear();
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
         const bool impOn = ImpOn();   // (test tooling: the forces by source, read only)
         if (impOn)
         {
            for (std::vector<V3> *const v : { &impG, &impX, &impFld, &impForm, &impBrake })
               v->assign((size_t)nb, V3());
            impIt.assign(4 * (size_t)ni, V3());
            impPF.assign((size_t)np * nb, V3());
            impPM.assign((size_t)np * nb, V3());
            for (int i = 0; i < nb; i++)
            {
               impG[i] = set.g * ball[i].m;
               impX[i] = ball[i].fExt;
            }
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
         if (impOn)
            for (int i = 0; i < nb; i++)
            {
               impBrake[i] = M[i];    // (M starts at zero: the brake is its only term so far)
               impFld[i] = F[i];      // (the fields' and the pair form's: by the change they make to F)
            }
         // the fields' stored energy: its discrete gradient over the ball's move (Gonzalez), so the work it does is the exact difference at any step length
         for (int i = 0; nf > 0 && i < nb; i++)
            if (ball[i].fieldScale != 0.0)
            {
               V3 gf;
               M3 hf;
               const double e1 = FieldEnergy(b1[i]);
               FieldEnergy(bm[i], &gf, J ? &hf : nullptr);
               const double gm[3] = { gf.x, gf.y, gf.z }, dq[3] = { dc[i].x, dc[i].y, dc[i].z }, iM[3] = { 1.0 / ball[i].m, 1.0 / ball[i].m, 1.0 / ball[i].m };
               double gd[3];
               Gonzalez(gm, e1 - fieldE0[i], 16.0 * 2.220446049250313e-16 * std::fabs(e1) + gf.Length() * noiseD, dq, iM, 3, set.gonzalezMass, gd);
               F[i] = F[i] - V3(gd[0], gd[1], gd[2]);
               if (J)
                  AddBlock(J, 6 * i, 6 * i, hf, 0.5 * h);   // (the midpoint moves half as far as the end)
            }
         if (impOn)
            for (int i = 0; i < nb; i++)
            {
               impFld[i] = F[i] - impFld[i];
               impForm[i] = F[i];
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
            const double iM[6] = { 1.0 / ball[T.a].m, 1.0 / ball[T.a].m, 1.0 / ball[T.a].m, two ? 1.0 / ball[T.b].m : 0.0, two ? 1.0 / ball[T.b].m : 0.0, two ? 1.0 / ball[T.b].m : 0.0 };
            double gd[6];
            Gonzalez(gm, q1 - q0, 16.0 * 2.220446049250313e-16 * std::fabs(q1) + noiseD, dq, iM, two ? 6 : 3, set.gonzalezMass, gd);
            F[T.a] = F[T.a] + V3(gd[0], gd[1], gd[2]) * T.g;
            if (two)
               F[T.b] = F[T.b] + V3(gd[3], gd[4], gd[5]) * T.g;
         }
         if (impOn)
            for (int i = 0; i < nb; i++)
               impForm[i] = F[i] - impForm[i];

         for (int k = 0; k < ni; k++)
         {
            const Item &it = item[k];
            // impOn: this contact's force and torque on its balls, by the change it makes to F and M (taken on every way out of the loop's body)
            struct ImpItem
            {
               const V3 *F, *M;
               V3 *o;
               int a, b;
               V3 f0a, m0a, f0b, m0b;
               ImpItem(const V3 *F_, const V3 *M_, V3 *o_, const int a_, const int b_) : F(F_), M(M_), o(o_), a(a_), b(b_)
               {
                  if (o == nullptr)
                     return;
                  f0a = F[a];
                  m0a = M[a];
                  if (b >= 0)
                  {
                     f0b = F[b];
                     m0b = M[b];
                  }
               }
               ~ImpItem()
               {
                  if (o == nullptr)
                     return;
                  o[0] = F[a] - f0a;
                  o[1] = M[a] - m0a;
                  o[2] = b >= 0 ? F[b] - f0b : V3();
                  o[3] = b >= 0 ? M[b] - m0b : V3();
               }
            } impItem(F, M, impOn ? &impIt[4 * (size_t)k] : nullptr, it.a, it.b);
            const Item &L = law[k];
            ItemEnd &e = ie[k];
            if (J)
               szRec[k].on = false;   // (stick_zone_trial: written below for a sliding patch)
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
            {
               if (clampDump && it.b < 0 && hm.nc > 0)   // VPX_GROUP_CLAMP_TRACE: a carried surface's contact that is off at this Eval
               {
                  char ln[200];
                  snprintf(ln, sizeof ln, "  c k %d ball %d part %d %s off d0 %.5g dm %.5g d1 %.5g um\n", k, it.a, shape[it.s] != nullptr ? shape[it.s]->part : -1, shape[it.s] != nullptr ? shape[it.s]->Name().c_str() : "-", d0 * 1e6, dm * 1e6,
                           d1 * 1e6);
                  clampCon += ln;
               }
               continue;
            }
            e.active = true;
            const bool pr = it.b >= 0;
            const int na = pr ? 6 : 3;
            const int nq = pr ? 0 : hm.nc;   // a carried surface's columns
            const int qc0 = nq > 0 ? coordOff[shape[it.s]->part] : 0;
            const V3 &nm = hm.n;
            // the compression's discrete gradient over the move (ball a's centre, then ball b's, or the columns' coordinates)
            double gm[kGradMax] = { -nm.x * hm.grad, -nm.y * hm.grad, -nm.z * hm.grad, nm.x, nm.y, nm.z }, dq[kGradMax], iM[kGradMax], gd[kGradMax];
            const Ball &A = ball[it.a];
            dq[0] = dc[it.a].x; dq[1] = dc[it.a].y; dq[2] = dc[it.a].z;
            iM[0] = iM[1] = iM[2] = 1.0 / A.m;
            if (pr)
            {
               dq[3] = dc[it.b].x; dq[4] = dc[it.b].y; dq[5] = dc[it.b].z;
               iM[3] = iM[4] = iM[5] = 1.0 / ball[it.b].m;
            }
            for (int j = 0; j < nq; j++)
            {
               const int cj = qc0 + hm.c[j].coord;
               gm[3 + j] = hm.c[j].slope * hm.grad;
               dq[3 + j] = qEnd[cj];
               iM[3 + j] = cUnk[cj] >= 0 ? ci[cj].invM : 0.0;   // (a move given to it takes none of the correction: only the unknowns do)
            }
            bool cf0U = false;
            uint64_t tCf0 = 0;   // (VPX_GROUP_CLAMP_TRACE: the terms whose correction was left out)
            Gonzalez(gm, d1 - d0, noiseD, dq, iM, na + nq, set.gonzalezMass, gd, clampDump ? &cf0U : nullptr);
            // [Contact] hinge_piece_sum: every term of the union by its own chord and its own exact gradient (each term's depth is smooth, so its midpoint gradient
            // misses by the third order only); the elastic force sum_t s_t fe_t gd_t does sum_t s_t dU_t = dU exactly. The contact's normal force N_el and its
            // direction G (the ball's part a unit normal) are that sum's; the damping acts along G on the rate G . dq / h: its loss cD (G . dq)^2 / h is exact for
            // any G, so no misfit reaches the part through it. (A massless part's spring is in its U*: not here.)
            const Hit &h0 = h0s[k];
            const bool terms = hm.nt > 0 && h1.nt == hm.nt && h0.nt == hm.nt && nq == 1;
            std::vector<std::array<double, 4>> &tg = tgS;
            std::vector<double> &tk = tkS;
            std::vector<char> &tOn = tOnS;
            std::vector<std::array<double, 4>> &tgE = tgES;   // each term's end gradient (endGradSlope; else = tg): the column side of its slope
            std::vector<std::array<double, 4>> &tgF = tgFS;   // each term's force direction gdT (termForceRows; else = tg): the row side of its slope
            double dd = d1 - d0, feT = 0.0;
            if (terms)
            {
               if (tk.size() < (size_t)hm.nt)
               {
                  tg.resize((size_t)hm.nt);
                  tgE.resize((size_t)hm.nt);
                  tgF.resize((size_t)hm.nt);
                  tk.resize((size_t)hm.nt);
                  tOn.resize((size_t)hm.nt);
               }
               double Fel[4] = { 0.0, 0.0, 0.0, 0.0 };
               for (int t = 0; t < hm.nt; t++)
               {
                  const double a0 = h0.t[t].d, a1 = h1.t[t].d, am = hm.t[t].d, s = hm.t[t].s;
                  tOn[t] = a0 > 0.0 || a1 > 0.0;
                  if (!tOn[t])
                     continue;
                  const V3 &tn = hm.t[t].n;
                  tg[t][0] = -tn.x; tg[t][1] = -tn.y; tg[t][2] = -tn.z; tg[t][3] = hm.t[t].slope;
                  if (set.endGradSlope)   // d a1/du: the end place's own gradient
                  {
                     const V3 &en = h1.t[t].n;
                     tgE[t][0] = -en.x; tgE[t][1] = -en.y; tgE[t][2] = -en.z; tgE[t][3] = h1.t[t].slope;
                  }
                  else
                     tgE[t] = tg[t];
                  double gdT[kGradMax];
                  bool cf0T = false;
                  Gonzalez(tg[t].data(), a1 - a0, noiseD, dq, iM, 4, set.gonzalezMass, gdT, clampDump ? &cf0T : nullptr);
                  if (cf0T && t < 64)
                     tCf0 |= 1ull << t;
                  for (int j = 0; j < 4; j++)   // the slope's row side: the direction this term's force acts along (termForceRows), or its mid gradient
                     tgF[t][j] = set.termForceRows ? gdT[j] : tg[t][j];
                  const double fT = L.Chord(a0, a1);
                  for (int j = 0; j < 4; j++)
                     Fel[j] += s * fT * gdT[j];
                  tk[t] = s * (std::fabs(a1 - a0) > 1e-13 ? (L.F(a1) - fT) / (a1 - a0) : 0.5 * L.K(std::max(std::max(am, a1), 0.0)));
               }
               const double nEl = std::sqrt(Fel[0] * Fel[0] + Fel[1] * Fel[1] + Fel[2] * Fel[2]);
               if (nEl > 1e-300)
               {
                  feT = nEl;
                  for (int j = 0; j < 4; j++)
                     gd[j] = Fel[j] / nEl;
               }
               else
                  for (int j = 0; j < 4; j++)   // (no elastic force: the union's own mid gradient carries the damping)
                     gd[j] = gm[j];
               dd = gd[0] * dq[0] + gd[1] * dq[1] + gd[2] * dq[2] + gd[3] * dq[3];
            }
            const double fe = terms ? feT : L.Chord(d0, d1);
            // [Contact] damping_integral (di): each damping by its integrated form over the step; cDs its exact slope d(c (d1 - d0))/d d1 (Hunt-Crossley without a bend: its
            // chord's, 1.5 alpha (fe + kc dd), taken in the slope below: hcChord). 0: the committed forms (KK by KKStepMid, the loss tangent and the bend factor at dm)
            const bool di = set.dampingIntegral;
            double cD, kkJac = 0.0, cDs = 0.0;
            bool hcChord = false;
            if (di)
            {
               double sl = 0.0;
               if (L.gammaKK > 0.0)
               {
                  cD = L.KKStep(d0, d1, sl);
                  cDs = sl;
               }
               else if (L.kb > 0.0)
               {
                  cD = L.HCBendStep(d0, d1, sl);
                  cDs = sl;
               }
               else
               {
                  cD = 1.5 * L.alpha * fe;
                  hcChord = true;
               }
               if (L.normalLoss && L.tanDelta > 0.0)
               {
                  cD += L.LossStep(d0, d1, sl);
                  cDs += sl;
               }
            }
            else if (L.gammaKK > 0.0)
            {
               cD = L.KKStepMid(d0, d1, dm, kkJac);
               if (dm > 0.0 && L.normalLoss && L.tanDelta > 0.0)
                  cD += L.tanDelta * std::sqrt(1.5 * L.F(dm) * L.mRed / dm);
            }
            else
               cD = L.DampC(dm, fe);
            // [Contact] damping_integral: the patch of the step's mean elastic force fe (F(dE) = fe), aE, and its slope daE/dfe: the twist limit takes it (no cusp where a
            // contact arrives or leaves within the step: fe and N are step means); [Contact] rolling_patch: rolling resistance takes it too
            const bool rpP = set.rollingPatch;
            double aE = 0.0, dadfE = 0.0;
            if (di || rpP)
               L.PatchOfForce(fe, aE, dadfE);
            const double fd = cD * dd / h;
            const double N = std::max(0.0, fe + fd);   // a contact pushes, it cannot pull (Hertz, no adhesion)
            e.normalLoss = (N - fe) * dd;              // (cD dd^2 / h while it pushes; all the spring's release once damping would make it pull)

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
            // the patch rule (rubbercontact.h ContactCore): a patch that grows holds the force, one that shrinks holds the stretch; the spring's fall is a loss.
            // [Contact] patch_path (pp): the step's spring kS and damper ct are their means over the depth path (u^0.5, u^0.25: Item::PatchMeans; kSp, ctSp their slopes in d1);
            // the stretch goes from k0 to kS by the rule here, and from kS to k1 after the move (below). 0: kS = k1 = Kt(d1), the damper from k1
            const bool pp = set.patchPath;
            Item::Pows pM, pS;
            if (pp)
               L.PatchMeans(d0, d1, pM, pS);
            const double ktU = pp ? L.ktPerA * std::sqrt(L.reff) : 0.0, ctU = pp ? L.shearTan * std::sqrt(ktU * L.mt) : 0.0;   // Kt and its damper per u^0.5, u^0.25
            const double k0 = L.Kt(d0), k1 = L.Kt(d1);
            double kS = pp ? ktU * pM.p05 : k1, kSp = pp ? ktU * pS.p05 : 0.0, ctSp = pp ? ctU * pS.p025 : 0.0;
            // a bend across in series (Item::ktb, a wire switch's bent wire): Kt is no power of u, so its mean over the path is by the quadrature or its integral (Item::MeanOver,
            // KtIntegral) and its slope by MeanSlope; the damper shearTan sqrt(kS mt) on that mean
            const bool bendT = pp && L.ktb > 0.0;
            if (bendT)
            {
               kS = Item::MeanOver(d0, d1, [&L](const double d) { return L.Kt(d); }, [&L](const double d) { return L.KtIntegral(d); });
               kSp = Item::MeanSlope(d0, d1, [&L](const double d) { return L.Kt(d); }, [&L](const double d) { return L.KtSlope(d); }, kS);
               ctSp = kS > 0.0 ? 0.5 * L.shearTan * std::sqrt(L.mt / kS) * kSp : 0.0;
            }
            const double held = 0.5 * k0 * xi0.Dot(xi0);
            if (kS > k0)
               xi0 = xi0 * (k0 / kS);
            e.slip = held - 0.5 * kS * xi0.Dot(xi0);
            // the patch: the spring kS and the damper ct side by side, in series with the Coulomb slip. The trial force (no slip) is the spring's midpoint force and the
            // damper's; past mu N the patch slips along the trial and the spring's end stretch is the one that carries the limit force with the damper. The force is
            // continuous across stick and slip.
            V3 Ft, xi1, FsTr;
            double ct = 0.0, limS = 0.0, flTr = 0.0;
            bool slipping = false;
            if (kS > 0.0)
            {
               const V3 xtr = xi0 + ut * h;
               ct = bendT ? L.shearTan * std::sqrt(kS * L.mt) : pp ? ctU * pM.p025 : L.shearTan * std::sqrt(k1 * L.mt);
               Ft = (xi0 + xtr) * (0.5 * kS) + ut * ct;
               FsTr = Ft;
               const double lim = L.mu * N, fl = Ft.Length();
               limS = lim;
               flTr = fl;
               xi1 = xtr;
               if (fl > lim)
               {
                  Ft = fl > 0.0 ? Ft * (lim / fl) : V3();
                  // Ft = kS (xi0 + xi1) / 2 + ct (xi1 - xi0) / h, so the slip xtr - xi1 = (trial - Ft) / (kS / 2 + ct / h) lies along Ft: its loss Ft . slip >= 0
                  const double cs = 0.5 * kS + ct / h;
                  xi1 = (Ft - xi0 * (0.5 * kS - ct / h)) * (1.0 / cs);
                  slipping = true;
                  e.slipping = true;
                  e.slip += Ft.Dot(xtr - xi1);
               }
            }
            {
               const V3 dxi = xi1 - xi0;   // the damper's loss over its own stretch (ct ut^2 h while it sticks)
               e.shearDamp = slipping ? ct * dxi.Dot(dxi) / h : ct * ut.Dot(ut) * h;
            }
            if (pp && kS > 0.0)   // patch_path: the stretch from kS to the end's patch k1 by the patch rule (a loss), stored on k1: the next step starts on it
            {
               const double heldS = 0.5 * kS * xi1.Dot(xi1);
               if (k1 > kS)
                  xi1 = xi1 * (kS / k1);
               e.slip += heldS - 0.5 * k1 * xi1.Dot(xi1);
            }
            e.xi1 = xi1;
            F[it.a] = F[it.a] - V3(gd[0], gd[1], gd[2]) * N - Ft;   // (a massless part's contact: its spring is in the part's U*, only the damping here)
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
            if (ContactTraceNow() || HistoryTraceOn())   // (diagnostics: copies only)
            {
               if ((int)cDiag.size() < ni)
                  cDiag.resize(ni);
               ContactDiag &cd = cDiag[k];
               cd.n = nm;
               cd.dm = dm;
               cd.N = N;
               cd.fe = fe;
               cd.fd = fd;
               cd.ft = Ft.Length();
               cd.ftTr = flTr;
               cd.lim = limS;
               cd.ut = ut.Length();
               cd.nq = nq;
               cd.qN = nq > 0 ? -gd[3] * N : 0.0;
               cd.qF = nq > 0 ? Ft.Dot(hm.c[0].V()) : 0.0;
            }
            // twist and rolling resistance on the mid-step spin (relative to the other ball, or the surface)
            V3 wr = pr ? wm[it.a] - wm[it.b] : wm[it.a];
            for (int j = 0; j < nq; j++)
               wr = wr - hm.c[j].W() * (dq[3 + j] / h);
            const double am = L.A(dm);
            double kTw = 0.0, sRo = 0.0, rhoR = 0.0, wlR = 0.0, kkTw = 0.0;
            double limN = 0.0;   // rolling at the limit (roll_model 1): d lim / dN = muRr R, the torque along mhR (leverRollSlope)
            double limF = 0.0;   // rolling_patch: rolling at the limit, d lim / dfe through the patch aE (muRr between its floor and cap), on h dfe/du; limN then R (muRr + N dmuRr/dN)
            const double aR = rpP ? aE : am;   // rolling resistance's patch: rolling_patch, the step's mean elastic force's aE; 0, the mid depth's a(dm)
            // twist (roll_model 1; twistSlope): sticking dMw/dd1 (twD1); at the limit dMw/d dm (twA) and dMw/dN (twN). twSt 0 none, 1 stick, 2 limit; twM, twLim the torque and
            // its limit (diagnostics)
            double twD1 = 0.0, twA = 0.0, twN = 0.0, twM = 0.0, twLim = 0.0;
            double twF = 0.0;   // damping_integral: dMw/dfe at the limit (the limit's patch aE), on h dfe/du (h kc gm; a terms contact: kvE)
            int twSt = 0;
            V3 wtR, mhR;
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
               // patch_path: the step's springs krS, kwS and dampers are their means over the depth path (u^1.5, u^0.75), the patch rule from the start's to them here and
               // from them to the end's after each move; 0: krS = kr1, kwS = kw1, the dampers from them
               const double kr0 = L.Kr(d0), kr1 = L.Kr(d1), kw0 = L.Ktw(d0), kw1 = L.Ktw(d1);
               const double Ir = L.mRed * A.R * A.R;   // with 1.4: the ball's inertia about its contact point; with 0.4: about its centre (the contact's mass)
               const double krU = pp ? L.kH * L.reff : 0.0, kwU = 16.0 / 3.0 * 0.35 * 0.75 * krU;   // Kr, Ktw per u^1.5 (a^3 = reff^1.5 u^1.5)
               const double krS = pp ? krU * pM.p15 : kr1, kwS = pp ? kwU * pM.p15 : kw1;
               const double held = 0.5 * kr0 * r0.th.Dot(r0.th) + 0.5 * kw0 * r0.tw * r0.tw;
               if (krS > kr0)
                  r0.th = r0.th * (kr0 / krS);
               if (kwS > kw0)
                  r0.tw *= kw0 / kwS;
               e.roll = held - 0.5 * krS * r0.th.Dot(r0.th) - 0.5 * kwS * r0.tw * r0.tw;   // the patch's change: a loss >= 0
               if (di ? L.twist && kwS > 0.0 && N > 0.0 : L.twist && kwS > 0.0 && N > 0.0 && am > 0.0)   // (damping_integral: no switch on dm; with no patch on the step's path, kwS = 0, there is no twist spring or damper)
               {
                  const double wn = wr.Dot(nm);
                  const double cwU = pp ? L.shearTan * std::sqrt(kwU * 0.4 * Ir) : 0.0;   // patch_path: the damper per u^0.75
                  const double cw = pp ? cwU * pM.p075 : L.shearTan * std::sqrt(kw1 * 0.4 * Ir);
                  const double twTr = r0.tw + wn * h;
                  double Mw = 0.5 * kwS * (r0.tw + twTr) + cw * wn;
                  const double limW = 3.0 * kPi / 16.0 * L.mu * (di ? aE : am) * N;
                  double tw1 = twTr;
                  kkTw = 0.5 * kwS * h + cw;   // the resisting torque's slope in wn while it sticks; 0 at the limit
                  if (std::fabs(Mw) > limW)
                  {
                     Mw = Mw > 0.0 ? limW : -limW;
                     const double cs = 0.5 * kwS + cw / h;
                     tw1 = (Mw - r0.tw * (0.5 * kwS - cw / h)) / cs;
                     e.twist = Mw * (twTr - tw1) + cw * (tw1 - r0.tw) * (tw1 - r0.tw) / h;
                     kkTw = 0.0;
                     twSt = 2;
                     if (set.twistSlope)   // Mw = s 3 pi/16 mu a(dm) N: its change with dm and with N (damping_integral: a(fe) N, with fe and N)
                     {
                        const double cl = (Mw > 0.0 ? 1.0 : -1.0) * 3.0 * kPi / 16.0 * L.mu;
                        if (di)
                        {
                           twF = cl * dadfE * N;
                           twN = cl * aE;
                        }
                        else
                        {
                           twA = cl * L.ASlope(dm) * N;
                           twN = cl * am;
                        }
                     }
                  }
                  else
                  {
                     e.twist = cw * wn * wn * h;
                     twSt = 1;
                     if (set.twistSlope)   // Mw = kwS tw0' + (kwS h / 2 + cw) wn: kwS = Ktw(d1) ~ a^3 (patch_path: its mean), cw ~ sqrt(kwS) (its mean); a growing patch holds kwS tw0' = kw0 tw0
                     {
                        if (pp)   // the means' exact slopes in d1
                           twD1 = kwU * pS.p15 * ((kwS > kw0 ? 0.0 : r0.tw) + 0.5 * h * wn) + cwU * pS.p075 * wn;
                        else
                        {
                           const double a1 = L.A(d1);
                           const double kw1p = a1 > 0.0 ? 3.0 * kw1 * L.ASlope(d1) / a1 : 0.0;
                           const double cwp = 0.5 * cw * kw1p / kw1;
                           twD1 = kw1p * ((kw1 > kw0 ? 0.0 : r0.tw) + 0.5 * h * wn) + cwp * wn;
                        }
                     }
                  }
                  if (pp)   // patch_path: the twist from kwS to the end's patch kw1 by the patch rule (a loss), stored on kw1
                  {
                     const double heldS = 0.5 * kwS * tw1 * tw1;
                     if (kw1 > kwS)
                        tw1 *= kwS / kw1;
                     e.roll += heldS - 0.5 * kw1 * tw1 * tw1;
                  }
                  twM = Mw;
                  twLim = limW;
                  e.ra1.tw = tw1;
                  kTw = kkTw > 0.0 ? 1.0 : 0.0;   // (the Jacobian's gate)
                  M[it.a] = M[it.a] - nm * Mw;
                  if (pr)
                     M[it.b] = M[it.b] + nm * Mw;
                  for (int j = 0; j < nq; j++)
                     qf[j] += Mw * nm.Dot(hm.c[j].W());
               }
               if ((L.tanDelta > 0.0 || L.rollingMu > 0.0) && N > 0.0 && aR > 0.0 && krS > 0.0)   // (rolling_patch: no switch on dm, the patch aE of the step's fe)
               {
                  const V3 wt = wr - nm * wr.Dot(nm);
                  const double pa = N / (kPi * aR * aR);
                  const double muP = 2.34 * (pa / L.E) * L.tanDelta;   // the law's own value, before its cap and floor
                  const double muRr = std::max(std::min(muP, RubberContact::kRollingMuMax), L.rollingMu);
                  const double lim = muRr * A.R * N;
                  const double cr = pp ? L.shearTan * std::sqrt(krU * 1.4 * Ir) * pM.p075 : L.shearTan * std::sqrt(kr1 * 1.4 * Ir);   // (patch_path: the damper's mean, u^0.75)
                  const V3 thTr = r0.th + wt * h;
                  V3 Mr = (r0.th + thTr) * (0.5 * krS) + wt * cr;
                  V3 th1 = thTr;
                  const double ml = Mr.Length(), sl = 0.5 * krS * h + cr;
                  const M3 Pt = M3::Tangent(nm);
                  if (ml > lim)
                  {
                     // rolls: the torque at the limit along the trial; the roll thTr - th1 lies along it, its loss Mr . roll >= 0
                     Mr = Mr * (lim / ml);
                     const double cs = 0.5 * krS + cr / h;
                     th1 = (Mr - r0.th * (0.5 * krS - cr / h)) * (1.0 / cs);
                     const V3 dth = th1 - r0.th;
                     e.roll += Mr.Dot(thTr - th1) + cr * dth.Dot(dth) / h;
                     const V3 mh = Mr * (1.0 / lim);
                     const M3 Q = M3::Outer(mh, mh);
                     for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                           RmS.a[i][j] = lim / ml * sl * (Pt.a[i][j] - Q.a[i][j]);
                     limN = muRr * A.R;
                     // rolling_patch: muRr = 2.34 tanDelta N / (pi aE^2 E) between floor and cap: dmuRr/dN = muRr / N, dmuRr/daE = -2 muRr / aE (0 at the cap or the floor)
                     if (rpP && muP < RubberContact::kRollingMuMax && muP > L.rollingMu)
                     {
                        limN = 2.0 * muRr * A.R;
                        limF = -2.0 * muRr * A.R * N / aR * dadfE;
                     }
                     mhR = mh;
                  }
                  else
                  {
                     e.roll += cr * wt.Dot(wt) * h;
                     for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                           RmS.a[i][j] = sl * Pt.a[i][j];
                  }
                  if (pp)   // patch_path: the tilt from krS to the end's patch kr1 by the patch rule (a loss), stored on kr1
                  {
                     const double heldS = 0.5 * krS * th1.Dot(th1);
                     if (kr1 > krS)
                        th1 = th1 * (krS / kr1);
                     e.roll += heldS - 0.5 * kr1 * th1.Dot(th1);
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
            if (set.rollModel != 1 && (L.tanDelta > 0.0 || L.rollingMu > 0.0) && N > 0.0 && aR > 0.0)   // (rolling_patch: the same patch aE; no slope of rho here, as with 0)
            {
               const V3 wt = wr - nm * wr.Dot(nm);
               const double wl = wt.Length();
               const double pa = N / (kPi * aR * aR);
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

            if (clampDump && nq > 0)   // VPX_GROUP_CLAMP_TRACE: this carried surface's contact at this Eval (read only): the residual's switches
            {
               const char kk = L.gammaKK > 0.0 ? 'K' : 'H';
               std::string on, cf;
               for (int t = 0; terms && t < hm.nt; t++)
               {
                  on += tOn[t] ? '1' : '0';
                  cf += t < 64 && (tCf0 >> t & 1u) ? '1' : '0';
               }
               // meff: the mass along the contact's normal, 1 / sum_j iM_j gd_j^2 (the inverse masses and Newton's gd), g (-1: no unknown takes the correction); ks: the law's
               // dF/dd at d1, N/m; cD: this round's damper, N s/m; vn0: the closing rate at the step's start, the mid gradient gm on the start's speeds and rates, m/s
               double wSum = 0.0, vn0 = 0.0;
               for (int j = 0; j < na + nq; j++)
                  wSum += iM[j] * gd[j] * gd[j];
               vn0 += gm[0] * A.v.x + gm[1] * A.v.y + gm[2] * A.v.z;
               if (pr)
                  vn0 += gm[3] * ball[it.b].v.x + gm[4] * ball[it.b].v.y + gm[5] * ball[it.b].v.z;
               for (int j = 0; j < nq; j++)
                  vn0 += gm[3 + j] * ci[qc0 + hm.c[j].coord].rate;
               char ln[600];
               snprintf(ln, sizeof ln, "  c k %d ball %d part %d %s on d0 %.5g dm %.5g d1 %.5g um fe %.6g fd %.6g N %.6g clamp %d slide %d kk %c cf0 %d terms %s tcf0 %s meff %.6g g ks %.6g N/m cD %.6g Ns/m vn0 %.6g m/s tw %c Mw %.4g limW %.4g N m\n", k, it.a,
                        shape[it.s] != nullptr ? shape[it.s]->part : -1, shape[it.s] != nullptr ? shape[it.s]->Name().c_str() : "-", d0 * 1e6, dm * 1e6, d1 * 1e6, fe, fd, N, fe + fd <= 0.0 ? 1 : 0, slipping ? 1 : 0, kk, cf0U ? 1 : 0,
                        terms ? on.c_str() : "-", terms ? cf.c_str() : "-", wSum > 0.0 ? 1e3 / wSum : -1.0, L.K(d1), cD, vn0, "-sl"[twSt], twM, twLim);
               clampCon += ln;
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
               const double dSl = cDs + (hcChord ? 1.5 * L.alpha * (fe + kc * dd) : 0.0);   // damping_integral: the damping's exact slope (all its laws)
               const double kn = fe + fd <= 0.0 ? 0.0 : di ? h * kc + dSl : L.gammaKK > 0.0 ? h * kc + kkJac : h * kc + cD + cDp * (d1 - d0);
               int ix[kGradMax];   // each gradient entry's unknown; -1 = a given move
               for (int p = 0; p < na; p++)
                  ix[p] = (p < 3 ? 6 * it.a : 6 * it.b) + p % 3;
               for (int j = 0; j < nq; j++)
                  ix[3 + j] = cUnk[qc0 + hm.c[j].coord];
               // dampDirSlope: G's turning, h dG/du = (I - G g^T) K / fe (K = sum_t h tk_t gdT_t tgE_t^T, termForceRows; else tg_t as the row; g = G's ball part), in Newton's slope: the damping's fd along it and
               // dd = G . dq's change through it, h dfd/du += (cD / (h fe)) (dq - dd g)^T K (gW, on kv too). The residual is unchanged.
               const bool dirOn = terms && set.dampDirSlope && feT > 0.0 && fe + fd > 0.0;
               double gK[4][4] = {}, gW[4] = {};
               if (dirOn)
               {
                  double K[4][4] = {};
                  for (int t = 0; t < hm.nt; t++)
                     if (tOn[t])
                        for (int p = 0; p < 4; p++)
                           for (int q = 0; q < 4; q++)
                              K[p][q] += h * tk[t] * tgF[t][p] * tgE[t][q];   // (row p: the force's direction, termForceRows; column q: d a1/du, endGradSlope)
                  for (int q = 0; q < 4; q++)
                  {
                     double gKq = 0.0, w = 0.0;
                     for (int r = 0; r < 4; r++)
                     {
                        if (r < 3)
                           gKq += gd[r] * K[r][q];
                        w += (dq[r] - (r < 3 ? dd * gd[r] : 0.0)) * K[r][q];
                     }
                     gW[q] = cD / (h * feT) * w;   // (fd = cD dd / h: dd's change through G takes the coefficient; cD's own depth slope is on G g^T)
                     for (int p = 0; p < 4; p++)
                        gK[p][q] = fd / feT * (K[p][q] - gd[p] * gKq);
                  }
                  for (int p = 0; p < 4; p++)
                     for (int q = 0; q < 4; q++)
                        gK[p][q] += gd[p] * gW[q];
               }
               if (terms)   // each term's own chord slope on its own gradient, the damping's on G
               {
                  if (fe + fd > 0.0)
                  {
                     const double kd = di ? dSl : L.gammaKK > 0.0 ? kkJac : cD + cDp * dd;
                     for (int p = 0; p < 4; p++)
                        if (ix[p] >= 0)
                           for (int q = 0; q < 4; q++)
                              if (ix[q] >= 0)
                              {
                                 double a = kd * gd[p] * gd[q];
                                 for (int t = 0; t < hm.nt; t++)
                                    if (tOn[t])
                                       a += h * tk[t] * tgF[t][p] * tgE[t][q];   // (row: the force's direction, termForceRows; column: d a1/du)
                                 if (dirOn)
                                    a += gK[p][q];
                                 J[ix[p]][ix[q]] += a;
                              }
                  }
               }
               else
                  for (int p = 0; p < na + nq; p++)
                     if (ix[p] >= 0)
                        for (int q = 0; q < na + nq; q++)
                           if (ix[q] >= 0)
                              J[ix[p]][ix[q]] += kn * gm[p] * gm[q];
               if (jacWinDump)   // VPX_GROUP_JAC_WINDOW: this contact's state at this Eval (read only)
               {
                  const int sp = !pr && shape[it.s] != nullptr ? shape[it.s]->part : -1;
                  char with[40], ln[780];
                  if (pr)
                     snprintf(with, sizeof with, "ball %d", it.b);
                  else if (sp >= 0)
                     snprintf(with, sizeof with, "part %d (shape %d)", sp, it.s);
                  else
                     snprintf(with, sizeof with, "surface (shape %d)", it.s);
                  snprintf(ln, sizeof ln, "  contact k %d ball %d with %s: d0 %.5g dm %.5g d1 %.5g um n (%.4f %.4f %.4f) N %.5g fe %.5g fd %.5g cD %.4g kc %.5g kn %.5g | %s ut (%.3g %.3g %.3g) m/s xi0 (%.3g %.3g %.3g) xi1 (%.3g %.3g %.3g) um Ft (%.4g %.4g %.4g) muN %.4g kS %.4g ct %.4g terms %d tw %c Mw %.4g limW %.4g\n",
                           k, it.a, with, d0 * 1e6, dm * 1e6, d1 * 1e6, nm.x, nm.y, nm.z, N, fe, fd, cD, kc, kn, slipping ? "slide" : "stick", ut.x, ut.y, ut.z, xi0.x * 1e6, xi0.y * 1e6, xi0.z * 1e6,
                           xi1.x * 1e6, xi1.y * 1e6, xi1.z * 1e6, Ft.x, Ft.y, Ft.z, limS, kS, ct, terms ? hm.nt : 0, "-sl"[twSt], twM, twLim);
                  jacWinCon += ln;
                  if (terms)
                     for (int t = 0; t < hm.nt; t++)
                     {
                        snprintf(ln, sizeof ln, "    term %d s %g d0 %.5g dm %.5g d1 %.5g um n (%.4f %.4f %.4f) slope %.5g on %d tk %.5g\n", t, hm.t[t].s, h0.t[t].d * 1e6, hm.t[t].d * 1e6, h1.t[t].d * 1e6,
                                 hm.t[t].n.x, hm.t[t].n.y, hm.t[t].n.z, hm.t[t].slope, tOn[t] ? 1 : 0, tOn[t] ? tk[t] : 0.0);
                        jacWinCon += ln;
                     }
               }
               // shear: the tangential force Ft changes with the slip rate ut (dFt/dut = Cm, in the contact plane) and with the depth d1 (dFt/dd1 = fd1).
               // ut = sp (dc_p / h + 1/2 Skew(r_p)^T dw_p) projected, d1 changes by -sq n . dc_q; on ball p the force is -sp Ft at r_p. So, for the residual of p:
               //    d/d dc_q: sp sq [Cm; Skew(r_p) Cm] - h sp sq [fd1; Skew(r_p) fd1] n^T,   d/d dw_q: 1/2 h sp sq [Cm; Skew(r_p) Cm] Skew(r_q)^T
               // The spring's force is kS xi0 + kS ut h / 2 (kS = k1, or patch_path: its mean over the path), and a growing patch holds kS xi0 (= k0 xi0 before the hold): g = kS'(d1) ((kS > k0 ? 0 : xi0) + ut h / 2) is its change with d1.
               // Sticking: Cm = (kS h / 2 + ct) Pt, fd1 = g + ct'(d1) ut (ct grows with the patch too: shearTan sqrt(mt / k1) k1' / 2, or patch_path its mean's slope).
               // Sliding (Ft = lim f^, f^ along the total trial force T = kS xi0 + (kS h / 2 + ct) ut): Cm = (lim / |T|) (kS h / 2 + ct) (Pt - f^ f^^T),
               // fd1 = mu N'(d1) f^ + (lim / |T|) (Pt - f^ f^^T) (g + gc), N' = kn / h, the normal row's own slope.
               const M3 Pt = M3::Tangent(nm);
               const int nbod = pr ? 2 : 1;
               const int bi[2] = { it.a, it.b };
               const V3 rr[2] = { rA, rB };
               const double sg[2] = { 1.0, -1.0 };
               // hinge_piece_sum: the normal force's own slope over the unknowns (ball centre, column), h dN/du: each pressing term's chord slope on its gradient,
               // projected on the force's direction, and the damping's on G. (The summed depth's slope kn along one normal is not N's: two terms at an angle.)
               double kv[4] = { 0.0, 0.0, 0.0, 0.0 };
               if (terms && fe + fd > 0.0)
               {
                  const double kd = di ? dSl : L.gammaKK > 0.0 ? kkJac : cD + cDp * dd;
                  for (int t = 0; t < hm.nt; t++)
                     if (tOn[t])
                     {
                        const double w = h * tk[t] * (gd[0] * tgF[t][0] + gd[1] * tgF[t][1] + gd[2] * tgF[t][2]);
                        for (int u = 0; u < 4; u++)
                           kv[u] += w * tgE[t][u];   // (d a1/du, endGradSlope)
                     }
                  for (int u = 0; u < 4; u++)
                     kv[u] += dirOn ? kd * gd[u] + gW[u] : kd * gd[u];
               }
               // damping_integral, a terms contact: the elastic force's own slope h dfe/du (each pressing term's chord slope projected on G; kv without the damping), the
               // column the twist limit's patch aE takes (and rolling_patch's rolling limit). Not a terms contact: h dfe/du = h kc gm
               double kvE[4] = { 0.0, 0.0, 0.0, 0.0 };
               if ((di || rpP) && terms)
                  for (int t = 0; t < hm.nt; t++)
                     if (tOn[t])
                     {
                        const double w = h * tk[t] * (gd[0] * tgF[t][0] + gd[1] * tgF[t][1] + gd[2] * tgF[t][2]);
                        for (int u = 0; u < 4; u++)
                           kvE[u] += w * tgE[t][u];
                     }
               if (kS > 0.0)
               {
                  M3 Cm;
                  V3 fd1, fN;   // fN: sliding, the limit's change per unit of N (hinge_piece_sum: taken on kv, not in fd1)
                  // the spring's and damper's change with d1: patch_path, the means' exact slopes (kSp, ctSp); 0: at d1, KtSlope and ct's through it
                  const double kp = pp ? kSp : L.KtSlope(d1);
                  const V3 g = ((kS > k0 ? V3() : xi0) + ut * (0.5 * h)) * kp;
                  const V3 gc = pp ? ut * ctSp : ut * (0.5 * L.shearTan * std::sqrt(L.mt / k1) * kp);
                  const double C = 0.5 * kS * h + ct;
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
                     if (terms)
                     {
                        fN = fh * L.mu;
                        fd1 = Pg * r;
                     }
                     else
                        fd1 = fh * (L.mu * dN) + Pg * r;
                  }
                  const M3 Dn = M3::Outer(fd1, nm);
                  if (set.stickZoneTrial && slipping && limS > 0.0 && flTr > 0.0)   // stick_zone_trial: the trial force's linear model along a move (kS, ct, stretch, geometry held)
                  {
                     StickRec &sr = szRec[k];
                     sr.on = true;
                     sr.a = FsTr;
                     sr.C = C;
                     sr.mu = L.mu;
                     sr.N0 = N;
                     sr.A = it.a;
                     sr.B = pr ? it.b : -1;
                     sr.rA = rA;
                     sr.rB = rB;
                     sr.nm = nm;
                     sr.nq = nq;
                     for (int j = 0; j < nq; j++)
                     {
                        sr.cV[j] = hm.c[j].V();
                        sr.qU[j] = ix[3 + j];
                     }
                     sr.nn = 0;   // N's slope over the unknowns: the normal row's own (h dN/du = kv, or kn gm), over h
                     if (terms)
                     {
                        for (int u = 0; u < 4; u++)
                           if (ix[u] >= 0 && kv[u] != 0.0)
                           {
                              sr.nG[sr.nn] = kv[u] / h;
                              sr.nIx[sr.nn++] = ix[u];
                           }
                     }
                     else if (kn != 0.0)
                        for (int p = 0; p < na + nq; p++)
                           if (ix[p] >= 0)
                           {
                              sr.nG[sr.nn] = kn / h * gm[p];
                              sr.nIx[sr.nn++] = ix[p];
                           }
                  }
                  if (terms)   // (only a ball on a part: nbod 1) the limit's change with N: mu f^ (h dN/du)
                  {
                     const M3 Nn = M3::Outer(fN, V3(kv[0], kv[1], kv[2]));
                     AddBlock(J, 6 * it.a, 6 * it.a, Nn, 1.0);
                     AddBlock(J, 6 * it.a + 3, 6 * it.a, M3::Skew(rA) * Nn, 1.0);
                  }
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
                        const V3 rc = cTv * -1.0 + nm * (h * vi.Dot(fd1)) - (terms ? V3(kv[0], kv[1], kv[2]) * vi.Dot(fN) : V3());
                        const V3 rw = MulV(M3::Skew(rA), cTv) * (-0.5 * h);
                        const int ia = 6 * it.a;
                        J[qi][ia] += rc.x; J[qi][ia + 1] += rc.y; J[qi][ia + 2] += rc.z;
                        J[qi][ia + 3] += rw.x; J[qi][ia + 4] += rw.y; J[qi][ia + 5] += rw.z;
                     }
                     const double li = hm.c[i].slope * hm.grad;
                     if (qi >= 0)   // the ball's rows against it: h dFt / dq = -Cm v + h slope fd1
                     {
                        const V3 dFq = MulV(Cm, vi) * -1.0 + fd1 * (h * li) + (terms ? fN * kv[3] : V3());
                        const V3 dMq = MulV(M3::Skew(rA), dFq);
                        const int ia = 6 * it.a;
                        J[ia][qi] += dFq.x; J[ia + 1][qi] += dFq.y; J[ia + 2][qi] += dFq.z;
                        J[ia + 3][qi] += dMq.x; J[ia + 4][qi] += dMq.y; J[ia + 5][qi] += dMq.z;
                     }
                     for (int j = 0; j < nq && qi >= 0; j++)
                     {
                        const int qj = cUnk[qc0 + hm.c[j].coord];
                        if (qj >= 0)
                           J[qi][qj] += vi.Dot(MulV(Cm, hm.c[j].V())) - h * vi.Dot(fd1) * hm.c[j].slope * hm.grad - (terms ? vi.Dot(fN) * kv[3] : 0.0);
                     }
                  }
                  // leverRollSlope: the contact point moves with the mid depth, d rr_p / du_j = s_p (gm_j / 4) n (dm moves half the depth's change), so the slip takes
                  // (gm_j / 4) wS, wS = sum_p wm_p x n. On ball p's force row (h s_p Ft) and spin row (h s_p rr_p x Ft), every unknown the depth moves:
                  //    h gm_j s_p Cm wS / 4,   h gm_j (n x Ft / 4 + s_p rr_p x Cm wS / 4)
                  if (set.leverRollSlope && dm > 0.0)
                  {
                     V3 wS;
                     for (int p = 0; p < nbod; p++)
                        wS = wS + wm[bi[p]].Cross(nm);
                     const V3 fw = MulV(Cm, wS) * 0.25, nF = nm.Cross(Ft) * 0.25;
                     for (int p = 0; p < nbod; p++)
                     {
                        const V3 fr = fw * sg[p], mr = nF + rr[p].Cross(fw) * sg[p];
                        const int ib = 6 * bi[p];
                        for (int j = 0; j < na + nq; j++)
                           if (ix[j] >= 0)
                           {
                              const double c = h * gm[j];
                              J[ib][ix[j]] += c * fr.x; J[ib + 1][ix[j]] += c * fr.y; J[ib + 2][ix[j]] += c * fr.z;
                              J[ib + 3][ix[j]] += c * mr.x; J[ib + 4][ix[j]] += c * mr.y; J[ib + 5][ix[j]] += c * mr.z;
                           }
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
               // twistSlope: the twist torque's change with the depth, h dMw/du_j = h gm_j (twD1 + twA / 2) + twN h dN/du_j (d d1/du = gm, d dm/du = gm / 2; h dN/du = kn gm,
               // a terms contact: kv). Ball p's spin row holds h s_p n Mw, a carried part's row -h (n . W_i) Mw
               if (set.twistSlope && (twD1 != 0.0 || twA != 0.0 || twN != 0.0 || twF != 0.0))
                  for (int j = 0; j < na + nq; j++)
                  {
                     if (ix[j] < 0)
                        continue;
                     double c = h * gm[j] * (twD1 + 0.5 * twA) + twN * (terms ? kv[j] : kn * gm[j]);
                     if (di)   // damping_integral: dMw/dfe on h dfe/du (h kc gm; a terms contact: kvE)
                        c += twF * (terms ? kvE[j] : h * kc * gm[j]);
                     for (int p = 0; p < nbod; p++)
                     {
                        const int iw = 6 * bi[p] + 3;
                        const double cs = c * sg[p];
                        J[iw][ix[j]] += cs * nm.x; J[iw + 1][ix[j]] += cs * nm.y; J[iw + 2][ix[j]] += cs * nm.z;
                     }
                     for (int i = 0; i < nq; i++)
                     {
                        const int qi = cUnk[qc0 + hm.c[i].coord];
                        if (qi >= 0)
                           J[qi][ix[j]] += -c * nm.Dot(hm.c[i].W());
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
               // leverRollSlope: rolling at the limit, ball p's spin row h s_p mhR lim changes with N by s_p mhR muRr R (h dN/du_j); h dN/du is the normal row's
               // own slope (kn gm; a terms contact: kv). rolling_patch: limN = R (muRr + N dmuRr/dN), and limF = d lim/dfe through the patch aE on h dfe/du (h kc gm; a terms
               // contact: kvE)
               if (set.leverRollSlope && limN > 0.0)
                  for (int p = 0; p < nbod; p++)
                  {
                     const V3 mr = mhR * (limN * sg[p]), mf = mhR * (limF * sg[p]);
                     const int iw = 6 * bi[p] + 3;
                     for (int j = 0; j < na + nq; j++)
                        if (ix[j] >= 0)
                        {
                           const double c = terms ? kv[j] : kn * gm[j];
                           J[iw][ix[j]] += c * mr.x; J[iw + 1][ix[j]] += c * mr.y; J[iw + 2][ix[j]] += c * mr.z;
                           if (limF != 0.0)
                           {
                              const double cF = terms ? kvE[j] : h * kc * gm[j];
                              J[iw][ix[j]] += cF * mf.x; J[iw + 1][ix[j]] += cF * mf.y; J[iw + 2][ix[j]] += cF * mf.z;
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
         {
            if (impOn)
               for (int i = 0; i < nb; i++)
               {
                  impPF[(size_t)p * nb + i] = F[i];
                  impPM[(size_t)p * nb + i] = M[i];
               }
            if (!part[p]->Eval(xv, partOff[p], h, F, M, Res, J))
            {
               partEvalFailed = true;
               return false;
            }
            if (impOn)
               for (int i = 0; i < nb; i++)
               {
                  impPF[(size_t)p * nb + i] = F[i] - impPF[(size_t)p * nb + i];
                  impPM[(size_t)p * nb + i] = M[i] - impPM[(size_t)p * nb + i];
               }
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
            const V3 rc = InertiaRow(dc[i], b.v, b.m, h) - F[i] * h;
            const V3 rw = dw[i] * b.I() - M[i] * h;
            Res[6 * i] = rc.x; Res[6 * i + 1] = rc.y; Res[6 * i + 2] = rc.z;
            Res[6 * i + 3] = rw.x; Res[6 * i + 4] = rw.y; Res[6 * i + 5] = rw.z;
            if (J)
               for (int a = 0; a < 3; a++)
               {
                  J[6 * i + a][6 * i + a] += InertiaSlope(b.m, h);
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
      double EnergyStepCore(const double hTry)
      {
         struct InEnergy { bool &f; explicit InEnergy(bool &b) : f(b) { f = true; } ~InEnergy() { f = false; } } inEnergy(histInEnergy);   // (diagnostics: the history's label)
         const NearOff nearOff(nearOn);
         nearGrow = INFINITY;   // (its explicit moves are not checked against a list: every piece)
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
         PhaseScope ps(phBegin);
         const int nU = nUnk = NumUnk();
         BuildNear(h);
         for (int k = 0; k < ni; k++)
         {
            Hit h0;
            if (nowStamp == stamp && item[k].nearK < 0)   // (a listed contact's start takes its list's terms: probed again)
               h0 = hNow[k];
            else
               Probe(item[k], ball, h0);
            d0s[k] = h0.delta;
            if (h0.nt > 0)
               h0s[k] = h0;
            else
               h0s[k].nt = 0;
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
         {
            part[p]->guessAccel = set.partGuessAccel;
            part[p]->gripInStep = set.gripInStep;
            part[p]->WipeAttemptWork();   // (an attempt's own values start afresh: PartRecords)
            if (!part[p]->BeginStep(ball, nb, h, x, partOff[p]))
            {
               failWhy = 'P';
               return false;
            }
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
         nearGrow = 1.0;
         for (;;)
         {
            farOff = false;
            if (!SolveOnce(h))
               return false;
            if (!FarStillFar(x))
            {
               farOff = true;
               const bool ok = SolveOnce(h);
               farOff = false;
               return ok;   // (its lists hold every piece)
            }
            if (NearStillNear(x))
               return true;
            nearGrow *= 4.0;   // a piece off a list was reached: the lists anew from a longer reach, the attempt again
         }
      }
      // diagnostics, off by default (VPX_GROUP_ROUND_TRACE=1): each round of an attempt, the residual by row in speed units and the parts' row terms; the first failed attempt of a step kept in failRoundBuf
      static bool RoundTraceOn() { static const bool on = []() { const char *const e = getenv("VPX_GROUP_ROUND_TRACE"); return e != nullptr && atoi(e) != 0; }(); return on; }
      bool searchTrace = false;
      std::string searchBuf;
      static bool SearchTraceOn() { static const bool on = []() { const char *const e = getenv("VPX_GROUP_SEARCH_TRACE"); return e != nullptr && atoi(e) != 0; }(); return on; }
      bool SolveOnce(double h)
      {
         if (ClampFile() == nullptr)
            return SolveOnceCore(h);
         clampAttempt++;
         clampRetake = false;
         const bool ok = SolveOnceCore(h);
         if (!ok)
            clampFails++;
         char ln[120];
         snprintf(ln, sizeof ln, "A att %d h %.4g us %s\n", clampAttempt, h * 1e6, ok ? "converged" : "failed");
         clampStep += ln;
         return ok;
      }
      // VPX_GROUP_CLAMP_TRACE: the step's rounds to the file if it took more than 10 rounds or an attempt failed, and it had a carried surface's contact
      void ClampFlush()
      {
         if (clampStep.find("  c k ") != std::string::npos && (clampRounds > 10 || clampFails > 0))
         {
            fprintf(ClampFile(), "S t %u step %ld attempts %d rounds %d failed %d\n%s", phaseMs, phaseStep + 1, clampAttempt, clampRounds, clampFails, clampStep.c_str());
            fflush(ClampFile());
         }
         clampStep.clear();
         clampAttempt = clampRounds = clampFails = 0;
      }
      bool SolveOnceCore(double h)
      {
         if (RoundTraceOn())
            roundBuf.clear();
         if (ContactTraceNow())
            contactBuf.clear();
         jacWinTry.clear();
         jacWinCon.clear();
         jacWinTryEntries = 0;
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
         double scaleFull = 0.0;           // stick_zone_trial: the move's full scale; its stick-zone candidates szCand[iCand ..] still to try, then the full move
         int nCand = 0, iCand = 0;
         bool fullTried = true;
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
            jacWinDump = set.jacMode != 0 && jacWinLines < 2000 && JacWindow(phaseMs);
            if (jacWinDump)   // VPX_GROUP_JAC_WINDOW: the balls at the start and this round's move (read only)
            {
               jacWinCon.clear();
               char ln[300];
               for (int i = 0; i < nb; i++)
               {
                  snprintf(ln, sizeof ln, "  ball %d v (%.5g %.5g %.5g) m/s w (%.5g %.5g %.5g) rad/s move (%.5g %.5g %.5g) um unknowns 3-5 (%.5g %.5g %.5g)\n", i, ball[i].v.x, ball[i].v.y, ball[i].v.z,
                           ball[i].w.x, ball[i].w.y, ball[i].w.z, x[6 * i] * 1e6, x[6 * i + 1] * 1e6, x[6 * i + 2] * 1e6, x[6 * i + 3], x[6 * i + 4], x[6 * i + 5]);
                  jacWinCon += ln;
               }
            }
            clampDump = ClampFile() != nullptr;
            clampCon.clear();
            const bool evOk = Eval(h, x, Res, J, trialCached);
            jacWinDump = false;
            clampDump = false;
            if (!evOk)
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
            if (set.jacMode != 0 && jacWinLines < 2000 && JacWindow(phaseMs))   // VPX_GROUP_JAC_WINDOW (diagnostics, read only)
            {
               char ln[240];
               snprintf(ln, sizeof ln, "JACWIN t %u ms h %.4g us round %d rn %.6g evals %ld%s\n", phaseMs, h * 1e6, round, rn, evalCount - evalCountRound,
                        round > 1 && !gripMoved && rn > rnPrev * (1.0 - 1e-4) && rn > set.tolV ? " (no drop: backtrack)" : "");
               jacWinTry += ln;
               evalCountRound = evalCount;
               const Jac Ja = JacAna(), Jf = JacFd();
               for (int r = 0; r < nU; r++)
                  for (int c = 0; c < nU; c++)
                  {
                     const double a = Ja[r][c], f = Jf[r][c], d = std::fabs(a - f);
                     if (d > 1e-6 && d > 0.1 * std::max(std::fabs(a), std::fabs(f)))
                     {
                        snprintf(ln, sizeof ln, "  row %s col %s a %.6g f %.6g\n", UnkName(r).c_str(), UnkName(c).c_str(), a, f);
                        jacWinTry += ln;
                        jacWinTryEntries++;
                     }
                  }
            }
            if (ClampFile() != nullptr)   // VPX_GROUP_CLAMP_TRACE: this round's point and verdict (the same test as the line search below), then its contacts
            {
               const bool rej = round > 1 && !gripMoved && rn > rnPrev * (1.0 - 1e-4) && rn > set.tolV && back < 30;
               char ln[200];
               snprintf(ln, sizeof ln, "R t %u step %ld att %d round %d %s %s rn %.6g h %.4g us\n", phaseMs, phaseStep + 1, clampAttempt, round, clampRetake ? "retake" : "move",
                        rej ? "rejected" : "kept", rn, h * 1e6);
               clampStep += ln;
               clampStep += clampCon;
               clampRounds++;
               clampRetake = false;
            }
            // (a part's grip changed since the last round: the residual moved for that too, so no comparison this round)
            if (round > 1 && !gripMoved && rn > rnPrev * (1.0 - 1e-4) && rn > set.tolV && back < 30)
            {
               bool ok = false;
               while (!ok && (iCand < nCand || !fullTried || scalePrev > 1e-9))
               {
                  if (iCand < nCand)   // stick_zone_trial: the next stick-zone candidate, then the full move, then halves of it
                     scalePrev = scaleFull * szCand[iCand++];
                  else if (!fullTried)
                  {
                     scalePrev = scaleFull;
                     fullTried = true;
                  }
                  else
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
                  clampRetake = true;   // (VPX_GROUP_CLAMP_TRACE: the next round evaluates the shorter move)
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
            // stick_zone_trial: the move's convergence and progress measures (ch) stay the full move's; the first trial is the first stick-zone candidate in range
            scaleFull = scale;
            nCand = iCand = 0;
            fullTried = true;
            if (set.stickZoneTrial && moved)
            {
               nCand = StickCandidates(stp, scale, h);
               while (iCand < nCand)
               {
                  const double sc = scale * szCand[iCand++];
                  for (int a = 0; a < nU; a++)
                     trial[a] = x[a] + stp[a] * sc;
                  if (InRange(trial))
                  {
                     scale = sc;
                     fullTried = false;
                     break;
                  }
               }
               if (fullTried && nCand > 0)   // (none in range: the full move, its range check last so its probes are the cached ones)
               {
                  for (int a = 0; a < nU; a++)
                     trial[a] = x[a] + stp[a] * scale;
                  InRange(trial);
               }
            }
            for (int a = 0; a < nU; a++)
            {
               xPrev[a] = x[a];
               stpPrev[a] = stp[a];
               if (moved)
                  x[a] = trial[a];
            }
            scalePrev = scale;
            // the parts' grip on the new move (held through the next round's Newton move); with grip_in_step each Eval works it out itself: nothing moves between rounds
            double gch = 0.0;
            if (!set.gripInStep)
               for (int p = 0; p < np; p++)
                  gch = std::max(gch, part[p]->UpdateGrip(x, h));
            gripMoved = gch > set.tolV;
            rd.scale = scale;
            rd.grip = gch;
            if (RoundTraceOn())
            {
               char line[400];
               int n = snprintf(line, sizeof line, "  r%d h %.4g us rn %.4g |", round, h * 1e6, rn);
               for (int i = 0; i < nb && n < 300; i++)
               {
                  const double ic = 1.0 / ball[i].m, iw = ball[i].R / ball[i].I();
                  n += snprintf(line + n, sizeof line - n, " b%d v %.3g %.3g %.3g w %.3g %.3g %.3g |", i, Res[6 * i] * ic, Res[6 * i + 1] * ic, Res[6 * i + 2] * ic, Res[6 * i + 3] * iw, Res[6 * i + 4] * iw, Res[6 * i + 5] * iw);
               }
               for (int c = 0; c < nCoord && n < 360; c++)
                  if (cUnk[c] >= 0)
                     n += snprintf(line + n, sizeof line - n, " q%d %.3g step %.3g |", c, Res[cUnk[c]] * ci[c].invM * ci[c].reach, stp[cUnk[c]] * scale);
               snprintf(line + n, sizeof line - n, " ch %.3g grip %.3g settled %d moved %d", ch, gch, settled ? 1 : 0, moved ? 1 : 0);
               roundBuf += line;
               for (int p = 0; p < np; p++)
                  roundBuf += " || " + part[p]->RowNote();
               roundBuf += "\n";
               roundBuf += probeTrace;
               if (set.jacMode != 0)   // VPX_GROUP_JAC_MODE: each part coordinate's row and column, the solve's own Jacobian (a) vs central differences (f), against itself and each ball it couples to (move x y z, spin x y z)
                  for (int c = 0; c < nCoord; c++)
                  {
                     const int qi = cUnk[c];
                     if (qi < 0)
                        continue;
                     const Jac Ja = JacAna(), Jf = JacFd();
                     char jl[200];
                     snprintf(jl, sizeof jl, "    jac q%d (round %d): qq a %.6g f %.6g\n", c, jacRound, Ja[qi][qi], Jf[qi][qi]);
                     roundBuf += jl;
                     const double floorF = 1e-9 * std::fabs(Jf[qi][qi]);
                     for (int i = 0; i < nb; i++)
                     {
                        bool any = false;
                        for (int a = 0; a < 6; a++)
                           any = any || Ja[qi][6 * i + a] != 0.0 || Ja[6 * i + a][qi] != 0.0 || std::fabs(Jf[qi][6 * i + a]) > floorF || std::fabs(Jf[6 * i + a][qi]) > floorF;
                        if (!any)
                           continue;
                        for (int side = 0; side < 2; side++)   // 0: the coordinate's row against the ball's unknowns; 1: the ball's rows against the coordinate
                        {
                           int n = snprintf(jl, sizeof jl, "      %s b%d:", side == 0 ? "row q vs" : "col q in", i);
                           for (int a = 0; a < 6 && n < 180; a++)
                           {
                              const double va = side == 0 ? Ja[qi][6 * i + a] : Ja[6 * i + a][qi], vf = side == 0 ? Jf[qi][6 * i + a] : Jf[6 * i + a][qi];
                              n += snprintf(jl + n, sizeof jl - n, " %.4g/%.4g", va, vf);
                           }
                           snprintf(jl + n, sizeof jl - n, "\n");
                           roundBuf += jl;
                        }
                     }
                  }
            }
            if (ContactTraceNow() && (int)cDiag.size() >= ni)
            {
               // the values of this round's evaluation at x (before its Newton move): forces N, depths um, slip speed mm/s, generalized forces on the part's coordinate (N m for a hinge)
               char line[600];
               snprintf(line, sizeof line, "  r%d h %.4g us rn %.4g ch %.3g scale %.3g grip %.3g moved %d\n", round, h * 1e6, rn, ch, scale, gch, moved ? 1 : 0);
               contactBuf += line;
               for (int k = 0; k < ni; k++)
               {
                  const Item &it = item[k];
                  const ItemEnd &e = ie[k];
                  const ContactDiag &cd = cDiag[k];
                  int n = snprintf(line, sizeof line, "    c%d b%d %s%d", k, it.a, it.b >= 0 ? "ball b" : "s", it.b >= 0 ? it.b : it.s);
                  if (!e.active)
                     snprintf(line + n, sizeof line - n, " off d0 %.4g d1 %.4g um\n", e.d0 * 1e6, e.d1 * 1e6);
                  else
                     snprintf(line + n, sizeof line - n, " d0 %.4g d1 %.4g dm %.4g um n %.4f %.4f %.4f (nz %.4f) N %.5g (spring %.5g damp %.4g) Ft %.5g lim muN %.5g trial %.5g %s slip %.4g mm/s | q-force N %.4g F %.4g\n",
                              e.d0 * 1e6, e.d1 * 1e6, cd.dm * 1e6, cd.n.x, cd.n.y, cd.n.z, cd.n.z, cd.N, cd.fe, cd.fd, cd.ft, cd.lim, cd.ftTr, e.slipping ? "SLIDE" : "stick", cd.ut * 1e3, cd.qN, cd.qF);
                  contactBuf += line;
               }
               for (int p = 0; p < np; p++)
               {
                  snprintf(line, sizeof line, "    p%d ", p);
                  contactBuf += line + part[p]->ContactNote(xPrev, partOff[p], h, Res) + "\n";   // (xPrev: the point this round evaluated)
               }
            }
            if (searchTrace)
            {
               char line[200];
               snprintf(line, sizeof line, "  round %d h %.6g us rn %.4g ch %.4g grip %.4g scale %.3g back %d", round, h * 1e6, rn, ch, gch, scale, back);
               searchBuf += line;
               for (int p = 0; p < np; p++)
                  searchBuf += " | part res " + std::to_string(Res[partOff[p]]) + " " + part[p]->SearchNote(x, partOff[p], h);
               searchBuf += "\n";
            }
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
         if (!jacWinTry.empty())   // (VPX_GROUP_JAC_WINDOW: a solve that converged)
         {
            if (jacWinLines < 2000)   // (the rounds' lines, then the contacts as the last round's Eval had them)
            {
               jacWinBuf += jacWinTry;
               jacWinBuf += "  contacts at the last round:\n" + jacWinCon;
               jacWinLines += jacWinTryEntries;
            }
            jacWinTry.clear();
         }
         if (endValuesHere)
            return true;
         return Eval(h, x, Res, nullptr, trialCached);   // the end values at the solution
      }

      // stick_zone_trial: the line search's first trials for Newton's move d x scale (the full move the line search tries first): per contact sliding at the point
      // the round evaluated (szRec), the middle of the part of (0, 1] where its trial force a + s b lies inside mu N(s), N(s) = N0 + s dN >= 0; b = C dut, dut the slip
      // rate's change over the move by ut's own form (ut = Pt (dc_a / h + wm_a x rA - dc_b / h - wm_b x rB - sum_j cV_j q_j / h), wm = w + dw / 2), the normal, levers,
      // kS, ct and stretch held at the point. A contact whose stick interval contains s = 1 (the full move already ends inside its zone) gives no candidate: only
      // one whose interval ends before s = 1 does. Into szCand, ascending, duplicates dropped; returns how many.
      int StickCandidates(const double *d, double scale, double h)
      {
         szCand.clear();
         for (int k = 0; k < ni; k++)
         {
            const StickRec &sr = szRec[k];
            if (!sr.on)
               continue;
            const int ia = 6 * sr.A;
            V3 du = V3(d[ia], d[ia + 1], d[ia + 2]) * (1.0 / h) + (V3(d[ia + 3], d[ia + 4], d[ia + 5]) * 0.5).Cross(sr.rA);
            if (sr.B >= 0)
            {
               const int ib = 6 * sr.B;
               du = du - (V3(d[ib], d[ib + 1], d[ib + 2]) * (1.0 / h) + (V3(d[ib + 3], d[ib + 4], d[ib + 5]) * 0.5).Cross(sr.rB));
            }
            for (int j = 0; j < sr.nq; j++)
               if (sr.qU[j] >= 0)
                  du = du - sr.cV[j] * (d[sr.qU[j]] / h);
            du = du - sr.nm * du.Dot(sr.nm);
            const V3 b = du * (sr.C * scale);
            double dN = 0.0;
            for (int i = 0; i < sr.nn; i++)
               dN += sr.nG[i] * d[sr.nIx[i]];
            dN *= scale;
            // |a + s b|^2 <= (mu (N0 + s dN))^2: qa s^2 + qb s + qc <= 0
            const double m2 = sr.mu * sr.mu;
            const double qa = b.Dot(b) - m2 * dN * dN, qb = 2.0 * sr.a.Dot(b) - 2.0 * m2 * sr.N0 * dN, qc = sr.a.Dot(sr.a) - m2 * sr.N0 * sr.N0;
            const double disc = qb * qb - 4.0 * qa * qc;
            if (disc < 0.0 || std::fabs(qa) < 1e-300)
               continue;
            const double sq = std::sqrt(disc);
            const double r0 = std::min((-qb - sq) / (2.0 * qa), (-qb + sq) / (2.0 * qa)), r1 = std::max((-qb - sq) / (2.0 * qa), (-qb + sq) / (2.0 * qa));
            const double ivLo[2] = { qa > 0.0 ? r0 : -INFINITY, r1 }, ivHi[2] = { qa > 0.0 ? r1 : r0, INFINITY };
            double best = INFINITY;
            bool endIn = false;
            for (int v = 0; v < (qa > 0.0 ? 1 : 2); v++)
            {
               double lo = ivLo[v], hi = ivHi[v];
               if (dN < 0.0)   // N(s) >= 0: the normal force cannot pull
                  hi = std::min(hi, -sr.N0 / dN);
               else if (dN > 0.0)
                  lo = std::max(lo, -sr.N0 / dN);
               if (lo <= 1.0 && hi >= 1.0)   // the full move already ends inside this contact's stick zone: no candidate from it
                  endIn = true;
               const double lo2 = std::max(lo, 0.0), hi2 = std::min(hi, 1.0);
               if (hi2 <= lo2 || hi2 <= 0.0)
                  continue;
               const double s = 0.5 * (lo2 + hi2);
               if (s > 0.0)
                  best = std::min(best, s);
            }
            if (!endIn && best < 1.0)
               szCand.push_back(best);
         }
         std::sort(szCand.begin(), szCand.end());
         szCand.erase(std::unique(szCand.begin(), szCand.end()), szCand.end());
         return (int)szCand.size();
      }

      void Accept(double h, bool onStop)
      {
         nearOn = false;   // (from here the balls are at the step's end: the probes take every piece)
         const bool hist = np > 0 && HistoryTraceOn();   // (diagnostics, read only)
         const HistE histS0 = hist ? HistEnergy() : HistE();
         const Books histBk0 = hist ? books : Books();
         if (ImpOn())   // (test tooling: the impulses of the Eval this step was accepted on)
            ImpAccept(h);
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
            if (b.brakeRate > 0.0)
            {
               const V3 wmid = b.w + dw * 0.5, wa = wmid - b.brakeAxis * wmid.Dot(b.brakeAxis);
               books.spinBrakeJ += b.brakeRate * b.I() * wa.Dot(wa) * h;
            }
            b.c = b.c + dc;
            b.v = EndRate(dc, b.v, h);
            b.w = b.w + dw;
            aLast[i] = evF[i] * (1.0 / b.m);         // the next step's first guess: the force and torque of the move taken (Eval's at it), not (v1 - v0) / h, which is
            alphaLast[i] = evM[i] * (1.0 / b.I());   // the solve's tolerance over h for a step near zero long (a stop search on a stop: 5.7e12 m/s^2, a guess 88 km off)
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
            {
               it.dZero = std::max(0.0, std::min(it.dZero, it.dZero + e.d1));
            }
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
         {
            PhaseScope ps(phPartAfter);
            for (int p = 0; p < np; p++)
               part[p]->AfterStep(*this);
         }
         if (hist)
            HistoryKeep(h, histS0, histBk0);
         if (hist)
            histClock += h;
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

   // A part with no ball this cycle: the same solve, no balls, on the engine's settings and step floor (EngineSettings); one function for every part. False: a step could not be solved (time is lost,
   // no energy can appear: the next cycle goes on from here). out: the group's books of this cycle, when wanted.
   inline bool AdvanceAlone(Part &p, const double dt, const int tMs, Books *const out = nullptr)
   {
      static Group g;
      g.BeginUse();
      g.set = EngineSettings();
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
   // A part's cycle alone, as every part's host runs it: the step, and the part's counters (a step not solved, the detector of rule 8; steps taken, cost; of them on energy alone, rule 8). The steps taken
   inline int AloneCycle(Part &p, const double dt, const int tMs, long &failedSteps, long &aloneSteps, long &energySteps)
   {
      Books bk;
      if (!AdvanceAlone(p, dt, tMs, &bk))
         failedSteps++;   // (time is lost, no energy can appear: the next cycle goes on from here)
      aloneSteps += bk.steps;
      energySteps += bk.energySteps;
      return (int)bk.steps;
   }
}
