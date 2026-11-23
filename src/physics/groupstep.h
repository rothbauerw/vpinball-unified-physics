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
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

#include "physics/rubbercontact.h"
#include "physics/pairform.h"
#include "physics/partrecords.h"

namespace ContactSystem
{
   void NoticeRetired(const std::string &msg);   // contactsystem.cpp: one WARN and one on-screen notice per table start for each text (link_windows: a link order net)
}

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
   // A force on the balls that depends on the ball's place but is not the slope of a stored energy the group can read (an air-core coil's pull from its force table): per unit of
   // the ball's Ball::forceScale, read at the step's midpoint in every evaluation; a drive, its work (midpoint force . move) counted. Fixed for the cycle; it owns no state of a step.
   // slope: dF/dc (N/m, row i = force component i), when wanted
   class BallForce
   {
   public:
      virtual ~BallForce() = default;
      virtual V3 At(const V3 &c, double R, M3 *slope) const = 0;
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
      int id = -1;   // the piece it comes from (the same piece at every place and in every list): its found-inside zero is kept under it (Item::tZero)
      const std::vector<int> *over = nullptr;   // a joint term: the pieces (ids) it corrects, its zero floored by theirs (jointZeroFloor); null for any other term
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
   struct HitPeaks   // peak_slope: a static surface's level-form nodes (RubberContact::PeakTerms), the first kInline in place
   {
      static constexpr int kInline = 8;
      RubberContact::PeakTerm &operator[](const int i) { return i < kInline ? in[i] : ex[(size_t)(i - kInline)]; }
      const RubberContact::PeakTerm &operator[](const int i) const { return i < kInline ? in[i] : ex[(size_t)(i - kInline)]; }
      void Fit(const int n) { if (n > kInline && ex.size() < (size_t)(n - kInline)) ex.resize((size_t)(n - kInline)); }
      RubberContact::PeakTerm in[kInline];
      std::vector<RubberContact::PeakTerm> ex;
   };
   // patch_friction (piece 2b-1): one pressed patch of a static surface: its depth (from the object's zero), normal, gradient and its pieces (Hit::pc[p0 .. p0 + pn - 1])
   struct HitPatch
   {
      double d, grad;
      V3 n;
      int p0, pn;
      int gk;      // its deepest piece's signed geometry (RubberContact::PieceGeo): 0 the plane through n, 1 the segment ga-gb, 2 the point ga
      V3 ga, gb;
   };
   template <class T, int K> struct InList   // the first K in place, the rest on the heap, grown as needed (storage, not a limit)
   {
      T &operator[](const int i) { return i < K ? in[i] : ex[(size_t)(i - K)]; }
      const T &operator[](const int i) const { return i < K ? in[i] : ex[(size_t)(i - K)]; }
      void Fit(const int n) { if (n > K && ex.size() < (size_t)(n - K)) ex.resize((size_t)(n - K)); }
      T in[K];
      std::vector<T> ex;
   };
   // step pieces (runs/2026-11-22-ringcap design 2.4, 2.7): a static object's other convex piece in reach (not the one its own answer comes from), or a ridge of its pressed
   // patches, found by the step-start probe (Shape::ProbeStart). For the step rule only: ChooseStep takes it as it takes an object (no force, no state kept)
   struct HitEntry
   {
      double d;      // its depth as Hit::delta (> 0 pressed, < 0 the gap): a piece's before the item's zero (ZeroHit takes it off), a ridge's final (zeroed)
      V3 n;          // its normal (surface -> ball centre): the closing speed is -v . n
      bool zeroed;   // true: d is already from the object's zero (a ridge: the pass's depth, or the edge's gap from R - zero)
   };
   // step pieces counters (always on, cheap; written to the VPX_GROUP_STATS file): step-start lookups of a static surface, entries they built, steps an entry set
   struct PieceStats
   {
      long lookups = 0, entries = 0, set = 0;
   };
   inline PieceStats pieceStats;
   // patch_zero counters (running totals logged at 1, 2, 4, ...): items on the per-patch path at a cycle's start, first-sight patches, rule-5 joins at Accept, start merges
   // (F4), nets (a probe's raw patch with different zeros, F2)
   inline void PzCount(std::atomic<long long> &n, const char *what)
   {
      const long long was = n.fetch_add(1), c = was + 1;
      if (was == 0 || (c ^ was) > was)   // a power of 2 crossed
         PLOGI << "patch_zero: " << c << " " << what;
   }
   inline std::atomic<long long> *PzStats()   // [0] items per cycle, [1] first-sight patches, [2] rule-5 joins, [3] F4 start merges, [4] F2 nets
   {
      static std::atomic<long long> c[5] { { 0 }, { 0 }, { 0 }, { 0 }, { 0 } };
      return c;
   }
   // touch_predict counters (always on, cheap; written to the VPX_GROUP_STATS file): steps a predicted touch set and ended (predicted), of them landed within d and closing
   // (taken) or more than d short (shortL); re-solves at the Hermite root of a step that ended past a predicted touch (hermite); re-solves at today's limit after that (net).
   // Each by kind (piece 2): [0] a fixed surface, [1] a moving part's surface (a union's term with it), [2] a ball pair; a re-solve counts under the kind whose touch set its length
   struct TouchStats
   {
      long predicted[3] = {}, taken[3] = {}, shortL[3] = {}, hermite[3] = {}, net[3] = {};
      long partEvals = 0;   // the parts' start accelerations formed (ChooseStep, PartAccel: one BeginStep and one Eval each)
      static long Sum(const long (&c)[3]) { return c[0] + c[1] + c[2]; }
   };
   inline TouchStats touchStats;
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
      int np = 0;              // peak_slope: > 0 a static surface's level form with 2 or more peaks, its peaks and passes (pk), for the elastic push's Newton slope per peak
      HitPeaks pk;
      int npat = 0;            // patch_friction: > 0 a static surface's pressed patches (pat), their pieces (pc, npc ids)
      InList<HitPatch, 4> pat;
      int npc = 0;
      InList<int, 16> pc;
      InList<RubberContact::PieceData, 16> pd;   // level_stiffness: each piece's depth (from the object's zero, as pat's), normal, geometry and neighbours (parallel to pc)
      int ne = 0;              // step pieces: the step-start probe's entries (en): the object's other convex pieces and ridges for the step rule; 0 from every other probe
      InList<HitEntry, 8> en;
      Hit() {}
      Hit(const Hit &o) { *this = o; }
      Hit &operator=(const Hit &o)
      {
         ne = o.ne;
         en.Fit(ne);
         for (int j = 0; j < ne; j++)
            en[j] = o.en[j];
         npat = o.npat;
         pat.Fit(npat);
         for (int j = 0; j < npat; j++)
            pat[j] = o.pat[j];
         npc = o.npc;
         pc.Fit(npc);
         pd.Fit(npc);
         for (int j = 0; j < npc; j++)
         {
            pc[j] = o.pc[j];
            pd[j] = o.pd[j];
         }
         np = o.np;
         pk.Fit(np);
         for (int j = 0; j < np; j++)
            pk[j] = o.pk[j];
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
   // peak_slope (plan runs/2026-11-22-per-peak 2.5a): the merged depth's Hessian for a level form of peaks and passes, e = s25^0.4, s25 = sum s_i d_i^2.5:
   // H = (M - 1.5 e^0.5 gm gm^T) / e^1.5, M = sum s_i (1.5 d_i^0.5 n_i n_i^T - d_i^1.5 (I - n_i n_i^T - t_i t_i^T) / r_i) (r_i = 0: n_i does not turn; t_i a line piece's axis), gm the merged depth's gradient on
   // the ball's centre (-n grad). One peak with no turn: M = 1.5 d^0.5 n n^T = 1.5 e^0.5 gm gm^T, H = 0 (to the rounding of e = (d^2.5)^0.4). False: no pressed sum.
   inline bool PeakHessian(const HitPeaks &pk, const int np, const double *gm, double H[3][3])
   {
      double s25 = 0.0, M[3][3] = {};
      for (int i = 0; i < np; i++)
      {
         const RubberContact::PeakTerm &t = pk[i];
         const double sd = std::sqrt(t.d), a = t.s * 1.5 * sd, b = t.r > 0.0 ? t.s * t.d * sd / t.r : 0.0;
         const double n[3] = { t.n.x, t.n.y, t.n.z }, ax[3] = { t.t.x, t.t.y, t.t.z };   // C4: a line piece turns by (I - n n^T - t t^T) / r (t = 0: a joint, vertex)
         s25 += t.s * t.d * t.d * sd;
         for (int p = 0; p < 3; p++)
            for (int q = 0; q < 3; q++)
               M[p][q] += (a + b) * n[p] * n[q] - b * ((p == q ? 1.0 : 0.0) - ax[p] * ax[q]);
      }
      if (s25 <= 0.0)
         return false;
      const double e = std::pow(s25, 0.4), se = std::sqrt(e), c = 1.5 * se, k = 1.0 / (e * se);
      for (int p = 0; p < 3; p++)
         for (int q = 0; q < 3; q++)
            H[p][q] = (M[p][q] - c * gm[p] * gm[q]) * k;
      return true;
   }
   inline bool PeakSlopeCheckEnv() { static const bool on = []() { const char *const e = getenv("VPX_PATCH_CHECK"); return e != nullptr && atoi(e) != 0; }(); return on; }

   class Shape
   {
   public:
      RubberParams rp;         // the material (rubbercontact.h)
      double layerM = 0.0;     // the rubber's thickness for the layer term, m (0 = no layer term)
      int onlyBall = -1;       // >= 0: only this ball meets it (a surface's list of near faces is one ball's)
      int part = -1;           // the part that carries it, its place among the group's parts; -1 = fixed. Its Probe fills the same columns, in the same order, every time
      bool patchZero = false;  // patch_zero: a static mesh or tube (a fixed surface's own shape): with Settings::patchZero its item's zero is kept per patch (Item::tZero)
      mutable const double *trial = nullptr;   // a carried surface: each of its part's coordinates moved by this for the probe; null = where they are now
      virtual ~Shape() = default;
      virtual void Probe(const V3 &c, double R, Hit &h) const = 0;
      // the step-start probe (ChooseStep's, only): Probe, and a static surface also lists its other convex pieces and ridges in reach as entries (Hit::en; shared_step_only)
      virtual void ProbeStart(const V3 &c, const double R, Hit &h) const { Probe(c, R, h); }
      virtual std::string Name() const { return std::string(); }   // diagnostics (VPX_GROUP_CLAMP_TRACE): its owner's name
      // link_windows (rule T, to-do 1.33): a static surface's links pressed on the step's straight move c0 -> c1, each with its window (ContactSurface::LinkWindows); none by default
      virtual void LinkWindows(const V3 & /*c0*/, const V3 & /*c1*/, std::vector<RubberContact::LinkWindow> &out) const { out.clear(); }
      // the path test: does the ball's centre moving straight from s to s + d come within reach of a piece (depth above 0 anywhere on the path, before the item's zero)?
      // 1 = yes, 0 = no (exact, the path kit in rubbercontact.h; a curved piece by its cover, PathRevolve), -1 = no exact test for this shape: FarStillFar probes the middle and the end as before.
      // (A surface its part carries: PathMotion below, the rigid-body law in FarStillFar.)
      virtual int PathReach(const V3 & /*s*/, const V3 & /*d*/, double /*R*/) const { return -1; }
      // the part finder, a surface its part carries (group A, runs/2026-11-13-finder-P6/brief_groupA.md): trial = its part's move over the step (the end trial). The most m (m) any
      // point of the surface moves relative to the point P over the step (the ball's centre moving dLen): a turn about an axis (rho(P) + dLen / 2) |dq|, a turn about a pivot point
      // (|P - pivot| + dLen / 2) |dtilt|, a translation |dx|, sums for combined motions (paper: runs/2026-11-13-finder-P6/paper_groupA.txt). -1 = none: the two probes
      virtual double PathMotion(const V3 & /*P*/, double /*dLen*/) const { return -1.0; }
      // the largest slope of its probe's distance (1 for a true distance; a tapered Revolve profile more: JetBumper::Revolve::Slope), built at load; the law multiplies the moves by it
      virtual double PathSlope() const { return 1.0; }
      // the self-check of a cover (VPX_FINDER_CHECK=1, diagnostics only): its largest sagitta (m) and kind (0 post, 1 bumper body, 2 disc);
      // -1 = not a cover
      virtual double PathCoverSag(int & /*kind*/) const { return -1.0; }
      // the part finder, a surface of pieces: the cycle's list (in the pieces' own order) for a ball whose centre was c0 at the cycle's start and can move reach (m) in the cycle,
      // its part turning within the range its CycleTurn set; pieces: how many the surface has. False: not a surface of pieces (its probe takes no list)
      virtual bool CycleNear(const V3 &, double, double, std::vector<int> &, int &) const { return false; }
      // check mode (VPX_FINDER_CHECK): the pieces within R of a ball at c, its part moved by dq (the rigid piece at the trial angle), that are not on list (exact distances)
      virtual int CycleMisses(const V3 &, double, double, const std::vector<int> &) const { return 0; }
   };
   // part finder's counters per kind of part (its type's name; diagnostics only: VPX_CHECK_PROFILE or VPX_FINDER_CHECK, set by the owner in partFinderProf)
   struct PartFinderCount
   {
      long long lists = 0, listLen = 0, pieces = 0;   // lists built (per cycle and after a re-query), their pieces, the surface's pieces
      long long requeries = 0, rangeRuns = 0;   // re-queries (a ball outside its query volume), range re-runs (group A)
      long long checks = 0, misses = 0;              // check mode: positions checked (start, middle, end of each solve), pieces within R off the list
      long long ovCalls = 0, ovProbes = 0, ovSkips = 0, ovMismatch = 0;   // CheckOverlaps: calls (each part counted once per call), probes, parts skipped by their region, skipped parts overlapping (check mode)
      // FarStillFar's skipped contacts on this kind's surfaces (P6a2): tested on the path by the swept region, cleared, reached; left to the two probes (no region: -1); fixed surfaces
      // of the part with their own test (a bumper's body, a disc: counted also in "Path test totals"); check mode: swept-region self-checks (path end and middle) and disagreements
      long long pathTested = 0, pathClear = 0, pathReached = 0, pathTwo = 0, pathFixed = 0, pathChecks = 0, pathBad = 0;
   };
   inline std::map<std::string, PartFinderCount> partFinderCount;
   // energy skip test: static items of an own ball (all its ball-pair and moving-part items skipped in the attempt) not skipped by its cycle speed, tested on its own bound; skipped by it (diagnostics only)
   inline long long skipOwnTested = 0, skipOwnSkipped = 0;
   // energy skip test: own balls whose stored static-contact energy E_static was summed (at most once per ball per attempt, only when an own test needs it; diagnostics only)
   inline long long skipEsComputed = 0;
   // energy skip test: in one-trip attempts (BeginAttempt), the items looked at, the static items left past the cycle's speed for the own test after the trip (deferred), and the attempts (diagnostics only)
   inline long long skipTripItems = 0, skipDeferred = 0, skipTripAttempts = 0;
   inline bool partFinderProf = false;
   inline std::string partFinderFirstMiss;   // check mode: the first miss, for the owner's report
   inline std::string partFinderFirstPathBad;   // check mode: the swept region's first disagreement with its probe (P6a2)
   inline bool FinderCheckEnv() { static const bool on = []() { const char *const e = getenv("VPX_FINDER_CHECK"); return e != nullptr && atoi(e) != 0; }(); return on; }
   inline std::string PartFinderReport()
   {
      std::string s = "Part finder totals:";
      char b[900];
      for (const auto &e : partFinderCount)
      {
         const PartFinderCount &c = e.second;
         snprintf(b, sizeof b, " | %s: lists %lld, pieces per list %.3f, pieces per surface %.3f, re-queries %lld, range re-runs %lld, checks %lld, misses %lld; overlap check calls %lld, probes %lld (%.3f per call), skipped %lld, mismatches %lld"
                  "; path: swept tested %lld, cleared %lld, reached %lld, two probes %lld, fixed surface %lld, self-checks %lld, disagreements %lld",
                  e.first.c_str(), c.lists, c.lists > 0 ? (double)c.listLen / (double)c.lists : 0.0, c.lists > 0 ? (double)c.pieces / (double)c.lists : 0.0, c.requeries, c.rangeRuns, c.checks, c.misses,
                  c.ovCalls, c.ovProbes, c.ovCalls > 0 ? (double)c.ovProbes / (double)c.ovCalls : 0.0, c.ovSkips, c.ovMismatch, c.pathTested, c.pathClear, c.pathReached, c.pathTwo, c.pathFixed,
                  c.pathChecks, c.pathBad);
         s += b;
      }
      snprintf(b, sizeof b, " | energy_skip own bound: static items tested %lld, skipped by it %lld, E_static sums %lld; one trip: attempts %lld, items %lld, deferred %lld", skipOwnTested, skipOwnSkipped,
               skipEsComputed, skipTripAttempts, skipTripItems, skipDeferred);
      s += b;
      if (!partFinderFirstMiss.empty())
         s += " | first miss: " + partFinderFirstMiss;
      if (!partFinderFirstPathBad.empty())
         s += " | first swept-region disagreement: " + partFinderFirstPathBad;
      return s;
   }
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
      int PathReach(const V3 &s, const V3 &d, const double R) const override   // the height over the plane is linear along the path: its least is at an end
      {
         const double hs = (s - p).Dot(n), he = (s + d - p).Dot(n);
         return R - std::min(hs, he) > 0.0 ? 1 : 0;
      }
   };
   struct Ball
   {
      V3 c, v, w;            // centre, velocity, spin (m, m/s, rad/s)
      int adds = 0;          // the additions to c since the finder's query centre c0 (each stored c = fl(c + dc) rounds: the inside check allows it, FinderInsideQ); saved and put back with c
      double m = 0.080, R = 25.0 * 0.53975e-3;
      V3 fExt;               // a held outside force (a magnet not on a field, the magnets' ball-pair force, a part not yet moved over): a drive, its work counted
      double fieldScale = 0.0;   // how much of the group's fields it feels (Group::field); 0 = none
      double forceScale = 0.0;   // how much of the group's ball forces it feels (Group::force); 0 = none (a held ball: none)
      double brakeRate = 0.0;    // a viscous brake on its spin across an axis (the eddy currents in a coil's field), 1/s; 0 = none. The caller's, held for the cycle
      V3 brakeAxis;              // the axis (unit): a spin about it is not braked
      // held (KickerDesign 1h): a ball a kicker holds (parked) or carries (the arm) has infinite mass: forces do not move it (inverse mass and inertia 0), its move is given: velocity v
      // and acceleration heldA (m/s2, the arm's path's; 0 parked), spin 0, no gravity. Its rows in the solve are that move; its own energy is the kicker's, not the group's; the work the
      // contacts' forces on it do over its move is a drive (books.driveJ: 0 for a parked ball). The caller puts held balls after the others, so a pair with one has it as its b.
      bool held = false;
      V3 heldA;
      // the energy reach: from the contact system's grouping this cycle (ReachLinks): the energy that can feed the ball's island (its balls' kinetic and stored energy plus the
      // drive work of its parts, J) and the ball's travel over the cycle (m, EnergyReach); -1 = not given (the island balls' kinetic energy and ContactReachM)
      double feedE = -1.0, reachS = -1.0;
      // energy skip test: its top speed over the cycle (m/s; a held ball's given speed included) and its group's energy at that speed (J), its energy group; -1 = not given (skip nothing)
      double vTop = -1.0, eTop = -1.0;
      int eGroup = -1;
      // energy skip test: the bound on its acceleration by forces other than its contacts (m/s2): |gravity + nudge| + |a magnet's pull| now (ContactSystem Sources::Begin; 0 held);
      // -1 = not given (no own bound: the cycle's vTop)
      double aOwn = -1.0;
      double I() const { return 0.4 * m * R * R; }
      double InvM() const { return held ? 0.0 : 1.0 / m; }      // inverse mass (0: held)
      V3 Given(const double h) const { return v * h + heldA * (0.5 * h * h); }   // a held ball's move over h
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
      // touchPredict ([Contact] touch_predict): a ball not touching a fixed surface, a moving part's surface (with the part's start acceleration, PartAccel) or another ball gets
      // no clock; its step ends at the predicted touch (ChooseStep) and the solved end is checked (Step, TouchEnd). false = today's arrival rule and clock. Not under VPX_GROUP_HFIX.
      bool touchPredict = true;
      // dampDirSlope ([Contact] damp_dir_slope): a contact made of terms damps along G = Fel / |Fel| (the pressing terms' summed elastic force); G turns as the terms'
      // forces change, at h dG/du = (I - G g^T) K / fe (K = sum_t h tk_t tg_t tg_t^T, g = G's ball part). Newton's slope then takes fd h dG/du and dd = G . dq's change
      // through G, (kd / h) G dq^T dG/du; N's slope (kv) the latter. The force is unchanged. true (default) = on; false = G held fixed in the slope.
      bool dampDirSlope = true;
      // peakSlope ([Contact] peak_slope): 1: the elastic push's Newton slope summed per peak, each with its own normal and turn (plan runs/2026-11-22-per-peak 2.5a); 0: today's
      // slope along the merged normal. A static surface with 2 or more peaks only (Hit::np); the push itself is unchanged.
      bool peakSlope = true;
      // levelStiffness ([Contact] level_stiffness, patch_friction only): 1: a friction group's shear, rolling and twist stiffness from the level form over its pieces, averaged over the
      // step with each piece's depth straight between its ends (B-lin, plan runs/2026-11-22-per-peak 2.4b; one peak: today's PatchMeans); 0: from the group's merged depth
      bool levelStiffness = true;
      // linkWindows ([Contact] link_windows, patch_friction only; rule T, to-do 1.33, runs/2026-11-23-link-rule brief_code.md): 1: a step whose start has 2 or more patches, or a start or
      // end patch sharing no piece with one at the other end, links its friction groups also through every link (shared edge, tube joint) pressed on the step's straight move (the
      // union), and a step whose time-ordered feed (rule T's forward sweep) differs from the union's is cut and solved again (at most 4 cuts); 0: start and end patches linked by
      // shared pieces only
      bool linkWindows = true;
      // patchZero ([Contact] patch_zero, multi_contact and link_windows; to-do 1.33 piece 2b-3, runs/2026-11-23-zero-patch brief_code.md): 1: a static item with a zero keeps it per patch
      // (Item::tZero: each raw-pressed piece of a zeroed patch with its patch's zero); its probe gets them (RubberContact::PatchZeros) and returns each patch shifted by its own zero
      bool patchZero = true;
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
      // termEnergy ([Contact] term_energy; TODO 48 a): a terms contact's stored energy (Item::Energy, every book) is sum_t s_t V(max(d_t, 0)), the potential the terms' force
      // (Elastic: sum_t s_t Chord(a0_t, a1_t)) conserves. false = V of the union's depth (they agree only for pure Hertz and a positive sum).
      bool termEnergy = true;
      // jointZeroFloor ([Contact] joint_zero_floor; TODO 48 b): a joint term's zero >= the largest zero of the pieces it corrects (HitTerm::over), at first sight, follow-out
      // and the raise between cycles (FloorJointZeros), so d_J <= min(d_pieces) and U >= 0. false = each term's zero on its own.
      bool jointZeroFloor = true;
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
      double arrivalShare = 0.5;               // a contact or a term closing at a steady rate is reached in about two steps, not overshot: a step is at most this share of the time to get there (rule 3; StepFactList::Arrival)
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

   // VPX_STEP_SCALE=<k> (a diagnostic, read once; runs/2026-11-10-step-scale/design.md): every limit that sets a step for accuracy is k times its own (k < 1 = shorter steps, the accuracy bench's
   // convergence runs): omegaH, omegaHRest, arrivalShare, turnMaxRad, tauShare and fieldShare, and the at-rest phase budget epsRest by k^2 (its step goes as sqrt(epsRest)), set in
   // ContactSettings::ToStep; the old island loop's omega * h limits (contactsystem.cpp SetupIsland, NextStep) likewise. Not scaled: exact event times, the 1 ms cycle, the step floor, vSteady.
   // Unset, empty or <= 0 = 1 (nothing changes). Use k >= 0.25: the cost grows about 1/k and a cycle may pass the step guard's 3000 steps. With VPX_GROUP_HFIX set the fixed step wins.
   inline double GroupStepScale()
   {
      static const double k = []() { const char *const e = getenv("VPX_STEP_SCALE"); const double v = e != nullptr && *e ? atof(e) : 0.0; return v > 0.0 ? v : 1.0; }();
      return k;
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
      void Arrival(const double gap, const double rate, const double shortest, const int why) { Lim(std::max(set.arrivalShare * gap / rate, shortest), why); }   // closing on something at this rate; no shorter than this
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

   // A massless coordinate's contacts as its first guess sees them (Group::MasslessGuess, Part::MasslessGuess): Push(dq) their push along it at the step's end (N) when it
   // moves dq, by the full end-push law (Group::EndNormal, the law the solve's row takes) along each contact's depth linearized at the step's start in dq, the balls' guessed
   // moves and every other unknown held; dPush its slope (N/m; <= 0 where the law's force rises with depth: each push falls as the part moves on). Law evaluations only, no
   // probe. start: the move that keeps the depth of the contacts touching at the step's start (-inf with none); tol: the solver's tolerance on depth as a move of it (m)
   struct MasslessRow
   {
      virtual double Push(double dq, double &dPush) const = 0;
      double start = -INFINITY, tol = 0.0;
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
      // the part finder: at the cycle's start, how far it can turn in the cycle of T s with the balls of its island holding E J of kinetic energy (design_reach.md section 1). None: no lists
      virtual void CycleTurn(double /*E*/, double /*T*/) {}
      double eOwnJ = -1.0;   // energy skip test: what it can hold this cycle before a ball touches it (its kinetic energy, drive work, squash; J); -1 = not given (the host sets it each cycle)
      // the part finder (group A): the range of its coordinates its cycle lists were sized for is a size only, made complete by this check after each solve: false = the solved
      // move qe (the end trial) left it; CycleRegrow then sets the range again from its state now (CycleTurn over T with E) joined to the old range and the solved place
      virtual bool CycleRangeHolds(const double * /*qe*/) const { return true; }
      virtual void CycleRegrow(const double * /*qe*/, double /*E*/, double /*T*/) {}
      // what it is
      virtual int NumDof() const = 0;
      virtual int MaxDof() const { return NumDof(); }   // the most its NumDof can be inside a cycle (the group's storage is sized between cycles)
      virtual int NumCoords() const = 0;                // constant between cycles
      virtual int DofOf(int coord) const = 0;           // the coordinate's place among its unknowns now; -1 = it has none now
      virtual CoordInfo Coord(int coord) const = 0;
      // the part's stored energy and its own bodies' kinetic energy at the balls' positions b (the balls' own energy is the group's)
      virtual double Energy(const Ball *b, int nb) const = 0;
      // the part's own force on each ball at the balls' positions b and its state now, added to F[i] (one per group ball): only a force its Eval puts on a ball directly, not
      // through an item, and no rate term (touch_predict's net force, ChooseStep; runs/2026-11-22-touch-predict/design.md 2.7). Writes no state. None: nothing
      virtual void BallForceNow(const Ball * /*b*/, int /*nb*/, V3 * /*F*/) const {}
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
      virtual double StopSquashJ() const { return 0.0; }   // the energy reach: what its compliant stops hold above their unloaded state, by each stop's own law (J); a stop that acts by restitution holds nothing
      virtual double PendingNetJ() const { return 0.0; }   // drive less loss already in its state but not yet booked (an impact taken before the cycle's book opens, booked at the next Accept): the book opens before it (J)
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
      // a part with a massless coordinate (CoordInfo invM 0 with an unknown: a stand-up target's face, standupface.h): its stored energy's change over the move x (J), which the
      // balls touching it pay (Group::MasslessEnergy; runs/2026-11-19-standup-p2/design.md section 3 b); slope (not null): its change with each of its unknowns, added at off + d
      virtual double MasslessStoredChange(const double *, int, double *) const { return 0.0; }
      // its massless coordinate d's first guess (Group::MasslessGuess, design section 1 a: a massless face follows the ball): its own row solved with its contacts' end push
      // (row: MasslessRow), the balls' guessed moves held; into x at off + d
      virtual void MasslessGuess(int, const MasslessRow &, double *, int) {}

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
      std::vector<RubberContact::PieceMem> pm;   // patch_friction (plan 2.1): the memory kept on the pieces pressed at the last accepted step's end (each its patch group's)
      bool pmMulti = false;                      // that step had 2 or more friction groups: the next step reads its memory from pm, not xi / ra
      double pmEs = 0.0, pmEr = 0.0;             // pmMulti: the springs' stored energy, the sum over its groups (shear; rolling and twist), J: the book's (ShearE, RollE)
      bool kStored = false;                      // plan 2.4's memory rule: xi, ra were stored at kt1, kr1, kw1 (the last accepted step's one friction group's end stiffness)
      double kt1 = 0.0, kr1 = 0.0, kw1 = 0.0;
      double kReff = 0.0;                        // (diagnostics: the contact's reff when kt1 .. were stored)
      bool pmLev = false;                        // level_stiffness: that step had one friction group on the level form (B-lin): pmEs, pmEr hold its springs at its end stiffness
      double ShearE(const double d) const { return pmMulti || pmLev ? pmEs : 0.5 * Kt(d) * xi.Dot(xi); }
      double RollE(const double d) const { return pmMulti || pmLev ? pmEr : 0.5 * Kr(d) * ra.th.Dot(ra.th) + 0.5 * Ktw(d) * ra.tw * ra.tw; }
      double dZero = 0.0;             // the depth its squash counts from: a ball found inside has stored nothing (no energy without the work that stored it; only ever down)
      // a union of terms (hinge_piece_sum) found inside: each term's own zero (piece id, its raw depth then), so force and book read one potential U = sum_t s_t V(a_t - z_t).
      // While it is not empty the terms carry the zeros and dZero is only their largest (the item's "has a zero" mark); each follows its term out (BackOut), dropped once clear
      std::vector<std::pair<int, double>> tZero;
      double TermZero(const int id) const
      {
         for (const std::pair<int, double> &z : tZero)
            if (z.first == id)
               return z.second;
         return 0.0;
      }
      // [Contact] joint_zero_floor: each joint term of h (HitTerm::over) takes a zero no smaller than the largest zero of the pieces it corrects, so its depth from its zero
      // is never above theirs (d_J <= min d_p: the joint never subtracts more than its pieces add, U >= 0). Raises a joint's zero only: a joint at or below its pieces' depths
      // stays where it was. The caller keeps dZero (the largest) in step
      void FloorJointZeros(const Hit &h)
      {
         for (int t = 0; t < h.nt; t++)
         {
            if (h.t[t].over == nullptr)
               continue;
            double f = 0.0;
            for (const int p : *h.t[t].over)
               f = std::max(f, TermZero(p));
            if (f <= 0.0)
               continue;
            size_t j = 0;
            while (j < tZero.size() && tZero[j].first != h.t[t].id)
               j++;
            if (j == tZero.size())
               tZero.emplace_back(h.t[t].id, 0.0);
            tZero[j].second = std::max(tZero[j].second, f);
            dZero = std::max(dZero, f);
         }
      }
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
      // [Contact] the elastic seam (multi-contact piece 1): the contact's stored energy at a place of the ball, V of its depth. The energy books read it here.
      // termSum ([Contact] term_energy): a contact made of terms reads the potential its force conserves, sum_t s_t V(max(d_t, 0)) over its zeroed terms (Elastic sums
      // s_t Chord(a0_t, a1_t) of the same V); a one-term contact (nt = 0) is V of its depth as before
      double Energy(const Hit &h, const bool termSum) const
      {
         if (!termSum || h.nt == 0)
            return V(h.delta);
         double e = 0.0;
         for (int t = 0; t < h.nt; t++)
            if (h.t[t].d > 0.0)
               e += h.t[t].s * V(h.t[t].d);
         return e;
      }
      // the stiffness slot k: the held kH (dk/dc = 0), read through here so a later piece can give it a law of its own. No reader yet. (An accessor, not a field:
      // kH is also written outside this file, contactsystem.cpp and cordfixed.h, and a second copy could drift from it.)
      double Stiff() const { return kH; }
      // what the elastic seam hands the step for one contact over one move
      struct ElasticEnd
      {
         bool terms = false;       // a union of hinge terms (hinge_piece_sum): the force is the terms' sum
         double fe = 0.0;          // the elastic normal force over the move, N: the chord (V(d1) - V(d0)) / (d1 - d0), or the terms' sum |Fel|
         double feT = 0.0;         // the terms' sum |Fel| (0 when not terms)
         double dd = 0.0;          // the depth change along gd: d1 - d0, or gd . dq for a terms contact
         double kc = 0.0;          // the chord's slope d fe / d d1 (K / 2 for a vanishing change): Newton's stiffness (only when asked for)
         double dE = 0.0;          // the elastic energy the move stores, fe dd (no reader yet)
         double Fel[4] = { 0.0, 0.0, 0.0, 0.0 };   // a terms contact's summed force s fe_t gd_t on the centre and the column
         bool cf0U = false;        // (VPX_GROUP_CLAMP_TRACE: the union's correction was left out)
         uint64_t tCf0 = 0;        // (VPX_GROUP_CLAMP_TRACE: the terms whose correction was left out)
      };
      // the elastic force of one contact over the move h0 -> hm -> h1 (depths d0, dm, d1): its mid gradient Gm, its discrete gradient gd over the move (Gonzalez on dq, iM),
      // the force fe along gd and, asked for, its slope kc. Defined after Group (it takes Group's Gonzalez). Scratch: the hinge terms' rows (Group's tgS ...).
      ElasticEnd Elastic(const Hit &h0, const Hit &h1, const Hit &hm, double d0, const double *dq, const double *iM, int na, int nq, double noiseD, const Settings &set, bool trace, bool wantKc,
                         double *gm, double *gd, std::vector<std::array<double, 4>> &tg, std::vector<double> &tk, std::vector<char> &tOn, std::vector<std::array<double, 4>> &tgE,
                         std::vector<std::array<double, 4>> &tgF) const;
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
      bool nearAll = false;      // true = every piece (an explicit move's step); false = the cycle's lists
      bool farOff = false;
      // the part finder: the cycle's lists per item (cycHas 0 = to build, 1 = built, 2 = not a surface of pieces; cycAS: the ball and shape it was built for), each ball's query
      // centre and reach, the cycle's length; nearCyc: this attempt's lists are the cycle's
      std::vector<std::vector<int>> cycL;
      std::vector<char> cycHas;
      std::vector<std::pair<int, int>> cycAS;
      std::vector<V3> cycC0;
      std::vector<double> cycReach;
      double cycT = 1e-3;
      double cycE = 0.0;   // the energy reach: the energy the cycle's part turns were given (FinderCycleStart)
      bool cycOn = false, nearCyc = false;
   };

   class Group : public PartRecords, public GroupUseWork
   {
      friend struct Item;   // the elastic seam (Item::Elastic) takes Gonzalez and kGradMax
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
         for (Buf<double> *b : { &d0s, &x, &resS, &stpS, &xPrevS, &stpPrevS, &trialS, &zeroS, &esRes, &xpS, &rpS, &rmS, &qEnd, &qMid, &qF, &fieldE0, &mlJii, &mlWork, &mlGm, &mlDq, &mlIm, &mlOut, &mlWg, &mlSnap, &mlDmis })
            WipeEachUse(b);
         WipeEachUse(&mlCt);
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
         nb = ns = ni = np = nf = nbf = nform = 0;
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
      int nb = 0, ns = 0, ni = 0, np = 0, nf = 0, nbf = 0;
      Buf<Ball> ball;
      AddedList<const BallField *> field;   // the position energies on the balls: one energy, each ball's difference over its move exact
      AddedList<const BallForce *> force;   // the place-dependent forces on the balls (BallForce): read at the midpoint, a drive
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
         paF.Grow(n);
         paM.Grow(n);
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
      void AddForce(const BallForce *const f)
      {
         force.b.Grow((size_t)nbf + 1);
         force.b[nbf++] = f;
      }
      // The ball's place-dependent forces at its centre c (per BallForce, times its forceScale); their slope added into slope when wanted
      V3 BallForceAt(const Ball &b, const V3 &c, M3 *const slope = nullptr) const
      {
         V3 f;
         for (int k = 0; k < nbf; k++)
         {
            M3 s;
            f = f + force[k]->At(c, b.R, slope != nullptr ? &s : nullptr) * b.forceScale;
            if (slope != nullptr)
               for (int i = 0; i < 3; i++)
                  for (int j = 0; j < 3; j++)
                     slope->a[i][j] += b.forceScale * s.a[i][j];
         }
         return f;
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
      // shared finder re-engage (runs/2026-11-13-finder-reengage/design.md): a surface the ball's inside check found inside an attempt joins the group in the
      // middle of a step. Its contacts are added after the others with nothing kept (as Configure makes them), every other contact keeps its state, the books go on; the owner sets each new
      // contact's zero (SetShapeItem) and the step's attempt runs again from the step's start state (Solve). The shape's storage is the owner's and must not move. Returns the first new item.
      int AddShapeMidStep(const Shape *const s)
      {
         AddShape(s);
         const int sIdx = ns - 1, k0 = ni;
         const size_t n = (size_t)ni + (s->onlyBall >= 0 ? 1 : (size_t)nb);
         const auto grow = [&](auto &b) {
            b.Grow(n);
            for (size_t k = (size_t)ni; k < n; k++)
               b[k] = std::decay_t<decltype(b[0])>();
         };
         grow(item); grow(law); grow(hNow); grow(h1c); grow(hmc); grow(d0s); grow(h0s); grow(ie); grow(farNow); grow(prevSlipS); grow(prevOnS); grow(szRec); grow(esI0);
         for (int a = 0; a < nb; a++)
         {
            if (s->onlyBall >= 0 && s->onlyBall != a)
               continue;
            farNow[ni] = 1;   // (a skipped contact until the next attempt's start classifies it: FarStillFar tests it)
            ie[ni].d1 = -1.0;   // (its end as a skipped contact's: Eval's far branch)
            Item &it = item[ni++];
            it = Item();
            it.a = a;
            it.s = sIdx;
         }
         stamp++;   // (every cached probe of the step's start is taken again)
         return k0;
      }
      // the owner's mid-step join, after each attempt: the number of shapes it added (> 0: the attempt runs again); allow false (past kJoinRounds): add none, -1 if one is waiting.
      // Unset: none (every use but the contact system's islands)
      std::function<int(bool)> midJoin;
      // the move of ball i in the current attempt (its centre's displacement over the step), m
      V3 MoveOf(const int i) const { return V3(x[6 * i], x[6 * i + 1], x[6 * i + 2]); }
      bool timeJoins = false;    // the owner's VPX_CHECK_PROFILE: time the attempts run again after a join (Solve)
      double joinRerunNs = 0.0;  // their time, ns, and how many (the owner reads and clears them)
      long joinReruns = 0;
      int joinRounds = 0;   // joins in the current step (Step)
      long joinCapHit = 0;  // attempts that found more after kJoinRounds joins in one step (the owner reads and clears it)
      static constexpr int kJoinRounds = 8;   // each round adds an object, so the rounds are bounded by the objects; more than this in one step is a fault, logged
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
            mlJii.Grow(nCoord);
            mlWork.Grow(np);
            mlGm.Grow(3 * (size_t)nb);
            mlDq.Grow(3 * (size_t)nb);
            mlIm.Grow(3 * (size_t)nb);
            mlOut.Grow(3 * (size_t)nb);
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
            mlWg.Grow((size_t)np * (size_t)jStride);
            mlSnap.Grow(7 * (size_t)jStride);
            mlDmis.Grow((size_t)jStride);
            paX.Grow(nU);   // (touch_predict piece 2: PartAccel's own move, rows and Jacobian, so the solve's are not touched)
            paRes.Grow(nU);
            paJ.Grow(nU * nU);
            paQ.Grow(nCoord);
            paAcc.Grow(nCoord);
            paF.Grow((size_t)nb);
            paM.Grow((size_t)nb);
            paDone.Grow((size_t)np);
            for (int p = 0; p < np; p++)
               paSave.Grow(part[p]->StateBytes());
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
            const double vn = it.touching ? it.Energy(hp, set.termEnergy) : 0.0, vs = it.touching ? it.ShearE(hp.delta) : 0.0;
            const double vr = it.touching && set.rollModel == 1 ? it.RollE(hp.delta) : 0.0;
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
            if (b.held)
               continue;   // (a held ball's motion and height are the kicker's: KickerDesign 1h)
            e += 0.5 * b.m * b.v.Dot(b.v) + 0.5 * b.I() * b.w.Dot(b.w) - b.m * set.g.Dot(b.c);
            if (nf > 0 && b.fieldScale != 0.0)
               e += FieldEnergy(b);
         }
         for (int k = 0; k < ni; k++)
            if (item[k].touching)
            {
               Hit hp;
               Probe(item[k], ball, hp);
               e += item[k].Energy(hp, set.termEnergy) + item[k].ShearE(hp.delta);
               if (set.rollModel == 1)
                  e += item[k].RollE(hp.delta);   // (rolling and twist as the step stored them: pmMulti / pmLev, pmEr)
            }
         for (int p = 0; p < np; p++)
            e += part[p]->Energy(ball, nb);
         e += FormEnergy(ball);
         return e;
      }
      // the energy reach: what part p's contacts hold now: each touching ball-part item's stored energy (the same terms as Energy) and its compliant stops' squash (J)
      double PartContactJ(const int p) const
      {
         double e = part[p]->StopSquashJ();
         for (int k = 0; k < ni; k++)
            if (item[k].touching && item[k].b < 0 && item[k].s >= 0 && shape[item[k].s] != nullptr && shape[item[k].s]->part == p)
               e += ItemHeldJ(k);
         return std::max(e, 0.0);
      }
      // the energy reach: what the balls' touching items on shapes s0 .. s1 - 1 hold now (a cord's posts and plate: CordBand::heldJ) (J)
      double ShapesContactJ(const int s0, const int s1) const
      {
         double e = 0.0;
         for (int k = 0; k < ni; k++)
            if (item[k].touching && item[k].b < 0 && item[k].s >= s0 && item[k].s < s1)
               e += ItemHeldJ(k);
         return std::max(e, 0.0);
      }
      // one item's stored energy now: its squash and its shear (and rolling and twist) springs, as Energy counts them
      double ItemHeldJ(const int k) const
      {
         Hit hp;
         Probe(item[k], ball, hp);
         double e = item[k].Energy(hp, set.termEnergy) + item[k].ShearE(hp.delta);
         if (set.rollModel == 1)
            e += item[k].RollE(hp.delta);   // (as Energy)
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
      // patch_zero (brief_code.md Z2, D2): a static mesh's or tube's item with its zero kept per patch: its probe returns every patch already shifted by its own zero
      bool PatchZeroItem(const Item &it) const { return set.patchZero && it.b < 0 && shape[it.s]->patchZero && !it.tZero.empty(); }
      bool islandRing = false;   // (diagnostics: this use is a ring island, else a cord island; set by the caller)
      const Hit *StartHit(const int k) const { return nowStamp == stamp && k >= 0 && k < ni ? &hNow[k] : nullptr; }   // patch_friction (ridge arrival): item k's probe at the positions now (ChooseStep's), null when stale
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
         struct CycOff { bool &f; explicit CycOff(bool &b) : f(b) {} ~CycOff() { f = false; } } cycOff(cycOn);   // part finder: the cycle's lists end with it
         FinderCycleStart(dt);
         SkipCycleStart();
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
      void SetShapeItem(const int k, const bool wasTouching, const double zero, const V3 &xi, const RubberContact::RollAngle &ra = RubberContact::RollAngle(),
                        const std::vector<std::pair<int, double>> *tZero = nullptr, const std::vector<RubberContact::PieceMem> *pm = nullptr, const bool pmMulti = false)
      {
         Item &it = item[k];
         Hit hp;
         it.dZero = 0.0;
         it.tZero.clear();
         if (!wasTouching && zero <= 0.0 && set.patchZero && shape[it.s]->patchZero)   // patch_zero (Z4): a zero per raw patch
            FirstSightPatches(it, hp);
         else if (!wasTouching && zero <= 0.0)
         {
            Probe(it, ball, hp);
            it.dZero = std::max(hp.delta, 0.0);
            if (it.dZero > 0.0 && hp.nt > 0)   // a union of terms: each term takes its own raw depth as its zero
            {
               it.dZero = 0.0;
               for (int t = 0; t < hp.nt; t++)
                  if (hp.t[t].d > 0.0)
                  {
                     it.tZero.emplace_back(hp.t[t].id, hp.t[t].d);
                     it.dZero = std::max(it.dZero, hp.t[t].d);
                  }
               if (set.jointZeroFloor)
                  it.FloorJointZeros(hp);
            }
         }
         else
         {
            it.dZero = zero;
            if (zero > 0.0 && tZero != nullptr)
               it.tZero = *tZero;
         }
         it.xi = xi;
         it.ra = ra;
         if (pm != nullptr)   // patch_friction: the piece memory kept between cycles (none: cleared)
            it.pm = *pm;
         else
            it.pm.clear();
         it.pmMulti = pm != nullptr && pmMulti;
         it.pmLev = false;   // (level_stiffness: the caller restores it with the piece memory)
         it.kStored = false;
         Probe(it, ball, hp);
         it.touching = hp.delta > 0.0;
         if (it.touching)
            TakeLaw(it, hp);
         else
         {
            it.xi = V3();
            it.ra = RubberContact::RollAngle();
            it.pm.clear();
            it.pmMulti = false;
            it.pmLev = false;
            it.kStored = false;
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
         int entry = -1;                     // >= 0: the contact's entry (Hit::en, step pieces) that set the step; rule is then its own + 10 (11 stiffness, 13 arrival, 17 steady clock); delta and dDot its
      };
      mutable StepDiag diag;

      // the depth an approach at closing speed v reaches on law L (all the closing energy into the squash): the stiffness while it lasts (ChooseStep; touch_predict's size d)
      double ReachDepth(const Item &L, const double v) const { return ReachDepth(L, v, L.mRed); }
      double ReachDepth(const Item &L, const double v, const double mRed) const   // (mRed: the mass the touch meets, a part's body with the ball's: ChooseStep's L.mRed)
      {
         const double eC = 0.5 * mRed * v * v;
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
      }

      // touch_predict (runs/2026-11-22-touch-predict/design.md 2.1, 2.3): a ball not touching a fixed surface (an item or a static piece entry), a moving part's surface (an item, or
      // a term of its union: 2.1b, 2.2) or another ball (2.1b, 2.6). ChooseStep writes one per such item, entry or term it predicted; the next Step reads them while the positions
      // are the same (tpStamp == stamp): nothing else is kept between steps
      enum TouchKind { tkFixed = 0, tkPart = 1, tkPair = 2, tkTerm = 3 };   // (a term counts as a part's surface: TouchStats)
      struct TouchPred
      {
         int k = -1, e = -1;      // the item and its entry (-1 the object); tkTerm: e is the term's index in the union
         int kind = tkFixed;
         double gap = 0.0, v = 0.0, a = 0.0, d = 0.0;   // at the step's start: gap (m), closing speed (m/s), closing acceleration (m/s^2), the law's squash d at the arrival speed (m)
         double t = 0.0, hToday = 0.0;   // the predicted touch time d/2 short (inf: none), and today's limit for the item (the net: max(0.5 gap / v, its clock))
         double mRed = 0.0;       // the mass the touch meets (ChooseStep's L.mRed): d at the end's closing speed
         V3 n;                    // its normal at the start (the end's closing speed is taken along it)
         int nc = 0;              // a part's surface: its columns at the start, each as the group's coordinate (coordOff + coord) and its point speed per unit rate along n
         int cc[7] = {};          // (a term: one column, the union's column 0, with the term's own slope)
         double cn[7] = {};
         bool past = false;       // (Step: touching at the end of the first solve)
      };
      mutable std::vector<TouchPred> tpList;
      mutable long tpStamp = -1;
      mutable int tpSet = -1;     // the prediction that set the step (-1 none) and its h
      mutable double tpH = 0.0;
      bool tpCounted = false;     // (Step: the step's prediction is counted once)

      // touch_predict piece 2 (design 2.2): a part's coordinates' accelerations at the state now, formed at most once per ChooseStep and only for a part whose surface an untouched
      // ball's item reaches within the longest step. One Newton correction of the part's own rows from the free move (each unknown at rate x h, the balls at v h with no spin change)
      // at h = kAccH: dx = -J^-1 R, R its rows from its Eval less h times the balls' loads on its coordinates now (each touching contact's start force, as fNet takes it, through
      // its columns' point speeds, as the group adds them), J its Jacobian block; a = 2 dx / h^2. The inertia row 2 m dq / h - 2 m v - h Q is -h Q on the free move: exact for a
      // constant mass, a pose-varying or linked inertia from the rows themselves (paper runs/2026-11-23-touch-p2paper P1-P2: <= 2.4e-6 relative at 1e-8 s, rounding ~1e-9 s).
      // A coordinate with no unknown (on a stop, a given move): 0 (P2 d); a massless one (no inertia row: its move is a balance, not an acceleration): 0.
      // The part's kept state is not written: its state records (PartRecords: state, lists, queue counts) are copied before and put back after, then AfterRestore re-makes what
      // derives from them (the bat's pivot, the cord's probe angle), as the net's retake does; its step work and attempt work are wiped at Step's start (WipeStepWork) and at each
      // attempt's start (BeginAttempt: WipeAttemptWork, then BeginStep makes its start values again) before the real solve reads them. Its Eval's diagnostics (last row terms)
      // are rewritten by every Eval of the real solve. Indexed by the part's coordinate; nullptr never (a failed Eval or a singular block: all 0)
      static constexpr double kAccH = 1e-8;   // (s) the free move's length (paper P1)
      mutable Buf<double> paX, paRes, paJ, paQ, paAcc;
      mutable Buf<V3> paF, paM;
      mutable Buf<unsigned char> paSave;
      mutable Buf<char> paDone;
      const double *PartAccel(const int p) const
      {
         double *const acc = (double *)paAcc + coordOff[p];
         if (paDone[p])
            return acc;
         paDone[p] = 1;
         Part &P = *part[p];
         const int nc = P.NumCoords();
         for (int c = 0; c < nc; c++)
            acc[c] = 0.0;
         int off = 6 * nb, nU = 6 * nb;   // (its unknowns' place and the unknowns' count, as NumUnk lays them out now)
         for (int q = 0; q < np; q++)
         {
            if (q == p)
               off = nU;
            nU += part[q]->NumDof();
         }
         const int nd = P.NumDof();
         if (nd == 0)
            return acc;
         touchStats.partEvals++;
         const double h = kAccH;
         const size_t sb = P.StateBytes();
         paSave.Grow(sb);   // (grows to the largest part's records once; no allocation after)
         P.SaveRecords(paSave);
         double *const X = paX, *const R = paRes, *const Jd = paJ, *const Q = paQ;
         V3 *const F = paF, *const M = paM;
         for (int i = 0; i < nb; i++)
         {
            const V3 g = ball[i].held ? ball[i].Given(h) : ball[i].v * h;
            X[6 * i] = g.x; X[6 * i + 1] = g.y; X[6 * i + 2] = g.z;
            X[6 * i + 3] = X[6 * i + 4] = X[6 * i + 5] = 0.0;
            F[i] = M[i] = V3();
         }
         for (int u = 6 * nb; u < nU; u++)
            X[u] = 0.0;
         for (int u = 0; u < nU; u++)
            R[u] = 0.0;
         std::fill(Jd, Jd + (size_t)nU * (size_t)nU, 0.0);
         // the balls' loads on its coordinates now: each touching contact on a surface it carries, its start force on the ball (fNet's) reversed, through each column's point speed
         V3 loadP;
         for (int c = 0; c < nc; c++)
            Q[c] = 0.0;
         for (int k = 0; k < ni; k++)
         {
            const Item &it = item[k];
            const Hit &hk = hNow[k];
            if (it.b >= 0 || !it.touching || hk.delta <= 0.0 || it.s < 0 || shape[it.s] == nullptr || shape[it.s]->part != p)
               continue;
            const V3 f = hk.n * it.F(hk.delta) - it.xi * it.Kt(hk.delta);
            loadP = loadP - f;
            for (int j = 0; j < hk.nc; j++)
               Q[hk.c[j].coord] -= f.Dot(hk.c[j].V());
         }
         P.guessAccel = set.partGuessAccel;
         P.gripInStep = set.gripInStep;
         P.WipeAttemptWork();
         bool ok = P.BeginStep(ball, nb, h, X, off);
         for (int c = 0; ok && c < nc; c++)   // the free move (BeginStep's guess may carry the last step's acceleration: part_guess_accel)
            if (const int d = P.DofOf(c); d >= 0)
               X[off + d] = P.Coord(c).rate * h;
         if (ok)
         {
            P.Loads(Q, loadP);
            ok = P.Eval(X, off, h, F, M, R, Jac(Jd, nU));
         }
         if (ok)
         {
            // its block: rows and columns off .. off + nd - 1, the loads' share of its rows (as the group's Eval: Res -= h qF), then Gaussian elimination with partial pivoting
            for (int c = 0; c < nc; c++)
               if (const int d = P.DofOf(c); d >= 0)
                  R[off + d] -= h * Q[c];
            const auto Jb = [&](const int r, const int c) -> double & { return Jd[(size_t)(off + r) * (size_t)nU + (size_t)(off + c)]; };
            double *const b = R + off;
            for (int r = 0; r < nd; r++)
               b[r] = -b[r];
            for (int col = 0; ok && col < nd; col++)
            {
               int piv = col;
               for (int r = col + 1; r < nd; r++)
                  if (std::fabs(Jb(r, col)) > std::fabs(Jb(piv, col)))
                     piv = r;
               if (!(std::fabs(Jb(piv, col)) > 0.0))
               {
                  ok = false;
                  break;
               }
               if (piv != col)
               {
                  for (int c = 0; c < nd; c++)
                     std::swap(Jb(piv, c), Jb(col, c));
                  std::swap(b[piv], b[col]);
               }
               for (int r = col + 1; r < nd; r++)
               {
                  const double m = Jb(r, col) / Jb(col, col);
                  for (int c = col; c < nd; c++)
                     Jb(r, c) -= m * Jb(col, c);
                  b[r] -= m * b[col];
               }
            }
            for (int r = nd - 1; ok && r >= 0; r--)
            {
               double s = b[r];
               for (int c = r + 1; c < nd; c++)
                  s -= Jb(r, c) * b[c];
               b[r] = s / Jb(r, r);
            }
            for (int c = 0; ok && c < nc; c++)
               if (const int d = P.DofOf(c); d >= 0 && P.Coord(c).invM != 0.0)
                  acc[c] = 2.0 * b[d] / (h * h);
         }
         P.RestoreRecords(paSave);
         P.AfterRestore();
         return acc;
      }

      // The step choice from the state now (accuracy only)
      double ChooseStep() const
      {
         PhaseScope ps(phChoose);
         double h = set.hMax;
         StepDiag dg;
         const bool tpOn = set.touchPredict && GroupFixedStepS() <= 0.0;   // (touch_predict; never under VPX_GROUP_HFIX)
         tpList.clear();
         tpSet = -1;
         for (int p = 0; tpOn && p < np; p++)   // (piece 2: no part's start acceleration formed yet)
            paDone[p] = 0;
         V3 *const fNet = fNetS;   // the net force on each ball now: gravity, each touching contact's spring and shear stretch; with touch_predict each part's own force (BallForceNow)
         for (int i = 0; i < nb; i++)
         {
            fNet[i] = (set.g * ball[i].m) + ball[i].fExt;
            if (nbf > 0 && ball[i].forceScale != 0.0)
               fNet[i] = fNet[i] + BallForceAt(ball[i], ball[i].c);
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
            Probe(it, ball, hp, true);   // (the step-start probe: a static surface's entries, step pieces)
            if (it.touching && hp.delta > 0.0)
            {
               const V3 f = hp.n * it.F(hp.delta) - it.xi * it.Kt(hp.delta);
               fNet[it.a] = fNet[it.a] + f;
               if (it.b >= 0)
                  fNet[it.b] = fNet[it.b] - f;
            }
         }
         if (tpOn)   // touch_predict: each part's own force on a ball now, not through an item (a cord's stretch and grip; design.md 2.7)
            for (int p = 0; p < np; p++)
               part[p]->BallForceNow(ball, nb, fNet);
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
            // one object, or one entry of it (step pieces: another convex piece or a ridge of a static surface, Hit::en): the same rule. delta, n, curvDist: its depth,
            // normal and turning distance; e: its entry, -1 the object (runs/2026-11-22-ringcap design 2.5)
            const auto one = [&](const double delta, const V3 &n, const double curvDist, const int e) {
               const V3 vr = it.b >= 0 ? A.v - ball[it.b].v : A.v - vs;
               const double dDot = -vr.Dot(n);   // closing speed
               const double vc = std::max(dDot, set.vFloor);
               // termArrival: each term's closing rate by its own gradient, -v . n_t + slope_t * rate (the part's one coordinate, column 0)
               const bool terms = e < 0 && set.termArrival && hp.nt > 0 && nq > 0;   // (an entry has no terms)
               const double qRate = terms ? part[shape[it.s]->part]->Coord(hp.c[0].coord).rate : 0.0;
               const auto termDot = [&](const HitTerm &ht) { return -A.v.Dot(ht.n) + ht.slope * qRate; };
               // the reach over the longest step along normal n: distance under a steady acceleration, s = v h + a h^2 / 2, a the closing part of the net acceleration
               // now (fNet: gravity with the nudge, external, fields, touching contacts, with touch_predict the parts' own forces; a surface's own acceleration is not in it). Only a closing a widens the reach,
               // an opening one is not taken off: the test only decides whether the item is probed; probing one more costs one probe, dropping one that arrives is a fault
               V3 aRelNow = fNet[it.a] * (1.0 / A.m);
               if (it.b >= 0)
                  aRelNow = aRelNow - (ball[it.b].held ? ball[it.b].heldA : fNet[it.b] * (1.0 / ball[it.b].m));
               const auto reach = [&](const double v, const V3 &n) { return v * set.hMax + 0.5 * std::max(-aRelNow.Dot(n), 0.0) * set.hMax * set.hMax; };
               bool termNear = false;   // a term the longest step can reach though the union's depth cannot be
               for (int t = 0; terms && !termNear && t < hp.nt; t++)
                  termNear = hp.t[t].d > 0.0 || -hp.t[t].d < reach(std::max(termDot(hp.t[t]), set.vFloor), hp.t[t].n) + 1e-9;
               if (delta <= 0.0 && -delta >= reach(vc, n) + 1e-9 && !termNear)
                  return;   // cannot be reached within the longest step
               Item L = it;
               if (!it.touching || e >= 0)   // (an entry: an untouched contact of the same material, its law from the item's probe)
                  TakeLaw(L, hp);
               if (invQ > 0.0)
                  L.mRed = 1.0 / (1.0 / A.m + invQ);
               // the depth an approach at closing speed v reaches (all the closing energy into the squash): the stiffness while it lasts
               const auto reachDepth = [&](const double v) { return ReachDepth(L, v); };
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
               const double dE = std::max(delta, lo);
               const double ktE = L.Kt(dE);
               const double omN = std::sqrt(L.K(dE) / L.mRed);
               // steady (Settings): by the motion, at the contact's own time scale. At the friction limit (Coulomb): a spring within mu N holds at most (mu N)^2 / (2 kt) of energy, so a slip with more than that
               // in 1/2 mt u^2 cannot be stopped by it: the contact slides, its force is mu N whatever the spring's stiffness, and the spring's time scale is not in the motion.
               bool steady = false;
               int why = 0;
               double whyVal = 0.0, whyLim = 0.0;
               if (set.vSteady > 0.0)
               {
                  if (delta > 0.0)
                  {
                     V3 u = A.v + A.w.Cross(n * (-(A.R - 0.5 * delta)));
                     if (it.b >= 0)
                        u = u - (ball[it.b].v + ball[it.b].w.Cross(n * (ball[it.b].R - 0.5 * delta)));
                     else
                        u = u - vs;   // (a fixed surface: 0)
                     const V3 ut = u - n * u.Dot(n);
                     const double slipNow = ut.Length();
                     V3 aRel = fNet[it.a] * (1.0 / A.m);
                     if (it.b >= 0)
                        aRel = aRel - (ball[it.b].held ? ball[it.b].heldA : fNet[it.b] * (1.0 / ball[it.b].m));
                     const double omNow = std::sqrt(L.K(delta) / L.mRed);
                     steady = std::fabs(dDot) <= set.vSteady && std::fabs(aRel.Dot(n)) <= set.vSteady * omNow && (slipNow <= set.vSteady || (set.slideFree && slipNow * std::sqrt(ktE * L.mt) > L.mu * L.F(dE)));
                     if (!steady)
                     {
                        if (std::fabs(dDot) > set.vSteady)
                           why = 1, whyVal = dDot, whyLim = set.vSteady;
                        else if (std::fabs(aRel.Dot(n)) > set.vSteady * omNow)
                           why = 2, whyVal = aRel.Dot(n), whyLim = set.vSteady * omNow;
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
               if (set.stepModel == 1 && steady && delta > 0.0)   // by the motion in its mode (Settings::stepModel)
               {
                  const double omR = std::sqrt(L.K(delta) / L.mRed);
                  V3 aRel = fNet[it.a] * (1.0 / A.m);
                  if (it.b >= 0)
                     aRel = aRel - (ball[it.b].held ? ball[it.b].heldA : fNet[it.b] * (1.0 / ball[it.b].m));
                  const double x = dDot / omR, y = aRel.Dot(n) / (omR * omR), amp = std::sqrt(x * x + y * y);
                  const double hm = amp > 0.0 ? std::sqrt(12.0 * set.epsRest / amp) / omR : set.hMax;
                  hs = std::max(hs, std::min(hm, set.hMax));
                  if (set.swingCap && 2.0 * amp <= set.epsRest)
                     hs = set.hMax;
               }
               const auto lim = [&](const double c, const int rule) {
                  if (c < h)
                  {
                     h = c;
                     dg.rule = e >= 0 ? rule + 10 : rule;
                     dg.item = k;
                     dg.entry = e;
                     dg.delta = delta;
                     dg.dDot = dDot;
                     dg.omN = omN;
                     dg.unsteady = why;
                     dg.unsteadyVal = whyVal;
                     dg.unsteadyLim = whyLim;
                     dg.term = -1;
                  }
               };
               // touch_predict: an item not touching a fixed surface gets no clock (2.1b, e): its closing speed v, closing acceleration a (the reach test's), arrival speed
               // va = sqrt(v^2 + 2 a gap) and d the law's squash at va; the step ends at the smallest root t of gap - d/2 = v t + a t^2 / 2, taken closing (the stable form
               // 2 s / (v + sqrt(v^2 + 2 a s)), or (sqrt(v^2 + 2 a s) - v) / a where that has no positive denominator: a receding item inside d/2 takes the later, closing
               // crossing; none: no limit). Not an event, today's rule instead: a closing item within d at the start (J5), an arrival slower than vSteady (J7).
               // Piece 2: the same rule on a moving part's surface and on a ball pair. a = -(aBall - aPoint) . n: aBall the reach test's (a pair: a's less b's, b's a held ball's
               // heldA), aPoint the surface point's acceleration, its columns' point speeds (the ones that give vs) times the part's coordinate accelerations now (PartAccel, 2.2);
               // the centripetal part of a turning point left out (2.2: 2.3's re-solve takes it)
               bool tpDone = false;
               const int tpKind = !tpOn || it.touching ? -1 : it.b >= 0 ? tkPair : it.s < 0 || shape[it.s] == nullptr ? -1 : shape[it.s]->part >= 0 ? tkPart : nq == 0 ? tkFixed : -1;
               const double *const qAcc = tpKind == tkPart && nq > 0 ? PartAccel(shape[it.s]->part) : nullptr;
               // the root rule for one touch (gap, closing speed v, closing acceleration a): false = not an event (J5, J7); d the law's squash at the arrival speed, t the root
               const auto predict = [&](const double gap, const double v, const double a, double &d, double &t) {
                  const double va2 = v * v + 2.0 * a * gap;
                  const bool arrives = va2 >= 0.0 && (v > 0.0 || a > 0.0);
                  const double va = arrives ? std::sqrt(va2) : 0.0;
                  d = reachDepth(va);
                  if (!((gap > d || v <= 0.0) && (!arrives || va >= set.vSteady)))
                     return false;
                  const double s = gap - 0.5 * d, disc = v * v + 2.0 * a * s;
                  t = std::numeric_limits<double>::infinity();
                  if (disc >= 0.0)
                  {
                     const double q = v + std::sqrt(disc);
                     const double tr = q > 0.0 ? 2.0 * s / q : a > 0.0 ? (std::sqrt(disc) - v) / a : -1.0;
                     if (tr > 0.0)
                        t = tr;
                  }
                  return true;
               };
               // its record (the columns of a part's surface for the end's closing speed), and the step ended at it when it is the earliest so far
               const auto record = [&](TouchPred &tp) {
                  tp.mRed = L.mRed;
                  tpList.push_back(tp);
                  if (tp.t < h)
                  {
                     lim(tp.t, 9);
                     tpSet = (int)tpList.size() - 1;
                     tpH = tp.t;
                     return true;
                  }
                  return false;
               };
               if (tpKind >= 0)
               {
                  V3 aPt;
                  for (int j = 0; qAcc != nullptr && j < nq; j++)
                     aPt = aPt + hp.c[j].V() * qAcc[hp.c[j].coord];
                  const double gap = -delta, a = qAcc != nullptr ? -aRelNow.Dot(n) + aPt.Dot(n) : -aRelNow.Dot(n);
                  double d = 0.0, t = 0.0;
                  if (predict(gap, dDot, a, d, t))
                  {
                     TouchPred tp;
                     tp.k = k;
                     tp.e = e;
                     tp.kind = tpKind;
                     tp.gap = gap, tp.v = dDot, tp.a = a, tp.d = d, tp.t = t;
                     tp.hToday = delta < 0.0 && dDot > 0.0 ? std::max(set.arrivalShare * -delta / dDot, hs) : hs;
                     tp.n = n;
                     for (int j = 0; tpKind == tkPart && j < nq; j++)
                     {
                        tp.cc[j] = coordOff[shape[it.s]->part] + hp.c[j].coord;
                        tp.cn[j] = hp.c[j].V().Dot(n);
                     }
                     tp.nc = tpKind == tkPart ? nq : 0;
                     record(tp);
                     tpDone = true;
                  }
               }
               if (!tpDone && delta < 0.0 && dDot > 0.0)   // reached in about two steps, not overshot
               {
                  const double ha = set.arrivalShare * -delta / dDot;
                  lim(std::max(ha, hs), ha > hs ? 3 : 1);
               }
               else if (!tpDone)
                  lim(hs, steady ? 7 : 1);
               // termArrival: each term not yet pressing is an approaching contact of its own: the arrival rule with the clock an untouched contact gets there
               // (unsteady, at the depth its own closing speed reaches), not the contact's (a resting contact's clock can be hMax)
               // touch_predict (piece 2): on an untouched item each term not pressing takes the root rule instead, by its own normal and slope: closing speed termDot, closing
               // acceleration -aBall . n_t + slope_t a_q (a_q the union's column-0 coordinate's, PartAccel); its record carries that one column with the term's slope
               for (int t = 0; terms && t < hp.nt; t++)
               {
                  const HitTerm &ht = hp.t[t];
                  const double dDotT = termDot(ht);
                  if (qAcc != nullptr && ht.d < 0.0)
                  {
                     const double aT = -aRelNow.Dot(ht.n) + ht.slope * qAcc[hp.c[0].coord];
                     double dT = 0.0, tT = 0.0;
                     if (predict(-ht.d, dDotT, aT, dT, tT))
                     {
                        double omNT = 0.0;
                        const double hsT = unsteadyClock(reachDepth(std::max(dDotT, set.vFloor)), omNT);
                        TouchPred tp;
                        tp.k = k;
                        tp.e = t;
                        tp.kind = tkTerm;
                        tp.gap = -ht.d, tp.v = dDotT, tp.a = aT, tp.d = dT, tp.t = tT;
                        tp.hToday = dDotT > 0.0 ? std::max(set.arrivalShare * -ht.d / dDotT, hsT) : hsT;
                        tp.n = ht.n;
                        tp.nc = 1;
                        tp.cc[0] = coordOff[shape[it.s]->part] + hp.c[0].coord;
                        tp.cn[0] = ht.slope;
                        if (record(tp))
                        {
                           dg.term = t;
                           dg.delta = ht.d;
                           dg.dDot = dDotT;
                        }
                        continue;
                     }
                  }
                  if (ht.d >= 0.0 || dDotT <= 0.0)
                     continue;
                  const double ha = set.arrivalShare * -ht.d / dDotT;
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
               if (curvDist > 0.0 && delta > -1e-4)   // rolling round a post or a ball
               {
                  const V3 vt = vr - n * vr.Dot(n);
                  const double sp = vt.Length();
                  if (sp > 1e-9)
                     lim(set.turnMaxRad * curvDist / sp, 4);
               }
            };
            one(hp.delta, hp.n, hp.curvDist, -1);
            for (int e = 0; e < hp.ne; e++)   // (an item whose own answer is out of reach is still read for its entries)
               one(hp.en[e].d, hp.en[e].n, 0.0, e);
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
         if (dg.rule != 5 && dg.rule != 8 && dg.entry >= 0)   // (step pieces counter: an entry set the step)
            pieceStats.set++;
         if (dg.rule != 9 && dg.rule != 19)   // (touch_predict: a later limit took the step)
            tpSet = -1;
         tpStamp = stamp;
         diag = dg;
         nowStamp = stamp;
         return std::max(h, 1e-7);
      }

      // touch_predict: the first s in (0, h) where the cubic Hermite of a gap's start and end values and rates (r = d gap / dt) reaches target (G(0) > target > G(h)):
      // its turning points split [0, h] into monotone pieces, the first piece that crosses is bisected to the last representable split. <= 0: none
      static double HermiteFirst(const double g0, const double r0, const double g1, const double r1, const double h, const double target)
      {
         const double c0 = g0 - target, c1 = h * r0, c2 = 3.0 * (g1 - g0) - h * (2.0 * r0 + r1), c3 = 2.0 * (g0 - g1) + h * (r0 + r1);   // on u = s / h
         const auto G = [&](const double u) { return c0 + u * (c1 + u * (c2 + u * c3)); };
         double br[4] = { 0.0, 1.0, 1.0, 1.0 };
         int nbr = 1;
         const auto turn = [&](const double u) {
            if (u > 0.0 && u < 1.0)
               br[nbr++] = u;
         };
         if (c3 != 0.0)
         {
            const double D = c2 * c2 - 3.0 * c1 * c3;
            if (D > 0.0)
            {
               turn((-c2 - std::sqrt(D)) / (3.0 * c3));
               turn((-c2 + std::sqrt(D)) / (3.0 * c3));
            }
         }
         else if (c2 != 0.0)
            turn(-c1 / (2.0 * c2));
         std::sort(br + 1, br + nbr);
         br[nbr++] = 1.0;
         for (int i = 0; i + 1 < nbr; i++)
         {
            double lo = br[i], hi = br[i + 1];
            if (!(G(lo) > 0.0 && G(hi) <= 0.0))
               continue;
            for (;;)
            {
               const double mid = 0.5 * (lo + hi);
               if (!(mid > lo && mid < hi))
                  break;
               (G(mid) > 0.0 ? lo : hi) = mid;
            }
            return lo * h;   // (the short side)
         }
         return -1.0;
      }

      // touch_predict, the step's end (design 2.3), after a solve of length h (the step was asked for hTry). stage 0: the first solve; 1: after the Hermite re-solve; 2: after the net.
      // A prediction past at the end (touching): one re-solve at the earliest Hermite root of the past ones, aimed d/2 short; still past: the net, a re-solve at today's limit for
      // that item (counted; WARN and notice: contactsystem.cpp). Not past, the prediction that set the step: taken (end gap in [0, d], d at the end's closing speed, and closing)
      // or short (end gap > d). Returns the h to solve again, 0 = none. The end's closing speed along the start normal; an entry's end gap the larger of the object's (its
      // nearest piece, from the solve) and its own face moved along its start normal. Piece 2 (design 2.3): a moving part's surface by its end probe at the part's solved end and
      // the relative end speed (its point's through the start columns), a ball pair by the pair's end depth and the two balls' end speeds; the same bands, re-solve and net
      double TouchEnd(const double h, const double hTry, int &stage)
      {
         // VPX_TOUCH_TRACE=<csv> (a diagnostic, read once; off: nothing changes): one row per predicted touch outcome
         static FILE *const tf = []() -> FILE * {
            const char *const e = getenv("VPX_TOUCH_TRACE");
            FILE *const o = e != nullptr && *e ? fopen(e, "w") : nullptr;
            if (o != nullptr)
               fprintf(o, "t_ms,object,entry,gap_um,v,a,d_um,tpred_us,h_us,endgap_um,endv,outcome,kind\n");
            return o;
         }();
         // (kind: fixed, part, term (entry = the term's index), pair (object: ball_pair))
         const auto row = [&](const TouchPred &p, const double gEnd, const double vEnd, const char *const out) {
            if (tf == nullptr)
               return;
            static const char *const kKind[4] = { "fixed", "part", "pair", "term" };
            const int s = item[p.k].s;
            fprintf(tf, "%u,%s,%d,%.4f,%.5f,%.5g,%.4f,%.4f,%.4f,%.4f,%.5f,%s,%s\n", phaseMs, p.kind == tkPair || s < 0 ? "ball_pair" : shape[s]->Name().c_str(), p.e, p.gap * 1e6, p.v,
                    p.a, p.d * 1e6, std::isfinite(p.t) ? p.t * 1e6 : -1.0, h * 1e6, gEnd * 1e6, vEnd, out, kKind[p.kind]);
         };
         if (stage >= 3)
            return 0.0;
         const bool counting = stage == 0 && !tpCounted && tpSet >= 0 && hTry == tpH;
         double hHerm = std::numeric_limits<double>::infinity(), hNet = hHerm, setG = 0.0, setV = 0.0, setD = 0.0;
         int kindHerm = tkFixed, kindNet = tkFixed;   // (the counters' kind: the touch that set the re-solve's length)
         const auto statKind = [](const int kind) { return kind == tkTerm ? (int)tkPart : kind; };
         bool anyPast = false;
         for (size_t j = 0; j < tpList.size(); j++)
         {
            TouchPred &p = tpList[j];
            if (p.k >= ni)
               continue;
            const Item &it = item[p.k];
            const V3 dc(x[6 * it.a], x[6 * it.a + 1], x[6 * it.a + 2]);
            double vEnd = -EndRate(dc, ball[it.a].v, h).Dot(p.n);
            // piece 2: a moving surface's point at the end, its start columns at each coordinate's end rate (the inertia row's trapezoid, 2 dq / h - rate, a given move's too;
            // a pose-varying inertia's p1 / I1 differs by its change over the step, a band test only); its own face moved (an entry's rule, for a term) takes them too.
            // A pair: the other ball's end rate
            double qMove = 0.0;   // (the surface's move along n over the step: a term's own end gap)
            for (int j = 0; j < p.nc; j++)
            {
               const int c = p.cc[j];
               const double dq = cUnk[c] >= 0 ? x[cUnk[c]] : ci[c].given;
               vEnd += p.cn[j] * EndRate(dq, ci[c].rate, h);
               qMove += p.cn[j] * dq;
            }
            if (p.kind == tkPair)
               vEnd += EndRate(V3(x[6 * it.b], x[6 * it.b + 1], x[6 * it.b + 2]), ball[it.b].v, h).Dot(p.n);
            const double d1 = farNow[p.k] ? -1.0 : ie[p.k].d1;
            // the end gap: the object's (fixed or pair: the solve's end depth; a part's surface: its end probe at the part's solved end), an entry's or a term's the larger of that
            // and its own face moved along its start normal
            const double gEnd = p.e < 0 ? -d1 : std::max(-d1, p.gap + dc.Dot(p.n) - qMove);
            if (gEnd < 0.0)   // touching at the end: past
            {
               anyPast = true;
               if (p.hToday < hNet)
               {
                  hNet = p.hToday;
                  kindNet = p.kind;
               }
               if (stage == 0)
               {
                  p.past = true;
                  const double s = HermiteFirst(p.gap, -p.v, gEnd, -vEnd, h, 0.5 * p.d);
                  if (s > 0.0 && s < hHerm)
                  {
                     hHerm = s;
                     kindHerm = p.kind;
                  }
               }
               row(p, gEnd, vEnd, stage == 0 ? "past" : stage == 1 ? "hermite_past" : "net_past");
               continue;
            }
            const double dEnd = farNow[p.k] ? p.d : ReachDepth(law[p.k], std::max(vEnd, 0.0), p.mRed);
            const char *const band = gEnd <= dEnd && vEnd > 0.0 ? "taken" : gEnd > dEnd ? "short" : "open";
            if (stage > 0 && p.past)
            {
               char out[24];
               snprintf(out, sizeof out, "%s_%s", stage == 1 ? "hermite" : "net", band);
               row(p, gEnd, vEnd, out);
            }
            if ((int)j == tpSet)
            {
               setG = gEnd;
               setV = vEnd;
               setD = dEnd;
            }
         }
         if (counting)
         {
            tpCounted = true;
            const int ks = statKind(tpList[(size_t)tpSet].kind);
            touchStats.predicted[ks]++;
            if (!anyPast)
            {
               const bool taken = setG <= setD && setV > 0.0, shortL = setG > setD;
               touchStats.taken[ks] += taken ? 1 : 0;
               touchStats.shortL[ks] += shortL ? 1 : 0;
               row(tpList[(size_t)tpSet], setG, setV, taken ? "taken" : shortL ? "short" : "open");
            }
         }
         if (stage == 2 || (!anyPast && stage == 1))
         {
            stage = 3;
            return 0.0;
         }
         if (!anyPast)
            return 0.0;
         if (stage == 0 && hHerm < h)
         {
            stage = 1;
            touchStats.hermite[statKind(kindHerm)]++;
            return hHerm;
         }
         stage = 2;   // (still past after the Hermite re-solve, or no root short of the end: the net)
         touchStats.net[statKind(kindNet)]++;
         return hNet < h ? hNet : 0.0;
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
         int tpStage = 0;   // (touch_predict: TouchEnd's stage)
         tpCounted = false;
         const bool tpOn = set.touchPredict && tpStamp == stamp && !tpList.empty() && GroupFixedStepS() <= 0.0;
         nHalveDiag = 0;
         joinRounds = 0;
         lwCuts = 0;
         pzMerged.clear();
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
               if (MergeBeforeAccept())   // patch_zero (F4, F7): the union joined start patches with different zeros: one zero for them at the start, the same h again
                  continue;
               if (set.linkWindows && LinkOrderCut(h))   // link_windows: the links' time order feeds an end patch less than the union: the step is too long
                  continue;
               if (tpOn && tpStage < 3)   // touch_predict: the step's end against its predicted touches (design 2.3)
               {
                  const double hRe = TouchEnd(h, hTry, tpStage);
                  if (hRe > 0.0)
                  {
                     h = hRe;
                     continue;
                  }
               }
               if (StopExcess() <= 0.0)
               {
                  Accept(h, false);
                  return h;
               }
               // the step that ends on the stop (Illinois false position on the step length), then solved again there (its end state is the one taken)
            stopSearch:   // (patch_zero, F8: a merge's re-solve at hBelow that passes a stop searches again from that solve)
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
               {
                  bool merged = false;
                  while (MergeBeforeAccept())   // patch_zero (F7): the stop search's final solve merged at its start: hBelow solved again (a patch pair merges once per step)
                  {
                     merged = true;
                     if (!Solve(hBelow))
                        return EnergyStep(hTry);
                  }
                  if (merged && StopExcess() > 0.0)   // patch_zero (F8): the re-solved step passes a stop: no Accept past it, the stop search again from this solve (h = hBelow)
                  {
                     h = hBelow;
                     goto stopSearch;
                  }
               }
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
               {   // the last round's residual: its three largest rows in speed units (as rn weighs them), and the contacts' states
                  const double *const R = resS;
                  for (int u = 0; u < nUnk; u++)
                  {
                     double v = 0.0;
                     if (u < 6 * nb)
                        v = R[u] * (u % 6 < 3 ? 1.0 / ball[u / 6].m : ball[u / 6].R / ball[u / 6].I());
                     else
                        for (int c = 0; c < nCoord; c++)
                           if (cUnk[c] == u)
                              v = Massless(c) ? (mlJii[c] != 0.0 ? R[u] / mlJii[c] * (2.0 / h) : 0.0) : R[u] * ci[c].invM * ci[c].reach;
                     for (int t = 0; t < 3; t++)
                        if (hd.topRow[t] < 0 || std::fabs(v) > std::fabs(hd.topVal[t]))
                        {
                           for (int m = 2; m > t; m--)
                           {
                              hd.topRow[m] = hd.topRow[m - 1];
                              hd.topVal[m] = hd.topVal[m - 1];
                           }
                           hd.topRow[t] = u;
                           hd.topVal[t] = (float)v;
                           break;
                        }
                  }
                  for (int k = 0; k < ni && hd.nCon < HalveDiag::kCon; k++)
                     if (ie[k].active)
                     {
                        HalveDiag::Con &c = hd.con[hd.nCon++];
                        c.k = k;
                        c.slide = ie[k].slipping;
                        c.ut = (float)std::sqrt(ie[k].ut2);
                        c.fRatio = ie[k].lim > 0.0 ? (float)(ie[k].ftTr / ie[k].lim) : -1.0f;
                        c.d1 = (float)ie[k].d1;
                     }
               }
               hd.nRound = std::min(lastRounds, kRounds);
               for (int r = 0; r < hd.nRound; r++)   // each round's move and residual (the halved line, AppendHalveNotes)
               {
                  hd.roundMv[r] = (float)roundDiag[r].mv;
                  hd.roundRn[r] = roundDiag[r].back ? -1.0f : (float)roundDiag[r].rn;
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
      // roundMv, roundRn: per round (the first kRounds) the round's move max(ch, grip) and the residual norm rn it started from (m/s); rn -1 marks a line-search retake
      static constexpr int kRounds = 64;
      // topRow, topVal: the last round's three largest residual rows (unknown index; speed units, m/s, as rn weighs them); con: the contacts in play at the last round
      // (item k, sliding at mu N or sticking, mid-step slip speed m/s, trial friction over mu N, end depth m)
      struct HalveDiag
      {
         char why = '?'; double h = 0.0; int rounds = 0; double moveC = 0.0, moveW = 0.0; int slipFlips = 0, touchFlips = 0, maxSliding = 0; double lastW[6] = {}; int nRound = 0; float roundMv[kRounds] = {}, roundRn[kRounds] = {};
         int topRow[3] = { -1, -1, -1 }; float topVal[3] = {};
         static constexpr int kCon = 6;
         struct Con { int k = 0; bool slide = false; float ut = 0.0f, fRatio = 0.0f, d1 = 0.0f; };
         int nCon = 0; Con con[kCon];
      };
      static constexpr int kHalveDiag = 8;
      HalveDiag halveDiag[kHalveDiag];
      // the step's failed attempts as the log lines print them (the energy-step line and the halved line, contactsystem.cpp and AdvanceAlone; one format): per attempt
      // [why, h, rounds, move, rim, slip and touch flips, sliding]; with rounds, each round's move/residual (m/s, "b" a line-search retake). Called only when a line is printed
      void AppendHalveNotes(std::string &msg, const bool rounds) const
      {
         char s[260];
         for (int a = 0; a < nHalveDiag; a++)
         {
            const HalveDiag &hd = halveDiag[a];
            snprintf(s, sizeof s, " [%c h %.3g us rounds %d move %.3g rim %.3g slipFlips %d touchFlips %d sliding %d", hd.why, hd.h * 1e6, hd.rounds, hd.moveC, hd.moveW, hd.slipFlips, hd.touchFlips, hd.maxSliding);
            msg += s;
            if (rounds)
            {
               msg += " | move/rn";
               for (int r = 0; r < hd.nRound; r++)
               {
                  if (hd.roundRn[r] < 0.0f)
                     snprintf(s, sizeof s, " b");
                  else
                     snprintf(s, sizeof s, " %.2g/%.2g", hd.roundMv[r], hd.roundRn[r]);
                  msg += s;
               }
               msg += " | last res";
               for (int t = 0; t < 3 && hd.topRow[t] >= 0; t++)
               {
                  const int u = hd.topRow[t];
                  std::string name = UnkName(u);
                  if (u >= 6 * nb)   // a part's coordinate: the part's kind
                     for (int p = np - 1; p >= 0; p--)
                        if (partOff[p] <= u)
                        {
                           name += std::string("(") + typeid(*part[p]).name() + ")";
                           break;
                        }
                  snprintf(s, sizeof s, " %s %.2g", name.c_str(), hd.topVal[t]);
                  msg += s;
               }
               msg += " | contacts";
               for (int c = 0; c < hd.nCon; c++)
               {
                  const HalveDiag::Con &cn = hd.con[c];
                  const Item &it = item[cn.k];
                  snprintf(s, sizeof s, " c%d b%d %s%d %s ut %.3g mm/s f/muN %.4g d1 %.2g um", cn.k, it.a, it.b >= 0 ? "ball" : "s", it.b >= 0 ? it.b : it.s, cn.slide ? "SLIDE" : "stick", cn.ut * 1e3, cn.fRatio, cn.d1 * 1e6);
                  msg += s;
               }
            }
            msg += "]";
         }
      }
      // the halved lines (one per step with a halved attempt, VPX_GROUP_ROUND_TRACE not needed) printed this table start, and those over the cap not printed (contactsystem.cpp: reset at
      // SetTable, the count of the rest at the table's end)
      static constexpr int kHalvedLines = 60;
      static int &HalvedLines() { static int n = 0; return n; }
      static long &HalvedNotPrinted() { static long n = 0; return n; }
      // Rule 8 diagnostic, always on, log only: after a step with a halved attempt, one INFO line (capped at kHalvedLines per table start); pre: the caller's own part (its names, its
      // balls) ahead of the attempts. Builds the string only for a line printed
      template <class Pre> void LogHalved(const Pre &pre) const
      {
         if (nHalveDiag == 0)
            return;
         if (HalvedLines() >= kHalvedLines)
         {
            HalvedNotPrinted()++;
            return;
         }
         HalvedLines()++;
         std::string msg = "Group step: halved,";
         pre(msg);
         msg += ", attempts";
         AppendHalveNotes(msg, true);
         for (int p = 0; p < np; p++)   // each part's own note of each failed attempt's last round (a bumper's skirt rows by term, a bat's rows)
            for (int a = 0; a < nHalveDiag; a++)
            {
               const std::string note = part[p]->FailNote(a);
               if (!note.empty())
                  msg += " | part " + std::to_string(p) + " (" + typeid(*part[p]).name() + ") attempt " + std::to_string(a) + ": " + note;
            }
         if (HalvedLines() == kHalvedLines)
            msg += " (no more of these lines this table start: the count of the rest comes at its end)";
         PLOGI << msg;
      }
      std::string roundBuf, failRoundBuf;   // VPX_GROUP_ROUND_TRACE: this attempt's rounds; the first failed attempt of the step
      mutable std::string probeTrace;       // VPX_GROUP_ROUND_TRACE: the carried surfaces' probes of the last Eval (end and mid, their pieces)
      // VPX_GROUP_JAC_WINDOW=<t0>-<t1> (ms of game time; with VPX_GROUP_JAC_MODE=1 or 2; a diagnostic, read once; unset = off): the difference Jacobian is formed (mode 2: and solved with) only in
      // solves inside that window (outside it a solve runs exactly as with VPX_GROUP_JAC_MODE=0); for solves that converge in that window, each round's
      // residual and evaluations, and every Jacobian entry (balls and parts) where the solve's own (a) and central differences (f) differ by more than 10 % and 1e-6
      // (jacWinBuf, written by the owner; 600 lines in all)
      static const double *JacWindowRange()
      {
         static const double *const w = []() -> const double * {
            static double v[2];
            const char *const e = getenv("VPX_GROUP_JAC_WINDOW");
            return e != nullptr && sscanf(e, "%lf-%lf", &v[0], &v[1]) == 2 ? v : nullptr;
         }();
         return w;
      }
      static bool JacWindow(const uint32_t ms)
      {
         const double *const w = JacWindowRange();
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
      unsigned pfMs = 0;          // patch_friction: the table time at the cycle's start (ms), for its counter's log
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
            if (b.held)
               continue;
            s.ke += 0.5 * b.m * b.v.Dot(b.v) + 0.5 * b.I() * b.w.Dot(b.w);
            s.pe += -b.m * set.g.Dot(b.c) + (nf > 0 && b.fieldScale != 0.0 ? FieldEnergy(b) : 0.0);
         }
         for (int k = 0; k < ni; k++)
            if (item[k].touching)
            {
               Hit hp;
               Probe(item[k], ball, hp);
               s.con += item[k].Energy(hp, set.termEnergy) + item[k].ShearE(hp.delta);
               if (set.rollModel == 1)
                  s.con += item[k].RollE(hp.delta);   // (as Energy)
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
            const double vn = it.touching ? it.Energy(hp, set.termEnergy) : 0.0, vs = it.touching ? it.ShearE(hp.delta) : 0.0;
            const double vr = it.touching && set.rollModel == 1 ? it.RollE(hp.delta) : 0.0;
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
      // mv: the round's move max(ch, grip), as the no-progress rule reads it; rn: the residual norm the round started from (m/s)
      struct RoundDiag { double moveC = 0.0, moveW = 0.0, scale = 1.0, grip = 0.0, mv = 0.0, rn = 0.0; int nActive = 0, nSlip = 0, nSlipFlip = 0, nOnFlip = 0, back = 0; };
      RoundDiag roundDiag[kRounds];
      // Settings::jacMode: the solve's own Jacobian and the central-difference one at the last round of the last solve (rows: per ball the momentum residual 3,
      // then the spin residual 3; columns: per ball the move 3, then the spin change 3)
      Jac JacAna() { return Jac(jacAnaS, jStride); }
      Jac JacFd() { return Jac(jacFdS, jStride); }
      int jacRound = 0;
      int jacNow = 0;   // this solve's Settings::jacMode (with VPX_GROUP_JAC_WINDOW: the mode inside the window, 0 outside)

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
         double pmEs = 0.0, pmEr = 0.0;   // patch_friction, 2 or more groups: the springs' end energy summed over the groups (Item::pmEs)
         bool kStored = false;            // plan 2.4: one friction group: its end stiffness kt1, kr1, kw1 (Item::kStored)
         double kt1 = 0.0, kr1 = 0.0, kw1 = 0.0, kReff = 0.0;
         double memLossS = 0.0, memLossR = 0.0;   // patch_friction: the memory rule's loss booked in e.slip, e.roll (shear; rolling and twist), J
         bool pmLev = false;              // level_stiffness, one group on the level form: pmEs, pmEr are its springs' end energy at the level form's end stiffness (Item::pmLev)
         double ut2 = 0.0, ftTr = 0.0, lim = 0.0;   // diagnostics, never read by the step: the mid-step slip speed squared (m2/s2), the trial friction force and mu N (N)
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
      // patch_friction (piece 2b-1): one friction group of a step (Eval): its inputs (normal at the middle n, end gradient direction nD, depths, its load fe + fd = N, aE,
      // its chord slope kc and normal stiffness kn, the end gradient gm, its memory xi, ra) and the friction's results the slope reads
      struct FG
      {
         V3 n, nD, xi;
         RubberContact::RollAngle ra;
         double d0 = 0.0, d1 = 0.0, dm = 0.0, fe = 0.0, fd = 0.0, N = 0.0, aE = 0.0, dadfE = 0.0, kc = 0.0, kn = 0.0, gm[3] = {};
         bool touching = false;
         int trCmn = -1, trFound = -1;   // (VPX_PATCH_TRACE: hm patches sharing a piece with the group (0: its reference plane); GroupMemory's start patches found in it.pm (-1 not called, -2 not touching))
         // results
         V3 Ft, xi0, xi1, ut, FsTr, rA, rB, wr, wtR, mhR;
         M3 RmS;
         RubberContact::RollAngle ra1;
         double kS = 0.0, kSp = 0.0, ctSp = 0.0, k0 = 0.0, k1 = 0.0, ct = 0.0, limS = 0.0, flTr = 0.0, kTw = 0.0, kkTw = 0.0, sRo = 0.0, rhoR = 0.0, wlR = 0.0;
         double twD1 = 0.0, twA = 0.0, twN = 0.0, twF = 0.0, limN = 0.0, limF = 0.0, twM = 0.0, twLim = 0.0;
         int twSt = 0;
         bool slipping = false, rollS = false;
         bool blin = false;               // level_stiffness: the group's stiffness from B-lin (blS): its slopes in the end position are the vectors below, kSp, ctSp, twD1 0
         V3 kSpV, ctSpV, twD1V;           // dkS/dc1, dct/dc1, dMw/dc1 (sticking twist) (B-lin: per piece, MeanSlope times its end gradient -n_i)
         double kr1L = 0.0, kw1L = 0.0;   // B-lin: the rolling and twist stiffness at the end (the level form's; k1 the shear's)
         double ks = 0.0, krs = 0.0, kws = 0.0;      // plan 2.4: the stiffness its memory was stored at (shared at a split, summed at a join), 0 none
         double k0p = 0.0, kr0p = 0.0, kw0p = 0.0;   // (its start stiffness, as the block takes it: for the split's shares)
         double heldS = 0.0, heldR = 0.0;   // (VPX_PATCH_CHECK: the shear and rolling + twist springs' held energy at the step's start, as the block takes it, J)
      };
      mutable std::vector<long> chkAtt;   // (VPX_PATCH_CHECK: per item, the attempt last compared held vs book)
      std::vector<FG> fgS;                                   // (Eval's scratch, capacity reused)
      std::vector<int> fgEnd;                                // each end patch's group (-1 not pressed)
      std::vector<std::vector<RubberContact::PieceMem>> pmEnd;   // per item: the last Eval's end pieces with their group's memory (Accept writes it to the item)
      std::vector<signed char> pmEndOn;                      // per item: 0 not this path, 1 one group, 2 two or more
      // diagnostics VPX_PATCH_TRACE=t0-t1 (engine ms, read once): inside the window, one log line per Eval of a patch_friction contact (PTEV: each friction group's start and end
      // pieces, depths, load, stick/slip, memory found, stretch) and one per start patch with none of its pieces in it.pm (PTMISS); 400 lines in all. Outside: a bool test only
      static bool PatchTraceWin(const uint32_t ms)
      {
         static const double *const w = []() -> const double * {
            static double v[2];
            const char *const e = getenv("VPX_PATCH_TRACE");
            return e != nullptr && sscanf(e, "%lf-%lf", &v[0], &v[1]) == 2 ? v : nullptr;
         }();
         return w != nullptr && ms >= w[0] && ms <= w[1];
      }
      mutable int trLines = 0;               // (VPX_PATCH_TRACE lines written)
      int trRound = 0, trFd = 0;             // (VPX_PATCH_TRACE: the solve's round; 1 inside NumericJacobian)
      long trAtt = 0;                        // (VPX_PATCH_TRACE: attempts started)
      static std::string TraceIds(const Hit &hh, const HitPatch &P) { std::string s; char b[24]; for (int q = 0; q < P.pn; q++) { snprintf(b, sizeof b, "%s%d", q ? "," : "", hh.pc[P.p0 + q]); s += b; } return s; }
      static std::string TracePmIds(const std::vector<RubberContact::PieceMem> &v) { std::string s; char b[24]; for (size_t q = 0; q < v.size(); q++) { snprintf(b, sizeof b, "%s%d", q ? "," : "", v[q].id); s += b; } return s; }
      std::vector<signed char> pmEndJoin, pmEndSplit;        // per item: the last Eval's groups with more start than end patches (join), more end than start (split)
      int fgJoins = 0, fgSplits = 0;                         // FrictionGroups' last call: its joins and splits
      // patch_friction counters (running totals logged at 1, 2, 4, ...; the bench never closes the table): Evals with 2+ groups, joins, splits, memory pieces dropped unpressed
      static void PfCount(std::atomic<long long> &n, const long long add, const char *what, const unsigned ms)
      {
         if (add <= 0)
            return;
         const long long was = n.fetch_add(add), c = was + add;
         if (was == 0 || (c ^ was) > was)   // a power of 2 crossed (the top bit moved)
            PLOGI << "patch_friction: " << c << " " << what << " (table time " << ms << " ms)";
      }
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
            warmA[i] = evF[i] * ball[i].InvM();
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
      // a massless coordinate (CoordInfo invM 0 with an unknown; runs/2026-11-19-standup-p2/design.md section 3): any this attempt (CoordsNow); each one's row slope at the last
      // Eval with the Jacobian (the line search's units, CoordRes2); per part the work the contacts' step-mean reaction does on it in an Eval (the energy correction); scratch
      bool anyMassless = false;
      Buf<double> mlJii, mlWork, mlGm, mlDq, mlIm, mlOut;
      // with the Jacobian: per part the slope of mlWork over every unknown (row p, stride jStride); each massless column's row as the contact's shear and spin terms
      // left it (mlSnap, 7 per contact: a Hit's columns); the mismatch's slope (mlDmis)
      Buf<double> mlWg, mlSnap, mlDmis;
      // the first guess (MasslessGuess): each contact a massless column reaches in the step, its coordinate c, its item k, its depth at the coordinate's move 0 (a: the start
      // depth d0 plus the balls' guessed moves and every other unknown's) and its slope s along the move
      struct MlContact
      {
         int c = 0, k = 0;
         double a = 0.0, s = 0.0, d0 = 0.0;
      };
      Buf<MlContact> mlCt;
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
      // [Contact] hinge_piece_sum: per contact, the pieces of its surface the ball can reach in this cycle (its cycle's list, Shape::CycleNear); every probe
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
      // the lists of this attempt: each contact's cycle list (CycleList). The far test's second solve (farOff), an explicit move (nearAll) and the
      // diagnostic VPX_HINGE_NEAR_ALL take every piece
      void BuildNear(const double h)
      {
         nearOn = false;
         if (nearL.size() < (size_t)ni)
            nearL.resize((size_t)ni);
         const bool all = farOff || NearAllOn() || nearAll;
         // the part finder: the cycle's lists (CycleList); a step outside a cycle's StepsOver, or one that takes every piece, takes all (never the per-step lists)
         nearCyc = !all;
         for (int k = 0; k < ni; k++)
         {
            Item &it = item[k];
            it.nearK = -1;
            if (it.b >= 0 || shape[it.s]->part < 0)
               continue;
            if (nearCyc && cycOn && CycleList(k))
            {
               nearL[(size_t)k] = cycL[(size_t)k];
               it.nearK = k;
            }
         }
         nearOn = true;
      }
      // the part finder: item k's cycle list, built once per cycle (and again after its ball's re-query) from the ball's centre at the cycle's start and its reach; false = not a
      // surface of pieces
      bool CycleList(const int k)
      {
         const Item &it = item[k];
         if (cycHas.size() <= (size_t)k)
         {
            cycHas.resize((size_t)k + 1, 0);
            cycL.resize((size_t)k + 1);
            cycAS.resize((size_t)k + 1, std::make_pair(-1, -1));
         }
         if (cycHas[(size_t)k] != 0 && cycAS[(size_t)k] == std::make_pair(it.a, it.s))
            return cycHas[(size_t)k] == 1;
         if (cycC0.size() <= (size_t)it.a)   // (a ball the cycle's start did not have: its query from where it is now)
         {
            const size_t n0 = cycC0.size();
            cycC0.resize((size_t)it.a + 1);
            cycReach.resize((size_t)it.a + 1);
            for (size_t i = n0; i < cycC0.size(); i++)
            {
               cycC0[i] = ball[i].c;
               cycReach[i] = ContactReachM(ball[i].v.Length(), 0.0, cycT);
            }
         }
         const Shape *const sh = shape[it.s];
         int pieces = 0;
         const bool ok = sh->CycleNear(cycC0[(size_t)it.a], ball[it.a].R, cycReach[(size_t)it.a], cycL[(size_t)k], pieces);
         cycHas[(size_t)k] = ok ? 1 : 2;
         cycAS[(size_t)k] = std::make_pair(it.a, it.s);
         if (ok && partFinderProf)
         {
            PartFinderCount &pc = partFinderCount[typeid(*part[sh->part]).name()];
            pc.lists++;
            pc.listLen += (long long)cycL[(size_t)k].size();
            pc.pieces += pieces;
         }
         return ok;
      }
      // the part finder, the cycle's start (StepsOver): each part's turn over the cycle from the island's balls' kinetic energy (spin included), each ball's query volume
      // (its centre now, its reach ContactReachM from its speed now; its acceleration is not known here: 0, the inside check covers it)
      void FinderCycleStart(const double T)
      {
         cycOn = true;
         cycT = T;
         double E = 0.0;
         for (int i = 0; i < nb; i++)
            E += 0.5 * ball[i].m * ball[i].v.Dot(ball[i].v) + 0.5 * ball[i].I() * ball[i].w.Dot(ball[i].w);
         for (int i = 0; i < nb; i++)   // the island's feed energy from the grouping (its balls' kinetic and stored energy and its parts' drive work), never less than its balls' energy now
            E = std::max(E, ball[i].feedE);
         cycE = E;
         for (int p = 0; p < np; p++)
            part[p]->CycleTurn(E, T);
         cycC0.resize((size_t)nb);
         cycReach.resize((size_t)nb);
         for (int i = 0; i < nb; i++)
         {
            cycC0[(size_t)i] = ball[i].c;
            cycReach[(size_t)i] = ContactReachM(ball[i].v.Length(), 0.0, T);
            cycReach[(size_t)i] = std::max(cycReach[(size_t)i], ball[i].reachS);   // (the grouping's travel: the feed energy given all to this ball, its forces as the drive)
         }
         std::fill(cycHas.begin(), cycHas.end(), (char)0);
      }
      // the part finder, after the solve: the inside check of every ball with a listed contact at the move's middle and end (its centre within its reach of where the cycle's
      // query was made). Outside: the query again from the ball's solved place, its reach ContactReachM over the cycle's length T (at least the time left from this step's start, this
      // step included, so its own middle is inside; the group does not keep the time done) from its solved speed and the acceleration over the step just solved (the end rate 2 dc / h - v); its lists rebuilt, the attempt again.
      // Then the parts' range check (group A): a part whose solved move left the range its lists were sized for takes a new range, every list on it rebuilt, the attempt again
      bool CycleStillInside(const double *xv, const double h)
      {
         bool inside = true;
         for (int k = 0; k < ni; k++)
         {
            const Item &it = item[k];
            if (it.nearK < 0)
               continue;
            const Ball &A = ball[it.a];
            const V3 dc(xv[6 * it.a], xv[6 * it.a + 1], xv[6 * it.a + 2]);
            const V3 &c0 = cycC0[(size_t)it.a];
            const double d = std::max((A.c + dc * 0.5 - c0).Length(), (A.c + dc - c0).Length());
            if (d <= cycReach[(size_t)it.a])
               continue;
            inside = false;
            const V3 v1 = EndRate(dc, A.v, h);
            cycC0[(size_t)it.a] = A.c + dc;
            cycReach[(size_t)it.a] = ContactReachM(v1.Length(), (v1 - A.v).Length() / h, cycT);
            for (int j = 0; j < ni; j++)
               if (item[j].a == it.a && (size_t)j < cycHas.size())
                  cycHas[(size_t)j] = 0;
            if (partFinderProf)
               partFinderCount[typeid(*part[shape[it.s]->part]).name()].requeries++;
         }
         double E = -1.0;
         for (int k = 0; k < ni; k++)
         {
            const Item &it = item[k];
            if (it.nearK < 0)
               continue;
            const int p = shape[it.s]->part;
            const double *const qe = (const double *)qEnd + coordOff[p];
            if (part[p]->CycleRangeHolds(qe))
               continue;
            inside = false;
            if (E < 0.0)   // the island's balls' kinetic energy now (spin included), as at the cycle's start
            {
               E = 0.0;
               for (int i = 0; i < nb; i++)
                  E += 0.5 * ball[i].m * ball[i].v.Dot(ball[i].v) + 0.5 * ball[i].I() * ball[i].w.Dot(ball[i].w);
               E = std::max(E, cycE);   // (and never less than the cycle's feed energy)
            }
            part[p]->CycleRegrow(qe, E, cycT);
            for (int j = 0; j < ni; j++)
               if (item[j].b < 0 && shape[item[j].s]->part == p && (size_t)j < cycHas.size())
                  cycHas[(size_t)j] = 0;
            if (partFinderProf)
               partFinderCount[typeid(*part[p]).name()].rangeRuns++;
         }
         if (!inside || !FinderCheckEnv())
            return inside;
         // check mode: every piece within R of its ball at the move's start, middle and end (the rigid piece at the trial angle) is on the list
         for (int k = 0; k < ni; k++)
         {
            const Item &it = item[k];
            if (it.nearK < 0)
               continue;
            const Shape *const sh = shape[it.s];
            const Ball &A = ball[it.a];
            const V3 dc(xv[6 * it.a], xv[6 * it.a + 1], xv[6 * it.a + 2]);
            const double *const qe = (const double *)qEnd + coordOff[sh->part], *const qm = (const double *)qMid + coordOff[sh->part];
            const std::vector<int> &l = cycL[(size_t)k];
            const int m = sh->CycleMisses(A.c, A.R, 0.0, l) + sh->CycleMisses(A.c + dc * 0.5, A.R, qm[0], l) + sh->CycleMisses(A.c + dc, A.R, qe[0], l);
            PartFinderCount &pc = partFinderCount[typeid(*part[sh->part]).name()];
            pc.checks += 3;
            pc.misses += m;
            if (m > 0 && partFinderFirstMiss.empty())
            {
               char b[300];
               snprintf(b, sizeof b, "%s %s: %d piece(s) within R off the cycle list, ball at (%.5f %.5f %.5f) m, move (%.3g %.3g %.3g) m, cycle start (%.5f %.5f %.5f) reach %.4g m", typeid(*part[sh->part]).name(),
                        sh->Name().c_str(), m, A.c.x, A.c.y, A.c.z, dc.x, dc.y, dc.z, cycC0[(size_t)it.a].x, cycC0[(size_t)it.a].y, cycC0[(size_t)it.a].z, cycReach[(size_t)it.a]);
               partFinderFirstMiss = b;
            }
         }
         return inside;
      }
      // after the solve: every piece off a list is off its ball at the move's middle and end (the start: by its building)
      bool NearStillNear(const double *xv, const double h)
      {
         if (carried)
            SetTrials(xv);
         return !nearCyc || CycleStillInside(xv, h);
      }
      Buf<double> fieldE0;       // what each ball's place stores at the start of the attempt
      // out of reach: the gap is more than both bodies' speeds can close in twice the step. A body that picks up speed inside the step is what the check after the solve is for
      // energy skip test (runs/2026-11-15-finder-P6b2/paper.md candidate c), at the cycle's start from the grouping's feed: each part's top point speed (m/s) over the cycle,
      // sum over its coordinates of reach x sqrt(2 E invM) (1/2 I w^2 <= E), with E its own (what it holds before a ball touches it: one ball in the group, the item's own) or the most any
      // ball's group can hold (more balls; any held ball: a given motion, no energy bound). A coordinate with no mass (no unknown and none free) or no energy given: no bound
      std::vector<double> skPOwn, skPAll;
      int skNb = -1;
      void SkipCycleStart()
      {
         skNb = nb;
         double E = 0.0;
         bool held = false, given = true;
         for (int i = 0; i < nb; i++)
         {
            E = std::max(E, ball[i].eTop);
            held = held || ball[i].held;
            given = given && ball[i].vTop >= 0.0;
         }
         if (!given)
            skNb = -1;   // (a ball with no feed: skip nothing this cycle)
         skPOwn.assign((size_t)np, INFINITY);
         skPAll.assign((size_t)np, INFINITY);
         for (int p = 0; p < np; p++)
         {
            double so = 0.0, sa = 0.0;
            for (int c = 0; c < part[p]->NumCoords(); c++)
            {
               const CoordInfo q = part[p]->Coord(c);
               const double im = q.invM > 0.0 ? q.invM : q.invMFree;
               if (im <= 0.0 || q.given != 0.0)
               {
                  so = sa = INFINITY;
                  break;
               }
               so += q.reach * std::sqrt(2.0 * std::max(part[p]->eOwnJ, 0.0) * im);
               sa += q.reach * std::sqrt(2.0 * E * im);
            }
            if (part[p]->eOwnJ >= 0.0)
               skPOwn[(size_t)p] = so;
            if (!held)
               skPAll[(size_t)p] = sa;
         }
      }
      // energy skip test: ball A's own bound on its speed over the step when all its ball-pair and moving-part items are skipped (paper_stage3.md section 1: before any of them
      // touches, only its static contacts, gravity, the nudge and magnets act on it, so 1/2 m v^2 <= its kinetic energy with spin + what its static contacts hold + the forces' work):
      // sqrt(v^2 + I w^2 / m + 2 E_static / m) + a h, a = aOwn + |fExt| / m; a held ball: its given speed |v| + |heldA| h. Taken only past the cycle's speed (both are bounds: the least holds)
      std::vector<char> skOwn;     // per ball: own in this attempt (set in the trip through the items)
      std::vector<int> skDef;      // the attempt's deferred items: static, not touching, not skipped by the cycle's speed (index only; own test and law after the trip)
      std::vector<Hit> skDefH;     // their start probes, for the law (as the trip took them)
      std::vector<double> skEs;    // per ball: what its touching static items hold at the attempt's start, J (-1: not summed yet; summed on the first own test that needs it)
      bool OwnFar(const Item &it, const double d0, const double h)
      {
         const Ball &A = ball[it.a];
         const double g = -d0;
         if (g <= 0.0 || A.aOwn < 0.0)
            return false;
         skipOwnTested++;
         if (!A.held && skEs[(size_t)it.a] < 0.0)   // (a touching static item: what it holds now, as ItemHeldJ from this start's depth d0s; the trip probed them all)
         {
            skipEsComputed++;
            double es = 0.0;
            for (int k = 0; k < ni; k++)
            {
               const Item &q = item[k];
               if (q.a != it.a || q.b >= 0 || q.s < 0 || shape[q.s] == nullptr || shape[q.s]->part >= 0 || !q.touching)
                  continue;
               const double d = d0s[k];
               double e = q.V(d) + 0.5 * q.Kt(d) * q.xi.Dot(q.xi);
               if (set.rollModel == 1)
                  e += q.RollE(d);   // (as Energy)
               es += std::max(e, 0.0);
            }
            skEs[(size_t)it.a] = es;
         }
         const double vo = A.held ? A.v.Length() + A.heldA.Length() * h
                                  : std::sqrt(A.v.Dot(A.v) + A.I() / A.m * A.w.Dot(A.w) + 2.0 * skEs[(size_t)it.a] / A.m) + (A.aOwn + A.fExt.Length() / A.m) * h;
         const bool beyond = g > vo * h;
         if (beyond)
            skipOwnSkipped++;
         return beyond;
      }
      bool OutOfReach(const Item &it, const double d0, const double h, const Hit &h0) const
      {
         if (it.touching || d0 > 0.0)
            return false;
         // gap > what both sides can close in h by the energy law (a gap equal to it is in reach)
         if (skNb != nb || it.a >= skNb || (it.b >= skNb))
            return false;
         const Ball &A = ball[it.a];
         const double g = -d0;
         if (g <= 0.0)
            return false;
         if (it.b >= 0)
         {
            const Ball &B = ball[it.b];
            if (!A.held && !B.held && A.eGroup == B.eGroup)   // (Cauchy-Schwarz: (|va| + |vb|)^2 <= 2 Etot (1/ma + 1/mb), squared: no sqrt)
               return g * g > 2.0 * std::max(A.eTop, B.eTop) * (1.0 / A.m + 1.0 / B.m) * h * h;
            return g > (A.vTop + B.vTop) * h;
         }
         double sp = A.vTop;
         if (h0.nc > 0)
         {
            const int p = shape[it.s]->part;
            sp += nb == 1 ? skPOwn[(size_t)p] : skPAll[(size_t)p];
         }
         return g > sp * h;
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
         anyMassless = false;
         for (int p = 0; p < np; p++)
         {
            const int n = part[p]->NumCoords();
            for (int c = 0; c < n; c++)
            {
               const int i = coordOff[p] + c, d = part[p]->DofOf(c);
               ci[i] = part[p]->Coord(c);
               cUnk[i] = d >= 0 ? partOff[p] + d : -1;
               anyMassless = anyMassless || (cUnk[i] >= 0 && ci[i].invM == 0.0);
            }
         }
      }
      bool Massless(const int c) const { return anyMassless && cUnk[c] >= 0 && ci[c].invM == 0.0; }   // a massless coordinate with an unknown (design section 3)
      // the parts' rows in the line search's speed units: over the mass, times the reach. A massless row (design section 3 c): its residual turned into a move by its own slope
      // (the last Jacobian's diagonal, mlJii), then into a speed (2 / h)
      double CoordRes2(const double *Res, const double h) const
      {
         double r = 0.0;
         for (int c = 0; c < nCoord; c++)
            if (cUnk[c] >= 0)
            {
               const double s = Massless(c) ? (mlJii[c] != 0.0 ? Res[cUnk[c]] / mlJii[c] * (2.0 / h) : 0.0) : Res[cUnk[c]] * ci[c].invM * ci[c].reach;
               r += s * s;
            }
         return r;
      }
      // After the solve: did a skipped contact get reached after all? Then the attempt runs again with none skipped.
      // the path test: a static shape (straight, or curved by its cover) and a ball pair are tested over the whole move (PathReach, the path kit); a moving shape at the move's middle and end
      bool FarStillFar(const double *xv)
      {
         if (!ContactProf::CheckOn())
            return FarStillFarRun(xv);
         const auto t0 = std::chrono::steady_clock::now();   // (VPX_CHECK_PROFILE: the call's time, diagnostics only)
         const bool r = FarStillFarRun(xv);
         bool any = false;
         for (int k = 0; k < ni && !any; k++)
            any = farNow[k] != 0;
         if (any)
         {
            const uint64_t ns = (uint64_t)std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count();
            RubberContact::PathStats &ps = RubberContact::pathStats;
            ps.calls++;
            ps.ns += ns;
            ps.nsWorst = std::max(ps.nsWorst, ns);
         }
         return r;
      }
      // VPX_FINDER_CHECK=1 (diagnostics only: the counters are put back, nothing else changes): a cover against its shape's own probe at the path's end e. The probe's gap there g = -delta
      // (R + g: the distance from the surface the probe measures). A zero-length path at e must reach with a ball of R + g + 1 um (the cover lies within its sagitta of the curve, and its
      // test adds the sagitta) and must not with R + g - 2 sagitta - 1 um (the cover can lie a sagitta nearer than the curve, and the test adds one more)
      void PathSelfCheck(const Shape &sh, const V3 &e, const double R, const int ia) const
      {
         int kind = -1;
         const double sag = sh.PathCoverSag(kind);
         if (sag < 0.0 || kind < 0 || kind > 2)
            return;
         const RubberContact::PathStats save = RubberContact::pathStats;
         Hit h;
         sh.Probe(e, R, h);
         if (h.delta <= -1.0)   // the probe's "not there" (no face in its reach, or the surface switched off): no gap to compare, as PathSelfCheckSwept
         {
            RubberContact::pathStats = save;
            return;
         }
         const double g = -h.delta, rIn = R + g + 1e-6, rOut = R + g - 2.0 * sag - 1e-6;
         const int a = sh.PathReach(e, V3(), rIn);
         const int b = rOut > 0.0 ? sh.PathReach(e, V3(), rOut) : 0;
         RubberContact::pathStats = save;
         RubberContact::PathCheck &pc = RubberContact::pathCheck;
         pc.n[kind]++;
         if (a == 1 && b == 0)
            return;
         pc.bad[kind]++;
         if (!pc.warned[kind])
         {
            pc.warned[kind] = true;
            char buf[400];
            snprintf(buf, sizeof(buf), "shape %s, ball %d at (%.9f, %.9f, %.9f) m, probe gap %.9g m, sagitta %.3g m: radius R + g + 1 um -> %d (want 1), R + g - 2 sag - 1 um -> %d (want 0)", sh.Name().c_str(), ia,
                     e.x, e.y, e.z, g, sag, a, b);
            pc.first[kind] = buf;
         }
      }
      // the part finder, the path test of a skipped contact on a surface its part carries (group A; paper: runs/2026-11-13-finder-P6/paper_groupA.txt). The ball's centre moves
      // straight over the step, its middle P; the part moves from where it is to its end trial, its middle pose the middle trial. Every point of the part moves relative to P by at most
      // m(P) over the step (Shape::PathMotion) and the ball's centre at most dLen / 2 from P; a distance changes no faster than its two points move (1-Lipschitz), so over the whole step
      // the depth is at most the probe's depth at P with the part at its middle pose (the exact point probe, this attempt's list, less the item's zero) plus L (dLen / 2 + m / 2), L the
      // probe's largest slope (PathSlope: 1 for a true distance).
      // 1 = reached (that bound above 0: the attempt again with it in), 0 = clear, -1 = no motion term for this surface (the two probes)
      mutable Hit lawHit;
      int SweptLaw(const Item &it, const Shape &sh, const V3 &P, const double dLen, const double Rr) const
      {
         Part *const pp = part[sh.part];
         sh.trial = (const double *)qEnd + coordOff[sh.part];
         const double m = sh.PathMotion(P, dLen);
         int r = -1;
         if (m >= 0.0)
         {
            sh.trial = (const double *)qMid + coordOff[sh.part];
            if (nearOn && it.nearK >= 0)
               pp->nearNow = &nearL[(size_t)it.nearK];
            sh.Probe(P, Rr, lawHit);
            pp->nearNow = nullptr;
            ZeroHit(it, lawHit);
            r = lawHit.delta + sh.PathSlope() * (0.5 * dLen + 0.5 * m) > 0.0 ? 1 : 0;
         }
         sh.trial = nullptr;
         return r;
      }
      // VPX_FINDER_CHECK=1 (diagnostics only: the counters are put back, nothing else changes), the law above against the probe: at the path's end with the part at its end trial, and at
      // the path's middle with the part at its middle trial, the probe's gap g there (R + g: its distance from the surface, this attempt's list). Both places are in the step, so the law
      // at a zero-length path there must reach with a ball of R + g + 1 um (g from the zeroed probe, as the law reads it; a one-way bound). A probe that answers "not there" (-1 or less: a part switched off) or a
      // centre inside the surface (R + g + 1 um <= 0) is not checked
      void PathSelfCheckSwept(const Item &it, const Shape &sh, const V3 &c, const V3 &d, const double R, PartFinderCount &pc) const
      {
         const RubberContact::PathStats save = RubberContact::pathStats;
         Part *const pp = part[sh.part];
         for (int w = 0; w < 2; w++)
         {
            const V3 e = w == 0 ? c + d : c + d * 0.5;
            Hit h;
            sh.trial = (const double *)(w == 0 ? qEnd : qMid) + coordOff[sh.part];
            if (nearOn && it.nearK >= 0)
               pp->nearNow = &nearL[(size_t)it.nearK];
            sh.Probe(e, R, h);
            pp->nearNow = nullptr;
            sh.trial = nullptr;
            if (h.delta <= -1.0)
               continue;
            ZeroHit(it, h);   // the gap the step reads: the probe less the item's zero (a union of terms: each term less its own zero)
            const double g = -h.delta, rIn = R + g + 1e-6;
            if (rIn <= 0.0)
               continue;
            const int a = SweptLaw(it, sh, e, 0.0, rIn);
            pc.pathChecks++;
            if (a == 1)
               continue;
            pc.pathBad++;
            if (partFinderFirstPathBad.empty())
            {
               char buf[400];
               snprintf(buf, sizeof(buf), "shape %s, ball %d at the path's %s (%.9f, %.9f, %.9f) m, probe gap %.9g m: radius R + g + 1 um -> %d (want 1)", sh.Name().c_str(), it.a,
                        w == 0 ? "end" : "middle", e.x, e.y, e.z, g, a);
               partFinderFirstPathBad = buf;
            }
         }
         RubberContact::pathStats = save;
      }
      bool FarStillFarRun(const double *xv)
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
         const bool prof = ContactProf::CheckOn();
         for (int k = 0; k < ni; k++)
            if (farNow[k])
            {
               {
                  const Item &it = item[k];
                  const Ball &A = ball[it.a];
                  const V3 dA(xv[6 * it.a], xv[6 * it.a + 1], xv[6 * it.a + 2]);
                  int r;
                  if (it.b >= 0)   // a ball pair: the other centre is fixed in the relative move (both move straight over the step)
                  {
                     const Ball &B = ball[it.b];
                     const V3 dB(xv[6 * it.b], xv[6 * it.b + 1], xv[6 * it.b + 2]);
                     r = RubberContact::PathWithin(RubberContact::PathPointD2(A.c - B.c, dA - dB, V3()), A.R + B.R - it.dZero) ? 1 : 0;
                     if (prof)
                        RubberContact::pathStats.kit++;
                  }
                  else if (shape[it.s]->part >= 0)   // a part's surface: a fixed one (a bumper's body, a disc) its own test; one that moves with it the rigid-body law (group A)
                  {
                     const Shape &sh = *shape[it.s];
                     r = sh.PathReach(A.c, dA, A.R);
                     const bool fixedS = r >= 0;
                     if (!fixedS)
                        r = SweptLaw(it, sh, A.c + dA * 0.5, dA.Length(), A.R);
                     if (partFinderProf || RubberContact::PathCheckOn())
                     {
                        PartFinderCount &pc = partFinderCount[typeid(*part[sh.part]).name()];
                        if (fixedS)
                           pc.pathFixed++;
                        else if (r < 0)
                           pc.pathTwo++;
                        else
                        {
                           pc.pathTested++;
                           (r > 0 ? pc.pathReached : pc.pathClear)++;
                           if (RubberContact::PathCheckOn())
                              PathSelfCheckSwept(it, sh, A.c, dA, A.R, pc);
                        }
                     }
                  }
                  else
                     r = shape[it.s]->PathReach(A.c, dA, A.R);
                  if (r >= 0 && it.b < 0 && RubberContact::PathCheckOn())
                     PathSelfCheck(*shape[it.s], A.c + dA, A.R, it.a);
                  if (r >= 0)
                  {
                     if (prof)
                     {
                        RubberContact::PathStats &ps = RubberContact::pathStats;
                        ps.tested++;
                        if (r > 0)
                           ps.reached++;
                        else if (it.b >= 0)
                           ps.clearC++;   // (a surface host counts its clears where it decided: a, b or c)
                     }
                     if (r > 0)
                        return false;
                     continue;
                  }
               }
               Hit h1, hm;
               ProbeMove(item[k], b1, bm, h1, hm);
               if (prof)
                  RubberContact::pathStats.twoPoint++;
               if (h1.delta > 0.0 || hm.delta > 0.0)
               {
                  if (prof)
                     RubberContact::pathStats.reachedTwo++;
                  return false;
               }
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

      // A shape's raw hit less the item's zero. A union of terms with term zeros: each term's depth less its own zero (a_t - z_t), the union rebuilt from the zeroed terms as
      // the probe builds it from raw ones (hingegroup ProbeTerms: (sum_t s_t d_t^2.5)^0.4 over the pressing terms, else the deepest + term), so force (the terms' chords),
      // book (V of the union) and gating read the one potential sum_t s_t V(a_t - z_t)
      // patch_zero: a static item's own probe with the given per-piece zeros (the item's or a trial list), every patch shifted by its own zero; nothing to take off after
      void ProbeWithZeros(const Item &it, const std::vector<std::pair<int, double>> &zl, Hit &h) const
      {
         h.ne = 0;
         h.kb = 0.0;
         h.ktb = 0.0;
         h.nc = 0;
         h.npat = 0;
         RubberContact::PatchZeros() = &zl;
         RubberContact::PatchZeroNet() = false;
         shape[it.s]->Probe(ball[it.a].c, ball[it.a].R, h);
         RubberContact::PatchZeros() = nullptr;
         PzNetCheck(it);
      }
      // patch_zero (F2, F6): a probe that found a raw patch whose pieces carry different zeros (the largest-entry rule is only a net: F4's start merge keeps a patch's zeros equal
      // within a step). Only an Eval's probes are kept (pzNetS: the first such shape of this Eval; Eval clears it), and only the accepted solve's final Eval is told (PzNetTell,
      // at Accept): counted, WARN and a notice once per table start; the next Accept gives that patch the parts' energy (rule 5). A solve's iterates are not counted
      mutable bool pzEvalOn = false;
      mutable int pzNetS = -1;
      void PzNetCheck(const Item &it) const
      {
         if (!RubberContact::PatchZeroNet())
            return;
         RubberContact::PatchZeroNet() = false;
         if (pzEvalOn && pzNetS < 0)
            pzNetS = it.s;
      }
      void PzNetTell()
      {
         if (pzNetS < 0)
            return;
         PzCount(PzStats()[4], "nets (an accepted solve's final Eval: a raw patch with different zeros, F2)");
         const Shape *const sh = shape[pzNetS];
         ContactSystem::NoticeRetired("patch zero net: " + (sh != nullptr ? sh->Name() : std::string("?")));
         pzNetS = -1;
      }
      // patch_zero (Z4): first sight of a static mesh or tube: each raw patch's zero = its raw level-form depth (its energy 0), an entry for each of its pieces, dZero the largest.
      // A probe that names no pieces (patch_friction 0, or a path without ids): the scalar zero as before
      void FirstSightPatches(Item &it, Hit &hp) const
      {
         static const std::vector<std::pair<int, double>> none;
         ProbeWithZeros(it, none, hp);
         it.tZero.clear();
         it.dZero = 0.0;
         for (int p = 0; p < hp.npat; p++)
            if (hp.pat[p].d > 0.0)
            {
               for (int q = 0; q < hp.pat[p].pn; q++)
                  if (hp.pc[hp.pat[p].p0 + q] != RubberContact::PatchPieceList::kNoPiece)
                     it.tZero.emplace_back(hp.pc[hp.pat[p].p0 + q], hp.pat[p].d);
               it.dZero = std::max(it.dZero, hp.pat[p].d);
               PzCount(PzStats()[1], "first-sight patches");
            }
         if (it.tZero.empty())
            it.dZero = std::max(hp.delta, 0.0);
      }
      // patch_zero (Z5, rule 5): at Accept, per patch of the end probe (raw patches, each shifted by its zero z, the largest entry of its pieces): d' its shifted depth (its level
      // form's; no node above z: its deepest raw node less z), z = max(0, min(z, z + d')); every raw-pressed piece of the patch takes it (pieces that joined in the step included),
      // entries of pieces not raw-pressed dropped. A patch whose entries differ (a join at the step's end): one zero by bisection, the joined patch's level-form sum (s25) equal
      // to its parts' with the zeros they carried (zero_rule.md rule 5): only a net with F4 (a raw join inside a step with mixed zeros, F2).
      void BackOutPatches(Item &it) const
      {
         Hit hr;
         ProbeWithZeros(it, it.tZero, hr);
         thread_local std::vector<std::pair<int, double>> nz, trial;
         nz.clear();
         const auto s25 = [](const Hit &h, const int k) { return k < h.npat && h.pat[k].d > 0.0 ? std::pow(h.pat[k].d, 2.5) : 0.0; };
         for (int p = 0; p < hr.npat; p++)
         {
            const HitPatch &P = hr.pat[p];
            double zlo = INFINITY, zhi = 0.0, deep = -INFINITY;
            for (int q = 0; q < P.pn; q++)
            {
               const int id = hr.pc[P.p0 + q];
               for (const std::pair<int, double> &e : it.tZero)
                  if (e.first == id)
                  {
                     zlo = std::min(zlo, e.second);
                     zhi = std::max(zhi, e.second);
                     break;
                  }
            }
            if (zhi <= 0.0)
               continue;   // a patch with no zero: none kept
            for (int q = 0; q < P.pn; q++)
               deep = std::max(deep, hr.pd[P.p0 + q].d + zhi);   // (each piece's raw depth: the probe's is from the patch's zero)
            double z = std::max(0.0, std::min(zhi, zhi + P.d));
            if (zlo < zhi)   // rule 5
            {
               Hit hj;
               RubberContact::PatchZerosPerPiece() = true;
               ProbeWithZeros(it, it.tZero, hj);
               RubberContact::PatchZerosPerPiece() = false;
               const double target = s25(hj, p);
               double lo = 0.0, hi = std::max(deep, 0.0);
               for (int r = 0; r < 60; r++)
               {
                  const double mid = 0.5 * (lo + hi);
                  trial.assign(it.tZero.begin(), it.tZero.end());
                  for (int q = 0; q < P.pn; q++)
                  {
                     const int id = hr.pc[P.p0 + q];
                     size_t j = 0;
                     while (j < trial.size() && trial[j].first != id)
                        j++;
                     if (j == trial.size())
                        trial.emplace_back(id, mid);
                     else
                        trial[j].second = mid;
                  }
                  Hit ht;
                  ProbeWithZeros(it, trial, ht);
                  (s25(ht, p) > target ? lo : hi) = mid;
               }
               z = hi;
               PzCount(PzStats()[2], "rule-5 joins at Accept");
            }
            if (z <= 0.0)
               continue;
            for (int q = 0; q < P.pn; q++)
               if (hr.pd[P.p0 + q].d + zhi > 0.0 && hr.pc[P.p0 + q] != RubberContact::PatchPieceList::kNoPiece)   // raw-pressed
                  nz.emplace_back(hr.pc[P.p0 + q], z);
         }
         it.tZero.assign(nz.begin(), nz.end());
         it.dZero = 0.0;
         for (const std::pair<int, double> &e : it.tZero)
            it.dZero = std::max(it.dZero, e.second);
      }
      // patch_zero (F4, start merge): the solve's last Eval (at its answer) found its union joining pressed start patches whose zeros differ (LinkUnion, pzMg). Before the same h
      // is solved again, each such group's patches take ONE zero at the step's start position: bisection (60 halvings) on z with their summed level-form sums (s25) on the
      // start probe equal to the sum with their own zeros (zero_rule.md rule 5: the start energy kept); every piece of them that has an entry or is raw-pressed gets it.
      // Zeros change only here, between solves; the join then happens with equal zeros. A patch pair (their first pieces) merges at most once per step (pzMerged): a group
      // whose differing pairs all merged already is left to the net (F2). True: a zero changed, solve again
      // patch_zero (F7): the one call of the start merge, on the final solve of every step before it is accepted (Step: the normal and a halved step's solve, before the link
      // order and stop checks, so no other solve comes between it and their Accept; the stop search's final solve, before its Accept). True: a zero changed, the same step
      // length is solved again. (The energy step has no merge: the net, F2, checks it)
      bool MergeBeforeAccept() { return set.patchZero && StartMerges(); }   // (patchZero implies linkWindows: ToStep, F9)
      bool StartMerges()
      {
         bool any = false;
         thread_local std::vector<int> jp;
         thread_local std::vector<std::pair<int, double>> trial;
         const auto s25 = [](const Hit &h, const int k) { return k < h.npat && h.pat[k].d > 0.0 ? std::pow(h.pat[k].d, 2.5) : 0.0; };
         for (int k = 0; k < (int)pzMg.size() && k < ni; k++)
         {
            if (pzMg[k].empty() || !PatchZeroItem(item[k]))
               continue;
            Item &it = item[k];
            const std::vector<std::pair<int, int>> &mg = pzMg[k];
            Hit hr;
            ProbeWithZeros(it, it.tZero, hr);   // (the balls are where the step starts)
            const auto patchZ = [&](const HitPatch &P) {   // a patch's zero in hr: the largest entry of its pieces
               double z = 0.0;
               for (int q = 0; q < P.pn; q++)
                  z = std::max(z, RubberContact::PieceZero(it.tZero, hr.pc[P.p0 + q]));
               return z;
            };
            bool changed = false;
            for (size_t g0 = 0; g0 < mg.size();)
            {
               size_t g1 = g0;
               while (g1 < mg.size() && mg[g1].first == mg[g0].first)
                  g1++;
               jp.clear();
               for (int p = 0; p < hr.npat; p++)   // the group's patches in hr: those holding one of its pieces
               {
                  const HitPatch &P = hr.pat[p];
                  bool in = false;
                  for (int q = 0; q < P.pn && !in; q++)
                     for (size_t m = g0; m < g1 && !in; m++)
                        in = hr.pc[P.p0 + q] == mg[m].second;
                  if (in && P.pn > 0)
                     jp.push_back(p);
               }
               g0 = g1;
               bool fresh = false;
               for (size_t a = 0; a < jp.size(); a++)
                  for (size_t b = a + 1; b < jp.size(); b++)
                  {
                     const HitPatch &A = hr.pat[jp[a]], &B = hr.pat[jp[b]];
                     const std::pair<int, int> pr(std::min(hr.pc[A.p0], hr.pc[B.p0]), std::max(hr.pc[A.p0], hr.pc[B.p0]));
                     fresh = fresh || (patchZ(A) != patchZ(B) && std::find(pzMerged.begin(), pzMerged.end(), pr) == pzMerged.end());
                  }
               if (!fresh)
                  continue;
               double target = 0.0, deep = 0.0;
               for (const int p : jp)
               {
                  const HitPatch &P = hr.pat[p];
                  const double zp = patchZ(P);
                  target += s25(hr, p);
                  for (int q = 0; q < P.pn; q++)
                     deep = std::max(deep, hr.pd[P.p0 + q].d + zp);   // (each piece's raw depth)
               }
               const auto withZero = [&](const double z) {   // trial = the item's zeros with z on every piece of the group that has an entry or is raw-pressed
                  trial.assign(it.tZero.begin(), it.tZero.end());
                  for (const int p : jp)
                  {
                     const HitPatch &P = hr.pat[p];
                     const double zp = patchZ(P);
                     for (int q = 0; q < P.pn; q++)
                     {
                        const int id = hr.pc[P.p0 + q];
                        if (id == RubberContact::PatchPieceList::kNoPiece)
                           continue;
                        size_t j = 0;
                        while (j < trial.size() && trial[j].first != id)
                           j++;
                        if (j < trial.size())
                           trial[j].second = z;
                        else if (hr.pd[P.p0 + q].d + zp > 0.0)
                           trial.emplace_back(id, z);
                     }
                  }
               };
               double lo = 0.0, hi = deep;
               for (int r = 0; r < 60; r++)
               {
                  const double mid = 0.5 * (lo + hi);
                  withZero(mid);
                  Hit ht;
                  ProbeWithZeros(it, trial, ht);
                  double s = 0.0;
                  for (const int p : jp)
                     s += s25(ht, p);
                  (s > target ? lo : hi) = mid;
               }
               withZero(hi);
               trial.erase(std::remove_if(trial.begin(), trial.end(), [](const std::pair<int, double> &e) { return e.second <= 0.0; }), trial.end());
               it.tZero.assign(trial.begin(), trial.end());
               for (size_t a = 0; a < jp.size(); a++)
                  for (size_t b = a + 1; b < jp.size(); b++)
                  {
                     const int pa = hr.pc[hr.pat[jp[a]].p0], pb = hr.pc[hr.pat[jp[b]].p0];
                     pzMerged.emplace_back(std::min(pa, pb), std::max(pa, pb));
                  }
               PzCount(PzStats()[3], "start merges (F4: start patches with different zeros joined by a step's union, one zero at the start)");
               changed = true;
               if (it.tZero.empty())
                  break;
               ProbeWithZeros(it, it.tZero, hr);   // (the next group's zeros from the item's new ones)
            }
            if (!changed)
               continue;
            it.dZero = 0.0;
            for (const std::pair<int, double> &e : it.tZero)
               it.dZero = std::max(it.dZero, e.second);
            if (nowStamp == stamp)   // the step's start probe kept for the next attempt (BeginAttempt) reads the new zeros
               Probe(it, ball, hNow[k], true);
            any = true;
         }
         return any;
      }
      void ZeroHit(const Item &it, Hit &hp) const
      {
         if (PatchZeroItem(it))   // patch_zero: nothing to take off (the step pieces' entries too: each from its own patch's zero, StepPieces, D3)
            return;
         for (int j = 0; j < hp.ne; j++)   // step pieces: each piece entry from the object's zero, as the merged depth (a ridge's is already)
            if (!hp.en[j].zeroed)
               hp.en[j].d -= it.dZero;
         if (hp.nt == 0 || it.tZero.empty())
         {
            hp.delta -= it.dZero;
            for (int j = 0; j < hp.npat; j++)   // patch_friction: each patch from the object's zero, as the merged depth
               hp.pat[j].d -= it.dZero;
            for (int j = 0; j < hp.npc; j++)   // (level_stiffness: each piece likewise)
               hp.pd[j].d -= it.dZero;
            return;
         }
         double s25 = 0.0, best = -INFINITY;
         V3 g, bestN = hp.n;
         for (int t = 0; t < hp.nt; t++)
         {
            HitTerm &ht = hp.t[t];
            ht.d -= it.TermZero(ht.id);
            if (ht.d > 0.0)
            {
               const double d15 = ht.d * std::sqrt(ht.d);
               s25 += ht.s * ht.d * d15;
               g = g + ht.n * (ht.s * d15);
            }
            if (ht.s > 0.0 && ht.d > best)
            {
               best = ht.d;
               bestN = ht.n;
            }
         }
         hp.delta = best;
         hp.n = bestN;
         hp.grad = 1.0;
         if (s25 > 0.0)
         {
            const double d = std::pow(s25, 0.4);
            const V3 gr = g / (d * std::sqrt(d));
            const double gl = gr.Length();
            hp.delta = d;
            if (gl > 1e-6)
            {
               hp.n = gr / gl;
               hp.grad = gl;
            }
         }
      }
      // start: the step-start probe (ChooseStep): a shape's ProbeStart, which may write entries; every other probe writes none
      void Probe(const Item &it, const Ball *bs, Hit &hp, const bool start = false) const
      {
         const Ball &A = bs[it.a];
         hp.ne = 0;   // (step pieces: hits are kept and reused; only ProbeStart writes entries)
         hp.kb = 0.0;
         hp.ktb = 0.0;
         hp.nc = 0;   // (hits are kept and reused: a fixed surface writes no columns, and must not inherit a carried one's)
         hp.npat = 0;   // (patch_friction: a static surface writes its patches)
         if (it.b < 0)
         {
            const Shape *const sh = shape[it.s];
            const bool listed = nearOn && it.nearK >= 0 && sh->part >= 0;
            if (listed)
               part[sh->part]->nearNow = &nearL[(size_t)it.nearK];
            const bool pz = PatchZeroItem(it);   // patch_zero (D6): the item's live per-patch zeros to the probe, around its call
            if (pz)
            {
               RubberContact::PatchZeros() = &it.tZero;
               RubberContact::PatchZeroNet() = false;
            }
            if (start)
               sh->ProbeStart(A.c, A.R, hp);
            else
               sh->Probe(A.c, A.R, hp);
            if (pz)
            {
               RubberContact::PatchZeros() = nullptr;
               PzNetCheck(it);
            }
            if (listed)
               part[sh->part]->nearNow = nullptr;
            ZeroHit(it, hp);
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
            it.mRed = B.held ? A.m : A.m * B.m / (A.m + B.m);   // (a held ball: infinite mass, KickerDesign 1h)
            invMt += B.held ? 0.0 : 1.0 / B.m + B.R * B.R / B.I();
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
      // patch_friction (piece 2b-1): a group's load from its own depths, the same laws as the merged contact's (Item::Elastic's chord, the damping above Eval's friction):
      // fe the chord of F over d0 .. d1, kc its slope, fd the damping on d1 - d0, kn Newton's normal stiffness (as Eval's kn)
      void GroupLoad(const Item &L, const double h, FG &G) const
      {
         const double d0 = G.d0, d1 = G.d1, dm = G.dm;
         const bool di = set.dampingIntegral;
         const double fe = L.Chord(d0, d1);
         const double kc = std::fabs(d1 - d0) > 1e-13 ? (L.F(d1) - fe) / (d1 - d0) : 0.5 * L.K(std::max(std::max(dm, d1), 0.0));
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
         const double dd = d1 - d0, fd = cD * dd / h;
         G.fe = fe;
         G.fd = fd;
         G.N = std::max(0.0, fe + fd);
         G.kc = kc;
         G.aE = G.dadfE = 0.0;
         if (di || set.rollingPatch)
            L.PatchOfForce(fe, G.aE, G.dadfE);
         const double cDp = L.gammaKK > 0.0 ? (dm > 0.0 ? 0.5 * L.gammaKK * 0.5 * std::sqrt(L.reff / dm) : 0.0) : 1.5 * L.alpha * kc;
         const double dSl = cDs + (hcChord ? 1.5 * L.alpha * (fe + kc * dd) : 0.0);
         G.kn = fe + fd <= 0.0 ? 0.0 : di ? h * kc + dSl : L.gammaKK > 0.0 ? h * kc + kkJac : h * kc + cD + cDp * dd;
      }
      static bool SharePiece(const Hit &a, const HitPatch &pa, const Hit &b, const HitPatch &pb)
      {
         for (int i = 0; i < pa.pn; i++)
            for (int j = 0; j < pb.pn; j++)
               if (a.pc[pa.p0 + i] == b.pc[pb.p0 + j])
                  return true;
         return false;
      }
      // patch_friction (plan 2.5 C1/C2): the step's friction groups: the start (h0) and end (h1) patches pressed there, linked by shared pieces. Per group: d0, d1 its start's and
      // end's depth (several patches: the level form's one depth (sum d^2.5)^0.4); at the middle the hm patches that share a piece with it, else its reference patch's plane (its
      // deepest end patch, else start patch; a face's depth is linear in the centre: d + grad n . (c - c_ref)); the end gradient from the reference. Returns the count; fgEnd[j],
      // gStart[i] the group of end patch j, start patch i (-1 not pressed). The memory is read by GroupMemory.
      // link_windows (rule T, to-do 1.33; runs/2026-11-23-link-rule brief_code.md, rule_t.py, check2_out.txt). The gate: windows are needed only when the start has 2 or more
      // pressed patches, or a pressed start or end patch shares no piece with a pressed patch at the other end
      static bool LinkGate(const Hit &h0, const Hit &h1)
      {
         int p0 = 0;
         for (int i = 0; i < h0.npat; i++)
            p0 += h0.pat[i].d > 0.0 ? 1 : 0;
         if (p0 >= 2)
            return true;
         for (int i = 0; i < h0.npat; i++)
            if (h0.pat[i].d > 0.0)
            {
               bool any = false;
               for (int j = 0; j < h1.npat && !any; j++)
                  any = h1.pat[j].d > 0.0 && SharePiece(h0, h0.pat[i], h1, h1.pat[j]);
               if (!any)
                  return true;
            }
         for (int j = 0; j < h1.npat; j++)
            if (h1.pat[j].d > 0.0)
            {
               bool any = false;
               for (int i = 0; i < h0.npat && !any; i++)
                  any = h0.pat[i].d > 0.0 && SharePiece(h0, h0.pat[i], h1, h1.pat[j]);
               if (!any)
                  return true;
            }
         return false;
      }
      // the step's pieces (the patches' of h0, h1, hm and the windows' ends; lwIds sorted, unique) and their components (lwPar: union-find); the windows (lwS) with their ends' places
      // (lwA, lwB). Per item, the last Eval's order check: lwCut[k] the fraction of the step to cut at (> 0), -1 none, -2 rule T's feed differs from the union's with no cut point;
      // lwAdd[k] the start-end patch links the windows added beyond shared pieces. Scratch of Eval, capacity reused
      std::vector<int> lwIds, lwPar, lwA, lwB, lwOrd, lwP2;
      std::vector<RubberContact::LinkWindow> lwS;
      std::vector<uint64_t> lwSet, lwU, lwFT, lwFU;
      std::vector<double> lwCut;
      std::vector<int> lwAdd;
      int lwCuts = 0;   // the order cuts of this step (Step)
      // patch_zero (F4): per item, the last Eval's start patches to merge (group number, piece id; LinkUnion); the patch pairs merged in this step (their first pieces; StartMerges)
      std::vector<std::vector<std::pair<int, int>>> pzMg;
      std::vector<std::pair<int, int>> pzMerged;
      static constexpr int kLinkCuts = 4;
      int LwIx(const int id) const { return (int)(std::lower_bound(lwIds.begin(), lwIds.end(), id) - lwIds.begin()); }
      int LwFind(int x) { while (lwPar[x] != x) x = lwPar[x] = lwPar[lwPar[x]]; return x; }
      int LwRoot(const Hit &h, const HitPatch &P) { return P.pn > 0 ? LwFind(LwIx(h.pc[P.p0])) : -1; }   // a patch's component (-1: no pieces)
      bool LwIn(const Hit &h, const HitPatch &P, const int root)   // a piece of P in the component root
      {
         for (int q = 0; root >= 0 && q < P.pn; q++)
            if (LwFind(LwIx(h.pc[P.p0 + q])) == root)
               return true;
         return false;
      }
      // the counters (always on; running totals logged at 1, 2, 4, ...: the bench never closes the table)
      static void LwCount(std::atomic<long long> &n, const long long add, const char *what, const unsigned ms)
      {
         if (add <= 0)
            return;
         const long long was = n.fetch_add(add), c = was + add;
         if (was == 0 || (c ^ was) > was)
            PLOGI << "link_windows: " << c << " " << what << " (table time " << ms << " ms)";
      }
      struct LinkStats
      {
         std::atomic<long long> gated { 0 }, windows { 0 }, added { 0 }, cuts { 0 }, nets { 0 };
      };
      static LinkStats &LinkCounters()
      {
         static LinkStats s;
         return s;
      }
      // the union rule (brief 3): the piece components through each pressed start and end patch's own pieces and every windowed link; a start and an end patch in one component are
      // one group. Then the order check (brief 4): rule T's forward sweep (rule_t.py rule_T: each start patch's pieces seeded with it; at each window's t_in in time order the
      // component of the links active then takes the union of its pieces' sets); an end patch fed by fewer start patches than the union's: the cut point in its component
      template <class Join> void LinkUnion(const Hit &h0, const Hit &h1, const Hit &hm, const V3 &c0, const V3 &c1, const int k, const Join &join)
      {
         const int n0 = h0.npat, n1 = h1.npat;
         lwS.clear();
         const Shape *const sh = shape[item[k].s];
         if (sh != nullptr)
         {
            static const std::vector<std::pair<int, double>> none;
            const bool pz = PatchZeroItem(item[k]);   // patch_zero (D4): the links' reach with no zero (the patches are found on raw depth)
            if (pz)
               RubberContact::PatchZeros() = &none;
            sh->LinkWindows(c0, c1, lwS);
            if (pz)
               RubberContact::PatchZeros() = nullptr;
         }
         LinkStats &ls = LinkCounters();
         LwCount(ls.gated, 1, "gated Evals", pfMs);
         LwCount(ls.windows, (long long)lwS.size(), "link windows computed", pfMs);
         lwIds.clear();
         for (const Hit *hh : { &h0, &h1, &hm })
            for (int q = 0; q < hh->npc; q++)
               lwIds.push_back(hh->pc[q]);
         for (const RubberContact::LinkWindow &w : lwS)
         {
            lwIds.push_back(w.a);
            lwIds.push_back(w.b);
         }
         std::sort(lwIds.begin(), lwIds.end());
         lwIds.erase(std::unique(lwIds.begin(), lwIds.end()), lwIds.end());
         const int nPc = (int)lwIds.size();
         lwPar.resize(nPc);
         for (int x = 0; x < nPc; x++)
            lwPar[x] = x;
         const auto unite = [&](const int a, const int b) { lwPar[LwFind(a)] = LwFind(b); };
         for (const Hit *hh : { &h0, &h1 })   // a pressed patch's pieces are joined (one patch)
            for (int p = 0; p < hh->npat; p++)
            {
               const HitPatch &P = hh->pat[p];
               if (P.d > 0.0)
                  for (int q = 1; q < P.pn; q++)
                     unite(LwIx(hh->pc[P.p0]), LwIx(hh->pc[P.p0 + q]));
            }
         const int nw = (int)lwS.size();
         lwA.resize(nw);
         lwB.resize(nw);
         for (int w = 0; w < nw; w++)
         {
            lwA[w] = LwIx(lwS[w].a);
            lwB[w] = LwIx(lwS[w].b);
            unite(lwA[w], lwB[w]);
         }
         if ((int)lwAdd.size() <= k)
            lwAdd.resize(k + 1, 0);
         if ((int)lwCut.size() <= k)
            lwCut.resize(k + 1, -1.0);
         lwAdd[k] = 0;
         for (int i = 0; i < n0; i++)
            if (h0.pat[i].d > 0.0)
            {
               const int ri = LwRoot(h0, h0.pat[i]);
               for (int j = 0; j < n1; j++)
                  if (h1.pat[j].d > 0.0 && ri >= 0 && LwRoot(h1, h1.pat[j]) == ri)
                  {
                     join(i, n0 + j);
                     if (!SharePiece(h0, h0.pat[i], h1, h1.pat[j]))
                        lwAdd[k]++;
                  }
            }
         if (PatchZeroItem(item[k]))   // patch_zero (F4, F5): the RAW-pressed start patches (a piece's raw depth > 0, whatever the patch's shifted depth) that share a union
         {                             // component through any of their pieces, zeros differing: one group (group number, each piece id) in pzMg[k], for StartMerges (one zero at the step's start)
            const Item &it = item[k];
            const auto zeroOf = [&](const HitPatch &P) {
               double z = 0.0;
               for (int q = 0; q < P.pn; q++)
                  z = std::max(z, it.TermZero(h0.pc[P.p0 + q]));
               return z;
            };
            const auto rawPressed = [&](const HitPatch &P) {   // (a piece's raw depth: its probe depth, from the patch's zero, plus that zero)
               const double z = zeroOf(P);
               for (int q = 0; q < P.pn; q++)
                  if (h0.pd[P.p0 + q].d + z > 0.0)
                     return true;
               return false;
            };
            if ((int)pzMg.size() <= k)
               pzMg.resize(k + 1);
            pzMg[k].clear();
            thread_local std::vector<char> taken;
            thread_local std::vector<int> grp, pp, owner;
            taken.assign((size_t)n0, 0);
            pp.resize((size_t)n0);
            owner.assign((size_t)nPc, -1);
            const auto findP = [&](int a) { while (pp[a] != a) a = pp[a] = pp[pp[a]]; return a; };
            for (int i = 0; i < n0; i++)   // the patches joined through any piece's component (a raw-pressed patch at its zero need not have its pieces united above)
            {
               pp[i] = i;
               if (!rawPressed(h0.pat[i]))
               {
                  taken[i] = 1;
                  continue;
               }
               for (int q = 0; q < h0.pat[i].pn; q++)
               {
                  if (h0.pc[h0.pat[i].p0 + q] == RubberContact::PatchPieceList::kNoPiece)
                     continue;
                  const int r = LwFind(LwIx(h0.pc[h0.pat[i].p0 + q]));
                  if (owner[r] < 0)
                     owner[r] = i;
                  else
                     pp[findP(i)] = findP(owner[r]);
               }
            }
            int ng = 0;
            for (int i = 0; i < n0; i++)
            {
               if (taken[i] || h0.pat[i].pn <= 0)
                  continue;
               const int ri = findP(i);
               const double zi = zeroOf(h0.pat[i]);
               bool differ = false;
               grp.assign(1, i);
               for (int i2 = i + 1; i2 < n0; i2++)
                  if (!taken[i2] && findP(i2) == ri)
                  {
                     grp.push_back(i2);
                     differ = differ || zeroOf(h0.pat[i2]) != zi;
                  }
               for (const int g : grp)
                  taken[g] = 1;
               if (!differ)
                  continue;
               for (const int g : grp)
                  for (int q = 0; q < h0.pat[g].pn; q++)
                     if (h0.pc[h0.pat[g].p0 + q] != RubberContact::PatchPieceList::kNoPiece)
                        pzMg[k].emplace_back(ng, h0.pc[h0.pat[g].p0 + q]);
               ng++;
            }
         }
         // the order check: rule T's feed per pressed end patch against the union's (bit i: start patch i)
         if (n0 == 0 || n1 == 0)
            return;
         const int nWd = (n0 + 63) / 64;
         lwSet.assign((size_t)nPc * nWd, 0u);
         for (int i = 0; i < n0; i++)
         {
            const HitPatch &P = h0.pat[i];
            if (P.d > 0.0)
               for (int q = 0; q < P.pn; q++)
                  lwSet[(size_t)LwIx(h0.pc[P.p0 + q]) * nWd + i / 64] |= 1ull << (i % 64);
         }
         lwOrd.resize(nw);
         for (int w = 0; w < nw; w++)
            lwOrd[w] = w;
         std::sort(lwOrd.begin(), lwOrd.end(), [&](const int a, const int b) { return lwS[a].t0 < lwS[b].t0 || (lwS[a].t0 == lwS[b].t0 && a < b); });
         lwP2.resize(nPc);
         lwU.resize(nWd);
         const auto find2 = [&](int x) { while (lwP2[x] != x) x = lwP2[x] = lwP2[lwP2[x]]; return x; };
         for (const int w : lwOrd)
         {
            const double t = lwS[w].t0;
            for (int x = 0; x < nPc; x++)
               lwP2[x] = x;
            for (int v = 0; v < nw; v++)   // the links active at t
               if (lwS[v].t0 <= t && t <= lwS[v].t1)
                  lwP2[find2(lwA[v])] = find2(lwB[v]);
            const int r = find2(lwA[w]);
            std::fill(lwU.begin(), lwU.end(), 0u);
            for (int x = 0; x < nPc; x++)
               if (find2(x) == r)
                  for (int b = 0; b < nWd; b++)
                     lwU[b] |= lwSet[(size_t)x * nWd + b];
            for (int x = 0; x < nPc; x++)
               if (find2(x) == r)
                  for (int b = 0; b < nWd; b++)
                     lwSet[(size_t)x * nWd + b] = lwU[b];
         }
         lwFT.resize(nWd);
         lwFU.resize(nWd);
         for (int j = 0; j < n1; j++)
         {
            const HitPatch &E = h1.pat[j];
            if (E.d <= 0.0)
               continue;
            std::fill(lwFT.begin(), lwFT.end(), 0u);
            std::fill(lwFU.begin(), lwFU.end(), 0u);
            for (int q = 0; q < E.pn; q++)
               for (int b = 0; b < nWd; b++)
                  lwFT[b] |= lwSet[(size_t)LwIx(h1.pc[E.p0 + q]) * nWd + b];
            const int rj = LwRoot(h1, E);
            for (int i = 0; i < n0; i++)
               if (h0.pat[i].d > 0.0 && rj >= 0 && LwRoot(h0, h0.pat[i]) == rj)
                  lwFU[i / 64] |= 1ull << (i % 64);
            if (lwFT == lwFU)
               continue;
            // the cut (rule_t.py): in the end patch's union component, the earliest t_out of a link before a later t_in there; the cut point is the middle between that t_out and the
            // next t_in in the component
            double tOut = 2.0;
            for (int w = 0; w < nw; w++)
               if (LwFind(lwA[w]) == rj)
                  for (int v = 0; v < nw; v++)
                     if (LwFind(lwA[v]) == rj && lwS[v].t0 > lwS[w].t1)
                     {
                        tOut = std::min(tOut, lwS[w].t1);
                        break;
                     }
            if (tOut > 1.0)
            {
               if (lwCut[k] < 0.0)
                  lwCut[k] = -2.0;
               continue;
            }
            double tNext = 2.0;
            for (int v = 0; v < nw; v++)
               if (LwFind(lwA[v]) == rj && lwS[v].t0 > tOut)
                  tNext = std::min(tNext, lwS[v].t0);
            const double cut = 0.5 * (tOut + tNext);
            if (cut > 0.0 && (lwCut[k] < 0.0 || cut < lwCut[k]))
               lwCut[k] = cut;
         }
      }
      // the order check after a step solves (brief 4; the solve's last Eval is at its answer): the earliest cut of any gated item; true = h cut to it, solve again. Past kLinkCuts cuts,
      // or a difference with no cut point: the step is taken with the union (its groups only merge: energy can only fall), counted, and told (WARN and a notice once per table start)
      bool LinkOrderCut(double &h)
      {
         double cut = 2.0;
         int kc = -1, kn = -1;
         for (int k = 0; k < (int)lwCut.size() && k < ni; k++)
            if (lwCut[k] > 0.0 && lwCut[k] < cut)
            {
               cut = lwCut[k];
               kc = k;
            }
            else if (lwCut[k] == -2.0 && kn < 0)
               kn = k;
         if (kc >= 0 && lwCuts < kLinkCuts && cut < 1.0)
         {
            h *= cut;
            lwCuts++;
            LwCount(LinkCounters().cuts, 1, "order cuts", pfMs);
            return true;
         }
         const int kw = kc >= 0 ? kc : kn;
         if (kw < 0)
            return false;
         LwCount(LinkCounters().nets, 1, "order nets (steps taken with the union)", pfMs);
         const Shape *const sh = shape[item[kw].s];
         ContactSystem::NoticeRetired("link order net: " + (sh != nullptr ? sh->Name() : std::string("?"))
                                      + ": a step's friction groups taken with the union of its links (rule T's order still differed after 4 cuts, or had no cut point; [Contact] link_windows)");
         return false;
      }
      int FrictionGroups(const Hit &h0, const Hit &h1, const Hit &hm, const V3 &c0, const V3 &c1, const V3 &cm, const int k)
      {
         const int n0 = h0.npat, n1 = h1.npat, nn = n0 + n1;
         fgEnd.assign(n1, -1);
         fgJoins = fgSplits = 0;
         gStart.assign(n0, -1);
         if (nn == 0)
            return 0;
         thread_local std::vector<int> par, gid;
         par.resize(nn);
         gid.assign(nn, -1);
         for (int x = 0; x < nn; x++)
            par[x] = x;
         const auto find = [&](int x) { while (par[x] != x) x = par[x] = par[par[x]]; return x; };
         const auto pat = [&](const int x) -> const HitPatch & { return x < n0 ? h0.pat[x] : h1.pat[x - n0]; };
         const bool lw = set.linkWindows && LinkGate(h0, h1);   // link_windows: the union rule (LinkUnion); else today's shared pieces
         if (lw)
            LinkUnion(h0, h1, hm, c0, c1, k, [&](const int a, const int b) { par[find(b)] = find(a); });
         else
            for (int i = 0; i < n0; i++)
               if (h0.pat[i].d > 0.0)
                  for (int j = 0; j < n1; j++)
                     if (h1.pat[j].d > 0.0 && SharePiece(h0, h0.pat[i], h1, h1.pat[j]))
                        par[find(n0 + j)] = find(i);
         int ng = 0;
         for (int x = 0; x < nn; x++)
            if (pat(x).d > 0.0)
            {
               const int r = find(x);
               if (gid[r] < 0)
                  gid[r] = ng++;
               gid[x] = gid[r];
            }
         if (ng == 0)
            return 0;
         if ((int)fgS.size() < ng)
            fgS.resize(ng);
         for (int i = 0; i < n0; i++)
            gStart[i] = gid[i];
         for (int j = 0; j < n1; j++)
            fgEnd[j] = gid[n0 + j];
         for (int g = 0; g < ng; g++)
         {
            FG &G = fgS[g];
            G = FG();
            double s0 = 0.0, s1 = 0.0;
            int c0n = 0, c1n = 0, last0 = -1, last1 = -1, ref = -1;
            for (int x = 0; x < nn; x++)
            {
               if (gid[x] != g)
                  continue;
               const double d = pat(x).d;
               if (x < n0)
               {
                  s0 += d * d * std::sqrt(d);
                  c0n++;
                  last0 = x;
               }
               else
               {
                  s1 += d * d * std::sqrt(d);
                  c1n++;
                  last1 = x;
               }
               if (ref < 0 || (x >= n0 && ref < n0) || ((x >= n0) == (ref >= n0) && d > pat(ref).d))   // the deepest end patch, else start patch
                  ref = x;
            }
            const HitPatch &R = pat(ref);
            const V3 cR = ref >= n0 ? c1 : c0;
            const auto dist = [&R](const V3 &c) { return R.gk == 1 ? (c - RubberContact::SegNearest(R.ga, R.gb, c)).Length() : (c - R.ga).Length(); };
            const auto plane = [&](const V3 &c) {   // plan 2.5 C1: the reference's deepest piece's signed depth at c (its depth falls as c moves along n)
               return R.gk == 0 ? R.d - R.n.Dot(c - cR) : R.d + dist(cR) - dist(c);
            };
            G.d0 = c0n == 1 ? pat(last0).d : c0n > 1 ? std::pow(s0, 0.4) : plane(c0);
            G.d1 = c1n == 1 ? pat(last1).d : c1n > 1 ? std::pow(s1, 0.4) : plane(c1);
            if (c0n > 0 && c1n > 0 && c0n != c1n)   // (counters: more start than end patches a join, more end than start a split)
               (c0n > c1n ? fgJoins : fgSplits)++;
            double sm = 0.0;
            V3 gn;
            int cmn = 0, lastM = -1;
            for (int p = 0; p < hm.npat; p++)
            {
               const HitPatch &P = hm.pat[p];
               if (P.d <= 0.0)
                  continue;
               bool hit = false;
               for (int x = 0; x < nn && !hit; x++)
                  if (gid[x] == g)
                     hit = lw ? LwIn(hm, P, LwRoot(x < n0 ? h0 : h1, pat(x))) : SharePiece(hm, P, x < n0 ? h0 : h1, pat(x));   // (link_windows: by the same components)
               if (!hit)
                  continue;
               const double d15 = P.d * std::sqrt(P.d);
               sm += P.d * d15;
               gn = gn + P.n * d15;
               cmn++;
               lastM = p;
            }
            G.trCmn = cmn;
            if (cmn == 1)
            {
               G.dm = hm.pat[lastM].d;
               G.n = hm.pat[lastM].n;
            }
            else if (cmn > 1)
            {
               G.dm = std::pow(sm, 0.4);
               const double l = gn.Length();
               G.n = l > 0.0 ? gn * (1.0 / l) : R.n;
            }
            else
            {
               G.dm = plane(cm);
               G.n = R.n;
            }
            G.nD = R.n * R.grad;
            G.gm[0] = -G.nD.x;
            G.gm[1] = -G.nD.y;
            G.gm[2] = -G.nD.z;
         }
         return ng;
      }
      // level_stiffness (B-lin, plan runs/2026-11-22-per-peak 2.4b; paper/p2b_avg2.py, p2b_avg2_out.txt; levelk.py): a friction group's shear, rolling and twist stiffness from
      // the level form (the same the elastic energy uses: sum over the peaks of u^p less their passes'), averaged over the step. B: one group's (blS[g]); on false: today's
      // PatchMeans on the group's d0, d1 (one piece, or one peak with no role change: bit for bit today's). m the means of u^p, lev0 / lev1 the level form at the start / end;
      // s05 .. s075 the means' slopes in the end position, vectors (each piece's MeanSlope times its end gradient -n_i)
      struct BLin
      {
         bool on = false;
         Item::Pows m, lev0, lev1;
         V3 s05, s025, s15, s075;
      };
      std::vector<BLin> blS;   // (Eval's scratch per friction group, capacity reused)
      bool LevelOn(const Item &L) const { return set.levelStiffness && set.patchPath && L.ktb <= 0.0; }   // (a bend across in series: Kt is no power of u, today's mean)
      // the stiffnesses of a level form lev (u^p summed with roles): Kt = ktPerA sqrt(reff) u^0.5, Kr = kH reff u^1.5, Ktw = 16/3 0.35 0.75 Kr (patch_path's ktU, krU, kwU)
      static void LevK(const Item &L, const Item::Pows &lev, double &kt, double &kr, double &kw)
      {
         const double ktU = L.ktPerA * std::sqrt(L.reff), krU = L.kH * L.reff, kwU = 16.0 / 3.0 * 0.35 * 0.75 * krU;
         kt = ktU * lev.p05;
         kr = krU * lev.p15;
         kw = kwU * lev.p15;
      }
      // MergeBars' level form on n pieces (paper/p2b_avg.py roles()): the pressed (d > 0) taken deepest first (ties on index); a piece that touches two groups joins them at its own
      // depth, the shallower group's peak passing there (a group's peak is its root). role[i]: +1 a peak, -1 per pass, 0 none; a peak passing at its own depth adds nothing.
      // adj(i, j): the pieces' merge-tree neighbours. Scratch kept per thread (no allocation per step once grown)
      template <class Adj> static void LevelRoles(const int n, const double *d, const Adj &adj, int *role)
      {
         thread_local std::vector<int> ord, parent, passat;
         thread_local std::vector<unsigned char> onf, passed;
         ord.clear();
         for (int i = 0; i < n; i++)
         {
            role[i] = 0;
            if (d[i] > 0.0)
               ord.push_back(i);
         }
         const auto deeper = [d](const int a, const int b) { return d[a] != d[b] ? d[a] > d[b] : a < b; };
         std::sort(ord.begin(), ord.end(), deeper);
         parent.resize(n);
         passat.assign(n, -1);
         onf.assign(n, 0);
         passed.assign(n, 0);
         for (int i = 0; i < n; i++)
            parent[i] = i;
         const auto find = [&](int x) {
            while (parent[x] != x)
            {
               parent[x] = parent[parent[x]];
               x = parent[x];
            }
            return x;
         };
         for (const int x : ord)
         {
            onf[x] = 1;
            for (int y = 0; y < n; y++)
            {
               if (y == x || !onf[y] || !adj(x, y))
                  continue;
               const int rx = find(x), ry = find(y);
               if (rx == ry)
                  continue;
               int keep = rx, lose = ry;
               if (deeper(ry, rx))
                  std::swap(keep, lose);
               passed[lose] = 1;
               passat[lose] = x;
               parent[lose] = keep;
            }
         }
         for (const int i : ord)
         {
            if (passed[i] && d[passat[i]] == d[i])
               continue;
            if (!passed[i] && find(i) != i)
               continue;
            role[i]++;
            if (passed[i])
               role[passat[i]]--;
         }
      }
      // the pieces' merge-tree neighbours from their ids (RubberContact::PieceData's pair)
      static bool PieceAdj(const int ia, const RubberContact::PieceData &a, const int ib, const RubberContact::PieceData &b)
      {
         return a.nbA == ib || a.nbB == ib || b.nbA == ia || b.nbB == ia;
      }
      // the level form of one patch P of h over its own pieces at their depths (GroupMemory's weights)
      Item::Pows PatchLevel(const Item &L, const Hit &h, const HitPatch &P) const
      {
         thread_local std::vector<double> dd;
         thread_local std::vector<int> role;
         const int n = P.pn;
         dd.resize(n);
         role.resize(n);
         for (int q = 0; q < n; q++)
            dd[q] = h.pd[P.p0 + q].d;
         LevelRoles(n, dd.data(), [&](const int x, const int y) { return PieceAdj(h.pc[P.p0 + x], h.pd[P.p0 + x], h.pc[P.p0 + y], h.pd[P.p0 + y]); }, role.data());
         Item::Pows lev;
         for (int q = 0; q < n; q++)
            if (role[q] != 0)
               lev = lev + (double)role[q] * L.PowRate(dd[q]);
         return lev;
      }
      // B-lin for friction group g (FrictionGroups' gStart, fgEnd): its pieces are those of its start and end patches; each piece's start and end depth (pressed there: its probe's;
      // not pressed at that end: its own signed geometry there (PieceGeo) from the end where it is pressed), straight between; the step cut at every zero crossing and pairwise
      // crossing in (0, 1) (paper cuts_of); per segment the roles at its middle (LevelRoles), each piece's MeanOver on its segment weighted by the segment's length; the slope from
      // MeanSlope at both segment ends (dl/dd1 = tb at the end, ta at the start: MeanOver is symmetric, so the start's is MeanSlope with the ends swapped) times the piece's end
      // gradient -n_i, no cut-time terms (they cancel: paper/p2b_avg2.py, finite differences to 2e-7). One piece, or every segment one peak on the same piece: off
      void BLinGroup(const Item &L, const Hit &h0, const Hit &h1, const int g, const V3 &c0, const V3 &c1, BLin &B) const
      {
         B = BLin();
         thread_local std::vector<int> pid, role;
         thread_local std::vector<double> d0, d1, la, lb, lm, ts;
         thread_local std::vector<V3> nE;
         thread_local std::vector<RubberContact::PieceData> pdS, pdE;
         thread_local std::vector<unsigned char> on0, on1;
         pid.clear();
         pdS.clear();
         pdE.clear();
         on0.clear();
         on1.clear();
         const auto add = [&](const Hit &h, const HitPatch &P, const bool end) {
            for (int q = 0; q < P.pn; q++)
            {
               const int id = h.pc[P.p0 + q];
               int x = -1;
               for (int y = 0; y < (int)pid.size() && x < 0; y++)
                  if (pid[y] == id)
                     x = y;
               if (x < 0)
               {
                  x = (int)pid.size();
                  pid.push_back(id);
                  pdS.push_back(RubberContact::PieceData());
                  pdE.push_back(RubberContact::PieceData());
                  on0.push_back(0);
                  on1.push_back(0);
               }
               (end ? pdE : pdS)[x] = h.pd[P.p0 + q];
               (end ? on1 : on0)[x] = 1;
            }
         };
         for (int i = 0; i < h0.npat; i++)
            if (gStart[i] == g)
               add(h0, h0.pat[i], false);
         for (int j = 0; j < h1.npat; j++)
            if (fgEnd[j] == g)
               add(h1, h1.pat[j], true);
         const int n = (int)pid.size();
         if (n < 2)
            return;
         // a piece's signed depth at c from its depth at cR (its geometry; grad: its normal at c, the direction its depth falls along)
         // a face (kind 3, B2): its own triangle's closest point (the probe's ClosestOnTriangle), exact also where that point is on its rim; kind 0 (a plane) only for a probe piece with no geometry
         const auto nearOf = [](const RubberContact::PieceGeo &g, const V3 &c) {
            if (g.kind == 1)
               return RubberContact::SegNearest(g.a, g.b, c);
            if (g.kind == 3)
            {
               RubberContact::SurfaceTri t {};
               t.a = g.a;
               t.b = g.b;
               t.c = g.c;
               int kind;
               return RubberContact::ClosestOnTriangle(t, c, kind);
            }
            return g.a;
         };
         const auto geoDepth = [&nearOf](const RubberContact::PieceData &p, const V3 &cR, const V3 &c, V3 &grad) {
            if (p.g.kind == 0)
            {
               grad = p.n;
               return p.d - p.n.Dot(c - cR);
            }
            const V3 qR = nearOf(p.g, cR), q = nearOf(p.g, c);
            const double lR = (cR - qR).Length(), l = (c - q).Length();
            grad = l > 0.0 ? (c - q) * (1.0 / l) : p.n;
            return p.d + lR - l;
         };
         d0.resize(n);
         d1.resize(n);
         nE.resize(n);
         for (int x = 0; x < n; x++)
         {
            V3 gs;
            d0[x] = on0[x] ? pdS[x].d : geoDepth(pdE[x], c1, c0, gs);
            if (on1[x])
            {
               d1[x] = pdE[x].d;
               nE[x] = pdE[x].n;
            }
            else
               d1[x] = geoDepth(pdS[x], c0, c1, nE[x]);
         }
         const double *dCur = nullptr;   // (the depths the roles are taken at: a joint-joint link holds only while the segment between has no pressed foot, as BarsTube's)
         const auto adj = [&](const int x, const int y) {
            const RubberContact::PieceData &a = on1[x] ? pdE[x] : pdS[x], &b = on1[y] ? pdE[y] : pdS[y];
            if (PieceAdj(pid[x], a, pid[y], b))
               return true;
            int seg = RubberContact::PatchPieceList::kNoPiece;   // two consecutive tube joints: linked when segment seg between them is not pressed now
            if (a.nbC == pid[y])
               seg = a.nbA;
            else if (a.nbD == pid[y])
               seg = a.nbB;
            else
               return false;
            for (int z = 0; z < n; z++)
               if (pid[z] == seg && dCur != nullptr && dCur[z] > 0.0)
                  return false;
            return true;
         };
         ts.clear();   // the cuts (paper cuts_of)
         ts.push_back(0.0);
         ts.push_back(1.0);
         for (int i = 0; i < n; i++)
         {
            const double si = d1[i] - d0[i];
            if (si != 0.0)
            {
               const double t = -d0[i] / si;
               if (t > 0.0 && t < 1.0)
                  ts.push_back(t);
            }
            for (int j = i + 1; j < n; j++)
            {
               const double sj = d1[j] - d0[j];
               if (si != sj)
               {
                  const double t = (d0[j] - d0[i]) / (si - sj);
                  if (t > 0.0 && t < 1.0)
                     ts.push_back(t);
               }
            }
         }
         std::sort(ts.begin(), ts.end());
         const auto lin = [&](std::vector<double> &l, const double t) {
            l.resize(n);
            for (int x = 0; x < n; x++)
               l[x] = t == 1.0 ? d1[x] : t == 0.0 ? d0[x] : d0[x] + t * (d1[x] - d0[x]);
         };
         role.resize(n);
         const auto rate = [&L](const double d) { return L.PowRate(d); };
         const auto integral = [&L](const double d) { return L.PowIntegral(d); };
         const auto rateSlope = [&L](const double d) { return L.PowSlope(d); };
         int onePeak = -1;   // the one peak's piece while every segment has just it (-2: more, or a change)
         for (size_t s = 0; s + 1 < ts.size(); s++)
         {
            const double ta = ts[s], tb = ts[s + 1];
            if (!(tb > ta))
               continue;
            lin(la, ta);
            lin(lb, tb);
            lin(lm, 0.5 * (ta + tb));
            dCur = lm.data();
            LevelRoles(n, lm.data(), adj, role.data());
            int cnt = 0, top = -1;
            for (int x = 0; x < n; x++)
            {
               if (role[x] == 0)
                  continue;
               cnt++;
               top = role[x] == 1 ? x : -2;
               const double v = (double)role[x], w = tb - ta;
               const Item::Pows mo = Item::MeanOver(la[x], lb[x], rate, integral);
               B.m = B.m + (w * v) * mo;
               const Item::Pows sl = tb * Item::MeanSlope(la[x], lb[x], rate, rateSlope, mo) + ta * Item::MeanSlope(lb[x], la[x], rate, rateSlope, mo);
               const V3 gE = nE[x] * (-w * v);
               B.s05 = B.s05 + gE * sl.p05;
               B.s025 = B.s025 + gE * sl.p025;
               B.s15 = B.s15 + gE * sl.p15;
               B.s075 = B.s075 + gE * sl.p075;
            }
            if (cnt == 0)   // (nothing pressed on this segment: no role)
               continue;
            if (cnt != 1 || top < 0 || (onePeak != -1 && onePeak != top))
               onePeak = -2;
            else if (onePeak == -1)
               onePeak = top;
         }
         if (trLines < 400 && PatchTraceWin(pfMs))   // (VPX_PATCH_TRACE: PTBL, each group with 2+ pieces: per piece id, pressed at start / end, d0, d1 um, role at t = 0 / 1,
         {                                          //  and a face unpressed at one end: the old plane depth there beside the exact one; blin 1 B-lin on, 0 today's PatchMeans)
            trLines++;
            thread_local std::vector<int> r0, r1;
            r0.resize(n);
            r1.resize(n);
            dCur = d0.data();
            LevelRoles(n, d0.data(), adj, r0.data());
            dCur = d1.data();
            LevelRoles(n, d1.data(), adj, r1.data());
            std::ostringstream o;
            o << "PTBL ms " << pfMs << " evc " << evalCount << " g " << g << " blin " << (onePeak == -2 ? 1 : 0) << " |";
            for (int x = 0; x < n; x++)
            {
               o << " [" << pid[x] << " on " << (int)on0[x] << (int)on1[x] << " d " << d0[x] * 1e6 << " " << d1[x] * 1e6 << " role " << r0[x] << " " << r1[x];
               const RubberContact::PieceData &p = on0[x] ? pdS[x] : pdE[x];
               if (p.g.kind == 3 && on0[x] != on1[x])
                  o << " plane " << (on0[x] ? p.d - p.n.Dot(c1 - c0) : p.d - p.n.Dot(c0 - c1)) * 1e6;
               o << "]";
            }
            PLOGI << o.str();
         }
         if (onePeak != -2)   // one peak, no role change (or none): today's PatchMeans (FrictionGroups' d0, d1)
            return;
         for (int e = 0; e < 2; e++)   // the level form at the start and at the end (t = 0, 1)
         {
            const std::vector<double> &l = e == 0 ? d0 : d1;
            dCur = l.data();
            LevelRoles(n, l.data(), adj, role.data());
            Item::Pows lev;
            for (int x = 0; x < n; x++)
               if (role[x] != 0)
                  lev = lev + (double)role[x] * L.PowRate(l[x]);
            (e == 0 ? B.lev0 : B.lev1) = lev;
         }
         B.on = true;
      }
      // patch_friction: group g's memory from its start patches' pieces (plan 2.1, 2.4): each patch's (the first of its pieces found in it.pm), turned into the group's plane
      // (length kept); several by the parallel rule (k1 x1 + k2 x2) / (k1 + k2) per spring, k at the patch's start depth. None found: no memory
      void GroupMemory(const Item &it, const Item &L, const Hit &h0, const int g, FG &G, const int k) const
      {
         G.touching = false;
         G.xi = V3();
         G.ra = RubberContact::RollAngle();
         G.trFound = -2;
         if (!it.touching)
            return;
         double wt = 0.0, wr = 0.0, ww = 0.0, stw = 0.0;
         V3 sx, sth;
         int found = 0;
         const RubberContact::PieceMem *one = nullptr;
         const auto turn = [&G](const V3 &v) {
            const double l0 = v.Length();
            const V3 p = v - G.n * v.Dot(G.n);
            const double lp = p.Length();
            return lp > 1e-300 ? p * (l0 / lp) : V3();
         };
         for (int i = 0; i < h0.npat; i++)
         {
            if (gStart[i] != g)
               continue;
            const HitPatch &P = h0.pat[i];
            const RubberContact::PieceMem *m = nullptr;
            for (int q = 0; q < P.pn && m == nullptr; q++)
               for (const RubberContact::PieceMem &pm : it.pm)
                  if (pm.id == h0.pc[P.p0 + q])
                  {
                     m = &pm;
                     break;
                  }
            if (m == nullptr)
            {
               if (trLines < 400 && PatchTraceWin(pfMs))   // (VPX_PATCH_TRACE: a start patch with none of its pieces in it.pm)
               {
                  trLines++;
                  PLOGI << "PTMISS ms " << pfMs << " evc " << evalCount << " att " << trAtt << " rnd " << trRound << " fd " << trFd << " a " << it.a << " s " << it.s << " k " << k
                        << " g " << g << " patch " << i << "/" << h0.npat << " d " << P.d << " ids [" << TraceIds(h0, P) << "] pm [" << TracePmIds(it.pm) << "] pmMulti " << it.pmMulti;
               }
               continue;
            }
            bool seen = false;   // (plan 2.4: one stored group read once, however many of its start patches find it)
            for (const GmRead &r : gmRead)
               seen = seen || (r.g == g && r.sg == m->sg);
            if (seen)
               continue;
            gmRead.push_back({ g, m->sg, m->kt, m->kr, m->kw });
            found++;
            one = m;
            const double kt = m->kt, kr = m->kr, kw = m->kw;   // plan 2.4: the parallel rule's weights are the stiffness each was stored at
            sx = sx + turn(m->xi) * kt;
            sth = sth + turn(m->ra.th) * kr;
            stw += m->ra.tw * kw;
            wt += kt;
            wr += kr;
            ww += kw;
         }
         G.trFound = found;
         if (found == 0)
            return;
         G.touching = true;
         if (found == 1)
         {
            G.xi = one->xi;
            G.ra = one->ra;
            return;
         }
         G.xi = wt > 0.0 ? sx * (1.0 / wt) : V3();
         G.ra.th = wr > 0.0 ? sth * (1.0 / wr) : V3();
         G.ra.tw = ww > 0.0 ? stw / ww : 0.0;
      }
      std::vector<int> gStart;   // patch_friction: each start patch's group (FrictionGroups)
      struct GmRead   // plan 2.4: a new friction group g reads stored group sg, stored at kt, kr, kw (GroupMemory)
      {
         int g, sg;
         double kt, kr, kw;
      };
      mutable std::vector<GmRead> gmRead;   // (Eval's scratch, capacity reused)

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
            dc[i] = ball[i].held ? ball[i].Given(h) : V3(xv[6 * i], xv[6 * i + 1], xv[6 * i + 2]);   // (a held ball: its given move, spin 0)
            dw[i] = ball[i].held ? V3() : V3(xv[6 * i + 3], xv[6 * i + 4], xv[6 * i + 5]);
            wm[i] = ball[i].w + dw[i] * 0.5;
            b1[i] = bm[i] = ball[i];
            b1[i].c = ball[i].c + dc[i];
            bm[i].c = ball[i].c + dc[i] * 0.5;
            F[i] = ball[i].held ? V3() : set.g * ball[i].m + ball[i].fExt;   // (a held ball: no gravity; F gathers the contacts' forces on it)
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
         struct PzEvalScope   // patch_zero (F6): this Eval's probes keep a net (PzNetCheck); none from an earlier Eval
         {
            const Group &g;
            explicit PzEvalScope(const Group &gs) : g(gs) { g.pzEvalOn = true; g.pzNetS = -1; }
            ~PzEvalScope() { g.pzEvalOn = false; }
         } pzScope(*this);
         if (set.linkWindows)   // link_windows: this Eval's order check and added links per item (LinkUnion writes a gated item's)
         {
            lwCut.assign((size_t)ni, -1.0);
            pzMg.resize((size_t)ni);   // (patch_zero, F4: LinkUnion fills a gated item's)
            for (std::vector<std::pair<int, int>> &g : pzMg)
               g.clear();
            lwAdd.assign((size_t)ni, 0);
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
            for (int p = 0; anyMassless && p < np; p++)
            {
               mlWork[p] = 0.0;
               for (int u = 0; J && u < nU; u++)
                  mlWg[(size_t)p * jStride + u] = 0.0;
            }
         }
         driveW = 0.0;
         for (int i = 0; i < nb; i++)
            driveW += ball[i].fExt.Dot(dc[i]);
         // the place-dependent forces (BallForce): read at the midpoint, a drive as fExt is (its work the midpoint force times the move); its slope in the Jacobian
         for (int i = 0; nbf > 0 && i < nb; i++)
            if (ball[i].forceScale != 0.0)
            {
               M3 sf;
               const V3 f = BallForceAt(ball[i], bm[i].c, J ? &sf : nullptr);
               F[i] = F[i] + f;
               driveW += f.Dot(dc[i]);
               if (impOn)
                  impX[i] = impX[i] + f;   // (test tooling: booked with fExt, the held drive it replaces)
               if (J)
                  AddBlock(J, 6 * i, 6 * i, sf, -0.5 * h);   // (R = inertia - F h; the midpoint moves half as far as the end)
            }
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
            const double iM[6] = { ball[T.a].InvM(), ball[T.a].InvM(), ball[T.a].InvM(), two ? ball[T.b].InvM() : 0.0, two ? ball[T.b].InvM() : 0.0, two ? ball[T.b].InvM() : 0.0 };
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
            e.ut2 = e.ftTr = e.lim = 0.0;   // (diagnostics: this Eval's own)
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
            // [Contact] the elastic seam (multi-contact piece 1): the elastic force, its gradients and its slope come from the law's Item::Elastic; here only the move's
            // coordinates (ball a's centre, then ball b's, or the columns' coordinates) and their inverse masses
            double gm[kGradMax], dq[kGradMax], iM[kGradMax], gd[kGradMax];
            const Ball &A = ball[it.a];
            dq[0] = dc[it.a].x; dq[1] = dc[it.a].y; dq[2] = dc[it.a].z;
            iM[0] = iM[1] = iM[2] = 1.0 / A.m;
            if (pr)
            {
               dq[3] = dc[it.b].x; dq[4] = dc[it.b].y; dq[5] = dc[it.b].z;
               iM[3] = iM[4] = iM[5] = ball[it.b].InvM();   // (a held ball takes none of the correction: infinite mass)
            }
            for (int j = 0; j < nq; j++)
            {
               const int cj = qc0 + hm.c[j].coord;
               dq[3 + j] = qEnd[cj];
               iM[3 + j] = cUnk[cj] >= 0 ? ci[cj].invM : 0.0;   // (a move given to it takes none of the correction: only the unknowns do)
            }
            const Hit &h0 = h0s[k];
            std::vector<std::array<double, 4>> &tg = tgS;
            std::vector<double> &tk = tkS;
            std::vector<char> &tOn = tOnS;
            std::vector<std::array<double, 4>> &tgE = tgES;   // each term's end gradient (endGradSlope; else = tg): the column side of its slope
            std::vector<std::array<double, 4>> &tgF = tgFS;   // each term's force direction gdT (termForceRows; else = tg): the row side of its slope
            const Item::ElasticEnd el = L.Elastic(h0, h1, hm, d0, dq, iM, na, nq, noiseD, set, clampDump, J ? true : false, gm, gd, tg, tk, tOn, tgE, tgF);
            const bool terms = el.terms, cf0U = el.cf0U;
            const uint64_t tCf0 = el.tCf0;
            const double dd = el.dd, feT = el.feT, fe = el.fe;
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
            // [Contact] patch_friction (piece 2b-1; plan runs/2026-11-22-per-peak 2.1, 2.4, 2.5b): a static mesh's or tube's friction, damping, rolling and twist per friction group
            // (FrictionGroups: the start and end patches linked by shared pieces), each with its own normal, depths, load fe + fd and memory (GroupMemory). One group: today's inputs,
            // the same arithmetic: the block below runs once per group and keeps its results in fgS for the slope. The push (gd N) stays the merged contact's. Not in this build:
            // the B-lin stiffness average (a group's stiffness is PatchMeans on its own one depth), ridge arrival, turn terms
            const bool pfPath = !pr && nq == 0 && !terms && (h0.npat > 0 || h1.npat > 0);
            const int ngF = pfPath ? FrictionGroups(h0, h1, hm, A.c, b1[it.a].c, bm[it.a].c, k) : 0;
            const bool multiF = ngF > 1;
            const int ngLoop = multiF ? ngF : 1;
            if (fgS.empty())
               fgS.resize(1);
            if ((int)blS.size() < std::max(ngF, 1))
               blS.resize(std::max(ngF, 1));
            for (int g = 0; g < std::max(ngF, 1); g++)   // level_stiffness: each group's B-lin (one group too: a split happens inside a one-group step)
               blS[g].on = false;
            if (ngF > 0 && LevelOn(L))
               for (int g = 0; g < ngF; g++)
                  BLinGroup(L, h0, h1, g, A.c, b1[it.a].c, blS[g]);
            gmRead.clear();
            if (!multiF)
            {
               FG &G = fgS[0];
               G = FG();
               G.n = nm;
               G.nD = nm;
               G.d0 = d0;
               G.d1 = d1;
               G.dm = dm;
               G.fe = fe;
               G.fd = fd;
               G.N = N;
               G.aE = aE;
               G.dadfE = dadfE;
               G.touching = it.touching;
               G.xi = it.xi;
               G.ra = it.ra;
               if (ngF == 1 && it.pmMulti)   // the last step had 2 or more groups: the memory from the pieces (a join: the parallel rule)
                  GroupMemory(it, L, h0, 0, G, k);
            }
            else
               for (int g = 0; g < ngF; g++)
               {
                  GroupLoad(L, h, fgS[g]);
                  GroupMemory(it, L, h0, g, fgS[g], k);
               }
            if (pfPath && ngF > 0)   // plan 2.4's memory rule: each group's stored stiffness ks (the patch rule from ks to its start stiffness k0 runs in the block below)
            {
               for (int g = 0; g < ngLoop; g++)
               {
                  FG &Gg = fgS[g];
                  if (blS[g].on)
                     LevK(L, blS[g].lev0, Gg.k0p, Gg.kr0p, Gg.kw0p);
                  else
                  {
                     Gg.k0p = L.Kt(Gg.d0);
                     Gg.kr0p = L.Kr(Gg.d0);
                     Gg.kw0p = L.Ktw(Gg.d0);
                  }
               }
               if (!multiF && !(ngF == 1 && it.pmMulti))   // the item's one stretch: stored at kt1, kr1, kw1
               {
                  if (it.kStored && fgS[0].touching)
                  {
                     fgS[0].ks = it.kt1;
                     fgS[0].krs = it.kr1;
                     fgS[0].kws = it.kw1;
                  }
               }
               else   // read through the pieces: each stored group's k shared among the new groups that read it, in proportion to their k0 (a split); summed over the stored groups (a join)
                  for (const GmRead &r : gmRead)
                  {
                     double s0 = 0.0, s1 = 0.0, s2 = 0.0;
                     for (const GmRead &q : gmRead)
                        if (q.sg == r.sg)
                        {
                           s0 += fgS[q.g].k0p;
                           s1 += fgS[q.g].kr0p;
                           s2 += fgS[q.g].kw0p;
                        }
                     FG &Gg = fgS[r.g];
                     Gg.ks += s0 > 0.0 ? r.kt * Gg.k0p / s0 : 0.0;
                     Gg.krs += s1 > 0.0 ? r.kr * Gg.kr0p / s1 : 0.0;
                     Gg.kws += s2 > 0.0 ? r.kw * Gg.kw0p / s2 : 0.0;
                  }
            }
            F[it.a] = F[it.a] - V3(gd[0], gd[1], gd[2]) * N;   // the push (a massless part's contact: its spring is in the part's U*, only the damping here)
            if (pr)
               F[it.b] = F[it.b] - V3(gd[3], gd[4], gd[5]) * N;
            double accSlip = 0.0, accShear = 0.0, accRoll = 0.0, accTwist = 0.0;
            double qf[7];   // (a carried surface's coordinates: one group only)
            for (int gF = 0; gF < ngLoop; gF++)
            {
            FG &GF = fgS[gF];
            const V3 &nm = GF.n;
            const double dm = GF.dm, d0 = GF.d0, d1 = GF.d1, fe = GF.fe, fd = GF.fd, N = GF.N, aE = GF.aE, dadfE = GF.dadfE;
            if (multiF)
               e.slip = e.shearDamp = e.roll = e.twist = 0.0;

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
            V3 xi0 = GF.touching ? GF.xi : V3();
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
            const BLin *const BL = blS[gF].on ? &blS[gF] : nullptr;   // level_stiffness: the group's means from B-lin (its slopes are vectors: BL->s05 .., kSp, ctSp, twD1 stay 0)
            Item::Pows pM, pS;
            if (BL != nullptr)
               pM = BL->m;
            else if (pp)
               L.PatchMeans(d0, d1, pM, pS);
            const double ktU = pp ? L.ktPerA * std::sqrt(L.reff) : 0.0, ctU = pp ? L.shearTan * std::sqrt(ktU * L.mt) : 0.0;   // Kt and its damper per u^0.5, u^0.25
            double k0 = 0.0, k1 = 0.0, kr0L = 0.0, kr1L = 0.0, kw0L = 0.0, kw1L = 0.0;   // (B-lin: the level form's at the start and the end, shear, rolling, twist)
            if (BL != nullptr)
            {
               LevK(L, BL->lev0, k0, kr0L, kw0L);
               LevK(L, BL->lev1, k1, kr1L, kw1L);
            }
            else
            {
               k0 = L.Kt(d0);
               k1 = L.Kt(d1);
            }
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
            if (GF.ks > 0.0 && k0 > GF.ks)   // plan 2.4's memory rule: from the stiffness it was stored at to k0 (grows: holds the force; shrinks: holds the stretch), a loss (booked below)
               xi0 = xi0 * (GF.ks / k0);
            const double held = 0.5 * k0 * xi0.Dot(xi0);
            GF.heldS = held;
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
               e.ut2 = ut.Dot(ut);   // (diagnostics, copies only: the halved line's contact states)
               e.ftTr = fl;
               e.lim = lim;
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
            F[it.a] = F[it.a] - Ft;
            M[it.a] = M[it.a] + rA.Cross(Ft * -1.0);
            if (pr)
            {
               F[it.b] = F[it.b] + Ft;
               M[it.b] = M[it.b] + rB.Cross(Ft);
            }
            // on a carried surface's coordinates: the same work with the other sign
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
            V3 twD1V;   // (B-lin: dMw/dc1 sticking, a vector; twD1 0)
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
               RubberContact::RollAngle r0 = GF.touching ? GF.ra : RubberContact::RollAngle();
               {
                  const double l0 = r0.th.Length();
                  const V3 tp = r0.th - nm * r0.th.Dot(nm);
                  const double lp = tp.Length();
                  r0.th = lp > 1e-300 ? tp * (l0 / lp) : V3();
               }
               // patch_path: the step's springs krS, kwS and dampers are their means over the depth path (u^1.5, u^0.75), the patch rule from the start's to them here and
               // from them to the end's after each move; 0: krS = kr1, kwS = kw1, the dampers from them
               const double kr0 = BL != nullptr ? kr0L : L.Kr(d0), kr1 = BL != nullptr ? kr1L : L.Kr(d1), kw0 = BL != nullptr ? kw0L : L.Ktw(d0), kw1 = BL != nullptr ? kw1L : L.Ktw(d1);   // (B-lin: the level form's)
               const double Ir = L.mRed * A.R * A.R;   // with 1.4: the ball's inertia about its contact point; with 0.4: about its centre (the contact's mass)
               const double krU = pp ? L.kH * L.reff : 0.0, kwU = 16.0 / 3.0 * 0.35 * 0.75 * krU;   // Kr, Ktw per u^1.5 (a^3 = reff^1.5 u^1.5)
               const double krS = pp ? krU * pM.p15 : kr1, kwS = pp ? kwU * pM.p15 : kw1;
               if (GF.krs > 0.0 && kr0 > GF.krs)   // plan 2.4's memory rule (rolling and twist), as the shear's above
                  r0.th = r0.th * (GF.krs / kr0);
               if (GF.kws > 0.0 && kw0 > GF.kws)
                  r0.tw *= GF.kws / kw0;
               const double held = 0.5 * kr0 * r0.th.Dot(r0.th) + 0.5 * kw0 * r0.tw * r0.tw;
               GF.heldR = held;
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
                        if (BL != nullptr)   // level_stiffness: the means' slopes in the end position (B-lin's vectors)
                           twD1V = BL->s15 * (kwU * ((kwS > kw0 ? 0.0 : r0.tw) + 0.5 * h * wn)) + BL->s075 * (cwU * wn);
                        else if (pp)   // the means' exact slopes in d1
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
            GF.Ft = Ft;
            GF.xi0 = xi0;
            GF.xi1 = xi1;
            GF.ut = ut;
            GF.FsTr = FsTr;
            GF.rA = rA;
            GF.rB = rB;
            GF.wr = wr;
            GF.wtR = wtR;
            GF.mhR = mhR;
            GF.RmS = RmS;
            GF.ra1 = e.ra1;
            GF.kS = kS;
            GF.kSp = kSp;
            GF.ctSp = ctSp;
            GF.k0 = k0;
            GF.k1 = k1;
            GF.ct = ct;
            GF.limS = limS;
            GF.flTr = flTr;
            GF.kTw = kTw;
            GF.kkTw = kkTw;
            GF.sRo = sRo;
            GF.rhoR = rhoR;
            GF.wlR = wlR;
            GF.twD1 = twD1;
            GF.blin = BL != nullptr;
            GF.kSpV = BL != nullptr ? BL->s05 * ktU : V3();
            GF.ctSpV = BL != nullptr ? BL->s025 * ctU : V3();
            GF.twD1V = twD1V;
            GF.kr1L = kr1L;
            GF.kw1L = kw1L;
            GF.twA = twA;
            GF.twN = twN;
            GF.twF = twF;
            GF.limN = limN;
            GF.limF = limF;
            GF.twM = twM;
            GF.twLim = twLim;
            GF.twSt = twSt;
            GF.slipping = slipping;
            GF.rollS = rollS;
            if (multiF)
            {
               accSlip += e.slip;
               accShear += e.shearDamp;
               accRoll += e.roll;
               accTwist += e.twist;
            }
            }
            if (pmEnd.size() < (size_t)ni)
            {
               pmEnd.resize(ni);
               pmEndOn.resize(ni, 0);
            }
            pmEnd[k].clear();
            pmEndOn[k] = (signed char)(ngF == 0 ? 0 : multiF ? 2 : 1);
            if (pmEndJoin.size() < (size_t)ni)
            {
               pmEndJoin.resize(ni, 0);
               pmEndSplit.resize(ni, 0);
            }
            pmEndJoin[k] = (signed char)(ngF > 0 ? std::min(fgJoins, 127) : 0);
            pmEndSplit[k] = (signed char)(ngF > 0 ? std::min(fgSplits, 127) : 0);
            if (ngF > 0 && trLines < 400 && PatchTraceWin(pfMs))   // (VPX_PATCH_TRACE: this Eval's friction groups)
            {
               trLines++;
               char b[512];
               snprintf(b, sizeof b, "PTEV ms %u evc %ld att %ld rnd %d fd %d h %.6g us k %d a %d s %d c %.9g %.9g %.9g pm0 %d pmMulti %d n0 %d n1 %d nm %d ngF %d",
                        pfMs, evalCount, trAtt, trRound, trFd, h * 1e6, k, it.a, it.s, b1[it.a].c.x, b1[it.a].c.y, b1[it.a].c.z, (int)it.pm.size(), (int)it.pmMulti,
                        h0.npat, h1.npat, hm.npat, ngF);
               std::string s = b;
               for (int g = 0; g < ngF; g++)
               {
                  const FG &Gg = fgS[g];
                  std::string si, se;
                  for (int i = 0; i < h0.npat; i++)
                     if (gStart[i] == g)
                        si += "[" + TraceIds(h0, h0.pat[i]) + "]";
                  for (int j = 0; j < h1.npat; j++)
                     if (fgEnd[j] == g)
                        se += "[" + TraceIds(h1, h1.pat[j]) + "]";
                  snprintf(b, sizeof b, " | g%d s%s e%s d0 %.6g dm %.6g d1 %.6g cmn %d N %.6g %s found %d xi0 %.6g xi1 %.6g n %.6f %.6f %.6f", g, si.c_str(), se.c_str(),
                           Gg.d0 * 1e6, Gg.dm * 1e6, Gg.d1 * 1e6, Gg.trCmn, Gg.N, Gg.slipping ? "slip" : "stick", ngF == 1 && !it.pmMulti ? 9 : Gg.trFound,   // (found 9: one group, memory from it.xi)
                           Gg.xi0.Length() * 1e6, Gg.xi1.Length() * 1e6, Gg.n.x, Gg.n.y, Gg.n.z);
                  s += b;
               }
               PLOGI << s;
            }
            if (ngF > 0 && it.touching)   // patch_friction (plan 2.4): the memory read (the item's or through the pieces, the stored-stiffness rule, a join's parallel rule, a group that is gone) (GroupMemory: a join's parallel rule, a group that is gone): the springs'
            {                                       // energy the book holds at the start (ShearE, RollE) less what the groups hold after the rule is the rule's loss, booked in this step (>= 0)
               double hS = 0.0, hR = 0.0;
               for (int g = 0; g < ngLoop; g++)
               {
                  hS += fgS[g].heldS;
                  hR += fgS[g].heldR;
               }
               const double lS = (it.touching ? it.ShearE(d0) : 0.0) - hS, lR = (it.touching && set.rollModel == 1 ? it.RollE(d0) : 0.0) - hR;
               e.memLossS = lS > 0.0 ? lS : 0.0;   // (a negative one is not booked: it shows in the held check)
               e.memLossR = lR > 0.0 ? lR : 0.0;
            }
            if (pfPath && RubberContact::PatchCheckEnv() && trFd == 0)   // VPX_PATCH_CHECK: the groups' held energy at the start vs the book's (Energy: ShearE, Kr / Ktw at the depth), once per attempt
            {
               if ((int)chkAtt.size() < ni)
                  chkAtt.resize(ni, -1);
               if (chkAtt[k] != trAtt)
               {
                  chkAtt[k] = trAtt;
                  double heldSum = 0.0;
                  for (int g = 0; g < ngLoop; g++)
                     heldSum += fgS[g].heldS + fgS[g].heldR;
                  const double bookS = it.touching ? it.ShearE(d0) : 0.0;
                  const double bookR = it.touching && set.rollModel == 1 ? it.RollE(d0) : 0.0;
                  const double diff = heldSum + e.memLossS + e.memLossR - (bookS + bookR);   // (held plus the memory loss booked in this step)
                  RubberContact::HeldCheck()[0]++;
                  if (ngF == 1 && blS[0].on)   // lev1 of the one B-lin group against the end probe's own level form (MergeBars' peaks and passes, Hit::pk; one peak: its depth)
                  {
                     double lp = 0.0;
                     if (h1.np > 0)
                        for (int j = 0; j < h1.np; j++)
                           lp += h1.pk[j].s * L.PowRate(h1.pk[j].d).p05;
                     else
                        lp = L.PowRate(h1.delta).p05;
                     const double lb = blS[0].lev1.p05;
                     RubberContact::HeldCheck()[4]++;
                     if (std::fabs(lb - lp) > 1e-12 * std::max(std::fabs(lb), std::fabs(lp)))
                     {
                        const long long cL = ++RubberContact::HeldCheck()[5];
                        if (cL <= 20)
                           PLOGI << "LEV1CHK ms " << pfMs << " item " << k << " s " << it.s << " B-lin lev1 " << lb << " probe " << lp << " (u^0.5, m^0.5) peaks/passes " << h1.np << " d1 " << h1.delta * 1e6 << " um";
                     }
                  }
                  for (int g = 0; g < ngLoop; g++)   // kind C: a group's start stiffness against the stiffness its memory was stored at (plan 2.4)
                  {
                     const FG &Gg = fgS[g];
                     if (Gg.ks <= 0.0 || std::fabs(Gg.k0 - Gg.ks) <= 1e-9 * std::max(Gg.k0, Gg.ks))
                        continue;
                     const long long cC = ++RubberContact::HeldCheck()[3];
                     if (cC > 20)
                        continue;
                     std::ostringstream o;
                     o << "KINDC ms " << pfMs << " item " << k << " s " << it.s << " ngF " << ngF << " g " << g << " blin " << Gg.blin << " pmLev " << it.pmLev << " pmMulti " << it.pmMulti
                       << " kStored " << it.kStored << " island " << (islandRing ? "ring" : "cord") << " reff " << L.reff << " storedReff " << it.kReff << " firstStep " << (books.steps == 0)
                       << " | k0 " << Gg.k0 << " ks " << Gg.ks << " kr0 " << Gg.kr0p << " krs " << Gg.krs << " kw0 " << Gg.kw0p << " kws " << Gg.kws << " | start pieces";
                     for (int i = 0; i < h0.npat; i++)
                     {
                        if (ngF > 0 && gStart[i] != g)
                           continue;
                        const HitPatch &P = h0.pat[i];
                        thread_local std::vector<double> dd;
                        thread_local std::vector<int> rl;
                        dd.resize(P.pn);
                        rl.resize(P.pn);
                        for (int q = 0; q < P.pn; q++)
                           dd[q] = h0.pd[P.p0 + q].d;
                        LevelRoles(P.pn, dd.data(), [&](const int x, const int y) { return PieceAdj(h0.pc[P.p0 + x], h0.pd[P.p0 + x], h0.pc[P.p0 + y], h0.pd[P.p0 + y]); }, rl.data());
                        for (int q = 0; q < P.pn; q++)
                           o << " [" << h0.pc[P.p0 + q] << " " << dd[q] * 1e6 << " r" << rl[q] << "]";
                     }
                     o << " | stored (id sg kt)";
                     for (const RubberContact::PieceMem &pm : it.pm)
                        o << " [" << pm.id << " " << pm.sg << " " << pm.kt << "]";
                     PLOGI << o.str();
                  }
                  if (std::fabs(diff) > 1e-12)
                  {
                     const long long c = ++RubberContact::HeldCheck()[1];
                     if (c <= 20)
                     {
                        std::ostringstream o;
                        o << "HELDCHK ms " << pfMs << " item " << k << " s " << it.s << " ball " << it.a << " ngF " << ngF << " pmLev " << it.pmLev << " pmMulti " << it.pmMulti << " memLoss " << (e.memLossS + e.memLossR) * 1e6 << " uJ"
                          << " d0 " << d0 * 1e6 << " um, Kt(d0) " << L.Kt(d0) << " held " << heldSum * 1e6 << " uJ (shear+roll per group below), book " << (bookS + bookR) * 1e6
                          << " uJ (ShearE " << bookS * 1e6 << ", rolling " << bookR * 1e6 << ", RollE " << (it.touching ? it.RollE(d0) * 1e6 : 0.0) << "), diff " << diff * 1e6 << " uJ |";
                        for (int g = 0; g < ngLoop; g++)
                        {
                           const FG &Gg = fgS[g];
                           o << " g" << g << " blin " << Gg.blin << " k0 " << Gg.k0 << " Kt(g d0) " << L.Kt(Gg.d0) << " d0 " << Gg.d0 * 1e6 << " xi0 " << Gg.xi0.Length() * 1e6
                             << " um heldS " << Gg.heldS * 1e6 << " heldR " << Gg.heldR * 1e6 << " pieces";
                           for (int i = 0; i < h0.npat; i++)   // (its start patches' pieces: id, depth um; t = 0 roles over each patch's own pieces)
                           {
                              if (ngF > 0 && gStart[i] != g)
                                 continue;
                              const HitPatch &P = h0.pat[i];
                              thread_local std::vector<double> dd;
                              thread_local std::vector<int> rl;
                              dd.resize(P.pn);
                              rl.resize(P.pn);
                              for (int q = 0; q < P.pn; q++)
                                 dd[q] = h0.pd[P.p0 + q].d;
                              LevelRoles(P.pn, dd.data(), [&](const int x, const int y) { return PieceAdj(h0.pc[P.p0 + x], h0.pd[P.p0 + x], h0.pc[P.p0 + y], h0.pd[P.p0 + y]); }, rl.data());
                              for (int q = 0; q < P.pn; q++)
                                 o << " [" << h0.pc[P.p0 + q] << " " << dd[q] * 1e6 << " r" << rl[q] << "]";
                           }
                        }
                        PLOGI << o.str();
                     }
                  }
               }
            }
            if (multiF)
            {
               static std::atomic<long long> nEvalMulti { 0 };
               PfCount(nEvalMulti, 1, "Evals with 2 or more friction groups", pfMs);
            }
            if (ngF > 0)   // patch_friction: the end pieces with their group's memory (Accept keeps them on the item)
            {
               for (int j = 0; j < h1.npat; j++)
               {
                  const int g = fgEnd[j];
                  if (g < 0)
                     continue;
                  const HitPatch &P = h1.pat[j];
                  for (int q = 0; q < P.pn; q++)
                     pmEnd[k].push_back({ h1.pc[P.p0 + q], fgS[g].xi1, fgS[g].ra1, fgS[g].k1, fgS[g].blin ? fgS[g].kr1L : L.Kr(fgS[g].d1), fgS[g].blin ? fgS[g].kw1L : L.Ktw(fgS[g].d1), g });
               }
               if (multiF)   // the item's one stretch and angles (the energy book, a following one-group step's start): the groups' by the parallel rule at their end stiffness
               {
                  e.slip = accSlip;
                  e.shearDamp = accShear;
                  e.roll = accRoll;
                  e.twist = accTwist;
                  V3 sx, sth;
                  double wt = 0.0, wrS = 0.0, ww = 0.0, stw = 0.0;
                  for (int g = 0; g < ngF; g++)
                  {
                     const FG &Gg = fgS[g];
                     if (Gg.d1 <= 0.0)
                        continue;
                     const double kt = Gg.blin ? Gg.k1 : L.Kt(Gg.d1), kr = Gg.blin ? Gg.kr1L : L.Kr(Gg.d1), kw = Gg.blin ? Gg.kw1L : L.Ktw(Gg.d1);   // (B-lin: the level form's end stiffness, the next step's start)
                     e.pmEs += 0.5 * kt * Gg.xi1.Dot(Gg.xi1);   // (plan_p2b1 A4: the book's springs, one per group)
                     e.pmEr += 0.5 * kr * Gg.ra1.th.Dot(Gg.ra1.th) + 0.5 * kw * Gg.ra1.tw * Gg.ra1.tw;
                     sx = sx + Gg.xi1 * kt;
                     sth = sth + Gg.ra1.th * kr;
                     stw += Gg.ra1.tw * kw;
                     wt += kt;
                     wrS += kr;
                     ww += kw;
                  }
                  e.xi1 = wt > 0.0 ? sx * (1.0 / wt) : V3();
                  e.ra1.th = wrS > 0.0 ? sth * (1.0 / wrS) : V3();
                  e.ra1.tw = ww > 0.0 ? stw / ww : 0.0;
               }
               e.slip += e.memLossS;   // (the memory rule's loss, above; after the groups' sum)
               e.roll += e.memLossR;
               e.kStored = !multiF;   // (plan 2.4: the one group's end stiffness, B-lin's or the merged depth's)
               e.kt1 = fgS[0].k1;
               e.kReff = L.reff;
               e.kr1 = fgS[0].blin ? fgS[0].kr1L : L.Kr(fgS[0].d1);
               e.kw1 = fgS[0].blin ? fgS[0].kw1L : L.Ktw(fgS[0].d1);
               e.pmLev = !multiF && fgS[0].blin;   // level_stiffness, one group on B-lin: the book's springs at its level-form end stiffness (Item::pmLev), not Kt(merged depth)
               if (e.pmLev)
               {
                  const FG &G0 = fgS[0];
                  e.pmEs = 0.5 * G0.k1 * G0.xi1.Dot(G0.xi1);
                  e.pmEr = 0.5 * G0.kr1L * G0.ra1.th.Dot(G0.ra1.th) + 0.5 * G0.kw1L * G0.ra1.tw * G0.ra1.tw;
               }
            }
            // (one group, or the first: the code below until the slope reads these; the slope runs per group)
            const V3 &Ft = fgS[0].Ft, &xi0 = fgS[0].xi0, &xi1 = fgS[0].xi1, &ut = fgS[0].ut, &rA = fgS[0].rA, &rB = fgS[0].rB;
            const bool slipping = fgS[0].slipping;
            const double limS = fgS[0].limS, kS = fgS[0].kS, ct = fgS[0].ct, twM = fgS[0].twM, twLim = fgS[0].twLim;
            const int twSt = fgS[0].twSt;
            bool mlCol[7] = {};     // a massless coordinate's column (design section 3 a): its row takes the push at the step's end
            double mlSlope[7] = {};  // its end normal force's slope dN1/dd1 (N/m)
            for (int j = 0; j < nq; j++)
            {
               const int cj = qc0 + hm.c[j].coord;
               if (Massless(cj))
               {
                  // rule C (runs/2026-11-19-standup-p2/design.md section 2 C): the force balance at the step's end, so the contacts' push on it is the end probe's; the
                  // work the step's mean reaction does on it is kept for the energy correction (MasslessEnergy)
                  mlCol[j] = true;
                  mlWork[shape[it.s]->part] += qf[j] * dq[3 + j];
                  qF[cj] += MasslessPush(L, h1, j, d0, d1, h, Ft, mlSlope[j]);
               }
               else
                  qF[cj] += qf[j];
            }

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
               const double kc = el.kc;   // (the elastic seam's)
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
               // peak_slope (P2a-3): a static surface with 2 or more peaks: the push P = N gd (N = fe + fd) by the product rule, in place of today's kn gm gm^T below:
               // h dP/du = kn gd g1^T + h N dgd/du, g1 the end depth's gradient on the ball's centre; Gonzalez gd = gm + lam W u, lam = (dd - gm . u) / (u . W u), so
               // dgd/du = H / 2 + lam W + (W u) dlam^T, dlam = (g1 - gm - (H / 2) u - 2 lam W u) / (u . W u); H / 2 alone where Gonzalez left its correction out
               // (H the merged depth's Hessian, PeakHessian; paper check runs/2026-11-22-per-peak/paper/dampslope.py form A3: exact to 2e-5 of |J|)
               double pkH[3][3];
               const bool pkOn = set.peakSlope && hm.np > 0 && nq == 0 && !pr && fe > 0.0 && fe + fd > 0.0 && PeakHessian(hm.pk, hm.np, gm, pkH);
               if (terms)   // each term's own chord slope on its own gradient, the damping's on G
               {
                  if (fe + fd > 0.0)
                  {
                     const double kd = di ? dSl : L.gammaKK > 0.0 ? kkJac : cD + cDp * dd;
                     for (int p = 0; p < 4; p++)
                        if (ix[p] >= 0 && !(p >= 3 && mlCol[p - 3]))   // (a massless row: its slope below)
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
               else if (!pkOn)
                  for (int p = 0; p < na + nq; p++)
                     if (ix[p] >= 0 && !(p >= na && mlCol[p - na]))
                        for (int q = 0; q < na + nq; q++)
                           if (ix[q] >= 0)
                              J[ix[p]][ix[q]] += kn * gm[p] * gm[q];
               // peak_slope: the push's product-rule slope (above); hm is the middle probe: it moves half the move. Friction, rolling and twist unchanged.
               if (pkOn)
               {
                  const double g1[3] = { -h1.n.x * h1.grad, -h1.n.y * h1.grad, -h1.n.z * h1.grad }, N = fe + fd;
                  double W[3], Wu[3], Hu[3], uWu = 0.0, gu = 0.0;
                  for (int a = 0; a < 3; a++)
                  {
                     W[a] = !set.gonzalezMass ? iM[a] : iM[a] > 0.0 ? iM[0] / iM[a] : 0.0;   // (Gonzalez's weights)
                     Wu[a] = W[a] * dq[a];
                     gu += gm[a] * dq[a];
                     uWu += dq[a] * W[a] * dq[a];
                     Hu[a] = pkH[a][0] * dq[0] + pkH[a][1] * dq[1] + pkH[a][2] * dq[2];
                  }
                  const double lam = uWu > 1e-300 && std::fabs(d1 - d0 - gu) > noiseD ? (d1 - d0 - gu) / uWu : 0.0;   // (Gonzalez's cf, its noise test)
                  double dlam[3] = {};
                  if (lam != 0.0)
                     for (int q = 0; q < 3; q++)
                        dlam[q] = (g1[q] - gm[q] - 0.5 * Hu[q] - 2.0 * lam * Wu[q]) / uWu;
                  for (int p = 0; p < 3; p++)
                     if (ix[p] >= 0)
                        for (int q = 0; q < 3; q++)
                           if (ix[q] >= 0)
                              J[ix[p]][ix[q]] += kn * gd[p] * g1[q] + h * N * (0.5 * pkH[p][q] + (p == q ? lam * W[p] : 0.0) + Wu[p] * dlam[q]);
                  if (PeakSlopeCheckEnv())   // check mode (VPX_PATCH_CHECK=1): the list's first peak alone, unturned, gives today's slope (H = 0 to rounding)
                  {
                     static std::atomic<long long> nChk { 0 }, nBad { 0 };
                     HitPeaks one;
                     one[0] = hm.pk[0];
                     one[0].s = 1.0;
                     one[0].r = 0.0;
                     const double g1[3] = { -one[0].n.x, -one[0].n.y, -one[0].n.z };
                     double H1[3][3], w = 0.0;
                     if (PeakHessian(one, 1, g1, H1))
                        for (int p = 0; p < 3; p++)
                           for (int q = 0; q < 3; q++)
                              w = std::max(w, std::abs(H1[p][q]) * one[0].d / 1.5);   // relative to the peak's own 1.5 d^0.5 / d^1.5
                     const long long c = ++nChk, b = w > 1e-12 ? ++nBad : nBad.load();
                     if ((c & (c - 1)) == 0 || (w > 1e-12 && b <= 20))   // running totals at 1, 2, 4, ... lists (the bench never closes the table), and the first 20 misses
                        PLOGI << "peak_slope check: " << c << " lists, " << b << " with one peak not today's slope (this one " << w << ", np " << hm.np << ")";
                  }
               }
               // a massless row (rule C): -h times the end push's slope, h s1 dN1/dd1 g1 over the unknowns (g1 the end depth's gradient: the ball's centre -n1 grad1, the
               // columns' slope x grad1); the shear's push keeps its slope below
               for (int j = 0; j < nq; j++)
                  if (mlCol[j] && ix[na + j] >= 0 && mlSlope[j] > 0.0)
                  {
                     const double s1 = h1.c[j].slope * h1.grad, a = h * s1 * mlSlope[j];
                     const int qi = ix[na + j];
                     for (int q = 0; q < 3; q++)
                        if (ix[q] >= 0)
                           J[qi][ix[q]] += a * -(q == 0 ? h1.n.x : q == 1 ? h1.n.y : h1.n.z) * h1.grad;
                     for (int jj = 0; jj < nq; jj++)
                        if (ix[na + jj] >= 0)
                           J[qi][ix[na + jj]] += a * h1.c[jj].slope * h1.grad;
                  }
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
               // a massless column (MasslessEnergy's slope): its row as it stands before the shear's and the spins' terms, so their share of the mean reaction's slope is the change below
               for (int j = 0; j < nq; j++)
                  if (mlCol[j] && ix[3 + j] >= 0)
                     for (int u = 0; u < nU; u++)
                        mlSnap[(size_t)j * jStride + u] = J[ix[3 + j]][u];
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
               const double knO = kn, kcO = kc;
               const double *const gmO = gm;
               for (int gF = 0; gF < ngLoop; gF++)   // patch_friction: the friction's slope per group (one: today's inputs)
               {
               const FG &GF = fgS[gF];
               const V3 &nm = GF.n, &nD = GF.nD, &Ft = GF.Ft, &xi0 = GF.xi0, &ut = GF.ut, &FsTr = GF.FsTr, &rA = GF.rA, &rB = GF.rB, &wr = GF.wr, &wtR = GF.wtR, &mhR = GF.mhR;
               const M3 &RmS = GF.RmS;
               const double dm = GF.dm, d1 = GF.d1, N = GF.N, kn = multiF ? GF.kn : knO, kc = multiF ? GF.kc : kcO;
               const bool pp = set.patchPath;
               const double *const gm = multiF ? GF.gm : gmO;
               const double kS = GF.kS, kSp = GF.kSp, ctSp = GF.ctSp, k0 = GF.k0, k1 = GF.k1, ct = GF.ct, limS = GF.limS, flTr = GF.flTr, kTw = GF.kTw, kkTw = GF.kkTw;
               const double sRo = GF.sRo, rhoR = GF.rhoR, wlR = GF.wlR, twD1 = GF.twD1, twA = GF.twA, twN = GF.twN, twF = GF.twF, limN = GF.limN, limF = GF.limF;
               const bool slipping = GF.slipping, rollS = GF.rollS;
               const M3 Pt = M3::Tangent(nm);
               const V3 rr[2] = { rA, rB };
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
                  M3 Dn = M3::Outer(fd1, nD);   // (patch_friction: the group's end gradient; one group: nm, today's)
                  if (GF.blin)   // level_stiffness (B-lin): kS and ct change with the end position through each piece (kSpV, ctSpV), h dF/dc1 = a1 kSpV^T + a2 ctSpV^T, on Dn's -h
                  {
                     const V3 base = (kS > k0 ? V3() : xi0) + ut * (0.5 * h);
                     V3 a1, a2;
                     if (!slipping)
                     {
                        a1 = base;
                        a2 = ut;
                     }
                     else if (flTr > 0.0)   // past the limit: the trial's change projected off its direction, times lim / |trial|
                     {
                        const V3 fh = FsTr * (1.0 / flTr);
                        const M3 Q = M3::Outer(fh, fh);
                        const double r = limS / flTr;
                        M3 Pf;
                        for (int i = 0; i < 3; i++)
                           for (int j = 0; j < 3; j++)
                              Pf.a[i][j] = Pt.a[i][j] - Q.a[i][j];
                        a1 = MulV(Pf, base) * r;
                        a2 = MulV(Pf, ut) * r;
                     }
                     const M3 D1 = M3::Outer(a1, GF.kSpV), D2 = M3::Outer(a2, GF.ctSpV);
                     for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                           Dn.a[i][j] -= D1.a[i][j] + D2.a[i][j];
                  }
                  if (set.stickZoneTrial && !multiF && slipping && limS > 0.0 && flTr > 0.0)   // (patch_friction: one item's record: one group only)   // stick_zone_trial: the trial force's linear model along a move (kS, ct, stretch, geometry held)
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
               if (set.twistSlope && (twD1 != 0.0 || twA != 0.0 || twN != 0.0 || twF != 0.0 || GF.blin))
                  for (int j = 0; j < na + nq; j++)
                  {
                     if (ix[j] < 0)
                        continue;
                     double c = h * gm[j] * (twD1 + 0.5 * twA) + twN * (terms ? kv[j] : kn * gm[j]);
                     if (GF.blin && j < 3)   // level_stiffness (B-lin): the sticking twist's change with the ball's end position, its vector slope (a patch_friction contact: one ball, no part)
                        c += h * (j == 0 ? GF.twD1V.x : j == 1 ? GF.twD1V.y : GF.twD1V.z);
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
               // a massless column: the slope of the work its contact's mean reaction does on it, qf q (MasslessEnergy), over every unknown: q dqf/du + qf on its own.
               // dqf/du: the normal push -gd N by the slope a coordinate with mass takes in its row (h d(gd N)/du, above), the shear's and the spins' by their change to the row
               for (int j = 0; j < nq; j++)
                  if (mlCol[j] && ix[3 + j] >= 0)
                  {
                     const int qi = ix[3 + j], p3 = 3 + j;
                     double *const wg = (double *)mlWg + (size_t)shape[it.s]->part * jStride;
                     const double *const sn = (const double *)mlSnap + (size_t)j * jStride;
                     const double wq = dq[p3];
                     for (int u = 0; u < nU; u++)
                        wg[u] -= wq * (J[qi][u] - sn[u]) / h;
                     if (fe + fd > 0.0)
                        for (int q = 0; q < na + nq; q++)
                           if (ix[q] >= 0)
                           {
                              double a;
                              if (terms)
                              {
                                 const double kd = di ? dSl : L.gammaKK > 0.0 ? kkJac : cD + cDp * dd;
                                 a = kd * gd[p3] * gd[q];
                                 for (int t = 0; t < hm.nt; t++)
                                    if (tOn[t])
                                       a += h * tk[t] * tgF[t][p3] * tgE[t][q];
                                 if (dirOn)
                                    a += gK[p3][q];
                              }
                              else
                                 a = kn * gm[p3] * gm[q];
                              wg[ix[q]] -= wq * a / h;
                           }
                     wg[qi] += qf[j];
                  }
            }
         }
         if (anyMassless)
            MasslessEnergy(xv, dc, F, J, h);
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
            if (b.held)   // infinite mass (KickerDesign 1h): its rows are its given move, met exactly (residual 0, a row and column of its own); what the contacts' forces on it do over
            {             // that move is the kicker's work on the group, a drive (0 for a parked ball)
               driveW -= F[i].Dot(dc[i]);
               for (int a = 0; a < 6; a++)
               {
                  Res[6 * i + a] = 0.0;
                  if (J)
                  {
                     for (int c = 0; c < nU; c++)
                        J[6 * i + a][c] = J[c][6 * i + a] = 0.0;
                     J[6 * i + a][6 * i + a] = a < 3 ? InertiaSlope(b.m, h) : b.I();
                  }
               }
               continue;
            }
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
         if (J && anyMassless)   // (each massless row's slope: the line search's units, CoordRes2)
            for (int c = 0; c < nCoord; c++)
               if (Massless(c))
                  mlJii[c] = J[cUnk[c]][cUnk[c]];
         for (int a = 0; a < nU; a++)
            if (!std::isfinite(Res[a]))
               return false;
         return true;
      }

      // the first guess of a massless coordinate (design section 1 a: a massless face follows the ball): its own row, the part's force against its contacts' push at the step's end
      // by the full end-push law (EndNormal, as the solve's row: MasslessPush), solved exactly by the part (Part::MasslessGuess) along the depths linearized at the step's
      // start in its own move dq, the balls' guessed moves (x, filled before the parts' BeginStep) and every other unknown held: depth d = a + s dq, a = d0 + g . dc (g = -n grad,
      // the ball's guessed move dc; plus a carried coordinate's slope times its guessed move), s = slope x grad. Every contact the step can reach counts, touching at the start or
      // not (a touch that starts inside the step: its start probe's terms carried on). Each contact's start probe is taken again here; nothing per round, and the part's solve
      // evaluates the law only. One massless column per contact (piece 2's face): with more, each coordinate's guess holds the others' at 0
      struct MlRow : MasslessRow
      {
         const Group *g = nullptr;
         int c = 0;
         double h = 0.0;
         double Push(const double dq, double &dPush) const override
         {
            double P = 0.0;
            dPush = 0.0;
            for (int e = 0; e < g->mlN; e++)
            {
               const MlContact &m = g->mlCt[e];
               if (m.c != c)
                  continue;
               double nS;
               const double N = EndNormal(g->law[m.k], m.d0, m.a + m.s * dq, h, nS);
               P -= m.s * N;
               dPush -= m.s * m.s * nS;
            }
            return P;
         }
      };
      int mlN = 0;   // (mlCt's entries this attempt)
      void MasslessGuess(const double h)
      {
         mlN = 0;
         for (int k = 0; k < ni; k++)
         {
            const Item &it = item[k];
            if (farNow[k] || it.b >= 0 || it.s < 0 || shape[it.s] == nullptr || shape[it.s]->part < 0)
               continue;
            const int p = shape[it.s]->part;
            Hit h0;
            if (nowStamp == stamp && it.nearK < 0)   // (as BeginAttempt takes the start)
               h0 = hNow[k];
            else
               Probe(it, ball, h0);
            double a = d0s[k] - h0.grad * h0.n.Dot(V3(x[6 * it.a], x[6 * it.a + 1], x[6 * it.a + 2]));
            for (int j = 0; j < h0.nc; j++)
            {
               const int c = coordOff[p] + h0.c[j].coord;
               if (!Massless(c) && cUnk[c] >= 0)
                  a += h0.c[j].slope * h0.grad * x[cUnk[c]];
            }
            for (int j = 0; j < h0.nc; j++)
            {
               const int c = coordOff[p] + h0.c[j].coord;
               if (!Massless(c))
                  continue;
               mlCt.Grow((size_t)mlN + 1);
               MlContact &m = mlCt[mlN++];
               m.c = c;
               m.k = k;
               m.a = a;
               m.s = h0.c[j].slope * h0.grad;
               m.d0 = d0s[k];
            }
         }
         for (int p = 0; p < np; p++)
            for (int c = coordOff[p]; c < coordOff[p] + part[p]->NumCoords(); c++)
            {
               if (!Massless(c))
                  continue;
               // where the part's solve starts: the move that keeps the touching contacts' depths (least squares over them; with one, its own: the part's first Newton step is
               // then the row linearized at the start depth). Its tolerance: the solver's on depth (a ball's move within tolV h / 2: Solve's convergence test), as a move of
               // c by its steepest contact
               MlRow row;
               row.g = this;
               row.c = c;
               row.h = h;
               double num = 0.0, den = 0.0, sMax = 0.0;
               for (int e = 0; e < mlN; e++)
                  if (mlCt[e].c == c)
                  {
                     const MlContact &m = mlCt[e];
                     sMax = std::max(sMax, std::fabs(m.s));
                     if (m.d0 > 0.0)
                     {
                        num += m.s * (m.d0 - m.a);
                        den += m.s * m.s;
                     }
                  }
               row.start = den > 0.0 ? num / den : -INFINITY;
               row.tol = sMax > 0.0 ? 0.5 * set.tolV * h / sMax : INFINITY;
               part[p]->MasslessGuess(c - coordOff[p], row, x, partOff[p]);
            }
      }

      // a contact's normal force at the step's end on a massless coordinate (rule C, design section 2 C): N1 = max(0, F(d1) + c(d1) (d1 - max(d0, 0)) / h), the law's elastic
      // force and its damping at the end depth on the step's depth rate (as the paper, runs/2026-11-19-standup-p2/paper_kk.py). nSlope: dN1/dd1 (0 when N1 is 0)
      static double EndNormal(const Item &L, const double d0, const double d1, const double h, double &nSlope)
      {
         double N1 = 0.0;
         nSlope = 0.0;
         if (d1 > 0.0)
         {
            const double f = L.F(d1), c = L.DampC(d1, f), rate = (d1 - std::max(d0, 0.0)) / h;
            N1 = f + c * rate;
            if (N1 > 0.0)
            {
               const double e = 1e-6 * d1;   // (the damping coefficient's slope by a central difference: Newton's only)
               const double cp = (L.DampC(d1 + e, L.F(d1 + e)) - L.DampC(d1 - e, L.F(d1 - e))) / (2.0 * e);
               nSlope = L.K(d1) + cp * rate + c / h;
            }
            else
               N1 = 0.0;
         }
         return N1;
      }
      // a contact's push on a massless coordinate at the step's end: its end normal force (EndNormal) on the end column (-slope1 N1), and the shear's push on its point
      // (Ft . v: the step's, the shear has no end value of its own). nSlope: dN1/dd1 (0 when N1 is 0)
      static double MasslessPush(const Item &L, const Hit &h1, const int j, const double d0, const double d1, const double h, const V3 &Ft, double &nSlope)
      {
         return -h1.c[j].slope * h1.grad * EndNormal(L, d0, d1, h, nSlope) + Ft.Dot(h1.c[j].V());
      }
      // the energy correction of a massless coordinate (design section 3 b): per part, its contacts' stored change, its own (MasslessStoredChange) and the step's losses on those
      // contacts, less the work the contacts' forces do on the balls; by the contacts' own energy books (every contact's forces on the balls and on the part, its mean reaction qf,
      // do minus its stored change and loss) that mismatch is the part's stored change less the work its contacts' mean reaction does on it (mlWork). It is added to the balls
      // touching the part as a Gonzalez correction in the mass metric, the part's weight 0: the ball, the only body with mass, pays the blade's energy. In 1-D the ball's force is
      // then (dU + dV + loss) / dc (paper_kk.py). Its full slope goes in Newton's Jacobian (the mismatch's through mlWg, the metric's through the balls' moves)
      void MasslessEnergy(const double *xv, const V3 *dc, V3 *F, const Jac J, const double h)
      {
         for (int p = 0; p < np; p++)
         {
            bool ml = false;
            for (int c = coordOff[p]; c < coordOff[p] + part[p]->NumCoords(); c++)
               ml = ml || Massless(c);
            if (!ml)
               continue;
            const double dV = part[p]->MasslessStoredChange(xv, partOff[p], nullptr), mis = dV - mlWork[p];
            int n = 0;
            for (int i = 0; i < nb; i++)
            {
               bool touches = false;
               for (int k = 0; k < ni && !touches; k++)
                  touches = ie[k].active && item[k].a == i && item[k].b < 0 && shape[item[k].s] != nullptr && shape[item[k].s]->part == p;
               if (!touches || ball[i].held)
                  continue;
               for (int a = 0; a < 3; a++)
               {
                  mlGm[n + a] = 0.0;
                  mlDq[n + a] = a == 0 ? dc[i].x : a == 1 ? dc[i].y : dc[i].z;
                  mlIm[n + a] = 1.0 / ball[i].m;
               }
               n += 3;
            }
            if (n == 0)
               continue;
            Gonzalez(mlGm, mis, 16.0 * 2.220446049250313e-16 * (std::fabs(dV) + std::fabs(mlWork[p])), mlDq, mlIm, n, set.gonzalezMass, mlOut);
            // its slope (with the Jacobian): out_a = cf W_a dq_a, cf = mis / qq, qq = sum_a W_a dq_a^2 (W the metric's weights, as Gonzalez's), so on ball row a over unknown u
            // h W_a (dq_a dcf/du + cf [u = a]), dcf/du = (dmis/du - cf dqq/du) / qq, dqq/du = 2 W_u dq_u on the balls' own moves; dmis/du = the part's stored change's slope less
            // mlWork's (mlWg). The smooth law's slope inside the rounding band Gonzalez leaves out too
            double cfS = 0.0, qqS = 0.0;
            const auto Wt = [&](const int a) { return !set.gonzalezMass ? mlIm[a] : mlIm[a] > 0.0 ? mlIm[0] / mlIm[a] : 0.0; };
            if (J)
            {
               for (int a = 0; a < n; a++)
                  qqS += mlDq[a] * Wt(a) * mlDq[a];
               if (qqS > 1e-300)
               {
                  cfS = mis / qqS;
                  double *const dm = mlDmis;
                  for (int u = 0; u < nUnk; u++)
                     dm[u] = 0.0;
                  part[p]->MasslessStoredChange(xv, partOff[p], dm);
                  const double *const wg = (const double *)mlWg + (size_t)p * jStride;
                  for (int u = 0; u < nUnk; u++)
                     dm[u] -= wg[u];
                  int s = 0;
                  for (int i = 0; i < nb; i++)
                  {
                     bool touches = false;
                     for (int k = 0; k < ni && !touches; k++)
                        touches = ie[k].active && item[k].a == i && item[k].b < 0 && shape[item[k].s] != nullptr && shape[item[k].s]->part == p;
                     if (!touches || ball[i].held)
                        continue;
                     for (int a = 0; a < 3; a++)
                        dm[6 * i + a] -= cfS * 2.0 * Wt(s + a) * mlDq[s + a];
                     s += 3;
                  }
                  for (int u = 0; u < nUnk; u++)
                     dm[u] /= qqS;   // (dcf/du)
               }
            }
            n = 0;
            for (int i = 0; i < nb; i++)
            {
               bool touches = false;
               for (int k = 0; k < ni && !touches; k++)
                  touches = ie[k].active && item[k].a == i && item[k].b < 0 && shape[item[k].s] != nullptr && shape[item[k].s]->part == p;
               if (!touches || ball[i].held)
                  continue;
               const V3 corr(mlOut[n], mlOut[n + 1], mlOut[n + 2]);
               F[i] = F[i] - corr;
               if (J && qqS > 1e-300)
                  for (int a = 0; a < 3; a++)
                  {
                     const int r = 6 * i + a;
                     const double w = h * Wt(n + a);
                     for (int u = 0; u < nUnk; u++)
                        J[r][u] += w * mlDq[n + a] * mlDmis[u];
                     J[r][r] += w * cfS;
                  }
               if (ImpOn())   // (test tooling: booked on the ball's first contact with the part, the force it corrects)
                  for (int k = 0; k < ni; k++)
                     if (ie[k].active && item[k].a == i && item[k].b < 0 && shape[item[k].s] != nullptr && shape[item[k].s]->part == p)
                     {
                        impIt[4 * (size_t)k] = impIt[4 * (size_t)k] - corr;
                        break;
                     }
               n += 3;
            }
         }
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
         nearAll = true;   // (its explicit moves are not checked against a list: every piece)
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
                  if (!ball[i].held)
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
                     if (!ball[i].held)   // (a held ball's speed is given)
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
                     if (!ball[i].held)
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
                     if (!ball[i].held)
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
         {
            if (ball[i].held)   // (its given move)
            {
               ball[i].c = ball[i].c + ball[i].Given(hTry);
               ball[i].v = ball[i].v + ball[i].heldA * hTry;
            }
            else
               ball[i].c = ball[i].c + ball[i].v * hTry;
            ball[i].adds++;
         }
         for (int p = 0; p < np; p++)
            part[p]->BallsMoved(ball, nb);
         stamp++;
         {
            const double ex = Energy() - E0;
            double K = 0.0;
            for (int i = 0; i < nb; i++)
               if (!ball[i].held)
                  K += 0.5 * ball[i].m * ball[i].v.Dot(ball[i].v) + 0.5 * ball[i].I() * ball[i].w.Dot(ball[i].w);
            for (int p = 0; p < np; p++)
               K += part[p]->KineticJ();
            if (!(ex < K))
            {
               Restore();
               for (int i = 0; i < nb; i++)
                  if (!ball[i].held)
                     ball[i].v = ball[i].v * -1.0;
               for (int p = 0; p < np; p++)
                  part[p]->ScaleSpeeds(-1.0);
               books.energyBounces++;
            }
            else if (ex > 0.0)
            {
               const double sf = std::sqrt(1.0 - ex / K);
               for (int i = 0; i < nb; i++)
                  if (!ball[i].held)
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
         // energy skip test: one trip through the items (runs/2026-11-15-finder-P6b2/paper_stage3.md section 1). A ball is own when all its ball-pair and moving-part items are
         // skipped by the cycle's speed. A ball's static items not touching and not skipped by it are deferred; after the trip, an own ball's are tested on its own bound (OwnFar)
         const bool oneTrip = !farOff && skNb == nb;
         if (oneTrip)
         {
            skOwn.assign((size_t)nb, 1);
            skEs.assign((size_t)nb, -1.0);
            skDef.clear();
            skDefH.clear();
            skipTripAttempts++;
            skipTripItems += ni;
         }
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
            {
               h0s[k].nt = 0;
               Hit &hs = h0s[k];   // patch_friction: the start's patches and their pieces (the friction groups' start side and memory; without them every group's d0 is a plane's)
               hs.ne = h0.ne;   // (step pieces: the entries, as every copy of a Hit)
               hs.en.Fit(hs.ne);
               for (int j = 0; j < hs.ne; j++)
                  hs.en[j] = h0.en[j];
               hs.npat = h0.npat;
               hs.pat.Fit(hs.npat);
               for (int j = 0; j < hs.npat; j++)
                  hs.pat[j] = h0.pat[j];
               hs.npc = h0.npc;
               hs.pc.Fit(hs.npc);
               hs.pd.Fit(hs.npc);
               for (int j = 0; j < hs.npc; j++)
               {
                  hs.pc[j] = h0.pc[j];
                  hs.pd[j] = h0.pd[j];   // (level_stiffness: each piece's data)
               }
            }
            farNow[k] = !farOff && OutOfReach(item[k], h0.delta, h, h0);
            if (farNow[k])
               continue;
            if (oneTrip)
            {
               if (item[k].b < 0 && item[k].s >= 0 && shape[item[k].s] != nullptr && shape[item[k].s]->part < 0)
               {
                  if (!item[k].touching)
                  {
                     skDef.push_back(k);
                     skDefH.push_back(h0);
                     continue;
                  }
               }
               else   // (a ball pair or a moving part not skipped: it can feed both its balls in the step)
               {
                  skOwn[(size_t)item[k].a] = 0;
                  if (item[k].b >= 0)
                     skOwn[(size_t)item[k].b] = 0;
               }
            }
            law[k] = item[k];
            if (!item[k].touching)
               TakeLaw(law[k], h0);
         }
         if (oneTrip)   // (the law is per item: taking it after the trip changes no sum)
         {
            skipDeferred += (long long)skDef.size();
            for (size_t j = 0; j < skDef.size(); j++)
            {
               const int k = skDef[j];
               if (skOwn[(size_t)item[k].a] && OwnFar(item[k], d0s[k], h))
               {
                  farNow[k] = true;
                  continue;
               }
               law[k] = item[k];
               TakeLaw(law[k], skDefH[j]);
            }
         }
         for (int i = 0; nf > 0 && i < nb; i++)
            fieldE0[i] = ball[i].fieldScale != 0.0 ? FieldEnergy(ball[i]) : 0.0;
         for (int i = 0; i < nb; i++)
         {
            // the first guess: the last step's accelerations carried on, or in the stop search this step's own from its last solve (a guess only: the solve's
            // answer does not depend on it)
            const V3 &ag = warmOn ? warmA[i] : aLast[i], &alg = warmOn ? warmAl[i] : alphaLast[i];
            const V3 g0 = ball[i].held ? ball[i].Given(h) : (ball[i].v + ag * (0.5 * h)) * h, w0 = ball[i].held ? V3() : alg * h;
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
         if (anyMassless)
            MasslessGuess(h);
         for (int it = 0; !InRange(x); it++)   // a first guess past a layer's range is shortened
         {
            if (it == 60)
            {
               failWhy = 'G';
               return false;
            }
            for (int a = 0; a < nU; a++)
               x[a] *= 0.5;
            PinHeld(h);
         }
         return true;
      }
      // a held ball's unknowns: its given move, spin 0 (its rows keep them: residual 0, a row and column of their own)
      void PinHeld(const double h)
      {
         for (int i = 0; i < nb; i++)
            if (ball[i].held)
            {
               const V3 g = ball[i].Given(h);
               x[6 * i] = g.x; x[6 * i + 1] = g.y; x[6 * i + 2] = g.z;
               x[6 * i + 3] = x[6 * i + 4] = x[6 * i + 5] = 0.0;
            }
      }

      // Solve one step of h from the state now (the end values stay in x and ie[])
      // shared_finder re-engage (build 3): a surface the ball's inside check found joins as a skipped contact (far) when the step's path can reach its box (the owner's test), with no re-run
      // of its own: the check after the solve (FarStillFar) runs the attempt again only when a skipped contact is reached. Finds made by the checks' own probes are joined the same way and
      // tested by FarStillFar again, so nothing found is left at the step's end.
      bool Solve(double h)
      {
         nearAll = false;
         // VPX_CHECK_PROFILE (timeJoins, set by the owner; diagnostics only): an attempt run again in a step with a join (FarStillFar reached a contact, or a join after an all-in attempt)
         const auto t = [&]() { return timeJoins ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point(); };
         const auto timedRerun = [&](const std::chrono::steady_clock::time_point t0) {
            if (timeJoins && joinRounds > 0)
            {
               joinRerunNs += std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count();
               joinReruns++;
            }
         };
         for (;;)
         {
            farOff = false;
            const bool ok0 = SolveOnce(h);
            Joined();   // (joined as skipped contacts; a failed attempt halves the step with them in)
            if (!ok0)
               return false;
            bool again = false;
            for (;;)
            {
               if (!FarStillFar(x))
               {
                  const auto t0 = t();
                  farOff = true;
                  const bool ok = SolveOnce(h);
                  farOff = false;
                  timedRerun(t0);
                  if (Joined())
                  {
                     again = true;   // found in the all-in attempt: the attempt again, its new contacts classified at its start
                     break;
                  }
                  return ok;   // (its lists hold every piece)
               }
               if (!NearStillNear(x, h))   // a ball centre outside its query volume: its re-query (CycleStillInside), the attempt again
               {
                  again = true;
                  break;
               }
               if (!Joined())
                  return true;
               // found by the checks' own probes: joined as skipped contacts, tested by FarStillFar again
            }
            if (again)
               continue;
         }
      }
      // shared_finder re-engage: the owner's join after an attempt (midJoin); true = a shape was added, run the attempt again. Bounded per step (kJoinRounds): past it, a fault logged once per table use
      bool Joined()
      {
         if (!midJoin)
            return false;
         if (joinRounds >= kJoinRounds)
         {
            if (midJoin(false) < 0)
               joinCapHit++;   // (the owner logs it as a fault after the cycle's steps)
            return false;
         }
         if (midJoin(true) <= 0)
            return false;
         joinRounds++;
         return true;
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
         jacNow = set.jacMode != 0 && JacWindowRange() != nullptr ? (JacWindow(phaseMs) ? set.jacMode : 0) : set.jacMode;
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
         trAtt++;   // (VPX_PATCH_TRACE's attempt, round and difference-Jacobian marks)
         for (; round <= set.maxRounds; round++)
         {
            trRound = round;
            if (jacNow != 0)
            {
               trFd = 1;
               NumericJacobian(h);   // into jacFd (the Eval calls it makes are undone by the Eval below)
               trFd = 0;
            }
            partEvalFailed = false;
            if (jacNow == 0 && round > 1 && lastCh < set.tolCheckMove && !gripMoved && Eval(h, x, Res, nullptr, trialCached))
            {
               double rc = 0.0;
               for (int i = 0; i < nb; i++)
               {
                  const double ic = 1.0 / ball[i].m, iw = ball[i].R / ball[i].I();
                  for (int a = 0; a < 3; a++)
                     rc += Res[6 * i + a] * Res[6 * i + a] * ic * ic + Res[6 * i + 3 + a] * Res[6 * i + 3 + a] * iw * iw;
               }
               rc += CoordRes2(Res, h);
               if (rc <= set.tolRes * set.tolRes)
               {
                  converged = true;
                  endValuesHere = true;
                  break;
               }
            }
            jacWinDump = jacNow != 0 && jacWinLines < 2000 && JacWindow(phaseMs);
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
            if (jacNow != 0)
            {
               for (int r = 0; r < nU; r++)
                  for (int c = 0; c < nU; c++)
                  {
                     JacAna()[r][c] = J[r][c];
                     if (jacNow == 2)
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
            rn += CoordRes2(Res, h);
            rn = std::sqrt(rn);
            if (jacNow != 0 && jacWinLines < 2000 && JacWindow(phaseMs))   // VPX_GROUP_JAC_WINDOW (diagnostics, read only)
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
            rd.mv = std::max(ch, gch);
            rd.rn = rn;
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
               if (jacNow != 0)   // VPX_GROUP_JAC_MODE: each part coordinate's row and column, the solve's own Jacobian (a) vs central differences (f), against itself and each ball it couples to (move x y z, spin x y z)
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
         PzNetTell();   // patch_zero (F6): the net of the Eval this step is accepted on (every Accept: Step's, the energy step's)
         nearOn = false;   // (from here the balls are at the step's end: the probes take every piece)
         const bool hist = np > 0 && HistoryTraceOn();   // (diagnostics, read only)
         const HistE histS0 = hist ? HistEnergy() : HistE();
         const Books histBk0 = hist ? books : Books();
         if (ImpOn())   // (test tooling: the impulses of the Eval this step was accepted on)
            ImpAccept(h);
         PinHeld(h);   // (a held ball's move: given)
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
            b.adds++;
            b.v = EndRate(dc, b.v, h);
            b.w = b.w + dw;
            aLast[i] = evF[i] * b.InvM();            // the next step's first guess: the force and torque of the move taken (Eval's at it), not (v1 - v0) / h, which is
            alphaLast[i] = b.held ? V3() : evM[i] * (1.0 / b.I());   // the solve's tolerance over h for a step near zero long (a stop search on a stop: 5.7e12 m/s^2, a guess 88 km off)
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
               if (k < (int)pmEndOn.size() && pmEndOn[k] > 0)   // patch_friction: the memory on the end pieces (plan 2.1); a piece not pressed now loses it
               {
                  {
                     static std::atomic<long long> nJoin { 0 }, nSplit { 0 }, nDrop { 0 };
                     long long dropped = 0;
                     for (const RubberContact::PieceMem &om : it.pm)
                     {
                        bool kept = false;
                        for (const RubberContact::PieceMem &nm2 : pmEnd[k])
                           if (nm2.id == om.id)
                           {
                              kept = true;
                              break;
                           }
                        dropped += kept ? 0 : 1;
                     }
                     PfCount(nJoin, pmEndJoin[k], "joins in accepted steps", pfMs);
                     PfCount(nSplit, pmEndSplit[k], "splits in accepted steps", pfMs);
                     PfCount(nDrop, dropped, "memory pieces dropped unpressed", pfMs);
                     if (set.linkWindows && k < (int)lwAdd.size())   // link_windows: the start-end patch links its windows added (the accepted Eval's)
                        LwCount(LinkCounters().added, lwAdd[k], "links added beyond shared pieces in accepted steps", pfMs);
                  }
                  it.pm = pmEnd[k];
                  it.pmMulti = pmEndOn[k] == 2;
                  it.pmEs = e.pmEs;
                  it.pmEr = e.pmEr;
                  it.pmLev = e.pmLev;
                  it.kStored = e.kStored;   // (plan 2.4: one group's end stiffness with its memory)
                  it.kt1 = e.kt1;
                  it.kReff = e.kReff;
                  it.kr1 = e.kr1;
                  it.kw1 = e.kw1;
                  if (it.pmMulti)   // (VPX_GROUP_STATS-free count: accepted steps with 2 or more friction groups, running totals at 1, 2, 4, ...; the bench never closes the table)
                  {
                     static std::atomic<long long> nMulti { 0 };
                     const long long c = ++nMulti;
                     if ((c & (c - 1)) == 0)
                        PLOGI << "patch_friction: " << c << " accepted contact steps with 2 or more friction groups (this one at table time " << pfMs << " ms, ball " << it.a << ")";
                  }
               }
               else
               {
                  it.pm.clear();
                  it.pmMulti = false;
                  it.pmLev = false;
                  it.kStored = false;
               }
            }
            else
            {
               it.touching = false;
               it.xi = V3();
               it.ra = RubberContact::RollAngle();
               it.pm.clear();
               it.pmMulti = false;
               it.pmLev = false;
               it.kStored = false;
            }
            if (PatchZeroItem(it))   // patch_zero (Z5): each patch's zero follows it out
               BackOutPatches(it);
            else if (!it.tZero.empty())   // a union of terms: each term's zero follows its term out (its raw depth at the step's end), dropped once the term is clear
            {
               Hit hr;
               Probe(it, ball, hr);
               it.dZero = 0.0;
               for (int t = 0; t < hr.nt; t++)
                  for (size_t j = 0; j < it.tZero.size(); j++)
                     if (it.tZero[j].first == hr.t[t].id)
                     {
                        double &z = it.tZero[j].second;
                        z = std::max(0.0, std::min(z, z + hr.t[t].d));
                        break;
                     }
               if (set.jointZeroFloor)   // (a joint's zero follows its own term out only down to its pieces' zeros: it stays off, no energy)
                  it.FloorJointZeros(hr);
               size_t n = 0;
               for (size_t j = 0; j < it.tZero.size(); j++)
                  if (it.tZero[j].second > 0.0)
                  {
                     it.tZero[n++] = it.tZero[j];
                     it.dZero = std::max(it.dZero, it.tZero[j].second);
                  }
               it.tZero.resize(n);
            }
            else if (it.dZero > 0.0)   // backing out, with nothing stored: the zero follows
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

   // [Contact] the elastic seam (multi-contact piece 1): today's elastic force, moved here unchanged (same operations, same order). The compression's discrete gradient
   // over the move (ball a's centre, then ball b's, or the columns' coordinates); a union of hinge terms: each term by its own chord and its own exact gradient (hinge_piece_sum).
   inline Item::ElasticEnd Item::Elastic(const Hit &h0, const Hit &h1, const Hit &hm, const double d0, const double *dq, const double *iM, const int na, const int nq, const double noiseD,
                                         const Settings &set, const bool trace, const bool wantKc, double *gm, double *gd, std::vector<std::array<double, 4>> &tg, std::vector<double> &tk,
                                         std::vector<char> &tOn, std::vector<std::array<double, 4>> &tgE, std::vector<std::array<double, 4>> &tgF) const
   {
      ElasticEnd o;
      const double d1 = h1.delta, dm = hm.delta;
      const V3 &nm = hm.n;
      gm[0] = -nm.x * hm.grad; gm[1] = -nm.y * hm.grad; gm[2] = -nm.z * hm.grad; gm[3] = nm.x; gm[4] = nm.y; gm[5] = nm.z;
      for (int j = 6; j < Group::kGradMax; j++)
         gm[j] = 0.0;
      for (int j = 0; j < nq; j++)
         gm[3 + j] = hm.c[j].slope * hm.grad;
      Group::Gonzalez(gm, d1 - d0, noiseD, dq, iM, na + nq, set.gonzalezMass, gd, trace ? &o.cf0U : nullptr);
      // the elastic force sum_t s_t fe_t gd_t does sum_t s_t dU_t = dU exactly. The contact's normal force N_el and its direction G (the ball's part a unit normal) are that
      // sum's; the damping acts along G on the rate G . dq / h: its loss cD (G . dq)^2 / h is exact for any G, so no misfit reaches the part through it. (A massless part's
      // spring is in its U*: not here.)
      o.terms = hm.nt > 0 && h1.nt == hm.nt && h0.nt == hm.nt && nq == 1;
      double dd = d1 - d0, feT = 0.0;
      if (o.terms)
      {
         if (tk.size() < (size_t)hm.nt)
         {
            tg.resize((size_t)hm.nt);
            tgE.resize((size_t)hm.nt);
            tgF.resize((size_t)hm.nt);
            tk.resize((size_t)hm.nt);
            tOn.resize((size_t)hm.nt);
         }
         double *Fel = o.Fel;
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
            double gdT[Group::kGradMax];
            bool cf0T = false;
            Group::Gonzalez(tg[t].data(), a1 - a0, noiseD, dq, iM, 4, set.gonzalezMass, gdT, trace ? &cf0T : nullptr);
            if (cf0T && t < 64)
               o.tCf0 |= 1ull << t;
            for (int j = 0; j < 4; j++)   // the slope's row side: the direction this term's force acts along (termForceRows), or its mid gradient
               tgF[t][j] = set.termForceRows ? gdT[j] : tg[t][j];
            const double fT = Chord(a0, a1);
            for (int j = 0; j < 4; j++)
               Fel[j] += s * fT * gdT[j];
            tk[t] = s * (std::fabs(a1 - a0) > 1e-13 ? (F(a1) - fT) / (a1 - a0) : 0.5 * K(std::max(std::max(am, a1), 0.0)));
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
      o.fe = o.terms ? feT : Chord(d0, d1);
      o.feT = feT;
      o.dd = dd;
      o.dE = o.fe * dd;
      // the chord's exact slope (F(d1) - chord) / (d1 - d0) (K / 2 for a vanishing change)
      if (wantKc)
      {
         if (std::fabs(d1 - d0) > 1e-13)
            o.kc = (F(d1) - o.fe) / (d1 - d0);
         else
            o.kc = 0.5 * K(std::max(std::max(dm, d1), 0.0));
      }
      return o;
   }

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
         for (int p = 0; p < g.np; p++)
            e0 -= g.part[p]->PendingNetJ();   // an impact taken at CycleBegin (a driven stop's touch) is in Energy() and in this cycle's drive and loss: open the book before it
         t0 = std::chrono::steady_clock::now();
      }
      void AfterStep()
      {
         if (f == nullptr || !perStep)
            return;
         const double now = g.Energy() - e0 + g.books.Losses() - g.books.driveJ, rise = now - last;
         if (rise > 1e-9 && RubberContact::PatchCheckEnv())   // (VPX_PATCH_CHECK: the first 20 step rises > 1e-9 J, with the touching items' memory kind)
         {
            const long long c = ++RubberContact::HeldCheck()[2];
            if (c <= 20)
            {
               std::ostringstream o;
               o << "BOOKRISE ms " << g.pfMs << " step " << g.books.steps << " rise " << rise * 1e6 << " uJ, items";
               for (int k = 0; k < g.ni; k++)
                  if (g.item[k].touching)
                     o << " [" << k << " s " << g.item[k].s << " pmMulti " << g.item[k].pmMulti << " pmLev " << g.item[k].pmLev << " pm " << g.item[k].pm.size() << "]";
               PLOGI << o.str();
            }
         }
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
      double stepStart = 0.0;
      const bool ok = g.StepsOver(dt, done, [&](const double h) { return hFix > 0.0 ? std::min(hFix, dt - done) : h; }, [&]() {
         book.AfterStep();
         g.LogHalved([&](std::string &msg) {   // (the halved line: Group::LogHalved)
            char s[200];   // (a part keeps no name: its kind, and its state after the step)
            snprintf(s, sizeof s, " t=%d ms + %.1f us, step %.3f us, part alone [%s], balls 0, state", tMs, stepStart * 1e6, (done - stepStart) * 1e6, typeid(p).name());
            msg += s + p.HistoryNote();
         });
         stepStart = done;
      });
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
