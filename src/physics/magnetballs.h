// license:GPLv3+
// Iron-core magnet: the forces balls exert on each other through their magnetisation (2026-10-19). Built at load per magnet from its MagnetCore::Model, evaluated
// per physics step. No VPX dependencies (tested offline in runs/2026-10-14-mirror: test_magnetballs.cpp; design and checks: REPORT.md, sections of 2026-10-19).
//
// Model ("engine recipe", optA_recipe.py PLUG_ORDER 8, PAIR_EXT_ORDER 8, TOUCH_ONLY 1):
//   drive    = coil + core seen by each ball (balls = linear spheres, mu 200, potential multipoles): order-5 tables over (offset, gap) built from the bare field
//              and from the field of the core holding ONE ball on the face (the plug), blended by the plug's gap with the on-axis reaction curve; the plug
//              itself (on the axis) gets its zonal order-8 drive from the on-axis solve;
//   balls    = order-4 global solve (Gauss-Seidel) of every ball given, then every pair within the pair-fix range (surface gap < Settings::pairFixGap, 12 mm)
//              re-solved at order 24 inside the pair from its order-8 drive (precomputed pair operators per gap, 0.1 .. 12 mm), the difference to the order-4
//              pair added to both balls, 3 passes (Option A); fields between pair-fixed balls at order 8, all others at order 4;
//   saturation = contact factors from the linear gap field (ball-ball: uniform-pair FEM table; face: magnetcore's ball-saturation rule), cut rule: each
//              contact's factor scales the total linear force on the balls beyond it (tree rooted at the plug);
//   output   = for every ball the CHANGE of its force caused by the other balls; the lone-ball pull stays MagnetCore::Model::Force.
// Not handled (docs/NewPhysicsReference.txt, MAGNETS): plug off the axis (order-4 drive, accuracy unknown), two balls on one face, rings of touching balls
// (extra contacts ignored), balls of another size (tables for the standard ball). Any magnet axis (each magnet in its own frame, 2026-10-20).
#pragma once
#include "physics/magnetcore.h"
#include <complex>
#include <cstdint>
#include <string>
#include <vector>

namespace MagnetBalls
{
   using cd = std::complex<double>;
   struct V3 { double x = 0, y = 0, z = 0; };

   // coil + core drive at a ball centred at offset x (phi = 0), gap h above the face: scaled order-5 coefficients b_bar_nm = b_nm R^n
   struct DriveTable
   {
      std::vector<double> x, h; std::vector<cd> b;    // b[(ix * h.size() + ih) * K5 + k]
      void Build(const MagnetCore::BareTable &T, double ballR);
      bool Eval(double rho, double gap, cd *b5) const; // false outside (no drive)
   };

   struct MagnetDrive
   {
      const MagnetCore::Model *core = nullptr;
      std::vector<DriveTable> bare, plug;              // per level of core
      std::string log;
      bool Build(const MagnetCore::Model &m);
   };

   // origin = centre of the pole face (m); axis = the magnet axis out of the face (unit vector, table frame). Each magnet's drive is looked up in its own frame
   // (offset from the axis, height above the face) and rotated to the table frame before the sources are summed (build step B4).
   struct Source { const MagnetDrive *drive = nullptr; double amps = 0; V3 origin; V3 axis { 0, 0, 1 }; };

   // One part of a ball's dF (build step B7): the solve runs once per physics step, the parts let the caller re-aim it at every substep.
   //   b >= 0: a PAIR part: F follows the line between balls a and b (it was computed with the line pos[a] - pos[b] of this Evaluate);
   //   b < 0 : a MAGNET part: F follows ball a round the axis of source src (a change of a pull from the magnet, here on ball `ball`).
   // The parts of a ball add up to its dF exactly (saturation included: the cut rule spreads a ball's force over the balls toward the plug).
   struct Part { int ball = 0, a = 0, b = -1, src = 0; V3 F; };
   // the re-aim of a part (build step B7): a pair part by the smallest rotation taking the step's line u0 to the current line u1 (unit vectors); a magnet part
   // round the magnet axis n (unit) by the angle dphi its ball moved round that axis
   V3 ReaimPair(const V3 &F, const V3 &u0, const V3 &u1);
   V3 ReaimMagnet(const V3 &F, const V3 &n, double dphi);
   V3 ReaimMagnetCS(const V3 &F, const V3 &n, double c, double s);   // the same with the angle's cosine and sine (no trig in the substeps)
   // pairFix: pairs within the pair-fix range; passes: pair-fix passes taken; change: the force the last sweep or pass could still have moved (N, see stopForce)
   struct Stats { int balls = 0, pairs = 0, pairFix = 0, sweeps = 0, passes = 0, offAxisPlug = 0; bool turned = false; double change = 0; };
   struct Settings
   {
      double pairFixGap = 12e-3;   // m, surface to surface: pairs closer than this get the pair fix (build step B2; the operator table covers 0.1 .. 12 mm)
      // The order-4 sweeps and the pair-fix passes each stop when their last round moved no ball's force by more than half of stopForce (2026-11-02): bound
      // mu0 |a|max |da|max over the balls (2n+1 weighted), the force between two balls' moments when one changes by da. Half: the rounds shrink the error by about
      // half or better, so what is left is at most about the last change again (measured offline: error up to 1.85 x the last change on arrivals).
      // 1 mN: an eighth of the smallest pair force the model counts at all (the gate, MagnetCore::B_GATE: ~8 mN, 1 % of the ball's weight). A ball with no
      // last answer (new, or a cold start) gets today's fixed rounds first (gsCold, passesCold), then the same stop. 0 = no stop (offline references).
      double stopForce = 1e-3;
      // The most rounds are today's counts, so the solve never costs more than before: with a poor start (a ball arriving, a second magnet, the current
      // changing) the stop may not be reached and the error is what these rounds leave, as before (offline: up to 6 sweeps and 6 passes reach 1 mN there).
      int gsRounds = 3;            // most Gauss-Seidel rounds of the order-4 global solve per Evaluate, from the last answer (build step B5; stops on stopForce)
      int passes = 2;              // most feedback passes of the pair fix per Evaluate, from the last answer (build step B5; stops on stopForce)
      int gsCold = 3, passesCold = 2;   // the least, when a ball has no last answer
      bool saturation = true;      // false = linear balls (contact factors 1): offline tests against the linear reference only
      bool linearTotal = false;    // with saturation = false: dF = the TOTAL linear force on each ball (lone part included), for offline tests only
   };

   class Solver
   {
   public:
      Settings set;
      Solver();
      ~Solver();
      // dF (N): per ball, the change of its magnetic force caused by the other balls (positions in m, the same frame as the sources); ids: one per ball
      // (stable across cycles), kept with the result for FieldFromBalls
      void Evaluate(const std::vector<Source> &src, const std::vector<V3> &pos, std::vector<V3> &dF, Stats *st = nullptr, const std::vector<uint64_t> *ids = nullptr,
                    std::vector<Part> *parts = nullptr);
      // field (T) at p from the balls of the LAST Evaluate (their positions and magnetisation then, order 4), leaving out the ball with id `exclude`;
      // false if there is no last result (gate (b) of build step B3)
      bool FieldFromBalls(const V3 &p, uint64_t exclude, V3 &B) const;
      void ClearLast();   // no ball-ball solve this cycle: the last result and the warm start are stale
      double BallRadius() const;
   private:
      struct Impl; Impl *p;
   };
}
