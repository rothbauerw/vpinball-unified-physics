// license:GPLv3+
//
// A band (a rubber stretched round posts) and a sling (a band with a kicker plate, its coil and switches behind it) as a moving part of the contact system
// (docs/UnifiedContactDesign.txt 12.3). One part per rubber: the band's free spans, the plate, its coil driver and the switch blades step together, because the plate's face
// and the band's nodes push on each other.
//
// Laws:
//  - the band: every free span is a tensioned cord of straight legs between lumped nodes (RubberContact::BandChain: Kelvin-Voigt tension from the real leg length, the
//    node's own damping, stop posts as velocity constraints). Two spans of a loop that meet stop closing, with no bounce (SpanConstraint: a velocity constraint
//    between the legs' closest points, momentum kept, so it can only remove energy; the spring between cords, BandChainContact, made energy). The cord bends where it is loaded: a ball that
//    touches a free span has a node within band_follow_mm of its place; idle nodes that relaxed onto the line are merged away (a merge never gains energy).
//  - ball contact: the shared tube law (RubberContact::EvaluateTube: Hertz on the cord + damping + Mindlin + twist + rolling resistance) on the deformed centreline, the
//    surface moving with the leg the ball touches. The reaction is shared between the leg's two ends (a share on an anchor is carried by the post).
//  - Newton's third law at every island substep: the reactions of all the balls act on the nodes over the whole substep (StepInIsland), in which the band takes its own
//    steps: each as long as the stiffness acting on every node at that moment allows (its legs' axial stiffness, the tension's geometric stiffness, the plate's face
//    pressing on it; omega h <= nodeOmegaH) and the damping's explicit limit.
//  - the plate: the model flipper's mechanics with one winding (SlingPlate, hitsling.h); its face presses on the nodes under it (a spring and damper from the rubber's
//    modulus, one tributary width each) and takes the reaction torque in the same band step.
//  - closing speed: the fastest node now, or what the band's stored deformation energy can give its lightest node, plus the plate face's present speed and what its coil
//    and spring can add in a cycle (SurfaceSpeedBound).
//  - the step a ball's contact needs: the reduced mass of the ball and the cord's mass at the contact (the leg's two ends) with the law's own stiffness at the depth the
//    closing speed gives (BallStep, as the bat's: omega h <= ballOmegaH).
// Energy: only the sling coil adds energy. Leg tension stores and returns it; the viscous tension, the nodes' damping, the face's damper and the contact remove it; a node
// stopped by a stop post gives its speed into the post; two spans that meet lose their closing speed.
// Pure SI maths with no VPX types: the same code runs in the offline tests (runs/2026-10-29-band). The table's geometry (the centreline, the spans, the plate's face, the
// switch blades) is measured by the adapter (VPX side) and handed in.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

#include "physics/contactmovingpart.h"
#include "physics/hitsling.h"

namespace BandPart
{
   using V3 = RubberContact::V3;
   using RubberContact::BandChain;

   // A straight run of the band between two consecutive posts: its centreline vertices (indices into the closed loop, last may run past its end)
   struct Span
   {
      int first = 0, last = 0;
      bool free = false;   // long enough to deflect
   };

   // The plate's face, measured by the adapter from the plate's mesh (hitrubber.cpp SlingPlace): at plate angle phi the face at band height is at
   // originN + a / cos(t) + h tan(t) from the span's rest centreline (+ toward the ball), t = phi + rest
   struct PlateGeom
   {
      bool valid = false;
      int span = -1;                 // the span the face looks at
      double nx = 0.0, ny = 0.0;     // unit normal of that span, from the plate toward the band
      double originN = 0.0, a = 0.0, h = 0.0;   // pivot axis place, face offset from the axis, band height above the axis, m
      double rest = 0.0;             // the plate's angle at rest from the touching pose, rad (<= 0: turned back for the gap)
      double f0 = 0.5, f1 = 0.5;     // the face's extent along the span, fractions of the chord
   };

   // A switch blade behind the band: where it senses the band and which way is behind
   struct SwitchSpec
   {
      int span = -1;
      double f = 0.0;                // fraction of the span's chord
      double bx = 0.0, by = 0.0;     // unit vector from the band toward the blade: the band pushed this way closes the switch
      bool coil = false;             // a sling switch (wired to the coil driver); else read by the script only
   };
   struct SwitchState
   {
      bool closed = false;
      double closedForS = 0.0;
      double dispM = 0.0;            // band displacement toward the blade at the switch, m
   };

   struct Settings
   {
      RubberContact::RubberParams rp;     // the band's rubber against a ball (band_shear_zeta applied)
      double hLayer = 0.0;                // layer thickness of the ball contact, m (layer_mm, or the cord's diameter)
      double radius = 0.0;                // cord radius, m
      RubberContact::BandParams bp;       // the cord (stops in SI, their radius including the cord's)
      bool follow = true;                 // band_follow
      double followM = 6e-3;              // band_follow_mm
      double nodeSpacingM = 0.0;          // band_node_mm: permanent nodes this far apart on every free span (a sling's too); 0 = nodes where the balls load it
      bool ringing = true;                // band_ringing
      double ringZeta = -1.0;             // band_ring_zeta (< 0 = the nodes' own)
      double restM = 5e-5;                // band_rest_mm
      double releaseGapM = 0.2e-3;        // a ball this far from the band counts as still touching (the BandHit event's hysteresis)
      SlingSettings sling;                // the plate, its coil and its driver
      double switchGapM = 1.5875e-3;      // a script switch blade closes at this displacement
      double switchMinClosedS = 2e-3;
      int maxNodes = 12;                  // per span (fixed nodes: FillNodes' own limit)
      double maxStretch = 0.0;            // the largest leg stretch the linear tension law holds for (the model's valid range); 0 = not set
      double nodeOmegaH = 1.0;            // the band's own step: omega h <= this for every node (found offline, runs/2026-10-29-band)
      double ballOmegaH = 0.05;           // a ball's contact: omega h <= this (the island's first-order coupling, as the bat)
      bool onePoint = true;               // a ball touches the band at its nearest point only (one contact); false = every leg it presses into (Patches; made energy
                                          // where the ball pressed 2 - 3 legs: runs/2026-10-29-band, 30 deg hit +61 mJ with no losses)
      double plateMaxStepS = 20e-6;       // the plate's own step (its link and stops), at most
   };

   // A ball touching the band: which nodes it loads (merge), the BandHit event and the script's hit values, and the shear state of each of its contact patches (a ball
   // pressed into a cord with nodes a few mm apart touches several legs at once: more than the three states SurfaceContact holds)
   struct Touch
   {
      static constexpr int kSlots = 8;
      RubberContact::RubberState slotState[kSlots];
      int ball = -1;
      int chain = -1;
      double f = 0.5;           // its place on the chain, fraction of the span
      double peak = 0.0;        // the largest deflection of the cord at the contact, m
      double speedIn = 0.0;     // the ball's speed when it touched, m/s
      int slotKey[kSlots] = {};      // the patch a shear state belongs to: the span (>= 0) of a leg, or -1 - r for rigid run r
      double slotF[kSlots] = {};     // a leg's patch: where it is along the span (fraction)
      bool slotUsed[kSlots] = {};
   };

   // The energy books and counters (diagnostics; not saved with the state)
   struct Books
   {
      double faceOnNodesJ = 0.0;   // work of the plate's face on the nodes
      double faceOnPlateJ = 0.0;   // work of the nodes' reaction on the plate (its torque times its turn)
      double ballOnNodesJ = 0.0;   // work of the balls' reactions on the nodes
      double maxStretch = 0.0;     // the largest leg stretch seen (prestretch included)
      long bandSteps = 0, islandCalls = 0;
      // detail (off unless the offline test sets it: it costs an energy sum per band step):
      //   chainContactJ  the nodes' kinetic energy change made by the spans meeting (SpanConstraint): never above 0
      //   residualJ      per band step, the band's energy change less the work of every force applied from outside the chain step (balls, face, the spans' contact):
      //                  what the chain step itself made or lost (its damping, the stop posts, the integrator)
      //   residualStopJ  the part of residualJ from band steps in which a node was within reach of a stop post (the stops' velocity rule acts there)
      //   addJ, mergeJ   the band's energy change when nodes are added (under a ball, at the plate's face, the band made deformable) or merged away (re-shared masses, a
      //                  merged node's own motion, the two legs becoming one): no force does this work
      //   resetJ         the deformation energy dropped when the band is put back at rest
      bool detail = false;
      double chainContactJ = 0.0, residualJ = 0.0, residualStopJ = 0.0, addJ = 0.0, mergeJ = 0.0, resetJ = 0.0;
      int lastPatches = 0;         // the patches of the last ContactItem call (the trace)
   };

   class Part : public ContactMovingPart
   {
   public:
      // ---- shape and settings: set by the adapter, then Configure ----
      Settings set;
      std::vector<V3> pts;                  // the rest centreline, closed: n + 1 points, the last = the first (m)
      std::vector<Span> spans;
      bool hasPlate = false;
      PlateGeom pg;
      std::vector<SwitchSpec> switchSpec;
      // ---- inputs from the table between cycles (not state) ----
      bool slingEnabled = true;
      bool triggerReq = false;              // the script fired one pulse (taken at the next CycleBegin)
      // ---- dynamic state ----
      bool active = false;                  // the band is deformable now (a ball loads it, it rings, the plate presses on it)
      std::vector<BandChain> chains;        // one per free span
      std::vector<int> chainSpan;           // per chain: its span
      SlingPlate plate;
      SlingDriver driver;
      bool plateIdle = false;
      SlingPlate::State plateIdleState;
      bool anyCoilClosed = false;
      std::vector<SwitchState> sw;
      std::vector<Touch> touches;
      Touch lastHit;                        // the last ball that left the band (the script's hit values)
      double lastHitSpeedOut = 0.0;
      // ---- derived from the state (Rebuild) ----
      std::vector<V3> dpts;                 // the deformed centreline: the ring with the span of every chain that has nodes replaced by A - P_1 - ... - P_M - B
      std::vector<int> base;                // per chain with nodes: index of its A in dpts
      std::vector<int> segChain, segLeg;    // per segment of dpts: its chain (-1 = a post wrap or a short span) and leg (-1 = the rest line of a chain with no nodes)
      std::vector<double> restE;            // per chain: its elastic energy at rest (prestretch, pinned corners)
      Books books;
      RubberContact::BallSI ball;

      bool FixedNodes() const { return set.nodeSpacingM > 0.0; }

      void Configure()
      {
         n = (int)pts.size() - 1;
         RubberContact::BuildTubePre(restPre, pts.data(), (int)pts.size());
         restSegSpan.assign(std::max(n, 0), -1);
         for (size_t k = 0; k < spans.size(); k++)
            if (spans[k].free)
               for (int m = spans[k].first; m < spans[k].last; m++)
                  restSegSpan[((m % n) + n) % n] = (int)k;
         origin = V3();
         for (int i = 0; i < n; i++)
            origin = origin + pts[i] * (1.0 / n);
         double extent = 0.0, longest = 0.0;
         for (int i = 0; i < n; i++)
            extent = std::max(extent, (pts[i] - origin).Length());
         for (const Span &sp : spans)
            if (sp.free)
               longest = std::max(longest, (pts[sp.last % n] - pts[sp.first % n]).Length());
         // how far the cord can bow out: a node in the middle with both legs at the valid range's stretch (half the chord when the range is not set)
         double bow = 0.5 * longest;
         if (set.maxStretch > 0.0)
         {
            const double half = 0.5 * longest * (1.0 + set.maxStretch) / (1.0 + set.bp.prestretch);
            bow = std::sqrt(std::max(0.0, half * half - 0.25 * longest * longest));
         }
         reach = extent + set.radius + bow;
         // the rigid runs: the centreline between one free span's B and the next free span's A (the post wraps and the short spans), each probed as a ring rubber
         runs.clear();
         std::vector<int> freeSpans;
         for (size_t k = 0; k < spans.size(); k++)
            if (spans[k].free && spans[k].last - spans[k].first >= 1 && spans[k].last - spans[k].first < n)
               freeSpans.push_back((int)k);
         if (freeSpans.empty())
         {
            Run r;
            r.pts = pts;
            runs.push_back(r);
         }
         for (size_t q = 0; q < freeSpans.size(); q++)
         {
            const Span &a = spans[freeSpans[q]], &b = spans[freeSpans[(q + 1) % freeSpans.size()]];
            Run r;
            r.spanBefore = freeSpans[q];
            r.spanAfter = freeSpans[(q + 1) % freeSpans.size()];
            int end = b.first;
            while (end < a.last)
               end += n;
            for (int m = a.last; m <= end; m++)
               r.pts.push_back(pts[((m % n) + n) % n]);
            runs.push_back(r);
         }
         for (Run &r : runs)
            RubberContact::BuildTubePre(r.pre, r.pts.data(), (int)r.pts.size());
         sw.assign(switchSpec.size(), SwitchState());
         touches.clear();
         Deactivate();
         plateAlphaMax = faceLeverMax = 0.0;
         if (hasPlate)
         {
            plate.Configure(set.sling, true);
            driver.SetPulse(set.sling.pulseMs * 1e-3);
            plateIdle = false;
            // the largest angular acceleration the coil and the spring can give the plate, and the face's largest lever over the stroke (SurfaceSpeedBound)
            const FlipperSettings &s = plate.set.mech;
            const FlipperDerived &p = plate.der;
            double leverMax = 0.0, springMax = 0.0;
            for (int i = 0; i <= 100; i++)
            {
               const double phi = s.stopDeg * FlipperModel::DEG * i / 100.0;
               const FlipperModel::Kin k = FlipperModel::Kinematics(s, p, p.theta0 + phi);
               leverMax = std::max(leverMax, std::fabs(k.lever) / 1000.0);
               springMax = std::max(springMax, std::fabs(FlipperModel::SpringTorqueK(s, p, k)));
               if (pg.valid)
                  faceLeverMax = std::max(faceLeverMax, std::fabs(FaceDN(phi)));
            }
            plateAlphaMax = (p.strengthP * leverMax + springMax) / std::max(p.iFlipper, 1e-9);
         }
      }

      // ---- the moving-part interface (contactmovingpart.h) ----
      void SetBall(const RubberContact::BallSI &b) override { ball = b; }
      void SetGravity(const V3 &) override {}   // the cord is held at its height by the posts (band_lock_z) and the plate turns in its own plane: gravity does not move either
      V3 Origin() const override { return origin; }
      double Reach() const override { return reach; }

      RubberContact::TubeProbe Probe(const V3 &c) const override
      {
         if (active)
            return RubberContact::ProbeTube(ball, dpts.data(), (int)dpts.size(), set.radius, c);
         return RubberContact::ProbeTube(ball, pts.data(), (int)pts.size(), set.radius, c, &restPre);
      }

      void SurfaceVelocity(const V3 &point, V3 &v, V3 &w) const override
      {
         w = V3();
         v = V3();
         if (!active)
            return;
         int bestK = -1;
         double bestD = 1e30, bestS = 0.0;
         for (size_t k = 0; k + 1 < dpts.size(); k++)
         {
            if (segChain[k] < 0 || segLeg[k] < 0)
               continue;
            const V3 ab = dpts[k + 1] - dpts[k];
            const double l2 = ab.Dot(ab);
            const double s = l2 > 1e-18 ? std::min(std::max((point - dpts[k]).Dot(ab) / l2, 0.0), 1.0) : 0.0;
            const double d = (point - (dpts[k] + ab * s)).Length();
            if (d < bestD)
            {
               bestD = d;
               bestK = (int)k;
               bestS = s;
            }
         }
         if (bestK < 0)
            return;
         const BandChain &ch = chains[segChain[bestK]];
         const int j = segLeg[bestK];
         v = ch.PtVel(j) * (1.0 - bestS) + ch.PtVel(j + 1) * bestS;
      }

      // One substep of one ball's contact. The cord bends where it is loaded: a ball that touches a free span, or can reach it within this substep, gets a node within
      // band_follow_mm of its place first (on the current shape, with the leg's velocity there: no energy is added).
      // The ball touches every part of the band it presses into, each a contact of its own with the shared tube law and its own shear state (Patches): every leg is a straight
      // cord (its stiffness depends on the depth alone, so the elastic force is the gradient of an energy: it cannot make energy), the post wraps and short spans are a ring
      // rubber as before. A ball in the bend under a node presses on both legs at once. Each leg's reaction is kept on its two ends for StepInIsland (a share on an anchor
      // is carried by the post). Measured 2026-10-29 (runs/2026-10-29-band, short3 elastic): the nearest point alone, with the bend at the leg's ends setting its stiffness,
      // gave back 46.8 mJ more than it took.
      RubberContact::ContactResult ContactItem(int ballIndex, SurfaceContact &it, const V3 &c, const V3 &v, const V3 &w, double h) override
      {
         const RubberContact::TubeProbe near = Probe(c);
         if (set.follow || !it.partTouching)
         {
            int hc, hj;
            double s;
            LegOf(near, hc, hj, s);
            const V3 vs0 = (hc >= 0 && hj >= 0) ? chains[hc].PtVel(hj) * (1.0 - s) + chains[hc].PtVel(hj + 1) * s : V3();
            const double gap = near.dist - (ball.R + set.radius);
            const double closing = std::max(0.0, -(v - vs0).Dot(near.normal));
            if (near.overlap || gap <= closing * h)
               NodeUnder(near);
         }
         Patch ps[kMaxPatches];
         const int np = set.onePoint ? NearestPatch(c, ps[0]) : Patches(c, ps);
         books.lastPatches = np;
         Touch *t = FindTouch(ballIndex);
         if (np > 0 && t == nullptr)
         {
            touches.push_back(Touch());
            t = &touches.back();
            t->ball = ballIndex;
            t->speedIn = v.Length();
         }
         RubberContact::RubberState *states[Touch::kSlots];
         for (int q = 0; q < Touch::kSlots; q++)
            states[q] = t != nullptr ? &t->slotState[q] : nullptr;
         RubberContact::ContactResult sum;
         sum.dist = near.dist;
         bool matched[Touch::kSlots] = {};
         double deepest = -1.0;
         for (int k = 0; k < np; k++)
         {
            Patch &p = ps[k];
            // its shear state: the slot that held the same patch last time (the same rigid run, or the same span within band_follow_mm), else a fresh one
            int slot = -1;
            for (int q = 0; q < Touch::kSlots && slot < 0; q++)
               if (t->slotUsed[q] && !matched[q] && t->slotKey[q] == p.key && (p.key < 0 || std::fabs(t->slotF[q] - p.f) * p.Ls <= set.followM))
                  slot = q;
            for (int q = 0; q < Touch::kSlots && slot < 0; q++)
               if (!t->slotUsed[q] && !matched[q])
               {
                  slot = q;
                  states[q]->Reset();
               }
            for (int q = 0; q < Touch::kSlots && slot < 0; q++)
               if (!matched[q])
               {
                  slot = q;
                  states[q]->Reset();
               }
            matched[slot] = true;
            t->slotUsed[slot] = true;
            t->slotKey[slot] = p.key;
            t->slotF[slot] = p.f;
            const RubberContact::ContactResult r = RubberContact::EvaluateTube(set.rp, ball, p.pr, set.hLayer, v, w, h, *states[slot], p.vs);
            if (!r.inContact)
               continue;
            if (p.chain >= 0)
            {
               const V3 reaction = r.force * -1.0;
               const BandChain &ch = chains[p.chain];
               const int n0 = ch.PtNode(p.leg), n1 = ch.PtNode(p.leg + 1);
               if (n0 >= 0)
                  fBall[p.chain][n0] = fBall[p.chain][n0] + reaction * (1.0 - p.s);
               if (n1 >= 0)
                  fBall[p.chain][n1] = fBall[p.chain][n1] + reaction * p.s;
            }
            sum.inContact = true;
            sum.force = sum.force + r.force;
            sum.torqueOnBall = sum.torqueOnBall + r.torqueOnBall;
            sum.fn += r.fn;
            if (r.delta > deepest)
            {
               // the summary of the contact for the engine: the deepest patch's place and depth (the force and torque are the sums)
               deepest = r.delta;
               sum.normal = r.normal;
               sum.point = r.point;
               sum.delta = r.delta;
               sum.patchA = r.patchA;
               sum.reff = r.reff;
               if (p.chain >= 0)
               {
                  t->chain = p.chain;
                  t->f = p.f;
               }
            }
         }
         // a shear state whose patch has ended starts again from rest
         if (t != nullptr)
            for (int q = 0; q < Touch::kSlots; q++)
               if (!matched[q])
               {
                  t->slotUsed[q] = false;
                  states[q]->Reset();
               }
         if (sum.inContact)
         {
            it.partTouching = true;
            if (t->chain >= 0 && t->chain < (int)chains.size())
               t->peak = std::max(t->peak, DeflectionAt(t->chain, t->f));
         }
         else if (it.partTouching && near.dist - (ball.R + set.radius) > set.releaseGapM)
         {
            it.partTouching = false;
            EndTouch(ballIndex, v.Length());
         }
         return sum;
      }

      // The longest substep this ball's contact allows: the law's own stiffness (Hertz with the thin layer's stiffening, as the bat's BallStep) at the depth the relative
      // closing speed reaches, on the reduced mass of the ball and the cord's mass at the contact (the leg's two ends; an anchor is the post).
      double BallStep(const V3 &c, const V3 &v, double hMax) const override
      {
         const RubberContact::TubeProbe pr = Probe(c);
         if (!pr.overlap && pr.dist - (ball.R + set.radius) > 0.002)
            return hMax;   // not near: the arrival rule of the contact system bounds the step
         int hc, hj;
         double s;
         LegOf(pr, hc, hj, s);
         double mCord = 0.0;
         V3 vs;
         if (hc >= 0 && hj >= 0)
         {
            const BandChain &ch = chains[hc];
            const int n0 = ch.PtNode(hj), n1 = ch.PtNode(hj + 1);
            const double inv = (n0 >= 0 ? (1.0 - s) * (1.0 - s) / ch.nodes[n0].m : 0.0) + (n1 >= 0 ? s * s / ch.nodes[n1].m : 0.0);
            mCord = inv > 0.0 ? 1.0 / inv : 0.0;
            vs = ch.PtVel(hj) * (1.0 - s) + ch.PtVel(hj + 1) * s;
         }
         else
         {
            // a free span with no node under the ball yet (one is made at the touch): the mass a single node of that span carries
            const int span = SpanUnder(pr);
            if (span >= 0)
            {
               const double Ls = (pts[spans[span].last % n] - pts[spans[span].first % n]).Length();
               const double area = M_PI * set.radius * set.radius;
               mCord = set.bp.nodeMassKg > 0.0 ? set.bp.nodeMassKg : set.bp.densityKgM3 * area * Ls * 0.5;
            }
         }
         const double mRed = mCord > 0.0 ? ball.m * mCord / (ball.m + mCord) : ball.m;
         const double closing = std::max(std::fabs((v - vs).Dot(pr.normal)), 0.02);
         const RubberContact::RubberParams &rp = set.rp;
         const double Estar = rp.eMpa * 1e6 / (1.0 - rp.nu * rp.nu);
         const double reff = std::max(pr.reff, 1e-4);
         const double dImpact = std::pow(15.0 * mRed * closing * closing / (16.0 * Estar * std::sqrt(reff)), 0.4);
         const double delta = std::max({ pr.overlap ? pr.delta : 0.0, dImpact, 1e-6 });
         const double ah = std::sqrt(reff * delta) / std::max(set.hLayer, 1e-6);
         const double k = 2.0 * Estar * std::sqrt(reff * delta) * (1.0 + 2.0 * rp.betaLayer * ah * ah * ah);
         const double omega = std::sqrt(k / mRed);
         return std::min(hMax, set.ballOmegaH * ContactOmegaHScale() / omega);
      }

      void StepInIsland(double h) override
      {
         books.islandCalls++;
         // the balls' reactions act over the whole substep: their work on the nodes is the reaction times each node's displacement
         pStart.resize(chains.size());
         for (size_t c = 0; c < chains.size(); c++)
         {
            pStart[c].resize(chains[c].nodes.size());
            for (size_t id = 0; id < chains[c].nodes.size(); id++)
               pStart[c][id] = chains[c].nodes[id].P;
         }
         Advance(h, -1.0);
         for (size_t c = 0; c < chains.size() && c < fBall.size(); c++)
            for (size_t id = 0; id < fBall[c].size() && id < pStart[c].size(); id++)
            {
               books.ballOnNodesJ += fBall[c][id].Dot(chains[c].nodes[id].P - pStart[c][id]);
               fBall[c][id] = V3();
            }
         AfterSteps();
         steppedThisCycle = true;
         steppedS += h;
      }

      // No ball near: the band rings down (band_ringing; band_ring_zeta) and the plate moves alone; both sleep when nothing changes
      int AdvanceAlone(double dt) override
      {
         for (auto &f : fBall)
            std::fill(f.begin(), f.end(), V3());
         if (!steppedThisCycle)
            while (!touches.empty())   // no island held the band this cycle: every ball has left it
               EndTouch(touches.back().ball, -1.0);
         if (dt <= 0.0 || Asleep())
            return 0;
         const SlingPlate::State s0 = plate.st;
         const double zRing = set.ringing ? set.ringZeta : 1.0;   // no ringing: critically damped
         const int steps = Advance(dt, zRing);
         AfterSteps();
         if (active && touches.empty())
         {
            bool rest = true;
            for (const BandChain &ch : chains)
               rest = rest && ch.AtRest(set.restM, set.restM * 1000.0);   // m and m/s: 0.05 mm and 0.05 m/s at the default
            if (rest)
            {
               if (books.detail)
                  books.resetJ -= DeformEnergy();
               Deactivate();
            }
         }
         if (hasPlate && !active)
         {
            const SlingPlate::State &s1 = plate.st;
            plateIdle = dt >= 0.9e-3 && !driver.Coil() && driver.remaining <= 0.0 && std::fabs(s1.phi - s0.phi) < 1e-12 && std::fabs(s1.omega - s0.omega) < 1e-9
                        && std::fabs(s1.xp - s0.xp) < 1e-12 && std::fabs(s1.vp - s0.vp) < 1e-9 && std::fabs(s1.iHigh - s0.iHigh) < 1e-12 && s1.rigid == s0.rigid;
            if (plateIdle)
               plateIdleState = s1;
         }
         return steps;
      }

      // ---- once per physics cycle, from the host (not inside the step) ----
      // Before the step: the script's trigger and the switches read at the end of the last cycle drive the coil (a closing starts one pulse)
      void CycleBegin()
      {
         if (!hasPlate)
            return;
         driver.SetPulse(set.sling.pulseMs * 1e-3);
         const int firedBefore = driver.fired;
         bool triggered = false;
         if (slingEnabled)
         {
            if (triggerReq)
               triggered = driver.Trigger();
            driver.Update(anyCoilClosed);
         }
         else
            driver.Cancel();
         triggerReq = false;
         if (driver.fired != firedBefore)
         {
            plateIdle = false;
            events.push_back({ Event::SlingFire, -1, triggered ? 1.0 : 0.0 });
         }
      }
      // After the step: every blade reads the band; a switch closes when the band is pushed back by its gap and opens when it no longer is, not before its minimum closed
      // time. Events: SwitchClosed / SwitchOpened with the blade's index in `ball` and the displacement (m) in `value`.
      void CycleEnd(double dt)
      {
         anyCoilClosed = false;
         for (size_t i = 0; i < switchSpec.size(); i++)
         {
            const SwitchSpec &ss = switchSpec[i];
            SwitchState &st = sw[i];
            const double gap = ss.coil ? set.sling.switchGapMm * 1e-3 : set.switchGapM;
            const double minS = ss.coil ? set.sling.minClosedMs * 1e-3 : set.switchMinClosedS;
            const bool was = st.closed;
            st.dispM = SwitchDisp(ss);
            const bool raw = st.dispM >= gap;
            if (raw && !st.closed)
            {
               st.closed = true;
               st.closedForS = 0.0;
            }
            else if (st.closed)
            {
               st.closedForS += dt;
               if (!raw && st.closedForS >= minS)
                  st.closed = false;
            }
            if (st.closed != was)
               events.push_back({ st.closed ? Event::SwitchClosed : Event::SwitchOpened, (int)i, st.dispM });
            if (ss.coil)
               anyCoilClosed = anyCoilClosed || st.closed;
         }
      }

      // The model's valid range (the safety net, UnifiedContactDesign 12.3 g): a node or the plate is not a number, or a leg is stretched past maxStretch
      bool OutOfRange() const
      {
         if (hasPlate && !std::isfinite(plate.st.phi + plate.st.omega + plate.st.xp + plate.st.vp + plate.st.iHigh))
            return true;
         for (const BandChain &ch : chains)
         {
            for (const int id : ch.order)
            {
               const BandChain::Node &nd = ch.nodes[id];
               if (!std::isfinite(nd.P.x + nd.P.y + nd.P.z + nd.vp.x + nd.vp.y + nd.vp.z))
                  return true;
            }
            if (set.maxStretch > 0.0 && MaxLegStretch(ch) > set.maxStretch)
               return true;
         }
         return false;
      }

      // The band's deformation energy: elastic energy above its rest state plus the nodes' kinetic energy (J)
      double DeformEnergy() const
      {
         double e = 0.0;
         for (size_t c = 0; c < chains.size(); c++)
         {
            e += chains[c].Energy() - restE[c];
            for (const int id : chains[c].order)
               e += 0.5 * chains[c].nodes[id].m * chains[c].nodes[id].vp.Dot(chains[c].nodes[id].vp);
         }
         return e;
      }

      static double MaxLegStretch(const BandChain &ch)
      {
         double e = 0.0;
         for (int j = 0; j <= ch.Count(); j++)
         {
            const double l0 = std::max((ch.PtF(j + 1) - ch.PtF(j)) * ch.L0, 1e-9);
            e = std::max(e, (ch.PtPos(j + 1) - ch.PtPos(j)).Length() / l0 - 1.0);
         }
         return e;
      }

      // The fastest a point of the band can move toward a ball within cycleS: the fastest node, or what the deformation energy stored now can give the lightest node,
      // plus the plate face's speed now and what its coil and spring can add in the cycle
      double SurfaceSpeedBound(double cycleS) const override
      {
         double vNode = 0.0, mMin = 1e30;
         for (const BandChain &ch : chains)
            for (const int id : ch.order)
            {
               vNode = std::max(vNode, ch.nodes[id].vp.Length());
               mMin = std::min(mMin, ch.nodes[id].m);
            }
         if (mMin < 1e29)
            vNode = std::max(vNode, std::sqrt(2.0 * std::max(DeformEnergy(), 0.0) / mMin));
         const double vFace = hasPlate ? (std::fabs(plate.st.omega) + plateAlphaMax * cycleS) * faceLeverMax : 0.0;
         return vNode + vFace;
      }

      // ---- rewind (contactmovingpart.h): everything that changes while the part runs; the shape, the settings and the inputs are not state ----
      void SaveState(State &out) const override
      {
         out.used = 0;
         Put(out, active);
         Put(out, steppedS);
         Put(out, steppedThisCycle);
         Put(out, plate.st);
         Put(out, driver.remaining);
         Put(out, driver.fired);
         Put(out, driver.armed);
         Put(out, plateIdle);
         Put(out, plateIdleState);
         Put(out, anyCoilClosed);
         for (const SwitchState &s : sw)
            Put(out, s);
         Put(out, (int)touches.size());
         for (const Touch &t : touches)
            Put(out, t);
         Put(out, lastHit);
         Put(out, lastHitSpeedOut);
         Put(out, (int)chains.size());
         for (const BandChain &ch : chains)
         {
            Put(out, (int)ch.nodes.size());
            for (const BandChain::Node &nd : ch.nodes)
               Put(out, nd);
            Put(out, (int)ch.order.size());
            for (const int id : ch.order)
               Put(out, id);
         }
      }
      void RestoreState(const State &in) override
      {
         StateReader rd(in);
         rd.Get(active);
         rd.Get(steppedS);
         rd.Get(steppedThisCycle);
         rd.Get(plate.st);
         rd.Get(driver.remaining);
         rd.Get(driver.fired);
         rd.Get(driver.armed);
         rd.Get(plateIdle);
         rd.Get(plateIdleState);
         rd.Get(anyCoilClosed);
         for (SwitchState &s : sw)
            rd.Get(s);
         int nt = 0;
         rd.Get(nt);
         touches.resize(nt);
         for (Touch &t : touches)
            rd.Get(t);
         rd.Get(lastHit);
         rd.Get(lastHitSpeedOut);
         int nc = 0;
         rd.Get(nc);
         if (active)
            InitChains(false);
         else
            Deactivate();
         for (int c = 0; c < nc; c++)
         {
            int nn = 0, no = 0;
            rd.Get(nn);
            BandChain &ch = chains[c];
            ch.nodes.resize(nn);
            for (BandChain::Node &nd : ch.nodes)
               rd.Get(nd);
            rd.Get(no);
            ch.order.resize(no);
            for (int &id : ch.order)
               rd.Get(id);
         }
         if (active)
            Rebuild();
         events.clear();
      }
      // The largest state this part can save: every span at its node limit and maxBalls balls touching (the snapshot reserves it once per part, UnifiedContactDesign 12.3 d)
      size_t MaxStateBytes(int maxBalls) const
      {
         size_t b = 64 + 2 * sizeof(SlingPlate::State) + sw.size() * sizeof(SwitchState) + (size_t)(maxBalls + 1) * sizeof(Touch) + 16;
         int perSpan = set.maxNodes;
         if (FixedNodes())
            perSpan = 60;   // (FillNodes' own limit)
         for (const Span &sp : spans)
            if (sp.free)
               b += 8 + (size_t)perSpan * (sizeof(BandChain::Node) + sizeof(int));
         return b;
      }

      std::string StateNote() const override
      {
         int nodes = 0;
         for (const BandChain &ch : chains)
            nodes += ch.Count();
         return std::string(active ? "deformed, " : "at rest, ") + std::to_string(nodes) + " nodes" + (hasPlate ? ", plate " + std::to_string(plate.PhiDeg()) + " deg, coil " + (driver.Coil() ? "on" : "off") : std::string());
      }

   private:
      int n = 0;                                 // centreline points (the loop without its repeat)
      RubberContact::TubePre restPre;
      std::vector<int> restSegSpan;              // per rest segment: its free span, -1 = a post wrap or a short span
      V3 origin;
      double reach = 0.0;
      double plateAlphaMax = 0.0, faceLeverMax = 0.0;
      std::vector<std::vector<V3>> fBall;        // the balls' reactions on the nodes in the present island substep, per chain by node id
      std::vector<std::vector<V3>> fWork;        // scratch: all the forces on the nodes in one band step
      std::vector<std::vector<V3>> pStart;       // scratch: node positions at the start of an island substep
      std::vector<std::vector<V3>> pPrev;        // scratch (detail books): node positions before a band step
      struct FaceLoad { int id; double F, xn0; };
      std::vector<FaceLoad> faceLoads;           // scratch: the face's push on each node in one band step
      // A rigid run of the centreline (post wraps and short spans between two free spans), probed as a ring rubber
      struct Run
      {
         std::vector<V3> pts;
         RubberContact::TubePre pre;
         int spanBefore = -1, spanAfter = -1;   // the free span whose B is its first point, and the one whose A is its last (-1: a band with no free span)
      };
      std::vector<Run> runs;
      // One patch of a ball's contact with the band: a leg (chain, leg, s along it), or a rigid run (chain -1)
      struct Patch
      {
         RubberContact::TubeProbe pr;
         int chain = -1, leg = -1;
         double s = 0.0;
         int key = 0;          // the span (>= 0) of a leg, or -1 - r for rigid run r
         double f = 0.0;       // a leg's patch: where it is along the span (fraction)
         double Ls = 1.0;      // that span's chord, m
         V3 vs;                // the band's velocity at the patch
      };
      static constexpr int kMaxPatches = Touch::kSlots;

      // A free span's legs, at rest or deformed: points 0 .. LegPoints() - 1 (0 = A, the last = B)
      int LegPoints(int span) const
      {
         const int ci = active ? ChainOfSpan(span) : -1;
         return ci >= 0 ? chains[ci].Count() + 2 : 2;
      }
      V3 LegPt(int span, int k) const
      {
         const int ci = active ? ChainOfSpan(span) : -1;
         if (ci >= 0)
            return chains[ci].PtPos(k);
         return k == 0 ? pts[spans[span].first % n] : pts[spans[span].last % n];
      }
      // Where along leg j of a span the foot of the perpendicular from c falls (0 .. 1)
      double FootS(int span, int j, const V3 &c) const
      {
         const V3 p = LegPt(span, j), e = LegPt(span, j + 1) - p;
         const double l2 = e.Dot(e);
         return l2 > 1e-18 ? std::min(std::max((c - p).Dot(e) / l2, 0.0), 1.0) : 0.0;
      }

      // The one contact of a ball at c with the band (onePoint): the nearest point of the whole band, if the ball presses into it. On a leg the cord is straight (the effective
      // radius of a sphere on a cylinder at every depth: the contact's energy depends on the distance to the band alone, which is continuous, so it cannot make energy); on
      // a rigid run (post wraps, short spans) the ring law with its curvature, as before. Returns 1, or 0 when the ball does not touch.
      int NearestPatch(const V3 &c, Patch &out) const
      {
         const double touch = ball.R + set.radius;
         double best = 1e30;
         int bestSpan = -1, bestLeg = -1, bestRun = -1;
         double bestS = 0.0;
         V3 bestFoot;
         for (size_t k = 0; k < spans.size(); k++)
         {
            const Span &sp = spans[k];
            if (!sp.free || sp.last - sp.first < 1 || sp.last - sp.first >= n)
               continue;
            const int span = (int)k, legs = LegPoints(span) - 1;
            for (int j = 0; j < legs; j++)
            {
               const V3 p = LegPt(span, j), e = LegPt(span, j + 1) - p;
               const double l2 = e.Dot(e);
               const double s = l2 > 1e-18 ? std::min(std::max((c - p).Dot(e) / l2, 0.0), 1.0) : 0.0;
               const V3 foot = p + e * s;
               const double d = (c - foot).Length();
               if (d < best)
               {
                  best = d;
                  bestSpan = span;
                  bestLeg = j;
                  bestS = s;
                  bestFoot = foot;
                  bestRun = -1;
               }
            }
         }
         RubberContact::TubeProbe runPr;
         for (size_t r = 0; r < runs.size(); r++)
         {
            const RubberContact::TubeProbe pr = RubberContact::ProbeTube(ball, runs[r].pts.data(), (int)runs[r].pts.size(), set.radius, c, &runs[r].pre);
            if (pr.seg >= 0 && pr.dist < best)
            {
               best = pr.dist;
               bestRun = (int)r;
               runPr = pr;
            }
         }
         if (best >= touch || best <= 1e-12)
            return 0;
         Patch p;
         if (bestRun >= 0)
         {
            p.pr = runPr;
            p.key = -1 - bestRun;
            out = p;
            return 1;
         }
         const double R1 = 1.0 / (1.0 / ball.R + 1.0 / set.radius);
         p.pr.point = bestFoot;
         p.pr.normal = (c - bestFoot) / best;
         p.pr.dist = best;
         p.pr.delta = touch - best;
         p.pr.overlap = true;
         p.pr.reff = std::sqrt(R1 * ball.R);
         p.pr.rho = set.radius;
         p.chain = active ? ChainOfSpan(bestSpan) : -1;
         p.leg = bestLeg;
         p.s = bestS;
         p.key = bestSpan;
         p.Ls = (LegPt(bestSpan, LegPoints(bestSpan) - 1) - LegPt(bestSpan, 0)).Length();
         if (p.chain >= 0)
         {
            const BandChain &ch = chains[p.chain];
            p.f = ch.PtF(bestLeg) + bestS * (ch.PtF(bestLeg + 1) - ch.PtF(bestLeg));
            p.vs = ch.PtVel(bestLeg) * (1.0 - bestS) + ch.PtVel(bestLeg + 1) * bestS;
         }
         else
            p.f = bestS;
         out = p;
         return 1;
      }

      // Every patch of the band a ball at c presses into (the deepest kMaxPatches of them):
      //  - a leg whose foot lies inside it (a straight cord: the effective radius of a sphere on a cylinder, the same at every depth and every bend)
      //  - a node where the feet of both its legs are (the ball on the outside of the bend): one patch, the node takes all of it
      //  - a rigid run's nearest point (the ring law, curvature included), unless it is the run's end at an anchor while the leg there has its foot inside it (the same
      //    place, counted by the leg)
      int Patches(const V3 &c, Patch *ps) const
      {
         constexpr double kEps = 1e-9;
         int np = 0;
         const double touch = ball.R + set.radius;
         auto add = [&](const Patch &p) {
            if (np < kMaxPatches)
            {
               ps[np++] = p;
               return;
            }
            int w = 0;
            for (int k = 1; k < np; k++)
               if (ps[k].pr.delta < ps[w].pr.delta)
                  w = k;
            if (p.pr.delta > ps[w].pr.delta)
               ps[w] = p;
         };
         const double R1 = 1.0 / (1.0 / ball.R + 1.0 / set.radius);
         const double reffStraight = std::sqrt(R1 * ball.R);
         auto legPatch = [&](int span, int j, double s, const V3 &foot, double d) {
            Patch p;
            p.pr.point = foot;
            p.pr.normal = (c - foot) / d;
            p.pr.dist = d;
            p.pr.delta = touch - d;
            p.pr.overlap = true;
            p.pr.reff = reffStraight;
            p.pr.rho = set.radius;
            p.chain = active ? ChainOfSpan(span) : -1;
            p.leg = j;
            p.s = s;
            p.key = span;
            p.Ls = (LegPt(span, LegPoints(span) - 1) - LegPt(span, 0)).Length();
            if (p.chain >= 0)
            {
               const BandChain &ch = chains[p.chain];
               p.f = ch.PtF(j) + s * (ch.PtF(j + 1) - ch.PtF(j));
               p.vs = ch.PtVel(j) * (1.0 - s) + ch.PtVel(j + 1) * s;
            }
            else
               p.f = s;
            add(p);
         };
         for (size_t k = 0; k < spans.size(); k++)
         {
            const Span &sp = spans[k];
            if (!sp.free || sp.last - sp.first < 1 || sp.last - sp.first >= n)
               continue;
            const int span = (int)k, legs = LegPoints(span) - 1;
            double sPrev = 0.0;
            for (int j = 0; j < legs; j++)
            {
               const V3 p = LegPt(span, j), e = LegPt(span, j + 1) - p;
               const double l2 = e.Dot(e);
               const double s = l2 > 1e-18 ? std::min(std::max((c - p).Dot(e) / l2, 0.0), 1.0) : 0.0;
               const V3 foot = p + e * s;
               const double d = (c - foot).Length();
               if (s > kEps && s < 1.0 - kEps && d < touch && d > 1e-12)
                  legPatch(span, j, s, foot, d);
               if (j > 0 && sPrev >= 1.0 - kEps && s <= kEps && d < touch && d > 1e-12)
                  legPatch(span, j, 0.0, p, d);   // the node between legs j - 1 and j, from outside the bend
               sPrev = s;
            }
         }
         for (size_t r = 0; r < runs.size(); r++)
         {
            const Run &run = runs[r];
            const int m = (int)run.pts.size();
            const RubberContact::TubeProbe pr = RubberContact::ProbeTube(ball, run.pts.data(), m, set.radius, c, &run.pre);
            if (!pr.overlap)
               continue;
            if (run.spanBefore >= 0 && (pr.point - run.pts[0]).Length() < 1e-12 && FootS(run.spanBefore, LegPoints(run.spanBefore) - 2, c) < 1.0 - kEps)
               continue;   // at the anchor B of the span before: that span's last leg has it
            if (run.spanAfter >= 0 && (pr.point - run.pts[m - 1]).Length() < 1e-12 && FootS(run.spanAfter, 0, c) > kEps)
               continue;   // at the anchor A of the span after: its first leg has it
            Patch p;
            p.pr = pr;
            p.key = -1 - (int)r;
            add(p);
         }
         return np;
      }

      Touch *FindTouch(int b)
      {
         for (Touch &t : touches)
            if (t.ball == b)
               return &t;
         return nullptr;
      }
      void EndTouch(int b, double speedOut)
      {
         for (size_t i = 0; i < touches.size(); i++)
            if (touches[i].ball == b)
            {
               lastHit = touches[i];
               lastHitSpeedOut = speedOut;
               events.push_back({ Event::BandHit, b, touches[i].peak });
               touches.erase(touches.begin() + (long)i);
               return;
            }
      }

      bool Asleep() const { return !active && (!hasPlate || (plateIdle && !driver.Coil() && plate.st == plateIdleState)); }

      void Deactivate()
      {
         active = false;
         chains.clear();
         chainSpan.clear();
         base.clear();
         dpts.clear();
         segChain.clear();
         segLeg.clear();
         restE.clear();
         fBall.clear();
      }

      int ChainOfSpan(int span) const
      {
         for (size_t k = 0; k < chainSpan.size(); k++)
            if (chainSpan[k] == span)
               return (int)k;
         return -1;
      }

      // One chain per free span (no nodes; withNodes: the corners pinned on stop posts and the permanent nodes too)
      void InitChains(bool withNodes)
      {
         chains.clear();
         chainSpan.clear();
         restE.clear();
         for (size_t k = 0; k < spans.size(); k++)
         {
            const Span &sp = spans[k];
            if (!sp.free || sp.last - sp.first < 1 || sp.last - sp.first >= n)
               continue;
            BandChain ch;
            ch.Init(pts[sp.first % n], pts[sp.last % n], set.radius, set.bp);
            if (withNodes)
            {
               ch.PinStops();   // a post the band lies against is a fixed corner of the cord
               if (FixedNodes())
                  ch.FillNodes(set.nodeSpacingM);
            }
            restE.push_back(ch.Energy());
            chains.push_back(ch);
            chainSpan.push_back((int)k);
         }
         if (!withNodes)
         {
            // (a restored state brings its own nodes: the rest energy is that of the pinned rest shape, worked out again)
            for (size_t c = 0; c < chains.size(); c++)
            {
               BandChain ch = chains[c];
               ch.PinStops();
               restE[c] = ch.Energy();
            }
         }
         fBall.assign(chains.size(), std::vector<V3>());
      }

      bool EnsureDyn()
      {
         if (active)
            return true;
         InitChains(true);
         if (chains.empty())
            return false;
         active = true;
         Rebuild();
         return true;
      }

      // The deformed centreline and what each of its segments is. Called whenever a node is added or merged away.
      void Rebuild()
      {
         std::vector<char> removed(n, 0);
         std::vector<int> insertAfter(n, -1), restToD(n, -1);
         for (size_t c = 0; c < chains.size(); c++)
            if (chains[c].Count() > 0)
            {
               const Span &sp = spans[chainSpan[c]];
               for (int m = sp.first + 1; m < sp.last; m++)
                  removed[m % n] = 1;
               insertAfter[sp.first % n] = (int)c;
            }
         dpts.clear();
         base.assign(chains.size(), -1);
         for (int i = 0; i < n; i++)
            if (!removed[i])
            {
               restToD[i] = (int)dpts.size();
               dpts.push_back(pts[i]);
               if (insertAfter[i] >= 0)
               {
                  const int c = insertAfter[i];
                  base[c] = (int)dpts.size() - 1;
                  for (int k = 1; k <= chains[c].Count(); k++)
                     dpts.push_back(chains[c].PtPos(k));
               }
            }
         dpts.push_back(dpts[0]);
         segChain.assign(dpts.size() - 1, -1);
         segLeg.assign(dpts.size() - 1, -1);
         for (size_t c = 0; c < chains.size(); c++)
         {
            if (chains[c].Count() > 0)
               for (int j = 0; j <= chains[c].Count(); j++)
               {
                  segChain[base[c] + j] = (int)c;
                  segLeg[base[c] + j] = j;
               }
            else
            {
               const Span &sp = spans[chainSpan[c]];
               for (int m = sp.first; m < sp.last; m++)
                  if (restToD[m % n] >= 0)
                     segChain[restToD[m % n]] = (int)c;
            }
            fBall[c].resize(chains[c].nodes.size());
         }
      }

      void SyncPts()
      {
         for (size_t c = 0; c < chains.size(); c++)
            for (int k = 1; k <= chains[c].Count(); k++)
               dpts[base[c] + k] = chains[c].PtPos(k);
      }

      // The free span under a probe's nearest point (-1 = a post wrap or a short span)
      int SpanUnder(const RubberContact::TubeProbe &pr) const
      {
         if (pr.seg < 0)
            return -1;
         if (!active)
            return pr.seg < (int)restSegSpan.size() ? restSegSpan[pr.seg] : -1;
         const int c = pr.seg < (int)segChain.size() ? segChain[pr.seg] : -1;
         return c >= 0 ? chainSpan[c] : -1;
      }

      // The chain and the leg a probe's nearest point is on, and where along the leg (s: 0 at its start). hc = -1: no chain there; hj = -1: the rest line of a chain with no nodes.
      void LegOf(const RubberContact::TubeProbe &pr, int &hc, int &hj, double &s) const
      {
         hc = hj = -1;
         s = 0.0;
         if (!active || pr.seg < 0 || pr.seg >= (int)segChain.size())
            return;
         hc = segChain[pr.seg];
         hj = segLeg[pr.seg];
         if (hc < 0 || hj < 0)
            return;
         const V3 p = chains[hc].PtPos(hj), q = chains[hc].PtPos(hj + 1);
         const double len = (q - p).Length();
         s = len > 1e-9 ? std::min(1.0, std::max(0.0, (pr.point - p).Length() / len)) : 0.0;
      }

      // The material fraction of the cord at a probe's nearest point (on a leg: between its ends' fractions; on a rest line: along the chord)
      double FractionOf(const RubberContact::TubeProbe &pr, int hc, int hj, double s) const
      {
         const BandChain &ch = chains[hc];
         if (hj >= 0)
            return ch.PtF(hj) + s * (ch.PtF(hj + 1) - ch.PtF(hj));
         const V3 ab = ch.B - ch.A;
         return std::min(std::max((pr.point - ch.A).Dot(ab) / std::max(ab.Dot(ab), 1e-18), 0.0), 1.0);
      }

      // How far the cord at fraction f is from its rest place, m
      double DeflectionAt(int c, double f) const
      {
         const BandChain &ch = chains[c];
         int j = 0;
         while (j < ch.Count() && ch.PtF(j + 1) < f)
            j++;
         const double f0 = ch.PtF(j), f1 = ch.PtF(j + 1);
         const double u = f1 > f0 ? std::min(std::max((f - f0) / (f1 - f0), 0.0), 1.0) : 0.0;
         const V3 P = ch.PtPos(j) + (ch.PtPos(j + 1) - ch.PtPos(j)) * u;
         return (P - (ch.A + (ch.B - ch.A) * f)).Length();
      }

      // A node within band_follow_mm of the probe's nearest point on a free span. True if the shape changed (the band became deformable, or a node was added).
      bool NodeUnder(const RubberContact::TubeProbe &pr)
      {
         const int span = SpanUnder(pr);
         if (span < 0)
            return false;
         const bool wasActive = active;
         double f;
         if (active)
         {
            int hc, hj;
            double s;
            LegOf(pr, hc, hj, s);
            f = FractionOf(pr, hc, hj, s);
         }
         else
         {
            const V3 a = pts[spans[span].first % n], ab = pts[spans[span].last % n] - a;
            f = (pr.point - a).Dot(ab) / std::max(ab.Dot(ab), 1e-18);
         }
         f = std::min(std::max(f, 0.05), 0.95);
         const double e0 = books.detail ? DeformEnergy() : 0.0;
         if (!EnsureDyn())
            return false;
         if (books.detail)
            books.addJ += DeformEnergy() - e0;
         if (FixedNodes())
            return !wasActive;
         const int ci = ChainOfSpan(span);
         if (ci < 0)
            return !wasActive;
         const BandChain &ch = chains[ci];
         for (const int id : ch.order)
            if (std::fabs(ch.nodes[id].f - f) * ch.Ls < set.followM)
               return !wasActive;   // a node is there
         if (ch.Count() >= set.maxNodes)
            return !wasActive;
         const double e1 = books.detail ? DeformEnergy() : 0.0;
         AddNode(ci, f);
         if (books.detail)
            books.addJ += DeformEnergy() - e1;
         return true;
      }

      // A node at fraction f of chain ci, and a passive one on every other chain at the point nearest to it (unless one is within 10 mm), so the loaded span can push on
      // the others (a loop's back span)
      int AddNode(int ci, double f)
      {
         constexpr double kCompanionSpot = 0.010;   // m
         const int nid = chains[ci].Add(f);
         const V3 P0 = chains[ci].nodes[nid].P0;
         for (size_t cj = 0; cj < chains.size(); cj++)
         {
            if ((int)cj == ci)
               continue;
            BandChain &oc = chains[cj];
            const V3 ab = oc.B - oc.A;
            const double l2 = ab.Dot(ab);
            if (l2 < 1e-12 || oc.Count() >= set.maxNodes)
               continue;
            const double fo = std::min(std::max((P0 - oc.A).Dot(ab) / l2, 0.05), 0.95);
            bool close = false;
            for (const int id : oc.order)
               close = close || std::fabs(oc.nodes[id].f - fo) * oc.Ls < kCompanionSpot;
            if (!close)
               oc.Add(fo);
         }
         Rebuild();
         return nid;
      }

      // ---- the plate ----
      double FaceN(double phi) const
      {
         const double t = phi + pg.rest;
         return pg.originN + pg.a / std::cos(t) + pg.h * std::tan(t);
      }
      double FaceDN(double phi) const
      {
         const double t = phi + pg.rest, c = std::cos(t);
         return (pg.a * std::sin(t) + pg.h) / (c * c);
      }
      int PlateChain() const { return hasPlate && pg.valid ? ChainOfSpan(pg.span) : -1; }
      bool InFootprint(int c, int id) const
      {
         const BandChain &ch = chains[c];
         const double tol = set.followM / std::max(ch.Ls, 1e-6);   // a node within a node spacing of the face carries its share
         const double f = ch.nodes[id].f;
         return f >= pg.f0 - tol && f <= pg.f1 + tol;
      }
      // Nodes at the two ends of the face, so the cord between them can lie flat on it (a face narrower than 6 mm gets one node in the middle)
      void EnsurePlateNodes()
      {
         const double e0 = books.detail ? DeformEnergy() : 0.0;
         EnsurePlateNodesImpl();
         if (books.detail)
            books.addJ += DeformEnergy() - e0;
      }
      void EnsurePlateNodesImpl()
      {
         if (!EnsureDyn() || FixedNodes())
            return;   // (fixed nodes: the face presses on the permanent nodes under it)
         const int ci = PlateChain();
         if (ci < 0)
            return;
         BandChain &ch = chains[ci];
         double fs[2] = { pg.f0, pg.f1 };
         int nf = 2;
         if ((pg.f1 - pg.f0) * ch.Ls < 6e-3)
         {
            fs[0] = 0.5 * (pg.f0 + pg.f1);
            nf = 1;
         }
         bool changed = false;
         for (int k = 0; k < nf; k++)
         {
            bool have = false;
            for (const int id : ch.order)
               have = have || std::fabs(ch.nodes[id].f - fs[k]) * ch.Ls < set.followM;
            if (!have && ch.Count() < set.maxNodes)
            {
               ch.Add(fs[k]);
               changed = true;
            }
         }
         if (changed)
            Rebuild();
      }
      // The face against the nodes under it (a spring from the rubber's modulus over one tributary width each, and a damper): the forces go into f, the reaction torque on
      // the plate (in the direction of increasing phi) is returned. The pushes are recorded for the energy books.
      double FaceForces(std::vector<std::vector<V3>> &f, double &kNode)
      {
         faceLoads.clear();
         kNode = 0.0;
         const int ci = PlateChain();
         if (ci < 0 || chains[ci].Count() == 0)
            return 0.0;
         const BandChain &ch = chains[ci];
         const double phi = plate.st.phi;
         const double face = FaceN(phi), dF = FaceDN(phi), vFace = dF * plate.st.omega;
         int count = 0;
         for (const int id : ch.order)
            count += InFootprint(ci, id) ? 1 : 0;
         if (count == 0)
            return 0.0;
         const double width = std::max((pg.f1 - pg.f0) * ch.Ls, 3e-3) / count;
         const double k = set.rp.eMpa * 1e6 * width;
         kNode = k;
         double torque = 0.0;
         for (const int id : ch.order)
         {
            if (!InFootprint(ci, id))
               continue;
            const BandChain::Node &nd = ch.nodes[id];
            const double xn = (nd.P.x - nd.P0.x) * pg.nx + (nd.P.y - nd.P0.y) * pg.ny;   // + toward the ball
            const double pen = face - xn + set.radius;                                   // the face against the back of the cord
            if (pen <= 0.0)
               continue;
            const double vn = nd.vp.x * pg.nx + nd.vp.y * pg.ny;
            const double c = 2.0 * set.sling.plateContactZeta * std::sqrt(k * nd.m);
            const double F = k * pen + c * (vFace - vn);
            if (F <= 0.0)
               continue;
            f[ci][id] = f[ci][id] + V3(pg.nx, pg.ny, 0.0) * F;
            torque -= F * dF;
            faceLoads.push_back({ id, F, xn });
         }
         return torque;
      }

      // ---- the band's own steps ----
      // The longest band step the nodes allow now: omega h <= nodeOmegaH for the stiffness acting on every node (its legs' axial stiffness, the tension's geometric
      // stiffness, the face pressing on it), and the explicit limit of its damping (c h / m <= 1.5; the scheme's limit is 2)
      double NodeStepLimit(double hMax, double zeta, double kFaceNode) const
      {
         constexpr double kDampLimit = 1.5;
         double lim = hMax;
         const int pc = PlateChain();
         double face = 0.0, vFace = 0.0;
         if (pc >= 0)
         {
            face = FaceN(plate.st.phi);
            vFace = FaceDN(plate.st.phi) * plate.st.omega;
         }
         for (size_t c = 0; c < chains.size(); c++)
         {
            const BandChain &ch = chains[c];
            const int M = ch.Count();
            const double z = zeta >= 0.0 ? zeta : ch.zeta;
            for (int i = 0; i < M; i++)
            {
               const BandChain::Node &nd = ch.nodes[ch.order[i]];
               const double lL0 = std::max((ch.PtF(i + 1) - ch.PtF(i)) * ch.L0, 1e-9), lR0 = std::max((ch.PtF(i + 2) - ch.PtF(i + 1)) * ch.L0, 1e-9);
               const double LL = std::max((ch.PtPos(i + 1) - ch.PtPos(i)).Length(), 1e-9), LR = std::max((ch.PtPos(i + 2) - ch.PtPos(i + 1)).Length(), 1e-9);
               const double kAx = ch.EA * (1.0 / lL0 + 1.0 / lR0);
               const double kGeo = ch.LegTension(i) / LL + ch.LegTension(i + 1) / LR;
               double kFace = 0.0;
               if ((int)c == pc && kFaceNode > 0.0 && InFootprint(pc, ch.order[i]))
               {
                  const double xn = (nd.P.x - nd.P0.x) * pg.nx + (nd.P.y - nd.P0.y) * pg.ny;
                  if (face - xn + set.radius + std::fabs(vFace) * hMax > 0.0)   // the face presses on it, or can within the step
                     kFace = kFaceNode;
               }
               const double m = nd.m;
               const double omega = std::sqrt((kAx + kGeo + kFace) / m);
               const double cDamp = ch.eta * (1.0 / lL0 + 1.0 / lR0) + 2.0 * z * std::sqrt(kGeo * m) + (kFace > 0.0 ? 2.0 * set.sling.plateContactZeta * std::sqrt(kFace * m) : 0.0);
               lim = std::min(lim, set.nodeOmegaH / omega);
               if (cDamp > 0.0)
                  lim = std::min(lim, kDampLimit * m / cDamp);
            }
         }
         return lim;
      }

      double NodeKinetic() const
      {
         double e = 0.0;
         for (const BandChain &ch : chains)
            for (const int id : ch.order)
               e += 0.5 * ch.nodes[id].m * ch.nodes[id].vp.Dot(ch.nodes[id].vp);
         return e;
      }

      // Two spans of a loop meeting (the front span pressed onto the back one), with no bounce: for every pair of legs of two spans that touch or would touch within the step h,
      // the part of their closing speed at the closest points beyond what closes the gap within the step is taken away by an equal and opposite impulse, shared between
      // the legs' ends by where the touch is and their masses (an anchor, or a node pinned to a stop, is held: it takes none). Momentum is kept and the closing speed
      // only falls, so each impulse removes energy. Legs already overlapping may not close further but are not pushed apart (the band's tension separates them).
      // Several pairs share nodes: a few passes. Replaces the spring between the cords (BandChainContact, Hertz), which made energy: measured 2026-10-29
      // (runs/2026-10-29-band, no losses: hit3_a0 +3,233 mJ, +69 mJ with its stiffness in the step rule; hit1_a0 +0.43 mJ at 1 us).
      void SpanConstraint(double h)
      {
         constexpr int kPasses = 4;
         const double rho = set.radius;
         auto invMass = [](const BandChain &ch, int k) {   // point k of a chain (0 = A, Count() + 1 = B)
            const int id = ch.PtNode(k);
            return (id < 0 || ch.nodes[id].pinned) ? 0.0 : 1.0 / ch.nodes[id].m;
         };
         for (int pass = 0; pass < kPasses; pass++)
         {
            bool any = false;
            for (size_t a = 0; a < chains.size(); a++)
               for (size_t b = a + 1; b < chains.size(); b++)
               {
                  BandChain &X = chains[a], &Y = chains[b];
                  if (X.Count() == 0 && Y.Count() == 0)
                     continue;   // (both at rest on their anchors)
                  double vMax = 0.0;   // (the boxes, widened by what the fastest node can cover in the step)
                  V3 xl = X.PtPos(0), xh = xl, yl = Y.PtPos(0), yh = yl;
                  for (int k = 1; k <= X.Count() + 1; k++)
                  {
                     const V3 q = X.PtPos(k);
                     xl = V3(std::min(xl.x, q.x), std::min(xl.y, q.y), std::min(xl.z, q.z));
                     xh = V3(std::max(xh.x, q.x), std::max(xh.y, q.y), std::max(xh.z, q.z));
                     vMax = std::max(vMax, X.PtVel(k).Length());
                  }
                  for (int k = 1; k <= Y.Count() + 1; k++)
                  {
                     const V3 q = Y.PtPos(k);
                     yl = V3(std::min(yl.x, q.x), std::min(yl.y, q.y), std::min(yl.z, q.z));
                     yh = V3(std::max(yh.x, q.x), std::max(yh.y, q.y), std::max(yh.z, q.z));
                     vMax = std::max(vMax, Y.PtVel(k).Length());
                  }
                  const double g = 2.0 * rho + 2.0 * vMax * h;
                  if (xl.x > yh.x + g || yl.x > xh.x + g || xl.y > yh.y + g || yl.y > xh.y + g || xl.z > yh.z + g || yl.z > xh.z + g)
                     continue;
                  for (int j = 0; j <= X.Count(); j++)
                     for (int k = 0; k <= Y.Count(); k++)
                     {
                        double s, t;
                        const V3 xa = X.PtPos(j), xb = X.PtPos(j + 1), ya = Y.PtPos(k), yb = Y.PtPos(k + 1);
                        const double d = RubberContact::SegSegClosest(xa, xb, ya, yb, s, t);
                        if (d < 1e-9)
                           continue;
                        V3 n = ((xa + (xb - xa) * s) - (ya + (yb - ya) * t)) / d;   // from Y toward X
                        if (X.lockZ)
                        {
                           n.z = 0.0;   // (the nodes keep their height)
                           const double l = n.Length();
                           if (l < 1e-9)
                              continue;
                           n = n / l;
                        }
                        const double gap = d - 2.0 * rho;
                        const V3 vx = X.PtVel(j) * (1.0 - s) + X.PtVel(j + 1) * s, vy = Y.PtVel(k) * (1.0 - t) + Y.PtVel(k + 1) * t;
                        const double vn = (vx - vy).Dot(n);                   // > 0: moving apart
                        const double vMin = gap > 0.0 ? -gap / h : 0.0;       // the most closing speed the step allows
                        if (vn >= vMin)
                           continue;
                        const double wx0 = invMass(X, j) * (1.0 - s) * (1.0 - s), wx1 = invMass(X, j + 1) * s * s;
                        const double wy0 = invMass(Y, k) * (1.0 - t) * (1.0 - t), wy1 = invMass(Y, k + 1) * t * t;
                        const double W = wx0 + wx1 + wy0 + wy1;
                        if (W <= 0.0)
                           continue;
                        const double J = (vMin - vn) / W;   // impulse on X along n (on Y the opposite)
                        auto push = [&](BandChain &ch, int pt, double share) {
                           const int id = ch.PtNode(pt);
                           if (id >= 0 && !ch.nodes[id].pinned)
                              ch.nodes[id].vp = ch.nodes[id].vp + n * (J * share / ch.nodes[id].m);
                        };
                        push(X, j, 1.0 - s);
                        push(X, j + 1, s);
                        push(Y, k, -(1.0 - t));
                        push(Y, k + 1, -t);
                        any = true;
                     }
               }
            if (!any)
               break;
         }
      }

      // Advance the band and the plate by h with the balls' reactions (fBall) held; zeta >= 0 replaces the nodes' damping ratio (the ring-down). Returns the band steps taken.
      int Advance(double h, double zeta)
      {
         int steps = 0;
         double left = h;
         while (left > 1e-15)
         {
            // the face about to reach the cord: the band needs the nodes at its ends
            if (hasPlate && pg.valid)
            {
               const double vFace = FaceDN(plate.st.phi) * plate.st.omega;
               if (FaceN(plate.st.phi) + set.radius + std::max(vFace, 0.0) * std::min(left, set.plateMaxStepS) >= 0.0)
                  EnsurePlateNodes();
            }
            if (!active && !hasPlate)
               break;
            double kFaceNode = 0.0;
            fWork.resize(chains.size());
            for (size_t c = 0; c < chains.size(); c++)
               fWork[c] = fBall[c];
            const double torque = hasPlate ? FaceForces(fWork, kFaceNode) : 0.0;
            double hb = left;
            if (hasPlate)
               hb = std::min(hb, set.plateMaxStepS);
            if (active)
               hb = std::min(hb, NodeStepLimit(hb, zeta, kFaceNode));
            if (left - hb < 1e-9 * h)
               hb = left;
            const double phi0 = plate.st.phi;
            if (hasPlate)
            {
               plate.Step(hb, driver.Coil(), torque);
               driver.Tick(hb);
            }
            if (active)
            {
               // (detail books: the positions and the band's energy before the step)
               double e0 = 0.0;
               bool nearStop = false;
               if (books.detail)
               {
                  e0 = DeformEnergy();
                  pPrev.resize(chains.size());
                  for (size_t c = 0; c < chains.size(); c++)
                  {
                     pPrev[c].resize(chains[c].nodes.size());
                     for (size_t id = 0; id < chains[c].nodes.size(); id++)
                        pPrev[c][id] = chains[c].nodes[id].P;
                     for (const int id : chains[c].order)
                        for (const RubberContact::StopCircle &s : chains[c].stops)
                           nearStop = nearStop || std::hypot(chains[c].nodes[id].P.x - s.x, chains[c].nodes[id].P.y - s.y) <= s.r + 0.02;
                  }
               }
               // the new velocities, the spans meeting (a constraint on those velocities), then the nodes move
               for (size_t c = 0; c < chains.size(); c++)
                  chains[c].StepVelocities(fWork[c], hb, zeta);
               const double keSpan0 = books.detail ? NodeKinetic() : 0.0;
               SpanConstraint(hb);
               const double wCC = books.detail ? NodeKinetic() - keSpan0 : 0.0;
               for (size_t c = 0; c < chains.size(); c++)
               {
                  chains[c].StepPositions(hb);
                  chains[c].ApplyStops();
               }
               const int pc = PlateChain();
               double wFace = 0.0;
               for (const FaceLoad &fl : faceLoads)
               {
                  const BandChain::Node &nd = chains[pc].nodes[fl.id];
                  const double xn = (nd.P.x - nd.P0.x) * pg.nx + (nd.P.y - nd.P0.y) * pg.ny;
                  wFace += fl.F * (xn - fl.xn0);
               }
               books.faceOnNodesJ += wFace;
               books.faceOnPlateJ += torque * (plate.st.phi - phi0);
               if (books.detail)
               {
                  double wBall = 0.0;
                  for (size_t c = 0; c < chains.size(); c++)
                     for (size_t id = 0; id < pPrev[c].size(); id++)
                        if (c < fBall.size() && id < fBall[c].size())
                           wBall += fBall[c][id].Dot(chains[c].nodes[id].P - pPrev[c][id]);
                  books.chainContactJ += wCC;
                  const double resid = DeformEnergy() - e0 - wBall - wCC - wFace;
                  books.residualJ += resid;
                  if (nearStop)
                     books.residualStopJ += resid;
               }
            }
            left -= hb;
            steps++;
         }
         books.bandSteps += steps;
         return steps;
      }

      // After a step: the drawn shape, the stretch seen, idle nodes merged away (those no ball loads and the face does not press, lying on the line between their neighbours)
      void AfterSteps()
      {
         if (!active)
            return;
         SyncPts();
         for (const BandChain &ch : chains)
            books.maxStretch = std::max(books.maxStretch, MaxLegStretch(ch));
         if (!set.follow || FixedNodes())
            return;
         const int pc = PlateChain();
         bool merged = false;
         const double e0 = books.detail ? DeformEnergy() : 0.0;
         for (size_t c = 0; c < chains.size(); c++)
         {
            BandChain &ch = chains[c];
            loadedBuf.assign(ch.nodes.size(), 0);
            guardBuf.clear();
            for (const Touch &t : touches)
            {
               if (t.chain < 0 || t.chain >= (int)chains.size())
                  continue;
               const BandChain &tc = chains[t.chain];
               int near = -1;
               double best = 1e30;
               for (const int id : tc.order)
                  if (std::fabs(tc.nodes[id].f - t.f) < best)
                  {
                     best = std::fabs(tc.nodes[id].f - t.f);
                     near = id;
                  }
               if (near < 0)
                  continue;
               if ((size_t)t.chain == c)
                  loadedBuf[near] = 1;
               else
                  guardBuf.push_back(tc.nodes[near].P);
            }
            if ((int)c == pc && FaceN(plate.st.phi) + set.radius > -set.followM)
               for (const int id : ch.order)
                  if (InFootprint(pc, id))
                     loadedBuf[id] = 1;
            merged = ch.MergeIdle(loadedBuf, std::max(set.restM, 0.3e-3), 3.0, guardBuf, 0.03) || merged;
         }
         if (merged)
            Rebuild();
         if (books.detail)
            books.mergeJ += DeformEnergy() - e0;
      }
      std::vector<char> loadedBuf;
      std::vector<V3> guardBuf;

      // How far the band is pushed toward a blade at a switch, m (0 while it is at rest)
      double SwitchDisp(const SwitchSpec &ss) const
      {
         if (!active || ss.span < 0)
            return 0.0;
         const int ci = ChainOfSpan(ss.span);
         if (ci < 0 || chains[ci].Count() == 0)
            return 0.0;
         const BandChain &ch = chains[ci];
         int j = 0;
         while (j < ch.Count() && ch.PtF(j + 1) < ss.f)
            j++;
         const double f0 = ch.PtF(j), f1 = ch.PtF(j + 1);
         const double u = f1 > f0 ? std::min(std::max((ss.f - f0) / (f1 - f0), 0.0), 1.0) : 0.0;
         const V3 P = ch.PtPos(j) + (ch.PtPos(j + 1) - ch.PtPos(j)) * u;
         const V3 rest = ch.A + (ch.B - ch.A) * ss.f;
         return (P.x - rest.x) * ss.bx + (P.y - rest.y) * ss.by;
      }
   };
}
