// license:GPLv3+
//
// A moving part of the contact system: the engine-facing interface of anything that has its own bodies and surfaces, is stepped inside the islands of the balls that touch it, and
// raises its own events. Gates and spinners (HingeModel::Engine, contacthingepart.h) and the jet bumper (JetBumper::EnginePart, contactbumperpart.h) implement it; the contact system
// (contactsystem.cpp) and the host (ContactPartHost, contactparthost.h) know only this interface. The turntable disc (DiscModel::Engine, contactdiscpart.h) implements it too. A new part (a kicker door, a drop target, a plunger) implements it and nothing else changes.
//
// What the engine does with a part, in the order of a physics cycle:
//   SetGravity       the gravity the part feels (gravity minus the cabinet acceleration, so a nudge acts on it), once per cycle
//   Probe            where is the nearest surface of the part to a ball (distance, normal, curvature): engagement, release, the arrival-time rule
//   AdvanceAlone     a cycle with no ball near: the part moves by itself (a part at rest sleeps)
//   events           what happened (Hit, switch closed / opened, stops, spins), in order; the host fires them on the table part after the step
// Pure SI maths with no VPX types: the same code runs in the offline tests.

#pragma once

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "physics/rubbercontact.h"
#include "physics/contactsurface.h"
#include "physics/groupstep.h"
#include "physics/partrecords.h"

// [Contact] target_patches (P3b piece 1, to-do 1.33; runs/2026-11-23-targets-patch): a mesh that translates along one unit axis a by s (a drop target's z, a stand-up's travel),
// probed in the mesh's rest frame: the ball's centre moved by -s a and the fixed mesh probed there. A pressed place: the shared static mesh path (BarsTriangles), so every pressed
// place of a non-convex mesh pushes (the merged contact over its patches, the level form's peaks for peak_slope). Nothing pressed (or target_patches 0): the nearest face
// (ProbeTriangles over every face), its true gap and normal, which approach, arrival and the touch prediction read. A pure translation leaves the normal as it is: one column, the
// depth's slope along s is n . a (times Hit::grad, the group's rule), the point's speed a per unit rate, no spin. Links and folds built once at load (RubberContact::LinkFaces). Pure SI
// patch_friction: a pressed probe's patch list (plist, the probe's BarsTriangles / BarsTube patchesOut) and each patch's pieces (RubberContact::PatchPieces) into the hit (npat, pat,
// npc, pc, pd). One copy: the static surfaces (ContactSystem's SurfaceShape::ProbeSeg) and the translating mesh (SlideMesh::Probe)
inline void HitPatches(const std::vector<RubberContact::TubeProbe> &plist, GroupStep::Hit &h)
{
   const RubberContact::PatchPieceList &pl = RubberContact::PatchPieces();
   const bool withIds = pl.off.size() == plist.size() + 1;   // (a path that names no pieces: the patches carry none)
   h.npat = (int)plist.size();
   h.pat.Fit(h.npat);
   h.npc = withIds ? (int)pl.id.size() : 0;
   h.pc.Fit(h.npc);
   h.pd.Fit(h.npc);
   const bool withPd = pl.pd.size() == pl.id.size();   // (level_stiffness: each piece's data; a path without it: none, depth -1)
   for (int j = 0; j < h.npc; j++)
   {
      h.pc[j] = pl.id[j];
      h.pd[j] = withPd ? pl.pd[j] : RubberContact::PieceData();
   }
   for (int j = 0; j < h.npat; j++)
   {
      GroupStep::HitPatch &q = h.pat[j];
      const RubberContact::TubeProbe &t = plist[j];
      q.d = t.overlap ? t.delta : -1.0;
      q.n = t.normal;
      q.grad = t.grad;
      q.p0 = withIds ? pl.off[j] : 0;
      q.pn = withIds ? pl.off[j + 1] - pl.off[j] : 0;
      const bool withGeo = withIds && pl.geo.size() == plist.size();
      q.gk = withGeo ? pl.geo[j].kind : 0;
      q.ga = withGeo ? pl.geo[j].a : GroupStep::V3();
      q.gb = withGeo ? pl.geo[j].b : GroupStep::V3();
   }
}

struct SlideMesh
{
   using V3 = RubberContact::V3;
   // [Contact] target_patches, written by ContactSettings' Publish (play start and every live reload), read by every probe
   static bool &PatchesOn()
   {
      static bool on = true;
      return on;
   }
   // [Contact] patch_friction (P3b piece 2) with target_patches: a pressed probe also writes its patch list (Hit::npat, pat, pc, pd), so the group step takes the target's friction per
   // patch group. Its geometry (HitPatch ga, gb; PieceData's points) is in the mesh's rest frame, as the probe's centre cr = c - s a: the group step hands its consumers (FrictionGroups'
   // plane, BLinGroup's piece depths, the link windows) the ball's centres in the same frame (GroupStep::Shape::RestCentre). Written by Publish beside PatchesOn
   static bool &PatchFrictionOn()
   {
      static bool on = true;
      return on;
   }
   // patch_zero (P3b piece 4): the mesh's probe takes per-patch zeros (RubberContact::PatchZeros, read by BarsTriangles) only on the patch path: target_patches and the links
   // built. The face ids are rest-frame ids, fixed in the mesh under its translation (report.md: patch depths at (c, s) = the static mesh's at c - s a)
   bool PatchZeroOn() const { return PatchesOn() && nb.size() == tris.size(); }
   std::vector<RubberContact::SurfaceTri> tris;    // the mesh at its rest pose (world, m)
   std::vector<std::array<int, 3>> vid;            // per face: its welded vertex ids (the static weld rule, RubberSurface::Add)
   std::vector<std::array<int, 3>> nb;             // per face: the face across each edge, -1 none (LinkFaces)
   std::vector<std::array<int8_t, 3>> fold;        // per face: each shared edge's fold (LinkFaces)
   std::vector<V3> lo, hi;                         // per face: its box (the near pass)
   std::vector<int> all;                           // every face (ProbeTriangles' list)
   V3 a;                                           // the translation axis, unit
   V3 nFar;                                        // the normal written with no face at all (as the target's own probe did)
   double edgeRho = 1e-3;                          // the mesh edges' radius for the contact, m

   void Build(const std::vector<RubberContact::SurfaceTri> &faces, const std::vector<std::array<int, 3>> &ids, const V3 &axis, const V3 &farNormal, const double edgeRhoM)
   {
      tris = faces;
      vid = ids;
      a = axis / axis.Length();
      nFar = farNormal;
      edgeRho = edgeRhoM;
      all.resize(tris.size());
      for (size_t i = 0; i < tris.size(); i++)
         all[i] = (int)i;
      if (vid.size() == tris.size())
         RubberContact::LinkFaces(tris, vid, nb, fold);
      lo.resize(tris.size());
      hi.resize(tris.size());
      for (size_t i = 0; i < tris.size(); i++)
      {
         const RubberContact::SurfaceTri &t = tris[i];
         lo[i] = V3(std::min({ t.a.x, t.b.x, t.c.x }), std::min({ t.a.y, t.b.y, t.c.y }), std::min({ t.a.z, t.b.z, t.c.z }));
         hi[i] = V3(std::max({ t.a.x, t.b.x, t.c.x }), std::max({ t.a.y, t.b.y, t.c.y }), std::max({ t.a.z, t.b.z, t.c.z }));
      }
   }

   // The probe at shift s. Near pass: the faces whose box is within R + edgeRho of the shifted centre cr. Exact: BarsTriangles reads only faces within R of the centre (kind 0 and
   // plane distance below R - zero, zero 0 here: the static path passes it only its near pass, contactsurfacehost.cpp ProbePatches); every point of a face lies in its box, so a face
   // left out is farther than R + edgeRho from cr. Pressed: the merged answer, the level form's peaks always (Hit::np, pk: the solver's peak_slope decides their use)
   void Probe(const V3 &c, const double R, const double s, GroupStep::Hit &h) const
   {
      const V3 cr = c - a * s;
      if (PatchesOn() && nb.size() == tris.size())
      {
         const double reach = R + edgeRho, reach2 = reach * reach;
         thread_local std::vector<int> cand;   // (capacity reused: no allocation in the probe)
         cand.clear();
         for (int i = 0; i < (int)tris.size(); i++)
         {
            const double dx = std::max({ lo[i].x - cr.x, 0.0, cr.x - hi[i].x }), dy = std::max({ lo[i].y - cr.y, 0.0, cr.y - hi[i].y }),
                         dz = std::max({ lo[i].z - cr.z, 0.0, cr.z - hi[i].z });
            if (dx * dx + dy * dy + dz * dz <= reach2)
               cand.push_back(i);
         }
         if (!cand.empty())
         {
            RubberContact::BallSI b;
            b.R = R;
            RubberContact::PeakTerms().clear();
            const bool pf = PatchFrictionOn();
            thread_local std::vector<RubberContact::TubeProbe> plist;   // (patch_friction: the patch list; capacity reused)
            if (pf)
               RubberContact::PatchPieces().Clear();   // (as the static ProbePatches: a path that writes no pieces leaves the list short; HitPatches checks)
            // (zero 0: no scalar zero in the probe, ZeroHit takes it off after. patch_zero: PatchZeros() set, BarsTriangles shifts each raw patch by its own zero; a place
            // raw-pressed with nothing past its zero still returns overlap with its shifted depth (<= 0) and its patches, so the answer below is the static ProbeSeg's; Nearest
            // only when nothing is raw-pressed)
            const RubberContact::TubeProbe pr = RubberContact::BarsTriangles(b, tris, vid, nb, cand, edgeRho, cr, 0.0, nullptr, &fold, pf ? &plist : nullptr);
            if (pr.seg >= 0 && pr.overlap)
            {
               h.nc = 1;
               h.npat = 0;
               h.npc = 0;
               if (pf)
                  HitPatches(plist, h);   // (rest-frame geometry: see PatchFrictionOn)
               h.delta = pr.delta;
               h.n = pr.normal;
               h.grad = pr.grad;
               h.reff = pr.reff > 0.0 ? pr.reff : R;
               h.curvDist = pr.reff > 0.0 && pr.reff < 0.999 * R ? pr.dist : 0.0;   // (as the fixed mesh: an edge the ball rolls round)
               h.c[0].Set(0, pr.normal.Dot(a), a, V3());
               const std::vector<RubberContact::PeakTerm> &pt = RubberContact::PeakTerms();   // (as ContactSystem's ProbeSeg for a static surface)
               h.np = (int)pt.size();
               h.pk.Fit(h.np);
               for (int j = 0; j < h.np; j++)
                  h.pk[j] = pt[j];
               if (h.np <= 1 && GroupStep::PeakSlopeCheckEnv() && RubberContact::PatchZeros() == nullptr)   // (raw answers only: Nearest takes no zeros)
                  CheckOnePeak(cr, R, h);
               return;
            }
         }
      }
      Nearest(cr, R, h);
   }

   // link_windows (rule T, P3b piece 3) on the relative segment: c0r -> c1r the ball's centre in the rest frame at the step's start and end (GroupStep::Shape::RestCentre; a straight
   // segment, the part's move uniform over the step: runs/2026-11-23-targets-patch report.md Q2). The shared mesh kit (RubberContact::LinkWindowsMesh) on the faces whose box is within
   // R + |c1r - c0r| / 2 of the segment's middle: exact (an edge point within R of a point of the segment is within R + half its length of the middle; the box holds the face), reach R
   // (no zero, as the probe). Window times are fractions of the step: no shift back. Scratch per thread
   void LinkWindows(const V3 &c0r, const V3 &c1r, const double R, std::vector<RubberContact::LinkWindow> &out) const
   {
      out.clear();
      if (!PatchesOn() || nb.size() != tris.size() || R <= 0.0)
         return;
      const V3 mid = (c0r + c1r) * 0.5;
      const double reach = R + 0.5 * (c1r - c0r).Length(), reach2 = reach * reach;
      thread_local std::vector<int> cand;
      cand.clear();
      for (int i = 0; i < (int)tris.size(); i++)
      {
         const double dx = std::max({ lo[i].x - mid.x, 0.0, mid.x - hi[i].x }), dy = std::max({ lo[i].y - mid.y, 0.0, mid.y - hi[i].y }),
                      dz = std::max({ lo[i].z - mid.z, 0.0, mid.z - hi[i].z });
         if (dx * dx + dy * dy + dz * dz <= reach2)
            cand.push_back(i);
      }
      RubberContact::LinkWindowsMesh(tris, nb, cand, c0r, c1r, R, out);
   }

   // Today's target probe (one copy): the nearest face over every face at cr
   void Nearest(const V3 &cr, const double R, GroupStep::Hit &h) const
   {
      RubberContact::BallSI b;
      b.R = R;
      const RubberContact::TubeProbe pr = RubberContact::ProbeTriangles(b, tris, all, edgeRho, cr);
      h.nc = 1;
      h.npat = 0;   // (nothing pressed: no patches)
      h.npc = 0;
      if (pr.seg < 0)
      {
         h.delta = -1.0;
         h.n = nFar;
         h.reff = R;
         h.curvDist = 0.0;
         h.grad = 1.0;
         h.c[0].Set(0, -1.0, a, V3());
         return;
      }
      h.delta = pr.delta;
      h.n = pr.normal;
      h.grad = pr.grad;
      h.reff = pr.reff > 0.0 ? pr.reff : R;
      h.curvDist = pr.reff > 0.0 && pr.reff < 0.999 * R ? pr.dist : 0.0;   // (as the fixed mesh: an edge the ball rolls round)
      h.c[0].Set(0, pr.normal.Dot(a), a, V3());
   }

   // check mode (VPX_PATCH_CHECK=1, diagnostics only; no result changes): a pressed answer with one peak or none against today's nearest face at the same cr. Mismatch: delta off by
   // more than 1e-12 R, a normal component or grad by more than 1e-12, reff by more than 1e-12 R. Running totals at 1, 2, 4, ... calls and mismatches, the first 5 mismatches' values
   void CheckOnePeak(const V3 &cr, const double R, const GroupStep::Hit &h) const
   {
      static std::atomic<long long> nChk { 0 }, nBad { 0 };
      GroupStep::Hit t;
      Nearest(cr, R, t);
      const bool bad = std::abs(h.delta - t.delta) > 1e-12 * R || std::abs(h.n.x - t.n.x) > 1e-12 || std::abs(h.n.y - t.n.y) > 1e-12 || std::abs(h.n.z - t.n.z) > 1e-12
                    || std::abs(h.grad - t.grad) > 1e-12 || std::abs(h.reff - t.reff) > 1e-12 * R;
      const long long c = ++nChk, b = bad ? ++nBad : nBad.load();
      if ((c & (c - 1)) == 0)
         PLOGI << "target_patches check: " << c << " pressed probes with one peak or none, " << b << " not today's nearest face";
      if (bad && b <= 5)
         PLOGI << "target_patches check miss " << b << ": delta " << h.delta * 1e6 << " / " << t.delta * 1e6 << " um, n (" << h.n.x << ", " << h.n.y << ", " << h.n.z << ") / ("
               << t.n.x << ", " << t.n.y << ", " << t.n.z << "), grad " << h.grad << " / " << t.grad << ", reff " << h.reff * 1e3 << " / " << t.reff * 1e3 << " mm";
   }
};

// What the engine keeps of every moving part between its calls (PartRecords: part of its state)
struct MovingPartState
{
   double steppedS = 0.0;           // how long the contact system stepped the part with its balls in this cycle (s); the host advances it alone for the rest
   bool steppedThisCycle = false;   // the contact system stepped the part with the balls of an island this cycle (the host's own advance skips it once)
   double eOwnJ = -1.0;             // energy skip test: what it can hold this cycle before a ball touches it (its kinetic energy and drive work, EnergyMove with no feed, plus squashJ), J; -1 = not set (ContactSystem ReachEnergy)
   double squashJ = 0.0;            // the energy reach: what its contacts hold at the cycle's end, J: each ball-part contact's stored energy by the contact law and its compliant stops' squash (GroupStep::Part::StopSquashJ); fed to the balls that can reach it (ContactSystem ReachEnergy)
   bool restTouched = false;        // a ball asleep in a pile touches it (the engine's pile sleep): the part does not sleep while one does (GroupStep::RestWatch)
};

class ContactMovingPart : public PartRecords, public MovingPartState
{
public:
   using V3 = RubberContact::V3;

   struct Event
   {
      enum Kind { Hit, SwitchClosed, SwitchOpened, LimitEOS, LimitBOS, Spin, BandHit, SlingFire, Dropped, Raised } kind;   // (Dropped, Raised: a drop target, droptargetsimple.h)
      int ball = -1;        // Hit, BandHit: the ball's index as the owner numbers them (a drop target's Hit: the ball's own id, Ball::m_id); a band's SwitchClosed / SwitchOpened: the blade's index
      double value = 0.0;   // Hit: normal approach speed (m/s); Limit: angular speed (deg/s); BandHit: the cord's largest deflection under the ball (m); SlingFire: 1 = the script's trigger
   };

   // A saved copy of a part's dynamic state, for rewinding an island (docs/unified-physics/archive/UnifiedContactDesign.txt section 11, point 2): plain bytes, no allocation. Its state records (PartRecords)
   // and its group part's, whole: the values that change while it runs, not its settings or shape.
   struct State
   {
      alignas(8) unsigned char bytes[2048];
      size_t used = 0;
   };

   std::vector<Event> events;

   ContactMovingPart()
   {
      KeepState<MovingPartState>(this);
      KeepQueue(&events);
   }
   virtual ~ContactMovingPart() = default;

   virtual void SetBall(const RubberContact::BallSI &b) = 0;   // the ball of the present call (balls may differ in size and mass)
   virtual void SetGravity(const V3 &g) = 0;
   virtual V3 Origin() const = 0;                              // the part's reference point (world, m): the host's box and the log
   virtual double Reach() const = 0;                           // the swept reach about it (m)
   virtual double MotionRange() const { return 2.0 * Reach(); }   // the farthest any point of its surface can move from where it is, over its whole range of motion (m); every moving point is within Reach of Origin
   virtual RubberContact::TubeProbe Probe(const V3 &c) const = 0;
   virtual void SurfaceVelocity(const V3 &point, V3 &v, V3 &w) const = 0;
   virtual int AdvanceAlone(double dt) = 0;
   // Save the dynamic state / put it back exactly: its records and its group part's, copied whole (one shared code for every part). Restoring cuts each event queue back to what it held at the save
   // (the replay makes the later ones again) and rebuilds what is derived from the state.
   void SaveState(State &out)
   {
      GroupStep::Part *const g = GroupPart();
      const size_t a = StateBytes(), b = g != nullptr ? g->StateBytes() : 0;
      assert(a + b <= sizeof out.bytes);
      SaveRecords(out.bytes);
      if (g != nullptr)
         g->SaveRecords(out.bytes + a);
      out.used = a + b;
   }
   void RestoreState(const State &in)
   {
      GroupStep::Part *const g = GroupPart();
      RestoreRecords(in.bytes);
      if (g != nullptr)
      {
         g->RestoreRecords(in.bytes + StateBytes());
         g->AfterRestore();
      }
   }
   // The fastest any point of the part's surfaces can move ALONG ITS OWN NORMAL (toward a ball) within the next cycle of cycleS seconds (m/s): its present speed plus what its own drive and springs
   // can add in a cycle. Balls are not counted (their contact is the contact). The engine adds it to the ball's closing speed when it decides to engage a ball and when it decides to let one go
   // (UnifiedContactDesign.txt 11 point 1, 12.1 d). A part may answer its whole point speed (an upper bound); a disc turning about its own axis answers 0 (its surfaces slide along themselves).
   virtual double SurfaceSpeedBound(double cycleS) const = 0;
   // the energy reach (design_reach.md sections 1, 4b; rubbercontact.h EnergyReach): the most any point of the part's surface can move relative to a point c over a cycle of T
   // seconds (m), its turn or travel bounded by energy: its own motion, its drives (gravity, springs, blades, cams, coil) and E, the kinetic energy (J) of the balls that can reach it
   // this cycle. work (if not null): what the part can hand to the balls in the cycle beyond E, its kinetic energy now plus its drives' work over that move (J). Each part its own law
   virtual double EnergyMove(const V3 &c, double E, double T, double *work) const = 0;
   // the one-pass energy bound (runs/2026-11-20-energy-bound/design.md, AGREED 2026-11-20 items 1-6, 4b, 5b): the most the part can hand to the balls over
   // a cycle of T seconds whatever they feed it (J), from its state now and its built data, no move: its energy in the group's book now (GroupStep::Part::Energy: its kinetic energy,
   // its springs, blades, cams and weight) less the least that book's stored energy takes over the poses it can reach between its stops (item 4b: one source, no per-spring formula),
   // + its drives over T (a coil's field energy now + V^2 T / R, a motor's largest power x T). Past its stops under the force that presses it there (weight, preload): a linear stop
   // adds aLin = F sqrt(2 / k) (J^0.5, item 5: E = E0 + aLin sqrt(E)), a Hertz stop (k d^1.5) aHertz = F (5 / (2 k))^0.4 (J^0.6, item 5b: E = E0 + aHertz E^0.4); the caller closes
   // them in one step (ContactSystem OnePassE). The squash its contacts and compliant stops hold now is squashJ, not here
   virtual double EnergyCap(double T, double &aLin, double &aHertz) const = 0;
   // between cycles (the host's advance, after the cycle): what starts at the next cycle's start (a jet's pending pulse: a jet's pulse at a cycle's start); cycleS: that cycle's length
   virtual void NextCycleStart(double /*cycleS*/) {}
   // the energy reach: a lower bound on the distance from c to the part at any pose it can take in the cycle of T seconds with E (J) fed to it (m; its region, design_reach.md
   // section 3 in its bounding form; a part whose region does not depend on its feed ignores E); -INFINITY = none (the move above alone). Grouping takes the larger of this less R
   // and the gap now less EnergyMove
   virtual double CycleRegionDist(const V3 &, double /*T*/, double /*E*/) const { return -INFINITY; }
   // the part cull (runs/2026-11-20-part-cull/design.md 1c, 2R, 2J): CullCycle once per cycle, before the grouping, with E the table's whole energy this cycle (J:
   // every ball's, every stored and every part's work, an upper bound on any group's feed) and the cycle T (s); then CullDist(c) is a lower bound on the distance (m) from c to any
   // touchable point of the part at any pose it can take in the cycle (the ball's radius not included). -INFINITY = no region: a candidate for every ball. Default: the
   // grouping's own region at the table's energy (CycleRegionDist: the bat, the hinge parts, the disc; the jet and the target precompute the same region once per cycle)
   virtual void CullCycle(const double E, const double T) { cullE = E; cullT = T; }
   virtual double CullDist(const V3 &c) const { return CycleRegionDist(c, cullT, cullE); }
   virtual std::string RegionNote(double /*E*/, double /*T*/) const { return std::string(); }   // (the log, once per table start: how its region is built; empty = nothing to tell)
   double cullT = 0.0, cullE = 0.0;   // (the cycle and the table's energy CullCycle was set for, s and J; not part of the state: set every cycle)
   // the part's state as a short text for the log (an engagement of a ball with a part that is moving is the case to look at)
   virtual std::string StateNote() const { return std::string(); }
   // A part on the group step (docs/unified-physics/GroupStepDesign.txt 18): the engine moves it with the balls in one solve and holds its shapes. Null: the part is on the old loop until it moves over. No part implements these yet.
   virtual GroupStep::Part *GroupPart() { return nullptr; }
   virtual int GroupSurfaces() const { return 1; }                                                   // its surfaces a ball can meet: one shape per ball and surface, each with its own contact and state
   virtual void GroupProbe(int, const V3 &, double, const double *, GroupStep::Hit &) const {}       // a surface for a ball centre and radius, its coordinates moved by trial (null = where they are); fills the hit's columns
   // a translating mesh (P3b pieces 2-3, SlideMesh): the ball's centre c in the surface's rest frame with the part's coordinates moved by move (null = where they are): c - s a; false =
   // no rest frame (the surface turns, or the part has none): its patches and windows stay world. And the link windows on the rest-frame segment c0r -> c1r (SlideMesh::LinkWindows)
   virtual bool GroupRestCentre(int, const V3 &, const double *, V3 &) const { return false; }
   virtual void GroupLinkWindows(int, const V3 &, const V3 &, double, std::vector<RubberContact::LinkWindow> &out) const { out.clear(); }
   // patch_zero on a translating mesh (P3b piece 4, SlideMesh::PatchZeroOn): true = the surface keeps its found-inside zero per patch, keyed by face id in its rest frame (the
   // patches move with the mesh); false = one scalar zero (every other part)
   virtual bool GroupPatchZero(int) const { return false; }
   virtual bool GroupCycleNear(int, const V3 &, double, double, std::vector<int> &, int &) const { return false; }              // the part finder: the cycle's list (GroupStep::Shape::CycleNear)
   virtual int GroupCycleMisses(int, const V3 &, double, double, const std::vector<int> &) const { return 0; }                // (GroupStep::Shape::CycleMisses, check mode)
   // the part finder, the overlap check of free balls (ContactSystem CheckOverlaps): a lower bound on the distance from c to the part where it is now; each part its bounding
   // radius from its built pieces (group A: not Reach(), which carries margins). Here: none, always probed
   virtual double RegionDist(const V3 &) const { return -INFINITY; }
   virtual const RubberContact::RubberParams *GroupLaw(int, double &) const { return nullptr; }      // the surface's material and layer; null = the host's
   // the path test (GroupStep::Shape::PathReach): a surface of the part that never moves (a bumper's body, a disc turning about its own axis), the ball's centre moving straight
   // from s to s + d: 1 = it can touch it on the way, 0 = no (exact), -1 = no exact test (a surface that moves with the part: build 3, piece 6), FarStillFar probes as before
   virtual int GroupPathReach(int, const V3 &, const V3 &, double) const { return -1; }
   // the part finder, a surface that moves with the part (group A; GroupStep::Shape::PathMotion): the most any point of it moves relative to P over the step, the ball's centre
   // moving dLen and the part by trial (the solve's end move; null = no move); its probe (GroupProbe) must give the exact distance at any range. -1 = none: the two probes
   virtual double GroupPathMotion(int, const V3 &, double, const double *) const { return -1.0; }
   virtual double GroupPathSlope(int) const { return 1.0; }   // its probe's largest slope (GroupStep::Shape::PathSlope), built at load
   virtual double GroupPathCoverSag(int, int &) const { return -1.0; }   // the cover's largest sagitta and kind (the self-check, VPX_FINDER_CHECK; GroupStep::Shape::PathCoverSag); -1 = none
   virtual void CycleBegin() {}                                                                      // a group holding it is being built this cycle
   virtual void TouchBegan(int, int, int, double, const V3 &) {}                                     // this surface, this ball of the group (its index; its own id, Ball::m_id), closing at vn (m/s), the touch's normal (toward the ball): a touch began. A part's Hit is what it decides
   virtual bool Sleeping() const { return false; }                                                   // no ball near and nothing to step: no group this cycle
   virtual bool RestNumbers(double &, double &) const { return false; }                             // pile sleep's log: what AtRest is decided on, the move over the last cycle in a group and the surface speed (m, m/s); false = not given
   virtual bool AtRest() const { return false; }                                                     // nothing driven and nothing moved over the last whole cycle: it stays out of a ball's solve while no ball can reach it (a part that does not say so is always in)
   virtual GroupStep::RestWatch *Rest() { return nullptr; }                                          // its rest watch (GroupStep::RestWatch): the group restarts it where it is built, with RestPose
   virtual GroupStep::RestWatch::Pose RestPose() const { return {}; }                                // its one rest pose, alone and in a group: its coordinates and internal states
};
