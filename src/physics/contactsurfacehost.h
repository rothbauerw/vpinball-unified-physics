// license:GPLv3+

#pragma once

#include <filesystem>
#include <chrono>
#include <string>
#include <vector>
#include <functional>
#include <set>
#include <memory>
#include <map>
#include <array>
#include <cstdint>

#include "physics/collide.h"
#include "physics/rubbercontact.h"
#include "physics/contactsurface.h"
#include "physics/contacthost.h"
#include "physics/hitsling.h"
struct CordBand;

// Static rubber (post ring, wall rubber) as ONE collider: a closed tube of circular cross-section swept along the rubber's centreline.
// Replaces the triangle/edge/point set of Rubber::PhysicSetup when rubber_contact is on for that rubber (vpx_physics.ini).
// The ball is detected by an ordinary swept test (first touch, exactly as the other colliders). Collide() then solves the whole
// contact in one call: it integrates the force model of physics/rubbercontact.h for the ball alone (SI, small local steps) until the
// ball has left, and hands the exit velocity, spin and position back. Nothing else acts on the ball during the dwell (about 5 ms at
// post speeds), which is the same simplification as the 1 ms engine step.

// One post of a band: the circle the band wraps (centreline radius = post radius + band radius) and the centreline vertices on it
struct BandPost
{
   double x = 0.0, y = 0.0, r = 0.0;   // VPX units
   int first = 0, last = 0;            // centreline vertex run (indices into the closed loop, may wrap past the end)
   bool fitted = true;                 // false = the run was too short to fit a circle
};

// A straight run of the band between two consecutive posts
struct BandSpan
{
   int postA = 0, postB = 0;           // indices into the post list (B follows A around the loop)
   int first = 0, last = 0;            // centreline vertices of the run (inclusive)
   double lengthMm = 0.0;
   bool free = false;                  // long enough to deflect
};

// What scripts can read about a band (Rubber.BandDeflection etc.), refreshed by the physics every cycle the band moves.
struct BandInfo
{
   bool active = false;          // a span is deflected or being loaded
   int span = -1;                // index of the span of the most deflected node
   float deflectionMm = 0.f;     // current deflection of the most deflected node
   float peakMm = 0.f;           // peak deflection of the last hit
   float posSpan = 0.f;          // where that node is along its span: 0 at the first post, 1 at the second
   float posRubber = 0.f;        // where it is along the whole rubber centreline, 0..1 from the first point
   float nodeX = 0.f, nodeY = 0.f, nodeZ = 0.f;   // its position, VPX units
   int hitTimeMs = 0;            // game time when the last ball left the band
   int hitSpan = -1;             // where the last ball loaded the band: span, place along the span, place along the whole rubber
   float hitPosSpan = 0.f, hitPosRubber = 0.f;
   float speedIn = 0.f, speedOut = 0.f;           // that ball's speed when it hit and when it left, m/s
};

// What scripts can read about a sling (Rubber.SlingPlateAngle etc.), refreshed by the physics every cycle.
struct SlingInfo
{
   bool enabled = false;         // the coil may fire
   float plateAngleDeg = 0.f;    // plate angle from rest toward the coil stop
   bool switchClosed = false;    // any switch is closed
   bool coilOn = false;          // a pulse is running
   float coilCurrentA = 0.f;
   int pulses = 0;               // pulses started
   int fireTimeMs = 0;           // game time when the last one started
};

// One ball on the band
struct BandBall
{
   HitBall *ball = nullptr;
   RubberContact::V3 vin;                // its velocity at engagement
};

// The faces of a wall or primitive flagged as rubber, in the world frame (metres). Owned by the ContactSurfaceHost that hosts the object's contact.
struct RubberSurface
{
   std::vector<RubberContact::SurfaceTri> tris;
   float lo[3] = { 0.f, 0.f, 0.f }, hi[3] = { 0.f, 0.f, 0.f };   // bounding box, VPX units
   void Add(const Vertex3Ds &a, const Vertex3Ds &b, const Vertex3Ds &c);   // one face (VPX units)
   // multi_contact 1: which faces share an edge (the ball's contacts on two faces joined by a pressed edge are one contact; RubberContact::BarsTriangles)
   std::vector<std::array<int, 3>> vid;   // per face: its welded vertices (vertices equal in all three float bits are one)
   std::vector<std::array<int, 3>> nb;    // per face: the face across edge e = (vertex e, vertex e+1), -1 none; empty until Link
   // cut B (runs/2026-11-12-mc-cutB): per face and shared edge e, which side of this face's plane (by its normal n) the neighbour's third vertex lies on: +1 / -1; 0 when the fold is within
   // the stored mesh's float32 rounding on both faces (flat); 2 unknown (a zero-area face, or a vertex exactly on the plane outside the band): always taken as concave. Built by Link.
   std::vector<std::array<int8_t, 3>> fold;
   std::map<std::array<uint32_t, 3>, int> weld;   // vertex bits -> id, while faces are added
   int Link();   // builds nb from vid (once, after the last Add); returns the number of edges on more than two faces (each links its first pair only)
   // shared_finder (to-do 1.34, piece P1; runs/2026-11-11-finder/design.md sections 1-2): the finder's index over these faces, built once before the first cycle (ContactSurfaceHost::FinderBuild).
   // Each face's tight box is SurfaceTri lo / hi (the triangle's corners; rounding zero for a mesh). In P1 only check mode (VPX_FINDER_CHECK=1) reads the index.
   struct FinderNode
   {
      RubberContact::V3 lo, hi;     // the box over its items' boxes, m
      int left = -1, right = -1;    // children (inner node); left < 0 = a leaf
      int first = 0, count = 0;     // a leaf's items: order[first .. first + count)
   };
   std::vector<FinderNode> fnodes;  // the box tree over the faces (median split on the longest side, 4 faces per leaf, ties by face index); fnodes[0] is the root
   std::vector<int> forder;         // face indices in the tree's leaf order
   std::vector<int> piece;          // per face: its mesh piece (faces joined through shared edges, nb), numbered by each piece's first face in load order
   int nPieces = 0;
   int fobj = -1;                   // this object's index in the finder's object list (load order), -1 = not in the index
};

// One ball in an engine-stepped contact with a static rubber (post ring, wall rubber): the rubbers it touches and the bookkeeping for the log and the release.
class FlipperMoverObject;

struct BallContact
{
   class ContactHost *host = nullptr;   // the surface or part that engaged the ball (its settings drive the step)
   class ContactSurfaceHost *band = nullptr;       // a band that holds this ball (its BandBall is in band->m_engaged); null = none
   HitBall *ball = nullptr;
   std::vector<SurfaceContact> items;
   RubberContact::V3 vin, win;           // velocity (m/s) and spin (rad/s) at engagement
   double dwellS = 0.0;                  // time engaged
   double peakFn = 0.0;                  // largest normal force, N
   int steps = 0;                        // substeps taken
   double delta0 = 0.0, rewound = 0.0;   // start overlap (m) and how far back the touch was found (s)
   Vertex3Ds posStart;                   // where the ball was at engagement, VPX units
   Vertex3Ds regatherPos;                // where the surfaces near the ball were last collected, and when (Regather is a Verlet list: it rescans after 5 U of travel or 10 ms)
   uint32_t regatherMs = 0;
   bool regathered = false;
   bool inHole = false;                  // the ball was over a non-legacy kicker's hole when the surfaces were collected: the playfield has no floor there and is left out
   double endStored = -1.0;              // the energy check: what its contacts held at the end of the last cycle, J (< 0 = no cycle yet)
   bool groupLoop = false;               // ring_model: the last cycle went through the group step, so the touches' state is in SurfaceContact::gs (GroupStepDesign 14.4)
   bool seatPending = false;             // ring_model: found inside a rubber at engagement and not seated: the group counts it from zero; the old loop seats it if the island goes there (14.5)
};

class ContactSurfaceHost;
namespace PlayfieldCuts { struct CutQuery; }

// Shared with the contact system (contactsystem.cpp)
bool ContactTraceOn();
bool ContactTraceDetail();
bool InPlayfieldHole(const Vertex3Ds &pos);
bool PlayfieldHolesExist();
bool PairEligible(const HitBall *x);          // a ball the pair model can take part in
// shear state of a ball-to-ball contact, kept between cycles (a few pairs at most)
struct BallPairState
{
   const HitBall *a = nullptr, *b = nullptr;   // a < b by address
   RubberContact::RubberState st;
   bool inContact = false;
   bool used = false;
   double zero = 0.0;                    // the depth the squash counts from: a pair found inside with no work done has stored nothing (only ever down, as a surface's gs.zero)
   double endDepth = -1.0, endHeld = 0.0;   // at the last cycle's end: how far inside (zero not taken off), and what the pair held
};
extern std::vector<BallPairState> g_ballPairs;
class Primitive;
// Lets the physics engine step a band with no ball on it, and finish up the one that had, once per cycle, before any displacement.
class BandMover final : public MoverObject
{
public:
   explicit BandMover(ContactSurfaceHost *owner) : m_owner(owner) {}
   bool AddToList() const override;
   void UpdateDisplacements(const float) override {}
   void UpdateVelocities() override {}
   void PreUpdateDisplacements(const float dtime) override;

private:
   ContactSurfaceHost *m_owner;
};

class ContactSurfaceHost final : public ContactHost
{
public:
   // pts: closed centreline (first vertex NOT repeated), VPX units, world frame; radius: tube radius in VPX units
   ContactSurfaceHost(IEditable *editable, const std::vector<Vertex3Ds> &pts, float radius, const RubberObjectSettings &settings, const std::filesystem::path &iniPath, const std::string &name);
   // The host of a wall or primitive flagged as rubber: no collider of its own (the object's stock colliders detect the ball and call SurfaceCollide);
   // it holds the faces, the settings and steps the contact. kind: "wall" or "primitive".
   ContactSurfaceHost(IEditable *editable, std::shared_ptr<RubberSurface> surface, const RubberObjectSettings &settings, const std::filesystem::path &iniPath, const std::string &name, const char *kind);
   ~ContactSurfaceHost() override;
   // The height (VPX units) of the highest upward-facing face of this surface at table point (x, y); false if there is none there. Used to seat a turntable in the playfield.
   bool SurfaceHeightAt(double xU, double yU, double &zU) const override;
   // ContactSurface
   void Refresh(SurfaceContact &it, double ballR, double reachM, const RubberContact::V3 &c) const override;   // the faces the ball can reach this cycle
   double ArrivalS(const RubberContact::BallSI &bs, const SurfaceContact &it, int nearestSeg, const RubberContact::V3 &c, const RubberContact::V3 &v, double a, double hMax, double touchM, double arriveM) const override;
   RubberContact::TubeProbe ProbePatch(const RubberContact::BallSI &bs, const SurfaceContact &it, const RubberContact::V3 &c) const override;
   RubberContact::TubeProbe Probe(const RubberContact::BallSI &bs, const SurfaceContact &it, const RubberContact::V3 &c) const override;   // the nearest point of one surface contact (tube or faces)
   double LayerM() const override;
   // shared_finder (to-do 1.34, piece P1): builds the index over every static mesh host once. atLoad: the call at the end of the table's load (PhysicsEngine), which drops
   // any index left from a previous table first; the call at each cycle's gather is a no-op once built (it builds only if the load call did not run)
   static void FinderBuild(bool atLoad = false);
   // check mode (VPX_FINDER_CHECK=1, diagnostics only): once per engaged ball and cycle, every enabled mesh object with a face within the ball's radius must be among the ball's
   // gathered surfaces (items) and found by the index's object tree; inPlayfieldHole: the ball is over a kicker's hole (the gather leaves the playfield out on purpose)
   static void FinderCheckGather(const HitBall *ball, const std::vector<SurfaceContact> &items, const RubberContact::V3 &c, double ballR, bool inPlayfieldHole);
   // shared_finder piece P2a (runs/2026-11-11-finder/design.md section 3): the per-cycle query of one engaged ball at the cycle start (centre c, radius R + reach, v the
   // ball's velocity at the cycle start, dt the cycle), kept for the cycle; the ball's surface items point at it (SurfaceContact::finderQ). In P2a nothing reads its result
   // except check mode and the timing counters, so no result changes.
   static void FinderQuery(const HitBall *ball, const RubberContact::V3 &c, const RubberContact::V3 &v, double ballR, double dt);
   // shared finder P2b: the query finds every static mesh contact. FinderQueryAt: a query of a ball around c, radius ballR + extM, with no
   // cycle bookkeeping (a free ball's engagement test, a new contact). FinderBind: the ball's items point at its query, and each static mesh item reads that object's faces from it.
   // FinderHosts: the static mesh hosts the ball's query found, in load order.
   static void FinderQueryAt(const HitBall *ball, const RubberContact::V3 &c, double ballR, double extM);
   static void FinderBind(const HitBall *ball, std::vector<SurfaceContact> &items, bool cycleItems);   // cycleItems: the engaged ball's own set (kept for the re-query's count)
   static void FinderBindOne(const HitBall *ball, SurfaceContact &it);
   static void FinderHosts(const HitBall *ball, std::vector<ContactHost *> &out);
   // shared_finder re-engage (runs/2026-11-13-finder-reengage/design.md): an object the inside check's re-query found that is not among the ball's surfaces is recorded for engagement;
   // the group step takes it (FinderTakeNew: the ball's recorded hosts, cleared) and joins it mid-step. FinderPendAny: some ball has one recorded. One not taken by the ball's next query is
   // a fault, logged and counted ("not handled"). FinderJoinNote: counters of the joins (kind: the host's surface kind; ns: the join's own time, 0 when not timed)
   static bool FinderPendAny();
   static void FinderTakeNew(const HitBall *ball, std::vector<ContactHost *> &out);
   static void FinderJoinNote(const std::string &kind, bool againAttempt, bool insideAtStart, double ns);
   static void FinderJoinNotHandled(const HitBall *ball, const std::string &host, const char *why);
   static void FinderJoinCapHit(long n);
   static void FinderJoinTouched(int items, int touched);   // the contacts joined in an island's cycle, and how many of them were touching at the end of some step of it
   // shared_finder P4b: the cut pieces around c0 within r (the cuts' box tree, built once per cut Setup), in index order, into cq; the item's ball's cut query (never null: an item with no query gets an empty one, queried around the ball)
   static void FinderCuts(PlayfieldCuts::CutQuery &cq, const RubberContact::V3 &c0, double r);
   static PlayfieldCuts::CutQuery *FinderCutsOf(const SurfaceContact &it);
   int FinderObjIndex() const override { return m_surface != nullptr ? m_surface->fobj : -1; }
   // the finder's totals (query cost, re-queries), once a second of table time and when the table closes
   static void FinderTotals(bool force);
   // shared_finder P3a: the ball's measurement table of this object at c (null: the item does not read the ball's query); the cost log's counters since its last row
   RubberContact::FaceTable *FinderTable(const SurfaceContact &it, const RubberContact::V3 &c) const;
   static void FinderTabTake(uint64_t &checks, uint64_t &fills, uint64_t &reads);
   // Non-legacy kickers are round holes in the playfield: the playfield host leaves no floor inside them (the ball rests on the rim, or falls in and is captured).
   // Filled by Kicker::PhysicSetup, cleared when the physics engine is built. x, y, r in VPX units.
   static void AddPlayfieldHole(float x, float y, float r);
   static void ClearPlayfieldHoles();
   double m_tableMu = 0.3;   // a playfield host: the table's friction (the default mu of its material)
   static bool ContactPairModel(const HitBall *a, const BallS &b);   // both balls are contact-held and [Contact] ball_contact is on: the force model handles the pair, not the stock collision

   float HitTest(const BallS &ball, const float dtime, CollisionEvent &coll) const override;
   BandInfo GetBandInfo() const { return m_info; }
   SlingInfo GetSlingInfo() const { return m_slingInfo; }
   int PtCount() const { return (int)m_ptsSI.size(); }
   bool IsBand() const { return m_set.band; }   // found from its shape at load (DetectBand)
   bool FeatureOff() const;                     // [Features] rubbers / bands is off for what this rubber turned out to be: stock collisions (a development switch)
   float BandSectionU() const;                  // band_section_mm in VPX units for the drawn band, 0 = the rubber's Thickness
   int BandSwitchCount() const { return (int)m_switches.size(); }
   bool BandSwitchClosedAt(int i) const { return i >= 0 && i < (int)m_switches.size() && m_switches[i].closed; }
   using StopResolver = std::function<bool(const std::string &name, bool hasRadius, double radiusIn, double &x, double &y, double &r)>;
   void SetStopResolver(StopResolver f) { m_stopResolver = std::move(f); }   // band_stops entries that name a table object are looked up by the rubber's owner
   void BandSwitchesAttach(const std::vector<Primitive *> &blades);   // switch blades behind the band that only the script reads (collision off on them)
   void SetSlingEnabled(bool on) { if (m_slingPlate != nullptr) m_slingEnabled = on; }
   bool SlingTriggerPending() const { return m_slingTriggerReq; }
   void SlingTrigger() { if (m_slingPlate != nullptr) m_slingTriggerReq = true; }   // script: fire one coil pulse (taken at the next physics cycle)
   // For the renderer: whether the band is deflected or ringing, and where the rendered ring centred at c (table units) goes: the shift of its centre
   // and the rotation (row-major 3x3) of the ring about it. False = this ring does not move.
   bool BandAnimating() const { return BandFixedActive(); }
   bool BandFixedActive() const;   // its cord has something to draw
   void BandFixedDrawPrepare() const;
   bool BandDeformAt(const Vertex3Ds &c, Vertex3Ds &shift, float rot[9]) const;
   int GetType() const override { return e3DLine; }
   void Collide(const CollisionEvent &coll) override;
   void Contact(CollisionEvent &coll, const float dtime) override;
   MoverObject *GetMoverObject() override { return &m_mover; }
   // Sling: the table's primitives for the plate and the switch blades (found by name by the Rubber). Turns their collision off and places the plate.
   void SlingAttach(Primitive *plate, const std::vector<Primitive *> &switches);
   void CalcHitBBox() override;
   void DrawUI(std::function<Vertex2D(Vertex3Ds)> project, ImDrawList *drawList, bool fill) const override;

   void SetEventCarrier(HitObject *carrier) { m_eventCarrier = carrier; }   // a wall / primitive / ramp host: its first stock collider carries the part's hit event settings
   HitObject *HitEventObject() const override;

private:
   void Reload() const override;
   void DetectBand();
   void SlingFinish();                     // the table's primitives follow the plate and the switches
   struct BandSwitch;
   void SwitchPlace(BandSwitch &sw) const;   // find the span a blade senses and the side it is on
   void SwitchesFinish();                   // the blades follow the band
   void SlingPlace() const;   // put the plate so its face is slingGapMm from the band (from its table pose, so it can be repeated when the gap changes)
   friend class BandMover;
   friend void HitRubberBandBallDeleted(HitBall *ball);
   friend void HitRubberContactRelease(HitBall *ball, const char *why);
   friend struct ContactSystemImpl;   // the contact system steps the contact-held balls (contactsystem.cpp)
   // A ball on a band: the contact system steps the cord with its balls (CordIsland), with the band's state kept in the cord between cycles.
   bool BandWants(const HitBall *ball) const;
   bool BandEngage(HitBall *ball);
   void BandRelease(HitBall *ball, const char *why);
   void BandReleaseAll(const char *why);
   void BandAdvance(const float dtime);
   // The cord in the group step (cordband.h; docs/unified-physics/GroupStepDesign.txt 8, 11): the contact system steps it with its balls
   ContactSurfaceHost *BandHost() override { return this; }
   bool CordOn() const { return m_cord != nullptr; }
   bool CordModel() const { return m_set.band && m_bandPosts.size() >= 2; }   // this band runs on the cord: the contact system takes its balls (CordNear), no stock collision
   CordBand *Cord() const { return m_cord.get(); }
   void CordEnsure();                  // the cord from the band's posts, plate and switches
   double CordRing(double loopM, double rcM, double EA) const;   // the ring's inner diameter, inches: band_ring_in, or the standard size for this loop
   void CordBegin();                   // once per cycle, before the step: the script's trigger, the switches, the driver
   void CordFinish(float dtime);       // once per cycle, after: the script's events, the table's primitives
   void CordSolveFailed(const char *why);
   bool m_cordStepped = false;         // the contact system stepped the cord with its balls this cycle
   std::shared_ptr<CordBand> m_cord;
   int m_cordFired = 0;

   // resolver for band_stops entries "Post7" or "Post7, r" (hasRadius / radiusIn: an explicit radius was given; false = not found or no radius of its own;
   // x, y, r in VPX units)
   StopResolver m_stopResolver;
   std::set<std::string> m_stopWarned;   // names already reported as not found (once each)

   // check mode (VPX_FINDER_CHECK=1): at one evaluation of this mesh, every face within the ball's radius (a full scan) against today's list (it.Faces(c)) and the index
   void FinderCheck(const SurfaceContact &it, const RubberContact::V3 &c, double ballR, int evalKind) const;
   // shared_finder P2a (design.md section 4): the inside check of one evaluation at ball centre c: outside the ball's query volume, the query runs again around its centre
   // with radius R + 2 x the stray (counted; nothing reads the result in P2a)
   void FinderInside(const SurfaceContact &it, const RubberContact::V3 &c) const;
   static void FinderRun(void *q, double r);   // the query itself (q: a FinderBallQ), radius r around its centre
   std::shared_ptr<RubberSurface> m_surface;   // a wall / primitive host: its faces (m_pts is then empty)
   int m_linkOver = 0;                         // multi_contact: edges on more than two faces, counted when the faces' shared edges are built at load (Link)
   mutable bool m_linkNoticed = false;         // that count reported (log and screen) at the first ball near the surface, once
   HitObject *m_eventCarrier = nullptr;        // see SetEventCarrier
   std::vector<Vertex3Ds> m_pts; // closed loop, last = first
   RubberContact::TubePre m_tubePre;       // what ProbeTube precomputes from m_ptsSI
   // shared_finder P5: each segment of the (static) centreline is a piece of this tube host, its box the segment's end points grown by the tube radius; one box tree over them
   // (4 per leaf), built once at load. A ball's query finds the segments within its volume (SegsFor); the tube readers read that list in place of the whole tube. Cord bands keep the whole tube (P7).
   std::vector<RubberSurface::FinderNode> m_segNodes;
   std::vector<int> m_segOrder;
   bool SegsOn() const;   // the finder is built and not a cord band
   void SegQuery(void *t, const RubberContact::V3 &c0, double cov) const;   // t: a TubeSegQ; every segment whose centreline lies within cov of c0, in index order
   const std::vector<int> &SegsFor(const SurfaceContact &it, const RubberContact::V3 &c, double need, double *covAt = nullptr) const;   // the ball's segments, every one within need of c among them (re-queried if not)
   RubberContact::TubeProbe ProbeSegs(const RubberContact::BallSI &bs, const SurfaceContact &it, const RubberContact::V3 &c, double need, const std::vector<int> *&segs) const;   // ProbeTube's exact answer from the list
   std::vector<RubberContact::V3> m_ptsSI; // the same in metres (pivot-free world frame), for the solve
   double m_radiusSI = 0.0;                // tube radius, m
   float m_radius;               // tube radius, VPX units
   std::vector<BandPost> m_bandPosts;
   std::vector<BandSpan> m_bandSpans;
   BandMover m_mover{ this };
   BandInfo m_info;
   std::vector<BandBall> m_engaged;       // the balls on the band, one per span
   mutable Primitive *m_slingPlate = nullptr;
   mutable std::vector<Primitive *> m_slingSwitches;
   mutable Vertex3Ds m_slingBase;          // the plate's position as the table has it
   mutable double m_slingGapTableMm = 0.0; // the face-to-band gap at that position, mm
   mutable double m_slingPlacedGapMm = -1.0;
   // the sling's plate body, driver and switches (see hitsling.h)
   struct BandSwitch
   {
      Primitive *prim = nullptr;
      float baseTransX = 0.f;      // the table's TransX: the blade's resting value; the engine slides it with the band
      int span = -1;               // the span of the band it senses (index into m_bandSpans)
      double f = 0.0;              // where it senses the band: fraction of the span chord
      double bx = 0.0, by = 0.0;   // unit vector from the band toward the blade (behind the band): the band pushed this way closes the switch
      double axisProj = 0.0;       // the blade's own X axis dotted with that vector: how far TransX must move for a given push
      bool coil = false;           // a sling switch: wired to the coil driver (the others are read by the script only)
      bool legacySign = false;     // a sling switch: TransX decreases with the push (the table authors orient the blades that way)
      bool closed = false;
      double dispM = 0.0;          // band displacement toward the blade at the switch, m
      float shownTransX = 0.f;
   };
   mutable std::vector<BandSwitch> m_switches;
   mutable int m_slingSpan = -1;                 // the span the plate faces
   mutable double m_slingNx = 0.0, m_slingNy = 0.0;   // unit normal of that span, from the plate toward the band
   SlingPlate m_plate;
   SlingDriver m_driver;
   bool m_plateConfigured = false;
   bool m_slingEnabled = true;
   bool m_slingTriggerReq = false;   // SlingTrigger() was called; consumed by SlingAdvance
   bool m_slingAnyClosed = false;
   SlingInfo m_slingInfo;
   // the plate's geometry (SI), measured when the plate is placed: the face at band height is at originN + a / cos(phi) + h tan(phi) from the rest centreline
   mutable bool m_plateGeomValid = false;
   mutable double m_plateOriginN = 0.0, m_plateA = 0.0, m_plateH = 0.0;   // pivot axis place, face offset from the axis, band height above the axis
   mutable double m_plateRest = 0.0;                                       // the plate's angle at rest from the touching pose (RotX = 0), rad, <= 0: turned back for the gap
   mutable double m_plateF0 = 0.5, m_plateF1 = 0.5;                        // the face's extent along the span, fractions of the chord
   mutable float m_slingShownRotX = 0.f;
};

// A ball being deleted must not stay engaged with a band
void HitRubberBandBallDeleted(HitBall *ball);
// A force-model contact (not a band) holds this ball / releases it: the flipper takes over a ball that touches its bat.
bool HitRubberContactHolds(const HitBall *ball);
void HitRubberContactRelease(HitBall *ball, const char *why);
