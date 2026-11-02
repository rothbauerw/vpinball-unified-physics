// license:GPLv3+

#pragma once

#include <filesystem>
#include <chrono>
#include <string>
#include <vector>
#include <functional>
#include <set>
#include <memory>

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

// The faces of a wall or primitive flagged as rubber, in the world frame (metres). Owned by the HitRubber that hosts the object's contact.
struct RubberSurface
{
   std::vector<RubberContact::SurfaceTri> tris;
   float lo[3] = { 0.f, 0.f, 0.f }, hi[3] = { 0.f, 0.f, 0.f };   // bounding box, VPX units
   void Add(const Vertex3Ds &a, const Vertex3Ds &b, const Vertex3Ds &c);   // one face (VPX units)
};

// One ball in an engine-stepped contact with a static rubber (post ring, wall rubber): the rubbers it touches and the bookkeeping for the log and the release.
class FlipperMoverObject;

struct BallContact
{
   class ContactHost *host = nullptr;   // the surface or part that engaged the ball (its settings drive the step)
   class HitRubber *band = nullptr;       // a band that holds this ball (its BandBall is in band->m_engaged); null = none
   HitBall *ball = nullptr;
   std::vector<SurfaceContact> items;
   double omega = 0.0;                   // natural frequency of the contact in the last cycle, rad/s (sets the largest stable step); 0 = unknown / no contact
   double lastSpeed = 1e9;               // speed at the end of the last cycle, m/s
   double lastVn = 1e9;                  // fastest approach speed along a contact normal at the end of the last cycle, m/s (1e9 = not in contact)
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

class HitRubber;

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
   double maxDelta = 0.0, maxFn = 0.0;   // this cycle (trace)
   double lastDelta = 0.0;               // the compression of the last cycle: the stiffness the next cycle's step is chosen on
   double lastSpeed = 0.0;               // the normal speed at the last cycle's start, either way: the step rule's "slow"
   double zero = 0.0;                    // the depth the squash counts from: a pair found inside with no work done has stored nothing (only ever down, as a surface's gs.zero)
   double endDepth = -1.0, endHeld = 0.0;   // at the last cycle's end: how far inside (zero not taken off), and what the pair held
};
extern std::vector<BallPairState> g_ballPairs;
class Primitive;
// Lets the physics engine step a band with no ball on it, and finish up the one that had, once per cycle, before any displacement.
class BandMover final : public MoverObject
{
public:
   explicit BandMover(HitRubber *owner) : m_owner(owner) {}
   bool AddToList() const override;
   void UpdateDisplacements(const float) override {}
   void UpdateVelocities() override {}
   void PreUpdateDisplacements(const float dtime) override;

private:
   HitRubber *m_owner;
};

class HitRubber final : public ContactHost
{
public:
   // pts: closed centreline (first vertex NOT repeated), VPX units, world frame; radius: tube radius in VPX units
   HitRubber(IEditable *editable, const std::vector<Vertex3Ds> &pts, float radius, const RubberObjectSettings &settings, const std::filesystem::path &iniPath, const std::string &name);
   // The host of a wall or primitive flagged as rubber: no collider of its own (the object's stock colliders detect the ball and call SurfaceCollide);
   // it holds the faces, the settings and steps the contact. kind: "wall" or "primitive".
   HitRubber(IEditable *editable, std::shared_ptr<RubberSurface> surface, const RubberObjectSettings &settings, const std::filesystem::path &iniPath, const std::string &name, const char *kind);
   ~HitRubber() override;
   // The height (VPX units) of the highest upward-facing face of this surface at table point (x, y); false if there is none there. Used to seat a turntable in the playfield.
   bool SurfaceHeightAt(double xU, double yU, double &zU) const override;
   // ContactSurface
   void Refresh(SurfaceContact &it, double ballR, double reachM, const RubberContact::V3 &c) const override;   // the faces the ball can reach this cycle
   double ArrivalS(const RubberContact::BallSI &bs, const SurfaceContact &it, int nearestSeg, const RubberContact::V3 &c, const RubberContact::V3 &v, double a, double hMax, double touchM, double arriveM) const override;
   RubberContact::TubeProbe ProbePatch(const RubberContact::BallSI &bs, const SurfaceContact &it, const RubberContact::V3 &c) const override;
   RubberContact::TubeProbe Probe(const RubberContact::BallSI &bs, const SurfaceContact &it, const RubberContact::V3 &c) const override;   // the nearest point of one surface contact (tube or faces)
   double LayerM() const override;
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
   // The cord in the group step (cordband.h; docs/GroupStepDesign.txt 8, 11): the contact system steps it with its balls
   HitRubber *BandHost() override { return this; }
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

   std::shared_ptr<RubberSurface> m_surface;   // a wall / primitive host: its faces (m_pts is then empty)
   HitObject *m_eventCarrier = nullptr;        // see SetEventCarrier
   std::vector<Vertex3Ds> m_pts; // closed loop, last = first
   RubberContact::TubePre m_tubePre;       // what ProbeTube precomputes from m_ptsSI
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
