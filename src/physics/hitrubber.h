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
#include "physics/hitsling.h"
namespace HingePart { class Part; }
class ContactMovingPart;
class Gate;
class Spinner;
class Bumper;
namespace JetBumper { class EnginePart; struct Params; }

// Static rubber (post ring, wall rubber) as ONE collider: a closed tube of circular cross-section swept along the rubber's centreline.
// Replaces the triangle/edge/point set of Rubber::PhysicSetup when rubber_contact is on for that rubber (vpx_physics.ini).
// The ball is detected by an ordinary swept test (first touch, exactly as the other colliders). Collide() then solves the whole
// contact in one call: it integrates the force model of physics/rubbercontact.h for the ball alone (SI, small local steps) until the
// ball has left, and hands the exit velocity, spin and position back. Nothing else acts on the ball during the dwell (about 5 ms at
// post speeds), which is the same simplification as the 1 ms engine step.

struct RubberObjectSettings
{
   bool contact = false;          // use the force model for this rubber (else the stock triangle/edge/point colliders are built)
   // Walls and primitives flagged as rubber ([Wall:<name>] / [Primitive:<name>] rubber = 1): the same model on the object's faces. Never set by [Rubbers].
   bool surfaceRubber = false;
   double surfaceLayerMm = 5.4;   // layer thickness of a flat rubber face when layer_mm is 0: a typical ring's diameter (5.4 mm on the Gorgar rings)
   double edgeRadiusMm = 1.0;     // convex edges and corners of a rubber face are treated as rounded to this radius (Hertz effective radius)
   double eMpa = 12.0;            // effective Young's modulus
   double layerMm = 0.0;          // rubber layer thickness for the thin-layer stiffening; 0 = the tube diameter
   double alphaSPerM = 0.3;       // normal damping; rebound ratio ~ 1 - alpha*v (0.3 fits the Houdini post data; unverified)
   double betaLayer = 0.33;       // thin-layer stiffening coefficient
   double plateArealKgM2 = 0.0;   // a gate or spinner plate: its areal density (kg/m^2), or its mass (g) when plate_mass_g > 0; 0 = the default of the kind (UNVERIFIED)
   double plateMassG = 0.0;
   double plateClearanceU = 8.0;  // the plate's axis is at least this far above the top of the ball (50 U): the plate swings clear of the ball instead of resting on its top
   bool viscoelastic = false;     // viscoelastic (Kuwabara-Kono) normal damping calibrated to `restitution` at `restitutionRefMs` instead of alpha_s_per_m: the rebound ratio falls only slowly with speed (hard surfaces)
   double restitution = 0.0;      // the rebound ratio at the reference speed; 0 = from the part's elasticity (with its falloff at that speed)
   double restitutionRefMs = 1.0; // the reference speed, m/s
   double gammaKK = 0.0;          // the calibrated damping coefficient (not an ini key)
   double mu = 0.9;               // rubber / steel friction
   double shearScale = 1.0;       // multiplier on the shear stiffness
   double shearZeta = 0.2;        // shear damping ratio
   bool twist = true;             // drilling friction
   double tanDelta = 0.05;        // loss tangent for the rolling resistance, 0 = off
   double rollingMu = 0.0;        // floor of the rolling-resistance coefficient (the playfield: [Contact] rolling_resistance, because the hysteresis law below gives almost nothing on hard wood); 0 = none
   bool normalLoss = true;        // the same tan_delta also damps the normal force of a post / surface rubber (a dashpot at the contact's own natural frequency)
   double substepUs = 5.0;        // local integration step inside Collide()
   double restOmegaH = 0.2;       // a resting ball steps at most restOmegaH / omega (omega = the contact's natural frequency); symplectic Euler is stable below 2, 0.2 is accurate for an impact
   double restSubstepUs = 50.0;   // largest step for a ball at rest on the rubber (the stiffness sets the real limit); 0 = always substep_us
   double restSpeedMmS = 30.0;    // a ball slower than this (start and end of the cycle) counts as at rest
   double rollTravelUm = 50.0;    // a ball rolling or sliding along a surface (normal speed below rest_speed_mm_s) steps at most this far per substep, micrometres; 0 = always substep_us
   double maxMs = 30.0;           // give up after this much simulated contact time
   bool slowContact = true;       // a ball resting on or rolling along a post / surface rubber (slower than a hit) is carried by the force model too, not the stock contact
   double releaseGapMm = 0.2;     // an engaged ball is released only once it is farther than this from every rubber (and receding): a resting ball must not engage and release over and over

   // Band (a rubber stretched between posts): the posts and free spans are found at load (logged); a ball hitting a free span deflects it
   // (a tensioned cord, see BandChain in rubbercontact.h). Other parts of the rubber (the post wraps, short spans) stay rigid.
   bool band = false;             // this rubber is a band
   double bandPrestretchPct = 7.0; // stretch of a span at rest, percent
   double bandTensionN = 0.0;      // rest tension of a span, N; 0 = from the rubber's modulus and the stretch (the cord's EA is then e_mpa x area), > 0 = this tension (EA = tension / stretch)
   double bandRefHz = 20.0;       // frequency at which tan_delta (the rubber's own loss tangent, shared with the rings) applies to the stretch
   double bandNodeMassG = 0.0;    // moving mass of the deflected span, g; 0 = half the span's mass
   double bandNodeZeta = 0.05;    // damping ratio of the deflected span's own vibration (internal friction, air): how fast the ringing dies out
   double bandShearZeta = 1.0;    // damping ratio of the sideways (shear) spring of the ball / band contact; < 0 = the same as shear_zeta. 1 = critical: a ball resting on the band does not rock
   double bandSubstepUs = 20.0;   // step of the band's integration, microseconds; 0 = the same as substep_us. Larger = faster, a little less exact (40 is too coarse)
   bool bandFollow = true;       // a ball rolling along the span loads a new node under itself once it is more than band_follow_mm from its node, and idle nodes that have relaxed onto the line are merged away; false = nodes stay where the ball engaged
   double bandNodeMm = 0.0;       // > 0: a plain band (no sling) gets permanent nodes this far apart along every span, a ball rolls on the polyline they form and loads the legs under it; 0 = nodes are made under the ball (band_follow)
   double bandFollowMm = 6.0;     // how far a ball may be from the node it loads along the span before it takes a new one, mm (nodes closer than about 5 mm are too stiff for the 20 us step: 4 mm already trips)
   bool bandRinging = true;       // false = no ringing: after the last ball has left the band returns without oscillating (critical damping)
   double bandRingZeta = -1.0;    // damping ratio of the ringing after the last ball has left; < 0 = the same as band_node_zeta
   double bandRestMm = 0.05;      // the band counts as at rest (stops being simulated and drawn deflected) below this deflection, mm
   bool bandLockZ = true;         // the deflected span keeps its height (the posts hold the band); false = it can ride up over the ball
   double bandMaxMs = 150.0;      // not used (a band contact lasts as long as the ball is on it); kept so old files load
   double bandFlipperMarginU = 1.0; // a ball engages the band only if no flipper is closer than this to its surface, VPX units
   double bandCurveRadiusMm = 40.0; // a centreline vertex whose local bend radius is below this belongs to a post wrap
   double bandMinSpanMm = 20.0;   // a straight run shorter than this between two posts is not a free span
   // Sling (a band with a kicker plate behind it): the plate and the switch blades are table primitives named here. The engine places the plate at load
   // (and again when the gap changes) so its face is sling_gap_mm from the band, and turns collision off on all of them. The table file is never changed.
   std::string slingPlate;        // name of the plate primitive; its position is the pivot. Empty = this band is not a sling
   std::string slingSwitches;     // names of the switch blade primitives, comma separated (wired in parallel)
   std::string switches;          // names of switch blade primitives behind the band that only the script reads (BandSwitchClosed, SwitchClosed / SwitchOpened), comma separated
   double switchGapMm = 1.5875;   // band displacement toward such a blade that closes its switch, mm (1/16 in)
   double switchMinClosedMs = 2.0; // such a switch counts as closed for at least this long, ms
   double slingGapMm = 4.7625;    // distance from the plate face to the back of the band at rest, mm (3/16 in)
   SlingSettings sling;           // the plate, its coil and its driver (sling_<flipper key>, sling_pulse_ms, sling_min_closed_ms, sling_switch_gap_mm)
   std::string bandStops;         // obstacles the band cannot enter: "x,y,r | x,y,r | ..." in VPX units (centre, radius of the post); empty = none
   std::string bandPosts;         // optional override: "x,y,r | x,y,r | ..." in VPX units (centre, centreline radius of the wrap); empty = fit from the centreline

   RubberContact::RubberParams ToParams() const;
};

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

// The band's state while a span is deflected: the deformable centreline (the ring with the span replaced by A - P - B) and the node P.
// Lives in the HitRubber, not in a contact solve, so the band keeps moving (and can be read by scripts) after the ball has left.
struct BandDyn
{
   bool active = false;                  // the band is being loaded or is ringing
   std::vector<RubberContact::V3> pts;   // closed ring with the span of every chain that has nodes replaced by A - P_1 - ... - P_M - B
   std::vector<RubberContact::BandChain> chains;   // one per free span
   std::vector<int> span;                // per chain: index into the span list
   std::vector<int> base;                // per chain with nodes: index of its A in pts (its nodes follow, then B)
   RubberContact::BandParams bp;         // what the chains were built with
   double peakFn = 0.0;                  // peak ball contact force, N
   double peakBB = 0.0;                  // deepest overlap of two spans, m
};

// One ball on the band: the node of the span it loads, its shear state, and the bookkeeping for the log and the release.
struct BandBall
{
   HitBall *ball = nullptr;
   int chain = 0;                        // the chain (span) this ball loads, index into BandDyn::chains
   int node = 0;                         // the node of that chain this ball loads (a node id)
   RubberContact::RubberState st;        // shear state of this ball's contact with the band
   double clearS = 0.0;                  // how long the ball has been clear of the band and receding
   double dwellS = 0.0;
   double peak = 0.0;                    // peak deflection of its node while it was on the band, m
   RubberContact::V3 vin;                // its velocity at engagement
   double side = 0.0;                    // +1 / -1: which side of the span's rest line the ball came from (the sign of (centre - A) . (ab.y, -ab.x)); 0 = not known
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
   class HitRubber *host = nullptr;   // the rubber that engaged the ball (its settings drive the step)
   FlipperMoverObject *bat = nullptr;     // a bat whose rubber is part of this contact (the ball is a flipper's as well as the surfaces'); null = none
   class HitRubber *band = nullptr;       // a band whose nodes this ball loads (its BandBall is in band->m_engaged); null = none
   RubberContact::RubberState batSt;     // the bat rubber's shear state
   bool batJustEngaged = false;          // the first step of this bat engagement seats the starting overlap
   int batClear = 0;                     // consecutive substeps clear of the bat and receding
   double batLastDist = 0.0;             // ball centre to bat surface at the previous substep, m
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
};

// The force-model rubbers around a ball that a flipper's rubber model drives (the cradle: a post and the bat both touch the ball). The flipper's step adds their
// force to its own, so the ball is integrated once; the shear state of each rubber is kept here between cycles.
struct RubberSideSet
{
   std::vector<SurfaceContact> items;
   std::vector<double> depthM;   // deepest overlap of each item since the last reset, m (log and safety net)
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
};
extern std::vector<BallPairState> g_ballPairs;
class Primitive;
// Lets the physics engine step a band that has a ball on it (and the band alone while it rings down), once per cycle, before any displacement.
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

class HitRubber final : public HitObject, public ContactSurface
{
public:
   // pts: closed centreline (first vertex NOT repeated), VPX units, world frame; radius: tube radius in VPX units
   HitRubber(IEditable *editable, const std::vector<Vertex3Ds> &pts, float radius, const RubberObjectSettings &settings, const std::filesystem::path &iniPath, const std::string &name);
   // The host of a wall or primitive flagged as rubber: no collider of its own (the object's stock colliders detect the ball and call SurfaceCollide);
   // it holds the faces, the settings and steps the contact. kind: "wall" or "primitive".
   HitRubber(IEditable *editable, std::shared_ptr<RubberSurface> surface, const RubberObjectSettings &settings, const std::filesystem::path &iniPath, const std::string &name, const char *kind);
   ~HitRubber() override;
   // A stock collider of the flagged object found the ball. True if the force model holds it from now on (the stock impulse is skipped).
   bool SurfaceCollide(HitBall *ball, const Vertex3Ds &n);
   // A resting or slow touch (the engine's contact pass) of a rubber post or a flagged object. True if the force model holds the ball from now on (the stock contact is skipped).
   bool SlowContact(HitBall *ball);
   // The playfield as a surface on the force model. s gets the material of the [Playfield] section (defaults: a steel ball on lacquered wood, all UNVERIFIED), tableMu the
   // table's friction. False if [Features] playfield_contact is off.
   static bool LoadPlayfieldSettings(RubberObjectSettings &s, const std::filesystem::path &file, double tableMu);
   // A turntable's disc: a hard textured surface (the playfield's material with the given friction). Never off.
   static RubberObjectSettings TurntableSettings(double mu, double rollingResistance);
   // The height (VPX units) of the highest upward-facing face of this surface at table point (x, y); false if there is none there. Used to seat a turntable in the playfield.
   bool SurfaceHeightAt(double xU, double yU, double &zU) const;
   int m_turntable = -1;     // >= 0: this host is the disc of that turntable (physics/turntable.h): its surface moves, and the balls' reaction turns the disc
   static std::vector<HitRubber *> &AllHosts() { return Registry(); }
   bool IsPlayfieldKind() const { return m_surfaceKind == "playfield"; }
   // ContactSurface
   bool IsPlayfield() const override { return IsPlayfieldKind(); }
   double ReleaseGapM() const override { return m_set.releaseGapMm * 1e-3; }
   void Refresh(SurfaceContact &it, double ballR, double reachM, const RubberContact::V3 &c) const override;   // the faces the ball can reach this cycle
   RubberContact::TubeProbe Probe(const RubberContact::BallSI &bs, const SurfaceContact &it, const RubberContact::V3 &c) const override;   // the nearest point of one surface contact (tube or faces)
   void VelocityAt(const RubberContact::TubeProbe &pr, RubberContact::V3 &v, RubberContact::V3 &w) const override;   // a turntable disc turns about its axis
   void ApplyReaction(const RubberContact::TubeProbe &pr, const RubberContact::ContactResult &r, double ballR, double h) const override;   // the balls turn the disc
   const RubberObjectSettings *Settings() const override { return m_part != nullptr ? nullptr : &m_set; }   // a gate or spinner sets its own step (Part::BallStep), not substep_us
   bool IsFloor() const override { return m_surfaceKind == "playfield" || m_surfaceKind == "turntable"; }   // a surface a ball rests on: gravity acts in the contact substeps while one is in the contact
   // Non-legacy kickers are round holes in the playfield: the playfield host leaves no floor inside them (the ball rests on the rim, or falls in and is captured).
   // Filled by Kicker::PhysicSetup, cleared when the physics engine is built. x, y, r in VPX units.
   static void AddPlayfieldHole(float x, float y, float r);
   static void ClearPlayfieldHoles();
   // Surface hosts that are still registered: called when a new physics engine starts, where there must be none (a host the previous table did not release is tested against every ball of this one).
   // Empty if there are none, else "N (kind name, ...)".
   static std::string LeftoverHosts();
   double m_tableMu = 0.3;   // a playfield host: the table's friction (the default mu of its material)
   static bool ContactPairModel(const HitBall *a, const BallS &b);   // both balls are contact-held and [Contact] ball_contact is on: the force model handles the pair, not the stock collision
   bool ContactOwns(const BallS &ball) const;   // this rubber is one of the rubbers of a force-model contact that holds the ball
   // A ball driven by a flipper's rubber model feels the rubbers around it too (RubberSideSet). SideGather keeps the state of the rubbers that stay near, adds new
   // ones and drops the rest (c, ballR: ball centre and radius in the world frame, m; reachM: how far the ball can travel this cycle). SideStep is one substep: it adds
   // the force (N) and the torque on the ball (N m) of every item to F and Mb and records the overlap in set.depthM.
   static void SideGather(const BallS &ball, const RubberContact::V3 &c, double ballR, double reachM, RubberSideSet &set);
   static void SideStep(RubberSideSet &set, const RubberContact::BallSI &bs, const RubberContact::V3 &c, const RubberContact::V3 &v, const RubberContact::V3 &w, double h,
                        RubberContact::V3 &F, RubberContact::V3 &Mb);
   const std::string &Name() const override { return m_name; }
   bool SideHolds(const HitBall *ball) const;   // a flipper's rubber model drives this ball and has this rubber in its set: the stock contact stays out

   float HitTest(const BallS &ball, const float dtime, CollisionEvent &coll) const override;
   BandInfo GetBandInfo() const { return m_info; }
   SlingInfo GetSlingInfo() const { return m_slingInfo; }
   int PtCount() const { return (int)m_ptsSI.size(); }
   int BandSwitchCount() const { return (int)m_switches.size(); }
   bool BandSwitchClosedAt(int i) const { return i >= 0 && i < (int)m_switches.size() && m_switches[i].closed; }
   using StopResolver = std::function<bool(const std::string &name, bool hasRadius, double radiusIn, double &x, double &y, double &r)>;
   void SetStopResolver(StopResolver f) { m_stopResolver = std::move(f); }   // band_stops entries that name a table object are looked up by the rubber's owner
   void BandSwitchesAttach(const std::vector<Primitive *> &blades);   // switch blades behind the band that only the script reads (collision off on them)
   void SetSlingEnabled(bool on) { if (m_slingPlate != nullptr) m_slingEnabled = on; }
   void SlingTrigger() { if (m_slingPlate != nullptr) m_slingTriggerReq = true; }   // script: fire one coil pulse (taken at the next physics cycle)
   // For the renderer: whether the band is deflected or ringing, and where the rendered ring centred at c (table units) goes: the shift of its centre
   // and the rotation (row-major 3x3) of the ring about it. False = this ring does not move.
   bool BandAnimating() const { return m_dyn.active; }
   bool BandDeformAt(const Vertex3Ds &c, Vertex3Ds &shift, float rot[9]) const;
   int GetType() const override { return e3DLine; }
   void Collide(const CollisionEvent &coll) override;
   void Contact(CollisionEvent &coll, const float dtime) override;
   MoverObject *GetMoverObject() override { return &m_mover; }
   // Sling: the table's primitives for the plate and the switch blades (found by name by the Rubber). Turns their collision off and places the plate.
   void SlingAttach(Primitive *plate, const std::vector<Primitive *> &switches);
   void CalcHitBBox() override;
   void DrawUI(std::function<Vertex2D(Vertex3Ds)> project, ImDrawList *drawList, bool fill) const override;

   // Reads [Rubbers], [Slings], [Rubber:<name>], [Sling:<name>] of vpx_physics.ini over the built-in defaults (no file: the defaults). False on the stock engine.
   static bool LoadSettings(RubberObjectSettings &s, const std::filesystem::path &file, const std::string &name);
   // A wall or primitive: [Rubbers] as the defaults, then [<kind>:<name>] (kind = "wall" or "primitive"). True if the object is flagged rubber = 1 and
   // [Features] surface_rubbers allows it; s.contact is then set. Only a unified vpx_physics.ini can flag surfaces.
   static bool LoadSurfaceSettings(RubberObjectSettings &s, const std::filesystem::path &file, const char *kind, const std::string &name);
   // Unified engine: any other wall or primitive is a hard surface of the contact model. HardSurfaceWanted: the feature is on and the object does not opt out ([Wall:<name>] rubber = 0).
   // LoadHardSurfaceSettings: the playfield's material (steel ball on lacquered wood) with the part's own friction mu and its restitution e turned into the normal damping
   // (rebound ratio ~ 1 - alpha v at 1 m/s), then [Surfaces] and the object's own section on top. PROVISIONAL, unverified: the mapping from the stock elasticity.
   // A gate or spinner of the unified engine as a host (docs/GateSpinnerDesign.txt): the part model (contactgate.h / contactspinner.h) is the surface and the body; the host engages the
   // balls near it, its mover advances the part when no ball is near, writes its angle back to the part's render state (angle, angleSpeed: the stock mover's fields, rad and rad per VPT)
   // and fires its events on the table part (gate or spinner).
   static HitRubber *MakePartHost(IEditable *part, const char *kind, const std::string &name, const std::filesystem::path &ini, std::shared_ptr<HingePart::Part> model, Gate *gate, Spinner *spinner,
                                  float *angle, float *angleSpeed, double elasticity, double friction);
   static bool HingeWanted(const std::filesystem::path &file, const char *kind, const std::string &name);
   // A jet bumper of the unified engine as a host (docs/JetBumperDesign.txt): one host with the bumper's three surfaces; its mover advances the model when no ball is near, writes the
   // pose back to the table part (ring drop, skirt tilt, switch, coil current) and fires Hit / Switch on the table part. box: x, y, z of the lowest corner and the highest, VPX units.
   static HitRubber *MakeBumperHost(IEditable *part, const std::string &name, const std::filesystem::path &ini, std::shared_ptr<ContactMovingPart> model, Bumper *bumper, const float lo[3], const float hi[3]);
   static bool BumperWanted(const std::filesystem::path &file, const std::string &name);   // [Features] bumpers on, VPX_BUMPERS not 0, and the bumper's own section does not say physical = 0   // like HardSurfaceWanted, for [Features] hinges (gates and spinners)
   static bool HardSurfaceWanted(const std::filesystem::path &file, const char *kind, const std::string &name);
   static bool LoadHardSurfaceSettings(RubberObjectSettings &s, const std::filesystem::path &file, const char *kind, const std::string &name, double elasticity, double friction);
   void SetEventCarrier(HitObject *carrier) { m_eventCarrier = carrier; }   // a wall / primitive / ramp host: its first stock collider carries the part's hit event settings
   HitObject *HitEventObject() const override;
   ContactMovingPart *MovingPart() const override { return m_part.get(); }
   double SurfaceSpeedBoundMS(double cycleS) const override;
   RubberContact::ContactResult EvaluateOwn(SurfaceContact &it, const RubberContact::BallSI &bs, const RubberContact::V3 &c, const RubberContact::V3 &v, const RubberContact::V3 &w, double h) const override;
   bool HitEventWanted(const RubberContact::V3 &normalToBall) const override;
   void PartAdvance(float dtime);           // once per physics interval (the mover): advance the part alone unless an island stepped it, write the angle back, fire its events
   static RubberContact::V3 PartGravitySI(); // the gravity a part feels (m/s^2, table frame): gravity minus the cabinet's nudge acceleration, as the balls feel it
   void SetHardSurface(double elasticity, double friction) { m_hardSurface = true; m_partE = elasticity; m_partMu = friction; }
   static double HardSurfaceE(const HitObject *first);   // the part's restitution at 1 m/s: its elasticity with the stock falloff (0.3 without a collider)
   static std::filesystem::path IniPathFor(const std::filesystem::path &tableFile);

private:
   // The contact solve for one ball (see the header comment). Returns false if no contact developed and the ball is untouched.
   void ReloadIfChanged() const;
   void DetectBand();
   void SlingAdvance(const float dtime);   // switches and driver, once per physics cycle; steps the plate itself unless it is working on the band
   void SlingFinish();                     // the table's primitives follow the plate and the switches
   struct BandSwitch;
   void SwitchPlace(BandSwitch &sw) const;   // find the span a blade senses and the side it is on
   bool SwitchesUpdate(double dt);          // once per cycle: read the band at every blade, close / open the switches, fire the events; true if a sling (coil) switch is closed
   void SwitchesFinish();                   // the blades follow the band
   void BandAdvanceCore(const float dtime);
   // the band's cycle in phases (see BandCycleBegin in hitrubber.cpp)
   struct BandCycle
   {
      double dt = 0.0, h = 0.0;        // the cycle and the substep, s
      int steps = 1;
      RubberContact::RubberParams rp;  // the band's rubber (band_shear_zeta applied)
      double hLayer = 0.0;
      std::vector<std::vector<RubberContact::V3>> fNode;   // contact reactions on each chain's nodes, by node id
      std::vector<std::vector<char>> loaded;               // by node id: loaded by a ball (not to be merged away or taken over)
      bool nodeTripped = false;
      bool active = false;
   };
   BandCycle m_bc;
   struct BandBallContact
   {
      bool inContact = false;
      RubberContact::V3 force, torque;   // on the ball
      RubberContact::V3 vs, normal;      // the band's velocity at the contact and the contact normal
      double delta = 0.0, gap = 0.0;
   };
   bool BandCycleBegin(float dtime);
   void BandZeroForces();
   void BandChainContacts();
   void BandMarkLoaded();
   void BandRingDown();
   void BandSubstepBegin();
   void BandBallForce(BandBall &e, const RubberContact::BallSI &bs, const RubberContact::V3 &c, const RubberContact::V3 &v, const RubberContact::V3 &w, BandBallContact &out);
   bool BandSubstepEnd();
   double BandBallPeak(BandBall &e);
   bool BandCycleEnd();
   bool BandBallClear(BandBall &e, const RubberContact::BallSI &bs, const RubberContact::V3 &c, bool contact, double lastGap, double relNormalVel);
   bool BandDiverged() const;              // a node of the deformable band is not a number or is farther from its rest place than the model allows
   bool BandEnsureDyn();                   // the deformable band exists (all free spans as chains); false if it cannot
   double SlingFaceN(double phi) const;    // where the plate face is at band height, m from the span's rest centreline (+ toward the ball), at plate angle phi
   double SlingFaceDN(double phi) const;   // its derivative with phi, m/rad
   int SlingChain() const;                 // index into m_dyn.chains of the span the plate faces, -1 = none
   void SlingEnsureNodes();                // band nodes at the two ends of the plate face
   bool SlingInFootprint(int chain, int node) const;
   void SlingBandStep(std::vector<std::vector<RubberContact::V3>> &fNode, double h);   // one substep: the face presses on the band nodes, the plate steps with the reaction
   void SlingStepPlate(double h, double torqueNm);
   void SlingPlace() const;   // put the plate so its face is slingGapMm from the band (from its table pose, so it can be repeated when the gap changes)
   void BandPublish();   // refresh m_info from the band state
   float BandRubberPos(int span, double f) const;   // place along the whole rubber, 0..1, of the point at fraction f of a span
   friend class BandMover;
   friend void HitRubberBandBallDeleted(HitBall *ball);
   friend void HitRubberContactRelease(HitBall *ball, const char *why);
   friend struct ContactSystemImpl;   // the contact system steps the contact-held balls (contactsystem.cpp)
   // Engine-stepped band contact: while a ball is engaged, this rubber (not a one-call solve) integrates the ball against the deformable band every
   // physics cycle, with the band's state kept between cycles. Everything else the engine does to the ball (gravity, playfield, other objects) is unchanged.
   bool BandWants(const HitBall *ball) const;
   bool BandFixedNodes() const { return m_set.bandNodeMm > 0.0 && m_slingPlate == nullptr; }   // permanent evenly spaced nodes (band_node_mm)
   bool BandBuild(const RubberContact::V3 &c, double ballR, int &chain, int &node, double sameSpot = 0.008);   // choose the free span under the ball, add its node (and the passive ones on the other spans)
   void BandRebuild();   // the deformed centreline m_dyn.pts after a node was added
   void BandSyncPts();   // node positions into m_dyn.pts
   bool BandEngage(HitBall *ball);
   void BandRelease(HitBall *ball, const char *why);
   void BandReleaseAll(const char *why);
   bool BandIsEngaged(const BallS &ball) const;
   // Engine-stepped contact (a static rubber, no deformable state): the same idea as the band. The ball is engaged and every physics cycle advances its
   // contact by one cycle's worth of substeps (a whole hit in one collision call would take 900 - 1500 substeps, 400 - 600 us).
   static bool FlipperBlocks(const HitBall *ball);   // a flipper owns the ball: no force-model contact
   bool ContactEngage(HitBall *ball, const char *how = "hit");
   void ContactReleaseAll(const char *why);
   double LayerM() const;   // thickness of the compressible layer, m
   void BandNote(const HitBall *pb) const;   // logs (twice a second at most) why a ball that touches the band is not carried by the band model
   mutable long long m_bandNoteMs = -1000000;
   void BandAdvance(const float dtime);   // fills m_bandPosts / m_bandSpans from the centreline (or the band_posts override) and logs them
   void BandBeforeCore(float dtime);   // once per cycle before the band is stepped: the sling's driver and switches
   void BandAfterCore(float dtime);    // once per cycle after: diverged check, switches, the sling's table primitives
   bool m_bandIslandStepped = false;   // unified engine: the contact system stepped the band (and its balls) this cycle; BandAdvance only finishes up

   static std::vector<HitRubber *> &Registry();   // every live HitRubber, so a solve can include the rubbers near the ball

   // resolver for band_stops entries "Post7" or "Post7, r" (hasRadius / radiusIn: an explicit radius was given; false = not found or no radius of its own;
   // x, y, r in VPX units)
   StopResolver m_stopResolver;
   std::set<std::string> m_stopWarned;   // names already reported as not found (once each)

   std::shared_ptr<RubberSurface> m_surface;   // a wall / primitive host: its faces (m_pts is then empty)
   std::string m_surfaceKind;                  // "wall" or "primitive": the ini section prefix
   Bumper *m_partBumper = nullptr;
   bool m_bumperWarned = false, m_bumperSpikeShown = false;
   JetBumper::EnginePart *m_bumperModel = nullptr;   // m_part, as the bumper model (null for a gate or spinner)
   std::shared_ptr<ContactMovingPart> m_part;   // a gate, spinner or bumper of the unified engine (MakePartHost, MakeBumperHost)
   Gate *m_partGate = nullptr;                // its table part, for the events
   Spinner *m_partSpinner = nullptr;
   float *m_partAngle = nullptr, *m_partSpeed = nullptr;
   bool m_partLogged = false;   // the part's render state (the stock mover's angle and angular speed)
   HitObject *m_eventCarrier = nullptr;        // see SetEventCarrier
   bool m_hardSurface = false;                 // a part that is a hard surface of the contact model without the rubber flag (unified engine)
   double m_partE = 0.3, m_partMu = 0.3;       // that part's elasticity and friction (for the reload of its settings)
   std::vector<Vertex3Ds> m_pts; // closed loop, last = first
   RubberContact::TubePre m_tubePre;       // what ProbeTube precomputes from m_ptsSI
   std::vector<RubberContact::V3> m_ptsSI; // the same in metres (pivot-free world frame), for the solve
   double m_radiusSI = 0.0;                // tube radius, m
   float m_radius;               // tube radius, VPX units
   std::vector<BandPost> m_bandPosts;
   std::vector<BandSpan> m_bandSpans;
   BandMover m_mover{ this };
   BandDyn m_dyn;
   BandInfo m_info;
   std::vector<double> m_cumSI;           // centreline length from the first vertex to each vertex, m (n+1 entries)
   std::vector<BandBall> m_engaged;       // the balls on the band, one per span
   mutable long long m_slowNoteMs = -1000000;   // trace: last "slow contact refused" note
   uint32_t m_contactCooldownUntilMs = 0;    // after a trip of the contact model, the stock collision for a while
   uint32_t m_bandCooldownUntilMs = 0;
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
      double closedForS = 0.0;
      double dispM = 0.0;          // band displacement toward the blade at the switch, m
      float shownTransX = 0.f;
   };
   mutable std::vector<BandSwitch> m_switches;
   mutable int m_slingSpan = -1;                 // the span the plate faces
   mutable double m_slingNx = 0.0, m_slingNy = 0.0;   // unit normal of that span, from the plate toward the band
   SlingPlate m_plate;
   SlingDriver m_driver;
   bool m_plateConfigured = false;
   bool m_slingIdle = false;
   bool m_slingEnabled = true;
   bool m_slingTriggerReq = false;   // SlingTrigger() was called; consumed by SlingAdvance
   bool m_slingTriggered = false;    // the pulse that has just started came from the script
   bool m_slingAnyClosed = false;
   SlingInfo m_slingInfo;
   bool m_slingInBand = false;                   // this cycle the plate is stepped together with the band (inside BandAdvanceCore)
   bool m_slingStopLogged = false;
   bool m_slingWasActive = false;                // log: the band was deflected last cycle
   double m_slingPkBackMm = 0.0, m_slingPkFwdMm = 0.0;   // log: the largest displacement of the plate's span since the band became active, toward the plate / toward the ball
   // the plate's geometry (SI), measured when the plate is placed: the face at band height is at originN + a / cos(phi) + h tan(phi) from the rest centreline
   mutable bool m_plateGeomValid = false;
   mutable double m_plateOriginN = 0.0, m_plateA = 0.0, m_plateH = 0.0;   // pivot axis place, face offset from the axis, band height above the axis
   mutable double m_plateRest = 0.0;                                       // the plate's angle at rest from the touching pose (RotX = 0), rad, <= 0: turned back for the gap
   mutable double m_plateF0 = 0.5, m_plateF1 = 0.5;                        // the face's extent along the span, fractions of the chord
   SlingPlate::State m_slingIdleState;
   mutable float m_slingShownRotX = 0.f;
   int m_slingFired = 0;
   double SwitchDisplacementM(const BandSwitch &sw) const;   // how far the band is pushed from the ball toward the plate at a switch, m (> 0 = toward the plate)
   mutable RubberObjectSettings m_set;
   std::filesystem::path m_ini;
   std::string m_name;
   mutable std::filesystem::file_time_type m_iniTime;
   mutable std::chrono::steady_clock::time_point m_lastCheck;
};

// A ball being deleted must not stay engaged with a band
void HitRubberBandBallDeleted(HitBall *ball);
// A force-model contact (not a band) holds this ball / releases it: the flipper takes over a ball that touches its bat.
bool HitRubberContactHolds(const HitBall *ball);
void HitRubberContactRelease(HitBall *ball, const char *why);
