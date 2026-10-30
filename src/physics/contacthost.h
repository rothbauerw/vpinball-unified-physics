// license:GPLv3+

#pragma once
// A host of the contact system: the VPX-facing object that every surface or part the balls touch is registered under (docs/UnifiedContactDesign.txt 12.5).
// It holds what every host needs and nothing about a kind of object: the collision box and the Collidable switch (HitObject), the name and kind, the contact
// settings and their reload from the table's physics ini, the list the contact system walks (Registry), and the entry points through which a stock collider or
// the engine's contact pass hands a ball to the contact system (SurfaceCollide, SlowContact).
// Two kinds derive from it: HitRubber (hitrubber.h: rubber posts and walls, the faces of walls, primitives, ramps and the playfield, bands and slings) and
// ContactPartHost (contactparthost.h: a moving part of the contact system, with the table part's own adapter for its pose and events).

#include <chrono>
#include <filesystem>
#include <string>
#include <vector>

#include "physics/collide.h"
#include "physics/rubbercontact.h"
#include "physics/contactsurface.h"
#include "physics/hitsling.h"

struct RubberObjectSettings
{
   bool contact = false;          // use the force model for this rubber (else the stock triangle/edge/point colliders are built)
   // Walls and primitives flagged as rubber ([Wall:<name>] / [Primitive:<name>] rubber = 1): the same model on the object's faces. Never set by [Rubbers].
   bool surfaceRubber = false;
   double surfaceLayerMm = 5.4;   // layer thickness of a flat rubber face when layer_mm is 0: a typical ring's diameter (5.4 mm on the Gorgar rings)
   double edgeRadiusMm = 1.0;     // convex edges and corners of a rubber face are treated as rounded to this radius (Hertz effective radius)
   double eMpa = 12.0;            // effective Young's modulus
   double layerMm = 0.0;          // rubber layer thickness for the thin-layer stiffening; 0 = the tube diameter
   double alphaSPerM = 0.18;      // normal damping; rebound ratio ~ 1 - alpha*v (0.18 fits the user's rebound targets 0.97 / 0.87 / 0.64 at 0.1 / 0.86 / 3 m/s, 2026-10-29; not measured)
   double betaLayer = 0.33;       // thin-layer stiffening coefficient
   double plateArealKgM2 = 0.0;   // a gate or spinner plate: its areal density (kg/m^2), or its mass (g) when plate_mass_g > 0; 0 = the default of the kind (UNVERIFIED)
   double plateMassG = 0.0;
   double plateClearanceU = 8.0;  // the plate's axis is at least this far above the top of the ball (50 U): the plate swings clear of the ball instead of resting on its top
   bool viscoelastic = false;     // viscoelastic (Kuwabara-Kono) normal damping calibrated to `restitution` at `restitutionRefMs` instead of alpha_s_per_m: the rebound ratio falls only slowly with speed (hard surfaces)
   double restitution = 0.0;      // the rebound ratio at the reference speed; 0 = from the part's elasticity (with its falloff at that speed)
   double restitutionRefMs = 1.0; // the reference speed, m/s
   double gammaKK = 0.0;          // the calibrated damping coefficient (not an ini key)
   double mu = 0.8;               // rubber / steel friction (black / white rubber; user, 2026-10-30)
   double shearScale = 1.0;       // multiplier on the shear stiffness
   double shearZeta = 0.2;        // shear damping ratio
   bool twist = true;             // drilling friction
   double tanDelta = 0.05;        // loss tangent for the rolling resistance, 0 = off
   double rollingMu = 0.0;        // floor of the rolling-resistance coefficient (the playfield: [Contact] rolling_resistance, because the hysteresis law below gives almost nothing on hard wood); 0 = none
   bool normalLoss = false;       // the same tan_delta also damps the normal force (a dashpot at the contact's own natural frequency); off for rubbers (2026-10-29: it caps a gentle bounce near 0.78), on for the hard surfaces, playfield and turntable (their loaders)
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
   double bandModel = 1.0;        // 1 = the cord in the group step (bandcord.h, GroupStepDesign 8, 11); 0 = the old node band, for comparison only
   double bandRingIn = 0.0;       // the ring's inner diameter, inches: the stretch at rest comes from it and the posts; 0 = band_prestretch_pct
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

class ContactHost : public HitObject, public ContactSurface
{
public:
   // kind: the ini section prefix of the object ("wall", "primitive", "ramp", "playfield", "gate", "spinner", "bumper", "turntable"; empty for a rubber post)
   ContactHost(IEditable *editable, const RubberObjectSettings &settings, const std::filesystem::path &iniPath, const std::string &name, const char *kind);
   ~ContactHost() override;

   // every live host, so a solve can include the surfaces near the ball
   static std::vector<ContactHost *> &Registry();
   static std::vector<ContactHost *> &AllHosts() { return Registry(); }
   // Hosts that are still registered: called when a new physics engine starts, where there must be none (a host the previous table did not release is tested against every
   // ball of this one). Empty if there are none, else "N (kind name, ...)".
   static std::string LeftoverHosts();

   // ---- settings (the table's vpx_physics.ini) ----
   static std::filesystem::path IniPathFor(const std::filesystem::path &tableFile);
   // Reads [Rubbers], [Slings], [Rubber:<name>], [Sling:<name>] of vpx_physics.ini over the built-in defaults (no file: the defaults). False on the stock engine.
   static bool LoadSettings(RubberObjectSettings &s, const std::filesystem::path &file, const std::string &name);
   // A wall or primitive: [Rubbers] as the defaults, then [<kind>:<name>] (kind = "wall" or "primitive"). True if the object is flagged rubber = 1 and
   // [Features] surface_rubbers allows it; s.contact is then set. Only a unified vpx_physics.ini can flag surfaces.
   static bool LoadSurfaceSettings(RubberObjectSettings &s, const std::filesystem::path &file, const char *kind, const std::string &name);
   // Unified engine: any other wall or primitive is a hard surface of the contact model. HardSurfaceWanted: the feature is on and the object does not opt out ([Wall:<name>] rubber = 0).
   // LoadHardSurfaceSettings: the playfield's material (steel ball on lacquered wood) with the part's own friction mu and its restitution e turned into the normal damping
   // (rebound ratio ~ 1 - alpha v at 1 m/s), then [Surfaces] and the object's own section on top. PROVISIONAL, unverified: the mapping from the stock elasticity.
   static bool HardSurfaceWanted(const std::filesystem::path &file, const char *kind, const std::string &name);
   static bool LoadHardSurfaceSettings(RubberObjectSettings &s, const std::filesystem::path &file, const char *kind, const std::string &name, double elasticity, double friction);
   static double HardSurfaceE(const HitObject *first);   // the part's restitution at 1 m/s: its elasticity with the stock falloff (0.3 without a collider)
   static bool HingeWanted(const std::filesystem::path &file, const char *kind, const std::string &name);   // like HardSurfaceWanted, for [Features] hinges (gates and spinners)
   static bool BumperWanted(const std::filesystem::path &file, const std::string &name);   // [Features] bumpers on, VPX_BUMPERS not 0, and the bumper's own section does not say physical = 0
   // The playfield as a surface on the force model. s gets the material of the [Playfield] section (defaults: a steel ball on lacquered wood, all UNVERIFIED), tableMu the
   // table's friction. False if [Features] playfield_contact is off.
   static bool LoadPlayfieldSettings(RubberObjectSettings &s, const std::filesystem::path &file, double tableMu);
   // A turntable's disc: a hard textured surface (the playfield's material with the given friction). Never off.
   static RubberObjectSettings TurntableSettings(double mu, double rollingResistance);
   void SetHardSurface(double elasticity, double friction) { m_hardSurface = true; m_partE = elasticity; m_partMu = friction; }
   const RubberObjectSettings &ContactSettings() const { return m_set; }
   // The settings file changed (checked at most once a second): read the settings again (Reload: what a host reloads depends on its kind)
   void ReloadIfChanged() const;

   // ---- the ball handed over by the stock engine ----
   // A stock collider of the object found the ball. True if the force model holds it from now on (the stock impulse is skipped).
   bool SurfaceCollide(HitBall *ball, const Vertex3Ds &n);
   // A resting or slow touch (the engine's contact pass). True if the force model holds the ball from now on (the stock contact is skipped).
   bool SlowContact(HitBall *ball);
   bool ContactOwns(const BallS &ball) const;   // this host is one of the surfaces of a force-model contact that holds the ball
   bool SideHolds(const HitBall *ball) const;   // a flipper's rubber model drives this ball and has this host in its set: the stock contact stays out
   // Once per physics interval (the host's mover): a moving part advances alone for what the contact system did not step of the interval (ContactPartHost); nothing for a surface.
   virtual void PartAdvance(float /*dtime*/) {}

   // ---- ContactSurface ----
   const std::string &Name() const override { return m_name; }
   bool IsPlayfieldKind() const { return m_surfaceKind == "playfield"; }
   bool IsPlayfield() const override { return IsPlayfieldKind(); }
   bool IsFloor() const override { return m_surfaceKind == "playfield" || m_surfaceKind == "turntable"; }   // a surface a ball rests on: gravity acts in the contact substeps while one is in the contact
   double ReleaseGapM() const override { return m_set.releaseGapMm * 1e-3; }
   const RubberObjectSettings *Settings() const override { return &m_set; }
   // The height (VPX units) of the highest upward-facing face of this host at table point (x, y); false if there is none there. Used to seat a turntable in the playfield.
   virtual bool SurfaceHeightAt(double /*xU*/, double /*yU*/, double & /*zU*/) const { return false; }
   virtual double LayerM() const;   // thickness of the compressible layer, m

protected:
   static bool FlipperBlocks(const HitBall *ball);   // a flipper owns the ball: no force-model contact
   bool ContactEngage(HitBall *ball, const char *how = "hit");
   void ContactReleaseAll(const char *why);
   virtual void Reload() const = 0;

   mutable RubberObjectSettings m_set;
   std::filesystem::path m_ini;
   std::string m_name;
   std::string m_surfaceKind;                  // the ini section prefix (see the constructor)
   bool m_hardSurface = false;                 // a part that is a hard surface of the contact model without the rubber flag (unified engine)
   double m_partE = 0.3, m_partMu = 0.3;       // that part's elasticity and friction (for the reload of its settings)
   uint32_t m_contactCooldownUntilMs = 0;      // after a trip of the contact model, the stock collision for a while
   mutable long long m_slowNoteMs = -1000000;  // trace: last "slow contact refused" note

private:
   mutable std::filesystem::file_time_type m_iniTime;
   mutable std::chrono::steady_clock::time_point m_lastCheck;
   friend struct ContactSystemImpl;   // the contact system steps the contact-held balls (contactsystem.cpp)
};
