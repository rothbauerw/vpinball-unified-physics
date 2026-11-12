// license:GPLv3+

#pragma once
// A host of the contact system: the VPX-facing object that every surface or part the balls touch is registered under (docs/unified-physics/archive/UnifiedContactDesign.txt 12.5).
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
   double shearScale = 0.1;       // sideways_stiffness: multiplier on Mindlin's 8 G* a. Mindlin assumes a thick solid body; a ring's section also shears and turns on its post, in series with the patch. Rough calculation 0.14 - 0.25 (runs/2026-10-30-rubber-review/notes.md, assumed geometry, not measured); 0.1 fits the user's rebound on a ring at 35 U with mu 0.8. A hard surface sets 1
   double shearTan = -1.0;        // the sideways spring's loss tangent: < 0 = tan_delta (a rubber: the same material). Not an ini key; hard surfaces, the playfield and the disc set their own (rubbercontact.h RubberParams)
   bool twist = true;             // drilling friction
   double tanDelta = 0.05;        // loss tangent for the rolling resistance, 0 = off
   double rollingMu = 0.0;        // floor of the rolling-resistance coefficient (the playfield: [Contact] rolling_resistance, because the hysteresis law below gives almost nothing on hard wood); 0 = none
   bool normalLoss = false;       // the same tan_delta also damps the normal force (a dashpot at the contact's own natural frequency); off for rubbers (2026-10-29: it caps a gentle bounce near 0.78), on for the hard surfaces, playfield and turntable (their loaders)
   double maxMs = 30.0;           // give up after this much simulated contact time
   double releaseGapMm = 0.2;     // an engaged ball is released only once it is farther than this from every rubber (and receding): a resting ball must not engage and release over and over

   // Band (a rubber stretched between posts): the posts and free spans are found at load (logged); a ball on a span deflects it (cordfixed.h: the fixed-height cord)
   bool band = false;             // this rubber is a band: found from its shape at load (DetectBand), not an ini key
   double bandRingIn = 0.0;       // the ring's inner diameter, inches: the stretch at rest comes from it and the posts; 0 = the standard size that gives a tension within the range (HitRubber::CordRing)
   double bandSectionMm = 0.0;     // the cord's section (its diameter), mm, for physics and drawing; 0 = the rubber's Thickness in the editor
   double bandStretchAlpha = 0.09; // Hunt-Crossley on the cord's stretch, s/m: half the ring's 0.18, for half its drop-off with speed. Fitted to the user's best guess (about 0.8 at 3 m/s), not measured
   double bandTwistScale = 1.0;    // the twist spring against its calculated stiffness (kT / 6)
   bool bandTwistLoss = true;      // the twist spring takes the rubber's tan_delta
   bool bandPlateReach = true;     // band_plate_reach (default 1 since 2026-11-05): 1 = a ball is taken into a sling's band within the plate's own stroke at the cord (CordFixed::PlateReachM), and the band's near search reaches as far as its engage rule needs; 0 = the fixed 30 mm and the search's 4 U pad
   double bandRefHz = 20.0;       // frequency at which tan_delta (the rubber's own loss tangent, shared with the rings) applies to the stretch
   double bandCurveRadiusMm = 40.0; // a centreline vertex whose local bend radius is below this belongs to a post wrap
   double bandMinSpanMm = 20.0;   // a straight run shorter than this between two posts is not a free span
   // Sling (a band with a kicker plate behind it): the plate and the switch blades are table primitives named here. The engine places the plate at load
   // (and again when the gap changes) so its face is sling_gap_mm from the band, and turns collision off on all of them. The table file is never changed.
   std::string slingPlate;        // name of the plate primitive; its position is the pivot. Empty = this band is not a sling
   std::string slingSwitches;     // names of the switch blade primitives, comma separated (wired in parallel)
   std::string switches;          // names of switch blade primitives behind the band that only the script reads (BandSwitchClosed, SwitchClosed / SwitchOpened), comma separated
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

   virtual class HitRubber *BandHost() { return nullptr; }   // a rubber: itself

   // every live host, so a solve can include the surfaces near the ball
   static std::vector<ContactHost *> &Registry();
   static std::vector<ContactHost *> &AllHosts() { return Registry(); }
   int m_regOrder = 0;   // shared_finder P2b: this host's place in the registry when the finder was built (a ball's surfaces are kept in this order, as the registry loops gave them)
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
   static bool IsSurfaceKey(const std::string &key);   // a key of the surface material (LoadHardSurfaceSettings): a part's own loader does not call it unknown
   static double HardSurfaceE(const HitObject *first);   // the part's restitution at 1 m/s: its elasticity with the stock falloff (0.3 without a collider)
   static bool HingeWanted(const std::filesystem::path &file, const char *kind, const std::string &name);   // like HardSurfaceWanted, for [Features] hinges (gates and spinners)
   // The playfield as a surface on the force model. s gets the material of the [Playfield] section (defaults: a steel ball on lacquered wood, all UNVERIFIED), tableMu the
   // table's friction. False if [Features] playfield_contact is off.
   static bool LoadPlayfieldSettings(RubberObjectSettings &s, const std::filesystem::path &file, double tableMu);
   // A turntable's disc: a hard textured surface (the playfield's material with the given friction). Never off.
   static RubberObjectSettings TurntableSettings(double mu, double rollingResistance);
   void SetHardSurface(double elasticity, double friction) { m_hardSurface = true; m_partE = elasticity; m_partMu = friction; }
   const RubberObjectSettings &ContactSettings() const { return m_set; }
   // The settings file changed (checked at most once a second): read the settings again (Reload: what a host reloads depends on its kind)
   void ReloadIfChanged() const;

   // ---- the ball found by a stock collider ----
   // A stock collider of the object found the ball. True if the force model holds it from now on (the stock impulse is skipped). False only for an object not converted (a setting) and the
   // playfield over a kicker's hole.
   bool SurfaceCollide(HitBall *ball, const Vertex3Ds &n);
   // A resting or slow touch (the engine's contact pass). True if the force model holds the ball from now on (the stock contact is skipped).
   bool SlowContact(HitBall *ball);
   bool ContactOwns(const BallS &ball) const;   // this host is one of the surfaces of a force-model contact that holds the ball
   // Once per physics interval (the host's mover): a moving part advances alone for what the contact system did not step of the interval (ContactPartHost); nothing for a surface.
   // The interval is the engine's current sub-interval in seconds, as ContactSystem::Step stored it (ContactSystem::SubIntervalS).
   virtual void PartAdvance() {}

   // ---- ContactSurface ----
   const std::string &Name() const override { return m_name; }
   bool IsPlayfieldKind() const { return m_surfaceKind == "playfield"; }
   bool IsPlayfield() const override { return IsPlayfieldKind(); }
   bool IsFloor() const override { return m_surfaceKind == "playfield" || m_surfaceKind == "turntable"; }   // a surface a ball rests on: gravity acts in the contact substeps while one is in the contact
   double ReleaseGapM() const override { return m_set.releaseGapMm * 1e-3; }
   bool GroupModel() const override;
   const RubberObjectSettings *Settings() const override { return &m_set; }
   // The height (VPX units) of the highest upward-facing face of this host at table point (x, y); false if there is none there. Used to seat a turntable in the playfield.
   virtual bool SurfaceHeightAt(double /*xU*/, double /*yU*/, double & /*zU*/) const { return false; }
   virtual double LayerM() const;   // thickness of the compressible layer, m

protected:
   bool ContactEngage(HitBall *ball, const char *how = "hit");
   void ContactReleaseAll(const char *why);
   virtual void Reload() const = 0;

   mutable RubberObjectSettings m_set;
   std::filesystem::path m_ini;
   std::string m_name;
   std::string m_surfaceKind;                  // the ini section prefix (see the constructor)
   bool m_hardSurface = false;                 // a part that is a hard surface of the contact model without the rubber flag (unified engine)
   double m_partE = 0.3, m_partMu = 0.3;       // that part's elasticity and friction (for the reload of its settings)

private:
   mutable std::filesystem::file_time_type m_iniTime;
   mutable std::chrono::steady_clock::time_point m_lastCheck;
   friend struct ContactSystemImpl;   // the contact system steps the contact-held balls (contactsystem.cpp)
};
