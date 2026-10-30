// license:GPLv3+
// The host every surface and part of the contact system is registered under: settings, registry, reload, and the hand-over of a ball from the stock engine (contacthost.h).

#include "core/stdafx.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

#include "physics/contacthost.h"
#include "physics/contactsystem.h"
#include "physics/contactpartsini.h"
#include "physics/hitball.h"
#include "physics/hitflipper.h"
#include "physics/hitrubber.h"
#include "physics/physicsini.h"
#include "physics/PhysicsEngine.h"

namespace
{
   constexpr double U = ContactSystem::kMetresPerUnit;
   constexpr double VS = ContactSystem::kVelocityScale;
   constexpr double kMassUnitKg = ContactSystem::kMassUnitKg;

   struct NumKey { const char *name; double RubberObjectSettings::*field; };
   struct BoolKey { const char *name; bool RubberObjectSettings::*field; };
   const NumKey NUM_KEYS[] = {
      { "sling_gap_mm", &RubberObjectSettings::slingGapMm },
      { "band_model", &RubberObjectSettings::bandModel }, { "band_ring_in", &RubberObjectSettings::bandRingIn },
      { "surface_layer_mm", &RubberObjectSettings::surfaceLayerMm }, { "edge_radius_mm", &RubberObjectSettings::edgeRadiusMm },
      { "e_mpa", &RubberObjectSettings::eMpa }, { "layer_mm", &RubberObjectSettings::layerMm },
      { "alpha_s_per_m", &RubberObjectSettings::alphaSPerM }, { "restitution", &RubberObjectSettings::restitution }, { "plate_areal_density_kg_m2", &RubberObjectSettings::plateArealKgM2 }, { "plate_mass_g", &RubberObjectSettings::plateMassG }, { "plate_clearance_u", &RubberObjectSettings::plateClearanceU }, { "restitution_ref_speed_ms", &RubberObjectSettings::restitutionRefMs }, { "beta_layer", &RubberObjectSettings::betaLayer },
      { "mu", &RubberObjectSettings::mu }, { "shear_scale", &RubberObjectSettings::shearScale },
      { "shear_zeta", &RubberObjectSettings::shearZeta }, { "tan_delta", &RubberObjectSettings::tanDelta },
      { "substep_us", &RubberObjectSettings::substepUs }, { "rest_substep_us", &RubberObjectSettings::restSubstepUs }, { "rest_omega_h", &RubberObjectSettings::restOmegaH }, { "rest_speed_mm_s", &RubberObjectSettings::restSpeedMmS }, { "roll_travel_um", &RubberObjectSettings::rollTravelUm }, { "max_ms", &RubberObjectSettings::maxMs }, { "release_gap_mm", &RubberObjectSettings::releaseGapMm },
      { "band_curve_radius_mm", &RubberObjectSettings::bandCurveRadiusMm }, { "band_min_span_mm", &RubberObjectSettings::bandMinSpanMm },
      { "band_prestretch_pct", &RubberObjectSettings::bandPrestretchPct }, { "band_tension_n", &RubberObjectSettings::bandTensionN }, { "band_node_zeta", &RubberObjectSettings::bandNodeZeta },
      { "band_ring_zeta", &RubberObjectSettings::bandRingZeta }, { "band_rest_mm", &RubberObjectSettings::bandRestMm },
      { "band_substep_us", &RubberObjectSettings::bandSubstepUs }, { "band_shear_zeta", &RubberObjectSettings::bandShearZeta }, { "band_ref_hz", &RubberObjectSettings::bandRefHz },
      { "band_node_mass_g", &RubberObjectSettings::bandNodeMassG }, { "band_max_ms", &RubberObjectSettings::bandMaxMs },
      { "band_flipper_margin_u", &RubberObjectSettings::bandFlipperMarginU }, { "band_follow_mm", &RubberObjectSettings::bandFollowMm }, { "band_node_mm", &RubberObjectSettings::bandNodeMm },
      { "switch_gap_mm", &RubberObjectSettings::switchGapMm }, { "switch_min_closed_ms", &RubberObjectSettings::switchMinClosedMs },
   };
   const BoolKey BOOL_KEYS[] = {
      { "contact", &RubberObjectSettings::contact }, { "rubber", &RubberObjectSettings::surfaceRubber }, { "normal_loss", &RubberObjectSettings::normalLoss }, { "viscoelastic_damping", &RubberObjectSettings::viscoelastic }, { "slow_contact", &RubberObjectSettings::slowContact }, { "twist", &RubberObjectSettings::twist }, { "band", &RubberObjectSettings::band },
      { "band_ringing", &RubberObjectSettings::bandRinging }, { "band_follow", &RubberObjectSettings::bandFollow },      { "band_lock_z", &RubberObjectSettings::bandLockZ },
   };

   using PhysicsIni::Trim;
   using PhysicsIni::Lower;

   bool ApplyKey(RubberObjectSettings &s, const std::string &key, const std::string &val)
   {
      if (key == "band_stops")
      {
         s.bandStops = val;
         return true;
      }
      if (key == "band_posts")
      {
         s.bandPosts = val;
         return true;
      }
      if (key == "sling_plate")
      {
         s.slingPlate = val;
         return true;
      }
      if (key == "sling_switches")
      {
         s.slingSwitches = val;
         return true;
      }
      if (key == "switches")
      {
         s.switches = val;
         return true;
      }
      const std::string v = Lower(val);
      for (const BoolKey &k : BOOL_KEYS)
         if (key == k.name)
         {
            s.*(k.field) = (v == "1" || v == "true" || v == "yes" || v == "on");
            return true;
         }
      for (const NumKey &k : NUM_KEYS)
         if (key == k.name)
         {
            char *end = nullptr;
            const double d = strtod(val.c_str(), &end);
            if (end == val.c_str())
               return false;
            s.*(k.field) = d;
            return true;
         }
      if (key.rfind("sling_", 0) == 0)   // the sling's own settings, then any [Flippers] key for its plate and coil with the prefix sling_
      {
         const std::string k = key.substr(6);
         auto num = [&](double &d) {
            char *end = nullptr;
            const double x = strtod(val.c_str(), &end);
            if (end == val.c_str())
               return false;
            d = x;
            return true;
         };
         if (k == "pulse_ms")
            return num(s.sling.pulseMs);
         if (k == "coil_delay_ms")
            return num(s.sling.coilDelayMs);
         if (k == "min_closed_ms")
            return num(s.sling.minClosedMs);
         if (k == "switch_gap_mm")
            return num(s.sling.switchGapMm);
         if (k == "plate_contact_zeta")
            return num(s.sling.plateContactZeta);
         return FlipperApplyKey(s.sling.mech, k, val);
      }
      return false;
   }

   bool ApplySection(RubberObjectSettings &s, const std::filesystem::path &file, const std::string &section, const std::string &who, const char *keyPrefix = "")
   {
      std::ifstream in(file);
      if (!in)
         return false;
      std::string line, cur;
      const std::string want = Lower(section);
      while (std::getline(in, line))
      {
         const size_t c = line.find_first_of("#;");
         if (c != std::string::npos)
            line.resize(c);
         line = Trim(line);
         if (line.empty())
            continue;
         if (line.front() == '[' && line.back() == ']')
         {
            cur = Lower(Trim(line.substr(1, line.size() - 2)));
            continue;
         }
         if (cur != want)
            continue;
         const size_t eq = line.find('=');
         std::string key = eq == std::string::npos ? std::string() : Lower(Trim(line.substr(0, eq)));
         if (*keyPrefix != 0 && key.rfind("sling_", 0) != 0)   // in a [Slings] / [Sling:<name>] section the keys are written without sling_
            key = keyPrefix + key;
         if (eq == std::string::npos || !ApplyKey(s, key, Trim(line.substr(eq + 1))))
            PLOGW << "Rubber settings (" << who << "): ignored line '" << line << "'";
      }
      return true;
   }

   void Clamp(RubberObjectSettings &s, const std::string &who, const double eMaxMpa = 1000.0)
   {
      auto clampKey = [&](const char *key, double &v, double lo, double hi) {
         const double was = v;
         v = std::min(std::max(v, lo), hi);
         if (v != was)
            PLOGW << "Rubber settings (" << who << "): " << key << " = " << was << " is out of range, using " << v;
      };
      clampKey("e_mpa", s.eMpa, 0.1, eMaxMpa);
      clampKey("layer_mm", s.layerMm, 0.0, 50.0);
      clampKey("surface_layer_mm", s.surfaceLayerMm, 0.2, 50.0);
      clampKey("edge_radius_mm", s.edgeRadiusMm, 0.05, 20.0);
      clampKey("alpha_s_per_m", s.alphaSPerM, 0.0, 2.0);
      clampKey("restitution", s.restitution, 0.0, 0.99);
      clampKey("restitution_ref_speed_ms", s.restitutionRefMs, 0.05, 20.0);
      clampKey("beta_layer", s.betaLayer, 0.0, 5.0);
      clampKey("mu", s.mu, 0.0, 3.0);
      clampKey("shear_scale", s.shearScale, 0.01, 100.0);
      clampKey("shear_zeta", s.shearZeta, 0.0, 10.0);
      clampKey("tan_delta", s.tanDelta, 0.0, 2.0);
      clampKey("substep_us", s.substepUs, 1.0, 100.0);
      clampKey("rest_substep_us", s.restSubstepUs, 0.0, 500.0);
      clampKey("rest_omega_h", s.restOmegaH, 0.05, 1.5);
      clampKey("rest_speed_mm_s", s.restSpeedMmS, 0.0, 500.0);
      clampKey("roll_travel_um", s.rollTravelUm, 0.0, 1000.0);
      clampKey("max_ms", s.maxMs, 1.0, 200.0);
      clampKey("release_gap_mm", s.releaseGapMm, 0.0, 5.0);
      clampKey("band_curve_radius_mm", s.bandCurveRadiusMm, 2.0, 500.0);
      clampKey("band_min_span_mm", s.bandMinSpanMm, 0.0, 500.0);
      clampKey("switch_gap_mm", s.switchGapMm, 0.05, 20.0);
      clampKey("switch_min_closed_ms", s.switchMinClosedMs, 0.0, 1000.0);
      clampKey("band_prestretch_pct", s.bandPrestretchPct, 0.0, 50.0);
      clampKey("band_tension_n", s.bandTensionN, 0.0, 2000.0);
      clampKey("band_ref_hz", s.bandRefHz, 1.0, 10000.0);
      clampKey("band_node_zeta", s.bandNodeZeta, 0.0, 2.0);
      if (s.bandRingZeta >= 0.0)
         clampKey("band_ring_zeta", s.bandRingZeta, 0.0, 2.0);
      clampKey("band_rest_mm", s.bandRestMm, 0.001, 5.0);
      clampKey("band_substep_us", s.bandSubstepUs, 0.0, 50.0);
      if (s.bandShearZeta >= 0.0)
         clampKey("band_shear_zeta", s.bandShearZeta, 0.0, 10.0);
      clampKey("band_node_mass_g", s.bandNodeMassG, 0.0, 100.0);
      clampKey("band_max_ms", s.bandMaxMs, 1.0, 1000.0);
      clampKey("band_flipper_margin_u", s.bandFlipperMarginU, 0.0, 200.0);
      clampKey("band_follow_mm", s.bandFollowMm, 5.0, 50.0);
      if (s.bandNodeMm > 0.0)
         clampKey("band_node_mm", s.bandNodeMm, 3.0, 30.0);
      clampKey("sling_gap_mm", s.slingGapMm, 0.0, 50.0);
      clampKey("sling_pulse_ms", s.sling.pulseMs, 1.0, 500.0);
      clampKey("sling_coil_delay_ms", s.sling.coilDelayMs, 0.0, 500.0);
      clampKey("sling_min_closed_ms", s.sling.minClosedMs, 0.0, 100.0);
      clampKey("sling_switch_gap_mm", s.sling.switchGapMm, 0.05, 20.0);
      clampKey("sling_plate_contact_zeta", s.sling.plateContactZeta, 0.0, 2.0);
   }

}


RubberContact::RubberParams RubberObjectSettings::ToParams() const
{
   RubberContact::RubberParams p;
   p.eMpa = eMpa;
   p.alphaSPerM = alphaSPerM;
   p.gammaKK = viscoelastic ? gammaKK : 0.0;
   p.betaLayer = betaLayer;
   p.mu = mu;
   p.shearScale = shearScale;
   p.shearZeta = shearZeta;
   p.twist = twist;
   p.tanDelta = tanDelta;
   p.rollingMu = rollingMu;
   p.normalLoss = !band && normalLoss;   // a band keeps its tuned rebound
   return p;
}

std::filesystem::path ContactHost::IniPathFor(const std::filesystem::path &tableFile)
{
   return PhysicsIni::FindUnified(tableFile);   // the table's vpx_physics.ini (or VPX_PHYSICS_INI); empty = none
}

bool ContactHost::LoadSettings(RubberObjectSettings &s, const std::filesystem::path &file, const std::string &name)
{
   if (!ContactSystem::Unified())
      return false;   // stock engine: stock rubbers
   // the rubber model is on for every rubber the [Features] allow (no file: every model on with its defaults; a rubber says contact = 0 to opt out);
   // a sling section makes its rubber a band
   const PhysicsIni::Features feat = PhysicsIni::Active(file);
   s.contact = true;
   ApplySection(s, file, "rubbers", name);
   ApplySection(s, file, "slings", name, "sling_");
   ApplySection(s, file, "rubber:" + name, name);
   ApplySection(s, file, "sling:" + name, name, "sling_");
   Clamp(s, name);
   if (!s.slingPlate.empty())
      s.band = true;
   if (s.band ? !feat.bands : !feat.rubbers)
      s.contact = false;
   if (!feat.slings)
      s.slingPlate.clear();
   return true;
}

bool ContactHost::LoadSurfaceSettings(RubberObjectSettings &s, const std::filesystem::path &file, const char *kind, const std::string &name)
{
   if (file.empty() || !PhysicsIni::Active(file).surfaceRubbers)   // the rubber = 1 flag is in the object's own section: needs the file
      return false;
   const std::string section = std::string(kind) + ":" + name;
   ApplySection(s, file, "rubbers", name);
   s.surfaceRubber = false;   // the flag is read from the object's own section only: [Rubbers] is the default for the numbers, not a switch for every wall
   s.band = false;
   ApplySection(s, file, section, name);
   s.band = false;
   Clamp(s, name);
   s.contact = s.surfaceRubber;
   return s.surfaceRubber;
}

double ContactHost::HardSurfaceE(const HitObject *first)
{
   return first != nullptr ? (double)ElasticityWithFalloff(first->m_elasticity, first->m_elasticityFalloff, (float)ContactSystem::kVelocityScale > 0.f ? (float)(1.0 / ContactSystem::kVelocityScale) : 18.53f) : 0.3;
}

bool ContactHost::BumperWanted(const std::filesystem::path &file, const std::string &name)
{
   static const bool envOn = []() { const char *const e = getenv("VPX_BUMPERS"); return !(e != nullptr && atoi(e) == 0); }();
   if (!envOn)
      return false;
   if (!PhysicsIni::Active(file).bumpers)   // no file: on with the defaults
      return false;
   bool physical = true;   // [Bumper:<name>] physical = 0 keeps the stock bumper for this one
   if (!file.empty())
      PartsIni::ReadSections(file, "", "bumper:" + name, [&](const std::string &k, double v, const std::string &) {
         if (k == "physical")
            physical = v != 0.0;
      });
   return physical;
}

bool ContactHost::HingeWanted(const std::filesystem::path &file, const char *kind, const std::string &name)
{
   static const bool envOn = []() { const char *const e = getenv("VPX_HINGES"); return !(e != nullptr && atoi(e) == 0); }();
   if (!envOn)
      return false;
   if (!PhysicsIni::Active(file).hinges)   // no file: on with the defaults
      return false;
   RubberObjectSettings probe;
   probe.surfaceRubber = true;
   ApplySection(probe, file, std::string(kind) + ":" + name, name);
   return probe.surfaceRubber;
}

bool ContactHost::HardSurfaceWanted(const std::filesystem::path &file, const char *kind, const std::string &name)
{
   if (!PhysicsIni::Active(file).hardSurfaces)   // no file: on with the defaults
      return false;
   RubberObjectSettings probe;
   probe.surfaceRubber = true;   // a hard surface unless the object's own section says rubber = 0
   ApplySection(probe, file, std::string(kind) + ":" + name, name);
   return probe.surfaceRubber;
}

bool ContactHost::LoadHardSurfaceSettings(RubberObjectSettings &s, const std::filesystem::path &file, const char *kind, const std::string &name, const double elasticity, const double friction)
{
   s = RubberObjectSettings();
   static const bool veDefault = []() { const char *const e = getenv("VPX_VISCOELASTIC"); return !(e != nullptr && atoi(e) == 0); }();
   s.viscoelastic = veDefault;
   s.eMpa = 7875.0;   // as the playfield: E* = 10.5 GPa
   s.betaLayer = 0.0;
   s.alphaSPerM = std::min(2.0, std::max(0.02, 1.1 * (1.0 - elasticity)));   // rebound ratio ~ 1 - alpha v at about 1 m/s: the part's restitution (the stock falloff with speed is not modelled)
   s.tanDelta = 0.05;
   s.normalLoss = true;   // (the rubber default is off; the right law for hard materials is not measured, to-do 1.21)
   s.mu = friction;
   s.rollingMu = ContactSettings::Get().enable ? ContactSettings::Get().rollingResistance : 0.0;
   s.restSubstepUs = 250.0;
   s.restOmegaH = 0.8;
   s.surfaceRubber = true;
   s.band = false;
   ApplySection(s, file, "surfaces", name);
   s.band = false;
   ApplySection(s, file, std::string(kind) + ":" + name, name);
   s.band = false;
   Clamp(s, name, 1.0e6);
   s.contact = s.surfaceRubber;
   // Viscoelastic normal damping (the default for hard surfaces; VPX_VISCOELASTIC=0 or viscoelastic_damping = 0 restores the alpha law): gamma is calibrated so that a ball of the default
   // size and mass rebounds with the part's restitution (its elasticity with the stock falloff at the reference speed) at the reference speed; at other speeds the law gives
   // e ~ v^(-1/5) relative to that, not the alpha law's fall to zero.
   if (s.viscoelastic && s.surfaceRubber)
   {
      const double eT = s.restitution > 0.0 ? s.restitution : elasticity;
      static std::map<std::pair<long long, long long>, double> memo;   // gamma by (e, parameters): most parts share a few elasticities
      const RubberContact::RubberParams rp = s.ToParams();
      const long long k1 = (long long)(eT * 1e6), k2 = (long long)(s.restitutionRefMs * 1e3) * 1000003LL + (long long)(s.eMpa * 10.0);
      auto it = memo.find({ k1, k2 });
      if (it == memo.end())
         it = memo.emplace(std::make_pair(k1, k2), RubberContact::CalibrateKK(rp, 25.0 * ContactSystem::kMetresPerUnit, ContactSystem::kMassUnitKg, s.restitutionRefMs, eT)).first;
      s.gammaKK = it->second;
      PLOGI << "Hard surface " << name << ": viscoelastic damping gamma " << s.gammaKK << " N s/m^2 for restitution " << eT << " at " << s.restitutionRefMs << " m/s";
   }
   return s.surfaceRubber;
}

// The playfield material. Steel ball on a lacquered wood playfield, all values UNVERIFIED (published figures, not measured on a table):
//   e_mpa 7875    effective modulus E/(1 - 0.5^2) of the model's fixed nu = 0.5, so that E* = E/(1 - nu^2) = 10.5 GPa (steel 200 GPa / 0.3, wood plus lacquer ~10 GPa / 0.3)
//   beta_layer 0  no thin-layer stiffening (a hard, thick support)
//   alpha_s_per_m 0.5, tan_delta 0.05   normal damping (rebound ratio ~ 1 - alpha v) and loss tangent (drop tests: first rebound e ~ 0.55 - 0.65 for 0.2 - 1.3 m/s)
//   mu            the table's friction (it is the static mu of every surface)
bool ContactHost::LoadPlayfieldSettings(RubberObjectSettings &s, const std::filesystem::path &file, const double tableMu)
{
   if (!PhysicsIni::Active(file).playfieldContact)   // no file: on with the defaults
      return false;
   s = RubberObjectSettings();
   s.eMpa = 7875.0;
   s.betaLayer = 0.0;
   s.alphaSPerM = 0.5;
   s.tanDelta = 0.05;
   s.normalLoss = true;   // (the rubber default is off; wood is not measured, to-do 1.21)
   s.mu = tableMu;
   s.rollingMu = ContactSettings::Get().enable ? ContactSettings::Get().rollingResistance : 0.0;   // a ball rolling on hard wood: the measured coefficient of [Contact], not the rubber hysteresis law
   s.restSubstepUs = 250.0;   // a hard contact's own frequency (about 1 kHz) is what limits the step of a resting ball
   s.restOmegaH = 0.8;
   ApplySection(s, file, "playfield", "playfield");
   s.band = false;
   Clamp(s, "playfield", 1.0e6);
   s.surfaceRubber = true;
   s.contact = true;
   return true;
}

RubberObjectSettings ContactHost::TurntableSettings(const double mu, const double rollingResistance)
{
   RubberObjectSettings s;
   s.eMpa = 7875.0;
   s.betaLayer = 0.0;
   s.alphaSPerM = 0.5;
   s.tanDelta = 0.05;
   s.normalLoss = true;   // (as the playfield)
   s.mu = mu;
   s.rollingMu = rollingResistance >= 0.0 ? rollingResistance : (ContactSettings::Get().enable ? ContactSettings::Get().rollingResistance : 0.0);
   s.restSubstepUs = 250.0;
   s.restOmegaH = 0.8;
   s.surfaceRubber = true;
   s.contact = true;
   return s;
}

// A ball that a flipper's rubber model drives is not taken by the contact model. The other way round, a contact-held ball that touches a bat is handed to the flipper at the
// moment of contact (HitFlipper::Collide / Contact release it: the flipper's step carries the nearby rubbers, so the post stays on the force model).
bool ContactHost::FlipperBlocks(const HitBall *pball)
{
   return pball->m_rubberOwner != nullptr;
}

// A stock collider of a flagged wall / primitive found the ball: take it into the force model, as Collide() does for a rubber post.
bool ContactHost::SurfaceCollide(HitBall *const pball, const Vertex3Ds &n)
{
   if (!m_set.contact)
      return false;
   if (ContactOwns(pball->m_d) || SideHolds(pball))
      return true;   // already held (by a force-model contact or by the flipper's model, which carries this rubber): another collider of the same hit
   const float dot = -(n.Dot(pball->m_d.m_vel));
   if (dot < C_LOWNORMVEL)
      return false;   // receding or sliding: nothing to solve (a resting ball is the stock contact's)
   ReloadIfChanged();
   if (!m_set.contact)
      return false;
   const bool flipperBusy = FlipperBlocks(pball);
   if (!flipperBusy && pball->m_bandOwner == nullptr && pball->m_contact == nullptr && g_pplayer->m_time_msec >= m_contactCooldownUntilMs && ContactEngage(pball))
   {
      g_pplayer->m_physics->OnBallWallHit(*pball, n, dot);   // the force-model contact is stepped by the engine from here on (no impulse here)
      return true;
   }
   return false;
}

bool ContactHost::SlowContact(HitBall *const pball)
{
   if (!m_set.contact || !m_set.slowContact)
      return false;
   if (ContactOwns(pball->m_d) || SideHolds(pball))
      return true;   // already held (by this or another collider of the same touch, or by the flipper's model that carries this rubber)
   ReloadIfChanged();
   if (!m_set.contact || !m_set.slowContact)
      return false;
   const char *why = nullptr;
   if (pball->m_rubberOwner != nullptr)
      why = "a flipper's rubber model owns the ball";
   else if (pball->m_bandOwner != nullptr || pball->m_contact != nullptr)
      why = "another band / force-model contact owns the ball";
   else if (g_pplayer->m_time_msec < m_contactCooldownUntilMs)
      why = "cooldown";
   else if (ContactEngage(pball, "slow"))
      return true;
   else
      why = "ContactEngage declined";
   if (ContactTraceOn() && (long long)g_pplayer->m_time_msec - m_slowNoteMs >= 500)
   {
      m_slowNoteMs = g_pplayer->m_time_msec;
      PLOGI << "SlowContact " << m_name << " t=" << g_pplayer->m_time_msec << " refused: " << why;
   }
   return false;
}

bool ContactHost::SideHolds(const HitBall *const pball) const
{
   return pball->m_rubberOwner != nullptr && pball->m_rubberOwner->RubberSideHas(this);
}

bool ContactHost::ContactOwns(const BallS &ball) const
{
   return ContactSystem::Owns(ball, this);
}

bool ContactHost::ContactEngage(HitBall *const pb, const char *how)
{
   return ContactSystem::Engage(this, pb, how);
}

void ContactHost::ContactReleaseAll(const char *)
{
   ContactSystem::ReleaseHost(this);
}

ContactHost::ContactHost(IEditable *editable, const RubberObjectSettings &settings, const std::filesystem::path &iniPath, const std::string &name, const char *kind)
   : HitObject(editable), m_set(settings), m_ini(iniPath), m_name(name), m_surfaceKind(kind != nullptr ? kind : "")
{
   Registry().push_back(this);
   std::error_code ec;
   if (!m_ini.empty())
      m_iniTime = std::filesystem::last_write_time(m_ini, ec);
   m_lastCheck = std::chrono::steady_clock::now();
}

ContactHost::~ContactHost()
{
   // A ball can outlive its physics engine (a script still references it when the table closes): it must not keep a pointer to a host that is gone, or ~HitBall
   // would release it from freed memory.
   ContactSystem::ReleaseHost(this);
   auto &reg = Registry();
   reg.erase(std::remove(reg.begin(), reg.end(), this), reg.end());
}

std::vector<ContactHost *> &ContactHost::Registry()
{
   static std::vector<ContactHost *> reg;
   return reg;
}

std::string ContactHost::LeftoverHosts()
{
   const std::vector<ContactHost *> &reg = Registry();
   if (reg.empty())
      return std::string();
   std::string s = std::to_string(reg.size()) + " (";
   for (size_t i = 0; i < reg.size() && i < 6; i++)
      s += (i ? ", " : "") + reg[i]->m_surfaceKind + " " + reg[i]->m_name;
   return s + (reg.size() > 6 ? ", ...)" : ")");
}

void ContactHost::ReloadIfChanged() const
{
   if (m_ini.empty())
      return;
   const auto now = std::chrono::steady_clock::now();
   if (now - m_lastCheck < std::chrono::seconds(1))
      return;
   m_lastCheck = now;
   std::error_code ec;
   const auto mtime = std::filesystem::last_write_time(m_ini, ec);
   if (ec || mtime == m_iniTime)
      return;
   m_iniTime = mtime;
   Reload();
}

double ContactHost::LayerM() const
{
   return m_set.layerMm > 0.0 ? m_set.layerMm * 1e-3 : m_set.surfaceLayerMm * 1e-3;
}
