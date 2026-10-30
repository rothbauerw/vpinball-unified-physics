// license:GPLv3+

// The bat's settings: defaults, then [Flippers], then [Flipper:<name>] of the table's vpx_physics.ini (VPX side: file, log). The bat itself is in contactbat.cpp (pure SI).

#include "core/stdafx.h"
#include "contactbat.h"
#include "physics/physicsini.h"

#include <cstdlib>
#include <fstream>

// ---------------------------------------------------------------------------
namespace FlipperModel
{
   // ---- settings: defaults, then [Flippers], then [Flipper:<name>] of the table's vpx_physics.ini ----

   struct NumKey { const char *name; double FlipperSettings::*field; };
   struct BoolKey { const char *name; bool FlipperSettings::*field; };

   static const NumKey NUM_KEYS[] = {
      { "l_rod_mm", &FlipperSettings::lRodMm }, { "r_crank_mm", &FlipperSettings::rCrankMm },
      { "center_deg", &FlipperSettings::centerDeg }, { "stop_deg", &FlipperSettings::stopDeg },
      { "i_bat", &FlipperSettings::iBat }, { "i_crank", &FlipperSettings::iCrank },
      { "plunger_kg", &FlipperSettings::plungerKg }, { "b_damping", &FlipperSettings::bDamping },
      { "v_supply", &FlipperSettings::vSupply }, { "bobbin_mm", &FlipperSettings::bobbinMm },
      { "eos_offset_deg", &FlipperSettings::eosOffsetDeg },
      { "r_power", &FlipperSettings::rP }, { "l_air_power", &FlipperSettings::lAirP },
      { "l_arm_power", &FlipperSettings::lArmP }, { "x0_power", &FlipperSettings::x0P },
      { "r_hold", &FlipperSettings::rH }, { "l_air_hold", &FlipperSettings::lAirH },
      { "l_arm_hold", &FlipperSettings::lArmH }, { "x0_hold", &FlipperSettings::x0H },
      { "preload_n", &FlipperSettings::preloadN }, { "return_ms", &FlipperSettings::returnMs },
      { "k_spring_n_per_mm", &FlipperSettings::kSpringNPerMm },
      { "rest_k_nm_per_rad", &FlipperSettings::restKNmPerRad }, { "rest_zeta", &FlipperSettings::restZeta }, { "rest_n", &FlipperSettings::restN },
      { "slop_deg", &FlipperSettings::slopDeg }, { "k_contact_n_per_mm", &FlipperSettings::kContactNPerMm },
      { "contact_zeta", &FlipperSettings::contactZeta }, { "stiff_zeta", &FlipperSettings::stiffZeta },
      { "rigid_max_rel_speed_mm_s", &FlipperSettings::rigidMaxRelSpeedMmS },
      { "rubber_e_mpa", &FlipperSettings::rubberEMpa },
      { "rubber_thickness_mm", &FlipperSettings::rubberThicknessMm }, { "rubber_alpha_s_per_m", &FlipperSettings::rubberAlphaSPerM },
      { "rubber_thin_layer_beta", &FlipperSettings::rubberBetaLayer }, { "rubber_mu", &FlipperSettings::rubberMu },
      { "rubber_shear_scale", &FlipperSettings::rubberShearScale }, { "rubber_shear_zeta", &FlipperSettings::rubberShearZeta },
      { "rubber_z_lo_mm", &FlipperSettings::rubberZLoMm }, { "rubber_height_mm", &FlipperSettings::rubberHeightMm },
      { "rubber_substep_us", &FlipperSettings::rubberSubstepUs }, { "rubber_release_steps", &FlipperSettings::rubberReleaseSteps }, { "rubber_release_gap_mm", &FlipperSettings::rubberReleaseGapMm },
      { "rubber_tan_delta", &FlipperSettings::rubberTanDelta },
      { "nudge_gain", &FlipperSettings::nudgeGain }, { "nudge_rise_ms", &FlipperSettings::nudgeRiseMs }, { "nudge_fall_ms", &FlipperSettings::nudgeFallMs },
   };
   static const BoolKey BOOL_KEYS[] = {
      { "free_play_damping", &FlipperSettings::freePlayDamping },
      { "rubber_contact", &FlipperSettings::rubberContact }, { "rubber_twist", &FlipperSettings::rubberTwist },
   };

   inline std::string Trim(const std::string& s)
   {
      const size_t a = s.find_first_not_of(" \t\r\n");
      if (a == std::string::npos)
         return std::string();
      return s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
   }
   inline std::string Lower(std::string s)
   {
      for (char& c : s)
         c = (char)tolower((unsigned char)c);
      return s;
   }

   inline bool ApplyKey(FlipperSettings& s, const std::string& key, const std::string& val)
   {
      const std::string v = Lower(val);
      if (key == "wiring")
      {
         if (v != "series" && v != "parallel" && v != "single")
            return false;
         s.series = (v == "series");
         s.single = (v == "single");
         return true;
      }
      if (key == "eos_open_rule")
      {
         if (v != "flux" && v != "continuity")
            return false;
         s.eosFlux = (v == "flux");
         return true;
      }
      if (key == "nudge_target")
      {
         s.nudgeTarget = val;
         return true;
      }
      for (const BoolKey& k : BOOL_KEYS)
         if (key == k.name)
         {
            s.*(k.field) = (v == "1" || v == "true" || v == "yes" || v == "on");
            return true;
         }
      for (const NumKey& k : NUM_KEYS)
         if (key == k.name)
         {
            char *end = nullptr;
            const double d = strtod(val.c_str(), &end);
            if (end == val.c_str())
               return false;
            s.*(k.field) = d;
            return true;
         }
      return false;
   }

   // applies one section ("" never matches); returns false if the file cannot be read
   inline bool ApplySection(FlipperSettings& s, const std::filesystem::path& file, const std::string& section, const std::string& who)
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
         if (eq == std::string::npos || !ApplyKey(s, Lower(Trim(line.substr(0, eq))), Trim(line.substr(eq + 1))))
            PLOGW << "Flipper settings (" << who << "): ignored line '" << line << "'";
      }
      return true;
   }

   // Keeps the rubber settings in a range where the contact model is defined and stable. Returns the number of values changed;
   // each change is reported through `warn`, so a typo in the ini is visible instead of silently ignored.
   template <class Warn>
   inline int ClampRubberSettings(FlipperSettings& s, Warn warn)
   {
      int changed = 0;
      auto clampKey = [&](const char *key, double& v, double lo, double hi) {
         double c = v < lo ? lo : (v > hi ? hi : v);
         if (std::isnan(v))
            c = lo;
         if (c != v)
         {
            warn(key, v, c);
            v = c;
            changed++;
         }
      };
      clampKey("rubber_e_mpa", s.rubberEMpa, 0.1, 1000.0);
      clampKey("rubber_thickness_mm", s.rubberThicknessMm, 0.2, 20.0);
      clampKey("rubber_alpha_s_per_m", s.rubberAlphaSPerM, 0.0, 2.0);
      clampKey("rubber_thin_layer_beta", s.rubberBetaLayer, 0.0, 5.0);
      clampKey("rubber_mu", s.rubberMu, 0.0, 3.0);
      clampKey("rubber_shear_scale", s.rubberShearScale, 0.01, 100.0);
      clampKey("rubber_shear_zeta", s.rubberShearZeta, 0.0, 10.0);
      clampKey("rubber_z_lo_mm", s.rubberZLoMm, 0.0, 60.0);
      clampKey("rubber_height_mm", s.rubberHeightMm, 1.0, 60.0);
      clampKey("rubber_substep_us", s.rubberSubstepUs, 2.0, 100.0);
      clampKey("rubber_release_steps", s.rubberReleaseSteps, 1.0, 1000.0);
      clampKey("rubber_release_gap_mm", s.rubberReleaseGapMm, 0.0, 5.0);
      clampKey("rubber_tan_delta", s.rubberTanDelta, 0.0, 2.0);
      clampKey("nudge_gain", s.nudgeGain, 0.0, 10000.0);
      clampKey("nudge_rise_ms", s.nudgeRiseMs, 0.1, 100.0);
      clampKey("nudge_fall_ms", s.nudgeFallMs, 0.1, 500.0);
      return changed;
   }

}

bool FlipperApplyKey(FlipperSettings &s, const std::string &key, const std::string &val)
{
   return FlipperModel::ApplyKey(s, FlipperModel::Lower(key), FlipperModel::Trim(val));
}

void ContactBat::LoadSettings(const std::string &who, const std::filesystem::path &tableFile, const bool force)
{
   using namespace FlipperModel;
   std::error_code ec;
   if (iniPath.empty())
   {
      if (!tableFile.empty())
         iniPath = PhysicsIni::FindUnified(tableFile);   // the table's own vpx_physics.ini (or VPX_PHYSICS_INI); none = the built-in defaults
   }
   if (iniPath.empty())
   {
      if (force)
         Derive(set, der);
      return;
   }
   const auto mtime = std::filesystem::last_write_time(iniPath, ec);
   if (!force && (ec || mtime == iniTime))
      return;
   iniTime = ec ? std::filesystem::file_time_type{} : mtime;
   FlipperSettings s;
   ApplySection(s, iniPath, "flippers", who);
   ApplySection(s, iniPath, "flipper:" + who, who);
   if (!PhysicsIni::Active(iniPath).flipperRubber)
      s.rubberContact = false;
   ClampRubberSettings(s, [&who](const char *key, double was, double now) { PLOGW << "Flipper settings (" << who << "): " << key << " = " << was << " is out of range, using " << now; });
   set = s;
   Derive(set, der);
   idle = false;   // new settings: the old equilibrium no longer holds
   PLOGI << "Flipper settings for " << who << " loaded from " << iniPath.string() << ": " << (s.single ? "single" : (s.series ? "series" : "parallel"))
         << ", K_spring " << der.kSpringNPerMm << " N/mm, preload " << s.preloadN << " N, rest stop K " << s.restKNmPerRad
         << " zeta " << s.restZeta << " n " << s.restN << ", rubber contact " << (s.rubberContact ? "on" : "off") << ", tan delta " << s.rubberTanDelta;
}
