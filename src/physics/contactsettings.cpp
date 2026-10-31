// license:GPLv3+

#include "core/stdafx.h"
#include "physics/hitball.h"
#include "parts/pintable.h"
#include "physics/physicsini.h"
#include <filesystem>
#include <fstream>

namespace
{
   std::string Trim(const std::string& s)
   {
      const size_t a = s.find_first_not_of(" \t\r\n");
      return a == std::string::npos ? std::string() : s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
   }

   std::string Lower(std::string s)
   {
      for (char& c : s)
         c = (char)tolower((unsigned char)c);
      return s;
   }

   void ClampKey(const char* key, double& v, const double lo, const double hi)
   {
      const double was = v;
      v = std::min(std::max(v, lo), hi);
      if (v != was)
         PLOGW << "Contact settings: " << key << " = " << was << " is out of range, using " << v;
   }

   void LoadContactSettings(ContactSettings& s, const std::filesystem::path& file)
   {
      std::ifstream in(file);
      if (!in)
         return;
      std::string line, section;
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
            section = Lower(Trim(line.substr(1, line.size() - 2)));
            continue;
         }
         if (section != "contact")
            continue; // table-wide settings only for now: [Contact] of vpx_physics.ini
         const size_t eq = line.find('=');
         const std::string key = eq == std::string::npos ? std::string() : Lower(Trim(line.substr(0, eq)));
         const std::string val = eq == std::string::npos ? std::string() : Trim(line.substr(eq + 1));
         const std::string v = Lower(val);
         const bool on = (v == "1" || v == "true" || v == "yes" || v == "on");
         char* end = nullptr;
         const double d = strtod(val.c_str(), &end);
         const bool num = end != val.c_str();
         if (key == "enable")
            s.enable = on;
         else if (key == "twist")
            s.twist = on;
         else if (key == "mu_dynamic_ratio" && num)
            s.muDynamicRatio = d;
         else if (key == "rolling_resistance" && num)
            s.rollingResistance = d;
         else if (key == "ball_contact")
            s.ballContact = on;
         else if (key == "ball_e_mpa" && num)
            s.ballEMpa = d;
         else if (key == "ball_nu" && num)
            s.ballNu = d;
         else if (key == "ball_alpha" && num)
            s.ballAlpha = d;
         else if (key == "ball_tan_delta" && num)
            s.ballTanDelta = d;
         else if (key == "ball_rolling_resistance" && num)
            s.ballRollingResistance = d;
         else if (key == "ball_mu" && num)
            s.ballMu = d;
         else if (key == "ball_substep_us" && num)
            s.ballSubstepUs = d;
         else if (key == "ball_rest_substep_us" && num)
            s.ballRestSubstepUs = d;
         else if (key == "arrival_step")
            s.arrivalStep = on;
         else if (key == "static_model")
            s.staticModel = on;
         else if (key == "arrival_depth_um" && num)
            s.arrivalDepthUm = d;
         else if (key == "ball_reach_mm" && num)
            s.ballReachMm = d;
         else if (key == "twist_patch_mm" && num)
            s.twistPatchMm = d;
         else
            PLOGW << "Contact settings: ignored line '" << line << "'";
      }
      if (!PhysicsIni::Active(file).contact)   // [Features] contact = 0 is the old enable = 0; stock engine: off
         s.enable = false;
      ClampKey("mu_dynamic_ratio", s.muDynamicRatio, 0.1, 1.0);
      ClampKey("rolling_resistance", s.rollingResistance, 0.0, 0.1);
      ClampKey("twist_patch_mm", s.twistPatchMm, 0.0, 5.0);
      ClampKey("ball_e_mpa", s.ballEMpa, 1000.0, 1.0e6);
      ClampKey("ball_nu", s.ballNu, 0.0, 0.49);
      ClampKey("ball_alpha", s.ballAlpha, 0.0, 1.0);
      ClampKey("ball_tan_delta", s.ballTanDelta, 0.0, 0.5);
      ClampKey("ball_mu", s.ballMu, 0.0, 1.5);
      ClampKey("ball_rolling_resistance", s.ballRollingResistance, 0.0, 0.1);
      ClampKey("ball_substep_us", s.ballSubstepUs, 0.2, 50.0);
      ClampKey("ball_rest_substep_us", s.ballRestSubstepUs, 0.5, 100.0);
      ClampKey("ball_reach_mm", s.ballReachMm, 0.0, 20.0);
      ClampKey("arrival_depth_um", s.arrivalDepthUm, 0.1, 100.0);
   }
}

const ContactSettings& ContactSettings::Get()
{
   static ContactSettings cur;
   static std::filesystem::path path;
   static std::filesystem::file_time_type stamp{};
   static bool pathKnown = false, loaded = false;
   static uint32_t nextCheckMs = 0;
   if (g_pplayer == nullptr)
      return cur;
   if (!pathKnown)
   {
      if (g_pplayer->m_ptable != nullptr && !g_pplayer->m_ptable->m_filename.empty())
         path = PhysicsIni::FindUnified(g_pplayer->m_ptable->m_filename);   // the table's vpx_physics.ini (or VPX_PHYSICS_INI)
      pathKnown = !path.empty();
      if (!pathKnown && !loaded)
      {
         loaded = true;   // no file: the defaults (off on the stock engine)
         cur = ContactSettings();
         cur.enable = PhysicsIni::Active(std::filesystem::path()).contact;
      }
   }
   if (pathKnown && (!loaded || g_pplayer->m_time_msec >= nextCheckMs))
   {
      nextCheckMs = g_pplayer->m_time_msec + 500;
      std::error_code ec;
      const auto mtime = std::filesystem::last_write_time(path, ec);
      if (!ec && (!loaded || mtime != stamp))
      {
         stamp = mtime;
         loaded = true;
         ContactSettings s;
         LoadContactSettings(s, path);
         cur = s;
         PLOGI << "Contact settings loaded from " << path.string() << ": " << (s.enable ? "on" : "off") << ", dynamic ratio " << s.muDynamicRatio
               << ", rolling resistance " << s.rollingResistance << ", twist " << (s.twist ? "on" : "off") << " patch " << s.twistPatchMm << " mm, ball contact " << (s.ballContact ? "on" : "off");
      }
      else if (ec)
      {
         loaded = true; // no file: the defaults stay (off on the stock engine)
         cur = ContactSettings();
         cur.enable = PhysicsIni::Active(std::filesystem::path()).contact;
      }
   }
   return cur;
}

// Ported at: VisualPinball.Unity/VisualPinball.Unity/VPT/Ball/BallCollider.cs
