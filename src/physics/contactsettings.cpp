// license:GPLv3+

#include "core/stdafx.h"
#include "physics/hitball.h"
#include "parts/pintable.h"
#include "physics/physicsini.h"
#include "physics/groupstep.h"
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
         else if (key == "roll_smooth_rad_s" && num)
            s.rollSmoothRadS = d;
         else if (key == "roll_model" && num)
         {
            double m = d;
            ClampKey("roll_model", m, 0.0, 1.0);
            s.rollModel = (int)std::lround(m);
         }
         else if (key == "eps_rest_um" && num)
         {
            double m = d;
            ClampKey("eps_rest_um", m, 1e-4, 100.0);
            s.epsRestUm = m;
         }
         else if (key == "step_swing_cap")
            s.stepSwingCap = on;
         else if (key == "term_arrival")
            s.termArrival = on;
         else if (key == "damp_dir_slope")
            s.dampDirSlope = on;
         else if (key == "lever_roll_slope")
            s.leverRollSlope = on;
         else if (key == "twist_slope")
            s.twistSlope = on;
         else if (key == "end_grad_slope")
            s.endGradSlope = on;
         else if (key == "damping_integral")
            s.dampingIntegral = on;
         else if (key == "rolling_patch")
            s.rollingPatch = on;
         else if (key == "term_force_rows")
            s.termForceRows = on;
         else if (key == "grip_in_step")
            s.gripInStep = on;
         else if (key == "patch_path")
            s.patchPath = on;
         else if (key == "hinge_piece_sum")
            s.hingePieceSum = on;
         else if (key == "part_guess_accel")
            s.partGuessAccel = on;
         else if (key == "step_model" && num)
         {
            double m = d;
            ClampKey("step_model", m, 0.0, 1.0);
            s.stepModel = (int)std::lround(m);
         }
         else if (key == "arrival_depth_um" && num)
            s.arrivalDepthUm = d;
         else if (key == "ball_reach_mm" && num)
            s.ballReachMm = d;
         else if (key == "island_reach")
            s.islandReach = on;
         else if (key == "ball_pair_reach")
            s.ballPairReach = on;
         else if (key == "sleep")
            s.sleep = on;
         else if (key == "sleep_speed_mm_s" && num)
            s.sleepSpeedMmS = d;
         else if (key == "sleep_hold_ms" && num)
            s.sleepHoldMs = d;
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
      ClampKey("roll_smooth_rad_s", s.rollSmoothRadS, 0.001, 100.0);
      ClampKey("sleep_speed_mm_s", s.sleepSpeedMmS, 0.001, 10.0);
      ClampKey("sleep_hold_ms", s.sleepHoldMs, 1.0, 10000.0);
   }
}

namespace
{
   ContactSettings cur;
   std::filesystem::path path;
   std::filesystem::file_time_type stamp{};
   bool pathKnown = false, loaded = false;
   uint32_t nextCheckMs = 0;

   // what the group step takes from the keys: renewed whenever they are (GroupStep::EngineSettings, the settings of a part stepping alone)
   void Publish() { cur.ToStep(GroupStep::EngineSettings()); }

   void LoadNow()
   {
      ContactSettings s;
      LoadContactSettings(s, path);
      if (PhysicsIni::Active(path).unified)
         s.ballContact = true;   // one contact path: on the unified engine a ball pair is never left to the stock collision (ball_contact no longer acts there)
      cur = s;
      Publish();
      PLOGI << "Contact settings loaded from " << path.string() << ": " << (s.enable ? "on" : "off") << ", dynamic ratio " << s.muDynamicRatio
            << ", rolling resistance " << s.rollingResistance << ", twist " << (s.twist ? "on" : "off") << " patch " << s.twistPatchMm << " mm, ball contact " << (s.ballContact ? "on" : "off") << ", roll_model " << s.rollModel << ", step_model " << s.stepModel << ", step_swing_cap " << (s.stepSwingCap ? 1 : 0) << ", term_arrival " << (s.termArrival ? 1 : 0) << ", eps_rest_um " << s.epsRestUm << ", hinge_piece_sum " << (s.hingePieceSum ? 1 : 0) << ", damp_dir_slope " << (s.dampDirSlope ? 1 : 0) << ", lever_roll_slope " << (s.leverRollSlope ? 1 : 0) << ", twist_slope " << (s.twistSlope ? 1 : 0) << ", end_grad_slope " << (s.endGradSlope ? 1 : 0) << ", term_force_rows " << (s.termForceRows ? 1 : 0) << ", damping_integral " << (s.dampingIntegral ? 1 : 0) << ", rolling_patch " << (s.rollingPatch ? 1 : 0) << ", patch_path " << (s.patchPath ? 1 : 0) << ", grip_in_step " << (s.gripInStep ? 1 : 0) << ", part_guess_accel " << (s.partGuessAccel ? 1 : 0);
   }
}

// The start of a play, before any part is built: [Contact] as the file is now, nothing kept from the last play or table. Game time starts again at 0, so the check in Get would not see a change until it passed the last play's
void ContactSettings::StartPlay(const std::filesystem::path& tableFile)
{
   cur = ContactSettings();
   stamp = {};
   nextCheckMs = 0;
   path = PhysicsIni::FindUnified(tableFile);
   pathKnown = !path.empty();
   loaded = true;
   std::error_code ec;
   if (pathKnown)
      stamp = std::filesystem::last_write_time(path, ec);
   if (pathKnown && !ec)
      LoadNow();
   else
      cur.enable = PhysicsIni::Active(std::filesystem::path()).contact;   // no file: the defaults (off on the stock engine)
   Publish();
}

// The group step's settings from [Contact], and its floor: a step that cannot be solved gives up at hMin (in real time it must not spend the frame halving; the fixed-step reference may)
void ContactSettings::ToStep(GroupStep::Settings& s) const
{
   s.rollSmooth = rollSmoothRadS;
   s.rollModel = rollModel;
   s.stepModel = stepModel;
   s.swingCap = stepSwingCap;
   s.termArrival = termArrival;
   s.dampDirSlope = dampDirSlope;
   s.leverRollSlope = leverRollSlope;
   s.twistSlope = twistSlope;
   s.endGradSlope = endGradSlope;
   s.termForceRows = termForceRows;
   s.dampingIntegral = dampingIntegral;
   s.rollingPatch = rollingPatch;
   s.patchPath = patchPath;
   s.gripInStep = gripInStep;
   s.partGuessAccel = partGuessAccel;
   s.epsRest = epsRestUm * 1e-6;   // (0.1 um gives 1e-7 m exactly, the Settings default)
   s.hMin = GroupStep::GroupFixedStepS() > 0.0 ? 1e-9 : 2e-6;
}

const ContactSettings& ContactSettings::Get()
{
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
         Publish();
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
         LoadNow();
      }
      else if (ec)
      {
         loaded = true; // no file: the defaults stay (off on the stock engine)
         cur = ContactSettings();
         cur.enable = PhysicsIni::Active(std::filesystem::path()).contact;
         Publish();
      }
   }
   return cur;
}

// Ported at: VisualPinball.Unity/VisualPinball.Unity/VPT/Ball/BallCollider.cs
