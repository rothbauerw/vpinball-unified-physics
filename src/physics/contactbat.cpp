// license:GPLv3+

// The bat as a body of the contact system: settings file, state, the substep and the pivot nudge. See contactbat.h.
// The model is a direct port of runs/2026-09-28-rampup0/flipper_model.py (rothbauerw's real-data model):
//  - plunger: translates along the coil axis, pulled toward closure by the power + hold coils, hard steel stopper at the closed end, hard stop at the open end.
//  - flipper (bat + crank + shaft, one rigid body): the return spring pulls on it continuously; a hysteretic rubber stop acts on it below phi = 0.
//  - link: nylon link with pin slop. Loaded contact is mismatch = 0 (coil and spring both press the pins against the same hole walls); free play for 0 < mismatch < SLOP;
//    stiff damped contact beyond.
//  - wiring: series (one supply feeds power + hold; the EOS switch shorts the hold winding during the stroke) or parallel.
//  - hand-off: two-body while the coil is energized; rigid one-body (plunger slaved to the flipper) when the coil is off and the link is taut; back to two-body when the
//    coil fires again.
// All dynamics are in SI units; the owner converts at the VPX boundary (1 U = 0.53975 mm, 1 VPT = 10 ms, 1 mass unit = 80 g). Every tunable is in FlipperSettings, per
// flipper, set from the table's vpx_physics.ini ([Flippers], [Flipper:<name>]), hot-reloaded on change.

#include "core/stdafx.h"
#include "contactbat.h"
#include "physics/physicsini.h"
#include "physics/rubbercontact.h"
#include "physics/physconst.h"

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

   RubberContact::RubberParams MakeRubberParams(const FlipperSettings& s)
   {
      RubberContact::RubberParams p;
      p.eMpa = s.rubberEMpa;
      p.thicknessMm = s.rubberThicknessMm;
      p.alphaSPerM = s.rubberAlphaSPerM;
      p.betaLayer = s.rubberBetaLayer;
      p.mu = s.rubberMu;
      p.shearScale = s.rubberShearScale;
      p.shearZeta = s.rubberShearZeta;
      p.twist = s.rubberTwist;
      p.tanDelta = s.rubberTanDelta;
      p.zLoMm = s.rubberZLoMm;
      p.heightMm = s.rubberHeightMm;
      return p;
   }
}

bool FlipperApplyKey(FlipperSettings &s, const std::string &key, const std::string &val)
{
   return FlipperModel::ApplyKey(s, FlipperModel::Lower(key), FlipperModel::Trim(val));
}

static constexpr double FLIPPER_MASS_UNIT_KG = 0.080;
static constexpr double FLIPPER_U_M = 1.0625 * 0.0254 / 50.0;

// (re)builds this flipper's settings from the ini file; the model state is left untouched
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

// one SI substep of length h seconds, state at the start of the step used throughout (explicit Euler,
// same ordering as flipper_model.py)
void ContactBat::Step(double h, const double extTorqueNm, const bool solenoid)
{
   using namespace FlipperModel;
   const FlipperSettings& s = set;
   const FlipperDerived& p = der;
   const double hMs = h * 1000.0;

   const double phiDeg = phi / DEG;
   const double theta = p.theta0 + phi;
   const Kin kin = Kinematics(s, p, theta);
   const double leverMm = kin.lever;
   const double leverM = leverMm / 1000.0;
   const double xRigid = kin.x;

   const double frac = (xp - p.xMin) / s.bobbinMm;
   const double gHigh = Sq(s.x0P / (frac + s.x0P));
   const double gLow = Sq(s.x0H / (frac + s.x0H));

   // back-EMF depends on the plunger's own closing speed
   const double vPos = fmax(-(vp / 1000.0), 0.0);
   const double dLdxP = p.bP / (s.bobbinMm * Sq(frac + s.x0P));   // mH/mm
   const double dLdxH = p.bH / (s.bobbinMm * Sq(frac + s.x0H));
   const double lP = p.aP + p.bP / (frac + s.x0P);                // mH
   const double lH = p.aH + p.bH / (frac + s.x0H);
   const double tauHigh = lP / s.rP, tauLow = lH / s.rH;          // ms

   const bool coilOff = !solenoid;
   const bool preEos = phiDeg < p.eosDeg;   // EOS switch rides on the crank/lever arm, not the plunger
   double iHighNext, iLowNext;
   if (s.single)
   {
      // one winding (power values only): the supply is across R + L(x) while the coil is on; when it is off the current decays through R (V = 0).
      // There is no EOS switch and no hold winding, so iLow stays 0 and the hold force term vanishes. Same circuit equation as the series
      // branch's loop, di/dt = (V - R i - i dL/dx v) / L (mH and ms give amps), for one winding.
      const double issP = s.vSupply / s.rP;
      const double iP = iHigh * issP;
      const double vNow = coilOff ? 0.0 : s.vSupply;
      const double iNextP = fmax(0.0, iP + ((vNow - s.rP * iP - iP * dLdxP * vPos) / lP) * hMs);
      iHighNext = iNextP / issP;
      iLowNext = 0.0;
      coilOffPrev = coilOff;
   }
   else if (s.series)
   {
      // currents in amps; iHigh/iLow stay normalized to V/R of each winding
      const double issP = s.vSupply / s.rP, issH = s.vSupply / s.rH;
      const double iP = iHigh * issP, iH = iLow * issH;
      const double dLdxTot = dLdxP + dLdxH;
      double iLoop, iNextP, iNextH;
      if (coilOff)
      {
         // supply removed: both windings in ONE loop with V = 0
         if (!coilOffPrev)
         {
            iLoop = (preEos && s.eosFlux) ? (lP * iP + lH * iH) / (lP + lH) : iP;
            iHigh = iLoop / issP;
            iLow = iLoop / issH;
         }
         else
            iLoop = iP;
         iNextP = iNextH = fmax(0.0, iLoop + ((0.0 - (s.rP + s.rH) * iLoop - iLoop * dLdxTot * vPos) / (lP + lH)) * hMs);
      }
      else if (preEos)
      {
         // EOS closed: hold winding shorted (decays), power winding on the supply
         iNextP = fmax(0.0, iP + ((s.vSupply - s.rP * iP - iP * dLdxP * vPos) / lP) * hMs);
         iNextH = iH * exp(-hMs / tauLow);
      }
      else
      {
         if (wasPreEos)
         {
            // EOS just opened: windings must now share one current
            iLoop = s.eosFlux ? (lP * iP + lH * iH) / (lP + lH) : iP;
            iHigh = iLoop / issP;
            iLow = iLoop / issH;
         }
         else
            iLoop = iP;
         iNextP = iNextH = fmax(0.0, iLoop + ((s.vSupply - (s.rP + s.rH) * iLoop - iLoop * dLdxTot * vPos) / (lP + lH)) * hMs);
      }
      iHighNext = iNextP / issP;
      iLowNext = iNextH / issH;
      wasPreEos = preEos;
      coilOffPrev = coilOff;
   }
   else
   {
      const double bemfHigh = (dLdxP / s.rP) * vPos;
      const double bemfLow = (dLdxH / s.rH) * vPos;
      const double tgt = coilOff ? 0.0 : 1.0;
      if (preEos)
         iHighNext = fmax(0.0, iHigh + ((tgt - iHigh - bemfHigh * iHigh) / tauHigh) * hMs);
      else
         iHighNext = iHigh * exp(-hMs / tauHigh);
      iLowNext = fmax(0.0, iLow + ((tgt - iLow - bemfLow * iLow) / tauLow) * hMs);
      coilOffPrev = coilOff;
   }

   const double fCoil = p.strengthP * gHigh * iHigh * iHigh + p.strengthH * gLow * iLow * iLow;

   // --- hand-off between the two-body model and the rigid one-body ---
   const double mismatch = xp - xRigid;
   const double rate = vp + leverMm * omega;
   const bool goRigid = mismatch < 0.0 && fabs(rate) < s.rigidMaxRelSpeedMmS;
   if (!rigid && coilOff && goRigid)
      rigid = true;
   else if (rigid && !coilOff)
      rigid = false;   // coil fired again; the plunger is slaved in rigid mode, so mismatch = 0 here and the state is continuous

   plungerHitStop = false;

   if (rigid)
   {
      // one rigid body: the coil's generalized force on phi is fCoil*lever (closing); the link carries no separate
      // force and the plunger is slaved to the flipper angle
      const double iEff = p.iFlipper + s.plungerKg * leverM * leverM;
      const double restT = RestStopTorque(s, phi, omega, iEff);
      const double net = fCoil * leverM - SpringTorqueK(s, p, kin) - s.bDamping * omega + restT + extTorqueNm;
      const double alpha = net / iEff;
      double omegaNext = omega + alpha * h;
      double phiNext = phi + omega * h;
      const double thetaN = p.theta0 + phiNext;
      const Kin kn = Kinematics(s, p, thetaN);
      xp = kn.x;
      vp = -kn.lever * omegaNext;
      if (phiNext > s.stopDeg * DEG)   // steel end stop still applies in rigid mode (inelastic)
      {
         if (omega > 0.0)
         {
            plungerHitStop = true;
            hitL = iEff * omega;
         }
         phiNext = s.stopDeg * DEG;
         omegaNext = fmin(omegaNext, 0.0);
      }
      lastAlpha = alpha;
      omega = omegaNext;
      phi = phiNext;
      iHigh = iHighNext;
      iLow = iLowNext;
      return;
   }

   const double fLink = LinkForceN(s, p, mismatch, rate);

   // --- plunger: coil pulls toward closure (-x); hard stops at both ends ---
   const double fPlunger = -fCoil + fLink;
   const bool atClosed = xp <= p.xMin;
   const bool atOpen = xp >= p.xMax;
   double xNext, vNext;
   if (atClosed && fPlunger <= 0.0)
   {
      xNext = p.xMin;
      vNext = 0.0;
   }
   else if (atOpen && fPlunger >= 0.0)
   {
      xNext = p.xMax;
      vNext = 0.0;
   }
   else
   {
      const double a = fPlunger / s.plungerKg;                 // m/s^2
      vNext = vp + a * 1000.0 * h;
      xNext = xp + vp * h;
      if (!atClosed && xNext <= p.xMin)
      {
         if (vp < 0.0)
         {
            plungerHitStop = true;
            hitL = (p.iFlipper + s.plungerKg * leverM * leverM) * fabs(omega);
         }
         xNext = p.xMin;
         vNext = 0.0;
      }
      else if (!atOpen && xNext >= p.xMax)
      {
         xNext = p.xMax;
         vNext = 0.0;
      }
   }

   // --- flipper: return spring acts directly and continuously; constrained only through the link ---
   const double restT = RestStopTorque(s, phi, omega, p.iFlipper);
   const double netTorque = -SpringTorqueK(s, p, kin) - s.bDamping * omega + restT + fLink * leverM + extTorqueNm;
   lastAlpha = netTorque / p.iFlipper;
   const double omegaNext = omega + lastAlpha * h;
   const double phiNext = phi + omega * h;

   xp = xNext;
   vp = vNext;
   omega = omegaNext;
   phi = phiNext;
   iHigh = iHighNext;
   iLow = iLowNext;
}

// ---- flipper nudge: one-shot x pulse of the pivot ----

// Pivot displacement (m) and velocity (m/s) t seconds into the pulse: half-cosine rise over nzRise, half-cosine return over nzFall.
void ContactBat::NudgeEval(const double t, double& dispM, double& velMs) const
{
   constexpr double PI = 3.14159265358979323846;
   if (!nzActive || t < 0.0 || t >= nzRiseS + nzFallS)
      dispM = velMs = 0.0;
   else if (t < nzRiseS)
   {
      dispM = nzA * 0.5 * (1.0 - cos(PI * t / nzRiseS));
      velMs = nzA * 0.5 * (PI / nzRiseS) * sin(PI * t / nzRiseS);
   }
   else
   {
      const double u = t - nzRiseS;
      dispM = nzA * 0.5 * (1.0 + cos(PI * u / nzFallS));
      velMs = -nzA * 0.5 * (PI / nzFallS) * sin(PI * u / nzFallS);
   }
}

void ContactBat::NudgeStart(const double baseXU, const double dispM, const double riseS, const double fallS)
{
   if (nzActive)   // a pulse is still running: let it finish
      return;
   nzBaseX = baseXU;
   nzA = dispM;
   nzRiseS = riseS;
   nzFallS = fallS;
   nzT = 0.0;
   nzActive = true;
   idle = false;
}

// Advances the pulse by dtS and moves the pivot to where it is at the end of that time.
void ContactBat::NudgeStep(const double dtS, float &pivotXU)
{
   if (!nzActive)
      return;
   constexpr double U = 0.53975e-3;   // m per VPX unit
   nzT += dtS;
   double d, v;
   if (nzT >= nzRiseS + nzFallS)
   {
      nzActive = false;
      d = v = 0.0;
   }
   else
      NudgeEval(nzT, d, v);
   pivotXU = (float)(nzBaseX + d / U);
   nzVelU = (float)(v / U * DEFAULT_STEPTIME_S);
}

void ContactBat::Reset(const double phi0, const bool sol)
{
   phi = phi0;
   omega = 0.0;
   xp = FlipperModel::SliderX(set, der, der.theta0 + phi);
   vp = 0.0;
   iHigh = 0.0;
   iLow = 0.0;
   lastAlpha = 0.0;
   rigid = !sol;      // a bat at rest with the coil off is one rigid body
   wasPreEos = true;
   coilOffPrev = true;
   init = true;
}

double ContactBat::InertiaSI() const
{
   const double theta = der.theta0 + phi;
   const double leverM = FlipperModel::LeverMm(set, der, theta) / 1000.0;
   double iSi = der.iFlipper;
   if (rigid || xp - FlipperModel::SliderX(set, der, theta) <= 0.0)
      iSi += set.plungerKg * leverM * leverM;
   return iSi;
}

bool ContactBat::IdleMatches(const bool sol) const
{
   return idle && sol == idleSol && rigid == idleRigid && phi == idlePhi && omega == idleOmega && xp == idleXp && vp == idleVp && iHigh == idleIHigh && iLow == idleILow;
}

void ContactBat::IdleCapture(const bool sol)
{
   idleSol = sol;
   idleRigid = rigid;
   idlePhi = phi;
   idleOmega = omega;
   idleXp = xp;
   idleVp = vp;
   idleIHigh = iHigh;
   idleILow = iLow;
}
