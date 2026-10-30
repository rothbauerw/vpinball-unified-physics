// license:GPLv3+

// The bat as a body of the contact system: state, the substep and the pivot nudge (pure SI: no VPX types; the settings file is read in contactbatini.cpp). See contactbat.h.
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


#include "contactbat.h"
#include "physics/rubbercontact.h"
#include "physics/physconst.h"

#include <cmath>
#include <cstdlib>

// ---------------------------------------------------------------------------
namespace FlipperModel
{
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
   coilWorkJ += fCoil * (-vp / 1000.0) * h;   // the coil pulls toward closure (-x); vp in mm/s (the step uses the start-of-step velocity throughout)

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
