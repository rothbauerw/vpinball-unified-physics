// license:GPLv3+

#pragma once

#include <cmath>
#include <cstdint>
#include <string>

// The coil / plunger / link mechanics, shared by the model flipper (hitflipper.cpp) and the sling plate (hitsling.cpp): the settings, the
// quantities derived from them, and the kinematics and force laws of the slider-crank link, the link's slop and stops, the return spring and the
// rest stop. Pure functions of the settings and the state; nothing here owns state or knows about flippers or slings.

// flipper model settings, one set per flipper (see FlipperLoadSettings in hitflipper.cpp: defaults, then
// [Flippers] of vpx_physics.ini, then the [Flipper:<name>] section). Units: mm, ohm, mH, N, kg, degrees unless named otherwise.
struct FlipperSettings
{
   bool series = true;               // wiring: series (EOS shorts the hold winding) or parallel
   bool single = false;              // wiring = single: one winding (the "power" values r_power, l_air_power, l_arm_power, x0_power), no hold winding, no EOS switch; series/parallel are then false
   bool eosFlux = true;              // series, EOS-open current rule: flux conserved (true) or power-winding current kept
   bool freePlayDamping = true;
   // linkage
   double lRodMm = 19.05, rCrankMm = 25.4, centerDeg = 90.0, stopDeg = 41.5;
   // inertia / mass
   double iBat = 2.299e-5, iCrank = 6.81e-6, plungerKg = 0.0413355, bDamping = 0.0;
   // coils
   double vSupply = 28.0, bobbinMm = 35.0, eosOffsetDeg = 6.0;
   double rP = 1.25, lAirP = 0.673, lArmP = 2.205, x0P = 0.2;
   double rH = 25.05, lAirH = 4.47, lArmH = 15.68, x0H = 0.0048;
   // return spring: K is solved so the return takes returnMs, unless kSpringNPerMm > 0
   double preloadN = 0.5, returnMs = 60.0, kSpringNPerMm = 0.0;
   // rest-position hysteretic rubber stop
   double restKNmPerRad = 15.0, restZeta = 0.45, restN = 1.0;   // rest stop: spring (N m/rad at 1.4 deg), damping ratio, force exponent (1 = linear)
   // link slop and contacts
   double slopDeg = 6.0, kContactNPerMm = 20.0, contactZeta = 0.1, stiffZeta = 0.12;
   // rigid one-body hand-off
   double rigidMaxRelSpeedMmS = 50.0;   // two-body -> rigid only when the link is taut and slower than this
   // ball / rubber contact (physics/rubbercontact.h). rubberContact = false keeps the impulse path. All values are
   // unverified defaults meant to be edited; see MakeRubberParams / ClampRubberSettings in hitflipper.cpp.
   bool rubberContact = false;          // use the rubber force model for ball contacts on this flipper
   double rubberEMpa = 12.0;            // effective Young's modulus
   double rubberThicknessMm = 3.2;      // rubber layer thickness
   double rubberAlphaSPerM = 0.05;      // normal damping; rebound ratio ~ 1 - alpha*v
   double rubberBetaLayer = 0.33;       // thin-layer stiffening coefficient
   double rubberMu = 0.9;               // rubber / steel friction
   double rubberShearScale = 1.0;       // multiplier on the shear stiffness
   double rubberShearZeta = 0.2;        // shear damping ratio
   bool rubberTwist = true;             // resistance to spin about the contact normal
   double rubberTanDelta = 0.05;        // rubber loss tangent for the rolling resistance, 0 = off. PROVISIONAL: starts low on purpose (typical elastomer 0.1-0.3), to be tuned
   double rubberZLoMm = 9.144;          // rubber bottom above the playfield surface (0.36 in)
   double rubberHeightMm = 11.176;      // rubber height (0.44 in)
   double rubberReleaseSteps = 20.0;    // substeps the ball must stay clear of the rubber and receding before the engagement ends
   double rubberReleaseGapMm = 0.2;     // ...and be farther than this from the rubber surface (release hysteresis for a resting ball)
   double rubberSubstepUs = 10.0;       // substep while a ball is engaged with the rubber (bat and ball step together on this)
   // flipper nudge: when this flipper's plunger hits the coil stop, the named flipper's pivot is displaced along x, away from this one, by a
   // single pulse (fast rise, slow return) of nudgeGain * (angular momentum at the stop, kg m^2/s) millimetres. nudgeGain 0 = off.
   std::string nudgeTarget;             // name of the flipper that is nudged ("" = none)
   double nudgeGain = 50.0;             // mm per kg m^2/s (inert until nudge_target is set)
   double nudgeRiseMs = 3.0, nudgeFallMs = 5.0;
};

// quantities derived from FlipperSettings (rebuilt whenever the settings change)
struct FlipperDerived
{
   double offsetMm = 0.0, theta0 = 0.0, xMin = 0.0, xMax = 0.0, xRest = 0.0;
   double aP = 0.0, bP = 0.0, aH = 0.0, bH = 0.0, strengthP = 0.0, strengthH = 0.0;
   double iFlipper = 0.0, slopMm = 0.0, cContactNsPerMm = 0.0, cStiffNsPerMm = 0.0;
   double kSpringNPerMm = 0.0, eosDeg = 0.0;
};

namespace FlipperModel
{
   constexpr double PI = 3.14159265358979323846;
   constexpr double DEG = PI / 180.0;

   inline double Sq(double x) { return x * x; }

   inline double SliderX(const FlipperSettings& s, const FlipperDerived& p, double theta)
   {
      const double b = p.offsetMm - s.rCrankMm * sin(theta);
      return s.rCrankMm * cos(theta) + sqrt(fmax(s.lRodMm * s.lRodMm - b * b, 1e-6));
   }
   inline double LeverMm(const FlipperSettings& s, const FlipperDerived& p, double theta)
   {
      const double sn = sin(theta), cs = cos(theta);
      const double b = p.offsetMm - s.rCrankMm * sn;
      const double root = sqrt(fmax(s.lRodMm * s.lRodMm - b * b, 1e-6));
      return fabs(-s.rCrankMm * sn + (b * s.rCrankMm * cs) / root);
   }
   // slider position and lever arm at the same crank angle: one sin, cos and sqrt for both (what SliderX and LeverMm compute separately)
   struct Kin { double x, lever; };
   inline Kin Kinematics(const FlipperSettings& s, const FlipperDerived& p, double theta)
   {
      const double sn = sin(theta), cs = cos(theta);
      const double b = p.offsetMm - s.rCrankMm * sn;
      const double root = sqrt(fmax(s.lRodMm * s.lRodMm - b * b, 1e-6));
      return Kin{ s.rCrankMm * cs + root, fabs(-s.rCrankMm * sn + (b * s.rCrankMm * cs) / root) };
   }
   inline void SolveAB(double larm, double lair, double x0, double& a, double& b)
   {
      b = (larm - lair) * x0 * (1.0 + x0);
      a = larm - b / x0;
   }

   // time the flipper takes to return from the end stop to rest under the spring alone (rigid one-body)
   inline double ReturnTimeS(const FlipperSettings& s, const FlipperDerived& p, double kSpring)
   {
      constexpr double dt = 1e-6;
      double phi = s.stopDeg * DEG, omega = 0.0, t = 0.0;
      while (phi > 0.0)
      {
         const double theta = p.theta0 + phi;
         const double leverM = LeverMm(s, p, theta) / 1000.0;
         const double force = s.preloadN + kSpring * (p.xRest - SliderX(s, p, theta));
         const double iEff = p.iFlipper + s.plungerKg * leverM * leverM;
         omega += (-force * leverM / iEff) * dt;
         phi += omega * dt;
         t += dt;
         if (t > 1.0)
            return t;
      }
      return t;
   }

   inline void Derive(const FlipperSettings& s, FlipperDerived& p)
   {
      p.iFlipper = s.iBat + s.iCrank;
      p.offsetMm = s.rCrankMm * cos(s.stopDeg * 0.5 * DEG);
      p.theta0 = (s.centerDeg - s.stopDeg * 0.5) * DEG;
      p.xMin = 1e30;
      p.xMax = -1e30;
      for (int i = 0; i <= 500; i++)
      {
         const double x = SliderX(s, p, p.theta0 + s.stopDeg * DEG * i / 500.0);
         p.xMin = fmin(p.xMin, x);
         p.xMax = fmax(p.xMax, x);
      }
      p.xRest = SliderX(s, p, p.theta0);
      SolveAB(s.lArmP, s.lAirP, s.x0P, p.aP, p.bP);
      SolveAB(s.lArmH, s.lAirH, s.x0H, p.aH, p.bH);
      const double issP = s.vSupply / s.rP, issH = s.vSupply / s.rH;
      p.strengthP = 0.5 * issP * issP * p.bP / (s.bobbinMm * s.x0P * s.x0P);
      p.strengthH = 0.5 * issH * issH * p.bH / (s.bobbinMm * s.x0H * s.x0H);
      p.eosDeg = s.stopDeg - s.eosOffsetDeg;
      const double leverAtStop = LeverMm(s, p, p.theta0 + s.stopDeg * DEG);
      p.slopMm = leverAtStop * s.slopDeg * DEG;
      const double mReflected = p.iFlipper / Sq(leverAtStop / 1000.0);   // flipper inertia reflected to the plunger axis, kg
      p.cContactNsPerMm = 2.0 * s.contactZeta * sqrt(s.kContactNPerMm * 1000.0 * mReflected) / 1000.0;
      p.cStiffNsPerMm = 2.0 * s.stiffZeta * sqrt(s.kContactNPerMm * 1000.0 * mReflected) / 1000.0;
      if (s.kSpringNPerMm > 0.0)
         p.kSpringNPerMm = s.kSpringNPerMm;
      else
      {
         double lo = 0.001, hi = 50.0;   // return time falls as K rises
         for (int i = 0; i < 50; i++)
         {
            const double mid = 0.5 * (lo + hi);
            if (ReturnTimeS(s, p, mid) > s.returnMs / 1000.0)
               lo = mid;
            else
               hi = mid;
         }
         p.kSpringNPerMm = 0.5 * (lo + hi);
      }
   }

   // force the link exerts ON THE PLUNGER (N); reaction torque on the flipper = this * lever.
   // Contacts are one-sided (never pull); damping acts in the stiff contacts and, if enabled, in free play.
   inline double LinkForceN(const FlipperSettings& s, const FlipperDerived& p, double mismatchMm, double rateMmPerS)
   {
      if (mismatchMm < 0.0)
         return fmax(0.0, -(s.kContactNPerMm * mismatchMm) - p.cStiffNsPerMm * rateMmPerS);                       // loaded side
      if (mismatchMm > p.slopMm)
         return fmin(0.0, -(s.kContactNPerMm * (mismatchMm - p.slopMm)) - p.cStiffNsPerMm * rateMmPerS);          // far side
      return s.freePlayDamping ? -(p.cContactNsPerMm * rateMmPerS) : 0.0;                                         // free play
   }

   // rest stop: the crank arm on its rubber stop, a one-sided spring and damper that act from phi = 0 (pen = -phi). Torque on the flipper
   // (closing direction) = K*pen + c*d(pen)/dt, never pulling. All energy loss is velocity-dependent, so a slow push stores its energy and
   // gives it all back (no permanent set), while a fast hit loses energy in the damper. c = 2*zeta*sqrt(K*I), I = the inertia the stop acts on.
   inline double RestStopTorque(const FlipperSettings& s, double phi, double omega, double inertia)
   {
      const double pen = fmax(0.0, -phi);
      if (pen <= 0.0)
         return 0.0;
      // progressive stop: torque ~ pen^n, K is the torque slope at the reference penetration (1.4 deg); n = 1 is the plain linear stop
      const double penRef = 1.4 * M_PI / 180.0;
      const double r = pen / penRef;
      const double c = 2.0 * s.restZeta * sqrt(s.restKNmPerRad * inertia) * pow(r, 0.5 * (s.restN - 1.0));
      return fmax(0.0, s.restKNmPerRad * penRef * pow(r, s.restN) - c * omega);
   }

   inline double SpringTorque(const FlipperSettings& s, const FlipperDerived& p, double theta)
   {
      const double force = s.preloadN + p.kSpringNPerMm * (p.xRest - SliderX(s, p, theta));
      return force * LeverMm(s, p, theta) / 1000.0;
   }
   inline double SpringTorqueK(const FlipperSettings& s, const FlipperDerived& p, const Kin& k)
   {
      const double force = s.preloadN + p.kSpringNPerMm * (p.xRest - k.x);
      return force * k.lever / 1000.0;
   }
}

// Applies one [Flippers] key (l_rod_mm, r_crank_mm, wiring, ...) to a settings set; false if the key or its value is not understood.
// Defined in hitflipper.cpp. The sling uses it for its sling_<key> lines.
bool FlipperApplyKey(FlipperSettings &s, const std::string &key, const std::string &val);
