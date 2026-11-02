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
   // ball / rubber contact: the shared rubber law in the group step. All values are unverified defaults meant to be edited; see MakeRubberParams / ClampRubberSettings.
   bool rubberContact = true;           // rubber_contact: read so a table that sets it to 0 is told it is ignored (a bat has no impulse path, MakePart)
   double rubberEMpa = 12.0;            // effective Young's modulus
   double rubberThicknessMm = 3.2;      // rubber layer thickness
   double rubberAlphaSPerM = 0.05;      // normal damping; rebound ratio ~ 1 - alpha*v
   double rubberBetaLayer = 0.33;       // thin-layer stiffening coefficient
   double rubberMu = 0.8;               // rubber / steel friction (black / white rubber; user, 2026-10-30)
   double rubberShearScale = 0.33;      // rubber_sideways_stiffness: multiplier on Mindlin's 8 G* a. rubber_e_mpa 12 is the normal stiffness of a thin confined strip; shear is not stiffened by that confinement, so G is the plain rubber's (E about 4, as the rings): 4 / 12. The strip lies on the bat, so nothing softer is in series (a backed layer is at least as stiff as a half-space). Cross-check G pi (a + h/2)^2 / h, a = 2.3 - 3.7 mm: 0.33 - 0.46. Calculated, not measured (runs/2026-10-30-rubber-review/notes.md)
   bool rubberTwist = true;             // resistance to spin about the contact normal
   double rubberTanDelta = 0.05;        // rubber loss tangent for the rolling resistance, 0 = off. PROVISIONAL: starts low on purpose (typical elastomer 0.1-0.3), to be tuned
   double rubberZLoMm = 9.144;          // rubber bottom above the playfield surface (0.36 in)
   double rubberHeightMm = 11.176;      // rubber height (0.44 in)
   double rubberReleaseGapMm = 0.2;     // ...and be farther than this from the rubber surface (release hysteresis for a resting ball)
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
   struct Kin { double x, lever, dLever; };   // dLever: d(lever) / d(theta), mm / rad
   inline Kin Kinematics(const FlipperSettings& s, const FlipperDerived& p, double theta)
   {
      const double sn = sin(theta), cs = cos(theta);
      const double b = p.offsetMm - s.rCrankMm * sn;
      const double root = sqrt(fmax(s.lRodMm * s.lRodMm - b * b, 1e-6));
      const double r = s.rCrankMm, db = -r * cs;
      const double g = -r * sn + (b * r * cs) / root;
      const double dg = -r * cs + (db * r * cs - b * r * sn) / root + (b * r * cs) * (b * db) / (root * root * root);
      return Kin{ r * cs + root, fabs(g), g >= 0.0 ? dg : -dg };
   }
   // One rigid body with the plunger slaved through the link: its inertia depends on the angle, so I w' + 1/2 I' w^2 = torque (Lagrange). This is the 1/2 I' w^2 term, N m; without it the
   // bat's energy changes with no work done.
   inline double SlavedInertiaTerm(const FlipperSettings& s, const Kin& k, double omega)
   {
      return s.plungerKg * (k.lever / 1000.0) * (k.dLever / 1000.0) * omega * omega;
   }
   // The windings over a step of h: the circuit's exponential with L and R + dL/dx v held (the plunger's place at mid-step, by its speed now). From the state at the step's start, so a step taken again is the same.
   // iHigh, iLow: the currents at the end, fractions of each winding's V / R; i2High, i2Low: their squares' exact means over the step, same units (the force goes as i^2, and its fall is twice as fast as the current's).
   struct CoilStep { double iHigh = 0.0, iLow = 0.0, i2High = 0.0, i2Low = 0.0; };
   inline CoilStep CoilAdvance(const FlipperSettings& s, const FlipperDerived& p, const double xM, const double vX, const double iHigh, const double iLow, const bool coilOn, const bool preEos, const bool coilOffPrev, const bool wasPreEos, const double h)
   {
      const bool coilOff = !coilOn;
      double frac = (xM + 0.5 * h * vX - p.xMin * 1e-3) / (s.bobbinMm * 1e-3);
      frac = fmax(frac, 0.0);
      const double vClose = fmax(-vX, 0.0);
      const double dP = p.bP / (s.bobbinMm * Sq(frac + s.x0P)), dH = p.bH / (s.bobbinMm * Sq(frac + s.x0H));
      const double lP = (p.aP + p.bP / (frac + s.x0P)) * 1e-3, lH = (p.aH + p.bH / (frac + s.x0H)) * 1e-3;
      const double issP = s.vSupply / s.rP, issH = s.vSupply / s.rH;
      double iP = iHigh * issP, iH = iLow * issH;
      const auto Exp = [&](const double i0, const double V, const double R, const double L, double& meanI2) {
         const double iInf = V / R, d = i0 - iInf, u = h * R / L;
         const double e1 = exp(-u), m1 = u > 1e-6 ? (1.0 - e1) / u : 1.0 - 0.5 * u, m2 = u > 1e-6 ? (1.0 - e1 * e1) / (2.0 * u) : 1.0 - u;
         meanI2 = iInf * iInf + 2.0 * iInf * d * m1 + d * d * m2;
         return iInf + d * e1;
      };
      double mP = 0.0, mH = 0.0, iP1, iH1;
      if (s.single)
      {
         iP1 = Exp(iP, coilOff ? 0.0 : s.vSupply, s.rP + dP * vClose, lP, mP);
         iH1 = 0.0;
      }
      else if (s.series)
      {
         if (coilOff)
         {
            if (!coilOffPrev)
            {
               const double iLoop = (preEos && s.eosFlux) ? (lP * iP + lH * iH) / (lP + lH) : iP;
               iP = iH = iLoop;
            }
            else
               iH = iP;
            iP1 = iH1 = Exp(iP, 0.0, s.rP + s.rH + (dP + dH) * vClose, lP + lH, mP);
            mH = mP;
         }
         else if (preEos)
         {
            iP1 = Exp(iP, s.vSupply, s.rP + dP * vClose, lP, mP);
            iH1 = Exp(iH, 0.0, s.rH, lH, mH);
         }
         else
         {
            if (wasPreEos)
            {
               const double iLoop = s.eosFlux ? (lP * iP + lH * iH) / (lP + lH) : iP;
               iP = iH = iLoop;
            }
            else
               iH = iP;
            iP1 = iH1 = Exp(iP, s.vSupply, s.rP + s.rH + (dP + dH) * vClose, lP + lH, mP);
            mH = mP;
         }
      }
      else
      {
         const double V = coilOff ? 0.0 : 1.0;   // parallel
         iP1 = Exp(iP, preEos ? V * s.vSupply : 0.0, s.rP + (preEos ? dP * vClose : 0.0), lP, mP);
         iH1 = Exp(iH, V * s.vSupply, s.rH + dH * vClose, lH, mH);
      }
      CoilStep o;
      o.iHigh = iP1 / issP;
      o.iLow = iH1 / issH;
      o.i2High = mP / (issP * issP);
      o.i2Low = mH / (issH * issH);
      return o;
   }
   // The time constant the windings are changing with now, s; infinite when they are at their targets (for a step's facts)
   inline double CoilTimeConstant(const FlipperSettings& s, const FlipperDerived& p, const double xM, const double vX, const double iHigh, const double iLow, const bool coilOn, const bool preEos)
   {
      const double frac = (xM - p.xMin * 1e-3) / (s.bobbinMm * 1e-3);
      const double vClose = fmax(-vX, 0.0);
      const double lP = (p.aP + p.bP / (frac + s.x0P)) * 1e-3, lH = (p.aH + p.bH / (frac + s.x0H)) * 1e-3;
      const double dP = p.bP / (s.bobbinMm * Sq(frac + s.x0P)), dH = p.bH / (s.bobbinMm * Sq(frac + s.x0H));
      const double issP = s.vSupply / s.rP, issH = s.vSupply / s.rH;
      double tau = INFINITY;
      if (s.single)
      {
         if (fabs((coilOn ? 1.0 : 0.0) - iHigh) > 1e-3)
            tau = lP / (s.rP + dP * vClose);
      }
      else if (!coilOn || !preEos)
      {
         const double tgt = coilOn ? s.vSupply / (s.rP + s.rH) : 0.0;
         if (fabs(iHigh * issP - tgt) > 1e-3 * issP * 0.1 || fabs(iLow * issH - tgt) > 1e-3 * issP * 0.1)
            tau = (lP + lH) / (s.rP + s.rH + (dP + dH) * vClose);
      }
      else
      {
         if (fabs(1.0 - iHigh) > 1e-3)
            tau = fmin(tau, lP / (s.rP + dP * vClose));
         if (iLow > 1e-3)
            tau = fmin(tau, lH / s.rH);
      }
      return tau;
   }
   // A given motion: one pulse along an axis, half-cosine rise then half-cosine return (the other bat's strike moving a pivot)
   struct PivotPulse
   {
      bool active = false;
      double t = 0.0, amp = 0.0, rise = 0.001, fall = 0.008;
      void Eval(const double at, double& disp, double& vel) const
      {
         if (!active || at < 0.0 || at >= rise + fall)
            disp = vel = 0.0;
         else if (at < rise)
         {
            disp = amp * 0.5 * (1.0 - cos(PI * at / rise));
            vel = amp * 0.5 * (PI / rise) * sin(PI * at / rise);
         }
         else
         {
            const double u = at - rise;
            disp = amp * 0.5 * (1.0 + cos(PI * u / fall));
            vel = -amp * 0.5 * (PI / fall) * sin(PI * u / fall);
         }
      }
      double Move(const double h) const   // its move over the next h
      {
         if (!active)
            return 0.0;
         double d0, d1, v;
         Eval(t, d0, v);
         Eval(t + h, d1, v);
         return d1 - d0;
      }
      void Start(const double a, const double riseS, const double fallS)
      {
         if (active)   // one still running: let it finish
            return;
         amp = a;
         rise = riseS;
         fall = fallS;
         t = 0.0;
         active = true;
      }
      void Advance(const double h)
      {
         t += h;
         if (t >= rise + fall)
            active = false;
      }
   };
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
         const Kin k = Kinematics(s, p, theta);
         const double leverM = k.lever / 1000.0;
         const double force = s.preloadN + kSpring * (p.xRest - k.x);
         const double iEff = p.iFlipper + s.plungerKg * leverM * leverM;
         omega += ((-force * leverM - SlavedInertiaTerm(s, k, omega)) / iEff) * dt;
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

   inline double SpringTorqueK(const FlipperSettings& s, const FlipperDerived& p, const Kin& k)
   {
      const double force = s.preloadN + p.kSpringNPerMm * (p.xRest - k.x);
      return force * k.lever / 1000.0;
   }
}

// A single-winding coil with a plunger: the circuit and force laws the sling plate uses (hitsling.cpp), as shared functions with the coil parameters as input.
// The plunger position is the fraction frac = (distance from the closed end) / bobbin length; inductance L = a + b / (frac + x0) (mH), force = strength * g * (i / iss)^2
// with g = (x0 / (frac + x0))^2 and iss = V / R. Pure functions; the owner holds the current.
struct CoilSettings
{
   double vSupply = 28.0, rOhm = 3.769, lAirMH = 5.26, lArmMH = 15.77, bobbinMm = 35.3, x0 = 0.3;   // AE-23-800
};
struct CoilDerived
{
   double a = 0.0, b = 0.0, iss = 0.0, strength = 0.0;   // mH, mH, A, N
};
namespace CoilModel
{
   inline void Derive(const CoilSettings &s, CoilDerived &d)
   {
      FlipperModel::SolveAB(s.lArmMH, s.lAirMH, s.x0, d.a, d.b);
      d.iss = s.vSupply / s.rOhm;
      d.strength = 0.5 * d.iss * d.iss * d.b / (s.bobbinMm * s.x0 * s.x0);
   }
   // force on the plunger toward the closed end (N) at current i (A)
   inline double Force(const CoilSettings &s, const CoilDerived &d, double frac, double i)
   {
      const double g = FlipperModel::Sq(s.x0 / (frac + s.x0)), ih = i / d.iss;
      return d.strength * g * ih * ih;
   }
   // The circuit and the force over a step of h in which the plunger moves dx toward the closed end, from frac0 to frac1, under volts vNow: a function of the start and the move only, so a retake is identical and the solve can hold it.
   // L di/dt + (R + C v) i = V over the step, with L the mean of its ends, C the chord of L along the move and v = dx / h (signed: a plunger drawn out gives its work back to the circuit), solved exactly. The force is what closes the energy exactly:
   // supplied - lost in R - change in 1/2 L i^2 = F dx at any step length, F = C (mean i^2 - (i0^2 + i1^2) / 4), which is C i^2 / 2 to second order. The coil's only energy law: nothing else books it.
   struct OverMove { double i1 = 0.0, meanI = 0.0, meanI2 = 0.0, force = 0.0; };
   inline OverMove Over(const CoilSettings &s, const CoilDerived &d, const double frac0, const double frac1, const double dx, const double vNow, const double i0, const double h)
   {
      const double a0 = frac0 + s.x0, a1 = frac1 + s.x0;
      const double C = d.b / (s.bobbinMm * a0 * a1);   // (mH/mm = H/m) (L1 - L0) / dx, in closed form: no difference of near-equal numbers
      const double lH = (d.a + 0.5 * d.b * (1.0 / a0 + 1.0 / a1)) * 1e-3;
      double R = s.rOhm + C * dx / h;
      if (R == 0.0)
         R = 1e-12;   // (the line below divides; a drawn-out plunger so fast it cancels R exactly)
      const double iInf = vNow / R, di = i0 - iInf, u = h * R / lH;
      const double e1 = std::exp(-u), au = std::fabs(u), m1 = au > 1e-6 ? (1.0 - e1) / u : 1.0 - 0.5 * u, m2 = au > 1e-6 ? (1.0 - e1 * e1) / (2.0 * u) : 1.0 - u;
      OverMove o;
      o.i1 = iInf + di * e1;
      o.meanI = iInf + di * m1;
      o.meanI2 = iInf * iInf + 2.0 * iInf * di * m1 + di * di * m2;
      o.force = C * (o.meanI2 - 0.25 * (i0 * i0 + o.i1 * o.i1));
      return o;
   }
}

// Applies one [Flippers] key (l_rod_mm, r_crank_mm, wiring, ...) to a settings set; false if the key or its value is not understood.
// Defined in hitflipper.cpp. The sling uses it for its sling_<key> lines.
bool FlipperApplyKey(FlipperSettings &s, const std::string &key, const std::string &val);
