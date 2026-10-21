// license:GPLv3+

#include "core/stdafx.h"

#include "physics/hitsling.h"

SlingSettings::SlingSettings()
{
   // the model flipper's mechanics and return spring rate, with the sling's lever arm, travel and coil. All unverified until measured.
   FlipperSettings &g = mech;
   g.series = false;
   g.single = true;               // one winding
   g.freePlayDamping = true;
   g.lRodMm = 19.05;              // the model flipper's link
   g.rCrankMm = 50.8;             // 2 in lever arm
   g.centerDeg = 90.0;
   g.stopDeg = 20.0;              // plate travel
   g.iBat = 1.383e-5;             // 2 x 16.1 g steel arms (2 x 1/2 x 1/8 in), each a rod about its end
   g.iCrank = 1.383e-5;
   g.plungerKg = 0.0413355;       // the model flipper's plunger
   g.bDamping = 0.001;            // pivot and link friction, N m s/rad (a guess: it settles the plate after a hit)
   g.vSupply = 28.0;
   g.bobbinMm = 35.3;             // AE-23-800 (MAYA data; 35.3 mm is the figure for that group of playfield coils)
   g.rP = 3.769;
   g.lAirP = 5.26;
   g.lArmP = 15.77;
   g.x0P = 0.3;                   // tried 0.1 - 1.0: 0.3 - 0.6 gives the strongest, most even kick (the model flipper's power coil has 0.2)
   g.preloadN = 0.3;              // the coil is vertical under the plate with the plunger above it, so the plunger's weight (0.405 N) offsets most of the spring; this small preload holds the plate on its rest stop
   g.kSpringNPerMm = 0.0619;
   g.restKNmPerRad = 15.0;        // rubber rest stop
   g.restZeta = 0.45;
   g.restN = 1.0;
   g.slopDeg = 0.0;               // no slop in the link
   g.kContactNPerMm = 20.0;       // the coil stop is the stiff contact
   g.contactZeta = 0.1;
   g.stiffZeta = 0.12;
   g.rigidMaxRelSpeedMmS = 50.0;
}

void SlingPlate::Configure(const SlingSettings &s, const bool reset)
{
   set = s;
   set.mech.series = false;
   set.mech.single = true;   // the sling coil is always one winding
   FlipperModel::Derive(set.mech, der);
   if (reset)
   {
      st = State();
      st.xp = FlipperModel::SliderX(set.mech, der, der.theta0);   // at rest, the link taut
      alpha = 0.0;
      hitStop = false;
   }
}

// One SI substep of h seconds, explicit Euler with the state at the start of the step used throughout. The same equations as FlipperMoverObject::FlipperStep
// (hitflipper.cpp) for a single winding with no EOS switch; the coil is on while the driver pulses.
void SlingPlate::Step(const double h, const bool coilOn, const double extTorqueNm)
{
   using namespace FlipperModel;
   const FlipperSettings &s = set.mech;
   const FlipperDerived &p = der;
   const double hMs = h * 1000.0;

   const Kin kin = Kinematics(s, p, p.theta0 + st.phi);
   const double leverMm = kin.lever;
   const double leverM = leverMm / 1000.0;

   const double frac = (st.xp - p.xMin) / s.bobbinMm;
   const double gHigh = Sq(s.x0P / (frac + s.x0P));
   const double vPos = fmax(-(st.vp / 1000.0), 0.0);                 // back-EMF depends on the plunger's closing speed
   const double dLdxP = p.bP / (s.bobbinMm * Sq(frac + s.x0P));      // mH/mm
   const double lP = p.aP + p.bP / (frac + s.x0P);                   // mH

   // circuit: di/dt = (V - R i - i dL/dx v) / L, mH and ms give amps; i is kept as a fraction of V/R
   const double issP = s.vSupply / s.rP;
   const double iP = st.iHigh * issP;
   const double vNow = coilOn ? s.vSupply : 0.0;
   const double iNext = fmax(0.0, iP + ((vNow - s.rP * iP - iP * dLdxP * vPos) / lP) * hMs);
   const double iHighNext = iNext / issP;

   const double fCoil = p.strengthP * gHigh * st.iHigh * st.iHigh;

   // hand-off between the two-body model and the rigid one-body
   const double mismatch = st.xp - kin.x;
   const double rate = st.vp + leverMm * st.omega;
   if (!st.rigid && !coilOn && mismatch < 0.0 && fabs(rate) < s.rigidMaxRelSpeedMmS)
      st.rigid = true;
   else if (st.rigid && coilOn)
      st.rigid = false;

   hitStop = false;

   if (st.rigid)
   {
      // one rigid body: the coil's generalized force on phi is fCoil*lever; the link carries no separate force and the plunger is slaved to the plate angle
      const double iEff = p.iFlipper + s.plungerKg * leverM * leverM;
      const double restT = RestStopTorque(s, st.phi, st.omega, iEff);
      const double net = fCoil * leverM - SpringTorqueK(s, p, kin) - s.bDamping * st.omega + restT + extTorqueNm;
      const double a = net / iEff;
      double omegaNext = st.omega + a * h;
      double phiNext = st.phi + st.omega * h;
      const Kin kn = Kinematics(s, p, p.theta0 + phiNext);
      st.xp = kn.x;
      st.vp = -kn.lever * omegaNext;
      if (phiNext > s.stopDeg * DEG)   // steel end stop (inelastic)
      {
         if (st.omega > 0.0)
            hitStop = true;
         phiNext = s.stopDeg * DEG;
         omegaNext = fmin(omegaNext, 0.0);
      }
      alpha = a;
      st.omega = omegaNext;
      st.phi = phiNext;
      st.iHigh = iHighNext;
      return;
   }

   const double fLink = LinkForceN(s, p, mismatch, rate);

   // plunger: the coil pulls toward closure (-x); hard stops at both ends
   const double fPlunger = -fCoil + fLink;
   const bool atClosed = st.xp <= p.xMin;
   const bool atOpen = st.xp >= p.xMax;
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
      const double a = fPlunger / s.plungerKg;   // m/s^2
      vNext = st.vp + a * 1000.0 * h;
      xNext = st.xp + st.vp * h;
      if (!atClosed && xNext <= p.xMin)
      {
         if (st.vp < 0.0)
            hitStop = true;
         xNext = p.xMin;
         vNext = 0.0;
      }
      else if (!atOpen && xNext >= p.xMax)
      {
         xNext = p.xMax;
         vNext = 0.0;
      }
   }

   // the plate: the return spring acts directly and continuously; it is constrained only through the link
   const double restT = RestStopTorque(s, st.phi, st.omega, p.iFlipper);
   const double netTorque = -SpringTorqueK(s, p, kin) - s.bDamping * st.omega + restT + fLink * leverM + extTorqueNm;
   alpha = netTorque / p.iFlipper;
   const double omegaNext = st.omega + alpha * h;
   st.phi += st.omega * h;   // explicit Euler on the state at the start of the step: phi uses the old omega
   st.omega = omegaNext;
   st.xp = xNext;
   st.vp = vNext;
   st.iHigh = iHighNext;
}
