// license:GPLv3+
//
// The mechanical plunger's rod (docs/unified-physics/PlungerDesign.txt section 1; piece 1a: the rod alone, no ball, no hand). One body on one coordinate (bodypart.h), a slide along
// the lane: s = the pull, 0 at rest, + toward the player (pulled), - toward the ball (pushed in). Rod, knob, tip, washer and clip move as one, with a third of the main spring's mass.
//   main spring   linear, F = k_m (x0 + s) toward the ball (-s); k_m = G d^4 / (8 (OD - d)^3 n), x0 = free length - gap (its preload compression at rest)
//   outer spring  between the knob and the housing, pushes the knob back (+s) only while touching: F = k_o (c0 - s) for s < c0, 0 beyond; c0 = k_m x0 / k_o (rest balance)
//   gravity       the part's gravity along the lane: m g.axis toward +s (the table's slope is in g)
//   stops         full pull: the main spring solid at s = travel (side +1); pushed in: the outer spring solid at s = -(outer free - outer solid - c0) (side -1).
//                 Steel on steel BodyStops (the impact law, the rest rule), as the gates' stops
//   losses        viscous drag F = -c v (the rod and spring rubbing and rattling in the sleeve); Coulomb friction in the bore (piece 1b, design 1f(i)): the bore's force across the
//                 axis N = |(m g + F_ball) less its part along the axis|, slip limit mu N, the stick spring Mindlin's 8 G* a of one steel bearing patch (a = cbrt(3 N r / (4 E*)), r the
//                 rod's radius; the hinge's law with lever 1). One patch: a simplification (the real rod bears on two ends of its bore)
//   tip           (piece 1b, design 3a/3c) the rubber tip as a Revolve about the rod's axis: a hemisphere of the rod's radius on a cylinder of the same radius tipCylMm back from its
//                 centre, then closed (ButtonCap::RoundedOutline). Apex at s = 0 is the point tipApex; the surface moves along +axis for +q. Ring rubber (RubberObjectSettings defaults)
//   hand          (piece 2, design 2a/2b/2g) two stops the player's hand moves (bodypart.h driven stops), off until used: HOLD (side -1: the rod cannot move toward the
//                 ball past it; the hand pulls) and PUSH (side +1: the rod cannot move away from the ball past it). A key uses HOLD alone; the cabinet sensor ties
//                 both together at its reading (Tie). Each placed at the rod when the hand takes hold (BodyPart::SetStop), then driven at the hand's rate; their work
//                 is the hand's drive. No launch rule: let go, the spring drives.
//                 A soft gripping hand, not steel (user 2026-11-18): restitution 0, the knob joins the hand on contact (a perfectly inelastic catch: the rod takes
//                 the stop's speed, the kinetic energy difference booked as stop loss; no bounce). The Hertz k at the rod's radius stays only for the rest test
// Stored energy V(s) = k_m s (x0 + s / 2) + k_o max(c0 - s, 0)^2 / 2 - m (g.axis) s; its chord across the outer spring's kink at c0 is exact (piecewise), so the book closes at any step.
// Pure SI inside; the settings in mm, g, N/mm as the hardware data is given (PartData.txt section 11). No VPX types.

#pragma once

#include <algorithm>
#include <cmath>

#include "physics/bodypart.h"
#include "physics/buttoncap.h"

namespace PlungerRod
{
   using GroupStep::V3;

   // Every value: MEASURED (user's hardware), DERIVED (lead, from measured), ESTIMATED (no data); source in PartData.txt section 11 / PlungerDesign.txt section 1
   struct Settings
   {
      // main (inner) spring: the default brown spring in the default gap
      double mainWireMm = 1.0;          // wire diameter d, mm: MEASURED (user 2026-11-17, brown spring)
      double mainOdMm = 12.0;           // outside diameter, mm: MEASURED (user 2026-11-17)
      double mainActiveCoils = 24.5;    // active coils n: MEASURED (user 2026-11-17: 24 - 25 counted, closed ends not counted)
      double mainClosedPerEnd = 4.0;    // closed (inactive) coils at each end: MEASURED (user 2026-11-17, brown spring)
      double mainFreeMm = 152.4;        // free length, mm: MEASURED (user 2026-11-17, 6 in)
      double gapMm = 93.3;              // the length the rod gives the main spring at rest, mm: MEASURED (user 2026-11-17, 3.675 in without the simulator's spacer and washer)
      double shearModulusGPa = 79.3;    // music wire G, GPa: handbook value
      double steelDensity = 7850.0;     // kg/m3, the spring's wire (its moving third): handbook value
      // outer (barrel) spring between the knob and the housing
      double outerRateNmm = 9.4;        // k_o, N/mm: MEASURED (user 2026-11-17: about 20 lbf at 0.375 in; linear ASSUMED)
      double outerFreeMm = 19.05;       // free length, mm: MEASURED (0.75 in)
      double outerSolidMm = 8.5;        // solid length, mm: DERIVED (5 coils x 1.7 mm wire, MEASURED geometry)
      double outerWireMm = 1.7;         // wire diameter, mm: MEASURED (user 2026-11-17); its radius is the solid stop's Hertz radius
      // moving mass, g (the main spring's own third is added from its wire)
      double rodG = 110.0;              // DERIVED (3/8 in x 7-3/4 in steel rod)
      double knobG = 10.0;              // MEASURED typical (ball knobs 15 - 25)
      double tipG = 3.0;                // ESTIMATED
      double washerClipG = 2.0;         // ESTIMATED
      // stops: steel on steel coil bind, restitution stopE at stopVRef m/s (the gates' steel stop value): ESTIMATED (coil bind not measured)
      double stopE = 0.5;
      double stopVRef = 0.1;
      // viscous drag c, N s/m: DERIVED (runs/2026-11-18-plunger-bounce/results_damping.txt: with rubN, the simulator plunger's empty full release in the user's
      // 960 fps video, 2026-11-18: about 40 % of the spring energy lost by rest, 4 bounces back past 0.5 mm, still at about 220 ms; rounded, not a close fit)
      double dragNsPerM = 1.0;
      // the rod and main spring rubbing in the sleeve: a constant Coulomb force, N, whatever the side load. DERIVED from the same video (2.0 N fits it; the bounce
      // ratio falling with amplitude is Coulomb's signature) less an ESTIMATED share for the simulator's optical slide sensor (user 2026-11-18): 1.5 N
      double rubN = 1.5;
      // Coulomb friction in the bore on the side load (piece 1b): mu ESTIMATED; the rod's radius for the bearing patch (3/8 in rod, DERIVED)
      double boreMu = 0.2;
      double rodRadiusMm = 4.7625;
      // the tip (design 3a): hemisphere radius MEASURED (3/8 in across); the cylinder behind it ESTIMATED (only has to exceed any squeeze); rubber from the apex to the rod end
      // for the thin-layer term ESTIMATED (0 = a thick body: no layer term)
      double tipRadiusMm = 4.7625;
      double tipCylMm = 6.0;
      double tipRubberMm = 4.7625;
      // the tip's rubber: the rings' values (contacthost.h RubberObjectSettings defaults: rubber_e_mpa 12, alpha 0.18, beta_layer 0.33, mu 0.8, sideways_stiffness 0.1, tan_delta 0.05)
      double tipEMpa = 12.0, tipAlpha = 0.18, tipBetaLayer = 0.33, tipMu = 0.8, tipShearScale = 0.1, tipTanDelta = 0.05;
      // the hand's pull speed on a key (key down to full pull, PullBackandRetract's return), mm/s: ESTIMATED (design 2b: a hand pulling the 61 mm in about half a second)
      double pullMmS = 120.0;
      // the cabinet sensor's release (design 2k MERGED, TODO 14.3): a release when a quadratic fit of the pulled reports gives the knob an acceleration toward rest over
      // release_accel, m/s^2 (or a fitted speed toward rest over 0.5 m/s). MEASURED boundary (the user's recording 20261118_231530, 21 ms fit, peak per move: held and
      // guided returns 3-9 m/s^2; free releases 15-83, the short ones of 212416 10-25): above the held and guided, at the bottom of the free releases (the two
      // shortest, 10 and 14, stay tied: staying tied is acceptable, a false release is not, user 2026-11-18)
      double releaseAccel = 15.0;
      // ... and a fitted speed toward rest over release_accel_min_speed, m/s: above a hand easing back from a hold (e2k run U: a guided return starting from a hold
      // fitted 5.7 m/s^2 at 0.081 m/s, a false release at 0.05). MEASURED: held and guided at most 0.39 m/s, free releases 0.27-2.6 m/s (peak per move); free
      // releases at detection 0.15 and up
      double releaseAccelMinSpeed = 0.15;
      // the fit's window W (s) from the board's step q (m, the smallest change seen between reports): W = sqrt(2 C q / release_window_noise), so a steady move's
      // quantization noise C q / W^2 stays under release_window_noise / 2; at least 4 reports. C = release_window_c, from the lead's simulation (20000 steady returns
      // 0.03-0.2 m/s, +-1 ms jitter, 3 ms reports: worst fake acceleration 4.7 m/s^2 at q 0.1 mm, W 15 ms: C = 4.7 x 0.015^2 / 1e-4). The PinOne (q 0.118 mm): W 22 ms
      double releaseWindowC = 10.6;
      // the window's noise allowance, m/s^2, independent of release_accel (TODO 14.3, user 2026-11-19): W tied to release_accel shrank to 13 ms when the threshold
      // went to 15 and held / guided returns then read 15-23 m/s^2 (19 false releases in the user's regression set); at 5 (W 22 ms) none. MEASURED on the PinOne
      double releaseWindowNoise = 5.0;
      // the window's floor, ms (user 2026-11-19): over a short window a hand's own quick speed changes read as acceleration (the user's guided moves read 15-23 m/s^2
      // over 13 ms where the PinOne's step noise is only about 7), so a finer board must not shorten it; the board's step only lengthens it. MEASURED on the user's set
      double releaseWindowMinMs = 22.0;
      // the cabinet sensor's source (design 2i a, user 2026-11-18): false (default) = the sensor's own reports (the position mapping's value, a new one about every 6 ms);
      // true = VPX's filtered position and speed (PlungerSensor's Kalman filter: lags 10-20 ms, overshoots 2-4 mm, MEASURED), for a noisy or irregular sensor
      bool sensorFiltered = false;
   };

   // its own state beside the body's (PartRecords)
   struct RodState
   {
      V3 g = V3(0.0, 0.0, -9.81);   // the gravity the part feels (the owner's, between cycles; kept so a replay runs on the same)
   };

   class Part : public GroupStep::BodyPart, public RodState
   {
   public:
      Part() { KeepState<RodState>(this); }

      V3 axis = V3(0.0, 1.0, 0.0);   // the unit direction of +s (toward the player) in the world
      double km = 0.0, x0 = 0.0, ko = 0.0, c0 = 0.0, c = 0.0;   // N/m, m, N/m, m, N s/m
      double solid = 0.0, travel = 0.0, outerSolidS = 0.0;      // m: main spring solid length, full pull, the outer spring's solid point (s < 0)
      double springMass = 0.0;                                  // kg: the main spring's whole mass (a third of it moves with the rod)
      int stopPull = -1, stopIn = -1, stopHold = -1, stopPush = -1;
      double pullRate = 0.0;                                    // m/s: the key's pull speed (Settings::pullMmS)
      double boreMu = 0.0, rRod = 0.0;                          // bore friction, the rod's radius (m)
      double rub = 0.0;                                         // N: the sleeve rub (Settings::rubN), added to the bore's Coulomb limit
      V3 tipApex;                                               // the tip's apex at s = 0 in the world, m
      JetBumper::Revolve prof;                                  // the tip in (rho from the axis, z along the ball-ward axis -axis, apex at z 0)
      RubberContact::RubberParams rp;                           // the tip's ball law
      double layerM = 0.0;                                      // the rubber's thickness for the layer term, m (the group's Shape::layerM)
      static constexpr double kSteelEStar = 115e9, kSteelGStar = 23.8e9;   // steel on steel, as the hinge's pin (Pa)
      static constexpr double kEdgeMinM = 0.2e-3;               // the smallest rounding a probe counts (the Button's)

      // axisDir: +s in the world (unit); gravity: the part's gravity now; apex: the tip's apex at s = 0 (the owner's placement)
      void Build(const Settings &s, const V3 &axisDir, const V3 &gravity, const V3 &apex = V3())
      {
         const double mm = 1e-3, d = s.mainWireMm * mm, D = (s.mainOdMm - s.mainWireMm) * mm, pi = GroupStep::kPi;
         axis = axisDir;
         g = gravity;
         km = s.shearModulusGPa * 1e9 * d * d * d * d / (8.0 * D * D * D * s.mainActiveCoils);
         x0 = (s.mainFreeMm - s.gapMm) * mm;
         solid = (s.mainActiveCoils + 2.0 * s.mainClosedPerEnd) * d;
         travel = s.gapMm * mm - solid;
         ko = s.outerRateNmm * 1e3;
         c0 = km * x0 / ko;
         outerSolidS = -((s.outerFreeMm - s.outerSolidMm) * mm - c0);
         springMass = (s.mainActiveCoils + 2.0 * s.mainClosedPerEnd) * pi * D * (0.25 * pi * d * d) * s.steelDensity;
         mass = (s.rodG + s.knobG + s.tipG + s.washerClipG) * 1e-3 + springMass / 3.0;
         c = s.dragNsPerM;
         rub = s.rubN;
         boreMu = s.boreMu;
         rRod = s.rodRadiusMm * mm;
         tipApex = apex;
         {
            // (0, 0) the apex, (rt, 0), (rt, -(2 rt + cyl)): the face rounded rt is the hemisphere centred at (0, -rt); the side's end, rounded as the Button's side, closes it at
            // the cylinder's end -(rt + cyl)
            const double rt = s.tipRadiusMm * mm, cyl = s.tipCylMm * mm;
            const ButtonCap::RoundedOutline o({ { 0.0, 0.0 }, { rt, 0.0 }, { rt, -(2.0 * rt + cyl) } });
            prof.v = o.Sweep({ 0.0, rt, 0.0 });
         }
         rp = RubberContact::RubberParams();
         rp.eMpa = s.tipEMpa;
         rp.alphaSPerM = s.tipAlpha;
         rp.gammaKK = 0.0;
         rp.betaLayer = s.tipBetaLayer;
         rp.mu = s.tipMu;
         rp.shearScale = s.tipShearScale;
         rp.shearTan = s.tipTanDelta;
         rp.tanDelta = s.tipTanDelta;
         rp.twist = true;
         rp.normalLoss = false;
         rp.thicknessMm = s.tipRubberMm;
         layerM = s.tipRubberMm * mm;
         reach = 1.0;   // a slide: the rate is the surface's speed
         maxMove = travel - outerSolidS;
         stops.clear();
         const auto stop = [&](const double at, const int side, const double reff) {
            GroupStep::BodyStopPlace st;
            st.at = at;
            st.side = side;
            st.lever = 1.0;
            st.Setup(mass, reff, s.stopE, s.stopVRef);
            return st;
         };
         stopPull = PlaceStop(stop(travel, 1, 0.5 * d));
         stopIn = PlaceStop(stop(outerSolidS, -1, 0.5 * s.outerWireMm * mm));
         GroupStep::BodyStopPlace hand = stop(0.0, -1, rRod);   // the hand's stops, off until the hand takes hold (KeyDown, Tie)
         hand.law = GroupStep::BodyStopLaw{};   // restitution 0 at every speed: the hand catches the knob, no rebound (Impact rests it on the stop)
         hand.enabled = false;
         stopHold = PlaceStop(hand);
         hand.side = 1;
         stopPush = PlaceStop(hand);
         pullRate = s.pullMmS * mm;
         Release(0.0);
      }
      // let go at pull s0 (m), at rest: no unknown held, the stops armed from their free side
      void Release(const double s0)
      {
         q = s0;
         rate = 0.0;
         pinStop = -1;
         for (GroupStep::BodyStop &st : stops)
         {
            st.touching = false;
            st.armed = st.side * (q - st.at) <= 0.0;
         }
         CycleBegin();
      }
      // ---- the hand (between cycles, before CycleBegin; rates in m/s) ----
      // key down: the hand takes the knob where it is and pulls it back to full pull; the spring presses the rod onto the hold stop, which carries it
      void KeyDown(const double rate)
      {
         SetStop(stopHold, q, true);
         DriveStop(stopHold, travel, rate);
      }
      // key up: the hand lets go; the spring alone drives the rod
      void KeyUp() { SetStop(stopHold, q, false); }
      // PullBackandRetract after its pull: the hand eases the rod back to rest on the hold stop
      void Retract(const double rate) { DriveStop(stopHold, 0.0, rate); }
      // the cabinet sensor's tie (design 2c, user U1 2026-11-18): the hold and push stops together at `at` hold the knob both ways (the rod is where the sensor says);
      // DriveTie moves both toward the reading at one rate, so the hand carries the rod either way (pull and palm alike); Untie: the hand lets go, the springs drive
      void Tie(const double at)
      {
         SetStop(stopHold, at, true);
         SetStop(stopPush, at, true);
      }
      void DriveTie(const double target, const double rate)
      {
         DriveStop(stopHold, target, rate);
         DriveStop(stopPush, target, rate);
      }
      void Untie()
      {
         SetStop(stopHold, q, false);
         SetStop(stopPush, q, false);
      }
      double GravityForce() const { return mass * g.Dot(axis); }   // along +s, N

      // ---- the body's terms ----
      double OuterEnergy(const double at) const
      {
         const double u = std::max(c0 - at, 0.0);
         return 0.5 * ko * u * u;
      }
      double StoredEnergy(const double at) const override { return km * at * (x0 + 0.5 * at) + OuterEnergy(at) - GravityForce() * at; }
      double StoredChord(const double at, const double dq, double, double &slope, double &lossF, double &lossSlope) const override
      {
         lossF = lossSlope = 0.0;
         // main spring and gravity: (V(at + dq) - V(at)) / dq exactly
         double ch = km * (x0 + at + 0.5 * dq) - GravityForce();
         slope = 0.5 * km;
         // outer spring: u = c0 - s, V = k_o u^2 / 2 where u > 0
         const double ua = c0 - at, ub = c0 - at - dq;
         if (ua >= 0.0 && ub >= 0.0)
         {
            ch += -0.5 * ko * (ua + ub);
            slope += 0.5 * ko;
         }
         else if (ua > 0.0 || ub > 0.0)   // across the kink (dq != 0): each side's energy, exactly
         {
            const double pa = std::max(ua, 0.0), pb = std::max(ub, 0.0), chO = 0.5 * ko * (pb * pb - pa * pa) / dq;
            ch += chO;
            slope += (-ko * pb - chO) / dq;   // d(chord)/d(dq) = (V'(at + dq) - chord) / dq
         }
         return ch;
      }
      double StoredForce(const double at) const override { return -km * (x0 + at) + ko * std::max(c0 - at, 0.0) + GravityForce(); }
      double Viscous() const override { return c; }
      // the bore's force on the rod across its axis, N: its weight and the balls' push, less their parts along the axis
      double BoreLoad(const V3 &ballForce) const
      {
         const V3 f = g * mass + ballForce;
         return (f - axis * f.Dot(axis)).Length();
      }
      double FrictionLimit(double, double, const V3 &ballForce) const override { return rub + boreMu * BoreLoad(ballForce); }   // the sleeve rub plus the bore on the side load
      // the rod on its bore: Mindlin's sideways spring 8 G* a of the Hertz patch (the hinge's law, lever 1)
      double PinStiffness(double, double, const V3 &ballForce) const override
      {
         if (boreMu <= 0.0 && rub <= 0.0)
            return 0.0;
         const double a = std::cbrt(3.0 * BoreLoad(ballForce) * rRod / (4.0 * kSteelEStar));
         return 8.0 * kSteelGStar * a;
      }
      // The tip's surface for a ball of radius Rb centred at c, the rod moved dq more (the trial): the ball's centre in the rod frame (z along the ball-ward axis from the apex,
      // rho across), the Revolve's closest point, its normal back to the world; the column: the rod's +s direction (the surface moves along +axis for +q)
      void ProbeSurface(const V3 &c, const double Rb, const double dq, GroupStep::Hit &h) const
      {
         const V3 u = axis * -1.0, rel = c - (tipApex + axis * (q + dq));
         const double z = rel.Dot(u);
         V3 side = rel - u * z;
         const double rho = side.Length();
         if (rho > 1e-12)
            side = side * (1.0 / rho);
         else   // on the axis: any direction across it (the normal has no part across the axis there)
         {
            side = std::fabs(u.z) < 0.9 ? u.Cross(V3(0.0, 0.0, 1.0)) : u.Cross(V3(1.0, 0.0, 0.0));
            side = side * (1.0 / side.Length());
         }
         const JetBumper::Revolve::Near nr = prof.Closest(rho, z, kEdgeMinM);
         h.delta = Rb - nr.g;
         h.n = side * nr.nr + u * nr.nz;
         h.reff = BumperGroup::KappaEff(Rb, nr);
         h.curvDist = nr.edge ? nr.g + nr.rr : 0.0;
         h.grad = 1.0;
         h.nc = 1;
         h.c[0].Set(0, h.n.Dot(axis), axis, V3());
      }
      // the part finder: |c - the apex now| less the profile's farthest point plus the travel both ways
      double RegionDist(const V3 &c) const { return (c - (tipApex + axis * q)).Length() - (BumperGroup::FarR(prof) + travel - outerSolidS); }
      // the springs' stiffness in the motion, and the outer spring's onset as an approach (accuracy only; the group holds the factors)
      void StepFacts(const GroupStep::Ball *b, const int nb, const GroupStep::Settings &gs, GroupStep::StepFactList &f) override
      {
         BodyPart::StepFacts(b, nb, gs, f);
         if (NumDof() == 0)
            return;
         const double kNow = km + (q < c0 ? ko : 0.0);
         f.Stiffness(std::sqrt(kNow / mass), std::fabs(rate) <= gs.vSteady, 2);
         if (q > c0 && rate < 0.0)
            f.Approach(q - c0, -rate, std::sqrt((km + ko) / mass), 3, [&] {   // (touch_predict: closing acceleration -a_q, PartAccel; sign as its rate)
               const double *const acc = f.Accel();
               return acc != nullptr ? -acc[0] : 0.0;
            }, GroupStep::SpringCols(0, -1.0));
      }
   };
}
