// license:GPLv3+
//
// The auto plunger (docs/unified-physics/PlungerDesign.txt section 4; piece 2, offline). One body on one coordinate (bodypart.h), a straight slide along the lane:
// s = the carriage's travel, 0 at rest, + toward the ball (up the lane). Carriage, head, the coil's plunger and linkage move as one (4a, 4f).
//   drive         on Fire a constant force F = autoEnergyJ / stroke along +s (BodyPart's DriveForce, booked as a drive), off once the carriage reaches the end stop:
//                 its work is F x stroke = autoEnergyJ (4d; a stand-in for the coil)
//   return spring a light linear spring toward the rest stop: F = -(P + k s) (4e)
//   gravity       the part's gravity along the lane: m g.axis toward +s (the table's slope is in g)
//   stops         rest (s = 0, side -1) and end of stroke (s = stroke, side +1, the coil's plunger bottoming): steel on steel BodyStops, the impact law (4c)
//   head          fork (default, 4b): two steel tines, each a rounded corner (a sphere of radius tineRadiusMm) placed so that a ball centred at the cradle point
//                 meets it with the derived normal (4b: sideways, down and along the lane); or tip (4h): the mechanical rod's revolved ring-rubber tip (plungerrod.h)
//   losses        none in this piece (no drag, no friction): the stops' impact law and the head's contact damping only
// Stored energy V(s) = P s + k s^2 / 2 - m (g.axis) s: its chord is exact, so the book closes at any step. Pure SI inside; settings in mm, g, N, N/mm, J.

#pragma once

#include <algorithm>
#include <cmath>

#include "physics/plungerrod.h"

namespace AutoPlunger
{
   using GroupStep::V3;

   enum class Head { Fork, Tip };

   // Every value: MEASURED (user's hardware), DERIVED (lead, from measured), ESTIMATED (no data); source PartData.txt section 11 / PlungerDesign.txt section 4
   struct Settings
   {
      double massG = 80.0;              // carriage and forks about 50 g, the coil's plunger about 20 g, linkage: ESTIMATED (4f)
      double strokeMm = 27.0;           // about one ball diameter: ESTIMATED (user, 4c)
      double autoEnergyJ = 1.0;         // the drive's work over the stroke: ESTIMATED (4d; per table, [Plunger:<name>] auto_energy_j)
      Head head = Head::Fork;           // fork (default) or tip (4h)
      // return spring (4e): ESTIMATED (the photo's spring on the coil plunger): holds the carriage on its rest stop, returns it from the end in well under 0.2 s
      double returnPreloadN = 1.0;
      double returnRateNmm = 0.1;
      // end and rest stops: steel on steel, restitution stopE at stopVRef m/s, Hertz radius stopRadiusMm: ESTIMATED as plungerrod.h's steel stops (coil plunger bottoming not measured)
      double stopE = 0.5;
      double stopVRef = 0.1;
      double stopRadiusMm = 5.0;        // ESTIMATED (the coil plunger's end; only sets the stop's contact stiffness)
      // fork (4b): user's estimates (PartData 11)
      double forkGapMm = 20.0;          // inside, plate to plate: ESTIMATED (user, 3/4 of the ball)
      double forkHeightMm = 18.0;       // contact height above the playfield: ESTIMATED (user, 2/3 up the ball)
      double tineRadiusMm = 1.0;        // the tine's rounded corner: ESTIMATED (a strip about 2 mm thick, its corner rounded to half that)
      // the tines' steel: as the jet bumper's steel ring (contactbumper.h: E* of steel on steel, mu 0.35, restitution 0.6 on a ball at 0.5 m/s): ESTIMATED
      double steelMu = 0.35;
      double steelE = 0.6;
      double steelVRef = 0.5;
      // tip head (4h): the mechanical rod's tip, its own values (plungerrod.h Settings: radius, cylinder, rubber)
      PlungerRod::Settings tip;
   };

   class Part : public GroupStep::BodyPart
   {
   public:
      V3 axis = V3(0.0, -1.0, 0.0);   // the unit direction of +s (toward the ball) in the world
      V3 g = V3(0.0, 0.0, -9.81);     // the gravity the part feels
      double kRet = 0.0, pRet = 0.0;  // return spring, N/m and N
      double stroke = 0.0, force = 0.0;
      bool firing = false;            // the drive is on (Fire until the end stop)
      int stopRest = -1, stopEnd = -1;
      Head head = Head::Fork;
      V3 cradle;                      // a ball centred here touches the head at s = 0
      V3 tine[2];                     // fork: the tines' corner centres at s = 0
      double rTine = 0.0;
      RubberContact::RubberParams forkRp;   // the tines' ball law
      // tip head
      V3 tipApex;
      JetBumper::Revolve prof;
      RubberContact::RubberParams tipRp;
      double tipLayerM = 0.0;
      static constexpr double kEdgeMinM = 0.2e-3;

      // axisDir: +s in the world (unit), toward the ball; up: the playfield's normal; gravity: the part's gravity now; cradleAt: the centre of a ball (radius Rb, mass mb)
      // that touches the head at s = 0
      void Build(const Settings &s, const V3 &axisDir, const V3 &up, const V3 &gravity, const V3 &cradleAt, const double Rb, const double mb)
      {
         const double mm = 1e-3;
         axis = axisDir;
         g = gravity;
         head = s.head;
         mass = s.massG * 1e-3;
         stroke = s.strokeMm * mm;
         force = s.autoEnergyJ / stroke;
         kRet = s.returnRateNmm * 1e3;
         pRet = s.returnPreloadN;
         cradle = cradleAt;
         // fork: the contact point on the ball, from its centre: half the gap sideways, the height above the centre up, the rest of the radius back along -axis (4b)
         {
            V3 side = axis.Cross(up);
            side = side * (1.0 / side.Length());
            const double a = 0.5 * s.forkGapMm * mm, b = s.forkHeightMm * mm - Rb, c = std::sqrt(std::max(Rb * Rb - a * a - b * b, 0.0));
            rTine = s.tineRadiusMm * mm;
            for (int i = 0; i < 2; i++)
            {
               const V3 u = (side * (i == 0 ? a : -a) + up * b - axis * c) * (1.0 / Rb);   // unit, ball centre to contact point
               tine[i] = cradle + u * (Rb + rTine);
            }
            const double steelEstarMpa = 200e9 / (2.0 * (1.0 - 0.3 * 0.3)) * 1e-6, reff = Rb * rTine / (Rb + rTine);
            forkRp = RubberContact::RubberParams();
            forkRp.eMpa = steelEstarMpa; forkRp.nu = 0.0;   // E* itself (as the bumper's ring)
            forkRp.thicknessMm = 100.0; forkRp.betaLayer = 0.0; forkRp.mu = s.steelMu; forkRp.twist = true; forkRp.normalLoss = false;
            forkRp.gammaKK = RubberContact::CalibrateKK(forkRp, reff, mb, s.steelVRef, s.steelE);
         }
         // tip: the rod's profile and rubber (plungerrod.h Build), apex on the axis Rb behind the cradle
         {
            const PlungerRod::Settings &t = s.tip;
            const double rt = t.tipRadiusMm * mm, cyl = t.tipCylMm * mm;
            const ButtonCap::RoundedOutline o({ { 0.0, 0.0 }, { rt, 0.0 }, { rt, -(2.0 * rt + cyl) } });
            prof.v = o.Sweep({ 0.0, rt, 0.0 });
            tipApex = cradle - axis * Rb;
            tipRp = RubberContact::RubberParams();
            tipRp.eMpa = t.tipEMpa; tipRp.alphaSPerM = t.tipAlpha; tipRp.gammaKK = 0.0; tipRp.betaLayer = t.tipBetaLayer; tipRp.mu = t.tipMu;
            tipRp.shearScale = t.tipShearScale; tipRp.shearTan = t.tipTanDelta; tipRp.tanDelta = t.tipTanDelta; tipRp.twist = true; tipRp.normalLoss = false;
            tipRp.thicknessMm = t.tipRubberMm;
            tipLayerM = t.tipRubberMm * mm;
         }
         reach = 1.0;   // a slide: the rate is the surface's speed
         maxMove = stroke;
         stops.clear();
         const auto stop = [&](const double at, const int side) {
            GroupStep::BodyStopPlace st;
            st.at = at;
            st.side = side;
            st.lever = 1.0;
            st.Setup(mass, s.stopRadiusMm * mm, s.stopE, s.stopVRef);
            return st;
         };
         stopRest = PlaceStop(stop(0.0, -1));
         stopEnd = PlaceStop(stop(stroke, 1));
         q = rate = 0.0;
         firing = false;
         CycleBegin();
      }
      // between cycles, before CycleBegin: the coil fires; the drive stays on until the carriage reaches the end stop
      void Fire() { firing = true; }
      double GravityForce() const { return mass * g.Dot(axis); }   // along +s, N

      // ---- the body's terms ----
      double StoredEnergy(const double at) const override { return pRet * at + 0.5 * kRet * at * at - GravityForce() * at; }
      double StoredChord(const double at, const double dq, double, double &slope, double &lossF, double &lossSlope) const override
      {
         lossF = lossSlope = 0.0;
         slope = 0.5 * kRet;
         return pRet + kRet * (at + 0.5 * dq) - GravityForce();   // (V(at + dq) - V(at)) / dq exactly
      }
      double StoredForce(const double at) const override { return -(pRet + kRet * at) + GravityForce(); }
      double DriveForce() const override { return firing ? force : 0.0; }
      double FrictionLimit(double, double, const V3 &) const override { return 0.0; }
      // the coil's plunger has bottomed: the drive is off (its work F x stroke is done)
      void AfterMove() override
      {
         if (firing && stops[stopEnd].touching)
            firing = false;
      }

      // ---- the head (a carried surface; the column: +s moves it along +axis) ----
      // fork tine i (0, 1): a steel sphere of radius rTine at tine[i] + axis (q + dq)
      void ProbeTine(const int i, const V3 &c, const double Rb, const double dq, GroupStep::Hit &h) const
      {
         const V3 rel = c - (tine[i] + axis * (q + dq));
         const double d = rel.Length();
         h.n = d > 1e-12 ? rel * (1.0 / d) : axis;
         h.delta = Rb + rTine - d;
         h.reff = Rb * rTine / (Rb + rTine);
         h.curvDist = d;
         h.grad = 1.0;
         h.nc = 1;
         h.c[0].Set(0, h.n.Dot(axis), axis, V3());
      }
      // tip: the revolved profile about the axis, apex at tipApex + axis (q + dq), z along +axis (toward the ball) from the apex (plungerrod.h ProbeSurface)
      void ProbeTip(const V3 &c, const double Rb, const double dq, GroupStep::Hit &h) const
      {
         const V3 u = axis, rel = c - (tipApex + axis * (q + dq));
         const double z = rel.Dot(u);
         V3 side = rel - u * z;
         const double rho = side.Length();
         if (rho > 1e-12)
            side = side * (1.0 / rho);
         else
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
      void StepFacts(const GroupStep::Ball *b, const int nb, const GroupStep::Settings &gs, GroupStep::StepFactList &f) override
      {
         BodyPart::StepFacts(b, nb, gs, f);
         if (NumDof() == 0)
            return;
         f.Stiffness(std::sqrt(kRet / mass), std::fabs(rate) <= gs.vSteady, 2);
      }
   };
}
