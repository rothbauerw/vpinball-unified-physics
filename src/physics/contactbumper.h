// license:GPLv3+
//
// Jet bumper (docs/JetBumperDesign.txt): a skirt that rocks on its base rim, a stem that pushes a spoon switch down, and a steel ring driven by the AE-23-800 coil.
//
// The settings, the shapes, the masses and the materials, pure SI with no VPX types. The moving is the group step's: BumperGroup::Part (bumpergroup.h), JetBumper::EnginePart (contactbumperpart.h).
//
// Frame: the bumper axis is vertical through the origin; z is up from the playfield top. Lengths are metres inside; Params are in mm where the name says Mm.

#pragma once

#include <cmath>
#include <vector>
#include <algorithm>

#include "physics/rubbercontact.h"
#include "physics/coilmodel.h"

namespace JetBumper
{
   using RubberContact::V3;
   using RubberContact::RubberParams;

   constexpr double MM = 1e-3;
   constexpr double DEG = M_PI / 180.0;

   enum class DriveMode { Held, Pulse, External };

   // ---- settings. [user] stated by the user, [derived] from user numbers, [guess] mine, unverified (design section 5) ----
   struct Params
   {
      // geometry, mm
      double bodyDiaMm = 38.0;             // [user: typically]
      double bodyTopMm = 40.0;             // [guess]
      double skirtOdMm = 63.5;             // [user]
      double skirtInnerDiaMm = 44.9;       // [skirt.obj at OD 63.5]
      double skirtTInnerMm = 2.27;         // [skirt.obj]
      double skirtTEdgeMm = 0.5;           // [user: finite rim]
      double skirtClearanceMm = 3.18;      // [user: 1/8 in "for now"]
      double lipDiaMm = 38.0;              // the skirt's inner lip meets the base rim here
      double stemLenMm = 45.0;             // [user: below the skirt]
      double stemDiaMm = 4.76;             // [user: 3/16 in]
      double tipRadiusMm = 2.38;           // [assumed: full round of the stem]
      double coneRadiusMm = 6.35;          // [user: 1/2 in across]
      double coneDepthMm = 6.35;           // [user: 1/4 in deep]
      double ringOdMm = 63.5, ringIdMm = 45.0;   // [user]
      double ringSheetMm = 1.9;            // vertical thickness [stock ring mesh]
      double ringRiseMm = 7.7;             // [stock ring mesh]
      double ringMidRestMm = 26.9875;      // [user: the mid-face of the sloped top is at the top of the ball at rest]
      double ringStrokeMm = 10.0;          // [user]
      // masses, kg [guess, accepted by the user]
      double skirtMass = 7e-3, ringMass = 105e-3, spoonMass = 1e-3;
      double plasticDensity = 1100.0;      // kg/m^3 [guess], splits the skirt mass between wafer and stem
      // springs [guess, tuned by feel]
      double skirtKNPerM = 200.0, skirtPreloadN = 0.4, springRadiusMm = 8.0, skirtZeta = 0.3;   // preloads 0.8 / 0.3 N made a ball rolling in from an inch up the 6 deg slope miss; 0.4 / 0.15 N fire it (T11)
      double spoonKNPerM = 100.0, spoonPreloadN = 0.15, spoonZeta = 0.3;
      double returnKNPerM = 100.0, returnPreloadN = 1.6, ringViscousNsPerM = 0.5;
      // stops [guess]
      double ringStopReffMm = 5.0, ringRestE = 0.4, ringCoilE = 0.3, ringRestEMpa = 30.0;   // ring stops: Hertz contact at this effective radius; rest stop rubber (E), coil stop steel on steel
      double skirtStopKNPerM = 5e5, skirtStopZeta = 0.3;
      // switch
      double switchCloseMm = 0.79375;      // [user 2026-11-02: 1/32 in]
      double switchHystMm = 0.1;           // [guess]
      // coil: AE-23-800 as the sling code has it
      DriveMode mode = DriveMode::Pulse;
      double pulseMs = 40.0;               // [user 2026-11-02: 20 was too short in play; the ring doesn't reach its stop on slow hits] (was 20 ms [user 2026-10-26: 15 ms first, 20 ms after the first play test])
      CoilSettings coil;                   // AE-23-800 defaults; vSupply 28 V [user], x0 0.3 [user accepted]
      double plungerClosedFrac = 0.0;      // plunger position fraction when the ring is fully down (0 = plunger fully in) [guess]
      // contact materials [guess]: Young's modulus of the soft partner, restitution at 0.5 m/s, friction
      double plasticEMpa = 2500.0, plasticNu = 0.35, plasticE = 0.5, plasticMu = 0.3;
      double steelRingE = 0.6, steelMu = 0.35;
      double ballRadiusM = 0.5 * 1.0625 * 0.0254, ballMassKg = 0.080;
      double edgeMinMm = 0.2;              // smallest rounding radius of a profile vertex
   };

   // ---- a surface of revolution: a polyline in (rho, z) swept by a radius rr that varies along it (a thin sheet with rounded edges) ----
   struct Revolve
   {
      struct Vtx { double r, z, rr; };
      std::vector<Vtx> v;
      struct Near
      {
         double g = 1e9;          // distance from the query point to the surface (negative inside)
         double nr = 0, nz = 1;   // unit normal in the meridian plane, surface -> point
         double rs = 0;           // radius of the surface point
         bool edge = false;       // the nearest point is a vertex (a rounded edge), not the flat of a segment
         double rr = 0;
         double gr = 0, gz = 1;   // the distance's exact gradient in the meridian plane: the normal, less the taper of the sweep's radius along the segment
      };
      Near Closest(double rho, double z, double edgeMin) const
      {
         Near best;
         for (size_t i = 0; i + 1 < v.size(); i++)
         {
            const Vtx &a = v[i], &b = v[i + 1];
            const double dr = b.r - a.r, dz = b.z - a.z, len2 = dr * dr + dz * dz;
            double t = len2 > 0.0 ? ((rho - a.r) * dr + (z - a.z) * dz) / len2 : 0.0;
            t = std::min(1.0, std::max(0.0, t));
            const double cr = a.r + t * dr, cz = a.z + t * dz, rr = a.rr + t * (b.rr - a.rr);
            const double qr = rho - cr, qz = z - cz, d = std::sqrt(qr * qr + qz * qz);
            const double g = d - rr;
            if (g < best.g)
            {
               best.g = g;
               const double len = std::sqrt(len2);
               if (d > 1e-12) { best.nr = qr / d; best.nz = qz / d; }
               else if (len > 0.0) { best.nr = -dz / len; best.nz = dr / len; }   // on the line: the normal to the left of the direction of travel
               best.rs = cr + rr * best.nr;
               best.edge = (t <= 0.0 || t >= 1.0);
               best.rr = std::max(rr, edgeMin);
               const double taper = len > 0.0 && t > 0.0 && t < 1.0 ? (b.rr - a.rr) / len : 0.0;
               best.gr = best.nr - taper * (len > 0.0 ? dr / len : 0.0);
               best.gz = best.nz - taper * (len > 0.0 ? dz / len : 0.0);
            }
         }
         return best;
      }
   };

   // What the settings give: the skirt's mass properties, the three surfaces, the materials, the coil's constants, the ring stops' stiffness. No state.
   class Bumper
   {
   public:
      Params p;
      double zRest = 0.0, zBase = 0.0;     // skirt centre height at rest; height of the base rim the lip rests against
      double cmZ = 0.0, mThetaEff = 0.0;   // skirt centre of mass (body frame z), inertia about x/y of the skirt about the rim-plane centre
      Revolve skirt0, ring0, body0;        // skirt in its own frame, ring at rest, body fixed
      RubberParams skirtMat, ringMat, bodyMat;
      CoilDerived coilD;
      double kRest = 0.0, kCoil = 0.0;     // ring stops: k (N/m^1.5)

      void Configure(const Params &prm)
      {
         p = prm;
         const double R = p.ballRadiusM;
         // --- skirt mass properties from the profile: flat bottom, top from tInner (to the inner dia) tapering to tEdge at the rim ---
         const double rOut = 0.5 * p.skirtOdMm * MM, rIn = 0.5 * p.skirtInnerDiaMm * MM, tIn = p.skirtTInnerMm * MM, tEd = p.skirtTEdgeMm * MM;
         const double rStem = 0.5 * p.stemDiaMm * MM, lStem = p.stemLenMm * MM;
         double mWafer = 0.0, iWaferDiam = 0.0, zWaferMoment = 0.0;
         const int n = 400;
         for (int i = 0; i < n; i++)
         {
            const double r = (i + 0.5) * rOut / n;
            const double th = r <= rIn ? tIn : tIn + (tEd - tIn) * (r - rIn) / (rOut - rIn);
            const double dm = p.plasticDensity * 2.0 * M_PI * r * th * (rOut / n);
            mWafer += dm;
            iWaferDiam += dm * 0.5 * r * r + dm * th * th / 12.0;
            zWaferMoment += dm * 0.5 * th;
         }
         // scale the densities so the total is the stated skirt mass (wafer + stem)
         const double stemFull = p.plasticDensity * M_PI * rStem * rStem * lStem;
         const double sc = p.skirtMass / (mWafer + stemFull);
         mWafer *= sc; iWaferDiam *= sc; zWaferMoment *= sc;
         const double mStem = stemFull * sc;
         const double zWaferCm = zWaferMoment / mWafer, zStemCm = -0.5 * lStem;
         cmZ = (mWafer * zWaferCm + mStem * zStemCm) / p.skirtMass;
         const double iStemAboutPlane = mStem * (lStem * lStem / 12.0 + zStemCm * zStemCm) + mStem * rStem * rStem / 4.0;
         const double iWaferAboutPlane = iWaferDiam + mWafer * zWaferCm * zWaferCm;
         mThetaEff = iWaferAboutPlane + iStemAboutPlane;                        // about the diameter through the centre of the bottom plane: the tilt coordinate turns about the plane centre
         zRest = p.skirtClearanceMm * MM;
         zBase = zRest + tIn;

         // --- profiles ---
         // skirt top surface (body frame, origin at the centre of the flat bottom): inner lip (r = lip radius) -> inner dia -> rim, as a variable-radius sweep ending in a round rim
         const double rLip = 0.5 * p.lipDiaMm * MM, rRim = rOut - 0.5 * tEd;
         skirt0.v = { { rLip, tIn, 0.0 }, { rIn, tIn, 0.0 }, { rRim, 0.5 * tEd, 0.5 * tEd } };
         // ring: the top face is the line from (ID/2, zTopIn) to (OD/2, zTopIn + rise); mid-face at ringMidRest at rest; the sheet is swept 0.5 * sheet (normal) below it
         {
            const double rI = 0.5 * p.ringIdMm * MM, rO = 0.5 * p.ringOdMm * MM, rise = p.ringRiseMm * MM;
            const double slope = std::atan2(rise, rO - rI);
            const double rrSheet = 0.5 * p.ringSheetMm * MM * std::cos(slope);   // half the sheet measured along the normal
            const double zTopIn = p.ringMidRestMm * MM - 0.5 * rise;
            const double down = 0.5 * p.ringSheetMm * MM;                         // the centreline is the top face moved down by half the vertical thickness (same ID and OD)
            ring0.v = { { rI, zTopIn - down, rrSheet }, { rO, zTopIn + rise - down, rrSheet } };
         }
         body0.v = { { 0.5 * p.bodyDiaMm * MM, zBase, 0.5 * MM }, { 0.5 * p.bodyDiaMm * MM, p.bodyTopMm * MM, 0.5 * MM }, { 0.0, p.bodyTopMm * MM, 0.5 * MM } };

         // --- materials ---
         const double steelEstarMpa = 200e9 / (2.0 * (1.0 - 0.3 * 0.3)) * 1e-6;   // steel on steel: E* = E / (2 (1 - nu^2))
         {
            RubberParams rp;
            rp.eMpa = p.plasticEMpa; rp.nu = p.plasticNu; rp.thicknessMm = 100.0; rp.betaLayer = 0.0; rp.mu = p.plasticMu; rp.twist = true; rp.normalLoss = false;
            rp.gammaKK = RubberContact::CalibrateKK(rp, R, p.ballMassKg, 0.5, p.plasticE);
            skirtMat = rp;
            bodyMat = rp;
         }
         {
            RubberParams rp;
            rp.eMpa = steelEstarMpa; rp.nu = 0.0;   // E* itself (as HingeBody::SetupStop)
            rp.thicknessMm = 100.0; rp.betaLayer = 0.0; rp.mu = p.steelMu; rp.twist = true; rp.normalLoss = false;
            rp.gammaKK = RubberContact::CalibrateKK(rp, R, p.ballMassKg, 0.5, p.steelRingE);
            ringMat = rp;
         }
         // ring stops: the Hertz stiffness of each, for the stop's contact time
         {
            const double sq = std::sqrt(p.ringStopReffMm * MM);
            kCoil = (4.0 / 3.0) * steelEstarMpa * 1e6 * sq;
            kRest = (4.0 / 3.0) * (p.ringRestEMpa * 1e6 / (1.0 - 0.25)) * sq;
         }
         // coil constants (the sling's law)
         CoilModel::Derive(p.coil, coilD);
      }
   };
}
