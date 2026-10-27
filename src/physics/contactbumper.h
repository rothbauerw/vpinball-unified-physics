// license:GPLv3+
//
// Jet bumper (docs/JetBumperDesign.txt): a skirt that rocks on its base rim, a stem that pushes a spoon switch down, and a steel ring driven by the AE-23-800 coil.
//
// Pure SI maths with no VPX types, so the same code runs inside VPX and in the offline tests (runs/2026-10-26-bumper). Step 1 of the build order: the parts, the
// switch, the coil and the three drive modes, plus the ball's contact with the skirt, the ring and the centre body (the same Hertz + viscoelastic + Mindlin law as
// every other surface, RubberContact::ContactCore). Nothing here is wired into the engine yet.
//
// Frame: the bumper axis is vertical through the origin; z is up from the playfield top. Lengths are metres inside; Params are in mm where the name says Mm.
//
// Bodies
//   skirt  3 degrees of freedom: the height z of the centre of its flat bottom and two tilt angles (thx about x, thy about y; R = Rx(thx) Ry(thy), exact rotation).
//          The lateral position is held (the base takes any sideways force at the lip; ASSUMPTION, design 2.2a). Equations of motion on the generalised coordinates
//          with the virtual-work Jacobian of each force point: m z'' = Qz, Mth th'' = Qth (the mass matrix is taken at the rest pose, where it is diagonal; the velocity
//          terms are second order in the angles, at most about 5 degrees). Forces: weight (minus the cabinet acceleration) at the centre of mass, 8 equal springs on a
//          circle under the skirt (axial rate k and tilt rate k a^2 / 2 from one set of numbers), a one-sided stop on the lip circle against the base rim, a one-sided
//          stop of the rim circle against the playfield, the stem tip on the spoon, the balls.
//   spoon  one vertical degree of freedom (ds, down positive), spring returned, with an upper stop (rest) and a lower stop. Its upper face is a cone.
//   ring   ring + rods + yoke + plunger as one body, one vertical degree of freedom (dr, down positive): coil force, weight, return spring (preload + rate),
//          viscous friction, rubber rest stop (up), steel coil stop (down), the balls.
// Switch   closed when ds >= close, open again when ds <= close - hysteresis. Hit is raised at the instant it closes, before the coil is commanded.
// Coil     di/dt = (V - R i - i dL/dx v) / L with the sling's single-winding law (L = a + b / (frac + x0)); force on the plunger = strength * g * i^2.
//          Drive: held (on while the switch is closed), pulse (one timed pulse per closure, re-armed when the switch has reopened), external (Solenoid()).
//          Off = decay through R (freewheel diode across the coil ASSUMED, as the sling).
// Integration: semi-implicit Euler with a caller-chosen step (the offline tests use 1 us; the engine will take the step from the contact stiffness, rule 5).

#pragma once

#include <cmath>
#include <vector>
#include <algorithm>

#include "physics/rubbercontact.h"
#include "physics/coilmodel.h"

namespace JetBumper
{
   using RubberContact::V3;
   using RubberContact::BallSI;
   using RubberContact::RubberParams;
   using RubberContact::RubberState;

   constexpr double MM = 1e-3;
   constexpr double DEG = M_PI / 180.0;
   constexpr double kSteelEstar = 200e9 / (2.0 * (1.0 - 0.3 * 0.3));   // steel on steel: E* = E / (2 (1 - nu^2))

   enum class DriveMode { Held, Pulse, External };

   enum EventKind { EvHit, EvSwitchClosed, EvSwitchOpened, EvPulseStart, EvPulseEnd, EvCoilStop, EvRestStop };
   struct Event { double t; int kind; double value; };

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
      double spoonKNPerM = 100.0, spoonPreloadN = 0.15, spoonZeta = 0.3, spoonTravelMm = 5.0;
      double returnKNPerM = 100.0, returnPreloadN = 1.6, ringViscousNsPerM = 0.5;
      // stops [guess]
      double ringStopReffMm = 5.0, ringRestE = 0.4, ringCoilE = 0.3, ringRestEMpa = 30.0;   // ring stops: Hertz contact at this effective radius; rest stop rubber (E), coil stop steel on steel
      double skirtStopKNPerM = 5e5, skirtStopZeta = 0.3, spoonStopKNPerM = 2e5;
      // switch
      double switchCloseMm = 1.59;         // [user: 1/16 in typical, up to 1/8 in]
      double switchHystMm = 0.1;           // [guess]
      // coil: AE-23-800 as the sling code has it
      DriveMode mode = DriveMode::Pulse;
      double pulseMs = 20.0;               // [user 2026-10-26: 15 ms first, 20 ms after the first play test]
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
            }
         }
         return best;
      }
   };

   struct BallContactState { RubberState skirt, ring, body; };
   struct BallLoad
   {
      V3 force, torque;            // on the ball, N and N m (the playfield is the caller's)
      double fnSkirt = 0, fnRing = 0, fnBody = 0;
      double gapSkirt = 1e9, gapRing = 1e9, gapBody = 1e9;   // distance of the ball surface from each surface (negative = touching)
      V3 normal;                   // the contact normal (surface -> ball) of the deepest contact
      double delta = 0.0;          // its compression, m
   };

   // the nearest of the three surfaces to a ball (centre relative to the axis foot)
   struct Nearest
   {
      double g = 1e9;              // gap between the ball's surface and the surface, m (negative = overlap)
      V3 n = V3(0.0, 0.0, 1.0);    // normal, surface -> ball
      V3 point;                    // the surface point (relative to the axis foot)
      double reff = 0.0;           // Hertz effective radius
      int surf = -1;               // 0 skirt, 1 ring, 2 body
   };

   class Bumper
   {
   public:
      Params p;

      // ---- state ----
      double zc = 0.0, thx = 0.0, thy = 0.0, vz = 0.0, wx = 0.0, wy = 0.0;   // skirt
      double ds = 0.0, vds = 0.0;                                            // spoon (m, m/s, down positive)
      double dr = 0.0, vr = 0.0;                                             // ring (m, m/s, down positive)
      double iCoil = 0.0;                                                    // A
      bool switchClosed = false;
      double t = 0.0;
      std::vector<Event> events;
      V3 gEff = V3(0.0, 0.0, -9.81);                                         // gravity minus the cabinet acceleration

      // detectors / audit (rule 8, design 9)
      double coilOnS = 0.0;                // time the coil has been commanded on, this episode
      double coilOnMaxS = 0.0;
      double wCoil = 0.0;                  // work of the coil force on the ring, J
      double eDissRing = 0.0;              // dissipated in the ring's dampers, J
      bool tipOutsideCone = false;

      // ---- derived ----
      double zRest = 0.0, zBase = 0.0;     // skirt centre height at rest; height of the base rim the lip rests against
      double cmZ = 0.0, mThetaEff = 0.0, iCmT = 0.0;   // skirt centre of mass (body frame z), inertia about x/y of the skirt about the rim-plane centre
      double stemMass = 0.0;
      Revolve skirt0, ring0, body0;        // skirt in its own frame, ring at rest, body fixed
      RubberParams skirtMat, ringMat, bodyMat, tipMat;
      double coneApexRest = 0.0;           // z of the cone apex at rest
      CoilDerived coilD;
      double kRest = 0.0, gRest = 0.0, kCoil = 0.0, gCoil = 0.0;   // ring stops: k (N/m^1.5), gamma (N s/m^2)

      // drive state
      bool pulseActive = false;
      double pulseLeftS = 0.0;
      bool pulseArmed = true;
      bool extOn = false;
      double extVolts = -1.0, extDuty = 1.0;
      bool coilCommand = false;
      double coilVolts = 0.0;
      bool prevClosed = false;
      bool asleep = false;                 // nothing moves, nothing is driven: StepHeld does nothing until something wakes it
      bool collidable = true;              // false: the bumper is not there at all (no contact, no switch, no coil, no events, parts at rest)

      void Configure(const Params &prm)
      {
         p = prm;
         const double R = p.ballRadiusM;
         // --- skirt mass properties from the profile: flat bottom, top from tInner (to the inner dia) tapering to tEdge at the rim ---
         const double rOut = 0.5 * p.skirtOdMm * MM, rIn = 0.5 * p.skirtInnerDiaMm * MM, tIn = p.skirtTInnerMm * MM, tEd = p.skirtTEdgeMm * MM;
         const double rStem = 0.5 * p.stemDiaMm * MM, lStem = p.stemLenMm * MM;
         stemMass = p.plasticDensity * M_PI * rStem * rStem * lStem;
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
         const double sc = p.skirtMass / (mWafer + stemMass);
         mWafer *= sc; iWaferDiam *= sc; zWaferMoment *= sc;
         const double mStem = stemMass * sc;
         stemMass = mStem;
         const double zWaferCm = zWaferMoment / mWafer, zStemCm = -0.5 * lStem;
         cmZ = (mWafer * zWaferCm + mStem * zStemCm) / p.skirtMass;
         const double iStemAboutPlane = mStem * (lStem * lStem / 12.0 + zStemCm * zStemCm) + mStem * rStem * rStem / 4.0;
         const double iWaferAboutPlane = iWaferDiam + mWafer * zWaferCm * zWaferCm;
         const double iPlane = iWaferAboutPlane + iStemAboutPlane;              // about the diameter through the centre of the bottom plane
         iCmT = iPlane - p.skirtMass * cmZ * cmZ;
         mThetaEff = iPlane;                                                    // the tilt coordinate turns about the plane centre
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
         auto plastic = [&](double e) {
            RubberParams rp;
            rp.eMpa = p.plasticEMpa; rp.nu = p.plasticNu; rp.thicknessMm = 100.0; rp.betaLayer = 0.0; rp.mu = p.plasticMu; rp.twist = true; rp.normalLoss = false;
            rp.gammaKK = RubberContact::CalibrateKK(rp, R, p.ballMassKg, 0.5, e);
            return rp;
         };
         skirtMat = plastic(p.plasticE);
         bodyMat = skirtMat;
         tipMat = skirtMat;
         tipMat.mu = p.plasticMu;
         {
            RubberParams rp;
            rp.eMpa = 200e9 / (2.0 * (1.0 - 0.3 * 0.3)) * 1e-6; rp.nu = 0.0;   // E* itself (as HingeBody::SetupStop)
            rp.thicknessMm = 100.0; rp.betaLayer = 0.0; rp.mu = p.steelMu; rp.twist = true; rp.normalLoss = false;
            rp.gammaKK = RubberContact::CalibrateKK(rp, R, p.ballMassKg, 0.5, p.steelRingE);
            ringMat = rp;
         }
         // spoon cone apex at rest: the tip (a sphere) just touches the 45 degree wall with its centre on the axis
         {
            const double alpha = std::atan2(p.coneDepthMm, p.coneRadiusMm);
            const double hSphere = p.tipRadiusMm * MM / std::cos(alpha);
            const double zTipCentre = zRest - (p.stemLenMm - p.tipRadiusMm) * MM;
            coneApexRest = zTipCentre - hSphere;
         }
         // ring stops, Hertz + viscoelastic (Kuwabara-Kono) as the hinge stops: restitution e at 0.5 m/s for the ring's mass
         {
            const double reff = p.ringStopReffMm * MM, sq = std::sqrt(reff);
            RubberParams rp;
            rp.thicknessMm = 100.0; rp.betaLayer = 0.0; rp.normalLoss = false;
            rp.eMpa = 200e9 / (2.0 * (1.0 - 0.3 * 0.3)) * 1e-6; rp.nu = 0.0;
            kCoil = (4.0 / 3.0) * rp.eMpa * 1e6 * sq;
            gCoil = RubberContact::CalibrateKK(rp, reff, p.ringMass, 0.5, p.ringCoilE);
            rp.eMpa = p.ringRestEMpa; rp.nu = 0.5;
            kRest = (4.0 / 3.0) * (p.ringRestEMpa * 1e6 / (1.0 - 0.25)) * sq;
            gRest = RubberContact::CalibrateKK(rp, reff, p.ringMass, 0.5, p.ringRestE);
         }
         // coil constants (the sling's law)
         CoilModel::Derive(p.coil, coilD);
         Reset();
      }

      void Reset()
      {
         zc = zRest; thx = thy = vz = wx = wy = 0.0;
         ds = vds = 0.0; dr = vr = 0.0; iCoil = 0.0;
         switchClosed = false; prevClosed = false; t = 0.0; events.clear();
         pulseActive = false; pulseLeftS = 0.0; pulseArmed = true; extOn = false; extVolts = -1.0; extDuty = 1.0; coilCommand = false; coilVolts = 0.0;
         coilOnS = coilOnMaxS = 0.0; wCoil = eDissRing = 0.0; tipOutsideCone = false;
         m_Qext[0] = m_Qext[1] = m_Qext[2] = 0.0; m_ringExt = 0.0;
         gEff = V3(0.0, 0.0, -9.81);
      }

      // ---- script surface ----
      // Collidable = False takes the bumper out of the world: parts return to rest, coil off, balls pass through, Solenoid ignored. True puts it back at rest.
      void SetCollidable(bool c)
      {
         if (c == collidable)
            return;
         const V3 g = gEff;
         Reset();
         gEff = g;
         collidable = c;
         asleep = false;
      }
      // Bumper.Solenoid enabled[, volts[, duty]]: external mode sets the coil on/off; pulse mode fires one pulse on True; held mode ignores it (the caller warns once)
      bool Solenoid(bool on, double volts = -1.0, double duty = 1.0)
      {
         if (!collidable || p.mode == DriveMode::Held)
            return false;
         asleep = false;
         if (p.mode == DriveMode::External)
         {
            extOn = on; extVolts = on ? volts : -1.0; extDuty = duty;
         }
         else if (on && !pulseActive)
            StartPulse();
         return true;
      }
      double RotXDeg() const { return thx / DEG; }
      double RotYDeg() const { return thy / DEG; }
      double RingDropM() const { return dr; }
      double SkirtZM() const { return zc; }
      double SpoonDropM() const { return ds; }

      // ---- skirt kinematics ----
      V3 SkirtPoint(const V3 &b) const
      {
         const double cx = std::cos(thx), sx = std::sin(thx), cy = std::cos(thy), sy = std::sin(thy);
         const V3 a(cy * b.x + sy * b.z, b.y, -sy * b.x + cy * b.z);
         return V3(a.x, cx * a.y - sx * a.z, sx * a.y + cx * a.z + zc);
      }
      V3 ToBody(const V3 &w) const
      {
         const double cx = std::cos(thx), sx = std::sin(thx), cy = std::cos(thy), sy = std::sin(thy);
         const V3 r(w.x, w.y, w.z - zc);
         const V3 a(r.x, cx * r.y + sx * r.z, -sx * r.y + cx * r.z);
         return V3(cy * a.x - sy * a.z, a.y, sy * a.x + cy * a.z);
      }
      V3 SkirtRotate(const V3 &b) const   // body direction -> world direction
      {
         const double cx = std::cos(thx), sx = std::sin(thx), cy = std::cos(thy), sy = std::sin(thy);
         const V3 a(cy * b.x + sy * b.z, b.y, -sy * b.x + cy * b.z);
         return V3(a.x, cx * a.y - sx * a.z, sx * a.y + cx * a.z);
      }
      // d(point)/d(z, thx, thy) at body point b
      void Jacobian(const V3 &b, V3 col[3]) const
      {
         const double cx = std::cos(thx), sx = std::sin(thx), cy = std::cos(thy), sy = std::sin(thy);
         const V3 a(cy * b.x + sy * b.z, b.y, -sy * b.x + cy * b.z);
         col[0] = V3(0.0, 0.0, 1.0);
         col[1] = V3(0.0, -sx * a.y - cx * a.z, cx * a.y - sx * a.z);
         const V3 da(-sy * b.x + cy * b.z, 0.0, -cy * b.x - sy * b.z);
         col[2] = V3(da.x, -sx * da.z, cx * da.z);
      }
      V3 SkirtPointVelocity(const V3 &b) const
      {
         V3 col[3];
         Jacobian(b, col);
         return col[0] * vz + col[1] * wx + col[2] * wy;
      }
      V3 SkirtOmega() const { return V3(wx, std::cos(thx) * wy, std::sin(thx) * wy); }

      // ---- the ball's contact with the skirt, the ring and the centre body ----
      // Adds the reactions to the bodies (applied at the next Step) and returns the force and torque on the ball. c, v, w: ball centre, velocity, spin (world).
      void BallForces(const BallSI &ball, const V3 &c, const V3 &v, const V3 &w, double h, BallContactState &st, BallLoad &out)
      {
         if (!collidable)
         {
            st.skirt.Reset(); st.ring.Reset(); st.body.Reset();
            return;
         }
         const double R = ball.R;
         const double rho = std::sqrt(c.x * c.x + c.y * c.y);
         const double cphi = rho > 1e-9 ? c.x / rho : 1.0, sphi = rho > 1e-9 ? c.y / rho : 0.0;
         const double edgeMin = p.edgeMinMm * MM;
         auto kappaEff = [&](const Revolve::Near &nr) { return KappaEff(R, nr); };
         // skirt (body frame)
         {
            const V3 b = ToBody(c);
            const double rb = std::sqrt(b.x * b.x + b.y * b.y);
            const double cb = rb > 1e-9 ? b.x / rb : 1.0, sb = rb > 1e-9 ? b.y / rb : 0.0;
            const Revolve::Near nr = skirt0.Closest(rb, b.z, edgeMin);
            out.gapSkirt = nr.g - R;
            if (nr.g < R)
            {
               const V3 nb(nr.nr * cb, nr.nr * sb, nr.nz);
               const V3 n = SkirtRotate(nb);
               const double delta = R - nr.g;
               const V3 pc = c + n * (-(R - 0.5 * delta));
               const V3 vs = SkirtPointVelocity(ToBody(pc));
               const auto res = RubberContact::ContactCore(skirtMat, ball, n, delta, kappaEff(nr), skirtMat.thicknessMm * MM, vs, v, w, h, st.skirt, SkirtOmega());
               out.force = out.force + res.force;
               out.torque = out.torque + res.torque;
               out.fnSkirt = res.fn;
               if (delta > out.delta) { out.delta = delta; out.normal = n; }
               asleep = false;
               AddSkirtForce(ToBody(pc), res.force * -1.0, m_Qext);
            }
            else
               st.skirt.Reset();
         }
         // ring (translated down by dr: query the profile at z + dr)
         {
            const Revolve::Near nr = ring0.Closest(rho, c.z + dr, edgeMin);
            out.gapRing = nr.g - R;
            if (nr.g < R)
            {
               const V3 n(nr.nr * cphi, nr.nr * sphi, nr.nz);
               const double delta = R - nr.g;
               const auto res = RubberContact::ContactCore(ringMat, ball, n, delta, kappaEff(nr), ringMat.thicknessMm * MM, V3(0.0, 0.0, -vr), v, w, h, st.ring);
               out.force = out.force + res.force;
               out.torque = out.torque + res.torque;
               out.fnRing = res.fn;
               if (delta > out.delta) { out.delta = delta; out.normal = n; }
               asleep = false;
               m_ringExt += res.force.z;   // force on the ring is -F; down positive: +F.z
            }
            else
               st.ring.Reset();
         }
         // centre body (fixed)
         {
            const Revolve::Near nr = body0.Closest(rho, c.z, edgeMin);
            out.gapBody = nr.g - R;
            if (nr.g < R)
            {
               const V3 n(nr.nr * cphi, nr.nr * sphi, nr.nz);
               const double delta = R - nr.g;
               const auto res = RubberContact::ContactCore(bodyMat, ball, n, delta, kappaEff(nr), bodyMat.thicknessMm * MM, V3(), v, w, h, st.body);
               out.force = out.force + res.force;
               out.torque = out.torque + res.torque;
               out.fnBody = res.fn;
               if (delta > out.delta) { out.delta = delta; out.normal = n; }
               asleep = false;
            }
            else
               st.body.Reset();
         }
      }

      // The nearest of the three surfaces to a ball of radius R centred at c (relative to the axis foot): the engine's probe (distance, normal, curvature)
      Nearest NearestSurface(const V3 &c, double R) const
      {
         Nearest out;
         const double edgeMin = p.edgeMinMm * MM;
         const double rho = std::sqrt(c.x * c.x + c.y * c.y);
         const double cphi = rho > 1e-9 ? c.x / rho : 1.0, sphi = rho > 1e-9 ? c.y / rho : 0.0;
         auto consider = [&](const Revolve::Near &nr, const V3 &n, const V3 &pt, int surf) {
            if (nr.g < out.g)
            {
               out.g = nr.g - R;
               out.n = n;
               out.point = pt;
               out.reff = KappaEff(R, nr);
               out.surf = surf;
            }
         };
         double best = 1e9;
         {
            const V3 b = ToBody(c);
            const double rb = std::sqrt(b.x * b.x + b.y * b.y);
            const double cb = rb > 1e-9 ? b.x / rb : 1.0, sb = rb > 1e-9 ? b.y / rb : 0.0;
            const Revolve::Near nr = skirt0.Closest(rb, b.z, edgeMin);
            const V3 n = SkirtRotate(V3(nr.nr * cb, nr.nr * sb, nr.nz));
            if (nr.g < best) { best = nr.g; out.g = 1e9; consider(nr, n, c - n * nr.g, 0); }
         }
         {
            const Revolve::Near nr = ring0.Closest(rho, c.z + dr, edgeMin);
            const V3 n(nr.nr * cphi, nr.nr * sphi, nr.nz);
            if (nr.g < best) { best = nr.g; out.g = 1e9; consider(nr, n, c - n * nr.g, 1); }
         }
         {
            const Revolve::Near nr = body0.Closest(rho, c.z, edgeMin);
            const V3 n(nr.nr * cphi, nr.nr * sphi, nr.nz);
            if (nr.g < best) { best = nr.g; out.g = 1e9; consider(nr, n, c - n * nr.g, 2); }
         }
         return out;
      }

      // ---- stepping for the engine ----
      // dt seconds in steps no longer than hMax, with the external forces accumulated since the last call (the balls' reactions of one island substep) held over all of them,
      // as a gate's body holds its torque. Returns the number of steps (0 while asleep).
      int StepHeld(double dt, double hMax)
      {
         if (!collidable)
            return 0;
         const bool hasExt = m_Qext[0] != 0.0 || m_Qext[1] != 0.0 || m_Qext[2] != 0.0 || m_ringExt != 0.0;
         if (asleep && !hasExt)
            return 0;
         asleep = false;
         const double Q[3] = { m_Qext[0], m_Qext[1], m_Qext[2] }, Rf = m_ringExt;
         const int n = std::max(1, (int)std::ceil(dt / hMax - 1e-9));
         const double hs = dt / n;
         for (int i = 0; i < n; i++)
         {
            m_Qext[0] = Q[0]; m_Qext[1] = Q[1]; m_Qext[2] = Q[2]; m_ringExt = Rf;
            Step(hs);
         }
         m_Qext[0] = m_Qext[1] = m_Qext[2] = 0.0;
         m_ringExt = 0.0;
         if (!hasExt && AtRest())
            asleep = true;
         return n;
      }
      // the step the bumper's own bodies need (s): fine while a ring stop is in contact or about to be, otherwise the skirt / spoon stops and the tip contact set it
      double OwnStep() const
      {
         const double stroke = p.ringStrokeMm * MM;
         const bool nearStop = m_inRest || m_inCoilStop || (std::abs(vr) > 0.02 && (dr < 0.3e-3 || dr > stroke - 0.3e-3));
         return nearStop ? kStepStop : kStepFree;
      }
      static constexpr double kStepFree = 4e-6, kStepStop = 1e-6;
      // resting: nothing is moving, driven or pressed
      bool AtRest() const
      {
         if (switchClosed || coilCommand || pulseActive || iCoil > 1e-3)
            return false;
         if (p.mode == DriveMode::Pulse && !pulseArmed)
            return false;
         if (std::abs(vz) > 1e-4 || std::abs(wx) > 1e-4 || std::abs(wy) > 1e-4 || std::abs(zc - zRest) > 5e-6 || std::abs(thx) > 1e-6 || std::abs(thy) > 1e-6)
            return false;
         if (std::abs(vds) > 1e-4 || std::abs(ds) > 5e-6 || std::abs(vr) > 1e-4 || std::abs(dr) > 5e-5)
            return false;
         return true;
      }

      // ---- the dynamic state, for rewinding an island (ContactMovingPart::State) ----
      struct Dyn
      {
         double zc, thx, thy, vz, wx, wy, ds, vds, dr, vr, iCoil, t;
         bool switchClosed, tipOutsideCone, asleep, collidable, pulseActive, pulseArmed, extOn, coilCommand, prevClosed, inRest, inCoilStop, restActive, coilActive;
         V3 gEff;
         double coilOnS, coilOnMaxS, wCoil, eDissRing, pulseLeftS, extVolts, extDuty, coilVolts, Qext[3], ringExt, restStored, coilStored;
         RubberState tip[2];
      };
      void SaveDyn(Dyn &d) const
      {
         d.zc = zc; d.thx = thx; d.thy = thy; d.vz = vz; d.wx = wx; d.wy = wy; d.ds = ds; d.vds = vds; d.dr = dr; d.vr = vr; d.iCoil = iCoil; d.t = t;
         d.switchClosed = switchClosed; d.tipOutsideCone = tipOutsideCone; d.asleep = asleep; d.collidable = collidable; d.pulseActive = pulseActive; d.pulseArmed = pulseArmed;
         d.extOn = extOn; d.coilCommand = coilCommand; d.prevClosed = prevClosed; d.inRest = m_inRest; d.inCoilStop = m_inCoilStop; d.restActive = m_restActive; d.coilActive = m_coilActive;
         d.gEff = gEff;
         d.coilOnS = coilOnS; d.coilOnMaxS = coilOnMaxS; d.wCoil = wCoil; d.eDissRing = eDissRing; d.pulseLeftS = pulseLeftS; d.extVolts = extVolts; d.extDuty = extDuty; d.coilVolts = coilVolts;
         d.Qext[0] = m_Qext[0]; d.Qext[1] = m_Qext[1]; d.Qext[2] = m_Qext[2]; d.ringExt = m_ringExt; d.restStored = m_restStored; d.coilStored = m_coilStored;
         d.tip[0] = m_tipState[0]; d.tip[1] = m_tipState[1];
      }
      void LoadDyn(const Dyn &d)
      {
         zc = d.zc; thx = d.thx; thy = d.thy; vz = d.vz; wx = d.wx; wy = d.wy; ds = d.ds; vds = d.vds; dr = d.dr; vr = d.vr; iCoil = d.iCoil; t = d.t;
         switchClosed = d.switchClosed; tipOutsideCone = d.tipOutsideCone; asleep = d.asleep; collidable = d.collidable; pulseActive = d.pulseActive; pulseArmed = d.pulseArmed;
         extOn = d.extOn; coilCommand = d.coilCommand; prevClosed = d.prevClosed; m_inRest = d.inRest; m_inCoilStop = d.inCoilStop; m_restActive = d.restActive; m_coilActive = d.coilActive;
         gEff = d.gEff;
         coilOnS = d.coilOnS; coilOnMaxS = d.coilOnMaxS; wCoil = d.wCoil; eDissRing = d.eDissRing; pulseLeftS = d.pulseLeftS; extVolts = d.extVolts; extDuty = d.extDuty; coilVolts = d.coilVolts;
         m_Qext[0] = d.Qext[0]; m_Qext[1] = d.Qext[1]; m_Qext[2] = d.Qext[2]; m_ringExt = d.ringExt; m_restStored = d.restStored; m_coilStored = d.coilStored;
         m_tipState[0] = d.tip[0]; m_tipState[1] = d.tip[1];
         events.clear();
      }

      // a force applied to the skirt at a body point (script, test rigs); applied at the next Step
      void PushSkirt(const V3 &bodyPoint, const V3 &F) { AddSkirtForce(bodyPoint, F, m_Qext); asleep = false; }
      void PushRing(double downN) { m_ringExt += downN; asleep = false; }

      // ---- one step of h seconds ----
      void Step(double h)
      {
         if (!collidable)
         {
            t += h;
            return;
         }
         StepSkirtAndSpoon(h);
         UpdateSwitch();
         DriveLogic(h);
         StepRingAndCoil(h);
         t += h;
      }

      // the cabinet acceleration (nudge): the skirt and spoon feel gravity minus it
      void SetCabinetAcceleration(const V3 &a) { SetGravity(V3(0.0, 0.0, -9.81) - a); }
      // the gravity the skirt and the spoon feel (gravity minus the cabinet acceleration, world frame, m/s^2); a change wakes the bumper
      void SetGravity(const V3 &g)
      {
         if (g.x != gEff.x || g.y != gEff.y || g.z != gEff.z)
            asleep = false;
         gEff = g;
      }

      double SkirtEdgeDropM(double azimuth) const   // drop of the rim point (top at the rim radius) below its rest height, along the given azimuth
      {
         const double rRim = 0.5 * p.skirtOdMm * MM - 0.5 * p.skirtTEdgeMm * MM;
         const V3 b(rRim * std::cos(azimuth), rRim * std::sin(azimuth), 0.0);
         return zRest - SkirtPoint(b).z;
      }

   private:
      double KappaEff(double R, const Revolve::Near &nr) const
      {
         double inv = 1.0 / R;
         if (nr.edge)
            inv += 1.0 / nr.rr;
         if (nr.rs > 1e-4)
            inv += nr.nr / nr.rs;
         return 1.0 / std::max(inv, 1.0 / (50.0 * R));
      }
      double m_Qext[3] = { 0.0, 0.0, 0.0 };
      double m_ringExt = 0.0;       // external force on the ring, down positive
      RubberState m_tipState[2];

      static void AddSkirtForceJ(const V3 col[3], const V3 &F, double Q[3])
      {
         Q[0] += col[0].Dot(F); Q[1] += col[1].Dot(F); Q[2] += col[2].Dot(F);
      }
      void AddSkirtForce(const V3 &b, const V3 &F, double Q[3]) const
      {
         V3 col[3];
         Jacobian(b, col);
         AddSkirtForceJ(col, F, Q);
      }

      void StepSkirtAndSpoon(double h)
      {
         double Q[3] = { m_Qext[0], m_Qext[1], m_Qext[2] };
         m_Qext[0] = m_Qext[1] = m_Qext[2] = 0.0;
         const double m = p.skirtMass;
         // weight (and the cabinet acceleration) at the centre of mass
         AddSkirtForce(V3(0.0, 0.0, cmZ), gEff * m, Q);
         // 8 springs on a circle under the skirt, each pushing up from a fixed rest height with the preload
         {
            const double a = p.springRadiusMm * MM, kk = p.skirtKNPerM / 8.0, f0 = p.skirtPreloadN / 8.0;
            for (int k = 0; k < 8; k++)
            {
               const V3 b(a * std::cos(k * M_PI / 4.0), a * std::sin(k * M_PI / 4.0), 0.0);
               const double pz = SkirtPoint(b).z;
               const double f = std::max(0.0, f0 + kk * (zRest - pz));
               AddSkirtForce(b, V3(0.0, 0.0, f), Q);
            }
         }
         // lip stop against the base rim, rim stop against the playfield (one-sided penalty spring and damper, 16 points each)
         {
            const double rLip = 0.5 * p.lipDiaMm * MM, tIn = p.skirtTInnerMm * MM, rRim = 0.5 * p.skirtOdMm * MM - 0.5 * p.skirtTEdgeMm * MM;
            const double kPt = p.skirtStopKNPerM / 16.0, cPt = 2.0 * p.skirtStopZeta * std::sqrt(p.skirtStopKNPerM * m) / 16.0;
            for (int k = 0; k < 16; k++)
            {
               const double cs = std::cos(k * M_PI / 8.0), sn = std::sin(k * M_PI / 8.0);
               {
                  const V3 b(rLip * cs, rLip * sn, tIn);
                  const double pen = SkirtPoint(b).z - zBase;
                  if (pen > 0.0)
                  {
                     const double f = std::max(0.0, kPt * pen + cPt * SkirtPointVelocity(b).z);
                     AddSkirtForce(b, V3(0.0, 0.0, -f), Q);
                  }
               }
               {
                  const V3 b(rRim * cs, rRim * sn, 0.0);
                  const double pz = SkirtPoint(b).z;
                  if (pz < 0.0)
                  {
                     const double f = std::max(0.0, kPt * (-pz) - cPt * SkirtPointVelocity(b).z);
                     AddSkirtForce(b, V3(0.0, 0.0, f), Q);
                  }
               }
            }
         }
         // stem tip on the spoon cone: two meridian contacts (the near and the far wall)
         double spoonDownN = 0.0;
         {
            const V3 bTip(0.0, 0.0, -(p.stemLenMm - p.tipRadiusMm) * MM);
            const V3 pt = SkirtPoint(bTip);
            const double rho = std::sqrt(pt.x * pt.x + pt.y * pt.y);
            const double ex = rho > 1e-9 ? pt.x / rho : 1.0, ey = rho > 1e-9 ? pt.y / rho : 0.0;
            const double rT = p.tipRadiusMm * MM, apex = coneApexRest - ds;
            const double coneR = p.coneRadiusMm * MM, coneD = p.coneDepthMm * MM;
            tipOutsideCone = tipOutsideCone || (rho > coneR - rT);
            BallSI tip; tip.R = rT; tip.m = m;
            const V3 vTip = SkirtPointVelocity(bTip), wS = SkirtOmega();
            for (int side = 0; side < 2; side++)
            {
               const double sg = side == 0 ? 1.0 : -1.0;
               // wall from the apex (0, apex) to the rim (sg * coneR, apex + coneD) in the (u, z) plane; the tip is at (rho, z)
               const double dr_ = sg * coneR, dz_ = coneD, len2 = dr_ * dr_ + dz_ * dz_;
               double tt = (rho * dr_ + (pt.z - apex) * dz_) / len2;
               tt = std::min(1.0, std::max(0.0, tt));
               const double cu = tt * dr_, cz = apex + tt * dz_;
               const double qu = rho - cu, qz = pt.z - cz, d = std::sqrt(qu * qu + qz * qz);
               if (d >= rT || d < 1e-12)
               {
                  m_tipState[side].Reset();
                  continue;
               }
               const double nu_ = qu / d, nz_ = qz / d;
               const V3 n(nu_ * ex, nu_ * ey, nz_);
               const double delta = rT - d;
               const V3 pc = pt + n * (-(rT - 0.5 * delta));
               const auto res = RubberContact::ContactCore(tipMat, tip, n, delta, rT, tipMat.thicknessMm * MM, V3(0.0, 0.0, -vds), vTip, wS, h, m_tipState[side], V3());
               AddSkirtForce(ToBody(pc), res.force, Q);
               spoonDownN += res.force.z;   // the spoon gets -F; down positive: +F.z
            }
         }
         // viscous damping on the generalised coordinates (ratio skirtZeta of the spring-set axial and tilt stiffness)
         {
            const double kTilt = p.skirtKNPerM * (p.springRadiusMm * MM) * (p.springRadiusMm * MM) * 0.5;
            Q[0] -= 2.0 * p.skirtZeta * std::sqrt(p.skirtKNPerM * m) * vz;
            Q[1] -= 2.0 * p.skirtZeta * std::sqrt(kTilt * mThetaEff) * wx;
            Q[2] -= 2.0 * p.skirtZeta * std::sqrt(kTilt * mThetaEff) * wy;
         }
         vz += Q[0] / m * h;
         wx += Q[1] / mThetaEff * h;
         wy += Q[2] / mThetaEff * h;
         zc += vz * h;
         thx += wx * h;
         thy += wy * h;

         // spoon (down positive): weight, spring, tip reaction, stops
         {
            const double ms = p.spoonMass;
            double f = ms * 9.81 + spoonDownN - (p.spoonPreloadN + p.spoonKNPerM * ds) - 2.0 * p.spoonZeta * std::sqrt(p.spoonKNPerM * ms) * vds;
            const double cStop = 2.0 * 0.3 * std::sqrt(p.spoonStopKNPerM * ms);
            if (ds < 0.0)
               f += std::max(0.0, p.spoonStopKNPerM * (-ds) - cStop * vds);
            const double dMax = p.spoonTravelMm * MM;
            if (ds > dMax)
               f -= std::max(0.0, p.spoonStopKNPerM * (ds - dMax) + cStop * vds);
            vds += f / ms * h;
            ds += vds * h;
         }
      }

      void UpdateSwitch()
      {
         const double close = p.switchCloseMm * MM, hyst = p.switchHystMm * MM;
         if (!switchClosed && ds >= close)
         {
            switchClosed = true;
            events.push_back({ t, EvHit, ds });   // Hit first, then the Switch event, then the drive reacts (DriveLogic)
            events.push_back({ t, EvSwitchClosed, ds });
         }
         else if (switchClosed && ds <= close - hyst)
         {
            switchClosed = false;
            events.push_back({ t, EvSwitchOpened, ds });
         }
      }

      void StartPulse()
      {
         pulseActive = true;
         pulseLeftS = p.pulseMs * 1e-3;
         pulseArmed = false;
         events.push_back({ t, EvPulseStart, 0.0 });
      }

      void DriveLogic(double h)
      {
         const bool closedEdge = switchClosed && !prevClosed;
         prevClosed = switchClosed;
         switch (p.mode)
         {
         case DriveMode::Held:
            coilCommand = switchClosed;
            coilVolts = p.coil.vSupply;
            break;
         case DriveMode::Pulse:
            if (pulseActive)
            {
               pulseLeftS -= h;
               if (pulseLeftS <= 0.0)
               {
                  pulseActive = false;
                  events.push_back({ t, EvPulseEnd, 0.0 });
               }
            }
            if (closedEdge && pulseArmed && !pulseActive)
               StartPulse();
            if (!pulseActive && !switchClosed)
               pulseArmed = true;
            coilCommand = pulseActive;
            coilVolts = p.coil.vSupply;
            break;
         case DriveMode::External:
            coilCommand = extOn;
            coilVolts = (extVolts > 0.0 ? extVolts : p.coil.vSupply) * extDuty;
            break;
         }
         if (coilCommand)
         {
            coilOnS += h;
            coilOnMaxS = std::max(coilOnMaxS, coilOnS);
         }
         else
            coilOnS = 0.0;
      }

      void StepRingAndCoil(double h)
      {
         const double stroke = p.ringStrokeMm * MM, bobbin = p.coil.bobbinMm * MM;
         const double frac = std::max(0.0, (stroke - dr) / bobbin + p.plungerClosedFrac);
         const double iNext = CoilModel::NextCurrent(p.coil, coilD, frac, std::max(vr, 0.0), coilCommand ? coilVolts : 0.0, iCoil, h);
         const double fCoil = CoilModel::Force(p.coil, coilD, frac, iCoil);   // down
         const double m = p.ringMass;
         double f = fCoil + m * 9.81 - (p.returnPreloadN + p.returnKNPerM * dr) - p.ringViscousNsPerM * vr + m_ringExt;
         double diss = p.ringViscousNsPerM * vr * vr;
         const double reff = p.ringStopReffMm * MM;
         if (dr < 0.0)   // rest stop (up): penetration -dr, closing rate -vr
         {
            const double dl = -dr, a = std::sqrt(reff * dl);
            const double fs = std::max(0.0, kRest * dl * std::sqrt(dl) + gRest * a * (-vr));
            f += fs;
            if (fs > 0.0) diss += gRest * a * vr * vr;
            // the viscoelastic force reaches zero before the compression does: the spring energy still stored then is not returned (lost, counted as dissipated)
            if (m_restActive && fs <= 0.0) eDissRing += 0.4 * kRest * dl * dl * std::sqrt(dl);
            m_restActive = fs > 0.0;
            m_restStored = m_restActive ? 0.4 * kRest * dl * dl * std::sqrt(dl) : 0.0;
            if (!m_inRest) { m_inRest = true; events.push_back({ t, EvRestStop, -vr }); }
         }
         else
         {
            m_inRest = false;
            m_restActive = false;
            m_restStored = 0.0;
         }
         if (dr > stroke)   // coil stop (down)
         {
            const double dl = dr - stroke, a = std::sqrt(reff * dl);
            const double fs = std::max(0.0, kCoil * dl * std::sqrt(dl) + gCoil * a * vr);
            f -= fs;
            if (fs > 0.0) diss += gCoil * a * vr * vr;
            if (m_coilActive && fs <= 0.0) eDissRing += 0.4 * kCoil * dl * dl * std::sqrt(dl);
            m_coilActive = fs > 0.0;
            m_coilStored = m_coilActive ? 0.4 * kCoil * dl * dl * std::sqrt(dl) : 0.0;
            if (!m_inCoilStop) { m_inCoilStop = true; events.push_back({ t, EvCoilStop, vr }); }
         }
         else
         {
            m_inCoilStop = false;
            m_coilActive = false;
            m_coilStored = 0.0;
         }
         wCoil += fCoil * vr * h;
         eDissRing += diss * h;
         m_ringExt = 0.0;
         vr += f / m * h;
         dr += vr * h;
         iCoil = iNext;
      }
      bool m_inRest = false, m_inCoilStop = false, m_restActive = false, m_coilActive = false;
      double m_restStored = 0.0, m_coilStored = 0.0;

   public:
      double StopStoredJ() const { return m_restStored + m_coilStored; }   // spring energy in the ring's stops while they are pushing (for the energy audit)
      // vertical gap between the ring's inner underside edge and the highest point of the skirt top under it (a check on demand, not part of Step)
      double RingSkirtGapM() const
      {
         const Revolve::Vtx &a = ring0.v[0];
         const double zUnder = a.z - dr - 0.5 * p.ringSheetMm * MM;
         double best = -1e9;
         for (int k = 0; k < 16; k++)
         {
            const double az = k * M_PI / 8.0;
            const V3 b = ToBody(V3(a.r * std::cos(az), a.r * std::sin(az), zc));
            const double rb = std::sqrt(b.x * b.x + b.y * b.y);
            const auto &v = skirt0.v;
            double zb;
            if (rb <= v[1].r) zb = v[1].z;
            else zb = v[1].z + std::min(1.0, (rb - v[1].r) / (v[2].r - v[1].r)) * (v[2].z - v[1].z);
            best = std::max(best, SkirtPoint(V3(b.x, b.y, zb)).z);
         }
         return zUnder - best;
      }
   };
}
