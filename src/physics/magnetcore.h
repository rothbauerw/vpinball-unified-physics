// license:GPLv3+
// Iron-core magnet: the core's answer to the ball (2026-10-18). Built at load per magnet, read per ball per cycle. No VPX dependencies (tested offline in runs/2026-10-14-mirror).
//
// Model (docs/NewPhysicsReference.txt, MAGNETS, core_model = cells):
//   core      = axisymmetric volume cells (4 x 30) of saturating steel (B-H M = Ms h / (h + Ms/chi)), solved by Newton for the coil alone and with one ball on the axis;
//   bare field = coil + coil-only core, read out smoothly into an (r, z) table above the face (bicubic Hermite);
//   no-reaction force = linear ball (mu 200) to order 4 in the bare field, measured on an inner sphere;
//   reaction  = on-axis reaction R(h) = M2 - M0 of the cell solve with the ball (order 12), carried off the axis by the equivalent height h* of the linear tube law
//               (P_axis(h*) = P(x, h)): down = R(h*) dh*/dz, centring = R(h*) dh*/dx;
//   ball saturation = sat/lin factor of the total from the linear gap field on the axis at h* (face_neck rule);
//   rim rule  = walking outward across the face edge (5.5 mm inside to 2.5 mm outside), the down force must lie between its corrected inner neighbour and the outer
//               trusted value, else linear interpolation (user rule, 2026-10-18).
// Everything is per operating level (coil current), interpolated linearly in the current at run time (0 A = no force).
#pragma once
#include <algorithm>
#include <string>
#include <vector>

namespace MagnetCore
{
   // Ball-ball magnetisation gate (REPORT.md "GATE", 2026-10-19): two touching balls in a field B along their line of centres pull 18080 N/T^2 x B^2 (linear FEM, the
   // largest orientation), so below 0.665 mT the pair force is under 8 mN (~1 % of the ball's weight). Each level's tables reach as far as its bare field is above it.
   constexpr double B_GATE = 0.665e-3;

   struct Geometry
   {
      double coreRadius = 9.5e-3;      // m
      double windLen = 0.038;          // winding length, m
      double exposed = 6e-3;           // core length above the winding (to the face), m
      double windIn = 10.5e-3;         // winding inner radius, m
      double windOut = 16.5e-3;        // winding outer radius, m
      double turns = 2000;
      double bSat = 1.6;               // core steel saturation, T
      double ballR = 25 * 0.5397e-3;   // ball radius the tables are built for, m
      bool rimRule = true;             // apply the rim rule to the down force (off: for checks)
   };

   struct BareTable   // H (A/m) of coil + coil-only core above the face plane: r from the axis, z above the face
   {
      std::vector<double> r, z, F[2], Fr[2], Fz[2], Frz[2];   // Hr, Hz and their spline slopes
      void Build();                                           // slopes from F
      void Eval(double rho, double zz, double v[2], double vr[2], double vz[2]) const;
      bool Inside(double rho, double zz) const { return !r.empty() && rho <= r.back() && zz >= z.front() && zz <= z.back(); }
   };

   struct Level
   {
      double amps = 0;
      double reach = 0;                               // m from the face centre: the bare field at this level is below B_GATE beyond it (rays at 0..90 deg from the axis)
      BareTable bare;
      std::vector<double> hR, R, sat;                 // on-axis: gap (m), reaction down force (N, linear ball), ball-saturation factor
      std::vector<double> down, side;                 // force table over the geometry grid (N): toward the face, toward the axis (rim rule applied to down)
      // for a ball-ball model (used by the multipole model, removed 2026-10-19; kept for its successor): the core with ONE ball on the face (plug, on the axis, gap hR.front()) and the plug's own drive per on-axis gap
      BareTable plug;                                 // H (A/m) of coil + core in the plug state (the plug's own field NOT included), same layout as bare
      std::vector<std::vector<double>> plugB;         // per gap (as hR): scaled zonal drive b_bar_n0 (orthonormal Y_n0, n = 1..13) of coil + core at the plug centre, core frozen
   };

   struct Model
   {
      Geometry g;
      std::vector<double> xg, hg, hs, dhx, dhz;       // geometry table: offset x, gap h (m), equivalent height h* and slopes; out to the furthest bare table (x 0.5 mm steps to 30 mm, 1 mm beyond; h log)
      std::vector<Level> levels;                      // ascending current
      std::string log;                                // build report (timings, checks)

      // Build for the given operating currents (A); returns false (with log) if nothing could be built
      bool Build(const Geometry &geom, const std::vector<double> &amps);
      double Reach() const { double r = 0; for (const Level &L : levels) r = std::max(r, L.reach); return r; }   // the top level's reach (m)
      // Force on a ball of radius g.ballR centred at offset x (m, >= 0) and gap h (m, above the face) at coil current I (A): down (toward the face) and side (toward the axis), N.
      // Returns false outside the bare table (no force).
      bool Force(double x, double h, double I, double &down, double &side) const;
      // Bare field B (T) at (x, height above the face): Br along the offset, Bz along the axis, interpolated in the current
      bool Field(double x, double z, double I, double &Br, double &Bz) const;
   };
}
