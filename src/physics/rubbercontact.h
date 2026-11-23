// license:GPLv3+
//
// Ball / rubber contact force model (unified engine).
//
// Pure SI maths with no VPX types, so the same code runs inside VPX and in the standalone test harness
// (runs/2026-10-04-rubber/rubber_harness.cpp). It is a direct port of the contact loop of
// runs/2026-10-04-rubber/rubber_contact.py; the two must stay in step. See DESIGN_engine_integration.md.
//
// Model (all constants are in RubberParams and meant to be edited; the defaults are unverified guesses):
//   geometry  the flipper outline (base circle, end circle, two tangent faces) is the OUTER SURFACE of the rubber, extruded
//             between zLo and zHi above the playfield. The contact point is the nearest point of that surface to the
//             ball centre, so the normal is fully 3D (it is horizontal while the ball centre is inside [zLo, zHi]).
//   patch     a = sqrt(Reff * delta); Reff = R*sqrt(r/(R+r)) on a curved surface of horizontal radius r, R on a flat face
//             (sphere on a vertical cylinder, circular-patch approximation of the elliptical contact)
//   normal    F = (4/3) E* a^3/Reff * (1 + beta (a/h)^3) * (1 + 1.5 alpha ddot), never pulling
//             Hertz + thin bonded layer stiffening + Hunt-Crossley damping (rebound ratio ~ 1 - alpha v)
//   shear     Mindlin spring kt = 8 G* a times shearScale (the sideways stiffness: < 1 where the rubber itself gives way, a tube on a post, a strip on a bat), a damper
//             shearTan sqrt(kt mt), spring and damper together within mu Fn. Mindlin-Deresiewicz: a patch that grows holds the force, one that shrinks holds the stretch.
//             The stretch xi is kept in RubberState and re-seated at the cap on the spring's force alone
//   twist     drilling-friction torque (3 pi/16) mu a Fn opposing spin about the contact normal
//   rolling   hysteresis rolling resistance: a ball rolling on rubber compresses the rubber ahead and the rubber behind
//             recovers with a lag, which moves the pressure centre ahead of the ball centre. Torque on the ball
//             T = mu_rr R Fn opposing the ball's rolling relative to the bat, with the equal and opposite torque on the bat.
//             mu_rr = 2.34 (p_a / E) tan(delta), p_a = Fn / (pi a^2)   (Persson, arXiv 1008.4493, hard sphere on a
//             viscoelastic solid at low speed). Evaluated every step from the contact state, so the only input is the
//             rubber's loss tangent. tan(delta) is PROVISIONAL (a typical elastomer, 0.1-0.3; not measured for flipper
//             rubber) and is taken constant; it really depends on frequency (~ v/a, 100-150 Hz here). 0 switches it off.
//             The roll axis is the vertical one (the way a ball rolls along a vertical rubber face), and the rolling
//             rate is the ball's spin relative to the bat's own rotation.
//
// Coordinates: the caller passes everything relative to the flipper pivot (pivot at the origin, z up from the playfield),
// angles in radians, forces in N, torques in N m, lengths in m.

#pragma once

#include <cmath>
#include <algorithm>
#include <vector>
#include <array>
#include <map>
#include <utility>
#include <cstring>
#include <cstdint>
#include <cstdlib>
#include <atomic>
#include <string>
#include "contactprof.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// The reach of a ball over one cycle (to-do N5, the one reach rule; shared finder design.md section 3): how far its centre can move in dt from its speed and its
// acceleration, reach = speed x dt + 1/2 x accel x dt^2 (m; speed m/s, accel m/s^2, dt s). Used by the shared finder's per-cycle query (to-do 1.34 P2a) and the
// group step's part finder lists (P6a). Here, not in contactsurface.h, so the group step needs no VPX types (ModelStandards 2e).
inline double ContactReachM(const double speed, const double accel, const double dt) { return speed * dt + 0.5 * accel * dt * dt; }
// the energy reach (runs/2026-11-13-finder-P6/design_reach.md sections 1, 4b; to-do N5's one reach rule): how far a body (a ball, a part's turn or travel) can get in T.
// Energy enters only through the drives and what feeds it: 1/2 m v^2 <= 1/2 m v0^2 + E + m a x, so dx/dt <= sqrt(c^2 + 2 a x), c^2 = v0^2 + 2 E / m:
// x <= c T + 1/2 a T^2, v_max = c + a T (m and v0 in the body's own units: kg and m/s, or kg m^2 and rad/s; E J; a the largest drive acceleration). E 0: the reach rule above
inline double EnergyReach(const double v0, const double E, const double m, const double a, const double T, double *const vMax = nullptr)
{
   const double c = std::sqrt(v0 * v0 + 2.0 * std::max(E, 0.0) / m);
   if (vMax != nullptr)
      *vMax = c + a * T;
   return c * T + 0.5 * a * T * T;
}
// the drive-capped move (design_reach.md section 1; one law for every part's EnergyMove): the energy reach above, and at most what the body makes in T with all its energy as speed:
// its kinetic energy, the feed E and its drives' work bound D (J, the drives' largest work over the move): sqrt(v0^2 + 2 (E + D) / m) T, plus add (a term the caller's law adds, e.g.
// gravity's 1/2 g T^2 on a part it pulls beyond its drives' bound; 0 none). a: the largest drive acceleration (+infinity: none known, the cap alone decides)
inline double DriveCappedReach(const double v0, const double E, const double m, const double a, const double D, const double T, const double add = 0.0)
{
   return std::min(EnergyReach(v0, E, m, a, T), std::sqrt(v0 * v0 + 2.0 * (std::max(E, 0.0) + D) / m) * T + add);
}

namespace RubberContact
{
   struct V3
   {
      double x = 0.0, y = 0.0, z = 0.0;
      V3() = default;
      V3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
      V3 operator+(const V3 &o) const { return V3(x + o.x, y + o.y, z + o.z); }
      V3 operator-(const V3 &o) const { return V3(x - o.x, y - o.y, z - o.z); }
      V3 operator*(double s) const { return V3(x * s, y * s, z * s); }
      V3 operator/(double s) const { return V3(x / s, y / s, z / s); }
      double Dot(const V3 &o) const { return x * o.x + y * o.y + z * o.z; }
      V3 Cross(const V3 &o) const { return V3(y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x); }
      double Length() const { return std::sqrt(Dot(*this)); }
   };

   struct RubberParams
   {
      double eMpa = 12.0;            // effective Young's modulus (rubber_e_mpa)
      double nu = 0.5;               // Poisson ratio: fixed, solid rubber is ~incompressible (not an ini setting)
      double thicknessMm = 3.2;      // rubber layer thickness
      double alphaSPerM = 0.05;      // normal damping; rebound ratio ~ 1 - alpha*v
      double gammaKK = 0.0;          // > 0: viscoelastic (Kuwabara-Kono) normal damping instead of alpha: Fd = gamma a v_n (a = contact radius), N s / m^2. The rebound ratio then falls only as v^(1/5), close to a constant restitution; see CalibrateKK
      double betaLayer = 0.33;       // thin-layer stiffening coefficient (approximation, unchecked)
      double mu = 0.8;               // rubber / steel friction (black / white rubber; user, 2026-10-30)
      double shearScale = 1.0;       // the sideways stiffness, a multiplier on Mindlin's 8 G* a. 1 = a thick solid body (steel, plastic, wood). Rubbers: [Rubbers] sideways_stiffness, [Flippers] rubber_sideways_stiffness
      double shearTan = 0.4;         // loss tangent of the sideways spring: damper c = shearTan sqrt(kt mt), a damping ratio of shearTan / 2. A rubber's is its own tan_delta (the same material; for nu 0.5, E = 3 G at every frequency). This default is for contacts that are not rubber (steel, plastic, wood; the old damping ratio 0.2, a guess): it stands in for their micro-slip
      bool twist = true;             // drilling friction on/off
      double tanDelta = 0.0;         // loss tangent for the rolling resistance (0 = off); see the header comment
      double rollingMu = 0.0;        // floor of the rolling-resistance coefficient of a static surface (the playfield: [Contact] rolling_resistance); 0 = none
      bool normalLoss = false;       // tanDelta also damps the normal force: a dashpot c = tanDelta sqrt(k m) (damping ratio tanDelta / 2 at the contact's own natural frequency), so a ball at rest settles; off for the flipper bat
      double zLoMm = 0.36 * 25.4;    // rubber bottom above the playfield surface (adjustable)
      double heightMm = 0.44 * 25.4; // rubber height
   };

   struct BallSI
   {
      double R = 25.0 * 0.53975e-3;  // radius, m
      double m = 0.080;              // mass, kg
      double Inertia() const { return 0.4 * m * R * R; }
   };

   // Flipper outline in the horizontal plane. rb: base circle radius (at the pivot), re: end circle radius, Lc: centre-to-centre.
   struct Outline
   {
      double rb = 0.0, re = 0.0, Lc = 0.0;

      struct Nearest
      {
         double px, py;       // nearest point on the outline
         double nx, ny;       // outward unit normal there
         double radius;       // horizontal radius of curvature: rb, re, or INFINITY on a face
      };

      Nearest NearestPoint(double qx, double qy, double phi) const
      {
         const double nu = (rb - re) / Lc;
         const double nv = std::sqrt(std::max(0.0, 1.0 - nu * nu));
         const double ux = std::cos(phi), uy = std::sin(phi);
         const double vx = -std::sin(phi), vy = std::cos(phi);
         const double ex = Lc * ux, ey = Lc * uy;
         Nearest best{};
         double bestD = INFINITY;
         auto consider = [&](double px, double py, double nx, double ny, double rad) {
            const double d = std::sqrt((qx - px) * (qx - px) + (qy - py) * (qy - py));
            if (d < bestD)
            {
               bestD = d;
               best = Nearest{ px, py, nx, ny, rad };
            }
         };

         // two faces
         for (int k = 0; k < 2; k++)
         {
            const double side = (k == 0) ? 1.0 : -1.0;
            const double nx = nu * ux + side * nv * vx, ny = nu * uy + side * nv * vy;
            const double ax = rb * nx, ay = rb * ny;
            const double bx = ex + re * nx, by = ey + re * ny;
            const double tx = bx - ax, ty = by - ay;
            const double L2 = tx * tx + ty * ty;
            const double s = std::min(std::max(((qx - ax) * tx + (qy - ay) * ty) / L2, 0.0), 1.0);
            consider(ax + s * tx, ay + s * ty, nx, ny, INFINITY);
         }
         // base circle: valid where the outward normal m satisfies m.u <= nu
         {
            const double dx = qx, dy = qy;
            const double dn = std::sqrt(dx * dx + dy * dy);
            if (dn > 1e-12)
            {
               const double mx = dx / dn, my = dy / dn;
               if (mx * ux + my * uy <= nu)
                  consider(rb * mx, rb * my, mx, my, rb);
               else
                  for (int k = 0; k < 2; k++)
                  {
                     const double side = (k == 0) ? 1.0 : -1.0;
                     const double nx = nu * ux + side * nv * vx, ny = nu * uy + side * nv * vy;
                     consider(rb * nx, rb * ny, nx, ny, rb);
                  }
            }
         }
         // end circle: valid where m.u >= nu
         {
            const double dx = qx - ex, dy = qy - ey;
            const double dn = std::sqrt(dx * dx + dy * dy);
            if (dn > 1e-12)
            {
               const double mx = dx / dn, my = dy / dn;
               if (mx * ux + my * uy >= nu)
                  consider(ex + re * mx, ey + re * my, mx, my, re);
               else
                  for (int k = 0; k < 2; k++)
                  {
                     const double side = (k == 0) ? 1.0 : -1.0;
                     const double nx = nu * ux + side * nv * vx, ny = nu * uy + side * nv * vy;
                     consider(ex + re * nx, ey + re * ny, nx, ny, re);
                  }
            }
         }
         return best;
      }
   };

   // Hertz effective radius of a sphere (radius R) on a vertical cylinder of horizontal radius r (INFINITY = flat face)
   inline double EffectiveRadius(double R, double r)
   {
      if (!std::isfinite(r))
         return R;
      const double Rx = 1.0 / (1.0 / R + 1.0 / r);
      return std::sqrt(Rx * R);
   }

   struct Probe
   {
      bool overlap = false;  // the ball (radius R) penetrates the rubber surface
      V3 point;              // nearest point on the rubber surface (outline extruded over [zLo, zHi])
      double dist = 0.0;     // distance from the ball centre to that point
      double delta = 0.0;    // penetration depth R - dist (positive = overlapping), m
      double reff = 0.0;     // Hertz effective radius at that point, m
   };

   // Pure geometry, no state: nearest point of the 3D rubber surface to the ball centre c (pivot frame, m).
   // Works for balls that already overlap (delta > 0); the outline is the outer surface of the rubber, extruded
   // between zLo and zHi, so the normal tilts above and below the rubber and is horizontal in between.
   inline Probe ProbeSurface(const RubberParams &rp, const BallSI &ball, const Outline &ol, double phi, const V3 &c)
   {
      const double zLo = rp.zLoMm * 1e-3, zHi = (rp.zLoMm + rp.heightMm) * 1e-3;
      const Outline::Nearest nn = ol.NearestPoint(c.x, c.y, phi);
      Probe pr;
      pr.point = V3(nn.px, nn.py, std::min(std::max(c.z, zLo), zHi));
      pr.reff = EffectiveRadius(ball.R, nn.radius);
      pr.dist = (c - pr.point).Length();
      pr.delta = ball.R - pr.dist;
      pr.overlap = pr.delta > 0.0;
      return pr;
   }

   // Upper bound on the rolling-resistance coefficient. Persson's expression is a small-pressure limit (p_a << E); this keeps
   // it from running away if a very soft rubber or a hard hit pushes p_a/E up.
   constexpr double kRollingMuMax = 0.3;

   // The rolling and twist stick springs' angles at a contact (group step, [Contact] roll_model 1): th in the contact plane (rolling), tw about the normal (twist), rad
   struct RollAngle
   {
      V3 th;
      double tw = 0.0;
   };

   struct RubberState
   {
      V3 xi;                 // shear spring state (tangent plane), m; zero when not in contact
      RollAngle ra;          // rolling and twist springs (group step, roll_model 1); zero when not in contact
      void Reset() { xi = V3(); ra = RollAngle(); }
   };

   struct ContactResult
   {
      bool inContact = false;
      V3 force;              // force on the ball, N
      V3 torqueOnBall;       // torque on the ball about its centre, N m
      double torqueOnBat = 0.0; // torque on the bat about the vertical axis through the pivot, N m (closing = +phi direction)
      V3 normal;             // contact normal, pointing from the rubber surface to the ball centre
      V3 point;              // nearest point on the rubber surface
      double delta = 0.0;    // compression, m
      double patchA = 0.0;   // contact patch radius, m
      double fn = 0.0;       // normal force magnitude, N
      double reff = 0.0;     // effective radius used, m
      double dist = 0.0;     // distance from the ball centre to the rubber surface point, m (set even when not in contact)
   };

   struct CoreResult
   {
      V3 force;              // force on the ball, N
      V3 torque;             // torque on the ball about its centre (normal, shear, twist), N m
      V3 rB;                 // contact point relative to the ball centre
      double fn = 0.0;       // normal force magnitude, N
      double a = 0.0;        // contact patch radius, m
      double E = 0.0;        // Young's modulus, Pa
   };

   // The surface-independent part of the contact law: normal (Hertz + thin layer + Hunt-Crossley), Mindlin shear with the Coulomb cap, twist.
   // The caller supplies the geometry of the touch: n = outward surface normal (surface -> ball centre), delta = compression, reff = Hertz
   // effective radius, h = rubber layer thickness (m), vs = velocity of the surface point. Shared by the flipper (Evaluate) and static rubbers.
   inline CoreResult ContactCore(const RubberParams &rp, const BallSI &ball, const V3 &n, double delta, double reff, double h,
                                 const V3 &vs, const V3 &v, const V3 &w, double dt, RubberState &st, const V3 &ws = V3())
   {
      CoreResult out;
      const double R = ball.R, m = ball.m, Ib = ball.Inertia();
      const double E = rp.eMpa * 1e6;
      const double Estar = E / (1.0 - rp.nu * rp.nu);
      const double G = E / (2.0 * (1.0 + rp.nu));
      const double Gstar = G / (2.0 - rp.nu);

      const V3 rB = n * (-(R - delta / 2.0));
      const V3 vrel = v + w.Cross(rB) - vs;
      const double vn = vrel.Dot(n);
      const V3 vt = vrel - n * vn;

      const double a = std::sqrt(reff * delta);
      const double ah = a / h;
      const double fl = 1.0 + rp.betaLayer * (ah * ah * ah);
      const double Fe = (4.0 / 3.0) * Estar * (a * a * a / reff) * fl;
      // Impact damping (alpha, Hunt-Crossley: grows with speed) plus the rubber's loss tangent: a constant loss factor means a damping ratio tanDelta / 2 at whatever
      // frequency the contact rings, i.e. a dashpot c = tanDelta k / omega_n = tanDelta sqrt(k m) with the tangent stiffness k = 1.5 Fe / delta. It damps small motions at
      // any speed and any ball mass (alpha alone does nothing for a ball at rest).
      double Fd = rp.gammaKK > 0.0 ? rp.gammaKK * a * (-vn) : Fe * 1.5 * rp.alphaSPerM * (-vn);
      if (rp.normalLoss && rp.tanDelta > 0.0)
         Fd += rp.tanDelta * std::sqrt(1.5 * Fe * m / std::max(delta, 1e-9)) * (-vn);
      const double Fn = std::max(0.0, Fe + Fd);

      const double kt = 8.0 * Gstar * a * rp.shearScale;
      // Mindlin-Deresiewicz: a patch that grows takes in unstressed rubber at its rim, so the force holds and the stretch falls; a patch that shrinks keeps the stretch. What
      // the spring holds can then only fall. kt goes as sqrt(delta), and the depth a step ago comes from the closing speed, so no state is kept for it.
      if (vn < 0.0)
         st.xi = st.xi * std::sqrt(std::max(delta + vn * dt, 0.0) / delta);
      V3 xiOld;
      {
         const double l0 = st.xi.Length();
         const V3 xp = st.xi - n * st.xi.Dot(n);
         const double lp = xp.Length();
         xiOld = lp > 1e-300 ? xp * (l0 / lp) : V3();   // turned into the contact plane, its length kept
         st.xi = xiOld + vt * dt;
      }
      const double mt = 1.0 / (1.0 / m + R * R / Ib);          // tangential effective mass (m/3.5 for a solid sphere)
      const double ct = rp.shearTan * std::sqrt(kt * mt);
      // the patch: the spring kt and the damper ct side by side, in series with the Coulomb slip (the group step's law): the trial force is both, past the limit the
      // patch slips along it and the spring's end stretch is the one that carries the limit force with the damper (Ft = -kt xi - ct (xi - xiOld) / dt)
      const double lim = rp.mu * Fn;
      V3 Ft = st.xi * (-kt) + vt * (-ct);
      const double nt = Ft.Length();
      if (nt > lim)
      {
         Ft = Ft * (lim / nt);
         st.xi = (xiOld * (ct / dt) - Ft) * (1.0 / (kt + ct / dt));
      }

      const V3 F = n * Fn + Ft;
      V3 Mb = rB.Cross(F);
      if (rp.twist)
      {
         const double wn = (w - ws).Dot(n);   // the spin about the normal relative to the surface (a turning disc drives the ball's spin)
         Mb = Mb + n * (-(3.0 * M_PI / 16.0) * rp.mu * a * Fn * std::tanh(wn / 2.0));
      }
      out.force = F;
      out.torque = Mb;
      out.rB = rB;
      out.fn = Fn;
      out.a = a;
      out.E = E;
      return out;
   }

   // The restitution of a normal impact at speed v0 of a ball (radius R, mass m) on a flat surface with the Hertz law of ContactCore and viscoelastic damping gamma (no layer, no friction).
   inline double ImpactRestitutionKK(const RubberParams &rp, double R, double m, double gamma, double v0)
   {
      const double Estar = rp.eMpa * 1e6 / (1.0 - rp.nu * rp.nu);
      const double tc = 2.9 * std::pow(m * m / (R * Estar * Estar * v0), 0.2);   // Hertz contact time
      const double dt = tc / 4000.0;
      double d = 0.0, vd = v0;   // compression and its rate (approach speed positive)
      for (int i = 0; i < 40000; i++)
      {
         const double a = std::sqrt(R * std::max(d, 0.0));
         const double Fe = (4.0 / 3.0) * Estar * std::pow(std::max(d, 0.0), 1.5) * std::sqrt(R);
         double Fd = gamma * a * vd;
         if (rp.normalLoss && rp.tanDelta > 0.0 && d > 0.0)   // the loss tangent's dashpot is part of the contact too
            Fd += rp.tanDelta * std::sqrt(1.5 * Fe * m / std::max(d, 1e-9)) * vd;
         const double Fn = std::max(0.0, Fe + Fd);
         vd -= Fn / m * dt;
         d += vd * dt;
         if (d <= 0.0 && i > 0)
            return std::max(0.0, -vd) / v0;
      }
      return 0.0;
   }

   // The gamma that gives the restitution eTarget at speed vRef (bisection in log gamma; the rebound falls with gamma). eTarget >= 0.995 gives 0 (an elastic impact).
   inline double CalibrateKK(const RubberParams &rp, double R, double m, double vRef, double eTarget)
   {
      if (eTarget >= 0.995)
         return 0.0;
      double lo = 1.0, hi = 1e9;
      for (int i = 0; i < 48; i++)
      {
         const double mid = std::sqrt(lo * hi);
         if (ImpactRestitutionKK(rp, R, m, mid, vRef) > eTarget)
            lo = mid;
         else
            hi = mid;
      }
      return std::sqrt(lo * hi);
   }

   // Evaluate the contact for one time step of length dt. Updates the shear state.
   //   phi, omBat  bat angle (rad, from the +x axis, counter-clockwise) and angular velocity about +z (rad/s)
   //   c, v, w     ball centre, velocity and angular velocity relative to the pivot frame
   inline ContactResult Evaluate(const RubberParams &rp, const BallSI &ball, const Outline &ol, double phi, double omBat,
                                 const V3 &c, const V3 &v, const V3 &w, double dt, RubberState &st)
   {
      ContactResult out;
      const double R = ball.R;
      const double h = rp.thicknessMm * 1e-3;

      const Probe pr = ProbeSurface(rp, ball, ol, phi, c);
      out.dist = pr.dist;
      if (!pr.overlap)
      {
         st.Reset();
         return out;
      }
      const V3 P = pr.point;
      const double reff = pr.reff, delta = pr.delta, dist = pr.dist;
      const V3 n = (c - P) / dist;
      const V3 ez(0.0, 0.0, 1.0);
      const V3 vs = ez.Cross(P) * omBat;

      const CoreResult core = ContactCore(rp, ball, n, delta, reff, h, vs, v, w, dt, st);
      const V3 F = core.force;
      V3 Mb = core.torque;
      const double Fn = core.fn, a = core.a, E = core.E;

      double rollTorqueZ = 0.0;      // rolling-resistance torque on the ball about z (the reaction acts on the bat)
      if (rp.tanDelta > 0.0)
      {
         // roll axis: vertical, projected into the tangent plane (skipped where the surface is nearly flat, a ball on top of the rubber)
         const V3 et = ez - n * ez.Dot(n);
         const double etLen = et.Length();
         if (etLen > 0.2)
         {
            const V3 e = et / etLen;
            const double wRel = (w - ez * omBat).Dot(e);                         // rolling rate relative to the bat
            const double pa = Fn / (M_PI * a * a);                               // mean contact pressure
            const double muRr = std::min(2.34 * (pa / E) * rp.tanDelta, kRollingMuMax);
            const double T = -muRr * R * Fn * std::tanh(wRel / 2.0);            // opposes the rolling; smooth through zero
            Mb = Mb + e * T;
            rollTorqueZ = T * e.z;
         }
      }

      out.inContact = true;
      out.force = F;
      out.torqueOnBall = Mb;
      out.torqueOnBat = ez.Dot(P.Cross(F * -1.0)) - rollTorqueZ;
      out.normal = n;
      out.point = P;
      out.delta = delta;
      out.patchA = a;
      out.fn = Fn;
      out.reff = reff;
      return out;
   }

   // ---- Static rubber objects (posts, wall rubbers): a tube of circular cross-section swept along a polyline centreline ----

   struct TubeProbe
   {
      bool overlap = false;
      V3 point;              // nearest point on the centreline
      V3 normal;             // surface -> ball centre (unit)
      double dist = 0.0;     // ball centre to centreline, m
      double delta = 0.0;    // R + rho - dist, m
      double reff = 0.0;     // Hertz effective radius, m
      double rho = 0.0;      // tube radius, m
      int seg = -1;          // segment index of the nearest point
      double grad = 1.0;     // ProbePatch: the length of the depth's gradient along normal (below 1 where the patch spans a fold)
   };

   // Two parts of one surface pressing a ball (the arms of a fold, a concave corner of a mesh) are two Hertz contacts: their stored energies add, as one depth d = (dA^2.5 + dB^2.5)^0.4 with gradient
   // (dA^1.5 nA + dB^1.5 nB) / d^1.5. The share of a part fades as its nearest point nears its end (fit: how much of the patch is on it), so it leaves with no jump.
   struct FoldSum
   {
      double s25 = 0.0;
      V3 g;
      void Add(const double d, const double fit, const V3 &n)
      {
         const double w = fit * fit * std::sqrt(fit);
         const double d15 = d * std::sqrt(d);
         s25 += w * d * d15;
         g = g + n * (w * d15);
      }
      void Finish(TubeProbe &pr, const double zero) const
      {
         const double d = std::pow(s25, 0.4);
         const V3 gr = g / (d * std::sqrt(d));
         const double gl = gr.Length();
         pr.delta = d + zero;
         if (gl > 1e-6)
         {
            pr.normal = gr / gl;
            pr.grad = gl;
         }
      }
   };

   // Surface curvature along the tube at the contact (signed, 1/m, > 0 convex), the torus's curvature along its parallel (Meusnier): for a centreline of curvature k,
   // cos/(1/k + rho cos), cosang the 3D cosine between the contact's unit normal and the centreline's in-plane normal pointing away from its centre of curvature
   // (+1 the outer equator, -1 the inner one, 0 the top: a normal along the binormal sees no bend along the tube).
   inline double TubeAlongCurvature(double k, double rho, double cosang)
   {
      if (k < 1e-9)
         return 0.0;
      return cosang / (1.0 / k + rho * cosang);
   }

   // What a static centreline gives for free on every probe: the bend at each vertex and the segment bounds. Built once (BuildTubePre); ProbeTube(..., &pre) then returns
   // exactly what ProbeTube without it does, minus the acos / sqrt per probe and with most segments rejected by their bounding sphere.
   struct TubePre
   {
      std::vector<V3> t;                 // segment i: pts[i+1] - pts[i]
      std::vector<double> L2;            // |t|^2 (0 = degenerate, skipped)
      std::vector<V3> mid;               // bounding sphere of the segment
      std::vector<double> rad;
      std::vector<V3> gmid;              // bounding sphere of each group of kGroup segments
      std::vector<double> grad;
      std::vector<double> kap;           // turning curvature at vertex j (0 at the ends)
      std::vector<double> dnx, dny;      // horizontal unit direction of the change of heading at vertex j
      std::vector<char> dnOk;            // dn long enough to normalise
      static constexpr int kGroup = 8;
      mutable int hint = 0;              // segment that was nearest last time: its distance is the first bound
   };

   inline double TubeKappa(const V3 *pts, int n, int j)
   {
      if (j < 1 || j > n - 2)
         return 0.0;
      const V3 t0 = pts[j] - pts[j - 1], t1 = pts[j + 1] - pts[j];
      const double l0 = t0.Length(), l1 = t1.Length();
      if (l0 <= 0.0 || l1 <= 0.0)
         return 0.0;
      return std::acos(std::min(std::max(t0.Dot(t1) / (l0 * l1), -1.0), 1.0)) / (0.5 * (l0 + l1));
   }

   inline void BuildTubePre(TubePre &pre, const V3 *pts, int n)
   {
      const int ns = std::max(0, n - 1);
      pre.t.assign(ns, V3());
      pre.L2.assign(ns, 0.0);
      pre.mid.assign(ns, V3());
      pre.rad.assign(ns, 0.0);
      for (int i = 0; i < ns; i++)
      {
         pre.t[i] = pts[i + 1] - pts[i];
         pre.L2[i] = pre.t[i].Dot(pre.t[i]);
         pre.mid[i] = pts[i] + pre.t[i] * 0.5;
         pre.rad[i] = 0.5 * std::sqrt(pre.L2[i]);
      }
      const int ng = (ns + TubePre::kGroup - 1) / TubePre::kGroup;
      pre.gmid.assign(ng, V3());
      pre.grad.assign(ng, 0.0);
      for (int g = 0; g < ng; g++)
      {
         const int i0 = g * TubePre::kGroup, i1 = std::min(ns, i0 + TubePre::kGroup);   // vertices i0 .. i1
         V3 sum;
         for (int i = i0; i <= i1; i++)
            sum = sum + pts[i];
         pre.gmid[g] = sum / (double)(i1 - i0 + 1);
         double r2 = 0.0;
         for (int i = i0; i <= i1; i++)
            r2 = std::max(r2, (pts[i] - pre.gmid[g]).Dot(pts[i] - pre.gmid[g]));
         pre.grad[g] = std::sqrt(r2);
      }
      pre.kap.assign(n, 0.0);
      pre.dnx.assign(n, 0.0);
      pre.dny.assign(n, 0.0);
      pre.dnOk.assign(n, 0);
      for (int j = 0; j < n; j++)
      {
         pre.kap[j] = TubeKappa(pts, n, j);
         if (j >= 1 && j <= n - 2 && pre.kap[j] > 1e-9)
         {
            const V3 t0 = pts[j] - pts[j - 1], t1 = pts[j + 1] - pts[j];
            V3 dn = t1 / t1.Length() - t0 / t0.Length();
            dn.z = 0.0;
            const double dl = dn.Length();
            if (dl > 1e-12)
            {
               pre.dnx[j] = dn.x / dl;
               pre.dny[j] = dn.y / dl;
               pre.dnOk[j] = 1;
            }
         }
      }
      pre.hint = 0;
   }

   // pts: n+1 centreline vertices (pivot-free world frame, m). Tilted rubbers work because the points are taken after the rubber transform.
   // pre: optional, built from the same points (a static centreline); the result is the same.
   // segs (shared_finder P5, needs pre): the segments to read, in index order (the ball's query); the nearest of them, the same arithmetic and tie rule (the first in index order)
   inline TubeProbe ProbeTube(const BallSI &ball, const V3 *pts, int n, double rho, const V3 &c, const TubePre *pre = nullptr, const std::vector<int> *segs = nullptr)
   {
      TubeProbe pr;
      pr.rho = rho;
      double bestD2 = INFINITY;   // squared: the nearest segment is found without a square root per segment
      V3 bestP;
      int bestI = -1;
      if (pre != nullptr && segs != nullptr)
      {
         for (const int i : *segs)
         {
            const double L2 = pre->L2[i];
            if (L2 <= 0.0)
               continue;
            const V3 dm = c - pre->mid[i];
            const double lim = std::sqrt(bestD2) * (1.0 + 1e-9) + 1e-12 + pre->rad[i];
            if (dm.Dot(dm) > lim * lim)
               continue;
            const V3 &t = pre->t[i];
            const double s = std::min(std::max((c - pts[i]).Dot(t) / L2, 0.0), 1.0);
            const V3 p = pts[i] + t * s;
            const V3 dv = c - p;
            const double d2 = dv.Dot(dv);
            if (d2 < bestD2)
            {
               bestD2 = d2;
               bestP = p;
               bestI = i;
            }
         }
      }
      else if (pre != nullptr && n >= 2)
      {
         // upper bound on the nearest distance from last time's nearest segment; a segment whose bounding sphere is farther cannot be the nearest
         double ubD = INFINITY;
         const int h = std::min(pre->hint, n - 2);
         if (pre->L2[h] > 0.0)
         {
            const double s = std::min(std::max((c - pts[h]).Dot(pre->t[h]) / pre->L2[h], 0.0), 1.0);
            const V3 dv = c - (pts[h] + pre->t[h] * s);
            ubD = std::sqrt(dv.Dot(dv)) * (1.0 + 1e-9) + 1e-12;
         }
         for (int i = 0; i + 1 < n; i++)
         {
            if (i % TubePre::kGroup == 0)
            {
               const int g = i / TubePre::kGroup;
               const V3 dg = c - pre->gmid[g];
               const double glim = std::min(ubD, std::sqrt(bestD2) * (1.0 + 1e-9) + 1e-12) + pre->grad[g] * (1.0 + 1e-9);
               if (dg.Dot(dg) > glim * glim)
               {
                  i += TubePre::kGroup - 1;   // the whole group is farther than the best so far
                  continue;
               }
            }
            const double L2 = pre->L2[i];
            if (L2 <= 0.0)
               continue;
            const V3 dm = c - pre->mid[i];
            const double lim = std::min(ubD, std::sqrt(bestD2) * (1.0 + 1e-9) + 1e-12) + pre->rad[i];
            if (dm.Dot(dm) > lim * lim)
               continue;
            const V3 &t = pre->t[i];
            const double s = std::min(std::max((c - pts[i]).Dot(t) / L2, 0.0), 1.0);
            const V3 p = pts[i] + t * s;
            const V3 dv = c - p;
            const double d2 = dv.Dot(dv);
            if (d2 < bestD2)
            {
               bestD2 = d2;
               bestP = p;
               bestI = i;
            }
         }
         if (bestI >= 0)
            pre->hint = bestI;
      }
      else
      for (int i = 0; i + 1 < n; i++)
      {
         const V3 t = pts[i + 1] - pts[i];
         const double L2 = t.Dot(t);
         if (L2 <= 0.0)
            continue;
         const double s = std::min(std::max((c - pts[i]).Dot(t) / L2, 0.0), 1.0);
         const V3 p = pts[i] + t * s;
         const V3 dv = c - p;
         const double d2 = dv.Dot(dv);
         if (d2 < bestD2)
         {
            bestD2 = d2;
            bestP = p;
            bestI = i;
         }
      }
      const double bestD = std::sqrt(bestD2);
      if (bestI < 0 || bestD <= 1e-12)
         return pr;
      pr.point = bestP;
      pr.seg = bestI;
      pr.dist = bestD;
      pr.normal = (c - bestP) / bestD;
      pr.delta = ball.R + rho - bestD;
      pr.overlap = pr.delta > 0.0;

      // turning curvature around the segment (mean of its two end vertices) and where the normal sits relative to the centre of curvature
      auto kap = [&](int j) { return pre != nullptr ? pre->kap[j] : TubeKappa(pts, n, j); };
      const double k = 0.5 * (kap(bestI) + kap(bestI + 1));
      double kLong = 0.0;
      if (k > 1e-9)
      {
         const int j = std::min(std::max(bestI, 1), n - 2);
         double dnx, dny;
         bool ok;
         if (pre != nullptr)
         {
            dnx = pre->dnx[j];
            dny = pre->dny[j];
            ok = pre->dnOk[j] != 0;
         }
         else
         {
            const V3 t0 = pts[j] - pts[j - 1], t1 = pts[j + 1] - pts[j];
            V3 dn = t1 / t1.Length() - t0 / t0.Length();
            dn.z = 0.0;
            const double dl = dn.Length();
            ok = dl > 1e-12;
            dnx = ok ? dn.x / dl : 0.0;
            dny = ok ? dn.y / dl : 0.0;
         }
         if (ok)   // (the 3D cosine between the unit normal and the horizontal unit dn: Meusnier, TubeAlongCurvature)
            kLong = TubeAlongCurvature(k, rho, -(pr.normal.x * dnx + pr.normal.y * dny));
      }
      const double R1 = 1.0 / (1.0 / ball.R + 1.0 / rho);
      // concave bends tighter than the ball would give a zero or negative radius; keep the along-tube radius within 20 ball radii
      const double R2 = 1.0 / std::max(1.0 / ball.R + kLong, 0.05 / ball.R);
      pr.reff = std::sqrt(R1 * R2);
      return pr;
   }


   // ProbeTube's result with the other parts of the same tube pressing the ball: a nearest point in the middle of a segment (a local minimum of the distance along the tube)
   // segs (shared_finder P5): the segments to read, in index order (null: all)
   inline void PatchTube(const BallSI &ball, const V3 *pts, const int n, const double rho, const V3 &c, const TubePre &pre, const double zero, TubeProbe &pr, const std::vector<int> *segs = nullptr)
   {
      if (!pr.overlap || pr.seg < 0 || pr.reff <= 0.0)
         return;
      if (pr.delta <= zero)
         return;
      const double S = ball.R + rho;
      const double a0 = std::max(std::sqrt(pr.reff * (pr.delta - zero)), 1e-6);
      const auto fitOn = [&](const int i, const double s) { return std::min(1.0, std::min(s, 1.0 - s) * std::sqrt(pre.L2[i]) / a0); };   // how much of the patch is on this segment
      const double sb = pre.L2[pr.seg] > 0.0 ? (c - pts[pr.seg]).Dot(pre.t[pr.seg]) / pre.L2[pr.seg] : 0.0;
      const double fitBest = fitOn(pr.seg, std::min(std::max(sb, 0.0), 1.0));
      FoldSum sum;
      bool any = false;
      const int cnt = segs != nullptr ? (int)segs->size() : n - 1;
      for (int k = 0; k < cnt; k++)
      {
         const int i = segs != nullptr ? (*segs)[k] : k;
         if (i == pr.seg || pre.L2[i] <= 0.0)
            continue;
         const V3 dm = c - pre.mid[i];
         const double lim = S + pre.rad[i];
         if (dm.Dot(dm) > lim * lim)
            continue;
         const double s = (c - pts[i]).Dot(pre.t[i]) / pre.L2[i];
         if (s <= 0.0 || s >= 1.0)
            continue;
         const V3 dv = c - (pts[i] + pre.t[i] * s);
         const double dist = dv.Length();
         if (dist <= 1e-12 || S - dist <= zero)
            continue;
         const double fit = std::min(fitBest, fitOn(i, s));
         if (fit <= 0.0 || (dv / dist - pr.normal).Length() < 1e-9)
            continue;
         if (!any)
         {
            sum.Add(pr.delta - zero, 1.0, pr.normal);
            any = true;
         }
         sum.Add(S - dist - zero, fit, dv / dist);
      }
      if (any)
         sum.Finish(pr, zero);
   }

   // ---- Rubber faces: walls and primitives flagged as rubber. A triangle soup (world frame, metres); the ball meets a face, an edge or a corner. ----

   struct SurfaceTri
   {
      V3 a, b, c;            // vertices
      V3 n;                  // unit normal of the face (either orientation: the ball side is taken from where the ball is)
      V3 lo, hi;             // bounding box
      bool area;             // false: zero area (three points on a line), no face of its own
   };

   inline SurfaceTri MakeSurfaceTri(const V3 &a, const V3 &b, const V3 &c)
   {
      SurfaceTri t;
      t.a = a;
      t.b = b;
      t.c = c;
      const V3 cr = (b - a).Cross(c - a);
      const double l = cr.Length();
      t.area = l > 0.0;
      t.n = t.area ? cr / l : V3(0.0, 0.0, 1.0);
      t.lo = V3(std::min({ a.x, b.x, c.x }), std::min({ a.y, b.y, c.y }), std::min({ a.z, b.z, c.z }));
      t.hi = V3(std::max({ a.x, b.x, c.x }), std::max({ a.y, b.y, c.y }), std::max({ a.z, b.z, c.z }));
      return t;
   }

   // Closest point of a triangle to p (Ericson, Real-Time Collision Detection 5.1.5). kind: 0 face, 1 edge, 2 vertex.
   inline V3 ClosestOnTriangle(const SurfaceTri &t, const V3 &p, int &kind)
   {
      const V3 ab = t.b - t.a, ac = t.c - t.a, ap = p - t.a;
      const double d1 = ab.Dot(ap), d2 = ac.Dot(ap);
      kind = 2;
      if (d1 <= 0.0 && d2 <= 0.0)
         return t.a;
      const V3 bp = p - t.b;
      const double d3 = ab.Dot(bp), d4 = ac.Dot(bp);
      if (d3 >= 0.0 && d4 <= d3)
         return t.b;
      const double vc = d1 * d4 - d3 * d2;
      if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
      {
         kind = 1;
         return t.a + ab * (d1 / (d1 - d3));
      }
      const V3 cp = p - t.c;
      const double d5 = ab.Dot(cp), d6 = ac.Dot(cp);
      if (d6 >= 0.0 && d5 <= d6)
         return t.c;
      const double vb = d5 * d2 - d1 * d6;
      if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
      {
         kind = 1;
         return t.a + ac * (d2 / (d2 - d6));
      }
      const double va = d3 * d6 - d5 * d4;
      if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
      {
         kind = 1;
         return t.b + (t.c - t.b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
      }
      kind = 0;
      const double denom = 1.0 / (va + vb + vc);
      return t.a + ab * (vb * denom) + ac * (vc * denom);
   }

   // ---- The path kit (the path test; runs/2026-11-13-finder-exactpath/design.md section 2): exact squared distances from a ball centre's straight move in a step,
   // the segment s + u d (u in [0, 1]; the solver's own model of the move), to a straight piece. One place, shared by every shape that has straight pieces.
   // to a point (a ball pair: the relative move against the other centre)
   inline double PathPointD2(const V3 &s, const V3 &d, const V3 &p)
   {
      const double dd = d.Dot(d);
      const double u = dd > 0.0 ? std::min(std::max((p - s).Dot(d) / dd, 0.0), 1.0) : 0.0;
      const V3 q = s + d * u - p;
      return q.Dot(q);
   }
   // to a segment p2 + w d2 (an edge, a tube's centreline segment): closest points of two segments (Ericson, Real-Time Collision Detection 5.1.9)
   inline double PathSegD2(const V3 &p1, const V3 &d1, const V3 &p2, const V3 &d2)
   {
      const V3 r = p1 - p2;
      const double a = d1.Dot(d1), e = d2.Dot(d2), f = d2.Dot(r);
      double s = 0.0, t = 0.0;
      if (a <= 0.0 && e <= 0.0)
         return r.Dot(r);
      if (a <= 0.0)
         t = std::min(std::max(f / e, 0.0), 1.0);
      else
      {
         const double c = d1.Dot(r);
         if (e <= 0.0)
            s = std::min(std::max(-c / a, 0.0), 1.0);
         else
         {
            const double b = d1.Dot(d2), den = a * e - b * b;
            s = den > 0.0 ? std::min(std::max((b * f - c * e) / den, 0.0), 1.0) : 0.0;   // (parallel: any s; 0, then t is clamped below)
            t = (b * s + f) / e;
            if (t < 0.0)
            {
               t = 0.0;
               s = std::min(std::max(-c / a, 0.0), 1.0);
            }
            else if (t > 1.0)
            {
               t = 1.0;
               s = std::min(std::max((b - c) / a, 0.0), 1.0);
            }
         }
      }
      const V3 q = p1 + d1 * s - (p2 + d2 * t);
      return q.Dot(q);
   }
   // to a triangle (a mesh face with its edges and corners, as ProbeTriangles measures it): the segment and the triangle are convex, so apart from a crossing the nearest
   // pair has the segment's end or a triangle edge in it: both ends (ClosestOnTriangle, the probe's own function), the three edges, and the point where the segment crosses the face's plane
   inline double PathTriD2(const V3 &s, const V3 &d, const SurfaceTri &t)
   {
      int kind = 0;
      const V3 e = s + d;
      V3 q = s - ClosestOnTriangle(t, s, kind);
      double best = q.Dot(q);
      q = e - ClosestOnTriangle(t, e, kind);
      best = std::min(best, q.Dot(q));
      best = std::min(best, PathSegD2(s, d, t.a, t.b - t.a));
      best = std::min(best, PathSegD2(s, d, t.b, t.c - t.b));
      best = std::min(best, PathSegD2(s, d, t.c, t.a - t.c));
      if (t.area)
      {
         const double ds = (s - t.a).Dot(t.n), de = (e - t.a).Dot(t.n);
         if (ds != de && ((ds <= 0.0 && de >= 0.0) || (ds >= 0.0 && de <= 0.0)))
         {
            const V3 x = s + d * (ds / (ds - de));
            q = x - ClosestOnTriangle(t, x, kind);
            best = std::min(best, q.Dot(q));
         }
      }
      return best;
   }
   // within reach: the piece's depth (reach - distance, the probe's form) is above 0 somewhere on the path. The band 1e-9 x reach only widens it (a rounding guard against the
   // probes' own arithmetic at the ends: a re-run more, never one less)
   inline bool PathWithin(const double d2, const double reach) { return reach - std::sqrt(d2) > -1e-9 * reach; }
   // a box against the path's box grown by reach (lo, hi of the path's ends)
   inline bool PathBoxApart(const V3 &lo, const V3 &hi, const V3 &plo, const V3 &phi, const double reach)
   {
      return lo.x > phi.x + reach || hi.x < plo.x - reach || lo.y > phi.y + reach || hi.y < plo.y - reach || lo.z > phi.z + reach || hi.z < plo.z - reach;
   }
   // counters of the path test (VPX_CHECK_PROFILE=1 only; diagnostics, logged with the finder's once-a-second totals): skipped contacts tested by the path, by the two probes (a curved or moving
   // shape); cleared at the object box (a), at every piece box (b), by the exact distance (c); exact kit tests; reached; FarStillFar calls with a skipped contact, their time (ns) and worst
   // Build 2, per kind of curved shape (0 post, 1 bumper body, 2 disc): tested, reached, reached only inside the cover's sagitta band (the cover's distance
   // within reach + sagitta but not within reach: an upper bound on the re-runs the cover adds, the true curve not measured);
   // pieces passed over by their box (boxSkip) against the kit tests above (kit)
   struct PathStats
   {
      uint64_t tested = 0, twoPoint = 0, clearA = 0, clearB = 0, clearC = 0, kit = 0, reached = 0, reachedTwo = 0;
      uint64_t calls = 0, ns = 0, nsWorst = 0;
      uint64_t kTested[3] = {}, kReached[3] = {}, kBand[3] = {}, boxSkip = 0;
   };
   inline PathStats pathStats;
   // The covers' self-check (VPX_FINDER_CHECK=1; diagnostics only, no bit changes; groupstep.h PathSelfCheck): per kind as PathStats, the checks made and the disagreements, and the first
   // disagreement's numbers for the log's WARN (logged with "Path test totals")
   inline bool PathCheckOn()
   {
      static const bool on = [] { const char *e = std::getenv("VPX_FINDER_CHECK"); return e != nullptr && e[0] == '1'; }();
      return on;
   }
   struct PathCheck
   {
      uint64_t n[3] = {}, bad[3] = {};
      bool warned[3] = {}, logged[3] = {};
      std::string first[3];
   };
   inline PathCheck pathCheck;

   // ---- Covers of curved pieces (the path test, build 2; design.md sections 2 and 2b) ----
   // A curved centre line or surface is replaced, for the path test only, by straight chords or flat facets inscribed in it, built once with the shape. Each lies within its sagitta s of the true
   // curve and the curve within s of it (a two-way bound), so the path's distance to the curve >= its distance to the cover - s: a cover farther than reach + s clears the piece exactly.
   // kPathSagittaM: the largest sagitta of a cover (m). A COST setting, not physics: fewer chords mean fewer kit tests and more re-runs for a path that passes inside the band (a re-run only
   // puts the contact in the solve; it adds no force unless the ball touches). 20 um: about 32 chords on a 4 mm ring (design.md 2b, paper/chords.txt).
   constexpr double kPathSagittaM = 20e-6;
   // the chords of a circle of radius r that keep the sagitta r (1 - cos(pi / N)) at or under kPathSagittaM (at least 8)
   inline int PathChords(const double r)
   {
      if (r <= kPathSagittaM)
         return 8;
      return std::max(8, (int)std::ceil(3.14159265358979323846 / std::acos(1.0 - kPathSagittaM / r)));
   }
   // A profile in (rho, z) about the vertical axis through o, revolved (a post ring: one vertex, its circle; a bumper body or a disc: a polyline, each segment a band of flat facets, a
   // frustum or a flat ring). Each vertex has a sweep radius rr (the surface lies rr from the revolved line; a post's tube, a bumper sheet's half thickness). Reach of a segment: R + its larger rr.
   struct PathRevolve
   {
      struct Seg
      {
         int first = 0, count = 0;   // its pieces: chords (pts[first + k] to pts[first + k + 1]) for a ring, else facets (tris[first + k])
         bool ring = false;
         double rr = 0.0, sag = 0.0;   // the larger sweep radius, the sagitta at the larger radius
         V3 lo, hi;                    // the pieces' box
      };
      std::vector<Seg> seg;
      std::vector<V3> pts;
      std::vector<SurfaceTri> tris;
      double MaxSag() const
      {
         double s = 0.0;
         for (const Seg &g : seg)
            s = std::max(s, g.sag);
         return s;
      }
      void Build(const V3 &o, const double *r, const double *z, const double *rr, const int n)
      {
         seg.clear();
         pts.clear();
         tris.clear();
         const auto at = [&](const double rv, const double zv, const double ph) { return V3(o.x + rv * std::cos(ph), o.y + rv * std::sin(ph), o.z + zv); };
         const auto boxOf = [](Seg &g, const V3 &p) {
            g.lo = V3(std::min(g.lo.x, p.x), std::min(g.lo.y, p.y), std::min(g.lo.z, p.z));
            g.hi = V3(std::max(g.hi.x, p.x), std::max(g.hi.y, p.y), std::max(g.hi.z, p.z));
         };
         const double twoPi = 2.0 * 3.14159265358979323846;
         if (n == 1)   // a ring: N chords of the circle
         {
            Seg g;
            g.ring = true;
            g.rr = rr[0];
            const int N = PathChords(r[0]);
            g.sag = r[0] * (1.0 - std::cos(0.5 * twoPi / N));
            g.first = 0;
            g.count = N;
            g.lo = g.hi = at(r[0], z[0], 0.0);
            for (int k = 0; k <= N; k++)
            {
               pts.push_back(at(r[0], z[0], twoPi * (k % N) / N));
               boxOf(g, pts.back());
            }
            seg.push_back(g);
            return;
         }
         for (int i = 0; i + 1 < n; i++)   // a band of N facets per segment, each a flat trapezoid (two parallel chords at the two ends' heights) as two triangles
         {
            Seg g;
            const double rMax = std::max(r[i], r[i + 1]);
            const int N = PathChords(rMax);
            g.rr = std::max(rr[i], rr[i + 1]);
            g.sag = rMax * (1.0 - std::cos(0.5 * twoPi / N));
            g.first = (int)tris.size();
            g.lo = g.hi = at(r[i], z[i], 0.0);
            for (int k = 0; k < N; k++)
            {
               const double p0 = twoPi * k / N, p1 = twoPi * ((k + 1) % N) / N;
               const V3 a0 = at(r[i], z[i], p0), a1 = at(r[i], z[i], p1), b0 = at(r[i + 1], z[i + 1], p0), b1 = at(r[i + 1], z[i + 1], p1);
               if (r[i] > 0.0)   // (an end on the axis: its chord is a point, one triangle)
                  tris.push_back(MakeSurfaceTri(a0, a1, b1));
               if (r[i + 1] > 0.0)
                  tris.push_back(MakeSurfaceTri(a0, b1, b0));
               for (const V3 &p : { a0, a1, b0, b1 })
                  boxOf(g, p);
            }
            g.count = (int)tris.size() - g.first;
            seg.push_back(g);
         }
      }
      // 1 = a piece of the cover comes within its reach + sagitta of the path s -> s + d (a re-run), 0 = none (exact: the true surface is clear). band: reached only inside the sagitta band.
      int Reach(const V3 &s, const V3 &d, const double R, const bool prof, bool &band) const
      {
         PathStats &ps = pathStats;
         band = false;
         const V3 e = s + d;
         const V3 plo(std::min(s.x, e.x), std::min(s.y, e.y), std::min(s.z, e.z)), phi(std::max(s.x, e.x), std::max(s.y, e.y), std::max(s.z, e.z));
         for (const Seg &g : seg)
         {
            const double reach = R + g.rr, lim = reach + g.sag;
            if (PathBoxApart(g.lo, g.hi, plo, phi, lim))
            {
               if (prof)
                  ps.boxSkip += (uint64_t)g.count;
               continue;
            }
            for (int k = 0; k < g.count; k++)
            {
               double d2;
               if (g.ring)
               {
                  const V3 &a = pts[g.first + k], &b = pts[g.first + k + 1];
                  if (PathBoxApart(V3(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)), V3(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)), plo, phi, lim))
                  {
                     if (prof)
                        ps.boxSkip++;
                     continue;
                  }
                  d2 = PathSegD2(s, d, a, b - a);
               }
               else
               {
                  const SurfaceTri &t = tris[g.first + k];
                  if (PathBoxApart(t.lo, t.hi, plo, phi, lim))
                  {
                     if (prof)
                        ps.boxSkip++;
                     continue;
                  }
                  d2 = PathTriD2(s, d, t);
               }
               if (prof)
                  ps.kit++;
               if (PathWithin(d2, reach))
               {
                  band = false;
                  return 1;
               }
               if (PathWithin(d2, lim))
               {
                  band = true;
                  if (!prof)
                     return 1;   // (profile only: look on for a piece within reach itself, so the band count is not inflated; the answer is 1 either way)
               }
            }
         }
         return band ? 1 : 0;
      }
   };

   // Which triangles can touch a ball of radius R centred within `reach` of c (indices into tris)
   inline void SurfaceCandidates(const std::vector<SurfaceTri> &tris, const V3 &c, double reach, std::vector<int> &out)
   {
      out.clear();
      for (size_t i = 0; i < tris.size(); i++)
      {
         const SurfaceTri &t = tris[i];
         if (c.x + reach < t.lo.x || c.x - reach > t.hi.x || c.y + reach < t.lo.y || c.y - reach > t.hi.y || c.z + reach < t.lo.z || c.z - reach > t.hi.z)
            continue;
         out.push_back((int)i);
      }
   }

   // The nearest point on the faces in `cand`. A contact on a face is flat (effective radius R). On an edge or a corner, where the direction to the ball is
   // more than ~5 deg off the face normal, the rubber counts as rounded to edgeRho (the effective radius falls toward that of a ball on a tube of radius
   // edgeRho, fully at 15 deg), so a knife edge has a contact patch. The depth is measured from the sharp geometry. The result has rho = 0.

   // One candidate face as the probe measured it: its nearest point to the ball centre, the kind of that point, the squared distance.
   struct FaceHit
   {
      int i;
      int kind;
      V3 p;
      double d2;
   };

   // shared_finder P3a (to-do 1.34; runs/2026-11-11-finder/design.md section 4, R6): one measurement per face per ball position. One table per ball and static mesh object, keyed by the
   // ball centre (exact bits) and the ball's query version (a re-query or a rebind starts a new table). A slot per position in the object's query list holds the face's closest point,
   // its kind and the squared distance, filled on first request, so no face loop measures a face it did not measure before. check: every read of a filled slot also measures fresh and
   // counts any bit difference (diagnostics only, VPX_FINDER_CHECK=1).
   struct FaceTable
   {
      struct Slot
      {
         V3 p;
         double d2;
         int kind;
         uint32_t gen;
      };
      std::vector<Slot> s;
      uint32_t gen = 0;   // the table's generation: a slot is filled when its gen equals it
      uint64_t ver = ~0ull;
      double cx = 0.0, cy = 0.0, cz = 0.0;
      bool check = false;
      uint64_t fills = 0, reads = 0, mism = 0;   // faces measured, table reads, check-mode mismatches (summed by the owner)
      int lastMismJ = -1;
      uint64_t cyc = ~0ull;   // the counters: the cycle of the last contact check and the checks in it (contactsurfacehost.cpp TabCheck)
      uint32_t nCyc = 0;
      bool count = false;    // P3a's per-read counters (fills, reads): check or profile mode only
      int tag = -1;          // the counters: 1 = a rubber corner timed apart (design.md section 7), 0 = not, -1 = not yet known
      // shared_finder P3b (design.md section 4; list pass runs/2026-11-12-listpass/brief.md): the faces within R, found once per table (ball position) by one pass over the query list.
      // inJ / inI: the faces within R of c, in load order (query-list positions, face indices); useNear: the loops read them (false: the whole query list, check mode's comparison).
      uint32_t nearGen = 0;
      bool useNear = false;
      std::vector<int> inJ, inI;
      // the query list's index for the edge slots (EdgeSlot), rebuilt once per query version: idx by face id gives the face's list position (stamped, never cleared)
      struct Idx
      {
         uint32_t stamp;
         int j;
      };
      std::vector<Idx> idx;
      uint32_t idxStamp = 0;
      bool idxOk = false;
      // P3c: the faces' edge measurements (BarsTriangles), slot 3 j + e for edge e of the face at list position j, filled on first request at this ball position. A shared edge has one slot,
      // its lower face id's (when that face is in the query list): both faces measure it from its lower welded-vertex end, and welded vertices are the same floats, so the bits are the same.
      // check: every read of a filled slot also measures fresh and counts any bit difference (emism).
      struct ESlot
      {
         V3 n;          // unit direction from the edge's nearest point to the ball centre (0 when the centre lies on the edge)
         double dist;   // the centre's distance to the edge
         uint32_t gen;
      };
      std::vector<ESlot> es;
      uint64_t efills = 0, ereads = 0, emism = 0;   // edges measured, edge reads (every request), check-mode mismatches (counted in check or profile mode)
      uint64_t cbMism = 0;   // cut B check mode: calls the rule would skip whose full law ended with more than one peak (summed by the owner)
      int ListPos(const int i) const   // a face's query-list position from the index (-1: no index for this query list, or not in it)
      {
         return idxOk && i >= 0 && i < (int)idx.size() && idx[i].stamp == idxStamp ? idx[i].j : -1;
      }
      ESlot *EdgeSlot(const int t, const int e, const std::vector<std::array<int, 3>> &nb)   // null: no slot (measure directly)
      {
         const int j = ListPos(t);
         if (j < 0 || 3 * (size_t)j + 2 >= es.size())
            return nullptr;
         const int o = nb[t][e];
         if (o >= 0 && o < t)   // shared with a lower face: its slot
         {
            const int jo = ListPos(o);
            if (jo >= 0)
               for (int f = 0; f < 3; f++)
                  if (nb[o][f] == t)
                     return &es[3 * (size_t)jo + f];
         }
         return &es[3 * (size_t)j + e];
      }
      void *objStats = nullptr;   // P3b2 counters (check or profile mode): the owner's per-object entry, looked up once
      uint64_t nearMiss = 0, nearBox = 0, nearVisit = 0;   // check mode: faces within R the pass did not hand over; box tests and faces past the box (measured through the table) (counted with fills)
      int lastMissJ = -1;
      static double BoxD2(const SurfaceTri &t, const V3 &c)   // squared distance from c to the face's box: a lower bound of its distance
      {
         const double dx = std::max(std::max(t.lo.x - c.x, c.x - t.hi.x), 0.0), dy = std::max(std::max(t.lo.y - c.y, c.y - t.hi.y), 0.0), dz = std::max(std::max(t.lo.z - c.z, c.z - t.hi.z), 0.0);
         return dx * dx + dy * dy + dz * dz;
      }
      static constexpr double kBoxSlack = 1e-9;   // relative: the box bound and a measured d2 are each rounded, so a face is dropped only when clearly farther (a rounding guard, never a cutoff)
      Slot &Get(const int j, const SurfaceTri &t, const V3 &c)   // the slot, measured on first request
      {
         Slot &x = s[j];
         if (x.gen != gen)
         {
            x.p = ClosestOnTriangle(t, c, x.kind);
            const V3 dv = c - x.p;
            x.d2 = dv.Dot(dv);
            x.gen = gen;
            if (count)
               fills++;
         }
         return x;
      }
      // the faces within R: one pass over the query list in list order (load order); a face is dropped when its box is farther than R (an exact lower bound of its distance), else
      // measured through the table. list: the object's query list (face indices, ascending). nb: only for the edge slots' index (built when the faces are linked).
      void NearPass(const std::vector<SurfaceTri> &tris, const std::vector<std::array<int, 3>> &nb, const std::vector<int> &list, const V3 &c, const double R2)
      {
         if (nearGen == gen)
            return;
         nearGen = gen;
         inJ.clear();
         if (!idxOk && nb.size() == tris.size())   // a new query list: index it (once per query version)
         {
            idxOk = true;
            if (++idxStamp == 0)
            {
               for (Idx &x : idx)
                  x.stamp = 0;
               idxStamp = 1;
            }
            if (idx.size() < tris.size())
               idx.resize(tris.size(), Idx { 0, -1 });
            for (size_t j = 0; j < list.size(); j++)
               idx[list[j]] = Idx { idxStamp, (int)j };
         }
         const double boxR2 = R2 * (1.0 + kBoxSlack);
         for (size_t j = 0; j < list.size(); j++)
         {
            const int i = list[j];
            if (count)
               nearBox++;
            if (BoxD2(tris[i], c) > boxR2)
               continue;
            if (count)
               nearVisit++;
            if (Get((int)j, tris[i], c).d2 < R2)
               inJ.push_back((int)j);
         }
         inI.resize(inJ.size());
         for (size_t k = 0; k < inJ.size(); k++)
            inI[k] = list[inJ[k]];
         if (check)   // every face of the list measured fresh: one within R that the pass did not hand over is a miss
            for (size_t j = 0; j < list.size(); j++)
            {
               int k = 0;
               const V3 dv = c - ClosestOnTriangle(tris[list[j]], c, k);
               if (dv.Dot(dv) < R2 && !std::binary_search(inJ.begin(), inJ.end(), (int)j))
               {
                  nearMiss++;
                  lastMissJ = (int)j;
               }
            }
      }
      void Begin(const V3 &c, const uint64_t v, const size_t n)
      {
         if (gen == 0 || v != ver || std::memcmp(&c.x, &cx, sizeof(double)) != 0 || std::memcmp(&c.y, &cy, sizeof(double)) != 0 || std::memcmp(&c.z, &cz, sizeof(double)) != 0)
         {
            if (++gen == 0)   // wrapped: every slot empty again
            {
               for (Slot &x : s)
                  x.gen = 0;
               for (ESlot &x : es)
                  x.gen = 0;
               gen = 1;
            }
            if (v != ver)
               idxOk = false;   // P3b2: a new query list
            ver = v;
            cx = c.x;
            cy = c.y;
            cz = c.z;
         }
         if (s.size() < n)
            s.resize(n, Slot { V3(), 0.0, 0, 0 });
         if (es.size() < 3 * n)
            es.resize(3 * n, ESlot { V3(), 0.0, 0 });
      }
      const Slot &At(const int j, const SurfaceTri &t, const V3 &c)
      {
         Slot &x = s[j];
         if (count)
            reads++;
         if (x.gen != gen)
            Get(j, t, c);
         else if (check)
         {
            int k = 0;
            const V3 p = ClosestOnTriangle(t, c, k);
            const V3 dv = c - p;
            const double d2 = dv.Dot(dv);
            if (std::memcmp(&p.x, &x.p.x, sizeof(double)) != 0 || std::memcmp(&p.y, &x.p.y, sizeof(double)) != 0 || std::memcmp(&p.z, &x.p.z, sizeof(double)) != 0 || k != x.kind
                || std::memcmp(&d2, &x.d2, sizeof(double)) != 0)
            {
               mism++;
               lastMismJ = j;
            }
         }
         return x;
      }
   };
   // What a face loop reads: the table (null: measure as before) and the query-list positions of the loop's faces (null: the loop's own index, the loop runs over the query list itself)
   struct FaceRead
   {
      FaceTable *t = nullptr;
      const int *pos = nullptr;
      bool prune = false;   // P3b: the probe with no face within R: the nearest of the whole list, a face skipped when its box is strictly farther than the best so far (exact lower bound)
      mutable FaceTable::Slot tmp {};
      const FaceTable::Slot &At(const size_t j, const SurfaceTri &tri, const V3 &c) const
      {
         if (t == nullptr)   // no table: measured directly
         {
            tmp.p = ClosestOnTriangle(tri, c, tmp.kind);
            const V3 dv = c - tmp.p;
            tmp.d2 = dv.Dot(dv);
            return tmp;
         }
         return t->At(pos != nullptr ? pos[j] : (int)j, tri, c);
      }
   };

   // halfThick > 0: the faces are the mid-surface of a sheet of thickness 2 halfThick (a gate plate): the surface is everything within halfThick of the faces
   // (flat on both sides, edges rounded to radius halfThick), so the depth is R + halfThick - distance. Pass edgeRho = halfThick for such a sheet.
   // hits (optional): every candidate's measurement, in candidate order, for a caller that needs the other faces too (BarsTriangles), so each face is measured once.
   inline TubeProbe ProbeTriangles(const BallSI &ball, const std::vector<SurfaceTri> &tris, const std::vector<int> &cand, double edgeRho, const V3 &c,
                                   double halfThick = 0.0, std::vector<FaceHit> *hits = nullptr, const FaceRead *fr = nullptr)
   {
      TubeProbe pr;
      double bestD2 = INFINITY;
      V3 bestP;
      int bestI = -1, bestKind = 0;
      { const ContactProf::LoopTicks profLoop(cand.size());   // P3a profile (diagnostics only)
      for (size_t j = 0; j < cand.size(); j++)
      {
         const int i = cand[j];
         if (fr != nullptr && fr->prune && bestI >= 0 && FaceTable::BoxD2(tris[i], c) > (bestD2 + 1e-18) * (1.0 + FaceTable::kBoxSlack))
            continue;   // P3b: cannot be nearer than the best so far nor tie with it (exact lower bound; the 1e-18 tie band kept)
         int kind = 0;
         V3 p;
         const FaceTable::Slot *sl = nullptr;
         if (fr != nullptr)   // shared_finder P3a: the table's measurement
         {
            sl = &fr->At(j, tris[i], c);
            p = sl->p;
            kind = sl->kind;
         }
         else
            p = ClosestOnTriangle(tris[i], c, kind);
         const V3 dv = c - p;
         const double d2 = sl != nullptr ? sl->d2 : dv.Dot(dv);
         if (hits != nullptr)
            hits->push_back({ i, kind, p, d2 });
         // on a tie prefer the face over an edge or a corner (the diagonal of a flat face is an "edge" of both its triangles)
         if (d2 < bestD2 - 1e-18 || (d2 < bestD2 + 1e-18 && kind < bestKind))
         {
            bestD2 = d2;
            bestP = p;
            bestI = i;
            bestKind = kind;
         }
      }
      }
      if (bestI < 0)
         return pr;
      const double d = std::sqrt(bestD2);
      pr.point = bestP;
      pr.seg = bestI;
      pr.dist = d;
      pr.delta = ball.R + halfThick - d;
      pr.overlap = pr.delta > 0.0;
      const SurfaceTri &t = tris[bestI];
      V3 nf = t.n;
      if ((c - bestP).Dot(nf) < 0.0)
         nf = nf * -1.0;
      if (bestKind == 0 || d <= 1e-12)
         pr.normal = nf;
      else
         pr.normal = (c - bestP) / d;
      double edgeW = 0.0;   // 0 flat .. 1 fully rounded
      if (bestKind != 0 && d > 1e-12)
      {
         const double ang = std::acos(std::min(1.0, std::max(-1.0, pr.normal.Dot(nf))));
         edgeW = std::min(1.0, std::max(0.0, (ang - 5.0 * M_PI / 180.0) / (10.0 * M_PI / 180.0)));
      }
      const double R1 = 1.0 / ((1.0 - edgeW) / ball.R + edgeW * (1.0 / ball.R + 1.0 / edgeRho));
      pr.reff = std::sqrt(R1 * ball.R);
      return pr;
   }

   // ProbeTriangles' result with the other faces pressing the ball: the nearest point inside a face of another plane (a concave corner of the mesh). A face in a plane already counted adds nothing.
   inline void PatchTriangles(const BallSI &ball, const std::vector<SurfaceTri> &tris, const std::vector<int> &cand, const V3 &c, const double zero, TubeProbe &pr, const FaceRead *fr = nullptr)
   {
      if (!pr.overlap || pr.seg < 0 || pr.reff <= 0.0 || pr.delta <= zero)
         return;
      const double a0 = std::max(std::sqrt(pr.reff * (pr.delta - zero)), 1e-6);
      const auto fitOn = [&](const SurfaceTri &t, const V3 &p) {   // how much of the patch is on this face
         const auto edge = [&](const V3 &u, const V3 &v) { const V3 e = v - u; const double l = e.Length(); return l > 0.0 ? (p - u).Cross(e).Length() / l : 0.0; };
         return std::min(1.0, std::min({ edge(t.a, t.b), edge(t.b, t.c), edge(t.c, t.a) }) / a0);
      };
      V3 seen[4];
      int nSeen = 0;
      FoldSum sum;
      bool any = false;
      double fitBest = -1.0;
      { const ContactProf::LoopTicks profLoop(cand.size());   // P3a profile (diagnostics only)
      for (size_t j = 0; j < cand.size(); j++)
      {
         const int i = cand[j];
         if (i == pr.seg)
            continue;
         const SurfaceTri &t = tris[i];
         const double h = (c - t.a).Dot(t.n);
         const double ah = std::abs(h);
         if (ball.R - ah <= zero || ah <= 1e-12)
            continue;
         const V3 ni = h < 0.0 ? t.n * -1.0 : t.n;
         bool same = (ni - pr.normal).Dot(ni - pr.normal) < 1e-12;
         for (int k = 0; !same && k < nSeen; k++)
            same = (ni - seen[k]).Dot(ni - seen[k]) < 1e-12;
         if (same)
            continue;
         int kind = 0;
         V3 p;
         if (fr != nullptr)   // shared_finder P3a: the table's measurement
         {
            const FaceTable::Slot &sl = fr->At(j, t, c);
            p = sl.p;
            kind = sl.kind;
         }
         else
            p = ClosestOnTriangle(t, c, kind);
         if (kind != 0)
            continue;
         if (fitBest < 0.0)
            fitBest = fitOn(tris[pr.seg], pr.point);
         const double fit = std::min(fitBest, fitOn(t, p));
         if (fit <= 0.0)
            continue;
         if (nSeen < 4)
            seen[nSeen++] = ni;
         if (!any)
         {
            sum.Add(pr.delta - zero, 1.0, pr.normal);
            any = true;
         }
         sum.Add(ball.R - ah - zero, fit, ni);
      }
      }
      if (any)
         sum.Finish(pr, zero);
   }

   // ---- multi_contact 1: every part of one surface that presses the ball is its own Hertz bar (to-do 1.33, piece 2) ----
   // The depth over the surface, d(x) = R + rho - dist(c, x) - zero, has a peak on each part pressing the ball. Two peaks joined over the surface by a path that stays pressed are one
   // contact up to the lowest point of the best such path (the hump, a shared edge or a joint): the stored energy is the level form E = sum over peaks of V(d) - V(pass), pass the
   // hump where the peak joins a deeper one (0 if it never does), V the Hertz potential ~ d^2.5. For Hertz it is one depth with the same energy, the generalised fold sum
   // s25 = sum (d^2.5 - pass^2.5), g = sum (d^1.5 n - pass^1.5 n'), delta = s25^0.4 + zero, n = g / |g|, grad = |g| / s25^0.6 (FoldSum::Finish), n' the hump's normal.
   // The peaks and humps come from a union-find over the pieces taken deepest first (a merge tree): a piece with no pressed neighbour starts a peak; a piece that touches two
   // groups joins them at its own depth, the shallower group's peak passing there. A piece that joins a deeper group passes at its own depth: d = pass, it adds exactly nothing.
   struct BarNode
   {
      double d;      // depth above the found-inside zero, m
      V3 n;          // the depth's gradient direction at the piece (surface -> ball centre)
      int order;     // on equal depths: pieces (0) before the humps that join them (1)
      double r = 0.0;   // peak_slope: the distance from the ball's centre to the piece's nearest point, about which n turns as the ball moves (0: a face's inside, n fixed)
      V3 t;             // patch_friction (C4): a line piece's unit axis (tube segment, mesh edge inside its ends): n turns by (I - n n^T - t t^T) / r; 0 a joint or vertex
   };

   // peak_slope (runs/2026-11-22-per-peak plan 2.5a): every node that enters the level form's sum (s = +1 a peak, -1 its pass), for the elastic push's Newton slope taken per
   // peak. MergeBars writes it when it finds 2 or more peaks and empties it otherwise; thread-local, capacity reused (no allocation per step once grown).
   struct PeakTerm
   {
      double s, d, r;   // sign, depth above the zero (m), turning distance (m, 0 none)
      V3 n;
      V3 t;             // a line piece's unit axis (BarNode::t), 0 none
   };
   inline std::vector<PeakTerm> &PeakTerms()
   {
      thread_local std::vector<PeakTerm> v;
      return v;
   }

   // patch_friction (plan 2.1, piece 2b-1): each patch's pressed pieces (permanent ids: a mesh face t, a shared edge -(1 + 3 t + e) of its lower face t; a tube segment i, a joint
   // -(j + 1)), patch p's ids id[off[p]] .. id[off[p + 1] - 1]. Written by BarsTriangles / BarsTube beside their patch list; thread-local, capacity reused.
   // Each patch also names its deepest piece's signed geometry (plan 2.5 C1: a patch not pressed at the step's middle is evaluated there from it): kind 0 a face's plane (its
   // depth d - n . (c' - c)), 1 a segment a-b (a shared edge, a tube segment), 2 a point a (a joint, a face's rim point): depth d + dist(c) - dist(c')
   struct PieceGeo
   {
      int kind = 0;
      V3 a, b;
      V3 c;   // level_stiffness (fault B2): kind 3, a mesh face's own triangle a-b-c (PieceData only; its depth elsewhere from its closest point, ClosestOnTriangle)
   };
   // level_stiffness (B-lin, plan 2.4b): one pressed piece beside its id: its depth (raw, as a probe's delta: the zero added back), its normal (the depth's gradient direction), its own
   // signed geometry (PieceGeo, as above for the piece itself; a mesh face: kind 3, its triangle) and its two neighbours in the merge tree (static from the ids: a mesh edge -(1 + 3 t + e) joins faces t and
   // nb[t][e]; a tube segment i joins joints -(i + 1), -(i + 2); a joint -(j + 1) joins segments j - 1 and j (a ring's joint 0: n - 2 and 0); kNoPiece none, a face has none)
   struct PieceData
   {
      double d = -1.0;
      V3 n;
      PieceGeo g;
      int nbA = 0x7fffffff, nbB = 0x7fffffff;
      int nbC = 0x7fffffff, nbD = 0x7fffffff;   // a tube joint: the joints before and after it (through segments nbA, nbB); MergeBars links two pressed ones when that segment has no pressed foot
   };
   struct PatchPieceList
   {
      std::vector<int> off, id;
      std::vector<PieceGeo> geo;
      std::vector<PieceData> pd;   // level_stiffness: per id (parallel to id)
      void Clear() { off.assign(1, 0); id.clear(); geo.clear(); pd.clear(); }
      void One(const int piece, const double d, const V3 &n)   // the probe's nearest piece alone (d raw, n its normal)
      {
         Clear();
         if (piece != kNoPiece)
         {
            id.push_back(piece);
            PieceData p;
            p.d = d;
            p.n = n;
            pd.push_back(p);
         }
         off.push_back((int)id.size());
         geo.push_back(PieceGeo());
      }
      static constexpr int kNoPiece = 0x7fffffff;
   };
   // patch_friction (plan 2.1): one piece's memory, its patch's shear stretch and rolling and twist angles at the last accepted step's end
   struct PieceMem
   {
      int id;
      V3 xi;
      RollAngle ra;
      double kt = 0.0, kr = 0.0, kw = 0.0;   // plan 2.4's memory rule: its group's end stiffness the stretch and angles were stored at (shear, rolling, twist), 0 none
      int sg = -1;                           // that group's index at that step (the pieces of one stored group share it)
   };
   // patch_friction (plan 3.6): the ridge arrival's input, the item's step-start probe (the group step's Hit, valid while its positions are): its patches' pieces (off, id) and its
   // level form's passes and peaks (PeakTerm); a copy the caller fills (thread-local, capacity reused)
   struct RidgeStart
   {
      PatchPieceList pieces;
      std::vector<PeakTerm> passes;
   };
   // step pieces (runs/2026-11-22-ringcap design 2.4, 2.7): one entry of a static surface's step-start lookup for the shared step rule: another convex piece in reach, or a ridge
   struct StepPiece
   {
      double d;      // its depth (> 0 pressed, < 0 the gap): a piece's R (+ a tube's radius) - distance before the item's zero; a ridge's already from the zero (zeroed)
      V3 n;          // its normal (surface -> ball centre): the closing speed is -v . n
      bool zeroed;
   };
   // link_windows (rule T, to-do 1.33; runs/2026-11-23-link-rule brief_code.md, rule_t.py): one link of a static surface pressed during the step's straight move c(t) = c0 + t (c1 - c0):
   // the pieces it joins (ids as the probe's: face t, tube segment i, tube joint -(j + 1)) and the part of the move [t0, t1] (0 <= t0 <= t1 <= 1) where it is within reach of the centre
   struct LinkWindow
   {
      int a, b;
      double t0, t1;
   };
   // VPX_PATCH_CHECK=1 (diagnostics, read once): cut B's pressed faces from the probe's hits are also compared with each face probed alone (BarsTriangles): [0] shortcut calls
   // compared, [1] calls whose pressed set differs; logged on the patch check's totals line
   inline bool PatchCheckEnv()
   {
      static const bool on = [] { const char *e = std::getenv("VPX_PATCH_CHECK"); return e != nullptr && e[0] == '1'; }();
      return on;
   }
   inline std::atomic<long long> *CutBCheck()   // [2] calls whose listed pieces (faces and pressed shared edges) give more than one peak through MergeBars (expected 0)
   {
      static std::atomic<long long> c[3] { { 0 }, { 0 }, { 0 } };
      return c;
   }
   // VPX_PATCH_CHECK=1 (diagnostics): the friction groups' held spring energy at a step's start (as Eval's friction block takes it) against the item's book energy there
   // (Group::Energy's ShearE + rolling terms): [0] attempts compared, [1] |difference| > 1e-12 J, [2] book per-step rises > 1e-9 J (VPX_GROUP_BOOK_STEPS); on the patch check's line
   inline std::atomic<long long> *HeldCheck()
   {
      static std::atomic<long long> c[6] { { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 } };   // [4] B-lin lev1 checked, [5] lev1 vs the end probe's level form, rel > 1e-12   // [3] kind C: a group's start stiffness k0 differs from its stored ks (relative > 1e-9)
      return c;
   }
   inline PatchPieceList &PatchPieces()
   {
      thread_local PatchPieceList v;
      return v;
   }
   // MergeBars' node -> its patch (entry in patches; -1 none), written when patches is asked for
   inline std::vector<int> &NodePatch()
   {
      thread_local std::vector<int> v;
      return v;
   }
   // level_stiffness (B2): a mesh face's own geometry, its triangle (kind 3)
   inline PieceGeo FaceGeo(const SurfaceTri &t)
   {
      PieceGeo g;
      g.kind = 3;
      g.a = t.a;
      g.b = t.b;
      g.c = t.c;
      return g;
   }
   // the pieces of each patch from MergeBars' NodePatch (np patches) and each node's piece id
   inline void FillPatchPieces(const std::vector<int> &nodePiece, const std::vector<PieceData> &nodePd, const int np)
   {
      PatchPieceList &pl = PatchPieces();
      const std::vector<int> &npch = NodePatch();
      pl.off.assign(np + 1, 0);
      pl.id.clear();
      pl.geo.clear();
      pl.pd.clear();
      for (int p = 0; p < np; p++)
      {
         for (size_t x = 0; x < nodePiece.size(); x++)
            if (npch[x] == p)
            {
               pl.id.push_back(nodePiece[x]);
               pl.pd.push_back(nodePd[x]);
            }
         pl.off[p + 1] = (int)pl.id.size();
      }
   }
   // one patch holding every node (one peak, or several peaks in one patch)
   inline void AllPatchPieces(const std::vector<int> &nodePiece, const std::vector<PieceData> &nodePd)
   {
      PatchPieceList &pl = PatchPieces();
      pl.Clear();
      pl.id.assign(nodePiece.begin(), nodePiece.end());
      pl.pd.assign(nodePd.begin(), nodePd.end());
      pl.off.push_back((int)pl.id.size());
   }
   // the deepest of the nodes (MergeBars' order: depth, then index)
   inline int DeepestNode(const std::vector<BarNode> &nodes)
   {
      int best = -1;
      for (int x = 0; x < (int)nodes.size(); x++)
         if (best < 0 || nodes[x].d > nodes[best].d)
            best = x;
      return best;
   }

   // patch_contacts (runs/2026-11-22-per-peak plan 3.1): one pressed patch (a merge tree's root: every piece joined to it by a path that stays pressed), its own fold sum over its own
   // peaks (pass terms included) and its deepest node (the root)
   struct PatchSum
   {
      FoldSum sum;
      int deep = -1;
   };

   // nodes: the pieces pressing the ball (d > 0); adj: pairs of nodes that touch on the surface. Returns the number of peaks; sum holds their s25 and g.
   // patches (optional): one entry per patch, in the order of its first peak; sum itself is unchanged (the same additions in the same order).
   inline int MergeBars(const std::vector<BarNode> &nodes, const std::vector<std::pair<int, int>> &adj, FoldSum &sum, std::vector<PatchSum> *patches = nullptr)
   {
      const int m = (int)nodes.size();
      thread_local std::vector<int> ord, parent, peak, first, nbr, fill;   // scratch kept per thread: capacity reused, no heap use per call
      ord.resize(m);
      parent.resize(m);
      peak.resize(m);
      first.assign(m + 1, 0);
      nbr.resize(adj.size() * 2);
      for (int i = 0; i < m; i++)
         ord[i] = parent[i] = peak[i] = i;
      std::sort(ord.begin(), ord.end(), [&](const int a, const int b) {
         if (nodes[a].d != nodes[b].d)
            return nodes[a].d > nodes[b].d;
         return nodes[a].order != nodes[b].order ? nodes[a].order < nodes[b].order : a < b;
      });
      for (const auto &e : adj)
      {
         first[e.first + 1]++;
         first[e.second + 1]++;
      }
      for (int i = 0; i < m; i++)
         first[i + 1] += first[i];
      {
         fill.assign(first.begin(), first.end() - 1);
         for (const auto &e : adj)
         {
            nbr[fill[e.first]++] = e.second;
            nbr[fill[e.second]++] = e.first;
         }
      }
      thread_local std::vector<char> on, passed;
      on.assign(m, 0);
      passed.assign(m, 0);
      const auto find = [&](int x) { while (parent[x] != x) x = parent[x] = parent[parent[x]]; return x; };
      thread_local std::vector<double> pass;
      thread_local std::vector<V3> passN;
      pass.assign(m, 0.0);
      passN.assign(m, V3());
      thread_local std::vector<double> passR;   // peak_slope: the pass node's turning distance
      passR.assign(m, 0.0);
      thread_local std::vector<V3> passT;   // its line axis (C4)
      passT.assign(m, V3());
      std::vector<PeakTerm> &terms = PeakTerms();
      terms.clear();
      const auto deeper = [&](const int pa, const int pb) { return nodes[pa].d != nodes[pb].d ? nodes[pa].d > nodes[pb].d : pa < pb; };
      for (const int x : ord)
      {
         on[x] = 1;
         for (int k = first[x]; k < first[x + 1]; k++)
         {
            const int y = nbr[k];
            if (!on[y])
               continue;
            const int rx = find(x), ry = find(y);
            if (rx == ry)
               continue;
            int keep = rx, lose = ry;   // the group whose peak is shallower passes at x
            if (deeper(peak[ry], peak[rx]))
               std::swap(keep, lose);
            const int pl = peak[lose];
            pass[pl] = nodes[x].d;
            passN[pl] = nodes[x].n;
            passR[pl] = nodes[x].r;
            passT[pl] = nodes[x].t;
            passed[pl] = 1;
            parent[lose] = keep;
         }
      }
      thread_local std::vector<int> slot;   // patch_contacts: root -> its entry in patches, -1 none yet
      if (patches != nullptr)
      {
         patches->clear();
         slot.assign(m, -1);
      }
      int peaks = 0;
      for (int i = 0; i < m; i++)
      {
         if (peak[i] != i || (passed[i] && pass[i] == nodes[i].d))
            continue;   // not a peak (every node starts as its own peak; one that joins at its own depth adds nothing)
         if (!passed[i] && find(i) != i)
            continue;
         const double d = nodes[i].d, d15 = d * std::sqrt(d);
         sum.s25 += d * d15;
         sum.g = sum.g + nodes[i].n * d15;
         terms.push_back({ 1.0, d, nodes[i].r, nodes[i].n, nodes[i].t });
         if (passed[i] && pass[i] > 0.0)
         {
            const double p = pass[i], p15 = p * std::sqrt(p);
            sum.s25 -= p * p15;
            sum.g = sum.g - passN[i] * p15;
            terms.push_back({ -1.0, p, passR[i], passN[i], passT[i] });
         }
         if (patches != nullptr)   // the same terms into the peak's own patch (a root is its group's deepest node: the deeper group always keeps its root)
         {
            const int r = find(i);
            if (slot[r] < 0)
            {
               slot[r] = (int)patches->size();
               patches->push_back(PatchSum());
               patches->back().deep = r;
            }
            FoldSum &ps = (*patches)[slot[r]].sum;
            ps.s25 += d * d15;
            ps.g = ps.g + nodes[i].n * d15;
            if (passed[i] && pass[i] > 0.0)
            {
               const double p = pass[i], p15 = p * std::sqrt(p);
               ps.s25 -= p * p15;
               ps.g = ps.g - passN[i] * p15;
            }
         }
         peaks++;
      }
      if (patches != nullptr)   // patch_friction: each node's patch
      {
         std::vector<int> &np = NodePatch();
         np.resize(m);
         for (int x = 0; x < m; x++)
            np[x] = slot[find(x)];
      }
      if (peaks <= 1)
         terms.clear();   // one peak: today's slope, no list
      return peaks;
   }

   // patch_zero (to-do 1.33 piece 2b-3; runs/2026-11-23-zero-patch brief_code.md Z1, Z2, D6): the live found-inside zeros of the static item being probed, (piece id, its patch's
   // zero) (GroupStep::Item::tZero), set by the group step around its probe of that item; null: the probe takes the scalar zero it is given. Set (even empty): the caller passes
   // zero 0, the patches are found on raw depth and each is shifted by its own zero (a patch with no entry: 0)
   inline const std::vector<std::pair<int, double>> *&PatchZeros()
   {
      thread_local const std::vector<std::pair<int, double>> *p = nullptr;
      return p;
   }
   // patch_zero (rule 5, brief_code.md Z5): true: each node takes its own piece's entry (none: 0) instead of its patch's largest: a joined patch's parts with the zeros they carried
   inline bool &PatchZerosPerPiece()
   {
      thread_local bool on = false;
      return on;
   }
   // patch_zero (F2): set by ShiftedBars when a raw patch's pieces carry different entries (its largest taken: a net); the group step clears and reads it around its probe
   inline bool &PatchZeroNet()
   {
      thread_local bool on = false;
      return on;
   }
   inline double PieceZero(const std::vector<std::pair<int, double>> &zl, const int id)
   {
      for (const std::pair<int, double> &z : zl)
         if (z.first == id)
            return z.second;
      return 0.0;
   }
   // patch_zero (Z2): one raw patch's shifted answer: its fold sum over its shifted nodes, its deepest node (raw and shifted alike: one zero per patch) and that node's shifted depth
   struct ZeroPatch
   {
      FoldSum sum;
      int deep = -1;
      double top = 0.0;
      double z = 0.0;
   };
   // patch_zero (Z2): nodes, adj RAW (the probe's zero 0), piece each node's id. The raw patches (MergeBars), each one's zero = the largest entry of zl among its pieces, then the
   // level form of every patch on its own nodes shifted by its zero (d - z_p; a node at or below 0 drops out), all in one MergeBars: no link joins two raw patches, so the sums add.
   // Returns the shifted peaks; sum the shifted fold sum; zp one entry per raw patch (MergeBars' order); nodeZ each node's patch zero. NodePatch() is left the RAW one; PeakTerms()
   // as MergeBars on the shifted nodes leaves it (2 or more shifted peaks: the list)
   inline int ShiftedBars(const std::vector<BarNode> &nodes, const std::vector<std::pair<int, int>> &adj, const std::vector<int> &piece, const std::vector<std::pair<int, double>> &zl,
                          FoldSum &sum, std::vector<ZeroPatch> &zp, std::vector<double> &nodeZ)
   {
      thread_local std::vector<PatchSum> rawP, shP;
      thread_local std::vector<int> npRaw, shOf, orig;
      thread_local std::vector<BarNode> sn;
      thread_local std::vector<std::pair<int, int>> sa;
      FoldSum rawSum;
      MergeBars(nodes, adj, rawSum, &rawP);
      npRaw = NodePatch();
      const int m = (int)nodes.size(), np = (int)rawP.size();
      zp.assign(np, ZeroPatch());
      for (int x = 0; x < m; x++)
         if (npRaw[x] >= 0)
            zp[npRaw[x]].z = std::max(zp[npRaw[x]].z, PieceZero(zl, piece[x]));
      if (!PatchZerosPerPiece())   // F2: an entry below its patch's largest (pieces with no entry are not counted)
         for (int x = 0; x < m; x++)
            if (npRaw[x] >= 0)
            {
               const double e = PieceZero(zl, piece[x]);
               if (e > 0.0 && e < zp[npRaw[x]].z)
                  PatchZeroNet() = true;
            }
      nodeZ.assign(m, 0.0);
      sn.clear();
      orig.clear();
      shOf.assign(m, -1);
      for (int x = 0; x < m; x++)
      {
         const double z = PatchZerosPerPiece() ? PieceZero(zl, piece[x]) : npRaw[x] >= 0 ? zp[npRaw[x]].z : 0.0;
         nodeZ[x] = z;
         if (nodes[x].d - z > 0.0)
         {
            shOf[x] = (int)sn.size();
            BarNode b = nodes[x];
            b.d -= z;
            sn.push_back(b);
            orig.push_back(x);
         }
      }
      for (int p = 0; p < np; p++)
      {
         zp[p].deep = rawP[p].deep;
         zp[p].top = nodes[rawP[p].deep].d - zp[p].z;
      }
      sa.clear();
      for (const std::pair<int, int> &e : adj)
         if (shOf[e.first] >= 0 && shOf[e.second] >= 0)
            sa.emplace_back(shOf[e.first], shOf[e.second]);
      int peaks = 0;
      if (sn.empty())
         PeakTerms().clear();
      else
      {
         peaks = MergeBars(sn, sa, sum, &shP);
         for (const PatchSum &q : shP)   // each shifted patch into the raw patch it lies in
         {
            const int p = npRaw[orig[q.deep]];
            zp[p].sum.s25 += q.sum.s25;
            zp[p].sum.g = zp[p].sum.g + q.sum.g;
         }
      }
      NodePatch() = npRaw;
      return peaks;
   }

   // The nearest point of segment p0 -> p1 (the ends exactly at a clamp, so two pieces that meet at a corner report the same bits there)
   inline V3 SegNearest(const V3 &p0, const V3 &p1, const V3 &c)
   {
      const V3 t = p1 - p0;
      const double L2 = t.Dot(t);
      if (L2 <= 0.0)
         return p0;
      const double s = (c - p0).Dot(t) / L2;
      return s <= 0.0 ? p0 : s >= 1.0 ? p1 : p0 + t * s;
   }

   // link_windows (rule T; runs/2026-11-23-link-rule rule_t.py win_seg, win_pt): the window [t0, t1] of t in [0, 1] where the centre c0 + t m is within r of the segment p -> q (of the
   // point p), closed form: the squared distance is a quadratic in t on each region (the foot inside, the end p, the end q), the distance is convex in t so the window is one interval.
   // False: none. A crossing the roots miss at a region's border (rounding) is found by bisection on the distance (both ends known on opposite sides of r)
   inline int LinkQuadRoots(const double a, const double b, const double c, const double aMin, double *r)
   {
      if (a <= aMin)   // (the move along the region's level lines: no crossing inside it; its ends are tested apart)
         return 0;
      const double D = b * b - 4.0 * a * c;
      if (D < 0.0)
         return 0;
      const double s = std::sqrt(D);
      r[0] = (-b - s) / (2.0 * a);
      r[1] = (-b + s) / (2.0 * a);
      return 2;
   }
   template <class Dist> inline bool LinkWinFinish(const Dist &dist, const double r, const double *cand, const int nc, double &t0, double &t1)
   {
      const bool in0 = dist(0.0) < r, in1 = dist(1.0) < r;
      double lo = 2.0, hi = -1.0;
      for (int i = 0; i < nc; i++)
      {
         lo = std::min(lo, cand[i]);
         hi = std::max(hi, cand[i]);
      }
      if (nc == 0 && !in0 && !in1)
         return false;
      const auto bisect = [&](double a, double b) {   // dist(a) < r <= dist(b) or the reverse: the crossing
         const bool ina = dist(a) < r;
         for (int it = 0; it < 60; it++)
         {
            const double m = 0.5 * (a + b);
            ((dist(m) < r) == ina ? a : b) = m;
         }
         return 0.5 * (a + b);
      };
      t0 = in0 ? 0.0 : nc > 0 ? lo : bisect(1.0, 0.0);
      t1 = in1 ? 1.0 : nc > 0 ? hi : bisect(0.0, 1.0);
      return t0 <= t1;
   }
   inline bool LinkWinPt(const V3 &c0, const V3 &m, const V3 &p, const double r, double &t0, double &t1)
   {
      const V3 w = c0 - p;
      const double mm = m.Dot(m);
      double rt[2], cand[2];
      int nc = 0;
      const int n = LinkQuadRoots(mm, 2.0 * w.Dot(m), w.Dot(w) - r * r, 0.0, rt);
      for (int i = 0; i < n; i++)
         if (rt[i] >= 0.0 && rt[i] <= 1.0)
            cand[nc++] = rt[i];
      return LinkWinFinish([&](const double t) { return (w + m * t).Length(); }, r, cand, nc, t0, t1);
   }
   inline bool LinkWinSeg(const V3 &c0, const V3 &m, const V3 &p, const V3 &q, const double r, double &t0, double &t1)
   {
      const V3 u = q - p;
      const double L2 = u.Dot(u);
      if (L2 <= 0.0)
         return LinkWinPt(c0, m, p, r, t0, t1);
      const V3 w0 = c0 - p, w1 = c0 - q;
      const double mu = m.Dot(u), w0u = w0.Dot(u), mm = m.Dot(m), aMin = 1e-12 * mm;
      const auto sOf = [&](const double t) { return (w0u + t * mu) / L2; };
      double rt[2], cand[6];
      int nc = 0;
      int n = LinkQuadRoots(mm - mu * mu / L2, 2.0 * (w0.Dot(m) - w0u * mu / L2), w0.Dot(w0) - w0u * w0u / L2 - r * r, aMin, rt);   // the foot inside
      for (int i = 0; i < n; i++)
         if (rt[i] >= 0.0 && rt[i] <= 1.0 && sOf(rt[i]) >= 0.0 && sOf(rt[i]) <= 1.0)
            cand[nc++] = rt[i];
      n = LinkQuadRoots(mm, 2.0 * w0.Dot(m), w0.Dot(w0) - r * r, 0.0, rt);   // the end p
      for (int i = 0; i < n; i++)
         if (rt[i] >= 0.0 && rt[i] <= 1.0 && sOf(rt[i]) <= 0.0)
            cand[nc++] = rt[i];
      n = LinkQuadRoots(mm, 2.0 * w1.Dot(m), w1.Dot(w1) - r * r, 0.0, rt);   // the end q
      for (int i = 0; i < n; i++)
         if (rt[i] >= 0.0 && rt[i] <= 1.0 && sOf(rt[i]) >= 1.0)
            cand[nc++] = rt[i];
      return LinkWinFinish([&](const double t) { const V3 c = c0 + m * t; return (c - SegNearest(p, q, c)).Length(); }, r, cand, nc, t0, t1);
   }

   // A tube (centreline pts, n points; a closed ring repeats its first point at the end, the two are one joint): the pieces are the segments' inside nearest points and the joints.
   // One peak: ProbeTube's result unchanged. More: the level form above, reff and the nearest point kept from the deepest (held reff, piece 2).
   // segs, pr0 (shared_finder P5): the segments to read, in index order, and their probe (ProbeTube's result, given by the caller); null: all, probed here
   inline TubeProbe BarsTube(const BallSI &ball, const V3 *pts, const int n, const double rho, const V3 &c, const TubePre &pre, const double zero, const std::vector<int> *segs = nullptr,
                             const TubeProbe *pr0 = nullptr, std::vector<TubeProbe> *patchesOut = nullptr)
   {
      TubeProbe pr = pr0 != nullptr ? *pr0 : ProbeTube(ball, pts, n, rho, c, &pre);
      bool pcs = false;   // patch_friction: the patch pieces written (else one() takes the probe's nearest piece)
      const auto one = [&]() {   // patch_contacts: one patch (or none): today's answer is its single entry
         if (patchesOut != nullptr)
         {
            patchesOut->clear();
            patchesOut->push_back(pr);
            if (!pcs)
               PatchPieces().One(pr.overlap && pr.seg >= 0 ? pr.seg : PatchPieceList::kNoPiece, pr.delta, pr.normal);
         }
         return pr;
      };
      const std::vector<std::pair<int, double>> *const zl = PatchZeros();   // patch_zero: set, the caller passed zero 0 (Z2)
      if (!pr.overlap || pr.seg < 0 || pr.reff <= 0.0 || pr.delta <= zero || n < 3)
      {
         if (zl != nullptr && pr.overlap && pr.seg >= 0)   // patch_zero: the one piece, shifted by its patch's zero
            pr.delta -= PieceZero(*zl, pr.seg);
         return one();
      }
      const double S = ball.R + rho;
      const bool ring = pts[0].x == pts[n - 1].x && pts[0].y == pts[n - 1].y && pts[0].z == pts[n - 1].z;
      const int nv = ring ? n - 1 : n;   // joints: point n-1 of a ring is joint 0
      thread_local std::vector<BarNode> nodes;   // scratch kept per thread
      thread_local std::vector<std::pair<int, int>> adj;
      thread_local std::vector<int> jointNode;
      thread_local std::vector<int> nodeSeg;   // patch_contacts: node -> its piece (segment i, or joint j as -(j + 1))
      nodes.clear();
      adj.clear();
      nodeSeg.clear();
      jointNode.assign(nv, -2);   // -2 not looked at, -1 not pressed
      const auto joint = [&](int j) {
         if (j >= nv)
            j = 0;
         if (jointNode[j] == -2)
         {
            const V3 dv = c - pts[j];
            const double dist = dv.Length();
            const double d = S - dist - zero;
            jointNode[j] = -1;
            if (d > 0.0 && dist > 1e-12)
            {
               jointNode[j] = (int)nodes.size();
               nodes.push_back({ d, dv / dist, 1, dist });
               nodeSeg.push_back(-(j + 1));
            }
         }
         return jointNode[j];
      };
      const int cnt = segs != nullptr ? (int)segs->size() : n - 1;
      for (int k = 0; k < cnt; k++)
      {
         const int i = segs != nullptr ? (*segs)[k] : k;
         const V3 dm = c - pre.mid[i];
         const double lim = S - zero + pre.rad[i];
         if (dm.Dot(dm) > lim * lim)
            continue;
         const int j0 = joint(i), j1 = joint(i + 1);
         const double s = pre.L2[i] > 0.0 ? (c - pts[i]).Dot(pre.t[i]) / pre.L2[i] : 0.0;
         if (s > 0.0 && s < 1.0)
         {
            const V3 dv = c - (pts[i] + pre.t[i] * s);
            const double dist = dv.Length();
            const double d = S - dist - zero;
            if (d > 0.0 && dist > 1e-12)
            {
               const int k = (int)nodes.size();
               nodes.push_back({ d, dv / dist, 0, dist, pre.t[i] * (1.0 / std::sqrt(pre.L2[i])) });
               nodeSeg.push_back(i);
               if (j0 >= 0)
                  adj.emplace_back(k, j0);
               if (j1 >= 0)
                  adj.emplace_back(k, j1);
               continue;
            }
         }
         if (j0 >= 0 && j1 >= 0)   // the segment's depth runs from one joint to the other without a peak inside: the joints touch along it
            adj.emplace_back(j0, j1);
      }
      FoldSum sum;
      thread_local std::vector<PatchSum> patches;
      const auto tubeGeo = [&](const int x) {   // patch_friction: a node's piece geometry (PieceGeo)
         PieceGeo g;
         const int s = x >= 0 ? nodeSeg[x] : 0;
         if (x < 0)
            return g;
         if (s >= 0)
         {
            g.kind = 1;
            g.a = pts[s];
            g.b = pts[s + 1];
         }
         else
         {
            g.kind = 2;
            g.a = pts[-s - 1];
         }
         return g;
      };
      thread_local std::vector<PieceData> nodePd;   // level_stiffness: each node's piece data (PieceData), filled where the pieces are written
      const auto tubePd = [&]() {
         nodePd.resize(nodes.size());
         for (size_t x = 0; x < nodes.size(); x++)
         {
            PieceData &p = nodePd[x];
            p.d = nodes[x].d + zero;
            p.n = nodes[x].n;
            p.g = tubeGeo((int)x);
            const int s = nodeSeg[x];
            if (s >= 0)   // a segment: its two joints (a ring's point n - 1 is joint 0)
            {
               p.nbA = -(s + 1);
               p.nbB = -((s + 1 >= nv ? 0 : s + 1) + 1);
            }
            else   // a joint: the segments that meet there
            {
               const int j = -s - 1;
               p.nbA = j > 0 ? j - 1 : ring ? n - 2 : PatchPieceList::kNoPiece;
               p.nbB = j < n - 1 ? j : PatchPieceList::kNoPiece;
               p.nbC = j > 0 ? -j : ring ? -nv : PatchPieceList::kNoPiece;                          // (the joint before: j - 1; a ring's joint 0: nv - 1)
               p.nbD = j < n - 1 ? -((j + 1 >= nv ? 0 : j + 1) + 1) : PatchPieceList::kNoPiece;   // (the joint after: j + 1; a ring's last: joint 0)
            }
         }
      };
      if (zl != nullptr)   // patch_zero (Z2): the raw patches, each shifted by its own zero; the merged answer and the patch list from the shifted level form
      {
         if (nodes.empty())
         {
            pr.delta -= PieceZero(*zl, pr.seg);
            return one();
         }
         thread_local std::vector<ZeroPatch> zp;
         thread_local std::vector<double> nodeZ;
         const int peaks = ShiftedBars(nodes, adj, nodeSeg, *zl, sum, zp, nodeZ);
         const int rawTop = DeepestNode(nodes);
         const TubeProbe base = pr;   // the raw deepest piece's probe
         const auto alone = [&](const int x) {   // a node's piece probed alone: a segment, or a joint as the segments that meet there (index order)
            if (x == rawTop)
               return base;
            thread_local std::vector<int> one2;
            const int s = nodeSeg[x];
            one2.clear();
            if (s >= 0)
               one2.push_back(s);
            else
            {
               const int j = -s - 1;
               if (j > 0)
                  one2.push_back(j - 1);
               if (j < n - 1)
                  one2.push_back(j);
               if (ring && j == 0)
                  one2.push_back(n - 2);
            }
            const int hint = pre.hint;
            const TubeProbe q = ProbeTube(ball, pts, n, rho, c, &pre, &one2);
            pre.hint = hint;
            return q;
         };
         int best = 0;   // the patch whose deepest node is deepest after the shift: the merged answer keeps its reff, point and seg
         for (int p = 1; p < (int)zp.size(); p++)
            if (zp[p].top > zp[best].top)
               best = p;
         pr = alone(zp[best].deep);
         if (peaks > 0)
            sum.Finish(pr, 0.0);
         else   // no node above its zero: the gap is the deepest shifted node's
         {
            pr.delta = zp[best].top;
            pr.normal = nodes[zp[best].deep].n;
            pr.grad = 1.0;
         }
         if (patchesOut != nullptr)
         {
            patchesOut->clear();
            for (const ZeroPatch &p : zp)
            {
               TubeProbe q = alone(p.deep);
               if (p.sum.s25 > 0.0)
                  p.sum.Finish(q, 0.0);
               else
               {
                  q.delta = p.top;
                  q.normal = nodes[p.deep].n;
                  q.grad = 1.0;
               }
               patchesOut->push_back(q);
            }
            tubePd();
            for (size_t x = 0; x < nodes.size(); x++)   // (B-lin: each piece's depth from its patch's zero)
               nodePd[x].d -= nodeZ[x];
            FillPatchPieces(nodeSeg, nodePd, (int)zp.size());
            for (const ZeroPatch &p : zp)
               PatchPieces().geo.push_back(tubeGeo(p.deep));
            pcs = true;
         }
         return pr;
      }
      if (MergeBars(nodes, adj, sum, patchesOut != nullptr ? &patches : nullptr) <= 1)
      {
         if (patchesOut != nullptr)
         {
            tubePd();
            AllPatchPieces(nodeSeg, nodePd);
            PatchPieces().geo.assign(1, tubeGeo(DeepestNode(nodes)));
            pcs = true;
         }
         return one();
      }
      const TubeProbe base = pr;   // the deepest piece's probe (reff, point, seg)
      sum.Finish(pr, zero);
      if (patchesOut != nullptr)
      {
         if (patches.size() <= 1)
         {
            tubePd();
            AllPatchPieces(nodeSeg, nodePd);   // (every node is in the one patch)
            PatchPieces().geo.assign(1, tubeGeo(DeepestNode(nodes)));
            pcs = true;
            return one();   // one patch: today's merged answer
         }
         patchesOut->clear();
         int top = patches[0].deep;   // the deepest node of all (MergeBars' tie rule): its patch keeps the probe's reff, point and seg
         for (const PatchSum &p : patches)
            if (nodes[p.deep].d > nodes[top].d || (nodes[p.deep].d == nodes[top].d && p.deep < top))
               top = p.deep;
         thread_local std::vector<int> one1;
         for (const PatchSum &p : patches)
         {
            TubeProbe q = base;
            if (p.deep != top)   // the patch's deepest piece probed alone: a segment, or a joint as the segments that meet there (index order)
            {
               const int s = nodeSeg[p.deep];
               one1.clear();
               if (s >= 0)
                  one1.push_back(s);
               else
               {
                  const int j = -s - 1;
                  if (j > 0)
                     one1.push_back(j - 1);
                  if (j < n - 1)
                     one1.push_back(j);
                  if (ring && j == 0)
                     one1.push_back(n - 2);
               }
               const int hint = pre.hint;
               q = ProbeTube(ball, pts, n, rho, c, &pre, &one1);
               pre.hint = hint;
            }
            p.sum.Finish(q, zero);
            patchesOut->push_back(q);
         }
         tubePd();
         FillPatchPieces(nodeSeg, nodePd, (int)patches.size());
         for (const PatchSum &p : patches)
            PatchPieces().geo.push_back(tubeGeo(p.deep));
      }
      return pr;
   }

   // A mesh's links for BarsTriangles, from its faces and their welded vertex ids (vid: vertices equal in all three float32 bits of the VPX-unit vertex are one): nb[t][e] the face
   // across edge e = (vertex e, vertex e+1), -1 none; fold[t][e] that edge's fold (0 flat, 1 / -1 the other face's third vertex above / below this face's plane, 2 degenerate).
   // Static meshes (RubberSurface::Link) and translating ones (SlideMesh, contactmovingpart.h) build them here. Returns the number of faces left unlinked on an edge of three or more
   inline int LinkFaces(const std::vector<SurfaceTri> &tris, const std::vector<std::array<int, 3>> &vid, std::vector<std::array<int, 3>> &nb, std::vector<std::array<int8_t, 3>> &fold)
   {
      nb.assign(tris.size(), { -1, -1, -1 });
      std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> edges;   // (lower id, higher id) -> every (face, edge) on it, in face order
      for (int t = 0; t < (int)tris.size(); t++)
         for (int e = 0; e < 3; e++)
         {
            const int a = vid[t][e], b = vid[t][(e + 1) % 3];
            if (a != b)
               edges[std::make_pair(std::min(a, b), std::max(a, b))].emplace_back(t, e);
         }
      // An edge on two faces joins them. On three or more (a flat surface running on past a wall's foot, user 2026-11-16): the two faces that continue each other most nearly
      // flat are joined, greedily (the in-face directions away from the edge most nearly opposite), so the surface a ball rolls over stays one; the others stay unlinked
      const auto away = [&tris](int t, int e) {   // unit direction in face t, square to its edge e, from the edge into the face
         const V3 *v[3] = { &tris[t].a, &tris[t].b, &tris[t].c };
         const V3 p = *v[e], d = *v[(e + 1) % 3] - p, w = *v[(e + 2) % 3] - p;
         const double dd = d.Dot(d);
         const V3 s = dd > 0.0 ? w - d * (w.Dot(d) / dd) : w;
         const double l = std::sqrt(s.Dot(s));
         return l > 0.0 ? s / l : s;
      };
      int over = 0;
      for (const auto &[key, list] : edges)
      {
         if (list.size() == 2)
         {
            if (list[0].first != list[1].first)
            {
               nb[list[0].first][list[0].second] = list[1].first;
               nb[list[1].first][list[1].second] = list[0].first;
            }
            else
               over++;   // one face twice
            continue;
         }
         if (list.size() < 2)
            continue;
         std::vector<bool> used(list.size(), false);
         for (;;)
         {
            int bi = -1, bj = -1;
            double best = 2.0;
            for (size_t i = 0; i < list.size(); i++)
               for (size_t j = i + 1; j < list.size(); j++)
                  if (!used[i] && !used[j] && list[i].first != list[j].first)
                     if (const double c = away(list[i].first, list[i].second).Dot(away(list[j].first, list[j].second)); c < best)
                     {
                        best = c;
                        bi = (int)i;
                        bj = (int)j;
                     }
            if (bi < 0)
               break;
            used[bi] = used[bj] = true;
            nb[list[bi].first][list[bi].second] = list[bj].first;
            nb[list[bj].first][list[bj].second] = list[bi].first;
         }
         for (size_t i = 0; i < list.size(); i++)
            over += !used[i];   // a face left unlinked on an edge of three or more
      }
      // cut B: each shared edge's fold, from both faces. Flat band (derived, not tuned): every stored coordinate is a float32, rounded by at most 2^-24 of its size, so each vertex is off by
      // at most dv = sqrt(3) 2^-24 M (M the largest coordinate size of the four vertices; the SI scale adds 2^-53, neglected); the height of a point above a plane through three such vertices
      // is then off by at most dv (1 + |l0| + |l1| + |l2|), l the point's barycentric coordinates in the face. A fold whose two heights are both inside their bounds is flat.
      fold.assign(tris.size(), { 0, 0, 0 });
      const auto height = [](const SurfaceTri &A, const V3 &p, double &band) {
         double M = 0.0;
         for (const V3 *q : { &A.a, &A.b, &A.c, &p })
            M = std::max({ M, std::abs(q->x), std::abs(q->y), std::abs(q->z) });
         const V3 v0 = A.b - A.a, v1 = A.c - A.a, v2 = p - A.a;
         const double d00 = v0.Dot(v0), d01 = v0.Dot(v1), d11 = v1.Dot(v1), d20 = v2.Dot(v0), d21 = v2.Dot(v1), den = d00 * d11 - d01 * d01;
         const double l1 = den != 0.0 ? (d11 * d20 - d01 * d21) / den : INFINITY, l2 = den != 0.0 ? (d00 * d21 - d01 * d20) / den : INFINITY, l0 = 1.0 - l1 - l2;
         band = std::sqrt(3.0) * std::ldexp(M, -24) * (1.0 + std::abs(l0) + std::abs(l1) + std::abs(l2));
         return (p - A.a).Dot(A.n);
      };
      for (int t = 0; t < (int)tris.size(); t++)
         for (int e = 0; e < 3; e++)
         {
            const int o = nb[t][e];
            if (o < t)   // each shared edge once, from its lower face
               continue;
            int f = -1;
            for (int k = 0; k < 3; k++)
               if (nb[o][k] == t)
                  f = k;
            if (f < 0)
               continue;
            const SurfaceTri &A = tris[t], &B = tris[o];
            if (!A.area || !B.area)
            {
               fold[t][e] = fold[o][f] = 2;
               continue;
            }
            const V3 *va[3] = { &A.a, &A.b, &A.c }, *vb[3] = { &B.a, &B.b, &B.c };
            double bandA, bandB;
            const double hB = height(A, *vb[(f + 2) % 3], bandB);   // B's third vertex above A's plane
            const double hA = height(B, *va[(e + 2) % 3], bandA);   // A's third vertex above B's plane
            if (std::abs(hB) <= bandB && std::abs(hA) <= bandA)
               continue;   // flat
            fold[t][e] = hB > 0.0 ? 1 : hB < 0.0 ? -1 : 2;
            fold[o][f] = hA > 0.0 ? 1 : hA < 0.0 ? -1 : 2;
         }
      return over;
   }

   // link_windows (rule T) on a triangle mesh, one copy (the static mesh, ContactSurfaceHost::LinkWindows; a translating mesh, SlideMesh::LinkWindows on the relative segment): each
   // shared edge between two faces of the list L (nb's links), from its lower face, both faces with area, with the t in [0, 1] of the straight move c0 -> c1 where the edge is
   // within r of the centre (LinkWinSeg, closed form). Pure SI; scratch per thread, capacity reused
   inline void LinkWindowsMesh(const std::vector<SurfaceTri> &tris, const std::vector<std::array<int, 3>> &nb, const std::vector<int> &L, const V3 &c0, const V3 &c1, const double r,
                               std::vector<LinkWindow> &out)
   {
      thread_local std::vector<uint32_t> mark;   // stamp: a face in the list
      thread_local uint32_t stamp = 0;
      if (mark.size() < tris.size())
         mark.resize(tris.size(), 0u);
      if (++stamp == 0)
      {
         std::fill(mark.begin(), mark.end(), 0u);
         stamp = 1;
      }
      const V3 m = c1 - c0;
      double t0, t1;
      for (const int i : L)
         mark[i] = stamp;
      for (const int t : L)
      {
         const SurfaceTri &T = tris[t];
         if (!T.area)
            continue;
         const V3 *v[3] = { &T.a, &T.b, &T.c };
         for (int e = 0; e < 3; e++)
         {
            const int o = nb[t][e];
            if (o <= t || mark[o] != stamp || !tris[o].area)
               continue;
            if (LinkWinSeg(c0, m, *v[e], *v[(e + 1) % 3], r, t0, t1))
               out.push_back({ t, o, t0, t1 });
         }
      }
   }

   // A triangle mesh with shared-edge neighbours (nb[t][e] the triangle across edge e = (v e, v e+1), -1 none; vid the welded vertex ids): the pieces are the pressed triangles
   // (a face's inside nearest point, or the deepest of its edges where the nearest point is on the rim) and the shared edges between two of them (the humps).
   // One peak: ProbeTriangles' result unchanged. More: the level form, reff and the nearest point kept from the deepest.
   inline TubeProbe BarsTriangles(const BallSI &ball, const std::vector<SurfaceTri> &tris, const std::vector<std::array<int, 3>> &vid, const std::vector<std::array<int, 3>> &nb,
                                  const std::vector<int> &cand, const double edgeRho, const V3 &c, const double zero, const FaceRead *fr = nullptr,
                                  const std::vector<std::array<int8_t, 3>> *fold = nullptr, std::vector<TubeProbe> *patchesOut = nullptr)
   {
      const ContactProf::Ticks profBars(fr != nullptr && fr->t != nullptr && fr->t->tag == 1 ? ContactProf::CD().barsCorner : ContactProf::CD().bars);   // P3c profile (diagnostics only): the whole call, the probe inside included; the rubber corners apart
      if (ContactProf::CheckOn())
         ContactProf::CD().barsN++;
      // mc-cost profile (diagnostics only): stamps between the parts; the object's entry (named once its table has made its near pass)
      const bool profS = ContactProf::CheckOn();
      ContactProf::ObjBars *const profObj = profS && fr != nullptr && fr->t != nullptr ? ContactProf::ObjBarsOf(fr->t->objStats) : nullptr;
      struct ProfEnd
      {
         const bool on;
         ContactProf::ObjBars *const obj;
         const uint64_t t0;
         uint64_t cand = 0, nodes = 0;
         ~ProfEnd()
         {
            if (!on)
               return;
            const uint64_t dt = __rdtsc() - t0;
            ContactProf::CD().rowBarsN++;
            ContactProf::CD().rowBars += dt;
            if (obj != nullptr)
            {
               obj->t += dt;
               obj->n++;
               obj->cand += cand;
               obj->nodes += nodes;
            }
         }
      } profEnd { profS, profObj, profS ? __rdtsc() : 0 };
      profEnd.cand = cand.size();
      thread_local std::vector<FaceHit> hits;   // each candidate measured once, by the probe (scratch kept per thread)
      hits.clear();
      TubeProbe pr = ProbeTriangles(ball, tris, cand, edgeRho, c, 0.0, &hits, fr);
      bool pcs = false;   // patch_friction: the patch pieces written (else one() takes the probe's nearest piece)
      const auto one = [&]() {   // patch_contacts: one patch (or none): today's answer is its single entry
         if (patchesOut != nullptr)
         {
            patchesOut->clear();
            patchesOut->push_back(pr);
            if (!pcs)
            {
               PatchPieces().One(pr.overlap && pr.seg >= 0 ? pr.seg : PatchPieceList::kNoPiece, pr.delta, pr.normal);
               PatchPieceList &pl1 = PatchPieces();
               if (!pl1.pd.empty() && pr.seg >= 0 && pr.seg < (int)tris.size())   // (level_stiffness, B2: the face's triangle)
                  pl1.pd[0].g = FaceGeo(tris[pr.seg]);
            }
         }
         return pr;
      };
      uint64_t profT = profS ? __rdtsc() : 0;
      if (profS)
         ContactProf::CD().bProbe += profT - profEnd.t0;
      const auto profStamp = [&](uint64_t &to) {
         if (!profS)
            return;
         const uint64_t t = __rdtsc();
         to += t - profT;
         profT = t;
      };
      const std::vector<std::pair<int, double>> *const zl = PatchZeros();   // patch_zero: set, the caller passed zero 0 (Z2)
      if (!pr.overlap || pr.seg < 0 || pr.reff <= 0.0 || pr.delta <= zero || cand.size() < 2 || nb.size() != tris.size())
      {
         if (profS)
            ContactProf::CD().bRetProbe++;
         if (zl != nullptr && pr.overlap && pr.seg >= 0)   // patch_zero: the one piece, shifted by its patch's zero
            pr.delta -= PieceZero(*zl, pr.seg);
         return one();
      }
      if (profS)
         ContactProf::CD().bCand += cand.size();
      // cut B (runs/2026-11-12-mc-cutB/brief_build.md): the probe's result is the law's answer (one peak) when the faces within R (the near pass's set, cand) are one set joined through shared
      // edges and no shared edge between two of them is concave toward c (b on c's side of plane A and a on c's side of plane B; fold built at load, RubberSurface::Link). Then the face
      // pass, edge pass and merge are skipped. Check mode (VPX_FINDER_CHECK=1) runs the full law on every skippable call and counts any that ends with more than one peak.
      bool cbSkip = false;
      if (zl == nullptr && fold != nullptr && fold->size() == tris.size() && fr != nullptr && fr->t != nullptr && fr->t->useNear && &cand == &fr->t->inI)   // (patch_zero: the full law, each patch's pieces needed for its zero)
      {
         thread_local std::vector<uint32_t> mark;   // stamp: in the set; stamp + 1: reached
         thread_local std::vector<int> st;
         thread_local uint32_t stamp = 0;
         if (mark.size() < tris.size())
            mark.resize(tris.size(), 0u);
         stamp += 2;
         if (stamp < 2)
         {
            std::fill(mark.begin(), mark.end(), 0u);
            stamp = 2;
         }
         for (const int i : cand)
            mark[i] = stamp;
         cbSkip = true;
         for (size_t k = 0; k < cand.size() && cbSkip; k++)
         {
            const int i = cand[k];
            const SurfaceTri &A = tris[i];
            const double sA = (c - A.a).Dot(A.n);
            for (int e = 0; e < 3; e++)
            {
               const int o = nb[i][e];
               if (o <= i || mark[o] != stamp)
                  continue;
               const int fi = (*fold)[i][e];
               if (fi == 0)
                  continue;   // flat
               int fo = 2;
               for (int f = 0; f < 3; f++)
                  if (nb[o][f] == i)
                     fo = (*fold)[o][f];
               const double sB = (c - tris[o].a).Dot(tris[o].n);
               if (fi == 2 || fo == 2 || (fi * sA >= 0.0 && fo * sB >= 0.0))
               {
                  cbSkip = false;   // concave toward c (or unknown): a second patch is possible
                  break;
               }
            }
         }
         if (cbSkip)   // one set: everything reached from the first face through shared edges
         {
            st.clear();
            st.push_back(cand[0]);
            mark[cand[0]] = stamp + 1;
            size_t reached = 1;
            while (!st.empty())
            {
               const int i = st.back();
               st.pop_back();
               for (int e = 0; e < 3; e++)
               {
                  const int o = nb[i][e];
                  if (o >= 0 && mark[o] == stamp)
                  {
                     mark[o] = stamp + 1;
                     reached++;
                     st.push_back(o);
                  }
               }
            }
            cbSkip = reached == cand.size();
         }
      }
      const bool cbCheck = cbSkip && fr->t->check;
      if (profS)
      {
         if (cbSkip && !cbCheck)
            ContactProf::CD().cbSkip++;
         else
         {
            ContactProf::CD().cbFull++;
            if (profObj != nullptr)
               profObj->full++;
         }
      }
      // an edge's nearest point from its lower-id end, so both triangles on it get the same bits
      const auto edgeMeasure = [&](const int t, const int e, V3 &nOut) {
         const SurfaceTri &T = tris[t];
         const V3 *v[3] = { &T.a, &T.b, &T.c };
         const int a = e, b = (e + 1) % 3;
         const bool fwd = vid[t][a] < vid[t][b];
         const V3 q = fwd ? SegNearest(*v[a], *v[b], c) : SegNearest(*v[b], *v[a], c);
         const V3 dv = c - q;
         const double dist = dv.Length();
         nOut = dist > 1e-12 ? dv / dist : V3();
         return dist;
      };
      const auto edgeAxis = [&](const int t, const int e) {   // patch_friction (C4): edge e's unit axis when its nearest point to c lies inside it (n turns about it), else 0 (a vertex)
         const SurfaceTri &T = tris[t];
         const V3 *v[3] = { &T.a, &T.b, &T.c };
         const V3 ab = *v[(e + 1) % 3] - *v[e];
         const double l2 = ab.Dot(ab), s = l2 > 0.0 ? (c - *v[e]).Dot(ab) / l2 : 0.0;
         return s > 0.0 && s < 1.0 ? ab * (1.0 / std::sqrt(l2)) : V3();
      };
      FaceTable *const ft = fr != nullptr && fr->t != nullptr ? fr->t : nullptr;   // P3c: the edges' slots (no table: measured directly)
      const bool prof = ContactProf::CheckOn();
      const auto edgeDepth = [&](const int t, const int e, V3 &nOut) {
         FaceTable::ESlot *const sl = ft != nullptr ? ft->EdgeSlot(t, e, nb) : nullptr;
         if (prof)
            ContactProf::CD().eReads++;
         double dist;
         if (sl == nullptr)
         {
            dist = edgeMeasure(t, e, nOut);
            if (prof)
               ContactProf::CD().eFills++;
         }
         else
         {
            if (ft->count)
               ft->ereads++;
            if (sl->gen != ft->gen)
            {
               sl->dist = edgeMeasure(t, e, sl->n);
               sl->gen = ft->gen;
               if (ft->count)
                  ft->efills++;
               if (prof)
                  ContactProf::CD().eFills++;
            }
            else if (ft->check)
            {
               V3 n2;
               const double d2 = edgeMeasure(t, e, n2);
               if (std::memcmp(&d2, &sl->dist, sizeof(double)) != 0 || std::memcmp(&n2.x, &sl->n.x, sizeof(double)) != 0 || std::memcmp(&n2.y, &sl->n.y, sizeof(double)) != 0
                   || std::memcmp(&n2.z, &sl->n.z, sizeof(double)) != 0)
                  ft->emism++;
            }
            dist = sl->dist;
            nOut = sl->n;
         }
         return dist > 1e-12 ? ball.R - dist - zero : -1.0;
      };
      // a probed face's node: its depth above the zero d (<= 0: not pressed), normal nn, and the rim's deepest edge rimE (-1 the face's inside); the face pass and cut B
      const auto faceNode = [&](const FaceHit &fh, double &d, V3 &nn, int &rimE) {
         const SurfaceTri &T = tris[fh.i];
         d = -1.0;
         rimE = -1;
         if (std::abs((c - T.a).Dot(T.n)) >= ball.R - zero)
            return;   // every point of the face and its edges is at least the plane distance away: depth <= 0
         if (!T.area)
            return;   // zero area: no face (MakeSurfaceTri)
         if (fh.kind == 0)
         {
            const V3 dv = c - fh.p;
            const double dist = std::sqrt(fh.d2);   // = dv.Length()
            if (dist <= 1e-12)
               return;
            d = ball.R - dist - zero;
            nn = dv.Dot(T.n) < 0.0 ? T.n * -1.0 : T.n;
         }
         else
         {
            d = -1.0;   // the nearest point is on the rim: the deepest of the three edges
            for (int e = 0; e < 3; e++)
            {
               V3 ne;
               const double de = edgeDepth(fh.i, e, ne);
               if (de > d)
               {
                  d = de;
                  nn = ne;
                  rimE = e;
               }
            }
         }
      };
      if (cbSkip && !cbCheck)
      {
         if (patchesOut != nullptr)   // patch_friction (plan 2.1, 3.1): the shortcut still reports its pieces: the near set's PRESSED faces, from the probe's hits (each candidate measured once) by the face pass's depth (one patch, cut B's own test)
         {
            thread_local std::vector<int> pressed;
            thread_local std::vector<PieceData> pressedPd;
            pressed.clear();
            pressedPd.clear();
            for (const FaceHit &fh : hits)
            {
               double d;
               V3 nn;
               int rimE;
               faceNode(fh, d, nn, rimE);
               if (d <= 0.0)
                  continue;
               PieceData p;
               p.d = d + zero;
               p.n = nn;
               p.g = FaceGeo(tris[fh.i]);   // (B2: its triangle)
               pressed.push_back(fh.i);
               pressedPd.push_back(p);
            }
            // the pressed shared edges between two pressed faces, by the edge pass's own rule and data (edgeDepth: the face table), from the lower face: the merge tree's humps,
            // so B-lin's adjacency joins the faces (level_stiffness); the same PieceData as meshPd's edge node
            const size_t nFace = pressed.size();
            for (size_t x = 0; x < nFace; x++)
            {
               const int t = pressed[x];
               for (int e = 0; e < 3; e++)
               {
                  const int o = nb[t][e];
                  if (o <= t)
                     continue;
                  bool oOn = false;
                  for (size_t y = 0; y < nFace && !oOn; y++)
                     oOn = pressed[y] == o;
                  if (!oOn)
                     continue;
                  V3 ne;
                  const double de = edgeDepth(t, e, ne);
                  if (de <= 0.0)
                     continue;
                  const V3 *v[3] = { &tris[t].a, &tris[t].b, &tris[t].c };
                  PieceData p;
                  p.d = de + zero;
                  p.n = ne;
                  p.g.kind = 1;
                  p.g.a = *v[e];
                  p.g.b = *v[(e + 1) % 3];
                  p.nbA = t;
                  p.nbB = o;
                  pressed.push_back(-(1 + 3 * t + e));
                  pressedPd.push_back(p);
               }
            }
            if (PatchCheckEnv())   // VPX_PATCH_CHECK=1: the listed pieces through MergeBars (their depths, adjacency from the neighbour pairs): more than one peak counted (expected 0)
            {
               thread_local std::vector<BarNode> cn;
               thread_local std::vector<std::pair<int, int>> ca;
               thread_local std::vector<PeakTerm> keepTerms;
               cn.clear();
               ca.clear();
               for (size_t x = 0; x < pressed.size(); x++)
               {
                  const PieceData &p = pressedPd[x];
                  cn.push_back({ p.d - zero, p.n, x < nFace ? 0 : 1 });
                  for (size_t y = 0; y < x; y++)
                     if (p.nbA == pressed[y] || p.nbB == pressed[y] || pressedPd[y].nbA == pressed[x] || pressedPd[y].nbB == pressed[x])
                        ca.emplace_back((int)x, (int)y);
               }
               keepTerms = PeakTerms();   // (MergeBars writes the thread's peak list: the probe's own stays as it was)
               FoldSum cs;
               if (MergeBars(cn, ca, cs) > 1)
                  CutBCheck()[2]++;
               PeakTerms() = keepTerms;
            }
            if (PatchCheckEnv())   // VPX_PATCH_CHECK=1: against each near face probed alone (the set, any order)
            {
               thread_local std::vector<int> old, one1, a, b;
               old.clear();
               for (const int i : cand)
               {
                  one1.assign(1, i);
                  const TubeProbe q = ProbeTriangles(ball, tris, one1, edgeRho, c);
                  if (q.overlap && q.delta > zero)
                     old.push_back(i);
               }
               a.assign(pressed.begin(), pressed.begin() + (std::ptrdiff_t)nFace);   // (the faces only)
               b.assign(old.begin(), old.end());
               std::sort(a.begin(), a.end());
               std::sort(b.begin(), b.end());
               CutBCheck()[0]++;
               if (a != b)
                  CutBCheck()[1]++;
            }
            AllPatchPieces(pressed, pressedPd);
            PatchPieces().geo.assign(1, PieceGeo());
            pcs = true;
         }
         return one();
      }
      thread_local std::vector<BarNode> nodes;   // scratch kept per thread
      thread_local std::vector<std::pair<int, int>> adj;
      thread_local std::vector<std::pair<int, int>> triNode;   // (triangle, node)
      thread_local std::vector<int> triIdx;   // triangle -> its node, -1 none; every entry set here is put back to -1 before return
      thread_local std::vector<int> nodeTri;   // patch_contacts: node -> its triangle (a shared edge's node: its lower triangle)
      thread_local std::vector<int> nodePc;    // patch_friction: node -> its piece id (PatchPieceList)
      nodePc.clear();
      nodes.clear();
      adj.clear();
      triNode.clear();
      nodeTri.clear();
      if (triIdx.size() < tris.size())
         triIdx.resize(tris.size(), -1);
      for (const FaceHit &fh : hits)
      {
         const int i = fh.i;
         double d;
         V3 nn;
         int rimE;   // C4: the rim's deepest edge (its axis)
         faceNode(fh, d, nn, rimE);
         if (d <= 0.0)
            continue;
         if (triIdx[i] < 0)
            triIdx[i] = (int)nodes.size();   // the first node of a triangle (as the search it replaces)
         triNode.emplace_back(i, (int)nodes.size());
         nodes.push_back({ d, nn, 0, fh.kind == 0 ? 0.0 : ball.R - zero - d, rimE >= 0 ? edgeAxis(i, rimE) : V3() });   // (the rim: the edge's distance, edgeDepth's d = R - dist - zero)
         nodeTri.push_back(i);
         nodePc.push_back(i);
      }
      profStamp(ContactProf::CD().bFace);
      if (profS)
      {
         ContactProf::CD().bFaceNodes += nodes.size();
         profEnd.nodes = nodes.size();
      }
      const auto meshGeo = [&](const int x) {   // patch_friction: a node's piece geometry (PieceGeo)
         PieceGeo g;
         if (x < 0)
            return g;
         const int p = nodePc[x];
         if (p >= 0)
         {
            if (nodes[x].r > 0.0)   // a face's rim: its nearest point
            {
               g.kind = 2;
               g.a = c - nodes[x].n * nodes[x].r;
            }
            return g;
         }
         const int t = (-p - 1) / 3, e = (-p - 1) % 3;
         const V3 *v[3] = { &tris[t].a, &tris[t].b, &tris[t].c };
         g.kind = 1;
         g.a = *v[e];
         g.b = *v[(e + 1) % 3];
         return g;
      };
      thread_local std::vector<PieceData> nodePd;   // level_stiffness: each node's piece data (PieceData), filled where the pieces are written
      const auto meshPd = [&]() {
         nodePd.resize(nodes.size());
         for (size_t x = 0; x < nodes.size(); x++)
         {
            PieceData &p = nodePd[x];
            p.d = nodes[x].d + zero;
            p.n = nodes[x].n;
            const int q = nodePc[x];
            p.g = q >= 0 ? FaceGeo(tris[q]) : meshGeo((int)x);   // (B2: a face its triangle; an edge its segment)
            p.nbA = p.nbB = PatchPieceList::kNoPiece;
            if (q < 0)   // a shared edge: its two faces
            {
               const int t = (-q - 1) / 3, e = (-q - 1) % 3;
               p.nbA = t;
               p.nbB = nb[t][e];
            }
         }
      };
      if (nodes.size() < 2 && (zl == nullptr || nodes.empty()))   // (patch_zero: one node goes on to the shifted law)
      {
         if (zl != nullptr)   // patch_zero: no face pressed: the probe's piece, shifted by its patch's zero
            pr.delta -= PieceZero(*zl, pr.seg);
         if (patchesOut != nullptr)
         {
            meshPd();
            AllPatchPieces(nodePc, nodePd);
            PatchPieces().geo.assign(1, meshGeo(DeepestNode(nodes)));
            pcs = true;
         }
         for (const auto &tn : triNode)
            triIdx[tn.first] = -1;
         if (profS)
            ContactProf::CD().bRetNodes++;
         return one();
      }
      const auto nodeOf = [&](const int t) { return triIdx[t]; };
      for (const auto &tn : triNode)
         for (int e = 0; e < 3; e++)
         {
            const int o = nb[tn.first][e];
            if (o <= tn.first)   // each shared edge once (from its lower triangle)
               continue;
            const int on = nodeOf(o);
            if (on < 0)
               continue;
            V3 ne;
            const double de = edgeDepth(tn.first, e, ne);
            if (de <= 0.0)
               continue;   // the edge between them is not pressed: two separate contacts
            const int k = (int)nodes.size();
            nodes.push_back({ de, ne, 1, ball.R - zero - de, edgeAxis(tn.first, e) });
            nodeTri.push_back(tn.first);
            nodePc.push_back(-(1 + 3 * tn.first + e));
            adj.emplace_back(k, tn.second);
            adj.emplace_back(k, on);
         }
      for (const auto &tn : triNode)
         triIdx[tn.first] = -1;
      profStamp(ContactProf::CD().bEdge);
      if (profS)
         ContactProf::CD().bNodes += nodes.size();
      FoldSum sum;
      if (zl != nullptr)   // patch_zero (Z2): the raw patches, each shifted by its own zero; the merged answer and the patch list from the shifted level form
      {
         thread_local std::vector<ZeroPatch> zp;
         thread_local std::vector<double> nodeZ;
         const int peaks = ShiftedBars(nodes, adj, nodePc, *zl, sum, zp, nodeZ);
         profStamp(ContactProf::CD().bMerge);
         const int rawTop = DeepestNode(nodes);
         const TubeProbe base = pr;   // the raw deepest piece's probe
         const auto alone = [&](const int x) {   // a node's piece probed alone (its triangle)
            if (x == rawTop)
               return base;
            thread_local std::vector<int> one2;
            one2.assign(1, nodeTri[x]);
            return ProbeTriangles(ball, tris, one2, edgeRho, c);
         };
         int best = 0;   // the patch whose deepest node is deepest after the shift: the merged answer keeps its reff, point and seg
         for (int p = 1; p < (int)zp.size(); p++)
            if (zp[p].top > zp[best].top)
               best = p;
         pr = alone(zp[best].deep);
         if (peaks > 0)
            sum.Finish(pr, 0.0);
         else   // no node above its zero: the gap is the deepest shifted node's
         {
            pr.delta = zp[best].top;
            pr.normal = nodes[zp[best].deep].n;
            pr.grad = 1.0;
         }
         if (patchesOut != nullptr)
         {
            patchesOut->clear();
            for (const ZeroPatch &p : zp)
            {
               TubeProbe q = alone(p.deep);
               if (p.sum.s25 > 0.0)
                  p.sum.Finish(q, 0.0);
               else
               {
                  q.delta = p.top;
                  q.normal = nodes[p.deep].n;
                  q.grad = 1.0;
               }
               patchesOut->push_back(q);
            }
            meshPd();
            for (size_t x = 0; x < nodes.size(); x++)   // (B-lin: each piece's depth from its patch's zero)
               nodePd[x].d -= nodeZ[x];
            FillPatchPieces(nodePc, nodePd, (int)zp.size());
            for (const ZeroPatch &p : zp)
               PatchPieces().geo.push_back(meshGeo(p.deep));
            pcs = true;
         }
         return pr;
      }
      thread_local std::vector<PatchSum> patches;
      const int profPeaks = MergeBars(nodes, adj, sum, patchesOut != nullptr ? &patches : nullptr);
      profStamp(ContactProf::CD().bMerge);
      if (profPeaks <= 1)
      {
         if (patchesOut != nullptr)
         {
            meshPd();
            AllPatchPieces(nodePc, nodePd);
            PatchPieces().geo.assign(1, meshGeo(DeepestNode(nodes)));
            pcs = true;
         }
         if (profS)
            ContactProf::CD().bOnePeak++;
         return one();
      }
      if (cbCheck)   // cut B check mode: the rule would have returned the probe's single answer here
         fr->t->cbMism++;
      const TubeProbe base = pr;   // the deepest piece's probe (reff, point, seg)
      sum.Finish(pr, zero);
      profStamp(ContactProf::CD().bFinish);
      if (profS)
         ContactProf::CD().bMulti++;
      if (patchesOut != nullptr)
      {
         if (patches.size() <= 1)
         {
            meshPd();
            AllPatchPieces(nodePc, nodePd);   // (every node is in the one patch)
            PatchPieces().geo.assign(1, meshGeo(DeepestNode(nodes)));
            pcs = true;
            return one();   // one patch: today's merged answer
         }
         patchesOut->clear();
         int top = patches[0].deep;   // the deepest node of all (MergeBars' tie rule): its patch keeps the probe's reff, point and seg
         for (const PatchSum &p : patches)
            if (nodes[p.deep].d > nodes[top].d || (nodes[p.deep].d == nodes[top].d && p.deep < top))
               top = p.deep;
         thread_local std::vector<int> one1;
         for (const PatchSum &p : patches)
         {
            TubeProbe q = base;
            if (p.deep != top)   // the patch's deepest piece probed alone (its triangle)
            {
               one1.assign(1, nodeTri[p.deep]);
               q = ProbeTriangles(ball, tris, one1, edgeRho, c);
            }
            p.sum.Finish(q, zero);
            patchesOut->push_back(q);
         }
         meshPd();
         FillPatchPieces(nodePc, nodePd, (int)patches.size());
         for (const PatchSum &p : patches)
            PatchPieces().geo.push_back(meshGeo(p.deep));
      }
      return pr;
   }


   // Static-surface contact for one step. Same law as Evaluate (ContactCore), plus rolling resistance on the tangential part of the
   // ball's spin (the roll axis of a ball running along a post is arbitrary, unlike the flipper's vertical face). h = rubber layer thickness.
   // torqueOnBat is left 0 (nothing moves).
   // vs = velocity of the rubber surface at the contact (zero for a rigid static rubber, the node's motion for a band).
   inline ContactResult EvaluateTube(const RubberParams &rp, const BallSI &ball, const TubeProbe &pr, double hLayer,
                                     const V3 &v, const V3 &w, double dt, RubberState &st, const V3 &vs = V3(), const V3 &ws = V3())
   {
      ContactResult out;
      out.dist = pr.dist;
      if (!pr.overlap)
      {
         st.Reset();
         return out;
      }
      const CoreResult core = ContactCore(rp, ball, pr.normal, pr.delta, pr.reff, hLayer, vs, v, w, dt, st, ws);
      V3 Mb = core.torque;
      if ((rp.tanDelta > 0.0 || rp.rollingMu > 0.0) && core.fn > 0.0)
      {
         const V3 wt = w - pr.normal * w.Dot(pr.normal);
         const double wl = wt.Length();
         if (wl > 1e-9)
         {
            const double pa = core.fn / (M_PI * core.a * core.a);
            const double muRr = std::max(std::min(2.34 * (pa / core.E) * rp.tanDelta, kRollingMuMax), rp.rollingMu);
            Mb = Mb - wt * (muRr * ball.R * core.fn * std::tanh(wl / 2.0) / wl);
         }
      }
      out.inContact = true;
      out.force = core.force;
      out.torqueOnBall = Mb;
      out.normal = pr.normal;
      out.point = pr.point;
      out.delta = pr.delta;
      out.patchA = core.a;
      out.fn = core.fn;
      out.reff = pr.reff;
      return out;
   }
}
