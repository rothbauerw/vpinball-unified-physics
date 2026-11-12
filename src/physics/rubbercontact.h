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
#include <utility>
#include <cstring>
#include <cstdint>
#include "contactprof.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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
   inline TubeProbe ProbeTube(const BallSI &ball, const V3 *pts, int n, double rho, const V3 &c, const TubePre *pre = nullptr)
   {
      TubeProbe pr;
      pr.rho = rho;
      double bestD2 = INFINITY;   // squared: the nearest segment is found without a square root per segment
      V3 bestP;
      int bestI = -1;
      if (pre != nullptr && n >= 2)
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
   inline void PatchTube(const BallSI &ball, const V3 *pts, const int n, const double rho, const V3 &c, const TubePre &pre, const double zero, TubeProbe &pr)
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
      for (int i = 0; i + 1 < n; i++)
      {
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

   // A round hole of a non-legacy kicker in the playfield (centre and radius in the surface's plane, metres): a ball inside one is handed to the kicker, see InPlayfieldHole.
   struct SurfaceHole
   {
      double x, y, r;
   };

   // The nearest point on the faces in `cand`. A contact on a face is flat (effective radius R). On an edge or a corner, where the direction to the ball is
   // more than ~5 deg off the face normal, the rubber counts as rounded to edgeRho (the effective radius falls toward that of a ball on a tube of radius
   // edgeRho, fully at 15 deg), so a knife edge has a contact patch. The depth is measured from the sharp geometry. The result has rho = 0.
   // A round hole cut in a flat surface (a turntable's, metres): the surface within r of (x, y) at height z is missing; the nearest point of what is left is on the hole's edge.
   struct CarvedHole
   {
      double x, y, r, z;
      bool Holds(const V3 &p) const   // a point of a face inside the hole isn't surface
      {
         const double dx = p.x - x, dy = p.y - y;
         return std::abs(p.z - z) < 0.003 && dx * dx + dy * dy < r * r;
      }
   };

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
         uint32_t wgen;   // P3b: the walk's mark (visited at this table generation)
      };
      std::vector<Slot> s;
      uint32_t gen = 0;   // the table's generation: a slot is filled when its gen equals it
      uint64_t ver = ~0ull;
      double cx = 0.0, cy = 0.0, cz = 0.0;
      bool check = false;
      uint64_t fills = 0, reads = 0, mism = 0;   // faces measured, table reads, check-mode mismatches (summed by the owner)
      int lastMismJ = -1;
      uint64_t cyc = ~0ull;   // the counters: the cycle of the last contact check and the checks in it (hitrubber.cpp TabCheck)
      uint32_t nCyc = 0;
      bool count = false;    // P3a's per-read counters (fills, reads): check or profile mode only
      bool direct = false;   // VPX_FINDER_TABLE=0 (diagnostics only, results equal): the loops measure every face directly (the walk still keeps its own measurements here)
      int tag = -1;          // the counters: 1 = a rubber corner timed apart (design.md section 7), 0 = not, -1 = not yet known
      // shared_finder P3b (design.md section 4): the shared-edge walk, made once per table (ball position). inJ / inI: the faces within R of c, in load order (query-list positions,
      // face indices); seed: the faces within R at this ball's previous evaluation on this object (empty: none, the list with the box drop); useWalk: the loops read the walk (false: the
      // whole query list, check mode's comparison).
      uint32_t walkGen = 0;
      bool useWalk = false;
      std::vector<int> inJ, inI, seed, stack;
      // P3b2: the query list's index, rebuilt once per query version (one pass over the list, as the query's own gather): idx by face id gives the face's list position and its mesh
      // piece's local number (stamped, never cleared); the list's pieces with their faces' positions (pcJ from pcStart[lp]), so the walk's unreached check visits only pieces it did not reach.
      struct Idx
      {
         uint32_t stamp;
         int j, lp;
      };
      std::vector<Idx> idx;
      std::vector<int> pcLoc, pcStart, pcJ;
      std::vector<uint32_t> pcLocStamp, pcReached;
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
      int ListPos(const int i) const   // a face's query-list position from the walk's index (-1: no index for this query list, or not in it)
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
      uint64_t walkMiss = 0, walkBox = 0, walkVisit = 0;   // check mode: faces within R the walk did not hand over; box tests of unreached pieces and faces visited by the walk (counted with fills)
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
      // the walk: the seed measured, then across shared edges (nb) to each neighbour while it is within R; then every face of the list in a mesh piece the walk did not reach, dropped when
      // its box is farther than R, else measured. list: the object's query list (face indices, ascending).
      void Walk(const std::vector<SurfaceTri> &tris, const std::vector<std::array<int, 3>> &nb, const std::vector<int> &piece, const std::vector<int> &list, const V3 &c, const double R2)
      {
         if (walkGen == gen)
            return;
         walkGen = gen;
         inJ.clear();
         stack.clear();
         const bool linked = nb.size() == tris.size() && piece.size() == tris.size();
         const auto unreachedFace = [&](const int j) {   // a face of a piece the walk did not reach: dropped when its box is farther than R, else measured
            if (s[j].wgen == gen)
               return;
            const int i = list[j];
            if (count)
               walkBox++;
            if (BoxD2(tris[i], c) > R2 * (1.0 + kBoxSlack))
               return;
            s[j].wgen = gen;
            if (Get(j, tris[i], c).d2 < R2)
               inJ.push_back(j);
         };
         if (!linked)
         {
            for (size_t j = 0; j < list.size(); j++)
               unreachedFace((int)j);
         }
         else
         {
            if (!idxOk)   // a new query list: index it (once per query version)
            {
               idxOk = true;
               if (++idxStamp == 0)
               {
                  for (Idx &x : idx)
                     x.stamp = 0;
                  std::fill(pcLocStamp.begin(), pcLocStamp.end(), 0u);
                  idxStamp = 1;
               }
               if (idx.size() < tris.size())
                  idx.resize(tris.size(), Idx { 0, -1, -1 });
               pcStart.clear();
               for (size_t j = 0; j < list.size(); j++)
               {
                  const int i = list[j], pc = piece[i];
                  if (pc >= (int)pcLoc.size())
                  {
                     pcLoc.resize(pc + 1, -1);
                     pcLocStamp.resize(pc + 1, 0u);
                  }
                  if (pcLocStamp[pc] != idxStamp)
                  {
                     pcLocStamp[pc] = idxStamp;
                     pcLoc[pc] = (int)pcStart.size();
                     pcStart.push_back(0);
                  }
                  idx[i] = Idx { idxStamp, (int)j, pcLoc[pc] };
                  pcStart[pcLoc[pc]]++;
               }
               const size_t np = pcStart.size();   // counts to starts, then the positions grouped by piece (ascending within a piece)
               int sum = 0;
               for (size_t p = 0; p < np; p++)
               {
                  const int n = pcStart[p];
                  pcStart[p] = sum;
                  sum += n;
               }
               pcStart.push_back(sum);
               pcJ.resize(list.size());
               stack.assign(pcStart.begin(), pcStart.end() - 1);   // fill cursors (stack is empty again before the walk)
               for (size_t j = 0; j < list.size(); j++)
                  pcJ[stack[idx[list[j]].lp]++] = (int)j;
               stack.clear();
               pcReached.assign(np, 0u);
            }
            const auto visit = [&](const int i) {
               if (count)
                  walkVisit++;
               const Idx &x = idx[i];
               if (x.stamp != idxStamp)   // not in the query list
                  return;
               const int j = x.j;
               if (s[j].wgen == gen)
                  return;
               s[j].wgen = gen;
               if (Get(j, tris[i], c).d2 < R2)
               {
                  inJ.push_back(j);
                  stack.push_back(i);
                  pcReached[x.lp] = gen;
               }
            };
            for (const int i : seed)
               visit(i);
            while (!stack.empty())
            {
               const int i = stack.back();
               stack.pop_back();
               for (int e = 0; e < 3; e++)
                  if (nb[i][e] >= 0)
                     visit(nb[i][e]);
            }
            for (size_t p = 0; p + 1 < pcStart.size(); p++)   // only the pieces the walk did not reach: a reached piece costs nothing per face
               if (pcReached[p] != gen)
                  for (int k = pcStart[p]; k < pcStart[p + 1]; k++)
                     unreachedFace(pcJ[k]);
         }
         std::sort(inJ.begin(), inJ.end());   // load order (R9): the loops' tie rules see the faces in today's order
         inI.resize(inJ.size());
         for (size_t k = 0; k < inJ.size(); k++)
            inI[k] = list[inJ[k]];
         seed = inI;
         if (check)   // every face of the list measured fresh: one within R that the walk did not hand over is a miss
            for (size_t j = 0; j < list.size(); j++)
            {
               int k = 0;
               const V3 dv = c - ClosestOnTriangle(tris[list[j]], c, k);
               if (dv.Dot(dv) < R2 && !std::binary_search(inJ.begin(), inJ.end(), (int)j))
               {
                  walkMiss++;
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
                  x.gen = x.wgen = 0;
               for (ESlot &x : es)
                  x.gen = 0;
               std::fill(pcReached.begin(), pcReached.end(), 0u);
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
            s.resize(n, Slot { V3(), 0.0, 0, 0, 0 });
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
         if (t == nullptr || t->direct)   // VPX_FINDER_TABLE=0: measured directly, as with no table
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
                                   const CarvedHole *holes = nullptr, int nHoles = 0, double halfThick = 0.0, std::vector<FaceHit> *hits = nullptr, const FaceRead *fr = nullptr)
   {
      TubeProbe pr;
      double bestD2 = INFINITY;
      V3 bestP;
      int bestI = -1, bestKind = 0;
      { const ContactProf::LoopTicks profLoop(cand.size());   // P3a profile (diagnostics only)
      for (size_t j = 0; j < cand.size(); j++)
      {
         const int i = cand[j];
         if (fr != nullptr && fr->prune && nHoles == 0 && bestI >= 0 && FaceTable::BoxD2(tris[i], c) > (bestD2 + 1e-18) * (1.0 + FaceTable::kBoxSlack))
            continue;   // P3b: cannot be nearer than the best so far nor tie with it (exact lower bound; the 1e-18 tie band kept)
         int kind = 0;
         V3 p;
         const FaceTable::Slot *sl = nullptr;
         if (fr != nullptr)   // shared_finder P3a: the table's measurement (raw; the carved hole below is applied to it as before)
         {
            sl = &fr->At(j, tris[i], c);
            p = sl->p;
            kind = sl->kind;
         }
         else
            p = ClosestOnTriangle(tris[i], c, kind);
         for (int k = 0; k < nHoles; k++)
         {
            const CarvedHole &h = holes[k];
            if (h.Holds(p))
            {
               double ux = c.x - h.x, uy = c.y - h.y;
               const double ul = std::sqrt(ux * ux + uy * uy);
               if (ul > 1e-9) { ux /= ul; uy /= ul; } else { ux = 1.0; uy = 0.0; }
               p = V3(h.x + ux * h.r, h.y + uy * h.r, p.z);   // the nearest point of the surface that is left: on the rim of the hole
               kind = 1;
               sl = nullptr;   // the point moved: its distance is measured from it
            }
         }
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
   inline void PatchTriangles(const BallSI &ball, const std::vector<SurfaceTri> &tris, const std::vector<int> &cand, const V3 &c, const double zero, TubeProbe &pr, const CarvedHole *holes = nullptr, int nHoles = 0,
                              const FaceRead *fr = nullptr)
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
         if (fr != nullptr)   // shared_finder P3a: the table's measurement (the Holds test below applied to it as before)
         {
            const FaceTable::Slot &sl = fr->At(j, t, c);
            p = sl.p;
            kind = sl.kind;
         }
         else
            p = ClosestOnTriangle(t, c, kind);
         if (kind != 0 || std::any_of(holes, holes + nHoles, [&](const CarvedHole &h) { return h.Holds(p); }))
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
   };

   // nodes: the pieces pressing the ball (d > 0); adj: pairs of nodes that touch on the surface. Returns the number of peaks; sum holds their s25 and g.
   inline int MergeBars(const std::vector<BarNode> &nodes, const std::vector<std::pair<int, int>> &adj, FoldSum &sum)
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
            passed[pl] = 1;
            parent[lose] = keep;
         }
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
         if (passed[i] && pass[i] > 0.0)
         {
            const double p = pass[i], p15 = p * std::sqrt(p);
            sum.s25 -= p * p15;
            sum.g = sum.g - passN[i] * p15;
         }
         peaks++;
      }
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

   // A tube (centreline pts, n points; a closed ring repeats its first point at the end, the two are one joint): the pieces are the segments' inside nearest points and the joints.
   // One peak: ProbeTube's result unchanged. More: the level form above, reff and the nearest point kept from the deepest (held reff, piece 2).
   inline TubeProbe BarsTube(const BallSI &ball, const V3 *pts, const int n, const double rho, const V3 &c, const TubePre &pre, const double zero)
   {
      TubeProbe pr = ProbeTube(ball, pts, n, rho, c, &pre);
      if (!pr.overlap || pr.seg < 0 || pr.reff <= 0.0 || pr.delta <= zero || n < 3)
         return pr;
      const double S = ball.R + rho;
      const bool ring = pts[0].x == pts[n - 1].x && pts[0].y == pts[n - 1].y && pts[0].z == pts[n - 1].z;
      const int nv = ring ? n - 1 : n;   // joints: point n-1 of a ring is joint 0
      thread_local std::vector<BarNode> nodes;   // scratch kept per thread
      thread_local std::vector<std::pair<int, int>> adj;
      thread_local std::vector<int> jointNode;
      nodes.clear();
      adj.clear();
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
               nodes.push_back({ d, dv / dist, 1 });
            }
         }
         return jointNode[j];
      };
      for (int i = 0; i + 1 < n; i++)
      {
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
               nodes.push_back({ d, dv / dist, 0 });
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
      if (MergeBars(nodes, adj, sum) <= 1)
         return pr;
      sum.Finish(pr, zero);
      return pr;
   }

   // A triangle mesh with shared-edge neighbours (nb[t][e] the triangle across edge e = (v e, v e+1), -1 none; vid the welded vertex ids): the pieces are the pressed triangles
   // (a face's inside nearest point, or the deepest of its edges where the nearest point is on the rim) and the shared edges between two of them (the humps).
   // One peak: ProbeTriangles' result unchanged. More: the level form, reff and the nearest point kept from the deepest.
   inline TubeProbe BarsTriangles(const BallSI &ball, const std::vector<SurfaceTri> &tris, const std::vector<std::array<int, 3>> &vid, const std::vector<std::array<int, 3>> &nb,
                                  const std::vector<int> &cand, const double edgeRho, const V3 &c, const double zero, const FaceRead *fr = nullptr,
                                  const std::vector<std::array<int8_t, 3>> *fold = nullptr)
   {
      const ContactProf::Ticks profBars(fr != nullptr && fr->t != nullptr && fr->t->tag == 1 ? ContactProf::CD().barsCorner : ContactProf::CD().bars);   // P3c profile (diagnostics only): the whole call, the probe inside included; the rubber corners apart
      if (ContactProf::CheckOn())
         ContactProf::CD().barsN++;
      // mc-cost profile (diagnostics only): stamps between the parts; the object's entry (named once its table has walked)
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
      TubeProbe pr = ProbeTriangles(ball, tris, cand, edgeRho, c, nullptr, 0, 0.0, &hits, fr);
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
      if (!pr.overlap || pr.seg < 0 || pr.reff <= 0.0 || pr.delta <= zero || cand.size() < 2 || nb.size() != tris.size())
      {
         if (profS)
            ContactProf::CD().bRetProbe++;
         return pr;
      }
      if (profS)
         ContactProf::CD().bCand += cand.size();
      // cut B (runs/2026-11-12-mc-cutB/brief_build.md): the probe's result is the law's answer (one peak) when the faces within R (the walk's set, cand) are one set joined through shared
      // edges and no shared edge between two of them is concave toward c (b on c's side of plane A and a on c's side of plane B; fold built at load, RubberSurface::Link). Then the face
      // pass, edge pass and merge are skipped. Check mode (VPX_FINDER_CHECK=1) runs the full law on every skippable call and counts any that ends with more than one peak.
      bool cbSkip = false;
      if (fold != nullptr && fold->size() == tris.size() && fr != nullptr && fr->t != nullptr && fr->t->useWalk && &cand == &fr->t->inI)
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
      if (cbSkip && !cbCheck)
         return pr;
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
      FaceTable *const ft = fr != nullptr && fr->t != nullptr && !fr->t->direct ? fr->t : nullptr;   // P3c: the edges' slots (VPX_FINDER_TABLE=0 or no table: measured directly)
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
      thread_local std::vector<BarNode> nodes;   // scratch kept per thread
      thread_local std::vector<std::pair<int, int>> adj;
      thread_local std::vector<std::pair<int, int>> triNode;   // (triangle, node)
      thread_local std::vector<int> triIdx;   // triangle -> its node, -1 none; every entry set here is put back to -1 before return
      nodes.clear();
      adj.clear();
      triNode.clear();
      if (triIdx.size() < tris.size())
         triIdx.resize(tris.size(), -1);
      for (const FaceHit &fh : hits)
      {
         const int i = fh.i;
         const SurfaceTri &T = tris[i];
         if (std::abs((c - T.a).Dot(T.n)) >= ball.R - zero)
            continue;   // every point of the face and its edges is at least the plane distance away: depth <= 0
         if (!T.area)
            continue;   // zero area: no face (MakeSurfaceTri)
         double d;
         V3 nn;
         if (fh.kind == 0)
         {
            const V3 dv = c - fh.p;
            const double dist = std::sqrt(fh.d2);   // = dv.Length()
            if (dist <= 1e-12)
               continue;
            d = ball.R - dist - zero;
            nn = dv.Dot(T.n) < 0.0 ? T.n * -1.0 : T.n;
         }
         else
         {
            d = -1.0;   // the nearest point is on the rim: the deepest of the three edges
            for (int e = 0; e < 3; e++)
            {
               V3 ne;
               const double de = edgeDepth(i, e, ne);
               if (de > d)
               {
                  d = de;
                  nn = ne;
               }
            }
         }
         if (d <= 0.0)
            continue;
         if (triIdx[i] < 0)
            triIdx[i] = (int)nodes.size();   // the first node of a triangle (as the search it replaces)
         triNode.emplace_back(i, (int)nodes.size());
         nodes.push_back({ d, nn, 0 });
      }
      profStamp(ContactProf::CD().bFace);
      if (profS)
      {
         ContactProf::CD().bFaceNodes += nodes.size();
         profEnd.nodes = nodes.size();
      }
      if (nodes.size() < 2)
      {
         for (const auto &tn : triNode)
            triIdx[tn.first] = -1;
         if (profS)
            ContactProf::CD().bRetNodes++;
         return pr;
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
            nodes.push_back({ de, ne, 1 });
            adj.emplace_back(k, tn.second);
            adj.emplace_back(k, on);
         }
      for (const auto &tn : triNode)
         triIdx[tn.first] = -1;
      profStamp(ContactProf::CD().bEdge);
      if (profS)
         ContactProf::CD().bNodes += nodes.size();
      FoldSum sum;
      const int profPeaks = MergeBars(nodes, adj, sum);
      profStamp(ContactProf::CD().bMerge);
      if (profPeaks <= 1)
      {
         if (profS)
            ContactProf::CD().bOnePeak++;
         return pr;
      }
      if (cbCheck)   // cut B check mode: the rule would have returned the probe's single answer here
         fr->t->cbMism++;
      sum.Finish(pr, zero);
      profStamp(ContactProf::CD().bFinish);
      if (profS)
         ContactProf::CD().bMulti++;
      return pr;
   }

   // A probe reports the surface's nearest point only, so a part of the surface that is not the nearest (a ramp's wall while the ball flies just over its floor) has no say in the step and the
   // ball can be carried deep into it, with stored energy no work paid for. These give the longest step, at most hMax, that brings the ball no closer than arriveM to any other part of the surface
   // (the arrival rule, part by part: d <= vc h + a h^2). A part within touchM counts as touching. The nearest part (nearest: the probe's seg) is the caller's.
   inline double ArrivalStep(const double d, const double vc, const double a, const double h)
   {
      if (d > vc * h + a * h * h)
         return h;
      const double t = a > 1e-9 ? (-vc + std::sqrt(vc * vc + 4.0 * a * d)) / (2.0 * a) : vc > 1e-12 ? d / vc : h;
      return std::min(h, t);
   }

   // a face in the plane of the nearest face has the same gap and closing speed: it is passed over; so is one whose plane the ball cannot reach within the step
   inline double ArrivalTriangles(const BallSI &ball, const std::vector<SurfaceTri> &tris, const std::vector<int> &cand, const int nearest, const V3 &c, const V3 &v, const double a, const double hMax, const double touchM, const double arriveM,
                                  const FaceRead *fr = nullptr)
   {
      double h = hMax;
      const double speed = v.Length();
      const SurfaceTri *const best = nearest >= 0 ? &tris[nearest] : nullptr;
      const ContactProf::LoopTicks profLoop(cand.size());   // P3a profile (diagnostics only)
      for (size_t j = 0; j < cand.size(); j++)
      {
         const int i = cand[j];
         if (i == nearest)
            continue;
         const SurfaceTri &t = tris[i];
         if (std::abs((c - t.a).Dot(t.n)) - ball.R - touchM > speed * h + a * h * h)
            continue;
         if (best != nullptr && std::abs(t.n.Dot(best->n)) > 0.9999 && std::abs((t.a - best->a).Dot(best->n)) < 1e-5)
            continue;
         int kind = 0;
         const V3 dv = c - (fr != nullptr ? fr->At(j, t, c).p : ClosestOnTriangle(t, c, kind));   // shared_finder P3a: the table's measurement
         const double dist = dv.Length();
         if (dist <= 1e-12 || dist - ball.R <= touchM)
            continue;
         h = ArrivalStep(dist - ball.R - arriveM, std::max(-v.Dot(dv) / dist, 0.0), a, h);
      }
      return h;
   }

   // the segments next to the nearest one are passed over: the nearest point runs onto them with no jump in the gap
   inline double ArrivalTube(const BallSI &ball, const V3 *pts, const int n, const double rho, const TubePre &pre, const int nearest, const V3 &c, const V3 &v, const double a, const double hMax, const double touchM, const double arriveM)
   {
      double h = hMax;
      const double speed = v.Length(), S = ball.R + rho;
      for (int i = 0; i + 1 < n; i++)
      {
         if (std::abs(i - nearest) <= 1 || pre.L2[i] <= 0.0)
            continue;
         if ((c - pre.mid[i]).Length() - pre.rad[i] - S - touchM > speed * h + a * h * h)
            continue;
         const double s = std::min(std::max((c - pts[i]).Dot(pre.t[i]) / pre.L2[i], 0.0), 1.0);
         const V3 dv = c - (pts[i] + pre.t[i] * s);
         const double dist = dv.Length();
         if (dist <= 1e-12 || dist - S <= touchM)
            continue;
         h = ArrivalStep(dist - S - arriveM, std::max(-v.Dot(dv) / dist, 0.0), a, h);
      }
      return h;
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
