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

   struct RubberState
   {
      V3 xi;                 // shear spring state (tangent plane), m; zero when not in contact
      void Reset() { xi = V3(); }
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

   // How much of a damper's force D goes on top of the spring's force Fs and stays within the limit lim (|Fs| <= lim): the largest s in [0, 1] with |Fs + s D| <= lim.
   // Taken off the damper, so the spring keeps its own force and the damper can only remove energy.
   inline double DamperShare(const V3 &Fs, const V3 &D, const double lim)
   {
      const V3 sum = Fs + D;
      if (sum.Dot(sum) <= lim * lim)
         return 1.0;
      const double qa = D.Dot(D), qb = 2.0 * Fs.Dot(D), qc = Fs.Dot(Fs) - lim * lim;
      return qa > 0.0 ? std::min(1.0, std::max(0.0, (-qb + std::sqrt(std::max(0.0, qb * qb - 4.0 * qa * qc))) / (2.0 * qa))) : 1.0;
   }

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
      {
         const double l0 = st.xi.Length();
         const V3 xp = st.xi - n * st.xi.Dot(n);
         const double lp = xp.Length();
         st.xi = (lp > 1e-300 ? xp * (l0 / lp) : V3()) + vt * dt;   // turned into the contact plane, its length kept
      }
      const double mt = 1.0 / (1.0 / m + R * R / Ib);          // tangential effective mass (m/3.5 for a solid sphere)
      const double ct = rp.shearTan * std::sqrt(kt * mt);
      const V3 Fs = st.xi * (-kt);
      V3 Ft;
      const double lim = rp.mu * Fn;
      const double ns = Fs.Length();
      if (ns > lim)
      {
         Ft = Fs * (lim / ns);
         st.xi = Ft * (-1.0 / kt);   // on the spring's force alone: a damper's share here would charge the spring with work nobody did
      }
      else
      {
         const V3 D = vt * (-ct);
         Ft = Fs + D * DamperShare(Fs, D, lim);
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

   // Surface curvature along the tube at the contact (signed, 1/m, > 0 convex): for a centreline of curvature k the outer equator gives
   // cos/(1/k + rho cos), the inner one is concave. cosang = +1 outer equator, -1 inner (in-plane direction of the normal against the centre of curvature).
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
         const double nh = std::sqrt(pr.normal.x * pr.normal.x + pr.normal.y * pr.normal.y);
         if (ok && nh > 1e-9)
         {
            const double cosang = -(pr.normal.x * dnx + pr.normal.y * dny) / nh;
            kLong = TubeAlongCurvature(k, rho, cosang);
         }
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
   };

   inline SurfaceTri MakeSurfaceTri(const V3 &a, const V3 &b, const V3 &c)
   {
      SurfaceTri t;
      t.a = a;
      t.b = b;
      t.c = c;
      const V3 cr = (b - a).Cross(c - a);
      const double l = cr.Length();
      t.n = l > 0.0 ? cr / l : V3(0.0, 0.0, 1.0);
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
   };

   // halfThick > 0: the faces are the mid-surface of a sheet of thickness 2 halfThick (a gate plate): the surface is everything within halfThick of the faces
   // (flat on both sides, edges rounded to radius halfThick), so the depth is R + halfThick - distance. Pass edgeRho = halfThick for such a sheet.
   inline TubeProbe ProbeTriangles(const BallSI &ball, const std::vector<SurfaceTri> &tris, const std::vector<int> &cand, double edgeRho, const V3 &c,
                                   const CarvedHole *holes = nullptr, int nHoles = 0, double halfThick = 0.0)
   {
      TubeProbe pr;
      double bestD2 = INFINITY;
      V3 bestP;
      int bestI = -1, bestKind = 0;
      for (const int i : cand)
      {
         int kind = 0;
         V3 p = ClosestOnTriangle(tris[i], c, kind);
         for (int k = 0; k < nHoles; k++)
         {
            const CarvedHole &h = holes[k];
            const double dx = p.x - h.x, dy = p.y - h.y;
            if (std::abs(p.z - h.z) < 0.003 && dx * dx + dy * dy < h.r * h.r)
            {
               double ux = c.x - h.x, uy = c.y - h.y;
               const double ul = std::sqrt(ux * ux + uy * uy);
               if (ul > 1e-9) { ux /= ul; uy /= ul; } else { ux = 1.0; uy = 0.0; }
               p = V3(h.x + ux * h.r, h.y + uy * h.r, p.z);   // the nearest point of the surface that is left: on the rim of the hole
               kind = 1;
            }
         }
         const V3 dv = c - p;
         const double d2 = dv.Dot(dv);
         // on a tie prefer the face over an edge or a corner (the diagonal of a flat face is an "edge" of both its triangles)
         if (d2 < bestD2 - 1e-18 || (d2 < bestD2 + 1e-18 && kind < bestKind))
         {
            bestD2 = d2;
            bestP = p;
            bestI = i;
            bestKind = kind;
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
   inline void PatchTriangles(const BallSI &ball, const std::vector<SurfaceTri> &tris, const std::vector<int> &cand, const V3 &c, const double zero, TubeProbe &pr)
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
      for (const int i : cand)
      {
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
         const V3 p = ClosestOnTriangle(t, c, kind);
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
      if (any)
         sum.Finish(pr, zero);
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
   inline double ArrivalTriangles(const BallSI &ball, const std::vector<SurfaceTri> &tris, const std::vector<int> &cand, const int nearest, const V3 &c, const V3 &v, const double a, const double hMax, const double touchM, const double arriveM)
   {
      double h = hMax;
      const double speed = v.Length();
      const SurfaceTri *const best = nearest >= 0 ? &tris[nearest] : nullptr;
      for (const int i : cand)
      {
         if (i == nearest)
            continue;
         const SurfaceTri &t = tris[i];
         if (std::abs((c - t.a).Dot(t.n)) - ball.R - touchM > speed * h + a * h * h)
            continue;
         if (best != nullptr && std::abs(t.n.Dot(best->n)) > 0.9999 && std::abs((t.a - best->a).Dot(best->n)) < 1e-5)
            continue;
         int kind = 0;
         const V3 dv = c - ClosestOnTriangle(t, c, kind);
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
