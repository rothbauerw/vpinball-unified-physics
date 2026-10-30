// license:GPLv3+
//
// Ball / rubber contact force model for the model flipper (unified engine, rubber_contact=1).
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
//   shear     Mindlin-style spring kt = 8 G* a with damping ratio zeta, Coulomb cap mu*Fn; the spring state xi is
//             kept in RubberState and re-seated at the cap
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
      double shearScale = 1.0;       // multiplier on the shear stiffness
      double shearZeta = 0.2;        // shear damping ratio
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
      st.xi = st.xi - n * st.xi.Dot(n) + vt * dt;
      const double mt = 1.0 / (1.0 / m + R * R / Ib);          // tangential effective mass (m/3.5 for a solid sphere)
      const double ct = 2.0 * rp.shearZeta * std::sqrt(kt * mt);
      V3 Ft = st.xi * (-kt) - vt * ct;
      const double lim = rp.mu * Fn;
      const double nt = Ft.Length();
      if (nt > lim)
      {
         Ft = Ft * (lim / nt);
         st.xi = Ft * (-1.0 / kt);
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

   // ---- Band: a tensioned elastic cord stretched between two posts (a free span), loaded at nodes ----
   //
   // The span between its two anchors A and B (the tangent points on the post wraps, fixed) is a cord of cross-section area pi rho^2,
   // Young's modulus E, stretched by `prestretch` over the posts: unstretched length L0 = |AB| / (1 + prestretch). Balls (and other spans) load it at
   // material points: each is a lumped node P_i at the fraction f_i of the unstretched length, carrying the mass of the cord around it. The cord is the
   // chain A - P_1 - ... - P_M - B of straight legs; a node can be added anywhere at any time (its place on the current shape, the leg it splits keeps
   // its strain). The band sticks to a ball (rubber on steel), so each leg keeps its own unstretched length (f_{j+1} - f_j) L0.
   //   leg tension   T_j = max(0, EA e_j + eta de_j/dt),  e_j = L_j / l0_j - 1,  eta = EA tan(delta) / (2 pi refHz)   (Kelvin-Voigt)
   //   force on P_i  T_{i+1} u_{i+1} - T_i u_i (u = unit vector along the leg, away from A), plus the contact reactions (a share of a contact at the
   //                 fraction s of a leg goes to each end of the leg; the part that lands on an anchor is carried by the post), plus a damping of the
   //                 node's own vibration (ratio zeta of the critical damping of its two legs)
   // Contact with a ball is the same tube contact as a rigid rubber (ContactCore) on the deformed centreline, with the leg's velocity as the surface
   // velocity. The shape is straight legs: no wave propagation along a leg and no slip over a ball, which is fine for the low speeds and the 100 mm spans
   // of a pinball band; the band cannot push (T >= 0).

   // A band stop: a post behind the band that the cord cannot enter. A circle in the band's plane (m); r already includes the cord's radius.
   struct StopCircle
   {
      double x = 0.0, y = 0.0, r = 0.0;
   };

   struct BandParams
   {
      std::vector<StopCircle> stops;   // band_stops; a stop that the span's rest line already touches is left out of that span
      double eMpa = 3.0;           // Young's modulus of the band material (the rubber's own e_mpa)
      double prestretch = 0.07;    // strain of the span at rest
      double tensionN = 0.0;       // rest tension of a span, N; > 0 sets the cord's axial stiffness EA = tensionN / prestretch (the modulus eMpa is then not used for the cord)
      double tanDelta = 0.1;       // loss tangent of the band material
      double refHz = 20.0;         // frequency at which tanDelta applies to the stretch (a band bounce lasts 25 - 50 ms)
      double nodeMassKg = 0.0;     // 0 = from the cord the node carries; else the mass of a node that carries the whole span (scaled by the share it carries)
      double densityKgM3 = 1100.0; // rubber
      bool lockZ = true;           // the nodes stay at the band's height (the posts hold the band there); else they move freely in 3D
      double zeta = 0.05;          // damping ratio of a node's own vibration (internal friction and air; the stretch viscosity only damps a taut cord at second order)
   };

   struct BandChain
   {
      struct Node
      {
         double f = 0.5;            // material place along the span, 0 at A, 1 at B
         V3 P0, P, vp;              // rest position (on the chord), position, velocity
         double m = 1e-3;           // kg
         double peak = 0.0;         // peak deflection since the chain was built, m
         bool alive = true;         // false = merged away (its slot is reused by Add)
         bool pinned = false;       // fixed to a stop (the band lies against a post): a corner of the cord, it does not move until the forces pull it off the post
         int pinStop = -1;          // which stop it is pinned to
         double Deflection() const { return (P - P0).Length(); }
      };

      V3 A, B;
      double Ls = 1.0;             // |AB|, m
      double L0 = 1.0;             // unstretched length of the whole span, m
      double EA = 0.0;             // axial stiffness, N
      double eta = 0.0;            // viscous tension coefficient, N s
      double rhoA = 0.0;           // mass of the cord per metre of its stretched length, kg/m
      double nodeMassKg = 0.0;
      double zeta = 0.05;
      bool lockZ = true;
      std::vector<StopCircle> stops;   // obstacles for this span (see ApplyStops)
      std::vector<Node> nodes;     // by id (the order they were added in)
      std::vector<int> order;      // node ids sorted by f
      std::vector<double> T;       // leg tensions of the last step, N (Count() + 1 legs)
      std::vector<V3> uBuf;        // scratch for Step: leg directions and lengths (kept between calls, so a step does not allocate)
      std::vector<double> lBuf;
      std::vector<double> mBuf;    // scratch for Add / Remove: the masses before the re-lumping

      void Init(const V3 &a, const V3 &b, double tubeRadius, const BandParams &bp)
      {
         A = a;
         B = b;
         Ls = (b - a).Length();
         L0 = Ls / (1.0 + bp.prestretch);
         const double area = M_PI * tubeRadius * tubeRadius;
         EA = (bp.tensionN > 0.0 && bp.prestretch > 1e-6) ? bp.tensionN / bp.prestretch : bp.eMpa * 1e6 * area;
         eta = EA * bp.tanDelta / (2.0 * M_PI * std::max(bp.refHz, 1.0));
         rhoA = bp.densityKgM3 * area;
         nodeMassKg = bp.nodeMassKg;
         zeta = bp.zeta;
         lockZ = bp.lockZ;
         nodes.clear();
         order.clear();
         T.clear();
         stops.clear();
         for (const StopCircle &s : bp.stops)
         {
            // a stop that contains one of this span's anchors is the anchor's own post, not an obstacle. A stop that the rest line merely touches or
            // crosses stays: the band lies against it (a post between the two sides of a loop), and a node placed on the rest line is put on its surface.
            if (std::hypot(s.x - a.x, s.y - a.y) >= s.r && std::hypot(s.x - b.x, s.y - b.y) >= s.r)
               stops.push_back(s);
         }
      }

      // A stop that the rest line touches or crosses is a post the band lies against. The band is fixed to it like a corner post: a pinned node on the post's
      // surface where the rest line meets it, so the span is two independent spans (anchor to post, post to anchor) with the cord bending at the post. A ball
      // that lands on the post loads this node (a plain fixed rubber), a ball anywhere else deflects its own side. Called once, when the chain is made.
      void PinStops()
      {
         for (size_t si = 0; si < stops.size() && Count() < 12; si++)
         {
            const StopCircle &s = stops[si];
            const double abx = B.x - A.x, aby = B.y - A.y, l2 = abx * abx + aby * aby;
            if (l2 < 1e-18)
               continue;
            const double t = ((s.x - A.x) * abx + (s.y - A.y) * aby) / l2;
            if (t < 0.05 || t > 0.95)
               continue;
            const double fx = A.x + abx * t, fy = A.y + aby * t;   // the point of the rest line nearest to the stop's centre
            const double dx = fx - s.x, dy = fy - s.y, d = std::hypot(dx, dy);
            if (d >= s.r)
               continue;   // the band does not touch it at rest: the stop acts only when the cord is pushed against it
            const double nx = d > 1e-9 ? dx / d : -aby / std::sqrt(l2), ny = d > 1e-9 ? dy / d : abx / std::sqrt(l2);
            const int id = Add(t);
            Node &n = nodes[id];
            n.P = V3(s.x + nx * s.r, s.y + ny * s.r, n.P.z);
            n.vp = V3();
            n.pinned = true;
            n.pinStop = (int)si;
         }
      }

      // Evenly spaced permanent nodes (about spacingM apart) between the anchors and the pinned corners. A ball then rolls on the polyline they form, which already
      // has the dip, instead of loading one material point that holds it where it landed. Called once, after PinStops.
      void FillNodes(double spacingM)
      {
         std::vector<double> cuts(1, 0.0);
         for (int i = 0; i < Count(); i++)
            cuts.push_back(nodes[order[i]].f);
         cuts.push_back(1.0);
         for (size_t k = 0; k + 1 < cuts.size(); k++)
         {
            const double len = (cuts[k + 1] - cuts[k]) * Ls;
            const int n = std::max(1, (int)std::lround(len / spacingM));
            for (int j = 1; j < n && Count() < 60; j++)
               Add(cuts[k] + (cuts[k + 1] - cuts[k]) * j / n);
         }
      }

      // Residual position fix after a step: a node that has ended up inside a stop (the starting state, or a hit by a force the velocity constraint in StepOnce
      // could not absorb) is put on its surface and loses its inward velocity. Normally there is nothing to do. No node is ever added here: a wrap node moved out
      // to the post in one jump would stretch its legs and put energy into the cord. Returns false (kept for the caller, which rebuilds the ring when it is true).
      bool ApplyStops()
      {
         for (int i = 0; i < Count() && !stops.empty(); i++)
         {
            Node &n = nodes[order[i]];
            for (const StopCircle &s : stops)
            {
               const double dx = n.P.x - s.x, dy = n.P.y - s.y, d = std::hypot(dx, dy);
               constexpr double kDeep = 0.003;   // m: only a node deeper than this is moved at once (the velocity constraint in StepOnce pushes out the rest)
               if (d >= s.r - kDeep || d < 1e-9)
                  continue;
               const double nx = dx / d, ny = dy / d;
               n.P.x = s.x + nx * (s.r - kDeep);
               n.P.y = s.y + ny * (s.r - kDeep);
               const double vn = n.vp.x * nx + n.vp.y * ny;
               if (vn < 0.0)
               {
                  n.vp.x -= vn * nx;
                  n.vp.y -= vn * ny;
               }
            }
         }
         return false;
      }

      int Count() const { return (int)order.size(); }
      // points along the cord: 0 = A, 1..Count() = the nodes in order, Count() + 1 = B
      V3 PtPos(int k) const { return k <= 0 ? A : (k > Count() ? B : nodes[order[k - 1]].P); }
      V3 PtVel(int k) const { return (k <= 0 || k > Count()) ? V3() : nodes[order[k - 1]].vp; }
      int PtNode(int k) const { return (k <= 0 || k > Count()) ? -1 : order[k - 1]; }   // node id, -1 = an anchor
      double PtF(int k) const { return k <= 0 ? 0.0 : (k > Count() ? 1.0 : nodes[order[k - 1]].f); }

      void UpdateMasses()
      {
         for (int i = 0; i < Count(); i++)
         {
            const double share = 0.5 * (PtF(i + 2) - PtF(i));   // half the material on each side
            Node &n = nodes[order[i]];
            n.m = std::max(1e-5, nodeMassKg > 0.0 ? nodeMassKg * 2.0 * share : rhoA * Ls * share);
         }
      }

      // Add a node at the material fraction f. Its place is on the current shape of the leg it splits, with the velocity of that leg there. Returns its id.
      int Add(double f)
      {
         f = std::min(std::max(f, 0.02), 0.98);
         int pos = 0;   // the new node goes before order[pos]
         while (pos < Count() && nodes[order[pos]].f < f)
            pos++;
         const double f0 = PtF(pos), f1 = PtF(pos + 1);
         const double u = f1 > f0 ? std::min(std::max((f - f0) / (f1 - f0), 0.0), 1.0) : 0.5;
         Node n;
         n.f = f;
         n.P0 = A + (B - A) * f;
         n.P = PtPos(pos) + (PtPos(pos + 1) - PtPos(pos)) * u;
         n.vp = PtVel(pos) + (PtVel(pos + 1) - PtVel(pos)) * u;
         int id = -1;
         for (size_t k = 0; k < nodes.size() && id < 0; k++)
            if (!nodes[k].alive)
               id = (int)k;
         mBuf.resize(nodes.size() + 1);
         for (size_t k = 0; k < nodes.size(); k++)
            mBuf[k] = nodes[k].m;
         if (id < 0)
         {
            nodes.push_back(n);
            id = (int)nodes.size() - 1;
         }
         else
            nodes[id] = n;
         order.insert(order.begin() + pos, id);
         UpdateMasses();
         // Re-lumping the cord's mass must not make energy: the nodes that gave mass keep their speed, and the new node moves with the momentum of the mass it took
         // (never faster than that mass did). Measured 2026-10-29 (runs/2026-10-29-band): re-lumping at the old speeds made energy at merges.
         double S = 0.0;
         V3 p;
         for (const int k : order)
            if (k != id && nodes[k].m < mBuf[k])
            {
               S += mBuf[k] - nodes[k].m;
               p = p + nodes[k].vp * (mBuf[k] - nodes[k].m);
            }
         if (S > 0.0)
            nodes[id].vp = p / std::max(nodes[id].m, S);
         return id;
      }

      // Elastic energy of leg j (a leg cannot push: nothing stored when it is slack)
      double LegEnergy(int j) const
      {
         const double l0 = std::max((PtF(j + 1) - PtF(j)) * L0, 1e-9);
         const double e = (PtPos(j + 1) - PtPos(j)).Length() / l0 - 1.0;
         return e > 0.0 ? 0.5 * EA * e * e * l0 : 0.0;
      }
      double Energy() const
      {
         double E = 0.0;
         for (int j = 0; j <= Count(); j++)
            E += LegEnergy(j);
         return E;
      }
      double LegTension(int j) const   // elastic part only, from the geometry
      {
         const double l0 = std::max((PtF(j + 1) - PtF(j)) * L0, 1e-9);
         return std::max(0.0, EA * ((PtPos(j + 1) - PtPos(j)).Length() / l0 - 1.0));
      }

      // Merge a node away: its two legs become one, the unstretched lengths add (the same sticking rule as Add). A merge never gains energy (the chord is not
      // longer than the two legs). Never the last node of a chain.
      void Remove(int id)
      {
         if (!nodes[id].alive || Count() <= 1)
            return;
         mBuf.resize(nodes.size());
         for (size_t k = 0; k < nodes.size(); k++)
            mBuf[k] = nodes[k].m;
         const double mr = nodes[id].m;
         const V3 vr = nodes[id].vp;
         for (size_t k = 0; k < order.size(); k++)
            if (order[k] == id)
            {
               order.erase(order.begin() + (long)k);
               break;
            }
         nodes[id].alive = false;
         nodes[id].vp = V3();
         UpdateMasses();
         // Re-lumping without making energy: the neighbours take the merged node's mass as a plastic collision (its momentum shared by the mass each takes, never
         // more than it had), so they slow down to keep momentum; a node that lost mass keeps its speed. (They used to keep their speed at the larger mass: measured
         // 2026-10-29, runs/2026-10-29-band, +19.7 mJ from merges in one 30 deg hit with no losses.)
         double S = 0.0;
         for (const int k : order)
            if (nodes[k].m > mBuf[k])
               S += nodes[k].m - mBuf[k];
         if (S > 0.0)
            for (const int k : order)
               if (nodes[k].m > mBuf[k])
               {
                  const double dm = nodes[k].m - mBuf[k];
                  nodes[k].vp = (nodes[k].vp * mBuf[k] + vr * (dm / std::max(S, mr) * mr)) / nodes[k].m;
               }
      }

      // Idle nodes (loaded[id] == 0, not within guardR of a guard point) that lie on the straight line between their neighbours are merged away, one per call.
      bool MergeIdle(const std::vector<char> &loaded, double devTol, double speedTol, const std::vector<V3> &guards, double guardR)
      {
         for (int i = 0; i < Count() && Count() > 1; i++)
         {
            const int id = order[i];
            if (id < (int)loaded.size() && loaded[id])
               continue;
            const Node &n = nodes[id];
            if (n.pinned)
               continue;   // a corner on a post stays
            bool guarded = false;   // a companion: another span is being pushed here by a ball, this node lets it bend where it is pressed
            for (const V3 &g : guards)
               guarded = guarded || (n.P - g).Length() < guardR;
            if (guarded)
               continue;
            const V3 p0 = PtPos(i), d = PtPos(i + 2) - p0;
            const double dl = d.Length();
            const double dev = dl > 1e-9 ? (n.P - p0).Cross(d).Length() / dl : (n.P - p0).Length();   // perpendicular distance from the line through its neighbours
            if (dev < devTol && n.vp.Length() < speedTol)
            {
               Remove(id);
               return true;
            }
         }
         return false;
      }

      // Advance every node one step under fContact (N, indexed by node id; the contact reactions, already shared out to the legs' ends)
      // zetaOverride >= 0 replaces the nodes' damping ratio for this step (the band's ring-down after the last ball has left)
      // One step of dt. The explicit integration is stable only while (damping coefficient of a node) * dt / mass < 2 and (stiffness / mass) * dt^2 < 4. A leg of
      // length l has stiffness EA / l and damping eta / l, a node's mass goes with l, so both limits get worse as 1 / l^2: nodes closer than about 8 mm (a ball's
      // node next to a plate node, a node near the end of a span) blew up at the 20 us step (the plate's angle ran away to infinity). Such a step is split into as
      // many equal substeps as the shortest leg needs, with the contact forces held; with legs of 8 mm and more it is one step, as before.
      void Step(const std::vector<V3> &fContact, double dt, double zetaOverride = -1.0)
      {
         const int M = Count();
         if (M == 0)
            return;
         double worst = 0.0;
         for (int i = 0; i < M; i++)
         {
            const double lL = std::max((PtF(i + 1) - PtF(i)) * L0, 1e-9), lR = std::max((PtF(i + 2) - PtF(i + 1)) * L0, 1e-9);
            const double m = nodes[order[i]].m;
            const double damp = eta * (1.0 / lL + 1.0 / lR) * dt / m / 1.5;           // damping criterion, 1 = at the limit used here
            const double stiff = std::sqrt(EA * (1.0 / lL + 1.0 / lR) / m) * dt / 1.0;   // omega * dt, 1 = at the limit used here
            worst = std::max(worst, std::max(damp, stiff));
         }
         if (worst > 1.0)
         {
            const int sub = std::min(64, (int)std::ceil(worst));
            const double h = dt / (double)sub;
            for (int k = 0; k < sub; k++)
               StepOnce(fContact, h, zetaOverride);
            return;
         }
         StepOnce(fContact, dt, zetaOverride);
      }

      void StepOnce(const std::vector<V3> &fContact, double dt, double zetaOverride = -1.0)
      {
         StepVelocities(fContact, dt, zetaOverride);
         StepPositions(dt);
      }

      // StepOnce in its two halves, so a constraint between chains can act on the new velocities before the nodes move (BandPart: two spans of a loop meeting)
      void StepVelocities(const std::vector<V3> &fContact, double dt, double zetaOverride = -1.0)
      {
         const int M = Count();
         if (M == 0)
            return;
         const double z = zetaOverride >= 0.0 ? zetaOverride : zeta;
         T.resize(M + 1);
         uBuf.resize(M + 1);
         lBuf.resize(M + 1);
         std::vector<V3> &u = uBuf;
         std::vector<double> &L = lBuf;
         for (int j = 0; j <= M; j++)
         {
            const V3 d = PtPos(j + 1) - PtPos(j);
            L[j] = std::max(d.Length(), 1e-9);
            u[j] = d / L[j];
            const double l0 = std::max((PtF(j + 1) - PtF(j)) * L0, 1e-9);
            const double e = L[j] / l0 - 1.0;
            const double rate = (PtVel(j + 1) - PtVel(j)).Dot(u[j]) / l0;   // strain rate, 1/s
            T[j] = std::max(0.0, EA * e + eta * rate);
         }
         for (int i = 0; i < M; i++)
         {
            Node &n = nodes[order[i]];
            V3 F = fContact[order[i]] + u[i + 1] * T[i + 1] - u[i] * T[i];
            if (n.pinned)
            {
               // fixed to the post until the forces on it point away from it
               const StopCircle &s = stops[std::min<size_t>((size_t)std::max(n.pinStop, 0), stops.size() - 1)];
               if (F.x * (n.P.x - s.x) + F.y * (n.P.y - s.y) > 0.0)
                  n.pinned = false;
               else
               {
                  n.vp = V3();
                  continue;
               }
            }
            const double kT = T[i] / L[i] + T[i + 1] / L[i + 1];
            F = F - n.vp * (2.0 * z * std::sqrt(kT * n.m));   // the node's own vibration dies out (critical damping of its two legs times zeta)
            if (lockZ)
               F.z = 0.0;   // the vertical part is carried by the posts
            n.vp = n.vp + F * (dt / n.m);
            if (lockZ)
               n.vp.z = 0.0;
            // a stop is a velocity constraint (applied before the position update, so the node never has to be teleported out of the post, which would
            // stretch its legs and put energy into the cord): the node may approach until it touches, and what is left of the inward velocity is dropped
            // The same for the two legs at the node: the leg's point nearest to the stop moves with (1 - t) of its far end's and t of this node's velocity,
            // so a leg cannot be pulled through a post either (no wrap node is inserted: that would jump a node out to the post and stretch its legs).
            for (const StopCircle &s : stops)
               for (int side = 0; side < 3; side++)
               {
                  const V3 m = side == 0 ? n.P : (side == 1 ? PtPos(i) : PtPos(i + 2));   // the far end of the leg (side 0: the node itself)
                  const V3 vm = side == 0 ? V3() : (side == 1 ? PtVel(i) : PtVel(i + 2));
                  double t = 1.0;
                  V3 q = n.P;
                  if (side > 0)
                  {
                     const double ex = n.P.x - m.x, ey = n.P.y - m.y, l2 = ex * ex + ey * ey;
                     if (l2 < 1e-12)
                        continue;
                     t = std::min(std::max(((s.x - m.x) * ex + (s.y - m.y) * ey) / l2, 0.0), 1.0);
                     if (t < 0.15)
                        continue;   // the nearest point is by the far end, which this node hardly moves
                     q = m + (n.P - m) * t;
                  }
                  const double dx = q.x - s.x, dy = q.y - s.y, d = std::hypot(dx, dy);
                  if (d > s.r + 0.02 || d < 1e-9)   // 20 mm: beyond one step's reach
                     continue;
                  const double nx = dx / d, ny = dy / d;
                  const double vq = (vm.x * (1.0 - t) + n.vp.x * t) * nx + (vm.y * (1.0 - t) + n.vp.y * t) * ny;   // side 0: t = 1, vm unused (weight 0)
                  // outside: may approach until it touches. Inside (a post that overlaps the band at rest): is pushed out at a bounded speed, so the band
                  // moves away from the post smoothly instead of jumping (a jump stretches the legs at once and kicks the ball on it)
                  const double vMin = d >= s.r ? -(d - s.r) / dt : std::min((s.r - d) / dt, 0.3);
                  if (vq < vMin)
                  {
                     const double dv = (vMin - vq) / t;
                     n.vp.x += dv * nx;
                     n.vp.y += dv * ny;
                  }
                  // the band sticks to the post like a corner post (rubber on a post does not slide): a node that is on the stop and pushed into it stays where
                  // it is, so the two legs run from the post to the neighbours; it is let go when the forces on it point away from the post
                  if (side == 0 && std::fabs(d - s.r) <= 1e-5 && F.x * nx + F.y * ny < 0.0)   // (a node still inside is on its way out, not stuck)
                  {
                     n.vp.x = 0.0;
                     n.vp.y = 0.0;
                  }
               }
         }
      }
      void StepPositions(double dt)
      {
         const int M = Count();
         for (int i = 0; i < M; i++)
         {
            Node &n = nodes[order[i]];
            n.P = n.P + n.vp * dt;
            n.peak = std::max(n.peak, n.Deflection());
         }
      }

      // posTol: the largest deflection that still counts as at rest, m; velTol: the largest node speed, m/s
      bool AtRest(double posTol = 5e-5, double velTol = 5e-2) const
      {
         for (const Node &n : nodes)
            if (n.alive && !n.pinned && (n.Deflection() >= posTol || n.vp.Length() >= velTol))
               return false;
         return true;
      }
   };

   // Closest points of the segments p1-q1 and p2-q2: parameters s, t in [0,1] along each, returns the distance (Ericson, Real-Time Collision Detection)
   inline double SegSegClosest(const V3 &p1, const V3 &q1, const V3 &p2, const V3 &q2, double &s, double &t)
   {
      const V3 d1 = q1 - p1, d2 = q2 - p2, r = p1 - p2;
      const double a = d1.Dot(d1), e = d2.Dot(d2), f = d2.Dot(r);
      const double eps = 1e-18;
      if (a <= eps && e <= eps)
         s = t = 0.0;
      else if (a <= eps)
      {
         s = 0.0;
         t = std::min(std::max(f / e, 0.0), 1.0);
      }
      else
      {
         const double c = d1.Dot(r);
         if (e <= eps)
         {
            t = 0.0;
            s = std::min(std::max(-c / a, 0.0), 1.0);
         }
         else
         {
            const double b = d1.Dot(d2), denom = a * e - b * b;
            s = denom > eps ? std::min(std::max((b * f - c * e) / denom, 0.0), 1.0) : 0.0;
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
      return ((p1 + d1 * s) - (p2 + d2 * t)).Length();
   }

   // Band against band: two spans of the same loop pushed together (a ball presses the front span of a loop onto the back one, or two balls press it at
   // two places). The legs of the two chains are tubes of radius rho; they cannot overlap. Hertz between two cylinders touching at a point (sphere-like,
   // Reff = rho/2) with the same Hunt-Crossley damping as the ball contact, no friction. Each contact (the nearest pair of legs of a cluster: a leg pair with
   // a neighbour pair that is closer is skipped, so a node between two legs does not count twice) is shared out to the ends of the two legs.
   // fX / fY: force on the nodes of X / Y, indexed by node id. Returns whether any contact was found; depth = the deepest overlap.
   inline bool BandChainContact(const RubberParams &rp, double rho, const BandChain &X, const BandChain &Y, std::vector<V3> &fX, std::vector<V3> &fY, double *depth = nullptr)
   {
      constexpr int kMaxLegs = 16;
      const int mx = std::min(X.Count() + 1, kMaxLegs), my = std::min(Y.Count() + 1, kMaxLegs);
      {
         // boxes round the two chains' points: further apart than 2 rho in any axis means no leg of one can touch a leg of the other
         V3 xl = X.PtPos(0), xh = xl, yl = Y.PtPos(0), yh = yl;
         for (int k = 1; k <= mx; k++)
         {
            const V3 q = X.PtPos(k);
            xl = V3(std::min(xl.x, q.x), std::min(xl.y, q.y), std::min(xl.z, q.z));
            xh = V3(std::max(xh.x, q.x), std::max(xh.y, q.y), std::max(xh.z, q.z));
         }
         for (int k = 1; k <= my; k++)
         {
            const V3 q = Y.PtPos(k);
            yl = V3(std::min(yl.x, q.x), std::min(yl.y, q.y), std::min(yl.z, q.z));
            yh = V3(std::max(yh.x, q.x), std::max(yh.y, q.y), std::max(yh.z, q.z));
         }
         const double g = 2.0 * rho;
         if (xl.x > yh.x + g || yl.x > xh.x + g || xl.y > yh.y + g || yl.y > xh.y + g || xl.z > yh.z + g || yl.z > xh.z + g)
            return false;
      }
      double D[kMaxLegs][kMaxLegs], S[kMaxLegs][kMaxLegs], Tt[kMaxLegs][kMaxLegs];
      bool any = false;
      for (int j = 0; j < mx; j++)
         for (int k = 0; k < my; k++)
            D[j][k] = SegSegClosest(X.PtPos(j), X.PtPos(j + 1), Y.PtPos(k), Y.PtPos(k + 1), S[j][k], Tt[j][k]);
      const double Estar = rp.eMpa * 1e6 / (1.0 - rp.nu * rp.nu);
      for (int j = 0; j < mx; j++)
         for (int k = 0; k < my; k++)
         {
            const double d = D[j][k];
            if (!(d < 2.0 * rho))
               continue;
            if ((j > 0 && D[j - 1][k] < d) || (j + 1 < mx && D[j + 1][k] < d) || (k > 0 && D[j][k - 1] < d) || (k + 1 < my && D[j][k + 1] < d))
               continue;
            const double s = S[j][k], t = Tt[j][k];
            const V3 xa = X.PtPos(j), xb = X.PtPos(j + 1), ya = Y.PtPos(k), yb = Y.PtPos(k + 1);
            const V3 qx = xa + (xb - xa) * s, qy = ya + (yb - ya) * t;
            V3 n = qx - qy;
            n = d > 1e-9 ? n / d : V3(0.0, 0.0, 1.0);
            const double delta = 2.0 * rho - d;
            if (depth)
               *depth = std::max(*depth, delta);
            const V3 vx = X.PtVel(j) * (1.0 - s) + X.PtVel(j + 1) * s, vy = Y.PtVel(k) * (1.0 - t) + Y.PtVel(k + 1) * t;
            const double vn = (vx - vy).Dot(n);
            const double Fe = (4.0 / 3.0) * Estar * std::sqrt(0.5 * rho) * delta * std::sqrt(delta);
            const double Fn = std::max(0.0, Fe * (1.0 + 1.5 * rp.alphaSPerM * (-vn)));
            const V3 F = n * Fn;   // on X at the contact; the reaction on Y
            const int nx0 = X.PtNode(j), nx1 = X.PtNode(j + 1), ny0 = Y.PtNode(k), ny1 = Y.PtNode(k + 1);
            if (nx0 >= 0) fX[nx0] = fX[nx0] + F * (1.0 - s);
            if (nx1 >= 0) fX[nx1] = fX[nx1] + F * s;
            if (ny0 >= 0) fY[ny0] = fY[ny0] - F * (1.0 - t);
            if (ny1 >= 0) fY[ny1] = fY[ny1] - F * t;
            any = true;
         }
      return any;
   }
}
