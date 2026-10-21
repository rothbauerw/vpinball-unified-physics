// license:GPLv3+
//
// A part that turns about one fixed horizontal axis (a gate, a spinner): the shared body of the rebuilt gate and spinner models (docs/GateSpinnerDesign.txt, section 3).
//
// Pure SI maths with no VPX types, so the same code runs inside VPX and in the offline tests (runs/2026-10-14-hinge). The part supplies its shape (wire segments and
// flat plates, from which the mass, the centre of mass and the inertia about the axis follow), its stops and its pin; the owner feeds the gravity it feels (gravity
// minus the cabinet acceleration, so a nudge swings it) and the torque and force of the balls touching it, and reads the angle back.
//
// Frame: the axis passes through A along E (unit, horizontal). N is the direction a ball moves when it passes from the front, Z = E x N is up, E = N x Z.
// theta = 0: the body as built (a gate hangs straight down); a positive theta swings the lower part toward N (a rotation about +E). A point of the body is given in
// body coordinates (e along E, d "down" at theta = 0, f "front" at theta = 0) and sits at A + E e + Down(theta) d + Front(theta) f,
// Down(theta) = -Z cos theta + N sin theta, Front(theta) = N cos theta + Z sin theta.
//
// Torques (all summed before the body moves):
//   gravity     E . (r_cm x m g)                                   g = the gravity vector the body feels (SI, world frame)
//   balls       the reactions the contact system adds (AddExternal)
//   viscous     -c omega                                           c >= 0 (a spinner's spin-down; 0 for a gate)
//   stops       Hertz + viscoelastic damping at each stop's lever arm (section 4.5 of the design): Fn = k delta^1.5 + gamma a delta_dot, a = sqrt(Reff delta),
//               k = 4/3 E* sqrt(Reff) (steel on steel); gamma from the stop's restitution at a reference speed (CalibrateKK). No angle clamp, no bounce factor.
//   pin         Coulomb friction of the pin wire turning in its hole: |tau| <= mu r_pin |N|, N = the load the pins carry (weight, balls, stops, the centripetal
//               force of the centre of mass). The body sticks when the other torques cannot overcome the limit, and slides with the full limit otherwise.
// Integration: semi-implicit Euler (omega, then theta). A body that is stuck, untouched and not driven sleeps: Advance does nothing until something changes.

#pragma once

#include <cmath>
#include <algorithm>
#include <vector>

#include "physics/rubbercontact.h"

namespace HingeBody
{
   using RubberContact::V3;

   // Mass properties of the body about its axis, built from its shape (body coordinates, metres: x = e along the axis, y = d down, z = f front)
   struct MassProps
   {
      double m = 0.0;     // kg
      V3 mc;              // sum of m * position (first moment), kg m
      double I = 0.0;     // about the axis, kg m^2
      V3 Cm() const { return m > 0.0 ? mc / m : V3(); }

      // A straight wire from p to q, linear density lambda (kg/m). Exact for a thin rod: the squared distance to the axis is quadratic along it.
      void AddWire(const V3 &p, const V3 &q, double lambda)
      {
         const double L = (q - p).Length();
         const double dm = lambda * L;
         m += dm;
         mc = mc + (p + q) * (0.5 * dm);
         // r_perp^2 = d^2 + f^2: integral over the segment of (a + t (b - a))^2 = L (a^2 + a b + b^2) / 3 per component
         auto sq = [](double a, double b) { return (a * a + a * b + b * b) / 3.0; };
         I += dm * (sq(p.y, q.y) + sq(p.z, q.z));
      }
      // A flat triangle (a, b, c), areal density sigma (kg/m^2). Exact for a thin plate (second moment of a triangle).
      void AddPlate(const V3 &a, const V3 &b, const V3 &c, double sigma)
      {
         const double A = 0.5 * (b - a).Cross(c - a).Length();
         const double dm = sigma * A;
         m += dm;
         mc = mc + (a + b + c) * (dm / 3.0);
         auto sq = [](double x, double y, double z) { return (x * x + y * y + z * z + x * y + y * z + z * x) / 6.0; };
         I += dm * (sq(a.y, b.y, c.y) + sq(a.z, b.z, c.z));   // integral of y^2 over a triangle = A/6 (sum of squares + cross terms of the vertices)
      }
      void AddQuad(const V3 &a, const V3 &b, const V3 &c, const V3 &d, double sigma)
      {
         AddPlate(a, b, c, sigma);
         AddPlate(a, c, d, sigma);
      }
   };

   // Steel on steel (E 200 GPa, nu 0.3, both bodies): E* = E / (2 (1 - nu^2))
   constexpr double kSteelEstar = 200e9 / (2.0 * (1.0 - 0.3 * 0.3));

   // A stop: a piece of the bracket the body touches at a lever arm from the axis. side = +1: the stop is reached by turning to larger theta (it blocks theta > angle);
   // side = -1: it blocks theta < angle. A driven stop (the hold actuator of a gate) moves its angle toward target at rate (rad/s); its motion is part of the
   // contact's approach speed.
   struct Stop
   {
      double angle = 0.0;      // rad
      int side = 1;
      double lever = 0.005;    // m
      double reff = 5e-4;      // Hertz effective radius of the touch, m (crossed cylinders: sqrt(r1 r2))
      double k = 0.0;          // N / m^1.5
      double gamma = 0.0;      // N s / m^2
      bool enabled = true;
      bool armed = true;       // a stop engages only from its free side: one placed while the body is past it waits until the body is on its free side
      bool touching = false;
      // driven stop
      double target = 0.0, rate = 0.0;   // rate 0 = fixed
      double speed = 0.0;                // the stop's own angular speed this substep (rad/s)
      double Delta(double theta) const { return side * (theta - angle) * lever; }
   };

   // k and gamma of a stop: steel on steel at Reff, restitution e at the reference impact speed vRef (m/s, at the lever arm), for the body's effective mass
   // at the lever arm (I / lever^2). The viscoelastic law makes the restitution fall slowly with speed (as real metal impacts do).
   inline void SetupStop(Stop &s, double inertia, double eStop, double vRef)
   {
      s.k = (4.0 / 3.0) * kSteelEstar * std::sqrt(s.reff);
      RubberContact::RubberParams rp;
      rp.nu = 0.0;                      // ImpactRestitutionKK uses E / (1 - nu^2) as E*: pass E* itself
      rp.eMpa = kSteelEstar * 1e-6;
      rp.normalLoss = false;
      const double mEff = inertia / (s.lever * s.lever);
      s.gamma = RubberContact::CalibrateKK(rp, s.reff, mEff, vRef, eStop);
   }

   // A one-sided torsion spring (a switch blade the body presses once it is past angle): torque -side k (theta - angle) while side (theta - angle) > 0.
   // Conservative: it takes energy while it is pressed and gives it back as it is released.
   struct Spring
   {
      double angle = 0.0;   // rad
      int side = 1;
      double k = 0.0;       // N m / rad
      bool enabled = true;
   };

   // A crank on the shaft pressing a switch blade (a spinner's leaf or micro switch): the crank of radius r points along the angle phase; the blade sits gap from
   // the shaft, so it is pushed in by delta = r cos(theta - phase) - gap while that is positive. The blade pushes back with F = F0 + k delta + c delta_dot (F0: a
   // micro switch's preload; 0 for a leaf), never pulling; the torque on the body is -F d(delta)/d(theta) = F r sin(theta - phase). The spring part is
   // conservative (it takes energy while the blade is pressed in and gives it back as it is released: the speed-up / slow-down of a spinner over its
   // switch); c is the blade's small loss.
   struct Cam
   {
      double phase = 0.0;   // rad
      double r = 3e-3;      // crank radius, m
      double gap = 0.5e-3;  // blade distance from the shaft, m
      double k = 0.0;       // N/m
      double F0 = 0.0;      // N (preload while pressed)
      double c = 0.0;       // N s/m
      bool enabled = true;
      double Delta(double th) const { return r * std::cos(th - phase) - gap; }
   };

   struct StopEvent
   {
      int stop;           // index into stops
      double speedDegS;   // |omega| at the touch, deg/s (the stock LimitEOS/LimitBOS parameter)
   };

   class Body
   {
   public:
      // ---- shape and settings ----
      V3 A, E, N, Z;              // axis point and frame (world, m)
      MassProps mp;               // from the shape
      double rPin = 0.6e-3;       // pin (wire) radius, m
      double muPivot = 0.3;       // pin in its hole (steel on steel)
      double c = 0.0;             // viscous loss, N m s
      std::vector<Stop> stops;
      std::vector<Spring> springs;
      std::vector<Cam> cams;

      // ---- state ----
      double theta = 0.0, omega = 0.0;
      bool stuck = false;         // held by the pin friction in the last substep
      bool asleep = false;
      std::vector<StopEvent> events;   // stop touches since the owner last cleared it

      void SetFrame(const V3 &axisPoint, const V3 &front)
      {
         A = axisPoint;
         N = front / front.Length();
         Z = V3(0.0, 0.0, 1.0);
         E = N.Cross(Z);
         E = E / E.Length();
         N = Z.Cross(E);   // horizontal, perpendicular to E (E = N x Z)
      }
      V3 Down(double th) const { return Z * -std::cos(th) + N * std::sin(th); }
      V3 Front(double th) const { return N * std::cos(th) + Z * std::sin(th); }
      // body coordinates (x = e, y = d, z = f) -> world
      V3 ToWorld(const V3 &b, double th) const { return A + E * b.x + Down(th) * b.y + Front(th) * b.z; }
      V3 ToWorld(const V3 &b) const { return ToWorld(b, theta); }
      // velocity of a body point (world, relative to A) and the angular velocity vector
      V3 VelocityAt(const V3 &worldPoint) const { return (E * omega).Cross(worldPoint - A); }

      // ---- per-substep input from the contact system ----
      void AddExternal(double torqueAboutE, const V3 &force)
      {
         m_extTau += torqueAboutE;
         m_extF = m_extF + force;
         m_hadExternal = true;
      }
      // takes the external torque and force accumulated so far (the balls' reactions of an island substep) and clears them
      void TakeExternal(double &tau, V3 &F)
      {
         tau = m_extTau;
         F = m_extF;
         m_extTau = 0.0;
         m_extF = V3();
      }
      // the gravity the body feels (world, m/s^2): gravity minus the cabinet acceleration
      void SetGravity(const V3 &g)
      {
         if (g.x != m_g.x || g.y != m_g.y || g.z != m_g.z)
            asleep = false;   // a nudge or a slope change: the body has a new equilibrium
         m_g = g;
      }
      const V3 &Gravity() const { return m_g; }

      // adds a stop; it engages only from its free side, so one placed while the body is past it waits until the body comes back (a one-way stop switched on
      // while a two-way gate is swung backward). The body at exactly the stop's angle (a one-way gate at its rest angle) is on the free side.
      int PlaceStop(Stop s)
      {
         s.touching = false;
         s.armed = s.Delta(theta) <= 0.0;
         stops.push_back(s);
         asleep = false;
         return (int)stops.size() - 1;
      }

      // the script moves a driven stop (a gate's hold actuator) toward a new angle
      void DriveStop(int i, double target, double rate)
      {
         stops[i].target = target;
         stops[i].rate = std::abs(rate);
         asleep = false;
      }
      void Wake() { asleep = false; }

      // ---- one substep of h seconds (the contact system calls it inside an island; Advance calls it for a body alone) ----
      void Step(double h)
      {
         SetDrivenStopSpeeds();   // a driven stop's speed enters the contact damping now; its angle moves with the body's at the end of the substep
         V3 pinLoad;   // force the body gets from everything but the pins (the pins carry the opposite, plus the centripetal force)
         double tauOther = m_extTau + GravityTorque(theta) - c * omega + StopTorques(theta, omega, pinLoad) + SpringTorques(theta) + CamTorques(theta, omega);
         pinLoad = pinLoad + m_g * mp.m + m_extF;
         // centripetal: the pins also hold the centre of mass on its circle (m a = F_pins + F_other, a = -omega^2 r_cm perp)
         const V3 rcm = ToWorld(mp.Cm()) - A;
         const V3 rcmPerp = rcm - E * rcm.Dot(E);
         const V3 pinForce = rcmPerp * (-mp.m * omega * omega) - pinLoad;
         const V3 pinPerp = pinForce - E * pinForce.Dot(E);
         const double tauLim = muPivot * rPin * pinPerp.Length();
         const double I = mp.I;

         double wNew;
         if (omega != 0.0)
         {
            const double sgn = omega > 0.0 ? 1.0 : -1.0;
            const double alpha = (tauOther - sgn * tauLim) / I;
            wNew = omega + alpha * h;
            if (wNew * omega < 0.0)   // the body stops inside the substep: Coulomb decides whether it stays or turns back
            {
               if (std::abs(tauOther) <= tauLim)
                  wNew = 0.0;
               else
               {
                  const double tStop = -omega / alpha;
                  const double sgn2 = tauOther > 0.0 ? 1.0 : -1.0;
                  wNew = (tauOther - sgn2 * tauLim) / I * (h - tStop);
               }
            }
         }
         else if (std::abs(tauOther) <= tauLim)
            wNew = 0.0;
         else
         {
            const double sgn = tauOther > 0.0 ? 1.0 : -1.0;
            wNew = (tauOther - sgn * tauLim) / I * h;
         }
         stuck = wNew == 0.0;
         omega = wNew;
         theta += omega * h;
         MoveDrivenStops(h);
         UpdateStopStates();
         m_extTau = 0.0;
         m_extF = V3();
      }

      // ---- a whole interval (a physics cycle) for a body no ball touches: adaptive substeps ----
      // holdTorque / holdForce act during the whole interval (added before every substep), as the bat's extTorque. Returns the number of substeps (0 while asleep).
      int Advance(double dt, double holdTorque = 0.0, const V3 &holdForce = V3())
      {
         return AdvanceWith(dt, [] {}, holdTorque, holdForce);
      }
      // the same, calling afterSubstep() after every substep (a part's switches must see every substep: a spinner turns 45 deg in 1 ms at 800 rad/s)
      template <class F> int AdvanceWith(double dt, F afterSubstep, double holdTorque = 0.0, const V3 &holdForce = V3())
      {
         const bool held = holdTorque != 0.0 || holdForce.Dot(holdForce) != 0.0;
         if (held)
            m_hadExternal = true;
         if (asleep && !m_hadExternal)
            return 0;
         asleep = false;
         // equal substeps over the interval (semi-implicit Euler keeps the energy only with a steady step; unequal pieces drift). The step is re-chosen
         // after each one in case a stop comes into reach, and the rest of the interval is again split evenly.
         int n = 0;
         double t = 0.0;
         while (t < dt)
         {
            const double rest = dt - t;
            if (held)   // the held torque is part of the acceleration the step must resolve
            {
               m_extTau += holdTorque;
               m_extF = m_extF + holdForce;
            }
            const int k = std::max(1, (int)std::ceil(rest / ChooseStep(rest) - 1e-9));
            if (held)
            {
               m_extTau -= holdTorque;
               m_extF = m_extF - holdForce;
            }
            const double h = rest / k;
            if (held)
               AddExternal(holdTorque, holdForce);
            Step(h);
            afterSubstep();
            t = (k == 1) ? dt : t + h;
            n++;
         }
         // asleep: held by the pin (or at rest on a stop with the pin holding), nothing driving it
         asleep = stuck && !m_hadExternal && !AnyStopMoving();
         m_hadExternal = false;
         return n;
      }

      // The largest substep that resolves the body's own forces (the island takes the smaller of this and the ball contacts' step)
      double ChooseStep(double hMax) const
      {
         double h = hMax;
         const double I = mp.I;
         // the angle step comes from the highest speed the body's energy allows (all of it kinetic at the bottom of the swing), not from the present speed: a
         // step that follows the state through every swing breaks the energy behaviour of semi-implicit Euler (measured: a frictionless 60 deg swing lost half
         // its energy in 8 s that way); this bound only changes when the energy does
         const double eKinMax = Energy() + mp.m * std::sqrt(m_g.Dot(m_g)) * PendulumArm();
         const double wMax = std::max(std::abs(omega), std::sqrt(std::max(0.0, 2.0 * eKinMax / I)));
         if (wMax > 0.0)
            h = std::min(h, kMaxAngleStep / wMax);
         if (c > 0.0)
            h = std::min(h, 0.2 * I / c);
         for (const Spring &sp : springs)
            if (sp.enabled && sp.k > 0.0)
               h = std::min(h, 0.1 / std::sqrt(sp.k / I));
         for (const Cam &cm : cams)
            if (cm.enabled && cm.k > 0.0)
               h = std::min(h, 0.1 / std::sqrt(cm.k * cm.r * cm.r / I));   // the blade as a torsion spring at the crank (its stiffest)
         const double d = PendulumArm();
         if (d > 0.0)
            h = std::min(h, 0.05 / std::sqrt(mp.m * std::sqrt(m_g.Dot(m_g)) * d / I));
         // the body's angular acceleration without the stops, bounded from above: gravity, the balls, the viscous loss, plus the most the pin friction can add
         // (it opposes the motion, so it can speed up the closing on a stop the body moves away from)
         const double alphaFree = std::abs(m_extTau + GravityTorque(theta) - c * omega + SpringTorques(theta) + CamTorques(theta, omega)) / I
            + muPivot * rPin * (mp.m * std::sqrt(m_g.Dot(m_g)) + m_extF.Length()) / I;
         for (const Stop &s : stops)
         {
            if (!s.enabled || !s.armed)
               continue;
            const double mEff = I / (s.lever * s.lever);
            const double delta = s.Delta(theta);
            // a driven stop moves at its rate during the coming substep (its speed from the last one may still be 0)
            const double sp = DrivenSpeed(s);
            const double vc = s.side * (omega - sp) * s.lever;   // closing speed (> 0: closing)
            const double ac = alphaFree * s.lever;               // bound on the closing acceleration
            auto hContact = [&](double vImpact) {
               // the deepest compression this approach can reach (Hertz: 1/2 m v^2 = 2/5 k d^2.5) sets the stiffness to resolve
               const double v = std::max(vImpact, 0.0);
               const double dMax = std::max(delta, 0.0) + std::pow(1.25 * mEff * v * v / s.k, 0.4);
               const double kt = 1.5 * s.k * std::sqrt(std::max(dMax, 1e-12));
               return kContactStepFraction / std::sqrt(kt / mEff);
            };
            if (delta > 0.0)
               h = std::min(h, hContact(vc));
            else
            {
               // earliest time the gap -delta can close. Semi-implicit Euler moves the body by (omega + alpha h) h in a step, so the closing over a step of
               // length t is vc t + ac t^2 (not 1/2 ac t^2); vc may be negative (moving apart, then pulled back)
               const double gap = -delta;
               const double tArr = ac > 0.0 ? (-vc + std::sqrt(vc * vc + 4.0 * ac * gap)) / (2.0 * ac) : (vc > 0.0 ? gap / vc : 1e30);
               if (tArr < h)   // reaches the stop in this step: step up to it, then resolve the touch
                  h = std::max(tArr, hContact(vc + ac * tArr));
            }
         }
         return std::max(h, kMinStep);
      }

      // ---- diagnostics: energy (J) = rotation + gravity potential + stop springs ----
      double Energy() const
      {
         double e = 0.5 * mp.I * omega * omega - mp.m * m_g.Dot(ToWorld(mp.Cm()) - A);
         for (const Stop &s : stops)
         {
            const double d = s.Delta(theta);
            if (s.enabled && s.armed && d > 0.0)
               e += 0.4 * s.k * std::pow(d, 2.5);
         }
         for (const Spring &sp : springs)
         {
            const double x = sp.side * (theta - sp.angle);
            if (sp.enabled && x > 0.0)
               e += 0.5 * sp.k * x * x;
         }
         for (const Cam &cm : cams)
         {
            const double d = cm.Delta(theta);
            if (cm.enabled && d > 0.0)
               e += cm.F0 * d + 0.5 * cm.k * d * d;
         }
         return e;
      }
      // distance of the centre of mass from the axis (m)
      double PendulumArm() const
      {
         const V3 cm = mp.Cm();
         return std::sqrt(cm.y * cm.y + cm.z * cm.z);
      }

      static constexpr double kMaxAngleStep = 0.02;          // rad per substep: the body's shape moves at most ~1 deg (the balls' contact sees it move smoothly)
      static constexpr double kContactStepFraction = 0.1;    // substep = 0.1 / (stop contact frequency)
      static constexpr double kMinStep = 1e-7;               // numerical floor (s), far below any contact time here (a stop touch lasts ~100 us)

   private:
      double m_extTau = 0.0;
      V3 m_extF;
      V3 m_g = V3(0.0, 0.0, -9.81);
      bool m_hadExternal = false;

      double CamTorques(double th, double w) const
      {
         double tau = 0.0;
         for (const Cam &cm : cams)
         {
            const double d = cm.Delta(th);
            if (!cm.enabled || d <= 0.0)
               continue;
            const double s = std::sin(th - cm.phase);
            const double ddot = -cm.r * s * w;
            const double F = std::max(0.0, cm.F0 + cm.k * d + cm.c * ddot);
            tau += F * cm.r * s;
         }
         return tau;
      }
      double SpringTorques(double th) const
      {
         double tau = 0.0;
         for (const Spring &sp : springs)
         {
            const double x = sp.side * (th - sp.angle);
            if (sp.enabled && x > 0.0)
               tau -= sp.side * sp.k * x;
         }
         return tau;
      }
      double GravityTorque(double th) const
      {
         const V3 r = ToWorld(mp.Cm(), th) - A;
         return E.Dot(r.Cross(m_g * mp.m));
      }
      // stop torques at (th, w); adds the stops' forces on the body to load (each acts across its lever arm, along the body's front direction: used only for
      // the pin load, i.e. the friction limit)
      double StopTorques(double th, double w, V3 &load) const
      {
         double tau = 0.0;
         for (const Stop &s : stops)
         {
            if (!s.enabled || !s.armed)
               continue;
            const double delta = s.Delta(th);
            if (delta <= 0.0)
               continue;
            const double ddot = s.side * (w - s.speed) * s.lever;
            const double a = std::sqrt(s.reff * delta);
            const double Fn = std::max(0.0, s.k * delta * std::sqrt(delta) + s.gamma * a * ddot);
            tau -= s.side * s.lever * Fn;
            load = load + Front(th) * (-s.side * Fn);
         }
         return tau;
      }
      static double DrivenSpeed(const Stop &s) { return (s.rate > 0.0 && s.angle != s.target) ? (s.target > s.angle ? s.rate : -s.rate) : 0.0; }
      void SetDrivenStopSpeeds()
      {
         for (Stop &s : stops)
            s.speed = DrivenSpeed(s);
      }
      void MoveDrivenStops(double h)
      {
         for (Stop &s : stops)
         {
            if (s.speed == 0.0)
               continue;
            const double diff = s.target - s.angle;
            const double move = s.speed * h;
            s.angle = std::abs(diff) <= std::abs(move) ? s.target : s.angle + move;
         }
      }
      bool AnyStopMoving() const
      {
         for (const Stop &s : stops)
            if (s.rate > 0.0 && s.angle != s.target)
               return true;
         return false;
      }
      void UpdateStopStates()
      {
         for (int i = 0; i < (int)stops.size(); i++)
         {
            Stop &s = stops[i];
            const double delta = s.Delta(theta);
            if (!s.armed)
            {
               if (delta <= 0.0)
                  s.armed = true;   // the body is on the free side now
               continue;
            }
            const bool t = s.enabled && delta > 0.0;
            if (t && !s.touching)
               events.push_back({ i, std::abs(omega) * 180.0 / M_PI });
            s.touching = t;
         }
      }
   };
}
