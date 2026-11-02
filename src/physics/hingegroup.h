// license:GPLv3+
//
// What every part on one horizontal axis shares as a part of the group step (a gate, a spinner; docs/GateSpinnerDesign.txt sections 3 and 6): its touchable shape (wire capsules and sheet triangles in body coordinates), its stored energy
// (gravity on the centre of mass, a switch blade), the pivot's friction. Its body is a BodyPart (bodypart.h); the moving is the group's. Pure SI, no VPX types.
//
// Frame: the axis through A along E (horizontal). N: the way a ball passes from the front; Z = E x N, up. Body coordinates: x along E, y down at theta = 0, z front at theta = 0; a positive theta swings the lower part toward N.
// Ball contact: an ordinary Item, the steel-on-steel Hertz + viscoelastic law, plus the bend: out of its plane the part bends like a cantilever from the axis, far softer than Hertz, a massless spring in series (Hit::kb), taken when a touch starts and held while it lasts.
// Pivot: Coulomb friction of the pin wire in its hole: |tau| <= mu r_pin |N|, N = the load on the pins (weight, balls, the centre of mass's centripetal force).

#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "physics/bodypart.h"
#include "physics/contacthingebody.h"
#include "physics/rubbercontact.h"

namespace HingeGroup
{
   using GroupStep::Hit;
   using GroupStep::V3;
   using RubberContact::BallSI;
   using RubberContact::RubberParams;
   using RubberContact::TubeProbe;

   class Part : public GroupStep::BodyPart
   {
   public:
      RubberParams rp;            // ball contact law
      V3 A, E, N, Z;              // axis point and frame
      V3 g = V3(0.0, 0.0, -9.81); // the gravity the part feels: gravity minus the cabinet acceleration
      HingeBody::MassProps mp;
      double rPin = 0.6e-3, muPivot = 0.3, cViscous = 0.0;
      std::vector<HingeBody::Cam> cams;
      double kBlade = 0.0, aBlade = 0.0;   // a one-sided torsion spring: a blade pressed once theta is past aBlade
      std::vector<char> swLog;    // what happened, in order, until the owner reads it: 1 switch closed, 0 opened, 2 a spin
      bool switchClosed = false;
      double prevWrapped = 0.0;   // (a spinner with no switch: where it was against its switch angle)

      // ---- what the derived part supplies ----
      virtual double BendKb(double y) const = 0;           // the bend stiffness out of the part's plane at the distance y from the axis
      virtual bool NearStart(double th) const = 0;         // near its start angle: a Hit fires only then
      virtual int LimitKind(int) const { return 0; }       // a touch of this stop as an event: +1 EOS, -1 BOS, 0 none
      virtual void UpdateSwitches() = 0;                   // after a step

      int Surfaces() const { return surfaces; }
      // Where gravity alone hangs it: on a sloped table not theta 0 (6 deg on a 6 deg table)
      double HangAngle() const { return std::atan2(GravB(), GravA()); }
      double PendulumArm() const { const V3 cm = mp.Cm(); return std::sqrt(cm.y * cm.y + cm.z * cm.z); }
      // The most the stored energy can push it with (the speed bound)
      double MaxStoredTorque() const
      {
         double tau = mp.m * g.Length() * PendulumArm() + kBlade * 0.5;
         for (const HingeBody::Cam &cm : cams)
            tau += (cm.F0 + cm.k * cm.r) * cm.r;
         return tau;
      }
      // One contact's probe: only that side's wires and faces, the part at its angle now + dq
      void ProbeSurface(const int surface, const V3 &c, const double R, const double dq, Hit &h) const
      {
         const double th = q + dq;
         const V3 dn = Down(th), fr = Front(th), rel = c - A, cb(rel.Dot(E), rel.Dot(dn), rel.Dot(fr));
         BallSI b;
         b.R = R;
         TubeProbe best;
         best.dist = INFINITY;
         best.delta = -INFINITY;
         for (const Seg &w : wires)
         {
            if (w.surface != surface)
               continue;
            const V3 pts[2] = { w.a, w.b };
            TubeProbe pr = RubberContact::ProbeTube(b, pts, 2, rho, cb);
            if (pr.seg < 0)
               continue;
            pr.delta = R + rho - pr.dist;
            if (pr.delta > best.delta)
               best = pr;
         }
         if (!tris.empty() && surface == 0)
         {
            const TubeProbe pr = RubberContact::ProbeTriangles(b, tris, cand, halfT, cb, nullptr, 0, halfT);
            if (pr.seg >= 0 && pr.delta > best.delta)
               best = pr;
         }
         const V3 nb = best.normal;
         h.delta = best.delta;
         h.n = E * nb.x + dn * nb.y + fr * nb.z;
         h.reff = best.reff > 0.0 ? best.reff : R;
         h.curvDist = 0.0;
         h.grad = 1.0;
         // the contact point (body coordinates), its speed per unit rate of theta (E x Down = Front, E x Front = -Down): the depth's slope with theta is the normal's part of it
         const V3 pb = cb - nb * (R - 0.5 * std::max(best.delta, 0.0));
         h.nc = 1;
         h.c[0].Set(0, nb.z * pb.y - nb.y * pb.z, fr * pb.y - dn * pb.z, E);
         // the bend out of the part's plane, along the normal
         const double oop = std::fabs(nb.z);
         h.kb = oop > 1e-6 ? BendKb(std::max(std::fabs(pb.y), 1e-4)) / (oop * oop) : 0.0;
      }

      // ---- the body's terms ----
      double Viscous() const override { return cViscous; }
      double StoredEnergy(const double th) const override
      {
         double e = GravV(th);
         const double x = th - aBlade;
         if (kBlade > 0.0 && x > 0.0)
            e += 0.5 * kBlade * x * x;
         for (const HingeBody::Cam &cm : cams)
            e += cm.V(cm.Delta(th));
         return e;
      }
      double StoredForce(const double th) const override
      {
         double tau = -GravVp(th);
         const double x = th - aBlade;
         if (kBlade > 0.0 && x > 0.0)
            tau -= kBlade * x;
         for (const HingeBody::Cam &cm : cams)
            tau += cm.F(cm.Delta(th)) * cm.r * std::sin(th - cm.phase);
         return tau;
      }
      double StoredChord(const double th, const double dq, const double h, double &slope, double &lossF, double &lossSlope) const override
      {
         // gravity: closed form, no cancellation (sin a - sin b = 2 cos((a + b) / 2) sin((a - b) / 2))
         const double half = 0.5 * dq, thm = th + half, sinc = std::fabs(half) > 1e-4 ? std::sin(half) / half : 1.0 - half * half / 6.0;
         double chord = GravVp(thm) * sinc;
         slope = 0.5 * GravVpp(thm);
         lossF = lossSlope = 0.0;
         if (kBlade > 0.0)
         {
            const double x0 = th - aBlade, x1 = x0 + dq;
            if (x0 > 0.0 || x1 > 0.0)
            {
               const double a0 = std::max(x0, 0.0), a1 = std::max(x1, 0.0);
               chord += kBlade * 0.5 * (a1 + a0) * (dq != 0.0 ? (a1 - a0) / dq : (x0 > 0.0 ? 1.0 : 0.0));
               slope += 0.5 * kBlade;
            }
         }
         for (const HingeBody::Cam &cm : cams)
         {
            const double d0 = cm.Delta(th), d1 = cm.Delta(th + dq);
            if (d0 <= 0.0 && d1 <= 0.0)
               continue;
            const double sn = std::sin(thm - cm.phase), ddt = -cm.r * sn;   // d delta / d theta at mid-step
            const double Fm = cm.F(std::max(0.5 * (d0 + d1), 0.0));
            if (std::fabs(dq) > 1e-7)
               chord += (cm.V(d1) - cm.V(d0)) / dq;
            else
               chord += Fm * ddt;
            slope += 0.5 * (cm.k + (cm.F0 > 0.0 ? cm.F0 / cm.rampM : 0.0)) * ddt * ddt;
            if (cm.c > 0.0)
            {
               // the damper on the part of the move inside the zone: continuous in the move, and its work has one sign; never pulling
               const double p0 = std::max(d0, 0.0), p1 = std::max(d1, 0.0);
               double Fd = cm.c * (p1 - p0) / h;
               if (Fd < -Fm)
                  Fd = -Fm;
               if (std::fabs(dq) > 1e-12)
               {
                  lossF += Fd * (p1 - p0) / dq;
                  lossSlope += cm.c * ddt * ddt / h;
               }
            }
         }
         return chord;
      }
      double FrictionLimit(const double th, const double om, const V3 &ballForce) const override
      {
         const V3 cm = mp.Cm(), rcm = Down(th) * cm.y + Front(th) * cm.z, load = g * mp.m + ballForce;
         V3 pin = rcm * (-mp.m * om * om) - load;
         pin = pin - E * pin.Dot(E);
         return muPivot * rPin * pin.Length();
      }
      void AfterMove() override { UpdateSwitches(); }
      void SaveState() override
      {
         BodyPart::SaveState();
         svSwitch = switchClosed;
         svWrapped = prevWrapped;
         svSwLog = swLog.size();
      }
      void RestoreState() override
      {
         BodyPart::RestoreState();
         switchClosed = svSwitch;
         prevWrapped = svWrapped;
         swLog.resize(svSwLog);
      }

   protected:
      struct Seg { V3 a, b; int surface; };
      void Begin(const V3 &axisPoint, const V3 &front, const V3 &gravity)
      {
         g = gravity;
         mp = HingeBody::MassProps();
         stops.clear();
         cams.clear();
         wires.clear();
         tris.clear();
         cand.clear();
         q = rate = 0.0;
         A = axisPoint;
         N = front / front.Length();
         Z = V3(0.0, 0.0, 1.0);
         E = N.Cross(Z);
         E = E / E.Length();
         N = Z.Cross(E);
         surfaces = 1;
         surfaceNow = 0;
         kBlade = aBlade = 0.0;
         switchClosed = false;
      }
      void AddWire(const V3 &p, const V3 &u, const double lambda, const bool touchable)
      {
         mp.AddWire(p, u, lambda);
         if (touchable)
            wires.push_back({ p, u, surfaceNow });
      }
      void AddPlate(const V3 &a, const V3 &b, const V3 &c, const V3 &d, const double sigma)
      {
         mp.AddQuad(a, b, c, d, sigma);
         tris.push_back(RubberContact::MakeSurfaceTri(a, b, c));
         tris.push_back(RubberContact::MakeSurfaceTri(a, c, d));
      }
      // the shape is complete: the faces' list, the body's inertia and reach
      void FinishShape(const double reachM)
      {
         cand.resize(tris.size());
         for (size_t i = 0; i < tris.size(); i++)
            cand[i] = (int)i;
         mass = mp.I;
         reach = reachM;
      }
      // the ball contact law: steel on steel (E* of the pair in the law's single-body form), restitution of the part at 1 m/s on a fixed surface
      void SetContactLaw(const double elasticity, const double friction, const BallSI &b)
      {
         rp = RubberParams();
         rp.nu = 0.3;
         rp.eMpa = GroupStep::kSteelEstar * (1.0 - 0.3 * 0.3) * 1e-6;
         rp.betaLayer = 0.0;
         rp.mu = friction;
         rp.twist = true;
         rp.tanDelta = 0.0;
         rp.normalLoss = false;
         rp.gammaKK = RubberContact::CalibrateKK(rp, b.R, b.m, 1.0, std::min(0.99, std::max(0.05, elasticity)));
      }
      std::vector<Seg> wires;
      std::vector<RubberContact::SurfaceTri> tris;
      std::vector<int> cand;
      double rho = 0.0, halfT = 0.0;
      int surfaces = 1, surfaceNow = 0;   // a wire shape's separate contacts: split where a ball can touch two wires at once (the V of a wire W): nearest-only has a force step there, no root

   private:
      V3 Down(const double th) const { return Z * -std::cos(th) + N * std::sin(th); }
      V3 Front(const double th) const { return N * std::cos(th) + Z * std::sin(th); }
      // gravity on the centre of mass, V = -m g . r(theta) = -m (A cos theta + B sin theta) + const
      double GravA() const { const V3 cm = mp.Cm(); return -cm.y * g.Dot(Z) + cm.z * g.Dot(N); }
      double GravB() const { const V3 cm = mp.Cm(); return cm.y * g.Dot(N) + cm.z * g.Dot(Z); }
      double GravV(const double th) const { return -mp.m * (GravA() * std::cos(th) + GravB() * std::sin(th)); }
      double GravVp(const double th) const { return -mp.m * (-GravA() * std::sin(th) + GravB() * std::cos(th)); }
      double GravVpp(const double th) const { return mp.m * (GravA() * std::cos(th) + GravB() * std::sin(th)); }
      bool svSwitch = false;
      double svWrapped = 0.0;
      size_t svSwLog = 0;
   };

   // One side of a part, a surface it carries (the offline bench's; the engine has its own, on the same probe)
   class Shape : public GroupStep::Shape
   {
   public:
      const Part *owner = nullptr;
      int surface = 0;
      void Probe(const V3 &c, const double R, Hit &h) const override { owner->ProbeSurface(surface, c, R, trial != nullptr ? trial[0] : 0.0, h); }
   };
}
