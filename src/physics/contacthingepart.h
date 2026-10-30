// license:GPLv3+
//
// What every part on one axis shares (a gate, a spinner; docs/GateSpinnerDesign.txt sections 3 and 6): the body (contacthingebody.h), its touchable shape (wire
// capsules and sheet triangles in body coordinates), and the ball contact with it.
//
// Pure SI maths with no VPX types (offline tests: runs/2026-10-14-hinge).
//
// Ball contact: the steel-on-steel Hertz + viscoelastic law of the contact model (EvaluateTube) on the part's present geometry, with the part's surface velocity;
// the reaction goes to the body. The part is not rigid: out of its plane it bends like a cantilever from the axis (the derived part gives the stiffness at the
// distance y of the contact from the axis), far softer than the Hertz contact. The bend is a massless spring in series with the contact: each substep the
// overlap splits so that both carry the same force (k_h delta_c^1.5 = k_b x_b), the contact's damping acts on its share of the approach rate, the bend stores
// and returns its energy. A load in the part's plane does not bend it: the compliance scales with the square of the normal's out-of-plane part.
// The bend's stiffness changes with where the ball touches; a bent cantilever is tilted at the load, so the force on the ball has the component
// -1/2 x^2 grad(k_eff) (the gradient of the stored bend energy at fixed bend, the envelope theorem), acting through the ball's centre. Without it the contact
// gains energy when the contact point slides to a stiffer place under load (measured on gates: +6 % / +16 % of the ball's energy).
// Step: the contact is resolved for the reduced mass of the ball and the part's effective mass at the contact (I / lever^2), with the series stiffness.
// Hit: a ball starts touching the part (after it was more than release_gap away) while the part is near its start angle.

#pragma once

#include <cmath>
#include <algorithm>
#include <vector>

#include "physics/rubbercontact.h"
#include "physics/contacthingebody.h"
#include "physics/contactsurface.h"
#include "physics/contactmovingpart.h"
#include "physics/contactprof.h"

namespace HingePart
{
   using RubberContact::V3;
   using RubberContact::BallSI;
   using RubberContact::TubeProbe;
   using RubberContact::RubberParams;
   using RubberContact::RubberState;
   using RubberContact::ContactResult;
   using RubberContact::SurfaceTri;

   using Event = ContactMovingPart::Event;

   // The per-ball state of a ball's contact with a part (owned by the caller with the ball)
   struct BallLink
   {
      RubberState st;
      bool touching = false;   // inside the release gap since the last touch began
   };

   // A part that turns about one fixed axis (a gate, a spinner): the body of contacthingebody.h and its wires and plates. Implements the engine's moving-part interface.
   class Part : public ContactMovingPart
   {
   public:
      HingeBody::Body body;
      RubberParams rp;            // ball contact law
      BallSI ball;                // the ball of the present call (SetBall: balls may differ in size and mass)

      void SetBall(const BallSI &b) override { ball = b; }

      // ---- the moving-part interface (contactmovingpart.h) ----
      // the body's dynamic state, then whatever the derived part adds (its switch state): SaveExtra / LoadExtra
      void SaveState(State &out) const override
      {
         out.used = 0;
         Put(out, steppedThisCycle);
         Put(out, steppedS);
         HingeBody::Body::Dyn d;
         std::memset(&d, 0, sizeof d);   // (the padding too: saved bytes must be the same for the same state)
         body.SaveDyn(d);
         Put(out, d);
         SaveExtra(out);
      }
      void RestoreState(const State &in) override
      {
         StateReader r(in);
         r.Get(steppedThisCycle);
         r.Get(steppedS);
         HingeBody::Body::Dyn d;
         r.Get(d);
         body.LoadDyn(d);
         events.clear();
         LoadExtra(r);
         UpdateWorld();   // the touchable geometry follows the angle
      }
      virtual void SaveExtra(State &) const {}
      virtual void LoadExtra(StateReader &) {}
      void SetGravity(const V3 &g) override { body.SetGravity(g); }
      V3 Origin() const override { return body.A; }
      // the part's own contact with a ball: the one-surface contact of a gate or spinner
      ContactResult ContactItem(int ballIndex, SurfaceContact &it, const V3 &c, const V3 &v, const V3 &w, double h) override { return Contact(ballIndex, it.st, it.partTouching, c, v, w, h); }
      // the velocity and angular velocity of the part's surface at a point
      void SurfaceVelocity(const V3 &point, V3 &v, V3 &w) const override
      {
         v = body.VelocityAt(point);
         w = body.E * body.omega;
      }

      // ---- shape: the part's distance from the axis (m) where the bend is computed; the derived part gives the bend stiffness out of its plane there ----
      virtual double BendKb(double y) const = 0;
      // whether the part is near its start angle (a Hit only fires then)
      virtual bool NearStart() const = 0;
      // the part's switches and stop events after a substep (the geometry is not needed for them)
      virtual void UpdateSwitches() = 0;

      // after body.Step inside an island: the geometry follows the new angle, then the switches
      void AfterStep()
      {
         UpdateWorld();
         UpdateSwitches();
      }
      // One substep h of an island: the balls' reactions of this substep (added by Contact) act during all of it, and the part advances in its OWN substeps within h (its stops
      // and switch spring set them: a steel stop impact lasts ~100 us and needs a few us, far finer than the ball contacts' step). Measured in VPX: a gate stepped with the
      // island's 77 us substep hit its open stop and came back at -7427 rad/s.
      void StepInIsland(double h) override
      {
         double tau;
         V3 F;
         body.TakeExternal(tau, F);
         body.AdvanceWith(h, [this] { CPROF(PartSwitch); UpdateSwitches(); }, tau, F);
         {
            CPROF(PartWorld);
            UpdateWorld();
         }
         steppedThisCycle = true;
         steppedS += h;
      }
      // a part no ball is near: one interval (a physics cycle) in its own substeps, the switches after every substep, the geometry once at the end
      int AdvanceAlone(double dt) override
      {
         const int n = body.AdvanceWith(dt, [this] { UpdateSwitches(); });
         if (n > 0)
            UpdateWorld();
         return n;
      }
      // the swept reach of the part from its axis (m)
      double Reach() const override { return m_reach; }
      double SurfaceSpeedBound(double cycleS) const override { return body.SurfaceSpeedBound(cycleS, m_reach); }

      // ---- geometry at the present angle ----
      void UpdateWorld()
      {
         m_wireW.resize(m_wire.size());
         for (size_t i = 0; i < m_wire.size(); i++)
            m_wireW[i] = { body.ToWorld(m_wire[i].a), body.ToWorld(m_wire[i].b) };
         m_trisW.resize(m_plate.size());
         m_cand.resize(m_plate.size());
         for (size_t i = 0; i < m_plate.size(); i++)
         {
            m_trisW[i] = RubberContact::MakeSurfaceTri(body.ToWorld(m_plate[i].a), body.ToWorld(m_plate[i].b), body.ToWorld(m_plate[i].c));
            m_cand[i] = (int)i;
         }
      }

      // nearest point of the part to a ball centred at c (the present geometry)
      TubeProbe Probe(const V3 &c) const override
      {
         TubeProbe best;
         best.dist = INFINITY;
         best.delta = -INFINITY;
         for (const auto &w : m_wireW)
         {
            const V3 pts[2] = { w.a, w.b };
            const TubeProbe pr = RubberContact::ProbeTube(ball, pts, 2, m_rho, c);
            if (pr.seg >= 0 && pr.delta > best.delta)
               best = pr;
         }
         if (!m_trisW.empty())
         {
            const TubeProbe pr = RubberContact::ProbeTriangles(ball, m_trisW, m_cand, m_halfT, c, nullptr, 0, m_halfT);
            if (pr.seg >= 0 && pr.delta > best.delta)
               best = pr;
         }
         return best;
      }

      // One substep of the contact of a ball (centre c, velocity v, spin w) with the part: the force and torque on the ball, the reaction added to the body.
      // Call it for every ball in reach before body.Step(h). ballIndex numbers the ball for the Hit event.
      ContactResult Contact(int ballIndex, BallLink &link, const V3 &c, const V3 &v, const V3 &w, double h) { return Contact(ballIndex, link.st, link.touching, c, v, w, h); }
      // the same on the engine's per-ball state of the contact (its shear spring and the touching flag of the release-gap hysteresis)
      ContactResult Contact(int ballIndex, RubberState &st, bool &touching, const V3 &c, const V3 &v, const V3 &w, double h)
      {
         ContactResult r;
         const TubeProbe pr = Probe(c);
         if (pr.seg < 0)
            return r;
         // the part's surface velocity at the contact point
         const V3 pc = c - pr.normal * (ball.R - 0.5 * std::max(pr.delta, 0.0));
         const V3 vs = body.VelocityAt(pc);
         const V3 ws = body.E * body.omega;
         // Hit: a new touch after the ball was farther than the release gap, on a part near its start angle (touches of a swinging part are not hits)
         if (pr.overlap && !touching)
         {
            touching = true;
            if (NearStart())
            {
               const double vn = -(v - vs).Dot(pr.normal);
               events.push_back({ Event::Hit, ballIndex, vn });
            }
         }
         else if (!pr.overlap && touching && -pr.delta > m_releaseGap)
            touching = false;
         // the bend in series with the contact: the contact gets its share of the overlap and of the approach rate
         TubeProbe prc = pr;
         RubberParams rpc = rp;
         if (pr.overlap)
         {
            const Split sp = SeriesSplit(pr, pc);
            prc.delta = sp.deltaC;
            rpc.gammaKK = rp.gammaKK * sp.rateShare;
         }
         r = RubberContact::EvaluateTube(rpc, ball, prc, 1.0, v, w, h, st, vs, ws);
         if (r.inContact && prc.delta < pr.delta)
         {
            // the bend's share of the force off the rigid normal: -1/2 x^2 grad k_eff(c) at fixed bend x (central differences over the ball position)
            const double x = pr.delta - prc.delta;
            const double eps = 1e-6;
            V3 grad;
            const V3 ex(eps, 0, 0), ey(0, eps, 0), ez(0, 0, eps);
            grad.x = (KeffAt(c + ex) - KeffAt(c - ex)) / (2.0 * eps);
            grad.y = (KeffAt(c + ey) - KeffAt(c - ey)) / (2.0 * eps);
            grad.z = (KeffAt(c + ez) - KeffAt(c - ez)) / (2.0 * eps);
            if (std::isfinite(grad.x) && std::isfinite(grad.y) && std::isfinite(grad.z))
            {
               // the tilt of the bent surface: still a frictionless normal force on a sphere, so it acts through the ball's centre (no torque on the ball;
               // measured: applying it at the rigid contact point spun the ball up and doubled its energy on a long plate hit from behind)
               const V3 Fb = grad * (-0.5 * x * x);
               r.force = r.force + Fb;
            }
         }
         if (r.inContact)
         {
            // reaction on the part: the force at the contact point plus the couple (twist) the contact puts on the ball, both opposite
            const V3 rB = r.normal * (-(ball.R - r.delta / 2.0));
            const V3 Mc = r.torqueOnBall - rB.Cross(r.force);
            const double tau = body.E.Dot((pc - body.A).Cross(r.force * -1.0)) - body.E.Dot(Mc);
            body.AddExternal(tau, r.force * -1.0);
         }
         return r;
      }

      // the bend stiffness along the normal for a ball centred at c (INFINITY where it does not bend)
      double KeffAt(const V3 &c) const
      {
         const TubeProbe pr = Probe(c);
         if (pr.seg < 0)
            return INFINITY;
         const V3 pc = c - pr.normal * (ball.R - 0.5 * std::max(pr.delta, 0.0));
         return BendStiffness(pc, pr.normal);
      }

      // Bending stiffness of the part at a contact point, along the contact normal (N/m); infinite for a load in the part's plane
      double BendStiffness(const V3 &pc, const V3 &n) const
      {
         const V3 rel = pc - body.A;
         const double y = std::max(std::abs(rel.Dot(body.Down(body.theta))), 1e-4);   // distance of the contact from the axis, in the part's plane
         const double outOfPlane = std::abs(n.Dot(body.Front(body.theta)));
         if (outOfPlane < 1e-6)
            return INFINITY;
         return BendKb(y) / (outOfPlane * outOfPlane);
      }

      struct Split
      {
         double deltaC = 0.0;     // contact compression
         double bend = 0.0;       // bend along the normal
         double kh = 0.0, kb = 0.0;
         double rateShare = 1.0;  // the contact's share of the approach rate
      };
      // Splits the overlap of a probe between the Hertz contact and the bend so that both carry the same force (Newton on k_h u^1.5 = k_b (delta - u))
      Split SeriesSplit(const TubeProbe &pr, const V3 &pc) const
      {
         Split sp;
         const double delta = std::max(pr.delta, 0.0);
         sp.kh = (4.0 / 3.0) * HingeBody::kSteelEstar * std::sqrt(pr.reff > 0.0 ? pr.reff : ball.R);
         sp.kb = BendStiffness(pc, pr.normal);
         if (delta <= 0.0 || !std::isfinite(sp.kb))
         {
            sp.deltaC = delta;
            return sp;
         }
         double u = std::min(delta, std::pow(sp.kb * delta / sp.kh, 2.0 / 3.0));
         for (int i = 0; i < 20; i++)
         {
            const double su = std::sqrt(std::max(u, 0.0));
            const double f = sp.kh * u * su - sp.kb * (delta - u);
            const double df = 1.5 * sp.kh * su + sp.kb;
            const double du = f / df;
            u = std::min(delta, std::max(0.0, u - du));
            if (std::abs(du) < 1e-12 * (delta + 1e-9))
               break;
         }
         sp.deltaC = u;
         sp.bend = delta - u;
         const double kt = 1.5 * sp.kh * std::sqrt(std::max(u, 0.0));   // the contact's tangent stiffness
         sp.rateShare = sp.kb / (sp.kb + kt);
         return sp;
      }

      // elastic energy stored in a ball's contact with the part (contact + bend), J
      double ContactEnergy(const V3 &c) const
      {
         const TubeProbe pr = Probe(c);
         if (pr.seg < 0 || pr.delta <= 0.0)
            return 0.0;
         const V3 pc = c - pr.normal * (ball.R - 0.5 * pr.delta);
         const Split sp = SeriesSplit(pr, pc);
         return 0.4 * sp.kh * std::pow(sp.deltaC, 2.5) + (std::isfinite(sp.kb) ? 0.5 * sp.kb * sp.bend * sp.bend : 0.0);
      }

      // the step the part needs for itself (the island takes the smaller of this and the ball contacts')
      double ChooseStep(double hMax) const { return body.ChooseStep(hMax); }

      // The step a ball's contact with the part needs (design 6.3), by the same rule as every other surface (contactsystem.cpp NextStep): the reduced mass of the ball and the part's effective mass
      // at the contact, the series stiffness of contact and bend; resting (closing speed under 0.02 m/s, compression under 3 x the ball's weight compression): omega * h <= 0.8, impact: 0.2 at the
      // compression the closing speed reaches. Before the touch the step ends at the arrival (closing speed plus a bound on the closing acceleration: gravity on the ball and the part's own
      // angular acceleration).
      double BallStep(const V3 &c, const V3 &v, double hMax) const override
      {
         if ((c - body.A).Length() > Reach() + ball.R + std::max(v.Length(), 0.1) * hMax)
            return hMax;
         const TubeProbe pr = Probe(c);
         if (pr.seg < 0)
            return hMax;
         const V3 n = pr.normal;
         const V3 pc = c - n * (ball.R - 0.5 * std::max(pr.delta, 0.0));
         const double lever = (pc - body.A).Cross(n).Dot(body.E);
         const double invMeff = lever * lever / body.mp.I;
         const double mu = 1.0 / (1.0 / ball.m + invMeff);
         const double vc = -(v - body.VelocityAt(pc)).Dot(n);   // closing speed
         // gravity on the ball, the part's own swing (pendulum) and, for a part that spins, the centripetal acceleration of the contact point
         const double ac = 9.81 + std::abs(lever) * (body.mp.m * 9.81 * body.PendulumArm() / body.mp.I) + body.omega * body.omega * (pc - body.A).Length();
         const double k = (4.0 / 3.0) * HingeBody::kSteelEstar * std::sqrt(pr.reff > 0.0 ? pr.reff : ball.R);
         const double kb = BendStiffness(pc, n);
         constexpr double kRestSpeed = 0.02, kImpactOmegaH = 0.2, kRestOmegaH = 0.8;
         const double dRest = std::pow(ball.m * 9.81 / k, 2.0 / 3.0);   // k delta^1.5 = the ball's weight
         const bool resting = std::abs(vc) <= kRestSpeed && pr.delta <= 3.0 * dRest;
         auto hContact = [&](double vImpact) {
            // the contact's stiffness at the compression it rests at, or at the deepest a Hertz-only impact would reach (an upper bound), in series with the bend
            const double vi = std::max(vImpact, 0.0);
            const double dMax = resting ? std::max(pr.delta, dRest) : std::max(pr.delta, 0.0) + std::pow(1.25 * mu * vi * vi / k, 0.4);
            const double kt = 1.5 * k * std::sqrt(std::max(dMax, 1e-12));
            const double ks = std::isfinite(kb) ? 1.0 / (1.0 / kt + 1.0 / kb) : kt;
            return (resting ? kRestOmegaH : kImpactOmegaH * ContactOmegaHScale()) / std::sqrt(ks / mu);
         };
         double h = hMax;
         if (pr.delta > 0.0)
            h = std::min(h, hContact(vc));
         else
         {
            const double gap = -pr.delta;
            const double tArr = (-vc + std::sqrt(vc * vc + 4.0 * ac * gap)) / (2.0 * ac);
            if (tArr < h)
               h = std::max(tArr, hContact(vc + ac * tArr));
         }
         return std::max(h, HingeBody::Body::kMinStep);
      }

      double AngleDeg() const { return body.theta * 180.0 / M_PI; }

   protected:
      struct Seg { V3 a, b; };
      struct Tri { V3 a, b, c; };
      std::vector<Seg> m_wire, m_wireW;
      std::vector<Tri> m_plate;
      std::vector<SurfaceTri> m_trisW;
      std::vector<int> m_cand;
      double m_rho = 0.0, m_halfT = 0.0, m_reach = 0.0, m_releaseGap = 0.2e-3;

      // the ball contact law: steel on steel (E* of the pair in the law's single-body form), restitution of the part at 1 m/s on a fixed surface
      void SetupContactLaw(double elasticity, double friction)
      {
         rp = RubberParams();
         rp.nu = 0.3;
         rp.eMpa = HingeBody::kSteelEstar * (1.0 - 0.3 * 0.3) * 1e-6;
         rp.betaLayer = 0.0;
         rp.mu = friction;
         rp.twist = true;
         rp.tanDelta = 0.0;
         rp.normalLoss = false;
         rp.gammaKK = RubberContact::CalibrateKK(rp, ball.R, ball.m, 1.0, std::min(0.99, std::max(0.05, elasticity)));
      }
      void ClearShape()
      {
         m_wire.clear();
         m_plate.clear();
      }
      void AddWire(const V3 &p, const V3 &q, double lambda, bool touchable)
      {
         body.mp.AddWire(p, q, lambda);
         if (touchable)
            m_wire.push_back({ p, q });
      }
      void AddPlate(const V3 &a, const V3 &b, const V3 &c, const V3 &d, double sigma)
      {
         body.mp.AddQuad(a, b, c, d, sigma);
         m_plate.push_back({ a, b, c });
         m_plate.push_back({ a, c, d });
      }
   };
}
