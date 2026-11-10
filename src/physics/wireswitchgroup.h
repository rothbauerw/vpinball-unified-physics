// license:GPLv3+
//
// A rollover wire switch as a part of the group step (runs/2026-11-03-switch-requirements/design-wire.md, topic 1): a rigid steel wire on the lever of a micro switch, one angle
// about the switch's hinge below the playfield (hingegroup.h, bodypart.h): a gate turned on its side.
//   spring   the micro switch's lever spring, a preloaded torsion spring tau = -k (theta - a), a < 0 its free angle beyond the rest stop (the hinge part's one-sided blade);
//            k and a set by the vertical force at the peak: spring_n at first touch and spring_full_n at full travel, the wire's own weight included
//   stops    rest (wire up, theta 0) and full travel (the peak stop_below_mm below the playfield): BodyStop, the Hertz steel stop with viscoelastic restitution, as gates
//   pivot    Coulomb pin friction, gravity with the cabinet's acceleration, as gates
//   switch   closes when the peak has travelled close_mm; opens when it rises back release_mm past that point (the micro switch's differential travel). Hit on close, Unhit on open
// Frame (hingegroup.h): the wire is fixed to the micro switch's lever and turns with it about the lever's pivot (PartData S7: the switch on a bracket flush under the playfield,
// the pivot 19.3 mm below the playfield top). The axis runs ACROSS the lane (perpendicular to the wire's plane), N points along the lane from the wire to the pivot, so a positive
// theta swings the wire down and away from the pivot. The shape, the pivot, the lever end and the peak are fixed per shape (wireshape.h MakeCentreline, placed by Place); the
// caller passes them in metres. Every capsule is touchable and carries mass; the lever's own mass is left out (user 2026-11-09: about 0.16 g within 20 mm of the pivot,
// 5 - 10 % of the inertia). Out of the wire's plane (along the axis) the wire bends as a cantilever from the lever end. Pure SI.

#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "physics/hingegroup.h"
#include "physics/wireshape.h"

namespace WireSwitchGroup
{
   using GroupStep::M3;
   using GroupStep::V3;
   using RubberContact::BallSI;

   struct Settings
   {
      double heightOffsetMm = 0.0;    // the switch's mount raised (+) or lowered (-) from the agreed bracket height (pivot 19.3 mm below the playfield top): the whole switch and wire move
      double stopBelowMm = 3.0;       // full travel: the peak this far below the playfield (user 2026-11-09: a ball cannot press it below the surface, so the wire never bottoms out in play)
      double closeMm = 3.0;           // the switch closes at this travel of the peak
      double releaseMm = 0.5;         // and opens when the peak rises back this far past the close point
      double springN = 0.02;          // net vertical force at the peak at first touch, the wire's weight included
      double springFullN = 0.03;      // at full travel
      double wireDMm = 1.0;           // the steel wire's diameter (PartData S1): the touched tube, the drawn tube and the mass; under 2.6 mm (the 1.3 mm bends fold above that)
      double bendStepDeg = 30.0;      // bend_step_deg: each 1.3 mm bend drawn and touched as straight pieces of at most this angle (runs/2026-11-09-bend-steps/design.md; user 2026-11-09: the
                                      // bends need only look right at play distance). 30: a 90 deg bend in 3 pieces of 0.68 mm, 0.044 mm off the true arc. 5 = the shapes before 2026-11-09
      double archStepDeg = 10.0;      // arch_step_deg: D's 14.7 mm arch the same way. 10: pieces of 2.6 mm, 0.056 mm off the true arc. 5 = the shape before 2026-11-09
      int wireBend = 1;               // wire_bend (default 1, on while tested; runs/2026-11-09-wire-bend/design.md): 1 = the wire bends in every direction as a round steel wire clamped at the lever end
                                      // (bending E I about both section axes, torsion G J, G = E / 2.6, J = 2 I; unit-load method over the wire from the touch to the clamp): along the normal in series with
                                      // the contact (kb) and across it in series with the shear spring (the stiffer tangential direction); 0 = out of the wire's plane only, a straight cantilever (in its plane rigid)
      double stopArmMm = 5.0, stopEdgeMm = 0.4;   // the lever on the switch body's stops, as the gate's
      double stopE = 0.2;             // both stops' restitution, ESTIMATE (unmeasured): the thin spring-steel wire loses most of an impact in its own bending (user 2026-11-04; 0.5 rebounded past the close point at 3 m/s)
      double muPivot = 0.3;
      double density = 7850.0;
      double steelEGpa = 200.0;
      double elasticity = 0.3, friction = 0.02;   // ball contact: the wire's restitution at 1 m/s, friction
   };

   class Part : public HingeGroup::Part
   {
   public:
      Settings set;

      // body: the placed wire in metres (wireshape.h Place: chain, pivot, lane from the switch toward the wire, lever end, peak); surfaceZ: the playfield top under the wire, m
      void Build(const Settings &s, const double surfaceZ, const WireShape::WireBody &body_, const BallSI &b, const V3 &gravity)
      {
         set = s;
         const double mm = 1e-3;
         const std::vector<WireShape::Capsule> &chain = body_.caps;
         const V3 lane(body_.lane.x, body_.lane.y, 0.0);
         Begin(V3(body_.pivot.x, body_.pivot.y, body_.pivot.z), lane * -1.0, gravity);
         muPivot = s.muPivot;
         rho = chain.empty() ? 0.5 * s.wireDMm * mm : chain[0].r;
         rPin = rho;
         const double lam = s.density * GroupStep::kPi * rho * rho;
         // body coordinates at theta 0: x along the axis (E), y = -height above the hinge, z = distance along N from the hinge (the wire is in front of it, z < 0)
         const auto body = [&](const WireShape::V3 &p) { const V3 r = V3(p.x, p.y, p.z) - A; return V3(r.Dot(E), -r.Dot(Z), r.Dot(N)); };
         double reachM = 0.0;
         for (const WireShape::Capsule &c : chain)
            for (const V3 &p : { body(c.a), body(c.b) })
               reachM = std::max(reachM, p.Length() + rho);
         const V3 pk = body(body_.peak);   // the named peak (the wire file's corner), not searched for
         peakY = pk.y;
         peakZ = pk.z;
         for (const WireShape::Capsule &c : chain)
            AddWire(body(c.a), body(c.b), lam, true);
         tip = V3(body_.leverEnd.x, body_.leverEnd.y, body_.leverEnd.z);
         tipBody = body(body_.leverEnd);
         FinishShape(reachM);
         BuildBend(chain, body);

         // the travel of the peak: rest stop at theta 0, full travel where it is stop_below under the playfield, the switch's make and break
         q = rate = 0.0;
         thFull = AngleAtDrop(body_.peak.z - surfaceZ + s.stopBelowMm * mm);
         switchMake = AngleAtDrop(s.closeMm * mm);
         switchBreak = AngleAtDrop(std::max(s.closeMm - s.releaseMm, 0.0) * mm);
         const auto stop = [&](const double at, const int side) {
            GroupStep::BodyStopPlace st;
            st.at = at;
            st.side = side;
            st.lever = s.stopArmMm * mm;
            st.Setup(mp.I, std::sqrt(rho * s.stopEdgeMm * mm), s.stopE, 0.1);
            return st;
         };
         stopRest = PlaceStop(stop(0.0, -1));
         stopFull = PlaceStop(stop(thFull, 1));

         // the spring from the net vertical force at the peak: k (theta - a) - tau_g(theta) = F(theta) * arm(theta) at theta 0 and at full travel (tau_g: gravity's torque toward +theta)
         kBlade = aBlade = 0.0;
         const double t0 = s.springN * PeakArm(0.0) + StoredForce(0.0), tf = s.springFullN * PeakArm(thFull) + StoredForce(thFull);
         if (tf > t0 && t0 > 0.0)
         {
            kBlade = (tf - t0) / thFull;
            aBlade = -t0 / kBlade;
         }
         SetContactLaw(s.elasticity, s.friction, b);
         CycleBegin();
      }
      double PeakDrop(const double th) const { return PeakHeight(0.0) - PeakHeight(th); }   // how far the peak has travelled down, m
      double FullAngle() const { return thFull; }
      double PeakLane(const double th) const { return -(peakY * std::sin(th) + peakZ * std::cos(th)); }   // the peak along the lane away from the pivot, m
      V3 TipWorld() const { return tip; }
      double MakeAngle() const { return switchMake; }
      double BreakAngle() const { return switchBreak; }


      // The trigger's Hit is the switch closing (swLog), not a touch
      bool NearStart(double) const override { return false; }
      void UpdateSwitches() override
      {
         if (!switchClosed && q >= switchMake)
         {
            switchClosed = true;
            swLog.push_back(1);
         }
         else if (switchClosed && q < switchBreak)
         {
            switchClosed = false;
            swLog.push_back(0);
         }
      }
      // bend out of the wire's plane at the distance y from where it is held: one wire, a cantilever, k = 3 E I / y^3, I = pi d^4 / 64
      double BendKb(const double y) const override
      {
         const double d = 2.0 * rho, I = GroupStep::kPi * d * d * d * d / 64.0;
         return 3.0 * set.steelEGpa * 1e9 * I / (y * y * y);
      }
      // the wire's plane holds N, out of it is the axis (x); held at the lever end
      double BendAlongNormal(const V3 &nb, const V3 &pb) const override
      {
         const double oop = std::fabs(nb.x);
         return oop > 1e-6 ? BendKb(std::max((pb - tipBody).Length(), 1e-4)) / (oop * oop) : 0.0;
      }
      // [Triggers] wire_bend: the bend in series with a touch at pb (body coordinates) along nb on wire piece w. Along the normal kb = 1 / (n^T C n); across it the stiffer
      // tangential direction, kt_b = 1 / (the smaller eigenvalue of C's block in the contact plane). The coupling between the directions is left out (design section 3).
      // A touch on the lever (before the clamp) or on no wire piece is rigid (0). wire_bend 0: out of the plane only (BendAlongNormal), no shear bend
      void Bend(const V3 &nb, const V3 &pb, const int w, double &kb, double &ktb) const override
      {
         if (set.wireBend == 0)
         {
            HingeGroup::Part::Bend(nb, pb, w, kb, ktb);
            return;
         }
         kb = ktb = 0.0;
         M3 C;
         if (!Compliance(w, pb, C))
            return;
         const double cn = Quad(C, nb, nb);
         kb = cn > kRigidC ? 1.0 / cn : 0.0;
         // the contact plane: e1 across the wire's plane (x) less its share along n (else the wire's own direction), e2 = n x e1
         V3 e1 = V3(1.0, 0.0, 0.0) - nb * nb.x;
         if (e1.Length() < 1e-6)
            e1 = bend[(size_t)w].t - nb * nb.Dot(bend[(size_t)w].t);
         const double l1 = e1.Length();
         if (l1 < 1e-12)
            return;
         e1 = e1 / l1;
         const V3 e2 = nb.Cross(e1);
         const double c11 = Quad(C, e1, e1), c22 = Quad(C, e2, e2), c12 = Quad(C, e1, e2);
         const double lo = 0.5 * (c11 + c22) - std::sqrt(0.25 * (c11 - c22) * (c11 - c22) + c12 * c12);
         ktb = lo > kRigidC ? 1.0 / lo : 0.0;
      }
      // the wire's compliance (3 x 3, m/N, body coordinates) at the centreline point of piece w nearest pb, clamped at the lever end: C(p) = integral of S^T W S ds from the clamp
      // to p, S = Skew(p - r(s)), W = (I - t t^T) / EI + t t^T / GJ. With the sums over the whole pieces from the clamp (BuildBend) and this piece's part up to p, closed form
      // (W and t constant on a piece, r linear): C = -P M0 P + X + X^T - M2, X = P M1, P = Skew(p), M0 = int W, M1 = int W Skew(r), M2 = int Skew(r) W Skew(r). False on the lever
      bool Compliance(const int w, const V3 &pb, M3 &C) const
      {
         if (bendFrom < 0 || w < bendFrom || w >= (int)bend.size())
            return false;
         const BendPiece &bp = bend[(size_t)w];
         const double l = std::min(std::max((pb - clampBody - bp.a).Dot(bp.t), 0.0), bp.L);
         M3 m0, m1, m2;
         bp.Part(l, m0, m1, m2);
         m0 = Add(m0, bp.s0);
         m1 = Add(m1, bp.s1);
         m2 = Add(m2, bp.s2);
         const M3 P = M3::Skew(bp.a + bp.t * l), X = P * m1;
         C = Add(Add(Scale(P * m0 * P, -1.0), Add(X, X.T())), Scale(m2, -1.0));
         return true;
      }

   private:
      // the peak's height above the hinge at theta (Down(th) . Z = -cos, Front(th) . Z = sin) and its vertical lever, |d height / d theta|
      double PeakHeight(const double th) const { return -peakY * std::cos(th) + peakZ * std::sin(th); }
      double PeakArm(const double th) const { return std::fabs(peakY * std::sin(th) + peakZ * std::cos(th)); }
      // the angle where the peak has dropped d: -y cos(th) + z sin(th) = h0 - d, the root nearest 0 on the pressing side
      double AngleAtDrop(const double d) const
      {
         const double r = std::sqrt(peakY * peakY + peakZ * peakZ), phi = std::atan2(peakZ, -peakY), c = (PeakHeight(0.0) - d) / r;   // height = r cos(th - phi)
         return phi + std::acos(std::max(-1.0, std::min(1.0, c)));
      }
      double peakY = 0.0, peakZ = 0.0, thFull = 0.0, switchMake = 0.0, switchBreak = 0.0;
      V3 tip, tipBody;   // the lever end (where the wire leaves the lever): world, body coordinates
      int stopRest = -1, stopFull = -1;

      // [Triggers] wire_bend: one wire piece from the clamp's side, coordinates from the clamp: its start a, direction t, length L, its W and the products of the closed form,
      // and the sums s0, s1, s2 (M0, M1, M2) over the whole pieces between the clamp and its start
      static constexpr double kRigidC = 1e-12;   // m/N: a compliance below this is a rigid wire (kb, kt_b above 1e12 N/m: the contact alone)
      static M3 Add(const M3 &x, const M3 &y) { M3 m; for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) m.a[i][j] = x.a[i][j] + y.a[i][j]; return m; }
      static M3 Scale(const M3 &x, const double s) { M3 m; for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) m.a[i][j] = s * x.a[i][j]; return m; }
      static double Quad(const M3 &C, const V3 &u, const V3 &v)
      {
         const double uu[3] = { u.x, u.y, u.z }, vv[3] = { v.x, v.y, v.z };
         double s = 0.0;
         for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
               s += uu[i] * C.a[i][j] * vv[j];
         return s;
      }
      struct BendPiece
      {
         V3 a, t;
         double L = 0.0;
         M3 W, AWA, AWT, TWT;   // W, Skew(a) W Skew(a), Skew(a) W Skew(t) + Skew(t) W Skew(a), Skew(t) W Skew(t)
         M3 s0, s1, s2;
         // the piece from its start over the length l: M0 = W l, M1 = W Skew(a l + t l^2 / 2), M2 = AWA l + AWT l^2 / 2 + TWT l^3 / 3
         void Part(const double l, M3 &m0, M3 &m1, M3 &m2) const
         {
            m0 = Scale(W, l);
            m1 = W * M3::Skew(a * l + t * (0.5 * l * l));
            m2 = Add(Add(Scale(AWA, l), Scale(AWT, 0.5 * l * l)), Scale(TWT, l * l * l / 3.0));
         }
      };
      std::vector<BendPiece> bend;   // per wire piece (the chain's order); those before bendFrom lie on the lever
      int bendFrom = -1;
      V3 clampBody;                  // the clamp: the chain's node nearest the lever end, body coordinates
      // at Build: the clamp, each piece's terms and the sums from the clamp (design section 2: E = steel_e_gpa, G = E / 2.6, I = pi d^4 / 64, J = 2 I, d = the wire's diameter)
      template <class ToBody> void BuildBend(const std::vector<WireShape::Capsule> &chain, const ToBody &body)
      {
         bend.clear();
         bendFrom = -1;
         const int n = (int)chain.size();
         if (n == 0)
            return;
         int c0 = 0;
         double best = INFINITY;
         for (int i = 0; i <= n; i++)
         {
            const double dist = (body(i < n ? chain[(size_t)i].a : chain[(size_t)n - 1].b) - tipBody).Length();
            if (dist < best)
            {
               best = dist;
               c0 = i;
            }
         }
         clampBody = body(c0 < n ? chain[(size_t)c0].a : chain[(size_t)n - 1].b);
         const double d = 2.0 * rho, I = GroupStep::kPi * d * d * d * d / 64.0, E = set.steelEGpa * 1e9, EI = E * I, GJ = E / 2.6 * 2.0 * I;
         bendFrom = c0;
         bend.resize((size_t)n);
         M3 s0, s1, s2;
         for (int k = c0; k < n; k++)
         {
            BendPiece &bp = bend[(size_t)k];
            const V3 a = body(chain[(size_t)k].a) - clampBody, b = body(chain[(size_t)k].b) - clampBody;
            bp.a = a;
            bp.L = (b - a).Length();
            bp.t = bp.L > 0.0 ? (b - a) / bp.L : V3(0.0, 1.0, 0.0);
            const M3 tt = M3::Outer(bp.t, bp.t), Wt = M3::Tangent(bp.t);
            bp.W = Add(Scale(Wt, 1.0 / EI), Scale(tt, 1.0 / GJ));
            const M3 A = M3::Skew(a), T = M3::Skew(bp.t);
            bp.AWA = A * bp.W * A;
            bp.AWT = Add(A * bp.W * T, T * bp.W * A);
            bp.TWT = T * bp.W * T;
            bp.s0 = s0;
            bp.s1 = s1;
            bp.s2 = s2;
            M3 m0, m1, m2;
            bp.Part(bp.L, m0, m1, m2);
            s0 = Add(s0, m0);
            s1 = Add(s1, m1);
            s2 = Add(s2, m2);
         }
      }
   };
}
