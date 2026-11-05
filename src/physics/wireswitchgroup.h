// license:GPLv3+
//
// A rollover wire switch as a part of the group step (runs/2026-11-03-switch-requirements/design-wire.md, topic 1): a rigid steel wire on the lever of a micro switch, one angle
// about the switch's hinge below the playfield (hingegroup.h, bodypart.h): a gate turned on its side.
//   spring   the micro switch's lever spring, a preloaded torsion spring tau = -k (theta - a), a < 0 its free angle beyond the rest stop (the hinge part's one-sided blade);
//            k and a set by the vertical force at the peak: spring_n at first touch and spring_full_n at full travel, the wire's own weight included
//   stops    rest (wire up, theta 0) and full travel (the peak stop_below_mm below the playfield): BodyStop, the Hertz steel stop with viscoelastic restitution, as gates
//   pivot    Coulomb pin friction, gravity with the cabinet's acceleration, as gates
//   switch   closes when the peak has travelled close_mm; opens when it rises back release_mm past that point (the micro switch's differential travel). Hit on close, Unhit on open
// Frame (hingegroup.h): the wire is fixed to the tip of the micro switch's lever and turns with it about the lever's pivot (user's photos, 2026-11-04). The switch hangs under the playfield below the near end
// of the slot, the axis runs ACROSS the lane (perpendicular to the wire's plane), N points along the lane from the wire to the pivot, so a positive theta swings the wire down and away from the pivot.
// The lever tip holds the chain node leverNode (user's layouts 2026-11-04: A, B, D a chain end, C the low end of its visible strand); the pivot sits hinge_back_mm along the lane behind it
// (away from the wire's top) and hinge_below_mm below the playfield (both UNVERIFIED, estimated from the product photo). Shape: the mesh's centreline (wireshape.h), every capsule touchable and carrying mass;
// the caller passes it in metres. Out of the wire's plane (along the axis) the wire bends as a cantilever from the lever tip. Pure SI.
// The Inder trigger is a flat steel blade standing on edge, not a wire (2026-11-04, bench runs/2026-11-04-inder): the mesh VPX draws (the caller passes it) is its collidable shape, each triangle a face
// (PoseBlade, the triangle Build) with a face skin edge_mm (ESTIMATE), mass from its two flat side faces, plate bending; the same stops, spring, pivot and switch as the wires.

#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "physics/hingegroup.h"
#include "physics/wireshape.h"

namespace WireSwitchGroup
{
   using GroupStep::V3;
   using RubberContact::BallSI;

   struct Settings
   {
      double peakMm = 6.0;            // wire top above the playfield at rest
      double stopBelowMm = 1.0;       // full travel: the peak this far below the playfield
      double closeMm = 3.0;           // the switch closes at this travel of the peak
      double releaseMm = 0.5;         // and opens when the peak rises back this far past the close point
      double springN = 0.02;          // net vertical force at the peak at first touch, the wire's weight included
      double springFullN = 0.03;      // at full travel
      double wireDMm = 1.4;           // 0: the mesh's tube + WireThickness (the caller's WireShape::Pose)
      double hingeBackMm = 18.0;      // UNVERIFIED (product photo): the pivot, along the lane behind the switch leg's end, away from the span
      double hingeBelowMm = 15.0;     // UNVERIFIED (product photo): the pivot below the playfield
      int leverNode = -1;             // the chain node the lever tip holds (0 = the chain's first point, i = the end of capsule i - 1); -1: the chain's lower end. Per shape, user's layouts 2026-11-04
      double stopArmMm = 5.0, stopEdgeMm = 0.4;   // the lever on the switch body's stops, as the gate's
      double stopE = 0.2;             // both stops' restitution, ESTIMATE (unmeasured): the thin spring-steel wire loses most of an impact in its own bending (user 2026-11-04; 0.5 rebounded past the close point at 3 m/s)
      double muPivot = 0.3;
      double density = 7850.0;
      double steelEGpa = 200.0;
      double elasticity = 0.3, friction = 0.02;   // ball contact: the wire's restitution at 1 m/s, friction
      double edgeMm = 0.05;           // a blade (Inder): the faces' skin and edge radius, ESTIMATE (the triangle probe needs one > 0; the blade acts 2 x edge_mm thicker)
      int wireModel = 1;              // [Triggers] wire_model: 1 = a massless lever (runs/2026-11-04-massless-wire/design.md), 0 = the hinge body with its mass. Wires only: a blade keeps the body
      int wireSum = 0;                // [Triggers] wire_sum: 1 = every wire piece pressing the ball adds its own contact (one summed force where two pieces meet), 0 = the deepest piece alone
   };

   // A blade trigger as VPX draws it (the caller passes the drawn vertices, Trigger::GenerateMesh: the author's rotation, x * Scale X, y * Scale Y, z unscaled on the trigger's
   // Surface, WireThickness along the normals), VPX units, converted to metres. lane: the blade's length direction (the mesh's y as drawn), horizontal. tri: 3 points per triangle;
   // thick: the extent across the lane; tip: the lever tip, the lowest point within 10 % of the length of the -lane end, on the blade's centre plane; top0/top1: the top edge's ends
   struct Blade
   {
      std::vector<V3> tri;
      double thick = 0.0;
      V3 tip, top0, top1;
   };
   template <class Index>
   inline Blade PoseBlade(const float *verts, const unsigned nVerts, const unsigned stride, const Index *idx, const unsigned nIdx, const V3 &laneDir, const double metresPerUnit)
   {
      Blade bl;
      V3 lane(laneDir.x, laneDir.y, 0.0);
      lane = lane / lane.Length();
      const V3 side = lane.Cross(V3(0.0, 0.0, 1.0));
      const auto P = [&](const unsigned i) { const float *v = verts + (size_t)i * stride; return V3(v[0] * metresPerUnit, v[1] * metresPerUnit, v[2] * metresPerUnit); };
      double s0 = 1e30, s1 = -1e30, c0 = 1e30, c1 = -1e30, z1 = -1e30;
      for (unsigned i = 0; i < nVerts; i++)
      {
         const V3 p = P(i);
         s0 = std::min(s0, p.Dot(lane)); s1 = std::max(s1, p.Dot(lane));
         c0 = std::min(c0, p.Dot(side)); c1 = std::max(c1, p.Dot(side));
         z1 = std::max(z1, p.z);
      }
      bl.thick = c1 - c0;
      for (unsigned i = 0; i < nIdx; i++)
         bl.tri.push_back(P((unsigned)idx[i]));
      const double cm = 0.5 * (c0 + c1);
      const auto onPlane = [&](const V3 &p) { return p - side * (p.Dot(side) - cm); };
      V3 tip(0.0, 0.0, 1e30);
      for (unsigned i = 0; i < nVerts; i++)
      {
         const V3 p = P(i);
         const double s = p.Dot(lane);
         if (s <= s0 + 0.1 * (s1 - s0) && (p.z < tip.z - 1e-9 || (std::fabs(p.z - tip.z) <= 1e-9 && s < tip.Dot(lane))))
            tip = p;
      }
      bl.tip = onPlane(tip);
      // the top edge's ends: on the centre plane, at the lane's ends, at the top's height (through the tip's point across the lane)
      const V3 base = onPlane(tip) - lane * (tip.Dot(lane));
      bl.top0 = V3(base.x + lane.x * s0, base.y + lane.y * s0, z1);
      bl.top1 = V3(base.x + lane.x * s1, base.y + lane.y * s1, z1);
      return bl;
   }

   class Part : public HingeGroup::Part
   {
   public:
      Settings set;

      // surfaceZ: the playfield's height under the wire. chain: the wire's traced centreline (wireshape.h Build, converted to metres, world frame), each capsule touchable and carrying its mass;
      // rho from the chain (the Pose's wire_d_mm). The lever tip: Settings::leverNode; the lane runs from it toward the wire's top
      void Build(const Settings &s, const double surfaceZ, const std::vector<WireShape::Capsule> &chain, const BallSI &b, const V3 &gravity)
      {
         set = s;
         const double mm = 1e-3;
         const V3 e0(chain.front().a.x, chain.front().a.y, chain.front().a.z), e1(chain.back().b.x, chain.back().b.y, chain.back().b.z);
         const int nNode = (int)chain.size() + 1;
         const auto node = [&](const int i) { const WireShape::V3 &p = i == 0 ? chain.front().a : chain[i - 1].b; return V3(p.x, p.y, p.z); };
         tip = s.leverNode >= 0 && s.leverNode < nNode ? node(s.leverNode) : e0.z <= e1.z ? e0 : e1;
         // the lane runs from the lever tip toward the wire's top (the pivot sits behind the tip, away from it)
         V3 top = node(0);
         for (int i = 1; i < nNode; i++)
            if (node(i).z > top.z)
               top = node(i);
         V3 lane(top.x - tip.x, top.y - tip.y, 0.0);
         lane = lane / lane.Length();
         Begin(V3(tip.x, tip.y, surfaceZ) - lane * (s.hingeBackMm * mm) - V3(0.0, 0.0, s.hingeBelowMm * mm), lane * -1.0, gravity);
         muPivot = s.muPivot;
         sumWires = s.wireSum != 0;
         rho = chain.empty() ? 0.5 * s.wireDMm * mm : chain[0].r;
         rPin = rho;
         const double lam = s.density * GroupStep::kPi * rho * rho;
         // body coordinates at theta 0: x along the axis (E), y = -height above the hinge, z = distance along N from the hinge (the wire is in front of it, z < 0); the peak is the highest centreline point
         const auto body = [&](const WireShape::V3 &p) { const V3 r = V3(p.x, p.y, p.z) - A; return V3(r.Dot(E), -r.Dot(Z), r.Dot(N)); };
         double reachM = 0.0;
         peakY = INFINITY;
         peakZ = 0.0;
         for (const WireShape::Capsule &c : chain)
            for (const V3 &p : { body(c.a), body(c.b) })
            {
               reachM = std::max(reachM, p.Length() + rho);
               if (p.y - rho < peakY)
               {
                  peakY = p.y - rho;
                  peakZ = p.z;
               }
            }
         for (const WireShape::Capsule &c : chain)
            AddWire(body(c.a), body(c.b), lam, true);
         tipBody = body(WireShape::V3(tip.x, tip.y, tip.z));
         FinishShape(reachM);
         FinishSwitch(s, b);
      }

      // A blade (Inder): bl the posed mesh (PoseBlade), every triangle touchable with the skin edge_mm; the mass from the two flat side faces (|n.x| > 0.9 in the mesh's frame,
      // the blade's thickness along mesh x), density x thick / 2 each. The lever tip: bl.tip; the lane runs from it toward the blade's top
      void Build(const Settings &s, const double surfaceZ, const Blade &bl, const BallSI &b, const V3 &gravity)
      {
         set = s;
         const double mm = 1e-3;
         const std::vector<V3> &tri = bl.tri;
         tip = bl.tip;
         double zLo = INFINITY, zHi = -INFINITY;
         for (const V3 &p : tri)
         {
            zLo = std::min(zLo, p.z);
            zHi = std::max(zHi, p.z);
         }
         bladeH = zHi - zLo;
         set.peakMm = (zHi - surfaceZ) / mm;   // the drawn blade's top above its surface (peak_mm is not used to place a blade)
         // the lane: along the blade (the mesh's y) from its -y end, on the centre plane (not toward the top vertex, which sits on one face: a skewed axis)
         V3 lane(bl.top1.x - bl.top0.x, bl.top1.y - bl.top0.y, 0.0);
         lane = lane / lane.Length();
         Begin(V3(tip.x, tip.y, surfaceZ) - lane * (s.hingeBackMm * mm) - V3(0.0, 0.0, s.hingeBelowMm * mm), lane * -1.0, gravity);
         muPivot = s.muPivot;
         rho = 0.5 * bl.thick;   // the blade's half thickness: pin, stop edge, bending
         rPin = rho;
         halfT = s.edgeMm * mm;
         const auto body = [&](const V3 &p) { const V3 r = p - A; return V3(r.Dot(E), -r.Dot(Z), r.Dot(N)); };
         const V3 side = lane.Cross(V3(0.0, 0.0, 1.0));   // the blade's thickness direction (the mesh's x, rotated)
         double reachM = 0.0;
         peakY = INFINITY;
         peakZ = 0.0;
         for (size_t i = 0; i + 2 < tri.size(); i += 3)
         {
            const V3 a = body(tri[i]), bb = body(tri[i + 1]), c = body(tri[i + 2]);
            for (const V3 &p : { a, bb, c })
            {
               reachM = std::max(reachM, p.Length() + halfT);
               if (p.y - halfT < peakY)
               {
                  peakY = p.y - halfT;
                  peakZ = p.z;
               }
            }
            tris.push_back(RubberContact::MakeSurfaceTri(a, bb, c));
            const V3 n = (tri[i + 1] - tri[i]).Cross(tri[i + 2] - tri[i]);
            const double an = n.Length();
            if (an > 0.0 && std::fabs(n.Dot(side)) / an > 0.9)
               mp.AddPlate(a, bb, c, s.density * bl.thick * 0.5);
         }
         tipBody = body(tip);
         FinishShape(reachM);
         FinishSwitch(set, b);
      }

      // both shapes: the stops, the switch points, the spring and the ball contact law
      void FinishSwitch(const Settings &s, const BallSI &b)
      {
         const double mm = 1e-3;
         // the travel of the peak: rest stop at theta 0, full travel where it has dropped peak + stop_below, the switch's make and break
         q = rate = 0.0;
         thFull = AngleAtDrop((s.peakMm + s.stopBelowMm) * mm);
         switchMake = AngleAtDrop(s.closeMm * mm);
         switchBreak = AngleAtDrop(std::max(s.closeMm - s.releaseMm, 0.0) * mm);
         const auto stop = [&](const double at, const int side) {
            GroupStep::BodyStop st;
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

      // ---- wire_model 1: a massless lever (runs/2026-11-04-massless-wire/design.md). No inertia, no unknown, no pivot friction, no stop impacts: the group solves its angle
      // from all the balls on it at once, the minimum of the spring + weight + their contacts' springs between the rest stop and full travel (Group::EnvSolve), and takes that
      // one energy U* over the balls' moves. Its angle after each step is the one solved for the balls where they end (AfterStep): the switch and the drawing read it.
      bool Massless() const { return set.wireModel == 1 && bladeH == 0.0; }
      bool Envelope() const override { return Massless(); }
      void EnvRange(double &q0, double &lo, double &hi) const override
      {
         q0 = q;
         lo = 0.0;
         hi = thFull;
      }
      double EnvPotential(const double th, double *d1, double *d2) const override
      {
         const double e = 1e-6;   // (the curvature only steers the solve and the step rule; the energy and slope are exact)
         *d1 = -StoredForce(th);
         *d2 = (StoredForce(th - e) - StoredForce(th + e)) / (2.0 * e);
         return StoredEnergy(th);
      }
      int NumDof() const override { return Massless() ? 0 : BodyPart::NumDof(); }
      int MaxDof() const override { return Massless() ? 0 : BodyPart::MaxDof(); }
      int NumCoords() const override { return Massless() ? 0 : BodyPart::NumCoords(); }
      int DofOf(const int c) const override { return Massless() ? -1 : BodyPart::DofOf(c); }
      double KineticJ() const override { return Massless() ? 0.0 : BodyPart::KineticJ(); }
      double Energy(const GroupStep::Ball *b, const int nb) const override { return Massless() ? 0.0 : BodyPart::Energy(b, nb); }   // (massless: all of it is in the group's U*)
      void StepFacts(const GroupStep::Ball *b, const int nb, const GroupStep::Settings &st, GroupStep::StepFactList &f) override
      {
         if (!Massless())
            BodyPart::StepFacts(b, nb, st, f);
      }
      bool BeginStep(const GroupStep::Ball *b, const int nb, const double h, double *x, const int off) override
      {
         hStep = h;
         return Massless() ? true : BodyPart::BeginStep(b, nb, h, x, off);
      }
      bool Eval(const double *x, const int off, const double h, V3 *F, V3 *M, double *Res, GroupStep::Jac J) override { return Massless() ? true : BodyPart::Eval(x, off, h, F, M, Res, J); }
      void Loads(const double *qf, const V3 &F) override
      {
         if (!Massless())
            BodyPart::Loads(qf, F);
      }
      double StopExcess(const double *x, const int off, const unsigned mask = ~0u) const override { return Massless() ? -INFINITY : BodyPart::StopExcess(x, off, mask); }
      unsigned CrossingStops(const double *x, const int off) const override { return Massless() ? 0u : BodyPart::CrossingStops(x, off); }
      bool PinOnStop(const double *x, const int off) override { return Massless() ? false : BodyPart::PinOnStop(x, off); }
      bool FreedByStep() override { return Massless() ? false : BodyPart::FreedByStep(); }
      void Accept(const double *x, const int off, const double h, const bool onStop, double &driveJ, double &lossJ, bool &through) override
      {
         if (!Massless())
         {
            BodyPart::Accept(x, off, h, onStop, driveJ, lossJ, through);
            return;
         }
         driveJ = lossJ = 0.0;
         through = false;
      }
      void AfterStep(const GroupStep::Group &g) override
      {
         if (!Massless())
            return;
         const double qn = g.EnvAngleNow(this);
         const double vBall = g.EnvBallSpeed(this);   // (its tip moves no faster than the balls on it: the speed bound and the rest test; it carries no energy)
         rate = qn == q || reach <= 0.0 ? 0.0 : std::copysign(vBall / reach, qn - q);
         q = qn;
         UpdateSwitches();
      }
      void ExplicitMove(const double *Res, const int off, const double h, double *x) const override
      {
         if (!Massless())
            BodyPart::ExplicitMove(Res, off, h, x);
      }
      double EndResidual() const override { return Massless() ? 0.0 : BodyPart::EndResidual(); }
      std::string ContactNote(const double *x, const int off, const double h, const double *Res) const override
      {
         if (!Massless())
            return HingeGroup::Part::ContactNote(x, off, h, Res);
         char t[120];
         snprintf(t, sizeof t, " massless wire q %.6g deg switch %d", q * 57.29577951308232, (int)switchClosed);
         return t;
      }
      std::string SearchNote(const double *x, const int off, const double h) const override { return Massless() ? std::string(" massless wire") : BodyPart::SearchNote(x, off, h); }
      std::string RowNote() const override { return Massless() ? std::string(" massless wire") : BodyPart::RowNote(); }

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
      // a blade (bladeH > 0): a plate cantilever, I = h t^3 / 12 (h the blade's height, t its thickness)
      double BendKb(const double y) const override
      {
         const double d = 2.0 * rho, I = bladeH > 0.0 ? bladeH * d * d * d / 12.0 : GroupStep::kPi * d * d * d * d / 64.0;
         return 3.0 * set.steelEGpa * 1e9 * I / (y * y * y);
      }
      // the wire's plane holds N, out of it is the axis (x); held at the lever tip
      double BendAlongNormal(const V3 &nb, const V3 &pb) const override
      {
         const double oop = std::fabs(nb.x);
         return oop > 1e-6 ? BendKb(std::max((pb - tipBody).Length(), 1e-4)) / (oop * oop) : 0.0;
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
      double bladeH = 0.0;   // a blade's height (0: a wire)
      double hStep = 0.0;    // the step being taken (massless: its angle's change over it as a rate)
      V3 tip, tipBody;   // the switch leg's end (the lever tip): world, body coordinates
      int stopRest = -1, stopFull = -1;
   };
}
