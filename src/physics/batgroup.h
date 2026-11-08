// license:GPLv3+
//
// The flipper bat as a part of the group step (runs/2026-10-31-bat-in-group/design.md, runs/2026-10-31-flipper-engine/notes.md). Settings, kinematics and force laws are coilmodel.h's; the moving is the group's.
// Modes as ContactBat::Step: rigid (coil off, link taut: one unknown, plunger slaved, inertia I(phi)), two-body (phi and the plunger's x, joined by the link with free play),
// a body pressed on a hard stop has no unknown. Stored energy with exact differences: return spring, rest stop, link. The coil is the only drive; its work is counted.
// The link's spring is C1 and convex in the mismatch, so the step keeps its root. The dampers act on the part of the move inside their zone (continuous in the move), which is what
// removes the force step at contact. Pure SI.
// Coordinates (GroupStep::Part): the bat's angle, the plunger's place, and the pivot's place along x, which has no unknown: its pulse (the other bat's strike) is a given move, a drive, carried inside the step.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstring>

#include "physics/coilmodel.h"
#include "physics/groupstep.h"

namespace BatGroup
{
   using GroupStep::V3;
   using GroupStep::Ball;
   using GroupStep::Hit;
   namespace FM = FlipperModel;

   // The bat's rubber: the surface its balls meet (BatPart::Part::GroupProbe)
   struct BatGeom
   {
      V3 pivotRest;                  // the pivot at rest, m; z = 0: the rubber's heights are from the playfield
      V3 pivot;                      // the pivot now (a nudge moves it along x), m
      double vPivotX = 0.0;          // its speed along x, m/s
      RubberContact::Outline ol;
      double angleStart = 0.0, dirSign = 1.0;
   };

   struct PartBooks
   {
      double coilJ = 0.0, restDampJ = 0.0, linkDampJ = 0.0, frictionJ = 0.0, stopJ = 0.0, captureJ = 0.0;
      long plungerStops = 0, phiStops = 0, captures = 0;
   };

   // The part's values by lifetime (PartRecords): its state and what one attempt keeps (the books are state: a retake books the same)
   struct BatState
   {
      PartBooks bk;
      bool rigid = true, wasPreEos = true, coilOffPrev = true, plungerPinned = false, phiPinned = false, hitStop = false;
      double phi = 0.0, omega = 0.0, xp = 0.0, vp = 0.0, iHigh = 0.0, iLow = 0.0;
      bool hitStopLatch = false;   // the plunger reached the coil stop since the table last read it (EOS, the nudge of the other bat)
      double hitL = 0.0;           // the angular momentum at the last such hit, kg m^2/s
      double aPhi = 0.0, pendingStopJ = 0.0;
      FM::PivotPulse pulse;   // the other bat's strike moves the pivot
   };
   struct BatAttemptWork
   {
      int nDof = 0;
      FM::Kin k0{};
      double gI0 = 1.0, gP0 = 0.0, gV0 = 0.0, gRest0 = 0.0, gCr = 0.0, gGivenX = 0.0, gLoad = 0.0;
      FM::CoilStep gCoil;
      bool gLoop = false, gDeadP = false, gDeadH = false, gCoilLive = false;
      double gX0 = 0.0, gFrac0 = 0.0, gVolts = 0.0, gIP0 = 0.0, gIH0 = 0.0;
      double gDrive = 0.0, gLossRest = 0.0, gLossLink = 0.0, gLossB = 0.0;
      bool gPreEos = true;
   };

   class Part : public GroupStep::Part, public BatState, protected BatAttemptWork
   {
   public:
      Part()
      {
         KeepState<BatState>(this);
         WipeEachAttempt<BatAttemptWork>(this);
      }
      FlipperSettings s;
      FlipperDerived p;
      BatGeom geo;            // (its pivot: made from the pulse, MovePivot)
      RubberContact::RubberParams rp;         // the bat's rubber
      double layerM = 0.0;
      bool solenoid = false;  // the table's coil command (between cycles)
      double reach = 0.1;     // m, the lever that turns the bat's rate into a speed
      static constexpr int kPhi = 0, kX = 1, kPivot = 2;   // the coordinates

      // the settings (a change while running keeps the state; the rubber is the caller's: FlipperModel::MakeRubberParams)
      void Configure(const FlipperSettings &st, const RubberContact::RubberParams &rubber)
      {
         s = st;
         rp = rubber;
         layerM = s.rubberThicknessMm * 1e-3;
         FM::Derive(s, p);
         coilP.vSupply = coilH.vSupply = s.vSupply;
         coilP.bobbinMm = coilH.bobbinMm = s.bobbinMm;
         coilP.rOhm = s.rP;
         coilP.lAirMH = s.lAirP;
         coilP.lArmMH = s.lArmP;
         coilP.x0 = s.x0P;
         coilH.rOhm = s.rH;
         coilH.lAirMH = s.lAirH;
         coilH.lArmMH = s.lArmH;
         coilH.x0 = s.x0H;
         CoilModel::Derive(coilP, derP);
         CoilModel::Derive(coilH, derH);
         kLink = s.kContactNPerMm * 1000.0;
         cStiff = p.cStiffNsPerMm * 1000.0;
         cFree = s.freePlayDamping ? p.cContactNsPerMm * 1000.0 : 0.0;
         slop = p.slopMm * 1e-3;
         kSpring = p.kSpringNPerMm * 1000.0;
         xRest = p.xRest * 1e-3;
         xMin = p.xMin * 1e-3;
         xMax = p.xMax * 1e-3;
         bobbin = s.bobbinMm * 1e-3;
         phiStop = s.stopDeg * FM::DEG;
         mp = s.plungerKg;
         If = p.iFlipper;
         reach = geo.ol.Lc + geo.ol.re;
      }
      // The bat at angle phi0 (rad from rest), still, plunger slaved: one rigid body with the coil off (as ContactBat::Reset)
      void Reset(const double phi0, const bool sol)
      {
         phi = phi0;
         omega = 0.0;
         vp = 0.0;
         xp = X(phi0);
         iHigh = iLow = 0.0;
         rigid = !sol;
         wasPreEos = coilOffPrev = true;
         plungerPinned = phiPinned = false;
         solenoid = sol;
         aPhi = pendingStopJ = 0.0;
      }
      // the table's coil, between cycles. A coil that fires ends the rigid mode: the plunger is slaved, so the mismatch is 0 and nothing jumps
      void SetCoil(const bool on)
      {
         solenoid = on;
         if (on && rigid)
         {
            rigid = false;
            phiPinned = false;
            plungerPinned = false;
            const FM::Kin k = Kin(phi);
            xp = k.x * 1e-3;
            vp = -k.lever * 1e-3 * omega;
            PinCheck(false, false);
         }
         else if (!on && !rigid)   // the coil let go with the link taut and slow: one body from here, not a step later: its windings' flux is at its highest now, and a pinned plunger takes none of it
            bk.captureJ += Capture();
      }
      // the pivot now, from the pulse
      void MovePivot()
      {
         double d, v;
         pulse.Eval(pulse.t, d, v);
         geo.pivot = V3(geo.pivotRest.x + d, geo.pivotRest.y, geo.pivotRest.z);
         geo.vPivotX = v;
      }
      // The rubber for a ball centre and radius, the coordinates moved by trial (null = where they are): depth, normal, radius, and how its point moves with each coordinate
      void ProbeRubber(const V3 &c, const double R, const double *const trial, Hit &h) const
      {
         const V3 cRel = c - V3(geo.pivot.x + (trial != nullptr ? trial[kPivot] : 0.0), geo.pivot.y, geo.pivot.z);
         RubberContact::BallSI b;
         b.R = R;
         const RubberContact::Probe pr = RubberContact::ProbeSurface(rp, b, geo.ol, geo.angleStart + geo.dirSign * (phi + (trial != nullptr ? trial[kPhi] : 0.0)) - FM::PI / 2.0, cRel);
         h.delta = pr.delta;
         h.n = pr.dist > 1e-12 ? (cRel - pr.point) * (1.0 / pr.dist) : V3(0.0, 0.0, 1.0);
         h.reff = pr.reff;
         h.curvDist = 0.0;
         h.grad = 1.0;
         const V3 wq(0.0, 0.0, geo.dirSign), vq = wq.Cross(V3(pr.point.x, pr.point.y, 0.0));
         h.nc = 2;
         h.c[0].Set(kPhi, h.n.Dot(vq), vq, wq);
         h.c[1].Set(kPivot, h.n.x, V3(1.0, 0.0, 0.0), V3());
      }

      // ---- the group's interface ----
      int NumDof() const override { return rigid ? (phiPinned ? 0 : 1) : (plungerPinned ? 1 : 2); }
      int MaxDof() const override { return 2; }
      int NumCoords() const override { return 3; }
      int DofOf(const int c) const override { return c == kPhi ? (NumDof() >= 1 ? 0 : -1) : c == kX && NumDof() == 2 ? 1 : -1; }
      GroupStep::CoordInfo Coord(const int c) const override
      {
         GroupStep::CoordInfo q;
         if (c == kPhi)
         {
            q.rate = omega;
            q.invM = 1.0 / InertiaNow();
            q.invMFree = q.invM;
            q.reach = reach;
            q.maxMove = 0.05;
            q.blockedSide = rigid && phiPinned ? 1 : 0;
         }
         else if (c == kX)
         {
            q.rate = vp;
            q.invM = 1.0 / mp;
            q.reach = 1.0;
            q.maxMove = 2e-3;
         }
         else
         {
            q.rate = geo.vPivotX;
            q.given = gGivenX;
            q.reach = 1.0;
         }
         return q;
      }
      double Energy(const Ball *, int) const override
      {
         return KineticJ() + Vbat(phi) + (rigid ? 0.0 : VLink(xp - X(phi)));
      }
      double KineticJ() const override
      {
         return rigid ? 0.5 * Inertia(phi) * omega * omega : 0.5 * If * omega * omega + 0.5 * mp * vp * vp;
      }
      // The step's facts, from the state now (accuracy only): the group holds the factors
      void StepFacts(const Ball *, int, const GroupStep::Settings &gs, GroupStep::StepFactList &f) override
      {
         const FM::Kin k = Kin(phi);
         const double lever = k.lever * 1e-3;
         const double I = rigid ? If + mp * lever * lever : If;
         const double m = rigid ? 0.0 : xp - k.x * 1e-3, rate = rigid ? 0.0 : vp + lever * omega;
         const bool still = std::fabs(omega) * reach <= gs.vSteady && std::fabs(vp) <= gs.vSteady;
         if (!phiPinned)
         {
            const double omR = std::sqrt(s.restKNmPerRad / I);
            if (phi <= 0.0)
               f.Stiffness(omR, still, 1);
            else if (omega < 0.0)
               f.Approach(phi, -omega, omR, 2);
            if (std::fabs(omega) > 1e-9)
               f.Turn(1.0, std::fabs(omega), 5);
         }
         if (!rigid)
         {
            const double invMu = lever * lever / If + (plungerPinned ? 0.0 : 1.0 / mp);
            const double omL = std::sqrt(kLink * invMu);
            if (m < 0.0 || m > slop)
               f.Stiffness(omL, still, 3);
            else if (std::fabs(rate) > 1e-9)
               f.Approach(rate > 0.0 ? slop - m : m, std::fabs(rate), omL, 4);
         }
         if (solenoid && phi < p.eosDeg * FM::DEG && omega > 1e-3)
            f.Arrival(p.eosDeg * FM::DEG - phi, omega, 0.0, 6);
         // the coil's own time
         const double vX = rigid ? -lever * omega : vp;
         const double tau = FM::CoilTimeConstant(s, p, xp, vX, iHigh, iLow, solenoid, phi < p.eosDeg * FM::DEG);
         if (std::isfinite(tau))
            f.TimeConstant(tau, 7);
         if (rigid ? !phiPinned : !plungerPinned)   // steep field: a winding that holds force. The plunger moves at most a share of its distance from the field's centre
         {
            const double frac = (xp - xMin) / bobbin;
            const double fP = p.strengthP * FM::Sq(s.x0P / (frac + s.x0P)) * iHigh * iHigh, fL = p.strengthH * FM::Sq(s.x0H / (frac + s.x0H)) * iLow * iLow;
            const double acc = (fP + fL) / (rigid ? I / std::max(lever * lever, 1e-12) : mp), vc = std::fabs(vX);
            if (fP > 0.0)
               f.Field((frac + s.x0P) * bobbin, vc, acc, 8);
            if (fL > 0.0)
               f.Field((frac + s.x0H) * bobbin, vc, acc, 9);
         }
         if (pulse.active)   // the pulse's turn and end
            f.TimeLeft(pulse.t < pulse.rise ? pulse.rise - pulse.t : pulse.rise + pulse.fall - pulse.t, 10);
      }
      bool BeginStep(const Ball *, int, const double h, double *x, const int off) override
      {
         nDof = NumDof();
         gPreEos = phi < p.eosDeg * FM::DEG;
         gGivenX = pulse.Move(h);
         gLoad = 0.0;
         k0 = Kin(phi);
         BeginCoil();
         gI0 = rigid ? If + mp * Sq(k0.lever * 1e-3) : If;
         gP0 = gI0 * omega;
         gV0 = Vbat(phi);
         gRest0 = VRest(phi);
         gCr = RestC(gI0);
         if (nDof >= 1)
            x[off] = omega * h;
         if (nDof == 2)
            x[off + 1] = vp * h;
         return true;
      }
      void Loads(const double *q, const V3 &) override { gLoad = q[kPhi]; }
      bool Eval(const double *x, const int off, const double h, V3 *, V3 *, double *Res, GroupStep::Jac J) override
      {
         gLossRest = gLossLink = gLossB = gDrive = 0.0;
         if (nDof == 0)
         {
            TakeCoil(CoilOver(0.0, h));
            nM1 = rigid ? 0.0 : xp - k0.x * 1e-3;
            nLink = rigid ? 'r' : 'h';
            nFrac1 = gFrac0;
            return true;
         }
         const double dphi = x[off], dx = nDof == 2 ? x[off + 1] : 0.0;
         const double phi1 = phi + dphi, phm = phi + 0.5 * dphi;
         const FM::Kin k1 = Kin(phi1), km = Kin(phm);
         const double l1 = k1.lever * 1e-3, lm = km.lever * 1e-3;
         const GroupStep::MeanInertia mi(gI0, rigid ? If + mp * l1 * l1 : If, gP0, dphi, h);   // (the shared inertia row: the bat's inertia varies with its angle)
         double dT = 0.0;
         if (rigid)
         {
            if (std::fabs(dphi) > 1e-9)
               dT = mi.KineticSlope((1.0 / mi.I1 - 1.0 / gI0) / dphi);
            else
            {
               // d(1/I)/dphi = -I'/I^2, I' = 2 mp l dl (the lever's slope along phi is k.dLever)
               const double Im = If + mp * lm * lm, Ip = 2.0 * mp * lm * km.dLever * 1e-3;
               dT = mi.KineticSlope(-Ip / (Im * Im));
            }
         }
         // spring and rest stop: exact differences
         const double V1 = Vbat(phi1), Vr1 = VRest(phi1);
         const double Vp = Chord(VbatP(phm, km), V1, gV0, dphi);
         const double tauSr = -Chord(VRestP(phm), Vr1, gRest0, dphi);
         // the rest stop's damper on the part of the move inside the stop; the stop pushes, it cannot pull (the shared contact's N = max(0, fe + fd))
         const double pen0 = std::max(0.0, -phi), pen1 = std::max(0.0, -phi1);
         const double tauN = std::max(0.0, tauSr + gCr * (pen1 - pen0) / h), tauD = tauN - tauSr;
         const bool restClamped = tauN <= 0.0;
         const double om = dphi / h;
         double Qphi = -Vp + tauD - s.bDamping * om - dT;
         gLossRest = -tauD * dphi;
         gLossB = s.bDamping * om * dphi;
         double jpp = mi.Slope(h) + s.bDamping;
         if (J)
         {
            double kc;
            if (std::fabs(dphi) > 1e-9)
               kc = (VbatP(phi1, k1) - Vp) / dphi;
            else
               kc = 0.5 * (kSpring * lm * lm + (phm < 0.0 ? s.restKNmPerRad : 0.0));
            jpp += h * std::max(kc, 0.0);
            if (!restClamped && phi1 < 0.0)
               jpp += gCr;
         }
         // the coil, every wiring: the circuit and its force over the solved move, in the solve (CoilModel), so what the coil gives is what the circuit supplied at any step length; rigid, through the link's chord, so the bat is given exactly the plunger's work
         double dF = 0.0;
         double Qx = 0.0;
         {
            const double dClose = rigid ? CloseOfBat(dphi, k1, lm) : -dx;
            nFrac1 = (gX0 - dClose - xMin) / bobbin;
            const CoilMove m = CoilOver(dClose, h);
            TakeCoil(m);
            if (rigid)
            {
               const double Tc = BatTorque(m.force, dphi, dClose, lm);
               Qphi += Tc;
               gDrive = Tc * dphi;
               if (J)
               {
                  const double e = 1e-9 + 1e-6 * std::fabs(dphi);
                  const auto TorqueAt = [&](const double d) { const FM::Kin kd = Kin(phi + d); const double ld = Kin(phi + 0.5 * d).lever * 1e-3, c = CloseOfBat(d, kd, ld); return BatTorque(CoilOver(c, h).force, d, c, ld); };
                  jpp -= h * (TorqueAt(dphi + e) - TorqueAt(dphi - e)) / (2.0 * e);
               }
            }
            else
            {
               Qx = -m.force;
               gDrive = m.force * dClose;
               if (J)
               {
                  const double e = 1e-9 + 1e-6 * std::fabs(dx);
                  dF = (CoilOver(-dx - e, h).force - CoilOver(-dx + e, h).force) / (2.0 * e);   // along the plunger's x, for the row below
               }
            }
         }
         double jxx = GroupStep::InertiaSlope(mp, h), jpx = 0.0;
         if (!rigid)
         {
            // the link: the mismatch's discrete gradient over (phi, x)
            const double m0 = xp - k0.x * 1e-3, m1 = xp + dx - k1.x * 1e-3, dm = m1 - m0;
            const double gP = lm, gX = nDof == 2 ? 1.0 : 0.0;
            const double W0 = 1.0 / If, W1 = 1.0 / mp;
            const double gdq = gP * dphi + gX * dx, qq = dphi * W0 * dphi + (nDof == 2 ? dx * W1 * dx : 0.0);
            const double cf = qq > 1e-300 && std::fabs(dm - gdq) > 16.0 * 2.220446049250313e-16 * 0.1 ? (dm - gdq) / qq : 0.0;   // (a mismatch within the rounding of m is left out, as Group::Gonzalez does: over a very small move it would be a large false slope)
            const double gdP = gP + cf * W0 * dphi, gdX = gX + (nDof == 2 ? cf * W1 * dx : 0.0);
            const double Fe = LinkChord(m0, m1);
            const double dGc = cStiff * (std::min(m1, 0.0) - std::min(m0, 0.0) + std::max(m1 - slop, 0.0) - std::max(m0 - slop, 0.0));
            const double dGf = cFree * (Clamp(m1) - Clamp(m0));
            double Tc = -Fe - dGc / h;
            bool clamped = false;
            if (Fe < 0.0 && Tc < 0.0)
            {
               Tc = 0.0;
               clamped = true;
            }
            else if (Fe > 0.0 && Tc > 0.0)
            {
               Tc = 0.0;
               clamped = true;
            }
            else if (Fe == 0.0)
               Tc = 0.0;
            const double T = Tc - dGf / h;
            Qphi += T * gdP;
            Qx += T * gdX;
            gLossLink = -(T + Fe) * dm;
            nM1 = m1;
            nLink = clamped ? 'c' : (m1 < 0.0 ? 'p' : (m1 > slop ? 'f' : 's'));
            if (J)
            {
               const bool in = m1 < 0.0 || m1 > slop;
               const double kc = std::fabs(dm) > 1e-13 ? (VLinkP(m1) - Fe) / dm : (in ? 0.5 * kLink : 0.0);
               double kn = 0.0;
               if (!clamped && Fe != 0.0)
                  kn += h * std::max(kc, 0.0) + (in ? cStiff : 0.0);
               if (!in)
                  kn += cFree;
               jpp += kn * gP * gP;
               jpx += kn * gP * gX;
               jxx += kn * gX * gX + h * dF;   // (the magnet's pull: a negative stiffness)
            }
         }
         if (rigid)
         {
            nM1 = 0.0;
            nLink = 'r';
         }
         Res[off] = mi.Row() - h * Qphi;
         if (nDof == 2)
            Res[off + 1] = GroupStep::InertiaRow(dx, vp, mp, h) - h * Qx;
         if (J)
         {
            J[off][off] += jpp;
            if (nDof == 2)
            {
               J[off][off + 1] += jpx;
               J[off + 1][off] += jpx;
               J[off + 1][off + 1] += jxx;
            }
         }
         return true;
      }
      // stops: bit 0 the bat's up-stop (rigid), bit 1 the plunger's closed end, bit 2 its open end (two-body)
      double StopExcess(const double *x, const int off, const unsigned mask = ~0u) const override
      {
         if (nDof == 0)
            return -INFINITY;
         if (rigid)
            return (mask & 1u) ? phi + x[off] - phiStop : -INFINITY;
         if (nDof == 2)
         {
            // (a plunger slaved to a bat pressed below its rest angle starts past the open end: its stop is where it is)
            const double x1 = xp + x[off + 1];
            double e = -INFINITY;
            if (mask & 2u)
               e = std::max(e, std::min(xMin, xp) - x1);
            if (mask & 4u)
               e = std::max(e, x1 - std::max(xMax, xp));
            return e;
         }
         return -INFINITY;
      }
      unsigned CrossingStops(const double *x, const int off) const override
      {
         unsigned m = 0u;
         for (unsigned b = 1u; b <= 4u; b <<= 1)
            if (StopExcess(x, off, b) > 0.0)
               m |= b;
         return m;
      }
      // The body is on its stop (the group's test) and the free step presses it in: it stays on the stop for this step (a contact's complementarity, decided per step: no time is lost and nothing chatters); inelastic, as the stop's law
      bool PinOnStop(const double *, int) override
      {
         if (rigid)
         {
            if (phiPinned || omega < 0.0)
               return false;
            pendingStopJ += 0.5 * Inertia(phi) * omega * omega;
            omega = 0.0;
            phiPinned = true;
            return true;
         }
         if (plungerPinned || (xp < 0.5 * (xMin + xMax) ? vp > 0.0 : vp < 0.0))
            return false;
         pendingStopJ += 0.5 * mp * vp * vp;
         vp = 0.0;
         plungerPinned = true;
         return true;
      }
      void Accept(const double *x, const int off, const double h, const bool onStop, double &driveJ, double &lossJ, bool &through) override
      {
         double stopJ = pendingStopJ, capJ = 0.0;
         pendingStopJ = 0.0;
         hitStop = false;
         if (nDof >= 1)
         {
            const double dphi = x[off], phi1 = phi + dphi;
            const GroupStep::MeanInertia mi(gI0, rigid ? Inertia(phi1) : If, gP0, dphi, h);
            const double I1 = mi.I1;
            double w1 = mi.EndRate();
            aPhi = (w1 - omega) / h;
            if (rigid)
            {
               if (onStop)
               {
                  stopJ += 0.5 * I1 * w1 * w1;
                  hitStop = w1 > 0.0;
                  if (hitStop)
                     hitL = I1 * w1;
                  w1 = 0.0;
                  phiPinned = true;
                  bk.phiStops++;
               }
               phi = phi1;
               omega = w1;
               const FM::Kin k = Kin(phi);
               xp = k.x * 1e-3;
               vp = -k.lever * 1e-3 * omega;
            }
            else
            {
               phi = phi1;
               omega = w1;
               if (nDof == 2)
               {
                  double v1 = GroupStep::EndRate(x[off + 1], vp, h);
                  xp += x[off + 1];
                  if (onStop)
                  {
                     stopJ += 0.5 * mp * v1 * v1;
                     v1 = 0.0;
                     plungerPinned = true;
                     bk.plungerStops++;
                     hitStop = xp < 0.5 * (xMin + xMax);
                     if (hitStop)
                        hitL = (If + mp * Sq(Kin(phi).lever * 1e-3)) * std::fabs(omega);
                  }
                  vp = v1;
               }
            }
         }
         hitStopLatch = hitStopLatch || hitStop;
         iHigh = gCoil.iHigh;
         iLow = gCoil.iLow;
         wasPreEos = gPreEos;
         coilOffPrev = !solenoid;
         PinCheck(onStop, true);
         stopJ += pendingStopJ;
         pendingStopJ = 0.0;
         capJ = Capture();
         bk.coilJ += gDrive;
         bk.restDampJ += gLossRest + gLossB;
         bk.linkDampJ += gLossLink;
         bk.stopJ += stopJ;
         bk.captureJ += capJ;
         driveJ = gDrive;
         lossJ = gLossRest + gLossB + gLossLink + stopJ + capJ;
         through = false;
         if (pulse.active)
         {
            pulse.Advance(h);
            MovePivot();
         }
      }
      void AfterRestore() override { MovePivot(); }
      std::string StartNote() const override
      {
         BatState sv;
         std::memcpy(&sv, SavedState(), sizeof sv);   // (its state record)
         char s[400];
         snprintf(s, sizeof s, "rigid %d preEos %d coilOffPrev %d plungerPinned %d phiPinned %d solenoid %d phi %.17g omega %.17g xp %.17g vp %.17g iHigh %.17g iLow %.17g pulse %d t %.17g amp %.17g",
            sv.rigid, sv.wasPreEos, sv.coilOffPrev, sv.plungerPinned, sv.phiPinned, solenoid, sv.phi, sv.omega, sv.xp, sv.vp, sv.iHigh, sv.iLow, sv.pulse.active, sv.pulse.t, sv.pulse.amp);
         return s;
      }
      std::string StopNote() const override
      {
         char s[260];
         snprintf(s, sizeof s, "rigid %d phiPinned %d plungerPinned %d solenoid %d dof %d phi %.9g omega %.6g xp %.9g vp %.6g xMin %.9g xMax %.9g", rigid, phiPinned, plungerPinned, solenoid, nDof, phi, omega, xp, vp, xMin, xMax);
         return s;
      }
      void ScaleSpeeds(const double f) override { omega *= f; vp *= f; }
      void ExplicitMove(const double *Res, const int off, const double h, double *x) const override
      {
         if (nDof >= 1)
            x[off] -= Res[off] * h / (2.0 * gI0);
         if (nDof == 2)
            x[off + 1] -= Res[off + 1] * h / (2.0 * mp);
      }

      void KeepFailNote(const int slot, const double *Res, const int off) override
      {
         if (slot < 0 || slot >= GroupStep::Group::kHalveDiag)
            return;
         FailNoteData &d = failNote[slot];
         d.nDof = nDof;
         d.resPhi = nDof >= 1 ? Res[off] : 0.0;
         d.resX = nDof == 2 ? Res[off + 1] : 0.0;
         d.m1 = nM1;
         d.link = nLink;
         d.iP = gCoil.iHigh * (s.vSupply / s.rP);
         d.iH = gCoil.iLow * (s.vSupply / s.rH);
         d.frac = nFrac1;
      }
      std::string FailNote(const int slot) const override
      {
         if (slot < 0 || slot >= GroupStep::Group::kHalveDiag)
            return std::string();
         const FailNoteData &d = failNote[slot];
         char t[200];
         snprintf(t, sizeof t, "dof %d res phi %.3g N m s x %.3g N s, link %c m %.4g mm, coil %.3g A / %.3g A, frac %.5f", d.nDof, d.resPhi, d.resX, d.link, d.m1 * 1e3, d.iP, d.iH, d.frac);
         return t;
      }

      double InertiaNow() const { return rigid ? Inertia(phi) : If; }   // what the balls feel
      double InertiaRead() const { return rigid || xp - X(phi) <= 0.0 ? Inertia(phi) : If; }   // what a script reads (Mass): with the plunger's reflected mass while the link is loaded, as before
      double Inertia(const double ph) const { return If + mp * Sq(Kin(ph).lever * 1e-3); }
      double X(const double ph) const { return FM::SliderX(s, p, p.theta0 + ph) * 1e-3; }
      double PhiStop() const { return phiStop; }

   private:
      static double Sq(const double a) { return a * a; }
      FM::Kin Kin(const double ph) const { return FM::Kinematics(s, p, p.theta0 + ph); }
      double VSpring(const double x) const
      {
         const double u = xRest - x;
         return s.preloadN * u + 0.5 * kSpring * u * u;
      }
      double VRest(const double ph) const
      {
         const double pen = -ph;
         if (pen <= 0.0)
            return 0.0;
         const double penRef = 1.4 * FM::PI / 180.0;
         return s.restKNmPerRad * penRef * penRef / (s.restN + 1.0) * std::pow(pen / penRef, s.restN + 1.0);
      }
      double VRestP(const double ph) const
      {
         const double pen = -ph;
         if (pen <= 0.0)
            return 0.0;
         const double penRef = 1.4 * FM::PI / 180.0;
         return -s.restKNmPerRad * penRef * std::pow(pen / penRef, s.restN);
      }
      double RestC(const double I) const { return 2.0 * s.restZeta * std::sqrt(s.restKNmPerRad * I); }
      double Vbat(const double ph) const { return VSpring(X(ph)) + VRest(ph); }
      double VbatP(const double ph, const FM::Kin &k) const
      {
         return (s.preloadN + kSpring * (xRest - k.x * 1e-3)) * k.lever * 1e-3 + VRestP(ph);
      }
      double Clamp(const double m) const { return std::min(std::max(m, 0.0), slop); }
      double VLink(const double m) const
      {
         const double o = m < 0.0 ? m : (m > slop ? m - slop : 0.0);
         return 0.5 * kLink * o * o;
      }
      double VLinkP(const double m) const { return kLink * (m < 0.0 ? m : (m > slop ? m - slop : 0.0)); }
      // (V1 - V0) / move, or the slope at mid-step over a move so small that the two differ only by rounding (V near 0.2 J rounds to 1e-17 J; the slope's error is V''' move^2 / 24, under 1e-11 of a N m below 1e-5 rad, which is what the rounding is above it)
      static double Chord(const double slopeMid, const double V1, const double V0, const double move) { return std::fabs(move) > 1e-5 ? (V1 - V0) / move : slopeMid; }
      // the link's chord, with no difference of energies inside one zone of the mismatch
      double LinkChord(const double m0, const double m1) const
      {
         const int z0 = m0 < 0.0 ? -1 : (m0 > slop ? 1 : 0), z1 = m1 < 0.0 ? -1 : (m1 > slop ? 1 : 0);
         if (z0 == z1)
            return z0 == 0 ? 0.0 : kLink * (0.5 * (m0 + m1) - (z0 > 0 ? slop : 0.0));
         return (VLink(m1) - VLink(m0)) / (m1 - m0);
      }
      double CoilForce(const double frac, const double i2P, const double i2H, double &dFdx) const
      {
         const double aP = frac + s.x0P, aH = frac + s.x0H;
         const double gP = FM::Sq(s.x0P / aP), gH = FM::Sq(s.x0H / aH);
         const double fP = p.strengthP * gP * i2P, fH = p.strengthH * gH * i2H;
         dFdx = (-2.0 * fP / aP - 2.0 * fH / aH) / bobbin;
         return fP + fH;
      }
      // Every wiring: one law, once per current loop, with no coupling between the windings.
      // Series: before EOS the power winding on the supply and the hold winding in its short are two loops; after EOS, or with the coil off, they are one loop and L is the two windings' sum.
      // Parallel: each winding is its own loop, on the supply while the coil is on and in the flyback diode's loop when off. The EOS switch is in the power winding's loop: open, that loop is gone, so its current is 0 and its field goes to the switch.
      // Single: the power winding alone, on the supply or in its diode's loop.
      // The switches' rule on the currents, at the step's start.
      void BeginCoil()
      {
         gX0 = rigid ? k0.x * 1e-3 : xp;
         gFrac0 = std::max((gX0 - xMin) / bobbin, FracFloor());
         double iP = iHigh * (s.vSupply / s.rP), iH = iLow * (s.vSupply / s.rH);
         const double lP = (p.aP + p.bP / (gFrac0 + s.x0P)) * 1e-3, lH = (p.aH + p.bH / (gFrac0 + s.x0H)) * 1e-3;
         gVolts = solenoid ? s.vSupply : 0.0;
         gLoop = s.series && (!solenoid || !gPreEos);
         gDeadP = gDeadH = false;
         if (s.single)
         {
            iH = 0.0;
            gDeadH = true;
         }
         else if (!s.series)
         {
            if (!gPreEos)
            {
               iP = 0.0;
               gDeadP = true;
            }
         }
         else if (!solenoid)
         {
            if (!coilOffPrev)
               iP = iH = (gPreEos && s.eosFlux) ? (lP * iP + lH * iH) / (lP + lH) : iP;
            else
               iH = iP;
         }
         else if (!gPreEos)
         {
            if (wasPreEos)
               iP = iH = s.eosFlux ? (lP * iP + lH * iH) / (lP + lH) : iP;
            else
               iH = iP;
         }
         gIP0 = iP;
         gIH0 = iH;
         gCoilLive = gVolts != 0.0 || iP != 0.0 || iH != 0.0;
      }
      struct CoilMove { CoilModel::OverMove a, b; double force = 0.0; };   // a: the power winding or the series loop; b: the hold winding in its own loop. force: over the move as asked, so force * dClose is the work
      // The windings over a trial move dClose toward the closed end, from the step's start, so a retake is identical. L's curve runs on a little past the closed end, so a trial there sees no kink (a step ends past the stop only by the stop's tolerance). Far past it, where the curve would blow up, a trial does no more work.
      CoilMove CoilOver(const double dClose, const double h) const
      {
         CoilMove m;
         if (!gCoilLive)
            return m;
         const double fTo = (gX0 - dClose - xMin) / bobbin, f1 = std::max(fTo, FracFloor());
         const bool past = fTo < FracFloor();
         const double dx = past ? (gFrac0 - f1) * bobbin : dClose;
         const double cP = CoilModel::ChordHPerM(coilP, derP, gFrac0, f1), cH = CoilModel::ChordHPerM(coilH, derH, gFrac0, f1);
         if (gLoop)
         {
            m.a = CoilModel::OverLoop(s.rP + s.rH, cP + cH, CoilModel::MeanH(coilP, derP, gFrac0, f1) + CoilModel::MeanH(coilH, derH, gFrac0, f1), dx, gVolts, gIP0, h);
            m.force = m.a.force;
         }
         else
         {
            if (!gDeadP)
               m.a = CoilModel::OverLoop(s.rP, cP, CoilModel::MeanH(coilP, derP, gFrac0, f1), dx, gVolts, gIP0, h);
            if (!gDeadH)
               m.b = CoilModel::OverLoop(s.rH, cH, CoilModel::MeanH(coilH, derH, gFrac0, f1), dx, s.series ? 0.0 : gVolts, gIH0, h);   // (series: the hold winding in the EOS switch's short)
            m.force = m.a.force + m.b.force;
         }
         if (past && dClose != 0.0)
            m.force *= dx / dClose;
         return m;
      }
      void TakeCoil(const CoilMove &m)
      {
         gCoil.iHigh = m.a.i1 / (s.vSupply / s.rP);
         gCoil.iLow = (gLoop ? m.a.i1 : m.b.i1) / (s.vSupply / s.rH);
      }
      double FracFloor() const { return -0.5 * (s.single ? s.x0P : std::min(s.x0P, s.x0H)); }   // (single: the power winding alone; there is no hold winding's curve to stop short of)
      // The plunger's move toward the closed end for a bat's move: the link's chord; over a move so small that it would be rounding, the lever
      double CloseOfBat(const double dphi, const FM::Kin &k1, const double lm) const { return std::fabs(dphi) > 1e-5 ? (k0.x - k1.x) * 1e-3 : lm * dphi; }
      static double BatTorque(const double f, const double dphi, const double dClose, const double lm) { return std::fabs(dphi) > 1e-5 ? f * dClose / dphi : f * lm; }
      // The coil's force on the plunger at its place x (m), toward the closed end, at the currents now (N)
      double CoilForceNow(const double x) const
      {
         double dF;
         return CoilForce(std::max((x - xMin) / bobbin, 0.0), iHigh * iHigh, iLow * iLow, dF);
      }
      // A body on a hard stop stays while the whole load at the end of the step presses into it (the stop pushes, never pulls): spring, ball and the coil (one body: through the link and crank)
      void PinCheck(const bool onStop, const bool afterStep)
      {
         (void)onStop;
         if (rigid)
         {
            if (phiPinned || phi >= phiStop)
            {
               const FM::Kin k = Kin(phi);
               const double net = -VbatP(phi, k) + CoilForceNow(k.x * 1e-3) * k.lever * 1e-3 + (afterStep ? gLoad : 0.0);
               phiPinned = net > 0.0;
            }
            return;
         }
         phiPinned = false;
         const double Fc = CoilForceNow(xp);
         const FM::Kin k = Kin(phi);
         const double m = xp - k.x * 1e-3, rate = vp + k.lever * 1e-3 * omega;
         double T;
         if (m < 0.0)
            T = std::max(0.0, -kLink * m - cStiff * rate);
         else if (m > slop)
            T = std::min(0.0, -kLink * (m - slop) - cStiff * rate);
         else
            T = -cFree * rate;
         const double net = -Fc + T;
         const double tol = GroupStep::Group::kOnStop;
         if ((xp <= xMin + tol && net <= 0.0 && vp <= 1e-9) || (xp >= xMax - tol && net >= 0.0 && vp >= -1e-9))   // a body moving off its stop is not held by it
         {
            plungerPinned = true;
            if (afterStep)
               pendingStopJ += 0.5 * mp * vp * vp;
            else
               bk.stopJ += 0.5 * mp * vp * vp;   // (the coil's switch, between cycles: no step's books)
            vp = 0.0;
         }
         else
            plungerPinned = false;
      }
      // Two-body to rigid (coil off, link taut and slow): an inelastic constraint. Generalized momentum is kept, so the kinetic energy can only fall; the link's held energy goes with it. Returns the loss.
      double Capture()
      {
         if (rigid || solenoid)
            return 0.0;
         const FM::Kin k = Kin(phi);
         const double l = k.lever * 1e-3, m = xp - k.x * 1e-3, rate = vp + l * omega;
         if (!(m < 0.0 && std::fabs(rate) * 1000.0 < s.rigidMaxRelSpeedMmS))
            return 0.0;
         const double E0 = 0.5 * If * omega * omega + 0.5 * mp * vp * vp + VLink(m);
         const double I = If + mp * l * l;
         omega = (If * omega - mp * l * vp) / I;
         vp = -l * omega;
         xp = k.x * 1e-3;
         rigid = true;
         plungerPinned = false;
         phiPinned = false;
         bk.captures++;
         return E0 - 0.5 * I * omega * omega;
      }

      double kLink = 0.0, cStiff = 0.0, cFree = 0.0, slop = 0.0, kSpring = 0.0, xRest = 0.0, xMin = 0.0, xMax = 0.0, bobbin = 0.0, phiStop = 0.0, mp = 0.0, If = 0.0;
      CoilSettings coilP, coilH;   // the windings as the coil's one energy law takes them
      CoilDerived derP, derH;
      // the last evaluation, for the failure notes: the link's mismatch at the move's end (m), its state (r rigid, h plunger pinned, s slack, p pushing,
      // f pressed at the slop's far end, c clamped: its damper would pull, so its force is 0) and the coil's place at the move's end (frac of the bobbin)
      double nM1 = 0.0, nFrac1 = 0.0;
      char nLink = '?';
      struct FailNoteData { int nDof = 0; double resPhi = 0.0, resX = 0.0, m1 = 0.0, iP = 0.0, iH = 0.0, frac = 0.0; char link = '?'; };
      FailNoteData failNote[GroupStep::Group::kHalveDiag];
   };
}
