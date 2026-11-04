// license:GPLv3+

// The fixed-height rubber band on the group step. The cord stays in the band's plane at its own height, no points of its own: its length is a closed form of the balls' positions and its
// stored energy one function of them. Bench: runs/2026-11-01-cord-fixed (design_two.md, the laws, what breaks), runs/2026-11-01-cord-100us (design_fixed.md, design_plate.md).
//
// V = sum over spans [ 1/2 (EA / l0) (L - l0)^2 + the squash at each body on a run ] + the grip springs.
// L = the taut length round the balls (tangents and arcs, in the band's plane) + what the runs gain bending round bodies.
// A body on a run is a contact force: at its station the run is behind it by a gap, shared between lifting the run and squashing the cord where the two forces balance; that sum is one function of the positions.
// A span's ends are notches: fixed points. A ball on a post or the plate is the group's own Item with the ring's law.

#pragma once
#include <algorithm>
#include <cmath>
#include <vector>
#include <memory>
#include <string>
#include "physics/groupstep.h"
#include "physics/hitsling.h"

namespace CordFixed
{
using GroupStep::Ball;
using GroupStep::V3;
using RubberContact::RubberParams;

// A body behind a span: at a cord's centre x, the signed clearance of the cord (radius r), and its gradient. Below zero: the cord is pressed in, however far (the free path is behind the body)
class CordBody
{
public:
   virtual ~CordBody() = default;
   virtual void Gap(const V3 &x, double r, double &g, V3 &m) const = 0;
};
// a round post behind a span, at the cord's height (band_stops); no top: its height isn't read from the table
class StopBody : public CordBody
{
public:
   V3 c0;
   double rCore = 0.0;
   void Gap(const V3 &x, const double r, double &g, V3 &m) const override
   {
      const double dx = x.x - c0.x, dy = x.y - c0.y, p = std::max(std::sqrt(dx * dx + dy * dy), 1e-12);
      m = V3(dx / p, dy / p, 0.0);
      g = p - rCore - r;
   }
};
struct Obstacle
{
   const CordBody *body = nullptr;   // null: an edge of the moving plate
   double x = 0.0;          // where along the span (m from A): a plate's edge, a stop post
   GroupStep::Item law;     // the shared squash law, cord on this body
   double bend = 0.0;       // how far this body holds the cord off its straight line at rest, signed along the span's yh (0: clear of it): part of the band's rest shape and stored energy, ball or no ball
};
struct Span
{
   V3 A, B, t, up = V3(0.0, 0.0, 1.0), yh;   // yh: across the span, in the band's plane
   double len = 0.0, l0 = 0.0, k = 0.0;
   std::vector<Obstacle> ob;
   bool bent = false;       // a body holds this span off its straight line at rest
   double bendMost = 0.0;
   double restV = 0.0, restT = 0.0, restLift = 0.0;   // with no ball: its stored energy and tension, set at load, and how far the cord's centre lifts at the first bend (the rest of the bend is squash)
   void Set(const V3 &a, const V3 &b, const double stretch, const double EA)
   {
      A = a;
      B = b;
      len = (B - A).Length();
      t = (B - A) * (1.0 / len);
      yh = up.Cross(t);
      l0 = len / (1.0 + stretch);
      k = EA / l0;
   }
};

// One side of a taut length: from an end (bx along the axis, dy off it, the ball's side) to where it leaves the ball. Value and slopes are C1 through s = 0: a closed form with no search
struct SideLen
{
   double s = 0.0, beta = 0.0, val = 0.0, Sb = 0.0, Sd = 0.0;
   void Set(const double bx, const double dy, const double R)
   {
      const double D2 = std::max(bx * bx + dy * dy, 1e-18), X = D2 - R * R;
      s = X > 0.0 ? std::sqrt(X) : 0.0;
      beta = std::atan2(R, s) - std::atan2(dy, bx);
      val = s + R * beta;
      Sb = (s * bx + R * dy) / D2;
      Sd = (s * dy - R * bx) / D2;
   }
};

// A ball on a span: its label and the grip's own state. The label changes only with the ball off the cord, where V does not read it
struct Touch
{
   bool on = false;     // touching at the state now
   V3 nh;               // the side the ball came from (+- yh)
   bool top = false;    // it came down onto the cord: the cord under it holds it up (the sphere's own circle); otherwise it came from the side, and the tension never lifts it
   double d0 = 0.0;     // the depth it counts from: nothing stored without the work (the Item's zero)
   double xiT = 0.0, xiB = 0.0, kT = 0.0, kB = 0.0;   // the grip: along the cord, and round it (the twist's spring)
};

class Part final : public GroupStep::Part
{
public:
   static constexpr int kMaxEdges = 4;
   double EA = 126.7, rc = 0.5 * 6.35e-3, mu = 0.8, tanD = 0.15, refHz = 20.0;
   double w = 0.0, T0 = 37.4;   // bending length sqrt(EI / T0): rounds the line's axis
   std::vector<Span> span;
   bool grip = true;
   bool chain = true;   // two or more balls on one side of a span: one taut path round them all (design_two.md); false: each ball's length added, for the comparison
   // the plate: one body (plate and plunger through the rigid link), its spring and rest stop as stored energy, the coil's torque as a drive. Its values and laws: the sling's (hitsling.cpp, coilmodel.h, as bandcord.h)
   struct PlateDef
   {
      bool on = false;
      int span = 0;
      double x0 = 0.0, x1 = 0.0, sgn = 1.0;   // its edges on the span, and which way along yh it pushes
      double originN = 0.0, a = 0.0, h = 0.0, rest = 0.0;
      GroupStep::Item law;
   } pl;
   SlingSettings sl;
   FlipperDerived der;
   CoilSettings coilS;   // the plate's winding (sl.mech's power values), for the coil's one energy law
   CoilDerived coilD;
   SlingDriver driver;
   double phi = 0.0, omega = 0.0, iCoil = 0.0;   // iCoil: A
   bool pinned = false, atRest = false;   // no unknown: on a stop; atRest: the rest stop, otherwise the coil stop
   double minPhi = 0.0;                   // how far it ever went behind its seat (a rigid seat: never)
   mutable double phiQ = 0.0;   // the angle V is asked at
   mutable double dVdP = 0.0;   // d V / d phi from the last Solve with forces
   // what fires it: the cord dented this far at the ball (the real switch's gap), or a time. Reading it does no work
   double swGapM = 0.0, fireAtS = -1.0, tNow = 0.0;
   double maxPhi = 0.0, coilJ = 0.0;
   // the book, by what took it (diagnostics: nothing reads them)
   GroupStep::Item stretchLaw;   // the cord's stretch loss: only its Hunt-Crossley coefficient is set
   double bkStretch = 0.0, bkTwist = 0.0;
   bool twistLoss = true;   // the twist spring takes the rubber's tan_delta
   bool patchShear = true;  // the ball's patch in series with the grip, on the ring's sideways stiffness, friction limit and tan_delta
   GroupStep::Item patchLaw;
   double patchKH = 0.0, patchBeta = 0.0;
   double twistScale = 1.0;   // the twist spring against its calculated value
   double bkDash = 0.0, bkSquash = 0.0, bkGrip = 0.0, bkPivot = 0.0, bkRest = 0.0, bkStop = 0.0, bkTube = 0.0, bkCord = 0.0;
   void Loads(const double *q, const V3 &) override { gLoad = q[0]; }
   // diagnostics
   double maxT = 0.0;
   long nFar = 0;

   int nb = 0, ns = 0;
   std::vector<Touch> tch;

   void SetupPlate()
   {
      FlipperModel::Derive(sl.mech, der);
      coilS.vSupply = sl.mech.vSupply;
      coilS.rOhm = sl.mech.rP;
      coilS.lAirMH = sl.mech.lAirP;
      coilS.lArmMH = sl.mech.lArmP;
      coilS.bobbinMm = sl.mech.bobbinMm;
      coilS.x0 = sl.mech.x0P;
      CoilModel::Derive(coilS, coilD);
      driver.SetPulse(sl.pulseMs * 1e-3);
      // seated on its rest stop: metal or wood, far stiffer than anything the ball meets on the way (user 2026-11-01), and the motion during the contact doesn't matter: a hard stop with an impact law, no stiff spring (ModelStandards 3a)
      phi = 0.0;
      omega = iCoil = 0.0;
      pinned = atRest = true;
      phiQ = phi;
   }
   void Setup(const int balls)
   {
      nb = balls;
      ns = (int)span.size();
      const double r2 = rc * rc, EI = (EA / (kPi_ * r2)) * 0.25 * kPi_ * r2 * r2;
      w = std::sqrt(EI / T0);
      tch.assign((size_t)nb * ns, Touch());
      st.assign((size_t)nb * ns, Scr());
      V0S.assign(ns, 0.0);
      dLS.assign(ns, 0.0);
      TdS.assign(ns, 0.0);
      TmS.assign(ns, 0.0);
      gv.assign(nb, V3());
      farS.assign(ns, 0);
      save = tch;
   }
   // With the ball off the cord (or joining): which way it comes, and where it counts from so that nothing is stored. V is zero before and after: the label can't alter it
   void Label(const Span &sp, Touch &u, const Ball &b) const
   {
      const V3 w0 = b.c - sp.A;
      phiQ = phi;
      const double aL = std::min(std::max(w0.Dot(sp.t), 0.0), sp.len);
      const double y = w0.Dot(sp.yh) - BaseY(sp, aL).y, dz = w0.Dot(sp.up), Rc = b.R + rc;
      u.nh = sp.yh * (y >= 0.0 ? 1.0 : -1.0);
      u.top = dz > 0.0 && -b.v.Dot(sp.up) > std::fabs(b.v.Dot(sp.yh));
      u.on = false;
      u.d0 = 0.0;
      Out o;
      Length(sp, u, b.c, b.R, false, o);
      // from the side above the plane the circle in the plane is smaller than Rc: it touches where the sphere does
      const double dPhys = !u.top && dz > 0.0 ? Rc - std::sqrt(std::max(Rc * Rc - dz * dz, 0.0)) : 0.0;
      u.d0 = std::max({ dPhys, o.D, 0.0 });
   }
   // a ball that joins: its side, and if it is found on the cord, nothing stored
   void Join(const Ball *b)
   {
      for (int i = 0; i < nb; i++)
         for (int j = 0; j < ns; j++)
         {
            Touch &u = tch[(size_t)i * ns + j];
            u = Touch();
            Label(span[j], u, b[i]);
         }
   }

   // GroupStep::Part
   int NumDof() const override { return pl.on && !pinned ? 1 : 0; }
   int MaxDof() const override { return pl.on ? 1 : 0; }
   int NumCoords() const override { return pl.on ? 1 : 0; }
   int DofOf(int) const override { return NumDof() - 1; }
   GroupStep::CoordInfo Coord(int) const override
   {
      GroupStep::CoordInfo q;
      q.rate = omega;
      q.invM = 1.0 / Inertia(phi);
      q.reach = kArm;
      q.maxMove = 0.05;
      q.settle = 1e-9;
      q.invMFree = 1.0 / Inertia(phi);
      q.blockedSide = pinned ? (atRest ? -1 : 1) : 0;
      return q;
   }
   double Energy(const Ball *b, const int nBalls) const override
   {
      phiQ = phi;
      double e = pl.on ? PlateV(phi) + 0.5 * Inertia(phi) * omega * omega : 0.0;
      thread_local std::vector<Out> oo;
      oo.resize(nBalls);
      for (int j = 0; j < ns; j++)
      {
         for (int i = 0; i < nBalls; i++)
         {
            Length(span[j], tch[(size_t)i * ns + j], b[i].c, b[i].R, false, oo[i]);
            const Touch &u = tch[(size_t)i * ns + j];
            e += 0.5 * u.kT * u.xiT * u.xiT + 0.5 * u.kB * u.xiB * u.xiB;
         }
         Chain(span[j], [&](const int i) -> Out & { return oo[i]; }, nBalls, false);
         double L, T;
         e += Solve(span[j], [&](const int i) -> const Out & { return oo[i]; }, nBalls, L, T, nullptr);
      }
      return e;
   }
   bool BeginStep(const Ball *b, const int, const double h, double *x, int off) override
   {
      b0 = b;
      phiQ = phi;
      if (pl.on)
      {
         gI0 = Inertia(phi);
         gp0 = gI0 * omega;
         gVolts = driver.Coil() ? coilS.vSupply : 0.0;
         gCoilLive = gVolts != 0.0 || iCoil != 0.0;
         gKin0 = FlipperModel::Kinematics(sl.mech, der, der.theta0 + phi);
         if (!pinned)
            x[off] = (omega + 0.5 * gAlpha * h) * h;
      }
      for (int j = 0; j < ns; j++)
         if (farS[j])
            V0S[j] = span[j].restV;
         else
            BeginSpan(b, j);
      return true;
   }
   // The patch's shear, as Group::TakeLaw takes it: the law of a ball on a straight tube of the ring's rubber
   void SetPatchLaw(const RubberContact::RubberParams &rp)
   {
      const double E = rp.eMpa * 1e6, Estar = E / (1.0 - rp.nu * rp.nu), G = E / (2.0 * (1.0 + rp.nu)), Gstar = G / (2.0 - rp.nu);
      patchLaw.ktPerA = 8.0 * Gstar * rp.shearScale;
      patchLaw.hL = 2.0 * rc;
      patchKH = 4.0 / 3.0 * Estar;
      patchBeta = rp.betaLayer;
   }
   // Mindlin's stiffness at the patch this force makes: the shared law's depth at the force, its Kt there
   double PatchKt(const double N, const double R) const
   {
      if (N <= 0.0)
         return 0.0;
      GroupStep::Item L = patchLaw;
      const double R1 = 1.0 / (1.0 / R + 1.0 / rc);
      L.reff = std::sqrt(R1 * R);
      const double sr = std::sqrt(L.reff);
      L.kH = patchKH * sr;
      L.cL = L.hL > 0.0 ? L.kH * patchBeta * L.reff * sr : 0.0;
      double d = std::cbrt(N / L.kH);
      d = d * d;
      for (int it = 0; it < 3; it++)
         d = std::min(std::max(d - (L.F(d) - N) / L.K(d), 0.25 * d), 0.999 * L.hL);
      return L.Kt(d);
   }
   void BeginSpan(const Ball *b, const int j)
   {
      const Span &sp = span[j];
      double L, T;
      if (startReady && nb == 1)   // the step's rule has just worked this out at the same positions
      {
         V0S[j] = startV[j];
         T = startT[j];
      }
      else
      {
         if (!startReady)
            for (int i = 0; i < nb; i++)
               Length(sp, tch[(size_t)i * ns + j], b[i].c, b[i].R, true, st[(size_t)i * ns + j].o0);
         Chain(sp, [&](const int i) -> Out & { return st[(size_t)i * ns + j].o0; }, nb, true);
         V0S[j] = Solve(sp, [&](const int i) -> const Out & { return st[(size_t)i * ns + j].o0; }, nb, L, T, nullptr);
      }
      for (int i = 0; i < nb; i++)
      {
         Touch &u = tch[(size_t)i * ns + j];
         Scr &q = st[(size_t)i * ns + j];
         const Out &o = q.o0;
         // the grip's stiffness for this step (the shared Mindlin-Deresiewicz rule: one that rises holds the force, one that falls holds the stretch; what it stored only falls: heat)
         const double aw = std::max(o.a, w), bw = std::max(sp.len - o.a, w), l0a = sp.l0 * aw / sp.len, l0b = sp.l0 * bw / sp.len;
         double kT = EA * (1.0 / l0a + 1.0 / l0b), kB = kT / 6.0 * twistScale;
         if (patchShear)   // the ball's patch on the cord, in series with the cord's own stretch and twist: the shared law's shear, at the patch the pull makes
         {
            const double kP = PatchKt(T * o.Gfree, b[i].R);
            kT = kP > 0.0 ? kT * kP / (kT + kP) : 0.0;
            kB = kP > 0.0 && kB > 0.0 ? kB * kP / (kB + kP) : 0.0;
         }
         q.kT = kT;
         q.kB = kB;
         q.xiT = (kT > u.kT && u.kT > 0.0) ? u.xiT * u.kT / kT : u.xiT;
         q.xiB = (kB > u.kB && u.kB > 0.0) ? u.xiB * u.kB / kB : u.xiB;
         q.lossK = 0.5 * u.kT * u.xiT * u.xiT + 0.5 * u.kB * u.xiB * u.xiB - 0.5 * q.kT * q.xiT * q.xiT - 0.5 * q.kB * q.xiB * q.xiB;
         q.f = grip ? mu * T * o.Gfree : 0.0;
      }
   }
   // A span no ball can touch this step: nothing to do but its constant V. Clear by more than a move the step can make (checked in Eval), not holding anything; the plate's span always counts as it carries the plate's own terms
   bool SpanFar(const Ball *b, const int j) const
   {
      if (pl.on && j == pl.span)
         return false;
      const Span &sp = span[j];
      for (int i = 0; i < nb; i++)
      {
         const Touch &u = tch[(size_t)i * ns + j];
         if (u.on || u.xiT != 0.0 || u.xiB != 0.0)
            return false;
         const V3 w0 = b[i].c - sp.A;
         const double aRaw = w0.Dot(sp.t), a = std::min(std::max(aRaw, 0.0), sp.len), over = aRaw - a, y = w0.Dot(sp.yh), lim = b[i].R + rc + kFarM + sp.bendMost;
         if (y * y + over * over <= lim * lim)
            return false;
      }
      return true;
   }
   bool Eval(const double *x, int off, const double h, V3 *F, V3 *M, double *Res, GroupStep::Jac J) override
   {
      const bool plDof = pl.on && !pinned;
      const double dphi = plDof ? x[off] : 0.0, phi1 = phi + dphi, phm = phi + 0.5 * dphi;
      double tauCord = 0.0;
      for (int j = 0; j < ns; j++)
      {
         const Span &sp = span[j];
         if (farS[j])
         {
            bool still = true;
            for (int i = 0; i < nb && still; i++)
               still = x[6 * i] * x[6 * i] + x[6 * i + 1] * x[6 * i + 1] + x[6 * i + 2] * x[6 * i + 2] < kFarMove * kFarMove;
            if (still)
            {
               maxT = std::max(maxT, sp.restT);
               continue;
            }
            farS[j] = 0;
            BeginSpan(b0, j);
         }
         double dn = 0.0, dLf = 0.0;
         for (int i = 0; i < nb; i++)
         {
            const Ball &bb = b0[i];
            const V3 dc(x[6 * i], x[6 * i + 1], x[6 * i + 2]);
            const Touch &u = tch[(size_t)i * ns + j];
            Scr &q = st[(size_t)i * ns + j];
            phiQ = phi1;
            Length(sp, u, bb.c + dc, bb.R, false, q.o1);
            phiQ = phm;
            Length(sp, u, bb.c + dc * 0.5, bb.R, true, q.om);
         }
         Chain(sp, [&](const int i) -> Out & { return st[(size_t)i * ns + j].o1; }, nb, false);
         Chain(sp, [&](const int i) -> Out & { return st[(size_t)i * ns + j].om; }, nb, true);
         double gdSum = 0.0;
         for (int i = 0; i < nb; i++)
         {
            const V3 dc(x[6 * i], x[6 * i + 1], x[6 * i + 2]);
            const Scr &q = st[(size_t)i * ns + j];
            dn += dc.Dot(dc);
            dLf += q.o1.L - q.o0.L;
            gdSum += q.om.G.Dot(dc) + q.om.LP * dphi;
         }
         dn += kArm * kArm * dphi * dphi;
         const double rsAll = dn > 1e-24 ? (dLf - gdSum) / dn : 0.0;   // the dashpot's path closed over the whole span, as V's
         double L1, T1, Lm, Tm;
         phiQ = phi1;
         const double V1 = Solve(sp, [&](const int i) -> const Out & { return st[(size_t)i * ns + j].o1; }, nb, L1, T1, nullptr);
         phiQ = phm;
         Solve(sp, [&](const int i) -> const Out & { return st[(size_t)i * ns + j].om; }, nb, Lm, Tm, gv.data());
         const double dVp = dVdP;
         TmS[j] = Tm;
         maxT = std::max(maxT, T1);
         // Gonzalez over the whole span: the mid-step force, closed along the move to the exact difference of V
         double mid = 0.0;
         for (int i = 0; i < nb; i++)
            mid += gv[i].Dot(V3(x[6 * i], x[6 * i + 1], x[6 * i + 2]));
         mid += dVp * dphi;
         // (a difference within the rounding of V and of the length it is made from is no difference: divided by the move, it would be a force that grows as the ball slows, in any step length)
         const double missV = (V1 - V0S[j]) - mid, roundV = std::numeric_limits<double>::epsilon() * (std::fabs(V1) + std::fabs(V0S[j]) + 2.0 * Tm * Lm);
         const double close = dn > 1e-24 && std::fabs(missV) > roundV ? missV / dn : 0.0;
         // the cord's own loss: a dashpot on the stretch rate (the shared law's form), along the free path
         const double eta = EA * tanD / (2.0 * kPi_ * refHz), Td = eta * dLf / (h * sp.l0);
         dLS[j] = dLf;
         TdS[j] = Td;
         tauCord -= dVp + kArm * kArm * dphi * close;
         for (int i = 0; i < nb; i++)
         {
            const Ball &bb = b0[i];
            const V3 dc(x[6 * i], x[6 * i + 1], x[6 * i + 2]), dw(x[6 * i + 3], x[6 * i + 4], x[6 * i + 5]);
            Scr &q = st[(size_t)i * ns + j];
            const double dc2 = dc.Dot(dc);
            F[i] = F[i] - gv[i] - dc * close;
            const double dn1 = dc2 + kArm * kArm * dphi * dphi;
            const V3 Gd = q.om.G + dc * rsAll;
            const double GdP = q.om.LP + kArm * kArm * dphi * rsAll;
            F[i] = F[i] - Gd * Td;
            tauCord -= GdP * Td;
            // Hunt-Crossley on the stretch: the cord's pull on the ball grows with how fast the ball deepens its dent, by the shared law's coefficient (Item::DampC) with the pull in place of the squash's force, push-only. On the dent's exact difference, closed along the move as the dashpot's is, so it only takes energy out; V doesn't read it
            q.lossS = 0.0;
            const double feS = Tm * std::max(-q.om.Ld, 0.0);
            if (stretchLaw.alpha > 0.0 && feS > 0.0)
            {
               const double dD = q.o1.Deff - q.o0.Deff;
               const V3 gD = q.om.gR - q.om.gC;
               const double gP = -q.om.cenP, rsD = dn1 > 1e-24 ? (dD - gD.Dot(dc) - gP * dphi) / dn1 : 0.0, cD = stretchLaw.DampC(q.om.Deff, feS), fS = std::max(cD * dD / h, -feS);
               const V3 gDd = gD + dc * rsD;
               F[i] = F[i] - gDd * fS;
               tauCord -= (gP + kArm * kArm * dphi * rsD) * fS;
               q.lossS = fS * dD;
               if (J)
               {
                  const double gg[3] = { gD.x, gD.y, gD.z };
                  for (int a = 0; a < 3; a++)
                     for (int c = 0; c < 3; c++)
                        J[6 * i + a][6 * i + c] += cD * gg[a] * gg[c];
               }
            }
            // the grip
            q.xiT1 = q.xiB1 = 0.0;
            q.lossTw = 0.0;
            if (q.f > 0.0 && (q.om.Deff > 0.0 || q.o1.Deff > 0.0))
            {
               const V3 &tt = sp.t;
               const V3 nm0 = q.o0.n + q.o1.n;
               const double nl = nm0.Length();
               const V3 nm = nl > 1e-9 ? nm0 * (1.0 / nl) : q.o1.n;
               const V3 bm = nm.Cross(tt);
               const V3 dth = (bb.w + dw * 0.5) * h;
               const V3 Rn = nm * bb.R;
               const double spT = dth.Cross(Rn).Dot(tt), spB = dth.Cross(Rn).Dot(bm);
               // the cord moves with the ball only across the span, in the plane; along it and up and down it stays
               const double St = dc.Dot(tt) + spT, upB = sp.up.Dot(bm), Sb = dc.Dot(sp.up) * upB + spB;
               double fT = 0.0, fB = 0.0;
               Chord(q.kT, q.f, q.xiT, St, fT, q.xiT1, q.pT0, q.pT1);
               Chord(q.kB, q.f, q.xiB, Sb, fB, q.xiB1, q.pB0, q.pB1);
               // the rubber's own loss on the twist spring: tan_delta, the shared law's shear damper (ContactCore: tan_delta sqrt(k mt), spring and damper within the limit), on the step's own slip, so its work is the loss booked; none while it slides
               q.lossTw = 0.0;
               if (twistLoss)
               {
                  const double mt = 1.0 / (1.0 / bb.m + bb.R * bb.R / bb.I());
                  const auto Damp = [&](const double k, const double xi1, const double S, double &fc) {
                     if (k <= 0.0 || std::fabs(xi1) * k > q.f)
                        return;
                     const double fD = std::min(std::max(fc + tanD * std::sqrt(k * mt) * S / h, -q.f), q.f) - fc;
                     fc += fD;
                     q.lossTw += fD * S;
                  };
                  Damp(q.kB, q.xiB1, Sb, fB);
                  if (patchShear)
                     Damp(q.kT, q.xiT1, St, fT);
               }
               const V3 Fs = (tt * fT + bm * fB) * -1.0;
               F[i] = F[i] - tt * fT - sp.up * (upB * fB);
               M[i] = M[i] + Rn.Cross(Fs);
            }
            else
            {
               q.pT0 = q.pT1 = q.pB0 = q.pB1 = 0.0;
            }
            if (J)
            {
               // Newton's slope: the cord is soft beside 2 m / h^2, so the outer product and the string's own pull are enough
               const double ga[3] = { q.om.G.x, q.om.G.y, q.om.G.z }, tv[3] = { sp.t.x, sp.t.y, sp.t.z };
               const double aw = std::max(q.om.a, w), bw = std::max(sp.len - q.om.a, w), kg = q.om.Deff > 0.0 && !q.om.shadow ? Tm * (1.0 / aw + 1.0 / bw) : 0.0;
               for (int a = 0; a < 3; a++)
                  for (int c = 0; c < 3; c++)
                     J[6 * i + a][6 * i + c] += 0.5 * h * kg * ((a == c ? 1.0 : 0.0) - tv[a] * tv[c]);
               for (int k = 0; k < nb; k++)
               {
                  const Out &ok = st[(size_t)k * ns + j].om;
                  const double gb[3] = { ok.G.x, ok.G.y, ok.G.z };
                  for (int a = 0; a < 3; a++)
                     for (int c = 0; c < 3; c++)
                        J[6 * i + a][6 * k + c] += 0.5 * h * sp.k * ga[a] * gb[c];
               }
            }
         }
      }
      phiQ = phi;
      gPlateLoss = gCoilW = 0.0;
      double tauC = 0.0;
      if (pl.on)
      {
         gOver = gCoilLive ? CoilOver(dphi, h) : CoilModel::OverMove{};
         tauC = gCoilLive ? CoilTorque(gOver, dphi) : 0.0;
      }
      if (pl.on && !plDof)
         gTauPinned = -PlateVp(phi) + tauC + tauCord;
      if (plDof)
      {
         // the plate's own row: its momentum, its spring and rest stop as exact differences, its losses on the mid-step rate, the coil's torque a drive (as bandcord.h)
         const double I1 = Inertia(phi1), Jb = 0.5 * (1.0 / gI0 + 1.0 / I1), pEnd = 2.0 * dphi / (h * Jb) - gp0;
         double dTdphi;
         if (std::fabs(dphi) > 1e-7)
            dTdphi = (1.0 / I1 - 1.0 / gI0) / dphi * 0.25 * (gp0 * gp0 + pEnd * pEnd);
         else
         {
            const double e = 1e-6;
            dTdphi = (1.0 / Inertia(phm + e) - 1.0 / Inertia(phm - e)) / (2.0 * e) * 0.25 * (gp0 * gp0 + pEnd * pEnd);
         }
         const double Vp = std::fabs(dphi) > 1e-7 ? (PlateV(phi1) - PlateV(phi)) / dphi : PlateVp(phm);
         double tauS, cR;
         RestStop(phm, tauS, cR);
         const double om = dphi / h, tauD = tauS > 0.0 ? std::max(-cR * om, -tauS) : 0.0;
         const double bd = sl.mech.bDamping;
         const double Fl = -Vp + tauD - bd * om + tauC - dTdphi + tauCord;
         gPlateLoss = -(tauD - bd * om) * dphi;
         gRestL = -tauD * dphi;
         gPivotL = bd * om * dphi;
         gCordTau = tauCord;
         gCoilW = tauC * dphi;
         Res[off] = pEnd - gp0 - h * Fl;
         if (J)
         {
            // (the coil's term by differences of the same closed form: the magnet's pull, a negative stiffness, and the closing speed's back EMF)
            const double e = 1e-9 + 1e-6 * std::fabs(dphi);
            const double dTau = gCoilLive ? (CoilTorque(CoilOver(dphi + e, h), dphi + e) - CoilTorque(CoilOver(dphi - e, h), dphi - e)) / (2.0 * e) : 0.0;
            J[off][off] += 2.0 / (h * Jb) + 0.5 * h * std::max(PlateVpp(phm), 0.0) + bd + cR - h * dTau;
         }
      }
      return true;
   }
   double StopExcess(const double *x, int off) const override { return pl.on && !pinned ? std::max(phi + x[off] - PhiStop(), -(phi + x[off])) : -INFINITY; }
   bool PinOnStop(const double *, int) override
   {
      const bool low = phi < 0.5 * PhiStop();
      if (pinned || (low ? omega > 0.0 : omega < 0.0))
         return false;
      pinned = true;
      atRest = low;
      gPinJ += 0.5 * Inertia(phi) * omega * omega;   // the plunger stops on the steel coil stop: inelastic
      omega = 0.0;
      return true;
   }
   // off the rest stop: the step's own forces take it off, the coil's or the cord's against the ball's reaction
   bool FreedByStep() override
   {
      if (!pl.on || !pinned || !atRest || gFreed || gTauPinned + gLoad <= 0.0)
         return false;
      pinned = atRest = false;
      gFreed = true;
      return true;
   }
   double KineticJ() const override { return pl.on ? 0.5 * Inertia(phi) * omega * omega : 0.0; }
   void ScaleSpeeds(double f) override { omega *= f; }
   // the limit from the middle of the step
   double UpdateGrip(const double *, double h) override
   {
      if (!grip)
         return 0.0;
      double ch = 0.0;
      for (int j = 0; j < ns; j++)
         for (int i = 0; i < nb && !farS[j]; i++)
         {
            Scr &q = st[(size_t)i * ns + j];
            const double f = mu * TmS[j] * q.om.Gfree;
            ch = std::max(ch, std::fabs(f - q.f) * h / b0[i].m);
            q.f = f;
         }
      return ch;
   }
   void StepFacts(const Ball *b, const int, const GroupStep::Settings &, GroupStep::StepFactList &f) override
   {
      if (pl.on)
      {
         phiQ = phi;
         double yl, dy;
         Limit(phi, yl, dy);
         if (!pinned)
         {
            f.Stiffness(std::sqrt(std::max(PlateVpp(phi), 0.0) / Inertia(phi)), false, 103);
            f.Field(w, std::fabs(omega) * dy, 0.0, 104);   // the plate moves the cord by no more than a share of its bending length
         }
         if (driver.Coil() || iCoil > 1e-6 * coilD.iss)
         {
            const FlipperSettings &s = sl.mech;
            const FlipperModel::Kin k = FlipperModel::Kinematics(s, der, der.theta0 + phi);
            const double frac = (k.x - der.xMin) / s.bobbinMm;
            f.TimeConstant((der.aP + der.bP / (frac + s.x0P)) * 1e-3 / s.rP, 105);
         }
      }
      startV.resize(ns);
      startT.resize(ns);
      for (int j = 0; j < ns; j++)
      {
         const Span &sp = span[j];
         for (int i = 0; i < nb; i++)
            Length(sp, tch[(size_t)i * ns + j], b[i].c, b[i].R, true, st[(size_t)i * ns + j].o0);
         double L, T;
         startV[j] = Solve(sp, [&](const int i) -> const Out & { return st[(size_t)i * ns + j].o0; }, nb, L, T, nullptr);
         startT[j] = T;
         for (int i = 0; i < nb; i++)
         {
            const Out &o = st[(size_t)i * ns + j].o0;
            const double aw = std::max(o.a, w), bw = std::max(sp.len - o.a, w);
            double kTens = sp.k * o.G.Dot(o.G) + T * (1.0 / aw + 1.0 / bw);
            for (int e = 0; e < o.nEd; e++)
               kTens += T / std::max(o.ed[e].a1, w);   // a run bent round a body: its pull rate
            const double omega = std::sqrt(kTens / b[i].m);
            if (o.Deff > 0.0)
            {
               f.Stiffness(omega, false, 100);
            }
            else
            {
               const double rate = -o.ohat.Dot(b[i].v);
               if (rate > 1e-6 && o.gap < 0.05)
                  f.Approach(std::max(o.gap, 0.0), rate, omega, 102);
            }
         }
      }
      startReady = true;
   }
   void Accept(const double *x, int off, double h, bool onStop, double &driveJ, double &lossJ, bool &) override
   {
      if (pl.on)
      {
         lossJ += gPinJ;
         bkStop += gPinJ;
         gPinJ = 0.0;
         if (!pinned)
         {
            const double dq = x[off], phi1 = phi + dq, I1 = Inertia(phi1), Jb = 0.5 * (1.0 / gI0 + 1.0 / I1);
            double om1 = (2.0 * dq / (h * Jb) - gp0) / I1;
            gAlpha = (om1 - omega) / h;
            if (onStop)
            {
               lossJ += 0.5 * I1 * om1 * om1;
               bkStop += 0.5 * I1 * om1 * om1;
               om1 = 0.0;
               pinned = true;
               atRest = phi1 < 0.5 * PhiStop();
            }
            phi = phi1;
            omega = om1;
            maxPhi = std::max(maxPhi, phi);
            minPhi = std::min(minPhi, phi);
            driveJ += gCoilW;
            coilJ += gCoilW;
            lossJ += gPlateLoss;
            bkRest += gRestL;
            bkPivot += gPivotL;
            bkTube += gLoad * dq;
            bkCord += gCordTau * dq;
         }
         iCoil = gOver.i1;
         gFreed = false;
         driver.Tick(h);
         tNow += h;
         phiQ = phi;
      }
      for (int j = 0; j < ns; j++)
      {
         if (farS[j])
            continue;
         lossJ += TdS[j] * dLS[j];
         bkDash += TdS[j] * dLS[j];
         for (int i = 0; i < nb; i++)
         {
            Touch &u = tch[(size_t)i * ns + j];
            Scr &q = st[(size_t)i * ns + j];
            // what the step's limit cut at the start, what stays beyond it at the end, and what the stiffness rule let go: all only falls
            lossJ += q.lossS + q.lossTw;
            bkStretch += q.lossS;
            bkTwist += q.lossTw;
            lossJ += q.lossK;
            const double lossBefore = lossJ;
            const auto Close = [&](const double k, const double f, const double xiStart, const double p0, const double p1, const double xi1, double &keep) {
               if (k <= 0.0 || f <= 0.0)
               {
                  keep = 0.0;
                  return 0.5 * k * xiStart * xiStart;
               }
               const double lim = f / k;
               keep = std::min(std::max(xi1, -lim), lim);
               const double cut = 0.5 * k * xiStart * xiStart - p0;
               return cut + p1 - 0.5 * k * keep * keep;
            };
            double keepT = 0.0, keepB = 0.0;
            const bool active = q.f > 0.0 && (q.om.Deff > 0.0 || q.o1.Deff > 0.0);
            if (active)
            {
               lossJ += Close(q.kT, q.f, q.xiT, q.pT0, q.pT1, q.xiT1, keepT);
               lossJ += Close(q.kB, q.f, q.xiB, q.pB0, q.pB1, q.xiB1, keepB);
               if (q.o1.Deff <= 0.0)   // off the cord at the step's end: what it held goes as heat
               {
                  lossJ += 0.5 * q.kT * keepT * keepT + 0.5 * q.kB * keepB * keepB;
                  keepT = keepB = 0.0;
               }
            }
            else
               lossJ += 0.5 * q.kT * q.xiT * q.xiT + 0.5 * q.kB * q.xiB * q.xiB;
            bkGrip += lossJ - lossBefore + q.lossK;
            u.kT = q.kT;
            u.kB = q.kB;
            u.xiT = keepT;
            u.xiB = keepB;
         }
      }
   }
   void AfterStep(const GroupStep::Group &g) override
   {
      phiQ = phi;
      // what each ball keeps on each post, stop and the plate between the cycles (a group is built afresh each cycle)
      for (int k = 0; k < g.ni; k++)
      {
         const GroupStep::Item &it = g.item[k];
         if (it.b >= 0 || it.s < postFirst || it.s >= postFirst + nShapes)
            continue;
         Keep &kp = keep[(size_t)it.a * nShapes + (it.s - postFirst)];
         kp.touching = it.touching;
         kp.xi = it.xi;
         kp.ra = it.ra;
         kp.zero = it.dZero;
      }
      for (int i = 0; i < nb; i++)
      {
         lastC[i] = g.ball[i].c;
         lastR[i] = g.ball[i].R;
      }
      if (pl.on)
      {
         bool sw = false;
         for (int i = 0; i < nb && swGapM > 0.0; i++)
         {
            Out o;
            Length(span[pl.span], tch[(size_t)i * ns + pl.span], g.ball[i].c, g.ball[i].R, false, o);
            sw = sw || o.Deff >= swGapM;
         }
         if (swGapM > 0.0 && enabled)
            driver.Update(sw);
         if (fireAtS >= 0.0 && tNow >= fireAtS)
         {
            driver.Trigger();
            fireAtS = -1.0;
         }
         if (pinned && !atRest)
         {
            // off the coil stop once the torque on it points away
            thread_local std::vector<Out> oo;
            thread_local std::vector<V3> gg;
            oo.resize(nb);
            gg.resize(nb);
            const Span &sp = span[pl.span];
            for (int i = 0; i < nb; i++)
               Length(sp, tch[(size_t)i * ns + pl.span], g.ball[i].c, g.ball[i].R, true, oo[i]);
            Chain(sp, [&](const int i) -> Out & { return oo[i]; }, nb, true);
            double L, T;
            Solve(sp, [&](const int i) -> const Out & { return oo[i]; }, nb, L, T, gg.data());
            if (-dVdP - PlateVp(phi) + CoilTorqueAt(phi) < 0.0)
               pinned = false;
         }
      }
      startReady = false;
      for (int j = 0; j < ns; j++)
         farS[j] = SpanFar(&g.ball[0], j) ? 1 : 0;
      for (int i = 0; i < nb; i++)
         for (int j = 0; j < ns; j++)
         {
            if (farS[j])
               continue;
            Touch &u = tch[(size_t)i * ns + j];
            Out o;
            Length(span[j], u, g.ball[i].c, g.ball[i].R, false, o);
            u.on = o.Deff > 0.0;
            if (o.pastLine)
               nFar++;
            if (!u.on)
               Label(span[j], u, g.ball[i]);
         }
   }
   void SaveState() override { save = tch; savePl[0] = phi; savePl[1] = omega; savePl[2] = iCoil; savePin = pinned; saveRest = atRest; }
   void RestoreState() override { tch = save; phi = phiQ = savePl[0]; omega = savePl[1]; iCoil = savePl[2]; pinned = savePin; atRest = saveRest; }
   // ---- what the engine's cord path asks of a band (cordband.h, contactsystem.cpp, hitrubber.cpp) ----
   struct Circle { double x = 0.0, y = 0.0, r = 0.0; bool twin = false; };   // r: the cord centre's radius round it
   struct SwitchSpec { double x = 0.0, y = 0.0, bx = 0.0, by = 0.0; bool coil = false; };
   static constexpr int kMaxSwitches = 8;
   struct Geo
   {
      bool hasPlate = false;
      int nSwitches = 0;
      SwitchSpec sw[kMaxSwitches];
   } geo;
   std::vector<Circle> circles;   // posts and stops, the Items a ball meets at the cord's height
   double zc = 0.0;               // the cord's height, world
   RubberParams ringRp;           // the table's rubber: the ring's law
   bool enabled = true;           // the sling is on (script): no coil from the switch otherwise
   double switchGapM = 0.0;
   long energySteps = 0;          // counted by the contact system: none expected
   void AttachGroup(const Ball *b, const int n)
   {
      Setup(n);
      keep.assign((size_t)n * nShapes, Keep());
      lastC.assign(n, V3());
      lastR.assign(n, 0.0);
      for (int i = 0; i < n; i++)
      {
         lastC[i] = b[i].c;
         lastR[i] = b[i].R;
      }
      if (n > 0)
         Join(b);
   }
   void RemapBalls(const int *oldOfNew, const int n, const Ball *b)
   {
      const std::vector<Touch> oldT = tch;
      const std::vector<Keep> oldK = keep;
      const int ob = nb;
      Setup(n);
      keep.assign((size_t)n * nShapes, Keep());
      lastC.assign(n, V3());
      lastR.assign(n, 0.0);
      for (int i = 0; i < n; i++)
      {
         lastC[i] = b[i].c;
         lastR[i] = b[i].R;
         const int o = oldOfNew[i];
         for (int j = 0; j < ns; j++)
         {
            Touch &u = tch[(size_t)i * ns + j];
            if (o >= 0 && o < ob)
               u = oldT[(size_t)o * ns + j];
            else
               Label(span[j], u, b[i]);
         }
         for (int s = 0; s < nShapes && o >= 0 && o < ob; s++)
            keep[(size_t)i * nShapes + s] = oldK[(size_t)o * nShapes + s];
      }
   }
   // the start of a cycle with its balls: what was worked out for the last cycle's end isn't kept; a ball off the cord may have moved since, its side and zero are read again (V doesn't read them)
   void BeginGroupCycle(const Ball *b)
   {
      startReady = false;
      phiQ = phi;
      for (int j = 0; j < ns; j++)
         farS[j] = 0;
      for (int i = 0; i < nb; i++)
      {
         lastC[i] = b[i].c;
         lastR[i] = b[i].R;
         for (int j = 0; j < ns; j++)
         {
            Touch &u = tch[(size_t)i * ns + j];
            if (!u.on)
               Label(span[j], u, b[i]);
         }
      }
   }
   class PostShape : public GroupStep::Shape
   {
   public:
      const Part *pt = nullptr;
      int j = 0;
      void Probe(const V3 &c, const double R, GroupStep::Hit &h) const override
      {
         const Circle &q = pt->circles[j];
         if (q.twin)
         {
            h.delta = -1.0;
            h.n = V3(0.0, 0.0, 1.0);
            h.reff = R;
            h.curvDist = 0.0;
            return;
         }
         const double dx = c.x - q.x, dy = c.y - q.y, D = std::sqrt(dx * dx + dy * dy), aa = D - q.r, zz = c.z - pt->zc, dist = std::sqrt(aa * aa + zz * zz);
         const double ex = D > 0.0 ? dx / D : 1.0, ey = D > 0.0 ? dy / D : 0.0, na = dist > 0.0 ? aa / dist : 1.0, nz = dist > 0.0 ? zz / dist : 0.0;
         h.delta = R + pt->rc - dist;
         h.n = V3(ex * na, ey * na, nz);
         const double R1 = 1.0 / (1.0 / R + 1.0 / pt->rc), kLong = RubberContact::TubeAlongCurvature(1.0 / q.r, pt->rc, 1.0);   // the effective radius as ProbeTube
         h.reff = std::sqrt(R1 / std::max(1.0 / R + kLong, 0.05 / R));
         h.curvDist = D;
      }
   };
   // Into the group, after AddPart, before its Configure
   template <class Tube> void AddPostShapes(GroupStep::Group &g, Tube &tube)
   {
      postFirst = g.ns;
      postShape.resize(circles.size());
      for (size_t k = 0; k < circles.size(); k++)
      {
         PostShape &q = postShape[k];
         q.pt = this;
         q.j = (int)k;
         q.rp = ringRp;
         q.layerM = 2.0 * rc;
         g.AddShape(&q);
      }
      if (pl.on)
      {
         tube.pt = this;
         tube.rp = ringRp;
         tube.layerM = 2.0 * rc;
         tube.part = g.np - 1;
         g.AddShape(&tube);
      }
   }
   // After its Configure: what each ball kept on each. One not touching that is found inside has stored nothing (the shared law's rule, SetShapeItem)
   void RestorePostItems(GroupStep::Group &g) const
   {
      for (int k = 0; k < g.ni; k++)
      {
         const GroupStep::Item &it = g.item[k];
         if (it.b >= 0 || it.s < postFirst || it.s >= postFirst + nShapes)
            continue;
         const Keep &kp = keep[(size_t)it.a * nShapes + (it.s - postFirst)];
         g.SetShapeItem(k, kp.touching, kp.zero, kp.xi, kp.ra);
      }
   }
   bool OnCord(const int i) const
   {
      for (int j = 0; j < ns; j++)
         if (tch[(size_t)i * ns + j].on)
            return true;
      return false;
   }
   // how fast a ball came on, by the span it is on: the script's Hit
   double ClosingOnCord(const int i, const V3 &v) const
   {
      phiQ = phi;
      double best = 0.0;
      for (int j = 0; j < ns; j++)
      {
         Out o;
         const Touch &u = tch[(size_t)i * ns + j];
         if (u.on)
         {
            Length(span[j], u, lastC[i], lastR[i], false, o);
            best = std::max(best, -o.ohat.Dot(v));
         }
      }
      return best;
   }
   // how far a ball is from the cord and from anything of its own it touches: it stays in the group until clear
   double CordGap(const int i) const
   {
      if (i < 0 || i >= nb)
         return 1.0;
      phiQ = phi;
      double gap = 1e30;
      for (int j = 0; j < ns; j++)
      {
         Out o;
         const Touch &u = tch[(size_t)i * ns + j];
         Length(span[j], u, lastC[i], lastR[i], false, o);
         gap = std::min(gap, u.on ? 0.0 : o.gap);
      }
      for (size_t k = 0; k < circles.size(); k++)
      {
         GroupStep::Hit hp;
         postShape.empty() ? (void)0 : postShape[k].Probe(lastC[i], lastR[i], hp);
         if (!postShape.empty())
            gap = std::min(gap, keep[(size_t)i * nShapes + k].touching ? 0.0 : -(hp.delta - keep[(size_t)i * nShapes + k].zero));
      }
      return gap;
   }
   bool PlateQuiet() const { return !geo.hasPlate || (std::fabs(phi) < 1e-5 && std::fabs(omega) < 1e-3 && !driver.Coil() && driver.pending <= 0.0); }
   // the cord's displacement at a switch blade, from the deepest ball's own runs (a blade reads where the cord is: nothing here acts on the ball)
   double SwitchPush(const int s) const
   {
      const SwitchSpec &sw = geo.sw[s];
      phiQ = phi;
      double best = 0.0;
      for (int j = 0; j < ns; j++)
      {
         const Span &sp = span[j];
         const V3 w0 = V3(sw.x, sw.y, sp.A.z) - sp.A;
         const double x = w0.Dot(sp.t);
         if (x < 0.0 || x > sp.len || std::fabs(w0.Dot(sp.yh)) > 1e-2)
            continue;
         const double blade = std::fabs(sp.yh.x * sw.bx + sp.yh.y * sw.by);
         for (int i = 0; i < nb; i++)
         {
            Out o;
            Length(sp, tch[(size_t)i * ns + j], lastC[i], lastR[i], false, o);
            if (o.Deff <= 0.0)
               continue;
            const double dx = x - o.a, right = dx >= 0.0 ? 1.0 : -1.0, bx = dx >= 0.0 ? sp.len - o.a : o.a;
            SideLen S;
            S.Set(bx, o.dt, o.Rr);
            const double sb = std::sin(S.beta), cb = std::cos(S.beta);
            double depth;
            if (std::fabs(dx) <= o.Rr * sb)
               depth = std::sqrt(std::max(o.Rr * o.Rr - dx * dx, 0.0)) - o.dt;   // under the ball's arc
            else
               depth = cb > 1e-6 ? (dx >= 0.0 ? sp.len - x : x) * sb / cb : 0.0;   // on its straight run, down to the post
            (void)right;
            best = std::max(best, depth * blade);
         }
      }
      return best;
   }
   bool SwitchIsClosed(const int s) const { return swClosed[s] != 0; }
   void CycleBegin()
   {
      bool any = false;
      for (int s = 0; s < geo.nSwitches; s++)
      {
         const bool closed = SwitchPush(s) >= switchGapM;
         swClosed[s] = closed ? 1 : 0;
         any = any || (closed && geo.sw[s].coil);
      }
      (void)any;
      if (pl.on && swGapM > 0.0)
      {
         bool sw = false;
         for (int i = 0; i < nb; i++)
         {
            Out o;
            Length(span[pl.span], tch[(size_t)i * ns + pl.span], lastC[i], lastR[i], false, o);
            sw = sw || o.Deff >= swGapM;
         }
         driver.Update(sw);
      }
   }
   // ---- the drawn band (requirement 11): the cord's path across its span, post, tangent, round the balls, tangent, post, read once a frame. No state, nothing the step reads ----
   struct DrawV { double x, y, r, thOut; };
   struct DrawSide { std::vector<DrawV> v; };   // from the first notch to the last, ends with radius 0
   mutable std::vector<DrawSide> drawS;         // two a span: the balls pressing from each side
   bool DrawActive() const
   {
      if (!PlateQuiet())
         return true;
      for (const Span &sp : span)
         if (sp.bent)
            return true;
      for (const Touch &u : tch)
         if (u.on)
            return true;
      return false;
   }
   void DrawPrepare() const
   {
      phiQ = phi;
      drawS.resize((size_t)ns * 2);
      thread_local std::vector<Out> oo;
      oo.resize(nb);
      for (int j = 0; j < ns; j++)
      {
         const Span &sp = span[j];
         for (int i = 0; i < nb; i++)
            Length(sp, tch[(size_t)i * ns + j], lastC[i], lastR[i], false, oo[i]);
         Chain(sp, [&](const int i) -> Out & { return oo[i]; }, nb, false);
         for (int side = 0; side < 2; side++)
         {
            std::vector<DrawV> &v = drawS[(size_t)j * 2 + side].v;
            v.clear();
            const double sg = side == 0 ? 1.0 : -1.0;
            v.push_back({ 0.0, 0.0, 0.0, 0.0 });
            for (int i = 0; i < nb; i++)
               if (oo[i].Deff > 0.0 && !oo[i].shadow && oo[i].sg == sg)
                  v.push_back({ oo[i].a, oo[i].dt, oo[i].Rr, 0.0 });
            std::sort(v.begin() + 1, v.end(), [](const DrawV &p, const DrawV &q) { return p.x < q.x; });
            v.push_back({ sp.len, 0.0, 0.0, 0.0 });
            for (size_t k = 0; k + 1 < v.size(); k++)
            {
               const double dx = v[k + 1].x - v[k].x, dy = v[k + 1].y - v[k].y, dr = v[k + 1].r - v[k].r, X = dx * dx + dy * dy - dr * dr, s = X > 0.0 ? std::sqrt(X) : 0.0;
               v[k].thOut = std::atan2(dy, dx) - std::atan2(dr, s);
            }
         }
      }
   }
   // how far the cord is off its straight line at a station, across the span, in the span's own +yh
   double DrawY(const int j, const double x) const
   {
      const Span &sp = span[j];
      double y = BaseY(sp, std::min(std::max(x, 0.0), sp.len)).y;
      for (int side = 0; side < 2; side++)
      {
         const std::vector<DrawV> &v = drawS[(size_t)j * 2 + side].v;
         if (v.size() <= 2)
            continue;
         double here = 0.0;
         bool found = false;
         for (size_t k = 0; k + 1 < v.size() && !found; k++)
         {
            // the straight run from this one to the next, between where it leaves one and meets the other
            const DrawV &p = v[k], &q = v[k + 1];
            const double sn = std::sin(p.thOut), cs = std::cos(p.thOut), xa = p.x + p.r * sn, ya = p.y - p.r * cs, xb = q.x + q.r * sn, yb = q.y - q.r * cs;
            if (x >= xa && x <= xb && xb > xa)
            {
               here = ya + (yb - ya) * (x - xa) / (xb - xa);
               found = true;
            }
         }
         for (size_t k = 1; k + 1 < v.size() && !found; k++)   // or on a ball's arc
            if (std::fabs(x - v[k].x) <= v[k].r)
            {
               here = v[k].y - std::sqrt(v[k].r * v[k].r - (x - v[k].x) * (x - v[k].x));
               found = true;
            }
         y += (side == 0 ? 1.0 : -1.0) * std::min(here, 0.0);
      }
      return y;
   }
   // a ring of the drawn mesh with its centre at (x, y), m: where it goes and how it turns. False: on no span
   bool DrawAt(const double x, const double y, double &sx, double &sy, double &cosA, double &sinA) const
   {
      int best = -1;
      double bestD = 5e-3, st = 0.0;
      for (int j = 0; j < ns; j++)
      {
         const Span &sp = span[j];
         const double wx = x - sp.A.x, wy = y - sp.A.y, a = wx * sp.t.x + wy * sp.t.y, d = std::fabs(wx * sp.yh.x + wy * sp.yh.y);
         if (a > 0.0 && a < sp.len && d < bestD)
         {
            best = j;
            bestD = d;
            st = a;
         }
      }
      if (best < 0)
         return false;
      const Span &sp = span[best];
      const double Y = DrawY(best, st), e = 1e-3, slope = (DrawY(best, st + e) - DrawY(best, st - e)) / (2.0 * e), turn = std::atan(slope) * (sp.t.x * sp.yh.y - sp.t.y * sp.yh.x >= 0.0 ? 1.0 : -1.0);
      if (std::fabs(Y) < 1e-6 && std::fabs(slope) < 1e-5)
         return false;
      sx = sp.yh.x * Y;
      sy = sp.yh.y * Y;
      cosA = std::cos(turn);
      sinA = std::sin(turn);
      return true;
   }
   // At load, after the spans, their bodies and the plate: a span's rest with no ball. Nothing here is worked out again while no ball is near
   void SetupRest()
   {
      phiQ = phi;
      for (Span &sp : span)
      {
         double L;
         sp.restV = Solve(sp, [](const int) -> const Out & { static const Out none; return none; }, 0, L, sp.restT, nullptr);
         for (const Obstacle &ob : sp.ob)
            if (ob.bend != 0.0 && sp.restLift == 0.0)
            {
               Edge e;
               e.law = &ob.law;
               e.gr = std::fabs(ob.bend);
               e.a1 = ob.x;
               e.a2 = sp.len - ob.x;
               sp.restLift = Share(e, sp.restT);
            }
      }
   }
   void SetNShapes() { nShapes = (int)circles.size() + (pl.on ? 1 : 0); swClosed.assign(kMaxSwitches, 0); }
   int postFirst = 0, nShapes = 0;
   struct Keep { bool touching = false; double zero = 0.0; V3 xi; RubberContact::RollAngle ra; };
   std::vector<Keep> keep;
   std::vector<PostShape> postShape;
   std::vector<V3> lastC;
   std::vector<double> lastR;
   std::vector<char> swClosed;
   std::vector<std::shared_ptr<StopBody>> stopBodies;   // what the obstacles point at
   // the plate's laws
   double PhiStop() const { return sl.mech.stopDeg * FlipperModel::DEG; }
   double Inertia(const double ph) const
   {
      const FlipperModel::Kin k = FlipperModel::Kinematics(sl.mech, der, der.theta0 + ph);
      return der.iFlipper + sl.mech.plungerKg * (k.lever * 1e-3) * (k.lever * 1e-3);
   }
   double PlateV(const double ph) const
   {
      const FlipperSettings &s = sl.mech;
      const double u = der.xRest - FlipperModel::SliderX(s, der, der.theta0 + ph);
      return 1e-3 * (s.preloadN * u + 0.5 * der.kSpringNPerMm * u * u);
   }
   double PlateVp(const double ph) const { const double e = 1e-7; return (PlateV(ph + e) - PlateV(ph - e)) / (2.0 * e); }
   double PlateVpp(const double ph) const { const double e = 1e-5; return (PlateV(ph + e) - 2.0 * PlateV(ph) + PlateV(ph - e)) / (e * e); }
   void RestStop(const double ph, double &tauS, double &cD) const
   {
      (void)ph;
      tauS = cD = 0.0;   // the rest stop is the hard seat
   }
   // The coil's torque on the plate at the current now, N m
   double CoilTorqueAt(const double ph) const
   {
      const FlipperModel::Kin k = FlipperModel::Kinematics(sl.mech, der, der.theta0 + ph);
      return CoilModel::Force(coilS, coilD, (k.x - der.xMin) / coilS.bobbinMm, iCoil) * k.lever * 1e-3;
   }
   // The circuit and the force over a trial move of the plate: the plunger's move through the link, from the step's start, so a retake is identical (CoilModel::Over, the coil's only energy law)
   CoilModel::OverMove CoilOver(const double dphi, const double h) const
   {
      const double x1 = FlipperModel::SliderX(sl.mech, der, der.theta0 + phi + dphi);
      return CoilModel::Over(coilS, coilD, (gKin0.x - der.xMin) / coilS.bobbinMm, (x1 - der.xMin) / coilS.bobbinMm, (gKin0.x - x1) * 1e-3, gVolts, iCoil, h);
   }
   // The plunger's force as a torque on the plate over the same move: F dx / dphi, the link's chord, so the plate is given exactly the plunger's work. No move: the lever there
   double CoilTorque(const CoilModel::OverMove &o, const double dphi) const
   {
      if (std::fabs(dphi) > 1e-7)
         return o.force * (gKin0.x - FlipperModel::SliderX(sl.mech, der, der.theta0 + phi + dphi)) * 1e-3 / dphi;
      return o.force * FlipperModel::LeverMm(sl.mech, der, der.theta0 + phi + 0.5 * dphi) * 1e-3;
   }
   double RestPhi() const
   {
      double lo = -0.2, hi = 0.0;
      for (int it = 0; it < 40; it++)
      {
         const double mid = 0.5 * (lo + hi);
         (PlateVp(mid) > 0.0 ? hi : lo) = mid;
      }
      return 0.5 * (lo + hi);
   }
   // where the cord's centre must be clear of the plate's face at the cord's height, along the push, from the rest line; and its rate with the angle
   void Limit(const double ph, double &y, double &dy) const
   {
      const double t = ph + pl.rest, ct = std::cos(t), st = std::sin(t);
      y = pl.originN + pl.a / ct + pl.h * st / ct + rc / ct;
      dy = (pl.a * st + pl.h) / (ct * ct) + rc * st / (ct * ct);
   }
   // The cord as the plate has bent it, with no ball: straight to the plate's edges from the posts, flat along its face. A function of the span's station and the plate's angle only, so no state and no choice. It follows the face with a bend rounded by the bending length, about the depth the cord sinks into the plate
   struct Base { double y = 0.0, ya = 0.0, yp = 0.0; };
   Base BaseY(const Span &sp, const double a) const
   {
      Base b;
      if (sp.bent)
         for (const Obstacle &ob : sp.ob)
            if (ob.bend != 0.0)   // straight from each post to the body, the same rounding as the plate's edges
            {
               const double p = a / ob.x, q = (sp.len - a) / (sp.len - ob.x), e1 = w / ob.x, e2 = w / (sp.len - ob.x);
               const double r1 = std::sqrt((p - 1.0) * (p - 1.0) + e1 * e1), s1 = 0.5 * (p + 1.0 - r1), ds1 = 0.5 * (1.0 - (p - 1.0) / r1) / ob.x;
               const double r2 = std::sqrt((s1 - q) * (s1 - q) + e2 * e2), s2 = 0.5 * (s1 + q - r2), dq = -1.0 / (sp.len - ob.x);
               b.y += ob.bend * s2;
               b.ya += ob.bend * 0.5 * (ds1 + dq - (s1 - q) * (ds1 - dq) / r2);
            }
      if (!pl.on || &sp != &span[pl.span])
         return b;
      double yl, dy;
      Limit(phiQ, yl, dy);
      if (yl <= 0.0)
         return b;
      const double yr = yl < w ? yl * yl / (2.0 * w) : yl - 0.5 * w, dyr = yl < w ? yl / w : 1.0;
      const double p = a / pl.x0, q = (sp.len - a) / (sp.len - pl.x1), e1 = w / pl.x0, e2 = w / (sp.len - pl.x1);
      // the least of p, 1, q, smoothly
      const double r1 = std::sqrt((p - 1.0) * (p - 1.0) + e1 * e1), s1 = 0.5 * (p + 1.0 - r1), ds1 = 0.5 * (1.0 - (p - 1.0) / r1) / pl.x0;
      const double r2 = std::sqrt((s1 - q) * (s1 - q) + e2 * e2), s2 = 0.5 * (s1 + q - r2);
      const double dq = -1.0 / (sp.len - pl.x1), ds2 = 0.5 * (ds1 + dq - (s1 - q) * (ds1 - dq) / r2);
      b.y += pl.sgn * yr * s2;
      b.ya += pl.sgn * yr * ds2;
      b.yp = pl.sgn * dyr * dy * s2;
      return b;
   }

   // A run blocked by a body: what it shares between lift and squash
   struct Edge
   {
      const GroupStep::Item *law = nullptr;
      double gr = 0.0, g = 0.0, ell = 0.0, a1 = 0.0, a2 = 0.0;   // the gap to share; the full gap; how far the run has left the ball there; the run's lengths either side of the body
      V3 dgr, da1;
      double dgrP = 0.0, da1P = 0.0;   // their rates with the plate's angle
      double x = 0.0;
      mutable double xLast = -1.0;   // this Solve's last root: the next tension's Newton starts here (Solve clears it, so nothing is kept between calls)
   };
   // One ball on one span
   struct Out
   {
      double L = 0.0, D = 0.0, Deff = 0.0, a = 0.0, qLen = 0.0, gap = 0.0, Rc = 0.0, Gfree = 0.0;
      bool pastLine = false;
      V3 G, n, ohat;
      // what the squash needs: the same wrap at a smaller radius
      double La = 0.0, Ld = 0.0, LR = 0.0, sg = 1.0;   // the free length's slopes along the span, across it, with the circle's radius; which side the ball presses from
      bool shadow = false;   // the taut path round the other balls is clear of this one: no force, whatever its own dent
      double Rr = 0.0, dt = 0.0, gradA = 0.0, cenP = 0.0, LP = 0.0;   // cenP: how the ball's offset from the cord moves with the plate's angle; LP: how the free length does
      V3 gC, gR;
      mutable double sq = 0.0;   // the squash Solve found
      int nEd = 0;
      Edge ed[kMaxEdges];
   };
   // The free length of one ball on one span, and the bodies its runs are blocked by. Everything is in the plane: a = along, y = across, dz = off the plane
   void Length(const Span &sp, const Touch &u, const V3 &c, const double R, const bool wantG, Out &o) const
   {
      const double Rc = R + rc;
      const V3 w0 = c - sp.A;
      const double aRaw = w0.Dot(sp.t), a = std::min(std::max(aRaw, 0.0), sp.len), over = aRaw - a, dz = w0.Dot(sp.up);
      const Base bs = BaseY(sp, a);
      const double y = w0.Dot(sp.yh) - bs.y;   // measured from the cord as the plate has it: a ball touching it stores nothing, wherever the plate is
      const bool inA = aRaw > 0.0 && aRaw < sp.len;
      o.cenP = o.LP = 0.0;
      const double sgn = u.on ? u.nh.Dot(sp.yh) : (y >= 0.0 ? 1.0 : -1.0);
      o.a = a;
      o.Rc = Rc;
      o.L = sp.len;
      o.nEd = 0;
      o.pastLine = false;
      o.La = o.Ld = o.LR = 0.0;
      o.sg = sgn >= 0.0 ? 1.0 : -1.0;
      o.shadow = false;
      o.ohat = sp.yh * (y >= 0.0 ? 1.0 : -1.0);
      // the circle the ball cuts in the plane. From the side: below the plane the sphere's own, which presses the ball down; above it, no smaller: the tension never lifts the ball. Coming down: above the plane the sphere's own, which holds it up
      const double Rs2 = std::max(Rc * Rc - dz * dz, 1e-12), Rs = std::sqrt(Rs2);
      const bool sph = u.top ? dz >= 0.0 : dz < 0.0;
      const double Rr = sph ? Rs : Rc, Rp = sph ? -dz / Rs : 0.0;
      const V3 gR = sp.up * Rp;
      // where its centre is across the line: from the side, signed by its side, so past the line the length keeps growing; coming down, the line rounded by the bending length, either side
      double cen;
      V3 gC;
      if (u.top)
      {
         const double e1 = std::sqrt(Rc * Rc + w * w) - Rc, r = std::sqrt(y * y + over * over + w * w);
         cen = r - e1;
         gC = (sp.yh * y + sp.t * over + sp.t * (inA ? -y * bs.ya : 0.0)) * (1.0 / r);
         o.cenP = -y * bs.yp / r;
      }
      else
      {
         const double ys = sgn * y;
         if (over == 0.0)
         {
            cen = ys;
            gC = sp.yh * sgn + sp.t * (inA ? -sgn * bs.ya : 0.0);
            o.cenP = -sgn * bs.yp;
         }
         else
         {
            const double yy = std::sqrt(ys * ys + over * over), sq = ys >= 0.0 ? 1.0 : -1.0;
            cen = sq * yy;
            gC = (sp.yh * (sgn * ys) + sp.t * over) * (sq / std::max(yy, 1e-12));
         }
         o.pastLine = u.on && ys < 0.0;
      }
      o.D = Rr - cen;
      o.Deff = o.D - u.d0;
      o.gap = -o.Deff;
      // the cord's place from the ball's centre
      const V3 q = sp.yh * (-sgn * Rr) + sp.up * -dz;
      o.qLen = q.Length();
      o.n = o.qLen > 1e-12 ? q * (1.0 / o.qLen) : sp.up * -1.0;
      o.sq = 0.0;
      if (o.Deff <= 0.0)
      {
         o.G = V3();
         o.Gfree = 0.0;
         return;
      }
      const double dt = cen + u.d0, b = sp.len - a, gradA = (aRaw > 0.0 && aRaw < sp.len) ? 1.0 : 0.0;
      o.Rr = Rr;
      o.dt = dt;
      o.gradA = gradA;
      o.gC = gC;
      o.gR = gR;
      SideLen SL, SR;
      SL.Set(a, dt, Rr);
      SR.Set(b, dt, Rr);
      o.L = SL.val + SR.val;
      o.La = SL.Sb - SR.Sb;
      o.Ld = SL.Sd + SR.Sd;
      o.LR = SL.beta + SR.beta;
      if (wantG)
      {
         o.G = sp.t * ((SL.Sb - SR.Sb) * gradA) + gC * (SL.Sd + SR.Sd) + gR * (SL.beta + SR.beta);
         o.Gfree = o.G.Length();
         o.LP = (SL.Sd + SR.Sd) * o.cenP;
      }
      const V3 ec = sp.yh * -sgn;   // the way the cord is pressed, in the plane
      for (const Obstacle &ob : sp.ob)
      {
         // where the ball's run is at the body's station: only a straight run can be blocked; under the ball's arc the cord lies on the ball, and that is the Item's squeeze
         const bool right = ob.x >= a;
         const SideLen &S = right ? SR : SL;
         if (S.s <= 1e-7)
            continue;
         const double sb = std::sin(S.beta), cb = std::cos(S.beta), Tp = right ? a + Rr * sb : a - Rr * sb, a1 = right ? ob.x - Tp : Tp - ob.x;
         if (a1 <= 0.0 || cb <= 1e-6)
            continue;
         const double far2 = right ? sp.len - ob.x : ob.x, hh = far2 * sb / cb;
         // how far the ball's run is behind the body at its station
         double raw, yl = 0.0, dyl = 0.0;
         V3 mv;
         if (ob.body)
         {
            double g;
            ob.body->Gap(sp.A + sp.t * ob.x + sp.yh * ob.bend + ec * hh, rc, g, mv);
            raw = -g;
         }
         else
         {
            if (sgn * pl.sgn < 0.0)   // the ball is behind the cord: it presses the cord away from the plate
               continue;
            Limit(phiQ, yl, dyl);
            raw = hh + std::min(yl, 0.0);
         }
         if (raw <= 0.0 || o.nEd >= kMaxEdges)
            continue;
         // the run can lift only until it lies on the ball: a straight run has left the ball by a1^2 / 2 R this far from where it let go. The smaller of the two, smoothly, zero at either's zero
         const double ell = a1 * a1 / (2.0 * Rr), r2g = std::sqrt(raw * raw + ell * ell);
         Edge &e = o.ed[o.nEd++];
         e.law = &ob.law;
         e.g = -raw;
         e.ell = ell;
         e.a1 = a1;
         e.a2 = far2;
         e.x = ob.x;
         e.gr = raw * ell / r2g;
         if (!wantG)
            continue;
         const double bx = right ? b : a, D2 = std::max(bx * bx + dt * dt, 1e-18);
         const double dBetaB = -Rr * bx / (D2 * S.s) + dt / D2, dBetaD = -Rr * dt / (D2 * S.s) - bx / D2, sec2 = 1.0 / (cb * cb);
         const V3 gA = sp.t * gradA;
         const V3 gBeta = gA * (right ? -dBetaB : dBetaB) + gC * dBetaD + gR * (1.0 / S.s);
         const double gBetaP = dBetaD * o.cenP;
         const V3 gH = gBeta * (far2 * sec2);                          // grad of the run's depth at the station, through its lever: the rest goes to the far post
         const double hP = gBetaP * far2 * sec2;
         const V3 gRaw = ob.body ? gH * -mv.Dot(ec) : gH;
         const double rawP = ob.body ? 0.0 : hP + (yl < 0.0 ? dyl : 0.0);
         const double sgnT = right ? 1.0 : -1.0;
         const V3 gTp = gA + (gR * sb + gBeta * (Rr * cb)) * sgnT;
         const double gTpP = gBetaP * Rr * cb * sgnT;
         e.da1 = right ? gTp * -1.0 : gTp;
         e.da1P = right ? -gTpP : gTpP;
         const V3 gEll = e.da1 * (a1 / Rr) + gR * (-a1 * a1 / (2.0 * Rr * Rr));
         const double ellP = e.da1P * a1 / Rr, r3 = r2g * r2g * r2g, wR = ell * ell * ell / r3, wE = raw * raw * raw / r3;
         e.dgr = gRaw * wR + gEll * wE;
         e.dgrP = rawP * wR + ellP * wE;
      }
   }
   // Two or more balls pressing a span from one side: the cord's free length is the one shortest path from notch to notch with their circles on its one side, closed form: tangents and arcs. Each run's angle from the same rule as SideLen; a ball whose arc would turn the wrong way is clear of the path and takes no part (a lower hull, ordered along the span).
   // V stays one function of the positions, and C1: a ball joins or leaves with its arc at zero, where its two runs are one straight line, so the length is the same either way and its force (the difference of the runs' directions, and the arc for the radius) starts from zero. With one ball on a side this is the single-ball form exactly. Balls on opposite sides add, group against group.
   template <class Get> void Chain(const Span &sp, const Get &get, const int n, const bool wantG) const
   {
      if (!chain || n < 2)
         return;
      struct Vx { double x, y, r, thIn, sIn; int i; };
      thread_local std::vector<int> id;
      thread_local std::vector<Vx> st_;
      for (const double side : { 1.0, -1.0 })
      {
         id.clear();
         for (int i = 0; i < n; i++)
            if (get(i).Deff > 0.0 && get(i).sg == side)
               id.push_back(i);
         if (id.size() < 2)
            continue;
         std::sort(id.begin(), id.end(), [&](const int p, const int q) { return get(p).a < get(q).a; });
         int endJ = -1, startJ = -1;
         double endE[2] = { 0.0, 0.0 }, startE[2] = { 0.0, 0.0 }, endRoD = 0.0, startRoD = 0.0, baseX = 0.0, baseY = 0.0, endX = sp.len, endY = 0.0;
         bool redo = false;
         const auto Push = [&](const double xx, const double yy, const double rr, const int who) {
            for (;;)
            {
               const Vx &p = st_.back();
               const double dx = xx - p.x, dy = yy - p.y, dr = rr - p.r, X = dx * dx + dy * dy - dr * dr, s = X > 0.0 ? std::sqrt(X) : 0.0, th = std::atan2(dy, dx) - std::atan2(dr, s);
               if (st_.size() >= 2 && th < p.thIn)
               {
                  // a ball clear of the run. If a notch is inside its circle (it sits on the post), the single-ball form ends the cord at the notch's radial point on the circle, so the run ends there instead: the same length where the ball would have left the path, so V doesn't jump. The path's end has moved, so what was clear of the old run is looked at again
                  const auto Radial = [&](const double ex, const double ey, double &px, double &py, double &D) {
                     D = std::max(std::hypot(ex - p.x, ey - p.y), 1e-12);
                     px = p.x + p.r * (ex - p.x) / D;
                     py = p.y + p.r * (ey - p.y) / D;
                     return D < p.r;
                  };
                  double px, py, D;
                  if (who == -1 && endJ < 0 && Radial(sp.len, 0.0, px, py, D))
                  {
                     endJ = p.i;
                     endE[0] = (sp.len - p.x) / D;
                     endE[1] = -p.y / D;
                     endRoD = p.r / D;
                     endX = px;
                     endY = py;
                     redo = true;
                     return;
                  }
                  if (st_.size() == 2 && startJ < 0 && Radial(0.0, 0.0, px, py, D))
                  {
                     startJ = p.i;
                     startE[0] = -p.x / D;
                     startE[1] = -p.y / D;
                     startRoD = p.r / D;
                     baseX = px;
                     baseY = py;
                     redo = true;
                     return;
                  }
                  st_.pop_back();
                  continue;
               }
               st_.push_back({ xx, yy, rr, th, s, who });
               return;
            }
         };
         for (const int i : id)
         {
            Out &o = get(i);
            o.L = sp.len;
            o.La = o.Ld = o.LR = o.Gfree = o.LP = 0.0;
            o.G = V3();
            o.shadow = true;
         }
         do
         {
            redo = false;
            st_.clear();
            st_.push_back({ baseX, baseY, 0.0, 0.0, 0.0, -1 });
            for (const int i : id)
            {
               if (i == endJ || i == startJ)
                  continue;
               const Out &o = get(i);
               Push(o.a, o.dt, o.Rr, i);
               if (redo)
                  break;
            }
            if (!redo)
               Push(endX, endY, 0.0, endJ >= 0 ? -2 : -1);
         } while (redo);
         const auto ArcSlope = [](const Vx &v, const double ex, const double ey, const double sg, double &ux, double &uy) {
            const double dx = ex - v.x, dy = ey - v.y, D2 = std::max(dx * dx + dy * dy, 1e-18);
            ux = -sg * v.r * dy / D2;
            uy = sg * v.r * dx / D2;
         };
         double L = st_.back().sIn;
         for (size_t k = 1; k + 1 < st_.size(); k++)
         {
            const Vx &v = st_[k];
            const double in = v.thIn, out = st_[k + 1].thIn;
            L += v.sIn + v.r * (out - in);
            Out &o = get(v.i);
            o.shadow = false;
            // the runs' directions; with the path's end inside the circle, no run: the slope of the arc to its radial point, as SideLen's
            double ci = std::cos(in), si = std::sin(in), co = std::cos(out), so = std::sin(out);
            if (k == 1 && v.sIn <= 0.0)
               ArcSlope(v, st_[0].x, st_[0].y, 1.0, ci, si);
            if (k + 2 == st_.size() && st_[k + 1].sIn <= 0.0)
               ArcSlope(v, st_[k + 1].x, st_[k + 1].y, 1.0, co, so);
            o.La = ci - co;
            o.Ld = si - so;
            o.LR = out - in;
            if (wantG)
            {
               o.G = sp.t * (o.La * o.gradA) + o.gC * o.Ld + o.gR * o.LR;
               o.Gfree = o.G.Length();
               o.LP = o.Ld * o.cenP;
            }
         }
         // a ball whose notch point the run ends on: the run's pull, through how that point moves with the ball
         const auto OnNotch = [&](const int j, const double ux, const double uy, const double *e, const double RoD, const double sg) {
            Out &o = get(j);
            const double eu = e[0] * ux + e[1] * uy;
            o.shadow = false;
            o.La = sg * (ux - RoD * (ux - e[0] * eu));
            o.Ld = sg * (uy - RoD * (uy - e[1] * eu));
            o.LR = sg * eu;
            if (wantG)
            {
               o.G = sp.t * (o.La * o.gradA) + o.gC * o.Ld + o.gR * o.LR;
               o.Gfree = o.G.Length();
               o.LP = o.Ld * o.cenP;
            }
         };
         // how the length moves with that point: the run's direction, or if the point is inside the next ball's circle too (two balls, one each side of a notch), that ball's arc
         if (endJ >= 0)
         {
            double ux = std::cos(st_.back().thIn), uy = std::sin(st_.back().thIn);
            if (st_.back().sIn <= 0.0 && st_.size() >= 3)
               ArcSlope(st_[st_.size() - 2], st_.back().x, st_.back().y, 1.0, ux, uy);
            OnNotch(endJ, ux, uy, endE, endRoD, 1.0);
         }
         if (startJ >= 0)
         {
            double ux = std::cos(st_[1].thIn), uy = std::sin(st_[1].thIn);
            if (st_[1].sIn <= 0.0 && st_.size() >= 3)
               ArcSlope(st_[1], st_[0].x, st_[0].y, 1.0, ux, uy);
            OnNotch(startJ, ux, uy, startE, startRoD, -1.0);
         }
         // the group's whole length goes on its first ball; only the sum is ever read
         get(id[0]).L = L;
      }
   }
   // The cord's force on each ball now, over all its spans, and the highest tension (diagnostics: nothing reads them)
   void CordForce(const Ball *b, std::vector<V3> &F, double &Tmax) const
   {
      phiQ = phi;
      thread_local std::vector<Out> oo;
      thread_local std::vector<V3> gg;
      oo.resize(nb);
      gg.resize(nb);
      F.assign(nb, V3());
      Tmax = 0.0;
      for (int j = 0; j < ns; j++)
      {
         for (int i = 0; i < nb; i++)
            Length(span[j], tch[(size_t)i * ns + j], b[i].c, b[i].R, true, oo[i]);
         Chain(span[j], [&](const int i) -> Out & { return oo[i]; }, nb, true);
         double L, T;
         Solve(span[j], [&](const int i) -> const Out & { return oo[i]; }, nb, L, T, gg.data());
         for (int i = 0; i < nb; i++)
            F[i] = F[i] - gg[i];
         Tmax = std::max(Tmax, T);
      }
   }
   // What a run gains lifted by e off its straight line between two points a1, a2 away, without the cancellation
   static double Lift(const double e, const double a1, const double a2)
   {
      const double r1 = std::sqrt(a1 * a1 + e * e), r2 = std::sqrt(a2 * a2 + e * e);
      return e * e / (r1 + a1) + e * e / (r2 + a2);
   }
   // The whole span at these balls: V, the length, the tension. At each ball the wrap's pull meets the shared law's squash of the ball on the cord, and at each body the lifted run's pull meets its squash: massless series springs, solved out, so V is one function of the positions
   template <class Get> double Solve(const Span &sp, const Get &get, const int n, double &Ltot, double &T, V3 *grad) const
   {
      double Lf = sp.len;
      int on = 0;
      for (int i = 0; i < n; i++)
      {
         const Out &o = get(i);
         Lf += o.L - sp.len;
         on += o.Deff > 0.0 ? 1 : 0;
      }
      T = std::max(sp.k * (Lf - sp.l0), 0.0);
      Ltot = Lf;
      double Vsq = 0.0;
      // the plate pressing the cord, ball or no ball: past the line, each edge lifts the run from it to its post, and the squash takes its share (the same rule as any body on a run). The face between them adds no length
      Edge base[2 + kMaxEdges];
      int nBase = 0, nPlate = 0;
      double dyl = 0.0;
      if (pl.on && &sp == &span[pl.span])
      {
         double yl;
         Limit(phiQ, yl, dyl);
         if (yl > 0.0)
         {
            base[0].law = base[1].law = &pl.law;
            base[0].gr = base[1].gr = yl;
            base[0].a1 = pl.x0;
            base[1].a1 = sp.len - pl.x1;
            base[0].a2 = base[1].a2 = 1e3;
            nBase = nPlate = 2;
         }
      }
      if (sp.bent)   // a body holding the cord at rest: its bend is stored whether a ball is on the cord or not, so a ball's touch changes nothing
         for (const Obstacle &ob : sp.ob)
            if (ob.bend != 0.0 && nBase < 2 + kMaxEdges)
            {
               Edge &e = base[nBase++];
               e.law = &ob.law;
               e.gr = std::fabs(ob.bend);
               e.a1 = ob.x;
               e.a2 = sp.len - ob.x;
            }
      for (int i = 0; i < n; i++)
      {
         const Out &o = get(i);
         for (int k = 0; k < o.nEd; k++)
            o.ed[k].xLast = -1.0;
      }
      if (on > 0 || nBase > 0)
      {
         for (int it = 0; it < 30; it++)
         {
            double add = sp.len;
            Vsq = 0.0;
            for (int k = 0; k < nBase; k++)
            {
               const double lift = Share(base[k], T);
               add += Lift(lift, base[k].a1, base[k].a2);
               Vsq += base[k].law->V(base[k].gr - lift);
            }
            for (int i = 0; i < n; i++)
            {
               const Out &o = get(i);
               if (o.Deff <= 0.0)
                  continue;
               add += o.L - sp.len;
               for (int k = 0; k < o.nEd; k++)
               {
                  const Edge &e = o.ed[k];
                  const double lift = Share(e, T);
                  add += Lift(lift, e.a1, e.a2);
                  Vsq += e.law->V(e.gr - lift);
               }
            }
            Ltot = add;
            const double Tn = std::max(sp.k * (Ltot - sp.l0), 0.0), dT = std::fabs(Tn - T);
            T = Tn;
            if (dT < 1e-11)
               break;
         }
      }
      if (grad)
         dVdP = 0.0;
      if (grad)
         for (int k = 0; k < nPlate; k++)
            dVdP += pl.law.F(base[k].gr - Share(base[k], T)) * dyl;
      if (grad)
         for (int i = 0; i < n; i++)
         {
            const Out &o = get(i);
            grad[i] = V3();
            if (o.Deff <= 0.0)
               continue;
            grad[i] = (sp.t * (o.La * o.gradA) + o.gC * o.Ld + o.gR * o.LR) * T;
            dVdP += T * o.Ld * o.cenP;
            for (int k = 0; k < o.nEd; k++)
            {
               const Edge &e = o.ed[k];
               const double lift = Share(e, T), r1 = std::sqrt(e.a1 * e.a1 + lift * lift);
               const double fe = e.law->F(e.gr - lift), ta = T * -(lift * lift) / (r1 * (r1 + e.a1));
               grad[i] = grad[i] + e.dgr * fe + e.da1 * ta;
               dVdP += e.dgrP * fe + e.da1P * ta;
            }
         }
      return Es(sp, Ltot) + Vsq;
   }
   // The lift at which the run's pull meets the squash (both monotonic: one root; a safeguarded Newton from the same start each time, no state)
   static double Share(const Edge &e, const double T)
   {
      const double hL = e.law->hL;
      double lo = hL > 0.0 ? std::max(0.0, e.gr - hL * (1.0 - 1e-9)) : 0.0, hi = e.gr, x = e.xLast > lo && e.xLast < hi ? e.xLast : 0.5 * (lo + hi);
      for (int it = 0; it < 60; it++)
      {
         const double r1 = std::sqrt(e.a1 * e.a1 + x * x), r2 = std::sqrt(e.a2 * e.a2 + x * x), d = e.gr - x;
         const double f = T * (x / r1 + x / r2) - e.law->F(d);
         if (f > 0.0)
            hi = x;
         else
            lo = x;
         const double fp = T * (e.a1 * e.a1 / (r1 * r1 * r1) + e.a2 * e.a2 / (r2 * r2 * r2)) + e.law->K(d);
         double nx = fp > 0.0 ? x - f / fp : 0.5 * (lo + hi);
         // solved before the bracket is looked at: at the root a Newton move can land on the bracket's own end, and that isn't a reason to start again from its middle
         if (std::fabs(nx - x) < 1e-14 && nx >= lo && nx <= hi)
         {
            x = nx;
            break;
         }
         if (!(nx > lo && nx < hi))
            nx = 0.5 * (lo + hi);
         const bool done = std::fabs(nx - x) < 1e-14;
         x = nx;
         if (done)
            break;
      }
      e.xLast = x;
      return x;
   }
   static double Es(const Span &sp, const double L)
   {
      const double e = L - sp.l0;
      return e > 0.0 ? 0.5 * sp.k * e * e : 0.0;
   }

private:
   static constexpr double kPi_ = 3.14159265358979323846;
   // A grip axis inside a step: its limit and the start held, Coulomb's elastic-plastic energy; the chord's force over the slip S; the energies at both ends
   static void Chord(const double k, const double f, const double xi0, const double S, double &fc, double &xi1, double &p0, double &p1)
   {
      if (k <= 0.0 || f <= 0.0)
      {
         fc = 0.0;
         xi1 = 0.0;
         p0 = p1 = 0.0;
         return;
      }
      const double lim = f / k, x0 = std::min(std::max(xi0, -lim), lim);
      xi1 = x0 + S;
      const auto Psi = [&](const double z) {
         const double az = std::fabs(z);
         return az <= lim ? 0.5 * k * z * z : f * az - 0.5 * f * lim;
      };
      p0 = 0.5 * k * x0 * x0;
      p1 = Psi(xi1);
      if (std::fabs(S) > 1e-14)
         fc = (p1 - p0) / S;
      else
         fc = k * x0;
   }
   struct Scr
   {
      Out o0, o1, om;
      double lossS = 0.0, lossTw = 0.0;
      double kT = 0.0, kB = 0.0, xiT = 0.0, xiB = 0.0, xiT1 = 0.0, xiB1 = 0.0, f = 0.0, lossK = 0.0, pT0 = 0.0, pT1 = 0.0, pB0 = 0.0, pB1 = 0.0;
   };
   static constexpr double kFarM = 5e-3, kFarMove = 4e-3;   // a span is skipped when clear of every ball by this, and while no ball's move this step passes the smaller figure (the rounding of a coming-down ball's line is 0.3 mm)
   std::vector<char> farS;
   bool startReady = false;   // StepFacts has left each ball's length at the positions now in st[].o0
   std::vector<double> startV, startT;
   std::vector<Scr> st;
   std::vector<double> V0S, dLS, TdS, TmS;
   std::vector<V3> gv;
   std::vector<Touch> save;
   double savePl[3] = { 0, 0, 0 };
   bool savePin = false, saveRest = false, gFreed = false;
   double gTauPinned = 0.0;
   static constexpr double kArm = 0.05;   // the distance that makes the plate's angle a length, in the solver's rows and in sharing the closure
   double gLoad = 0.0, gRestL = 0.0, gPivotL = 0.0, gCordTau = 0.0;
   double gI0 = 1.0, gp0 = 0.0, gVolts = 0.0, gAlpha = 0.0, gPlateLoss = 0.0, gCoilW = 0.0, gPinJ = 0.0;
   bool gCoilLive = false;   // volts on or current flowing: otherwise the coil costs nothing
   FlipperModel::Kin gKin0{ 0.0, 0.0, 0.0 };   // the link at the step's start
   CoilModel::OverMove gOver;
   const Ball *b0 = nullptr;
};
// What the ball meets on the plate: the cord lying on the face at the cord's height, a straight tube carried by the plate. The ring's Item, the same code; the group books its reaction on the plate
class PlateTube : public GroupStep::Shape
{
public:
   const Part *pt = nullptr;
   void Probe(const V3 &c, const double R, GroupStep::Hit &h) const override
   {
      const Span &sp = pt->span[pt->pl.span];
      double yl, dy;
      pt->Limit(pt->phi + (trial ? trial[0] : 0.0), yl, dy);
      const V3 up = sp.yh * pt->pl.sgn, p0 = sp.A + sp.t * pt->pl.x0 + up * yl, dp = sp.t * (pt->pl.x1 - pt->pl.x0);
      const double ll = dp.Dot(dp), u = std::min(std::max((c - p0).Dot(dp) / ll, 0.0), 1.0);
      const V3 d = c - (p0 + dp * u);
      const double dist = std::max(d.Length(), 1e-12), Rc = R + pt->rc, R1 = R * pt->rc / (R + pt->rc);
      h.delta = Rc - dist;
      h.n = d * (1.0 / dist);
      h.reff = std::sqrt(R1 * R);
      h.curvDist = 0.0;
      h.grad = 1.0;
      h.nc = 1;
      const V3 vq = up * dy;
      h.c[0].Set(0, h.n.Dot(vq), vq, V3());
   }
};
}
