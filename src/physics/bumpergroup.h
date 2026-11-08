// license:GPLv3+
//
// The jet bumper as a part of the group step (runs/2026-10-31-bumper-design/design.md sections 5, 6, 7: skirt, ring and coil, in its three drive modes). Geometry, masses and materials are contactbumper.h's; the moving is the group's.
// Coordinates: the skirt's height z and its two tilts thx, thy (R = Rx(thx) Ry(thy)), always unknowns; the ring's drop dr (down positive; with its rods, yoke and plunger), an unknown except on a stop. The spoon has none: a massless spring on the follower distance D(q).
// The ring's stops are impacts (step ends on the stop, e from the settings, the e(v) table is open), and a rebound shorter than the stop's own contact time never left it. Its coil is the only drive: the circuit and its force over the solved move, in the solve (CoilModel::Over), so what the coil gives the ring is what the circuit supplied, at any step length.
// Stored energy, each in closed form with its own discrete gradient, so the step keeps its energy at any length: gravity from the shared g, the conical spring, the seat at the lip, the seat at the rim on the playfield, the spoon.
// Losses, each only removing energy: the seats' and the spoon's dampers act on the part of the move inside them and push only; the skirt's viscous loss on its mid-step rates.
// A seat's damper works on its equivalent depth s, V = 1/2 k s^2: the same form as a ball's contact (N = max(0, fe + fd), fe = (V1 - V0) / (s1 - s0), fd = c (s1 - s0) / h). Its force along the discrete gradient of V is that gradient times 1 + fd / fe. Pure SI.

#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "physics/contactbumper.h"
#include "physics/groupstep.h"

namespace BumperGroup
{
   using GroupStep::V3;
   using GroupStep::Ball;
   using GroupStep::Hit;
   namespace JB = JetBumper;

   // The effective radius of a ball on a surface of revolution
   inline double KappaEff(const double R, const JB::Revolve::Near &nr)
   {
      double inv = 1.0 / R;
      if (nr.edge)
         inv += 1.0 / nr.rr;
      if (nr.rs > 1e-4)
         inv += nr.nr / nr.rs;
      return 1.0 / std::max(inv, 1.0 / (50.0 * R));
   }

   struct PartBooks
   {
      double lipJ = 0.0, floorJ = 0.0, spoonJ = 0.0, springJ = 0.0, viscousJ = 0.0;
      double ringStopJ = 0.0, ringViscousJ = 0.0;                 // the ring's losses
      double coilJ = 0.0, elecJ = 0.0, resistJ = 0.0;              // coil: work on the ring (the drive), electrical energy delivered, lost in its resistance
      long restHits = 0, coilHits = 0;                             // ring impacts
      double tFirstCoilHit = -1.0, tLastCoilHit = -1.0, tFirstRestHit = -1.0, tLastRestHit = -1.0, vFirstCoilHit = 0.0, vFirstRestHit = 0.0;
      double tPulseStart = -1.0, tPulseEnd = -1.0;
      long pulses = 0;
   };

   class Part : public GroupStep::Part
   {
   public:
      static constexpr int kZ = 0, kThx = 1, kThy = 2, kR = 3;   // the coordinates
      JB::Bumper cfg;           // its geometry, masses and materials
      PartBooks bk;
      // state (SI)
      double z = 0.0, thx = 0.0, thy = 0.0, vz = 0.0, wx = 0.0, wy = 0.0;
      double t = 0.0;           // s, the time the part has been stepped
      bool closed = false;      // the switch
      int closures = 0;
      double tClose = -1.0;     // the first closure
      double reach = 0.05;      // m, the lever that turns a tilt rate into a speed
      double dr = 0.0, vr = 0.0;   // the ring's drop and rate, down positive
      int ringStop = -1;        // -1 on the rest stop (no unknown), +1 on the coil stop (no unknown), 0 free
      double iCoil = 0.0;       // A
      bool coilEnabled = true, pulseActive = false, pulseArmed = true;
      double pulseLeft = 0.0;
      bool extOn = false;       // external drive: the script's command
      double extVolts = -1.0, extDuty = 1.0;
      double coilOnS = 0.0, coilOnMaxS = 0.0;   // how long the coil has been commanded on, this episode and the longest (a detector's)
      std::vector<char> swLog;  // the switch's changes in order, until the owner reads them: 1 closed, 0 opened
      struct CoilStepBook { bool free = false; double dx = 0.0, v0 = 0.0, v1 = 0.0, i0 = 0.0, i1 = 0.0, h = 0.0, load = 0.0, elec = 0.0, res = 0.0, dField = 0.0, work = 0.0; };
      CoilStepBook cb;          // the last step's coil balance (a diagnostic)

      // The settings. coneRoundMm: the spoon cone's rounded tip (estimated by eye, design section 5). gs: the group's settings, for the shared g
      void Configure(const JB::Params &prm, const GroupStep::Settings &gs, const double coneRoundMm = 0.25)
      {
         cfg.Configure(prm);
         set = &gs;
         const JB::Params &p = cfg.p;
         m = p.skirtMass;
         I = cfg.mThetaEff;
         cmZ = cfg.cmZ;
         tIn = p.skirtTInnerMm * JB::MM;
         zRest = cfg.zRest;
         zBase = cfg.zBase;
         rLip = 0.5 * p.lipDiaMm * JB::MM;
         rRim = 0.5 * p.skirtOdMm * JB::MM - 0.5 * p.skirtTEdgeMm * JB::MM;
         Lt = (p.stemLenMm - p.tipRadiusMm) * JB::MM;
         kSpring = p.skirtKNPerM;
         fSpring = p.skirtPreloadN;
         aSpring = p.springRadiusMm * JB::MM;
         kSeat = p.skirtStopKNPerM;
         cSeat = 2.0 * p.skirtStopZeta * std::sqrt(kSeat * m);
         kSpoon = p.spoonKNPerM;
         fSpoon = p.spoonPreloadN;
         cSpoon = 2.0 * p.spoonZeta * std::sqrt(kSpoon * p.spoonMass);   // the old spoon's damper value; its mass is not a mass here (design section 5)
         cZ = 2.0 * p.skirtZeta * std::sqrt(kSpring * m);
         taWall = p.coneDepthMm / p.coneRadiusMm;
         const double wall = std::atan(taWall);
         w = coneRoundMm * JB::MM * (1.0 / std::cos(wall) - 1.0) / taWall;   // puts D's wall line where a round of this radius puts it (0.414 x the round at 45 deg)
         dClose = p.switchCloseMm * JB::MM;
         dOpen = (p.switchCloseMm - p.switchHystMm) * JB::MM;
         // The tilt's damper, zeta of a viscous damper on the tilt mode's own stiffness: the spring's 0.5 k a^2 plus the spoon's kSpoon s^2, s = D's slope per tilt where the tilt alone brings D to the closing distance
         // (the cone's rounded tip makes the slope 0 at its apex, so the closing point; in the wall the slope stays within 1 % of it, runs/2026-11-02-bumper-kick/notes.md)
         {
            const auto DOfTilt = [&](const double th) {
               const double dq[4] = { 0.0, 0.0, th, 0.0 };
               const Pose q = PoseAt(dq);
               return D(q) - zRest + q.z;
            };
            double lo = 1e-6, hi = 0.3;
            for (int i = 0; i < 80; i++)
            {
               const double mid = 0.5 * (lo + hi);
               (DOfTilt(mid) < dClose ? lo : hi) = mid;
            }
            const double thClose = 0.5 * (lo + hi), e = 1e-7;
            const double s = (DOfTilt(thClose + e) - DOfTilt(thClose - e)) / (2.0 * e);
            const double kTilt = kSpring * aSpring * aSpring * 0.5 + kSpoon * s * s;
            cTilt = 2.0 * p.skirtZeta * std::sqrt(kTilt * I);
         }
         invM[0] = 1.0 / m;
         invM[1] = invM[2] = 1.0 / I;
         omLip = std::max(std::sqrt(kSeat / m), std::sqrt(kSeat * rLip * rLip * 0.5 / I));
         omRim = std::max(std::sqrt(kSeat / m), std::sqrt(kSeat * rRim * rRim * 0.5 / I));
         reach = std::max(rRim, Lt);
         mR = p.ringMass;
         invM[3] = 1.0 / mR;
         stroke = p.ringStrokeMm * JB::MM;
         bobbin = p.coil.bobbinMm * JB::MM;
         fReturn = p.returnPreloadN;
         kReturn = p.returnKNPerM;
         cRing = p.ringViscousNsPerM;
         Reset();
      }
      // At rest on its seat, on the shared g: the lip pressed by the springs less the weight, no tilt. Nothing to settle
      void Reset()
      {
         const double a = (fSpring + fSpoon + m * set->g.z) / (kSeat + kSpring + kSpoon);
         z = zRest + a;
         thx = thy = vz = wx = wy = t = 0.0;
         closed = false;
         closures = 0;
         tClose = -1.0;
         bk = PartBooks();
         dr = vr = iCoil = pulseLeft = 0.0;
         ringStop = -1;
         pulseActive = false;
         pulseArmed = true;
         extOn = false;
         extVolts = -1.0;
         extDuty = 1.0;
         coilOnS = coilOnMaxS = 0.0;
         swLog.clear();
         gTriedFree = false;
      }
      // One pulse now, as the table's or a test's solenoid command
      void FirePulse()
      {
         pulseActive = true;
         pulseLeft = cfg.p.pulseMs * 1e-3;
         pulseArmed = false;
         if (bk.tPulseStart < 0.0)
            bk.tPulseStart = t;
         bk.pulses++;
      }
      // The coil's command now, by drive mode: held follows the switch, pulse its timer, external the script. It changes between steps only
      bool CoilOn() const { return cfg.p.mode == JB::DriveMode::Held ? coilEnabled && closed : cfg.p.mode == JB::DriveMode::External ? extOn : pulseActive; }
      double CoilVolts() const { return !CoilOn() ? 0.0 : cfg.p.mode == JB::DriveMode::External ? (extVolts > 0.0 ? extVolts : cfg.p.coil.vSupply) * extDuty : cfg.p.coil.vSupply; }
      double RingDrop(const double *const trial) const { return dr + (trial != nullptr ? trial[kR] : 0.0); }
      double FieldJ() const { return 0.5 * InductanceH(Frac(dr)) * iCoil * iCoil; }
      double FollowerD() const { return D(PoseAt(nullptr)); }
      double TiltDeg() const { return std::sqrt(thx * thx + thy * thy) / JB::DEG; }

      // The skirt's top for a ball centre and radius, the coordinates moved by trial (null = where they are): depth, normal, radius, and how its point moves with each coordinate
      void ProbeSkirt(const V3 &c, const double R, const double *const trial, Hit &h) const
      {
         const Pose q = PoseAt(trial);
         const V3 b = q.ToBody(c);
         const double rb = std::sqrt(b.x * b.x + b.y * b.y);
         const double cb = rb > 1e-9 ? b.x / rb : 1.0, sb = rb > 1e-9 ? b.y / rb : 0.0;
         const JB::Revolve::Near nr = cfg.skirt0.Closest(rb, b.z, cfg.p.edgeMinMm * JB::MM);
         h.delta = R - nr.g;
         // the depth's exact gradient, so the step's correction stays second order in the move (the sweep's radius tapers: the normal alone is out by that taper, and the correction would load the light tilt with it)
         h.grad = std::sqrt(nr.gr * nr.gr + nr.gz * nr.gz);
         h.n = q.Rotate(V3(nr.gr * cb, nr.gr * sb, nr.gz) * (1.0 / h.grad));
         h.reff = KappaEff(R, nr);
         h.curvDist = rb;
         // how a body point moves with each coordinate: the surface's point for its speed, the ball's centre for the depth's slope (the body turning moves the centre in its frame; the two agree where the gradient is the normal)
         const auto cols = [&q](const V3 &bq, V3 &vx, V3 &vy) {
            const V3 a(q.cy * bq.x + q.sy * bq.z, bq.y, -q.sy * bq.x + q.cy * bq.z);
            const V3 da(-q.sy * bq.x + q.cy * bq.z, 0.0, -q.cy * bq.x - q.sy * bq.z);
            vx = V3(0.0, -q.sx * a.y - q.cx * a.z, q.cx * a.y - q.sx * a.z);
            vy = V3(da.x, -q.sx * da.z, q.cx * da.z);
         };
         V3 vx, vy, cx, cy;
         cols(b - V3(nr.nr * cb, nr.nr * sb, nr.nz) * nr.g, vx, vy);
         cols(b, cx, cy);
         const V3 up(0.0, 0.0, 1.0);
         h.nc = 3;
         h.c[0].Set(kZ, h.n.Dot(up), up, V3());
         h.c[1].Set(kThx, h.n.Dot(cx), vx, V3(1.0, 0.0, 0.0));
         h.c[2].Set(kThy, h.n.Dot(cy), vy, V3(0.0, q.cx, q.sx));
      }

      // ---- the group's interface ----
      int NumDof() const override { return ringStop == 0 ? 4 : 3; }
      int MaxDof() const override { return 4; }
      int NumCoords() const override { return 4; }
      int DofOf(const int c) const override { return c != kR || ringStop == 0 ? c : -1; }
      GroupStep::CoordInfo Coord(const int c) const override
      {
         GroupStep::CoordInfo q;
         if (c == kR)
         {
            q.rate = vr;
            q.invM = invM[kR];
            q.invMFree = q.invM;
            q.reach = 1.0;
            q.maxMove = 2e-3;
            q.blockedSide = ringStop;   // the rest stop blocks a smaller drop, the coil stop a larger one
            return q;
         }
         q.rate = c == kZ ? vz : c == kThx ? wx : wy;
         q.invM = invM[c];
         q.invMFree = q.invM;
         q.reach = c == kZ ? 1.0 : reach;
         q.maxMove = c == kZ ? 2e-3 : 0.05;
         return q;
      }
      double KineticJ() const override { return 0.5 * m * vz * vz + 0.5 * I * (wx * wx + wy * wy) + 0.5 * mR * vr * vr; }
      double Energy(const Ball *, int) const override
      {
         const Pose q = PoseAt(nullptr);
         const double dq = D(q), s = zRest - z;
         return KineticJ() + VSmooth(q) + fSpring * s + 0.5 * kSpring * s * s + Seat(ALip(q), rLip * std::sqrt(q.A2)).V + Seat(-z, rRim * std::sqrt(q.A2)).V + fSpoon * dq + 0.5 * kSpoon * dq * dq + mR * set->g.z * dr + fReturn * dr + 0.5 * kReturn * dr * dr;
      }
      // The step's facts, from the state now (accuracy only): the group holds the factors
      void StepFacts(const Ball *, int, const GroupStep::Settings &gs, GroupStep::StepFactList &f) override
      {
         const Pose q = PoseAt(nullptr);
         const double tw = std::sqrt(wx * wx + wy * wy), b = std::sqrt(q.A2);
         const auto seat = [&](const double a, const double r, const double om, const int whyIn, const int whyGap) {
            const double rate = std::fabs(vz) + r * tw;
            if (a > -r * b)
               f.Stiffness(om, rate <= gs.vSteady, whyIn);
            else if (rate > 1e-9)
               f.Approach(-a - r * b, rate, om, whyGap);
         };
         seat(ALip(q), rLip, omLip, 1, 2);
         seat(-z, rRim, omRim, 3, 4);
         if (tw > 1e-6)
            f.Turn(kTurnRadius, tw, 5);
         // the coil: the pulse's end, the current's own time, the steep field on a moving ring
         if (pulseActive)
            f.TimeLeft(pulseLeft, 6);
         const JB::Params &p = cfg.p;
         const double fr = Frac(dr), dLdx = cfg.coilD.b / (p.coil.bobbinMm * (fr + p.coil.x0) * (fr + p.coil.x0));
         const double target = CoilVolts() / p.coil.rOhm;
         if (std::fabs(iCoil - target) > 1e-3 * cfg.coilD.iss)
            f.TimeConstant(InductanceH(fr) / (p.coil.rOhm + dLdx * std::fabs(vr)), 7);
         if (ringStop == 0 && iCoil > 0.0)
            f.Field((fr + p.coil.x0) * bobbin, std::fabs(vr), CoilModel::Force(p.coil, cfg.coilD, fr, iCoil) / mR, 8);
      }
      bool BeginStep(const Ball *, int, const double h, double *x, const int off) override
      {
         x[off + kZ] = vz * h;
         x[off + kThx] = wx * h;
         x[off + kThy] = wy * h;
         if (ringStop == 0)
            x[off + kR] = vr * h;
         gH = h;
         gCoilOn = CoilOn();
         gVolts = CoilVolts();
         const Pose q = PoseAt(nullptr);
         gD0 = D(q);
         gVs0 = VSmooth(q);
         gLip0 = Seat(ALip(q), rLip * std::sqrt(q.A2)).V;
         gRim0 = Seat(-z, rRim * std::sqrt(q.A2)).V;
         gALip0 = ALip(q);
         gCz0 = q.Cz;
         gA20 = q.A2;
         return true;
      }
      bool Eval(const double *x, const int off, const double h, V3 *, V3 *, double *Res, GroupStep::Jac J) override
      {
         double r[3];
         if (J)   // its own block by differences of its own rows: closed forms, so they are cheap, and the seats are C1 only. (A shared analytic form is a cost matter, rule 9)
         {
            double xp[3] = { x[off], x[off + 1], x[off + 2] }, rp[3], rm[3];
            StepLosses ignore;
            for (int c = 0; c < 3; c++)
            {
               const double e = (c == kZ ? 1e-9 : 1e-8) + 1e-6 * std::fabs(x[off + c]);
               xp[c] = x[off + c] + e;
               Rows(xp, h, rp, ignore);
               xp[c] = x[off + c] - e;
               Rows(xp, h, rm, ignore);
               xp[c] = x[off + c];
               for (int a = 0; a < 3; a++)
                  J[off + a][off + c] += (rp[a] - rm[a]) / (2.0 * e);
            }
         }
         Rows(x + off, h, r, gLoss);
         for (int a = 0; a < 3; a++)
            Res[off + a] = r[a];
         // the ring: weight, return spring and viscous loss over the move; the coil and its force over the same move
         gDrive = gRingLoss = 0.0;
         const double dx = ringStop == 0 ? x[off + kR] : 0.0;
         gCoil = CoilOver(dx, h);
         gFCoil = gCoil.force;
         if (ringStop != 0)
            return true;
         const double fe = fReturn + kReturn * (dr + 0.5 * dx) + mR * set->g.z;
         Res[off + kR] = GroupStep::InertiaRow(dx, vr, mR, h) + h * (fe - gFCoil) + cRing * dx;
         gDrive = gFCoil * dx;
         gRingLoss = cRing * dx * dx / h;
         if (J)
         {
            // (the coil's term by differences of the same closed form: the magnet's pull, a negative stiffness, and the closing speed's back EMF)
            const double e = 1e-9 + 1e-6 * std::fabs(dx);
            const double dF = gCoilOn || iCoil != 0.0 ? (CoilOver(dx + e, h).force - CoilOver(dx - e, h).force) / (2.0 * e) : 0.0;
            J[off + kR][off + kR] += GroupStep::InertiaSlope(mR, h) + 0.5 * h * kReturn + cRing - h * dF;
         }
         return true;
      }
      void Loads(const double *q, const V3 &) override { gLoad = q[kR]; }
      // The ring past a stop over the move. With no move (the group's test of where the step starts): a ring leaving a stop isn't on it, and how far off it counts by its rate, so the search for the next arrival doesn't start at the stop it left
      // stops: bit 0 the ring's rest, bit 1 the end of its stroke
      double StopExcess(const double *x, const int off, const unsigned mask = ~0u) const override
      {
         if (ringStop != 0)
            return -INFINITY;
         const double dx = x[off + kR], x1 = dr + dx;
         const bool nearRest = x1 < 0.5 * stroke;
         gArriving = nearRest ? GroupStep::EndRate(dx, vr, gH) < 0.0 : GroupStep::EndRate(dx, vr, gH) > 0.0;
         double e = -INFINITY;
         if (mask & 1u)
            e = std::max(e, -x1);
         if (mask & 2u)
            e = std::max(e, x1 - stroke);
         return e;
      }
      unsigned CrossingStops(const double *x, const int off) const override
      {
         const bool arriving = gArriving;   // (what the last StopExcess found stays for ArrivingAtStop)
         unsigned m = 0u;
         for (unsigned b = 1u; b <= 2u; b <<= 1)
            if (StopExcess(x, off, b) > 0.0)
               m |= b;
         gArriving = arriving;
         return m;
      }
      bool ArrivingAtStop() const override { return gArriving; }
      // On its stop and the move just solved would take it past: it stays, the small speed it had goes
      bool PinOnStop(const double *x, const int off) override
      {
         if (ringStop != 0)
            return false;
         const bool rest = dr + x[off + kR] < 0.5 * stroke;
         if (!gTriedFree && (rest ? vr > 0.0 : vr < 0.0))
            return false;
         gPendingStopJ += 0.5 * mR * vr * vr;
         vr = 0.0;
         ringStop = rest ? -1 : 1;
         return true;
      }
      // On a stop, and the step's forces (coil, weight, spring, the ball) take it off: it has an unknown, and the step is solved again. Once a step: if that move presses it back, it stays for this step (PinOnStop)
      bool FreedByStep() override
      {
         if (ringStop == 0 || gTriedFree)
            return false;
         const double net = RingNetDown();
         if (ringStop < 0 ? net <= 0.0 : net >= 0.0)
            return false;
         ringStop = 0;
         vr = 0.0;
         gTriedFree = true;
         return true;
      }
      void Accept(const double *x, const int off, const double h, const bool onStop, double &driveJ, double &lossJ, bool &through) override
      {
         vz = GroupStep::EndRate(x[off + kZ], vz, h);
         wx = GroupStep::EndRate(x[off + kThx], wx, h);
         wy = GroupStep::EndRate(x[off + kThy], wy, h);
         z += x[off + kZ];
         thx += x[off + kThx];
         thy += x[off + kThy];
         t += h;
         cb.free = ringStop == 0;
         cb.dx = cb.free ? x[off + kR] : 0.0;
         cb.v0 = vr;
         cb.i0 = iCoil;
         cb.h = h;
         cb.load = gLoad;
         cb.dField = -FieldJ();
         double stopJ = gPendingStopJ;
         gPendingStopJ = 0.0;
         if (ringStop == 0)
         {
            const double dx = x[off + kR];
            double v1 = GroupStep::EndRate(dx, vr, h);
            dr += dx;
            const bool rest = dr < 0.5 * stroke;
            if (onStop && (rest ? dr : stroke - dr) < 1e-7)   // the impact: the rate reverses by the restitution; a rebound that would be back within the stop's own contact time never left it
            {
               const double vimp = std::fabs(v1), e = rest ? cfg.p.ringRestE : cfg.p.ringCoilE, vAfter = e * vimp;
               const double net = RingNetDown(), aIn = (rest ? -net : net) / mR;
               const double kStop = rest ? cfg.kRest : cfg.kCoil;
               const bool stays = aIn > 0.0 && (vAfter <= 1e-12 || 2.0 * vAfter / aIn < 2.94 * std::pow(1.25 * mR / kStop, 0.4) * std::pow(vAfter, -0.2));
               stopJ += 0.5 * mR * (stays ? vimp * vimp : vimp * vimp - vAfter * vAfter);
               if (rest)
               {
                  bk.restHits++;
                  if (bk.tFirstRestHit < 0.0)
                  {
                     bk.tFirstRestHit = t;
                     bk.vFirstRestHit = vimp;
                  }
                  bk.tLastRestHit = t;
               }
               else
               {
                  bk.coilHits++;
                  if (bk.tFirstCoilHit < 0.0)
                  {
                     bk.tFirstCoilHit = t;
                     bk.vFirstCoilHit = vimp;
                  }
                  bk.tLastCoilHit = t;
               }
               v1 = stays ? 0.0 : (rest ? vAfter : -vAfter);
               if (stays)
                  ringStop = rest ? -1 : 1;
            }
            vr = v1;
         }
         gTriedFree = false;
         iCoil = gCoil.i1;
         // the step's coil balance, a diagnostic: what the circuit supplied, the field's change and the work
         {
            const JB::Params &pc = cfg.p;
            cb.dField += FieldJ();
            cb.elec = gVolts * gCoil.meanI * h;
            cb.res = pc.coil.rOhm * gCoil.meanI2 * h;
            cb.work = gDrive;
            cb.i1 = iCoil;
            cb.v1 = vr;
         }
         bk.elecJ += gVolts * gCoil.meanI * h;
         bk.resistJ += cfg.p.coil.rOhm * gCoil.meanI2 * h;
         bk.coilJ += gDrive;
         bk.ringStopJ += stopJ;
         bk.ringViscousJ += gRingLoss;
         bk.lipJ += gLoss.lip;
         bk.floorJ += gLoss.rim;
         bk.spoonJ += gLoss.spoon;
         bk.springJ += gLoss.spring;
         bk.viscousJ += gLoss.viscous;
         driveJ = gDrive;
         lossJ = gLoss.lip + gLoss.rim + gLoss.spoon + gLoss.spring + gLoss.viscous + stopJ + gRingLoss;
         through = false;
         coilOnS = gCoilOn ? coilOnS + h : 0.0;
         coilOnMaxS = std::max(coilOnMaxS, coilOnS);
         // the pulse's clock, then the switch at the step's end; in pulse mode a closure fires one pulse, re-armed when the pulse is over and the switch has opened
         if (gCoilOn && pulseActive)
         {
            pulseLeft -= h;
            if (pulseLeft <= 1e-12)
            {
               pulseActive = false;
               bk.tPulseEnd = t;
            }
         }
         const double d = FollowerD();
         if (!closed && d >= dClose)
         {
            closed = true;
            closures++;
            if (tClose < 0.0)
               tClose = t;
            swLog.push_back(1);
            if (cfg.p.mode == JB::DriveMode::Pulse && coilEnabled && pulseArmed && !pulseActive)
               FirePulse();
         }
         else if (closed && d <= dOpen)
         {
            closed = false;
            swLog.push_back(0);
         }
         if (!pulseActive && !closed)
            pulseArmed = true;
      }
      // ---- the net (Group::EnergyStep): counted and reported, not expected to fire ----
      struct Saved { double z, thx, thy, vz, wx, wy, t, tClose, dr, vr, iCoil, pulseLeft, extVolts, extDuty, coilOnS, coilOnMaxS; int closures, ringStop; bool closed, pulseActive, pulseArmed, extOn; size_t nLog; PartBooks bk; };
      Saved sv{};
      void Fill(Saved &o) const
      {
         o.z = z; o.thx = thx; o.thy = thy; o.vz = vz; o.wx = wx; o.wy = wy; o.t = t; o.tClose = tClose; o.dr = dr; o.vr = vr; o.iCoil = iCoil; o.pulseLeft = pulseLeft; o.extVolts = extVolts; o.extDuty = extDuty; o.coilOnS = coilOnS; o.coilOnMaxS = coilOnMaxS;
         o.closures = closures; o.ringStop = ringStop; o.closed = closed; o.pulseActive = pulseActive; o.pulseArmed = pulseArmed; o.extOn = extOn; o.nLog = swLog.size(); o.bk = bk;
      }
      void Load(const Saved &o)
      {
         z = o.z; thx = o.thx; thy = o.thy; vz = o.vz; wx = o.wx; wy = o.wy; t = o.t; tClose = o.tClose; dr = o.dr; vr = o.vr; iCoil = o.iCoil; pulseLeft = o.pulseLeft; extVolts = o.extVolts; extDuty = o.extDuty; coilOnS = o.coilOnS; coilOnMaxS = o.coilOnMaxS;
         closures = o.closures; ringStop = o.ringStop; closed = o.closed; pulseActive = o.pulseActive; pulseArmed = o.pulseArmed; extOn = o.extOn; bk = o.bk;
         gTriedFree = false;
         gPendingStopJ = 0.0;
      }
      void SaveState() override { Fill(sv); }
      void RestoreState() override
      {
         Load(sv);
         swLog.resize(sv.nLog);
      }
      void ScaleSpeeds(const double f) override { vz *= f; wx *= f; wy *= f; vr *= f; }
      void ExplicitMove(const double *Res, const int off, const double h, double *x) const override
      {
         for (int a = 0; a < NumDof(); a++)
            x[off + a] -= Res[off + a] * h * invM[a] * 0.5;
      }

   private:
      // the coil's place: 0 with the ring on its coil stop
      double Frac(const double drop) const { return std::max(0.0, (stroke - drop) / bobbin + cfg.p.plungerClosedFrac); }
      double InductanceH(const double fr) const { return (cfg.coilD.a + cfg.coilD.b / (fr + cfg.p.coil.x0)) * 1e-3; }
      // the circuit and the force over the move dx, from the step's start
      CoilModel::OverMove CoilOver(const double dx, const double h) const
      {
         return CoilModel::Over(cfg.p.coil, cfg.coilD, Frac(dr), Frac(dr + dx), dx, gVolts, iCoil, h);
      }
      // the net force on the ring without its stops, this step's coil and ball, down positive
      double RingNetDown() const { return gFCoil - mR * set->g.z - (fReturn + kReturn * dr) + gLoad; }
      static constexpr double kTurnRadius = 0.125;   // with the group's turn limit of 2 deg: the skirt turns at most 0.25 deg a step (the design's first value)
      struct Pose
      {
         double z, cx, sx, cy, sy, Cz, A2;   // A2: the tilt's sine squared; Cz its cosine
         V3 ToBody(const V3 &p) const
         {
            const V3 r(p.x, p.y, p.z - z), a(r.x, cx * r.y + sx * r.z, -sx * r.y + cx * r.z);
            return V3(cy * a.x - sy * a.z, a.y, sy * a.x + cy * a.z);
         }
         V3 Rotate(const V3 &b) const
         {
            const V3 a(cy * b.x + sy * b.z, b.y, -sy * b.x + cy * b.z);
            return V3(a.x, cx * a.y - sx * a.z, sx * a.y + cx * a.z);
         }
      };
      // The pose with the coordinates moved by dq (null = where they are); s scales the move (0.5: mid-step)
      Pose PoseAt(const double *const dq, const double s = 1.0) const
      {
         Pose q;
         const double tx = thx + (dq != nullptr ? s * dq[kThx] : 0.0), ty = thy + (dq != nullptr ? s * dq[kThy] : 0.0);
         q.z = z + (dq != nullptr ? s * dq[kZ] : 0.0);
         q.cx = std::cos(tx);
         q.sx = std::sin(tx);
         q.cy = std::cos(ty);
         q.sy = std::sin(ty);
         q.Cz = q.cx * q.cy;
         q.A2 = q.cx * q.cx * q.sy * q.sy + q.sx * q.sx;
         return q;
      }
      double ALip(const Pose &q) const { return q.z + q.Cz * tIn - zBase; }
      // The seat: a flat ring on an elastic rim, integrated exactly round the circle (design section 5). a: pressing depth at the centre, b: radius x tilt. Va = dV/da, VbB = (dV/db) / b
      struct SeatE { double V, Va, VbB; };
      SeatE Seat(const double a, const double b) const
      {
         if (a >= b)
            return { 0.5 * kSeat * (a * a + 0.5 * b * b), kSeat * a, 0.5 * kSeat };
         if (a <= -b)
            return { 0.0, 0.0, 0.0 };
         const double r = a / b, p = std::acos(-r), s = std::sqrt(std::max(1.0 - r * r, 0.0)), c = -r;
         return { kSeat / (2.0 * GroupStep::kPi) * b * b * (r * r * p + 2.0 * r * s + 0.5 * p + 0.5 * s * c), kSeat / GroupStep::kPi * b * (r * p + s), kSeat / GroupStep::kPi * (r * s + 0.5 * (p + s * c)) };
      }
      // the gravity and the spring's tilt part; the spring's axial part goes on its own, exact
      double VSmooth(const Pose &q) const
      {
         const V3 g = set->g;
         return 0.25 * kSpring * aSpring * aSpring * q.A2 - m * (g.x * q.sy * cmZ - g.y * q.sx * q.cy * cmZ + g.z * (q.z + q.Cz * cmZ));
      }
      // The follower: the spoon's place if the tip stays on the cone wall (design section 5), with no cancellation at a small tilt
      double D(const Pose &q) const
      {
         const double rho2 = Lt * Lt * q.A2;
         return zRest - q.z - Lt * q.A2 / (1.0 + q.Cz) + taWall * rho2 / (std::sqrt(rho2 + w * w) + w);
      }
      // Gonzalez over the three coordinates: g + (df - g . dq) / (dq^T W dq) W dq, so that the result . dq = df exactly; a mismatch within the rounding of f is left out
      void Discrete(const double *g, const double df, const double noise, const double *dq, double *o) const
      {
         double gdq = 0.0, qq = 0.0;
         for (int a = 0; a < 3; a++)
         {
            gdq += g[a] * dq[a];
            qq += dq[a] * invM[a] * dq[a];
         }
         const double cf = qq > 1e-300 && std::fabs(df - gdq) > noise ? (df - gdq) / qq : 0.0;
         for (int a = 0; a < 3; a++)
            o[a] = g[a] + cf * invM[a] * dq[a];
      }
      static constexpr double kEps16 = 16.0 * 2.220446049250313e-16;
      struct StepLosses { double lip = 0.0, rim = 0.0, spoon = 0.0, spring = 0.0, viscous = 0.0; };
      // The chord of a function of the tilts over their move: g + (df - g . dq) dq / (dq . dq), so that the result . dq = df exactly. The two tilts have one weight, so nothing comes in from another coordinate; a mismatch within the rounding of f is left out
      static void TiltChord(const double g[2], const double df, const double noise, const double *dq, double o[2])
      {
         const double dx = dq[kThx], dy = dq[kThy], qq = dx * dx + dy * dy, mis = df - (g[0] * dx + g[1] * dy);
         const double cf = qq > 1e-300 && std::fabs(mis) > noise ? mis / qq : 0.0;
         o[0] = g[0] + cf * dx;
         o[1] = g[1] + cf * dy;
      }
      // A seat's force over the move, added into Q: the discrete gradient of its energy times the damper's multiplier. Returns the loss, positive or zero.
      // V is a function of two scalars, the depth a and the tilt's A2 (b = r sqrt(A2)), so its gradient is taken in those: dV = Abar da + Cbar dA2 exactly, the mean of the two orders of moving them (reversing the move gives the same force).
      // Each then goes onto the coordinates by its own chord. A 3-coordinate Gonzalez here spreads V's mismatch by the mass weight, and near the seat's onset, where the move is as large as the overlap, the tilt, 3,700 times lighter, took a steep false force from a z move.
      double SeatForce(const double am, const double a0, const double a1, const double dadz, const double tiltA, const double r, const double V0, const Pose &qm, const Pose &q1, const double cz[2], const double ca[2], const double h, double *Q) const
      {
         const double b0 = r * std::sqrt(gA20), b1 = r * std::sqrt(q1.A2);
         const double V1 = Seat(a1, b1).V, V10 = Seat(a1, b0).V, V01 = Seat(a0, b1).V;
         if (V0 <= 0.0 && V1 <= 0.0 && V10 <= 0.0 && V01 <= 0.0)
            return 0.0;
         const double bm = r * std::sqrt(qm.A2);
         const SeatE em = Seat(am, bm);
         const double da = a1 - a0, dA = q1.A2 - gA20, dVa = 0.5 * (V10 - V0 + V1 - V01), dVb = 0.5 * (V01 - V0 + V1 - V10);
         const double noise = kEps16 * (std::max(std::max(V0, V1), std::max(V10, V01)) + std::fabs(em.Va) * zBase + em.VbB * bm * bm);   // (the rounding of V is the rounding of the depth a, whose scale is the height, and of b, times their forces; far above V's own)
         const double Abar = std::fabs(dVa) > noise ? dVa / da : em.Va, Cbar = std::fabs(dVb) > noise ? dVb / dA : em.VbB * 0.5 * r * r;
         const double s0 = std::sqrt(2.0 * V0 / kSeat), s1 = std::sqrt(2.0 * V1 / kSeat), sum = s0 + s1;
         // 1 + fd / fe with fe = k (s0 + s1) / 2, fd = c (s1 - s0) / h: pushes only
         const double mu = sum > 0.0 ? std::max(0.0, 1.0 + 2.0 * cSeat * (s1 - s0) / (h * kSeat * sum)) : 1.0;
         Q[kZ] -= mu * Abar * dadz;
         Q[kThx] -= mu * (Abar * tiltA * cz[0] + Cbar * ca[0]);
         Q[kThy] -= mu * (Abar * tiltA * cz[1] + Cbar * ca[1]);
         return (mu - 1.0) * (V1 - V0);
      }
      // The part's own rows over the move dq: the shared inertia row (2 M dq / h - 2 M v) - h Q, and the losses
      void Rows(const double *dq, const double h, double *r, StepLosses &l) const
      {
         const Pose qm = PoseAt(dq, 0.5), q1 = PoseAt(dq);
         double Q[3] = { 0.0, 0.0, 0.0 };
         const double dCz[2] = { -qm.sx * qm.cy, -qm.cx * qm.sy };
         const double dA2[2] = { 2.0 * qm.sx * qm.cx * qm.cy * qm.cy, 2.0 * qm.cx * qm.cx * qm.sy * qm.cy };
         // gravity and the spring's tilt
         {
            const V3 g = set->g;
            const double k4 = 0.25 * kSpring * aSpring * aSpring, gm = m * cmZ;
            const double gr[3] = { -m * g.z, k4 * dA2[0] - gm * (-g.y * qm.cx * qm.cy + g.z * dCz[0]), k4 * dA2[1] - gm * (g.x * qm.cy + g.y * qm.sx * qm.sy + g.z * dCz[1]) };
            const double V1 = VSmooth(q1);
            double gd[3];
            Discrete(gr, V1 - gVs0, kEps16 * std::max(std::fabs(V1), std::fabs(gVs0)), dq, gd);
            for (int a = 0; a < 3; a++)
               Q[a] -= gd[a];
         }
         // the spring's axial part: exact, along z only; it pushes, with its damper
         {
            const double ds = -dq[kZ], fe = fSpring + kSpring * (zRest - z + 0.5 * ds), N = std::max(0.0, fe + cZ * ds / h);
            Q[kZ] += N;
            l.spring = (N - fe) * ds;
         }
         // the seats, on the chords of the tilt's two functions
         {
            double cz[2], ca[2];
            TiltChord(dCz, q1.Cz - gCz0, kEps16, dq, cz);
            TiltChord(dA2, q1.A2 - gA20, kEps16 * std::max(q1.A2, gA20), dq, ca);
            l.lip = SeatForce(ALip(qm), gALip0, ALip(q1), 1.0, tIn, rLip, gLip0, qm, q1, cz, ca, h, Q);
            l.rim = SeatForce(-qm.z, -z, -q1.z, -1.0, 0.0, rRim, gRim0, qm, q1, cz, ca, h, Q);
         }
         // the spoon: a spring on D with its damper, along D's discrete gradient
         {
            const double D1 = D(q1), dD = D1 - gD0, rho2 = Lt * Lt * qm.A2, root = std::sqrt(rho2 + w * w), wallK = taWall * Lt * Lt * 0.5 / root;
            const double gr[3] = { -1.0, Lt * dCz[0] + wallK * dA2[0], Lt * dCz[1] + wallK * dA2[1] };
            double gd[3];
            Discrete(gr, dD, kEps16 * Lt, dq, gd);
            const double fe = fSpoon + kSpoon * (gD0 + 0.5 * dD), N = std::max(0.0, fe + cSpoon * dD / h);
            for (int a = 0; a < 3; a++)
               Q[a] -= N * gd[a];
            l.spoon = (N - fe) * dD;
         }
         // the skirt's viscous loss on its tilts' mid-step rates (the axial one is with the spring)
         Q[kThx] -= cTilt * dq[kThx] / h;
         Q[kThy] -= cTilt * dq[kThy] / h;
         l.viscous = cTilt * (dq[kThx] * dq[kThx] + dq[kThy] * dq[kThy]) / h;
         const double v[3] = { vz, wx, wy }, M[3] = { m, I, I };
         for (int a = 0; a < 3; a++)
            r[a] = GroupStep::InertiaRow(dq[a], v[a], M[a], h) - h * Q[a];
      }

      const GroupStep::Settings *set = nullptr;
      double invM[4] = { 0.0, 0.0, 0.0, 0.0 };
      double mR = 0.0, stroke = 0.0, bobbin = 0.0, fReturn = 0.0, kReturn = 0.0, cRing = 0.0;
      CoilModel::OverMove gCoil;
      bool gCoilOn = false, gTriedFree = false;
      mutable bool gArriving = true;
      double gH = 1e-6, gVolts = 0.0, gFCoil = 0.0, gLoad = 0.0, gDrive = 0.0, gRingLoss = 0.0, gPendingStopJ = 0.0;
      double m = 0.0, I = 0.0, cmZ = 0.0, tIn = 0.0, zRest = 0.0, zBase = 0.0, rLip = 0.0, rRim = 0.0, Lt = 0.0, kSpring = 0.0, fSpring = 0.0, aSpring = 0.0, kSeat = 0.0, cSeat = 0.0, kSpoon = 0.0, fSpoon = 0.0, cSpoon = 0.0, cZ = 0.0, cTilt = 0.0, taWall = 1.0, w = 0.0, dClose = 0.0, dOpen = 0.0, omLip = 0.0, omRim = 0.0;
      double gD0 = 0.0, gVs0 = 0.0, gLip0 = 0.0, gRim0 = 0.0, gALip0 = 0.0, gCz0 = 1.0, gA20 = 0.0;
      StepLosses gLoss;
   };

   // The skirt's top: a surface the part carries
   class SkirtShape : public GroupStep::Shape
   {
   public:
      const Part *owner = nullptr;
      void Probe(const V3 &c, const double R, Hit &h) const override { owner->ProbeSkirt(c, R, trial, h); }
   };
   // A surface of revolution round the bumper's axis: the centre body, fixed; the ring, carried by the part along its drop (one column, down)
   class RevolveShape : public GroupStep::Shape
   {
   public:
      const JB::Revolve *prof = nullptr;
      double edgeMin = 0.0;
      const Part *owner = nullptr;   // null = fixed
      void Probe(const V3 &c, const double R, Hit &h) const override
      {
         const double rho = std::sqrt(c.x * c.x + c.y * c.y);
         const double cphi = rho > 1e-9 ? c.x / rho : 1.0, sphi = rho > 1e-9 ? c.y / rho : 0.0;
         const JB::Revolve::Near nr = prof->Closest(rho, c.z + (owner != nullptr ? owner->RingDrop(trial) : 0.0), edgeMin);
         h.delta = R - nr.g;
         h.n = V3(nr.nr * cphi, nr.nr * sphi, nr.nz);
         h.reff = KappaEff(R, nr);
         h.curvDist = rho;
         h.grad = 1.0;
         if (owner != nullptr)
         {
            const V3 down(0.0, 0.0, -1.0);
            h.nc = 1;
            h.c[0].Set(Part::kR, h.n.Dot(down), down, V3());
         }
      }
   };
}
