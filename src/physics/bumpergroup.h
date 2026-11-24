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

   struct StepLosses { double lip = 0.0, rim = 0.0, spoon = 0.0, spring = 0.0, viscous = 0.0; };

   // The part's values by lifetime (PartRecords): its state, what one step keeps over its attempts, what one attempt keeps
   struct BumperState
   {
      PartBooks bk;
      double z = 0.0, thx = 0.0, thy = 0.0, vz = 0.0, wx = 0.0, wy = 0.0;
      double t = 0.0;           // s, the time the part has been stepped
      bool closed = false;      // the switch
      int closures = 0;
      double tClose = -1.0;     // the first closure
      double dr = 0.0, vr = 0.0;   // the ring's drop and rate, down positive
      int ringStop = -1;        // -1 on the rest stop (no unknown), +1 on the coil stop (no unknown), 0 free
      double iCoil = 0.0;       // A
      bool pulseActive = false, pulseArmed = true;
      double pulseLeft = 0.0;
      bool pulsePending = false;   // a jet's pulse at a cycle's start: a closure armed a pulse that starts at a cycle's start (CycleStart), pendingS after the next cycle's start
      double pendingS = 0.0;
      bool extOn = false;       // external drive: the script's command
      double extVolts = -1.0, extDuty = 1.0;
      double coilOnS = 0.0, coilOnMaxS = 0.0;   // how long the coil has been commanded on, this episode and the longest (a detector's)
      double gPendingStopJ = 0.0;   // the ring's energy a stop took in PinOnStop, booked by the next Accept (a retake's Accept books it too)
   };
   struct BumperStepWork
   {
      bool gTriedFree = false;   // FreedByStep freed the ring in this step: once a step
      bool gSeatOut[2] = { false, false };   // touch_predict design 3.8f: the seat (0 lip, 1 rim) left out of this step's solve (SpringsOut, from the start state; SpringIn puts it back)
   };
   struct BumperAttemptWork
   {
      CoilModel::OverMove gCoil;
      bool gCoilOn = false;
      mutable bool gArriving = true;
      double gLoad = 0.0;   // the balls' load on the ring at this attempt's evaluation (Loads); a step alone has no ball, so none (as a hinge's ClearLoads)
      double gH = 1e-6, gVolts = 0.0, gFCoil = 0.0, gDrive = 0.0, gRingLoss = 0.0;
      double gD0 = 0.0, gLip0 = 0.0, gRim0 = 0.0, gALip0 = 0.0, gCz0 = 1.0, gA20 = 0.0;
      double gCx0 = 1.0, gSx0 = 0.0, gCy0 = 1.0;   // the tilts' cosines and sine at the step's start (MoveOf)
      double gCtcQ[3] = { 0.0, 0.0, 0.0 };         // diagnostics (VPX_SKIRT_ROW_TRACE): the balls' generalised forces on z, thx, thy at this evaluation (Loads)
      StepLosses gLoss;
   };

   class Part : public GroupStep::Part, public BumperState, protected BumperStepWork, protected BumperAttemptWork
   {
   public:
      Part()
      {
         KeepState<BumperState>(this);
         WipeEachStep<BumperStepWork>(this);
         WipeEachAttempt<BumperAttemptWork>(this);
         KeepQueue(&swLog);
      }
      static constexpr int kZ = 0, kThx = 1, kThy = 2, kR = 3;   // the coordinates
      JB::Bumper cfg;           // its geometry, masses and materials
      double reach = 0.05;      // m, the lever that turns a tilt rate into a speed
      double TiltInertia() const { return I; }   // the skirt's tilt inertia about its pivot, kg m^2 (the energy reach)
      bool coilEnabled = true;  // the table's (between cycles)
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
         pulsePending = false;
         pendingS = 0.0;
         extOn = false;
         extVolts = -1.0;
         extDuty = 1.0;
         coilOnS = coilOnMaxS = 0.0;
         swLog.clear();
         gTriedFree = false;
         gPendingStopJ = 0.0;
         restJ = Energy(nullptr, 0);   // (one-pass bound: RestJ)
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
      // a jet's pulse at a cycle's start (the sling's form, hitsling.h): between cycles, a pending pulse starts at this cycle's start once the bumper's delay (jet_coil_delay_ms) has
      // passed; cycleS: the cycle about to start. So the coil's state for a cycle is known at its start (ContactBumperPart::EnergyMove)
      void CycleStart(const double cycleS)
      {
         if (!pulsePending)
            return;
         if (pendingS > 1e-12)
         {
            pendingS -= cycleS;
            return;
         }
         pulsePending = false;
         FirePulse();
      }
      // The coil's command now, by drive mode: held follows the switch, pulse its timer, external the script. It changes between steps only
      bool CoilOn() const { return cfg.p.mode == JB::DriveMode::Held ? coilEnabled && closed : cfg.p.mode == JB::DriveMode::External ? extOn : pulseActive; }
      double CoilVolts() const { return !CoilOn() ? 0.0 : cfg.p.mode == JB::DriveMode::External ? (extVolts > 0.0 ? extVolts : cfg.p.coil.vSupply) * extDuty : cfg.p.coil.vSupply; }
      double RingDrop(const double *const trial) const { return dr + (trial != nullptr ? trial[kR] : 0.0); }
      double FieldJ() const { return 0.5 * InductanceH(Frac(dr)) * iCoil * iCoil; }
      double FollowerD() const { return D(PoseAt(nullptr)); }
      double TiltDeg() const { return std::sqrt(thx * thx + thy * thy) / JB::DEG; }

      // the part finder, the path test's motion term of the skirt (group A; GroupStep::Shape::PathMotion), P in the bumper's frame: it turns about its pivot (0, 0, z) by at most
      // sqrt(dthx^2 + dthy^2) (Ry Rx: two poses differ by one turn of at most that, paper_groupA.txt) and the pivot slides |dz|: (|P - pivot at the middle| + dLen / 2) |dtilt| + |dz|.
      // Its probe (ProbeSkirt: JB::Revolve::Closest in its own frame) is the profile's distance at any range
      double SkirtMotion(const V3 &P, const double dLen, const double *const trial) const
      {
         if (trial == nullptr)
            return 0.0;
         const double zm = z + 0.5 * trial[kZ];
         return (std::sqrt(P.x * P.x + P.y * P.y + (P.z - zm) * (P.z - zm)) + 0.5 * dLen) * std::hypot(trial[kThx], trial[kThy]) + std::fabs(trial[kZ]);
      }
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
      // the one-pass bound (ContactMovingPart::EnergyCap, ContactBumperPart, design.md AGREED item 4b): RestJ: the book's Energy at the rest Reset places it in (on its seat, springs
      // against the weight, no tilt, the ring up, still), the least of its stored energy (INFERRED: the rest is that least; the table's slope tilts the true least by micro-joules);
      // SeatK: the skirt seat's linear stiffness (N/m) for the squash term past it
      double RestJ() const { return restJ; }
      double SeatK() const { return kSeat; }
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
         const auto seat = [&](const bool lip, const double a, const double r, const double om, const int whyIn, const int whyGap) {
            const double rate = std::fabs(vz) + r * tw;
            if (a > -r * b)
               f.Stiffness(om, rate <= gs.vSteady, whyIn);
            else if (f.predictOn)   // (touch_predict design 3.8f: within d at the start engaged; else left out of the step's solve, SpringsOut: no limit)
            {
               double vc, gap;
               SeatNow(q, lip, gap, vc);
               f.Spring(gap, vc, om, whyGap);
            }
            else if (rate > 1e-9)
               f.Approach(-a - r * b, rate, om, whyGap);   // (touch_predict 0: today's bound)
         };
         seat(true, ALip(q), rLip, omLip, 1, 2);
         seat(false, -z, rRim, omRim, 3, 4);
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
         gLip0 = Seat(ALip(q), rLip * std::sqrt(q.A2)).V;
         gRim0 = Seat(-z, rRim * std::sqrt(q.A2)).V;
         gALip0 = ALip(q);
         gCz0 = q.Cz;
         gA20 = q.A2;
         gCx0 = q.cx;
         gSx0 = q.sx;
         gCy0 = q.cy;
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
               RowsTerms(xp, h, rp, ignore, &rtP);   // (the terms traced: copies only, the halved line's slopes)
               xp[c] = x[off + c] - e;
               RowsTerms(xp, h, rm, ignore, &rtM);
               xp[c] = x[off + c];
               for (int a = 0; a < 3; a++)
                  J[off + a][off + c] += (rp[a] - rm[a]) / (2.0 * e);
               eNow[c] = e;
               for (int k = 0; k < 6; k++)
                  for (int a = 0; a < 3; a++)
                     slopeNow[c][k][a] = (TermOf(rtP, k)[a] - TermOf(rtM, k)[a]) / (2.0 * e);
            }
            slopeSet = true;
         }
         Rows(x + off, h, r, gLoss);
         {   // (the halved line's history: this evaluation, copies only)
            HistE &he = hist[histN % kHist];
            he.t = rt;
            he.hasSlope = slopeSet;
            std::memcpy(he.slope, slopeNow, sizeof slopeNow);
            std::memcpy(he.e, eNow, sizeof eNow);
            if (++histN > (1 << 30))
               histN = kHist + histN % kHist;
         }
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
      void Loads(const double *q, const V3 &) override
      {
         gLoad = q[kR];
         for (int a = 0; a < 3; a++)   // (diagnostics: the row terms' ctc)
            gCtcQ[a] = q[a];
      }
      // diagnostics (VPX_SKIRT_ROW_TRACE=1, with VPX_GROUP_ROUND_TRACE): the skirt's three rows term by term at the round's last evaluation, in the round trace's units (m/s: the row x invM x reach),
      // with the move, the step-start tilt, and each function's chord: its change df (MoveOf) and the mismatch df - g . dq
      static bool SkirtRowTraceOn() { static const bool on = []() { const char *const e = getenv("VPX_SKIRT_ROW_TRACE"); return e != nullptr && atoi(e) != 0; }(); return on; }
      std::string RowNote() const override
      {
         if (!SkirtRowTraceOn())
            return std::string();
         return TermsNote(-1);
      }
      // the halved line (GroupStep::Group::LogHalved): each failed attempt's skirt rows term by term at its last round, the format of RowNote
      void KeepFailNote(const int slot, const double *, const int) override
      {
         if (slot >= 0 && slot < GroupStep::Group::kHalveDiag)
         {
            failRt[slot] = rt;
            const int n = std::min(histN, kHist);
            failHistN[slot] = n;
            for (int i = 0; i < n; i++)   // oldest first
               failHist[slot][i] = hist[(histN - n + i) % kHist];
         }
         histN = 0;
      }
      std::string FailNote(const int slot) const override { return slot >= 0 && slot < GroupStep::Group::kHalveDiag ? TermsNote(slot) + HistNote(slot) : std::string(); }
      // the failed attempt's last evaluations, each against the one before: the move per coordinate (m, rad), the branches (lip on/Abar/Cbar/mu0, rim the same, spring, spoon),
      // and per row and term its change d against its slope x the move p (m/s as RowNote; the slope of the latest Jacobian before; ctc: its change only)
      std::string HistNote(const int slot) const
      {
         std::string o;
         char s[200];
         static const char *const tn[6] = { "inert", "grav", "spring", "seat", "spoon", "damp" }, *const rn[3] = { "z", "thx", "thy" };
         const HistE *const H = failHist[slot];
         if (failHistN[slot] > 0)
         {
            snprintf(s, sizeof s, " || jac e z %.3g m thx %.3g thy %.3g rad", H[failHistN[slot] - 1].e[0], H[failHistN[slot] - 1].e[1], H[failHistN[slot] - 1].e[2]);
            o += s;
         }
         for (int i = 1; i < failHistN[slot]; i++)
         {
            const HistE &e0 = H[i - 1], &e1 = H[i];
            double d[3];
            for (int c = 0; c < 3; c++)
               d[c] = e1.t.dq[c] - e0.t.dq[c];
            char st0[11], st1[11];
            for (int k = 0; k < 10; k++)
            {
               st0[k] = (char)('0' + e0.t.st[k]);
               st1[k] = (char)('0' + e1.t.st[k]);
            }
            st0[10] = st1[10] = 0;
            snprintf(s, sizeof s, " || ev%d move z %.3g thx %.3g thy %.3g st %s>%s", i, d[0], d[1], d[2], st0, st1);
            o += s;
            for (int a = 0; a < 3; a++)
            {
               const double sc = invM[a] * (a == kZ ? 1.0 : reach);
               o += std::string(" | ") + rn[a] + ":";
               for (int k = 0; k < 6; k++)
               {
                  const double dt = (TermOf(e1.t, k)[a] - TermOf(e0.t, k)[a]) * sc;
                  double p = 0.0;
                  for (int c = 0; c < 3; c++)
                     p += e0.slope[c][k][a] * d[c];
                  p *= sc;
                  if (dt == 0.0 && p == 0.0)
                     continue;
                  snprintf(s, sizeof s, " %s d %.4g p %.4g", tn[k], dt, e0.hasSlope ? p : 0.0);
                  o += s;
               }
               snprintf(s, sizeof s, " ctc d %.4g own %.4g", (e1.t.ctc[a] - e0.t.ctc[a]) * sc, (e1.t.own[a] + e1.t.ctc[a]) * sc);
               o += s;
            }
            // each seat's force (SeatForceAt, sd): V0 and V1 - V0; the branch taken (ch chord, md midpoint) with the chord's and the midpoint's estimated errors; da, dA;
            // dVa and dVb; Abar and Cbar with the midpoint's; mu and s1 - s0; then each one's change from the evaluation before and its share of the seat term's thx change (m/s as
            // above: h mu dAbar tiltA cz, h mu dCbar ca, h dmu (Abar tiltA cz + Cbar ca)) and the seat's whole thx change
            static const char *const sn[2] = { "lip", "rim" };
            for (int j = 0; j < 2; j++)
            {
               const double *const a0 = e0.t.sd[j], *const a1 = e1.t.sd[j];
               if (e1.t.st[4 * j] == 0 && e0.t.st[4 * j] == 0)
                  continue;
               char w[700];
               const double sc = invM[kThx] * reach;
               const auto X = [](const double *v) { return v[20] * v[13] * (v[9] * v[15] * v[16] + v[11] * v[18]); };
               snprintf(w, sizeof w, " | %s V0 %.17g V1-V0 %.4g %s ech %.3g emd %.3g da %.4g dA %.4g dVa %.6g dVb %.6g A %.9g Am %.9g C %.9g Cm %.9g mu %.12g ds %.4g"
                                     " | ch dV0 %.3g ddVa %.3g ddVb %.3g dAbar %.3g dCbar %.3g dmu %.3g thx A %.3g C %.3g mu %.3g seat %.3g",
                        sn[j], a1[0], a1[1] - a1[0], a1[8] != 0.0 ? "ch" : "md", a1[2], a1[3], a1[4], a1[5], a1[6], a1[7],
                        a1[9], a1[10], a1[11], a1[12], a1[13], a1[14],
                        a1[0] - a0[0], a1[6] - a0[6], a1[7] - a0[7], a1[9] - a0[9], a1[11] - a0[11], a1[13] - a0[13],
                        a1[20] * a1[13] * (a1[9] - a0[9]) * a1[15] * a1[16] * sc, a1[20] * a1[13] * (a1[11] - a0[11]) * a1[18] * sc,
                        a1[20] * (a1[13] - a0[13]) * (a0[9] * a0[15] * a0[16] + a0[11] * a0[18]) * sc, (X(a1) - X(a0)) * sc);
               o += w;
            }
         }
         return o;
      }
      std::string TermsNote(const int slot) const   // (slot -1: the last evaluation's terms)
      {
         const auto &t = slot < 0 ? rt : failRt[slot];
         char line[2400];
         int n = snprintf(line, sizeof line, "skirt dq z %.6g thx %.6g thy %.6g th0 %.9g %.9g |", t.dq[0], t.dq[1], t.dq[2], thx, thy);
         static const char *const gn[4] = { "Vs", "Cz", "A2", "D" };
         for (int k = 0; k < 4; k++)
            n += snprintf(line + n, sizeof line - n, " %s df %.6g mis %.4g;", gn[k], t.df[k], t.mis[k]);
         static const char *const rn[3] = { "z", "thx", "thy" };
         for (int a = 0; a < 3; a++)
         {
            const double s = invM[a] * (a == kZ ? 1.0 : reach);
            n += snprintf(line + n, sizeof line - n, " | %s: inert %.15g grav %.15g spring %.15g seat %.15g spoon %.15g damp %.15g ctc %.15g own %.6g own+ctc %.6g", rn[a], t.inert[a] * s, t.grav[a] * s,
                          t.spring[a] * s, t.seat[a] * s, t.spoon[a] * s, t.damp[a] * s, t.ctc[a] * s, t.own[a] * s, (t.own[a] + t.ctc[a]) * s);
         }
         return line;
      }
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
      // touch_predict piece 3 (design 3.2a): the same stops at the step's start (bit 0 the rest, bit 1 the stroke's end) on dr, the ring's mass, each stop's own Hertz k
      // (contactbumper.h: cfg.kRest the rest's rubber, cfg.kCoil the stroke end's steel; the values its contact time takes) for d. None while on a stop (no unknown)
      int StopsNow(const double *acc, GroupStep::StopNow *out, const int max) const override
      {
         if (ringStop != 0 || max < 2)
            return 0;
         int n = 0;
         const double ar = acc != nullptr ? acc[kR] : 0.0;
         const auto add = [&](const double gap, const double v, const double a, const double k, const int side, const unsigned id) {
            if (!(gap > 0.0))
               return;
            GroupStep::StopNow &o = out[n++];
            o.gap = gap, o.v = v, o.a = a, o.vs = 0.0, o.m = mR, o.k = k, o.coord = kR, o.side = side, o.id = id;
         };
         add(dr, -vr, -ar, cfg.kRest, -1, 1u);
         add(stroke - dr, vr, ar, cfg.kCoil, 1, 2u);
         return n;
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
      // On its stop and the move just solved would take it past: it stays, the small speed it had goes. Stop rule 1 (bounce), runs/2026-11-22-touch-predict/design.md 3.2d:
      // within d at the start (touch_predict: stopNear > 0) it is on the stop whatever its start rate
      bool PinOnStop(const double *x, const int off) override
      {
         if (ringStop != 0)
            return false;
         const bool rest = dr + x[off + kR] < 0.5 * stroke;
         if (!(stopNear > 0.0) && !gTriedFree && (rest ? vr > 0.0 : vr < 0.0))
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
            if (onStop && (rest ? dr : stroke - dr) < (stopNear > 0.0 ? stopNear : 1e-7))   // (touch_predict: the taken stop's d, design 2.4)   // the impact: the rate reverses by the restitution; a rebound that would be back within the stop's own contact time never left it
            {
               const double vimp = std::fabs(v1), e = rest ? cfg.p.ringRestE : cfg.p.ringCoilE, vAfter = e * vimp;
               const double net = RingNetDown(), aIn = (rest ? -net : net) / mR;
               const double kStop = rest ? cfg.kRest : cfg.kCoil;
               // stop rule 1 (bounce), runs/2026-11-22-touch-predict/design.md 3.2d (touch_predict: stopNear > 0, the taken stop's d): a rebound whose rise against the
               // ring's net force toward the stop, vAfter^2 / (2 aIn), is within d stays on it
               const bool stays = aIn > 0.0 && (stopNear > 0.0 ? vAfter * vAfter <= 2.0 * aIn * stopNear
                                                               : vAfter <= 1e-12 || 2.0 * vAfter / aIn < 2.94 * std::pow(1.25 * mR / kStop, 0.4) * std::pow(vAfter, -0.2));
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
            if (cfg.p.mode == JB::DriveMode::Pulse && coilEnabled && pulseArmed && !pulseActive && !pulsePending)
            {
               pulsePending = true;   // armed here, started at a cycle's start (CycleStart)
               pulseArmed = false;
               pendingS = cfg.p.coilDelayMs * 1e-3;
            }
         }
         else if (closed && d <= dOpen)
         {
            closed = false;
            swLog.push_back(0);
         }
         if (!pulseActive && !pulsePending && !closed)
            pulseArmed = true;
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
      // touch_predict (K2-K4, runs/2026-11-23-touch-p3paper/s_paper.py): a seat's gap (m, > 0 open), lip: -ALip - rLip b, rim: z - rRim b, b = sqrt(A2) (StepFacts' own form)
      double SeatGap(const Pose &q, const bool lip) const { return lip ? -ALip(q) - rLip * std::sqrt(q.A2) : q.z - rRim * std::sqrt(q.A2); }
      // its signed closing rate vc = -(grad g . qdot) and closing acceleration ac = -(grad g . qddot + qdot' H qdot) in the coordinates (z, thx, thy); v, a: qdot, qddot.
      // g = -sz z - tL Cz - r b (+ const), sz +1 lip / -1 rim, tL = tIn lip / 0 rim; Cz = cx cy, b = sqrt(1 - Cz^2) (A2 = 1 - Cz^2 identically). With c_k = dCz/dth_k:
      // b_k = -Cz c_k / b, b_kl = -(c_k c_l + Cz c_kl) / b - Cz^2 c_k c_l / b^3 (the tilt cone's curvature r (|w|^2 - (u.w)^2) / b, paper S1)
      void SeatClosing(const Pose &q, const bool lip, const double v[3], const double a[3], double &vc, double &ac) const
      {
         const double sz = lip ? 1.0 : -1.0, tL = lip ? tIn : 0.0, r = lip ? rLip : rRim, b = std::sqrt(q.A2), c = q.Cz;
         if (!(b > 0.0))   // the cone point (no tilt): one-sided, b' = |w| (never the central 0, paper S2), b'' = w.alpha / |w| (|alpha| at w = 0), Cz'' = -|w|^2
         {
            const double tw = std::sqrt(v[1] * v[1] + v[2] * v[2]), ta = std::sqrt(a[1] * a[1] + a[2] * a[2]);
            vc = sz * v[0] + r * tw;
            ac = sz * a[0] - tL * tw * tw + r * (tw > 0.0 ? (v[1] * a[1] + v[2] * a[2]) / tw : ta);
            return;
         }
         const double ck[2] = { -q.sx * q.cy, -q.cx * q.sy };                               // dCz/dthx, dCz/dthy
         const double ckl[2][2] = { { -c, q.sx * q.sy }, { q.sx * q.sy, -c } };              // d2Cz
         double gk[2], vHv = 0.0;
         for (int k = 0; k < 2; k++)
            gk[k] = -tL * ck[k] + r * c * ck[k] / b;                                         // dg/dth_k = -tL c_k - r b_k
         for (int k = 0; k < 2; k++)
            for (int l = 0; l < 2; l++)
            {
               const double bkl = -(ck[k] * ck[l] + c * ckl[k][l]) / b - c * c * ck[k] * ck[l] / (b * b * b);
               vHv += v[1 + k] * (-tL * ckl[k][l] - r * bkl) * v[1 + l];                    // H_kl = -tL c_kl - r b_kl (z enters g linearly: no H row)
            }
         vc = sz * v[0] - gk[0] * v[1] - gk[1] * v[2];                                       // -(grad g . qdot), dg/dz = -sz
         ac = sz * a[0] - gk[0] * a[1] - gk[1] * a[2] - vHv;                                 // -(grad g . qddot + qdot' H qdot)
      }
      // a seat's gap and signed closing rate at the step's start (the state now)
      void SeatNow(const Pose &q, const bool lip, double &gap, double &vc) const
      {
         const double v3[3] = { vz, wx, wy }, a3[3] = {};
         double ac;
         gap = SeatGap(q, lip);
         SeatClosing(q, lip, v3, a3, vc, ac);
      }
      // touch_predict design 3.8f: a seat (why 2 lip, 4 rim: StepFacts' whyGap) open by more than d at the step's start (d = SpringD(gap, rate, 0, omega), as StepFacts'
      // Spring) is left out of this step's solve (gSeatOut: its force and Jacobian zero; step work)
      int SpringsOut(GroupStep::SpringOut *out, const int max) override
      {
         if (NumDof() == 0)
            return 0;
         const Pose q = PoseAt(nullptr);
         int n = 0;
         for (int s = 0; s < 2 && n < max; s++)
         {
            const bool lip = s == 0;
            const double om = lip ? omLip : omRim;
            double gap, vc;
            SeatNow(q, lip, gap, vc);
            if (!(om > 0.0) || !(gap > GroupStep::StepFactList::SpringD(gap, vc, 0.0, om)))
               continue;
            gSeatOut[s] = true;
            out[n].id = lip ? 2 : 4;
            out[n].gap = gap;
            out[n].v = vc;
            out[n].omega = om;
            n++;
         }
         return n;
      }
      void SpringIn(const int id) override
      {
         for (int s = 0; s < 2; s++)
            if (id < 0 || id == (s == 0 ? 2 : 4))
               gSeatOut[s] = false;
      }
      // the seat's exact gap and signed closing rate at the solved end (why 2 lip, 4 rim: StepFacts' whyGap)
      bool SpringGapEnd(const int why, const double *const dq, const double *const rate, double &gap, double &closing) const override
      {
         if (why != 2 && why != 4)
            return false;
         const Pose q = PoseAt(dq);
         const double v3[3] = { rate[kZ], rate[kThx], rate[kThy] }, a3[3] = {};
         double ac;
         gap = SeatGap(q, why == 2);
         SeatClosing(q, why == 2, v3, a3, closing, ac);
         return true;
      }
      // The seat: a flat ring on an elastic rim, integrated exactly round the circle (design section 5). a: pressing depth at the centre, b: radius x tilt. Va = dV/da, VbB = (dV/db) / b
      struct SeatE { double V, Va, VbB; };
      SeatE Seat(const double a, const double b) const { return SeatAt(kSeat, a, a + b, b); }
      // The partial seat in theta = acos(-a / b), the half-angle of the rim's pressed arc (r = a / b = -cos theta, so acos(-r) = theta, sqrt(1 - r^2) = sin theta):
      // V = k b^2 / (2 pi) F, Va = k b / pi G, VbB = k / pi H, with F = theta + theta/2 cos 2theta - 3/4 sin 2theta, G = sin theta - theta cos theta, H = theta/2 - sin(2 theta)/4
      // (the forms in r, rewritten). Near lift-off (theta -> 0) the closed forms cancel down to theta^5, theta^3, theta^3: relative error about eps / theta^2 (at e = 1 + r = 1e-6
      // V had the wrong sign). Below kThSeries the three come from their Taylor series (entire functions; coefficients exact fractions, correctly rounded):
      // F = sum_{n>=2} (-1)^n 4^n (n-1) / (2n+1)! theta^(2n+1), G = sum_{n>=1} (-1)^(n+1) 2n / (2n+1)! theta^(2n+1), H = sum_{n>=1} (-1)^(n+1) 2^(2n-1) / (2n+1)! theta^(2n+1),
      // 13 terms each (truncation below 1e-20 relative at the switch). The switch is where both forms are within 2 eps of the true value (offline test,
      // runs/2026-11-21-halvings/fix_results.md). At the full side (theta -> pi) nothing cancels.
      struct SeatShape { double F, G, H; };
      static constexpr double kThSeries = 1.5;
      static SeatShape ShapeOf(const double th)
      {
         if (th >= kThSeries)
         {
            const double s2 = std::sin(2.0 * th), c2 = std::cos(2.0 * th);
            return { th + 0.5 * th * c2 - 0.75 * s2, std::sin(th) - th * std::cos(th), 0.5 * th - 0.25 * s2 };
         }
         static constexpr double cF[13] = { 0.13333333333333333, -0.025396825396825397, 0.0021164021164021165, -0.00010261343594676929, 3.2888921777810667e-06, -7.517467834928152e-08,
            1.2897616383455163e-09, -1.7239921648728706e-11, 1.84713446236379e-13, -1.622428161935696e-15, 1.1897806520861771e-17, -7.395683929051606e-20, 3.94679355819667e-22 };
         static constexpr double cG[13] = { 0.3333333333333333, -0.03333333333333333, 0.0011904761904761906, -2.2045855379188714e-05, 2.505210838544172e-07, -1.9270852604185937e-09,
            1.0706029224547743e-11, -4.498331606952833e-14, 1.4797143443923793e-16, -3.9145882126782523e-19, 8.509974375387505e-22, -1.5472680682522736e-24, 2.387759364586842e-27 };
         static constexpr double cH[13] = { 0.3333333333333333, -0.06666666666666667, 0.006349206349206349, -0.0003527336860670194, 1.282667949334616e-05, -3.2888921777810664e-07,
            6.264556529106794e-09, -9.212583131039402e-11, 1.0774951030455441e-12, -1.0261858124243278e-14, 8.11214080967848e-17, -5.408093873118987e-19, 3.081534970438169e-21 };
         const double x = th * th;
         double f = cF[12], g = cG[12], hh = cH[12];
         for (int n = 11; n >= 0; n--)
         {
            f = f * x + cF[n];
            g = g * x + cG[n];
            hh = hh * x + cH[n];
         }
         const double t3 = th * x;
         return { f * t3 * x, g * t3, hh * t3 };
      }
      // theta from the depth a, the sum a + b (formed by the caller from its own parts, where it cancels) and b > 0: on the lift-off side by 2 asin(sqrt(e / 2)),
      // e = (a + b) / b = 1 - cos theta, which keeps e's digits (acos(-r) loses them in r = a / b). As sqrt(a + b) / sqrt(2 b): /fp:fast formed (a + b) / b as a / b + 1,
      // which threw those digits away again (offline test: V off by 6e14 eps at e = 1e-15)
      static double ThetaOf(const double a, const double apb, const double b) { return a < 0.0 ? 2.0 * std::asin(std::sqrt(apb) / std::sqrt(2.0 * b)) : std::acos(-a / b); }
      // The seat (k its stiffness): full (a >= b), partial, or off (a + b <= 0, the rim lifted)
      static SeatE SeatAt(const double k, const double a, const double apb, const double b)
      {
         if (a >= b)
            return { 0.5 * k * (a * a + 0.5 * b * b), k * a, 0.5 * k };
         if (apb <= 0.0)
            return { 0.0, 0.0, 0.0 };
         const SeatShape f = ShapeOf(ThetaOf(a, apb, b));
         return { k / (2.0 * GroupStep::kPi) * b * b * f.F, k / GroupStep::kPi * b * f.G, k / GroupStep::kPi * f.H };
      }
      // The seat energy's change over the straight move (a0, b0) -> (a0 + da, b0 + db) as Ia da + Ib db, (Ia, Ib) its gradient in the two lengths at the rim, the depth a and
      // the tilt's lift b. Two ways to form it, one law: each call takes the one whose estimated error is the smaller (no threshold).
      // - Midpoint: (Ia, Ib) = (Va, VbB b) at the move's middle. Its error comes from the force's curvature (order move^3) and is estimated from the values at hand as
      //   |trapezoid - midpoint| / 3, the trapezoid being the mean of the end forces times the move (the midpoint's leading error is -1/2 the trapezoid's). Exact for the full seat,
      //   whose force is linear in a and b, and free of V's rounding: over a small move it gives the true force where V1 - V0 does not (design.md 1.2, the bottom jet).
      // - Chord: Ia da + Ib db = V1 - V0 exactly, the midpoint force corrected along the move (Gonzalez's discrete gradient, a and b with one weight, both lengths).
      //   Its error is V's rounding, eps max(V0, V1) (V from the theta form is a few eps, also near lift-off), plus the end's: a1 and b1 are known to eps of themselves,
      //   which moves V1 by eps (|Va1 a1| + VbB1 b1^2); near lift-off, where Va a / V ~ 1 / theta^2, that term is the larger (offline test: with eps max(V0, V1) alone
      //   the chord was taken at 1e-20 m moves and its force was 45 % off). Over a move it is a force error of that over the move: the chord wins over a large move,
      //   where the midpoint's curvature error is the larger.
      // e1: the seat at the move's end (the caller's a1, b1, whose V1 also feeds the damper). a + b at the middle is formed from a0 + b0 (exact near lift-off) and da + db.
      struct SeatGrad { double Ia, Ib, Iam, Ibm, errChord, errMid; bool chord; };
      static SeatGrad SeatGradOf(const double k, const double a0, const double b0, const double da, const double db, const double V0, const SeatE &e1, const double b1)
      {
         const SeatE e0 = SeatAt(k, a0, a0 + b0, b0), em = SeatAt(k, a0 + 0.5 * da, (a0 + b0) + 0.5 * (da + db), b0 + 0.5 * db);
         SeatGrad g;
         g.Iam = em.Va;
         g.Ibm = em.VbB * (b0 + 0.5 * db);
         const double mid = g.Iam * da + g.Ibm * db, trap = 0.5 * ((e0.Va + e1.Va) * da + (e0.VbB * b0 + e1.VbB * b1) * db), qq = da * da + db * db;
         g.errMid = std::fabs(trap - mid) / 3.0;
         g.errChord = 2.220446049250313e-16 * (std::max(V0, e1.V) + std::fabs(e1.Va * (a0 + da)) + e1.VbB * b1 * b1);
         g.chord = g.errChord < g.errMid && qq > 1e-300;
         const double cf = g.chord ? (e1.V - V0 - mid) / qq : 0.0;
         g.Ia = g.Iam + cf * da;
         g.Ib = g.Ibm + cf * db;
         return g;
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
      // The chord of a function of the tilts over their move: g + (df - g . dq) dq / (dq . dq), so that the result . dq = df exactly. The two tilts have one weight, so nothing comes in from another coordinate; a mismatch within the rounding of f is left out
      static void TiltChord(const double g[2], const double df, const double noise, const double *dq, double o[2])
      {
         const double dx = dq[kThx], dy = dq[kThy], qq = dx * dx + dy * dy, mis = df - (g[0] * dx + g[1] * dy);
         const double cf = qq > 1e-300 && std::fabs(mis) > noise ? mis / qq : 0.0;
         o[0] = g[0] + cf * dx;
         o[1] = g[1] + cf * dy;
      }
      // A seat's force over the move, added into Q: the discrete gradient of its energy times the damper's multiplier. Returns the loss, positive or zero.
      // V is a function of two scalars, the depth a and the tilt's A2 (b = r sqrt(A2)), so its gradient is taken in those: dV = Abar da + Cbar dA2, from the gradient in (a, b)
      // along the straight move (SeatGradOf: the midpoint force or the chord V1 - V0, whichever has the smaller estimated error).
      // Each then goes onto the coordinates by its own chord. A 3-coordinate Gonzalez here spreads V's mismatch by the mass weight, and near the seat's onset, where the move is as large as the overlap, the tilt, 3,700 times lighter, took a steep false force from a z move.
      // a1, A21: the depth and A2 at the move's end; da, dA: their changes over the move (from MoveOf)
      double SeatForceAt(const double a0, const double a1, const double da, const double dadz, const double tiltA, const double r, const double V0, const double A21, const double dA,
                       const double cz[2], const double ca[2], const double h, double *Q, int *const st = nullptr, double *const sd = nullptr) const
      {
         const double sA0 = std::sqrt(gA20), sA1 = std::sqrt(A21), sAs = sA0 + sA1, b0 = r * sA0, b1 = r * sA1;
         const SeatE e1 = Seat(a1, b1);
         const double V1 = e1.V;
         if (st != nullptr)   // (diagnostics, the halved line: the seat's branches, on / chord taken / - / mu clamped at 0)
            st[0] = st[1] = st[2] = st[3] = 0;
         if (sd != nullptr)   // (diagnostics, the halved line: the seat's inputs and outputs, RowTrace::sd; copies only)
            std::fill(sd, sd + kSd, 0.0);
         // both ends off (a + b <= 0, a half-plane): the whole move is off
         if (V0 <= 0.0 && V1 <= 0.0)
            return 0.0;
         // the move in b from dA by the identity b1 - b0 = r (A21 - A20) / (sqrt A21 + sqrt A20); with no tilt at either end b stays 0
         const double db = sAs > 0.0 ? r * dA / sAs : 0.0;
         const SeatGrad g = SeatGradOf(kSeat, a0, b0, da, db, V0, e1, b1);
         // dV = Abar da + Cbar dA: Abar = Ia; Cbar dA = Ib db, so Cbar = Ib r / (sqrt A21 + sqrt A20). With no tilt at either end, dV/dA2 at b = 0: VbB r^2 / 2.
         const double Cb0 = sAs > 0.0 ? 0.0 : 0.5 * r * r * Seat(a0 + 0.5 * da, 0.0).VbB;
         const double Abar = g.Ia, Cbar = sAs > 0.0 ? g.Ib * r / sAs : Cb0, dV = Abar * da + Cbar * dA;
         const double s0 = std::sqrt(2.0 * V0 / kSeat), s1 = std::sqrt(2.0 * V1 / kSeat), sum = s0 + s1;
         // 1 + fd / fe with fe = k (s0 + s1) / 2, fd = c (s1 - s0) / h: pushes only
         const double mu = sum > 0.0 ? std::max(0.0, 1.0 + 2.0 * cSeat * (s1 - s0) / (h * kSeat * sum)) : 1.0;
         if (st != nullptr)
         {
            st[0] = 1;
            st[1] = g.chord ? 1 : 0;
            st[3] = mu == 0.0 ? 1 : 0;
         }
         if (sd != nullptr)   // slots 2, 3: the chord's and the midpoint's estimated errors; 8: 1 if the chord was taken; 10, 12: the midpoint's Abar, Cbar
         {
            const double v[kSd] = { V0, V1, g.errChord, g.errMid, da, dA, Abar * da, Cbar * dA, g.chord ? 1.0 : 0.0, Abar, g.Iam, Cbar, sAs > 0.0 ? g.Ibm * r / sAs : Cb0, mu, s1 - s0,
                                    tiltA, cz[0], cz[1], ca[0], ca[1], h, dadz };
            std::copy(v, v + kSd, sd);
         }
         Q[kZ] -= mu * Abar * dadz;
         Q[kThx] -= mu * (Abar * tiltA * cz[0] + Cbar * ca[0]);
         Q[kThy] -= mu * (Abar * tiltA * cz[1] + Cbar * ca[1]);
         return (mu - 1.0) * dV;
      }
      // The part's own rows over the move dq (RowsTerms), with VPX_SKIRT_ROW_TRACE its terms for the trace
      // (the terms are kept at every main evaluation, not the Jacobian's differences: a failed attempt's last round goes to the halved line, KeepFailNote; copies only)
      void Rows(const double *dq, const double h, double *r, StepLosses &l) const { RowsTerms(dq, h, r, l, &l == &gLoss ? &rt : nullptr); }
      // the step's change of each function of the tilt the rows take a difference of, from the move by exact identities, so each carries the rounding of its own
      // size and not of the near-1 value it is the change of (Cz = cos thx cos thy is 1 - 9e-5 at 0.78 deg: one ulp of it is a tilt move of 8e-15 rad; bumper2/results.md W1):
      // cos a - cos b = -2 sin((a + b) / 2) sin((a - b) / 2), sin a - sin b = 2 cos((a + b) / 2) sin((a - b) / 2) (the mean angle is the mid-step pose's), a1 b1 - a0 b0 =
      // ((a1 - a0)(b1 + b0) + (b1 - b0)(a1 + a0)) / 2, and A2 = 1 - Cz^2 so dA2 = -dCz (Cz1 + Cz0)
      struct TiltMove { double dCz = 0.0, dA2 = 0.0, dSy = 0.0, dSxCy = 0.0; };
      TiltMove MoveOf(const double *dq, const Pose &qm, const Pose &q1) const
      {
         const double shx = std::sin(0.5 * dq[kThx]), shy = std::sin(0.5 * dq[kThy]);
         const double dCx = -2.0 * qm.sx * shx, dSx = 2.0 * qm.cx * shx, dCy = -2.0 * qm.sy * shy;
         TiltMove t;
         t.dCz = 0.5 * (dCx * (q1.cy + gCy0) + dCy * (q1.cx + gCx0));
         t.dA2 = -t.dCz * (q1.Cz + gCz0);
         t.dSy = 2.0 * qm.cy * shy;
         t.dSxCy = 0.5 * (dSx * (q1.cy + gCy0) + dCy * (q1.sx + gSx0));
         return t;
      }
      // diagnostics (VPX_SKIRT_ROW_TRACE): the rows' terms at the last evaluation (row units: inertia row, -h x each force), the move, the four chords (Vs, Cz, A2, D)
      // st: the branches (lip on, chord taken, -, mu at 0; rim the same; spring N > 0; spoon N > 0)
      // sd: per seat (0 lip, 1 rim) SeatForceAt's inputs and outputs (V0, V1, chord's est. error, midpoint's est. error, da, dA, dVa, dVb, chord taken, Abar, midpoint Abar, Cbar, midpoint Cbar, mu, s1 - s0,
      // tiltA, cz, ca, h, dadz), for the halved line's seat chord note (HistNote; to-do 1.42, which of Abar / Cbar / mu carries the seat term's noise)
      static constexpr int kSd = 22;
      struct RowTrace { double inert[3], grav[3], spring[3], seat[3], spoon[3], damp[3], ctc[3], own[3], dq[3], df[4], mis[4]; int st[10]; double sd[2][kSd]; };
      mutable RowTrace rt = {}, rtP = {}, rtM = {};
      // the halved line's history (diagnostics, copies only): the last kHist main evaluations' terms, move and branches, and the terms' slopes (per coordinate,
      // central differences as the Jacobian takes them) of the latest Jacobian before each
      static constexpr int kHist = 4;
      struct HistE { RowTrace t; double slope[3][6][3]; double e[3]; bool hasSlope; };
      mutable double slopeNow[3][6][3] = {}, eNow[3] = {};
      mutable bool slopeSet = false;
      mutable HistE hist[kHist] = {};
      mutable int histN = 0;
      HistE failHist[GroupStep::Group::kHalveDiag][kHist] = {};
      int failHistN[GroupStep::Group::kHalveDiag] = {};
      static const double *TermOf(const RowTrace &t, const int k) { const double *const p[6] = { t.inert, t.grav, t.spring, t.seat, t.spoon, t.damp }; return p[k]; }
      RowTrace failRt[GroupStep::Group::kHalveDiag] = {};   // each failed attempt's rt at its last round (KeepFailNote)
      // The part's own rows over the move dq: the shared inertia row (2 M dq / h - 2 M v) - h Q, and the losses
      void RowsTerms(const double *dq, const double h, double *r, StepLosses &l, RowTrace *const trTo) const
      {
         const bool tr = trTo != nullptr;   // (trace: the terms go to *trTo, the member rt or the Jacobian's differences rtP / rtM)
         RowTrace &rt = tr ? *trTo : this->rt;
         const Pose qm = PoseAt(dq, 0.5), q1 = PoseAt(dq);
         double Q[3] = { 0.0, 0.0, 0.0 };
         const double dCz[2] = { -qm.sx * qm.cy, -qm.cx * qm.sy };
         const double dA2[2] = { 2.0 * qm.sx * qm.cx * qm.cy * qm.cy, 2.0 * qm.cx * qm.cx * qm.sy * qm.cy };
         const TiltMove tm = MoveOf(dq, qm, q1);
         const double A21 = gA20 + tm.dA2;   // A2 at the move's end, as the seats take it
         double Qt[3] = { 0.0, 0.0, 0.0 };   // (trace: Q before a term)
         const auto mark = [&]() { for (int a = 0; a < 3; a++) Qt[a] = Q[a]; };
         const auto took = [&](double *o) { for (int a = 0; a < 3; a++) o[a] = -h * (Q[a] - Qt[a]); };
         const auto chord = [&](const int k, const double *g, const int n0, const double df) {   // (trace: a chord's change df and its mismatch df - g . dq)
            double gdq = 0.0;
            for (int a = n0; a < 3; a++)
               gdq += g[a - n0] * dq[a];
            rt.df[k] = df;
            rt.mis[k] = df - gdq;
         };
         if (tr)
            mark();
         // gravity and the spring's tilt
         {
            const V3 g = set->g;
            const double k4 = 0.25 * kSpring * aSpring * aSpring, gm = m * cmZ;
            const double gr[3] = { -m * g.z, k4 * dA2[0] - gm * (-g.y * qm.cx * qm.cy + g.z * dCz[0]), k4 * dA2[1] - gm * (g.x * qm.cy + g.y * qm.sx * qm.sy + g.z * dCz[1]) };
            double gd[3];
            const double dVs = k4 * tm.dA2 - m * (g.x * cmZ * tm.dSy - g.y * cmZ * tm.dSxCy + g.z * (dq[kZ] + cmZ * tm.dCz));
            Discrete(gr, dVs, 0.0, dq, gd);
            if (tr)
               chord(0, gr, 0, dVs);
            for (int a = 0; a < 3; a++)
               Q[a] -= gd[a];
         }
         if (tr)
         {
            took(rt.grav);
            mark();
         }
         // the spring's axial part: exact, along z only; it pushes, with its damper
         {
            const double ds = -dq[kZ], fe = fSpring + kSpring * (zRest - z + 0.5 * ds), N = std::max(0.0, fe + cZ * ds / h);
            if (tr)
               rt.st[8] = fe + cZ * ds / h > 0.0 ? 1 : 0;
            Q[kZ] += N;
            l.spring = (N - fe) * ds;
         }
         if (tr)
         {
            took(rt.spring);
            mark();
         }
         // the seats, on the chords of the tilt's two functions
         {
            double cz[2], ca[2];
            TiltChord(dCz, tm.dCz, 0.0, dq, cz);
            TiltChord(dA2, tm.dA2, 0.0, dq, ca);
            const double daLip = dq[kZ] + tIn * tm.dCz;
            l.lip = gSeatOut[0] ? 0.0 : SeatForceAt(gALip0, gALip0 + daLip, daLip, 1.0, tIn, rLip, gLip0, A21, tm.dA2, cz, ca, h, Q, tr ? rt.st : nullptr, tr ? rt.sd[0] : nullptr);   // (touch_predict design 3.8f: a left-out seat: no force)
            l.rim = gSeatOut[1] ? 0.0 : SeatForceAt(-z, -z - dq[kZ], -dq[kZ], -1.0, 0.0, rRim, gRim0, A21, tm.dA2, cz, ca, h, Q, tr ? rt.st + 4 : nullptr, tr ? rt.sd[1] : nullptr);
            if (tr)
            {
               chord(1, dCz, 1, tm.dCz);
               chord(2, dA2, 1, tm.dA2);
            }
         }
         if (tr)
         {
            took(rt.seat);
            mark();
         }
         // the spoon: a spring on D with its damper, along D's discrete gradient (D = zRest - z - Lt (1 - Cz) + taWall (sqrt(rho2 + w^2) - w), rho2 = Lt^2 A2, so its change is
         // -dz + Lt dCz + taWall Lt^2 dA2 / (root1 + root0))
         {
            const double rho2 = Lt * Lt * qm.A2, root = std::sqrt(rho2 + w * w), wallK = taWall * Lt * Lt * 0.5 / root;
            const double dD = -dq[kZ] + Lt * tm.dCz + taWall * Lt * Lt * tm.dA2 / (std::sqrt(Lt * Lt * gA20 + w * w) + std::sqrt(Lt * Lt * A21 + w * w));
            const double gr[3] = { -1.0, Lt * dCz[0] + wallK * dA2[0], Lt * dCz[1] + wallK * dA2[1] };
            double gd[3];
            Discrete(gr, dD, 0.0, dq, gd);
            if (tr)
               chord(3, gr, 0, dD);
            const double fe = fSpoon + kSpoon * (gD0 + 0.5 * dD), N = std::max(0.0, fe + cSpoon * dD / h);
            if (tr)
               rt.st[9] = fe + cSpoon * dD / h > 0.0 ? 1 : 0;
            for (int a = 0; a < 3; a++)
               Q[a] -= N * gd[a];
            l.spoon = (N - fe) * dD;
         }
         if (tr)
         {
            took(rt.spoon);
            mark();
         }
         // the skirt's viscous loss on its tilts' mid-step rates (the axial one is with the spring)
         Q[kThx] -= cTilt * dq[kThx] / h;
         Q[kThy] -= cTilt * dq[kThy] / h;
         l.viscous = cTilt * (dq[kThx] * dq[kThx] + dq[kThy] * dq[kThy]) / h;
         const double v[3] = { vz, wx, wy }, M[3] = { m, I, I };
         for (int a = 0; a < 3; a++)
            r[a] = GroupStep::InertiaRow(dq[a], v[a], M[a], h) - h * Q[a];
         if (tr)
         {
            took(rt.damp);
            for (int a = 0; a < 3; a++)
            {
               rt.inert[a] = GroupStep::InertiaRow(dq[a], v[a], M[a], h);
               rt.ctc[a] = -h * gCtcQ[a];
               rt.own[a] = r[a];
               rt.dq[a] = dq[a];
            }
         }
      }

      const GroupStep::Settings *set = nullptr;
      double invM[4] = { 0.0, 0.0, 0.0, 0.0 };
      double mR = 0.0, stroke = 0.0, bobbin = 0.0, fReturn = 0.0, kReturn = 0.0, cRing = 0.0;
      double restJ = 0.0;   // Energy at Reset's rest (one-pass bound: RestJ)
      double m = 0.0, I = 0.0, cmZ = 0.0, tIn = 0.0, zRest = 0.0, zBase = 0.0, rLip = 0.0, rRim = 0.0, Lt = 0.0, kSpring = 0.0, fSpring = 0.0, aSpring = 0.0, kSeat = 0.0, cSeat = 0.0, kSpoon = 0.0, fSpoon = 0.0, cSpoon = 0.0, cZ = 0.0, cTilt = 0.0, taWall = 1.0, w = 0.0, dClose = 0.0, dOpen = 0.0, omLip = 0.0, omRim = 0.0;
   };

   // The skirt's top: a surface the part carries
   class SkirtShape : public GroupStep::Shape
   {
   public:
      const Part *owner = nullptr;
      void Probe(const V3 &c, const double R, Hit &h) const override { owner->ProbeSkirt(c, R, trial, h); }
   };
   // the part finder (group A): the farthest point of a solid of revolution's profile from its frame's origin (each vertex's |(r, z)| plus its sweep radius: along a segment
   // |core| + rr is convex, so its largest is at a vertex), for the overlap check's bounding radius
   inline double FarR(const JB::Revolve &p)
   {
      double f = 0.0;
      for (const JB::Revolve::Vtx &v : p.v)
         f = std::max(f, std::hypot(v.r, v.z) + v.rr);
      return f;
   }
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
