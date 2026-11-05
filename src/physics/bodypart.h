// license:GPLv3+
//
// One body on one coordinate as a part of the group step: its mass, its speed, its hard stops and its pivot's friction. What is shared by a body that a ball can push (a gate today; a spinner, a disc next). The derived part says only
// what is its own: stored energy, viscous loss, a held drive, the friction's limit, the surfaces it carries. Pure SI, no VPX types.
//
// One unknown: the move dq over a step. Its row is m (2 dq / h - 2 rate) = h Q, so with the solve done the energy closes at any step:
//   stored     V(q), exact difference (the derived part)
//   losses     a viscous c on the mid-step rate; Coulomb friction at its limit, against the step's direction; the derived part's own damper on the move
//   drives     a held force whose work is counted, the push of a driven stop
// Hard things are events of the step (Group::Step), not terms of the equation:
//   stops      the step ends on the stop; the impact law of the stop's own Hertz + damping contact gives the rebound, a loss. A body whose rebound would not leave the squash its pressing force holds (1/2 m v^2 <= 2/5 F d(F)) rests: no unknown while there
//   stick      a body the friction stops: the step ends where the speed is zero; while the other forces stay within the friction limit it has no unknown. So friction never pushes: its work has one sign
//   leaving    a body with no unknown that the step's forces take off has one, and the step is solved again
//   driven     a stop that moves: a body resting on it moves with it, CoordInfo::given; the stop pays for the stored energy and the friction, the group books its work on the balls

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <vector>

#include "physics/groupstep.h"
#include "physics/rubbercontact.h"

namespace GroupStep
{
   // Steel on steel (E 200 GPa, nu 0.3, both bodies): E* = E / (2 (1 - nu^2))
   constexpr double kSteelEstar = 200e9 / (2.0 * (1.0 - 0.3 * 0.3));

   // A stop: side +1 blocks q above at, side -1 blocks q below it; touched at a lever from the axis (1 for a body that slides). A driven stop moves at toward target at rate.
   struct BodyStop
   {
      double at = 0.0;
      int side = 1;
      double lever = 1.0;
      double k = 0.0;          // Hertz stiffness at the lever, N / m^1.5
      bool enabled = true;
      bool armed = true;       // engages only from its free side: one placed while the body is past it waits until the body is on its free side
      bool touching = false;
      double target = 0.0, rate = 0.0;   // rate 0 = fixed
      double speed = 0.0;                // its own rate over the step being solved
      // the restitution of the stop's contact by the closing speed at the lever (the viscoelastic law: it falls slowly with speed); log-spaced from kEV0, a factor kEFac apart
      static constexpr int kE = 12;
      static constexpr double kEV0 = 1e-4, kEFac = 3.0;
      double eTab[kE] = {};
      double Restitution(const double v) const
      {
         const double u = v > kEV0 ? std::log(v / kEV0) / std::log(kEFac) : 0.0;
         const int i = std::min((int)u, kE - 2);
         const double f = std::min(u - i, 1.0);
         return eTab[i] + (eTab[i + 1] - eTab[i]) * f;
      }
      double EffMass(const double m) const { return m / (lever * lever); }
      double HeldEnergy(const double F) const { return F > 0.0 && k > 0.0 ? 0.4 * F * std::pow(F / k, 2.0 / 3.0) : 0.0; }   // the squash a pressing force holds: 2/5 F d, with k d^1.5 = F
      double SpeedNow() const { return rate > 0.0 && at != target ? (target > at ? rate : -rate) : 0.0; }
      double AfterStep(const double h) const
      {
         if (speed == 0.0)
            return at;
         const double a = at + speed * h;
         return speed > 0.0 ? std::min(a, target) : std::max(a, target);
      }
      // k and the restitution table: steel on steel at the effective radius reff, restitution e at the reference speed vRef (m/s at the lever), for the body's mass at the lever
      void Setup(const double mass, const double reff, const double e, const double vRef)
      {
         k = (4.0 / 3.0) * kSteelEstar * std::sqrt(reff);
         RubberContact::RubberParams rp;
         rp.nu = 0.0;   // ImpactRestitutionKK uses E / (1 - nu^2) as E*: pass E* itself
         rp.eMpa = kSteelEstar * 1e-6;
         rp.normalLoss = false;
         const double mEff = mass / (lever * lever), gamma = RubberContact::CalibrateKK(rp, reff, mEff, vRef, e);
         double v = kEV0;
         for (int i = 0; i < kE; i++, v *= kEFac)
            eTab[i] = RubberContact::ImpactRestitutionKK(rp, reff, mEff, gamma, v);
      }
   };

   struct BodyStopEvent
   {
      int stop;
      double speed;   // |rate| at the touch
   };

   struct BodyBooks
   {
      double driveJ = 0.0, viscousJ = 0.0, frictionJ = 0.0, stopJ = 0.0;
      long stopHits = 0, rests = 0, sticks = 0, leaves = 0;
   };

   class BodyPart : public Part
   {
   public:
      double q = 0.0, rate = 0.0, mass = 1.0;
      double reach = 0.05;      // m: the lever that turns the rate into a speed, for the solver's units
      double maxMove = 0.05;    // a Newton round moves it at most this
      double qAtStepStart = 0.0;
      std::vector<BodyStop> stops;
      std::vector<BodyStopEvent> stopEvents;   // until the owner reads them
      BodyBooks bk;
      // [Contact] pivot_model, set by the owner when the part is made: 1 = the pin's stick spring (it holds to the static limit, slides at dynRatio times it);
      // 0 = the slip direction, stuck state and turn event (the part keeps 0 until its owner sets it)
      int pivotModel = 1;   // [Contact] pivot_model (default 1, as the ini); the engine sets it per part
      double dynRatio = 0.85;   // [Contact] mu_dynamic_ratio
      virtual void SetPivotModel(const int m) { pivotModel = m; }   // (a part whose losses depend on the law sets them here)

      // ---- what the derived part supplies ----
      virtual double StoredEnergy(double at) const = 0;
      // (V(at + dq) - V(at)) / dq and its slope with dq; lossF: its own dampers over the move (they only remove energy), and their slope
      virtual double StoredChord(double at, double dq, double h, double &slope, double &lossF, double &lossSlope) const = 0;
      virtual double StoredForce(double at) const = 0;              // -dV / dq
      virtual double Viscous() const { return 0.0; }
      virtual double DriveForce() const { return 0.0; }             // a held force, its work counted
      virtual double FrictionLimit(double at, double r, const V3 &load) const = 0;   // load: the balls' force on the part
      virtual double PinStiffness(double, double, const V3 &) const { return 0.0; }   // the pin's stick spring in its hole, Nm/rad (pivot_model 1)
      virtual void AfterMove() {}                                   // a step was taken: switches

      // Between cycles, after the owner moved a driven stop or changed the gravity: who has an unknown now
      void CycleBegin()
      {
         for (BodyStop &s : stops)
            s.speed = s.SpeedNow();
         Decide();
         DrivenTouch(pendDrive, pendLoss);
         gGiven = 0.0;
      }
      int PlaceStop(BodyStop s)
      {
         s.touching = false;
         s.armed = s.side * (q - s.at) <= 0.0;
         stops.push_back(s);
         return (int)stops.size() - 1;
      }
      // the owner moves a driven stop toward a new place. A body resting on it no longer does: if the stop comes back at it, that is an impact (CycleBegin)
      void DriveStop(const int i, const double target, const double r)
      {
         stops[i].target = target;
         stops[i].rate = std::abs(r);
         if (pinStop == i)
            pinStop = -1;
      }
      void ClearLoads()   // no ball now: what the last one pressed with is gone
      {
         gBallQ = 0.0;
         gBallF = V3();
      }
      bool AnyStopMoving() const
      {
         for (const BodyStop &s : stops)
            if (s.rate > 0.0 && s.at != s.target)
               return true;
         return false;
      }
      bool Pinned() const { return pinStop >= 0; }
      bool Stuck() const { return pivotModel == 1 ? pinHolds : stuck; }   // (pivot_model 1: the pin's spring holds)
      int PinStop() const { return pinStop; }
      double PinSpringJ() const { return 0.5 * pinK * pinS * pinS; }

      // ---- the group's interface ----
      int NumDof() const override { return pinStop >= 0 || stuck ? 0 : 1; }
      int MaxDof() const override { return 1; }
      int NumCoords() const override { return 1; }
      int DofOf(int) const override { return NumDof() - 1; }
      CoordInfo Coord(int) const override
      {
         CoordInfo c;
         c.rate = rate;
         c.given = NumDof() == 0 ? gGiven : 0.0;
         c.invM = 1.0 / mass;
         c.invMFree = 1.0 / mass;
         c.reach = reach;
         c.maxMove = maxMove;
         c.blockedSide = pinStop >= 0 ? stops[pinStop].side : 0;
         return c;
      }
      double KineticJ() const override { return 0.5 * mass * rate * rate; }
      double Energy(const Ball *, int) const override { return pivotModel == 1 ? KineticJ() + StoredEnergy(q) + PinSpringJ() : KineticJ() + StoredEnergy(q); }
      void StepFacts(const Ball *, int, const Settings &, StepFactList &f) override
      {
         if (std::fabs(rate) > 1e-9)
            f.Turn(1.0, std::fabs(rate), 1);   // at most the group's turn limit per step
      }
      bool BeginStep(const Ball *, int, const double h, double *x, const int off) override
      {
         gH = h;
         for (BodyStop &s : stops)
            s.speed = s.SpeedNow();
         gGiven = pinStop >= 0 ? stops[pinStop].AfterStep(h) - stops[pinStop].at : 0.0;
         if (NumDof() == 1)
            x[off] = rate * h;
         gDrive = gViscous = gFriction = 0.0;
         gBallQ = 0.0;
         gBallF = V3();
         return true;
      }
      void Loads(const double *qf, const V3 &F) override
      {
         gBallQ = qf[0];
         gBallF = F;
      }
      bool Eval(const double *x, const int off, const double h, V3 *, V3 *, double *Res, Jac J) override
      {
         gDrive = gViscous = gFriction = 0.0;
         if (NumDof() == 0)
            return true;
         const double dq = x[off], r = dq / h, qm = q + 0.5 * dq;
         double slope = 0.0, lossF = 0.0, lossSlope = 0.0;
         const double Vp = StoredChord(q, dq, h, slope, lossF, lossSlope);
         const double cv = Viscous(), Fd = DriveForce();
         if (pivotModel == 1)
         {
            const PinMove pm = PinStep(dq, h);
            Res[off] = mass * (2.0 * dq / h - 2.0 * rate) - h * (-Vp - lossF - cv * r + Fd + pm.M);
            rowT[0] = mass * (2.0 * dq / h - 2.0 * rate); rowT[1] = h * (Vp + lossF); rowT[2] = h * cv * r; rowT[3] = -h * pm.M; rowT[4] = -h * gBallQ; rowT[5] = pm.K; rowT[6] = FrictionLimit(qm, r, gBallF); rowT[7] = pm.holds ? 1.0 : 0.0; rowT[8] = pm.s1; rowT[9] = dq;
            if (J)
               J[off][off] += 2.0 * mass / h + h * slope + h * lossSlope + cv + (pm.holds ? h * 0.5 * pm.K + pm.c : 0.0);
            gViscous = cv * r * dq + lossF * dq;
            gFriction = pm.loss;
            gDrive = Fd * dq;
            return true;
         }
         const double lim = dirF != 0 ? FrictionLimit(qm, r, gBallF) : 0.0;
         const double Q = -Vp - lossF - cv * r + Fd - dirF * lim;
         Res[off] = mass * (2.0 * dq / h - 2.0 * rate) - h * Q;
         if (J)
            J[off][off] += 2.0 * mass / h + h * slope + h * lossSlope + cv;
         gViscous = cv * r * dq + lossF * dq;
         gFriction = dirF * lim * dq;
         gDrive = Fd * dq;
         return true;
      }
      // > 0: the move carries the body past a stop (how far), or past the speed's zero with friction acting (the rate past it).
      // With no move (the group's test of where the step starts): a body leaving a stop isn't on it, and how far off it counts by its rate, so the search for the next arrival does not start at the stop it left (the jet bumper's ring does the same)
      // stops: bit i its stop i (the last bit for any beyond), bit 31 the speed's zero (pivot_model 0's friction turn)
      static unsigned StopBit(const int i) { return 1u << std::min(i, 30); }
      double StopExcess(const double *x, const int off, const unsigned mask = ~0u) const override
      {
         const bool wentBack = gExcessTurn;   // (the solved move, tested last, went back past the speed's zero)
         gExcessStop = -1;
         gExcessTurn = false;
         if (NumDof() == 0)
            return -INFINITY;
         const double dq = x[off], q1 = q + dq;
         double e = -INFINITY;
         for (int i = 0; i < (int)stops.size(); i++)
         {
            const BodyStop &s = stops[i];
            if (!s.enabled || !s.armed || !(mask & StopBit(i)))
               continue;
            const double ex = s.side * (q1 - (dq == 0.0 ? s.at : s.AfterStep(gH)));
            if (ex > e)
            {
               e = ex;
               gExcessStop = i;
            }
         }
         if (pivotModel != 1 && dirF != 0 && (mask & (1u << 31)))
         {
            // From rest the speed's zero is where it starts: an event only if the solved move goes back against the friction's way (it stays: PinOnStop). Going the friction's way there is no zero ahead, and the search for a stop must not take the start for one
            const double back = dq == 0.0 ? (rate != 0.0 ? -std::fabs(rate) : wentBack ? 0.0 : -INFINITY) : -dirF * (2.0 * dq / gH - rate);
            if (back > e && (rate != 0.0 || dq == 0.0 || back > 0.0))
            {
               e = back;
               gExcessStop = -1;
               gExcessTurn = true;
            }
         }
         gDqLast = dq;
         return e;
      }
      unsigned CrossingStops(const double *x, const int off) const override
      {
         const int stop = gExcessStop;   // (what the last StopExcess found stays for ArrivingAtStop and PinOnStop)
         const bool turn = gExcessTurn;
         const double dqLast = gDqLast;
         unsigned m = 0u;
         for (int i = 0; i < (int)stops.size() && i <= 30; i++)
         {
            gExcessTurn = turn;
            if (StopExcess(x, off, StopBit(i)) > 0.0)
               m |= StopBit(i);
         }
         gExcessTurn = turn;
         if (StopExcess(x, off, 1u << 31) > 0.0)
            m |= 1u << 31;
         gExcessStop = stop;
         gExcessTurn = turn;
         gDqLast = dqLast;
         return m;
      }
      bool ArrivingAtStop() const override
      {
         if (gExcessTurn || gExcessStop < 0)
            return true;
         const BodyStop &s = stops[gExcessStop];
         return s.side * (2.0 * gDqLast / gH - rate - s.speed) > 0.0;
      }
      // The free step passed a stop or the speed's zero, and the body is already there: it stays for this step, and the step is solved again as it is
      bool PinOnStop(const double *, int) override
      {
         if (NumDof() == 0)
            return false;
         if (gExcessTurn)
         {
            if (rate != 0.0)
               return false;
            stuck = true;
            dirF = 0;
            bk.sticks++;
            gPinnedInStep = true;
            return true;
         }
         if (gExcessStop < 0)
            return false;
         BodyStop &s = stops[gExcessStop];
         if (s.side * (q - s.at) < -kNearStop)
            return false;
         if (!s.touching)
            stopEvents.push_back({ gExcessStop, std::fabs(rate) });
         s.touching = true;
         const double r1 = 2.0 * gDqLast / gH - rate, f = mass * (r1 - rate) / gH;
         const double vRel = std::fabs(rate - s.speed) * s.lever, F = s.side * f > 0.0 ? std::fabs(f) / s.lever : 0.0;
         if (0.5 * s.EffMass(mass) * vRel * vRel > s.HeldEnergy(F))
            return false;
         pendLoss += 0.5 * mass * (rate * rate - s.speed * s.speed);
         bk.stopJ += 0.5 * mass * (rate * rate - s.speed * s.speed);
         if (s.speed != 0.0)   // what the stop gave it is the stop's work: nothing here is lost beyond the book above, whose sign the drive's share closes
         {
            pendDrive += mass * (s.speed - rate) * s.speed;
            pendLoss += mass * (s.speed - rate) * s.speed;
            bk.stopJ += mass * (s.speed - rate) * s.speed;
         }
         rate = s.speed;
         pinStop = gExcessStop;
         ReleasePin(pendLoss);
         stuck = false;
         dirF = 0;
         bk.rests++;
         gGiven = stops[pinStop].AfterStep(gH) - stops[pinStop].at;
         gPinnedInStep = true;
         return true;
      }
      // The step just solved left the body with no unknown, and its forces take it off: it moves in this step. (Not a body that the free step has just put there: that answer is already in.)
      bool FreedByStep() override
      {
         if (NumDof() == 1 || gPinnedInStep)
            return false;
         const double net = NetForce(), lim = FrictionLimit(q, rate, gBallF);
         const bool leaves = pinStop >= 0 ? stops[pinStop].side * net < -lim : std::fabs(net) > lim;
         if (!leaves)
            return false;
         pinStop = -1;
         stuck = false;
         dirF = rate != 0.0 ? (rate > 0.0 ? 1 : -1) : (lim > 0.0 ? (net > 0.0 ? 1 : -1) : 0);
         bk.leaves++;
         gGiven = 0.0;
         return true;
      }
      void Accept(const double *x, const int off, const double h, const bool onStop, double &driveJ, double &lossJ, bool &through) override
      {
         gPinnedInStep = false;
         through = false;
         PinMove pm;
         if (pivotModel == 1 && NumDof() == 1)
         {
            pm = PinStep(x[off], h);   // (the solved move's: the spring's end state and its loss)
            gFriction = pm.loss;
         }
         driveJ = pendDrive + gDrive;
         lossJ = pendLoss + gViscous + gFriction;
         pendDrive = pendLoss = 0.0;
         bk.viscousJ += gViscous;
         bk.frictionJ += gFriction;
         bk.driveJ += gDrive;
         const double r0 = rate;
         qAtStepStart = q;
         if (NumDof() == 1)
         {
            const double dq = x[off];
            rate = 2.0 * dq / h - rate;
            q += dq;
            if (pivotModel == 1)
            {
               pinS = pm.s1;
               pinK = pm.K;
               pinHolds = pm.holds;
            }
         }
         else if (gGiven != 0.0)
         {
            // carried by its driven stop: the stop pays for the stored energy and the friction; its work on the balls is the group's book
            const double dq = gGiven, V0 = StoredEnergy(q);
            const double lim = FrictionLimit(q + 0.5 * dq, rate, gBallF), cv = Viscous();
            q += dq;
            const double loss = lim * std::fabs(dq) + cv * rate * dq, d = StoredEnergy(q) - V0 + loss;
            driveJ += d;
            lossJ += loss;
            bk.driveJ += d;
            bk.frictionJ += loss;
         }
         for (BodyStop &s : stops)
         {
            s.at = s.AfterStep(h);
            s.speed = s.SpeedNow();
         }
         if (onStop)
            StopOrTurn(r0, h, driveJ, lossJ);
         UpdateStops();
         Decide();
         DrivenTouch(driveJ, lossJ);
         AfterMove();
      }
      // The dynamic state as plain bytes: for the engine's snapshot and the net's retake (the settings are not state)
      static constexpr int kMaxStops = 4;
      struct Dyn
      {
         double q, rate, pendDrive, pendLoss, qAtStepStart, gBallQ, gGiven, pinS, pinK;
         V3 gBallF;
         BodyBooks bk;
         struct { double at, target, rate, speed; bool enabled, armed, touching; } stop[kMaxStops];
         int pin, dir;
         bool stuck, pinnedInStep, pinHolds;
      };
      void Fill(Dyn &d) const
      {
         assert(stops.size() <= kMaxStops);
         d.q = q; d.rate = rate; d.pendDrive = pendDrive; d.pendLoss = pendLoss; d.qAtStepStart = qAtStepStart; d.gBallQ = gBallQ; d.gGiven = gGiven; d.gBallF = gBallF; d.bk = bk; d.pin = pinStop; d.dir = dirF; d.stuck = stuck; d.pinnedInStep = gPinnedInStep; d.pinS = pinS; d.pinK = pinK; d.pinHolds = pinHolds;
         for (size_t i = 0; i < stops.size(); i++)
         {
            const BodyStop &s = stops[i];
            d.stop[i].at = s.at; d.stop[i].target = s.target; d.stop[i].rate = s.rate; d.stop[i].speed = s.speed; d.stop[i].enabled = s.enabled; d.stop[i].armed = s.armed; d.stop[i].touching = s.touching;
         }
      }
      void Load(const Dyn &d)
      {
         q = d.q; rate = d.rate; pendDrive = d.pendDrive; pendLoss = d.pendLoss; qAtStepStart = d.qAtStepStart; gBallQ = d.gBallQ; gGiven = d.gGiven; gBallF = d.gBallF; bk = d.bk; pinStop = d.pin; dirF = d.dir; stuck = d.stuck; gPinnedInStep = d.pinnedInStep; pinS = d.pinS; pinK = d.pinK; pinHolds = d.pinHolds;
         for (size_t i = 0; i < stops.size(); i++)
         {
            BodyStop &s = stops[i];
            s.at = d.stop[i].at; s.target = d.stop[i].target; s.rate = d.stop[i].rate; s.speed = d.stop[i].speed; s.enabled = d.stop[i].enabled; s.armed = d.stop[i].armed; s.touching = d.stop[i].touching;
         }
      }
      // ---- the net (Group::EnergyStep): counted and reported, not expected to fire ----
      void SaveState() override
      {
         std::memset(&sv, 0, sizeof sv);
         Fill(sv);
         svEvents = stopEvents.size();
      }
      void RestoreState() override
      {
         Load(sv);
         stopEvents.resize(svEvents);
      }
      // after an energy step: the state it saved at the step's start, full precision (for an exact replay; rad, rad/s, the ball load in N m and N). The books are left out (counted, never read by a step)
      std::string StartNote() const override
      {
         char s[700];
         snprintf(s, sizeof s, "q %.17g rate %.17g qStart %.17g pendDrive %.17g pendLoss %.17g ballQ %.17g ballF %.17g %.17g %.17g given %.17g pin %d dir %d stuck %d pinnedInStep %d pinHolds %d pinS %.17g pinK %.17g pivotModel %d stops %d",
            sv.q, sv.rate, sv.qAtStepStart, sv.pendDrive, sv.pendLoss, sv.gBallQ, sv.gBallF.x, sv.gBallF.y, sv.gBallF.z, sv.gGiven, sv.pin, sv.dir, (int)sv.stuck, (int)sv.pinnedInStep, (int)sv.pinHolds,
            sv.pinS, sv.pinK, pivotModel, (int)stops.size());
         std::string msg = s;
         for (size_t i = 0; i < stops.size() && i < kMaxStops; i++)
         {
            snprintf(s, sizeof s, " {stop %d side %d at %.17g target %.17g rate %.17g speed %.17g enabled %d armed %d touching %d}", (int)i, stops[i].side,
               sv.stop[i].at, sv.stop[i].target, sv.stop[i].rate, sv.stop[i].speed, (int)sv.stop[i].enabled, (int)sv.stop[i].armed, (int)sv.stop[i].touching);
            msg += s;
         }
         return msg;
      }
      double EndResidual() const override { return 1e-8; }   // a noise stop can land past a steep onset (a light body's first touch): 0.45 uJ, runs/2026-10-31-gate-spinner-disc
      void ScaleSpeeds(const double f) override { rate *= f; }
      void ExplicitMove(const double *Res, const int off, const double h, double *x) const override
      {
         if (NumDof() == 1)
            x[off] -= Res[off] * h / (2.0 * mass);
      }

   protected:
      static constexpr double kNearStop = 1e-9;   // a body this close to a stop's place is at it
      double NetForce() const { return StoredForce(q) - Viscous() * rate + DriveForce() + gBallQ; }
      // The pin's stick spring over the move dq (pivot_model 1): the pin's twist s in its hole, Mindlin's patch spring K (the part's PinStiffness) with the contact's
      // damper c = shearTan sqrt(K I) beside it (steel 0.4, as the ball's shear spring), stepped as the rolling spring is (trapezoid torque on the trial s0 + dq, the
      // damper on the mid-step rate). Spring and damper hold while their torque is within the static limit, else they slide at dynRatio times it.
      // A linear spring under the trapezoid rule gives back what it took; the damper loses c ds^2 / h, slip |M| |trial - s1|, a softer spring (K falls, or rises with
      // s scaled by K0 / K1) the stored energy it can no longer hold.
      static constexpr double kPinShearTan = 0.4;
      struct PinMove { double M = 0.0, s1 = 0.0, K = 0.0, c = 0.0, loss = 0.0; bool holds = false; };
      // the wire row's terms at the last evaluation (pivot_model 1), each as it enters Res (N m s): inertia, stored (spring + gravity), viscous, pin, balls; the pin's K, static limit, holds, s1; dq
      double rowT[10] = {};
      std::string RowNote() const override
      {
         char line[300];
         snprintf(line, sizeof line, "row: inertia %.4g stored %.4g visc %.3g pin %.4g balls %.4g (sum %.3g) | pin K %.4g lim %.4g holds %d s1 %.6g dq %.4g gBallF %.5g",
                  rowT[0], rowT[1], rowT[2], rowT[3], rowT[4], rowT[0] + rowT[1] + rowT[2] + rowT[3] + rowT[4], rowT[5], rowT[6], (int)rowT[7], rowT[8], rowT[9], gBallF.Length());
         return line;
      }
      // diagnostics (VPX_GROUP_HISTORY_TRACE): after an accepted step (energies in uJ)
      std::string HistoryNote() const override
      {
         char line[300];
         snprintf(line, sizeof line, "q %.6f deg rate %.5g rad/s KE %.5g stored %.9g pinJ %.5g uJ pinS %.4g holds %d stuck %d pinStop %d", q * 57.29577951308232, rate, KineticJ() * 1e6, StoredEnergy(q) * 1e6,
                  PinSpringJ() * 1e6, pinS, pinHolds ? 1 : 0, stuck ? 1 : 0, pinStop);
         return line;
      }
      // diagnostics (VPX_GROUP_CONTACT_TRACE): the move x over h as this round evaluated it (read only; torques in N m, the ball load as the last Loads gave it), its stops, its residual
      std::string ContactNote(const double *x, const int off, const double h, const double *Res) const override
      {
         const int dof = NumDof();
         const double dq = dof == 1 ? x[off] : 0.0, r = dq / h, qm = q + 0.5 * dq;
         char line[600];
         int n = snprintf(line, sizeof line, "q %.6f deg rate %.4g rad/s dq %.4g rad dof %d stuck %d pinStop %d pivotModel %d dirF %d | ballQ %.5g ballF %.4g N", q * 57.29577951308232, rate, dq, dof,
                          Stuck() ? 1 : 0, pinStop, pivotModel, dirF, gBallQ, gBallF.Length());
         if (dof == 1)
         {
            double slope = 0.0, lossF = 0.0, lossSlope = 0.0;
            const double Vp = StoredChord(q, dq, h, slope, lossF, lossSlope), lim = FrictionLimit(qm, r, gBallF);
            n += snprintf(line + n, sizeof line - n, " stored %.5g visc %.4g drive %.4g", -Vp - lossF, -Viscous() * r, DriveForce());
            if (pivotModel == 1)
            {
               const PinMove pm = PinStep(dq, h);
               n += snprintf(line + n, sizeof line - n, " pin %.5g lim %.5g holds %d (K %.4g s0 %.4g s1 %.4g)", pm.M, lim, pm.holds ? 1 : 0, pm.K, pinS, pm.s1);
            }
            else
               n += snprintf(line + n, sizeof line - n, " pin %.5g lim %.5g", -dirF * lim, lim);
            n += snprintf(line + n, sizeof line - n, " | res %.5g N m s (%.4g m/s)", Res[off], Res[off] / mass * reach);
         }
         std::string note = line;
         for (size_t i = 0; i < stops.size(); i++)
         {
            const BodyStop &s = stops[i];
            snprintf(line, sizeof line, " {stop %d side %d at %.5f deg excess %.4g deg en %d armed %d touching %d}", (int)i, s.side, s.at * 57.29577951308232, s.side * (q + dq - s.at) * 57.29577951308232,
                     (int)s.enabled, (int)s.armed, (int)s.touching);
            note += line;
         }
         return note;
      }
      std::string SearchNote(const double *x, const int off, const double h) const override
      {
         if (NumDof() != 1)
            return std::string();
         const double dq = x[off];
         const PinMove pm = pivotModel == 1 ? PinStep(dq, h) : PinMove();
         char line[300];
         snprintf(line, sizeof line, "q0 %.6g deg rate0 %.4g dq %.4g q1 %.6g deg rate1 %.4g pin %d | pin spring s0 %.4g s1 %.4g K %.4g M %.4g holds %d | stop excess %.4g | ball force %.4g N",
                  q * 57.29577951308232, rate, dq, (q + dq) * 57.29577951308232, 2.0 * dq / h - rate, pinStop, pinS, pm.s1, pm.K, pm.M, pm.holds ? 1 : 0, StopExcess(x, off), gBallF.Length());
         return line;
      }
      PinMove PinStep(const double dq, const double h) const
      {
         PinMove pm;
         const double qm = q + 0.5 * dq, r = dq / h, held0 = PinSpringJ();
         pm.K = PinStiffness(qm, r, gBallF);
         if (pm.K <= 0.0)
         {
            pm.loss = held0;
            return pm;
         }
         const double lim = FrictionLimit(qm, r, gBallF);
         pm.c = kPinShearTan * std::sqrt(pm.K * mass);
         const double s0 = pm.K > pinK ? pinS * (pinK / pm.K) : pinS, sTr = s0 + dq, mTr = -0.5 * pm.K * (s0 + sTr) - pm.c * dq / h;
         pm.loss = held0 - 0.5 * pm.K * s0 * s0;
         if (std::fabs(mTr) <= lim)
         {
            pm.M = mTr;
            pm.s1 = sTr;
            pm.holds = true;
            pm.loss += pm.c * dq * dq / h;
         }
         else
         {
            pm.M = mTr > 0.0 ? dynRatio * lim : -dynRatio * lim;
            pm.s1 = (-pm.M - s0 * (0.5 * pm.K - pm.c / h)) / (0.5 * pm.K + pm.c / h);   // -M = 0.5 K (s0 + s1) + c (s1 - s0) / h
            pm.loss += std::fabs(pm.M * (sTr - pm.s1)) + pm.c * (pm.s1 - s0) * (pm.s1 - s0) / h;
         }
         return pm;
      }
      // the body rests on a stop: the pin's spring lets go, its energy lost
      void ReleasePin(double &lossJ)
      {
         if (pivotModel != 1)
            return;
         lossJ += PinSpringJ();
         bk.frictionJ += PinSpringJ();
         pinS = pinK = 0.0;
         pinHolds = false;
      }

   private:
      // a stop placed with the body past it waits for the free side; a touch is one event
      void UpdateStops()
      {
         for (BodyStop &s : stops)
         {
            const double gap = s.side * (q - s.at);
            if (!s.armed && gap <= 0.0)
               s.armed = true;
            if (s.touching && gap < -kNearStop)
               s.touching = false;
         }
      }
      // The step ended on a stop (the impact law), or where the friction brought the speed to zero
      void StopOrTurn(const double r0, const double h, double &driveJ, double &lossJ)
      {
         int hit = -1;
         double best = -1e-8;
         for (int i = 0; i < (int)stops.size(); i++)
         {
            const BodyStop &s = stops[i];
            if (!s.enabled || !s.armed)
               continue;
            const double gap = s.side * (q - s.at);
            if (gap >= best && s.side * (rate - s.speed) > 0.0)
            {
               best = gap;
               hit = i;
            }
         }
         if (hit < 0)
         {
            if (gExcessTurn)
            {
               lossJ += 0.5 * mass * rate * rate;   // (the speed left at the zero the step ended on)
               bk.frictionJ += 0.5 * mass * rate * rate;
               rate = 0.0;
            }
            return;
         }
         const double f = mass * (rate - r0) / h;
         const BodyStop &s = stops[hit];
         Impact(hit, s.side * f > 0.0 ? std::fabs(f) / s.lever : 0.0, driveJ, lossJ);
      }
      // The impact on stop i at the closing speed now: rebound by the stop's own contact, or rest if it would not leave the squash its pressing force F holds
      void Impact(const int i, const double F, double &driveJ, double &lossJ)
      {
         BodyStop &s = stops[i];
         const double vr = rate - s.speed, vClose = std::fabs(vr) * s.lever;
         if (!s.touching)
            stopEvents.push_back({ i, std::fabs(rate) });
         s.touching = true;
         const double e = s.Restitution(vClose), vOut = e * vClose;
         double rNew = s.speed - e * vr;
         if (0.5 * s.EffMass(mass) * vOut * vOut <= s.HeldEnergy(F))
         {
            rNew = s.speed;
            pinStop = i;
            ReleasePin(lossJ);
            stuck = false;
            dirF = 0;
            bk.rests++;
         }
         const double work = mass * (rNew - rate) * s.speed;   // the stop's push, against its own speed
         const double loss = 0.5 * mass * rate * rate + work - 0.5 * mass * rNew * rNew;
         driveJ += work;
         lossJ += loss;
         bk.driveJ += work;
         bk.stopJ += loss;
         bk.stopHits++;
         rate = rNew;
      }
      // A driven stop that has come up against a body not moving away: the impact at once, between steps
      void DrivenTouch(double &driveJ, double &lossJ)
      {
         for (int i = 0; i < (int)stops.size(); i++)
         {
            BodyStop &s = stops[i];
            if (i == pinStop || s.speed == 0.0 || !s.enabled || !s.armed || s.side * (q - s.at) < -1e-12 || s.side * (rate - s.speed) <= 0.0)
               continue;
            pinStop = -1;
            stuck = false;
            Impact(i, 0.0, driveJ, lossJ);
            dirF = rate > 0.0 ? 1 : rate < 0.0 ? -1 : 0;
         }
      }
      // Who has an unknown in the next step, from the forces at the end of this one
      void Decide()
      {
         const double net = NetForce(), lim = FrictionLimit(q, rate, gBallF);
         const bool wasPinned = pinStop >= 0;
         if (pinStop >= 0)
         {
            const BodyStop &s = stops[pinStop];
            // its stop has arrived (the body keeps the speed it had: it flies on, and lands if it is pressed), or the force takes it off the stop
            const bool arrived = rate != 0.0 && s.SpeedNow() == 0.0;
            if (arrived || s.side * net < -lim)
            {
               pinStop = -1;
               bk.leaves++;
            }
         }
         if (pivotModel == 1)   // the pin's spring holds or slides inside the step: no state to choose here
            return;
         if (!wasPinned && stuck && std::fabs(net) > lim)
         {
            stuck = false;
            bk.leaves++;
         }
         if (pinStop < 0 && !stuck)
         {
            if (rate != 0.0)
               dirF = rate > 0.0 ? 1 : -1;
            else if (lim > 0.0 && std::fabs(net) <= lim)
            {
               stuck = true;
               dirF = 0;
               bk.sticks++;
            }
            else
               dirF = lim > 0.0 ? (net > 0.0 ? 1 : -1) : 0;
         }
         else
            dirF = 0;
      }

      int pinStop = -1, dirF = 0;
      bool stuck = false, gPinnedInStep = false, pinHolds = false;
      double pinS = 0.0, pinK = 0.0;   // the pin's spring: its twist (rad) and stiffness at the end of the last step (pivot_model 1)
      double gH = 1e-3, gGiven = 0.0, gDrive = 0.0, gViscous = 0.0, gFriction = 0.0, pendDrive = 0.0, pendLoss = 0.0;
      double gBallQ = 0.0;
      V3 gBallF;
      mutable int gExcessStop = -1;
      mutable bool gExcessTurn = false;
      mutable double gDqLast = 0.0;
      Dyn sv{};
      size_t svEvents = 0;
   };
}
