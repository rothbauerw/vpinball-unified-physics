// license:GPLv3+

// Two-body flipper model of the unified engine (plunger, coil, link, return spring, rest stop), the flipper nudge and the ball/bat rubber contact: members of
// FlipperMoverObject (declared in hitflipper.h). Kept apart from hitflipper.cpp so the stock flipper file only holds the redirects into it.

#include "core/stdafx.h"
#include "hitflipper.h"
#include "parts/flipper.h"

#include "physics/PhysicsEngine.h"
#include "parts/ball.h"
#include "physics/rubbercontact.h"
#include "physics/hitrubber.h"
#include "physics/contactbat.h"
#include "physics/contactbatpart.h"
#include "physics/contactparthost.h"
#include "physics/contactsystem.h"
#include "physics/physicsini.h"

#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>

// ---------------------------------------------------------------------------
// Two-body flipper model (unified engine).
// A direct port of runs/2026-09-28-rampup0/flipper_model.py (rothbauerw's real-data model):
//  - plunger: translates along the coil axis, pulled toward closure by the power + hold coils,
//    hard steel stopper at the closed end, hard stop at the open end.
//  - flipper (bat + crank + shaft, one rigid body): the return spring pulls on it continuously; a
//    hysteretic rubber stop acts on it below phi = 0.
//  - link: nylon link with pin slop. Loaded contact is mismatch = 0 (coil and spring both press the
//    pins against the same hole walls); free play for 0 < mismatch < SLOP; stiff damped contact beyond.
//  - wiring: series (one supply feeds power + hold; the EOS switch shorts the hold winding during the
//    stroke) or parallel.
//  - hand-off: two-body while the coil is energized; rigid one-body (plunger slaved to the flipper) when
//    the coil is off and the link is taut; back to two-body when the coil fires again.
// All dynamics are in SI units inside; converted at the VPX boundary (1 U = 0.53975 mm,
// 1 VPT = 10 ms, 1 mass unit = 80 g). Every tunable is in FlipperSettings, per flipper, set from
// the table's vpx_physics.ini ([Flippers], [Flipper:<name>]; none = the defaults), hot-reloaded on change.
// The flipper model replaces the stock flipper torque model on the unified engine (ContactSystem::Unified: [Features] unified = 1 in the table's physics ini; no ini or key = stock)
// unless [Features] flippers = 0 in the table's vpx_physics.ini; read once, when the table has loaded.
bool FlipperModelEnabled()
{
   if (!ContactSystem::Unified())
      return false;
   static int featureOn = -1;
   if (featureOn < 0)
   {
      if (g_pplayer == nullptr || g_pplayer->m_ptable == nullptr)
         return true;
      const std::filesystem::path u = PhysicsIni::FindUnified(g_pplayer->m_ptable->m_filename);
      featureOn = PhysicsIni::Active(u).flippers ? 1 : 0;
      if (featureOn == 0)
         PLOGI << "Physics settings: [Features] flippers = 0 in " << u.string() << ": stock flippers";
   }
   return featureOn == 1;
}

static constexpr double FLIPPER_MASS_UNIT_KG = 0.080;
static constexpr double FLIPPER_U_M = 1.0625 * 0.0254 / 50.0;

// (re)builds this flipper's settings from the ini file; the model state is left untouched
void FlipperMoverObject::FlipperReloadSettings(const bool force)
{
   m_bat.LoadSettings(m_pflipper->GetName(), m_pflipper->m_ptable != nullptr ? m_pflipper->m_ptable->m_filename : std::filesystem::path(), force);
}

void FlipperMoverObject::FlipperEnsureInit()
{
   if (m_bat.init)
      return;
   FlipperReloadSettings(true);
   const double dirSign = m_direction ? 1.0 : -1.0;
   m_bat.Reset(dirSign * (double)(m_angleCur - m_angleStart), m_solState);
}

// effective inertia seen by ball impacts, in VPX units: the flipper alone while the link is loose,
// flipper + the plunger's reflected mass while the link is loaded or rigid (plunger dragged along)
void FlipperMoverObject::FlipperUpdateInertia()
{
   m_inertia = (float)(m_bat.InertiaSI() / (FLIPPER_MASS_UNIT_KG * FLIPPER_U_M * FLIPPER_U_M));
}

// one SI substep of length h seconds (the bat's state and equations are in contactbat.cpp)
void FlipperMoverObject::FlipperStep(double h, const double extTorqueNm)
{
   m_bat.Step(h, extTorqueNm, m_solState);
}

void FlipperMoverObject::UpdateDisplacementsFlipper(const float dtime)
{
   FlipperEnsureInit();
   if (++m_bat.iniCheckCounter >= 500)   // ~0.5 s of physics steps
   {
      m_bat.iniCheckCounter = 0;
      FlipperReloadSettings(false);
   }
   const double dt = (double)dtime * DEFAULT_STEPTIME_S;   // VPT -> s
   if (dt <= 0.0)
      return;
   bool hitStop = false;
   if (!m_gAdvanced)   // else RubberAdvance moved the pivot with the ball
      NudgeStep(dt);
   if (m_gAdvanced)   // RubberAdvance already stepped the bat together with its ball
   {
      m_gAdvanced = false;
      m_bat.idle = false;
      hitStop = m_gAdvHitStop;
   }
   else if (m_bat.idle && m_solState == m_bat.idleSol && m_bat.rigid == m_bat.idleRigid && m_bat.phi == m_bat.idlePhi && m_bat.omega == m_bat.idleOmega && m_bat.xp == m_bat.idleXp && m_bat.vp == m_bat.idleVp
            && m_bat.iHigh == m_bat.idleIHigh && m_bat.iLow == m_bat.idleILow)
   {
      // at rest and nothing has touched it since: the substeps would only repeat the same state
   }
   else
   {
      const double phi0 = m_bat.phi, omega0 = m_bat.omega, xp0 = m_bat.xp, vp0 = m_bat.vp, iHigh0 = m_bat.iHigh, iLow0 = m_bat.iLow;
      const bool sol0 = m_solState;
      const int n = max(1, (int)ceil(dt / 2e-5));
      const double h = dt / n;
      for (int i = 0; i < n; i++)
      {
         FlipperStep(h);
         hitStop |= m_bat.plungerHitStop;
      }
      // a whole physics step in which nothing moved (angle, speed, plunger, both coil currents) and the coil did not change: the flipper is at rest
      m_bat.idle = dt >= 0.9 * PHYSICS_STEPTIME_S && !hitStop && sol0 == m_solState && fabs(m_bat.phi - phi0) < 1e-12 && fabs(m_bat.omega - omega0) < 1e-9 && fabs(m_bat.xp - xp0) < 1e-12
         && fabs(m_bat.vp - vp0) < 1e-9 && fabs(m_bat.iHigh - iHigh0) < 1e-12 && fabs(m_bat.iLow - iLow0) < 1e-12;
      if (m_bat.idle)
      {
         m_bat.idleSol = m_solState;
         m_bat.idleRigid = m_bat.rigid;
         m_bat.idlePhi = m_bat.phi;
         m_bat.idleOmega = m_bat.omega;
         m_bat.idleXp = m_bat.xp;
         m_bat.idleVp = m_bat.vp;
         m_bat.idleIHigh = m_bat.iHigh;
         m_bat.idleILow = m_bat.iLow;
      }
   }

   const double dirSign = m_direction ? 1.0 : -1.0;
   m_angleCur = m_angleStart + (float)(dirSign * m_bat.phi);
   m_angleSpeed = (float)(dirSign * m_bat.omega * DEFAULT_STEPTIME_S);                           // rad / VPT
   m_angularAcceleration = (float)(dirSign * m_bat.lastAlpha * DEFAULT_STEPTIME_S * DEFAULT_STEPTIME_S);
   FlipperUpdateInertia();
   m_angularMomentum = m_inertia * m_angleSpeed;

   // EOS = plunger reaches the coil stopper; BOS = flipper returns to rest (same events the stock model fires)
   const float anglespd = fabsf(RADTOANG(m_angleSpeed));
   if (hitStop && m_bat.set.nudgeGain > 0.0 && !m_bat.set.nudgeTarget.empty())
   {
      FlipperMoverObject *const t = NudgeFindTarget();
      if (t != nullptr && t != this)
      {
         const double sgn = t->m_hitcircleBase.center.x >= m_hitcircleBase.center.x ? 1.0 : -1.0;
         const double ampMm = m_bat.set.nudgeGain * m_bat.hitL;
         t->NudgeStart(sgn * ampMm * 1e-3, m_bat.set.nudgeRiseMs * 1e-3, m_bat.set.nudgeFallMs * 1e-3);
         PLOGI << "Nudge: " << m_pflipper->GetName() << " struck its stop at t=" << g_pplayer->m_time_msec << " ms, L " << m_bat.hitL << " kg m2/s, pivot of " << m_bat.set.nudgeTarget << " moves " << sgn * ampMm << " mm";
      }
   }
   if (hitStop && m_enableRotateEvent > 0)
   {
      m_pflipper->FireVoidEventParm(DISPID_LimitEvents_EOS, anglespd);
      m_enableRotateEvent = 0;
   }
   else if (m_enableRotateEvent < 0 && m_bat.phi <= 0.0 && m_bat.omega <= 0.0)
   {
      m_pflipper->FireVoidEventParm(DISPID_LimitEvents_BOS, anglespd);
      m_enableRotateEvent = 0;
   }
}

// ---- flipper nudge: one-shot x pulse of the pivot (the pulse is in contactbat.cpp) ----

void FlipperMoverObject::NudgeEval(const double t, double& dispM, double& velMs) const
{
   m_bat.NudgeEval(t, dispM, velMs);
}

void FlipperMoverObject::NudgeStart(const double dispM, const double riseS, const double fallS)
{
   if (m_part != nullptr)
      m_part->bat.NudgeStart(m_hitcircleBase.center.x, dispM, riseS, fallS);   // (the part moves its pivot: BatPart::Part::PivotNow)
   else
      m_bat.NudgeStart(m_hitcircleBase.center.x, dispM, riseS, fallS);
}

// Advances the pulse by dtS and moves the pivot to where it is at the end of that time.
void FlipperMoverObject::NudgeStep(const double dtS)
{
   m_bat.NudgeStep(dtS, m_hitcircleBase.center.x);
}

FlipperMoverObject *FlipperMoverObject::NudgeFindTarget()
{
   if (m_nzTarget != nullptr && m_nzTargetName == m_bat.set.nudgeTarget)
      return m_nzTarget;
   m_nzTarget = nullptr;
   m_nzTargetName = m_bat.set.nudgeTarget;
   if (g_pplayer == nullptr || g_pplayer->m_ptable == nullptr)
      return nullptr;
   m_nzTarget = g_pplayer->m_physics->FindFlipperMover(m_nzTargetName);
   if (m_nzTarget == nullptr)
      PLOGW << "Nudge: " << m_pflipper->GetName() << " nudge_target '" << m_nzTargetName << "' is not a flipper in this table";
   return m_nzTarget;
}

// Unit boundary for the rubber model: VPX units -> SI, pivot frame. The bat angle in the engine is the angle of the
// centre line measured so that angle 0 points along (0,-1) (see HitTestFlipperEnd), i.e. phi = angle - 90 deg from +x.
RubberContact::Probe FlipperMoverObject::RubberProbe(const BallS& ball) const
{
   constexpr double U = 0.53975e-3;   // m per VPX unit
   RubberContact::Outline ol;
   ol.rb = (double)m_hitcircleBase.radius * U;
   ol.re = (double)m_endradius * U;
   ol.Lc = (double)m_flipperradius * U;
   RubberContact::BallSI b;
   b.R = (double)ball.m_radius * U;
   const RubberContact::V3 c(((double)ball.m_pos.x - m_hitcircleBase.center.x) * U,
                             ((double)ball.m_pos.y - m_hitcircleBase.center.y) * U,
                             (double)ball.m_pos.z * U);   // playfield surface is z = 0
   return RubberContact::ProbeSurface(FlipperModel::MakeRubberParams(m_bat.set), b, ol, (double)m_angleCur - M_PI / 2.0, c);
}

// A ball may engage when rubber_contact is on for this flipper, no other ball is engaged here, and the ball is not
// already driven by another flipper. Anything else falls back to the impulse path.
bool FlipperMoverObject::RubberWants(const HitBall *ball) const
{
   return FlipperModelEnabled() && m_bat.set.rubberContact
      && g_pplayer->m_time_msec >= m_rbCooldownUntilMs
      && (MultiBall() || m_rbBall == nullptr || m_rbBall == ball)
      && ball->m_bandOwner == nullptr   // a force-model contact is handed over at the moment of contact (Collide / Contact release it)
      && (ball->m_rubberOwner == nullptr || ball->m_rubberOwner == this);
}

bool FlipperMoverObject::RubberSideHas(const ContactHost *r) const
{
   if (m_rbSide)
      for (const SurfaceContact &it : m_rbSide->items)
         if (it.r == r)
            return true;
   return false;
}

// Unified engine: any number of balls can be on one bat (the contact system steps the bat once with all their torques); VPX_MULTIBAT=0 restores one ball per bat
bool FlipperMoverObject::MultiBall()
{
   return ContactSystem::Unified();
}

bool FlipperMoverObject::RubberIsEngaged(const BallS& ball) const
{
   if (MultiBall())
      return ball.m_owner != nullptr && ball.m_owner->m_rubberOwner == this;
   return m_rbBall != nullptr && &m_rbBall->m_d == &ball;
}

void FlipperMoverObject::RubberEngage(HitBall *ball)
{
   m_rbBall = ball;
   m_rbClearSteps = 0;
   m_rbState.Reset();
   m_rbLastDist = 0.0;
   m_rbJustEngaged = true;
   ball->m_rubberOwner = this;
   m_rbSide = std::make_shared<RubberSideSet>();   // the rubbers around the ball are part of the step too (a post beside the base)
   {
      constexpr double U = 0.53975e-3;
      HitRubber::SideGather(ball->m_d, RubberContact::V3((double)ball->m_d.m_pos.x * U, (double)ball->m_d.m_pos.y * U, (double)ball->m_d.m_pos.z * U), (double)ball->m_d.m_radius * U, 0.005, *m_rbSide);
   }
   PLOGI << "Rubber: " << m_pflipper->GetName() << " engaged ball at t=" << g_pplayer->m_time_msec << " ms, pos (" << ball->m_d.m_pos.x << ',' << ball->m_d.m_pos.y << ',' << ball->m_d.m_pos.z << ") vel (" << ball->m_d.m_vel.x << ',' << ball->m_d.m_vel.y << ',' << ball->m_d.m_vel.z << ") U/T";
}

FlipperMoverObject::~FlipperMoverObject()
{
   RubberRelease();   // (a bat that is a moving part: its host lets go of the balls when it is destroyed, ~ContactHost)
}

void FlipperMoverObject::RubberRelease()
{
   if (m_rbBall != nullptr)
   {
      PLOGI << "Rubber: " << m_pflipper->GetName() << " released ball at t=" << g_pplayer->m_time_msec << " ms, vel (" << m_rbBall->m_d.m_vel.x << ',' << m_rbBall->m_d.m_vel.y << ',' << m_rbBall->m_d.m_vel.z << ") U/T";
      m_rbBall->m_rubberOwner = nullptr;
   }
   m_rbSide.reset();
   m_rbBall = nullptr;
   m_rbClearSteps = 0;
}

bool FlipperMoverObject::RubberReleaseStep(const bool clearAndReceding)
{
   if (m_rbBall == nullptr)
      return false;
   m_rbClearSteps = clearAndReceding ? m_rbClearSteps + 1 : 0;
   if (m_rbClearSteps < (int)m_bat.set.rubberReleaseSteps)
      return false;
   RubberRelease();
   return true;
}

// Diagnostics: when VPX_RUBBER_LOG names a file, the LeftFlipper (or the one named by VPX_RUBBER_LOG_FLIPPER) writes one CSV row per physics step with its own state, the first
// ball's state and what the rubber model did since the previous row (engaged or not, so rubber_contact 0 and 1 compare directly).
void FlipperMoverObject::RubberLogRow()
{
   static FILE *const f = []() -> FILE * {
      const char *const p = getenv("VPX_RUBBER_LOG");
      FILE *const fp = (p != nullptr && *p) ? fopen(p, "w") : nullptr;
      if (fp != nullptr)
         fprintf(fp, "t_ms,lf_angle_deg,lf_omega_rad_s,ball_x,ball_y,ball_z,ball_vx,ball_vy,ball_vz,ball_Lx,ball_Ly,ball_Lz,engaged,contact_steps,peak_fn_n,peak_delta_mm\n");
      return fp;
   }();
   static const std::string logFlipper = []() { const char *const n = getenv("VPX_RUBBER_LOG_FLIPPER"); return std::string((n != nullptr && *n) ? n : "LeftFlipper"); }();
   if (f == nullptr || g_pplayer->m_vball.empty() || m_pflipper->GetName() != logFlipper)
      return;
   const HitBall *const b = &g_pplayer->m_vball[0]->m_hitBall;
   fprintf(f, "%u,%.5f,%.5f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%d,%d,%.5f,%.6f\n", g_pplayer->m_time_msec,
      (double)m_angleCur * 180.0 / M_PI, m_bat.omega * (m_direction ? 1.0 : -1.0),
      b->m_d.m_pos.x, b->m_d.m_pos.y, b->m_d.m_pos.z, b->m_d.m_vel.x, b->m_d.m_vel.y, b->m_d.m_vel.z,
      b->m_angularmomentum.x, b->m_angularmomentum.y, b->m_angularmomentum.z,
      (m_rbBall == b) ? 1 : 0, m_rbLogSteps, m_rbLogPeakFn, m_rbLogPeakDelta * 1e3);
   fflush(f);
   m_rbLogSteps = 0;
   m_rbLogPeakFn = m_rbLogPeakDelta = 0.0;
}

// Coupled bat + ball step while a ball is engaged (called before any mover is displaced in this cycle, so the ball's own
// straight-line move is skipped via HitBall::m_rubberAdvanced whatever the mover order). Gravity, playfield and everything
// else the engine does to the ball is unchanged; only the contact with this flipper's rubber is integrated here, in SI:
// per substep the rubber force acts on the ball, the equal and opposite torque goes into the bat (FlipperStep), then both move.
// Returns true when it advanced the pair.
bool FlipperMoverObject::RubberAdvance(const float dtime)
{
   if (m_rbBall == nullptr)
      return false;
   FlipperEnsureInit();
   if (!m_bat.set.rubberContact)   // switched off in the ini while a ball is engaged
   {
      RubberRelease();
      return false;
   }
   const double dt = (double)dtime * DEFAULT_STEPTIME_S;   // VPT -> s
   if (dt <= 0.0)
      return false;

   constexpr double U = 0.53975e-3;                          // m per VPX unit
   constexpr double VS = U / DEFAULT_STEPTIME_S;             // VPX velocity (U/VPT) -> m/s
   HitBall *const pb = m_rbBall;
   const double dirSign = m_direction ? 1.0 : -1.0;
   const RubberContact::RubberParams rp = FlipperModel::MakeRubberParams(m_bat.set);
   RubberContact::Outline ol;
   ol.rb = (double)m_hitcircleBase.radius * U;
   ol.re = (double)m_endradius * U;
   ol.Lc = (double)m_flipperradius * U;
   RubberContact::BallSI bs;
   bs.R = (double)pb->m_d.m_radius * U;
   bs.m = FLIPPER_MASS_UNIT_KG * (double)pb->m_d.m_mass;
   const double inertiaU = (double)pb->Inertia();            // VPX units; spin rate (rad/VPT) = angular momentum / inertiaU

   const double px = m_hitcircleBase.center.x, py = m_hitcircleBase.center.y;
   const Vertex3Ds posStart = pb->m_d.m_pos;
   RubberContact::V3 c(((double)posStart.x - px) * U, ((double)posStart.y - py) * U, (double)posStart.z * U);
   RubberContact::V3 v((double)pb->m_d.m_vel.x * VS, (double)pb->m_d.m_vel.y * VS, (double)pb->m_d.m_vel.z * VS);
   const double wScale = 1.0 / (inertiaU * DEFAULT_STEPTIME_S);   // angular momentum -> rad/s
   RubberContact::V3 w((double)pb->m_angularmomentum.x * wScale, (double)pb->m_angularmomentum.y * wScale, (double)pb->m_angularmomentum.z * wScale);

   const double hMax = std::min(2e-5, m_bat.set.rubberSubstepUs * 1e-6);
   const int n = max(1, (int)ceil(dt / hMax));
   const double h = dt / n;
   bool hitStop = false;
   bool live = true;   // false once the ball is released mid-cycle: it then just flies on while the bat keeps stepping

   // Safety net. A ball that starts the engagement already buried in the rubber would release the stored compression as speed; a non-finite state
   // or a runaway speed means the step left the model's range. Either way the ball is left untouched (the engine moves it as usual), the
   // engagement ends, and the impulse path takes over for a short while. The bat keeps stepping whatever happens.
   bool tripped = false;
   const double vCap = std::max(50.0, 4.0 * v.Length() + 5.0);   // m/s
   RubberContact::V3 absCorr;   // position correction for the ball at the end of the step when the start was rewound against a moving bat
   if (m_rbJustEngaged)
   {
      m_rbJustEngaged = false;
      // The ball starts the engagement inside the rubber: the hit was found late. Go back to the instant it first touched: in the bat's frame, along
      // the relative velocity, by depth / closing speed, at most one physics step. Moving the ball along its relative velocity puts it where it was
      // when it touched the bat's earlier pose; the ball's own absolute path is then restored at the end by the bat's displacement over that time.
      // If it is not closing (pushed in by something else) it is seated at touching distance along the normal. Nothing is stored as compression.
      constexpr double kMaxRewindS = PHYSICS_STEPTIME_S;
      const double phi0 = (double)m_angleStart + dirSign * m_bat.phi - M_PI / 2.0;
      const RubberContact::V3 ez(0.0, 0.0, 1.0);
      RubberContact::V3 vsLast;
      double rewoundS = 0.0;
      for (int pass = 0; pass < 3; pass++)
      {
         const RubberContact::Probe pr = RubberContact::ProbeSurface(rp, bs, ol, phi0, c);
         if (!pr.overlap || pr.delta <= 1e-6 || pr.dist <= 1e-12)
            break;
         const RubberContact::V3 nrm = (c - pr.point) / pr.dist;
         const RubberContact::V3 vs = ez.Cross(pr.point) * (dirSign * m_bat.omega);   // bat surface velocity at the contact point
         const RubberContact::V3 vrel = v - vs;
         const double closing = -vrel.Dot(nrm);
         if (closing > 0.01)
         {
            const double t = std::min(pr.delta / closing, kMaxRewindS);
            c = c - vrel * t;
            absCorr = absCorr - vs * t;
            vsLast = vs;
            rewoundS += t;
         }
         else
            c = c + nrm * pr.delta;   // seat
      }
      if (rewoundS > 0.0)
         PLOGI << "Rubber: " << m_pflipper->GetName() << " rewound the ball " << rewoundS * 1e3 << " ms to its first touch";
   }
   const RubberContact::V3 pivotW(px * U, py * U, 0.0);
   if (m_rbSide)   // the force-model rubbers within reach of the ball this cycle
      HitRubber::SideGather(pb->m_d, c + pivotW, bs.R, 2.0 * v.Length() * dt + 0.005, *m_rbSide);
   for (int i = 0; i < n; i++)
   {
      double ctq = 0.0;
      // pivot velocity (flipper nudge): c is relative to the pivot, so the contact sees the ball's velocity relative to it
      double nzD = 0.0, nzV = 0.0;
      if (m_bat.nzActive)
         NudgeEval(m_bat.nzT + (i + 0.5) * h, nzD, nzV);
      const RubberContact::V3 vRel = v - RubberContact::V3(nzV, 0.0, 0.0);
      if (live)
      {
         const double phiMath = (double)m_angleStart + dirSign * m_bat.phi - M_PI / 2.0;   // engine angle 0 points along (0,-1)
         const RubberContact::ContactResult r = RubberContact::Evaluate(rp, bs, ol, phiMath, dirSign * m_bat.omega, c, vRel, w, h, m_rbState);
         if (r.inContact)
         {
            m_rbLogSteps++;
            m_rbLogPeakFn = std::max(m_rbLogPeakFn, r.fn);
            m_rbLogPeakDelta = std::max(m_rbLogPeakDelta, r.delta);
            v = v + r.force * (1.0 / bs.m) * h;
            w = w + r.torqueOnBall * (1.0 / bs.Inertia()) * h;
            ctq = dirSign * r.torqueOnBat;   // FlipperStep wants the torque in the closing direction
            if (r.delta > 0.8 * bs.R)        // sunk deep into the rubber
               tripped = true;
         }
         if (m_rbSide && !m_rbSide->items.empty())
         {
            // rubbers beside the bat (a post by the base): same ball, same substep. c is relative to the pivot, which the nudge may be moving.
            RubberContact::V3 sF, sM;
            const double pxW = m_bat.nzActive ? m_bat.nzBaseX * U + nzD : px * U;
            HitRubber::SideStep(*m_rbSide, bs, RubberContact::V3(c.x + pxW, c.y + py * U, c.z), v, w, h, sF, sM);
            v = v + sF * (1.0 / bs.m) * h;
            w = w + sM * (1.0 / bs.Inertia()) * h;
            for (const double d : m_rbSide->depthM)
               if (d > 0.8 * bs.R)   // sunk deep into a rubber
                  tripped = true;
         }
         const bool receding = r.dist > m_rbLastDist;
         m_rbLastDist = r.dist;
         // Resting regime: a ball hovering a few micrometres off the rubber is still engaged. It is let go only once it is clear by more than
         // rubber_release_gap_mm as well as receding (release hysteresis), else it is released and re-hit by the engine dozens of times a second.
         const bool clearOfRubber = !r.inContact && r.dist - bs.R > m_bat.set.rubberReleaseGapMm * 1e-3;
         if (RubberReleaseStep(clearOfRubber && receding))
            live = false;
      }
      c = c + (v - RubberContact::V3(nzV, 0.0, 0.0)) * h;
      if (live && (!std::isfinite(c.x + c.y + c.z + v.x + v.y + v.z + w.x + w.y + w.z) || v.Length() > vCap))
         tripped = true;
      if (tripped)
      {
         live = false;
         ctq = 0.0;   // do not feed a bad contact torque into the bat
      }
      FlipperStep(h, ctq);
      hitStop |= m_bat.plungerHitStop;
   }
   NudgeStep(dt);   // the pivot ends the cycle where the substeps above assumed

   if (tripped)
   {
      PLOGW << "Rubber: " << m_pflipper->GetName() << " contact step left the model's range, using the impulse path (ball untouched this cycle)";
      RubberRelease();
      m_rbCooldownUntilMs = g_pplayer->m_time_msec + 20;
      m_gAdvanced = true;   // the bat has been stepped; the ball is left for the engine to move
      m_gAdvHitStop = hitStop;
      return true;
   }

   c = c + absCorr;
   if (m_rbSide && ContactSystem::EnvVar("VPX_CONTACT_TRACE", "VPX_RING_TRACE") && !m_rbSide->items.empty())
   {
      std::ostringstream o;
      for (size_t i = 0; i < m_rbSide->items.size(); i++)
         o << ' ' << m_rbSide->items[i].r->Name() << ':' << m_rbSide->depthM[i] * 1e6 << "um";
      PLOGI << "FlipperSide " << m_pflipper->GetName() << " t=" << g_pplayer->m_time_msec << " pos " << c.x / U + m_hitcircleBase.center.x << ' ' << c.y / U + py << ' ' << c.z / U << " v " << v.x << ' ' << v.y << ' ' << v.z << " depth" << o.str();
   }
   pb->m_d.m_pos =Vertex3Ds((float)(c.x / U + m_hitcircleBase.center.x), (float)(c.y / U + py), (float)(c.z / U));
   pb->m_d.m_vel = Vertex3Ds((float)(v.x / VS), (float)(v.y / VS), (float)(v.z / VS));
   const double lScale = inertiaU * DEFAULT_STEPTIME_S;
   pb->m_angularmomentum = Vertex3Ds((float)(w.x * lScale), (float)(w.y * lScale), (float)(w.z * lScale));
   const Vertex3Ds ds = pb->m_d.m_pos - posStart;
   pb->m_lastEventSqrDist += ds.LengthSquared();
#ifdef C_DYNAMIC
   pb->m_drsq = ds.LengthSquared();
#endif
   pb->m_rubberAdvanced = true;   // the ball has been moved for the whole cycle, released or not
   m_gAdvanced = true;
   m_gAdvHitStop = hitStop;
   return true;
}

// ---------------------------------------------------------------------------
// The bat as a moving part of the contact system (unified engine, rubber_contact on; docs/UnifiedContactDesign.txt 12.2). The part (contactbatpart.h) steps the bat with its
// balls; this adapter is the table side: the state the script and the renderer read, the stroke events, the nudge of the other bat, the flipper's Hit event.

namespace
{
   // The host's contact settings from the bat's rubber keys: the release gap and the material the contact system's own step rule reads (the contact itself is the part's)
   RubberObjectSettings FlipperHostSettings(const FlipperSettings &f)
   {
      RubberObjectSettings s;
      s.contact = true;
      s.surfaceRubber = true;
      s.eMpa = f.rubberEMpa;
      s.layerMm = f.rubberThicknessMm;
      s.betaLayer = f.rubberBetaLayer;
      s.alphaSPerM = f.rubberAlphaSPerM;
      s.mu = f.rubberMu;
      s.shearScale = f.rubberShearScale;
      s.twist = f.rubberTwist;
      s.tanDelta = f.rubberTanDelta;
      s.normalLoss = false;   // (MakeRubberParams: the bat rubber keeps its tuned rebound)
      s.releaseGapMm = f.rubberReleaseGapMm;
      return s;
   }

   class FlipperPartAdapter final : public PartAdapter
   {
   public:
      FlipperPartAdapter(FlipperMoverObject *mover, HitFlipper *flipper) : m_mover(mover), m_flipper(flipper) {}
      ContactPartHost *m_host = nullptr;

      void Publish() override
      {
         FlipperMoverObject &m = *m_mover;
         BatPart::Part &p = *m.m_part;
         constexpr double U = ContactSystem::kMetresPerUnit;
         // the bat's ini keys are read again about every 0.5 s (as the flipper's own stepping did); a change reconfigures the part and the host
         if (++m_reloadCounter >= 500)
         {
            m_reloadCounter = 0;
            const auto before = p.bat.iniTime;
            p.bat.LoadSettings(m.m_pflipper->GetName(), m.m_pflipper->m_ptable != nullptr ? m.m_pflipper->m_ptable->m_filename : std::filesystem::path(), false);
            if (p.bat.iniTime != before)
            {
               p.Configure();
               m.m_bat.set = p.bat.set;   // (the mover's copy: the nudge target, RubberProbe)
               m_host->SetContactSettings(FlipperHostSettings(p.bat.set));
            }
         }
         m_host->m_enabled = m.m_enabled;   // the flipper's Enabled: the balls touch the bat or not
         // the state the script and the renderer read
         const double dirSign = p.dirSign;
         m.m_angleCur = m.m_angleStart + (float)(dirSign * p.bat.phi);
         m.m_angleSpeed = (float)(dirSign * p.bat.omega * DEFAULT_STEPTIME_S);   // rad / VPT
         m.m_angularAcceleration = (float)(dirSign * p.bat.lastAlpha * DEFAULT_STEPTIME_S * DEFAULT_STEPTIME_S);
         m.m_inertia = (float)(p.bat.InertiaSI() / (FLIPPER_MASS_UNIT_KG * FLIPPER_U_M * FLIPPER_U_M));
         m.m_angularMomentum = m.m_inertia * m.m_angleSpeed;
         {
            double vx;
            m.m_hitcircleBase.center.x = (float)(p.PivotNow(vx).x / U);   // the pivot during a nudge pulse (the renderer draws the bat there)
         }
         // EOS = the plunger reached the coil stopper; BOS = the bat is back at rest; the striking bat's stop hit nudges the other bat's pivot
         const bool hitStop = p.hitStopLatch;
         p.hitStopLatch = false;
         const float anglespd = fabsf(RADTOANG(m.m_angleSpeed));
         if (hitStop && p.bat.set.nudgeGain > 0.0 && !p.bat.set.nudgeTarget.empty())
         {
            FlipperMoverObject *const t = m.NudgeFindTarget();
            if (t != nullptr && t != &m)
            {
               const double sgn = t->m_hitcircleBase.center.x >= m.m_hitcircleBase.center.x ? 1.0 : -1.0;
               const double ampMm = p.bat.set.nudgeGain * p.bat.hitL;
               t->NudgeStart(sgn * ampMm * 1e-3, p.bat.set.nudgeRiseMs * 1e-3, p.bat.set.nudgeFallMs * 1e-3);
               PLOGI << "Nudge: " << m.m_pflipper->GetName() << " struck its stop at t=" << g_pplayer->m_time_msec << " ms, L " << p.bat.hitL << " kg m2/s, pivot of " << p.bat.set.nudgeTarget << " moves " << sgn * ampMm << " mm";
            }
         }
         if (hitStop && m.m_enableRotateEvent > 0)
         {
            m.m_pflipper->FireVoidEventParm(DISPID_LimitEvents_EOS, anglespd);
            m.m_enableRotateEvent = 0;
         }
         else if (m.m_enableRotateEvent < 0 && p.bat.phi <= 0.0 && p.bat.omega <= 0.0)
         {
            m.m_pflipper->FireVoidEventParm(DISPID_LimitEvents_BOS, anglespd);
            m.m_enableRotateEvent = 0;
         }
      }

      // the bat's Hit (a ball's contact began: its closing speed, m/s) as the flipper's script Hit / Collide event and the haptics
      void Fire(const std::vector<ContactMovingPart::Event> &events) override
      {
         for (const ContactMovingPart::Event &e : events)
            if (e.kind == ContactMovingPart::Event::Hit)
            {
               const float bnv = -(float)(e.value / ContactSystem::kVelocityScale);   // engine velocity units, < 0 approaching
               m_flipper->FireHitEventsAndRumble(bnv, false, -bnv);
            }
      }

      bool FixedSettings() const override { return true; }   // (the host has no ini of its own: Publish reloads the bat's keys and passes them on)

   private:
      FlipperMoverObject *m_mover;
      HitFlipper *m_flipper;
      int m_reloadCounter = 0;
   };
}

void FlipperMoverObject::SetPartSolenoid(const bool s)
{
   m_part->solenoid = s;
}

ContactPartHost *FlipperMoverObject::MakePart(HitFlipper *const owner, const float zlow, const float zhigh)
{
   if (!FlipperModelEnabled())
      return nullptr;
   constexpr double U = ContactSystem::kMetresPerUnit;
   auto part = std::make_shared<BatPart::Part>();
   part->bat.LoadSettings(m_pflipper->GetName(), m_pflipper->m_ptable != nullptr ? m_pflipper->m_ptable->m_filename : std::filesystem::path(), true);
   if (!part->bat.set.rubberContact)
      return nullptr;   // the impulse path: the flipper keeps its own stepping
   const double dirSign = m_direction ? 1.0 : -1.0;
   part->bat.Reset(dirSign * (double)(m_angleCur - m_angleStart), m_solState);
   part->pivot = RubberContact::V3((double)m_hitcircleBase.center.x * U, (double)m_hitcircleBase.center.y * U, 0.0);   // (z = 0: the rubber heights are from the playfield)
   part->ol.rb = (double)m_hitcircleBase.radius * U;
   part->ol.re = (double)m_endradius * U;
   part->ol.Lc = (double)m_flipperradius * U;
   part->angleStart = (double)m_angleStart;
   part->dirSign = dirSign;
   part->solenoid = m_solState;
   part->Configure();
   m_bat.set = part->bat.set;   // (the mover's copy: the nudge target, RubberProbe)
   m_part = part;
   // the box the bat can sweep: the circle of its reach about the pivot (VPX units)
   const float reachU = (float)(part->TipReach() / U) + 1.f;
   const float lo[3] = { m_hitcircleBase.center.x - reachU, m_hitcircleBase.center.y - reachU, zlow };
   const float hi[3] = { m_hitcircleBase.center.x + reachU, m_hitcircleBase.center.y + reachU, zhigh };
   auto adapter = std::make_unique<FlipperPartAdapter>(this, owner);
   FlipperPartAdapter *const a = adapter.get();
   // (no ini path for the host: the adapter reads the bat's keys and passes the host its settings)
   ContactPartHost *const host = new ContactPartHost(m_pflipper, part, std::move(adapter), FlipperHostSettings(part->bat.set), std::filesystem::path(), m_pflipper->GetName(), "flipper", lo, hi);
   a->m_host = host;
   host->m_enabled = m_enabled;
   PLOGI << "Flipper " << m_pflipper->GetName() << ": unified part (bat inertia " << part->bat.der.iFlipper << " kg m^2, reach " << part->TipReach() * 1e3 << " mm, rubber e_mpa "
         << part->bat.set.rubberEMpa << ", largest acceleration " << part->alphaMax << " rad/s^2)";
   return host;
}
