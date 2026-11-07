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
#include "ui/live/LiveUI.h"

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
// The flipper model is the unified engine's only flipper (ContactSystem::Unified: [Features] unified = 1 in the table's physics ini; no ini or key = a stock table). The model itself is batgroup.h.

static constexpr double FLIPPER_MASS_UNIT_KG = 0.080;
static constexpr double FLIPPER_U_M = 1.0625 * 0.0254 / 50.0;

// ---- flipper nudge: one-shot x pulse of the pivot (the pulse is in coilmodel.h) ----

void FlipperMoverObject::NudgeStart(const double dispM, const double riseS, const double fallS)
{
   if (m_part != nullptr)
      m_part->NudgeStart(dispM, riseS, fallS);   // (the part moves its pivot)
}

FlipperMoverObject *FlipperMoverObject::NudgeFindTarget()
{
   if (m_nzTarget != nullptr && m_nzTargetName == m_part->bat.set.nudgeTarget)
      return m_nzTarget;
   m_nzTarget = nullptr;
   m_nzTargetName = m_part->bat.set.nudgeTarget;
   if (g_pplayer == nullptr || g_pplayer->m_ptable == nullptr)
      return nullptr;
   m_nzTarget = g_pplayer->m_physics->FindFlipperMover(m_nzTargetName);
   if (m_nzTarget == nullptr)
      PLOGW << "Nudge: " << m_pflipper->GetName() << " nudge_target '" << m_nzTargetName << "' is not a flipper in this table";
   return m_nzTarget;
}

FlipperMoverObject::~FlipperMoverObject() = default;

// ---------------------------------------------------------------------------
// The bat as a moving part of the contact system, on the group step (unified engine; docs/unified-physics/archive/UnifiedContactDesign.txt 12.2). The part (contactbatpart.h) steps the bat with its
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
      ~FlipperPartAdapter() override { PLOGI << "Group step: energy steps alone, table end: " << m_name << " " << m_esSeen; }
      ContactPartHost *m_host = nullptr;
      std::string m_name;   // (kept: the flipper may be gone when the adapter goes)
      long m_esSeen = 0;
      int m_esLines = 0;

      void Publish() override
      {
         FlipperMoverObject &m = *m_mover;
         BatPart::Part &p = *m.m_part;
         constexpr double U = ContactSystem::kMetresPerUnit;
         p.bookMs = (int)g_pplayer->m_time_msec;
         // the bat's ini keys are read again about every 0.5 s (as the flipper's own stepping did); a change reconfigures the part and the host
         if (++m_reloadCounter >= 500)
         {
            m_reloadCounter = 0;
            const auto before = p.bat.iniTime;
            p.bat.LoadSettings(m.m_pflipper->GetName(), m.m_pflipper->m_ptable != nullptr ? m.m_pflipper->m_ptable->m_filename : std::filesystem::path(), false);
            if (p.bat.iniTime != before)
            {
               p.Configure();   // (the state is kept)
               m_host->SetContactSettings(FlipperHostSettings(p.bat.set));
            }
         }
         m_host->m_enabled = m.m_enabled;   // the flipper's Enabled: the balls touch the bat or not
         // the state the script and the renderer read
         const BatGroup::Part &g = p.gp;
         const double dirSign = m.m_direction ? 1.0 : -1.0;
         m.m_angleCur = m.m_angleStart + (float)(dirSign * g.phi);
         m.m_angleSpeed = (float)(dirSign * g.omega * DEFAULT_STEPTIME_S);   // rad / VPT
         m.m_angularAcceleration = (float)(dirSign * g.aPhi * DEFAULT_STEPTIME_S * DEFAULT_STEPTIME_S);
         m.m_inertia = (float)(g.InertiaRead() / (FLIPPER_MASS_UNIT_KG * FLIPPER_U_M * FLIPPER_U_M));
         m.m_angularMomentum = m.m_inertia * m.m_angleSpeed;
         m.m_hitcircleBase.center.x = (float)(g.geo.pivot.x / U);   // the pivot during a nudge pulse (the renderer draws the bat there)
         // EOS = the plunger reached the coil stopper; BOS = the bat is back at rest; the striking bat's stop hit nudges the other bat's pivot
         const bool hitStop = g.hitStopLatch;
         p.gp.hitStopLatch = false;
         const float anglespd = fabsf(RADTOANG(m.m_angleSpeed));
         if (hitStop && p.bat.set.nudgeGain > 0.0 && !p.bat.set.nudgeTarget.empty())
         {
            FlipperMoverObject *const t = m.NudgeFindTarget();
            if (t != nullptr && t != &m)
            {
               const double sgn = t->m_hitcircleBase.center.x >= m.m_hitcircleBase.center.x ? 1.0 : -1.0;
               const double ampMm = p.bat.set.nudgeGain * g.hitL;
               t->NudgeStart(sgn * ampMm * 1e-3, p.bat.set.nudgeRiseMs * 1e-3, p.bat.set.nudgeFallMs * 1e-3);
               PLOGI << "Nudge: " << m.m_pflipper->GetName() << " struck its stop at t=" << g_pplayer->m_time_msec << " ms, L " << g.hitL << " kg m2/s, pivot of " << p.bat.set.nudgeTarget << " moves " << sgn * ampMm << " mm";
            }
         }
         if (hitStop && m.m_enableRotateEvent > 0)
         {
            m.m_pflipper->FireVoidEventParm(DISPID_LimitEvents_EOS, anglespd);
            m.m_enableRotateEvent = 0;
         }
         else if (m.m_enableRotateEvent < 0 && g.phi <= 0.0 && g.omega <= 0.0)
         {
            m.m_pflipper->FireVoidEventParm(DISPID_LimitEvents_BOS, anglespd);
            m.m_enableRotateEvent = 0;
         }
         if (p.energySteps > m_esSeen)   // rule 8: a step of the bat alone was taken on energy (no ball: nothing but the bat in the group)
         {
            m_esSeen = p.energySteps;
            if (m_esLines < 5)
            {
               m_esLines++;
               PLOGW << "Group step: energy step alone, part " << m.m_pflipper->GetName() << ", balls 0, t=" << g_pplayer->m_time_msec << " ms, " << m_esSeen << " so far, why not known for a part alone; " << p.StateNote() << ". A net, not expected to fire: keep vpinball.log and report it";
            }
         }
         m_name = m.m_pflipper->GetName();
         if (p.failedSteps > 0 && !m_warned)   // rule 8
         {
            m_warned = true;
            const std::string msg = std::string("Flipper ") + m.m_pflipper->GetName() + ": a step of the bat alone could not be solved: the rest of that cycle was lost";
            PLOGW << msg;
            if (g_pplayer->m_liveUI != nullptr)
               g_pplayer->m_liveUI->PushNotification(msg, 6000);
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
      bool m_warned = false;
   };
}

void FlipperMoverObject::SetPartSolenoid(const bool s)
{
   m_part->SetCoil(s);
}

ContactPartHost *FlipperMoverObject::MakePart(HitFlipper *const owner, const float zlow, const float zhigh)
{
   constexpr double U = ContactSystem::kMetresPerUnit;
   auto part = std::make_shared<BatPart::Part>();
   part->bat.LoadSettings(m_pflipper->GetName(), m_pflipper->m_ptable != nullptr ? m_pflipper->m_ptable->m_filename : std::filesystem::path(), true);
   // A bat has one path. A table that asks for another is told it is ignored (rule 8).
   const std::filesystem::path iniU = m_pflipper->m_ptable != nullptr ? PhysicsIni::FindUnified(m_pflipper->m_ptable->m_filename) : std::filesystem::path();
   if (!part->bat.set.rubberContact || !PhysicsIni::Active(iniU).flippers)
   {
      const std::string msg = std::string("Flipper ") + m_pflipper->GetName() + ": rubber_contact = 0, flipper_rubber = 0 or flippers = 0 is ignored: a bat has no impulse path and no stock path on the unified engine";
      PLOGW << msg;
      if (g_pplayer != nullptr && g_pplayer->m_liveUI != nullptr)
         g_pplayer->m_liveUI->PushNotification(msg, 6000);
   }
   const double dirSign = m_direction ? 1.0 : -1.0;
   part->pivot = RubberContact::V3((double)m_hitcircleBase.center.x * U, (double)m_hitcircleBase.center.y * U, 0.0);   // (z = 0: the rubber heights are from the playfield)
   part->gp.geo.ol.rb = (double)m_hitcircleBase.radius * U;
   part->gp.geo.ol.re = (double)m_endradius * U;
   part->gp.geo.ol.Lc = (double)m_flipperradius * U;
   part->gp.geo.angleStart = (double)m_angleStart;
   part->gp.geo.dirSign = dirSign;
   part->Configure();
   part->gp.Reset(dirSign * (double)(m_angleCur - m_angleStart), m_solState);
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
   PLOGI << "Flipper " << m_pflipper->GetName() << " geometry for tests: pivot " << m_hitcircleBase.center.x << " " << m_hitcircleBase.center.y << " start " << RADTOANG(m_angleStart) << " end " << RADTOANG(m_angleEnd) << " L " << m_flipperradius << " rb " << m_hitcircleBase.radius << " re " << m_endradius << " dir " << dirSign;
   host->m_enabled = m_enabled;
   PLOGI << "Flipper " << m_pflipper->GetName() << ": unified part (bat inertia " << part->bat.der.iFlipper << " kg m^2, reach " << part->TipReach() * 1e3 << " mm, rubber e_mpa "
         << part->bat.set.rubberEMpa << ", largest acceleration " << part->alphaMax << " rad/s^2)";
   return host;
}
