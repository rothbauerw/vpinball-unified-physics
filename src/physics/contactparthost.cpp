// license:GPLv3+
// The host of a moving part of the contact system (contactparthost.h).

#include "core/stdafx.h"

#include "parts/pintable.h"
#include "core/player.h"

#include "physics/contactparthost.h"
#include "physics/contactsystem.h"
#include "physics/hitrubber.h"
#include "physics/PhysicsEngine.h"
#include "physics/cabinet/NudgeHandler.h"

ContactPartHost::ContactPartHost(IEditable *editable, std::shared_ptr<ContactMovingPart> part, std::unique_ptr<PartAdapter> adapter, const RubberObjectSettings &settings,
                                 const std::filesystem::path &iniPath, const std::string &name, const char *kind, const float lo[3], const float hi[3])
   : ContactHost(editable, settings, iniPath, name, kind), m_part(std::move(part)), m_adapter(std::move(adapter))
{
   for (int i = 0; i < 3; i++)
   {
      m_lo[i] = lo[i];
      m_hi[i] = hi[i];
   }
}

RubberContact::V3 ContactPartHost::PartGravitySI()
{
   // the balls' ExternalVelocityIncrement (gravity, nudge) without the magnets, as an acceleration: per physics step in VPX units -> m/s^2
   if (g_pplayer == nullptr || g_pplayer->m_physics == nullptr || g_pplayer->m_pininput.m_nudgeHandler == nullptr)
      return RubberContact::V3(0.0, 0.0, -9.81);
   const Vertex3Ds g = g_pplayer->m_physics->GetGravity();
   const float slope = ANGTORAD(g_pplayer->m_ptable->GetPlayfieldSlope());
   const Vertex2D &ca = g_pplayer->m_pininput.m_nudgeHandler->GetCabinetAcceleration();
   Vertex3Ds dv = (float)PHYS_FACTOR * g;
   dv.x -= (float)PHYS_FACTOR * MS2TOVPUVPT2(ca.x);
   dv.y -= (float)PHYS_FACTOR * MS2TOVPUVPT2(ca.y) * cosf(slope);
   dv.z -= (float)PHYS_FACTOR * MS2TOVPUVPT2(ca.y) * sinf(slope);
   const double k = ContactSystem::kVelocityScale / PHYSICS_STEPTIME_S;
   return RubberContact::V3(dv.x * k, dv.y * k, dv.z * k);
}

void ContactPartHost::PartAdvance(const float dtime)
{
   const RubberContact::V3 g = PartGravitySI();
   if (!m_gravityLogged)
   {
      m_gravityLogged = true;
      PLOGI << "Part " << m_name << ": gravity felt " << g.x << ' ' << g.y << ' ' << g.z << " m/s^2 (|g| " << g.Length() << ")";
   }
   m_part->SetGravity(g);   // (a nudge wakes a sleeping part)
   {
      // the contact system stepped it with its balls for steppedS; the rest of the interval it moves alone
      const double rest = (double)dtime * DEFAULT_STEPTIME_S - m_part->steppedS;
      if (rest > 1e-9)
         m_part->AdvanceAlone(rest);
      m_part->steppedS = 0.0;
      m_part->steppedThisCycle = false;
   }
   m_adapter->Publish();
   // the part's events, in the order they happened: queued, the engine fires them after the free balls have moved (ContactSystem::Dispatch)
   std::vector<ContactMovingPart::Event> q;
   q.swap(m_part->events);
   if (q.empty())
      return;
   PartAdapter *const adapter = m_adapter.get();
   ContactSystem::DeferEvent([adapter, q]() { adapter->Fire(q); });
}

void ContactPartHost::Reload() const
{
   if (m_adapter->FixedSettings())
      return;   // (a turntable disc keeps the settings it was built with: its keys are read when the table loads)
   RubberObjectSettings s;
   if (m_hardSurface)
      LoadHardSurfaceSettings(s, m_ini, m_surfaceKind.c_str(), m_name, m_partE, m_partMu);
   else
      LoadSurfaceSettings(s, m_ini, m_surfaceKind.c_str(), m_name);
   m_set = s;
   PLOGI << "Rubber surface settings for " << m_name << " reloaded from " << m_ini.string() << ": rubber " << (s.contact ? 1 : 0) << ", alpha " << s.alphaSPerM << ", mu " << s.mu;
}

void ContactPartHost::Refresh(SurfaceContact &it, double, double, const RubberContact::V3 &) const
{
   // the part probes its own few wires and faces: no candidate list
   it.candR = -1.0;
   it.cand.clear();
   it.candV.clear();
}

RubberContact::TubeProbe ContactPartHost::Probe(const RubberContact::BallSI &bs, const SurfaceContact &, const RubberContact::V3 &c) const
{
   m_part->SetBall(bs);
   return m_part->Probe(c);
}

void ContactPartHost::VelocityAt(const RubberContact::TubeProbe &pr, RubberContact::V3 &v, RubberContact::V3 &w) const
{
   m_part->SurfaceVelocity(pr.point, v, w);
}

RubberContact::ContactResult ContactPartHost::EvaluateOwn(SurfaceContact &it, const RubberContact::BallSI &bs, const RubberContact::V3 &c, const RubberContact::V3 &v, const RubberContact::V3 &w,
                                                          const double h) const
{
   m_part->SetBall(bs);
   return m_part->ContactItem(0, it, c, v, w, h);
}

double ContactPartHost::SurfaceSpeedBoundMS(const double cycleS) const
{
   return m_part->SurfaceSpeedBound(cycleS);
}

void ContactPartHost::CalcHitBBox()
{
   const float m = 1.0f;
   m_hitBBox = FRect3D(m_lo[0] - m, m_hi[0] + m, m_lo[1] - m, m_hi[1] + m, m_lo[2] - m, m_hi[2] + m);
   if (ContactTraceOn() && m_set.contact)
      PLOGI << "RubberGeom " << m_name << " " << m_surfaceKind << " x " << m_lo[0] << " " << m_hi[0] << " y " << m_lo[1] << " " << m_hi[1] << " z " << m_lo[2] << " " << m_hi[2];
}
