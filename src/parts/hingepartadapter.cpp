// license:GPLv3+
// The table side of a gate's or spinner's moving part (hingepartadapter.h).

#include "core/stdafx.h"

#include "parts/hingepartadapter.h"
#include "parts/gate.h"
#include "parts/spinner.h"
#include "physics/contacthingepart.h"
#include "physics/contactsettings.h"
#include "physics/contactsystem.h"
#include "physics/hitrubber.h"
#include "ui/live/LiveUI.h"

HingePartAdapter::HingePartAdapter(std::shared_ptr<HingeModel::Engine> model, Gate *gate, Spinner *spinner, float *angle, float *angleSpeed, const std::string &name)
   : m_model(std::move(model)), m_gate(gate), m_spinner(spinner), m_angle(angle), m_angleSpeed(angleSpeed), m_name(name)
{
   const ContactSettings &cs = ContactSettings::Get();
   m_model->gp.dynRatio = cs.muDynamicRatio;
   m_model->gp.SetPieceSum(cs.hingePieceSum);
   if (cs.hingePieceSum && !m_model->gp.OverlapWarn().empty())
      PLOGW << "hinge_piece_sum: part " << m_name << ": " << m_model->gp.OverlapWarn();
}

void HingePartAdapter::Publish()
{
   // the stock sign of the part's state: the stock renderer turns a two-way gate's mesh by +angle and a one-way gate's by -angle (Gate::Render), so a two-way
   // gate's angle is written with the opposite sign of the body's (which is the same for both: positive swings toward the front)
   m_model->Collect();
   m_model->bookMs = (int)g_pplayer->m_time_msec;
   const double sign = (m_gate != nullptr && m_gate->m_d.m_twoWay) ? -1.0 : 1.0;
   if (m_angle != nullptr)
      *m_angle = (float)(sign * m_model->gp.q);
   if (m_angleSpeed != nullptr)
      *m_angleSpeed = (float)(sign * m_model->gp.rate * DEFAULT_STEPTIME_S);   // rad per VPT
   // detector (rule 8): a step of the part alone was taken on energy, or could not be solved
   if (m_model->energySteps > m_esSeen)
   {
      m_esSeen = m_model->energySteps;
      if (m_esLines < 5)
      {
         m_esLines++;
         PLOGW << "Group step: energy step alone, part " << m_name << ", balls 0, t=" << g_pplayer->m_time_msec << " ms, " << m_esSeen << " so far; " << m_model->StateNote() << ". A net, not expected to fire: keep vpinball.log and report it";
      }
   }
   if (m_model->failedSteps > 0 && !m_failShown)
   {
      m_failShown = true;
      const std::string msg = m_name + ": a step of the part alone could not be solved: the rest of that cycle was lost";
      PLOGW << msg;
      if (g_pplayer->m_liveUI != nullptr)
         g_pplayer->m_liveUI->PushNotification(msg, 6000);
   }
}

void HingePartAdapter::Fire(const std::vector<ContactMovingPart::Event> &events)
{
   for (const ContactMovingPart::Event &e : events)
   {
      if (m_gate == nullptr && m_spinner == nullptr)
         break;
      switch (e.kind)
      {
      case ContactMovingPart::Event::Hit:
         if (m_gate != nullptr)
            m_gate->FireGroupEvent(DISPID_HitEvents_Hit);
         break;
      case ContactMovingPart::Event::Spin:
         if (m_spinner != nullptr)
            m_spinner->FireGroupEvent(DISPID_SpinnerEvents_Spin);
         break;
      case ContactMovingPart::Event::LimitEOS:
         if (m_gate != nullptr)
            m_gate->FireVoidEventParm(DISPID_LimitEvents_EOS, (float)e.value);
         else
            m_spinner->FireVoidEventParm(DISPID_LimitEvents_EOS, (float)e.value);
         break;
      case ContactMovingPart::Event::LimitBOS:
         if (m_gate != nullptr)
            m_gate->FireVoidEventParm(DISPID_LimitEvents_BOS, (float)e.value);
         else
            m_spinner->FireVoidEventParm(DISPID_LimitEvents_BOS, (float)e.value);
         break;
      case ContactMovingPart::Event::SwitchClosed:
      case ContactMovingPart::Event::SwitchOpened:
         if (ContactTraceOn())   // (no script event yet: SwitchClosed / SwitchOpened for gates and spinners come with the IDL change, design build step 7)
            PLOGI << "PartSwitch " << m_name << (e.kind == ContactMovingPart::Event::SwitchClosed ? " closed" : " opened") << " t=" << g_pplayer->m_time_msec;
         break;
      default:
         break;
      }
   }
}

ContactPartHost *HingePartAdapter::MakeHost(IEditable *part, const char *kind, const std::string &name, const std::filesystem::path &ini, std::shared_ptr<HingeModel::Engine> model, Gate *gate,
                                            Spinner *spinner, float *angle, float *angleSpeed, const double elasticity, const double friction)
{
   constexpr double U = ContactSystem::kMetresPerUnit;
   RubberObjectSettings s;
   ContactHost::LoadHardSurfaceSettings(s, ini, kind, name, elasticity, friction);
   // the box the part can sweep (its reach about the axis), VPX units: the proximity engagement and the stock tree see it
   const RubberContact::V3 &A = model->gp.A;
   const double r = model->Reach();
   const float lo[3] = { (float)((A.x - r) / U), (float)((A.y - r) / U), (float)((A.z - r) / U) };
   const float hi[3] = { (float)((A.x + r) / U), (float)((A.y + r) / U), (float)((A.z + r) / U) };
   ContactPartHost *const host = new ContactPartHost(part, model, std::make_unique<HingePartAdapter>(model, gate, spinner, angle, angleSpeed, name), s, ini, name, kind, lo, hi);
   host->SetHardSurface(elasticity, friction);
   // (the gravity it feels is set at every advance: the physics engine is still being built here)
   PLOGI << (gate != nullptr ? "Gate " : "Spinner ") << name << ": unified part on the group step (mass " << model->gp.mp.m * 1e3 << " g, inertia " << model->gp.mp.I << " kg m^2, centre of mass "
         << model->gp.PendulumArm() * 1e3 << " mm from the axis, reach " << r * 1e3 << " mm)";
   return host;
}
