// license:GPLv3+
// The table side of a gate's or spinner's moving part (hingepartadapter.h).

#include "core/stdafx.h"

#include "parts/hingepartadapter.h"
#include "parts/gate.h"
#include "parts/spinner.h"
#include "physics/contacthingepart.h"
#include "physics/contactsystem.h"
#include "physics/hitrubber.h"

HingePartAdapter::HingePartAdapter(std::shared_ptr<HingePart::Part> model, Gate *gate, Spinner *spinner, float *angle, float *angleSpeed, const std::string &name)
   : m_model(std::move(model)), m_gate(gate), m_spinner(spinner), m_angle(angle), m_angleSpeed(angleSpeed), m_name(name)
{
}

void HingePartAdapter::Publish()
{
   // the stock sign of the part's state: the stock renderer turns a two-way gate's mesh by +angle and a one-way gate's by -angle (Gate::Render), so a two-way
   // gate's angle is written with the opposite sign of the body's (which is the same for both: positive swings toward the front)
   const double sign = (m_gate != nullptr && m_gate->m_d.m_twoWay) ? -1.0 : 1.0;
   if (m_angle != nullptr)
      *m_angle = (float)(sign * m_model->body.theta);
   if (m_angleSpeed != nullptr)
      *m_angleSpeed = (float)(sign * m_model->body.omega * DEFAULT_STEPTIME_S);   // rad per VPT
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
      }
   }
}

ContactPartHost *HingePartAdapter::MakeHost(IEditable *part, const char *kind, const std::string &name, const std::filesystem::path &ini, std::shared_ptr<HingePart::Part> model, Gate *gate,
                                            Spinner *spinner, float *angle, float *angleSpeed, const double elasticity, const double friction)
{
   constexpr double U = ContactSystem::kMetresPerUnit;
   RubberObjectSettings s;
   ContactHost::LoadHardSurfaceSettings(s, ini, kind, name, elasticity, friction);
   // the box the part can sweep (its reach about the axis), VPX units: the proximity engagement and the stock tree see it
   const RubberContact::V3 &A = model->body.A;
   const double r = model->Reach();
   const float lo[3] = { (float)((A.x - r) / U), (float)((A.y - r) / U), (float)((A.z - r) / U) };
   const float hi[3] = { (float)((A.x + r) / U), (float)((A.y + r) / U), (float)((A.z + r) / U) };
   ContactPartHost *const host = new ContactPartHost(part, model, std::make_unique<HingePartAdapter>(model, gate, spinner, angle, angleSpeed, name), s, ini, name, kind, lo, hi);
   host->SetHardSurface(elasticity, friction);
   // (the gravity it feels is set at every advance: the physics engine is still being built here)
   PLOGI << (gate != nullptr ? "Gate " : "Spinner ") << name << ": unified part (mass " << model->body.mp.m * 1e3 << " g, inertia " << model->body.mp.I << " kg m^2, centre of mass "
         << model->body.PendulumArm() * 1e3 << " mm from the axis, reach " << r * 1e3 << " mm)";
   return host;
}
