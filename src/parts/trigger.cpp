// license:GPLv3+

#include "core/stdafx.h"
#include "trigger.h"

#include "core/VPApp.h"
#include "math/MeshUtils.h"
#include "math/matrix.h"
#include "meshes/triggerButtonMesh.h"
#include "meshes/triggerInderMesh.h"
#include "meshes/triggerSimpleMesh.h"
#include "meshes/triggerStarMesh.h"
#include "meshes/triggerWireDMesh.h"
#include "parts/ball.h"
#include "parts/Collection.h"
#include "physics/collideex.h"
#include "physics/contacthost.h"
#include "physics/contactparthost.h"
#include "physics/contactpartsini.h"
#include "physics/contactsettings.h"
#include "physics/contactsystem.h"
#include "physics/contactwireswitch.h"
#include "physics/hitrubber.h"
#include "physics/physicsini.h"
#include "physics/testdriver.h"
#include "core/player.h"
#include "ui/live/LiveUI.h"
#include "physics/PhysicsEngine.h"
#include "renderer/Renderer.h"
#include "renderer/Shader.h"
#include "renderer/trace.h"
#include "utils/objloader.h"


Trigger::~Trigger()
{
   assert(m_renderer == nullptr);
}

Trigger *Trigger::CopyForPlay() const
{
   STANDARD_EDITABLE_WITH_DRAGPOINT_COPY_FOR_PLAY_IMPL(Trigger, m_curve)
   return dst;
}

Trigger::StaticMeshData Trigger::SetupMeshData(const TriggerShape shape)
{
   switch (shape)
   {
   case TriggerNone: return {};
   case TriggerWireA:
   case TriggerWireB:
   case TriggerWireC: return { triggerSimple, triggerSimpleIndices };
   case TriggerWireD: return { triggerDWireMesh, triggerDWireIndices };
   case TriggerInder: return { triggerInderMesh, triggerInderIndices };
   case TriggerButton: return { triggerButtonMesh, triggerButtonIndices };
   case TriggerStar: return { triggerStar, triggerStarIndices };
   default: return {};
   }
}

namespace
{
   WireShape::Shape ToWireShape(const TriggerShape shape)
   {
      return shape == TriggerWireB ? WireShape::Shape::WireB : shape == TriggerWireC ? WireShape::Shape::WireC : shape == TriggerWireD ? WireShape::Shape::WireD
         : shape == TriggerWireE ? WireShape::Shape::WireE : WireShape::Shape::WireA;
   }
}

void Trigger::GetWireOutline(vector<Vertex2D> &outline) const
{
   if (CentrelineShape(m_d.m_shape))
   {
      // the plan view of the wire's tube (the same centreline as the physics, wireshape.h): its plane is vertical, so a strip one wire wide from end to end
      const WireSwitchGroup::Settings ds;   // the editor's outline: the default steps and wire (the table's ini is read at play)
      const WireShape::Centreline cl = WireShape::MakeCentreline(ToWireShape(m_d.m_shape), ds.bendStepDeg, ds.archStepDeg);
      WireShape::Placement pl;
      pl.x = m_d.m_vCenter.x;
      pl.y = m_d.m_vCenter.y;
      pl.rotDeg = m_d.m_rotation;
      const WireShape::WireBody wb = WireShape::Place(cl, pl, WireSwitchGroup::Settings().wireDMm);
      double s0 = 1e30, s1 = -1e30;
      for (const WireShape::V3 &p : wb.nodes)
      {
         const double s = (p - wb.pivot).Dot(wb.lane);
         s0 = std::min(s0, s - wb.r);
         s1 = std::max(s1, s + wb.r);
      }
      const WireShape::V3 e(wb.lane.y, -wb.lane.x, 0.);
      const auto corner = [&](const double s, const double w) { const WireShape::V3 p = wb.pivot + wb.lane * s + e * w; return Vertex2D((float)p.x, (float)p.y); };
      outline = { corner(s0, -wb.r), corner(s1, -wb.r), corner(s1, wb.r), corner(s0, wb.r), corner(s0, -wb.r) };
      return;
   }
   const StaticMeshData meshData = SetupMeshData(m_d.m_shape);
   if (meshData.indices.empty())
      return;

   const float scaleX = (m_d.m_shape == TriggerStar || m_d.m_shape == TriggerButton) ? m_d.m_radius : m_d.m_scaleX;
   const float scaleY = (m_d.m_shape == TriggerStar || m_d.m_shape == TriggerButton) ? m_d.m_radius : m_d.m_scaleY;
   const Matrix3D fullMatrix
      = Matrix3D::MatrixRotateZ(ANGTORAD(m_d.m_rotation)) * Matrix3D::MatrixScale(scaleX, scaleY, 1.0f) * Matrix3D::MatrixTranslate(m_d.m_vCenter.x, m_d.m_vCenter.y, 0.0f);

   outline.reserve(meshData.indices.size() / 3 + 1);
   const Vertex3Ds A = fullMatrix * meshData.vertices[meshData.indices[0]];
   outline.emplace_back(A.x, A.y);
   for (size_t i = 0; i < meshData.indices.size(); i += 3)
   {
      const Vertex3Ds B = fullMatrix * meshData.vertices[meshData.indices[i + 1]];
      outline.emplace_back(B.x, B.y);
   }
}

void Trigger::InitShape(float x, float y)
{
   constexpr float lengthX = 30.0f;
   constexpr float lengthY = 30.0f;

   m_curve.ClearPoints();

   // First time shape has been set to custom - set up some points
   m_curve.PushPoint(std::make_unique<DragPoint>(&m_curve, x - lengthX, y - lengthY, 0.f, false));
   m_curve.PushPoint(std::make_unique<DragPoint>(&m_curve, x - lengthX, y + lengthY, 0.f, false));
   m_curve.PushPoint(std::make_unique<DragPoint>(&m_curve, x + lengthX, y + lengthY, 0.f, false));
   m_curve.PushPoint(std::make_unique<DragPoint>(&m_curve, x + lengthX, y - lengthY, 0.f, false));
}

HRESULT Trigger::Init(const float x, const float y, const bool fromMouseClick, const bool forPlay)
{
   SetDefaults(fromMouseClick);
   m_d.m_vCenter.x = x;
   m_d.m_vCenter.y = y;
   if (m_curve.GetPoints().empty())
      InitShape(x, y);
   return S_OK;
}

void Trigger::SetDefaults(const bool fromMouseClick)
{
#define LinkProp(field, prop) field = fromMouseClick ? g_settingsService.GetAppSettings().GetDefaultPropsTrigger_##prop() : Settings::GetDefaultPropsTrigger_##prop##_Default()
   LinkProp(m_d.m_radius, Radius);
   LinkProp(m_d.m_rotation, Rotation);
   LinkProp(m_d.m_wireThickness, WireThickness);
   LinkProp(m_d.m_scaleX, ScaleX);
   LinkProp(m_d.m_scaleY, ScaleY);
   LinkProp(m_d.m_enabled, Enabled);
   LinkProp(m_d.m_visible, Visible);
   LinkProp(m_d.m_hit_height, HitHeight);
   LinkProp(m_d.m_shape, Shape);
   LinkProp(m_d.m_szSurface, Surface);
   LinkProp(m_d.m_animSpeed, AnimSpeed);
   LinkProp(m_d.m_reflectionEnabled, ReflectionEnabled);
   LinkProp(m_timerEnabled, TimerEnabled);
   LinkProp(m_timerInterval, TimerInterval);
#undef LinkProp
}

void Trigger::WriteRegDefaults()
{
#define LinkProp(field, prop) g_settingsService.GetAppSettings().SetDefaultPropsTrigger_##prop(field, false)
   LinkProp(m_d.m_radius, Radius);
   LinkProp(m_d.m_rotation, Rotation);
   LinkProp(m_d.m_wireThickness, WireThickness);
   LinkProp(m_d.m_scaleX, ScaleX);
   LinkProp(m_d.m_scaleY, ScaleY);
   LinkProp(m_d.m_enabled, Enabled);
   LinkProp(m_d.m_visible, Visible);
   LinkProp(m_d.m_hit_height, HitHeight);
   LinkProp(m_d.m_shape, Shape);
   LinkProp(m_d.m_szSurface, Surface);
   LinkProp(m_d.m_animSpeed, AnimSpeed);
   LinkProp(m_d.m_reflectionEnabled, ReflectionEnabled);
   LinkProp(m_timerEnabled, TimerEnabled);
   LinkProp(m_timerInterval, TimerInterval);
#undef LinkProp
}


#pragma region Physics

// Ported at: VisualPinball.Engine/VPT/Trigger/TriggerHitGenerator.cs

namespace
{
   // The table side of a wire trigger's rollover wire switch (contactwireswitch.h): the switch's close is the trigger's Hit, its open the Unhit; ActiveBall is the ball
   // nearest the wire (one ball presses a wire at a time, design-wire.md topic 6). It publishes the part's angle, which turns the drawn switch (Trigger::Render).
   class WireSwitchAdapter final : public PartAdapter
   {
   public:
      WireSwitchAdapter(Trigger *trigger, std::shared_ptr<WireSwitchModel::WireSwitch> model, std::vector<WireShape::Capsule> chainVpu)
         : m_trigger(trigger), m_model(std::move(model)), m_chain(std::move(chainVpu)), m_name(trigger->GetName())
      {
         const ContactSettings &cs = ContactSettings::Get();
         m_model->gp.dynRatio = cs.muDynamicRatio;
         m_model->gp.SetPieceSum(cs.hingePieceSum);
      }
      void Publish() override
      {
         m_model->Collect();
         m_model->bookMs = (int)g_pplayer->m_time_msec;
         m_trigger->m_partAngle = (float)m_model->part.q;   // (the angle the renderer turns the drawn switch by)
         if (!m_notices.empty() && g_pplayer->m_liveUI != nullptr)   // the load-time detectors, on screen once the player runs
         {
            for (const std::string &msg : m_notices)
               g_pplayer->m_liveUI->PushNotification(msg, 10000);
            m_notices.clear();
         }
         static const bool gapLog = []() { const char *const e = getenv("VPX_WIRE_GAP"); return e != nullptr && atoi(e) != 0; }();
         if (gapLog && m_host != nullptr && m_gapLines < 600)   // diagnostic (off by default): each ball near the wire, its probe gap and whether its contact holds the wire
         {
            constexpr double U = ContactSystem::kMetresPerUnit;
            if (m_gapLines == 0)
            {
               const FRect3D &bb = m_host->m_hitBBox;
               PLOGI << "WireGap " << m_name << " box x " << bb.left << " " << bb.right << " y " << bb.top << " " << bb.bottom << " z " << bb.zlow << " " << bb.zhigh << " enabled " << m_host->m_enabled << " contact " << m_host->ContactSettings().contact;
               m_gapLines++;
            }
            for (Ball *b : g_pplayer->m_vball)
            {
               const HitBall &hb = b->m_hitBall;
               const Vertex3Ds &p = hb.m_d.m_pos;
               if (WireShape::DistToChain(WireShape::V3(p.x, p.y, p.z), m_chain) > hb.m_d.m_radius + 40.0)
                  continue;
               RubberContact::BallSI bs;
               bs.R = hb.m_d.m_radius * U;
               m_model->SetBall(bs);
               const RubberContact::TubeProbe pr = m_model->Probe(RubberContact::V3(p.x * U, p.y * U, p.z * U));
               int held = -1;
               if (hb.m_contact != nullptr)
                  held = std::any_of(hb.m_contact->items.begin(), hb.m_contact->items.end(), [this](const SurfaceContact &it) { return it.r == m_host; }) ? 1 : 0;
               PLOGI << "WireGap " << m_name << " t=" << g_pplayer->m_time_msec << " ball (" << p.x << ", " << p.y << ", " << p.z << ") probe gap " << -pr.delta * 1e3 << " mm, chain gap "
                     << (WireShape::DistToChain(WireShape::V3(p.x, p.y, p.z), m_chain) - hb.m_d.m_radius - m_chain[0].r) * U * 1e3 << " mm, in a contact " << (hb.m_contact != nullptr) << ", holds the wire " << held
                     << ", items " << (hb.m_contact != nullptr ? (int)hb.m_contact->items.size() : 0) << ", q " << m_model->gp.q;
               m_gapLines++;
            }
         }
         if (m_model->energySteps > m_esSeen)   // detector (rule 8)
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
            const std::string msg = m_name + ": a step of the wire alone could not be solved: the rest of that cycle was lost";
            PLOGW << msg;
            if (g_pplayer->m_liveUI != nullptr)
               g_pplayer->m_liveUI->PushNotification(msg, 6000);
         }
      }
      void Fire(const std::vector<ContactMovingPart::Event> &events) override
      {
         for (const ContactMovingPart::Event &e : events)
         {
            if (e.kind != ContactMovingPart::Event::SwitchClosed && e.kind != ContactMovingPart::Event::SwitchOpened)
               continue;
            const bool closed = e.kind == ContactMovingPart::Event::SwitchClosed;
            if (closed)
            {
               if (!m_trigger->m_d.m_enabled)
                  continue;
               HitBall *best = nullptr;
               double bestD = INFINITY;
               for (Ball *b : g_pplayer->m_vball)
               {
                  const Vertex3Ds &p = b->m_hitBall.m_d.m_pos;
                  const double d = WireShape::DistToChain(WireShape::V3(p.x, p.y, p.z), m_chain) - b->m_hitBall.m_d.m_radius;
                  if (d < bestD)
                  {
                     bestD = d;
                     best = &b->m_hitBall;
                  }
               }
               if (best == nullptr)
                  continue;
               m_holder = best;
               g_pplayer->m_pactiveball = best->m_pBall;
               if (best->m_d.m_vpVolObjs != nullptr)
                  best->m_d.m_vpVolObjs->push_back(m_trigger);
               m_trigger->TriggerAnimationHit();
            }
            else
            {
               if (m_holder == nullptr)
                  continue;   // (closed while disabled)
               const auto &vb = g_pplayer->m_vball;
               if (std::any_of(vb.begin(), vb.end(), [&](const Ball *b) { return &b->m_hitBall == m_holder; }))
               {
                  g_pplayer->m_pactiveball = m_holder->m_pBall;
                  if (m_holder->m_d.m_vpVolObjs != nullptr)
                     std::erase(*m_holder->m_d.m_vpVolObjs, (IFireEvents *)m_trigger);
               }
               m_holder = nullptr;
               m_trigger->TriggerAnimationUnhit();
            }
            (closed ? m_hits : m_unhits)++;
            if (m_hits + m_unhits <= 40)
               PLOGI << "Trigger " << m_name << ": wire switch " << (closed ? "closed: Hit" : "opened: Unhit") << " t=" << g_pplayer->m_time_msec << " ms (" << m_hits << " Hit, " << m_unhits << " Unhit so far)";
            m_trigger->OnHit(closed);
         }
      }

   private:
      Trigger *m_trigger;
      std::shared_ptr<WireSwitchModel::WireSwitch> m_model;
      std::vector<WireShape::Capsule> m_chain;   // the wire at rest, VPX units: ActiveBall is the ball nearest it
      std::string m_name;
      HitBall *m_holder = nullptr;   // the ball that closed the switch
      long m_esSeen = 0;
      int m_esLines = 0, m_hits = 0, m_unhits = 0;
      bool m_failShown = false;
      int m_gapLines = 0;

   public:
      const ContactHost *m_host = nullptr;   // (the gap log)
      std::vector<std::string> m_notices;    // the load-time detectors' messages, shown at the first cycle
   };

   // The table side of a Button trigger's cap (buttoncap.h, runs/2026-11-10-rollover/design.md): the switch's close is the trigger's Hit, its open the Unhit, as the wire switch;
   // ActiveBall is the ball nearest the cap's axis. A disabled trigger: no events, its cap still moves (design 4e). The ball that closed it destroyed: no script Unhit, the
   // state link opens (the None rule, design 4d). It publishes the drawn cap's height offset (VPU) in m_partAngle: the stock mesh follows the cap (UpdateAnimation)
   class ButtonCapAdapter final : public PartAdapter
   {
   public:
      ButtonCapAdapter(Trigger *trigger, std::shared_ptr<ButtonCap::Cap> model, const double drawRestVpu)
         : m_trigger(trigger), m_model(std::move(model)), m_name(trigger->GetName()), m_drawRest(drawRestVpu)
      {
      }
      void Publish() override
      {
         m_model->Collect();
         m_model->bookMs = (int)g_pplayer->m_time_msec;
         m_trigger->m_partAngle = (float)(m_drawRest - m_model->part.Sink() / ContactSystem::kMetresPerUnit);   // (the axial sink; the tilt is not drawn)
         if (!m_notices.empty() && g_pplayer->m_liveUI != nullptr)   // the load-time detectors, on screen once the player runs
         {
            for (const std::string &msg : m_notices)
               g_pplayer->m_liveUI->PushNotification(msg, 10000);
            m_notices.clear();
         }
         if (m_model->energySteps > m_esSeen)   // detector (rule 8)
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
            const std::string msg = m_name + ": a step of the button cap alone could not be solved: the rest of that cycle was lost";
            PLOGW << msg;
            if (g_pplayer->m_liveUI != nullptr)
               g_pplayer->m_liveUI->PushNotification(msg, 6000);
         }
      }
      void Fire(const std::vector<ContactMovingPart::Event> &events) override
      {
         for (const ContactMovingPart::Event &e : events)
         {
            if (e.kind != ContactMovingPart::Event::SwitchClosed && e.kind != ContactMovingPart::Event::SwitchOpened)
               continue;
            const bool closed = e.kind == ContactMovingPart::Event::SwitchClosed;
            if (closed)
            {
               if (!m_trigger->m_d.m_enabled)
                  continue;
               HitBall *best = nullptr;
               double bestD = INFINITY;
               for (Ball *b : g_pplayer->m_vball)
               {
                  const Vertex3Ds &p = b->m_hitBall.m_d.m_pos;
                  const double dx = p.x - m_trigger->m_d.m_vCenter.x, dy = p.y - m_trigger->m_d.m_vCenter.y, d = dx * dx + dy * dy;
                  if (d < bestD)
                  {
                     bestD = d;
                     best = &b->m_hitBall;
                  }
               }
               if (best == nullptr)
                  continue;
               m_holder = best;
               g_pplayer->m_pactiveball = best->m_pBall;
               if (best->m_d.m_vpVolObjs != nullptr)
                  best->m_d.m_vpVolObjs->push_back(m_trigger);
               m_trigger->TriggerAnimationHit();
            }
            else
            {
               if (m_holder == nullptr)
                  continue;   // (closed while disabled)
               const auto &vb = g_pplayer->m_vball;
               HitBall *const holder = m_holder;
               m_holder = nullptr;
               m_trigger->TriggerAnimationUnhit();
               if (!std::any_of(vb.begin(), vb.end(), [&](const Ball *b) { return &b->m_hitBall == holder; }))
               {
                  if (!m_trigger->m_d.m_stateSrcLink.empty())
                     g_pplayer->m_resURIResolver.SetBoolState(m_trigger->m_d.m_stateSrcLink, false);
                  PLOGI << "Trigger " << m_name << ": button opened after its ball was destroyed: no Unhit, state link off, t=" << g_pplayer->m_time_msec << " ms";
                  continue;
               }
               g_pplayer->m_pactiveball = holder->m_pBall;
               if (holder->m_d.m_vpVolObjs != nullptr)
                  std::erase(*holder->m_d.m_vpVolObjs, (IFireEvents *)m_trigger);
            }
            (closed ? m_hits : m_unhits)++;
            if (m_hits + m_unhits <= 40)
               PLOGI << "Trigger " << m_name << ": button switch " << (closed ? "closed: Hit" : "opened: Unhit") << " t=" << g_pplayer->m_time_msec << " ms (" << m_hits << " Hit, " << m_unhits << " Unhit so far)";
            m_trigger->OnHit(closed);
         }
      }

   private:
      Trigger *m_trigger;
      std::shared_ptr<ButtonCap::Cap> m_model;
      std::string m_name;
      double m_drawRest;             // the drawn mesh's height offset at rest (VPU): its top at the cap's physical top (looks only)
      HitBall *m_holder = nullptr;   // the ball that closed the switch
      long m_esSeen = 0;
      int m_esLines = 0, m_hits = 0, m_unhits = 0;
      bool m_failShown = false;
   public:
      std::vector<std::string> m_notices;    // the load-time detectors' messages, shown at the first cycle
   };
}

// None trigger on the unified engine (to-do 11.2, runs/2026-11-10-none-trigger/design.md): an invisible sensing zone (an opto beam, a proximity sensor, a script zone). It only reads the ball:
// no force, no position change. Zone: the ball centre inside the polygon (even-odd, half-open on the faces) and height <= z <= height + HitHeight (read live). The ball's volume list
// (m_vpVolObjs) is the state, as stock: BallCntOver, DestroyBall and the kickers read it. The zone list is the engine's (PhysicsEngine::m_noneZones, contactsystem.h).
namespace
{
   const bool g_noneTimed = getenv("VPX_CONTACT_TIMING") != nullptr && *getenv("VPX_CONTACT_TIMING") != 0;   // VPX_CONTACT_TIMING: the zones' cost, logged when the table ends

   int64_t NoneNs() { return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

   // a ball other than this one still listed in the trigger's zone
   bool OtherBallListed(const Trigger *const trigger, const HitBall *const except)
   {
      for (const Ball *const ball : g_pplayer->m_vball)
      {
         const HitBall *const b = &ball->m_hitBall;
         if (b != except && b->m_d.m_vpVolObjs != nullptr && FindIndexOf(*b->m_d.m_vpVolObjs, (IFireEvents *)trigger) >= 0)
            return true;
      }
      return false;
   }
}

bool ContactSystem::NoneZones::Inside(const Zone &z, const double x, const double y, const double zc)
{
   if (zc < z.zlow || zc > (double)z.zlow + z.trigger->m_d.m_hit_height || x < z.x0 || x > z.x1 || y < z.y0 || y > z.y1)
      return false;
   bool in = false;   // even-odd crossing count, half-open: a point on a face belongs to one side only
   for (size_t i = 0, j = z.poly.size() - 1; i < z.poly.size(); j = i++)
   {
      const Vertex2D &a = z.poly[i], &b = z.poly[j];
      if ((a.y > y) != (b.y > y) && x < (double)(b.x - a.x) * (y - a.y) / (double)(b.y - a.y) + a.x)
         in = !in;
   }
   return in;
}

// the ball's list against where its centre now is: Hit for a centre inside and not listed, Unhit for one outside and listed
void ContactSystem::NoneZones::Judge(const Zone &z, HitBall *const b, const bool in, const double t)
{
   vector<IFireEvents *> &vol = *b->m_d.m_vpVolObjs;
   const int i = FindIndexOf(vol, (IFireEvents *)z.trigger);
   if (in && i < 0)
   {
      vol.push_back((IFireEvents *)z.trigger);
      m_events.push_back({ t, z.trigger, b, true });
   }
   else if (!in && i >= 0)
   {
      vol.erase(vol.begin() + i);
      m_events.push_back({ t, z.trigger, b, false });
   }
}

// the placement rule (design 3, 4a, 4d, 4e): a ball created or moved by the script or the driver, and every ball for a zone enabled again or given a new HitHeight, judged where it is
void ContactSystem::NoneZones::ApplyPlacements(const double t)
{
   bool any = !m_placed.empty();
   for (const Zone &z : m_zones)
      any |= z.recheck;
   if (!any)
      return;
   for (Ball *const ball : g_pplayer->m_vball)
   {
      HitBall *const b = &ball->m_hitBall;
      if (b->m_d.m_vpVolObjs == nullptr)
         continue;
      const bool placed = std::find(m_placed.begin(), m_placed.end(), b) != m_placed.end();
      for (const Zone &z : m_zones)
         if (z.trigger->m_d.m_enabled && (placed || z.recheck))
            Judge(z, b, Inside(z, b->m_d.m_pos.x, b->m_d.m_pos.y, b->m_d.m_pos.z), t);
   }
   m_placed.clear();
   for (Zone &z : m_zones)
      z.recheck = false;
}

// the chord s -> e against one zone: the inside state along it changes only where the chord crosses a face or a z limit; each change is judged at its time
void ContactSystem::NoneZones::Chord(const Zone &z, HitBall *const b, const Vertex3Ds &s, const Vertex3Ds &e, const double startS, const double hS)
{
   const double zhigh = (double)z.zlow + z.trigger->m_d.m_hit_height;
   if (max(s.x, e.x) < z.x0 || min(s.x, e.x) > z.x1 || max(s.y, e.y) < z.y0 || min(s.y, e.y) > z.y1 || max(s.z, e.z) < z.zlow || min(s.z, e.z) > zhigh)
      return;   // the chord never meets the zone's box: outside all along, nothing changes
   const double dx = (double)e.x - s.x, dy = (double)e.y - s.y, dz = (double)e.z - s.z;
   vector<double> &cut = m_cut;
   cut.clear();
   for (size_t i = 0, j = z.poly.size() - 1; i < z.poly.size(); j = i++)
   {
      const double ax = z.poly[j].x, ay = z.poly[j].y, rx = z.poly[i].x - ax, ry = z.poly[i].y - ay;
      const double den = dx * ry - dy * rx;
      if (den == 0.0)
         continue;   // parallel to the face: no crossing point (the z cuts and the end states cover the rest)
      const double qx = ax - s.x, qy = ay - s.y;
      const double t = (qx * ry - qy * rx) / den, u = (qx * dy - qy * dx) / den;
      if (t > 0.0 && t < 1.0 && u >= 0.0 && u <= 1.0)
         cut.push_back(t);
   }
   if (dz != 0.0)
   {
      const double t0 = ((double)z.zlow - s.z) / dz, t1 = (zhigh - s.z) / dz;
      if (t0 > 0.0 && t0 < 1.0)
         cut.push_back(t0);
      if (t1 > 0.0 && t1 < 1.0)
         cut.push_back(t1);
   }
   std::sort(cut.begin(), cut.end());
   const int n = (int)cut.size();
   // the state at the start, on each piece between cuts (its middle), at the end; a change between neighbours is judged at the cut between them
   bool was = Inside(z, s.x, s.y, s.z);
   double tPrev = 0.0;
   for (int k = 0; k <= n; k++)
   {
      const double tNext = k < n ? cut[k] : 1.0, tm = 0.5 * (tPrev + tNext);
      const bool in = Inside(z, s.x + tm * dx, s.y + tm * dy, s.z + tm * dz);
      if (in != was)
         Judge(z, b, in, startS + tPrev * hS);
      was = in;
      tPrev = tNext;
   }
   const bool in = Inside(z, e.x, e.y, e.z);
   if (in != was)
      Judge(z, b, in, startS + hS);
}

void ContactSystem::NoneZones::SubStart(const double startS)
{
   if (m_zones.empty())
      return;
   const int64_t t0 = g_noneTimed ? NoneNs() : 0;
   ApplyPlacements(startS);
   const auto &vb = g_pplayer->m_vball;
   m_start.resize(vb.size());
   for (size_t i = 0; i < vb.size(); i++)
      m_start[i] = { &vb[i]->m_hitBall, vb[i]->m_hitBall.m_d.m_pos };
   if (g_noneTimed)
      m_ns += NoneNs() - t0;
}

void ContactSystem::NoneZones::SubTest(const double startS, const double hS)
{
   if (m_zones.empty())
      return;
   const int64_t t0 = g_noneTimed ? NoneNs() : 0;
   const auto &vb = g_pplayer->m_vball;
   for (size_t i = 0; i < vb.size() && i < m_start.size(); i++)
   {
      HitBall *const b = &vb[i]->m_hitBall;
      if (m_start[i].first != b || b->m_d.m_vpVolObjs == nullptr)
         continue;   // (no ball is made or destroyed between SubStart and here)
      for (const Zone &z : m_zones)
         if (z.trigger->m_d.m_enabled)
            Chord(z, b, m_start[i].second, b->m_d.m_pos, startS, hS);
   }
   if (g_noneTimed)
      m_ns += NoneNs() - t0;
}

void ContactSystem::NoneZones::CycleEnd(const double cycleS)
{
   if (m_zones.empty())
      return;
   const int64_t t0 = g_noneTimed ? NoneNs() : 0;
   ApplyPlacements(cycleS);   // a ball the script moved after the last test
   m_cycles++;
   m_ballCycles += g_pplayer->m_vball.size();
   if (g_noneTimed && m_cycles % 1000 == 0)   // (also every 1000 cycles: the bench ends without a table teardown)
      PLOGI << "None triggers: zone test cost " << m_ns * 1e-3 / (double)m_cycles << " us per cycle, " << (m_ballCycles > 0 ? m_ns * 1e-3 / (double)m_ballCycles : 0.0) << " us per ball per cycle (" << m_cycles << " cycles, " << m_ballCycles << " ball cycles)";
   if (m_events.empty())
   {
      if (g_noneTimed)
         m_ns += NoneNs() - t0;
      return;
   }
   m_firing.swap(m_events);
   m_events.clear();
   std::stable_sort(m_firing.begin(), m_firing.end(), [](const Event &a, const Event &b) { return a.t < b.t; });   // crossing order; one ball's Hit before its Unhit at the same time stays first
   if (g_noneTimed)
      m_ns += NoneNs() - t0;
   for (size_t i = 0; i < m_firing.size(); i++)   // (a script may destroy a ball: BallGone clears its entries here)
   {
      const Event ev = m_firing[i];
      if (ev.ball == nullptr)
         continue;
      g_pplayer->m_pactiveball = ev.ball->m_pBall;   // the ball for the script's ActiveBall
      if (ev.hit)
         ev.trigger->TriggerAnimationHit();
      else
         ev.trigger->TriggerAnimationUnhit();
      ev.trigger->OnHit(ev.hit);
   }
   m_firing.clear();
}

void ContactSystem::NoneZones::Placed(const HitBall *const ball)
{
   if (!m_zones.empty())
      m_placed.push_back(ball);
}

void ContactSystem::NoneZones::BallGone(const HitBall *const ball)
{
   if (m_zones.empty())
      return;
   std::erase(m_placed, ball);
   for (Event &ev : m_events)
      if (ev.ball == ball)
         ev.ball = nullptr;
   for (Event &ev : m_firing)
      if (ev.ball == ball)
         ev.ball = nullptr;
   if (ball->m_d.m_vpVolObjs != nullptr)
      for (const Zone &z : m_zones)
         if (FindIndexOf(*ball->m_d.m_vpVolObjs, (IFireEvents *)z.trigger) >= 0 && !z.trigger->m_d.m_stateSrcLink.empty() && !OtherBallListed(z.trigger, ball))
            g_pplayer->m_resURIResolver.SetBoolState(z.trigger->m_d.m_stateSrcLink, false);   // no Unhit (design 4b): the switch opens when no other ball is still in the zone, no script call
}

void ContactSystem::NoneZones::Add(Trigger *const trigger, vector<Vertex2D> &&poly, const float zlow)
{
   Zone z { trigger, std::move(poly), FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX, zlow, false };
   for (const Vertex2D &v : z.poly)
   {
      z.x0 = min(z.x0, v.x);
      z.x1 = max(z.x1, v.x);
      z.y0 = min(z.y0, v.y);
      z.y1 = max(z.y1, v.y);
   }
   m_zones.push_back(std::move(z));
}

void ContactSystem::NoneZones::Remove(const Trigger *const trigger)
{
   std::erase_if(m_zones, [trigger](const Zone &z) { return z.trigger == trigger; });
   if (!m_zones.empty())
      return;
   if (g_noneTimed && m_cycles > 0)
      PLOGI << "None triggers: zone test cost " << m_ns * 1e-3 / (double)m_cycles << " us per cycle, " << (m_ballCycles > 0 ? m_ns * 1e-3 / (double)m_ballCycles : 0.0) << " us per ball per cycle (" << m_cycles << " cycles, " << m_ballCycles << " ball cycles)";
   m_start.clear();
   m_placed.clear();
   m_events.clear();
   m_firing.clear();
   m_ns = m_cycles = m_ballCycles = 0;
}

void ContactSystem::NoneZones::Recheck(const Trigger *const trigger)
{
   for (Zone &z : m_zones)
      if (z.trigger == trigger)
         z.recheck = true;
}

void Trigger::PhysicSetup(PhysicsEngine* physics, const bool isUI)
{
   m_hitEvent = false;

   if (!isUI && GetPartGroup() != nullptr && GetPartGroup()->GetReferenceSpace() != PartGroupData::SpaceReference::SR_PLAYFIELD)
      return;

   const float height = m_ptable->GetSurfaceHeight(m_d.m_szSurface, m_d.m_vCenter.x, m_d.m_vCenter.y);
   const bool wire = CentrelineShape(m_d.m_shape);
   if (!isUI && ContactSystem::Unified())
   {
      const std::filesystem::path ini = ContactHost::IniPathFor(m_ptable->m_filename);
      const bool on = PhysicsIni::Active(ini).triggers;
      static const void *s_noted = nullptr;   // once per play
      if (s_noted != g_pplayer)
      {
         s_noted = g_pplayer;
         PLOGI << "Triggers: [Features] triggers = " << (on ? 1 : 0) << (on ? ": WireA-E are rollover switches of the contact system (group step); None is a sensing zone of the unified engine; Button is a rollover cap of the group step; Star and Inder stay stock sensors (not converted)" : ": None is a sensing zone of the unified engine; Button is a rollover cap of the group step; every other trigger is a stock sensor");
      }
      if (on && wire)
      {
         // the real wire (PartData S7): a steel wire on a micro switch lever, the switch on a bracket under the playfield (GroupStepDesign 23.4); no stock sensor, mover or volume
         constexpr double U = ContactSystem::kMetresPerUnit;
         WireSwitchModel::Settings s;
         std::string ignoredKeys;
         PartsIni::LoadWireSwitch(s, ini, GetName(), &ignoredKeys);
         std::vector<std::string> notices;   // load-time detectors (rule 8): on screen at the first cycle and in the log as WARN
         if (!ignoredKeys.empty())
            notices.push_back("Wire trigger " + GetName() + ": ini key(s) " + ignoredKeys + " not used (removed or unknown); remove them from the table's physics ini");
         if (s.wireDMm <= 0.0)
         {
            notices.push_back("Wire trigger " + GetName() + ": wire_d_mm " + std::to_string(s.wireDMm) + " is not a wire; the default 1.0 mm is used. Set wire_d_mm above 0");
            s.wireDMm = WireSwitchGroup::Settings().wireDMm;
         }
         if (s.wireDMm >= WireShape::kMaxWireDMm)
            notices.push_back("Wire trigger " + GetName() + ": wire_d_mm " + std::to_string(s.wireDMm) + " folds on the inside of its 1.3 mm bends (needs under 2.6 mm); check wire_d_mm");
         if (!(s.bendStepDeg > 0.0))
         {
            notices.push_back("Wire trigger " + GetName() + ": bend_step_deg " + std::to_string(s.bendStepDeg) + " cuts no pieces; the default is used. Set bend_step_deg above 0");
            s.bendStepDeg = WireSwitchGroup::Settings().bendStepDeg;
         }
         if (!(s.archStepDeg > 0.0))
         {
            notices.push_back("Wire trigger " + GetName() + ": arch_step_deg " + std::to_string(s.archStepDeg) + " cuts no pieces; the default is used. Set arch_step_deg above 0");
            s.archStepDeg = WireSwitchGroup::Settings().archStepDeg;
         }
         if (m_d.m_scaleX != 1.f || m_d.m_scaleY != 1.f)
            notices.push_back("Wire trigger " + GetName() + ": Scale X / Y (" + std::to_string(m_d.m_scaleX) + ", " + std::to_string(m_d.m_scaleY) + ") are not used: a real wire has one size. Set them to 1");
         // the playfield top = the trigger's Surface (user 2026-11-09); the trigger's x, y = the middle of the wire's part above it; the pivot 19.3 mm below it plus heightoffset_mm
         const float floorZ = height;
         const float playfieldZ = m_ptable->GetSurfaceHeight(string(), m_d.m_vCenter.x, m_d.m_vCenter.y);
         const WireShape::Shape shape = ToWireShape(m_d.m_shape);
         const WireShape::Centreline cl = WireShape::MakeCentreline(shape, s.bendStepDeg, s.archStepDeg);
         WireShape::Placement pl;
         pl.x = m_d.m_vCenter.x;
         pl.y = m_d.m_vCenter.y;
         pl.rotDeg = m_d.m_rotation;
         pl.surfaceZ = floorZ;
         pl.heightOffsetMm = s.heightOffsetMm;
         const WireShape::WireBody wb = WireShape::Place(cl, pl, s.wireDMm);
         WireShape::WireBody wm = wb;   // the part's wire, metres (the adapter takes wb.caps, VPU: ActiveBall is the ball nearest it)
         for (WireShape::Capsule &c : wm.caps)
         {
            c.a = c.a * U;
            c.b = c.b * U;
            c.r = c.r * U;
         }
         for (WireShape::V3 *p : { &wm.pivot, &wm.leverEnd, &wm.peak, &wm.knee })
            *p = *p * U;
         wm.r = wb.r * U;
         m_partModel = std::make_shared<WireSwitchModel::WireSwitch>();
         m_partModel->Build(s, (double)floorZ * U, wm, RubberContact::BallSI());
         // the drawn wire is the physics' wire (requirements.md section 4): a tube of the same radius around the same centreline (wireshape.h BuildTube)
         {
            std::vector<float> v;
            std::vector<uint16_t> idx;
            WireShape::BuildTube(wb.nodes, wb.lane, wb.r, 12, v, idx);
            m_tubeVerts.resize(v.size() / 8);
            memcpy(m_tubeVerts.data(), v.data(), v.size() * sizeof(float));
            m_tubeIdx.assign(idx.begin(), idx.end());
            m_partDrawWorst = 0.;
            for (const Vertex3D_NoTex2 &tv : m_tubeVerts)
               m_partDrawWorst = std::max(m_partDrawWorst, std::fabs(WireShape::DistToChain(WireShape::V3(tv.x, tv.y, tv.z), wb.caps) - wb.r));
         }
         {
            // the drawn switch turns with the part about its hinge (Render); the check (rule 8 log): rest vertices against the physics shape, the rest pose's named points
            const HingeGroup::Part &g = m_partModel->part;
            m_hingeA = Vertex3Ds((float)(g.A.x / U), (float)(g.A.y / U), (float)(g.A.z / U));
            m_hingeE = Vertex3Ds((float)g.E.x, (float)g.E.y, (float)g.E.z);
            m_hingeN = Vertex3Ds((float)g.N.x, (float)g.N.y, (float)g.N.z);
            m_hingeZ = Vertex3Ds((float)g.Z.x, (float)g.Z.y, (float)g.Z.z);
            m_partAngle = (float)g.q;
            m_partDraw = true;
            double topZ = -1e30, nodeTopZ = -1e30;
            for (const Vertex3D_NoTex2 &v : m_tubeVerts)
               topZ = std::max(topZ, (double)v.z);
            for (const WireShape::V3 &n : wb.nodes)
               nodeTopZ = std::max(nodeTopZ, n.z);
            const auto mmAbove = [&](const double z) { return (z - floorZ) * WireShape::kMmPerVpu; };
            PLOGI << "Trigger " << GetName() << ": drawn tube " << m_tubeVerts.size() << " vertices, off the physics shape by at most " << m_partDrawWorst * WireShape::kMmPerVpu * 1e3 << " um; drawn top "
                  << (topZ - playfieldZ) * WireShape::kMmPerVpu << " mm above the playfield (" << mmAbove(topZ) << " mm above its surface); centreline top " << mmAbove(nodeTopZ) << " mm";
            PLOGI << "Trigger " << GetName() << ": rest pose above its surface, mm: pivot " << mmAbove(wb.pivot.z) << ", lever end " << mmAbove(wb.leverEnd.z) << ", peak " << mmAbove(wb.peak.z)
                  << (cl.hasKnee ? ", knee " + std::to_string(mmAbove(wb.knee.z)) : std::string()) << "; pivot at (" << wb.pivot.x << ", " << wb.pivot.y << ") VPU, the wire toward (" << wb.lane.x << ", " << wb.lane.y << ")";
         }
         RubberObjectSettings hs;
         ContactHost::LoadHardSurfaceSettings(hs, ini, "trigger", GetName(), s.elasticity, s.friction);
         const RubberContact::V3 &A = m_partModel->gp.A;
         const double reach = m_partModel->Reach();
         const float lo[3] = { (float)((A.x - reach) / U), (float)((A.y - reach) / U), (float)((A.z - reach) / U) };
         const float hi[3] = { (float)((A.x + reach) / U), (float)((A.y + reach) / U), (float)((A.z + reach) / U) };
         auto adapter = std::make_unique<WireSwitchAdapter>(this, m_partModel, wb.caps);
         WireSwitchAdapter *const adapterRaw = adapter.get();
         m_partHost = new ContactPartHost(this, m_partModel, std::move(adapter), hs, ini, GetName(), "trigger", lo, hi);
         adapterRaw->m_host = m_partHost;
         // a wire whose peak ends below the playfield under it can never be reached (Gorgar's old triggerWall at -25, a workaround that sank the drawn wire)
         if (wb.peak.z < playfieldZ)
            notices.push_back("Wire trigger " + GetName() + " sits below the playfield: no ball can reach it; check its Surface and heightoffset_mm");
         for (const std::string &n : notices)
            PLOGW << n;
         adapterRaw->m_notices = std::move(notices);
         m_partHost->SetHardSurface(s.elasticity, s.friction);
         physics->AddCollider(m_partHost, isUI);
         const WireSwitchGroup::Part &w = m_partModel->part;
         PLOGI << "Trigger " << GetName() << ": wire switch part on the group step (shape Wire" << "ABCDE"[(int)shape] << ", " << wm.caps.size() << " capsules, at most " << cl.sagittaMm * 1e3 << " um off the true bends, wire " << 2e3 * wm.r << " mm: hinge body with its mass)";
         PLOGI << "Trigger " << GetName() << ": surface '" << m_d.m_szSurface << "' at z " << floorZ << " VPU (playfield z " << playfieldZ << "), heightoffset " << s.heightOffsetMm << " mm; spring " << s.springN
               << " N at first touch, " << s.springFullN << " N at full travel; close " << s.closeMm << " mm, release " << s.releaseMm << " mm, full travel " << s.stopBelowMm << " mm under the playfield; wire_d_mm " << s.wireDMm << ", bend_step_deg " << s.bendStepDeg << ", arch_step_deg " << s.archStepDeg << "; mass "
               << w.mp.m * 1e3 << " g, I " << w.mp.I << " kg m^2, reach " << reach * 1e3 << " mm)";
         return;
      }
   }
   if (!isUI && ContactSystem::Unified() && m_d.m_shape == TriggerButton)
   {
      // the rollover button (runs/2026-11-10-rollover/design.md): a plastic cap on its leaf switch, a part of the group step; no stock sensor (every Button on a unified table)
      constexpr double U = ContactSystem::kMetresPerUnit;
      const std::filesystem::path ini = ContactHost::IniPathFor(m_ptable->m_filename);
      ButtonCap::Settings s;
      std::string ignoredKeys;
      PartsIni::LoadButtonCap(s, ini, GetName(), &ignoredKeys);
      std::vector<std::string> notices;   // load-time detectors (rule 8): on screen at the first cycle and in the log as WARN
      if (!ignoredKeys.empty())
         notices.push_back("Button trigger " + GetName() + ": ini key(s) " + ignoredKeys + " not used (removed or unknown); remove them from the table's physics ini");
      const double radiusMm = m_d.m_radius * U * 1e3, fitMm = std::min(radiusMm, s.restMm);   // design 4f: the corner fits under the cap's radius and its height above the wood
      if (s.cornerMm > fitMm)
      {
         notices.push_back("Button trigger " + GetName() + ": corner_mm " + std::to_string(s.cornerMm) + " does not fit the cap (radius " + std::to_string(radiusMm) + " mm, rest_mm " + std::to_string(s.restMm) + "); " + std::to_string(fitMm) + " is used");
         s.cornerMm = fitMm;
      }
      if (s.capHMm <= s.restMm)
      {
         notices.push_back("Button trigger " + GetName() + ": cap_h_mm " + std::to_string(s.capHMm) + " does not reach into the hole at rest (it must exceed rest_mm); the default 3.97 mm is used");
         s.capHMm = ButtonCap::Settings().capHMm;
      }
      if (s.springN <= s.massG * 1e-3 * 9.81 || s.springFullN < s.springN)
      {
         notices.push_back("Button trigger " + GetName() + ": spring_n " + std::to_string(s.springN) + " / spring_full_n " + std::to_string(s.springFullN) + " do not hold the cap up (the preload must exceed its weight, the full force the preload); the defaults 0.34 / 0.69 N are used");
         s.springN = ButtonCap::Settings().springN;
         s.springFullN = ButtonCap::Settings().springFullN;
      }
      const GroupStep::V3 axis((double)m_d.m_vCenter.x * U, (double)m_d.m_vCenter.y * U, (double)height * U);
      auto model = std::make_shared<ButtonCap::Cap>();
      model->Build(s, axis, m_d.m_radius * U, RubberContact::BallSI());
      // the drawn cap: the stock mesh (taller than the real cap), its top put at the physical top at rest, moved by the travel (looks only)
      float meshTop = -FLT_MAX;
      for (const Vertex3D_NoTex2 &v : SetupMeshData(TriggerButton).vertices)
         meshTop = max(meshTop, v.z * m_d.m_radius + 5.0f);   // (GenerateMesh's zoffset for a Button)
      const double drawRest = s.restMm * 1e-3 / U - meshTop;
      RubberObjectSettings hs;
      const JetBumper::Params pp;
      ContactHost::LoadHardSurfaceSettings(hs, ini, "trigger", GetName(), pp.plasticE, pp.plasticMu);
      const double reach = model->Reach();
      const float lo[3] = { (float)((axis.x - reach) / U), (float)((axis.y - reach) / U), (float)((axis.z - reach) / U) };
      const float hi[3] = { (float)((axis.x + reach) / U), (float)((axis.y + reach) / U), (float)((axis.z + reach) / U) };
      auto adapter = std::make_unique<ButtonCapAdapter>(this, model, drawRest);
      ButtonCapAdapter *const adapterRaw = adapter.get();
      m_partHost = new ContactPartHost(this, model, std::move(adapter), hs, ini, GetName(), "trigger", lo, hi);
      for (const std::string &n : notices)
         PLOGW << n;
      adapterRaw->m_notices = std::move(notices);
      m_partHost->SetHardSurface(pp.plasticE, pp.plasticMu);
      physics->AddCollider(m_partHost, isUI);
      m_partAngle = (float)drawRest;
      const ButtonCap::Part &c = model->part;
      PLOGI << "Trigger " << GetName() << ": button cap on the group step at (" << m_d.m_vCenter.x << ", " << m_d.m_vCenter.y << ") VPU, surface z " << height << " VPU; radius " << c.R * 1e3 << " mm, rest " << s.restMm
            << " mm, corner " << s.cornerMm << " mm, close " << s.closeMm << " mm; spring " << s.springN << " N at rest, " << s.springFullN << " N flush; mass " << s.massG << " g; stops e " << s.stopE << " (sink only, one coordinate; stop k " << c.stops[c.stopUp].k << " N/m^1.5); cap " << s.capHMm << " mm tall; plastic E " << c.rp.eMpa
            << " MPa, mu " << c.rp.mu << "; enabled " << m_d.m_enabled << "; drawn offset at rest " << drawRest << " VPU";
      return;
   }
   if (!isUI && ContactSystem::Unified() && m_d.m_shape == TriggerNone)
   {
      // the sensing zone (design.md 1, 3): no stock hit objects, no ball move; tested on the unified engine's chords
      vector<RenderVertex> vvertex;
      m_curve.GetRgVertex(vvertex);
      vector<Vertex2D> poly(vvertex.size());
      for (size_t i = 0; i < vvertex.size(); i++)
         poly[i] = Vertex2D(vvertex[i].x, vvertex[i].y);
      PLOGI << "Trigger " << GetName() << ": None sensing zone on the unified engine, " << poly.size() << " corners, z " << height << " to " << height << " + HitHeight (" << m_d.m_hit_height << ") VPU, centre test, enabled " << m_d.m_enabled;
      physics->m_noneZones.Add(this, std::move(poly), height);
      return;
   }
   if (m_d.m_shape == TriggerStar || m_d.m_shape == TriggerButton)
   {
      TriggerHitCircle* hitcircle = new TriggerHitCircle(this, m_d.m_vCenter, m_d.m_radius, height, height + m_d.m_hit_height);
      hitcircle->m_enabled = isUI ? true : m_d.m_enabled;
      hitcircle->m_ObjType = eTrigger;
      hitcircle->m_obj = (IFireEvents *)this;
      if (!isUI) // Keep reference to update enable flag
      {
         assert(m_ptriggerhitcircle == nullptr);
         m_ptriggerhitcircle = hitcircle;
      }
      physics->AddCollider(hitcircle, isUI);
   }
   else
   {
      vector<RenderVertex> vvertex;
      m_curve.GetRgVertex(vvertex);
      // The volume hit objects derive inside/outside from the winding dependent
      // segment and polygon normals, so restore the canonical winding: a reversed
      // point order would swap the Hit and UnHit events
      NormalizeWindingOrder(vvertex);

      const int count = (int)vvertex.size();
      for (int i = 0; i < count; i++)
      {
         const RenderVertex &pv1 = vvertex[(i < count - 1) ? (i + 1) : 0];
         const RenderVertex &pv2 = vvertex[(i < count - 2) ? (i + 2) : (i + 2 - count)];

         LineSeg *const plineseg = isUI ? 
              new LineSeg(this, pv1, pv2, height, height + max(m_d.m_hit_height - 8.0f, 0.f)) //adjust for same hit height as circular
            : new TriggerLineSeg(this, pv1, pv2, height, height + max(m_d.m_hit_height - 8.0f, 0.f));
         plineseg->m_ObjType = eTrigger;
         plineseg->m_obj = (IFireEvents *)this;
         physics->AddCollider(plineseg, isUI);
      }

      Vertex3Ds *const rgv3D = new Vertex3Ds[count];
      for (int i = 0; i < count; i++)
      {
         rgv3D[i].x = vvertex[i].x;
         rgv3D[i].y = vvertex[i].y;
         rgv3D[i].z = height + (float)(PHYS_SKIN * 2.0);
      }
      Hit3DPoly *const ph3dpoly = new Hit3DPoly(this, rgv3D, count);
      ph3dpoly->m_ObjType = eTrigger;
      ph3dpoly->m_obj = (IFireEvents *)this;
      if (!isUI) // Keep reference to update enable flag
      {
         assert(m_ptriggerhitpoly == nullptr);
         m_ptriggerhitpoly = ph3dpoly;
      }
      physics->AddCollider(ph3dpoly, isUI);
   }
}

void Trigger::PhysicRelease(PhysicsEngine* physics, const bool isUI)
{
   if (!isUI)
   {
      physics->m_noneZones.Remove(this);
      m_ptriggerhitcircle = nullptr;
      m_ptriggerhitpoly = nullptr;
      m_partModel = nullptr;   // (the host is owned by the physics engine)
      m_partHost = nullptr;
      m_partDraw = false;
      m_tubeVerts.clear();
      m_tubeIdx.clear();
      m_partAngle = 0.f;
   }
}

void Trigger::TriggerAnimationHit()
{
   m_hitEvent = true;
}

void Trigger::TriggerAnimationUnhit()
{
   m_unhitEvent = true;
}

void Trigger::OnHit(const bool hit)
{
   TestDriver::TriggerEvent(GetName(), hit);   // (test driver: events.csv, every trigger type, both engines)
   if (!m_d.m_stateSrcLink.empty() && g_pplayer != nullptr)
      g_pplayer->m_resURIResolver.SetBoolState(m_d.m_stateSrcLink, hit);
   FireGroupEvent(hit ? DISPID_HitEvents_Hit : DISPID_HitEvents_Unhit);
}

#pragma endregion


#pragma region Rendering

void Trigger::RenderSetup(Renderer *renderer)
{
   assert(m_renderer == nullptr);
   m_renderer = renderer;
   m_hitEvent = false;
   m_unhitEvent = false;
   m_doAnimation = false;
   m_moveDown = false;
   m_animHeightOffset = 0.0f;
   UpdateMesh();
}

void Trigger::UpdateMesh()
{
   m_meshBuffer = nullptr;
   m_triggerVertices = nullptr;
   m_meshShape = m_d.m_shape;
   m_vertexBuffer_animHeightOffset = -FLT_MAX;

   if (m_d.m_shape == TriggerNone)
      return;

   if (!m_partDraw)   // no switch part: the stock mesh, or for WireE (no stock mesh) its rest tube
   {
      m_tubeVerts.clear();
      m_tubeIdx.clear();
      if (m_d.m_shape == TriggerWireE)
         MakeRestTube();
   }
   const std::span<const WORD> indices = m_tubeVerts.empty() ? SetupMeshData(m_d.m_shape).indices : std::span<const WORD>(m_tubeIdx);
   if (indices.empty())
      return;

   m_triggerVertices = GenerateMesh(m_boundingSphereCenter);
   if (m_triggerVertices == nullptr)
      return;
   std::shared_ptr<IndexBuffer> triggerIndexBuffer = std::make_shared<IndexBuffer>(m_renderer->m_renderDevice, (unsigned int)indices.size(), indices.data());
   std::shared_ptr<VertexBuffer> vertexBuffer = std::make_shared<VertexBuffer>(m_renderer->m_renderDevice, (unsigned int)m_triggerVertices->size(), (float *)m_triggerVertices->data(), true);
   m_meshBuffer = std::make_shared<MeshBuffer>(GetName(), vertexBuffer, triggerIndexBuffer, true);
}

void Trigger::RenderRelease()
{
   assert(m_renderer != nullptr);
   m_renderer = nullptr;
   m_meshBuffer = nullptr;
   m_triggerVertices = nullptr;
}

// Ported at: VisualPinball.Unity/VisualPinball.Unity/VPT/Trigger/TriggerAnimationSystem.cs

void Trigger::UpdateAnimation(const float diff_time_msec)
{
   assert(m_renderer != nullptr);
   float animLimit;
   switch (m_d.m_shape)
   {
      case TriggerStar:   animLimit = m_d.m_radius * (float)(1.0 /  5.0); break;
      case TriggerButton: animLimit = m_d.m_radius * (float)(1.0 / 10.0); break;
      case TriggerWireC:  animLimit = 60.0f; break;
      case TriggerWireD:  animLimit = 25.0f; break;
      case TriggerInder:  animLimit = 25.0f; break;
      default:            animLimit = 32.0f; break;
   }

   if (m_partDraw)   // a converted switch moves by its hinge angle (Render), not the stock slide
   {
      m_hitEvent = m_unhitEvent = m_doAnimation = false;
      m_animHeightOffset = 0.0f;
      return;
   }
   if (m_partHost != nullptr && m_d.m_shape == TriggerButton)   // a button cap: the stock mesh at the cap's travel (its adapter publishes the drawn offset in m_partAngle), not the timed slide
   {
      m_hitEvent = m_unhitEvent = m_doAnimation = false;
      m_animHeightOffset = m_partAngle;
      return;
   }

   if (m_hitEvent)
   {
      m_doAnimation = true;
      m_hitEvent = false;
      // unhitEvent = false;   // Bugfix: If HitEvent and unhitEvent happen at the same time, you want to favor the unhit, otherwise the switch gets stuck down.
      m_animHeightOffset = 0.0f;
      m_moveDown = true;
   }
   if (m_unhitEvent)
   {
      m_doAnimation = true;
      m_unhitEvent = false;
      m_hitEvent = false;
      m_animHeightOffset = animLimit;
      m_moveDown = false;
   }

   if (m_doAnimation)
   {
      float step = diff_time_msec*m_d.m_animSpeed;
      if (m_moveDown)
         step = -step;
      m_animHeightOffset += step;

      if (m_moveDown)
      {
         if (m_animHeightOffset <= -animLimit)
         {
            m_animHeightOffset = -animLimit;
            m_doAnimation = false;
            m_moveDown = false;
         }
      }
      else
      {
         if (m_animHeightOffset >= 0.0f)
         {
            m_animHeightOffset = 0.0f;
            m_doAnimation = false;
            m_moveDown = true;
         }
      }
      FireGroupEvent(DISPID_AnimateEvents_Animate);
   }
}

void Trigger::Render(const unsigned int renderMask)
{
   assert(m_renderer != nullptr);
   assert(!m_desktopBackdrop);
   const bool isStaticOnly = renderMask & Renderer::STATIC_ONLY;
   const bool isDynamicOnly = renderMask & Renderer::DYNAMIC_ONLY;
   const bool isReflectionPass = renderMask & Renderer::REFLECTION_PASS;
   TRACE_FUNCTION();

   if (!m_d.m_visible 
   || m_d.m_shape == TriggerNone 
   || isStaticOnly 
   || (isReflectionPass && !m_d.m_reflectionEnabled))
      return;

   if (m_meshBuffer == nullptr || m_meshShape != m_d.m_shape)
      UpdateMesh();
   if (m_meshBuffer == nullptr)
      return;

   if (m_animHeightOffset != m_vertexBuffer_animHeightOffset)
   {
      m_vertexBuffer_animHeightOffset = m_animHeightOffset;
      const Vertex3D_NoTex2 *const verts = m_triggerVertices->data();
      Vertex3D_NoTex2 *buf;
      m_meshBuffer->m_vb->Lock(buf);
      for (unsigned int i = 0; i < m_meshBuffer->m_vb->m_count; i++)
      {
         buf[i].x = verts[i].x;
         buf[i].y = verts[i].y;
         buf[i].z = verts[i].z + m_animHeightOffset;
         buf[i].nx = verts[i].nx;
         buf[i].ny = verts[i].ny;
         buf[i].nz = verts[i].nz;
         buf[i].tu = verts[i].tu;
         buf[i].tv = verts[i].tv;
      }
      m_meshBuffer->m_vb->Unlock();
   }

   m_renderer->m_renderDevice->ResetRenderState();
   if (CentrelineShape(m_d.m_shape) || m_d.m_shape == TriggerInder)
      m_renderer->m_renderDevice->SetRenderState(RenderState::CULLMODE, RenderState::CULL_NONE);
   m_renderer->m_renderDevice->m_basicShader->SetBasic(m_ptable->GetMaterial(m_d.m_szMaterial), nullptr);
   if (m_partDraw)
   {
      // a converted switch: the rest mesh turned rigidly by the part's angle about its hinge, as hingegroup.h turns the body (Down(th) = -Z cos + N sin, Front(th) = N cos + Z sin):
      // world = A + E (r.E) + Down (-r.Z) + Front (r.N), one matrix per frame (row vectors), no vertex rewrite
      const float c = cosf(m_partAngle), s = sinf(m_partAngle);
      const Vertex3Ds &E = m_hingeE, &N = m_hingeN, &Z = m_hingeZ, &A = m_hingeA;
      const Vertex3Ds D = Z * -c + N * s, F = N * c + Z * s;
      float m[3][3];
      for (int i = 0; i < 3; i++)
         for (int j = 0; j < 3; j++)
            m[i][j] = (&E.x)[i] * (&E.x)[j] - (&Z.x)[i] * (&D.x)[j] + (&N.x)[i] * (&F.x)[j];
      float t[3];
      for (int j = 0; j < 3; j++)
         t[j] = (&A.x)[j] - (A.x * m[0][j] + A.y * m[1][j] + A.z * m[2][j]);
      m_renderer->UpdateBasicShaderMatrix(Matrix3D(m[0][0], m[0][1], m[0][2], 0.f, m[1][0], m[1][1], m[1][2], 0.f, m[2][0], m[2][1], m[2][2], 0.f, t[0], t[1], t[2], 1.f));
   }
   m_renderer->m_renderDevice->DrawMesh(
      m_renderer->m_renderDevice->m_basicShader, false, m_boundingSphereCenter, 0.f, m_meshBuffer, RenderDevice::TRIANGLELIST, 0, m_meshBuffer->m_ib->m_count);
   if (m_partDraw)
      m_renderer->UpdateBasicShaderMatrix();
}

#pragma endregion


void Trigger::ExportMesh(ObjLoader& loader)
{
   if (!m_d.m_visible || m_d.m_shape == TriggerNone)
      return;

   const string& name = m_name;
   Vertex3Ds boundingSphereCenter;
   const auto triggerVertices = GenerateMesh(boundingSphereCenter);
   if (!triggerVertices)
      return;
   const std::span<const WORD> indices = m_tubeVerts.empty() ? SetupMeshData(m_d.m_shape).indices : std::span<const WORD>(m_tubeIdx);
   loader.WriteObjectName(name);
   loader.WriteVertexInfo(triggerVertices->data(), (unsigned int)triggerVertices->size());
   const Material * const mat = m_ptable->GetMaterial(m_d.m_szMaterial);
   loader.WriteMaterial(m_d.m_szMaterial, string(), mat);
   loader.UseTexture(m_d.m_szMaterial);

   loader.WriteFaceInfoList(indices.data(), (unsigned int)indices.size());
   loader.UpdateFaceOffset((unsigned int)triggerVertices->size());
}

// Ported at: VisualPinball.Engine/VPT/Trigger/TriggerMeshGenerator.cs

std::unique_ptr<std::vector<Vertex3D_NoTex2>> Trigger::GenerateMesh(Vertex3Ds &boundingSphereCenter) const
{
   const float baseHeight = m_ptable->GetSurfaceHeight(m_d.m_szSurface, m_d.m_vCenter.x, m_d.m_vCenter.y);
   float zoffset = (m_d.m_shape == TriggerButton) ? 5.0f : 0.0f;
   if (m_d.m_shape == TriggerWireC) zoffset = -19.0f;

   boundingSphereCenter.Set(m_d.m_vCenter.x, m_d.m_vCenter.y, baseHeight);

   if (!m_tubeVerts.empty())   // WireA-E: the tube around the physics' centreline, already placed (PhysicSetup, MakeRestTube)
      return std::make_unique<std::vector<Vertex3D_NoTex2>>(m_tubeVerts);

   const StaticMeshData meshData = SetupMeshData(m_d.m_shape);
   const std::span<const Vertex3D_NoTex2> &verts = meshData.vertices;
   if (verts.empty())
   {
      ShowError("Unknown Trigger type");
      return nullptr;
   }

   auto triggerVertices = std::make_unique<std::vector<Vertex3D_NoTex2>>(verts.size());
   Vertex3D_NoTex2 *const outVerts = triggerVertices->data();

   Matrix3D fullMatrix;
   if (m_d.m_shape == TriggerWireB)
   {
      fullMatrix = Matrix3D::MatrixRotateX(ANGTORAD(-23.f))
                 * Matrix3D::MatrixRotateZ(ANGTORAD(m_d.m_rotation));
   }
   else if (m_d.m_shape == TriggerWireC)
   {
      fullMatrix = Matrix3D::MatrixRotateX(ANGTORAD(140.f))
                 * Matrix3D::MatrixRotateZ(ANGTORAD(m_d.m_rotation));
   }
   else
      fullMatrix = Matrix3D::MatrixRotateZ(ANGTORAD(m_d.m_rotation));

   for (size_t i = 0; i < verts.size(); i++)
   {
      Vertex3Ds vert = fullMatrix * Vertex3Ds{verts[i].x, verts[i].y, verts[i].z};

      if (m_d.m_shape == TriggerButton || m_d.m_shape == TriggerStar)
      {
         outVerts[i].x = vert.x * m_d.m_radius + m_d.m_vCenter.x;
         outVerts[i].y = vert.y * m_d.m_radius + m_d.m_vCenter.y;
         outVerts[i].z = vert.z * m_d.m_radius + baseHeight + zoffset;
      }
      else
      {
         outVerts[i].x = vert.x * m_d.m_scaleX + m_d.m_vCenter.x;
         outVerts[i].y = vert.y * m_d.m_scaleY + m_d.m_vCenter.y;
         outVerts[i].z = vert.z * 1.0f + baseHeight + zoffset;
      }

      vert = Vertex3Ds(verts[i].nx, verts[i].ny, verts[i].nz);
      vert = fullMatrix.MultiplyVectorNoTranslate(vert);
      outVerts[i].nx = vert.x;
      outVerts[i].ny = vert.y;
      outVerts[i].nz = vert.z;
      outVerts[i].tu = verts[i].tu;
      outVerts[i].tv = verts[i].tv;

      if (m_d.m_shape == TriggerWireA || m_d.m_shape == TriggerWireB || m_d.m_shape == TriggerWireC || m_d.m_shape == TriggerWireD || m_d.m_shape == TriggerInder)
      {
         outVerts[i].x += outVerts[i].nx * m_d.m_wireThickness;
         outVerts[i].y += outVerts[i].ny * m_d.m_wireThickness;
         outVerts[i].z += outVerts[i].nz * m_d.m_wireThickness;
      }
   }
   return triggerVertices;
}

void Trigger::MakeRestTube()
{
   const WireSwitchGroup::Settings ds;   // the default steps and wire (the table's ini is read at play)
   const WireShape::Centreline cl = WireShape::MakeCentreline(ToWireShape(m_d.m_shape), ds.bendStepDeg, ds.archStepDeg);
   WireShape::Placement pl;
   pl.x = m_d.m_vCenter.x;
   pl.y = m_d.m_vCenter.y;
   pl.rotDeg = m_d.m_rotation;
   pl.surfaceZ = m_ptable->GetSurfaceHeight(m_d.m_szSurface, m_d.m_vCenter.x, m_d.m_vCenter.y);
   const WireShape::WireBody wb = WireShape::Place(cl, pl, WireSwitchGroup::Settings().wireDMm);
   std::vector<float> v;
   std::vector<uint16_t> idx;
   WireShape::BuildTube(wb.nodes, wb.lane, wb.r, 12, v, idx);
   m_tubeVerts.resize(v.size() / 8);
   memcpy(m_tubeVerts.data(), v.data(), v.size() * sizeof(float));
   m_tubeIdx.assign(idx.begin(), idx.end());
}

void Trigger::FlipX(const Vertex2D &pvCenter)
{
   m_curve.FlipPointX(pvCenter);
   const float deltax = m_d.m_vCenter.x - pvCenter.x;
   m_d.m_vCenter.x -= deltax * 2.0f;
}

void Trigger::FlipY(const Vertex2D &pvCenter)
{
   m_curve.FlipPointY(pvCenter);
   const float deltay = m_d.m_vCenter.y - pvCenter.y;
   m_d.m_vCenter.y -= deltay * 2.0f;
}

void Trigger::Rotate(const float ang, const Vertex2D &center, const bool useElementCenter)
{
   m_curve.RotatePoints(ang, useElementCenter ? GetCenter() : center);
   if (!useElementCenter)
   {
      const float sn = sinf(ANGTORAD(ang));
      const float cs = cosf(ANGTORAD(ang));
      const float dx = m_d.m_vCenter.x - center.x;
      const float dy = m_d.m_vCenter.y - center.y;
      const float dx2 = cs * dx - sn * dy;
      const float dy2 = cs * dy + sn * dx;
      m_d.m_vCenter.x = center.x + dx2;
      m_d.m_vCenter.y = center.y + dy2;
   }
   m_d.m_rotation += ang;
}

void Trigger::Scale(const float scalex, const float scaley, const Vertex2D &center, const bool useElementCenter)
{
   m_curve.ScalePoints(scalex, scaley, useElementCenter ? GetCenter() : center);
   if (!useElementCenter)
   {
      const float dx = (m_d.m_vCenter.x - center.x) * scalex;
      const float dy = (m_d.m_vCenter.y - center.y) * scaley;
      m_d.m_vCenter.x = center.x + dx;
      m_d.m_vCenter.y = center.y + dy;
   }
   m_d.m_scaleX *= scalex;
   m_d.m_scaleY *= scaley;
}

void Trigger::Translate(const Vertex2D &offset)
{
   m_curve.TranslatePoints(offset);
   m_d.m_vCenter.x += offset.x;
   m_d.m_vCenter.y += offset.y;
}

void Trigger::Save(IObjectWriter& writer, const bool saveForUndo)
{
   writer.WriteVector2(FID(VCEN), m_d.m_vCenter);
   writer.WriteFloat(FID(RADI), m_d.m_radius);
   writer.WriteFloat(FID(ROTA), m_d.m_rotation);
   writer.WriteFloat(FID(WITI), m_d.m_wireThickness);
   writer.WriteFloat(FID(SCAX), m_d.m_scaleX);
   writer.WriteFloat(FID(SCAY), m_d.m_scaleY);
   writer.WriteBool(FID(TMON), m_timerEnabled);
   writer.WriteInt(FID(TMIN), m_timerInterval);
   writer.WriteString(FID(SURF), m_d.m_szSurface);
   writer.WriteString(FID(MATR), m_d.m_szMaterial);
   writer.WriteBool(FID(EBLD), m_d.m_enabled);
   writer.WriteString(FID(LINK), m_d.m_stateSrcLink);
   writer.WriteBool(FID(VSBL), m_d.m_visible);
   writer.WriteFloat(FID(THOT), m_d.m_hit_height);
   writer.WriteWideString(FID(NAME), MakeWString(m_name));
   writer.WriteInt(FID(SHAP), m_d.m_shape);
   writer.WriteFloat(FID(ANSP), m_d.m_animSpeed);
   writer.WriteBool(FID(REEN), m_d.m_reflectionEnabled);
   SaveSharedEditableFields(writer);
   m_curve.SavePoints(writer);
   writer.EndObject();
}

void Trigger::ClearForOverwrite() { m_curve.ClearPoints(); }

void Trigger::Load(IObjectReader& reader)
{
   SetDefaults(false);
   reader.AsObject(
      [this](int tag, IObjectReader& reader)
      {
         switch (tag)
         {
         case FID(PIID): reader.AsInt(); break;
         case FID(VCEN): m_d.m_vCenter = reader.AsVector2(); break;
         case FID(RADI): m_d.m_radius = reader.AsFloat(); break;
         case FID(ROTA): m_d.m_rotation = reader.AsFloat(); break;
         case FID(WITI): m_d.m_wireThickness = reader.AsFloat(); break;
         case FID(SCAX): m_d.m_scaleX = reader.AsFloat(); break;
         case FID(SCAY): m_d.m_scaleY = reader.AsFloat(); break;
         case FID(MATR): m_d.m_szMaterial = reader.AsString(); break;
         case FID(TMON): m_timerEnabled = reader.AsBool(); break;
         case FID(TMIN): m_timerInterval = reader.AsInt(); break;
         case FID(SURF): m_d.m_szSurface = reader.AsString(); break;
         case FID(EBLD): m_d.m_enabled = reader.AsBool(); break;
         case FID(LINK): m_d.m_stateSrcLink = reader.AsString(); break;
         case FID(THOT): m_d.m_hit_height = reader.AsFloat(); break;
         case FID(VSBL): m_d.m_visible = reader.AsBool(); break;
         case FID(REEN): m_d.m_reflectionEnabled = reader.AsBool(); break;
         case FID(SHAP): m_d.m_shape = static_cast<TriggerShape>(reader.AsInt()); break;
         case FID(ANSP): m_d.m_animSpeed = reader.AsFloat(); break;
         case FID(NAME): m_name = MakeString(reader.AsWideString()); break;
         case FID(DPNT): m_curve.LoadPointToken(reader); break;
         default: LoadSharedEditableField(tag, reader); break;
         }
         return true;
      });

   // Seed the default shape for tables saved without drag points
   if (m_curve.GetPoints().empty())
      InitShape(m_d.m_vCenter.x, m_d.m_vCenter.y);
}

STDMETHODIMP Trigger::InterfaceSupportsErrorInfo(REFIID riid)
{
   static const IID* arr[] =
   {
      &IID_ITrigger,
   };

   for (size_t i = 0; i < std::size(arr); i++)
      if (InlineIsEqualGUID(*arr[i], riid))
         return S_OK;

   return S_FALSE;
}

STDMETHODIMP Trigger::get_Radius(float *pVal)
{
   *pVal = m_d.m_radius;
   return S_OK;
}

STDMETHODIMP Trigger::put_Radius(float newVal)
{
   m_d.m_radius = newVal;
   return S_OK;
}

STDMETHODIMP Trigger::get_X(float *pVal)
{
   *pVal = m_d.m_vCenter.x;
   return S_OK;
}

STDMETHODIMP Trigger::put_X(float newVal)
{
   m_d.m_vCenter.x = newVal;
   return S_OK;
}

STDMETHODIMP Trigger::get_Y(float *pVal)
{
   *pVal = m_d.m_vCenter.y;
   return S_OK;
}

STDMETHODIMP Trigger::put_Y(float newVal)
{
   m_d.m_vCenter.y = newVal;
   return S_OK;
}

STDMETHODIMP Trigger::get_Surface(BSTR *pVal)
{
   *pVal = MakeWideBSTR(m_d.m_szSurface);
   return S_OK;
}

STDMETHODIMP Trigger::put_Surface(BSTR newVal)
{
   m_d.m_szSurface = MakeString(newVal);
   return S_OK;
}

STDMETHODIMP Trigger::get_Enabled(VARIANT_BOOL *pVal)
{
   *pVal = FTOVB(m_d.m_enabled);
   return S_OK;
}

STDMETHODIMP Trigger::put_Enabled(VARIANT_BOOL newVal)
{
   const bool was = m_d.m_enabled;
   m_d.m_enabled = VBTOb(newVal);
   if (m_d.m_enabled && !was && g_pplayer != nullptr && g_pplayer->m_physics != nullptr)
      g_pplayer->m_physics->m_noneZones.Recheck(this);   // a None zone enabled again: every ball judged where it is (design 4d)
   if (m_ptriggerhitcircle)
      m_ptriggerhitcircle->m_enabled = m_d.m_enabled;
   if (m_ptriggerhitpoly)
      m_ptriggerhitpoly->m_enabled = m_d.m_enabled;

   return S_OK;
}

STDMETHODIMP Trigger::get_Visible(VARIANT_BOOL *pVal)
{
   *pVal = FTOVB(m_d.m_visible);
   return S_OK;
}

STDMETHODIMP Trigger::put_Visible(VARIANT_BOOL newVal)
{
   m_d.m_visible = VBTOb(newVal);
   return S_OK;
}

STDMETHODIMP Trigger::BallCntOver(int *pVal)
{
   int cnt = 0;

   if (g_pplayer)
   {
      for (size_t i = 0; i < g_pplayer->m_vball.size(); i++)
      {
         Ball *const pball = g_pplayer->m_vball[i];

         if (pball->m_hitBall.m_d.m_vpVolObjs && FindIndexOf(*(pball->m_hitBall.m_d.m_vpVolObjs), (IFireEvents*)this) >= 0) // cast to IFireEvents necessary, as it is stored like this in HitObject.m_obj
         {
            g_pplayer->m_pactiveball = pball; // set active ball for scriptor
            ++cnt;
         }
      }
   }

   *pVal = cnt;

   return S_OK;
}

STDMETHODIMP Trigger::DestroyBall(int *pVal)
{
   if (pVal)
      *pVal = 0;
   
   if (!g_pplayer)
      return S_OK;

   for (Ball *ball : g_pplayer->m_vball)
   {
      if (ball->m_hitBall.m_d.m_vpVolObjs)
      {
         const auto it = std::ranges::find(*(ball->m_hitBall.m_d.m_vpVolObjs), (IFireEvents *)this); // cast to IFireEvents necessary, as it is stored like this in HitObject.m_obj
         if (it != ball->m_hitBall.m_d.m_vpVolObjs->end())
         {
            if (pVal)
               (*pVal) = (*pVal) + 1;
            g_pplayer->m_physics->m_noneZones.BallGone(&ball->m_hitBall);   // (a None zone: the state link opens while the list still names it)
            ball->m_hitBall.m_d.m_vpVolObjs->erase(it);
            g_pplayer->DestroyBall(ball); // inside trigger volume?
         }
      }
   }

   return S_OK;
}

STDMETHODIMP Trigger::get_HitHeight(float *pVal)
{
   *pVal = m_d.m_hit_height;
   return S_OK;
}

STDMETHODIMP Trigger::put_HitHeight(float newVal)
{
   const float was = m_d.m_hit_height;
   m_d.m_hit_height = newVal;
   if (newVal != was && g_pplayer != nullptr && g_pplayer->m_physics != nullptr)
      g_pplayer->m_physics->m_noneZones.Recheck(this);   // a None zone grown or shrunk round a ball: every ball judged where it is (design 4e)
   return S_OK;
}

STDMETHODIMP Trigger::get_Rotation(float *pVal)
{
   *pVal = m_d.m_rotation;
   return S_OK;
}

STDMETHODIMP Trigger::put_Rotation(float newVal)
{
   m_d.m_rotation = newVal;

   return S_OK;
}

STDMETHODIMP Trigger::get_WireThickness(float *pVal)
{
   *pVal = m_d.m_wireThickness;
   return S_OK;
}

STDMETHODIMP Trigger::put_WireThickness(float newVal)
{
   m_d.m_wireThickness = newVal;
   return S_OK;
}

STDMETHODIMP Trigger::get_AnimSpeed(float *pVal)
{
   *pVal = m_d.m_animSpeed;
   return S_OK;
}

STDMETHODIMP Trigger::put_AnimSpeed(float newVal)
{
   m_d.m_animSpeed = newVal;
   return S_OK;
}

STDMETHODIMP Trigger::get_CurrentAnimOffset(float *pVal)
{
   *pVal = m_animHeightOffset;
   return S_OK;
}

STDMETHODIMP Trigger::get_Material(BSTR *pVal)
{
   *pVal = MakeWideBSTR(m_d.m_szMaterial);
   return S_OK;
}

STDMETHODIMP Trigger::put_Material(BSTR newVal)
{
   m_d.m_szMaterial = MakeString(newVal);
   return S_OK;
}

STDMETHODIMP Trigger::get_TriggerShape(TriggerShape *pVal)
{
   *pVal = m_d.m_shape;
   return S_OK;
}

STDMETHODIMP Trigger::put_TriggerShape(TriggerShape newVal)
{
   m_d.m_shape = newVal;

   return S_OK;
}

STDMETHODIMP Trigger::get_ReflectionEnabled(VARIANT_BOOL *pVal)
{
   *pVal = FTOVB(m_d.m_reflectionEnabled);
   return S_OK;
}

STDMETHODIMP Trigger::put_ReflectionEnabled(VARIANT_BOOL newVal)
{
   m_d.m_reflectionEnabled = VBTOb(newVal);
   return S_OK;
}
