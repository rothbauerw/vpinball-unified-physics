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

void Trigger::GetWireOutline(vector<Vertex2D> &outline) const
{
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
         m_model->gp.SetPivotModel(cs.pivotModel);
         if (!m_model->gp.SetPieceSum(cs.hingePieceSum) && cs.hingePieceSum)
            PLOGW << "hinge_piece_sum: part " << m_name << " has more pieces than a contact's terms hold (" << GroupStep::kTermMax << "): its contact stays the deepest piece";
      }
      void Publish() override
      {
         m_model->Collect();
         m_model->bookMs = (int)g_pplayer->m_time_msec;
         m_trigger->m_partAngle = (float)m_model->part.q;   // (the angle the renderer turns the drawn switch by)
         if (!m_notice.empty() && g_pplayer->m_liveUI != nullptr)   // the load-time detector, on screen once the player runs
         {
            g_pplayer->m_liveUI->PushNotification(m_notice, 10000);
            m_notice.clear();
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
      std::string m_notice;                  // a load-time detector's message, shown at the first cycle
   };
}

void Trigger::PhysicSetup(PhysicsEngine* physics, const bool isUI)
{
   m_hitEvent = false;

   if (!isUI && GetPartGroup() != nullptr && GetPartGroup()->GetReferenceSpace() != PartGroupData::SpaceReference::SR_PLAYFIELD)
      return;

   const float height = m_ptable->GetSurfaceHeight(m_d.m_szSurface, m_d.m_vCenter.x, m_d.m_vCenter.y);
   const bool inder = m_d.m_shape == TriggerInder;
   const bool wire = m_d.m_shape == TriggerWireA || m_d.m_shape == TriggerWireB || m_d.m_shape == TriggerWireC || m_d.m_shape == TriggerWireD || inder;
   if (!isUI && ContactSystem::Unified())
   {
      const std::filesystem::path ini = ContactHost::IniPathFor(m_ptable->m_filename);
      const bool on = PhysicsIni::Active(ini).triggers;
      static const void *s_noted = nullptr;   // once per play
      if (s_noted != g_pplayer)
      {
         s_noted = g_pplayer;
         PLOGI << "Triggers: [Features] triggers = " << (on ? 1 : 0) << (on ? ": WireA-D and Inder are rollover switches of the contact system (group step; Inder a flat blade on its mesh triangles); None, Star and Button stay stock sensors (not converted yet)" : ": every trigger is a stock sensor");
      }
      if (on && wire)
      {
         // the real wire: a steel wire on a micro switch lever, hinged below the playfield (design-wire.md); no stock sensor, mover or volume
         constexpr double U = ContactSystem::kMetresPerUnit;
         WireSwitchModel::Settings s;
         if (m_d.m_shape == TriggerWireB)
            s.hingeBackMm = 27.0;   // the user's layouts 2026-11-04 (bench runs/2026-11-04-wire-bcd): B a standard-size micro switch; C held at the low end of its visible strand
         if (m_d.m_shape == TriggerWireC)
            s.leverNode = 1;
         PartsIni::LoadWireSwitch(s, ini, GetName());
         WireShape::Pose pose;
         pose.x = m_d.m_vCenter.x;
         pose.y = m_d.m_vCenter.y;
         pose.rotDeg = m_d.m_rotation;
         pose.scaleX = m_d.m_scaleX;
         // on the trigger's own surface (user 2026-11-04): the arch top peak_mm above it. Detector (rule 8): a wire whose top ends below the playfield under it can never be reached
         // (Gorgar's old triggerWall at -25, a workaround that sank the drawn wire)
         const float floorZ = height;
         const float playfieldZ = m_ptable->GetSurfaceHeight(string(), m_d.m_vCenter.x, m_d.m_vCenter.y);
         pose.surfaceZ = floorZ;
         pose.peakMm = s.peakMm;
         pose.wireDmm = s.wireDMm;
         pose.thicknessVpu = m_d.m_wireThickness;
         const WireShape::Shape shape = m_d.m_shape == TriggerWireB ? WireShape::Shape::WireB : m_d.m_shape == TriggerWireC ? WireShape::Shape::WireC : m_d.m_shape == TriggerWireD ? WireShape::Shape::WireD : WireShape::Shape::WireA;
         std::vector<WireShape::Capsule> chain, adapterChain;   // the part's wire (metres); the adapter's (VPU: ActiveBall is the ball nearest it)
         WireSwitchGroup::Blade blade;
         m_partModel = std::make_shared<WireSwitchModel::WireSwitch>();
         if (inder)
         {
            // a flat steel blade on edge (runs/2026-11-04-inder): the drawn mesh is its collidable shape, the very vertices GenerateMesh draws (rotation, Scale X / Y, Surface, WireThickness)
            Vertex3Ds sphere;
            const std::unique_ptr<std::vector<Vertex3D_NoTex2>> drawn = GenerateMesh(sphere);
            const Vertex3Ds ly = Matrix3D::MatrixRotateZ(ANGTORAD(m_d.m_rotation)).MultiplyVectorNoTranslate(Vertex3Ds(0.f, 1.f, 0.f));   // the mesh's y as drawn
            blade = WireSwitchGroup::PoseBlade(&(*drawn)[0].x, (unsigned)drawn->size(), 8, triggerInderIndices, triggerInderNumIndices, RubberContact::V3(ly.x * m_d.m_scaleX, ly.y * m_d.m_scaleY, 0.0), U);
            WireShape::Capsule top;   // the adapter's: the blade's top edge
            top.a = WireShape::V3(blade.top0.x / U, blade.top0.y / U, blade.top0.z / U);
            top.b = WireShape::V3(blade.top1.x / U, blade.top1.y / U, blade.top1.z / U);
            top.r = 0.5 * blade.thick / U;
            adapterChain.push_back(top);
            m_partModel->Build(s, (double)floorZ * U, blade, RubberContact::BallSI());
            s.peakMm = m_partModel->part.set.peakMm;   // the drawn top above the surface (the detector and the log)
         }
         else
         {
            const WireShape::WireBody wb = shape == WireShape::Shape::WireD ? WireShape::Build(shape, &triggerDWireMesh[0].x, triggerDWireNumVertices, 8, pose)
                                                                             : WireShape::Build(shape, &triggerSimple[0].x, triggerSimpleNumVertices, 8, pose);
            if (wb.caps.empty())
            {
               PLOGW << "Trigger " << GetName() << ": the wire's centreline could not be traced from its mesh: no switch";
               m_partModel.reset();
               return;
            }
            for (const WireShape::Capsule &c : wb.caps)
            {
               WireShape::Capsule m;
               m.a = c.a * U;
               m.b = c.b * U;
               m.r = c.r * U;
               chain.push_back(m);
            }
            adapterChain = wb.caps;
            m_partModel->Build(s, (double)floorZ * U, chain, RubberContact::BallSI());
            // the drawn wire is the physics' wire (requirements.md section 4): the same body transform and tube (wireshape.h PoseDrawn)
            const StaticMeshData md = SetupMeshData(m_d.m_shape);
            m_partDrawn.assign(md.vertices.size() * 8, 0.f);
            m_partDrawWorst = WireShape::PoseDrawn(wb, &md.vertices[0].x, md.vertices.size(), 8, m_partDrawn.data(), 8);
            for (size_t i = 0; i < md.vertices.size(); i++)
            {
               m_partDrawn[i * 8 + 6] = md.vertices[i].tu;
               m_partDrawn[i * 8 + 7] = md.vertices[i].tv;
            }
         }
         {
            // the drawn switch turns with the part about its hinge (Render); the check (rule 8 log): rest vertices against the physics shape, the drawn top above the playfield
            const HingeGroup::Part &g = m_partModel->part;
            m_hingeA = Vertex3Ds((float)(g.A.x / U), (float)(g.A.y / U), (float)(g.A.z / U));
            m_hingeE = Vertex3Ds((float)g.E.x, (float)g.E.y, (float)g.E.z);
            m_hingeN = Vertex3Ds((float)g.N.x, (float)g.N.y, (float)g.N.z);
            m_hingeZ = Vertex3Ds((float)g.Z.x, (float)g.Z.y, (float)g.Z.z);
            m_partAngle = (float)g.q;
            m_partDraw = true;
            Vertex3Ds sphere;
            const std::unique_ptr<std::vector<Vertex3D_NoTex2>> drawn = GenerateMesh(sphere);
            double worst = m_partDrawWorst, topZ = -1e30;
            for (const Vertex3D_NoTex2 &v : *drawn)
            {
               topZ = std::max(topZ, (double)v.z);
               if (inder)   // the blade's faces are the drawn triangles: each drawn vertex against its nearest face corner
               {
                  double d = 1e30;
                  for (const RubberContact::V3 &p : blade.tri)
                     d = std::min(d, (RubberContact::V3(v.x, v.y, v.z) - p / U).Length());
                  worst = std::max(worst, d);
               }
            }
            PLOGI << "Trigger " << GetName() << ": drawn on the physics pose: rest vertices off the physics shape by at most " << worst * WireShape::kMmPerVpu * 1e3 << " um, drawn top "
                  << (topZ - playfieldZ) * WireShape::kMmPerVpu << " mm above the playfield (" << (topZ - floorZ) * WireShape::kMmPerVpu << " mm above its surface)";
         }
         RubberObjectSettings hs;
         ContactHost::LoadHardSurfaceSettings(hs, ini, "trigger", GetName(), s.elasticity, s.friction);
         const RubberContact::V3 &A = m_partModel->gp.A;
         const double reach = m_partModel->Reach();
         const float lo[3] = { (float)((A.x - reach) / U), (float)((A.y - reach) / U), (float)((A.z - reach) / U) };
         const float hi[3] = { (float)((A.x + reach) / U), (float)((A.y + reach) / U), (float)((A.z + reach) / U) };
         auto adapter = std::make_unique<WireSwitchAdapter>(this, m_partModel, adapterChain);
         WireSwitchAdapter *const adapterRaw = adapter.get();
         m_partHost = new ContactPartHost(this, m_partModel, std::move(adapter), hs, ini, GetName(), "trigger", lo, hi);
         adapterRaw->m_host = m_partHost;
         if (floorZ + s.peakMm / WireShape::kMmPerVpu < playfieldZ)
         {
            adapterRaw->m_notice = "Wire trigger " + GetName() + " sits below the playfield: no ball can reach it; check its Surface";
            PLOGW << adapterRaw->m_notice << " (its top at z " << floorZ + s.peakMm / WireShape::kMmPerVpu << " VPU on '" << m_d.m_szSurface << "', the playfield at z " << playfieldZ << ")";
         }
         m_partHost->SetHardSurface(s.elasticity, s.friction);
         physics->AddCollider(m_partHost, isUI);
         const WireSwitchGroup::Part &w = m_partModel->part;
         const WireShape::V3 tip(w.TipWorld().x / U, w.TipWorld().y / U, w.TipWorld().z / U);
         if (inder)
            PLOGI << "Trigger " << GetName() << ": blade switch part on the group step (shape Inder, " << blade.tri.size() / 3 << " mesh triangles as faces, blade " << blade.thick * 1e3 << " mm thick, face skin "
                  << s.edgeMm << " mm";
         else
            PLOGI << "Trigger " << GetName() << ": wire switch part on the group step (shape Wire" << "ABCD"[(int)shape] << ", " << chain.size() << " capsules, wire " << 2e3 * chain[0].r << " mm, wire_model " << s.wireModel
                  << (s.wireModel == 1 ? ": massless lever, angle solved from the balls on it" : ": hinge body with its mass") << ", wire_sum " << s.wireSum << ")";
         PLOGI << "Trigger " << GetName() << ": scale " << m_d.m_scaleX << ", peak " << s.peakMm << " mm above its surface '" << m_d.m_szSurface << "' at z " << floorZ << " VPU (playfield z " << playfieldZ << "); lever end node " << s.leverNode << " at (" << tip.x << ", " << tip.y << ") VPU; pivot " << s.hingeBackMm << " mm back, "
               << s.hingeBelowMm << " mm below; spring " << s.springN << " N at first touch, " << s.springFullN << " N at full travel; close " << s.closeMm << " mm, release " << s.releaseMm << " mm; mass "
               << w.mp.m * 1e3 << " g, reach " << reach * 1e3 << " mm)";
         return;
      }
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
      m_ptriggerhitcircle = nullptr;
      m_ptriggerhitpoly = nullptr;
      m_partModel = nullptr;   // (the host is owned by the physics engine)
      m_partHost = nullptr;
      m_partDraw = false;
      m_partDrawn.clear();
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

   const StaticMeshData meshData = SetupMeshData(m_d.m_shape);
   if (meshData.indices.empty())
      return;

   m_triggerVertices = GenerateMesh(m_boundingSphereCenter);
   if (m_triggerVertices == nullptr)
      return;
   std::shared_ptr<IndexBuffer> triggerIndexBuffer = std::make_shared<IndexBuffer>(m_renderer->m_renderDevice, (unsigned int)meshData.indices.size(), meshData.indices.data());
   std::shared_ptr<VertexBuffer> vertexBuffer = std::make_shared<VertexBuffer>(m_renderer->m_renderDevice, (unsigned int)meshData.vertices.size(), (float *)m_triggerVertices->data(), true);
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
   if (m_d.m_shape == TriggerWireA || m_d.m_shape == TriggerWireB || m_d.m_shape == TriggerWireC || m_d.m_shape == TriggerWireD || m_d.m_shape == TriggerInder)
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
   const StaticMeshData meshData = SetupMeshData(m_d.m_shape);
   loader.WriteObjectName(name);
   loader.WriteVertexInfo(triggerVertices->data(), (unsigned int)meshData.vertices.size());
   const Material * const mat = m_ptable->GetMaterial(m_d.m_szMaterial);
   loader.WriteMaterial(m_d.m_szMaterial, string(), mat);
   loader.UseTexture(m_d.m_szMaterial);

   loader.WriteFaceInfoList(meshData.indices.data(), (unsigned int)meshData.indices.size());
   loader.UpdateFaceOffset((unsigned int)meshData.vertices.size());
}

// Ported at: VisualPinball.Engine/VPT/Trigger/TriggerMeshGenerator.cs

std::unique_ptr<std::vector<Vertex3D_NoTex2>> Trigger::GenerateMesh(Vertex3Ds &boundingSphereCenter) const
{
   const float baseHeight = m_ptable->GetSurfaceHeight(m_d.m_szSurface, m_d.m_vCenter.x, m_d.m_vCenter.y);
   float zoffset = (m_d.m_shape == TriggerButton) ? 5.0f : 0.0f;
   if (m_d.m_shape == TriggerWireC) zoffset = -19.0f;

   boundingSphereCenter.Set(m_d.m_vCenter.x, m_d.m_vCenter.y, baseHeight);

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
      if (m_partDrawn.size() == verts.size() * 8)   // a converted wire (WireA-D): the physics' pose (PhysicSetup, wireshape.h PoseDrawn)
         memcpy(&outVerts[i].x, &m_partDrawn[i * 8], 6 * sizeof(float));
   }
   return triggerVertices;
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
   m_d.m_enabled = VBTOb(newVal);
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
   m_d.m_hit_height = newVal;
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
