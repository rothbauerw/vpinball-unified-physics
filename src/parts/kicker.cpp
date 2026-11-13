// license:GPLv3+

#include "core/stdafx.h"
#include "kicker.h"
#include "physics/contactsurfacehost.h"

#include "core/VPApp.h"
#include "math/matrix.h"
#include "meshes/kickerCupMesh.h"
#include "meshes/kickerCupSys7.h"
#include "meshes/kickerGottlieb.h"
#include "meshes/kickerHitMesh.h"
#include "meshes/kickerHoleMesh.h"
#include "meshes/kickerSimpleHoleMesh.h"
#include "meshes/kickerT1Mesh.h"
#include "meshes/kickerWilliams.h"
#include "renderer/IndexBuffer.h"
#include "renderer/MeshBuffer.h"
#include "renderer/RenderDevice.h"
#include "renderer/Shader.h"
#include "renderer/trace.h"
#include "renderer/VertexBuffer.h"
#include "parts/ball.h"
#include "parts/PartGroup.h"
#include "physics/collideex.h"
#include "physics/PhysicsEngine.h"
#include "physics/carvedholes.h"
#include "physics/contactpartsini.h"
#include "physics/contactsystem.h"
#include "physics/kickerarm.h"
#include "physics/physicsini.h"
#include "physics/testdriver.h"
#include "ui/live/LiveUI.h"
#include <set>
#include "utils/objloader.h"


Kicker::~Kicker()
{
   assert(m_renderer == nullptr);
}

Kicker *Kicker::CopyForPlay() const
{
   STANDARD_EDITABLE_COPY_FOR_PLAY_IMPL(Kicker)
   return dst;
}

HRESULT Kicker::Init(const float x, const float y, const bool fromMouseClick, const bool forPlay)
{
   SetDefaults(fromMouseClick);
   m_d.m_vCenter.x = x;
   m_d.m_vCenter.y = y;
   return S_OK;
}

#define LinkProp(field, prop) field = fromMouseClick ? g_settingsService.GetAppSettings().GetDefaultPropsKicker_##prop() : Settings::GetDefaultPropsKicker_##prop##_Default()
void Kicker::SetDefaults(const bool fromMouseClick)
{
   LinkProp(m_d.m_enabled, Enabled);
   LinkProp(m_d.m_hitAccuracy, HitAccuracy);
   LinkProp(m_d.m_hit_height, HitHeight);
   LinkProp(m_d.m_orientation, Orientation);
   LinkProp(m_d.m_radius, Radius);
   LinkProp(m_d.m_kickertype, KickerType);
   LinkProp(m_d.m_szSurface, Surface);
   LinkProp(m_d.m_fallThrough, FallThrough);
   LinkProp(m_d.m_legacyMode, Legacy);
   LinkProp(m_timerEnabled, TimerEnabled);
   LinkProp(m_timerInterval, TimerInterval);
   SetDefaultPhysics(fromMouseClick);
}

void Kicker::SetDefaultPhysics(const bool fromMouseClick)
{
   LinkProp(m_d.m_scatter, Scatter);
}
#undef LinkProp

void Kicker::WriteRegDefaults()
{
#define LinkProp(field, prop) g_settingsService.GetAppSettings().SetDefaultPropsKicker_##prop(field, false)
   LinkProp(m_d.m_enabled, Enabled);
   LinkProp(m_d.m_hitAccuracy, HitAccuracy);
   LinkProp(m_d.m_hit_height, HitHeight);
   LinkProp(m_d.m_orientation, Orientation);
   LinkProp(m_d.m_radius, Radius);
   LinkProp(m_d.m_kickertype, KickerType);
   LinkProp(m_d.m_szSurface, Surface);
   LinkProp(m_d.m_fallThrough, FallThrough);
   LinkProp(m_d.m_legacyMode, Legacy);
   LinkProp(m_d.m_scatter, Scatter);
   LinkProp(m_timerEnabled, TimerEnabled);
   LinkProp(m_timerInterval, TimerInterval);
#undef LinkProp
}


#pragma region Physics

// Ported at: VisualPinball.Engine/VPT/Kicker/KickerHit.cs

namespace
{
   void KickerNotice(const std::string &msg)   // rule 8: on screen and in the log, once per cause per table start (each caller says it once)
   {
      PLOGW << msg;
      if (g_pplayer != nullptr && g_pplayer->m_liveUI != nullptr)
         g_pplayer->m_liveUI->PushNotification(msg, 10000);
   }
}

// Physical kicker (docs/unified-physics/KickerDesign.txt sections 1, 2, 5, 7; Requirements F 1-8), on the unified engine with [Features] kickers = 1 and type physical: no KickerHitCircle and no
// AddPlayfieldHole. The ball meets the cup, a static hard-surface host built from the cup data compiled in (meshes/kickerCupSys7.h: the cup and the hole's wall below the rim, mm, kicker
// frame), under a round hole carved in the playfield (CarvedHoles: the playfield's faces do not press a ball inside the real cutout's radius, from the cup data; the rim is the hole's
// edge). Gravity, the cup and friction bring the ball to rest; nothing else holds it (1d). Held: the cup's switch, a switch with travel (Requirements F 6): it closes when the ball's centre
// is inside the hole and its bottom reaches the close depth, and opens only when the bottom rises above the open depth (cup data), judged on each step's chord, Hit / Unhit queued and fired
// after the step (ContactSystem::NoneZones with the kicker as owner, 1f, 7e). The kick (KickXYZ) is the arm stroke (4c: the held ball carried along a quarter circle, KickerArm). No Enabled: the cup and its switch are
// always there (F 9). Type virtual (the default) is the virtual kicker (VirtualSetup). Returns false for today's kicker (kickers = 0).
bool Kicker::PhysicalSetup(PhysicsEngine *const physics)
{
   KickerArm::Remove(this);   // (a stroke left from the last table start)
   KickerHold::Remove(this);
   m_physical = m_physicalRefused = m_virtual = m_parkMovedNoted = false;
   m_heldBall = m_lastHeldBall = nullptr;
   m_offsetNoted = m_enabledNoted = m_createNoted = m_inclNoted = false;
   const std::filesystem::path ini = ContactHost::IniPathFor(m_ptable->m_filename);
   const PhysicsIni::Features feat = PhysicsIni::Active(ini);
   if (!feat.kickers)
      return false;
   // [Kicker:<name>] in the table's physics ini, every key read here or by the hard-surface loader below; unknown or repeated keys are said on screen (rule 8)
   bool physical = false;               // type = physical | virtual: default virtual; physical only where the table's physics ini says so (design 7b, Requirements F 9, user 2026-11-12)
   std::string cup = "sys7";            // cup = sys7: the collision cup and its switch depths, compiled in, one per cup type (System 7 the only one so far)
   std::set<std::string> seen;
   std::string unknown, repeated, badType;
   PartsIni::ReadSections(ini, std::string(), "kicker:" + GetName(), [&](const std::string &k, const double v, const std::string &lv) {
      if (!seen.insert(k).second)
         repeated += (repeated.empty() ? "" : ", ") + k;
      if (k == "type")
      {
         if (lv == "physical")
            physical = true;
         else if (lv == "virtual")
            physical = false;
         else
            badType = lv;
      }
      else if (k == "cup")
         cup = lv;
      else if (ContactHost::IsSurfaceKey(k)) {}   // the hard surface's material (ContactHost::LoadHardSurfaceSettings reads it)
      else
         unknown += (unknown.empty() ? "" : ", ") + k;
   });
   const std::string name = GetName();
   if (!unknown.empty())
      KickerNotice("Kicker " + name + ": [Kicker:" + name + "] key(s) " + unknown + " not used (unknown); remove them from the table's physics ini");
   if (!repeated.empty())
      KickerNotice("Kicker " + name + ": [Kicker:" + name + "] key(s) " + repeated + " given more than once; the last one is used. Keep one of each");
   if (!badType.empty())
      KickerNotice("Kicker " + name + ": type '" + badType + "' is neither physical nor virtual; the default (" + (physical ? "physical" : "virtual") + ") is used. Set type = physical or virtual");
   if (!physical)
   {
      VirtualSetup(physics, seen.count("type") != 0);
      return true;
   }
   if (cup != "sys7")
      KickerNotice("Kicker " + name + ": cup '" + cup + "' is not a known cup; the System 7 cup (sys7) is used. Set cup = sys7");
   m_physical = true;
   const float height = m_ptable->GetSurfaceHeight(m_d.m_szSurface, m_d.m_vCenter.x, m_d.m_vCenter.y);
   if (!m_d.m_szSurface.empty() || !feat.playfieldContact)
   {
      // design 6a: the hole is carved in the playfield's force model only; elsewhere there is nothing to carve (stated limit until the cut-hole work). Nothing is built: no stock path
      m_physical = false;
      m_physicalRefused = true;
      KickerNotice("Kicker " + name + ": a physical kicker needs the playfield" + (feat.playfieldContact ? " (it is on " + m_d.m_szSurface + ")" : " on the force model ([Features] playfield_contact = 0)")
                   + "; it does nothing. Set type = virtual in [Kicker:" + name + "]");
      return true;
   }
   constexpr double U = ContactSystem::kMetresPerUnit;
   const double mmU = 1e-3 / U;   // mm -> VPX units
   // the cup host: the cup's faces turned by the kicker's orientation about its centre, at the surface height (as a wall's hard-surface host, parts/surface.cpp)
   const Matrix3D rot = Matrix3D::MatrixRotateZ(ANGTORAD(m_d.m_orientation));
   std::vector<Vertex3Ds> v(kickerCupSys7NumVertices);
   for (unsigned int i = 0; i < kickerCupSys7NumVertices; i++)
   {
      const Vertex3Ds q = rot.MultiplyVectorNoTranslate(Vertex3Ds { (float)(kickerCupSys7Vertices[i][0] * mmU), (float)(kickerCupSys7Vertices[i][1] * mmU), (float)(kickerCupSys7Vertices[i][2] * mmU) });
      v[i] = Vertex3Ds { m_d.m_vCenter.x + q.x, m_d.m_vCenter.y + q.y, height + q.z };
   }
   auto faces = std::make_shared<RubberSurface>();
   for (unsigned int t = 0; t < kickerCupSys7NumTriangles; t++)
      faces->Add(v[kickerCupSys7Indices[t][0]], v[kickerCupSys7Indices[t][1]], v[kickerCupSys7Indices[t][2]]);
   // the material: a hard surface as a wall's, [Surfaces] then [Kicker:<name>] over the stock wall defaults (elasticity 0.3, friction 0.3: a kicker has no physics material; NOT measured)
   constexpr double kE = 0.3, kMu = 0.3;
   RubberObjectSettings rs;
   ContactHost::LoadHardSurfaceSettings(rs, ini, "kicker", name, kE, kMu);
   ContactSurfaceHost *const host = new ContactSurfaceHost(this, faces, rs, ini, name, "kicker");
   host->SetHardSurface(kE, kMu);
   host->m_elasticity = (float)kE;
   host->SetFriction((float)kMu);
   host->m_enabled = true;   // the cup is always there (a physical kicker has no Enabled, F 9)
   physics->AddCollider(host, false);
   // the hole carved in the playfield: the real cutout's radius at the playfield top (part of the cup data), not the VPX kicker radius
   const double holeMm = kickerCupSys7HoleRadiusMm;
   CarvedHoles::Add({ (double)m_d.m_vCenter.x * U, (double)m_d.m_vCenter.y * U, holeMm * 1e-3, (double)height * U });
   // the cup's switch (cup data: close and open depth of the ball's bottom, TUNED, kickerCupSys7.h)
   const double zClose = (double)height - kickerCupSys7SwitchCloseDepthMm * mmU, zOpen = (double)height - kickerCupSys7SwitchOpenDepthMm * mmU;
   physics->m_noneZones.AddKicker(this, m_d.m_vCenter.x, m_d.m_vCenter.y, holeMm * mmU, zClose, zOpen);
   PLOGI << "Kicker " << name << ": physical (" << (seen.count("type") ? "ini" : "default") << "), cup " << cup << " (" << faces->tris.size() << " faces, box z " << faces->lo[2] << ".." << faces->hi[2]
         << " VPU), centre (" << m_d.m_vCenter.x << ", " << m_d.m_vCenter.y << ") VPU, surface z " << height << " VPU, orientation " << m_d.m_orientation << " deg; hole radius " << holeMm
         << " mm (VPX radius " << m_d.m_radius * U * 1e3 << " mm, drawing only); switch closes at a bottom depth of " << kickerCupSys7SwitchCloseDepthMm << " mm (z " << zClose << " VPU), opens above "
         << kickerCupSys7SwitchOpenDepthMm << " mm (z " << zOpen << " VPU), rest " << kickerCupSys7RestDepthMm << " mm; material e_mpa " << rs.eMpa
         << ", mu " << rs.mu << ", restitution " << (rs.restitution > 0.0 ? rs.restitution : kE) << ", edge_radius_mm " << rs.edgeRadiusMm << "; scatter " << m_d.m_scatter << " deg (table difficulty " << m_ptable->m_globalDifficulty << "); Enabled ignored";
   return true;
}

// The virtual kicker (KickerDesign.txt 1g-1i, 4b, 6c; Requirements F 9-10), on the unified engine with [Features] kickers = 1 and type virtual (the default): a ball store. No KickerHitCircle,
// no playfield hole, no m_lockedInKicker. A sensing zone, the stock hit circle's (radius with the legacy factor, height from the surface to the surface + HitHeight), tests the ball's centre on
// each step's chord as a None trigger's zone does (ContactSystem::NoneZones). An empty, enabled kicker parks the entering ball where it is at the end of that step (Park: the capture,
// speed and spin to 0, energy out, booked as a loss, 1i); from then on it is a held ball (KickerHold, design 1h): infinite mass in the shared pair law, not stepped. Hit is fired after the step.
// A ball entering a full kicker is not caught (not listed, no Hit, as stock). Kick (4b): today's kick, angle and inclination in degrees; Unhit when the kicked ball's centre leaves the zone.
void Kicker::VirtualSetup(PhysicsEngine *const physics, const bool fromIni)
{
   const float height = m_ptable->GetSurfaceHeight(m_d.m_szSurface, m_d.m_vCenter.x, m_d.m_vCenter.y);
   const float radius = m_d.m_radius * (m_d.m_legacyMode ? (m_d.m_fallThrough ? 0.75f : 0.6f) : 1.f);   // (the stock hit circle's radius: only the inner circle of a legacy kicker starts a hit)
   physics->m_noneZones.AddKicker(this, m_d.m_vCenter.x, m_d.m_vCenter.y, radius, height, height + m_d.m_hit_height, true);
   KickerHold::Add(this);
   m_virtual = true;
   PLOGI << "Kicker " << GetName() << ": virtual (" << (fromIni ? "its ini" : "the default; type = physical in [Kicker:" + GetName() + "] makes it physical") << "): sensing zone centre ("
         << m_d.m_vCenter.x << ", " << m_d.m_vCenter.y << ") VPU, radius " << radius << " VPU (" << radius * ContactSystem::kMetresPerUnit * 1e3 << " mm; VPX radius " << m_d.m_radius
         << (m_d.m_legacyMode ? " x legacy factor" : "") << "), ball centre height " << height << " .. " << height + m_d.m_hit_height << " VPU; legacy " << m_d.m_legacyMode << ", fall-through "
         << m_d.m_fallThrough << " (not used: the ball is parked), enabled " << m_d.m_enabled << "; a parked ball is held (infinite mass, not stepped)";
}

bool Kicker::Park(HitBall *const ball, const Vertex3Ds &at)
{
   if (!m_virtual || !m_d.m_enabled || m_heldBall != nullptr)
      return false;   // full, or disabled (as stock's hit circle): not caught
   m_heldBall = m_lastHeldBall = ball;
   m_parkPlace = at;
   m_parkMovedNoted = false;
   Capture(ball);
   return true;
}

namespace
{
   // a ball's mechanical energy, J: kinetic (motion and spin) plus gravity's potential relative to ref (VPX units); the table's gravity (GetGravity), as the arm's work in KickXYZ
   double BallEnergyJ(const HitBall *const b, const Vertex3Ds &ref)
   {
      constexpr double U = ContactSystem::kMetresPerUnit, VS = ContactSystem::kVelocityScale;
      const double m = b->m_d.m_mass * ContactSystem::kMassUnitKg;
      const Vertex3Ds &v = b->m_d.m_vel, &L = b->m_angularmomentum, &p = b->m_d.m_pos;
      const double ke = 0.5 * m * VS * VS * ((double)v.x * v.x + (double)v.y * v.y + (double)v.z * v.z)
         + 0.5 * ContactSystem::kMassUnitKg * VS * VS * ((double)L.x * L.x + (double)L.y * L.y + (double)L.z * L.z) / b->Inertia();
      const Vertex3Ds &g = g_pplayer->m_physics->GetGravity();   // U/VPT^2
      const double gSI = U / (0.01 * 0.01);
      const double pe = -m * gSI * U * ((double)g.x * ((double)p.x - ref.x) + (double)g.y * ((double)p.y - ref.y) + (double)g.z * ((double)p.z - ref.z));
      return ke + pe;
   }
   constexpr double kResetLogJ = 10e-6;    // a reset that gives or takes more than this is logged on its own line (a hit: 0.5 m/s into a parked ball = 10000 uJ); smaller ones (a resting ball's
                                           // share of gravity in a cycle, 0.3 - 1.5 uJ measured, K4 run A) only add to the totals
}

// The capture (design 1i): the ball is parked where it is at the end of the step in which its centre entered the zone (no move back to the crossing: no snap), its speed and spin to 0, once.
// The energy it had (kinetic and spin) leaves: booked between cycles as a loss (ContactSystem::NoteDriveWork, negative). Its surface contact ends (a held ball has none). Nothing moves, so
// nothing is a placement. From then on nothing moves it: it is a held ball (KickerHold, design 1h). The Hit keeps the chord's crossing time (NoneZones::Chord).
void Kicker::Capture(HitBall *const b)
{
   const double e = BallEnergyJ(b, m_parkPlace);   // (the place at rest has 0)
   const char *const name = TestDriver::BallName(b);
   PLOGI << "Kicker " << GetName() << ": parks ball " << (name != nullptr ? name : "") << " where it is at the end of the step its centre entered the zone, (" << m_parkPlace.x << ", " << m_parkPlace.y
         << ", " << m_parkPlace.z << ") (" << (b->m_d.m_pos - m_parkPlace).Length() * ContactSystem::kMetresPerUnit * 1e3 << " mm from the ball: 0 expected), speed " << b->m_d.m_vel.Length() * ContactSystem::kVelocityScale
         << " m/s and spin " << b->m_angularmomentum.Length() / b->Inertia() / 0.01 << " rad/s set to 0, t=" << (g_pplayer != nullptr ? g_pplayer->m_time_msec : 0u);
   ContactSystem::Release(b, "kicker (held ball)");
   b->m_d.m_vel.SetZero();
   b->m_angularmomentum.SetZero();
   ContactSystem::DropKeptDoubles(b);   // the kicker wrote its state: the group's kept copy of the moving ball no longer describes it
   ContactSystem::NoteDriveWork(b, -e, "capture (virtual kicker, speed and spin to 0)");
}

void Kicker::CheckParked()
{
   const HitBall *const b = m_heldBall;
   if (m_parkMovedNoted || (b->m_d.m_pos == m_parkPlace && b->m_d.m_vel.IsZero() && b->m_angularmomentum.IsZero()))
      return;
   m_parkMovedNoted = true;
   KickerNotice("Kicker " + GetName() + ": its parked ball moved (" + std::to_string((b->m_d.m_pos - m_parkPlace).Length() * ContactSystem::kMetresPerUnit * 1e3) + " mm from its place, speed "
                + std::to_string(b->m_d.m_vel.Length() * ContactSystem::kVelocityScale) + " m/s): nothing should move a held ball (a fault; said once)");
}

// The arm stroke (physics/kickerarm.h; design 4c). Each stroke: the carried ball, its circle and how far along it the ball is. The path is prescribed (an idealised arm: the ball is carried,
// not riding, design 4c); the carried ball is a held ball (design 1h): not stepped; within a cycle the shared pair law moves it with the path's velocity and acceleration at the cycle's
// start (GroupStep::Ball::held), and at the cycle's end it is put on the path at that time (KickerHold::CycleEnd).
namespace
{
   struct ArmStroke
   {
      Kicker *kicker;
      HitBall *ball;
      Vertex3Ds rest;         // the ball's centre where the stroke began, U
      double dx, dy;          // the kick direction, unit (VPX x, y)
      double L;               // the circle's radius, U; its centre is L ahead of rest along (dx, dy), at rest's height
      double phi;             // how far along the circle, rad: 0 = rest (moving straight up), pi/2 = the top (moving along (dx, dy))
      double phiEnd;          // where the arm lets go, rad
      double speed;           // along the path, U/VPT (the kick's speed, held from the first cycle)
      uint32_t t0;            // the kick's physics time, ms (log)
   };
   std::vector<ArmStroke> g_strokes;   // (at most one per physical kicker)
   std::vector<Kicker *> g_virtualKickers;   // KickerHold: the virtual kickers, whose parked balls are held

   // the stroke's place and velocity at phi (U, U/VPT)
   void ArmPath(const ArmStroke &k, const double phi, Vertex3Ds &pos, Vertex3Ds &vel)
   {
      const double c = std::cos(phi), s = std::sin(phi);
      vel = Vertex3Ds((float)(k.speed * s * k.dx), (float)(k.speed * s * k.dy), (float)(k.speed * c));
      pos = Vertex3Ds((float)(k.rest.x + k.L * (1.0 - c) * k.dx), (float)(k.rest.y + k.L * (1.0 - c) * k.dy), (float)(k.rest.z + k.L * s));
   }

   // the path's acceleration at phi, U/VPT^2: the circle's centripetal pull, speed^2 / L towards the pivot (the speed along the path is constant)
   Vertex3Ds ArmAccel(const ArmStroke &k, const double phi)
   {
      const double c = std::cos(phi), s = std::sin(phi), a = k.speed * k.speed / k.L;
      return Vertex3Ds((float)(a * c * k.dx), (float)(a * c * k.dy), (float)(-a * s));
   }

   // The carried ball at the cycle's end (design 1h, 4c): on the path at this time, with the path's velocity and no spin. Within the cycle the pair law moved it with the path's velocity and
   // acceleration at the cycle's start (a held ball: GroupStep::Ball::held), or nothing moved it (no ball near): the place here differs from the pair law's by the path's third-order term
   // (at 1 m/s and L = 17 mm, under 1 um along the path). Nothing is booked: the arm's own work was booked at the kick, its work on other balls by the pair law (the group's drive work).
   // At the release point the ball is let go: it flies the rest of the cycle under gravity (no drive work after the release point) and is an ordinary ball from then on.
   bool ArmCycleEnd(ArmStroke &k, const double cycleS)
   {
      HitBall *const ball = k.ball;
      constexpr double U = ContactSystem::kMetresPerUnit;
      const double dtime = cycleS / 0.01;   // VPT
      const double ds = k.speed * dtime, left = (k.phiEnd - k.phi) * k.L;   // U along the circle: this cycle's, and what is left to the release point
      const bool letGo = ds >= left;
      k.phi = letGo ? k.phiEnd : k.phi + ds / k.L;
      Vertex3Ds pos, vel;
      ArmPath(k, k.phi, pos, vel);
      if (letGo)
      {
         const char *const name = TestDriver::BallName(ball);
         const double c = std::cos(k.phi), s = std::sin(k.phi);
         PLOGI << "Kicker " << k.kicker->GetName() << ": arm lets go of ball " << (name != nullptr ? name : "") << " at (" << pos.x << ", " << pos.y << ", " << pos.z << ") = " << k.L * (1.0 - c) * U * 1e3
               << " mm forward, " << k.L * s * U * 1e3 << " mm up from rest; velocity (" << vel.x << ", " << vel.y << ", " << vel.z << ") U/VPT = " << k.speed * ContactSystem::kVelocityScale << " m/s at "
               << std::atan2(c, s) * (180.0 / M_PI) << " deg up; carried " << (g_pplayer != nullptr ? g_pplayer->m_time_msec - k.t0 : 0u) << " cycles (arc " << k.L * k.phiEnd * U * 1e3 << " mm), t="
               << (g_pplayer != nullptr ? g_pplayer->m_time_msec : 0u);
         const float tFree = (float)((ds - left) / k.speed);   // VPT of this cycle after the release: free flight under gravity
         const Vertex3Ds &gU = g_pplayer->m_physics->GetGravity();
         pos += vel * tFree + gU * (0.5f * tFree * tFree);
         vel += gU * tFree;
      }
      static double worstUm = 0.0;   // (log: the largest put-on-path correction so far, against where the pair law or nothing left it)
      const double corrUm = (ball->m_d.m_pos - pos).Length() * U * 1e6;
      if (!letGo && corrUm > worstUm + 0.01)
      {
         worstUm = corrUm;
         PLOGI << "Kicker " << k.kicker->GetName() << ": carried ball put on the arm's path, largest correction so far " << corrUm << " um, t=" << (g_pplayer != nullptr ? g_pplayer->m_time_msec : 0u);
      }
      ball->m_d.m_pos = pos;
      ball->m_d.m_vel = vel;
      ball->m_angularmomentum.SetZero();
      ContactSystem::DropKeptDoubles(ball);   // the arm wrote its state (on the path, or let go): the group's kept copy no longer describes it
      return letGo;
   }
}

void KickerArm::Start(Kicker *const kicker, HitBall *const ball, const double dirX, const double dirY, const double radiusU, const double releaseRad, const double speed)
{
   BallGone(ball);
   g_strokes.push_back({ kicker, ball, ball->m_d.m_pos, dirX, dirY, radiusU, 0.0, releaseRad, speed, g_pplayer != nullptr ? g_pplayer->m_time_msec : 0u });
}

bool KickerArm::Carried(const HitBall *const ball)
{
   for (const ArmStroke &k : g_strokes)
      if (k.ball == ball)
         return true;
   return false;
}

void KickerArm::BallGone(const HitBall *const ball)
{
   g_strokes.erase(std::remove_if(g_strokes.begin(), g_strokes.end(), [ball](const ArmStroke &k) { return k.ball == ball; }), g_strokes.end());
}

void KickerArm::Remove(const Kicker *const kicker)
{
   g_strokes.erase(std::remove_if(g_strokes.begin(), g_strokes.end(), [kicker](const ArmStroke &k) { return k.kicker == kicker; }), g_strokes.end());
}

void KickerHold::Add(Kicker *const kicker)
{
   if (std::find(g_virtualKickers.begin(), g_virtualKickers.end(), kicker) == g_virtualKickers.end())
      g_virtualKickers.push_back(kicker);
}

void KickerHold::Remove(const Kicker *const kicker)
{
   std::erase(g_virtualKickers, kicker);
}

bool KickerHold::Held(const HitBall *const ball, Vertex3Ds *const accelU)
{
   if (accelU != nullptr)
      accelU->SetZero();
   for (const ArmStroke &k : g_strokes)
      if (k.ball == ball)
      {
         if (accelU != nullptr)
            *accelU = ArmAccel(k, k.phi);
         return true;
      }
   for (const Kicker *const k : g_virtualKickers)
      if (k->Parked(ball))
         return true;
   return false;
}

void KickerHold::CycleEnd(const double cycleS)
{
   for (Kicker *const k : g_virtualKickers)
      if (k->Parked(k->HeldBall()))
         k->CheckParked();
   for (size_t i = 0; i < g_strokes.size();)
      if (ArmCycleEnd(g_strokes[i], cycleS))
         g_strokes.erase(g_strokes.begin() + i);
      else
         i++;
}

void Kicker::HeldChange(HitBall *const ball, const bool in)
{
   if (in)
   {
      if (m_heldBall == nullptr)
         m_heldBall = ball;
      m_lastHeldBall = ball;
      return;
   }
   if (m_heldBall != ball)
      return;
   m_heldBall = nullptr;   // another ball still held takes its place (two balls in one cup: not designed for, design 5c)
   for (Ball *const b : g_pplayer->m_vball)
      if (&b->m_hitBall != ball && b->m_hitBall.m_d.m_vpVolObjs != nullptr && FindIndexOf(*b->m_hitBall.m_d.m_vpVolObjs, (IFireEvents *)this) >= 0)
      {
         m_heldBall = &b->m_hitBall;
         break;
      }
}

void Kicker::HeldEvent(const bool hit)
{
   TestDriver::TriggerEvent(GetName(), hit);   // (test driver: events.csv)
   FireGroupEvent(hit ? DISPID_HitEvents_Hit : DISPID_HitEvents_Unhit);
}

void Kicker::PhysicSetup(PhysicsEngine* physics, const bool isUI)
{
   if (!isUI && GetPartGroup() != nullptr && GetPartGroup()->GetReferenceSpace() != PartGroupData::SpaceReference::SR_PLAYFIELD)
      return;
   if (!isUI && ContactSystem::Unified() && PhysicalSetup(physics))
      return;

   const float height = m_ptable->GetSurfaceHeight(m_d.m_szSurface, m_d.m_vCenter.x, m_d.m_vCenter.y);
   const float radius = m_d.m_radius * (m_d.m_legacyMode ? (m_d.m_fallThrough ? 0.75f : 0.6f) // reduce the hit circle radius because only the inner circle of the kicker should start a hit event
                                                         : 1.f);
   if (isUI)
   {
      HitCircle *const phitcircle = new HitCircle(this, m_d.m_vCenter, radius, height, height + m_d.m_hit_height); // height of kicker hit cylinder
      phitcircle->m_ObjType = eKicker;
      phitcircle->m_obj = (IFireEvents *)this;
      physics->AddCollider(phitcircle, isUI);
   }
   else
   {
      KickerHitCircle *const phitcircle = new KickerHitCircle(this, m_d.m_vCenter, radius, height, height + m_d.m_hit_height); // height of kicker hit cylinder
      if (!m_d.m_legacyMode)
         ContactSurfaceHost::AddPlayfieldHole(m_d.m_vCenter.x, m_d.m_vCenter.y, m_d.m_radius);   // a non-legacy kicker is a hole in the playfield (a playfield on the force model has no floor in it)
      if (!m_d.m_legacyMode)
      {
         const float rad = phitcircle->radius * 0.8f;
         m_hitMesh.resize(kickerHitNumVertices);
         for (unsigned int t = 0; t < kickerHitNumVertices; t++)
         {
            // find the right normal by calculating the distance from current ball position to vertex of the kicker mesh
            Vertex3Ds vpos = Vertex3Ds(kickerHitMesh[t].x, kickerHitMesh[t].y, kickerHitMesh[t].z);
            vpos.x = vpos.x * rad + m_d.m_vCenter.x;
            vpos.y = vpos.y * rad + m_d.m_vCenter.y;
            vpos.z = vpos.z * rad + height;
            m_hitMesh[t] = vpos;
         }
      }
      phitcircle->m_enabled = m_d.m_enabled;
      phitcircle->m_ObjType = eKicker;
      phitcircle->m_obj = (IFireEvents *)this;
      phitcircle->m_pkicker = this;
      m_phitkickercircle = phitcircle;
      physics->AddCollider(phitcircle, isUI);
   }
}

void Kicker::PhysicRelease(PhysicsEngine* physics, const bool isUI)
{
   if (!isUI)
   {
      m_phitkickercircle = nullptr;
      m_hitMesh.clear();
      if (m_physical || m_virtual)
         physics->m_noneZones.Remove(this);   // (the cup host is a hit object of the engine, deleted with the others)
      KickerArm::Remove(this);
      KickerHold::Remove(this);
      m_physical = m_physicalRefused = m_virtual = m_parkMovedNoted = false;
      m_heldBall = m_lastHeldBall = nullptr;
   }
}

#pragma endregion



#pragma region Rendering

void Kicker::RenderSetup(Renderer *renderer)
{
   assert(m_renderer == nullptr);
   m_renderer = renderer;

   if (m_d.m_kickertype == KickerInvisible)
      return;

   m_baseHeight = m_ptable->GetSurfaceHeight(m_d.m_szSurface, m_d.m_vCenter.x, m_d.m_vCenter.y);

   // Ported at: VisualPinball.Engine/VPT/Kicker/KickerMeshGenerator.cs

   {
      Vertex3D_NoTex2 *buf = new Vertex3D_NoTex2[kickerPlateNumVertices];
      float rad = m_d.m_radius; 

      if (m_d.m_kickertype == KickerWilliams || m_d.m_kickertype == KickerGottlieb)
          rad = m_d.m_radius * 0.88f;
      else if (m_d.m_kickertype == KickerCup2)
          rad = m_d.m_radius * 0.87f;
      else if (m_d.m_kickertype != KickerCup && m_d.m_kickertype != KickerWilliams)
          rad = m_d.m_radius * 0.82f;

      for (unsigned int i = 0; i < kickerPlateNumVertices; i++)
      {
         buf[i].x = kickerPlate[i].x*rad + m_d.m_vCenter.x;
         buf[i].y = kickerPlate[i].y*rad + m_d.m_vCenter.y;
         buf[i].z = kickerPlate[i].z*rad + m_baseHeight;
         buf[i].nx = kickerPlate[i].nx;
         buf[i].ny = kickerPlate[i].ny;
         buf[i].nz = kickerPlate[i].nz;
         buf[i].tu = 0.0f;
         buf[i].tv = 0.0f;
      }

      std::shared_ptr<VertexBuffer> plateVertexBuffer = std::make_shared<VertexBuffer>(m_renderer->m_renderDevice, kickerPlateNumVertices, (float *)buf);
      std::shared_ptr<IndexBuffer> plateIndexBuffer = std::make_shared<IndexBuffer>(m_renderer->m_renderDevice, kickerPlateNumIndices, kickerPlateIndices);
      m_plateMeshBuffer = std::make_shared<MeshBuffer>(GetName() + ".Plate"s, plateVertexBuffer, plateIndexBuffer, true);

      delete[] buf;
   }

   const WORD * indices;
   switch (m_d.m_kickertype)
   {
       case KickerInvisible:
          assert(false);
       break;
       case KickerCup:
       {
          m_texture.reset(Texture::CreateFromFile(g_app->m_fileLocator.GetAppPath(FileLocator::AppSubFolder::Assets, "KickerCup.webp")));
          m_numIndices = kickerCupNumIndices;
          m_numVertices = kickerCupNumVertices;
          indices = kickerCupIndices;
       }
       break;
       case KickerWilliams:
       {
          m_texture.reset(Texture::CreateFromFile(g_app->m_fileLocator.GetAppPath(FileLocator::AppSubFolder::Assets, "KickerWilliams.webp")));
          m_numIndices = kickerWilliamsNumIndices;
          m_numVertices = kickerWilliamsNumVertices;
          indices = kickerWilliamsIndices;
       }
       break;
       case KickerGottlieb:
       {
          m_texture.reset(Texture::CreateFromFile(g_app->m_fileLocator.GetAppPath(FileLocator::AppSubFolder::Assets, "KickerGottlieb.webp")));
          m_numIndices = kickerGottliebNumIndices;
          m_numVertices = kickerGottliebNumVertices;
          indices = kickerGottliebIndices;
       }
       break;
       case KickerCup2:
       {
          m_texture.reset(Texture::CreateFromFile(g_app->m_fileLocator.GetAppPath(FileLocator::AppSubFolder::Assets, "KickerT1.webp")));
          m_numIndices = kickerT1NumIndices;
          m_numVertices = kickerT1NumVertices;
          indices = kickerT1Indices;
       }
       break;
       case KickerHole:
       {
          m_texture.reset(Texture::CreateFromFile(g_app->m_fileLocator.GetAppPath(FileLocator::AppSubFolder::Assets, "KickerHoleWood.webp")));
          m_numIndices = kickerHoleNumIndices;
          m_numVertices = kickerHoleNumVertices;
          indices = kickerHoleIndices;
       }
       break;
       default:
       case KickerHoleSimple:
       {
          m_texture.reset(Texture::CreateFromFile(g_app->m_fileLocator.GetAppPath(FileLocator::AppSubFolder::Assets, "KickerHoleWood.webp")));
          m_numIndices = kickerSimpleHoleNumIndices;
          m_numVertices = kickerSimpleHoleNumVertices;
          indices = kickerSimpleHoleIndices;
       }
       break;
   }

   //

   std::shared_ptr<VertexBuffer> vertexBuffer = std::make_shared<VertexBuffer>(m_renderer->m_renderDevice, m_numVertices);
   Vertex3D_NoTex2 *buf;
   vertexBuffer->Lock(buf);
   GenerateMesh(buf);
   vertexBuffer->Unlock();
   std::shared_ptr<IndexBuffer> indexBuffer = std::make_shared<IndexBuffer>(m_renderer->m_renderDevice, m_numIndices, indices);
   m_meshBuffer = std::make_shared<MeshBuffer>(GetName() + ".Kicker"s, vertexBuffer, indexBuffer, true);
}

void Kicker::RenderRelease()
{
   assert(m_renderer != nullptr);
   m_meshBuffer = nullptr;
   m_plateMeshBuffer = nullptr;
   m_meshEdgeBuffer = nullptr;
   m_plateMeshEdgeBuffer = nullptr;
   m_texture = nullptr;
   m_renderer = nullptr;
}

void Kicker::UpdateAnimation(const float diff_time_msec)
{
   assert(m_renderer != nullptr);
}

void Kicker::Render(const unsigned int renderMask)
{
   assert(m_renderer != nullptr);
   assert(!m_desktopBackdrop);
   const bool isStaticOnly = renderMask & Renderer::STATIC_ONLY;
   const bool isDynamicOnly = renderMask & Renderer::DYNAMIC_ONLY;
   const bool isReflectionPass = renderMask & Renderer::REFLECTION_PASS;
   const bool isUIPass = renderMask & Renderer::UI_EDGES || renderMask & Renderer::UI_FILL;
   TRACE_FUNCTION();

   if (isStaticOnly 
   || isReflectionPass
   || m_d.m_kickertype == KickerInvisible)
      return;

   const Vertex3Ds pos(m_d.m_vCenter.x, m_d.m_vCenter.y, m_baseHeight);
   if (isUIPass)
   {
      if (m_meshBuffer == nullptr || m_plateMeshBuffer == nullptr)
         return;
      if (renderMask & Renderer::UI_FILL)
      {
         m_renderer->m_renderDevice->DrawMesh(m_renderer->m_renderDevice->m_basicShader, true, pos, 0.f, m_plateMeshBuffer, RenderDevice::TRIANGLELIST, 0, kickerPlateNumIndices);
         m_renderer->m_renderDevice->DrawMesh(m_renderer->m_renderDevice->m_basicShader, true, pos, 0.f, m_meshBuffer, RenderDevice::TRIANGLELIST, 0, m_numIndices);
      }
      if (renderMask & Renderer::UI_EDGES)
      {
         if (m_meshEdgeBuffer == nullptr)
         {
            const WORD *indices;
            switch (m_d.m_kickertype)
            {
            case KickerInvisible: assert(false); break;
            case KickerCup: indices = kickerCupIndices; break;
            case KickerWilliams: indices = kickerWilliamsIndices; break;
            case KickerGottlieb: indices = kickerGottliebIndices; break;
            case KickerCup2: indices = kickerT1Indices; break;
            case KickerHole: indices = kickerHoleIndices; break;
            default:
            case KickerHoleSimple: indices = kickerSimpleHoleIndices; break;
            }
            vector<unsigned int> indices2(indices, indices + m_numIndices);
            m_meshEdgeBuffer = m_meshBuffer->CreateEdgeMeshBuffer(indices2);
         }
         m_renderer->m_renderDevice->DrawMesh(m_renderer->m_renderDevice->m_basicShader, false, pos, 0.f, m_meshEdgeBuffer, RenderDevice::LINELIST, 0, m_meshEdgeBuffer->m_ib->m_count);
         if (m_plateMeshEdgeBuffer == nullptr)
         {
            vector<unsigned int> indices(kickerPlateIndices, kickerPlateIndices + kickerPlateNumIndices);
            m_plateMeshEdgeBuffer = m_plateMeshBuffer->CreateEdgeMeshBuffer(indices);
         }
         m_renderer->m_renderDevice->DrawMesh(m_renderer->m_renderDevice->m_basicShader, false, pos, 0.f, m_plateMeshEdgeBuffer, RenderDevice::LINELIST, 0, m_plateMeshEdgeBuffer->m_ib->m_count);
      }
   }
   else if (m_d.m_kickertype != KickerInvisible)
   {
      m_renderer->m_renderDevice->ResetRenderState();
      if (m_d.m_kickertype == KickerHoleSimple)
         m_renderer->m_renderDevice->SetRenderState(RenderState::CULLMODE, RenderState::CULL_NONE);

      const Material *const mat = m_ptable->GetMaterial(m_d.m_szMaterial);
      m_renderer->m_renderDevice->m_basicShader->SetMaterial(mat);
      m_renderer->m_renderDevice->m_basicShader->SetTechniqueMaterial(ShaderTechnique::kickerBoolean, *mat);
      m_renderer->m_renderDevice->SetRenderState(RenderState::ZFUNC, RenderState::Z_ALWAYS);
      m_renderer->m_renderDevice->DrawMesh(m_renderer->m_renderDevice->m_basicShader, false, pos, 0.f, m_plateMeshBuffer, RenderDevice::TRIANGLELIST, 0, kickerPlateNumIndices);

      m_renderer->m_renderDevice->SetRenderState(RenderState::ZFUNC, RenderState::Z_LESSEQUAL);
      m_renderer->m_renderDevice->m_basicShader->SetBasic(mat, m_d.m_kickertype == KickerHoleSimple ? nullptr : m_texture.get());
      m_renderer->m_renderDevice->EnableAlphaBlend(false);
      m_renderer->m_renderDevice->DrawMesh(m_renderer->m_renderDevice->m_basicShader, false, pos, 0.f, m_meshBuffer, RenderDevice::TRIANGLELIST, 0, m_numIndices);
   }
}

#pragma endregion

void Kicker::ExportMesh(ObjLoader& loader)
{
   if (m_d.m_kickertype == KickerInvisible)
      return;

   const string& name = m_name;
   m_baseHeight = m_ptable->GetSurfaceHeight(m_d.m_szSurface, m_d.m_vCenter.x, m_d.m_vCenter.y);

   int num_vertices;
   int num_indices;
   const WORD *indices;
   switch (m_d.m_kickertype)
   {
       case KickerCup:
       {
           num_vertices = kickerCupNumVertices;
           num_indices = kickerCupNumIndices;
           indices = kickerCupIndices;
       }
       break;
       case KickerWilliams:
       {
           num_vertices = kickerWilliamsNumVertices;
           num_indices = kickerWilliamsNumIndices;
           indices = kickerWilliamsIndices;
       }
       break;
       case KickerGottlieb:
       {
           num_vertices = kickerGottliebNumVertices;
           num_indices = kickerGottliebNumIndices;
           indices = kickerGottliebIndices;
       }
       break;
       case KickerCup2:
       {
           num_vertices = kickerT1NumVertices;
           num_indices = kickerT1NumIndices;
           indices = kickerT1Indices;
       }
       break;
       case KickerHole:
       {
           num_vertices = kickerHoleNumVertices;
           num_indices = kickerHoleNumIndices;
           indices = kickerHoleIndices;
       }
       break;
       case KickerInvisible:
           assert(false);
       default:
       case KickerHoleSimple:
       {
           num_vertices = kickerSimpleHoleNumVertices;
           num_indices = kickerSimpleHoleNumIndices;
           indices = kickerSimpleHoleIndices;
       }
       break;
   }

   Vertex3D_NoTex2 *vertices = new Vertex3D_NoTex2[num_vertices];
   GenerateMesh(vertices);
   loader.WriteObjectName(name);
   loader.WriteVertexInfo(vertices, num_vertices);
   const Material * const mat = m_ptable->GetMaterial(m_d.m_szMaterial);
   loader.WriteMaterial(m_d.m_szMaterial, string(), mat);
   loader.UseTexture(m_d.m_szMaterial);
   loader.WriteFaceInfoList(indices, num_indices);
   loader.UpdateFaceOffset(num_vertices);
   delete[] vertices;
}

void Kicker::GenerateMesh(Vertex3D_NoTex2 *const buf) const
{
   int num_vertices;
   const Vertex3D_NoTex2 *vertices;
   float zoffset = 0.f;
   float zrot = m_d.m_orientation;
   switch (m_d.m_kickertype)
   {
       case KickerCup:
       {
           num_vertices = kickerCupNumVertices;
           vertices = kickerCupMesh;
           zoffset = -0.18f;
       }
       break;
       case KickerWilliams:
       {
           num_vertices = kickerWilliamsNumVertices;
           vertices = kickerWilliamsMesh;
           zrot = m_d.m_orientation + 90.f;
       }
       break;
       case KickerGottlieb:
       {
           num_vertices = kickerGottliebNumVertices;
           vertices = kickerGottliebMesh;
       }
       break;
       case KickerCup2:
       {
           num_vertices = kickerT1NumVertices;
           vertices = kickerT1Mesh;
       }
       break;
       case KickerHole:
       {
           num_vertices = kickerHoleNumVertices;
           vertices = kickerHoleMesh;
           zrot = 0.f;
       }
       break;
       case KickerInvisible:
           assert(false);
       default:
       case KickerHoleSimple:
       {
           num_vertices = kickerSimpleHoleNumVertices;
           vertices = kickerSimpleHoleMesh;
           zrot = 0.f;
       }
       break;
   }

   const Matrix3D fullMatrix = Matrix3D::MatrixRotateZ(ANGTORAD(zrot));

   for (int i = 0; i < num_vertices; i++)
   {
      Vertex3Ds vert = fullMatrix * Vertex3Ds{vertices[i].x, vertices[i].y, vertices[i].z + zoffset};

      buf[i].x = vert.x*m_d.m_radius + m_d.m_vCenter.x;
      buf[i].y = vert.y*m_d.m_radius + m_d.m_vCenter.y;
      buf[i].z = vert.z*m_d.m_radius + m_baseHeight;

      vert = fullMatrix.MultiplyVectorNoTranslate(Vertex3Ds{vertices[i].nx, vertices[i].ny, vertices[i].nz});
      buf[i].nx = vert.x;
      buf[i].ny = vert.y;
      buf[i].nz = vert.z;
      buf[i].tu = vertices[i].tu;
      buf[i].tv = vertices[i].tv;
   }
}

void Kicker::Translate(const Vertex2D &offset)
{
   m_d.m_vCenter.x += offset.x;
   m_d.m_vCenter.y += offset.y;
}

Vertex2D Kicker::GetCenter() const
{
   return m_d.m_vCenter;
}


void Kicker::Save(IObjectWriter& writer, const bool saveForUndo)
{
   writer.WriteVector2(FID(VCEN), m_d.m_vCenter);
   writer.WriteFloat(FID(RADI), m_d.m_radius);
   writer.WriteBool(FID(TMON), m_timerEnabled);
   writer.WriteInt(FID(TMIN), m_timerInterval);
   writer.WriteString(FID(MATR), m_d.m_szMaterial);
   writer.WriteString(FID(SURF), m_d.m_szSurface);
   writer.WriteBool(FID(EBLD), m_d.m_enabled);
   writer.WriteWideString(FID(NAME), MakeWString(m_name));
   writer.WriteInt(FID(TYPE), m_d.m_kickertype);
   writer.WriteFloat(FID(KSCT), m_d.m_scatter);
   writer.WriteFloat(FID(KHAC), m_d.m_hitAccuracy);
   writer.WriteFloat(FID(KHHI), m_d.m_hit_height);
   writer.WriteFloat(FID(KORI), m_d.m_orientation);
   writer.WriteBool(FID(FATH), m_d.m_fallThrough);
   writer.WriteBool(FID(LEMO), m_d.m_legacyMode);
   SaveSharedEditableFields(writer);
   writer.EndObject();
}

void Kicker::Load(IObjectReader& reader)
{
   SetDefaults(false);
   reader.AsObject(
      [this](int tag, IObjectReader& reader)
      {
         switch (tag)
         {
         case FID(PIID): reader.AsInt(); break;
         case FID(VCEN):
         {
            auto v = reader.AsVector2();
            m_d.m_vCenter.x = v.x;
            m_d.m_vCenter.y = v.y;
            break;
         }
         case FID(RADI): m_d.m_radius = reader.AsFloat(); break;
         case FID(KSCT): m_d.m_scatter = reader.AsFloat(); break;
         case FID(KHAC): m_d.m_hitAccuracy = reader.AsFloat(); break;
         case FID(KHHI): m_d.m_hit_height = reader.AsFloat(); break;
         case FID(KORI): m_d.m_orientation = reader.AsFloat(); break;
         case FID(MATR): m_d.m_szMaterial = reader.AsString(); break;
         case FID(TMON): m_timerEnabled = reader.AsBool(); break;
         case FID(TMIN): m_timerInterval = reader.AsInt(); break;
         case FID(EBLD): m_d.m_enabled = reader.AsBool(); break;
         case FID(TYPE):
         {
            m_d.m_kickertype = static_cast<KickerType>(reader.AsInt());
            // legacy handling:
            if (m_d.m_kickertype > KickerCup2)
               m_d.m_kickertype = KickerInvisible;
            break;
         }
         case FID(SURF): m_d.m_szSurface = reader.AsString(); break;
         case FID(NAME): m_name = MakeString(reader.AsWideString()); break;
         case FID(FATH): m_d.m_fallThrough = reader.AsBool(); break;
         case FID(LEMO): m_d.m_legacyMode = reader.AsBool(); break;
         default: LoadSharedEditableField(tag, reader); break;
         }
         return true;
      });
   m_phitkickercircle = nullptr;
}

STDMETHODIMP Kicker::InterfaceSupportsErrorInfo(REFIID riid)
{
   static const IID* arr[] =
   {
      &IID_IKicker,
   };

   for (size_t i = 0; i < std::size(arr); i++)
      if (InlineIsEqualGUID(*arr[i], riid))
         return S_OK;

   return S_FALSE;
}

STDMETHODIMP Kicker::CreateSizedBallWithMass(/*[in]*/ float radius, /*[in]*/ float mass, /*out, retval]*/ IBall **pResult)
{
   if (m_physical || m_virtual)
   {
      // a physical kicker: the ball is made at the hole's centre with its bottom at the surface, at rest, and falls into the cup (gravity the drive). A virtual kicker: the same place, and the
      // ball is parked there (listed, no Hit, as stock; a held ball, design 1h). Not when the kicker holds a ball or a ball overlaps that place (design 6c): no ball, said once
      const float height = m_ptable->GetSurfaceHeight(m_d.m_szSurface, m_d.m_vCenter.x, m_d.m_vCenter.y);
      const Vertex3Ds at { m_d.m_vCenter.x, m_d.m_vCenter.y, height + radius };
      bool occupied = m_heldBall != nullptr;
      for (const Ball *const b : g_pplayer->m_vball)
         occupied |= (b->m_hitBall.m_d.m_pos - at).Length() < b->m_hitBall.m_d.m_radius + radius;
      if (occupied)
      {
         if (!m_createNoted)
            KickerNotice("Kicker " + GetName() + ": CreateBall refused, the kicker holds a ball or a ball is in the way; no ball was made");
         m_createNoted = true;
         *pResult = nullptr;
         return S_OK;
      }
      Ball *const pball = g_pplayer->CreateBall(at.x, at.y, height, 0.f, 0.f, 0.f, radius, mass);
      if (m_virtual)
      {
         HitBall *const hb = &pball->m_hitBall;
         if (hb->m_d.m_vpVolObjs != nullptr)
            hb->m_d.m_vpVolObjs->push_back((IFireEvents *)this);   // listed (BallCntOver), no Hit
         m_heldBall = m_lastHeldBall = hb;
         m_parkPlace = hb->m_d.m_pos;
         m_parkMovedNoted = false;
         ContactSystem::DropKeptDoubles(hb);   // a new ball: a copy kept under its address belongs to a ball gone before
         PLOGI << "Kicker " << GetName() << ": CreateBall parks a new ball at (" << m_parkPlace.x << ", " << m_parkPlace.y << ", " << m_parkPlace.z << "), t=" << g_pplayer->m_time_msec;
      }
      pball->AddRef();
      *pResult = pball;
      return S_OK;
   }
   if (m_phitkickercircle)
   {
      const float height = m_ptable->GetSurfaceHeight(m_d.m_szSurface, m_d.m_vCenter.x, m_d.m_vCenter.y);

      Ball *const pball = g_pplayer->CreateBall(m_phitkickercircle->center.x, m_phitkickercircle->center.y, height, 0.1f, 0, 0, radius, mass);
      
      pball->AddRef();
      *pResult = pball;

      pball->m_hitBall.m_coll.m_hitflag = true; // HACK: avoid capture leaving kicker
      static constexpr Vertex3Ds hitnormal { FLT_MAX, FLT_MAX, FLT_MAX }; // unused due to newBall being true
      m_phitkickercircle->DoCollide(&(pball->m_hitBall), hitnormal, false, true);
   }

   return S_OK;
}

STDMETHODIMP Kicker::CreateSizedBall(/*[in]*/ float radius, /*out, retval]*/ IBall **pResult)
{
   return CreateSizedBallWithMass(radius, 1.f, pResult);
}

STDMETHODIMP Kicker::CreateBall(IBall **pResult)
{
   return CreateSizedBallWithMass(DEFAULT_BALL_SIZE, 1.f, pResult);
}


void Kicker::BallGone(const HitBall *const ball)
{
   if (m_phitkickercircle != nullptr && m_phitkickercircle->m_pHitBall == ball)
      m_phitkickercircle->m_pHitBall = nullptr;
   if (m_virtual && m_heldBall == ball)
   {
      PLOGI << "Kicker " << GetName() << ": parked ball destroyed";
      m_heldBall = nullptr;   // (no Unhit for a destroyed ball, as None zones)
   }
   else if (m_heldBall == ball)
      HeldChange(const_cast<HitBall *>(ball), false);   // (no Unhit for a destroyed ball, as None zones)
   KickerArm::BallGone(ball);
}

STDMETHODIMP Kicker::DestroyBall(int *pVal)
{
   int cnt = 0;

   if ((m_physical || m_virtual) && m_heldBall != nullptr)
   {
      ++cnt;
      g_pplayer->DestroyBall(m_heldBall->m_pBall);   // (BallGone moves the held state)
   }
   else if (m_phitkickercircle && m_phitkickercircle->m_pHitBall)
   {
      _ASSERTE(g_pplayer);
      ++cnt;
      HitBall *const b = m_phitkickercircle->m_pHitBall;
      m_phitkickercircle->m_pHitBall = nullptr;
      g_pplayer->DestroyBall(b->m_pBall);
   }

   if (pVal) *pVal = cnt;
   return S_OK;
}

// Ported at: VisualPinball.Unity/VisualPinball.Unity/VPT/Ball/BallManager.cs

STDMETHODIMP Kicker::KickXYZ(float angle, float speed, float inclination, float x, float y, float z)
{
   if (m_physicalRefused)
      return S_OK;   // a physical kicker that could not be built does nothing (design 6a)
   if (m_physical)
   {
      // Design 4c, 3c (approved by the user 2026-11-13; replaces 4a's velocity set): the arm stroke. The held ball is handed to an idealised kicker arm (KickerArm, below): its centre is
      // carried along a quarter circle of radius L = the cup's resting depth + kArmClearMm, in the vertical plane of the kick direction (angle + scatter, degrees), its pivot L ahead of the
      // resting centre at the centre's height, from rest (moving straight up) to where the path's direction equals the inclination (degrees: 0 = the full quarter, the ball's bottom then
      // kArmClearMm above the playfield, moving flat; 90 = at once, straight up), at the kick's speed from the first cycle on, then let go with that velocity and no spin. The arm's work
      // (kinetic energy at release less the kinetic energy before, plus the lift against gravity to the release point) is the drive's, reported to the energy book now (3c); the ball is
      // noted as placed now and again when it is let go (its first-sight rules, 3d). No held ball (the switch open), a stroke already under way, or speed 0: nothing (6d).
      constexpr double kArmClearMm = 6.35;   // the arm's top: the ball's bottom 0.25 in above the playfield (design 4c; chosen for feel by the user 2026-11-13, NOT measured)
      if ((x != 0.f || y != 0.f || z != 0.f) && !m_offsetNoted)
      {
         m_offsetNoted = true;
         KickerNotice("Kicker " + GetName() + ": KickXYZ's x, y, z offset and KickZ's height are ignored on a physical kicker (the ball is kicked where it rests); remove them from the script");
      }
      if (g_pplayer == nullptr || m_heldBall == nullptr || KickerArm::Carried(m_heldBall) || speed <= 0.f)
         return S_OK;
      float inclDeg = inclination;
      if (inclDeg < 0.f || inclDeg > 90.f)
      {
         // the quarter circle's direction runs from straight up (90) to flat (0) only; outside that the design gives no path (question in runs/2026-11-12-kicker-cup/K3c/results.txt)
         if (!m_inclNoted)
            KickerNotice("Kicker " + GetName() + ": inclination " + std::to_string(inclination) + " deg is outside the arm's 0..90 deg on a physical kicker; " + (inclDeg < 0.f ? "0" : "90") + " is used");
         m_inclNoted = true;
         inclDeg = inclDeg < 0.f ? 0.f : 90.f;
      }
      HitBall *const hb = m_heldBall;
      BallControl &ballControl = g_pplayer->m_liveUI->m_ballControl;   // (as today: the debug ball control follows the kicked ball)
      if (ballControl.GetSelectedBall() == nullptr)
         ballControl.SetDraggedBall(hb->m_pBall);
      if (ballControl.GetSelectedBall() == hb->m_pBall)
         ballControl.ReleaseDragTarget();
      // kinetic energy, motion and spin, J (VPX mass, U and 10 ms units -> kg, m, s)
      const double before = ContactSystem::kMassUnitKg * ContactSystem::kVelocityScale * ContactSystem::kVelocityScale
         * (0.5 * hb->m_d.m_mass * (double)hb->m_d.m_vel.Dot(hb->m_d.m_vel) + 0.5 * (double)hb->m_angularmomentum.Dot(hb->m_angularmomentum) / hb->Inertia());
      float anglerad = ANGTORAD(angle);   // yaw angle, zero is along -Y
      float scatterAngle = (m_d.m_scatter < 0.0f) ? c_hardScatter : ANGTORAD(m_d.m_scatter);   // if < 0 use the global value (as today)
      scatterAngle *= m_ptable->m_globalDifficulty;
      if (scatterAngle > 1.0e-5f)
      {
         float scatter = rand_mt_m11();   // -1..1, shaped quadratic and scaled (as today)
         scatter *= (1.0f - scatter * scatter) * 2.59808f * scatterAngle;
         anglerad += scatter;
      }
      constexpr double U = ContactSystem::kMetresPerUnit;
      const double dirX = std::sin((double)anglerad), dirY = -std::cos((double)anglerad);   // the kick direction, unit, VPX x y
      const double L = (kickerCupSys7RestDepthMm + kArmClearMm) * 1e-3 / U;                    // the arm's radius, U (Sys7: 17.0 mm)
      const double phiEnd = (double)ANGTORAD(90.0f - inclDeg);                                   // the release point along the quarter circle, rad
      hb->m_angularmomentum.SetZero();
      hb->m_coll.m_hitdistance = 0.0f;   // (the stock hit search's state, cleared as today)
      hb->m_coll.m_hittime = -1.0f;
      hb->m_coll.m_hitnormal.SetZero();
      hb->m_coll.m_hitvel.SetZero();
      hb->m_coll.m_hitflag = false;
      hb->m_coll.m_isContact = false;
      hb->m_coll.m_hitmoment_bit = true;
      hb->m_d.m_vel = Vertex3Ds(0.f, 0.f, speed);   // the path's start: straight up at the kick's speed (with 90 deg: let go here)
      ContactSystem::DropKeptDoubles(hb);   // the kick wrote its state: the group's kept copy of the resting ball no longer describes it
#ifdef C_DYNAMIC
      hb->m_dynamic = C_DYNAMIC;
#endif
      // the arm's work, J: 1/2 m v^2 at release less the kinetic energy before, plus m g h, h the rise to the release point along gravity (the lift normal to the playfield, g_n, and the
      // playfield slope's share along the forward move)
      const Vertex3Ds &gU = g_pplayer->m_physics->GetGravity();   // U/VPT^2
      const double gSI = U / (0.01 * 0.01), mKg = hb->m_d.m_mass * ContactSystem::kMassUnitKg, vSI = (double)speed * ContactSystem::kVelocityScale;
      const double fwdM = L * (1.0 - std::cos(phiEnd)) * U, upM = L * std::sin(phiEnd) * U;
      const double keJ = 0.5 * mKg * vSI * vSI, liftJ = -mKg * (double)gU.z * gSI * upM, slopeJ = -mKg * gSI * ((double)gU.x * dirX + (double)gU.y * dirY) * fwdM;
      PLOGI << "Kicker " << GetName() << ": arm stroke, ball " << (TestDriver::BallName(hb) != nullptr ? TestDriver::BallName(hb) : "") << " from (" << hb->m_d.m_pos.x << ", " << hb->m_d.m_pos.y << ", "
            << hb->m_d.m_pos.z << "), heading " << RADTOANG(anglerad) << " deg, " << vSI << " m/s, inclination " << inclDeg << " deg: L " << L * U * 1e3 << " mm, release " << fwdM * 1e3 << " mm forward / "
            << upM * 1e3 << " mm up, arc " << L * phiEnd * U * 1e3 << " mm (" << (vSI > 0.0 ? L * phiEnd * U / vSI * 1e3 : 0.0) << " ms); work: 1/2 m v^2 " << keJ * 1e6 << " uJ - before " << before * 1e6
            << " + m g_n h " << liftJ * 1e6 << " + slope " << slopeJ * 1e6 << " uJ, t=" << g_pplayer->m_time_msec;
      ContactSystem::NoteDriveWork(hb, keJ - before + liftJ + slopeJ, "kick (arm stroke)");
      if (phiEnd > 0.0)
      {
         KickerArm::Start(this, hb, dirX, dirY, L, phiEnd, (double)speed);
         ContactSystem::Release(hb, "kicker arm (held ball)");   // a carried ball is held (design 1h): no surface contacts while the arm carries it
      }
      ContactSystem::NotePlaced(hb);
      return S_OK;
   }
   if (m_virtual)
   {
      // Design 4b (user 2026-11-12): today's kick for compatibility, the parked ball released where it is held plus the script's offset (a teleport; a ball started inside something
      // falls under the first-sight rule, 3d), its velocity set with stock's speed split (an upward inclination takes its share of the speed, a downward one keeps the full horizontal speed)
      // and scatter; but angle and inclination in degrees only (no radians guess). The kick's work (energy after less before, the offset's height included) is the drive's, booked now (3c).
      // The ball stays listed: Unhit when its centre leaves the zone. No parked ball: nothing (6d).
      if (g_pplayer == nullptr || m_heldBall == nullptr)
         return S_OK;
      HitBall *const hb = m_heldBall;
      BallControl &ballControl = g_pplayer->m_liveUI->m_ballControl;   // (as today: the debug ball control follows the kicked ball)
      if (ballControl.GetSelectedBall() == nullptr)
         ballControl.SetDraggedBall(hb->m_pBall);
      if (ballControl.GetSelectedBall() == hb->m_pBall)
         ballControl.ReleaseDragTarget();
      const Vertex3Ds from = hb->m_d.m_pos;
      const double before = BallEnergyJ(hb, from);
      float anglerad = ANGTORAD(angle);                  // yaw angle, zero is along -Y
      const float inclrad = ANGTORAD(inclination);       // degrees only (design 4b)
      float scatterAngle = (m_d.m_scatter < 0.0f) ? c_hardScatter : ANGTORAD(m_d.m_scatter);   // if < 0 use the global value (as today)
      scatterAngle *= m_ptable->m_globalDifficulty;
      if (scatterAngle > 1.0e-5f)
      {
         float scatter = rand_mt_m11();   // -1..1, shaped quadratic and scaled (as today)
         scatter *= (1.0f - scatter * scatter) * 2.59808f * scatterAngle;
         anglerad += scatter;
      }
      const float speedz = sinf(inclrad) * speed;
      float speedxy = speed;
      if (speedz > 0.0f)
         speedxy *= cosf(inclrad);
      hb->m_angularmomentum.SetZero();
      hb->m_coll.m_hitdistance = 0.0f;   // (the stock hit search's state, cleared as today)
      hb->m_coll.m_hittime = -1.0f;
      hb->m_coll.m_hitnormal.SetZero();
      hb->m_coll.m_hitvel.SetZero();
      hb->m_coll.m_hitflag = false;
      hb->m_coll.m_isContact = false;
      hb->m_coll.m_hitmoment_bit = true;
      hb->m_d.m_pos.x += x;
      hb->m_d.m_pos.y += y;
      hb->m_d.m_pos.z += z;
      hb->m_d.m_vel = Vertex3Ds(sinf(anglerad) * speedxy, -cosf(anglerad) * speedxy, speedz);
      ContactSystem::DropKeptDoubles(hb);   // the kick wrote its state: the group's kept copy of the parked ball no longer describes it
#ifdef C_DYNAMIC
      hb->m_dynamic = C_DYNAMIC;
#endif
      m_heldBall = nullptr;   // no longer parked: an ordinary ball from here
      const double work = BallEnergyJ(hb, from) - before;
      const char *const name = TestDriver::BallName(hb);
      PLOGI << "Kicker " << GetName() << ": virtual kick of ball " << (name != nullptr ? name : "") << " from (" << from.x << ", " << from.y << ", " << from.z << "), offset (" << x << ", " << y << ", "
            << z << "), heading " << RADTOANG(anglerad) << " deg, inclination " << inclination << " deg: velocity (" << hb->m_d.m_vel.x << ", " << hb->m_d.m_vel.y << ", " << hb->m_d.m_vel.z << ") U/VPT = "
            << hb->m_d.m_vel.Length() * ContactSystem::kVelocityScale << " m/s; t=" << g_pplayer->m_time_msec;
      ContactSystem::NoteDriveWork(hb, work, "kick (virtual kicker)");
      ContactSystem::NotePlaced(hb);
      return S_OK;
   }
   if (g_pplayer && m_phitkickercircle && m_phitkickercircle->m_pHitBall)
   {
      BallControl& ballControl = g_pplayer->m_liveUI->m_ballControl;
      Ball* const kickedBall = m_phitkickercircle->m_pHitBall->m_pBall;
      if (ballControl.GetSelectedBall() == nullptr)
      {
         // Ball control most recently kicked if none currently.
         ballControl.SetDraggedBall(kickedBall);
      }
      if (ballControl.GetSelectedBall() == kickedBall)
      {
         // Clear any existing ball control target to allow kickout to work correctly.
         ballControl.ReleaseDragTarget();
      }
      float anglerad = ANGTORAD(angle);					// yaw angle, zero is along -Y axis

      if (fabsf(inclination) > (float)(M_PI / 2.0))		// radians or degrees?  if greater PI/2 assume degrees
         inclination *= (float)(M_PI / 180.0);			// convert to radians

      float scatterAngle = (m_d.m_scatter < 0.0f) ? c_hardScatter : ANGTORAD(m_d.m_scatter); // if < 0 use global value
      scatterAngle *= m_ptable->m_globalDifficulty;		// apply dificulty weighting

      if (scatterAngle > 1.0e-5f)						// ignore near zero angles
      {
         float scatter = rand_mt_m11();					// -1.0f..1.0f
         scatter *= (1.0f - scatter*scatter)*2.59808f * scatterAngle; // shape quadratic distribution and scale
         anglerad += scatter;
      }

      const float speedz = sinf(inclination) * speed;
      if (speedz > 0.0f)
         speed *= cosf(inclination);

      m_phitkickercircle->m_pHitBall->m_angularmomentum.SetZero();
      m_phitkickercircle->m_pHitBall->m_coll.m_hitdistance = 0.0f;
      m_phitkickercircle->m_pHitBall->m_coll.m_hittime = -1.0f;
      m_phitkickercircle->m_pHitBall->m_coll.m_hitnormal.SetZero();
      m_phitkickercircle->m_pHitBall->m_coll.m_hitvel.SetZero();
      m_phitkickercircle->m_pHitBall->m_coll.m_hitflag = false;
      m_phitkickercircle->m_pHitBall->m_coll.m_isContact = false;
      m_phitkickercircle->m_pHitBall->m_coll.m_hitmoment_bit = true;
      m_phitkickercircle->m_pHitBall->m_d.m_pos.x += x; // brian's suggestion
      m_phitkickercircle->m_pHitBall->m_d.m_pos.y += y;
      m_phitkickercircle->m_pHitBall->m_d.m_pos.z += z;
      m_phitkickercircle->m_pHitBall->m_d.m_vel.x =  sinf(anglerad) * speed;
      m_phitkickercircle->m_pHitBall->m_d.m_vel.y = -cosf(anglerad) * speed;
      m_phitkickercircle->m_pHitBall->m_d.m_vel.z = speedz;
      m_phitkickercircle->m_pHitBall->m_d.m_lockedInKicker = false;
#ifdef C_DYNAMIC
      m_phitkickercircle->m_pHitBall->m_dynamic = C_DYNAMIC;
#endif
      m_phitkickercircle->m_pHitBall = nullptr;
   }

   return S_OK;
}

STDMETHODIMP Kicker::KickZ(float angle, float speed, float inclination, float heightz)
{
   KickXYZ(angle, speed, inclination, 0, 0, heightz);
   return S_OK;
}

STDMETHODIMP Kicker::Kick(float angle, float speed, float inclination)
{
   KickXYZ(angle, speed, inclination, 0, 0, 0);
   return S_OK;
}

STDMETHODIMP Kicker::get_X(float *pVal)
{
   *pVal = m_d.m_vCenter.x;
   return S_OK;
}

STDMETHODIMP Kicker::put_X(float newVal)
{
   m_d.m_vCenter.x = newVal;

   return S_OK;
}

STDMETHODIMP Kicker::get_Y(float *pVal)
{
   *pVal = m_d.m_vCenter.y;
   return S_OK;
}

STDMETHODIMP Kicker::put_Y(float newVal)
{
   m_d.m_vCenter.y = newVal;

   return S_OK;
}

STDMETHODIMP Kicker::get_Surface(BSTR *pVal)
{
   *pVal = MakeWideBSTR(m_d.m_szSurface);
   return S_OK;
}

STDMETHODIMP Kicker::put_Surface(BSTR newVal)
{
   m_d.m_szSurface = MakeString(newVal);

   return S_OK;
}

STDMETHODIMP Kicker::get_Enabled(VARIANT_BOOL *pVal)
{
   *pVal = FTOVB(m_d.m_enabled);
   return S_OK;
}

STDMETHODIMP Kicker::put_Enabled(VARIANT_BOOL newVal)
{
   m_d.m_enabled = VBTOb(newVal);   // (kept for the script to read back, e.g. a saucer class's HasBall)
   if (m_phitkickercircle)
      m_phitkickercircle->m_enabled = m_d.m_enabled;
   if ((m_physical || m_physicalRefused) && !m_enabledNoted)
   {
      // a physical kicker has no Enabled (Requirements F 9): the cup and its switch stay; said once per table start
      m_enabledNoted = true;
      KickerNotice("Kicker " + GetName() + ": Enabled is ignored on a physical kicker (the cup and its switch are always there); its value is only kept for the script to read");
   }

   return S_OK;
}

STDMETHODIMP Kicker::get_Scatter(float *pVal)
{
   *pVal = m_d.m_scatter;
   return S_OK;
}

STDMETHODIMP Kicker::put_Scatter(float newVal)
{
   m_d.m_scatter = newVal;

   return S_OK;
}

STDMETHODIMP Kicker::get_HitAccuracy(float *pVal)
{
   *pVal = m_d.m_hitAccuracy;
   return S_OK;
}

STDMETHODIMP Kicker::put_HitAccuracy(float newVal)
{
   m_d.m_hitAccuracy = saturate(newVal);

   return S_OK;
}

STDMETHODIMP Kicker::get_HitHeight(float *pVal)
{
   *pVal = m_d.m_hit_height;
   return S_OK;
}

STDMETHODIMP Kicker::put_HitHeight(float newVal)
{
   m_d.m_hit_height = newVal;

   return S_OK;
}

STDMETHODIMP Kicker::get_Orientation(float *pVal)
{
   *pVal = m_d.m_orientation;
   return S_OK;
}

STDMETHODIMP Kicker::put_Orientation(float newVal)
{
   m_d.m_orientation = newVal;

   return S_OK;
}

STDMETHODIMP Kicker::get_Radius(float *pVal)
{
   *pVal = m_d.m_radius;
   return S_OK;
}

STDMETHODIMP Kicker::put_Radius(float newVal)
{
   m_d.m_radius = newVal;

   return S_OK;
}

STDMETHODIMP Kicker::get_FallThrough(VARIANT_BOOL *pVal)
{
   *pVal = FTOVB(m_d.m_fallThrough);
   return S_OK;
}

STDMETHODIMP Kicker::put_FallThrough(VARIANT_BOOL newVal)
{
   m_d.m_fallThrough = VBTOb(newVal);

   return S_OK;
}

STDMETHODIMP Kicker::get_Legacy(VARIANT_BOOL *pVal)
{
   *pVal = FTOVB(m_d.m_legacyMode);
   return S_OK;
}

STDMETHODIMP Kicker::put_Legacy(VARIANT_BOOL newVal)
{
   m_d.m_legacyMode = VBTOb(newVal);

   return S_OK;
}


STDMETHODIMP Kicker::get_DrawStyle(KickerType *pVal)
{
   *pVal = m_d.m_kickertype;
   return S_OK;
}

STDMETHODIMP Kicker::put_DrawStyle(KickerType newVal)
{
   m_d.m_kickertype = newVal;
   // legacy handling:
   if (m_d.m_kickertype > KickerCup2)
	   m_d.m_kickertype = KickerInvisible;

   return S_OK;
}

STDMETHODIMP Kicker::get_Material(BSTR *pVal)
{
   *pVal = MakeWideBSTR(m_d.m_szMaterial);
   return S_OK;
}

STDMETHODIMP Kicker::put_Material(BSTR newVal)
{
   m_d.m_szMaterial = MakeString(newVal);

   return S_OK;
}

STDMETHODIMP Kicker::BallCntOver(int *pVal)
{
   int cnt = 0;

   if (g_pplayer)
   {
      for (const auto pball : g_pplayer->m_vball)
      {
         if (pball->m_hitBall.m_d.m_vpVolObjs && FindIndexOf(*(pball->m_hitBall.m_d.m_vpVolObjs), (IFireEvents*)this) >= 0) // cast to IFireEvents necessary, as it is stored like this in HitObject.m_obj
         {
            ++cnt;
            g_pplayer->m_pactiveball = pball; // set active ball for scriptor
         }
      }
   }

   *pVal = cnt;
   return S_OK;
}

STDMETHODIMP Kicker::get_LastCapturedBall(IBall **pVal)
{
    if (pVal && g_pplayer && (m_physical || m_virtual))   // a physical kicker: the last ball its held test held; a virtual kicker: the last it parked; if it still exists
    {
        for (const auto ball : g_pplayer->m_vball)
            if (m_lastHeldBall != nullptr && &ball->m_hitBall == m_lastHeldBall)
            {
                ball->QueryInterface(IID_IBall, (void **)pVal);
                return S_OK;
            }
        ShowError(m_lastHeldBall == nullptr ? "LastCapturedBall was called but no ball was captured!" : "LastCapturedBall was called but ball is already destroyed!");
        return E_POINTER;
    }
    if (!pVal || !g_pplayer || !m_phitkickercircle)
    {
        return E_POINTER;
    }

    if (!m_phitkickercircle->m_lastCapturedBall)
    {
        ShowError("LastCapturedBall was called but no ball was captured!");
        return E_POINTER;
    }

    bool ballFound = false;
    for (const auto ball : g_pplayer->m_vball)
    {
        if (ball == m_phitkickercircle->m_lastCapturedBall->m_pBall)
        {
            ballFound = true;
            break;
        }
    }

    if (!ballFound)
    {
        ShowError("LastCapturedBall was called but ball is already destroyed!");
        return E_POINTER;
    }

    Ball * const pBall = m_phitkickercircle->m_lastCapturedBall->m_pBall;

    if (!pBall)
    {
        return E_POINTER;
    }
    pBall->QueryInterface(IID_IBall, (void **)pVal);

    return S_OK;
}

float KickerHitCircle::HitTest(const BallS& ball, const float dtime, CollisionEvent& coll) const
{
   return HitTestBasicRadius(ball, dtime, coll, false, false, false); //any face, not-lateral, non-rigid
}

// Ported at: VisualPinball.Unity/VisualPinball.Unity/VPT/Kicker/KickerCollider.cs

void KickerHitCircle::DoChangeBallVelocity(HitBall *const pball, const Vertex3Ds &hitnormal) const
{
    float minDist_sqr = FLT_MAX;
    unsigned int idx = ~0u;
    for (size_t t = 0; t < m_pkicker->m_hitMesh.size(); t++)
    {
        // find the right normal by calculating the distance from current ball position to vertex of the kicker mesh               
        const float length_sqr = (pball->m_d.m_pos - m_pkicker->m_hitMesh[t]).LengthSquared();
        if (length_sqr < minDist_sqr)
        {
            minDist_sqr = length_sqr;
            idx = (unsigned int)t;
        }
    }
    //minDist_sqr = sqrtf(minDist_sqr);

    if (idx != ~0u)
    {
        // we have the nearest vertex now use the normal and damp it so it doesn't speed up the ball velocity too much
        const Vertex3Ds hitnorm(kickerHitMesh[idx].nx, kickerHitMesh[idx].ny, kickerHitMesh[idx].nz);
        Vertex3Ds surfVel, tangent, surfP;
        const float dot = -pball->m_d.m_vel.Dot(hitnorm);
        const float reactionImpulse = pball->m_d.m_mass * fabsf(dot);
        /*
        if (pball->m_pos.z > pball->m_radius)
        {
        // ball is on a surface(e.g. upper playfield) use a high friction and a different calculation to compensate surface collision
        friction = 1.0f;
        surfP = -pball->m_radius * hitnorm;      // surface contact point relative to center of mass

        surfVel = pball->SurfaceVelocity(surfP); // velocity at impact point

        tangent = surfVel - surfVel.Dot(hitnorm) * hitnorm; // calc the tangential velocity
        }
        else
        */
        surfP = -pball->m_d.m_radius * hitnormal;    // surface contact point relative to center of mass

        surfVel = pball->SurfaceVelocity(surfP);         // velocity at impact point

        tangent = surfVel - surfVel.Dot(hitnormal) * hitnorm; // calc the tangential velocity

        pball->m_d.m_vel += dot * hitnorm; // apply collision impulse (along normal, so no torque)
#ifdef C_DYNAMIC
        pball->m_dynamic = C_DYNAMIC;
#endif

        constexpr float friction = 0.3f;
        const float tangentSpSq = tangent.LengthSquared();

        if (tangentSpSq > 1e-6f)
        {
            tangent /= sqrtf(tangentSpSq);           // normalize to get tangent direction
            const float vt = surfVel.Dot(tangent);   // get speed in tangential direction

            // compute friction impulse
            const Vertex3Ds cross = CrossProduct(surfP, tangent);
            const float kt = 1.0f/pball->m_d.m_mass + tangent.Dot(CrossProduct(cross / pball->Inertia(), surfP));

            // friction impulse can't be greater than coefficient of friction times collision impulse (Coulomb friction cone)
            const float maxFric = friction * reactionImpulse;
            const float jt = clamp(-vt / kt, -maxFric, maxFric);

            pball->ApplySurfaceImpulse(jt * cross, jt * tangent);
        }
    }
}

void KickerHitCircle::DoCollide(HitBall *const pball, const Vertex3Ds &hitnormal, const bool hitbit, const bool newBall)
{
   if (m_pHitBall) return;                              // a previous ball already in kicker

   const int i = FindIndexOf(*(pball->m_d.m_vpVolObjs), m_obj); // check if kicker in ball's volume set

   if (newBall || ((!hitbit) == (i < 0)))            // New or (Hit && !Vol || UnHit && Vol)
   {
      if (m_pkicker->m_d.m_legacyMode || newBall)
      {
         pball->m_d.m_pos += STATICTIME * pball->m_d.m_vel; // move ball slightly forward
      }

      if (i < 0) // entering Kickers volume
      {
         const float grabHeight = (m_hitBBox.zlow + pball->m_d.m_radius) * m_pkicker->m_d.m_hitAccuracy;
         // early out here if the ball is slow and we are near the kicker center
         const bool hitEvent = (pball->m_d.m_pos.z < grabHeight || m_pkicker->m_d.m_legacyMode || newBall);

         if(!hitEvent)
         {
            DoChangeBallVelocity(pball, hitnormal);

            // this is an ugly hack to prevent the ball stopping rapidly at the kicker bevel
            // something with the friction calculation is wrong in the physics engine
            // so we monitor the ball velocity if it drop under a length value of 0.2
            // if so we take the last "good" velocity to help the ball moving over the critical spot at the kicker bevel
            // this hack seems to work only if the kicker is on the playfield, a kicker attached to a wall has still problems
            // because the friction calculation for a wall is also different
            if (pball->m_d.m_vel.LengthSquared() < (float)(0.2*0.2))
               pball->m_d.m_vel = pball->m_oldVel;

            pball->m_oldVel = pball->m_d.m_vel;
         }

         if (hitEvent)
         {
            pball->m_d.m_lockedInKicker = !m_pkicker->m_d.m_fallThrough;
            if (pball->m_d.m_lockedInKicker)
            {
               pball->m_d.m_vpVolObjs->push_back(m_obj);		// add kicker to ball's volume set
               m_pHitBall = pball;
               m_lastCapturedBall = pball;
               if (BallControl* const ballControl = m_physics->GetBallControl(); ballControl && pball->m_pBall == ballControl->GetSelectedBall())
                  ballControl->ReleaseDragTarget();
            }

            // Don't fire the hit event if the ball was just created
            // Fire the event before changing ball attributes, so scripters can get a useful ball state
            if (!newBall)
               m_pkicker->FireGroupEvent(DISPID_HitEvents_Hit);

            if (pball->m_d.m_lockedInKicker || m_pkicker->m_d.m_fallThrough)	// script may have unfrozen the ball
            {
               // if ball falls through hole, we fake the collision algo by changing the ball height
               // in HitTestBasicRadius() the z-position of the ball is checked if it is >= to the hit cylinder
               // if we don't change the height of the ball we get a lot of hit events while the ball is falling!!

               // Only mess with variables if ball was not kicked during event
               pball->m_d.m_vel.SetZero();
               pball->m_angularmomentum.SetZero();
               pball->m_d.m_pos.x = center.x;
               pball->m_d.m_pos.y = center.y;
#ifdef C_DYNAMIC
               pball->m_dynamic = 0;
#endif
               if (m_pkicker->m_d.m_fallThrough)
                  pball->m_d.m_pos.z = m_hitBBox.zlow - pball->m_d.m_radius - 5.0f;
               else
                  pball->m_d.m_pos.z = m_hitBBox.zlow + pball->m_d.m_radius/**pball->m_radius/radius*/;
            }
            else
               m_pHitBall = nullptr; // make sure
         }
      }
      else // exiting kickers volume
      {
         pball->m_d.m_vpVolObjs->erase(pball->m_d.m_vpVolObjs->begin() + i); // remove kicker to ball's volume set
         m_pkicker->FireGroupEvent(DISPID_HitEvents_Unhit);
      }
   }
}
