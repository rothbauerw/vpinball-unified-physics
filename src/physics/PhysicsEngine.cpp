// license:GPLv3+

#include "core/stdafx.h"
#include <array>

#include <iomanip>
#include <optional>

#include "PhysicsEngine.h"

#include "hitflipper.h"
#include "hitplunger.h"
#include "contactsystem.h"
#include "contactprof.h"

#include "input/PlungerHandler.h"
#include "physics/cabinet/NudgeHandler.h"
#include "plugins/MsgPlugin.h"
#include "plugins/VPXPlugin.h"
#include "core/VPXPluginAPIImpl.h"
#include "ui/live/LiveUI.h"

#include "utils/ushock_output.h"

#include "physics/hitrubber.h"
#include "physics/physicsini.h"
#include "physics/magnet.h"
#include "physics/turntable.h"
#include "physics/testdriver.h"
#include "parts/ball.h"
#include "parts/flipper.h"
#include "parts/primitive.h"
static unsigned s_ctReplays = 0, s_ctReplayUs = 0;   // VPX_CONTACT_TIMING: overlap replays and their time since the last cost-log row
static uint64_t s_ctMagUs = 0;   // VPX_CONTACT_TIMING: time of Magnets::Step since the last cost-log row, us

PhysicsEngine::PhysicsEngine(PinTable *const table)
   : m_plumbHandler(table->GetSettings())
   , m_table(table)
   , m_hitPlayfield(table)
   , m_hitTopGlass(table)
   , m_hitoctree(this)
   , m_hitoctree_dynamic(this)
{
   ContactSystem::SetTable(table->m_filename);   // unified or stock for this table: from its physics ini, before any part's PhysicSetup asks
   m_physicsMaxLoops = table->m_PhysicsMaxLoops == 0xFFFFFFFFu ? 0 : table->m_PhysicsMaxLoops * (10000 / PHYSICS_STEPTIME) /*2*/;
   m_contacts.reserve(8);

   // a host the previous table did not release is tested against every ball of this one (slowdown that grows with each restart in one session)
   m_leftoverHosts = ContactHost::LeftoverHosts();
   if (!m_leftoverHosts.empty())
      PLOGW << "Physics: surface hosts left over from the previous table (a leak, the game will slow down with each restart): " << m_leftoverHosts;
   else
      PLOGI << "Physics: no surface hosts left over from the previous table";
   HitRubber::ClearPlayfieldHoles();   // the kickers add theirs in PhysicSetup, below

   // Collect all hit objects
   const FRect3D tableBounds = table->GetBoundingBox();
   m_hitoctree.SetBounds(FRect(tableBounds.left, tableBounds.right, tableBounds.top, tableBounds.bottom)); // Limit to table bounds as we don't expect to play outside of it
   m_pendingHitObjects = &m_hitoctree.BeginReset();
   m_pendingHitObjects->clear();
   for (IEditable *const pe : table->GetParts())
      if (IHitable *const ph = pe->GetIHitable(); ph)
         ph->PhysicSetup(this, false);

   AddCabinetBoundingHitShapes(table);
   for (HitObject *const pho : *m_pendingHitObjects)
      RegisterHitObject(pho);

   PLOGI << "Initializing octree"; // For profiling
   m_hitoctree.EndReset();
   m_hitoctree.Finalize();
   m_pendingHitObjects = nullptr;
   #if !defined(NDEBUG) && defined(PRINT_DEBUG_COLLISION_TREE)
      m_hitoctree.DumpTree(0);
   #endif

   m_onUpdatePhysicsMsgId = g_pplayer ? g_pplayer->m_pluginAPI.GetMsgID(VPXPI_NAMESPACE, VPXPI_EVT_ON_UPDATE_PHYSICS) : 0;

#ifdef DEBUGPHYSICS
   c_hitcnts = 0;
   c_collisioncnt = 0;
   c_contactcnt = 0;
   #ifdef C_DYNAMIC
   c_staticcnt = 0;
   #endif
   c_embedcnts = 0;
   c_timesearch = 0;

   c_traversed = 0;
   c_tested = 0;
   c_deepTested = 0;
#endif
}

PhysicsEngine::~PhysicsEngine()
{
   delete m_UIQuadTtree;
   delete m_playfieldHost;
   for (ContactHost *const h : m_ownedHosts)
      delete h;

   if (m_pendingHitObjects)
      ReleaseVHO(*m_pendingHitObjects, false);
   ReleaseVHO(m_hitoctree.GetHitObjects(), false);

   if (g_pplayer)
      g_pplayer->m_pluginAPI.ReleaseMsgID(m_onUpdatePhysicsMsgId);

   // We should release objects from the dynamic tree except HitBall (but there are only HitBall...)
}

void PhysicsEngine::ReleaseVHO(const vector<HitObject *> &vho, bool isUI)
{
   vector<IEditable *> editables;
   for (size_t i = 0; i < vho.size(); i++)
   {
      if (FindIndexOf(editables, vho[i]->m_editable) == -1)
      {
         editables.push_back(vho[i]->m_editable);
         if (auto ph = vho[i]->m_editable->GetIHitable(); ph)
            ph->PhysicRelease(this, isUI);
      }
      if (vho[i]->GetType() != eBall) // As balls own their HitBall hit object
         delete vho[i];
   }
}

void PhysicsEngine::SetGravity(float slopeDeg, float strength)
{
   m_gravity.x = 0;
   m_gravity.y = sinf(ANGTORAD(slopeDeg)) * strength;
   m_gravity.z = -cosf(ANGTORAD(slopeDeg)) * strength;
}

void PhysicsEngine::Update(IEditable *editable)
{
   // Editable parts without physics (timers, light sequencers, part groups,...)
   if (editable->GetIHitable() == nullptr)
      return;

   // Modifying physics suspends the simulation, allowing for interactive edit with a deferred quadtree rebuild
   if (g_pplayer)
      g_pplayer->SetPlayState(false);

   // The UI quadtree is updated immediately as it supports interactive picking
   GetUIQuadTree()->Update(editable);

   // The gameplay colliders are also updated immediately (they may be displayed in the editor)
   // but the static quadtree rebuild is deferred as it is a lengthy operation (see FlushStaticQuadTree)
   if (editable->GetItemType() != eItemBall) // Balls manage their own HitBall in the dynamic quadtree
   {
      ReleaseStaticColliders(editable);
      AddStaticColliders(editable);
      m_staticQuadTreeDirty = true;
   }
}

void PhysicsEngine::Add(IEditable *editable)
{
   assert(editable->GetIHitable() != nullptr);

   // Modifying physics suspends the simulation, allowing for interactive edit with a deferred quadtree rebuild
   if (g_pplayer)
      g_pplayer->SetPlayState(false);

   AddStaticColliders(editable);
   m_staticQuadTreeDirty = true;

   if (m_UIQuadTtree)
      m_UIQuadTtree->AddEditable(editable);
}

void PhysicsEngine::Remove(IEditable* editable)
{
   assert(editable->GetIHitable() != nullptr);

   // Modifying physics suspends the simulation, allowing for interactive edit with a deferred quadtree rebuild
   if (g_pplayer)
      g_pplayer->SetPlayState(false);

   ReleaseStaticColliders(editable);
   m_staticQuadTreeDirty = true;

   if (m_UIQuadTtree)
      m_UIQuadTtree->Remove(editable);
}

void PhysicsEngine::AddStaticColliders(IEditable *editable)
{
   // Collect the editable's colliders and add them to the gameplay quadtree (balls insert and own their HitBall through AddCollider)
   vector<HitObject *> hitObjects;
   CollectColliders(editable, &hitObjects, false);
   for (HitObject *const pho : hitObjects)
      RegisterHitObject(pho);
   vector<HitObject *> &vho = m_hitoctree.BeginReset();
   vho.insert(vho.end(), hitObjects.begin(), hitObjects.end());
}

void PhysicsEngine::ReleaseStaticColliders(IEditable *editable)
{
   editable->GetIHitable()->PhysicRelease(this, false); // Balls remove their HitBall from the dynamic quadtree through RemoveCollider
   vector<HitObject *> &vho = m_hitoctree.BeginReset();
   std::erase_if(vho,
      [this, editable](HitObject *ho)
      {
         if (ho->m_editable == editable)
         {
            UnregisterHitObject(ho);
            delete ho;
            return true;
         }
         return false;
      });
}

void PhysicsEngine::RegisterHitObject(HitObject *hitObject)
{
   if (hitObject->GetType() == eFlipper)
      m_vFlippers.push_back(static_cast<HitFlipper *>(hitObject));
   else if (hitObject->GetType() == ePlunger)
      m_vPlungers.push_back(static_cast<HitPlunger *>(hitObject));
   if (MoverObject *const pmo = hitObject->GetMoverObject(); pmo && pmo->AddToList()) // Spinner, Gate, Flipper, Plunger (ball is added separately on each create ball)
      m_vmover.push_back(pmo);
}

void PhysicsEngine::UnregisterHitObject(HitObject *hitObject)
{
   if (hitObject->GetType() == eFlipper)
      RemoveFromVectorSingle(m_vFlippers, static_cast<HitFlipper *>(hitObject));
   else if (hitObject->GetType() == ePlunger)
      RemoveFromVectorSingle(m_vPlungers, static_cast<HitPlunger *>(hitObject));
   if (MoverObject *const pmo = hitObject->GetMoverObject(); pmo && pmo->AddToList())
      RemoveFromVectorSingle(m_vmover, pmo);
}

void PhysicsEngine::FlushStaticQuadTree()
{
   if (m_staticQuadTreeDirty)
   {
      m_staticQuadTreeDirty = false;
      m_hitoctree.EndReset();
   }
}

void PhysicsEngine::AddCollider(HitObject *collider, const bool isUI)
{
   assert(collider->m_editable != nullptr);
   collider->m_physics = this;
   if (MoverObject *const mover = collider->GetMoverObject())
      mover->m_physics = this;
   collider->CalcHitBBox();
   if (!isUI && (collider->GetType() == eBall))
   {
      HitBall *const hitBall = static_cast<HitBall *>(collider);
      m_vmover.push_back(&hitBall->m_mover); // balls are always added separately to this list!
      m_vballs.push_back(hitBall);
      m_hitoctree_dynamic.Insert(collider);
   }
   else
   {
      assert(m_pendingHitObjects != nullptr);
      m_pendingHitObjects->push_back(collider);
      //assert(collider->m_hitBBox.left >= -100000); // Invalid bound
   }
}

void PhysicsEngine::RemoveCollider(HitObject * collider, const bool isUI)
{
   assert(collider->GetType() == eBall); // Ball are the only object owning its HitObject, therefore the only one allowed to call here
   if (!isUI)
   {
      HitBall *const ball = static_cast<HitBall *>(collider);
      RemoveFromVectorSingle<MoverObject *>(m_vmover, &ball->m_mover);
      RemoveFromVectorSingle(m_vballs, ball);
      m_hitoctree_dynamic.Remove(collider);
      // the ball object may outlive its removal (a destroyed ball is freed later, it may still be rendered): a band must let go of it now, not when it is freed
      if (ball->m_bandOwner != nullptr || ball->m_contact != nullptr)
         HitRubberBandBallDeleted(ball);
   }
}

void PhysicsEngine::CollectColliders(IEditable *editable, vector<HitObject *> *hitObjects, bool isUI)
{
   assert(m_pendingHitObjects == nullptr);
   m_pendingHitObjects = hitObjects;
   if (auto ph = editable->GetIHitable(); ph)
      ph->PhysicSetup(this, isUI);
   m_pendingHitObjects = nullptr;
}

void PhysicsEngine::AddCabinetBoundingHitShapes(PinTable *const table)
{
   // TODO these are the only colliders that are added without an associated editable. We could enforce one here to guarantee m_editable != nullptr and simplify the code
   
   // simple outer borders:
   AddCollider(new LineSeg(table, Vertex2D(table->m_right, table->m_top), Vertex2D(table->m_right, table->m_bottom), 0.f, table->m_glassTopHeight), false);
   AddCollider(new LineSeg(table, Vertex2D(table->m_left, table->m_bottom), Vertex2D(table->m_left, table->m_top), 0.f, table->m_glassBottomHeight), false);
   AddCollider(new LineSeg(table, Vertex2D(table->m_right, table->m_bottom), Vertex2D(table->m_left, table->m_bottom), 0.f, table->m_glassBottomHeight),  false);
   AddCollider(new LineSeg(table, Vertex2D(table->m_left, table->m_top), Vertex2D(table->m_right, table->m_top), 0.f, table->m_glassTopHeight), false);

   // glass:
   Vertex3Ds * const rgv3D = new Vertex3Ds[4];
   rgv3D[0] = Vertex3Ds(table->m_left, table->m_top, table->m_glassTopHeight);
   rgv3D[1] = Vertex3Ds(table->m_right, table->m_top, table->m_glassTopHeight);
   rgv3D[2] = Vertex3Ds(table->m_right, table->m_bottom, table->m_glassBottomHeight);
   rgv3D[3] = Vertex3Ds(table->m_left, table->m_bottom, table->m_glassBottomHeight);
   AddCollider(new Hit3DPoly(table, rgv3D, 4), false);

   /*
   // playfield:
   Vertex3Ds * const rgv3D = new Vertex3Ds[4];
   rgv3D[3] = Vertex3Ds(table->m_left, table->m_top, 0.f);
   rgv3D[2] = Vertex3Ds(table->m_right, table->m_top, 0.f);
   rgv3D[1] = Vertex3Ds(table->m_right, table->m_bottom, 0.f);
   rgv3D[0] = Vertex3Ds(table->m_left, table->m_bottom, 0.f);
   Hit3DPoly * const ph3dpoly = new Hit3DPoly(table, rgv3D, 4); //!!
   ph3dpoly->SetFriction(table->m_overridePhysics ? table->m_fOverrideContactFriction : table->m_friction);
   ph3dpoly->m_elasticity = table->m_overridePhysics ? table->m_fOverrideElasticity : table->m_elasticity;
   ph3dpoly->m_elasticityFalloff = table->m_overridePhysics ? table->m_fOverrideElasticityFalloff : table->m_elasticityFalloff;
   ph3dpoly->m_scatter = ANGTORAD(table->m_overridePhysics ? table->m_fOverrideScatterAngle : table->m_scatter);
   AddCollider(ph3dpoly, false);
   */

   // playfield:
   m_hitPlayfield = HitPlane(table, Vertex3Ds(0, 0, 1), 0.f);
   m_hitPlayfield.m_physics = this; // copy assignment above does not preserve it
   m_hitPlayfield.SetFriction(table->m_overridePhysics ? table->m_fOverrideContactFriction : table->m_friction);
   m_hitPlayfield.m_elasticity = table->m_overridePhysics ? table->m_fOverrideElasticity : table->m_elasticity;
   m_hitPlayfield.m_elasticityFalloff = table->m_overridePhysics ? table->m_fOverrideElasticityFalloff : table->m_elasticityFalloff;
   m_hitPlayfield.m_scatter = ANGTORAD(table->m_overridePhysics ? table->m_fOverrideScatterAngle : table->m_scatter);

   // the playfield on the force model (vpx_physics.ini, [Features] playfield_contact = 1, material in [Playfield]): a surface host the plane's collide and contact hand the ball to
   if (m_playfieldHost == nullptr && g_pplayer != nullptr && g_pplayer->m_implicitPlayfieldMesh != nullptr) // (the engine can run without a player) an explicit playfield primitive is its own host (Primitive::PhysicSetup)
   {
      const std::filesystem::path ini = HitRubber::IniPathFor(table->m_filename);
      RubberObjectSettings ps;
      if (HitRubber::LoadPlayfieldSettings(ps, ini, (double)m_hitPlayfield.m_friction))
      {
         // as far as a ball can fly: stock's plane has no limit, and a ball outside these faces would be engaged and released every cycle with no force and fall through
         const float x0 = table->m_left - 100000.f, x1 = table->m_right + 100000.f, y0 = table->m_top - 100000.f, y1 = table->m_bottom + 100000.f;
         auto faces = std::make_shared<RubberSurface>();
         faces->Add(Vertex3Ds(x0, y0, 0.f), Vertex3Ds(x1, y0, 0.f), Vertex3Ds(x1, y1, 0.f));
         faces->Add(Vertex3Ds(x0, y0, 0.f), Vertex3Ds(x1, y1, 0.f), Vertex3Ds(x0, y1, 0.f));
         m_playfieldHost = new HitRubber(table, faces, ps, ini, "playfield", "playfield");
         m_playfieldHost->m_tableMu = (double)m_hitPlayfield.m_friction;
         m_playfieldHost->CalcHitBBox();
         m_hitPlayfield.m_surfaceHost = m_playfieldHost;
         // the engine steps a ball the host holds through the host's mover (a primitive's host gets it from its hit objects; this one has none, so add it here)
         if (MoverObject *const pmo = m_playfieldHost->GetMoverObject(); pmo != nullptr && pmo->AddToList())
            m_vmover.push_back(pmo);
         PLOGI << "Playfield: force model on (e_mpa " << ps.eMpa << ", alpha " << ps.alphaSPerM << ", tan_delta " << ps.tanDelta << ", mu " << ps.mu << "; UNVERIFIED material)";
      }
   }
   m_hitPlayfield.m_surfaceHost = m_playfieldHost;
   Turntables::Setup(this, table, m_vmover, m_ownedHosts);   // discs set into the playfield (the holes are read by the playfield host)
   {
      const std::filesystem::path ini = HitRubber::IniPathFor(table->m_filename);
      const PhysicsIni::Features feat = PhysicsIni::Active(ini);
      m_hitSpeedEvents = feat.hitSpeed;
      m_hitSpeedMin = std::max(0.f, feat.hitSpeedMin);
      if (m_hitSpeedEvents)
         PLOGI << "Physics: hit_speed on (Ball.HitSpeed, playfield Hit event above " << m_hitSpeedMin << ")";
   }

   // glass:
   Vertex3Ds glassNormal(0, table->m_bottom - table->m_top, table->m_glassBottomHeight - table->m_glassTopHeight);
   glassNormal.Normalize();
   m_hitTopGlass = HitPlane(table, Vertex3Ds(0, glassNormal.z, -glassNormal.y), -table->m_glassTopHeight);
   m_hitTopGlass.m_physics = this; // copy assignment above does not preserve it
   m_hitTopGlass.m_elasticity = 0.2f;
}

// Diagnostics: VPX_CONTACT_LOG=<file> writes one row per collision (C) and per resting contact (T) of the first ball: game time, kind, object type, object name,
// hit distance (U). Which objects touch the ball while it rests, and which of them fire every step.
static void LogBallContact(const char kind, const CollisionEvent &ce)
{
   static FILE *const f = []() -> FILE * {
      const char *const p = getenv("VPX_CONTACT_LOG");
      FILE *const fp = (p != nullptr && *p) ? fopen(p, "w") : nullptr;
      if (fp != nullptr)
         fprintf(fp, "t_ms,kind,obj_type,obj_name,hitdist_u,hitnormal_x,hitnormal_y,hitnormal_z,vx,vy,vz,Lx,Ly,Lz\n");
      return fp;
   }();
   if (f == nullptr || ce.m_obj == nullptr || g_pplayer->m_vball.empty() || ce.m_ball != &g_pplayer->m_vball[0]->m_hitBall)
      return;
   fprintf(f, "%u,%c,%d,%s,%.5f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", g_pplayer->m_time_msec, kind, (int)ce.m_obj->m_ObjType,
      ce.m_obj->m_editable != nullptr ? ce.m_obj->m_editable->GetName().c_str() : "-", ce.m_hitdistance, ce.m_hitnormal.x, ce.m_hitnormal.y, ce.m_hitnormal.z,
      ce.m_ball->m_d.m_vel.x, ce.m_ball->m_d.m_vel.y, ce.m_ball->m_d.m_vel.z, ce.m_ball->m_angularmomentum.x, ce.m_ball->m_angularmomentum.y, ce.m_ball->m_angularmomentum.z);
}

bool PhysicsEngine::RecordContact(const CollisionEvent& newColl)
{
   if (m_recordContacts) // remember all contacts?
   {
      #ifdef USE_EMBREE
      const std::lock_guard<std::mutex> lg(mtx); // multiple threads may end up here and call push_back
      #endif
      m_contacts.push_back(newColl);
      return true;
   }
   return false;
}

AsyncDynamicQuadTree *PhysicsEngine::GetUIQuadTree()
{
   if (m_UIQuadTtree == nullptr)
      m_UIQuadTtree = new AsyncDynamicQuadTree(this, m_table, true);
   return m_UIQuadTtree;
}

vector<HitObject *> PhysicsEngine::GetUIHitObjects(IEditable *editable)
{
   return GetUIQuadTree()->GetHitObjects(editable);
}

void PhysicsEngine::RayCast(const Vertex3Ds &source, const Vertex3Ds &target, const bool uiCast, vector<HitTestResult> &vhoHit)
{
   // Create a ray (ball) that travels in 3D space along the given ray, and find what it intersects with.
   HitBall ballT;
   ballT.m_d.m_pos = source;
   ballT.m_d.m_vel = target - source;
   ballT.m_d.m_radius = 0.01f;
   ballT.m_coll.m_hittime = 1.0f;
   ballT.CalcHitBBox(); // need to update here, as only done lazily

   if (uiCast)
   {
      GetUIQuadTree()->HitTestXRay(&ballT, vhoHit, ballT.m_coll);
   }
   else
   {
      FlushStaticQuadTree();
      m_hitoctree_dynamic.HitTestXRay(&ballT, vhoHit, ballT.m_coll);
      m_hitoctree.HitTestXRay(&ballT, vhoHit, ballT.m_coll);
   }

   // Sort result by distance from viewer
   std::ranges::sort(vhoHit.begin(), vhoHit.end(), [](const HitTestResult& a, const HitTestResult& b) { return a.m_time < b.m_time; });
}


FlipperMoverObject *PhysicsEngine::FindFlipperMover(const std::string& name) const
{
   auto lower = [](std::string t) { for (char& c : t) c = (char)tolower((unsigned char)c); return t; };
   const std::string want = lower(name);
   for (HitFlipper *const f : m_vFlippers)
      if (lower(f->m_flipperMover.m_pflipper->GetName()) == want)
         return &f->m_flipperMover;
   return nullptr;
}

bool PhysicsEngine::FlipperNear(const BallS &ball, const float margin) const
{
   for (const HitFlipper *const f : m_vFlippers)
      if (f->GapToBall(ball) < margin)
         return true;
   return false;
}

// A bat the ball could meet within `horizon` physics steps: the gap to the bat is less than `margin` plus what the nearest bat point and the ball travel in that time.
// Unlike FlipperNear this ignores a bat that is still or a ball beside the base (which hardly moves), so a rubber next to the base keeps its force model.
float PhysicsEngine::FlipperClearance(const BallS &ball, const float horizon) const
{
   const float ballReach = ball.m_vel.Length() * horizon;
   float clearance = FLT_MAX;
   for (const HitFlipper *const f : m_vFlippers)
      clearance = std::min(clearance, f->GapToBall(ball) - ballReach - f->ReachToBall(ball, horizon));
   return clearance;
}


// Diagnostics: VPX_PHYS_TIMING=<file> times the 1 ms physics iterations and groups them per RENDERED frame (PhysicsEngine::OnFinishFrame). One row is written for each
// frame that is slow, plus a summary row every 5 s. A frame is slow when the time since the previous frame exceeds VPX_PHYS_TIMING_GAP_US (default 25000: a dropped frame at
// 60 Hz), its physics (all iterations together) took longer than VPX_PHYS_TIMING_FRAME_US (default 3000), or any single iteration took longer than VPX_PHYS_TIMING_ITER_US
// (default 1000). The clock is read only inside real iterations, so measuring does not load the idle physics loop. Columns:
//   kind          F (slow frame) or S (5 s summary)
//   t_ms          game time
//   gap_us        time since the previous rendered frame (S: the longest in the window)
//   phys_us       time spent in physics iterations during this frame (S: the longest frame)
//   iters         1 ms iterations run this frame (S: frames in the window)
//   sim_sum_us / sim_max_us      PhysicsSimulateCycle, the collision physics itself (S: sum over the window / longest single iteration)
//   other_sum_us / other_max_us  the rest of each iteration: script timers, input, nudge, plunger, mover velocities (S: same)
//   balls, engaged               balls in play, balls currently driven by a flipper rubber (S: engaged = number of slow frames in the window)
namespace
{
   struct PhysTiming
   {
      FILE *f = nullptr;
      uint32_t frameThrUs = 3000, iterThrUs = 1000, gapThrUs = 25000;
      uint64_t lastFrameUs = 0;
      uint32_t nextSummaryMs = 5000;
      // this frame
      uint32_t iters = 0;
      uint64_t simSum = 0, simMax = 0, otherSum = 0, otherMax = 0;
      // summary window
      uint32_t frames = 0, slow = 0;
      uint64_t wGap = 0, wPhys = 0, wSimSum = 0, wSimMax = 0, wOtherSum = 0, wOtherMax = 0;

      PhysTiming()
      {
         const char *const p = getenv("VPX_PHYS_TIMING");
         if (p == nullptr || *p == '\0')
            return;
         f = fopen(p, "w");
         if (f == nullptr)
            return;
         auto env = [](const char *name, uint32_t def) { const char *const v = getenv(name); return (v != nullptr && *v) ? (uint32_t)strtoul(v, nullptr, 10) : def; };
         frameThrUs = env("VPX_PHYS_TIMING_FRAME_US", frameThrUs);
         iterThrUs = env("VPX_PHYS_TIMING_ITER_US", iterThrUs);
         gapThrUs = env("VPX_PHYS_TIMING_GAP_US", gapThrUs);
         fprintf(f, "kind,t_ms,gap_us,phys_us,iters,sim_sum_us,sim_max_us,other_sum_us,other_max_us,balls,engaged\n");
      }
      static PhysTiming &Get() { static PhysTiming t; return t; }
      bool On() const { return f != nullptr; }

      void AddIter(const uint64_t simUs, const uint64_t otherUs)
      {
         iters++;
         simSum += simUs;
         otherSum += otherUs;
         simMax = std::max(simMax, simUs);
         otherMax = std::max(otherMax, otherUs);
      }

      void EndFrame()
      {
         const uint64_t now = usec();
         const uint64_t gap = lastFrameUs ? now - lastFrameUs : 0;
         lastFrameUs = now;
         const uint64_t phys = simSum + otherSum;
         uint32_t engaged = 0;
         for (const Ball *const b : g_pplayer->m_vball)
            engaged += (b->m_hitBall.m_rubberOwner != nullptr) ? 1u : 0u;
         const bool isSlow = gap > gapThrUs || phys > frameThrUs || simMax > iterThrUs || otherMax > iterThrUs;
         frames++;
         slow += isSlow ? 1u : 0u;
         wGap = std::max(wGap, gap);
         wPhys = std::max(wPhys, phys);
         wSimSum += simSum;
         wOtherSum += otherSum;
         wSimMax = std::max(wSimMax, simMax);
         wOtherMax = std::max(wOtherMax, otherMax);
         if (isSlow)
            fprintf(f, "F,%u,%llu,%llu,%u,%llu,%llu,%llu,%llu,%u,%u\n", g_pplayer->m_time_msec, (unsigned long long)gap, (unsigned long long)phys, iters, (unsigned long long)simSum,
               (unsigned long long)simMax, (unsigned long long)otherSum, (unsigned long long)otherMax, (uint32_t)g_pplayer->m_vball.size(), engaged);
         iters = 0;
         simSum = simMax = otherSum = otherMax = 0;
         if (g_pplayer->m_time_msec >= nextSummaryMs)
         {
            fprintf(f, "S,%u,%llu,%llu,%u,%llu,%llu,%llu,%llu,%u,%u\n", g_pplayer->m_time_msec, (unsigned long long)wGap, (unsigned long long)wPhys, frames, (unsigned long long)wSimSum,
               (unsigned long long)wSimMax, (unsigned long long)wOtherSum, (unsigned long long)wOtherMax, (uint32_t)g_pplayer->m_vball.size(), slow);
            fflush(f);
            nextSummaryMs = g_pplayer->m_time_msec + 5000;
            frames = slow = 0;
            wGap = wPhys = wSimSum = wSimMax = wOtherSum = wOtherMax = 0;
         }
      }
   };
}

void PhysicsEngine::ResetPerFrameStats()
{
   // Reset per frame debug counters
   m_phys_iterations = 0;

#ifdef DEBUGPHYSICS
   c_hitcnts = 0;
   c_collisioncnt = 0;
   c_contactcnt = 0;
   #ifdef C_DYNAMIC
   c_staticcnt = 0;
   #endif
   c_embedcnts = 0;
   c_timesearch = 0;

   c_traversed = 0;
   c_tested = 0;
   c_deepTested = 0;
#endif
}

void PhysicsEngine::OnFinishFrame()
{
   m_lastFlipTime = usec();
   if (PhysTiming::Get().On())
      PhysTiming::Get().EndFrame();
}

void PhysicsEngine::StartPhysics()
{
   m_startTime_usec = usec();
   m_curPhysicsFrameTime = m_startTime_usec;
   m_nextPhysicsFrameTime = m_curPhysicsFrameTime + PHYSICS_STEPTIME;
   m_physicsStarted = true;

#ifdef PLAYBACK
   if (m_playback)
   {
      float physicsStepTime;
      ParseLog((LARGE_INTEGER*)&physicsStepTime, (LARGE_INTEGER*)&m_startTime_usec);
   }
#endif

#ifdef LOG
   PLOGD.printf("Step Time %llu", m_startTime_usec);
   PLOGD.printf("End Frame");
#endif
}

void PhysicsEngine::UpdatePhysics(uint64_t targetTimeUs)
{
   // Simulation can only run once it has been started (timestamps are set by StartPhysics)
   if (!m_physicsStarted)
      return;

   if (g_pplayer)
      g_pplayer->m_logicProfiler.EnterProfileSection(FrameProfiler::PROFILE_PHYSICS);
   uint64_t initial_time_usec = targetTimeUs;
   PhysTiming &pt = PhysTiming::Get();

   // DJRobX's crazy latency-reduction code
   uint64_t delta_frame = 0;
   #if !defined(ENABLE_BGFX)
   if (g_pplayer && g_pplayer->m_minphyslooptime > 0 && m_lastFlipTime > 0)
   {
      // We want the physics loops to sync up to the frames, not
      // the post-render period, as that can cause some judder.
      delta_frame = initial_time_usec - m_lastFlipTime;
      initial_time_usec -= delta_frame;
   }
   #endif

   // When paused or after debugging, shift whole game forward in time
   // TODO not sure why we would need noTimeCorrect, as pause should already have shifted the timings
   if (g_pplayer && (!g_pplayer->IsPlaying() || g_pplayer->m_noTimeCorrect))
   {
      const uint64_t curPhysicsFrameTime = m_startTime_usec + (uint64_t)(m_time_sec * 1000000.0);
      if (initial_time_usec > curPhysicsFrameTime)
      {
         // m_time_sec can lag m_curPhysicsFrameTime by up to one step: cap the shift so the frame boundary never overshoots the target
         const uint64_t timeShift = min(initial_time_usec - curPhysicsFrameTime, initial_time_usec > m_curPhysicsFrameTime ? initial_time_usec - m_curPhysicsFrameTime : (uint64_t)0);
         m_startTime_usec += timeShift;
         m_nextPhysicsFrameTime += timeShift;
         m_curPhysicsFrameTime += timeShift;
         g_pplayer->m_noTimeCorrect = false;
      }
   }

   // Walk a single physics step forward
   if (g_pplayer && g_pplayer->m_step)
   {
      m_curPhysicsFrameTime -= PHYSICS_STEPTIME;
      g_pplayer->m_step = false;
   }

#ifdef LOG
   const double timepassed = (double)(initial_time_usec - m_curPhysicsFrameTime) / 1000000.0;
   const float frametime =
   #ifdef PLAYBACK
      (!m_playback) ? (float)(timepassed * 100.0) : ParseLog((LARGE_INTEGER*)&initial_time_usec, (LARGE_INTEGER*)&m_nextPhysicsFrameTime);
   #else
      #define TIMECORRECT 1
      #ifdef TIMECORRECT
         (float)(timepassed * 100.0);
         // 1.456927f;
      #else
         0.45f;
      #endif
   #endif //PLAYBACK
   PLOGD.printf("Frame Time %.20f %u %u %u %u", frametime, initial_time_usec >> 32, initial_time_usec, m_nextPhysicsFrameTime >> 32, m_nextPhysicsFrameTime);
   PLOGD.printf("End Frame");
#endif

   if (g_pplayer && m_nextPhysicsFrameTime < initial_time_usec)
      g_pplayer->m_pluginAPI.BroadcastVPXMsg(m_onUpdatePhysicsMsgId, nullptr);

   // Rebuild the static quadtree if physics was edited while the simulation was suspended (live edit)
   FlushStaticQuadTree();

   while (m_nextPhysicsFrameTime < initial_time_usec) // loop here until physics (=simulated) time catches up to current real time, still staying behind real time by up to one physics emulation step
   {
      TestDriver::LoopStart();   // cycle timing (docs/TestMatrixDesign.txt section 3 point 9); nothing unless VPX_TEST_INPUT is set
      const uint64_t ptIter0 = pt.On() ? usec() : 0;
      m_time_sec = max(m_time_sec, (double)(m_curPhysicsFrameTime - m_startTime_usec) / 1000000.0); // First iteration is done before precise time
      m_time_msec = (uint32_t)((m_curPhysicsFrameTime - m_startTime_usec) / 1000); // Get time in milliseconds for timers
      if (g_pplayer)
      {
         g_pplayer->m_timeUpdateTimeStamp = usec();
         g_pplayer->m_time_sec = m_time_sec;
         g_pplayer->m_time_msec = m_time_msec;
      }

      m_phys_iterations++;

      // Get the time until the next physics tick is done, and get the time
      // until the next frame is done
      // If the frame is the next thing to happen, update physics to that
      // point next update acceleration, and continue loop

      const float physics_diff_time = (float)((double)(m_nextPhysicsFrameTime - m_curPhysicsFrameTime)*(1.0 / DEFAULT_STEPTIME));
      //const float physics_to_graphic_diff_time = (float)((double)(initial_time_usec - m_curPhysicsFrameTime)*(1.0 / DEFAULT_STEPTIME));

      //if (physics_to_graphic_diff_time < physics_diff_time)          // is graphic frame time next???
      //{
      //      PhysicsSimulateCycle(physics_to_graphic_diff_time);      // advance physics to this time
      //      m_curPhysicsFrameTime = initial_time_usec;               // now current to the wall clock
      //      break;  //this is the common exit from the loop          // exit skipping accelerate
      //}                     // some rare cases will exit from while()


      // DJRobX's crazy latency-reduction code: Artificially lengthen the execution of the physics loop by X usecs, to give more opportunities to read changes from input(s) (try values in the multiple 100s up to maximum 1000 range, in general: the more, the faster the CPU is)
      //                                        Intended mainly to be used if vsync is enabled (e.g. most idle time is shifted from vsync-waiting to here)
      // FIXME the initial idea of this implementation is somewhat defeated by the fact that in single threaded mode, the main thread is mostly stalled waiting for GPU (solved in multithreaded mode) => remove ?
      #if !defined(ENABLE_BGFX)
      if (g_pplayer && g_pplayer->m_minphyslooptime > 0)
      {
         const uint64_t basetime = usec();
         const uint64_t targettime = ((uint64_t)g_pplayer->m_minphyslooptime * m_phys_iterations) + m_lastFlipTime;
         // If we're 3/4 of the way through the loop, fire a "controller sync" timer (timers with an interval set to -2) event so VPM can react to input.
         if (m_phys_iterations == 750 / ((int)g_pplayer->m_fps + 1))
            g_pplayer->FireTimers(-2);
         if (basetime < targettime)
         {
            g_pplayer->m_renderProfiler->EnterProfileSection(FrameProfiler::PROFILE_SLEEP);
            uSleep(targettime - basetime);
            g_pplayer->m_renderProfiler->ExitProfileSection();
         }
      }
      #endif
      // end DJRobX's crazy code

      // Anti hung mechanism
      if (g_pplayer && g_pplayer->m_playMode != Player::PlayMode::CaptureAttract)
      {
         const TestDriver::Timed tdInput(TestDriver::kInput);
         const uint64_t cur_time_usec = usec()
            - delta_frame; //!! one could also do this directly in the while loop condition instead (so that the while loop will really match with the current time), but that leads to some stuttering on some heavy frames

         // hung in the physics loop over 200 milliseconds or the number of physics iterations to catch up on is high (i.e. very low/unplayable FPS)
         // The test harness (VPX_PLAY_UNFOCUSED) runs slow frames on purpose: it must not skip simulated time, or the physics only runs for a part of the game time.
         if (!Player::PlayUnfocused() && ((cur_time_usec - initial_time_usec > 200000) || (m_physicsMaxLoops != 0 && m_phys_iterations > m_physicsMaxLoops)))
         { // can not keep up to real time
            m_curPhysicsFrameTime = initial_time_usec; // skip physics forward ... slip-cycles -> 'slowed' down physics
            m_nextPhysicsFrameTime = initial_time_usec + PHYSICS_STEPTIME;
            break; // go draw frame
         }

         //update keys, hid, plumb, nudge, timers, etc
         //const uint32_t sim_msec = (uint32_t)(m_curPhysicsFrameTime / 1000);
         const uint32_t cur_time_msec = (uint32_t)(cur_time_usec / 1000);

#if !defined(ENABLE_BGFX)
         // FIXME remove ? To be done correctly, we should process OS messages and sync back controller
         if (g_pplayer)
            g_pplayer->m_pininput.ProcessInput();
#endif

         // FIXME remove or at least move legacy ushock to a plugin
         ushock_output_update(/*sim_msec*/ cur_time_msec);
      }

#ifdef ACCURATETIMERS
      if (g_pplayer
#if !defined(ENABLE_BGFX)
         && (g_pplayer->GetVideoSyncMode() == VideoSyncMode::VSM_FRAME_PACING
            || g_pplayer->m_logicProfiler.Get(FrameProfiler::PROFILE_SCRIPT) <= 1000 * MAX_TIMERS_MSEC_OVERALL) // if overall script time per frame exceeded, skip
#endif
      )
      {
         const TestDriver::Timed tdTimers(TestDriver::kTimers);
         g_pplayer->FireTimers(0);
      }
#endif

      {
         const TestDriver::Timed tdDriver(TestDriver::kDriver);
         TestDriver::BeginCycle(m_time_msec);   // harness inputs by physics time (docs/TestMatrixDesign.txt section 2); nothing unless VPX_TEST_INPUT is set
      }

      // Since we are deriving forces/accelerations from velocities by doing a simple substract without scaling by delta time, we need dtime to be constant
      assert(fabs(physics_diff_time - ((double)PHYSICS_STEPTIME / (double)DEFAULT_STEPTIME)) < 1e-5);

      {
         const TestDriver::Timed tdNudge(TestDriver::kNudge);
         if (g_pplayer)
         {
            g_pplayer->m_pininput.m_nudgeHandler->StepOneMillisecond();
            g_pplayer->m_pininput.m_plungerHandler->StepOneMillisecond();
         }

         const Vertex2D cabinetAcceleration = GetCabinetAcceleration();
         m_plumbHandler.StepOneMillisecond(cabinetAcceleration);
         if (g_pplayer)
            g_pplayer->m_pininput.PlayNudgeRumble(cabinetAcceleration);
      }

      {
         const TestDriver::Timed tdMag(TestDriver::kMagnets);
         static const bool ctOn = getenv("VPX_CONTACT_TIMING") != nullptr && *getenv("VPX_CONTACT_TIMING") != 0;
         const uint64_t m0 = ctOn ? usec() : 0;
         Magnets::Step();   // coil currents and ball-ball forces: its time goes to the next cost-log row (mag_us)
         if (ctOn)
            s_ctMagUs += usec() - m0;
      }

      {
         const TestDriver::Timed tdVel(TestDriver::kVelocities);
         for (size_t i = 0; i < m_vmover.size(); i++)
            m_vmover[i]->UpdateVelocities();      // always on integral physics frame boundary (spinner, gate, flipper, plunger, ball)
      }

      //primary physics loop
      const uint64_t ptSim0 = pt.On() ? usec() : 0;
      {
         const TestDriver::Timed tdStep(TestDriver::kStep);
         PhysicsSimulateCycle(physics_diff_time); // main simulator call
      }
      const uint64_t ptSim1 = pt.On() ? usec() : 0;
      {
         const TestDriver::Timed tdDriver(TestDriver::kDriver);
         TestDriver::EndCycle(m_time_msec);
      }

      // Fast-ball watchdog: a ball above VPX_SPEED_WATCH m/s (default 8; 0 = off) is logged once per second with its last 8 cycles (position, velocity, who holds it), so a rare launch
      // can be traced afterwards. Costs a copy of a few floats per ball per cycle.
      {
         const TestDriver::Timed tdWatch(TestDriver::kWatchdog);
         static const double watch = []() { const char *const e = getenv("VPX_SPEED_WATCH"); return e != nullptr && *e ? atof(e) : 8.0; }();
         if (watch > 0.0 && g_pplayer != nullptr)
         {
            struct Rec { float x, y, z, vx, vy, vz; unsigned flags; uint32_t t; };
            static std::vector<std::array<Rec, 8>> hist;
            static std::vector<uint32_t> lastWarn;
            static unsigned head = 0;
            const size_t nb = g_pplayer->m_vball.size();
            if (hist.size() != nb)
            {
               hist.assign(nb, std::array<Rec, 8>());
               lastWarn.assign(nb, 0);
            }
            head = (head + 1) & 7;
            const float lim = (float)(watch * 18.53);
            for (size_t i = 0; i < nb; i++)
            {
               const HitBall &hb = g_pplayer->m_vball[i]->m_hitBall;
               const unsigned fl = (hb.m_contact != nullptr ? 1u : 0u) | (hb.m_rubberOwner != nullptr ? 2u : 0u) | (hb.m_bandOwner != nullptr ? 4u : 0u) | (hb.m_d.m_lockedInKicker ? 8u : 0u);
               hist[i][head] = { hb.m_d.m_pos.x, hb.m_d.m_pos.y, hb.m_d.m_pos.z, hb.m_d.m_vel.x, hb.m_d.m_vel.y, hb.m_d.m_vel.z, fl, g_pplayer->m_time_msec };
               const Vertex3Ds &v = hb.m_d.m_vel;
               if (!hb.m_d.m_lockedInKicker && v.x * v.x + v.y * v.y + v.z * v.z > lim * lim && g_pplayer->m_time_msec - lastWarn[i] > 1000)
               {
                  lastWarn[i] = g_pplayer->m_time_msec;
                  PLOGW << "Fast ball " << i << " at t=" << g_pplayer->m_time_msec << " ms: " << sqrt(v.x * v.x + v.y * v.y + v.z * v.z) / 18.53 << " m/s (VPX_SPEED_WATCH " << watch << "), last cycles (t, pos U, vel U/T, flags 1=contact 2=bat 4=band 8=kicker):";
                  for (unsigned k = 1; k <= 8; k++)
                  {
                     const Rec &r = hist[i][(head + k) & 7];
                     PLOGW << "   t=" << r.t << " pos (" << r.x << ", " << r.y << ", " << r.z << ") vel (" << r.vx << ", " << r.vy << ", " << r.vz << ") flags " << r.flags;
                  }
                  for (size_t j = 0; j < nb; j++)
                     if (j != i)
                     {
                        const HitBall &o = g_pplayer->m_vball[j]->m_hitBall;
                        PLOGW << "   other ball " << j << " pos (" << o.m_d.m_pos.x << ", " << o.m_d.m_pos.y << ", " << o.m_d.m_pos.z << ") vel (" << o.m_d.m_vel.x << ", " << o.m_d.m_vel.y << ", " << o.m_d.m_vel.z << ")";
                     }
               }
            }
         }
      }

      if (g_pplayer != nullptr && Magnets::HasGuides()) // air-core guide tube and stop
      {
         const TestDriver::Timed tdGuides(TestDriver::kGuides);
         for (const auto ball : g_pplayer->m_vball)
         {
            HitBall &hb = ball->m_hitBall;
            if (hb.m_d.m_lockedInKicker)
               continue;
            float p[3] = { hb.m_d.m_pos.x, hb.m_d.m_pos.y, hb.m_d.m_pos.z }, v[3] = { hb.m_d.m_vel.x, hb.m_d.m_vel.y, hb.m_d.m_vel.z };
            float L[3] = { hb.m_angularmomentum.x, hb.m_angularmomentum.y, hb.m_angularmomentum.z };
            if (Magnets::ConstrainBall(p, v, L, hb.m_d.m_radius, hb.m_d.m_mass, hb.Inertia()))
            {
               hb.m_d.m_pos = Vertex3Ds(p[0], p[1], p[2]);
               hb.m_d.m_vel = Vertex3Ds(v[0], v[1], v[2]);
               hb.m_angularmomentum = Vertex3Ds(L[0], L[1], L[2]);
            }
         }
      }

      // Store new position of balls for ball trails
      {
         const TestDriver::Timed tdTrails(TestDriver::kTrails);
         for (HitBall *const ball : m_vballs)
            ball->OnPhysicStepProcessed(m_time_msec);
      }

      //PLOGD << "PT: " << physics_diff_time << ' ' << physics_to_graphic_diff_time << ' ' << (uint32_t)(m_curPhysicsFrameTime/1000) << ' ' << (uint32_t)(initial_time_usec/1000) << ' ' << cur_time_msec;

      m_curPhysicsFrameTime = m_nextPhysicsFrameTime; // new cycle, on physics frame boundary
      m_nextPhysicsFrameTime += PHYSICS_STEPTIME;     // advance physics position

      if (pt.On())
      {
         const uint64_t sim = ptSim1 - ptSim0;
         pt.AddIter(sim, (usec() - ptIter0) - sim);
      }
      TestDriver::LoopEnd();
   } // end while (m_curPhysicsFrameTime < initial_time_usec)

   assert(m_curPhysicsFrameTime < m_nextPhysicsFrameTime);
   assert(m_curPhysicsFrameTime <= initial_time_usec);
   assert(initial_time_usec <= m_nextPhysicsFrameTime);

   // The physics is emulated by PHYSICS_STEPTIME, but the overall emulation time is more precise
   m_time_sec = (double)(initial_time_usec - m_startTime_usec) / 1000000.0;
   //m_time_sec = (double)(max(initial_time_usec, m_curPhysicsFrameTime) - m_startTime_usec) / 1000000.0;
   // m_time_msec = (uint32_t)((max(initial_time_usec, m_curPhysicsFrameTime) - m_startTime_usec) / 1000); // Not needed since PHYSICS_STEPTIME happens to be 1ms
   if (g_pplayer)
      g_pplayer->m_time_sec = m_time_sec;

   if (g_pplayer)
      g_pplayer->m_logicProfiler.ExitProfileSection();
}

void PhysicsEngine::PhysicsSimulateCycle(float dtime) // move physics forward to this time
{
   // PLOGD << "Cycle " << dtime;

   int StaticCnts = STATICCNTS; // maximum number of static counts
   // it's okay to have this code outside of the inner loop, as the ball hitrects already include the maximum distance they can travel in that timespan
   {
      const TestDriver::Timed tdSearch(TestDriver::kStockSearch);
      m_hitoctree_dynamic.Update();
   }

   static FILE *const ctLog = []() -> FILE * { const char *const p = getenv("VPX_CONTACT_TIMING"); FILE *const f = (p != nullptr && *p) ? fopen(p, "w") : nullptr; if (f != nullptr) fputs("t_ms,step_us,pre_us,near_us,bodies,islands,substeps,bats,band_bodies,rest_bodies,max_steps,items,probes,gather_us,setup_us,substep_us,write_us,rewinds,mag_us,snap_us,replays,replay_us\n", f); return f; }();   // cost log: one row per call that did work
   uint64_t ctNear = 0;
   if (ContactSystem::Unified())
   {
      const TestDriver::Timed tdContact(TestDriver::kContact);
      const uint64_t n0 = ctLog != nullptr ? usec() : 0;
      ContactSystem::EngageNear(dtime);
      ctNear = ctLog != nullptr ? usec() - n0 : 0;
   }

   while (dtime > 0.f)
   {
      // first find hits, if any +++++++++++++++++++++ 
      #ifdef DEBUGPHYSICS
         c_timesearch++;
      #endif
      float hittime = dtime; // begin time search from now ...  until delta ends
      const float hittimeStart = hittime;
      int splitCause = -1;
      const HitObject *splitObj = nullptr;
      const HitBall *splitBall = nullptr;
      bool splitStatic = false;
      std::optional<TestDriver::Timed> tdSearch(std::in_place, TestDriver::kStockSearch);   // the stock hit search, to m_recordContacts = false below

      // find earliest time where a flipper collides with its stop
      for (size_t i = 0; i < m_vFlippers.size(); ++i)
      {
         const float fliphit = m_vFlippers[i]->GetHitTime();
         //if ((fliphit >= 0.f) && !sign(fliphit) && (fliphit <= hittime))
         if ((fliphit > 0.f) && (fliphit <= hittime)) //!! >= 0.f causes infinite loop
            hittime = fliphit;
            splitCause = 100;
            splitObj = nullptr;
      }

      m_recordContacts = true;
      m_contacts.clear();

      #ifdef USE_EMBREE
            for (size_t i = 0; i < m_vballs.size(); i++)
               if (!m_vballs[i]->m_d.m_lockedInKicker
#ifdef C_DYNAMIC
                  && m_vballs[i]->m_dynamic > 0
#endif
                  ) // don't play with frozen balls
               {
                  m_vballs[i]->m_coll.m_hittime = hittime; // search upto current hittime
                  m_vballs[i]->m_coll.m_obj = nullptr;
               }

            if (!m_vballs.empty())
            {
               m_hitoctree.HitTestBall(m_vballs); // find the hit objects hit times
               m_hitoctree_dynamic.HitTestBall(m_vballs); // dynamic objects !! should reuse the same embree scene created already in m_hitoctree.HitTestBall!
            }
      #endif

      for (size_t i = 0; i < m_vballs.size(); i++)
      {
         HitBall *const pball = m_vballs[i];

         if (!pball->m_d.m_lockedInKicker
         #ifdef C_DYNAMIC
             && pball->m_dynamic > 0
         #endif
            ) // don't play with frozen balls
         {
            #ifndef USE_EMBREE
               pball->m_coll.m_hittime = hittime;          // search upto current hittime
               pball->m_coll.m_obj = nullptr;
            #endif

            // always check for playfield and top glass
            if (m_implicitPlayfieldMesh && !m_hitPlayfield.SurfaceOwns(pball->m_d)) // a ball the playfield force model holds is not tested
               DoHitTest(pball, &m_hitPlayfield, pball->m_coll);
            DoHitTest(pball, &m_hitTopGlass, pball->m_coll);

            #ifndef USE_EMBREE
               if (rand_mt_01() < 0.5f) // swap order of dynamic and static obj checks randomly
               {
                  m_hitoctree_dynamic.HitTestBall(pball, pball->m_coll); // dynamic objects
                  m_hitoctree.HitTestBall(pball, pball->m_coll);         // find the static hit objects hit times
               }
               else
               {
                  m_hitoctree.HitTestBall(pball, pball->m_coll);         // find the static hit objects hit times
                  m_hitoctree_dynamic.HitTestBall(pball, pball->m_coll); // dynamic objects
               }
            #endif
            const float htz = pball->m_coll.m_hittime; // this ball's hit time
            if (htz < 0.f) pball->m_coll.m_obj = nullptr; // no negative time allowed

            if (pball->m_coll.m_obj)                   // hit object
            {
               #ifdef DEBUGPHYSICS
                  ++c_hitcnts;                            // stats for display
                  if (/*pball->m_coll.m_hitRigid &&*/ pball->m_coll.m_hitdistance < -0.0875f) //rigid and embedded
                     ++c_embedcnts;
               #endif
               ///////////////////////////////////////////////////////////////////////////

               if (htz <= hittime)                     // smaller hit time??
               {
                  hittime = htz;                       // record actual event time
                  splitCause = (int)pball->m_coll.m_obj->m_ObjType;
                  splitObj = pball->m_coll.m_obj;
                  splitBall = pball;

                  if (hittime < STATICTIME)            // less than static time interval
                  {
                     /*if (!pball->m_coll.m_hitRigid) hittime = STATICTIME; // non-rigid ... set Static time
                     else*/ if (--StaticCnts < 0)
                     {
                        StaticCnts = 0;                // keep from wrapping
                        hittime = STATICTIME;
                        splitStatic = true;
                     }
                  }
               }
            }
         }
      } // end loop over all balls

      m_recordContacts = false;
      tdSearch.reset();

      // hittime now set ... or full frame if no hit
      // now update displacements to collide-contact or end of physics frame
      // !!!!! 2) move objects to hittime

      if (hittime > STATICTIME) StaticCnts = STATICCNTS; // allow more zeros next round

      std::optional<TestDriver::Timed> tdContact(std::in_place, TestDriver::kContact);   // the contact system's step, to the snapshot test below
      if (ContactSystem::Unified())
         ContactSystem::SaveStart();   // the state of the start of this sub-interval, for a replay (docs/UnifiedContactDesign.txt 11, point 2)
      const uint64_t ctT0 = ctLog != nullptr ? usec() : 0;
      if (hittime < hittimeStart && ContactProf::On())
      {
         std::string what = splitObj != nullptr ? typeid(*splitObj).name() : "flipper stop";
         if (splitObj != nullptr && splitObj->m_editable != nullptr)
            what += std::string(" part ") + splitObj->m_editable->GetName();
         if (splitBall != nullptr)
            what += splitBall->m_contact != nullptr ? " ballheld" : " ballfree";
         if (splitBall != nullptr && splitBall->m_contact == nullptr)
            what += " at " + std::to_string((int)(splitBall->m_d.m_pos.x / 100) * 100) + "," + std::to_string((int)(splitBall->m_d.m_pos.y / 100) * 100) + " z" + std::to_string((int)(splitBall->m_d.m_pos.z / 20) * 20)
                    + (splitBall->m_d.m_lockedInKicker ? " kickerlocked" : "") + (InPlayfieldHole(splitBall->m_d.m_pos) ? " inhole" : "") + (splitBall->m_bandOwner != nullptr ? " band" : "") + (splitBall->m_rubberOwner != nullptr ? " batowned" : "");
         ContactProf::Split(splitCause, splitStatic, what.c_str());
      }
      ContactSystem::Step(hittime); // contact-held balls: surface contacts and, for touching balls, ball-to-ball contact, stepped together
      if (ContactSystem::Unified() && ContactSystem::SnapshotTestOn())
         ContactSystem::SnapshotTest(hittime);   // test mode: the same step again from the snapshot must give the same result
      tdContact.reset();
      const uint64_t ctT1 = ctLog != nullptr ? usec() : 0;
      {
         const TestDriver::Timed tdMove(TestDriver::kMoveParts);
         for (size_t i = 0; i < m_vmover.size(); i++)
            m_vmover[i]->PreUpdateDisplacements(hittime); // flippers with a rubber-engaged ball move bat and ball together here
      }
      if (ctLog != nullptr)
      {
         const uint64_t ctT2 = usec();
         const ContactSystem::CycleStats &cs = ContactSystem::LastCycle();
         if (cs.bodies > 0 || ctT2 - ctT1 >= 3 || ctNear >= 3 || ContactSystem::LastSaveUs() >= 3)
         {
            static int rows = 0;
            fprintf(ctLog, "%u,%llu,%llu,%llu,%d,%d,%d,%d,%d,%d,%d,%d,%lld,%d,%d,%d,%d,%d,%llu,%u,%u,%u\n", g_pplayer->m_time_msec, (unsigned long long)(ctT1 - ctT0), (unsigned long long)(ctT2 - ctT1), (unsigned long long)ctNear, cs.bodies, cs.islands, cs.substeps, cs.bats, cs.bandBodies, cs.restBodies, cs.maxSteps, cs.items, cs.probes, cs.gatherUs, cs.setupUs, cs.substepUs, cs.writeUs, cs.rewinds, (unsigned long long)s_ctMagUs, ContactSystem::LastSaveUs(), s_ctReplays, s_ctReplayUs);
            s_ctReplays = s_ctReplayUs = 0;
            s_ctMagUs = 0;
            if ((++rows & 255) == 0)
               fflush(ctLog);
         }
      }
      {
         const TestDriver::Timed tdMove(TestDriver::kMoveParts);
         for (size_t i = 0; i < m_vmover.size(); i++)
            m_vmover[i]->UpdateDisplacements(hittime); // step 2: move the objects about according to velocities (spinner, gate, flipper, plunger, ball)
      }

      if (ContactSystem::Unified())
      {
         const TestDriver::Timed tdContactOverlaps(TestDriver::kContact);
         ContactSystem::CheckOverlaps(hittime);   // a free ball found inside a moving part: the sub-interval is replayed with the ball engaged at its first contact (docs/UnifiedContactDesign.txt 11, point 2)
         unsigned nRep = 0, usRep = 0;
         ContactSystem::TakeReplayStats(nRep, usRep);
         s_ctReplays += nRep;
         s_ctReplayUs += usRep;
      }

      // the events the contact step and the parts queued in this sub-interval: after the free balls have moved, as stock fires its collision events (docs/UnifiedContactDesign.txt 11, point 2)
      {
         const TestDriver::Timed tdDispatch(TestDriver::kContact);
         ContactSystem::Dispatch();
      }

      // find balls that need to be collided and script'ed (generally there will be one, but more are possible)

      for (size_t i = 0; i < m_vballs.size(); i++) // use m_vballs.size(), in case script deletes a ball
      {
         HitBall *const pball = m_vballs[i];

         if (
         #ifdef C_DYNAMIC
             pball->m_dynamic > 0 &&
         #endif
             pball->m_coll.m_obj && pball->m_coll.m_hittime <= hittime) // find balls with hit objects and minimum time
         {
            // now collision, contact and script reactions on active ball (object)+++++++++
            HitObject * const pho = pball->m_coll.m_obj; // object that ball hit in trials
            if (g_pplayer)
               g_pplayer->m_pactiveball = pball->m_pBall; // For script that wants the ball doing the collision
#ifdef DEBUGPHYSICS
            c_collisioncnt++;
#endif
            LogBallContact('C', pball->m_coll);
            if (m_hitSpeedEvents)
               pball->m_hitSpeed = std::max(0.f, -pball->m_coll.m_hitnormal.Dot(pball->m_d.m_vel)); // normal speed at the instant before the impact, whichever model resolves it
            pho->Collide(pball->m_coll);                 //!!!!! 3) collision on active ball
            // The playfield (the plane, the implicit mesh or an explicit playfield primitive) has no hit event of its own: the table gets one, for impacts only (not resting or rolling contact)
            if (m_hitSpeedEvents && g_pplayer != nullptr && pball->m_hitSpeed >= m_hitSpeedMin && i < g_pplayer->m_vball.size() && g_pplayer->m_vball[i] == pball->m_pBall
                && (pho == &m_hitPlayfield || (pho->m_ObjType == ePrimitive && pho->m_obj != nullptr && static_cast<Primitive *>(pho->m_obj)->IsPlayfield())))
               g_pplayer->m_ptable->FireVoidEvent(DISPID_HitEvents_Hit);
            pball->m_coll.m_obj = nullptr;               // remove trial hit object pointer

            // Collide may have changed the velocity of the ball, 
            // and therefore the bounding box for the next hit cycle
            if (i >= m_vballs.size() || m_vballs[i] != pball) // Ball still exists? may have been deleted from list
            {
               // collision script deleted the ball, back up one count
               --i;
               continue;
            }
            else
            {
               #ifdef C_DYNAMIC
                  // is this ball static? .. set static and quench
                  if (/*pball->m_coll.m_hitRigid &&*/ (pball->m_coll.m_hitdistance < (float)PHYS_TOUCH)) //rigid and close distance contacts //!! rather test isContact??
                  {
                     const float mag = pball->m_vel.x*pball->m_vel.x + pball->m_vel.y*pball->m_vel.y; // values below are taken from simulation
                     if (pball->m_drsq < 8.0e-5f && mag < 1.0e-3f*m_gravity*m_gravity / GRAVITYCONST / GRAVITYCONST && fabsf(pball->m_vel.z) < 0.2f*m_gravity / GRAVITYCONST)
                     {
                        if (--pball->m_dynamic <= 0)            //... ball static, cancels next gravity increment
                        {                                       // m_dynamic is cleared in ball gravity section
                           pball->m_dynamic = 0;
                           #ifdef DEBUGPHYSICS
                              c_staticcnt++;
                           #endif
                           pball->m_vel.x = pball->m_vel.y = pball->m_vel.z = 0.f; //quench the remaining velocity and set ...
                        }
                     }
                  }
               #endif
            }
         }
      }

      #ifdef DEBUGPHYSICS
         c_contactcnt = (uint32_t)m_contacts.size();
      #endif

      // Now handle contacts.

      // At this point UpdateDisplacements() was already called, so the state is different
      // from that at HitTest(). However, contacts have zero relative velocity, so
      // hopefully nothing catastrophic has happened in the meanwhile.

      // Maybe a two-phase setup where we first process only contacts, then only collisions
      // could also work.
      for (const CollisionEvent &ce : m_contacts)
         LogBallContact('T', ce);
#ifdef FIX_PHYSICS
      // Ball/ball contacts first: contacts can change a ball's normal velocity on its other supports (the friction stopping the
      // spin of a ball pressed against another ball or a wall pushes it up), which the static contacts resolved after them
      // (see HitBall::HandleStaticContact)
      for (HitBall *const pball : m_vballs)
         pball->m_contactDeltaV.SetZero();
      const bool forward = rand_mt_01() < 0.5f; // swap order of contact handling randomly
      for (int pass = 0; pass < 2; ++pass)
         for (size_t i = 0; i < m_contacts.size(); ++i)
         {
            CollisionEvent& contact = m_contacts[forward ? i : m_contacts.size() - 1 - i];
            if ((contact.m_obj->GetType() == eBall) == (pass == 0))
               contact.m_obj->Contact(contact, hittime);
         }
#else
      if (rand_mt_01() < 0.5f) // swap order of contact handling randomly
         for (size_t i = 0; i < m_contacts.size(); ++i)
            //if (m_contacts[i].m_hittime <= hittime) // does not happen often, and values then look sane, so do this check //!! why does this break some collisions (MM NZ&TT Reloaded Skitso, also CCC (Saloon))? maybe due to ball colliding with multiple things and then some sideeffect?
               m_contacts[i].m_obj->Contact(m_contacts[i], hittime);
      else
         for (int i = (int)m_contacts.size() - 1; i != -1; --i)
            //if (m_contacts[i].m_hittime <= hittime) // does not happen often, and values then look sane, so do this check //!! why does this break some collisions (MM NZ&TT Reloaded Skitso, also CCC (Saloon))? maybe due to ball colliding with multiple things and then some sideeffect?
               m_contacts[i].m_obj->Contact(m_contacts[i], hittime);
#endif

      m_contacts.clear();

      #ifdef C_BALL_SPIN_HACK
         // hacky killing of ball spin on resting balls (very low and very high spinning)
         for (HitBall *const hitBall : m_vballs)
         {
            if (hitBall->m_coll.m_hitdistance >= (float)PHYS_TOUCH)
               continue;

            const Vertex3Ds oldPos0 = hitBall->GetOldPosition(m_time_msec - 90); // Position 90ms ago
            if (oldPos0.x == FLT_MAX)
               continue;

            const Vertex3Ds oldPos1 = hitBall->GetOldPosition(m_time_msec - 80); // Position 80ms ago
            if (oldPos1.x == FLT_MAX)
               continue;

            /*
            const float mag = hitBall->m_vel.x*hitBall->m_vel.x + hitBall->m_vel.y*hitBall->m_vel.y; // values below are copy pasted from above
            if (hitBall->m_drsq < 8.0e-5f
                  && mag < 1.0e-3f*m_gravity*m_gravity / GRAVITYCONST / GRAVITYCONST
                  && fabsf(hitBall->m_vel.z) < 0.2f*m_gravity / GRAVITYCONST
                  && hitBall->m_angularmomentum.Length() < 0.9f*m_gravity / GRAVITYCONST)
                  //&& rand_mt_01() < 0.95f)
            {
               hitBall->m_angularmomentum *= 0.05f; // do not kill spin completely, otherwise stuck balls will happen during regular gameplay
            }
            */

            const Vertex3Ds diff_pos = oldPos0 - hitBall->m_d.m_pos;
            const float mag = diff_pos.x * diff_pos.x + diff_pos.y * diff_pos.y; // squared distance in the playfield plane, since 90ms
            const Vertex3Ds diff_pos2 = oldPos1 - hitBall->m_d.m_pos;
            const float mag2 = diff_pos2.x * diff_pos2.x + diff_pos2.y * diff_pos2.y; // squared distance in the playfield plane, since 80ms

            const float threshold = (hitBall->m_angularmomentum.x*hitBall->m_angularmomentum.x + hitBall->m_angularmomentum.y*hitBall->m_angularmomentum.y) / max(mag, mag2);

            if (!infNaN(threshold) && threshold > 666.f)
            {
               const float damp = clamp(1.0f - (threshold - 666.f) / 10000.f, 0.23f, 1.f); // do not kill spin completely, otherwise stuck balls will happen during regular gameplay
               hitBall->m_angularmomentum *= damp;
            }
         }
      #endif

      dtime -= hittime; //new delta .. i.e. time remaining

      m_swap_ball_collision_handling = !m_swap_ball_collision_handling; // swap order of ball-ball collisions

   } // end physics loop
}

string PhysicsEngine::GetPerfInfo(bool resetMax)
{
   const uint32_t lastPhysTime = g_pplayer ? g_pplayer->m_logicProfiler.GetPrev(FrameProfiler::PROFILE_PHYSICS) : 0;
   if (resetMax || lastPhysTime > m_phys_max)
      m_phys_max = lastPhysTime;

   if (resetMax || m_phys_iterations > m_phys_max_iterations)
      m_phys_max_iterations = m_phys_iterations;

   if (m_count == 0)
   {
      m_phys_total_iterations = m_phys_iterations;
      m_count = 1;
   }
   else
   {
      m_phys_total_iterations += m_phys_iterations;
      m_count++;
   }

   std::ostringstream info;
   info << std::fixed << std::setprecision(1);

   info << "Physics: " << m_phys_iterations << " iterations per frame (" << ((uint32_t)(m_phys_total_iterations / m_count)) << " avg " << m_phys_max_iterations
        << " max)\n";
#ifdef DEBUGPHYSICS
   info << std::setprecision(5);
   info << "Hits:" << c_hitcnts << " Collide:" << c_collisioncnt << " Ctacs:" << c_contactcnt;
#ifdef C_DYNAMIC
   info << " Static:" << c_staticcnt;
#endif
   info << " Embed:" << c_embedcnts << " TimeSearch:" << c_timesearch << '\n';
   info << " kDObjects:" << m_hitoctree_dynamic.GetObjectCount() << " kD:" << m_hitoctree_dynamic.GetNLevels() << '\n'
        << " QuadObjects:" << m_hitoctree.GetObjectCount() << " Quadtree:" << m_hitoctree.GetNLevels() << '\n'
        << " Traversed:" << c_traversed << " Tested:" << c_tested << " DeepTested:" << c_deepTested << '\n';
   info << std::setprecision(1);
#endif

   return info.str();
}

void PhysicsEngine::OnBallWallHit(const HitBall& ball, const Vertex3Ds& hitNormal, const float impactSpeed)
{
   for (HitPlunger* const plunger : m_vPlungers)
      plunger->OnBallWallHit(ball, hitNormal, impactSpeed);
}

Vertex2D PhysicsEngine::GetCabinetAcceleration() const { return g_pplayer ? g_pplayer->m_pininput.m_nudgeHandler->GetCabinetAcceleration() : Vertex2D(0.f, 0.f); }

BallControl *PhysicsEngine::GetBallControl() const { return (g_pplayer && g_pplayer->m_liveUI) ? &g_pplayer->m_liveUI->m_ballControl : nullptr; }

PlungerHandler *PhysicsEngine::GetPlungerHandler() const { return g_pplayer ? g_pplayer->m_pininput.m_plungerHandler.get() : nullptr; }

InputAction *PhysicsEngine::GetLaunchBallAction() const { return g_pplayer ? g_pplayer->m_pininput.GetInputActions()[g_pplayer->m_pininput.GetLaunchBallActionId()].get() : nullptr; }

void PhysicsEngine::PlayBallBallRumble(const float speed) const
{
   if (g_pplayer)
      g_pplayer->m_pininput.PlayBallBallRumble(speed);
}

void PhysicsEngine::PlaySlingshotRumble() const
{
   if (g_pplayer)
      g_pplayer->m_pininput.PlaySlingshotRumble();
}

void PhysicsEngine::PlayFlipperContactRumble(const float speed) const
{
   if (g_pplayer)
      g_pplayer->m_pininput.PlayFlipperContactRumble(speed);
}

void PhysicsEngine::PlayPlungerRumble(const float amplitude) const
{
   if (g_pplayer)
      g_pplayer->m_pininput.PlayPlungerRumble(amplitude);
}

void PhysicsEngine::PlayPlungerLaunchRumble(const float amplitude) const
{
   if (g_pplayer)
      g_pplayer->m_pininput.PlayPlungerLaunchRumble(amplitude);
}

// Diagnostics: one line per flipper with the numbers FlipperClearance uses.
std::string PhysicsEngine::FlipperDebug(const BallS &ball, const float horizon) const
{
   std::ostringstream o;
   for (const HitFlipper *const f : m_vFlippers)
      o << " [flipper base " << f->m_flipperMover.m_hitcircleBase.center.x << ',' << f->m_flipperMover.m_hitcircleBase.center.y << " angle " << f->m_flipperMover.m_angleCur
        << " speed " << f->m_flipperMover.m_angleSpeed << " gap " << f->GapToBall(ball) << " reach " << f->ReachToBall(ball, horizon) << ']';
   return o.str();
}
