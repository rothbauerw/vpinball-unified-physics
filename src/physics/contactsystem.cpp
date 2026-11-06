// license:GPLv3+

#include "core/stdafx.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <list>
#include <map>
#include <sstream>

#include "parts/ball.h"
#include "parts/flipper.h"
#include "contactsystem.h"
#include "hitrubber.h"
#include "contactparthost.h"
#include "hitball.h"
#include "hitflipper.h"
#include "magnet.h"
#include "magnetballs.h"
#include "magnetcore.h"
#include "turntable.h"
#include "contactprof.h"
#include "contacthingepart.h"
#include "contactmovingpart.h"
#include <map>
#include <unordered_map>
#include "physicsini.h"
#include "testdriver.h"
void ContactTraceBall(const char *stage, const HitBall *pb, bool items, int island);   // VPX_TRACE_BALL (end of file)
#include "ui/live/LiveUI.h"
#include "PhysicsEngine.h"
#include "physics/cordband.h"
#include "input/InputManager.h"

namespace
{
   constexpr double U = ContactSystem::kMetresPerUnit;
   constexpr double VS = ContactSystem::kVelocityScale;
   constexpr double kMassUnitKg = ContactSystem::kMassUnitKg;
   // Energy steps with a moving part in the group, this table start: per part, and the lines logged so far (rule 8; CordIsland, ReleaseHost)
   std::map<std::string, long> g_partEnergy;
   int g_partEnergyLines = 0;
   bool g_partEnergyTold = false;
}

namespace ContactSystem
{
   static bool s_unified = false, s_unifiedSet = false;

   void SetTable(const std::filesystem::path &tableFile)
   {
      g_partEnergy.clear();
      g_partEnergyLines = 0;
      g_partEnergyTold = false;
      const std::filesystem::path file = PhysicsIni::FindUnified(tableFile);
      PhysicsIni::Features f;
      const bool found = !file.empty() && PhysicsIni::ReadFeatures(file, f);
      s_unified = found && f.unified;   // no file, no [Features] or no  unified = 1  = stock VPX physics
      s_unifiedSet = true;
      ContactSettings::StartPlay(tableFile);   // (after the switch above: the settings ask it)
      if (s_unified)
         PLOGI << "Contact system: unified engine on (" << file.string() << "; dev fallback: stock hit search still runs for object types not yet converted)";
      else
         PLOGI << "Contact system: stock VPX physics, every new model off (" << (file.empty() ? std::string("no physics ini") : file.string() + " has no [Features] unified = 1") << ")";
   }

   bool Unified()
   {
      if (!s_unifiedSet)
         PLOGE << "Contact system: Unified() asked before the table's physics ini was read (SetTable); answering stock";
      return s_unified;
   }

   int Islands::Root(int i)
   {
      while (m_parent[i] != i)
         i = m_parent[i] = m_parent[m_parent[i]];
      return i;
   }

   void Islands::Build(const int n, const std::vector<std::pair<int, int>> &edges, const std::vector<std::pair<int, int>> &links)
   {
      m_parent.resize(n);
      for (int i = 0; i < n; i++)
         m_parent[i] = i;
      for (const auto &e : edges)
         m_parent[Root(e.first)] = Root(e.second);
      for (const auto &e : links)   // (balls on one band: connected without interacting)
         m_parent[Root(e.first)] = Root(e.second);
      count = 0;
      for (int root = 0; root < n; root++)
      {
         if (Root(root) != root)
            continue;
         if (count == list.size())
            list.emplace_back();
         Island &isl = list[count++];
         isl.members.clear();
         isl.edges.clear();
         for (int i = 0; i < n; i++)
            if (Root(i) == root)
               isl.members.push_back(i);
         for (const auto &e : edges)
            if (Root(e.first) == root)
               isl.edges.push_back(e);
      }
   }

}

// ---- The cycle's step -------------------------------------------------------------------------------------------------------------------------------------------------
// One physics cycle of every ball held by a force-model contact (the playfield, a post, a disc, ...). Each ball is advanced through the same substeps as before (its own surface
// contacts, gravity, magnets), but balls that touch each other are stepped TOGETHER, in lockstep: at every substep each ball's surface forces and the ball-to-ball Hertz forces
// are summed, then all of them are integrated. A stack of balls (balls pulled together by a magnet, a ball resting on another) is therefore one coupled system: the force through
// the stack is the weight plus the pull, with no stock impulse or position push in between.
//
// Phases: GatherBalls + GatherFree (SI state of every body) -> FindIslands -> per island: SetupIsland, Substeps (SurfaceForces + PairForces, then Integrate with the Sources),
// StepFailed or WriteBack -> Finish (releases, pair bookkeeping, deferred events).
using namespace RubberContact;

// The balls the contact system holds in an engine-stepped contact (a post, a wall rubber, the playfield, a disc): the rubbers each one touches, the shear state and the
// bookkeeping for the log and the release. HitBall::m_contact points at the ball's entry.
static std::list<BallContact> g_engaged;
static double g_frameUsedS = 0.0;   // how much of this physics step (1 ms) the sub-intervals before this one took (the engine splits a step at a hit)
static std::vector<const HitBall *> g_placed;   // the balls a script or the driver put where they are since the last cycle (NotePlaced)
static ContactSystem::CycleStats g_last;
static std::vector<std::function<void()>> g_deferred;      // events raised inside the step, fired after it (Finish)
// The test driver's cost split of this cycle (docs/TestMatrixDesign.txt section 3 point 10 b-d): null unless the driver times the cycle in detail. Per island: what set each
// substep (TestDriver::ContactCost::Why) and the smallest substep.
static TestDriver::ContactCost *s_cc = nullptr;
static int s_islWhy[TestDriver::ContactCost::kWhys];
static double s_islMinStep = 1.0;
void ContactSystem::DeferEvent(std::function<void()> f) { g_deferred.push_back(std::move(f)); }
static constexpr float kFlipperHorizonVpt = 1.0f;   // trace only: VPX time units (1 VPT = 10 ms) a bat is looked ahead for in the flip_clear numbers

struct ContactSystemImpl
{
   static void StepBalls(float dtime);
   static void Release(HitBall *ball, const char *why);
   static void Regather(BallContact &e, const Vertex3Ds &pos, float radiusU);
   static bool EngageImpl(ContactHost *host, HitBall *pb, const char *how);
   static void EngageNear(float dtime);
   static void CordNear(float dtime);

   // A ball in the step: its SI state, the substep's forces, and what the step needs to know about it.
   struct Body
   {
      ContactHost *o = nullptr;   // the rubber that engaged the ball (its settings drive the step)
      BallContact *e = nullptr;    // its engagement; null = a free ball that only joined an island through a pair
      HitBall *pb = nullptr;
      BallSI bs;
      V3 c, v, w;
      V3 extBaseVS;             // m/s per physics step: gravity and nudge (the part that does not depend on the position)
      Vertex3Ds posStart;
      double lScale = 1.0;
      double vCap = 50.0;
      int steps = 1;
      int stepsWhy = TestDriver::ContactCost::kEnd;   // the rule that set steps (test driver cost split, TestMatrixDesign section 3 point 10c)
      bool magnetsOn = false;
      bool extTaken = false;    // Sources::Begin took this cycle's gravity and nudge out of the velocity (they act in the substeps)
      V3 extOwe;                // the step's gravity and nudge for the sub-intervals after this one: out of the velocity while this one is stepped, back at the write-back
      double endPair = 0.0;     // the energy check: its half of what its ball pairs hold at the island's end
      bool magParts = false;    // the Option A ball-ball force is applied re-aimed at every substep (MagnetPartForces), not as the held increment
      double magAcc = 0.0;      // m/s2: the magnets' acceleration of this ball at the start of the step (the substep bound)
      double eddyRate = -1.0;   // eddy-current braking of the spin across the coil field, 1/s: evaluated once per cycle (negative = not yet); eddyDir = the field direction
      V3 eddyDir;
      V3 F, Mb;                 // forces of the substep
      V3 vS;                    // velocity at the start of the last substep, and that substep's length (0 = unknown: first substep of the cycle, or taken again): the ball's net acceleration
      double hS = 0.0;
      std::vector<TubeProbe> probes;   // the probes of this substep's start state, one per surface of the contact (arrival-time substeps: probed once, used for the step and the forces)
      bool probesReady = false;
      double omega2 = 0.0, deepest = 0.0, vnLast = 0.0;
      double eStart = 0.0, eScale = 0.0;   // the energy check: the ball's energy at the cycle's start, and what the allowance is a share of
      bool tripped = false;
      const char *tripWhy = "";
      int tripStep = 0;
      HitRubber *band = nullptr;   // a band that holds the ball (the group step takes the whole island)
   };

   struct BallHit
   {
      HitBall *a, *b;
      float speed;
   };

   // A ball started touching a surface of a part: the part's script hit event (threshold, duplicate filter) and, for the playfield, the table's impact event
   struct PartHit
   {
      const ContactSurface *s;
      HitBall *pb;
      float speed;        // along the surface normal, engine velocity units
      Vertex3Ds normal;   // from the surface to the ball
   };

   struct Cycle
   {
      double dt = 0.0;   // s
      const ContactSettings *cs = nullptr;
      std::vector<Body> bodies;
      std::vector<std::pair<HitBall *, const char *>> releases;
      std::vector<BallHit> hitEvents;    // deferred events: dispatched after the write-back
      std::vector<PartHit> partHits;     // likewise
      float hitFloor = 0.5f;             // slowest impact anything but a part's own hit event looks at (the plungers' 0.5, the table's hit_speed_min)
      std::vector<std::pair<int, int>> pairs;
      ContactSystem::Islands islands;
      RubberParams rpBall;               // steel on steel: Hertz with the pair's effective modulus (half of one ball's), no thin layer, no loss tangent
      std::vector<BallPairState *> pst;  // the current island's pair states
      int steps = 1;                     // the current island's substeps
      int stepsWhy = TestDriver::ContactCost::kEnd;   // the rule that set them (fixed-step islands; test driver cost split)
      int stepsPair = 1;                 // the substeps its ball pairs need alone (the arrival-time loop adds the surfaces' own)
      float dtimeF = 0.0f;               // the cycle in engine time units
      std::vector<std::pair<int, int>> links;
      std::vector<ContactMovingPart *> islandParts;   // the moving parts (gates, spinners) the current island steps (contactmovingpart.h)
      std::vector<std::pair<HitBall *, const char *>> bandReleases;   // balls the band lets go of (after the write-back: the release can fire the script)
      // the current island's Option A parts, resolved to its bodies once (MagnetPartForces): ia / ib = body index, -1 = a ball outside the island (held where it began)
      struct MagPartRT { int body = -1, ia = -1, ib = -1; double scale = 0; MagnetBalls::V3 F, u0, n, o, r0, pa0, pb0; bool pair = false; };
      std::vector<MagPartRT> magParts;
      const void *magIsl = nullptr;
   };

   // Sources: anything that acts on one body with no partner. Gravity and the nudge act in the substeps, not as one kick per physics step (a contact stepped without gravity
   // and kicked once per step rattles at the step rate); a magnet's pull follows the ball, so it is evaluated at every substep. Both are velocity increments.
   struct Sources
   {
      // Gravity and nudge: the engine has added this step's increment to the ball already (the stock contacts of everything else rely on it): this cycle's share of it is taken
      // back here and given back in small pieces, so the step ends with the same total.
      // Only while the playfield is one of the surfaces of the contact, or the ball has no contact at all: with the stock plane contact under it the ball needs the engine's kick to stay on the plane.
      static void Begin(Body &w, const double dt, const bool always = false)
      {
         HitBall *const pb = w.pb;
         // A free ball (no contact of its own, in an island only through a ball pair) has nothing of the stock engine under it: its gravity, nudge and magnet pull act in the substeps like an engaged ball's.
         const bool onPlayfieldModel = always || w.e == nullptr || std::any_of(w.e->items.begin(), w.e->items.end(), [](const SurfaceContact &it) { return it.r->IsFloor(); });
         w.extTaken = onPlayfieldModel;
         const Vertex3Ds extInc = onPlayfieldModel ? pb->ExternalVelocityIncrement() : Vertex3Ds(0.f, 0.f, 0.f);   // what the engine added this step: gravity, nudge, magnets at the start position
         const Vertex3Ds extBase = onPlayfieldModel ? pb->ExternalVelocityIncrement(false) : Vertex3Ds(0.f, 0.f, 0.f);   // the part that does not depend on the position
         w.magnetsOn = onPlayfieldModel && pb->m_magnetic && Magnets::Active();   // a magnet's pull follows the ball: evaluated at every substep, not once per cycle
         w.magParts = w.magnetsOn && Magnets::HasBallBallParts(pb);
         w.magAcc = w.magnetsOn ? V3((extInc.x - extBase.x) * VS, (extInc.y - extBase.y) * VS, (extInc.z - extBase.z) * VS).Length() / PHYSICS_STEPTIME_S : 0.0;
         const double cycleShare = dt / PHYSICS_STEPTIME_S;
         w.extBaseVS = V3(extBase.x * VS, extBase.y * VS, extBase.z * VS);
         // the step's increment is all in the velocity: the sub-intervals after this one (a step split at a hit) get theirs back at the write-back, not while this one is stepped
         const double later = std::max(PHYSICS_STEPTIME_S - g_frameUsedS - dt, 0.0) / PHYSICS_STEPTIME_S;
         w.extOwe = V3(extInc.x * VS * later, extInc.y * VS * later, extInc.z * VS * later);
         w.v = w.v - V3(extInc.x * VS * cycleShare, extInc.y * VS * cycleShare, extInc.z * VS * cycleShare) - w.extOwe;
      }

      // Eddy-current braking of the ball's spin across the coil field, the one form for the old loop and the group: a steel sphere turning at w in a field B feels a torque (2 pi / 15) sigma R^5 B^2 w
      // against the part of w across B (a spin about the field lines induces nothing). Its rate and direction into the body. The field changes little within a cycle: one evaluation per cycle, at the ball's place then.
      static void EddyBrake(Body &w)
      {
         if (w.eddyRate >= 0.0)
            return;
         w.eddyRate = 0.0;
         if (!w.magnetsOn)
            return;
         const float mp[3] = { (float)(w.c.x / U), (float)(w.c.y / U), (float)(w.c.z / U) };
         const double sigma = Magnets::EddyConductivity();
         double Bf[3];
         if (sigma > 0.0 && Magnets::FieldAt(mp, w.pb->m_d.m_radius, Bf))
         {
            const double b2 = Bf[0] * Bf[0] + Bf[1] * Bf[1] + Bf[2] * Bf[2];
            if (b2 > 1e-8)
            {
               w.eddyDir = V3(Bf[0], Bf[1], Bf[2]) / std::sqrt(b2);
               const double R2 = w.bs.R * w.bs.R;
               w.eddyRate = (2.0 * M_PI / 15.0) * sigma * R2 * R2 * w.bs.R * b2 / w.bs.Inertia();
            }
         }
      }

      // One substep's share (aScale = h / one physics step), at the body's current position.
      static void Add(Body &w, const double aScale)
      {
         w.v = w.v + w.extBaseVS * aScale;
         if (w.magnetsOn)
         {
            const float mp[3] = { (float)(w.c.x / U), (float)(w.c.y / U), (float)(w.c.z / U) };
            float mdv[3] = { 0.f, 0.f, 0.f };   // the other balls' part: held for the cycle (dipoles), or re-aimed at every substep (Option A: MagnetPartForces)
            if (!w.magParts) { mdv[0] = w.pb->m_magPair[0]; mdv[1] = w.pb->m_magPair[1]; mdv[2] = w.pb->m_magPair[2]; }
            Magnets::AddVelocityIncrement(mp, w.pb->m_d.m_radius, w.pb->m_d.m_mass, mdv);
            const double k = VS * aScale;
            w.v = w.v + V3(mdv[0] * k, mdv[1] * k, mdv[2] * k);
            // Exact exponential decay over the substep.
            EddyBrake(w);
            if (w.eddyRate > 0.0)
            {
               const V3 wPerp = w.w - w.eddyDir * w.w.Dot(w.eddyDir);
               const double x = w.eddyRate * aScale * PHYSICS_STEPTIME_S;
               w.w = w.w - wPerp * (x < 0.05 ? x * (1.0 - 0.5 * x) : 1.0 - std::exp(-x));
            }
         }
      }
   };

   // The constants of a body come from the ball's radius, mass and inertia each cycle
   static void LoadBody(Body &w, HitBall *const pb)
   {
      w.pb = pb;
      w.bs.R = (double)pb->m_d.m_radius * U;
      w.bs.m = kMassUnitKg * (double)pb->m_d.m_mass;
      w.lScale = (double)pb->Inertia() * DEFAULT_STEPTIME_S;
      w.posStart = pb->m_d.m_pos;
   }

   // SI state of a ball
   static void LoadState(Body &w)
   {
      HitBall *const pb = w.pb;
      const double wScale = 1.0 / w.lScale;
      w.c = V3((double)w.posStart.x * U, (double)w.posStart.y * U, (double)w.posStart.z * U);
      w.v = V3((double)pb->m_d.m_vel.x * VS, (double)pb->m_d.m_vel.y * VS, (double)pb->m_d.m_vel.z * VS);
      w.w = V3((double)pb->m_angularmomentum.x * wScale, (double)pb->m_angularmomentum.y * wScale, (double)pb->m_angularmomentum.z * wScale);
   }

   // How fast the ball approaches the surfaces next to it (m/s, along their normals), for the choice of the step. When the last cycle had a contact the engagement
   // remembers it (lastVn); a ball hopping a micrometre above the surface has none, and then this says whether it is rolling or about to hit: the nearest surface within a few
   // release gaps decides, a ball farther out (in flight) keeps the impact step.
   static double ApproachSpeed(Body &w, BallContact &e, const double dt)
   {
      CPROF(ApproachSpeed);
      // A surface other than the one the ball rests on (a wall, a gate or spinner plate that moves) that the ball, or the surface, can reach within this cycle counts with its closing speed
      // (relative to the surface's own velocity): a ball rolling on the playfield must meet it with the impact step, not the rolling step (a light plate takes the ball's momentum in a few
      // microseconds).
      double vnNear = 0.0;
      for (SurfaceContact &it : e.items)
      {
         if (it.r->IsFloor() || it.r->MovingPart() != nullptr)   // (a gate or spinner bounds the step by its own arrival rule: Part::BallStep in SetupIsland)
            continue;
         it.r->Refresh(it, w.bs.R, 2.0 * w.v.Length() * dt + 0.005, w.c);
         if (it.NoFaces(w.c))
            continue;
         const TubeProbe pr = it.r->Probe(w.bs, it, w.c);
         if (pr.dist <= 0.0)
            continue;
         if (it.partOut)   // (island_reach: a part out of this ball's reach this cycle)
            continue;
         V3 vs, ws;
         if (it.r->MovingPart() != nullptr)
            it.r->VelocityAt(pr, vs, ws);
         const V3 vrel = w.v - vs;
         const double gap = pr.dist - (w.bs.R + pr.rho);
         if (gap < 2.0 * vrel.Length() * dt + 3.0 * it.r->ReleaseGapM() + 1e-3)
            vnNear = std::max(vnNear, std::max(0.0, -vrel.Dot(pr.normal)));
      }
      if (e.lastVn < 1e8)
         return std::max(e.lastVn, vnNear);
      double vn = 1e9;
      for (SurfaceContact &it : e.items)
      {
         if (it.r->MovingPart() != nullptr)
            continue;
         it.r->Refresh(it, w.bs.R, 2.0 * w.v.Length() * dt + 0.005, w.c);
         const TubeProbe pr = it.r->Probe(w.bs, it, w.c);
         if (pr.dist <= 0.0)
            continue;
         const double gap = pr.dist - (w.bs.R + pr.rho);
         if (gap > 3.0 * it.r->ReleaseGapM() + 1e-6)
            continue;   // not close to this one
         const double a = std::max(0.0, -w.v.Dot(pr.normal));
         vn = vn >= 1e8 ? a : std::max(vn, a);
      }
      return vn >= 1e8 ? vn : std::max(vn, vnNear);
   }

   // The balls held by a force-model contact
   static void GatherBalls(Cycle &cy)
   {
      CPROF(GatherBalls);
      const double dt = cy.dt;
      for (auto eit = g_engaged.begin(); eit != g_engaged.end();)
      {
         BallContact &e = *eit++;   // (a release below erases only this entry)
         if (e.ball->m_sleepGroup != 0)   // pile sleep: held where it is, its contact state kept (SleepInStep)
            continue;
         ContactHost *const o = e.host;   // null for a ball only a band holds
         if (o != nullptr && !o->m_set.contact)   // switched off in the ini while a ball is engaged
         {
            Release(e.ball, "rubber contact switched off");
            continue;
         }
         {
            HitBall *const pb = e.ball;
            if (pb->m_d.m_lockedInKicker)   // captured by a kicker: the kicker holds it, not the contact
            {
               cy.releases.emplace_back(pb, "kicker");
               continue;
            }
            Body w;
            w.o = o;
            w.e = &e;
            LoadBody(w, pb);
            {
               CPROF(Regather);
               Regather(e, w.posStart, pb->m_d.m_radius);   // (over a kicker's hole the playfield is left out)
            }
            if (PlayfieldHolesExist() && InPlayfieldHole(pb->m_d.m_pos) && std::any_of(e.items.begin(), e.items.end(), [](const SurfaceContact &it) { return it.r->IsPlayfield(); }))
            {
               cy.releases.emplace_back(pb, "kicker hole");   // rolled into a non-legacy kicker's circle with nothing else to hold it
               continue;
            }
            LoadState(w);
            // The settings that set the substep are those of the finest surface the ball touches or can reach: the one that engaged it and every surface in its contact set
            // (measured: replay of the 2026-10-13 launch, a ball just let go by a bat stepped once per cycle against the playfield, 28 um deep at the start of the next, 247 N held
            // for 1 ms, +3 m/s: it is the contact's own stiffness that must limit the step, whatever engaged the ball).
            const RubberObjectSettings *so = o != nullptr ? o->Settings() : nullptr;   // (a gate or spinner has none: its step is its own, SetupIsland)
            for (const SurfaceContact &it : e.items)
               if (const RubberObjectSettings *const s = it.r->Settings())
                  if (so == nullptr || s->substepUs < so->substepUs)
                     so = s;
            const int stepsN = so != nullptr ? std::max(1, (int)std::ceil(dt / (so->substepUs * 1e-6))) : 1;
            // A ball at rest on the rubber (slow at both ends of the cycle) does not need the impact step: the contact's own frequency sets the stable step (omega * h <= 0.2),
            // up to rest_substep_us. The rest depth is the force balance, not a property of the step.
            int steps = stepsN;
            w.stepsWhy = stepsN > 1 ? TestDriver::ContactCost::kSetting : TestDriver::ContactCost::kEnd;
            if (so != nullptr)
            {
               const RubberObjectSettings &set = *so;
               if (set.restSubstepUs > set.substepUs && w.v.Length() < set.restSpeedMmS * 1e-3 && e.lastSpeed < set.restSpeedMmS * 1e-3 && e.omega > 0.0)
               {
               const double hRest = std::min(set.restSubstepUs * 1e-6, set.restOmegaH / e.omega);
               steps = std::max(1, std::min(stepsN, (int)std::ceil(dt / hRest)));
               w.stepsWhy = steps > 1 ? TestDriver::ContactCost::kSurfaceStiffness : TestDriver::ContactCost::kEnd;
               }
               else if (set.restSubstepUs > set.substepUs && set.rollTravelUm > 0.0 && e.omega > 0.0 && ApproachSpeed(w, e, dt) < set.restSpeedMmS * 1e-3)
               {
               // Rolling or sliding along a surface (nothing approaches it): the contact's stiffness and the distance covered per step limit the step, not the impact step.
               const double hRoll = std::min({ set.restSubstepUs * 1e-6, set.restOmegaH / e.omega, set.rollTravelUm * 1e-6 / std::max(w.v.Length(), 1e-6) });
               steps = std::max(1, std::min(stepsN, (int)std::ceil(dt / hRoll)));
               w.stepsWhy = steps > 1 ? TestDriver::ContactCost::kSurfaceStiffness : TestDriver::ContactCost::kEnd;
               }
            }
            w.band = e.band;
            w.steps = steps;
            Sources::Begin(w, dt);
            w.vCap = std::max(50.0, 4.0 * w.v.Length() + 5.0);   // m/s
            {
               CPROF(RefreshCycle);
               for (SurfaceContact &it : e.items)
                  it.r->Refresh(it, w.bs.R, 2.0 * w.v.Length() * dt + 0.005, w.c);   // a surface: the faces the ball can reach this cycle
            }
            cy.bodies.push_back(std::move(w));
         }
      }
   }

   // A free ball (not held by any contact) that can reach another ball, held or free, this cycle joins the island with no surface contacts of its own: the engine has given it
   // gravity already, and the stock collision is skipped for the pair, so the two are solved together here. Chains grow until no free ball is within reach of a member.
   static void GatherFree(Cycle &cy)
   {
      CPROF(GatherFree);
      const ContactSettings &cs = *cy.cs;
      if (!cs.ballContact)
         return;
      std::vector<Body> cand;
      for (auto *const ball : g_pplayer->m_vball)
      {
         HitBall *const pb = &ball->m_hitBall;
         if (pb->m_bandOwner != nullptr || pb->m_contact != nullptr || pb->m_sleepGroup != 0 || !PairEligible(pb))
            continue;
         Body w;
         LoadBody(w, pb);
         LoadState(w);
         w.vCap = std::max(50.0, 4.0 * w.v.Length() + 5.0);
         cand.push_back(std::move(w));
      }
      if (cand.empty() || (cy.bodies.empty() && cand.size() < 2))
         return;
      const auto reaches = [&](const Body &x, const Body &y) {
         const double gap = (x.c - y.c).Length() - x.bs.R - y.bs.R;
         return gap < cs.ballReachMm * 1e-3 + 2.0 * (x.v - y.v).Length() * cy.dt;
      };
      std::vector<char> taken(cand.size(), 0);
      for (bool grew = true; grew;)
      {
         grew = false;
         const size_t end = cy.bodies.size();
         for (size_t k = 0; k < cand.size(); k++)
         {
            if (taken[k])
               continue;
            bool close = false;
            for (size_t i = 0; i < end && !close; i++)
               close = reaches(cy.bodies[i], cand[k]);
            if (!close)   // (a pair of candidates that are close to each other but to nothing else starts a cluster of free balls)
               for (size_t m = 0; m < cand.size() && !close; m++)
                  close = m != k && reaches(cand[m], cand[k]);
            if (close)
            {
               taken[k] = 1;
               Sources::Begin(cand[k], cy.dt);   // (its gravity, nudge and magnet pull: taken back, then given in the substeps)
               cy.bodies.push_back(cand[k]);
               grew = true;
            }
         }
      }
   }

   // island_reach: a ball and a moving part share an island only if the ball can reach the part's surface within this cycle (the reachable set over the cycle, runs/2026-11-05-island-reach/design.md).
   // A ball's travel: the speeds of every ball of its group (a hit hands speed on), what their forces add, and twice the speed of each part in the group (a ball leaves a moving surface at up to twice
   // its speed). A part's travel: its whole range of motion once a ball can reach it (a hit can move it anywhere it can go), else what its own drives give it (SurfaceSpeedBound). Grown to a fixed point:
   // the bounds only grow. A part no ball can reach is left out of every ball's solve (partOut) and moves alone; one a ball reaches is in that ball's island only. The gap is to its real surface (Probe).
   // ball_pair_reach: the one ball-ball reach test (island pairing in ReachLinks, the sleep wake in SleepInStep). V = a ball's speed bound over the cycle: the kinetic energy
   // that can reach it (spin included) given all to it, sqrt(2 E / m), plus what forces and moving parts add; a held (sleeping) ball has V = 0 at its held position
   static double BallSpeedBound(const double e, const double m, const double add) { return std::sqrt(2.0 * e / m) + add; }
   static bool BallsMayMeet(const double gap, const double vi, const double vj, const double T) { return gap < (vi + vj) * T + 100e-6; }   // (100 um: the contact law's approach margin)

   // ball_pair_reach rule 8: two balls in different islands last cycle that overlap at this cycle's start were missed by the bound
   inline static std::vector<std::pair<const HitBall *, uintptr_t>> s_ballGroupLast;   // (ball, its group) of the last cycle: its island's root ball, or its sleep group (held balls)
   static void BallPairMissCheck(const Cycle &cy, const std::vector<double> &bgap)
   {
      const int nw = (int)cy.bodies.size();
      const auto rootOf = [](const HitBall *pb) -> uintptr_t { for (const auto &g : s_ballGroupLast) if (g.first == pb) return g.second; return 0; };
      for (int i = 0; i < nw; i++)
         for (int j = i + 1; j < nw; j++)
            if (bgap[(size_t)i * nw + j] < 0.0)
            {
               const uintptr_t ri = rootOf(cy.bodies[i].pb), rj = rootOf(cy.bodies[j].pb);
               if (ri == 0 || rj == 0 || ri == rj)
                  continue;
               static int warned = 0;
               const char *const ni = TestDriver::BallName(cy.bodies[i].pb), *const nj = TestDriver::BallName(cy.bodies[j].pb);
               std::ostringstream o;
               o << "Ball pair reach: balls " << (ni != nullptr ? ni : "") << " and " << (nj != nullptr ? nj : "") << " overlap " << -bgap[(size_t)i * nw + j] * 1e3 << " mm, in different islands last cycle (a missed contact)";
               PLOGW << o.str();
               if (warned++ == 0 && g_pplayer->m_liveUI != nullptr)
                  g_pplayer->m_liveUI->PushNotification(o.str(), 6000);
            }
   }
   // ball_pair_reach diagnostic: a pair only the closure groups (the pairwise rule would not), logged when that grouping starts, with the terms of each V (200 lines at most)
   inline static std::vector<std::pair<const HitBall *, const HitBall *>> s_pairReachLast, s_pairReachNow;
   template <class Terms> static void BallPairReachLog(const Cycle &cy, const int i, const int j, const double gap, Terms terms)
   {
      const std::pair<const HitBall *, const HitBall *> key(cy.bodies[i].pb, cy.bodies[j].pb);
      s_pairReachNow.push_back(key);
      static int lines = 0;
      if (lines >= 200 || std::find(s_pairReachLast.begin(), s_pairReachLast.end(), key) != s_pairReachLast.end())
         return;
      lines++;
      const char *const ni = TestDriver::BallName(cy.bodies[i].pb), *const nj = TestDriver::BallName(cy.bodies[j].pb);
      const V3 &ci = cy.bodies[i].c, &cj = cy.bodies[j].c;
      PLOGI << "Ball pair reach t=" << g_pplayer->m_time_msec << ": " << (ni != nullptr ? ni : "?") << " (" << ci.x << " " << ci.y << " " << ci.z << ") and " << (nj != nullptr ? nj : "?") << " (" << cj.x << " "
            << cj.y << " " << cj.z << ") gap " << gap * 1e3 << " mm < (Vi + Vj) T: i " << terms(i) << "; j " << terms(j);
   }
   template <class Find> static void BallPairMissKeep(const Cycle &cy, Find &find)
   {
      s_ballGroupLast.clear();
      for (int i = 0; i < (int)cy.bodies.size(); i++)
         s_ballGroupLast.emplace_back(cy.bodies[i].pb, (uintptr_t)cy.bodies[find(i)].pb);
      for (auto *const b : g_pplayer->m_vball)   // held (asleep) this cycle: grouped by sleep group, checked at its held position when it wakes
         if (b->m_hitBall.m_sleepGroup != 0)
            s_ballGroupLast.emplace_back(&b->m_hitBall, (uintptr_t)b->m_hitBall.m_sleepGroup);
   }

   static void ReachLinks(Cycle &cy)
   {
      constexpr double kApproachM = 100e-6;   // m: the contact law's approach margin (SetupIsland kApproachMarginM): a gap under it is stepped as a contact
      const int nw = (int)cy.bodies.size();
      const double T = cy.dt;
      struct Edge { int i, p; SurfaceContact *it; double gap; bool on; };
      std::vector<Edge> edges;
      std::vector<ContactMovingPart *> parts;
      for (int i = 0; i < nw; i++)
         if (cy.bodies[i].e != nullptr)
            for (SurfaceContact &it : cy.bodies[i].e->items)
               if (ContactMovingPart *const mp = it.r->MovingPart())
               {
                  it.partOut = false;
                  const Body &w = cy.bodies[i];
                  mp->SetBall(w.bs);
                  const double gap = -mp->Probe(w.c).delta;
                  if (it.partOutLast && gap < 0.0)   // left out last cycle, and the ball is inside it now: the bound missed a contact (rule 8)
                  {
                     static int warned = 0;
                     const char *const name = TestDriver::BallName(w.pb);
                     std::ostringstream o;
                     o << "Island reach: ball " << (name != nullptr ? name : "") << " is " << -gap * 1e3 << " mm inside " << it.r->Name() << ", which it was left out of last cycle (a missed contact)";
                     PLOGW << o.str();
                     if (warned++ == 0 && g_pplayer->m_liveUI != nullptr)
                        g_pplayer->m_liveUI->PushNotification(o.str(), 6000);
                  }
                  it.partOutLast = false;
                  const int p = (int)(std::find(parts.begin(), parts.end(), mp) - parts.begin());
                  if (p == (int)parts.size())
                     parts.push_back(mp);
                  edges.push_back({ i, p, &it, gap, false });
               }
      const bool ballPairs = cy.cs->ballPairReach && cy.cs->ballContact && nw >= 2;
      if (edges.empty() && !ballPairs)
         return;
      std::vector<int> root(nw);
      for (int i = 0; i < nw; i++)
         root[i] = i;
      const auto find = [&](int i) { while (root[i] != i) i = root[i] = root[root[i]]; return i; };
      const auto unite = [&](const int a, const int b) { root[find(a)] = find(b); };
      for (const auto &pr : cy.pairs)
         unite(pr.first, pr.second);
      for (const auto &lk : cy.links)
         unite(lk.first, lk.second);
      const int np = (int)parts.size();
      std::vector<double> sPart(np, -1.0), range(np, 0.0);
      std::vector<int> partBall(np, -1);   // the first ball that can reach it
      std::vector<const std::string *> partName(np, nullptr);
      for (const Edge &ed : edges)
         if (partName[ed.p] == nullptr)
            partName[ed.p] = &ed.it->r->Name();
      for (const Edge &ed : edges)
         if (sPart[ed.p] < 0.0)
         {
            parts[ed.p]->SetBall(cy.bodies[ed.i].bs);
            sPart[ed.p] = parts[ed.p]->SurfaceSpeedBound(T);
            range[ed.p] = parts[ed.p]->MotionRange();
         }
      // ball_pair_reach: balls i, j share an island when gap_ij < (V_i + V_j) T + margin, V = the group speed bound below (plain sum: walls and hits can turn a ball, so |vi - vj| is no bound)
      std::vector<char> paired;
      std::vector<double> bgap;
      if (ballPairs)
      {
         paired.assign((size_t)nw * nw, 0);
         bgap.resize((size_t)nw * nw);
         for (const auto &pr : cy.pairs)
            paired[(size_t)pr.first * nw + pr.second] = 1;
         for (int i = 0; i < nw; i++)
            for (int j = i + 1; j < nw; j++)
               bgap[(size_t)i * nw + j] = (cy.bodies[i].c - cy.bodies[j].c).Length() - cy.bodies[i].bs.R - cy.bodies[j].bs.R;
         BallPairMissCheck(cy, bgap);
      }
      // ball_pair_reach: the balls' part of a ball's speed bound is the group's kinetic energy (spin included) given all to that ball, sqrt(2 E / m): a hit only shares energy out
      // (a light ball hit by a heavy one leaves at up to twice its speed, beyond the sum of speeds); equal masses: sqrt(sum v^2) <= sum |v|. Without it: the sum of the group's speeds.
      std::vector<double> vc(nw), eg(ballPairs ? nw : 0);
      const auto V = [&](const int i) { const int g = find(i); return ballPairs ? BallSpeedBound(eg[g], cy.bodies[i].bs.m, vc[g]) : vc[g]; };
      for (bool grew = true; grew;)
      {
         grew = false;
         std::fill(vc.begin(), vc.end(), 0.0);
         std::fill(eg.begin(), eg.end(), 0.0);
         for (int i = 0; i < nw; i++)
         {
            const Body &w = cy.bodies[i];
            vc[find(i)] += (ballPairs ? 0.0 : w.v.Length()) + (w.extBaseVS.Length() / PHYSICS_STEPTIME_S + w.magAcc + w.F.Length() / w.bs.m) * T;
            if (ballPairs)
               eg[find(i)] += 0.5 * w.bs.m * w.v.Dot(w.v) + 0.5 * w.bs.Inertia() * w.w.Dot(w.w);
         }
         for (int p = 0; p < np; p++)
            if (partBall[p] >= 0)
               vc[find(partBall[p])] += 2.0 * sPart[p];
         for (Edge &ed : edges)
         {
            if (ed.on)
               continue;
            const double travelPart = partBall[ed.p] >= 0 ? range[ed.p] : sPart[ed.p] * T;
            if (ed.gap < V(ed.i) * T + travelPart + kApproachM)
            {
               ed.on = grew = true;
               if (partBall[ed.p] < 0)
                  partBall[ed.p] = ed.i;
               else if (find(ed.i) != find(partBall[ed.p]))
                  unite(ed.i, partBall[ed.p]);
            }
         }
         if (ballPairs)
            for (int i = 0; i < nw; i++)
               for (int j = i + 1; j < nw; j++)
                  if (!paired[(size_t)i * nw + j] && BallsMayMeet(bgap[(size_t)i * nw + j], V(i), V(j), T))
                  {
                     BallPairReachLog(cy, i, j, bgap[(size_t)i * nw + j], [&](const int k) {   // (the terms of V, read before the pair joins the groups)
                        const int g = find(k);
                        std::ostringstream o;
                        double parts2S = 0.0;
                        o << "V " << V(k) << " = energy " << std::sqrt(2.0 * eg[g] / cy.bodies[k].bs.m);
                        for (int p = 0; p < np; p++)
                           if (partBall[p] >= 0 && find(partBall[p]) == g)
                           {
                              parts2S += 2.0 * sPart[p];
                              o << " + 2S " << 2.0 * sPart[p] << " (" << *partName[p] << ")";
                           }
                        o << " + forces " << vc[g] - parts2S << ", |v| " << cy.bodies[k].v.Length();
                        return o.str(); });
                     paired[(size_t)i * nw + j] = 1;
                     cy.pairs.emplace_back(i, j);
                     unite(i, j);
                     grew = true;
                  }
      }
      if (ballPairs)
      {
         BallPairMissKeep(cy, find);
         std::swap(s_pairReachLast, s_pairReachNow);
         s_pairReachNow.clear();
      }
      for (const Edge &ed : edges)
      {
         const int pb = partBall[ed.p];
         if (pb < 0 || find(ed.i) != find(pb))
            ed.it->partOut = true;   // solved in another island, or alone
         else if (ed.on && ed.i != pb)
            cy.links.emplace_back(pb, ed.i);
      }
   }

   // The cycle's island_reach marks end with it (other searches between cycles see every part); a part left out is checked at the next cycle's start (ReachLinks)
   static void ReachClear()
   {
      for (BallContact &e : g_engaged)
         for (SurfaceContact &it : e.items)
         {
            it.partOutLast = it.partOut;
            it.partOut = false;
         }
   }

   // Bodies whose surfaces are within reach of each other (plus what they can close in this cycle) interact: they form islands
   static void FindIslands(Cycle &cy)
   {
      CPROF(FindIslands);
      const ContactSettings &cs = *cy.cs;
      const int nw = (int)cy.bodies.size();
      cy.pairs.clear();
      if (cs.ballContact && nw >= 2)
         for (int i = 0; i < nw; i++)
            for (int j = i + 1; j < nw; j++)
            {
               const double gap = (cy.bodies[i].c - cy.bodies[j].c).Length() - cy.bodies[i].bs.R - cy.bodies[j].bs.R;
               const double reach = cs.ballReachMm * 1e-3 + 2.0 * (cy.bodies[i].v - cy.bodies[j].v).Length() * cy.dt;
               if (gap < reach)
                  cy.pairs.emplace_back(i, j);
            }
      cy.links.clear();
      for (int i = 0; i < nw; i++)
         if (cy.bodies[i].band != nullptr)
            for (int j = 0; j < i; j++)
               if (cy.bodies[j].band == cy.bodies[i].band)
               {
                  cy.links.emplace_back(j, i);
                  break;
               }
      if (cs.islandReach)
         ReachLinks(cy);
      else
      for (int i = 0; i < nw; i++)   // balls touching the same gate or spinner are one island: the part takes all their reactions
         if (cy.bodies[i].e != nullptr)
            for (const SurfaceContact &it : cy.bodies[i].e->items)
               if (ContactMovingPart *const hg = it.r->MovingPart())
               {
                  for (int j = 0; j < i; j++)
                     if (cy.bodies[j].e != nullptr && std::any_of(cy.bodies[j].e->items.begin(), cy.bodies[j].e->items.end(), [hg](const SurfaceContact &o) { return o.r->MovingPart() == hg; }))
                     {
                        cy.links.emplace_back(j, i);
                        break;
                     }
               }
      cy.islands.Build(nw, cy.pairs, cy.links);
   }

   // The pair states of an island (shear state kept between cycles), the substep count (the finest any member or pair needs) and the zeroed per-cycle counters
   static void SetupIsland(Cycle &cy, const ContactSystem::Islands::Island &isl)
   {
      CPROF(SetupIsland);
      const ContactSettings &cs = *cy.cs;
      const double dt = cy.dt;
      const std::vector<std::pair<int, int>> &cpairs = isl.edges;
      constexpr double kRestSpeedFree = 0.02;   // m/s: closing speed above which a pair is stepped for the accuracy of its impact, not only for stability
      constexpr double kImpactOmegaH = 0.2, kRestOmegaH = 0.8, kMaxPairStepS = 1e-4;   // omega * h limits: an impact resolved with about 30 steps, a resting contact only stable
      constexpr double kApproachMarginM = 100e-6;   // m: gap beyond the closing distance under which a not-yet-touching pair still gets its contact step

      // the shear state of each pair: found or made (two passes, because making one can move the vector)
      cy.pst.assign(cpairs.size(), nullptr);
      for (int pass = 0; pass < 2; pass++)
         for (size_t k = 0; k < cpairs.size(); k++)
         {
            const HitBall *a = cy.bodies[cpairs[k].first].pb, *b = cy.bodies[cpairs[k].second].pb;
            if (b < a)
               std::swap(a, b);
            BallPairState *found = nullptr;
            for (BallPairState &p : g_ballPairs)
               if (p.a == a && p.b == b)
                  found = &p;
            if (found == nullptr && pass == 0)
            {
               g_ballPairs.push_back(BallPairState());
               g_ballPairs.back().a = a;
               g_ballPairs.back().b = b;
            }
            if (pass == 1)
            {
               found->used = true;
               cy.pst[k] = found;
            }
         }

      int steps = 1, stepsWhy = TestDriver::ContactCost::kEnd;
      for (const int i : isl.members)
         if (cy.bodies[i].steps > steps)
         {
            steps = cy.bodies[i].steps;
            stepsWhy = cy.bodies[i].stepsWhy;
         }
      int stepsPair = 1;
      // The pair's step comes from its stiffness at the compression it will reach (the same criterion as a ball on a surface: omega * h limited), not from a constant per
      // situation: Hertz k = 2 E* sqrt(Reff delta), delta the depth it rests at or the depth an impact at the closing speed reaches (Hertz energy). An impact is resolved with
      // omega * h <= kImpactOmegaH, a resting contact only kept stable (kRestOmegaH).
      {
         const RubberContact::RubberParams &rp = cy.rpBall;
         const double Estar = rp.eMpa * 1e6 / (1.0 - rp.nu * rp.nu);
         for (const auto &p : cpairs)
         {
            const Body &A = cy.bodies[p.first], &B = cy.bodies[p.second];
            const V3 d = A.c - B.c;
            const double dist = d.Length();
            if (dist < 1e-9)
               continue;
            BallPairState &ps = *cy.pst[(size_t)(&p - cpairs.data())];
            const double vn = (A.v - B.v).Dot(d / dist), closing = std::max(0.0, -vn);
            // resting = slow now and slow a cycle ago: a hit's rebound half recedes and its deepest point falls on a cycle's start now and then; neither is a rest (a long step there gave back up to 14 % more speed than came in)
            const bool slow = std::max(std::abs(vn), ps.lastSpeed) <= kRestSpeedFree;
            ps.lastSpeed = std::abs(vn);
            // a pair that cannot touch within this cycle (the gap is more than the closing distance plus a margin for the force's own acceleration) needs no step of its own:
            // it joined the cluster only because it is within reach
            if (!ps.inContact && dist - A.bs.R - B.bs.R > closing * dt * 1.5 + kApproachMarginM)
               continue;
            const double reff = A.bs.R * B.bs.R / (A.bs.R + B.bs.R), mRed = A.bs.m * B.bs.m / (A.bs.m + B.bs.m);
            const double dImpact = std::pow(15.0 * mRed * closing * closing / (16.0 * Estar * std::sqrt(reff)), 0.4);
            const double delta = std::max({ ps.lastDelta, dImpact, 1e-6 });
            const double omega = std::sqrt(2.0 * Estar * std::sqrt(reff * delta) / mRed);
            const double omegaH = slow ? kRestOmegaH : kImpactOmegaH;
            const double h = std::min(kMaxPairStepS, std::max(cs.ballSubstepUs * 1e-6, omegaH / omega));
            if ((int)std::ceil(dt / h) > steps)
            {
               steps = (int)std::ceil(dt / h);
               stepsWhy = TestDriver::ContactCost::kPairStiffness;
            }
            stepsPair = std::max(stepsPair, (int)std::ceil(dt / h));
         }
      }
      cy.steps = steps;
      cy.stepsPair = stepsPair;
      // The gates and spinners of the island. Their step is the reduced-mass rule of a moving body (design 6.3): every ball's contact with a part is resolved for the reduced mass of
      // the ball and the part's effective mass at the contact, with the series stiffness of contact and bend, and ends at the arrival (Part::BallStep); the part's own forces
      // (pendulum, stops, switch spring) are resolved in its own substeps inside each island substep (Part::StepInIsland). A turntable disc is held along its axis by its bearing, so no
      // normal moves (its BallStep sets no limit).
      cy.islandParts.clear();
      for (const int i : isl.members)
         if (cy.bodies[i].e != nullptr)
            for (const SurfaceContact &it : cy.bodies[i].e->items)
               if (ContactMovingPart *const hg = it.partOut ? nullptr : it.r->MovingPart())   // (island_reach: a part out of this ball's reach is solved elsewhere)
               {
                  if (std::find(cy.islandParts.begin(), cy.islandParts.end(), hg) == cy.islandParts.end())
                  {
                     cy.islandParts.push_back(hg);
                     hg->SetGravity(ContactPartHost::PartGravitySI());   // (its own forces are resolved in its own substeps: Part::StepInIsland)
                  }
                  const Body &w = cy.bodies[i];
                  hg->SetBall(w.bs);
                  const int sp = (int)std::ceil(dt / hg->BallStep(w.c, w.v, dt) - 1e-9);
                  if (sp > steps)
                  {
                     steps = sp;
                     stepsWhy = TestDriver::ContactCost::kPartStep;
                  }
               }
      cy.steps = steps;
      cy.stepsWhy = stepsWhy;
      for (BallPairState *const ps : cy.pst)
      {
         ps->lastDelta = ps->maxDelta;
         ps->maxDelta = ps->maxFn = 0.0;
      }
      for (const int i : isl.members)
      {
         Body &w = cy.bodies[i];
         w.omega2 = w.deepest = w.vnLast = 0.0;
         w.tripped = false;
      }
   }

   // Forces of the ball's own surfaces (rubbers, the playfield, the disc) at the current state; the disc gets the reaction
   static void SurfaceForces(Cycle &cy, Body &w, const double h, const int s)
   {
      w.F = V3();
      w.Mb = V3();
      w.vnLast = 0.0;
      if (w.e == nullptr)
         return;   // a free ball: only the pair forces act on it
      CPROF(SurfaceForces);
      BallContact &e = *w.e;
      size_t itemIdx = 0;
      for (SurfaceContact &it : e.items)
      {
         const size_t k = itemIdx++;
         if (it.partOut || it.NoFaces(w.c))   // a surface of faces with none within reach of the ball this cycle (or a part out of its reach, island_reach): no contact
         {
            it.st.Reset();
            it.touching = false;
            it.fnLast = 0.0;
            continue;
         }
         TubeProbe pr;
         if (w.probesReady)
            pr = w.probes[k];
         else
         {
            CPROF(SurfProbe);
            pr = it.r->Probe(w.bs, it, w.c);
         }
         if (it.r->MovingPart() != nullptr)
         {
            // a gate or spinner: its own contact (series bend, tilt force, its Hit event near the start angle, the reaction on its body)
            const ContactResult r = it.r->EvaluateOwn(it, w.bs, w.c, w.v, w.w, h);
            it.touching = r.inContact;
            it.fnLast = r.inContact ? r.fn : 0.0;
            if (r.inContact)
            {
               w.vnLast = std::max(w.vnLast, std::abs(w.v.Dot(r.normal)));
               w.deepest = std::max(w.deepest, r.delta);
               e.peakFn = std::max(e.peakFn, r.fn);
               w.F = w.F + r.force;
               w.Mb = w.Mb + r.torqueOnBall;
            }
            continue;
         }
         V3 vsurf, wsurf;   // velocity of the surface at the contact and its angular velocity (zero for a static surface: the moving parts took the branch above)
         if (pr.overlap)
            it.r->VelocityAt(pr, vsurf, wsurf);
         ContactResult r;
         {
            CPROF(SurfEval);
            r = EvaluateTube(it.rp, w.bs, pr, it.hLayer, w.v, w.w, h, it.st, vsurf, wsurf);
         }
         it.fnLast = r.inContact ? r.fn : 0.0;
         if (r.inContact != it.touching)
         {
            it.touching = r.inContact;
            if (r.inContact)
            {
               const float vn = (float)(-w.v.Dot(r.normal) / VS);   // the approach speed, engine velocity units (the stock hit events take it along the face normal)
               const HitObject *const ho = it.r->HitEventObject();
               if ((vn >= cy.hitFloor || (ho != nullptr && ho->m_fe && vn >= ho->m_threshold)) && it.r->HitEventWanted(r.normal))   // (below the floor, only a part with a hit event that asks for it)
                  cy.partHits.push_back({ it.r, w.pb, vn, Vertex3Ds((float)r.normal.x, (float)r.normal.y, (float)r.normal.z) });
            }
         }
         if (r.inContact)
         {
            CPROF(SurfReact);
            it.r->ApplyReaction(pr, r, w.bs.R, h);   // a moving surface gets the equal and opposite
         }
         if (r.inContact)
         {
            w.vnLast = std::max(w.vnLast, std::abs(w.v.x * r.normal.x + w.v.y * r.normal.y + w.v.z * r.normal.z));
            if (s == 0 && r.fn > 0.0 && r.delta > 1e-9)
            {
               // contact stiffness (Hertz 1.5 Fn / delta, doubled for the layer's stiffening) and the sideways spring, on the ball's mass
               const double E = it.rp.eMpa * 1e6, G = E / (2.0 * (1.0 + it.rp.nu)), Gs = G / (2.0 - it.rp.nu);
               const double kt = 8.0 * Gs * r.patchA * it.rp.shearScale, mt = 1.0 / (1.0 / w.bs.m + w.bs.R * w.bs.R / w.bs.Inertia());
               w.omega2 += 3.0 * r.fn / r.delta / w.bs.m + kt / mt;
            }
            e.peakFn = std::max(e.peakFn, r.fn);
            w.deepest = std::max(w.deepest, r.delta);
            w.F = w.F + r.force;
            w.Mb = w.Mb + r.torqueOnBall;
         }
      }
   }

   // The Option A ball-ball magnet force (Magnets::BallBallForce, solved once per physics step), re-aimed at every substep: each pair part turns with the CURRENT
   // line between its two balls, each magnet part with its ball round the magnet axis, so a ball moving round its partner gets no work from a held direction
   // (the force held for the whole 1 ms step does F r (w dt)^2 / 2 per step: orbiting balls spun up, s19). A ball outside this island counts where the step began.
   // Sized as the held increment: the force on the steel mass the magnets use, times the body's mass.
   static void MagnetPartForces(Cycle &cy, const ContactSystem::Islands::Island &isl, const int s)
   {
      using M3 = MagnetBalls::V3;
      if (s == 0 || cy.magIsl != &isl)   // resolve this island's parts once: each part's balls become body indices
      {
         cy.magIsl = &isl; cy.magParts.clear();
         const auto bodyOf = [&](const HitBall *x) { for (const int k : isl.members) if (cy.bodies[k].pb == x) return k; return -1; };
         for (const int i : isl.members)
         {
            const Body &B = cy.bodies[i]; double massKg = 0; const std::vector<Magnets::BallBallPart> *pp;
            if (!B.magParts || (pp = Magnets::BallBallParts(B.pb, massKg)) == nullptr || massKg <= 0)
               continue;
            for (const Magnets::BallBallPart &p : *pp)
            {
               Cycle::MagPartRT r; r.body = i; r.scale = B.bs.m / massKg; r.pair = p.b != nullptr;
               r.ia = bodyOf(p.a); r.ib = r.pair ? bodyOf(p.b) : -1;
               r.F = { p.F[0], p.F[1], p.F[2] }; r.u0 = { p.u0[0], p.u0[1], p.u0[2] }; r.n = { p.n[0], p.n[1], p.n[2] }; r.o = { p.o[0], p.o[1], p.o[2] };
               r.r0 = { p.r0[0], p.r0[1], p.r0[2] }; r.pa0 = { p.pa0[0], p.pa0[1], p.pa0[2] }; r.pb0 = { p.pb0[0], p.pb0[1], p.pb0[2] };
               cy.magParts.push_back(r);
            }
         }
      }
      for (const Cycle::MagPartRT &r : cy.magParts)
      {
         const auto at = [&](int k, const M3 &p0) { if (k < 0) return p0; const V3 &c = cy.bodies[k].c; return M3 { c.x, c.y, c.z }; };
         const M3 pa = at(r.ia, r.pa0);
         M3 f = r.F;
         if (r.pair)   // turned with the current line between its balls
         {
            const M3 pb = at(r.ib, r.pb0), d = { pa.x - pb.x, pa.y - pb.y, pa.z - pb.z }; const double l = std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
            if (l > 0) f = MagnetBalls::ReaimPair(r.F, r.u0, { d.x / l, d.y / l, d.z / l });
         }
         else if (r.n.x != 0 || r.n.y != 0 || r.n.z != 0)   // turned round the magnet axis as its ball moved round it: cos and sin from the radial directions
         {
            M3 d = { pa.x - r.o.x, pa.y - r.o.y, pa.z - r.o.z }; const double dn = d.x * r.n.x + d.y * r.n.y + d.z * r.n.z;
            d = { d.x - dn * r.n.x, d.y - dn * r.n.y, d.z - dn * r.n.z }; const double l = std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
            if (l > 1e-6)
            {
               const M3 r1 = { d.x / l, d.y / l, d.z / l }, x = { r.r0.y * r1.z - r.r0.z * r1.y, r.r0.z * r1.x - r.r0.x * r1.z, r.r0.x * r1.y - r.r0.y * r1.x };
               f = MagnetBalls::ReaimMagnetCS(r.F, r.n, r.r0.x * r1.x + r.r0.y * r1.y + r.r0.z * r1.z, x.x * r.n.x + x.y * r.n.y + x.z * r.n.z);
            }
         }
         Body &B = cy.bodies[r.body];
         B.F = B.F + V3(f.x, f.y, f.z) * r.scale;
      }
   }

   // Ball to ball: the same contact law, equal and opposite on the two balls
   static void PairForces(Cycle &cy, const ContactSystem::Islands::Island &isl, const double h, const int s)
   {
      MagnetPartForces(cy, isl, s);
      CPROF(PairForces);
      const std::vector<std::pair<int, int>> &cpairs = isl.edges;
      for (size_t k = 0; k < cpairs.size(); k++)
      {
         Body &A = cy.bodies[cpairs[k].first], &B = cy.bodies[cpairs[k].second];
         BallPairState &ps = *cy.pst[k];
         const V3 d = A.c - B.c;
         const double dist = d.Length();
         ps.zero = std::max(0.0, std::min(ps.zero, A.bs.R + B.bs.R - dist));   // (backing out, with nothing stored: the zero follows)
         const double delta = A.bs.R + B.bs.R - dist - ps.zero;
         if (delta <= 0.0 || dist < 1e-9)
         {
            ps.st.Reset();
            ps.inContact = false;
            continue;
         }
         const V3 n = d / dist;   // from B to A
         const V3 rbB = n * (B.bs.R - delta / 2.0);   // the contact point relative to B's centre
         const V3 vsB = B.v + B.w.Cross(rbB);
         const double reff = A.bs.R * B.bs.R / (A.bs.R + B.bs.R);
         const CoreResult core = ContactCore(cy.rpBall, A.bs, n, delta, reff, 1.0, vsB, A.v, A.w, h, ps.st, B.w);
         const V3 couple = core.torque - core.rB.Cross(core.force);
         A.F = A.F + core.force;
         A.Mb = A.Mb + core.torque;
         B.F = B.F - core.force;
         B.Mb = B.Mb + rbB.Cross(core.force * -1.0) - couple;
         if (cy.cs->ballRollingResistance > 0.0 && core.fn > 0.0)   // rolling resistance, as every other contact has: the balls' relative spin in the tangent plane is resisted
         {
            const V3 wr = A.w - B.w;
            const V3 wt = wr - n * wr.Dot(n);
            const double wl = wt.Length();
            if (wl > 1e-9)
            {
               const V3 T = wt * (-cy.cs->ballRollingResistance * reff * core.fn * std::tanh(wl / 2.0) / wl);
               A.Mb = A.Mb + T;
               B.Mb = B.Mb - T;
            }
         }
         if (!ps.inContact)
         {
            ps.inContact = true;
            const double closing = -(A.v - B.v).Dot(n);
            if (closing / VS > 0.25)   // as the stock collision: only impacts, not resting contact, reach the script
               cy.hitEvents.push_back({ A.pb, B.pb, (float)(closing / VS) });
         }
         ps.maxDelta = std::max(ps.maxDelta, delta);
         ps.maxFn = std::max(ps.maxFn, core.fn);
         if (A.e) A.e->peakFn = std::max(A.e->peakFn, core.fn);
         if (B.e) B.e->peakFn = std::max(B.e->peakFn, core.fn);
         if (s == 0 && core.fn > 0.0 && delta > 1e-9)
         {
            const double mRed = A.bs.m * B.bs.m / (A.bs.m + B.bs.m);
            const double o2 = 3.0 * core.fn / delta / mRed;
            A.omega2 += o2;
            B.omega2 += o2;
         }
      }
   }

   // Semi-implicit Euler for every member: the summed forces, then the sources; returns true if a member left the model's range
   static bool Integrate(Cycle &cy, const ContactSystem::Islands::Island &isl, const double h, const int s)
   {
      CPROF(Integrate);
      const double aScale = h / PHYSICS_STEPTIME_S;
      bool tripped = false;
      for (const int i : isl.members)
      {
         Body &w = cy.bodies[i];
         w.v = w.v + w.F * (1.0 / w.bs.m) * h;
         {
            CPROF(Sources);
            Sources::Add(w, aScale);
         }
         w.w = w.w + w.Mb * (1.0 / w.bs.Inertia()) * h;
         w.c = w.c + w.v * h;
         // a non-finite state, a runaway speed or a ball sunk deep into a surface means the step left the model's range (a ball wedged in a loop tighter than itself, a ball placed
         // inside a solid). The step is taken again in two; one that still trips stops the ball (StepFailed).
         if (!std::isfinite(w.c.x + w.c.y + w.c.z + w.v.x + w.v.y + w.v.z + w.w.x + w.w.y + w.w.z) || w.v.Length() > w.vCap || w.deepest > 0.8 * w.bs.R)
         {
            w.tripped = true;
            w.tripWhy = !std::isfinite(w.c.x + w.c.y + w.c.z + w.v.x + w.v.y + w.v.z + w.w.x + w.w.y + w.w.z) ? "not finite" : w.v.Length() > w.vCap ? "speed" : "depth";
            w.tripStep = s;
            tripped = true;
         }
      }
      return tripped;
   }

   // ---- Arrival-time substeps (islands of balls and static surfaces) ----
   // A step never carries a ball into a surface. At the start of every substep each ball is probed against the surfaces of its contact; a surface that is not touching yet bounds the step
   // by the time the ball needs to get within kArriveM of it (distance, closing speed, and the acceleration it can pick up), and a touching one by its own stiffness at the depth the
   // contact will reach (omega * h limited, Hertz). A ball with nothing near takes the rest of the cycle as one step. The ball pairs keep the step their own rule gives (cy.stepsPair).

   // Islands the arrival loop (with the rewind) can step: balls, static surfaces and moving parts (gates, spinners, turntable discs: they step inside the loop and are saved
   // with the snapshot, UnifiedContactDesign.txt 11 point 7, 12.1). A bat is on the group step, never here. A band keeps the uniform path until it can be saved and restored (12.3).
   static bool IslandRewindable(const Cycle &cy, const ContactSystem::Islands::Island &isl)
   {
      for (const int i : isl.members)
      {
         const Body &w = cy.bodies[i];
         if (w.band != nullptr)
            return false;
      }
      return true;
   }

   // t: the time of this substep into the cycle (s). At the cycle's first substep every surface is probed; after that a static surface whose gap, less what the ball has moved
   // since its last probe, is still more than the ball can close in the rest of the cycle (its speed and the acceleration bound of the arrival rule, NextStep) is not probed: it
   // cannot touch the ball in this cycle, and NextStep would have passed it over (d > vc h + a h^2 for every h up to the rest of the cycle). Moving parts are always probed.
   static void ProbeBody(Cycle &cy, Body &w, const double t)
   {
      CPROF(ProbeBody);
      constexpr double kTouchM = 1.0e-6;   // as NextStep: within this gap a surface counts as touching
      w.probesReady = false;
      if (w.e == nullptr)
         return;
      w.probes.resize(w.e->items.size());
      const bool first = t <= 0.0;
      const double rem = std::max(cy.dt - t, 0.0);
      const double aBall = w.extBaseVS.Length() / PHYSICS_STEPTIME_S + w.magAcc + w.F.Length() / w.bs.m;
      const double reach = w.v.Length() * rem + aBall * rem * rem;
      size_t k = 0;
      for (SurfaceContact &it : w.e->items)
      {
         TubeProbe &pr = w.probes[k++];
         if (it.NoFaces(w.c))   // no face within reach this cycle: far
         {
            pr = TubeProbe();
            pr.dist = -1.0;
            it.probeT = -1.0;
            continue;
         }
         if (first)
            it.probeT = -1.0;
         else if (it.probeT >= 0.0 && it.r->MovingPart() == nullptr && it.probeGap - (w.c - it.probeC).Length() - kTouchM > reach)
         {
            pr = TubeProbe();   // cannot be reached in the rest of the cycle: far
            pr.dist = -1.0;
            it.hcStep = 0.0;    // (as NextStep does for a surface that is not touching)
            continue;
         }
         if (ContactProf::On())
         {
            const uint64_t t0 = __rdtsc();
            pr = it.r->Probe(w.bs, it, w.c);
            if (ContactProf::ProbeTime(it.r, __rdtsc() - t0))
               ContactProf::ProbeName(it.r, it.r->Name());
         }
         else
            pr = it.r->Probe(w.bs, it, w.c);
         it.probeGap = pr.overlap ? -pr.delta : pr.dist - (w.bs.R + pr.rho);
         it.probeC = w.c;
         it.probeT = t;
         // VPX_PROBE_CHECK=1: the short face list must give the same nearest face as the whole Verlet list wherever it matters (within the ball's radius plus 3 mm); mismatches are counted in the profile notes
         static const bool check = getenv("VPX_PROBE_CHECK") != nullptr;
         if (check && it.candV.size() > it.cand.size() && it.candR >= 0.0)
         {
            SurfaceContact full = it;
            full.cand = it.candV;
            full.tightC = w.c;
            const TubeProbe p2 = it.r->Probe(w.bs, full, w.c);
            const double nearM = w.bs.R + 0.003;
            const bool relevant = std::min(pr.dist, p2.dist) <= nearM;
            if (relevant && (std::abs(pr.dist - p2.dist) > 1e-12 || pr.seg != p2.seg))
               ContactProf::Note(std::string("PROBE MISMATCH ") + it.r->Name());
            else if (relevant)
               ContactProf::Note("probe check ok");
         }
      }
      w.probesReady = true;
   }

   // The next substep of an island: at most `remaining`
   // why: the rule that set the step (TestDriver::ContactCost::Why, for the test driver's cost split)
   static double NextStep(Cycle &cy, const ContactSystem::Islands::Island &isl, const double remaining, int &arrivals, int &why)
   {
      using CC = TestDriver::ContactCost;
      CPROF(NextStep);
      constexpr double kTouchM = 1.0e-6;          // within this gap a surface counts as touching
      constexpr double kArriveM = 0.5e-6;         // an arrival step ends this far above the surface
      constexpr double kMinStepS = 0.2e-6, kMaxContactStepS = 100e-6;
      const double kImpactOmegaH = 0.2 * ContactOmegaHScale();
      constexpr double kRestOmegaH = 0.8;
      constexpr double kRestSpeed = 0.02;   // as the ball pairs' rule
      // resting = slow AND the forces on the ball balance: what is left over (m times its acceleration over the last substep) is at most this share of the contact's own force.
      // At rest the contact carries everything that presses the ball into it (weight, magnets, balls on top); at the deepest point of an impact it far exceeds the rest.
      constexpr double kBalance = 0.2;   // measured 2026-10-21 (6 balls on the magnet): the strongest contact of a resting ball shows net force / contact force 0.005 median and 0.01 at most; in the bursts after a landing 99 % of cycles are under 0.08. 0.2 is 20x the resting maximum
      double h = remaining;
      why = CC::kEnd;
      if (cy.stepsPair > 1 && cy.dt / cy.stepsPair < h)
      {
         h = cy.dt / cy.stepsPair;
         why = CC::kPairStiffness;
      }
      for (const int i : isl.members)
      {
         Body &w = cy.bodies[i];
         if (w.e == nullptr || !w.probesReady)
            continue;
         // the acceleration the ball can pick up: gravity and the nudge, a magnet's pull, and what the contact forces of the last substep gave it
         const double fLeft = w.hS > 0.0 ? (w.v - w.vS).Length() / w.hS * w.bs.m : -1.0;   // N, the ball's net force over the last substep (< 0 = unknown)
         const double aBallBound = w.extBaseVS.Length() / PHYSICS_STEPTIME_S + w.magAcc + w.F.Length() / w.bs.m;
         size_t k = 0;
         for (SurfaceContact &it : w.e->items)
         {
            const TubeProbe &pr = w.probes[k++];
            if (it.partOut || (!pr.overlap && pr.dist < 0.0))
               continue;   // far (or a part out of its reach, island_reach)
            // a moving part closes on the ball too: its present speed and what its drive can add (SurfaceSpeedBound, linear in the time) join the ball's closing speed and acceleration
            ContactMovingPart *const mp = it.r->MovingPart();
            double vSurf = 0.0, aSurf = 0.0;
            if (mp != nullptr)
            {
               mp->SetBall(w.bs);
               vSurf = mp->SurfaceSpeedBound(0.0);
               aSurf = (mp->SurfaceSpeedBound(cy.dt) - vSurf) / cy.dt;
            }
            const double aBound = aBallBound + aSurf;
            const double gap = pr.overlap ? -pr.delta : pr.dist - (w.bs.R + pr.rho);
            const double vn = -w.v.Dot(pr.normal);
            // The probe sees the nearest part of the surface only: every other part the ball can reach in this step bounds it too (no ball is carried into a part no probe looked at).
            // Distance cannot shrink faster than the ball moves, so a nearest gap beyond that travel rules out the rest.
            if (mp == nullptr && gap <= w.v.Length() * h + aBound * h * h)
            {
               const double ha = it.r->ArrivalS(w.bs, it, pr.seg, w.c, w.v, aBound, h, kTouchM, kArriveM);
               if (ha < h)
               {
                  h = std::max(ha, kMinStepS);
                  why = ha < kMinStepS ? CC::kMinimum : CC::kArrival;
                  arrivals++;
               }
            }
            if (gap > kTouchM)
            {
               it.hcStep = 0.0;
               const double d = gap - kArriveM, vc = std::max(vn, 0.0) + vSurf;
               if (d > vc * h + aBound * h * h)   // (semi-implicit Euler moves the ball by the new velocity: v h + a h^2 for a constant acceleration)
                  continue;   // cannot get there within the step
               { static const bool noSplit = getenv("VPX_ARRIVAL_NOSPLIT") != nullptr; if (noSplit) continue; }   // diagnostic: no arrival split, the rewind has to catch every late contact
               double t = 1e9;
               if (aBound > 1e-9)
                  t = (-vc + std::sqrt(vc * vc + 4.0 * aBound * d)) / (2.0 * aBound);
               else if (vc > 1e-12)
                  t = d / vc;
               if (t < h)
               {
                  h = std::max(t, kMinStepS);
                  why = t < kMinStepS ? CC::kMinimum : CC::kArrival;
                  arrivals++;
               }
            }
            else
            {
               static const bool noLand = getenv("VPX_ARRIVAL_NOSPLIT") != nullptr;
               const double vl = std::max(vn, 0.0) + vSurf;
               if (!noLand && !pr.overlap && vl > 1e-9 && gap > 0.0)   // inside the touch gap but not touching yet: land on the surface exactly, so the force starts at the touch
               {
                  const double tl = aBound > 1e-9 ? (-vl + std::sqrt(vl * vl + 4.0 * aBound * gap)) / (2.0 * aBound) : gap / vl;
                  if (tl < h)
                  {
                     h = std::max(tl, kMinStepS);
                     why = tl < kMinStepS ? CC::kMinimum : CC::kArrival;
                     arrivals++;
                  }
               }
               // touching: the contact's stiffness at the compression it will reach (resting, or an impact at the closing speed), k = 2 E* sqrt(R delta)
               const double Estar = it.rp.eMpa * 1e6 / (1.0 - it.rp.nu * it.rp.nu);
               // (the speed of the contact motion, either way: a compressed contact leaving is as stiff as one arriving, and its energy is still to be given back)
               const double sp = std::fabs(vn), dNow = pr.overlap ? pr.delta : 0.0;
               const bool balanced = fLeft >= 0.0 && it.fnLast > 0.0 && fLeft <= kBalance * it.fnLast;
               double hc;
               if (it.hcStep > 0.0 && sp <= it.hcSpeed * 1.25 + 0.005 && dNow <= it.hcDelta * 1.1 + 1e-7 && !(sp <= kRestSpeed && balanced && !it.hcRest))
                  hc = it.hcStep;   // covered by the step worked out for this contact (a bigger speed or compression recomputes it, and so does a contact that has come to rest)
               else
               {
                  const double dImpact = sp > kRestSpeed ? std::pow(15.0 * w.bs.m * sp * sp / (16.0 * Estar * std::sqrt(w.bs.R)), 0.4) : 0.0;
                  // the depth the ball's weight (and what pushes it) presses it to: a contact that is not moving and not deeper than a few times that is quasi-static and only has to be
                  // stable (omega * h <= 0.8); anything else, a speed through zero at the deepest point of an impact included, is stepped to resolve its dynamics (omega * h <= 0.2)
                  const double aG = std::max(w.extBaseVS.Length() / PHYSICS_STEPTIME_S, 1.0);
                  const double dRest = std::pow(w.bs.m * aG / ((4.0 / 3.0) * Estar * std::sqrt(w.bs.R)), 2.0 / 3.0);
                  const bool quasiStatic = sp <= kRestSpeed && (balanced || dNow <= 3.0 * dRest);
                  const double delta = std::max({ dNow, dImpact, quasiStatic ? dRest : 0.0, 1e-6 });
                  const double omega = std::sqrt(2.0 * Estar * std::sqrt(w.bs.R * delta) / w.bs.m);
                  hc = std::min(kMaxContactStepS, (quasiStatic ? kRestOmegaH : kImpactOmegaH) / omega);
                  it.hcStep = hc;
                  it.hcRest = quasiStatic;
                  it.hcSpeed = std::max(sp, kRestSpeed);
                  it.hcDelta = delta;
               }
               if (const RubberObjectSettings *const set = it.r->Settings())   // rolling along a surface: the distance covered per step is limited too
                  if (set->rollTravelUm > 0.0)
                  {
                     const V3 vt = w.v - pr.normal * w.v.Dot(pr.normal);
                     hc = std::min(hc, set->rollTravelUm * 1e-6 / std::max(vt.Length(), 1e-6));
                  }
               int hcWhy = CC::kSurfaceStiffness;
               if (mp != nullptr)   // the reduced mass of the ball and the part, the part's own stiffness and the surface speed it can reach (Part::BallStep)
               {
                  const double hp = mp->BallStep(w.c, w.v, h);
                  if (hp < hc)
                  {
                     hc = hp;
                     hcWhy = CC::kPartStep;
                  }
               }
               if (std::max(hc, kMinStepS) < h)
               {
                  h = std::max(hc, kMinStepS);
                  why = hc < kMinStepS ? CC::kMinimum : hcWhy;
               }
            }
         }
      }
      if (h < kMinStepS && kMinStepS < remaining)
         why = CC::kMinimum;
      else if (remaining <= std::max(h, kMinStepS))
         why = CC::kEnd;
      return std::min(remaining, std::max(h, kMinStepS));
   }

   // The state a substep changes, kept at its start so the substep can be taken again (a contact that started too deep, or a trip): the balls, the shear springs and touching flags
   // of their contacts, the ball pairs, the events queued. The gap and closing speed of each surface at that point say whether the next substep started a contact late.
   struct RwSnap
   {
      struct Bd { V3 c, v, w; double omega2, deepest, vnLast, peakFn; size_t itemBase; size_t nItems; };
      struct It { V3 xi; bool touching; RubberContact::RollAngle ra; };
      struct Pr { V3 xi; bool inContact; double maxDelta, maxFn, zero; RubberContact::RollAngle ra; };
      std::vector<Bd> b;
      std::vector<It> it;
      std::vector<Pr> p;
      std::vector<double> gap, vn;   // per surface of each ball (flat, in the order of b): the gap and the closing speed at the snapshot (1e9 / 0 = far)
      std::vector<double> pairGap;
      std::vector<ContactMovingPart::State> parts;   // the island's moving parts (cy.islandParts), saved with the balls
      size_t partHits = 0, hitEvents = 0;
      double t = 0.0, h = 0.0;
      int s = 0;
      bool approaching = false;   // some surface or pair was within reach but not touching at the snapshot: the substep must be checked at its end
   };

   static void RwTake(Cycle &cy, const ContactSystem::Islands::Island &isl, RwSnap &sn, const double t, const int s)
   {
      sn.b.clear();
      sn.it.clear();
      sn.gap.clear();
      sn.vn.clear();
      for (const int i : isl.members)
      {
         Body &w = cy.bodies[i];
         RwSnap::Bd d = { w.c, w.v, w.w, w.omega2, w.deepest, w.vnLast, w.e != nullptr ? w.e->peakFn : 0.0, sn.it.size(), 0 };
         if (w.e != nullptr)
         {
            size_t k = 0;
            for (SurfaceContact &it : w.e->items)
            {
               sn.it.push_back({ it.st.xi, it.touching, it.st.ra });
               const TubeProbe &pr = w.probes[k++];
               if (!pr.overlap && pr.dist < 0.0)
               {
                  sn.gap.push_back(1e9);
                  sn.vn.push_back(0.0);
               }
               else
               {
                  sn.gap.push_back(pr.overlap ? -pr.delta : pr.dist - (w.bs.R + pr.rho));
                  sn.vn.push_back(-w.v.Dot(pr.normal));
               }
               d.nItems++;
            }
         }
         sn.b.push_back(d);
      }
      sn.p.clear();
      sn.pairGap.clear();
      for (size_t k = 0; k < isl.edges.size(); k++)
      {
         const BallPairState &ps = *cy.pst[k];
         sn.p.push_back({ ps.st.xi, ps.inContact, ps.maxDelta, ps.maxFn, ps.zero, ps.st.ra });
         const Body &A = cy.bodies[isl.edges[k].first], &B = cy.bodies[isl.edges[k].second];
         sn.pairGap.push_back((A.c - B.c).Length() - A.bs.R - B.bs.R + ps.zero);
      }
      sn.partHits = cy.partHits.size();
      sn.hitEvents = cy.hitEvents.size();
      sn.parts.resize(cy.islandParts.size());
      for (size_t k = 0; k < cy.islandParts.size(); k++)
         cy.islandParts[k]->SaveState(sn.parts[k]);
      sn.t = t;
      sn.s = s;
      sn.approaching = false;
      for (const double g : sn.gap)
         if (g > 1.0e-6 && g < 1.0e8)
            sn.approaching = true;
      for (const double g : sn.pairGap)
         if (g > 1.0e-6)
            sn.approaching = true;
   }

   static void RwRestore(Cycle &cy, const ContactSystem::Islands::Island &isl, const RwSnap &sn)
   {
      for (const int i : isl.members)
         cy.bodies[i].hS = 0.0;   // the substep is taken again: no net acceleration known for it
      size_t bi = 0;
      for (const int i : isl.members)
      {
         Body &w = cy.bodies[i];
         const RwSnap::Bd &d = sn.b[bi++];
         w.c = d.c;
         w.v = d.v;
         w.w = d.w;
         w.omega2 = d.omega2;
         w.deepest = d.deepest;
         w.vnLast = d.vnLast;
         w.tripped = false;
         if (w.e != nullptr)
         {
            w.e->peakFn = d.peakFn;
            size_t k = d.itemBase;
            for (SurfaceContact &it : w.e->items)
            {
               it.st.xi = sn.it[k].xi;
               it.st.ra = sn.it[k].ra;
               it.touching = sn.it[k].touching;
               k++;
            }
         }
      }
      for (size_t k = 0; k < isl.edges.size(); k++)
      {
         BallPairState &ps = *cy.pst[k];
         ps.st.xi = sn.p[k].xi;
         ps.st.ra = sn.p[k].ra;
         ps.inContact = sn.p[k].inContact;
         ps.maxDelta = sn.p[k].maxDelta;
         ps.maxFn = sn.p[k].maxFn;
         ps.zero = sn.p[k].zero;
      }
      cy.partHits.resize(sn.partHits);
      cy.hitEvents.resize(sn.hitEvents);
      for (size_t k = 0; k < cy.islandParts.size() && k < sn.parts.size(); k++)
         cy.islandParts[k]->RestoreState(sn.parts[k]);
   }

   // Did the substep just taken start a contact deeper than it should (the surface was not touching at its start)? The allowed first depth is the larger of the setting and a quarter of
   // the compression the impact will reach. On a hit: `f` is the share of the substep before the touch, and `what` / `depth` / `prevGap` say which contact (for the log).
   static bool RwLate(Cycle &cy, const ContactSystem::Islands::Island &isl, const RwSnap &sn, const double allowAbs, double &f, std::string &what, double &depthOut, double &prevGapOut, double &vnOut)
   {
      constexpr double kTouchM = 1.0e-6, kRestSpeed = 0.02;
      size_t bi = 0;
      for (const int i : isl.members)
      {
         const Body &w = cy.bodies[i];
         const RwSnap::Bd &d = sn.b[bi++];
         if (w.e == nullptr || !w.probesReady)
            continue;
         size_t k = 0;
         for (const SurfaceContact &it : w.e->items)
         {
            const TubeProbe &pr = w.probes[k];
            const double prevGap = sn.gap[d.itemBase + k];
            // the closing speed of ball and surface: a moving part adds its own speed and what its drive can add in the substep (SurfaceSpeedBound)
            const double vnPrev = std::max(sn.vn[d.itemBase + k], 0.0) + (!it.partOut && it.r->MovingPart() != nullptr ? it.r->MovingPart()->SurfaceSpeedBound(sn.h) : 0.0);
            k++;
            if (!pr.overlap)
               continue;
            if (prevGap <= kTouchM)
            {
               // a contact that was already touching: the depth cannot grow by more than the closing speed carries it in the substep (a margin of 3 for the contact's own motion), or the step was too long
               const double prevDepth = std::max(0.0, -prevGap);
               if (pr.delta > prevDepth + 3.0 * vnPrev * sn.h + 2.0 * allowAbs)
               {
                  f = 0.55;
                  what = std::string("deepening surface ") + it.r->Name();
                  depthOut = pr.delta;
                  prevGapOut = prevGap;
                  vnOut = vnPrev;
                  return true;
               }
               continue;
            }
            const double depth = pr.delta;
            const double Estar = it.rp.eMpa * 1e6 / (1.0 - it.rp.nu * it.rp.nu);
            const double dImp = vnPrev > kRestSpeed ? std::pow(15.0 * w.bs.m * vnPrev * vnPrev / (16.0 * Estar * std::sqrt(w.bs.R)), 0.4) : 0.0;
            if (depth > std::max(allowAbs, 0.25 * dImp))
            {
               f = prevGap / (prevGap + depth);
               what = std::string("surface ") + it.r->Name();
               depthOut = depth;
               prevGapOut = prevGap;
               vnOut = vnPrev;
               return true;
            }
         }
      }
      for (size_t k = 0; k < isl.edges.size(); k++)
      {
         const Body &A = cy.bodies[isl.edges[k].first], &B = cy.bodies[isl.edges[k].second];
         const double gapNow = (A.c - B.c).Length() - A.bs.R - B.bs.R + cy.pst[k]->zero;
         if (gapNow >= -allowAbs || sn.pairGap[k] <= kTouchM)
            continue;
         // the closing speed of the pair at the snapshot
         size_t ia = 0, ib = 0, bi2 = 0;
         for (const int i : isl.members)
         {
            if (i == isl.edges[k].first)
               ia = bi2;
            if (i == isl.edges[k].second)
               ib = bi2;
            bi2++;
         }
         const V3 dd = sn.b[ia].c - sn.b[ib].c;
         const double dist = dd.Length();
         const double closing = dist > 1e-12 ? std::max(0.0, -(sn.b[ia].v - sn.b[ib].v).Dot(dd / dist)) : 0.0;
         const double Estar = cy.rpBall.eMpa * 1e6 / (1.0 - cy.rpBall.nu * cy.rpBall.nu);
         const double reff = A.bs.R * B.bs.R / (A.bs.R + B.bs.R), mRed = A.bs.m * B.bs.m / (A.bs.m + B.bs.m);
         const double dImp = closing > kRestSpeed ? std::pow(15.0 * mRed * closing * closing / (16.0 * Estar * std::sqrt(reff)), 0.4) : 0.0;
         if (-gapNow > std::max(allowAbs, 0.25 * dImp))
         {
            f = sn.pairGap[k] / (sn.pairGap[k] - gapNow);
            what = "ball pair " + std::to_string(isl.edges[k].first) + "," + std::to_string(isl.edges[k].second);
            depthOut = -gapNow;
            prevGapOut = sn.pairGap[k];
            vnOut = closing;
            return true;
         }
      }
      return false;
   }

   static bool SubstepsArrival(Cycle &cy, const ContactSystem::Islands::Island &isl)
   {
      CPROF(Substeps);
      static RwSnap sn;
      static uint32_t lastLogMs = 0;
      bool tripped = false, havePrev = false;
      double t = 0.0, hCap = 0.0;
      int s = 0, arrivals = 0, rewinds = 0;
      constexpr int kMaxSubsteps = 20000, kMaxRewinds = 16;
      constexpr double kMinStepS = 0.2e-6;
      const double allowAbs = cy.cs->arrivalDepthUm * 1e-6;
      auto log = [&](const char *kind, const std::string &what, const double depth, const double prevGap, const double vn, const double hPrev, const double hNew)
      {
         if (g_pplayer->m_time_msec - lastLogMs < 500)
            return;
         lastLogMs = g_pplayer->m_time_msec;
         PLOGW << "Rewind (" << kind << ") t=" << g_pplayer->m_time_msec << " ms, island of " << isl.members.size() << ": " << what << ", depth " << depth * 1e6 << " um after a substep that started " << prevGap * 1e6
               << " um away, closing " << vn << " m/s, substep " << hPrev * 1e6 << " -> " << hNew * 1e6 << " us (" << rewinds << " in this cycle)";
      };
      for (;;)
      {
         const bool finished = !(t < cy.dt * (1.0 - 1e-9) && !tripped && s < kMaxSubsteps);
         if (finished && !(havePrev && sn.approaching && !tripped && rewinds < kMaxRewinds))
            break;   // (the last substep of the cycle is checked too when a surface was approaching at its start)
         {
            const TestDriver::KindTimer kt(s_cc, TestDriver::ContactCost::kSurfaces);
            for (const int i : isl.members)
               ProbeBody(cy, cy.bodies[i], t);
         }
         if (havePrev && rewinds < kMaxRewinds)
         {
            double f = 1.0, depth = 0.0, prevGap = 0.0, vn = 0.0;
            std::string what;
            if (RwLate(cy, isl, sn, allowAbs, f, what, depth, prevGap, vn))
            {
               const double hNew = std::max(kMinStepS, 0.9 * f * sn.h);
               rewinds++;
               log("late contact", what, depth, prevGap, vn, sn.h, hNew);
               RwRestore(cy, isl, sn);
               t = sn.t;
               s = sn.s;
               hCap = hNew;
               havePrev = false;
               continue;
            }
         }
         if (finished)
            break;
         RwTake(cy, isl, sn, t, s);
         havePrev = true;
         int why = TestDriver::ContactCost::kEnd;
         double h = (s == kMaxSubsteps - 1) ? cy.dt - t : NextStep(cy, isl, cy.dt - t, arrivals, why);
         for (const int i : isl.members)
         {
            cy.bodies[i].vS = cy.bodies[i].v;
            cy.bodies[i].hS = h;
         }
         if (hCap > 0.0)
         {
            if (hCap < h)
               why = TestDriver::ContactCost::kRetake;
            h = std::min(h, hCap);
            hCap = 0.0;
         }
         sn.h = h;
         if (s_cc != nullptr)
         {
            s_islWhy[why]++;
            s_islMinStep = std::min(s_islMinStep, h);
         }
         {
            const TestDriver::KindTimer kt(s_cc, TestDriver::ContactCost::kSurfaces);
            for (const int i : isl.members)
               SurfaceForces(cy, cy.bodies[i], h, s);
         }
         {
            const TestDriver::KindTimer kt(s_cc, TestDriver::ContactCost::kPairs);
            PairForces(cy, isl, h, s);
         }
         tripped = Integrate(cy, isl, h, s);
         {
            CPROF(PartStep);
            const TestDriver::KindTimer kt(s_cc, TestDriver::ContactCost::kMovingParts);
            for (ContactMovingPart *const hg : cy.islandParts)   // the part moves with the reactions of all its balls summed, in its own substeps within h
               hg->StepInIsland(h);
         }
         if (tripped && rewinds < kMaxRewinds && h > 4.0 * kMinStepS)
         {
            // the model left its range in this substep: take it again in two
            rewinds++;
            log("trip", "substep tripped (" + std::string(cy.bodies[isl.members[0]].tripWhy) + ")", cy.bodies[isl.members[0]].deepest, 0.0, 0.0, h, 0.5 * h);
            RwRestore(cy, isl, sn);
            t = sn.t;
            s = sn.s;
            hCap = 0.5 * h;
            tripped = false;
            havePrev = false;
            continue;
         }
         t += h;
         s++;
      }
      for (const int i : isl.members)
         cy.bodies[i].probesReady = false;
      // the cost counters were taken from the plan (the old plan: uniform steps); here is the real count
      g_last.substeps += s - cy.steps;
      for (const int i : isl.members)
         if (cy.bodies[i].e != nullptr)
            g_last.probes += (long long)cy.bodies[i].e->items.size() * (s - cy.steps);
      g_last.maxSteps = std::max(g_last.maxSteps, s);
      g_last.rewinds += rewinds;
      cy.steps = s;
      if (ContactTraceDetail())
         for (size_t k = 0; k < isl.edges.size(); k++)
         {
            const Body &A = cy.bodies[isl.edges[k].first], &B = cy.bodies[isl.edges[k].second];
            PLOGI << "ContactPair t=" << g_pplayer->m_time_msec << " balls " << isl.edges[k].first << ',' << isl.edges[k].second << " steps " << cy.steps << " max_delta_um " << cy.pst[k]->maxDelta * 1e6 << " max_fn " << cy.pst[k]->maxFn << " gap_end_um "
                  << ((A.c - B.c).Length() - A.bs.R - B.bs.R) * 1e6 << (tripped ? " TRIPPED" : "");
         }
      return tripped;
   }

   static bool ArrivalOn(const Cycle &cy)
   {
      static const int env = []() { const char *const e = getenv("VPX_ARRIVAL"); return e != nullptr && *e ? atoi(e) : -1; }();
      return env >= 0 ? env != 0 : cy.cs->arrivalStep;
   }

   // The island's substeps: at every one, all forces are summed from the current states before any body moves
   static bool Substeps(Cycle &cy, const ContactSystem::Islands::Island &isl)
   {
      if (ArrivalOn(cy) && IslandRewindable(cy, isl))
         return SubstepsArrival(cy, isl);
      CPROF(Substeps);
      using CC = TestDriver::ContactCost;
      const double h = cy.dt / cy.steps;
      if (s_cc != nullptr)
      {
         s_islWhy[cy.stepsWhy] += cy.steps;
         s_islMinStep = std::min(s_islMinStep, h);
      }
      bool tripped = false;
      int s = 0;
      for (; s < cy.steps && !tripped; s++)
      {
         {
            const TestDriver::KindTimer kt(s_cc, CC::kSurfaces);
            for (const int i : isl.members)
               SurfaceForces(cy, cy.bodies[i], h, s);
         }
         {
            const TestDriver::KindTimer kt(s_cc, CC::kPairs);
            PairForces(cy, isl, h, s);
         }
         tripped = Integrate(cy, isl, h, s);
         {
            const TestDriver::KindTimer kt(s_cc, CC::kMovingParts);
            for (ContactMovingPart *const hg : cy.islandParts)   // the part moves with the reactions of all its balls summed, in its own substeps within h
               hg->StepInIsland(h);
         }
      }
      // The parts keep stepping whatever happens to the ball: the rest of the cycle without the balls' reactions
      for (; s < cy.steps; s++)
      {
         const TestDriver::KindTimer kt(s_cc, CC::kMovingParts);
         for (ContactMovingPart *const hg : cy.islandParts)
            hg->StepInIsland(h);
      }
      if (ContactTraceDetail())
         for (size_t k = 0; k < isl.edges.size(); k++)
         {
            const Body &A = cy.bodies[isl.edges[k].first], &B = cy.bodies[isl.edges[k].second];
            PLOGI << "ContactPair t=" << g_pplayer->m_time_msec << " balls " << isl.edges[k].first << ',' << isl.edges[k].second << " steps " << cy.steps << " max_delta_um " << cy.pst[k]->maxDelta * 1e6 << " max_fn " << cy.pst[k]->maxFn << " gap_end_um "
                  << ((A.c - B.c).Length() - A.bs.R - B.bs.R) * 1e6 << (tripped ? " TRIPPED" : "");
         }
      return tripped;
   }

   // A step that left the model's range and could not be taken again. Nothing was written, so the island stays at the start of the cycle, and the same state would trip again. A ball that tripped
   // takes an inelastic stop there, which only takes energy out: pressed too deep, it loses its speed into the surfaces that press it (it goes on along them); after a runaway or a non-finite
   // step the state cannot be trusted and it loses all its motion. The contacts are let go and start again from zero (nothing stored): the near scan takes the balls back next cycle, on the
   // same law. Never stock. A run that gets here fails (TestMatrixDesign limit 5).
   static void StepFailed(Cycle &cy, const ContactSystem::Islands::Island &isl)
   {
      for (const int i : isl.members)
      {
         Body &w = cy.bodies[i];
         if (w.tripped)
         {
            const std::string name = w.o != nullptr ? w.o->m_name : w.band != nullptr ? w.band->m_name : std::string("a ball pair");
            static uint32_t lastMs = 0;
            if (lastMs == 0 || g_pplayer->m_time_msec - lastMs >= 500)
            {
               lastMs = g_pplayer->m_time_msec;
               PLOGW << "Contact step on " << name << " left the model's range (" << w.tripWhy << " at substep " << w.tripStep << " of " << cy.steps << ", depth " << w.deepest * 1e3 << " mm, speed " << w.v.Length() << " m/s, cluster " << isl.members.size() << "): the ball was stopped at (" << w.posStart.x << ", " << w.posStart.y << ", " << w.posStart.z << ")";
            }
            static std::vector<std::string> shown;   // one notification per surface per table start
            if (g_pplayer->m_liveUI != nullptr && std::find(shown.begin(), shown.end(), name) == shown.end())
            {
               shown.push_back(name);
               g_pplayer->m_liveUI->PushNotification("A contact step on " + name + " left the model's range: the ball was stopped (see vpinball.log)", 8000);
            }
            {
               V3 v((double)w.pb->m_d.m_vel.x, (double)w.pb->m_d.m_vel.y, (double)w.pb->m_d.m_vel.z);
               bool pressed = false;
               if (std::strcmp(w.tripWhy, "depth") == 0 && w.e != nullptr && w.probesReady)
                  for (const TubeProbe &pr : w.probes)
                     if (pr.overlap)
                     {
                        pressed = true;
                        const double vn = v.Dot(pr.normal);
                        if (vn < 0.0)
                           v = v - pr.normal * vn;
                     }
               if (pressed)
                  w.pb->m_d.m_vel = Vertex3Ds((float)v.x, (float)v.y, (float)v.z);
               else
               {
                  w.pb->m_d.m_vel = Vertex3Ds(0.f, 0.f, 0.f);
                  w.pb->m_angularmomentum = Vertex3Ds(0.f, 0.f, 0.f);
               }
            }
         }
         if (w.e != nullptr)
            cy.releases.emplace_back(w.pb, "tripped");
      }
   }

   // State back into the balls, the engagement bookkeeping and the release decision
   static void WriteBack(Cycle &cy, const ContactSystem::Islands::Island &isl)
   {
      CPROF(WriteBack);
      const double dt = cy.dt;
      for (const int i : isl.members)
      {
         Body &w = cy.bodies[i];
         HitBall *const pb = w.pb;
         BallContact *const ep = w.e;
         if (ep != nullptr)
         {
            if (w.omega2 > 0.0)   // a cycle without a contact (a ball hopping a micrometre above the surface) keeps the last known stiffness: the surface's stiffness has not changed, and without it the ball would be stepped at the impact step
               ep->omega = std::sqrt(w.omega2);
            ep->lastSpeed = w.v.Length();
            ep->lastVn = w.vnLast > 0.0 ? w.vnLast : 1e9;
            ep->steps += cy.steps;
            ep->dwellS += dt;
         }
         if (ContactTraceDetail() && ep != nullptr)
         {
            std::ostringstream o;
            for (SurfaceContact &it : ep->items)
            {
               const TubeProbe pr = it.r->Probe(w.bs, it, w.c);
               o << ' ' << it.r->Name() << ':' << (pr.overlap ? pr.delta * 1e6 : 0.0) << "um";
            }
            PLOGI << "ContactTrace t=" << g_pplayer->m_time_msec << " pos " << w.c.x / U << ' ' << w.c.y / U << ' ' << w.c.z / U << " v " << w.v.x << ' ' << w.v.y << ' ' << w.v.z << " w " << w.w.x << ' ' << w.w.y << ' ' << w.w.z
                  << " depth" << o.str() << " cluster " << isl.members.size();
         }
         pb->m_d.m_pos = Vertex3Ds((float)(w.c.x / U), (float)(w.c.y / U), (float)(w.c.z / U));
         pb->m_d.m_vel = Vertex3Ds((float)((w.v.x + w.extOwe.x) / VS), (float)((w.v.y + w.extOwe.y) / VS), (float)((w.v.z + w.extOwe.z) / VS));
         w.extOwe = V3();
         pb->m_angularmomentum = Vertex3Ds((float)(w.w.x * w.lScale), (float)(w.w.y * w.lScale), (float)(w.w.z * w.lScale));
         const Vertex3Ds ds = pb->m_d.m_pos - w.posStart;
         pb->m_lastEventSqrDist += ds.LengthSquared();
#ifdef C_DYNAMIC
         pb->m_drsq = ds.LengthSquared();
#endif
         pb->m_rubberAdvanced = true;   // the ball has been moved for the whole cycle
         // release: out of contact with every rubber and either moving away or more than 2 mm clear (the end of a hit), or the engagement has lasted max_ms
         // A ball within release_gap of a rubber stays engaged (hysteresis: the engine's contact pass would take it straight back).
         if (ep == nullptr)
            continue;   // a free ball has nothing to release
         BallContact &e = *ep;
         if (e.band != nullptr)   // the band holds the ball until its own release rule lets go (the ball stays in its other contacts)
         {
            // off the cord and clear of it; a plate on its way out can still catch a ball close by, so it stays in the group until the plate is quiet or the ball is far
            CordBand &cb = *e.band->Cord();
            int slot = -1;
            for (int k = 0; k < (int)cb.slots.size(); k++)
               if (cb.slots[k] == pb)
                  slot = k;
            const double gap = slot >= 0 ? cb.Gap(slot) : 1.0;
            if (gap > 0.03 || (gap > 0.001 && cb.PlateQuiet()))
               cy.bandReleases.emplace_back(pb, "clear");
            continue;
         }
         bool free = true, touching = false, closeBy = false;
         CPROF(ReleaseTest);
         for (const SurfaceContact &it : e.items)
         {
            const TubeProbe pr = it.r->Probe(w.bs, it, w.c);
            const double bound = PreventOn() ? it.r->SurfaceSpeedBoundMS(dt) * dt : 0.0;   // what the surface itself can cover in a cycle: a ball it can reach stays engaged
            const double releaseGap = it.r->ReleaseGapM() + bound;
            if (pr.overlap && pr.delta > 0.0)
            {
               free = false;
               touching = true;
            }
            else if (pr.dist > 0.0)
            {
               const double gap = pr.dist - (w.bs.R + pr.rho);
               if (gap <= releaseGap)
                  closeBy = true;
               if (gap <= releaseGap || !(w.v.Dot(pr.normal) > 0.0 || gap > 0.002 + bound))
                  free = false;
            }
         }
         if (free)
            cy.releases.emplace_back(pb, "clear");
         else if (w.o != nullptr && e.dwellS >= w.o->m_set.maxMs * 1e-3 && !touching && !closeBy)   // a ball in contact stays in the model however long (a resting ball, a pinch): re-engaging it would seat it and kick it
            cy.releases.emplace_back(pb, "timeout");
      }
   }

   // ---- The whole-ball energy check (ModelStandards 2b, rule 8): always on ----
   // A ball's energy: motion, spin, height, and what each of its contacts holds (normal, layer, shear stretch; a moving part's or a cord's own books are theirs). Every loss only lowers it, and height is in it, so:
   //  - between one cycle's end and the next one's start nothing may raise what the contacts hold: nothing but a step does work on them (a seat, a teleport, a face change, a law change);
   //  - across a cycle the island's energy may not rise, unless a drive is in it (a moving part, a band or cord, a magnet, the nudge).
   // A rise over 1 % of the energy present or 1 mJ (the test matrix's allowance) is a WARN naming ball, surface and depth, and on screen once per cause per table start.
   static void EnergyRaise(const char *const kind, const HitBall *const pb, const std::string &surface, const double depth0, const double depth1, const double rise)
   {
      static const void *owner = nullptr;
      static std::map<std::string, int> seen;
      if (owner != g_pplayer)
      {
         owner = g_pplayer;
         seen.clear();
      }
      int &count = seen[std::string(kind) + surface];
      if (count++ >= 20)
         return;
      const char *const name = TestDriver::BallName(pb);
      const Vertex3Ds &p = pb->m_d.m_pos;
      std::ostringstream o;
      o << "Energy check: ball " << (name != nullptr ? name : "") << " at (" << p.x << ", " << p.y << ", " << p.z << ") gained " << rise * 1e3 << " mJ " << kind << ", " << surface << " depth " << depth0 * 1e3 << " -> " << depth1 * 1e3 << " mm, t=" << g_pplayer->m_time_msec;
      PLOGW << o.str();
      if (count == 1 && g_pplayer->m_liveUI != nullptr)
         g_pplayer->m_liveUI->PushNotification(o.str(), 6000);
   }

   // what one contact holds, by the law its loop uses; depth: the depth it counts
   static double ContactHolds(const SurfaceContact &it, const BallSI &bs, const V3 &c, const bool group, double &depth)
   {
      depth = 0.0;
      if (it.r->MovingPart() != nullptr || it.NoFaces(c))
         return 0.0;
      const TubeProbe pr = group ? it.r->ProbePatch(bs, it, c) : it.r->Probe(bs, it, c);
      if (!pr.overlap)
         return 0.0;
      depth = group ? pr.delta - it.gs.zero : pr.delta;
      if (depth <= 0.0)
         return 0.0;
      const double reff = group && it.gs.reff > 0.0 ? it.gs.reff : pr.reff > 0.0 ? pr.reff : bs.R;
      const RubberParams &rp = it.rp;
      const double E = rp.eMpa * 1e6, Es = E / (1.0 - rp.nu * rp.nu), Gs = E / (2.0 * (1.0 + rp.nu)) / (2.0 - rp.nu);
      const double kH = 4.0 / 3.0 * Es * std::sqrt(reff), cL = kH * rp.betaLayer * reff * std::sqrt(reff);
      double held;
      if (group && cL > 0.0 && it.hLayer > 0.0)   // the group's law: the layer's thickness left
      {
         GroupStep::Item law;
         law.kH = kH;
         law.cL = cL;
         law.hL = it.hLayer;
         held = law.V(std::min(depth, 0.999 * it.hLayer));
      }
      else   // the old loop's: on the full thickness
         held = 0.4 * kH * depth * depth * std::sqrt(depth) + (it.hLayer > 0.0 ? cL * depth * depth * depth * depth / (4.0 * it.hLayer * it.hLayer * it.hLayer) : 0.0);
      const V3 &xi = group ? it.gs.xi : it.st.xi;
      return held + 4.0 * Gs * std::sqrt(reff * depth) * rp.shearScale * xi.Dot(xi);
   }

   // what a ball pair holds, by the law of both loops: Hertz on the pair's modulus (cy.rpBall holds it halved), and the shear stretch; its depth counts from its zero
   static double PairHolds(const Cycle &cy, const Body &A, const Body &B, const BallPairState &ps)
   {
      const double depth = A.bs.R + B.bs.R - (A.c - B.c).Length() - ps.zero;
      if (depth <= 0.0)
         return 0.0;
      const double reff = A.bs.R * B.bs.R / (A.bs.R + B.bs.R);
      const RubberParams &rp = cy.rpBall;
      const double E = rp.eMpa * 1e6, Es = E / (1.0 - rp.nu * rp.nu), Gs = E / (2.0 * (1.0 + rp.nu)) / (2.0 - rp.nu);
      return 0.4 * (4.0 / 3.0 * Es * std::sqrt(reff)) * depth * depth * std::sqrt(depth) + 4.0 * Gs * std::sqrt(reff * depth) * rp.shearScale * ps.st.xi.Dot(ps.st.xi);
   }

   // Between cycles nothing does work on a pair (no step has run), so a pair may not hold more than it did: one found inside that was not inside at the last cycle's end (a ball created, placed or
   // moved onto another) has stored nothing. Its depth is its zero; the zero only ever comes down as the balls back out, so the pair can give back only what a step stored. No position moves.
   static void PairFirstSight(Cycle &cy)
   {
      for (const auto &p : cy.pairs)
      {
         const Body &A = cy.bodies[p.first], &B = cy.bodies[p.second];
         const double depth = A.bs.R + B.bs.R - (A.c - B.c).Length();
         BallPairState *const ps = FindPair(A.pb, B.pb);
         const bool placed = std::find(g_placed.begin(), g_placed.end(), A.pb) != g_placed.end() || std::find(g_placed.begin(), g_placed.end(), B.pb) != g_placed.end();
         if (ps == nullptr)
         {
            if (depth > 0.0)
            {
               MakePair(A.pb, B.pb).zero = depth;
               if (!placed && depth > 1e-9)
                  DepthFault(DepthCause::FirstSight, A.pb, "a ball pair", depth);
            }
         }
         else if (ps->endDepth <= ps->zero && depth > ps->zero)
         {
            if (!placed && depth - ps->zero > 1e-9)
               DepthFault(DepthCause::StockMoved, A.pb, "a ball pair", depth - ps->zero);
            ps->zero = depth;
         }
      }
   }

   static double BodyMotion(const Body &w, const V3 &g)
   {
      return 0.5 * w.bs.m * w.v.Dot(w.v) + 0.5 * w.bs.Inertia() * w.w.Dot(w.w) - w.bs.m * g.Dot(w.c);
   }

   static V3 CheckGravity()   // the table's gravity, m/s2, without the nudge
   {
      const Vertex3Ds &gv = g_pplayer->m_physics->GetGravity();
      return V3(gv.x, gv.y, gv.z) * (VS / DEFAULT_STEPTIME_S);
   }

   // The cycle's start, every held ball: what its contacts hold against what they held at the last cycle's end
   static void EnergyCheckStart(Cycle &cy)
   {
      const V3 g = CheckGravity();
      for (size_t i = 0; i < cy.bodies.size(); i++)   // VPX_TRACE_BALL (off by default)
      {
         int island = -1;
         for (size_t k = 0; k < cy.islands.count && island < 0; k++)
            if (std::find(cy.islands.list[k].members.begin(), cy.islands.list[k].members.end(), (int)i) != cy.islands.list[k].members.end())
               island = (int)k;
         ContactTraceBall("EnergyCheckStart", cy.bodies[i].pb, true, island);
      }
      for (Body &w : cy.bodies)
      {
         double held = 0.0, worst = 0.0, depth0 = 0.0, depth1 = 0.0;
         const SurfaceContact *where = nullptr;
         if (w.e != nullptr)
            for (const SurfaceContact &it : w.e->items)
            {
               double d = 0.0;
               const double v = w.e->seatPending ? 0.0 : ContactHolds(it, w.bs, w.c, w.e->groupLoop, d);   // (found inside and not seated: it counts from zero)
               held += v;
               const double was = w.e->endStored >= 0.0 ? it.endDepth : 0.0;
               if (d - was > worst)
               {
                  worst = d - was;
                  where = &it;
                  depth0 = was;
                  depth1 = d;
               }
            }
         w.eStart = BodyMotion(w, g) + held;
         w.eScale = 0.5 * w.bs.m * w.v.Dot(w.v) + 0.5 * w.bs.Inertia() * w.w.Dot(w.w) + held;
         if (w.e != nullptr)
         {
            const double rise = held - std::max(w.e->endStored, 0.0);
            if (rise > std::max(1e-3, 0.01 * w.eScale))
               EnergyRaise(w.e->endStored >= 0.0 ? "between cycles" : "at first sight", w.pb, where != nullptr ? where->r->Name() : std::string("?"), depth0, depth1, rise);
         }
      }
      for (const auto &p : cy.pairs)   // what the ball pairs hold, half to each ball
      {
         Body &A = cy.bodies[p.first], &B = cy.bodies[p.second];
         const BallPairState *const ps = FindPair(A.pb, B.pb);
         if (ps == nullptr)
            continue;
         const double held = PairHolds(cy, A, B, *ps);
         A.eStart += 0.5 * held;
         B.eStart += 0.5 * held;
         A.eScale += 0.5 * held;
         B.eScale += 0.5 * held;
         const double rise = held - ps->endHeld;
         if (rise > std::max(1e-3, 0.01 * (A.eScale + B.eScale)))
         {
            const double depth = A.bs.R + B.bs.R - (A.c - B.c).Length() - ps->zero;
            EnergyRaise("between cycles", A.pb, "ball pair", std::max(ps->endDepth - ps->zero, 0.0), depth, rise);
         }
      }
   }

   // The island's end: the rise across the cycle, and what the contacts hold now for the next cycle's start
   static void EnergyCheckEnd(Cycle &cy, const ContactSystem::Islands::Island &isl)
   {
      const V3 g = CheckGravity();
      bool driven = false;
      double e0 = 0.0, e1 = 0.0, scale = 0.0, top = -1e300;
      const Body *rose = nullptr;
      for (const int i : isl.members)
         cy.bodies[i].endPair = 0.0;
      for (const auto &p : isl.edges)   // what the ball pairs hold now, half to each ball; and for the next cycle's start
      {
         Body &A = cy.bodies[p.first], &B = cy.bodies[p.second];
         BallPairState *const ps = FindPair(A.pb, B.pb);
         if (ps == nullptr)
            continue;
         ps->endDepth = A.bs.R + B.bs.R - (A.c - B.c).Length();
         ps->endHeld = PairHolds(cy, A, B, *ps);
         A.endPair += 0.5 * ps->endHeld;
         B.endPair += 0.5 * ps->endHeld;
      }
      for (const int i : isl.members)
      {
         Body &w = cy.bodies[i];
         driven = driven || w.band != nullptr || w.magAcc > 0.0
            || (w.extBaseVS.Length() > 0.0 && (w.extBaseVS * (1.0 / PHYSICS_STEPTIME_S) - g).Length() > 0.02 * g.Length());   // (the nudge)
         double held = 0.0;
         if (w.e != nullptr)
         {
            for (SurfaceContact &it : w.e->items)
            {
               driven = driven || (!it.partOut && it.r->MovingPart() != nullptr);
               held += ContactHolds(it, w.bs, w.c, w.e->groupLoop, it.endDepth);
            }
            w.e->endStored = held;
         }
         const double e = BodyMotion(w, g) + held + w.endPair;
         e0 += w.eStart;
         e1 += e;
         scale += w.eScale;
         if (e - w.eStart > top)
         {
            top = e - w.eStart;
            rose = &w;
         }
      }
      if (driven || rose == nullptr || e1 - e0 <= std::max(1e-3, 0.01 * scale))
         return;
      // the surface it went deepest into, the floor last: the likely cause, for the message
      const SurfaceContact *where = nullptr;
      if (rose->e != nullptr)
         for (const SurfaceContact &it : rose->e->items)
            if (it.endDepth > 0.0 && (where == nullptr || where->r->IsFloor() || (!it.r->IsFloor() && it.endDepth > where->endDepth)))
               where = &it;
      EnergyRaise("across a cycle", rose->pb, where != nullptr ? where->r->Name() : std::string("no surface"), 0.0, where != nullptr ? where->endDepth : 0.0, e1 - e0);
   }

   // The bats that took part: finish the cycle (the pivot's nudge, then UpdateDisplacements only publishes the angle and fires the events)
   // The island's members for the test driver's cost split: its balls (the driver's names), then the surfaces, flippers and bands they touch at the end of the step
   // ---- An island with a cord: one group step for the whole island (docs/GroupStepDesign.txt 8, 11) ----
   // A surface of a ball's contact as a shape of the group (its own list of near faces, so one ball's)
   class SurfaceShape : public GroupStep::Shape
   {
   public:
      const ContactSurface *r = nullptr;
      const SurfaceContact *it = nullptr;
      BallSI bs;
      double heldReff = 0.0;   // > 0: the radius this touch started with (the group takes a shape's law anew each cycle; a polygon's radius jumps at a corner, the energy must not)
      double heldKb = 0.0;     // the bend this touch started with, for the same reason
      ContactMovingPart *mp = nullptr;   // a part on the group: its own surface, which moves with it
      int surface = 0;
      bool NearPieces(const V3 &c, const double R, const double reach, const double swing, const double minReach, std::vector<int> &list) const override
      {
         return mp != nullptr && mp->GroupNearPieces(surface, c, R, reach, swing, minReach, list);
      }
      bool NearHolds(const V3 &c, const double R, const double dq, const std::vector<int> &list) const override { return mp == nullptr || mp->GroupNearHolds(surface, c, R, dq, list); }
      void Probe(const V3 &c, double R, GroupStep::Hit &h) const override
      {
         if (mp != nullptr)
         {
            mp->GroupProbe(surface, c, R, trial, h);
            if (heldReff > 0.0)
            {
               h.reff = heldReff;
               h.kb = heldKb;
            }
            return;
         }
         const TubeProbe pr = r->ProbePatch(bs, *it, c);
         h.delta = pr.overlap ? pr.delta : (pr.dist > 0.0 ? R + pr.rho - pr.dist : -1.0);
         h.n = pr.normal;
         h.grad = pr.grad;
         h.reff = heldReff > 0.0 ? heldReff : pr.reff > 0.0 ? pr.reff : R;
         h.curvDist = pr.reff > 0.0 && pr.reff < 0.999 * R ? pr.dist : 0.0;
      }
   };

   // a ball pair's state between cycles (a < b by address)
   static BallPairState *FindPair(const HitBall *a, const HitBall *b)
   {
      if (b < a)
         std::swap(a, b);
      for (BallPairState &p : g_ballPairs)
         if (p.a == a && p.b == b)
            return &p;
      return nullptr;
   }
   static BallPairState &MakePair(const HitBall *const a, const HitBall *const b)
   {
      if (BallPairState *const p = FindPair(a, b))
         return *p;
      g_ballPairs.push_back(BallPairState());
      BallPairState &p = g_ballPairs.back();
      p.a = std::min(a, b);
      p.b = std::max(a, b);
      return p;
   }

   static bool IslandHasCord(const Cycle &cy, const ContactSystem::Islands::Island &isl)
   {
      for (const int i : isl.members)
         if (cy.bodies[i].band != nullptr)
            return true;
      return false;
   }

   // An island holding a part on the group step (a bat): always the group, whatever else is in it. Such a part has no other path: no old loop, no stock.
   static bool IslandHasGroupPart(const Cycle &cy, const ContactSystem::Islands::Island &isl)
   {
      for (const int i : isl.members)
         if (cy.bodies[i].e != nullptr)
            for (const SurfaceContact &it : cy.bodies[i].e->items)
               if (ContactMovingPart *const mp = it.r->MovingPart(); mp != nullptr && !it.partOut && mp->GroupPart() != nullptr)
                  return true;
      return false;
   }

   // An island goes to the group step without a cord ([Rubbers] ring_model, GroupStepDesign 14.1) when a ball can touch a rubber surface within this cycle, or still counts a squash from a
   // zero only the group knows. An island with a moving part stays on the old loop until that part moves over (14.2). The reach is the proximity engagement's, so a ball joins before it touches.
   static bool IslandOnRing(const Cycle &cy, const ContactSystem::Islands::Island &isl)
   {
      bool want = false;
      for (const int i : isl.members)
      {
         const Body &w = cy.bodies[i];
         if (w.band != nullptr)
            return false;
         if (w.e == nullptr)
            continue;
         for (const SurfaceContact &it : w.e->items)
         {
            if (it.partOut)   // (island_reach: a part out of this ball's reach is not in its solve; the ball keeps the loop it has without it)
               continue;
            if (it.r->MovingPart() != nullptr)
               return false;
            if (want)
               continue;
            if (w.e->groupLoop && it.gs.zero > 0.0)
               want = true;
            else if (it.r->GroupModel() && !it.NoFaces(w.c))
            {
               const TubeProbe pr = it.r->Probe(w.bs, it, w.c);
               want = pr.overlap || (pr.dist > 0.0 && pr.dist - (w.bs.R + pr.rho) < 2.0 * w.v.Length() * cy.dt + it.r->ReleaseGapM());
            }
         }
      }
      return want;
   }

   static bool g_contactKeptEnd(const HitBall *const pb)   // the ball's contact has a last cycle's end to compare with
   {
      return pb->m_contact != nullptr && pb->m_contact->endStored >= 0.0;
   }

   // Why a contact of the group went deeper with no step. Placed by a script or the driver, or handed over from the old loop (to go with that loop), is expected. Anything else is a fault: its zero has
   // taken the depth, so nothing is stored and the ball goes on, and it is shown once per cause and surface (rule 8).
   enum class DepthCause { None, Placed, HandedOver, StockMoved, Faces, FirstSight };
   static void DepthFault(const DepthCause why, const HitBall *const pb, const std::string &surface, const double depth)
   {
      static const void *owner = nullptr;
      static std::map<std::string, int> seen;
      if (owner != g_pplayer)
      {
         owner = g_pplayer;
         seen.clear();
      }
      const char *const text = why == DepthCause::StockMoved ? "the ball was moved outside any step (a stock part)" : why == DepthCause::FirstSight ? "the ball was already inside at its first touch" : "the surface's faces in reach changed";
      int &count = seen[std::string(text) + surface];
      if (count++ >= 20)
         return;
      const char *const name = TestDriver::BallName(pb);
      const Vertex3Ds &p = pb->m_d.m_pos;
      std::ostringstream o;
      o << "Contact: ball " << (name != nullptr ? name : "") << " at (" << p.x << ", " << p.y << ", " << p.z << ") is " << depth * 1e6 << " um deeper in " << surface << " with no step: " << text << ". It counts from there, nothing stored. t=" << g_pplayer->m_time_msec;
      PLOGW << o.str();
      if (count == 1 && g_pplayer->m_liveUI != nullptr)
         g_pplayer->m_liveUI->PushNotification(o.str(), 6000);
   }

   // A ball coming from the old loop: the same touch, the same stretch, so the same energy; a surface it is found inside has no history and counts from zero (14.4). Goes with the old loop.
   static void RingHandIn(BallContact &e, const BallSI &bs, const V3 &c)
   {
      if (e.groupLoop)
         return;
      for (SurfaceContact &it : e.items)
      {
         it.gs.touching = it.touching && e.dwellS > 0.0;
         it.gs.xi = it.st.xi;
         it.gs.ra = it.st.ra;
         it.gs.zero = it.gs.reff = 0.0;
         if (!it.gs.touching || it.hLayer <= 0.0 || it.r->MovingPart() != nullptr || it.NoFaces(c))
            continue;
         // the group's law at this depth may not hold more than the old loop's did: its zero takes up the rest (a change of loop does no work)
         const TubeProbe pr = it.r->Probe(bs, it, c);
         if (!pr.overlap)
            continue;
         const double reff = pr.reff > 0.0 ? pr.reff : bs.R;
         const RubberParams &rp = it.rp;
         const double Es = rp.eMpa * 1e6 / (1.0 - rp.nu * rp.nu), kH = 4.0 / 3.0 * Es * std::sqrt(reff), cL = kH * rp.betaLayer * reff * std::sqrt(reff);
         if (cL <= 0.0)
            continue;
         const double d = pr.delta, hL = it.hLayer;
         const double heldOld = 0.4 * kH * d * d * std::sqrt(d) + cL * d * d * d * d / (4.0 * hL * hL * hL);
         GroupStep::Item law;
         law.kH = kH;
         law.cL = cL;
         law.hL = hL;
         if (d < 0.999 * hL && law.V(d) <= heldOld)
            continue;
         double lo = 0.0, hi = std::min(d, 0.999 * hL);
         for (int r = 0; r < 40; r++)
         {
            const double mid = 0.5 * (lo + hi);
            (law.V(mid) < heldOld ? lo : hi) = mid;
         }
         it.gs.zero = d - lo;
      }
      e.groupLoop = true;
   }

   // The island stays on (or falls back to) the old loop. The stretch goes back to its copy there. The old loop has no zero: a ball still inside a surface is seated, as it would have been at
   // engagement; that moves it out and stores nothing.
   static void RingHandOut(Cycle &cy, const ContactSystem::Islands::Island &isl)
   {
      for (const int i : isl.members)
      {
         Body &w = cy.bodies[i];
         if (w.e == nullptr || !(w.e->groupLoop || w.e->seatPending))
            continue;
         BallContact &e = *w.e;
         bool seat = e.seatPending;
         for (SurfaceContact &it : e.items)
         {
            if (e.groupLoop)
            {
               it.st.xi = it.gs.xi;
               it.st.ra = it.gs.ra;
            }
            seat = seat || it.gs.zero > 0.0;
            it.gs.touching = false;
            it.gs.zero = it.gs.reff = 0.0;
         }
         e.groupLoop = e.seatPending = false;
         if (seat)
         {
            V3 c = w.c;
            double delta0 = 0.0, rewound = 0.0;
            SeatAtTouch(e.items, w.bs, c, w.v, delta0, rewound);
            w.c = c;
         }
      }
   }

   // the group's balls change between cycles: each cord's own per-ball state follows them
   // What a ball's state was to double precision when it left a group step, with the floats that state was written as (the rest of the engine keeps floats). While its floats are still those, the doubles go on: a rounded place
   // would make the cord re-solve its path, and a pulley could change with nothing paying. Keyed by the ball; a ball gone leaves an entry no other ball's floats match.
   struct KeptDoubles { Vertex3Ds pos, vel, am; V3 c, v, w; };
   inline static std::unordered_map<const HitBall *, KeptDoubles> g_keptDoubles;

   static void CordFollow(CordBand &cb, const GroupStep::Ball *const groupBalls, HitBall *const *balls, const int n)
   {
      if (!cb.attached)
      {
         cb.fixed->AttachGroup(groupBalls, n);
         cb.attached = true;
      }
      else
      {
         bool same = (int)cb.slots.size() == n;
         for (int i = 0; same && i < n; i++)
            same = cb.slots[i] == balls[i];
         if (!same)
         {
            std::vector<int> oldOfNew(n, -1);
            for (int i = 0; i < n; i++)
               for (int k = 0; k < (int)cb.slots.size(); k++)
                  if (cb.slots[k] == balls[i])
                     oldOfNew[i] = k;
            cb.fixed->RemapBalls(oldOfNew.data(), n, groupBalls);
         }
      }
      if (n > 0)
         cb.fixed->BeginGroupCycle(groupBalls);
      cb.slots.assign(balls, balls + n);
   }

   // The group parts in this island's solve. A part at rest stays out while no ball of the island can reach it in this cycle: Group::OutOfReach's test over the whole cycle. A ball can pick up speed inside the cycle, so the
   // speed is the bound on what it can have: all the island's balls' (a ball's hit passes it on), what its forces can add, and twice the fastest surface in the island that is not at rest (a ball leaves a moving surface at
   // up to twice its speed). Out of the solve it moves alone (AdvanceAlone), so nothing about it stops; the script and a nudge can only act between cycles.
   static void PartsInSolve(Cycle &cy, const ContactSystem::Islands::Island &isl, std::vector<const ContactMovingPart *> &in)
   {
      in.clear();
      double sp = 0.0, acc = 0.0, vPart = 0.0;
      for (const int i : isl.members)
      {
         const Body &w = cy.bodies[i];
         sp += w.v.Length();
         acc = std::max(acc, w.extBaseVS.Length() / PHYSICS_STEPTIME_S + w.magAcc + w.F.Length() / w.bs.m);
         if (w.e != nullptr)
            for (const SurfaceContact &it : w.e->items)
               if (ContactMovingPart *const mp = it.r->MovingPart(); mp != nullptr && !it.partOut && !mp->AtRest())
               {
                  mp->SetBall(w.bs);
                  vPart = std::max(vPart, mp->SurfaceSpeedBound(cy.dt));
               }
      }
      const double reach = 2.0 * (sp + acc * cy.dt + 2.0 * vPart) * cy.dt + 1e-4;
      for (const int i : isl.members)
      {
         Body &w = cy.bodies[i];
         if (w.e != nullptr)
            for (SurfaceContact &it : w.e->items)
            {
               ContactMovingPart *const mp = it.r->MovingPart();
               if (mp == nullptr || it.partOut || mp->GroupPart() == nullptr || std::find(in.begin(), in.end(), mp) != in.end())
                  continue;
               bool joins = !mp->AtRest();
               const int surfaces = mp->GroupSurfaces();
               if ((int)it.gsMore.size() < surfaces - 1)
                  it.gsMore.resize(surfaces - 1);
               for (int sf = 0; sf < surfaces && !joins; sf++)
                  joins = it.Surf(sf).touching;
               if (!joins)
               {
                  mp->SetBall(w.bs);
                  joins = mp->Probe(w.c).delta > -reach;
               }
               if (joins)
                  in.push_back(mp);
            }
      }
   }

   // A cells magnet's pull on the group's balls: what the ball's place stores at the coil's current (magnetcore.h Model::Energy), the current held for the cycle. Behind the pole face's plane the tables have nothing: what the plane's points store goes on, unchanged in height.
   class PullField final : public GroupStep::BallField
   {
   public:
      Magnets::CellPull p;
      double Eval(const V3 &c, const double R, V3 *const grad, GroupStep::M3 *const hess) const override { return At(p.amps, c, R, grad, hess); }
      double At(const double amps, const V3 &c, const double R, V3 *const grad, GroupStep::M3 *const hess) const
      {
         const V3 n(p.n[0], p.n[1], p.n[2]), d = c - V3(p.o[0], p.o[1], p.o[2]);
         const double zc = d.Dot(n);
         const V3 rv = d - n * zc;
         const double x = rv.Length(), h = zc > 0.0 ? zc - R : -R;
         MagnetCore::Model::Pull e;
         p.core->Energy(x, h, amps, e);
         if (zc <= 0.0)
            e.Vh = e.Vxh = e.Vhh = 0.0;
         const V3 ur = x > 1e-12 ? rv * (1.0 / x) : V3();
         if (grad != nullptr)
            *grad = n * e.Vh + ur * e.Vx;
         if (hess != nullptr)
         {
            using M3 = GroupStep::M3;
            const M3 nn = M3::Outer(n, n), rr = M3::Outer(ur, ur), nr = M3::Outer(n, ur);
            const double across = x > 1e-9 ? e.Vx / x : e.Vxx;   // the turn of the radial direction round the axis
            for (int i = 0; i < 3; i++)
               for (int j = 0; j < 3; j++)
                  hess->a[i][j] = e.Vhh * nn.a[i][j] + e.Vxx * rr.a[i][j] + e.Vxh * (nr.a[i][j] + nr.a[j][i]) + across * ((i == j ? 1.0 : 0.0) - nn.a[i][j] - rr.a[i][j]);
         }
         return e.V;
      }
   };

   // ring: the island has no cord (IslandOnRing). Returns false when a step of such an island could not be solved: nothing has been written, and the cycle is retaken on the old loop, the contact
   // system's own (GroupStepDesign 14.3; never stock, no hold). A cord island keeps its own rule (15.3) and returns true. So does an island with a part on the group step (a bat), which has no old
   // loop to go to: what is left of the cycle is lost, no energy can appear.
   static bool CordIsland(Cycle &cy, const ContactSystem::Islands::Island &isl, const bool ring = false)
   {
      static constexpr int kPartEnergyLines = 20;
      static GroupStep::Group g;
      static std::deque<SurfaceShape> shapes;
      shapes.clear();
      g.nb = g.ns = g.np = g.nf = 0;
      g.nform = 0;
      static std::deque<PullField> pulls;
      pulls.clear();
      for (const Magnets::CellPull &cp : Magnets::CellPulls())
      {
         pulls.emplace_back();
         pulls.back().p = cp;
         g.AddField(&pulls.back());
      }
      g.set.pair = cy.rpBall;
      g.set.rollSmooth = cy.cs->rollSmoothRadS;
      g.set.rollModel = cy.cs->rollModel;
      g.set.stepModel = cy.cs->stepModel;
      g.set.swingCap = cy.cs->stepSwingCap;
      g.set.shearMotion = cy.cs->shearStepMotion;
      g.set.termArrival = cy.cs->termArrival;
      g.set.dampDirSlope = cy.cs->dampDirSlope;
      g.set.epsRest = cy.cs->epsRestUm * 1e-6;   // (0.01 um gives 1e-8 m exactly, the Settings default)
      if (ring)   // (cy.rpBall holds the old loop's halved modulus, and the group halves it for a pair itself; its rolling resistance is on the pair's radius, R / 2)
      {
         g.set.pair.eMpa = cy.cs->ballEMpa;
         g.set.pair.rollingMu = 0.5 * cy.cs->ballRollingResistance;
      }
      g.set.hMin = GroupStep::GroupFixedStepS() > 0.0 ? 1e-9 : 2e-6;   // a step that cannot be solved gives up here: in real time it must not spend the frame halving (a fixed-step reference may)
      const bool hasPart = IslandHasGroupPart(cy, isl);
      static std::vector<HitBall *> balls;
      static std::vector<DepthCause> cause;
      balls.resize(isl.members.size());
      cause.resize(isl.members.size());
      g.Size((int)isl.members.size());
      std::vector<HitRubber *> cords;
      static std::vector<ContactMovingPart *> parts;   // the moving parts on the group in this island
      parts.clear();
      static std::vector<const ContactMovingPart *> partsIn;
      PartsInSolve(cy, isl, partsIn);
      bool gSet = false;
      for (const int i : isl.members)
      {
         Body &w = cy.bodies[i];
         if (!ring && g.nb >= CordBand::kMaxBalls)   // (the cord's own limit)
         {
            PLOGW << "Cord group: more than " << CordBand::kMaxBalls << " balls in one island: the others are left out this cycle";
            break;
         }
         if (!w.extTaken)   // the group gives every ball its gravity in its steps, wherever the ball is (a ball high on the cord has no floor in its contact)
            Sources::Begin(w, cy.dt, true);
         const int ib = g.nb++;
         balls[ib] = w.pb;
         GroupStep::Ball &b = g.ball[ib];
         b.m = w.bs.m;
         b.R = w.bs.R;
         b.c = w.c;
         b.v = w.v;
         b.w = w.w;
         if (ring)
         {
            const auto kd = g_keptDoubles.find(w.pb);
            if (kd != g_keptDoubles.end() && kd->second.pos == w.pb->m_d.m_pos)
            {
               const KeptDoubles &k = kd->second;
               b.c = k.c;
               b.v = w.v + (k.v - V3((double)k.vel.x * VS, (double)k.vel.y * VS, (double)k.vel.z * VS));
               b.w = w.w + (k.w - V3((double)k.am.x / w.lScale, (double)k.am.y / w.lScale, (double)k.am.z / w.lScale));
            }
         }
         else
         {
            HitBall *const pb = w.pb;
            const auto kd = g_keptDoubles.find(pb);
            if (kd != g_keptDoubles.end() && kd->second.pos == pb->m_d.m_pos && kd->second.vel == pb->m_d.m_vel && kd->second.am == pb->m_angularmomentum
               && w.c.x == (double)pb->m_d.m_pos.x * U && w.c.y == (double)pb->m_d.m_pos.y * U && w.c.z == (double)pb->m_d.m_pos.z * U   // (nothing has moved it this cycle before here)
               && w.v.x == (double)pb->m_d.m_vel.x * VS && w.v.y == (double)pb->m_d.m_vel.y * VS && w.v.z == (double)pb->m_d.m_vel.z * VS)
            {
               b.c = kd->second.c;
               b.v = kd->second.v;
               b.w = kd->second.w;
            }
         }
         b.fExt = V3();
         if (!gSet && w.extBaseVS.Length() > 0.0)
         {
            g.set.g = w.extBaseVS * (1.0 / PHYSICS_STEPTIME_S);   // gravity along the table's slope, and the nudge
            gSet = true;
         }
         Sources::EddyBrake(w);   // (the old loop's own: its rate and direction, held for the cycle as the coil's current is)
         b.brakeRate = w.eddyRate;
         b.brakeAxis = w.eddyDir;
         // a cells magnet's pull: a stored energy of the ball's place; on the body it is the same acceleration as AddVelocityIncrement's (m / the steel mass it uses)
         b.fieldScale = g.nf > 0 && w.pb->m_magnetic ? b.m / (7800.0 * 4.0 / 3.0 * M_PI * b.R * b.R * b.R * std::max(0.001f, w.pb->m_d.m_mass)) : 0.0;
         if (w.magAcc > 0.0)   // what stays held for the cycle: the other balls' part, the magnets not on a field: a drive
         {
            const float mp[3] = { w.pb->m_d.m_pos.x, w.pb->m_d.m_pos.y, w.pb->m_d.m_pos.z };
            // the other balls' part: a stored energy of place when the solve gave a form (pairform.h, 2026-11-03, added below), else held as before (the dipoles)
            const bool form = Magnets::HasPairForm(w.pb);
            float mdv[3] = { form ? 0.f : w.pb->m_magPair[0], form ? 0.f : w.pb->m_magPair[1], form ? 0.f : w.pb->m_magPair[2] };
            Magnets::AddVelocityIncrement(mp, w.pb->m_d.m_radius, w.pb->m_d.m_mass, mdv, true);
            b.fExt = V3(mdv[0] * VS, mdv[1] * VS, mdv[2] * VS) * (b.m / PHYSICS_STEPTIME_S);
         }
         if (w.band != nullptr && std::find(cords.begin(), cords.end(), w.band) == cords.end())
            cords.push_back(w.band);
         if (ring)
         {
            const auto kd = g_keptDoubles.find(w.pb);
            const bool same = kd != g_keptDoubles.end() && kd->second.pos == w.pb->m_d.m_pos;
            cause[ib] = std::find(g_placed.begin(), g_placed.end(), w.pb) != g_placed.end() ? DepthCause::Placed
               : w.e != nullptr && !w.e->groupLoop && w.e->endStored >= 0.0 ? DepthCause::HandedOver
               : w.e != nullptr && w.e->endStored < 0.0 ? DepthCause::FirstSight
               : same ? DepthCause::Faces : DepthCause::StockMoved;
         }
         if (ring && w.e != nullptr)
            RingHandIn(*w.e, w.bs, w.c);
         if (w.e != nullptr)
            for (SurfaceContact &it : w.e->items)
            {
               if (ring && it.NoFaces(w.c))
               {
                  it.gs.touching = it.touching = false;
                  it.gs.xi = V3();
                  it.gs.ra = RubberContact::RollAngle();
                  it.gs.zero = it.gs.reff = 0.0;
                  continue;
               }
               if (it.partOut)   // (island_reach: a part out of this ball's reach is solved elsewhere)
                  continue;
               ContactMovingPart *const mp = it.r->MovingPart() != nullptr && it.r->MovingPart()->GroupPart() != nullptr ? it.r->MovingPart() : nullptr;
               if (mp != nullptr && std::find(partsIn.begin(), partsIn.end(), mp) == partsIn.end())
                  continue;
               const int surfaces = mp != nullptr ? mp->GroupSurfaces() : 1;
               if ((int)it.gsMore.size() < surfaces - 1)
                  it.gsMore.resize(surfaces - 1);
               for (int sf = 0; sf < surfaces; sf++)
               {
                  shapes.emplace_back();
                  SurfaceShape &sh = shapes.back();
                  sh.r = it.r;
                  sh.it = &it;
                  sh.bs = w.bs;
                  sh.rp = it.rp;
                  sh.layerM = it.hLayer;
                  sh.onlyBall = ib;
                  sh.heldReff = (ring || mp != nullptr) && it.Surf(sf).touching ? it.Surf(sf).reff : 0.0;
                  sh.heldKb = mp != nullptr && it.Surf(sf).touching ? it.Surf(sf).kb : 0.0;
                  if (hasPart && mp == nullptr && it.r->MovingPart() != nullptr)   // a part still on the old loop is a surface at its present pose in this cycle, its reaction left out, until it moves over (GroupStepDesign 3)
                  {
                     static const void *warned = nullptr;
                     if (warned != g_pplayer)
                     {
                        warned = g_pplayer;
                        const std::string msg = "Group step: a part on the group step and " + it.r->Name() + " share an island: " + it.r->Name() + " counts as a fixed surface in those cycles until it moves to the group step";
                        PLOGW << msg;
                        g_pplayer->m_liveUI->PushNotification(msg, 6000);
                     }
                  }
                  if (mp != nullptr)
                  {
                     sh.mp = mp;
                     sh.surface = sf;
                     if (const RubberParams *const law = mp->GroupLaw(sf, sh.layerM))
                        sh.rp = *law;
                     const auto at = std::find(parts.begin(), parts.end(), mp);
                     sh.part = (int)(at - parts.begin());
                     if (at == parts.end())
                     {
                        mp->CycleBegin();
                        parts.push_back(mp);
                        g.AddPart(mp->GroupPart());
                     }
                  }
                  g.AddShape(&sh);
               }
            }
      }
      // the magnets' ball-pair form: its terms on this group's balls. A pair with one ball outside the group: that ball held where it is now (its own island steps it)
      for (const Magnets::FormTerm &ft : Magnets::PairFormTerms())
      {
         int ia = -1, ib = -1;
         for (int k = 0; k < g.nb; k++)
         {
            if (balls[k] == ft.a) ia = k;
            if (ft.b != nullptr && balls[k] == ft.b) ib = k;
         }
         const HitBall *out = nullptr;
         if (ia < 0 && ib >= 0) { ia = ib; ib = -1; out = ft.a; }   // (a distance: either end may lead)
         else if (ia >= 0 && ib < 0 && ft.b != nullptr) out = ft.b;
         if (ia < 0)
            continue;
         PairForm::Term t = ft.t;
         t.a = ia;
         t.b = ib;
         if (out != nullptr)
         {
            V3 c((double)out->m_d.m_pos.x * U, (double)out->m_d.m_pos.y * U, (double)out->m_d.m_pos.z * U);
            for (const Body &w : cy.bodies)
               if (w.pb == out) { c = w.c; break; }
            t.o[0] = c.x; t.o[1] = c.y; t.o[2] = c.z;
         }
         g.AddForm(t);
      }
      const int nSurf = g.ns;   // the surfaces' shapes; after them each cord's posts (rubbers on posts: GroupStepDesign 13)
      for (HitRubber *const c : cords)
      {
         CordFollow(*c->Cord(), g.ball, balls.data(), g.nb);
         g.AddPart(c->Cord()->fixed.get());
         c->Cord()->fixed->AddPostShapes(g, c->Cord()->fixedTube);
         c->CordBegin();
         c->m_cordStepped = true;
      }
      g.Configure();
      // a coil's current changed since the last cycle: what the balls' places store changed with it, with nothing moved. That is the coil's work, a drive. Once per ball per physics step (a split cycle comes back here)
      for (const PullField &pf : pulls)
         if (pf.p.amps != pf.p.ampsBefore)
         {
            static std::unordered_map<const HitBall *, long long> booked;
            const long long stepNo = Magnets::StepNumber();
            for (int i = 0; i < g.nb; i++)
            {
               const GroupStep::Ball &b = g.ball[i];
               if (b.fieldScale == 0.0)
                  continue;
               long long &last = booked[balls[i]];
               if (last == stepNo)
                  continue;
               last = stepNo;
               for (const PullField &q : pulls)
                  g.books.coilJ += b.fieldScale * (q.At(q.p.amps, b.c, b.R, nullptr, nullptr) - q.At(q.p.ampsBefore, b.c, b.R, nullptr, nullptr));
            }
            break;
         }
      for (HitRubber *const c : cords)
         c->Cord()->fixed->RestorePostItems(g);
      // what a touch keeps between cycles: the shear spring's stretch
      for (int k = 0; k < g.ni; k++)
      {
         GroupStep::Item &it = g.item[k];
         if (it.b >= 0 || it.s >= nSurf)
            continue;
         const SurfaceShape *const sh = static_cast<const SurfaceShape *>((const GroupStep::Shape *)g.shape[it.s]);
         const SurfaceState &kept = sh->it->Surf(sh->surface);
         g.SetShapeItem(k, kept.touching, kept.zero, kept.xi, kept.ra);
      }
      // ... then each kept touch's depth now against its end depth, with every contact set (a massless part's angle is solved from all of them, as it was at the step's end)
      for (int k = 0; k < g.ni; k++)
      {
         GroupStep::Item &it = g.item[k];
         if (it.b >= 0 || it.s >= nSurf)
            continue;
         const SurfaceShape *const sh = static_cast<const SurfaceShape *>((const GroupStep::Shape *)g.shape[it.s]);
         const SurfaceState &kept = sh->it->Surf(sh->surface);
         if (ring && it.touching && g_contactKeptEnd(balls[it.a]))
         {
            GroupStep::Hit hp;
            g.ProbeItem(k, hp);
            const double was = sh->mp != nullptr ? kept.endDepth : sh->it->endDepth;
            if (hp.delta > was + 1e-9)
               it.dZero += hp.delta - was;
         }
         // what the zero took this cycle, by who moved the ball
         if (ring && it.dZero > kept.zero + 1e-9 && cause[it.a] != DepthCause::Placed && cause[it.a] != DepthCause::HandedOver)
            DepthFault(cause[it.a], balls[it.a], sh->r->Name(), it.dZero - kept.zero);
      }
      for (int k = 0; k < g.ni; k++)   // ... and a ball pair's zero
         if (g.item[k].b >= 0)
            if (const BallPairState *const ps = FindPair(balls[g.item[k].a], balls[g.item[k].b]))
               if (ps->zero > 0.0)
                  g.SetPairItem(k, ps->zero);
      static std::vector<V3> v0, vStep;
      v0.resize(g.nb);
      vStep.resize(g.nb);
      for (int i = 0; i < g.nb; i++)
         v0[i] = g.ball[i].v;
      // ring_model: a touch's hit event at the step it starts in, by the old loop's rule (a hit shorter than a cycle fires too); a ball pair keeps its stretch and fires its hit (14.8). The events wait for
      // the cycle to be solved: a failed one is retaken.
      static std::vector<char> hitOn;
      hitOn.resize(g.ni);
      static std::vector<PartHit> ringParts;
      static std::vector<BallHit> ringPairs;
      ringParts.clear();
      ringPairs.clear();
      if (ring || hasPart)
      {
         g.ClearGuess();
         vStep = v0;
         for (int k = 0; k < g.ni; k++)
         {
            GroupStep::Item &it = g.item[k];
            if (it.b < 0)
               hitOn[k] = it.s < nSurf && (static_cast<const SurfaceShape *>(g.shape[it.s])->mp != nullptr ? it.touching : static_cast<const SurfaceShape *>(g.shape[it.s])->it->touching);
            else if (ring)
            {
               const BallPairState *const ps = FindPair(balls[it.a], balls[it.b]);
               hitOn[k] = ps != nullptr && ps->inContact;
               if (hitOn[k] && it.touching)
               {
                  it.xi = ps->st.xi;
                  it.ra = ps->st.ra;
               }
            }
         }
      }
      const auto ringStep = [&]() {
         if (!ring && !hasPart)
            return;
         for (int k = 0; k < g.ni; k++)
         {
            const GroupStep::Item &gi = g.item[k];
            if (!ring && (gi.b >= 0 || gi.s >= nSurf || static_cast<const SurfaceShape *>(g.shape[gi.s])->mp == nullptr))   // (a cord island: a part's touches only)
               continue;
            if (gi.touching == hitOn[k])
               continue;
            hitOn[k] = gi.touching;
            if (!gi.touching)
               continue;
            if (gi.b >= 0)
            {
               const V3 d = g.ball[gi.a].c - g.ball[gi.b].c;
               const double dist = d.Length();
               const double closing = dist > 1e-12 ? -(vStep[gi.a] - vStep[gi.b]).Dot(d / dist) : 0.0;
               if (closing / VS > 0.25)   // as the stock collision: only impacts, not resting contact, reach the script
                  ringPairs.push_back({ balls[gi.a], balls[gi.b], (float)(closing / VS) });
            }
            else if (gi.s < nSurf)
            {
               const SurfaceShape *const sh = static_cast<const SurfaceShape *>(g.shape[gi.s]);
               GroupStep::Hit hp;
               sh->Probe(g.ball[gi.a].c, g.ball[gi.a].R, hp);
               if (sh->mp != nullptr)   // what its Hit is, the part decides; a touch that starts again within the release gap is the same touch
               {
                  SurfaceContact &sc = *const_cast<SurfaceContact *>(sh->it);
                  if (!sc.partTouching)
                  {
                     sc.partTouching = true;
                     sh->mp->TouchBegan(sh->surface, gi.a, -(vStep[gi.a] - g.SurfaceSpeed(*sh, hp)).Dot(hp.n));
                  }
                  continue;
               }
               const float vn = (float)(-vStep[gi.a].Dot(hp.n) / VS);
               const HitObject *const ho = sh->r->HitEventObject();
               if ((vn >= cy.hitFloor || (ho != nullptr && ho->m_fe && vn >= ho->m_threshold)) && sh->r->HitEventWanted(hp.n))
                  ringParts.push_back({ sh->r, balls[gi.a], vn, Vertex3Ds((float)hp.n.x, (float)hp.n.y, (float)hp.n.z) });
            }
         }
         for (int i = 0; i < g.nb; i++)
            vStep[i] = g.ball[i].v;
      };
      // VPX_CORD_TRACE=<file>: per cycle the group's energy at its start and end, and what happened to the ball on the cord
      static FILE *const cordTrace = []() -> FILE * {
         const char *const e = getenv("VPX_CORD_TRACE");
         FILE *const f = e != nullptr ? fopen(e, "w") : nullptr;
         if (f != nullptr)
            fprintf(f, "t_ms,dt_us,nb,np,E_start_mJ,E_end_mJ,loss_uJ,drive_uJ,book_uJ,x_U,y_U,z_U,vx,vy,vz,wx,wy,wz,on_cord,gap_mm,L_normal,items\n");
         return f;
      }();
      const double traceE0 = cordTrace != nullptr ? g.Energy() : 0.0;
      GroupStep::CycleBook book(g);
      double done = 0.0;
      // The group's step rule sees a surface's nearest face only, as a probe does: the other faces a ball can reach bound the step too (NextStep's rule). Not for energy, which the solve keeps
      // at any step: so that a hit on a hidden face is resolved, not one long step into it.
      const auto ringCap = [&](double h) {
         if (!ring)
            return h;
         const double a = g.set.g.Length();
         for (int s = 0; s < nSurf; s++)
         {
            const SurfaceShape &sh = *static_cast<const SurfaceShape *>(g.shape[s]);
            if (sh.it->candR < 0.0 || sh.it->NoFaces(g.ball[sh.onlyBall].c))
               continue;
            const GroupStep::Ball &b = g.ball[sh.onlyBall];
            const TubeProbe pr = sh.r->Probe(sh.bs, *sh.it, b.c);
            const double gap = pr.overlap ? -pr.delta : pr.dist - (b.R + pr.rho);
            if (gap <= b.v.Length() * h + a * h * h)
               h = std::max(sh.r->ArrivalS(sh.bs, *sh.it, pr.seg, b.c, b.v, a, h, 1.0e-6, 0.5e-6), std::min(h, 1e-6));
         }
         return h;
      };
      // for a test run (VPX_CORD_ENERGY): the group's own books over the cycle, stored energy included, told at each new largest rise
      static const bool cordEnergy = std::getenv("VPX_CORD_ENERGY") != nullptr;
      const bool cordEnergyOn = cordEnergy && !cords.empty();
      for (HitRubber *const c : cords)
         c->ReloadIfChanged();   // as a ring: the file is looked at while a ball is on it
      const double cordE0 = cordEnergyOn ? g.Energy() : 0.0;
      // a part's books (always on, rule 8): the group's energy across the cycle against what its drives did and its losses took
      const double partE0 = hasPart ? g.Energy() : 0.0;
      const double hFix = GroupStep::GroupFixedStepS();
      // the group's moving parts and cords by name
      const auto groupParts = [&]() {
         std::vector<std::string> names;
         for (const ContactMovingPart *const mp : parts)
            for (int s = 0; s < nSurf; s++)
               if (static_cast<const SurfaceShape *>(g.shape[s])->mp == mp)
               {
                  names.push_back(static_cast<const SurfaceShape *>(g.shape[s])->r->Name());
                  break;
               }
         for (HitRubber *const c : cords)
            names.push_back(c->Name());
         return names;
      };
      // Rule 8: each energy step with a moving part in the group: which parts, where, when, how long and why (why: GroupStep::Group's letters, each attempt of the step). The lines are capped per table start; the totals go at the table's end (ReleaseHost)
      long esSeen = g.books.energySteps, esBounces = g.books.energyBounces, esFree = g.books.energyFree;
      double esDone = 0.0;
      const auto partEnergy = [&]() {
         const double hTaken = done - esDone;
         esDone = done;
         if (!hasPart || g.books.energySteps == esSeen)
            return;
         const char *const how = g.books.energyFree > esFree ? "moved with no path" : g.books.energyBounces > esBounces ? "speeds turned back" : "move closed on energy";
         esSeen = g.books.energySteps;
         esBounces = g.books.energyBounces;
         esFree = g.books.energyFree;
         const std::vector<std::string> names = groupParts();
         for (const std::string &n : names)
            g_partEnergy[n]++;
         if (g_partEnergyLines >= kPartEnergyLines)
            return;
         g_partEnergyLines++;
         char s[260];
         std::string msg = "Group step: energy step, parts [";
         for (size_t a = 0; a < names.size(); a++)
            msg += (a > 0 ? ", " : "") + names[a];
         msg += "], balls " + std::to_string(g.nb);
         for (int i = 0; i < g.nb; i++)
         {
            snprintf(s, sizeof s, " (%.2f, %.2f, %.2f U, %.3f m/s)", g.ball[i].c.x / U, g.ball[i].c.y / U, g.ball[i].c.z / U, g.ball[i].v.Length());
            msg += s;
         }
         snprintf(s, sizeof s, ", t=%u ms + %.1f us, step %.3f us, %s, why", g_pplayer->m_time_msec, (done - hTaken) * 1e6, hTaken * 1e6, how);
         msg += s;
         for (int a = 0; a < g.nHalveDiag; a++)
         {
            const GroupStep::Group::HalveDiag &hd = g.halveDiag[a];
            snprintf(s, sizeof s, " [%c h %.3g us rounds %d move %.3g rim %.3g slipFlips %d touchFlips %d sliding %d]", hd.why, hd.h * 1e6, hd.rounds, hd.moveC, hd.moveW, hd.slipFlips, hd.touchFlips, hd.maxSliding);
            msg += s;
         }
         snprintf(s, sizeof s, " last %c, touching", g.failWhy);
         msg += s;
         for (int k = 0; k < g.ni; k++)
            if (g.item[k].touching)
            {
               GroupStep::Hit hp;
               if (g.item[k].b >= 0)
                  snprintf(s, sizeof s, " [ball %d on ball %d]", g.item[k].a, g.item[k].b);
               else
               {
                  g.ProbeItem(k, hp);
                  const SurfaceShape *const sh = g.item[k].s < nSurf ? static_cast<const SurfaceShape *>(g.shape[g.item[k].s]) : nullptr;
                  snprintf(s, sizeof s, " [ball %d %s/%d depth %.1f um]", g.item[k].a, sh != nullptr ? sh->r->Name().c_str() : "post", sh != nullptr ? sh->surface : 0, hp.delta * 1e6);
               }
               msg += s;
            }
         for (ContactMovingPart *const mp : parts)
         {
            msg += "; " + mp->StateNote();
            // each failed attempt's last round, as the part kept it (the halvings above, in order)
            if (GroupStep::Part *const gp = mp->GroupPart())
               for (int a = 0; a < g.nHalveDiag; a++)
               {
                  const std::string note = gp->FailNote(a);
                  if (!note.empty())
                  {
                     snprintf(s, sizeof s, " [%.3g us: ", g.halveDiag[a].h * 1e6);
                     msg += s + note + "]";
                  }
               }
         }
         msg += ". A net, not expected to fire: keep vpinball.log and report it";
         if (g_partEnergyLines == kPartEnergyLines)
            msg += " (no more of these lines this table start: the totals come at its end)";
         PLOGW << msg;
         // VPX_GROUP_HISTORY_TRACE=1: the accepted steps before the run's first energy step (the last ones marked ENERGY-STEP are its own tries, maybe undone)
         static bool historyDumped = false;
         if (GroupStep::Group::HistoryTraceOn() && !historyDumped)
         {
            historyDumped = true;
            PLOGW << "Group step: history, the last " << std::min<long>(g.historyN, GroupStep::Group::kHistory) << " accepted steps before the run's first energy step (energies uJ):\n" << g.HistoryDump();
         }
         // the state the energy step saved at the step's start, full precision, for an exact replay (balls in m, m/s, rad/s; contacts: their springs; parts: their own)
         snprintf(s, sizeof s, "Group step: energy step start state, t=%u ms + %.17g s, h %.17g s", g_pplayer->m_time_msec, done - hTaken, hTaken);
         msg = s + g.StartNote();
         for (size_t p = 0; p < parts.size(); p++)   // in the order of the names above
            if (GroupStep::Part *const gp = parts[p]->GroupPart())
               msg += " [part " + (p < names.size() ? names[p] : std::to_string(p)) + " " + gp->StartNote() + "]";
         PLOGW << msg;
         if (!g.failRoundBuf.empty())   // VPX_GROUP_ROUND_TRACE=1: the rounds of the step's first failed attempt
            PLOGW << "Group step: energy step, first failed attempt's rounds:\n" << g.failRoundBuf;
         // VPX_GROUP_CONTACT_TRACE=1: the first failed attempt's rounds, contact by contact and part by part, for the first kContactTraceSteps energy steps of the run; the contacts named first
         if (!g.failContactBuf.empty() && GroupStep::Group::ContactTraceLeft() > 0)
         {
            const int left = --GroupStep::Group::ContactTraceLeft();
            msg = "Group step: contact trace, t=" + std::to_string(g_pplayer->m_time_msec) + " ms, first failed attempt's rounds (" + std::to_string(left) + " more steps to print). Contacts:";
            for (int k = 0; k < g.ni; k++)
            {
               const GroupStep::Item &it = g.item[k];
               if (it.b >= 0)
                  snprintf(s, sizeof s, " [c%d ball %d on ball %d]", k, it.a, it.b);
               else
               {
                  const SurfaceShape *const sh = it.s < nSurf ? static_cast<const SurfaceShape *>(g.shape[it.s]) : nullptr;
                  snprintf(s, sizeof s, " [c%d ball %d on s%d %s/%d part p%d]", k, it.a, it.s, sh != nullptr ? sh->r->Name().c_str() : "post", sh != nullptr ? sh->surface : 0, g.shape[it.s]->part);
               }
               msg += s;
            }
            PLOGW << msg << "\n" << g.failContactBuf;
            g.failContactBuf.clear();
         }
      };
      // VPX_GROUP_STATS=<file> (a diagnostic, read once): the group step's totals over the run, rewritten every 200 ms of table time: cycles, time per cycle, steps, solve rounds, halvings,
      // energy steps, failed cycles (the behaviour check's step statistics)
      static const char *const statsPath = []() -> const char * { const char *const e = getenv("VPX_GROUP_STATS"); return e != nullptr && *e ? e : nullptr; }();
      const uint64_t statsT0 = statsPath != nullptr ? usec() : 0;
      if (!g.jacWinBuf.empty())   // VPX_GROUP_JAC_WINDOW: the last cycle's converged solves' rounds
      {
         PLOGI << "Group step: Jacobian window:\n" << g.jacWinBuf;
         g.jacWinBuf.clear();
      }
      if (GroupStep::Group::HistoryTraceOn())   // (diagnostics: the history's step times)
      {
         static long cycles = 0;
         g.histMs = g_pplayer->m_time_msec;
         g.histClock = done;
         g.histCycle = ++cycles;
      }
      // VPX_STEP_RULE=<file> (a diagnostic, read once; off: nothing changes): one row per step: what ChooseStep's rule was (StepDiag), its h and the h after ringCap, the solve rounds and the wall time of the step before
      static FILE *const ruleF = []() -> FILE * {
         const char *const e = getenv("VPX_STEP_RULE");
         FILE *const o = e != nullptr && *e ? fopen(e, "w") : nullptr;
         if (o != nullptr)
            fprintf(o, "t_ms,done_us,h_us,hcap_us,rule,partRule,item,itemPart,touching,b,delta_um,dDot,omN,prevRounds,prevUs,envSolves,envIters,capsules,envHits,envMisses,envUs,unsteady,unsteadyVal,unsteadyLim,A_um,missProbe,missEnd,missMid,missStart,missRate,missAngle,missEnergy,term\n");
         return o;
      }();
      long ruleRounds = g.books.rounds;
      auto ruleT = std::chrono::steady_clock::now();
      const auto ruleCap = [&](const double h) {
         const double hc = hFix > 0.0 ? std::min(hFix, cy.dt - done) : ringCap(h);
         if (ruleF != nullptr)
         {
            const auto now = std::chrono::steady_clock::now();
            const GroupStep::Group::StepDiag &dg = g.diag;
            const bool hasItem = dg.rule != 8 && dg.item >= 0 && dg.item < g.ni;
            const GroupStep::EnvStats es = GroupStep::envStats;   // (the massless parts' solves since the last row: the step before and this choice)
            GroupStep::envStats = GroupStep::EnvStats();
            fprintf(ruleF, "%u,%.3f,%.4f,%.4f,%d,%d,%d,%d,%d,%d,%.4f,%.5f,%.1f,%ld,%.1f,%ld,%ld,%ld,%ld,%ld,%.1f,%d,%.6g,%.6g,%.4g,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%d\n", g_pplayer->m_time_msec, done * 1e6, h * 1e6, hc * 1e6, dg.rule, dg.partRule, dg.item,
                    hasItem && g.item[dg.item].b < 0 ? g.shape[g.item[dg.item].s]->part : -9, hasItem ? (int)g.item[dg.item].touching : -1, hasItem ? g.item[dg.item].b : -9, dg.delta * 1e6, dg.dDot, dg.omN,
                    g.books.rounds - ruleRounds, std::chrono::duration<double, std::micro>(now - ruleT).count(), es.solves, es.iters, es.capsules, es.hits, es.misses, es.us, dg.unsteady, dg.unsteadyVal, dg.unsteadyLim, dg.aShear * 1e6,
                    es.missBy[0], es.missBy[1], es.missBy[2], es.missBy[3], es.missBy[4], es.missBy[5], es.missBy[6], dg.term);
            ruleRounds = g.books.rounds;
            ruleT = now;
         }
         return hc;
      };
      g.phaseMs = g_pplayer->m_time_msec;   // (VPX_GROUP_PHASE_TIMING rows)
      const bool ok = g.StepsOver(cy.dt, done, ruleCap, [&]() { ringStep(); book.AfterStep(); partEnergy(); });
      if (statsPath != nullptr)
      {
         static long n = 0, steps = 0, rounds = 0, eSteps = 0, fails = 0, halv = 0;
         static double us = 0.0, usMax = 0.0;
         static uint32_t lastDump = 0;
         const double t = (double)(usec() - statsT0);
         n++;
         us += t;
         usMax = std::max(usMax, t);
         steps += g.books.steps;
         rounds += g.books.rounds;
         eSteps += g.books.energySteps;
         halv += g.books.halvings;
         fails += ok ? 0 : 1;
         if (g_pplayer->m_time_msec - lastDump > 200)
         {
            lastDump = g_pplayer->m_time_msec;
            if (FILE *const f = fopen(statsPath, "w"))
            {
               fprintf(f, "cycles %ld us/cycle %.2f max %.0f steps %ld rounds %ld halvings %ld energysteps %ld failed %ld balls %d\n", n, us / n, usMax, steps, rounds, halv, eSteps, fails, g.nb);
               fclose(f);
            }
         }
      }
      if (cordEnergyOn)
      {
         static double worstRise = 1e-9;
         const double riseE = g.Energy() - cordE0 + g.books.Losses() - g.books.driveJ;
         if (riseE > worstRise)
         {
            worstRise = riseE;
            PLOGI << "Cord energy: largest rise over a cycle so far " << riseE * 1e6 << " uJ, t=" << g_pplayer->m_time_msec;
         }
      }
      book.End((int)g_pplayer->m_time_msec, false, cy.dt, !ok);
      for (ContactMovingPart *const mp : parts)
      {
         mp->steppedS += done;
         mp->steppedThisCycle = true;
      }
      if (cordTrace != nullptr)
      {
         const double E1 = g.Energy(), loss = g.books.Losses(), drive = g.books.driveJ;
         const CordBand *const pc = cords.empty() ? nullptr : cords[0]->Cord();
         fprintf(cordTrace, "%d,%.2f,%d,%d,%.6f,%.6f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.4f,%.4f,%.4f,%.2f,%.2f,%.2f,%d,%.3f,%.2f,%s\n", (int)g_pplayer->m_time_msec, cy.dt * 1e6, g.nb, g.np,
                 traceE0 * 1e3, E1 * 1e3, loss * 1e6, drive * 1e6, (E1 - traceE0 + loss - drive) * 1e6, g.ball[0].c.x / U, g.ball[0].c.y / U, g.ball[0].c.z / U, g.ball[0].v.x, g.ball[0].v.y, g.ball[0].v.z,
                 g.ball[0].w.x, g.ball[0].w.y, g.ball[0].w.z, pc != nullptr && pc->OnCord(0) ? 1 : 0, pc != nullptr ? pc->Gap(0) * 1e3 : 0.0, g.books.normalLossJ * 1e6, [&]() {
                    static std::string names;
                    names.clear();
                    for (int k = 0; k < g.ni; k++)
                       if (g.item[k].touching && g.item[k].b < 0)
                          names += g.item[k].s < nSurf ? static_cast<const SurfaceShape *>(g.shape[g.item[k].s])->r->Name() + "|" : std::string("post|");
                    return names.c_str();
                 }());
         fflush(cordTrace);
      }
      if (ring && !ok && !hasPart)   // a net, not the fix: on screen once (rule 8)
      {
         static uint32_t lastMs = 0;
         static long failed = 0;
         static const void *warned = nullptr;
         failed++;
         if (g_pplayer->m_time_msec - lastMs > 500)
         {
            lastMs = g_pplayer->m_time_msec;
            PLOGW << "Ring group failed (" << failed << " in all): why " << g.failWhy << " at " << done * 1e6 << " us of " << cy.dt * 1e6 << " us, balls " << g.nb << ", ball 0 at (" << g.ball[0].c.x / U << ", " << g.ball[0].c.y / U << ", " << g.ball[0].c.z / U << ") v " << g.ball[0].v.Length() << " m/s: the cycle is retaken on the old loop";
         }
         if (warned != g_pplayer)
         {
            warned = g_pplayer;
            const std::string msg = "Rubbers: a group step could not be solved: the cycle was retaken on the old loop (see vpinball.log)";
            PLOGW << msg;
            g_pplayer->m_liveUI->PushNotification(msg, 6000);
         }
         return false;
      }
      int ib = 0;
      for (const int i : isl.members)
      {
         if (ib >= g.nb)
            break;
         Body &w = cy.bodies[i];
         if (ring && w.e != nullptr)
            w.e->seatPending = false;   // (a surface it was found inside counted from zero)
         w.c = g.ball[ib].c;
         w.v = g.ball[ib].v;
         w.w = g.ball[ib].w;   // (a step that could not be solved: what is left of the cycle is not stepped, so nothing is gained: below)
         {
            KeptDoubles &kd = g_keptDoubles[w.pb];
            kd.pos = Vertex3Ds((float)(w.c.x / U), (float)(w.c.y / U), (float)(w.c.z / U));   // (as the cycle's end writes them)
            kd.vel = Vertex3Ds((float)(w.v.x / VS), (float)(w.v.y / VS), (float)(w.v.z / VS));
            kd.am = Vertex3Ds((float)(w.w.x * w.lScale), (float)(w.w.y * w.lScale), (float)(w.w.z * w.lScale));
            kd.c = w.c;
            kd.v = w.v;
            kd.w = w.w;
         }
         ib++;
      }
      // the part's hit events, and what the touches keep
      for (int k = 0; k < g.ni; k++)
      {
         const GroupStep::Item &gi = g.item[k];
         if (gi.b >= 0 && (gi.touching || gi.dZero > 0.0))
         {
            BallPairState &ps = MakePair(balls[gi.a], balls[gi.b]);
            if (ring)
            {
               ps.st.xi = gi.xi;
               ps.st.ra = gi.ra;
            }
            ps.zero = gi.dZero;
            ps.inContact = gi.touching;
            ps.used = true;
         }
         if (gi.b >= 0 || gi.s >= nSurf)
            continue;
         const SurfaceShape *const sh = static_cast<const SurfaceShape *>(g.shape[gi.s]);
         SurfaceContact &it = *const_cast<SurfaceContact *>(sh->it);
         SurfaceState &kept = it.Surf(sh->surface);
         if (sh->mp != nullptr)
         {
            it.touching = sh->surface == 0 ? gi.touching : it.touching || gi.touching;
            kept.reff = gi.touching ? gi.reff : 0.0;
            kept.kb = gi.touching ? gi.kb : 0.0;
            GroupStep::Hit hp;
            g.ProbeItem(k, hp);
            kept.endDepth = hp.delta;
            if (!gi.touching && it.partTouching && hp.delta < -sh->r->ReleaseGapM())
               it.partTouching = false;
         }
         else if (ring)
         {
            it.gs.reff = gi.touching ? gi.reff : 0.0;
            it.touching = hitOn[k];
         }
         else if (gi.touching && !it.touching)
         {
            GroupStep::Hit hp;
            sh->Probe(g.ball[gi.a].c, g.ball[gi.a].R, hp);
            const float vn = (float)(-v0[gi.a].Dot(hp.n) / VS);
            const HitObject *const ho = it.r->HitEventObject();
            if ((vn >= cy.hitFloor || (ho != nullptr && ho->m_fe && vn >= ho->m_threshold)) && it.r->HitEventWanted(hp.n))
               cy.partHits.push_back({ it.r, balls[gi.a], vn, Vertex3Ds((float)hp.n.x, (float)hp.n.y, (float)hp.n.z) });
         }
         if (!ring && sh->mp == nullptr)
            it.touching = gi.touching;
         kept.touching = gi.touching;
         kept.xi = gi.xi;
         kept.ra = gi.ra;
         kept.zero = gi.dZero;
      }
      for (HitRubber *const c : cords)
      {
         CordBand &cb = *c->Cord();
         std::vector<HitBall *> &on = cb.onCord;
         for (int i = 0; i < g.nb; i++)
            if (cb.OnCord(i) && std::find(on.begin(), on.end(), balls[i]) == on.end())
               cb.hits.emplace_back(balls[i], (float)(cb.Closing(i, v0[i]) / VS));
         on.clear();
         for (int i = 0; i < g.nb; i++)
            if (cb.OnCord(i))
               on.push_back(balls[i]);
         if (g.books.energySteps > 0)   // none expected: told each time (rule 8)
         {
            cb.fixed->energySteps += g.books.energySteps;
            PLOGW << "Rubber " << c->Name() << ": energy step on the fixed-height cord, " << cb.fixed->energySteps << " so far, t=" << g_pplayer->m_time_msec;
            if (g_pplayer->m_liveUI != nullptr)
               g_pplayer->m_liveUI->PushNotification("Rubber " + c->Name() + ": energy step on the fixed-height cord", 6000);
         }
      }
      if (ring)
      {
         cy.partHits.insert(cy.partHits.end(), ringParts.begin(), ringParts.end());
         cy.hitEvents.insert(cy.hitEvents.end(), ringPairs.begin(), ringPairs.end());
         cy.steps = (int)g.books.steps;
         g_last.substeps += cy.steps;
         g_last.maxSteps = std::max(g_last.maxSteps, cy.steps);
      }
      if (!ok)
      {
         PLOGW << "Cord group failed: why " << g.failWhy << " at " << done * 1e6 << " us of " << cy.dt * 1e6 << " us, balls " << g.nb << ", ball 0 at (" << g.ball[0].c.x / U << ", " << g.ball[0].c.y / U << ", " << g.ball[0].c.z / U << ") v " << g.ball[0].v.Length() << " m/s, halvings " << g.books.halvings << ", plate " << (cords.empty() ? 0.0 : cords[0]->Cord()->fixed->phi * 57.29578) << " deg";
         {
            std::string how;
            char s[220];
            for (int a = 0; a < g.nHalveDiag; a++)
            {
               const GroupStep::Group::HalveDiag &hd = g.halveDiag[a];
               snprintf(s, sizeof s, " [%c h %.3g us rounds %d move %.3g rim %.3g slipFlips %d touchFlips %d sliding %d]", hd.why, hd.h * 1e6, hd.rounds, hd.moveC, hd.moveW, hd.slipFlips, hd.touchFlips, hd.maxSliding);
               how += s;
            }
            for (int k = 0; k < g.ni; k++)
               if (g.item[k].touching)
               {
                  snprintf(s, sizeof s, " item a%d b%d s%d%s zero %.3g", g.item[k].a, g.item[k].b, g.item[k].s, g.item[k].b < 0 && g.item[k].s >= nSurf ? " post" : "", g.item[k].dZero);
                  how += s;
               }
            for (int i = 0; i < g.nb; i++)
            {
               snprintf(s, sizeof s, " ball%d v (%.4f %.4f %.4f) w (%.2f %.2f %.2f)", i, g.ball[i].v.x, g.ball[i].v.y, g.ball[i].v.z, g.ball[i].w.x, g.ball[i].w.y, g.ball[i].w.z);
               how += s;
            }
            for (HitRubber *const c : cords)
            {
               const CordFixed::Part &p = *c->Cord()->fixed;
               snprintf(s, sizeof s, " cord %s plate %.3f omega %.3f on0 %d", c->Name().c_str(), p.phi * 57.29578, p.omega, p.OnCord(0) ? 1 : 0);
               how += s;
            }
            PLOGW << "Cord group how:" << how << " (cycle lost: " << g.overWhy << ")";
            if (!g.stopBuf.empty())   // VPX_GROUP_STOP_TRACE=1
               PLOGW << "Cord group stop checks this cycle:\n" << g.stopBuf;
         }
         // A hard step (GroupStepDesign 15): the balls stay as they are for the rest of the cycle and the next cycle tries again from the same state: time is lost, no energy can appear,
         // and nothing goes to stock. If it will not solve for kHeldCycles in a row the cord lets the balls go, which only takes energy away.
         constexpr int kHeldCycles = 20;
         for (HitRubber *const c : cords)
         {
            CordBand &cb = *c->Cord();
            if (!cb.warnedHeld)
            {
               cb.warnedHeld = true;
               const std::string msg = "Rubber " + c->Name() + ": a step could not be solved: the ball was held for the rest of the cycle";
               PLOGW << msg;
               g_pplayer->m_liveUI->PushNotification(msg, 6000);
            }
            if (++cb.heldCycles >= kHeldCycles)
            {
               cb.heldCycles = 0;
               c->CordSolveFailed("the group step could not solve a step for 20 cycles: the ball was let go");
            }
         }
      }
      else
         for (HitRubber *const c : cords)
         {
            c->Cord()->heldCycles = 0;
         }
      WriteBack(cy, isl);
      if (hasPart)   // rule 8
      {
         static const void *warnedFail = nullptr, *warnedBook = nullptr, *warnedGain = nullptr;
         const double gain = g.Energy() - partE0 + g.books.Losses() - g.books.driveJ;
         if (!ok && warnedFail != g_pplayer)
         {
            warnedFail = g_pplayer;
            const std::string msg = "Group step: a step with a moving part could not be solved: the rest of that cycle was lost (see vpinball.log)";
            PLOGW << msg;
            g_pplayer->m_liveUI->PushNotification(msg, 6000);
         }
         if (g.books.energySteps > 0 && cords.empty() && warnedBook != g_pplayer)
         {
            warnedBook = g_pplayer;
            std::string msg = "Group step: a step with";
            for (const std::string &n : groupParts())
               msg += " " + n;
            msg += " was taken on energy alone: a net, not expected to fire. Keep vpinball.log and report it";
            g_pplayer->m_liveUI->PushNotification(msg, 6000);
         }
         if (gain > 1e-4 && warnedGain != g_pplayer)
         {
            warnedGain = g_pplayer;
            const std::string msg = "Group step: a moving part's group gained " + std::to_string(gain * 1e3) + " mJ beyond its drives in one cycle, t=" + std::to_string(g_pplayer->m_time_msec);
            PLOGW << msg;
            g_pplayer->m_liveUI->PushNotification(msg, 6000);
         }
      }
      return true;
   }

   // a cord with no ball this cycle: a plate moving by itself, or the path letting go of a ball just gone (GroupStepDesign 8 point 5)
   static void CordAlone(HitRubber *const band, const float dtime)
   {
      CordBand &cb = *band->Cord();
      const double dt = (double)dtime * DEFAULT_STEPTIME_S;
      if (dt <= 0.0)
         return;
      if (cb.attached && cb.slots.empty() && cb.PlateQuiet() && !band->SlingTriggerPending())
      {
         band->CordBegin();   // (the switches, the driver: nothing to step unless they start something)
         if (cb.PlateQuiet())
            return;
      }
      CordFollow(cb, nullptr, nullptr, 0);
      if (!(cb.attached && cb.slots.empty() && !cb.PlateQuiet()))
         band->CordBegin();
      if (!GroupStep::AdvanceAlone(*cb.fixed, dt, 2e-6, (int)g_pplayer->m_time_msec))
         band->CordSolveFailed("the plate's step could not be solved");
   }

   static void IslandMembers(const Cycle &cy, const ContactSystem::Islands::Island &isl, std::string &out)
   {
      std::vector<std::string> balls, parts;
      auto add = [&parts](const std::string &n) { if (parts.size() < 10 && std::find(parts.begin(), parts.end(), n) == parts.end()) parts.push_back(n); };
      for (const int i : isl.members)
      {
         const Body &w = cy.bodies[i];
         const char *const n = TestDriver::BallName(w.pb);
         balls.push_back(n != nullptr ? n : "ball");
         if (w.e != nullptr)
            for (const SurfaceContact &it : w.e->items)
               if (it.touching)
                  add(it.r->Name());
         if (w.band != nullptr)
            add("band " + w.band->Name());
      }
      out = "balls";
      for (const std::string &b : balls)
         out += ' ' + b;
      out += " | touching:";
      for (const std::string &p : parts)
         out += ' ' + p + ';';
      if (parts.empty())
         out += " nothing";
   }

   // Releases and the pair clean-up of the step, after the free balls have moved and the overlap check has passed
   static void ApplyReleases(Cycle &cy)
   {
      for (const auto &r : cy.releases)
         Release(r.first, r.second);
      for (const auto &r : cy.bandReleases)
         if (r.first->m_bandOwner != nullptr)
            r.first->m_bandOwner->BandRelease(r.first, r.second);
      cy.releases.clear();
      cy.bandReleases.clear();
      for (size_t i = 0; i < g_ballPairs.size();)
         if (!g_ballPairs[i].used)
            g_ballPairs.erase(g_ballPairs.begin() + (long)i);
         else
            i++;
   }

   // The events of the step, then those the parts queued, in the order they happened (a script may move or delete balls: the queue is taken first)
   static void DispatchEvents(Cycle &cy)
   {
      CPROF(Finish);
      for (const BallHit &hv : cy.hitEvents)
      {
         g_pplayer->m_ptable->InvokeBallBallCollisionCallback(hv.a, hv.b, hv.speed);
         g_pplayer->m_pininput.PlayBallBallRumble(hv.speed);
      }
      for (const PartHit &ph : cy.partHits)
      {
         const auto &vb = g_pplayer->m_vball;
         if (std::none_of(vb.begin(), vb.end(), [&](const Ball *b) { return &b->m_hitBall == ph.pb; }))
            continue;   // an earlier event deleted the ball
         g_pplayer->m_pactiveball = ph.pb->m_pBall;
         g_pplayer->m_physics->OnBallWallHit(*ph.pb, ph.normal, ph.speed);
         if (ph.s->IsPlayfield())
         {
            const PhysicsEngine &pe = *g_pplayer->m_physics;
            if (pe.HitSpeedEvents())
            {
               ph.pb->m_hitSpeed = ph.speed;
               if (ph.speed >= pe.HitSpeedMin())
                  g_pplayer->m_ptable->FireVoidEvent(DISPID_HitEvents_Hit);
            }
         }
         if (ContactTraceOn())
            PLOGI << "PartHit " << ph.s->Name() << " t=" << g_pplayer->m_time_msec << " speed_ms " << ph.speed * VS << " event " << (ph.s->HitEventObject() != nullptr && ph.s->HitEventObject()->m_fe ? 1 : 0)
                  << " thr " << (ph.s->HitEventObject() != nullptr ? ph.s->HitEventObject()->m_threshold : 0.f);
         if (HitObject *const ho = ph.s->HitEventObject(); ho != nullptr && ho->m_obj != nullptr && ho->m_fe && ph.speed >= ho->m_threshold)
         {
            ho->m_obj->m_currentHitThreshold = ph.speed;
            ho->FireHitEvent(ph.pb);
         }
      }
      cy.hitEvents.clear();
      cy.partHits.clear();
      DrainDeferred();
   }
   static void DrainDeferred()
   {
      if (!g_deferred.empty())
      {
         std::vector<std::function<void()>> q;
         q.swap(g_deferred);
         for (auto &f : q)
            f();
      }
   }
   static Cycle *s_pendingCy;   // the step whose events have not been dispatched yet
   static std::vector<ContactHost *> &Hosts() { return ContactHost::Registry(); }   // (for the snapshot, outside the struct)
   // ---- Pile sleep ([Contact] sleep = 1; runs/2026-11-03-pile-sleep/design.md) ----
   // A group (an island) whose balls have all stayed under sleep_speed_mm_s, in travel and in the surface speed of their spin, for sleep_hold_ms, touching only parts at rest, with every
   // coil current steady and no contact deeper than 0.05 mm, sleeps: its balls keep their place and every contact state (surface shear, ball pairs), their speed and spin go to 0 (a loss
   // of at most the kinetic energy under the line; nothing stored changes), and nothing steps them until a wake. The line is requirement 4.1's 'still' (under 0.1 mm/s).
   struct SleepGroup
   {
      int id = 0;
      std::vector<HitBall *> balls;
      std::vector<ContactMovingPart *> parts;   // the moving parts its balls touch
      std::vector<char> partIdle;              // the part has been seen asleep since: a drive that wakes it then wakes the group
   };
   inline static std::vector<SleepGroup> g_sleep;
   inline static int g_sleepNextId = 1;
   inline static Vertex3Ds g_sleepBase;          // gravity and nudge over one physics step at the last cycle start (U/VPT)
   inline static bool g_sleepBaseKnown = false;
   inline static double g_sleepCycleDt = 0.0;    // the share of this physics cycle the contact steps have covered (s): holds count at the cycle's end only
   inline static int g_sleepLog = 0;             // sleep and wake lines this table start
   inline static bool g_sleepCreepWarned = false, g_sleepDepthWarned = false;
   inline static const void *g_sleepOwner = nullptr;
   static constexpr double kSleepDepthM = 0.05e-3;   // requirement 4.2: a contact deeper than this is a fault to see, not a pile to put to sleep
   static constexpr double kSleepCreepMS = 1e-3;     // the slow-creep notice: a group between the sleep line and 1 mm/s ...
   static constexpr double kSleepCreepS = 5.0;       // ... for this long

   static void SleepTable()
   {
      if (g_sleepOwner == g_pplayer)
         return;
      g_sleepOwner = g_pplayer;
      g_sleep.clear();
      g_sleepNextId = 1;
      g_sleepBaseKnown = false;
      g_sleepCycleDt = 0.0;
      g_sleepLog = 0;
      g_sleepCreepWarned = g_sleepDepthWarned = false;
   }
   static bool SleepLogOk()
   {
      if (g_sleepLog < 1000)
      {
         g_sleepLog++;
         return true;
      }
      if (g_sleepLog++ == 1000)
         PLOGI << "Pile sleep: further sleep and wake lines this table start are not logged";
      return false;
   }
   static void SleepNotice(const std::string &msg)
   {
      PLOGW << msg;
      if (g_pplayer != nullptr && g_pplayer->m_liveUI != nullptr)
         g_pplayer->m_liveUI->PushNotification(msg, 8000);
   }
   static bool SameV(const Vertex3Ds &a, const Vertex3Ds &b) { return a.x == b.x && a.y == b.y && a.z == b.z; }
   static bool PartQuiet(const ContactMovingPart *mp) { return mp->Sleeping() || mp->AtRest(); }

   // owe: woken after this cycle's velocity update (which gave a sleeping ball no gravity): the contact steps of the rest of the cycle give the share back first (SleepInStep)
   static void SleepWakeGroup(const size_t gi, const char *cause, const bool owe)
   {
      SleepGroup &g = g_sleep[gi];
      if (SleepLogOk())
         PLOGI << "Pile sleep: wake group " << g.id << " (" << g.balls.size() << " balls) at t=" << g_pplayer->m_time_msec << " ms: " << cause;
      for (HitBall *const pb : g.balls)
      {
         pb->m_sleepGroup = 0;
         pb->m_sleepHoldS = 0.0;
         pb->m_sleepCreepS = 0.0;
         pb->m_sleepOwe = owe;
      }
      g_sleep.erase(g_sleep.begin() + (long)gi);
   }
   static void SleepWakeBall(const HitBall *pb, const char *cause, const bool owe)
   {
      if (pb == nullptr || pb->m_sleepGroup == 0)
         return;
      for (size_t gi = 0; gi < g_sleep.size(); gi++)
         if (g_sleep[gi].id == pb->m_sleepGroup)
         {
            SleepWakeGroup(gi, cause, owe);
            return;
         }
   }
   // Something no step did: a ball moved, given a speed or spin, captured, or put somewhere (a script, the driver, a kicker, a stock collision)
   static const char *SleepMoved(const SleepGroup &g)
   {
      for (const HitBall *const pb : g.balls)
      {
         if (pb->m_d.m_lockedInKicker)
            return "ball captured by a kicker";
         if (!SameV(pb->m_d.m_pos, pb->m_sleepPos))
            return "ball moved outside the step";
         if (!SameV(pb->m_d.m_vel, Vertex3Ds(0.f, 0.f, 0.f)) || !SameV(pb->m_angularmomentum, Vertex3Ds(0.f, 0.f, 0.f)))
            return "ball given a speed or spin outside the step";
         if (std::find(g_placed.begin(), g_placed.end(), pb) != g_placed.end())
            return "ball placed";
      }
      return nullptr;
   }
   // A part the balls touch that a drive woke or that moved; a part not at rest that can reach a ball of the group within the cycle (the rule of PartsInSolve)
   static const char *SleepParts(SleepGroup &g, const double dt)
   {
      for (size_t k = 0; k < g.parts.size(); k++)
      {
         ContactMovingPart *const mp = g.parts[k];
         if (mp->Sleeping())
            g.partIdle[k] = 1;
         else if (g.partIdle[k])
            return "a part it touches was driven";
         else if (!mp->AtRest())
            return "a part it touches moved";
      }
      for (ContactHost *const o : ContactHost::Registry())
      {
         ContactMovingPart *const mp = o->MovingPart();
         if (mp == nullptr || PartQuiet(mp))
            continue;
         for (HitBall *const pb : g.balls)
         {
            Body s;
            LoadBody(s, pb);
            LoadState(s);
            mp->SetBall(s.bs);
            const double reach = 2.0 * mp->SurfaceSpeedBound(dt) * dt + 1e-4;
            if ((s.c - mp->Origin()).Length() > mp->Reach() + s.bs.R + reach)
               continue;
            if (mp->Probe(s.c).delta > -reach)
               return "a moving part can reach it";
         }
      }
      return nullptr;
   }
   // Once per cycle, before the magnets and the balls' velocity update: gravity and nudge, coil currents, moves outside the step, parts
   static void SleepStart()
   {
      if (!ContactSettings::Get().sleep || g_pplayer == nullptr)
         return;
      SleepTable();
      g_sleepCycleDt = 0.0;
      for (auto *const ball : g_pplayer->m_vball)
         ball->m_hitBall.m_sleepOwe = false;
      if (g_pplayer->m_vball.empty())
         return;
      const Vertex3Ds base = g_pplayer->m_vball[0]->m_hitBall.ExternalVelocityIncrement(false);
      const bool baseChanged = g_sleepBaseKnown && !SameV(base, g_sleepBase);
      g_sleepBase = base;
      g_sleepBaseKnown = true;
      const bool steady = Magnets::CurrentSteady();
      if (baseChanged || !steady)   // a hold counts only under the same gravity, nudge and currents
         for (auto *const ball : g_pplayer->m_vball)
            if (baseChanged || ball->m_hitBall.m_magnetic)
               ball->m_hitBall.m_sleepHoldS = 0.0;
      for (size_t gi = 0; gi < g_sleep.size();)
      {
         SleepGroup &g = g_sleep[gi];
         const char *cause = baseChanged ? "gravity or nudge changed" : nullptr;
         if (cause == nullptr && !steady)
            for (const HitBall *const pb : g.balls)
               if (pb->m_magnetic)
                  cause = "a coil current changes";
         if (cause == nullptr)
            cause = SleepMoved(g);
         if (cause == nullptr)
            cause = SleepParts(g, PHYSICS_STEPTIME_S);
         if (cause != nullptr)
            SleepWakeGroup(gi, cause, false);
         else
            gi++;
      }
   }
   // At each contact step, before the gather: moves outside the step since the last one, an awake ball within reach (the islands' own test), then the gravity share a ball woken
   // within this cycle is owed; the pairs of sleeping balls are kept
   static void SleepInStep(Cycle &cy)
   {
      SleepTable();
      if (!s_append)
         g_sleepCycleDt += cy.dt;
      for (size_t gi = 0; gi < g_sleep.size();)
         if (const char *const cause = SleepMoved(g_sleep[gi]))
            SleepWakeGroup(gi, cause, true);
         else
            gi++;
      for (bool woke = !g_sleep.empty(); woke;)
      {
         woke = false;
         // ball_pair_reach: the islands' test (BallsMayMeet) with the groups not built yet: the energy of every awake ball can reach an awake ball, J (spin included)
         const bool pairReach = cy.cs->ballPairReach && cy.cs->islandReach;
         double eAwake = 0.0;
         if (pairReach)
            for (auto *const a : g_pplayer->m_vball)
            {
               const HitBall &h = a->m_hitBall;
               if (h.m_sleepGroup == 0 && !h.m_d.m_lockedInKicker)
                  eAwake += 0.5 * (h.m_d.m_mass * ContactSystem::kMassUnitKg * h.m_d.m_vel.Dot(h.m_d.m_vel) * VS * VS + h.m_angularmomentum.Dot(h.m_angularmomentum) / h.Inertia() * ContactSystem::kMassUnitKg * VS * VS);
            }
         for (auto *const a : g_pplayer->m_vball)
         {
            const HitBall *const pa = &a->m_hitBall;
            if (pa->m_sleepGroup != 0 || pa->m_d.m_lockedInKicker)
               continue;
            const double ra = (double)pa->m_d.m_radius * U, va = std::sqrt((double)pa->m_d.m_vel.Dot(pa->m_d.m_vel)) * VS;
            const double va2 = pairReach ? BallSpeedBound(eAwake, pa->m_d.m_mass * ContactSystem::kMassUnitKg, 0.0) : 0.0;
            for (auto *const b : g_pplayer->m_vball)
            {
               const HitBall *const pb = &b->m_hitBall;
               if (pb->m_sleepGroup == 0)
                  continue;
               const Vertex3Ds d = pa->m_d.m_pos - pb->m_d.m_pos;
               const double gap = std::sqrt((double)d.Dot(d)) * U - ra - (double)pb->m_d.m_radius * U;
               if (gap < cy.cs->ballReachMm * 1e-3 + 2.0 * va * cy.dt || (pairReach && BallsMayMeet(gap, va2, 0.0, cy.dt)))   // the islands' pairing: the pairwise rule, or the closure (sleeper V = 0)
               {
                  SleepWakeBall(pb, "an awake ball within reach", true);
                  woke = true;
               }
            }
         }
      }
      const double share = cy.dt / PHYSICS_STEPTIME_S;
      for (auto *const ball : g_pplayer->m_vball)
      {
         HitBall &hb = ball->m_hitBall;
         if (hb.m_sleepOwe && !s_append)
            hb.m_d.m_vel += hb.ExternalVelocityIncrement() * (float)share;
      }
      for (BallPairState &p : g_ballPairs)
         if (p.a != nullptr && p.b != nullptr && p.a->m_sleepGroup != 0 && p.b->m_sleepGroup != 0)
            p.used = true;
   }
   // After an island's step: its balls' holds at the cycle's end, the slow-creep notice, and the sleep itself
   static void SleepConsider(Cycle &cy, const ContactSystem::Islands::Island &isl)
   {
      const ContactSettings &cs = *cy.cs;
      if (!cs.sleep || s_append || g_sleepCycleDt < PHYSICS_STEPTIME_S * (1.0 - 1e-6))   // (between the sub-intervals of a split cycle a ball still carries the rest of the cycle's gravity)
         return;
      const double vs = cs.sleepSpeedMmS * 1e-3, hold = cs.sleepHoldMs * 1e-3;
      bool ok = true, bandKick = false, partBlock = false;
      double vmax = 0.0, spinMax = 0.0, minHold = 1e9;
      for (const int i : isl.members)
      {
         const Body &b = cy.bodies[i];
         HitBall *const pb = b.pb;
         Body s;
         LoadBody(s, pb);
         LoadState(s);
         const double v = s.v.Length(), spin = s.w.Length() * s.bs.R;
         vmax = std::max(vmax, v);
         spinMax = std::max(spinMax, spin);
         pb->m_sleepHoldS = std::max(v, spin) < vs ? pb->m_sleepHoldS + PHYSICS_STEPTIME_S : 0.0;
         minHold = std::min(minHold, pb->m_sleepHoldS);
         if (b.band != nullptr || pb->m_bandOwner != nullptr || pb->m_d.m_lockedInKicker)
            bandKick = true;
         if (pb->m_contact != nullptr)
            for (const SurfaceContact &it : pb->m_contact->items)
               if (const ContactMovingPart *const mp = it.r->MovingPart(); mp != nullptr && !PartQuiet(mp))
                  partBlock = true;
      }
      ok = minHold >= hold - 1e-9 && !bandKick && !partBlock;
      const double vtop = std::max(vmax, spinMax);
      for (const int i : isl.members)
      {
         HitBall *const pb = cy.bodies[i].pb;
         pb->m_sleepCreepS = vtop >= vs && vtop < kSleepCreepMS ? pb->m_sleepCreepS + PHYSICS_STEPTIME_S : 0.0;
         if (pb->m_sleepCreepS > kSleepCreepS && !g_sleepCreepWarned)
         {
            g_sleepCreepWarned = true;
            std::ostringstream os;
            os << "Pile sleep: a group of " << isl.members.size() << " balls has crept at " << vtop * 1e3 << " mm/s (between " << cs.sleepSpeedMmS << " and 1 mm/s) for over 5 s at t=" << g_pplayer->m_time_msec
               << " ms: it never comes to rest (requirement 4.1)";
            SleepNotice(os.str());
         }
      }
      const auto whyNot = [&](const char *stage, const double deepM)
      {
         if (vtop >= vs)
            return;
         static std::vector<std::pair<const HitBall *, uint32_t>> last;   // per group (its first ball), when it was last logged
         const HitBall *const key = cy.bodies[isl.members[0]].pb;
         const uint32_t now = g_pplayer->m_time_msec;
         auto at = std::find_if(last.begin(), last.end(), [key](const std::pair<const HitBall *, uint32_t> &p) { return p.first == key; });
         if (at != last.end() && now - at->second < 1000)
            return;
         if (at != last.end())
            at->second = now;
         else if (last.size() < 64)
            last.emplace_back(key, now);
         std::ostringstream os;
         os << "Pile sleep: a group of " << isl.members.size() << " balls under the speed line is not asleep at t=" << now << " ms, held by: " << stage << ". Worst speed " << vmax * 1e3
            << " mm/s, worst spin surface speed " << spinMax * 1e3 << " mm/s (line " << cs.sleepSpeedMmS << "), shortest hold " << minHold * 1e3 << " ms (need " << cs.sleepHoldMs << "), band or kicker "
            << (int)bandKick << ", coil currents steady " << (int)Magnets::CurrentSteady() << ", deepest contact " << deepM * 1e6 << " um (-1 = not checked)";
         std::vector<const ContactMovingPart *> seen;
         for (const int i : isl.members)
            if (const BallContact *const e = cy.bodies[i].pb->m_contact)
               for (const SurfaceContact &it : e->items)
                  if (const ContactMovingPart *const mp = it.r->MovingPart(); mp != nullptr && !PartQuiet(mp) && std::find(seen.begin(), seen.end(), mp) == seen.end())
                  {
                     seen.push_back(mp);
                     double mv = 0.0, sp = 0.0;
                     os << "; part not at rest at (" << mp->Origin().x << ", " << mp->Origin().y << ") m: ";
                     if (mp->RestNumbers(mv, sp))
                        os << "moved " << mv * 1e9 << " nm over its last cycle (AtRest limit 1 nm), surface speed " << sp * 1e6 << " um/s (limit 1 um/s), ";
                     os << mp->StateNote();
                  }
         PLOGI << os.str();
      };
      const bool steady = Magnets::CurrentSteady();
      if (!ok || !steady)
      {
         whyNot(minHold < hold - 1e-9 ? "the hold" : bandKick ? "a band or kicker" : partBlock ? "a moving part not at rest" : "a coil current not steady", -1.0);
         return;
      }
      // no contact deeper than requirement 4.2 allows: a pile held that deep is a fault to see, not to hide
      double deep = 0.0;
      for (size_t a = 0; a < isl.members.size(); a++)
      {
         Body s;
         LoadBody(s, cy.bodies[isl.members[a]].pb);
         LoadState(s);
         for (size_t b = a + 1; b < isl.members.size(); b++)
         {
            Body t;
            LoadBody(t, cy.bodies[isl.members[b]].pb);
            LoadState(t);
            deep = std::max(deep, s.bs.R + t.bs.R - (s.c - t.c).Length());
         }
         if (s.pb->m_contact != nullptr)
            for (const SurfaceContact &it : s.pb->m_contact->items)
            {
               const TubeProbe pr = it.r->Probe(s.bs, it, s.c);
               if (pr.overlap)
                  deep = std::max(deep, pr.delta);
            }
      }
      if (deep > kSleepDepthM)
      {
         whyNot("the depth guard", deep);
         if (!g_sleepDepthWarned)
         {
            g_sleepDepthWarned = true;
            std::ostringstream os;
            os << "Pile sleep: a group of " << isl.members.size() << " balls at rest holds a contact " << deep * 1e3 << " mm deep (limit 0.05 mm) at t=" << g_pplayer->m_time_msec << " ms: it is not put to sleep";
            SleepNotice(os.str());
         }
         return;
      }
      SleepGroup g;
      g.id = g_sleepNextId++;
      for (const int i : isl.members)
      {
         HitBall *const pb = cy.bodies[i].pb;
         g.balls.push_back(pb);
         if (pb->m_contact != nullptr)
            for (const SurfaceContact &it : pb->m_contact->items)
               if (ContactMovingPart *const mp = it.r->MovingPart(); mp != nullptr && std::find(g.parts.begin(), g.parts.end(), mp) == g.parts.end())
                  g.parts.push_back(mp);
         pb->m_d.m_vel = Vertex3Ds(0.f, 0.f, 0.f);
         pb->m_angularmomentum = Vertex3Ds(0.f, 0.f, 0.f);
         pb->m_sleepPos = pb->m_d.m_pos;
         pb->m_sleepGroup = g.id;
         pb->m_sleepOwe = false;
         pb->m_sleepCreepS = 0.0;
      }
      g.partIdle.assign(g.parts.size(), 0);
      if (SleepLogOk())
         PLOGI << "Pile sleep: group " << g.id << " (" << g.balls.size() << " balls, " << g.parts.size() << " moving parts) sleeps at t=" << g_pplayer->m_time_msec << " ms, max speed " << vmax * 1e3
               << " mm/s, max spin surface speed " << spinMax * 1e3 << " mm/s, deepest contact " << deep * 1e6 << " um";
      g_sleep.push_back(std::move(g));
   }

   static bool s_append;   // StepBalls: do not flush or clear the queued events: the steps of one replayed sub-interval chain
   static bool PreventOn()   // point 1 of the overlap standard; VPX_PREVENT=0 switches it off (test only)
   {
      static const bool on = []() { const char *const e = getenv("VPX_PREVENT"); return !(e != nullptr && atoi(e) == 0); }();
      return on;
   }
   static void CheckOverlaps(float hittime);
   static void SeatAtTouch(const std::vector<SurfaceContact> &items, const BallSI &bs, V3 &c, const V3 &v, double &delta0, double &rewound);
   static void SeatKeepsOthers(const std::vector<SurfaceContact> &items, const ContactSurface *self, const BallSI &bs, const V3 &c0, V3 &c);
   static void DispatchPending()
   {
      if (s_pendingCy != nullptr)
      {
         Cycle &cy = *s_pendingCy;
         s_pendingCy = nullptr;
         ApplyReleases(cy);
         DispatchEvents(cy);
      }
      else
         DrainDeferred();
   }
};
ContactSystemImpl::Cycle *ContactSystemImpl::s_pendingCy = nullptr;
bool ContactSystemImpl::s_append = false;

void ContactProfCensus(const ContactSurface *s, bool empty, size_t cand);

void ContactSystemImpl::StepBalls(const float dtime)
{
   if (!s_append)
      DispatchPending();   // (the engine dispatches after every step; this keeps a step that was not followed by one from losing its events)
   g_last = ContactSystem::CycleStats();
   if (g_engaged.empty() && !(ContactSettings::Get().ballContact && g_pplayer->m_vball.size() >= 2))   // (two free balls are solved here too)
      return;
   const double dt = (double)dtime * DEFAULT_STEPTIME_S;
   if (dt <= 0.0)
      return;
   static Cycle cy;   // buffers are reused between cycles
   cy.dt = dt;
   cy.cs = &ContactSettings::Get();
   cy.bodies.clear();
   if (!s_append)   // (the steps of a replayed sub-interval add to the queue)
   {
      cy.releases.clear();
      cy.hitEvents.clear();
      cy.partHits.clear();
   }
   cy.hitFloor = std::min(0.5f, g_pplayer->m_physics->HitSpeedMin());
   cy.dtimeF = dtime;
   if (!s_append)
      cy.bandReleases.clear();
   for (BallPairState &p : g_ballPairs)
      p.used = false;
   cy.rpBall = RubberParams();
   cy.rpBall.eMpa = cy.cs->ballEMpa * 0.5;
   cy.rpBall.nu = cy.cs->ballNu;
   cy.rpBall.alphaSPerM = cy.cs->ballAlpha;
   cy.rpBall.mu = cy.cs->ballMu;
   cy.rpBall.tanDelta = cy.cs->ballTanDelta;
   cy.rpBall.rollingMu = 0.0;
   cy.rpBall.normalLoss = cy.cs->ballTanDelta > 0.0;

   CPROF(StepBalls);
   if (cy.cs->sleep)
      SleepInStep(cy);
   s_cc = TestDriver::Contact();
   static const bool timed = getenv("VPX_CONTACT_TIMING") != nullptr && *getenv("VPX_CONTACT_TIMING") != 0;
   uint64_t tp = timed ? usec() : 0;
   auto lap = [&](int &acc) { if (timed) { const uint64_t t = usec(); acc += (int)(t - tp); tp = t; } };
   GatherBalls(cy);
   GatherFree(cy);
   FindIslands(cy);
   PairFirstSight(cy);
   EnergyCheckStart(cy);
   lap(g_last.gatherUs);
   g_last.bodies = (int)cy.bodies.size();
   g_last.islands = (int)cy.islands.count;
   for (size_t ii = 0; ii < cy.islands.count; ii++)
   {
      const ContactSystem::Islands::Island &isl = cy.islands.list[ii];
      const uint64_t islT0 = s_cc != nullptr ? TestDriver::Ticks() : 0;
      if (s_cc != nullptr)
      {
         std::fill(std::begin(s_islWhy), std::end(s_islWhy), 0);
         s_islMinStep = 1.0;
      }
      if (IslandHasCord(cy, isl))
      {
         CordIsland(cy, isl);
         lap(g_last.substepUs);
         EnergyCheckEnd(cy, isl);
         continue;
      }
      if (IslandHasGroupPart(cy, isl))
      {
         CordIsland(cy, isl, true);
         lap(g_last.substepUs);
         EnergyCheckEnd(cy, isl);
         SleepConsider(cy, isl);
         continue;
      }
      if (IslandOnRing(cy, isl) && CordIsland(cy, isl, true))
      {
         lap(g_last.substepUs);
         EnergyCheckEnd(cy, isl);
         SleepConsider(cy, isl);
         continue;
      }
      RingHandOut(cy, isl);
      SetupIsland(cy, isl);
      lap(g_last.setupUs);
      if (ContactProf::On())
         for (const int i : isl.members)
            if (cy.bodies[i].e != nullptr)
               for (const SurfaceContact &it : cy.bodies[i].e->items)
                  ContactProfCensus(it.r, it.candR >= 0.0 && it.cand.empty(), it.cand.size());
      g_last.substeps += cy.steps;
      for (const int i : isl.members)
         if (cy.bodies[i].e != nullptr)
         {
            g_last.items += (int)cy.bodies[i].e->items.size();
            g_last.probes += (long long)cy.bodies[i].e->items.size() * cy.steps;
         }
      g_last.maxSteps = std::max(g_last.maxSteps, cy.steps);
      for (const int i : isl.members)
      {
         g_last.bandBodies += cy.bodies[i].band != nullptr ? 1 : 0;
         g_last.restBodies += cy.bodies[i].steps < 25 ? 1 : 0;
      }
      const bool tripped = Substeps(cy, isl);
      lap(g_last.substepUs);
      const bool top = s_cc != nullptr && cy.steps > s_cc->topSteps;
      if (top)
         IslandMembers(cy, isl, s_cc->topMembers);   // (before the write-back, which may release a ball's contact)
      if (tripped)
         StepFailed(cy, isl);
      else
      {
         WriteBack(cy, isl);
         EnergyCheckEnd(cy, isl);
         SleepConsider(cy, isl);
      }
      lap(g_last.writeUs);
      if (s_cc != nullptr)
      {
         s_cc->islands++;
         if (top)
         {
            s_cc->topSteps = cy.steps;
            s_cc->topTicks = TestDriver::Ticks() - islT0;
            s_cc->topMinStepUs = s_islMinStep * 1e6;
            std::copy(std::begin(s_islWhy), std::end(s_islWhy), s_cc->topWhy);
         }
      }
   }
   if (cy.cs->islandReach)
      ReachClear();
   s_pendingCy = &cy;   // its events are dispatched by ContactSystem::Dispatch, after the free balls have moved
   lap(g_last.writeUs);
   if (ContactProf::On())
   {
      static int calls = 0;
      static const int period = []() { const char *const e = getenv("VPX_CONTACT_PROFILE"); const int n = e != nullptr ? atoi(e) : 0; return n > 1 ? n : 5000; }();   // VPX_CONTACT_PROFILE=<n>: a dump every n step calls
      if (++calls % period == 0)
         ContactProf::Dump("periodic");
   }
}


// ---- Proximity engagement ---------------------------------------------------------------------------------------------------------------------------------------------
// A free ball that can touch a surface this cycle is taken into the contact before the stock hit search runs: the surface is engaged when the ball is within its release gap, or
// closing on it faster than it can cross the gap in two cycles (reach 2 vn dt + the release gap; a ball released at the gap, receding, is not taken back). The first touch is
// then found by the substeps, and the part's hit event is queued from there (PartHit).
// A cord takes a ball before it touches, and the group steps it from there (GroupStepDesign 15): the band is never found by the stock collision, and the ball does not start inside the cord.
// Every ball, one already held for the playfield included: its contact is handed over to the band.
void ContactSystemImpl::CordNear(const float dtime)
{
   const double dt = (double)dtime * DEFAULT_STEPTIME_S;
   for (auto *const ball : g_pplayer->m_vball)
   {
      HitBall *const pb = &ball->m_hitBall;
      if (pb->m_bandOwner != nullptr || pb->m_d.m_lockedInKicker)
         continue;
      const Vertex3Ds &pos = pb->m_d.m_pos;
      const double speed = std::sqrt((double)pb->m_d.m_vel.Dot(pb->m_d.m_vel)) * VS;
      BallSI bs;
      bs.R = (double)pb->m_d.m_radius * U;
      const V3 c((double)pos.x * U, (double)pos.y * U, (double)pos.z * U);
      const V3 v((double)pb->m_d.m_vel.x * VS, (double)pb->m_d.m_vel.y * VS, (double)pb->m_d.m_vel.z * VS);
      for (ContactHost *const o : ContactHost::Registry())
      {
         HitRubber *const band = o->m_set.band && o->m_enabled ? o->BandHost() : nullptr;
         if (band == nullptr || !band->CordModel())
            continue;
         // a plate on its way out carries the cord this far past its rest line (as the release rule)
         const double plateReach = band->CordOn() && !band->Cord()->PlateQuiet() ? 0.03 : 0.0;
         // a real section thicker than the table's: the cord's surface stands out from the table's tube by twice the difference (its centre moves out by it, and its radius grows by it); the probe below reads the table's tube
         const double fixedRc = band->CordOn() ? band->Cord()->fixed->rc : 0.0;
         const float reachU = pb->m_d.m_radius + (float)((2.0 * speed * dt + plateReach + (fixedRc > 0.0 ? 2.0 * fixedRc : 0.0)) / U) + 4.0f;
         const FRect3D &bb = o->m_hitBBox;
         if (pos.x < bb.left - reachU || pos.x > bb.right + reachU || pos.y < bb.top - reachU || pos.y > bb.bottom + reachU || pos.z < bb.zlow - reachU || pos.z > bb.zhigh + reachU)
            continue;
         o->ReloadIfChanged();
         SurfaceContact it;
         it.r = o;
         const TubeProbe pr = o->Probe(bs, it, c);
         if (pr.dist <= 0.0)
            continue;
         const double gap = pr.dist - (bs.R + pr.rho) - (fixedRc > 0.0 ? 2.0 * (fixedRc - pr.rho) : 0.0), rel = o->ReleaseGapM(), closing = -v.Dot(pr.normal);
         if ((gap <= rel + plateReach || (closing > 0.0 && gap < closing * 2.0 * dt + rel)) && (pb->m_sleepGroup == 0 || (SleepWakeBall(pb, "a band can take it", true), true)) && band->BandEngage(pb))
            break;
      }
   }
}

void ContactSystemImpl::EngageNear(const float dtime)
{
   CordNear(dtime);
   CPROF(EngageNear);
   ContactProf::Cycle();
   static const bool on = []()
   {
      const char *const e = getenv("VPX_PROXIMITY");
      if (e != nullptr)
         return atoi(e) != 0;
      const std::filesystem::path ini = ContactHost::IniPathFor(g_pplayer->m_ptable->m_filename);
      return PhysicsIni::Active(ini).proximityEngage;
   }();
   if (!on)
      return;
   // The faces near each free ball, per surface, kept between cycles (a Verlet list: Refresh only rescans when the ball has moved out of the list): scanning a mesh playfield
   // afresh for every free ball every cycle would cost milliseconds.
   struct NearEntry
   {
      const HitBall *ball;
      const ContactHost *host;
      SurfaceContact it;
      uint32_t seenMs;
   };
   static std::vector<NearEntry> cache;
   const uint32_t nowMs = g_pplayer->m_time_msec;
   for (size_t i = 0; i < cache.size();)
      if (nowMs - cache[i].seenMs > 200)
      {
         cache[i] = std::move(cache.back());
         cache.pop_back();
      }
      else
         i++;
   const double dt = (double)dtime * DEFAULT_STEPTIME_S;
   for (auto *const ball : g_pplayer->m_vball)
   {
      HitBall *const pb = &ball->m_hitBall;
      if (pb->m_contact != nullptr || pb->m_bandOwner != nullptr || pb->m_d.m_lockedInKicker || pb->m_sleepGroup != 0)
         continue;
      const Vertex3Ds &pos = pb->m_d.m_pos;
      const double speed = std::sqrt((double)pb->m_d.m_vel.Dot(pb->m_d.m_vel)) * VS;
      const double travelM = 2.0 * speed * dt;
      const float reachU = pb->m_d.m_radius + (float)(travelM / U) + 20.0f;
      // Verlet: the last scan found every surface farther than m_nearSafeU from where the ball was then (box distance less the radius and a margin); until the ball has
      // travelled that far, less what it can cross in two cycles, there is nothing to find
      if (pb->m_nearSafeU >= 0.f && nowMs - pb->m_nearMs < 10 && std::sqrt((double)(pos - pb->m_nearPos).LengthSquared()) + travelM / U < (double)pb->m_nearSafeU)
         continue;
      float minBoxGapU = 1e9f;
      BallSI bs;
      bs.R = (double)pb->m_d.m_radius * U;
      const V3 c((double)pos.x * U, (double)pos.y * U, (double)pos.z * U);
      const V3 v((double)pb->m_d.m_vel.x * VS, (double)pb->m_d.m_vel.y * VS, (double)pb->m_d.m_vel.z * VS);
      ContactHost *best = nullptr;
      double bestGap = 1e9, pfGap = 1e9, pfClosing = 0.0;
      for (ContactHost *const o : ContactHost::Registry())
      {
         if (!o->m_enabled || !o->m_set.contact || o->m_set.band)
            continue;
         const double vb = PreventOn() ? o->SurfaceSpeedBoundMS(dt) : 0.0;   // the surface's own closing-speed bound (a moving part; 0 for a static surface)
         const float reachHost = reachU + (float)(2.0 * vb * dt / U);
         const FRect3D &bb = o->m_hitBBox;
         {
            const float dx = std::max({ bb.left - pos.x, 0.f, pos.x - bb.right }), dy = std::max({ bb.top - pos.y, 0.f, pos.y - bb.bottom }), dz = std::max({ bb.zlow - pos.z, 0.f, pos.z - bb.zhigh });
            minBoxGapU = std::min(minBoxGapU, std::sqrt(dx * dx + dy * dy + dz * dz) - pb->m_d.m_radius);
         }
         if (pos.x < bb.left - reachHost || pos.x > bb.right + reachHost || pos.y < bb.top - reachHost || pos.y > bb.bottom + reachHost || pos.z < bb.zlow - reachHost || pos.z > bb.zhigh + reachHost)
            continue;
         if (o->m_surfaceKind == "playfield" && InPlayfieldHole(pos))
            continue;   // the kicker's
         o->ReloadIfChanged();
         auto ce = std::find_if(cache.begin(), cache.end(), [&](const NearEntry &n) { return n.ball == pb && n.host == o; });
         if (ce == cache.end())
         {
            cache.push_back({ pb, o, SurfaceContact(), nowMs });
            ce = cache.end() - 1;
            ce->it.r = o;
         }
         ce->seenMs = nowMs;
         SurfaceContact &it = ce->it;
         o->Refresh(it, bs.R, travelM + 2.0 * vb * dt + 0.005, c);
         const TubeProbe pr = o->Probe(bs, it, c);
         if (pr.dist <= 0.0)
            continue;
         const double gap = pr.dist - (bs.R + pr.rho);
         const double rel = o->ReleaseGapM();
         const double closing = -v.Dot(pr.normal) + vb;   // relative to the surface: its bound counts as closing speed
         if (o->m_surfaceKind == "playfield")
         {
            pfGap = gap;
            pfClosing = closing;
         }
         if ((gap <= rel || (closing > 0.0 && gap < closing * 2.0 * dt + rel)) && gap < bestGap)
         {
            best = o;
            bestGap = gap;
         }
      }
      if (best != nullptr)
         EngageImpl(best, pb, "near");
      else
      {
         if (ContactProf::On() && pb->m_d.m_pos.y > 1500.f)
         {
            const double mm = pfGap * 1e3;
            ContactProf::Note(std::string("free ball not engaged (y>1500): playfield gap ") + (pfGap > 1e8 ? "none(box/hole/cooldown)" : mm <= 0.0 ? "<=0mm" : mm < 0.5 ? "<0.5mm" : mm < 1.0 ? "<1mm" : mm < 3.0 ? "<3mm" : mm < 10.0 ? "<10mm" : ">=10mm") + (pfClosing > 0.0 ? " closing" : " receding"));
         }
         pb->m_nearPos = pos;
         pb->m_nearMs = nowMs;
         pb->m_nearSafeU = minBoxGapU - 6.0f;   // (a surface within its release gap, about 4 U, qualifies; the margin keeps the box test conservative)
      }
   }
}

// ---- Engagement: which balls the contact system holds ----------------------------------------------------------------------------------------------------------------
// Go back to the instant the ball first touched: along its own velocity, by depth / closing speed, at most one physics step; if it is not closing, to touching distance along the normal.
void ContactSystemImpl::SeatAtTouch(const std::vector<SurfaceContact> &items, const BallSI &bs, V3 &c, const V3 &v, double &delta0, double &rewound)
{
   constexpr double kMaxRewindS = PHYSICS_STEPTIME_S;
   const V3 c0 = c;
   for (int pass = 0; pass < 3; pass++)
      for (const SurfaceContact &it : items)
      {
         const TubeProbe p0 = it.r->Probe(bs, it, c);
         if (!p0.overlap || p0.delta <= 1e-6)
            continue;
         delta0 = std::max(delta0, p0.delta);
         if (ContactTraceOn())
            PLOGI << "ContactSeat " << it.r->Name() << " pass " << pass << " overlap_um " << p0.delta * 1e6 << " normal " << p0.normal.x << ' ' << p0.normal.y << ' ' << p0.normal.z << " ball " << c.x / U << ' ' << c.y / U << ' ' << c.z / U;
         const double closing = -v.Dot(p0.normal);
         if (closing > 0.01)
         {
            const double t = std::min(p0.delta / closing, kMaxRewindS);
            c = c - v * t;
            rewound = std::max(rewound, t);
            const TubeProbe p1 = it.r->Probe(bs, it, c);
            if (p1.overlap && p1.delta > 1e-6)
               c = c + p1.normal * p1.delta;
         }
         else
            c = c + p0.normal * p0.delta;
         SeatKeepsOthers(items, it.r, bs, c0, c);
      }
}

// A seat is a position move: it does no work, so it may store nothing (no stored energy without the work that stored it). It may not press the ball deeper into any other surface of its contact than
// it was at c0: that depth is taken back along that surface's normal, and the next pass of the seat goes on along the surface.
void ContactSystemImpl::SeatKeepsOthers(const std::vector<SurfaceContact> &items, const ContactSurface *const self, const BallSI &bs, const V3 &c0, V3 &c)
{
   for (const SurfaceContact &it : items)
   {
      if (it.r == self)
         continue;
      const TubeProbe p1 = it.r->Probe(bs, it, c);
      if (!p1.overlap)
         continue;
      const TubeProbe p0 = it.r->Probe(bs, it, c0);
      const double more = p1.delta - (p0.overlap ? p0.delta : 0.0);
      if (more > 1e-6)
         c = c + p1.normal * more;
   }
}

void ContactSystem::ShoveKeepsOthers(const HitBall *const pb, Vertex3Ds &delta)
{
   if (pb->m_contact == nullptr)
      return;
   BallSI bs;
   bs.R = (double)pb->m_d.m_radius * U;
   const Vertex3Ds &p = pb->m_d.m_pos;
   const V3 c0((double)p.x * U, (double)p.y * U, (double)p.z * U);
   V3 c = c0 + V3((double)delta.x * U, (double)delta.y * U, (double)delta.z * U);
   for (int pass = 0; pass < 3; pass++)
      ContactSystemImpl::SeatKeepsOthers(pb->m_contact->items, nullptr, bs, c0, c);
   delta = Vertex3Ds((float)((c.x - c0.x) / U), (float)((c.y - c0.y) / U), (float)((c.z - c0.z) / U));
}

void ContactSystem::SeatKeepsOthers(const HitBall *const pb, const ContactSurface *const self, const RubberContact::BallSI &bs, const RubberContact::V3 &c0, RubberContact::V3 &c)
{
   if (pb->m_contact != nullptr)
      ContactSystemImpl::SeatKeepsOthers(pb->m_contact->items, self, bs, c0, c);
}

// A ball enters when a surface's stock collider finds it (the stock impulse is skipped) and leaves when it is clear of every surface and receding.
bool ContactSystemImpl::EngageImpl(ContactHost *const host, HitBall *const pb, const char *how)
{
   CPROF(EngageImpl);
   using namespace RubberContact;
   if ((host != nullptr && !host->m_set.contact) || pb->m_bandOwner != nullptr || pb->m_contact != nullptr || pb->m_d.m_lockedInKicker)
      return false;   // (a ball held in a kicker belongs to the kicker)
   if (host != nullptr && host->m_surfaceKind == "playfield" && InPlayfieldHole(pb->m_d.m_pos))
      return false;   // the kicker's
   BallContact rb;
   rb.ball = pb;
   BallSI bs;
   bs.R = (double)pb->m_d.m_radius * U;
   const Vertex3Ds posStart = pb->m_d.m_pos;
   V3 c((double)posStart.x * U, (double)posStart.y * U, (double)posStart.z * U);
   const V3 v((double)pb->m_d.m_vel.x * VS, (double)pb->m_d.m_vel.y * VS, (double)pb->m_d.m_vel.z * VS);
   {
      const float reach = pb->m_d.m_radius + 20.0f;   // VPX units: the ball radius plus what it can travel during a hit
      for (ContactHost *const o : ContactHost::Registry())
      {
         if (o != host)
         {
            if (!o->m_enabled || !o->m_set.contact)
               continue;
            if (o->IsPlayfield() && InPlayfieldHole(posStart))
               continue;   // no floor over a kicker's hole
            if (o->m_set.band)
               continue;   // the band rubbers are the band model's: as a rigid surface contact they would stop the ball before the band could deflect (as Regather)
            const FRect3D &bb = o->m_hitBBox;
            if (posStart.x < bb.left - reach || posStart.x > bb.right + reach || posStart.y < bb.top - reach || posStart.y > bb.bottom + reach
                || posStart.z < bb.zlow - reach || posStart.z > bb.zhigh + reach)
               continue;
            o->ReloadIfChanged();
         }
         SurfaceContact it;
         it.r = o;
         it.rp = o->m_set.ToParams();
         it.hLayer = o->LayerM();
         it.st.Reset();
         it.touching = o == host && std::strcmp(how, "near") != 0;   // the stock collider that found the ball fires its own hit event; a surface found by proximity gets its event from the first touch
         rb.items.push_back(it);
      }
   }
   // The ball starts the engagement touching or a little inside the rubber (the hit was found within a step). Go back to the instant it first touched: along its
   // own velocity, by depth / closing speed, at most one physics step; if it is not closing (pushed in by something else), seat it at touching distance along the normal.
   // Nothing is stored as compression.
   // Found inside a part with no earlier clean state (docs/UnifiedContactDesign.txt 11, point 6): CheckOverlaps found no clean start of the sub-interval, or the engagement by proximity finds the ball already inside a
   // moving part (point 1 engages a ball that has a past before it touches, so a ball inside at that point was created or moved there)
   bool noPast = std::strcmp(how, "no past") == 0 && !(host != nullptr && host->MovingPart() != nullptr && host->MovingPart()->GroupPart() != nullptr);   // (a part on the group step: its zero, below)
   if (!noPast && host != nullptr && std::strcmp(how, "near") == 0)
      for (SurfaceContact &it : rb.items)
         if (it.r->MovingPart() != nullptr && it.r->MovingPart()->GroupPart() == nullptr)
         {
            it.r->Refresh(it, bs.R, 0.03, c);
            const TubeProbe p0 = it.r->Probe(bs, it, c);
            if (p0.overlap && p0.delta > 1e-6)
               noPast = true;
         }
   double noPastMoved = 0.0, noPastLeft = 0.0;
   int noPastPasses = 0;
   if (noPast && host != nullptr)
   {
      // Minimum translation: along the normal of each overlapped surface (the playfield is one of them) to touching distance, repeated up to 8 passes. Velocities unchanged. A ball that is squeezed and does not fit keeps
      // the overlap that is left, and the force law carries it.
      const V3 c0 = c;
      for (SurfaceContact &it : rb.items)
         it.r->Refresh(it, bs.R, 0.03, c);
      for (; noPastPasses < 8; noPastPasses++)
      {
         bool moved = false;
         for (const SurfaceContact &it : rb.items)
         {
            const TubeProbe p0 = it.r->Probe(bs, it, c);
            if (!p0.overlap || p0.delta <= 1e-6)
               continue;
            if (noPastPasses == 0)
               rb.delta0 = std::max(rb.delta0, p0.delta);
            c = c + p0.normal * p0.delta;
            moved = true;
         }
         if (!moved)
            break;
      }
      noPastMoved = (c - c0).Length();
      for (const SurfaceContact &it : rb.items)
      {
         const TubeProbe p0 = it.r->Probe(bs, it, c);
         if (p0.overlap)
            noPastLeft = std::max(noPastLeft, p0.delta);
      }
   }
   else if (host != nullptr)
   {
   for (SurfaceContact &it : rb.items)
      it.r->Refresh(it, bs.R, 0.03, c);   // a surface: the faces within 30 mm of the ball
   // Found inside a surface going to the group: zero replaces the seat, no position snap (GroupStepDesign 14.4, 14.5); not with a moving part in the set, which stays on the old loop
   bool zeroNotSeat = false;
   for (const SurfaceContact &it : rb.items)
      if (it.r->GroupModel() || (it.r->MovingPart() != nullptr && it.r->MovingPart()->GroupPart() != nullptr))
      {
         const TubeProbe p0 = it.r->Probe(bs, it, c);
         zeroNotSeat = zeroNotSeat || (p0.overlap && p0.delta > 1e-6);
      }
   if (std::any_of(rb.items.begin(), rb.items.end(), [](const SurfaceContact &it) { return it.r->MovingPart() != nullptr && it.r->MovingPart()->GroupPart() == nullptr; }))
      zeroNotSeat = false;
   rb.seatPending = zeroNotSeat;
   if (!zeroNotSeat)
      SeatAtTouch(rb.items, bs, c, v, rb.delta0, rb.rewound);
   }
   pb->m_d.m_pos = Vertex3Ds((float)(c.x / U), (float)(c.y / U), (float)(c.z / U));
   rb.vin = v;
   const double inertiaU = (double)pb->Inertia();
   const double wScale = 1.0 / (inertiaU * DEFAULT_STEPTIME_S);
   rb.win = V3((double)pb->m_angularmomentum.x * wScale, (double)pb->m_angularmomentum.y * wScale, (double)pb->m_angularmomentum.z * wScale);
   rb.posStart = posStart;
   if (ContactTraceOn())
      PLOGI << "ContactEngage " << (host != nullptr ? host->m_name : std::string("band")) << " t=" << g_pplayer->m_time_msec << " how " << how << " speed " << v.Length() << " vn_in " << rb.vin.Length() << " pos " << posStart.x << ' ' << posStart.y << ' ' << posStart.z << " delta0_um " << rb.delta0 * 1e6
            << " rewind_us " << rb.rewound * 1e6 << [&]() { std::string s; for (const SurfaceContact &it : rb.items) if (ContactMovingPart *const mp = it.r->MovingPart()) s += " [" + it.r->Name() + " " + mp->StateNote() + "]"; return s; }() << " flip_clear " << g_pplayer->m_physics->FlipperClearance(pb->m_d, kFlipperHorizonVpt) << g_pplayer->m_physics->FlipperDebug(pb->m_d, kFlipperHorizonVpt) << " rubbers " << rb.items.size() << " pos " << pb->m_d.m_pos.x << ' ' << pb->m_d.m_pos.y << ' ' << pb->m_d.m_pos.z;
   if (noPast)
   {
      static uint32_t lastMs = 0;
      if (g_pplayer->m_time_msec - lastMs > 500)
      {
         lastMs = g_pplayer->m_time_msec;
         PLOGW << "Ball with no past found inside " << (host != nullptr ? host->m_name : std::string("?")) << " (" << rb.delta0 * 1e6 << " um deep at first sight): seated by minimum translation, moved " << noPastMoved * 1e6 << " um in " << noPastPasses << " passes, " << noPastLeft * 1e6 << " um overlap left; ball at " << posStart.x << ' ' << posStart.y << ' ' << posStart.z;
      }
      static std::vector<std::string> shown;   // one notification per part per table start
      const std::string nm = host != nullptr ? host->m_name : std::string();
      if (g_pplayer->m_liveUI != nullptr && std::find(shown.begin(), shown.end(), nm) == shown.end())
      {
         shown.push_back(nm);
         g_pplayer->m_liveUI->PushNotification("A ball was found inside " + nm + " with no earlier position (created or moved there by a script, or Collidable turned on): moved out by the shortest way (see vpinball.log)", 8000);
      }
   }
   rb.host = host;
   g_engaged.push_back(rb);
   pb->m_contact = &g_engaged.back();
   return true;
}

void ContactSystemImpl::Release(HitBall *const ball, const char *why)
{
   if (ball->m_contact == nullptr)
      return;
   const auto itE = std::find_if(g_engaged.begin(), g_engaged.end(), [ball](const BallContact &r) { return &r == ball->m_contact; });
   if (itE == g_engaged.end())
      return;
   if (itE->band != nullptr && ball->m_bandOwner == itE->band)
      itE->band->BandRelease(ball, why);   // clears itE->band and the ball's band owner
   itE->band = nullptr;
   const BallContact e = *itE;
   ContactHost *const host = e.host;
   if (getenv("VPX_RUBBERS_LOG"))
   {
      using namespace RubberContact;
      const double inertiaU = (double)ball->Inertia();
      const double wScale = 1.0 / (inertiaU * DEFAULT_STEPTIME_S);
      const V3 v((double)ball->m_d.m_vel.x * VS, (double)ball->m_d.m_vel.y * VS, (double)ball->m_d.m_vel.z * VS);
      const V3 w((double)ball->m_angularmomentum.x * wScale, (double)ball->m_angularmomentum.y * wScale, (double)ball->m_angularmomentum.z * wScale);
      PLOGI << "RubberSolve " << (host != nullptr ? host->m_name : std::string("band")) << " t=" << g_pplayer->m_time_msec << " steps " << e.steps << " ms " << e.dwellS * 1e3 << " peakFn " << e.peakFn << " why " << why
            << " vin " << e.vin.x << ' ' << e.vin.y << ' ' << e.vin.z << " vout " << v.x << ' ' << v.y << ' ' << v.z
            << " win " << e.win.x << ' ' << e.win.y << ' ' << e.win.z << " wout " << w.x << ' ' << w.y << ' ' << w.z
            << " delta0 " << e.delta0 * 1e3 << " rewind_ms " << e.rewound * 1e3 << " rubbers " << e.items.size()
            << " pos " << e.posStart.x << ' ' << e.posStart.y << ' ' << e.posStart.z << " (m/s, rad/s, VPX units)";
   }
   if (ContactTraceOn())
      PLOGI << "ContactRelease " << (host != nullptr ? host->m_name : std::string("band")) << " t=" << g_pplayer->m_time_msec << " why " << why << " dwell_ms " << e.dwellS * 1e3 << " peakFn " << e.peakFn
            << " speed " << std::sqrt((double)ball->m_d.m_vel.Dot(ball->m_d.m_vel)) * 0.53975e-3 / DEFAULT_STEPTIME_S << " pos " << ball->m_d.m_pos.x << ' ' << ball->m_d.m_pos.y << ' ' << ball->m_d.m_pos.z;
   if (ContactProf::On())
      ContactProf::Note(std::string("release ") + why);
   ball->m_contact = nullptr;
   ball->m_nearSafeU = -1.f;   // scan again at the next cycle
   g_engaged.erase(itE);
}

void ContactSystem::StepCordAlone(HitRubber *const band, const float dtime)
{
   ContactSystemImpl::CordAlone(band, dtime);
}

bool ContactSystem::EngageBand(HitRubber *const band, HitBall *const pb)
{
   BallContact *e = pb->m_contact;
   if (e != nullptr && e->band != nullptr && e->band != band)
      return false;   // another band has the ball
   if (e == nullptr)
   {
      if (!ContactSystemImpl::EngageImpl(nullptr, pb, "band"))
         return false;
      e = pb->m_contact;
   }
   e->band = band;
   return true;
}

void ContactSystemImpl::Regather(BallContact &e, const Vertex3Ds &pos, const float radiusU)
{
   const float reach = radiusU + 20.0f;   // as at engagement
   // A Verlet list: the surfaces within reach + 10 U are kept and collected again only after the ball has moved 5 U (or 10 ms have passed: a part can become collidable)
   constexpr float kMarginU = 10.0f, kMoveU = 5.0f;
   const uint32_t nowMs = g_pplayer->m_time_msec;
   const bool hole = PlayfieldHolesExist() && InPlayfieldHole(pos);
   if (e.regathered && hole == e.inHole && nowMs - e.regatherMs < 10 && (pos - e.regatherPos).LengthSquared() < kMoveU * kMoveU)
      return;
   e.inHole = hole;
   e.regathered = true;
   e.regatherPos = pos;
   e.regatherMs = nowMs;
   std::vector<SurfaceContact> next;
   for (ContactHost *const o : ContactHost::Registry())
   {
      if (!o->m_enabled || !o->m_set.contact)
         continue;
      if (hole && o->IsPlayfield())
         continue;   // no floor over a kicker's hole
      const auto old = std::find_if(e.items.begin(), e.items.end(), [o](const SurfaceContact &it) { return it.r == o; });
      if (old == e.items.end() && o->m_set.band)
         continue;   // a band (and a sling) is the band model's: as a rigid surface contact it would stop the ball before the band could deflect, and hide the band from its own collision
      const FRect3D &bb = o->m_hitBBox;
      const float rm = reach + kMarginU;
      if (pos.x < bb.left - rm || pos.x > bb.right + rm || pos.y < bb.top - rm || pos.y > bb.bottom + rm || pos.z < bb.zlow - rm || pos.z > bb.zhigh + rm)
         continue;
      if (old != e.items.end())
      {
         next.push_back(std::move(*old));
         continue;
      }
      o->ReloadIfChanged();
      SurfaceContact it;
      it.r = o;
      it.rp = o->m_set.ToParams();
      it.hLayer = o->LayerM();
      it.st.Reset();
      next.push_back(std::move(it));
   }
   if (!next.empty())   // never drop the whole set: the release test below decides when the ball is free
      e.items = std::move(next);
}

void ContactSystem::EngageNear(const float dtime)
{
   ContactSystemImpl::EngageNear(dtime);
}

// A surface's collider, or the slow pass, found the ball: the ball is on this surface's law from now on, whatever holds it. False only for what no surface law covers yet: the playfield
// over a kicker's hole (the kicker, not converted).
bool ContactSystem::Engage(ContactHost *const host, HitBall *const pb, const char *how)
{
   if (pb->m_d.m_lockedInKicker)
      return true;   // the kicker holds it: no surface acts on it, nothing to solve
   if (pb->m_contact == nullptr)
      return ContactSystemImpl::EngageImpl(host, pb, how);
   if (host->IsPlayfield() && InPlayfieldHole(pb->m_d.m_pos))
      return false;
   // Already held, and this surface is not in its set (made collidable, or the ball was moved, since the last gather): it joins the set, as Regather would add it
   std::vector<SurfaceContact> &items = pb->m_contact->items;
   if (std::none_of(items.begin(), items.end(), [host](const SurfaceContact &it) { return it.r == host; }))
   {
      host->ReloadIfChanged();
      SurfaceContact it;
      it.r = host;
      it.rp = host->ContactSettings().ToParams();
      it.hLayer = host->LayerM();
      it.st.Reset();
      it.touching = true;   // the collider that found the ball fires its own hit event
      items.push_back(std::move(it));
   }
   return true;
}

void ContactSystem::Release(HitBall *const ball, const char *why)
{
   ContactSystemImpl::Release(ball, why);
}

void ContactSystem::ReleaseHost(ContactHost *const host)
{
   // the first host to go is the table's end: the energy steps with a moving part, per part
   if (!g_partEnergyTold)
   {
      g_partEnergyTold = true;
      std::string totals;
      for (const auto &p : g_partEnergy)
         totals += " " + p.first + " " + std::to_string(p.second);
      if (totals.empty())
         PLOGI << "Group step: energy steps with a ball, table end: none";
      else
         PLOGW << "Group step: energy steps with a ball, table end, per part in the group:" << totals;
   }
   // the host is going away: its engagements end without a log line (a ball can outlive the physics engine: it must not keep a pointer to the entry)
   for (BallContact &e : g_engaged)
      if (e.band == host)
         e.band = nullptr;
   for (auto it = g_engaged.begin(); it != g_engaged.end();)
      if (it->host == host)
      {
         if (it->ball != nullptr && it->ball->m_contact == &*it)
            it->ball->m_contact = nullptr;
         it = g_engaged.erase(it);
      }
      else
         ++it;
}

bool ContactSystem::Holds(const HitBall *ball)
{
   return ball->m_contact != nullptr;
}

bool ContactSystem::Owns(const BallS &ball, const ContactSurface *surface)
{
   if (ball.m_owner != nullptr)   // the ball knows its contact: only its own surfaces are looked at
   {
      if (const BallContact *const e = ball.m_owner->m_contact)
         for (const SurfaceContact &it : e->items)
            if (it.r == surface)
               return true;
      return false;
   }
   for (const BallContact &e : g_engaged)
      if (&e.ball->m_d == &ball)
         for (const SurfaceContact &it : e.items)
            if (it.r == surface)
               return true;
   return false;
}

size_t ContactSystem::Count()
{
   return g_engaged.size();
}
// ---- Profile output (VPX_CONTACT_PROFILE) ----
namespace
{
   struct Census
   {
      uint64_t itemCycles = 0, emptyCycles = 0, candSum = 0;
      size_t faces = 0;
      std::string name;   // taken when first seen: the surface can be gone by the time the profile is written
   };
   std::unordered_map<const ContactSurface *, Census> g_census;
   std::map<std::string, uint64_t> g_notes;
   std::map<int, uint64_t> g_splits;      // cause -> count
   std::map<std::string, uint64_t> g_splitWhat;   // class of the object that ended the sub-step -> count
   uint64_t g_cyclesSeen = 0, g_stepCalls = 0, g_staticSplits = 0, g_bodyCycles = 0;
   uint64_t g_t0Us = 0, g_tsc0 = 0;
}

void ContactProfCensus(const ContactSurface *s, bool empty, size_t cand)
{
   Census &c = g_census[s];
   if (c.itemCycles == 0)
      c.name = s->Name();
   c.candSum += cand;
   c.itemCycles++;
   c.emptyCycles += empty ? 1 : 0;
   g_bodyCycles++;
}

void ContactProf::Note(const std::string &what)
{
   g_notes[what]++;
}

void ContactProf::Cycle()
{
   if (g_t0Us == 0)
   {
      g_t0Us = usec();
      g_tsc0 = __rdtsc();
   }
   g_cyclesSeen++;
}

void ContactProf::Split(int cause, bool staticTime, const char *what)
{
   g_splits[cause]++;
   g_splitWhat[std::string(what) + (staticTime ? " [static clamp]" : "")]++;
   g_staticSplits += staticTime ? 1 : 0;
}

namespace
{
   struct ProbeCost { std::string name; uint64_t ticks = 0, calls = 0; };
   std::unordered_map<const void *, ProbeCost> &ProbeCosts()
   {
      static std::unordered_map<const void *, ProbeCost> m;
      return m;
   }
}
bool ContactProf::ProbeTime(const void *surface, uint64_t ticks)
{
   auto it = ProbeCosts().find(surface);
   const bool fresh = it == ProbeCosts().end();
   ProbeCost &pc = fresh ? ProbeCosts()[surface] : it->second;
   pc.ticks += ticks;
   pc.calls++;
   return fresh;
}
void ContactProf::ProbeName(const void *surface, const std::string &name) { ProbeCosts()[surface].name = name; }

void ContactProf::Dump(const char *why)
{
   const double us = (double)(usec() - g_t0Us), cyc = (double)(__rdtsc() - g_tsc0);
   const double usPerTick = cyc > 0 ? us / cyc : 0.0;
   const Data &d = D();
   PLOGI << "PROF begin (" << why << ") cycles " << g_cyclesSeen << " wall_us " << (uint64_t)us;
   for (int i = 0; i < Count; i++)
      if (d.n[i] > 0)
         PLOGI << "PROF " << Name(i) << " total_us " << (uint64_t)(d.t[i] * usPerTick) << " calls " << d.n[i] << " us_per_call " << (d.t[i] * usPerTick) / (double)d.n[i] << " us_per_cycle " << (d.t[i] * usPerTick) / std::max<uint64_t>(1, g_cyclesSeen);
   for (const auto &kv : g_notes)
      PLOGI << "PROF note " << kv.first << " " << kv.second;
   std::ostringstream o;
   uint64_t tot = 0;
   for (const auto &kv : g_splits)
      tot += kv.second;
   o << "PROF splits total " << tot << " static_clamp " << g_staticSplits << " by_cause(objtype:count)";
   for (const auto &kv : g_splits)
      o << ' ' << kv.first << ':' << kv.second;
   PLOGI << o.str();
   {
      std::vector<std::pair<uint64_t, std::string>> sw;
      for (const auto &kv : g_splitWhat)
         sw.emplace_back(kv.second, kv.first);
      std::sort(sw.rbegin(), sw.rend());
      for (size_t i = 0; i < sw.size() && i < 25; i++)
         PLOGI << "PROF split by " << sw[i].second << " " << sw[i].first;
   }
   {
      std::vector<const ProbeCost *> v;
      for (const auto &kv : ProbeCosts())
         v.push_back(&kv.second);
      std::sort(v.begin(), v.end(), [](const ProbeCost *a, const ProbeCost *b) { return a->ticks > b->ticks; });
      for (size_t q = 0; q < v.size() && q < 8; q++)
         PLOGI << "PROF probe " << v[q]->name << " total_us " << (uint64_t)(v[q]->ticks * usPerTick) << " calls " << v[q]->calls << " us_per_call " << (v[q]->ticks * usPerTick) / (double)std::max<uint64_t>(1, v[q]->calls);
   }
   std::vector<std::pair<uint64_t, std::string>> top;
   for (const auto &kv : g_census)
      top.emplace_back(kv.second.itemCycles, kv.second.name + " (no face in reach " + std::to_string((int)(100.0 * kv.second.emptyCycles / std::max<uint64_t>(1, kv.second.itemCycles))) + "%, avg faces on list " + std::to_string((double)kv.second.candSum / std::max<uint64_t>(1, kv.second.itemCycles)) + ")");
   std::sort(top.rbegin(), top.rend());
   for (size_t i = 0; i < top.size() && i < 14; i++)
      PLOGI << "PROF surface " << top[i].second << " body-cycles " << top[i].first;
}

void ContactSystem::Dispatch()
{
   ContactSystemImpl::DispatchPending();
}

void ContactSystem::NotePlaced(const HitBall *const ball)
{
   g_placed.push_back(ball);
}

void ContactSystem::SleepCycleStart()
{
   g_frameUsedS = 0.0;   // (a new physics step: called before the balls' velocities)
   ContactSystemImpl::SleepStart();
}

void ContactSystem::SleepWake(const HitBall *const ball, const char *const cause)
{
   ContactSystemImpl::SleepWakeBall(ball, cause, true);
}

void ContactSystem::SleepWakeAtStart(const HitBall *const ball, const char *const cause)
{
   ContactSystemImpl::SleepWakeBall(ball, cause, false);
}

void ContactSystem::Step(float dtime)
{
   ContactSystemImpl::StepBalls(dtime);
   if (dtime > 0.f)
   {
      g_placed.clear();
      g_frameUsedS += (double)dtime * DEFAULT_STEPTIME_S;
   }
}

const ContactSystem::CycleStats &ContactSystem::LastCycle()
{
   return g_last;
}


// ---- The snapshot of the start of a sub-interval (docs/UnifiedContactDesign.txt section 11, point 2) ---------------------------------------------------------------------------------
// Everything a replay of the step would have to put back: every ball's position, velocity and spin, the engaged balls' contacts (the whole entry, copy-assigned in place so the list and the
// pointers to its entries stay as they are; the vectors inside keep their capacity), the ball pairs, the saved state of every moving part, the size of the event queue. The engaged list and the pair
// list only change in ApplyReleases, after the overlap check, so a snapshot of their contents is enough. Bats and bands are not saved yet: islands with them are not replayed.
namespace
{
   struct PartSnap
   {
      ContactMovingPart *p = nullptr;
      ContactMovingPart::State st;
   };
   struct BallSnap
   {
      HitBall *pb = nullptr;
      Vertex3Ds pos, vel, angm;
      float lastEventSqrDist = 0.f;
      bool advanced = false;
      bool hasKept = false;
      ContactSystemImpl::KeptDoubles kept;
      int sleepGroup = 0;   // pile sleep
      bool sleepOwe = false;
      double sleepHoldS = 0.0, sleepCreepS = 0.0;
      Vertex3Ds sleepPos;
   };
   struct StartSnap
   {
      std::vector<BallSnap> balls;
      std::vector<BallContact> contacts;
      std::vector<BallPairState> pairs;
      std::vector<PartSnap> parts;
      std::vector<ContactSystemImpl::SleepGroup> sleep;   // pile sleep
      double sleepCycleDt = 0.0;
      size_t deferred = 0;
      bool valid = false;
      unsigned us = 0;
   };
   StartSnap g_snap;
   bool g_snapSkipParts = false;   // test control (VPX_SNAPSHOT_TEST=2): the parts are not put back; the test must then see a difference wherever a part moved
}

void ContactSystem::SaveStart()
{
   const uint64_t t0 = usec();
   StartSnap &s = g_snap;
   const auto &vb = g_pplayer->m_vball;
   s.balls.resize(vb.size());
   for (size_t i = 0; i < vb.size(); i++)
   {
      HitBall &hb = vb[i]->m_hitBall;
      s.balls[i] = { &hb, hb.m_d.m_pos, hb.m_d.m_vel, hb.m_angularmomentum, hb.m_lastEventSqrDist, hb.m_rubberAdvanced };
      const auto kd = ContactSystemImpl::g_keptDoubles.find(&hb);
      s.balls[i].hasKept = kd != ContactSystemImpl::g_keptDoubles.end();
      if (s.balls[i].hasKept)
         s.balls[i].kept = kd->second;
      s.balls[i].sleepGroup = hb.m_sleepGroup;
      s.balls[i].sleepOwe = hb.m_sleepOwe;
      s.balls[i].sleepHoldS = hb.m_sleepHoldS;
      s.balls[i].sleepCreepS = hb.m_sleepCreepS;
      s.balls[i].sleepPos = hb.m_sleepPos;
   }
   s.sleep = ContactSystemImpl::g_sleep;
   s.sleepCycleDt = ContactSystemImpl::g_sleepCycleDt;
   s.contacts.resize(g_engaged.size());
   size_t k = 0;
   for (const BallContact &e : g_engaged)
      s.contacts[k++] = e;
   s.pairs = g_ballPairs;
   size_t np = 0;
   for (ContactHost *const o : ContactSystemImpl::Hosts())
      if (ContactMovingPart *const mp = o->MovingPart())
      {
         if (np >= s.parts.size())
            s.parts.emplace_back();
         s.parts[np].p = mp;
         mp->SaveState(s.parts[np].st);
         np++;
      }
   s.parts.resize(np);
   s.deferred = g_deferred.size();
   s.valid = true;
   s.us = (unsigned)(usec() - t0);
}

unsigned ContactSystem::LastSaveUs()
{
   return g_snap.us;
}

// Puts the snapshot back. Balls: all of them (only == nullptr), or the engaged ones and those in `only` (a replay: the free balls have been moved on, and stay)
static void RestoreSnap(const std::vector<HitBall *> *only)
{
   StartSnap &s = g_snap;
   if (!s.valid)
      return;
   ContactSystemImpl::g_sleep = s.sleep;   // pile sleep: every ball's, whatever the replay restores of the rest
   ContactSystemImpl::g_sleepCycleDt = s.sleepCycleDt;
   for (const BallSnap &b : s.balls)
   {
      b.pb->m_sleepGroup = b.sleepGroup;
      b.pb->m_sleepOwe = b.sleepOwe;
      b.pb->m_sleepHoldS = b.sleepHoldS;
      b.pb->m_sleepCreepS = b.sleepCreepS;
      b.pb->m_sleepPos = b.sleepPos;
      if (only != nullptr && b.pb->m_contact == nullptr && std::find(only->begin(), only->end(), b.pb) == only->end())
         continue;
      b.pb->m_d.m_pos = b.pos;
      b.pb->m_d.m_vel = b.vel;
      b.pb->m_angularmomentum = b.angm;
      b.pb->m_lastEventSqrDist = b.lastEventSqrDist;
      b.pb->m_rubberAdvanced = b.advanced;
      if (b.hasKept)
         ContactSystemImpl::g_keptDoubles[b.pb] = b.kept;
      else
         ContactSystemImpl::g_keptDoubles.erase(b.pb);
   }
   if (g_engaged.size() == s.contacts.size())
   {
      size_t k = 0;
      for (BallContact &e : g_engaged)
         e = s.contacts[k++];
   }
   else
      PLOGW << "Snapshot: the engaged list changed since it was taken (" << s.contacts.size() << " -> " << g_engaged.size() << "): contacts not restored";
   g_ballPairs = s.pairs;
   if (!g_snapSkipParts)
      for (const PartSnap &p : s.parts)
         p.p->RestoreState(p.st);
   if (g_deferred.size() > s.deferred)
      g_deferred.erase(g_deferred.begin() + (long)s.deferred, g_deferred.end());
   if (ContactSystemImpl::s_pendingCy != nullptr)   // the events of the step being undone
   {
      ContactSystemImpl::s_pendingCy->hitEvents.clear();
      ContactSystemImpl::s_pendingCy->partHits.clear();
      ContactSystemImpl::s_pendingCy->releases.clear();
      ContactSystemImpl::s_pendingCy->bandReleases.clear();
      ContactSystemImpl::s_pendingCy = nullptr;
   }
}

void ContactSystem::RestoreStart()
{
   RestoreSnap(nullptr);
}

bool ContactSystem::SnapshotTestOn()
{
   static const int level = []() { const char *const e = getenv("VPX_SNAPSHOT_TEST"); return e != nullptr ? atoi(e) : 0; }();
   g_snapSkipParts = level == 2;
   return level != 0;
}

namespace
{
   struct Fp
   {
      uint64_t balls = 14695981039346656037ull, contacts = 14695981039346656037ull, pairs = 14695981039346656037ull, parts = 14695981039346656037ull;
   };
   uint64_t Fnv(uint64_t h, const void *p, const size_t n)
   {
      const unsigned char *b = static_cast<const unsigned char *>(p);
      for (size_t i = 0; i < n; i++)
      {
         h ^= b[i];
         h *= 1099511628211ull;
      }
      return h;
   }
   template <class T> uint64_t FnvV(uint64_t h, const T &v) { return Fnv(h, &v, sizeof(T)); }
   Fp Fingerprint()
   {
      Fp f;
      for (const Ball *const b : g_pplayer->m_vball)
      {
         const HitBall &hb = b->m_hitBall;
         f.balls = FnvV(f.balls, hb.m_d.m_pos);
         f.balls = FnvV(f.balls, hb.m_d.m_vel);
         f.balls = FnvV(f.balls, hb.m_angularmomentum);
         f.balls = FnvV(f.balls, hb.m_lastEventSqrDist);
      }
      for (const BallContact &e : g_engaged)
      {
         f.contacts = FnvV(f.contacts, e.omega);
         f.contacts = FnvV(f.contacts, e.lastSpeed);
         f.contacts = FnvV(f.contacts, e.lastVn);
         f.contacts = FnvV(f.contacts, e.dwellS);
         f.contacts = FnvV(f.contacts, e.peakFn);
         f.contacts = FnvV(f.contacts, e.steps);
         for (const SurfaceContact &it : e.items)
         {
            f.contacts = FnvV(f.contacts, it.st.xi);
            f.contacts = FnvV(f.contacts, it.gs.touching);
            f.contacts = FnvV(f.contacts, it.gs.xi);
            if (it.gs.ra.th.Dot(it.gs.ra.th) + it.gs.ra.tw * it.gs.ra.tw > 0.0 || it.st.ra.th.Dot(it.st.ra.th) + it.st.ra.tw * it.st.ra.tw > 0.0)
            {
               f.contacts = FnvV(f.contacts, it.gs.ra.th);
               f.contacts = FnvV(f.contacts, it.gs.ra.tw);
               f.contacts = FnvV(f.contacts, it.st.ra.th);
               f.contacts = FnvV(f.contacts, it.st.ra.tw);
            }
            f.contacts = FnvV(f.contacts, it.gs.zero);
            f.contacts = FnvV(f.contacts, it.gs.reff);
            f.contacts = FnvV(f.contacts, it.gs.kb);
            f.contacts = FnvV(f.contacts, it.touching);
            f.contacts = FnvV(f.contacts, it.partTouching);
            f.contacts = FnvV(f.contacts, it.fnLast);
         }
      }
      for (const BallPairState &p : g_ballPairs)
      {
         f.pairs = FnvV(f.pairs, p.st.xi);
         if (p.st.ra.th.Dot(p.st.ra.th) + p.st.ra.tw * p.st.ra.tw > 0.0)
         {
            f.pairs = FnvV(f.pairs, p.st.ra.th);
            f.pairs = FnvV(f.pairs, p.st.ra.tw);
         }
         f.pairs = FnvV(f.pairs, p.inContact);
         f.pairs = FnvV(f.pairs, p.lastDelta);
         f.pairs = FnvV(f.pairs, p.zero);
      }
      ContactMovingPart::State st;
      for (ContactHost *const o : ContactSystemImpl::Hosts())
         if (ContactMovingPart *const mp = o->MovingPart())
         {
            mp->SaveState(st);
            f.parts = Fnv(f.parts, st.bytes, st.used);
         }
      return f;
   }
}

void ContactSystem::SnapshotTest(const float dtime)
{
   static long compared = 0, mismatches = 0, skipped = 0;
   if (!g_snap.valid)
      return;
   for (const BallContact &e : g_engaged)
      if (e.band != nullptr)
      {
         skipped++;
         return;
      }
   if (g_engaged.empty())
      return;   // nothing the step touched
   const Fp a = Fingerprint();
   RestoreStart();
   Step(dtime);
   const Fp b = Fingerprint();
   compared++;
   if (a.balls != b.balls || a.contacts != b.contacts || a.pairs != b.pairs || a.parts != b.parts)
   {
      mismatches++;
      PLOGW << "SnapshotTest t=" << g_pplayer->m_time_msec << ": the step taken again from the snapshot gave a different result in:" << (a.balls != b.balls ? " balls" : "") << (a.contacts != b.contacts ? " contacts" : "")
            << (a.pairs != b.pairs ? " pairs" : "") << (a.parts != b.parts ? " parts" : "") << " (" << g_engaged.size() << " engaged balls)";
   }
   if (compared % 2000 == 1)
      PLOGI << "SnapshotTest: " << compared << " steps compared, " << mismatches << " different, " << skipped << " skipped (band)";
}


// ---- The overlap check and the replay of the sub-interval (docs/UnifiedContactDesign.txt section 11, point 2) ------------------------------------------------------------------------
// Run after UpdateDisplacements, before the events are dispatched. A free ball found overlapping a moving part was born inside this sub-interval: it moved on a straight line (stock free flight)
// while the part moved by its own motion. Find the time t* in the sub-interval at which the gap was zero with BOTH of them where they really were, put the snapshot of the start back, step the
// system to t* (the engaged balls and the parts), engage the ball there with its true position and velocity, and step on to the end. The contact starts at zero depth with the real velocities.
// Several balls: the earliest first contact first, the others are found again on the new end state (chained, at most the free balls and parts near + 1 times).
namespace
{
   unsigned g_replayCount = 0, g_replayUs = 0;   // since the cost log last asked
}

void ContactSystem::TakeReplayStats(unsigned &count, unsigned &us)
{
   count = g_replayCount;
   us = g_replayUs;
   g_replayCount = g_replayUs = 0;
}

void ContactSystem::CheckOverlaps(const float hittime)
{
   ContactSystemImpl::CheckOverlaps(hittime);
}

void ContactSystemImpl::CheckOverlaps(const float hittime)
{
   using namespace RubberContact;
   static const bool on = []() { const char *const e = getenv("VPX_OVERLAP_REPLAY"); return !(e != nullptr && atoi(e) == 0); }();
   if (!on || !g_snap.valid || hittime <= 0.f)
      return;
   std::vector<ContactHost *> parts;   // the part hosts that can be overlapped (collidable); `allParts` below for the advance
   for (ContactHost *const o : Hosts())
      if (o->m_enabled && o->MovingPart() != nullptr)
         parts.push_back(o);
   if (parts.empty())
      return;
   const double tol = ContactSettings::Get().arrivalDepthUm * 1e-6;
   const double dtS = (double)hittime * DEFAULT_STEPTIME_S;
   auto isFree = [](const HitBall *pb) { return pb->m_contact == nullptr && pb->m_bandOwner == nullptr && !pb->m_d.m_lockedInKicker; };
   struct Cand
   {
      HitBall *pb;
      ContactHost *host;
      double delta;   // the overlap at the end of the sub-interval, m
   };
   std::vector<HitBall *> noPast;   // balls found overlapping with no earlier clean state (point 6): seated by minimum translation and engaged once the parts have advanced
   std::vector<std::pair<HitBall *, ContactHost *>> noPastHost;
   auto scan = [&](std::vector<Cand> &out) {
      out.clear();
      for (Ball *const b : g_pplayer->m_vball)
      {
         HitBall *const pb = &b->m_hitBall;
         if (!isFree(pb) || std::find(noPast.begin(), noPast.end(), pb) != noPast.end())
            continue;
         BallSI bs;
         bs.R = (double)pb->m_d.m_radius * U;
         const V3 c((double)pb->m_d.m_pos.x * U, (double)pb->m_d.m_pos.y * U, (double)pb->m_d.m_pos.z * U);
         for (ContactHost *const o : parts)
         {
            ContactMovingPart *const mp = o->MovingPart();
            mp->SetBall(bs);
            const TubeProbe pr = mp->Probe(c);
            if (pr.overlap && pr.delta > tol)
               out.push_back({ pb, o, pr.delta });
         }
      }
   };
   std::vector<Cand> cands;
   scan(cands);
   if (cands.empty())
      return;
   for (const BallContact &e : g_engaged)
      if (e.band != nullptr)
      {
         static uint32_t lastMs = 0;
         if (g_pplayer->m_time_msec - lastMs > 2000)
         {
            lastMs = g_pplayer->m_time_msec;
            PLOGW << "Overlap with a moving part found, but an island holds a band, whose state is not saved yet: no replay (the engagement of the next cycle seats the ball)";
         }
         return;
      }

   const uint64_t tStart = usec();
   struct Rp
   {
      HitBall *pb;
      ContactHost *host;
      double t;   // the share of the sub-interval at which the ball is engaged
   };
   std::vector<Rp> list;
   auto snapBall = [&](const HitBall *pb) -> const BallSnap * {
      for (const BallSnap &b : g_snap.balls)
         if (b.pb == pb)
            return &b;
      return nullptr;
   };
   // The system at share tEnd of the sub-interval: the snapshot put back, the engaged balls and the parts stepped to each engagement, the balls of `list` engaged there, then on to tEnd
   auto simulate = [&](const double tEnd) {
      std::vector<HitBall *> only;
      for (const Rp &e : list)
      {
         if (e.pb->m_contact != nullptr)
            Release(e.pb, "replay");   // engaged by an earlier simulation: the snapshot has no entry for it
         only.push_back(e.pb);
      }
      RestoreSnap(&only);
      double t = 0.0;
      bool first = true;
      for (Rp &e : list)
      {
         if (e.t > tEnd)
            break;
         if (e.t > t)
         {
            s_append = !first;
            ContactSystem::Step((float)((e.t - t) * (double)hittime));
            s_append = false;
            first = false;
            t = e.t;
         }
         for (ContactHost *const o : Hosts())   // a part no engaged ball held has not moved yet: it is at the start of the sub-interval. It moves on to the engagement, as it would have alone (gapAt measures it so)
            if (ContactMovingPart *const mp = o->MovingPart())
            {
               const double idle = dtS * e.t - mp->steppedS;
               if (idle > 1e-9)
               {
                  mp->AdvanceAlone(idle);
                  mp->steppedS += idle;   // (PartAdvance advances it for the rest of the sub-interval only)
               }
            }
         if (const BallSnap *const b0 = snapBall(e.pb))
            e.pb->m_d.m_pos = b0->pos + b0->vel * (float)(e.t * (double)hittime);   // its free flight to the engagement
         EngageImpl(e.host, e.pb, "overlap");
      }
      if (tEnd > t)
      {
         s_append = !first;
         ContactSystem::Step((float)((tEnd - t) * (double)hittime));
         s_append = false;
      }
   };
   // The gap between a free ball (on its straight flight from the snapshot) and the nearest part, at share t, in the system simulate(t) left. A part that was stepped for less than t is advanced alone
   // for the rest on a scratch copy of its state.
   auto gapAt = [&](HitBall *pb, const double t) {
      const BallSnap *const b0 = snapBall(pb);
      BallSI bs;
      bs.R = (double)pb->m_d.m_radius * U;
      const Vertex3Ds p = b0->pos + b0->vel * (float)(t * (double)hittime);
      const V3 c((double)p.x * U, (double)p.y * U, (double)p.z * U);
      double gmin = 1e9;
      for (ContactHost *const o : parts)
      {
         ContactMovingPart *const mp = o->MovingPart();
         const double idle = dtS * t - mp->steppedS;
         ContactMovingPart::State tmp;
         const bool adv = idle > 1e-9;
         if (adv)
         {
            mp->SaveState(tmp);
            mp->AdvanceAlone(idle);
         }
         mp->SetBall(bs);
         const TubeProbe pr = mp->Probe(c);
         if (adv)
            mp->RestoreState(tmp);
         gmin = std::min(gmin, pr.overlap ? -pr.delta : pr.dist - bs.R);
      }
      return gmin;
   };
   int evaluations = 0;
   // The share of the sub-interval at which the ball first touched: regula falsi on the gap with the system simulated at each trial time (the contact never starts deeper than the arrival depth)
   auto findT = [&](const Cand &cd, double &tOut) {
      simulate(0.0);
      evaluations++;
      const double g0 = gapAt(cd.pb, 0.0);
      if (g0 < -tol)
         return false;   // overlapping at the start of the sub-interval: no past
      double lo = 0.0, glo = g0, hi = 1.0, ghi = -cd.delta;
      for (int it = 0; it < 8; it++)
      {
         double t = (glo - ghi) > 1e-12 ? lo + (hi - lo) * glo / (glo - ghi) : 0.5 * (lo + hi);
         t = std::min(hi - 0.02 * (hi - lo), std::max(lo + 0.02 * (hi - lo), t));   // (keeps the bracket shrinking)
         simulate(t);
         evaluations++;
         const double g = gapAt(cd.pb, t);
         if (std::abs(g) <= tol)
         {
            tOut = t;
            return true;
         }
         if (g > 0.0)
         {
            lo = t;
            glo = g;
         }
         else
         {
            hi = t;
            ghi = g;
         }
      }
      tOut = lo;   // the side that does not overlap
      return true;
   };

   const int maxChain = (int)(g_pplayer->m_vball.size() + parts.size()) + 1;
   std::string what;
   for (int chain = 0; chain < maxChain && !cands.empty(); chain++)
   {
      Rp best = { nullptr, nullptr, 2.0 };
      double bestDelta = 0.0;
      for (const Cand &cd : cands)
      {
         double t = 0.0;
         if (!findT(cd, t))
         {
            noPast.push_back(cd.pb);
            noPastHost.emplace_back(cd.pb, cd.host);
            static uint32_t lastMs = 0;
            if (g_pplayer->m_time_msec - lastMs > 500)
            {
               lastMs = g_pplayer->m_time_msec;
               PLOGW << "Overlap with " << cd.host->m_name << " at the START of the sub-interval (" << cd.delta * 1e6 << " um at its end): no past to rewind to; the ball is seated by minimum translation (point 6)";
            }
            continue;
         }
         if (t < best.t)
         {
            best = { cd.pb, cd.host, t };
            bestDelta = cd.delta;
         }
      }
      if (best.pb == nullptr)
         break;
      list.push_back(best);
      std::sort(list.begin(), list.end(), [](const Rp &a, const Rp &b) { return a.t < b.t; });
      g_replayCount++;
      what += (what.empty() ? "" : "; ") + best.host->m_name + " at " + std::to_string((int)(best.t * 100.0)) + " % (" + std::to_string((int)(bestDelta * 1e6)) + " um deep at the end)";
      simulate(1.0);   // the end state with the new ball in: the other free balls are checked against it
      scan(cands);
   }
   if (list.empty())
   {
      simulate(1.0);   // the system as it was before (the same step again; nothing was engaged)
   }
   // the parts: advanced alone for what the contact system did not step them, as the movers do after a step
   for (ContactHost *const o : Hosts())
      if (o->MovingPart() != nullptr)
         o->PartAdvance(hittime);
   for (Ball *const b : g_pplayer->m_vball)
      b->m_hitBall.m_rubberAdvanced = false;   // every ball has been moved for this sub-interval
   for (const auto &nh : noPastHost)
      if (nh.first->m_contact == nullptr)
         EngageImpl(nh.second, nh.first, "no past");
   g_replayUs += (unsigned)(usec() - tStart);
   if (!list.empty())
   {
      static uint32_t lastMs = 0;
      if (g_pplayer->m_time_msec - lastMs > 500)
      {
         lastMs = g_pplayer->m_time_msec;
         PLOGW << "Overlap replay t=" << g_pplayer->m_time_msec << " ms: " << what << "; " << evaluations << " evaluations, " << (usec() - tStart) << " us";
      }
      if (list.size() > 2 && g_pplayer->m_liveUI != nullptr)
      {
         static bool shown = false;
         if (!shown)
         {
            shown = true;
            g_pplayer->m_liveUI->PushNotification("Several balls overlapped moving parts in one sub-interval and were replayed (see vpinball.log, 'Overlap replay')", 8000);
         }
      }
   }
}

// VPX_TRACE_BALL=<driver ball name> (a diagnostic, off by default): that ball's position (U) and velocity at a named stage of the cycle, and, when asked, its contact's items and its
// island, one log line each. Nothing is printed or changed without the variable.
void ContactTraceBall(const char *const stage, const HitBall *const pb, const bool items, const int island)
{
   static const std::string want = []() { const char *const e = getenv("VPX_TRACE_BALL"); return std::string(e != nullptr ? e : ""); }();
   if (want.empty() || pb == nullptr || g_pplayer == nullptr)
      return;
   const char *const name = TestDriver::BallName(pb);
   if (name == nullptr || want != name)
      return;
   std::ostringstream o;
   o << std::fixed << std::setprecision(4) << "TraceBall " << name << " t=" << g_pplayer->m_time_msec << ' ' << stage << " pos " << pb->m_d.m_pos.x << ' ' << pb->m_d.m_pos.y << ' ' << pb->m_d.m_pos.z
     << " vel " << pb->m_d.m_vel.x << ' ' << pb->m_d.m_vel.y << ' ' << pb->m_d.m_vel.z;
   if (items)
   {
      o << " held " << (pb->m_contact != nullptr ? 1 : 0) << " band " << (pb->m_bandOwner != nullptr ? 1 : 0) << " island " << island << " items";
      if (pb->m_contact != nullptr)
         for (const SurfaceContact &it : pb->m_contact->items)
            o << ' ' << it.r->Name();
   }
   PLOGI << o.str();
}
