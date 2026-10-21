// license:GPLv3+

#include "core/stdafx.h"
#include <algorithm>
#include <cmath>
#include <list>
#include <sstream>

#include "parts/ball.h"
#include "parts/flipper.h"
#include "contactsystem.h"
#include "hitrubber.h"
#include "hitball.h"
#include "hitflipper.h"
#include "magnet.h"
#include "magnetballs.h"
#include "turntable.h"
#include "contactprof.h"
#include "contacthingepart.h"
#include <map>
#include <unordered_map>
#include "physicsini.h"
#include "PhysicsEngine.h"
#include "input/InputManager.h"

namespace
{
   constexpr double U = ContactSystem::kMetresPerUnit;
   constexpr double VS = ContactSystem::kVelocityScale;
   constexpr double kMassUnitKg = ContactSystem::kMassUnitKg;
}

namespace ContactSystem
{
   static bool s_unified = false, s_unifiedSet = false;

   void SetTable(const std::filesystem::path &tableFile)
   {
      const std::filesystem::path file = PhysicsIni::FindUnified(tableFile);
      PhysicsIni::Features f;
      const bool found = !file.empty() && PhysicsIni::ReadFeatures(file, f);
      s_unified = found && f.unified;   // no file, no [Features] or no  unified = 1  = stock VPX physics
      s_unifiedSet = true;
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
// HandBack or WriteBack -> Finish (releases, pair bookkeeping, deferred events).
using namespace RubberContact;

// The balls the contact system holds in an engine-stepped contact (a post, a wall rubber, the playfield, a disc): the rubbers each one touches, the shear state and the
// bookkeeping for the log and the release. HitBall::m_contact points at the ball's entry.
static std::list<BallContact> g_engaged;
static ContactSystem::CycleStats g_last;
static std::vector<std::function<void()>> g_deferred;      // events raised inside the step, fired after it (Finish)
void ContactSystem::DeferEvent(std::function<void()> f) { g_deferred.push_back(std::move(f)); }
static constexpr float kFlipperHorizonVpt = 1.0f;   // trace only: VPX time units (1 VPT = 10 ms) a bat is looked ahead for in the flip_clear numbers

struct ContactSystemImpl
{
   static void StepBalls(float dtime);
   static void Release(HitBall *ball, const char *why);
   static void Regather(BallContact &e, const Vertex3Ds &pos, float radiusU);
   static bool EngageImpl(HitRubber *host, FlipperMoverObject *bat, HitBall *pb, const char *how);
   static void EngageNear(float dtime);
   static void DetachBat(BallContact &e);

   // A ball in the step: its SI state, the substep's forces, and what the step needs to know about it.
   struct Body
   {
      HitRubber *o = nullptr;   // the rubber that engaged the ball (its settings drive the step)
      BallContact *e = nullptr;    // its engagement; null = a free ball that only joined an island through a pair
      HitBall *pb = nullptr;
      BallSI bs;
      V3 c, v, w;
      V3 extBaseVS;             // m/s per physics step: gravity and nudge (the part that does not depend on the position)
      Vertex3Ds posStart;
      double lScale = 1.0;
      double vCap = 50.0;
      int steps = 1;
      bool magnetsOn = false;
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
      bool tripped = false;
      const char *tripWhy = "";
      int tripStep = 0;
      // a bat engaged with the ball (its rubber is one more surface; the bat is stepped with the island)
      FlipperMoverObject *bat = nullptr;
      bool batLive = false;     // the bat's rubber is still in contact bookkeeping (false once released mid-cycle)
      RubberParams batRp;
      Outline batOl;
      V3 absCorr;               // position correction at the end when the start was rewound against a moving bat
      double pivOff = 0.0;      // how far the pivot has moved (nudge) since the start of the cycle, m
      double ctq = 0.0;         // torque on the bat from the ball in this substep, closing direction
      bool hitStop = false;     // the bat's plunger reached the coil stop during this cycle
      bool noPairs = false;     // a flipper's ball: ball-to-ball stays with the stock collision
      // a band whose nodes the ball loads (the band's nodes are stepped with the island)
      HitRubber *band = nullptr;
      BandBall *bb = nullptr;
      bool bandContact = false; // the ball touched the band in this cycle
      double lastGap = 0.0, relNormalVel = 0.0;   // for the band release test
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
      int stepsPair = 1;                 // the substeps its ball pairs need alone (the arrival-time loop adds the surfaces' own)
      float dtimeF = 0.0f;               // the cycle in engine time units
      std::vector<HitRubber *> bands;    // bands begun this cycle
      std::vector<HitRubber *> islandBands;   // the current island's bands
      std::vector<std::pair<int, int>> links;
      std::vector<HingePart::Part *> islandParts;   // the gates and spinners the current island steps (contacthingepart.h)
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
      static void Begin(Body &w, const double dt)
      {
         HitBall *const pb = w.pb;
         // A free ball (no contact of its own, in an island only through a ball pair) has nothing of the stock engine under it: its gravity, nudge and magnet pull act in the substeps like an engaged ball's.
         const bool onPlayfieldModel = w.e == nullptr || std::any_of(w.e->items.begin(), w.e->items.end(), [](const SurfaceContact &it) { return it.r->IsFloor(); });
         const Vertex3Ds extInc = onPlayfieldModel ? pb->ExternalVelocityIncrement() : Vertex3Ds(0.f, 0.f, 0.f);   // what the engine added this step: gravity, nudge, magnets at the start position
         const Vertex3Ds extBase = onPlayfieldModel ? pb->ExternalVelocityIncrement(false) : Vertex3Ds(0.f, 0.f, 0.f);   // the part that does not depend on the position
         w.magnetsOn = onPlayfieldModel && pb->m_magnetic && Magnets::Active();   // a magnet's pull follows the ball: evaluated at every substep, not once per cycle
         w.magParts = w.magnetsOn && Magnets::HasBallBallParts(pb);
         w.magAcc = w.magnetsOn ? V3((extInc.x - extBase.x) * VS, (extInc.y - extBase.y) * VS, (extInc.z - extBase.z) * VS).Length() / PHYSICS_STEPTIME_S : 0.0;
         const double cycleShare = dt / PHYSICS_STEPTIME_S;
         w.extBaseVS = V3(extBase.x * VS, extBase.y * VS, extBase.z * VS);
         w.v = w.v - V3(extInc.x * VS * cycleShare, extInc.y * VS * cycleShare, extInc.z * VS * cycleShare);
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
            // Eddy-current braking of the ball's spin across the coil field: a steel sphere turning at w in a field B feels a torque (2 pi / 15) sigma R^5 B^2 w
            // against the part of w across B (a spin about the field lines induces nothing). Exact exponential decay over the substep.
            if (w.eddyRate < 0.0)   // the field changes little within a cycle: one evaluation, the substeps reuse its rate and direction
            {
               w.eddyRate = 0.0;
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
         if (it.candR >= 0.0 && it.cand.empty())
            continue;
         const TubeProbe pr = it.r->Probe(w.bs, it, w.c);
         if (pr.dist <= 0.0)
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
         HitRubber *const o = e.host;   // null for a ball only a bat holds
         if (o != nullptr && !o->m_set.contact)   // switched off in the ini while a ball is engaged
         {
            Release(e.ball, "rubber contact switched off");
            continue;
         }
         if (e.bat != nullptr)
         {
            e.bat->FlipperEnsureInit();
            if (!e.bat->m_bat.set.rubberContact)   // switched off in the ini while a ball is engaged
               DetachBat(e);
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
            // The settings that set the substep are those of the finest surface the ball touches or can reach: the one that engaged it and every surface in its contact set. A ball only
            // a bat engaged has no engaging surface; stepped once per cycle after the bat let go (measured: replay of the 2026-10-13 launch, 1 ms step against the playfield, 28 um deep
            // at the start of the next, 247 N held for 1 ms, +3 m/s) it is the contact's own stiffness that must limit the step, whatever engaged the ball.
            const RubberObjectSettings *so = o != nullptr ? o->Settings() : nullptr;   // (a gate or spinner has none: its step is its own, SetupIsland)
            for (const SurfaceContact &it : e.items)
               if (const RubberObjectSettings *const s = it.r->Settings())
                  if (so == nullptr || s->substepUs < so->substepUs)
                     so = s;
            const int stepsN = so != nullptr ? std::max(1, (int)std::ceil(dt / (so->substepUs * 1e-6))) : 1;
            // A ball at rest on the rubber (slow at both ends of the cycle) does not need the impact step: the contact's own frequency sets the stable step (omega * h <= 0.2),
            // up to rest_substep_us. The rest depth is the force balance, not a property of the step.
            int steps = stepsN;
            if (e.bat != nullptr)
               SetupBat(w, e, dt);   // (rewinds the start against a moving bat: before the step count is chosen)
            if (so != nullptr)
            {
               const RubberObjectSettings &set = *so;
               if (set.restSubstepUs > set.substepUs && w.v.Length() < set.restSpeedMmS * 1e-3 && e.lastSpeed < set.restSpeedMmS * 1e-3 && e.omega > 0.0)
               {
               const double hRest = std::min(set.restSubstepUs * 1e-6, set.restOmegaH / e.omega);
               steps = std::max(1, std::min(stepsN, (int)std::ceil(dt / hRest)));
               }
               else if (set.restSubstepUs > set.substepUs && set.rollTravelUm > 0.0 && e.omega > 0.0 && ApproachSpeed(w, e, dt) < set.restSpeedMmS * 1e-3)
               {
               // Rolling or sliding along a surface (nothing approaches it): the contact's stiffness and the distance covered per step limit the step, not the impact step.
               const double hRoll = std::min({ set.restSubstepUs * 1e-6, set.restOmegaH / e.omega, set.rollTravelUm * 1e-6 / std::max(w.v.Length(), 1e-6) });
               steps = std::max(1, std::min(stepsN, (int)std::ceil(dt / hRoll)));
               }
            }
            if (w.bat != nullptr)
            {
               // The bat's own substep while a ball is on it; a ball at rest on a still bat steps by the contact's stiffness, as on any surface (the rest limits are the ball's other
               // surface's, else the usual ones)
               const double hBat = std::min(2e-5, w.bat->m_bat.set.rubberSubstepUs * 1e-6);
               const double restSpeed = (so != nullptr ? so->restSpeedMmS : 30.0) * 1e-3, restSub = (so != nullptr ? so->restSubstepUs : 50.0) * 1e-6, restOmegaH = so != nullptr ? so->restOmegaH : 0.8;
               const double batSpeed = std::abs((double)w.bat->m_bat.omega) * (double)w.bat->m_flipperradius * U;
               double h = hBat;
               if (restSub > hBat && w.v.Length() < restSpeed && e.lastSpeed < restSpeed && batSpeed < restSpeed && e.omega > 0.0)
                  h = std::max(hBat, std::min(restSub, restOmegaH / e.omega));
               steps = std::max(steps, (int)std::ceil(dt / h));
            }
            if (e.band != nullptr)
            {
               SetupBand(cy, w, e);
               if (w.band != nullptr)   // the band's substep
                  steps = std::max(steps, w.band->m_bc.steps);
            }
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
         if (pb->m_bandOwner != nullptr || pb->m_contact != nullptr || !PairEligible(pb))
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
               if (gap < reach && !cy.bodies[i].noPairs && !cy.bodies[j].noPairs)
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
      for (int i = 0; i < nw; i++)   // balls touching the same gate or spinner are one island: the part takes all their reactions
         if (cy.bodies[i].e != nullptr)
            for (const SurfaceContact &it : cy.bodies[i].e->items)
               if (HingePart::Part *const hg = it.r->MovingPart())
               {
                  for (int j = 0; j < i; j++)
                     if (cy.bodies[j].e != nullptr && std::any_of(cy.bodies[j].e->items.begin(), cy.bodies[j].e->items.end(), [hg](const SurfaceContact &o) { return o.r->MovingPart() == hg; }))
                     {
                        cy.links.emplace_back(j, i);
                        break;
                     }
               }
      for (int i = 0; i < nw; i++)   // balls on the same bat are one island: the bat takes all their torques
         if (cy.bodies[i].bat != nullptr)
            for (int j = 0; j < i; j++)
               if (cy.bodies[j].bat == cy.bodies[i].bat)
               {
                  cy.links.emplace_back(j, i);
                  break;
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

      int steps = 1;
      for (const int i : isl.members)
         steps = std::max(steps, cy.bodies[i].steps);
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
            const BallPairState &ps = *cy.pst[(size_t)(&p - cpairs.data())];
            const double closing = std::max(0.0, -(A.v - B.v).Dot(d / dist));
            // a pair that cannot touch within this cycle (the gap is more than the closing distance plus a margin for the force's own acceleration) needs no step of its own:
            // it joined the cluster only because it is within reach
            if (!ps.inContact && dist - A.bs.R - B.bs.R > closing * dt * 1.5 + kApproachMarginM)
               continue;
            const double reff = A.bs.R * B.bs.R / (A.bs.R + B.bs.R), mRed = A.bs.m * B.bs.m / (A.bs.m + B.bs.m);
            const double dImpact = std::pow(15.0 * mRed * closing * closing / (16.0 * Estar * std::sqrt(reff)), 0.4);
            const double delta = std::max({ ps.lastDelta, dImpact, 1e-6 });
            const double omega = std::sqrt(2.0 * Estar * std::sqrt(reff * delta) / mRed);
            const double omegaH = closing > kRestSpeedFree ? kImpactOmegaH : kRestOmegaH;
            const double h = std::min(kMaxPairStepS, std::max(cs.ballSubstepUs * 1e-6, omegaH / omega));
            steps = std::max(steps, (int)std::ceil(dt / h));
            stepsPair = std::max(stepsPair, (int)std::ceil(dt / h));
         }
      }
      cy.steps = steps;
      cy.stepsPair = stepsPair;
      // The gates and spinners of the island. Their step is the reduced-mass rule of a moving body (design 6.3): every ball's contact with a part is resolved for the reduced mass of
      // the ball and the part's effective mass at the contact, with the series stiffness of contact and bend, and ends at the arrival (Part::BallStep); the part's own forces
      // (pendulum, stops, switch spring) are resolved in its own substeps inside each island substep (Part::StepInIsland). A heavy moving body (the turntable's disc: the normal of a ball on it is along its axis, so the normal force
      // cannot turn it; a bat: hundreds of grams at the contact against an 80 g ball) has an effective mass far above the ball's and is not affected.
      cy.islandParts.clear();
      for (const int i : isl.members)
         if (cy.bodies[i].e != nullptr)
            for (const SurfaceContact &it : cy.bodies[i].e->items)
               if (HingePart::Part *const hg = it.r->MovingPart())
               {
                  if (std::find(cy.islandParts.begin(), cy.islandParts.end(), hg) == cy.islandParts.end())
                  {
                     cy.islandParts.push_back(hg);
                     hg->body.SetGravity(HitRubber::PartGravitySI());   // (its own forces are resolved in its own substeps: Part::StepInIsland)
                  }
                  const Body &w = cy.bodies[i];
                  hg->SetBall(w.bs);
                  steps = std::max(steps, (int)std::ceil(dt / hg->BallStep(w.c, w.v, dt) - 1e-9));
               }
      cy.steps = steps;
      cy.islandBands.clear();
      for (const int i : isl.members)
      {
         HitRubber *const b = cy.bodies[i].band;
         if (b != nullptr && std::find(cy.islandBands.begin(), cy.islandBands.end(), b) == cy.islandBands.end())
         {
            cy.islandBands.push_back(b);
            b->m_bc.steps = steps;   // the band steps with the island
            b->m_bc.h = cy.dt / steps;
         }
      }
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

   // ---- A bat engaged with a ball: its rubber is one more surface of the ball's contact, and the bat is stepped with the island ----
   // (the contact law is the bat's own, RubberContact::Evaluate: the bat is a moving surface whose reaction is a torque on it)

   // Per cycle: the bat's outline and rubber, the pivot, and, for the first step of an engagement, the seat of the starting overlap.
   static void SetupBat(Body &w, BallContact &e, const double dt)
   {
      FlipperMoverObject &m = *e.bat;
      w.bat = &m;
      w.batLive = true;
      w.noPairs = !ContactSystem::Unified();   // (the bat is stepped with the island, so its balls take part in pair contact)
      w.batRp = FlipperModel::MakeRubberParams(m.m_bat.set);
      w.batOl.rb = (double)m.m_hitcircleBase.radius * U;
      w.batOl.re = (double)m.m_endradius * U;
      w.batOl.Lc = (double)m.m_flipperradius * U;
      if (!e.batJustEngaged)
         return;
      e.batJustEngaged = false;
      // The ball starts the engagement inside the rubber: the hit was found late. Go back to the instant it first touched: in the bat's frame, along
      // the relative velocity, by depth / closing speed, at most one physics step. Moving the ball along its relative velocity puts it where it was
      // when it touched the bat's earlier pose; the ball's own absolute path is then restored at the end by the bat's displacement over that time.
      // If it is not closing (pushed in by something else) it is seated at touching distance along the normal. Nothing is stored as compression.
      (void)dt;
      constexpr double kMaxRewindS = PHYSICS_STEPTIME_S;
      const double dirSign = m.m_direction ? 1.0 : -1.0;
      const V3 pivotW((double)m.m_hitcircleBase.center.x * U, (double)m.m_hitcircleBase.center.y * U, 0.0);
      V3 c = w.c - pivotW;
      const double phi0 = (double)m.m_angleStart + dirSign * m.m_bat.phi - M_PI / 2.0;
      const V3 ez(0.0, 0.0, 1.0);
      double rewoundS = 0.0;
      for (int pass = 0; pass < 3; pass++)
      {
         const Probe pr = ProbeSurface(w.batRp, w.bs, w.batOl, phi0, c);
         if (!pr.overlap || pr.delta <= 1e-6 || pr.dist <= 1e-12)
            break;
         const V3 nrm = (c - pr.point) / pr.dist;
         const V3 vs = ez.Cross(pr.point) * (dirSign * m.m_bat.omega);   // bat surface velocity at the contact point
         const V3 vrel = w.v - vs;
         const double closing = -vrel.Dot(nrm);
         if (closing > 0.01)
         {
            const double t = std::min(pr.delta / closing, kMaxRewindS);
            c = c - vrel * t;
            w.absCorr = w.absCorr - vs * t;
            rewoundS += t;
         }
         else
            c = c + nrm * pr.delta;   // seat
      }
      w.c = c + pivotW;
      if (rewoundS > 0.0)
         PLOGI << "Rubber: " << m.m_pflipper->GetName() << " rewound the ball " << rewoundS * 1e3 << " ms to its first touch";
   }

   // The bat rubber's force on the ball at substep s (and the torque on the bat), the release bookkeeping, and the pivot's nudge motion
   static void BatForces(Body &w, const double h, const int s)
   {
      CPROF(BatForces);
      w.ctq = 0.0;
      if (w.bat == nullptr || !w.batLive || w.e == nullptr)
         return;
      FlipperMoverObject &m = *w.bat;
      BallContact &e = *w.e;
      const double dirSign = m.m_direction ? 1.0 : -1.0;
      const FlipperSettings &fs = m.m_bat.set;
      // pivot velocity (flipper nudge): the contact sees the ball's velocity relative to it
      double nzD = 0.0, nzV = 0.0;
      if (m.m_bat.nzActive)
         m.NudgeEval(m.m_bat.nzT + (s + 0.5) * h, nzD, nzV);
      const V3 pivotW((double)m.m_hitcircleBase.center.x * U + w.pivOff, (double)m.m_hitcircleBase.center.y * U, 0.0);
      const V3 cRel = w.c - pivotW;
      const V3 vRel = w.v - V3(nzV, 0.0, 0.0);
      const double phiMath = (double)m.m_angleStart + dirSign * m.m_bat.phi - M_PI / 2.0;   // engine angle 0 points along (0,-1)
      const ContactResult r = Evaluate(w.batRp, w.bs, w.batOl, phiMath, dirSign * m.m_bat.omega, cRel, vRel, w.w, h, e.batSt);
      if (r.inContact)
      {
         m.m_rbLogSteps++;
         m.m_rbLogPeakFn = std::max(m.m_rbLogPeakFn, r.fn);
         m.m_rbLogPeakDelta = std::max(m.m_rbLogPeakDelta, r.delta);
         w.F = w.F + r.force;
         w.Mb = w.Mb + r.torqueOnBall;
         w.ctq = dirSign * r.torqueOnBat;   // the bat step wants the torque in the closing direction
         w.deepest = std::max(w.deepest, r.delta);
         w.vnLast = std::max(w.vnLast, std::abs(w.v.x * r.normal.x + w.v.y * r.normal.y + w.v.z * r.normal.z));
         if (s == 0 && r.fn > 0.0 && r.delta > 1e-9)
            w.omega2 += 3.0 * r.fn / r.delta / w.bs.m;   // the bat rubber's stiffness (Hertz 1.5 Fn / delta, doubled for the layer's stiffening) is part of what limits the step, like a surface's
      }
      // Resting regime: a ball hovering a few micrometres off the rubber is still engaged. It is let go only once it is clear by more than
      // rubber_release_gap_mm as well as receding (release hysteresis), else it is released and re-hit by the engine dozens of times a second.
      const bool receding = r.dist > e.batLastDist;
      e.batLastDist = r.dist;
      const bool clearOfRubber = !r.inContact && r.dist - w.bs.R > fs.rubberReleaseGapMm * 1e-3;
      e.batClear = (clearOfRubber && receding) ? e.batClear + 1 : 0;
      if (e.batClear >= (int)fs.rubberReleaseSteps)
      {
         DetachBat(e);
         w.batLive = false;
      }
      w.pivOff += nzV * h;   // the pivot has moved on by the end of the substep
   }

   // ---- A ball on a band: the band's nodes are stepped with the island, the ball's contact with the band is one more force on the ball ----

   // Per cycle, the first ball of a band begins the band's cycle (the sling's driver and switches, the substep, the follow nodes); then this ball is attached to its BandBall.
   static void SetupBand(Cycle &cy, Body &w, BallContact &e)
   {
      HitRubber *const band = e.band;
      if (std::find(cy.bands.begin(), cy.bands.end(), band) == cy.bands.end())
      {
         cy.bands.push_back(band);
         band->BandBeforeCore(cy.dtimeF);
         band->m_bandIslandStepped = true;   // BandAdvance later this cycle only finishes up
         band->BandCycleBegin(cy.dtimeF);    // may let go of the band's balls (switched off in the ini)
      }
      if (e.band == nullptr)
         return;
      for (BandBall &b : band->m_engaged)
         if (b.ball == e.ball)
         {
            w.band = band;
            w.bb = &b;
            w.noPairs = !ContactSystem::Unified();   // (the band is stepped with the island, so its balls take part in pair contact)
            return;
         }
   }

   // The band's force on the ball at substep s; the reaction goes into the band's nodes
   static void BandForces(Body &w)
   {
      CPROF(BandForces);
      if (w.band == nullptr || w.bb == nullptr)
         return;
      HitRubber::BandBallContact out;
      w.band->BandBallForce(*w.bb, w.bs, w.c, w.v, w.w, out);
      if (out.inContact)
      {
         w.bandContact = true;
         w.F = w.F + out.force;
         w.Mb = w.Mb + out.torque;
         w.deepest = std::max(w.deepest, out.delta);
      }
      w.lastGap = out.gap;
      w.relNormalVel = (w.v - out.vs).Dot(out.normal);
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
         if (it.candR >= 0.0 && it.cand.empty())   // a surface of faces with none within reach of the ball this cycle: no contact
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
         V3 vsurf, wsurf;   // velocity of the surface at the contact and its angular velocity (a moving surface: a turntable's disc)
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
         const double delta = A.bs.R + B.bs.R - dist;
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
         // safety net: a non-finite state, a runaway speed or a ball sunk deep into a rubber means the model left its range (a ball wedged in a loop tighter than
         // itself, a degenerate rubber). Hand the ball back untouched; the stock collision takes over for a few ms.
         if (!std::isfinite(w.c.x + w.c.y + w.c.z + w.v.x + w.v.y + w.v.z + w.w.x + w.w.y + w.w.z) || w.v.Length() > w.vCap || w.deepest > 0.8 * w.bs.R)
         {
            w.tripped = true;
            w.tripWhy = !std::isfinite(w.c.x + w.c.y + w.c.z + w.v.x + w.v.y + w.v.z + w.w.x + w.w.y + w.w.z) ? "not finite" : w.v.Length() > w.vCap ? "speed" : "depth";
            w.tripStep = s;
            tripped = true;
         }
      }
      StepBats(cy, isl, h, true);   // the bat takes the balls' reaction in the same substep (no contact torque from a bad contact)
      return tripped;
   }

   // Every bat of the island steps ONCE per substep, with the torque of all the balls on it (several balls can be on one bat); withTorque false = the rest of the cycle after the balls were handed back
   static void StepBats(Cycle &cy, const ContactSystem::Islands::Island &isl, const double h, const bool withTorque)
   {
      for (const int i : isl.members)
      {
         Body &w = cy.bodies[i];
         if (w.bat == nullptr)
            continue;
         bool seen = false;
         for (const int j : isl.members)
            if (j == i)
               break;
            else if (cy.bodies[j].bat == w.bat)
            {
               seen = true;
               break;
            }
         if (seen)
            continue;
         double tq = 0.0;
         if (withTorque)
            for (const int j : isl.members)
               if (cy.bodies[j].bat == w.bat && !cy.bodies[j].tripped)
                  tq += cy.bodies[j].ctq;
         w.bat->FlipperStep(h, tq);
         const bool stop = w.bat->m_bat.plungerHitStop;
         for (const int j : isl.members)
            if (cy.bodies[j].bat == w.bat)
               cy.bodies[j].hitStop |= stop;
      }
   }

   // ---- Arrival-time substeps (islands of balls and static surfaces) ----
   // A step never carries a ball into a surface. At the start of every substep each ball is probed against the surfaces of its contact; a surface that is not touching yet bounds the step
   // by the time the ball needs to get within kArriveM of it (distance, closing speed, and the acceleration it can pick up), and a touching one by its own stiffness at the depth the
   // contact will reach (omega * h limited, Hertz). A ball with nothing near takes the rest of the cycle as one step. The ball pairs keep the step their own rule gives (cy.stepsPair).

   static bool IslandStatic(const Cycle &cy, const ContactSystem::Islands::Island &isl)
   {
      for (const int i : isl.members)
      {
         const Body &w = cy.bodies[i];
         if (w.bat != nullptr || w.band != nullptr)
            return false;
         if (w.e != nullptr)
            for (const SurfaceContact &it : w.e->items)
            {
               if (it.r->MovingPart() != nullptr)
                  return false;
               if (const HitRubber *const hr = dynamic_cast<const HitRubber *>(it.r))
                  if (hr->m_turntable >= 0)
                     return false;
            }
      }
      return true;
   }

   static void ProbeBody(Body &w)
   {
      CPROF(ProbeBody);
      w.probesReady = false;
      if (w.e == nullptr)
         return;
      w.probes.resize(w.e->items.size());
      size_t k = 0;
      for (SurfaceContact &it : w.e->items)
      {
         TubeProbe &pr = w.probes[k++];
         if (it.candR >= 0.0 && it.cand.empty())   // no face within reach this cycle: far
         {
            pr = TubeProbe();
            pr.dist = -1.0;
            continue;
         }
         pr = it.r->Probe(w.bs, it, w.c);
      }
      w.probesReady = true;
   }

   // The next substep of an island: at most `remaining`
   static double NextStep(Cycle &cy, const ContactSystem::Islands::Island &isl, const double remaining, int &arrivals)
   {
      CPROF(NextStep);
      constexpr double kTouchM = 1.0e-6;          // within this gap a surface counts as touching
      constexpr double kArriveM = 0.5e-6;         // an arrival step ends this far above the surface
      constexpr double kMinStepS = 0.2e-6, kMaxContactStepS = 100e-6;
      constexpr double kImpactOmegaH = 0.2, kRestOmegaH = 0.8;
      constexpr double kRestSpeed = 0.02;   // as the ball pairs' rule
      // resting = slow AND the forces on the ball balance: what is left over (m times its acceleration over the last substep) is at most this share of the contact's own force.
      // At rest the contact carries everything that presses the ball into it (weight, magnets, balls on top); at the deepest point of an impact it far exceeds the rest.
      constexpr double kBalance = 0.2;   // measured 2026-10-21 (6 balls on the magnet): the strongest contact of a resting ball shows net force / contact force 0.005 median and 0.01 at most; in the bursts after a landing 99 % of cycles are under 0.08. 0.2 is 20x the resting maximum
      double h = remaining;
      if (cy.stepsPair > 1)
      {
         h = std::min(h, cy.dt / cy.stepsPair);
      }
      for (const int i : isl.members)
      {
         Body &w = cy.bodies[i];
         if (w.e == nullptr || !w.probesReady)
            continue;
         // the acceleration the ball can pick up: gravity and the nudge, a magnet's pull, and what the contact forces of the last substep gave it
         const double fLeft = w.hS > 0.0 ? (w.v - w.vS).Length() / w.hS * w.bs.m : -1.0;   // N, the ball's net force over the last substep (< 0 = unknown)
         const double aBound = w.extBaseVS.Length() / PHYSICS_STEPTIME_S + w.magAcc + w.F.Length() / w.bs.m;
         size_t k = 0;
         for (SurfaceContact &it : w.e->items)
         {
            const TubeProbe &pr = w.probes[k++];
            if (!pr.overlap && pr.dist < 0.0)
               continue;   // far
            const double gap = pr.overlap ? -pr.delta : pr.dist - (w.bs.R + pr.rho);
            const double vn = -w.v.Dot(pr.normal);   // (static surface)
            if (gap > kTouchM)
            {
               it.hcStep = 0.0;
               const double d = gap - kArriveM, vc = std::max(vn, 0.0);
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
                  arrivals++;
               }
            }
            else
            {
               static const bool noLand = getenv("VPX_ARRIVAL_NOSPLIT") != nullptr;
               if (!noLand && !pr.overlap && vn > 1e-9 && gap > 0.0)   // inside the touch gap but not touching yet: land on the surface exactly, so the force starts at the touch
               {
                  const double tl = aBound > 1e-9 ? (-vn + std::sqrt(vn * vn + 4.0 * aBound * gap)) / (2.0 * aBound) : gap / vn;
                  if (tl < h)
                  {
                     h = std::max(tl, kMinStepS);
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
               h = std::min(h, std::max(hc, kMinStepS));
            }
         }
      }
      return std::min(remaining, std::max(h, kMinStepS));
   }

   // The state a substep changes, kept at its start so the substep can be taken again (a contact that started too deep, or a trip): the balls, the shear springs and touching flags
   // of their contacts, the ball pairs, the events queued. The gap and closing speed of each surface at that point say whether the next substep started a contact late.
   struct RwSnap
   {
      struct Bd { V3 c, v, w; double omega2, deepest, vnLast, peakFn; size_t itemBase; size_t nItems; };
      struct It { V3 xi; bool touching; };
      struct Pr { V3 xi; bool inContact; double maxDelta, maxFn; };
      std::vector<Bd> b;
      std::vector<It> it;
      std::vector<Pr> p;
      std::vector<double> gap, vn;   // per surface of each ball (flat, in the order of b): the gap and the closing speed at the snapshot (1e9 / 0 = far)
      std::vector<double> pairGap;
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
               sn.it.push_back({ it.st.xi, it.touching });
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
         sn.p.push_back({ ps.st.xi, ps.inContact, ps.maxDelta, ps.maxFn });
         const Body &A = cy.bodies[isl.edges[k].first], &B = cy.bodies[isl.edges[k].second];
         sn.pairGap.push_back((A.c - B.c).Length() - A.bs.R - B.bs.R);
      }
      sn.partHits = cy.partHits.size();
      sn.hitEvents = cy.hitEvents.size();
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
               it.touching = sn.it[k].touching;
               k++;
            }
         }
      }
      for (size_t k = 0; k < isl.edges.size(); k++)
      {
         BallPairState &ps = *cy.pst[k];
         ps.st.xi = sn.p[k].xi;
         ps.inContact = sn.p[k].inContact;
         ps.maxDelta = sn.p[k].maxDelta;
         ps.maxFn = sn.p[k].maxFn;
      }
      cy.partHits.resize(sn.partHits);
      cy.hitEvents.resize(sn.hitEvents);
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
            const double prevGap = sn.gap[d.itemBase + k], vnPrev = sn.vn[d.itemBase + k];
            k++;
            if (!pr.overlap || prevGap <= kTouchM)
               continue;   // not touching, or it was touching at the start of the substep: not an arrival
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
         const double gapNow = (A.c - B.c).Length() - A.bs.R - B.bs.R;
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
         for (const int i : isl.members)
            ProbeBody(cy.bodies[i]);
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
         double h = (s == kMaxSubsteps - 1) ? cy.dt - t : NextStep(cy, isl, cy.dt - t, arrivals);
         for (const int i : isl.members)
         {
            cy.bodies[i].vS = cy.bodies[i].v;
            cy.bodies[i].hS = h;
         }
         if (hCap > 0.0)
         {
            h = std::min(h, hCap);
            hCap = 0.0;
         }
         sn.h = h;
         for (const int i : isl.members)
            SurfaceForces(cy, cy.bodies[i], h, s);
         PairForces(cy, isl, h, s);
         tripped = Integrate(cy, isl, h, s);
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
      if (ArrivalOn(cy) && IslandStatic(cy, isl))
         return SubstepsArrival(cy, isl);
      CPROF(Substeps);
      const double h = cy.dt / cy.steps;
      bool tripped = false;
      int s = 0;
      for (; s < cy.steps && !tripped; s++)
      {
         for (HitRubber *const b : cy.islandBands)
            b->BandSubstepBegin();
         for (const int i : isl.members)
            SurfaceForces(cy, cy.bodies[i], h, s);
         for (const int i : isl.members)
            BatForces(cy.bodies[i], h, s);
         for (const int i : isl.members)
            BandForces(cy.bodies[i]);
         PairForces(cy, isl, h, s);
         tripped = Integrate(cy, isl, h, s);
         for (HingePart::Part *const hg : cy.islandParts)   // the part moves with the reactions of all its balls summed, in its own substeps within h
            hg->StepInIsland(h);
         for (HitRubber *const b : cy.islandBands)   // the nodes step after the forces on them were summed
            tripped = b->BandSubstepEnd() || tripped;
         for (const int i : isl.members)
            if (cy.bodies[i].band != nullptr && cy.bodies[i].bb != nullptr && !std::isfinite(cy.bodies[i].band->BandBallPeak(*cy.bodies[i].bb)))
               tripped = true;
      }
      // The bat keeps stepping whatever happens to the ball: the rest of the cycle without a contact torque
      for (; s < cy.steps; s++)
      {
         StepBats(cy, isl, h, false);
         for (HingePart::Part *const hg : cy.islandParts)
            hg->StepInIsland(h);
      }
      for (HitRubber *const b : cy.islandBands)   // idle nodes merged away; a band that left the model's range is reset and its balls let go
         if (!b->BandCycleEnd())
            tripped = true;
      if (ContactTraceDetail())
         for (size_t k = 0; k < isl.edges.size(); k++)
         {
            const Body &A = cy.bodies[isl.edges[k].first], &B = cy.bodies[isl.edges[k].second];
            PLOGI << "ContactPair t=" << g_pplayer->m_time_msec << " balls " << isl.edges[k].first << ',' << isl.edges[k].second << " steps " << cy.steps << " max_delta_um " << cy.pst[k]->maxDelta * 1e6 << " max_fn " << cy.pst[k]->maxFn << " gap_end_um "
                  << ((A.c - B.c).Length() - A.bs.R - B.bs.R) * 1e6 << (tripped ? " TRIPPED" : "");
         }
      return tripped;
   }

   // The whole island is handed back untouched (nothing was written to the balls)
   static void HandBack(Cycle &cy, const ContactSystem::Islands::Island &isl)
   {
      for (const int i : isl.members)
      {
         Body &w = cy.bodies[i];
         if (w.tripped && w.o != nullptr)
         {
            PLOGW << "Rubber " << w.o->m_name << ": force-model contact left the model's range (" << w.tripWhy << " at substep " << w.tripStep << " of " << cy.steps << ", depth " << w.deepest * 1e3 << " mm, speed " << w.v.Length() << " m/s, cluster " << isl.members.size() << "), using the stock collision";
            w.o->m_contactCooldownUntilMs = g_pplayer->m_time_msec + 20;
         }
         if (w.tripped && w.o == nullptr && w.bat != nullptr)
         {
            PLOGW << "Rubber: " << w.bat->m_pflipper->GetName() << " contact step left the model's range (" << w.tripWhy << "), using the impulse path (ball untouched this cycle)";
            w.bat->m_rbCooldownUntilMs = g_pplayer->m_time_msec + 20;
         }
         else if (w.tripped && w.bat != nullptr)
            w.bat->m_rbCooldownUntilMs = g_pplayer->m_time_msec + 20;
         if (w.tripped && w.band != nullptr)
            w.band->m_bandCooldownUntilMs = g_pplayer->m_time_msec + 20;
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
         w.c = w.c + w.absCorr;
         pb->m_d.m_pos = Vertex3Ds((float)(w.c.x / U), (float)(w.c.y / U), (float)(w.c.z / U));
         pb->m_d.m_vel = Vertex3Ds((float)(w.v.x / VS), (float)(w.v.y / VS), (float)(w.v.z / VS));
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
         if (e.bat != nullptr)
            continue;   // the bat holds the ball until its own release rule lets go
         if (e.band != nullptr)   // the band holds the ball until its own release rule lets go (the ball stays in its other contacts)
         {
            for (BandBall &b : e.band->m_engaged)
               if (b.ball == pb && e.band->BandBallClear(b, w.bs, w.c, w.bandContact, w.lastGap, w.relNormalVel))
                  cy.bandReleases.emplace_back(pb, "clear");
            continue;
         }
         bool free = true, touching = false, closeBy = false;
         CPROF(ReleaseTest);
         for (const SurfaceContact &it : e.items)
         {
            const TubeProbe pr = it.r->Probe(w.bs, it, w.c);
            const double releaseGap = it.r->ReleaseGapM();
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
               if (gap <= releaseGap || !(w.v.Dot(pr.normal) > 0.0 || gap > 0.002))
                  free = false;
            }
         }
         if (free)
            cy.releases.emplace_back(pb, "clear");
         else if (w.o != nullptr && e.dwellS >= w.o->m_set.maxMs * 1e-3 && !touching && !closeBy)   // a ball in contact stays in the model however long (a resting ball, a pinch): re-engaging it would seat it and kick it
            cy.releases.emplace_back(pb, "timeout");
      }
   }

   // The bats that took part: finish the cycle (the pivot's nudge, then UpdateDisplacements only publishes the angle and fires the events)
   static void FinishBats(Cycle &cy, const ContactSystem::Islands::Island &isl)
   {
      for (const int i : isl.members)
      {
         Body &w = cy.bodies[i];
         if (w.bat == nullptr)
            continue;
         bool seen = false;
         for (const int j : isl.members)
            if (j == i)
               break;
            else if (cy.bodies[j].bat == w.bat)
               seen = true;
         if (seen)
            continue;   // (a bat with several balls on it finishes once)
         w.bat->NudgeStep(cy.dt);
         w.bat->m_gAdvanced = true;
         w.bat->m_gAdvHitStop = w.hitStop;
      }
   }

   // Releases, forgetting the pairs that no longer touch, and the deferred events (the script never runs inside the step)
   static void Finish(Cycle &cy)
   {
      CPROF(Finish);
      for (const auto &r : cy.releases)
         Release(r.first, r.second);
      for (const auto &r : cy.bandReleases)
         if (r.first->m_bandOwner != nullptr)
            r.first->m_bandOwner->BandRelease(r.first, r.second);
      for (HitRubber *const b : cy.bands)
         b->BandPublish();
      for (size_t i = 0; i < g_ballPairs.size();)
         if (!g_ballPairs[i].used)
            g_ballPairs.erase(g_ballPairs.begin() + (long)i);
         else
            i++;
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
      if (!g_deferred.empty())   // events raised inside the step, in the order they happened (a script may run another step's worth of code: the queue is taken first)
      {
         std::vector<std::function<void()>> q;
         q.swap(g_deferred);
         for (auto &f : q)
            f();
      }
   }
};

void ContactProfCensus(const ContactSurface *s, bool empty, size_t cand);

void ContactSystemImpl::StepBalls(const float dtime)
{
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
   cy.releases.clear();
   cy.hitEvents.clear();
   cy.partHits.clear();
   cy.hitFloor = std::min(0.5f, g_pplayer->m_physics->HitSpeedMin());
   cy.dtimeF = dtime;
   cy.bands.clear();
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
   static const bool timed = getenv("VPX_CONTACT_TIMING") != nullptr && *getenv("VPX_CONTACT_TIMING") != 0;
   uint64_t tp = timed ? usec() : 0;
   auto lap = [&](int &acc) { if (timed) { const uint64_t t = usec(); acc += (int)(t - tp); tp = t; } };
   GatherBalls(cy);
   GatherFree(cy);
   FindIslands(cy);
   lap(g_last.gatherUs);
   g_last.bodies = (int)cy.bodies.size();
   g_last.islands = (int)cy.islands.count;
   for (size_t ii = 0; ii < cy.islands.count; ii++)
   {
      const ContactSystem::Islands::Island &isl = cy.islands.list[ii];
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
      for (const int i : isl.members)
         g_last.bats += cy.bodies[i].bat != nullptr ? 1 : 0;
      const bool tripped = Substeps(cy, isl);
      lap(g_last.substepUs);
      FinishBats(cy, isl);
      if (tripped)
         HandBack(cy, isl);
      else
         WriteBack(cy, isl);
      lap(g_last.writeUs);
   }
   Finish(cy);
   lap(g_last.writeUs);
   if (ContactProf::On())
   {
      static int calls = 0;
      if (++calls % 5000 == 0)
         ContactProf::Dump("periodic");
   }
}


// ---- Proximity engagement ---------------------------------------------------------------------------------------------------------------------------------------------
// A free ball that can touch a surface this cycle is taken into the contact before the stock hit search runs: the surface is engaged when the ball is within its release gap, or
// closing on it faster than it can cross the gap in two cycles (reach 2 vn dt + the release gap; a ball released at the gap, receding, is not taken back). The first touch is
// then found by the substeps, and the part's hit event is queued from there (PartHit).
void ContactSystemImpl::EngageNear(const float dtime)
{
   CPROF(EngageNear);
   ContactProf::Cycle();
   static const bool on = []()
   {
      const char *const e = getenv("VPX_PROXIMITY");
      if (e != nullptr)
         return atoi(e) != 0;
      const std::filesystem::path ini = HitRubber::IniPathFor(g_pplayer->m_ptable->m_filename);
      return PhysicsIni::Active(ini).proximityEngage;
   }();
   if (!on)
      return;
   // The faces near each free ball, per surface, kept between cycles (a Verlet list: Refresh only rescans when the ball has moved out of the list): scanning a mesh playfield
   // afresh for every free ball every cycle would cost milliseconds.
   struct NearEntry
   {
      const HitBall *ball;
      const HitRubber *host;
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
      if (pb->m_contact != nullptr || pb->m_bandOwner != nullptr || pb->m_rubberOwner != nullptr || pb->m_d.m_lockedInKicker)
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
      HitRubber *best = nullptr;
      double bestGap = 1e9, pfGap = 1e9, pfClosing = 0.0;
      for (HitRubber *const o : HitRubber::Registry())
      {
         if (!o->m_enabled || !o->m_set.contact || o->m_set.band || g_pplayer->m_time_msec < o->m_contactCooldownUntilMs)
            continue;
         const FRect3D &bb = o->m_hitBBox;
         {
            const float dx = std::max({ bb.left - pos.x, 0.f, pos.x - bb.right }), dy = std::max({ bb.top - pos.y, 0.f, pos.y - bb.bottom }), dz = std::max({ bb.zlow - pos.z, 0.f, pos.z - bb.zhigh });
            minBoxGapU = std::min(minBoxGapU, std::sqrt(dx * dx + dy * dy + dz * dz) - pb->m_d.m_radius);
         }
         if (pos.x < bb.left - reachU || pos.x > bb.right + reachU || pos.y < bb.top - reachU || pos.y > bb.bottom + reachU || pos.z < bb.zlow - reachU || pos.z > bb.zhigh + reachU)
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
         o->Refresh(it, bs.R, travelM + 0.005, c);
         const TubeProbe pr = o->Probe(bs, it, c);
         if (pr.dist <= 0.0)
            continue;
         const double gap = pr.dist - (bs.R + pr.rho);
         const double rel = o->ReleaseGapM();
         const double closing = -v.Dot(pr.normal);
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
         EngageImpl(best, nullptr, pb, "near");
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
// A ball enters when a surface's stock collider finds it (the stock impulse is skipped) and leaves when it is clear of every surface and receding.
bool ContactSystemImpl::EngageImpl(HitRubber *const host, FlipperMoverObject *const bat, HitBall *const pb, const char *how)
{
   CPROF(EngageImpl);
   using namespace RubberContact;
   if ((host != nullptr && !host->m_set.contact) || pb->m_bandOwner != nullptr || pb->m_contact != nullptr || pb->m_rubberOwner != nullptr || pb->m_d.m_lockedInKicker)
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
      for (HitRubber *const o : HitRubber::Registry())
      {
         if (o != host)
         {
            if (!o->m_enabled || !o->m_set.contact)
               continue;
            if (o->IsPlayfield() && InPlayfieldHole(posStart))
               continue;   // no floor over a kicker's hole
            if (host == nullptr && bat == nullptr && o->m_set.band)
               continue;   // a band engagement: the band rubbers are the band model's; as a rigid surface contact they would stop the ball before the band could deflect
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
   if (bat == nullptr && host != nullptr)   // (a bat engagement seats its start against the moving bat in the first step: SetupBat)
   {
   constexpr double kMaxRewindS = PHYSICS_STEPTIME_S;
   const double zBeforeSeat = c.z;
   for (SurfaceContact &it : rb.items)
      it.r->Refresh(it, bs.R, 0.03, c);   // a surface: the faces within 30 mm of the ball
   for (int pass = 0; pass < 3; pass++)
      for (const SurfaceContact &it : rb.items)
      {
         const TubeProbe p0 = it.r->Probe(bs, it, c);
         if (!p0.overlap || p0.delta <= 1e-6)
            continue;
         rb.delta0 = std::max(rb.delta0, p0.delta);
         if (ContactTraceOn())
            PLOGI << "ContactSeat " << it.r->Name() << " pass " << pass << " overlap_um " << p0.delta * 1e6 << " normal " << p0.normal.x << ' ' << p0.normal.y << ' ' << p0.normal.z << " ball " << c.x / U << ' ' << c.y / U << ' ' << c.z / U;
         const double closing = -v.Dot(p0.normal);
         if (closing > 0.01)
         {
            const double t = std::min(p0.delta / closing, kMaxRewindS);
            c = c - v * t;
            rb.rewound = std::max(rb.rewound, t);
            const TubeProbe p1 = it.r->Probe(bs, it, c);
            if (p1.overlap && p1.delta > 1e-6)
               c = c + p1.normal * p1.delta;
         }
         else
            c = c + p0.normal * p0.delta;
      }
   // A ball deep inside a rubber (a hit the stock collision missed, or a script placing it there) can be seated along a normal that points down. With the playfield in the
   // contact that would put it under the floor, beyond what the stock plane can still catch (more than two radii below it): it would fall out of the table. Do not engage then;
   // the stock collision resolves the overlap the way it always did.
   if (c.z < zBeforeSeat - 1.0 * U && std::any_of(rb.items.begin(), rb.items.end(), [](const SurfaceContact &it) { return it.r->IsFloor(); }))
   {
      PLOGW << "Rubber " << host->m_name << ": the seat of a " << rb.delta0 * 1e3 << " mm overlap would put the ball under the playfield, not engaging (stock collision)";
      return false;
   }
   }
   pb->m_d.m_pos = Vertex3Ds((float)(c.x / U), (float)(c.y / U), (float)(c.z / U));
   rb.vin = v;
   const double inertiaU = (double)pb->Inertia();
   const double wScale = 1.0 / (inertiaU * DEFAULT_STEPTIME_S);
   rb.win = V3((double)pb->m_angularmomentum.x * wScale, (double)pb->m_angularmomentum.y * wScale, (double)pb->m_angularmomentum.z * wScale);
   rb.posStart = posStart;
   if (ContactTraceOn())
      PLOGI << "ContactEngage " << (host != nullptr ? host->m_name : std::string("bat")) << " t=" << g_pplayer->m_time_msec << " how " << how << " speed " << v.Length() << " vn_in " << rb.vin.Length() << " delta0_um " << rb.delta0 * 1e6
            << " rewind_us " << rb.rewound * 1e6 << " flip_clear " << g_pplayer->m_physics->FlipperClearance(pb->m_d, kFlipperHorizonVpt) << g_pplayer->m_physics->FlipperDebug(pb->m_d, kFlipperHorizonVpt) << " rubbers " << rb.items.size() << " pos " << pb->m_d.m_pos.x << ' ' << pb->m_d.m_pos.y << ' ' << pb->m_d.m_pos.z;
   rb.host = host;
   rb.bat = bat;
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
   DetachBat(*itE);
   if (itE->band != nullptr && ball->m_bandOwner == itE->band)
      itE->band->BandRelease(ball, why);   // clears itE->band and the ball's band owner
   itE->band = nullptr;
   const BallContact e = *itE;
   HitRubber *const host = e.host;
   if (getenv("VPX_RUBBERS_LOG"))
   {
      using namespace RubberContact;
      const double inertiaU = (double)ball->Inertia();
      const double wScale = 1.0 / (inertiaU * DEFAULT_STEPTIME_S);
      const V3 v((double)ball->m_d.m_vel.x * VS, (double)ball->m_d.m_vel.y * VS, (double)ball->m_d.m_vel.z * VS);
      const V3 w((double)ball->m_angularmomentum.x * wScale, (double)ball->m_angularmomentum.y * wScale, (double)ball->m_angularmomentum.z * wScale);
      PLOGI << "RubberSolve " << (host != nullptr ? host->m_name : std::string("bat")) << " t=" << g_pplayer->m_time_msec << " steps " << e.steps << " ms " << e.dwellS * 1e3 << " peakFn " << e.peakFn << " why " << why
            << " vin " << e.vin.x << ' ' << e.vin.y << ' ' << e.vin.z << " vout " << v.x << ' ' << v.y << ' ' << v.z
            << " win " << e.win.x << ' ' << e.win.y << ' ' << e.win.z << " wout " << w.x << ' ' << w.y << ' ' << w.z
            << " delta0 " << e.delta0 * 1e3 << " rewind_ms " << e.rewound * 1e3 << " rubbers " << e.items.size()
            << " pos " << e.posStart.x << ' ' << e.posStart.y << ' ' << e.posStart.z << " (m/s, rad/s, VPX units)";
   }
   if (ContactTraceOn())
      PLOGI << "ContactRelease " << (host != nullptr ? host->m_name : std::string("bat")) << " t=" << g_pplayer->m_time_msec << " why " << why << " dwell_ms " << e.dwellS * 1e3 << " peakFn " << e.peakFn
            << " speed " << std::sqrt((double)ball->m_d.m_vel.Dot(ball->m_d.m_vel)) * 0.53975e-3 / DEFAULT_STEPTIME_S << " pos " << ball->m_d.m_pos.x << ' ' << ball->m_d.m_pos.y << ' ' << ball->m_d.m_pos.z;
   if (ContactProf::On())
      ContactProf::Note(std::string("release ") + why);
   ball->m_contact = nullptr;
   ball->m_nearSafeU = -1.f;   // scan again at the next cycle
   g_engaged.erase(itE);
}

// The bat lets go of the ball (clear and receding, or the contact is dropped): the ball stays in whatever other contact it has
void ContactSystemImpl::DetachBat(BallContact &e)
{
   FlipperMoverObject *const m = e.bat;
   if (m == nullptr)
      return;
   if (e.ball != nullptr)
      PLOGI << "Rubber: " << m->m_pflipper->GetName() << " released ball at t=" << g_pplayer->m_time_msec << " ms, vel (" << e.ball->m_d.m_vel.x << ',' << e.ball->m_d.m_vel.y << ',' << e.ball->m_d.m_vel.z << ") U/T";
   if (e.ball != nullptr && e.ball->m_rubberOwner == m)
      e.ball->m_rubberOwner = nullptr;
   if (m->m_rbBall == e.ball)
      m->m_rbBall = nullptr;
   m->m_rbClearSteps = 0;
   e.bat = nullptr;
}

bool ContactSystem::EngageBand(HitRubber *const band, HitBall *const pb)
{
   BallContact *e = pb->m_contact;
   if (e != nullptr && e->band != nullptr && e->band != band)
      return false;   // another band has the ball
   if (e == nullptr)
   {
      if (!ContactSystemImpl::EngageImpl(nullptr, nullptr, pb, "band"))
         return false;
      e = pb->m_contact;
   }
   e->band = band;
   return true;
}

bool ContactSystem::EngageBat(FlipperMoverObject *const bat, HitBall *const pb)
{
   BallContact *e = pb->m_contact;
   if (e != nullptr && e->bat != nullptr && e->bat != bat)
      return false;   // another bat has the ball
   if (e == nullptr)
   {
      if (!ContactSystemImpl::EngageImpl(nullptr, bat, pb, "bat"))
         return false;
      e = pb->m_contact;
   }
   e->bat = bat;
   e->batSt.Reset();
   e->batLastDist = 0.0;
   e->batClear = 0;
   e->batJustEngaged = true;
   bat->m_rbBall = pb;
   bat->m_rbClearSteps = 0;
   pb->m_rubberOwner = bat;
   PLOGI << "Rubber: " << bat->m_pflipper->GetName() << " engaged ball at t=" << g_pplayer->m_time_msec << " ms, pos (" << pb->m_d.m_pos.x << ',' << pb->m_d.m_pos.y << ',' << pb->m_d.m_pos.z << ") vel (" << pb->m_d.m_vel.x << ',' << pb->m_d.m_vel.y << ',' << pb->m_d.m_vel.z << ") U/T";
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
   for (HitRubber *const o : HitRubber::Registry())
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

bool ContactSystem::Engage(HitRubber *const host, HitBall *const pb, const char *how)
{
   return ContactSystemImpl::EngageImpl(host, nullptr, pb, how);
}

void ContactSystem::Release(HitBall *const ball, const char *why)
{
   ContactSystemImpl::Release(ball, why);
}

void ContactSystem::ReleaseHost(HitRubber *const host)
{
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

void ContactSystem::ReleaseBat(FlipperMoverObject *const bat)
{
   // the bat is going away: every ball on it lets go without a log line. With several balls on one bat the bat's m_rbBall names only the last one, so the
   // others are found here (a ball can outlive the physics engine: a pointer left to the bat crashed ~HitBall at table close, 2026-10-20)
   for (BallContact &e : g_engaged)
      if (e.bat == bat)
      {
         if (e.ball != nullptr && e.ball->m_rubberOwner == bat)
            e.ball->m_rubberOwner = nullptr;
         e.bat = nullptr;
      }
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
   std::vector<std::pair<uint64_t, std::string>> top;
   for (const auto &kv : g_census)
      top.emplace_back(kv.second.itemCycles, kv.second.name + " (no face in reach " + std::to_string((int)(100.0 * kv.second.emptyCycles / std::max<uint64_t>(1, kv.second.itemCycles))) + "%, avg faces on list " + std::to_string((double)kv.second.candSum / std::max<uint64_t>(1, kv.second.itemCycles)) + ")");
   std::sort(top.rbegin(), top.rend());
   for (size_t i = 0; i < top.size() && i < 14; i++)
      PLOGI << "PROF surface " << top[i].second << " body-cycles " << top[i].first;
}

void ContactSystem::Step(const float dtime)
{
   ContactSystemImpl::StepBalls(dtime);
}

const ContactSystem::CycleStats &ContactSystem::LastCycle()
{
   return g_last;
}
