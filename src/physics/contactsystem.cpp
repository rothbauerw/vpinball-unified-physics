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
static double g_subS = 0.0;         // the engine's current sub-interval in seconds, as ContactSystem::Step was given it (the parts and cords advanced alone after it read it)
static std::vector<const HitBall *> g_placed;   // the balls a script or the driver put where they are since the last cycle (NotePlaced)
static ContactSystem::CycleStats g_last;
static std::vector<std::function<void()>> g_deferred;      // events raised inside the step, fired after it (Finish)
// The test driver's cost split of this cycle (docs/unified-physics/TestMatrixDesign.txt section 3 point 10 b-d): null unless the driver times the cycle in detail. Per island: what set each
// substep (TestDriver::ContactCost::Why) and the smallest substep.
void ContactSystem::DeferEvent(std::function<void()> f) { g_deferred.push_back(std::move(f)); }
static constexpr float kFlipperHorizonVpt = 1.0f;   // trace only: VPX time units (1 VPT = 10 ms) a bat is looked ahead for in the flip_clear numbers

struct ContactSystemImpl
{
   static void StepBalls(double dt);   // dt: seconds
   static void Release(HitBall *ball, const char *why);
   static void Regather(BallContact &e, const Vertex3Ds &pos, float radiusU, bool finder);   // finder: shared_finder P2b, the static mesh hosts come from the ball's query
   static bool EngageImpl(ContactHost *host, HitBall *pb, const char *how);
   static void EngageNear(double dt);   // dt: the cycle in seconds
   static void CordNear(double dt);     // likewise

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
      bool magnetsOn = false;
      bool extTaken = false;    // Sources::Begin took this cycle's gravity and nudge out of the velocity (they act in the substeps)
      V3 extOwe;                // the step's gravity and nudge for the sub-intervals after this one: out of the velocity while this one is stepped, back at the write-back
      double endPair = 0.0;     // the energy check: its half of what its ball pairs hold at the island's end
      double magAcc = 0.0;      // m/s2: the magnets' acceleration of this ball at the start of the step (the substep bound)
      double eddyRate = -1.0;   // eddy-current braking of the spin across the coil field, 1/s: evaluated once per cycle (negative = not yet); eddyDir = the field direction
      V3 eddyDir;
      double eStart = 0.0, eScale = 0.0;   // the energy check: the ball's energy at the cycle's start, and what the allowance is a share of
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
      int steps = 1;                     // the current island's substeps
      std::vector<std::pair<int, int>> links;
      std::vector<std::pair<HitBall *, const char *>> bandReleases;   // balls the band lets go of (after the write-back: the release can fire the script)
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

   // The balls held by a force-model contact
   static void GatherBalls(Cycle &cy)
   {
      CPROF(GatherBalls);
      HitRubber::FinderBuild();   // shared_finder: the index over every static mesh, once before the first cycle (P1: only check mode reads it)
      HitRubber::FinderTotals(false);   // shared_finder P2a: counts the cycle; the query's totals once a second of table time
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
            const bool finder = HitRubber::FinderDrives();   // shared_finder P2b (key on): the cycle's query finds the ball's static mesh surfaces
            if (finder)
            {
               LoadState(w);
               HitRubber::FinderQuery(pb, w.c, w.v, w.bs.R, dt);   // the per-cycle query (design.md section 3), before the gather that reads it
            }
            {
               CPROF(Regather);
               Regather(e, w.posStart, pb->m_d.m_radius, finder);   // (over a kicker's hole the playfield is left out)
            }
            if (PlayfieldHolesExist() && InPlayfieldHole(pb->m_d.m_pos) && std::any_of(e.items.begin(), e.items.end(), [](const SurfaceContact &it) { return it.r->IsPlayfield(); }))
            {
               cy.releases.emplace_back(pb, "kicker hole");   // rolled into a non-legacy kicker's circle with nothing else to hold it
               continue;
            }
            if (!finder)
               LoadState(w);
            HitRubber::FinderBind(pb, e.items, true);   // shared_finder (key on): each static mesh item reads its faces from the ball's query (its Refresh below does nothing)
            w.band = e.band;
            Sources::Begin(w, dt);
            {
               CPROF(RefreshCycle);
               for (SurfaceContact &it : e.items)
                  it.r->Refresh(it, w.bs.R, 2.0 * w.v.Length() * dt + 0.005, w.c);   // a surface: the faces the ball can reach this cycle
            }
            HitRubber::FinderCheckGather(pb, e.items, w.c, w.bs.R, PlayfieldHolesExist() && InPlayfieldHole(pb->m_d.m_pos));   // check mode only (VPX_FINDER_CHECK=1)
            cy.bodies.push_back(std::move(w));
         }
      }
   }

   // A free ball (not held by any contact) that can reach another ball, held or free, this cycle joins the island with no surface contacts of its own: the engine has given it
   // gravity already, and the stock collision is skipped for the pair, so the two are solved together here. Chains grow until no free ball is within reach of a member.
   // [Features] free_flight (physicsini.h): every free ball is a body of the step, not only one within reach of another ball, so no ball is moved by the stock flight
   static bool FreeFlightOn()
   {
      static const bool on = []() { return PhysicsIni::Active(ContactHost::IniPathFor(g_pplayer->m_ptable->m_filename)).freeFlight; }();
      return on;
   }

   static void GatherFree(Cycle &cy)
   {
      CPROF(GatherFree);
      const ContactSettings &cs = *cy.cs;
      const bool every = FreeFlightOn();
      if (!cs.ballContact && !every)
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
         cand.push_back(std::move(w));
      }
      if (every)   // free_flight: each one a body with no contacts (its pair contacts, if any, are found as for any body: FindIslands)
      {
         for (Body &w : cand)
         {
            Sources::Begin(w, cy.dt);   // (its gravity, nudge and magnet pull: taken back, then given in the substeps)
            cy.bodies.push_back(std::move(w));
         }
         return;
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
            vc[find(i)] += (ballPairs ? 0.0 : w.v.Length()) + (w.extBaseVS.Length() / PHYSICS_STEPTIME_S + w.magAcc) * T;
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
            // off the cord and clear of it; a plate on its way out can still catch a ball close by, so it stays in the group until the plate is back at its rest angle and still, or the ball is far
            // (band_plate_reach: beyond the plate's stroke at the cord, CordFixed::PlateReachM, and this rule's 1 mm; otherwise 30 mm)
            CordBand &cb = *e.band->Cord();
            int slot = -1;
            for (int k = 0; k < (int)cb.slots.size(); k++)
               if (cb.slots[k] == pb)
                  slot = k;
            const double gap = slot >= 0 ? cb.Gap(slot) : 1.0;
            const double clearM = e.band->m_set.bandPlateReach ? 0.001 + cb.fixed->PlateReachM() : 0.03;
            if (gap > clearM || (gap > 0.001 && cb.PlateAtRest()))
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
   // ---- An island with a cord: one group step for the whole island (docs/unified-physics/GroupStepDesign.txt 8, 11) ----
   // A surface of a ball's contact as a shape of the group (its own list of near faces, so one ball's)
   class SurfaceShape : public GroupStep::Shape
   {
   public:
      const ContactSurface *r = nullptr;
      const SurfaceContact *it = nullptr;
      BallSI bs;
      double heldReff = 0.0;   // > 0: the radius this touch started with (the group takes a shape's law anew each cycle; a polygon's radius jumps at a corner, the energy must not)
      double heldKb = 0.0;     // the bend this touch started with, for the same reason
      double heldKtb = 0.0;    // and its bend across the normal
      ContactMovingPart *mp = nullptr;   // a part on the group: its own surface, which moves with it
      int surface = 0;
      std::string Name() const override { return r != nullptr ? r->Name() : std::string(); }
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
               h.ktb = heldKtb;
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
         acc = std::max(acc, w.extBaseVS.Length() / PHYSICS_STEPTIME_S + w.magAcc);
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
      g.BeginUse();
      g.ImpClear();   // (test tooling: a use whose result was not kept leaves nothing)
      static std::map<int, std::string> impShapeNames, impPartNames;   // (test tooling: names of the shapes and parts with no Name() of their own, by index in this use)
      impShapeNames.clear();
      impPartNames.clear();
      static std::deque<PullField> pulls;
      pulls.clear();
      for (const Magnets::CellPull &cp : Magnets::CellPulls())
      {
         pulls.emplace_back();
         pulls.back().p = cp;
         g.AddField(&pulls.back());
      }
      cy.cs->ToStep(g.set);   // the engine's step settings and floor, as a part alone gets them (GroupStep::EngineSettings)
      g.set.pair = cy.rpBall;
      if (ring)   // (cy.rpBall holds the old loop's halved modulus, and the group halves it for a pair itself; its rolling resistance is on the pair's radius, R / 2)
      {
         g.set.pair.eMpa = cy.cs->ballEMpa;
         g.set.pair.rollingMu = 0.5 * cy.cs->ballRollingResistance;
      }
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
                  sh.heldKtb = mp != nullptr && it.Surf(sf).touching ? it.Surf(sf).ktb : 0.0;
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
                        if (GroupStep::RestWatch *const rw = mp->Rest())   // the one restart of the rest rule: a ball touches the part this cycle, so its watch starts again from its pose now (GroupStep::RestWatch)
                           rw->Restart(mp->RestPose());
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
         c->Cord()->fixed->plateRest.Reset();   // the same restart for the cord's plate (it joins the group here, not above): a ball touches it, its watch alone starts again when it is alone
         c->Cord()->fixed->plateIdle = false;
         const int impS0 = g.ns;   // (test tooling: the cord's shapes' names below)
         c->Cord()->fixed->AddPostShapes(g, c->Cord()->fixedTube);
         if (GroupStep::Group::ImpOn())   // test tooling (VPX_TEST_IMPULSES=1): the cord's posts and plate, and its own force on a ball, by its rubber's name
         {
            for (int s = impS0; s < g.ns; s++)
               impShapeNames[s] = c->Name() + (g.shape[s] == &c->Cord()->fixedTube ? std::string(".plate") : ".post" + std::to_string(s - impS0));
            impPartNames[g.np - 1] = c->Name() + ".cord";
         }
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
      // the cycle's end (a failed step: those of the part of the cycle that was stepped).
      static std::vector<char> hitOn;
      hitOn.resize(g.ni);
      static std::vector<PartHit> ringParts;
      static std::vector<BallHit> ringPairs;
      ringParts.clear();
      ringPairs.clear();
      if (ring || hasPart)
      {
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
         // rule 8: a tripped net shows on the playfield; every energy step, not once per table start (user, 2026-11-11: a break, to be seen each time it happens, for now)
         if (g_pplayer != nullptr && g_pplayer->m_liveUI != nullptr)
         {
            std::string note = "Energy step (a fault): parts [";
            for (size_t a = 0; a < names.size(); a++)
               note += (a > 0 ? ", " : "") + names[a];
            note += "] at " + std::to_string(g_pplayer->m_time_msec) + " ms: " + how + " (see vpinball.log, 'Group step: energy step')";
            g_pplayer->m_liveUI->PushNotification(note, 6000);
         }
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
            fprintf(o, "t_ms,done_us,h_us,hcap_us,rule,partRule,item,itemPart,touching,b,delta_um,dDot,omN,prevRounds,prevUs,envSolves,envIters,capsules,envHits,envMisses,envUs,unsteady,unsteadyVal,unsteadyLim,missProbe,missEnd,missMid,missStart,missRate,missAngle,missEnergy,term\n");
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
            fprintf(ruleF, "%u,%.3f,%.4f,%.4f,%d,%d,%d,%d,%d,%d,%.4f,%.5f,%.1f,%ld,%.1f,%ld,%ld,%ld,%ld,%ld,%.1f,%d,%.6g,%.6g,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%d\n", g_pplayer->m_time_msec, done * 1e6, h * 1e6, hc * 1e6, dg.rule, dg.partRule, dg.item,
                    hasItem && g.item[dg.item].b < 0 ? g.shape[g.item[dg.item].s]->part : -9, hasItem ? (int)g.item[dg.item].touching : -1, hasItem ? g.item[dg.item].b : -9, dg.delta * 1e6, dg.dDot, dg.omN,
                    g.books.rounds - ruleRounds, std::chrono::duration<double, std::micro>(now - ruleT).count(), es.solves, es.iters, es.capsules, es.hits, es.misses, es.us, dg.unsteady, dg.unsteadyVal, dg.unsteadyLim,
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
               const long *const rc = GroupStep::RestWatch::Counts();   // (parts alone: bumper, hinge part, disc, flipper, sling plate)
               fprintf(f, "alone calls %ld %ld %ld %ld %ld asleep after %ld %ld %ld %ld %ld\n", rc[0], rc[1], rc[2], rc[3], rc[4], rc[5], rc[6], rc[7], rc[8], rc[9]);
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
      if (ring && !ok && !hasPart)   // a net, not the fix: on screen once (rule 8). As an island with a moving part (below): the balls keep the state reached and the rest of the
      {                              // cycle is lost, so no energy can appear; nothing retakes it (the old island loop is no longer reached)
         static uint32_t lastMs = 0;
         static long failed = 0;
         static const void *warned = nullptr;
         failed++;
         if (g_pplayer->m_time_msec - lastMs > 500)
         {
            lastMs = g_pplayer->m_time_msec;
            PLOGW << "Ring group failed (" << failed << " in all): why " << g.failWhy << " at " << done * 1e6 << " us of " << cy.dt * 1e6 << " us, balls " << g.nb << ", ball 0 at (" << g.ball[0].c.x / U << ", " << g.ball[0].c.y / U << ", " << g.ball[0].c.z / U << ") v " << g.ball[0].v.Length() << " m/s: the rest of the cycle is lost";
         }
         if (warned != g_pplayer)
         {
            warned = g_pplayer;
            const std::string msg = "Group step: a ball's step could not be solved: the rest of that cycle was lost (see vpinball.log)";
            PLOGW << msg;
            g_pplayer->m_liveUI->PushNotification(msg, 6000);
         }
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
      if (GroupStep::Group::ImpOn())   // test tooling (VPX_TEST_IMPULSES=1): the cycle's impulses by ball and source to the test driver, now that its result is kept
      {
         for (const auto &e : g.imp)
         {
            const int i = e.first.first;
            if (i < 0 || i >= g.nb)
               continue;
            std::string src = e.first.second;
            if (src.size() > 1 && src[0] == '#')   // a ball pair: "ball" and the other ball's driver name
            {
               const int j = std::atoi(src.c_str() + 1);
               const char *const nj = j >= 0 && j < g.nb ? TestDriver::BallName(balls[j]) : nullptr;
               src = std::string("ball") + (nj != nullptr ? nj : src.c_str() + 1);
            }
            else if ((src.compare(0, 5, "shape") == 0 || src.compare(0, 4, "part") == 0) && src.find_first_not_of("0123456789", src[0] == 's' ? 5 : 4) == std::string::npos
                     && src.size() > (src[0] == 's' ? 5u : 4u))   // the group's own "shape<s>" / "part<p>" (a cord's post, plate or own force): its rubber's name
            {
               const bool isShape = src[0] == 's';
               const std::map<int, std::string> &nm = isShape ? impShapeNames : impPartNames;
               const auto f = nm.find(std::atoi(src.c_str() + (isShape ? 5 : 4)));
               if (f != nm.end())
                  src = f->second;
            }
            const double J[3] = { e.second.J.x, e.second.J.y, e.second.J.z }, L[3] = { e.second.L.x, e.second.L.y, e.second.L.z };
            TestDriver::AddImpulse(balls[i], src, J, L);
         }
         g.ImpClear();
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
            kept.ktb = gi.touching ? gi.ktb : 0.0;
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

   // a cord with no ball this cycle: a plate moving by itself, or the path letting go of a ball just gone (GroupStepDesign 8 point 5); dt: seconds
   static void CordAlone(HitRubber *const band, const double dt)
   {
      CordBand &cb = *band->Cord();
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
      const bool was = cb.PlateQuiet();
      if (!GroupStep::AdvanceAlone(*cb.fixed, dt, (int)g_pplayer->m_time_msec))
         band->CordSolveFailed("the plate's step could not be solved");
      cb.fixed->PlateRestStep(dt);
      GroupStep::RestWatch::Book(4, was, cb.PlateQuiet());
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
         if (const PhysicsEngine &pe = *g_pplayer->m_physics; pe.HitSpeedEvents())
         {
            ph.pb->m_hitSpeed = ph.speed;   // every surface's first touch (shared_finder P9: a converted part's stock copy no longer sets it), before the part's Hit event below
            if (ph.s->IsPlayfield() && ph.speed >= pe.HitSpeedMin())
               g_pplayer->m_ptable->FireVoidEvent(DISPID_HitEvents_Hit);   // the playfield has no Hit event of its own: the table's
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
      for (ContactMovingPart *const mp : g.parts)
         mp->restTouched = false;
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
         mp->restTouched = true;   // (a ball asleep on it: the part does not sleep, and the pile stays asleep while the part does not move)
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
   static void CheckOverlaps(float hittime, double hS);   // hittime: engine time units (the free flight only); hS: the same in seconds
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

// dt: the time to step, seconds
void ContactSystemImpl::StepBalls(const double dt)
{
   if (!s_append)
      DispatchPending();   // (the engine dispatches after every step; this keeps a step that was not followed by one from losing its events)
   g_last = ContactSystem::CycleStats();
   if (g_engaged.empty() && !(ContactSettings::Get().ballContact && g_pplayer->m_vball.size() >= 2) && !(FreeFlightOn() && !g_pplayer->m_vball.empty()))   // (two free balls are solved here too; with free_flight every ball)
      return;
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
      if (IslandHasCord(cy, isl))
      {
         CordIsland(cy, isl);
         lap(g_last.substepUs);
         EnergyCheckEnd(cy, isl);
         continue;
      }
      // Every other island goes to the group step (a moving part, a rubber in reach, a held ball with nothing in reach, free balls near each other): one contact path for
      // every ball. A step it cannot solve loses the rest of that cycle (CordIsland), as with a moving part; there is no other path to retake it on.
      CordIsland(cy, isl, true);
      lap(g_last.substepUs);
      EnergyCheckEnd(cy, isl);
      SleepConsider(cy, isl);
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
// Every ball, one already held for the playfield included: its contact is handed over to the band. dt: the cycle in seconds (CordNear and EngageNear).
void ContactSystemImpl::CordNear(const double dt)
{
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
         // a plate on its way out carries the cord this far past its rest line (as the release rule): band_plate_reach, the plate's stroke at the cord (CordFixed::PlateReachM); otherwise 30 mm
         const bool bound = band->m_set.bandPlateReach;
         const double plateReach = band->CordOn() && !band->Cord()->PlateAtRest() ? (bound ? band->Cord()->fixed->PlateReachM() : 0.03) : 0.0;
         // a real section thicker than the table's: the cord's surface stands out from the table's tube by twice the difference (its centre moves out by it, and its radius grows by it); the probe below reads the table's tube
         const double fixedRc = band->CordOn() ? band->Cord()->fixed->rc : 0.0;
         // band_plate_reach: as far past the band's box (its points and the table's tube radius) as the engage rule below can take a ball: the real section standing out, the release gap, and the farther of the plate's stroke and two cycles' closing
         const float reachU = bound ? pb->m_d.m_radius + (float)((std::max(2.0 * (fixedRc - band->m_radiusSI), 0.0) + o->ReleaseGapM() + std::max(plateReach, 2.0 * speed * dt)) / U)
                                    : pb->m_d.m_radius + (float)((2.0 * speed * dt + plateReach + (fixedRc > 0.0 ? 2.0 * fixedRc : 0.0)) / U) + 4.0f;
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

void ContactSystemImpl::EngageNear(const double dt)
{
   CordNear(dt);
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
      // shared_finder P2b (key on): the static mesh hosts are tested every cycle from the ball's query (radius R + the reach below), the Verlet skip covers only the others
      const bool finder = HitRubber::FinderDrives();
      const bool skipOthers = pb->m_nearSafeU >= 0.f && nowMs - pb->m_nearMs < 10 && std::sqrt((double)(pos - pb->m_nearPos).LengthSquared()) + travelM / U < (double)pb->m_nearSafeU;
      if (skipOthers && !finder)
         continue;
      float minBoxGapU = 1e9f;
      BallSI bs;
      bs.R = (double)pb->m_d.m_radius * U;
      const V3 c((double)pos.x * U, (double)pos.y * U, (double)pos.z * U);
      const V3 v((double)pb->m_d.m_vel.x * VS, (double)pb->m_d.m_vel.y * VS, (double)pb->m_d.m_vel.z * VS);
      ContactHost *best = nullptr;
      double bestGap = 1e9, pfGap = 1e9, pfClosing = 0.0;
      thread_local std::vector<ContactHost *> hosts;
      hosts.clear();
      if (finder)
      {
         HitRubber::FinderQueryAt(pb, c, bs.R, travelM + 0.005);   // as the Refresh below asked of each mesh: the faces within R + 2 x speed x 1 ms + 5 mm (a static surface: no bound of its own)
         HitRubber::FinderHosts(pb, hosts);
      }
      if (!skipOthers)
         for (ContactHost *const o : ContactHost::Registry())
            if (!finder || o->FinderObjIndex() < 0)
               hosts.push_back(o);
      if (finder)   // the registry's order, as the one loop gave them (the nearest wins; a tie goes to the first)
         std::stable_sort(hosts.begin(), hosts.end(), [](const ContactHost *a, const ContactHost *b) { return a->m_regOrder < b->m_regOrder; });
      for (ContactHost *const o : hosts)
      {
         if (!o->m_enabled || !o->m_set.contact || o->m_set.band)
            continue;
         const bool mesh = finder && o->FinderObjIndex() >= 0;
         const double vb = PreventOn() ? o->SurfaceSpeedBoundMS(dt) : 0.0;   // the surface's own closing-speed bound (a moving part; 0 for a static surface)
         const float reachHost = reachU + (float)(2.0 * vb * dt / U);
         const FRect3D &bb = o->m_hitBBox;
         if (!mesh)   // (a mesh host found by the query: it has a face within the reach; the Verlet skip's bound is the others')
         {
            const float dx = std::max({ bb.left - pos.x, 0.f, pos.x - bb.right }), dy = std::max({ bb.top - pos.y, 0.f, pos.y - bb.bottom }), dz = std::max({ bb.zlow - pos.z, 0.f, pos.z - bb.zhigh });
            minBoxGapU = std::min(minBoxGapU, std::sqrt(dx * dx + dy * dy + dz * dz) - pb->m_d.m_radius);
         }
         if (!mesh && (pos.x < bb.left - reachHost || pos.x > bb.right + reachHost || pos.y < bb.top - reachHost || pos.y > bb.bottom + reachHost || pos.z < bb.zlow - reachHost || pos.z > bb.zhigh + reachHost))
            continue;
         if (o->m_surfaceKind == "playfield" && InPlayfieldHole(pos))
            continue;   // the kicker's
         o->ReloadIfChanged();
         SurfaceContact itq;   // a mesh host found by the query: its faces from the query (no cached list)
         if (mesh)
         {
            itq.r = o;
            HitRubber::FinderBindOne(pb, itq);
         }
         auto ce = cache.end();
         if (!mesh)
         {
         ce = std::find_if(cache.begin(), cache.end(), [&](const NearEntry &n) { return n.ball == pb && n.host == o; });
         if (ce == cache.end())
         {
            cache.push_back({ pb, o, SurfaceContact(), nowMs });
            ce = cache.end() - 1;
            ce->it.r = o;
         }
         ce->seenMs = nowMs;
         }
         SurfaceContact &it = mesh ? itq : ce->it;
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
         if (!skipOthers)   // (key on with the skip: the others were not scanned, their bound stands)
         {
         pb->m_nearPos = pos;
         pb->m_nearMs = nowMs;
         pb->m_nearSafeU = minBoxGapU - 6.0f;
         }   // (a surface within its release gap, about 4 U, qualifies; the margin keeps the box test conservative)
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
      // shared_finder P2b (key on): the static mesh hosts other than the one that found the ball come from a query at the engagement, radius R + 20 U (the reach above, to a face
      // instead of the host's box); kept in the registry's order
      const bool finder = HitRubber::FinderDrives();
      thread_local std::vector<ContactHost *> hosts;
      hosts.clear();
      if (finder)
      {
         HitRubber::FinderQueryAt(pb, c, bs.R, 20.0 * U);
         HitRubber::FinderHosts(pb, hosts);
      }
      for (ContactHost *const o : ContactHost::Registry())
         if (!finder || o->FinderObjIndex() < 0 || o == host)
            hosts.push_back(o);
      if (finder)
      {
         std::stable_sort(hosts.begin(), hosts.end(), [](const ContactHost *a, const ContactHost *b) { return a->m_regOrder < b->m_regOrder; });
         hosts.erase(std::unique(hosts.begin(), hosts.end()), hosts.end());   // (the host found by both)
      }
      for (ContactHost *const o : hosts)
      {
         const bool mesh = finder && o->FinderObjIndex() >= 0;
         if (o != host)
         {
            if (!o->m_enabled || !o->m_set.contact)
               continue;
            if (o->IsPlayfield() && InPlayfieldHole(posStart))
               continue;   // no floor over a kicker's hole
            if (o->m_set.band)
               continue;   // the band rubbers are the band model's: as a rigid surface contact they would stop the ball before the band could deflect (as Regather)
            const FRect3D &bb = o->m_hitBBox;
            if (!mesh && (posStart.x < bb.left - reach || posStart.x > bb.right + reach || posStart.y < bb.top - reach || posStart.y > bb.bottom + reach
                || posStart.z < bb.zlow - reach || posStart.z > bb.zhigh + reach))
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
      HitRubber::FinderBind(pb, rb.items, false);   // key on: the mesh items read the engagement query's faces (their Refresh below does nothing; the inside check re-queries as the seat moves the ball)
   }
   // The ball starts the engagement touching or a little inside the rubber (the hit was found within a step). Go back to the instant it first touched: along its
   // own velocity, by depth / closing speed, at most one physics step; if it is not closing (pushed in by something else), seat it at touching distance along the normal.
   // Nothing is stored as compression.
   // Found inside a part with no earlier clean state (docs/unified-physics/archive/UnifiedContactDesign.txt 11, point 6): CheckOverlaps found no clean start of the sub-interval, or the engagement by proximity finds the ball already inside a
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

void ContactSystem::StepCordAlone(HitRubber *const band)
{
   ContactSystemImpl::CordAlone(band, g_subS);
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

void ContactSystemImpl::Regather(BallContact &e, const Vertex3Ds &pos, const float radiusU, const bool finder)
{
   const float reach = radiusU + 20.0f;   // as at engagement
   // A Verlet list: the surfaces within reach + 10 U are kept and collected again only after the ball has moved 5 U (or 10 ms have passed: a part can become collidable)
   constexpr float kMarginU = 10.0f, kMoveU = 5.0f;
   const uint32_t nowMs = g_pplayer->m_time_msec;
   const bool hole = PlayfieldHolesExist() && InPlayfieldHole(pos);
   const bool rescan = !(e.regathered && hole == e.inHole && nowMs - e.regatherMs < 10 && (pos - e.regatherPos).LengthSquared() < kMoveU * kMoveU);
   if (!rescan && !finder)
      return;
   std::vector<SurfaceContact> next;
   if (finder)
   {
      // shared_finder P2b (design.md section 3): the static mesh hosts every cycle from the ball's query (a face within R + reach), the others from the Verlet list below as today;
      // kept in the registry's order, as the registry loop gave them
      thread_local std::vector<ContactHost *> found;
      HitRubber::FinderHosts(e.ball, found);
      for (ContactHost *const o : found)
      {
         if (!o->m_set.contact)
            continue;
         if (hole && o->IsPlayfield())
            continue;   // no floor over a kicker's hole
         const auto old = std::find_if(e.items.begin(), e.items.end(), [o](const SurfaceContact &it) { return it.r == o; });
         if (old == e.items.end() && o->m_set.band)
            continue;   // (as below)
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
      if (!rescan)   // the other surfaces: the set as it is
         for (SurfaceContact &it : e.items)
            if (it.r != nullptr && it.r->FinderObjIndex() < 0)
               next.push_back(std::move(it));
   }
   if (rescan)
   {
   e.inHole = hole;
   e.regathered = true;
   e.regatherPos = pos;
   e.regatherMs = nowMs;
   for (ContactHost *const o : ContactHost::Registry())
   {
      if (finder && o->FinderObjIndex() >= 0)
         continue;   // a static mesh host: from the query, above
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
   }
   if (finder)
      std::stable_sort(next.begin(), next.end(), [](const SurfaceContact &a, const SurfaceContact &b) { return static_cast<const ContactHost *>(a.r)->m_regOrder < static_cast<const ContactHost *>(b.r)->m_regOrder; });
   if (!next.empty())   // never drop the whole set: the release test below decides when the ball is free
      e.items = std::move(next);
}

void ContactSystem::EngageNear(const double cycleS)
{
   ContactSystemImpl::EngageNear(cycleS);
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
      HitRubber::FinderBindOne(pb, it);   // shared_finder (key on): its faces from the ball's query of this cycle (the inside check re-queries if the ball is outside it)
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
   if (g_pplayer != nullptr && g_pplayer->m_physics != nullptr)
      g_pplayer->m_physics->m_noneZones.Placed(ball);   // the None triggers judge it where it is (placement rule)
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

// A step of dt seconds: the engine's whole sub-interval (Step), or a piece of it in an overlap replay (CheckOverlaps)
static void StepSeconds(const double dt)
{
   ContactSystemImpl::StepBalls(dt);
   if (dt > 0.0)
   {
      g_placed.clear();
      g_frameUsedS += dt;
   }
}

void ContactSystem::Step(const double hS)
{
   g_subS = hS;
   StepSeconds(hS);
}

double ContactSystem::SubIntervalS()
{
   return g_subS;
}

const ContactSystem::CycleStats &ContactSystem::LastCycle()
{
   return g_last;
}


// ---- The snapshot of the start of a sub-interval (docs/unified-physics/archive/UnifiedContactDesign.txt section 11, point 2) ---------------------------------------------------------------------------------
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
      int sleepNextId = 1;                                 // (the next pile's number: a pile formed in the step taken again gets the same one)
      std::vector<const HitBall *> placed;                 // g_placed: the balls a script or the driver put since the last cycle (the step reads it; Step clears it after)
      double sleepCycleDt = 0.0;
      double frameUsedS = 0.0;   // g_frameUsedS: how much of the physics step the sub-intervals before this one took (Body::Begin's share of gravity and nudge owed to the later ones reads it, and Step advances it)
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
   s.sleepNextId = ContactSystemImpl::g_sleepNextId;
   s.placed = g_placed;
   s.sleepCycleDt = ContactSystemImpl::g_sleepCycleDt;
   s.frameUsedS = g_frameUsedS;
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
   ContactSystemImpl::g_sleepNextId = s.sleepNextId;
   g_placed = s.placed;
   ContactSystemImpl::g_sleepCycleDt = s.sleepCycleDt;
   g_frameUsedS = s.frameUsedS;   // (the step being taken again starts where this sub-interval started in the physics step)
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
            f.contacts = FnvV(f.contacts, it.gs.ktb);
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

void ContactSystem::SnapshotTest(const double hS)
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
   Step(hS);
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


// ---- The overlap check and the replay of the sub-interval (docs/unified-physics/archive/UnifiedContactDesign.txt section 11, point 2) ------------------------------------------------------------------------
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

void ContactSystem::CheckOverlaps(const float hittime, const double hS)
{
   ContactSystemImpl::CheckOverlaps(hittime, hS);
}

// hittime: the sub-interval in the engine's time units, only for the free balls' straight flight (velocity x stock time); hS: the same sub-interval in seconds, for everything else
void ContactSystemImpl::CheckOverlaps(const float hittime, const double hS)
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
      RestoreSnap(FreeFlightOn() ? nullptr : &only);   // free_flight: every free ball is stepped again by StepSeconds below, so each starts from the snapshot (else one outside `only` would move twice)
      double t = 0.0;
      bool first = true;
      for (Rp &e : list)
      {
         if (e.t > tEnd)
            break;
         if (e.t > t)
         {
            s_append = !first;
            StepSeconds((e.t - t) * hS);
            s_append = false;
            first = false;
            t = e.t;
         }
         for (ContactHost *const o : Hosts())   // a part no engaged ball held has not moved yet: it is at the start of the sub-interval. It moves on to the engagement, as it would have alone (gapAt measures it so)
            if (ContactMovingPart *const mp = o->MovingPart())
            {
               const double idle = hS * e.t - mp->steppedS;
               if (idle > 1e-9)
               {
                  mp->AdvanceAlone(idle);
                  mp->steppedS += idle;   // (PartAdvance advances it for the rest of the sub-interval only)
               }
            }
         if (const BallSnap *const b0 = snapBall(e.pb); b0 != nullptr && !FreeFlightOn())   // (free_flight: the step above has already flown it to the engagement)
            e.pb->m_d.m_pos = b0->pos + b0->vel * (float)(e.t * (double)hittime);   // its free flight to the engagement
         EngageImpl(e.host, e.pb, "overlap");
      }
      if (tEnd > t)
      {
         s_append = !first;
         StepSeconds((tEnd - t) * hS);
         s_append = false;
      }
   };
   // The gap between a free ball (on its straight flight from the snapshot) and the nearest part, at share t, in the system simulate(t) left. A part that was stepped for less than t is advanced alone
   // for the rest on a scratch copy of its state.
   auto gapAt = [&](HitBall *pb, const double t) {
      const BallSnap *const b0 = snapBall(pb);
      BallSI bs;
      bs.R = (double)pb->m_d.m_radius * U;
      const Vertex3Ds p = FreeFlightOn() ? pb->m_d.m_pos : b0->pos + b0->vel * (float)(t * (double)hittime);   // free_flight: where simulate(t) stepped it
      const V3 c((double)p.x * U, (double)p.y * U, (double)p.z * U);
      double gmin = 1e9;
      for (ContactHost *const o : parts)
      {
         ContactMovingPart *const mp = o->MovingPart();
         const double idle = hS * t - mp->steppedS;
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
         o->PartAdvance();
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
