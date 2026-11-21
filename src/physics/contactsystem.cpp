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
#include "contactsurfacehost.h"
#include "kickerarm.h"
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

// VPX_RELOOK_TRACE=1 (diagnostics, TestingPlaybook 8; read once): each ball's energy-reach bound terms (ReachEnergy, every cycle), printed by the finder at an energy-reach re-look
// (ContactSurfaceHost::FinderInsideQ, declared there). Off: no string is built
bool RelookTraceOn()
{
   static const bool on = [] { const char *const e = std::getenv("VPX_RELOOK_TRACE"); return e != nullptr && e[0] == '1'; }();
   return on;
}
std::unordered_map<const HitBall *, std::string> &RelookTerms()
{
   static std::unordered_map<const HitBall *, std::string> m;
   return m;
}

// the part cull (runs/2026-11-20-part-cull/design.md): the table's moving parts and cord bands in the registry's order, listed at the first cycle of a table start
// (ContactSystemImpl::CullLists); each band's outline box (world, m) fixed when its cord is made; each body's travel bound this cycle (m, by its index in the cycle's bodies; set by
// CullPrepare for g_cullCy only); the misses told (check mode, once per part or band per table start); the work counters (diagnostics, logged every 1000 cycles)
namespace
{
   struct CullBand { ContactHost *host = nullptr; ContactSurfaceHost *band = nullptr; double lo[3] = { 0.0, 0.0, 0.0 }, hi[3] = { 0.0, 0.0, 0.0 }; bool box = false; };
   std::vector<ContactHost *> g_cullParts;
   std::vector<CullBand> g_cullBands;
   std::vector<char> g_cullPartTold, g_cullBandTold;
   bool g_cullListed = false, g_cullNamed = false;
   size_t g_cullRegSize = 0;
   std::vector<double> g_cullTravel;
   const void *g_cullCy = nullptr;
   double g_cullE = 0.0;
   std::vector<double> g_cullCap, g_cullA, g_cullAH;   // one-pass bound: each cull part's EnergyCap this cycle (J; -1 = off), its aLin (J^0.5) and aHertz (J^0.6); sized once per table start (CullPrepare)
   uint64_t g_cullLeftOut = 0, g_cullCapTests = 0;   // part caps left out of E (no body reaches the part at Efull) of those tested
   double g_cullBoundUs = 0.0;   // check mode: the bound's own time (us, summed)
   uint64_t g_hertzNoBound = 0;   // one-pass bound: closed forms with no bound (the Hertz tangent's slope 0.4 aHertz c^-0.6 >= 1) this table start; the notice is raised at the first (CullPrepare)
   bool g_hertzNoticed = false;
   // one-pass bound: E = E0 + aL sqrt(E) closed in one step (the larger root): sqrt(E) = (aL + sqrt(aL^2 + 4 E0)) / 2; aL = 0: E0
   inline double SqrtClosedE(const double E0, const double aL)
   {
      if (!(aL > 0.0))
         return E0;
      const double r = 0.5 * (aL + std::sqrt(aL * aL + 4.0 * std::max(E0, 0.0)));
      return r * r;
   }
   // one-pass bound (design.md AGREED items 5, 5b): E <= E0 + aL sqrt(E) + aH E^0.4 (aL: linear stops' squash, aH: Hertz stops'), one pass. E^0.4 is concave, so its tangent at c lies
   // above it: E^0.4 <= 0.6 c^0.4 + 0.4 c^-0.6 E; with s = 0.4 aH c^-0.6 < 1, (1 - s) E <= E0 + 0.6 aH c^0.4 + aL sqrt(E), closed as above. c = the sqrt form's E (aH = 0: E0, as
   // item 5b states; it is <= E, so a nearer tangent). s >= 1 (or c 0): no bound (INFINITY: the caller's no-bound path), counted in g_hertzNoBound
   inline double OnePassE(const double E0, const double aL, const double aH)
   {
      const double c = SqrtClosedE(E0, aL);
      if (!(aH > 0.0))
         return c;
      const double s = c > 0.0 ? 0.4 * aH * std::pow(c, -0.6) : INFINITY;
      if (!(s < 1.0))
      {
         g_hertzNoBound++;
         return INFINITY;
      }
      return SqrtClosedE((E0 + 0.6 * aH * std::pow(c, 0.4)) / (1.0 - s), aL / (1.0 - s));
   }
   // a ball's travel bound over T (the part cull's, CullPrepare, and the grouping's, ReachEnergy): EnergyReach(0, E, m, a, T) + hv T = sqrt(0*0 + 2E/m) T + 0.5 a T T + hv T, rounded up by
   // its own arithmetic's worst case. Every term >= 0 (no cancellation); 0*0, 2E, 0.5a and + 0 are exact; the longest chain of rounded operations is 5 (/m, sqrt, *T, +, + hv T; the a term
   // has 4: *T, *T, +, +; hv T has 2), so the computed value is at least the exact one times (1 - 5u) (u = eps/2, Higham's gamma_5); times (1 + 5 eps) = (1 + 10u), its own rounding
   // (1 - u) included, it is not below the exact value of this formula for its inputs E, a, hv as given
   inline double TravelBound(const double E, const double m, const double a, const double hv, const double T)
   {
      return (EnergyReach(0.0, E, m, a, T) + hv * T) * (1.0 + 5.0 * DBL_EPSILON);
   }
   uint64_t g_cullCycles = 0, g_cullTests = 0, g_cullSkipped = 0, g_cullBandTests = 0, g_cullBandSkipped = 0, g_cullMisses = 0, g_cullNoBound = 0, g_cullGatherDiff = 0;
}

namespace ContactSystem
{
   static bool s_unified = false, s_unifiedSet = false;
   static bool s_seatEditorBalls = false;   // [Features] seat_editor_balls of this table (physicsini.h), read likewise; off with the unified engine off
   static std::filesystem::path s_iniFile;   // the table's physics ini (PhysicsIni::FindUnified at SetTable; empty = none): every reader reads this one file
   static std::filesystem::file_time_type s_iniTime{};   // its write time when the readers last read it (SetTable, ReloadChangedIni)
   static bool s_aircoreInStep = false;     // [Features] aircore_in_step of this table (physicsini.h), read likewise: 1 = an air-core coil's pull is read at each step's midpoint (AircorePullForce, CordIsland); 0 = held for the cycle as before
   static std::set<std::string> s_retiredTold;   // the retired-setting notices told since this table started (NoticeRetired: one WARN and one notice each per table start)

   void SetTable(const std::filesystem::path &tableFile)
   {
      g_partEnergy.clear();
      g_partEnergyLines = 0;
      g_partEnergyTold = false;
      const std::filesystem::path file = PhysicsIni::FindUnified(tableFile);
      s_iniFile = file;
      s_iniTime = {};
      if (!file.empty())
      {
         std::error_code ec;
         s_iniTime = std::filesystem::last_write_time(file, ec);   // (before the readers read it: a save after this is seen at the next focus gain)
      }
      PhysicsIni::Features f;
      const bool found = !file.empty() && PhysicsIni::ReadFeatures(file, f);
      s_unified = found && f.unified;   // no file, no [Features] or no  unified = 1  = stock VPX physics
      s_unifiedSet = true;
      s_seatEditorBalls = s_unified && f.seatEditorBalls;
      s_aircoreInStep = s_unified && f.aircoreInStep;
      s_retiredTold.clear();
      g_cullListed = g_cullNamed = false;   // (the part cull: its lists, told misses and counters per table start)
      g_cullParts.clear();
      g_cullBands.clear();
      g_cullTravel.clear();
      g_cullCy = nullptr;
      g_cullCycles = g_cullTests = g_cullSkipped = g_cullBandTests = g_cullBandSkipped = g_cullMisses = g_cullNoBound = g_cullGatherDiff = 0;
      g_cullCap.clear();
      g_cullA.clear();
      g_cullAH.clear();
      g_hertzNoBound = 0;
      g_hertzNoticed = false;
      g_cullLeftOut = g_cullCapTests = 0;
      g_cullBoundUs = 0.0;
      GroupStep::partFinderProf = ContactProf::CheckOn() || GroupStep::FinderCheckEnv();   // (its counters: diagnostics only)
      GroupStep::partFinderCount.clear();
      GroupStep::skipOwnTested = GroupStep::skipOwnSkipped = GroupStep::skipEsComputed = GroupStep::skipTripItems = GroupStep::skipDeferred = GroupStep::skipTripAttempts = 0;
      GroupStep::partFinderFirstMiss.clear();
      ContactSettings::StartPlay(tableFile);   // (after the switch above: the settings ask it)
      if (s_unified)
      {
         PLOGI << "Contact system: unified engine on (" << file.string() << "; the unified cycle: no stock search, Collide or Contact)";
         for (const std::string &key : f.retiredOff)   // a setting whose 0 once sent an object or a law to stock (PhysicsIni::RETIRED)
            NoticeRetired(key + " = 0 is ignored: no stock path on the unified engine");
         PLOGI << "Contact system: [Features] aircore_in_step = " << (s_aircoreInStep ? 1 : 0) << (s_aircoreInStep ? " (an air-core coil's pull read at each step's midpoint)" : " (an air-core coil's pull held for the cycle)");
      }
      else
         PLOGI << "Contact system: stock VPX physics, every new model off (" << (file.empty() ? std::string("no physics ini") : file.string() + " has no [Features] unified = 1") << ")";
   }

   void ReloadChangedIni()
   {
      if (!s_unifiedSet || s_iniFile.empty())
         return;   // (no table read yet, or no physics ini at table load: nothing to look at)
      std::error_code ec;
      const auto t = std::filesystem::last_write_time(s_iniFile, ec);
      if (ec || t == s_iniTime)
         return;   // (unchanged; or gone or locked: everything keeps its settings, the next focus gain looks again)
      s_iniTime = t;
      PLOGI << "Physics ini " << s_iniFile.string() << " changed: reloading [Contact], the magnets and every contact host (the bats included); [Features] take the next play";
      ContactSettings::Reload();
      Magnets::Reload();
      ContactHost::ReloadAll();
   }

   bool Unified()
   {
      if (!s_unifiedSet)
         PLOGE << "Contact system: Unified() asked before the table's physics ini was read (SetTable); answering stock";
      return s_unified;
   }

   void NoticeRetired(const std::string &msg)
   {
      if (!s_retiredTold.insert(msg).second)
         return;
      PLOGW << msg;
      if (g_pplayer != nullptr && g_pplayer->m_liveUI != nullptr)
         g_pplayer->m_liveUI->PushNotification(msg, 8000);
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
// one entry per ball counters (diagnostics only, deterministic; running totals logged every 10 s of table time and at the table's end, ReleaseHost): the nets that should never fire (stock: a stock collider found a ball
// on a converted surface; noPast: a ball seated by minimum translation; overlap: the overlap replay's scan found a ball; snap: a snapshot restore with a changed entry list),
// items joined / left between cycles, and the tightness: per kind, the ball-surface pairs kept in a ball's solve at a cycle's end and those that touched in that cycle
struct QueryEngageCounts
{
   uint64_t made = 0, cleared = 0, noPast = 0, overlap = 0, snap = 0, joined = 0, left = 0, cycles = 0, firstSightZero = 0;   // firstSightZero: items new in the query found inside (zero = start depth)
   std::map<std::string, std::pair<uint64_t, uint64_t>> kept;   // kind -> (kept, touched)
};
static QueryEngageCounts g_qe;
static std::map<std::string, std::pair<double, int>> g_partWorkMax;   // the energy reach: per part type, the largest work of one part in one cycle, J, and the table time (ms) of it (diagnostics only)
static void QueryEngageLog(const char *const when)   // the one entry per ball counters' line (running totals; diagnostics only)
{
   std::ostringstream o;
   o << "Query engage (" << when << "): cycles " << g_qe.cycles << ", entries made " << g_qe.made << ", items cleared " << g_qe.cleared << "; nets: no past " << g_qe.noPast << ", overlap " << g_qe.overlap
     << ", snapshot " << g_qe.snap << "; items joined " << g_qe.joined << ", left " << g_qe.left << ", found inside at first sight (zero = depth) " << g_qe.firstSightZero << "; tightness (kept in solve at cycle end / touched in the cycle):";
   for (const auto &k : g_qe.kept)
      o << " [" << k.first << " " << k.second.first << " / " << k.second.second << "]";
   PLOGI << o.str();
   std::ostringstream w;
   w << "Energy reach part work (" << when << "): largest per cycle per part type, J (table ms):";
   for (const auto &k : g_partWorkMax)
      w << " [" << k.first << " " << k.second.first << " (" << k.second.second << ")]";
   PLOGI << w.str();
}
static double g_cycleS = 0.0;       // the cycle in seconds, as ContactSystem::Step was given it (the parts and cords advanced alone after it read it)
static std::vector<const HitBall *> g_placed;   // the balls placed since the last cycle: the load saw their position changed outside the engine, or loaded them first (NotePlaced)
// a placement (GroupStepDesign 17.4): a position move no step did, and no fault; the None triggers judge the ball where it is. Called by the load only (OwnLoad).
static void NotePlaced(const HitBall *const ball)
{
   g_placed.push_back(ball);
   if (g_pplayer != nullptr && g_pplayer->m_physics != nullptr)
      g_pplayer->m_physics->m_noneZones.Placed(ball);   // the None triggers judge it where it is (placement rule)
}
static ContactSystem::CycleStats g_last;
static std::vector<std::function<void()>> g_deferred;      // events raised inside the step, fired after it (Finish)
// The test driver's cost split of this cycle (docs/unified-physics/TestMatrixDesign.txt section 3 point 10 b-d): null unless the driver times the cycle in detail. Per island: what set each
// substep (TestDriver::ContactCost::Why) and the smallest substep.
void ContactSystem::DeferEvent(std::function<void()> f) { g_deferred.push_back(std::move(f)); }
static constexpr float kFlipperHorizonVpt = 1.0f;   // trace only: VPX time units (1 VPT = 10 ms) a bat is looked ahead for in the flip_clear numbers

struct ContactSystemImpl
{
   static void StepBalls(double dt);   // dt: seconds
   static void Release(HitBall *ball, const char *why, bool remove = false);   // remove: the ball leaves play (one entry per ball: only then is its entry erased; otherwise its items are cleared)
   static void Regather(BallContact &e, const Vertex3Ds &pos, float radiusU, int part = 0, double cullTravel = -1.0);   // the static mesh hosts come from the ball's query (shared finder), the others from a Verlet list
   // (the energy reach, part: 1 = the others only (before the grouping), 2 = the static mesh hosts only (from the query at the ball's energy travel, after it); 0 = both)
   // (cullTravel, part 1: the ball's travel bound this cycle, m (the part cull, CullPrepare): only the moving parts whose region it reaches are candidates; < 0 = every one)
   static bool EngageImpl(ContactHost *host, HitBall *pb, const char *how);
   static void EngageNear(double dt);   // dt: the cycle in seconds
   // one entry per ball: the ball's one entry, made when it has none (its first cycle, or the first collider to find it before that); no seat, no items: its query fills them
   static BallContact &EnsureOne(HitBall *const pb)
   {
      if (pb->m_contact != nullptr)
         return *pb->m_contact;
      BallContact rb;
      rb.ball = pb;
      rb.posStart = pb->m_d.m_pos;
      g_engaged.push_back(std::move(rb));
      pb->m_contact = &g_engaged.back();
      g_qe.made++;
      return g_engaged.back();
   }
   // every ball in play (not locked in a kicker) has its entry from its first cycle
   static void EnsureEntries()
   {
      for (auto *const ball : g_pplayer->m_vball)
      {
         HitBall *const pb = &ball->m_hitBall;
         if (pb->m_contact == nullptr && PairEligible(pb))
            EnsureOne(pb);
      }
   }

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
      bool dragged = false;     // the debug mouse drags the ball: held against gravity and the nudge (Ball::fExt), no magnets
      bool extTaken = false;    // Sources::Begin took this cycle's gravity and nudge out of the velocity (they act in the substeps)
      double endPair = 0.0;     // the energy check: its half of what its ball pairs hold at the island's end
      double magAcc = 0.0;      // m/s2: the magnets' acceleration of this ball at the start of the step (the substep bound)
      double eddyRate = -1.0;   // eddy-current braking of the spin across the coil field, 1/s: evaluated once per cycle (negative = not yet); eddyDir = the field direction
      V3 eddyDir;
      double eStart = 0.0, eScale = 0.0;   // the energy check: the ball's energy at the cycle's start, and what the allowance is a share of
      ContactSurfaceHost *band = nullptr;   // a band that holds the ball (the group step takes the whole island)
      bool held = false;           // a ball a kicker holds (parked) or carries (KickerDesign 1h): infinite mass in the pair law, its velocity given (v), no gravity, spin 0
      V3 heldA;                    // its given acceleration, m/s2 (the arm's path's centripetal pull; 0 parked)
      double feedE = -1.0, reachS = -1.0;   // the energy reach (ReachEnergy): its group's feed energy (J) and its travel over the cycle (m); -1 = not set
      double vTop = -1.0, eTop = -1.0;      // energy skip test (ReachEnergy): its top speed over the cycle (m/s, EnergyReach's v_max plus a held ball's given speed) and its group's energy at that speed, 1/2 m v_max^2 (J); -1 = not set
      int eGroup = -1;                      // energy skip test: its energy group (ReachEnergy's root)
      double magReach = 0.0;  // the energy reach: |a magnet's pull| on it now, m/s2 (Sources::Begin; on or off the playfield model)
      double gravAcc = 0.0;   // the energy reach: |gravity + nudge| on it, m/s2 (Sources::Begin; 0 for a held ball, whose motion is given), whether or not the substeps take it
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
      std::vector<int> ballOrder;        // the bodies sorted along the table's x axis by their centre at the cycle start (SortBallsX; ball pair candidates)
      std::vector<std::pair<int, int>> near0;   // the ball pairs whose x intervals [x - R, x + R] overlap (i < j, in i, j order): every pair that may touch now (FindIslands)
   };

   // Sources: anything that acts on one body with no partner. Gravity and the nudge act in the substeps, not as one kick per physics step (a contact stepped without gravity
   // and kicked once per step rattles at the step rate); a magnet's pull follows the ball, so it is evaluated at every substep. Both are velocity increments.
   struct Sources
   {
      // Gravity and nudge: the engine has added this step's increment to the ball already (the stock contacts of everything else rely on it): this cycle's share of it is taken
      // back here and given back in small pieces, so the step ends with the same total.
      // Only while the playfield is one of the surfaces of the contact, or the ball has no contact at all: with the stock plane contact under it the ball needs the engine's kick to stay on the plane.
      // the energy reach: gravity and the nudge act on every ball every cycle, whichever model gives them (the stock engine's kick or the substeps): the drive's acceleration
      // (ReachEnergy's forces term; also a new ball's early query, EngageNear). Reads only: the ball's velocity is not changed here
      static void ReachAcc(Body &w)
      {
         HitBall *const pb = w.pb;
         const Vertex3Ds g = pb->ExternalVelocityIncrement(false);
         w.gravAcc = V3(g.x * VS, g.y * VS, g.z * VS).Length() / PHYSICS_STEPTIME_S;
         // a magnet's pull on a magnetic ball while a coil carries current, whether or not the substeps give it (the stock increment gives it to a ball off the playfield model);
         // at the ball's place now (the field's change over the cycle's travel is not bounded here)
         w.magReach = 0.0;
         if (pb->m_magnetic && Magnets::Active())
         {
            const Vertex3Ds a = pb->ExternalVelocityIncrement(true);
            w.magReach = V3((a.x - g.x) * VS, (a.y - g.y) * VS, (a.z - g.z) * VS).Length() / PHYSICS_STEPTIME_S;
         }
      }

      static void Begin(Body &w, const double dt, const bool always = false)
      {
         HitBall *const pb = w.pb;
         if (w.held)   // a ball a kicker holds is not stepped (KickerDesign 1h): no gravity, nudge or magnets were given, none are taken or owed
         {
            w.extTaken = true;
            w.magnetsOn = false;
            w.magAcc = 0.0;
            w.extBaseVS = V3();
            return;
         }
         ReachAcc(w);
         // the engine added nothing (HitBall::UpdateVelocities): every ball's gravity, nudge and magnets act in the substeps; the velocity is the true one, left as it is
         BallControl *const ballControl = g_pplayer->m_physics->GetBallControl();
         w.dragged = ballControl != nullptr && pb->m_pBall == ballControl->GetDraggedBall();
         const Vertex3Ds extInc = pb->ExternalVelocityIncrement();   // gravity, nudge, magnets at the start position
         const Vertex3Ds extBase = pb->ExternalVelocityIncrement(false);   // the part that does not depend on the position
         w.extTaken = true;
         w.extBaseVS = V3(extBase.x * VS, extBase.y * VS, extBase.z * VS);
         w.magnetsOn = !w.dragged && pb->m_magnetic && Magnets::Active();
         w.magAcc = w.magnetsOn ? V3((extInc.x - extBase.x) * VS, (extInc.y - extBase.y) * VS, (extInc.z - extBase.z) * VS).Length() / PHYSICS_STEPTIME_S : 0.0;
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
   // Always on (DevelopmentRules 8): a ball whose centre, velocity or spin is not a finite number reached the step: a fault upstream (a freed ball, a script, a solver). Reported only,
   // never repaired: one WARN and one on-screen notice per table start (EndTable resets it), naming the ball, its place in the player's ball list (-1 = not a ball in play) and its state.
   inline static bool s_nonFiniteTold = false;
   static void NonFiniteBall(const HitBall *const pb)
   {
      s_nonFiniteTold = true;
      const char *const name = TestDriver::BallName(pb);
      int index = -1;
      for (size_t i = 0; i < g_pplayer->m_vball.size(); i++)
         if (&g_pplayer->m_vball[i]->m_hitBall == pb)
            index = (int)i;
      const Vertex3Ds &p = pb->m_d.m_pos, &v = pb->m_d.m_vel, &am = pb->m_angularmomentum;
      std::ostringstream o;
      o << "Contact: ball " << (name != nullptr ? name : "") << " (index " << index << " of " << g_pplayer->m_vball.size() << " balls in play) has a non-finite state at t=" << g_pplayer->m_time_msec
        << " ms: centre (" << p.x << ", " << p.y << ", " << p.z << ") velocity (" << v.x << ", " << v.y << ", " << v.z << ") spin (" << am.x << ", " << am.y << ", " << am.z
        << "): a fault upstream, reported only (once per table start)";
      PLOGW << o.str();
      if (g_pplayer->m_liveUI != nullptr)
         g_pplayer->m_liveUI->PushNotification(o.str(), 8000);
   }
   static void LoadBody(Body &w, HitBall *const pb)
   {
      if (!s_nonFiniteTold)
      {
         const Vertex3Ds &p = pb->m_d.m_pos, &v = pb->m_d.m_vel, &am = pb->m_angularmomentum;
         if (!(std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z) && std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z) && std::isfinite(am.x) && std::isfinite(am.y) && std::isfinite(am.z)))
            NonFiniteBall(pb);
      }
      w.pb = pb;
      w.bs.R = (double)pb->m_d.m_radius * U;
      w.bs.m = kMassUnitKg * (double)pb->m_d.m_mass;
      w.lScale = (double)pb->Inertia() * DEFAULT_STEPTIME_S;
      w.posStart = pb->m_d.m_pos;
      Vertex3Ds aU;
      w.held = KickerHold::Held(pb, &aU);
      w.heldA = V3((double)aU.x, (double)aU.y, (double)aU.z) * (U / (DEFAULT_STEPTIME_S * DEFAULT_STEPTIME_S));   // U/VPT^2 -> m/s2
   }

   // ---- the ball's owned state (to-do 1.51, runs/2026-11-21-ball-state/design.md; HitBall::m_own) ----
   // The load's change test, part by part (design 3.1-3.4): position, velocity and spin each keep the exact double while their floats are what the engine last wrote for them
   // (a script writing back the very float it read changes nothing); a part whose floats differ was changed outside the engine and takes the float, which replaces the owned
   // value (a script that sets only the velocity leaves the exact position). A changed position, or a ball's first load (no state yet), is a placement (NotePlaced).
   // Idempotent: the floats it takes become the recorded ones, so every later load in the cycle sees no change. A change smaller than the float's step is not seen (design 5.2: no work).
   static void OwnLoad(HitBall *const pb)
   {
      HitBall::OwnState &s = pb->m_own;
      const bool first = !s.set;
      if (first || !(pb->m_d.m_pos == s.pos))
      {
         s.pos = pb->m_d.m_pos;
         s.c[0] = (double)s.pos.x * U; s.c[1] = (double)s.pos.y * U; s.c[2] = (double)s.pos.z * U;
         NotePlaced(pb);
      }
      if (first || !(pb->m_d.m_vel == s.vel))
      {
         s.vel = pb->m_d.m_vel;
         s.v[0] = (double)s.vel.x * VS; s.v[1] = (double)s.vel.y * VS; s.v[2] = (double)s.vel.z * VS;
      }
      if (first || !(pb->m_angularmomentum == s.am))
      {
         s.am = pb->m_angularmomentum;
         const double wScale = 1.0 / ((double)pb->Inertia() * DEFAULT_STEPTIME_S);
         s.w[0] = (double)s.am.x * wScale; s.w[1] = (double)s.am.y * wScale; s.w[2] = (double)s.am.z * wScale;
      }
      s.set = true;
   }
   static V3 OwnC(const HitBall *const pb) { return V3(pb->m_own.c[0], pb->m_own.c[1], pb->m_own.c[2]); }   // (after the load: every ball has a state)
   static V3 OwnV(const HitBall *const pb) { return V3(pb->m_own.v[0], pb->m_own.v[1], pb->m_own.v[2]); }
   static V3 OwnW(const HitBall *const pb) { return V3(pb->m_own.w[0], pb->m_own.w[1], pb->m_own.w[2]); }
   // The engine's one write of a ball's state (design 4.1: the step's write-back, the sleep's zero, an engagement's seat): the owned doubles and the floats together, the floats
   // recorded as written, so the next load sees no outside change. Writing the floats does no work (17.5): the ball goes on from the doubles.
   static void WriteBall(HitBall *const pb, const V3 &c, const V3 &v, const V3 &w)
   {
      HitBall::OwnState &s = pb->m_own;
      const double lScale = (double)pb->Inertia() * DEFAULT_STEPTIME_S;
      s.c[0] = c.x; s.c[1] = c.y; s.c[2] = c.z;
      s.v[0] = v.x; s.v[1] = v.y; s.v[2] = v.z;
      s.w[0] = w.x; s.w[1] = w.y; s.w[2] = w.z;
      pb->m_d.m_pos = Vertex3Ds((float)(c.x / U), (float)(c.y / U), (float)(c.z / U));
      pb->m_d.m_vel = Vertex3Ds((float)(v.x / VS), (float)(v.y / VS), (float)(v.z / VS));
      pb->m_angularmomentum = Vertex3Ds((float)(w.x * lScale), (float)(w.y * lScale), (float)(w.z * lScale));
      s.pos = pb->m_d.m_pos;
      s.vel = pb->m_d.m_vel;
      s.am = pb->m_angularmomentum;
      s.set = true;
   }

   // SI state of a ball: its owned state, after the change test (the floats only where something outside the engine changed them)
   static void LoadState(Body &w)
   {
      HitBall *const pb = w.pb;
      OwnLoad(pb);
      w.c = OwnC(pb);
      w.v = OwnV(pb);
      w.w = OwnW(pb);
      if (w.held)
         w.w = V3();   // (spin zero: KickerDesign 1h)
   }

   // ---- the part cull (runs/2026-11-20-part-cull/design.md) ----
   // the table's moving parts and cord bands, in the registry's order, once per table start (and again if the registry grows)
   static void CullLists()
   {
      const std::vector<ContactHost *> &reg = ContactHost::Registry();
      if (g_cullListed && g_cullRegSize == reg.size())
         return;
      g_cullListed = true;
      g_cullNamed = false;
      g_cullRegSize = reg.size();
      g_cullParts.clear();
      g_cullBands.clear();
      for (ContactHost *const o : reg)
      {
         if (o->FinderObjIndex() < 0 && o->MovingPart() != nullptr)
            g_cullParts.push_back(o);
         if (ContactSurfaceHost *const b = o->BandHost(); b != nullptr && b->CordModel())
         {
            CullBand cb;
            cb.host = o;
            cb.band = b;
            g_cullBands.push_back(cb);
         }
      }
      g_cullPartTold.assign(g_cullParts.size(), 0);
      g_cullBandTold.assign(g_cullBands.size(), 0);
   }
   // a band's region, fixed (design 2R, amended: a ball only pushes a band inward, so its cord stays inside its outline): a box round its outline with the sling plate over its whole
   // stroke, the cord's thickness included. The full test with no ball held (BandQueryGap: each span's rest line and each post's cord circle less rc and the span's bend, and on the
   // plate the face now less rc and its move this cycle) lies inside it: each span's ends and each post's circle grown by rc + |bendMost| + the plate's whole stroke move, the
   // plate's face edges at the two ends of its stroke grown by rc + that move (a box round the outline, not the outline itself). The plate's stroke: the face's offset y(phi) for phi
   // in [0, PhiStop], at both ends and where dy = 0 ((a + rc) sin t = -h, t = phi + rest; CordFixed::Part::Limit)
   static void CullBandBox(CullBand &cb)
   {
      const CordFixed::Part &p = *cb.band->Cord()->fixed;
      double yMin = 0.0, yMax = 0.0, plate = 0.0;
      if (p.geo.hasPlate)
      {
         const double stop = p.PhiStop();
         const auto y = [&p](const double ph) { double yy = 0.0, dy = 0.0; p.Limit(ph, yy, dy); return yy; };
         yMin = std::min(y(0.0), y(stop));
         yMax = std::max(y(0.0), y(stop));
         const double ar = p.pl.a + p.rc;
         if (std::fabs(ar) > std::fabs(p.pl.h))
            for (const double t : { std::asin(-p.pl.h / ar), GroupStep::kPi - std::asin(-p.pl.h / ar) })
               if (t - p.pl.rest > 0.0 && t - p.pl.rest < stop)
               {
                  yMin = std::min(yMin, y(t - p.pl.rest));
                  yMax = std::max(yMax, y(t - p.pl.rest));
               }
         plate = yMax - yMin;
      }
      for (int k = 0; k < 3; k++)
      {
         cb.lo[k] = INFINITY;
         cb.hi[k] = -INFINITY;
      }
      const auto add = [&cb](const RubberContact::V3 &q, const double g) {
         const double v[3] = { q.x, q.y, q.z };
         for (int k = 0; k < 3; k++)
         {
            cb.lo[k] = std::min(cb.lo[k], v[k] - g);
            cb.hi[k] = std::max(cb.hi[k], v[k] + g);
         }
      };
      for (const CordFixed::Span &sp : p.span)
      {
         add(sp.A, p.rc + std::fabs(sp.bendMost) + plate);
         add(sp.B, p.rc + std::fabs(sp.bendMost) + plate);
      }
      for (const CordFixed::Part::Circle &q : p.circles)
         if (!q.twin)
            add(RubberContact::V3(q.x, q.y, p.zc), q.r + p.rc + plate);
      if (p.geo.hasPlate)
      {
         const CordFixed::Span &sp = p.span[(size_t)p.pl.span];
         for (const double yy : { yMin, yMax })
            for (const double x : { p.pl.x0, p.pl.x1 })
               add(sp.A + sp.t * x + sp.yh * (p.pl.sgn * yy), p.rc + plate);
      }
      cb.box = true;
      PLOGI << "Part cull: band " << cb.host->Name() << " outline box x " << cb.lo[0] << " .. " << cb.hi[0] << ", y " << cb.lo[1] << " .. " << cb.hi[1] << ", z " << cb.lo[2] << " .. " << cb.hi[2]
            << " m (plate stroke move " << plate * 1e3 << " mm)";
   }
   static double CullBandDist(const CullBand &cb, const RubberContact::V3 &c)
   {
      const double dx = std::max({ cb.lo[0] - c.x, 0.0, c.x - cb.hi[0] }), dy = std::max({ cb.lo[1] - c.y, 0.0, c.y - cb.hi[1] }), dz = std::max({ cb.lo[2] - c.z, 0.0, c.z - cb.hi[2] });
      return std::sqrt(dx * dx + dy * dy + dz * dz);
   }
   // check mode (VPX_FINDER_CHECK=1): a part or band the cull left out that the full test takes: a missed candidate (DevelopmentRules 8: WARN and on screen, once per part or band per table start)
   static void CullMiss(const char *const kind, const std::string &name, char &told, const Body &w, const double full, const double bound, const double travel)
   {
      g_cullMisses++;
      if (told)
         return;
      told = 1;
      const char *const ball = TestDriver::BallName(w.pb);
      std::ostringstream o;
      o << "Part cull: ball " << (ball != nullptr ? ball : "") << " can reach " << kind << " " << name << " this cycle by the full test (" << full * 1e3 << " mm, within its travel "
        << w.reachS * 1e3 << " mm) but the cull left it out (its region less R " << bound * 1e3 << " mm, the ball's travel bound " << travel * 1e3 << " mm) at t=" << g_pplayer->m_time_msec
        << " ms: a missed candidate (check mode; first for this " << kind << ")";
      PLOGW << o.str();
      if (g_pplayer->m_liveUI != nullptr)
         g_pplayer->m_liveUI->PushNotification(o.str(), 6000);
   }
   // once per cycle before the grouping (design 1b): an upper bound on any group's feed in ReachEnergy, E = every body's kinetic energy and spin and what its contacts held (endStored),
   // every live ball pair's stored energy (PairHolds), every cord's (heldJ) and its sling coil's work bound (CoilWorkBound), every part's squash (squashJ) and the caps of the parts some
   // ball can reach (EnergyCap) with their stops' squash closed in one pass (OnePassE; runs/2026-11-20-energy-bound/design.md, AGREED 2026-11-20). No group can be fed more (energy
   // conservation). Each body's travel bound: EnergyReach(0, E, its mass, every body's forces (gravity, nudge, magnets) summed, T) plus every held ball's given speed x T. No bound (a
   // closed form with none): no cull this cycle (logged)
   static void CullPrepare(Cycle &cy, const size_t first)
   {
      CPROF(CullPrepare);
      CullLists();
      g_cullCycles++;
      g_cullCy = &cy;
      const double T = cy.dt;
      const size_t nw = cy.bodies.size();
      double Eb = 0.0, af = 0.0, hv = 0.0;
      for (size_t i = 0; i < nw; i++)
      {
         const Body &w = cy.bodies[i];
         if (w.held)
            hv += 2.0 * (w.v.Length() + w.heldA.Length() * T);
         else
            Eb += 0.5 * w.bs.m * w.v.Dot(w.v) + 0.5 * w.bs.Inertia() * w.w.Dot(w.w) + (w.e != nullptr ? std::max(w.e->endStored, 0.0) : 0.0);
         af += w.gravAcc + w.magReach;
      }
      for (const BallPairState &ps : g_ballPairs)
      {
         const Body *a = nullptr, *b = nullptr;
         for (size_t i = 0; i < nw; i++)
         {
            a = cy.bodies[i].pb == ps.a ? &cy.bodies[i] : a;
            b = cy.bodies[i].pb == ps.b ? &cy.bodies[i] : b;
         }
         if (a != nullptr && b != nullptr)
            Eb += std::max(PairHolds(cy, *a, *b, ps), 0.0);
      }
      for (CullBand &cb : g_cullBands)
      {
         if (cb.host->m_set.band && cb.host->m_enabled)
            cb.band->CordEnsure();   // (as BandQuery does this cycle: the cord is made once)
         if (!cb.band->CordOn())
            continue;
         Eb += cb.band->Cord()->heldJ + cb.band->Cord()->fixed->CoilWorkBound(T);
         if (!cb.box)
            CullBandBox(cb);
      }
      for (ContactHost *const o : g_cullParts)
         if (o->m_enabled && o->m_set.contact)
            Eb += o->MovingPart()->squashJ;
      const bool check = GroupStep::FinderCheckEnv();
      const auto t0 = check ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point();
      double E = -1.0;
      {
         // E = Eb + the caps of the parts some ball can reach (EnergyCap: each from its state now, no move)
         // + the squash of stops past their touch, aLin sqrt(E) (linear) + aHertz E^0.4 (Hertz), closed in one pass (OnePassE). Locality (item 2) with no circle:
         // first Efull with EVERY part's cap (an upper bound on any E below), each body's travel and each part's region at Efull (the cull's own test, Regather part 1); a part no body
         // reaches at Efull is reached by none at any smaller E (travel and region only grow with E), so its cap is left out of E. E <= Efull, so the cull at E is inside the test at Efull
         const size_t np = g_cullParts.size();
         g_cullCap.resize(np);
         g_cullA.resize(np);
         g_cullAH.resize(np);
         double capAll = 0.0, aAll = 0.0, hAll = 0.0;
         {
            CPROF(BoundCaps);
            for (size_t k = 0; k < np; k++)
            {
               ContactHost *const o = g_cullParts[k];
               g_cullCap[k] = -1.0;
               g_cullA[k] = g_cullAH[k] = 0.0;
               if (o->m_enabled && o->m_set.contact)
               {
                  g_cullCap[k] = o->MovingPart()->EnergyCap(T, g_cullA[k], g_cullAH[k]);
                  capAll += g_cullCap[k];
                  aAll += g_cullA[k];
                  hAll += g_cullAH[k];
               }
            }
         }
         double Efull;
         {
            CPROF(BoundClosed);
            Efull = OnePassE(Eb + capAll, aAll, hAll);
         }
         E = Efull;
         if (std::isfinite(Efull))
         {
            CPROF(BoundLocality);
            for (ContactHost *const o : g_cullParts)
               o->MovingPart()->CullCycle(Efull, T);
            g_cullTravel.assign(nw, -1.0);
            for (size_t i = 0; i < nw; i++)
               g_cullTravel[i] = TravelBound(Efull, cy.bodies[i].bs.m, af, hv, T);
            double E0 = Eb, a = 0.0, h = 0.0;
            for (size_t k = 0; k < np; k++)
            {
               if (g_cullCap[k] < 0.0)
                  continue;
               ContactHost *const o = g_cullParts[k];
               bool reached = o->m_set.band;   // (a band's host: the band's own rule, BandQuery; always counted)
               for (size_t i = 0; i < nw && !reached; i++)
               {
                  const Body &w = cy.bodies[i];
                  const RubberContact::V3 c((double)w.posStart.x * U, (double)w.posStart.y * U, (double)w.posStart.z * U);   // (as Regather part 1)
                  reached = o->MovingPart()->CullDist(c) - (double)w.pb->m_d.m_radius * U <= g_cullTravel[i];
               }
               if (reached)
               {
                  E0 += g_cullCap[k];
                  a += g_cullA[k];
                  h += g_cullAH[k];
               }
               else
                  g_cullLeftOut++;
               g_cullCapTests++;
            }
            CPROF(BoundClosed);
            E = OnePassE(E0, a, h);
         }
         if (check)
            g_cullBoundUs += std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count();
      }
      if (!g_cullNamed && E >= 0.0)   // (once per table start: the bound's terms)
      {
         std::ostringstream o;
         o << "Part cull: the table's energy bound at t=" << g_pplayer->m_time_msec << " ms: " << E * 1e3 << " mJ; balls, stored, cords and squash " << Eb * 1e3
           << " mJ; each part's cap (mJ), aLin and aHertz (EnergyCap):";
         for (ContactHost *const o2 : g_cullParts)
            if (o2->m_enabled && o2->m_set.contact)
            {
               double a = 0.0, ah = 0.0;
               const double cap = o2->MovingPart()->EnergyCap(T, a, ah);
               o << " " << o2->Name() << " " << cap * 1e3 << " (aLin " << a << ", aHertz " << ah << ")";
            }
         PLOGI << o.str();
      }
      g_cullTravel.clear();
      if (g_hertzNoBound > 0 && !g_hertzNoticed)   // (one-pass bound: the first closed form with no bound this table start, here or in last cycle's grouping; rule 8, once)
      {
         g_hertzNoticed = true;
         std::ostringstream o;
         o << "Energy bound: a Hertz stop's squash term (plunger, auto plunger or button cap pressed on its stop) has no one-pass bound at t=" << g_pplayer->m_time_msec
           << " ms (0.4 aHertz E0^-0.6 >= 1): every moving part is a candidate that cycle. Report the table and vpinball.log";
         PLOGW << o.str();
         g_pplayer->m_liveUI->PushNotification(o.str(), 6000);
      }
      if (!(E >= 0.0) || !std::isfinite(E))
      {
         if (g_cullNoBound++ == 0)
            PLOGW << "Part cull: no bound on the table's energy (a closed form with none) at t=" << g_pplayer->m_time_msec << " ms: every moving part and band is a candidate this cycle (first of this table start)";
         return;
      }
      g_cullE = E;
      {
         CPROF(CullRegions);
         g_cullTravel.assign(nw, -1.0);
         for (size_t i = first; i < nw; i++)
            g_cullTravel[i] = TravelBound(E, cy.bodies[i].bs.m, af, hv, T);
         for (ContactHost *const o : g_cullParts)
            o->MovingPart()->CullCycle(E, T);
      }
      if (!g_cullNamed)   // (once per table start: what the cull covers; a part with no region is a candidate for every ball)
      {
         g_cullNamed = true;
         std::string none;
         for (ContactHost *const o : g_cullParts)
            if (o->MovingPart()->CullDist(o->MovingPart()->Origin()) == -INFINITY)
               none += " " + o->Name();
         size_t boxes = 0;
         for (const CullBand &cb : g_cullBands)
            boxes += cb.box ? 1 : 0;
         for (ContactHost *const o : g_cullParts)
            if (const std::string note = o->MovingPart()->RegionNote(E, T); !note.empty())
               PLOGI << "Part cull: " << o->Name() << " region (mm): " << note;
         PLOGI << "Part cull: " << g_cullParts.size() << " moving parts, " << g_cullBands.size() << " cord bands (" << boxes << " with an outline box); parts with no region (a candidate for every ball):"
               << (none.empty() ? std::string(" none") : none);
      }
      if (g_cullCycles % 1000 == 0)   // (diagnostics: the running totals)
         PLOGI << "Part cull totals: " << g_cullCycles << " cycles; ball-part tests " << g_cullTests << ", culled " << g_cullSkipped << "; ball-band tests " << g_cullBandTests << ", bands skipped "
               << g_cullBandSkipped << "; misses (check mode) " << g_cullMisses << ", gather test differences (check mode) " << g_cullGatherDiff << ", cycles with no bound " << g_cullNoBound
               << "; E now " << g_cullE * 1e3 << " mJ; part caps left out (no ball reaches the part) " << g_cullLeftOut << " of " << g_cullCapTests << "; check mode: the bound's time " << g_cullBoundUs / (double)g_cullCycles
               << " us per cycle";
   }
   // check mode (VPX_FINDER_CHECK=1), after the grouping: every moving part the cull left out of a ball's list gets the full test (ReachEnergy's: its gap now less its move at the
   // ball's group feed, or its region less R, within the ball's travel); one it takes is a miss (CullMiss)
   static void CullCheck(Cycle &cy)
   {
      if (g_cullCy != &cy || g_cullTravel.empty() || !GroupStep::FinderCheckEnv())
         return;
      for (size_t i = 0; i < cy.bodies.size() && i < g_cullTravel.size(); i++)
      {
         const Body &w = cy.bodies[i];
         if (w.e == nullptr || g_cullTravel[i] < 0.0)
            continue;
         for (size_t k = 0; k < g_cullParts.size(); k++)
         {
            ContactHost *const o = g_cullParts[k];
            if (!o->m_enabled || !o->m_set.contact || std::any_of(w.e->items.begin(), w.e->items.end(), [o](const SurfaceContact &it) { return it.r == o; }))
               continue;
            ContactMovingPart *const mp = o->MovingPart();
            mp->SetBall(w.bs);
            const double gap = -mp->Probe(w.c).delta, full = std::max(gap - mp->EnergyMove(w.c, w.feedE, cy.dt, nullptr), mp->CycleRegionDist(w.c, cy.dt, w.feedE) - w.bs.R);
            if (full <= w.reachS)
               CullMiss("part", o->Name(), g_cullPartTold[k], w, full, mp->CullDist(w.c) - w.bs.R, g_cullTravel[i]);
         }
      }
   }

   // The balls held by a force-model contact
   static void GatherBalls(Cycle &cy)
   {
      CPROF(GatherBalls);
      ContactSurfaceHost::FinderBuild();   // shared_finder: the index over every static mesh, once before the first cycle (P1: only check mode reads it)
      ContactSurfaceHost::FinderTotals(false);   // shared_finder P2a: counts the cycle; the query's totals once a second of table time
      const double dt = cy.dt;
      const size_t first = cy.bodies.size();
      for (auto eit = g_engaged.begin(); eit != g_engaged.end();)
      {
         BallContact &e = *eit++;   // (a release below erases only this entry)
         if (e.ball->m_sleepGroup != 0)   // pile sleep: held where it is, its contact state kept (SleepInStep)
            continue;
         ContactHost *const o = e.host;   // null for a ball only a band holds
         if (o != nullptr && !o->m_set.contact)   // switched off in the ini while a ball is engaged
         {
            e.host = nullptr;   // one entry per ball: the entry stays; the gather below drops the host's items (a host switched off is no candidate)
         }
         {
            HitBall *const pb = e.ball;
            if (KickerHold::Held(pb) || pb->m_d.m_lockedInKicker)   // held by a kicker (KickerDesign 1h: no surface contacts): the kicker holds it, not the contact
            {
               if (!e.items.empty())
                  g_qe.cleared++;
               e.items.clear();   // one entry per ball: its items cleared, its entry kept
               continue;
            }
            e.lastHosts.clear();   // (counters: the items before this cycle's gather)
            for (const SurfaceContact &it : e.items)
               e.lastHosts.push_back(it.r);
            Body w;
            w.o = e.host;
            w.e = &e;
            LoadBody(w, pb);
            LoadState(w);
            w.band = e.band;
            Sources::Begin(w, dt);   // (it reads only the ball's floor items, which Regather part 1 keeps as they are)
            cy.bodies.push_back(std::move(w));
         }
      }
      // the part cull (runs/2026-11-20-part-cull/design.md): every body's travel bound from all the energy on the table first (it needs every body), then each ball's moving parts
      // (Regather part 1; the static mesh query waits for the ball's energy travel: GatherStatic) and the refresh
      CullPrepare(cy, first);
      for (size_t k = first; k < cy.bodies.size(); k++)
      {
         Body &w = cy.bodies[k];
         BallContact &e = *w.e;
         {
            CPROF(Regather);
            Regather(e, w.posStart, w.pb->m_d.m_radius, 1, g_cullTravel.empty() ? -1.0 : g_cullTravel[k]);   // (no bound on the table's energy: every part, as a new ball's first query)
         }
         {
            CPROF(RefreshCycle);
            for (SurfaceContact &it : e.items)
               if (it.r->FinderObjIndex() < 0)   // (the static mesh items: GatherStatic, after the grouping)
                  it.r->Refresh(it, w.bs.R, 2.0 * w.v.Length() * dt + 0.005, w.c);   // a surface: the faces the ball can reach this cycle
         }
      }
   }

   // the energy reach (runs/2026-11-14-reach-p1b/brief.md; design.md section 3, "result per ball"): after the grouping (ReachEnergy: each ball's travel from its energy, grown
   // to a fixed point by what it can reach and be fed by), each engaged ball's one query of the cycle, out to that travel: its static mesh surfaces come from it, before any solve.
   // Nothing joins later except by the inside check (a ball outside its query: a bound too small, logged as a fault)
   // the energy gather (runs/2026-11-15-finder-P6b3/brief_B2-1.md): after the grouping, a ball's non-mesh surfaces are those its travel this cycle (reachS) can reach.
   // A moving part (a grouping candidate when the part cull leaves it in, Regather part 1): kept when grouped with the ball. A fixed analytic surface: its bounding box less R, a
   // true lower bound of its distance. A band joins only by the band model (kept, never added). No margin, no rescan: every cycle from the grouping.
   static void GatherNonMesh(Body &w, const double T)
   {
      BallContact &e = *w.e;
      const double reach = std::max(w.reachS, 0.0);
      const Vertex3Ds &pos = w.posStart;
      const float boxU = w.pb->m_d.m_radius + 20.0f + 10.0f;   // (counters only: the old box)
      static std::map<std::string, std::pair<uint64_t, uint64_t>> s_kind;   // per host kind: kept that the box would not hold; held by the box that cannot be reached (diagnostics)
      static uint64_t s_calls = 0;
      std::vector<SurfaceContact> next;
      next.reserve(e.items.size());
      for (ContactHost *const o : ContactHost::Registry())
      {
         if (o->FinderObjIndex() >= 0)
            continue;   // a static mesh host: the query's (Regather part 2)
         const auto old = std::find_if(e.items.begin(), e.items.end(), [o](const SurfaceContact &it) { return it.r == o; });
         const FRect3D &bb = o->m_hitBBox;
         const bool usable = o->m_enabled && o->m_set.contact && (old != e.items.end() || !o->m_set.band);
         const bool box = usable && !(pos.x < bb.left - boxU || pos.x > bb.right + boxU || pos.y < bb.top - boxU || pos.y > bb.bottom + boxU || pos.z < bb.zlow - boxU || pos.z > bb.zhigh + boxU);
         bool keep = false;
         if (usable)
         {
            if (ContactMovingPart *const mp = o->MovingPart())
            {
               if (old != e.items.end())
               {
                  // the part cull (design 3): a part is kept when it is grouped with the ball (not partOut); ReachEnergy's test is not repeated: a part left out failed it in its last
                  // round at a feed at least the ball's group's (feedE) and the same travel, so it fails here too. Check mode: that test beside it, a difference counted
                  mp->SetBall(w.bs);
                  keep = !old->partOut;
                  if (GroupStep::FinderCheckEnv() && keep != (keep || std::max(old->reachGap - mp->EnergyMove(w.c, w.feedE, T, nullptr), mp->CycleRegionDist(w.c, T, w.feedE) - w.bs.R) <= reach))
                  {
                     if (g_cullGatherDiff++ == 0)
                        PLOGW << "Part cull: the energy gather's repeated part test keeps " << o->Name() << " (left out of the ball's group) at t=" << g_pplayer->m_time_msec
                              << " ms: the repeat is not a copy of ReachEnergy's test (check mode; first of this table start)";
                  }
               }
            }
            else
            {
               const double dx = std::max({ (double)bb.left - pos.x, 0.0, (double)pos.x - bb.right }), dy = std::max({ (double)bb.top - pos.y, 0.0, (double)pos.y - bb.bottom }),
                            dz = std::max({ (double)bb.zlow - pos.z, 0.0, (double)pos.z - bb.zhigh });
               keep = std::sqrt(dx * dx + dy * dy + dz * dz) * U - w.bs.R <= reach;
            }
         }
         if (keep != box)
         {
            auto &k = s_kind[std::string(o->MovingPart() != nullptr ? "part " : "fixed ") + typeid(*o).name()];
            (keep ? k.first : k.second)++;
         }
         if (!keep)
            continue;
         if (old != e.items.end())
         {
            next.push_back(std::move(*old));
            continue;
         }
         SurfaceContact it;
         it.r = o;
         it.rp = o->m_set.ToParams();
         it.hLayer = o->LayerM();
         it.st.Reset();
         next.push_back(std::move(it));
      }
      for (SurfaceContact &it : e.items)
         if (it.r != nullptr && it.r->FinderObjIndex() >= 0)
            next.push_back(std::move(it));   // the static mesh items as Regather part 2 left them
      std::stable_sort(next.begin(), next.end(), [](const SurfaceContact &a, const SurfaceContact &b) { return static_cast<const ContactHost *>(a.r)->m_regOrder < static_cast<const ContactHost *>(b.r)->m_regOrder; });
      e.items = std::move(next);   // one entry per ball: an empty set is a ball touching nothing
      if (++s_calls % 10000 == 0)
      {
         std::ostringstream o;
         o << "Energy gather: " << s_calls << " ball queries; per kind (added beyond the box / box items out of reach):";
         for (const auto &k : s_kind)
            o << " [" << k.first << " " << k.second.first << " / " << k.second.second << "]";
         PLOGI << o.str();
      }
   }

   static void GatherStatic(Cycle &cy)
   {
      for (Body &w : cy.bodies)
      {
         if (w.e == nullptr)
            continue;
         HitBall *const pb = w.pb;
         BallContact &e = *w.e;
         // the query's centre is where the step starts the ball: its owned state
         ContactSurfaceHost::FinderQuery(pb, w.c, w.v, w.bs.R, cy.dt, std::max(w.reachS, 0.0));
         e.reachS = w.reachS;   // (the stock net's log reads it)
         e.reachMs = g_pplayer->m_time_msec;
         {
            CPROF(Regather);
            Regather(e, w.posStart, pb->m_d.m_radius, 2);
            GatherNonMesh(w, cy.dt);
         }
         ContactSurfaceHost::FinderBind(pb, e.items, true);
         {
            CPROF(RefreshCycle);
            for (SurfaceContact &it : e.items)
               if (it.r->FinderObjIndex() >= 0)
                  it.r->Refresh(it, w.bs.R, 2.0 * w.v.Length() * cy.dt + 0.005, w.c);
         }
         ContactSurfaceHost::FinderCheckGather(pb, e.items, w.c, w.bs.R);   // check mode only (VPX_FINDER_CHECK=1)
         // (counters: items joined and left since the last cycle's gather)
         for (const SurfaceContact &it : e.items)
            if (std::find(e.lastHosts.begin(), e.lastHosts.end(), it.r) == e.lastHosts.end())
               g_qe.joined++;
         for (const void *const h : e.lastHosts)
            if (std::none_of(e.items.begin(), e.items.end(), [h](const SurfaceContact &it) { return it.r == h; }))
               g_qe.left++;
      }
   }


   // island reach: a ball and a moving part share an island only if the ball can reach the part's surface within this cycle (the reachable set over the cycle, runs/2026-11-05-island-reach/design.md).
   // A ball's travel: the speeds of every ball of its group (a hit hands speed on), what their forces add, and twice the speed of each part in the group (a ball leaves a moving surface at up to twice
   // its speed). A part's travel: its whole range of motion once a ball can reach it (a hit can move it anywhere it can go), else what its own drives give it (SurfaceSpeedBound). Grown to a fixed point:
   // the bounds only grow. A part no ball can reach is left out of every ball's solve (partOut) and moves alone; one a ball reaches is in that ball's island only. The gap is to its real surface (Probe).

   // ball pairs by a sorted list (to-do 1.34 P8; finder design.md section 3): the bodies kept sorted along the table's x axis by their centre at the cycle start, from last cycle's
   // order (insertion sort: the list is nearly sorted). Two balls can touch in the cycle only if their x intervals [x - R - r, x + R + r] overlap (r: each one's travel over the cycle):
   // |dx| <= |ci - cj|, so gap = |ci - cj| - Ri - Rj <= ri + rj needs |dx| <= Ri + ri + Rj + rj. The list only finds the candidates; the pair's own test decides.
   // Diagnostic (VPX_PAIR_CHECK=1): the old every-pair tests run beside the candidate ones and count each pair the candidates missed (logged with the running counters)
   inline static std::vector<const HitBall *> s_ballOrderLast;
   static inline uint64_t s_pairCheckTests = 0, s_pairCheckMissed = 0;
   static bool PairCheckOn() { static const bool on = [] { const char *const e = std::getenv("VPX_PAIR_CHECK"); return e != nullptr && e[0] == '1'; }(); return on; }
   static void PairCheckMiss(const Cycle &cy, const int i, const int j, const char *const what, const double gap, const double reach)
   {
      s_pairCheckMissed++;
      static int lines = 0;
      if (lines++ >= 30)
         return;
      const char *const ni = TestDriver::BallName(cy.bodies[i].pb), *const nj = TestDriver::BallName(cy.bodies[j].pb);
      PLOGW << "Pair check: " << what << " pair " << (ni != nullptr ? ni : "?") << " and " << (nj != nullptr ? nj : "?") << " gap " << gap * 1e3 << " mm, reach " << reach * 1e3
            << " mm, not among the sorted list's candidates (t " << g_pplayer->m_time_msec << " ms)";
   }
   static void SortBallsX(Cycle &cy)
   {
      const int nw = (int)cy.bodies.size();
      std::vector<int> &ord = cy.ballOrder;
      ord.clear();
      std::vector<char> in((size_t)nw, 0);
      for (const HitBall *const pb : s_ballOrderLast)   // last cycle's order, the balls still bodies; then the new ones
         for (int i = 0; i < nw; i++)
            if (!in[(size_t)i] && cy.bodies[i].pb == pb)
            {
               ord.push_back(i);
               in[(size_t)i] = 1;
               break;
            }
      for (int i = 0; i < nw; i++)
         if (!in[(size_t)i])
            ord.push_back(i);
      for (int a = 1; a < nw; a++)   // insertion sort by x (ties keep their order)
      {
         const int k = ord[(size_t)a];
         int b = a - 1;
         for (; b >= 0 && cy.bodies[ord[(size_t)b]].c.x > cy.bodies[k].c.x; b--)
            ord[(size_t)b + 1] = ord[(size_t)b];
         ord[(size_t)b + 1] = k;
      }
      s_ballOrderLast.resize((size_t)nw);
      for (int a = 0; a < nw; a++)
         s_ballOrderLast[(size_t)a] = cy.bodies[ord[(size_t)a]].pb;
   }
   // the candidate pairs (i < j), sorted in i, j order (the order the every-pair loops used): x intervals overlap with r[i] each ball's travel (null: 0)
   static void BallPairCandidates(const Cycle &cy, const double *const r, std::vector<std::pair<int, int>> &out)
   {
      out.clear();
      const std::vector<int> &ord = cy.ballOrder;
      const int nw = (int)ord.size();
      const auto hw = [&](const int i) { return cy.bodies[i].bs.R + (r != nullptr ? r[i] : 0.0); };
      double hwMax = 0.0;
      for (int i = 0; i < nw; i++)
         hwMax = std::max(hwMax, hw(i));
      for (int a = 0; a < nw; a++)
      {
         const int i = ord[(size_t)a];
         const double hi = hw(i);
         for (int b = a + 1; b < nw; b++)
         {
            const int j = ord[(size_t)b];
            const double dx = cy.bodies[j].c.x - cy.bodies[i].c.x;   // (>= 0: sorted)
            if (dx > hi + hwMax)
               break;   // every later ball is farther along x
            if (dx <= hi + hw(j))
               out.emplace_back(std::min(i, j), std::max(i, j));
         }
      }
      std::sort(out.begin(), out.end());
   }

   // ball pair reach rule 8: two balls in different islands last cycle that overlap at this cycle's start were missed by the bound
   inline static std::vector<std::pair<const HitBall *, uintptr_t>> s_ballGroupLast;   // (ball, its group) of the last cycle: its island's root ball, or its sleep group (held balls)
   static void BallPairMissCheck(const Cycle &cy, const std::vector<double> &bgap)
   {
      const int nw = (int)cy.bodies.size();
      const auto rootOf = [](const HitBall *pb) -> uintptr_t { for (const auto &g : s_ballGroupLast) if (g.first == pb) return g.second; return 0; };
      for (const auto &pr : cy.near0)   // (an overlap is among the touch candidates)
         if (const int i = pr.first, j = pr.second; bgap[(size_t)i * nw + j] < 0.0)
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
   // ball pair reach diagnostic: a pair the energy test groups (not in contact at the cycle start), logged when that grouping starts, with the terms of each V (200 lines at most)
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

   // PartsInSolve counters (deterministic, diagnostics only; printed on the "Energy reach rounds" line): parts in island solves summed over the island solves, of them at rest,
   // a part at rest PartsInSolve left out whose ball is inside it at the next cycle's start (a miss: rule 8, WARN once per table start), and the release test's outcomes by reason
   static inline uint64_t s_solveParts = 0, s_solveIslands = 0, s_solveRest = 0, s_solveMiss = 0, s_relClear = 0, s_relTimeout = 0, s_relHeld = 0;
   struct Edge { int i, p; SurfaceContact *it; double gap; bool on; };   // island reach: a ball and a moving part in its list: the gap from its probe now; on = grouped
   // the energy reach (runs/2026-11-14-reach-p1/brief.md; design_reach.md sections 1, 4a, 4b): ReachLinks with the one energy reach (rubbercontact.h EnergyReach). Grown to a
   // fixed point (the bounds only grow), each round from the groups as they are:
   //  - a group's energy: its balls' kinetic energy (spin included) and what their contacts hold (the stored energy at the last cycle's end, endStored, and the ball pairs' now);
   //    a held ball's given motion as a speed (a ball leaves it at up to twice its speed); its balls' forces (gravity and nudge, magnets) as the drive
   //  - each grouped part's work: its kinetic energy and its drives' work over its move, the move from its group's balls' energy (EnergyMove); the group's feed energy is the
   //    group's energy plus its parts' work, all given to one ball: its travel EnergyReach(0, feed, m, forces, T) plus the held balls' given speed x T
   //  - a ball joins a part when gap <= its travel + the part's move at the ball (the feed energy of the ball's group and of the part's);
   //    two balls when gap <= the sum of their travels (ball pair reach)
   // Each body keeps its travel and its group's feed energy (Body::reachS, feedE: the group step's cycle turns and lists); the ball's one static mesh query of the cycle reaches its
   // travel (GatherStatic, after this). Static surfaces feed no energy (what a contact holds is in endStored), so they do not change the fixed point
   // NOTE: the fixed point stops at the first round that adds nothing: each round adds a part edge or ball pair or stops (no round cap)
   static void ReachEnergy(Cycle &cy, std::vector<Edge> &edges, const std::vector<ContactMovingPart *> &parts, const bool ballPairs)
   {
      const int nw = (int)cy.bodies.size(), np = (int)parts.size();
      const double T = cy.dt;
      if (nw == 0)
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
      std::vector<int> partBall(np, -1);   // the first ball that can reach it
      std::vector<char> paired;
      std::vector<double> bgap;
      std::vector<std::pair<int, int>> cand;   // a round's candidate ball pairs (BallPairCandidates)
      if (ballPairs)
      {
         paired.assign((size_t)nw * nw, 0);
         bgap.assign((size_t)nw * nw, std::numeric_limits<double>::quiet_NaN());   // (each pair's gap, measured when it is first a candidate)
         for (const auto &pr : cy.pairs)
            paired[(size_t)pr.first * nw + pr.second] = 1;
         for (const auto &pr : cy.near0)
            bgap[(size_t)pr.first * nw + pr.second] = (cy.bodies[pr.first].c - cy.bodies[pr.second].c).Length() - cy.bodies[pr.first].bs.R - cy.bodies[pr.second].bs.R;
         BallPairMissCheck(cy, bgap);
         if (PairCheckOn())   // (diagnostic: every overlapping pair among the candidates)
            for (int i = 0; i < nw; i++)
               for (int j = i + 1; j < nw; j++)
               {
                  const double gap = (cy.bodies[i].c - cy.bodies[j].c).Length() - cy.bodies[i].bs.R - cy.bodies[j].bs.R;
                  s_pairCheckTests++;
                  if (gap < 0.0 && std::isnan(bgap[(size_t)i * nw + j]))
                     PairCheckMiss(cy, i, j, "overlapping (rule 8 check)", gap, 0.0);
               }
      }
      std::vector<double> pairHeld(cy.pairs.size(), 0.0);   // what each ball pair holds now (its stored energy)
      for (size_t k = 0; k < cy.pairs.size(); k++)
      {
         const Body &A = cy.bodies[cy.pairs[k].first], &B = cy.bodies[cy.pairs[k].second];
         if (const BallPairState *const ps = FindPair(A.pb, B.pb))
            pairHeld[k] = std::max(PairHolds(cy, A, B, *ps), 0.0);
      }
      const size_t nPairs0 = cy.pairs.size();
      std::vector<double> eb(nw), af(nw), hv(nw), ef(nw), sb(nw), wp(np), gw(nw), gh(nw);
      std::vector<int> np2(nw);   // the full feed: the grouped parts of each group (by its root ball); gw: their work and squash
      int rounds = 0;
      // the one-pass energy bound (runs/2026-11-20-energy-bound/design.md; the cull's law, CullPrepare): each part's work is its EnergyCap, from its state now with no feed, so a group's feed
      // is its balls' energy + its parts' caps and squash, with its stops' squash closed in one pass (OnePassE: its parts' aLin and aHertz): no sweeps. Caps once per call
      static std::vector<double> s_cap, s_capA, s_capH;   // (each part's cap, aLin and aHertz this call; grown only, no allocation once sized)
      if (s_cap.size() < (size_t)np)
      {
         s_cap.resize(np);
         s_capA.resize(np);
         s_capH.resize(np);
      }
      for (int p = 0; p < np; p++)
      {
         s_capA[p] = s_capH[p] = 0.0;
         s_cap[p] = parts[p]->EnergyCap(T, s_capA[p], s_capH[p]);
      }
      for (bool grew = true; grew;)
      {
         grew = false;
         rounds++;
         std::fill(eb.begin(), eb.end(), 0.0);
         std::fill(af.begin(), af.end(), 0.0);
         std::fill(hv.begin(), hv.end(), 0.0);
         for (int i = 0; i < nw; i++)
         {
            const Body &w = cy.bodies[i];
            const int g = find(i);
            if (w.held)
               hv[g] += 2.0 * (w.v.Length() + w.heldA.Length() * T);
            else
               eb[g] += 0.5 * w.bs.m * w.v.Dot(w.v) + 0.5 * w.bs.Inertia() * w.w.Dot(w.w) + (w.e != nullptr ? std::max(w.e->endStored, 0.0) : 0.0);
            af[g] += w.gravAcc + w.magReach;
            if (w.band != nullptr && w.band->CordOn())   // a band that holds the ball: what its cord holds (CordBand::heldJ), once per group
            {
               bool first = true;
               for (int j = 0; j < i && first; j++)
                  first = !HoldsBand(cy.bodies[j], w.band);
               if (first)   // (and its sling coil's work bound for this cycle: CordFixed::Part::CoilWorkBound)
                  eb[g] += w.band->Cord()->heldJ + w.band->Cord()->fixed->CoilWorkBound(T);
            }
            if (w.e != nullptr)
               for (ContactSurfaceHost *const b : w.e->bandsMore)   // each other band that holds it, likewise once per group
                  if (b->CordOn())
                  {
                     bool first = true;
                     for (int j = 0; j < i && first; j++)
                        first = !HoldsBand(cy.bodies[j], b);
                     if (first)
                        eb[g] += b->Cord()->heldJ + b->Cord()->fixed->CoilWorkBound(T);
                  }
         }
         for (size_t k = 0; k < nPairs0; k++)
            eb[find(cy.pairs[k].first)] += pairHeld[k];
         ef = eb;
         for (int p = 0; p < np; p++)
            wp[p] = partBall[p] >= 0 ? s_cap[p] : 0.0;
         for (int p = 0; p < np; p++)
            if (partBall[p] >= 0)
            {
               auto &k = g_partWorkMax[typeid(*parts[p]).name()];   // (diagnostics only)
               if (wp[p] > k.first)
                  k = { wp[p], g_pplayer->m_time_msec };
            }
         for (int p = 0; p < np; p++)
            if (partBall[p] >= 0)
               ef[find(partBall[p])] += wp[p] + parts[p]->squashJ;   // (its kinetic energy is in its work; what its contacts held at the last cycle's end, MovingPartState::squashJ)
         {   // (each group's stops' squash past their touch: E = E0 + aLin sqrt(E) + aHertz E^0.4 over its grouped parts, closed in one pass; no bound: INFINITY)
            std::fill(gw.begin(), gw.end(), 0.0);
            std::fill(gh.begin(), gh.end(), 0.0);
            for (int p = 0; p < np; p++)
               if (partBall[p] >= 0)
               {
                  gw[find(partBall[p])] += s_capA[p];
                  gh[find(partBall[p])] += s_capH[p];
               }
            for (int i = 0; i < nw; i++)
               if (find(i) == i && (gw[i] > 0.0 || gh[i] > 0.0))
                  ef[i] = OnePassE(ef[i], gw[i], gh[i]);
         }
         for (int i = 0; i < nw; i++)
         {
            const int g = find(i);
            sb[i] = TravelBound(ef[g], cy.bodies[i].bs.m, af[g], hv[g], T);
         }
         for (Edge &ed : edges)
         {
            if (ed.on)
               continue;
            const int gi = find(ed.i), pb = partBall[ed.p];
            const double E = ef[gi] + (pb >= 0 && find(pb) != gi ? ef[find(pb)] : 0.0);
            const Body &w = cy.bodies[ed.i];
            if (!std::isfinite(E) || std::max(ed.gap - parts[ed.p]->EnergyMove(w.c, E, T, nullptr), parts[ed.p]->CycleRegionDist(w.c, T, E) - w.bs.R) <= sb[ed.i])   // (no bound: joined)
            {
               ed.on = grew = true;
               if (pb < 0)
                  partBall[ed.p] = ed.i;
               else if (find(ed.i) != find(pb))
                  unite(ed.i, pb);
            }
         }
         if (ballPairs)
         {
            BallPairCandidates(cy, sb.data(), cand);
            if (PairCheckOn())   // (diagnostic: the every-pair energy test, before this round joins anything; sb is fixed for the round)
               for (int i = 0; i < nw; i++)
                  for (int j = i + 1; j < nw; j++)
                  {
                     const double gap = (cy.bodies[i].c - cy.bodies[j].c).Length() - cy.bodies[i].bs.R - cy.bodies[j].bs.R;
                     s_pairCheckTests++;
                     if (!paired[(size_t)i * nw + j] && gap <= sb[i] + sb[j] && !std::binary_search(cand.begin(), cand.end(), std::make_pair(i, j)))
                        PairCheckMiss(cy, i, j, "energy", gap, sb[i] + sb[j]);
                  }
            for (const auto &pr : cand)
               if (const int i = pr.first, j = pr.second; !paired[(size_t)i * nw + j])
               {
                  double &gij = bgap[(size_t)i * nw + j];
                  if (std::isnan(gij))
                     gij = (cy.bodies[i].c - cy.bodies[j].c).Length() - cy.bodies[i].bs.R - cy.bodies[j].bs.R;
                  if (gij <= sb[i] + sb[j])
                  {
                     BallPairReachLog(cy, i, j, bgap[(size_t)i * nw + j], [&](const int k) {   // (the terms of the travel, read before the pair joins the groups)
                        const int g = find(k);
                        std::ostringstream o;
                        o << "travel " << sb[k] * 1e3 << " mm = energy " << eb[g] * 1e3 << " mJ + parts' work " << (ef[g] - eb[g]) * 1e3 << " mJ, forces " << af[g] << " m/s2, held " << hv[g] << " m/s, |v| " << cy.bodies[k].v.Length();
                        return o.str(); });
                     paired[(size_t)i * nw + j] = 1;
                     cy.pairs.emplace_back(i, j);
                     unite(i, j);
                     grew = true;
                  }
               }
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
      for (int i = 0; i < nw; i++)
      {
         Body &w = cy.bodies[i];
         w.reachS = sb[i];
         w.feedE = ef[find(i)];
         double vm = 0.0;   // energy skip test: the top speed (rubbercontact.h EnergyReach) and the group's energy at it, Etot = feed + the forces' work over the cycle
         EnergyReach(0.0, ef[find(i)], w.bs.m, af[find(i)], T, &vm);
         w.vTop = vm + hv[find(i)];
         w.eTop = 0.5 * w.bs.m * vm * vm;
         w.eGroup = find(i);
      }
      if (RelookTraceOn())   // VPX_RELOOK_TRACE: each ball's bound term by term, from this call's last round (the finder prints it at a re-look)
      {
         const auto pname = [&](const int p) {   // (its kind, origin in mm and state note)
            const V3 a = parts[p]->Origin();
            std::ostringstream n;
            n << typeid(*parts[p]).name() << " at (" << a.x * 1e3 << ", " << a.y * 1e3 << ", " << a.z * 1e3 << ") mm " << parts[p]->StateNote();
            return n.str(); };
         for (int i = 0; i < nw; i++)
         {
            const Body &w = cy.bodies[i];
            if (w.pb == nullptr)
               continue;
            const int g = find(i);
            int nb = 0;
            for (int j = 0; j < nw; j++)
               nb += find(j) == g ? 1 : 0;
            std::ostringstream o;
            o << "T " << T * 1e3 << " ms, m " << w.bs.m << " kg, |v| " << w.v.Length() << " m/s, own KE " << 0.5 * w.bs.m * w.v.Dot(w.v) * 1e3 << " mJ, spin " << 0.5 * w.bs.Inertia() * w.w.Dot(w.w) * 1e3
              << " mJ, endStored " << (w.e != nullptr ? w.e->endStored * 1e3 : 0.0) << " mJ, held " << w.held << "; group balls " << nb << ", E0 (balls, stored, pairs, bands) " << eb[g] * 1e3 << " mJ; parts in E1:";
            for (int p = 0; p < np; p++)
               if (partBall[p] >= 0 && find(partBall[p]) == g)
                  o << " [" << pname(p) << " work " << wp[p] * 1e3 << " mJ, squash " << parts[p]->squashJ * 1e3 << " mJ, aLin " << s_capA[p] << ", aHertz " << s_capH[p] << "]";
            o << "; aLin sum " << gw[g] << " J^0.5, aHertz sum " << gh[g] << " J^0.6; E1 " << ef[g] * 1e3 << " mJ; travel inputs: v from E1 " << std::sqrt(2.0 * std::max(ef[g], 0.0) / w.bs.m)
              << " m/s, acceleration " << af[g] << " m/s2 (own gravity " << w.gravAcc << ", magnets " << w.magReach << "), hv " << hv[g] << " m/s; travel " << sb[i] * 1e3 << " mm; parts left out:";
            for (int p = 0; p < np; p++)
               if (partBall[p] < 0 || find(partBall[p]) != g)
               {
                  o << " [" << pname(p) << (partBall[p] < 0 ? " (no ball)" : " (other group)");
                  o << " cap " << s_cap[p] * 1e3 << " mJ";
                  for (const Edge &ed : edges)
                     if (ed.i == i && ed.p == p)
                     {
                        const int pb = partBall[p];
                        const double E = ef[g] + (pb >= 0 && find(pb) != g ? ef[find(pb)] : 0.0);
                        o << ", in list: gap " << ed.gap * 1e3 << " mm, part move at E " << parts[p]->EnergyMove(w.c, E, T, nullptr) * 1e3 << " mm, region dist " << (parts[p]->CycleRegionDist(w.c, T, E) - w.bs.R) * 1e3 << " mm";
                     }
                  o << "]";
               }
            RelookTerms()[w.pb] = o.str();
         }
      }
      for (int p = 0; p < np; p++)   // energy skip test: what each part can hold before any ball touches it: its cap (EnergyCap: its own kinetic energy and drive work) and its contacts' squash, its stops' squash closed in one pass
      {
         parts[p]->eOwnJ = OnePassE(s_cap[p] + parts[p]->squashJ, s_capA[p], s_capH[p]);
      }
      static uint64_t nCyc = 0, nRounds = 0;   // the fixed point's rounds per cycle (diagnostics: logged every 10000 cycles)
      static int worst = 0;
      static uint64_t nMulti = 0, nGrown = 0;   // the full feed diagnostics: groups with two or more grouped parts, part moves fed by another part's work or squash
      std::fill(np2.begin(), np2.end(), 0);
      std::fill(gw.begin(), gw.end(), 0.0);
      for (int p = 0; p < np; p++)
         if (partBall[p] >= 0)
         {
            np2[find(partBall[p])]++;
            gw[find(partBall[p])] += wp[p] + parts[p]->squashJ;
         }
      for (int i = 0; i < nw; i++)
         nMulti += np2[i] >= 2 ? 1 : 0;
      for (int p = 0; p < np; p++)
         if (partBall[p] >= 0 && np2[find(partBall[p])] >= 2 && gw[find(partBall[p])] - wp[p] - parts[p]->squashJ > 0.0)
            nGrown++;
      nCyc++;
      nRounds += (uint64_t)rounds;
      worst = std::max(worst, rounds);
      if (nCyc % 10000 == 0)
         PLOGI << "Energy reach rounds: " << nCyc << " cycles, rounds per cycle mean " << (double)nRounds / (double)nCyc << ", worst " << worst << "; groups with two or more parts " << nMulti << ", part moves fed by another part " << nGrown << "; parts in solve " << s_solveParts << " over " << s_solveIslands << " island solves (at rest " << s_solveRest << "), PartsInSolve misses " << s_solveMiss << ", releases clear " << s_relClear << ", timeout " << s_relTimeout << ", kept within the surface.s bound " << s_relHeld;
   }

   static void ReachLinks(Cycle &cy)
   {
      const int nw = (int)cy.bodies.size();
      const double T = cy.dt;
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
                  if (it.solveOutLast && gap < 0.0)   // PartsInSolve left this part at rest out last cycle, and the ball is inside it now: a missed contact (rule 8)
                  {
                     static const void *owner = nullptr;
                     if (owner != g_pplayer)
                     {
                        owner = g_pplayer;
                        const char *const name = TestDriver::BallName(w.pb);
                        PLOGW << "PartsInSolve: ball " << (name != nullptr ? name : "") << " is " << -gap * 1e3 << " mm inside " << it.r->Name() << ", which PartsInSolve left out of its solve last cycle (a missed contact; first of this table start)";
                     }
                     s_solveMiss++;
                  }
                  it.solveOutLast = false;
                  it.reachGap = gap;
                  const int p = (int)(std::find(parts.begin(), parts.end(), mp) - parts.begin());
                  if (p == (int)parts.size())
                     parts.push_back(mp);
                  edges.push_back({ i, p, &it, gap, false });
               }
      ReachEnergy(cy, edges, parts, nw >= 2);
   }

   // The energy reach of one ball alone: ReachLinks / ReachEnergy on a cycle of this one ball (its own kinetic energy and spin, what its contacts hold, gravity + nudge + magnets
   // (ReachAcc), the drive work of every moving part in its list; no ball pairs). cy gets the ball as its one body, its travel in cy.bodies[0].reachS. Used by a new ball's early
   // query (EngageNear) and a sleeping group's wake test (SleepParts)
   static void LoneReach(Cycle &cy, BallContact &e, HitBall *const pb, const double dt)
   {
      cy.dt = dt;
      cy.cs = &ContactSettings::Get();
      Body w;
      w.o = e.host;
      w.e = &e;
      LoadBody(w, pb);
      LoadState(w);
      w.band = e.band;
      Sources::ReachAcc(w);
      cy.bodies.push_back(std::move(w));
      ReachLinks(cy);
   }

   // The cycle's island reach marks end with it (other searches between cycles see every part); a part left out is checked at the next cycle's start (ReachLinks)
   static void ReachClear()
   {
      for (BallContact &e : g_engaged)
         for (SurfaceContact &it : e.items)
         {
            it.partOutLast = it.partOut;
            it.partOut = false;
         }
   }

   // (runs/2026-11-15-finder-P6b3/paper_P7.md 2) a band is a part of every ball's energy query. How far each span of the cord can move in the coming cycle T:
   // its balls' deepest deflection now on that span plus the largest move of what can push that span, its balls' travels that reach it (SpanPush, m) and, on the plate's span, the
   // plate's energy move (PlateMoveM); the cord is massless, its span ends fixed on their posts and its stops hard, so it is nowhere else. The ball's distance less that move
   // (QueryGap) is what the ball's own travel must reach to touch it. m
   // whether a band holds the body's ball: its first band, or one of its others
   static bool HoldsBand(const Body &w, const ContactSurfaceHost *const b)
   {
      return w.band == b || (w.e != nullptr && !w.e->bandsMore.empty() && std::find(w.e->bandsMore.begin(), w.e->bandsMore.end(), b) != w.e->bandsMore.end());
   }
   static double BandQueryGap(ContactSurfaceHost *const band, const V3 &c, const double R, const std::vector<double> &dNow, const std::vector<double> &push, const double T)
   {
      const CordFixed::Part &p = *band->Cord()->fixed;
      return p.QueryGap(c, R, dNow, push, p.PlateMoveM(T));
   }
   // per span, the cord's deflection now and the travels of its balls that reach that span this cycle (CordFixed::Part::SpanPush), from the bodies it holds (a band-held ball never
   // sleeps, so all of them are bodies)
   static void BandHeld(const Cycle &cy, ContactSurfaceHost *const band, std::vector<double> &dNow, std::vector<double> &push)
   {
      dNow.clear();
      push.clear();
      for (const Body &w : cy.bodies)
         if (HoldsBand(w, band))
            band->Cord()->fixed->SpanPush(w.c, w.bs.R, w.reachS, dNow, push);
   }
   static inline uint64_t s_bandBallCycles = 0;   // ball-cycles held by a band (FindIslands), either key: the take rule's cost (diagnostics)
   static inline uint64_t s_bandSecond = 0;   // takes of a ball by a second band while another holds it (BandQuery; diagnostics)
   // after the energy reach (ReachLinks), each ball's band membership is its query's: a ball whose travel (reachS) reaches the band (BandQueryGap <= reachS, the grouping's
   // own test) is taken, a held one that no longer reaches it is let go (no hit: BandHit fires at the first touch, ContactSurfaceHost::CordFinish). A take links the band's balls and feeds
   // the cord's energy to the group (ReachEnergy), so the reach is grown again until no ball joins: a fixed point (the releases on the first pass only, so it ends). Returns whether the
   // links changed
   static bool BandQuery(Cycle &cy)
   {
      bool any = false;
      for (int pass = 0; pass < 16; pass++)
      {
         bool changed = false;
         // the part cull (design 2R, bands): the table's cord bands only (CullLists, the registry's order); a band no body holds and whose outline box no body's travel bound
         // reaches (CullBandDist less R) is skipped this pass (check mode: the full test for each body, a take = a miss). A cycle with no travel bounds (a new ball's own cycle, LoneReach;
         // no bound on the table's energy, CullPrepare) tests every host in the registry
         const bool cull = g_cullCy == &cy && g_cullTravel.size() == cy.bodies.size();
         const size_t nHosts = cull ? g_cullBands.size() : ContactHost::Registry().size();
         for (size_t oi = 0; oi < nHosts; oi++)
         {
            ContactHost *const o = cull ? g_cullBands[oi].host : ContactHost::Registry()[oi];
            ContactSurfaceHost *const band = o->m_set.band && o->m_enabled ? o->BandHost() : nullptr;
            if (band == nullptr || !band->CordModel())
               continue;
            band->CordEnsure();   // (the cord's geometry from the band's posts and plate, made once)
            if (!band->CordOn())
               continue;
            if (cull && g_cullBands[oi].box)
            {
               const CullBand &cb = g_cullBands[oi];
               bool reached = false;
               for (size_t i = 0; i < cy.bodies.size() && !reached; i++)
               {
                  const Body &w = cy.bodies[i];
                  if (HoldsBand(w, band))
                     reached = true;
                  else if (w.e != nullptr && g_cullTravel[i] >= 0.0)
                  {
                     g_cullBandTests++;
                     reached = CullBandDist(cb, w.c) - w.bs.R <= g_cullTravel[i];
                  }
                  else
                     reached = true;   // (a body with no bound: the full test)
               }
               if (!reached)
               {
                  g_cullBandSkipped++;
                  if (GroupStep::FinderCheckEnv())
                  {
                     static const std::vector<double> none;
                     for (size_t i = 0; i < cy.bodies.size(); i++)
                     {
                        const Body &w = cy.bodies[i];
                        if (w.e == nullptr || KickerHold::Held(w.pb) || w.pb->m_d.m_lockedInKicker)
                           continue;
                        const double full = BandQueryGap(band, w.c, w.bs.R, none, none, cy.dt);
                        if (full <= std::max(w.reachS, 0.0))
                           CullMiss("band", o->Name(), g_cullBandTold[oi], w, full, CullBandDist(cb, w.c) - w.bs.R, g_cullTravel[i]);
                     }
                  }
                  continue;
               }
            }
            std::vector<double> dNow, push;
            BandHeld(cy, band, dNow, push);
            for (Body &w : cy.bodies)
            {
               if (w.e == nullptr || KickerHold::Held(w.pb) || w.pb->m_d.m_lockedInKicker)
                  continue;
               const bool reaches = BandQueryGap(band, w.c, w.bs.R, dNow, push, cy.dt) <= std::max(w.reachS, 0.0);
               const bool holds = HoldsBand(w, band);
               if (holds && !reaches && pass == 0)
               {
                  band->BandRelease(w.pb, "clear");   // (another band that holds it becomes its first)
                  w.band = w.e->band;
                  changed = true;
               }
               else if (!holds && reaches)
               {
                  if (w.band == nullptr)
                  {
                     if (band->BandEngage(w.pb))
                     {
                        w.band = band;
                        changed = true;
                     }
                  }
                  else   // a ball within reach of a second cord: held by both (each cord steps it, CordIsland; CordFollow gives every cord the island's balls)
                  {
                     band->CordEnsure();
                     BandBall bb;
                     bb.ball = w.pb;
                     bb.vin = w.v;
                     band->m_engaged.push_back(bb);
                     w.e->bandsMore.push_back(band);
                     if (++s_bandSecond % 100 == 1)   // (counted: how often a ball is held by two bands at once)
                        PLOGI << "Band " << band->Name() << ": a ball held by band " << w.band->Name() << " is within reach of it too, held by both, " << s_bandSecond << " so far, t=" << g_pplayer->m_time_msec;
                     changed = true;
                  }
               }
            }
         }
         if (!changed)
            break;
         any = true;
         BandLinks(cy);
         ReachLinks(cy);
      }
      return any;
   }
   // a band's balls step together (CordIsland): each one linked, for each band that holds it, to the first body that band holds
   static void BandLinks(Cycle &cy)
   {
      const int nw = (int)cy.bodies.size();
      cy.links.clear();
      for (int i = 0; i < nw; i++)
         if (cy.bodies[i].band != nullptr)
         {
            for (int j = 0; j < i; j++)
               if (HoldsBand(cy.bodies[j], cy.bodies[i].band))
               {
                  cy.links.emplace_back(j, i);
                  break;
               }
            if (cy.bodies[i].e != nullptr)
               for (ContactSurfaceHost *const b : cy.bodies[i].e->bandsMore)   // (its other bands)
                  for (int j = 0; j < i; j++)
                     if (HoldsBand(cy.bodies[j], b))
                     {
                        cy.links.emplace_back(j, i);
                        break;
                     }
         }
   }

   // Bodies whose surfaces are within reach of each other (plus what they can close in this cycle) interact: they form islands
   static void FindIslands(Cycle &cy)
   {
      CPROF(FindIslands);
      const int nw = (int)cy.bodies.size();
      cy.pairs.clear();
      cy.near0.clear();
      // the seed pairs: the pairs in contact now (gap <= 0, or a live pair state: it may hold stored energy), so ReachEnergy reads what they hold before the reach is known;
      // the energy test (ReachEnergy) decides every other pair
      SortBallsX(cy);
      if (nw >= 2)
      {
         BallPairCandidates(cy, nullptr, cy.near0);
         for (const auto &pr : cy.near0)
            if ((cy.bodies[pr.first].c - cy.bodies[pr.second].c).Length() - cy.bodies[pr.first].bs.R - cy.bodies[pr.second].bs.R <= 0.0)
               cy.pairs.push_back(pr);
         for (const BallPairState &ps : g_ballPairs)   // a live pair state (its balls may be apart): both balls bodies this cycle
         {
            int ia = -1, ib = -1;
            for (int i = 0; i < nw; i++)
            {
               ia = cy.bodies[i].pb == ps.a ? i : ia;
               ib = cy.bodies[i].pb == ps.b ? i : ib;
            }
            if (ia >= 0 && ib >= 0)
               cy.pairs.emplace_back(std::min(ia, ib), std::max(ia, ib));
         }
         std::sort(cy.pairs.begin(), cy.pairs.end());   // (the every-pair loop's order)
         cy.pairs.erase(std::unique(cy.pairs.begin(), cy.pairs.end()), cy.pairs.end());
         if (PairCheckOn())   // (diagnostic: every pair in contact among the seeds)
            for (int i = 0; i < nw; i++)
               for (int j = i + 1; j < nw; j++)
               {
                  const double gap = (cy.bodies[i].c - cy.bodies[j].c).Length() - cy.bodies[i].bs.R - cy.bodies[j].bs.R;
                  s_pairCheckTests++;
                  if ((gap <= 0.0 || FindPair(cy.bodies[i].pb, cy.bodies[j].pb) != nullptr) && !std::binary_search(cy.pairs.begin(), cy.pairs.end(), std::make_pair(i, j)))
                     PairCheckMiss(cy, i, j, "seed (in contact)", gap, 0.0);
               }
      }
      BandLinks(cy);
      ReachLinks(cy);
      BandQuery(cy);   // (a band's balls by their query; the links and the reach again when they change)
      for (const Body &w : cy.bodies)
         s_bandBallCycles += w.band != nullptr ? 1 : 0;   // (the cost of the band's take rule: ball-cycles stepped with a cord, told with the running counters)
      cy.islands.Build(nw, cy.pairs, cy.links);
   }

   // State back into the balls, the engagement bookkeeping and the release decision
   static void WriteBack(Cycle &cy, const ContactSystem::Islands::Island &isl)
   {
      CPROF(WriteBack);
      const ContactProf::Stage profStage(ContactProf::StWriteBack);   // P3a profile: the calling stage of the checks here (diagnostics only)
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
         WriteBall(pb, w.c, w.v, w.w);
         const Vertex3Ds ds = pb->m_d.m_pos - w.posStart;
         pb->m_lastEventSqrDist += ds.LengthSquared();
#ifdef C_DYNAMIC
         pb->m_drsq = ds.LengthSquared();
#endif
         pb->m_rubberAdvanced = true;   // the ball has been moved for the whole cycle
         // release: a band lets go of a ball when the ball's query no longer reaches it (BandQuery, next cycle)
      }
   }

   // one entry per ball, rule 8: a ball found inside a surface where its group state is (re)built (zero = that depth, nothing stored): a placement fault or a contact missed earlier.
   // Shown once per surface per table start; every one is counted (g_qe.firstSightZero).
   static void FoundInside(const HitBall *const pb, const std::string &surface, const double depth)
   {
      static const void *owner = nullptr;
      static std::vector<std::string> seen;
      if (owner != g_pplayer)
      {
         owner = g_pplayer;
         seen.clear();
      }
      if (std::find(seen.begin(), seen.end(), surface) != seen.end())
         return;
      seen.push_back(surface);
      const char *const name = TestDriver::BallName(pb);
      const Vertex3Ds &p = pb->m_d.m_pos;
      std::ostringstream o;
      o << "Contact: ball " << (name != nullptr ? name : "") << " at (" << p.x << ", " << p.y << ", " << p.z << ") found " << depth * 1e3 << " mm inside " << surface
        << ": a ball started or was moved inside " << surface << ": a placement fault or a contact missed earlier. It counts from there, nothing stored. t=" << g_pplayer->m_time_msec;
      PLOGW << o.str();
      if (g_pplayer->m_liveUI != nullptr)
         g_pplayer->m_liveUI->PushNotification(o.str(), 6000);
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
      const ContactProf::Stage profStage(ContactProf::StHolds);   // P3a profile: the calling stage of the checks here (diagnostics only)
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
      // one entry per ball: a ball handed to the group step this cycle (a new entry in an island without a cord: RingHandIn runs in its CordIsland, after this check)
      auto handedIn = [&cy](const Body &w)
      {
         if (w.e == nullptr || w.e->groupLoop)
            return false;
         const int i = (int)(&w - &cy.bodies[0]);
         for (size_t k = 0; k < cy.islands.count; k++)
         {
            const std::vector<int> &m = cy.islands.list[k].members;
            if (std::find(m.begin(), m.end(), i) == m.end())
               continue;
            for (const int j : m)
               if (cy.bodies[j].band != nullptr)
                  return false;
            return true;
         }
         return false;
      };
      for (Body &w : cy.bodies)
      {
         double held = 0.0, worst = 0.0, depth0 = 0.0, depth1 = 0.0, heldEnd = 0.0;
         const SurfaceContact *where = nullptr;
         // one entry per ball: first sight per item (a surface new in the items), not per entry
         const bool handIn = handedIn(w);
         if (w.e != nullptr)
            for (SurfaceContact &it : w.e->items)
            {
               // one entry per ball, first sight of an item (new in the query) already overlapping: the zero rule (GroupStepDesign 14.4): its zero is its start depth, nothing stored,
               // no position snap (rule 10b: the overlap came from no drive). A moving part's surfaces take the same rule in the group step (SetShapeItem: not touching, zero 0)
               if (it.endHeld < 0.0 && w.e->groupLoop && !it.gs.touching && it.gs.zero <= 0.0 && it.r->MovingPart() == nullptr && !it.NoFaces(w.c))
               {
                  const TubeProbe pr = it.r->ProbePatch(w.bs, it, w.c);
                  if (pr.overlap && pr.delta > 0.0)
                  {
                     it.gs.zero = pr.delta;
                     g_qe.firstSightZero++;
                     FoundInside(w.pb, it.r->Name(), pr.delta);
                  }
               }
               // a ball handed in this cycle: an item RingHandIn will find inside takes zero = its depth there, so it holds nothing at this reading either (same rule, same probe)
               bool zeroAtHandIn = false;
               if (handIn && !(it.touching && w.e->dwellS > 0.0) && it.r->MovingPart() == nullptr && !it.NoFaces(w.c))
               {
                  const TubeProbe pr = it.r->ProbePatch(w.bs, it, w.c);
                  zeroAtHandIn = pr.overlap && pr.delta > 0.0;
               }
               double d = 0.0;
               const double v = w.e->seatPending || zeroAtHandIn ? 0.0 : ContactHolds(it, w.bs, w.c, w.e->groupLoop, d);   // (found inside and not seated: it counts from zero)
               held += v;
               heldEnd += std::max(it.endHeld, 0.0);
               const double was = it.endHeld >= 0.0 ? it.endDepth : 0.0;
               if (d - was > worst)
               {
                  worst = d - was;
                  where = &it;
                  depth0 = was;
                  depth1 = d;
               }
            }
         w.eStart = (w.held ? 0.0 : BodyMotion(w, g)) + held;   // (a held ball's motion and height are the kicker's: KickerDesign 1h)
         w.eScale = (w.held ? 0.0 : 0.5 * w.bs.m * w.v.Dot(w.v) + 0.5 * w.bs.Inertia() * w.w.Dot(w.w)) + held;
         if (w.e != nullptr)
         {
            const double rise = held - heldEnd;
            const bool first = where != nullptr && where->endHeld < 0.0;
            if (rise > std::max(1e-3, 0.01 * w.eScale))
               EnergyRaise(!first ? "between cycles" : "at first sight", w.pb, where != nullptr ? where->r->Name() : std::string("?"), depth0, depth1, rise);
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
         driven = driven || w.band != nullptr || w.magAcc > 0.0 || (w.held && w.v.Length() > 0.0)   // (a carried ball: the arm, a drive)
            || (w.extBaseVS.Length() > 0.0 && (w.extBaseVS * (1.0 / PHYSICS_STEPTIME_S) - g).Length() > 0.02 * g.Length());   // (the nudge)
         double held = 0.0;
         if (w.e != nullptr)
         {
            for (SurfaceContact &it : w.e->items)
            {
               driven = driven || (!it.partOut && it.r->MovingPart() != nullptr);
               const double v = ContactHolds(it, w.bs, w.c, w.e->groupLoop, it.endDepth);
               held += v;
               it.endHeld = v;
               if (it.partOut)   // (out of the ball's solve: not touching; its zero is 0, as any item clear at a cycle's end)
               {
                  it.gs.zero = 0.0;
                  for (SurfaceState &m : it.gsMore)
                     m.zero = 0.0;
               }
               // tightness counters: a pair kept in the ball's solve, and whether it touched in this cycle (per kind; diagnostics only)
               if (!it.partOut)
               {
                  const ContactMovingPart *const mp = it.r->MovingPart();
                  auto &k = g_qe.kept[mp != nullptr ? std::string("part ") + typeid(*it.r).name() : it.r->FinderObjIndex() >= 0 ? std::string("static mesh") : std::string("static other")];
                  k.first++;
                  if (it.cycleTouched || it.touching)
                     k.second++;
               }
               it.cycleTouched = false;
            }
            w.e->endStored = held;
         }
         const double e = (w.held ? 0.0 : BodyMotion(w, g)) + held + w.endPair;
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
      bool CycleNear(const V3 &c, const double R, const double reach, std::vector<int> &list, int &n) const override { return mp != nullptr && mp->GroupCycleNear(surface, c, R, reach, list, n); }
      int CycleMisses(const V3 &c, const double R, const double dq, const std::vector<int> &list) const override { return mp != nullptr ? mp->GroupCycleMisses(surface, c, R, dq, list) : 0; }
      int PathReach(const V3 &s, const V3 &d, const double R) const override { return mp != nullptr ? mp->GroupPathReach(surface, s, d, R) : r->PathReach(*it, s, d, R); }
      double PathMotion(const V3 &P, const double dLen) const override { return mp != nullptr ? mp->GroupPathMotion(surface, P, dLen, trial) : -1.0; }   // (group A)
      double PathSlope() const override { return mp != nullptr ? mp->GroupPathSlope(surface) : 1.0; }
      double PathCoverSag(int &kind) const override { return mp != nullptr ? mp->GroupPathCoverSag(surface, kind) : -1.0; }   // (the self-check, VPX_FINDER_CHECK)   // (the path test; a part's surface: its own answer, -1 for one that moves: build 3)
      void Probe(const V3 &c, double R, GroupStep::Hit &h) const override
      {
         const ContactProf::Stage profStage(ContactProf::StGroupShape);   // P3a profile: the calling stage of the checks here (diagnostics only)
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
      const ContactProf::Stage profStage(ContactProf::StRingHandIn);   // P3a profile: the calling stage of the checks here (diagnostics only)
      if (e.groupLoop)
         return;
      for (SurfaceContact &it : e.items)
      {
         it.gs.touching = it.touching && e.dwellS > 0.0;
         it.gs.xi = it.st.xi;
         it.gs.ra = it.st.ra;
         it.gs.zero = it.gs.reff = 0.0;
         // one entry per ball: an item not touching that the ball is found inside has stored nothing (rule 10b: the overlap came from no drive): its zero is its depth, by the
         // group's own probe of a fixed surface (SurfaceShape::Probe, ProbePatch); a moving part's surfaces take the same rule in SetShapeItem
         if (!it.gs.touching && it.r->MovingPart() == nullptr && !it.NoFaces(c))
         {
            const TubeProbe pr = it.r->ProbePatch(bs, it, c);
            if (pr.overlap && pr.delta > 0.0)
            {
               it.gs.zero = pr.delta;
               g_qe.firstSightZero++;
               if (e.ball != nullptr)
                  FoundInside(e.ball, it.r->Name(), pr.delta);
            }
            continue;
         }
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
   inline static std::unordered_map<const HitBall *, long long> g_coilBooked;   // a coil's current change: the physics step its work was last booked for each ball (once per ball per step); cleared at play end (EndTable)

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

   // The group parts in this island's solve. A part not at rest, or touching a ball of the island, joins. Before its first touch a part at rest moves only by its own drives, and every ball
   // that can touch it this cycle is in this island: it joins when the ball's travel (reachS, the energy reach) can close the gap ReachLinks probed less the part's own move (EnergyMove with
   // no feed; no reach set: it joins). Out of the solve it moves alone (AdvanceAlone), so nothing about it stops; the script and a nudge can only act between cycles.
   static void PartsInSolve(Cycle &cy, const ContactSystem::Islands::Island &isl, std::vector<const ContactMovingPart *> &in)
   {
      in.clear();
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
                  joins = w.reachS < 0.0 || it.reachGap - mp->EnergyMove(w.c, 0.0, cy.dt, nullptr) <= w.reachS;
                  it.solveOutLast = !joins;
               }
               if (joins && mp->AtRest())
                  s_solveRest++;
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

   // [Features] aircore_in_step: the air-core coils' pull on the group's balls, read where the ball is in the step (Magnets::AircoreForce, the law of AddVelocityIncrement's
   // air-core branch), the coils' current held for the cycle. Per unit of Ball::forceScale (the ball's mass over the steel mass the magnets use), in the units fExt's held
   // increment is turned into (VPU/VPT per physics step -> m/s2): the same force as the held one at the same place.
   class AircorePullForce final : public GroupStep::BallForce
   {
   public:
      V3 At(const V3 &c, const double R, GroupStep::M3 *const slope) const override
      {
         static const double kUnits = (double)MS2TOVPUVPT2(1.0f) * PHYS_FACTOR * VS / PHYSICS_STEPTIME_S;   // (fExt's chain: the increment's units back to SI, 1 to float rounding)
         const double pos[3] = { c.x / U, c.y / U, c.z / U };
         double F[3], S[3][3];
         Magnets::AircoreForce(pos, R / U, F, slope != nullptr ? S : nullptr);
         if (slope != nullptr)
            for (int i = 0; i < 3; i++)
               for (int j = 0; j < 3; j++)
                  slope->a[i][j] = kUnits * S[i][j];
         return V3(F[0], F[1], F[2]) * kUnits;
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
      static const AircorePullForce aircorePull;
      if (ContactSystem::s_aircoreInStep && Magnets::Active())
         g.AddForce(&aircorePull);
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
      static std::vector<Body *> bodyOf;   // re-engage: each group ball's body (its surface items)
      bodyOf.assign(isl.members.size(), nullptr);
      g.Size((int)isl.members.size());
      std::vector<ContactSurfaceHost *> cords;
      static std::vector<ContactMovingPart *> parts;   // the moving parts on the group in this island
      parts.clear();
      static std::vector<const ContactMovingPart *> partsIn;
      PartsInSolve(cy, isl, partsIn);
      s_solveParts += partsIn.size();   // (PartsInSolve counters)
      s_solveIslands++;
      bool gSet = false;
      static std::vector<int> order;   // the island's members, held balls last (KickerDesign 1h): a pair with a held ball has it as its b, the body forces do not move (GroupStep::Ball::held)
      order.clear();
      for (const int i : isl.members)
         if (!cy.bodies[i].held)
            order.push_back(i);
      for (const int i : isl.members)
         if (cy.bodies[i].held)
            order.push_back(i);
      for (const int i : order)
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
         bodyOf[ib] = &w;
         GroupStep::Ball &b = g.ball[ib];
         b.m = w.bs.m;
         b.R = w.bs.R;
         b.c = w.c;
         b.adds = ContactSurfaceHost::FinderAddsBind(w.pb, &b.adds);   // the inside check's rounding room: the additions to the ball's position since its query (0 unless stepped before in the cycle)
         b.v = w.v;
         b.w = w.w;   // (c, v, w: the ball's owned state, as loaded: one state, nothing to pick)
         b.fExt = w.dragged ? (w.extBaseVS * (1.0 / PHYSICS_STEPTIME_S)) * (-b.m) : V3();   // (a dragged ball: the hand holds it against gravity and the nudge, every ball's g.set.g)
         b.held = w.held;
         b.heldA = w.heldA;
         b.feedE = w.feedE;   // (the energy reach: -1 when not set)
         b.reachS = w.reachS;
         b.vTop = w.vTop;   // (energy skip test: -1 when not set)
         b.eTop = w.eTop;
         b.eGroup = w.eGroup;
         b.aOwn = w.vTop >= 0.0 ? w.gravAcc + w.magReach : -1.0;   // (energy skip test: its own bound's acceleration; -1 when ReachEnergy did not set it)
         if (w.held)
            b.w = V3();   // (spin zero)
         if (!gSet && w.extBaseVS.Length() > 0.0)
         {
            g.set.g = w.extBaseVS * (1.0 / PHYSICS_STEPTIME_S);   // gravity along the table's slope, and the nudge
            gSet = true;
         }
         Sources::EddyBrake(w);   // (the old loop's own: its rate and direction, held for the cycle as the coil's current is)
         b.brakeRate = w.eddyRate;
         b.brakeAxis = w.eddyDir;
         // a cells magnet's pull: a stored energy of the ball's place; on the body it is the same acceleration as AddVelocityIncrement's (m / the steel mass it uses)
         b.fieldScale = g.nf > 0 && w.pb->m_magnetic && !w.held && !w.dragged ? b.m / (7800.0 * 4.0 / 3.0 * M_PI * b.R * b.R * b.R * std::max(0.001f, w.pb->m_d.m_mass)) : 0.0;   // (a held ball: no fields)
         // [Features] aircore_in_step: an air-core coil's pull read at the step's midpoint (AircorePullForce), the same scale (a held ball: none; a ball whose magnets are off this cycle: none)
         b.forceScale = g.nbf > 0 && w.magnetsOn && !w.held ? b.m / (7800.0 * 4.0 / 3.0 * M_PI * b.R * b.R * b.R * std::max(0.001f, w.pb->m_d.m_mass)) : 0.0;
         if (w.magAcc > 0.0)   // what stays held for the cycle: the other balls' part, the magnets not on a field (nor, with aircore_in_step, the air-core coils): a drive
         {
            const float mp[3] = { (float)(w.c.x / U), (float)(w.c.y / U), (float)(w.c.z / U) };   // (the magnets' float interface, from the owned centre)
            // the other balls' part: a stored energy of place when the solve gave a form (pairform.h, 2026-11-03, added below), else held as before (the dipoles)
            const bool form = Magnets::HasPairForm(w.pb);
            float mdv[3] = { form ? 0.f : w.pb->m_magPair[0], form ? 0.f : w.pb->m_magPair[1], form ? 0.f : w.pb->m_magPair[2] };
            Magnets::AddVelocityIncrement(mp, w.pb->m_d.m_radius, w.pb->m_d.m_mass, mdv, true, ContactSystem::s_aircoreInStep);
            b.fExt = V3(mdv[0] * VS, mdv[1] * VS, mdv[2] * VS) * (b.m / PHYSICS_STEPTIME_S);
         }
         if (w.band != nullptr && std::find(cords.begin(), cords.end(), w.band) == cords.end())
            cords.push_back(w.band);
         if (w.e != nullptr)
            for (ContactSurfaceHost *const b : w.e->bandsMore)   // (each other band that holds it steps with this island)
               if (std::find(cords.begin(), cords.end(), b) == cords.end())
                  cords.push_back(b);
         if (ring)
         {
            // (a ball the load did not see moved outside the engine is where the engine left it: a depth found there is the faces in reach changing)
            cause[ib] = std::find(g_placed.begin(), g_placed.end(), w.pb) != g_placed.end() ? DepthCause::Placed
               : w.e != nullptr && !w.e->groupLoop && w.e->endStored >= 0.0 ? DepthCause::HandedOver
               : DepthCause::Faces;
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
               if (it.partOut)   // (island reach: a part out of this ball's reach is solved elsewhere)
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
                        mp->GroupPart()->eOwnJ = mp->eOwnJ;   // (energy skip test: this cycle's only; -1 when ReachEnergy did not set it this cycle)
                        mp->eOwnJ = -1.0;
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
            V3 c = OwnC(out);   // (a partner outside the group: its owned centre)
            for (const Body &w : cy.bodies)
               if (w.pb == out) { c = w.c; break; }
            t.o[0] = c.x; t.o[1] = c.y; t.o[2] = c.z;
         }
         g.AddForm(t);
      }
      // which shapes are a ball's surface (SurfaceShape: a static surface or a part's) and which a cord's posts (rubbers on posts: GroupStepDesign 13). The surfaces come first, then each cord's
      // posts; a surface joined mid-step (re-engage) comes after them, so the kind is told by this flag, never by the shape's index
      static std::vector<char> surfShape;
      surfShape.assign((size_t)g.ns, 1);
      const auto isSurf = [&](const int s) { return s >= 0 && s < (int)surfShape.size() && surfShape[(size_t)s] != 0; };
      for (ContactSurfaceHost *const c : cords)
      {
         c->Cord()->fixed->found.clear();
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
      surfShape.resize((size_t)g.ns, 0);   // (the posts)
      g.Configure();
      // a coil's current changed since the last cycle: what the balls' places store changed with it, with nothing moved. That is the coil's work, a drive. Once per ball per physics step (a split cycle comes back here)
      for (const PullField &pf : pulls)
         if (pf.p.amps != pf.p.ampsBefore)
         {
            std::unordered_map<const HitBall *, long long> &booked = g_coilBooked;
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
      // a band found inside (the cord, its posts, its plate): no drive stored that overlap, so each piece took the depth it was found at as its zero (Label, SetShapeItem), released
      // as the ball backs out; rule 8: one WARN per band per table start, every one counted
      for (ContactSurfaceHost *const c : cords)
      {
         c->Cord()->fixed->RestorePostItems(g);
         for (const std::pair<int, double> &f : c->Cord()->fixed->found)
         {
            g_qe.firstSightZero++;
            FoundInside(balls[f.first], "band " + c->Name(), f.second);
         }
         c->Cord()->fixed->found.clear();
      }
      // what a touch keeps between cycles: the shear spring's stretch
      for (int k = 0; k < g.ni; k++)
      {
         GroupStep::Item &it = g.item[k];
         if (it.b >= 0 || !isSurf(it.s))
            continue;
         const SurfaceShape *const sh = static_cast<const SurfaceShape *>((const GroupStep::Shape *)g.shape[it.s]);
         const SurfaceState &kept = sh->it->Surf(sh->surface);
         g.SetShapeItem(k, kept.touching, kept.zero, kept.xi, kept.ra, &kept.tZero);
         if (sh->mp != nullptr && !kept.touching && kept.zero <= 0.0 && it.dZero > 0.0)   // rule 8 for a moving part: its item took a found-inside zero here
         {
            g_qe.firstSightZero++;
            FoundInside(balls[it.a], sh->r->Name(), it.dZero);
         }
      }
      // ... then each kept touch's depth now against its end depth, with every contact set (a massless part's angle is solved from all of them, as it was at the step's end)
      for (int k = 0; k < g.ni; k++)
      {
         GroupStep::Item &it = g.item[k];
         if (it.b >= 0 || !isSurf(it.s))
            continue;
         const SurfaceShape *const sh = static_cast<const SurfaceShape *>((const GroupStep::Shape *)g.shape[it.s]);
         const SurfaceState &kept = sh->it->Surf(sh->surface);
         if (ring && it.touching && g_contactKeptEnd(balls[it.a]))
         {
            GroupStep::Hit hp;
            g.ProbeItem(k, hp);
            const double was = sh->mp != nullptr ? kept.endDepth : sh->it->endDepth;
            if (hp.nt > 0 && sh->mp != nullptr)   // a union of terms: each term by the same rule, its depth now (raw less its own zero) against its own end depth; the excess goes
            {                                     // into that term's zero, so force (the terms) and book (the union of the zeroed terms) read one U = sum_t s_t V(a_t - z_t)
               for (int t = 0; t < hp.nt; t++)
               {
                  double end = 0.0;   // a term with no end depth kept was not pressing then (V = 0 there)
                  for (const std::pair<int, double> &e : kept.tEnd)
                     if (e.first == hp.t[t].id)
                     {
                        end = std::max(e.second, 0.0);
                        break;
                     }
                  if (hp.t[t].d > end + 1e-9)
                  {
                     size_t j = 0;
                     while (j < it.tZero.size() && it.tZero[j].first != hp.t[t].id)
                        j++;
                     if (j == it.tZero.size())
                        it.tZero.emplace_back(hp.t[t].id, 0.0);
                     it.tZero[j].second += hp.t[t].d - end;
                     it.dZero = std::max(it.dZero, it.tZero[j].second);
                  }
               }
               if (g.set.jointZeroFloor)   // [Contact] joint_zero_floor: a raised piece raises its joints' zeros with it (d_J <= min d_p, U >= 0)
                  it.FloorJointZeros(hp);
            }
            else if (hp.delta > was + 1e-9)
               it.dZero += hp.delta - was;
         }
         // what the zero took this cycle, by who moved the ball
         // (one entry per ball: first sight per item: the zero rule took its start depth, SetShapeItem / EnergyCheckStart, as a placed ball: nothing to report)
         const DepthCause why = cause[it.a] != DepthCause::Placed && sh->it->endHeld < 0.0 ? DepthCause::FirstSight : cause[it.a];
         if (ring && it.dZero > kept.zero + 1e-9 && why != DepthCause::Placed && why != DepthCause::HandedOver && why != DepthCause::FirstSight)
            DepthFault(why, balls[it.a], sh->r->Name(), it.dZero - kept.zero);
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
               hitOn[k] = isSurf(it.s) && (static_cast<const SurfaceShape *>(g.shape[it.s])->mp != nullptr ? it.touching : static_cast<const SurfaceShape *>(g.shape[it.s])->it->touching);
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
            if (!ring && (gi.b >= 0 || !isSurf(gi.s) || static_cast<const SurfaceShape *>(g.shape[gi.s])->mp == nullptr))   // (a cord island: a part's touches only)
               continue;
            if (gi.touching && gi.b < 0 && isSurf(gi.s))
               const_cast<SurfaceContact *>(static_cast<const SurfaceShape *>(g.shape[gi.s])->it)->cycleTouched = true;   // (one entry per ball counters)
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
            else if (isSurf(gi.s))
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
                     int ballId = -1;   // (the ball's own id, Ball::m_id)
                     if (balls[gi.a]->m_pBall != nullptr)
                        balls[gi.a]->m_pBall->get_ID(&ballId);
                     sh->mp->TouchBegan(sh->surface, gi.a, ballId, -(vStep[gi.a] - g.SurfaceSpeed(*sh, hp)).Dot(hp.n), hp.n);
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
         const ContactProf::Stage profStage(ContactProf::StRingCap);   // P3a profile: the calling stage of the checks here (diagnostics only)
         if (!ring)
            return h;
         const double a = g.set.g.Length();
         for (int s = 0; s < g.ns; s++)
         {
            if (!isSurf(s))
               continue;
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
      const double cordE0 = cordEnergyOn ? g.Energy() : 0.0;
      // a part's books (always on, rule 8): the group's energy across the cycle against what its drives did and its losses took
      const double partE0 = hasPart ? [&]() {
         double e = g.Energy();
         for (int p = 0; p < g.np; p++)
            e -= g.part[p]->PendingNetJ();   // an impact taken at CycleBegin (a driven stop's touch) is in Energy() and in this cycle's books: start before it (as GroupStep::CycleBook)
         return e;
      }() : 0.0;
      const double hFix = GroupStep::GroupFixedStepS();
      // the group's moving parts and cords by name
      const auto groupParts = [&]() {
         std::vector<std::string> names;
         for (const ContactMovingPart *const mp : parts)
            for (int s = 0; s < g.ns; s++)
               if (isSurf(s) && static_cast<const SurfaceShape *>(g.shape[s])->mp == mp)
               {
                  names.push_back(static_cast<const SurfaceShape *>(g.shape[s])->r->Name());
                  break;
               }
         for (ContactSurfaceHost *const c : cords)
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
                  const SurfaceShape *const sh = isSurf(g.item[k].s) ? static_cast<const SurfaceShape *>(g.shape[g.item[k].s]) : nullptr;
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
                  const SurfaceShape *const sh = isSurf(it.s) ? static_cast<const SurfaceShape *>(g.shape[it.s]) : nullptr;
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
      // re-engage (runs/2026-11-13-finder-reengage/design.md). Law: non-penetration, a surface acts on a ball from its first touch. An object the ball's inside check
      // found beyond its predicted reach, inside an attempt, joins the ball's surfaces (as the cycle's gather makes them: Regather) and the group (its shape after the others), with nothing kept
      // and its zero from the step's start (SetShapeItem, as the cycle's start sets a new touch), when this step's path can reach its box (FinderTakeReachable; the rest stays pending for the
      // cycle's later steps). It joins as a skipped contact: the check after the solve (GroupStep::Group::Solve, FarStillFar) runs the attempt again only if it is reached.
      // Not joined here (a fault, counted and logged: FinderJoinNotHandled): a ball with no surface items of its own. In an island with a cord the joined shape comes after the posts (isSurf tells them apart).
      static std::vector<ContactHost *> joinHosts;
      static std::vector<int> joinedK;   // the items joined this cycle; joinedTouch: touching at the end of some step of the cycle (counted: FinderJoinTouched)
      static std::vector<char> joinedTouch;
      joinedK.clear();
      joinedTouch.clear();
      double joinSetupNs = 0.0;   // VPX_CHECK_PROFILE: this cycle's join setup time (the lambda below), ns
      g.timeJoins = ContactProf::CheckOn();
      g.joinRerunNs = 0.0;
      g.joinReruns = 0;
      g.midJoin = [&](const bool allow) -> int {
            // the ball's inside check: every ball of the group is checked against its own cycle query at the attempt's end, whatever its items (design.md section 4: the
            // inside check first). The attempt moves the ball's centre on the straight segment c -> c + MoveOf, whose farthest point from the query's centre is one of its ends;
            // the start passed the last check, so the end covers the move. A ball with no BallContact entry has no cycle query (FinderQuery runs per entry): nothing to check.
            for (int ib = 0; ib < g.nb; ib++)
               if (bodyOf[ib] != nullptr && bodyOf[ib]->e != nullptr)
                  ContactSurfaceHost::FinderInsideBall(balls[ib], g.ball[ib].c, g.MoveOf(ib));
            if (!ContactSurfaceHost::FinderPendAny())
               return 0;
            if (!allow)
               return -1;
            const bool timed = ContactProf::CheckOn();
            const auto t0 = timed ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point();
            int added = 0;
            bool inside = false;
            const int nbNow = g.nb;
            for (int ib = 0; ib < nbNow; ib++)
            {
               ContactSurfaceHost::FinderTakeReachable(balls[ib], g.ball[ib].c, g.ball[ib].c + g.MoveOf(ib), g.ball[ib].R, joinHosts);   // (build 3: only what this step's path can reach; the rest stays pending)
               if (joinHosts.empty())
                  continue;
               Body &w = *bodyOf[ib];
               HitBall *const pb = balls[ib];
               for (ContactHost *const o : joinHosts)
               {
                  if (w.e == nullptr)
                  {
                     ContactSurfaceHost::FinderJoinNotHandled(pb, o->m_name, "the ball has no surface items");
                     continue;
                  }
                  std::vector<SurfaceContact> &items = w.e->items;
                  if (!o->m_set.contact || o->m_set.band || std::any_of(items.begin(), items.end(), [o](const SurfaceContact &it) { return it.r == o; }))
                     continue;   // (as Regather leaves them out)
                  SurfaceContact nit;
                  nit.r = o;
                  nit.rp = o->m_set.ToParams();
                  nit.hLayer = o->LayerM();
                  nit.st.Reset();
                  const uintptr_t old0 = (uintptr_t)items.data(), old1 = (uintptr_t)(items.data() + items.size());
                  items.push_back(std::move(nit));
                  if ((uintptr_t)items.data() != old0)   // the shapes point at the ball's items: moved with them
                     for (SurfaceShape &s : shapes)
                        if ((uintptr_t)s.it >= old0 && (uintptr_t)s.it < old1)
                           s.it = items.data() + ((uintptr_t)s.it - old0) / sizeof(SurfaceContact);
                  SurfaceContact &it = items.back();
                  ContactSurfaceHost::FinderBindOne(pb, it);   // (a static mesh: its faces come from the ball's query; no Refresh)
                  shapes.emplace_back();
                  SurfaceShape &sh = shapes.back();
                  sh.r = it.r;
                  sh.it = &it;
                  sh.bs = w.bs;
                  sh.rp = it.rp;
                  sh.layerM = it.hLayer;
                  sh.onlyBall = ib;
                  const int k0 = g.AddShapeMidStep(&sh);
                  surfShape.resize((size_t)g.ns, 1);   // (the joined surface: after any cord's posts)
                  hitOn.resize(g.ni, 0);
                  for (int k = k0; k < g.ni; k++)
                  {
                     g.SetShapeItem(k, false, 0.0, V3());
                     hitOn[k] = 0;
                     joinedK.push_back(k);
                     joinedTouch.push_back(g.item[k].dZero > 0.0 ? 1 : 0);
                     if (g.item[k].dZero > 0.0)
                     {
                        inside = true;
                        PLOGW << "Contact finder: ball " << (const void *)pb << " joined " << o->m_name << " mid-step already " << g.item[k].dZero * 1e6 << " um inside at the step's start (its zero takes it, nothing stored), t=" << g_pplayer->m_time_msec;
                     }
                  }
                  if (ContactTraceOn())
                     PLOGI << "ContactJoin " << o->m_name << " t=" << g_pplayer->m_time_msec << " mid-step (found by the inside check beyond the ball's predicted reach), ball " << (const void *)pb << (cords.empty() ? "" : " (island with a cord)");
                  if (timed)   // the per-join bookkeeping (a label and a map entry): profile only
                     ContactSurfaceHost::FinderJoinNote((cords.empty() ? std::string() : std::string("cord island:")) + (o->m_surfaceKind.empty() ? std::string("?") : o->m_surfaceKind), false, false, 0.0);
                  added++;
               }
            }
            if (added > 0)
            {
               const double ns = timed ? std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count() : 0.0;
               joinSetupNs += ns;
               ContactSurfaceHost::FinderJoinNote(std::string(), true, inside, ns);
            }
            return added;
         };
      const auto tSteps0 = g.timeJoins ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point();
      const bool ok = g.StepsOver(cy.dt, done, ruleCap, [&]() {
         ringStep();
         book.AfterStep();
         partEnergy();
         for (size_t j = 0; j < joinedK.size(); j++)   // re-engage: a joined contact touching at a step's end
            joinedTouch[j] = joinedTouch[j] || g.item[joinedK[j]].touching ? 1 : 0;
      });
      g.midJoin = nullptr;
      if (g.timeJoins && (joinSetupNs > 0.0 || g.joinReruns > 0))   // VPX_CHECK_PROFILE: a cycle with a join, its setup and re-run attempts against the cycle's whole step time
         ContactSurfaceHost::FinderJoinTiming(joinSetupNs, g.joinRerunNs, g.joinReruns, std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - tSteps0).count());
      if (!joinedK.empty())
         ContactSurfaceHost::FinderJoinTouched((int)joinedK.size(), (int)std::count(joinedTouch.begin(), joinedTouch.end(), (char)1));
      if (g.joinCapHit > 0)
      {
         ContactSurfaceHost::FinderJoinCapHit(g.joinCapHit);
         g.joinCapHit = 0;
      }
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
      // each cord's energy above its no-ball state and its plate's motion, for the next cycle's feed (CordBand::heldJ)
         for (ContactSurfaceHost *const c : cords)
         {
            CordBand &cb = *c->Cord();
            const int s0 = cb.fixed->postFirst, s1 = s0 + (int)cb.fixed->postShape.size() + (cb.fixed->pl.on ? 1 : 0);   // (its posts, then its plate: AddPostShapes)
            cb.heldJ = std::max(cb.fixed->Energy(&g.ball[0], g.nb) - cb.fixed->Energy(&g.ball[0], 0) + cb.fixed->KineticJ() + g.ShapesContactJ(s0, s1), 0.0);   // (and each ball's squash on its posts and plate)
         }
      // what each moving part's contacts hold (its ball contacts' and compliant stops' stored energy), for the next cycle's feed
         for (size_t p = 0; p < parts.size(); p++)   // (parts[p] is the group's part p: the cords are added after them)
            parts[p]->squashJ = g.PartContactJ((int)p);
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
                          names += isSurf(g.item[k].s) ? static_cast<const SurfaceShape *>(g.shape[g.item[k].s])->r->Name() + "|" : std::string("post|");
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
      for (const int i : order)   // (the group's order: held balls last)
      {
         if (ib >= g.nb)
            break;
         Body &w = cy.bodies[i];
         if (ring && w.e != nullptr)
            w.e->seatPending = false;   // (a surface it was found inside counted from zero)
         w.c = g.ball[ib].c;
         ContactSurfaceHost::FinderAddsBind(w.pb, nullptr);   // (the count kept for the cycle's rest; the group's slot is the next island's)
         w.v = g.ball[ib].v;
         w.w = g.ball[ib].w;   // (a step that could not be solved: what is left of the cycle is not stepped, so nothing is gained: below)
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
         if (gi.b >= 0 || !isSurf(gi.s))
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
            kept.tEnd.clear();
            for (int t = 0; t < hp.nt; t++)
               kept.tEnd.emplace_back(hp.t[t].id, hp.t[t].d);
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
         // one entry per ball: an item not touching keeps its zero while the ball is still inside it (the zero follows the ball out, groupstep BackOut); zero 0 only once the ball is
         // clear: its raw depth at the cycle's end <= 0. Raw depth: ProbeItem is the shape's own probe less the item's zero, so adding gi.dZero back gives the geometric depth
         double zeroKept = gi.dZero;
         if (!gi.touching && gi.dZero > 0.0)   // (a zero of 0 stays 0)
         {
            GroupStep::Hit hr;
            g.ProbeItem(k, hr);
            bool termsInside = false;
            for (int t = 0; t < hr.nt; t++)
               if (gi.TermZero(hr.t[t].id) > 0.0 && hr.t[t].d + gi.TermZero(hr.t[t].id) > 0.0)
                  termsInside = true;
            if (gi.tZero.empty() ? hr.delta + gi.dZero <= 0.0 : !termsInside)   // a union of terms: clear when every zeroed term's raw depth is <= 0
               zeroKept = 0.0;
         }
         kept.zero = zeroKept;
         kept.tZero.clear();
         if (zeroKept > 0.0)
            kept.tZero = gi.tZero;
      }
      for (ContactSurfaceHost *const c : cords)
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
                  snprintf(s, sizeof s, " item a%d b%d s%d%s zero %.3g", g.item[k].a, g.item[k].b, g.item[k].s, g.item[k].b < 0 && !isSurf(g.item[k].s) ? " post" : "", g.item[k].dZero);
                  how += s;
               }
            for (int i = 0; i < g.nb; i++)
            {
               snprintf(s, sizeof s, " ball%d v (%.4f %.4f %.4f) w (%.2f %.2f %.2f)", i, g.ball[i].v.x, g.ball[i].v.y, g.ball[i].v.z, g.ball[i].w.x, g.ball[i].w.y, g.ball[i].w.z);
               how += s;
            }
            for (ContactSurfaceHost *const c : cords)
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
         for (ContactSurfaceHost *const c : cords)
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
         for (ContactSurfaceHost *const c : cords)
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
   static void CordAlone(ContactSurfaceHost *const band, const double dt)
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
      cy.releases.clear();
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

   static void SleepWakeGroup(const size_t gi, const char *cause)
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
      }
      g_sleep.erase(g_sleep.begin() + (long)gi);
   }
   static void SleepWakeBall(const HitBall *pb, const char *cause)
   {
      if (pb == nullptr || pb->m_sleepGroup == 0)
         return;
      for (size_t gi = 0; gi < g_sleep.size(); gi++)
         if (g_sleep[gi].id == pb->m_sleepGroup)
         {
            SleepWakeGroup(gi, cause);
            return;
         }
   }
   // Something no step did: a ball moved, given a speed or spin, captured, or put somewhere (a script, the driver, a kicker, a stock collision)
   static const char *SleepMoved(const SleepGroup &g)
   {
      for (const HitBall *const pb : g.balls)
      {
         if (KickerHold::Held(pb))
            return "ball held by a kicker";
         if (!SameV(pb->m_d.m_pos, pb->m_sleepPos))
            return "ball moved outside the step";
         if (!SameV(pb->m_d.m_vel, Vertex3Ds(0.f, 0.f, 0.f)) || !SameV(pb->m_angularmomentum, Vertex3Ds(0.f, 0.f, 0.f)))
            return "ball given a speed or spin outside the step";
         if (std::find(g_placed.begin(), g_placed.end(), pb) != g_placed.end())
            return "ball placed";
      }
      return nullptr;
   }
   // A part the balls touch that a drive woke or that moved; a part not at rest that can reach a ball of the group within the cycle by the energy reach (the rule of PartsInSolve):
   // its gap to the ball less its own move (EnergyMove with no feed: its kinetic energy and drives) within the ball's own travel (LoneReach) (<=)
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
      std::vector<double> travel(g.balls.size(), -1.0);   // each ball's travel over the cycle (LoneReach), found when a part needs it
      for (ContactHost *const o : ContactHost::Registry())
      {
         ContactMovingPart *const mp = o->MovingPart();
         if (mp == nullptr || PartQuiet(mp))
            continue;
         for (size_t b = 0; b < g.balls.size(); b++)
         {
            HitBall *const pb = g.balls[b];
            if (travel[b] < 0.0)
            {
               Cycle cy;
               LoneReach(cy, EnsureOne(pb), pb, dt);
               travel[b] = std::max(cy.bodies[0].reachS, 0.0);
            }
            Body s;
            LoadBody(s, pb);
            LoadState(s);
            mp->SetBall(s.bs);
            if (-mp->Probe(s.c).delta - mp->EnergyMove(s.c, 0.0, dt, nullptr) <= travel[b])
               return "a moving part can reach it";
         }
      }
      // a cord that can move this cycle (its plate can, or balls hold it) wakes a sleeping ball it can reach, by the take's own test (BandQueryGap <= the ball's travel); the
      // held balls' travels are their last cycle's (BallContact::reachS)
      for (ContactHost *const o : ContactHost::Registry())
      {
         ContactSurfaceHost *const band = o->m_set.band && o->m_enabled ? o->BandHost() : nullptr;
         if (band == nullptr || !band->CordModel() || !band->CordOn())
            continue;
         std::vector<double> dNow, push;
         bool held = false;
         for (auto *const ball : g_pplayer->m_vball)
            if (ball->m_hitBall.m_bandOwner == band || (ball->m_hitBall.m_contact != nullptr && std::find(ball->m_hitBall.m_contact->bandsMore.begin(), ball->m_hitBall.m_contact->bandsMore.end(), band) != ball->m_hitBall.m_contact->bandsMore.end()))
            {
               const HitBall &h = ball->m_hitBall;
               held = true;
               band->Cord()->fixed->SpanPush(OwnC(&h), (double)h.m_d.m_radius * U, h.m_contact != nullptr ? h.m_contact->reachS : 0.0, dNow, push);
            }
         if (!held && band->Cord()->fixed->PlateMoveM(dt) <= 0.0)
            continue;   // a cord at rest with nothing to move it: the sleeping balls are where they were when they slept
         for (size_t b = 0; b < g.balls.size(); b++)
         {
            HitBall *const pb = g.balls[b];
            if (travel[b] < 0.0)
            {
               Cycle cy;
               LoneReach(cy, EnsureOne(pb), pb, dt);
               travel[b] = std::max(cy.bodies[0].reachS, 0.0);
            }
            Body s;
            LoadBody(s, pb);
            LoadState(s);
            if (BandQueryGap(band, s.c, s.bs.R, dNow, push, dt) <= travel[b])
               return "a band can reach it";
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
            SleepWakeGroup(gi, cause);
         else
            gi++;
      }
   }
   // At each contact step, before the gather: moves outside the step since the last one, an awake ball within reach (the islands' own test); the pairs of sleeping
   // balls are kept
   static void SleepInStep(Cycle &cy)
   {
      SleepTable();
      if (!s_append)
         g_sleepCycleDt += cy.dt;
      for (size_t gi = 0; gi < g_sleep.size();)
         if (const char *const cause = SleepMoved(g_sleep[gi]))
            SleepWakeGroup(gi, cause);
         else
            gi++;
      for (bool woke = !g_sleep.empty(); woke;)
      {
         woke = false;
         // ball pair reach: the islands' reach (ReachLinks, EnergyReach) with the groups not built yet, all awake balls as one group: the energy of every awake ball can reach an
         // awake ball, J (spin included), and their forces, m/s2 (gravity, nudge and a magnet's pull, as Sources::ReachAcc); a sleeping ball's travel is 0
         double eAwake = 0.0, aAwake = 0.0;
         for (auto *const a : g_pplayer->m_vball)
         {
            const HitBall &h = a->m_hitBall;
            if (h.m_sleepGroup == 0 && !KickerHold::Held(&h))
            {
               const V3 v = OwnV(&h), am = OwnW(&h) * ((double)h.Inertia() * DEFAULT_STEPTIME_S);   // (the owned state; am: the spin as angular momentum, VPX units)
               eAwake += 0.5 * (h.m_d.m_mass * ContactSystem::kMassUnitKg * v.Dot(v) + am.Dot(am) / h.Inertia() * ContactSystem::kMassUnitKg * VS * VS);
               const Vertex3Ds g = h.ExternalVelocityIncrement(false);
               aAwake += V3(g.x * VS, g.y * VS, g.z * VS).Length() / PHYSICS_STEPTIME_S;
               if (h.m_magnetic && Magnets::Active())
               {
                  const Vertex3Ds am = h.ExternalVelocityIncrement(true);
                  aAwake += V3((am.x - g.x) * VS, (am.y - g.y) * VS, (am.z - g.z) * VS).Length() / PHYSICS_STEPTIME_S;
               }
            }
         }
         for (auto *const a : g_pplayer->m_vball)
         {
            const HitBall *const pa = &a->m_hitBall;
            if (pa->m_sleepGroup != 0 || KickerHold::Held(pa))
               continue;
            const double ra = (double)pa->m_d.m_radius * U;
            const double reach = EnergyReach(0.0, eAwake, pa->m_d.m_mass * ContactSystem::kMassUnitKg, aAwake, cy.dt);
            for (auto *const b : g_pplayer->m_vball)
            {
               const HitBall *const pb = &b->m_hitBall;
               if (pb->m_sleepGroup == 0)
                  continue;
               const double gap = (OwnC(pa) - OwnC(pb)).Length() - ra - (double)pb->m_d.m_radius * U;
               if (gap <= reach)   // the islands' energy pairing (a gap equal to the reach can meet)
               {
                  SleepWakeBall(pb, "an awake ball within reach");
                  woke = true;
               }
            }
         }
      }
      for (BallPairState &p : g_ballPairs)
         if (p.a != nullptr && p.b != nullptr && p.a->m_sleepGroup != 0 && p.b->m_sleepGroup != 0)
            p.used = true;
   }
   // After an island's step: its balls' holds at the cycle's end, the slow-creep notice, and the sleep itself
   static void SleepConsider(Cycle &cy, const ContactSystem::Islands::Island &isl)
   {
      const ContactSettings &cs = *cy.cs;
      if (!cs.sleep || s_append || g_sleepCycleDt < PHYSICS_STEPTIME_S * (1.0 - 1e-6))   // (the pieces of a replayed cycle: a ball still carries the rest of the cycle's gravity)
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
         if (b.band != nullptr || pb->m_bandOwner != nullptr || KickerHold::Held(pb))
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
               const ContactProf::Stage profStage(ContactProf::StLoadState);   // P3a profile: the calling stage of the checks here (diagnostics only)
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
         WriteBall(pb, OwnC(pb), V3(), V3());   // (the engine's write: the owned state stops with its floats)
         pb->m_sleepPos = pb->m_d.m_pos;
         pb->m_sleepGroup = g.id;
         pb->m_sleepCreepS = 0.0;
      }
      g.partIdle.assign(g.parts.size(), 0);
      if (SleepLogOk())
         PLOGI << "Pile sleep: group " << g.id << " (" << g.balls.size() << " balls, " << g.parts.size() << " moving parts) sleeps at t=" << g_pplayer->m_time_msec << " ms, max speed " << vmax * 1e3
               << " mm/s, max spin surface speed " << spinMax * 1e3 << " mm/s, deepest contact " << deep * 1e6 << " um";
      g_sleep.push_back(std::move(g));
   }

   static bool s_append;   // StepBalls: do not flush or clear the queued events: the steps of one replayed cycle chain
   static bool PreventOn()   // point 1 of the overlap standard; VPX_PREVENT=0 switches it off (test only)
   {
      static const bool on = []() { const char *const e = getenv("VPX_PREVENT"); return !(e != nullptr && atoi(e) == 0); }();
      return on;
   }
   static void CheckOverlaps(float hittime, double hS);   // hittime: engine time units; hS: the same in seconds
   static void SeatAtTouch(const std::vector<SurfaceContact> &items, const BallSI &bs, V3 &c, const V3 &v, double &delta0, double &rewound);
   static void SeatKeepsOthers(const std::vector<SurfaceContact> &items, const ContactSurface *self, const BallSI &bs, const V3 &c0, V3 &c);
   // The surfaces a ball at its position may touch, as an engagement takes them: host, the surface whose collider found the ball (taken without the checks below; null = none, the
   // start-up pass), how: as EngageImpl. c: the ball's centre (m), for the finder's query.
   static void GatherItems(ContactHost *host, const HitBall *pb, const BallSI &bs, const V3 &c, const char *how, std::vector<SurfaceContact> &items);
   // A ball found inside with no earlier clean state: minimum translation along the normal of each overlapped surface to touching distance, up to 8 passes (UnifiedContactDesign 11, point 6).
   // Velocities unchanged, nothing stored. delta0: the deepest overlap at the first pass (m, kept if larger); passes: passes taken; left: the deepest overlap left (m); inside (may be null):
   // the surfaces it moved the ball out of. Returns the distance moved (m).
   static double SeatNoPast(std::vector<SurfaceContact> &items, const BallSI &bs, V3 &c, double &delta0, int &passes, double &left, std::vector<const ContactSurface *> *inside);
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
   if (g_engaged.empty() && g_pplayer->m_vball.empty())   // (every ball is stepped here, held or free)
      return;
   if (dt <= 0.0)
      return;
   static Cycle cy;   // buffers are reused between cycles
   cy.dt = dt;
   cy.cs = &ContactSettings::Get();
   cy.bodies.clear();
   if (!s_append)   // (the steps of a replayed cycle add to the queue)
   {
      cy.releases.clear();
      cy.hitEvents.clear();
      cy.partHits.clear();
   }
   cy.hitFloor = std::min(0.5f, g_pplayer->m_physics->HitSpeedMin());
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
   EnsureEntries();   // one entry per ball: every ball in play has its entry before the gather
   g_qe.cycles++;
   static uint32_t nextLogMs = 10000;
   if (g_pplayer->m_time_msec < nextLogMs - 10000)   // (a new table)
      nextLogMs = 10000;
   if (g_pplayer->m_time_msec >= nextLogMs)   // the counters' running totals every 10 s of table time (the bench stops before the table's end)
   {
      nextLogMs = g_pplayer->m_time_msec - g_pplayer->m_time_msec % 10000 + 10000;
      QueryEngageLog("running");
      PLOGI << "Bands: " << s_bandBallCycles << " ball-cycles held by a band, " << s_bandSecond << " second-band takes";
      if (PairCheckOn())
         PLOGI << "Pair check (VPX_PAIR_CHECK): " << s_pairCheckTests << " every-pair tests, " << s_pairCheckMissed << " pairs the sorted list missed";
   }
   GatherBalls(cy);
   FindIslands(cy);
   CullCheck(cy);   // (the part cull, check mode only: the parts it left out against the full test)
   GatherStatic(cy);   // (the energy reach: each ball's static mesh query at its travel from the grouping)
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
// A free ball that can touch a surface this cycle is taken into the contact before the cycle's Step: the surface is engaged when the ball is within its release gap, or
// closing on it faster than it can cross the gap in two cycles (reach 2 vn dt + the release gap; a ball released at the gap, receding, is not taken back). The first touch is
// then found by the substeps, and the part's hit event is queued from there (PartHit).
// A band takes its balls by their energy query (BandQuery). dt: the cycle in seconds.
void ContactSystemImpl::EngageNear(const double dt)
{
   CPROF(EngageNear);
   const ContactProf::Stage profStage(ContactProf::StEngageNear);   // P3a profile: the calling stage of the checks here (diagnostics only)
   ContactProf::Cycle();
   {
      // one entry per ball: a ball that entered play since the last step (made or moved in by the driver, the script or a kicker) has its entry and its first query here, before
      // the Step of its first cycle: one contact path for every surface that touches a ball, from its first instant. Its reach is the cycle's own energy reach
      // (ReachLinks / ReachEnergy on a cycle of this one ball: its own kinetic energy and spin, gravity + nudge + magnets (ReachAcc), the drive work of every moving part in reach
      // (Regather part 1: each one a candidate); a new ball holds no stored energy and is in no ball pair yet), then the cycle's own query and gather (GatherStatic).
      // Only this ball's entry is written; the cycle's Step groups it with the others and regathers it as every ball.
      for (auto *const ball : g_pplayer->m_vball)
      {
         HitBall *const pb = &ball->m_hitBall;
         if (pb->m_contact != nullptr || !PairEligible(pb) || KickerHold::Held(pb))
            continue;
         BallContact &e = EnsureOne(pb);
         Regather(e, pb->m_d.m_pos, pb->m_d.m_radius, 1);
         Cycle cy;
         LoneReach(cy, e, pb, dt);
         GatherStatic(cy);
      }
   }
}

// ---- Engagement: which balls the contact system holds ----------------------------------------------------------------------------------------------------------------
// Go back to the instant the ball first touched: along its own velocity, by depth / closing speed, at most one physics step; if it is not closing, to touching distance along the normal.
void ContactSystemImpl::SeatAtTouch(const std::vector<SurfaceContact> &items, const BallSI &bs, V3 &c, const V3 &v, double &delta0, double &rewound)
{
   const ContactProf::Stage profStage(ContactProf::StSeatAtTouch);   // P3a profile: the calling stage of the checks here (diagnostics only)
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

void ContactSystemImpl::GatherItems(ContactHost *const host, const HitBall *const pb, const BallSI &bs, const V3 &c, const char *const how, std::vector<SurfaceContact> &items)
{
   using namespace RubberContact;
   const Vertex3Ds posStart = pb->m_d.m_pos;
   const float reach = pb->m_d.m_radius + 20.0f;   // VPX units: the ball radius plus what it can travel during a hit
   // shared finder: the static mesh hosts other than the one that found the ball come from a query at the engagement, radius R + 20 U (the reach above, to a face
   // instead of the host's box); kept in the registry's order
   thread_local std::vector<ContactHost *> hosts;
   hosts.clear();
   // the path test: the query's size is the one reach rule (ContactReachM; the acceleration unknown at engagement: 0, the inside check re-queries past it); 0 = 20 U as before
   ContactSurfaceHost::FinderQueryAt(pb, c, bs.R, ContactReachM(V3((double)pb->m_d.m_vel.x * VS, (double)pb->m_d.m_vel.y * VS, (double)pb->m_d.m_vel.z * VS).Length(), 0.0, PHYSICS_STEPTIME_S));
   ContactSurfaceHost::FinderHosts(pb, hosts);
   for (ContactHost *const o : ContactHost::Registry())
      if (o->FinderObjIndex() < 0 || o == host)
         hosts.push_back(o);
   std::stable_sort(hosts.begin(), hosts.end(), [](const ContactHost *a, const ContactHost *b) { return a->m_regOrder < b->m_regOrder; });
   hosts.erase(std::unique(hosts.begin(), hosts.end()), hosts.end());   // (the host found by both)
   for (ContactHost *const o : hosts)
   {
      const bool mesh = o->FinderObjIndex() >= 0;
      if (o != host)
      {
         if (!o->m_enabled || !o->m_set.contact)
            continue;
         if (o->m_set.band)
            continue;   // the band rubbers are the band model's: as a rigid surface contact they would stop the ball before the band could deflect (as Regather)
         const FRect3D &bb = o->m_hitBBox;
         if (!mesh && (posStart.x < bb.left - reach || posStart.x > bb.right + reach || posStart.y < bb.top - reach || posStart.y > bb.bottom + reach
             || posStart.z < bb.zlow - reach || posStart.z > bb.zhigh + reach))
            continue;
      }
      SurfaceContact it;
      it.r = o;
      it.rp = o->m_set.ToParams();
      it.hLayer = o->LayerM();
      it.st.Reset();
      it.touching = o == host && std::strcmp(how, "near") != 0;   // the stock collider that found the ball fires its own hit event; a surface found by proximity gets its event from the first touch
      items.push_back(it);
   }
   ContactSurfaceHost::FinderBind(pb, items, false);   // the mesh items read the engagement query's faces (their Refresh below does nothing; the inside check re-queries as the seat moves the ball)
}

double ContactSystemImpl::SeatNoPast(std::vector<SurfaceContact> &items, const BallSI &bs, V3 &c, double &delta0, int &passes, double &left, std::vector<const ContactSurface *> *const inside)
{
   // along the normal of each overlapped surface (the playfield is one of them) to touching distance, repeated up to 8 passes
   using namespace RubberContact;
   const V3 c0 = c;
   for (SurfaceContact &it : items)
      it.r->Refresh(it, bs.R, 0.03, c);
   for (; passes < 8; passes++)
   {
      bool moved = false;
      for (const SurfaceContact &it : items)
      {
         const TubeProbe p0 = it.r->Probe(bs, it, c);
         if (!p0.overlap || p0.delta <= 1e-6)
            continue;
         if (passes == 0)
            delta0 = std::max(delta0, p0.delta);
         c = c + p0.normal * p0.delta;
         moved = true;
         if (inside != nullptr && std::find(inside->begin(), inside->end(), it.r) == inside->end())
            inside->push_back(it.r);
      }
      if (!moved)
         break;
   }
   for (const SurfaceContact &it : items)
   {
      const TubeProbe p0 = it.r->Probe(bs, it, c);
      if (p0.overlap)
         left = std::max(left, p0.delta);
   }
   return (c - c0).Length();
}

// A ball enters when a surface's stock collider finds it (the stock impulse is skipped) and leaves when it is clear of every surface and receding.
bool ContactSystemImpl::EngageImpl(ContactHost *const host, HitBall *const pb, const char *how)
{
   CPROF(EngageImpl);
   const ContactProf::Stage profStage(ContactProf::StEngageImpl);   // P3a profile: the calling stage of the checks here (diagnostics only)
   using namespace RubberContact;
   if ((host != nullptr && !host->m_set.contact) || pb->m_bandOwner != nullptr || pb->m_contact != nullptr || KickerHold::Held(pb))
      return false;   // (a ball held in a kicker, parked or carried, belongs to the kicker: KickerDesign 1h)
   BallContact rb;
   rb.ball = pb;
   BallSI bs;
   bs.R = (double)pb->m_d.m_radius * U;
   const Vertex3Ds posStart = pb->m_d.m_pos;
   OwnLoad(pb);   // (the owned state: its floats only for a ball with no state yet)
   V3 c = OwnC(pb);
   const V3 v = OwnV(pb);
   GatherItems(host, pb, bs, c, how, rb.items);
   // The ball starts the engagement touching or a little inside the rubber (the hit was found within a step). Go back to the instant it first touched: along its
   // own velocity, by depth / closing speed, at most one physics step; if it is not closing (pushed in by something else), seat it at touching distance along the normal.
   // Nothing is stored as compression.
   // Found inside a part with no earlier clean state (docs/unified-physics/archive/UnifiedContactDesign.txt 11, point 6): CheckOverlaps found no clean start of the cycle, or the engagement by proximity finds the ball already inside a
   // moving part (point 1 engages a ball that has a past before it touches, so a ball inside at that point was created or moved there)
   bool noPast = std::strcmp(how, "no past") == 0 && !(host != nullptr && host->MovingPart() != nullptr && host->MovingPart()->GroupPart() != nullptr);   // (a part on the group step: its zero, below)
   if (ContactSystem::s_unified)   // one entry per ball nets (0 expected): a ball the overlap replay found, a ball with no past
      (std::strcmp(how, "no past") == 0 ? g_qe.noPast : g_qe.overlap)++;
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
      // Minimum translation (SeatNoPast). A ball that is squeezed and does not fit keeps the overlap that is left, and the force law carries it.
      noPastMoved = SeatNoPast(rb.items, bs, c, rb.delta0, noPastPasses, noPastLeft, nullptr);
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
   if (ContactSystem::s_unified)   // one entry per ball: no rewind to the touch; a ball found inside counts from zero (the group step's zero rule)
   {
      zeroNotSeat = false;
      for (const SurfaceContact &it : rb.items)
      {
         const TubeProbe p0 = it.r->Probe(bs, it, c);
         zeroNotSeat = zeroNotSeat || (p0.overlap && p0.delta > 1e-6);
      }
   }
   rb.seatPending = zeroNotSeat;
   if (!zeroNotSeat && !ContactSystem::s_unified)
      SeatAtTouch(rb.items, bs, c, v, rb.delta0, rb.rewound);
   }
   WriteBall(pb, c, v, OwnW(pb));   // (the seat: the engine's write)
   rb.vin = v;
   rb.win = OwnW(pb);
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

void ContactSystemImpl::Release(HitBall *const ball, const char *why, const bool remove)
{
   if (ball->m_contact == nullptr)
      return;
   const auto itE = std::find_if(g_engaged.begin(), g_engaged.end(), [ball](const BallContact &r) { return &r == ball->m_contact; });
   if (itE == g_engaged.end())
      return;
   while (itE->band != nullptr && ball->m_bandOwner == itE->band)   // (each band that holds it, the next one becoming the first)
   {
      ContactSurfaceHost *const first = itE->band;
      first->BandRelease(ball, why);   // clears itE->band and the ball's band owner
      if (itE->band == first)
         break;
   }
   itE->band = nullptr;
   itE->bandsMore.clear();
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
   if (ContactSystem::s_unified && !remove)   // one entry per ball: the entry stays for the ball's life; what it held of the old position is gone (a script move, a kicker)
   {
      if (!itE->items.empty())
         g_qe.cleared++;
      BallContact fresh;
      fresh.ball = ball;
      fresh.posStart = ball->m_d.m_pos;
      *itE = std::move(fresh);
      return;
   }
   ball->m_contact = nullptr;
   g_engaged.erase(itE);
}

void ContactSystem::StepCordAlone(ContactSurfaceHost *const band)
{
   ContactSystemImpl::CordAlone(band, g_cycleS);
}

bool ContactSystem::EngageBand(ContactSurfaceHost *const band, HitBall *const pb)
{
   BallContact *e = pb->m_contact;
   if (e != nullptr && e->band != nullptr && e->band != band)
      return false;   // another band has the ball
   if (e == nullptr && ContactSystem::s_unified && PairEligible(pb) && !KickerHold::Held(pb))
      e = &ContactSystemImpl::EnsureOne(pb);   // one entry per ball: the ball's own entry (not yet made this cycle), the band joins it
   if (e == nullptr)
   {
      if (!ContactSystemImpl::EngageImpl(nullptr, pb, "band"))
         return false;
      e = pb->m_contact;
   }
   e->band = band;
   return true;
}

void ContactSystemImpl::Regather(BallContact &e, const Vertex3Ds &pos, const float radiusU, const int part, const double cullTravel)
{
   const float reach = radiusU + 20.0f;   // as at engagement
   // A Verlet list: the surfaces within reach + 10 U are kept and collected again only after the ball has moved 5 U (or 10 ms have passed: a part can become collidable)
   constexpr float kMarginU = 10.0f, kMoveU = 5.0f;
   const uint32_t nowMs = g_pplayer->m_time_msec;
   const bool rescan = part != 2 && (part == 1 || !(e.regathered && nowMs - e.regatherMs < 10 && (pos - e.regatherPos).LengthSquared() < kMoveU * kMoveU));
   std::vector<SurfaceContact> next;
   // shared finder (design.md section 3): the static mesh hosts every cycle from the ball's query (a face within R + reach), the others from the Verlet list below;
   // kept in the registry's order, as the registry loop gave them
   thread_local std::vector<ContactHost *> found;
   if (part == 1)   // the others only: the static mesh items as they are (the query at the ball's energy travel replaces them after the grouping)
      found.clear();
   else
      ContactSurfaceHost::FinderHosts(e.ball, found);
   for (SurfaceContact &it : e.items)
      if (part != 0 && it.r != nullptr && (it.r->FinderObjIndex() >= 0) == (part == 1))
         next.push_back(std::move(it));   // (the part not gathered now, kept)
   for (ContactHost *const o : found)
   {
      if (!o->m_set.contact)
         continue;
      const auto old = std::find_if(e.items.begin(), e.items.end(), [o](const SurfaceContact &it) { return it.r == o; });
      if (old == e.items.end() && o->m_set.band)
         continue;   // (as below)
      if (old != e.items.end())
      {
         next.push_back(std::move(*old));
         continue;
      }
      SurfaceContact it;
      it.r = o;
      it.rp = o->m_set.ToParams();
      it.hLayer = o->LayerM();
      it.st.Reset();
      next.push_back(std::move(it));
   }
   if (!rescan && part != 2)   // the other surfaces: the set as it is (part 2: kept above)
      for (SurfaceContact &it : e.items)
         if (it.r != nullptr && it.r->FinderObjIndex() < 0)
            next.push_back(std::move(it));
   if (rescan && part == 1)   // the moving parts below are the grouping candidates (ReachLinks probes them, ReachEnergy groups them);
   {                          // GatherNonMesh keeps after the grouping those grouped with the ball, and the fixed ones by their box
      for (SurfaceContact &it : e.items)
         if (it.r != nullptr && it.r->FinderObjIndex() < 0 && it.r->MovingPart() == nullptr)
            next.push_back(std::move(it));
      if (cullTravel >= 0.0)   // the part cull (design 1c): a moving part is a candidate only when the ball's travel bound reaches its region over the cycle (CullDist less R); the
      {                        // others get no entry, no probe and no energy test this cycle (the table's moving parts only, in the registry's order: CullLists)
         const RubberContact::V3 c((double)pos.x * U, (double)pos.y * U, (double)pos.z * U);
         const double R = (double)radiusU * U;
         for (ContactHost *const o : g_cullParts)
         {
            if (!o->m_enabled || !o->m_set.contact)
               continue;
            const auto old = std::find_if(e.items.begin(), e.items.end(), [o](const SurfaceContact &it) { return it.r == o; });
            if (old == e.items.end() && o->m_set.band)
               continue;   // (a band's own rule, below)
            g_cullTests++;
            if (o->MovingPart()->CullDist(c) - R > cullTravel)
            {
               g_cullSkipped++;
               continue;
            }
            if (old != e.items.end())
            {
               next.push_back(std::move(*old));
               continue;
            }
            SurfaceContact it;
            it.r = o;
            it.rp = o->m_set.ToParams();
            it.hLayer = o->LayerM();
            it.st.Reset();
            next.push_back(std::move(it));
         }
      }
      else   // no travel bound (a new ball's first query, EngageNear; no bound on the table's energy, CullPrepare): every moving part
      for (ContactHost *const o : ContactHost::Registry())
      {
         if (o->FinderObjIndex() >= 0 || o->MovingPart() == nullptr || !o->m_enabled || !o->m_set.contact)
            continue;
         const auto old = std::find_if(e.items.begin(), e.items.end(), [o](const SurfaceContact &it) { return it.r == o; });
         if (old == e.items.end() && o->m_set.band)
            continue;   // (a band's own rule, below)
         if (old != e.items.end())
         {
            next.push_back(std::move(*old));
            continue;
         }
         SurfaceContact it;
         it.r = o;
         it.rp = o->m_set.ToParams();
         it.hLayer = o->LayerM();
         it.st.Reset();
         next.push_back(std::move(it));
      }
   }
   else if (rescan)
   {
   e.regathered = true;
   e.regatherPos = pos;
   e.regatherMs = nowMs;
   for (ContactHost *const o : ContactHost::Registry())
   {
      if (o->FinderObjIndex() >= 0)
         continue;   // a static mesh host: from the query, above
      if (!o->m_enabled || !o->m_set.contact)
         continue;
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
      SurfaceContact it;
      it.r = o;
      it.rp = o->m_set.ToParams();
      it.hLayer = o->LayerM();
      it.st.Reset();
      next.push_back(std::move(it));
   }
   }
   std::stable_sort(next.begin(), next.end(), [](const SurfaceContact &a, const SurfaceContact &b) { return static_cast<const ContactHost *>(a.r)->m_regOrder < static_cast<const ContactHost *>(b.r)->m_regOrder; });
   if (!next.empty() || ContactSystem::s_unified)   // never drop the whole set: the release test below decides when the ball is free (one entry per ball: an empty set is legal)
      e.items = std::move(next);
}

void ContactSystem::EngageNear(const double cycleS)
{
   ContactSystemImpl::EngageNear(cycleS);
}

void ContactSystem::Release(HitBall *const ball, const char *why)
{
   ContactSystemImpl::Release(ball, why);
}

void ContactSystem::Remove(HitBall *const ball)
{
   ContactSystemImpl::Release(ball, "ball deleted", true);
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
   if (ContactSystem::s_unified && g_qe.cycles > 0)   // one entry per ball counters at the table's end (diagnostics only)
   {
      QueryEngageLog("table end");
      g_qe = QueryEngageCounts();
   }
   // the host is going away: its engagements end without a log line (a ball can outlive the physics engine: it must not keep a pointer to the entry)
   for (BallContact &e : g_engaged)
   {
      e.bandsMore.erase(std::remove(e.bandsMore.begin(), e.bandsMore.end(), host), e.bandsMore.end());   // (a second band that holds it)
      if (e.band == host)
      {
         e.band = e.bandsMore.empty() ? nullptr : e.bandsMore.front();   // (another band that holds it becomes its first)
         if (!e.bandsMore.empty())
            e.bandsMore.erase(e.bandsMore.begin());
      }
   }
   if (ContactSystem::s_unified)   // one entry per ball: the host's items go, every entry stays (its ball's life)
   {
      for (BallContact &e : g_engaged)
      {
         if (e.host == host)
            e.host = nullptr;
         e.items.erase(std::remove_if(e.items.begin(), e.items.end(), [host](const SurfaceContact &it) { return it.r == host; }), e.items.end());
         e.lastHosts.erase(std::remove(e.lastHosts.begin(), e.lastHosts.end(), (const void *)static_cast<const ContactSurface *>(host)), e.lastHosts.end());
      }
      return;
   }
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
   uint64_t g_cyclesSeen = 0, g_stepCalls = 0, g_bodyCycles = 0;
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

void ContactSystem::NoteDriveWork(const HitBall *const ball, const double joules, const char *const what)
{
   static const void *owner = nullptr;
   static double total = 0.0;   // the run's drive work between cycles, J (from this table start)
   if (owner != g_pplayer)
   {
      owner = g_pplayer;
      total = 0.0;
   }
   total += joules;
   const char *const name = TestDriver::BallName(ball);
   const Vertex3Ds &p = ball->m_d.m_pos;
   PLOGI << "Drive between cycles: " << what << ", ball " << (name != nullptr ? name : "") << " at (" << p.x << ", " << p.y << ", " << p.z << "): work " << joules * 1e6 << " uJ (drive total "
         << total * 1e6 << " uJ), t=" << (g_pplayer != nullptr ? g_pplayer->m_time_msec : 0);
}

void ContactSystem::LoadBalls()
{
   if (g_pplayer == nullptr)
      return;
   for (Ball *const b : g_pplayer->m_vball)
      ContactSystemImpl::OwnLoad(&b->m_hitBall);
}

void ContactSystem::BallCentreU(const HitBall *const ball, double p[3])
{
   using namespace RubberContact;
   for (int i = 0; i < 3; i++)
      p[i] = ball->m_own.set ? ball->m_own.c[i] / U : (double)(&ball->m_d.m_pos.x)[i];
}

void ContactSystem::PushOutAtStart(HitBall *const ball, const std::string &name)
{
   if (!Unified() || !s_seatEditorBalls)
      return;
   using namespace RubberContact;
   BallSI bs;
   bs.R = (double)ball->m_d.m_radius * U;
   const Vertex3Ds p = ball->m_d.m_pos;
   const V3 c0((double)p.x * U, (double)p.y * U, (double)p.z * U);
   V3 c = c0;
   std::vector<SurfaceContact> items;
   ContactSystemImpl::GatherItems(nullptr, ball, bs, c, "start", items);
   double delta0 = 0.0, left = 0.0;
   int passes = 0;
   std::vector<const ContactSurface *> inside;
   const double moved = ContactSystemImpl::SeatNoPast(items, bs, c, delta0, passes, left, &inside);
   if (inside.empty())
      return;
   ball->m_d.m_pos = Vertex3Ds((float)(c.x / U), (float)(c.y / U), (float)(c.z / U));
   std::string names;
   for (const ContactSurface *const h : inside)
      names += (names.empty() ? "" : ", ") + h->Name();
   PLOGI << "Editor ball " << name << " was inside " << names << " at play start (" << delta0 * 1e6 << " um deep): moved out by the shortest way, " << moved * 1e3 << " mm in " << passes << " passes, "
         << left * 1e6 << " um overlap left; from " << p.x << ' ' << p.y << ' ' << p.z << " to " << ball->m_d.m_pos.x << ' ' << ball->m_d.m_pos.y << ' ' << ball->m_d.m_pos.z;
}

void ContactSystem::SleepCycleStart()
{
   ContactSystemImpl::SleepStart();
}

void ContactSystem::SleepWake(const HitBall *const ball, const char *const cause)
{
   ContactSystemImpl::SleepWakeBall(ball, cause);
}

// A step of dt seconds: the engine's whole cycle (Step), or a piece of it in an overlap replay (CheckOverlaps)
static void StepSeconds(const double dt)
{
   ContactSystemImpl::StepBalls(dt);
   if (dt > 0.0)
      g_placed.clear();
}

void ContactSystem::Step(const double cycleS)
{
   g_cycleS = cycleS;
   StepSeconds(cycleS);
}

double ContactSystem::CycleS()
{
   return g_cycleS;
}

const ContactSystem::CycleStats &ContactSystem::LastCycle()
{
   return g_last;
}


// ---- The snapshot of the start of a cycle (docs/unified-physics/archive/UnifiedContactDesign.txt section 11, point 2) ---------------------------------------------------------------------------------
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
      HitBall::OwnState own;   // the engine's owned state, with the floats it wrote (to-do 1.51): the pair goes back together
      int sleepGroup = 0;   // pile sleep
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
      s.balls[i].own = hb.m_own;
      s.balls[i].sleepGroup = hb.m_sleepGroup;
      s.balls[i].sleepHoldS = hb.m_sleepHoldS;
      s.balls[i].sleepCreepS = hb.m_sleepCreepS;
      s.balls[i].sleepPos = hb.m_sleepPos;
   }
   s.sleep = ContactSystemImpl::g_sleep;
   s.sleepNextId = ContactSystemImpl::g_sleepNextId;
   s.placed = g_placed;
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

// Puts the snapshot back, every ball
static void RestoreSnap()
{
   StartSnap &s = g_snap;
   if (!s.valid)
      return;
   ContactSystemImpl::g_sleep = s.sleep;   // pile sleep: every ball's, whatever the replay restores of the rest
   ContactSystemImpl::g_sleepNextId = s.sleepNextId;
   g_placed = s.placed;
   ContactSystemImpl::g_sleepCycleDt = s.sleepCycleDt;
   for (const BallSnap &b : s.balls)
   {
      b.pb->m_sleepGroup = b.sleepGroup;
      b.pb->m_sleepHoldS = b.sleepHoldS;
      b.pb->m_sleepCreepS = b.sleepCreepS;
      b.pb->m_sleepPos = b.sleepPos;
      b.pb->m_d.m_pos = b.pos;
      b.pb->m_d.m_vel = b.vel;
      b.pb->m_angularmomentum = b.angm;
      b.pb->m_lastEventSqrDist = b.lastEventSqrDist;
      b.pb->m_rubberAdvanced = b.advanced;
      b.pb->m_own = b.own;   // (with the floats above: as the engine left them, or as changed outside since; the next load tells which)
   }
   if (g_engaged.size() == s.contacts.size())
   {
      size_t k = 0;
      for (BallContact &e : g_engaged)
         e = s.contacts[k++];
   }
   else
   {
      if (ContactSystem::s_unified)
         g_qe.snap++;   // one entry per ball net (0 expected): the list is one entry per ball, it changes only with a ball added or removed
      PLOGW << "Snapshot: the engaged list changed since it was taken (" << s.contacts.size() << " -> " << g_engaged.size() << "): contacts not restored";
   }
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
      ContactSystemImpl::s_pendingCy = nullptr;
   }
}

// The play ends: one entry per ball keeps a ball's entry for its life and the player frees its balls without Remove, so without this a later play in the same process steps freed balls
// (runs/2026-11-16-crash: NaN, energy steps, a crash). Called once per play, at the end of PhysicsEngine's destructor: its hosts have gone (ReleaseHost has run for each), the player's balls are
// still alive. Only bookkeeping goes; no ball's place, velocity or spin is touched.
void ContactSystem::EndTable()
{
   size_t notInPlay = 0;   // entries whose ball is not in the player's list (should be 0: a ball leaving play is removed); compared by address only, never read
   for (const BallContact &e : g_engaged)
      if (g_pplayer == nullptr || std::none_of(g_pplayer->m_vball.begin(), g_pplayer->m_vball.end(), [&e](const Ball *b) { return &b->m_hitBall == e.ball; }))
         notInPlay++;
   PLOGI << "Contact system: play end: " << g_engaged.size() << " ball entries (" << notInPlay << " of balls not in play), " << g_ballPairs.size() << " ball pairs, "
         << ContactSystemImpl::g_sleep.size() << " sleep groups: cleared";
   if (g_pplayer != nullptr)
      for (Ball *const b : g_pplayer->m_vball)   // a ball that outlives the play (an editor ball): no pointer into the stores, no sleep group
      {
         b->m_hitBall.m_contact = nullptr;
         b->m_hitBall.m_sleepGroup = 0;
         b->m_hitBall.m_own.set = false;   // (an editor ball: the next play loads it from its floats)
      }
   g_engaged.clear();
   g_ballPairs.clear();
   g_placed.clear();
   g_deferred.clear();   // (events of the last step: their parts are gone)
   ContactSystemImpl::g_coilBooked.clear();
   ContactSystemImpl::s_ballOrderLast.clear();
   ContactSystemImpl::s_ballGroupLast.clear();
   ContactSystemImpl::s_pairReachLast.clear();
   ContactSystemImpl::s_pairReachNow.clear();
   ContactSystemImpl::g_sleep.clear();
   ContactSystemImpl::g_sleepOwner = nullptr;   // (SleepTable resets the rest at the next play, even if the new player has the old one's address)
   ContactSystemImpl::s_pendingCy = nullptr;
   ContactSystemImpl::s_nonFiniteTold = false;
   g_snap = StartSnap();   // (ball and part pointers)
}

void ContactSystem::RestoreStart()
{
   RestoreSnap();
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
         f.balls = FnvV(f.balls, hb.m_own.c);   // the owned state (to-do 1.51): what the next cycle steps from
         f.balls = FnvV(f.balls, hb.m_own.v);
         f.balls = FnvV(f.balls, hb.m_own.w);
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


// ---- The overlap check and the replay of the cycle (docs/unified-physics/archive/UnifiedContactDesign.txt section 11, point 2) ------------------------------------------------------------------------
// Run after UpdateDisplacements, before the events are dispatched. A free ball found overlapping a moving part was born inside this cycle: it moved as a body with no contacts
// while the part moved by its own motion. Find the time t* in the cycle at which the gap was zero with BOTH of them where they really were, put the snapshot of the start back, step the
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

// hittime: the cycle in the engine's time units (only to skip an empty one); hS: the same cycle in seconds
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
   auto isFree = [](const HitBall *pb) { return pb->m_contact == nullptr && pb->m_bandOwner == nullptr && !KickerHold::Held(pb); };
   struct Cand
   {
      HitBall *pb;
      ContactHost *host;
      double delta;   // the overlap at the end of the cycle, m
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
         const V3 c = ContactSystemImpl::OwnC(pb);
         for (ContactHost *const o : parts)
         {
            ContactMovingPart *const mp = o->MovingPart();
            // the part finder: only a part whose region where it is now is within the ball's radius (plus the tolerance) can overlap it
            const bool skip = mp->RegionDist(c) > bs.R + tol;
            GroupStep::PartFinderCount *const pc = GroupStep::partFinderProf ? &GroupStep::partFinderCount[typeid(*mp).name()] : nullptr;
            if (pc != nullptr)
            {
               pc->ovCalls++;
               pc->ovProbes += skip ? 0 : 1;
               pc->ovSkips += skip ? 1 : 0;
            }
            if (skip && !GroupStep::FinderCheckEnv())
               continue;
            mp->SetBall(bs);
            const TubeProbe pr = mp->Probe(c);
            if (skip)   // check mode: a skipped part must not overlap
            {
               if (pr.overlap && pr.delta > tol && pc != nullptr)
                  pc->ovMismatch++;
               continue;
            }
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
      double t;   // the share of the cycle at which the ball is engaged
   };
   std::vector<Rp> list;
   // The system at share tEnd of the cycle: the snapshot put back, the engaged balls and the parts stepped to each engagement, the balls of `list` engaged there, then on to tEnd
   auto simulate = [&](const double tEnd) {
      for (const Rp &e : list)
         if (e.pb->m_contact != nullptr)
            Release(e.pb, "replay");   // engaged by an earlier simulation: the snapshot has no entry for it
      RestoreSnap();   // every free ball is stepped again by StepSeconds below, so each starts from the snapshot
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
         for (ContactHost *const o : Hosts())   // a part no engaged ball held has not moved yet: it is at the start of the cycle. It moves on to the engagement, as it would have alone (gapAt measures it so)
            if (ContactMovingPart *const mp = o->MovingPart())
            {
               const double idle = hS * e.t - mp->steppedS;
               if (idle > 1e-9)
               {
                  mp->AdvanceAlone(idle);
                  mp->steppedS += idle;   // (PartAdvance advances it for the rest of the cycle only)
               }
            }
         // (the step above has already flown the ball to the engagement)
         EngageImpl(e.host, e.pb, "overlap");
      }
      if (tEnd > t)
      {
         s_append = !first;
         StepSeconds((tEnd - t) * hS);
         s_append = false;
      }
   };
   // The gap between a free ball (where simulate(t) stepped it) and the nearest part, at share t, in the system simulate(t) left. A part that was stepped for less than t is advanced alone
   // for the rest on a scratch copy of its state.
   auto gapAt = [&](HitBall *pb, const double t) {
      BallSI bs;
      bs.R = (double)pb->m_d.m_radius * U;
      const V3 c = ContactSystemImpl::OwnC(pb);
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
   // The share of the cycle at which the ball first touched: regula falsi on the gap with the system simulated at each trial time (the contact never starts deeper than the arrival depth)
   auto findT = [&](const Cand &cd, double &tOut) {
      simulate(0.0);
      evaluations++;
      const double g0 = gapAt(cd.pb, 0.0);
      if (g0 < -tol)
         return false;   // overlapping at the start of the cycle: no past
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
               PLOGW << "Overlap with " << cd.host->m_name << " at the START of the cycle (" << cd.delta * 1e6 << " um at its end): no past to rewind to; the ball is seated by minimum translation (point 6)";
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
      b->m_hitBall.m_rubberAdvanced = false;   // every ball has been moved for this cycle
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
            g_pplayer->m_liveUI->PushNotification("Several balls overlapped moving parts in one cycle and were replayed (see vpinball.log, 'Overlap replay')", 8000);
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
