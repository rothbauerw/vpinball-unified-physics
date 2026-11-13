// license:GPLv3+

#include "core/stdafx.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <sstream>
#include <array>
#include <chrono>
#include <map>
#include <mutex>
#include <unordered_map>
#include <intrin.h>

#include "imgui/imgui.h"
#include "parts/ball.h"
#include "physics/hitrubber.h"
#include "physics/cordband.h"
#include "physics/turntable.h"
#include "physics/playfieldcut.h"
#include "physics/contactsystem.h"
#include "physics/hitball.h"
#include "physics/magnet.h"
#include "physics/hitflipper.h"
#include "physics/physicsini.h"
#include "ui/live/LiveUI.h"
#include "physics/cabinet/NudgeHandler.h"
#include "physics/PhysicsEngine.h"
#include "parts/rubber.h"
#include "parts/primitive.h"

bool ContactTraceDetail()   // VPX_CONTACT_TRACE=2: also the per-cycle lines (ContactTrace, ContactPair): huge logs, the log file rotates; 1 = engage / release / part hit lines only
{
   static const bool on = []() { const char *const e = ContactSystem::EnvVar("VPX_CONTACT_TRACE", "VPX_RING_TRACE"); return e != nullptr && atoi(e) >= 2; }();
   return on;
}

bool ContactTraceOn()   // read once: getenv in a per-cycle path costs microseconds
{
   static const bool on = ContactSystem::EnvVar("VPX_CONTACT_TRACE", "VPX_RING_TRACE") != nullptr;
   return on;
}

namespace
{
   constexpr double U = ContactSystem::kMetresPerUnit;
   constexpr double VS = ContactSystem::kVelocityScale;
   constexpr double kMassUnitKg = ContactSystem::kMassUnitKg;
   using PhysicsIni::Trim;

   // Earliest time t >= 0 at which a point at p0 moving with velocity v is at distance r from the point q (-1 = none)
   double SphereTime(const RubberContact::V3 &p0, const RubberContact::V3 &v, const RubberContact::V3 &q, double r)
   {
      const RubberContact::V3 w = p0 - q;
      const double a = v.Dot(v);
      if (a < 1e-18)
         return -1.0;
      const double b = w.Dot(v);
      const double c = w.Dot(w) - r * r;
      const double disc = b * b - a * c;
      if (disc < 0.0)
         return -1.0;
      const double t = (-b - std::sqrt(disc)) / a;
      return t >= 0.0 ? t : -1.0;
   }
}

bool HitRubber::SurfaceHeightAt(const double xU, const double yU, double &zU) const
{
   if (m_surface == nullptr)
      return false;
   const double x = xU * U, y = yU * U;
   bool found = false;
   double best = -1e30;
   for (const RubberContact::SurfaceTri &t : m_surface->tris)
   {
      if (std::abs(t.n.z) < 0.5 || x < t.lo.x || x > t.hi.x || y < t.lo.y || y > t.hi.y)
         continue;
      // point in triangle (xy), then the height of its plane there
      const double d = (t.b.y - t.c.y) * (t.a.x - t.c.x) + (t.c.x - t.b.x) * (t.a.y - t.c.y);
      if (std::abs(d) < 1e-18)
         continue;
      const double l1 = ((t.b.y - t.c.y) * (x - t.c.x) + (t.c.x - t.b.x) * (y - t.c.y)) / d;
      const double l2 = ((t.c.y - t.a.y) * (x - t.c.x) + (t.a.x - t.c.x) * (y - t.c.y)) / d;
      const double l3 = 1.0 - l1 - l2;
      if (l1 < -1e-9 || l2 < -1e-9 || l3 < -1e-9)
         continue;
      const double z = l1 * t.a.z + l2 * t.b.z + l3 * t.c.z;
      if (z > best)
      {
         best = z;
         found = true;
      }
   }
   zU = best / U;
   return found;
}


// ---- The shared finder, piece P1 (to-do 1.34; runs/2026-11-11-finder/design.md sections 1, 2 and 5) ----
//
// The index is built once, before the first cycle, over every static mesh host (a wall, ramp, primitive or the playfield with faces): per object a box tree over its faces'
// tight boxes and the faces' mesh pieces (joined through shared edges), and one box tree over the objects' boxes. In P1 nothing in the engine reads it: only check mode
// (VPX_FINDER_CHECK=1) queries it and compares it, and today's face lists, with a full scan. Collidable is the host's own m_enabled, read at the query (a flag on the leaf).

namespace
{
   using FinderNode = RubberSurface::FinderNode;

   struct FinderObj
   {
      const HitRubber *host;
      RubberContact::V3 lo, hi;   // the object's box: its face tree's root, m
   };

   // P2a: one engaged ball's per-cycle query (design.md section 3), kept for the cycle; only the ball's own evaluations (one thread) write it after the gather
   struct FinderBallQ
   {
      RubberContact::V3 c0;      // the query's centre: the ball's centre at the cycle start, m
      double R = 0.0, reach = 0.0, ext = 0.0;   // ball radius; reach (ContactReachM); ext: how far the ball may stray from c0 inside the volume (reach, 2 x stray after a re-query)
      std::vector<int> objs;     // the objects found (index into the finder's objs), in load order
      std::vector<int> first;    // per found object k: its faces are faces[first[k] .. first[k + 1])
      std::vector<int> faces;    // the faces found, per object in load order
      std::vector<std::vector<int>> per;   // P2b: per finder object (index into the finder's objs), its faces found (empty: not found); items read it in place (SurfaceContact::qFaces)
      const std::vector<SurfaceContact> *items = nullptr;   // P2b: the ball's surfaces this cycle (FinderBind): a re-query counts the objects it finds that are not among them
      RubberContact::V3 vPrev;   // the ball's velocity at the last query (m/s): the acceleration is its change over one cycle
      uint64_t vPrevCycle = 0;
      bool queried = false;
      int reqThis = 0;           // re-queries in the current cycle
      double reqNs = 0.0;        // time spent in re-queries, ns (summed at the totals)
      uint64_t reqNew = 0;       // P2b: objects a re-query found that were not among the ball's surfaces (summed at the totals)
      uint64_t ver = 0;          // P3a: the query's version: each query, re-query and bind starts every table of the ball anew
      std::vector<RubberContact::FaceTable> tab;   // P3a: per finder object, the ball's measurement table of its faces in per[k] (one per ball position)
      const HitBall *ball = nullptr;   // P3a: the ball (check mode's log)
      PlayfieldCuts::CutQuery cq;   // P4b: the playfield cuts found by the query, pieces of their playfield host (NearCuts reads them; playfieldcut.h)
   };

   // the timing counters' histogram: 10 ns bins up to 200 us
   struct FinderHist
   {
      std::vector<uint32_t> b = std::vector<uint32_t>(20001, 0);
      uint64_t n = 0;
      double sum = 0.0, worst = 0.0;
      void Add(const double ns)
      {
         b[std::min((size_t)std::max(ns / 10.0, 0.0), b.size() - 1)]++;
         n++;
         sum += ns;
         worst = std::max(worst, ns);
      }
      double P(const double q) const
      {
         uint64_t acc = 0;
         const uint64_t need = (uint64_t)std::ceil(q * (double)n);
         for (size_t i = 0; i < b.size(); i++)
            if ((acc += b[i]) >= need)
               return (double)(i + 1) * 10.0;
         return worst;
      }
   };

   struct FinderIndex
   {
      bool built = false, on = false;
      std::vector<FinderObj> objs;     // in the host registry's order (load order)
      std::vector<FinderNode> nodes;   // the box tree over the objects
      std::vector<int> order;
      // P2a: the per-cycle queries and their counters
      std::unordered_map<const HitBall *, FinderBallQ> balls;   // node-based: an item's pointer to its ball's query stays valid as balls are added
      uint64_t cycle = 0, queries = 0, noPrev = 0, objsFound = 0, facesFound = 0, qCycles = 0, reqCycles = 0, reqTotal = 0;
      uint64_t queriesAt = 0, reqNewObj = 0;   // P2b: queries with no cycle bookkeeping (FinderQueryAt); objects a re-query found that were not among the ball's surfaces
      int reqWorst = 0;
      FinderHist qNs, qObjNs;          // query time per ball-cycle; per ball-cycle and object found
      uint32_t nextTotalsMs = 0;
      // check mode
      std::mutex mtx;
      uint64_t evals = 0, facesNear = 0, missA = 0, missB = 0, gathers = 0, gMissA = 0, gMissB = 0, missQ = 0, gMissQ = 0, gWider = 0, noQ = 0;
      std::map<std::string, std::array<uint64_t, 5>> byObj;   // object -> misses: list, index, gather, object tree, per-cycle query
      int logged = 0;
      bool notifiedA = false, notifiedB = false, notifiedG = false, notifiedQ = false;
      uint32_t nextReportMs = 0;   // check mode: the totals are also written once a second of table time (a run may end without closing the table)
      std::vector<FinderNode> cutNodes;   // P4b: the box tree over the playfield cuts' boxes (built at the first query after the cuts' Setup)
      std::vector<int> cutOrder;
      uint64_t cutVer = ~0ull;
      void Reset()   // the table closes (its first mesh host is deleted): the next table builds its own
      {
         built = on = false;
         objs.clear();
         nodes.clear();
         order.clear();
         balls.clear();
         cycle = queries = noPrev = objsFound = facesFound = qCycles = reqCycles = reqTotal = 0;
         queriesAt = reqNewObj = 0;
         reqWorst = 0;
         qNs = FinderHist();
         qObjNs = FinderHist();
         nextTotalsMs = 0;
         evals = facesNear = missA = missB = gathers = gMissA = gMissB = missQ = gMissQ = gWider = noQ = 0;
         byObj.clear();
         logged = 0;
         notifiedA = notifiedB = notifiedG = notifiedQ = false;
         nextReportMs = 0;
         cutNodes.clear();
         cutOrder.clear();
         cutVer = ~0ull;
      }
   };
   FinderIndex g_finder;

   // shared_finder P3a counters: one contact check = one outermost Probe, ProbePatch or ArrivalS of a static mesh host with a measurement table (key on)
   struct FinderTabStats
   {
      uint64_t checks = 0, fills = 0, reads = 0, mism = 0;          // the whole run
      uint64_t rowChecks = 0, rowFills = 0, rowReads = 0;           // since the cost log's last row (FinderTabTake)
      std::vector<uint64_t> perCyc = std::vector<uint64_t>(4097, 0);   // contact checks per ball, object and cycle (last bin: 4096 and more)
      uint64_t perCycWorst = 0;
      std::vector<uint64_t> ticks = std::vector<uint64_t>(65537, 0);   // time per contact check, 16-tick bins (TSC)
      uint64_t tickSum = 0, tickWorst = 0;
      uint64_t tsc0 = 0;
      std::chrono::steady_clock::time_point t0;
      int logged = 0;
      uint64_t passes = 0, nearVisit = 0, nearFills = 0, nearIn = 0, nearList = 0;   // P3b counters (check or profile mode): near passes made, faces box-tested and measured by them, faces handed over, query-list length
      uint64_t efills = 0, ereads = 0, emism = 0;   // P3c (check or profile mode): edges measured, edge reads, check-mode edge mismatches
      uint64_t cbMism = 0;   // cut B check mode: skippable calls whose full law ended with more than one peak
      uint64_t nearMiss = 0, probeDiff = 0, nearBox = 0;   // P3b check mode: faces within R the near pass did not hand over; probe answers that differ from the whole list's; box tests (counters)
      int loggedNear = 0;
      std::vector<uint64_t> cornerTicks = std::vector<uint64_t>(65537, 0);   // P3b: time per contact check on the rubber corners (design.md section 7), 16-tick bins
      uint64_t cornerSum = 0, cornerN = 0;
      struct ObjNear   // P3b2 (check or profile mode): per object, near passes made, query-list faces, faces past the box, the pass's ticks
      {
         uint64_t passes = 0, list = 0, visit = 0, ticks = 0;
         ContactProf::ObjBars bars;   // mc-cost profile (diagnostics only)
      };
      std::map<std::string, ObjNear> perObj;
   };
   FinderTabStats g_tab;
   thread_local int g_tabDepth = 0;
   thread_local bool g_nearOff = false;   // P3b check mode: the comparison's pass over the whole query list
   thread_local bool g_nearCmp = false;   // P3b check mode: inside one comparison (ProbePatch calls Probe)
   bool FinderTableOn()   // P3a's table on; VPX_FINDER_TABLE=0 (diagnostics only, results equal): the face loops measure every face directly
   {
      static const bool on = [] { const char *e = std::getenv("VPX_FINDER_TABLE"); return e == nullptr || e[0] != '0'; }();
      return on;
   }
   bool SameProbe(const RubberContact::TubeProbe &a, const RubberContact::TubeProbe &b)   // P3b check mode: the same bits
   {
      const auto eq = [](const double x, const double y) { return std::memcmp(&x, &y, sizeof(double)) == 0; };
      return a.overlap == b.overlap && a.seg == b.seg && eq(a.point.x, b.point.x) && eq(a.point.y, b.point.y) && eq(a.point.z, b.point.z) && eq(a.normal.x, b.normal.x) && eq(a.normal.y, b.normal.y)
         && eq(a.normal.z, b.normal.z) && eq(a.dist, b.dist) && eq(a.delta, b.delta) && eq(a.reff, b.reff) && eq(a.rho, b.rho) && eq(a.grad, b.grad);
   }
   // P3b: the faces a loop reads. in: the near pass's faces within R (load order, with their query-list positions); probe: the same, or with none within R the whole query list with the box bound
   // (the probe's nearest wherever it is, as today). No near pass (key off, check mode's whole-list pass): both the query list as P3a.
   void NearLists(RubberContact::FaceTable *const ft, const std::vector<int> &all, const bool nearOn, const std::vector<int> *&inL, RubberContact::FaceRead &inR, const std::vector<int> *&prL,
                  RubberContact::FaceRead &prR)
   {
      inR.t = prR.t = ft;
      inL = prL = &all;
      if (!nearOn || ft == nullptr || !ft->useNear)
         return;
      inL = &ft->inI;
      inR.pos = ft->inJ.data();
      if (!ft->inI.empty())
      {
         prL = inL;
         prR.pos = inR.pos;
      }
      else
         prR.prune = true;
   }

   uint64_t HistP(const std::vector<uint64_t> &h, const double q)   // the bin holding quantile q
   {
      uint64_t n = 0;
      for (const uint64_t b : h)
         n += b;
      const uint64_t need = std::max<uint64_t>((uint64_t)std::ceil(q * (double)n), 1);
      uint64_t acc = 0;
      for (size_t i = 0; i < h.size(); i++)
         if ((acc += h[i]) >= need)
            return i;
      return h.size() - 1;
   }

   thread_local int g_finderDepth = 0;   // check mode: one check per evaluation (ProbePatch calls Probe)
   struct FinderGuard
   {
      const bool on;
      explicit FinderGuard(const bool o) : on(o) { if (on) g_finderDepth++; }
      ~FinderGuard() { if (on) g_finderDepth--; }
   };

   bool FinderCheckOn()
   {
      static const bool on = [] { const char *e = std::getenv("VPX_FINDER_CHECK"); return e != nullptr && e[0] == '1'; }();
      return on;
   }

   inline double Comp(const RubberContact::V3 &v, const int a) { return a == 0 ? v.x : a == 1 ? v.y : v.z; }

   struct CheckProfScope   // P3a profile (diagnostics only, VPX_CHECK_PROFILE=1): one contact check's time, split in contactprof.h CheckData
   {
      const bool on;
      const int kind;
      bool outer = false;
      uint64_t t0 = 0, l0 = 0, f0 = 0, r0 = 0;
      const RubberContact::FaceTable *ft = nullptr;
      explicit CheckProfScope(const int k) : on(ContactProf::CheckOn()), kind(k)
      {
         if (!on)
            return;
         ContactProf::CheckData &d = ContactProf::CD();
         outer = d.depth++ == 0;
         if (outer)
         {
            l0 = d.loop;
            t0 = __rdtsc();
         }
      }
      void Tab(const RubberContact::FaceTable *t)
      {
         if (on && outer && t != nullptr)
         {
            ft = t;
            f0 = t->fills;
            r0 = t->reads;
         }
      }
      ~CheckProfScope()
      {
         if (!on)
            return;
         ContactProf::CheckData &d = ContactProf::CD();
         d.depth--;
         if (!outer)
            return;
         const uint64_t dt = __rdtsc() - t0;
         const int st = ContactProf::StageCur();
         d.st[st][kind]++;
         d.stT[st][kind] += dt;
         if (ft == nullptr)
         {
            d.nNoTab++;
            d.wholeNoTab += dt;
            return;
         }
         const double T = (double)dt, L = (double)(d.loop - l0), F = (double)(ft->fills - f0), R = (double)(ft->reads - r0);
         d.n++;
         d.whole += dt;
         d.fills += ft->fills - f0;
         d.reads += ft->reads - r0;
         d.sT += T; d.sL += L; d.sF += F; d.sR += R; d.sFF += F * F; d.sRR += R * R; d.sFR += F * R; d.sTF += T * F; d.sTR += T * R; d.sLF += L * F; d.sLR += L * R;
      }
   };

   struct TabCheck   // P3a: one contact check's counters (the outermost of nested calls: ProbePatch calls Probe)
   {
      RubberContact::FaceTable *const ft;
      const HitRubber *const host;
      const SurfaceContact &it;
      const RubberContact::V3 c;
      bool outer = false;
      uint64_t tk = 0, f0 = 0, r0 = 0, m0 = 0, ef0 = 0, er0 = 0, em0 = 0, cb0 = 0;
      TabCheck(RubberContact::FaceTable *t, const HitRubber *h, const SurfaceContact &i, const RubberContact::V3 &cc) : ft(t), host(h), it(i), c(cc)
      {
         if (ft == nullptr)
            return;
         outer = g_tabDepth++ == 0;
         if (!outer)
            return;
         f0 = ft->fills;
         r0 = ft->reads;
         m0 = ft->mism;
         ef0 = ft->efills;
         er0 = ft->ereads;
         em0 = ft->emism;
         cb0 = ft->cbMism;
         if (g_tab.tsc0 == 0)
         {
            ContactProf::CD().cyc0 = g_finder.cycle;   // P3a profile: the cycles the totals cover
            g_tab.tsc0 = __rdtsc();
            g_tab.t0 = std::chrono::steady_clock::now();
         }
         tk = __rdtsc();
      }
      ~TabCheck()
      {
         if (ft == nullptr)
            return;
         g_tabDepth--;
         if (!outer)
            return;
         const uint64_t dt = __rdtsc() - tk;
         if (ft->tag == 1)   // P3b: the rubber corners apart
         {
            g_tab.cornerTicks[std::min<size_t>((size_t)(dt / 16), g_tab.cornerTicks.size() - 1)]++;
            g_tab.cornerSum += dt;
            g_tab.cornerN++;
         }
         g_tab.ticks[std::min<size_t>((size_t)(dt / 16), g_tab.ticks.size() - 1)]++;
         g_tab.tickSum += dt;
         g_tab.tickWorst = std::max(g_tab.tickWorst, dt);
         const uint64_t df = ft->fills - f0, dr = ft->reads - r0;
         g_tab.checks++;
         g_tab.fills += df;
         g_tab.reads += dr;
         g_tab.rowChecks++;
         g_tab.rowFills += df;
         g_tab.rowReads += dr;
         if (ft->cyc != g_finder.cycle)
         {
            if (ft->nCyc > 0)
            {
               g_tab.perCyc[std::min<size_t>(ft->nCyc, g_tab.perCyc.size() - 1)]++;
               g_tab.perCycWorst = std::max<uint64_t>(g_tab.perCycWorst, ft->nCyc);
            }
            ft->cyc = g_finder.cycle;
            ft->nCyc = 0;
         }
         ft->nCyc++;
         g_tab.efills += ft->efills - ef0;   // P3c
         g_tab.ereads += ft->ereads - er0;
         if (ft->emism != em0)
         {
            g_tab.emism += ft->emism - em0;
            if (g_tab.logged < 50)
            {
               g_tab.logged++;
               PLOGW << "Finder edge slot mismatch: object " << host->Name() << ", cycle " << g_finder.cycle << ", c " << c.x << " " << c.y << " " << c.z;
            }
         }
         if (ft->cbMism != cb0)   // cut B check mode: a call the rule would skip whose full law gave more than one peak
         {
            g_tab.cbMism += ft->cbMism - cb0;
            if (g_tab.logged < 50)
            {
               g_tab.logged++;
               const FinderBallQ *const q = static_cast<const FinderBallQ *>(it.finderQ);
               PLOGW << "Cut B mismatch: object " << host->Name() << ", ball " << (const void *)(q != nullptr ? q->ball : nullptr) << ", cycle " << g_finder.cycle << ", c " << c.x << " " << c.y << " " << c.z;
            }
         }
         if (ft->mism != m0)   // check mode: a table read whose fresh measurement differs (bits of point, kind or d2)
         {
            g_tab.mism += ft->mism - m0;
            if (g_tab.logged < 50)
            {
               g_tab.logged++;
               const FinderBallQ *const q = static_cast<const FinderBallQ *>(it.finderQ);
               const int k = host->FinderObjIndex();
               const int face = (q != nullptr && k >= 0 && k < (int)q->per.size() && ft->lastMismJ >= 0 && ft->lastMismJ < (int)q->per[k].size()) ? q->per[k][ft->lastMismJ] : -1;
               PLOGW << "Finder table mismatch: object " << host->Name() << ", face " << face << ", ball " << (const void *)(q != nullptr ? q->ball : nullptr) << ", cycle " << g_finder.cycle << ", c "
                     << c.x << " " << c.y << " " << c.z;
            }
         }
      }
   };

   inline double BoxD2(const RubberContact::V3 &lo, const RubberContact::V3 &hi, const RubberContact::V3 &c)
   {
      const double dx = std::max(std::max(lo.x - c.x, c.x - hi.x), 0.0), dy = std::max(std::max(lo.y - c.y, c.y - hi.y), 0.0), dz = std::max(std::max(lo.z - c.z, c.z - hi.z), 0.0);
      return dx * dx + dy * dy + dz * dz;
   }

   // A box tree over items with boxes (box(i, lo, hi)): split at the median of the box centres along the longest side of the node's box, ties by item index; leafSize items per leaf.
   template <class BoxOf> int BuildTree(std::vector<FinderNode> &nodes, std::vector<int> &order, const BoxOf &box, const int first, const int count, const int leafSize)
   {
      FinderNode n;
      n.lo = RubberContact::V3(INFINITY, INFINITY, INFINITY);
      n.hi = RubberContact::V3(-INFINITY, -INFINITY, -INFINITY);
      for (int k = first; k < first + count; k++)
      {
         RubberContact::V3 lo, hi;
         box(order[k], lo, hi);
         n.lo = RubberContact::V3(std::min(n.lo.x, lo.x), std::min(n.lo.y, lo.y), std::min(n.lo.z, lo.z));
         n.hi = RubberContact::V3(std::max(n.hi.x, hi.x), std::max(n.hi.y, hi.y), std::max(n.hi.z, hi.z));
      }
      const int self = (int)nodes.size();
      nodes.push_back(n);
      if (count <= leafSize)
      {
         nodes[self].first = first;
         nodes[self].count = count;
         return self;
      }
      const RubberContact::V3 ext = n.hi - n.lo;
      const int axis = ext.x >= ext.y && ext.x >= ext.z ? 0 : ext.y >= ext.z ? 1 : 2;
      auto key = [&](const int i) { RubberContact::V3 lo, hi; box(i, lo, hi); return Comp(lo, axis) + Comp(hi, axis); };
      const int mid = first + count / 2;
      std::nth_element(order.begin() + first, order.begin() + mid, order.begin() + first + count, [&](const int a, const int b) { const double ka = key(a), kb = key(b); return ka < kb || (ka == kb && a < b); });
      const int l = BuildTree(nodes, order, box, first, mid - first, leafSize);
      const int r = BuildTree(nodes, order, box, mid, first + count - mid, leafSize);
      nodes[self].left = l;
      nodes[self].right = r;
      return self;
   }

   // The items whose box lies within r of c (r2 = r^2: an exact box distance, a true lower bound on the item's distance), appended to out in tree order
   template <class BoxOf> void QueryTree(const std::vector<FinderNode> &nodes, const std::vector<int> &order, const BoxOf &box, const RubberContact::V3 &c, const double r2, std::vector<int> &out)
   {
      if (nodes.empty())
         return;
      int stack[128];
      int sp = 0;
      stack[sp++] = 0;
      while (sp > 0)
      {
         const FinderNode &n = nodes[stack[--sp]];
         if (BoxD2(n.lo, n.hi, c) > r2)
            continue;
         if (n.left < 0)
         {
            for (int k = n.first; k < n.first + n.count; k++)
            {
               RubberContact::V3 lo, hi;
               box(order[k], lo, hi);
               if (BoxD2(lo, hi, c) <= r2)
                  out.push_back(order[k]);
            }
         }
         else if (sp + 2 <= 128)
         {
            stack[sp++] = n.right;
            stack[sp++] = n.left;
         }
      }
   }

   void FinderReport()
   {
      FinderIndex &g = g_finder;
      if (!g.on || !FinderCheckOn())
         return;
      std::ostringstream o;
      o << "Finder check totals: " << g.evals << " mesh evaluations (" << g.facesNear << " faces within the ball's radius), misses of today's face lists " << g.missA << ", of the index "
        << g.missB << ", of the per-cycle query " << g.missQ << " (evaluations with no query " << g.noQ << "); " << g.gathers << " ball gathers, objects with a face within the radius missing from the gather "
        << g.gMissA << ", from the index's object tree " << g.gMissB << "; gathered objects missing from the query " << g.gMissQ << " with a face within R + reach, " << g.gWider << " with none (the gather is wider)";
      PLOGI << o.str();
      for (const auto &e : g.byObj)
         PLOGI << "Finder check totals: " << e.first << ": misses list " << e.second[0] << ", index " << e.second[1] << ", gather " << e.second[2] << ", object tree " << e.second[3] << ", query " << e.second[4];
   }

   void FinderPeriodic()   // (the caller holds g_finder.mtx)
   {
      FinderIndex &g = g_finder;
      const uint32_t now = g_pplayer != nullptr ? g_pplayer->m_time_msec : 0;
      if (now < g.nextReportMs)
         return;
      if (g.nextReportMs != 0)
         FinderReport();
      g.nextReportMs = now + 1000;
   }

   void FinderMiss(const int what, const std::string &obj, const int face, const int piece, const double dM, const RubberContact::V3 &c, const char *where)
   {
      FinderIndex &g = g_finder;
      static const char *const names[5] = { "today's face list", "the shared finder's index", "the ball's gathered objects", "the index's object tree", "the per-cycle query" };
      g.byObj[obj][what]++;
      if (g.logged < 200)
      {
         g.logged++;
         PLOGW << "Finder check: " << names[what] << " missed " << obj << " face " << face << " (piece " << piece << ") at " << dM * 1e3 << " mm from the ball centre (" << c.x << ", " << c.y << ", " << c.z
               << ") m, t=" << (g_pplayer != nullptr ? g_pplayer->m_time_msec : 0) << " ms, " << where;
      }
      bool &said = what == 0 ? g.notifiedA : what == 1 ? g.notifiedB : what == 4 ? g.notifiedQ : g.notifiedG;
      if (!said)
      {
         said = true;
         const std::string msg = "Finder check: " + std::string(names[what]) + " missed a face of " + obj + " within the ball's radius (VPX_FINDER_CHECK; see vpinball.log 'Finder check')";
         PLOGW << msg;
         if (g_pplayer != nullptr && g_pplayer->m_liveUI != nullptr)
            g_pplayer->m_liveUI->PushNotification(msg, 5000);
      }
   }

   void FinderClose(const bool meshHost)
   {
      if (meshHost && g_finder.built)
      {
         HitRubber::FinderTotals(true);
         FinderReport();
         g_finder.Reset();
      }
   }
}

HitRubber::HitRubber(IEditable *editable, const std::vector<Vertex3Ds> &pts, const float radius, const RubberObjectSettings &settings, const std::filesystem::path &iniPath, const std::string &name)
   : ContactHost(editable, settings, iniPath, name, ""), m_pts(pts), m_radius(radius)
{
   if (!m_pts.empty())
      m_pts.push_back(m_pts.front()); // close the loop
   for (const Vertex3Ds &p : m_pts)
      m_ptsSI.emplace_back((double)p.x * U, (double)p.y * U, (double)p.z * U);
   m_radiusSI = (double)m_radius * U;
   RubberContact::BuildTubePre(m_tubePre, m_ptsSI.data(), (int)m_ptsSI.size());
   // A straight span makes it a band. A rubber without one wraps its post all the way round: it is a ring, and a ball on it takes the ring's law (ModelStandards 2d), never stock.
   // Collidable off at load: a ring too. A decorative rubber or a script's animation pose: no ball touches it, so it pays for no cord or mover. Turned on later by a script, it has the ring's law.
   if (static_cast<Rubber *>(editable)->m_d.m_collidable)
      DetectBand();
   m_set.band = !m_bandSpans.empty();
   if (!m_set.band && !m_set.slingPlate.empty())
   {
      const std::string msg = "Rubber " + m_name + ": the physics ini gives it a sling plate, but it has no straight span or its Collidable is off. Running it as a ring: check the [Sling:] section";
      PLOGW << msg;
      if (g_pplayer != nullptr && g_pplayer->m_liveUI != nullptr)
         g_pplayer->m_liveUI->PushNotification(msg, 5000);
      m_set.slingPlate.clear();
   }
}

bool HitRubber::FeatureOff() const
{
   const PhysicsIni::Features feat = PhysicsIni::Active(m_ini);
   return m_set.band ? !feat.bands : !feat.rubbers;
}

float HitRubber::BandSectionU() const
{
   return m_set.band && m_set.bandSectionMm > 0.0 ? (float)(m_set.bandSectionMm * 1e-3 / U) : 0.f;
}

void RubberSurface::Add(const Vertex3Ds &a, const Vertex3Ds &b, const Vertex3Ds &c)
{
   using RubberContact::V3;
   const Vertex3Ds *p[3] = { &a, &b, &c };
   if (tris.empty())
      for (int k = 0; k < 3; k++)
         lo[k] = hi[k] = (&a.x)[k];
   for (const Vertex3Ds *v : p)
      for (int k = 0; k < 3; k++)
      {
         lo[k] = std::min(lo[k], (&v->x)[k]);
         hi[k] = std::max(hi[k], (&v->x)[k]);
      }
   tris.push_back(RubberContact::MakeSurfaceTri(V3((double)a.x * U, (double)a.y * U, (double)a.z * U), V3((double)b.x * U, (double)b.y * U, (double)b.z * U),
                                                V3((double)c.x * U, (double)c.y * U, (double)c.z * U)));
   std::array<int, 3> ids;
   for (int k = 0; k < 3; k++)
   {
      std::array<uint32_t, 3> key;
      std::memcpy(key.data(), &p[k]->x, sizeof(float));
      std::memcpy(key.data() + 1, &p[k]->y, sizeof(float));
      std::memcpy(key.data() + 2, &p[k]->z, sizeof(float));
      ids[k] = weld.emplace(key, (int)weld.size()).first->second;
   }
   vid.push_back(ids);
   nb.clear();
}

int RubberSurface::Link()
{
   nb.assign(tris.size(), { -1, -1, -1 });
   std::map<std::pair<int, int>, std::pair<int, int>> edges;   // (lower id, higher id) -> (face, edge) of the first face on it
   int over = 0;
   for (int t = 0; t < (int)tris.size(); t++)
      for (int e = 0; e < 3; e++)
      {
         const int a = vid[t][e], b = vid[t][(e + 1) % 3];
         if (a == b)
            continue;
         const auto key = std::make_pair(std::min(a, b), std::max(a, b));
         const auto it = edges.find(key);
         if (it == edges.end())
            edges.emplace(key, std::make_pair(t, e));
         else if (it->second.first >= 0 && nb[it->second.first][it->second.second] < 0 && it->second.first != t)
         {
            nb[it->second.first][it->second.second] = t;
            nb[t][e] = it->second.first;
         }
         else
            over++;   // a third face on the edge (or one face twice): left unlinked
      }
   // cut B: each shared edge's fold, from both faces. Flat band (derived, not tuned): every stored coordinate is a float32, rounded by at most 2^-24 of its size, so each vertex is off by
   // at most dv = sqrt(3) 2^-24 M (M the largest coordinate size of the four vertices; the SI scale adds 2^-53, neglected); the height of a point above a plane through three such vertices
   // is then off by at most dv (1 + |l0| + |l1| + |l2|), l the point's barycentric coordinates in the face. A fold whose two heights are both inside their bounds is flat.
   using RubberContact::V3;
   fold.assign(tris.size(), { 0, 0, 0 });
   const auto height = [](const RubberContact::SurfaceTri &A, const V3 &p, double &band) {
      double M = 0.0;
      for (const V3 *q : { &A.a, &A.b, &A.c, &p })
         M = std::max({ M, std::abs(q->x), std::abs(q->y), std::abs(q->z) });
      const V3 v0 = A.b - A.a, v1 = A.c - A.a, v2 = p - A.a;
      const double d00 = v0.Dot(v0), d01 = v0.Dot(v1), d11 = v1.Dot(v1), d20 = v2.Dot(v0), d21 = v2.Dot(v1), den = d00 * d11 - d01 * d01;
      const double l1 = den != 0.0 ? (d11 * d20 - d01 * d21) / den : INFINITY, l2 = den != 0.0 ? (d00 * d21 - d01 * d20) / den : INFINITY, l0 = 1.0 - l1 - l2;
      band = std::sqrt(3.0) * std::ldexp(M, -24) * (1.0 + std::abs(l0) + std::abs(l1) + std::abs(l2));
      return (p - A.a).Dot(A.n);
   };
   for (int t = 0; t < (int)tris.size(); t++)
      for (int e = 0; e < 3; e++)
      {
         const int o = nb[t][e];
         if (o < t)   // each shared edge once, from its lower face
            continue;
         int f = -1;
         for (int k = 0; k < 3; k++)
            if (nb[o][k] == t)
               f = k;
         if (f < 0)
            continue;
         const RubberContact::SurfaceTri &A = tris[t], &B = tris[o];
         if (!A.area || !B.area)
         {
            fold[t][e] = fold[o][f] = 2;
            continue;
         }
         const V3 *va[3] = { &A.a, &A.b, &A.c }, *vb[3] = { &B.a, &B.b, &B.c };
         double bandA, bandB;
         const double hB = height(A, *vb[(f + 2) % 3], bandB);   // B's third vertex above A's plane
         const double hA = height(B, *va[(e + 2) % 3], bandA);   // A's third vertex above B's plane
         if (std::abs(hB) <= bandB && std::abs(hA) <= bandA)
            continue;   // flat
         fold[t][e] = hB > 0.0 ? 1 : hB < 0.0 ? -1 : 2;
         fold[o][f] = hA > 0.0 ? 1 : hA < 0.0 ? -1 : 2;
      }
   return over;
}

HitRubber::HitRubber(IEditable *editable, std::shared_ptr<RubberSurface> surface, const RubberObjectSettings &settings, const std::filesystem::path &iniPath, const std::string &name,
                     const char *kind)
   : ContactHost(editable, settings, iniPath, name, kind), m_surface(std::move(surface)), m_radius(0.f)
{
   // multi_contact: the faces' shared edges, built here at load (every face is added before the host is made), not at the first ball near the surface, where a large
   // mesh cost a cycle of milliseconds; built whatever the key, since nothing reads them with multi_contact 0
   if (m_surface != nullptr && m_surface->nb.size() != m_surface->tris.size())
      m_linkOver = m_surface->Link();
}

HitRubber::~HitRubber()
{
   FinderClose(m_surface != nullptr);   // the shared finder's index belongs to the table being closed
   // A ball can outlive its physics engine (a script still references it when the table closes): it must not keep a pointer to a band that is gone (the contact system lets
   // go of the host's balls in ~ContactHost)
   for (const BandBall &e : m_engaged)
      if (e.ball != nullptr && e.ball->m_bandOwner == this)
         e.ball->m_bandOwner = nullptr;
}

static std::vector<RubberContact::SurfaceHole> g_playfieldHoles;   // metres
bool PlayfieldHolesExist() { return !g_playfieldHoles.empty(); }

void HitRubber::AddPlayfieldHole(const float x, const float y, const float r)
{
   g_playfieldHoles.push_back({ (double)x * U, (double)y * U, (double)r * U });
}

void HitRubber::ClearPlayfieldHoles()
{
   g_playfieldHoles.clear();
}

// A ball whose centre is inside a non-legacy kicker's circle belongs to the kicker: VPX implements the kicker's funnel, capture and release together with the stock floor
// contact (the kicker's own impulses against the floor's embed correction), so the playfield model leaves it to them.
bool InPlayfieldHole(const Vertex3Ds &pos)
{
   for (const RubberContact::SurfaceHole &h : g_playfieldHoles)
   {
      const double dx = (double)pos.x * U - h.x, dy = (double)pos.y * U - h.y;
      if (dx * dx + dy * dy < h.r * h.r)
         return true;
   }
   return false;
}

// Band detection, for every rubber: find the posts and the free spans. With at least one free span the rubber is a band (the constructor), otherwise a ring, which is normal and says nothing.
// A VPX rubber stores only a closed centreline, so the posts are recovered from it: a vertex whose local bend radius is below
// band_curve_radius_mm is on a post wrap, a run of such vertices is one post (a circle is fitted to it), and the straight runs between
// consecutive posts are the spans. band_posts = "x,y,r | ..." replaces the fit with given circles (vertices within 15% of r of a circle are on it).
void HitRubber::DetectBand()
{
   m_bandPosts.clear();
   m_bandSpans.clear();
   const int n = (int)m_pts.size() - 1;   // the loop closed by repeating the first vertex
   if (n < 4)
      return;
   const double mmPerU = U * 1000.0;
   auto at = [&](int i) -> const Vertex3Ds & { return m_pts[((i % n) + n) % n]; };

   // posts given by the user, if any
   struct Circle { double x, y, r; };
   std::vector<Circle> given;
   {
      const std::string rest = m_set.bandPosts;
      size_t pos = 0;
      while (pos < rest.size())
      {
         size_t bar = rest.find('|', pos);
         if (bar == std::string::npos)
            bar = rest.size();
         const std::string item = rest.substr(pos, bar - pos);
         double a = 0, b = 0, c = 0;
         if (sscanf(item.c_str(), " %lf , %lf , %lf", &a, &b, &c) == 3 && c > 0.0)
            given.push_back({ a, b, c });
         else if (!Trim(item).empty())
            PLOGW << "Band " << m_name << ": ignored band_posts entry '" << Trim(item) << "' (want x,y,r)";
         pos = bar + 1;
      }
   }

   // which vertices are on a post wrap; owner = index into `given` (override mode), or -1
   std::vector<char> onPost(n, 0);
   std::vector<int> owner(n, -1);
   const double rCurve = m_set.bandCurveRadiusMm / mmPerU;   // VPX units
   for (int i = 0; i < n; i++)
   {
      if (!given.empty())
      {
         for (size_t k = 0; k < given.size(); k++)
         {
            const double d = std::hypot((double)at(i).x - given[k].x, (double)at(i).y - given[k].y);
            if (std::fabs(d - given[k].r) <= std::max(0.15 * given[k].r, 1.0))
            {
               onPost[i] = 1;
               owner[i] = (int)k;
               break;
            }
         }
      }
      else
      {
         const double ax = at(i).x - at(i - 1).x, ay = at(i).y - at(i - 1).y;
         const double bx = at(i + 1).x - at(i).x, by = at(i + 1).y - at(i).y;
         const double la = std::hypot(ax, ay), lb = std::hypot(bx, by);
         if (la < 1e-6 || lb < 1e-6)
            continue;
         const double turn = std::atan2(ax * by - ay * bx, ax * bx + ay * by);   // signed turning angle at the vertex
         const double k = std::fabs(turn) / (0.5 * (la + lb));                    // curvature, 1/U
         onPost[i] = (k > 1.0 / rCurve) ? 1 : 0;
      }
   }

   // runs of post vertices around the loop, starting just after a free vertex
   int start = -1;
   for (int i = 0; i < n; i++)
      if (!onPost[i])
      {
         start = i;
         break;
      }
   if (start < 0)
      return;   // every centreline vertex bends tighter than band_curve_radius_mm: no straight span
   for (int k = 1; k < n;)
   {
      const int i = start + k;
      if (!onPost[i % n])
      {
         k++;
         continue;
      }
      BandPost bp;
      bp.first = i;
      int j = i;
      while (j - i < n && onPost[(j + 1) % n] && (given.empty() || owner[(j + 1) % n] == owner[i % n]))
         j++;
      bp.last = j;
      const int cnt = j - i + 1;
      if (!given.empty())
      {
         bp.x = given[owner[i % n]].x;
         bp.y = given[owner[i % n]].y;
         bp.r = given[owner[i % n]].r;
      }
      else
      {
         // algebraic circle fit (Kasa) to the run, with the vertex on each side (the tangent points) included
         std::vector<std::pair<double, double>> q;
         for (int m = i - 1; m <= j + 1; m++)
            q.emplace_back(at(m).x, at(m).y);
         double mx = 0, my = 0;
         for (const auto &p : q) { mx += p.first; my += p.second; }
         mx /= (double)q.size();
         my /= (double)q.size();
         double suu = 0, suv = 0, svv = 0, suuu = 0, svvv = 0, suvv = 0, svuu = 0;
         for (const auto &p : q)
         {
            const double u = p.first - mx, v = p.second - my;
            suu += u * u; suv += u * v; svv += v * v; suuu += u * u * u; svvv += v * v * v; suvv += u * v * v; svuu += v * u * u;
         }
         const double det = suu * svv - suv * suv;
         bp.fitted = cnt >= 3 && std::fabs(det) > 1e-9 * (suu + svv) * (suu + svv);
         if (bp.fitted)
         {
            const double b1 = 0.5 * (suuu + suvv), b2 = 0.5 * (svvv + svuu);
            const double uc = (b1 * svv - b2 * suv) / det, vc = (b2 * suu - b1 * suv) / det;
            bp.x = mx + uc;
            bp.y = my + vc;
            bp.r = std::sqrt(uc * uc + vc * vc + (suu + svv) / (double)q.size());
         }
         else
         {
            bp.x = mx;
            bp.y = my;
            double rr = 0;
            for (const auto &p : q)
               rr += std::hypot(p.first - mx, p.second - my);
            bp.r = rr / (double)q.size();
         }
      }
      m_bandPosts.push_back(bp);
      k = j - start + 1;
   }
   if (m_bandPosts.size() < 2)   // a band needs at least 2
   {
      m_bandPosts.clear();
      return;
   }

   // spans: the vertices between one post run and the next
   const size_t np = m_bandPosts.size();
   for (size_t a = 0; a < np; a++)
   {
      const BandPost &pa = m_bandPosts[a];
      const BandPost &pb = m_bandPosts[(a + 1) % np];
      BandSpan sp;
      sp.postA = (int)a;
      sp.postB = (int)((a + 1) % np);
      sp.first = pa.last;                                         // the last post vertex is where the straight run leaves the wrap
      sp.last = pb.first + (pb.first <= pa.last ? n : 0);         // and the first post vertex of the next post is where it arrives
      double len = 0.0;
      for (int m = sp.first; m < sp.last; m++)
         len += std::hypot((double)at(m + 1).x - at(m).x, (double)at(m + 1).y - at(m).y);
      sp.lengthMm = len * mmPerU;
      sp.free = sp.lengthMm >= m_set.bandMinSpanMm;
      m_bandSpans.push_back(sp);
   }
   if (std::none_of(m_bandSpans.begin(), m_bandSpans.end(), [](const BandSpan &sp) { return sp.free; }))   // only short runs between tight bends: not a band
   {
      m_bandPosts.clear();
      m_bandSpans.clear();
      return;
   }

   PLOGI << "Band " << m_name << ": " << (given.empty() ? "fitted" : "given") << ' ' << np << " posts, " << m_bandSpans.size() << " spans, tube radius "
         << m_radius * mmPerU << " mm, " << n << " centreline points";
   for (size_t a = 0; a < np; a++)
      PLOGI << "  post " << a << ": centre (" << m_bandPosts[a].x << ", " << m_bandPosts[a].y << ") VPX units, centreline radius " << m_bandPosts[a].r << " U = "
            << m_bandPosts[a].r * mmPerU << " mm, post radius ~" << (m_bandPosts[a].r - m_radius) * mmPerU << " mm, vertices " << m_bandPosts[a].first << ".." << m_bandPosts[a].last
            << (m_bandPosts[a].fitted ? "" : " (too few vertices to fit a circle: centre is the mean)");
   for (size_t a = 0; a < m_bandSpans.size(); a++)
   {
      // how far the drawn run strays from its chord: a span the table draws curved shows here
      const BandSpan &sp = m_bandSpans[a];
      const double ax = at(sp.first).x, ay = at(sp.first).y, bx = at(sp.last).x - ax, by = at(sp.last).y - ay, chord = std::max(std::hypot(bx, by), 1e-6);
      double stray = 0.0;
      for (int m = sp.first; m <= sp.last; m++)
         stray = std::max(stray, std::fabs(((double)at(m).x - ax) * by - ((double)at(m).y - ay) * bx) / chord);
      PLOGI << "  span " << a << ": post " << sp.postA << " -> post " << sp.postB << ", vertices " << sp.first << ".." << sp.last << ", length " << sp.lengthMm << " mm, " << (sp.free ? "free" : "short") << ", off its chord by " << stray * mmPerU << " mm";
   }
}

// ---------------------------------------------------------------------------
// Sling: the plate and the switch blades are the table's own primitives. Nothing is written to the table file; the plate's pose is changed in memory.
void HitRubber::SlingAttach(Primitive *const plate, const std::vector<Primitive *> &switches)
{
   if (m_bandSpans.empty())
      return;   // a ring (see the constructor): nothing to put a plate behind
   m_slingPlate = plate;
   m_slingSwitches = switches;
   if (plate == nullptr)
      return;
   m_slingBase = plate->m_d.m_vPosition;
   m_slingShownRotX = plate->m_d.m_aRotAndTra[0];
   // the ball must pass through these: only the band touches it
   plate->put_Collidable(VARIANT_FALSE);
   m_switches.clear();
   for (Primitive *const sw : switches)
   {
      sw->put_Collidable(VARIANT_FALSE);
      BandSwitch ss;
      ss.prim = sw;
      ss.baseTransX = sw->m_d.m_aRotAndTra[3];   // the table's value is the blade's resting place; the engine slides it with the band
      ss.shownTransX = ss.baseTransX;
      ss.coil = true;
      ss.legacySign = true;
      m_switches.push_back(ss);
   }
   m_plate.Configure(m_set.sling, true);
   m_plateConfigured = true;
   m_driver.SetPulse(m_set.sling.pulseMs * 1e-3);
   m_driver.SetDelay(m_set.sling.coilDelayMs * 1e-3);
   PLOGI << "Sling " << m_name << ": plate " << plate->GetName() << " (pivot " << m_slingBase.x << ", " << m_slingBase.y << ", " << m_slingBase.z << "), " << switches.size()
         << " switch blades, collision off on all of them";
   SlingPlace();
}

void HitRubber::SlingPlace() const
{
   Primitive *const pl = m_slingPlate;
   if (pl == nullptr || m_bandSpans.empty())
      return;
   const int n = (int)m_pts.size() - 1;
   if (n < 4)
      return;
   // measure from the table's own pose, with the plate at rest (the engine may have it swung)
   const float shownRotX = pl->m_d.m_aRotAndTra[0];
   pl->m_d.m_aRotAndTra[0] = 0.f;
   pl->m_d.m_vPosition = m_slingBase;
   const Matrix3D &M = pl->RecalculateMatrices();
   std::vector<Vertex3Ds> w(pl->m_mesh.NumVertices());
   double cx = 0.0, cy = 0.0;
   for (size_t i = 0; i < w.size(); i++)
   {
      w[i] = M * pl->m_mesh.m_vertices[i];
      cx += w[i].x;
      cy += w[i].y;
   }
   if (w.empty())
      return;
   cx /= (double)w.size();
   cy /= (double)w.size();

   // the span the face looks at: for every free span, the face is the plate's surface at the band's height that lies along the span; its distance from the back of the tube is the gap. The span with the smallest gap wins.
   double bestGap = 1e30, nfx = 0.0, nfy = 0.0;
   int bestSpan = -1;
   std::vector<std::pair<double, double>> facePts, bestPts;   // (place along the chord, distance from the rest centreline) of the plate's surface at band height
   for (size_t k = 0; k < m_bandSpans.size(); k++)
   {
      const BandSpan &sp = m_bandSpans[k];
      if (!sp.free)
         continue;
      const Vertex3Ds &A = m_pts[((sp.first % n) + n) % n];
      const Vertex3Ds &B = m_pts[((sp.last % n) + n) % n];
      double tx = (double)B.x - A.x, ty = (double)B.y - A.y;
      const double len = std::sqrt(tx * tx + ty * ty);
      if (len < 1e-6)
         continue;
      tx /= len;
      ty /= len;
      double fx = ty, fy = -tx;                                   // a normal of the span
      if (((cx - A.x) * fx + (cy - A.y) * fy) > 0.0)              // make it point from the plate toward the band
      {
         fx = -fx;
         fy = -fy;
      }
      const double bandZ = 0.5 * ((double)A.z + B.z);
      double face = -1e30;
      // the plate's surface at the band's height: where the mesh edges cross the plane z = bandZ (the mesh may have vertices only at the ends of the face,
      // and the hood, above the band, never reaches this plane)
      facePts.clear();
      auto consider = [&](const Vertex3Ds &v) {
         const double u = ((double)v.x - A.x) * tx + ((double)v.y - A.y) * ty;
         if (u < 0.0 || u > len)
            return;
         const double dn = ((double)v.x - A.x) * fx + ((double)v.y - A.y) * fy;
         facePts.emplace_back(u, dn);
         face = std::max(face, dn);
      };
      const size_t ni = pl->m_mesh.NumIndices();
      for (size_t i = 0; i + 2 < ni; i += 3)
         for (int e = 0; e < 3; e++)
         {
            const Vertex3Ds &p0 = w[pl->m_mesh.m_indices[i + e]];
            const Vertex3Ds &p1 = w[pl->m_mesh.m_indices[i + (e + 1) % 3]];
            const double d0 = (double)p0.z - bandZ, d1 = (double)p1.z - bandZ;
            if (d0 * d1 > 0.0)
               continue;
            if (d0 == d1)
            {
               consider(p0);
               consider(p1);
               continue;
            }
            const double t = d0 / (d0 - d1);
            consider(Vertex3Ds(p0.x + (float)(t * (p1.x - p0.x)), p0.y + (float)(t * (p1.y - p0.y)), (float)bandZ));
         }
      if (face < -1e29 || face > 0.0)
         continue;
      const double gap = -face - m_radius;
      if (gap < bestGap)
      {
         bestGap = gap;
         bestSpan = (int)k;
         nfx = fx;
         nfy = fy;
         bestPts = facePts;
      }
   }
   if (bestSpan < 0)
   {
      PLOGW << "Sling " << m_name << ": plate " << pl->GetName() << " has no face behind a free span of the band at the band's height; the plate is left where the table has it";
      pl->m_d.m_aRotAndTra[0] = shownRotX;
      pl->RecalculateMatrices();
      return;
   }
   constexpr double mmPerU = U * 1000.0;
   m_slingGapTableMm = bestGap * mmPerU;
   const double shiftU = bestGap;   // the plate goes where its face touches the band at RotX = 0; the gap at rest is made by turning the plate back (below)
   pl->m_d.m_vPosition.x = m_slingBase.x + (float)(nfx * shiftU);
   pl->m_d.m_vPosition.y = m_slingBase.y + (float)(nfy * shiftU);
   pl->m_d.m_aRotAndTra[0] = 0.f;
   pl->RecalculateMatrices();
   // the switches sense the band where their blades are: the place along the span's chord
   {
      const BandSpan &sp = m_bandSpans[bestSpan];
      const Vertex3Ds &A = m_pts[((sp.first % n) + n) % n];
      const Vertex3Ds &B = m_pts[((sp.last % n) + n) % n];
      const double abx = (double)B.x - A.x, aby = (double)B.y - A.y, l2 = abx * abx + aby * aby;
      m_slingSpan = bestSpan;
      m_slingNx = nfx;
      m_slingNy = nfy;
      // the plate's geometry for the contact: the face (its surface within half a unit of the front) and the pivot axis, which runs through the plate's position
      const double faceRest = -m_radius;   // VPX units, from the rest centreline: after the move the face touches the back of the band at RotX = 0
      double uLo = 1e30, uHi = -1e30;
      for (const auto &pt : bestPts)
         if (pt.second + shiftU >= faceRest - 0.6)
         {
            uLo = std::min(uLo, pt.first);
            uHi = std::max(uHi, pt.first);
         }
      if (uHi > uLo && std::sqrt(abx * abx + aby * aby) > 1e-6)
      {
         const double len = std::sqrt(abx * abx + aby * aby);
         const double originN = (((double)pl->m_d.m_vPosition.x - A.x) * nfx + ((double)pl->m_d.m_vPosition.y - A.y) * nfy);
         m_plateOriginN = originN * U;
         m_plateA = (faceRest - originN) * U;
         m_plateH = (0.5 * ((double)A.z + B.z) - (double)pl->m_d.m_vPosition.z) * U;
         m_plateF0 = uLo / len;
         m_plateF1 = uHi / len;
         m_plateGeomValid = m_plateA > 0.0 && m_plateH > 0.0;
         // the gap at rest: the plate stands turned back from the touching pose by the angle at which the face is slingGapMm from the band
         //   a (1 / cos t - 1) + h tan t = -gap   for t < 0   (monotonic in t as h > a)
         m_plateRest = 0.0;
         if (m_plateGeomValid && m_set.slingGapMm > 0.0)
         {
            const double gapM = m_set.slingGapMm * 1e-3;
            double lo = -1.2, hi = 0.0;
            for (int it = 0; it < 60; it++)
            {
               const double mid = 0.5 * (lo + hi);
               const double g = m_plateA * (1.0 / std::cos(mid) - 1.0) + m_plateH * std::tan(mid);
               if (g > -gapM)
                  hi = mid;
               else
                  lo = mid;
            }
            m_plateRest = 0.5 * (lo + hi);
         }
         const float rx = (float)((m_plate.PhiDeg() + m_plateRest * 180.0 / M_PI));
         pl->m_d.m_aRotAndTra[0] = rx;
         m_slingShownRotX = rx;
         pl->RecalculateMatrices();
      }
      for (BandSwitch &sw : m_switches)
      {
         if (!sw.coil)
            continue;
         sw.span = bestSpan;
         sw.bx = -nfx;   // behind the band is toward the plate
         sw.by = -nfy;
         // the blade's place in the table: its position moved by the table's TransX along its own X axis (turned by its RotZ)
         const double rz = (double)sw.prim->m_d.m_aRotAndTra[2] * M_PI / 180.0;
         const double px = (double)sw.prim->m_d.m_vPosition.x + std::cos(rz) * sw.baseTransX, py = (double)sw.prim->m_d.m_vPosition.y + std::sin(rz) * sw.baseTransX;
         sw.f = std::min(std::max(((px - A.x) * abx + (py - A.y) * aby) / std::max(l2, 1e-12), 0.0), 1.0);
      }
   }
   if (std::fabs(m_slingPlacedGapMm - m_set.slingGapMm) > 1e-9)
   {
      m_slingPlacedGapMm = m_set.slingGapMm;
      PLOGI << "Sling " << m_name << ": plate " << pl->GetName() << " faces span " << bestSpan << ", face-to-band gap in the table " << m_slingGapTableMm << " mm, set to " << m_set.slingGapMm
            << " mm: the plate is moved " << std::fabs(shiftU) * mmPerU << (shiftU >= 0.0 ? " mm toward" : " mm away from") << " the band (" << nfx << ", " << nfy << ") so that its face touches the band at RotX = 0, pivot now ("
            << pl->m_d.m_vPosition.x << ", " << pl->m_d.m_vPosition.y << "), and it rests turned back to RotX = " << m_plateRest * 180.0 / M_PI << " deg for the gap";
   }
}

// Switch blades behind the band that only the script reads: collision off, the span they sense and the side they are on found from where they stand.
void HitRubber::BandSwitchesAttach(const std::vector<Primitive *> &blades)
{
   if (m_bandSpans.empty())
      return;
   for (Primitive *const prim : blades)
   {
      prim->put_Collidable(VARIANT_FALSE);
      BandSwitch ss;
      ss.prim = prim;
      ss.baseTransX = prim->m_d.m_aRotAndTra[3];
      ss.shownTransX = ss.baseTransX;
      SwitchPlace(ss);
      m_switches.push_back(ss);
   }
   PLOGI << "Band " << m_name << ": " << blades.size() << " switch blades for the script, " << m_switches.size() << " in all, collision off on them";
}

void HitRubber::SwitchPlace(BandSwitch &sw) const
{
   const int n = (int)m_pts.size() - 1;
   if (n < 4 || m_bandSpans.empty())
      return;
   // the blade's place in the table: its position moved by the table's TransX along its own X axis (turned by its RotZ)
   const double rz = (double)sw.prim->m_d.m_aRotAndTra[2] * M_PI / 180.0;
   const double px = (double)sw.prim->m_d.m_vPosition.x + std::cos(rz) * sw.baseTransX, py = (double)sw.prim->m_d.m_vPosition.y + std::sin(rz) * sw.baseTransX;
   int best = -1;
   double bestD = 1e30, bestT = 0.0, bnx = 0.0, bny = 0.0, bestS = 0.0;
   for (size_t k = 0; k < m_bandSpans.size(); k++)
   {
      const BandSpan &sp = m_bandSpans[k];
      if (!sp.free)
         continue;
      const Vertex3Ds &A = m_pts[((sp.first % n) + n) % n];
      const Vertex3Ds &B = m_pts[((sp.last % n) + n) % n];
      const double abx = (double)B.x - A.x, aby = (double)B.y - A.y, l2 = abx * abx + aby * aby;
      if (l2 < 1e-9)
         continue;
      const double t = std::min(std::max(((px - A.x) * abx + (py - A.y) * aby) / l2, 0.0), 1.0);
      const double dx = px - (A.x + abx * t), dy = py - (A.y + aby * t);
      const double d = std::sqrt(dx * dx + dy * dy);
      if (d < bestD)
      {
         const double len = std::sqrt(l2);
         best = (int)k;
         bestD = d;
         bestT = t;
         bnx = aby / len;
         bny = -abx / len;
         bestS = (px - A.x) * bnx + (py - A.y) * bny;
      }
   }
   if (best < 0)
   {
      PLOGW << "Band " << m_name << ": switch blade " << sw.prim->GetName() << " has no free span to sense";
      return;
   }
   const double side = bestS >= 0.0 ? 1.0 : -1.0;
   if (std::fabs(bestS) < 1.0)
      PLOGW << "Band " << m_name << ": switch blade " << sw.prim->GetName() << " stands on the band's chord (" << bestS * U * 1e3 << " mm off it): which side is behind is a guess";
   sw.span = best;
   sw.f = bestT;
   sw.bx = bnx * side;
   sw.by = bny * side;
   sw.axisProj = std::cos(rz) * sw.bx + std::sin(rz) * sw.by;
   PLOGI << "Band " << m_name << ": switch blade " << sw.prim->GetName() << " senses span " << best << " at " << bestT << " of its chord, " << std::fabs(bestS) * U * 1e3
         << " mm behind it, back direction (" << sw.bx << ", " << sw.by << "), blade X axis along it " << sw.axisProj;
}

// The blades follow the band: pushed back toward a blade, the blade slides with it; when the band stretches toward the ball the blade stays at its resting place.
void HitRubber::SwitchesFinish()
{
   for (BandSwitch &sw : m_switches)
   {
      const double push = std::max(sw.dispM, 0.0) / U;   // VPX units
      const float tx = (!sw.legacySign && std::fabs(sw.axisProj) >= 0.2) ? sw.baseTransX + (float)(push / sw.axisProj) : sw.baseTransX - (float)push;
      if (std::fabs(tx - sw.shownTransX) > 1e-3f)
      {
         sw.shownTransX = tx;
         sw.prim->m_d.m_aRotAndTra[3] = tx;
         sw.prim->RecalculateMatrices();
      }
   }
}

// the table's primitives follow the plate and the switches
void HitRubber::SlingFinish()
{
   m_slingInfo.enabled = m_slingEnabled;
   m_slingInfo.plateAngleDeg = (float)m_plate.PhiDeg();
   m_slingInfo.switchClosed = m_slingAnyClosed;
   m_slingInfo.coilOn = m_driver.Coil();
   m_slingInfo.coilCurrentA = (float)m_plate.CurrentA();
   m_slingInfo.pulses = m_driver.fired;
   const float rotX = (float)(m_plate.PhiDeg() + m_plateRest * 180.0 / M_PI);   // RotX = 0 is the pose in which the face touches the band
   if (std::fabs(rotX - m_slingShownRotX) > 1e-4f)
   {
      m_slingShownRotX = rotX;
      m_slingPlate->m_d.m_aRotAndTra[0] = rotX;
      m_slingPlate->RecalculateMatrices();
   }
   SwitchesFinish();
}

void HitRubber::CalcHitBBox()
{
   if (m_surface != nullptr)
   {
      const float m = 1.0f;
      m_hitBBox = FRect3D(m_surface->lo[0] - m, m_surface->hi[0] + m, m_surface->lo[1] - m, m_surface->hi[1] + m, m_surface->lo[2] - m, m_surface->hi[2] + m);
      if (const char *nearArg = ContactSystem::EnvVar("VPX_CONTACT_NEAR", "VPX_RING_NEAR"); nearArg != nullptr && m_set.contact)   // "x,y,r" (VPX units): the faces within r, as a 4 U grid of cells (count, z range)
      {
         double cx = 0, cy = 0, cr = 0;
         if (sscanf(nearArg, "%lf,%lf,%lf", &cx, &cy, &cr) == 3)
         {
            std::map<std::pair<int, int>, std::array<double, 3>> cells;   // count, zmin, zmax
            constexpr double U = 0.53975e-3;
            for (const RubberContact::SurfaceTri &t : m_surface->tris)
            {
               const double x = (t.a.x + t.b.x + t.c.x) / (3.0 * U), y = (t.a.y + t.b.y + t.c.y) / (3.0 * U);
               if (std::hypot(x - cx, y - cy) > cr)
                  continue;
               auto &c = cells.try_emplace({ (int)std::floor(x / 4.0), (int)std::floor(y / 4.0) }, std::array<double, 3>{ 0.0, 1e9, -1e9 }).first->second;
               c[0] += 1.0;
               c[1] = std::min(c[1], std::min({ t.a.z, t.b.z, t.c.z }) / U);
               c[2] = std::max(c[2], std::max({ t.a.z, t.b.z, t.c.z }) / U);
            }
            for (const RubberContact::SurfaceTri &t : m_surface->tris)   // DIAG: the vertices (VPX units) of every face whose bounding box reaches the circle
               if (t.lo.x / U <= cx + cr && t.hi.x / U >= cx - cr && t.lo.y / U <= cy + cr && t.hi.y / U >= cy - cr)
                  PLOGI << "RubberNearTri " << m_name << " a " << t.a.x / U << ' ' << t.a.y / U << ' ' << t.a.z / U << " b " << t.b.x / U << ' ' << t.b.y / U << ' ' << t.b.z / U << " c " << t.c.x / U << ' ' << t.c.y / U << ' ' << t.c.z / U;
            for (const auto &kv : cells)
               PLOGI << "RubberNear " << m_name << " cell " << kv.first.first * 4 << ' ' << kv.first.second * 4 << " faces " << kv.second[0] << " z " << kv.second[1] << ' ' << kv.second[2];
         }
      }
      if (ContactTraceOn() && m_set.contact)
         PLOGI << "RubberGeom " << m_name << " " << m_surfaceKind << " x " << m_surface->lo[0] << " " << m_surface->hi[0] << " y " << m_surface->lo[1] << " " << m_surface->hi[1] << " z " << m_surface->lo[2] << " " << m_surface->hi[2];
      return;
   }
   if (m_pts.empty())
      return;
   float x0 = FLT_MAX, x1 = -FLT_MAX, y0 = FLT_MAX, y1 = -FLT_MAX, z0 = FLT_MAX, z1 = -FLT_MAX;
   for (const Vertex3Ds &p : m_pts)
   {
      x0 = std::min(x0, p.x); x1 = std::max(x1, p.x);
      y0 = std::min(y0, p.y); y1 = std::max(y1, p.y);
      z0 = std::min(z0, p.z); z1 = std::max(z1, p.z);
   }
   m_hitBBox = FRect3D(x0 - m_radius, x1 + m_radius, y0 - m_radius, y1 + m_radius, z0 - m_radius, z1 + m_radius);
   if (ContactTraceOn() && m_set.contact)
      PLOGI << "RubberGeom " << m_name << " ring x " << x0 << " " << x1 << " y " << y0 << " " << y1 << " z " << z0 << " " << z1 << " radius " << m_radius << " points " << m_pts.size();
}

float HitRubber::HitTest(const BallS &ball, const float dtime, CollisionEvent &coll) const
{
   if (!m_enabled || m_pts.size() < 2)
      return -1.0f;
   if (CordModel())
      return -1.0f;   // the cord never meets the stock collision: the contact system takes the ball before it touches (CordNear; GroupStepDesign 15)
   if (ContactOwns(ball))
      return -1.0f;   // the contact model owns this ball until it is released

   using RubberContact::V3;
   const size_t np = m_pts.size();
   auto pt = [&](size_t k) -> V3 { return V3(m_pts[k].x, m_pts[k].y, m_pts[k].z); };
   const double rt = (double)ball.m_radius + (double)m_radius;      // ball centre to centreline at touch, VPX units
   const V3 c(ball.m_pos.x, ball.m_pos.y, ball.m_pos.z);
   const V3 v(ball.m_vel.x, ball.m_vel.y, ball.m_vel.z);

   double bestT = -1.0;
   bool bestContact = false;
   double bestBnv = 0.0, bestBnd = 0.0;

   for (size_t i = 0; i + 1 < np; i++)
   {
      const V3 A = pt(i), B = pt(i + 1);
      const V3 seg = B - A;
      const double L2 = seg.Dot(seg);
      if (L2 <= 1e-12)
         continue;
      const double L = std::sqrt(L2);
      const V3 u = seg / L;

      // distance now, to the nearest point of this segment
      const double s0 = std::min(std::max((c - A).Dot(seg) / L2, 0.0), 1.0);
      const V3 P0 = A + seg * s0;
      const double d0 = (c - P0).Length();
      if (d0 <= 1.0e-6)
         continue;
      const double bnv = (c - P0).Dot(v) / d0;                    // ball normal velocity
      if (bnv > (double)C_CONTACTVEL)
         continue;                                                // clearly receding
      const double bnd = d0 - rt;                                 // ball distance to the rubber surface

      double t;
      bool isContact = false;
      if (bnd < (double)PHYS_TOUCH)                               // already touching
      {
         if (std::fabs(bnv) <= (double)C_CONTACTVEL)
         {
            isContact = true;
            t = 0.0;
         }
         else
            t = std::max(0.0, -bnd / bnv);
      }
      else
      {
         // swept test of the ball centre against the capsule of radius rt around the segment: side of the cylinder, then both end caps
         t = -1.0;
         const V3 w = c - A;
         const V3 wp = w - u * w.Dot(u), vp = v - u * v.Dot(u);
         const double a = vp.Dot(vp);
         if (a > 1e-18)
         {
            const double b = wp.Dot(vp), cc = wp.Dot(wp) - rt * rt, disc = b * b - a * cc;
            if (disc >= 0.0)
            {
               const double tc = (-b - std::sqrt(disc)) / a;
               if (tc >= 0.0)
               {
                  const double s = ((w + v * tc).Dot(u)) / L;
                  if (s >= 0.0 && s <= 1.0)
                     t = tc;
               }
            }
         }
         for (const V3 *q : { &A, &B })
         {
            const double te = SphereTime(c, v, *q, rt);
            if (te >= 0.0 && (t < 0.0 || te < t))
               t = te;
         }
         if (t < 0.0)
            continue;
      }
      if (!std::isfinite(t) || t > (double)dtime)
         continue;
      if (bestT < 0.0 || t < bestT)
      {
         bestT = t;
         bestContact = isContact;
         bestBnv = bnv;
         bestBnd = bnd;
      }
   }
   if (bestT < 0.0)
      return -1.0f;

   // normal at the hit position: from the nearest centreline point
   const V3 hp = c + v * bestT;
   double bd = 1e30;
   V3 bp;
   for (size_t i = 0; i + 1 < np; i++)
   {
      const V3 A = pt(i), B = pt(i + 1);
      const V3 seg = B - A;
      const double L2 = seg.Dot(seg);
      if (L2 <= 1e-12)
         continue;
      const double s = std::min(std::max((hp - A).Dot(seg) / L2, 0.0), 1.0);
      const V3 P = A + seg * s;
      const double d = (hp - P).Length();
      if (d < bd)
      {
         bd = d;
         bp = P;
      }
   }
   if (bd <= 1e-9)
      return -1.0f;
   const V3 n = (hp - bp) / bd;
   coll.m_hitnormal = Vertex3Ds((float)n.x, (float)n.y, (float)n.z);
   coll.m_isContact = bestContact;
   if (bestContact)
      coll.m_hit_org_normalvelocity = (float)bestBnv;
   coll.m_hitdistance = (float)bestBnd;
   return (float)bestT;
}

// What the cord reads from the rubber's settings that can change while it plays, at load and at a reload. Its geometry (section, ring, stops, e_mpa, the sling) waits for the next play
static void FixedCordValues(CordFixed::Part &p, const RubberObjectSettings &s)
{
   p.mu = s.mu;
   p.tanD = s.tanDelta;
   p.refHz = s.bandRefHz;
   p.stretchLaw.alpha = s.bandStretchAlpha;
   p.twistLoss = s.bandTwistLoss;
   p.twistScale = s.bandTwistScale;
   p.ringRp = s.ToParams();
   p.SetPatchLaw(p.ringRp);
}

void HitRubber::Reload() const
{
   RubberObjectSettings s;
   if (m_surface != nullptr && m_surfaceKind == "playfield")
   {
      if (!LoadPlayfieldSettings(s, m_ini, m_tableMu))
         s.contact = false;
      m_set = s;
      PLOGI << "Playfield settings reloaded from " << m_ini.string() << ": contact " << (s.contact ? 1 : 0) << ", e_mpa " << s.eMpa << ", alpha " << s.alphaSPerM << ", tan_delta " << s.tanDelta << ", mu " << s.mu;
      return;
   }
   if (m_surface != nullptr)
   {
      if (m_hardSurface)
         LoadHardSurfaceSettings(s, m_ini, m_surfaceKind.c_str(), m_name, m_partE, m_partMu);
      else
         LoadSurfaceSettings(s, m_ini, m_surfaceKind.c_str(), m_name);   // contact follows the object's  rubber  flag, so it can be switched off while running
      m_set = s;
      PLOGI << "Rubber surface settings for " << m_name << " reloaded from " << m_ini.string() << ": rubber " << (s.contact ? 1 : 0) << ", alpha " << s.alphaSPerM << ", mu " << s.mu;
      return;
   }
   LoadSettings(s, m_ini, m_name);
   s.band = m_set.band;   // from the shape, at load
   if (!s.band)   // a rubber with no straight span stays a ring whatever a reloaded section says
      s.slingPlate.clear();
   m_set = s;   // contact on/off is decided when the table is loaded; the values apply immediately
   m_set.contact = true;
   SlingPlace();
   if (m_plateConfigured)
      const_cast<HitRubber *>(this)->m_plate.Configure(m_set.sling, false);
   PLOGI << "Rubber settings for " << m_name << " reloaded from " << m_ini.string() << ": alpha " << s.alphaSPerM << ", mu " << s.mu << ", tan delta " << s.tanDelta;
   if (m_cord != nullptr)
   {
      FixedCordValues(*m_cord->fixed, m_set);
      PLOGI << "Band " << m_name << ": values now: stretch alpha " << m_set.bandStretchAlpha << " s/m, twist scale " << m_set.bandTwistScale << ", twist loss " << (m_set.bandTwistLoss ? 1 : 0) << ", mu " << m_set.mu << ", tan delta " << m_set.tanDelta << " (geometry takes the next play)";
   }
}

static constexpr float kFlipperHorizonVpt = 1.0f;   // trace only: VPX time units (1 VPT = 10 ms) a bat is looked ahead for in the flip_clear numbers

void HitRubber::Collide(const CollisionEvent &coll)
{
   HitBall *const pball = coll.m_ball;
   const Vertex3Ds &n = coll.m_hitnormal;
   const float dot = -(n.Dot(pball->m_d.m_vel));
   if (dot < C_LOWNORMVEL)
      return;   // receding or sliding: nothing to solve (a resting ball is handled by Contact)

   ReloadIfChanged();
   if (ContactEngage(pball))
      g_pplayer->m_physics->OnBallWallHit(*pball, n, dot);   // the force-model contact is stepped by the engine from here on (no impulse here)
   else
   {
      // A post or ring never gets here (a band never gets to Collide: HitTest)
      if (ContactTraceOn())
         PLOGI << "StockImpulse " << m_name << " t=" << g_pplayer->m_time_msec << " held " << (pball->m_contact != nullptr ? 1 : 0) << " owns " << (ContactOwns(pball->m_d) ? 1 : 0) << " dot " << dot;
      pball->Collide3DWall(n, m_elasticity, m_elasticityFalloff, m_friction, m_scatter);
   }

   if (m_obj && m_fe && dot >= m_threshold)
   {
      m_obj->m_currentHitThreshold = dot;
      FireHitEvent(pball);
   }
}

HitObject *HitRubber::HitEventObject() const
{
   if (m_surface != nullptr)
      return m_eventCarrier;
   return m_set.band ? nullptr : const_cast<HitRubber *>(this);   // a rubber post / ring fires its own hit event (Collide above)
}

double HitRubber::LayerM() const
{
   if (m_set.layerMm > 0.0 || m_surface != nullptr)
      return ContactHost::LayerM();
   return 2.0 * m_radiusSI;
}

void HitRubber::Refresh(SurfaceContact &it, const double ballR, const double reachM, const RubberContact::V3 &c) const
{
   if (m_surface == nullptr)
      return;
   if (ContactSettings::Get().multiContact && m_linkOver > 0 && !m_linkNoticed)   // multi_contact: the load's count of edges on more than two faces, reported once
   {
      m_linkNoticed = true;
      const std::string msg = "Rubber surface " + m_name + ": " + std::to_string(m_linkOver) + " edges are on more than two faces; each joins its first two faces only (multi_contact)";
      PLOGW << msg;
      if (g_pplayer != nullptr && g_pplayer->m_liveUI != nullptr)
         g_pplayer->m_liveUI->PushNotification(msg, 5000);
   }
   if (it.qFaces != nullptr)
      return;   // shared_finder P2b (key on): the ball's query holds this object's faces (FinderBind)
   // A list built with some margin stays good until the ball has moved farther than the margin allows: scanning every face of a playfield mesh for every ball on every
   // call was the main cost of a table with several balls on it.
   const double need = ballR + reachM;
   if (!(it.candR >= 0.0 && (c - it.candC).Length() + need <= it.candR))
   {
      constexpr double kMargin = 0.008;   // m
      RubberContact::SurfaceCandidates(m_surface->tris, c, need + kMargin, it.candV);
      it.candC = c;
      it.candR = need + kMargin;
   }
   // The faces the ball can touch in this cycle: those whose box is within `need` of the ball (a face farther than the ball's radius plus its reach cannot touch it, and the nearest-face
   // query only has to find faces that can). The probes of the cycle's substeps test this short list, not the whole Verlet list. A surface with carved holes keeps the whole list (the hole moves the nearest point).
   const RubberContact::CarvedHole *holes = nullptr;
   int nHoles = 0;
   if (m_surfaceKind == "playfield")
      Turntables::Holes(holes, nHoles);
   it.cand.clear();
   it.tightC = c;
   if (nHoles > 0)
   {
      it.cand = it.candV;
      return;
   }
   it.tightC = c;
   const double need2 = (need + SurfaceContact::kTightSkin) * (need + SurfaceContact::kTightSkin);
   for (const int i : it.candV)
   {
      const RubberContact::SurfaceTri &t = m_surface->tris[i];
      const double dx = std::max(std::max(t.lo.x - c.x, c.x - t.hi.x), 0.0), dy = std::max(std::max(t.lo.y - c.y, c.y - t.hi.y), 0.0), dz = std::max(std::max(t.lo.z - c.z, c.z - t.hi.z), 0.0);
      if (dx * dx + dy * dy + dz * dz <= need2)
         it.cand.push_back(i);
   }
}

RubberContact::TubeProbe HitRubber::Probe(const RubberContact::BallSI &bs, const SurfaceContact &it, const RubberContact::V3 &c) const
{
   if (m_surface != nullptr && !g_nearCmp && it.qFaces != nullptr && FinderCheckOn())   // P3b check mode: the near pass's answer against the whole query list's (diagnostics only)
   {
      g_nearCmp = true;
      g_nearOff = true;
      const RubberContact::TubeProbe r0 = Probe(bs, it, c);
      g_nearOff = false;
      const RubberContact::TubeProbe r1 = Probe(bs, it, c);
      g_nearCmp = false;
      if (!SameProbe(r0, r1))
      {
         g_tab.probeDiff++;
         if (g_tab.loggedNear++ < 50)
            PLOGW << "Finder near pass probe differs: object " << Name() << ", seg " << r1.seg << " vs " << r0.seg << ", dist " << r1.dist << " vs " << r0.dist << ", cycle " << g_finder.cycle;
      }
      return r1;
   }
   CheckProfScope checkProf(0);   // P3a profile (diagnostics only)
   const bool finderCheck = m_surface != nullptr && g_finderDepth == 0 && FinderCheckOn();   // check mode (VPX_FINDER_CHECK=1): diagnostics only
   if (m_surface != nullptr)
      FinderInside(it, c);   // shared_finder P2a: the inside check of the ball's query volume (re-query counted; nothing reads the query yet)
   if (finderCheck)
      FinderCheck(it, c, bs.R, 0);
   const FinderGuard finderGuard(finderCheck);
   RubberContact::FaceTable *const ft = FinderTable(it, c);   // P3a (key on): the ball's measurement table of this object at c (after the inside check, which may re-query)
   const RubberContact::FaceRead frd { ft, nullptr };
   const RubberContact::FaceRead *const fr = ft != nullptr ? &frd : nullptr;
   checkProf.Tab(ft);
   const TabCheck tabCheck(ft, this, it, c);
   if (m_surface != nullptr)
   {
      const std::vector<int> *inL, *prL;   // P3b: the near pass's faces (the probe: the whole list with the box bound when none is within R)
      RubberContact::FaceRead frIn, frPr;
      NearLists(ft, it.Faces(c), true, inL, frIn, prL, frPr);
      const RubberContact::FaceRead *const pPr = ft != nullptr ? &frPr : nullptr;
      if (m_surfaceKind == "playfield")
      {
         const RubberContact::CarvedHole *holes = nullptr;
         int nHoles = 0;
         Turntables::Holes(holes, nHoles);   // the playfield has a round hole under every turntable disc
         if (PlayfieldCuts::Any(this))   // holes cut in the playfield (playfieldcut.h): the faces inside a cut left out, its exact rim's contacts added. The cuts within the face list's reach are tested
            return PlayfieldCuts::Probe(this, bs, m_surface->tris, it.Faces(c), m_set.edgeRadiusMm * 1e-3, c, it.candR >= 0.0 ? it.candR : 1e30, holes, nHoles, false, 0.0, ft, FinderCutsOf(it));
         if (nHoles > 0)   // P3b: a turntable's hole moves a face's nearest point to its rim, so the near pass's distances are not the probe's: the whole list as P3a
            return RubberContact::ProbeTriangles(bs, m_surface->tris, it.Faces(c), m_set.edgeRadiusMm * 1e-3, c, holes, nHoles, 0.0, nullptr, fr);
         return RubberContact::ProbeTriangles(bs, m_surface->tris, *prL, m_set.edgeRadiusMm * 1e-3, c, holes, nHoles, 0.0, nullptr, pPr);
      }
      return RubberContact::ProbeTriangles(bs, m_surface->tris, *prL, m_set.edgeRadiusMm * 1e-3, c, nullptr, 0, 0.0, nullptr, pPr);   // (no carved hole: a ball inside a non-legacy kicker's circle is handed to the stock contact instead, see InPlayfieldHole; a carved hole let balls fall through a plane playfield)
   }
   return RubberContact::ProbeTube(bs, m_ptsSI.data(), (int)m_ptsSI.size(), m_radiusSI, c, &m_tubePre);
}

RubberContact::TubeProbe HitRubber::ProbePatch(const RubberContact::BallSI &bs, const SurfaceContact &it, const RubberContact::V3 &c) const
{
   if (m_surface != nullptr && !g_nearCmp && it.qFaces != nullptr && FinderCheckOn())   // P3b check mode: as Probe
   {
      g_nearCmp = true;
      g_nearOff = true;
      const RubberContact::TubeProbe r0 = ProbePatch(bs, it, c);
      g_nearOff = false;
      const RubberContact::TubeProbe r1 = ProbePatch(bs, it, c);
      g_nearCmp = false;
      if (!SameProbe(r0, r1))
      {
         g_tab.probeDiff++;
         if (g_tab.loggedNear++ < 50)
            PLOGW << "Finder near pass patch differs: object " << Name() << ", seg " << r1.seg << " vs " << r0.seg << ", delta " << r1.delta << " vs " << r0.delta << ", cycle " << g_finder.cycle;
      }
      return r1;
   }
   CheckProfScope checkProf(1);   // P3a profile (diagnostics only)
   const bool finderCheck = m_surface != nullptr && g_finderDepth == 0 && FinderCheckOn();   // check mode (as Probe; the Probe inside is not checked again)
   if (m_surface != nullptr)
      FinderInside(it, c);   // shared_finder P2a: the inside check of the ball's query volume (re-query counted; nothing reads the query yet)
   if (finderCheck)
      FinderCheck(it, c, bs.R, 1);
   const FinderGuard finderGuard(finderCheck);
   RubberContact::FaceTable *const ft = FinderTable(it, c);   // P3a (key on): the ball's measurement table of this object at c (after the inside check, which may re-query)
   const RubberContact::FaceRead frd { ft, nullptr };
   const RubberContact::FaceRead *const fr = ft != nullptr ? &frd : nullptr;
   checkProf.Tab(ft);
   const TabCheck tabCheck(ft, this, it, c);
   const std::vector<int> *inL = nullptr, *prL = nullptr;   // P3b: the near pass's faces (as Probe); the patch reads only faces within R (kind 0 and plane distance below R - zero)
   RubberContact::FaceRead frIn, frPr;
   if (m_surface != nullptr)
      NearLists(ft, it.Faces(c), true, inL, frIn, prL, frPr);
   const RubberContact::FaceRead *const pIn = ft != nullptr ? &frIn : nullptr, *const pPr = ft != nullptr ? &frPr : nullptr;
   if (m_surface != nullptr && m_surfaceKind == "playfield" && PlayfieldCuts::Any(this))   // (as Probe, with the patch sum of the faces' folds or of the cut's rim contacts)
   {
      const RubberContact::CarvedHole *holes = nullptr;
      int nHoles = 0;
      Turntables::Holes(holes, nHoles);
      return PlayfieldCuts::Probe(this, bs, m_surface->tris, it.Faces(c), m_set.edgeRadiusMm * 1e-3, c, it.candR >= 0.0 ? it.candR : 1e30, holes, nHoles, true, it.gs.zero, ft, FinderCutsOf(it));
   }
   if (ContactSettings::Get().multiContact && !(m_surface != nullptr && m_surfaceKind == "playfield"))   // every pressed part its own Hertz bar (rubbercontact.h BarsTube / BarsTriangles); the playfield (turntable holes) keeps the patch sum
   {
      if (m_surface != nullptr)
         return RubberContact::BarsTriangles(bs, m_surface->tris, m_surface->vid, m_surface->nb, *prL, m_set.edgeRadiusMm * 1e-3, c, it.gs.zero, pPr, &m_surface->fold);   // P3b: the near pass's faces (bars only use faces within R)
      return RubberContact::BarsTube(bs, m_ptsSI.data(), (int)m_ptsSI.size(), m_radiusSI, c, m_tubePre, it.gs.zero);
   }
   RubberContact::TubeProbe pr = Probe(bs, it, c);
   if (m_surface != nullptr)
   {
      const RubberContact::CarvedHole *holes = nullptr;
      int nHoles = 0;
      if (m_surfaceKind == "playfield")
         Turntables::Holes(holes, nHoles);   // (as Probe)
      if (nHoles > 0)   // P3b: a turntable's hole: the whole list as P3a (as Probe)
         RubberContact::PatchTriangles(bs, m_surface->tris, it.Faces(c), c, it.gs.zero, pr, holes, nHoles, fr);
      else
         RubberContact::PatchTriangles(bs, m_surface->tris, *inL, c, it.gs.zero, pr, holes, nHoles, pIn);
   }
   else
      RubberContact::PatchTube(bs, m_ptsSI.data(), (int)m_ptsSI.size(), m_radiusSI, c, m_tubePre, it.gs.zero, pr);
   return pr;
}

double HitRubber::ArrivalS(const RubberContact::BallSI &bs, const SurfaceContact &it, const int nearestSeg, const RubberContact::V3 &c, const RubberContact::V3 &v, const double a, const double hMax, const double touchM, const double arriveM) const
{
   CheckProfScope checkProf(2);   // P3a profile (diagnostics only)
   const bool finderCheck = m_surface != nullptr && g_finderDepth == 0 && FinderCheckOn();   // check mode (as Probe)
   if (m_surface != nullptr)
      FinderInside(it, c);   // shared_finder P2a: the inside check of the ball's query volume (re-query counted; nothing reads the query yet)
   if (finderCheck)
      FinderCheck(it, c, bs.R, 2);
   const FinderGuard finderGuard(finderCheck);
   RubberContact::FaceTable *const ft = FinderTable(it, c);   // P3a (key on): the ball's measurement table of this object at c (after the inside check, which may re-query)
   const RubberContact::FaceRead frd { ft, nullptr };
   const RubberContact::FaceRead *const fr = ft != nullptr ? &frd : nullptr;
   checkProf.Tab(ft);
   const TabCheck tabCheck(ft, this, it, c);
   if (m_surface != nullptr && m_surfaceKind == "playfield" && PlayfieldCuts::Any(this))   // (as Probe: the faces not inside a cut, and the cut's rim contacts other than the probe's nearest)
      return PlayfieldCuts::ArrivalS(this, bs, m_surface->tris, it.Faces(c), nearestSeg, c, v, a, hMax, touchM, arriveM, it.candR >= 0.0 ? it.candR : 1e30, ft, FinderCutsOf(it));
   if (m_surface != nullptr)
      return RubberContact::ArrivalTriangles(bs, m_surface->tris, it.Faces(c), nearestSeg, c, v, a, hMax, touchM, arriveM, fr);
   return RubberContact::ArrivalTube(bs, m_ptsSI.data(), (int)m_ptsSI.size(), m_radiusSI, m_tubePre, nearestSeg, c, v, a, hMax, touchM, arriveM);
}

// ---- Engine-stepped contact ----
//
// A ball hitting a static rubber (a post ring, a wall rubber) used to be solved to the end inside the one Collide() call: 900 - 1500 substeps of 5 us
// (400 - 600 us, half a physics cycle, for every hit; two balls in the same cycle would not fit). Now the ball is engaged: every physics cycle RingAdvance()
// advances the same force model by that cycle's substeps (about 80 us), and the ball is handed back to the engine when it is clear of every rubber and moving away.
// Only the rubber's own contact force and torque act in the substeps; gravity, the playfield and every other object are left to the engine (as for a band).
// Every force-model rubber within reach of the ball is part of the contact, each with its own parameters and shear state, so a ball between two rubbers feels both.

bool HitRubberContactHolds(const HitBall *ball);

std::vector<BallPairState> g_ballPairs;

// A ball the pair model can take part in: not a flipper's, not in a kicker, not held by a band; a contact-held ball or a free one.
bool PairEligible(const HitBall *x)
{
   // Every ball takes part in pair contact except one held by a kicker (a state of its own, until the kicker is redesigned). Unified: the bat and the band are stepped
   // with the island, so their balls are members like any other; the legacy engine steps them apart and keeps its balls out.
   return !x->m_d.m_lockedInKicker && (ContactSystem::Unified() || x->m_bandOwner == nullptr);
}

bool HitRubber::ContactPairModel(const HitBall *a, const BallS &b)
{
   if (!ContactSettings::Get().ballContact || !PairEligible(a))
      return false;
   const HitBall *hb = nullptr;
   for (const auto *const ball : g_pplayer->m_vball)
      if (&ball->m_hitBall.m_d == &b)
         hb = &ball->m_hitBall;
   // Unified: the contact system owns every pair (two free balls too). The legacy engine only solves pairs with a held ball.
   return hb != nullptr && PairEligible(hb) && (ContactSystem::Unified() || HitRubberContactHolds(a) || HitRubberContactHolds(hb));
}

bool HitRubberContactPairModel(const HitBall *a, const BallS &b)
{
   return HitRubber::ContactPairModel(a, b);
}

// ---- A band's mover: the contact system steps the cord with its balls; once per cycle the cord steps alone if no ball did, and finishes up ----

bool BandMover::AddToList() const
{
   return m_owner->m_set.band;
}

void BandMover::PreUpdateDisplacements(const float dtime)
{
   m_owner->BandAdvance(dtime);
}

bool HitRubberContactHolds(const HitBall *ball)
{
   return ContactSystem::Holds(ball);
}

void HitRubberContactRelease(HitBall *ball, const char *why)
{
   ContactSystem::Release(ball, why);
}

void HitRubberBandBallDeleted(HitBall *ball)
{
   if (ball->m_bandOwner != nullptr)
      ball->m_bandOwner->BandRelease(ball, "ball deleted");
   ContactSystem::Release(ball, "ball deleted");
}

bool HitRubber::BandWants(const HitBall *ball) const
{
   return m_set.band && ball->m_bandOwner == nullptr;   // a contact (the playfield's, above all) is handed over to the band. No cooling down and no giving way to a flipper: there is no stock collision to leave the ball to
}

bool HitRubber::BandEngage(HitBall *const pb)
{
   if (!BandWants(pb))
      return false;
   CordEnsure();
   // the group step takes it from here: the cord finds the ball itself
   if (!ContactSystem::EngageBand(this, pb))
      return false;
   BandBall e;
   e.ball = pb;
   e.vin = RubberContact::V3((double)pb->m_d.m_vel.x * VS, (double)pb->m_d.m_vel.y * VS, (double)pb->m_d.m_vel.z * VS);
   m_engaged.push_back(e);
   pb->m_bandOwner = this;
   return true;
}

void HitRubber::BandRelease(HitBall *const ball, const char *why)
{
   size_t i = 0;
   while (i < m_engaged.size() && m_engaged[i].ball != ball)
      i++;
   if (i == m_engaged.size())
      return;
   const BandBall e = m_engaged[i];
   const Vertex3Ds &vel = ball->m_d.m_vel;
   const double speedOut = std::hypot(std::hypot(vel.x, vel.y), vel.z) * VS;
   PLOGI << "Band " << m_name << ": released ball (" << why << ") at t=" << g_pplayer->m_time_msec << " ms, speed in " << e.vin.Length() << " m/s, out " << speedOut << " m/s, vel out ("
         << vel.x * VS << ", " << vel.y * VS << ", " << vel.z * VS << ") m/s, balls left on the band " << m_engaged.size() - 1;
   if (ball->m_contact != nullptr && ball->m_contact->band == this)
      ball->m_contact->band = nullptr;   // (the ball may stay in its other contacts)
   ball->m_bandOwner = nullptr;
   m_engaged.erase(m_engaged.begin() + (long)i);
   m_info.hitTimeMs = (int)g_pplayer->m_time_msec;
   m_info.speedIn = (float)e.vin.Length();
   m_info.speedOut = (float)speedOut;
   if (std::strcmp(why, "ball deleted") != 0 && std::strcmp(why, "band switched off") != 0 && m_editable != nullptr)
      static_cast<Rubber *>(m_editable)->FireVoidEventParm(DISPID_RubberEvents_BandHit, m_info.peakMm);   // the script handler may do anything: nothing of ours is in flight
}

void HitRubber::BandReleaseAll(const char *why)
{
   while (!m_engaged.empty())
      BandRelease(m_engaged.back().ball, why);
}

void HitRubber::BandAdvance(const float dtime)
{
   CordEnsure();
   if (!CordOn())
      return;
   if (!m_cordStepped)
      ContactSystem::StepCordAlone(this);   // (for the seconds the contact system's Step stored)
   m_cordStepped = false;
   CordFinish(dtime);
}

// ---- The cord in the group step (docs/unified-physics/GroupStepDesign.txt 8, 11) ----
void HitRubber::CordEnsure()
{
   if (m_cord != nullptr || !m_set.band || m_bandPosts.size() < 2)
      return;
   // the fixed-height cord (cordfixed.h): no copy of the rubber's values
   {
      using V3 = GroupStep::V3;
      auto cf = std::make_shared<CordBand>();
      cf->fixed = std::make_unique<CordFixed::Part>();
      CordFixed::Part &p = *cf->fixed;
      const double rc = m_set.bandSectionMm > 0.0 ? 0.5e-3 * m_set.bandSectionMm : m_radiusSI, rcDelta = rc - m_radiusSI;
      p.rc = rc;
      p.EA = m_set.eMpa * 1e6 * M_PI * rc * rc;
      p.rc = rc;
      FixedCordValues(p, m_set);
      double zc = 0.0;
      for (const Vertex3Ds &v : m_pts)
         zc += v.z;
      zc = m_pts.empty() ? 0.0 : zc / (double)m_pts.size() * U;
      p.zc = zc;
      // the posts, in the loop's order, anticlockwise: a span each, end points in the notches at the cord's height
      const int np = (int)m_bandPosts.size();
      double area = 0.0;
      for (int i = 0; i < np; i++)
         area += m_bandPosts[i].x * m_bandPosts[(i + 1) % np].y - m_bandPosts[(i + 1) % np].x * m_bandPosts[i].y;
      p.circles.clear();
      for (int i = 0; i < np; i++)
      {
         const BandPost &bp = m_bandPosts[area < 0.0 ? np - 1 - i : i];
         CordFixed::Part::Circle q;
         q.x = bp.x * U;
         q.y = bp.y * U;
         q.r = bp.r * U + rcDelta;
         for (const auto &e : p.circles)
            q.twin = q.twin || (std::fabs(e.x - q.x) < 1e-6 && std::fabs(e.y - q.y) < 1e-6);   // a post the loop passes twice is one post
         p.circles.push_back(q);
      }
      // the stretch at rest: the ring on the loop (tangents and arcs), or the percentage
      struct Tan { V3 a, b; };
      std::vector<Tan> tans(np);
      double loop = 0.0, prevAng = 0.0, firstAng = 0.0;
      std::vector<double> ang(np);
      for (int i = 0; i < np; i++)
      {
         const auto &a = p.circles[i], &b = p.circles[(i + 1) % np];
         const double dx = b.x - a.x, dy = b.y - a.y, L = std::hypot(dx, dy), sn = (a.r - b.r) / L, cs = std::sqrt(std::max(1.0 - sn * sn, 0.0)), tx = dx / L, ty = dy / L;
         const double nx = ty * cs + tx * sn, ny = -tx * cs + ty * sn;   // the outside is on the right of the travel
         tans[i].a = V3(a.x + nx * a.r, a.y + ny * a.r, zc);
         tans[i].b = V3(b.x + nx * b.r, b.y + ny * b.r, zc);
         loop += (tans[i].b - tans[i].a).Length();
         ang[i] = std::atan2(tans[i].b.y - tans[i].a.y, tans[i].b.x - tans[i].a.x);
      }
      (void)prevAng;
      (void)firstAng;
      for (int i = 0; i < np; i++)
      {
         double turn = ang[i] - ang[(i + np - 1) % np];
         while (turn < 0.0)
            turn += 2.0 * M_PI;
         while (turn >= 2.0 * M_PI)
            turn -= 2.0 * M_PI;
         loop += p.circles[i].r * turn;
      }
      const double ringIn = CordRing(loop, rc, p.EA);
      double stretch = loop / (M_PI * (ringIn * 25.4e-3 + 2.0 * rc)) - 1.0;
      if (stretch <= 0.0)
      {
         PLOGW << "Band " << m_name << ": a " << ringIn << " in ring is longer than the loop on its posts (" << loop * 1e3 << " mm): no stretch";
         stretch = 0.001;
      }
      p.T0 = p.EA * stretch;
      for (int i = 0; i < np; i++)
      {
         CordFixed::Span sp;
         sp.Set(tans[i].a, tans[i].b, stretch, p.EA);
         p.span.push_back(sp);
      }
      const auto LawOf = [&](GroupStep::Item &it, const double reff) {
         const RubberContact::RubberParams rp = p.ringRp;
         const double Estar = rp.eMpa * 1e6 / (1.0 - rp.nu * rp.nu);
         it.reff = reff;
         it.kH = 4.0 / 3.0 * Estar * std::sqrt(reff);
         it.hL = 2.0 * rc;
         it.cL = it.hL > 0.0 ? it.kH * rp.betaLayer * reff * std::sqrt(reff) : 0.0;
      };
      // the posts behind the band (band_stops, as below): each goes on the span it stands behind, by its nearest line; an Item for the ball
      std::vector<std::shared_ptr<CordFixed::StopBody>> stopBodies;
      {
         const std::string &rest = m_set.bandStops;
         size_t pos = 0;
         while (pos < rest.size())
         {
            size_t bar = rest.find('|', pos);
            if (bar == std::string::npos)
               bar = rest.size();
            const std::string item = rest.substr(pos, bar - pos);
            double a = 0, b = 0, cc = 0, x = 0, y = 0, r = 0;
            bool have = false;
            if (sscanf(item.c_str(), " %lf , %lf , %lf", &a, &b, &cc) == 3 && cc > 0.0)
            {
               x = a;
               y = b;
               r = cc;
               have = true;
            }
            else if (!Trim(item).empty() && m_stopResolver)
            {
               std::string nm = Trim(item);
               bool hasR = false;
               double rIn = 0.0;
               const size_t comma = nm.find(',');
               if (comma != std::string::npos)
               {
                  hasR = sscanf(nm.c_str() + comma + 1, " %lf", &rIn) == 1;
                  nm = Trim(nm.substr(0, comma));
               }
               have = m_stopResolver(nm, hasR, rIn, x, y, r) && r > 0.0;   // (false for a part with Collidable off: it stops nothing)
            }
            if (have)
            {
               CordFixed::Part::Circle q;
               q.x = x * U;
               q.y = y * U;
               q.r = r * U + rc;
               p.circles.push_back(q);
               auto body = std::make_shared<CordFixed::StopBody>();
               body->c0 = V3(q.x, q.y, zc);
               body->rCore = r * U;
               stopBodies.push_back(body);
               int best = -1;
               double bestD = 0.03, bestX = 0.0;
               for (int j = 0; j < np; j++)
               {
                  const CordFixed::Span &sp = p.span[j];
                  const V3 w0 = body->c0 - sp.A;
                  const double st = w0.Dot(sp.t), dl = std::fabs(w0.Dot(sp.yh)) - q.r;
                  if (st > 0.0 && st < sp.len && dl < bestD)
                  {
                     bestD = dl;
                     best = j;
                     bestX = st;
                  }
               }
               if (best >= 0)
               {
                  CordFixed::Obstacle ob;
                  ob.body = body.get();
                  ob.x = bestX;
                  LawOf(ob.law, rc / std::sqrt(2.0));
                  {
                     CordFixed::Span &sb = p.span[best];
                     double g0;
                     V3 mv0;
                     body->Gap(sb.A + sb.t * bestX, rc, g0, mv0);
                     if (g0 < 0.0)   // the post is in the cord's straight line: the cord lies bent round it at rest
                     {
                        ob.bend = (mv0.Dot(sb.yh) >= 0.0 ? -g0 : g0);
                        sb.bent = true;
                        sb.bendMost = std::max(sb.bendMost, -g0);
                     }
                  }
                  p.span[best].ob.push_back(ob);
               }
            }
            pos = bar + 1;
         }
      }
      p.stopBodies = stopBodies;
      // the plate, as the sling's own setup measured it (SlingPlace); seated on a hard stop (the part's rule, no soft stop)
      const int n = (int)m_pts.size() - 1;
      if (m_slingPlate != nullptr && m_plateGeomValid && m_slingSpan >= 0 && n >= 4)
      {
         const BandSpan &bs = m_bandSpans[m_slingSpan];
         const Vertex3Ds &A = m_pts[((bs.first % n) + n) % n], &B = m_pts[((bs.last % n) + n) % n];
         const V3 aw(A.x * U, A.y * U, zc), bw(B.x * U, B.y * U, zc);
         int best = 0;
         double bestE = 1e30;
         for (int j = 0; j < np; j++)
         {
            const double e1 = (p.span[j].A - aw).Length() + (p.span[j].B - bw).Length(), e2 = (p.span[j].A - bw).Length() + (p.span[j].B - aw).Length(), e = std::min(e1, e2);
            if (e < bestE)
            {
               bestE = e;
               best = j;
            }
         }
         CordFixed::Span &hs = p.span[best];
         const V3 e0 = aw + (bw - aw) * m_plateF0, e1 = aw + (bw - aw) * m_plateF1, nv(m_slingNx, m_slingNy, 0.0);
         p.pl.on = true;
         p.pl.span = best;
         p.pl.x0 = std::min((e0 - hs.A).Dot(hs.t), (e1 - hs.A).Dot(hs.t));
         p.pl.x1 = std::max((e0 - hs.A).Dot(hs.t), (e1 - hs.A).Dot(hs.t));
         p.pl.sgn = hs.yh.Dot(nv) >= 0.0 ? 1.0 : -1.0;
         p.pl.originN = m_plateOriginN - rcDelta;   // measured from the table's own centreline; the span is the cord's, moved out by the real section with the plate where it is
         p.pl.a = m_plateA;
         p.pl.h = m_plateH;
         p.pl.rest = m_plateRest;
         LawOf(p.pl.law, rc / std::sqrt(2.0));
         for (const double xs : { p.pl.x0, p.pl.x1 })
         {
            CordFixed::Obstacle ob;
            ob.body = nullptr;
            ob.x = xs;
            LawOf(ob.law, rc / std::sqrt(2.0));
            hs.ob.push_back(ob);
         }
         p.sl = m_set.sling;
         p.SetupPlate();
         p.driver.SetDelay(m_set.sling.coilDelayMs * 1e-3);
         p.swGapM = m_set.sling.switchGapMm * 1e-3;
         p.geo.hasPlate = true;
      }
      p.switchGapM = m_set.sling.switchGapMm * 1e-3;
      for (const BandSwitch &sw : m_switches)
         if (sw.span >= 0 && sw.span < (int)m_bandSpans.size() && p.geo.nSwitches < CordFixed::Part::kMaxSwitches && n >= 4)
         {
            const BandSpan &sp = m_bandSpans[sw.span];
            const Vertex3Ds &A = m_pts[((sp.first % n) + n) % n], &B = m_pts[((sp.last % n) + n) % n];
            CordFixed::Part::SwitchSpec &ss = p.geo.sw[p.geo.nSwitches++];
            ss.x = (A.x + (B.x - A.x) * sw.f) * U;
            ss.y = (A.y + (B.y - A.y) * sw.f) * U;
            ss.bx = sw.bx;
            ss.by = sw.by;
            ss.coil = sw.coil;
         }
      p.SetNShapes();
      p.Setup(0);
      p.SetupRest();
      PLOGI << "Band " << m_name << ": values at load: stretch alpha " << p.stretchLaw.alpha << " s/m, twist scale " << p.twistScale << ", twist loss " << (p.twistLoss ? 1 : 0) << ", mu " << p.mu << ", tan delta " << p.tanD << ", roll smooth " << ContactSettings::Get().rollSmoothRadS << " rad/s, section " << 2e3 * rc << " mm" << (m_set.bandSectionMm > 0.0 ? " (band_section_mm)" : " (the rubber's thickness)") << ", ring " << ringIn << " in" << (m_set.bandRingIn > 0.0 ? " (band_ring_in)" : " (automatic)") << ", tension " << p.T0 << " N";
      for (int j = 0; j < np; j++)
         if (p.span[j].bent)
            PLOGI << "Band " << m_name << ": a stop holds span " << j << " " << p.span[j].bendMost * 1e3 << " mm off its straight line at rest: the cord's centre lifts " << p.span[j].restLift * 1e3 << " mm, the rest is squash; tension " << p.span[j].restT << " N, " << p.T0 << " N straight";
      m_cord = cf;
      m_cordFired = 0;
      for (int j = 0; j < np; j++)
         PLOGI << "Band " << m_name << ": span " << j << " from (" << p.span[j].A.x / U << ", " << p.span[j].A.y / U << ") to (" << p.span[j].B.x / U << ", " << p.span[j].B.y / U << ") U";
      PLOGI << "Band " << m_name << ": fixed-height cord, " << np << " posts, " << stopBodies.size() << " stops, " << (p.pl.on ? "plate, " : "") << p.geo.nSwitches << " switches, cord " << 2e3 * rc << " mm at " << zc * 1e3
            << " mm, E " << m_set.eMpa << " MPa, EA " << p.EA << " N, tan delta " << p.tanD << ", loop " << loop * 1e3 << " mm, stretch " << stretch * 100.0 << " %, tension " << p.T0 << " N"
            << (p.pl.on ? ", plate edges " : "") << (p.pl.on ? std::to_string(p.pl.x0 * 1e3) + " - " + std::to_string(p.pl.x1 * 1e3) + " mm of " + std::to_string(p.span[p.pl.span].len * 1e3) + " mm, originN " + std::to_string(p.pl.originN * 1e3) + " a " + std::to_string(p.pl.a * 1e3) + " h " + std::to_string(p.pl.h * 1e3) + " rest " + std::to_string(p.pl.rest * 57.29578) + " deg, " + std::to_string(p.sl.mech.vSupply) + " V, stroke at the cord " + std::to_string(p.PlateReachM() * 1e3) + " mm (band_plate_reach " + (m_set.bandPlateReach ? "1" : "0") + ")" : "");
   }
}

// The ring on this loop: band_ring_in, or the standard size (3/4 in, then 1 in and up in 1/4 in steps) whose tension on the loop is nearest the middle of 35 - 60 N, among the sizes within that range.
// None within it: the size whose tension is nearest to the range, and the log says so. Tension = EA (loop / ring's centreline - 1), the cord's own law; the ring's centreline is pi (inner diameter + section).
double HitRubber::CordRing(const double loopM, const double rcM, const double EA) const
{
   if (m_set.bandRingIn > 0.0)
      return m_set.bandRingIn;
   constexpr double kLow = 35.0, kHigh = 60.0, kMid = 0.5 * (kLow + kHigh), kLargest = 20.0;   // N, N, N, in (the largest band_ring_in takes)
   const auto TensionOf = [&](const double in) { return EA * (loopM / (M_PI * (in * 25.4e-3 + 2.0 * rcM)) - 1.0); };
   double best = 0.0, bestOff = 1e30, nearest = 0.0, nearestOff = 1e30;
   for (double in = 0.75; in <= kLargest; in = in < 1.0 ? 1.0 : in + 0.25)
   {
      const double T = TensionOf(in), off = T < kLow ? kLow - T : T > kHigh ? T - kHigh : 0.0;
      if (off == 0.0 && std::fabs(T - kMid) < bestOff)
      {
         best = in;
         bestOff = std::fabs(T - kMid);
      }
      if (off < nearestOff)
      {
         nearest = in;
         nearestOff = off;
      }
   }
   if (best > 0.0)
      return best;
   PLOGW << "Band " << m_name << ": no standard ring gives a tension within " << kLow << " - " << kHigh << " N on this loop (" << loopM * 1e3 << " mm, section " << 2e3 * rcM << " mm): taking the nearest, " << nearest << " in, " << TensionOf(nearest) << " N";
   return nearest;
}

void HitRubber::CordBegin()
{
   CordFixed::Part &p = *m_cord->fixed;
   if (!m_slingEnabled)
      p.driver.Cancel();
   if (p.geo.hasPlate)
   {
      if (m_slingEnabled && m_slingTriggerReq)
         p.driver.Trigger();
      m_slingTriggerReq = false;
   }
   if (m_slingEnabled)
      p.CycleBegin();
}

void HitRubber::CordFinish(const float /*dtime*/)
{
   for (const auto &h : m_cord->hits)   // the script's Hit event, as Collide fires it for a stock hit
      if (m_obj && m_fe && h.second >= m_threshold)
      {
         m_obj->m_currentHitThreshold = h.second;
         FireHitEvent(h.first);
      }
   m_cord->hits.clear();
   CordFixed::Part &p = *m_cord->fixed;
   std::vector<std::pair<int, bool>> changes;
   const int ns = std::min((int)m_switches.size(), p.geo.nSwitches);
   bool anyClosed = false;
   for (int i = 0; i < ns; i++)
   {
      BandSwitch &sw = m_switches[i];
      sw.dispM = m_cord->attached ? p.SwitchPush(i) : 0.0;
      const bool closed = p.SwitchIsClosed(i);
      if (closed != sw.closed)
         changes.emplace_back(i, closed);
      sw.closed = closed;
      anyClosed = anyClosed || (closed && sw.coil);
   }
   if (m_editable != nullptr)
      for (const auto &c : changes)
         static_cast<Rubber *>(m_editable)->FireVoidEventParm(c.second ? DISPID_RubberEvents_SwitchClosed : DISPID_RubberEvents_SwitchOpened, c.first);
   if (m_slingPlate != nullptr && p.geo.hasPlate)
   {
      m_plate.st.phi = p.phi;
      m_plate.st.omega = p.omega;
      m_slingAnyClosed = anyClosed;
      if (p.driver.fired != m_cordFired)
      {
         m_cordFired = p.driver.fired;
         m_slingInfo.fireTimeMs = (int)g_pplayer->m_time_msec;
         static_cast<Rubber *>(m_editable)->FireGroupEvent(DISPID_RubberEvents_SlingFire);
      }
      SlingFinish();
   }
   else if (!m_switches.empty())
      SwitchesFinish();
}

void HitRubber::CordSolveFailed(const char *why)
{
   PLOGW << "Band " << m_name << ": " << why;
   if (g_pplayer != nullptr && g_pplayer->m_liveUI != nullptr)
      g_pplayer->m_liveUI->PushNotification("Band " + m_name + ": " + why, 5000);
   m_cord->attached = false;
   m_cord->slots.clear();
   BandReleaseAll("solve failed");
}

// The mesh animation: where the rendered tube goes when the band is deflected: the cord's own path, at its own height (the ring of the mesh centred at c shifts and turns).
bool HitRubber::BandFixedActive() const { return m_cord != nullptr && m_cord->fixed->DrawActive(); }
void HitRubber::BandFixedDrawPrepare() const
{
   if (m_cord != nullptr)
      m_cord->fixed->DrawPrepare();
}
bool HitRubber::BandDeformAt(const Vertex3Ds &c, Vertex3Ds &shift, float rot[9]) const
{
   double sx, sy, ca, sa;
   if (m_cord == nullptr || !m_cord->fixed->DrawAt((double)c.x * U, (double)c.y * U, sx, sy, ca, sa))
      return false;
   shift = Vertex3Ds((float)(sx / U), (float)(sy / U), 0.f);
   rot[0] = (float)ca; rot[1] = (float)-sa; rot[2] = 0.f;
   rot[3] = (float)sa; rot[4] = (float)ca; rot[5] = 0.f;
   rot[6] = 0.f; rot[7] = 0.f; rot[8] = 1.f;
   return true;
}

void HitRubber::Contact(CollisionEvent &coll, const float dtime)
{
   if (SlowContact(coll.m_ball))
      return;   // a post rubber: the force model carries the resting or slow touch
   HitObject::Contact(coll, dtime);
}

void HitRubber::DrawUI(std::function<Vertex2D(Vertex3Ds)> project, ImDrawList *drawList, bool fill) const
{
   if (!m_enabled)
      return;
   const ImU32 lCol = ImGui::GetColorU32(ImGuiCol_PlotLines);
   for (size_t i = 0; i + 1 < m_pts.size(); i++)
   {
      const Vertex2D p0 = project(m_pts[i]);
      const Vertex2D p1 = project(m_pts[i + 1]);
      if (p0.x != FLT_MAX && p1.x != FLT_MAX)
         drawList->AddLine(ImVec2(p0.x, p0.y), ImVec2(p1.x, p1.y), lCol);
   }
}
void HitRubber::FinderBuild(const bool atLoad)
{
   FinderIndex &g = g_finder;
   if (atLoad && g.built)
   {
      FinderReport();
      g.Reset();
   }
   if (g.built)
      return;
   g.built = true;
   std::vector<HitRubber *> hosts;
   int skipped = 0;
   for (ContactHost *const o : ContactHost::Registry())
   {
      HitRubber *const h = dynamic_cast<HitRubber *>(o);
      if (h == nullptr || h->m_surface == nullptr || h->m_surface->tris.empty())
         continue;
      if (h->m_surfaceKind != "playfield" && h->m_surfaceKind != "wall" && h->m_surfaceKind != "ramp" && h->m_surfaceKind != "primitive")
      {
         skipped++;
         continue;
      }
      hosts.push_back(h);
   }
   if (hosts.empty())
      return;
   if (!PhysicsIni::Active(hosts[0]->m_ini).sharedFinder)
   {
      PLOGI << "Shared finder: off ([Features] shared_finder = 0)";
      return;
   }
   g.on = true;
   if (!PhysicsIni::Active(hosts[0]->m_ini).proximityEngage)
   {
      // with shared_finder on, stock's copies of converted meshes are out of stock's trees (P9), and proximity_engage = 0 leaves the contact system nothing that engages a free ball
      // with them: balls pass through walls, ramps and the playfield (old-law switch, to-do 10.13)
      const std::string msg = "[Features] proximity_engage = 0 with shared_finder = 1: no engagement path for converted surfaces, balls will fall through (set proximity_engage = 1)";
      PLOGW << msg;
      if (g_pplayer != nullptr && g_pplayer->m_liveUI != nullptr)
         g_pplayer->m_liveUI->PushNotification(msg, 10000);
   }
   {
      int i = 0;   // P2b: the registry's order, kept for a ball's surfaces (Regather and EngageImpl sort them by it)
      for (ContactHost *const o : ContactHost::Registry())
         o->m_regOrder = i++;
   }
   const auto t0 = std::chrono::steady_clock::now();
   double slowestMs = 0.0, l2Ms = 0.0;
   std::string slowest;
   size_t faces = 0, l2Nodes = 0, bytes = 0, slowestFaces = 0;
   int pieces = 0, disabled = 0;
   for (HitRubber *const h : hosts)
   {
      const auto a = std::chrono::steady_clock::now();
      RubberSurface &s = *h->m_surface;
      const int n = (int)s.tris.size();
      // mesh pieces: faces joined through shared edges (nb, built by Link at load), numbered in load order of each piece's first face
      s.piece.assign(n, -1);
      s.nPieces = 0;
      std::vector<int> stack;
      for (int i = 0; i < n; i++)
      {
         if (s.piece[i] >= 0)
            continue;
         s.piece[i] = s.nPieces;
         stack.push_back(i);
         while (!stack.empty())
         {
            const int f = stack.back();
            stack.pop_back();
            if (s.nb.size() != s.tris.size())
               continue;
            for (int e = 0; e < 3; e++)
            {
               const int m = s.nb[f][e];
               if (m >= 0 && s.piece[m] < 0)
               {
                  s.piece[m] = s.nPieces;
                  stack.push_back(m);
               }
            }
         }
         s.nPieces++;
      }
      s.fnodes.clear();
      s.forder.resize(n);
      for (int i = 0; i < n; i++)
         s.forder[i] = i;
      BuildTree(s.fnodes, s.forder, [&s](const int i, RubberContact::V3 &lo, RubberContact::V3 &hi) { lo = s.tris[i].lo; hi = s.tris[i].hi; }, 0, n, 4);
      g.objs.push_back({ h, s.fnodes[0].lo, s.fnodes[0].hi });
      s.fobj = (int)g.objs.size() - 1;
      const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
      l2Ms += ms;
      if (ms > slowestMs)
      {
         slowestMs = ms;
         slowest = h->m_name;
         slowestFaces = n;
      }
      faces += n;
      pieces += s.nPieces;
      l2Nodes += s.fnodes.size();
      bytes += s.fnodes.size() * sizeof(FinderNode) + s.forder.size() * sizeof(int) + s.piece.size() * sizeof(int);
      if (!h->m_enabled)
         disabled++;
   }
   const auto b = std::chrono::steady_clock::now();
   g.order.resize(g.objs.size());
   for (size_t i = 0; i < g.objs.size(); i++)
      g.order[i] = (int)i;
   BuildTree(g.nodes, g.order, [&g](const int i, RubberContact::V3 &lo, RubberContact::V3 &hi) { lo = g.objs[i].lo; hi = g.objs[i].hi; }, 0, (int)g.objs.size(), 4);
   bytes += g.nodes.size() * sizeof(FinderNode) + g.order.size() * sizeof(int) + g.objs.size() * sizeof(FinderObj);
   const auto t1 = std::chrono::steady_clock::now();
   const double l1Ms = std::chrono::duration<double, std::milli>(t1 - b).count(), totalMs = std::chrono::duration<double, std::milli>(t1 - t0).count();
   PLOGI << "Shared finder: index built for " << g.objs.size() << " static mesh objects (" << disabled << " not collidable at load; " << skipped << " other surface hosts left out), " << faces << " faces in "
         << pieces << " mesh pieces; nodes: " << l2Nodes << " in the objects' face trees, " << g.nodes.size() << " in the object tree; memory " << bytes << " bytes (" << sizeof(FinderNode)
         << " per node); build " << totalMs << " ms (face trees and pieces " << l2Ms << " ms, slowest " << slowest << " " << slowestMs << " ms for " << slowestFaces << " faces; object tree " << l1Ms << " ms)"
         << (atLoad ? "; at load" : "; at the first cycle (the load call did not run)") << (FinderCheckOn() ? "; check mode on (VPX_FINDER_CHECK)" : "");
}

void HitRubber::FinderCheck(const SurfaceContact &it, const RubberContact::V3 &c, const double ballR, const int evalKind) const
{
   if (!g_finder.built)
      FinderBuild();   // (an evaluation before the first cycle's gather)
   FinderIndex &g = g_finder;
   if (!g.on || m_surface == nullptr)
      return;
   static const char *const kinds[3] = { "probe", "probe patch", "arrival" };
   const std::vector<RubberContact::SurfaceTri> &tris = m_surface->tris;
   const double r2 = ballR * ballR;
   thread_local std::vector<int> full, listA, listB, objHit;
   thread_local std::vector<double> dist;
   full.clear();
   dist.clear();
   for (int i = 0; i < (int)tris.size(); i++)   // the slow way: every face, an exact box reject then the exact distance
   {
      if (BoxD2(tris[i].lo, tris[i].hi, c) > r2)
         continue;
      int kind = 0;
      const RubberContact::V3 p = RubberContact::ClosestOnTriangle(tris[i], c, kind);
      const RubberContact::V3 dv = c - p;
      const double d2 = dv.Dot(dv);
      if (d2 <= r2)
      {
         full.push_back(i);
         dist.push_back(std::sqrt(d2));
      }
   }
   std::lock_guard<std::mutex> lk(g.mtx);
   FinderPeriodic();
   g.evals++;
   if (full.empty())
      return;
   g.facesNear += full.size();
   listA = it.Faces(c);
   std::sort(listA.begin(), listA.end());
   objHit.clear();
   QueryTree(g.nodes, g.order, [&g](const int i, RubberContact::V3 &lo, RubberContact::V3 &hi) { lo = g.objs[i].lo; hi = g.objs[i].hi; }, c, r2, objHit);
   const bool found = std::any_of(objHit.begin(), objHit.end(), [&](const int k) { return g.objs[k].host == this; });
   listB.clear();
   if (found)
   {
      const RubberSurface &s = *m_surface;
      QueryTree(s.fnodes, s.forder, [&s](const int i, RubberContact::V3 &lo, RubberContact::V3 &hi) { lo = s.tris[i].lo; hi = s.tris[i].hi; }, c, r2, listB);
      std::sort(listB.begin(), listB.end());
   }
   // P2a: the ball's per-cycle query (after the inside check) must hold every face within the ball's radius
   const FinderBallQ *const q = static_cast<const FinderBallQ *>(it.finderQ);
   int qa = 0, qb = 0;
   if (q == nullptr)
      g.noQ++;
   else
   {
      const auto o = std::lower_bound(q->objs.begin(), q->objs.end(), m_surface->fobj);
      if (o != q->objs.end() && *o == m_surface->fobj)
      {
         qa = q->first[o - q->objs.begin()];
         qb = q->first[o - q->objs.begin() + 1];
      }
   }
   for (size_t k = 0; k < full.size(); k++)
   {
      const int f = full[k];
      const int pc = (int)m_surface->piece.size() > f ? m_surface->piece[f] : -1;
      if (q != nullptr && !std::binary_search(q->faces.begin() + qa, q->faces.begin() + qb, f))
      {
         g.missQ++;
         FinderMiss(4, m_name, f, pc, dist[k], c, kinds[evalKind]);
      }
      if (!std::binary_search(listA.begin(), listA.end(), f))
      {
         g.missA++;
         FinderMiss(0, m_name, f, pc, dist[k], c, kinds[evalKind]);
      }
      if (!std::binary_search(listB.begin(), listB.end(), f))
      {
         g.missB++;
         FinderMiss(1, m_name, f, pc, dist[k], c, kinds[evalKind]);
      }
   }
}

void HitRubber::FinderCheckGather(const HitBall *ball, const std::vector<SurfaceContact> &items, const RubberContact::V3 &c, const double ballR, const bool inPlayfieldHole)
{
   if (!g_finder.built)
      FinderBuild();   // (an evaluation before the first cycle's gather)
   FinderIndex &g = g_finder;
   if (!g.on || !FinderCheckOn())
      return;
   const double r2 = ballR * ballR;
   thread_local std::vector<int> objHit;
   objHit.clear();
   QueryTree(g.nodes, g.order, [&g](const int i, RubberContact::V3 &lo, RubberContact::V3 &hi) { lo = g.objs[i].lo; hi = g.objs[i].hi; }, c, r2, objHit);
   std::lock_guard<std::mutex> lk(g.mtx);
   FinderPeriodic();
   g.gathers++;
   for (int k = 0; k < (int)g.objs.size(); k++)
   {
      const FinderObj &ob = g.objs[k];
      const HitRubber *const h = ob.host;
      if (!h->m_enabled || !h->m_set.contact || h->m_set.band || (inPlayfieldHole && h->m_surfaceKind == "playfield"))
         continue;   // (as the gather: a disabled object, one off the contact model, and the playfield over a kicker's hole are left out there on purpose)
      if (BoxD2(ob.lo, ob.hi, c) > r2)
         continue;
      const std::vector<RubberContact::SurfaceTri> &tris = h->m_surface->tris;
      int nearest = -1;
      double nearestD2 = r2;
      for (int i = 0; i < (int)tris.size(); i++)
      {
         if (BoxD2(tris[i].lo, tris[i].hi, c) > nearestD2)
            continue;
         int kind = 0;
         const RubberContact::V3 p = RubberContact::ClosestOnTriangle(tris[i], c, kind);
         const RubberContact::V3 dv = c - p;
         const double d2 = dv.Dot(dv);
         if (d2 <= nearestD2)
         {
            nearestD2 = d2;
            nearest = i;
         }
      }
      if (nearest < 0)
         continue;
      const int pc = (int)h->m_surface->piece.size() > nearest ? h->m_surface->piece[nearest] : -1;
      if (std::none_of(items.begin(), items.end(), [h](const SurfaceContact &it) { return it.r == static_cast<const ContactSurface *>(h); }))
      {
         g.gMissA++;
         FinderMiss(2, h->m_name, nearest, pc, std::sqrt(nearestD2), c, "gather");
      }
      if (std::find(objHit.begin(), objHit.end(), k) == objHit.end())
      {
         g.gMissB++;
         FinderMiss(3, h->m_name, nearest, pc, std::sqrt(nearestD2), c, "gather");
      }
   }
   // P2a: every mesh object today's gather engages must be among the per-cycle query's objects (one with no face within R + reach of the query's centre: the gather is wider)
   const FinderBallQ *const q = items.empty() ? nullptr : static_cast<const FinderBallQ *>(items[0].finderQ);
   if (q != nullptr)
      for (int k = 0; k < (int)g.objs.size(); k++)
      {
         const HitRubber *const h = g.objs[k].host;
         if (std::none_of(items.begin(), items.end(), [h](const SurfaceContact &it) { return it.r == static_cast<const ContactSurface *>(h); }) || std::binary_search(q->objs.begin(), q->objs.end(), k))
            continue;
         const std::vector<RubberContact::SurfaceTri> &tris = h->m_surface->tris;
         const double rq = q->R + q->ext;
         int nearest = -1;
         double nearestD2 = rq * rq;
         for (int i = 0; i < (int)tris.size(); i++)
         {
            if (BoxD2(tris[i].lo, tris[i].hi, q->c0) > nearestD2)
               continue;
            int kind = 0;
            const RubberContact::V3 dv = q->c0 - RubberContact::ClosestOnTriangle(tris[i], q->c0, kind);
            if (dv.Dot(dv) <= nearestD2)
            {
               nearestD2 = dv.Dot(dv);
               nearest = i;
            }
         }
         if (nearest < 0)
            g.gWider++;
         else
         {
            g.gMissQ++;
            FinderMiss(4, h->m_name, nearest, (int)h->m_surface->piece.size() > nearest ? h->m_surface->piece[nearest] : -1, std::sqrt(nearestD2), q->c0, "gather");
         }
      }
   (void)ball;
}

// ---- The shared finder, piece P2a: the per-cycle query (design.md sections 3 and 4), run beside today's gathers; nothing in the engine reads it yet ----

void HitRubber::FinderRun(void *const qv, const double r)
{
   FinderIndex &g = g_finder;
   FinderBallQ &q = *static_cast<FinderBallQ *>(qv);
   const double r2 = r * r;
   if (q.per.size() != g.objs.size())
      q.per.assign(g.objs.size(), std::vector<int>());
   for (const int k : q.objs)
      q.per[k].clear();   // the last result's
   q.ver++;   // P3a: the lists change: every table of the ball starts anew
   q.objs.clear();
   q.first.clear();
   q.faces.clear();
   thread_local std::vector<int> hit;
   hit.clear();
   QueryTree(g.nodes, g.order, [&g](const int i, RubberContact::V3 &lo, RubberContact::V3 &hi) { lo = g.objs[i].lo; hi = g.objs[i].hi; }, q.c0, r2, hit);
   std::sort(hit.begin(), hit.end());   // load order (R9)
   for (const int k : hit)
   {
      const HitRubber *const h = g.objs[k].host;
      if (!h->m_enabled)
         continue;   // not collidable now: skipped by its flag
      const RubberSurface &s = *h->m_surface;
      const size_t b = q.faces.size();
      QueryTree(s.fnodes, s.forder, [&s](const int i, RubberContact::V3 &lo, RubberContact::V3 &hi) { lo = s.tris[i].lo; hi = s.tris[i].hi; }, q.c0, r2, q.faces);
      if (q.faces.size() == b)
         continue;
      std::sort(q.faces.begin() + b, q.faces.end());   // load order (R9)
      q.objs.push_back(k);
      q.first.push_back((int)b);
      q.per[k].assign(q.faces.begin() + b, q.faces.end());
   }
   q.first.push_back((int)q.faces.size());
   // P4b: the playfield's cuts in the same query, pieces of the playfield host with their own boxes. The cut readers reach R + reach from a ball anywhere inside the
   // volume (within ext of c0, r = R + ext), so the cuts are found within r + reach: no re-query in the usual case.
   if (PlayfieldCuts::Count() > 0)
      FinderCuts(q.cq, q.c0, r + q.reach);
}

static bool CutProfOn() { return ContactProf::CheckOn() || FinderCheckOn(); }   // P4b counters (diagnostics only)

void HitRubber::FinderCuts(PlayfieldCuts::CutQuery &cq, const RubberContact::V3 &c0, const double r)
{
   FinderIndex &g = g_finder;
   const bool prof = CutProfOn();
   const uint64_t t0 = prof ? __rdtsc() : 0;
   const uint64_t ver = PlayfieldCuts::Version();
   const auto box = [](const int i, RubberContact::V3 &lo, RubberContact::V3 &hi) { PlayfieldCuts::Box(i, lo, hi); };
   if (g.cutVer != ver)   // the cuts' box tree, once per cut Setup (4 cuts per leaf, as the faces)
   {
      const int n = PlayfieldCuts::Count();
      g.cutNodes.clear();
      g.cutOrder.resize(n);
      for (int i = 0; i < n; i++)
         g.cutOrder[i] = i;
      if (n > 0)
         BuildTree(g.cutNodes, g.cutOrder, box, 0, n, 4);
      g.cutVer = ver;
   }
   cq.c0 = c0;
   cq.r = r;
   cq.ver = ver;
   cq.cuts.clear();
   QueryTree(g.cutNodes, g.cutOrder, box, c0, r * r, cq.cuts);
   std::sort(cq.cuts.begin(), cq.cuts.end());   // index order (R9)
   if (prof)
   {
      PlayfieldCuts::CutStats &st = PlayfieldCuts::Stats();
      st.queries++;
      st.queryTicks += __rdtsc() - t0;
   }
}

PlayfieldCuts::CutQuery *HitRubber::FinderCutsOf(const SurfaceContact &it)
{
   if (!g_finder.on)
      return nullptr;   // key off: the scan
   FinderBallQ *const q = static_cast<FinderBallQ *>(it.finderQ);
   if (q != nullptr)
      return &q->cq;
   thread_local PlayfieldCuts::CutQuery none;   // no query for this item: NearCuts queries around the ball (counted as a re-query)
   none.r = -1.0;
   return &none;
}

// P4b counters (VPX_CHECK_PROFILE=1 or VPX_FINDER_CHECK=1; diagnostics only): once a second of table time, key on or off
static void FinderP4Log(const uint64_t cycles)
{
   static uint32_t next = 0;
   const uint32_t now = g_pplayer != nullptr ? g_pplayer->m_time_msec : 0;
   if (!CutProfOn() || now < next || cycles == 0)
      return;
   next = now + 1000;
   static const uint64_t tsc0 = __rdtsc();
   static const auto c0 = std::chrono::steady_clock::now();
   const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - c0).count();
   const uint64_t dt = __rdtsc() - tsc0;
   if (dt == 0 || s <= 0.0)
      return;
   const double nsPerTick = s * 1e9 / (double)dt, cyc = (double)cycles;
   const PlayfieldCuts::CutStats &c = PlayfieldCuts::Stats();
   PLOGI << "Finder P4b totals: " << cycles << " cycles (key " << (g_finder.on ? 1 : 0) << "); NearCuts calls " << c.calls << " (" << (double)c.calls / cyc << " per cycle), cuts tested per call "
         << (c.calls > 0 ? (double)c.tested / (double)c.calls : 0.0) << ", us per cycle " << (double)c.ticks * nsPerTick * 1e-3 / cyc << ", re-queries " << c.requery << "; cut queries " << c.queries
         << ", us per cycle " << (double)c.queryTicks * nsPerTick * 1e-3 / cyc << "; check calls " << c.checks << ", mismatches " << c.mismatch << ", misses " << c.miss;
}

void HitRubber::FinderQuery(const HitBall *const ball, const RubberContact::V3 &c, const RubberContact::V3 &v, const double ballR, const double dt)
{
   FinderIndex &g = g_finder;
   if (!g.on)
      return;
   const auto t0 = std::chrono::steady_clock::now();
   FinderBallQ &q = g.balls[ball];
   if (q.queried)   // the last query's re-queries, per ball-cycle
   {
      g.qCycles++;
      g.reqCycles += q.reqThis > 0 ? 1 : 0;
      g.reqTotal += q.reqThis;
      g.reqWorst = std::max(g.reqWorst, q.reqThis);
   }
   // a: the ball's acceleration from its net force in the previous cycle, taken as its velocity change over that cycle (the engine keeps no net force per ball);
   // none known (the ball was not engaged last cycle): 0, and the inside check covers it
   double accel = 0.0;
   if (q.queried && q.vPrevCycle + 1 == g.cycle)
      accel = (v - q.vPrev).Length() / dt;
   else
      g.noPrev++;
   q.vPrev = v;
   q.vPrevCycle = g.cycle;
   q.queried = true;
   q.c0 = c;
   q.R = ballR;
   q.reach = ContactReachM(v.Length(), accel, dt);
   q.ext = q.reach;
   q.reqThis = 0;
   q.items = nullptr;
   FinderRun(&q, ballR + q.reach);
   const double ns = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count();
   g.queries++;
   g.objsFound += q.objs.size();
   g.facesFound += q.faces.size();
   g.qNs.Add(ns);
   if (!q.objs.empty())
      g.qObjNs.Add(ns / (double)q.objs.size());
}

void HitRubber::FinderInside(const SurfaceContact &it, const RubberContact::V3 &c) const
{
   const ContactProf::Ticks profInside(ContactProf::CD().inside);   // P3a profile (diagnostics only)
   FinderBallQ *const q = static_cast<FinderBallQ *>(it.finderQ);
   if (q == nullptr || !g_finder.on)
      return;
   const double stray = (c - q->c0).Length();
   if (stray <= q->ext)
      return;
   const auto t0 = std::chrono::steady_clock::now();
   q->ext = 2.0 * stray;   // the factor 2 bounds the repeats; completeness comes from this check (design.md section 4)
   FinderRun(q, q->R + q->ext);   // (P2b: every item of the ball reads its object's faces in place, so the lists are refreshed before this evaluation measures)
   q->reqThis++;
   if (q->items != nullptr)   // P2b: an object the re-query found that is not among the ball's surfaces (engagement comes from the cycle's query; counted, reported in the totals)
      for (const int k : q->objs)
         if (std::none_of(q->items->begin(), q->items->end(), [k](const SurfaceContact &o) { return o.r->FinderObjIndex() == k; }))
            q->reqNew++;
   q->reqNs += std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count();
}

RubberContact::FaceTable *HitRubber::FinderTable(const SurfaceContact &it, const RubberContact::V3 &c) const
{
   const ContactProf::Ticks profBegin(ContactProf::CD().begin);   // P3a profile (diagnostics only)
   FinderBallQ *const q = static_cast<FinderBallQ *>(it.finderQ);
   if (q == nullptr || !g_finder.on || it.qFaces == nullptr || m_surface == nullptr)
      return nullptr;
   const int k = FinderObjIndex();
   if (k < 0 || k >= (int)q->per.size() || it.qFaces != &q->per[k])
      return nullptr;
   if (q->tab.size() != q->per.size())
      q->tab.resize(q->per.size());
   RubberContact::FaceTable &t = q->tab[k];
   t.check = FinderCheckOn();
   t.count = t.check || ContactProf::CheckOn();   // P3a's per-read counters: check or profile mode only
   t.direct = !FinderTableOn();
   if (t.tag < 0)
      t.tag = Name() == "zCol_Rubber_Corner_15" || Name() == "zCol_Rubber_Corner_16" || Name() == "zCol_Rubber_Corner_21" ? 1 : 0;
   t.Begin(c, q->ver, q->per[k].size());
   t.useNear = !g_nearOff;   // P3b: the near pass, once per ball position (check mode's whole-list pass reads the list as P3a)
   if (t.useNear)
   {
      const uint64_t m0 = t.nearMiss, b0 = t.nearBox, v0 = t.nearVisit, f0 = t.fills;
      const uint32_t w0 = t.nearGen;
      const uint64_t tw0 = t.count ? __rdtsc() : 0;
      t.NearPass(m_surface->tris, m_surface->nb, q->per[k], c, q->R * q->R);
      if (t.count && t.nearGen != w0)   // the near pass's counters: check or profile mode only
      {
         const uint64_t tw1 = __rdtsc();
         if (t.objStats == nullptr)
            t.objStats = &g_tab.perObj[Name()];   // std::map: the entry's address stays valid
         FinderTabStats::ObjNear &o = *static_cast<FinderTabStats::ObjNear *>(t.objStats);
         o.ticks += tw1 - tw0;
         o.passes++;
         o.list += q->per[k].size();
         o.visit += t.nearVisit - v0;
         g_tab.nearBox += t.nearBox - b0;
         g_tab.nearVisit += t.nearVisit - v0;
         g_tab.nearFills += t.fills - f0;
         g_tab.passes++;
         g_tab.nearIn += t.inI.size();
         g_tab.nearList += q->per[k].size();
      }
      if (t.nearMiss != m0)
      {
         g_tab.nearMiss += t.nearMiss - m0;
         if (g_tab.loggedNear < 50)
         {
            g_tab.loggedNear++;
            const int face = t.lastMissJ >= 0 && t.lastMissJ < (int)q->per[k].size() ? q->per[k][t.lastMissJ] : -1;
            PLOGW << "Finder near pass miss: object " << Name() << ", face " << face << ", piece " << (face >= 0 && face < (int)m_surface->piece.size() ? m_surface->piece[face] : -1) << ", ball "
                  << (const void *)q->ball << ", cycle " << g_finder.cycle << ", c " << c.x << " " << c.y << " " << c.z;
         }
      }
   }
   return &t;
}

void HitRubber::FinderTabTake(uint64_t &checks, uint64_t &fills, uint64_t &reads)
{
   checks = g_tab.rowChecks;
   fills = g_tab.rowFills;
   reads = g_tab.rowReads;
   g_tab.rowChecks = g_tab.rowFills = g_tab.rowReads = 0;
}

void HitRubber::FinderTotals(const bool force)
{
   FinderIndex &g = g_finder;
   if (!force)
      g.cycle++;
   FinderP4Log(g.cycle);   // P4b counters (profile or check mode)
   if (!g.on)
      return;
   const uint32_t now = g_pplayer != nullptr ? g_pplayer->m_time_msec : 0;
   if (!force && now < g.nextTotalsMs)
      return;
   g.nextTotalsMs = now + 1000;
   if (g.queries == 0)
      return;
   double reqNs = 0.0;
   g.reqNewObj = 0;
   for (const auto &e : g.balls)
   {
      reqNs += e.second.reqNs;
      g.reqNewObj += e.second.reqNew;
   }
   PLOGI << "Finder query totals: " << g.queries << " ball-cycle queries (" << g.noPrev << " with no previous velocity), objects found " << (double)g.objsFound / (double)g.queries << " and faces "
         << (double)g.facesFound / (double)g.queries << " per query; time per query mean " << g.qNs.sum / (double)g.qNs.n << " ns, p99 " << g.qNs.P(0.99) << ", worst " << g.qNs.worst
         << "; per query and object found mean " << (g.qObjNs.n > 0 ? g.qObjNs.sum / (double)g.qObjNs.n : 0.0) << " ns, p99 " << g.qObjNs.P(0.99) << ", worst " << g.qObjNs.worst << "; re-queries "
         << g.reqTotal << " in " << g.reqCycles << " of " << g.qCycles << " ball-cycles (" << (g.qCycles > 0 ? 100.0 * (double)g.reqCycles / (double)g.qCycles : 0.0) << " %), worst " << g.reqWorst
         << " in one ball-cycle, time " << reqNs << " ns in all; P2b: queries at engagement " << g.queriesAt << ", objects a re-query found that were not among the ball's surfaces " << g.reqNewObj;
   if (g_tab.checks > 0)   // P3a: the measurement tables
   {
      const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - g_tab.t0).count();
      const double nsPerTick = s > 0.0 ? s * 1e9 / (double)(__rdtsc() - g_tab.tsc0) : 0.0;
      PLOGI << "Finder table totals: " << g_tab.checks << " contact checks; faces measured " << (double)g_tab.fills / (double)g_tab.checks << " and table reads " << (double)g_tab.reads / (double)g_tab.checks
            << " per check (" << g_tab.fills << ", " << g_tab.reads << "); checks per ball, object and cycle median " << HistP(g_tab.perCyc, 0.5) << ", p90 " << HistP(g_tab.perCyc, 0.9) << ", worst "
            << g_tab.perCycWorst << "; time per check mean " << nsPerTick * (double)g_tab.tickSum / (double)g_tab.checks << " ns, p99 " << nsPerTick * 16.0 * (double)(HistP(g_tab.ticks, 0.99) + 1)
            << ", worst " << nsPerTick * (double)g_tab.tickWorst << " (" << nsPerTick << " ns per tick); check-mode table mismatches " << g_tab.mism << "; P3b near pass misses " << g_tab.nearMiss
            << ", probe answers differing " << g_tab.probeDiff << ", box tests " << g_tab.nearBox << "; near passes " << g_tab.passes << ", per pass: faces past the box "
            << (g_tab.passes > 0 ? (double)g_tab.nearVisit / (double)g_tab.passes : 0.0) << ", measured " << (g_tab.passes > 0 ? (double)g_tab.nearFills / (double)g_tab.passes : 0.0) << ", box tests "
            << (g_tab.passes > 0 ? (double)g_tab.nearBox / (double)g_tab.passes : 0.0) << ", handed over " << (g_tab.passes > 0 ? (double)g_tab.nearIn / (double)g_tab.passes : 0.0) << ", query list "
            << (g_tab.passes > 0 ? (double)g_tab.nearList / (double)g_tab.passes : 0.0) << "; rubber corners " << g_tab.cornerN << " checks, mean "
            << (g_tab.cornerN > 0 ? nsPerTick * (double)g_tab.cornerSum / (double)g_tab.cornerN : 0.0) << " ns, p99 " << nsPerTick * 16.0 * (double)(HistP(g_tab.cornerTicks, 0.99) + 1)
            << ", time " << nsPerTick * (double)g_tab.cornerSum * 1e-3 << " us in all; P3c edge slots: measured " << g_tab.efills << ", reads " << g_tab.ereads << ", check-mode edge mismatches " << g_tab.emism << "; cut B check-mode mismatches " << g_tab.cbMism;
      if (ContactProf::CheckOn())   // P3a profile (diagnostics only)
         ContactProf::CheckDump(g.cycle - ContactProf::CD().cyc0, nsPerTick);
      if (!g_tab.perObj.empty())   // P3b2: the 5 objects with the longest query list per pass (check or profile mode)
      {
         std::vector<std::pair<double, std::string>> v;
         for (const auto &e : g_tab.perObj)
            v.emplace_back(-(double)e.second.list / (double)e.second.passes, e.first);
         std::sort(v.begin(), v.end());
         std::ostringstream os;
         for (size_t i = 0; i < v.size() && i < 5; i++)
         {
            const FinderTabStats::ObjNear &o = g_tab.perObj[v[i].second];
            os << "; " << v[i].second << " passes " << o.passes << ", list " << (double)o.list / (double)o.passes << ", past the box " << (double)o.visit / (double)o.passes << ", ns per pass "
               << nsPerTick * (double)o.ticks / (double)o.passes;
         }
         PLOGI << "Finder near pass densest (" << g_tab.perObj.size() << " objects)" << os.str();
      }
   }
}

// ---- The shared finder, piece P2b: the query finds every static mesh contact ([Features] shared_finder = 1) ----

bool HitRubber::FinderDrives()
{
   return g_finder.on;
}

void HitRubber::FinderQueryAt(const HitBall *const ball, const RubberContact::V3 &c, const double ballR, const double extM)
{
   FinderIndex &g = g_finder;
   if (!g.on)
      return;
   FinderBallQ &q = g.balls[ball];
   q.c0 = c;
   q.R = ballR;
   q.reach = q.ext = extM;
   q.items = nullptr;
   FinderRun(&q, ballR + extM);
   g.queriesAt++;
}

void HitRubber::FinderBindOne(const HitBall *const ball, SurfaceContact &it)
{
   FinderIndex &g = g_finder;
   if (!g.on)
      return;
   FinderBallQ &q = g.balls[ball];
   it.finderQ = &q;
   q.ball = ball;
   q.ver++;   // P3a: a rebind starts the ball's tables anew
   const int k = it.r->FinderObjIndex();
   if (k < 0)
      return;
   if (q.per.size() != g.objs.size())
      q.per.assign(g.objs.size(), std::vector<int>());
   it.qFaces = &q.per[k];
   it.candR = q.R + q.reach;   // how far the list reaches from the query's centre (>= 0: a list is held; the cuts read it as their reach)
}

void HitRubber::FinderBind(const HitBall *const ball, std::vector<SurfaceContact> &items, const bool cycleItems)
{
   if (!g_finder.on)
      return;
   for (SurfaceContact &it : items)
      FinderBindOne(ball, it);
   g_finder.balls[ball].items = cycleItems ? &items : nullptr;
}

void HitRubber::FinderHosts(const HitBall *const ball, std::vector<ContactHost *> &out)
{
   out.clear();
   FinderIndex &g = g_finder;
   if (!g.on)
      return;
   const FinderBallQ &q = g.balls[ball];
   for (const int k : q.objs)
      out.push_back(const_cast<HitRubber *>(g.objs[k].host));
}

// P3a profile (diagnostics only, VPX_CHECK_PROFILE=1): one log line "Check profile" with the split of a contact check's time and the fit time = a + b fills + c reads
void ContactProf::CheckDump(const uint64_t cycles, const double nsPerTick)
{
   const CheckData &d = CD();
   if (d.n == 0 || cycles == 0)
      return;
   const double cyc = (double)cycles, k = nsPerTick;
   const auto fit = [&](const double sY, const double sYF, const double sYR, double &a, double &b, double &c) {   // least squares y = a + b F + c R (normal equations, Cramer)
      const double n = (double)d.n;
      const double m[3][3] = { { n, d.sF, d.sR }, { d.sF, d.sFF, d.sFR }, { d.sR, d.sFR, d.sRR } }, y[3] = { sY, sYF, sYR };
      const auto det = [](const double q[3][3]) { return q[0][0] * (q[1][1] * q[2][2] - q[1][2] * q[2][1]) - q[0][1] * (q[1][0] * q[2][2] - q[1][2] * q[2][0]) + q[0][2] * (q[1][0] * q[2][1] - q[1][1] * q[2][0]); };
      const double D = det(m);
      double r[3];
      for (int col = 0; col < 3; col++)
      {
         double q[3][3];
         for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
               q[i][j] = j == col ? y[i] : m[i][j];
         r[col] = D != 0.0 ? det(q) / D : 0.0;
      }
      a = r[0] * k; b = r[1] * k; c = r[2] * k;
   };
   double ta, tb, tc, la, lb, lc;
   fit(d.sT, d.sTF, d.sTR, ta, tb, tc);
   fit(d.sL, d.sLF, d.sLR, la, lb, lc);
   std::ostringstream o;
   o << "Check profile: " << cycles << " cycles; table checks " << (double)d.n / cyc << " per cycle, whole " << (double)d.whole * k / cyc * 1e-3 << " us/cycle (" << (double)d.whole * k / (double)d.n
     << " ns/check); no-table checks " << (double)d.nNoTab / cyc << " per cycle, " << (double)d.wholeNoTab * k / cyc * 1e-3 << " us/cycle; split us/cycle: inside " << (double)d.inside * k / cyc * 1e-3
     << ", table begin " << (double)d.begin * k / cyc * 1e-3 << ", face loops in checks " << (double)d.loop * k / cyc * 1e-3 << " (" << (double)d.loops / cyc << " loops, " << (double)d.loopFaces / cyc
     << " faces looped per cycle), loops outside checks " << (double)d.loopOut * k / cyc * 1e-3 << " (" << (double)d.loopFacesOut / cyc << " faces); fills " << (double)d.fills / cyc << ", reads "
     << (double)d.reads / cyc << " per cycle; fit whole ns = " << ta << " + " << tb << " fills + " << tc << " reads; fit loops ns = " << la << " + " << lb << " fills + " << lc << " reads; by stage (probe/patch/arrival checks per cycle, us/cycle):";
   static const char *const nm[StCount] = { "other", "ContactHolds", "GroupShape", "RingHandIn", "RingCap", "LoadState", "CordNear", "EngageNear", "SeatAtTouch", "EngageImpl", "WriteBack", "ReleaseTest" };
   for (int st = 0; st < StCount; st++)
      if (d.st[st][0] + d.st[st][1] + d.st[st][2] > 0)
         o << " " << nm[st] << " " << (double)d.st[st][0] / cyc << "/" << (double)d.st[st][1] / cyc << "/" << (double)d.st[st][2] / cyc << " "
           << (double)(d.stT[st][0] + d.stT[st][1] + d.stT[st][2]) * k / cyc * 1e-3;
   o << "; P3c BarsTriangles " << (double)d.barsN / cyc << " calls, " << (double)(d.bars + d.barsCorner) * k / cyc * 1e-3 << " us/cycle (rubber corners " << (double)d.barsCorner * k / cyc * 1e-3
     << "), edges measured " << (double)d.eFills / cyc << " of " << (double)d.eReads / cyc << " requested per cycle";
   PLOGI << o.str();
   // mc-cost profile (diagnostics only): BarsTriangles split and per object
   const auto us = [&](const uint64_t t) { return (double)t * k / cyc * 1e-3; };
   std::ostringstream b;
   b << "Bars split us/cycle: probe " << us(d.bProbe) << ", face pass " << us(d.bFace) << ", edge pass " << us(d.bEdge) << ", MergeBars " << us(d.bMerge) << ", Finish " << us(d.bFinish)
     << ", whole " << us(d.bars + d.barsCorner) << "; calls per cycle " << (double)d.barsN / cyc << ": end after probe " << (double)d.bRetProbe / cyc << ", fewer than 2 nodes " << (double)d.bRetNodes / cyc
     << ", one peak " << (double)d.bOnePeak / cyc << ", more peaks " << (double)d.bMulti / cyc << "; per call past the probe: candidates "
     << (double)d.bCand / (double)std::max<uint64_t>(1, d.barsN - d.bRetProbe) << ", face nodes " << (double)d.bFaceNodes / (double)std::max<uint64_t>(1, d.barsN - d.bRetProbe)
     << ", nodes at merge " << (double)d.bNodes / (double)std::max<uint64_t>(1, d.bOnePeak + d.bMulti) << "; cut B per cycle: skipped " << (double)d.cbSkip / cyc << ", full law " << (double)d.cbFull / cyc;
   std::vector<std::pair<uint64_t, std::string>> v;
   uint64_t tObj = 0;
   for (const auto &e : g_tab.perObj)
   {
      v.emplace_back(e.second.bars.t, e.first);
      tObj += e.second.bars.t;
   }
   std::sort(v.rbegin(), v.rend());
   b << "; objects named " << us(tObj) << " us/cycle of it; top:";
   for (size_t i = 0; i < v.size() && i < 10; i++)
   {
      const ObjBars &ob = g_tab.perObj[v[i].second].bars;
      b << " " << v[i].second << " " << us(ob.t) << " us/cycle, " << (double)ob.n / cyc << " calls, " << (double)ob.cand / (double)std::max<uint64_t>(1, ob.n) << " cand/call, "
        << (double)ob.nodes / (double)std::max<uint64_t>(1, ob.n) << " face nodes/call, " << (double)ob.full / cyc << " full-law calls;";
   }
   PLOGI << b.str();
}

ContactProf::ObjBars *ContactProf::ObjBarsOf(void *objStats)   // mc-cost profile (diagnostics only)
{
   return objStats != nullptr ? &static_cast<FinderTabStats::ObjNear *>(objStats)->bars : nullptr;
}
