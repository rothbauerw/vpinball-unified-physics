#pragma once
// A scoped-timer profiler for the contact system (VPX_CONTACT_PROFILE=1). Time is the CPU cycle counter, converted to microseconds against the wall clock when the table is
// dumped; scopes nest, so a time is INCLUSIVE of the scopes inside it. Every 5000 step calls a table goes to the log: lines "PROF ...".
#include <intrin.h>
#include <cstdint>
#include <cstdlib>
#include <string>

namespace ContactProf
{
   enum Id
   {
      StepBalls, GatherBalls, Regather, RefreshCycle, GatherFree, FindIslands, BatForces, BandForces, WriteBack, ReleaseTest,
      Finish, EngageNear, EngageImpl, PartWorld, PartSwitch, Count
   };
   inline const char *Name(int i)
   {
      static const char *const n[Count] = { "StepBalls", "GatherBalls", "  Regather", "  RefreshCycle", "GatherFree", "FindIslands", "  BatForces", "  BandForces", "WriteBack", "  ReleaseTest", "Finish", "EngageNear", "EngageImpl", "    PartWorld (a gate / spinner's geometry update, once per island substep)", "    PartSwitch (a gate / spinner's switch update, once per own substep)" };
      return n[i];
   }
   inline bool On()
   {
      static const bool on = getenv("VPX_CONTACT_PROFILE") != nullptr;
      return on;
   }
   struct Data
   {
      uint64_t t[Count] = {}, n[Count] = {};
   };
   inline Data &D()
   {
      static Data d;
      return d;
   }
   struct Scope
   {
      int id;
      uint64_t t0;
      explicit Scope(int i) : id(i), t0(On() ? __rdtsc() : 0) {}
      ~Scope()
      {
         if (t0 != 0)
         {
            D().t[id] += __rdtsc() - t0;
            D().n[id]++;
         }
      }
   };
   void Dump(const char *why);               // contactsystem.cpp
   void Note(const std::string &what);      // count an event by name
   bool ProbeTime(const void *surface, uint64_t ticks);   // add the time of one probe of a surface; true the first time a surface is seen
   void ProbeName(const void *surface, const std::string &name);
   void Cycle();                             // one physics cycle began
   void Split(int cause, bool staticTime, const char *what);   // a cycle was cut short at a hit: the type of the object that ended the sub-step (100 = a flipper stop, 200 = the static-time clamp)
}
// P3a profile (diagnostics only, VPX_CHECK_PROFILE=1; runs/2026-11-12-finder-P3a-profile): where one contact check's time goes. A check is the outermost of HitRubber's Probe / ProbePatch / ArrivalS;
// its time splits into the inside check (FinderInside), the table lookup and Begin (FinderTable), the face loops (each loop timed whole: reads and fills together) and the rest. Per-check sums
// of time, faces measured (fills) and table reads give a least-squares fit time = a + b fills + c reads. Checks are counted by the calling stage (Stage, set by the callers).
namespace ContactProf
{
   inline bool CheckOn()
   {
      static const bool on = getenv("VPX_CHECK_PROFILE") != nullptr;
      return on;
   }
   enum StageId { StOther, StHolds, StGroupShape, StRingHandIn, StRingCap, StLoadState, StCordNear, StEngageNear, StSeatAtTouch, StEngageImpl, StWriteBack, StReleaseTest, StCount };
   struct CheckData
   {
      uint64_t inside = 0, begin = 0, loop = 0, loopOut = 0, loopFaces = 0, loopFacesOut = 0, loops = 0;
      uint64_t n = 0, whole = 0, nNoTab = 0, wholeNoTab = 0, fills = 0, reads = 0;
      uint64_t bars = 0, barsCorner = 0, barsN = 0, eFills = 0, eReads = 0;   // P3c: time in BarsTriangles (probe inside included; the rubber corners apart) and its calls; edges measured and edge requests
      // mc-cost profile (diagnostics only): BarsTriangles split in ticks (probe, face pass, edge pass, MergeBars, Finish); calls by where they end (after the probe, fewer than 2 nodes,
      // one peak, more peaks); candidates, face nodes and all nodes summed; bars calls and ticks since the last cost-log row (VPX_CONTACT_TIMING)
      uint64_t bProbe = 0, bFace = 0, bEdge = 0, bMerge = 0, bFinish = 0, bRetProbe = 0, bRetNodes = 0, bOnePeak = 0, bMulti = 0, bCand = 0, bFaceNodes = 0, bNodes = 0, rowBarsN = 0, rowBars = 0;
      uint64_t cbSkip = 0, cbFull = 0;   // cut B (runs/2026-11-12-mc-cutB, profile only): calls past the probe that the rule skipped, and that ran the full law
      double sT = 0, sL = 0, sF = 0, sR = 0, sFF = 0, sRR = 0, sFR = 0, sTF = 0, sTR = 0, sLF = 0, sLR = 0;
      uint64_t st[StCount][3] = {}, stT[StCount][3] = {};
      int depth = 0;
      uint64_t cyc0 = ~0ull;
   };
   inline CheckData &CD()
   {
      static CheckData d;
      return d;
   }
   inline int &StageCur()
   {
      static int s = StOther;
      return s;
   }
   struct Stage   // the calling stage of the checks inside this scope
   {
      const int prev;
      explicit Stage(const int s) : prev(StageCur()) { StageCur() = s; }
      ~Stage() { StageCur() = prev; }
   };
   struct Ticks   // add a scope's ticks to a counter (only while a check is open, else to 'out')
   {
      uint64_t *const to;
      const uint64_t t0;
      explicit Ticks(uint64_t &dst) : to(CheckOn() ? &dst : nullptr), t0(to != nullptr ? __rdtsc() : 0) {}
      ~Ticks() { if (to != nullptr) *to += __rdtsc() - t0; }
   };
   struct LoopTicks   // one face loop, timed whole
   {
      const bool on;
      const uint64_t t0;
      const size_t nf;
      explicit LoopTicks(const size_t faces) : on(CheckOn()), t0(on ? __rdtsc() : 0), nf(faces) {}
      ~LoopTicks()
      {
         if (!on)
            return;
         const uint64_t dt = __rdtsc() - t0;
         CheckData &d = CD();
         if (d.depth > 0) { d.loop += dt; d.loopFaces += nf; d.loops++; }
         else { d.loopOut += dt; d.loopFacesOut += nf; }
      }
   };
   void CheckDump(uint64_t cycles, double nsPerTick);   // hitrubber.cpp
   struct ObjBars   // mc-cost profile (diagnostics only): one object's BarsTriangles ticks, calls, candidates and nodes
   {
      uint64_t t = 0, n = 0, cand = 0, nodes = 0, full = 0;   // full: cut B, calls that ran the full law
   };
   ObjBars *ObjBarsOf(void *objStats);   // hitrubber.cpp: the per-object entry behind a FaceTable's objStats (nullptr if none)
}
#define CPROF(id) ContactProf::Scope cprofScope_##id(ContactProf::id)
