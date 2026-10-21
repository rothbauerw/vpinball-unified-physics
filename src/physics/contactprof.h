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
      StepBalls, GatherBalls, Regather, RefreshCycle, ApproachSpeed, GatherFree, FindIslands, SetupIsland, Substeps,
      SurfaceForces, SurfProbe, SurfEval, SurfReact, BatForces, BandForces, PairForces, Integrate,
      WriteBack, ReleaseTest, Finish, EngageNear, EngageImpl, ProbeBody, NextStep, Sources, Count
   };
   inline const char *Name(int i)
   {
      static const char *const n[Count] = { "StepBalls", "GatherBalls", "  Regather", "  RefreshCycle", "  ApproachSpeed", "GatherFree", "FindIslands", "SetupIsland", "Substeps",
                                            "  SurfaceForces", "    Probe", "    EvaluateTube", "    ApplyReaction", "  BatForces", "  BandForces", "  PairForces", "  Integrate",
                                            "WriteBack", "  ReleaseTest", "Finish", "EngageNear", "EngageImpl", "  ProbeBody", "  NextStep", "    Sources (gravity, magnet pull, eddy)" };
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
   void Cycle();                             // one physics cycle began
   void Split(int cause, bool staticTime, const char *what);   // a cycle was cut short at a hit: the type of the object that ended the sub-step (100 = a flipper stop, 200 = the static-time clamp)
}
#define CPROF(id) ContactProf::Scope cprofScope_##id(ContactProf::id)
