#pragma once
// Unified contact system (docs/UnifiedContactDesign.txt). One call per 1 ms physics cycle steps every ball it owns.
// Step 1: the entry point, the engine gate and the cycle's step (gather, islands, per-island substeps, write-back) live here.

#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <string>
#include <utility>
#include <vector>

class HitRubber;
class ContactHost;
class HitBall;
struct BallS;
class ContactSurface;
class FlipperMoverObject;

namespace ContactSystem
{
   // An environment variable by its current name, else by the name it had before the rename (VPX_RING_TRACE -> VPX_CONTACT_TRACE, VPX_RING_NEAR -> VPX_CONTACT_NEAR, VPX_RING_NOSIDE -> VPX_CONTACT_NOSIDE)
   inline const char *EnvVar(const char *name, const char *oldName)
   {
      const char *const e = getenv(name);
      return e != nullptr ? e : getenv(oldName);
   }

   constexpr double kMetresPerUnit = 0.53975e-3;                  // m per VPX unit
   constexpr double kVelocityScale = kMetresPerUnit / 0.01;       // VPX velocity (U/VPT, 1 VPT = 10 ms) -> m/s
   constexpr double kMassUnitKg = 0.080;                          // ball mass 1.0 in VPX = 80 g

   // Reads the table's physics ini (PhysicsIni::FindUnified) once per table, before any part's PhysicSetup; [Features] unified = 1 switches the unified engine on.
   void SetTable(const std::filesystem::path &tableFile);

   // The table runs on the unified engine: only with  unified = 1  in the table's physics ini; no ini or no key = stock VPX physics.
   bool Unified();

   // One cycle's step of every ball the contact system holds (dtime in the engine's step units). Its events (ball-ball callbacks, part hits, the parts' own events) are only queued here.
   void Step(float dtime);

   // The state of everything an overlap replay would have to put back (docs/UnifiedContactDesign.txt section 11, point 2), taken at the start of a sub-interval, before Step: every ball's position,
   // velocity and spin, the engaged balls' contacts (shear springs, touching flags, counters), the ball pairs, the saved state of every moving part, the size of the event queue. No allocation after the
   // first few calls. RestoreStart puts it all back and discards the events queued since. LastSaveUs: what SaveStart cost, for the cost log.
   void SaveStart();
   void RestoreStart();
   unsigned LastSaveUs();
   // Test mode (env VPX_SNAPSHOT_TEST=1; =2 is the negative control that leaves the parts out of the restore): after every Step the snapshot is restored and the step is taken again; any difference in the result is logged. Islands with a bat or a band are skipped
   // (their state is not saved yet).
   bool SnapshotTestOn();
   void SnapshotTest(float dtime);

   // The overlap check and the replay of the sub-interval (UnifiedContactDesign.txt section 11, point 2): after UpdateDisplacements, before Dispatch. A free ball found inside a moving part is engaged at
   // the time its gap was zero, with the part where it really was, and the sub-interval is stepped again from there. VPX_OVERLAP_REPLAY=0 switches it off (test). TakeReplayStats: the replays and their cost
   // since the last call, for the cost log.
   void CheckOverlaps(float hittime);
   void TakeReplayStats(unsigned &count, unsigned &us);

   // Fires the events queued by Step and by the parts, in the order they happened. Called by the engine after the free balls have moved (UpdateDisplacements), the order of stock VPX: displacements,
   // then collisions and script events (docs/UnifiedContactDesign.txt section 11, point 2). Nothing has reached a script when a late engagement has to rewind the interval.
   void Dispatch();

   // What the last Step did (for the cost log, VPX_CONTACT_TIMING): bodies stepped, islands, substeps summed over the islands, bats stepped
   struct CycleStats
   {
      int bodies = 0, islands = 0, substeps = 0, bats = 0, bandBodies = 0, restBodies = 0, maxSteps = 0;
      int rewinds = 0;                  // substeps taken again because a contact started deeper than it should (arrival-time loop)
      int items = 0;                    // surfaces summed over the bodies of the cycle
      long long probes = 0;             // surface probes of the cycle: items x substeps summed over the bodies
      int gatherUs = 0, setupUs = 0, substepUs = 0, writeUs = 0;
   };
   const CycleStats &LastCycle();

   // An event raised inside the step: fired after it (the script never runs inside the step)
   void DeferEvent(std::function<void()> f);

   // Takes the free balls that can touch a surface this cycle into the contact (before the stock hit search); [Features] proximity_engage, VPX_PROXIMITY=0|1
   void EngageNear(float dtime);

   // Engagement: a surface (the host) takes a ball into the force model; the ball stays until it is clear of every surface and receding
   bool Engage(ContactHost *host, HitBall *pb, const char *how);   // false = the stock collision keeps the ball
   bool EngageBand(HitRubber *band, HitBall *pb);               // a band takes the ball: its nodes join the ball's contact (the ball stays in the other contacts it has)
   void StepCordAlone(HitRubber *band, float dtime);   // a cord with no ball on it this cycle (its plate moving by itself, a ball just gone)
   void Release(HitBall *ball, const char *why);                  // no-op if the contact system does not hold the ball
   void ReleaseHost(ContactHost *host);                            // the host is destroyed
   bool Holds(const HitBall *ball);
   bool Owns(const BallS &ball, const ContactSurface *surface);  // the contact that holds the ball includes this surface
   size_t Count();                                               // balls held

   // Islands: the connected components of the interaction graph (bodies are indices 0..n-1, an interaction is an edge). Each island is stepped on its own substeps.
   // Islands are listed in order of their root index; members ascend; edges keep the order they were given in.
   class Islands
   {
   public:
      struct Island
      {
         std::vector<int> members;
         std::vector<std::pair<int, int>> edges;
      };
      std::vector<Island> list;   // only the first count entries are live this cycle (the rest keep their buffers)
      size_t count = 0;
      void Build(int n, const std::vector<std::pair<int, int>> &edges, const std::vector<std::pair<int, int>> &links = {});   // buffers are reused between cycles
   private:
      std::vector<int> m_parent;
      int Root(int i);
   };
}
