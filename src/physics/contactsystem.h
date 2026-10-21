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

   // One cycle's step of every ball the contact system holds (dtime in the engine's step units).
   void Step(float dtime);

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
   bool Engage(HitRubber *host, HitBall *pb, const char *how);   // false = the stock collision keeps the ball
   bool EngageBat(FlipperMoverObject *bat, HitBall *pb);         // a bat takes the ball (contact-held or free): its rubber joins the ball's contact
   bool EngageBand(HitRubber *band, HitBall *pb);               // a band takes the ball: its nodes join the ball's contact (the ball stays in the other contacts it has)
   void Release(HitBall *ball, const char *why);                  // no-op if the contact system does not hold the ball
   void ReleaseHost(HitRubber *host);                            // the host is destroyed
   void ReleaseBat(FlipperMoverObject *bat);                     // the bat is destroyed: every ball on it lets go
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
