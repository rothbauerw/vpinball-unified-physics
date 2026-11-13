#pragma once
// A surface a ball can touch (docs/unified-physics/archive/UnifiedContactDesign.txt, "Surface"): static geometry or a driven one (the turntable's disc). The contact system asks it for the nearest
// point, normal and curvature of the part of the surface near the ball (Probe), and, for a surface that moves, the velocity of the surface at the contact and the reaction
// on it. The contact law (Hertz + damping + Mindlin + twist + rolling resistance) and the integration are the same for every surface and live outside.
//
// Implemented by the hosts (contacthost.h): ContactSurfaceHost (a rubber tube, the faces of a wall, primitive, ramp or the playfield) and ContactPartHost (a moving part). The unified engine's own surfaces
// (playfield mesh, ramps, wall and primitive meshes) will implement this interface directly.

#include <string>
#include <vector>
#include "physics/rubbercontact.h"

class ContactSurface;
class HitObject;
class ContactMovingPart;
struct RubberObjectSettings;

// What the group step keeps of a touch between cycles, per surface: touching at the cycle's end, the shear spring's stretch then, the depth its squash counts from (GroupStep Item dZero), and the effective radius it
// took when it started: held for the life of the touch (ring_model; the law's energy is one function throughout).
struct SurfaceState
{
   bool touching = false;
   RubberContact::V3 xi;
   RubberContact::RollAngle ra;   // the rolling and twist springs' angles then (roll_model 1)
   double zero = 0.0, reff = 0.0;
   double kb = 0.0;         // the bend in series this touch started with (a part's surface)
   double ktb = 0.0;        // and its bend across the normal, in series with the shear spring (a wire switch's bent wire)
   double endDepth = 0.0;   // a part's surface: its depth at the cycle's end, from its zero (a fixed surface's is SurfaceContact::endDepth, the energy check's)
};

// What one ball has with ONE surface: the material parameters, the shear state kept between cycles, and the faces near the ball this cycle.
// The reach of a ball over one cycle (to-do N5, the one reach rule; shared finder design.md section 3): how far its centre can move in dt from its speed and its
// acceleration, reach = speed x dt + 1/2 x accel x dt^2 (m; speed m/s, accel m/s^2, dt s). Used first by the shared finder's per-cycle query (to-do 1.34 P2a).
inline double ContactReachM(const double speed, const double accel, const double dt) { return speed * dt + 0.5 * accel * dt * dt; }

struct SurfaceContact
{
   const ContactSurface *r = nullptr;
   RubberContact::RubberParams rp;
   double hLayer = 0.0;
   RubberContact::RubberState st;
   std::vector<int> cand;   // a surface of faces: the faces the ball can touch this cycle (those of candV whose box is within the ball's radius plus its reach of the ball)
   std::vector<int> candV;  // the Verlet list: the faces within the reach plus a margin of where the list was built; cand is filtered from it every cycle
   RubberContact::V3 candC;       // where the list was built and how far it reaches: it stays valid while the ball moves less than candR - the reach asked for (a Verlet list)
   RubberContact::V3 tightC;      // where cand was filtered (it holds every face within the reach plus kTightSkin of there: valid while the ball stays within kTightSkin of it)
   double candR = -1.0;           // < 0 = no list
   static constexpr double kTightSkin = 0.002;   // m
   // the faces to test for a ball at c: the short list while the ball is within the skin of where it was filtered, the whole Verlet list beyond it
   // shared finder P2b (a static mesh host): the ball's per-cycle query's faces of this object (qFaces, set by ContactSurfaceHost::FinderBind; a re-query refreshes it in place)
   const std::vector<int> &Faces(const RubberContact::V3 &c) const { return qFaces != nullptr ? *qFaces : (c - tightC).Length() <= kTightSkin ? cand : candV; }
   // no face within reach: the surface has nothing to probe this cycle
   bool NoFaces(const RubberContact::V3 &c) const { return candR >= 0.0 && Faces(c).empty(); }
   bool touching = false;         // in contact at the last substep (a hit event is a transition to touching)
   SurfaceState gs;               // the group step's state of this surface; a moving part on the group has one per surface: the first here, the others in gsMore (sized between cycles, GroupSurfaces)
   std::vector<SurfaceState> gsMore;
   SurfaceState &Surf(const int s) { return s == 0 ? gs : gsMore[s - 1]; }
   const SurfaceState &Surf(const int s) const { return s == 0 ? gs : gsMore[s - 1]; }
   double endDepth = 0.0;         // the energy check: its depth at the end of the last cycle (m)
   double hcStep = 0.0, hcSpeed = 0.0, hcDelta = 0.0;   // the arrival-time loop's step for this contact and the speed and compression it was computed for (it is reused while they hold)
   double fnLast = 0.0;           // normal force of this contact at the last substep (N): the resting test compares the ball's leftover force with it
   bool hcRest = false;           // hcStep came from the resting rule   // 6-ball stall: hcStep came from the resting (quasi-static) rule
   bool partOut = false;          // island_reach: a moving part this ball cannot reach this cycle, solved in another island or alone: left out of this ball's solve (cleared at the cycle's end)
   bool partOutLast = false;      // it was left out last cycle (the next cycle's check: the ball must not be inside it)
   bool partTouching = false;     // a moving part's own release-gap hysteresis for its Hit event (HingePart::Part::Contact)
   // The last probe of this surface in the arrival loop: the gap (m, < 0 = overlap; probeT < 0 = none yet this cycle), where the ball was and when (s into the cycle). A static
   // surface's gap can shrink by no more than the ball has moved since, so a surface the ball cannot reach in the rest of the cycle is not probed again (ProbeBody).
   double probeGap = 0.0, probeT = -1.0;
   void *finderQ = nullptr;       // shared_finder: the ball's per-cycle query (contactsurfacehost.cpp FinderBallQ), set at each cycle's gather and at engagement (ContactSurfaceHost::FinderBind); null = none
   const std::vector<int> *qFaces = nullptr;   // shared finder P2b: this static mesh object's faces in the ball's query (Faces reads it in place of cand / candV); null = today's lists
   RubberContact::V3 probeC;
};

class ContactSurface
{
public:
   virtual ~ContactSurface() = default;
   virtual const std::string &Name() const = 0;
   // The faces (or the part of the geometry) the ball can reach this cycle: within ballR + reachM of c (metres)
   virtual void Refresh(SurfaceContact &it, double ballR, double reachM, const RubberContact::V3 &c) const = 0;
   // shared_finder: this surface's index in the finder's objects (a static mesh host in the index), -1 = not in it (its own finder)
   virtual int FinderObjIndex() const { return -1; }
   // The nearest point of the surface to a ball at c
   virtual RubberContact::TubeProbe Probe(const RubberContact::BallSI &bs, const SurfaceContact &it, const RubberContact::V3 &c) const = 0;
   // The same, as the group step's law reads it: the ball's contact patch has a size, so a fold under it pushes along the patch's average (rubbercontact.h PatchTube, PatchTriangles)
   virtual RubberContact::TubeProbe ProbePatch(const RubberContact::BallSI &bs, const SurfaceContact &it, const RubberContact::V3 &c) const { return Probe(bs, it, c); }
   // The longest step, at most hMax, that brings the ball at most to arriveM of a part of the surface other than the probe's nearest (nearestSeg); a nearest-only surface: hMax. a: the acceleration bound.
   virtual double ArrivalS(const RubberContact::BallSI & /*bs*/, const SurfaceContact & /*it*/, int /*nearestSeg*/, const RubberContact::V3 & /*c*/, const RubberContact::V3 & /*v*/, double /*a*/, double hMax, double /*touchM*/, double /*arriveM*/) const { return hMax; }
   virtual const RubberObjectSettings *Settings() const { return nullptr; }   // the step and rest settings of a surface that has them (rubber, wall, primitive, ramp, playfield, disc); null = none
   virtual bool IsFloor() const = 0;        // a surface a ball rests on: gravity acts in the contact's substeps while one is in the contact
   virtual bool GroupModel() const { return false; }   // a ball on it goes through the group step ([Rubbers] ring_model, GroupStepDesign 14)
   virtual bool IsPlayfield() const = 0;
   virtual double ReleaseGapM() const = 0;  // a ball this close to the surface stays engaged (hysteresis)
   // The fastest the surface itself can move toward a ball within a cycle (m/s, ContactMovingPart::SurfaceSpeedBound); 0 for a surface that does not move
   virtual double SurfaceSpeedBoundMS(double /*cycleS*/) const { return 0.0; }
   // The part object whose script hit event, threshold and duplicate filter a first touch of this surface fires (null = none)
   virtual HitObject *HitEventObject() const { return nullptr; }
   // A part whose motion the contact system integrates (contactmovingpart.h: a gate, a spinner, a jet bumper). It fires its own events.
   virtual ContactMovingPart *MovingPart() const { return nullptr; }
   virtual bool HitEventWanted(const RubberContact::V3 & /*normalToBall*/) const { return true; }
   // A surface that moves: its velocity and angular velocity at a probed contact point (default: at rest), and the reaction of one substep's contact on it
   virtual void VelocityAt(const RubberContact::TubeProbe &, RubberContact::V3 &, RubberContact::V3 &) const {}
   virtual void ApplyReaction(const RubberContact::TubeProbe &, const RubberContact::ContactResult &, double /*ballR*/, double /*h*/) const {}
};
