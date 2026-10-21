#pragma once
// A surface a ball can touch (docs/UnifiedContactDesign.txt, "Surface"): static geometry or a driven one (the turntable's disc). The contact system asks it for the nearest
// point, normal and curvature of the part of the surface near the ball (Probe), and, for a surface that moves, the velocity of the surface at the contact and the reaction
// on it. The contact law (Hertz + damping + Mindlin + twist + rolling resistance) and the integration are the same for every surface and live outside.
//
// Today's implementation is HitRubber (a rubber tube, the faces of a wall or primitive, the playfield, a turntable disc). The unified engine's own surfaces
// (playfield mesh, ramps, wall and primitive meshes) will implement this interface directly.

#include <string>
#include <vector>
#include "physics/rubbercontact.h"

class ContactSurface;
class HitObject;
namespace HingePart { class Part; }
struct RubberObjectSettings;

// What one ball has with ONE surface: the material parameters, the shear state kept between cycles, and the faces near the ball this cycle.
struct SurfaceContact
{
   const ContactSurface *r = nullptr;
   RubberContact::RubberParams rp;
   double hLayer = 0.0;
   RubberContact::RubberState st;
   std::vector<int> cand;   // a surface of faces: the faces near the ball this cycle
   RubberContact::V3 candC;       // where the list was built and how far it reaches: it stays valid while the ball moves less than candR - the reach asked for (a Verlet list)
   double candR = -1.0;           // < 0 = no list
   bool touching = false;         // in contact at the last substep (a hit event is a transition to touching)
   double hcStep = 0.0, hcSpeed = 0.0, hcDelta = 0.0;   // the arrival-time loop's step for this contact and the speed and compression it was computed for (it is reused while they hold)
   double fnLast = 0.0;           // normal force of this contact at the last substep (N): the resting test compares the ball's leftover force with it
   bool hcRest = false;           // hcStep came from the resting rule   // 6-ball stall: hcStep came from the resting (quasi-static) rule
   bool partTouching = false;     // a moving part's own release-gap hysteresis for its Hit event (HingePart::Part::Contact)
};

class ContactSurface
{
public:
   virtual ~ContactSurface() = default;
   virtual const std::string &Name() const = 0;
   // The faces (or the part of the geometry) the ball can reach this cycle: within ballR + reachM of c (metres)
   virtual void Refresh(SurfaceContact &it, double ballR, double reachM, const RubberContact::V3 &c) const = 0;
   // The nearest point of the surface to a ball at c
   virtual RubberContact::TubeProbe Probe(const RubberContact::BallSI &bs, const SurfaceContact &it, const RubberContact::V3 &c) const = 0;
   virtual const RubberObjectSettings *Settings() const { return nullptr; }   // the step and rest settings of a surface that has them (rubber, wall, primitive, ramp, playfield, disc); null = none
   virtual bool IsFloor() const = 0;        // a surface a ball rests on: gravity acts in the contact's substeps while one is in the contact
   virtual bool IsPlayfield() const = 0;
   virtual double ReleaseGapM() const = 0;  // a ball this close to the surface stays engaged (hysteresis)
   // The part object whose script hit event, threshold and duplicate filter a first touch of this surface fires (null = none)
   virtual HitObject *HitEventObject() const { return nullptr; }
   // A part on one axis whose motion the contact system integrates (a gate, a spinner: contacthingepart.h). It evaluates its own contact with a ball (EvaluateOwn: the series bend,
   // the tilt force, its Hit event, the reaction on its body) and fires its own events.
   virtual HingePart::Part *MovingPart() const { return nullptr; }
   virtual RubberContact::ContactResult EvaluateOwn(SurfaceContact & /*it*/, const RubberContact::BallSI & /*bs*/, const RubberContact::V3 & /*c*/, const RubberContact::V3 & /*v*/,
                                                    const RubberContact::V3 & /*w*/, double /*h*/) const { return RubberContact::ContactResult(); }
   virtual bool HitEventWanted(const RubberContact::V3 & /*normalToBall*/) const { return true; }
   // A surface that moves: its velocity and angular velocity at a probed contact point (default: at rest), and the reaction of one substep's contact on it
   virtual void VelocityAt(const RubberContact::TubeProbe &, RubberContact::V3 &, RubberContact::V3 &) const {}
   virtual void ApplyReaction(const RubberContact::TubeProbe &, const RubberContact::ContactResult &, double /*ballR*/, double /*h*/) const {}
};
