// license:GPLv3+

#pragma once

#include "collide.h"
#include "physics/rubbercontact.h"

//#define DEBUG_FLIPPERS

class Flipper;

#include <filesystem>
#include <memory>
#include <string>

class ContactSurfaceHost;
class ContactHost;
class ContactPartHost;
class HitFlipper;
namespace BatPart { class Part; }

#include "physics/coilmodel.h"
#include "physics/contactbat.h"

class FlipperMoverObject : public MoverObject
{
public:
   FlipperMoverObject(const Vertex2D& center, const float baser, const float endr, const float flipr, const float angleStart, float angleEnd,
      const float zlow, const float zhigh, Flipper* const pflipper);

   void UpdateDisplacements(const float dtime) override;
   void UpdateVelocities() override;

   bool AddToList() const override { return m_part == nullptr; }   // (a bat that is a moving part is advanced by its host, ContactPartHost)

   void SetSolenoidState(const bool s);
   float GetStrokeRatio() const;

#ifdef FIX_PHYSICS
   // Recomputes m_zeroAngNorm, m_inertiaShape, m_torqueScale and m_inertia
   // from current geometry (m_hitcircleBase.radius, m_endradius, m_flipperradius)
   // and current mass (m_d.m_OverrideMass or m_d.m_mass).
   // Called from the constructor and from HitFlipper::UpdatePhysicsFromFlipper
   // so live edits to mass propagate. Geometry edits via the COM API today
   // do not propagate to the runtime mover, but this helper is ready for the wiring
   void UpdateInertia();
#endif
   // Angle range (rad) the ball collision tests may assume the bat reaches. Stock: the table's Start/End angle. With the Flipper
   // model on the bat has no collider to test (HitTest), so this is the stock flipper's.
   void HitTestAngleLimits(float& angleMin, float& angleMax) const;

   void SetStartAngle(const float r);
   void SetEndAngle(const float r);
   float GetReturnRatio() const;
   float GetMass() const;
   void SetMass(const float m);
   float GetStrength() const;

   // rigid body functions
   Vertex3Ds SurfaceVelocity(const Vertex3Ds& surfP) const;
   Vertex3Ds SurfaceAcceleration(const Vertex3Ds& surfP) const;

   float GetHitTime() const;

   void ApplyImpulse(const Vertex3Ds& rotI);


   ~FlipperMoverObject();

   // Flipper nudge. The striking flipper (settings nudge_target / nudge_gain) starts a one-shot pulse on the target flipper's pivot: x displacement
   // away from the striker, half-cosine rise then half-cosine return. The pivot (m_hitcircleBase.center.x) follows it, ball contacts on the bat
   // see the pivot velocity (BatGroup::Part).
   void NudgeStart(double dispM, double riseS, double fallS);
   FlipperMoverObject *NudgeFindTarget();
   std::string m_nzTargetName;
   FlipperMoverObject *m_nzTarget = nullptr;


   // Unified engine: the bat is a moving part of the contact system, on the group step (contactbatpart.h). This mover and the stock collider then only hold the state the script and the
   // renderer read (angle, speed, solenoid); the part's adapter writes it (hitflippermodel.cpp). There is no other path for the bat on the unified engine.
   std::shared_ptr<BatPart::Part> m_part;
   // Builds the part and its host. zlow, zhigh: the bat's height range (VPX units).
   ContactPartHost *MakePart(HitFlipper *owner, float zlow, float zhigh);
   void SetPartSolenoid(bool s);   // the coil input of the part (the table set it between cycles)

   Flipper *m_pflipper;

   //float m_faceLength;
   HitCircle m_hitcircleBase;
   float m_endradius;
   float m_flipperradius;

   // kinematic state
   float m_angularMomentum;
   float m_angularAcceleration;
   float m_angleSpeed;
   float m_angleCur;

   float m_curTorque;
   float m_contactTorque;

   float m_angleStart, m_angleEnd;

   float m_inertia;         // moment of inertia
#ifdef FIX_PHYSICS
   float m_inertiaShape;    // shape-only factor k^2 (length^2), so m_inertia = mass * m_inertiaShape
   float m_torqueScale;     // multiplier applied to coil torque; 1.0 unless m_backwards_compatibility is set
#endif

   Vertex2D m_zeroAngNorm;  // base norms at zero degrees

   short m_enableRotateEvent; // -1,0,1

   bool m_direction;

   bool m_solState;         // is solenoid enabled?
   bool m_isInContact;

   bool m_enabled;
   bool m_lastHitFace;

#ifdef FIX_PHYSICS
   // When true (legacy model/default), coil torque is scaled by
   // m_inertiaShape / (flipr^2 / 3) so that the per-Strength angular
   // acceleration matches the old rod-approximation inertia and existing
   // tables keep their feel. Set to false to drive the flipper with the
   // raw Strength value against the more accurate inertia
   bool m_backwards_compatibility;
#endif

#ifdef DEBUG_FLIPPERS
   uint32_t m_startTime;
#endif
};

class HitFlipper : public HitObject
{
public:
   HitFlipper(const Vertex2D& center, const float baser, const float endr, const float flipr, const float angleStart, const float angleEnd,
              const float zlow, const float zhigh, Flipper* const pflipper);
   ~HitFlipper() override { /*m_pflipper->m_phitflipper = nullptr;*/ }

   float HitTest(const BallS& ball, const float dtime, CollisionEvent& coll) const override;
   int GetType() const override { return eFlipper; }
   void Collide(const CollisionEvent& coll) override;
   void Contact(CollisionEvent& coll, const float dtime) override;
   void CalcHitBBox() override;
   MoverObject *GetMoverObject() override { return &m_flipperMover; }

   void UpdatePhysicsFromFlipper();

   float HitTestFlipperFace(const BallS& ball, const float dtime, CollisionEvent& coll, const bool face1) const;
   float HitTestFlipperEnd(const BallS& ball, const float dtime, CollisionEvent& coll) const;

   float GetHitTime() const { return m_flipperMover.GetHitTime(); }
   float ReachToBall(const BallS& ball, const float horizon) const; // how far (VPX units) the bat point nearest the ball travels in `horizon` physics steps at the bat's current angular speed (little near the base, most at the tip)
   float GapToBall(const BallS& ball) const; // distance (VPX units, plan view) from the ball surface to the bat outline at its current angle; negative = overlapping

   void DrawUI(std::function<Vertex2D(Vertex3Ds)> project, ImDrawList* drawList, bool fill) const override;

   FlipperMoverObject m_flipperMover;

   // the same from scalars: bnv the normal speed (engine units, < 0 approaching), hitMomentBit as CollisionEvent::m_hitmoment_bit, arrival the ball's own speed toward the face
   // (the bat's part adapter fires its Hit through this, hitflippermodel.cpp)
   void FireHitEventsAndRumble(float bnv, bool hitMomentBit, float arrival);
private:
   void FireHitEventsAndRumble(float bnv, const CollisionEvent& coll, const Vertex3Ds& normal, const Vertex3Ds& vB); // script Hit/Collide event and haptics, shared by the impulse and rubber paths

   uint32_t m_last_hittime;
   // Contacts of the last RUMBLE_WINDOW_MS, summed for the contact rumble, see Collide()
   static constexpr uint32_t RUMBLE_WINDOW_MS = 80;
   static constexpr int RUMBLE_CONTACTS = 32; // more than a slap delivers within a window
   static constexpr float RUMBLE_FULL_IMPACT = 30.f; // the sum at which the contact rumble saturates in level and length, see InputManager::PlayFlipperContactRumble
   uint32_t m_rumbleContactMs[RUMBLE_CONTACTS];
   float m_rumbleContactImpact[RUMBLE_CONTACTS];
   int m_rumbleContactTail = 0;
   int m_rumbleContactCount = 0;
   float m_rumblePlayed = 0.f; // sum last played, reset when the window runs empty
};

