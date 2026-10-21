// license:GPLv3+

#pragma once

#include "collide.h"
#include "physics/rubbercontact.h"

//#define DEBUG_FLIPPERS

class Flipper;

#include <filesystem>
#include <memory>
#include <string>

class HitRubber;
struct RubberSideSet;

#include "physics/coilmodel.h"
#include "physics/contactbat.h"

class FlipperMoverObject : public MoverObject
{
public:
   FlipperMoverObject(const Vertex2D& center, const float baser, const float endr, const float flipr, const float angleStart, float angleEnd,
      const float zlow, const float zhigh, Flipper* const pflipper);

   void UpdateDisplacements(const float dtime) override;
   void UpdateVelocities() override;

   bool AddToList() const override { return true; }

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
   // model on there is no fixed limit (its stroke is stop_deg and it overshoots by a load dependent amount), so no limit.
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

   // Two-body flipper model (unified engine): a plunger body and a
   // flipper (bat + crank + shaft, one rigid body) coupled through the link's pin slop, with a
   // power coil, a hold coil and the return spring. Replaces the stock torque model; table
   // strength/return/mass/sweep are ignored (rothbauerw's real-data values are used instead).
   void FlipperEnsureInit();
   void FlipperStep(double h, double extTorqueNm = 0.0);   // extTorqueNm: external torque on the flipper in the closing direction (rubber contact with a ball), 0 = none
   void FlipperUpdateInertia();
   void UpdateDisplacementsFlipper(const float dtime);
   void FlipperReloadSettings(bool force);

   // 3D ball/rubber geometry test (physics/rubbercontact.h): does the ball overlap the rubber surface right now?
   // Allows overlap, so it also works for a ball that is already engaged. Pure geometry, no state change.
   RubberContact::Probe RubberProbe(const BallS& ball) const;

   // Engagement of a ball with this flipper's rubber (one ball per flipper, one flipper per ball). While engaged the
   // ball is driven by the rubber force model and the impulse path is skipped.
   static bool MultiBall();
   bool RubberWants(const HitBall *ball) const;   // rubber_contact is on and this ball may engage (else the impulse path runs)
   bool RubberIsEngaged(const BallS& ball) const;
   bool RubberSideHas(const HitRubber *r) const;   // the force-model rubbers around the engaged ball (a post beside the bat) are part of this flipper's step
   std::shared_ptr<RubberSideSet> m_rbSide;        // those rubbers with their shear state, null = none engaged
   void RubberEngage(HitBall *ball);
   void RubberRelease();
   // One substep of release bookkeeping: pass whether the ball is clear of the rubber and moving away. Returns true
   // (and releases) once that has held for rubber_release_steps substeps in a row.
   bool RubberReleaseStep(bool clearAndReceding);
   HitBall *m_rbBall = nullptr;   // engaged ball, null = none
   int m_rbClearSteps = 0;        // consecutive substeps clear and receding
   RubberContact::RubberState m_rbState;   // shear spring state of the current engagement
   double m_rbLastDist = 0.0;     // ball centre to rubber surface distance at the previous substep, m
   bool m_rbJustEngaged = false;  // the next RubberAdvance is the first of this engagement (checks and seats the starting overlap)
   uint32_t m_rbCooldownUntilMs = 0;   // after the safety net trips, no new engagement before this game time (the impulse path runs meanwhile)
   int m_rbLogSteps = 0;          // log only: substeps in contact since the last log row
   double m_rbLogPeakFn = 0.0, m_rbLogPeakDelta = 0.0;   // log only: peak normal force (N) and compression (m) since the last row
   void RubberLogRow();           // VPX_RUBBER_LOG: one row per physics step (LeftFlipper only), engaged or not
   bool m_gAdvanced = false;      // RubberAdvance already moved bat and ball for this cycle; UpdateDisplacements only finishes up
   bool m_gAdvHitStop = false;    // plunger reached the coil stopper during RubberAdvance

   // Moves the bat and its engaged ball together in substeps (rubber force on the ball, opposite torque on the bat).
   void PreUpdateDisplacements(const float dtime) override;
   bool RubberAdvance(const float dtime);
   ~FlipperMoverObject();

   // Flipper nudge. The striking flipper (settings nudge_target / nudge_gain) starts a one-shot pulse on the target flipper's pivot: x displacement
   // away from the striker, half-cosine rise then half-cosine return. The pivot (m_hitcircleBase.center.x) follows it, ball contacts on the bat
   // see the pivot velocity (rubber path in RubberAdvance, impulse path through SurfaceVelocity).
   void NudgeStart(double dispM, double riseS, double fallS);
   void NudgeEval(double t, double& dispM, double& velMs) const;
   void NudgeStep(double dtS);
   FlipperMoverObject *NudgeFindTarget();
   std::string m_nzTargetName;
   FlipperMoverObject *m_nzTarget = nullptr;

   ContactBat m_bat;                 // the bat body: settings, plunger/coil/spring state, the substep and the pivot nudge (contactbat.h)

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

// The two-body flipper model is on (unified engine; hitflippermodel.cpp)
bool FlipperModelEnabled();
