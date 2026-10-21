// license:GPLv3+

#pragma once

#include "math/math.h"
#include "math/matrix.h"
#include "physics/collide.h"
#include "physics/contactsettings.h"


// Fake/artistic ball motion trail length, but also used for spin kill, so may not be lower than 10 (to keep 100ms of ball history at 10ms interval)
#define MAX_BALL_TRAIL_POS 10


class FlipperMoverObject;

class BallMoverObject : public MoverObject
{
public:
   bool AddToList() const override { return false; } // We add ourselves to the mover list.
                                                     // If we allow the table to do that, we might get added twice, if we get created in the player Init code
   void UpdateDisplacements(const float dtime) override;
   void UpdateVelocities() override;

   HitBall* m_pHitBall;
};

struct BallS
{
   Vertex3Ds m_pos = Vertex3Ds(0.f, 0.f, 0.f);
   Vertex3Ds m_vel = Vertex3Ds(0.f, 0.f, 0.f); // ball velocity
   float m_radius = DEFAULT_BALL_SIZE;
   float m_mass   = 1.f;
   vector<IFireEvents*>* m_vpVolObjs; // vector of triggers and kickers we are now inside (stored as IFireEvents* though, as HitObject.m_obj stores it like that!)
   bool m_lockedInKicker = false;
   class HitBall *m_owner = nullptr;   // the ball this state belongs to (set by HitBall; lets a check with only the state find the ball's contact in O(1))
};

class HitBall final : public HitObject
{
public:
   HitBall();
   ~HitBall() override;

   void UpdateDisplacements(const float dtime);
   void UpdateVelocities();
   Vertex3Ds ExternalVelocityIncrement(bool withMagnets = true) const;   // gravity, nudge and (optionally) the magnets over one physics step, U/VPT (hitballcontact.cpp)

   // From HitObject
   float HitTest(const BallS& ball, const float dtime, CollisionEvent& coll) const override;
   int GetType() const override { return eBall; }
   void Collide(const CollisionEvent& coll) override;
   void Contact(CollisionEvent& coll, const float dtime) override; // two body contact handling for coll.m_ball against this ball
   void CalcHitBBox() override;

   float HitRadiusSqr() const { return sqrf((m_hitBBox.right - m_hitBBox.left)*0.5f); } // this returns the extended (by m_vel + magic) squared radius, as needed to be used in the collision detection
   void Collide3DWall(const Vertex3Ds& hitNormal, float elasticity, const float elastFalloff, const float friction, float scatter_angle);

   // other: the contacted ball for ball/ball contacts (relative approach; relative friction only with FIX_PHYSICS), nullptr for static surfaces
   void ApplyFriction(const Vertex3Ds& hitnormal, const float dtime, const float fricCoeff,
      const float normalImpulse, const HitBall* const other = nullptr); // normalImpulse: normal delta-v applied by the contact this step (only used with FIX_PHYSICS)
   void HandleStaticContact(const CollisionEvent& coll, const float friction, const float dtime, const HitBall* const other = nullptr);
#ifdef NEW_PHYSICS
   void ApplyRotationalFriction(const Vertex3Ds& hitnormal, const float fricCoeff, const float normalImpulse, const HitBall* const other = nullptr); // as for ApplyFriction
#endif

   // Rolling resistance and twist friction of a resting contact (unified engine, ContactSettings): spin about the normal and in the tangent plane each lose up to what the
   // normal impulse can resist, never reversing (hitballcontact.cpp)
   void ApplyRestingSpinLoss(const Vertex3Ds& hitnormal, const float fricCoeff, const float normalImpulse);
   // A position push (overlap correction) that does not shove a ball down into a surface the contact model holds it on (hitballcontact.cpp)
   void Shove(Vertex3Ds delta);

   Vertex3Ds SurfaceVelocity(const Vertex3Ds& surfP) const;
   Vertex3Ds SurfaceAcceleration(const Vertex3Ds& surfP) const;
   float Inertia() const { return (float)(2.0/5.0) * m_d.m_radius*m_d.m_radius * m_d.m_mass; }

   void ApplySurfaceImpulse(const Vertex3Ds& rotI, const Vertex3Ds& impulse);

   void DrawUI(std::function<Vertex2D(Vertex3Ds)> project, ImDrawList* drawList, bool fill) const override;

   // Per frame info
   Ball *m_pBall = nullptr;  // Object model version of the ball

   CollisionEvent m_coll;    // collision information, may not be a actual hit if something else happens first

   BallMoverObject m_mover;

   BallS m_d;

   Vertex3Ds m_oldVel;       // hack for kicker hole handling only

   bool m_rubberAdvanced = false;               // the rubber model already moved this ball for the current cycle (UpdateDisplacements skips the straight-line move)
   class HitRubber *m_bandOwner = nullptr;      // band rubber whose engine-stepped contact currently drives this ball (null = none)
   Vertex3Ds m_nearPos;                                 // proximity engagement: where the registry was last scanned for this ball, and how far it may move before the next scan (VPX units, the nearest surface box less a margin); < 0 = scan
   float m_nearSafeU = -1.f;
   uint32_t m_nearMs = 0;
   struct BallContact *m_contact = nullptr;             // the contact system's engagement that holds this ball (a post, a wall rubber, the playfield, a disc); null = none
   FlipperMoverObject *m_rubberOwner = nullptr; // model flipper whose rubber contact model currently drives this ball (null = none)
   bool m_magnetic = true;                      // Ball.Magnetic: false = non-magnetic (ceramic): no pull from the magnets, no field of its own
   float m_magPair[3] = { 0.f, 0.f, 0.f };      // magnets: velocity change this physics step from the other balls' induced dipoles (VPX units), set by Magnets::Step
   float m_hitSpeed = 0.f;   // speed along the surface normal just before the latest impact (same units as the velocity), for scripts (Ball.HitSpeed) and the playfield hit event

   Vertex3Ds m_lastEventPos; // last hit event position (to filter hit 'equal' hit events)
   float m_lastEventSqrDist = 0.f; // distance travelled since last event

   Vertex3Ds m_angularmomentum;

   Matrix3 m_orientation;

#ifdef FIX_PHYSICS
   // velocity change the contacts resolved so far in the current contact pass applied to the ball (reset by PhysicsEngine):
   // all of it for ball/ball contacts, only the tangential (friction) part for static contacts. Static contacts add it to
   // the velocity they saw at hit test time (see HandleStaticContact)
   Vertex3Ds m_contactDeltaV = Vertex3Ds(0.f, 0.f, 0.f);
#endif

#ifdef C_DYNAMIC
   int m_dynamic = C_DYNAMIC; // used to determine static ball conditions and velocity quenching
   float m_drsq = 0.0f;       // square of distance moved
#endif

   void OnPhysicStepProcessed(uint64_t physicsTimeUs);
   const Vertex3Ds& GetOldPosition(uint64_t physicsTimeUs) const;

private:
   Vertex3Ds m_oldpos[MAX_BALL_TRAIL_POS]; // used for killing spin (see C_BALL_SPIN_HACK) and for ball trails
};

// Unified engine (hitrubber.cpp)
bool HitRubberContactPairModel(const HitBall *a, const BallS &b); // both balls are contact-held and ball contact is on
bool HitRubberContactHolds(const HitBall *ball);                  // a force-model contact (a post, the playfield) is stepping this ball
