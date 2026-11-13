// license:GPLv3+

#pragma once

class HitBall;
class Kicker;
class Vertex3Ds;

// The physical kicker's arm stroke (docs/unified-physics/KickerDesign.txt 4c, approved by the user 2026-11-13; code in parts/kicker.cpp). A Kick on a physical kicker hands the held
// ball to an idealised arm: its centre is carried along a quarter circle in the vertical plane of the kick direction, from where it rests (moving straight up) to where the path's direction
// equals the kick's inclination, at the kick's speed, then let go with that velocity and no spin. The arm's work is booked when the kick starts (ContactSystem::NoteDriveWork).
// The carried ball is a held ball (KickerHold, design 1h): not stepped, on the arm's path; another ball meets it through the shared pair law as a body of infinite mass.
namespace KickerArm
{
   void Start(Kicker *kicker, HitBall *ball, double dirX, double dirY, double radiusU, double releaseRad, double speed);   // Kicker::KickXYZ: dir unit (VPX x, y), speed U/VPT
   bool Carried(const HitBall *ball);                // the arm carries this ball (a second Kick does nothing while it does)
   void BallGone(const HitBall *ball);               // Kicker::BallGone: a destroyed ball is no longer carried
   void Remove(const Kicker *kicker);                // Kicker::PhysicSetup / PhysicRelease: this kicker's strokes end (table start and end)
}

// Held balls (KickerDesign.txt 1h, user 2026-11-13; Requirements F 10): a ball a kicker holds still (a virtual kicker's parked ball) or carries (a physical kicker's arm) has infinite mass.
// It is not stepped: no gravity, no surface contacts, no stock hit tests, no magnets, spin zero. Its place is fixed (parked) or the arm's path at each moment (carried), and its velocity is
// the given one (0, or the path's). In the shared ball-ball pair law (GroupStep, contactsystem.cpp) it is a body forces do not move (GroupStep::Ball::held): another ball meets it with the
// same Hertz + damping + friction, with the held ball's velocity in the relative velocity. A parked ball does no work; a carried ball's work on another ball is the arm's drive work
// (the group's books.driveJ). The capture (Kicker::Park) zeroes speed and spin once, a loss booked between cycles (1i); nothing is put back afterwards.
namespace KickerHold
{
   void Add(Kicker *kicker);                         // Kicker::PhysicSetup: a virtual kicker (it may hold a parked ball)
   void Remove(const Kicker *kicker);                // Kicker::PhysicRelease
   void CycleEnd(double cycleS);                     // PhysicsEngine::PhysicsSimulateCycle, unified engine, before the cycle's events: each carried ball on the arm's path at the cycle's
                                                     // end time (let go at the release point); each parked ball checked at its place at rest (a fault said once if not: nothing should move it)
   bool Held(const HitBall *ball, Vertex3Ds *accelU = nullptr);   // a ball a kicker holds (parked) or carries. accelU: its given acceleration, U/VPT^2 (the arm's path's at this
                                                     // moment, the centripetal pull of the circle; 0 parked); its given velocity is the ball's own m_d.m_vel (the kicker sets it)
}
