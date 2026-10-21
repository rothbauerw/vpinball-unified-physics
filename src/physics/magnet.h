#pragma once
#include <vector>
class HitBall;
// Electromagnets (coil + iron core) that pull the steel ball. Configured in [Magnet:<name>] sections of vpx_physics.ini;
// the table script switches a coil with Magnet "<name>", volts and reads its current with MagnetCurrent("<name>").
//
// Physics (SI inside, converted at the ball boundary):
//  - coil: N turns of round wire wound on the core -> resistance R from the wire length; inductance L = N^2 / (unloaded R).
//    I follows dI/dt = (V - R I) / L (exact exponential per 1 ms step); with V = 0 it decays through R (freewheel diode).
//  - flux: Phi = N I / R, R = far pole (1 - return_fraction) / (8 mu0 a) + near pole 1 / (8 mu0 a + mu0 A / gap)  (a = core radius, gap = face to ball surface); Bc = min(Bsat, Phi / A).
//  - field: the near pole face is a uniformly charged disc (24 point charges, surface density Bc); the opposite pole is one charge behind the coil (bare rod),
//    or, with return_pole_distance_mm, point charges on the yoke flanges in a plane parallel to the face. B and its gradient are the point-charge sums: finite at the face.
//  - ball: induced dipole of a steel sphere, M = min(3 |B| / mu0, Ms); force F = V_ball M grad|B| (along the field gradient).
//  - acceleration = F / (steel mass of the ball * Ball.Mass).
// Limitations: the flux sees the ball only through the gap term; the force is evaluated at the ball centre; one iron core, no saturation of the yoke.
// Keys: see the [Magnet:<name>] section of vpx_physics.ini.

namespace Magnets
{
   void Step();                                              // once per 1 ms physics frame: advance every coil current
   void RenderGuides(class Renderer *renderer);              // draws the guide tubes and stops that have tube_visible = 1
   void ReleaseRender();                                     // frees the drawn guide meshes (called when the renderer shuts down)
   bool HasGuides();                                         // true if an air-core magnet has a guide tube or stop (checked once per physics cycle)
   // air-core guide tube and stop: keeps the ball (pos U, vel U/VPT, angular momentum, radius U, mass, inertia) inside the tube and short of the stop; true if it moved the ball
   bool ConstrainBall(float pos[3], float vel[3], float angmom[3], float radiusU, float mass, float inertia);
   bool Active();                                            // true if at least one coil carries current (the force loop is skipped otherwise)
   // Option A ball-ball force (build step B7): solved once per physics step, re-aimed at every substep by the caller. One part of a ball's force (N, table frame, m):
   //   b != nullptr: a PAIR part: turn F with the line from b to a (u0 = that line's unit vector at the step start; pa0 / pb0 = the centres then);
   //   b == nullptr: a MAGNET part: turn F round the magnet axis n (o = face centre) as ball a moves round it (r0 = a's unit radial direction at the step start).
   struct BallBallPart { const HitBall *a = nullptr, *b = nullptr; double F[3] = { 0, 0, 0 }, u0[3] = { 0, 0, 1 }, n[3] = { 0, 0, 1 }, o[3] = { 0, 0, 0 }, r0[3] = { 1, 0, 0 }, pa0[3] = { 0, 0, 0 }, pb0[3] = { 0, 0, 0 }; };
   bool HasBallBallParts(const HitBall *hb);
   // the parts of hb this step (nullptr: none); massKg = the steel mass the magnets use for hb (7800 kg/m3 x volume x Ball.Mass): F x (body mass / massKg) gives
   // the same acceleration as the held increment (m_magPair)
   const std::vector<BallBallPart> *BallBallParts(const HitBall *hb, double &massKg);
   // adds the velocity change over one physics step (VPU/VPT) of a ball at pos (U) with radius (U) and mass (VPX mass units) to dv
   void AddVelocityIncrement(const float pos[3], float radius, float mass, float dv[3]);
   // field (T) inside the steel ball at pos (U), radius (U): the coil field times up to 3 (a permeable sphere gathers the field), for the eddy currents of a spinning ball; false if no coil carries current
   bool FieldAt(const float pos[3], float radiusU, double B[3]);
   // electrical conductivity of the ball steel (S/m) times the ball_eddy scale: 0 = no eddy-current braking of the ball's spin
   double EddyConductivity();
   bool SetVoltage(const char *name, float volts);           // false = no such magnet
   bool SetDuty(const char *name, float duty);               // PWM: mean voltage = duty x rated_volts (-1..1)
   bool GetCurrent(const char *name, float &amps);
}
