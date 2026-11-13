#pragma once
#include <vector>
#include "physics/pairform.h"
class HitBall;
namespace MagnetCore { struct Model; }
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
   void StartPlay();                                         // table load (Player::InitTableSession, after the script's Init, before the first physics step): reads the magnets from the physics ini
   void Step();                                              // once per 1 ms physics frame: advance every coil current (never reads the file)
   bool CurrentSteady();                                     // pile sleep: no coil current would change in the next Step (each at its own fixed point)
   void Reload();                                            // live tuning: the physics ini changed (ContactSystem::ReloadChangedIni, on playfield focus gain): the magnets read again, coils keep volts and amps
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
   // The ball-pair force as a stored energy (physics/pairform.h, 2026-11-03), fitted at this step's Option A solve; the group step takes it instead of the held m_magPair.
   // a, b: the balls of a term (b == nullptr: a one-ball term); t.a and t.b are the solve's own indices, unused by the caller. Empty: no form this step.
   struct FormTerm { const HitBall *a = nullptr, *b = nullptr; PairForm::Term t; };
   const std::vector<FormTerm> &PairFormTerms();
   bool HasPairForm(const HitBall *hb);
   // adds the velocity change over one physics step (VPU/VPT) of a ball at pos (U) with radius (U) and mass (VPX mass units) to dv. heldOnly: the magnets whose pull the group step holds for the cycle (not core_model = cells: CellPulls)
   // aircoreInStep ([Features] aircore_in_step, with heldOnly): the air-core coils are left out too, the group step reads them at the step's midpoint (AircoreForce)
   void AddVelocityIncrement(const float pos[3], float radius, float mass, float dv[3], bool heldOnly = false, bool aircoreInStep = false);
   // The air-core coils' pull on a solid steel ball (7800 kg/m3) of radius radiusU (U) centred at pos (U), summed over the coils with current: F (N, table frame) and, when S is
   // given, its slope S[i][j] = dF_i / dpos_j per metre of the centre's move. The law of AddVelocityIncrement's air-core branch (one shared function); a ball of Ball.Mass m
   // and steel mass M_s = 7800 V m feels F x (its mass / M_s), as there.
   void AircoreForce(const double pos[3], double radiusU, double F[3], double S[3][3]);
   // A cells magnet's pull as the group step takes it, a stored energy of the ball's place (magnetcore.h Model::Energy): its tables, its pole face (m), its axis, the coil's current held for this physics step and the one
   // before it (a change between them is the coil's work). Every cells magnet with a current in either. StepNumber: counts the physics steps.
   struct CellPull { const MagnetCore::Model *core = nullptr; double o[3] = { 0, 0, 0 }, n[3] = { 0, 0, 1 }, amps = 0, ampsBefore = 0; };
   const std::vector<CellPull> &CellPulls();
   long long StepNumber();
   // field (T) inside the steel ball at pos (U), radius (U): the coil field times up to 3 (a permeable sphere gathers the field), for the eddy currents of a spinning ball; false if no coil carries current
   bool FieldAt(const float pos[3], float radiusU, double B[3]);
   // electrical conductivity of the ball steel (S/m) times the ball_eddy scale: 0 = no eddy-current braking of the ball's spin
   double EddyConductivity();
   bool SetVoltage(const char *name, float volts);           // false = no such magnet
   bool SetDuty(const char *name, float duty);               // PWM: mean voltage = duty x rated_volts (-1..1)
   bool GetCurrent(const char *name, float &amps);
}
