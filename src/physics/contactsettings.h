// license:GPLv3+
#pragma once

#include <cmath>
#include <filesystem>

// Table-wide contact settings of the unified engine (see docs/NewPhysicsReference.txt, [Contact]); implemented in contactsettings.cpp.

// Table-wide contact friction settings beyond the per-object friction coefficient (the surface's static mu). Read from [Contact] of the table's
// vpx_physics.ini and reloaded when it changes. Without a file these defaults apply (off on the stock engine).
struct ContactSettings
{
   bool enable = true;               // false = plain Coulomb friction with one mu (the object's friction), no rolling resistance, no twist
   double muDynamicRatio = 0.85;     // dynamic (sliding) mu = ratio * the object's friction, which is the static mu (typically 0.7 - 0.95 for hard dry surfaces)
   double rollingResistance = 0.003; // rolling resistance coefficient (dimensionless, not a friction mu): torque = rollingResistance * R * N (steel ball on hard wood / plastic is ~0.002 - 0.005)
   bool twist = true;                // drilling friction: spin about the contact normal is resisted
   double twistPatchMm = 0.5;        // contact patch radius used for the twist torque (3 pi / 16) * mu * a * N; steel on hard wood ~ 0.3 - 0.6 mm
   // Ball-to-ball contact in the force model (two contact-held balls are stepped together; the stock impulse and position push are skipped for the pair)
   bool ballContact = false;         // ball_contact
   double ballEMpa = 200000.0;       // ball_e_mpa: Young's modulus of the ball steel (the pair's effective modulus is half of it)
   double ballNu = 0.3;              // ball_nu
   double ballAlpha = 0.05;          // ball_alpha: Hunt-Crossley damping, s/m (rebound ratio about 1 - alpha * v)
   double ballTanDelta = 0.0;        // ball_tan_delta: loss tangent of the steel pair (0 = elastic; a damping ratio tan_delta / 2 at the contact's own frequency, so a slow bounce loses energy too)
   double ballRollingResistance = 0.003; // ball_rolling_resistance: one ball rolling on another, torque = coefficient * reff * N (steel on steel, same class as rolling_resistance)
   double ballMu = 0.15;             // ball_mu: steel on steel friction
   double ballSubstepUs = 2.0;       // ball_substep_us: the smallest substep a pair may need (its step is omega_H / omega of the contact, see SetupIsland)
   double ballRestSubstepUs = 6.0;   // ball_rest_substep_us: no longer used, still read
   double ballReachMm = 1.0;         // ball_reach_mm: balls closer than this (surface gap) are stepped as one cluster
   bool sleep = false;               // sleep: 1 = a group of touching balls at rest sleeps (not stepped) until something can change it (runs/2026-11-03-pile-sleep/design.md)
   double sleepSpeedMmS = 0.1;       // sleep_speed_mm_s: every ball of the group under this speed and surface speed of spin (requirement 4.1, 'still')
   double sleepHoldMs = 100.0;       // sleep_hold_ms: for this long without a break (a ball only passing through zero speed does not sleep)
   // Arrival-time substeps for islands of static surfaces and balls: a step never carries a ball into a surface (it ends at a small gap when an approach is predicted within it), a contact
   // is stepped by its own stiffness at the depth the impact will reach, and a ball far from every surface takes the rest of the cycle as one step.
   bool arrivalStep = true;          // arrival_step (VPX_ARRIVAL=0/1 overrides it)
   bool staticModel = true;          // static_model (default 1 since 2026-11-03): 1 = every static surface (walls, primitives, ramps, playfield) takes a ball to the group step as a rubber does (ring_model, GroupStepDesign 14); 0 = only rubbers
   int rollModel = 1;                // roll_model (default 1): 0 = the rolling and twist torques smoothed through zero spin (roll_smooth_rad_s); 1 = rolling and twist stick springs at the contact (group step)
   int pivotModel = 1;               // pivot_model (default 1): 1 = a 1-dof part's pin (gate, spinner, wire switch) holds on a stick spring and slides at mu_dynamic_ratio times its static limit; 0 = the slip direction, stuck state and turn event
   double epsRestUm = 0.01;          // eps_rest_um (default 0.01): the position error a steady contact's step may carry (step_model 1, step_swing_cap), um
   bool stepSwingCap = false;        // step_swing_cap (default 0; with step_model 1): 1 = a steady contact whose whole swing (2 x its mode's amplitude) is within the step's error allowance sets no step limit
   bool hingePieceSum = false;       // hinge_piece_sum (default 0): 1 = a hinge part's contact (wire switch, gate, spinner) is the union of its pieces, each stepped by its own exact difference (runs/2026-11-05-hinge-piece-sum/design.md)
   int stepModel = 0;                // step_model (default 0): 0 = a steady contact is stepped at omegaHRest / its frequency at the expected depth; 1 = by the motion in its mode: h <= sqrt(12 eps / A) / omega at its real depth (A the mode's amplitude), never shorter than rule 0
   double rollSmoothRadS = 0.05;     // roll_smooth_rad_s: the rolling torque's smoothing through zero spin in the group step, rad/s
   double arrivalDepthUm = 2.0;      // arrival_depth_um: the deepest a ball may start a contact before the step is taken again (rewind), um
   static const ContactSettings& Get();   // checks the file for changes about twice a second of game time
   static void StartPlay(const std::filesystem::path& tableFile);   // the start of a play, before any part is built: the file as it is now
   // The friction impulse a contact applies for the impulse `x` that would cancel the slip (signed): inside the static limit it sticks (x); beyond it, it slides at the dynamic limit
   // (muDynamicRatio times the static one). With the settings off (stock engine) it is the plain clamp to the static limit.
   float LimitSlip(const float x, const float maxImpulse) const
   {
      if (!enable || fabsf(x) <= maxImpulse)
         return enable ? x : fmaxf(-maxImpulse, fminf(x, maxImpulse));
      return copysignf((float)muDynamicRatio * maxImpulse, x);
   }
};
