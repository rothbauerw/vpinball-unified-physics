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
   // Arrival-time substeps for islands of static surfaces and balls: a step never carries a ball into a surface (it ends at a small gap when an approach is predicted within it), a contact
   // is stepped by its own stiffness at the depth the impact will reach, and a ball far from every surface takes the rest of the cycle as one step.
   bool arrivalStep = true;          // arrival_step (VPX_ARRIVAL=0/1 overrides it)
   bool staticModel = false;         // static_model: 1 = every static surface (walls, primitives, ramps, playfield) takes a ball to the group step as a rubber does (ring_model, GroupStepDesign 14); 0 = only rubbers
   double rollSmoothRadS = 2.0;      // roll_smooth_rad_s: the rolling torque's smoothing through zero spin in the group step, rad/s
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
