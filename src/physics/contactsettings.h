// license:GPLv3+
#pragma once

#include <cmath>
#include <filesystem>
#include <string>

namespace GroupStep { struct Settings; }

// Table-wide contact settings of the unified engine (see docs/unified-physics/archive/NewPhysicsReference.txt, [Contact]); implemented in contactsettings.cpp.

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
   double ballReachMm = 1.0;         // ball_reach_mm: balls closer than this (surface gap) are stepped as one cluster
   bool islandReach = true;          // island_reach (default 1 since 2026-11-05): 1 = a ball and a moving part share an island only if the ball can reach the part's surface within the cycle (runs/2026-11-05-island-reach/design.md); 0 = any part in the ball's near list
   bool ballPairReach = true;        // ball_pair_reach (default 1 since 2026-11-05; needs island_reach): 1 = two balls share an island when the gap can close within the cycle at their speed bounds (the island_reach closure, energy bound, BallsMayMeet); 0 = the pairwise ball_reach_mm + 2|vi-vj|T rule only
   bool sleep = false;               // sleep: 1 = a group of touching balls at rest sleeps (not stepped) until something can change it (runs/2026-11-03-pile-sleep/design.md)
   double sleepSpeedMmS = 0.1;       // sleep_speed_mm_s: every ball of the group under this speed and surface speed of spin (requirement 4.1, 'still')
   double sleepHoldMs = 100.0;       // sleep_hold_ms: for this long without a break (a ball only passing through zero speed does not sleep)
   int rollModel = 1;                // roll_model (default 1): 0 = the rolling and twist torques smoothed through zero spin (roll_smooth_rad_s); 1 = rolling and twist stick springs at the contact (group step)
   double epsRestUm = 0.1;           // eps_rest_um (default 0.1 since 2026-11-05): the position error a steady contact's step may carry (step_model 1, step_swing_cap), um
   bool stepSwingCap = true;         // step_swing_cap (default 1 since 2026-11-05; with step_model 1): 1 = a steady contact whose whole swing (2 x its mode's amplitude) is within the step's error allowance sets no step limit
   bool termArrival = true;          // term_arrival (default 1 since 2026-11-05): 1 = each term of a contact made of terms (hinge_piece_sum) not yet pressing gets the arrival rule (rule 3) by its own depth and closing rate
   bool dampDirSlope = true;         // damp_dir_slope (default 1 since 2026-11-05): 1 = Newton's slope of a terms contact (hinge_piece_sum) takes the turning of the damping's direction G = Fel / |Fel| and dd's change through it (runs/2026-11-05-damp-dir-slope/design.md); the force is unchanged
   bool leverRollSlope = true;       // lever_roll_slope (default 1 since 2026-11-08, on while tested): 1 = Newton's slope takes the two laws' terms it held fixed: the contact point's move with the depth (lever arm rA = -n (R - dm / 2), in the slip and the torque) and the rolling limit's change with N (lim = muRr R N, roll_model 1); the force is unchanged (runs/2026-11-08-lead/slope_terms_design.md); 0 = today's slope
   bool twistSlope = true;           // twist_slope (default 1 since 2026-11-08, on while tested): 1 = Newton's slope takes the twist torque's change with the depth (roll_model 1): while it sticks through the torsion spring kw1 = Ktw(d1), its damper cw and the patch's hold of tw0; at the limit 3 pi / 16 mu a(dm) N through a(dm) and N; the torque is unchanged (runs/2026-11-08-lead/twist_slope_design.md); 0 = today's slope (d Mw / d wn only)
   bool endGradSlope = true;         // end_grad_slope (default 1 since 2026-11-09, on while tested): 1 = Newton's slope of a terms contact (hinge_piece_sum) takes each pressing term's change of its END depth over the unknowns, d a1 / du, from the end place's own gradient (its normal and column slope at the step's end), in the terms' slope, G's turning and the normal force's slope kv (friction, lever-roll and twist limits); the row side (the force's direction) stays the mid gradient; the force is unchanged (runs/2026-11-08-lead/end_gradient_design.md); 0 = today's slope (the mid gradient on both sides)
   bool dampingIntegral = true;      // damping_integral (default 1 since 2026-11-09, on while tested): 1 = every depth-dependent damping of the group step in its integrated form, the law's own work over the step at a steady rate, c = (W(d1) - W(d0)) / (d1 - d0), with its exact slope c(d1) in Newton's (Kuwabara-Kono, the loss tangent's dashpot, Hunt-Crossley on the contact's own share of a bend), and the twist limit from the patch of the step's mean elastic force a(fe) instead of the mid depth's (runs/2026-11-08-lead/damping_integral_design.md); 0 = the committed forms (point values at the mid depth, the Kuwabara-Kono chord/midpoint switch, the mid depth's patch for the twist limit)
   bool rollingPatch = true;         // rolling_patch (default 1 since 2026-11-09, on while tested): 1 = rolling resistance takes the patch of the step's mean elastic force a(fe) (as the twist limit under damping_integral) for its pressure N / (pi a^2) and its gate (no switch on the mid depth), with the rolling limit's full slope in Newton's (through N and a(fe); 0 at the floor or the cap of muRr) (runs/2026-11-08-lead/rolling_patch_design.md); 0 = the mid depth's patch a(dm), the limit's slope with muRr held
   bool stickZoneTrial = true;       // stick_zone_trial (default 1 since 2026-11-09, on while tested): 1 = the group step's line search first tries, for each contact whose shear patch slides at the current point with N > 0, the fraction of Newton's move at the middle of the part of the move where its trial force lies inside mu N (a quadratic along the move, the patch's spring, damper, stretch and geometry held at the point; none from a contact whose part contains the full move's end, s = 1), ascending, then the full move, 1/2, 1/4 ...; acceptance by the exact residual as every trial (runs/2026-11-08-lead/slip_reversal_design.md section 2b); 0 = the full move, 1/2, 1/4 ... only
   bool physicsDefaults = true;      // physics_defaults (default 1 since 2026-11-10, on while tested): 1 = a wall, rubber, primitive, ramp or hit target with override physics off and no physics material on the table takes its type's VPX defaults (elasticity, falloff, friction: Settings_properties.inl), with a notice, not the blank material's 0s (a dead, frictionless object; runs/2026-11-10-physics-defaults/design.md); scatter stays 0 (the unified contact never reads it); 0 = the blank material, as stock
   bool gonzalezMass = true;         // gonzalez_mass (default 1 since 2026-11-10, on while tested): 1 = each stored energy's exact-energy correction (Gonzalez) is split over the bodies by mass, the kinetic-energy metric, so each body's speed changes by its own move and a light part takes no false force from a heavy ball's move (runs/2026-11-10-gonzalez-metric/design.md; removed the 15441 ms halvings on a wire switch); 0 = split by inverse mass (the old rule)
   bool termForceRows = true;        // term_force_rows (default 1 since 2026-11-09, on while tested): 1 = Newton's slope of a terms contact (hinge_piece_sum) takes each pressing term's force direction gdT (its mid gradient Gonzalez-corrected so that gdT . dq = a1 - a0, the direction its force acts along) as the row side of the term's slope, in the terms' slope, G's turning (gK, and dd's change through G, gW: 0 for one term) and the normal force's slope kv (friction, lever-roll and twist limits); the column side stays d a1 / du; the force is unchanged (runs/2026-11-08-lead/term_rows_design.md); 0 = the mid gradient tg as the row (today's slope)
   bool gripInStep = true;           // grip_in_step (default 1 since 2026-11-09, on while tested): 1 = a part's grip limit (the cord's mu Tm Gfree) is worked out inside each of the group step's evaluations from that evaluation's own mid tension and mid geometry, as a contact's mu N is, with Newton's slope through it (runs/2026-11-08-lead/cord_grip_design.md); 0 = the limit refreshed between the solve's rounds from the last round (the line search is then off in any round it moved)
   bool patchPath = true;            // patch_path (default 1 since 2026-11-09, on while tested): 1 = every patch-sized coefficient of the group step's contact (the shear spring Kt ~ a and its damper ~ sqrt(Kt), the twist and rolling springs ~ a^3 and their dampers) is its mean over the depth the step travels, (P(d1) - P(d0)) / (d1 - d0), P its integral from 0 (short moves by quadrature), with the patch rule (growing holds the force, shrinking holds the stretch) from the start's patch to the mean before the step's move and from the mean to the end's patch after it, the end stretch stored on the end's patch; Newton's slope takes the means' exact slopes (runs/2026-11-08-lead/patch_state_design.md); 0 = the end depth's coefficients Kt(d1), Ktw(d1), Kr(d1) and their dampers for the whole step
   bool hingePieceSum = true;        // hinge_piece_sum (default 1 since 2026-11-05): 1 = a hinge part's contact (wire switch, gate, spinner) is the union of its pieces, each stepped by its own exact difference (runs/2026-11-05-hinge-piece-sum/design.md)
   bool partGuessAccel = false;      // part_guess_accel (default 0, 2026-11-06): 1 = a moving part's (wire switch, gate, spinner) first guess carries its last step's acceleration, (rate + a h / 2) h, as a ball's does; Newton's start only
   int stepModel = 1;                // step_model (default 1 since 2026-11-05): 0 = a steady contact is stepped at omegaHRest / its frequency at the expected depth; 1 = by the motion in its mode: h <= sqrt(12 eps / A) / omega at its real depth (A the mode's amplitude), never shorter than rule 0
   double rollSmoothRadS = 0.05;     // roll_smooth_rad_s: the rolling torque's smoothing through zero spin in the group step, rad/s
   double arrivalDepthUm = 2.0;      // arrival_depth_um: the deepest a ball may start a contact before the step is taken again (rewind), um
   static const ContactSettings& Get();   // checks the file for changes about twice a second of game time
   static void StartPlay(const std::filesystem::path& tableFile);   // the start of a play, before any part is built: the file as it is now
   void ToStep(GroupStep::Settings& s) const;   // the group step's settings from these keys and the step floor: every group the engine steps, an island's and a part's alone (GroupStep::EngineSettings)
   // The friction impulse a contact applies for the impulse `x` that would cancel the slip (signed): inside the static limit it sticks (x); beyond it, it slides at the dynamic limit
   // (muDynamicRatio times the static one). With the settings off (stock engine) it is the plain clamp to the static limit.
   float LimitSlip(const float x, const float maxImpulse) const
   {
      if (!enable || fabsf(x) <= maxImpulse)
         return enable ? x : fmaxf(-maxImpulse, fminf(x, maxImpulse));
      return copysignf((float)muDynamicRatio * maxImpulse, x);
   }
};

class Material;
class HitObject;
class PinTable;
// An object whose physics come from its physics material (override physics off): with no such material on the table (the blank one, all 0) and the unified physics on,
// its type's VPX defaults (ContactSettings::physicsDefaults); falloff < 0: the type has none. A notice once per object and play. True when they were used
bool UsePhysicsDefaults(const PinTable *table, const Material *mat, const char *kind, const std::string &name, float elasticity, float falloff, float friction, HitObject *obj);
