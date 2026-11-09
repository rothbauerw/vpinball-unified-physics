# Unified physics to-do

Open items only, one line each (300 characters at most), numbers kept. Finished = deleted; git history is the archive.
No progress notes, results or time stamps here: those go in the lead log, runs/<date>-lead/log.md.
History before 2026-11-07: docs/unified-physics/archive/ (three old to-do files).
A hook refuses an edit that breaks these limits (.claude/hooks/check_todo.py).

## Goals (shared engine; what must be finished, in order; changes only with the user)
G1. One shared path: the old loop and the old-law switches deleted; one shared step rule (today 4 drifted copies).
G2. Energy by construction: energy enters only through drives; EnergyStep removed once no state needs it.
G3. No stock in converted objects: the 13 stock routes closed (stock_routes.md); stock detection off for them.
G4. Contact system runs once per 1 ms cycle, not per stock sub-interval.

## Now
In flight: nothing running. grip_in_step committed (1cec677ab). energy_decay deleted (user, 11-09; snapshot backups/2026-11-09-before-delete-energy-decay). Uncommitted: clamp-trace twist fields (groupstep.h, diagnostic). CLAUDE.md, hooks and memory are outside git.
Next decision: trace Gorgar's ball-only halvings (a ball touching down on the playfield mid-step, 3897 ms; balls against Wall6 at the Drain kicker).
N1. Wire switch: switch bench now 0 failed attempts (11-09 commits c54f837ef, eef728db8, 08e00900d, 07b40b3c9). Cost: Gorgar's not a concern (user, 11-09); the switches' own cost vs the 30 us/ms limit still open, not yet measured per ms.
N12. Gorgar load: cord-grip halvings fixed (1cec677ab); 3-4 ball-only halvings left (playfield touch-down mid-step; Wall6 at the Drain kicker); ~60 depth warnings < 20 um in the bottom pile and at the Drain kicker (stock kicker).
N2. Wires B, C, D as the real switches (B V wire, C angled straight wire, D arch through playfield); lever ends wrong. Derive sizes from photos (ball 1-1/16 in as scale) and parts drawings; user checks before build.
N3. Wire switches: option for a playfield hole under each wire.
N4. Inder blade: hinge sizes and skin are estimates; 10-17 us per step; crawls from rest. Derive sizes; user checks before build.
N5. One reach rule from physics (speed plus force) for grouping and step length. Grouping: Gate001's 81 mm swing as one-cycle travel joins two lanes. Step length: speed only, no force. Check all parts; ball pairs of unequal mass and resting-ball wake (read only).
N6. Bounce too lively on long steps (bench: rebound 0.83 vs law 0.51). 11-09: the deleted energy_decay fixed 7 wall hits on the 750 bench (one 7.2 -> 1.0 deg): long-step wall ringing is real; needs a physical fix. Includes steel impacts (VPX_OMEGAH_SCALE).
N7. Stress v2 bench: ball travel energy +2 % at 57 ms; spin not counted. Check total energy in the saved run (no new run).
N8. Test driver must refuse a ball placed inside any object, unless the file declares the overlap (rule 10b: then judged on recovery, no energy gained).
N9. Branch pushed: every commit is on the fork. Engine work (part_guess_accel, clamp trace) uncommitted.
N10. Performance tests: not built; review TestMatrixDesign with the user when building them. Includes the stress table driver clock after start-up and the 144 Hz frame budget. Any safety-net trip (e.g. the cord's hold and let-go) fails the run.

## 1. Contact system (tasks)
10. Sleep (default off): a ball on a flipper never sleeps (bats have no rest check). Whether sleep becomes default: user.
13. Moving script/animated primitives: details from the user.
17. Shared-loop cost per cycle (waits until models work).
26. Ball started inside: the solve must never keep a runaway answer. May be handled already: check.
27. Shared stop piece (gate, spinner, disc, ring, bat, plunger).
29. Counter columns filled only on the old path.
32. Env switches that change behaviour (rule: diagnostics only): VPX_HINGES (=0 gives the stock kick), VPX_PREVENT, VPX_PROXIMITY, VPX_VISCOELASTIC, VPX_OVERLAP_REPLAY, VPX_ARRIVAL_NOSPLIT; VPX_BUMPERS only warns. Remove each or make it an ini key.
31. After G1: check every converted part (flipper, gates, spinners, bumper, disc, rubbers, cord) against the ModelStandards definition of done (top of the file).
30. part_guess_accel (default off): design evaluation for all moving parts, then a test plan if kept. After the shared-engine items, unless one of them needs it.

## 2. Flippers
1. Rubber separation at the end stop: details from the user.
2. End stops as impact law.
3. Coil winding coupling.
4. Ball control after a flipper hit (later).
5. Binding that limits top speed: details from the user.
6. Shaft torque in the bushing: details from the user.
7. Solenoids for every era: coil data from the user.
9. Rest angle -79 vs -70 deg, crank 25.4 vs 22.2 mm: user to say.
11. Two balls on a raised bat; table close with balls on a bat (check).
12. Real Gorgar flip data: which clock is time zero for the 0.675 ms first movement (fits put the first 0.1 deg at 6-7 ms); coil delay and rise time are fitted against it.
13. Flipper inertia computed per flipper from its parts (bat from length, base and tip radius; crank, plunger, link), each mass a setting; methods in PartData section 1. Replaces the fixed i_bat and i_crank; then recheck the empty stroke.

## 3. Bands and slings
2. One Switch event instead of SwitchClosed/Opened.
3. SlingTrigger to the jet bumper Solenoid form (check).

## 5. Gates and spinners
1. Cases not covered: gate_hold, gap_back, overlap, spinner switches on a table, two-way and plate gates.
2. Cost over 30 us; stops as impact law.
5. Redraw to match the real parts.
6. Deep late contacts on gates (Gate002 317 um at 1.08 m/s).

## 6. Magnets (rules: ModelStandards section 8)
2. Energy by construction for magnets (goal G2): pair force not energy exact; iron-core and plug cases unchecked.
7. Long rod core, moving magnets.
9. Guide tube ConstrainBall sets held balls' positions outside the group step (PhysicsEngine.cpp:866): breaks the shared path. Code read only.

## 7. Turntable
1. Constant-speed AC motor option.

## 8. Testing
1. Accuracy bench exists (750 hits vs 1 us); 36-39 hits always off the 5 % limit: fix them or set the limit (user).
5. MCP server (later).
9. 750 bench: 4 energy lines (Hgo1_2, Hgo1_29 87.1 mJ; Hgo2_55 24.9; Hgo2_66 53.4) = ball 2.3-3.8 mm in Rubber4 at first sight; likely replay starting mid-hit (inferred). Confirm, then start earlier or declare (rule 10b).
8. Build-compare script compiles /fp:precise by default, so a reported difference is logic, not rounding (ModelStandards trap: a byte difference is not a logic difference).

## 9. Still stock (ask before converting)
1. Stand-up and drop targets.
2. Kickers.
3. Switches and triggers: wire switch work in progress (N1-N4).
4. Plunger: other tables, auto plunger, mechanical input.
5. Ramps: wire-ramp types, rebound, cost.
6. Playfield material values unverified: compare with the user's video of a ball rolling on a playfield.

## 10. Branch and housekeeping (release path 3, 5, 10: see ReleasePlan.txt)
3. unified = 0 against a plain upstream build; plunger fix.
5. Stock-file hooks to separate.
6. Editor changes to replace the ini files (panels, settings saved in the .vpx); unified off per table for release.
7. Old keys in 14 inis, then the unknown-key notice.
10. Review ReleasePlan with the user, then PR to the official repo after 3-5.
13. Old-law switches to delete once accepted (CLAUDE.md): static/roll/pivot/step_model 0, island_reach 0, ball_pair_reach 0, ball_rest_substep_us.
