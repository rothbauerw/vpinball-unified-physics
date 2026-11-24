# Unified physics to-do

Open items only, one line each (300 characters at most), numbers kept. Finished = deleted; git history is the archive.
No progress notes, results or time stamps here: those go in the lead log, runs/<date>-lead/log.md.
History before 2026-11-07: docs/unified-physics/archive/ (three old to-do files).
A hook refuses an edit that breaks these limits (.claude/hooks/check_todo.py).

## Goals (shared engine; what must be finished, in order; changes only with the user)
G1. One shared path: the old-law switches deleted (10.13).
G2. Energy by construction: energy enters only through drives; EnergyStep removed once no state needs it.
G3. No stock in converted objects (stock_routes.md).
G4. Contact system runs once per 1 ms cycle, not per stock sub-interval.

## Now
Touch prediction (1.52): design 3.8 + ball touches on the go-back pushed (6988fb3e1; a9ee0bd50 comment pushed with bca6c9925); Gorgar load clean. Worst Gorgar cycle = BottomJet ring kick on its contact clock, on and off (1.40, parked by the user).
Open (lead 14): BodyAxis held energy vs step (stops lead).
Step halvings (1.42, 1.50): averaged push (push_average, bca6c9925) cleared stress v2 76690 and 249784 (tested on 2df841ee4, not on 6988fb3e1). Snapshot test and cost not run. Rest: 1.50 below. runs/2026-11-23-lead19/log.md.
Drop targets: branch drop-axis (wt-drop-axis, not merged). PAUSED (user 11-22): hits 15-30x over 30 us. To be rewritten on the shared path to ModelStandards; drop-axis code not reused (user 11-22). DropTargetDesign 6 STATUS.
N7. Stress v2 bench: ball travel energy +2 % at 57 ms; spin not counted. Check total energy in the saved run (no new run).
N8. Test driver must refuse a ball placed inside any object, unless the file declares the overlap (rule 10b: then judged on recovery, no energy gained).
N12. Physical kicker hole sizes (user 2026-11-17, 'extend at some point'): today fixed at the Sys7 18.93 mm; another size needs a cup that follows it. Paper first.
N10. Performance tests: not built; review TestMatrixDesign with the user when building them. Includes the stress table driver clock after start-up and the 144 Hz frame budget. Any safety-net trip (e.g. the cord's hold and let-go) fails the run.

## 1. Contact system (tasks)
38. Overlap replay (CheckOverlaps, PhysicsEngine.cpp:1099): count how often it trips on stress v2 and Gorgar; never: delete it; else judge a replayed hit against the same hit without.
39. Stress v2 cost over rule 5: worst 7-cycle window 6.95-7.69 ms focused (limit 3 ms), worst cycle 2.4-2.6 ms (runs/2026-11-22-ringcap/cost). Cause not traced.
10. Sleep (default off): a ball on a flipper never sleeps (bats have no rest check). Whether sleep becomes default: user.
13. Moving script/animated primitives: details from the user.
27. Shared stop piece (runs/2026-11-22-stops/design.md): (a) done d408b70ac; (b) on branch stops-wip (ad9578d64), not merged: fix the released-row search first (s4_results 5-6). Then the skirt (1.40), then (c) = 1.52 (user 11-22).
32. Env switches that change behaviour (rule: diagnostics only): VPX_PREVENT, VPX_VISCOELASTIC, VPX_OVERLAP_REPLAY. Remove each or make it an ini key.
31. After G1: check every converted part (flipper, gates, spinners, bumper, disc, rubbers, cord) against the ModelStandards definition of done (top of the file).
30. part_guess_accel (default off): design evaluation for all moving parts, then a test plan if kept. After the shared-engine items, unless one of them needs it.
40. Jet bumpers PARKED (user 11-23). Skirt seat round 4 on branch stops-wip (0ee3cd69d). Struck cycles: ring kick ~45 % of steps, seat ~20 % (runs/2026-11-23-struck-skirt). Halving at the lip seat's cone point, skirt rocking at rest (stress gb4, Bumper2, t 20120).
42. Step halvings: stress v2 0 to 260 s and Gorgar 0 with push_average (bca6c9925). Recheck the 11-15 repro: slow ball rocking on a slope below the flippers, 1,269 halvings (runs/2026-11-15-finder-P6b3/drive_burst.txt).
43. After the finder work (1.34): tables/turn table/turn table test 4 (user's worst-case table) runs slower on today's build than on the 11-09 build (user's play). Measure both on one machine and find the cause.
48. Joint zero floor (built f1397fb66): paper check that a joint of 3+ wires never adds energy at the raise between cycles; shared-line negative terms get no floor.
41. Magnet pull in the energy reach is taken at the ball's place, not the largest over its travel: a true bound needs current, saturation, in-winding and per-type terms. Re-looks while a magnet is on (stress C1).
49. Jacobian: rolling resistance's torque has no slope with depth (Kr ~ a^3, its damper, its turn at the limit); every contact. Measured 1e-4 of the row (runs/2026-11-19-standup-p2/fix2).
50. Mid-probe slopes left after push_average (static surfaces fixed, bca6c9925): moving parts (0.36 of the row at a stand-up face crossing, fix2), ball pairs, terms path; friction plane still the mid normal, flips at a corner tie (runs/2026-11-23-lead19/endslope).
52. Stop search -> event prediction: one predictor for part stops and every ball touch on the shared path (runs/2026-11-22-touch-predict). Design of 11-21 reused, drop-axis code not.
54. Playfield and turntable damp the bounce by their own alpha law (contacthost.cpp LoadPlayfieldSettings, TurntableSettings); every other hard surface uses the viscoelastic law: use the shared one (user 2026-11-21). Needs a playfield restitution (see 9.7).
53. One 'on the stop' rule for the engine (design 1.52 Q7): 35 checks and 5+ nearness constants decide it today (runs/2026-11-21-onstop-read/onstop.md). Paper first.
56. Jet ring stops (user 11-23): pinned dead or bounce? Today elasticity 0.4 rest / 0.3 coil stop, code defaults, not measured (contactbumper.h:58). Kept on the bounce rule (piece 3, design 3.2d rule 1) until investigated.
55. Cost of a ball: 21 us/cycle on the bench (stock 2.5), ~90 in play; worst-cycle steps: bounces 49 %, seats 21 % (1.40), approach 17 % (1.52), steady 11 % (step-inventory). Bounce: one impulse at touch. Evaluations per object per cycle up to 598 (Hgo2_19, time3).

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
6. Deep late contacts on gates (Gate002 317 um at 1.08 m/s); check the reach test's missing part acceleration (runs/2026-11-10-reach-accel/design.md section 4).
7. Plate gates: bottom edge at or above the ball's centre at rest, the wire-gate rule (user 2026-11-19). Today 20 mm from a 31.75 mm pivot puts it 1.1 mm below at the 15 deg rest: a back hit lifts the ball over Gorgar Ramp2 (runs/2026-11-19-laptop-lead/gate).
8. Ball-gate restitution comes from the stock Elasticity (Gorgar 1.1 gives damping 0, no loss in the contact); needs a physical steel-on-steel value (user 2026-11-19).

## 6. Magnets (rules: ModelStandards section 8)
2. Energy by construction for magnets (goal G2): pair force not energy exact; iron-core and plug cases unchecked.
7. Long rod core, moving magnets.

## 7. Turntable
1. Constant-speed AC motor option.

## 8. Testing
1. Accuracy bench: add speed error (arrival timing). 24 bounce step faults, all under 1 mm, parked by the user (judge_RA3).
5. MCP server (later).
9. 750 bench: 4 energy lines (Hgo1_2, Hgo1_29 87.1 mJ; Hgo2_55 24.9; Hgo2_66 53.4) = ball 2.3-3.8 mm in Rubber4 at first sight; likely replay starting mid-hit (inferred). Confirm, then start earlier or declare (rule 10b).
8. Build-compare script compiles /fp:precise by default, so a reported difference is logic, not rounding (ModelStandards trap: a byte difference is not a logic difference).
10. A state carries between replays in one chunk: Hgo5_77 differs only after earlier balls change, identical when run alone; find what carries over (runs/2026-11-10-old-loop/results_piece1.txt).
12. Snapshot self-test fails on Gorgar: 30 of 8001 steps differ (19 with free_flight), "contacts (3 engaged balls)", before free flight too. Find what the snapshot misses. Lead: a retake dropped JewelRamp from a ball's list (aircore-ff/snap_trace).
11. Accuracy bench for moving parts: real hits on flipper bats, jet bumpers, gates, spinners, wire switches and sling plates, judged by impulse against the 1 us reference as the static bench is (today's 750 hits are static surfaces only).
13. 750 replays: ~390 halvings (408 on 11-19, 389 at 6244dbce1), mostly wire switches sw33 / sw32 with a ball; 0 in full Gorgar play on the same hits: find what the replay start misses.

## 9. Still stock (ask before converting)
1. Drop targets (stand-up targets converted: StandupTargetDesign).
6. Playfield material values unverified: compare with the user's video of a ball rolling on a playfield.
7. Playfield bounce too high (user 2026-11-19): a falling ball rebounds at 0.77-0.98 of its fall speed at 0.28-0.56 m/s (Gorgar, runs/2026-11-19-laptop-lead/gate run_F1, t 3559-3871). No measured value in PartData; find where the restitution comes from.

## 10. Branch and housekeeping (release path 3, 5, 10: see ReleasePlan.txt)
3. unified = 0 against a plain upstream build; plunger fix.
5. Stock-file hooks to separate.
6. Editor changes to replace the ini files (panels, settings saved in the .vpx); unified off per table for release.
7. Old keys in 14 inis.
10. Review ReleasePlan with the user, then PR to the official repo after 3-5.
13. Old-law switches to delete once accepted (CLAUDE.md): roll/pivot/step_model 0 (stock's copies of converted parts are out of the path: balls fall through; load warns).

## 11. Switches and triggers
1. Light-part contact cost (limit 30 us): wire switch 54-181 us, Button 271 us (runs/2026-11-10-rollover/S1); drop target 3 m/s ~180 us, its pinch ring is physical (runs/2026-11-19-drop-rattle). Midpoint solver follows the ball-on-part clock: needs another solver class. Parked to N10.

## 12. Kickers (when kickers are extended)
1. Physical cup settles too slowly (user's play): find the missing damping, maybe the switch blade. Resting cost 102 us per cycle is the contact system's (runs/2026-11-12-lead).

## 13. Ramps (converted: floor and side walls as a hard surface, ramphardsurface.cpp)
1. Performance: cost with a ball rolling and resting on a ramp never measured; measure against rule 5 (worst cycle, worst 7-cycle window).
2. Wire ramps as wires (enhancement, not a priority; user 2026-11-17): today a flat floor with walls, as stock; the ball should ride on the rods, a tube contact along each wire.
3. Holes in ramps: refused with a notice (user 2026-11-17). Physics cuts in ramps left unexplained stray edges (mesh holes P2); a curved ramp's drawing needs per-face clipping.

## 14. Plungers (PlungerDesign.txt; steps 1-3 built)
1. Launch cost FAILS H9: rod worst cycle 305 us vs 30 us, 17 steps/cycle on the tip's stiffness clock (e1-why); auto fork 147 us, 24 steps. Same root as 11.1. User 2026-11-18: does not hold up the model.
2. LimitBOS / LimitEOS: which end each means, and re-firing on the rod's bounces (10-16 per launch); no table in tables/ uses them. Left as built (BOS rest, EOS full pull); user 2026-11-18: later.
3. Release rule (2k MERGED) passes the user's regression set on one board (PinOne). Other boards, a weaker cabinet spring and other players' moves untested: add their recordings to runs/2026-11-18-regression when available.
4. The same 2 mm key launch gives 0.513-0.580 m/s between bench runs (e2d, e2f, e2g); INFERRED driver start offset, unproven. Understand before any release design.
5. Extension (user 2026-11-18, later): design 2h, a cabinet plunger on an auto-plunger table acts as the Launch Ball button (pull-and-release or a press past rest).
6. Auto plunger not rechecked with a clean ball: the tip head, a lower energy (0.64 J), the fork drawing on screen (a3c-clean ran the fork at 1 J only).
7. Plunger button lost control once (user 2026-11-19): Enter went down/up every 308 ms (tables/logs 183700). VPX drops OS repeats and acts on changes only: source outside VPX (runs/2026-11-19-laptop-lead/plunger_input.md). Rod log lines capped at 40 per load (plunger.cpp:163).
8. Other tables, auto plunger, mechanical input (moved from 9.4: the plunger is converted).