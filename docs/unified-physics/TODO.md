# Unified physics to-do

Open items only, one line each (300 characters at most), numbers kept. Finished = deleted; git history is the archive.
No progress notes, results or time stamps here: those go in the lead log, runs/<date>-lead/log.md.
History before 2026-11-07: docs/unified-physics/archive/ (three old to-do files).
A hook refuses an edit that breaks these limits (.claude/hooks/check_todo.py).

## Goals (shared engine; what must be finished, in order; changes only with the user)
G1. One shared path: the old-law switches deleted (10.13).
G2. Energy by construction: energy enters only through drives; EnergyStep removed once no state needs it.
G3. No stock in converted objects: the 13 stock routes closed (stock_routes.md); stock detection off for them.
G4. Contact system runs once per 1 ms cycle, not per stock sub-interval.

## Now
In flight: 1.34 first (shared finder; P1-P3, P9 static meshes, P4b cuts, P5 static tubes). 1.33: piece 2 kept with cut B (stress v2 cost 52 to 1.7 us per cycle), piece 3 dropped (cost). Owed: snapshot hashes refresh. Lead log runs/2026-11-10-lead.
Next decision: 1.34 P6 (rigid moving parts); P4c (turntable and kicker holes in the near pass) after the kicker lead's hole change. Refs since P5 (/fp:fast moved last digits): runs/2026-11-13-finder-P5. Separate: N3 ramp edge search (343 of 417 us per cycle, slot in groove).
Kickers: physical Sys7 cup and virtual store in the branch ([Features] kickers; KickerDesign.txt; hole list in carvedholes.cpp for P4c). Next: 10.14. Physical resting cost 102 us per cycle is the contact system's. Lead log runs/2026-11-12-lead.
Parked with playfield holes (Button and Star committed): a ball on a flush Button lands on the wood through the cap (no hole), bounce e 0.86; user expects much more loss on a Button.
N12. Gorgar load: 0 halvings (cord grip 1cec677ab, rolling patch 87ec36642, stick-zone trial 49155d9e2: Drain-kicker 2-cycle gone, rounds per step -10 %); ~60 depth warnings < 20 um in the bottom pile and at the Drain kicker (stock kicker).
N3. Playfield cuts (slot, groove, rectangle, meets: bbdee9278). Open: slot-in-groove 256 us/cycle; normal jump at slot tips (cut-junction/results_ramp.txt); ball held where a slot crosses a rect's edge gains energy (finder-P4b S20, S06). Later: draw new shapes, triangle, holes.
N5. One reach rule from physics (speed plus force) for grouping and step length; first built in 1.34 P2. Grouping: Gate001's 81 mm swing as one-cycle travel joins two lanes. Step length: speed only, no force. Check all parts; ball pairs of unequal mass and resting-ball wake (read only).
N7. Stress v2 bench: ball travel energy +2 % at 57 ms; spin not counted. Check total energy in the saved run (no new run).
N8. Test driver must refuse a ball placed inside any object, unless the file declares the overlap (rule 10b: then judged on recovery, no energy gained).
N9. Branch pushed: every commit is on the fork. Engine work (part_guess_accel, clamp trace) uncommitted.
N10. Performance tests: not built; review TestMatrixDesign with the user when building them. Includes the stress table driver clock after start-up and the 144 Hz frame budget. Any safety-net trip (e.g. the cord's hold and let-go) fails the run.

## 1. Contact system (tasks)
34. One shared finder for every contact (static, moving, tubes, ball pairs): complete, short, local, one index; replaces the per-part finders. Requirements accepted: runs/2026-11-11-finder/requirements.md. Design accepted: runs/2026-11-11-finder/design.md (build P1-P9).
36. Release check's fixed 2 mm hold for a closing ball (contactsystem.cpp:707, from the port): no physics behind it; release at the ball's reach instead (mesh hosts already do with shared_finder).
37. A no-op stock hit that splits a cycle moved S00's bounce spin 3.5 % (stress v2, 2894 ms; finder-P2b S00 trace): step-split sensitivity; check against the 1 us reference. Parked.
38. Overlap replay with free_flight (every free ball stepped by the contact system, bb3898626): the replay now restores and re-steps every ball; seen only not to crash (19 replays). Judge one replayed hit against the same hit without a replay.
39. Stress v2 cost over rule 5: worst 7-cycle window 9.6-9.8 ms (limit 3 ms), worst cycle 4.8-5.3 ms, with free_flight 0 and 1 alike (runs/2026-11-12-freeflight timing csv). Cause not traced.
35. Stress v2 diverter: both diverter walls collide (the script's switch is not run), S09 jams between them and climbs ~3x what its velocity gives (cy 49168; finder/gap/cy49168.py). Check energy (rule 10b).
10. Sleep (default off): a ball on a flipper never sleeps (bats have no rest check). Whether sleep becomes default: user.
13. Moving script/animated primitives: details from the user.
17. Shared-loop cost per cycle (waits until models work). Measured: evaluations per object per cycle median 4, p90 16, worst 598 over 144 steps (750 bench, Hgo2_19; time3/results.txt); cause unknown.
26. Ball started inside: the solve must never keep a runaway answer. May be handled already: check.
27. Shared stop piece (gate, spinner, disc, ring, bat, plunger).
32. Env switches that change behaviour (rule: diagnostics only): VPX_HINGES (=0 gives the stock kick), VPX_PREVENT, VPX_PROXIMITY, VPX_VISCOELASTIC, VPX_OVERLAP_REPLAY; VPX_BUMPERS only warns. Remove each or make it an ini key.
33. Several contacts on one object (mesh faces and edges, tube folds, cut rims): each true local minimum of the ball-surface distance its own contact and friction, merged only where patches overlap; replaces the fold-sum fade. Also: sling post cos (cordfixed.h:949). Check: 750-bench rubber hits.
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
6. Deep late contacts on gates (Gate002 317 um at 1.08 m/s); check the reach test's missing part acceleration (runs/2026-11-10-reach-accel/design.md section 4).

## 6. Magnets (rules: ModelStandards section 8)
2. Energy by construction for magnets (goal G2): pair force not energy exact; iron-core and plug cases unchecked.
7. Long rod core, moving magnets.
9. Guide tube ConstrainBall sets held balls' positions outside the group step (PhysicsEngine.cpp:866): breaks the shared path. Code read only.

## 7. Turntable
1. Constant-speed AC motor option.
2. Cost: with a turntable hole the playfield's probe tests every face on every call (hitrubber.cpp:981, the short list is dropped). Keep the short list and leave out points inside the hole, as the cuts do.

## 8. Testing
1. Accuracy bench, judged by impulse: all 750 pass (worst 0.68 mm at the next contact, limit 1 mm). Step or fault by path error over VPX_STEP_SCALE 1/0.5/0.25: 24 bounce faults, all under 1 mm, parked by the user (judge_RA3); add speed error (arrival timing).
5. MCP server (later).
9. 750 bench: 4 energy lines (Hgo1_2, Hgo1_29 87.1 mJ; Hgo2_55 24.9; Hgo2_66 53.4) = ball 2.3-3.8 mm in Rubber4 at first sight; likely replay starting mid-hit (inferred). Confirm, then start earlier or declare (rule 10b).
8. Build-compare script compiles /fp:precise by default, so a reported difference is logic, not rounding (ModelStandards trap: a byte difference is not a logic difference).
10. A state carries between replays in one chunk: Hgo5_77 differs only after earlier balls change, identical when run alone; find what carries over (runs/2026-11-10-old-loop/results_piece1.txt).
12. Snapshot self-test fails on Gorgar: 30 of 8001 steps differ (19 with free_flight), "contacts (3 engaged balls)", before free flight too. Find what the snapshot misses.
11. Accuracy bench for moving parts: real hits on flipper bats, jet bumpers, gates, spinners, wire switches and sling plates, judged by impulse against the 1 us reference as the static bench is (today's 750 hits are static surfaces only).

## 9. Still stock (ask before converting)
1. Stand-up and drop targets.
4. Plunger: other tables, auto plunger, mechanical input.
5. Ramps: wire-ramp types, rebound, cost.
6. Playfield material values unverified: compare with the user's video of a ball rolling on a playfield.

## 10. Branch and housekeeping (release path 3, 5, 10: see ReleasePlan.txt)
3. unified = 0 against a plain upstream build; plunger fix.
5. Stock-file hooks to separate.
6. Editor changes to replace the ini files (panels, settings saved in the .vpx); unified off per table for release.
7. Old keys in 14 inis, then the unknown-key notice.
10. Review ReleasePlan with the user, then PR to the official repo after 3-5.
13. Old-law switches to delete once accepted (CLAUDE.md): roll/pivot/step_model 0, island_reach 0, ball_pair_reach 0, free_flight 0, proximity_engage 0 (with shared_finder on, stock's copies are out of the path: balls fall through; load warns).
14. Rename HitRubber: it hosts every static mesh, the playfield included. Kicker lead, right after the kicker work is committed; agree timing with the finder lead (1.34 edits hitrubber.cpp).

## 11. Switches and triggers
1. Light-part contact cost (ball on a part of a few grams): worst cycle wire switch 54-181 us, Button 271 us (runs/2026-11-10-rollover/S1); limit 30 us. Midpoint solver must follow the ball-on-part contact clock: needs another solver class (rigid contact + impact law, or L-stable). Parked to N10.