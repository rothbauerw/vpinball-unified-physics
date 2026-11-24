# Testing playbook

How to build, run, measure and diagnose today. When a practice changes, edit its line in place; the reason goes in the commit message.
No dates, no "learned", no superseded text. What a test is for, its cap and hand-over: DevelopmentRules 4, 10, 11.

## 1. Which tool
- **Real bench** (section 3): VPX itself with no frames. Every behaviour and accuracy question goes here first. Deterministic: same exe, table, ini and driver file give the same bytes.
- **VPX driver run** (section 4): only what the bench does not cover (list in section 3), kept to one named question.
- **The user's play** (section 6): only when the user asks.
- Build and test BGFX only (physics is identical on GL); GL only for a renderer question.

## 2. Build
- Worktree `C:\Pinball\VPX-Physics\unified-physics-branch`; `runs\` and `tables\` are in the parent folder, outside git. Build folder `build-bgfx` -> `build-bgfx\Release\VPinballX_BGFX64.exe` (`build` is the GL one).
- From the build folder: `"C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" --build . --config Release --target vpinball -- -m 2>&1 | grep -E "error C|error LNK"`
  (`-- -m`, not /m, in Git Bash; grep "error" alone matches the /errorReport flag.) LNK1104 = a VPX is running on that exe. Check the exe's time with `ls -la` after every build.
- Temporary instrumentation is marked `/* TEMP */` and grepped for before finishing; a `//` comment after code can swallow a closing brace.
- Python source edits: keep CRLF; exact-string replace only.

## 3. Real bench
- Folder `runs\2026-10-31-real-bench` (README.md, RESULTS.md). Run from PowerShell or a script file, never Bash with `C:\...` paths (it strips the backslashes):
  `powershell -ExecutionPolicy Bypass -File runs\2026-10-31-real-bench\bench.ps1 -Name x -Table <full> -Ini <full> -Drive <full> -Dir <run folder>` (default `-Exe` is `main_exe\VPXBench_BGFX64.exe`;
  bench.ps1 refuses it when it is older than the branch's BGFX build: run `python bench_snapshot.py`).
  Options: `-Trace 1`, `-Vbs file`, `-StopMs N` (engine time = drive t + 3003), `-PreRollMs N` (default 3000). Output: `run_x\out\balls.csv, report.txt, cycles.csv, done.txt`, `run_x\vpinball.log`.
- The `VPXBench_*` name keeps run.ps1 and wait loops from taking it for a VPX. Many runs at once is fine, except on a table new to `runs\_tables`: run one first (the first runs race on INDEX.csv).
- Confirm a run started: the log says `Test driver: N actions` with N above 0. Cost on vs off: compare work counts (VPX_GROUP_STATS, finder totals), not wall time.
- **Patch not in the branch**: a reusable worktree per lead, `C:\Pinball\VPX-Physics\wt-bench` or `wt-bench-2` (own `build-bgfx`, incremental). Claim it first: read
  `IN_USE.txt` in its root; if absent, write session, task, time; delete it when done. Only the claimant resets or builds there. Check out the commit, patch, build,
  exe to `Release\VPXBench_<tag>.exe` for `-Exe`. Never compile in a run folder.
- **Accuracy bench**: `runs\2026-11-01-roll-smooth`, 750 real hits each replayed alone against a 1 us reference, judged by impulse; how to run and judge: BenchTools.md.
- **Not covered by the bench** (use a VPX run): cost against the frame limit, drawing, real input, table restart and teardown, scripts that act on the game, script and frame timers,
  drop targets (they fall per frame), ROM tables. A pile with address-order splits (Gorgar pile) does not repeat on the bench, nor in VPX.

## 4. VPX driver run
- `runs\2026-10-09-slow-contact\run.ps1` with `-Drive`, every path full (a run with no driver says why with `-NoDriver "<reason>"`); command and options: BenchTools.md.
- Behaviour runs are unfocused (run.ps1 sets VPX_PLAY_UNFOCUSED=1, which also turns off the anti-hung time skip). Cost runs: `-Focused -Timing` (prefs\stress-timing: vsync off, 144 fps cap); tell the user first so the machine is idle.
- run.ps1 writes `[Features] unified = <VPX_UNIFIED, default 1>` into the run's ini and points VPX_PHYSICS_INI at it. Check every run: the ini has one `[Features]`, the log says `loaded from` (not `wrote default settings`) and `unified engine on`.
- Never copy a table: run.ps1 plays `runs\_tables\<hash>.vpx` itself; a -Vbs or bench run links it as t.vpx (a full link limit adds `<hash>-2.vpx`). Never save these in the editor.
- Another session's VPX may be open: wait for it, never close it. A leftover VPX may be closed unless the user said they are using it. Never start a run that can show a window while the user plays.
- Ini variants: copy the table's ini and change one key; a key name can repeat in other sections (`contact = 1` is in [Features] and [Rubbers]), so match the section.

## 5. Test driver input (VPX_TEST_INPUT)
- File of actions by physics millisecond, `<t_ms> <action> <args>`, `#` comment, same-time actions in file order, applied at the start of the cycle (same result at any frame rate). Code: src/physics/testdriver.cpp; design TestMatrixDesign section 2.
- Actions: `ball NAME x y z vx vy vz [radius mass [wx wy wz]]` (x, y, z in U, z = the ball's bottom over the playfield; m/s; default radius 25 U, mass 1), `destroy NAME`, `drop TARGET 1|0`, `flip FLIPPER 1|0`,
  `enabled FLIPPER|TRIGGER 1|0`, `kick KICKER angle speed`, `sling RUBBER`, `magnet NAME V`, `turntable NAME V`, `respawn Y`, `respawn_z Z`, `speed_limit V` (default 8 m/s), `floor Z` (default -200 U), `record N`, `timing detail|total`, `end`.
- Checks every cycle: over the speed limit, below the floor, outside the table, removed by the table, skipped cycles. Unknown actions are a WARN and a notification. It resets the random seed at table start.
- Writes into VPX_TEST_OUT: balls.csv (t_ms, ball, position U, velocity U/10 ms, angular momentum), events.csv (t_ms, trigger, ball, Hit/Unhit), cycles.csv (cost split, us), frames.csv,
  report.txt (pass/fail, cost from 1 s: mean, p99, worst cycle, worst 7 cycles, 5 worst windows by part), done.txt.
  Time 0 is the first cycle after start-up (first frame under 100 ms). A restarted table writes balls_2.csv, ...
- A start spot inside or touching an object is a broken test (DevelopmentRules 4): check start spots against the table layout and the first 5 ms of each shot.
- A gate or spinner plate swings 3-8 s: space shots 5 s (gate) / 8 s (spinner).

## 6. The user's play and reports
- Launcher: `tables\VPX play - NEW BGFX.bat` (BGFX exe, prefs\rubber-contact, VPX_CONTACT_TRACE=1, copies the log to `tables\logs\vpinball_BGFX_<stamp>.log`). It does not set the engine: the table's own ini needs `[Features] unified = 1`.
- `prefs\rubber-contact\vpinball.log` is shared with every harness run and overwritten: copy it before any test when diagnosing a user report.
- Order for a report (DevelopmentRules 1e): the user's log and ini first, which model and values ran, then reproduce on that table unchanged. Crash: crash.txt and crash.dmp next to the table; the top frame names the function.
## 7. Traps
- Anti-hung time skip: a frame over 200 ms makes VPX skip simulated time. Off under VPX_PLAY_UNFOCUSED and on the bench; a focused run with a slow script still skips. The driver's report counts skipped cycles: must be 0.
- Unfocused frames are slow, so a table script's 1 ms timer fires once per frame: never take rebounds or timings from a script's own log; use the driver's balls.csv.
- The trace log rotates at ~1.7 MB (~the last 1 s at VPX_CONTACT_TRACE=2): end the drive ~1 s after the contact; check balls.csv matches the full run.
- VBScript is case-insensitive: a variable and a Sub with the same name fail to compile with no log error and no output.
- PowerShell: a function named `R` is shadowed by the alias `r`; an array of arrays is built with `@(,@(...))` and `+= ,@(...)`; pass table and ini to run.ps1 by hashtable splat, not array.
- Offline C++ benches from Bash: `timeout 250 cmd //c '<full path>\b3.bat <name>'`; PowerShell Start-Process -Wait can hang on the compiler's child processes. Hinge benches need `/FIcassert`.
- `dumpvpx.py` (olefile, runs\2026-10-31-parts-check) reads object positions and rotations from a .vpx.
- Bench group book t = balls.csv t + 4004 ms; log t = drive t + 3003 ms. Drives, stock runs, overlap tests: runs\2026-11-15-starts-fix\README.md.

## 8. Environment switches
Every switch the engine reads, one line each. A switch added, changed or removed in the code is changed here in the same commit. Settings that are ini keys are described beside the key in the code.

**Harness** (run mode and input; the laws are unchanged)
```
VPX_PHYSICS_INI=<file>     use this physics ini instead of <table>.physics.ini / vpx_physics.ini next to the table. The engine switch is the ini's [Features] unified key.
VPX_UNIFIED=0|1            run.ps1's own variable (default 1), written into the run's ini; the engine does not read it.
VPX_PLAY_UNFOCUSED=1       the game keeps running without window focus; also turns off the anti-hung time skip.
VPX_TEST_INPUT=<file>      test driver input (section 5).            VPX_TEST_OUT=<folder>   where the driver writes.
VPX_TEST_IMPULSES=1        the driver also writes impulses.csv: per ball and source (contact surface, gravity, stock object, driver), cumulative J and L in SI.
VPX_BENCH=1                the real bench: no frames, no input, virtual 7 ms frames, physics as fast as the CPU goes, one run per process; acts as VPX_PLAY_UNFOCUSED.
VPX_BENCH_PREROLL=<ms>     bench physics time before the driver's time 0 (default 3000).     VPX_BENCH_MS=<ms>   stop at this time on a run with no input file.
```

**Diagnostics** (log or file only; a ball's bytes do not change)
```
VPX_SPEED_WATCH=<m/s>      fast-ball warning with the ball's last 8 cycles and every other ball (default 8, 0 = off; at most once a second per ball).
VPX_CONTACT_TRACE=1|2      contact engage / release / hit lines (1); also per-cycle lines (2, huge).
VPX_CONTACT_LOG=<csv>      every stock collision and resting contact of the first ball: object type and name (stock engine, unified = 0, only; the playfield plane is named after the table, type 0).
VPX_CONTACT_NEAR=x,y,r     rubber faces within r of (x, y), U, when the surface is built.
VPX_CONTACT_TIMING=<csv>   per-cycle cost of the contact system by stage. One file per process.
VPX_CONTACT_PROFILE=<n>    scoped timers, a log table every n step calls (default 5000).     VPX_CHECK_PROFILE=1   one check's time split ("Check profile"), cut, tube counters ("Finder P4b/P5 totals").
VPX_PHYS_TIMING=<csv>      frame and iteration timing rows past a threshold (us): VPX_PHYS_TIMING_FRAME_US, _GAP_US, _ITER_US.
VPX_RUBBERS_LOG=1          a RubberSolve line when a rubber contact ends: in / out velocity, spin, dwell, peak force.
VPX_BAND_MESH_LOG=<csv>    band ring centres before and after deformation.
VPX_BUMPER_TRACE=<csv>     one row per cycle for each bumper not at rest.
VPX_CORD_TRACE=<csv>       per-cycle group energy and cord-ball state.     VPX_CORD_ENERGY=1   the largest per-cycle energy rise of a group with a cord.
VPX_WIRE_GAP=1             probe gap and hold state of balls near a wire switch (first 600 lines).
VPX_TRACE_BALL=<name>      the driver ball's position, velocity, contact items and every stock Collide / Contact, per stage of the cycle.
VPX_GROUP_BOOK=<csv>       one row per group per cycle: steps, rounds, halvings, stop steps, energy steps, failures, energy books, cost.
VPX_GROUP_BOOK_STEPS=1|2   with GROUP_BOOK: energy rise per step (1); also every stored energy term at cycle start and end (2, slow).
VPX_PLUNGER_TRACE=<file>   the cabinet plunger sensor per ms while active: raw reading, last report time, filtered position and speed.
VPX_GROUP_STATS=<file>     run totals (cycles, us/cycle, steps, rounds, halvings, energy steps, failures) every 200 ms of table time; line 2: steps alone and sleeps per part kind.
VPX_GROUP_PHASE_TIMING=<csv>  per-step phase timing and counts (attempts, evals, linear solves, halvings, capsules, terms).
VPX_STEP_RULE=<csv>        one row per group step: what set the step, h, solve rounds.
VPX_TOUCH_TRACE=<csv>      touch_predict: one row per touch, stop or seat end check (gap, v, a, d, step go-backs so far, h, end gap, outcome past/dip/taken/short/open/engage/net_*, kind).
VPX_GROUP_CLAMP_TRACE=<file>  every Newton round of steps with a wire contact.
VPX_GROUP_CONTACT_TRACE=1  per-round depth, force, friction (first 3 energy steps).     VPX_GROUP_HISTORY_TRACE=1   last 50 steps, at the first energy-step warning.
VPX_GROUP_ROUND_TRACE=1    per-round Newton trace with probe notes (group and hinge solves), printed only with an energy-step warning.
VPX_GROUP_SEARCH_TRACE=1   the failing stop-search solve, per round, to stderr.     VPX_GROUP_STOP_TRACE=1   the first 12 stop checks of a cycle.
VPX_GROUP_JAC_MODE=1       also a central-difference Jacobian, logged; changes last digits.  VPX_GROUP_JAC_WINDOW=t0-t1  mode 1 or 2 only in that ms (engine) window; else as mode 0.
VPX_SKIRT_ROW_TRACE=1      with ROUND_TRACE: bumper skirt rows by term per round (m/s), tilt move, each chord's change and mismatch.
VPX_FINDER_CHECK=1         finder vs the full scan: missed pieces within reach, face-table mismatches, playfield cut lists, static tube readers (log, screen).
VPX_RELOOK_TRACE=1         energy-reach re-look: each bound term, the move, solver trial point or accepted position; no cap.
VPX_PATCH_CHECK=1          per-patch checks: patch lists vs merged, cut B faces vs each alone, held friction vs book (HELDCHK, KINDC), B-lin vs probe (LEV1CHK).
VPX_PATCH_TRACE=t0-t1      engine ms: per Eval each friction group (PTEV), start patches not in memory (PTMISS); 400 lines.
VPX_PAIR_CHECK=1           ball pairs: the every-pair test beside the sorted list; counts any pair the list missed (log).
VPX_HINGE_NEAR_ALL=1       every near-piece list takes all pieces (A/B of the short list, same result expected).
VPX_SNAPSHOT_TEST=1        after each contact step: restore the snapshot, step again, log any difference (=2, =3: see below).
```

**Change the result** (test only; each is removed or made an ini key: to-do 1.32)
```
VPX_GROUP_HFIX=<us>        every group step this long; 1 = the accuracy bench's reference. A step may halve to 1 ns and a cycle may pass 3000 steps.
VPX_GROUP_JAC_MODE=2       solves with the difference Jacobian instead of the analytic one.
VPX_SNAPSHOT_TEST=2        snapshot test control: parts not put back (=3: the per-patch memory the step left kept); a difference must show where it changed.
VPX_STEP_SCALE=<k>         every accuracy step limit times k (epsRest k^2); bench convergence k 0.5 / 0.25. k >= 0.25; HFIX wins.
VPX_VISCOELASTIC=0         hard-surface damping by the old alpha law.
VPX_OVERLAP_REPLAY=0       no overlap replay.        VPX_PREVENT=0   no overlap prevention (overlap standard point 1).
```
