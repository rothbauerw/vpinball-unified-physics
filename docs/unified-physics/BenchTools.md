# Bench tools: how to run them

Detail for two tools the testing playbook (TestingPlaybook.md) points to. Current practice only: when a practice changes, edit its line in place; the reason goes in the commit message.

## Accuracy bench (TestingPlaybook section 3)
- Folder `runs\2026-11-01-roll-smooth`. 750 real hits (stress p1-p3, Gorgar), each replayed alone: `run_replay.py <tag> <exe> <mode>`; the reference is the same with `-env VPX_GROUP_HFIX=1` (1 us steps);
  runs with `-env VPX_TEST_IMPULSES=1`.
- `judge_hit.py <ref> <ref 0.5 us> <build> <build k 0.5> <build k 0.25>` (VPX_STEP_SCALE) judges each hit by the impulse its surface gave (bounce, slide, spin; runs\2026-11-10-hit-judge\design.md):
  noise floor from the 0.5 us twin, step error or fault by convergence, and pass = the exit direction error moves the ball under 1 mm by its next contact.
- `check_imp.py <run out folder>` checks the log adds up to each ball's momentum change.

## VPX driver run, run.ps1 (TestingPlaybook section 4)
- `powershell -ExecutionPolicy Bypass -File runs\2026-10-09-slow-contact\run.ps1 -Name N -Drive <full path> -Table <full> -Ini <full>`, every path full, no helper function wrapping it.
- run.ps1 defaults to the BGFX exe. It refuses a run without `-Drive`; a run that cannot use the driver (load check, restart, real input, the table script acting on the game)
  gives `-NoDriver "<reason>"`, kept in `run_N\nodriver.txt` and `runs\nodriver_runs.txt`.
- `-Vbs` takes a FILE NAME inside `runs\2026-10-09-slow-contact` (a full path breaks it); with `-Drive` the script only logs events.
- Other options: `-NoTrace` (trace is on by default, a huge log), `-Mass`, `-TimeoutS`, `-Seed`, `-Focused`, `-Timing`.
