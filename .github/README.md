# Unified physics for Visual Pinball X (prototype)

> This page is the README of the **unified-physics** fork. The original Visual Pinball X README is [here](../README.md); build instructions for VPX itself are in [make/README.md](../make/README.md).

This branch (`unified-physics`) is a work-in-progress, **opt-in** rebuild of the ball contact physics for [Visual Pinball X](https://github.com/vpinball/vpinball). It sits on top of the official `vpinball/vpinball` master (`da8b594ae`, 2026-10-21, which includes Toxie's physics fixes `8c4cf05d`) and adds new physics in separate files.

**Status: prototype for discussion and review. Not a release, not endorsed by the VPX maintainers, not ready to merge.** It was developed with AI assistance (Claude Code and ChatGPT for testing); it was developed against measurements and the checks listed below, but it has not had independent review yet.

## What it does
Physical models in place of the stock impulse-and-correction contact rules, all on one contact system that owns every ball it touches:

- **Flippers:** a two-body model (coil plunger + bat/crank, link, return spring, rest stop) with real-data parameters, and a rubber contact model on the bat.
- **Rubbers, posts, bands and slings:** force-model rubber rings and stretched bands (with sling plate, coil and switches), surface rubbers on walls and primitives.
- **Surfaces:** the playfield, walls, primitives and ramps as hard surfaces with Hertz contact plus Hunt-Crossley damping; resting-contact friction, rolling resistance and twist friction.
- **Gates and spinners:** hinged plates of the same contact system.
- **Magnets:** coil + core models including guide tubes, ball-ball induced forces.
- **Turntables:** a motor-driven disc set into the playfield.

Design and rationale: `docs/UnifiedContactDesign.txt`, `docs/UnifiedContactAudit.txt`, `docs/GateSpinnerDesign.txt`. Every key and feature: `docs/NewPhysicsReference.txt`. Ground rules used for the work: `docs/DevelopmentRules.txt` (no change without measured data, no hacks, list what breaks).

## Off by default
A table uses the new physics **only** if its physics ini says so:

1. Put a file named `<table name>.physics.ini` next to the table (for `MyTable.vpx`: `MyTable.physics.ini`). A shared `vpx_physics.ini` in the table's folder is read as a testing fallback only.
2. In it, under `[Features]`, set `unified = 1`. Every model is then on with its defaults; other `[Features]` keys switch single models off, and the other sections hold per-part settings (see the reference).

No ini, or no `unified = 1`, means stock VPX physics for that table. (A stock-mode comparison against plain upstream has **not** been done yet.)

## Where the code is
- New files (the bulk of the work): `src/physics/` (`contactsystem`, `hitrubber`, `contactbat`, `hitflippermodel`, `hitsling`, `hitballcontact`, `contactsettings`, `magnet*`, `turntable`, `rubbercontact.h`, `physicsini.h`, ...), `src/parts/ramphardsurface.cpp`, `src/parts/rubberband.cpp`.
- Edits to existing upstream files are small redirects into those files; `docs/StockFilesTouched.txt` lists every stock file touched and why.
- Where upstream already has its own fix for the same problem, upstream's version is used.

## Build
Same as upstream (see `make/README.md`). Windows, x64: build the external libraries with `platforms/windows-x64/external.sh`, then
```
cmake -DRENDERER=GL -DARCH=x64 -DPLATFORM=windows -B build
cmake --build build --config Release --target vpinball
```
The GL renderer is what the checks below ran on. A BGFX build (`-DRENDERER=BGFX`) builds and starts; the drawing hooks for magnet guide tubes and turntable discs have not yet been checked on screen on BGFX. DX9 has not been tried.

## What has been checked (and what has not)
Compared against the previous branch (`rubber-contact`) on the GL build with `unified = 1`, scripted ball/flipper scenes: ball drop, resting ball, rubber post scenes, flipper reflip and spin scenes are byte-identical; a flipper-cradle scene and a ball-near-post scene show small differences (the latter vanishes when the post is also a hard surface of the contact model). Known open items:

- The one stock fix not in upstream, a plunger fix (`HitPlunger`: ignore a hit that `Collide()` would not act on, which otherwise repeats at zero time), is not re-applied yet.
- `unified = 0` against a plain-upstream build: not yet checked.
- BGFX drawing hooks: not yet checked on screen. Other BGFX backends (D3D12, OpenGL, Vulkan) not tried.
- Some parts of the stock files are still hooks rather than separated code (`PhysicsEngine.cpp` loop, `hitflipper.cpp` hit tests, the force-model setup in `Rubber::PhysicSetup`).
- Editor support (UI panels) and storing the settings in the `.vpx` are not done; settings live in the ini.

## Feedback
Issues and pull requests on this fork are welcome. For anything affecting VPX itself, discussion with the VPX maintainers comes first.

License: as upstream (GPLv3+; see `LICENSE`).
