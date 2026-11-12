# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

THERMINATOR 2 (this fork: **hboltz**) is a Monte Carlo event generator for statistical thermal
models of heavy-ion collisions. Reference: [arXiv:1102.0273](https://arxiv.org/abs/1102.0273).
This fork modernizes the build system and C++ style while preserving the original physics.

## Build

C++17, CMake (>=3.14), ROOT 6 required. **Out-of-source build is mandatory, and the CMake source
root is `src/`, not the repo root:**

```bash
mkdir -p build && cd build
cmake -DCMAKE_CXX_COMPILER=clang++ ../src
make -j$(nproc)
```

- If ROOT isn't in default search paths: add `-DROOTSYS=/path/to/root` (or export `$ROOTSYS`
  beforehand).
- Build artifacts: `build/src/libtherminator.{so,dylib}`, `build/work/therm2_events.exe`, plus one
  `.exe` per example program in `src/work/`.
- The build generates `build/compile_commands.json`. Symlink it to the repo root for clangd:
  `ln -sf build/compile_commands.json compile_commands.json`

## Formatting

- Config at `src/.clang-format` (LLVM base, 2-space indent, 100-col limit).
- Run: `clang-format -i src/include/*.h src/src/*.cxx src/work/*.cxx`

## No CI or automated tests

Verify changes by rebuilding (`make -j$(nproc)`) and running `therm2_events.exe` or the relevant
`ex_*.exe` example against the template ini files (see below).

## Runtime prerequisites

- The SHARE particle database (`share/particles.data`, `share/decays.data`) must be on disk.
- Every binary resolves it in this order: the `ShareDir` key in `events.ini`, then the
  `$THERMINATOR_SHARE` environment variable. If neither is set, the program exits with an error.
  ```bash
  export THERMINATOR_SHARE=$(pwd)/share
  ```

## Running the event generator

Each run needs its own working directory with ini files copied from `ini/` and an empty `events/`
subdirectory:

```bash
mkdir myrun && cd myrun
cp ../ini/events.ini .
cp ../ini/blastwave.ini .
mkdir events
# Edit events.ini → FreezeOutModel=BlastWave, FreezeOutModelINI=blastwave.ini
../build/work/therm2_events.exe events.ini
```

Output: `events/event000.root`, `event001.root`, …

`events.ini` required keys:

| Key | Description |
|---|---|
| `FreezeOutModel` | `BlastWave`, `BlastWaveAniso` (anisotropic, Retiere & Lisa) or `HRG` |
| `FreezeOutModelINI` | path to the model parameter ini file (unused by `HRG`) |
| `NumberOfEvents` | number of events to generate |
| `IntegrateSamples` | Monte-Carlo samples for the Cooper-Frye integration |
| `EventDir` | output directory for ROOT files (must already exist) |
| `EventFileType` | `root`, `text`, or `root&text` |

Optional: `ShareDir` (overrides `$THERMINATOR_SHARE`).

Some `work/` examples (e.g. `ex_o4_pion_spectra.exe`) take their own CLI flags instead of an ini
file (`--share`, `--events`, `--output`) and cache Cooper-Frye integration results in a
`fmultiplicity_<hash>.txt` file next to where they're run — delete it to force recomputation after
changing model parameters.

## Architecture

The core pipeline, shared by `therm2_events.cxx` and every `src/work/ex_*.cxx` example:

```
Configurator  →        Model         →    Integrator     →   EventGenerator
(parses .ini      (abstract; GetIntegrand()   (Monte-Carlo         (drives the loop: generate
 key=value        implements the Cooper-Frye   integrates the       primordial particles per
 files)           formula on a freeze-out      integrand per         type, decay resonances,
                  hypersurface — concrete      particle type;         assemble an Event, write
                  models: Model_BlastWave,     caches per-type         ROOT/text output)
                  Model_HRG)                   max-integrand +
                                               multiplicity in
                                               fmultiplicity_<hash>.txt,
                                               keyed by Model::GetHash(),
                                               a CRC32 of the model
                                               parameters)
```

Supporting classes:
- **ParticleDB** — loads `share/particles.data` + `share/decays.data` into `ParticleType` records.
- **`DecayTable` / `DecayChannel` / `ParticleDecayer`** — resonance decay branching used by
  `EventGenerator` after primordial particles are generated.
- **`Particle` / `ParticleCoor` / `Event`** — per-particle state and the particle collection for one
  collision event.
- **Thermodynamics, Temperature, Pressure, Entropy, SoundVelocity, Viscosity, Chemistry, Energy** —
  equation-of-state helpers consumed by the freeze-out models.
- **Crc32** — hashing used for the multiplicity cache filename.
- **Parser** — low-level ini/XML file parsing underlying `Configurator`.

**Adding a new freeze-out model:** subclass `Model`, implement `GetIntegrand()` /
`GetHyperCubeVolume()`, add a `Model_<Name>.h`/`.cxx` pair, wire it into `therm2_events.cxx`'s
model-selection logic, and add the source file to `src/src/CMakeLists.txt`'s `SOURCES` list.

**Adding a new `work/` example program:** add `src/work/ex_name.cxx`, then add matching
`add_executable(ex_name.exe ex_name.cxx)` and `target_link_libraries(ex_name.exe therminator
${ROOT_LIBRARIES})` lines to `src/work/CMakeLists.txt` (see `ex_balance_function.exe` /
`ex_connected_correlations_2d.exe` for the current pattern). Examples that produce ROOT histograms
typically ship a companion `ex_plot_<name>.C` ROOT macro, registered via `configure_file(...)` in
the same CMakeLists so it's copied into the build directory.

### Directory layout (git-tracked)

| Path | Role |
|---|---|
| `src/include/` | Public headers |
| `src/src/` | Library sources (`.cxx`) → builds into `libtherminator` |
| `src/work/` | Executable sources: main generator + `ex_*` physics examples + `ex_plot_*.C` macros |
| `src/include/old/`, `src/src/old/`, `src/work/old/`, `ini/old/` | Legacy upstream THERMINATOR2 code (e.g. `Model_Example`, Lhyquid/Cracow/BWA models, femto/hbtfit). Committed to git but **not referenced by any CMakeLists.txt** — not part of the compiled build; treat as reference only. |
| `ini/` | Template configuration files (`events.ini`, `blastwave.ini`, `blastwaveaniso.ini`) copied into per-run working directories |
| `share/` | Read-only SHARE particle/decay data tables |

Root CMakeLists is at `src/CMakeLists.txt`, which adds subdirectories `src/src/` (library) and
`src/work/` (executables).

Note: `build/`, and any per-run working directories you create (e.g. `runBW/`, `runHRG/`), are
local build/output artifacts, not part of the tracked repo structure.


## Conventions

assume always mostly plus metric
Notation: K^\mu = (p_1^\mu+p_2^\mu), where p1 and p2 are the daughters' momenta. A third daughter will have p3 but this doesn't change K
Q^\mu = (p_1^\mu-p_2^\mu)
M mass of parent, m1, m2 etc masses of daughters
we work with constant matrix elements and no finite resonance width, unless specified otherwise

