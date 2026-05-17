# Therminator 2 (hboltz)

THERMINATOR 2 is a Monte Carlo event generator for statistical thermal models of
heavy-ion collisions.  This fork (hboltz) modernizes the build system and C++ style
while preserving the physics.

Reference: [arXiv:1102.0273](https://arxiv.org/abs/1102.0273)

---

## Building

Requires [ROOT 6](https://root.cern) and a C++17 compiler.  An out-of-source build
is required.

```bash
mkdir build && cd build
cmake -DCMAKE_CXX_COMPILER=clang++ ../src
make -j4
```

On macOS with ROOT installed via MacPorts the cmake step finds ROOT automatically.
If ROOT is installed elsewhere, pass `-DROOTSYS=/path/to/root` or set `$ROOTSYS`
before running cmake.

The build produces:
- `build/src/libtherminator.dylib` — shared library
- `build/work/therm2_events.exe` — main event generator
- `build/work/t1.exe`, `ex_o4_pion_spectra.exe` — example programs

---

## The SHARE particle database

All runs need access to the SHARE particle and decay tables (`particles.data`,
`decays.data`).  These live in the `share/` directory of this repository and are
read-only reference data shared across every run.

Set the environment variable once in your shell profile:

```bash
export THERMINATOR_SHARE=/path/to/hboltz/share/
```

Every program looks for `ShareDir` in the run's `events.ini` first.  If that key is
absent it falls back to `$THERMINATOR_SHARE`.  If neither is set the program exits
with a clear error message.

---

## Running therm2_events

Each run lives in its own working directory.  Set one up by copying the template ini
files and creating the output directory:

```bash
mkdir myrun && cd myrun
cp /path/to/hboltz/ini/events.ini .
cp -r /path/to/hboltz/ini/fomodel .
mkdir events
```

Edit `events.ini` to choose a freeze-out model and set run parameters, then edit the
corresponding model file in `fomodel/` to set the physics parameters.

Run the generator:

```bash
/path/to/build/work/therm2_events.exe events.ini
```

Output ROOT files appear in `events/` as `event000.root`, `event001.root`, etc.

### events.ini — required keys

| Key | Description |
|---|---|
| `FreezeOutModel` | model name: `BlastWave` or `HRG` |
| `FreezeOutModelINI` | path to the model parameter ini file (not used by `HRG`) |
| `NumberOfEvents` | number of events to generate |
| `IntegrateSamples` | Monte Carlo samples for Cooper-Frye integration |
| `EventDir` | directory to write output ROOT files (must already exist) |
| `EventFileType` | `root`, `text`, or `root&text` |

Optional key: `ShareDir` (overrides `$THERMINATOR_SHARE`).

### Freeze-out models

| `FreezeOutModel` | model ini file |
|---|---|
| `BlastWave` | specified by `FreezeOutModelINI` |
| `HRG` | *(no model ini needed)* |

---

## Quick-start example (ex_o4_pion_spectra.exe)

`ex_o4_pion_spectra.exe` generates a thermal HRG gas, runs decays, and fills pion
momentum spectra separated by origin (primordial vs. decay) and decay lifetime.
It also estimates the effect of a modified O(4) pion dispersion relation near the
chiral critical point by reweighting primordial pions by `f(p)/fvac(p)`, and
estimates stimulated emission from short-lived resonance decays via a `(1 + f(p))`
enhancement factor.

Run from any directory — the share path is taken from `$THERMINATOR_SHARE` or passed
explicitly:

```bash
../build/work/ex_o4_pion_spectra.exe [--share /path/to/share] [--events 500] [--output o4_pion_spectra.root]
```

On the first run it computes the Cooper-Frye integrals and writes a multiplicity cache
file (`fmultiplicity_<hash>.txt`); subsequent runs read the cache and start
generating immediately.

Output histograms are written to `o4_pion_spectra.root`.  To plot them:

```bash
root -b -q 'build/work/ex_plot_o4_pion_spectra.C("o4_pion_spectra.root")'
```

This produces `o4_pion_spectra.pdf` showing the pion `dN/dp` distributions broken
down by type: all pions, decay pions, thermal (primordial) pions, and pions from
short-lived decays (τ < 2.5 fm/c).

---

## VSCode / clangd integration

The build generates `build/compile_commands.json`.  A symlink at the repository root
points to it so clangd finds it automatically:

```bash
ln -s build/compile_commands.json compile_commands.json
```

This eliminates IntelliSense squiggles for ROOT headers and project includes.
