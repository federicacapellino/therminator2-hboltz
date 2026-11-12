# Therminator 2 (hboltz) — Agent Notes

## Build

- C++17, CMake (>=3.14), ROOT 6 required. Out-of-source build mandatory.
- Source root is `src/`, **not** the repo root:
  ```bash
  mkdir -p build && cd build
  cmake -DCMAKE_CXX_COMPILER=clang++ ../src
  make -j$(nproc)
  ```
- If ROOT isn't in default search paths: add `-DROOTSYS=/path/to/root` to cmake.
- Artifacts: `build/src/libtherminator.{so,dylib}`, `build/work/therm2_events.exe`, plus example binaries.

## Runtime Prerequisites

- SHARE particle database (`share/particles.data`, `share/decays.data`) must be on disk.
- Export before running any binary: `export THERMINATOR_SHARE=$(pwd)/share`
- Alternatively set `ShareDir` key inside `events.ini`.

## Running the Event Generator

Each run needs its own directory with ini files copied from `ini/` and an empty `events/` subdirectory. Example setup for BlastWave model:
```bash
mkdir myrun && cd myrun
cp ../ini/events.ini .
cp ../ini/blastwave.ini .
mkdir events
# Edit events.ini → FreezeOutModel=BlastWave, FreezeOutModelINI=blastwave.ini
../build/work/therm2_events.exe events.ini
```
Output: `events/event000.root`, `event001.root`, …

## Code Structure

| Path | Role |
|---|---|
| `src/include/` | All public headers (.h) |
| `src/src/` | Library sources (.cxx) → builds into `libtherminator` |
| `src/work/` | Executable sources: main generator + examples |
| `share/` | Read-only particle/decay data tables |
| `ini/` | Template configuration files |

Root CMakeLists at `src/CMakeLists.txt`. Subdirs `src/src/` (library) and `src/work/` (executables).

## Style & Formatting

- Clang-format config: `src/.clang-format` (LLVM base, 2-space indent, 100-col limit).
- Run as: `clang-format -i src/include/*.h src/src/*.cxx src/work/*.cxx`
- Compile commands: `build/compile_commands.json`. Symlink from repo root for clangd: `ln -sf build/compile_commands.json compile_commands.json`

## No CI or Tests

Nothing automated. Verify changes by rebuilding (`make -j$(nproc)`) and running `therm2_events.exe` with the example ini files from `myrun/`.
