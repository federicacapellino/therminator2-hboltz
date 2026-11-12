#!/usr/bin/env bash
# run_con_generators.sh — run ex_mass_distribution_connected.exe and
# ex_phi_distribution_connected.exe inside several run directories with the same
# settings, so their outputs can be compared with run_mass_con_plots.sh --compare.
#
# The model is chosen per run directory: a name starting with "hrg" runs the HRG
# model; otherwise a directory with blastwaveaniso.ini runs the anisotropic
# blast-wave model with it, and one with blastwave.ini runs the blast-wave model.  Inside one directory the mass and Delta-phi generators run one
# after the other (they share the fmultiplicity_<hash>.txt cache); different
# directories run in parallel.  Outputs use the generators' default names
# (mass_con/mass_con_<model>_<ymax>_<deltay>_<events>.root, phi_con/phi_con_...),
# overwriting earlier files with the same settings.  Logs go to <run>/logs/.
#
# Usage (from the folder holding the run directories, e.g. build/work):
#   bash run_con_generators.sh [--events N] [--ymax Y] [--deltay D] [--samples N]
#                              [--share DIR] [--plot] [run_dir ...]
#   defaults: --events 500 --ymax 3 --deltay 1.5, every *_run*/ directory
#   --plot:   afterwards run run_mass_con_plots.sh --compare on the same directories,
#             using exactly the files just generated (same --events)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
shopt -s nullglob

EVENTS=500
YMAX=3
DELTAY=1.5
SAMPLES=""
SHARE=""
PLOT=0
runs=()
while (( $# > 0 )); do
  case "$1" in
    --events)  EVENTS="$2";  shift 2 ;;
    --ymax)    YMAX="$2";    shift 2 ;;
    --deltay)  DELTAY="$2";  shift 2 ;;
    --samples) SAMPLES="$2"; shift 2 ;;
    --share)   SHARE="$2";   shift 2 ;;
    --plot)    PLOT=1;       shift ;;
    -h|--help) sed -n '2,19p' "${BASH_SOURCE[0]}"; exit 0 ;;
    -*)        echo "error: unknown option $1" >&2; exit 1 ;;
    *)         runs+=("${1%/}"); shift ;;
  esac
done
if (( ${#runs[@]} == 0 )); then
  for d in *_run*/; do runs+=("${d%/}"); done
fi
if (( ${#runs[@]} == 0 )); then
  echo "error: no run directories given and no *_run*/ directories here" >&2
  exit 1
fi

# SHARE particle database: --share, then $THERMINATOR_SHARE, then the repo's share/.
if [[ -z "$SHARE" ]]; then
  if [[ -n "${THERMINATOR_SHARE:-}" ]]; then
    SHARE="$THERMINATOR_SHARE"
  elif [[ -d "$SCRIPT_DIR/../../share" ]]; then
    SHARE="$(cd "$SCRIPT_DIR/../../share" && pwd)"
  fi
fi
if [[ -z "$SHARE" || ! -f "$SHARE/particles.data" ]]; then
  echo "error: SHARE database not found (use --share DIR or export THERMINATOR_SHARE)" >&2
  exit 1
fi
SHARE="$(cd "$SHARE" && pwd)"

for kind in mass phi; do
  exe="$SCRIPT_DIR/ex_${kind}_distribution_connected.exe"
  if [[ ! -x "$exe" ]]; then
    echo "error: cannot find $exe (build it with make)" >&2
    exit 1
  fi
done

common=(--events "$EVENTS" --ymax "$YMAX" --deltay "$DELTAY" --share "$SHARE")
[[ -n "$SAMPLES" ]] && common+=(--samples "$SAMPLES")
tag="${YMAX}_${DELTAY}_${EVENTS}"
statusDir="$(mktemp -d)"
trap 'rm -rf "$statusDir"' EXIT

# Run both generators in one run directory; record "<run> <kind> ok|FAILED" lines.
run_one() {
  local run="$1"; shift
  local model=("$@")
  mkdir -p "$run/logs"
  for kind in mass phi; do
    local log="logs/${kind}_con_${tag}.log"
    if (cd "$run" && "$SCRIPT_DIR/ex_${kind}_distribution_connected.exe" "${model[@]}" \
          "${common[@]}" > "$log" 2>&1); then
      echo "$run $kind ok" >> "$statusDir/status"
      echo "    $run: $kind ok"
    else
      echo "$run $kind FAILED" >> "$statusDir/status"
      echo "    $run: $kind FAILED (see $run/$log)"
    fi
  done
}

started=()
for run in "${runs[@]}"; do
  if [[ ! -d "$run" ]]; then
    echo "skip: $run (not a directory)"
    continue
  fi
  name="$(basename "$run")"
  if [[ "$name" == hrg* ]]; then
    model=(--model hrg)
  elif [[ -f "$run/blastwaveaniso.ini" ]]; then
    model=(--model blastwaveaniso --ini blastwaveaniso.ini)
  elif [[ -f "$run/blastwave.ini" ]]; then
    model=(--model blastwave --ini blastwave.ini)
  else
    echo "skip: $run (no blastwave.ini or blastwaveaniso.ini and not hrg*)"
    continue
  fi
  echo "==> $run: started (${model[*]}, events=$EVENTS, Ymax=$YMAX, DeltaY=$DELTAY)"
  run_one "$run" "${model[@]}" &
  started+=("$run")
done

if (( ${#started[@]} == 0 )); then
  echo "error: no run directory could be started" >&2
  exit 1
fi
wait

# Summary: key numbers from each log.
echo
for run in "${started[@]}"; do
  for kind in mass phi; do
    log="$run/logs/${kind}_con_${tag}.log"
    [[ -f "$log" ]] || continue
    echo "--- $run $kind"
    grep -E "yield|count|in window" "$log" | sed 's/^/    /' || true
  done
done

failed="$(grep FAILED "$statusDir/status" || true)"
if [[ -n "$failed" ]]; then
  echo
  echo "error: some generators failed:" >&2
  echo "$failed" | sed 's/^/    /' >&2
  exit 1
fi
echo
echo "done: $(wc -l < "$statusDir/status") generator run(s) finished"

if (( PLOT )); then
  echo
  bash "$SCRIPT_DIR/run_mass_con_plots.sh" --compare --ymax "$YMAX" --deltay "$DELTAY" \
    --events "$EVENTS" "${started[@]}"
fi
