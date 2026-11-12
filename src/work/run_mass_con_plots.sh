#!/usr/bin/env bash
# run_mass_con_plots.sh — batch plotting for the mass / Delta-phi examples.
#
# Per-folder mode (default): plot every ROOT file in one or more output folders,
# producing one PDF per file next to it.  The plotting macro is chosen per file
# from its filename prefix (the default output names of the generators):
#   phi_con_*   -> ex_plot_phi_distribution_connected.C
#   mass_con_*  -> ex_plot_mass_distribution_connected.C
#   mass_me_*   -> ex_plot_mass_distribution_mixed_events.C
# Files with any other prefix are skipped.  A missing folder is not an error; the
# script fails only if no file at all was plotted.
#
# Compare mode (--compare): overlay several run directories, once for mass_con and
# once for phi_con.  For each run the <kind>/<kind>_*_<ymax>_<deltay>_<events>.root
# file with the most events is used (or exactly <events> with --events N); runs
# without one are skipped.  Three plots are
# made per kind, written to the current folder:
#   full_compare_<mass|phi>_<ymax>_<deltay>.pdf     — uncut hMass  (ex_plot_compare_full.C)
#   conn_compare_<mass|phi>_<ymax>_<deltay>.pdf     — hMassConn    (ex_plot_compare_conn.C)
#   connwin_compare_<mass|phi>_<ymax>_<deltay>.pdf  — hMassConnWin (ex_plot_compare_conn.C),
#     normalized per pid1 inside |y| < Ymax/2; only files that contain it are used,
#     and the plot is skipped if none do (older outputs predate this histogram).
#
# Usage:
#   # from inside a run folder, e.g. build/work/hrg_run
#   bash ../run_mass_con_plots.sh [folder ...]        # default: mass_con mass_me phi_con
#   # from the folder holding the run folders, e.g. build/work
#   bash run_mass_con_plots.sh --compare [--ymax Y] [--deltay D] [--events N] [run_dir ...]
#                                                     # default: --ymax 3 --deltay 1.5,
#                                                     # every *_run*/ directory

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
shopt -s nullglob

# ---------------------------------------------------------------------------
# Compare mode
# ---------------------------------------------------------------------------
if [[ "${1:-}" == "--compare" ]]; then
  shift
  YMAX=3
  DELTAY=1.5
  EVENTS="*"
  runs=()
  while (( $# > 0 )); do
    case "$1" in
      --ymax)   YMAX="$2";   shift 2 ;;
      --deltay) DELTAY="$2"; shift 2 ;;
      --events) EVENTS="$2"; shift 2 ;;
      *)        runs+=("${1%/}"); shift ;;
    esac
  done
  if (( ${#runs[@]} == 0 )); then
    for d in *_run*/; do runs+=("${d%/}"); done
  fi
  if (( ${#runs[@]} == 0 )); then
    echo "error: no run directories given and no *_run*/ directories here" >&2
    exit 1
  fi

  FULL_MACRO="$SCRIPT_DIR/ex_plot_compare_full.C"
  CONN_MACRO="$SCRIPT_DIR/ex_plot_compare_conn.C"
  for m in "$FULL_MACRO" "$CONN_MACRO"; do
    if [[ ! -f "$m" ]]; then
      echo "error: cannot find $m" >&2
      exit 1
    fi
  done

  nPlots=0
  for kind in mass_con phi_con; do
    files=()
    for run in "${runs[@]}"; do
      best=""
      bestN=-1
      for f in "$run/$kind/${kind}_"*"_${YMAX}_${DELTAY}_"$EVENTS.root; do
        n="$(basename "$f" .root)"
        n="${n##*_}"
        [[ "$n" =~ ^[0-9]+$ ]] || continue
        if (( n > bestN )); then best="$f"; bestN=$n; fi
      done
      if [[ -z "$best" ]]; then
        echo "note: $run has no $kind file for Ymax=$YMAX, DeltaY=$DELTAY"
      else
        files+=("$best")
      fi
    done

    if (( ${#files[@]} == 0 )); then
      echo "note: no $kind files to compare — skipping"
      continue
    fi
    list="$(IFS=,; echo "${files[*]}")"
    echo "==> comparing $kind: ${files[*]}"
    for which in full conn connwin; do
      out="${which}_compare_${kind%_con}_${YMAX}_${DELTAY}.pdf"
      case "$which" in
        full)    call="${FULL_MACRO}(\"${list}\",\"\",\"${out}\")" ;;
        conn)    call="${CONN_MACRO}(\"${list}\",\"\",\"${out}\")" ;;
        connwin) call="${CONN_MACRO}(\"${list}\",\"\",\"${out}\",\"hMassConnWin\")" ;;
      esac
      rm -f "$out"
      root -l -b -q "$call"
      if [[ -f "$out" ]]; then
        nPlots=$((nPlots + 1))
      else
        echo "note: $out not written (no input file has the histogram)"
      fi
    done
  done

  if (( nPlots == 0 )); then
    echo "error: nothing to compare (runs: ${runs[*]})" >&2
    exit 1
  fi
  echo "done: $nPlots comparison plot(s) written"
  exit 0
fi

# ---------------------------------------------------------------------------
# Per-folder mode
# ---------------------------------------------------------------------------
if (( $# == 0 )); then
  set -- mass_con mass_me phi_con
fi

nPlotted=0

for folder in "$@"; do
  if [[ ! -d "$folder" ]]; then
    echo "note: folder '$folder' does not exist — skipping"
    continue
  fi

  roots=("$folder"/*.root)
  if (( ${#roots[@]} == 0 )); then
    echo "note: no .root files in '$folder' — skipping"
    continue
  fi

  for f in "${roots[@]}"; do
    case "$(basename "$f")" in
      phi_con_*)  macro="ex_plot_phi_distribution_connected.C" ;;
      mass_con_*) macro="ex_plot_mass_distribution_connected.C" ;;
      mass_me_*)  macro="ex_plot_mass_distribution_mixed_events.C" ;;
      *)
        echo "skip: $f (unknown prefix)"
        continue
        ;;
    esac
    if [[ ! -f "$SCRIPT_DIR/$macro" ]]; then
      echo "skip: $f (cannot find $SCRIPT_DIR/$macro)"
      continue
    fi
    echo "==> plotting $f with $macro"
    root -l -b -q "$SCRIPT_DIR/${macro}(\"${f}\")"
    nPlotted=$((nPlotted + 1))
  done
done

if (( nPlotted == 0 )); then
  echo "error: no ROOT files were plotted (folders: $*)" >&2
  exit 1
fi
echo "done: $nPlotted file(s) plotted"
