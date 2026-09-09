#!/usr/bin/env bash
# run_mass_con_plots.sh — run ex_plot_mass_distribution_connected.C on every ROOT
# file in a folder (default: ./mass_con), producing one PDF per file in that same
# folder.  Then, if a mixed-event folder exists (default: ./mass_me), run
# ex_plot_mass_distribution_mixed_events.C on every ROOT file there too; a missing
# folder or no ROOT files is not an error — it is skipped.
#
# Usage (from the directory where the generator was run, e.g. build/work):
#   bash run_mass_con_plots.sh [folder] [me_folder]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MACRO="$SCRIPT_DIR/ex_plot_mass_distribution_connected.C"
ME_MACRO="$SCRIPT_DIR/ex_plot_mass_distribution_mixed_events.C"
FOLDER="${1:-mass_con}"
ME_FOLDER="${2:-mass_me}"

if [[ ! -f "$MACRO" ]]; then
  echo "error: cannot find $MACRO" >&2
  exit 1
fi
if [[ ! -d "$FOLDER" ]]; then
  echo "error: folder '$FOLDER' does not exist" >&2
  exit 1
fi

shopt -s nullglob
roots=("$FOLDER"/*.root)
if (( ${#roots[@]} == 0 )); then
  echo "error: no .root files in '$FOLDER'" >&2
  exit 1
fi

for f in "${roots[@]}"; do
  echo "==> plotting $f"
  root -l -b -q "${MACRO}(\"${f}\")"
done

echo "done: ${#roots[@]} file(s) processed, PDFs written to '$FOLDER'"

# Mixed-event pass — optional. Skip silently (no error) if the macro, the folder,
# or its ROOT files are missing.
me_roots=()
if [[ -f "$ME_MACRO" && -d "$ME_FOLDER" ]]; then
  me_roots=("$ME_FOLDER"/*.root)
fi
if (( ${#me_roots[@]} == 0 )); then
  echo "note: no mixed-event ROOT files (looked in '$ME_FOLDER') — skipping"
else
  for f in "${me_roots[@]}"; do
    echo "==> plotting $f"
    root -l -b -q "${ME_MACRO}(\"${f}\")"
  done
  echo "done: ${#me_roots[@]} mixed-event file(s) processed, PDFs written to '$ME_FOLDER'"
fi
