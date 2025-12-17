#!/bin/bash
set -euo pipefail

PROJECT_ROOT="$(git rev-parse --show-toplevel)"
BUILD_DIR="$PROJECT_ROOT/build"
SCRIPTS_DIR="PROJECT_ROOT/scripts"
EXEC="$BUILD_DIR"/1D-FW-CADIS

# Record the start time
START_TIME=$SECONDS

# --- Your commands go here ---
echo "Executing some tasks..."

"$EXEC $1"

echo "Calculating runtime..."
# Calculate the duration
DURATION=$((SECONDS - START_TIME))
echo "Done."

echo "Beginning post-processing..."
base_name=$(basename "$1" .xml)
python3 "$SCRIPTS_DIR/plot_pathlength $base_name"

echo "Post-processing complete."
# -----------------------------

echo ""
echo "Execution time was $DURATION seconds."
