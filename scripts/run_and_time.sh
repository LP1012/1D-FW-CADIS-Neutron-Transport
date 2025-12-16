#!/bin/bash

# Record the start time
START_TIME=$SECONDS

# --- Your commands go here ---
echo "Executing some tasks..."
sleep 2 # Example task
echo "Tasks finished."
# -----------------------------

# Calculate the duration
DURATION=$((SECONDS - START_TIME))

echo "Execution time was $DURATION seconds."
