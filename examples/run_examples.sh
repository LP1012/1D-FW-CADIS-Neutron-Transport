#!/bin/bash

set -euo pipefail

PROJECT_ROOT="$(git rev-parse --show-toplevel)"
BUILD_DIR="$PROJECT_ROOT/build"
SCRIPTS_DIR="PROJECT_ROOT/scripts"
EXEC="$BUILD_DIR"/1D-FW-CADIS

# Directory containing XML files
XML_DIR="." # change this if needed

# Loop over all .xml files in the directory
for xml_file in "$XML_DIR"/*.xml; do
	# Get the base filename without extension
	base_name=$(basename "$xml_file" .xml)

	echo ""
	echo "Running MCSlab on $xml_file ..."
	"$EXEC $1"

	echo "Plotting results ..."
	python3 "$SCRIPTS_DIR/plot_collision.py $base_name"
	python3 "$SCRIPTS_DIR/plot_pathlength $base_name"

	echo "Done with $xml_file"
	echo "---------------------------"
	echo ""
done
