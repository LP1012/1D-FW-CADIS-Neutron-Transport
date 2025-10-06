#!/bin/bash
# run_simulations.sh
# This script runs MCSlab on all XML files and then plots the results.

# Directory containing XML files
XML_DIR="."  # change this if needed


# Loop over all .xml files in the directory
for xml_file in "$XML_DIR"/*.xml; do
    # Get the base filename without extension
    base_name=$(basename "$xml_file" .xml)
    
    echo ""
    echo "Running MCSlab on $xml_file ..."
    ../build/MCSlab "$xml_file"
    
    # Generate corresponding output CSV filename
    csv_file="${base_name}_out.csv"
    
    echo "Plotting results to $csv_file ..."
    python3 ../scripts/plot_simulation.py "$csv_file"
    
    echo "Done with $xml_file"
    echo "---------------------------"
    echo ""
done

