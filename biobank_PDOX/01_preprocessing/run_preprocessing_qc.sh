#!/bin/bash

# Environment to activate
# dir_env="/omics/groups/OE0596/internal/maria/analysis/conda/env_scRNA"

# Define the path you want to search for subdirectories
dir_data="../data/crg_multi/"
analysis="preprocessing"
script="preprocessing_qc.py"

########## Pre-flight checks

# Check for required directories
if [ ! -d "$dir_data" ]; then
    echo "Error: Raw data directory $dir_data does not exist."
    exit 1
fi

# Script: exit if non-existent
if [ ! -f "$script" ]; then
    echo "Error: Script $script does not exist."
    exit 1
fi

########## Actual code

echo -e "Preprocessing and QC at sample level\n"

# List files in input directory
samples=$(find "$dir_data" -type f -printf '%f\n')
length=$(echo "$samples" | wc -l)
echo -e "Number of samples in $dir_data: $length\n"

# Iterate over all samples
for sample in $samples; do
    # Run python script
    echo "Process sample $sample"
    python $script -i ${dir_data} -f ${sample} -an ${analysis}
    echo -e "\n"
done