#!/bin/bash
#SBATCH --job-name=merging-pipeline
#SBATCH --output=logs/merging-pipeline_%j.log
#SBATCH --error=logs/merging-pipeline_%j.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=fn_long,fn_medium,gpu8_medium,gpu8_long
#SBATCH --mem=500G

set -e
date
echo "Starting final scRNA-seq analysis"

SAMPLESHEET="samplesheet.csv"
PARAMS_FILE="params.yaml"
WORKING_DIR=$(pwd)

SCRIPTS_DIR="$WORKING_DIR/scripts"

# Load R module
module load r/4.3.2

cd $WORKING_DIR

# Run final analysis

Rscript $SCRIPTS_DIR/scRNA-processing.R "$SAMPLESHEET" "$PARAMS_FILE"

date
echo "Completed final analysis"

