#!/bin/bash
#SBATCH --job-name=scRNA-pipeline
#SBATCH --output=logs/scRNA-pipeline_%j.log
#SBATCH --error=logs/scRNA-pipeline_%j.err
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --partition=fn_long,fn_medium,gpu8_medium,gpu8_long
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G

set -e
date
echo "Starting LUAD sample processing pipeline"

# Load R module
module load r/4.3.2

# Define paths
SAMPLESHEET="samplesheet.csv"
PARAMS_FILE="params.yaml"
WORKING_DIR=$(pwd)

SCRIPTS_DIR="$WORKING_DIR/scripts"

# Check if sample sheet exists
if [ ! -f "$SAMPLESHEET" ]; then
    echo "Error: Sample sheet not found at $SAMPLESHEET"
    exit 1
fi

# Check if params file exists
if [ ! -f "$PARAMS_FILE" ]; then
    echo "Error: Parameters file not found at $PARAMS_FILE"
    exit 1
fi

# Check if scripts exist
if [ ! -f "$SCRIPTS_DIR/individual_samples.R" ]; then
    echo "Error: Process script not found at $SCRIPTS_DIR/individual_samples.R"
    exit 1
fi

if [ ! -f "$SCRIPTS_DIR/scRNA-processing.R" ]; then
    echo "Error: Final processing script not found at $SCRIPTS_DIR/scRNA-processing.R"
    exit 1
fi

# Count number of samples (excluding header)
NB_SAMPLES=$(wc -l < "$SAMPLESHEET")
NB_SAMPLES=$((NB_SAMPLES - 1))  # Subtract header line
echo "Processing $NB_SAMPLES samples from $SAMPLESHEET"

WORKING_DIR=$(pwd)
# Directory to store job IDs
JOB_DIR="$WORKING_DIR/logs/job_ids_$$"
mkdir -p $JOB_DIR

# Array to store job IDs
JOB_IDS=()

# Process each sample (skip header)
for ((n=2; n<=NB_SAMPLES+1; n++)); do
    # Extract the sample from the CSV file
    SAMPLE=$(sed -n "${n}p" "$SAMPLESHEET" | cut -d',' -f1)
    echo "Processing sample: $SAMPLE"
    
    # Create a temporary job script for this sample
    TEMP_SCRIPT="$JOB_DIR/j${n}-${SAMPLE}.sh"
    JOB_INDEX=$((n - 1))
    cat > "$TEMP_SCRIPT" << EOL
#!/bin/bash
#SBATCH --job-name=${SAMPLE}_${JOB_INDEX}
#SBATCH --output=${JOB_DIR}/${SAMPLE}_J${JOB_INDEX}_%j.log
#SBATCH --error=${JOB_DIR}/${SAMPLE}_J${JOB_INDEX}_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=fn_long,fn_medium,gpu8_medium,gpu8_long
#SBATCH --mem=72G

set -e
date
echo "Processing sample: $SAMPLE"

# Load R module
module load r/4.3.2

# Run R script for individual sample processing
cd $WORKING_DIR
Rscript $SCRIPTS_DIR/individual_samples.R "$SAMPLE" "$SAMPLESHEET" "$PARAMS_FILE"

date
echo "Completed processing sample: $SAMPLE"
EOL
    
    chmod +x "$TEMP_SCRIPT"
    
    # Submit the job and capture its ID
    JOB_ID=$(sbatch --parsable "$TEMP_SCRIPT")
    JOB_IDS+=($JOB_ID)
    echo "Submitted job $JOB_ID for sample: $SAMPLE"
    echo "$JOB_ID" >> "$JOB_DIR/all_job_ids.txt"
done

# Create dependency string from all job IDs
DEPENDENCY_LIST=$(IFS=:; echo "${JOB_IDS[*]}")


# Create the final analysis job
FINAL_SCRIPT="$JOB_DIR/final_analysis.sh"

cat > "$FINAL_SCRIPT" << EOL
#!/bin/bash
#SBATCH --job-name=merging-pipeline
#SBATCH --output=${JOB_DIR}/merging-pipeline_%j.log
#SBATCH --error=${JOB_DIR}/merging-pipeline_%j.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=fn_long,fn_medium,gpu8_medium,gpu8_long
#SBATCH --mem=500G

set -e
date
echo "Starting final scRNA-seq analysis"

# Load R module
module load r/4.3.2

cd $WORKING_DIR
# Extract output directory from params.yaml
OUTPUT_DIR=$(grep '^output:' "$PARAMS_FILE" | awk '{gsub(/"/, "", $2); print $2}')

echo "\$OUTPUT_DIR"

# Run QC metrics script
Rscript $SCRIPTS_DIR/qc_metrics.R "\$OUTPUT_DIR"

# Run final analysis

Rscript $SCRIPTS_DIR/scRNA-processing.R "$SAMPLESHEET" "$PARAMS_FILE"

date
echo "Completed final analysis"

# Cleanup
EOL

chmod +x "$FINAL_SCRIPT"

# Submit the final job with dependency on all sample jobs
FINAL_JOB_ID=$(sbatch --parsable --dependency=afterok:$DEPENDENCY_LIST "$FINAL_SCRIPT")
echo "Submitted final analysis job $FINAL_JOB_ID that will run after all sample processing jobs complete"

echo "All jobs submitted. Monitor progress with 'squeue -u \$USER'"
date
