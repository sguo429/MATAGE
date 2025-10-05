#!/bin/bash

#
# MATAGE Workflow Runner
#
# This script orchestrates the MATAGE analysis pipeline by first fitting a null model
# and then processing a genotype file in parallel chunks. All parameters are
# controlled via command-line flags.
#

# --- 1. Usage Function ---
# Displays help information and exits.
usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

This script runs the MATAGE workflow.

Required Arguments:
  --pheno_infile      Path to the phenotype data file.
  --geno_infile       Path to the genotype data file.
  --outfile_dir       Directory to save result files.
  --outfile_name      Basename for output files.
  --ID_name           Column name for subject IDs in the phenotype file.
  --formula           Model formula for the null model (e.g., "Phenotype~Age+Sex"). Must be quoted.
  --family            Model family type: 'linear' or 'logistic'.
  --geno_start_col    Integer specifying the starting column of genotype data in geno_infile.

Optional Arguments:
  --kins_infile       Path to the kinship matrix file. (Default: NULL)
  --total_rows        Total number of variants (rows) to process in geno_infile. (Default: line count of geno_infile, minus header)
  --num_rows          Number of rows to process per job. (Default: same as --total_rows)
  --max_jobs          Maximum number of parallel Rscript jobs to run. (Default: 1)
  -h, --help          Display this help message and exit.

EOF
    exit 1
}

# --- 2. Initialize Variables and Set Defaults ---
kins_infile="NULL"
num_rows=""
max_jobs=1
family_r_arg=""
outfile_name=""

# --- 3. Parse Command-Line Arguments ---
if [ "$#" -eq 0 ]; then
    usage
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --pheno_infile) pheno_infile="$2"; shift ;;
        --geno_infile) geno_infile="$2"; shift ;;
        --kins_infile) kins_infile="$2"; shift ;;
        --outfile_dir) outfile_dir="$2"; shift ;;
        --outfile_name) outfile_name="$2"; shift ;;
        --ID_name) ID_name="$2"; shift ;;
        --formula) formula="$2"; shift ;;
        --family) family="$2"; shift ;;
        --geno_start_col) geno_start_col="$2"; shift ;;
        --total_rows) total_rows="$2"; shift ;;
        --num_rows) num_rows="$2"; shift ;;
        --max_jobs) max_jobs="$2"; shift ;;
        -h|--help) usage ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# --- 4. Validate Inputs and Finalize Variables ---

# Check for all required arguments
if [ -z "$pheno_infile" ] || [ -z "$geno_infile" ] || [ -z "$outfile_dir" ] || [ -z "$outfile_name" ] || [ -z "$ID_name" ] || [ -z "$formula" ] || [ -z "$family" ] || [ -z "$geno_start_col" ]; then
    echo "Error: Missing one or more required arguments."
    usage
fi

# Set default for total_rows if not provided by counting lines in genotype file (minus header)
if [ -z "$total_rows" ]; then
    echo "Info: --total_rows not provided. Calculating from genotype file..."
    if [ ! -r "$geno_infile" ]; then
        echo "Error: --geno_infile (${geno_infile}) is not readable or does not exist. Cannot determine --total_rows."
        exit 1
    fi
    total_rows=$(($(wc -l < "$geno_infile") - 1))
    echo "Info: Calculated --total_rows to be ${total_rows}."
fi

# Set default for num_rows if not provided
if [ -z "$num_rows" ]; then
    num_rows=$total_rows
fi

# Validate that numeric inputs are integers
int_regex='^[0-9]+$'
if ! [[ $geno_start_col =~ $int_regex ]] || ! [[ $total_rows =~ $int_regex ]] || ! [[ $num_rows =~ $int_regex ]] || ! [[ $max_jobs =~ $int_regex ]]; then
    echo "Error: --geno_start_col, --total_rows, --num_rows, and --max_jobs must all be integers."
    exit 1
fi

# Translate --family argument to R syntax
if [ "$family" == "linear" ]; then
    family_r_arg="gaussian(link=\"identity\")"
elif [ "$family" == "logistic" ]; then
    family_r_arg="binomial(link = \"logit\")"
else
    echo "Error: --family must be either 'linear' or 'logistic'."
    usage
fi

# Get the absolute path to the directory containing this script
# This allows the R scripts to be found reliably
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
R_SCRIPT_DIR="${SCRIPT_DIR}"

# --- 5. Define Helper Function ---
# Function to check the number of running background jobs and wait if it exceeds max_jobs
wait_for_jobs() {
    while [ $(jobs -p | wc -l) -ge "$max_jobs" ]; do
        sleep 1
    done
}

# --- 6. Execute Workflow ---
echo "Starting MATAGE workflow..."

# Define and create output subdirectories
null_model_dir="${outfile_dir}/null_model_files"
results_subdir="${outfile_dir}/result_files"
log_dir="${outfile_dir}/log_files"
mkdir -p "$null_model_dir" "$results_subdir" "$log_dir"

# Construct the path for the null model file automatically
null_model_file="${null_model_dir}/${outfile_name}_null_model.RData"
fit_null_log_file="${log_dir}/${outfile_name}_fit_null_model.log"

echo "Step 1: Fitting the null model..."
echo "  - Null model will be saved to: ${null_model_file}"
echo "  - Log file will be saved to: ${fit_null_log_file}"

# Run the null model fitting script and capture time/memory usage to a log file.
# The command is wrapped in parentheses so that redirection captures the output of 'time'.
( /usr/bin/time -v Rscript "${R_SCRIPT_DIR}/Fit_null_model.R" \
    "$pheno_infile" \
    "$kins_infile" \
    "$ID_name" \
    "$formula" \
    "$family_r_arg" \
    "$null_model_file" ) > "$fit_null_log_file" 2>&1

if [ $? -ne 0 ]; then
    echo "Error: Null model fitting failed. Check the log for details: ${fit_null_log_file}"
    exit 1
fi

echo "Step 2: Processing genotype file in chunks..."
# Loop through the start_row values (from 1 to total_rows with steps of num_rows)
for start_row in $(seq 1 "$num_rows" "$total_rows"); do
    wait_for_jobs

    end_row=$((start_row + num_rows - 1))
    # Ensure end_row does not exceed total_rows
    if [ "$end_row" -gt "$total_rows" ]; then
        end_row=$total_rows
    fi

    outfile="${results_subdir}/${outfile_name}_${start_row}_${end_row}.txt"
    job_log_file="${log_dir}/${outfile_name}_${start_row}_${end_row}.log"
    echo "  - Dispatching job for rows ${start_row}-${end_row}. Log: ${job_log_file}"

    # Run the main analysis R script in the background, capturing time/memory to a log.
    ( /usr/bin/time -v Rscript "${R_SCRIPT_DIR}/MATAGE.R" \
        "$geno_infile" \
        "$start_row" \
        "$num_rows" \
        "$outfile" \
        "$null_model_file" \
        "$geno_start_col" ) > "$job_log_file" 2>&1 &

done

# Wait for all background jobs to complete
echo "Waiting for all dispatched jobs to finish..."
wait
echo "MATAGE workflow completed successfully."
