#!/bin/bash
#PBS -N jobname
#PBS -l nodes=node:ppn=cores
#PBS -l mem=ram
#PBS -o job.log

# Self-submission: if not inside a PBS job, parse flags and submit via qsub
if [ -z "$PBS_JOBID" ]; then
  JOB_NAME="numders-fw"
  NODE="nv81"
  PPN=32
  FILE_LIST="FILES"
  CORES_PER_CALC=8
  MAX_PARALLEL_CALCS=4

  while getopts "N:n:p:f:c:m:h" opt; do
    case $opt in
      N) JOB_NAME="$OPTARG" ;;         # Job name          (-N myjob)
      n) NODE="$OPTARG" ;;             # Node name         (-n nv82)
      p) PPN="$OPTARG" ;;             # Cores/node        (-p 16)
      f) FILE_LIST="$OPTARG" ;;        # File list         (-f MY_FILES)
      c) CORES_PER_CALC="$OPTARG" ;;   # Cores per calc    (-c 4)
      m) MAX_PARALLEL_CALCS="$OPTARG" ;; # Max parallel    (-m 8)
      h) echo "Usage: $0 [-N job_name] [-n node] [-p ppn] [-f file_list] [-c cores_per_calc] [-m max_parallel_calcs]"
         echo ""
         echo "Options:"
         echo "  -N  Job name (default: numders-fw)"
         echo "  -n  Node name (default: nv81)"
         echo "  -p  Cores per node / ppn (default: 32)"
         echo "  -f  File containing list of input .com files (default: FILES)"
         echo "  -c  Cores per Gaussian calculation (default: 8)"
         echo "  -m  Max parallel calculations (default: 4)"
         exit 0 ;;
      *) echo "Usage: $0 [-N job_name] [-n node] [-p ppn] [-f file_list] [-c cores_per_calc] [-m max_parallel_calcs]"
         echo "Run '$0 -h' for more information."
         exit 1 ;;
    esac
  done

  qsub -N "$JOB_NAME" \
       -l "nodes=${NODE}:ppn=${PPN}" \
       -v "FILE_LIST=$FILE_LIST,CORES_PER_CALC=$CORES_PER_CALC,MAX_PARALLEL_CALCS=$MAX_PARALLEL_CALCS" \
       "$0"
  exit $?
fi

# User-defined parameters (may be passed via -v from self-submission)
FILE_LIST="${FILE_LIST:-FILES}"
CORES_PER_CALC="${CORES_PER_CALC:-8}"
MAX_PARALLEL_CALCS="${MAX_PARALLEL_CALCS:-4}"

# Create a scratch directory for the job
if [ ! -d /scratch/users/$USER ]; then
  mkdir -p /scratch/users/$USER
fi
WORK_DIR=/scratch/users/$USER/$PBS_JOBID
mkdir -p "$WORK_DIR"

# Initialize output
echo "Job $PBS_JOBID running on host $HOSTNAME"

# Gaussian environment setup
GAUSSIAN_FOLDER=/home/apps/g16
source $GAUSSIAN_FOLDER/g16/bsd/g16.profile
export GAUSS_EXEDIR=$GAUSSIAN_FOLDER/g16
export GAUSS_SCRDIR=$WORK_DIR

# Go to the working directory where the input files are located
cd "$PBS_O_WORKDIR"
echo "Current working directory: $PBS_O_WORKDIR"

# Validate the input file list
if [ ! -f "$FILE_LIST" ]; then
  echo "Error: File list $FILE_LIST not found."
  exit 1
fi

LIST=$(cat "$FILE_LIST")

if [ -z "$LIST" ]; then
  echo "Error: No input files specified in $FILE_LIST."
  exit 1
fi

# Function to run Gaussian calculations
run_calculation() {
  local input_file=$1
  local log_file="${input_file%.com}.log"

  echo "Running: $input_file" | tee -a progress.log
  g16 "$input_file"

  if (( $? != 0 )); then
    echo "Failed: $input_file (g16 error)" | tee -a progress.log
    return 1
  fi

  echo "Transforming and compressing $input_file..." | tee -a progress.log
  formchk -3 "${input_file%.com}.chk" "${input_file%.com}.fchk"

  if (( $? == 0 )); then
    rm -f "${input_file%.com}.chk"
    gzip -f "${input_file%.com}.fchk"
    echo "Completed: $input_file" | tee -a progress.log
  else
    echo "Failed: $input_file (formchk error, .chk preserved)" | tee -a progress.log
    return 1
  fi
}

# Export the function for parallel execution
export -f run_calculation
export CORES_PER_CALC
export GAUSS_SCRDIR
export GAUSS_EXEDIR
export WORK_DIR

# Run calculations in parallel
echo "$LIST" | xargs -n 1 -P "$MAX_PARALLEL_CALCS" -I {} bash -c 'run_calculation "{}"'

# Clean up the scratch directory
echo "Cleaning up scratch directory: $WORK_DIR"
#rm -rf "$WORK_DIR"

# Job completed
echo "Job $PBS_JOBID completed on host $HOSTNAME."

