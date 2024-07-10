#!/bin/bash

#
#

#SBATCH --job-name=R_samplewise
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=80gb
#SBATCH -a 1-YYYY
#SBATCH --mail-user=XXXX@csiro.au
#SBATCH --mail-type=END

#set up an abort trap to exit and advise if something goes horribly wrong
abort(){
    echo >&2 '
*****************
*** CANCELLED ***
*****************
'
    echo "An error occurred. Exiting..." >&2
    exit 1
}

trap 'abort' 0

set -e

#end trap

start=$SECONDS

module load R/4.1.3	

amplicon=$1

echo $amplicon

param=$(sed -n ${SLURM_ARRAY_TASK_ID}p  parameters.txt)

echo $param

Rscript SSI_3c_CSIRO_sample_indexes.R "$amplicon" "$param"

duration=$(( SECONDS - start ))

echo "Time (s) for individual array script to run: "$duration

#reset the trap upon successful completion
trap : 0

echo >&2 '
************
*** DONE ***
************
'
