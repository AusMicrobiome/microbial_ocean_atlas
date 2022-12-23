#!/bin/bash

#
#

#SBATCH --job-name=R_bootstrap
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=30gb
#SBATCH -a 1-1000
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

##Get the amplicon passed from the setup script
amplicon=$1

echo $amplicon

rep_num=${SLURM_ARRAY_TASK_ID}

echo $rep_num

Rscript SSI_2c_CSIRO_MOA_bootstraps.R "$amplicon" "$rep_num"

duration=$(( SECONDS - start ))
echo "Time (s) for individual array script to run: "$duration

#reset the trap upon successful completion
trap : 0

echo >&2 '
************
*** DONE ***
************
'
