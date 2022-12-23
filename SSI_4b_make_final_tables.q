#!/bin/bash

#
#

#SBATCH --job-name=final_table
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100gb
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
module load python/3.9.4

Rscript SSI_4c_CSIRO_diversity.R

python SSI_4c_generate_final_table.py

duration=$(( SECONDS - start ))

echo "Time (s) for individual array script to run: "$duration

#reset the trap upon successful completion
trap : 0

echo >&2 '
************
*** DONE ***
************
'
