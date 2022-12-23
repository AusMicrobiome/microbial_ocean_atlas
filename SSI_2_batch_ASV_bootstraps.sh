#!/bin/bash

#make a copy of the batch script - we will modify it to contain the number of parameters to analyse - it will be deleted at the end
cp SSI_2b_bootstrap_ASV_array.q SSI_2b_bootstrap_ASV_array_mod.q

#fill the user info
echo $USER
sed -i 's/XXXX/'${USER}'/g' SSI_2b_bootstrap_ASV_array_mod.q

#batch the scripts for each amplicon
for amplicon in 16S A16S 18Sv4
do
echo $amplicon
sbatch SSI_2b_bootstrap_ASV_array_mod.q $amplicon
done
