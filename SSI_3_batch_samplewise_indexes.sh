#!/bin/bash

#make a copy of the batch script - we will modify it to contain the number of parameters to analyse - it will be deleted at the end
cp SSI_3b_combine_samplewise_indexes.q SSI_3b_combine_samplewise_indexes_mod.q
#make the parameters file based on the 16S array physchem parameters 
ls -1 16S*subset1.csv | cut -f4 -d"_" > parameters.txt;
count=$(wc -l < parameters.txt);
echo $count;
sed -i 's/YYYY/'${count}'/g' SSI_3b_combine_samplewise_indexes_mod.q

#fill the user info
echo $USER
sed -i 's/XXXX/'${USER}'/g' SSI_3b_combine_samplewise_indexes_mod.q

#batch the scripts for each amplicon
for amplicon in 16S A16S 18Sv4
do
echo $amplicon
sbatch SSI_3b_combine_samplewise_indexes_mod.q $amplicon
done
