#!/bin/bash

#make a copy of the batch script - we will modify it to contain the number of parameters to analyse - it will be deleted at the end
cp SSI_1b_input_file.q SSI_1b_input_file_mod.q

#fill the user info
echo $USER
sed -i 's/XXXX/'${USER}'/g' SSI_1b_input_file_mod.q

sbatch SSI_1b_input_file_mod.q

