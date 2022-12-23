#!/bin/bash

#make a copy of the batch script - we will modify it to contain the number of parameters to analyse - it will be deleted at the end
cp SSI_4b_make_final_tables.q SSI_4b_make_final_tables_mod.q

#fill the user info
echo $USER
sed -i 's/XXXX/'${USER}'/g' SSI_4b_make_final_tables_mod.q

sbatch SSI_4b_make_final_tables_mod.q

