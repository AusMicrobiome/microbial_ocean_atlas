Scripts were executed on CSIRO petrichor x86-64 server using the following software:

`SUSE Linux Enterprise Server 15 SP3`

`slurm 20.11.9`

`python/3.9.4`

`R/4.1.3`

The below workflow was used to produce the data outputs

1) clone the git repo to a local directory

`git clone https://github.com/AusMicrobiome/microbial_ocean_atlas.git ./`

2) Run the input file creation script to create the input files for each amplicon

`./SSI_1_batch_input_file_creation.sh`

The script will call the following scripts:
        
`SSI_1b_input_file.q`
            
`SSI_1c_make_atlas_input.py`

2a) Check for failed slurms:
`grep "CANCELLED" *.out`

Re-run if script fails

3) Run the script to batch the ASV bootstrap script for each amplicon/parameter

`./SSI_2_batch_ASV_bootstraps.sh`

The script will call the following scripts:

`SSI_2b_bootstrap_ASV_array.q`

`SSI_2c_CSIRO_MOA_bootstraps.R`

3a) Check for failed slurms:

`grep "CANCELLED" *.out`

Failed arrays can be re run by directly calling the script and passing appropriate amplicon and array number to the script, as below:

`module load R/4.1.3`

`Rscript SSI_CSIRO_MOA_bootstraps.R "16S" "579"`

If no failed slurms, remove the slurm out files

`rm *.out`

4) Run the script to generate the samplewise statistics files for each amplicon/environmental parameter

`./SSI_3_batch_samplewise_indexes.sh`

script will call the following scripts:

`SSI_3b_combine_samplewise_indexes.q`

`SSI_3c_CSIRO_sample_indexes.R`

4a) Check for failed slurms:

`grep "CANCELLED" *.out`

Failed arrays can be re run by directly calling the script and passing appropriate amplicon and array number to the script, as below:

`module load R/4.1.3`

`Rscript SSI_CSIRO_sample_indexes.R "16S" "Temperature"`

If no failed slurms, remove the slurm out files

`rm *.out`

5) Run the script to generate the final table for IMOS visualisation. The script will also clean up the directory and save appropriate run files.

`./SSI_4_batch_create_final_tables.sh`

script will call the following scripts:
`SSI_4b_make_final_tables.q`

`SSI_4c_CSIRO_diversity.R`

`SSI_4c_generate_final_table.py`

5a) Check for failed slurms:
`grep "CANCELLED" *.out`

Re-run if script fails

6) Clean up directory

remove modified batch scripts

`rm *_mod.q`

remove slurm out files

`rm *.out`

