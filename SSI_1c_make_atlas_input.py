#!/usr/bin/env python

'''
make_ocean_atlas_input.py

A script generate a input file to calculate community indexes for Marine samples held in the AM
The script will:
1) Query the Australian Microbiome database for samples designated as 'Pelagic' and 'Coastal water' and retrieve specific physiochemical parameters  
2) Merge pysiochemical information with ASV, abundance information for each sample 
'''

import pandas as pd
import glob
import sqlite3
import os
from pandas.io import sql
from datetime import datetime

#####
# Start Database functions
#
#Function to connect to the current DB
def connect_DB():
    global db_conn
    global c
    global db_name
    db_name = glob.glob("/datasets/work/oa-amd/work/amd-work/SQlite/*.db")[0]
    print(db_name)
    #connect to the sqlite database 
    db_conn = sqlite3.connect(db_name)
    #create a cursor object to communicate with SQL
    c = db_conn.cursor()

#Function to close the DB cursor
def close_DB():
    c.close()
    #Close the database connections
    db_conn.close()

'''
Function to return desired fields based on a search, variable formats are as below
returnfields=['*'] to return all fields; returnfields=['sample_id','lat_lon',...] to return selected fields
searchfield should be defined as a string e.g.: searchfield = 'sample_id'
var_list is a list to search for: var_list=['102.100.100/7031','102.100.100/7032',....]
usage selectColumns_by_fieldValues(returnfields, searchfield, var_list)
'''
def selectColumns_by_fieldValues(returnfields, searchfield, var_list):
    global selectFields
    connect_DB()
    getFields = ','.join(returnfields)
    search = '\',\''.join(var_list)
    search = "\'" + search + "\'"
    select_cmd = "SELECT " + getFields + " FROM AM_metadata WHERE " + searchfield + " IN (" + search + ")" 
    #print(select_cmd)
    selectFields = pd.read_sql(select_cmd, db_conn)
    print("output df = selectFields")
    print(selectFields.shape)
    close_DB()
# End Database Functions

#function io write filenames of files used to a metadata file 
def write_files_used(filename):
    with open('AM_input_file_metadata.txt', mode='a') as f:
        f.write(f'{filename}\n')

#########
#Query the AMDB for the samples we need
#########
returnfields = ['sample_id','temp', 'salinity', 'nitrate_nitrite', 'phosphate', 'silicate','oxygen','sample_integrity_warnings']
searchfield = 'sample_type'
var_list = ['Pelagic','Coastal water']
selectColumns_by_fieldValues(returnfields, searchfield, var_list)

#########
# Format the df for use by the R script
#########
print("Removing samples with integrity warnings")
selectFields = selectFields[~selectFields['sample_integrity_warnings'].notnull()]
drops = ['sample_integrity_warnings']
selectFields.drop(drops, axis=1,inplace=True)
selectFields.rename({'sample_id': 'Sample_only'}, axis=1, inplace=True)
#Format the sample ID to short form by removing part of the handle
selectFields['Sample_only'] = selectFields['Sample_only'].str.replace('102.100.100/','',regex=True)
#Remove NCBI samples as they cause problems for the Rscript
selectFields = selectFields[~selectFields['Sample_only'].str.contains('SAMN', regex=True)]
#write the DB file used to the metadata file
db_file = os.path.basename(db_name)
write_files_used(db_file)
print(selectFields.shape)

All_Amplicons = ['16S','A16S','18Sv4']
for amplicon in All_Amplicons:
    #########
    # Get the ASV abundance table and read it
    #########
    amplicon_path = glob.glob('/datasets/work/oa-amd/work/amd-work/zotutabs/' + amplicon + '/short/*_all_20K_combined.txt')[0]
    abund_file = glob.glob(amplicon_path)[0]
    print("reading abundance file: " + abund_file)
    OTUs = pd.read_csv(abund_file, sep=('\t'), low_memory=False)
    print("shape of OTUs:")
    print(OTUs.shape)
    #save the input abundance used to the metadata file
    abund_f = os.path.basename(abund_file)
    write_files_used(abund_f)

    #########
    # Merge abundance and metadata
    #########
    merged = pd.merge(selectFields,OTUs,left_on='Sample_only', right_on='Sample_only', how='left')
    #drop any IDs that dont have ASV's - this can happen if the sample hasnt been sequenced for a particular amplicon
    merged.dropna(subset=['#OTU ID'],inplace=True)
    #Order the columns 
    merged = merged[['Sample_only', '#OTU ID', 'Abundance', 'Abundance_20K', 'temp', 'salinity', 'nitrate_nitrite', 'phosphate', 'silicate','oxygen']]
    print("shape of merged abundance and metadata:")
    print(merged.shape)

    #########
    # Save time stamped file
    #########
    dt = datetime.now().strftime("%Y%m%d%H%M")
    outfile_all = amplicon + "_" + dt + "_infile_multi_params_all.csv"
    print("saving merged abundance and metadata: " + outfile_all)
    merged.to_csv(outfile_all, index=False)

