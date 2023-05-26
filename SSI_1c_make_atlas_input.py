#!/usr/bin/env python

'''
SSI_1c_make_atlas_input.py

A script generates an input file to calculate community indexes for Marine samples held in the AM
The script will:
1) Query the Australian Microbiome database for samples designated as 'Pelagic' and 'Coastal water' and retrieve specific physiochemical parameters
2) Add the sample specific physiochemical metadata to ASV, abundance and abundance_20K information for each sample
3) Script will generate a file containing the abundance (subsampled to 20K reads) of taxa of interest.
    - Taxa of interest are defined in the file `NRS_taxon_list.csv` and follow Silva138 taxonomy formatted to include taxonomic level descriptions (e.g., d__<name>,p__<name>,c__<name>,o__<name>,f__<name>,g__<name>,s__<name>).
    - Classifications from SKlearn are used.
'''

import pandas as pd
import numpy as np
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
    
'''
Function to return desired fields based on partial match, variable formats are as below
returnfields=['*'] to return all fields; returnfields=['sample_id','lat_lon',...] to return selected fields
searchfield should be defined as a string: e.g., searchfield = 'sample_id'
wild_card is a string to search for : e.g., wild_card='7031%' or wild_card='%7031%' (% sign will indicate trailing or wrapping chars)
usage: 
selectColumns_by_partialValues(returnfields, searchfield, wild_card)
'''
def selectColumns_by_partialValues(returnfields, searchfield, wild_card):
    global stringMatch
    connect_DB()
    getFields = ','.join(returnfields)
    select_cmd = "SELECT " + getFields + " FROM AM_metadata WHERE " + searchfield + " LIKE \'" + wild_card + "\'" 
    print(select_cmd)
    stringMatch = pd.read_sql(select_cmd, db_conn)
    print("output df = stringMatch")
    print(stringMatch.shape)
    close_DB()

# End Database Functions

def write_files_used(filename):
    with open('files_analysed.txt', mode='a') as f:
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
selectFields.fillna(np.nan, inplace=True)
print("Removing samples with integrity warnings")
selectFields = selectFields[~selectFields['sample_integrity_warnings'].notnull()]
drops = ['sample_integrity_warnings']
selectFields.drop(drops, axis=1,inplace=True)
selectFields.rename({'sample_id': 'Sample_only'}, axis=1, inplace=True)
#Format the sample ID to short form by removing part of the handle
selectFields['Sample_only'] = selectFields['Sample_only'].str.replace('102.100.100/','',regex=True)
#Remove NCBI samples as they cause problems for the Rscript
selectFields = selectFields[~selectFields['Sample_only'].str.contains('SAMN', regex=True)]
#change the index to the sample only colunm - we will use the index to create the key for the dict
selectFields.set_index('Sample_only', inplace=True)
#grab the fields we are working with
df_cols = selectFields.columns.tolist()
#make a dictionary of each row of the metadata dataframe
meta_dict = selectFields.to_dict(orient='index')

All_Amplicons = ['16S','A16S','18Sv4']

#get all NRS samples to crete the taxonomy abundance table
returnfields = ['sample_id','depth', 'nrs_trip_code', 'nrs_sample_code','sample_integrity_warnings']
searchfield = 'imos_site_code'
wild_card = 'NRS%'
selectColumns_by_partialValues(returnfields, searchfield, wild_card)

print("Removing samples with integrity warnings")
stringMatch = stringMatch[~stringMatch['sample_integrity_warnings'].notnull()]
stringMatch['Sample_only'] = stringMatch['sample_id'].apply(lambda x: x.replace('102.100.100/','') if 'SAMN' not in x else x)
print("shape of NRS samples df")
print(stringMatch.shape)

#get a list of the NRS sample IDs 
NRS_IDs = set(stringMatch['Sample_only'].values.tolist())
# make an empty list to hold a list of each line matching a NRS ID 
NRS_list = []

dt = datetime.now().strftime("%Y%m%d%H%M")


for amplicon in All_Amplicons:
    print(amplicon)
    #########
    # Get the ASV abundance table and read it
    #########
    
    #Sample_only,#OTU ID,Abundance,Abundance_20K,temp,salinity,nitrate_nitrite,phosphate,silicate,oxygen
    amplicon_path = glob.glob('/datasets/work/oa-amd/work/amd-work/zotutabs/' + amplicon + '/short/*_all_20K_combined.txt')[0]
    abund_file = glob.glob(amplicon_path)[0]

    outfile_all = amplicon + "_" + dt + "_infile_multi_params_all.csv"
    
    print("saving merged abundance and metadata: " + outfile_all)

    #########
    # Open the ASV abundance table and output file to write the combined data data
    #########

    with open(abund_file, mode='r') as abund,\
    open(outfile_all, mode='w') as out:
        for line in abund:
            if '#OTU ID' in line:
                header = "Sample_only,#OTU ID,Abundance,Abundance_20K,"
                #add the metadata columns used
                meta_cols = ','.join(df_cols)
                #join the abundance and meta headers
                header = header + meta_cols
                out.write(f'{header}\n')
                NRS_header = []
                line = line.rstrip('\n')
                for ele in line.split('\t'):
                    NRS_header.append(ele)
                continue
            otu,ID,Abund,Abund_20K = line.split('\t')
            Abund_20K = Abund_20K.rstrip()
            if ID in meta_dict:
                tmp_string = ''
                for col in meta_dict[ID]:
                    tmp_string = tmp_string + str(meta_dict[ID].get(col)) + ","
                tmp_string = tmp_string.rstrip(',')
                tmp = f'{ID},{otu},{Abund},{Abund_20K},{tmp_string}\n'
                out.write(f'{ID},{otu},{Abund},{Abund_20K},{tmp_string}\n')
            if set(line.split('\t')[1:-2]) & NRS_IDs:
                tmp = []
                line = line.rstrip('\n')
                for ele in line.split('\t'):
                    tmp.append(ele)
                NRS_list.append(tmp)

#make a pandas dataframe and format the columns ready for merging with taxonomy
NRS_df = pd.DataFrame(NRS_list, columns=NRS_header).fillna(np.nan)
NRS_df.drop(['Abundance'],axis=1,inplace=True)
NRS_df['Abundance_20K'].replace('', np.nan, inplace=True)
print(NRS_df.shape)
NRS_df.dropna(subset=['Abundance_20K'], inplace=True)
print("Shape of requested NRS amplicon abndance table: " )
print(NRS_df.shape)
                
#process the NRS samples to get the taxa required
Taxonomies= []
with open('NRS_taxon_list.csv', mode='r') as f:
    for line in f:
        line = line.strip()
        #skip the header
        if '#' in line:
            continue
        #skip any blank lines
        if len(line.strip()) == 0:
            continue
        tmp = []
        for ele in line.split(','):
            tmp.append(ele.rstrip('\n'))
        Taxonomies.append(tmp)
print(Taxonomies)

get_amplicon_tax = set()
for taxon in Taxonomies:
    if taxon[0] == 'd__Bacteria':
        get_amplicon_tax.add('16S')
    if taxon[0] == 'd__Eukaryota':
        get_amplicon_tax.add('18Sv4')
    if taxon[0] == 'd__Archaea':
        get_amplicon_tax.add('A16S')
print(get_amplicon_tax)

# read in and combine the Taxonomy tables
combined_tax_table = pd.DataFrame()
for amplicon in get_amplicon_tax:
    print("Getting taxonomy for: " + amplicon)
    tax = pd.read_csv(glob.glob('/datasets/work/oa-amd/work/amd-work/zotutabs/' + amplicon + '/short/*.silva138.SKlearn.taxonomy')[0], sep='\t')
    if "Traits" in tax.columns:
        tax.drop(['Traits'], axis=1,inplace=True)
    tax.drop(['confidence'], axis=1,inplace=True)
    drops = ['d__Unassigned', 'd__Unclassified']
    tax = tax[~tax['kingdom'].isin(drops)]
    print("amplicon " + amplicon + " table shape")
    print(tax.shape)
    combined_tax_table = pd.concat([combined_tax_table, tax], ignore_index=True)
print("Shape of combined amplicon taxonomy tables")
print(combined_tax_table.shape)
    
#merge the dataframes
NRS_tax = pd.merge(NRS_df, combined_tax_table, left_on='#OTU ID', right_on='#OTU ID', how='left')
NRS_tax['Abundance_20K'] = NRS_tax['Abundance_20K'].astype('int')
NRS_tax.drop(['#OTU ID'],axis=1,inplace=True)
NRS_tax = NRS_tax.groupby(['Sample_only','kingdom','phylum','class','order','family','genus','species','amplicon'])['Abundance_20K'].sum().reset_index()

print("Shape of merged abundance taxonomy table")
print(NRS_tax.shape)

#search for the requested taxonomies
NRS_abund_tax = pd.DataFrame()
taxon_columns =  []
for taxon in Taxonomies:
    #make the string of all but the last taxon level to filter the final results
    taxon_string = (',').join(taxon[:-1])
    print(taxon_string)
    #get the number of levels of requested taxonomies and make it relative to the positions in the table table
    taxTable_pos = len(taxon)
    print(taxTable_pos)
    print(taxon[-1])
    tmp = NRS_tax.copy()
    #filter by the final taxonomic level provided
    tmp = tmp[tmp.iloc[:, taxTable_pos].str.contains(taxon[-1])]
    #make a string of the taxonomic levels above the final level
    tmp['tax_string'] = tmp.iloc[:, 1:taxTable_pos].apply(lambda row: ','.join(row.values.astype(str)), axis=1)
    #filter by the remainder of the tax string as a final check
    tmp = tmp[tmp['tax_string'].isin(list(taxon_string.split(" ")))]
    #if the taxonomy string terminates before species we will replace the lower taxonomies with the same text and regroup 
    if taxTable_pos < 7:
        tmp.iloc[:,taxTable_pos+1:8] = taxon[-1]
        tmp = tmp.groupby(['Sample_only','kingdom','phylum','class','order','family','genus','species','amplicon'])['Abundance_20K'].sum().reset_index()
    #add a column descriptor for describing the taxon requested ( lowest level taxonomy (last level requested))
    tmp['taxon'] = taxon[-1] + "_abundance_20K"
    taxon_columns.append(taxon[-1] + "_abundance_20K")
    NRS_abund_tax = pd.concat([NRS_abund_tax, tmp], ignore_index=True)
    print(NRS_abund_tax.shape)
print('final taxonomy table shape')
print(NRS_abund_tax.shape)

#merge the abundance and taxonomy tables
NRS_abund_tax = NRS_abund_tax[['Sample_only','Abundance_20K','amplicon','taxon']]
NRSabund_meta = pd.merge(NRS_abund_tax,stringMatch, left_on='Sample_only', right_on='Sample_only', how='left')
print("Merged abundance and taxonomy shape")
print(NRSabund_meta.shape)

#pivot the table to get the desired format and format the numeric format for export
NRS_pivoted = pd.pivot_table(NRSabund_meta,index=['sample_id','depth','nrs_trip_code','nrs_sample_code','amplicon'], columns= ['taxon'], values=['Abundance_20K'])
print(NRS_pivoted.shape)
NRS_pivoted.columns = NRS_pivoted.columns.droplevel(0)
for i in NRS_pivoted.index:
    for col in taxon_columns:
        try:
            NRS_pivoted.loc[i,col] = str(int(NRS_pivoted.loc[i,col]))#.astype('Int64')
        except:
            pass

#sort the df on the NRS sample code to make it easier for a human to read
NRS_pivoted.sort_values(by=['nrs_sample_code'],ascending=True,inplace=True)

#write the outut to csv
final_NRSabund_file = os.path.join(os.path.abspath(os.getcwd()),'data/NRS_taxon_abundance.csv')
print("writing NRS taxonomy abundances to output file: " + final_NRSabund_file)
NRS_pivoted.to_csv(final_NRSabund_file)
NRS_pivoted.head()
