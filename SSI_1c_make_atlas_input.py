#!/usr/bin/env python

'''
make_ocean_atlas_input.py

A script generate a input file to calculate community indexes for Marine samples held in the AM
The script will:
1) Query the Australian Microbiome database for samples designated as 'Pelagic' and 'Coastal water' and retrieve specific physiochemical parameters  
2) Merge pysiochemical information with ASV, abundance information for each sample
3) Generate an ASV abundance table of requested taxa and traits 
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

def split_lines(line):
    global tmp
    tmp = []
    if trait_col == True:
        for ele in line.split('\t'):
            tmp.append(ele)
    else:
        for ele in line.split('\t'):
            tmp.append(ele)
        tmp.append('')
    return tmp

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


Taxonomies = []
lowest_tax = set()
with open('taxa_request.list', mode='r') as f:
    for line in f:
        #skip the header
        if '#' in line:
            continue
        #skip any blank lines
        if len(line.strip()) == 0:
            continue
        tmp = []
        for ele in line.split(';'):
            if '\n' in ele:
                lowest_tax.add(ele.strip())
            tmp.append(ele.strip())
        Taxonomies.append(tmp)
print(Taxonomies)
print(lowest_tax)

#Define the amplicon we will use to grab the taxa from
taxon_dict = {}
for taxon in Taxonomies:
    if taxon[0] == 'd__Bacteria':
        taxon_dict.update({taxon[0]:'16S'})
    if taxon[0] == 'd__Eukaryota':
        taxon_dict.update({taxon[0]:'18Sv4'})
    if taxon[0] == 'd__Archaea':
        taxon_dict.update({taxon[0]:'A16S'})
for domain, amplicon in taxon_dict.items():
    print(domain)
    print(amplicon)


All_Amplicons = ['16S','A16S','18Sv4']
ASV_list = []
for amplicon in All_Amplicons:
    abund_file = glob.glob('/datasets/work/oa-amd/work/amd-work/zotutabs/' + amplicon + '/short/*_all_20K_combined.txt')[0]
    outfile_all = amplicon + "_infile_multi_params_all.csv"
    
    print("saving merged abundance and metadata: " + outfile_all)

    #########
    # Open the ASV abundance table and output files for each aplicon to alculate indicies we wil also generate a lits to use to look at taxa and traits combined data
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
                ASV_header = []
                line = line.rstrip('\n')
                for ele in line.split('\t'):
                    ASV_header.append(ele)
                continue
            otu,ID,Abund,Abund_20K = line.split('\t')
            Abund_20K = Abund_20K.rstrip()
            if ID in meta_dict:
                tmp_string = ''
                for col in meta_dict[ID]:
                    tmp_string = tmp_string + str(meta_dict[ID].get(col)) + ","
                tmp_string = tmp_string.rstrip(',')
                #tmp = f'{ID},{otu},{Abund},{Abund_20K},{tmp_string}\n'
                out.write(f'{ID},{otu},{Abund},{Abund_20K},{tmp_string}\n')
                tmp = []
                line = line.rstrip('\n')
                for ele in line.split('\t'):
                    tmp.append(ele)
                tmp.append(amplicon)
                ASV_list.append(tmp)
print(len(ASV_list))

ASV_columns = ASV_header + ['amplicon_name']
ASV_df = pd.DataFrame(ASV_list, columns=ASV_columns)
#convert the abundance columns to mumbers
ASV_df['Abundance'] = ASV_df['Abundance'].astype('int')
#Abudance 20K will have trouble converting to int, we will leave it as a float for the moment and change to 'Int64' later
ASV_df['Abundance_20K'] = pd.to_numeric(ASV_df['Abundance_20K'])
print("Shape of combined ASV table")
print(ASV_df.shape)

#get the total number and counts of zotus for each sample/amplicon
print("Getting ASV abundance info")
for amplicon in All_Amplicons:
    amplicon_abund_col =  amplicon + "_total_abundance"
    amplicon_count_col =  amplicon + "_total_count"
    amplicon_abund20K_col =  amplicon + "_abundance_20K"
    amplicon_count20K_col =  amplicon + "_count_20K"
    #amplicon_count_order.extend([amplicon_abund_col,amplicon_count_col,amplicon_abund20K_col,amplicon_count20K_col])
    tmp = ASV_df[ASV_df['amplicon_name'] == amplicon].groupby(['Sample_only','amplicon_name'])['Abundance'].agg(['sum', 'count']).reset_index().rename(columns={'sum':amplicon_abund_col,'count':amplicon_count_col})
    ASV_df = pd.merge(ASV_df, tmp, left_on=['Sample_only','amplicon_name'], right_on=['Sample_only','amplicon_name'], how='left')
    tmp = ASV_df[ASV_df['amplicon_name'] == amplicon].groupby(['Sample_only','amplicon_name'])['Abundance_20K'].agg(['sum', 'count']).reset_index().rename(columns={'sum':amplicon_abund20K_col,'count':amplicon_count20K_col})
    ASV_df = pd.merge(ASV_df, tmp, left_on=['Sample_only','amplicon_name'], right_on=['Sample_only','amplicon_name'], how='left')
print(ASV_df.shape)

print("Getting requested traits and Taxa")
traits_wanted = set()
with open('trait_request.list', mode='r') as f:
    for line in f:
        if "#" in line:
            continue
        #skip any blank lines
        if len(line.strip()) == 0:
            continue
        traits_wanted.add(line.strip())
print("Requested Traits")
print(traits_wanted)

tax_list = []
trait_list = []
tax_columns = []
for amplicon in taxon_dict.values():
    print("Getting taxonomy for: " + amplicon)
    tax_file = glob.glob('/datasets/work/oa-amd/work/amd-work/zotutabs/' + amplicon + '/short/*.silva138.SKlearn.taxonomy')[0]
    print(tax_file)
    with open(tax_file, mode='r') as f:
        count = 0
        trait_col = False
        
        for line in f:
            line = line.rstrip('\n')
            if '#' in line:
                if 'traits' not in line:
                    trait_col = False
                else:
                    trait_col = True
                if tax_columns:
                    continue
                else:
                    for ele in line.split('\t'):
                        if '#' in ele:
                            tax_columns.append(ele)
                        else:
                            tax_columns.append(ele.title())
                    if 'Traits' not in tax_columns:
                        tax_columns.append('Traits')
                print(tax_columns)
            else:
                tmp = []
                for tax in lowest_tax:
                    if tax in line:
                        #use the function to generate a list from the line
                        split_lines(line)
                        tax_list.append(tmp)
                if amplicon == '18Sv4':
                    continue
                else:
                    tmp = []
                    for trait in traits_wanted:
                        if trait in line:
                            #use the function to generate a list from the line
                            split_lines(line)
                            trait_list.append(tmp)
print("Number of lines holding requested taxonomies")
print(len(tax_list))
print("Number of lines holding requested traits")
print(len(trait_list))

print("Calculating trait counts")
#convert the traits retrieved into a dataframe
trait_df = pd.DataFrame(trait_list, columns=tax_columns)
print("Shape of dataframe holding requested traits")
print(trait_df.shape)

#Merge the ASV abundance table withthe trait information requested so we can get the abundance of the identified traits
trait_count_order = []
trait_counts = pd.merge(ASV_df, trait_df, left_on='#OTU ID', right_on='#OTU ID', how = 'inner')

ASV_sample_abund_cols = ['Sample_only',
 'amplicon_name', 
 '16S_total_abundance',
 '16S_total_count',
 '16S_abundance_20K',
 '16S_count_20K',
 'A16S_total_abundance',
 'A16S_total_count',
 'A16S_abundance_20K',
 'A16S_count_20K',
 '18Sv4_total_abundance',
 '18Sv4_total_count',
 '18Sv4_abundance_20K',
 '18Sv4_count_20K']

#we will use the ASV_col list later so we will make a copy to extend
trait_count_order = ASV_sample_abund_cols.copy()
for trait in traits_wanted:
    trait_abund_col =  trait + "_abundance"
    trait_count_col =  trait + "_count"
    trait_abund20K_col =  trait + "_abundance_20K"
    trait_count20K_col =  trait + "_count_20K"
    trait_count_order.extend([trait_abund_col,trait_count_col,trait_abund20K_col,trait_count20K_col])
    tmp = trait_counts[trait_counts['Traits'].str.contains(trait)].groupby('Sample_only')['Abundance'].agg(['sum', 'count']).reset_index().rename(columns={'sum':trait_abund_col,'count':trait_count_col})
    trait_counts = pd.merge(trait_counts, tmp, left_on='Sample_only', right_on='Sample_only', how='left')
    tmp = trait_counts[trait_counts['Traits'].str.contains(trait)].groupby('Sample_only')['Abundance_20K'].agg(['sum', 'count']).reset_index().rename(columns={'sum':trait_abund20K_col,'count':trait_count20K_col})
    trait_counts = pd.merge(trait_counts, tmp, left_on='Sample_only', right_on='Sample_only', how='left')
trait_counts = trait_counts[['#OTU ID'] + trait_count_order]
trait_counts = trait_counts.drop_duplicates(keep='first')
print("Shape of dataframe holding trait counts")
print(trait_counts.shape)

#convert the taxonomies retrieved into a dataframe
tax_df = pd.DataFrame(tax_list, columns=tax_columns)
print("Shape of dataframe holding all requested taxonomies")
print(tax_df.shape)
#depending on the taxonomy levels searched (e.g. searching by genus and searching by species within that genus) we may have duplicate entries for the same ASV/taxonomy so we will drop them
tax_df = tax_df.drop_duplicates(keep='first')
print("Shape of dataframe holding unique requested taxonomies")
print(tax_df.shape)

print("Merging taxonomy and trait info")
# outer merge the taxonomy and traits to generate a table for all sample ID that had awanted trait and/or taxonomy  
tax_traits_df = pd.merge(tax_df,trait_counts, left_on=['#OTU ID'], right_on=['#OTU ID'], how='outer')
print(tax_traits_df.shape)
tax_traits_df = tax_traits_df[trait_count_order]
tax_traits_df = tax_traits_df.drop_duplicates(keep='first')
print("Shape of dataframe merged taxonomy trait counts")
print(tax_traits_df.shape)

print("Merging taxonomy and ASV abundances")
#merge the taxonomies and abundance tables
tax_df = pd.merge(tax_df,ASV_df,left_on=['#OTU ID'],right_on=['#OTU ID'], how='inner')

print("Merging Abundance and taxonomies with trait counts")
tax_df = pd.merge(tax_df,tax_traits_df, left_on=ASV_sample_abund_cols,right_on=ASV_sample_abund_cols,how='outer')
print(tax_df.shape)

print("Formatting final table")
#drop any columns that are all null
tax_df = tax_df.dropna(axis=1, how='all')
#add in the full AM sampleID
tax_df_cols = tax_df.columns.tolist()
tax_df['sample_id'] = '102.100.100/' + tax_df['Sample_only']
tax_df = tax_df[['sample_id'] + tax_df_cols ]
print(tax_df.shape)

for col in tax_df.columns:
    if "20K" in col:
        print(col)
        tax_df[col] = tax_df[col].apply(lambda x: np.nan if x == 0 else x)
        tax_df[col] = pd.to_numeric(tax_df[col]).astype('Int32')
print("Writing final table")
tax_df.to_csv('data/taxa_traits_ASV_request.csv', index=False)
print("DONE!")
