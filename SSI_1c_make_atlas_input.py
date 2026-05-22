import pandas as pd
import numpy as np
import glob
import sqlite3
import os
from pandas.io import sql
from datetime import datetime
import sys
from dotenv import load_dotenv

#####
# Start Database functions
#
#Function to connect to the current DB
def connect_DB():
    global db_conn
    global c
    global db_name
    db_name = glob.glob(DB)[0]
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
    connect_DB()
    getFields = ','.join(returnfields)
    search = '\',\''.join(var_list)
    search = "\'" + search + "\'"
    select_cmd = "SELECT " + getFields + " FROM AM_metadata WHERE " + searchfield + " IN (" + search + ")"
    #print(select_cmd)
    selectFields = pd.read_sql(select_cmd, db_conn)
    print(f'Shape of meatadata df: {selectFields.shape}')
    close_DB()
    return selectFields
# End Database Functions

def write_files_used(filename):
    with open('AM_input_file_metadata.txt', mode='a') as f:
        f.write(f'{filename}\n')
'''
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
'''
#A function to find dulicates in a list
def find_dups(inlist):
    seen = set()
    dups = []
    for x in inlist:
        if x in seen:
            dups.append(x)
        else:
            seen.add(x)
    return list(set(dups))

#load the environmental variables
env_path = os.getenv("MICRO_OCEAN_ATLAS_ENV")
load_dotenv(env_path)

if not env_path:
    raise EnvironmentError(
        "Environment variable MICRO_OCEAN_ATLAS_ENV is not set. "
        "Please add `export path_to_microbial_ocean_atlas.env` to your ~/.bashrc"
    )

if not os.path.isfile(env_path):
    raise FileNotFoundError(
        f"The .env file was not found at: {env_path}. "
        "Please check the path set in ENV_VARS."
    )

ZOTUS =  os.getenv("ZOTUS")
print(ZOTUS)

DB = os.getenv("DB_PATH")
print(DB)

#########
#Query the AMDB for the samples we need
#########
returnfields = ['sample_id',
                'temp', 
                'salinity', 
                'nitrate_nitrite', 
                'phosphate', 
                'silicate',
                'oxygen',
                'sample_integrity_warnings', 
                'nrs_trip_code',
                'synonyms',
                'imos_site_code',
                'sample_site_location_description']
searchfield = 'sample_type'
var_list = ['Pelagic','Coastal water']
selectFields = selectColumns_by_fieldValues(returnfields, searchfield, var_list)
AMDB_file = os.path.basename(db_name)
write_files_used(AMDB_file)

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

tax_req_df = pd.read_csv('taxa_request.list', comment="#", names=['taxonomy'])
tax_req_df['lowest'] = tax_req_df['taxonomy'].str.split(';').str[-1]
print(f'Number of requested Taxa: {len(tax_req_df.index)}')

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
print(f'Requested Traits: {traits_wanted}')

#Define the amplicon we will use to grab the taxa from
domain_amp_d = {'16S': 'd__Bacteria', '18Sv4': 'd__Eukaryota', 'A16S': 'd__Archaea'}

All_Amplicons = ['16S', '18Sv4', 'A16S']

#make some empty dfs to hold combined outputs for  all amplicons
ASV_df = pd.DataFrame()
trait_counts = pd.DataFrame()
tax_df = pd.DataFrame()
trait_df = pd.DataFrame()
for amplicon in All_Amplicons:
    if '16S' in amplicon:
        tax_file = glob.glob(os.path.join(ZOTUS, f'{amplicon}/short/*.silva138.SKlearn.taxonomy'))[0]
        #dict of taxonomic level prefixes
        level_prefix_d = {'domain':'d__',
                          'phylum': 'p__',
                          'class':'c__',
                          'order':'o__',
                          'family':'f__',
                          'genus':'g__',
                          'species':'s__'}
        
    if amplicon == '18Sv4':
        tax_file = glob.glob(os.path.join(ZOTUS, f'{amplicon}/short/*.PR2v500.wang.taxonomy'))[0]
        #dict of taxonomic level prefixes
        level_prefix_d = {'domain':'d__',
                          'supergroup':'sg__',
                          'division':'div__',
                          'subdivision':'sdiv__',
                          'class':'c__',
                          'order':'o__',
                          'family':'f__',
                          'genus':'g__',
                          'species':'s__'}
    taxonomy_file = os.path.basename(tax_file)
    write_files_used(taxonomy_file)
    tax = pd.read_csv(tax_file, sep='\t', dtype=str)
    print(f'Shape of taxonomy file with all OTUs: {tax.shape}')
    print(f'Filtering taxonomy for target domain: {domain_amp_d[amplicon]}')
    tax = tax[tax['domain'].str.contains(domain_amp_d[amplicon])].reset_index(drop=True)
    print(f'Shape of taxonomy file filtered to target domain: {tax.shape}')
    targetOTUs = tax['#OTU ID'].values.tolist()
    if amplicon != '18Sv4': # 18Sv4 as it doesnt have that traits analysed
        print("Getting taxonomy for: " + amplicon)
        tmp = tax[tax['traits'].fillna('').str.contains('|'.join(list(traits_wanted)))].reset_index(drop=True)
        trait_df = pd.concat([trait_df, tmp], ignore_index=True)
        print(f'Shape of dataframe holding requested traits: {trait_df.shape}')
    prefix_level_d = {v: k for k, v in level_prefix_d.items()}
    tax_req_target = tax_req_df[tax_req_df['taxonomy'].str.contains(domain_amp_d[amplicon])].reset_index(drop=True)
    if len(tax_req_target.index) > 0:
        for i, (k, v) in enumerate(prefix_level_d.items()):
            level_mask = tax_req_target['lowest'].str.contains(k)
            if level_mask.sum() > 0:
                taxa_names = tax_req_target.loc[level_mask, 'lowest'].values.tolist()
                tmp = tax[tax[v].isin(taxa_names)].reset_index(drop=True)
                if len(tmp.index) > 0:
                    print(f'Shape of dataframe matching requested taxa at {v}: {tmp.shape}')
                    #if the genus is not in the species name we will overwrite it to the genus_sp.
                    tmp['genus_only'] = tmp['genus'].replace('g__','', regex=True).str.split('_').str[0]
                    species_mask = tmp['species'].str.contains('s__')
                    #tax_df.loc[species_mask, 'species_only'] = 
                    tmp.loc[species_mask, 'species_only'] = tmp.loc[species_mask, 'species'].replace('s__','', regex=True).replace('_', ' ', regex=True)
                    sp_mask = (tmp['species_only'].notna()) & (tmp['genus_only'] != tmp['species_only'].str.split(' ').str[0])
                    if sp_mask.sum() > 0:
                        tmp.loc[sp_mask, 'species'] = 's__' + tmp.loc[sp_mask, 'genus'].replace('g__', '', regex=True) + '_sp.'
                    drops = ['genus_only','species_only']
                    tmp = tmp.drop(drops, axis=1)
                    tax_df = pd.concat([tax_df, tmp], ignore_index=True)
                    tax_df = tax_df.drop_duplicates().reset_index(drop=True)
                    #test the df to make sure there are no duplicate OTU IDs 
                    dups = find_dups(tax_df['#OTU ID'].values.tolist())
                    if len(dups) > 0:
                        print(f'{len(dups)} Duplicate OTUs in taxonomy file: {dups}')
                        sys.exit(1)
                del tmp
                tax_cols_order = ['#OTU ID',
                                  'domain',
                                  'supergroup',
                                  'division',
                                  'subdivision',
                                  'phylum',
                                  'class',
                                  'order',
                                  'family',
                                  'genus',
                                  'species',
                                  'amplicon',
                                  'traits']
                reorder = []
                for col in tax_cols_order:
                    if col in tax_df:
                        reorder.append(col)
                tax_df = tax_df[reorder] 
    abund_file = glob.glob(os.path.join(ZOTUS, f'{amplicon}/short/*_all_20K_combined.txt'))[0]
    print(abund_file)
    abundance_file = os.path.basename(abund_file)
    write_files_used(abundance_file)
    abund = pd.read_csv(abund_file, sep='\t', dtype=str)
    abund['amplicon_name'] = amplicon
    abund = pd.merge(abund, selectFields, left_on='Sample_only', right_on='Sample_only', how='inner')
    print(f'"Shape of combined ASV metadata table for {amplicon}": {abund.shape}')
    abund['Abundance'] = abund['Abundance'].astype(int)
    abund['Abundance_20K'] = abund['Abundance_20K'].astype(float)
    #filter the df to target 
    print(f'filtering OTU table to target domain: {domain_amp_d[amplicon]}')
    print(f'shape of OTU table before screening for target OTUs: {abund.shape}')
    abund = abund[abund['#OTU ID'].isin(targetOTUs)].reset_index(drop=True)
    print(f'shape of OTU table after screening for target OTUs: {abund.shape}')
    ASV_df = pd.concat([ASV_df, abund], ignore_index=True)
    outfile_all = amplicon + "_infile_multi_params_all.csv"
    abund_drops = ['nrs_trip_code',
                   'synonyms',
                   'imos_site_code',
                   'sample_site_location_description']
    abund = abund.drop(abund_drops, axis=1)
    abund.to_csv(outfile_all, index=False)
    del abund

#get the total number and counts of zotus for each sample/amplicon
print("Getting ASV abundance info")
for amplicon in All_Amplicons:
    amplicon_abund_col =  amplicon + "_total_abundance"
    amplicon_count_col =  amplicon + "_total_count"
    amplicon_abund20K_col =  amplicon + "_abundance_20K"
    amplicon_count20K_col =  amplicon + "_count_20K"
    tmp = ASV_df[ASV_df['amplicon_name'] == amplicon].groupby(['Sample_only','amplicon_name'])['Abundance'].agg(['sum', 'count']).reset_index().rename(columns={'sum':amplicon_abund_col,'count':amplicon_count_col})
    ASV_df = pd.merge(ASV_df, tmp, left_on=['Sample_only','amplicon_name'], right_on=['Sample_only','amplicon_name'], how='left')
    tmp = ASV_df[ASV_df['amplicon_name'] == amplicon].groupby(['Sample_only','amplicon_name'])['Abundance_20K'].agg(['sum', 'count']).reset_index().rename(columns={'sum':amplicon_abund20K_col,'count':amplicon_count20K_col})
    ASV_df = pd.merge(ASV_df, tmp, left_on=['Sample_only','amplicon_name'], right_on=['Sample_only','amplicon_name'], how='left')

print("Calculating trait counts")
trait_counts = pd.merge(ASV_df, trait_df, left_on='#OTU ID', right_on='#OTU ID', how = 'inner')
print(f'Shape of traits counts df: {trait_counts.shape}')
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
    tmp = trait_counts[trait_counts['traits'].str.contains(trait)].groupby('Sample_only')['Abundance'].agg(['sum', 'count']).reset_index().rename(columns={'sum':trait_abund_col,'count':trait_count_col})
    trait_counts = pd.merge(trait_counts, tmp, left_on='Sample_only', right_on='Sample_only', how='left')
    tmp = trait_counts[trait_counts['traits'].str.contains(trait)].groupby('Sample_only')['Abundance_20K'].agg(['sum', 'count']).reset_index().rename(columns={'sum':trait_abund20K_col,'count':trait_count20K_col})
    trait_counts = pd.merge(trait_counts, tmp, left_on='Sample_only', right_on='Sample_only', how='left')
trait_counts = trait_counts[['#OTU ID'] + trait_count_order]
trait_counts = trait_counts.drop_duplicates(keep='first')
print(f'Shape of dataframe holding trait counts: {trait_counts.shape}')

print('Proceesing requested taxonomies')
tax_df = tax_df.drop_duplicates()
print(f'Shape of dataframe holding all requested taxonomies: {tax_df.shape}')
print("Merging taxonomy and trait info")
# outer merge the taxonomy and traits to generate a table for all sample ID that had awanted trait and/or taxonomy  
tax_traits_df = pd.merge(tax_df,trait_counts, left_on=['#OTU ID'], right_on=['#OTU ID'], how='outer')
print(f'Shape of merged taxonomy and traits df: {tax_traits_df.shape}')
tax_traits_df = tax_traits_df[trait_count_order]
tax_traits_df = tax_traits_df.drop_duplicates(keep='first')
print(f'Shape of dataframe final taxonomy trait counts: {tax_traits_df.shape}')

print("Merging taxonomy and ASV abundances")
#merge the taxonomies and abundance tables
tax_df = pd.merge(tax_df,ASV_df,left_on=['#OTU ID'],right_on=['#OTU ID'], how='inner')

print("Merging Abundance and taxonomies with trait counts")
tax_df = pd.merge(tax_df,tax_traits_df, left_on=ASV_sample_abund_cols, right_on=ASV_sample_abund_cols,how='outer')

print("Formatting final table")
#drop any columns that are all null
tax_df = tax_df.dropna(axis=1, how='all')
#add in the full AM sampleID
tax_df_cols = tax_df.columns.tolist()
tax_df['sample_id'] = tax_df.apply(lambda x: '102.100.100/' + (str(x['Sample_only'])), axis=1)
tax_df = tax_df[['sample_id'] + tax_df_cols ]
print(f'Shape of final requested taxonomy and traits Table: {tax_df.shape}')
for col in tax_df.columns:
    if "20K" in col:
        tax_df[col] = tax_df[col].apply(lambda x: np.nan if x == 0 else x)
        tax_df[col] = pd.to_numeric(tax_df[col]).astype('Int32')

print("Writing final table")
tax_df.to_csv('data/taxa_traits_ASV_request.csv', index=False)
print("DONE!")