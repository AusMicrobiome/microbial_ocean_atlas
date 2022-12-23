#!/usr/bin/env python

'''
SSI_generate_final_table.py
'''

import pandas as pd
import numpy as np
import tarfile
import glob
import sqlite3
from pandas.io import sql
import os
import shutil

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
usage: selectColumns_by_fieldValues(returnfields, searchfield, var_list)
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
# Delete the subset bootstrap files
#########
for f in glob.glob("*_subset*.csv"):
    os.remove(f)

#########
# Set up: define parameters and make empty dataframe to hold the final table
#########
param_file = glob.glob(os.path.join(os.getcwd(),'parameters.txt'))[0]
parameters = pd.read_csv(param_file, names=['param'])
parameters_list = parameters['param'].values.tolist()
#Make an empty dataframe to hold the combined cti tables containing each amplicon
merge_cols = ['Sample_only'] + parameters_list
all_cti = pd.DataFrame(columns=merge_cols)
amplicons = ['16S', 'A16S', '18Sv4']
for amplicon in amplicons:
    #########
    # Define file paths
    #########
    tax_f = glob.glob('/datasets/work/oa-amd/work/amd-work/zotutabs/' + amplicon + '/short/*.silva132.SKlearn.taxonomy')[0]
    abund_f = glob.glob('/datasets/work/oa-amd/work/amd-work/zotutabs/' + amplicon + '/short/*_all_20K_combined.txt')[0]
    #write the taxonomy file used to the metadata file
    taxonomy_file = os.path.basename(tax_f)
    write_files_used(taxonomy_file)
    
    #########
    # Merge all the samplewise indexes
    #########
    indexes = glob.glob(os.path.join(os.getcwd(),amplicon + '*I.csv'))
    #create empty dataframe to merge the results into
    cti_merge = pd.DataFrame(columns=['Sample_only'])
    for f in indexes:
        df = pd.read_csv(f)
        df.drop('Unnamed: 0',axis=1, inplace=True)
        print(df.shape)
        cti_merge =  pd.merge(cti_merge, df, left_on='Sample_only',right_on='Sample_only',how='outer')

    # Get the sample IDs we are dealing with
    IDs = []
    for ID in cti_merge['Sample_only']:
        val = pd.isnull(ID)
        if val == True:
            continue
        else:
            IDs.append(str(ID))

    print("Reading abundance table")
    abund = pd.read_csv(abund_f, sep ='\t', low_memory=False)
    print(abund.shape)
    print("Selecting sampleIDs")
    abund = abund[abund['Sample_only'].isin(IDs)].astype(str)
    print(abund.shape)
    if amplicon == '18Sv4':
        #Currently no trait data for 18Sv4
        cti_merge['Sample_only'] = "102.100.100/" +  cti_merge['Sample_only'].astype(str)
    else:
        print("Reading taxonomy table")
        tax = pd.read_csv(tax_f,sep='\t')
        print(tax.shape)
        print("Selecting traits")
        tax['traits'].fillna('nan',inplace=True)
        tax = tax[tax['traits'].str.contains('fish_parasites|nitrogen_fixation')].reset_index(drop=True)
        print(tax.shape)
        tax = tax.reset_index(drop=True).astype(str)
        traits_merged = pd.merge(tax,abund,left_on='#OTU ID',right_on='#OTU ID', how='left')
        print(traits_merged.shape)
        traits_merged = traits_merged[traits_merged['Sample_only'].notna()].reset_index(drop=True)
        print(traits_merged.shape)
        traits_list = ['fish_parasites','nitrogen_fixation']
        trait_abund = []
        for trait in traits_list:
            col_name = amplicon + "_" + trait
            trait_abund.append(col_name)
            traits_merged[col_name] =  np.nan
            for i in traits_merged.index:
                val = traits_merged.loc[i,'traits']
                if trait in val:
                    traits_merged.loc[i,col_name] = traits_merged.loc[i,'Abundance_20K']
        traits_merged = traits_merged[['Sample_only'] + trait_abund]
        traits_merged = traits_merged.replace('nan', np.nan)
        traits_merged[trait_abund] = traits_merged[trait_abund].astype(float)
        traits_merged = traits_merged.groupby(['Sample_only'])[trait_abund].sum().reset_index()
        print("traits merged")
        print(traits_merged.shape)
        cti_merge['Sample_only'] = "102.100.100/" + cti_merge['Sample_only'].astype(str)
        traits_merged['Sample_only'] = "102.100.100/" + traits_merged['Sample_only'].astype(str)
        cti_merge = pd.merge(cti_merge,traits_merged, left_on='Sample_only', right_on='Sample_only', how='left')
        print("cti_merge and traits merged")
        print(cti_merge.shape)
    print("All_CTI shape")
    all_cti = pd.merge(all_cti,cti_merge, left_on=merge_cols, right_on=merge_cols, how='outer')
    print(all_cti.shape)
#sum the traits for each amplicon into one total column
for trait in traits_list:
    trait_columns = [s for s in all_cti.columns if trait in s]
    all_cti[trait] = all_cti[trait_columns].sum(axis=1)

#########
# Add the diversity outputs
#########
diversity = pd.read_csv('MOA_Diversity_stats.csv')
diversity['Sample_only'] = "102.100.100/" + diversity['Sample_only'].astype(str)
diversity.drop('Unnamed: 0',axis=1, inplace=True)
all_cti = pd.merge(all_cti,diversity, left_on='Sample_only', right_on='Sample_only', how='left')
print(diversity.shape)

#########
# Get the sample IDs in long and short format from the merged df that we are dealing with
#########
IDs = []
var_list = []
for ID in all_cti['Sample_only']:
    val = pd.isnull(ID)
    if val == True:
        continue
    else:
        #make a list of all full IDs to query the AMDB
        var_list.append(ID)
        ID = ID.replace("102.100.100/","")
        #make a list of short IDs to query the KEGG files
        IDs.append(str(ID))

#########
# Query the AMDB for the site metadata needed and merge it with the CI table
#########
returnfields=['sample_id','nrs_trip_code','nrs_sample_code','sample_site_location_description','latitude','longitude','collection_date','utc_time_sampled','utc_date_sampled','imos_site_code','depth','density','tot_depth_water_col','chlorophyll_a','ammonium','total_co2','turbidity','alkalinity','secchi_depth','picoeukaryotes','prochlorococcus','synechococcus']
searchfield = 'sample_id'
selectColumns_by_fieldValues(returnfields, searchfield, var_list)
#make the NRS TripCode_depth from the nrs_sample_code
for i in selectFields.index:
    null_val = pd.isnull(selectFields.loc[i,'nrs_sample_code'])
    if null_val == False:
        selectFields.loc[i,'TripCode_depth'] = selectFields.loc[i,'nrs_sample_code'].rsplit('_', 1)[1]
#Make the individual Date columns
selectFields['utc_date_sampled'] = pd.to_datetime(selectFields['utc_date_sampled'], format='%Y-%m-%d')
selectFields['Year'] = selectFields['utc_date_sampled'].dt.year.astype("Int64")
selectFields['Month'] = selectFields['utc_date_sampled'].dt.month.astype("Int64")
selectFields['Day'] = selectFields['utc_date_sampled'].dt.day.astype("Int64")
#make the StationCode
selectFields['imos_site_code'].fillna(np.nan,inplace=True)
selectFields['StationCode'] = selectFields['imos_site_code'].str.replace('NRS','')
#merge the metadata with the indecies
all_cti = pd.merge(all_cti,selectFields, left_on='Sample_only', right_on='sample_id', how='left')

#########
# Find the Kegg Terms defined in the gene_request file
#########
pathway_terms_df = pd.read_csv('gene_request.csv')
pathway_terms_df.columns= pathway_terms_df.columns.str.strip().str.lower()
pathway_terms = dict(zip(pathway_terms_df.kegg,pathway_terms_df.gene))
#Create a minimal dictionary to hold all the keg terms we find along with sample info - we will append to it as we go along
KegTerms = {'Sample_only':[]}
for key in pathway_terms:
    print(key)
    KegTerms[key] = []
for ID in IDs:
    print(ID)
    # the Rscript has cast the sample IDs as floats so we will drop the decimal
    ID = str(ID).replace('.0','')
    #read a local tar file
    tar_file_name = "/datasets/work/oa-amd/work/amd-work/SQM_READS/" + str(ID) + ".tar.gz"
    if os.path.exists(tar_file_name):
        #write the kegg file used to the metadata file
        tar_f = os.path.basename(tar_file_name)
        write_files_used(tar_f)
        map_stat = str(ID) + "/" + str(ID) +"_MGSD_CSIRO.sqmreads.out.mappingstat"
        #read the Mappingstat file within the tarball to get the toal number of reads proceesed
        with tarfile.open(tar_file_name, "r:gz") as tar:
            with tar.extractfile(map_stat) as f:
                for line in f:
                    line = line.decode("utf-8")
                    if '#' in line:
                        continue
                    else:
                        Sample,File,TotalReads,Reads_hit_nr =  line.split('\t')
                        total_reads = str(TotalReads)
        data_file_name = str(ID) + "/" + str(ID) +"_MGSD_CSIRO.sqmreads.out.allreads.funkegg"
        #Read the samples kegg file within the tarball to get the requested patways
        with tarfile.open(tar_file_name, "r:gz") as tar:
            with tar.extractfile(data_file_name) as f:
                #set up a tmp dictionary for the sample we will add to it if it finds a kegg term in the dict
                tmp = dict()
                for line in f:
                    line = line.decode("utf-8")
                    #skip the headers
                    if line.startswith("#"):
                        continue
                    keggT,Total,sampleAbund,Function,Class = line.split('\t')
                    #check if the kegg terms exist per line
                    for key in pathway_terms:
                        if key == keggT:
                            tmp.setdefault(key, []).append(Total)
                            continue
                #if we have found kegg terms we will add em to the master dict
                if len(tmp) > 0:
                    print(tmp)
                    KegTerms.setdefault('Sample_only',[]).append(ID)
                    KegTerms.setdefault('Total_reads',[]).append(total_reads)
                    #account for keggs that wernt found else write them to the master dict
                    for key in pathway_terms:
                        if key not in tmp:
                            val = np.nan
                        else:
                            val = ''.join(tmp[key])
                        KegTerms.setdefault(key,[]).append(val)
    else:
        print(str(ID) + ' No Metagenome')
Keg_terms_df = pd.DataFrame.from_dict(KegTerms)
for i in Keg_terms_df.index:
    for key in pathway_terms:
        val = pd.isnull(Keg_terms_df.loc[i, key])
        if val == False:
            standardized_name = key + "_per_mil_reads"
            Keg_terms_df.loc[i, standardized_name] = int(Keg_terms_df.loc[i, key]) / (int(Keg_terms_df.loc[i, 'Total_reads']) / 1000000)
print(Keg_terms_df.shape)
Keg_terms_df['Sample_only'] = "102.100.100/" + Keg_terms_df['Sample_only'].astype(str)
all_cti = pd.merge(all_cti, Keg_terms_df, left_on='Sample_only', right_on='Sample_only', how='left')
print(all_cti.shape)
all_cti.head()

#########
# Combine and zip the ASV statistics files
#########
for amplicon in amplicons:
    stats = glob.glob(os.path.join(os.getcwd(),amplicon + '*_statistics.csv'))
    #create empty dataframe to merge the results into
    stats_merge = pd.DataFrame(columns=['ASV'])
    for f in stats:
        df = pd.read_csv(f)
        df.drop('Unnamed: 0',axis=1, inplace=True)
        print(df.shape)
        stats_merge = pd.merge(stats_merge, df, left_on='ASV',right_on='ASV',how='outer')
    merged_stats_file = "data/" + amplicon + "_Summary_combined_statistics.csv.gz"
    stats_merge.to_csv(merged_stats_file, index=False, compression="gzip")

#########
# Move the metadata file containing the AM files used to the data folder
#########
shutil.move('AM_input_file_metadata.txt', 'data/AM_input_file_metadata.txt')

#########
# change the AM column names to the IMOS_visualisation columns 
#########
#Dictionary to map names 
#- dictionary should be structured with the column order that IMOS wants
#- dict structure should be <AM_column_name>:<IMOS_column name> - known mappings so far are below 
#- all columns to be retained in final table should be listed - even if the key:val are the same
column_dict = {'Sample_only':'code','TripCode_depth':'TripCode_depth','nrs_trip_code':'TripCode',\
'sample_site_location_description':'StationName','StationCode':'StationCode',\
'longitude':'Longitude','latitude':'Latitude','collection_date':'SampleDateUTC',\
'utc_time_sampled':'Time_24hr','Year':'Year','Month':'Month','Day':'Day',\
'Silicate':'Silicate_umolL','Nitrogen':'Nitrate_umolL','Phosphate':'Phosphate_umolL',\
'ammonium':'Ammonium_umolL','chlorophyll_a':'Chla_mgm3',\
'Temperature':'Temperature_degC','Salinity':'Salinity_psu',\
'Oxygen':'Oxygen_umolL','depth':'depth_m','turbidity':'Turbidity_NTU',\
'density':'Density_kgm3','tot_depth_water_col':'bottom_depth',\
'secchi_depth':'Secchi_m','total_co2':'DIC_umolkg',\
'alkalinity':'TAlkalinity_umolkg','prochlorococcus':'Prochlorococcus_cells_ml',\
'synechococcus':'Synechococcus_cells_ml','picoeukaryotes':'Picoeukaryotes_cells_ml',\
'K02588_per_mil_reads':'NifH_genes_per_mil_reads',\
'K01601_per_mil_reads':'RuBisCo_genes_per_mil_reads',\
'fish_parasites':'fish_parasites',\
'nitrogen_fixation':'nitrogen_fixation_organisms',\
'16S CSI_kernaldensity':'Bacterial_Salinity_Index_KD',\
'16S CSI_mean':'Bacterial_Salinity_Index_mean',\
'16S CSI_range':'Bacterial_Salinity_Range',\
'16S CSI_proportion':'Bacterial_Salinity_Proportion',\
'16S CSI_bias':'Bacterial_Salinity_Bias',\
'16S CSI_diversity':'Bacterial_Salinity_Diversity',\
'A16S CSI_kernaldensity':'Archaeal_Salinity_Index_KD',\
'A16S CSI_mean':'Archaeal_Salinity_Index_mean',\
'A16S CSI_range':'Archaeal_Salinity_Range',\
'A16S CSI_proportion':'Archaeal_Salinity_Proportion',\
'A16S CSI_bias':'Archaeal_Salinity_Bias',\
'A16S CSI_diversity':'Archaeal_Salinity_Diversity',\
'18Sv4 CSI_kernaldensity':'Eukaryote_Salinity_Index_KD',\
'18Sv4 CSI_mean':'Eukaryote_Salinity_Index_mean',\
'18Sv4 CSI_range':'Eukaryote_Salinity_Range',\
'18Sv4 CSI_proportion':'Eukaryote_Salinity_Proportion',\
'18Sv4 CSI_bias':'Eukaryote_Salinity_Bias',\
'18Sv4 CSI_diversity':'Eukaryote_Salinity_Diversity',\
'16S CNI_kernaldensity':'Bacterial_Nitrogen_Index_KD',\
'16S CNI_mean':'Bacterial_Nitrogen_Index_mean',\
'16S CNI_range':'Bacterial_Nitrogen_Range',\
'16S CNI_proportion':'Bacterial_Nitrogen_Proportion',\
'16S CNI_bias':'Bacterial_Nitrogen_Bias',\
'16S CNI_diversity':'Bacterial_Nitrogen_Diversity',\
'A16S CNI_kernaldensity':'Archaeal_Nitrogen_Index_KD',\
'A16S CNI_mean':'Archaeal_Nitrogen_Index_mean',\
'A16S CNI_range':'Archaeal_Nitrogen_Range',\
'A16S CNI_proportion':'Archaeal_Nitrogen_Proportion',\
'A16S CNI_bias':'Archaeal_Nitrogen_Bias',\
'A16S CNI_diversity':'Archaeal_Nitrogen_Diversity',\
'18Sv4 CNI_kernaldensity':'Eukaryote_Nitrogen_Index_KD',\
'18Sv4 CNI_mean':'Eukaryote_Nitrogen_Index_mean',\
'18Sv4 CNI_range':'Eukaryote_Nitrogen_Range',\
'18Sv4 CNI_proportion':'Eukaryote_Nitrogen_Proportion',\
'18Sv4 CNI_bias':'Eukaryote_Nitrogen_Bias',\
'18Sv4 CNI_diversity':'Eukaryote_Nitrogen_Diversity',\
'16S CPI_kernaldensity':'Bacterial_Phosphate_Index_KD',\
'16S CPI_mean':'Bacterial_Phosphate_Index_mean',\
'16S CPI_range':'Bacterial_Phosphate_Range',\
'16S CPI_proportion':'Bacterial_Phosphate_Proportion',\
'16S CPI_bias':'Bacterial_Phosphate_Bias',\
'16S CPI_diversity':'Bacterial_Phosphate_Diversity',\
'A16S CPI_kernaldensity':'Archaeal_Phosphate_Index_KD',\
'A16S CPI_mean':'Archaeal_Phosphate_Index_mean',\
'A16S CPI_range':'Archaeal_Phosphate_Range',\
'A16S CPI_proportion':'Archaeal_Phosphate_Proportion',\
'A16S CPI_bias':'Archaeal_Phosphate_Bias',\
'A16S CPI_diversity':'Archaeal_Phosphate_Diversity',\
'18Sv4 CPI_kernaldensity':'Eukaryote_Phosphate_Index_KD',\
'18Sv4 CPI_mean':'Eukaryote_Phosphate_Index_mean',\
'18Sv4 CPI_range':'Eukaryote_Phosphate_Range',\
'18Sv4 CPI_proportion':'Eukaryote_Phosphate_Proportion',\
'18Sv4 CPI_bias':'Eukaryote_Phosphate_Bias',\
'18Sv4 CPI_diversity':'Eukaryote_Phosphate_Diversity',\
'16S CSiI_kernaldensity':'Bacterial_Silicate_Index_KD',\
'16S CSiI_mean':'Bacterial_Silicate_Index_mean',\
'16S CSiI_range':'Bacterial_Silicate_Range',\
'16S CSiI_proportion':'Bacterial_Silicate_Proportion',\
'16S CSiI_bias':'Bacterial_Silicate_Bias',\
'16S CSiI_diversity':'Bacterial_Silicate_Diversity',\
'A16S CSiI_kernaldensity':'Archaeal_Silicate_Index_KD',\
'A16S CSiI_mean':'Archaeal_Silicate_Index_mean',\
'A16S CSiI_range':'Archaeal_Silicate_Range',\
'A16S CSiI_proportion':'Archaeal_Silicate_Proportion',\
'A16S CSiI_bias':'Archaeal_Silicate_Bias',\
'A16S CSiI_diversity':'Archaeal_Silicate_Diversity',\
'18Sv4 CSiI_kernaldensity':'Eukaryote_Silicate_Index_KD',\
'18Sv4 CSiI_mean':'Eukaryote_Silicate_Index_mean',\
'18Sv4 CSiI_range':'Eukaryote_Silicate_Range',\
'18Sv4 CSiI_proportion':'Eukaryote_Silicate_Proportion',\
'18Sv4 CSiI_bias':'Eukaryote_Silicate_Bias',\
'18Sv4 CSiI_diversity':'Eukaryote_Silicate_Diversity',\
'16S COI_kernaldensity':'Bacterial_Oxygen_Index_KD',\
'16S COI_mean':'Bacterial_Oxygen_Index_mean',\
'16S COI_range':'Bacterial_Oxygen_Range',\
'16S COI_proportion':'Bacterial_Oxygen_Proportion',\
'16S COI_bias':'Bacterial_Oxygen_Bias',\
'16S COI_diversity':'Bacterial_Oxygen_Diversity',\
'A16S COI_kernaldensity':'Archaeal_Oxygen_Index_KD',\
'A16S COI_mean':'Archaeal_Oxygen_Index_mean',\
'A16S COI_range':'Archaeal_Oxygen_Range',\
'A16S COI_proportion':'Archaeal_Oxygen_Proportion',\
'A16S COI_bias':'Archaeal_Oxygen_Bias',\
'A16S COI_diversity':'Archaeal_Oxygen_Diversity',\
'18Sv4 COI_kernaldensity':'Eukaryote_Oxygen_Index_KD',\
'18Sv4 COI_mean':'Eukaryote_Oxygen_Index_mean',\
'18Sv4 COI_range':'Eukaryote_Oxygen_Range',\
'18Sv4 COI_proportion':'Eukaryote_Oxygen_Proportion',\
'18Sv4 COI_bias':'Eukaryote_Oxygen_Bias',\
'18Sv4 COI_diversity':'Eukaryote_Oxygen_Diversity',\
'16S CTI_kernaldensity':'Bacterial_Temperature_Index_KD',\
'16S CTI_mean':'Bacterial_Temperature_Index_mean',\
'16S CTI_range':'Bacterial_Temperature_Range',\
'16S CTI_proportion':'Bacterial_Temperature_Proportion',\
'16S CTI_bias':'Bacterial_Temperature_Bias',\
'16S CTI_diversity':'Bacterial_Temperature_Diversity',\
'A16S CTI_kernaldensity':'Archaeal_Temperature_Index_KD',\
'A16S CTI_mean':'Archaeal_Temperature_Index_mean',\
'A16S CTI_range':'Archaeal_Temperature_Range',\
'A16S CTI_proportion':'Archaeal_Temperature_Proportion',\
'A16S CTI_bias':'Archaeal_Temperature_Bias',\
'A16S CTI_diversity':'Archaeal_Temperature_Diversity',\
'18Sv4 CTI_kernaldensity':'Eukaryote_Temperature_Index_KD',\
'18Sv4 CTI_mean':'Eukaryote_Temperature_Index_mean',\
'18Sv4 CTI_range':'Eukaryote_Temperature_Range',\
'18Sv4 CTI_proportion':'Eukaryote_Temperature_Proportion',\
'18Sv4 CTI_bias':'Eukaryote_Temperature_Bias',\
'18Sv4 CTI_diversity':'Eukaryote_Temperature_Diversity',\
'Bacteria_unique_ASVs':'Bacteria_unique_ASVs',\
'Bacteria_shannon_index':'Bacteria_shannon_index',\
'Bacteria_simpsons_index':'Bacteria_simpsons_index',\
'Bacteria_invsimpson_index':'Bacteria_invsimpson_index',\
'Bacteria_total_observations':'Bacteria_total_observations',\
'Archaea_unique_ASVs':'Archaea_unique_ASVs',\
'Archaea_shannon_index':'Archaea_shannon_index',\
'Archaea_simpsons_index':'Archaea_simpsons_index',\
'Archaea_invsimpson_index':'Archaea_invsimpson_index',\
'Archaea_total_observations':'Archaea_total_observations',\
'Eukaryote_unique_ASVs':'Eukaryote_unique_ASVs',\
'Eukaryote_shannon_index':'Eukaryote_shannon_index',\
'Eukaryote_simpsons_index':'Eukaryote_simpsons_index',\
'Eukaryote_invsimpson_index':'Eukaryote_invsimpson_index',\
'Eukaryote_total_observations':'Eukaryote_total_observations'}
#drop any columns not recognised by the dictionary
#make a list of known columns that map to the IMOS columns
known_col_variants = column_dict.keys()
#make a list of the df columns
incoming_cols = all_cti.columns.tolist()
#Make a list columns not in the known variants 
drops = [col for col in incoming_cols if col not in known_col_variants]
all_cti.drop(drops,axis=1, inplace=True)
#reorder the columns
all_cti = all_cti[known_col_variants]
#map the column names to the IMOS format
all_cti.columns = all_cti.columns.to_series().map(column_dict)

#########
# Write the final output file- we will overwrite the file in the data dir
#########
all_cti.to_csv('data/oceanViz_AM_data.csv', index=False)

#########
# Delete unwanted files
#########
unwanted_files = ['*_infile_multi_params_all.csv', '*I.csv','MOA_Diversity_stats.csv','*_statistics.csv','parameters.txt']
for file in unwanted_files:
    for f in glob.glob(file):
        os.remove(f)
