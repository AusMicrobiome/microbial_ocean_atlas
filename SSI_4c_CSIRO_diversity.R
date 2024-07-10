#!/usr/bin/env Rscript

library(tidyverse)
write("defined libs", stdout())

write("read files", stdout())
Arch_input <- list.files(pattern=glob2rx("A16S*_infile_multi_params_all.csv"),full.names = T)
print(Arch_input)
Archaea_rawdata <- read_csv(Arch_input) #read in table containing pelagic archaeal 20K samples x abundance with temperature
Archaea_rawdata  <- rename(Archaea_rawdata, ASV = "#OTU ID")

Bact_input <- list.files(pattern=glob2rx("16S*_infile_multi_params_all.csv"),full.names = T)
print(Bact_input)
Bacteria_rawdata <- read_csv(Bact_input)
Bacteria_rawdata <- rename(Bacteria_rawdata, ASV = "#OTU ID")

Euk_input <- list.files(pattern=glob2rx("18Sv4*_infile_multi_params_all.csv"),full.names = T)
print(Euk_input)
Eukaryote_rawdata <- read_csv(Euk_input)
Eukaryote_rawdata <- rename(Eukaryote_rawdata, ASV = "#OTU ID")


write("define functions", stdout())
richness <- function(x){
  sum(x > 0)
}

shannon <- function(x){
  rabund  <-  x[x>0]/sum(x)
  -sum(rabund  * log(rabund))
}

simpson <- function(x){
  n <-  sum(x)
  sum(x *  (x-1) / (n * (n -1)))
}




write("calculate diversity", stdout())
Archaea_Diversity <- Archaea_rawdata %>%  group_by(Sample_only) %>% summarize(Archaea_unique_ASVs =  richness(Abundance), 
                                                                             Archaea_shannon_index  =  shannon(Abundance),
                                                                             Archaea_simpsons_index  =  simpson(Abundance),
                                                                             Archaea_invsimpson_index  =  1/(Archaea_simpsons_index),
                                                                             Archaea_total_observations = sum(Abundance))

Bacteria_Diversity <- Bacteria_rawdata %>%  group_by(Sample_only) %>% summarize(Bacteria_unique_ASVs =  richness(Abundance), 
                                                                                Bacteria_shannon_index  =  shannon(Abundance),
                                                                                Bacteria_simpsons_index  =  simpson(Abundance),
                                                                                Bacteria_invsimpson_index  =  1/(Bacteria_simpsons_index),
                                                                                Bacteria_total_observations = sum(Abundance))

Eukaryote_Diversity <- Eukaryote_rawdata %>%  group_by(Sample_only) %>% summarize(Eukaryote_unique_ASVs =  richness(Abundance), 
                                                                                  Eukaryote_shannon_index  =  shannon(Abundance),
                                                                                  Eukaryote_simpsons_index  =  simpson(Abundance),
                                                                                  Eukaryote_invsimpson_index  =  1/(Eukaryote_simpsons_index),
                                                                                  Eukaryote_total_observations = sum(Abundance))


All_Diversity <-  Bacteria_Diversity %>% full_join(Archaea_Diversity, by= "Sample_only")
All_Diversity <- All_Diversity %>%  full_join(Eukaryote_Diversity, by = "Sample_only")

write.csv(All_Diversity, "MOA_Diversity_stats.csv")


