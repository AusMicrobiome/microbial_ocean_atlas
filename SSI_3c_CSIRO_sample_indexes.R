#!/usr/bin/env Rscript

########################
# Part2 Browns R script combine Tables from batch array
########################

## -----------------------------------------------------------------------------
library(tidyverse)
library(data.table)
write("defined libs", stdout())

## -----------------------------------------------------------------------------

print("using amplicon")
amplicon = commandArgs(trailingOnly=TRUE)[1]
print(amplicon)

input_file <-  paste(amplicon,"*_multi_params_all.csv", sep = "")

#read in the ASV table to get the ASV and amplicon info
input <- list.files(pattern=glob2rx(input_file),full.names = T) #read in table containing  20K samples x abundance with metadata
print(input)

rawdata <- read_csv(input) 
rawdata  <- rename(rawdata, ASV = "#OTU ID") # rename #OTU ID without the #
rawdata$Sample_only <- as.integer(rawdata$Sample_only)
distinctASV <- rawdata %>% select(ASV) %>% distinct() #for later use building table
#set up files for each index with ever ASV

tmp_d <- distinctASV
tmp_m <- distinctASV  
tmp_R <- distinctASV 
tmp_max <- distinctASV 
tmp_min <- distinctASV 

#grab the physchem passed to the script, we will use it to append to the output files
write("param",stdout())
param = commandArgs(trailingOnly=TRUE)[2]
print(param)

#params_list <- list("Nitrogen", "Oxygen", "Phosphate", "Salinity", "Silicate", "Temperature")
#for (param in params_list) {
if (param == 'Nitrogen'){
  short <- "SNI"
  cols <- c("ASV", "SNI_na", "SNI_count", "SNI_d_mean", "SNI_d_SD", "SNI_d_CV", "SNI_m_mean", "SNI_m_SD", "SNI_m_CV", "SNIR_mean", "SNIR_SD", "SNIR_CV", "SNImax_mean", "SNImax_SD", "SNImax_CV", "SNImin_mean", "SNImin_SD", "SNImin_CV")
  sub_cols <- c("ASV","SNI_d_mean","SNI_m_mean","SNIR_mean","SNImax_mean","SNImin_mean")
  commCI <- "CNI"
  env_col <- as.name("nitrate_nitrite")
}
if (param == 'Oxygen'){
  short <- "SOI"
    cols <- c("ASV", "SOI_na", "SOI_count", "SOI_d_mean", "SOI_d_SD", "SOI_d_CV", "SOI_m_mean", "SOI_m_SD", "SOI_m_CV", "SOIR_mean", "SOIR_SD", "SOIR_CV", "SOImax_mean", "SOImax_SD", "SOImax_CV", "SOImin_mean", "SOImin_SD", "SOImin_CV")
    sub_cols <- c("ASV","SOI_d_mean","SOI_m_mean","SOIR_mean","SOImax_mean","SOImin_mean")
    commCI <- "COI"
    env_col <- as.name("oxygen")    
}
if (param == 'Phosphate'){
  short <- "SPI"
  cols <- c("ASV", "SPI_na", "SPI_count", "SPI_d_mean", "SPI_d_SD", "SPI_d_CV", "SPI_m_mean", "SPI_m_SD", "SPI_m_CV", "SPIR_mean", "SPIR_SD", "SPIR_CV", "SPImax_mean", "SPImax_SD", "SPImax_CV", "SPImin_mean", "SPImin_SD", "SPImin_CV")
  sub_cols <- c("ASV","SPI_d_mean","SPI_m_mean","SPIR_mean","SPImax_mean","SPImin_mean")
  commCI <- "CPI"
  env_col <- as.name("phosphate") 
}
if (param == "Salinity"){
  short <- "SSI"
  cols <- c("ASV", "SPSI_na", "SSI_count", "SSI_d_mean", "SSI_d_SD", "SSI_d_CV", "SSI_m_mean", "SSI_m_SD", "SSI_m_CV", "SSIR_mean", "SSIOR_SD", "SSIR_CV", "SSImax_mean", "SSImax_SD", "SSImax_CV", "SSImin_mean", "SSImin_SD", "SSImin_CV")
  sub_cols <- c("ASV","SSI_d_mean","SSI_m_mean","SSIR_mean","SSImax_mean","SSImin_mean")
  commCI <- "CSI"
  env_col <- as.name("salinity")
}
if (param == "Silicate"){
  short = "SSiI"
  cols <- c("ASV", "SPSiI_na", "SSiI_count", "SSiI_d_mean", "SSiI_d_SD", "SSiI_d_CV", "SSiI_m_mean", "SSiI_m_SD", "SSiI_m_CV", "SSiIR_mean", "SSiIR_SD", "SSiIR_CV", "SSiImax_mean", "SSiImax_SD", "SSiImax_CV", "SSiImin_mean", "SSiImin_SD", "SSiImin_CV")
  sub_cols <- c("ASV","SSiI_d_mean","SSiI_m_mean","SSiIR_mean","SSiImax_mean","SSiImin_mean")
  commCI <- "CSiI"
  env_col <- as.name("silicate")
}
if (param == "Temperature"){
  short = "STI"
  cols <- c("ASV", "SPTI_na", "STI_count", "STI_d_mean", "STI_d_SD", "STI_d_CV", "STI_m_mean", "STI_m_SD", "STI_m_CV", "STIR_mean", "STIR_SD", "STIR_CV", "STImax_mean", "STImax_SD", "STImax_CV", "STImin_mean", "STImin_SD", "STImin_CV")
  sub_cols <- c("ASV","STI_d_mean","STI_m_mean","STIR_mean","STImax_mean","STImin_mean")
  commCI <- "CTI"
  env_col <- as.name("temp")
}
#read files from SNI directory and populate  indices
write("reading SNI bins", stdout())

files <- list.files(path = "./", pattern=glob2rx(paste(amplicon,"*", param,"_subset*.csv",sep="")), full.names = T)
    
#define the number of subsamples
r <- length(files)

write(files, stdout())
write("bootstrapping", stdout())
    
for (i in 1:length(files)){
  boot <- fread(files[i], header=T)
  boot <- as_tibble(boot)
      
  #define variables to hold contents of table columns
  boot_d <- boot %>% select(ASV,paste(short,"_density", sep=""))
  boot_m <- boot %>% select(ASV,paste(short,"_mean", sep=""))
  boot_R <- boot %>% select(ASV,paste(short,"Range", sep=""))
  boot_max <- boot %>% select(ASV,paste(short,"max", sep=""))
  boot_min <- boot %>% select(ASV,paste(short,"min", sep=""))
  
  #generate new column names to hold the incoming the concatenated table
  colnames(boot_d)[2] <- paste(colnames(boot_d[2]), "sub", i, sep = "_")
  colnames(boot_m)[2] <- paste(colnames(boot_m[2]), "sub", i, sep = "_")
  colnames(boot_R)[2] <- paste(colnames(boot_R[2]), "sub", i, sep = "_")
  colnames(boot_max)[2] <- paste(colnames(boot_max[2]), "sub", i, sep = "_")
  colnames(boot_min)[2] <- paste(colnames(boot_min[2]), "sub", i, sep = "_")
      
  #merging into the master table 
  tmp_d <- left_join(tmp_d,boot_d,by="ASV")
  tmp_m <- left_join(tmp_m,boot_m,by="ASV")
  tmp_R <- left_join(tmp_R,boot_R,by="ASV")
  tmp_max <- left_join(tmp_max,boot_max,by="ASV")
  tmp_min <- left_join(tmp_min,boot_min,by="ASV")
      
  rm(boot)
}

#calculate mean, standard deviation, number of observations (count) and coefficient of variation (CV) across all subsampled replicates  
#r = number of replicate subsets run
#calc St.dev
write("calculate stats", stdout())
Summary <- distinctASV
Summary$S_d_na <- apply(is.na(tmp_d[,-1]), 1, sum)
Summary$S_d_count <- r - Summary$S_d_na
Summary$S_d_mean <- rowMeans(tmp_d[2:r + 1], na.rm=TRUE)
Summary$S_d_SD <- apply(tmp_d[,-1],1,sd, na.rm=TRUE)
Summary$S_d_CV <- Summary$S_d_SD / Summary$S_d_mean

#calc mean
Summary$S_m_mean <- rowMeans(tmp_m[2:r + 1], na.rm=TRUE)
Summary$S_m_SD <- apply(tmp_m[,-1],1,sd, na.rm=TRUE)
Summary$S_m_CV <- Summary$S_m_SD / Summary$S_m_mean
    
#calc range
Summary$SR_mean <- rowMeans(tmp_R[2:r +1], na.rm=TRUE)
Summary$SR_SD <- apply(tmp_R[,-1],1,sd, na.rm=TRUE)
Summary$SR_CV <- Summary$SR_SD / Summary$SR_mean
    
#calc max
Summary$Smax_mean <- rowMeans(tmp_max[2:r + 1], na.rm=TRUE)
Summary$Smax_SD <- apply(tmp_max[,-1],1,sd, na.rm=TRUE)
Summary$Smax_CV <- Summary$Smax_SD / Summary$Smax_mean
    
#calc min
Summary$Smin_mean <- rowMeans(tmp_min[2:r + 1], na.rm=TRUE)
Summary$Smin_SD <- apply(tmp_min[,-1],1,sd, na.rm=TRUE)
Summary$Smin_CV <- Summary$Smin_SD / Summary$Smin_mean

#write to the output file
out_stats <- paste(amplicon,"_Summary_",short,"_statistics.csv", sep="")
colnames(Summary) <- cols
write.csv(Summary, out_stats)

## -----------------------------------------------------------------------------
write("create community indecies", stdout())
Summary_sub <- Summary %>% select(all_of(sub_cols))   

#use raw aundance rather than  Abundance_20k as it  becomes  proportional and we can use more data 
all <- rawdata %>% select(ASV, Sample_only, Abundance)
all <-  rename(all, abund = Abundance)
all$abund <- as.integer(all$abund)

all <- all %>% group_by(Sample_only) %>% mutate(total = sum(abund))  #get the total sequence count in order to calculate proportion of sample that had SNI used in the CTI calcualtion

all <- all %>% inner_join(Summary_sub, by="ASV") %>% drop_na(paste(short,"_d_mean",sep="")) %>% droplevels() #remove ASVs with no  SNI

head1 <- as.name(paste(short,"_d_mean",sep=""))
head2 <- as.name(paste(short,"_m_mean",sep=""))
head3 <- as.name(paste(short,"R_mean",sep=""))

#calculate community indices

ci_all <- all %>% group_by(Sample_only) %>% summarise(CI_d = sum(abund*!!head1)/sum(abund), 
                                                      CI_m = sum(abund*!!head2)/sum(abund),
                                                      CR = sum(abund*!!head3)/sum(abund),
                                                      proportion = sum(abund)/total) %>% distinct()
    

#calculate community nitrogen bias = difference between CNI and sample nitrite+nitrate
metadata <- rawdata %>% select(Sample_only, all_of(env_col)) %>% distinct()
ci_all <- inner_join(ci_all, metadata, by = c("Sample_only"))

ci_all  <- ci_all %>% mutate(CI_bias = CI_d - !!env_col)

#calculate community  thermal diversity
all <- inner_join(all, ci_all, by = c("Sample_only"))

ci_div <- all %>% group_by(Sample_only) %>% summarise(C_Div = sqrt(sum(((!!head1-CI_d)*(!!head1-CI_d)*abund))/sum(abund)))

ci_all <- left_join(ci_all, ci_div, by = c("Sample_only"))


a <- paste(amplicon, paste(commCI,"_kernaldensity", sep= ""))
b <- paste(amplicon, paste(commCI,"_mean", sep= ""))
c <- paste(amplicon, paste(commCI,"_range", sep= ""))
d <- paste(amplicon, paste(commCI,"_proportion", sep= ""))
e <- param
f <- paste(amplicon, paste(commCI,"_bias", sep= ""))
g <- paste(amplicon, paste(commCI,"_diversity", sep= ""))

colnames(ci_all) <- c("Sample_only", a, b, c, d, e,f,g)

asv_ci <-  paste(amplicon,"_",commCI,".csv", sep = "")
write.csv(ci_all, asv_ci)


