#!/usr/bin/env Rscript

#' ---
#' title: "R Notebook"
#' output: html_notebook
#' ---
#' 
#' This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 
#' 
#' Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Cmd+Shift+Enter*. 
#' 
## -----------------------------------------------------------------------------
library(tidyverse)
library(data.table)
write("defined libs", stdout())
#' 
#' 
#' 
#' 
## -----------------------------------------------------------------------------

#' 
#' 
## -----------------------------------------------------------------------------
write("read files", stdout())
temps_only_input <- list.files(pattern=glob2rx("*_temp.csv"),
         full.names = T)

amplicon <- gsub("./", "",temps_only_input)
amplicon <- strsplit(amplicon, "_")[[1]][1]
amplicon

rawdata <- read_csv(temps_only_input) #read in table containing pelagic archaeal 20K samples x abundance with temperature
rawdata  <- rename(rawdata, ASV = "#OTU ID") # rename #OTU ID without the #
distinctASV <- rawdata %>% select(ASV) %>% distinct() #for later use building table

data <- na.omit(rawdata) #working with Abundance_20K data for this example so na.omit removes samples with <20k
data <- data %>% filter(data$temp > -3) #remove samples where temperature has been entered as -9999

data <- data %>% mutate(SST = round(temp/0.2) * 0.2) #round temperatures to 0.2C increments to be used in kernal plots

data$Sample_only <- as.integer(data$Sample_only)




#' 
## -----------------------------------------------------------------------------
#bin samples into 3C bins and find the number of samnples in the smallest bin
write("bin into 3 degs", stdout())
data_sub <- data %>% select(Sample_only, temp) %>% distinct() #get sample_only and temperature and collapse table : results in 2545 samples with tempeature data
data_sub <- data_sub %>% mutate(SST_bin = round(`temp`/3) * 3) #this actually works to create temperature bins of 3C....not sure why now...
data_sub <- data_sub %>% mutate(SST_bin = case_when(SST_bin < 0 ~ 0, SST_bin >= 0 & SST_bin < 33 ~ SST_bin, SST_bin >= 33 ~ 30)) #move any tmperatures <0 into 0C bin and any temperature above 30 into 30C bin

sample_count <- as.data.frame(table(data_sub$SST_bin))

#sample_count
#   Var1 Freq
#1     0  148
#2     3  247
#3     6  149
#4     9  145
#5    12  216
#6    15  290
#7    18  362
#8    21  454
#9    24  231
#10   27  174
#11   30  129

smin <- min(sample_count$Freq) # select smallest bin as variable for resampling

#smin
#[1] 129

#' 
## -----------------------------------------------------------------------------
#for loop to randomly resample smin samples from each bin and calculate kernal density for each OTU
#set number of  replicates
write("resample bins", stdout())
r <- 100
for (p in 1:r) {
  
  
data_subset <- data_sub  %>% group_by(SST_bin) %>% sample_n(size = smin, replace=TRUE) #random sample from each bin 
select <- data_subset$Sample_only #read codes to selection variable

data_select <- data %>% filter(Sample_only %in% select) # select samples from main table based on select Sample Id 

x <- select(data_select, c(SST, ASV, Abundance_20K)) #pass columns temp, ASV and abund to x
x <- rename(x,abund = "Abundance_20K")

dfn <- x %>% group_by(ASV) %>% summarize(ntax = n()) # count the number of  samples each ASV is in in the subset

dfx <- inner_join(x, dfn, by="ASV") #join the count data back to table x

write("Setting abundance limits", stdout())
notax <- 20 # set number of occurrences - choose the threshold for how many samples an ASV must appear in before kernal density is  calculated
dfx <- subset(dfx, ntax>notax) #select only ASVs in > 20 samples


df <- dfx %>% group_by(SST, ASV) %>%
  summarize(abund = sum(abund), freq = n(), a = sum(abund)/n()) %>% droplevels() 

df$ASV <- factor(df$ASV, levels=unique(df$ASV))
  
sti_id <- matrix(0, nrow=nlevels(df$ASV), ncol=6) #set up a matrix to read results to - with the number of rows as ASVs and the number of columns as the indices being caluclated
colnames(sti_id)=c("ASV", "STI_density","STI_mean", "STRange","STmax","STmin") # set column names in sti_id

write("set up kernal parameters", stdout())
#set up kernal density parameters
kernStep <- 0.1
kernMin <- min(df$SST) - 3 
kernMax <- max(df$SST) + 3
kernN <- round((kernMax - kernMin) / kernStep + 1)
kernTemps <- seq(kernMin, kernMax, length.out=kernN)


#i=1
for (i in 1:nlevels(df$ASV)) {
  taxon <- levels(df$ASV)[i]
  kernData <- subset(df[,c(1:3,5)], ASV == taxon & abund > 0)
  kernData$ASV <- factor(kernData$ASV)
  kernData$weight <- with(kernData, abs(abund) / sum(abund))
  kernOut <- with(kernData,
                  density(SST, weight=weight,
                          bw="nrd0",
                          from=kernMin,
                          to=kernMax,
                          n=kernN))
  
  z <- cbind(kernTemps, kernOut$y)
  m <- as.data.frame(z)
  z <- z[which.max(z[,2]),] 
  write("averaging temps", stdout())
  #mean approach avareageing the temperatures of the 4 samples where ASV is most abdunant as sanity check
  meanData <- kernData %>% arrange(desc(abund))
  bins <- head(meanData,4)
  meanOut <- mean(bins$SST)
  
  #determine range based on level at which kernal desnioty curve is quater peak height in min and  max  direction 
  b = which(m$V2 > (as.numeric(z[2])/4))
  Tmin <- m$kernTemps[(head(b, n=1))]
  Tmax <- m$kernTemps[(tail(b, n=1))]
  rangeq4 <- (m$kernTemps[(tail(b, n=1))]) - (m$kernTemps[(head(b, n=1))])
   

  
  sti_id[i,1] <- paste0(taxon)
  sti_id[i,2] <- as.numeric(z[1])
  sti_id[i,3] <- as.numeric(meanOut)
  sti_id[i,4] <- as.numeric(rangeq4)
  sti_id[i,5] <- as.numeric(Tmax)
  sti_id[i,6] <- as.numeric(Tmin)
  
}

sti_id <- as.data.frame(sti_id)
#define output filename
out_bootstrap <- paste(amplicon, "_sti_id_Temperature_subset",p,".csv", sep="")
colnames(sti_id) <- c("ASV", "STI_density","STI_mean", "STRange","STmax","STmin")
write.csv(sti_id, file = out_bootstrap)
}

#' 
## -----------------------------------------------------------------------------
#I moved the .csv files into folder /Users/mvb815/MarineMicrobes Dropbox/AM_indices/STI/
#set up files for each index with ever ASV

temp_STI_d <- distinctASV
temp_STI_m <- distinctASV  
temp_STRq4 <- distinctASV 
temp_STmax <- distinctASV 
temp_STmin <- distinctASV 

#read files from STI directory and populate  indices
write("reading STI bins", stdout())
#files <- list.files(path = "./",
#         pattern="*subset*.csv", 
#         full.names = T)
         
files <- list.files(path = "./",
         pattern=glob2rx("*subset*.csv"),
         full.names = T)

write(files, stdout())
write("bootstrapping", stdout())
#i =1
for (i in 1:length(files))
{
  boot <- fread(files[i], header=T)
  boot <- as_tibble(boot)
  
boot_STI_d <- boot %>% select(ASV,STI_density)
boot_STI_m <- boot %>% select(ASV,STI_mean)
boot_STRq4 <- boot %>% select(ASV,STRange)
boot_STmax <- boot %>% select(ASV,STmax)
boot_STmin <- boot %>% select(ASV,STmin)
  
  
colnames(boot_STI_d)[2] <- paste(colnames(boot_STI_d[2]), "sub", i, sep = "_")
colnames(boot_STI_m)[2] <- paste(colnames(boot_STI_m[2]), "sub", i, sep = "_")
colnames(boot_STRq4)[2] <- paste(colnames(boot_STRq4[2]), "sub", i, sep = "_")
colnames(boot_STmax)[2] <- paste(colnames(boot_STmax[2]), "sub", i, sep = "_")
colnames(boot_STmin)[2] <- paste(colnames(boot_STmin[2]), "sub", i, sep = "_")

  
temp_STI_d <- left_join(temp_STI_d,boot_STI_d,by="ASV")
temp_STI_m <- left_join(temp_STI_m,boot_STI_m,by="ASV")
temp_STRq4 <- left_join(temp_STRq4,boot_STRq4,by="ASV")
temp_STmax <- left_join(temp_STmax,boot_STmax,by="ASV")
temp_STmin <- left_join(temp_STmin,boot_STmin,by="ASV")
  
  rm(boot)
}

#write.csv(temp_STI_d, "Archaea_Bootstraped_STI_density.csv")
#write.csv(temp_STI_m, "Archaea_Bootstraped_STI_mean.csv")
#write.csv(temp_STRq4, "Archaea_Bootstraped_STRange.csv")
#write.csv(temp_Tmax, "Archaea_Bootstraped_STmax.csv")
#write.csv(temp_Tmin, "Archaea_Bootstraped_STmin.csv")

#' 
#' 
#' 
#' 
## -----------------------------------------------------------------------------
#STI
#calculate mean, standard deviation, number of observations (count) and coefficient of variation (CV) across all subsampled replicates  
#r = number of replicate subsets run
write("calculate stats", stdout())
Summary <- distinctASV
Summary$STI_d_na <- apply(is.na(temp_STI_d[,-1]), 1, sum)
Summary$STI_d_count <- r - Summary$STI_d_na
Summary$STI_d_mean <- rowMeans(temp_STI_d[2:r + 1], na.rm=TRUE)
Summary$STI_d_SD <- apply(temp_STI_d[,-1],1,sd, na.rm=TRUE)
Summary$STI_d_CV <- Summary$STI_d_SD / Summary$STI_d_mean


Summary$STI_m_na <- apply(is.na(temp_STI_m[,-1]), 1, sum)
Summary$STI_m_count <- r - Summary$STI_m_na
Summary$STI_m_mean <- rowMeans(temp_STI_m[2:r + 1], na.rm=TRUE)
Summary$STI_m_SD <- apply(temp_STI_m[,-1],1,sd, na.rm=TRUE)
Summary$STI_m_CV <- Summary$STI_m_SD / Summary$STI_m_mean


Summary$STRq4_na <- apply(is.na(temp_STRq4[,-1]), 1, sum)
Summary$STRq4_count <- r - Summary$STRq4_na
Summary$STRq4_mean <- rowMeans(temp_STRq4[2:r +1], na.rm=TRUE)
Summary$STRq4_SD <- apply(temp_STRq4[,-1],1,sd, na.rm=TRUE)
Summary$STRq4_CV <- Summary$STRq4_SD / Summary$STRq4_mean

Summary$STmax_na <- apply(is.na(temp_STmax[,-1]), 1, sum)
Summary$STmax_count <- r - Summary$STmax_na
Summary$STmax_mean <- rowMeans(temp_STmax[2:r + 1], na.rm=TRUE)
Summary$STmax_SD <- apply(temp_STmax[,-1],1,sd, na.rm=TRUE)
Summary$STmax_CV <- Summary$STmax_SD / Summary$STmax_mean

Summary$STmin_na <- apply(is.na(temp_STmin[,-1]), 1, sum)
Summary$STmin_count <- r - Summary$STmax_na
Summary$STmin_mean <- rowMeans(temp_STmin[2:r + 1], na.rm=TRUE)
Summary$STmin_SD <- apply(temp_STmin[,-1],1,sd, na.rm=TRUE)
Summary$STmin_CV <- Summary$STmin_SD / Summary$STmin_mean

out_stats <- paste(amplicon,"_Summary_STI_statistics.csv", sep="")
write.csv(Summary, out_stats)
#write.csv(Summary, "Archaea_Summary_STI_statistics.csv")


#' 
#' Create community incices 
#' 
## -----------------------------------------------------------------------------
write("create community indecies", stdout())
Summary_sub <- Summary %>% select(ASV,STI_d_mean,STI_m_mean,STRq4_mean,STmax_mean,STmin_mean)          

#use raw aundance rather than  Abundance_20k as it  becomes  proportional and we can use more data 
temp_notemp_input <- list.files(pattern=glob2rx("*_temp_notemp.csv"),
         full.names = T)
alldata <- read_csv(temp_notemp_input) #read in table containing pelagic archaeal 20K samples x abundance with temperature
alldata  <- rename(alldata, ASV = "#OTU ID") # rename #OTU ID without the #
all <- alldata %>% select(ASV, Sample_only, Abundance)
all <-  rename(all, abund = Abundance)

all <- all %>% group_by(Sample_only) %>% mutate(total = sum(abund))  #get the total sequence count in order to calculate proportion of sample that had STI used in the CTI calcualtion

all <- all %>% inner_join(Summary_sub, by="ASV") %>% drop_na(STI_d_mean) %>% droplevels() #remove ASVs with no  STI

#calculate community indices
cti_all <- all %>% group_by(Sample_only) %>% summarise(CTI_d = sum(abund*STI_d_mean)/sum(abund), 
                                                CTI_m = sum(abund*STI_m_mean)/sum(abund),
                                                CTRq4 = sum(abund*STRq4_mean)/sum(abund),
                                                proportion_T = sum(abund)/total)

cti_all <- cti_all %>% distinct()

#calculate community temperature bias = difference between CTI and sample temperature
metadata <- alldata %>% select(Sample_only, temp) %>% distinct()
cti_all <- inner_join(cti_all, metadata, by = c("Sample_only"))
cti_all$CT_bias <- cti_all$CTI_d - cti_all$temp

#calculate community  thermal diversity
all <- inner_join(all, cti_all, by = c("Sample_only"))

cti_div <- all %>% group_by(Sample_only) %>% summarise(CTDiv = sqrt(sum(((STI_d_mean-CTI_d)*(STI_d_mean-CTI_d)*abund))/sum(abund)))

cti_all <- left_join(cti_all, cti_div, by = c("Sample_only"))

asv_cti <- paste(amplicon,"_CTI.csv", sep="")
write.csv(cti_all, asv_cti)


#' 
#' 
## -----------------------------------------------------------------------------


#' 
#' Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Cmd+Option+I*.
#' 
#' When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Cmd+Shift+K* to preview the HTML file). 
#' 
#' The preview shows you a rendered HTML copy of the contents of the editor. Consequently, unlike *Knit*, *Preview* does not run any R code chunks. Instead, the output of the chunk when it was last run in the editor is displayed.
#' 
#' 
