#!/usr/bin/env Rscript

#TO COME .... Script description.....

## -----------------------------------------------------------------------------
library(tidyverse)
library(data.table)
write("defined libs", stdout())
## -----------------------------------------------------------------------------
#read in table containing  20K samples x abundance with metadata
print("using amplicon")
amplicon = commandArgs(trailingOnly=TRUE)[1]
print(amplicon)

write("read files", stdout())
input_file <-  paste(amplicon,"*_multi_params_all.csv", sep = "")
#input <- list.files(pattern=glob2rx("*_multi_params_all.csv"),full.names = T) 
input <- list.files(pattern=glob2rx(input_file),full.names = T) 
print(input)
#split the inputfile name on the first uderscore delimiter to get the amplicon being analysed

#grab the array number passed to the script, we will use it to append to the output files
write("array_num",stdout())
array_num = commandArgs(trailingOnly=TRUE)[2]
print(array_num)

rawdata <- read_csv(input)
rawdata  <- rename(rawdata, ASV = "#OTU ID") # rename #OTU ID without the #
rawdata$Sample_only <- as.integer(rawdata$Sample_only)
distinctASV <- rawdata %>% select(ASV) %>% distinct() #for later use building table

data <- rawdata %>% drop_na(Abundance_20K) #working with Abundance_20K data for this example so na.omit removes samples with <20k

#define the list of environmental parameters we will use, we will loop through them to produce the output files
params_list <- list("Nitrogen", "Oxygen", "Phosphate", "Salinity", "Silicate", "Temperature")
for (param in params_list){
    print(paste("Processing ", param))
    #Perform specific actions for each parameter
    if (param == 'Nitrogen'){
        data <- data %>% drop_na(`nitrate_nitrite`) %>% filter(!nitrate_nitrite < 0) %>% mutate(SST = round(`nitrate_nitrite`/0.01) * 0.01)
        data <- data %>% filter(data$nitrate_nitrite > -3) #remove samples where the AM sentry value has been entered as -9999
        #bin samples into  and find the number of samnples in the smallest bin
        write("bin", stdout())
        data_sub <- data %>% select(Sample_only, nitrate_nitrite) %>% distinct()
        data_sub <- data_sub %>% mutate(SST_bin = cut_interval(nitrate_nitrite, n = 5,labels = FALSE))
        cols <- c("ASV", "SNI_density","SNI_mean", "SNIRange","SNImax","SNImin")
    }
    if (param == 'Oxygen'){
        data <- data %>% drop_na(`oxygen`) %>% filter(!oxygen < 0) %>% mutate(SST = round(`oxygen`/0.01) * 0.01)
        data <- data %>% filter(data$oxygen > -3) #remove samples where the AM sentry value has been entered as -9999
        #bin samples into  and find the number of samnples in the smallest bin
        write("bin", stdout())
        data_sub <- data %>% select(Sample_only, oxygen) %>% distinct() %>% mutate(oxygen_group = case_when( `oxygen` < 100  ~ 100,`oxygen` >= 100 ~ `oxygen`)) #make new variable for sample selection to include low salniities into first group as there arent currently many
        data_sub <- data_sub %>% mutate(SST_bin = cut_interval(oxygen_group, n = 5,labels = FALSE))
        cols <- c("ASV", "SOI_density","SOI_mean", "SOIRange","SOImax","SOImin")
    }
    if (param == 'Phosphate'){
        data <- data %>% drop_na(`phosphate`) %>% filter(!phosphate < 0) %>% mutate(SST = round(`phosphate`/0.01) * 0.01)
        data <- data %>% filter(data$phosphate > -3) #remove samples where the AM sentry value has been entered as -9999
        
        #NOTE the removal of phosphate values above 3 is being performed due to high values in ~20 samples from Long Reef Port Phillip Bay
        #These were causing binning issues resulting in the script to crash
        #will need to sort out a better way to deal with outliers but this will require more resourcing 
        data <- data %>% filter(data$phosphate < 2.9) #remove phosphate values >3.2 as they are causing issues with binning
        #bin samples into  and find the number of samnples in the smallest bin
        write("bin", stdout())
        data_sub <- data %>% select(Sample_only, phosphate) %>% distinct()
        data_sub <- data_sub %>% mutate(SST_bin = cut_interval(phosphate, n = 5,labels = FALSE))
        cols <- c("ASV", "SPI_density","SPI_mean", "SPIRange","SPImax","SPImin")
    }
    if (param == "Salinity"){
        data <- data %>% drop_na(`salinity`) %>% filter(!salinity < 0) %>% mutate(SST = round(`salinity`/0.01) * 0.01)
        data <- data %>% filter(data$salinity > -3) #remove samples where the AM sentry value has been entered as -9999
        #bin samples into  and find the number of samnples in the smallest bin
        write("bin", stdout())
        data_sub <- data %>% select(Sample_only, salinity) %>% distinct() %>% mutate(salinity_group = case_when( `salinity` < 34  ~ 34,`salinity` >= 34 ~ `salinity`)) #make new variable for sample selection to include low salniities into first group
        data_sub <- data_sub %>% mutate(SST_bin = cut_interval(salinity_group, n = 5,labels = FALSE))
        cols <- c("ASV", "SSI_density","SSI_mean", "SSIRange","SSImax","SSImin")
    }
    if (param == "Silicate"){
        data <- data %>% drop_na(`silicate`) %>% filter(!silicate < 0) %>% mutate(SST = round(`silicate`/0.01) * 0.01)
        data <- data %>% filter(data$silicate > -3) #remove samples where the AM sentry value has been entered as -9999
        #bin samples  and find the number of samnples in the smallest bin
        write("bin", stdout())
        data_sub <- data %>% select(Sample_only, silicate) %>% distinct()
        data_sub <- data_sub %>% mutate(SST_bin = cut_interval(silicate, n = 5,labels = FALSE))
        cols <- c("ASV", "SSiI_density","SSiI_mean", "SSiIRange","SSiImax","SSiImin")
    }
    if (param == "Temperature"){
        data <- rawdata %>% drop_na(Abundance_20K) #working with Abundance_20K data for this example so na.omit removes samples with <20k
        data <- data %>% filter(data$temp > -3) #remove samples where the AM sentry value has been entered as -9999
        data <- data %>% mutate(SST = round(temp/0.2) * 0.2) #round temperatures to 0.2C increments to be used in kernal plots
        #bin samples into 3C bins and find the number of samnples in the smallest bin
        write("bin into 3 degs", stdout())
        data_sub <- data %>% select(Sample_only, temp) %>% distinct() #get sample_only and temperature and collapse table : results in 2545 samples with tempeature data
        data_sub <- data_sub %>% mutate(SST_bin = round(`temp`/3) * 3) #this actually works to create temperature bins of 3C....not sure why now...
        data_sub <- data_sub %>% mutate(SST_bin = case_when(SST_bin < 0 ~ 0, SST_bin >= 0 & SST_bin < 33 ~ SST_bin, SST_bin >= 33 ~ 30)) #move any tmperatures <0 into 0C bin and any temperature above 30 into 30C bin
        cols <- c("ASV", "STI_density","STI_mean", "STIRange","STImax","STImin")
    }
    #CREATE 10 binsv
    sample_count <- as.data.frame(table(data_sub$SST_bin))

    smin <- min(sample_count$Freq) # select smallest bin as variable for resampling

    #for loop to randomly resample smin samples from each bin and calculate kernal density for each OTU
    #set number of  replicates
    write("Samplng bins", stdout())
    r <- 1
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
        #colnames(sti_id)=c("ASV", "SNI_density","SNI_mean", "SNRange","SNmax","SNmin") # set column names in sti_id

        write("Set up kernal parameters", stdout())
      
        #set up kernal density parameters
        kernStep <- 0.1
        kernMin <- min(df$SST) - 3
        kernMax <- max(df$SST) + 3
        kernN <- round((kernMax - kernMin) / kernStep + 1)
        kernTemps <- seq(kernMin, kernMax, length.out=kernN)

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
            #mean approach averaging the temperatures of the 4 samples where ASV is most abdunant as sanity check
            meanData <- kernData %>% arrange(desc(abund))
            bins <- head(meanData,4)
            meanOut <- mean(bins$SST)

            #determine range based on level at which kernal density curve is quater peak height in min and max direction
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
        print(paste("Number of ASVs: ", nlevels(df$ASV)))
        sti_id <- as.data.frame(sti_id)

        #define output filename and write the file
        out_bootstrap <- paste(amplicon, "_sti_id_",param,"_subset",array_num,".csv", sep="")
        print(paste("Writing output file: ", out_bootstrap))
        #write(out_bootstrap, stdout())
        #Use the column names defined for each parameter to name the columns
        colnames(sti_id) <- cols
        write.csv(sti_id, file = out_bootstrap)
    }
}