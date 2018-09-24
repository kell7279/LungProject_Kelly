###
# The purpose of this script is to create more user-friendly data files for 
# conducting analyses of the lung CTC data
###

library(readxl)
library(stringr)


dataDir <- "~/Dropbox (Partners HealthCare)/Lung CTC Project/Lung CTC Samples"
#Load my data + preprocess ---- (Code from Mark! Thanks Mark)
outDir <- "~/Dropbox (Partners HealthCare)/Lung CTC Project/Lung CTC Samples/Kelly/Data"

setwd(dataDir)
txt_files_ls = list.files(path=dataDir, pattern="*.txt") 
txt_files_df <- lapply(txt_files_ls, function(x) {read.table(file = x,row.names = 1)})
counts.df <- do.call("cbind", lapply(txt_files_df, as.data.frame))
counts.df<- counts.df[,!duplicated(colnames(counts.df))]
colnames(counts.df) <- str_split(colnames(counts.df),"_" ,n = 2,simplify = TRUE)[,1]


setwd(outDir)
lookup.df <- read.csv('NSCLC_Lookupedit.csv') #sample info sheet 
lookup.df$filename <- gsub('-','.',lookup.df$filename)

lookup.df <- lookup.df[order(match(lookup.df$filename,colnames(counts.df))),]
lookup.df$Smoker <- grepl("Healthy smoker", lookup.df$sample.type, ignore.case = T)
lookup.df$Smoker_HiRisk <- grepl("high risk", lookup.df$sample.type, ignore.case = T)
lookup.df$Cancer <- grepl("Lung Cancer", lookup.df$sample.type, ignore.case = T)
lookup.df$Metastatic <- grepl("metastatic", lookup.df$sample.type, ignore.case = T)
lookup.df$Adeno <- grepl("adeno", lookup.df$diagnosis, ignore.case=T)
sum(lookup.df$filename %in%  colnames(counts.df))

# Mark samples from same patient as "exclude"
colnames(counts.df)[grep("LS93", colnames(counts.df))]
colnames(counts.df)[grep("LS.110", colnames(counts.df))]
colnames(counts.df)[grep("S25", colnames(counts.df))]

lookup.df[lookup.df$filename %in% c("LS93.05.13.2016", "LS.110.06.02.2017", "S25.2"), ]$EXCLUDE <- TRUE

# Mark samples with low detected transcripts or total reads as "exclude"
source("~/Dropbox (Partners HealthCare)/Lung CTC Project/Lung CTC Samples/Kelly/Analyses/Supporting Scripts/RPM.R")
calcs <- RPM(counts.df)
detected_trans <- calcs$detected_trans
total_reads <- calcs$total_reads
log2RPM <- calcs$log2RPM

plot(total_reads, detected_trans)
# exlude samples with < 8000 detected transcripts
lookup.df[which(detected_trans<8000), "EXCLUDE"] <- TRUE
levels(lookup.df$Notes) <- c(levels(lookup.df$Notes), "low detected transcripts")
lookup.df[which(detected_trans<8000), ]$Notes <- "low detected transcripts"

# Mark squamous lung cancer samples as "exclude" since most are adeno 
# & they will have different gene expression profiles

lookup.df[grep("squamous", lookup.df$diagnosis), "EXCLUDE"] <- TRUE

### Save lookup, raw counts, and log2(1+RPM)

saveRDS(lookup.df, file = "lookup_k.RDS")
saveRDS(log2RPM, file = "log2RPM.RDS")
saveRDS(counts.df, file = "RawCounts.RDS")


