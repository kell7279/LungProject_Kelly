### Lung Analyses Without Utilizing TCGA Data for Feature Selection ###

dataDir <- "~/Dropbox (Partners HealthCare)/Lung CTC Project/Lung CTC Samples/LungProject_Kelly/Data"
#Save figures here:
figDir <- "~/Dropbox (Partners HealthCare)/Lung CTC Project/Lung CTC Samples/LungProject_Kelly/Output/Figures"
#Save R objects (Random forest CV output) here:
outDir <- "~/Dropbox (Partners HealthCare)/Lung CTC Project/Lung CTC Samples/LungProject_Kelly/Output/R objects"


### Load data and limit to healthy smokers (NOT high risk) and adeno lung cancer patents (NOT metastatic or Squamous)
setwd(dataDir)
lookup <- readRDS("lookup_k.RDS")
logRPM <- readRDS("log2RPM.RDS")

samps <- logRPM[lookup$EXCLUDE == FALSE,]
adeno <- lookup$Adeno[lookup$EXCLUDE==FALSE]

## RF classifier without PCA or TCGA
source("~/Documents/Martin Aryee/HCC:CLD/Analysis/CV.R")
set.seed(7279)
noPCA_noTCGA_ttest <- HCC_cv_ttest(data = samps, cancer=adeno, reps=10)
#mean AUC=0.70

## RF classifier with PCA, no TCGA
set.seed(7279)
PCA_noTCGA_ttest_3 <- HCC_PCA_ttest(data = samps, cancer = adeno, reps = 10, folds = 10, num_PC=3)

#Plot average ROC curves
png(paste0(figDir,"/Lung_noTCGA.png"), width = 6, height = 4, units = 'in', res = 300)
plot.roc(adeno, rowMeans(noPCA_noTCGA_ttest$votes), col="black", asp=NA, legacy.axes=TRUE, main="Lung - No TCGA")
lines.roc(adeno, rowMeans(PCA_noTCGA_ttest_3$votes), col="red")
lines.roc(adeno, rowMeans(PCA_noTCGA_ttest_3$votes_adj), col="green")
legend("bottomright",  title="Analysis",
       c(paste("No PCA",round(mean(noPCA_noTCGA_ttest$auc),3)),paste("PCA Feature Selection", round(mean(PCA_noTCGA_ttest_3$auc),3)), paste("PCA Adjustment (3 PCs)", round(mean(PCA_noTCGA_ttest_3$auc_adj),3))), col=c("black", "red", "green"), lty = 1, lwd=2)      
dev.off()

## Save CV results
saveRDS(noPCA_noTCGA_ttest, paste0(outDir,"/noPCA_noTCGA_RFcv.RDS"))
saveRDS(PCA_noTCGA_ttest_3, paste0(outDir, "/PCA3_noTCGA_RFcv.RDS"))

