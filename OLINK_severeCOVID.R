library(tidyr)
library(OlinkAnalyze)
#---Read in the OLINK NPX file
olinkdata_npx <- read.table("/MGH_COVID_OLINK_NPX.txt",
                         sep = ";", h = TRUE)
olinkdata_npx <- as.data.frame(olinkdata_npx)

#---Read in the clinical info
clinical <- read.table("/MGH_COVID_Clinical_Info.txt",
                       sep = ";", h = TRUE)

#---Remove ctrl samples
olinkdata_npx <- subset(olinkdata_npx, olinkdata_npx$subject_id %in% clinical$subject_id)

#---Retain only samples collected on D0
olinkdata_npx <- subset(olinkdata_npx, olinkdata_npx$Timepoint=="D0")

#---PCA biplot to examine outliers
olink_pca_plot(olinkdata_npx, outlierDefX = 3, outlierDefY = 3, 
               outlierLines = TRUE, label_outliers = TRUE)

#--Retain only the Inflammation panel for CXCL8, TNF, and IL8
olinkdata_npx <- subset(olinkdata_npx, !(Assay == "CXCL8" & Panel != "Inflammation"))
olinkdata_npx <- subset(olinkdata_npx, !(Assay == "IL6" & Panel != "Inflammation"))
olinkdata_npx <- subset(olinkdata_npx, !(Assay == "TNF" & Panel != "Inflammation"))

#---Format the NPX dataframe
olinkdata_npx <- subset(olinkdata_npx, select=c("subject_id", "Assay", "NPX"))
olinkdata_npx <- pivot_wider(olinkdata_npx, names_from=Assay, values_from=NPX)
colnames(olinkdata_npx) <- gsub("-", "_", colnames(olinkdata_npx))
#---Format the sample table dataframe TO SELECT COVARIATES
clinical <- subset(clinical, select =c("subject_id", "COVID", "Age_cat", "HEART",
                                       "DIABETES", "HTN", "KIDNEY", "IMMUNO", "LUNG"))

#---Merge the two dataframes
merged_df <- merge(olinkdata_npx, clinical, by="subject_id")
merged_df$COVID <- as.factor(merged_df$COVID)

#---Linear model to assess abundance
proteins <- colnames(olinkdata_npx)
proteins <- proteins[-1]
results_lm <- as.data.frame(proteins)
results_lm$pvalue <- rep("0", nrow(results_lm))
for(i in 1:length(proteins))
{
  model_lm <- summary(lm(paste(proteins[i],"~ COVID + Age_cat + HEART + DIABETES + HTN + KIDNEY + IMMUNO + LUNG"), data = merged_df))
  p_value <- model_lm$coefficients[2,4]
  results_lm$pvalue[i] <- as.numeric(p_value)
}
results_lm$pvalue <- as.numeric(results_lm$pvalue)
results_lm$FDR <- p.adjust(results_lm$pvalue, "BH")
