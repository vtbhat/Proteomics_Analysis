library(tidyr)
library(OlinkAnalyze)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggrepel)

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

#---Figure: Heatmap of top 200 differentially expressed proteins
results_lm <- results_lm[order(results_lm$FDR), ]
top200proteins <- results_lm$proteins[1:200]
top200npx_matrix <- olinkdata_npx
top200npx_matrix <- as.data.frame(top200npx_matrix)
rownames(top200npx_matrix) <- top200npx_matrix$subject_id
top200npx_matrix$subject_id <- NULL
top200npx_matrix <- as.matrix(subset(top200npx_matrix, select=top200proteins))
top200npx_matrix <- scale(top200npx_matrix)
top200npx_matrixt <- t(top200npx_matrix)
colnames(top200npx_matrixt) <- as.vector(olinkdata_npx$subject_id)
top200npx_matrixt <- as.data.frame(top200npx_matrixt)
top200npx_matrixt <- as.matrix(top200npx_matrixt)

#---Truncate values 4 or 5 SDs above mean
top200npx_matrixt <- pmax(pmin(top200npx_matrixt, 4), -4)

heatmap_ann <- subset(clinical, select=c("subject_id", "COVID"))
df_meta <- heatmap_ann[match(colnames(top200npx_matrixt), heatmap_ann$subject_id), ]

df_meta$COVID[df_meta$COVID==1] <- "COVID positive"
df_meta$COVID[df_meta$COVID==0] <- "COVID negative"

ha_col <- HeatmapAnnotation(
  COVID = df_meta$COVID,
  col = list(
    COVID = c("COVID negative" = "blue", "COVID positive" = "red")   # color mapping
  )
)

#---Heatmap
Heatmap(top200npx_matrixt, top_annotation = ha_col, 
        heatmap_legend_param=list(title="Z Scores"), show_row_names=FALSE, 
        show_column_names = FALSE)

#---Volcano Plot
res_volcano <- results_lm
res_volcano$DiffExpressed <- "Not Significant"
for(i in 1:nrow(res_volcano))
{
  if(res_volcano$npx_diff[i]<(0) && res_volcano$FDR[i]<0.05)
  {
    res_volcano$DiffExpressed[i] <- "Significant"
  }
  else if(res_volcano$npx_diff[i]>(0) && res_volcano$FDR[i]<0.05)
  {
    res_volcano$DiffExpressed[i] <- "Significant"
  }
  
}

for(i in 1:nrow(res_volcano))
{
  if(abs(res_volcano$pvalue[i])>(2.5e-10))
  {res_volcano$proteins[i] <- NA}
  
}
for(i in 1:nrow(res_volcano))
{
  if(res_volcano$DiffExpressed[i]=="Not Significant")
  {res_volcano$proteins[i] <- NA}
  
}

#--Plotting
ggplot(data=res_volcano, aes(x=npx_diff, y=-log10(pvalue), col=DiffExpressed, label=proteins)) + 
  scale_color_manual(breaks = c("Significant", "Not Significant"),
    values=c("#01016f", "grey"),name="Differential Expression")+
  theme_minimal() + geom_point() + 
  geom_text_repel() + xlim(-3, 3) + 
  xlab("NPX difference (day 0)") + 
  theme(plot.title = element_text(hjust = 0.5), panel.grid.minor = element_blank(), 
        panel.grid.major=element_blank(),
        axis.line = element_line(colour = "black"))
