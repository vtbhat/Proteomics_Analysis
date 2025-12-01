# Proteomics_Analysis
Analysis scripts for different types of proteomics platforms, including mass spectrometry and OLINK

### Project 1: Longitudinal Analysis of Severe COVID-19
Analysis script: OLINK_severeCOVID.R

In a paper by [Filbin et al.,](https://pmc.ncbi.nlm.nih.gov/articles/PMC8091031/#sec2), the authors investigated 1463 proteins longitudinally across 306 COVID-19 patients and 78 symptomatic controls. The authors identified multiple protein signatures associated with the severity of COVID-19. Here, I implement a linear model to reproduce the results from Figure 1 of the paper:

Heatmap of the top 200 differentially expressed proteins

<img width="765" height="547" alt="image" src="https://github.com/user-attachments/assets/73068591-9ce1-4ff0-81d5-0608b6dfbac7" />

Volcano Plot of Differentially Expressed Proteins

<img width="1037" height="591" alt="image" src="https://github.com/user-attachments/assets/2473fbe6-4cbd-4900-b3f6-cc8fb61b5ef3" />

Boxplot of Selected Differentially Expressed Proteins (Viral and Interferon Signaling Pathways)

<img width="1920" height="992" alt="image" src="https://github.com/user-attachments/assets/ae8e0b17-60dc-4e81-9536-89290d2adcd7" />

As shown in the paper, proteins associated with viral response and inflammatory signaling, such as DDX58, IFNG, CCL7, and CXCL10, were upregulated in COVID+ patients compared to symptomatic controls.

### Project 2: Analysis of Long COVID vs Recovered Samples at Multiple Timepoints
Analysis script: OLINK_analysis.R

In a paper by [Hamlin et al.,](https://pmc.ncbi.nlm.nih.gov/articles/PMC12148066/#S2), the authors analyzed sex-specific differences in long COVID development using proteomic data assayed using the OLINK Inflammation and Immune Response panels. The data was collected across three timepoints: during acute infection, 3 months after infection, and 12 months after infection. The authors observed that TGF-beta-1 levels (also called LAP TGF-beta-1) were higher in participants with long COVID after 3 months when compared to individuals who had fully recovered at the same timepoint (unadjusted p-value < 0.05). Here, I implement a limma model to reproduce the same results:

Table of all genes with p-value < 0.05. None of the genes passed the same Benjamini-Hochberg FDR threshold due to small sample size, but 11 genes had p-value < 0.05.
| Gene            | logFC   | P.Value  | adj.P.Val |
|-----------------|---------|----------|-----------|
| IL13        | -0.907  | 0.0003    | 0.061     |
| EIF5A       | -0.872  | 0.001    | 0.065     |
| KPNA1       | -0.868  | 0.002    | 0.122     |
| MCP-4       | -0.697  | 0.018    | 0.765     |
| FCRL3       | -0.403  | 0.030    | 0.765     |
| MCP-1       | -0.273  | 0.035    | 0.765     |
| NF2         | -0.756  | 0.035    | 0.765     |
| ARNT        | -0.463  | 0.036    | 0.765     |
| LIF-R       | 0.297   | 0.039    | 0.765     |
| MILR1       | 0.520   | 0.044    | 0.765     |
| **LAP TGF-beta-1** | **0.355** | **0.047** | **0.765** |

Furthermore, the authors found that the TGF-beta=1 levels were found to be significantly different during acute infection between participants who eventually recovered vs those who went on to develop long COVID (Wilcoxon p < 0.05). I've reproduced the same plot here:

<img width="1021" height="470" alt="image" src="https://github.com/user-attachments/assets/49af55fd-f424-4f82-91af-8f8860cc154a" />


