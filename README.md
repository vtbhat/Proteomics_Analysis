# Proteomics_Analysis
Analysis scripts for different types of proteomics platforms, including mass spectrometry and OLINK

### Project 1: Analysis of Long COVID vs Recovered Samples at Multiple Timepoints
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

### Project 2: Longitudinal Analysis of Severe COVID-19
