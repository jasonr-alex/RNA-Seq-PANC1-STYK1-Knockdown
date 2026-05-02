# Bulk RNA-seq Analysis: STYK1 Knockdown in Panc-1 Pancreatic Cancer Cells

## Overview

This repository contains a complete end-to-end bulk RNA-seq analysis pipeline investigating the transcriptional effects of **STYK1 (Serine Threonine Tyrosine Kinase 1) knockdown** in Panc-1 human pancreatic cancer cells. STYK1 is a receptor tyrosine kinase implicated in cancer cell proliferation and survival, making it a compelling target in pancreatic cancer research — a disease with one of the lowest 5-year survival rates in oncology (~12%).

The analysis was performed on publicly available paired-end sequencing data comprising **6 samples** (3 control, 3 knockdown) downloaded from the NCBI Sequence Read Archive (SRA).

---

## Biological Question

What are the transcriptional consequences of STYK1 knockdown in Panc-1 pancreatic cancer cells, and what biological pathways are most affected?

---

## Table of Contents

- [Pipeline Overview](#pipeline-overview)
- [Tools & Dependencies](#tools--dependencies)
- [Data](#data)
- [Methods](#methods)
- [Results](#results)
- [Repository Structure](#repository-structure)
- [How to Run](#how-to-run)
- [References](#references)
- [Author](#author)
- [Acknowledgements](#acknowledgements)

---

## Pipeline Overview

```
Raw FASTQ files
      ↓
Quality Control (FastQC + MultiQC)
      ↓
Adapter Trimming (fastp)
      ↓
Transcript Quantification (Salmon)
      ↓
Import & Gene-level Summarization (tximport)
      ↓
Differential Expression Analysis (DESeq2)
      ↓
Visualization & Annotation (ggplot2, pheatmap, org.Hs.eg.db)
      ↓
Annotated DE Gene List (CSV)
```

---

## Tools & Dependencies

### Command Line
| Tool | Version | Purpose |
|---|---|---|
| SRA Toolkit (fasterq-dump) | Latest | Raw data download |
| FastQC | 0.12.1 | Quality control |
| MultiQC | 1.33 | QC report aggregation |
| fastp | Latest | Adapter trimming |
| Salmon | 1.11.4 | Transcript quantification |

### R Packages
| Package | Purpose |
|---|---|
| DESeq2 | Differential expression analysis |
| tximport | Import Salmon output |
| biomaRt | tx2gene mapping |
| org.Hs.eg.db | Gene annotation |
| ggplot2 | Visualization |
| pheatmap | Heatmap generation |
| apeglm | LFC shrinkage |
| RColorBrewer | Color palettes |
| tidyverse | Data manipulation |

---

## Data

Data were retrieved from the NCBI Sequence Read Archive (SRA), hosted on AWS. The study accession is [SRP587732](https://www.ncbi.nlm.nih.gov/sra?term=SRP587732).

- **Organism:** Homo sapiens
- **Cell line:** Panc-1 (human pancreatic cancer)
- **Condition:** STYK1 shRNA knockdown vs. scramble control
- **Sequencing:** Paired-end, 151bp reads
- **Instrument:** Illumina NovaSeq X Plus
- **Samples:** 6 total (3 control, 3 knockdown)

| SRA Accession | Condition | Mapping Rate |
|---|---|---|
| SRR33705621 | shSTYK1 Knockdown | 87.16% |
| SRR33705622 | shSTYK1 Knockdown | 86.37% |
| SRR33705623 | Control | 87.33% |
| SRR33705624 | shSTYK1 Knockdown | 87.71% |
| SRR33705625 | Control | 88.38% |
| SRR33705626 | Control | 88.30% |

**Average mapping rate: 87.54%** across all samples.

---

## Methods

### 1. Quality Control
FastQC was run on all 12 FASTQ files (6 samples × R1 + R2). MultiQC aggregated the results into a single report. All samples showed high per-base sequence quality (Phred > 35 across all 151bp positions) with minimal adapter contamination. RQC was not used due to memory limitations with files of this size.

### 2. Adapter Trimming
fastp was used to trim adapter sequences with `--detect_adapter_for_pe` for automatic paired-end adapter detection.

### 3. Quantification
Salmon (v1.11.4) was used for quasi-mapping quantification against the human transcriptome (Ensembl GRCh38 release 111). The `--validateMappings` flag was applied for improved accuracy.

### 4. Differential Expression Analysis
Transcript-level counts were imported and summarized to gene level using `tximport` with a tx2gene mapping table generated via `biomaRt`. A DESeq2 object was constructed without pre-filtering to retain maximum statistical power. Normalization was performed using `rlog()` transformation (`blind = TRUE`) for visualization. Differential expression was tested using the Wald test with a Benjamini-Hochberg adjusted p-value threshold of **α = 0.05**. Log2 fold changes were shrunk using the `apeglm` method (Zhu et al., 2018) to reduce noise from low-count genes.

---

## Results

### Quality Metrics
- All 12 files passed per-base sequence quality checks
- Uniform read length of 151bp across all samples
- Overrepresented sequences < 1% across all samples
- Adapter content detected at 3' ends — corrected with fastp trimming

### Sample-Level QC
- PCA plot confirmed clean separation between control and knockdown groups along PC1 (75% variance)
- Sample correlation heatmap showed high within-group similarity (r > 0.999) and clear between-group differences

### Differential Expression
At FDR < 5%, **2,808 genes** were significantly differentially expressed across 22,248 genes with nonzero total read count:

| Category | Gene Count |
|---|---|
| Total genes with nonzero counts | 22,248 |
| Upregulated in knockdown (LFC > 0) | 1,382 (6.2%) |
| Downregulated in knockdown (LFC < 0) | 1,426 (6.4%) |
| Outlier genes | 28 (0.13%) |
| Low-count genes filtered | 6,223 (28%) |

With α = 0.05 applied across 22,248 genes, an estimated **~140 genes** are expected to be false positives, meaning the vast majority of the 2,808 significant genes represent real biological signals.

### Key Genes of Biological Interest

**ALCAM (Activated Leukocyte Cell Adhesion Molecule)**
Downregulated following STYK1 knockdown. ALCAM overexpression is directly linked to cancer cell migration and metastatic potential. Its downregulation suggests STYK1 knockdown may suppress the invasive capacity of Panc-1 cells — a therapeutically relevant finding in pancreatic cancer, where metastasis is a primary driver of mortality.

**HLA-A (Human Leukocyte Antigen A)**
Downregulated following STYK1 knockdown. HLA-A is a MHC class I molecule responsible for presenting tumor antigens to cytotoxic T cells. Reduced HLA-A expression is associated with decreased response to PD-L1 checkpoint inhibitor therapy, as T cells require antigen presentation to recognize and destroy tumor cells. This finding has direct implications for understanding how STYK1 knockdown may influence immunotherapy efficacy in pancreatic cancer — specifically suggesting that STYK1 loss may reduce the enhanced efficacy observed with anti-PD-L1 therapy by impairing antigen presentation.

**CDH2 (N-cadherin)**
Downregulated in knockdown — involved in epithelial-mesenchymal transition (EMT), a key process in cancer invasion and metastasis.

**INHBB & BMP5**
Differentially expressed genes belonging to the TGF-β signaling pathway, a major driver of pancreatic cancer progression.

**CLDN1 (Claudin-1)**
Downregulated — tight junction protein frequently dysregulated in pancreatic cancer, involved in tumor invasion.

**CYP1B1**
Differentially expressed — a drug metabolism enzyme with known implications in chemotherapy resistance in pancreatic cancer.

**MXRA5**
Downregulated — matrix remodeling-associated protein involved in extracellular matrix and tumor microenvironment remodeling.

These findings collectively suggest STYK1 plays a broad role in regulating **cell adhesion, epithelial-mesenchymal transition, immune recognition, and drug metabolism** pathways in Panc-1 pancreatic cancer cells — with direct implications for drug target prioritization, metastasis suppression, and immunotherapy response.

### Figures

| Figure | Description |
|---|---|
| Fig 1 | Heatmap for Exploratory Data Analysis (EDA) of Samples |
| Fig 2 | PCA — Variance Between PC1 and PC2 |
| Fig 3 | Dispersion Graph of DESeq2 Fitted Line |
| Fig 4 | MA Plot (unshrunken) |
| Fig 5 | MA Plot with LFC Shrinkage (apeglm) |
| Fig 6 | Heatmap with Hierarchical Clustering of Significant Genes |
| Fig 7 | Volcano Plot for DE Analysis |
| Fig 8 | Top 20 Genes by Significance Level |

---

## Repository Structure

```
├── RNA_seq_analysis.Rmd                              # Full annotated R analysis pipeline
├── DE_sig_results_annotated_STYK1_knockdown.csv      # Significant DE gene list with annotations
├── plots/
│   ├── correlation_heatmap.pdf                       # Fig 1
│   ├── PCA_plot.pdf                                  # Fig 2
│   ├── Dispersion_Graph_of_DESeq2_Fitted_Line.pdf    # Fig 3
│   ├── MA_plot_before_shrinkage.pdf                  # Fig 4
│   ├── MA_plot_after_shrinkage.pdf                   # Fig 5
│   ├── Heatmap_Hierarchical_Clustering.pdf           # Fig 6
│   ├── Volcano_Plot_for_DE_Analysis.pdf              # Fig 7
│   └── Top_20_Genes_by_Significance_Level.pdf        # Fig 8
└── README.md
```

---

## How to Run

### Prerequisites
Install required R packages:
```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "tximport", "org.Hs.eg.db", "biomaRt", "apeglm"))
install.packages(c("tidyverse", "pheatmap", "RColorBrewer"))
```

Install command-line tools via conda:
```bash
conda install -c bioconda fastqc multiqc fastp salmon sra-tools
```

### Steps
1. Download raw data using `fasterq-dump` with SRR accession numbers listed above
2. Run FastQC and MultiQC for quality control
3. Trim adapters with fastp
4. Build Salmon index from human transcriptome (Ensembl GRCh38)
5. Run Salmon quantification for each sample
6. Open and run `RNA_seq_analysis.Rmd` in RStudio

---

## References

1. Ma, Z., Liu, D., Li, W. et al. STYK1 promotes tumor growth and metastasis by reducing SPINT2/HAI-2 expression in non-small cell lung cancer. *Cell Death Dis* **10**, 435 (2019). https://doi.org/10.1038/s41419-019-1659-1

2. Mohammed, R., Ismaeel, A., Alshaikh, S. et al. STYK1 expression in breast cancer and its association with vascular invasion and clinicopathological features. *Sci Rep* **16**, 7775 (2026). https://doi.org/10.1038/s41598-026-39385-8

3. Wang, X., Zhang, Y., Li, A. et al. STYK1 as a targetable vulnerability enhances the efficacy of anti-PD-L1 therapy in pancreatic cancer. *J Gastroenterol* **60**, 1449–1469 (2025). https://doi.org/10.1007/s00535-025-02291-3

4. Zhu, A., Ibrahim, J.G., Love, M.I. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics* (2018). https://doi.org/10.1093/bioinformatics/bty895

---

## Acknowledgements

Raw sequencing data were obtained from the NCBI Sequence Read Archive (SRA) and used here for educational and research purposes. Analysis was performed using open-source bioinformatics tools. Gene annotations sourced from Ensembl GRCh38 and the org.Hs.eg.db Bioconductor package.

---

## Session Info

```
R version 4.5.3 (2026-03-11)
Platform: aarch64-apple-darwin20
Running under: macOS Tahoe 26.3.1

Key attached packages:
  DESeq2 1.50.2           tximport 1.38.2         apeglm 1.32.0
  biomaRt 2.66.1          org.Hs.eg.db 3.22.0     AnnotationDbi 1.72.0
  ggplot2 4.0.3           pheatmap 1.0.13          RColorBrewer 1.1-3
  tidyverse 2.0.0         dplyr 1.2.1              tibble 3.3.1
  GenomicRanges 1.62.1    SummarizedExperiment 1.40.0  BiocManager 1.30.27
```

Full session info available in `session_info.txt`.
