# STYK1 Knockdown RNA-Seq Analysis in Pancreatic Cancer

Bulk RNA-seq differential expression analysis to uncover biochemical regulatory changes between STYK1 knockdown and control samples in pancreatic cancer.

---

## Table of Contents

- [Background](#background)
- [Data](#data)
- [Pipeline Overview](#pipeline-overview)
  - [1. Data Retrieval](#1-data-retrieval)
  - [2. Quality Control](#2-quality-control)
  - [3. Adapter Trimming](#3-adapter-trimming)
  - [4. Transcript Quantification](#4-transcript-quantification)
  - [5. Differential Expression Analysis](#5-differential-expression-analysis)
- [Results](#results)
- [References](#references)
- [Session Info](#session-info)

---

## Background

STYK1 (Serine/Threonine/Tyrosine Kinase 1), also known as NOK (Novel Oncogene with Kinase-domain), is a receptor tyrosine kinase overexpressed across multiple cancer types, including non-small cell lung cancer, breast cancer, and pancreatic cancer. Recent work has demonstrated STYK1 as a targetable vulnerability that enhances anti-PD-L1 immunotherapy efficacy in pancreatic ductal adenocarcinoma (PDAC).

This project performs bulk RNA-seq differential expression analysis comparing STYK1 shRNA knockdown samples against scramble controls, with the goal of characterizing the downstream transcriptional consequences of STYK1 loss in pancreatic cancer cells.

---

## Data

Data were retrieved from the NCBI Sequence Read Archive (SRA), hosted on AWS. The study accession is [SRP587732](https://www.ncbi.nlm.nih.gov/sra?term=SRP587732).

| SRA Accession | Condition |
|---|---|
| SRR33705621 | shSTYK1 Knockdown |
| SRR33705622 | shSTYK1 Knockdown |
| SRR33705624 | shSTYK1 Knockdown |
| SRR33705623 | Control |
| SRR33705625 | Control |
| SRR33705626 | Control |

**Instrument:** Illumina NovaSeq X Plus  
**Library type:** Paired-end

---

## Pipeline Overview

### 1. Data Retrieval

Data are stored on AWS and require the [NCBI SRA Toolkit](https://github.com/ncbi/sra-tools/wiki) and an AWS account to access.

```bash
# Download the file from the cloud provider (AWS) to a directory
prefetch SRR33705621

# Run in the same folder where the file was prefetched
fasterq-dump [options] SRR33705621
```

### 2. Quality Control

Quality was assessed using **FastQC** and **MultiQC**. RQC was not used due to memory limitations with files of this size (~18 GB). FastQC processes files in parallel via the `-t` flag (one file per thread) and does not encounter memory constraints regardless of file size.

The average Phred score across all 6 samples (12 total FASTQ files) was **≥ 35**, exceeding the Q30 threshold (99.9% base call accuracy).

```bash
mkdir -p fastqc_out/R1 fastqc_out/R2 multiqc_out

# Run FastQC on forward reads
fastqc /path/to/sra_data/*_1.fastq \
  -t 4 \
  -o fastqc_out/R1/

# Run FastQC on reverse reads
fastqc /path/to/sra_data/*_2.fastq \
  -t 4 \
  -o fastqc_out/R2/

# Aggregate reports with MultiQC
multiqc fastqc_out/R1/ fastqc_out/R2/ \
  -o multiqc_out/ \
  --title "Combined QA for Forward and Reverse Reads (R1 & R2)"
```

### 3. Adapter Trimming

Adapter sequences — synthetic components ligated to reads during library preparation to anchor them to the flow cell surface — were removed using **fastp** prior to mapping. Adapter removal improves downstream mapping rates.

```bash
for sample in SRR33705621 SRR33705622 SRR33705623 SRR33705624 SRR33705625 SRR33705626
do
  fastp \
    -i ${sample}_1.fastq \
    -I ${sample}_2.fastq \
    -o trimmed/${sample}_1_trimmed.fastq \
    -O trimmed/${sample}_2_trimmed.fastq \
    --detect_adapter_for_pe \
    --thread 10 \
    -h trimmed/${sample}_fastp_report.html \
    -j trimmed/${sample}_fastp_report.json
done
```

### 4. Transcript Quantification

Trimmed reads were mapped to the human reference transcriptome using **Salmon** — a lightweight, alignment-free quantification tool designed for transcript-level expression estimation. A decoy-aware Salmon index was built prior to quantification.

```bash
salmon quant \
  -i ~/salmon_index \
  -l A \
  -1 SRR33705621_1_trimmed.fastq \
  -2 SRR33705621_2_trimmed.fastq \
  -p 9 \
  --validateMappings \
  -o SRR33705621_quant
```

> The output folder contains `quant.sf`, which holds the raw count data required for DESeq2 import via `tximport`.

**Mapping rates across samples:**

| Sample | Condition | Mapping Rate |
|---|---|---|
| SRR33705621 | shSTYK1 KD | 87.16% |
| SRR33705622 | shSTYK1 KD | 86.37% |
| SRR33705624 | shSTYK1 KD | 87.71% |
| SRR33705623 | Control | 87.33% |
| SRR33705625 | Control | 88.38% |
| SRR33705626 | Control | 88.30% |

### 5. Differential Expression Analysis

Downstream differential expression (DE) analysis was performed in **R** using **DESeq2**, with counts imported via `tximport`. LFC shrinkage was applied using the `apeglm` method (Zhu et al., 2018).

**Filtering note:** No more than 30% of genes were filtered prior to modeling. Over-filtering reduces DESeq2's ability to accurately fit gene dispersions across the mean, leading to incomplete dispersion estimates and decreased statistical power.

---

## Results

DE analysis was run at **α = 0.05** across 11,960 mapped genes, with an estimated ~21 expected false positives.

| Category | Gene Count |
|---|---|
| Total mapped genes | 11,960 |
| Upregulated (LFC > 0) | 211 |
| Downregulated (LFC < 0) | 214 |
| Outlier genes | 2 |
| Low-count genes filtered | 2,383 |

**PCA:** PC1 explained 89% of variance, driven by the knockdown vs. control condition. PC2 captured within-group biological variation.

**Notable differentially expressed genes:**

- **SHH** — Downregulated upon STYK1 knockdown. SHH (Sonic Hedgehog) is a developmental morphogen that, when dysregulated, drives the occurrence of multiple cancers including pancreatic cancer.
- **EFNB2** — Downregulated upon STYK1 knockdown. EFNB2 regulates endothelial cell survival and angiogenesis, suggesting STYK1 may promote tumor vascularization via this axis.

These results establish STYK1 as a valid therapeutic target and reveal the downstream transcriptional circuitry disrupted by its loss in pancreatic cancer.

**Figures:**

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

## References

1. Ma, Z., Liu, D., Li, W. et al. STYK1 promotes tumor growth and metastasis by reducing SPINT2/HAI-2 expression in non-small cell lung cancer. *Cell Death Dis* **10**, 435 (2019). https://doi.org/10.1038/s41419-019-1659-1

2. Mohammed, R., Ismaeel, A., Alshaikh, S. et al. STYK1 expression in breast cancer and its association with vascular invasion and clinicopathological features. *Sci Rep* **16**, 7775 (2026). https://doi.org/10.1038/s41598-026-39385-8

3. Wang, X., Zhang, Y., Li, A. et al. STYK1 as a targetable vulnerability enhances the efficacy of anti-PD-L1 therapy in pancreatic cancer. *J Gastroenterol* **60**, 1449–1469 (2025). https://doi.org/10.1007/s00535-025-02291-3

4. Zhu, A., Ibrahim, J.G., Love, M.I. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics* (2018). https://doi.org/10.1093/bioinformatics/bty895

---

## Acknowledgements

Raw sequencing data were provided by the authors of reference [1] and used here for educational purposes.

---

## Session Info

```
R version 4.5.3 (2026-03-11)
Platform: aarch64-apple-darwin20
Running under: macOS Tahoe 26.3.1

Key packages:
  DESeq2 1.50.2 | tximport 1.38.2 | apeglm 1.32.0
  ggplot2 4.0.3 | pheatmap 1.0.13 | biomaRt 2.66.1
  tidyverse 2.0.0 | org.Hs.eg.db 3.22.0
```

Full session info available in `session_info.txt`.
