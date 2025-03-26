*ATAC-seq Analysis in SMARCB1-Altered Bladder Cancer Cells*

##Overview

This repository contains the analysis workflow and results for studying chromatin accessibility changes in bladder cancer cells with SMARCB1 alterations using ATAC-seq data. The study aims to understand the impact of SMARCB1 on chromatin dynamics and its potential implications for cancer research and treatment.

##Introduction

Bladder cancer is the second most common urological cancer in the U.S., with advanced stages being difficult to treat. The SWI/SNF chromatin remodeling complex, frequently dysregulated in bladder cancer, plays a key role in epigenetic regulation. SMARCB1, a core component of this complex, is linked to aggressive tumor behavior. This study leverages ATAC-seq to analyze chromatin accessibility changes resulting from SMARCB1 knockout, rescue, and control conditions, potentially unveiling new therapeutic targets.

##Dataset

The dataset used in this study is publicly available from the Gene Expression Omnibus (GEO) under accession number GSE213964. The raw sequencing files were obtained from the Sequence Read Archive (SRA):

SMARCB1 Rescue Samples: SRR21677760, SRR21677761

SMARCB1 Knockout Samples: SRR21677762, SRR21677763

Control Samples: SRR21677764, SRR21677765

##Methods

###Workflow

1. Raw Data Acquisition: Downloaded using SRA Toolkit.

2. Quality Control: Assessed using FastQC.

3. Alignment & Indexing: Performed with Bowtie2 using the human genome reference (hg38).

4. Peak Calling: Conducted using MACS2 to identify open chromatin regions.

5. Differential Accessibility Analysis: Performed using DiffBind to compare chromatin accessibility across conditions.

##Results

1. Quality Control: High-quality data confirmed, no trimming required.

2. Alignment: Achieved a 97% alignment rate with over 88 million reads.

3. Peak Calling: Identified differential peaks across conditions.

4. Differential Accessibility Analysis:

         Control vs. Knockout: 80,347 differentially accessible peaks.

         Control vs. Rescue: 9,250 differentially accessible peaks.

         Notable differences in chromatin accessibility indicate that SMARCB1 deletion significantly alters chromatin structure.

##Conclusion

Our analysis highlights the critical role of SMARCB1 in modulating chromatin accessibility in bladder cancer. The results suggest that SMARCB1 loss may drive significant changes in chromatin structure, which could have implications for tumor progression and therapeutic strategies. Further research is needed to identify the genomic regions and genes affected by SMARCB1 alterations.

##Dependencies

    FastQC (for quality control)

    SRA Toolkit (for downloading raw data)

    Bowtie2 (for alignment)
 
    Samtools (for indexing)

    MACS2 (for peak calling)

    DiffBind (R package) (for differential accessibility analysis)
