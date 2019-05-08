## Human iPSC-derived microglia are genetically relevant to Alzheimer’s disease

### Introduction
Welcome! in [this tarball](_data/iMG_analysis.tar.gz) (4.8MB) you will find a mostly empty file tree containing all the scripts necessary to recreate the data generated for the paper.

```sh
[01;34miMG_analysis/[00m
├── [01;34mbam_counts_pe[00m
│   ├── command.txt
│   ├── featureCounts_pe_counts.txt
│   ├── featureCounts_pe_counts.txt.summary
│   ├── featureCounts_se_counts.txt
│   ├── featureCounts_se_counts.txt.summary
│   ├── featureCounts_unpr_counts.txt
│   └── featureCounts_unpr_counts.txt.summary
├── [01;34mbam_counts_se[00m
│   └── command.txt
├── [01;34mbam_counts_unpair[00m
│   └── command.txt
├── [01;34mldsc[00m
│   ├── 2019-03-07_generate_genesets.R
│   ├── 2019-03-08_generate_genesets_top20.R
│   ├── 2019-03-08_generate_genesets_top30.R
│   ├── 2019-03-08_generate_genesets_top40.R
│   ├── 2019-03-08_generate_genesets_top50.R
│   ├── 2019-03-08_ldsc_ranges.R
│   ├── 2019-03-12_ldsc_20kb_ranges.R
│   ├── 2019-03-12_ldsc_50kb_ranges.R
│   ├── command_Kozlova_20kb_top20.txt
│   ├── command_Kozlova_20kb_top30.txt
│   ├── command_Kozlova_20kb_top40.txt
│   ├── command_Kozlova_20kb_top50.txt
│   ├── command_Kozlova_20kb.txt
│   ├── command_Kozlova_50kb_top20.txt
│   ├── command_Kozlova_50kb_top30.txt
│   ├── command_Kozlova_50kb_top40.txt
│   ├── command_Kozlova_50kb_top50.txt
│   ├── command_Kozlova_50kb.txt
│   ├── command_Kozlova_top20.txt
│   ├── command_Kozlova_top30.txt
│   ├── command_Kozlova_top40.txt
│   ├── command_Kozlova_top50.txt
│   ├── command_Kozlova.txt
│   ├── Kozlova_20kb.ldcts
│   ├── Kozlova_20kb_top20.ldcts
│   ├── Kozlova_20kb_top30.ldcts
│   ├── Kozlova_20kb_top40.ldcts
│   ├── Kozlova_20kb_top50.ldcts
│   ├── Kozlova_50kb.ldcts
│   ├── Kozlova_50kb_top20.ldcts
│   ├── Kozlova_50kb_top30.ldcts
│   ├── Kozlova_50kb_top40.ldcts
│   ├── Kozlova_50kb_top50.ldcts
│   ├── Kozlova.ldcts
│   ├── Kozlova_top20.ldcts
│   ├── Kozlova_top30.ldcts
│   ├── Kozlova_top40.ldcts
│   ├── Kozlova_top50.ldcts
│   └── README
├── [01;34mMagma_analysis[00m
│   ├── 2019-02-06_sampleTable.txt
│   ├── 2019-03-03_Specificity_merged.R
│   ├── 2019-03-04_Magma_aMGL_merged_Jansen.R
│   ├── 2019-03-04_Magma_fMGL_merged_Jansen.R
│   ├── 2019-03-04_Magma_iMGL-Abud_merged_Jansen.R
│   ├── 2019-03-04_Magma_iMGL-Brownjohn_merged_Jansen.R
│   ├── 2019-03-04_Magma_iMGL-Kozlova_merged_Jansen.R
│   ├── 2019-03-04_Magma_iPMP-Kozlova_merged_Jansen.R
│   ├── 2019-03-04_Magma_scMGL_merged_Jansen.R
│   ├── 2019-03-05_Magma50-7_aMGL_merged_Jansen.R
│   ├── 2019-03-05_Magma50-7_fMGL_merged_Jansen.R
│   ├── 2019-03-05_Magma50-7_iMGL-Abud_merged_Jansen.R
│   ├── 2019-03-05_Magma50-7_iMGL-Brownjohn_merged_Jansen.R
│   ├── 2019-03-05_Magma50-7_iMGL-Kozlova_merged_Jansen.R
│   ├── 2019-03-05_Magma50-7_iPMP-Kozlova_merged_Jansen.R
│   ├── 2019-03-05_Magma50-7_scMGL_merged_Jansen.R
│   ├── 2019-03-12_Magma_aMGL_merged_BMI.R
│   ├── 2019-03-12_Magma_aMGL_merged_Nalls_PD.R
│   ├── 2019-03-12_Magma_aMGL_merged_PGC_SCZ.R
│   ├── 2019-03-12_Magma_fMGL_merged_BMI.R
│   ├── 2019-03-12_Magma_fMGL_merged_Nalls_PD.R
│   ├── 2019-03-12_Magma_fMGL_merged_PGC_SCZ.R
│   ├── 2019-03-12_Magma_iMGL-Abud_merged_BMI.R
│   ├── 2019-03-12_Magma_iMGL-Abud_merged_Nalls_PD.R
│   ├── 2019-03-12_Magma_iMGL-Abud_merged_PGC_SCZ.R
│   ├── 2019-03-12_Magma_iMGL-Brownjohn_merged_BMI.R
│   ├── 2019-03-12_Magma_iMGL-Brownjohn_merged_Nalls_PD.R
│   ├── 2019-03-12_Magma_iMGL-Brownjohn_merged_PGC_SCZ.R
│   ├── 2019-03-12_Magma_iMGL-Kozlova_merged_BMI.R
│   ├── 2019-03-12_Magma_iMGL-Kozlova_merged_Nalls_PD.R
│   ├── 2019-03-12_Magma_iMGL-Kozlova_merged_PGC_SCZ.R
│   ├── 2019-03-12_Magma_iPMP-Kozlova_merged_BMI.R
│   ├── 2019-03-12_Magma_iPMP-Kozlova_merged_Nalls_PD.R
│   ├── 2019-03-12_Magma_iPMP-Kozlova_merged_PGC_SCZ.R
│   ├── 2019-03-12_Magma_scMGL_merged_BMI.R
│   ├── 2019-03-12_Magma_scMGL_merged_Nalls_PD.R
│   ├── 2019-03-12_Magma_scMGL_merged_PGC_SCZ.R
│   ├── 2019-03-14_merged_pVals_all.R
│   ├── command.txt
│   └── README
├── [01;34mR_analysis[00m
│   ├── 2018-12-14_SraRunTable.xlsx
│   ├── CQN_expression_analysis.R
│   └── new_names.txt
├── [01;34mR_analysis_redux[00m
│   ├── 2019-01-03_SraRunTable.xlsx
│   ├── 2019-01-08_sampleTable.txt
│   ├── 2019-03-14_merged_specificity_ctdall.R
│   ├── 2019-04-06_enrichment_heatmap.R
│   ├── 2019-04-18_CQN_deseq_pe_se.R
│   ├── command.txt
│   ├── Haenseler_genes.txt
│   ├── MGL_enrich_Bennett.txt
│   └── new_names.txt
├── [01;34mread_files[00m
│   ├── command.txt
│   └── read_counts_fastq.sh
├── [01;34mreferences[00m
│   └── command.txt
├── run_pipeline.sh
├── [01;34mtrimmed_pe[00m
│   └── command.txt
└── [01;34mtrimmed_se[00m
    └── command.txt
```


### Read Processing
Given the sra data for each of the samples listed in Supplementary Table S2 are downloaded and placed in the read_files folder, the processing of the raw reads can be completed using the accompanying bash scripts, run from the master bash script `run_pipeline.sh` This will run through the basic workflow of read trimming, STAR mapping and featureCounts for the paired-end and single-end reads generating fragment counts. For examples of each step:

Trimmomatic:
```sh
```

STAR + featureCounts:
```sh
```

### General Transcriptomic Processing
Initial general transcriptomic analysis was done using the script in R_analysis, but the more complete analysis was conducted in R_analysis_redux. GWAS summary datasets will need to be downloaded into this folder before LDSC and MAGMA can be run(see command.txt). After downloading, the `2019-04-18_CQN_deseq_pe_se.R` script will do the bulk of the processing. The other two Rscripts in this folder require specificity and enrichment data from LDSC/MAGMA, so return to them after running those programs.

### MAGMA Celltyping
MAGMA should be run first, as it also generates specificity values. The scripts for various settings, cell types and GWAS traits. As described in the paper, each microglia or microglia-like cell type is run against the other non-micgroglia-like cell types. So for instance aMGL would be run against Ast, DC, Endo, Fib, MC, Neu, NSC(NPC), and Olig. A full tutorial of this protocol is available at [Nathan Skene's Guide](https://github.com/NathanSkene/MAGMA_Celltyping), but here are the basic command sets (see also the README in the folder).

Specificity (EWCE):
```R
```

Regression using MAGMA_Celltyping:
```R
```

### LDSC
LDSC also has excellent [documentation](https://github.com/bulik/ldsc/wiki/Cell-type-specific-analyses), and the instructions in the README in this folder should be sufficient. Here are the example options used for each command:

Partition heritability:
```sh
```

Munge_stats:
```sh
```

Regression:
```sh
```



[Home](README.md)