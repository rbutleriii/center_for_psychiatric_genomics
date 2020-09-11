## Sex-specific transcriptomic profiles and imprinting drive nicotine sensitization and self-administration in rats and inform the genetic basis of human smoking behavior

[Home](README.md)
### Introduction
Welcome! in <a href="https://github.com/rbutleriii/center_for_psychiatric_genomics/raw/master/_data/nic_analysis.tar.gz" download="download">this tarball</a> you will find a mostly empty file tree containing all the scripts necessary to recreate the data generated for the paper.

```bash
nic_analysis
├── Duan_Project_010
├── Liu_etal_genes.R
├── Liu_et_al_genes.txt
├── refs
│   └── command.txt
├── run_pipeline.sh
├── self_admin
│   ├── 20190815_sampleTable.txt
│   ├── allelic_imbalance
│   │   ├── ABgroups.R
│   │   ├── AB_samples.txt
│   │   ├── all_imb_calling.sh
│   │   ├── david
│   │   │   ├── 20200414_AI_human_david_plots.R
│   │   │   └── 20200414_david_table.sh
│   │   ├── masked_gvcfs
│   │   │   └── AI_proptest_plots.R
│   │   ├── merge_list.txt
│   │   ├── prop_test
│   │   │   ├── 20200808_AI_genesets.R
│   │   │   ├── 20200808_AI_human_genes.R
│   │   │   ├── 20200808_AI_self_admin_gene_overlaps.R
│   │   │   ├── 20200808_geneloc_generator.R
│   │   │   ├── AI_differentials.R
│   │   │   └── command.txt
│   │   └── samples_list.txt
│   ├── david_split_DE
│   │   ├── 20200416_david_plots.R
│   │   └── 20200416_david_table.sh
│   ├── david_split_downgenes
│   │   ├── 20200514_david_plots.R
│   │   └── 20200514_david_table.sh
│   ├── david_split_upgenes
│   │   ├── 20200514_david_plots.R
│   │   └── 20200514_david_table.sh
│   ├── dsq_broad_analysis
│   │   └── 20191004_Deseq_broad_analysis.R
│   ├── dsq_split_dropvta_DE
│   │   ├── 20200129_dsq_split_DE.R
│   │   ├── 20200513_dsq_human_genes.R
│   │   ├── 20200513_dsq_split_genes.R
│   │   ├── 20200515_SA_enrichment_sets.R
│   │   └── 20200527_SA_dsq_split_annotations.R
│   ├── dsq_split_dropvta_downgenes
│   │   ├── 20200513_dsq_human_downgenes.R
│   │   └── 20200513_dsq_split_downgenes.R
│   ├── dsq_split_dropvta_upgenes
│   │   ├── 20200513_dsq_human_upgenes.R
│   │   └── 20200513_dsq_split_upgenes.R
│   ├── magma_allelic_imbalance
│   │   ├── command.txt
│   │   ├── magma_geneset.py
│   │   ├── magma_posthoc_QC.R
│   │   ├── magma_ranges.R
│   │   └── posthoc_qc_107a.r
│   ├── magma_split_DE
│   │   ├── command.txt
│   │   ├── magma_geneset.py
│   │   ├── magma_posthoc_QC.R
│   │   ├── magma_ranges.R
│   │   └── posthoc_qc_107a.r
│   ├── magma_split_downgenes
│   │   ├── command.txt
│   │   ├── magma_geneset.py
│   │   ├── magma_posthoc_QC.R
│   │   ├── magma_ranges.R
│   │   └── posthoc_qc_107a.r
│   └── magma_split_upgenes
│       ├── command.txt
│       ├── magma_geneset.py
│       ├── magma_posthoc_QC.R
│       ├── magma_ranges.R
│       └── posthoc_qc_107a.r
├── sensitization
│   ├── 20190706_sampleTable.txt
│   ├── david_split_DE
│   │   ├── 20200416_david_plots.R
│   │   └── 20200416_david_table.sh
│   ├── david_split_downgenes
│   │   ├── 20200514_david_plots.R
│   │   └── 20200514_david_table.sh
│   ├── david_split_upgenes
│   │   ├── 20200514_david_plots.R
│   │   └── 20200514_david_table.sh
│   ├── dsq_broad_analysis
│   │   └── 20190930_Deseq_broad_analysis.R
│   ├── dsq_split_DE_analysis
│   │   ├── 20191211_dsq_split_DE.R
│   │   ├── 20200511_dsq_split_genes.R
│   │   ├── 20200511_dsq_split_human_genes.R
│   │   ├── 20200515_SST_enrichment_sets.R
│   │   └── 20200527_dsq_split_annotations.R
│   ├── dsq_split_downgenes
│   │   ├── 20200513_dsq_human_downgenes.R
│   │   └── 20200513_dsq_split_downgenes.R
│   ├── dsq_split_upgenes
│   │   ├── 20200513_dsq_human_upgenes.R
│   │   └── 20200513_dsq_split_upgenes.R
│   ├── magma_split_DE
│   │   ├── command.txt
│   │   ├── magma_geneset.py
│   │   ├── magma_posthoc_QC.R
│   │   ├── magma_ranges.R
│   │   └── posthoc_qc_107a.r
│   ├── magma_split_downgenes
│   │   ├── command.txt
│   │   ├── magma_geneset.py
│   │   ├── magma_posthoc_QC.R
│   │   ├── magma_ranges.R
│   │   └── posthoc_qc_107a.r
│   ├── magma_split_interaction
│   │   ├── command.txt
│   │   ├── magma_geneset.py
│   │   ├── magma_interact_ranges.R
│   │   ├── magma_posthoc_QC.R
│   │   └── posthoc_qc_107a.r
│   └── magma_split_upgenes
│       ├── command.txt
│       ├── magma_geneset.py
│       ├── magma_posthoc_QC.R
│       ├── magma_ranges.R
│       └── posthoc_qc_107a.r
├── SST-SA_comparison
│   ├── correlations
│   │   └── 20200204_dsq_gene_plots.R
│   ├── magma_interaction_males
│   │   ├── command.txt
│   │   ├── magma_geneset.py
│   │   ├── magma_interact_ranges.R
│   │   ├── magma_posthoc_QC.R
│   │   └── posthoc_qc_107a.r
│   ├── males_overlaps
│   │   └── 20200508_SST-SA_male_overlaps_sets.R
│   ├── overlaps
│   │   └── 20200203_upset_ultrawide.R
│   └── SST_overlaps
│       └── 20200511_SST_overlap_sets.R
├── Walker_AdIn_genes.txt
├── Walker_SA_genes.R
└── Walker_SA_genes.txt
```


### Read Processing
Given the sra data for each of the samples listed in GSE157726 and GSE157683 are downloaded and placed in the Duan_Project_010 folder, the processing of the raw reads can be completed using the accompanying bash scripts, run from the master bash script `run_pipeline.sh` This will run through the basic workflow of read trimming, salmon & STAR mapping and featureCounts for the paired-end reads. For this paper, salmon->tximport->deseq2 was used for DE, and STAR for allelic imbalance mapping. For examples of each step in the initial mapping/counting:

Fastp:
```bash
fastp \
  --thread 2\
  -i $1 \
  -I $R2 \
  -o ${out}_R1_trim.fastq.gz \
  -O ${out}_R2_trim.fastq.gz \
  --adapter_sequence AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
  --adapter_sequence_r2 GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG \
  -3 \
  -5 \
  -l 20 \
  -h ${out}.html \
  -j ${out}.json
```

Salmon:
```bash
salmon --no-version-check quant \
  --seqBias \
  --gcBias \
  --validateMappings \
  -i $2 \
  -l A \
  -p 20 \
  -1 <(pigz -dc "$1") \
  -2 <(pigz -dc "$R2") \
  -o ${out}
```

STAR:
```bash
STAR \
  --runThreadN 20 \
  --genomeDir $2 \
  --genomeLoad LoadAndKeep \
  --readFilesCommand zcat \
  --outFileNamePrefix ${out}_star \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMmapqUnique 60 \
  --outSAMattrRGline ID:$out LB:lib1 PL:illumina PU:unit1 SM:$out \
  --outBAMsortingThreadN 10 \
  --outBAMsortingBinsN 10 \
  --limitBAMsortRAM 10000000000 \
  --outBAMcompression 10 \
  --bamRemoveDuplicatesType UniqueIdentical \
  --outFilterIntronMotifs RemoveNoncanonical \
  --outFilterMismatchNmax 2 \
  --outFilterScoreMinOverLread 0.30 \
  --outFilterMatchNminOverLread 0.30 \
  --alignSoftClipAtReferenceEnds No\
  --readFilesIn $1 $R2
```

### General Transcriptomic Processing
Initial general transcriptomic analysis was done separately for SST and SA datasets using the scripts their respective folders. Individual gene sets need to be run in the DAVID and STRING web interfaces manually. 

### MAGMA Enrichment of disease risk
GWAS summary datasets will need to be downloaded into this folder before MAGMA can be run(see <a href="https://github.com/rbutleriii/center_for_psychiatric_genomics/raw/master/_data/nic_gwas_refs.tar.gz" download="download">this tarball</a>). The scripts for various settings and gene sets are in their associated folders. The custom script magma_geneset.py requires specification of gwas sets be hardcoded in for now, and multithreads by number of gwas sets (which can be then run in parallel). 

Example magma with multiple snp annotation windows (and 6 GWAS traits; 9x6 = 54 threads):
```bash
parallel -j9 --delay 2 --link python magma_geneset.py \
  --grch38 \
  -u {1} -d {2} \
  -b ai-UW \
  -c hla protein_coding expressed \
  -o magma_geneset_{1}_{2} \
  ::: '0' '10' '20' '50' '100' '10' '20' '50' '100'\
  ::: '0' '1.5' '3' '7.5' '15' '10' '20' '50' '100'
Rscript magma_ranges.R
Rscript magma_posthoc_QC.R
```

### In summary
So, there are the general steps to recreate the raw data used in our publication. As a disclaimer, this information is provided as is and likely requires modification to run correctly on other machines. Understand the commands before running them. I am not responsible if you try to run 100 threads on a laptop. 

[Home](README.md)

