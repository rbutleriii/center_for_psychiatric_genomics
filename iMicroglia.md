## Human iPSC-derived microglia are genetically relevant to Alzheimer’s disease

[Home](README.md)
### Introduction
Welcome! in <a href="https://github.com/rbutleriii/center_for_psychiatric_genomics/blob/master/_data/iMG_analysis.tar.gz" download="download">this tarball</a> you will find a mostly empty file tree containing all the scripts necessary to recreate the data generated for the paper.

```bash
iMG_analysis/
├── bam_counts_pe
│   └── command.txt
├── bam_counts_se
│   └── command.txt
├── bam_counts_unpair
│   └── command.txt
├── ldsc
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
├── Magma_analysis
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
├── R_analysis
│   ├── 2018-12-14_SraRunTable.xlsx
│   ├── CQN_expression_analysis.R
│   └── new_names.txt
├── R_analysis_redux
│   ├── 2019-01-03_SraRunTable.xlsx
│   ├── 2019-01-08_sampleTable.txt
│   ├── 2019-03-14_merged_specificity_ctdall.R
│   ├── 2019-04-06_enrichment_heatmap.R
│   ├── 2019-04-18_CQN_deseq_pe_se.R
│   ├── command.txt
│   ├── Haenseler_genes.txt
│   ├── MGL_enrich_Bennett.txt
│   └── new_names.txt
├── read_files
│   ├── command.txt
│   └── read_counts_fastq.sh
├── references
│   └── command.txt
├── run_pipeline.sh
├── trimmed_pe
│   └── command.txt
└── trimmed_se
    └── command.txt
```


### Read Processing
Given the sra data for each of the samples listed in Supplementary Table S2 are downloaded and placed in the read_files folder, the processing of the raw reads can be completed using the accompanying bash scripts, run from the master bash script `run_pipeline.sh` This will run through the basic workflow of read trimming, STAR mapping and featureCounts for the paired-end and single-end reads generating fragment counts. For examples of each step:

Trimmomatic:
```bash
java -jar $trimmopath PE \
    -threads 12 \
    -basein "$1" \
    -baseout "${out}_trim.fastq.gz" \
    ILLUMINACLIP:${adapterpath}:2:30:10 \
    LEADING:30 \
    TRAILING:30 \
    SLIDINGWINDOW:4:15 \
    MINLEN:20
```

STAR + featureCounts:
```bash
STAR \
    --runThreadN 20 \
    --genomeDir $DB \
    --genomeLoad LoadAndKeep \
    --readFilesCommand zcat \
    --outFileNamePrefix ${out}_star \
    --outSAMtype BAM SortedByCoordinate \
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
    --readFilesIn $1 $2

featureCounts \
    -T 64 \
    -p \
    -B \
    -C \
    -D 2000 \
    -t gene \
    -a $GTF \
    -o ${outfile}_pe_counts.txt *.bam
```

### General Transcriptomic Processing
Initial general transcriptomic analysis was done using the script in R_analysis, but the more complete analysis was conducted in R_analysis_redux. GWAS summary datasets will need to be downloaded into this folder before LDSC and MAGMA can be run(see command.txt). After downloading, the `2019-04-18_CQN_deseq_pe_se.R` script will do the bulk of the processing. The other two Rscripts in this folder require specificity and enrichment data from LDSC/MAGMA, so return to them after running those programs.

### MAGMA Celltyping
MAGMA should be run first, as it also generates specificity values. The scripts for various settings, cell types and GWAS traits. As described in the paper, each microglia or microglia-like cell type is run against the other non-micgroglia-like cell types. So for instance aMGL would be run against Ast, DC, Endo, Fib, MC, Neu, NSC(NPC), and Olig. A full tutorial of the MAGMA_Celltyping protocol is available at [Nathan Skene's Github](https://github.com/NathanSkene/MAGMA_Celltyping), but here are the basic command sets (see also the README in the folder).

Specificity (EWCE):
```R
# ctd object
input_ctd <- list(exp=exp, annot=annot)

# drop uninformative genes
exp_dropped <- drop.uninformative.genes(exp=input_ctd$exp, level2annot=input_ctd$annot$level1class)

# generate specificities
annotLevels <- list(level1class=input_ctd$annot$level1class)
fNames_out = generate.celltype.data(exp=exp_dropped,annotLevels=annotLevels,
                                  groupName=paste0(MGL_name, "_merged"))
```

Regression using MAGMA_Celltyping:
```R
# Map SNPs to Genes
genesOutPath = map.snps.to.genes(gwas_sumstats_path, genome_ref_path=genome_ref_path)

# linear enrichment
ctAssocsLinear = calculate_celltype_associations(ctd, gwas_sumstats_path, genome_ref_path=genome_ref_path,
                                                 specificity_species="human", analysis_name=mgl_type)
FigsLinear = plot_celltype_associations(ctAssocsLinear, ctd=ctd, fileTag=mgl_type)

# top 10%
ctAssocsTop = calculate_celltype_associations(ctd, gwas_sumstats_path, genome_ref_path=genome_ref_path,
                                              specificity_species="human", EnrichmentMode="Top 10%",
                                              analysis_name=paste0(mgl_type, "top10"))
FigsTopDecile = plot_celltype_associations(ctAssocsTop, ctd=ctd, fileTag=mgl_type)

# plot both together
ctAssocMerged = merge_magma_results(ctAssocsLinear, ctAssocsTop)
FigsMerged = plot_celltype_associations(ctAssocMerged, ctd=ctd, fileTag=mgl_type)
```

### LDSC
LDSC also has excellent [documentation](https://github.com/bulik/ldsc/wiki/Cell-type-specific-analyses), and the instructions in the README in this folder should be sufficient. Here are the example options used for each command:

Generate annotations:
```bash
parallel -j7 make_annot.py \
  --gene-set-file {1}_1000Gv3_ldscores/{1}.{2}.GeneSet \
  --gene-coord-file ENSG_coord.txt \
  --windowsize 100000 \
  --bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.{3}.bim \
  --annot-file {1}_1000Gv3_ldscores/{1}.{2}.{3}.annot.gz \
  ::: $cts_name ::: {1..15} control ::: {1..22}
```

Partition LD Scores:
```bash
parallel -j7 ldsc.py \
  --l2 \
  --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.{3} \
  --ld-wind-cm 1 \
  --annot {1}_1000Gv3_ldscores/{1}.{2}.{3}.annot.gz \
  --thin-annot \
  --out {1}_1000Gv3_ldscores/{1}.{2}.{3} \
  --print-snps hapmap3_snps/hm.{3}.snp \
  ::: $cts_name ::: {1..15} control ::: {1..22}
```

Munge_stats:
```bash
munge_sumstats.py \
  --sumstats ../R_analysis_redux/body_BMIz.sumstats.gz \
  --merge-alleles w_hm3.snplist \
  --out BMI
```

Regression:
```bash
ldsc.py \
  --h2-cts BMI.sumstats.gz \
  --ref-ld-chr 1000G_EUR_Phase3_baseline_v2.2/baselineLD. \
  --overlap-annot \
  --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. \
  --out BMI_${cts_name} \
  --ref-ld-chr-cts ${cts_name}.ldcts \
  --w-ld-chr weights_hm3_no_hla/weights.
```

### In summary
So, there are the general steps to recreate the raw data used in our publication. As a disclaimer, this information is provided as is and likely requires modification to run correctly on other machines. Understand the commands before running them. I am not responsible if you try to run 100 threads on a laptop. 

[Home](README.md)
