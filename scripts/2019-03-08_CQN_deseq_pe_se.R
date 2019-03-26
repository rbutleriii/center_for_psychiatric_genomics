#!/usr/bin/env Rscript

library(cqn)
library(edgeR)
library(dplyr)
library(biomaRt)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(PoiClaClu)


#####################################
# Kozlova data prep
# load featureCounts tables for Kozlova data
setwd('/Volumes/botweiser/microglia_iPSC/R_analysis_redux/')
# setwd('/media/sf_E_DRIVE/microglia_iPSC/R_analysis_redux/')
fcounts_pe <- read.table("featureCounts_pe_counts.txt", header = T)
fcounts_se <- read.table("featureCounts_se_counts.txt", header = T)

#unused column indices
uci <- -c(2:6)

# merge pe and se datasets for all fragment counts.
all_frag_count <- merge(fcounts_pe[ , uci], fcounts_se[ , uci], by = "Geneid")

# grabbing gene info for analysis from featurecounts, biomaRt
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
BMdat <- getBM(attributes = c("ensembl_gene_id",
                              "external_gene_name",
                              "percentage_gene_gc_content",
                              "gene_biotype",
                              "chromosome_name"),
               filters = "ensembl_gene_id",
               values = all_frag_count$Geneid,
               mart = ensembl)

# merging into features table for length, gc, and gene name
gene_features <- merge(fcounts_pe[order(fcounts_pe$Geneid), c(1,6)], 
                       BMdat,
                       by.x = "Geneid",
                       by.y = "ensembl_gene_id",
                       all.x = T)

# no missing gc values (0 genes)
nogc <- as.vector(gene_features[is.na(gene_features$percentage_gene_gc_content), "Geneid"])
nogc_counts <- all_frag_count[all_frag_count$Geneid %in% nogc,]

# keep only protein coding genes (19,951)
coding <- as.vector(gene_features[gene_features$gene_biotype == "protein_coding", "Geneid"])
coding_counts <- all_frag_count[all_frag_count$Geneid %in% coding,]

# 19,951 coding genes (of 58,735), keep intersection of coding and !(nogc) for final counts
final_frag_count <- all_frag_count[!(all_frag_count$Geneid %in% nogc) &
                                     all_frag_count$Geneid %in% coding, ]
final_gene_features <- gene_features[!(gene_features$Geneid %in% nogc) &
                                       gene_features$Geneid %in% coding, ]

# Cleaning up table
sum(final_frag_count$Geneid == final_gene_features$Geneid) == length(final_frag_count$Geneid)
row.names(final_frag_count) <- final_gene_features$Geneid
row.names(final_gene_features) <- final_gene_features$Geneid
final_frag_count <- final_frag_count[, -1]
final_frag_count <- final_frag_count %>%
  rename_at(.vars = vars(ends_with("_starAligned.sortedByCoord.out.bam")),
            .funs = list(~sub("_starAligned.sortedByCoord.out.bam$", "", .)))

# Rename samples
name <- read.table("new_names.txt", header = T)
# Loop through replacing sample names
for(i in seq_along(name$Run)) {
  final_frag_count <- final_frag_count %>%
    rename_at(.vars = vars(starts_with("SRR")),
              .funs = list(~sub(name$Run[i], name$short_name[i], .)))
}
for(i in seq_along(name$Run)) {
  final_frag_count <- final_frag_count %>%
    rename_at(.vars = vars(starts_with("ERR")),
              .funs = list(~sub(name$Run[i], name$short_name[i], .)))
}

# freeing up some memory
remove(fcounts_pe, fcounts_se, all_frag_count, BMdat, ensembl, gene_features,
       nogc, nogc_counts, uci, name, coding, coding_counts, i)


#####################################
# Filtering and ordering Kozlova datasets
# load sampleTable
sampleTable <- read.table("2019-01-08_sampleTable.txt", header = T)

# reorder columns like sampleTable
final_frag_count <- final_frag_count[ , as.character(sampleTable$Name) ]

#check all there
names(final_frag_count) == as.character(sampleTable$Name)

#####################################
# Cqn
# scale factor as sums of each column (read depth per sample)
# final_frag_count <- final_frag_count + 1
read_depth <- as.array(colSums(final_frag_count))
row.names(read_depth) <- colnames(final_frag_count)

# normalizing with cqn
reads_cqn_source <- cqn(final_frag_count, 
                        lengths = final_gene_features$Length,
                        x = final_gene_features$percentage_gene_gc_content, 
                        sizeFactors = read_depth,
                        verbose = T)

# FPKMs cqn normalized 
FPKM_cqn_table <- reads_cqn_source$y + reads_cqn_source$offset

# Raw FPKMs
FPKM_uncorrected_table <- rpkm(final_frag_count,
                               gene.length = final_gene_features$Length,
                               log = T,
                               prior.count = 1)


#####################################
# write out raw and cqn adjusted super tables (compressed)

# uncorrected FPKM
f1 <- gzfile("FPKM_uncorrected_all.txt.gz", open = "wb")
write.table(merge(final_gene_features[, c(1,3)],
                  FPKM_uncorrected_table,
                  by.x = "Geneid",
                  by.y = 0),
            file = f1,
            eol = "\n",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)
close(f1) 
remove(f1)

# cqn FPKM
f1 <- gzfile("FPKM_cqn_all.txt.gz", open = "wb")
write.table(merge(final_gene_features[, c(1,3)],
                  FPKM_cqn_table,
                  by.x = "Geneid",
                  by.y = 0),
            file = f1,
            eol = "\n",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)
close(f1) 
remove(f1)


#####################################
# check normalization to see issues
pdf(file = paste0(Sys.Date(), "_cqn_systematic_plots.pdf"), 
    onefile = F, 
    paper = "US",
    width = 11,
    height = 8.5)
par(mfrow=c(3,2))
cqnplot(reads_cqn_source, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
cqnplot(reads_cqn_source, n = 2, xlab = "length", lty = 1, ylim = c(1,7))

# MA plots to check normalization
# just split down the middle
grp1 <- c(colnames(final_frag_count[1:26]))
grp2 <- c(colnames(final_frag_count[27:52]))

# eliminate genes FPKM less than 2
whGenes <- which(rowMeans(FPKM_uncorrected_table) >= 2)
M.std <- rowMeans(FPKM_uncorrected_table[whGenes, grp1]) - rowMeans(FPKM_uncorrected_table[whGenes, grp2])
A.std <- rowMeans(FPKM_uncorrected_table[whGenes,])
M.cqn <- rowMeans(FPKM_cqn_table[whGenes, grp1]) - rowMeans(FPKM_cqn_table[whGenes, grp2])
A.cqn <- rowMeans(FPKM_cqn_table[whGenes,])

# plotting MA for overall correction
plot(A.std, M.std, cex = 0.5, pch = 16, xlab = "A", ylab = "M",
  main = "Standard RPKM", ylim = c(-4,4), xlim = c(0,12),
  col = alpha("black", 0.25))
plot(A.cqn, M.cqn, cex = 0.5, pch = 16, xlab = "A", ylab = "M",
  main = "CQN normalized RPKM", ylim = c(-4,4), xlim = c(0,12),
  col = alpha("black", 0.25))

# plotting MA for gc extremes
gccontent <- final_gene_features$percentage_gene_gc_content[whGenes]
whHigh <- which(gccontent > quantile(gccontent, 0.9))
whLow <- which(gccontent < quantile(gccontent, 0.1))
plot(A.std[whHigh], M.std[whHigh], cex = 0.5, pch = 16, xlab = "A",
  ylab = "M", main = "Standard RPKM",
  ylim = c(-4,4), xlim = c(0,12), col = "red")
points(A.std[whLow], M.std[whLow], cex = 0.5, pch = 16, col = "blue")
plot(A.cqn[whHigh], M.cqn[whHigh], cex = 0.5, pch = 16, xlab = "A",
  ylab = "M", main = "CQN normalized RPKM",
  ylim = c(-4,4), xlim = c(0,12), col = "red")
points(A.cqn[whLow], M.cqn[whLow], cex = 0.5, pch = 16, col = "blue")
dev.off()

# cleanup
remove(A.cqn, A.std, M.cqn, M.std, whGenes, whHigh, whLow, grp1, grp2,
       gccontent, FPKM_uncorrected_table, final_frag_count)


#####################################
# Poisson pairwise heatmap analysis 
poisd <- PoissonDistance(t(reads_cqn_source$counts))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( sampleTable$Name, sampleTable$cell_type, sep=" - " )
colnames(samplePoisDistMatrix) <- paste( sampleTable$Name, sampleTable$cell_type, sep=" - " )

# set colors and draw
colors <- colorRampPalette( brewer.pal(9, "RdYlBu") )(255)
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors,
         filename = paste0(Sys.Date(), "_featureCounts_poisson.pdf"),
         width = 11,
         height = 8.5)

remove(poisd)
remove(samplePoisDistMatrix)


#####################################
# AD gene_set heatmap for cqn_FPKMs
# selecting only AD genes form haenseler et al
haenseler <- scan(file = "Haenseler_genes.txt", what = "", sep = "\n")

# only first 31 are AD
haenseler_AD <- haenseler[1:31]

# check in our genes and fetch subset use_genes
haenseler_AD %in% final_gene_features$external_gene_name
use_genes <- final_gene_features[ final_gene_features$external_gene_name %in% haenseler_AD, ]

# setting up gene matrix
mat  <- FPKM_cqn_table[ as.character(use_genes$Geneid) , ]
colnames(mat) == as.character(sampleTable$Name)
colnames(mat) <- paste(sampleTable$Name, sampleTable$Run, sep=" - ")
rownames(mat) == as.character(use_genes$Geneid)
rownames(mat) <- paste(use_genes$Geneid, use_genes$external_gene_name, sep = " - ")
# order by rowSum expression
mat <- mat[ order(rowSums(mat[, c("aMGL_1 - SRR4450453",
                                  "aMGL_2 - SRR4450454",
                                  "aMGL_3 - SRR4450455")]), decreasing = T), ]

# cell type annotations
anno <- as.data.frame(sampleTable$cell_type)
colnames(anno) <- 'Cell Type'
rownames(anno) <- colnames(mat)

# regular heatmap
pheatmap(mat, annotation_col = anno, color = rev(colors), show_rownames = T,
         filename = paste0(Sys.Date(), "_featureCounts_cqn_gene_clust.pdf"),
         width = 11, height = 8.5, cluster_rows = F)

# row centered heatmap
pheatmap(mat, annotation_col = anno, color = rev(colors), show_rownames = T,
         filename = paste0(Sys.Date(), "_featureCounts_cqn_rowcentered_gene_clust.pdf"),
         width = 11, height = 8.5, scale = "row", cluster_rows = F)

remove(haenseler, haenseler_AD, use_genes, mat, anno, colors)


#####################################
# CQN normalized FPKM PCA analysis

# use rowMin value as filter (1 cpm in at least 3 col), 15,469 genes
FPKM_cqn_filtered <- FPKM_cqn_table[rowSums(cpm(reads_cqn_source$counts,
                                                lib.size = reads_cqn_source$sizeFactors)>1) >= 3, ]

# sample batches for coloring
colnames(FPKM_cqn_filtered) == as.character(sampleTable$Name)
batch = sampleTable$cell_type

# starting pdf
num_genes <- nrow(FPKM_cqn_filtered)
pdf(file = paste0(Sys.Date(), "_FPKM_cqn_filtered_pca.pdf"), onefile = T, 
    paper = "USr", width = 11, height = 8.5,
    title = paste0(num_genes, " genes"))

# FPKM_cqn_filtered
pca <- prcomp(t(FPKM_cqn_filtered), scale. = T)
percent_pc <- round((pca$sdev^2 / sum(pca$sdev^2))*100, 2)

# PC1 v PC2
ggPCA1 <- ggplot(data.frame(pca$x), aes(x = PC1, y = PC2, col=batch)) +
  labs(x = paste0("PC1 - ", percent_pc[1], "%"),
       y = paste0("PC2 - ", percent_pc[2], "%"), col = "Study") +
  geom_point(size=3, alpha=0.5) +
  geom_text_repel(aes(label=rownames(data.frame(pca$x)))) +
  scale_color_brewer(type = "qual", palette = "Paired") +
  ggtitle(paste0(num_genes, " genes")) + theme_classic()

# PC1 v PC3
ggPCA2 <- ggplot(data.frame(pca$x), aes(x = PC1, y = PC3, col=batch)) +
  labs(x = paste0("PC1 - ", percent_pc[1], "%"),
       y = paste0("PC3 - ", percent_pc[3], "%"), col = "Study") +
  geom_point(size=3, alpha=0.5) +
  geom_text_repel(aes(label=rownames(data.frame(pca$x)))) +
  scale_color_brewer(type = "qual", palette = "Paired") +
  ggtitle(paste0(num_genes, " genes")) + theme_classic()
plot(ggPCA1)
plot(ggPCA2)
dev.off()


#####################################
# Save FPKM data for later
save(FPKM_cqn_filtered, file = paste(Sys.Date(), "FPKM_cqn_filtered.txt.gz", sep = "_"), compress = T)




