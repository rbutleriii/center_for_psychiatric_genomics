#!/usr/bin/env Rscript

library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)

setwd('/media/sf_E_DRIVE/microglia_iPSC/R_analysis_redux/')
# setwd('/Volumes/botweiser/microglia_iPSC/R_analysis_redux/')
# setwd('/data/butlerr/microglia_iPSC/R_analysis_redux/')
gwas_names <- c("AD", "SCZ", "PD", "BMI")
# file root names for ldsc
gwas_list <- c("Jan_AD", "PGC_SCZ", "Nalls_PD", "BMI")


###########################################
# LDSC
# get files (from 3-15 sumstats)
files <- list.files(path="../ldsc", pattern="2019-03-15_.*[^kb]_stats\\.txt", 
                    full.names=T)
all(file.exists(files))
lapply(files, function(x) {
  assign(gsub(pattern = "../ldsc/2019-03-15_(.*)_stats.txt", 
              replacement = "\\1", x), 
         read.table(file = x, header = T), inherits = T)
})

# merge frames per GWAS (dropping height)
df_list <- lapply(gwas_list, function(x) get(x))
ldsc_result <- df_list[[1]][, c("Name", "Coefficient_P_value_top20") ]
for (i in head(seq_along(df_list), -1)) {
  ldsc_result <- merge(x = ldsc_result, 
                       y = df_list[[i+1]][, 
                           c("Name", "Coefficient_P_value_top20") ],
                       by = "Name", suffixes = gwas_list[i:(i+1)])
}

# cleanup
names(ldsc_result) <- c("VARIABLE", gwas_names)
remove(df_list, files, Height, BMI, Jan_AD, Nalls_PD, PGC_SCZ, i, gwas_list)

###########################################
# Magma
# get files (from 3-15 analysis)
files <- list.files(path="../Magma_analysis",
                    pattern="2019-03-15_.*_stats\\.txt", 
                    full.names=T)
all(file.exists(files))
lapply(files, function(x) {
  assign(gsub(pattern = "../Magma_analysis/2019-03-15_(.*)_stats.txt", 
              replacement = "\\1", x), 
         read.table(file = x, header = T), inherits = T)
})

# restructure stats for each GWAS
lapply(gwas_names, function(y) {
  df <- get(y)
  # Pvals for non-MGL from aMGL column
  var_names <- as.character(df$VARIABLE[ !is.na(df$P_aMGL) ])
  df_short <- df[ df$VARIABLE %in% var_names, c("VARIABLE", "P_aMGL") ]
  # add pvals for MGL-like to x_short list
  mgl_names <- as.character(df$VARIABLE[ is.na(df$P_aMGL) ])
  for (x in mgl_names) {
    df_short[ nrow(df_short) + 1, ] <-
      list(x, df[ df$VARIABLE == x, paste0("P_", x) ])
  }
  # write to _short df
  assign(sprintf("%s_short", y), df_short, inherits = T)
})

# merge dataframe columns
df_list <- lapply(gwas_names, function(x) get(sprintf("%s_short", x)))
magma_result <- df_list[[1]]
for (i in head(seq_along(df_list), -1)) {
  magma_result <- merge(x = magma_result, y = df_list[[i+1]], by = "VARIABLE",
                        suffixes = gwas_names[i:(i+1)])
}

# cleanup
magma_result$VARIABLE <- gsub("\\.", "_", magma_result$VARIABLE)
names(magma_result) <- c("VARIABLE", gwas_names)
remove(df_list, files, Height, BMI, AD, PD, SCZ, AD_short, PD_short, 
       SCZ_short, i)


###########################################
# graph heatmap
# merge datasets
heat_result <- merge(x = ldsc_result, y = magma_result, by = "VARIABLE",
                     suffixes = c("_ldsc", "_magma"))
rownames(heat_result) <- heat_result$VARIABLE
heat_result <- heat_result[ , -1 ]
log_heat <- -log10(heat_result)

# set colors
colors <- colorRampPalette(c("#FFFFFF", brewer.pal(6, "Blues")),
                           bias = 1.5)(255)
colors_bmi <- colorRampPalette(c("#FFFFFF", brewer.pal(6, "Blues")))(255)

# with BMI
pheatmap <- pheatmap(log_heat, col = colors, cluster_rows = T, 
                     cluster_cols = F, gaps_col = c(rep(length(gwas_names),2)), 
                     filename = paste0(Sys.Date(), 
                                       "_enrichment_heatmap_pheatmap.pdf"),
                     width = 11, height = 8.5)

# without BMI
no_bmi <- log_heat[ , c(1:3, 5:7)]
pheatmap(no_bmi, col = colors_bmi, cluster_rows = T, cluster_cols = F,
         filename = paste0(Sys.Date(), "_enrichment_heatmap_noBMI.pdf"),
         width = 11, height = 8.5, gaps_col = c(rep(length(no_bmi)/2,2)))

#Heatmap(log_heat, col = colors, 

# P-value Bar Graph
pdf(file = paste0(Sys.Date(), "AD_Pval_ldsc-magma.pdf"), 
    onefile = T, 
    paper = "USr",
    width = 11,
    height = 8.5)

pval_mat <- as.matrix(log_heat[ rev(pheatmap$tree_row[["order"]]),
                                c('AD_ldsc', 'AD_magma') ])
ggmelt <- melt(pval_mat)
ggmelt$Var2 <- gsub(pattern = "AD_", replacement = "", ggmelt$Var2)
ggAD <- ggplot(ggmelt, aes(x=Var1, y=value, fill=Var2)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  scale_fill_brewer(type = "qual", palette = "Paired") +
  ggtitle("AD P-values") + coord_flip() + scale_y_reverse() +
  labs(x="Cell type", y="-log10(P)", fill="Analysis") + 
  theme_classic()
plot(ggAD)
dev.off()



