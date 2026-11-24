rm(list=ls())

library(readxl)
library(tidyr)
library(dplyr)
library(matrixStats)
library(pheatmap)

# Immune filtering of metastatic castration-resistant
# prostate cancer (mCRPC) tumor samples.
#
# Update in_dir to specify the location of the TPM expression
# matrices and the immune cell marker gene list.
#
# All results and plots are saved in the current working directory.


#########################################
#### Set directories

# Input directory containing TPM expression files (TXT)
# and the immune cell marker gene list (XLSX)
in_dir <- "D:/SharedDesktopFiles/SharedDesktopFiles"

# Output directory where results and plots are saved
# Default is current working die.  update to folder where resullts should be saved

out_dir <- getwd()


#########################################
#### Load needed data

# Two RNA-seq gene expression matrices comprising 210 mCRPC
# tumor biopsies from the WCDT and PROMOTE cohorts,
# as analyzed by Lundberg et al. (PMID: 37289025).

# WCDT gene expression data (TPM)
xa <- read.table(file.path(in_dir, "Lundberg_at_al/counts_and_TPMs/Counts_TPM_Genes.txt"),
                 header=TRUE, check.names=FALSE)

# PROMOTE gene expression data (TPM)
xb <- read.table(file.path(in_dir, "Lundberg_at_al/counts_and_TPMs/Counts_TPM_Genes_PROMOTE.txt"),
                 header=TRUE, check.names=FALSE)

# Merge gene expression datasets by FEATURE_ID
x1 <- merge(xa[,colnames(xa) != "gene_name"], xb, by="FEATURE_ID")

rm(xa, xb)


#########################################
#### Load immune marker genes

# Immune cell marker gene list from PMID: 29230012 (XLSX)
g <- read_excel(file.path(in_dir, "ImmuneGenes_PMID_29230012.xlsx"))

g <- separate_rows(g, Required_Genes, sep=",")
g$Required_Genes <- trimws(g$Required_Genes)


#########################################
#### Subset immune marker expression

x2 <- x1[x1$gene_name %in% g$Required_Genes,]
x2 <- x2[!duplicated(x2$gene_name),]

sample_cols <- setdiff(colnames(x2), c("FEATURE_ID", "gene_name"))

t1 <- data.frame(t(x2[,sample_cols]), check.names=FALSE)
colnames(t1) <- x2$gene_name

# Verify sample count and expression values
stopifnot(nrow(t1) == 210,
          all(is.finite(as.matrix(t1))),
          all(t1 >= 0))


#########################################
#### Remove genes with >=50% zero expression

t1 <- t1[,colMeans(t1 == 0) < 0.5, drop=FALSE]


#########################################
#### Log2 transformation and Z-scoring

# Log2 transformation of TPM expression values
n.exp <- log2(t1 + 0.001)

# Gene-wise Z-score normalization across all samples
z.exp <- n.exp

for (i in 1:ncol(n.exp)) {
  
  m <- mean(n.exp[,i])
  s <- sd(n.exp[,i])
  
  if (s == 0) {
    s <- 0.001
  }
  
  z.exp[,i] <- (n.exp[,i] - m)/s
}


#########################################
#### Median Z-score for each immune cell type

cells <- c("Tcells", "CD4+ Tcells", "CD8a+ Tcells",
           "CD8b+ Tcells", "Tregs", "B cell",
           "Macrophage/Monocyte", "Dendritic Cell",
           "Natural Killer Cell")

NewDF <- data.frame(row.names=rownames(z.exp))

for (i in cells) {
  
  genes <- g$Required_Genes[g$Cell_type == i]
  
  tdf <- z.exp[,colnames(z.exp) %in% genes, drop=FALSE]
  
  stopifnot(ncol(tdf) > 0)
  
  NewDF[,i] <- rowMedians(as.matrix(tdf))
}


#########################################
#### Remove samples with any immune score >1

# Identify samples with a median Z-score >1
# for at least one immune cell type
test <- rowSums(NewDF > 1)

# Retain samples with all immune cell scores <=1
filtered_rows <- NewDF[test == 0,,drop=FALSE]

# Dimensions before and after immune filtering
print(dim(NewDF))
print(dim(filtered_rows))


#########################################
#### Save immune scores and filtered samples

# RDS file used for downstream stemness analysis
saveRDS(filtered_rows,
        file.path(out_dir, "ImmuneFiltered_Final.rds"))

# Immune cell scores before filtering (CSV)
write.csv(NewDF,
          file.path(out_dir, "ImmuneScores_AllSamples.csv"))

# Immune cell scores after filtering (CSV)
write.csv(filtered_rows,
          file.path(out_dir, "ImmuneScores_Filtered.csv"))


#########################################
#### Heatmap before immune filtering

pheatmap(t(NewDF),
         cluster_rows=FALSE,
         cluster_cols=TRUE,
         show_colnames=FALSE,
         scale="none",
         color=colorRampPalette(c("blue","white","red"))(100),
         filename=file.path(out_dir, "ImmuneScores_AllSamples.pdf"),
         width=12,
         height=5)


#########################################
#### Heatmap after immune filtering

pheatmap(t(filtered_rows),
         cluster_rows=FALSE,
         cluster_cols=TRUE,
         show_colnames=FALSE,
         scale="none",
         color=colorRampPalette(c("blue","white","red"))(100),
         filename=file.path(out_dir, "ImmuneScores_Filtered.pdf"),
         width=12,
         height=5)