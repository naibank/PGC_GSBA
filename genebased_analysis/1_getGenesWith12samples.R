#########################################################
#### Install GSBurden ####
### only run a line below when doing this in Rstudio
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
# devtools::install_github("naibank/GSBurden")
# Rscript 1_getGenesWith12samples.R <sample_path> <cnv_path> <disorder>
args = commandArgs(trailingOnly=TRUE)

sample_path <- args[1]
cnv_path <- args[2]
disorder <- args[3]

library(data.table)
library(GenomicRanges)

#### manifest exclude duplicated missing samples
exclude.samples <- readLines("duplicated_PGC.txt")

#### getsampleinfo
sampleinfo <- fread(sample_path, data.table = F)
sampleinfo <- sampleinfo[which(sampleinfo$AFF %in% c(1,2) & sampleinfo$SEX %in% c(1, 2)), ]
sampleinfo$AFF <- ifelse(sampleinfo$AFF == 2, 1, 0)
sampleinfo <- sampleinfo[!sampleinfo$IID %in% exclude.samples, ]
write.table(sampleinfo, sprintf("sampleinfo.%s.tsv", disorder), sep="\t", row.names=F, quote = F, col.names=T)

##### cnv data
cnvs <- fread(cnv_path, data.table = F)
cnvs <- cnvs[cnvs$TYPE %in% c(1, 3), ]
cnvs$CHR <- paste0("chr", cnvs$CHR)
cnvs$CHR <- gsub("chr23", "chrX", cnvs$CHR)
cnvs$CHR <- gsub("chr24", "chrY", cnvs$CHR)
cnvs$TYPE <- ifelse(cnvs$TYPE == 1, "DEL", "DUP")
cnvs <- cnvs[cnvs$IID %in% sampleinfo$IID, ]
cnvs <- cnvs[cnvs$CHR %in% paste0("chr", c(1:22, "X")),]
del <- cnvs[cnvs$TYPE == "DEL", ]
dup <- cnvs[cnvs$TYPE == "DUP", ]

del.g <- GRanges(del$CHR, IRanges(del$BP1, del$BP2), "*")
dup.g <- GRanges(dup$CHR, IRanges(dup$BP1, dup$BP2), "*")

genes <- fread("hg38_refGene_20200708.exon.txt", data.table = F)
genes.g <- GRanges(genes$V1, IRanges(genes$V2, genes$V3), "*")

if(!dir.exists("intermediateFiles"))
  dir.create("intermediateFiles")

del.olap <- data.frame(findOverlaps(del.g, genes.g))
del.olap$gene <- genes$V5[del.olap$subjectHits]
del.olap$sample <- del$IID[del.olap$queryHits]
del.olap <- unique(del.olap[, c("sample", "gene")])
del.genes <- names(which(table(del.olap$gene) >= 12))
writeLines(del.genes, "intermediateFiles/Genes_with12DELsOrMore.txt")

dup.olap <- data.frame(findOverlaps(dup.g, genes.g))
dup.olap$gene <- genes$V5[dup.olap$subjectHits]
dup.olap$sample <- dup$IID[dup.olap$queryHits]
dup.olap <- unique(dup.olap[, c("sample", "gene")])
dup.genes <- names(which(table(dup.olap$gene) >= 12))
writeLines(dup.genes, "intermediateFiles/Genes_with12DUPsOrMore.txt")
