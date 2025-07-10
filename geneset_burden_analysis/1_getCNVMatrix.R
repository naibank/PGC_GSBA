#########################################################
#### Install GSBurden ####
### only run a line below when doing this in Rstudio
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
# devtools::install_github("naibank/GSBurden")
# Rscript 1_getCNVMatrix_ASD.R <geneset_path> <sample_info_path> <cnv_path>
args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(GenomicRanges)
library(GSBurden)

gs_path <- args[1]
sample_path <- args[2]
cnv_path <- args[3]
gene_level_sumstat <- args[4]

disorder <- gsub("sampleinfo|\\.|tsv", "", basename(sample_path))

### read in gene definition
gene.in <- read.delim("hg38_refGene_20200708.exon.txt", header = F, 
                      col.names = c("chr", "start", "end", "isoid", "genesymbol", "enzid"), stringsAsFactors = F)
genes <- GeneAnnotation(gene.in$enzid, gene.in$chr, gene.in$start, gene.in$end, gene.in$genesymbol)
genes.g <- GRanges(genes$chr, IRanges(genes$start, genes$end), "*")

gs_variable <- ls()

load(gs_path)
gs_variable <- setdiff(ls(), c(gs_variable, "gs_variable"))
gsData.s <- get(gs_variable)

if(basename(gs_path) == "human_gsData_goincludingIEA100-1000_pathways15-500_20211125.RData"){
  gsData.s <- gsData.s$gs2gene
  gsData.s <- gsData.s[sapply(gsData.s, length) >= 50 & sapply(gsData.s, length) <= 500]
  
  names(gsData.s) <- gsub(":", "_", names(gsData.s))
  names(gsData.s) <- gsub("-", "_", names(gsData.s))
  names(gsData.s) <- gsub("\\.", "_", names(gsData.s))
}
names(gsData.s) <- gsub("\\ |-", ".", names(gsData.s))

sampleinfo <- data.table::fread(sample_path, data.table = F)
names(sampleinfo) <- c("FID", "sample", "status", "sex", "ancestry", "PLATFORM", "dataset", "phenotype", paste0("PC", 1:10))

cnvs.in <- fread(cnv_path, data.table = F)
names(cnvs.in) <- c("FID", "sample", "chr", "start", "end", "type", "score", "site")
cnvs.in <- cnvs.in[which(cnvs.in$type %in% c(1, 3)), ]
cnvs.in$chr <- gsub("23", "X", cnvs.in$chr)
cnvs.in$chr <- gsub("24", "Y", cnvs.in$chr)
cnvs.in$chr <- paste0("chr", cnvs.in$chr)
cnvs.in$type <- ifelse(cnvs.in$type == 1, "DEL", "DUP")

cnvs.in <- cnvs.in[cnvs.in$sample %in% sampleinfo$sample, ]
cnvs.in <- cnvs.in[cnvs.in$chr %in% paste0("chr", c(1:22, "X")),]
cnvs.in$cnvid <- paste(cnvs.in$chr, cnvs.in$start, cnvs.in$end, cnvs.in$type, sep="#")

cnvs.in.g <- GRanges(cnvs.in$chr, IRanges(cnvs.in$start, cnvs.in$end), "*")
olap <- data.frame(findOverlaps(cnvs.in.g, genes.g))
olap$gene <- genes$enzid[olap$subjectHits]
nrxn1_del <- unique(olap$queryHits[genes$gsymbol[olap$subjectHits] == "NRXN1" &
                              cnvs.in$type[olap$queryHits] == "DEL"])
nrxn1_del <- unique(cnvs.in$sample[nrxn1_del])
olap <- unique(olap[, c("queryHits", "gene")])


olap <- data.frame(table(olap$queryHits))

sumstat <- fread(gene_level_sumstat, data.table = F)
del.sumstat <- sumstat[sumstat$type == "DEL", ]
dup.sumstat <- sumstat[sumstat$type == "DUP", ]   

del.retained.genes <- unique(strsplit(paste(del.sumstat$genesymbol[!del.sumstat$filtered_call], collapse = ","), ",")[[1]])
del.excluded.genes <- unique(strsplit(paste(del.sumstat$genesymbol[del.sumstat$filtered_call], collapse = ","), ",")[[1]])
  
dup.retained.genes <-   unique(strsplit(paste(dup.sumstat$genesymbol[!dup.sumstat$filtered_call], collapse = ","), ",")[[1]])
dup.excluded.genes <- unique(strsplit(paste(dup.sumstat$genesymbol[dup.sumstat$filtered_call], collapse = ","), ",")[[1]])

del.exclude <- del.excluded.genes[!del.excluded.genes %in% del.retained.genes]
dup.exclude <- dup.excluded.genes[!dup.excluded.genes %in% dup.retained.genes]

del.matrix <- getCNVGSMatrix(cnvs.in[which(cnvs.in$type == "DEL"),], 
                             genes[which(!genes$gsymbol %in% del.exclude), ], gsData.s)
dup.matrix <- getCNVGSMatrix(cnvs.in[which(cnvs.in$type == "DUP"),], 
                             genes[which(!genes$gsymbol %in% dup.exclude), ], gsData.s)


cnv.matrix <- merge(del.matrix, dup.matrix, by = "sample", all = T)
cnv.matrix[is.na(cnv.matrix)] <- 0

rm(list = c("del.matrix", "dup.matrix"))

platform.case <- table(sampleinfo$PLATFORM[sampleinfo$status == 1])
platform.case <- names(which(platform.case >= 50))
platform.ctrl <- table(sampleinfo$PLATFORM[sampleinfo$status == 0])
platform.ctrl <- names(which(platform.ctrl >= 50))
platform.both <- intersect(platform.case, platform.ctrl)

sampleinfo <- sampleinfo[sampleinfo$PLATFORM %in% platform.both,]

cnv.matrix <- cnv.matrix[cnv.matrix$sample %in% sampleinfo$sample, ]
cnv.matrix <- merge(sampleinfo, cnv.matrix, by = "sample", 
                    all.x = T)
cnv.matrix[is.na(cnv.matrix)] <- 0

write.table(cnv.matrix, sprintf("%s.matrix.all.tsv", disorder), sep="\t", row.names=F, quote=F, col.names=T)

### get cnv.matrix excluding recurrent CNVs + genome-wide significant genes
recurrent_loci <- fread("/RecurrentCNVLoci_AnalysisPlan_Feb2023 - Sheet1.tsv", data.table = F)
recurrent_loci <- recurrent_loci[recurrent_loci$start_hg38 != "exonic deletions", ]
recurrent_loci.g <- GRanges(recurrent_loci$chr, IRanges(as.numeric(recurrent_loci$start_hg38), as.numeric(recurrent_loci$end_hg38)), "*")

recurrent_olap <- data.frame(findOverlaps(cnvs.in.g, recurrent_loci.g))
recurrent_olap$coverage <- width(pintersect(cnvs.in.g[recurrent_olap$queryHits],
                                            recurrent_loci.g[recurrent_olap$subjectHits]))/width(recurrent_loci.g[recurrent_olap$subjectHits])
recurrent_olap$sample <- cnvs.in$sample[recurrent_olap$queryHits]
recurrent_olap$type <- cnvs.in$type[recurrent_olap$queryHits]
recurrent_coverage <- aggregate(coverage ~sample + type + subjectHits, recurrent_olap, sum)
recurrent_coverage <- recurrent_coverage[which(recurrent_coverage$coverage > 0.75), ]
recurrent_coverage$locus <- recurrent_loci$locus_name[recurrent_coverage$subjectHits]
recurrent_coverage <- merge(recurrent_coverage, cnv.matrix[, c("sample", "status")], by = "sample", all.x = T)

excluding.samples <- unique(c(nrxn1_del, recurrent_coverage$sample, strsplit(paste(sumstat$sample[sumstat$genome_wide_significant], collapse = ","), ",")[[1]]  ))
writeLines(excluding.samples, "samples.excluded.with.known.significant.variants.txt")

cnv.matrix <- cnv.matrix[which(!cnv.matrix$sample %in% excluding.samples), ]
write.table(cnv.matrix, sprintf("%s.matrix.smaller.tsv", disorder), sep="\t", row.names=F, quote=F, col.names=T)

