library(data.table)
library(GSBurden)
library(GenomicRanges)
library(metafor)

args = commandArgs(trailingOnly=TRUE)

analysis <- args[1]
disorder <- args[2]
gs_path <- args[3]

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


ans <- toupper(strsplit(disorder, "_")[[1]][2])
if(is.na(ans)){
  ans <- "ALL"
}

disease <- toupper(strsplit(disorder, "_")[[1]][1])
if(length(list.files(sprintf("BurdenIndividualPlatform/%s/%s_FAM/results/", analysis, toupper(disorder)),
                     full.names = T)) > 0 &
   length(list.files(sprintf("BurdenIndividualPlatform/%s/%s/results/", analysis, toupper(disorder)),
                     full.names = T)) > 0){
  
  metadata <- fread(sprintf("BurdenIndividualPlatform/genebased_mega_analysis/%s/sampleinfo.%s.tsv", toupper(disorder), disorder), data.table = F)
  metadata_fam <- fread(sprintf("BurdenIndividualPlatform/genebased_mega_analysis/%s_FAM/sampleinfo.%s_fam.tsv", toupper(disorder), disorder), data.table = F)
  
  cnvs <- rbind(fread(sprintf("Data/CNVPlinkFormat/%s/%s/Fam/all_platforms_10kb_10probes_Fams_%s_%s.cnv",
                              ans, disease, disease, ans), data.table = F),
                fread(sprintf("Data/CNVPlinkFormat/%s/%s/noFam/all_platforms_10kb_10probes_noFams_%s_%s.cnv", 
                              ans, disease, disease, ans), data.table = F))
  
  explained_samples <- c(readLines(sprintf("BurdenIndividualPlatform/%s/%s/samples.excluded.with.known.significant.variants.txt",
                                analysis, disease)),
                        readLines(sprintf("BurdenIndividualPlatform/%s/%s_FAM/samples.excluded.with.known.significant.variants.txt",
                                analysis, disease)))
  
  gs_sum_path <- c(list.files(sprintf("BurdenIndividualPlatform/%s/%s/results/", analysis, toupper(disorder)),
                              full.names = T),
                   list.files(sprintf("BurdenIndividualPlatform/%s/%s_FAM/results/", analysis, toupper(disorder)),
                              full.names = T))
  
  
}else if(length(list.files(sprintf("BurdenIndividualPlatform/%s/%s/results/", analysis, toupper(disorder)),
                           full.names = T)) > 0){
  
  metadata <- fread(sprintf("BurdenIndividualPlatform/genebased_mega_analysis/%s/sampleinfo.%s.tsv", toupper(disorder), disorder), data.table = F)

  cnvs <- fread(sprintf("Data/CNVPlinkFormat/%s/%s/noFam/all_platforms_10kb_10probes_noFams_%s_%s.cnv",
                        ans, disease, disease, ans), data.table = F)
  
  explained_samples <- readLines(sprintf("BurdenIndividualPlatform/%s/%s/samples.excluded.with.known.significant.variants.txt",
                                analysis, disease))
  
  gs_sum_path <- list.files(sprintf("BurdenIndividualPlatform/%s/%s/results/", analysis, toupper(disorder)),
                              full.names = T)
}else{
  metadata_fam <- fread(sprintf("BurdenIndividualPlatform/genebased_mega_analysis/%s_FAM/sampleinfo.%s_fam.tsv", toupper(disorder), disorder), data.table = F)
  
  cnvs <- fread(sprintf("Data/CNVPlinkFormat/%s/%s/Fam/all_platforms_10kb_10probes_Fams_%s_%s.cnv",
                        ans, disease, disease, ans), data.table = F)
  
  explained_samples <- readLines(sprintf("BurdenIndividualPlatform/%s/%s_FAM/samples.excluded.with.known.significant.variants.txt",
                                         analysis, disease))
  
  gs_sum_path <- list.files(sprintf("BurdenIndividualPlatform/%s/%s_FAM/results/", analysis, toupper(disorder)),
                            full.names = T)
}

gene_level_sumstat <- fread(sprintf("MetaAnalysis/genebased_mega_analysis/%s/%s_final_summ_stat.tsv", toupper(disorder),
                                    toupper(disorder)), data.table = F)

del.sumstat <- gene_level_sumstat[gene_level_sumstat$type == "DEL", ]
dup.sumstat <- gene_level_sumstat[gene_level_sumstat$type == "DUP", ]   

del.retained.genes <- unique(strsplit(paste(del.sumstat$genesymbol[!del.sumstat$filtered_call], collapse = ","), ",")[[1]])
del.excluded.genes <- unique(strsplit(paste(del.sumstat$genesymbol[del.sumstat$filtered_call], collapse = ","), ",")[[1]])

dup.retained.genes <-   unique(strsplit(paste(dup.sumstat$genesymbol[!dup.sumstat$filtered_call], collapse = ","), ",")[[1]])
dup.excluded.genes <- unique(strsplit(paste(dup.sumstat$genesymbol[dup.sumstat$filtered_call], collapse = ","), ",")[[1]])

del.exclude <- del.excluded.genes[!del.excluded.genes %in% del.retained.genes]
dup.exclude <- dup.excluded.genes[!dup.excluded.genes %in% dup.retained.genes]


names(cnvs) <- c("FID", "sample", "chr", "start", "end", "type", "score", "site")
cnvs <- cnvs[which(cnvs$type %in% c(1,3)), ]
cnvs$chr <- gsub("23", "X", cnvs$chr)
cnvs$chr <- gsub("24", "Y", cnvs$chr)
cnvs$chr <- paste0("chr", cnvs$chr)
cnvs$type <- ifelse(cnvs$type == 1, "DEL", "DUP")

#### meta analysis
final_stats <- data.frame()
for(set in c("all", "smaller")){
  message(set)
  
  cnvs.here <- cnvs
  if(set == "smaller"){
    cnvs.here <- cnvs.here[!cnvs.here$sample %in% explained_samples, ]
  }
  cnvs.g <- GRanges(cnvs.here$chr, IRanges(cnvs.here$start, cnvs.here$end), "*")
  olap <- data.frame(findOverlaps(cnvs.g, genes.g))
  olap$sample <- cnvs$sample[olap$queryHits]
  olap$type <- cnvs$type[olap$queryHits]
  olap$enzid <- genes$enzid[olap$subjectHits]
  olap$gsymbol <- genes$gsymbol[olap$subjectHits]
  olap <- unique(olap[, c("sample", "type", "enzid", "gsymbol")])
  
  olap <- olap[(olap$type == "DEL" & (!olap$gsymbol %in% del.exclude)) |
                 (olap$type == "DUP" & (!olap$gsymbol %in% dup.exclude)), ]
  
  if(length(gs_sum_path) > 1){ ### when there is multiple platform
    gs_sum_stats <- data.frame()
    gs_sum_weights <- list()
    gs_se_dat <- data.frame()
    for(path in gs_sum_path){
      
      platform.folder <- basename(path)
      platform <- gsub("\\.|\\-", "_", basename(path))
      platform <- gsub("\\+", "PlusSign", platform)
      
      if(length(grep("FAM", path)) > 0){
        metadata.here <- metadata_fam
        disorder.here <- paste0(disorder, "_fam")
        platform <- paste0(platform, "_fam")
      }else{
        metadata.here <- metadata
        disorder.here <- disorder
      }
     
      sumstat <- data.frame()
      for(file in list.files(path, full.names = T, pattern = sprintf("\\.%s\\.", set))){
        message(file)
        sumstat <- rbind(sumstat, fread(file, data.table = F))
      } 
      
      sumstat$platform <- platform
      gs_se_dat <- rbind(gs_se_dat, sumstat[, c("geneset", "type", "platform", "coefficient", "stderr")])
      names(sumstat)[c(3, 6)] <- paste(names(sumstat)[c(3, 6)], platform, sep="_")
      sumstat <- sumstat[, c(1:3, 6)]
      if(nrow(sumstat) == 0){
        if(ncol(gs_sum_stats) == 0){
          gs_sum_stats <- sumstat
        }else{
          gs_sum_stats[, c(paste0(c("coefficient_", "pvalue_"), platform))] <- NA
        }
      }else{
        if(ncol(gs_sum_stats) == 0){
          gs_sum_stats <- sumstat
        }else{
          gs_sum_stats <- merge(gs_sum_stats, sumstat, by=c("geneset", "type"), all = T)
        }
      }
      
      ncase <- sum(metadata.here$AFF == 1 & metadata.here$CNV_platform == platform.folder) 
      nctrl <- sum(metadata.here$AFF == 0 & metadata.here$CNV_platform == platform.folder)
      N = 4/((1/ncase)+(1/nctrl))
      w = sqrt(N)
      
      gs_sum_weights[[platform]] <- w
    }
    z <- gs_sum_weights 
    for(platform in names(gs_sum_weights)){
      z.tmp <- sign(gs_sum_stats[, paste0("coefficient_", platform)]) * -1 * qnorm(gs_sum_stats[, paste0("pvalue_", platform)]/2) * 
        gs_sum_weights[[platform]]
      
      gs_sum_stats[, sprintf("z_%s", platform)] <- z.tmp
    }
    
    num = sapply(as.data.frame.matrix(t(gs_sum_stats[, paste0("z_", names(gs_sum_weights))])), sum, na.rm = T)
    
    for(i in 1:length(num)){
      if(length(which(!is.na(gs_sum_stats[i, paste0("z_", names(gs_sum_weights))]))) > 0){
        den = sqrt(sum(sapply(gs_sum_weights[which(!is.na(gs_sum_stats[i, paste0("z_", names(gs_sum_weights))]))], "^", 2)))
        
      }else{
        den = sqrt(sum(sapply(gs_sum_weights, "^", 2)))
        
      }
      num[i] <- num[i]/den
    }
    
    gs_sum_stats$z = num
    gs_sum_stats$pvalue = 2*pnorm(pmin(gs_sum_stats$z,-gs_sum_stats$z))
    gs_sum_stats$set <- set
    
    gs_sum_stats$se_test <- NA
    gs_sum_stats$se_coeff <- NA
    gs_sum_stats$se_ci_ub <- NA
    gs_sum_stats$se_ci_lb <- NA
    
    gs_sum_stats$genes <- ""
    for(i in 1:nrow(gs_sum_stats)){
      tmp.gs <- gsData.s[[gs_sum_stats$geneset[i]]]
      genes.here <- olap$gsymbol[olap$enzid %in% tmp.gs & olap$type == gs_sum_stats$type[i]]
      
      gs_sum_stats$genes[i] <- paste(unique(genes.here), collapse = ",")
      
      tmp.dat <- gs_se_dat[gs_se_dat$type == gs_sum_stats$type[i] &
                             gs_se_dat$geneset == gs_sum_stats$geneset[i], ]
      if(nrow(na.omit(tmp.dat)) > 0){
        tmp <- tryCatch(
          {res <- rma(yi = coefficient, sei = stderr, data = tmp.dat,
                      control=list(stepadj=0.5, maxiter=10000))},
          error = function(cond){
            res <- list()
            res$zval <- res$beta <- res$ci.ub <- res$ci.lb <- NA
            message("error")
            return(res)
          }
        )
        
        gs_sum_stats$se_test[i] <- tmp$zval
        gs_sum_stats$se_coeff[i] <- tmp$beta
        gs_sum_stats$se_ci_ub[i] <- tmp$ci.ub
        gs_sum_stats$se_ci_lb[i] <- tmp$ci.lb
      }else{
        gs_sum_stats$se_test[i] <- NA
        gs_sum_stats$se_coeff[i] <- NA
        gs_sum_stats$se_ci_ub[i] <- NA
        gs_sum_stats$se_ci_lb[i] <- NA
      }
      
    }
    
    final_stats <- rbind(final_stats, gs_sum_stats)
    
  }else{ ### when there is only one platform
    path = gs_sum_path
    
    platform.folder <- basename(gs_sum_path[1])
    platform <- gsub("\\.|\\-", "_", basename(gs_sum_path[1]))
    platform <- gsub("\\+", "Plus", platform)
    
    sumstat <- data.frame()
    for(file in list.files(path, full.names = T, pattern = sprintf("\\.%s\\.", set))){
      message(file)
      sumstat <- rbind(sumstat, fread(file, data.table = F))
    } 
    
    if(length(grep("FAM", path)) > 0){
      metadata.here <- metadata_fam
      disorder.here <- paste0(disorder, "_fam")
      platform <- paste0(platform, "_fam")
    }else{
      metadata.here <- metadata
      disorder.here <- disorder
    }
    
    names(sumstat)[c(3, 6)] <- paste(names(sumstat)[c(3, 6)], platform, sep="_")
    gs_sum_stats <- sumstat[, c(1:3, 6)]
    
    ncase <- sum(metadata.here$AFF == 1 & metadata.here$CNV_platform == platform.folder) 
    nctrl <- sum(metadata.here$AFF == 0 & metadata.here$CNV_platform == platform.folder)
    N = 4/((1/ncase)+(1/nctrl))
    w = sqrt(N)
    
    z.tmp <- sign(gs_sum_stats[, paste0("coefficient_", platform)]) * -1 * qnorm(gs_sum_stats[, paste0("pvalue_", platform)]/2) * w
    den <- w
    gs_sum_stats$z = z.tmp/w
    gs_sum_stats$pvalue = 2*pnorm(pmin(gs_sum_stats$z,-gs_sum_stats$z))
    gs_sum_stats$set <- set
    
    gs_sum_stats$se_test <- sumstat$testval
    gs_sum_stats$se_coeff <- sumstat$coefficient
    gs_sum_stats$se_ci_ub <- sumstat$coefficient.upper
    gs_sum_stats$se_ci_lb <- sumstat$coefficient.lower
    
    gs_sum_stats$genes <- ""
    for(i in 1:nrow(gs_sum_stats)){
      tmp.gs <- gsData.s[[gs_sum_stats$geneset[i]]]
      genes.here <- olap$gsymbol[olap$enzid %in% tmp.gs & olap$type == gs_sum_stats$type[i]]
      
      gs_sum_stats$genes[i] <- paste(unique(genes.here), collapse = ",")
    }
    
    final_stats <- rbind(final_stats, gs_sum_stats)
  }
}

final_stats$qvalue <- NA
final_stats$qvalue[final_stats$set == "all" & final_stats$type == "DEL"] <- p.adjust(final_stats$pvalue[final_stats$set == "all" & final_stats$type == "DEL"], method = "BH")
final_stats$qvalue[final_stats$set == "all" & final_stats$type == "DUP"] <- p.adjust(final_stats$pvalue[final_stats$set == "all" & final_stats$type == "DUP"], method = "BH")
final_stats$qvalue[final_stats$set == "smaller" & final_stats$type == "DEL"] <- p.adjust(final_stats$pvalue[final_stats$set == "smaller" & final_stats$type == "DEL"], method = "BH")
final_stats$qvalue[final_stats$set == "smaller" & final_stats$type == "DUP"] <- p.adjust(final_stats$pvalue[final_stats$set == "smaller" & final_stats$type == "DUP"], method = "BH")

write.table(final_stats, sprintf("MetaAnalysis/%s/%s/meta_analysis_all_tests.tsv", 
                                 analysis, toupper(disorder)), sep="\t", row.names=F, quote=F, col.names=T)
