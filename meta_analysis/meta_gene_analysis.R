library(data.table)
library(metafor)

args = commandArgs(trailingOnly=TRUE)

analysis <- args[1]
disorder <- args[2]

gene_path <- "hg38_refGene_20200708.transcript.txt"


if(length(list.files(sprintf("BurdenIndividualPlatform/%s/%s_FAM/results/", analysis, toupper(disorder)),
                     full.names = T)) > 0 &
   length(list.files(sprintf("BurdenIndividualPlatform/%s/%s/results/", analysis, toupper(disorder)),
                     full.names = T)) > 0){
  metadata <- fread(sprintf("BurdenIndividualPlatform/genebased_mega_analysis/%s/sampleinfo.%s.tsv", toupper(disorder), disorder), data.table = F)
  metadata_fam <- fread(sprintf("BurdenIndividualPlatform/genebased_mega_analysis/%s_FAM/sampleinfo.%s_fam.tsv", toupper(disorder), disorder), data.table = F)
  
  gene_sum_path <- c(list.files(sprintf("BurdenIndividualPlatform/%s/%s/results/", analysis, toupper(disorder)),
                                full.names = T),
                     list.files(sprintf("BurdenIndividualPlatform/%s/%s_FAM/results/", analysis, toupper(disorder)),
                                full.names = T))
}else if(length(list.files(sprintf("BurdenIndividualPlatform/%s/%s/results/", analysis, toupper(disorder)),
                           full.names = T)) > 0){
  metadata <- fread(sprintf("BurdenIndividualPlatform/genebased_mega_analysis/%s/sampleinfo.%s.tsv", toupper(disorder), disorder), data.table = F)

  gene_sum_path <- list.files(sprintf("BurdenIndividualPlatform/%s/%s/results/", analysis, toupper(disorder)),
                                full.names = T)
}else{
  metadata_fam <- fread(sprintf("BurdenIndividualPlatform/genebased_mega_analysis/%s_FAM/sampleinfo.%s_fam.tsv", toupper(disorder), disorder), data.table = F)
  
  gene_sum_path <- list.files(sprintf("BurdenIndividualPlatform/%s/%s_FAM/results/", analysis, toupper(disorder)),
                                full.names = T)
}


genes <- fread(gene_path, data.table = F)
names(genes) <- c("chr", "start", "end", "transcript", "gene", "enzid")

#### meta analysis
final_stats <- data.frame()
for(cnv in c("DEL", "DUP")){
  for(chr in paste0("chr", c(1:22, "X"))){
    message(cnv, chr)
    
    if(length(gene_sum_path) > 1){ ### when there is multiple platform
      gene_sum_stats <- data.frame()
      gene_sum_weights <- list()
      gs_se_dat <- data.frame()
      
      for(path in gene_sum_path){
        
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
        
        
        if(file.exists(sprintf("BurdenIndividualPlatform/%s/%s/results/%s/loci.test12sample.%s.%s.%s.tsv",
                                analysis, toupper(disorder.here), platform.folder, disorder.here, cnv, chr))){
          sumstat <- fread(sprintf("BurdenIndividualPlatform/%s/%s/results/%s/loci.test12sample.%s.%s.%s.tsv",
                                   analysis, toupper(disorder.here), platform.folder, disorder.here, cnv, chr), data.table = F, 
                           colClasses = c("character", "character", "double", "double", "double", "double", "character", "character", "double", "character",
                                          "double", "double", "double", "double"))
          
          if(nrow(sumstat) == 0){
            sumstat <- data.frame("enzid", "chr", "start", "end", "gene.start", "gene.end", "gsymbol", "type", "coefficient", "pvalue", "sampleid")
            names(sumstat) <- sumstat[1, ]
            sumstat <- sumstat[-1, ]
            names(sumstat)[c(3:4, 9:11)] <- c("start_tmp", "end_tmp", paste0(c("coefficient_", "pvalue_"), platform), "samples")
          }else{
            sumstat$platform <- platform
            
            if(nrow(sumstat) == 1){
              sumstat$gsymbol[1] <- paste(sort(strsplit(sumstat$gsymbol[1], ",")[[1]]), collapse = ",")
              sumstat$enzid[1] <- paste(sort(strsplit(as.character(sumstat$enzid[1]), ",")[[1]]), collapse = ",")
              
            }else{
              sumstat$gsymbol <- sapply(lapply(sapply(sumstat$gsymbol, strsplit, ","), sort), paste, collapse = ",")
              sumstat$enzid <- sapply(lapply(sapply(as.character(sumstat$enzid), strsplit, ","), sort), paste, collapse = ",")
            }
            
            gs_se_dat <- rbind(gs_se_dat, sumstat[, c("enzid", "chr",  "gene.start", "gene.end", "gsymbol", "type", "platform",
                                                      "coefficient", "stderr")])
            
            sumstat <- sumstat[, c("enzid", "chr", "start", "end", "gene.start", "gene.end", "gsymbol", "type", "coefficient", "pvalue", "sampleid")]
            names(sumstat)[c(3:4, 9:11)] <- c("start_tmp", "end_tmp", paste0(c("coefficient_", "pvalue_"), platform), "samples")
            
          }
        }else{
          sumstat <- data.frame("enzid", "chr", "start", "end", "gene.start", "gene.end", "gsymbol", "type", "coefficient", "pvalue", "sampleid")
          names(sumstat) <- sumstat[1, ]
          sumstat <- sumstat[-1, ]
          names(sumstat)[c(3:4, 9:11)] <- c("start_tmp", "end_tmp", paste0(c("coefficient_", "pvalue_"), platform), "samples")
        }
        
        if((!file.exists(sprintf("BurdenIndividualPlatform/%s/%s/results/%s/loci.test12sample.%s.%s.%s.tsv",
                               analysis, toupper(disorder.here), platform.folder, disorder.here, cnv, chr))) | nrow(sumstat) == 0){
          if(ncol(gene_sum_stats) == 0){
            gene_sum_stats <- sumstat
            names(gene_sum_stats)[c(3:4, 11)] <- c("start", "end", "sampleid")
          }else{
            gene_sum_stats <- merge(gene_sum_stats, sumstat[, -c(3:4, 11)], by=c("enzid", "chr", "gene.start", "gene.end", "gsymbol", "type"), all = T)
          }
        }else{
          if(ncol(gene_sum_stats) == 0){
            gene_sum_stats <- sumstat
            names(gene_sum_stats)[c(3:4, 11)] <- c("start", "end", "sampleid")
          }else{
            gene_sum_stats <- merge(gene_sum_stats, sumstat, by=c("enzid", "chr", "gene.start", "gene.end", "gsymbol", "type"), all = T)
            ##### combine sample list from the two columns by ',', then strsplit them to get unique samples and collapse them back with ','
            if(nrow(gene_sum_stats) == 1){
              gene_sum_stats$sampleid <- paste(unique(strsplit(paste(gene_sum_stats$sampleid, gene_sum_stats$samples, sep=","), ",")[[1]]), collapse=",")
              gene_sum_stats$start <- min(gene_sum_stats$start, gene_sum_stats$start_tmp, na.rm = T)
              gene_sum_stats$end <- max(gene_sum_stats$end, gene_sum_stats$end_tmp, na.rm = T)
            }else{
              gene_sum_stats$sampleid <- sapply(lapply(sapply(paste(gene_sum_stats$sampleid, gene_sum_stats$samples, sep=","), strsplit, ","), unique), paste, collapse=",")
              gene_sum_stats$start <- pmin(gene_sum_stats$start, gene_sum_stats$start_tmp, na.rm = T)
              gene_sum_stats$end <- pmax(gene_sum_stats$end, gene_sum_stats$end_tmp, na.rm = T)
            }
            gene_sum_stats <- gene_sum_stats[, -which(names(gene_sum_stats) %in% c("start_tmp", "end_tmp", "samples"))]
          }
        }
        
        ncase <- sum(metadata.here$AFF == 1 & metadata.here$CNV_platform == platform.folder) 
        nctrl <- sum(metadata.here$AFF == 0 & metadata.here$CNV_platform == platform.folder)
        N = 4/((1/ncase)+(1/nctrl))
        w = sqrt(N)
        
        gene_sum_weights[[platform]] <- w
      }
      
      if(nrow(gene_sum_stats) > 0){
        
        
        
        z <- gene_sum_weights 
        for(platform in names(gene_sum_weights)){
          if(sum(!is.na(gene_sum_stats[, paste0("coefficient_", platform)])) > 0){
            z.tmp <- sign(gene_sum_stats[, paste0("coefficient_", platform)]) * -1 * qnorm(gene_sum_stats[, paste0("pvalue_", platform)]/2) * 
              gene_sum_weights[[platform]]
          }else{
            z.tmp <- NA
          }
          
          gene_sum_stats[, sprintf("z_%s", platform)] <- z.tmp
        }
        
        num = sapply(as.data.frame.matrix(t(gene_sum_stats[, paste0("z_", names(gene_sum_weights))])), sum, na.rm = T)
        
        message("get p")
        for(i in 1:length(num)){
          if(length(which(!is.na(gene_sum_stats[i, paste0("z_", names(gene_sum_weights))]))) > 0){
            den = sqrt(sum(sapply(gene_sum_weights[which(!is.na(gene_sum_stats[i, paste0("z_", names(gene_sum_weights))]))], "^", 2)))
            
          }else{
            den = sqrt(sum(sapply(gene_sum_weights, "^", 2)))
            
          }
          num[i] <- num[i]/den
        }
        
        gene_sum_stats$z = num
        gene_sum_stats$pvalue = 2*pnorm(pmin(gene_sum_stats$z,-gene_sum_stats$z))
        
        for(i in 1:nrow(gene_sum_stats)){
          tmp.dat <- gs_se_dat[gs_se_dat$type == gene_sum_stats$type[i] &
                                 gs_se_dat$enzid == gene_sum_stats$enzid[i] &
                                 gs_se_dat$chr == gene_sum_stats$chr[i] &
                                 gs_se_dat$gene.start == gene_sum_stats$gene.start[i] &
                                 gs_se_dat$gene.end == gene_sum_stats$gene.end[i], ]
          
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
            
            gene_sum_stats$se_test[i] <- tmp$zval
            gene_sum_stats$se_coeff[i] <- tmp$beta
            gene_sum_stats$se_ci_ub[i] <- tmp$ci.ub
            gene_sum_stats$se_ci_lb[i] <- tmp$ci.lb
          }else{
            
            gene_sum_stats$se_test[i] <- NA
            gene_sum_stats$se_coeff[i] <- NA
            gene_sum_stats$se_ci_ub[i] <- NA
            gene_sum_stats$se_ci_lb[i] <- NA
          }
          
        }
        
        final_stats <- rbind(final_stats, gene_sum_stats)
      }
    }else{ ### when there is only one platform
      platform.folder <- basename(gene_sum_path[1])
      platform <- gsub("\\.|\\-", "_", basename(gene_sum_path[1]))
      platform <- gsub("\\+", "Plus", platform)
      
      sumstat <- fread(sprintf("BurdenIndividualPlatform/%s/%s/results/%s/loci.test12sample.%s.%s.%s.%s.tsv",
                               analysis, toupper(disorder.here), platform.folder, disorder.here, cnv, chr, date), data.table = F)
      sumstat <- sumstat[, c("enzid", "chr", "start", "end", "gene.start", "gene.end", "gsymbol", "type", "coefficient", "pvalue", "sampleid",
                             "stderr", "testval")]
      names(sumstat)[9:11] <- c(paste0(c("coefficient_", "pvalue_"), platform), "sampleid")
      
      gene_sum_stats <- sumstat[, -c(12:13)]
      
      if(length(grep("FAM", gene_sum_path[1])) > 0){
        metadata.here <- metadata_fam
        disorder.here <- paste0(disorder, "_fam")
        platform <- paste0(platform, "_fam")
      }else{
        metadata.here <- metadata
        disorder.here <- disorder
      }
      
      ncase <- sum(metadata.here$AFF == 1 & metadata.here$CNV_platform == platform.folder) 
      nctrl <- sum(metadata.here$AFF == 0 & metadata.here$CNV_platform == platform.folder)
      N = 4/((1/ncase)+(1/nctrl))
      w = sqrt(N)
      
      if(sum(!is.na(gene_sum_stats[, paste0("coefficient_", platform)])) > 0){
        z.tmp <- sign(gene_sum_stats[, paste0("coefficient_", platform)]) * -1 * qnorm(gene_sum_stats[, paste0("pvalue_", platform)]/2) * w
      }else{
        z.tmp <- NA
      }
      
      den <- w
      gene_sum_stats$z = z.tmp/w
      gene_sum_stats$pvalue = 2*pnorm(pmin(gene_sum_stats$z,-gene_sum_stats$z))
      
      gs_sum_stats$se_test <- sumstat$testval
      gs_sum_stats$se_coeff <- sumstat$coefficient
      gs_sum_stats$se_ci_ub <- NA
      gs_sum_stats$se_ci_lb <- NA
      
      final_stats <- rbind(final_stats, gene_sum_stats)
    }
  }
}

write.table(final_stats, sprintf("MetaAnalysis/%s/%s/meta_analysis_all_tests.tsv", 
                                 analysis, toupper(disorder)), sep="\t", row.names=F, quote=F, col.names=T)
#### merge loci to be done to reduce multiple testing correction
test.out <- final_stats[order(final_stats$pvalue), ]
mergable <- T
while (mergable) {
  test.out <- test.out[order(test.out$pvalue), ]
  t.g <- GenomicRanges::GRanges(test.out$chr, IRanges::IRanges(as.numeric(as.character(test.out$start)), 
                                                               as.numeric(as.character(test.out$end))), "*")
  olap.t <- data.frame(IRanges::findOverlaps(t.g, t.g))
  olap.t <- olap.t[olap.t$queryHits != olap.t$subjectHits, ]
  olap.t <- olap.t[test.out$type[olap.t$queryHits] == test.out$type[olap.t$subjectHits], ]
  olap.t$pair <- paste(pmin(olap.t$queryHits, olap.t$subjectHits), 
                       pmax(olap.t$queryHits, olap.t$subjectHits), sep = ",")
  olap.t <- olap.t[!duplicated(olap.t$pair), ]
  mergable <- F
  if (nrow(olap.t) > 0) {
    message(sprintf("%s loci pair/s can potentially be merged", nrow(olap.t)))
    mergable <- T
    olap.t$pvalue <- pmin(test.out$pvalue[olap.t$queryHits], 
                          test.out$pvalue[olap.t$subjectHits])
    remove.test <- c()
    for (i in 1:nrow(olap.t)) {
      olap.rec <- olap.t[i, ]
      subject1 <- unique(strsplit(as.character(test.out$sampleid[olap.rec$queryHits]), 
                                  ",")[[1]])
      subject2 <- unique(strsplit(as.character(test.out$sampleid[olap.rec$subjectHits]), 
                                  ",")[[1]])
      if (!(olap.rec$queryHits %in% remove.test | olap.rec$queryHits %in% 
            remove.test)) {
        if (length(intersect(subject1, subject2))/length(union(subject1, 
                                                               subject2)) >= 0.8) {
          enzid <- paste(unique(strsplit(paste(test.out$enzid[olap.rec$queryHits], 
                                               test.out$enzid[olap.rec$subjectHits], sep = ","), 
                                         ",")[[1]]), collapse = ",")
          start <- min(test.out$start[olap.rec$queryHits], 
                       test.out$start[olap.rec$subjectHits])
          end <- max(test.out$end[olap.rec$queryHits], 
                     test.out$end[olap.rec$subjectHits])
          gene.start <- min(test.out$gene.start[olap.rec$queryHits], 
                            test.out$gene.start[olap.rec$subjectHits])
          gene.end <- max(test.out$gene.end[olap.rec$queryHits], 
                          test.out$gene.end[olap.rec$subjectHits])
          gsymbol <- paste(unique(strsplit(paste(test.out$gsymbol[olap.rec$queryHits], 
                                                 test.out$gsymbol[olap.rec$subjectHits], 
                                                 sep = ","), ",")[[1]]), collapse = ",")
          
          sampleid <- paste(unique(strsplit(paste(test.out$sampleid[olap.rec$queryHits], 
                                                  test.out$sampleid[olap.rec$subjectHits], 
                                                  sep = ","), ",")[[1]]), collapse = ",")
          
          if(sum(is.na(test.out$pvalue[c(olap.rec$queryHits, 
                               olap.rec$subjectHits)])) == 2){
            merge.test <- test.out[olap.rec$queryHits, ]
          }else{
            merge.test <- test.out[c(olap.rec$queryHits, 
                                     olap.rec$subjectHits)[which.min(test.out$pvalue[c(olap.rec$queryHits, 
                                                                                       olap.rec$subjectHits)])], ]
          }
          
          merge.test[, c("enzid", "gsymbol", "sampleid")] <- c(enzid, gsymbol, sampleid)
          merge.test[,  c("start", "end", "gene.start", "gene.end")] <- c(start, end, gene.start, gene.end)
          
          
          test.out <- rbind(test.out, 
                            merge.test)
          
          remove.test <- c(remove.test, olap.rec$queryHits, 
                           olap.rec$subjectHits)
        }
      }
    }
    remove.test <- unique(remove.test)
    if (length(remove.test) == 0) {
      mergable = F
    }
    else {
      test.out <- unique(test.out[-c(remove.test), ])
    }
  }
}

test.out$BHFDR <- NA
test.out$FWER <- NA
test.out$BHFDR[test.out$type == "DEL"] <- p.adjust(test.out$pvalue[test.out$type == "DEL"], method = "BH")
test.out$FWER[test.out$type == "DEL"] <- p.adjust(test.out$pvalue[test.out$type == "DEL"], method = "bonferroni")
test.out$BHFDR[test.out$type == "DUP"] <- p.adjust(test.out$pvalue[test.out$type == "DUP"], method = "BH")
test.out$FWER[test.out$type == "DUP"] <- p.adjust(test.out$pvalue[test.out$type == "DUP"], method = "bonferroni")


write.table(test.out, sprintf("MetaAnalysis/%s/%s/meta_analysis_reduced_tests.tsv", 
                              analysis, toupper(disorder)), sep="\t", row.names=F, quote=F, col.names=T)
