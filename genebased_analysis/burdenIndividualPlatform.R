#########################################################
#### Install GSBurden ####
### only run a line below when doing this in Rstudio
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
# devtools::install_github("naibank/GSBurden")
library(GSBurden)
library(data.table)
library(parallel)
library(doSNOW)
library(survival)

#### read arguments from command line
args = commandArgs(trailingOnly=TRUE)
if(length(args) < 2){
  stop(call. = T, "Require more argument(s)")
}

paramNames <- grep("--", args)
paramValues <- paramNames + 1

if(length(grep("--", args[paramValues])) > 0){
  stop(call. = T, "Failed in reading arguments, '--' found in argument value")
}

if(length(args) != length(c(paramNames, paramValues))){
  message(sprintf("Ignore argument(s):%s", paste(args[-union(paramNames, paramValues)], collapse=",")))
}

params <- list()
paramNames <- gsub("--", "", args[paramNames])
for(i in 1:length(paramNames)){
  params[paramNames[i]] <- args[paramValues[i]]
}

#####
coreNumber <- 14
coreNumber <- ifelse(detectCores() < coreNumber, detectCores(), coreNumber)
cl <- makeCluster(coreNumber-1)
registerDoSNOW(cl)

# loci.test.out <- CNVLociTest(cnvs.in, cnv.matrix, genes, "status", covariates, permutation = F, nsubject = 10)

CNVLociTest <- function (cnv.table, cnv.matrix, annotation.table, label, covariates, 
                         geneset = list(), permutation = T, nperm = 1000, nsubject = 5, 
                         BiasedUrn = F, cnvtypes = "") {
  if(cnvtypes == "")
    cnvtypes <- unique(cnv.table$type)
  
  all.cnv.table <- cnv.table
  model = "lm"
  
  if (length(geneset) != 0) {
    annotation.table <- annotation.table[annotation.table$enzid %in% 
                                           unlist(geneset), ]
  }
  
  final.out <- data.frame()
  for (cnvtype in cnvtypes) {
    message(sprintf("Testing %s CNVs...", cnvtype))
    cnv.table <- all.cnv.table[all.cnv.table$type == cnvtype, ]
    cnv.g <- GenomicRanges::GRanges(cnv.table$chr, IRanges::IRanges(cnv.table$start, 
                                                                    cnv.table$end), "*")
    annotation.g <- GenomicRanges::GRanges(annotation.table$chr, 
                                           IRanges::IRanges(annotation.table$start, annotation.table$end), 
                                           "*")
    ###find CNVs overlap gene
    
    olap <- data.frame(IRanges::findOverlaps(cnv.g, annotation.g))
    olap$gsymbol <- annotation.table$gsymbol[olap$subjectHits]
    olap$chr <- annotation.table$chr[olap$subjectHits]
    olap <- unique(olap[, c("queryHits", "gsymbol", "chr")])
    
    ### collapse into single gene overlapped by multiple CNVs
    table <- aggregate(queryHits ~ gsymbol, olap, paste, 
                       collapse = ",")
    
    ### get unique loci overlapped by same set of CNVs (100% common)
    new.loci <- data.frame()
    for (i in unique(table$queryHits)) {
      temp.loci <- table[table$queryHits == i, ]
      this.samples <- unique(cnv.table$sample[as.numeric(as.character(strsplit(i, 
                                                                               ",")[[1]]))])
      temp.cnv <- cnv.table[unique(as.numeric(as.character(strsplit(i, 
                                                                    ",")[[1]]))), ]
      if (length(this.samples) >= nsubject) {
        sample <- paste(na.omit(this.samples), collapse = ",")
        temp.loci <- annotation.table[annotation.table$gsymbol %in% 
                                        temp.loci$gsymbol, ]
        temp.loci$gene.start <- min(temp.loci$start)
        temp.loci$gene.end <- max(temp.loci$end)
        temp.loci$enzid <- paste(unique(temp.loci$enzid), 
                                 collapse = ",")
        temp.loci$gsymbol <- paste(unique(temp.loci$gsymbol), 
                                   collapse = ",")
        temp.loci$start <- min(temp.cnv$start)
        temp.loci$end <- max(temp.cnv$end)
        temp.loci$sample <- sample
        temp.loci <- unique(temp.loci)
        new.loci <- rbind(new.loci, temp.loci)
      }
    }
    if (nrow(new.loci) == 0) {
      stop(sprintf("No loci having at least %s subjects. Please reduce nsubject param (default: nsubject=3)", 
                   nsubject))
    }
    current.count <- 1
    all.count <- nrow(new.loci)
    dt.out <- data.frame()
    dt.temp <- cnv.matrix[, c("sample", label, covariates)]
    
    for(covariate in covariates){
      if(length(unique(dt.temp[, covariate])) < 2){
        covariates <- covariates[covariates != covariate]
      }
    }
    
    ref.term <- sprintf("%s ~ %s", label, paste(covariates, 
                                                collapse = " + "))
    add.term <- sprintf("%s + %s", ref.term, "gene_count")
    sample.with.loci <- sapply(new.loci$sample, strsplit, ",")
    
    ref.model <- glm(ref.term, dt.temp, family = binomial(link = "logit"))
    
    if (permutation) {
      ref.perm.term <- sprintf("outcome.perm ~ %s", 
                               paste(covariates, collapse = " + "))
      ref.perm.models <- list()
      for(i in 1:nperm){
        dt.temp$outcome.perm <- perm.hg[, iperm]
        
        ref.perm.model <- glm(ref.perm.term, dt.temp,
                              family = binomial(link = "logit"))
        ref.perm.models[[length(ref.perm.models) + 1]] <- ref.perm.model
      }
    }
    
    for (iloci in 1:nrow(new.loci)) {
      temp.out <- data.frame()
      this.loci <- new.loci[iloci, ]
      
      dt.temp$gene_count <- 0
      # sample.with.loci <- strsplit(this.loci$sample, ",")[[1]]
      dt.temp$gene_count[which(dt.temp$sample %in% sample.with.loci[[iloci]])] <- 1
      
      add.model <- glm(add.term, dt.temp, family = binomial(link = "logit"))
      
      ano <- anova(ref.model, add.model, test = "Chisq")
      
      names(ano)[length(names(ano))] <- "pvalue"
      pvalue <- ano$pvalue[2]
      coefficient <- add.model$coefficients["gene_count"]
      
      sm <- summary(add.model)
      
      if("gene_count" %in% rownames(sm$coefficients)){
        waldp <- sm$coefficients["gene_count", 4]
        stderr <- sm$coefficients["gene_count", 2]
        testval <- sm$coefficients["gene_count", 3]
      }else{
        waldp <- stderr <- testval <- NA
      }
      
      temp.out <- data.frame(enzid = this.loci$enzid, chr = this.loci$chr, 
                             start = this.loci$start, end = this.loci$end, 
                             gsymbol = this.loci$gsymbol, type = cnvtype, 
                             coefficient = coefficient, pvalue = pvalue, stderr, testval, waldp, sampleid = this.loci$sample, 
                             gene.start = this.loci$gene.start, gene.end = this.loci$gene.end, 
                             sampleclass = paste(dt.temp$status[dt.temp$gene_count == 
                                                                  1], collapse = ","), stringsAsFactors = F)
      
      
      if (permutation) {
        ref.perm.term <- sprintf("outcome.perm ~ %s", 
                                 paste(covariates, collapse = " + "))
        add.perm.term <- sprintf("%s + %s", ref.perm.term, 
                                 "gene_count")
        
        perm.pvals <- foreach(iperm = 1:nperm, .combine=c) %dopar% {
          dt.temp$outcome.perm <- perm.hg[, iperm]
          
          ref.perm.model <- ref.perm.models[[iperm]]
          add.perm.model <- glm(add.perm.term, dt.temp,
                                family = binomial(link = "logit"))
          ano.perm <- anova(ref.perm.model, add.perm.model,
                            test = "Chisq")
          
          names(ano.perm)[length(names(ano.perm))] <- "pvalue"
          signif(ano.perm$pvalue[2], digits = 3)
        }
        
        # 
        temp.out[, sprintf("perm.pvalue.n%s", 1:nperm)] <- perm.pvals        # for (iperm in 1:nperm) {
        #   dt.temp$outcome.perm <- perm.hg[, iperm]
        #   
        #   ref.perm.model <- glm(ref.perm.term, dt.temp, 
        #                         family = binomial(link = "logit"))
        #   add.perm.model <- glm(add.perm.term, dt.temp, 
        #                         family = binomial(link = "logit"))
        #   ano.perm <- anova(ref.perm.model, add.perm.model, 
        #                     test = "Chisq")
        #   
        #   names(ano.perm)[length(names(ano.perm))] <- "pvalue"
        #   temp.out[, sprintf("perm.pvalue.n%s", iperm)] <- ano.perm$pvalue[2]
        # }
      }
      
      dt.out <- rbind(dt.out, temp.out)
      current.count <- current.count + 1
      if ((current.count%%10) == 0) {
        message(sprintf("%s Loci testing - %s out of %s tests were done %s", 
                        cnvtype, current.count, all.count, Sys.time()))
        flush.console()
      }
    }
    message("Loci testing done")
    # write.table(dt.out, sprintf("intermidiate.%s.dt.out.%s.tsv", cnvtype, Sys.Date()), sep="\t",row.names=F, quote=F, col.names=T)
    # message("Calculate FDR ...")
    dt.out <- dt.out[!is.na(dt.out$pvalue), ]
    # dt.out.merge <- mergeLoci(dt.out, "pvalue")
    # write.table(dt.out.merge, sprintf("intermidiate.%s.bd.merged.dt.out.%s.tsv", cnvtype, Sys.Date()), sep="\t",row.names=F, quote=F, col.names=T)
    
    final.out <- rbind(final.out, dt.out[, c("enzid", "chr", "start", "end", "gene.start", "gene.end", "gsymbol", "type", 
                                             "coefficient", "sampleid", "pvalue", "stderr", "testval", "waldp")])
    # if (permutation) {
    #   perm.pvalue <- c()
    #   for (i in 1:nperm) {
    #     message(sprintf("Merging permution #%s", i))
    #     perm.merge <- mergeLoci(dt.out, sprintf("perm.pvalue.n%s", 
    #                                             i))
    #     perm.pvalue <- c(perm.pvalue, perm.merge$pvalue)
    #   }
    #   dt.out.merge$permFDR <- 1
    #   for (i in 1:nrow(dt.out.merge)) {
    #     actual <- sum(dt.out.merge$pvalue <= dt.out.merge$pvalue[i])/nrow(dt.out.merge)
    #     perm <- sum(perm.pvalue <= dt.out.merge$pvalue[i])/length(perm.pvalue)
    #     dt.out.merge$permFDR[i] <- perm/actual
    #     dt.out.merge$permFDR[i] <- ifelse(dt.out.merge$permFDR[i] > 
    #                                         1, 1, dt.out.merge$permFDR[i])
    #   }
    #   rownames(dt.out.merge) <- NULL
    #   dt.out.merge <- dt.out.merge[, c(2:4, 9:10, 6, 7, 
    #                                    11:12, 1, 5, 8)]
    #   dt.out.merge <- dt.out.merge[order(dt.out.merge$pvalue), 
    #   ]
    #   for (i in 1:nrow(dt.out.merge)) {
    #     dt.out.merge$permFDR[i] <- min(dt.out.merge$permFDR[i:nrow(dt.out.merge)])
    #   }
    # }
    # final.out <- rbind(final.out, dt.out.merge)
  }
  #names(final.out)[2:3] <- c("loci.start", "loci.end")
  return(final.out)
}


gene.in <- read.delim("hg38_refGene_20200708.exon.txt", header = F, 
                      col.names = c("chr", "start", "end", "isoid", "genesymbol", "enzid"), stringsAsFactors = F)
genes <- GeneAnnotation(gene.in$enzid, gene.in$chr, gene.in$start, gene.in$end, gene.in$genesymbol)
genes <- genes[genes$chr == params$chr, ]

if(params$cnvtypes == "DEL"){
  genes <- genes[genes$gsymbol %in% readLines(sprintf("%s/intermediateFiles/Genes_with12DELsOrMore.txt",
                                                      toupper(params$disorder))), ]
}else{
  genes <- genes[genes$gsymbol %in% readLines(sprintf("%s/intermediateFiles/Genes_with12DUPsOrMore.txt",
                                                      toupper(params$disorder))), ]
}

#### read gene set ####
#### text file format (table with at least 2 columns; enzid and gene set name)
#### default separator is "\t" with header. If other separator is used, please define using parameter "sep".
load("gsMain_PGC_2021.RData")
gsMain <- gsMain[-grep("scRNA|CompleteLetha|GOP_less|SynKEGG|gnomAD", names(gsMain))]
# gs <- readGeneset("material/gs.table.tsv")
# outliers <- read.delim("QC_outlier_rare_cnvs.tsv", stringsAsFactors = F)
cnvs.in <- fread(params$cnv_path, data.table = F)
names(cnvs.in) <- c("FID", "sample", "chr", "start", "end", "type", "score", "site")
cnvs.in <- cnvs.in[which(cnvs.in$type %in% c(1, 3)), ]
cnvs.in$chr <- gsub("23", "X", cnvs.in$chr)
cnvs.in$chr <- gsub("24", "Y", cnvs.in$chr)
cnvs.in$chr <- paste0("chr", cnvs.in$chr)
cnvs.in$type <- ifelse(cnvs.in$type == 1, "DEL", "DUP")

cnvs.in <- cnvs.in[which(cnvs.in$type %in% params$cnvtypes & cnvs.in$chr %in% params$chr), ]

sampleinfo <- data.table::fread(params$sample_path, data.table = F)
names(sampleinfo) <- c("FID", "sample", "status", "sex", "ancestry", "PLATFORM", "dataset", "phenotype", paste0("PC", 1:10))

cnvs.in <- cnvs.in[cnvs.in$sample %in% sampleinfo$sample, ]

#### get gene sets matrix ####
#### transform CNV data into gene-set matrix
cnv.matrix <- getCNVGSMatrix(cnvs.in, genes, gsMain[1:2])

platform.case <- table(sampleinfo$PLATFORM[sampleinfo$status == 1])
platform.case <- names(which(platform.case >= 50))
platform.ctrl <- table(sampleinfo$PLATFORM[sampleinfo$status == 0])
platform.ctrl <- names(which(platform.ctrl >= 50))
platform.both <- intersect(platform.case, platform.ctrl)

sampleinfo <- sampleinfo[sampleinfo$PLATFORM %in% platform.both,]

cnvs.in <- cnvs.in[cnvs.in$sample %in% sampleinfo$sample, ]

### only have 62066 samples with cnvs > 5KB < 3MB
cnv.matrix <- merge(sampleinfo, cnv.matrix, by = "sample",
                    all.x = T)
cnv.matrix[is.na(cnv.matrix)] <- 0
covariates <- c("sex", paste0("PC", c(1:10)))

#### perform global burden and gene set burden test ####
# global.test.out <- CNVGlobalTest(cnv.matrix, "status", covariates, standardizeCoefficient = F)
# burden.test.out <- CNVBurdenTest(cnv.matrix, gsASD, "status", covariates, permutation = F)$Test
# write.table(global.test.out, "global.burden.tsv", sep="\t", row.names=F, quote=F, col.names=T)
# write.table(burden.test.out$Test, "gs.burden.tsv", sep="\t", row.names=F, quote=F, col.names=T)
for(platform in unique(cnv.matrix$PLATFORM)){
  loci.test.out <- CNVLociTest(cnvs.in, cnv.matrix[which(cnv.matrix$PLATFORM == platform), ], genes, "status", 
                               covariates, permutation = F, nsubject = 1, cnvtypes = params$cnvtypes)
  
  if(!dir.exists(sprintf("%s/results/%s", 
                         toupper(params$disorder), platform))){
    dir.create(sprintf("%s/results/%s", 
                       toupper(params$disorder), platform), recursive = T)
  }
  
  #### write output to files ####
  write.table(loci.test.out, sprintf("%s/results/%s/loci.test12sample.%s.%s.%s.tsv",
                                     toupper(params$disorder), platform, params$disorder, params$cnvtypes, params$chr), 
              sep="\t", row.names=F, quote=F, col.names=T)
}
