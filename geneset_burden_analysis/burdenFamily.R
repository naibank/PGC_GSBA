#########################################################
#### Install GSBurden ####
### only run a line below when doing this in Rstudio
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
# devtools::install_github("naibank/GSBurden")
library(GSBurden)
library(data.table)
library(survival)

CNVBurdenTest <- function (cnv.matrix, geneset, label, covariates, correctGlobalBurden = T, 
                           standardizeCoefficient = T, permutation = T, nperm = 100, 
                           BiasedUrn = F) 
{
  distinct.prefixes <- c("DEL", "DUP")
  model = "lm"
  if (length(unique(cnv.matrix[, label])) == 2) {
    message("dichotomous outcome variable detected. Logistic regression is being done ...")
    model = "glm"
  }
  else if (is.numeric(cnv.matrix[, label])) {
    message("continuous outcome variable detected. Linear regression is being done ...")
    model = "lm"
  }
  else if (is.factor(cnv.matrix[, label])) {
    message("ordinal outcome variable detected. Ordinal regression is being done ...")
    model = "clm"
  }
  else {
    stop("Non dichotomous or continuous outcome variable detected. The burden test cannot be run")
  }
  if (permutation) {
    ref.term <- sprintf("%s ~ %s", label, paste(covariates, 
                                                collapse = " + "))
    message(sprintf("Permuting sample labels for %s times", 
                    nperm))
    if (BiasedUrn) {
      lm.odds <- glm(ref.term, cnv.matrix, family = binomial(link = "logit"))
      d.odds <- exp(lm.odds$linear.predictors)
      n.case <- sum(cnv.matrix[, label])
      n.all <- length(cnv.matrix[, label])
      perm.hg <- BiasedUrn::rMFNCHypergeo(nran = nperm, 
                                          m = rep(1, n.all), n = n.case, odds = d.odds)
    }
    else {
      perm.hg <- data.frame(1:nrow(cnv.matrix), nrow(cnv.matrix):1)
      for (i in 1:nperm) {
        permuted <- sample(cnv.matrix[, label])
        perm.hg <- cbind(perm.hg, permuted)
      }
      perm.hg <- perm.hg[, -c(1:2)]
    }
  }
  perm.test.pvalues <- data.frame()
  test.out <- data.frame()
  for (cnvtype in distinct.prefixes) {
    for (this.gs in names(geneset)) {
      out.message <- sprintf("Testing %s", this.gs)
      if (length(distinct.prefixes) > 1) {
        out.message <- sprintf("%s for %s CNVs", out.message, 
                               cnvtype)
      }
      message(out.message)
      feature <- sprintf("%s_%s", this.gs, cnvtype)
      global <- sprintf("gene_count_%s", cnvtype)
      this.covariates <- covariates
      
      cnv.matrix$feature <- cnv.matrix[, feature]
      cnv.matrix$global <- cnv.matrix[, global] -  cnv.matrix[, feature]

      if(sd(cnv.matrix$feature) != 0)
        cnv.matrix$feature <- scale(cnv.matrix$feature)
      
      ref.model <- clogit(status ~ sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                            strata(FID) + global, cnv.matrix)
      add.model <- clogit(status ~ sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                            strata(FID) + global + feature, cnv.matrix)
      ano <- anova(ref.model, add.model, test = "Chisq")
      
      names(ano)[length(names(ano))] <- "pvalue"
      pvalue <- ano$pvalue[2]
      coefficient <- add.model$coefficients["feature"]
      sm <- summary(add.model)
      if ("feature" %in% rownames(sm$coefficients)) {
        waldp <- sm$coefficients["feature", 4]
        stderr <- sm$coefficients["feature", 2]
        testval <- sm$coefficients["feature", 3]
      }
      else {
        waldp <- stderr <- testval <- NA
      }
      
      conf <- tryCatch({
        confint(add.model)
      }, error = function(e) {
        return(NA)
      })
      if (is.null(nrow(conf))) {
        coeff.l <- 0
        coeff.u <- 0
      }
      else {
        coeff.l <- conf["feature", 1]
        coeff.u <- conf["feature", 2]
      }
      temp.out <- data.frame(geneset = this.gs, type = cnvtype, 
                             coefficient = coefficient, coeff.upper = coeff.u, 
                             coeff.lower = coeff.l, pvalue = pvalue, stderr, 
                             testval, waldp)
      test.out <- rbind(test.out, temp.out)
      if (permutation) {
        for (iperm in 1:nperm) {
          cnv.matrix$outcome.perm <- perm.hg[, iperm]
          ref.perm.term <- sprintf("outcome.perm ~ %s", 
                                   paste(this.covariates, collapse = " + "))
          add.perm.term <- sprintf("%s + %s", ref.perm.term, 
                                   feature)
          if (model == "lm") {
            ref.perm.model <- lm(ref.perm.term, cnv.matrix)
            add.perm.model <- lm(add.perm.term, cnv.matrix)
            ano.perm <- anova(ref.perm.model, add.perm.model, 
                              test = "Chisq")
          }
          else if (model == "glm") {
            ref.perm.model <- glm(ref.perm.term, cnv.matrix, 
                                  family = binomial(link = "logit"))
            add.perm.model <- glm(add.perm.term, cnv.matrix, 
                                  family = binomial(link = "logit"))
            ano.perm <- anova(ref.perm.model, add.perm.model, 
                              test = "Chisq")
          }
          else {
            ref.perm.model <- ordinal::clm(ref.perm.term, 
                                           data = cnv.matrix)
            add.perm.model <- ordinal::clm(add.perm.term, 
                                           data = cnv.matrix)
            ano.perm <- anova(ref.perm.model, add.perm.model)
          }
          names(ano.perm)[length(names(ano.perm))] <- "pvalue"
          coeff <- add.perm.model$coefficients[feature]
          perm.test.pvalues <- rbind(perm.test.pvalues, 
                                     data.frame(cnvtype = cnvtype, pvalue = ano.perm$pvalue[2], 
                                                coeff))
        }
      }
    }
  }
  if (permutation) {
    message("Calculating permutation-based FDR")
    test.out$permFDR <- 1
    for (i in 1:nrow(test.out)) {
      rec <- test.out[i, ]
      this.perm <- perm.test.pvalues[perm.test.pvalues$cnvtype == 
                                       rec$type, ]
      actual <- sum(test.out$pvalue[test.out$type == rec$type] <= 
                      rec$pvalue, na.rm = T)/nrow(na.omit(test.out[test.out$type == 
                                                                     rec$type, ]))
      perm <- sum(this.perm$pvalue <= rec$pvalue, na.rm = T)/nrow(na.omit(this.perm))
      fdr <- ifelse(perm/actual > 1, 1, perm/actual)
      test.out$permFDR[i] <- fdr
    }
  }
  test.out <- test.out[order(test.out$pvalue), ]
  for (i in 1:nrow(test.out)) {
    test.out$permFDR[i] <- min(test.out$permFDR[i:nrow(test.out)], 
                               na.rm = T)
  }
  list.out <- list(test.out, perm.test.pvalues)
  names(list.out) <- c("Test", "Permutation.Test")
  return(list.out)
}

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
# --set all --from 1 --end 12 --gs_path /hpf/largeprojects/tcagstor/users/worrawat/PGC/freeze_2024/BurdenAnalysis/geneset_celltype/gsCelltype.RData --disorder asd_fam
### params set, from, end
set <- params$set
from <- params$from
end <- params$end
gs_path <- params$gs_path
disorder <- params$disorder

#############################
#############################
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

cnv.matrix <- fread(sprintf("../%s.matrix.%s.tsv", disorder, set), data.table = F)

#nongs <- names(cnv.matrix)[1:(grep("DEL", names(cnv.matrix))[1]-1)]
#gs_del <- paste(c("cnv_count", "cnv_size", "gene_count", names(gsData.s)[from:end]), "DEL", sep="_")
#gs_dup <- paste(c("cnv_count", "cnv_size", "gene_count", names(gsData.s)[from:end]), "DUP", sep="_")

#count_combine <- cnv.matrix[, gs_del] + cnv.matrix[, gs_dup]

#cnv.matrix[, paste(c("cnv_count", "cnv_size", "gene_count", names(gsData.s)[from:end]), "All", sep="_")] <- count_combine

covariates <- c("PLATFORM", "sex", paste0("PC", c(1:10)))
if(length(unique(cnv.matrix$PLATFORM)) == 1){
  covariates <- c("sex", paste0("PC", c(1:10)))
}
#### perform global burden and gene set burden test ####
gsData.s <- gsData.s[from:end]
# gsData.s <- gsData.s[as.numeric(which(sapply(sapply(gsData.s, "%in%", "23191"), sum) > 0))]

if(length(gsData.s) > 0){
  for(platform in unique(cnv.matrix$PLATFORM)){
    dt <- CNVBurdenTest(cnv.matrix[which(cnv.matrix$PLATFORM == platform), ], gsData.s, "status", covariates, permutation = F)$Test
    
    if(!dir.exists(sprintf("../results/%s", 
                           platform))){
      dir.create(sprintf("../results/%s", 
                         platform), recursive = T)
    }
    
    #### write output to files ####
    
    write.table(dt, sprintf("../results/%s/%s.%s.%s.%s.burden.tsv", platform, disorder, set, from, end), sep="\t", row.names=F, quote=F, col.names=T)
  }
}