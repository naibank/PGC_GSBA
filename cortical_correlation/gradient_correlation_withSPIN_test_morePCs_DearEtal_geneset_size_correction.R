library(data.table)
library(ggplot2)
library(reshape2)
library(ggnewscale)
library(ggrepel)
library(scales)
library(tidyverse)
library(matrixStats)
library(ggpubr)
library(ggseg)
library(ggsegGlasser)
library(viridis)
# library(doSNOW)
# library(parallel)


rm(list=ls())
source("functions.R")

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
# --set all --from 1 --end 12 --gs_path /hpf/largeprojects/tcagstor/users/worrawat/PGC/freeze_2024/BurdenAnalysis/geneset_celltype/gsCelltype.RData --disorder asd
### params set, from, end
disorder <- params$disorder
type <- params$type
# time.start <- Sys.time()
# coreNumber <- 12
# cl <- makeCluster(coreNumber-1)
# registerDoSNOW(cl)

load("df_CorticalOrganizationHierarchy_Maps_Glasser180.RData")
gradient <- df_CorticalOrganizationHierarchy_Maps_Glasser180
gradient$region_id <- paste0("Brain_", 1:180)

newPCs <- fread("AHBA_PCs.csv", data.table = F)
newPCs$id <- paste("Brain", newPCs$id, sep="_")
newPCs <- newPCs[, -2]

gradient <- merge(gradient, newPCs, by.x = "region_id", by.y = "id", all = T)

dt.plot <- fread("all_test_results.tsv", data.table = F) #### summary statistics for brain regions
dt.plot <- dt.plot[dt.plot$ancestry == "ALL" & dt.plot$set == "all" & dt.plot$disorder != "XD", ]
dt.plot <- dt.plot[which(dt.plot$geneset_cat == "geneset_glasser"), ]

  # for(type in c("DEL", "DUP")){
dt <- dt.plot[dt.plot$disorder == disorder & dt.plot$type == type, ]
dt$geneset <- factor(dt$geneset, levels = paste0("Brain_", 1:180))
dt <- merge(dt, gradient, by.x = "geneset", by.y = "region_id", all.x = T)
dt <- dt[order(dt$geneset), ]

load("genesets/glasser.RData")
gsSize <- sapply(gsData.s, length)
gsSize <- data.frame("name" = names(gsData.s), "size" = gsSize)
dt <- merge(dt, gsSize, by.x = "geneset_lab", by.y = "name", all = T)

lm <- lm(z ~ size, dt)
dt$z <- lm$residuals

c1_cor <- cor.test(dt$z, dt$C1, method = "kendall")
# c2_cor <- cor.test(dt$z, dt$C2, method = "kendall")
# c3_cor <- cor.test(dt$z, dt$C3, method = "kendall")

c1_spin <- fGet_Corr_pval_SpinTest_two_maps_selectROIs_CorrArray(dt$z, dt$C1, 1:180, 10000, corr.type='kendall')[[2]]
# c2_spin <- fGet_Corr_pval_SpinTest_two_maps_selectROIs_CorrArray(dt$z, dt$C2, 1:180, 10000, corr.type='kendall')[[2]]
# c3_spin <- fGet_Corr_pval_SpinTest_two_maps_selectROIs_CorrArray(dt$z, dt$C3, 1:180, 10000, corr.type='kendall')[[2]]

c1_p <- c1_cor$p.value
# c2_p <- c2_cor$p.value
# c3_p <- c3_cor$p.value

c1_cor <- c1_cor$estimate
# c2_cor <- c2_cor$estimate
# c3_cor <- c3_cor$estimate

# color <- ifelse(type == "DEL", "red", "blue")
# 
# p1 <- ggplot(dt, aes(x = G1.fMRI, y = z)) + geom_point(color = color, size = 1) +
#   geom_smooth(method = "lm", color = "black") + theme_bw() + xlab("fMRI") +
#   ggtitle(sprintf("%s - \nR^2=%s, p=%s", paste(disorder, type), signif(f_MRI_cor, digits = 2), signif(f_MRI_spin, digits = 2)))
# 
# p2 <- ggplot(dt, aes(x = T1T2ratio, y = z)) + geom_point(color = color, size = 1) +
#     geom_smooth(method = "lm", color = "black") + theme_bw() + xlab("T1T2ratio") +
#   ggtitle(sprintf("%s - \nR^2=%s, p=%s", paste(disorder, type), signif(t1t2_cor, digits = 2), signif(t1t2_spin, digits = 2)))
#   
# p3 <- ggplot(dt, aes(x = PC1.AHBA, y = z)) + geom_point(color = color, size = 1) +
#   geom_smooth(method = "lm", color = "black") + theme_bw() + xlab("PC1.AHBA") +
#   ggtitle(sprintf("%s - \nR^2=%s, p=%s", paste(disorder, type), signif(trans_cor, digits = 2), signif(trans_spin, digits = 2)))
# 
# cowplot::plot_grid(p1, p2, p3, nrow = 3)
# ggsave(sprintf("%s_%s_spin.pdf", disorder, type), width = 3, height = 9)
# cor.test <- rbind(cor.test, data.frame(disorder, type, c1_cor, c1_p, c1_spin)) #, c2_cor, c2_p, c2_spin, c3_cor, c3_p, c3_spin))
tmp <- data.frame(disorder, type, c1_cor, c1_p, c1_spin)
write.table(tmp, paste0(disorder, "_", type, "_correlation_dear_size_corrected.tsv"), sep="\t", row.names=F, quote=F, col.names=T)
