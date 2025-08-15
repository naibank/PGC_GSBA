library(psych)
library(ggplot2)
library(pheatmap)
library(corrplot)
library(reshape2)
library(data.table)
library(ggdendro)
library(ggrepel)
library(scales)
library(Hmisc)
library(corrplot)

geneset_labels <- c("AxonalTrans_A"="AxonalTrans",
                    "AxonalTrans_ExNeuPre"="AxonalTrans",
                    "AxonalTrans_ExNeuPst"="AxonalTrans",
                    "AxonalTrans_GPCPre"="AxonalTrans",
                    "AxonalTrans_MGPst"="AxonalTrans",
                    "AxonalTrans_OligoPst"="AxonalTrans",
                    "AxonalTrans_S"="AxonalTrans",
                    "AxonGuid_A"="AxonGuid",
                    "AxonGuid_ASTPre"="AxonGuid",
                    "AxonGuid_ASTPst"="AxonGuid",
                    "AxonGuid_ExNeuPre"="AxonGuid",
                    "AxonGuid_ExNeuPre_A"="AxonGuid",
                    "AxonGuid_ExNeuPst"="AxonGuid",
                    "AxonGuid_GPCPre"="AxonGuid",
                    "AxonGuid_InNeuPre"="AxonGuid",
                    "AxonGuid_InNeuPst"="AxonGuid",
                    "AxonGuid_MGPst"="AxonGuid",
                    "AxonGuid_MGPst_S"="AxonGuid",
                    "AxonGuid_OligoPst"="AxonGuid",
                    "AxonGuid_S"="AxonGuid",
                    "AxonGuid_VASCPst"="AxonGuid",
                    "AxonGuid_VASCPst_S"="AxonGuid",
                    "Calcium_A"="Calcium",
                    "Calcium_ASTPst"="Calcium",
                    "Calcium_ExNeuPst"="Calcium",
                    "Calcium_InNeuPst"="Calcium",
                    "Calcium_S"="Calcium",
                    "Calcium_VASCPst"="Calcium",
                    "CellCycle_A"="CellCycle",
                    "CellCycle_ExNeuPre"="CellCycle",
                    "CellCycle_ExNeuPst"="CellCycle",
                    "CellCycle_GPCPre"="CellCycle",
                    "CellCycle_GPCPre_S"="CellCycle",
                    "CellCycle_MGPst"="CellCycle",
                    "CellCycle_OligoPst"="CellCycle",
                    "CellCycle_S"="CellCycle",
                    "CellCycle_VASCPst"="CellCycle",
                    "CellSign_A"="CellSign",
                    "CellSign_ASTPre"="CellSign",
                    "CellSign_ASTPst"="CellSign",
                    "CellSign_ExNeuPre"="CellSign",
                    "CellSign_ExNeuPre_A"="CellSign",
                    "CellSign_ExNeuPst"="CellSign",
                    "CellSign_ExNeuPst_A"="CellSign",
                    "CellSign_ExNeuPst_S"="CellSign",
                    "CellSign_GPCPre"="CellSign",
                    "CellSign_GPCPre_S"="CellSign",
                    "CellSign_InNeuPre"="CellSign",
                    "CellSign_InNeuPst"="CellSign",
                    "CellSign_MGPst"="CellSign",
                    "CellSign_MGPst_A"="CellSign",
                    "CellSign_MGPst_S"="CellSign",
                    "CellSign_OligoPst"="CellSign",
                    "CellSign_OligoPst_S"="CellSign",
                    "CellSign_OPCPst"="CellSign",
                    "CellSign_S"="CellSign",
                    "CellSign_VASCPst"="CellSign",
                    "CellSign_VASCPst_A"="CellSign",
                    "CellSign_VASCPst_S"="CellSign",
                    "Chromatin_A"="Chromatin",
                    "Chromatin_ASTPst"="Chromatin",
                    "Chromatin_ExNeuPre"="Chromatin",
                    "Chromatin_ExNeuPst"="Chromatin",
                    "Chromatin_GPCPre"="Chromatin",
                    "Chromatin_InNeuPre"="Chromatin",
                    "Chromatin_MGPst"="Chromatin",
                    "Chromatin_OligoPst"="Chromatin",
                    "Chromatin_S"="Chromatin",
                    "Chromatin_VASCPst"="Chromatin",
                    "CytoMetabolism_A"="CytoMetabolism",
                    "CytoMetabolism_ASTPst"="CytoMetabolism",
                    "CytoMetabolism_ExNeuPre"="CytoMetabolism",
                    "CytoMetabolism_ExNeuPst"="CytoMetabolism",
                    "CytoMetabolism_GPCPre"="CytoMetabolism",
                    "CytoMetabolism_MGPst"="CytoMetabolism",
                    "CytoMetabolism_OligoPst"="CytoMetabolism",
                    "CytoMetabolism_S"="CytoMetabolism",
                    "CytoMetabolism_VASCPst"="CytoMetabolism",
                    "Endosome_A"="Endosome",
                    "Endosome_ASTPst"="Endosome",
                    "Endosome_ExNeuPre"="Endosome",
                    "Endosome_ExNeuPst"="Endosome",
                    "Endosome_MGPst"="Endosome",
                    "Endosome_MGPst_A"="Endosome",
                    "Endosome_MGPst_S"="Endosome",
                    "Endosome_OligoPst"="Endosome",
                    "Endosome_OligoPst_S"="Endosome",
                    "Endosome_S"="Endosome",
                    "Endosome_VASCPst"="Endosome",
                    "Endosome_VASCPst_S"="Endosome",
                    "GTPase_A"="GTPase",
                    "GTPase_ASTPst"="GTPase",
                    "GTPase_ExNeuPre"="GTPase",
                    "GTPase_ExNeuPre_A"="GTPase",
                    "GTPase_ExNeuPst"="GTPase",
                    "GTPase_ExNeuPst_A"="GTPase",
                    "GTPase_GPCPre"="GTPase",
                    "GTPase_InNeuPre"="GTPase",
                    "GTPase_InNeuPst"="GTPase",
                    "GTPase_MGPst"="GTPase",
                    "GTPase_MGPst_A"="GTPase",
                    "GTPase_MGPst_S"="GTPase",
                    "GTPase_OligoPst"="GTPase",
                    "GTPase_OligoPst_S"="GTPase",
                    "GTPase_S"="GTPase",
                    "GTPase_VASCPst"="GTPase",
                    "GTPase_VASCPst_A"="GTPase",
                    "GTPase_VASCPst_S"="GTPase",
                    "Hormone_A"="Hormone",
                    "Hormone_GPCPre"="Hormone",
                    "Hormone_MGPst"="Hormone",
                    "Hormone_S"="Hormone",
                    "Hormone_VASCPst"="Hormone",
                    "MAPK_A"="MAPK",
                    "MAPK_ASTPre"="MAPK",
                    "MAPK_ASTPst"="MAPK",
                    "MAPK_ASTPst_A"="MAPK",
                    "MAPK_ASTPst_S"="MAPK",
                    "MAPK_ExNeuPre"="MAPK",
                    "MAPK_ExNeuPre_A"="MAPK",
                    "MAPK_ExNeuPre_S"="MAPK",
                    "MAPK_ExNeuPst"="MAPK",
                    "MAPK_ExNeuPst_A"="MAPK",
                    "MAPK_ExNeuPst_S"="MAPK",
                    "MAPK_GPCPre"="MAPK",
                    "MAPK_GPCPre_A"="MAPK",
                    "MAPK_GPCPre_S"="MAPK",
                    "MAPK_InNeuPre"="MAPK",
                    "MAPK_InNeuPst"="MAPK",
                    "MAPK_InNeuPst_S"="MAPK",
                    "MAPK_MGPst"="MAPK",
                    "MAPK_MGPst_A"="MAPK",
                    "MAPK_MGPst_S"="MAPK",
                    "MAPK_OligoPst"="MAPK",
                    "MAPK_OligoPst_S"="MAPK",
                    "MAPK_OPCPst"="MAPK",
                    "MAPK_S"="MAPK",
                    "MAPK_VASCPst"="MAPK",
                    "MAPK_VASCPst_A"="MAPK",
                    "MAPK_VASCPst_S"="MAPK",
                    "MitoMetabolism_A"="MitoMetabolism",
                    "MitoMetabolism_ASTPst"="MitoMetabolism",
                    "MitoMetabolism_ExNeuPre"="MitoMetabolism",
                    "MitoMetabolism_ExNeuPst"="MitoMetabolism",
                    "MitoMetabolism_MGPst"="MitoMetabolism",
                    "MitoMetabolism_S"="MitoMetabolism",
                    "MitoMetabolism_VASCPst"="MitoMetabolism",
                    "MitoSpindle_A"="MitoSpindle",
                    "MitoSpindle_ASTPst"="MitoSpindle",
                    "MitoSpindle_ExNeuPre"="MitoSpindle",
                    "MitoSpindle_ExNeuPst"="MitoSpindle",
                    "MitoSpindle_ExNeuPst_A"="MitoSpindle",
                    "MitoSpindle_GPCPre"="MitoSpindle",
                    "MitoSpindle_GPCPre_S"="MitoSpindle",
                    "MitoSpindle_InNeuPre"="MitoSpindle",
                    "MitoSpindle_InNeuPst"="MitoSpindle",
                    "MitoSpindle_MGPst"="MitoSpindle",
                    "MitoSpindle_OligoPst"="MitoSpindle",
                    "MitoSpindle_OligoPst_S"="MitoSpindle",
                    "MitoSpindle_S"="MitoSpindle",
                    "MitoSpindle_VASCPst"="MitoSpindle",
                    "NervousSys_A"="NervousSys",
                    "NervousSys_ASTPre"="NervousSys",
                    "NervousSys_ASTPst"="NervousSys",
                    "NervousSys_ASTPst_A"="NervousSys",
                    "NervousSys_ExNeuPre"="NervousSys",
                    "NervousSys_ExNeuPre_A"="NervousSys",
                    "NervousSys_ExNeuPre_S"="NervousSys",
                    "NervousSys_ExNeuPst"="NervousSys",
                    "NervousSys_ExNeuPst_A"="NervousSys",
                    "NervousSys_ExNeuPst_S"="NervousSys",
                    "NervousSys_GPCPre"="NervousSys",
                    "NervousSys_InNeuPre"="NervousSys",
                    "NervousSys_InNeuPst"="NervousSys",
                    "NervousSys_InNeuPst_A"="NervousSys",
                    "NervousSys_InNeuPst_S"="NervousSys",
                    "NervousSys_MGPst"="NervousSys",
                    "NervousSys_MGPst_A"="NervousSys",
                    "NervousSys_MGPst_S"="NervousSys",
                    "NervousSys_OligoPst"="NervousSys",
                    "NervousSys_OligoPst_S"="NervousSys",
                    "NervousSys_OPCPre"="NervousSys",
                    "NervousSys_OPCPst"="NervousSys",
                    "NervousSys_S"="NervousSys",
                    "NervousSys_VASCPst"="NervousSys",
                    "NervousSys_VASCPst_A"="NervousSys",
                    "NervousSys_VASCPst_S"="NervousSys",
                    "SynapPlasticity_A"="SynapPlasticity",
                    "SynapPlasticity_ASTPre"="SynapPlasticity",
                    "SynapPlasticity_ASTPst"="SynapPlasticity",
                    "SynapPlasticity_ExNeuPre"="SynapPlasticity",
                    "SynapPlasticity_ExNeuPre_A"="SynapPlasticity",
                    "SynapPlasticity_ExNeuPst"="SynapPlasticity",
                    "SynapPlasticity_ExNeuPst_A"="SynapPlasticity",
                    "SynapPlasticity_ExNeuPst_S"="SynapPlasticity",
                    "SynapPlasticity_GPCPre"="SynapPlasticity",
                    "SynapPlasticity_InNeuPst"="SynapPlasticity",
                    "SynapPlasticity_MGPst"="SynapPlasticity",
                    "SynapPlasticity_OligoPst"="SynapPlasticity",
                    "SynapPlasticity_S"="SynapPlasticity",
                    "SynapPlasticity_VASCPst"="SynapPlasticity",
                    "SynapPlasticity_VASCPst_S"="SynapPlasticity",
                    "SynapTransmiss_A"="SynapTransmiss",
                    "SynapTransmiss_ASTPst"="SynapTransmiss",
                    "SynapTransmiss_ASTPst_A"="SynapTransmiss",
                    "SynapTransmiss_ExNeuPre"="SynapTransmiss",
                    "SynapTransmiss_ExNeuPre_A"="SynapTransmiss",
                    "SynapTransmiss_ExNeuPst"="SynapTransmiss",
                    "SynapTransmiss_ExNeuPst_A"="SynapTransmiss",
                    "SynapTransmiss_ExNeuPst_S"="SynapTransmiss",
                    "SynapTransmiss_GPCPre"="SynapTransmiss",
                    "SynapTransmiss_InNeuPre"="SynapTransmiss",
                    "SynapTransmiss_InNeuPst"="SynapTransmiss",
                    "SynapTransmiss_InNeuPst_A"="SynapTransmiss",
                    "SynapTransmiss_InNeuPst_S"="SynapTransmiss",
                    "SynapTransmiss_MGPst"="SynapTransmiss",
                    "SynapTransmiss_OligoPst"="SynapTransmiss",
                    "SynapTransmiss_OligoPst_S"="SynapTransmiss",
                    "SynapTransmiss_OPCPst"="SynapTransmiss",
                    "SynapTransmiss_S"="SynapTransmiss",
                    "SynapTransmiss_VASCPst"="SynapTransmiss",
                    "SynapTransmiss_VASCPst_S"="SynapTransmiss",
                    "Translation_A"="Translation",
                    "Translation_ASTPre"="Translation",
                    "Translation_ASTPst"="Translation",
                    "Translation_ExNeuPre"="Translation",
                    "Translation_ExNeuPre_A"="Translation",
                    "Translation_ExNeuPre_S"="Translation",
                    "Translation_ExNeuPst"="Translation",
                    "Translation_GPCPre"="Translation",
                    "Translation_GPCPre_A"="Translation",
                    "Translation_GPCPre_S"="Translation",
                    "Translation_InNeuPst"="Translation",
                    "Translation_MGPst"="Translation",
                    "Translation_OligoPst"="Translation",
                    "Translation_S"="Translation",
                    "Translation_VASCPst"="Translation",
                    "Ubiquitin_A"="Ubiquitin",
                    "Ubiquitin_ASTPst"="Ubiquitin",
                    "Ubiquitin_ExNeuPre"="Ubiquitin",
                    "Ubiquitin_ExNeuPst"="Ubiquitin",
                    "Ubiquitin_GPCPre"="Ubiquitin",
                    "Ubiquitin_InNeuPre"="Ubiquitin",
                    "Ubiquitin_InNeuPst"="Ubiquitin",
                    "Ubiquitin_MGPst"="Ubiquitin",
                    "Ubiquitin_MGPst_S"="Ubiquitin",
                    "Ubiquitin_OligoPst"="Ubiquitin",
                    "Ubiquitin_S"="Ubiquitin",
                    "Ubiquitin_VASCPst"="Ubiquitin",
                    "VesicleTraffick_A"="VesicleTraffick",
                    "VesicleTraffick_ASTPst"="VesicleTraffick",
                    "VesicleTraffick_ASTPst_A"="VesicleTraffick",
                    "VesicleTraffick_ExNeuPre"="VesicleTraffick",
                    "VesicleTraffick_ExNeuPre_A"="VesicleTraffick",
                    "VesicleTraffick_ExNeuPst"="VesicleTraffick",
                    "VesicleTraffick_ExNeuPst_A"="VesicleTraffick",
                    "VesicleTraffick_ExNeuPst_S"="VesicleTraffick",
                    "VesicleTraffick_GPCPre"="VesicleTraffick",
                    "VesicleTraffick_InNeuPre"="VesicleTraffick",
                    "VesicleTraffick_InNeuPst"="VesicleTraffick",
                    "VesicleTraffick_MGPst"="VesicleTraffick",
                    "VesicleTraffick_MGPst_S"="VesicleTraffick",
                    "VesicleTraffick_OligoPst"="VesicleTraffick",
                    "VesicleTraffick_OPCPst"="VesicleTraffick",
                    "VesicleTraffick_S"="VesicleTraffick",
                    "VesicleTraffick_VASCPst"="VesicleTraffick",
                    "VesicleTraffick_VASCPst_S"="VesicleTraffick")

dt <- fread("association_results_2_3_ways_withcontributegenes_atp0.05.tsv", data.table = F)
dt <- dt[!dt$stratification %in% "celltype-sa-axis" & dt$set == "all", ]

#### write ylab gene set name for pathway loadings
label <- c()
for(i in 1:nrow(dt)){
  geneset <- dt$geneset[i]
  pathway <- dt$pathway[i]
  label <- c(label, paste0("\"", geneset, "\"", "=", "\"", pathway, "\","))
}
label <- unique(label)
writeLines(label, "ylab_label_all.txt")

stratification <- "pathway-celltype_pathway-sa-axis"
if(stratification == "all"){
  dt.plot <- dt
}else if(stratification == "pathway-celltype_pathway-sa-axis"){
  dt.plot <- dt[dt$stratification %in% c("pathway-sa-axis", "pathway-celltype"),]
}else{
  dt.plot <- dt[dt$stratification == stratification, ]
}

dt.plot$qvalue[dt.plot$type == "DEL"] <- p.adjust(dt.plot$pvalue[dt.plot$type == "DEL"], method = "BH")
dt.plot$qvalue[dt.plot$type == "DUP"] <- p.adjust(dt.plot$pvalue[dt.plot$type == "DUP"], method = "BH")

dt.plot$fwer[dt.plot$type == "DEL"] <- p.adjust(dt.plot$pvalue[dt.plot$type == "DEL"], method = "bonferroni")
dt.plot$fwer[dt.plot$type == "DUP"] <- p.adjust(dt.plot$pvalue[dt.plot$type == "DUP"], method = "bonferroni")

asd.del <- dt.plot[dt.plot$disorder == "ASD" & dt.plot$type == "DEL", ]
asd.dup <- dt.plot[dt.plot$disorder == "ASD" & dt.plot$type == "DUP", ]

scz.del <- dt.plot[dt.plot$disorder == "SCZ" & dt.plot$type == "DEL", ]
scz.dup <- dt.plot[dt.plot$disorder == "SCZ" & dt.plot$type == "DUP", ]

bd.del <- dt.plot[dt.plot$disorder == "BD" & dt.plot$type == "DEL", ]
bd.dup <- dt.plot[dt.plot$disorder == "BD" & dt.plot$type == "DUP", ]

adhd.del <- dt.plot[dt.plot$disorder == "ADHD" & dt.plot$type == "DEL", ]
adhd.dup <- dt.plot[dt.plot$disorder == "ADHD" & dt.plot$type == "DUP", ]

mdd.del <- dt.plot[dt.plot$disorder == "MDD" & dt.plot$type == "DEL", ]
mdd.dup <- dt.plot[dt.plot$disorder == "MDD" & dt.plot$type == "DUP", ]

ptsd.del <- dt.plot[dt.plot$disorder == "PTSD" & dt.plot$type == "DEL", ]
ptsd.dup <- dt.plot[dt.plot$disorder == "PTSD" & dt.plot$type == "DUP", ]

zmat <- data.frame("ASD_DEL" = asd.del$z,
                   "ASD_DUP" = asd.dup$z,
                   "SCZ_DEL" = scz.del$z,
                   "SCZ_DUP" = scz.dup$z,
                   "BD_DEL" = bd.del$z,
                   "BD_DUP" = bd.dup$z,
                   "ADHD_DEL" = adhd.del$z,
                   "ADHD_DUP" = adhd.dup$z,
                   "MDD_DEL" = mdd.del$z,
                   "MDD_DUP" = mdd.dup$z,
                   "PTSD_DEL" = ptsd.del$z,
                   "PTSD_DUP" = ptsd.dup$z)

rownames(zmat) <- asd.del$geneset
 
zmat_deldup <- data.frame("ASD" = c(asd.del$z, asd.dup$z),
                   "SCZ" = c(scz.del$z, scz.dup$z),
                   "BD" = c(bd.del$z, bd.dup$z),
                   "ADHD" = c(adhd.del$z, adhd.dup$z),
                   "MDD" = c(mdd.del$z, mdd.dup$z),
                   "PTSD" = c(ptsd.del$z, ptsd.dup$z))
rownames(zmat_deldup) <- c(paste(asd.del$geneset, "DEL", sep="_"), paste(asd.del$geneset, "DUP", sep="_"))
diag.clust <- hclust(dist(t(zmat_deldup)))

ggdendrogram(diag.clust, rotate = T, leaf_labels = F, label = T) + coord_flip() + scale_x_reverse()
ggsave("ggdendro.diag.all_factor.pdf", width = 3, height = 4.5)

# Standardize Z-scores by row (pathways)
zmat_scaled <- scale(zmat, center = TRUE, scale = TRUE)
pca_result <- prcomp(zmat_scaled)
pdf(sprintf("pca_result_%s", stratification))
plot(pca_result, type = "lines")  # Scree plot
dev.off()
s_pca <- summary(pca_result)

find_elbow_threshold <- function(eig_vals, drop_thresh = 0.5) {
  diffs <- -diff(eig_vals)  # drop in eigenvalues
  elbow <- which(diffs < drop_thresh)[1] + 1
  if (is.na(elbow)) elbow <- length(eig_vals)  # fallback
  return(elbow)
}

# Apply to eigenvalues
eig_vals <- (pca_result$sdev)^2
num_factors <- find_elbow_threshold(eig_vals)

print(paste(stratification, num_factors))
# print(min(which(s_pca$importance["Cumulative Proportion",] >= 0.5)))
# Compute correlation matrix across disorder-CNVs
# cor_mat <- cor(zmat_scaled, method = "kendall")
# 
# # Visualize correlation matrix
# cor_mat_plot <- cor_mat
# cor_mat_plot[cor_mat_plot > 0.3] <- 0.3
# cor_mat_plot[cor_mat_plot < -0.3] <- -0.3
# pheatmap(cor_mat_plot, main = paste("Kendall's Tau Correlation Across Disorder-CNVs\n",
#                                     gsub("Sa", "SA", tools::toTitleCase(stratification)), "Stratification"), 
#          filename = sprintf("%s.pdf", gsub("\\-", "_", stratification)), width = 5.5, height = 5.5)
# pheatmap(cor_mat_plot)
# dev.off()
# Compute Kendall's correlation and p-values manually
kendall_test <- function(x, y) {
  res <- cor.test(x, y, method = "kendall")
  return(c(estimate = res$estimate, p.value = res$p.value))
}

# Compute correlation matrix and p-values
n <- ncol(zmat_scaled)
cor_mat <- matrix(NA, n, n)
p_mat <- matrix(NA, n, n)
colnames(cor_mat) <- rownames(cor_mat) <- colnames(zmat_scaled)
colnames(p_mat) <- rownames(p_mat) <- colnames(zmat_scaled)

for (i in 1:n) {
  for (j in 1:n) {
    out <- kendall_test(zmat_scaled[, i], zmat_scaled[, j])
    cor_mat[i, j] <- out["estimate.tau"]
    p_mat[i, j] <- out["p.value"]
    if(i == j)
      p_mat[i, j] <- 1
  }
}

write.table(cor_mat, sprintf("cor_matrix_%s.tsv", gsub("\\-", "_", stratification)), sep="\t", row.names=F, quote=F, col.names=T)

# Cap correlations for visualization
cor_mat_plot <- cor_mat
cor_mat_plot[cor_mat_plot > 0.3] <- 0.3
cor_mat_plot[cor_mat_plot < -0.3] <- -0.3

# Add significance stars
signif_matrix <- matrix("", n, n)
signif_matrix[p_mat < 0.05] <- "*"
signif_matrix[p_mat*132 < 0.05] <- "**"

# Plot heatmap with annotations
pheatmap(cor_mat_plot,
         display_numbers = signif_matrix, fontsize_number = 14,
         main = "Correlation of disorders based on gene-set effect sizes",
         filename = sprintf("%s.pdf", gsub("\\-", "_", stratification)),
         width = 5.5,
         height = 5.5)

# Also show it in R session
pheatmap(cor_mat_plot)
dev.off()
# tmp <- fa.parallel(cor_mat, fa = "fa", n.obs = nrow(asd.del))
# vss_res <- vss(cor_mat, n = 10, n.obs = nrow(asd.del))

# if(tmp$nfact <= 1){
#   tmp <- fa.parallel(cor_mat, fa = "fa")
# }
# 
# if(tmp$nfact <= 1){
#   tmp$nfact <- 2
# }
# 
# if(stratification == "celltype-sa-axis"){
#   tmp$nfact <- which.min(vss_res$vss.stats$SABIC)
# }
# 
# tmp$nfact <- num_factors

efa_result <- fa(cor_mat, nfactors = num_factors, rotate = "oblimin", fm = "ml")
loadings_matrix <- as.matrix(efa_result$loadings[, 1:num_factors])
loadings_matrix <- loadings_matrix[, c("ML1", "ML2", "ML3")]

sem_df <- as.data.table(loadings_matrix)
sem_df$Trait <- rownames(loadings_matrix)
sem_long <- melt(sem_df, id.vars = "Trait", variable.name = "Factor", value.name = "Loading")
sem_long$Factor <- gsub("ML", "F", sem_long$Factor)
write.table(sem_long, sprintf("sem_loading_%s.tsv", gsub("\\-", "_", stratification)), sep="\t", row.names=F, quote=F, col.names=T)

diags <- rev(diag.clust$labels[diag.clust$order])
rep.diags <- c()
for(diag in diags){
  rep.diags <- c(rep.diags, rep(diag, 2))
}      
# Create SEM-style plot of disorder-CNVs and factors
sem_long$Trait <- factor(sem_long$Trait, levels = paste(rep.diags, rep(c("DEL", "DUP"), 6), sep="_"))
ggplot(sem_long, aes(x = Factor, y = Trait, fill = Loading)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() + ylab("Diagnosis-CNV") + theme(legend.position = "right", title = element_text(face = "bold")) +
  ggtitle("Factor loadings of disorders")
ggsave(sprintf("sem_loading_%s.pdf", gsub("\\-", "_", stratification)), width = num_factors*1.15, height = 5)

# Correlate each pathway with factor loadings
pathway_factor_corr <- zmat_scaled %*% loadings_matrix
write.csv(pathway_factor_corr, sprintf("pathway_factor_correlations_%s.csv",  gsub("\\-", "_", stratification)))

loadings <- as.data.frame(pathway_factor_corr)
loadings$color <- "red"
loadings$V1 <- rownames(loadings)
loadings$pathway <- sapply(sapply(loadings$V1, strsplit, "_"), "[", 1)
loadings$celltype <- sapply(sapply(loadings$V1, strsplit, "_"), "[", 2)
loadings$brain <- sapply(sapply(loadings$V1, strsplit, "_"), "[", 3)

# for(i in 1:nrow(loadings)){
#   if(loadings$celltype[i] %in% c("A", "S")){
#     loadings$brain[i] <- loadings$celltype[i]
#     loadings$celltype[i] <- NA
#   }
#   # if(is.na(loadings$brain[i])){
#   #   loadings$brain[i] <- ""
#   # }
# }
#### cluster all factor together
# tmp_melt <- reshape::melt(loadings[, c("ML1", "ML2", "ML3", "pathway", "celltype")], id.vars = c("pathway", "celltype"))
# tmp_cast <- reshape::cast(tmp_melt, celltype ~ pathway+variable, mean, na.rm=T)
# tmp_cast[is.na(tmp_cast)] <- 0
# 
# gs.clust <- hclust(dist(tmp_cast))
# 
# ggdendrogram(gs.clust, rotate = F, leaf_labels = F, label = T)
# ggsave(sprintf("ggdendro.%s.all_factor.pdf", stratification), width = 8, height = 4)

# celltype.out$geneset_lab <- factor(celltype.out$geneset_lab, levels=gs.clust$labels[gs.clust$order])
#### cluster all factor together end here
sign_bicluster_order <- function(mat) {
  stopifnot(is.matrix(mat))
  
  # Column grouping
  col_means <- colMeans(mat, na.rm = TRUE)
  col_A <- which(col_means > 0)
  col_B <- which(col_means < 0)
  
  # Row grouping
  row_means <- rowMeans(mat, na.rm = TRUE)
  row_A <- which(row_means > 0)
  row_B <- which(row_means < 0)
  
  # Column sorting
  if (length(row_A) > 0) {
    col_A_sorted <- col_A[order(colMeans(mat[row_A, col_A, drop = FALSE]), decreasing = TRUE)]
  } else {
    col_A_sorted <- integer(0)
  }
  
  if (length(row_B) > 0) {
    col_B_sorted <- col_B[order(colMeans(mat[row_B, col_B, drop = FALSE]), decreasing = TRUE)]
  } else {
    col_B_sorted <- integer(0)
  }
  
  # Row sorting
  if (length(col_A) > 0) {
    row_A_sorted <- row_A[order(rowMeans(mat[row_A, col_A, drop = FALSE]), decreasing = TRUE)]
  } else {
    row_A_sorted <- integer(0)
  }
  
  if (length(col_B) > 0) {
    row_B_sorted <- row_B[order(rowMeans(mat[row_B, col_B, drop = FALSE]), decreasing = TRUE)]
  } else {
    row_B_sorted <- integer(0)
  }
  
  # Final order
  new_row_order <- c(row_A_sorted, row_B_sorted)
  new_col_order <- c(col_A_sorted, col_B_sorted)
  
  mat_reordered <- mat[new_row_order, new_col_order, drop = FALSE]
  return(list(matrix = mat_reordered, row_order = new_row_order, col_order = new_col_order))
}


p <- list()
# p_nolabel <- list()
for(i in 1:num_factors){
  cols <- c("pathway", "celltype", "brain", names(loadings)[i], "color")
  cols <- cols[cols %in% names(loadings)]
  
  f <- loadings[, cols]    
  names(f)[4] <- "loading"
  
  f$color[f$loading < 0] <- "#3366FF"
  
  min <- -max(abs(f$loading))
  max <- max(abs(f$loading))

  if(max > 4){
    breaks = c(-4, 0, 4)
  }else{
    breaks = c(-2, 0 , 2)
  }
  # min <- -3
  # max <- 3
  # breaks <- c(-2, 0, 2)
  
  color_scale <- c("#053061", "#2166AC", "#4393C3", "#72A5CE", "#D1E5F0", "#FFFFFF", "#FFF0A0", "#FDAE61", "#F47D43", "#D75027", "#A50026")
  color_scale <- scale_fill_gradientn(colors = color_scale, limits = c(min, max), name = sprintf("Factor %s\nscore",i))
  
  f <- f[order(abs(f$loading), decreasing = T), ]
  # sel_pathway <- unique(f$pathway[1:round(nrow(f)/10)])
  # sel_celltype <- unique(f$celltype[1:round(nrow(f)/10)])
  # # f <- f[f$pathway %in% sel_pathway & f$celltype %in% sel_celltype, ]
  
  tmp_loadings <- f #[f$brain == tmp_brain, ]
  tmp_loadings <- tmp_loadings[order(tmp_loadings$loading, decreasing = T), ]
  
  for(celltype in unique(tmp_loadings$celltype)){
    for(pathway in unique(tmp_loadings$pathway)){
      if(sum(tmp_loadings$celltype == celltype & tmp_loadings$pathway == pathway) == 0){
        tmp_rec <- tmp_loadings[1, ]
        tmp_rec$pathway <- pathway
        tmp_rec$celltype <- celltype
        tmp_rec$loading <- 0
        
        tmp_loadings <- rbind(tmp_loadings, tmp_rec)
      }
    }
  }
  
  tmp_loadings$celltype <- gsub("VASC", "Vasc", tmp_loadings$celltype)
  tmp_loadings$celltype <- gsub("AST", "Ast", tmp_loadings$celltype)
  tmp_loadings$celltype <- gsub("Pst", "Post", tmp_loadings$celltype)
  tmp_loadings$celltype <- gsub("GPC", "Gpc", tmp_loadings$celltype)
  tmp_loadings$celltype <- gsub("OPC", "Opc", tmp_loadings$celltype)
  tmp_loadings$celltype <- gsub("MG", "Mg", tmp_loadings$celltype)
  tmp_loadings$pathway <- gsub("CytoMetabolism", "CytoMetab", tmp_loadings$pathway)
  tmp_loadings$pathway <- gsub("MitoMetabolism", "MitoMetab", tmp_loadings$pathway)
  tmp_loadings$pathway <- gsub("MitoSpindle", "MitotSpind", tmp_loadings$pathway)
  tmp_loadings$pathway <- gsub("NervousSys", "NervSysDev", tmp_loadings$pathway)
  tmp_loadings$pathway <- gsub("SynapPlasticity", "SynapPlast", tmp_loadings$pathway)
  tmp_loadings$pathway <- gsub("SynapTransmiss", "SynapTrans", tmp_loadings$pathway)
  tmp_loadings$pathway <- gsub("VesicleTraffick", "VesiclTraf", tmp_loadings$pathway)
  tmp_loadings$pathway <- gsub("VesicleTraffick", "VesiclTraf", tmp_loadings$pathway)
  tmp_loadings$celltype[tmp_loadings$celltype == "S"] <- "     Sensori."
  tmp_loadings$celltype[tmp_loadings$celltype == "A"] <- "Assoc."
  
  cast_pathway_celltype <- reshape::cast(tmp_loadings[!tmp_loadings$celltype %in% c("     Sensori.", "Assoc."), c("pathway", "celltype", "loading")], pathway~celltype, mean)
  cast_pathway_celltype <- as.matrix(cast_pathway_celltype)
  res_celltype <- sign_bicluster_order(cast_pathway_celltype)
  
  #### brain
  tmp_loadings_brain <- tmp_loadings[tmp_loadings$celltype %in% c("     Sensori.", "Assoc."), ]
  tmp_loadings_brain$sign <- sign(tmp_loadings_brain$loading)
  tmp_loadings_brain$sign[abs(tmp_loadings_brain$loading) < 0.1] <- -1
  tmp_loadings_brain <- tmp_loadings_brain[order(tmp_loadings_brain$celltype), ]
  
  
  
  pathway_brain_sign_loading <- aggregate(sign ~ pathway, tmp_loadings_brain, paste, collapse = ",")
  cast_pathway_brain <- reshape::cast(tmp_loadings_brain[, c("pathway", "celltype", "loading")], pathway~celltype, mean)
  cast_pathway_brain <- as.matrix(cast_pathway_brain)
  res_brain <- sign_bicluster_order(cast_pathway_brain)
  
  findMax <- function(v){
    return(max(v))
  }
  pathway_brain_order <- aggregate(loading ~ pathway, tmp_loadings_brain, findMax)
  pathway_brain_order <- merge(pathway_brain_order, pathway_brain_sign_loading, by = "pathway", all.x = T)
  pathway_brain_order <- pathway_brain_order$pathway[order(pathway_brain_order$sign, pathway_brain_order$loading)]
  
  pathway_brain_order <- rev(names(res_celltype$row_order))
  
  tmp_loadings_brain$pathway <- factor(tmp_loadings_brain$pathway, levels = pathway_brain_order)
  writeLines(pathway_brain_order, sprintf("pathway_brain_order_factor_%s", i))
  
  #### celltype
  tmp_loadings_celltype <- tmp_loadings[!tmp_loadings$celltype %in% c("     Sensori.", "Assoc."), ]
  
  pathway_celltype_order <- aggregate(loading ~ pathway, tmp_loadings_celltype, mean)
  celltype_order <- aggregate(loading ~ celltype, tmp_loadings_celltype, mean)
  
  #### order by mean/ cluster
  # pathway_celltype_order <- pathway_celltype_order$pathway[order(pathway_celltype_order$loading)]
  # celltype_order <- celltype_order$celltype[order(celltype_order$loading, decreasing = T)]
  
  #### order by bicluster of positive/negative mean factor score
  pathway_celltype_order <- names(res_celltype$row_order)
  celltype_order <- names(res_celltype$col_order)
  
  tmp_loadings_celltype$pathway <- factor(tmp_loadings_celltype$pathway, levels = pathway_brain_order)
  tmp_loadings_celltype$celltype <- factor(tmp_loadings_celltype$celltype, levels = celltype_order)
  
  writeLines(pathway_celltype_order, sprintf("pathway_celltype_order_factor_%s", i))
  writeLines(celltype_order, sprintf("celltype_order_factor_%s", i))
  writeLines(as.character(c(min, max)), sprintf("minmax_factor_%s", i))
  
  
  
  # tmp_loadings$loading[tmp_loadings$loading < min] <- min
  # tmp_loadings$loading[tmp_loadings$loading > max] <- max
  

  
  p_celltype <- ggplot(tmp_loadings_celltype, aes(y = pathway, x = celltype)) + geom_tile(fill = "white", color = "grey") +
    geom_tile(aes(fill = loading), color = "black", size = .3, na.rm = F) +
    color_scale + 
    theme_bw() + theme(axis.text = element_text(size = 12),
                       axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none", 
                       plot.margin = margin(10, 5, 0, 0)) +
    xlab("") + ylab("")# + ggtitle("brain regions")
  
  p_brain <- ggplot(tmp_loadings_brain, aes(y = pathway, x = celltype)) + geom_tile(fill = "white", color = "grey") +
    geom_tile(aes(fill = loading), color = "black", size = .3, na.rm = F) +
    color_scale + 
    theme_bw() + theme(axis.text = element_text(size = 12),
                       axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.text.y = element_blank(), 
                       axis.ticks.y = element_blank(),
                       legend.position = "right", plot.margin = margin(10, 5, 0, -5)) +
    xlab("") + ylab("")# + ggtitle("cell types")
  
  p[[length(p)+1]] <- cowplot::plot_grid(p_celltype, p_brain, nrow = 1, rel_widths = c(1, 0.38))
}
  
ggpubr::ggarrange(plotlist = p, ncol = 1)

height <- 13 #main
# height <- 15 #supplement

width <- 6

# ggpubr::ggarrange(plotlist = p, nrow = 1)
ggsave(sprintf("factor_loadings_%s_mean.pdf", gsub("\\-", "_", stratification)), width = width, height = height)


#### scatter plot for factors
plot_scatter <- function(gs2, gs1, title, column, y, x, factor_loadings){
  
  load("pathway_2_3_ways_stratification.RData")
  gsLength <- sapply(gsData.s, length)
  gsLength <- data.frame("geneset" = names(gsLength), "length" = gsLength)
  
  gs1 <- merge(gs1, gsLength, by = "geneset", all.x = T)
  lm <- lm(z~length, gs1)
  gs1$z <- lm$residuals
  
  gs2 <- merge(gs2, gsLength, by = "geneset", all.x = T)
  lm <- lm(z~length, gs2)
  gs2$z <- lm$residuals
  
  gs1$fwer <- p.adjust(gs1$pvalue, method = "bonferroni")
  gs2$fwer <- p.adjust(gs2$pvalue, method = "bonferroni")
  cutoff <- 0.1
  gs <- merge(gs1[, c("geneset", "type", column, "pvalue", "qvalue", "fwer")],
              gs2[, c("geneset", "type", column, "pvalue", "qvalue", "fwer")], by = "geneset", all = T)
  gs <- na.omit(gs)
  
  names(gs) <- c("geneset", paste(c("type", "coefficient", "pvalue", "qvalue", "fwer"), x, sep="_"),
                 paste(c("type", "coefficient", "pvalue", "qvalue", "fwer"), y, sep="_"))
  
  names(factor_loadings) <- c("geneset", "loading")
  # gs <- gs[which(gs[sprintf("qvalue_%s", x)] < cutoff | gs[sprintf("qvalue_%s", y)] < cutoff), ]
  gs <- merge(gs, factor_loadings, by = "geneset", all.x = T)
  
  min <- -max(abs(gs$loading))
  max <- max(abs(gs$loading))
  # min <- -3
  # max <- 3
  gs$loading[gs$loading < min] <- min
  gs$loading[gs$loading > max] <- max
  
  color_scale <- c("#053061", "#2166AC", "#4393C3", "#72A5CE", "#D1E5F0", "#FFFFFF", "#FFF0A0", "#FDAE61", "#F47D43", "#D75027", "#A50026")
  color_scale <- scale_fill_gradientn(colors = color_scale, limits = c(min, max))
  
  if(nrow(gs) > 1){
    gs$signal <- "NS"
    gs$signal[which(gs[sprintf("qvalue_%s", x)] < cutoff)] <- "FDR10%"
    gs$signal[which(gs[sprintf("qvalue_%s", y)] < cutoff)] <- "FDR10%"
    gs$signal[which(gs[sprintf("fwer_%s", x)] < cutoff)] <- "FWER10%"
    gs$signal[which(gs[sprintf("fwer_%s", y)] < cutoff)] <- "FWER10%"
    gs$signal <- factor(gs$signal, levels = c("FWER10%", "FDR10%", "NS"))
    
    gs[sprintf("%s_sign", x)] <- sign(gs[sprintf("coefficient_%s", x)])
    gs[sprintf("%s_sign", y)] <- sign(gs[sprintf("coefficient_%s", y)])
    
    dt <- gs
    # sign.test <- signif(sum(dt[sprintf("%s_sign", x)] == dt[sprintf("%s_sign", y)], na.rm = T)/nrow(na.omit(dt)), digits = 2)
    # cor.test <- cor.test(dt[,sprintf("OR_nocovariates_%s", x)],
    #                      dt[,sprintf("OR_nocovariates_%s", y)],
    #                      method = "spearman")
    cor.test <- cor.test(dt[,sprintf("coefficient_%s", x)],
                         dt[,sprintf("coefficient_%s", y)],
                         method = "kendall")
    rho.perm <- c()
    for(perm in 1:1000){
      perm.test <- cor.test(dt[,sprintf("coefficient_%s", x)],
                            sample(dt[,sprintf("coefficient_%s", y)]),
                            method = "kendall")
      rho.perm <- c(rho.perm, perm.test$estimate)
    }
    
    rho <- signif(cor.test$estimate, digits = 2)
    pval <- signif((sum(abs(rho.perm) >= abs(rho))/length(rho.perm)), digits = 2)
    pval <- ifelse(pval > 1, 1, pval)
    pval <- ifelse(pval <= 0, "<2.2e-16", pval)
    # gs$x <- log(gs[,sprintf("OR_nocovariates_%s", x)])
    # gs$y <- log(gs[,sprintf("OR_nocovariates_%s", y)])
    
    gs$x <- gs[,sprintf("coefficient_%s", x)]
    gs$y <- gs[,sprintf("coefficient_%s", y)]
    
    maxx <- max(c(gs$x, gs$y), na.rm = T) # 
    minx <- min(c(gs$x, gs$y), na.rm = T) # 
    
    maxy <- max(c(gs$x, gs$y), na.rm = T) #
    miny <- min(c(gs$x, gs$y), na.rm = T) #
    
    if(nrow(gs)<2){
      return(NULL)
    }else{
      line_color <- ifelse(pval < 0.05, "black", NA)
      
      leg.position <- "bottom"
      gs$label <- ""
      gs$label[gs$signal == "FWER10%"] <- gs$geneset[gs$signal == "FWER10%"]
      
      p <- ggplot(gs, aes(x = x, y = y)) + 
        geom_hline(yintercept = 0, lty = 1, lwd = .2) + geom_vline(xintercept = 0, lty = 1, lwd = .2) + 
        geom_smooth(data = gs,
                    method = "lm", size = .7, color = line_color, fill="lightgrey", alpha=.5, fullrange=TRUE) +
        geom_point(aes(fill = loading, size = signal, alpha= signal, color = signal), shape=21) + theme_bw() +
        ggtitle(sprintf("%s\ntau=%s, p=%s", title, rho, pval)) + 
        xlab(x) + ylab(y) + theme(legend.position = leg.position, title = element_text(size = 10), 
                                  legend.title = element_blank(), axis.text = element_text(size = 10)) +
        scale_y_continuous(labels = label_number(accuracy = 0.01), limits = c(miny-1, maxy+1)) +
        scale_x_continuous(labels = label_number(accuracy = 0.01), limits = c(miny-1, maxy+1)) +
        coord_cartesian(xlim=c(minx-0.01, maxx+0.01), ylim=c(miny-0.01, maxy+0.01)) +
        color_scale + 
        # geom_text_repel(aes(label = label), alpha=.7,  max.overlaps = 100, size = 3.5, box.padding = 0.5,
        #                 direction    = "y",
        #                 hjust        = -1, nudge_y = asd.del$nudge_y,
        #                 segment.size = 0.4, force_pull = 0,
        #                 point.padding = 0.2, 
        #                 nudge_x = asd.del$nudge_x,
        #                 segment.curvature = -1e-20,
        #                 max.iter = 1e4, max.time = 1) +
        scale_size_manual(values=c("FWER10%"=2, "FDR10%"=2, "NS"=2), drop = F)+ 
        scale_alpha_manual(values=c("FWER10%"=1, "FDR10%"=1, "NS"=.9), drop = F)+
        scale_color_manual(values=c("FWER10%"="black", "FDR10%"="darkgrey", "NS"="white"), drop = F)+
        guides(
               size = guide_legend(ncol = 3, drop = F), 
               alpha = guide_legend(ncol = 3, override.aes = c("fill" = "darkgrey"), drop = F))
      
      return(p)
    }
  }else{
    return(NULL)
  }
}

##factor 1
# gs2 <- dt.plot[dt.plot$disorder == "SCZ" & dt.plot$type == "DEL" & !is.na(dt.plot$celltype), ]
# gs1 <- dt.plot[dt.plot$disorder == "SCZ" & dt.plot$type == "DUP", ]
# y <- "SCZ_DEL"
# x <- "SCZ_DUP"
# factor_loadings <- loadings[, c("V1", "ML1")]
# column <- "z"
# title <- "SCZDELvSCZDUP"

p1 <- plot_scatter(dt.plot[dt.plot$disorder == "SCZ" & dt.plot$type == "DUP", ],
             dt.plot[dt.plot$disorder == "SCZ" & dt.plot$type == "DEL", ],
             "SCZ DUP vs SCZ DEL", "z",
             "SCZ DUP effect size (Z-score)", "SCZ DEL effect size (Z-score)", loadings[, c("V1", "ML1")])

p2 <- plot_scatter(dt.plot[dt.plot$disorder == "MDD" & dt.plot$type == "DUP", ],
             dt.plot[dt.plot$disorder == "BD" & dt.plot$type == "DUP", ],
             "MDD DUP vs BD DUP", "z",
             "MDD DUP effect size (Z-score)", "BD DUP effect size (Z-score)", loadings[, c("V1", "ML2")])

p3 <- plot_scatter(dt.plot[dt.plot$disorder == "MDD" & dt.plot$type == "DEL", ],
             dt.plot[dt.plot$disorder == "ADHD" & dt.plot$type == "DEL", ],
             "MDD DEL vs ADHD DEL", "z",
             "MDD DEL effect size (Z-score)", "ADHD DEL effect size (Z-score)", loadings[, c("V1", "ML3")])

cowplot::plot_grid(p1, p2, p3, nrow = 1, ncol = 3)
ggsave("scatter_plots_2way.pdf", width = 9, height = 4)
