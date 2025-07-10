library(ggplot2)
library(data.table)

# Read in your data
data.all <- fread("association_results_2_3_ways_withcontributegenes_atp0.05.tsv", data.table = F)
data.all <- data.all[data.all$set == "all", ]
data.all <- data.all[, c("geneset", "type", "z", "pvalue", "qvalue", "stratification", "disorder",
                                                                      "pathway", "celltype", "brain")]

data.all$brain <- factor(data.all$brain)
data.all$dosage <- factor(data.all$type)
data.all$celltype <- factor(data.all$celltype)

tests <-  list(c("pathway", "brain", "dosage"), 
               c("pathway", "celltype", "dosage"))

### by likelihood ratio test
res <- data.frame()
for(disorder in c("ASD", "SCZ", "BD", "ADHD", "MDD", "PTSD")){
  # for(stratification in c("pathway-sa-axis", "pathway-celltype")){
  for(test in tests){
    exclude <- c()
    for(subtest in test){
      exclude <- union(exclude, which(is.na(data.all[, c(subtest)])))
    }
    
    if(length(exclude) > 0){
      data <- data.all[-exclude, ]
    }else{
      data <- data.all
    }
    
    
    data <- data[which(data$disorder == disorder), ]
    
    m1.formula <- paste(test, collapse = "*")
    m2.formula <- paste(test[-2], collapse = "*")
    m3.formula <- paste(test[2], paste(test[-2], collapse = "*"), sep="+")
    
    m1.model <- lm(sprintf("z ~ %s", m1.formula), data)
    m2.model <- lm(sprintf("z ~ %s", m2.formula), data)
    m3.model <- lm(sprintf("z ~ %s", m3.formula), data)
    
    ano.res <- anova(m2.model, m3.model, m1.model)
    
    res <- rbind(res, data.frame(disorder, "model" = sprintf("%s v %s", m2.formula, m3.formula),
                                 "test" = paste(test[2], "main effect"),
                                 "p" = signif(ano.res$`Pr(>F)`[2], digits=2)))
    res <- rbind(res, data.frame(disorder, "model" = sprintf("%s v %s", m3.formula, m1.formula), 
                                 "test" = paste(test[2], "interaction effect"),
                                 "p" = signif(ano.res$`Pr(>F)`[3], digits = 2)))
    
  }
}

res$sig <- F
res$sig[res$p < 0.05] <- T

write.table(res, "model_comparison.tsv", sep="\t", row.names = F ,quote=F, col.names=T)

