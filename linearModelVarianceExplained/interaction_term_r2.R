# Load packages
library(lme4)
library(lmerTest)
library(car)
library(ggplot2)
library(emmeans)
library(Biobase)
library(data.table)
library(MuMIn)

# Read in your data
data.all <- fread("association_results_2_3_ways_withcontributegenes_atp0.05.tsv", data.table = F)
data.all <- data.all[data.all$set == "all", ]
# data.all <- data.all[data.all$stratification != "celltype-sa-axis", c("geneset", "type", "z", "pvalue", "qvalue", "stratification", "disorder",
#                                                                       "pathway", "celltype", "brain")]
data.all <- data.all[, c("geneset", "type", "z", "pvalue", "qvalue", "stratification", "disorder",
                                                                      "pathway", "celltype", "brain")]
# data.all <- data.all[data.all$stratification == "three-way", ]

data.all$brain <- factor(data.all$brain)
data.all$dosage <- factor(data.all$type)
data.all$celltype <- factor(data.all$celltype)

tests <-  list(#"pathway", "celltype", "dosage",
               c("pathway", "dosage"), 
               c("pathway", "brain"), 
               c("pathway", "celltype"), 
               c("brain", "dosage"),
               c("celltype", "dosage"),
               c("celltype", "brain"),
               c("pathway", "brain", "dosage"), 
               c("pathway", "celltype", "dosage"),
               c("celltype", "brain", "dosage"),
               c("pathway", "celltype", "brain"))

### by likelihood ratio test
dt.lh <- data.frame()
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
    
    p_factor <- 1
    if(length(test) == 1){
      
    }else{
      for(i in 2:length(test)){
        p_factor <- p_factor * length(unique(data[, test[i]]))
      }
    }
      
    if(length(test) == 3){
      test_additive <- paste(test[1], "*", test[2], "+", test[1], "*", test[3], "+", test[3], "*", test[2], "+", paste(test, collapse = "+"))
    }else if(length(test) == 2){
      test_additive <- paste(test, collapse = "+")
    }else{
      test_additive <- "1"
    }
    
    test <- paste(test, collapse = "*")
    
    
    formula_additive <- sprintf("z ~ %s", test_additive)
    formula <- sprintf("z ~ %s", test)
    
    set.seed(123)
    r2.bootstrap <- c()
    for(i in 1:1000){
      data.bootstrap <- data[sample(1:nrow(data), replace = T),]
      model_additive <- lm(formula = formula_additive, data.bootstrap)
      model <- lm(formula = formula, data.bootstrap)
      
      r2.bootstrap <- c(r2.bootstrap, summary(model)$r.squared - summary(model_additive)$r.squared)
    }
    
    
    # if(!is.na(addition)){
    #   formula <- paste(formula, " + ", addition)
    # }
    model <- lm(formula = formula, data)
    model_additive <- lm(formula = formula_additive, data)
    
    
    
    # adj.model.r2 <- 1 - ((1-summary(model)$r.squared)*nrow(data))/(nrow(data) - p_factor - 1)
    # adj.model_additive.r2 <- 1 - ((1-summary(model_additive)$r.squared)*nrow(data))/(nrow(data) - p_factor - 1)
    
    r2 <- summary(model)$r.squared - summary(model_additive)$r.squared
    # adj.r2 <- adj.model.r2 - adj.model_additive.r2
    
    se <- sd (r2.bootstrap)
    ci <- 1.96*se
    
    p <- anova(model_additive, model, test = "Chisq")$`Pr(>Chi)`[2]
    
    dt.lh <- rbind(dt.lh, data.frame(disorder, test, r2, "r2.lower" = r2-ci, "r2.upper" = r2+ci, p))
  }
}

write.table(dt.lh, "variance.tsv", sep="\t", row.names = F ,quote=F, col.names=T)
# ggplot(dt.lh, aes(x = disorder, y = adj.r2)) + geom_bar(stat = "identity", aes(fill = test), position = "dodge") +
#   theme_bw()
dt.lh$test <- factor(dt.lh$test, levels = unique(dt.lh$test))
dt.lh$FWER <- p.adjust(dt.lh$p, method = "bonferroni")
dt.lh$label <- ""
dt.lh$label[dt.lh$p < 0.05] <- "*"
dt.lh$label[dt.lh$FWER < 0.05] <- "**"
dt.lh$disorder <- factor(dt.lh$disorder, levels = c("ASD", "SCZ", "MDD", "PTSD", "ADHD", "BD"))

col.scheme <- c("pathway*dosage" = "#dc96c8",
                "pathway*brain" = "#508ab0",
                "pathway*celltype" = "#a487bd",
                "celltype*brain" = "#57a456",
                "brain*dosage" = "red",
                "celltype*dosage" = "grey",
                "pathway*brain*dosage" = "#b5b656",
                "pathway*celltype*dosage" = "#b1dbe1",
                "celltype*brain*dosage" = "#d28a5c",
                "pathway*celltype*brain" = "#f4C261")

# dt.individual <- dt.lh[dt.lh$test %in% c("pathway", "celltype", "dosage"), ]

dt.lh$test <- factor(dt.lh$test, levels = c("pathway*dosage", "pathway*celltype", "pathway*brain",
                                            "celltype*dosage", "celltype*brain", "brain*dosage",
                                            "pathway*celltype*dosage", "pathway*brain*dosage", "celltype*brain*dosage", "pathway*celltype*brain"))
dt.lh$r2.lower[dt.lh$r2.lower < 0] <- 0
write.table(dt.lh, "interaction_r2.tsv", sep="\t", row.names=F, quote=F, col.names=T)

ggplot(dt.lh, aes(x = disorder, y = r2)) + geom_bar(stat = "identity", aes(fill = test), position = "dodge", color = "white", width = .8) +
  geom_errorbar(aes(ymin = r2.lower, ymax = r2.upper, fill = test), position = position_dodge(width = .8), width =.2, lwd = .2) +
  theme_bw() + #geom_text(aes(label = label, fill = test, y = r2.upper), position = position_dodge(width = .8)) +
  scale_fill_manual(values = col.scheme, name = "", 
                    guide=guide_legend(nrow = 3, order = T, byrow = T)) + theme(legend.position = "bottom")
ggsave("pathway_stratification_interaction.pdf", width = 8, height = 3)

