# Load packages
library(lme4)
library(lmerTest)
library(car)
library(ggplot2)
library(emmeans)
library(Biobase)
library(data.table)
# Read in your data
data.all <- fread("association_results_2_3_ways_withcontributegenes_atp0.05.tsv", data.table = F)
data.all <- data.all[data.all$set == "all", ]

data.all <- data.all[, c("geneset", "type", "z", "pvalue", "qvalue", "stratification", "disorder",
                         "pathway", "celltype", "brain")]

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
  c("pathway", "celltype", "brain"),
  c("pathway", "celltype", "brain", "dosage"))

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
    
    test <- paste(test, collapse = "*")
    
    formula <- sprintf("z ~ %s", test)
    
    model <- lm(formula = formula, data)

    r2 <- summary(model)$r.squared
    
    if(test == "pathway*celltype*brain*dosage"){
      ci <- c(1, 1)
    }else{
      # Perform bootstrap
      set.seed(123)
      r2.bootstrap <- c()
      for(i in 1:1000){
        data.bootstrap <- data[sample(1:nrow(data), replace = T),]
        model <- lm(formula = formula, data.bootstrap)
        r2.bootstrap <- c(r2.bootstrap, summary(model)$r.squared)
      }
      
      # Get 95% CI
      se <- sd (r2.bootstrap)
      ci <- 1.96*se
      ci <- c(r2-ci, r2+ci)
    }
    
    
    dt.lh <- rbind(dt.lh, data.frame(disorder, test, r2, "r2.lower" = ci[1], "r2.upper" = ci[2]))
  }
}

write.table(dt.lh, "r2_all_models.tsv", sep="\t", row.names = F ,quote=F, col.names=T)

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
                "pathway*celltype*brain" = "#f4C261",
                "pathway*celltype*brain*dosage" = "#2ec4b6")

# dt.individual <- dt.lh[dt.lh$test %in% c("pathway", "celltype", "dosage"), ]

dt.lh$test <- factor(dt.lh$test, levels = c("pathway*dosage", "pathway*celltype", "pathway*brain",
                                            "celltype*dosage", "celltype*brain", "brain*dosage",
                                            "pathway*celltype*dosage", "pathway*brain*dosage", "celltype*brain*dosage", 
                                            "pathway*celltype*brain", "pathway*celltype*brain*dosage"))
dt.lh$r2.lower[dt.lh$r2.lower < 0] <- 0
dt.lh$r2.upper[dt.lh$r2.lower > 1] <- 1

write.table(dt.lh, "full_model_r2.tsv", sep="\t", row.names=F, quote=F, col.names=T)


ggplot(dt.lh, aes(x = disorder, y = r2)) + geom_bar(stat = "identity", aes(fill = test), position = "dodge", color = "white", width = .8) +
  theme_bw() + 
  geom_errorbar(aes(ymin = r2.lower, ymax = r2.upper, fill = test), position = position_dodge(width = .8), width =.2, lwd = .2) +
  scale_fill_manual(values = col.scheme, name = "", 
                    guide=guide_legend(nrow = 3, order = T, byrow = F)) + theme(legend.position = "bottom")

ggsave("pathway_stratification_fullmodel.pdf", width = 8, height = 3)

