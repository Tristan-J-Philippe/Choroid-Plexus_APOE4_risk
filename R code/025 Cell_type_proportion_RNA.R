###
### 031 Cell_type_proportion_RNA.R
###
# Purpose: Cell type proportion analysis in single nuclei RNA dataset
# Dependencies:
curr_dir <- "/ExVivo" #output results to specific folder
set.seed(1234)
##Cell proportion
####Per subtype

library(limma)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(viridis)

#load up to date meta
mdat <- read_csv(paste0(curr_dir, "/metadata/ex vivo ChP metadata (analytical) newAPOE RNAonly.csv", col_type=cols(.default="c")))

#cell counts
md <- data.frame(Ex_ChP@meta.data$Subcelltype, Ex_ChP@meta.data$Lane_projid, Ex_ChP@meta.data$projid)
colnames(md) <- c("active.ident", "Lane_projid", "projid")
prop <- tabyl(md, active.ident, Lane_projid)
prop <- prop %>% column_to_rownames("active.ident")
parenchyma <- prop[c("Neuro", "Astro"),]
parenchyma <- colSums(parenchyma)
prop <- prop[!rownames(prop)%in%c("Neuro", "Astro", "Oligo"),]

#discrete metadata
Ex_ChP_md <- left_join(md, mdat, by="projid")
Ex_ChP_md$active.ident <- NULL
Ex_ChP_md <- distinct(Ex_ChP_md, .keep_all = T)

#Factors
msex = factor(paste(Ex_ChP_md$msex))
X10X_batch = factor(paste(Ex_ChP_md$X10X_batch))
library_batch = factor(paste(Ex_ChP_md$library_batch))
apoe = factor(paste(Ex_ChP_md$apoe_genotype))
cogdx = factor(paste(Ex_ChP_md$cogdx))
Cog_Path = factor(paste(Ex_ChP_md$Cog_Path))
ad_adnc <- factor(paste(Ex_ChP_md$ad_adnc))

#continuous
age = as.numeric(Ex_ChP_md$age_death)
pmi = as.numeric(Ex_ChP_md$pmi)
educ = as.numeric(Ex_ChP_md$educ)
cogn_global = as.numeric(Ex_ChP_md$cogn_global_last)
parenchyma = as.numeric(parenchyma)
tdp_st4 = as.numeric(Ex_ChP_md$tdp_st4)

Lane <- factor(paste(sapply(str_split(colnames(prop), "_"), '[', 2), sapply(str_split(colnames(prop), "_"), '[', 3), sep = "_"))

####APOE

library(lme4)
library(lmerTest)
library(emmeans)

Ex_ChP_md <- Ex_ChP_md[order(match(Ex_ChP_md$Lane_projid, colnames(prop))),]
Ex_ChP_md$dx_group <- factor(ifelse(Ex_ChP_md$cogdx == 1, "NCI",
                                    ifelse(Ex_ChP_md$cogdx %in% c(2, 3), "MCI",
                                           ifelse(Ex_ChP_md$cogdx %in% c(4, 5), "AD", NA))),
                             levels = c("NCI", "MCI", "AD"))
Ex_ChP_md$apoe_genotype[is.na(Ex_ChP_md$apoe_genotype)] <- "23"
Ex_ChP_md$apoe_genotype <- relevel(factor(Ex_ChP_md$apoe_genotype), ref = "33")

cols <- c("contrast", "estimate", "SE", "df", "CI.L", "CI.R", "t.ratio", "p.value")
c_apoe <- data.frame(matrix(ncol = length(cols), nrow = 0))
colnames(c_apoe) <- cols
c_cogdx <- c_apoe

cols <- c("apoe_pairwise", "dx_group_identity", "estimate", "SE", "df", "CI.L", "CI.R", "t.ratio", "p.value")
c_apoeBYcogdx <- data.frame(matrix(ncol = length(cols), nrow = 0))
colnames(c_apoeBYcogdx) <- cols

prop_percent <- as.data.frame(proportions(as.matrix(prop), margin=2))*100

for(celltype in rownames(prop_percent)){
  dat <- data.frame(proportion = as.numeric(prop_percent[celltype, ]),
                    projid = Ex_ChP_md$projid,
                    Lane = Lane,
                    apoe = Ex_ChP_md$apoe_genotype,
                    dx_group = Ex_ChP_md$dx_group,
                    msex = Ex_ChP_md$msex,
                    age = as.numeric(Ex_ChP_md$age_death),
                    pmi = as.numeric(Ex_ChP_md$pmi),
                    parenchyma = parenchyma)
  
  fit <- lmer(proportion ~ apoe * dx_group + msex + age + pmi + parenchyma + (1 | projid) + (1 | Lane), data=dat)
  
  emm <- emmeans(fit, ~ apoe)
  int <- contrast(emm, method = list(apoe = c(-1, 0, 1)))
  c_apoe[celltype,] <- data.frame(summary(int, infer = c(TRUE, TRUE)))
  
  
  emm <- emmeans(fit, ~ dx_group)
  int <- contrast(emm, method = list(AD_vs_NCI = c(-1, 0, 1)))
  c_cogdx[celltype,] <- data.frame(summary(int, infer = c(TRUE, TRUE)))
  
  
  emm <- emmeans(fit, ~ apoe * dx_group, at = list(apoe = c("33", "34"), dx_group=c("AD", "NCI")))
  int <- contrast(emm, interaction = c("pairwise", "pairwise"))
  c_apoeBYcogdx[celltype,] <- data.frame(summary(int, infer = c(TRUE, TRUE)))
}

c_apoe$FDR <- ave(c_apoe$p.value, c_apoe$contrast, FUN = p.adjust, method = "BH")
c_cogdx$FDR <- ave(c_cogdx$p.value, c_cogdx$contrast, FUN = p.adjust, method = "BH")
c_apoeBYcogdx$FDR <- ave(c_apoeBYcogdx$p.value, c_apoeBYcogdx$apoe_pairwise, FUN = p.adjust, method = "BH")
c_apoe
c_cogdx
c_apoeBYcogdx

library(ggplot2)
library(patchwork)
plts1 <- ggplot(c_apoe, aes(x = estimate, y = rownames(c_apoe))) + geom_pointrange(aes(xmin = CI.L, xmax = CI.R)) + theme_minimal() + 
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") + labs(x="apoe34v33", y="Sub Cell Type", title="Percent proportion")
plts2 <- ggplot(c_cogdx, aes(x = estimate, y = rownames(c_cogdx))) + geom_pointrange(aes(xmin = CI.L, xmax = CI.R)) + theme_minimal() + geom_vline(xintercept = 0, linetype = "dotted", color = "black") + labs(x="ADvNCI", y="Sub Cell Type", title="Percent proportion")
plts3 <- ggplot(c_apoeBYcogdx, aes(x = -estimate, y = rownames(c_apoeBYcogdx))) + geom_pointrange(aes(xmin = -CI.L, xmax = -CI.R)) + theme_minimal() + geom_vline(xintercept = 0, linetype = "dotted", color = "black") + labs(x="apoe by NCI", y="Sub Cell Type", title="Percent proportion")

plts1+plts2+plts3

ggsave2(paste0(curr_dir, "/CellType/Publish/Proportions/apoe_cogdx_cogdxBYapoe.pdf"), plts1+plts2+plts3)
write.csv(c_apoe, paste0(curr_dir, "/CellType/Publish/Proportions/apoe_cogdx_cogdxBYapoe_apoe.cvs"))
write.csv(c_cogdx, paste0(curr_dir, "/CellType/Publish/Proportions/apoe_cogdx_cogdxBYapoe_cogdx.csv"))
write.csv(c_apoeBYcogdx, paste0(curr_dir, "/CellType/Publish/Proportions/apoe_cogdx_cogdxBYapoe_apoeBYcogdx.csv"))
