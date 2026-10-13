###
### 031 Cell_type_proportion_RNA.R
###
# Purpose: Cell type proportion analysis in single nuclei RNA dataset
# Dependencies:
curr_dir <- "/Tsai" #output results to specific folder similarly "/Jackson", "Jco_cult"
set.seed(1234)
##Cell proportion
####Per subtype

library(limma)
library(ggplot2)
library(ComplexHeatmap)
library(viridis)
#cell counts
md <- data.frame(Org_ChP@meta.data$Subcelltype, Org_ChP@meta.data$sample, Org_ChP@meta.data$Genotype)
colnames(md) <- c("Subcelltype", "sample", "Genotype")
prop <- tabyl(md, Subcelltype, sample)
prop <- prop %>% column_to_rownames("Subcelltype")

#metadata
mdat <- data.frame(samples = colnames(prop), row.names = colnames(prop))
mdat <- mdat %>% mutate(Genotype = ifelse(grepl("E3", samples), "33", "44")) #Extract APOE33 and APOE44 information
mdat <- mdat %>% mutate(Differentiation = t(data.frame(regmatches(samples, gregexpr("[[:digit:]]+", samples)))[2,])) #Extract differentiation number
mdat <- mdat %>% mutate(Genotype_Differentiation = paste0(Genotype, Differentiation)) #merge Genotype and differentiation number to specify single differentiation.

#limma
D = factor(paste(mdat$Differentiation))
G = factor(paste(mdat$Genotype))
GD = factor(paste(mdat$Genotype_Differentiation))

Maindm <- model.matrix(~0 + D + G) #!!!!
rownames(Maindm) <- mdat$samples
colnames(Maindm)
Maincm <- makeContrasts(levels=Maindm, MainGenotype = (G44))

vdat <- voom(prop, design = Maindm)
corfit <- duplicateCorrelation(vdat, block=GD)

Mainf <- lmFit(vdat, Maindm, weights = vdat$weights, block=GD, correlation=corfit$consensus.correlation)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

df <- topTable(Mainf, coef = "MainGenotype", p.value=1, number=Inf) #BAM_2 
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions_MainGenotype_EFI.csv", sep=""))

Idents(Org_ChP) <- Org_ChP$Majorcelltype
Ncellfun(Org_ChP, "CellType/Publish/Major_Genotype.csv", group.by="Genotype")

Idents(Org_ChP) <- Org_ChP$Majorcelltype
Ncellfun(Org_ChP, "CellType/Publish/Major_Sample.csv", group.by="sample")

Idents(Org_ChP) <- Org_ChP$Subcelltype
Ncellfun(Org_ChP, "CellType/Publish/Sub_Genotype.csv", group.by="Genotype")


###Imm

Imm_prop <- data.frame(prop[c(8:10),])

#metadata
mdat <- data.frame(samples = colnames(Imm_prop), row.names = colnames(Imm_prop))
mdat <- mdat %>% mutate(Genotype = ifelse(grepl("E3", samples), "33", "44")) #Extract APOE33 and APOE44 information
mdat <- mdat %>% mutate(Differentiation = t(data.frame(regmatches(samples, gregexpr("[[:digit:]]+", samples)))[2,])) #Extract differentiation number
mdat <- mdat %>% mutate(Genotype_Differentiation = paste0(Genotype, Differentiation)) #merge Genotype and differentiation number to specify single differentiation.

#limma
D = factor(paste(mdat$Differentiation))
G = factor(paste(mdat$Genotype))
GD = factor(paste(mdat$Genotype_Differentiation))

Maindm <- model.matrix(~0 + D + G) #!!!!
rownames(Maindm) <- mdat$samples
colnames(Maindm)
Maincm <- makeContrasts(levels=Maindm, MainGenotype = (G44))

vdat <- voom(Imm_prop, design = Maindm)
heatmap(cor(vdata), Colv=NA, Rowv=NA, symm=T, main="Sample correlations")

corfit <- duplicateCorrelation(vdat, block=GD)

Mainf <- lmFit(vdat, Maindm, weights = vdat$weights, block=GD, correlation=corfit$consensus.correlation)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

df <- topTable(Mainf, coef = "MainGenotype", p.value=1, number=Inf) #BAM_2 
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions_Imm_MainGenotype.csv", sep=""))


###Epi

Epi_prop <- data.frame(prop[c(1:4),])

#metadata
mdat <- data.frame(samples = colnames(Epi_prop), row.names = colnames(Epi_prop))
mdat <- mdat %>% mutate(Genotype = ifelse(grepl("E3", samples), "33", "44")) #Extract APOE33 and APOE44 information
mdat <- mdat %>% mutate(Differentiation = t(data.frame(regmatches(samples, gregexpr("[[:digit:]]+", samples)))[2,])) #Extract differentiation number
mdat <- mdat %>% mutate(Genotype_Differentiation = paste0(Genotype, Differentiation)) #merge Genotype and differentiation number to specify single differentiation.

#limma
D = factor(paste(mdat$Differentiation))
G = factor(paste(mdat$Genotype))
GD = factor(paste(mdat$Genotype_Differentiation))

Maindm <- model.matrix(~0 + D + G) #!!!!
rownames(Maindm) <- mdat$samples
colnames(Maindm)
Maincm <- makeContrasts(levels=Maindm, MainGenotype = (G44))

vdat <- voom(Epi_prop, design = Maindm)
heatmap(cor(vdata), Colv=NA, Rowv=NA, symm=T, main="Sample correlations")

corfit <- duplicateCorrelation(vdat, block=GD)

Mainf <- lmFit(vdat, Maindm, weights = vdat$weights, block=GD, correlation=corfit$consensus.correlation)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

df <- topTable(Mainf, coef = "MainGenotype", p.value=1, number=Inf) #BAM_2 
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions_Epi_MainGenotype.csv", sep=""))


###Fib

Fib_prop <- data.frame(prop[c(5:7),])

#metadata
mdat <- data.frame(samples = colnames(Fib_prop), row.names = colnames(Fib_prop))
mdat <- mdat %>% mutate(Genotype = ifelse(grepl("E3", samples), "33", "44")) #Extract APOE33 and APOE44 information
mdat <- mdat %>% mutate(Differentiation = t(data.frame(regmatches(samples, gregexpr("[[:digit:]]+", samples)))[2,])) #Extract differentiation number
mdat <- mdat %>% mutate(Genotype_Differentiation = paste0(Genotype, Differentiation)) #merge Genotype and differentiation number to specify single differentiation.

#limma
D = factor(paste(mdat$Differentiation))
G = factor(paste(mdat$Genotype))
GD = factor(paste(mdat$Genotype_Differentiation))

Maindm <- model.matrix(~0 + D + G) #!!!!
rownames(Maindm) <- mdat$samples
colnames(Maindm)
Maincm <- makeContrasts(levels=Maindm, MainGenotype = (G44))

vdat <- voom(Fib_prop, design = Maindm)
heatmap(cor(vdata), Colv=NA, Rowv=NA, symm=T, main="Sample correlations")

corfit <- duplicateCorrelation(vdat, block=GD)

Mainf <- lmFit(vdat, Maindm, weights = vdat$weights, block=GD, correlation=corfit$consensus.correlation)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

df <- topTable(Mainf, coef = "MainGenotype", p.value=1, number=Inf) #BAM_2 
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions_Fib_MainGenotype.csv", sep=""))


###Prog

Prog_prop <- data.frame(prop[c(11:15),])

#metadata
mdat <- data.frame(samples = colnames(Prog_prop), row.names = colnames(Prog_prop))
mdat <- mdat %>% mutate(Genotype = ifelse(grepl("E3", samples), "33", "44")) #Extract APOE33 and APOE44 information
mdat <- mdat %>% mutate(Differentiation = t(data.frame(regmatches(samples, gregexpr("[[:digit:]]+", samples)))[2,])) #Extract differentiation number
mdat <- mdat %>% mutate(Genotype_Differentiation = paste0(Genotype, Differentiation)) #merge Genotype and differentiation number to specify single differentiation.

#limma
D = factor(paste(mdat$Differentiation))
G = factor(paste(mdat$Genotype))
GD = factor(paste(mdat$Genotype_Differentiation))

Maindm <- model.matrix(~0 + D + G) #!!!!
rownames(Maindm) <- mdat$samples
colnames(Maindm)
Maincm <- makeContrasts(levels=Maindm, MainGenotype = (G44))

vdat <- voom(Prog_prop, design = Maindm)
heatmap(cor(vdata), Colv=NA, Rowv=NA, symm=T, main="Sample correlations")

corfit <- duplicateCorrelation(vdat, block=GD)

Mainf <- lmFit(vdat, Maindm, weights = vdat$weights, block=GD, correlation=corfit$consensus.correlation)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

df <- topTable(Mainf, coef = "MainGenotype", p.value=1, number=Inf) #BAM_2 
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions_Prog_MainGenotype.csv", sep=""))


###Dividing

Divi_prop <- data.frame(prop[c(16:19),])

#metadata
mdat <- data.frame(samples = colnames(Divi_prop), row.names = colnames(Divi_prop))
mdat <- mdat %>% mutate(Genotype = ifelse(grepl("E3", samples), "33", "44")) #Extract APOE33 and APOE44 information
mdat <- mdat %>% mutate(Differentiation = t(data.frame(regmatches(samples, gregexpr("[[:digit:]]+", samples)))[2,])) #Extract differentiation number
mdat <- mdat %>% mutate(Genotype_Differentiation = paste0(Genotype, Differentiation)) #merge Genotype and differentiation number to specify single differentiation.

#limma
D = factor(paste(mdat$Differentiation))
G = factor(paste(mdat$Genotype))
GD = factor(paste(mdat$Genotype_Differentiation))

Maindm <- model.matrix(~0 + D + G) #!!!!
rownames(Maindm) <- mdat$samples
colnames(Maindm)
Maincm <- makeContrasts(levels=Maindm, MainGenotype = (G44))

vdat <- voom(Divi_prop, design = Maindm)
heatmap(cor(vdata), Colv=NA, Rowv=NA, symm=T, main="Sample correlations")

corfit <- duplicateCorrelation(vdat, block=GD)

Mainf <- lmFit(vdat, Maindm, weights = vdat$weights, block=GD, correlation=corfit$consensus.correlation)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

df <- topTable(Mainf, coef = "MainGenotype", p.value=1, number=Inf) #BAM_2 
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions_Divi_MainGenotype.csv", sep=""))


####Per major celltype

md <- data.frame(Org_ChP@meta.data$Majorcelltype, Org_ChP@meta.data$sample, Org_ChP@meta.data$Genotype)
colnames(md) <- c("Majorcelltype", "sample", "Genotype")
prop <- tabyl(md, Majorcelltype, sample)
prop <- prop %>% column_to_rownames("Majorcelltype")

#metadata
mdat <- data.frame(samples = colnames(prop), row.names = colnames(prop))
mdat <- mdat %>% mutate(Genotype = ifelse(grepl("E3", samples), "33", "44")) #Extract APOE33 and APOE44 information
mdat <- mdat %>% mutate(Differentiation = t(data.frame(regmatches(samples, gregexpr("[[:digit:]]+", samples)))[2,])) #Extract differentiation number
mdat <- mdat %>% mutate(Genotype_Differentiation = paste0(Genotype, Differentiation)) #merge Genotype and differentiation number to specify single differentiation.

#limma
D = factor(paste(mdat$Differentiation))
G = factor(paste(mdat$Genotype))
GD = factor(paste(mdat$Genotype_Differentiation))

Maindm <- model.matrix(~0 + D + G) #!!!!
rownames(Maindm) <- mdat$samples
colnames(Maindm)
Maincm <- makeContrasts(levels=Maindm, MainGenotype = (G44))

vdat <- voom(prop, design = Maindm)
heatmap(cor(vdata), Colv=NA, Rowv=NA, symm=T, main="Sample correlations")

corfit <- duplicateCorrelation(vdat, block=GD)

Mainf <- lmFit(vdat, Maindm, weights = vdat$weights, block=GD, correlation=corfit$consensus.correlation)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

df <- topTable(Mainf, coef = "MainGenotype", p.value=1, number=Inf) #BAM_2 
df
write.csv(df, paste(curr_dir, "/CellType/Publish/Proportions_Major_MainGenotype.csv", sep=""))

