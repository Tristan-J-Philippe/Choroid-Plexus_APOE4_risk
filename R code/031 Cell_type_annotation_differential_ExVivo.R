###
### 032 Cell_type_annotation_differential_RNA.R 
###
# Purpose: Cell type annotation differential analysis and gsea of single nuclei RNA dataset.
# Dependencies:
source("/030 limma_pattern_fgsea_util.R", verbose = FALSE)
curr_dir <- "/ExVivo" #output results to specific folder
set.seed(42)


###Subcelltype
cts_fun <- function(celltype){
  sub <- subset(Ex_ChP, Subcelltype==celltype)
  cts <- AggregateExpression(sub, group.by = "library_batch", assay = "soupXcounts", slot = "counts", return.seurat = FALSE)
  cts <- cts[["soupXcounts"]]
  assign(c(paste("cts_soupX", celltype, sep = "_")), cts, envir = .GlobalEnv)
  pct <- pctincell(sub)
  assign(c(paste("pct", celltype, sep = "_")), pct, envir = .GlobalEnv)
}

cts_fun("Epi_1")
gc()
cts_fun("Epi_2")
gc()
cts_fun("Epi_3")
gc()
cts_fun("Epi_4")
gc()
cts_fun("Fib_1")
cts_fun("Fib_2")
cts_fun("Fib_3")
cts_fun("Endo_1")
cts_fun("Endo_3")
cts_fun("Endo_2")
cts_fun("Mural_1")
cts_fun("Mural_2")
cts_fun("Mural_3")
cts_fun("BAM_1")
cts_fun("BAM_2")
cts_fun("BAM_3")
cts_fun("T_Cell")
cts_fun("Oligo")
cts_fun("Astro")
cts_fun("Neuro")

dim(cts_soupX_Epi_1)
dim(cts_soupX_Epi_2)
dim(cts_soupX_Epi_3)
dim(cts_soupX_Epi_4)
dim(cts_soupX_Fib_1)
dim(cts_soupX_Fib_2)
dim(cts_soupX_Fib_3)
dim(cts_soupX_Endo_1)
dim(cts_soupX_Endo_3)
dim(cts_soupX_Endo_2)
dim(cts_soupX_Mural)
dim(cts_soupX_BAM_1)
dim(cts_soupX_BAM_2)
dim(cts_soupX_BAM_3)
dim(cts_soupX_T_Cell)
dim(cts_soupX_Oligo)
dim(cts_soupX_Astro)
dim(cts_soupX_Neuro)

cts_sub_paren <- cbind(cts_soupX_Epi_1, cts_soupX_Epi_2, cts_soupX_Epi_3, cts_soupX_Epi_4, cts_soupX_Fib_1, cts_soupX_Fib_2, cts_soupX_Fib_3, cts_soupX_Fib_4, cts_soupX_Endo_1, cts_soupX_Endo_3, cts_soupX_Endo_2, cts_soupX_Mural_1, cts_soupX_Pericyte, cts_soupX_BAM_1, cts_soupX_BAM_2, cts_soupX_BAM_3, cts_soupX_T_Cell, cts_soupX_Oligo, cts_soupX_Astro, cts_soupX_Neuro)
colnames(cts_sub_paren) <- c(rep("Epi_1", 7), rep("Epi_2", 7), rep("Epi_3", 7), rep("Epi_4", 7), rep("Fib_1", 7), rep("Fib_2", 7), rep("Fib_3", 7), rep("Endo_1", 7), rep("Capilary", 7), rep("Endo_2", 7), rep("Mural_1", 7), rep("Pericyte", 7), rep("BAM_1", 7), rep("BAM_2", 7), rep("BAM_3", 7), rep("T_Cell", 7), rep("Oligo", 7), rep("Astro", 7), rep("Neuro", 7))

celltype <- factor(paste(colnames(cts_sub_paren)))

options(na.action='na.pass')
Maindm <- model.matrix(~0 + celltype)

rownames(Maindm) <- colnames(cts_sub_paren)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
Maincm <- makeContrasts(levels=Maindm, 
                        Epi_1=celltypeEpi_1-(celltypeEpi_2+celltypeEpi_3+celltypeEpi_4+celltypeFib_1+celltypeFib_2+celltypeFib_3+celltypeFib_4+celltypeEndo+celltypeMural+celltypeBAM_1+celltypeBAM_2+celltypeBAM_3+celltypeT_Cell+celltypeOligo+celltypeAstro+celltypeNeuro)/16, 
                        Epi_2=celltypeEpi_2-(celltypeEpi_1+celltypeEpi_3+celltypeEpi_4+celltypeFib_1+celltypeFib_2+celltypeFib_3+celltypeFib_4+celltypeEndo+celltypeMural+celltypeBAM_1+celltypeBAM_2+celltypeBAM_3+celltypeT_Cell+celltypeOligo+celltypeAstro+celltypeNeuro)/16, 
                        Epi_3=celltypeEpi_3-(celltypeEpi_1+celltypeEpi_2+celltypeEpi_4+celltypeFib_1+celltypeFib_2+celltypeFib_3+celltypeFib_4+celltypeEndo+celltypeMural+celltypeBAM_1+celltypeBAM_2+celltypeBAM_3+celltypeT_Cell+celltypeOligo+celltypeAstro+celltypeNeuro)/16, 
                        Epi_4=celltypeEpi_4-(celltypeEpi_1+celltypeEpi_2+celltypeEpi_3+celltypeFib_1+celltypeFib_2+celltypeFib_3+celltypeFib_4+celltypeEndo+celltypeMural+celltypeBAM_1+celltypeBAM_2+celltypeBAM_3+celltypeT_Cell+celltypeOligo+celltypeAstro+celltypeNeuro)/16, 
                        Fib_1=celltypeFib_1-(celltypeFib_2+celltypeFib_3+celltypeFib_4+celltypeEpi_1+celltypeEpi_2+celltypeEpi_3+celltypeEpi_4+celltypeEndo+celltypeMural+celltypeBAM_1+celltypeBAM_2+celltypeBAM_3+celltypeT_Cell+celltypeOligo+celltypeAstro+celltypeNeuro)/16, 
                        Fib_2=celltypeFib_2-(celltypeFib_1+celltypeFib_3+celltypeFib_4+celltypeEpi_1+celltypeEpi_2+celltypeEpi_3+celltypeEpi_4+celltypeEndo+celltypeMural+celltypeBAM_1+celltypeBAM_2+celltypeBAM_3+celltypeT_Cell+celltypeOligo+celltypeAstro+celltypeNeuro)/16, 
                        Fib_3=celltypeFib_3-(celltypeFib_1+celltypeFib_2+celltypeFib_4+celltypeEpi_1+celltypeEpi_2+celltypeEpi_3+celltypeEpi_4+celltypeEndo+celltypeMural+celltypeBAM_1+celltypeBAM_2+celltypeBAM_3+celltypeT_Cell+celltypeOligo+celltypeAstro+celltypeNeuro)/16, 
                        Fib_4=celltypeFib_4-(celltypeFib_1+celltypeFib_2+celltypeFib_3+celltypeEpi_1+celltypeEpi_2+celltypeEpi_3+celltypeEpi_4+celltypeEndo+celltypeMural+celltypeBAM_1+celltypeBAM_2+celltypeBAM_3+celltypeT_Cell+celltypeOligo+celltypeAstro+celltypeNeuro)/16, 
                        Endo=celltypeEndo-(celltypeEpi_1+celltypeEpi_2+celltypeEpi_3+celltypeEpi_4+celltypeFib_1+celltypeFib_2+celltypeFib_3+celltypeFib_4+celltypeMural+celltypeBAM_1+celltypeBAM_2+celltypeBAM_3+celltypeT_Cell+celltypeOligo+celltypeAstro+celltypeNeuro)/16, 
                        Mural=celltypeMural-(celltypeEpi_1+celltypeEpi_2+celltypeEpi_3+celltypeEpi_4+celltypeFib_1+celltypeFib_2+celltypeFib_3+celltypeFib_4+celltypeEndo+celltypeBAM_1+celltypeBAM_2+celltypeBAM_3+celltypeT_Cell+celltypeOligo+celltypeAstro+celltypeNeuro)/16, 
                        BAM_1=celltypeBAM_1-(celltypeEpi_1+celltypeEpi_2+celltypeEpi_3+celltypeEpi_4+celltypeFib_1+celltypeFib_2+celltypeFib_3+celltypeFib_4+celltypeEndo+celltypeMural+celltypeBAM_2+celltypeBAM_3+celltypeT_Cell+celltypeOligo+celltypeAstro+celltypeNeuro)/16, 
                        BAM_2=celltypeBAM_2-(celltypeEpi_1+celltypeEpi_2+celltypeEpi_3+celltypeEpi_4+celltypeFib_1+celltypeFib_2+celltypeFib_3+celltypeFib_4+celltypeEndo+celltypeMural+celltypeBAM_1+celltypeBAM_3+celltypeT_Cell+celltypeOligo+celltypeAstro+celltypeNeuro)/16, 
                        BAM_3=celltypeBAM_3-(celltypeEpi_1+celltypeEpi_2+celltypeEpi_3+celltypeEpi_4+celltypeFib_1+celltypeFib_2+celltypeFib_3+celltypeFib_4+celltypeEndo+celltypeMural+celltypeBAM_1+celltypeBAM_2+celltypeT_Cell+celltypeOligo+celltypeAstro+celltypeNeuro)/16, 
                        T_Cell=celltypeT_Cell-c(celltypeEpi_1+celltypeEpi_2+celltypeEpi_3+celltypeEpi_4+celltypeFib_1+celltypeFib_2+celltypeFib_3+celltypeFib_4+celltypeEndo+celltypeMural+celltypeBAM_1+celltypeBAM_2+celltypeBAM_3+celltypeOligo+celltypeAstro+celltypeNeuro)/16, 
                        Oligo=celltypeOligo-(celltypeEpi_1+celltypeEpi_2+celltypeEpi_3+celltypeEpi_4+celltypeFib_1+celltypeFib_2+celltypeFib_3+celltypeFib_4+celltypeEndo+celltypeMural+celltypeBAM_1+celltypeBAM_2+celltypeBAM_3+celltypeT_Cell+celltypeAstro+celltypeNeuro)/16, 
                        Neuro=celltypeNeuro-(celltypeEpi_1+celltypeEpi_2+celltypeEpi_3+celltypeEpi_4+celltypeFib_1+celltypeFib_2+celltypeFib_3+celltypeFib_4+celltypeEndo+celltypeMural+celltypeBAM_1+celltypeBAM_2+celltypeBAM_3+celltypeT_Cell+celltypeOligo+celltypeAstro)/16, 
                        Astro=celltypeAstro-(celltypeEpi_1+celltypeEpi_2+celltypeEpi_3+celltypeEpi_4+celltypeFib_1+celltypeFib_2+celltypeFib_3+celltypeFib_4+celltypeEndo+celltypeMural+celltypeBAM_1+celltypeBAM_2+celltypeBAM_3+celltypeT_Cell+celltypeOligo+celltypeNeuro)/16
)

cts_sub_paren <- cts_sub_paren[rowSums(cts_sub_paren)>2*ncol(cts_sub_paren),]
vdat <- voom(cts_sub_paren, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

for(coef in colnames(Mainf)){
  limma_resfun(coef)
}


###Within Epi

cts_sub_paren <- cbind(cts_soupX_Epi_1, cts_soupX_Epi_2, cts_soupX_Epi_3, cts_soupX_Epi_4)
colnames(cts_sub_paren) <- c(rep("Epi_1", 7), rep("Epi_2", 7), rep("Epi_3", 7), rep("Epi_4", 7))

celltype <- factor(paste(colnames(cts_sub_paren)))

options(na.action='na.pass')
Maindm <- model.matrix(~0 + celltype)

rownames(Maindm) <- colnames(cts_sub_paren)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
Maincm <- makeContrasts(levels=Maindm, 
                        Epi_1=celltypeEpi_1-(celltypeEpi_2+celltypeEpi_3+celltypeEpi_4)/3, 
                        Epi_2=celltypeEpi_2-(celltypeEpi_1+celltypeEpi_3+celltypeEpi_4)/3, 
                        Epi_3=celltypeEpi_3-(celltypeEpi_4+celltypeEpi_2+celltypeEpi_4)/3, 
                        Epi_4=celltypeEpi_4-(celltypeEpi_1+celltypeEpi_2+celltypeEpi_3)/3)

cts_sub_paren <- cts_sub_paren[rowSums(cts_sub_paren)>2*ncol(cts_sub_paren),]
vdat <- voom(cts_sub_paren, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

for(coef in colnames(Mainf)){
  limma_resfun(coef)
}


###Within Fib

cts_sub_paren <- cbind(cts_soupX_Fib_1, cts_soupX_Fib_2, cts_soupX_Fib_3)
colnames(cts_sub_paren) <- c(rep("Fib_1", 7), rep("Fib_2", 7), rep("Fib_3", 7))

celltype <- factor(paste(colnames(cts_sub_paren)))

options(na.action='na.pass')
Maindm <- model.matrix(~0 + celltype)

rownames(Maindm) <- colnames(cts_sub_paren)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
Maincm <- makeContrasts(levels=Maindm, 
                        Fib_1=celltypeFib_1-(celltypeFib_2+celltypeFib_3)/2, 
                        Fib_2=celltypeFib_2-(celltypeFib_1+celltypeFib_3)/2, 
                        Fib_3=celltypeFib_3-(celltypeFib_1+celltypeFib_2)/2)

cts_sub_paren <- cts_sub_paren[rowSums(cts_sub_paren)>2*ncol(cts_sub_paren),]
vdat <- voom(cts_sub_paren, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

for(coef in colnames(Mainf)){
  limma_resfun(coef)
}


###Within Endo

cts_sub_paren <- cbind(cts_soupX_Endo_1, cts_soupX_Endo_3, cts_soupX_Endo_2)
colnames(cts_sub_paren) <- c(rep("Endo_1", 7), rep("Endo_3", 7), rep("Endo_2", 7))

celltype <- factor(paste(colnames(cts_sub_paren)))

options(na.action='na.pass')
Maindm <- model.matrix(~0 + celltype)

rownames(Maindm) <- colnames(cts_sub_paren)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
Maincm <- makeContrasts(levels=Maindm, 
                        Endo_1=celltypeEndo_1-(celltypeEndo_3+celltypeEndo_2)/2, 
                        Endo_3=celltypeEndo_3-(celltypeEndo_1+celltypeEndo_2)/2, 
                        Endo_2=celltypeEndo_2-(celltypeEndo_1+celltypeEndo_3)/2)

cts_sub_paren <- cts_sub_paren[rowSums(cts_sub_paren)>2*ncol(cts_sub_paren),]
vdat <- voom(cts_sub_paren, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

for(coef in colnames(Mainf)){
  limma_resfun(coef)
}


###Within Mural

cts_sub_paren <- cbind(cts_soupX_Mural_1, cts_soupX_Mural_2, cts_soupX_Mural_3)
colnames(cts_sub_paren) <- c(rep("Mural_1", 7), rep("Mural_2", 7),  rep("Mural_3", 7))

celltype <- factor(paste(colnames(cts_sub_paren)))

options(na.action='na.pass')
Maindm <- model.matrix(~0 + celltype)

rownames(Maindm) <- colnames(cts_sub_paren)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
Maincm <- makeContrasts(levels=Maindm, 
                        Mural_1=celltypeMural_1-(celltypeMural_2+celltypeMural_3)/2,
                        Mural_2=celltypeMural_2-(celltypeMural_1+celltypeMural_3)/2,
                        Mural_3=celltypeMural_3-(celltypeMural_1+celltypeMural_2)/2)

cts_sub_paren <- cts_sub_paren[rowSums(cts_sub_paren)>2*ncol(cts_sub_paren),]
vdat <- voom(cts_sub_paren, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

for(coef in colnames(Mainf)){
  limma_resfun(coef)
}


###Within Immune

cts_sub_paren <- cbind(cts_soupX_BAM_1, cts_soupX_BAM_2, cts_soupX_BAM_3, cts_soupX_T_Cell)
colnames(cts_sub_paren) <- c(rep("BAM_1", 7), rep("BAM_2", 7), rep("BAM_3", 7), rep("T_Cell", 7))

celltype <- factor(paste(colnames(cts_sub_paren)))

options(na.action='na.pass')
Maindm <- model.matrix(~0 + celltype)

rownames(Maindm) <- colnames(cts_sub_paren)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
Maincm <- makeContrasts(levels=Maindm, 
                        BAM_1=celltypeBAM_1-(celltypeBAM_2+celltypeBAM_3+celltypeT_Cell)/3, 
                        BAM_2=celltypeBAM_2-(celltypeBAM_1+celltypeBAM_3+celltypeT_Cell)/3, 
                        BAM_3=celltypeBAM_3-(celltypeBAM_1+celltypeBAM_2+celltypeT_Cell)/3, 
                        T_Cell=celltypeT_Cell-(celltypeBAM_1+celltypeBAM_2+celltypeBAM_3)/3)

cts_sub_paren <- cts_sub_paren[rowSums(cts_sub_paren)>2*ncol(cts_sub_paren),]
vdat <- voom(cts_sub_paren, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

for(coef in colnames(Mainf)){
  limma_resfun(coef)
}


###Majorcelltype

###Heatmap by pseudobulk
cts_fun <- function(celltype){
  sub <- subset(Ex_ChP, Majorcelltype==celltype)
  cts <- AggregateExpression(sub, group.by = "library", assay = "soupXcounts", slot = "counts", return.seurat = FALSE)
  cts <- cts[["soupXcounts"]]
  assign(c(paste("cts_soupX", celltype, sep = "_")), cts, envir = .GlobalEnv)
  pct <- pctincell(sub)
  assign(c(paste("pct", celltype, sep = "_")), pct, envir = .GlobalEnv)
}

cts_fun("Epithelial")
gc()
cts_fun("Fibroblast")
gc()
cts_fun("Endothelial")
cts_fun("Mural")
cts_fun("Immune")
cts_fun("Parenchyma")

dim(cts_soupX_Epithelial)
dim(cts_soupX_Fibroblast)
dim(cts_soupX_Endothelial)
dim(cts_soupX_Mural)
dim(cts_soupX_Immune)
dim(cts_soupX_Parenchyma)

cts_major_paren <- cbind(cts_soupX_Epithelial, cts_soupX_Fibroblast, cts_soupX_Endothelial, cts_soupX_Mural, cts_soupX_Immune, cts_soupX_Parenchyma)
colnames(cts_major_paren) <- c(rep("Epithelial", 7), rep("Fibroblast", 7), rep("Endothelial", 7), rep("Mural", 7), rep("Immune", 7), rep("Parenchyma", 7))

celltype <- factor(paste(colnames(cts_major_paren)))

options(na.action='na.pass')
Maindm <- model.matrix(~0 + celltype)

rownames(Maindm) <- colnames(cts_major_paren)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
cts_major_paren <- cts_major_paren[,colnames(cts_major_paren)%in%rownames(Maindm)]
cts_major_paren <- cts_major_paren[,colnames(cts_major_paren)%in%rownames(Maindm)]
Maincm <- makeContrasts(levels=Maindm, Epithelial=celltypeEpithelial-(celltypeFibroblast+celltypeEndothelial+celltypeMural+celltypeImmune+celltypeParenchyma)/5, Fibroblast=celltypeFibroblast-(celltypeEpithelial+celltypeEndothelial+celltypeMural+celltypeImmune+celltypeParenchyma)/5, Endothelial=celltypeEndothelial-(celltypeEpithelial+celltypeFibroblast+celltypeMural+celltypeImmune+celltypeParenchyma)/5, Mural=celltypeMural-(celltypeEpithelial+celltypeFibroblast+celltypeEndothelial+celltypeImmune+celltypeParenchyma)/5, Immune=celltypeImmune-(celltypeEpithelial+celltypeFibroblast+celltypeEndothelial+celltypeMural+celltypeParenchyma)/5, Parenchyma=celltypeParenchyma-(celltypeEpithelial+celltypeFibroblast+celltypeEndothelial+celltypeMural+celltypeImmune)/5)

cts_major_paren <- cts_major_paren[rowSums(cts_major_paren)>2*ncol(cts_major_paren),]
vdat <- voom(cts_major_paren, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

for(coef in colnames(Mainf)){
  limma_resfun(coef)
}



###In sn
####Load ST

cos_Epi_1_sig <- read_csv("/Ex_ST/Epi_1/Spreadsheet__Epi_1_corrected.csv")
cos_Epi_2_sig <- read_csv("/Ex_ST/Epi_2/Spreadsheet__Epi_2_corrected.csv")
cos_Epi_3_sig <- read_csv("/Ex_ST/Epi_3/Spreadsheet__Epi_3_corrected.csv")
cos_Epi_4_sig <- read_csv("/Ex_ST/Epi_4/Spreadsheet__Epi_4_corrected.csv")
cos_Fib_1_sig <- read_csv("/Ex_ST/Fib_1/Spreadsheet__Fib_1_corrected.csv")
cos_Fib_2_sig <- read_csv("/Ex_ST/Fib_2/Spreadsheet__Fib_2_corrected.csv")
cos_Fib_3_sig <- read_csv("/Ex_ST/Fib_3/Spreadsheet__Fib_3_corrected.csv")
cos_Fib_4_sig <- read_csv("/Ex_ST/Fib_4/Spreadsheet__Fib_4_corrected.csv")
cos_Mural_sig <- read_csv("/Ex_ST/Mural/Spreadsheet__Mural_corrected.csv")
cos_Endo_sig <- read_csv("/Ex_ST/Endo/Spreadsheet__Endo_corrected.csv")
cos_Oligo_sig <- read_csv("/Ex_ST/Olig/Spreadsheet__Olig_corrected.csv")
cos_BAM_1_sig <- read_csv("/Ex_ST/BAM_1/Spreadsheet__BAM_1_corrected.csv")
cos_BAM_2_sig <- read_csv("/Ex_ST/BAM_2/Spreadsheet__BAM_2_corrected.csv")

cos_Epithelial_sig <- read_csv("/Ex_ST/Epithelial/Spreadsheet__Epithelial_corrected.csv")
cos_Fibroblast_sig <- read_csv("/Ex_ST/Fibroblast/Spreadsheet__Fibroblast_corrected.csv")
cos_Endothelial_sig <- read_csv("/Ex_ST/Endothelial/Spreadsheet__Endothelial_corrected.csv")
cos_Immune_sig <- read_csv("/Ex_ST/Immune/Spreadsheet__Immune_corrected.csv")
cos_Parenchyma_sig <- read_csv("/Ex_ST/Parenchyma/Spreadsheet__Parenchyma_corrected.csv")

notinCosMX <- readRDS("/Ex_ST/genes_not_in_CosMX.RDS")


####Overlap venns

pairs <- list(HG38 = c("AARDS", "ACPP", "ATP5MPL", "C9orf16", "CLECL1", "MT-CO1", "MT-CO2", "DPAGT1", "H2AFZ", "HIST1H2BD", "H3F3A", "H3F3C", "HIST1H3H", "HIST1H3B", "HIST1H3G", "HIST1H4C", "HLA-DRB1", "LRMP", "KLRC1", "MYLPF", "CD3EAP", "DDX58", "WDR61", "TMEM173", "WARS", "YARS", "CCDC130"),
              CosMX2 = c("AARS1", "ACP3", "ATP5MJ", "BBLN", "CLECL1P", "COX1", "COX2", "H2AX", "H2AZ1", "H2BC5", "H3-3A", "H3-5", "H3C10", "H3C2", "H3C8", "H4C3", "HLA-DRB", "IRAG2", "KLRC1/2", "MYL11", "POLR1G", "RIGI", "SKIC8", "STING1", "WARS1", "YARS1", "YJU2B"))
sn_cosmx_overlap_fun(celltypeDEGs_Epi_1_sig, cos_Epi_1_sig)
sn_cosmx_overlap_fun(celltypeDEGs_Epi_2_sig, cos_Epi_2_sig)
sn_cosmx_overlap_fun(celltypeDEGs_Epi_3_sig, cos_Epi_3_sig)
sn_cosmx_overlap_fun(celltypeDEGs_Epi_4_sig, cos_Epi_4_sig)
sn_cosmx_overlap_fun(celltypeDEGs_Fib_1_sig, cos_Fib_1_sig)
sn_cosmx_overlap_fun(celltypeDEGs_Fib_2_sig, cos_Fib_2_sig)
sn_cosmx_overlap_fun(celltypeDEGs_Fib_3_sig, cos_Fib_3_sig)
sn_cosmx_overlap_fun(celltypeDEGs_Fib_4_sig, cos_Fib_4_sig)
sn_cosmx_overlap_fun(celltypeDEGs_Endo_sig, cos_Endo_sig)
sn_cosmx_overlap_fun(celltypeDEGs_BAM_1_sig, cos_BAM_1_sig)
sn_cosmx_overlap_fun(celltypeDEGs_BAM_2_sig, cos_BAM_2_sig)
sn_cosmx_overlap_fun(celltypeDEGs_Oligo_sig, cos_Oligo_sig)

sn_cosmx_overlap_fun(celltypeDEGs_Epithelial_sig, cos_Epithelial_sig)
sn_cosmx_overlap_fun(celltypeDEGs_Fibroblast_sig, cos_Fibroblast_sig)
sn_cosmx_overlap_fun(celltypeDEGs_Mural_sig, cos_Mural_sig)
sn_cosmx_overlap_fun(celltypeDEGs_Endothelial_sig, cos_Endothelial_sig)
sn_cosmx_overlap_fun(celltypeDEGs_Immune_sig, cos_Immune_sig)
sn_cosmx_overlap_fun(celltypeDEGs_Parenchyma_sig, cos_Parenchyma_sig)


###FGSEA celltype
Path_fun(celltypeDEGs_Epi_1)
Path_fun(celltypeDEGs_Epi_2)
Path_fun(celltypeDEGs_Epi_3)
Path_fun(celltypeDEGs_Epi_4)
Path_fun(celltypeDEGs_Fib_1)
Path_fun(celltypeDEGs_Fib_2)
Path_fun(celltypeDEGs_Fib_3)
Path_fun(celltypeDEGs_Mural)
Path_fun(celltypeDEGs_Endo_1)
Path_fun(celltypeDEGs_Endo_3)
Path_fun(celltypeDEGs_Endo_2)
Path_fun(celltypeDEGs_Mural_2)
Path_fun(celltypeDEGs_Mural_3)
Path_fun(celltypeDEGs_Mural_1)
Path_fun(celltypeDEGs_BAM_1)
Path_fun(celltypeDEGs_BAM_2)
Path_fun(celltypeDEGs_BAM_3)
Path_fun(celltypeDEGs_T_Cell)
Path_fun(celltypeDEGs_Oligo)
Path_fun(celltypeDEGs_Neuro)
Path_fun(celltypeDEGs_Astro)

Path_fun(celltypeDEGs_Epithelial)
Path_fun(celltypeDEGs_Fibroblast)
Path_fun(celltypeDEGs_Endothelial)
Path_fun(celltypeDEGs_Mural)
Path_fun(celltypeDEGs_Immune)
Path_fun(celltypeDEGs_Parenchyma)

##Figures Celltype
####UMAP

#Major types
Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
unique(Ex_ChP@meta.data$Subcelltype)
Ex_ChP <- RenameIdents(Ex_ChP, "Epi_1a"="Epithelial", "Epi_1b"="Epithelial", "Epi_2a"="Epithelial", "Epi_2b"="Epithelial", "Fib_1"="Fibroblast", "Fib_2"="Fibroblast", "Fib_3"="Fibroblast", "Endo_1"="Endothelial", "Endo_2"="Endothelial", "Endo_3"="Endothelial", "Mural_2"="Mural", "Mural_3"="Mural", "Mural_1"="Mural", "Oligo"="Parenchyma", "Astro"="Parenchyma", "Neuro"="Parenchyma", "BAM_1"="Immune", "BAM_2"="Immune", "BAM_3"="Immune", "T_Cell"="Immune")
levels(Ex_ChP) <- c("Epithelial", "Fibroblast", "Endothelial", "Mural", "Immune", "Parenchyma")
Ex_ChP@meta.data[["Majorcelltype"]] <- Idents(Ex_ChP)

plts <- DimPlot(Ex_ChP, raster=F,raster.dpi=c(600,600), cols=c("Epithelial"="#C48F45", "Fibroblast"="#3249A6", "Endothelial"="#926392", "Mural"="#261132", "Immune"="#26532B", "Parenchyma"="sienna")) + theme(text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Ex_ChP_UMAP_major.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Ex_ChP_UMAP_major.png", sep=""), width=10, height=5)
remove(plts)

#Subcelltypes
Idents(Ex_ChP) <- Ex_ChP@meta.data[["Subcelltype"]] 
levels(Ex_ChP) <- c("Epi_1a", "Epi_1b", "Epi_2a", "Epi_2b", "Fib_1", "Fib_2", "Fib_3", "Endo_1", "Endo_2", "Endo_3", "Mural_1", "Mural_2", "Mural_3", "BAM_1", "BAM_2", "BAM_3", "T_Cell", "Oligo", "Astro", "Neuro")
Ex_ChP@meta.data[["Subcelltype"]] <- Idents(Ex_ChP)

plts <- DimPlot(Ex_ChP, raster=F, raster.dpi=c(600,600), cols=c("Epi_1a"="#CD8500", "Epi_1b"="#FFA500",  "Epi_2a"="#FF4500", "Epi_2b"="#CD3700", "Fib_1"="#050382", "Fib_2"="#27D7EF", "Fib_3"="#2a91ea", "Endo_1"="#c875c4", "Endo_2"="#976697", "Endo_3"="#c5aac5", "Mural_1"="#261132", "Mural_2"="#8f2aa2", "Mural_3"="#6d207b", "Oligo"="lightsalmon", "Astro"="salmon2", "Neuro"="sienna", "BAM_1"="#70C25B", "BAM_2"="#0E9554", "BAM_3"="#266417", "T_Cell"="#3ded97"))  + theme(text=element_text(size=7))#sanity check and make publishable figs
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Ex_ChP_UMAP.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Ex_ChP_UMAP.png", sep=""), width=10, height=5)
remove(plts)

Ex_ChP <- subset(Ex_ChP, Majorcelltype!="Parenchyma")
Ex_ChP <- RunUMAP(object = Ex_ChP, dims = 1:40, assay = "soupXcounts", reduction =  "harmony", n.neighbors=50, min.dist=0.8)
Idents(Ex_ChP) <- Ex_ChP@meta.data[["Subcelltype"]] 
plts <- DimPlot(Ex_ChP, raster=T,raster.dpi=c(600,600), cols=c("Epi_1a"="#CD8500", "Epi_1b"="#FFA500",  "Epi_2a"="#FF4500", "Epi_2b"="#CD3700", "Fib_1"="#050382", "Fib_2"="#27D7EF", "Fib_3"="#2a91ea", "Endo_1"="#c875c4", "Endo_2"="#976697", "Endo_3"="#c5aac5", "Mural_1"="#261132", "Mural_2"="#8f2aa2", "Mural_3"="#6d207b", "BAM_1"="#70C25B", "BAM_2"="#0E9554", "BAM_3"="#266417", "T_Cell"="#3ded97"))  + theme(text=element_text(size=7))#sanity check and make publishable figs
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Ex_ChP_UMAP_noParenchyma.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Ex_ChP_UMAP_noParenchyma.png", sep=""), width=10, height=5)
remove(plts)

#Major types
Idents(Ex_ChP) <- Ex_ChP@meta.data$Majorcelltype
plts <- DimPlot(Ex_ChP, raster=T,raster.dpi=c(600,600), cols=c("Epithelial"="#C48F45", "Fibroblast"="#3249A6", "Endothelial"="#926392", "Mural"="#261132", "Immune"="#26532B")) + theme(text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Ex_ChP_UMAP_major_noParenchyma.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Ex_ChP_UMAP_major_noParenchyma.png", sep=""), width=10, height=5)
remove(plts)

#Sub type individual figures
Idents(Epi) <- Epi@meta.data[["Subcelltype"]]
plts <- DimPlot(Epi, raster=T, cols=c("Epi_1a"="#CD8500", "Epi_1b"="#FFA500",  "Epi_2a"="#FF4500", "Epi_2b"="#CD3700")) + theme(text=element_text(size=7)) #make publishable
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Epi_UMAP.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Epi_UMAP.png", sep=""), width=10, height=5)
remove(plts)

Idents(Fibro) <- Fibro@meta.data[["Subcelltype"]]
plts <- DimPlot(Fibro, raster=T, cols=c("Fib_1"="#050382", "Fib_2"="#27D7EF", "Fib_3"="#2a91ea")) + theme(text=element_text(size=7)) #make publishable figs
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Fibro_UMAP.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Fibro_UMAP.png", sep=""), width=10, height=5)
remove(plts)

Idents(Imm) <- Imm@meta.data[["Subcelltype"]]
levels(Imm) <- c("BAM_1", "BAM_2", "BAM_3", "T_Cell")
Imm@meta.data[["Subcelltype"]] <- Idents(Imm)
#Imm <- RunUMAP(Imm, dims=1:20, assay="soupXcounts", reduction="harmony", min.dist = .8)
plts <- DimPlot(Imm, raster=T, cols=c("BAM_1"="#70C25B", "BAM_2"="#0E9554", "BAM_3"="#266417", "T_Cell"="#3ded97")) + theme(text=element_text(size=7)) # make publishable figs
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Imm_UMAP.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Imm_UMAP.png", sep=""), width=10, height=5)
remove(plts)

Idents(Endo) <- Endo@meta.data[["Subcelltype"]]
levels(Endo) <- c("Endo_1", "Endo_3", "Endo_2")
Endo@meta.data[["Subcelltype"]] <- Idents(Endo)
#Endo <- RunUMAP(Endo, dims=1:20, assay="soupXcounts", reduction="harmony", min.dist = 1)
plts <- DimPlot(Endo, raster=F, cols=c("Endo_1"="#c875c4", "Endo_2"="#976697", "Endo_3"="#c5aac5")) + theme(text=element_text(size=7)) # make publishable figs
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Endo_UMAP.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Endo_UMAP.png", sep=""), width=10, height=5)
remove(plts)

Idents(Mural) <- Mural@meta.data[["Subcelltype"]]
levels(Mural) <- c("Mural_2", "Mural_3", "Mural_1")
Mural@meta.data[["Subcelltype"]] <- Idents(Mural)
#Mural <- RunUMAP(Mural, dims=1:20, assay="soupXcounts", reduction="harmony", min.dist = 1.4)
plts <- DimPlot(Mural, raster=F, cols=c("Mural_1"="#261132", "Mural_2"="#8f2aa2", "Mural_3"="#6d207b")) + theme(text=element_text(size=7)) # make publishable figs
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Mural_UMAP.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Mural_UMAP.png", sep=""), width=10, height=5)
remove(plts)

rm(Fibro_meta, Endo_meta, Neuro_meta, Epi_1_meta, Epi_2_meta, BAM_meta, Astro_meta, T_meta, Oligo_meta, Imm_meta)
gc()


###Heatmap

library(ComplexHeatmap)
library(circlize)
library(DESeq2)
library(viridis)
gc()

#select list
Epi1a_genes <- c("CDH12", "CDH1", "SLC26A3", "ADAM10", "MTUS2", "TMEM108", "MBD5", "GALNT17", "ARHGEF3", "ZBTB20", "NRG3", "MAGI2", "MAGI1", "MAPK10", "TBCK", "PTPRM", "KIAA1217", "BBS9", "WDPCP", "AFF3", "MVB12B", "FTO", "DOCK3", "DNM3", "ZNF521", "NBEA", "PLCB4", "DANT2", "ITGA6", "TMEM232", "ESRRG", "PPP1R9A", "TMEM117", "PHACTR2", "ATXN1", "SNTB1", "PDZRN3", "PCAT1", "ATP13A5", "DYM", "XYLT1", "LDLRAD4", "STK33", "GLIS3", "EYA4", "NFATC2", "TLN2", "ROR2") 
Epi1b_genes <- c("CLIC6", "FOXJ1", "DNAH12", "DNAH9", "DNAH6", "DNAH7", "DNAH11", "RMST", "ROR2", "RFX3", "ARL13B", "SPAG17", "TMEM72-AS1", "TMEM117", "RYR3", "SLC24A4", "SNTB1", "SLC12A2", "SLC39A12", "EYA4", "F5", "RIMBP2", "FYB2", "SLC4A10", "CFAP43", "CRPPA", "OCA2", "GNA14", "NEK11", "CFAP54", "MTUS2", "CFAP44", "DGKI", "GPR155", "DOCK8", "POU2AF1", "PACRG", "HYDIN", "PRDM16", "SHANK2", "CFAP100", "CFAP69", "TUBB4B", "CFAP61", "PTPRG", "DNAI1", "CFAP299", "DNAH9", "DNAH7", "DYNC1I1", "FAM120B", "FAM160A1", "ESRRG", "EPHA4", "LRRIQ1", "CHST11", "COL9A1", "CTNND2", "CUX2", "DCDC1", "LRRC9", "MRPS6", "SDK1", "MYO1E", "SLC15A4", "SPAG16", "HTR2C", "OTX2", "PRLR", "ABCA4", "ABCA5", "GALNT11", "TPD52L1", "TRPV4", "SGK1", "NWD1", "CDH12", "FBLN1", "OTX2-AS1", "TRMT9B", "VAT1L", "WDR49", "ZBBX", "LMX1A")
Epi2_genes <- c("CLDN1", "OCLN", "COX2", "MT-ND4", "MT-CO3", "MT-CO2", "MT-CO1", "MAGI1", "MT-ND3", "BBS9", "SPATA17", "GPX3", "MT-ATP6", "MT-ATP8", "MT-CYB", "MT-ND1", "MT-ND2", "SLC26A3", "MT-ND3", "MT-ND4L", "MT-ND5", "MTRNR2L8", "MT-ND1", "GALNT17", "AQP1", "ATP1A1", "SLC25A6", "SLC2A12", "SLC25A23", "SLC4A2", "FOLR1", "IGFBP2", "EEF1A1", "TMBIM6", "RPL13", "MTRNR2L12", "PDK4", "TTR")
Fib_genes <- c("PPM1H", "PDE7B", "PBX1", "NEBL", "ITGA9", "KLF9", "LRP1", "LRP1B", "PLCB1", "SESTD1", "SETBP1", "MID1", "GPC6", "MIR100HG", "GMDS", "PHLDB2", "EGFR", "EYA2", "FAM20A", "FOXP2", "ADAMTSL3", "ARHGAP10", "EGFR", "ATP13A3", "ARHGAP21", "CDH11", "ADK", "SLC22A23", "BCL6", "THSD4","BICC1", "BNC2", "CHSY3", "SLC7A2", "ATP1B3", "BMPR1B", "AKAP12", "ADAM12", "ADD3",  "ALPL", "ARHGAP6", "ATP13A3", "WWTR1", "TMOD1", "TMTC1", "USP53", "CDK14","EPS8", "GLCE", "GNAL", "LPAR1", "MAP3K5", "MAST4", "NID1", "NR2F1", "NR2F1-AS1", "SLC1A1",  "SLC38A2", "SMOC1", "ACADL", "JAK2", "ALDH1A2", "LEPR", "SLC24A3", "STXBP6", "TBC1D8")
Fib1_genes <- c("PLXNA4", "KCNMA1", "CACNA2D3", "SPTLC3", "SLC4A4", "ALPL", "ALCAM", "UNC5C")
Fib2_genes <- c("TGFBR3", "ABCA9", "APOD", "EDA", "EBF1", "KCNT2", "ITGA8", "ROBO1", "TRPM3", "MYRIP", "GPM6A", "CLDN11", "TJP1", "TJP2", "ITGBL1", "GJB6", "PTGDS", "ECM2", "SLC26A20", "SLC13A1", "LAMA4", "COCH", "NTM", "COLEC12", "HIVEP2", "LAMA4", "LDLRAD3", "SLC22A3")
Fib3_genes <- c("COL5A1", "COL5A2", "GRK5", "ABCA8", "SPAG9", "DCN", "FN1", "ABC9A9", "LAMA2", "COL4A2", "COL3A1", "COL6A3")
Endo_genes <- c("CCNY", "CLIC4", "ARHGAP26", "ARHGAP29", "ATP13A3", "PODXL", "PTPRB", "SHANK3", "RPGR", "CMIP", "EPAS1", "IQSEC1", "MSN", "IVNS1ABP", "GALNT18", "PLEKHG1", "PLPP1", "CRIM1", "FLT1", "BMPR2", "GRB10", "DOCK9", "UBC", "INSR", "LDB2", "NFIB", "ARL15", "PIK3R3", "STGALNAC3", "UGCG", "ADGRL4", 
                "PECAM1", "VWF", "INSR", "CLDN5", "CDH13", "ARL15", "MECOM", "VEGFC", "COL8A1", "FBN1", "THSD7A", 
                "TSHZ2", "CEMIP2", "IL1R1", "IL4R", "FLI1", "SORBS2", 
                "PLVAP", "ADGRF5", "KCNQ3", "SLC2A3", "ITGA1")
Mural_genes <- c("RBPMS", "RCAN2", "RGS6", "CARMN", "CALD1", "SYNPO2", "DLC1", "SGIP1", "MCAM", "LTBP1", "NOTCH3", "ADAMTS1", "FRY", "JAG1", "COL4A1", "COL4A2", "EBF1", "EPS8", "LPP", "SLIT3",
                 "SOX6", "MYO1D", "MYH11", "ACTA2", "PRUNE2", "RYR2", "TPM1", "LMOD1", "COL18A1", "SCL22A3",
                 "ENPP2", "CEMIP", "GLIS3", "PDZD2", 
                 "COL4A4", "CRISPLD2", "ADAMTS9", "TLE1", "PIPNC1", "ADCY3", "HIPK2")
BAM_genes <- c("STAB1", "MERTK", "CD163", "AIF1", "CSF1R", "APBB1IP", "C1QA", "C1QB", "C1QC", "CLEC7A", "ENPP2", "MERTK", "MS4A4E", "MS4A6A", "MS4A7", "PTPRC")
BAM1_genes <- c("ABCA1", "AOAH", "STARD13", "ADAM28", "MEF2A", "MERTK", "MGAT4A", "MRC1", "PDE4B", "PLCG2", "RAB31", "SFMBT2", "SLC1A3", "SPRED1", "SRGN", "F13A1", "LYVE1", "HCK", "IL4R", "CD74", "HLA-DRA", "CD86", "PRKN", "SNX6", "ARHGAP24", "RHBDF2", "WDFY3", "RIN2", "MSR1", "MEF2C", "STK3", "SRGAP2", "EPS8", "ITSN1", "EIF4G3", "FOXO3", "NUMB", "ATF6", "IGSF21", "MAN1A1", "FTX", "FER", "TLR2", "PTPRE", "ZSWIM6", "ZFP36L1", "MAML3", "ATP8B4", "FCHSD2", "SAT1", "BACH1", "TBC1D14", "ZNF710", "AGFG1", "LPAR6", "LYN", "ANKS1A", "SPRED1", "SRGAP1", "GAB2", "SUSD6", "RB1", "ARHGAP18", "CSGALNACT1", "TNRC18", "ST3GAL1", "VMP1")
BAM2_genes <- c("DNAJB6", "HSP90AA1", "HSPH1", "HSPB1", "DNAJB1", "HSPD1", "HSPE1", "DNAJA4", "HSP90AB1", "HSPA1A", "HSPA1B", "HSPA5", "HSPA6", "HSPA8", "HSPB1", "ATP2C1", "CD83", "CHORDC1", "P4HA1", "PTGES3", "SPECC1", "SPP1", "TRA2B", "UBC", "UGCG", "USP36", "USPL1")
BAM3_genes <- c("MKI67", "CENF", "TOP2A", "ITGAM", "CCR2", "ANLN", "ANXA2", "APOLD1", "ARHGAP11", "ARHGAP1B", "ASPM", "ATAD2", "ATAD5", "BARD1", "BRCA1", "BRCA2", "BRIP1", "BUB1", "BUB1B", "CCDC18", "CIP2A", "CIT", "CKAP2", "CKAP2L", "CLSPN", "DEPDC1B", "DIAPH3", "DTL", "ECT2", "EZH2", "FAM111B", "FANCI", "GINS1", "GPSM2", "GTSE1", "HJURP", "HMMR", "KIF4A", "KIF11", "KIF14", "KIF15", "KIF18A", "KIF18B", "KIF20B", "KIF23", "KIFC1", "KNL1", "KNTC1", "LMNB1", "LNCAROD", "MELK", "NCAPG", "NCAPG2", "MCAPH", "NUF2", "NUSAP1", "PARPBP", "POLQ", "PRR1", "PVT1", "RFC3", "RRM2", "SCLT1", "SGO2", "SMC4", "TPX2", "VRK1",  "CEP128", "ATAD2", "CEP152", "CENPP")
Tcell_genes <- c("THEMIS", "IL7R", "CANK4", "CBLB", "BCL11B", "ANK3", "ETS1", "FYN", "GLCCI1", "INPP4B", "IQGAP2", "ITGA4", "ITK", "KLF12", "PARP8", "PRKCH", "PPP1R16B", "PPP2R2B", "PRKACB", "PTPN22", "RASA3", "RESF1", "SKAP1", "SLC38A1", "SLFN12L", "STAT4", "SYNE2", "SYTL3", "TC2N", "TOX", "TRERF1", "TSPAN5", "ZEB1", "ZNF831")
paren_genes <- c("FGF14", "ADGRB3", "ARPP21", "ANKS1B", "GRIK2", "MEG3", "FGF12", "GRIA2", "RIMS2", "KCND2", "SNTG1", "TRIM9", "RALYL", "NRXN1", "PPFIA2", "RBFOX1", "LRRC7", "MYT1L", "SEZ6L", "TMEM178B", "GRIA4", "GRIA1", "KCNIP4", "SNAP25", "PHACTR1", "RYR2", "SLC35F3", "CACNA1B", "STXBP5L", "ANK3", "LRRTM4", "NEGR1", "ATP8A2", "CCSER1", "NRCAM", "TMEM132B", "PCSK2", "GRM5", "DLGAP1", "IQCJ-SCHIP1", "DGKB", "CELF4", "SYN3", "CADPS", "SRRM3", "CNTN1", "MAP2", "CSMD3", "KCNQ5", "XKR4", "GRIN2A", "RBFOX3", "MAPT", "ADGRL3", "UNC79", "PCLO", "CNTNAP2", "SYT1", "ERC2", "GABRB3", "CNTNAP5", "HECW1", "GALNT13", "DLGAP2", "PLCB1", "RIMS1", "FAM155A", "CACNA1D", "DCLK1", "CSMD1", "ASIC2", "PTPRD", "DSCAM", "MACROD2", "RGS7", "GRIN2B", "SPTBN4", "KIF5C", "NKAIN2", "LRRC4C", "GRM7", "SNHG14", "NCAM2", "GRID2", "KIRREL3", "DLG2", "OPCML", "MAP1B", "PRKCB", "KSR2", "EPHA6", "DAB1", "FRY", "PCDH9", "NRXN3", "RBFOX1", "NRG3", "CADM2")


###Heatmap by pseudobulk
cts_fun <- function(celltype){
  sub <- subset(Ex_ChP, Subcelltype==celltype)
  cts <- AggregateExpression(sub, group.by = "projid", assay = "soupXcounts", slot = "counts", return.seurat = FALSE)
  cts <- cts[["soupXcounts"]]
  assign(c(paste("cts", "soupX", celltype, sep = "_")), cts, envir = .GlobalEnv)
  gc()
}

cts_fun("Epi_1a")
cts_fun("Epi_1b")
cts_fun("Epi_2a")
cts_fun("Epi_2b")
cts_fun("Fib_1")
cts_fun("Fib_2")
cts_fun("Fib_3")
cts_fun("Endo_1")
cts_fun("Endo_2")
cts_fun("Endo_3")
cts_fun("Mural_1")
cts_fun("Mural_2")
cts_fun("Mural_3")
cts_fun("BAM_1")
cts_fun("BAM_2")
cts_fun("BAM_3")
cts_fun("T_Cell")


cts_fun <- function(celltype){
  sub <- subset(Ex_ChP, Majorcelltype==celltype)
  cts <- AggregateExpression(sub, group.by = "projid", assay = "soupXcounts", slot = "counts", return.seurat = FALSE)
  cts <- cts[["soupXcounts"]]
  assign(c(paste("cts", "soupX", celltype, sep = "_")), cts, envir = .GlobalEnv)
  gc()
}

cts_fun("Epithelial")
cts_fun("Fibroblast")
cts_fun("Endothelial")
cts_fun("Mural")
cts_fun("Immune")
cts_fun("Parenchyma")

#sub no paren
my_genes <- c(Epi1a_genes, Epi1b_genes, Epi_genes, Epi2_genes, Fib_genes, Fib1_genes, Fib2_genes, Fib3_genes, Endo_genes, Mural_genes, BAM_genes, BAM1_genes, BAM2_genes, BAM3_genes, Tcell_genes) %>% unique()
#my_genes <- my_genes[my_genes %in% rownames(Ex_ChP)]
cts_sub <- cbind(cts_soupX_Epi_1a, cts_soupX_Epi_1b, cts_soupX_Epi_2a, cts_soupX_Epi_2b, cts_soupX_Fib_1, cts_soupX_Fib_2, cts_soupX_Fib_3, cts_soupX_Endo_1, cts_soupX_Endo_2, cts_soupX_Endo_3, cts_soupX_Mural_1, cts_soupX_Mural_2, cts_soupX_Mural_3, cts_soupX_BAM_1, cts_soupX_BAM_2, cts_soupX_BAM_3, cts_soupX_T_Cell)

cts_sub <- data.frame(scale(cts_sub))
cts_sub <- cts_sub[rownames(cts_sub) %in% my_genes,]
cts_sub <- cts_sub[order(match(rownames(cts_sub), my_genes)),]
cts_sub <- t(scale(t(cts_sub)))
htg <- Heatmap(cts_sub, heatmap_legend_param = list(title = "Expression"), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, show_row_names = TRUE, col = colorRamp2(c(-2, 0, 2, 4), c("#000004FF","#721F81FF", "#F1605DFF", "#FCFDBFFF")))
draw(htg)
graph2png(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_batch_noparen.png", sep = ""), width=7, height=7)
graph2pdf(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_batch_noparen.pdf", sep = ""), width=7, height=7)

#major no paren
my_genes <- c(Epi1a_genes, Epi1b_genes, Epi_genes, Epi2_genes, Fib_genes, Fib1_genes, Fib2_genes, Fib3_genes, Endo_genes, Mural_genes, BAM_genes, BAM1_genes) %>% unique()
#my_genes <- my_genes[my_genes %in% rownames(Ex_ChP)]
cts_mjr <- cbind(cts_soupX_Epithelial, cts_soupX_Fibroblast, cts_soupX_Endothelial, cts_soupX_Mural, cts_soupX_Immune)
cts_mjr <- data.frame(scale(cts_mjr))
cts_mjr <- cts_mjr[rownames(cts_mjr) %in% my_genes,] #!!!!no parenchyma in my genes
cts_mjr <- cts_mjr[order(match(rownames(cts_mjr), my_genes)),]
cts_mjr <- t(scale(t(cts_mjr)))
htg <- Heatmap(cts_mjr, heatmap_legend_param = list(title = "Expression"), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, show_row_names = TRUE, col = colorRamp2(c(-2, 0, 2, 4), c("#000004FF","#721F81FF", "#F1605DFF", "#FCFDBFFF")))
draw(htg)
graph2png(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_batch_major_noparen.png", sep = ""), width=7, height=7)
graph2pdf(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_batch_major_noparen.pdf", sep = ""), width=7, height=7)


#sub with paren
my_genes <- c(Epi1a_genes, Epi1b_genes, Epi_genes, Epi2_genes, Fib_genes, Fib1_genes, Fib2_genes, Fib3_genes, Endo_genes, Mural_genes, BAM_genes, BAM1_genes, BAM2_genes, BAM3_genes, Tcell_genes, paren_genes) %>% unique()
#my_genes <- my_genes[my_genes %in% rownames(Ex_ChP)]
cts_sub <- cbind(cts_soupX_Epi_1a, cts_soupX_Epi_1b, cts_soupX_Epi_2a, cts_soupX_Epi_2b, cts_soupX_Fib_1, cts_soupX_Fib_2, cts_soupX_Fib_3, cts_soupX_Endo_1, cts_soupX_Endo_2, cts_soupX_Endo_3, cts_soupX_Mural_1, cts_soupX_Mural_2, cts_soupX_Mural_3, cts_soupX_BAM_1, cts_soupX_BAM_2, cts_soupX_BAM_3, cts_soupX_T_Cell, cts_soupX_Parenchyma)
cts_sub <- data.frame(scale(cts_sub))
cts_sub <- cts_sub[rownames(cts_sub) %in% my_genes,] #!!!!no parenchyma in my genes
cts_sub <- cts_sub[order(match(rownames(cts_sub), my_genes)),]
cts_sub <- t(scale(t(cts_sub)))
htg <- Heatmap(cts_sub, heatmap_legend_param = list(title = "Expression"), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, show_row_names = TRUE, col = colorRamp2(c(-2, 0, 2, 4), c("#000004FF","#721F81FF", "#F1605DFF", "#FCFDBFFF")))
draw(htg)
graph2png(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_batch.png", sep = ""), width=7, height=7)
graph2pdf(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_batch.pdf", sep = ""), width=7, height=7)

#major with paren
my_genes <- c(Epi1a_genes, Epi1b_genes, Epi_genes, Epi2_genes, Fib_genes, Fib1_genes, Fib2_genes, Fib3_genes, Endo_genes, Mural_genes, BAM_genes, BAM1_genes, paren_genes) %>% unique()
#my_genes <- my_genes[my_genes %in% rownames(Ex_ChP)]
cts_mjr <- cbind(cts_soupX_Epithelial, cts_soupX_Fibroblast, cts_soupX_Endothelial, cts_soupX_Mural, cts_soupX_Immune, cts_soupX_Parenchyma)
cts_mjr <- data.frame(scale(cts_mjr))
cts_mjr <- cts_mjr[rownames(cts_mjr) %in% my_genes,] #!!!!no parenchyma in my genes
cts_mjr <- cts_mjr[order(match(rownames(cts_mjr), my_genes)),]
cts_mjr <- t(scale(t(cts_mjr)))
htg <- Heatmap(cts_mjr, heatmap_legend_param = list(title = "Expression"), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, show_row_names = TRUE, col = colorRamp2(c(-2, 0, 2, 4), c("#000004FF","#721F81FF", "#F1605DFF", "#FCFDBFFF")))
draw(htg)
graph2png(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_batch_major.png", sep = ""), width=7, height=7)
graph2pdf(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_batch_major.pdf", sep = ""), width=7, height=7)



###Bubble VLN Feature QC plots
####Main

highlights_genes <- c("TTR", "HTR2C", "CFAP45", "DNAH7", "CLIC6", "COL1A1", "COL1A2", "COL4A2", "COL4A1", "DCN", "LUM", "FN1", "AIF1", "HLA-DMB", "STAB1", "CD163", "MSR1", "LYVE1", "F13A1", "MS4A7", "ABCA1", "AOAH", "IL10RA")
Idents(Ex_ChP) <- Ex_ChP@meta.data$Majorcelltype
plts <- DotPlot(Ex_ChP, features = highlights_genes, split.by="Majorcelltype", cols=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "brown")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Ex_ChP_bubble.pdf", sep=""), width=5+length(highlights_genes)*.15, height=1.5+length(unique(Ex_ChP@meta.data[["Majorcelltype"]]))*.4)

Idents(Ex_ChP) <- Ex_ChP@meta.data$Majorcelltype
plts <- DotPlot(Ex_ChP, features = highlights_genes) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Ex_ChP_bubble_scale.pdf", sep=""), width=5+length(highlights_genes)*.15, height=1.5+length(unique(Ex_ChP@meta.data[["Majorcelltype"]]))*.4)

plts <- VlnPlot(Ex_ChP, features = highlights_genes, split.by='Majorcelltype', stack=T, raster=T, cols=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "brown"), same.y.lims=T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Ex_ChP_vln.pdf", sep=""), width=3+length(highlights_genes)*.4, height=1.5+length(unique(Ex_ChP@meta.data[["Majorcelltype"]]))*1)

Featurefun(Ex_ChP, highlights_genes, "highlights_genes",  assay="soupXcounts", max.cutoff=c(3,4))

CNDRdimfun(Ex_ChP, reduction="harmony", res="Majorcelltype", split.by = c("X10X_batch", "projid", "Lane"), cols=c("#C48F45", "#3249A6", "#926392", "#261132", "#26532B", "brown"))
Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
CNDRdimfun(Ex_ChP, reduction="harmony", res="Subcelltype", split.by = c("X10X_batch", "projid", "Lane"), cols=c("#CD8500", "#FFA500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#c875c4", "#976697", "#c5aac5", "#261132", "#6d207b", "#ac32c3", "#70C25B", "#0E9554", "#266417", "#3ded97","lightsalmon", "salmon2", "sienna"))

for(var in c("X10X_batch", "projid", "Lane")){
  path <- paste(curr_dir, "/CellType/Ex_ChP/harmonySubcelltype/", sep="")
  Idents(Ex_ChP) <- Ex_ChP@meta.data[[var]]
  plts <- VlnPlot(Ex_ChP, "nCount_RNA", pt.size=F, y.max=10000)
  plot(plts)
  ggsave2(paste(path, "Ex_ChP_cluster_Vln_nCount_RNA_", var, ".pdf", sep=""), plot = plts, width=length(unique(Ex_ChP@meta.data[[var]]))*.2, height=4)
  
  plts <- VlnPlot(Ex_ChP, "nFeature_RNA", pt.size=F)
  plot(plts)
  ggsave2(paste(path, "Ex_ChP_cluster_Vln_nFeature_RNA_", var, ".pdf", sep=""), plot = plts, width=length(unique(Ex_ChP@meta.data[[var]]))*.2, height=4)
  
  plts <- VlnPlot(Ex_ChP, "percent.mt", pt.size=F)
  plot(plts)
  ggsave2(paste(path, "Ex_ChP_cluster_Vln_percent_mt_", var, ".pdf", sep=""), plot = plts, width=length(unique(Ex_ChP@meta.data[[var]]))*.2, height=4)
}


####Epi sub

Epi <- readRDS(paste0(curr_dir, "/CellType/Finalcelltype/Epi.RDS"))

Epi_genes <- c("TTR", "HTR2C", "CDH1", "HSP90AA1", "SLC26A3", "TCF12", "BNC2", "AUTS2", "ADAM10", "FOXJ1", "DNAH6", "CFAP100", "CFAP69", "SPAG17", "TUBB4B", "CLDN1", "CLDN5", "OCLN", "TJP1", "GPX3", "MT-CO1", "MT-CO2", "AQP1", "ATP1A1", "SLC25A6", "SLC25A23", "SLC4A2", "FOLR1", "IGFBP2")

plts <- DotPlot(Epi, features = Epi_genes, split.by="Subcelltype", cols=c("#CD8500", "#FFA500", "#FF4500", "#CD3700")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Epi_bubble.pdf", sep=""), width=3+length(Epi_genes)*.15, height=1.5+length(unique(Epi@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Epi, features = Epi_genes, cols=c("grey", "#C48F45")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Epi_bubble_scale.pdf", sep=""), width=5+length(Epi_genes)*.15, height=1.5+length(unique(Epi@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Epi, features = Epi_genes) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Epi_bubble_default.pdf", sep=""), width=5+length(Epi_genes)*.15, height=1.5+length(unique(Epi@meta.data[["Subcelltype"]]))*.4)

plts <- VlnPlot(Epi, features = Epi_genes, split.by='Subcelltype', stack=T, raster=T, cols=c("#CD8500", "#FFA500", "#FF4500", "#CD3700"), same.y.lims=F) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Epi_vln.pdf", sep=""), width=3+length(Epi_genes)*.4, height=1.5+length(unique(Epi@meta.data[["Subcelltype"]]))*1)
remove(plts)

plts <- VlnPlot(Epi, features = Epi_genes, split.by='Subcelltype', stack=T, raster=T, cols=c("#CD8500", "#FFA500", "#FF4500", "#CD3700"), same.y.lims=T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))  + geom_boxplot(outliers=F)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Epi_vln_box.pdf", sep=""), width=3+length(Epi_genes)*.4, height=1.5+length(unique(Epi@meta.data[["Subcelltype"]]))*1)
remove(plts)

Featurefun(Epi, Epi_genes, "highlights_genes",  assay="soupXcounts", max.cutoff=c(3,4))

Idents(Epi) <- Epi@meta.data$Subcelltype
CNDRdimfun(Epi, reduction="harmony", res="Subcelltype", split.by = c("X10X_batch", "projid", "Lane"), cols=c("#CD8500", "#FFA500", "#FF4500", "#CD3700"))


####Fibro sub

Fibro <- readRDS(paste0(curr_dir, "/CellType/Finalcelltype/Fibro.RDS"))

Fibro_genes <- c("LEPR", "CDH11", "STXBP6", "KCNMA1", "PLXNA4", "ALPL", "DPEP1", "DLK1", "NRG1", "ROBO1", "TRPM3", "MYRIP", "GPM6A", "CLDN11", "TJP1", "ECM2", "SLC6A20", "SLC13A3", "ABCA9", "LAMA2", "LAMA4", "COL5A2", "COL4A2", "COL1A2", "COL4A4", "COL3A1", "COL6A3", "DCN", "FN1", "PDGFRA", "SLC4A7", "SLC2A3", "LUM")
plts <- DotPlot(Fibro, features = Fibro_genes, split.by='Subcelltype', cols=c("#050382", "#27D7EF", "#2a91ea")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Fibro_bubble.pdf", sep=""), width=3+length(Fibro_genes)*.15, height=1.5+length(unique(Fibro@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Fibro, features = Fibro_genes, cols=c("grey", "#3249A6")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Fibro_bubble_scale.pdf", sep=""), width=5+length(Fibro_genes)*.15, height=1.5+length(unique(Fibro@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Fibro, features = Fibro_genes) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Fibro_bubble_default.pdf", sep=""), width=5+length(Fibro_genes)*.15, height=1.5+length(unique(Fibro@meta.data[["Subcelltype"]]))*.4)

plts <- VlnPlot(Fibro, features = Fibro_genes, split.by='Subcelltype', stack=T, raster=T, cols=c("#050382", "#27D7EF", "#2a91ea"), same.y.lims = T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Fibro_vln.pdf", sep=""), width=3+length(Fibro_genes)*.4, height=1.5+length(unique(Fibro@meta.data[["Subcelltype"]]))*1)
remove(plts)

Featurefun(Fibro, Fibro_genes, "highlights_genes",  assay="soupXcounts", max.cutoff=c(3,4))

Idents(Fibro) <- Fibro@meta.data$Subcelltype
CNDRdimfun(Fibro, reduction="harmony", res="Subcelltype", split.by = c("X10X_batch", "projid", "Lane"), cols=c("#050382", "#27D7EF", "#2a91ea"))


####Imm sub

Imm <- readRDS(paste0(curr_dir, "/CellType/Finalcelltype/Imm.RDS"))

Imm_genes <- c("F13A1", "LYVE1", "MRC1", "STAB1", "CD163", "CD68", "MERTK", "HCK", "IL4R", "CD74", "HLA-DRA", "S100A9", "CSF1R", "JAK2", "CD86", "VCAM1", "CD83", "LAMP1", "DNAJA1", "HSP90AA1", "HSPE1", "HSPD1", "HSPA8", "HSPB1", "MKI67", "CENPF", "TOP2A", "ITGAM", "S100A8", "CCR2", "PCNA", "EZH2", "GTSE1", "ANLN", "KIF11", "THEMIS", "CD3E", "PDE7A", "CAMK4", "HLA-B", "CD8A") %>% unique()

Idents(Imm) <- Imm@meta.data$Subcelltype
plts <- DotPlot(Imm, features = Imm_genes, split.by="Subcelltype", cols=c("#70C25B", "#0E9554", "#266417", "#3ded97")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Imm_bubble.pdf", sep=""), width=5+length(Imm_genes)*.15, height=1.5+length(unique(Imm@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Imm, features = Imm_genes, cols=c("grey", "#26532B")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Imm_bubble_scale.pdf", sep=""), width=5+length(Imm_genes)*.15, height=1.5+length(unique(Imm@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Imm, features = Imm_genes) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Imm_bubble_default.pdf", sep=""), width=5+length(Imm_genes)*.15, height=1.5+length(unique(Imm@meta.data[["Subcelltype"]]))*.4)

plts <- VlnPlot(Imm, features = Imm_genes, split.by='Subcelltype', stack=T, raster=T, cols=c("#70C25B", "#0E9554", "#266417", "#3ded97"), same.y.lims = T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Imm_vln.pdf", sep=""), width=3+length(Imm_genes)*.4, height=1.5+length(unique(Imm@meta.data[["Subcelltype"]]))*1)
remove(plts)

Featurefun(Imm, Imm_genes, "highlights_genes",  assay="soupXcounts", max.cutoff=c(3,4))

Idents(Imm) <- Imm@meta.data$Subcelltype
CNDRdimfun(Imm, reduction="harmony", res="Subcelltype", split.by = c("X10X_batch", "projid", "Lane"), cols=c("#70C25B", "#0E9554", "#266417", "#3ded97"))



####Endo sub

Endo <- readRDS(paste0(curr_dir, "/CellType/Finalcelltype/Endo.RDS"))

Endo_genes <- c("PECAM1", "VWF", "INSR", "CLDN5", "CDH13", "ARL15", "MECOM", "VEGFC", "COL8A1", "FBN1", "THSD7A", 
                "TSHZ2", "CEMIP2", "IL1R1", "IL4R", "FLI1", "SORBS2", 
                "PLVAP", "FLT1", "ADGRF5", "KCNQ3", "SLC2A3", "ITGA1") %>% unique()

Idents(Endo) <- Endo@meta.data$Subcelltype
plts <- DotPlot(Endo, features = Endo_genes, split.by="Subcelltype", cols=c("#c875c4", "#c5aac5", "#976697")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Endo_bubble.pdf", sep=""), width=5+length(Endo_genes)*.15, height=1.5+length(unique(Endo@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Endo, features = Endo_genes, cols=c("grey", "#926392")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Endo_bubble_scale.pdf", sep=""), width=5+length(Endo_genes)*.15, height=1.5+length(unique(Endo@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Endo, features = Endo_genes) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Endo_bubble_default.pdf", sep=""), width=5+length(Endo_genes)*.15, height=1.5+length(unique(Endo@meta.data[["Subcelltype"]]))*.4)

plts <- VlnPlot(Endo, features = Endo_genes, split.by='Subcelltype', stack=T, raster=T, cols=c("#c875c4", "#c5aac5", "#976697"), same.y.lims = T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Endo_vln.pdf", sep=""), width=3+length(Endo_genes)*.4, height=1.5+length(unique(Endo@meta.data[["Subcelltype"]]))*1)
remove(plts)

Featurefun(Endo, Endo_genes, "highlights_genes",  assay="soupXcounts", max.cutoff=c(3,4))

Idents(Endo) <- Endo@meta.data$Subcelltype
CNDRdimfun(Endo, reduction="harmony", res="Subcelltype", split.by = c("X10X_batch", "projid", "Lane"), cols=c("#c875c4", "#c5aac5", "#976697"))


####Mural sub

Mural <- readRDS(paste0(curr_dir, "/CellType/Finalcelltype/Mural.RDS"))

Mural_genes <- c("ENPP2","PLXDC2","PALMD", "MID1", "CEMIP", "GLIS3", "PDZD2", 
                 "COL4A4", "CRISPLD2", "ADAMTS9", "TLE1", "PIPNC1", "ADCY3", "HIPK2", "COL18A1", "SCL22A3",
                 "SOX6", "MYO1D", "MYH11", "ACTA2", "PRUNE2", "RYR2", "TPM1", "LMOD1") %>% unique()

Idents(Mural) <- Mural@meta.data$Subcelltype
plts <- DotPlot(Mural, features = Mural_genes, split.by="Subcelltype", cols=c("#8f2aa2", "#6d207b", "#261132")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Mural_bubble.pdf", sep=""), width=5+length(Mural_genes)*.15, height=1.5+length(unique(Mural@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Mural, features = Mural_genes, cols=c("grey", "#261132")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Mural_bubble_scale.pdf", sep=""), width=5+length(Mural_genes)*.15, height=1.5+length(unique(Mural@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Mural, features = Mural_genes) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Mural_bubble_default.pdf", sep=""), width=5+length(Mural_genes)*.15, height=1.5+length(unique(Mural@meta.data[["Subcelltype"]]))*.4)

plts <- VlnPlot(Mural, features = Mural_genes, split.by='Subcelltype', stack=T, raster=T, cols=c("#8f2aa2", "#6d207b", "#261132"), same.y.lims = T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Mural_vln.pdf", sep=""), width=3+length(Mural_genes)*.4, height=1.5+length(unique(Mural@meta.data[["Subcelltype"]]))*1)
remove(plts)

Featurefun(Mural, Mural_genes, "highlights_genes",  assay="soupXcounts", max.cutoff=c(3,4))

Idents(Mural) <- Mural@meta.data$Subcelltype
CNDRdimfun(Mural, reduction="harmony", res="Subcelltype", split.by = c("X10X_batch", "projid", "Lane"), cols=c("#8f2aa2", "#6d207b", "#261132"))



####All sub

Ex_ChP <- readRDS(paste0(curr_dir, "/CellType/Finalcelltype/Ex_ChP_noparen.RDS"))

Epi_genes <- c("TTR", "HTR2C", "CLIC6", 
               "CDH1", "HSP90AA1", "SLC26A3", "TCF12", "BNC2", "AUTS2", "ADAM10",
               "FOXJ1", "DNAH6", "CFAP100", "CFAP69", "SPAG17", "TUBB4B",
               "CLDN1", "CLDN5", "TJP1", "MT-CO1", "MT-CO2", "ATP1A1", "OCLN")
Fibro_genes <- c("LEPR", "CDH11", "STXBP6",
                 "KCNMA1", "SLC4A4", "FOXO1", "KLHL13", "PDGFC", "SLC24A3", "JAM3", "JAM2", 
                 "LAMA4", "CLDN11", "NRG1", "PDE1A", "ADAMTS12", "ADAMTSL1", "ROBO2", "ABCA9", "BICC1", "THBS2",
                 "COL6A3", "DCN", "FN1", "COL1A2", "COL27A1", "COL4A2", "LUM")
Endo_genes <- c("PECAM1", "VWF", "CLDN5", "CDH13", "ARL15", "MECOM", "VEGFC", "COL8A1", "FBN1", "THSD7A", 
                "TSHZ2", "CEMIP2", "IL1R1", "IL4R", "FLI1",
                "PLVAP", "FLT1", "ADGRF5", "SLC2A3", "INSR", "ITGA1") %>% unique()
Mural_genes <- c("SOX6", "MYO1D", "MYH11", "ACTA2", "PRUNE2", "RYR2", "TPM1", "LMOD1",
                 "COL4A4", "CRISPLD2", "ADAMTS9", "TLE1", "PIPNC1", "ADCY3", "HIPK2", "COL18A1", "SCL22A3",
                 "ENPP2","PLXDC2","PALMD", "MID1", "CEMIP", "GLIS3", "PDZD2") %>% unique()
Imm_genes <- c("F13A1", "MRC1", "STAB1", "CD163",
               "MERTK", "HCK", "IL4R",
               "DNAJA1", "HSP90AA1", "HSPE1", "HSPD1", "HSPA8", "HSPB1",
               "MKI67", "CENPF", "TOP2A", "PCNA", "EZH2", "GTSE1", "ANLN", "KIF11",
               "THEMIS", "CD3E", "PDE7A", "CAMK4", "HLA-B", "CD8A") %>% unique()
All_genes <- unique(c(Epi_genes, Fibro_genes, Endo_genes, Mural_genes, Imm_genes))

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
plts <- DotPlot(Ex_ChP, features = All_genes, split.by="Subcelltype", cols=c("#CD8500", "#FFA500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#c875c4", "#c5aac5", "#976697", "#8f2aa2", "#6d207b", "#261132", "#70C25B", "#0E9554", "#266417", "#3ded97")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Ex_ChP_sub_bubble.pdf", sep=""), width=5+length(All_genes)*.15, height=1.5+length(unique(Ex_ChP@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Ex_ChP, features = All_genes, split.by="Subcelltype", cols=c("#C48F45", "#C48F45", "#C48F45", "#C48F45", "#3249A6", "#3249A6", "#3249A6", "#926392", "#926392", "#926392", "#261132", "#261132", "#261132", "#26532B", "#26532B", "#26532B", "#26532B")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Ex_ChP_sub_bubble_limit cols.pdf", sep=""), width=5+length(All_genes)*.15, height=1.5+length(unique(Ex_ChP@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Ex_ChP, features = All_genes) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Ex_ChP_sub_bubble_default.pdf", sep=""), width=5+length(All_genes)*.15, height=1.5+length(unique(Ex_ChP@meta.data[["Subcelltype"]]))*.4)

plts <- VlnPlot(Ex_ChP, features = All_genes, split.by='Subcelltype', stack=T, raster=T, cols=c("#C48F45", "#C48F45", "#C48F45", "#C48F45", "#3249A6", "#3249A6", "#3249A6", "#926392", "#926392", "#926392", "#261132", "#261132", "#261132", "#26532B", "#26532B", "#26532B", "#26532B"), same.y.lims = T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Ex_ChP_sub_vln.pdf", sep=""), width=3+length(All_genes)*.4, height=1.5+length(unique(Ex_ChP@meta.data[["Subcelltype"]]))*1)
remove(plts)

CNDRdimfun(Ex_ChP, reduction="harmony", res="Subcelltype", split.by = c("X10X_batch", "projid", "Lane"), cols=c("#CD8500", "#FFA500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#c875c4", "#c5aac5", "#976697", "#8f2aa2", "#6d207b", "#261132", "#70C25B", "#0E9554", "#266417", "#3ded97"))

