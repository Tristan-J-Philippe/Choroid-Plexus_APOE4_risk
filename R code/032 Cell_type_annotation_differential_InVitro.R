###
### 032 Cell_type_annotation_differential_RNA.R 
###
# Purpose: Cell type annotation differential analysis and gsea of single nuclei RNA dataset.
# Dependencies:
source("/030 limma_pattern_fgsea_util.R", verbose = FALSE)
curr_dir <- "/Tsai" #output results to specific folder similarly "/Jackson", "Jco_cult"
set.seed(42)


##limma on celltypes
###Housekeeping

library(lattice)
library(limma)
library(edgeR)
library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(ggrepel)
library(splines)
library(stringr)
library(tidyverse)
library(ggVennDiagram)
library(dplyr)

set.seed(42)

#Set and load dependent variables for limma/DESeq functions
Gene.search <- read.csv("/home/tristan-j-philippe/DEG to GO.csv")
FDRcutoff <- 0.01
LFcutoff <- 0.25
pct_cutoff <- 0.1

limma_resfun <- function(coef){
  #Small function to extract coefficients (and save as csv) and display a p value distribution
  #Fill in coef = whatever the coefficient is
  #Expects to receive Mainf from global environment
  cutoff = 1
  df.1a <- topTable(Mainf, coef = coef, p.value=cutoff, n=length(Mainf[["Amean"]]))
  df.1a <- df.1a[order(row.names(df.1a)), ]
  vdata <- vdata[order(row.names(vdata)), ]
  df.1a <- cbind(df.1a, vdata)
  
  path <- paste(curr_dir, "/CellType/limma/", sep="") #specify dir
  if(dir.exists(path)==TRUE){print("Directory exists")}else{dir.create(file.path(path), recursive = TRUE)} #check dir exists
  
  hist(df.1a[,"P.Value"], main = paste(coef), ylab = "Number of genes", xlab = "P value")
  
  limma::plotMA(Mainf, coef = coef, main = coef)
  
  df.1b <-df.1a
  df.1b$diffexpressed <- "NS"
  df.1b$diffexpressed[df.1b$logFC > LFcutoff & df.1b$adj.P.Val < FDRcutoff] <- "UP"
  df.1b$diffexpressed[df.1b$logFC < -LFcutoff & df.1b$adj.P.Val < FDRcutoff] <- "DOWN"
  mycolors <- c("blue", "red", "black")
  names(mycolors) <- c("DOWN", "UP", "NS")
  
  plts <- ggplot(data=df.1b, aes(x=logFC, y=-log10(adj.P.Val), col = diffexpressed)) +
    geom_point() + theme_minimal() +
    scale_color_manual(values=mycolors) +
    geom_vline(xintercept=c(-LFcutoff, LFcutoff), col="red") +
    geom_hline(yintercept=-log10(FDRcutoff), col="red")
  plot(plts)
  graph2png(plts, paste(path, coef, "_Volcano.png", sep=""), width=4, height=2.5)
  
  df.2a <- df.1a %>% rownames_to_column("gene")
  df.2a <- left_join(df.2a, Gene.search, by = "gene")
  #cts <- cts %>% rownames_to_column("gene")
  #df.2a <- merge(df.2a, cts, by = "gene")
  df.2a <- df.2a[order(df.2a$adj.P.Val),, drop = FALSE]
  df.2a <- subset(df.2a, logFC > 0)
  write.csv(df.2a, file = paste(path, "Spreadsheet_", coef, "_uncorrected.csv", sep = ""))
  assign(c(paste("limma", coef, sep = "_")), df.2a, envir = .GlobalEnv) #Save in Global environment
  
  df.2c <- subset(df.2a, adj.P.Val < FDRcutoff)
  df.2c <- subset(df.2c, logFC > LFcutoff)
  if(nrow(df.2c)>1){
    write.csv(df.2c, file = paste(path, "Spreadsheet_",  coef, "_corrected.csv", sep = ""))
    assign(c(paste("limma", coef, "sig", sep = "_")), df.2c, envir = .GlobalEnv) #Save in Global environment
  }
}


###Subcelltype

###Heatmap by pseudobulk
cts_fun <- function(celltype){
  sub <- subset(Org_ChP, Subcelltype==celltype)
  cts <- AggregateExpression(sub, group.by = "X10X_batch", assay = "soupXcounts", slot = "counts", return.seurat = FALSE)
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

cts_fun("BAM_1")
cts_fun("BAM_2")
cts_fun("BAM_3")

cts_fun("Neuron")

cts_fun("Prog_Vasc")
cts_fun("Prog_Fib")
cts_fun("Prog_Epi")
cts_fun("Prog_Neuro_Epi")
cts_fun("Prog_Neuro")

cts_fun("Dividing_Fib")
cts_fun("Dividing_Epi")
cts_fun("Dividing_Neuro_Epi")

dim(cts_soupX_Epi_1)
dim(cts_soupX_Epi_2)
dim(cts_soupX_Epi_3)
dim(cts_soupX_Epi_4)

dim(cts_soupX_Fib_1)
dim(cts_soupX_Fib_2)
dim(cts_soupX_Fib_3)

dim(cts_soupX_BAM_1)
dim(cts_soupX_BAM_2)
dim(cts_soupX_BAM_3)

dim(cts_soupX_Prog_Vasc)
dim(cts_soupX_Prog_Fib)
dim(cts_soupX_Prog_Epi)
dim(cts_soupX_Prog_Neuro_Epi)
dim(cts_soupX_Prog_Neuro)

dim(cts_soupX_Dividing_Fib)
dim(cts_soupX_Dividing_Epi)
dim(cts_soupX_Dividing_Neuro_Epi)


###Within Epi

cts_sub_paren <- cbind(cts_soupX_Epi_1, cts_soupX_Epi_2, cts_soupX_Epi_3, cts_soupX_Epi_4)
colnames(cts_sub_paren) <- c(rep("Epi_1", 9), rep("Epi_2", 9), rep("Epi_3", 9), rep("Epi_4", 9))

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
colnames(cts_sub_paren) <- c(rep("Fib_1", 9), rep("Fib_2", 9), rep("Fib_3", 9))

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


###Within Immune

cts_sub_paren <- cbind(cts_soupX_BAM_1, cts_soupX_BAM_2, cts_soupX_BAM_3)
colnames(cts_sub_paren) <- c(rep("BAM_1", 9), rep("BAM_2", 7), rep("BAM_3", 8))

celltype <- factor(paste(colnames(cts_sub_paren)))

options(na.action='na.pass')
Maindm <- model.matrix(~0 + celltype)

rownames(Maindm) <- colnames(cts_sub_paren)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
Maincm <- makeContrasts(levels=Maindm, 
                        BAM_1=celltypeBAM_1-(celltypeBAM_2+celltypeBAM_3)/2, 
                        BAM_2=celltypeBAM_2-(celltypeBAM_1+celltypeBAM_3)/2, 
                        BAM_3=celltypeBAM_3-(celltypeBAM_1+celltypeBAM_2)/2)

cts_sub_paren <- cts_sub_paren[rowSums(cts_sub_paren)>2*ncol(cts_sub_paren),]
vdat <- voom(cts_sub_paren, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

for(coef in colnames(Mainf)){
  limma_resfun(coef)
}


###Within Prog

cts_sub_paren <- cbind(cts_soupX_Prog_Vasc, cts_soupX_Prog_Fib, cts_soupX_Prog_Epi, cts_soupX_Prog_Neuro_Epi, cts_soupX_Prog_Neuro)
colnames(cts_sub_paren) <- c(rep("Prog_Vasc", 9), rep("Prog_Fib", 9), rep("Prog_Epi", 9), rep("Prog_Neuro_Epi", 9), rep("Prog_Neuro", 9))

celltype <- factor(paste(colnames(cts_sub_paren)))

options(na.action='na.pass')
Maindm <- model.matrix(~0 + celltype)

rownames(Maindm) <- colnames(cts_sub_paren)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
Maincm <- makeContrasts(levels=Maindm, 
                        Prog_Vasc=celltypeProg_Vasc-(celltypeProg_Fib+celltypeProg_Epi+celltypeProg_Neuro_Epi+celltypeProg_Neuro)/4, 
                        Prog_Fib=celltypeProg_Fib-(celltypeProg_Vasc+celltypeProg_Epi+celltypeProg_Neuro_Epi+celltypeProg_Neuro)/4, 
                        Prog_Epi=celltypeProg_Epi-(celltypeProg_Vasc+celltypeProg_Fib+celltypeProg_Neuro_Epi+celltypeProg_Neuro)/4, 
                        Prog_Neuro_Epi=celltypeProg_Neuro_Epi-(celltypeProg_Vasc+celltypeProg_Fib+celltypeProg_Epi+celltypeProg_Neuro)/4, 
                        Prog_Neuro=celltypeProg_Neuro-(celltypeProg_Vasc+celltypeProg_Fib+celltypeProg_Epi+celltypeProg_Neuro_Epi)/4 
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


###Within Dividing

cts_sub_paren <- cbind(cts_soupX_Dividing_Fib, cts_soupX_Dividing_Epi, cts_soupX_Dividing_Neuro_Epi)
colnames(cts_sub_paren) <- c(rep("Dividing_Fib", 9), rep("Dividing_Epi", 9),  rep("Dividing_Neuro_Epi", 9))

celltype <- factor(paste(colnames(cts_sub_paren)))

options(na.action='na.pass')
Maindm <- model.matrix(~0 + celltype)

rownames(Maindm) <- colnames(cts_sub_paren)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
cts_sub_paren <- cts_sub_paren[,colnames(cts_sub_paren)%in%rownames(Maindm)]
Maincm <- makeContrasts(levels=Maindm, 
                        Dividing_Fib=celltypeDividing_Fib-(celltypeDividing_Epi+celltypeDividing_Neuro_Epi)/2,
                        Dividing_Epi=celltypeDividing_Epi-(celltypeDividing_Fib+celltypeDividing_Neuro_Epi)/2,
                        Dividing_Neuro_Epi=celltypeDividing_Neuro_Epi-(celltypeDividing_Fib+celltypeDividing_Epi)/2)

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
  sub <- subset(Org_ChP, Majorcelltype==celltype)
  cts <- AggregateExpression(sub, group.by = "X10X_batch", assay = "soupXcounts", slot = "counts", return.seurat = FALSE)
  cts <- cts[["soupXcounts"]]
  assign(c(paste("cts_soupX", celltype, sep = "_")), cts, envir = .GlobalEnv)
  pct <- pctincell(sub)
  assign(c(paste("pct", celltype, sep = "_")), pct, envir = .GlobalEnv)
}

cts_fun("Epithelial")
gc()
cts_fun("Fibroblast")
gc()
cts_fun("Immune")
gc()
cts_fun("Neuron")
gc()
cts_fun("Progenitor")
gc()
cts_fun("Dividing")

dim(cts_soupX_Epithelial)
dim(cts_soupX_Fibroblast)
dim(cts_soupX_Immune)
dim(cts_soupX_Progenitor)
dim(cts_soupX_Dividing)
dim(cts_soupX_Neuron)

cts_major_paren <- cbind(cts_soupX_Epithelial, cts_soupX_Fibroblast, cts_soupX_Immune, cts_soupX_Progenitor, cts_soupX_Dividing, cts_soupX_Neuron)
colnames(cts_major_paren) <- c(rep("Epithelial", 9), rep("Fibroblast", 9), rep("Immune", 9), rep("Progenitor", 9), rep("Dividing", 9), rep("Neuron", 9))

celltype <- factor(paste(colnames(cts_major_paren)))

options(na.action='na.pass')
Maindm <- model.matrix(~0 + celltype)

rownames(Maindm) <- colnames(cts_major_paren)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
cts_major_paren <- cts_major_paren[,colnames(cts_major_paren)%in%rownames(Maindm)]
cts_major_paren <- cts_major_paren[,colnames(cts_major_paren)%in%rownames(Maindm)]
Maincm <- makeContrasts(levels=Maindm, Epithelial=celltypeEpithelial-(celltypeFibroblast+celltypeImmune+celltypeProgenitor+celltypeDividing+celltypeNeuron)/5, Fibroblast=celltypeFibroblast-(celltypeEpithelial+celltypeImmune+celltypeProgenitor+celltypeDividing+celltypeNeuron)/5, Immune=celltypeImmune-(celltypeEpithelial+celltypeFibroblast+celltypeProgenitor+celltypeDividing+celltypeNeuron)/5, Progenitor=celltypeProgenitor-(celltypeEpithelial+celltypeFibroblast+celltypeImmune+celltypeDividing+celltypeNeuron)/5, Dividing=celltypeDividing-(celltypeEpithelial+celltypeFibroblast+celltypeImmune+celltypeProgenitor+celltypeNeuron)/5, Neuron=celltypeNeuron-(celltypeEpithelial+celltypeFibroblast+celltypeImmune+celltypeProgenitor+celltypeDividing)/5)

cts_major_paren <- cts_major_paren[rowSums(cts_major_paren)>2*ncol(cts_major_paren),]
vdat <- voom(cts_major_paren, design = Maindm, plot=T)
vdata <- data.frame(vdat)#limma
Mainf <- lmFit(vdat, weights=vdat$weights, Maindm)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf)

for(coef in colnames(Mainf)){
  limma_resfun(coef)
}



####Go run

#cogn_global
Path_fun(limma_Epi_1a)
Path_fun(limma_Epi_1b)
Path_fun(limma_Epi_2a)
Path_fun(limma_Epi_2b)

Path_fun(limma_Fib_1)
Path_fun(limma_Fib_2)
Path_fun(limma_Fib_3)

Path_fun(limma_BAM_1)
Path_fun(limma_BAM_2)
Path_fun(limma_BAM_3)

Path_fun(limma_Prog_Vasc)
Path_fun(limma_Prog_Fib)
Path_fun(limma_Prog_Epi)
Path_fun(limma_Prog_Neuro_Epi)
Path_fun(limma_Prog_Neuro)

Path_fun(limma_Dividing_Fib)
Path_fun(limma_Dividing_Epi)
Path_fun(limma_Dividing_Neuro_Epi)

Path_fun(limma_Epithelial)
Path_fun(limma_Fibroblast)
Path_fun(limma_Immune)
Path_fun(limma_Progenitor)
Path_fun(limma_Dividing)
Path_fun(limma_Neuron)




##Figures Celltype
###UMAP

#Major types
Idents(Org_ChP) <- Org_ChP@meta.data[["Majorcelltype"]]
plts <- DimPlot(Org_ChP, raster=F,raster.dpi=c(600,600), cols=c("Epithelial"="#C48F45", "Fibroblast"="#3249A6",  "Immune"="#26532B", "Progenitor"= "#926392", "Dividing"="#261132", "Neuron"="#784e00")) + theme(text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Org_ChP_UMAP_major.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Org_ChP_UMAP_major.png", sep=""), width=10, height=5)
remove(plts)

#Subcelltypes
Idents(Org_ChP) <- Org_ChP@meta.data[["Subcelltype"]] 
plts <- DimPlot(Org_ChP, raster=F, raster.dpi=c(600,600), cols=c("Epi_1a"="#FFA500", "Epi_1b"="#CD8500", "Epi_2a"="#FF4500", "Epi_2b"="#CD3700", "Fib_1"="#050382", "Fib_2"="#27D7EF", "Fib_3"="#2a91ea", "BAM_1"="#70C25B", "BAM_2"="#0E9554", "BAM_3"="#266417", "Prog_Vasc"="#6d207b", "Prog_Fib"="blue", "Prog_Epi"="#ff7450", "Prog_Neuro_Epi"="#d16120", "Prog_Neuro"="#a55d2a", "Dividing_Fib"="lightblue", "Dividing_Epi"="#ffc274", "Dividing_Neuro_Epi"="#e8a281", "Dividing_Neuro"="#bba486", "Neuron"="#784e00"))  + theme(text=element_text(size=7))#sanity check and make publishable figs
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Org_ChP_UMAP.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Org_ChP_UMAP.png", sep=""), width=10, height=5)
remove(plts)

#Epi
Idents(Epi) <- Epi@meta.data[["Subcelltype"]] 
plts <- DimPlot(Epi, raster=F, raster.dpi=c(600,600), cols=c("Epi_1a"="#FFA500", "Epi_1b"="#CD8500", "Epi_2a"="#FF4500", "Epi_2b"="#CD3700"))  + theme(text=element_text(size=7))#sanity check and make publishable figs
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Epi_UMAP.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Epi_UMAP.png", sep=""), width=10, height=5)
remove(plts)

#Fib
Idents(Fib) <- Fib@meta.data[["Subcelltype"]] 
plts <- DimPlot(Fib, raster=F, raster.dpi=c(600,600), cols=c("Fib_1"="#050382", "Fib_2"="#27D7EF", "Fib_3"="#2a91ea"))  + theme(text=element_text(size=7))#sanity check and make publishable figs
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Fib_UMAP.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Fib_UMAP.png", sep=""), width=10, height=5)
remove(plts)

#Imm
Idents(Imm) <- Imm@meta.data[["Subcelltype"]] 
plts <- DimPlot(Imm, raster=F, raster.dpi=c(600,600), cols=c("BAM_1"="#70C25B", "BAM_2"="#0E9554", "BAM_3"="#266417"))  + theme(text=element_text(size=7))#sanity check and make publishable figs
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Imm_UMAP.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Imm_UMAP.png", sep=""), width=10, height=5)
remove(plts)

#Prog
Idents(Prog) <- Prog@meta.data[["Subcelltype"]] 
plts <- DimPlot(Prog, raster=F, raster.dpi=c(600,600), cols=c("Prog_Vasc"="#6d207b", "Prog_Fib"="blue", "Prog_Epi"="#ff7450", "Prog_Neuro_Epi"="#d16120", "Prog_Neuro"="#a55d2a"))  + theme(text=element_text(size=7))#sanity check and make publishable figs
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Prog_UMAP.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Prog_UMAP.png", sep=""), width=10, height=5)
remove(plts)

#Div
Idents(Dividing) <- Dividing@meta.data[["Subcelltype"]] 
plts <- DimPlot(Dividing, raster=F, raster.dpi=c(600,600), cols=c("Dividing_Fib"="lightblue", "Dividing_Epi"="#ffc274", "Dividing_Neuro_Epi"="#e8a281", "Dividing_Neuro"="#bba486"))  + theme(text=element_text(size=7))#sanity check and make publishable figs
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Dividing_UMAP.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Dividing_UMAP.png", sep=""), width=10, height=5)
remove(plts)

##no prog/dividing/neuro
Org_ChP <- subset(Org_ChP, Majorcelltype=="Progenitor" | Majorcelltype=="Dividing" | Majorcelltype=="Neuron", invert=T)
Org_ChP <- RunUMAP(Org_ChP, dims=1:20, assay="soupXcounts", reduction="harmony", min.dist=1)
#Major types
Idents(Org_ChP) <- Org_ChP@meta.data[["Majorcelltype"]]
plts <- DimPlot(Org_ChP, raster=F,raster.dpi=c(600,600), cols=c("Epithelial"="#C48F45", "Fibroblast"="#3249A6",  "Immune"="#26532B")) + theme(text=element_text(size=7))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Org_ChP_UMAP_major_noextra.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Org_ChP_UMAP_major_noextra.png", sep=""), width=10, height=5)
remove(plts)

#Subcelltypes
Idents(Org_ChP) <- Org_ChP@meta.data[["Subcelltype"]] 
plts <- DimPlot(Org_ChP, raster=F, raster.dpi=c(600,600), cols=c("Epi_1a"="#FFA500", "Epi_1b"="#CD8500", "Epi_2a"="#FF4500", "Epi_2b"="#CD3700", "Fib_1"="#050382", "Fib_2"="#27D7EF", "Fib_3"="#2a91ea", "BAM_1"="#70C25B", "BAM_2"="#0E9554", "BAM_3"="#266417"))  + theme(text=element_text(size=7))#sanity check and make publishable figs
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Org_ChP_UMAP_noextra.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Org_ChP_UMAP_noextra.png", sep=""), width=10, height=5)
remove(plts)

Org_ChP33 <- subset(Org_ChP, Genotype=="APOE33", invert=T)
#Subcelltypes
Idents(Org_ChP33) <- Org_ChP33@meta.data[["Subcelltype"]] 
plts <- DimPlot(Org_ChP33, raster=F, raster.dpi=c(600,600), cols=c("Epi_1a"="#FFA500", "Epi_1b"="#CD8500", "Epi_2a"="#FF4500", "Epi_2b"="#CD3700", "Fib_1"="#050382", "Fib_2"="#27D7EF", "Fib_3"="#2a91ea", "BAM_1"="#70C25B", "BAM_2"="#0E9554", "BAM_3"="#266417"))  + theme(text=element_text(size=7))#sanity check and make publishable figs
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Org_ChP_UMAP_noextra_APOE33.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Org_ChP_UMAP_noextra_APOE33.png", sep=""), width=10, height=5)
remove(plts)

Org_ChP44 <- subset(Org_ChP, Genotype=="APOE44", invert=T)
#Subcelltypes
Idents(Org_ChP44) <- Org_ChP44@meta.data[["Subcelltype"]] 
plts <- DimPlot(Org_ChP44, raster=F, raster.dpi=c(600,600), cols=c("Epi_1a"="#FFA500", "Epi_1b"="#CD8500", "Epi_2a"="#FF4500", "Epi_2b"="#CD3700", "Fib_1"="#050382", "Fib_2"="#27D7EF", "Fib_3"="#2a91ea", "BAM_1"="#70C25B", "BAM_2"="#0E9554", "BAM_3"="#266417"))  + theme(text=element_text(size=7))#sanity check and make publishable figs
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Org_ChP_UMAP_noextra_APOE44.pdf", sep=""), width=10, height=5)
graph2png(plts, paste(curr_dir, "/CellType/Publish/Org_ChP_UMAP_noextra_APOE44.png", sep=""), width=10, height=5)
remove(plts)


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
  sub <- subset(Org_ChP, Subcelltype==celltype)
  cts <- AggregateExpression(sub, group.by = "Differentiation", assay = "soupXcounts", slot = "counts", return.seurat = FALSE)
  cts <- cts[["soupXcounts"]]
  assign(c(paste("cts", celltype, sep = "_")), cts, envir = .GlobalEnv)
  gc()
}

cts_fun("Epi_1a")
gc()
cts_fun("Epi_1b")
gc()
cts_fun("Epi_2a")
gc()
cts_fun("Epi_2b")
gc()
cts_fun("Fib_1")
cts_fun("Fib_2")
cts_fun("Fib_3")
cts_fun("BAM_1")
cts_fun("BAM_2")
cts_fun("BAM_3")


cts_fun <- function(celltype){
  sub <- subset(Org_ChP, Majorcelltype==celltype)
  cts <- AggregateExpression(sub, group.by = "Differentiation", assay = "soupXcounts", slot = "counts", return.seurat = FALSE)
  cts <- cts[["soupXcounts"]]
  assign(c(paste("cts", celltype, sep = "_")), cts, envir = .GlobalEnv)
  gc()
}

cts_fun("Epithelial")
cts_fun("Fibroblast")
cts_fun("Immune")
cts_fun("Neuron")

#sub no paren
my_genes <- c(Epi1a_genes, Epi1b_genes, Epi_genes, Epi2_genes, Fib_genes, Fib1_genes, Fib2_genes, Fib3_genes, BAM_genes, BAM1_genes, BAM2_genes, BAM3_genes) %>% unique()
my_genes <- my_genes[my_genes %in% rownames(Org_ChP)]

cts_sub <- cbind(cts_Epi_1a, cts_Epi_1b, cts_Epi_2a, cts_Epi_2b, cts_Fib_1, cts_Fib_2, cts_Fib_3, cts_BAM_1, cts_BAM_2, cts_BAM_3)

cts_sub <- data.frame(scale(cts_sub))
cts_sub <- cts_sub[rownames(cts_sub) %in% my_genes,]
cts_sub <- cts_sub[order(match(rownames(cts_sub), my_genes)),]
cts_sub <- t(scale(t(cts_sub)))
htg <- Heatmap(cts_sub, heatmap_legend_param = list(title = "Expression"), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, show_row_names = TRUE, col = colorRamp2(c(-2, 0, 2, 4), c("#000004FF","#721F81FF", "#F1605DFF", "#FCFDBFFF")))
draw(htg)
graph2png(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_batch_noparen.png", sep = ""), width=7, height=7)
graph2pdf(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_batch_noparen.pdf", sep = ""), width=7, height=7)

#major no paren
cts_mjr <- cbind(cts_Epithelial, cts_Fibroblast, cts_Immune)
cts_mjr <- data.frame(scale(cts_mjr))
cts_mjr <- cts_mjr[rownames(cts_mjr) %in% my_genes,] #!!!!no parenchyma in my genes
cts_mjr <- cts_mjr[order(match(rownames(cts_mjr), my_genes)),]
cts_mjr <- t(scale(t(cts_mjr)))
htg <- Heatmap(cts_mjr, heatmap_legend_param = list(title = "Expression"), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, show_row_names = TRUE, col = colorRamp2(c(-2, 0, 2, 4), c("#000004FF","#721F81FF", "#F1605DFF", "#FCFDBFFF")))
draw(htg)
graph2png(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_batch_major_noparen.png", sep = ""), width=7, height=7)
graph2pdf(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_batch_major_noparen.pdf", sep = ""), width=7, height=7)


#sub
my_genes <- c(Epi1a_genes, Epi1b_genes, Epi_genes, Epi2_genes, Fib_genes, Fib1_genes, Fib2_genes, Fib3_genes, BAM_genes, BAM1_genes, BAM2_genes, BAM3_genes, paren_genes) %>% unique()
my_genes <- my_genes[my_genes %in% rownames(Org_ChP)]
cts_sub <- cbind(cts_Epi_1a, cts_Epi_1b, cts_Epi_2a, cts_Epi_2b, cts_Fib_1, cts_Fib_2, cts_Fib_3, cts_BAM_1, cts_BAM_2, cts_BAM_3, cts_Neuron)
cts_sub <- data.frame(scale(cts_sub))
cts_sub <- cts_sub[rownames(cts_sub) %in% my_genes,] #!!!!no parenchyma in my genes
cts_sub <- cts_sub[order(match(rownames(cts_sub), my_genes)),]
cts_sub <- t(scale(t(cts_sub)))
htg <- Heatmap(cts_sub, heatmap_legend_param = list(title = "Expression"), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, show_row_names = TRUE, col = colorRamp2(c(-2, 0, 2, 4), c("#000004FF","#721F81FF", "#F1605DFF", "#FCFDBFFF")))
draw(htg)
graph2png(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_batch.png", sep = ""), width=7, height=7)
graph2pdf(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_batch.pdf", sep = ""), width=7, height=7)

#major with paren
cts_mjr <- cbind(cts_Epithelial, cts_Fibroblast, cts_Immune, cts_Neuron)
cts_mjr <- data.frame(scale(cts_mjr))
cts_mjr <- cts_mjr[rownames(cts_mjr) %in% my_genes,] #!!!!no parenchyma in my genes
cts_mjr <- cts_mjr[order(match(rownames(cts_mjr), my_genes)),]
cts_mjr <- t(scale(t(cts_mjr)))
htg <- Heatmap(cts_mjr, heatmap_legend_param = list(title = "Expression"), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, show_row_names = TRUE, col = colorRamp2(c(-2, 0, 2, 4), c("#000004FF","#721F81FF", "#F1605DFF", "#FCFDBFFF")))
draw(htg)
graph2png(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_batch_major.png", sep = ""), width=7, height=7)
graph2pdf(htg, file=paste(curr_dir, "/CellType/Publish/Heatmap_batch_major.pdf", sep = ""), width=7, height=7)



###Bubble VLN Feature QC plots
####Major

highlights_genes <- c("TTR", "HTR2C", "CFAP45", "DNAH7", "CLIC6", "COL1A1", "COL1A2", "COL4A2", "COL4A1", "DCN", "LUM", "FN1", "AIF1", "HLA-DMB", "STAB1", "CD163", "MSR1", "LYVE1", "F13A1", "MS4A7", "ABCA1", "AOAH", "IL10RA")
Idents(Org_ChP) <- (Org_ChP$Majorcelltype)

plts <- DotPlot(Org_ChP, features = highlights_genes, split.by="Majorcelltype", cols=c("#C48F45", "#3249A6", "#26532B", "#926392", "#261132", "brown")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Org_ChP_bubble.pdf", sep=""), width=5+length(highlights_genes)*.15, height=1.5+length(unique(Org_ChP@meta.data[["Majorcelltype"]]))*.4)

plts <- VlnPlot(Org_ChP, features = highlights_genes, split.by='Majorcelltype', stack=T, raster=T, cols=c("#C48F45", "#3249A6", "#26532B", "#926392", "#261132", "brown"), same.y.lims=T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Org_ChP_vln.pdf", sep=""), width=3+length(highlights_genes)*.4, height=1.5+length(unique(Org_ChP@meta.data[["Majorcelltype"]]))*1)

plts <- VlnPlot(Org_ChP, features = highlights_genes, split.by='Majorcelltype', stack=T, raster=T, cols=c("#C48F45", "#3249A6", "#26532B", "#926392", "#261132", "brown"), same.y.lims=T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5)) + geom_boxplot(outliers=F)
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Org_ChP_vln_box.pdf", sep=""), width=3+length(highlights_genes)*.4, height=1.5+length(unique(Org_ChP@meta.data[["Majorcelltype"]]))*1)


####Epi sub

#Epi <- readRDS(paste0(curr_dir, "/CellType/Finalcelltype/Epi.RDS"))
Idents(Epi) <- Epi@meta.data$Subcelltype
levels(Epi) <- c("Epi_1a", "Epi_1b", "Epi_2a", "Epi_2b")
Idents(Epi) -> Epi@meta.data$Subcelltype

Epi_genes <- c("TTR", "HTR2C", "CLIC6", 
               "FRMD4B", "BMPR2", "RNF38", "ERC1",
               "FOXJ1", "DNAH6", "CFAP100", "CFAP69", "SPAG17", "TUBB4B",
               "CLDN1", "CLDN5", "TJP1", "MT-CO1", "MT-CO2", "ATP1A1", "OCLN")

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


####Fibro sub

#Fib <- readRDS(paste0(curr_dir, "/CellType/Finalcelltype/Fib.RDS"))

Idents(Fib) <- Fib$Subcelltype
Fib_genes <- c("LEPR", "CDH11",
               "KCNMA1", "SLC4A4", "SLC24A3", "FOXO1", "KLHL13", "PDGFC", "JAM3", "CAV3",
               "LAMA4", "LAMB1", "CLDN11", "ADAMTS12", "ADAMTSL1", "ROBO2", "THBS2", "COL6A6",
               "COL4A2", "COL1A2", "COL5A2", "COL3A1", "SLC16A3", "SLC2A3", "LUM")

plts <- DotPlot(Fib, features = Fib_genes, split.by='Subcelltype', cols=c("#050382", "#27D7EF", "#2a91ea")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Fib_bubble.pdf", sep=""), width=3+length(Fib_genes)*.15, height=1.5+length(unique(Fib@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Fib, features = Fib_genes, cols=c("grey", "#3249A6")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Fib_bubble_scale.pdf", sep=""), width=5+length(Fib_genes)*.15, height=1.5+length(unique(Fib@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Fib, features = Fib_genes) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Fib_bubble_default.pdf", sep=""), width=5+length(Fib_genes)*.15, height=1.5+length(unique(Fib@meta.data[["Subcelltype"]]))*.4)

plts <- VlnPlot(Fib, features = Fib_genes, split.by='Subcelltype', stack=T, raster=T, cols=c("#050382", "#27D7EF", "#2a91ea"), same.y.lims = T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Fib_vln.pdf", sep=""), width=3+length(Fib_genes)*.4, height=1.5+length(unique(Fib@meta.data[["Subcelltype"]]))*1)
remove(plts)


####Imm sub

#Imm <- readRDS(paste0(curr_dir, "/CellType/Finalcelltype/Imm.RDS"))

Imm_genes <- c("F13A1", "MRC1", "STAB1", "CD163",
               "MERTK", "HCK", "IL4R",
               "DNAJA1", "HSP90AA1", "HSPE1", "HSPD1", "HSPA8", "HSPB1",
               "MKI67", "CENPF", "TOP2A", "PCNA", "EZH2", "GTSE1", "ANLN", "KIF11") %>% unique()

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



####All sub

Org_ChP <- readRDS("/Ex_ChP/Organoid9-5/CellType/Finalcelltype/Organoid_Tsai_EFI.RDS")

Epi_genes <- c("TTR", "HTR2C", "CLIC6", 
               "FRMD4B", "BMPR2", "RNF38", "ERC1",
               "FOXJ1", "DNAH6", "CFAP100", "CFAP69", "SPAG17", "TUBB4B",
               "CLDN1", "CLDN5", "TJP1", "MT-CO1", "MT-CO2", "ATP1A1", "OCLN")
Fibro_genes <- c("LEPR", "CDH11", "STXBP6",
                 "KCNMA1", "SLC4A4", "FOXO1", "KLHL13", "PDGFC", "SLC24A3", "JAM3", "JAM2", 
                 "LAMA4", "CLDN11", "NRG1", "PDE1A", "ABCA9", "BICC1", "ADAMTS12", "ADAMTSL1", "ROBO2", "THBS2",
                 "DCN", "FN1", "COL1A2", "COL27A1", "COL4A2", "COL6A3", "LUM")
Imm_genes <- c("F13A1", "MRC1", "STAB1", "CD163",
               "MERTK", "HCK", "IL4R",
               "DNAJA1", "HSP90AA1", "HSPE1", "HSPD1", "HSPA8", "HSPB1",
               "MKI67", "CENPF", "TOP2A", "PCNA", "EZH2", "GTSE1", "ANLN", "KIF11") %>% unique()
All_genes <- unique(c(Epi_genes, Fibro_genes, Imm_genes))

Idents(Org_ChP) <- Org_ChP@meta.data$Subcelltype
plts <- DotPlot(Org_ChP, features = All_genes, split.by="Subcelltype", cols=c("#CD8500", "#FFA500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#70C25B", "#0E9554", "#266417")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Org_ChP_sub_bubble.pdf", sep=""), width=5+length(All_genes)*.15, height=1.5+length(unique(Org_ChP@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Org_ChP, features = All_genes, split.by="Subcelltype", cols=c("#C48F45", "#C48F45", "#C48F45", "#C48F45", "#3249A6", "#3249A6", "#3249A6", "#26532B", "#26532B", "#26532B")) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Org_ChP_sub_bubble_limit cols.pdf", sep=""), width=5+length(All_genes)*.15, height=1.5+length(unique(Org_ChP@meta.data[["Subcelltype"]]))*.4)

plts <- DotPlot(Org_ChP, features = All_genes) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Org_ChP_sub_bubble_default.pdf", sep=""), width=5+length(All_genes)*.15, height=1.5+length(unique(Org_ChP@meta.data[["Subcelltype"]]))*.4)

plts <- VlnPlot(Org_ChP, features = All_genes, split.by='Subcelltype', stack=T, raster=T, cols=c("#C48F45", "#C48F45", "#C48F45", "#C48F45", "#3249A6", "#3249A6", "#3249A6", "#26532B", "#26532B", "#26532B"), same.y.lims = T) + theme(axis.text.x.bottom = element_text(angle=90), text=element_text(size=7)) + RotatedAxis() + scale_size(range = c(2,5))
plot(plts)
graph2pdf(plts, paste(curr_dir, "/CellType/Publish/Org_ChP_sub_vln.pdf", sep=""), width=3+length(All_genes)*.4, height=1.5+length(unique(Org_ChP@meta.data[["Subcelltype"]]))*1)
remove(plts)



####QC

Idents(Org_ChP) <- Org_ChP@meta.data$Subcelltype
cols <- c("#CD8500", "#FFA500", "#FF4500", "#CD3700", "#050382", "#27D7EF", "#2a91ea", "#70C25B", "#0E9554", "#266417", "#6d207b", "blue", "#ff7450", "#d16120", "#a55d2a", "lightblue", "#ffc274", "#e8a281", "#bba486", "brown")
plts <- VlnPlot(Org_ChP, "nFeature_RNA", pt.size=F, cols=cols)
plot(plts)
ggsave2(paste(curr_dir, "/CellType/Publish/QC_Vln_nFeature_RNA", ".pdf", sep=""), plot = plts, width=4, height=4)
plts <- VlnPlot(Org_ChP, "nCount_RNA", pt.size=F, cols=cols)
plot(plts)
ggsave2(paste(curr_dir, "/CellType/Publish/QC_Vln_nCount_RNA", ".pdf", sep=""), plot = plts, width=4, height=4)
plts <- VlnPlot(Org_ChP, "percent.mt", pt.size=F, cols=cols)
plot(plts)
ggsave2(paste(curr_dir, "/CellType/Publish/percent.mt", ".pdf", sep=""), plot = plts, width=4, height=4)

