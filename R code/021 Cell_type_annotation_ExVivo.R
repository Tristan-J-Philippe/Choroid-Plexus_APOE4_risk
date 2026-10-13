###
### 021 Cell_type_annotation_RNA.R
###
# Purpose: Cell type annotation and cleaning of single nuclei RNA dataset.
# Dependencies:
source("/ChP_RNA/020 Cell_type_annotation_util.R", verbose = FALSE)
curr_dir <- "/ExVivo" #output results to specific folder
set.seed(1234)


##Load and merge samples (note all projIDs replaced with fkID)
ChP_101 <- readRDS(file=paste(curr_dir, "/Sample_preprocess/ChP_101.RDS", sep=""))
...
ChP_112 <- readRDS(file=paste(curr_dir, "/Sample_preprocess/ChP_112.RDS", sep=""))

#merge first half due to limits on total number of cells and computational resources
Ex_ChP <- merge(x=ChP_101, y=c(ChP_102, ... ChP_112))
rm(list = grep("ChP_", ls(), value=TRUE))
gc()

##Run ##Wholeset major celltypes code and determine major celltypes

saveRDS(Ex_ChP, paste(curr_dir, "/CellType/Ex_ChPC.RDS", sep=""))

for(i in unique(Ex_ChP@meta.data$Firstcelltype)){
  sub <- subset(Ex_ChP, Firstcelltype==i)
  saveRDS(sub, paste(curr_dir, "/CellType/Firstcelltype/", i, "C.RDS", sep=""))
  gc()
}

Ex_ChP@meta.data$Lane_projid <- paste(Ex_ChP@meta.data$projid, Ex_ChP@meta.data$Lane, sep="_")
Ex_ChP@meta.data$sample <- Ex_ChP@meta.data$projid
gc()
Ex_ChP@meta.data[is.na(Ex_ChP@meta.data)==TRUE] <- "unk"

L_p <- unique(Ex_ChP@meta.data$Lane_projid)
L_p[order(L_p)] #Check this when making metadata for DEG later



#Celltyping
##Wholeset for each part
###NFSP, Harmonize, CNDR
ChP_genes <- c(known_arach, known_barrier, known_Fibro, known_Mesen, known_Mural, known_periv, known_Epithelial, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, known_dividing, known_epi_neuro_prog, known_epi_prog, known_fibro_prog, known_mglia_prog, known_neuro_prog, known_Neuron, known_neuro_Exc, known_neuro_Inh, known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B)
NFSP_fun(Ex_ChP, assay="soupXcounts", spikein = ChP_genes)
use_dims <- 1:15 #!!!Modify based on elbow plot!!!!
set.seed(1234)

Ex_ChP <- RunHarmony(object = Ex_ChP, group.by.vars = c("Lane"), dims.use = use_dims, assay="soupXcounts", max.iter.harmony = 50, plot_convergence = TRUE) #Correct for batches !!!!
#CNDRfun_quick(Ex_ChP, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony", resolution=0.1)
#Ex_ChP <- FindClusters(Ex_ChP, resolution = c(0.01), assay=assay)
#plts <- clustree(Ex_ChP)
#  plot(plts)
#  graph2png(plts, paste(curr_dir, "/CellType/", "Ex_ChP/Ex_ChP_cluster_tree.png", sep=""), width=14, height=6)
CNDRfun(Ex_ChP, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")


###Choose resolution
res="0.1" #Pick res
Idents(Ex_ChP) <- Ex_ChP@meta.data[[paste("soupXcounts_snn_res.", res, sep="")]]
CNDRdimfun(Ex_ChP, reduction = "harmony", res=res, split.by = c("X10X_batch", "projid", "Lane", "CMO_classification.global", "DoubletFinder_call"))
Ncellfun(Ex_ChP, paste0(curr_dir, "/CellType/Ex_ChP/harmony", res, "/res per Lane.csv", sep=""), "Lane")
Ncellfun(Ex_ChP, paste0(curr_dir, "/CellType/Ex_ChP/harmony", res, "/res per 10X_batch.csv", sep=""), "X10X_batch")
Ncellfun(Ex_ChP, paste0(curr_dir, "/CellType/Ex_ChP/harmony", res, "/res per projid.csv", sep=""), "projid")

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
DimPlot(Ex_ChP, raster=T)

return_markersfun(Ex_ChP, assay="soupXcounts", features = Ex_ChP_soupXcounts_var_feat)###Find markers

###Dividing cells spend more time in S phase
Ex_ChP <- CellCycleScoring(Ex_ChP, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes, set.ident=FALSE)
DimPlot(Ex_ChP, group.by = "Phase", raster=F)
DimPlot(Ex_ChP, split.by = "Phase", raster=F)
Ncellfun(Ex_ChP, paste0(curr_dir, "/CellType/Ex_ChP/CMO per phase.csv", "Phase"))



curr_clust <- 0 #iterate over top 30 markers, number of markers per set and which markers
Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic), na.rm=TRUE)
sum(Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_Mural_1), na.rm=TRUE)
sum(Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Ex_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

Top20_Featuresfun(Ex_ChP, curr_clust)

search_markersfun(Ex_ChP, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Ex_ChP, curr_clust, known_fibro_prog)

search_markersfun(Ex_ChP, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic) %>% unique())
search_markersfun(Ex_ChP, curr_clust, c(known_Mural, known_peric, known_Mural_1) %>% unique())

search_markersfun(Ex_ChP, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Ex_ChP, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Ex_ChP, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Ex_ChP, curr_clust, c(known_light) %>% unique())
search_markersfun(Ex_ChP, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Ex_ChP, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Ex_ChP, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Ex_ChP, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Ex_ChP, curr_clust, known_dividing %>% unique())

search_markersfun(Ex_ChP, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Ex_ChP, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Ex_ChP, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())

search_markersfun(Ex_ChP, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Ex_ChP, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Ex_ChP, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Ex_ChP, "LAMA4", raster=F)


#Ex_ChPC
#"0"="Epithelial", "1"="Fibroblast", "2"="Epithelial", "3"="Epithelial", "4"="Epithelial", "5"="Immune", "6"="Neuro", "7"="Doublet"

###Rename all identities
Ex_ChP <- RenameIdents(Ex_ChP, "0"="Epi_2", "1"="Epi_1", "2"="Fibro", "3"="Epi_1", "4"="Imm", "5"="Endo", "6"="Neuro", "7"="Oligo", "8"="Astro", "9"="Neuro", "10"="Oligo")
Ex_ChP@meta.data$Firstcelltype <- Idents(Ex_ChP)
saveRDS(Ex_ChP, paste(curr_dir, "/CellType/Ex_ChP.RDS", sep=""))


###Subtypes
for(i in unique(Ex_ChP@meta.data$Firstcelltype)){
  sub <- subset(Ex_ChP, Firstcelltype==i)
  saveRDS(sub, paste(curr_dir, "/CellType/Firstcelltype/", i, "C.RDS", sep=""))
}



##Epi_1
###PCA
Epi_1A <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/Epi_1A.RDS"))
Epi_1B <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/Epi_1B.RDS"))
Epi_1C <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/Epi_1C.RDS"))
Epi_1 <- merge(Epi_1A, Epi_1B, Epi_1C) #merge  batches
rm(Epi_1A, Epi_1B, Epi_1C)
gc()


ChP_genes <- c(known_arach, known_barrier, known_Fibro, known_Mesen, known_Mural, known_periv, known_Epithelial, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, known_dividing, known_epi_neuro_prog, known_epi_prog, known_fibro_prog, known_mglia_prog, known_neuro_prog, known_Neuron, known_neuro_Exc, known_neuro_Inh, known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B)

NFSP_fun(Epi, assay="soupXcounts")#, spikein = ChP_genes)

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Epi <- RunHarmony(object = Epi, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun_quick(Epi, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")
#Epi <- FindClusters(Epi, resolution = c(0.07), assay="soupXcounts")
#Epi <- FindClusters(Epi, resolution = c(0.09), assay="soupXcounts")
#Epi <- FindClusters(Epi, resolution = c(0.3), assay="soupXcounts")
#plts <- clustree(Epi)
#plot(plts)
#graph2png(plts, paste(curr_dir, "/CellType/", "Epi/Ex_ChP_cluster_tree.png", sep=""), width=14, height=6)

#Pick res
Idents(Epi) <- Epi@meta.data[["soupXcounts_snn_res.0.1"]]
CNDRdimfun(Epi, reduction = "harmony", res="_0.1", split.by = c("X10X_batch", "library_batch", "CMO_classification.global", "Subcelltype", "Majorcelltype"))
Ncellfun(Epi, paste0(curr_dir, "/CellType/Epi/Final per AD.csv", sep=""), "X10X_batch")
Ncellfun(Epi, paste0(curr_dir, "/CellType/Epi/Final per AD.csv", sep=""), "Cog_Path")
#return markers
return_markersfun(Epi, assay="soupXcounts", features = Epi_soupXcounts_var_feat, num_clust=0:3)

Featurefun(Epi, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "GPX3", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS", "XIST", "XACT"), "KeyMarkers",  assay="soupXcounts")


###remove doublets

Featurefun(Epi, c("HTR2C", "TTR", "CLIC6", "AQP1", "KRT8", "KRT18", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="soupXcounts")

#remove batch specific clusters with few cells on first pass
Epi <- subset(Epi, soupXcounts_snn_res.0.03=="0" | soupXcounts_snn_res.0.03=="1")
Epi <- subset(Epi, soupXcounts_snn_res.0.3=="13" | soupXcounts_snn_res.0.3=="14" | soupXcounts_snn_res.0.3=="15" | soupXcounts_snn_res.0.3=="16", invert=T)

Epi <- subset(Epi, CMO_classification.global=="Singlet")
gc()

#Remove cells that don't belong
#Epi <- subset(Epi, HTR2C < 4)
#Epi <- subset(Epi, TTR < 4)
#Epi <- subset(Epi, CLIC6 < 4)
#Epi <- subset(Epi, AQP1 < 3)
#Epi <- subset(Epi, CFAP54 < 3)
#Epi <- subset(Epi, FAM227A < 3)
#Epi <- subset(Epi, ARL13B < 3)
#Epi <- subset(Epi, TRPM3 < 3)
Epi <- subset(Epi, SLC4A4 < 3)
Epi <- subset(Epi, LAMA2 < 3)
Epi <- subset(Epi, LAMA4 < 3)
#Epi <- subset(Epi, CLDN1 < 3)
Epi <- subset(Epi, CLDN11 < 3)
Epi <- subset(Epi, PRICKLE1 < 3)
Epi <- subset(Epi, CEMIP < 3)
Epi <- subset(Epi, FOXP2 < 3)
Epi <- subset(Epi, COL6A3 < 3)
Epi <- subset(Epi, SLC6A13 < 3)
Epi <- subset(Epi, MECOM < 3)
Epi <- subset(Epi, VWF < 3)
Epi <- subset(Epi, TAGLN < 3)
Epi <- subset(Epi, MYH11 < 3)
Epi <- subset(Epi, NOTCH3 < 3)
Epi <- subset(Epi, CADM2 < 3)
Epi <- subset(Epi, MAP2 < 3)
Epi <- subset(Epi, GRIK2 < 3)
#Epi <- subset(Epi, CLU < 3)
#Epi <- subset(Epi, AQP4 < 3)
Epi <- subset(Epi, MOBP < 3)
Epi <- subset(Epi, MBP < 3)
Epi <- subset(Epi, CDK19 < 3)
Epi <- subset(Epi, CD68 < 3)
Epi <- subset(Epi, F13A1 < 3)
Epi <- subset(Epi, STAB1 < 3)
Epi <- subset(Epi, CD163 < 3)
Epi <- subset(Epi, LYVE1 < 3)
Epi <- subset(Epi, P2RY12 < 3)
Epi <- subset(Epi, THEMIS < 3)


###Label

#Second pass
Idents(Epi) <- Epi@meta.data[["soupXcounts_snn_res.0.3"]]
Epi <- RenameIdents(Epi, "0"="Epi_1a", "1"="Epi_1b", "2"="Epi_2a", "3"="Epi_2b")
Epi@meta.data$Subcelltype <- Idents(Epi)
Epi <- subset(Epi, Subcelltype=="Epi_1a" | Subcelltype=="Epi_1b" | Subcelltype=="Epi_2a" | Subcelltype=="Epi_2b")
saveRDS(Epi, paste0(curr_dir, "/CellType/Thirdcelltype/Epi_old.RDS"))

#Last pass
Idents(Epi) <- Epi@meta.data[["soupXcounts_snn_res.0.1"]]
Epi <- RenameIdents(Epi, "0"="Epi_1b", "1"="Epi_1a", "2"="Epi_2a", "3"="Epi_2b", "4"="Doublet", "5"="Doublet", "6"="Doublet", "7"="Doublet", "8"="Doublet", "9"="Doublet", "10"="Doublet")
Epi@meta.data$Subcelltype <- Idents(Epi)
DimPlot(Epi, raster=F, split.by="Subcelltype")
DimPlot(Epi, raster=F)


###Remove interdoublets

Idents(Epi) <- Epi@meta.data$Subcelltype
Epi <- CellSelector(DimPlot(Epi), Epi, ident = "Epi1")
Epi_1 <- subset(Epi, ident="Epi1")
Epi_1 <- subset(Epi_1, Subcelltype=="Epi_1a")


Idents(Epi) <- Epi@meta.data$Subcelltype
Epi <- CellSelector(DimPlot(Epi), Epi, ident = "Epi2")
Epi_2 <- subset(Epi, Subcelltype=="Epi_1b")
Epi_2 <- subset(Epi_2, ident="Epi2")


Idents(Epi) <- Epi@meta.data$Subcelltype
Epi <- CellSelector(DimPlot(Epi), Epi, ident = "Epi3")
Epi_3 <- subset(Epi, Subcelltype=="Epi_2a")
Epi_3 <- subset(Epi_3, ident="Epi3")


Idents(Epi) <- Epi@meta.data$Subcelltype
Epi <- CellSelector(DimPlot(Epi), Epi, ident = "Epi4")
Epi_4 <- subset(Epi, Subcelltype=="Epi_2b")
Epi_4 <- subset(Epi_4, ident="Epi4")

Epi_1_meta <- Epi_1@meta.data
Epi_2_meta <- Epi_2@meta.data
Epi_3_meta <- Epi_3@meta.data
Epi_4_meta <- Epi_4@meta.data
merged_meta <- rbind(Epi_1_meta, Epi_2_meta, Epi_3_meta, Epi_4_meta)

Epi@meta.data$Subcelltype <- NULL
Epi_meta <- Epi@meta.data

Epi_meta <- rownames_to_column(Epi_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$Subcelltype)
colnames(merged_meta) <- c("rowname", "Subcelltype")
Epi_meta <- left_join(Epi_meta, merged_meta, by="rowname")
Epi_meta <- column_to_rownames(Epi_meta, "rowname")
Epi_meta$Subcelltype <- as.character(Epi_meta$Subcelltype)
Epi_meta[is.na(Epi_meta)==TRUE] <- "Doublet"
Epi@meta.data <- Epi_meta
Epi <- subset(Epi, Subcelltype=="Doublet", invert=T)

#Reorder
Idents(Epi) <- Epi@meta.data[["Subcelltype"]] 
levels(Epi) <- c("Epi_1a", "Epi_1b", "Epi_2a", "Epi_2b")
Epi@meta.data[["Subcelltype"]] <- Idents(Epi)

Epi <- RunUMAP(object = Epi, dims = 1:20, assay = "soupXcounts", reduction =  "harmony")
DimPlot(Epi, raster=F)

saveRDS(Epi, paste0(curr_dir, "/CellType/Finalcelltype/Epi.RDS"))



curr_clust <- 1
Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:50]
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_Mural_1), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Epi, curr_clust)

search_markersfun(Epi, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Epi, curr_clust, known_fibro_prog)

search_markersfun(Epi, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic) %>% unique())
search_markersfun(Epi, curr_clust, c(known_Mural, known_peric, known_VSM) %>% unique())

search_markersfun(Epi, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Epi, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Epi, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Epi, curr_clust, c(known_light) %>% unique())
search_markersfun(Epi, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Epi, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Epi, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Epi, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Epi, curr_clust, known_dividing %>% unique())

search_markersfun(Epi, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Epi, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Epi, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Epi, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Epi, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Epi, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Epi, c("TTR", "CFAP46", "MT-CO2"), raster=F)

##Fibro
###PCA

FibroC <- readRDS("/Ex_ChP/Ex_ChP7/CellType/Firstcelltype/FibroblastA.RDS")
FibroAB <- readRDS("/Ex_ChP/Ex_ChP7/CellType/Secondcelltype_og/Fibro.RDS")
FibroAB$Lane <- FibroAB$LANE
FibroAB@meta.data[["DoubletFinder_DFclusters"]] <- NULL
MuralAB <- readRDS("/Ex_ChP/Ex_ChP7/CellType/Secondcelltype_og/Mural.RDS")
MuralAB$Lane <- MuralAB$LANE
MuralAB@meta.data[["DoubletFinder_DFclusters"]] <- NULL

Fibro <- merge(FibroAB, y=c(FibroC, MuralAB))
rm(FibroAB, FibroC, MuralAB)
gc()

ChP_genes <- c(known_arach, known_barrier, known_Fibro, known_Mesen, known_Mural, known_periv, known_Epithelial, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, known_dividing, known_epi_neuro_prog, known_epi_prog, known_fibro_prog, known_mglia_prog, known_neuro_prog, known_Neuron, known_neuro_Exc, known_neuro_Inh, known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B)

NFSP_fun(Fibro, assay="soupXcounts")#, spikein= ChP_genes)

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Fibro <- RunHarmony(object = Fibro, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun_quick(Fibro, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")  #!!!!pca or harmony

#Pick res
#Idents(Fibro) <- Fibro@meta.data$Subcelltype
Idents(Fibro) <- Fibro@meta.data[["soupXcounts_snn_res.0.1"]]
CNDRdimfun(Fibro, reduction = "harmony", res="0.1", split.by = c("X10X_batch", "library_batch", "Subcelltype", "Majorcelltype", "Secondcelltype", "CMO_classification.global"))
#return markers
return_markersfun(Fibro, assay="soupXcounts", features = Fibro_soupXcounts_var_feat)
saveRDS(Fibro, paste0(curr_dir, "/CellType/Firstcelltype/Fibro_combined_spiked.RDS"))



curr_clust <- 2
Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_VSM), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Fibro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Fibro, curr_clust)

search_markersfun(Fibro, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, "CLDN11", "TJP1", "TJP2") %>% unique())
search_markersfun(Fibro, curr_clust, known_fibro_prog)

search_markersfun(Fibro, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic) %>% unique())
search_markersfun(Fibro, curr_clust, c(known_Mural, known_peric, known_VSM) %>% unique())

search_markersfun(Fibro, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Fibro, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Fibro, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Fibro, curr_clust, c(known_light) %>% unique())
search_markersfun(Fibro, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Fibro, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Fibro, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Fibro, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Fibro, curr_clust, c(known_dividing, dividing_Org_Wang) %>% unique())

search_markersfun(Fibro, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog, prog_Org_Wang) %>% unique())

search_markersfun(Fibro, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Fibro, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Fibro, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Fibro, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Fibro, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Fibro, "CLDN11", raster=F)

Fib_2: LAMA2 TRPM3: 1 keep 7
Fib_1: LAMA2: 3, 11
Fib_3: SLC4A4: 0, 4 keep 2, 5 10
Fib_4: SLC4A4 TRPM3: 6, 9
Mural: 8

remove doublets and Mural res 0.3
Mural: 8
Doublets: 7, 2, 5, 10

2nd pass
Fib_2: 0
Fib_4: 7
Fib_3: 1,2,3,6
Fib_1: 4,5
8 Doublet

###remove doublets

Featurefun(Fibro, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS", "XIST", "XACT"), "Test2",  assay="soupXcounts")

Fibro <- subset(Fibro, CMO_classification.global=="Singlet")

#Remove cells that don't belong
Fibro <- subset(Fibro, HTR2C < 3)
Fibro <- subset(Fibro, TTR < 3)
Fibro <- subset(Fibro, CLIC6 < 3)
#Fibro <- subset(Fibro, AQP1 < 3)
Fibro <- subset(Fibro, CFAP54 < 3)
Fibro <- subset(Fibro, FAM227A < 3)
Fibro <- subset(Fibro, ARL13B < 3)
#Fibro <- subset(Fibro, TRPM3 < 3)
#Fibro <- subset(Fibro, SLC4A4 < 3)
#Fibro <- subset(Fibro, LAMA2 < 3)
#Fibro <- subset(Fibro, LAMA4 < 3)
#Fibro <- subset(Fibro, CLDN1 < 3)
#Fibro <- subset(Fibro, CLDN11 < 3)
#Fibro <- subset(Fibro, PRICKLE1 < 3)
#Fibro <- subset(Fibro, CEMIP < 3)
#Fibro <- subset(Fibro, FOXP2 < 3)
#Fibro <- subset(Fibro, COL6A3 < 3)
#Fibro <- subset(Fibro, SLC6A13 < 3)
Fibro <- subset(Fibro, MECOM < 3)
Fibro <- subset(Fibro, VWF < 3)
Fibro <- subset(Fibro, TAGLN < 3)
Fibro <- subset(Fibro, MYH11 < 3)
Fibro <- subset(Fibro, NOTCH3 < 3)
Fibro <- subset(Fibro, CADM2 < 3)
Fibro <- subset(Fibro, MAP2 < 3)
Fibro <- subset(Fibro, GRIK2 < 3)
#Fibro <- subset(Fibro, CLU < 3)
#Fibro <- subset(Fibro, AQP4 < 3)
Fibro <- subset(Fibro, MOBP < 3)
Fibro <- subset(Fibro, MBP < 3)
Fibro <- subset(Fibro, CDK19 < 3)
Fibro <- subset(Fibro, CD68 < 3)
Fibro <- subset(Fibro, F13A1 < 3)
Fibro <- subset(Fibro, STAB1 < 3)
Fibro <- subset(Fibro, CD163 < 3)
Fibro <- subset(Fibro, LYVE1 < 3)
Fibro <- subset(Fibro, P2RY12 < 3)
Fibro <- subset(Fibro, THEMIS < 3)


###Label

#Firstpass
Fibro <- subset(Fibro, percent.mt<5)

Idents(Fibro) <- Fibro@meta.data[["soupXcounts_snn_res.0.1"]]
Fibro <- RenameIdents(Fibro, "0"="Fib_1", "1"="Fib_2", "2"="Fib_3", "3"="Mural")
Fibro@meta.data$Secondcelltype <- Fibro@active.ident
Mural <- subset(Fibro, Secondcelltype=="Mural")
Fibro <- subset(Fibro, Secondcelltype=="Fib_1" | Secondcelltype=="Fib_2" | Secondcelltype=="Fib_3")

saveRDS(Fibro, paste0(curr_dir, "/CellType/Secondcelltype/Fibro.RDS"))
saveRDS(Mural, paste0(curr_dir, "/CellType/Secondcelltype/Mural.RDS"))

#Second pass
Idents(Fibro) <- Fibro@meta.data[["soupXcounts_snn_res.0.1"]]
Fibro <- RenameIdents(Fibro, "0"="Fib_1", "1"="Fib_2", "2"="Fib_3", "3"="Fib_1")
Idents(Fibro) -> Fibro@meta.data$Subcelltype
DimPlot(Fibro, raster=F)


###Remove interdoublets

Idents(Fibro) <- Fibro@meta.data$Subcelltype
Fibro <- CellSelector(DimPlot(Fibro), Fibro, ident = "Fib1")
Fib_1 <- subset(Fibro, Subcelltype=="Fib_1")
Fib_1 <- subset(Fib_1, ident="Fib1")

Idents(Fibro) <- Fibro@meta.data$Subcelltype
Fibro <- CellSelector(DimPlot(Fibro), Fibro, ident = "Fib2")
Fib_2 <- subset(Fibro, Subcelltype=="Fib_2")
Fib_2 <- subset(Fib_2, ident="Fib2")

Idents(Fibro) <- Fibro@meta.data$Subcelltype
Fibro <- CellSelector(DimPlot(Fibro), Fibro, ident = "Fib3")
Fib_3 <- subset(Fibro, Subcelltype=="Fib_3")
Fib_3 <- subset(Fib_3, ident="Fib3")

Fib_1_meta <- Fib_1@meta.data
Fib_2_meta <- Fib_2@meta.data
Fib_3_meta <- Fib_3@meta.data
merged_meta <- rbind(Fib_1_meta, Fib_2_meta, Fib_3_meta)
rm(Fib_1, Fib_2, Fib_3)
gc()

Fibro@meta.data$Subcelltype <- NULL
Fibro_meta <- Fibro@meta.data

Fibro_meta <- rownames_to_column(Fibro_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$Subcelltype)
colnames(merged_meta) <- c("rowname", "Subcelltype")
Fibro_meta <- left_join(Fibro_meta, merged_meta, by="rowname")
Fibro_meta <- column_to_rownames(Fibro_meta, "rowname")
Fibro_meta$Subcelltype <- as.character(Fibro_meta$Subcelltype)
Fibro_meta[is.na(Fibro_meta)==TRUE] <- "Doublet"
Fibro@meta.data <- Fibro_meta
Fibro <- subset(Fibro, Subcelltype=="Doublet", invert=T)

#Reorder
Idents(Fibro) <- Fibro@meta.data[["Subcelltype"]] 
levels(Fibro) <- c("Fib_1", "Fib_2", "Fib_3")
Fibro@meta.data[["Subcelltype"]] <- Idents(Fibro)


Fibro <- RunUMAP(object = Fibro, dims = 1:20, assay = "soupXcounts", reduction =  "harmony", min.dist = 1.4)
DimPlot(Fibro, raster=F)

saveRDS(Fibro, paste0(curr_dir, "/CellType/Finalcelltype/Fibro.RDS"))



##Endo
###PCA

EndoC <- readRDS("/Ex_ChP/Ex_ChP7/CellType/Firstcelltype/EndothelialA.RDS")
EndoAB <- readRDS("/Ex_ChP/Ex_ChP7/CellType/Secondcelltype_og/Endo.RDS")
EndoAB$Lane <- EndoAB$LANE
EndoAB@meta.data[["DoubletFinder_DFclusters"]] <- NULL

Endo <- merge(EndoC, Endo)

Endo <- readRDS("/Ex_ChP/Ex_ChP7/CellType/Secondcelltype/Endo.RDS")

#remove doublets
Endo <- subset(Endo, CMO_classification.global=="Singlet")
Endo <- subset(Endo, percent.mt<5)

#Featurefun(Endo, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS", "XIST", "XACT"), "Test2",  assay="soupXcounts")

#Remove cells that don't belong
Endo <- subset(Endo, HTR2C < 3)
Endo <- subset(Endo, TTR < 3)
Endo <- subset(Endo, CLIC6 < 3)
#Endo <- subset(Endo, AQP1 < 3)
Endo <- subset(Endo, CFAP54 < 3)
Endo <- subset(Endo, FAM227A < 3)
Endo <- subset(Endo, ARL13B < 3)
#Endo <- subset(Endo, TRPM3 < 3)
Endo <- subset(Endo, SLC4A4 < 3)
Endo <- subset(Endo, LAMA2 < 3)
Endo <- subset(Endo, LAMA4 < 3)
Endo <- subset(Endo, CLDN1 < 3)
Endo <- subset(Endo, CLDN11 < 3)
Endo <- subset(Endo, PRICKLE1 < 3)
Endo <- subset(Endo, CEMIP < 3)
Endo <- subset(Endo, FOXP2 < 3)
Endo <- subset(Endo, COL6A3 < 3)
Endo <- subset(Endo, SLC6A13 < 3)
#Endo <- subset(Endo, MECOM < 3)
#Endo <- subset(Endo, VWF < 3)
Endo <- subset(Endo, TAGLN < 3)
Endo <- subset(Endo, MYH11 < 3)
Endo <- subset(Endo, NOTCH3 < 3)
Endo <- subset(Endo, CADM2 < 3)
Endo <- subset(Endo, MAP2 < 3)
Endo <- subset(Endo, GRIK2 < 3)
#Endo <- subset(Endo, CLU < 3)
#Endo <- subset(Endo, AQP4 < 3)
Endo <- subset(Endo, MOBP < 3)
Endo <- subset(Endo, MBP < 3)
Endo <- subset(Endo, CDK19 < 3)
Endo <- subset(Endo, CD68 < 3)
Endo <- subset(Endo, F13A1 < 3)
Endo <- subset(Endo, STAB1 < 3)
Endo <- subset(Endo, CD163 < 3)
Endo <- subset(Endo, LYVE1 < 3)
Endo <- subset(Endo, P2RY12 < 3)
Endo <- subset(Endo, THEMIS < 3)



NFSP_fun(Endo, assay="soupXcounts")

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Endo <- RunHarmony(object = Endo, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Endo, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")  #!!!!pca or harmony

#Pick res
Idents(Endo) <- Endo@meta.data[["soupXcounts_snn_res.0.07"]]
CNDRdimfun(Endo, reduction = "harmony", res="_0.07", split.by = c("X10X_batch", "library_batch", "Subcelltype", "Majorcelltype", "CMO_classification.global"))
#return markers
return_markersfun(Endo, assay="soupXcounts", features = Endo_soupXcounts_var_feat)
saveRDS(Endo, paste0(curr_dir, "/CellType/Secondcelltype/Endo.RDS"))

Featurefun(Endo, c("HTR2C", "TTR", "CLIC6", "AQP1", "TUBB4B", "FOXJ1", "ARL13B", "KRT18", "COX2", "TRPM3", "VIM", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "FOXP2", "COL6A3", "SLC6A13", "LUM", "DCN", "FN1", "ACTA2", "KCNMA1", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "CD3E", "HSPE1", "HSPH1"), "Key Markers",  assay="soupXcounts")



curr_clust <- 2
Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic), na.rm=TRUE)
sum(Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_VSM), na.rm=TRUE)
sum(Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Endo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Endo, curr_clust)

search_markersfun(Endo, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
#search_markersfun(Endo, curr_clust, known_Fibro_prog)

search_markersfun(Endo, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic) %>% unique())
search_markersfun(Endo, curr_clust, c(known_veinous) %>% unique())
search_markersfun(Endo, curr_clust, c(known_capilary) %>% unique())
search_markersfun(Endo, curr_clust, c(known_arterial) %>% unique())
search_markersfun(Endo, curr_clust, c(known_lymphatic) %>% unique())

search_markersfun(Endo, curr_clust, c(known_Mural, known_peric, known_VSM) %>% unique())


search_markersfun(Endo, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Endo, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Endo, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Endo, curr_clust, c(known_light) %>% unique())
search_markersfun(Endo, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Endo, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Endo, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Endo, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Endo, curr_clust, known_dividing %>% unique())

search_markersfun(Endo, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Endo, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Endo, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Endo, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Endo, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Endo, curr_clust, c(Cell_Stress, whole_cells) %>% unique())

search_markersfun(Endo, 5, c("CLDN5", "IL1R1", "VEGFC", "RGCC", "PLVAP", "SGK1"))
FeaturePlot(Endo, c("CLDN5", "IL1R1"), raster=F)
FeaturePlot(Endo, c("VEGFC", "RGCC"), raster=F)
FeaturePlot(Endo, c("PLVAP", "SGK1"), raster=F)

VlnPlot(Endo, c("CLDN5", "IL1R1", "VEGFC", "RGCC", "PLVAP", "SGK1"))



###Label

Endo <- RenameIdents(Endo, "0"="Endo_3", "1"="Endo_1", "2"="Endo_3", "3"="Endo_3", "4"="Endo_2") #First pass
Endo@meta.data$Subcelltype <- Endo@active.ident

DimPlot(Endo)

saveRDS(Endo, paste0(curr_dir, "/CellType/Finalcelltype/Endo.RDS"))


###Remove interdoublets

Endo <- RenameIdents(Endo, "0"="Endo_3", "1"="Endo_1", "2"="Endo_2") #Last pass
Endo@meta.data$Subcelltype <- Endo@active.ident

Idents(Endo) <- Endo@meta.data$Subcelltype
Endo <- CellSelector(DimPlot(Endo), Endo, ident = "Cap")
Cap <- subset(Endo, ident="Cap")
Cap <- subset(Cap, Subcelltype=="Endo_3")


Idents(Endo) <- Endo@meta.data$Subcelltype
Endo <- CellSelector(DimPlot(Endo), Endo, ident = "Art")
Art <- subset(Endo, Subcelltype=="Endo_1")
Art <- subset(Art, ident="Art")


Idents(Endo) <- Endo@meta.data$Subcelltype
Endo <- CellSelector(DimPlot(Endo), Endo, ident = "Ven")
Ven <- subset(Endo, Subcelltype=="Endo_2")
Ven <- subset(Ven, ident="Ven")


Cap_meta <- Cap@meta.data
Art_meta <- Art@meta.data
Ven_meta <- Ven@meta.data
merged_meta <- rbind(Cap_meta, Art_meta, Ven_meta)

Endo@meta.data$Subcelltype <- NULL
Endo_meta <- Endo@meta.data

Endo_meta <- rownames_to_column(Endo_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$Subcelltype)
colnames(merged_meta) <- c("rowname", "Subcelltype")
Endo_meta <- left_join(Endo_meta, merged_meta, by="rowname")
Endo_meta <- column_to_rownames(Endo_meta, "rowname")
Endo_meta$Subcelltype <- as.character(Endo_meta$Subcelltype)
Endo_meta[is.na(Endo_meta)==TRUE] <- "Doublet"
Endo@meta.data <- Endo_meta
Endo <- subset(Endo, Subcelltype=="Doublet", invert=T)

#Reorder
Idents(Endo) <- Endo@meta.data[["Subcelltype"]] 
levels(Endo) <- c("Endo_1", "Endo_2", "Endo_3")
Endo@meta.data[["Subcelltype"]] <- Idents(Endo)
DimPlot(Endo, raster=F)

Endo <- RunUMAP(object = Endo, dims = 1:20, assay = "soupXcounts", reduction =  "harmony", min.dist = 1.2)
saveRDS(Endo, paste0(curr_dir, "/CellType/Finalcelltype/Endo.RDS"))

Ncellfun(Endo, paste0(curr_dir, "/CellType/Endo/Final per APOE.csv", sep=""), "apoe_genotype")
Ncellfun(Endo, paste0(curr_dir, "/CellType/Endo/Final per AD.csv", sep=""), "Cog_Path")
Ncellfun(Endo, paste0(curr_dir, "/CellType/Endo/Final per projid.csv", sep=""), "projid")





##Mural
###remove doublets

Mural <- readRDS("/Ex_ChP/Ex_ChP7/CellType/Secondcelltype/Mural.RDS")

#remove doublets
Mural <- subset(Mural, CMO_classification.global=="Singlet")
Mural <- subset(Mural, percent.mt<5)

Featurefun(Mural, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS", "XIST", "XACT"), "Test2",  assay="soupXcounts")

#Remove cells that don't belong
Mural <- subset(Mural, HTR2C < 3)
Mural <- subset(Mural, TTR < 3)
Mural <- subset(Mural, CLIC6 < 3)
#Mural <- subset(Mural, AQP1 < 3)
Mural <- subset(Mural, CFAP54 < 3)
Mural <- subset(Mural, FAM227A < 3)
Mural <- subset(Mural, ARL13B < 3)
#Mural <- subset(Mural, TRPM3 < 3)
Mural <- subset(Mural, SLC4A4 < 3)
#Mural <- subset(Mural, LAMA2 < 3)
#Mural <- subset(Mural, LAMA4 < 3)
Mural <- subset(Mural, CLDN1 < 3)
Mural <- subset(Mural, CLDN11 < 3)
Mural <- subset(Mural, PRICKLE1 < 3)
Mural <- subset(Mural, CEMIP < 3)
Mural <- subset(Mural, FOXP2 < 3)
Mural <- subset(Mural, COL6A3 < 3)
Mural <- subset(Mural, SLC6A13 < 3)
Mural <- subset(Mural, MECOM < 3)
Mural <- subset(Mural, VWF < 3)
#Mural <- subset(Mural, TAGLN < 3)
#Mural <- subset(Mural, MYH11 < 3)
#Mural <- subset(Mural, NOTCH3 < 3)
Mural <- subset(Mural, CADM2 < 3)
Mural <- subset(Mural, MAP2 < 3)
Mural <- subset(Mural, GRIK2 < 3)
#Mural <- subset(Mural, CLU < 3)
#Mural <- subset(Mural, AQP4 < 3)
Mural <- subset(Mural, MOBP < 3)
Mural <- subset(Mural, MBP < 3)
Mural <- subset(Mural, CDK19 < 3)
Mural <- subset(Mural, CD68 < 3)
Mural <- subset(Mural, F13A1 < 3)
Mural <- subset(Mural, STAB1 < 3)
Mural <- subset(Mural, CD163 < 3)
Mural <- subset(Mural, LYVE1 < 3)
Mural <- subset(Mural, P2RY12 < 3)
Mural <- subset(Mural, THEMIS < 3)


###PCA

NFSP_fun(Mural, assay="soupXcounts")

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Mural <- RunHarmony(object = Mural, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Mural, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")  #!!!!pca or harmony

#Pick res
Idents(Mural) <- Mural@meta.data[["soupXcounts_snn_res.0.1"]]
CNDRdimfun(Mural, reduction = "harmony", res="_0.1", split.by = c("X10X_batch", "library_batch", "Subcelltype", "Majorcelltype", "Firstcelltype", "CMO_classification.global", "Cog_Path"))
#return markers
return_markersfun(Mural, assay="soupXcounts", features = Mural_soupXcounts_var_feat)



curr_clust <- 2
Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic), na.rm=TRUE)
sum(Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_VSM), na.rm=TRUE)
sum(Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Mural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Mural, curr_clust)

search_markersfun(Mural, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Mural, curr_clust, known_Fibro_prog)

search_markersfun(Mural, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_Endo_1, known_lymphatic) %>% unique())
search_markersfun(Mural, curr_clust, c(known_Mural) %>% unique())
search_markersfun(Mural, curr_clust, c(known_peric) %>% unique())
search_markersfun(Mural, curr_clust, c(known_VSM) %>% unique())


search_markersfun(Mural, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Mural, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Mural, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Mural, curr_clust, c(known_light) %>% unique())
search_markersfun(Mural, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Mural, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Mural, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Mural, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Mural, curr_clust, known_dividing %>% unique())

search_markersfun(Mural, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Mural, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Mural, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Mural, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Mural, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Mural, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
Featurefun(Mural, c("ENPP2","PLXDC2","PALMD", "MID1", "CEMIP", "GLIS3", "PDZD2", 
                    "COL4A4", "CRISPLD2", "ADAMTS9", "TLE1", "PIPNC1", "ADCY3", "HIPK2", "COL18A1", "SCL22A3",
                    "SOX6", "MYO1D", "MYH11", "ACTA2", "PRUNE2", "RYR2", "TPM1", "LMOD1"), "Mural")

0 : Mural_2
1 : Mural_1
2 : Mural_3
3 : Mural_2
4 : Mural_2


###Label

Mural <- RenameIdents(Mural, "0"="Mural_3", "1"="Mural_1", "2"="Mural_2") #last pass
Mural@meta.data$Subcelltype <- Mural@active.ident
DimPlot(Mural)
saveRDS(Mural, paste0(curr_dir, "/CellType/Finalcelltype/Mural.RDS"))



###Remove interdoublets

Idents(Mural) <- Mural@meta.data$Subcelltype
Mural <- CellSelector(DimPlot(Mural), Mural, ident = "Peri1")
Peri1 <- subset(Mural, ident="Peri1")
Peri1 <- subset(Peri1, Subcelltype=="Mural_2")

Idents(Mural) <- Mural@meta.data$Subcelltype
Mural <- CellSelector(DimPlot(Mural), Mural, ident = "Peri2")
Peri2 <- subset(Mural, ident="Peri2")
Peri2 <- subset(Peri2, Subcelltype=="Mural_3")

Idents(Mural) <- Mural@meta.data$Subcelltype
Mural <- CellSelector(DimPlot(Mural), Mural, ident = "Mural_1s")
Mural_1 <- subset(Mural, ident="Mural_1s")
Mural_1 <- subset(Mural_1, Subcelltype=="Mural_1")


Peri1_meta <- Peri1@meta.data
Peri2_meta <- Peri2@meta.data
Mural_1_meta <- Mural_1@meta.data
merged_meta <- rbind(Peri1_meta, Peri2_meta, Mural_1_meta)

Mural@meta.data$Subcelltype <- NULL
Mural_meta <- Mural@meta.data

Mural_meta <- rownames_to_column(Mural_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$Subcelltype)
colnames(merged_meta) <- c("rowname", "Subcelltype")
Mural_meta <- left_join(Mural_meta, merged_meta, by="rowname")
Mural_meta <- column_to_rownames(Mural_meta, "rowname")
Mural_meta$Subcelltype <- as.character(Mural_meta$Subcelltype)
Mural_meta[is.na(Mural_meta)==TRUE] <- "Doublet"
Mural@meta.data <- Mural_meta
Mural <- subset(Mural, Subcelltype=="Doublet", invert=T)

#Reorder
Idents(Mural) <- Mural@meta.data[["Subcelltype"]] 
levels(Mural) <- c("Mural_1", "Mural_2", "Mural_3")
Mural@meta.data[["Subcelltype"]] <- Idents(Mural)
DimPlot(Mural, raster=F)

Mural <- RunUMAP(object = Mural, dims = 1:20, assay = "soupXcounts", reduction =  "harmony", min.dist = 1.5)
saveRDS(Mural, paste0(curr_dir, "/CellType/Finalcelltype/Mural.RDS"))

Ncellfun(Mural, paste0(curr_dir, "/CellType/Mural/Final per APOE.csv", sep=""), "apoe_genotype")
Ncellfun(Mural, paste0(curr_dir, "/CellType/Mural/Final per AD.csv", sep=""), "Cog_Path")
Ncellfun(Mural, paste0(curr_dir, "/CellType/Mural/Final per projid.csv", sep=""), "projid")



##Imm
###remove doublets

#remove doublets
Imm <- subset(Imm, CMO_classification.global=="Singlet")
Imm <- subset(Imm, percent.mt<5)

Featurefun(Imm, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS", "XIST", "XACT"), "Imm_markers",  assay="soupXcounts")
Featurefun(Imm, c("HSPH1"), "Imm_markers",  assay="soupXcounts")
#Remove cells that don't belong
Imm <- subset(Imm, HTR2C < 4)
Imm <- subset(Imm, TTR < 4)
Imm <- subset(Imm, CLIC6 < 4)
#Imm <- subset(Imm, AQP1 < 3)
Imm <- subset(Imm, CFAP54 < 3)
Imm <- subset(Imm, FAM227A < 3)
Imm <- subset(Imm, ARL13B < 3)
Imm <- subset(Imm, TRPM3 < 3)
Imm <- subset(Imm, SLC4A4 < 3)
Imm <- subset(Imm, LAMA2 < 3)
Imm <- subset(Imm, LAMA4 < 3)
Imm <- subset(Imm, CLDN1 < 3)
Imm <- subset(Imm, CLDN11 < 3)
Imm <- subset(Imm, PRICKLE1 < 3)
Imm <- subset(Imm, CEMIP < 3)
Imm <- subset(Imm, FOXP2 < 3)
Imm <- subset(Imm, COL6A3 < 3)
Imm <- subset(Imm, SLC6A13 < 3)
Imm <- subset(Imm, MECOM < 3)
Imm <- subset(Imm, VWF < 3)
Imm <- subset(Imm, TAGLN < 3)
Imm <- subset(Imm, MYH11 < 3)
Imm <- subset(Imm, NOTCH3 < 3)
Imm <- subset(Imm, CADM2 < 3)
Imm <- subset(Imm, MAP2 < 3)
Imm <- subset(Imm, GRIK2 < 3)
Imm <- subset(Imm, CLU < 3)
Imm <- subset(Imm, AQP4 < 3)
Imm <- subset(Imm, MOBP < 3)
#Imm <- subset(Imm, MBP < 3)
#Imm <- subset(Imm, CDK19 < 3)
#Imm <- subset(Imm, CD68 < 3)
#Imm <- subset(Imm, F13A1 < 3)
#Imm <- subset(Imm, STAB1 < 3)
#Imm <- subset(Imm, CD163 < 3)
#Imm <- subset(Imm, LYVE1 < 3)
#Imm <- subset(Imm, P2RY12 < 3)
#Imm <- subset(Imm, THEMIS < 3)

DimPlot(Imm)


###PCA

ImmC <- readRDS("/Ex_ChP/Ex_ChP7/CellType/Firstcelltype/ImmuneA.RDS")
ImmB <- readRDS("/Ex_ChP/Ex_ChP7/CellType/Secondcelltype_og/BAM.RDS")
ImmB$Lane <- ImmB$LANE
ImmA <- readRDS("/Ex_ChP/Ex_ChP7/CellType/Secondcelltype_og/T_cell.RDS")
ImmA$Lane <- ImmA$LANE

Imm <- merge(ImmC, y=c(ImmB, ImmA))
rm(ImmC, ImmB, ImmA)
gc()

ChP_genes <- c(known_arach, known_barrier, known_Fibro, known_Mesen, known_Mural, known_periv, known_Epithelial, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, known_dividing, known_epi_neuro_prog, known_epi_prog, known_fibro_prog, known_mglia_prog, known_neuro_prog, known_Neuron, known_neuro_Exc, known_neuro_Inh, known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B)

NFSP_fun(Imm, assay="soupXcounts")#, spikein = ChP_genes) #!!!!Using aggregated variable features from all samples instead of the default

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Imm <- RunHarmony(object = Imm, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Imm, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")  #!!!!pca or harmony
#Imm <- RunUMAP(Imm, dims=1:20, assay="soupXcounts", reduction="harmony", min.dist=1, n.neighbors=50L)
#Pick res
Idents(Imm) <- Imm@meta.data[["soupXcounts_snn_res.0.1"]]
CNDRdimfun(Imm, reduction = "harmony", res="_0.1", split.by = c("X10X_batch", "library_batch", "Subcelltype", "Majorcelltype", "Firstcelltype", "CMO_classification.global"))

Ncellfun(Imm, paste0(curr_dir, "/CellType/Imm_Cogdx"), group.by="cogdx")
#return markers
return_markersfun(Imm, assay="soupXcounts", features = Imm_soupXcounts_var_feat)
saveRDS(Imm, paste0(curr_dir, "/CellType/Firstcelltype/Imm_combined.RDS"))



curr_clust <- 3
Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_VSM), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Imm_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Imm, curr_clust)

search_markersfun(Imm, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Imm, curr_clust, known_Fibro_prog)

search_markersfun(Imm, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic) %>% unique())
search_markersfun(Imm, curr_clust, c(known_Mural, known_peric, known_VSM) %>% unique())

search_markersfun(Imm, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Imm, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Imm, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Imm, curr_clust, c(known_light) %>% unique())
search_markersfun(Imm, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Imm, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Imm, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Imm, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Imm, curr_clust, known_dividing %>% unique())

search_markersfun(Imm, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Imm, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Imm, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Imm, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Imm, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Imm, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Imm, "HSPH1", raster=F)

0: Fib_3 (                    ;; LAMA2+, SLC4A4+)
1: Fib_2 (                    ;; LAMA2+, TRPM3+)
2: Fib_3 (                    ;; LAMA2+)

###Label

Imm <- RenameIdents(Imm, "4"="T_Cell", "6"="BAM_3", "2"="BAM_2", "0"="BAM_1", "1"="BAM_1", "3"="BAM_1", "5"="BAM_1")
Imm@meta.data$Subcelltype <- Imm@active.ident
T_cell <- subset(Imm, Finalcelltype=="T_Cell")

DimPlot(Imm)


###Remove inter doublets

Idents(Imm) <- Imm@meta.data$Subcelltype
Imm <- CellSelector(DimPlot(Imm), Imm, ident = "BAM1")
BAM1 <- subset(Imm, ident="BAM1")
BAM1 <- subset(BAM1, Subcelltype=="BAM_1")

Idents(Imm) <- Imm@meta.data$Subcelltype
Imm <- CellSelector(DimPlot(Imm), Imm, ident = "BAM2")
BAM2 <- subset(Imm, ident="BAM2")
BAM2 <- subset(BAM2, Subcelltype=="BAM_2")

Idents(Imm) <- Imm@meta.data$Subcelltype
Imm <- CellSelector(DimPlot(Imm), Imm, ident = "BAM3")
BAM3 <- subset(Imm, ident="BAM3")
BAM3 <- subset(BAM3, Subcelltype=="BAM_3")


Idents(Imm) <- Imm@meta.data$Subcelltype
Imm <- CellSelector(DimPlot(Imm), Imm, ident = "T_c")
T_cell <- subset(Imm, ident="T_c")
T_cell <- subset(T_cell, Subcelltype=="T_Cell")

BAM1_meta <- BAM1@meta.data
BAM2_meta <- BAM2@meta.data
BAM3_meta <- BAM3@meta.data
T_cell_meta <- T_cell@meta.data
merged_meta <- rbind(BAM1_meta, BAM2_meta, BAM3_meta, T_cell_meta)

Imm@meta.data$Subcelltype <- NULL
Imm_meta <- Imm@meta.data

Imm_meta <- rownames_to_column(Imm_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$Subcelltype)
colnames(merged_meta) <- c("rowname", "Subcelltype")
Imm_meta <- left_join(Imm_meta, merged_meta, by="rowname")
Imm_meta <- column_to_rownames(Imm_meta, "rowname")
Imm_meta$Subcelltype <- as.character(Imm_meta$Subcelltype)
Imm_meta[is.na(Imm_meta)==TRUE] <- "Doublet"
Imm@meta.data <- Imm_meta
Imm <- subset(Imm, Subcelltype=="Doublet", invert=T)

#Reorder
Idents(Imm) <- Imm@meta.data[["Subcelltype"]] 
levels(Imm) <- c("BAM_1", "BAM_2", "BAM_3", "T_Cell")
Imm@meta.data[["Subcelltype"]] <- Idents(Imm)

Imm <- RunUMAP(object = Imm, dims = 1:20, assay = "soupXcounts", reduction =  "harmony", min.dist=0.8)
DimPlot(Imm, raster=F)

saveRDS(Imm, paste0(curr_dir, "/CellType/Finalcelltype/Imm.RDS"))


##Parenchyma

NeuralC <- readRDS("/Ex_ChP/Ex_ChP7/CellType/Firstcelltype/ParenchymaA.RDS")
NeuralB <- readRDS("/Ex_ChP/Ex_ChP4-5/CellType/Astro_first.RDS")
NeuralB$Lane <- NeuralB$LANE
NeuralA <- readRDS("/Ex_ChP/Ex_ChP7/CellType/Secondcelltype_og/Neuro.RDS")
NeuralA$Lane <- NeuralA$LANE
Neural <- readRDS("/Ex_ChP/Ex_ChP7/CellType/Secondcelltype_og/Oligo.RDS")
Neural$Lane <- Neural$LANE

Neural <- merge(NeuralC, y=c(NeuralB, NeuralA, Neural))

ChP_genes <- c(known_arach, known_barrier, known_Fibro, known_Mesen, known_Mural, known_periv, known_Epithelial, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, known_dividing, known_epi_neuro_prog, known_epi_prog, known_fibro_prog, known_mglia_prog, known_neuro_prog, known_Neuron, known_neuro_Exc, known_neuro_Inh, known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B)

NFSP_fun(Neural, assay="soupXcounts")#, spikein= ChP_genes)

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Neural <- RunHarmony(object = Neural, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Neural, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")  #!!!!pca or harmony

#Pick res
Idents(Neural) <- Neural@meta.data[["soupXcounts_snn_res.0.05"]]
CNDRdimfun(Neural, reduction = "harmony", res="_0.05", split.by = c("X10X_batch", "library_batch", "DoubletFinder_call", "CMO_classification.global", "Firstcelltype"))
#return markers
return_markersfun(Neural, assay="soupXcounts", features = Neural_soupXcounts_var_feat)
saveRDS(Neural, paste0(curr_dir, "/CellType/Firstcelltype/Neural_combined_spiked.RDS"))

Featurefun(Neural, c("KRT8", "KRT18"), "Test2",  assay="soupXcounts")



curr_clust <- 0
Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic), na.rm=TRUE)
sum(Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_VSM), na.rm=TRUE)
sum(Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Neural_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Neural, curr_clust)

search_markersfun(Neural, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Neural, curr_clust, known_fibro_prog)

search_markersfun(Neural, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic) %>% unique())
search_markersfun(Neural, curr_clust, c(known_Mural, known_peric, known_VSM) %>% unique())

search_markersfun(Neural, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Neural, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Neural, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Neural, curr_clust, c(known_light) %>% unique())
search_markersfun(Neural, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Neural, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Neural, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Neural, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Neural, curr_clust, known_dividing %>% unique())

search_markersfun(Neural, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Neural, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Neural, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Neural, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Neural, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Neural, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Neural, "TTR", raster=F)

0: Fib_3 (                    ;; LAMA2+, SLC4A4+)
1: Fib_2 (                    ;; LAMA2+, TRPM3+)
2: Fib_3 (                    ;; LAMA2+)

###Label

Neural <- RenameIdents(Neural, "0"="Oligo", "1"="Oligo", "2"="Astro", "3"="Astro", "4"="Neuron", "5"="Neuron", "6"="Neuron")
Neural@meta.data$Subcelltype <- Neural@active.ident
Neuro <- subset(Neural, Subcelltype=="Neuron")
Oligo <- subset(Neural, Subcelltype=="Oligo")
Astro <- subset(Neural, Subcelltype=="Astro")


###Venn of agreement

plts <- ggVennDiagram(list(colnames(Neural)[Neural@meta.data$DoubletFinder_call=="Doublet"], colnames(Neural)[Neural@meta.data$CMO_classification.global=="Doublet"], colnames(Neural)[Neural@meta.data$Doublet=="Doublet"]), category.names = c("Doublet Finder", "CMO", "Manual")) + theme(legend.position = "none") +labs(caption=length(colnames(Neural)))
plot(plts)
graph2png(plts, paste(curr_dir, "/CellType/Neural/Neural DoubletFinder-CMO-Manual venn.png", sep=""))



###Neuro

Neuro <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/Neuro.RDS"))
Neuro_all <- Neuro

Neuro <- subset(Neuro, CMO_classification.global=="Singlet")
Neuro <- subset(Neuro, percent.mt<5)

Featurefun(Neuro, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "VCAN", "GFAP", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="soupXcounts")
#Remove cells that don't belong
Neuro <- subset(Neuro, HTR2C < 4)
Neuro <- subset(Neuro, TTR < 4)
Neuro <- subset(Neuro, CLIC6 < 4)
#Neuro <- subset(Neuro, AQP1 < 3)
Neuro <- subset(Neuro, CFAP54 < 3)
Neuro <- subset(Neuro, FAM227A < 3)
Neuro <- subset(Neuro, ARL13B < 3)
#Neuro <- subset(Neuro, TRPM3 < 3)
#Neuro <- subset(Neuro, SLC4A4 < 3)
#Neuro <- subset(Neuro, LAMA2 < 3)
#Neuro <- subset(Neuro, LAMA4 < 3)
#Neuro <- subset(Neuro, CLDN1 < 3)
#Neuro <- subset(Neuro, CLDN11 < 3)
#Neuro <- subset(Neuro, PRICKLE1 < 3)
#Neuro <- subset(Neuro, CEMIP < 3)
#Neuro <- subset(Neuro, FOXP2 < 3)
#Neuro <- subset(Neuro, COL6A3 < 3)
#Neuro <- subset(Neuro, SLC6A13 < 3)
Neuro <- subset(Neuro, MECOM < 3)
Neuro <- subset(Neuro, VWF < 3)
Neuro <- subset(Neuro, TAGLN < 3)
Neuro <- subset(Neuro, MYH11 < 3)
Neuro <- subset(Neuro, NOTCH3 < 3)
#Neuro <- subset(Neuro, CADM2 < 3)
#Neuro <- subset(Neuro, MAP2 < 3)
#Neuro <- subset(Neuro, GRIK2 < 3)
Neuro <- subset(Neuro, CLU < 3)
Neuro <- subset(Neuro, AQP4 < 3)
Neuro <- subset(Neuro, MOBP < 3)
Neuro <- subset(Neuro, MBP < 3)
Neuro <- subset(Neuro, CDK19 < 3)
Neuro <- subset(Neuro, CD68 < 3)
Neuro <- subset(Neuro, F13A1 < 3)
Neuro <- subset(Neuro, STAB1 < 3)
Neuro <- subset(Neuro, CD163 < 3)
Neuro <- subset(Neuro, LYVE1 < 3)
Neuro <- subset(Neuro, P2RY12 < 3)
Neuro <- subset(Neuro, THEMIS < 3)


####PCA

NFSP_fun(Neuro, assay="soupXcounts", spikein=ChP_genes) #!!!!Using aggregated variable features from all samples instead of the default

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Neuro <- RunHarmony(object = Neuro, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Neuro, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")  #!!!!pca or harmony

#Pick res
Idents(Neuro) <- Neuro@meta.data[["soupXcounts_snn_res.0.03"]]
CNDRdimfun(Neuro, reduction = "harmony", res="_0.03", split.by = c("X10X_batch", "library_batch", "DoubletFinder_call", "CMO_classification.global"))
#return markers
return_markersfun(Neuro, assay="soupXcounts", features = Neuro_soupXcounts_var_feat)
saveRDS(Neuro, paste0(curr_dir, "/CellType/Secondcelltype/Neuro.RDS"))



curr_clust <- 2
Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic), na.rm=TRUE)
sum(Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_VSM), na.rm=TRUE)
sum(Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Neuro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Neuro, curr_clust)

search_markersfun(Neuro, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Neuro, curr_clust, known_Fibro_prog)

search_markersfun(Neuro, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic) %>% unique())
search_markersfun(Neuro, curr_clust, c(known_Mural, known_peric, known_VSM) %>% unique())

search_markersfun(Neuro, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Neuro, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Neuro, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Neuro, curr_clust, c(known_light) %>% unique())
search_markersfun(Neuro, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Neuro, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Neuro, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Neuro, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Neuro, curr_clust, known_dividing %>% unique())

search_markersfun(Neuro, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Neuro, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Neuro, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Neuro, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Neuro, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Neuro, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Neuro, "TTR", raster=F)


####Label

Neuro@meta.data$Subcelltype <- "Neuro"
saveRDS(Neuro, paste0(curr_dir, "/CellType/Finalcelltype/Neuro.RDS"))



###Astro

Astro <- readRDS(paste0(curr_dir, "/CellType/Secondcelltype/Astro.RDS"))
Astro_all <- Astro

Featurefun(Astro, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="soupXcounts")
#Remove cells that don't belong
Astro <- subset(Astro, HTR2C < 4)
Astro <- subset(Astro, TTR < 4)
Astro <- subset(Astro, CLIC6 < 4)
#Astro <- subset(Astro, AQP1 < 3)
#Astro <- subset(Astro, CFAP54 < 3)
#Astro <- subset(Astro, FAM227A < 3)
#Astro <- subset(Astro, ARL13B < 3)
#Astro <- subset(Astro, TRPM3 < 3)
#Astro <- subset(Astro, SLC4A4 < 3)
#Astro <- subset(Astro, LAMA2 < 3)
#Astro <- subset(Astro, LAMA4 < 3)
#Astro <- subset(Astro, CLDN1 < 3)
#Astro <- subset(Astro, CLDN11 < 3)
#Astro <- subset(Astro, PRICKLE1 < 3)
#Astro <- subset(Astro, CEMIP < 3)
#Astro <- subset(Astro, FOXP2 < 3)
#Astro <- subset(Astro, COL6A3 < 3)
#Astro <- subset(Astro, SLC6A13 < 3)
Astro <- subset(Astro, MECOM < 3)
Astro <- subset(Astro, VWF < 3)
Astro <- subset(Astro, TAGLN < 3)
Astro <- subset(Astro, MYH11 < 3)
Astro <- subset(Astro, NOTCH3 < 3)
Astro <- subset(Astro, CADM2 < 3)
Astro <- subset(Astro, MAP2 < 3)
Astro <- subset(Astro, GRIK2 < 3)
#Astro <- subset(Astro, CLU < 3)
#Astro <- subset(Astro, AQP4 < 3)
Astro <- subset(Astro, MOBP < 3)
Astro <- subset(Astro, MBP < 3)
Astro <- subset(Astro, CDK19 < 3)
Astro <- subset(Astro, CD68 < 3)
Astro <- subset(Astro, F13A1 < 3)
Astro <- subset(Astro, STAB1 < 3)
Astro <- subset(Astro, CD163 < 3)
Astro <- subset(Astro, LYVE1 < 3)
Astro <- subset(Astro, P2RY12 < 3)
Astro <- subset(Astro, THEMIS < 3)


####PCA

NFSP_fun(Astro, assay="soupXcounts") #!!!!Using aggregated variable features from all samples instead of the default

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Astro <- RunHarmony(object = Astro, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Astro, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")  #!!!!pca or harmony

#Pick res
Idents(Astro) <- Astro@meta.data[["soupXcounts_snn_res.0.05"]]
CNDRdimfun(Astro, reduction = "harmony", res="_0.05", split.by = c("X10X_batch", "library_batch", "CMO_classification.global"))
#return markers
return_markersfun(Astro, assay="soupXcounts", features = Astro_soupXcounts_var_feat)
saveRDS(Astro, paste0(curr_dir, "/CellType/Firstcelltype/Astro.RDS"))



curr_clust <- 0
Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic), na.rm=TRUE)
sum(Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_VSM), na.rm=TRUE)
sum(Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Astro_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Astro, curr_clust)

search_markersfun(Astro, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Astro, curr_clust, known_Fibro_prog)

search_markersfun(Astro, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic) %>% unique())
search_markersfun(Astro, curr_clust, c(known_Mural, known_peric, known_VSM) %>% unique())

search_markersfun(Astro, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Astro, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Astro, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Astro, curr_clust, c(known_light) %>% unique())
search_markersfun(Astro, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Astro, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Astro, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Astro, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Astro, curr_clust, known_dividing %>% unique())

search_markersfun(Astro, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Astro, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Astro, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Astro, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Astro, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Astro, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Astro, "VCAN", raster=F)

0: Astro

####Label

Astro <- RenameIdents(Astro, "2"="Astro", "0"="Doublet", "1"="Astro", "3"="Doublet", "4"="Astro", "5"="Astro", "6"="Doublet", "7"="Astro")
Astro@meta.data$Subcelltype <- Astro@active.ident
Astro <- subset(Astro, Subcelltype=="Astro")
#remove doublets
Astro <- subset(Astro, CMO_classification.global=="Singlet")
Astro <- subset(Astro, percent.mt<5)

saveRDS(Astro, paste0(curr_dir, "/CellType/Finalcelltype/Astro.RDS"))



###Oligo

Oligo <- readRDS(paste0(curr_dir, "/CellType/Firstcelltype/Oligo.RDS"))
Oligo_all <- Oligo

Featurefun(Oligo, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="soupXcounts")
#Remove cells that don't belong
Oligo <- subset(Oligo, HTR2C < 4)
Oligo <- subset(Oligo, TTR < 4)
Oligo <- subset(Oligo, CLIC6 < 4)
Oligo <- subset(Oligo, AQP1 < 3)
Oligo <- subset(Oligo, CFAP54 < 3)
Oligo <- subset(Oligo, FAM227A < 3)
Oligo <- subset(Oligo, ARL13B < 3)
Oligo <- subset(Oligo, TRPM3 < 3)
Oligo <- subset(Oligo, SLC4A4 < 3)
Oligo <- subset(Oligo, LAMA2 < 3)
Oligo <- subset(Oligo, LAMA4 < 3)
Oligo <- subset(Oligo, CLDN1 < 3)
Oligo <- subset(Oligo, CLDN11 < 3)
Oligo <- subset(Oligo, PRICKLE1 < 3)
Oligo <- subset(Oligo, CEMIP < 3)
Oligo <- subset(Oligo, FOXP2 < 3)
Oligo <- subset(Oligo, COL6A3 < 3)
Oligo <- subset(Oligo, SLC6A13 < 3)
Oligo <- subset(Oligo, MECOM < 3)
Oligo <- subset(Oligo, VWF < 3)
Oligo <- subset(Oligo, TAGLN < 3)
Oligo <- subset(Oligo, MYH11 < 3)
Oligo <- subset(Oligo, NOTCH3 < 3)
Oligo <- subset(Oligo, CADM2 < 3)
Oligo <- subset(Oligo, MAP2 < 3)
Oligo <- subset(Oligo, GRIK2 < 3)
Oligo <- subset(Oligo, CLU < 3)
Oligo <- subset(Oligo, AQP4 < 3)
#Oligo <- subset(Oligo, MOBP < 3)
#Oligo <- subset(Oligo, MBP < 3)
Oligo <- subset(Oligo, CDK19 < 3)
Oligo <- subset(Oligo, CD68 < 3)
Oligo <- subset(Oligo, F13A1 < 3)
Oligo <- subset(Oligo, STAB1 < 3)
Oligo <- subset(Oligo, CD163 < 3)
Oligo <- subset(Oligo, LYVE1 < 3)
Oligo <- subset(Oligo, P2RY12 < 3)
Oligo <- subset(Oligo, THEMIS < 3)


####PCA

NFSP_fun(Oligo, assay="soupXcounts") #!!!!Using aggregated variable features from all samples instead of the default

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Oligo <- RunHarmony(object = Oligo, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Oligo, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")  #!!!!pca or harmony

#Pick res
Idents(Oligo) <- Oligo@meta.data[["soupXcounts_snn_res.0.05"]]
CNDRdimfun(Oligo, reduction = "harmony", res="_0.05", split.by = c("X10X_batch", "library_batch", "CMO_classification.global"))
#return markers
return_markersfun(Oligo, assay="soupXcounts", features = Oligo_soupXcounts_var_feat)
saveRDS(Oligo, paste0(curr_dir, "/CellType/Firstcelltype/Oligo.RDS"))



curr_clust <- 0
Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic), na.rm=TRUE)
sum(Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_VSM), na.rm=TRUE)
sum(Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Oligo_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Oligo, curr_clust)

search_markersfun(Oligo, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Oligo, curr_clust, known_Fibro_prog)

search_markersfun(Oligo, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic) %>% unique())
search_markersfun(Oligo, curr_clust, c(known_Mural, known_peric, known_VSM) %>% unique())

search_markersfun(Oligo, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Oligo, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Oligo, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Oligo, curr_clust, c(known_light) %>% unique())
search_markersfun(Oligo, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Oligo, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Oligo, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Oligo, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Oligo, curr_clust, known_dividing %>% unique())

search_markersfun(Oligo, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Oligo, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Oligo, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Oligo, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Oligo, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Oligo, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Oligo, "TTR", raster=F)

0: Fib_3 (                    ;; LAMA2+, SLC4A4+)
1: Fib_2 (                    ;; LAMA2+, TRPM3+)
2: Fib_3 (                    ;; LAMA2+)

####Label

Oligo <- RenameIdents(Oligo, "0"="Oligo", "1"="Doublet", "2"="Oligo")
Oligo@meta.data$Subcelltype <- Oligo@active.ident
Oligo <- subset(Oligo, Subcelltype=="Oligo")
#remove doublets
Oligo <- subset(Oligo, CMO_classification.global=="Singlet")
Oligo <- subset(Oligo, percent.mt<5)

saveRDS(Oligo, paste0(curr_dir, "/CellType/Finalcelltype/Oligo.RDS"))




##Whole dataset reprocess
###load and merge

Ex_ChP_old <- readRDS("/Ex_ChP/Ex_ChP4-5/CellType/Finalcelltype/Ex_ChP_EMF.RDS")

Epi_old <- readRDS(paste(curr_dir, "/CellType/Thirdcelltype/Epi_old.RDS", sep=""))
Epi_new <- readRDS(paste(curr_dir, "/CellType/Thirdcelltype/Epi_new.RDS", sep=""))
Fibro <- readRDS(paste(curr_dir, "/CellType/Thirdcelltype/Fibro.RDS", sep=""))
Mural <- readRDS(paste(curr_dir, "/CellType/Thirdcelltype/Mural.RDS", sep=""))
Endo <- readRDS(paste(curr_dir, "/CellType/Thirdcelltype/Endo.RDS", sep=""))
Oligo <- readRDS(paste(curr_dir, "/CellType/Thirdcelltype/Oligo.RDS", sep=""))
Astro <- readRDS(paste(curr_dir, "/CellType/Thirdcelltype/Astro.RDS", sep=""))
Neuro <- readRDS(paste(curr_dir, "/CellType/Thirdcelltype/Neuro.RDS", sep=""))
Imm <- readRDS(paste(curr_dir, "/CellType/Thirdcelltype/Imm.RDS", sep=""))

Ex_ChP_merge_alt <- merge(Epi_new, c(Epi_old, Fibro, Mural, Endo, Oligo, Astro, Neuro, Imm))
rm (Epi_new, Epi_old, Fibro, Endo, Mural, Oligo, Astro, Neuro, Imm)

Ex_ChP_meta_new <- Ex_ChP_merge_alt@meta.data
Ex_ChP_meta_new <- rownames_to_column(Ex_ChP_meta_new, "rowname")
Ex_ChP_meta_new$Cell <- stringr::str_extract(Ex_ChP_meta_new$rowname, "[A-Za-z0-9]+-[0-9]")
Ex_ChP_meta_new$Lane[is.na(Ex_ChP_meta_new$Lane)] <- Ex_ChP_meta_new$LANE[is.na(Ex_ChP_meta_new$Lane)]
Ex_ChP_meta_new$Cell_Lane <- paste(Ex_ChP_meta_new$Cell, Ex_ChP_meta_new$Lane, sep="_")
rownames(Ex_ChP_meta_new) <- Ex_ChP_meta_new$rowname
Ex_ChP_merge_alt@meta.data <- Ex_ChP_meta_new


#Ex_ChP_old$NewSubcelltype <- NULL
Ex_ChP_meta_old <- Ex_ChP_old@meta.data
Ex_ChP_meta_old <- rownames_to_column(Ex_ChP_meta_old, "rowname")
Ex_ChP_meta_old$Cell <- stringr::str_extract(Ex_ChP_meta_old$rowname, "[A-Za-z0-9]+-[0-9]")
Ex_ChP_meta_old$Cell_Lane <- paste(Ex_ChP_meta_old$Cell, Ex_ChP_meta_old$Lane, sep="_")

Ex_ChP_meta_new <- data.frame(Ex_ChP_meta_new$Cell_Lane, Ex_ChP_meta_new$Subcelltype)
colnames(Ex_ChP_meta_new) <- c("Cell_Lane", "NewSubcelltype")
Ex_ChP_meta_old <- left_join(Ex_ChP_meta_old, Ex_ChP_meta_new, by="Cell_Lane")
Ex_ChP_meta_old <- Ex_ChP_meta_old[!(duplicated(Ex_ChP_meta_old$rowname)),]
rownames(Ex_ChP_meta_old) <- Ex_ChP_meta_old$rowname

Ex_ChP_meta_old$NewSubcelltype <- as.character(Ex_ChP_meta_old$NewSubcelltype)
Ex_ChP_old@meta.data <- Ex_ChP_meta_old

#Sanity check old vs new annotation
Idents(Ex_ChP_old) <- Ex_ChP_old$Subcelltype
Ncellfun(Ex_ChP_old, paste0(curr_dir, "/CellType/test.csv"), "NewSubcelltype")
#Nearly the same


Ex_ChP_new <- subset(Ex_ChP_merge_alt, c(Lane=="C39-C42_1_GEX" | Lane=="C39-C42_2_GEX" | Lane=="C39-C42_3_GEX" | Lane=="C39-C42_4_GEX" | Lane=="C43-46_1_GEX" | Lane=="C43-46_2_GEX" | Lane=="C43-46_3_GEX" | Lane=="C43-46_4_GEX" | Lane=="C47-50_1_GEX" | Lane=="C47-50_2_GEX" | Lane=="C47-50_3_GEX" | Lane=="C47-50_4_GEX"))
Ex_ChP_new@meta.data$Lane_projid <- paste(Ex_ChP_new@meta.data$projid, Ex_ChP_new@meta.data$Lane, sep="_") #I forgot this earlier
Ex_ChP_new@meta.data$sample <- Ex_ChP_new@meta.data$projid
gc()
Ex_ChP_new@meta.data[is.na(Ex_ChP_new@meta.data)==TRUE] <- "unk"

L_p <- unique(Ex_ChP_new@meta.data$Lane_projid)
L_p[order(L_p)] #Check this when making metadata for DEG later

#Cells to review
Ex_ChP_rev <- subset(Ex_ChP_merge_alt, c(Lane=="C39-C42_1_GEX" | Lane=="C39-C42_2_GEX" | Lane=="C39-C42_3_GEX" | Lane=="C39-C42_4_GEX" | Lane=="C43-46_1_GEX" | Lane=="C43-46_2_GEX" | Lane=="C43-46_3_GEX" | Lane=="C43-46_4_GEX" | Lane=="C47-50_1_GEX" | Lane=="C47-50_2_GEX" | Lane=="C47-50_3_GEX" | Lane=="C47-50_4_GEX"), invert=T)

cells_in_og <- c(Ex_ChP_old$Cell_Lane)
Ex_ChP_rev <- subset(Ex_ChP_rev, Cell_Lane==cells_in_og, invert=T)

Ex_ChP <- merge(Ex_ChP_old, c(Ex_ChP_new, Ex_ChP_rev))
saveRDS(Ex_ChP, paste0(curr_dir, "/CellType/Finalcelltype/Ex_ChP.RDS"))

rm(Ex_ChP_old, Ex_ChP_new, Ex_ChP_rev)
gc()



###process

NFSP_fun(Ex_ChP, assay="soupXcounts")
use_dims <- 1:15 #!!!Modify based on elbow plot!!!!
set.seed(1234)

Ex_ChP <- RunHarmony(object = Ex_ChP, group.by.vars = c("Lane"), dims.use = use_dims, assay="soupXcounts", max.iter.harmony = 50, plot_convergence = TRUE) #Correct for batches !!!!
CNDRfun_quick(Ex_ChP, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")
#Ex_ChP <- FindClusters(Ex_ChP, resolution = c(0.01), assay=assay)
#plts <- clustree(Ex_ChP)
#  plot(plts)
#  graph2png(plts, paste(curr_dir, "/CellType/", "Ex_ChP/Ex_ChP_cluster_tree.png", sep=""), width=14, height=6)

#Pick res
res="0.1"
Idents(Ex_ChP) <- Ex_ChP@meta.data[[paste0("soupXcounts_snn_res.", res)]]
CNDRdimfun(Ex_ChP, reduction = "harmony", res=res, split.by = c("X10X_batch", "projid", "Lane", "Subcelltype", "NewSubcelltype", "Majorcelltype"))
Ncellfun(Ex_ChP, paste0(curr_dir, "/CellType/Ex_ChP/harmony", res, "/res per cog_path.csv"), "Cog_Path")
Ncellfun(Ex_ChP, paste0(curr_dir, "/CellType/Ex_ChP/harmony", res, "/res per Lane.csv"), "Lane")
Ncellfun(Ex_ChP, paste0(curr_dir, "/CellType/Ex_ChP/harmony", res, "/res per 10X_batch.csv"), "X10X_batch")
Ncellfun(Ex_ChP, paste0(curr_dir, "/CellType/Ex_ChP/harmony", res, "/res per projid.csv"), "projid")

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
DimPlot(Ex_ChP, raster=F)
ggsave2(paste0(curr_dir, "/CellType/Ex_ChP/Finalmarkers.png", plot=plts))
#return_markersfun(Ex_ChP, assay="soupXcounts", features = Ex_ChP_soupXcounts_var_feat)

Ex_ChP <- CellCycleScoring(Ex_ChP, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes, set.ident=FALSE)
#DimPlot(Ex_ChP, group.by = "Phase", raster=F)
#DimPlot(Ex_ChP, split.by = "Phase", raster=F)
Ncellfun(Ex_ChP, paste0(curr_dir, "/CellType/Ex_ChP/CMO per phase.csv", "Phase"))
saveRDS(Ex_ChP, paste0(curr_dir, "/CellType/Finalcelltype/Ex_ChP.RDS"))


###Remove stragglers

Ex_ChP <- RunUMAP(Ex_ChP, dims=1:20, assay="soupXcounts", reduction="harmony")

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
Ex_ChP <- CellSelector(DimPlot(Ex_ChP), Ex_ChP, ident = "Epi")
Epi <- subset(Ex_ChP, ident="Epi")
Epi <- subset(Epi, Subcelltype=="Epi_1a" | Subcelltype=="Epi_1b" | Subcelltype=="Epi_2a" | Subcelltype=="Epi_2b")

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
Ex_ChP <- CellSelector(DimPlot(Ex_ChP), Ex_ChP, ident = "Fib")
Fibro <- subset(Ex_ChP, ident="Fib")
Fibro <- subset(Fibro, Subcelltype=="Fib_1" | Subcelltype=="Fib_2" | Subcelltype=="Fib_3")

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
Ex_ChP <- CellSelector(DimPlot(Ex_ChP), Ex_ChP, ident = "End")
Endo <- subset(Ex_ChP, ident="End")
Endo <- subset(Endo, Subcelltype=="Endo_1" | Subcelltype=="Endo_2" | Subcelltype=="Endo_3")

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
Ex_ChP <- CellSelector(DimPlot(Ex_ChP), Ex_ChP, ident = "Mur")
Mural <- subset(Ex_ChP, ident="Mur")
Mural <- subset(Mural, Subcelltype=="Mural_1" | Subcelltype=="Mural_2" | Subcelltype=="Mural_3")

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
Ex_ChP <- CellSelector(DimPlot(Ex_ChP), Ex_ChP, ident = "Neur")
Neuro <- subset(Ex_ChP, ident="Neur")
Neuro <- subset(Neuro, Subcelltype=="Neuro")

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
Ex_ChP <- CellSelector(DimPlot(Ex_ChP), Ex_ChP, ident = "Ast")
Astro <- subset(Ex_ChP, ident="Ast")
Astro <- subset(Astro, Subcelltype=="Astro")

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
Ex_ChP <- CellSelector(DimPlot(Ex_ChP), Ex_ChP, ident = "Oli")
Oligo <- subset(Ex_ChP, ident="Oli")
Oligo <- subset(Oligo, Subcelltype=="Oligo")

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
Ex_ChP <- CellSelector(DimPlot(Ex_ChP), Ex_ChP, ident = "Imm")
Imm <- subset(Ex_ChP, ident="Imm")
Imm <- subset(Imm, Subcelltype=="BAM_1" | Subcelltype=="BAM_2" | Subcelltype=="BAM_3" | Subcelltype=="T_Cell")

Epi_meta <- Epi@meta.data
Fibro_meta <- Fibro@meta.data
Endo_meta <- Endo@meta.data
Mural_meta <- Mural@meta.data
Neuro_meta <- Neuro@meta.data
Astro_meta <- Astro@meta.data
Oligo_meta <- Oligo@meta.data
Imm_meta <- Imm@meta.data
merged_meta <- rbind(Epi_meta, Fibro_meta, Endo_meta, Mural_meta, Neuro_meta, Astro_meta, Oligo_meta, Imm_meta)

Ex_ChP@meta.data$Subcelltype <- NULL
Ex_ChP_meta <- Ex_ChP@meta.data

Ex_ChP_meta <- rownames_to_column(Ex_ChP_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$Subcelltype)
colnames(merged_meta) <- c("rowname", "Subcelltype")
Ex_ChP_meta <- left_join(Ex_ChP_meta, merged_meta, by="rowname")
Ex_ChP_meta <- column_to_rownames(Ex_ChP_meta, "rowname")
Ex_ChP_meta$Subcelltype <- as.character(Ex_ChP_meta$Subcelltype)
Ex_ChP_meta[is.na(Ex_ChP_meta)==TRUE] <- "Doublet"
Ex_ChP@meta.data <- Ex_ChP_meta

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
DimPlot(Ex_ChP, label=TRUE, raster=F) #sanity check 

Ex_ChP <- subset(Ex_ChP, subset=Subcelltype=="Doublet", invert=TRUE)


###save and remove srt

saveRDS(Epi, paste(curr_dir, "/CellType/Finalcelltype/Epi.RDS", sep=""))
saveRDS(Fibro, paste(curr_dir, "/CellType/Finalcelltype/Fibro.RDS", sep=""))
saveRDS(Mural, paste(curr_dir, "/CellType/Finalcelltype/Mural.RDS", sep=""))
saveRDS(Endo, paste(curr_dir, "/CellType/Finalcelltype/Endo.RDS", sep=""))
saveRDS(Oligo, paste(curr_dir, "/CellType/Finalcelltype/Oligo.RDS", sep=""))
saveRDS(Astro, paste(curr_dir, "/CellType/Finalcelltype/Astro.RDS", sep=""))
saveRDS(Neuro, paste(curr_dir, "/CellType/Finalcelltype/Neuro.RDS", sep=""))
saveRDS(Imm, paste(curr_dir, "/CellType/Finalcelltype/Imm.RDS", sep=""))


Epi <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Epi.RDS", sep=""))
Fibro <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Fibro.RDS", sep=""))
Mural <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Mural.RDS", sep=""))
Endo <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Endo.RDS", sep=""))
Oligo <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Oligo.RDS", sep=""))
Astro <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Astro.RDS", sep=""))
Neuro <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Neuro.RDS", sep=""))
Imm <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Imm.RDS", sep=""))

Ex_ChP <- merge(Epi, c(Fibro, Mural, Endo, Oligo, Astro, Neuro, Imm))
rm(Epi, Fibro, Mural, Endo, Neuro, Astro, Oligo, Imm)
gc()

saveRDS(Ex_ChP, paste(curr_dir, "/CellType/Finalcelltype/Ex_ChP.RDS", sep=""))



###Ncells

Idents(Ex_ChP) <- Ex_ChP@meta.data$Subcelltype
#Cell numbers
Ncellfun(Ex_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Ex_ChP_projid_Lane.csv", sep=""), group.by="Lane_projid")
Ncellfun(Ex_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Ex_ChP_projid.csv", sep=""), group.by="projid")
Ncellfun(Ex_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Ex_ChP_cogdx.csv", sep=""), group.by="cogdx")
Ncellfun(Ex_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Ex_ChP_Cog_Path.csv", sep=""), group.by="Cog_Path")
Ncellfun(Ex_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Ex_ChP_apoe.csv", sep=""), group.by="apoe_genotype")

Idents(Ex_ChP) <- Ex_ChP@meta.data$Majorcelltype
#Cell numbers
Ncellfun(Ex_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Ex_ChP_Majorcelltype_projid.csv", sep=""), group.by="projid")
Ncellfun(Ex_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Ex_ChP_Majorcelltype_cogdx.csv", sep=""), group.by="cogdx")
Ncellfun(Ex_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Ex_ChP_Majorcelltype_Cog_Path.csv", sep=""), group.by="Cog_Path")
Ncellfun(Ex_ChP, paste(curr_dir, "/CellType/Publish/Proportions/Ex_ChP_Majorcelltype_apoe.csv", sep=""), group.by="apoe_genotype")

