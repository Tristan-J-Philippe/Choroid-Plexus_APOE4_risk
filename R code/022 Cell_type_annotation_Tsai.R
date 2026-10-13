###
### 022 Cell_type_annotation_ATAC.R
###
# Purpose: Cell type annotation and cleaning of single nuclei ATAC dataset.
# Dependencies:
source("/ChP_RNA/020 Cell_type_annotation_util.R", verbose = FALSE)
curr_dir <- "/Tsai" #output results to specific folder
set.seed(1234)
##merge
###load data

Org_ChP_ctrl_E3_06A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_06A.RDS", sep=""))
Org_ChP_ctrl_E4_06A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_06A.RDS", sep=""))
Org_ChP_ctrl_E3_15A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_15A.RDS", sep=""))
Org_ChP_ctrl_E4_13A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_13A.RDS", sep=""))
Org_ChP_ctrl_E4_15A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_15A.RDS", sep=""))
Org_ChP_ctrl_E3_36A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_36A.RDS", sep=""))
Org_ChP_ctrl_E3_36B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_36B.RDS", sep=""))
Org_ChP_ctrl_E3_36C <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_36C.RDS", sep=""))
Org_ChP_ctrl_E4_36A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_36A.RDS", sep=""))
Org_ChP_ctrl_E4_36B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_36B.RDS", sep=""))
Org_ChP_ctrl_E4_36C <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_36C.RDS", sep=""))
Org_ChP_ctrl_E3_38A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_38A.RDS", sep=""))
Org_ChP_ctrl_E3_38B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_38B.RDS", sep=""))
Org_ChP_ctrl_E3_38C <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_38C.RDS", sep=""))
Org_ChP_ctrl_E4_38A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_38A.RDS", sep=""))
Org_ChP_ctrl_E4_38B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_38B.RDS", sep=""))
Org_ChP_ctrl_E3_39A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_39A.RDS", sep=""))
Org_ChP_ctrl_E3_39B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_39B.RDS", sep=""))
Org_ChP_ctrl_E3_39C <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_39C.RDS", sep=""))
Org_ChP_ctrl_E4_39A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_39A.RDS", sep=""))
Org_ChP_ctrl_E4_39B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_39B.RDS", sep=""))
Org_ChP_ctrl_E4_39C <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_39C.RDS", sep=""))
Org_ChP_ctrl_E3_40A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_40A.RDS", sep=""))
Org_ChP_ctrl_E3_40B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_40B.RDS", sep=""))
Org_ChP_ctrl_E3_40C <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_40C.RDS", sep=""))
Org_ChP_ctrl_E4_40A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_40A.RDS", sep=""))
Org_ChP_ctrl_E4_40B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_40B.RDS", sep=""))
Org_ChP_ctrl_E4_40C <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_40C.RDS", sep=""))

Org_ChP_ctrl_E3_41A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_41A.RDS", sep=""))
Org_ChP_ctrl_E3_41B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_41B.RDS", sep=""))
Org_ChP_ctrl_E3_41C <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_41C.RDS", sep=""))
Org_ChP_ctrl_E3_44A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_44A.RDS", sep=""))
Org_ChP_ctrl_E3_44B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_44B.RDS", sep=""))
Org_ChP_ctrl_E3_44C <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_44C.RDS", sep=""))
Org_ChP_ctrl_E3_46A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_46A.RDS", sep=""))
Org_ChP_ctrl_E3_46B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_46B.RDS", sep=""))
Org_ChP_ctrl_E3_46C <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_46C.RDS", sep=""))
Org_ChP_ctrl_E3_47A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_47A.RDS", sep=""))
Org_ChP_ctrl_E3_47B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_47B.RDS", sep=""))
Org_ChP_ctrl_E3_47C <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_47C.RDS", sep=""))
Org_ChP_ctrl_E3_49A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_49A.RDS", sep=""))
Org_ChP_ctrl_E3_49B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_49B.RDS", sep=""))
Org_ChP_ctrl_E3_49C <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E3_49C.RDS", sep=""))
Org_ChP_ctrl_E4_41A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_41A.RDS", sep=""))
Org_ChP_ctrl_E4_41B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_41B.RDS", sep=""))
Org_ChP_ctrl_E4_41C <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_41C.RDS", sep=""))
Org_ChP_ctrl_E4_44A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_44A.RDS", sep=""))
Org_ChP_ctrl_E4_44B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_44B.RDS", sep=""))
Org_ChP_ctrl_E4_44C <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_44C.RDS", sep=""))
Org_ChP_ctrl_E4_46A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_46A.RDS", sep=""))
Org_ChP_ctrl_E4_46B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_46B.RDS", sep=""))
Org_ChP_ctrl_E4_46C <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_46C.RDS", sep=""))
Org_ChP_ctrl_E4_47A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_47A.RDS", sep=""))
Org_ChP_ctrl_E4_47B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_47B.RDS", sep=""))
Org_ChP_ctrl_E4_47C <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_47C.RDS", sep=""))
Org_ChP_ctrl_E4_49A <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_49A.RDS", sep=""))
Org_ChP_ctrl_E4_49B <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_49B.RDS", sep=""))
Org_ChP_ctrl_E4_49C <- readRDS(file=paste(curr_dir, "/Sample_preprocess/Org_ChP_ctrl_E4_49C.RDS", sep=""))


###merge

grep("Org_ChP_", ls(), value=TRUE) #list Org_ChP data
#merge all
#merge nuclei
Org_ChP <- merge(x=Org_ChP_ctrl_E3_06A, y=c(Org_ChP_ctrl_E4_06A, Org_ChP_ctrl_E3_15A, Org_ChP_ctrl_E4_15A, Org_ChP_ctrl_E4_13A, Org_ChP_ctrl_E3_36A, Org_ChP_ctrl_E3_36B, Org_ChP_ctrl_E3_36C, Org_ChP_ctrl_E4_36A, Org_ChP_ctrl_E4_36B, Org_ChP_ctrl_E4_36C, Org_ChP_ctrl_E3_38A, Org_ChP_ctrl_E3_38B, Org_ChP_ctrl_E3_38C, Org_ChP_ctrl_E4_38A, Org_ChP_ctrl_E4_38B, Org_ChP_ctrl_E3_39A, Org_ChP_ctrl_E3_39B, Org_ChP_ctrl_E3_39C, Org_ChP_ctrl_E4_39A, Org_ChP_ctrl_E4_39B, Org_ChP_ctrl_E4_39C, Org_ChP_ctrl_E3_40A, Org_ChP_ctrl_E3_40B, Org_ChP_ctrl_E3_40C, Org_ChP_ctrl_E4_40A, Org_ChP_ctrl_E4_40B, Org_ChP_ctrl_E4_40C))
gc()

Org_ChP <- merge(x=Org_ChP, y=c(Org_ChP_ctrl_E3_41A, Org_ChP_ctrl_E3_41B, Org_ChP_ctrl_E3_41C, Org_ChP_ctrl_E3_44A, Org_ChP_ctrl_E3_44B, Org_ChP_ctrl_E3_44C, Org_ChP_ctrl_E3_46A, Org_ChP_ctrl_E3_46B, Org_ChP_ctrl_E3_46C, Org_ChP_ctrl_E3_47A, Org_ChP_ctrl_E3_47B, Org_ChP_ctrl_E3_47C, Org_ChP_ctrl_E3_49A, Org_ChP_ctrl_E3_49B, Org_ChP_ctrl_E3_49C, Org_ChP_ctrl_E4_38A, Org_ChP_ctrl_E4_38B, Org_ChP_ctrl_E4_39A, Org_ChP_ctrl_E4_39B, Org_ChP_ctrl_E4_39C, Org_ChP_ctrl_E4_41A, Org_ChP_ctrl_E4_41B, Org_ChP_ctrl_E4_41C, Org_ChP_ctrl_E4_44A, Org_ChP_ctrl_E4_44B, Org_ChP_ctrl_E4_44C, Org_ChP_ctrl_E4_46A, Org_ChP_ctrl_E4_46B, Org_ChP_ctrl_E4_46C, Org_ChP_ctrl_E4_47A, Org_ChP_ctrl_E4_47B, Org_ChP_ctrl_E4_47C, Org_ChP_ctrl_E4_49A, Org_ChP_ctrl_E4_49B, Org_ChP_ctrl_E4_49C))
gc()

Org_ChP@meta.data$Lane <- Org_ChP@meta.data$LANE
Org_ChP@meta.data$LANE <- NULL



#Celltyping
##Main dataset
###Harmony

Org_ChP <- readRDS("/Org_ChP/Org_ChP9-5/CellType/Firstcelltype/Org_ChP.RDS")
Org_ChP <- subset(Org_ChP, CMO_classification.global=="Singlet")
gc()

NFSP_fun(Org_ChP, assay="soupXcounts")
use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
set.seed(1234)
Org_ChP <- RunHarmony(object = Org_ChP, group.by.vars = c("LANE"), dims.use = use_dims, assay="soupXcounts", max.iter.harmony = 50, plot_convergence = TRUE) #Correct for batches !!!!

CNDRfun(Org_ChP, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")
#Pick res
Idents(Org_ChP) <- Org_ChP@meta.data[["soupXcounts_snn_res.0.3"]]
CNDRdimfun(Org_ChP, reduction = "harmony", split.by = c("X10X_batch", "LANE", "sample", "CMO_classification.global", "DoubletFinder_call"), res = "0.3")

#return markers
return_markersfun(Org_ChP, assay="soupXcounts", features = Org_ChP_soupXcounts_var_feat)



Ncellfun(Org_ChP, "CellType/Org_ChP/cluster per sample.csv", "sample")
#Ncellfun(Org_ChP, "cluster per day.csv", "time")


###Dividing cells spend more time in S phase

Org_ChP <- CellCycleScoring(Org_ChP, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes, set.ident=FALSE)
DimPlot(Org_ChP, group.by = "Phase", raster = FALSE)
Ncellfun(Org_ChP, "CMO per phase.csv", "Phase")


###Resolution 0.3

curr_clust <- 6
Org_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Org_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% unique(c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach)), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_VSM), na.rm=TRUE)
sum(Org_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% unique(c(known_Epithelial, known_ChPMature, known_dark, other_dark, known_myoepi, known_light, known_ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3)), na.rm=TRUE)
sum(Org_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Org_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% unique(c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog)), na.rm=TRUE)
sum(Org_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% unique(c(known_epi_neuro_prog, known_neuro_prog)), na.rm=TRUE)
sum(Org_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% unique(c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS)), na.rm=TRUE)
sum(Org_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% unique(c(known_astro, known_oligo, known_OPC)), na.rm=TRUE)
sum(Org_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% unique(c(Cell_Stress, whole_cells)), na.rm=TRUE)
sum(Org_ChP_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% unique(c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B)), na.rm=TRUE)

#Top20_Featuresfun(Org_ChP, 1)

#FIBRO
search_markersfun(Org_ChP, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_arach)%>% unique())
search_markersfun(Org_ChP, curr_clust, known_periv)
search_markersfun(Org_ChP, curr_clust, known_menin)
search_markersfun(Org_ChP, curr_clust, known_fibro_prog)

#VASC
search_markersfun(Epi, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic) %>% unique())
search_markersfun(Epi, curr_clust, c(known_Mural, known_peric, known_VSM) %>% unique())

#EPI
search_markersfun(Org_ChP, curr_clust, c(known_dark, other_dark, Whyss_Correy_Epithelial_highMT) %>% unique())
search_markersfun(Org_ChP, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Org_ChP, curr_clust, c(known_light) %>% unique())
search_markersfun(Org_ChP, curr_clust, c(known_ependymal) %>% unique())
search_markersfun(Org_ChP, curr_clust, c(known_Epithelial))

#PROG
search_markersfun(Org_ChP, curr_clust, known_dividing)
search_markersfun(Org_ChP, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())
search_markersfun(Org_ChP, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

#NEURO
search_markersfun(Org_ChP, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())

#IMMUNE
search_markersfun(Org_ChP, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

#GLIA
search_markersfun(Org_ChP, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())

#low QC
search_markersfun(Org_ChP, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Org_ChP, "KRT17")

Featurefun(Org_ChP, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "MT-CO2", "RGN", "SI", "CFTR", "TRPM3", "KCNMA1", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "PECAM1", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS", "MT-CO2", "RSPO2", "RSPO3", "MKI67", "TUJ1"), "Test2",  assay="soupXcounts")



###First pass Isolate #4-5 Epis and keep high MT here, remove high MT everywhere else and rerun above

MT <- subset(Org_ChP, soupXcounts_snn_res.0.3=="4" | soupXcounts_snn_res.0.3=="5")
MT$temp <- "keep"
MT_less <- subset(Org_ChP, soupXcounts_snn_res.0.3=="4" | soupXcounts_snn_res.0.3=="5", invert=T)
MT_less <- subset(MT_less, percent.mt<5)
MT_less$temp <- "keep"

MT_meta <- MT@meta.data
MT_less_meta <- MT_less@meta.data
merged_meta <- rbind(MT_meta, MT_less_meta)

Org_ChP@meta.data$temp <- NULL
Org_ChP_meta <- Org_ChP@meta.data

Org_ChP_meta <- rownames_to_column(Org_ChP_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$temp)
colnames(merged_meta) <- c("rowname", "temp")
Org_ChP_meta <- left_join(Org_ChP_meta, merged_meta, by="rowname")
Org_ChP_meta <- column_to_rownames(Org_ChP_meta, "rowname")
Org_ChP_meta$temp <- as.character(Org_ChP_meta$temp)
Org_ChP_meta[is.na(Org_ChP_meta)==TRUE] <- "Doublet"

Org_ChP@meta.data <- Org_ChP_meta
Org_ChP <- subset(Org_ChP, subset=temp=="Doublet", invert=TRUE)
rm(MT, MT_less)
gc()
Org_ChP$temp <- NULL


###Second pass
0: Prog_Epi_Neuro, but some Epi
1: Neuro
2: Fib_3
3: Fib_3 (part)
4: Epi_1 (some prog)
5: Epi_4 (or intestinal)
6: Fib_2
7: Epi_23
8: Prog_Neuro
9: Prog_Fib
10: Dividing
11: BAM
12: Fib_1
13: Prog_Epi
14: Doublet



Featurefun(Org_ChP, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "MT-CO2", "RGN", "SI", "CFTR", "TRPM3", "KCNMA1", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "PECAM1", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS", "RSPO2", "RSPO3", "MKI67", "TUJ1"), "Test2",  assay="soupXcounts")


###Rename all identities

#Idents(Org_ChP) <- Org_ChP@meta.data[["soupXcounts_snn_res.0.3"]]
Org_ChP <- RenameIdents(Org_ChP, '0'='Prog_Epi_Neuro', '1'='Neuron', '2'='Fib_3', '3'='Fib_3', '4'='Epi_1', '5'='Epi_4', '6'='Fib_2', '7'='Epi_23', '8'='Prog_Neuro', '9'='Prog_Fib', '10'='Dividing', '11'='BAM', '12'='Fib_1', '13'='Prog_Epi', '14'='Doublet')
Org_ChP@meta.data$Firstcelltype <- Idents(Org_ChP)

Org_ChP <- RenameIdents(Org_ChP, 'Prog_Epi_Neuro'='Epi', 'Fib_1'='Fib', 'Fib_3'='Fib', 'Epi_1'='Epi', 'Epi_4'='Epi', 'Fib_3'='Fib', 'Epi_23'='Epi', 'Prog_Neuro'='Neuron', 'Prog_Fib'='Prog', 'Fib_2'='Fib', 'Prog_Epi'='Epi')
Org_ChP@meta.data$Majorcelltype <- Idents(Org_ChP)
saveRDS(Org_ChP, paste(curr_dir, "/CellType/Firstcelltype/Org_ChP.RDS", sep=""))



##Subtypes

Org_ChP <- readRDS(paste(curr_dir, "/CellType/Firstcelltype/Org_ChP.RDS", sep=""))
for(i in unique(Org_ChP@meta.data$Majorcelltype)){
  sub <- subset(Org_ChP, Majorcelltype==i)
  saveRDS(sub, paste(curr_dir, "/CellType/Firstcelltype/", i, ".RDS", sep=""))
}



###Epi
####Remove doublets

Epi <- readRDS(paste(curr_dir, "/CellType/Firstcelltype/Epi.RDS", sep=""))
#Epi <- subset(Org_ChP, Majorcelltype=="Epi")

Featurefun(Epi, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS", "MT-CO2"), "Test2",  assay="soupXcounts")

#"Second" pass
Idents(Epi) <- Epi@meta.data[["soupXcounts_snn_res.0.3"]]
Epi <- RenameIdents(Epi, "0"="Prog_Epi_Neuro", "5"="Prog_Epi_Neuro", "2"="Epi_2", "3"="Epi_2", "4"="Epi_3", "6"="Epi_3", "1"="Epi_4", "8"="Epi_4", "7"="Epi_Prog")
Epi@meta.data$Finalcelltype <- Idents(Epi)
Prog_EN <- subset(Epi, Finalcelltype=="Prog_Epi_Neuro")
saveRDS(Prog_EN, paste(curr_dir, "/CellType/Secondcelltype/Prog_Epi_Neuro.RDS", sep=""))

Prog_E <- subset(Epi, Finalcelltype=="Prog_Epi")
saveRDS(Prog_E, paste(curr_dir, "/CellType/Secondcelltype/Prog_Epi.RDS", sep=""))

Epi <- subset(Epi, Finalcelltype=="Prog_Epi_Neuro", invert=T)
Epi <- subset(Epi, Finalcelltype=="Prog_Epi", invert=T)

#Remove cells that don't belong
#Epi <- subset(Epi, HTR2C < 3)
#Epi <- subset(Epi, TTR < 3)
#Epi <- subset(Epi, CLIC6 < 3)
#Epi <- subset(Epi, AQP1 < 3)
#Epi <- subset(Epi, CFAP54 < 3)
#Epi <- subset(Epi, FAM227A < 3)
#Epi <- subset(Epi, ARL13B < 3)
#Epi <- subset(Epi, TRPM3 < 3)
Epi <- subset(Epi, SLC4A4 < 3)
Epi <- subset(Epi, LAMA2 < 3)
Epi <- subset(Epi, LAMA4 < 3)
Epi <- subset(Epi, CLDN1 < 3)
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


####PCA

NFSP_fun(Epi, assay="soupXcounts") #, spikein=known_dark) #only for cleanup

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Epi <- RunHarmony(object = Epi, group.by.vars = c("LANE"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Epi, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")

#Pick res
Idents(Epi) <- Epi@meta.data[["soupXcounts_snn_res.0.3"]]
CNDRdimfun(Epi, reduction = "harmony", res="_0.3", split.by = c("X10X_batch", "library_batch", "DoubletFinder_call", "CMO_classification.global"))
#return markers
return_markersfun(Epi, assay="soupXcounts", features = Epi_soupXcounts_var_feat)


####Dividing cells spend more time in S phase

DimPlot(Epi, group.by = "Phase", raster = FALSE)
Ncellfun(Epi, "CMO per phase.csv", "Phase")

Idents(Epi) <- Epi@meta.data$Firstcelltype
DimPlot(Epi)



curr_clust <- 3
Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic), na.rm=TRUE)
sum(Epi_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_VSM), na.rm=TRUE)
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
FeaturePlot(Epi, c("MT-CO2", "HTR2C", "TTR", "CFAP54"), raster=F)


####Label

#Third pass
Idents(Epi) <- Epi@meta.data[["soupXcounts_snn_res.0.3"]]
Epi <- RenameIdents(Epi, "1"="Epi_1", "2"="Epi_1", "3"="Epi_3", "4"="Epi_3", "5"="Epi_2", "6"="Epi_2", "7"="Epi_2", "8"="Dividing", "0"="Epi_4")
Epi@meta.data$Finalcelltype <- Epi@active.ident

Dividing_E <- subset(Epi, Finalcelltype=="Dividing")
saveRDS(Epi, paste(curr_dir, "/CellType/Secondcelltype/Dividing_E.RDS", sep=""))
Epi <- subset(Epi, Finalcelltype=="Dividing", invert=T)

#Fourth pass
Idents(Epi) <- Epi@meta.data[["soupXcounts_snn_res.0.3"]]
Epi <- RenameIdents(Epi, "1"="Epi_1", "2"="Epi_1", "3"="Epi_3", "4"="Epi_3", "5"="Epi_2", "6"="Epi_2", "7"="Epi_2", "8"="Epi_2", "0"="Epi_4")
Epi@meta.data$Finalcelltype <- Epi@active.ident
DimPlot(Epi)


####Remove interdoublets

Idents(Epi) <- Epi@meta.data$Finalcelltype
Epi <- CellSelector(DimPlot(Epi), Epi, ident = "Epi1")
Epi_1 <- subset(Epi, ident="Epi1")
Epi_1 <- subset(Epi_1, Finalcelltype=="Epi_1")

Idents(Epi) <- Epi@meta.data$Finalcelltype
Epi <- CellSelector(DimPlot(Epi), Epi, ident = "Epi2")
Epi_2 <- subset(Epi, ident="Epi2")
Epi_2 <- subset(Epi_2, Finalcelltype=="Epi_2")

Idents(Epi) <- Epi@meta.data$Finalcelltype
Epi <- CellSelector(DimPlot(Epi), Epi, ident = "Epi3")
Epi_3 <- subset(Epi, ident="Epi3")
Epi_3 <- subset(Epi_3, Finalcelltype=="Epi_3")

Idents(Epi) <- Epi@meta.data$Finalcelltype
Epi <- CellSelector(DimPlot(Epi), Epi, ident = "Epi4")
Epi_4 <- subset(Epi, ident="Epi4")
Epi_4 <- subset(Epi_4, Finalcelltype=="Epi_4")

Epi_1_meta <- Epi_1@meta.data
Epi_2_meta <- Epi_2@meta.data
Epi_3_meta <- Epi_3@meta.data
Epi_4_meta <- Epi_4@meta.data
merged_meta <- rbind(Epi_1_meta, Epi_2_meta, Epi_3_meta, Epi_4_meta)

Epi@meta.data$Finalcelltype <- NULL
Epi_meta <- Epi@meta.data

Epi_meta <- rownames_to_column(Epi_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$Finalcelltype)
colnames(merged_meta) <- c("rowname", "Finalcelltype")
Epi_meta <- left_join(Epi_meta, merged_meta, by="rowname")
Epi_meta <- column_to_rownames(Epi_meta, "rowname")
Epi_meta$Finalcelltype <- as.character(Epi_meta$Finalcelltype)
Epi_meta[is.na(Epi_meta)==TRUE] <- "Doublet"
Epi@meta.data <- Epi_meta
Epi <- subset(Epi, Finalcelltype=="Doublet", invert=T)

#Reorder
Idents(Epi) <- Epi@meta.data[["Finalcelltype"]] 
levels(Epi) <- c("Epi_1", "Epi_2", "Epi_3", "Epi_4")
Epi@meta.data[["Finalcelltype"]] <- Idents(Epi)
Epi <- RunUMAP(object = Epi, dims = 1:20, assay = "soupXcounts", reduction =  "harmony", min.dist=1.4)
DimPlot(Epi, raster=F)

saveRDS(Epi, paste(curr_dir, "/CellType/Secondcelltype/Epi.RDS", sep=""))



##Fib

Fib <- readRDS(paste(curr_dir, "/CellType/Firstcelltype/Fib.RDS", sep=""))

Featurefun(Fib, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="soupXcounts")


#Remove cells that don't belong
Fib <- subset(Fib, HTR2C < 3)
Fib <- subset(Fib, TTR < 3)
Fib <- subset(Fib, CLIC6 < 3)
Fib <- subset(Fib, AQP1 < 3)
Fib <- subset(Fib, CFAP54 < 3)
Fib <- subset(Fib, FAM227A < 3)
Fib <- subset(Fib, ARL13B < 3)
#Fib <- subset(Fib, TRPM3 < 3)
#Fib <- subset(Fib, SLC4A4 < 3)
#Fib <- subset(Fib, LAMA2 < 3)
#Fib <- subset(Fib, LAMA4 < 3)
#Fib <- subset(Fib, CLDN1 < 3)
#Fib <- subset(Fib, CLDN11 < 3)
#Fib <- subset(Fib, PRICKLE1 < 3)
#Fib <- subset(Fib, CEMIP < 3)
#Fib <- subset(Fib, FOXP2 < 3)
#Fib <- subset(Fib, COL6A3 < 3)
#Fib <- subset(Fib, SLC6A13 < 3)
Fib <- subset(Fib, MECOM < 3)
Fib <- subset(Fib, VWF < 3)
Fib <- subset(Fib, TAGLN < 3)
Fib <- subset(Fib, MYH11 < 3)
Fib <- subset(Fib, NOTCH3 < 3)
Fib <- subset(Fib, CADM2 < 3)
Fib <- subset(Fib, MAP2 < 3)
Fib <- subset(Fib, GRIK2 < 3)
Fib <- subset(Fib, CLU < 3)
Fib <- subset(Fib, AQP4 < 3)
Fib <- subset(Fib, MOBP < 3)
Fib <- subset(Fib, MBP < 3)
Fib <- subset(Fib, CDK19 < 3)
Fib <- subset(Fib, CD68 < 3)
Fib <- subset(Fib, F13A1 < 3)
Fib <- subset(Fib, STAB1 < 3)
Fib <- subset(Fib, CD163 < 3)
Fib <- subset(Fib, LYVE1 < 3)
Fib <- subset(Fib, P2RY12 < 3)
Fib <- subset(Fib, THEMIS < 3)


###PCA

NFSP_fun(Fib, assay="soupXcounts")

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Fib <- RunHarmony(object = Fib, group.by.vars = c("LANE"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Fib, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")

#Pick res
Idents(Fib) <- Fib@meta.data[["soupXcounts_snn_res.0.1"]]
CNDRdimfun(Fib, reduction = "harmony", res="_0.1", split.by = c("X10X_batch", "library_batch", "DoubletFinder_call", "CMO_classification.global"))
#return markers
return_markersfun(Fib, assay="soupXcounts", features = Fib_soupXcounts_var_feat)

#saveRDS(Fib, "/Org_ChP/Org_ChP4/CellType/Firstcelltype/Fib.RDS")


####Dividing cells spend more time in S phase

Idents(Fib) <- Fib@meta.data$Firstcelltype
DimPlot(Fib)

DimPlot(Fib, group.by = "Phase", raster = FALSE)
Ncellfun(Fib, "CMO per phase.csv", "Phase")



curr_clust <- 1
Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic), na.rm=TRUE)
sum(Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_VSM), na.rm=TRUE)
sum(Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Fib_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Fib, curr_clust)

search_markersfun(Fib, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Fib, curr_clust, known_fibro_prog)

search_markersfun(Fib, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic) %>% unique())
search_markersfun(Fib, curr_clust, c(known_Mural, known_peric, known_VSM) %>% unique())

search_markersfun(Fib, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Fib, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Fib, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Fib, curr_clust, c(known_light) %>% unique())
search_markersfun(Fib, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Fib, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Fib, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Fib, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Fib, curr_clust, known_dividing %>% unique())

search_markersfun(Fib, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Fib, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Fib, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Fib, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Fib, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Fib, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Fib, "TTR", raster=F)

0 Fib_4
1 Fib_2
2 Fib_2
3 Fib_2
4 Fib_3

###Label

#Second pass
Fib <- RenameIdents(Fib, "0"="Fib_3", "1"="Fib_3", "2"="Fib_3", "3"="Fib_3", "4"="Fib_2", "5"="Fib_3", "6"="Fib_1", "7"="Fib_1")
Fib@meta.data$Finalcelltype <- Fib@active.ident

DimPlot(Fib)


###Remove interdoublets

Idents(Fib) <- Fib@meta.data$Finalcelltype
Fib <- CellSelector(DimPlot(Fib), Fib, ident = "Fib1")
Fib_1 <- subset(Fib, ident="Fib1")
Fib_1 <- subset(Fib_1, Finalcelltype=="Fib_1")

Idents(Fib) <- Fib@meta.data$Finalcelltype
Fib <- CellSelector(DimPlot(Fib), Fib, ident = "Fib2")
Fib_2 <- subset(Fib, ident="Fib2")
Fib_2 <- subset(Fib_2, Finalcelltype=="Fib_2")

Idents(Fib) <- Fib@meta.data$Finalcelltype
Fib <- CellSelector(DimPlot(Fib), Fib, ident = "Fib3")
Fib_3 <- subset(Fib, ident="Fib3")
Fib_3 <- subset(Fib_3, Finalcelltype=="Fib_3")

Fib_1_meta <- Fib_1@meta.data
Fib_2_meta <- Fib_2@meta.data
Fib_3_meta <- Fib_3@meta.data
merged_meta <- rbind(Fib_1_meta, Fib_2_meta, Fib_3_meta)

Fib@meta.data$Finalcelltype <- NULL
Fib_meta <- Fib@meta.data

Fib_meta <- rownames_to_column(Fib_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$Finalcelltype)
colnames(merged_meta) <- c("rowname", "Finalcelltype")
Fib_meta <- left_join(Fib_meta, merged_meta, by="rowname")
Fib_meta <- column_to_rownames(Fib_meta, "rowname")
Fib_meta$Finalcelltype <- as.character(Fib_meta$Finalcelltype)
Fib_meta[is.na(Fib_meta)==TRUE] <- "Doublet"
Fib@meta.data <- Fib_meta
Fib <- subset(Fib, Finalcelltype=="Doublet", invert=T)

#Reorder
Idents(Fib) <- Fib@meta.data[["Finalcelltype"]] 
levels(Fib) <- c("Fib_1", "Fib_2", "Fib_3")
Fib@meta.data[["Finalcelltype"]] <- Idents(Fib)
Fib <- RunUMAP(object = Fib, dims = 1:20, assay = "soupXcounts", reduction =  "harmony", min.dist=0.8)
DimPlot(Fib, raster=F)

saveRDS(Fib, paste(curr_dir, "/CellType/Secondcelltype/Fib.RDS", sep=""))




##Neuro

Neuro <- readRDS(paste(curr_dir, "/CellType/Firstcelltype/Neuron.RDS", sep=""))

Featurefun(Neuro, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="soupXcounts")

#Remove cells that don't belong
Neuro <- subset(Neuro, HTR2C < 3)
Neuro <- subset(Neuro, TTR < 3)
Neuro <- subset(Neuro, CLIC6 < 3)
#Neuro <- subset(Neuro, AQP1 < 3)
Neuro <- subset(Neuro, CFAP54 < 3)
Neuro <- subset(Neuro, FAM227A < 3)
Neuro <- subset(Neuro, ARL13B < 3)
Neuro <- subset(Neuro, TRPM3 < 3)
#Neuro <- subset(Neuro, SLC4A4 < 3)
Neuro <- subset(Neuro, LAMA2 < 3)
Neuro <- subset(Neuro, LAMA4 < 3)
Neuro <- subset(Neuro, CLDN1 < 3)
Neuro <- subset(Neuro, CLDN11 < 3)
Neuro <- subset(Neuro, PRICKLE1 < 3)
Neuro <- subset(Neuro, CEMIP < 3)
Neuro <- subset(Neuro, FOXP2 < 3)
Neuro <- subset(Neuro, COL6A3 < 3)
#Neuro <- subset(Neuro, SLC6A13 < 3)
Neuro <- subset(Neuro, MECOM < 3)
Neuro <- subset(Neuro, VWF < 3)
Neuro <- subset(Neuro, TAGLN < 3)
Neuro <- subset(Neuro, MYH11 < 3)
Neuro <- subset(Neuro, NOTCH3 < 3)
#Neuro <- subset(Neuro, CADM2 < 3)
#Neuro <- subset(Neuro, MAP2 < 3)
#Neuro <- subset(Neuro, GRIK2 < 3)
#Neuro <- subset(Neuro, CLU < 3)
#Neuro <- subset(Neuro, AQP4 < 3)
#Neuro <- subset(Neuro, MOBP < 3)
#Neuro <- subset(Neuro, MBP < 3)
#Neuro <- subset(Neuro, CDK19 < 3)
Neuro <- subset(Neuro, CD68 < 3)
Neuro <- subset(Neuro, F13A1 < 3)
Neuro <- subset(Neuro, STAB1 < 3)
Neuro <- subset(Neuro, CD163 < 3)
Neuro <- subset(Neuro, LYVE1 < 3)
Neuro <- subset(Neuro, P2RY12 < 3)
Neuro <- subset(Neuro, THEMIS < 3)


###PCA

NFSP_fun(Neuro, assay="soupXcounts")

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Neuro <- RunHarmony(object = Neuro, group.by.vars = c("LANE"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Neuro, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")

#Pick res
Idents(Neuro) <- Neuro@meta.data[["soupXcounts_snn_res.0.03"]]
CNDRdimfun(Neuro, reduction = "harmony", res="_0.03", split.by = c("X10X_batch", "library_batch", "DoubletFinder_call", "CMO_classification.global"))
#return markers
return_markersfun(Neuro, assay="soupXcounts", features = Neuro_var_feat)

#saveRDS(Neuro, "/Org_ChP/Org_ChP4/CellType/Firstcelltype/Neuro.RDS")


###Dividing cells spend more time in S phase

DimPlot(Neuro, group.by = "Phase", raster = FALSE)
Ncellfun(Neuro, "CMO per phase.csv", "Phase")



curr_clust <- 1
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
search_markersfun(Neuro, curr_clust, known_fibro_prog)

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
FeaturePlot(Neuro, "MAP2", raster=F)

1 Doublet
2 Doublet


###Label

Neuro@meta.data$Finalcelltype <- Neuro@meta.data$Firstcelltype

saveRDS(Neuro, paste(curr_dir, "/CellType/Secondcelltype/Neuro.RDS", sep=""))



##BAM

Imm <- readRDS(paste(curr_dir, "/CellType/Firstcelltype/Imm.RDS", sep=""))


Featurefun(Imm, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="soupXcounts")

#"Second" pass
Idents(Imm) <- Imm@meta.data[["soupXcounts_snn_res.0.05"]]
Imm <- RenameIdents(Imm, "0"="BAM_1", "1"="BAM_2", "2"="BAM_1", "3"="BAM_2", "4"="Prog_Vasc", "5"="BAM_3")
Imm@meta.data$Finalcelltype <- Idents(Imm)
Vasc <- subset(Imm, Finalcelltype=="Prog_Vasc")
saveRDS(Vasc, paste(curr_dir, "/CellType/Finalcelltype/Vasc.RDS", sep=""))

Imm <- subset(Imm, Finalcelltype=="Prog_Vasc", invert=T)


###Remove cells that express genes that don't belong

Imm <- subset(Imm, HTR2C < 4)
Imm <- subset(Imm, TTR < 4)
Imm <- subset(Imm, CLIC6 < 4)
Imm <- subset(Imm, AQP1 < 3)
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
Imm <- subset(Imm, MBP < 3)
Imm <- subset(Imm, CDK19 < 3)
#Imm <- subset(Imm, CD68 < 3)
#Imm <- subset(Imm, F13A1 < 3)
#Imm <- subset(Imm, STAB1 < 3)
#Imm <- subset(Imm, CD163 < 3)
#Imm <- subset(Imm, LYVE1 < 3)
#Imm <- subset(Imm, P2RY12 < 3)
#Imm <- subset(Imm, THEMIS < 3)


###PCA

NFSP_fun(Imm, assay="soupXcounts")

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Imm <- RunHarmony(object = Imm, group.by.vars = c("LANE"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Imm, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")

#Pick res
Idents(Imm) <- Imm@meta.data[["soupXcounts_snn_res.0.09"]]
CNDRdimfun(Imm, reduction = "harmony", res="_0.09", split.by = c("X10X_batch", "library_batch", "DoubletFinder_call", "CMO_classification.global"))
#return markers
return_markersfun(Imm, assay="soupXcounts", features = Imm_soupXcounts_var_feat)


####Dividing cells spend more time in S phase

Idents(Imm) <- Imm@meta.data$Firstcelltype
DimPlot(Imm)

DimPlot(Imm, group.by = "Phase", raster = FALSE)
Ncellfun(Imm, "CMO per phase.csv", "Phase")



curr_clust <- 0
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
search_markersfun(Imm, curr_clust, known_fibro_prog)

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
FeaturePlot(Imm, "TTR", raster=F)

2? MECOM??? microglia? or macro prog
3?
  4 doublet?
  
  ###Label
  
#Third pass
Idents(Imm) <- Imm@meta.data[["soupXcounts_snn_res.0.09"]]
Imm <- RenameIdents(Imm, "0"="BAM_1", "1"="BAM_2", "2"="BAM_1", "3"="BAM_2", "4"="BAM_3")
Imm@meta.data$Finalcelltype <- Idents(Imm)
Imm <- RunUMAP(object = Imm, dims = 1:20, assay = "soupXcounts", reduction =  "harmony", min.dist=2)
DimPlot(Imm)


###Remove interdoublets

Idents(Imm) <- Imm@meta.data$Finalcelltype
Imm <- CellSelector(DimPlot(Imm), Imm, ident = "BAM1")
BAM_1 <- subset(Imm, ident="BAM1")
BAM_1 <- subset(BAM_1, Finalcelltype=="BAM_1")

Idents(Imm) <- Imm@meta.data$Finalcelltype
Imm <- CellSelector(DimPlot(Imm), Imm, ident = "BAM2")
BAM_2 <- subset(Imm, ident="BAM2")
BAM_2 <- subset(BAM_2, Finalcelltype=="BAM_2")

Idents(Imm) <- Imm@meta.data$Finalcelltype
Imm <- CellSelector(DimPlot(Imm), Imm, ident = "BAM3")
BAM_3 <- subset(Imm, ident="BAM3")
BAM_3 <- subset(BAM_3, Finalcelltype=="BAM_3")

BAM_1_meta <- BAM_1@meta.data
BAM_2_meta <- BAM_2@meta.data
BAM_3_meta <- BAM_3@meta.data
merged_meta <- rbind(BAM_1_meta, BAM_2_meta, BAM_3_meta)

Imm@meta.data$Finalcelltype <- NULL
Imm_meta <- Imm@meta.data

Imm_meta <- rownames_to_column(Imm_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$Finalcelltype)
colnames(merged_meta) <- c("rowname", "Finalcelltype")
Imm_meta <- left_join(Imm_meta, merged_meta, by="rowname")
Imm_meta <- column_to_rownames(Imm_meta, "rowname")
Imm_meta$Finalcelltype <- as.character(Imm_meta$Finalcelltype)
Imm_meta[is.na(Imm_meta)==TRUE] <- "Doublet"
Imm@meta.data <- Imm_meta
Imm <- subset(Imm, Finalcelltype=="Doublet", invert=T)

#Reorder
Idents(Imm) <- Imm@meta.data[["Finalcelltype"]] 
levels(Imm) <- c("Imm_1", "Imm_2", "Imm_3")
Imm@meta.data[["Finalcelltype"]] <- Idents(Imm)
Imm <- RunUMAP(object = Imm, dims = 1:20, assay = "soupXcounts", reduction =  "harmony", min.dist=2.5)
DimPlot(Imm, raster=F)

saveRDS(Imm, paste(curr_dir, "/CellType/Secondcelltype/Imm.RDS", sep=""))


##Dividing

Dividing <- readRDS(paste(curr_dir, "/CellType/Firstcelltype/Dividing.RDS", sep=""))
Dividing_E <- readRDS(paste(curr_dir, "/CellType/Firstcelltype/Dividing_E.RDS", sep=""))
Dividing <- merge(Dividing, Dividing_E)

Featurefun(Dividing, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="soupXcounts")


###PCA

NFSP_fun(Dividing, assay="soupXcounts")

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Dividing <- RunHarmony(object = Dividing, group.by.vars = c("LANE"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Dividing, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")

#Pick res
Idents(Dividing) <- Dividing@meta.data[["soupXcounts_snn_res.0.07"]]
CNDRdimfun(Dividing, reduction = "harmony", res="_0.07", split.by = c("X10X_batch", "library_batch", "DoubletFinder_call", "CMO_classification.global"))
#return markers
return_markersfun(Dividing, assay="soupXcounts", features = Dividing_soupXcounts_var_feat)

#saveRDS(Dividing, "/Org_ChP/Org_ChP4/CellType/Firstcelltype/Dividing.RDS")


###Dividing cells spend more time in S phase

DimPlot(Dividing, group.by = "Phase", raster = FALSE)
Ncellfun(Dividing, "CMO per phase.csv", "Phase")



curr_clust <- 6
Dividing_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Dividing_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Dividing_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic), na.rm=TRUE)
sum(Dividing_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_VSM), na.rm=TRUE)
sum(Dividing_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Dividing_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Dividing_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Dividing_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Dividing_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Dividing_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Dividing_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Dividing_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Dividing, curr_clust)

search_markersfun(Dividing, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Dividing, curr_clust, known_fibro_prog)

search_markersfun(Dividing, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic) %>% unique())
search_markersfun(Dividing, curr_clust, c(known_Mural, known_peric, known_VSM) %>% unique())

search_markersfun(Dividing, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Dividing, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Dividing, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Dividing, curr_clust, c(known_light) %>% unique())
search_markersfun(Dividing, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Dividing, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Dividing, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Dividing, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Dividing, curr_clust, known_dividing %>% unique())

search_markersfun(Dividing, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Dividing, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Dividing, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Dividing, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Dividing, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Dividing, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Dividing, "TTR", raster=F)

0-5 Epi_Neuro
1,2,4 Fib
3 Epi
6 Neuro

###Label

#First pass
Dividing <- RenameIdents(Dividing, "0"="Dividing_Neuro_Epi", "1"="Dividing_Fib", "2"="Dividing_Fib", "3"="Dividing_Epi", "4"="Dividing_Fib", "5"="Dividing_Neuro_Epi", "6"="Dividing_Neuro")
Dividing@meta.data$Finalcelltype <- Dividing@active.ident
Dividing <- RunUMAP(object = Dividing, dims = 1:20, assay = "soupXcounts", reduction =  "harmony", min.dist=1.5)
DimPlot(Dividing)
saveRDS(Dividing, paste(curr_dir, "/CellType/Secondcelltype/Dividing.RDS", sep=""))



##Prog

Prog_EN <- readRDS(paste(curr_dir, "/CellType/Secondcelltype/Prog_Epi_Neuro.RDS", sep=""))
Prog_E <- readRDS(paste(curr_dir, "/CellType/Secondcelltype/Prog_Epi.RDS", sep=""))
Prog_F <- readRDS(paste(curr_dir, "/CellType/Secondcelltype/Prog_Fib.RDS", sep=""))

Featurefun(Prog, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="soupXcounts")


###PCA

NFSP_fun(Prog, assay="soupXcounts")

use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
Prog <- RunHarmony(object = Prog, group.by.vars = c("LANE"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Prog, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")

#Pick res
Idents(Prog) <- Prog@meta.data[["soupXcounts_snn_res.0.3"]]
CNDRdimfun(Prog, reduction = "harmony", res="_0.3", split.by = c("X10X_batch", "library_batch", "DoubletFinder_call", "CMO_classification.global"))
#return markers
return_markersfun(Prog, assay="soupXcounts", features = Prog_soupXcounts_var_feat)

#saveRDS(Prog, "/Org_ChP/Org_ChP4/CellType/Firstcelltype/Prog.RDS")


###Prog cells spend more time in S phase

DimPlot(Prog, group.by = "Phase", raster = FALSE)
Ncellfun(Prog, "CMO per phase.csv", "Phase")



curr_clust <- 2
Prog_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]][0:30]
sum(Prog_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach, known_fibro_prog), na.rm=TRUE)
sum(Prog_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic), na.rm=TRUE)
sum(Prog_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Mural, known_peric, known_VSM), na.rm=TRUE)
sum(Prog_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Epithelial, known_myoepi, known_ependymal, known_ChPMature, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, Whyss_Correy_Epithelial_highGPX3, Whyss_Correy_Epithelial_highMT, Whyss_Correy_Epithelial_highSLC26A3), na.rm=TRUE)
sum(Prog_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% known_dividing, na.rm=TRUE)
sum(Prog_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog), na.rm=TRUE)
sum(Prog_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_epi_neuro_prog, known_neuro_prog), na.rm=TRUE)
sum(Prog_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS), na.rm=TRUE)
sum(Prog_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_astro, known_oligo, known_OPC), na.rm=TRUE)
sum(Prog_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(Cell_Stress, whole_cells), na.rm=TRUE)
sum(Prog_markers_list[[as.character(curr_clust)]][[2]][[1]][["gene"]] %in% c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B), na.rm=TRUE)

#Top20_Featuresfun(Prog, curr_clust)

search_markersfun(Prog, curr_clust, c(known_Mesen, known_Fibro, known_barrier, known_periv, known_menin, known_arach) %>% unique())
search_markersfun(Prog, curr_clust, known_fibro_prog)

search_markersfun(Prog, curr_clust, c(known_Endothelial, known_veinous, known_capilary, known_arterial, known_lymphatic) %>% unique())
search_markersfun(Prog, curr_clust, c(known_Mural, known_peric, known_VSM) %>% unique())

search_markersfun(Prog, curr_clust, c(known_Epithelial, Whyss_Correy_Epithelial) %>% unique())
search_markersfun(Prog, curr_clust, c(known_dark, other_dark) %>% unique())
search_markersfun(Prog, curr_clust, c(known_myoepi) %>% unique())
search_markersfun(Prog, curr_clust, c(known_light) %>% unique())
search_markersfun(Prog, curr_clust, c(known_ependymal, Whyss_Correy_Ependymal) %>% unique())
search_markersfun(Prog, curr_clust, Whyss_Correy_Epithelial_highGPX3 %>% unique())
search_markersfun(Prog, curr_clust, Whyss_Correy_Epithelial_highMT %>% unique())
search_markersfun(Prog, curr_clust, Whyss_Correy_Epithelial_highSLC26A3 %>% unique())

search_markersfun(Prog, curr_clust, known_dividing %>% unique())

search_markersfun(Prog, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Prog, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Prog, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Prog, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Prog, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Prog, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Prog, "TTR", raster=F)

0,1,3,5, Prog_Neuro_Epi
2 Prog_Neuro
4,6,9 Prog_Fib
7 Prog_Epi
8 Prog_Vasc

###Label

#First pass
Prog <- RenameIdents(Prog, "0"="Prog_Neuro_Epi", "1"="Prog_Neuro_Epi", "2"="Prog_Neuro", "3"="Prog_Neuro_Epi", "4"="Prog_Fib", "5"="Prog_Neuro_Epi", "6"="Prog_Fib", "7"="Prog_Epi", "8"="Prog_Vasc", "9"="Prog_Fib")
Prog@meta.data$Finalcelltype <- Prog@active.ident
Prog <- RunUMAP(object = Prog, dims = 1:20, assay = "soupXcounts", reduction =  "harmony", min.dist=.4)
DimPlot(Prog)

saveRDS(Prog, paste(curr_dir, "/CellType/Secondcelltype/Prog.RDS", sep=""))





##whole dataset

Prog <- readRDS(paste(curr_dir, "/CellType/Secondcelltype/Prog.RDS", sep=""))
Epi <- readRDS(paste(curr_dir, "/CellType/Secondcelltype/Epi.RDS", sep=""))
Fib <- readRDS(paste(curr_dir, "/CellType/Secondcelltype/Fibro.RDS", sep=""))
Neuro <- readRDS(paste(curr_dir, "/CellType/Secondcelltype/Neuro.RDS", sep=""))
Imm <- readRDS(paste(curr_dir, "/CellType/Secondcelltype/Imm.RDS", sep=""))
Dividing <- readRDS(paste(curr_dir, "/CellType/Secondcelltype/Dividing.RDS", sep=""))
Org_ChP <- readRDS(paste(curr_dir, "/CellType/Firstcelltype/Org_ChP.RDS", sep=""))

Prog_meta <- Prog@meta.data
Epi_meta <- Epi@meta.data
Fibro_meta <- Fib@meta.data
Neuro_meta <- Neuro@meta.data
Imm_meta <- Imm@meta.data
Dividing_meta <- Dividing@meta.data
merged_meta <- rbind(Prog_meta, Epi_meta, Fibro_meta, Neuro_meta, Imm_meta, Dividing_meta)

Org_ChP@meta.data$Finalcelltype <- NULL
Org_ChP_meta <- Org_ChP@meta.data

Org_ChP_meta <- rownames_to_column(Org_ChP_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$Finalcelltype)
colnames(merged_meta) <- c("rowname", "Finalcelltype")
Org_ChP_meta <- left_join(Org_ChP_meta, merged_meta, by="rowname")
Org_ChP_meta <- column_to_rownames(Org_ChP_meta, "rowname")
Org_ChP_meta$Subcelltype <- as.character(Org_ChP_meta$Finalcelltype)
Org_ChP_meta[is.na(Org_ChP_meta)==TRUE] <- "Doublet"
Org_ChP@meta.data <- Org_ChP_meta

Idents(Org_ChP) <- Org_ChP$Subcelltype
DimPlot(Org_ChP)

Org_ChP <- subset(Org_ChP, subset=Subcelltype=="Doublet", invert=TRUE)

rm(Prog_EN, Prog_E, Prog_F, Epi, Fib, Neuro, Imm, Vasc, Dividing)
gc()



NFSP_fun(Org_ChP, assay="soupXcounts")
use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
set.seed(1234)
Org_ChP <- RunHarmony(object = Org_ChP, group.by.vars = c("LANE"), dims.use = use_dims, assay="soupXcounts", max.iter.harmony = 50, plot_convergence = TRUE) #Correct for batches !!!!

CNDRfun(Org_ChP, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")
#Pick res
Idents(Org_ChP) <- Org_ChP@meta.data[["soupXcounts_snn_res.0.3"]]
CNDRdimfun(Org_ChP, reduction = "harmony", split.by = c("X10X_batch", "LANE", "sample", "CMO_classification.global", "DoubletFinder_call"), res = "0.3")

#Organize Sub and Major celltypes
Idents(Org_ChP) <- Org_ChP@meta.data[["Subcelltype"]]
levels(Org_ChP) <- c('Epi_1', 'Epi_2', 'Epi_3', 'Epi_4', 'Fib_1', 'Fib_2', 'Fib_3', 'BAM_1', 'BAM_2', 'BAM_3', "Prog_Vasc", "Prog_Fib", "Prog_Epi", "Prog_Neuro_Epi", "Prog_Neuro", "Dividing_Fib", "Dividing_Epi", "Dividing_Neuro_Epi", "Dividing_Neuro", 'Neuron')
Idents(Org_ChP) -> Org_ChP@meta.data[["Subcelltype"]]
DimPlot(Org_ChP, group.by = "Subcelltype", raster = FALSE)

Idents(Org_ChP) <- Org_ChP@meta.data[["Subcelltype"]]
Org_ChP <- RenameIdents(Org_ChP, 'Epi_1'="Epithelial", 'Epi_2'="Epithelial", 'Epi_3'="Epithelial", 'Epi_4'="Epithelial", 'Fib_1'="Fibroblast", 'Fib_2'="Fibroblast", 'Fib_3'="Fibroblast", 'BAM_1'="Immune", 'BAM_2'="Immune", 'BAM_3'="Immune", "Prog_Vasc"="Progenitor", "Prog_Fib"="Progenitor", "Prog_Epi"="Progenitor", "Prog_Neuro_Epi"="Progenitor", "Prog_Neuro"="Progenitor", "Dividing_Fib"="Dividing", "Dividing_Epi"="Dividing", "Dividing_Neuro_Epi"="Dividing", "Dividing_Neuro"="Dividing")
Idents(Org_ChP) -> Org_ChP@meta.data[["Majorcelltype"]]
DimPlot(Org_ChP, group.by = "Majorcelltype", raster = FALSE)

Idents(Org_ChP) <- Org_ChP@meta.data[["Subcelltype"]]
Org_ChP <- RenameIdents(Org_ChP, 'Epi_1'="Epithelial", 'Epi_2'="Epithelial", 'Epi_3'="Epithelial", 'Epi_4'="Epithelial", 'Fib_1'="Fibroblast", 'Fib_2'="Fibroblast", 'Fib_3'="Fibroblast", 'BAM_1'="Immune", 'BAM_2'="Immune", 'BAM_3'="Immune", "Prog_Vasc"="Immune", "Prog_Fib"="Fibroblast", "Prog_Epi"="Epithelial", "Prog_Neuro_Epi"="Epithelial", "Prog_Neuro"="Neuron", "Dividing_Fib"="Fibroblast", "Dividing_Epi"="Epithelial", "Dividing_Neuro_Epi"="Epithelial", "Dividing_Neuro"="Neuron")
Idents(Org_ChP) -> Org_ChP@meta.data[["Topcelltype"]]
DimPlot(Org_ChP, group.by = "Topcelltype", raster = FALSE)


##Remove stragglers

Idents(Org_ChP) <- Org_ChP@meta.data$Majorcelltype
Org_ChP <- CellSelector(DimPlot(Org_ChP), Org_ChP, ident = "Epi")
Epi <- subset(Org_ChP, ident="Epi")
Epi <- subset(Epi, Majorcelltype=="Epithelial")
Idents(Org_ChP) <- Org_ChP@meta.data$Majorcelltype
Org_ChP <- CellSelector(DimPlot(Org_ChP), Org_ChP, ident = "Fib")
Fib <- subset(Org_ChP, ident="Fib")
Fib <- subset(Fib, Majorcelltype=="Fibroblast")

Idents(Org_ChP) <- Org_ChP@meta.data$Majorcelltype
Prog <- subset(Org_ChP, Majorcelltype=="Progenitors" | Majorcelltype=="Dividing")

Idents(Org_ChP) <- Org_ChP@meta.data$Majorcelltype
Org_ChP <- CellSelector(DimPlot(Org_ChP), Org_ChP, ident = "Neur")
Neuro <- subset(Org_ChP, ident="Neur")
Neuro <- subset(Neuro, Majorcelltype=="Neuron")

Idents(Org_ChP) <- Org_ChP@meta.data$Majorcelltype
Org_ChP <- CellSelector(DimPlot(Org_ChP), Org_ChP, ident = "Imms")
Imm <- subset(Org_ChP, ident="Imms")
Imm <- subset(Imm, Majorcelltype=="Immune")


###label main dataset

Epi_meta <- Epi@meta.data
Fib_meta <- Fib@meta.data
Prog_meta <- Prog@meta.data
Neuro_meta <- Neuro@meta.data
Imm_meta <- Imm@meta.data
merged_meta <- rbind(Epi_meta, Fib_meta, Prog_meta, Neuro_meta, Imm_meta)

Org_ChP@meta.data$Majorcelltype <- NULL
Org_ChP_meta <- Org_ChP@meta.data

Org_ChP_meta <- rownames_to_column(Org_ChP_meta, "rowname")
merged_meta <- rownames_to_column(merged_meta, "rowname")
merged_meta <- data.frame(merged_meta$rowname, merged_meta$Majorcelltype)
colnames(merged_meta) <- c("rowname", "Majorcelltype")
Org_ChP_meta <- left_join(Org_ChP_meta, merged_meta, by="rowname")
Org_ChP_meta <- column_to_rownames(Org_ChP_meta, "rowname")
Org_ChP_meta$Majorcelltype <- as.character(Org_ChP_meta$Majorcelltype)
Org_ChP_meta[is.na(Org_ChP_meta)==TRUE] <- "Doublet"

Org_ChP@meta.data <- Org_ChP_meta
DimPlot(Org_ChP, group.by = "Majorcelltype", raster = FALSE)

Org_ChP <- subset(Org_ChP, subset=Majorcelltype=="Doublet", invert=TRUE)

rm(Fib, Epi, Prog, Neuro, Imm)
gc()


###Save and split subs

saveRDS(Org_ChP, paste(curr_dir, "/CellType/Finalcelltype/Org_ChP.RDS", sep=""))

Epi <- subset(Org_ChP, Majorcelltype=="Epithelial")
Fib <- subset(Org_ChP, Majorcelltype=="Fibroblast")
Imm <- subset(Org_ChP, Majorcelltype=="Immune")
Neuro <- subset(Org_ChP, Majorcelltype=="Neuron")






##Correlate to ex vivo
###load genesets

Epithelial <- read_csv("/Ex_ChP/sn-limma/Epithelial/Spreadsheet__Epithelial_corrected.csv")
Epithelial <- Epithelial[Epithelial$Percent_in_cells>.4,]
Epithelial <- Epithelial$gene
Epi_1a_genes <- read_csv("/Ex_ChP/sn-limma/Epi_1a/Spreadsheet__Epi_1a_corrected.csv")
Epi_1a_genes <- Epi_1a_genes[Epi_1a_genes$Percent_in_cells>.4,]
Epi_1a_genes <- Epi_1a_genes$gene
Epi_1b_genes <- read_csv("/Ex_ChP/sn-limma/Epi_1b/Spreadsheet__Epi_1b_corrected.csv")
Epi_1b_genes <- Epi_1b_genes[Epi_1b_genes$Percent_in_cells>.4,]
Epi_1b_genes <- Epi_1b_genes$gene
Epi_2a_genes <- read_csv("/Ex_ChP/sn-limma/Epi_2a/Spreadsheet__Epi_2a_corrected.csv")
Epi_2a_genes <- Epi_2a_genes[Epi_2a_genes$Percent_in_cells>.4,]
Epi_2a_genes <- Epi_2a_genes$gene
Epi_2b_genes <- read_csv("/Ex_ChP/sn-limma/Epi_2b/Spreadsheet__Epi_2b_corrected.csv")
Epi_2b_genes <- Epi_2b_genes[Epi_2b_genes$Percent_in_cells>.4,]
Epi_2b_genes <- Epi_2b_genes$gene


Fibroblast <- read_csv("/Ex_ChP/sn-limma/Fibroblast/Spreadsheet__Fibroblast_corrected.csv")
Fibroblast <- Fibroblast[Fibroblast$Percent_in_cells>.4,]
Fibroblast <- Fibroblast$gene
Fib_1_genes <- read_csv("/Ex_ChP/sn-limma/Fib_1/Spreadsheet__Fib_1_corrected.csv")
Fib_1_genes <- Fib_1_genes[Fib_1_genes$Percent_in_cells>.4,]
Fib_1_genes <- Fib_1_genes$gene
Fib_2_genes <- read_csv("/Ex_ChP/sn-limma/Fib_2/Spreadsheet__Fib_2_corrected.csv")
Fib_2_genes <- Fib_2_genes[Fib_2_genes$Percent_in_cells>.4,]
Fib_2_genes <- Fib_2_genes$gene
Fib_3_genes <- read_csv("/Ex_ChP/sn-limma/Fib_3/Spreadsheet__Fib_3_corrected.csv")
Fib_3_genes <- Fib_3_genes[Fib_3_genes$Percent_in_cells>.4,]
Fib_3_genes <- Fib_3_genes$gene

Immune_genes <- read_csv("/Ex_ChP/sn-limma/Immune/Spreadsheet__Immune_corrected.csv")
Immune_genes <- Immune_genes[Immune_genes$Percent_in_cells>.4,]
Immune_genes <- Immune_genes$gene
BAM_1_genes <- read_csv("/Ex_ChP/sn-limma/BAM_1/Spreadsheet__BAM_1_corrected.csv")
BAM_1_genes <- BAM_1_genes[BAM_1_genes$Percent_in_cells>.4,]
BAM_1_genes <- BAM_1_genes$gene
BAM_2_genes <- read_csv("/Ex_ChP/sn-limma/BAM_2/Spreadsheet__BAM_2_corrected.csv")
BAM_2_genes <- BAM_2_genes[BAM_2_genes$Percent_in_cells>.4,]
BAM_2_genes <- BAM_2_genes$gene
BAM_3_genes <- read_csv("/Ex_ChP/sn-limma/BAM_3/Spreadsheet__BAM_3_corrected.csv")
BAM_3_genes <- BAM_3_genes[BAM_3_genes$Percent_in_cells>.4,]
BAM_3_genes <- BAM_3_genes$gene
T_genes <- read_csv("/Ex_ChP/sn-limma/T_Cell/Spreadsheet__T_Cell_corrected.csv")
T_genes <- T_genes[T_genes$Percent_in_cells>.4,]
T_genes <- T_genes$gene

Neuro <- read_csv("/Ex_ChP/sn-limma/Neuro/Spreadsheet__Neuro_corrected.csv")
Neuro <- Neuro[Neuro$Percent_in_cells>.4,]
Neuro <- Neuro$gene

library(ComplexHeatmap)
library(colorRamp2)


####Major

Ex_ChP <- readRDS("/Ex_ChP/ref_sets/Ex_ChP.RDS")
Ex_ChP <- subset(Ex_ChP, Majorcelltype=="Epithelial" | Majorcelltype=="Fibroblast" | Majorcelltype=="Immune"| Subcelltype=="Neuro")
Ex_ChP <- subset(Ex_ChP, Subcelltype=="T_Cell", invert=T)
my_genes <- unique(c(Epithelial, Epi_1a_genes, Epi_1b_genes, Epi_2a_genes, Epi_2b_genes, Fibroblast, Fib_1_genes, Fib_2_genes, Fib_3_genes, Immune_genes, BAM_1_genes, BAM_2_genes, BAM_3_genes, Neuro))

#within overlapping genes
Ex_cts <- data.frame(AggregateExpression(Ex_ChP, group.by = "Majorcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Ex_cts) <- c("Epi_ExV", "Fib_ExV", "BAM_ExV", "Neuron_ExV")
Ex_cts <- data.frame(t(scale(t(Ex_cts))))
Ex_cts <- Ex_cts %>% rownames_to_column("gene")
Ex_cts[is.na(Ex_cts)] <- 0
Ex_cts <- Ex_cts[Ex_cts$gene %in% my_genes,]

Org_cts_T <- data.frame(AggregateExpression(Org_ChP, group.by = "Majorcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Org_cts_T) <- c("Epi_Tsai", "Fib_Tsai", "BAM_Tsai", "Prog", "Div", "Neuron_Tsai")
Org_cts_T <- Org_cts_T[,c("Epi_Tsai", "Fib_Tsai", "BAM_Tsai", "Neuron_Tsai")]
Org_cts_T <- data.frame(t(scale(t(Org_cts_T))))
Org_cts_T <- Org_cts_T %>% rownames_to_column("gene")
Org_cts_T <- Org_cts_T[Org_cts_T$gene %in% my_genes,]

cts <- left_join(Ex_cts, Org_cts_T, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c("Epi_Tsai", "Epi_ExV", "Fib_Tsai", "Fib_ExV", "BAM_Tsai", "BAM_ExV", "Neuron_Tsai", "Neuron_ExV")]
cts[is.na(cts)] <- 0

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("Tsai", rownames(cts_cor)), grepl("ExV", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-1, 0, 1)), width = ncol(cts_cor)*unit(2.7, "mm"), height = nrow(cts_cor)*unit(2.7, "mm"), col = colorRamp2(c(-.5, 0, .5), c("blue","white", "red")))
plot(plts)
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_Major_ExtoTsai.pdf", sep = ""))


####Sub

Ex_ChP <- readRDS("/Ex_ChP/ref_sets/Ex_ChP.RDS")
Ex_ChP <- subset(Ex_ChP, Majorcelltype=="Epithelial" | Majorcelltype=="Fibroblast" | Majorcelltype=="Immune")
Ex_ChP <- subset(Ex_ChP, Subcelltype=="T_Cell", invert=T)
my_genes <- unique(c(Epithelial, Epi_1a_genes, Epi_1b_genes, Epi_2a_genes, Epi_2b_genes, Fibroblast, Fib_1_genes, Fib_2_genes, Fib_3_genes, Immune_genes, BAM_1_genes, BAM_2_genes, BAM_3_genes))

#within overlapping genes
Ex_cts <- data.frame(AggregateExpression(Ex_ChP, group.by = "Subcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Ex_cts) <- c("BAM_1_ExV", "BAM_2_ExV", "BAM_3_ExV","Epi_1a_ExV", "Epi_1b_ExV", "Epi_2a_ExV", "Epi_2b_ExV", "Fib_1_ExV", "Fib_2_ExV", "Fib_3_ExV")
Ex_cts <- Ex_cts[,c("Epi_1a_ExV", "Epi_1b_ExV", "Epi_2a_ExV", "Epi_2b_ExV", "Fib_1_ExV", "Fib_2_ExV", "Fib_3_ExV", "BAM_1_ExV", "BAM_2_ExV", "BAM_3_ExV")]
Ex_cts <- data.frame(t(scale(t(Ex_cts))))
Ex_cts <- Ex_cts %>% rownames_to_column("gene")
Ex_cts[is.na(Ex_cts)] <- 0
Ex_cts <- Ex_cts[Ex_cts$gene %in% my_genes,]

Org_cts_T <- data.frame(AggregateExpression(Org_ChP, group.by = "Subcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Org_cts_T) <- c("Epi_1a_Tsai", "Epi_1b_Tsai", "Epi_2a_Tsai", "Epi_2b_Tsai", "Fib_1_Tsai", "Fib_2_Tsai", "Fib_3_Tsai", "BAM_1_Tsai", "BAM_2_Tsai", "BAM_3_Tsai")
Org_cts_T <- Org_cts_T[,c("Epi_1a_Tsai", "Epi_1b_Tsai", "Epi_2a_Tsai", "Epi_2b_Tsai", "Fib_1_Tsai", "Fib_2_Tsai", "Fib_3_Tsai", "BAM_1_Tsai", "BAM_2_Tsai", "BAM_3_Tsai")]
Org_cts_T <- data.frame(t(scale(t(Org_cts_T))))
Org_cts_T <- Org_cts_T %>% rownames_to_column("gene")
Org_cts_T <- Org_cts_T[Org_cts_T$gene %in% my_genes,]

cts <- left_join(Ex_cts, Org_cts_T, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c("Epi_1a_Tsai", "Epi_1a_ExV", "Epi_1b_Tsai", "Epi_1b_ExV", "Epi_2a_Tsai", "Epi_2a_ExV", "Epi_2b_Tsai", "Epi_2b_ExV", "Fib_1_Tsai", "Fib_1_ExV", "Fib_2_Tsai", "Fib_2_ExV", "Fib_3_Tsai", "Fib_3_ExV", "BAM_1_Tsai", "BAM_1_ExV", "BAM_2_Tsai", "BAM_2_ExV", "BAM_3_Tsai", "BAM_3_ExV")]
cts[is.na(cts)] <- 0

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("Tsai", rownames(cts_cor)), grepl("ExV", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-1, 0, 1)), width = ncol(cts_cor)*unit(2.7, "mm"), height = nrow(cts_cor)*unit(2.7, "mm"), col = colorRamp2(c(-.5, 0, .5), c("blue","white", "red")))
plot(plts)
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_Major_ExtoTsai.pdf", sep = ""))


####Epi

Epi_Ex <- readRDS("/Ex_ChP/ref_sets/Epi.rds")
my_genes <- unique(c(Epi_1a_genes, Epi_2b_genes, Epi_2a_genes, Epi_2b_genes))

#within overlapping genes
Ex_cts <- data.frame(AggregateExpression(Epi_Ex, group.by = "Subcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Ex_cts) <- c("Epi_1a_ExV", "Epi_1b_ExV", "Epi_2a_ExV", "Epi_2b_ExV")
Ex_cts <- data.frame(t(scale(t(Ex_cts))))
Ex_cts <- Ex_cts %>% rownames_to_column("gene")
Ex_cts[is.na(Ex_cts)] <- 0
Ex_cts <- Ex_cts[Ex_cts$gene %in% my_genes,]

Org_cts_T <- data.frame(AggregateExpression(Epi_T, group.by = "Subcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Org_cts_T) <- c("Epi_1a_Tsai", "Epi_1b_Tsai", "Epi_2a_Tsai", "Epi_2b_Tsai")
Org_cts_T <- Org_cts_T[,c("Epi_1a_Tsai", "Epi_1b_Tsai", "Epi_2a_Tsai", "Epi_2b_Tsai")]
Org_cts_T <- data.frame(t(scale(t(Org_cts_T))))
Org_cts_T <- Org_cts_T %>% rownames_to_column("gene")
Org_cts_T <- Org_cts_T[Org_cts_T$gene %in% my_genes,]

cts <- left_join(Ex_cts, Org_cts_T, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c("Epi_1a_Tsai", "Epi_1a_ExV", "Epi_1b_Tsai", "Epi_1b_ExV", "Epi_2a_Tsai", "Epi_2a_ExV", "Epi_2b_Tsai", "Epi_2b_ExV", )]
cts[is.na(cts)] <- 0

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("Tsai", rownames(cts_cor)), grepl("ExV", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-1, 0, 1)), width = ncol(cts_cor)*unit(2.7, "mm"), height = nrow(cts_cor)*unit(2.7, "mm"), col = colorRamp2(c(-.5, 0, .5), c("blue","white", "red")))
plot(plts)
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_Epi_ExtoTsai.pdf", sep = ""))


####Fib

Fibro_Ex <- readRDS("/Ex_ChP/ref_sets/Fibro.rds")
my_genes <- unique(c(Fib_1_genes, Fib_2_genes, Fib_3_genes))

#within overlapping genes
Ex_cts <- data.frame(AggregateExpression(Fib_Ex, group.by = "Subcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Ex_cts) <- c("Fib_1_ExV", "Fib_2_ExV", "Fib_3_ExV")
Ex_cts <- data.frame(t(scale(t(Ex_cts))))
Ex_cts <- Ex_cts %>% rownames_to_column("gene")
Ex_cts[is.na(Ex_cts)] <- 0
Ex_cts <- Ex_cts[Ex_cts$gene %in% my_genes,]

Org_cts_T <- data.frame(AggregateExpression(Fib_T, group.by = "Subcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Org_cts_T) <- c("Fib_1_Tsai", "Fib_2_Tsai", "Fib_3_Tsai")
Org_cts_T <- Org_cts_T[,c("Fib_1_Tsai", "Fib_2_Tsai", "Fib_3_Tsai")]
Org_cts_T <- data.frame(t(scale(t(Org_cts_T))))
Org_cts_T <- Org_cts_T %>% rownames_to_column("gene")
Org_cts_T <- Org_cts_T[Org_cts_T$gene %in% my_genes,]

cts <- left_join(Ex_cts, Org_cts_T, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c("Fib_1_Tsai", "Fib_1_ExV", "Fib_2_Tsai", "Fib_2_ExV", "Fib_3_Tsai", "Fib_3_ExV")]
cts[is.na(cts)] <- 0

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("Tsai", rownames(cts_cor)), grepl("ExV", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-1, 0, 1)), width = ncol(cts_cor)*unit(2.7, "mm"), height = nrow(cts_cor)*unit(2.7, "mm"), col = colorRamp2(c(-.5, 0, .5), c("blue","white", "red")))
plot(plts)
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_Fib_ExtoTsai.pdf", sep = ""))


####Imm

Imm_Ex <- readRDS("/Ex_ChP/ref_sets/Imm.rds")
my_genes <- unique(c(Immune_genes, BAM_1_genes, BAM_2_genes, BAM_3_genes, T_genes))

#within overlapping genes
Ex_cts <- data.frame(AggregateExpression(Imm_Ex, group.by = "Subcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Ex_cts) <- c("BAM_1_ExV", "BAM_2_ExV", "BAM_3_ExV", "T_Cell_ExV")
Ex_cts <- data.frame(t(scale(t(Ex_cts))))
Ex_cts <- Ex_cts %>% rownames_to_column("gene")
Ex_cts[is.na(Ex_cts)] <- 0
Ex_cts <- Ex_cts[Ex_cts$gene %in% my_genes,]

Org_cts_T <- data.frame(AggregateExpression(Imm_T, group.by = "Subcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Org_cts_T) <- c("BAM_1_Tsai", "BAM_2_Tsai", "BAM_3_Tsai")
Org_cts_T <- Org_cts_T[,c("BAM_1_Tsai", "BAM_2_Tsai", "BAM_3_Tsai")]
Org_cts_T <- data.frame(t(scale(t(Org_cts_T))))
Org_cts_T <- Org_cts_T %>% rownames_to_column("gene")
Org_cts_T <- Org_cts_T[Org_cts_T$gene %in% my_genes,]

cts <- left_join(Ex_cts, Org_cts_T, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c("BAM_1_Tsai", "BAM_1_ExV", "BAM_2_Tsai", "BAM_2_ExV", "BAM_3_Tsai", "BAM_3_ExV", "T_Cell_ExV")]
cts[is.na(cts)] <- 0

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("Tsai", rownames(cts_cor)), grepl("ExV", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-1, 0, 1)), width = ncol(cts_cor)*unit(2.7, "mm"), height = nrow(cts_cor)*unit(2.7, "mm"), col = colorRamp2(c(-.5, 0, .5), c("blue","white", "red")))
plot(plts)
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_Imm_ExtoTsai.pdf", sep = ""))



##Correlate to Lancaster
###load genesets

Epithelial <- read_csv("/Ex_ChP/sn-limma/Epithelial/Spreadsheet__Epithelial_corrected.csv")
Epithelial <- Epithelial[Epithelial$Percent_in_cells>.4,]
Epithelial <- Epithelial$gene
Epi_1a_genes <- read_csv("/Ex_ChP/sn-limma/Epi_1a/Spreadsheet__Epi_1a_corrected.csv")
Epi_1a_genes <- Epi_1a_genes[Epi_1a_genes$Percent_in_cells>.4,]
Epi_1a_genes <- Epi_1a_genes$gene
Epi_1b_genes <- read_csv("/Ex_ChP/sn-limma/Epi_1b/Spreadsheet__Epi_1b_corrected.csv")
Epi_1b_genes <- Epi_1b_genes[Epi_1b_genes$Percent_in_cells>.4,]
Epi_1b_genes <- Epi_1b_genes$gene
Epi_2a_genes <- read_csv("/Ex_ChP/sn-limma/Epi_2a/Spreadsheet__Epi_2a_corrected.csv")
Epi_2a_genes <- Epi_2a_genes[Epi_2a_genes$Percent_in_cells>.4,]
Epi_2a_genes <- Epi_2a_genes$gene
Epi_2b_genes <- read_csv("/Ex_ChP/sn-limma/Epi_2b/Spreadsheet__Epi_2b_corrected.csv")
Epi_2b_genes <- Epi_2b_genes[Epi_2b_genes$Percent_in_cells>.4,]
Epi_2b_genes <- Epi_2b_genes$gene


Fibroblast <- read_csv("/Ex_ChP/sn-limma/Fibroblast/Spreadsheet__Fibroblast_corrected.csv")
Fibroblast <- Fibroblast[Fibroblast$Percent_in_cells>.4,]
Fibroblast <- Fibroblast$gene
Fib_1_genes <- read_csv("/Ex_ChP/sn-limma/Fib_1/Spreadsheet__Fib_1_corrected.csv")
Fib_1_genes <- Fib_1_genes[Fib_1_genes$Percent_in_cells>.4,]
Fib_1_genes <- Fib_1_genes$gene
Fib_2_genes <- read_csv("/Ex_ChP/sn-limma/Fib_2/Spreadsheet__Fib_2_corrected.csv")
Fib_2_genes <- Fib_2_genes[Fib_2_genes$Percent_in_cells>.4,]
Fib_2_genes <- Fib_2_genes$gene
Fib_3_genes <- read_csv("/Ex_ChP/sn-limma/Fib_3/Spreadsheet__Fib_3_corrected.csv")
Fib_3_genes <- Fib_3_genes[Fib_3_genes$Percent_in_cells>.4,]
Fib_3_genes <- Fib_3_genes$gene

Immune_genes <- read_csv("/Ex_ChP/sn-limma/Immune/Spreadsheet__Immune_corrected.csv")
Immune_genes <- Immune_genes[Immune_genes$Percent_in_cells>.4,]
Immune_genes <- Immune_genes$gene
BAM_1_genes <- read_csv("/Ex_ChP/sn-limma/BAM_1/Spreadsheet__BAM_1_corrected.csv")
BAM_1_genes <- BAM_1_genes[BAM_1_genes$Percent_in_cells>.4,]
BAM_1_genes <- BAM_1_genes$gene
BAM_2_genes <- read_csv("/Ex_ChP/sn-limma/BAM_2/Spreadsheet__BAM_2_corrected.csv")
BAM_2_genes <- BAM_2_genes[BAM_2_genes$Percent_in_cells>.4,]
BAM_2_genes <- BAM_2_genes$gene
BAM_3_genes <- read_csv("/Ex_ChP/sn-limma/BAM_3/Spreadsheet__BAM_3_corrected.csv")
BAM_3_genes <- BAM_3_genes[BAM_3_genes$Percent_in_cells>.4,]
BAM_3_genes <- BAM_3_genes$gene
T_genes <- read_csv("/Ex_ChP/sn-limma/T_Cell/Spreadsheet__T_Cell_corrected.csv")
T_genes <- T_genes[T_genes$Percent_in_cells>.4,]
T_genes <- T_genes$gene

Neuro <- read_csv("/Ex_ChP/sn-limma/Neuro/Spreadsheet__Neuro_corrected.csv")
Neuro <- Neuro[Neuro$Percent_in_cells>.4,]
Neuro <- Neuro$gene

library(ComplexHeatmap)
library(colorRamp2)


####Major

Org_Ln <- readRDS("/Ex_ChP/ref_sets/Lancaster_ChPOrg.RDS")
my_genes <- unique(c(Epithelial, Epi_1a_genes, Epi_1b_genes, Epi_2a_genes, Epi_2b_genes, Fibroblast, Fib_1_genes, Fib_2_genes, Fib_3_genes, Immune_genes, BAM_1_genes, BAM_2_genes, BAM_3_genes, Neuro))

#within overlapping genes
Ex_cts <- data.frame(AggregateExpression(Org_Ln, group.by = "Lan_FirstCelltype", assay = "RNA", slot = "data", return.seurat = FALSE))
colnames(Ex_cts) <- c("Imm_ChP_Hem_Lan", "Neural_prog_Lan", "Mat_ChP_Lan", "Neurons_Lan", "ChP_stroma_Lan")
Ex_cts <- data.frame(t(scale(t(Ex_cts))))
Ex_cts <- Ex_cts %>% rownames_to_column("gene")
Ex_cts[is.na(Ex_cts)] <- 0
Ex_cts <- Ex_cts[Ex_cts$gene %in% my_genes,]

Org_cts_T <- data.frame(AggregateExpression(Org_ChP, group.by = "Majorcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Org_cts_T) <- c("Epi_Tsai", "Fib_Tsai", "Prog_Tsai", "Div_Tsai", "Neuron_Tsai")
Org_cts_T <- data.frame(t(scale(t(Org_cts_T))))
Org_cts_T <- Org_cts_T %>% rownames_to_column("gene")
Org_cts_T <- Org_cts_T[Org_cts_T$gene %in% my_genes,]

cts <- left_join(Ex_cts, Org_cts_T, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c("Epi_Tsai", "Mat_ChP_Lan", "Fib_Tsai", "ChP_stroma_Lan", "Neuron_Tsai", "Neurons_Lan", "Neural_prog_Lan", "Div_Tsai", "Imm_ChP_Hem_Lan", "Prog_Tsai")]
cts[is.na(cts)] <- 0

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("Tsai", rownames(cts_cor)), grepl("Lan", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-1, 0, 1)), width = ncol(cts_cor)*unit(2.7, "mm"), height = nrow(cts_cor)*unit(2.7, "mm"), col = colorRamp2(c(-.5, 0, .5), c("blue","white", "red")))
plot(plts)
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_Major_LancastertoTsai.pdf", sep = ""))



