###
### 023 Cell_type_annotation_ST.R
###
# Purpose:  Cell type annotation and cleaning of spatial transcriptomics data
# Dependencies:
source("/ChP_RNA/020 Cell_type_annotation_util.R", verbose = FALSE)
curr_dir <- "/Jackson" #output results to specific folder
set.seed(1234)
#Celltyping
##Main dataset
###Harmony

Org_ChP <- readRDS("/Ex_ChP/Org_ChP_BAM1/CellType/Firstcelltype/Organoid_Jackson.RDS")
Org_ChP <- subset(Org_ChP, CMO_classification.global=="Singlet")
gc()

ChP_genes <- c(known_arach, known_barrier, known_Fibro, known_Mesen, known_Mural, known_periv, known_Epithelial, known_dark, known_light, Whyss_Correy_Ependymal, Whyss_Correy_Epithelial, known_dividing, known_epi_neuro_prog, known_epi_prog, known_fibro_prog, known_mglia_prog, known_neuro_prog, known_Neuron, known_neuro_Exc, known_neuro_Inh, known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B)

NFSP_fun(Org_ChP, assay="soupXcounts", spikein= ChP_genes)
use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
set.seed(1234)
Org_ChP <- RunHarmony(object = Org_ChP, group.by.vars = c("Lane"), dims.use = use_dims, assay="soupXcounts", max.iter.harmony = 50, plot_convergence = TRUE) #Correct for batches !!!!

CNDRfun(Org_ChP, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")
#Pick res
Idents(Org_ChP) <- Org_ChP@meta.data[["soupXcounts_snn_res.0.3"]]
CNDRdimfun(Org_ChP, reduction = "harmony", split.by = c("X10X_batch", "Lane", "sample"), res = "0.3")

#return markers
return_markersfun(Org_ChP, assay="soupXcounts", features = Org_ChP_soupXcounts_var_feat, num_clust=0:8)



Ncellfun(Org_ChP, "CellType/Org_ChP/cluster per sample.csv", "sample")


###Dividing cells spend more time in S phase

Org_ChP <- CellCycleScoring(Org_ChP, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes, set.ident=FALSE)
plts <- DimPlot(Org_ChP, group.by = "Phase", raster = FALSE)
ggsave2(paste0(curr_dir, "/CellType/Org_ChP/Dividing.png"), plot = plts, width=14, height=6)
plot(plots)
Ncellfun(Org_ChP, "CMO per phase.csv", "Phase")


###Resolution 0.3

curr_clust <- 4
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

Featurefun(Org_ChP, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "MT-CO2", "KRT17", "KRT18", "KRT8", "RGN", "SI", "CFTR", "TRPM3", "KCNMA1", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "PECAM1", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS", "MT-CO2", "RSPO2", "RSPO3", "MKI67", "TUJ1"), "Test2",  assay="soupXcounts")


####Second pass
0: Prog
1: Fibro
2: 
  3: 
  4: Epi
5: 
  6: 
  7: 
  8: 
  9: 
  10: 
  11: 
  12: 
  13: 
  14: 
  
  ###First pass Isolate #5-6 Epis and keep high MT here, remove high MT everywhere else and rerun above
  
MT <- subset(Org_ChP, soupXcounts_snn_res.0.3=="5" | soupXcounts_snn_res.0.3=="6")
MT$temp <- "keep"
MT_less <- subset(Org_ChP, soupXcounts_snn_res.0.3=="5" | soupXcounts_snn_res.0.3=="6", invert=T)
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
MajorCelltype no MT no doublets
"0"="Fib"
"1"="Fib"
"2"="Epi"
"3"="Fib"
"4"="Neuro
"5"="Epi
"6"="Imm"
"7"="Mural"
"8"="Epi"

Progs probably with Fib.


Featurefun(Org_ChP, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "MT-CO2", "RGN", "SI", "CFTR", "TRPM3", "KCNMA1", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "PECAM1", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS", "RSPO2", "RSPO3", "MKI67", "TUJ1"), "Test2",  assay="soupXcounts")


###Rename all identities

#Idents(Org_ChP) <- Org_ChP@meta.data[["soupXcounts_snn_res.0.3"]]
Org_ChP <- RenameIdents(Org_ChP, '0'='Fib', '1'='Fib', '2'='Epi', '3'='Fib', '4'='Neuro', '5'='Epi', '6'='Imm', '7'='Mural', '8'='Epi')
Org_ChP@meta.data$Firstcelltype <- Idents(Org_ChP)

saveRDS(Org_ChP, paste(curr_dir, "/CellType/Firstcelltype/Org_ChP.RDS", sep=""))



##Subtypes

Org_ChP <- readRDS(paste(curr_dir, "/CellType/Firstcelltype/Org_ChP.RDS", sep=""))
for(i in unique(Org_ChP@meta.data$Firstcelltype)){
  sub <- subset(Org_ChP, Firstcelltype==i)
  saveRDS(sub, paste(curr_dir, "/CellType/Firstcelltype/", i, ".RDS", sep=""))
}



###Epi
####Remove doublets

Epi <- readRDS(paste(curr_dir, "/CellType/Firstcelltype/Epi.RDS", sep=""))
#Epi <- subset(Org_ChP, Majorcelltype=="Epi")

#Featurefun(Epi, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS", "MT-CO2"), "Test2",  assay="soupXcounts")

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
Epi <- RunHarmony(object = Epi, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Epi, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")

#Pick res
Idents(Epi) <- Epi@meta.data[["soupXcounts_snn_res.0.07"]]
CNDRdimfun(Epi, reduction = "harmony", res="_0.07", split.by = c("X10X_batch", "library_batch",  "CMO_classification.global", "Subcelltype"))
#return markers
return_markersfun(Epi, assay="soupXcounts", features = Epi_soupXcounts_var_feat)


####Dividing cells spend more time in S phase

DimPlot(Epi, group.by = "Phase", raster = FALSE)
Ncellfun(Epi, "CMO per phase.csv", "Phase")

Idents(Epi) <- Epi@meta.data$Firstcelltype
DimPlot(Epi)



curr_clust <- 5
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

search_markersfun(Epi, curr_clust, c(known_dividing, dividing_Org_Wang) %>% unique())

search_markersfun(Epi, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog) %>% unique())

search_markersfun(Epi, curr_clust, c(known_epi_neuro_prog, known_neuro_prog, prog_Org_Wang) %>% unique())

search_markersfun(Epi, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Epi, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Epi, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Epi, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Epi, c("MT-CO2", "HTR2C", "TTR", "CFAP54"), raster=F)


####Label

#First pass
Idents(Epi) <- Epi@meta.data[["soupXcounts_snn_res.0.3"]]
Epi <- RenameIdents(Epi, "0"="Progenitor_Epi", "1"="Progenitor_Epi", "2"="Epi_1b", "3"="Epi_1a", "4"="Epi_1a", "5"="Doublet", "6"="Progenitor", "7"="Epi_2b")
Epi@meta.data$Subcelltype <- Epi@active.ident

Epi <- subset(Epi, Subcelltype=="Doublet", invert=T)
saveRDS(Epi, paste(curr_dir, "/CellType/Secondcelltype/Epi.RDS", sep=""))

#Second pass
Idents(Epi) <- Epi@meta.data[["soupXcounts_snn_res.0.1"]]
Epi <- RenameIdents(Epi, "0"="Progenitor_Epi", "1"="Epi_1a", "2"="Epi_1b", "3"="Epi_2b", "4"="Epi_2a")
Epi@meta.data$Subcelltype <- Epi@active.ident
Prog_Epi <- subset(Epi, Subcelltype=="Progenitor_Epi")
Epi <- subset(Epi, Subcelltype=="Progenitor_Epi", invert=T)
saveRDS(Prog_Epi, paste(curr_dir, "/CellType/Finalcelltype/Prog_Epi.RDS", sep=""))

#Third pass
Idents(Epi) <- Epi@meta.data[["soupXcounts_snn_res.0.3"]]
Epi <- RenameIdents(Epi, "0"="Epi_1a", "1"="Epi_2b", "2"="Epi_1a", "3"="Epi_2a", "4"="Epi_1b", "5"="Epi_2b", "6"="Epi_1b", "7"="Epi_2a")
Epi@meta.data$Subcelltype <- Epi@active.ident
saveRDS(Epi, paste(curr_dir, "/CellType/Finalcelltype/Epi.RDS", sep=""))

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
Epi_2a <- subset(Epi, ident="Epi3")
Epi_2a <- subset(Epi_2a, Finalcelltype=="Epi_2a")

Idents(Epi) <- Epi@meta.data$Finalcelltype
Epi <- CellSelector(DimPlot(Epi), Epi, ident = "Epi4")
Epi_2b <- subset(Epi, ident="Epi4")
Epi_2b <- subset(Epi_2b, Finalcelltype=="Epi_2b")

Epi_1_meta <- Epi_1@meta.data
Epi_2_meta <- Epi_2@meta.data
Epi_2a_meta <- Epi_2a@meta.data
Epi_2b_meta <- Epi_2b@meta.data
merged_meta <- rbind(Epi_1_meta, Epi_2_meta, Epi_2a_meta, Epi_2b_meta)

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
levels(Epi) <- c("Epi_1", "Epi_2", "Epi_2a", "Epi_2b")
Epi@meta.data[["Finalcelltype"]] <- Idents(Epi)
Epi <- RunUMAP(object = Epi, dims = 1:20, assay = "soupXcounts", reduction =  "harmony", min.dist=1.4)
DimPlot(Epi, raster=F)

saveRDS(Epi, paste(curr_dir, "/CellType/Secondcelltype/Epi.RDS", sep=""))



##Fib

Fib <- readRDS(paste(curr_dir, "/CellType/Firstcelltype/Fib.RDS", sep=""))

#Featurefun(Fib, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="soupXcounts")


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
Fib <- RunHarmony(object = Fib, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Fib, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")

#Pick res
Idents(Fib) <- Fib@meta.data[["soupXcounts_snn_res.0.1"]]
CNDRdimfun(Fib, reduction = "harmony", res="_0.3", split.by = c("X10X_batch", "library_batch",  "CMO_classification.global"))
#return markers
return_markersfun(Fib, assay="soupXcounts", features = Fib_soupXcounts_var_feat)


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

search_markersfun(Fib, curr_clust, c(known_ChPImmature, known_epi_prog, known_epi_neuro_prog, prog_Org_Wang) %>% unique())

search_markersfun(Fib, curr_clust, c(known_epi_neuro_prog, known_neuro_prog) %>% unique())

search_markersfun(Fib, curr_clust, c(known_Neuron, known_neuro_Inh, known_neuro_Exc, known_neuro_Exc_subS, known_neuro_inh_CGE, known_neuro_inh_MGE, known_neuro_inh_subS) %>% unique())
search_markersfun(Fib, curr_clust, c(known_Immune, known_Macrophage, known_mglia_M0, known_mglia_prog, known_mglia_M1, known_mglia_M2, known_mglia_Dis, known_BAM, known_BAM_prog, known_BAM_strom, known_BAM_epi, known_BAM_MHCIIlo, known_BAM_MHCIIhi, known_Macro_M0, known_Macro_M1, known_Macro_M2, known_Macro_Dis, known_plasma_dendritic, known_conv_dendritic, known_mast, known_memory_T, known_effector_T, known_resident_T, known_NTK, known_plasma_B) %>% unique())

search_markersfun(Fib, curr_clust, c(known_astro, known_oligo, known_OPC) %>% unique())
search_markersfun(Fib, curr_clust, c(Cell_Stress, whole_cells) %>% unique())
FeaturePlot(Fib, c("LAMA2", "ADAM12", "PRICKLE1", "KCNMA1"), raster=F)


###Label

#First pass
Idents(Fib) <- Fib@meta.data[["soupXcounts_snn_res.0.07"]]
Fib <- RenameIdents(Fib, "0"="Prog2", "1"="Fib", "2"="Prog")
Fib@meta.data$Subcelltype <- Fib@active.ident
Prog <- subset(Fib, Subcelltype=="Prog" | Subcelltype=="Prog2")
Fib <- subset(Fib, Subcelltype=="Fib")

#Second pass
Fib@meta.data$Subcelltype <- "Fibroblast"
Fib@meta.data$Majorcelltype <- "Fibroblast"
saveRDS(Fib, paste0(curr_dir, "/CellType/Secondcelltype/Fib.RDS"))

DimPlot(Fib)


##Neuro

Neuro <- readRDS(paste(curr_dir, "/CellType/Firstcelltype/Neuro.RDS", sep=""))

#Featurefun(Neuro, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="soupXcounts")

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
Neuro <- RunHarmony(object = Neuro, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Neuro, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")

#Pick res
Idents(Neuro) <- Neuro@meta.data[["soupXcounts_snn_res.0.03"]]
CNDRdimfun(Neuro, reduction = "harmony", res="_0.03", split.by = c("X10X_batch", "library_batch",  "CMO_classification.global"))
#return markers
return_markersfun(Neuro, assay="soupXcounts", features = Neuro_soupXcounts_var_feat)

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

Neuro@meta.data$Subcelltype <- "Neuron"
Neuro@meta.data$Majorcelltype <- "Neuron"

saveRDS(Neuro, paste(curr_dir, "/CellType/Secondcelltype/Neuro.RDS", sep=""))



##BAM

Imm <- readRDS(paste(curr_dir, "/CellType/Firstcelltype/Imm.RDS", sep=""))

#Featurefun(Imm, c("HTR2C", "TTR", "CLIC6", "AQP1", "CFAP54", "FAM227A", "ARL13B", "TRPM3", "SLC4A4", "LAMA2", "LAMA4", "PRICKLE1", "CLDN1", "CLDN11", "CEMIP", "FOXP2", "COL6A3", "SLC6A13", "MECOM", "VWF", "TAGLN", "MYH11", "NOTCH3", "CADM2", "MAP2", "SNAP25", "GRIK2", "CLU", "AQP4", "MOBP", "MBP", "CDK19", "CD68", "F13A1", "P2RY12", "LYVE1", "STAB1", "CD163", "THEMIS"), "Test2",  assay="soupXcounts")

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
Imm <- RunHarmony(object = Imm, group.by.vars = c("Lane"), dims.use = use_dims, max.iter.harmony = 50, plot_convergence = TRUE, assay.use = "soupXcounts")
CNDRfun(Imm, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")

#Pick res
Idents(Imm) <- Imm@meta.data[["soupXcounts_snn_res.0.09"]]
CNDRdimfun(Imm, reduction = "harmony", res="_0.09", split.by = c("X10X_batch", "library_batch",  "CMO_classification.global"))
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
  
#First pass
Idents(Imm) <- Imm@meta.data[["soupXcounts_snn_res.0.09"]]
Imm <- RenameIdents(Imm, "0"="BAM_2", "1"="BAM_3", "2"="BAM_2")
Imm@meta.data$Subcelltype <- Idents(Imm)
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


##Whole dataset

Prog1 <- readRDS(paste(curr_dir, "/CellType/Secondcelltype/Prog_Fib1.RDS", sep=""))
Prog2 <- readRDS(paste(curr_dir, "/CellType/Secondcelltype/Prog_Fib2.RDS", sep=""))
Prog3 <- readRDS(paste(curr_dir, "/CellType/Secondcelltype/Prog_Epi.RDS", sep=""))
Epi <- readRDS(paste(curr_dir, "/CellType/Secondcelltype/Epi.RDS", sep=""))
Fib <- readRDS(paste(curr_dir, "/CellType/Secondcelltype/Fib.RDS", sep=""))
Neuro <- readRDS(paste(curr_dir, "/CellType/Secondcelltype/Neuro.RDS", sep=""))
Imm <- readRDS(paste(curr_dir, "/CellType/Secondcelltype/Imm.RDS", sep=""))
Org_ChP <- merge(Epi, c(Fib, Imm, Neuro, Prog1, Prog2, Prog3))



Org_ChP <- readRDS(paste(curr_dir, "/CellType/Finalcelltype/Org_ChP.RDS", sep=""))
NFSP_fun(Org_ChP, assay="soupXcounts")
use_dims <- 1:20 #!!!Modify based on elbow plot!!!!
set.seed(1234)
Org_ChP <- RunHarmony(object = Org_ChP, group.by.vars = c("Lane"), dims.use = use_dims, assay="soupXcounts", max.iter.harmony = 50, plot_convergence = TRUE) #Correct for batches !!!!

CNDRfun(Org_ChP, PCA_dims = use_dims, assay="soupXcounts", reduction = "harmony")
#Pick res
Idents(Org_ChP) <- Org_ChP@meta.data[["soupXcounts_snn_res.0.3"]]
CNDRdimfun(Org_ChP, reduction = "harmony", split.by = c("X10X_batch", "Lane", "sample", "CMO_classification.global", "Subcelltype", "Majorcelltype"), res = "0.3")

#Organize Sub and Major celltypes
Idents(Org_ChP) <- Org_ChP@meta.data[["Subcelltype"]]
levels(Org_ChP) <- c("Epi_1a", "Epi_1b", "Epi_2a", "Epi_2b", "Fibroblast", "BAM_1", "BAM_2", "BAM_3", "Neuron", "Progenitor", "Progenitor_Epi")
Idents(Org_ChP) -> Org_ChP@meta.data[["Subcelltype"]]
DimPlot(Org_ChP, group.by = "Subcelltype", raster = FALSE)

Idents(Org_ChP) <- Org_ChP@meta.data[["Subcelltype"]]
Org_ChP <- RenameIdents(Org_ChP, "Epi_1a"="Epithelial", "Epi_1b"="Epithelial", "Epi_2a"="Epithelial", "Epi_2b"="Epithelial", "BAM_1"="Immune", "BAM_2"="Immune", "BAM_3"="Immune", "Progenitor_Epi"="Progenitor")
levels(Org_ChP) <- c("Epithelial", "Fibroblast", "Immune", "Progenitor", "Neuron")
Idents(Org_ChP) -> Org_ChP@meta.data[["Majorcelltype"]]
DimPlot(Org_ChP, raster = FALSE)

saveRDS(Org_ChP, paste(curr_dir, "/CellType/Secondcelltype/Org_ChP.RDS", sep=""))


####Remove stragglers

Idents(Org_ChP) <- Org_ChP@meta.data$Majorcelltype
Org_ChP <- CellSelector(DimPlot(Org_ChP), Org_ChP, ident = "Epi")
Epi <- subset(Org_ChP, ident="Epi")
Epi <- subset(Epi, Majorcelltype=="Epithelial")

Idents(Org_ChP) <- Org_ChP@meta.data$Majorcelltype
Org_ChP <- CellSelector(DimPlot(Org_ChP), Org_ChP, ident = "Fib")
Fib <- subset(Org_ChP, ident="Fib")
Fib <- subset(Fib, Majorcelltype=="Fibroblast")

Idents(Org_ChP) <- Org_ChP@meta.data$Majorcelltype
Prog <- subset(Org_ChP, Majorcelltype=="Progenitor")

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

Idents(Org_ChP) <- Org_ChP@meta.data[["Lane"]]
Org_ChP <- RenameIdents(Org_ChP, "J1A_GEX"="Non_cocultured", "J1B_GEX"="Non_cocultured",  "J13A_GEX"="Non_cocultured", "J13B_GEX"="Non_cocultured", "J14A_GEX"="Non_cocultured", "J14B_GEX"="Non_cocultured", "J17A_GEX"="Non_cocultured", "J17B_GEX"="Non_cocultured", "B13A_GEX"="BAM_cocultured", "B13B_GEX"="BAM_cocultured", "B14A_GEX"="BAM_cocultured", "B14B_GEX"="BAM_cocultured", "B14C_GEX"="BAM_cocultured", "B17A_GEX"="BAM_cocultured", "B17B_GEX"="BAM_cocultured", "B17C_GEX"="BAM_cocultured")
Org_ChP_JBAM <- subset(Org_ChP, cocultured=="BAM_cocultured")
Org_ChP_J <- subset(Org_ChP_J, cocultured=="BAM_cocultured", invert=T)

saveRDS(Org_ChP, paste(curr_dir, "/CellType/Finalcelltype/Organoid_Jackson.RDS", sep=""))
saveRDS(Org_ChP_J, paste(curr_dir, "/CellType/Finalcelltype/Organoid_Jackson_non_cocultured.RDS", sep=""))
saveRDS(Org_ChP_JBAM, paste(curr_dir, "/CellType/Finalcelltype/Organoid_Jackson_BAM_cocultured.RDS", sep=""))

Epi_JBAM <- subset(Org_ChP_JBAM, Majorcelltype=="Epithelial")
Fib_JBAM <- subset(Org_ChP_JBAM, Majorcelltype=="Fibroblast")
Imm_JBAM <- subset(Org_ChP_JBAM, Majorcelltype=="Immune")



##Correlate to ex vivo
###load genesets

Epithelial <- read_csv("/Ex_ChP/sn-limma/Epithelial/Spreadsheet__Epithelial_corrected.csv")
Epithelial <- Epithelial[Epithelial$Percent_in_cells>.3,]
Epithelial <- Epithelial$gene
Epi_1a_genes <- read_csv("/Ex_ChP/sn-limma/Epi_1a/Spreadsheet__Epi_1a_corrected.csv")
Epi_1a_genes <- Epi_1a_genes[Epi_1a_genes$Percent_in_cells>.3,]
Epi_1a_genes <- Epi_1a_genes$gene
Epi_1b_genes <- read_csv("/Ex_ChP/sn-limma/Epi_1b/Spreadsheet__Epi_1b_corrected.csv")
Epi_1b_genes <- Epi_1b_genes[Epi_1b_genes$Percent_in_cells>.3,]
Epi_1b_genes <- Epi_1b_genes$gene
Epi_2a_genes <- read_csv("/Ex_ChP/sn-limma/Epi_2a/Spreadsheet__Epi_2a_corrected.csv")
Epi_2a_genes <- Epi_2a_genes[Epi_2a_genes$Percent_in_cells>.3,]
Epi_2a_genes <- Epi_2a_genes$gene
Epi_2b_genes <- read_csv("/Ex_ChP/sn-limma/Epi_2b/Spreadsheet__Epi_2b_corrected.csv")
Epi_2b_genes <- Epi_2b_genes[Epi_2b_genes$Percent_in_cells>.3,]
Epi_2b_genes <- Epi_2b_genes$gene


Fibroblast <- read_csv("/Ex_ChP/sn-limma/Fibroblast/Spreadsheet__Fibroblast_corrected.csv")
Fibroblast <- Fibroblast[Fibroblast$Percent_in_cells>.3,]
Fibroblast <- Fibroblast$gene
Fib_1_genes <- read_csv("/Ex_ChP/sn-limma/Fib_1/Spreadsheet__Fib_1_corrected.csv")
Fib_1_genes <- Fib_1_genes[Fib_1_genes$Percent_in_cells>.3,]
Fib_1_genes <- Fib_1_genes$gene
Fib_2_genes <- read_csv("/Ex_ChP/sn-limma/Fib_2/Spreadsheet__Fib_2_corrected.csv")
Fib_2_genes <- Fib_2_genes[Fib_2_genes$Percent_in_cells>.3,]
Fib_2_genes <- Fib_2_genes$gene
Fib_3_genes <- read_csv("/Ex_ChP/sn-limma/Fib_3/Spreadsheet__Fib_3_corrected.csv")
Fib_3_genes <- Fib_3_genes[Fib_3_genes$Percent_in_cells>.3,]
Fib_3_genes <- Fib_3_genes$gene

Immune_genes <- read_csv("/Ex_ChP/sn-limma/Immune/Spreadsheet__Immune_corrected.csv")
Immune_genes <- Immune_genes[Immune_genes$Percent_in_cells>.3,]
Immune_genes <- Immune_genes$gene
BAM_1_genes <- read_csv("/Ex_ChP/sn-limma/BAM_1/Spreadsheet__BAM_1_corrected.csv")
BAM_1_genes <- BAM_1_genes[BAM_1_genes$Percent_in_cells>.3,]
BAM_1_genes <- BAM_1_genes$gene
BAM_2_genes <- read_csv("/Ex_ChP/sn-limma/BAM_2/Spreadsheet__BAM_2_corrected.csv")
BAM_2_genes <- BAM_2_genes[BAM_2_genes$Percent_in_cells>.3,]
BAM_2_genes <- BAM_2_genes$gene
BAM_3_genes <- read_csv("/Ex_ChP/sn-limma/BAM_3/Spreadsheet__BAM_3_corrected.csv")
BAM_3_genes <- BAM_3_genes[BAM_3_genes$Percent_in_cells>.3,]
BAM_3_genes <- BAM_3_genes$gene

Neuro <- read_csv("/Ex_ChP/sn-limma/Neuro/Spreadsheet__Neuro_corrected.csv")
Neuro <- Neuro[Neuro$Percent_in_cells>.3,]
Neuro <- Neuro$gene

library(ComplexHeatmap)
library(colorRamp2)



####Major

Ex_ChP <- readRDS("/Ex_ChP/ref_sets/Ex_ChP.rds")
Ex_ChP <- subset(Ex_ChP, Majorcelltype=="Epithelial" | Majorcelltype=="Fibroblast" | Majorcelltype=="Immune"| Subcelltype=="Neuro")
Ex_ChP <- subset(Ex_ChP, Subcelltype=="T_Cell", invert=T)
my_genes <- unique(c(Epithelial, Epi_1a_genes, Epi_1b_genes, Epi_2a_genes, Epi_2b_genes, Fibroblast, Fib_1_genes, Fib_2_genes, Fib_3_genes, Immune_genes, BAM_1_genes, BAM_2_genes, BAM_3_genes))

#within overlapping genes
Ex_cts <- data.frame(AggregateExpression(Ex_ChP, group.by = "Majorcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Ex_cts) <- c("Epi_ExV", "Fib_ExV", "BAM_ExV", "Neuron_ExV")
Ex_cts <- data.frame(t(scale(t(Ex_cts))))
Ex_cts <- Ex_cts %>% rownames_to_column("gene")
Ex_cts[is.na(Ex_cts)] <- 0
Ex_cts <- Ex_cts[Ex_cts$gene %in% my_genes,]

Org_cts_J <- data.frame(AggregateExpression(Org_ChP_J, group.by = "Majorcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Org_cts_J) <- c("Epi_Jackson", "Fib_Jackson", "BAM_Jackson", "Neuron_Jackson", "Prog")
Org_cts_J <- Org_cts_J[,c("Epi_Jackson", "Fib_Jackson", "BAM_Jackson", "Neuron_Jackson")]
Org_cts_J <- data.frame(t(scale(t(Org_cts_J))))
Org_cts_J <- Org_cts_J %>% rownames_to_column("gene")
Org_cts_J <- Org_cts_J[Org_cts_J$gene %in% my_genes,]

cts <- left_join(Ex_cts, Org_cts_J, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c("Epi_Jackson", "Epi_ExV", "Fib_Jackson", "Fib_ExV", "BAM_Jackson", "BAM_ExV", "Neuron_Jackson", "Neuron_ExV")]
cts[is.na(cts)] <- 0

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("Jackson", rownames(cts_cor)), grepl("ExV", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-1, 0, 1)), width = ncol(cts_cor)*unit(2.7, "mm"), height = nrow(cts_cor)*unit(2.7, "mm"), col = colorRamp2(c(-.5, 0, .5), c("blue","white", "red")))
plot(plts)
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_Major_ExtoJ.pdf", sep = ""))


####Epi

Epi_Ex <- readRDS("/Ex_ChP/ref_sets/Epi.rds")
my_genes <- unique(c(Epi_1a_genes, Epi_1b_genes, Epi_2a_genes, Epi_2b_genes))

#within overlapping genes
Ex_cts <- data.frame(AggregateExpression(Epi_Ex, group.by = "Subcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Ex_cts) <- c("Epi_1a_ExV", "Epi_1b_ExV", "Epi_2a_ExV", "Epi_2b_ExV")
Ex_cts <- data.frame(t(scale(t(Ex_cts))))
Ex_cts <- Ex_cts %>% rownames_to_column("gene")
Ex_cts[is.na(Ex_cts)] <- 0
Ex_cts <- Ex_cts[Ex_cts$gene %in% my_genes,]

Org_cts_J <- data.frame(AggregateExpression(Epi_J, group.by = "Subcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Org_cts_J) <- c("Epi_1b_Jackson", "Epi_2b_Jackson")
Org_cts_J <- Org_cts_J[,c("Epi_1b_Jackson", "Epi_2b_Jackson")]
Org_cts_J <- data.frame(t(scale(t(Org_cts_J))))
Org_cts_J <- Org_cts_J %>% rownames_to_column("gene")
Org_cts_J <- Org_cts_J[Org_cts_J$gene %in% my_genes,]

cts <- left_join(Ex_cts, Org_cts_J, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c("Epi_1a_ExV", "Epi_1b_Jackson", "Epi_1b_ExV", "Epi_2a_ExV", "Epi_2b_Jackson", "Epi_2b_ExV")]
cts[is.na(cts)] <- 0

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("Jackson", rownames(cts_cor)), grepl("ExV", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-1, 0, 1)), width = ncol(cts_cor)*unit(2.7, "mm"), height = nrow(cts_cor)*unit(2.7, "mm"), col = colorRamp2(c(-.5, 0, .5), c("blue","white", "red")))
plot(plts)
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_Epi_ExtoJ.pdf", sep = ""))

####Imm

Imm_Ex <- readRDS("/Ex_ChP/ref_sets/Imm.rds")
my_genes <- unique(c(Immune_genes, BAM_1_genes, BAM_2_genes, BAM_3_genes))

#within overlapping genes
Ex_cts <- data.frame(AggregateExpression(Imm_Ex, group.by = "Subcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Ex_cts) <- c("BAM_1_ExV", "BAM_2_ExV", "BAM_3_ExV", "T_Cell_ExV")
Ex_cts <- data.frame(t(scale(t(Ex_cts))))
Ex_cts <- Ex_cts %>% rownames_to_column("gene")
Ex_cts[is.na(Ex_cts)] <- 0
Ex_cts <- Ex_cts[Ex_cts$gene %in% my_genes,]

Org_cts_J <- data.frame(AggregateExpression(Imm_J, group.by = "Subcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Org_cts_J) <- c("BAM_2_Jackson", "BAM_3_Jackson")
Org_cts_J <- Org_cts_J[,c("BAM_2_Jackson", "BAM_3_Jackson")]
Org_cts_J <- data.frame(t(scale(t(Org_cts_J))))
Org_cts_J <- Org_cts_J %>% rownames_to_column("gene")
Org_cts_J <- Org_cts_J[Org_cts_J$gene %in% my_genes,]

cts <- left_join(Ex_cts, Org_cts_J, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c("BAM_1_ExV", "BAM_2_Jackson", "BAM_2_ExV", "BAM_3_Jackson", "BAM_3_ExV", "T_Cell_ExV")]
cts[is.na(cts)] <- 0

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("Jackson", rownames(cts_cor)), grepl("ExV", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-1, 0, 1)), width = ncol(cts_cor)*unit(2.7, "mm"), height = nrow(cts_cor)*unit(2.7, "mm"), col = colorRamp2(c(-.5, 0, .5), c("blue","white", "red")))
plot(plts)
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_Imm_ExtoJ.pdf", sep = ""))


###Correlate Tsai to Jackson
####Major

my_genes <- unique(c(Epithelial, Epi_1_genes, Epi_2_genes, Epi_2a_genes, Epi_2b_genes, Fibroblast, Fib_1_genes, Fib_2_genes, Fib_3_genes, Immune_genes, BAM_1_genes, BAM_2_genes, BAM_3_genes))

Org_cts_J <- data.frame(AggregateExpression(Org_ChP_J, group.by = "Majorcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Org_cts_J) <- c("Epi_Jackson", "Fib_Jackson", "BAM_Jackson", "Neuron_Jackson")
Org_cts_J <- Org_cts_J[,c("Epi_Jackson", "Fib_Jackson", "BAM_Jackson", "Neuron_Jackson")]
Org_cts_J <- data.frame(t(scale(t(Org_cts_J))))
Org_cts_J <- Org_cts_J %>% rownames_to_column("gene")
Org_cts_J <- Org_cts_J[Org_cts_J$gene %in% my_genes,]

Org_cts_T <- data.frame(AggregateExpression(Org_ChP_T, group.by = "Majorcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Org_cts_T) <- c("Epi_Tsai", "Fib_Tsai", "BAM_Tsai", "Prog", "Div", "Neuron_Tsai")
Org_cts_T <- Org_cts_T[,c("Epi_Tsai", "Fib_Tsai", "BAM_Tsai", "Neuron_Tsai")]
Org_cts_T <- data.frame(t(scale(t(Org_cts_T))))
Org_cts_T <- Org_cts_T %>% rownames_to_column("gene")
Org_cts_T <- Org_cts_T[Org_cts_T$gene %in% my_genes,]

cts <- left_join(Org_cts_J, Org_cts_T, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c("Epi_Tsai", "Epi_Jackson", "Fib_Tsai", "Fib_Jackson", "BAM_Tsai", "BAM_Jackson", "Neuron_Tsai", "Neuron_Jackson")]
cts[is.na(cts)] <- 0

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("Tsai", rownames(cts_cor)), grepl("Jackson", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-1, 0, 1)), width = ncol(cts_cor)*unit(2.7, "mm"), height = nrow(cts_cor)*unit(2.7, "mm"), col = colorRamp2(c(-.5, 0, .5), c("blue","white", "red")))
plot(plts)
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_Major_TtoJ.pdf", sep = ""), width=7, height=7)


####Epi

my_genes <- unique(c(Epi_1a_genes, Epi_1b_genes, Epi_2a_genes, Epi_2b_genes))
#within overlapping genes
Org_cts_J <- data.frame(AggregateExpression(Epi_J, group.by = "Subcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Org_cts_J) <- c("Epi_1b_Jackson", "Epi_2b_Jackson")
Org_cts_J <- data.frame(t(scale(t(Org_cts_J))))
Org_cts_J <- Org_cts_J %>% rownames_to_column("gene")
Org_cts_J[is.na(Org_cts_J)] <- 0
Org_cts_J <- Org_cts_J[Org_cts_J$gene %in% my_genes,]

Org_cts_T <- data.frame(AggregateExpression(Epi_T, group.by = "Subcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Org_cts_T) <- c("Epi_1a_Tsai", "Epi_1b_Tsai", "Epi_2a_Tsai", "Epi_2b_Tsai")
Org_cts_T <- Org_cts_T[,c("Epi_1a_Tsai", "Epi_1b_Tsai", "Epi_2a_Tsai", "Epi_2b_Tsai")]
Org_cts_T <- data.frame(t(scale(t(Org_cts_T))))
Org_cts_T <- Org_cts_T %>% rownames_to_column("gene")
Org_cts_T <- Org_cts_T[Org_cts_T$gene %in% my_genes,]

cts <- left_join(Org_cts_J, Org_cts_T, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c("Epi_1a_Tsai", "Epi_1b_Tsai", "Epi_1b_Jackson", "Epi_2a_Tsai", "Epi_2b_Tsai", "Epi_2b_Jackson")]
cts[is.na(cts)] <- 0

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("Tsai", rownames(cts_cor)), grepl("Jackson", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-1, 0, 1)), width = ncol(cts_cor)*unit(2.7, "mm"), height = nrow(cts_cor)*unit(2.7, "mm"), col = colorRamp2(c(-.5, 0, .5), c("blue","white", "red")))
plot(plts)
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_Epi_TtoJ.pdf", sep = ""), width=7, height=7)


####Imm

my_genes <- unique(c(Immune_genes, BAM_1_genes, BAM_2_genes, BAM_3_genes))

#within overlapping genes
Org_cts_J <- data.frame(AggregateExpression(Imm_J, group.by = "Subcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Org_cts_J) <- c("BAM_2_Jackson", "BAM_3_Jackson")
Org_cts_J <- data.frame(t(scale(t(Org_cts_J))))
Org_cts_J <- Org_cts_J %>% rownames_to_column("gene")
Org_cts_J[is.na(Org_cts_J)] <- 0
Org_cts_J <- Org_cts_J[Org_cts_J$gene %in% my_genes,]

Org_cts_T <- data.frame(AggregateExpression(Imm, group.by = "Subcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Org_cts_T) <- c("BAM_1_Tsai", "BAM_3_Tsai", "BAM_2_Tsai")
Org_cts_T <- Org_cts_T[,c("BAM_2_Tsai", "BAM_3_Tsai")]
Org_cts_T <- data.frame(t(scale(t(Org_cts_T))))
Org_cts_T <- Org_cts_T %>% rownames_to_column("gene")
Org_cts_T <- Org_cts_T[Org_cts_T$gene %in% my_genes,]

cts <- left_join(Org_cts_J, Org_cts_T, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c("BAM_1_Tsai", "BAM_2_Tsai", "BAM_2_Jackson", "BAM_3_Tsai", "BAM_3_Jackson")]
cts[is.na(cts)] <- 0

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("Tsai", rownames(cts_cor)), grepl("Jackson", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-1, 0, 1)), width = ncol(cts_cor)*unit(2.7, "mm"), height = nrow(cts_cor)*unit(2.7, "mm"), col = colorRamp2(c(-.5, 0, .5), c("blue","white", "red")))
plot(plts)
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_Imm_TtoJ.pdf", sep = ""), width=7, height=7)



###Correlate to Lancaster

Org_Ln <- readRDS("/Ex_ChP/ref_sets/Lancaster_ChPOrg.RDS")
my_genes <- unique(c(Epithelial, Epi_1a_genes, Epi_1b_genes, Epi_2a_genes, Epi_2b_genes, Fibroblast, Fib_1_genes, Fib_2_genes, Fib_3_genes, Immune_genes, BAM_1_genes, BAM_2_genes, BAM_3_genes, Neuro))

#within overlapping genes
Ex_cts <- data.frame(AggregateExpression(Org_Ln, group.by = "Lan_FirstCelltype", assay = "RNA", slot = "data", return.seurat = FALSE))
colnames(Ex_cts) <- c("Imm_ChP_Hem_Lan", "Neural_prog_Lan", "Mat_ChP_Lan", "Neurons_Lan", "ChP_stroma_Lan")
Ex_cts <- data.frame(t(scale(t(Ex_cts))))
Ex_cts <- Ex_cts %>% rownames_to_column("gene")
Ex_cts[is.na(Ex_cts)] <- 0
Ex_cts <- Ex_cts[Ex_cts$gene %in% my_genes,]

Org_cts <- data.frame(AggregateExpression(Org_ChP_J, group.by = "Majorcelltype", assay = "soupXcounts", slot = "data", return.seurat = FALSE))
colnames(Org_cts) <- c("Epi_Jackson", "Fib_Jackson", "BAM_Jackson", "Neuron_Jackson", "Prog_Jackson")
Org_cts <- data.frame(t(scale(t(Org_cts))))
Org_cts <- Org_cts %>% rownames_to_column("gene")
Org_cts <- Org_cts[Org_cts$gene %in% my_genes,]

cts <- left_join(Ex_cts, Org_cts, by = "gene")
cts <- cts %>% column_to_rownames("gene")
colnames(cts)
cts <- cts[,c("Epi_Jackson", "Mat_ChP_Lan", "Fib_Jackson", "ChP_stroma_Lan", "Neuron_Jackson", "Neurons_Lan", "Neural_prog_Lan", "Imm_ChP_Hem_Lan", "Prog_Jackson")]
cts[is.na(cts)] <- 0

cts_cor <- cor(cts)
cts_cor <- cts_cor[grepl("Jackson", rownames(cts_cor)), grepl("Lan", colnames(cts_cor))]
plts <- Heatmap(cts_cor, cluster_rows = F, cluster_columns = F, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-1, 0, 1)), width = ncol(cts_cor)*unit(2.7, "mm"), height = nrow(cts_cor)*unit(2.7, "mm"), col = colorRamp2(c(-.5, 0, .5), c("blue","white", "red")))
plot(plts)
graph2pdf(plts, file=paste(curr_dir, "/CellType/Publish/Correlation_heatmap_Major_LancastertoJackson.pdf", sep = ""))

