###
### 034 differential_RNA.R
###
# Purpose: differential in single nuclei RNA dataset
# Dependencies:
source("/030 limma_pattern_fgsea_util.R", verbose = FALSE)
curr_dir <- "/ChP_RNA" #output results to specific folder
set.seed(42)


#Set and load dependent variables for limma/DESeq functions
Gene.search <- read.csv("/DEG to GO.csv")
FDRcutoff <- 0.05
LFcutoff <- 0.25
pct_cutoff <- 0.05

###limma run
####Sub

#!!!!Modify Majorcelltype/Subcelltype in limmafun above
unique(Org_ChP@meta.data$Finalcelltype)

limmafun_Org(srt=Org_ChP, celltype="Epi_1a")
limmafun_Org(srt=Org_ChP, celltype="Epi_1b")
limmafun_Org(srt=Org_ChP, celltype="Epi_2a")
limmafun_Org(srt=Org_ChP, celltype="Epi_2b")

limmafun_Org(srt=Org_ChP, celltype="Prog_Epi")
limmafun_Org(srt=Org_ChP, celltype="Prog_Neuro_Epi")

#limmafun_Org(srt=Org_ChP, celltype="Dividing_Epi")
limmafun_Org(srt=Org_ChP, celltype="Dividing_Neuro_Epi")

limmafun_Org(srt=Org_ChP, celltype="Fib_1")
limmafun_Org(srt=Org_ChP, celltype="Fib_2")
limmafun_Org(srt=Org_ChP, celltype="Fib_3")

limmafun_Org(srt=Org_ChP, celltype="Prog_Fib")
limmafun_Org(srt=Org_ChP, celltype="Dividing_Fib")

limmafun_Org(srt=Org_ChP, celltype="BAM_1")
#limmafun_Org(srt=Org_ChP, celltype="BAM_2")
#limmafun_Org(srt=Org_ChP, celltype="BAM_3")

#limmafun_Org(srt=Org_ChP, celltype="Prog_Vasc")

limmafun_Org(srt=Org_ChP, celltype="Neuron")
limmafun_Org(srt=Org_ChP, celltype="Prog_Neuro")
#limmafun_Org(srt=Org_ChP, celltype="Dividing_Neuro")


####Major

#!!!!Modify Majorcelltype/Subcelltype in limmafun_Org above
unique(Org_ChP@meta.data$Majorcelltype)
limmafun_Org(srt=Org_ChP, celltype="Epithelial")
limmafun_Org(srt=Org_ChP, celltype="Fibroblast")
limmafun_Org(srt=Org_ChP, celltype="Immune")
limmafun_Org(srt=Org_ChP, celltype="Progenitor")
limmafun_Org(srt=Org_ChP, celltype="Dividing")



###GO run

Path_fun(limma_MainGenotype_Epi_1a)
Path_fun(limma_MainGenotype_Epi_1b)
Path_fun(limma_MainGenotype_Epi_2a)
Path_fun(limma_MainGenotype_Epi_2b)

Path_fun(limma_MainGenotype_Prog_Epi)
Path_fun(limma_MainGenotype_Prog_Neuro_Epi)

#Path_fun(limma_MainGenotype_Dividing_Epi)
Path_fun(limma_MainGenotype_Dividing_Neuro_Epi)

Path_fun(limma_MainGenotype_Fib_1)
Path_fun(limma_MainGenotype_Fib_2)
Path_fun(limma_MainGenotype_Fib_3)
Path_fun(limma_MainGenotype_Prog_Fib)
Path_fun(limma_MainGenotype_Dividing_Fib)

Path_fun(limma_MainGenotype_BAM_1)
#Path_fun(limma_MainGenotype_BAM_2)
#Path_fun(limma_MainGenotype_BAM_3)

#Path_fun(limma_MainGenotype_Prog_Vasc)

Path_fun(limma_MainGenotype_Neuron)
Path_fun(limma_MainGenotype_Prog_Neuro)
#Path_fun(limma_MainGenotype_Dividing_Neuro)

Path_fun(limma_MainGenotype_Epithelial)
Path_fun(limma_MainGenotype_Fibroblast)
Path_fun(limma_MainGenotype_Immune)
Path_fun(limma_MainGenotype_Progenitor)
Path_fun(limma_MainGenotype_Dividing)



#Figures DEGs

library(ggplot2)
library(ggrepel)
FDRcutoff <- 0.01
LFcutoff <- 0.25
Epi_genes <- c("TUBB4B", "CFAP54", "DYNC2H1", "ARL13B", #cilia
               "INSIG1", "LAMP2", "LDLRAD4", "LDLR", "LRP1", "SORL1", #lipid
               "PFKL", "SLC2A1") #glucose and glut1
Fib_genes <- c("LUM", "COL7A1", "CDH11", "CNTN5", "SMAD3", "ITGA111", #adhesion
               "VLDLR-AS", "LRP1", "LRPPRC", "LRP1B", #lipid
               "GRIA2", "PIEZO2", "CACHD1", "KCNT2") #transport
BAM_genes <- c("PSEN1", "IL4R", "IL2R", "CD28", "MERTK", "CD36", "CLEC7A", "ITGAM")


###custom Volcano

Volcano_fun <- function(TT_res, genes){
  name <- deparse(substitute(TT_res))
  name <- str_sub(name, 7)
  df.1b <- TT_res
  df.1b <- df.1b[!duplicated(df.1b$gene),]
  rownames(df.1b) <- NULL
  df.1b <- df.1b %>% column_to_rownames("gene")
  df.1b$diffexpressed <- "NS"
  df.1b$diffexpressed[df.1b$logFC > LFcutoff & df.1b$adj.P.Val < FDRcutoff] <- "UP"
  df.1b$diffexpressed[df.1b$logFC < -LFcutoff & df.1b$adj.P.Val < FDRcutoff] <- "DOWN"
  mycolors <- c("blue", "red", "black")
  names(mycolors) <- c("DOWN", "UP", "NS")
  
  lbl_gene <- df.1b[rownames(df.1b) %in% genes,]
  
  plts <- ggplot(data=df.1b, aes(x=logFC, y=-log10(adj.P.Val), col = diffexpressed)) +
    ggrastr::rasterise(geom_point(size=.1), dpi=600) + theme_bw() + xlim(-2,2) + ylim(0,8) + 
    scale_color_manual(values=mycolors) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.position = "none", axis.text = element_text(size = 7), axis.title=element_text(size=7)) +
    geom_vline(xintercept=c(-LFcutoff, LFcutoff), col="red") +
    geom_hline(yintercept=-log10(FDRcutoff), col="red") + geom_text_repel(data=lbl_gene, mapping=aes(label=rownames(lbl_gene)), size=2, max.overlaps=1000)
  plot(plts)
  graph2pdf(plts, paste(curr_dir, "/limma/volcano_", name, ".pdf", sep=""), width=1.574, height=1.574)
}

Volcano_fun(limma_MainGenotype_Epithelial, Epi_genes)
Volcano_fun(limma_MainGenotype_Epi_1a, Epi_genes)
Volcano_fun(limma_MainGenotype_Epi_1b, Epi_genes)
Volcano_fun(limma_MainGenotype_Epi_2a, Epi_genes)
Volcano_fun(limma_MainGenotype_Epi_2b, Epi_genes)

Volcano_fun(limma_MainGenotype_Fib_1, Fib_genes)
Volcano_fun(limma_MainGenotype_Fib_2, Fib_genes)
Volcano_fun(limma_MainGenotype_Fib_3, Fib_genes)
Volcano_fun(limma_MainGenotype_Fibroblast, Fib_genes)

Volcano_fun(limma_MainGenotype_Immune, BAM_genes)

Volcano_fun(limma_MainGenotype_iCSF, c("EVPL", "BAIAP2", "PPL", "CGN", "LAYN", "KRT18", "NECTIN4", "ADAMTS4", "MUC16", "IL7", "ADGRE5", "TNFRSF10C", "PADI2"))


##upset plots

library(dplyr)
library(ggpmisc)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(DESeq2)
library(viridis)
library(export)
col_fun <- col_fun


###Epithelial

cm <- make_comb_mat(list_to_matrix(list(Ex_vivo_ADemvNCI = c(limma_cogdx_33_ADvNCI_Epithelial_sig$gene), Ex_vivo_NCI_APOE34v33 = c(limma_apoe_ADNC034v33_Epithelial_sig$gene), In_vitro_ctrl_APOE44v33 = c(MainGenotype_Epithelial_sig$gene))))
plts <- UpSet(cm, set_order = rownames(cm), pt_size = unit(2, "mm"), top_annotation = upset_top_annotation(cm, height = unit(2.5, "inch"), annotation_name_rot = 90))
plts
pdf(file=paste("/Org_ChP/Org_ChP9/limma/upset_Epithelial.pdf"), width=4, height=4)
plts
dev.off()


###Fibroblast

cm <- make_comb_mat(list_to_matrix(list(Ex_vivo_ADemvNCI = c(limma_cogdx_33_ADvNCI_Fibroblast_sig$gene), Ex_vivo_NCI_APOE34v33 = c(limma_apoe_ADNC034v33_Fibroblast_sig$gene), In_vitro_ctrl_APOE44v33 = c(MainGenotype_Fibroblast_sig$gene))))
plts <- UpSet(cm, set_order = rownames(cm), pt_size = unit(2, "mm"), top_annotation = upset_top_annotation(cm, height = unit(2.5, "inch"), annotation_name_rot = 90))
plts
pdf(file=paste("/Org_ChP/Org_ChP9/limma/upset_Fibroblast.pdf"), width=4, height=4)
plts
dev.off()





##GO heatmap of overlap

library(ComplexHeatmap)
library(circlize)
library(DESeq2)
library(viridis)
library(export)
col_fun <- colorRamp2(c(-2, 0, 2), c("blue","white", "red"))


###Overall GO

selection <- c(      #IMMUNE
  "GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS", "GOBP_ADAPTIVE_IMMUNE_RESPONSE", "GOBP_REGULATION_OF_IMMUNE_RESPONSE", "GOMF_IMMUNE_RECEPTOR_ACTIVITY", "GOBP_INFLAMMATORY_RESPONSE", #general immune
  "GOBP_CELL_CHEMOTAXIS", "GOBP_RESPONSE_TO_CYTOKINE", "GOBP_REGULATION_OF_RESPONSE_TO_STRESS", "GOBP_RESPONSE_TO_WOUNDING", "GOMF_CYTOKINE_RECEPTOR_ACTIVITY", #response immune
  "GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE", "GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE", "GOBP_REGULATION_OF_T_CELL_ACTIVATION", "GOBP_T_CELL_DIFFERENTIATION_INVOLVED_IN_IMMUNE_RESPONSE", "GOBP_T_CELL_ACTIVATION", "GOBP_MYELOID_LEUKOCYTE_ACTIVATION", "GOBP_B_CELL_ACTIVATION", "GOBP_B_CELL_MEDIATED_IMMUNITY", "GOBP_LEUKOCYTE_PROLIFERATION", "GOBP_REGULATION_OF_LEUKOCYTE_PROLIFERATION", #Tcell, Bcell, 
  "GOCC_MHC_PROTEIN_COMPLEX",  "GOMF_MHC_CLASS_II_PROTEIN_COMPLEX_BINDING", "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I", "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION", "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_ANTIGEN", "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN",  "GOMF_PEPTIDE_ANTIGEN_BINDING", "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN", #MHC to antigen
  "GOBP_INTERLEUKIN_8_PRODUCTION", "GOBP_INTERLEUKIN_6_PRODUCTION", "GOBP_PHAGOCYTOSIS",
  "GOMF_UNFOLDED_PROTEIN_BINDING", "GOMF_HSP90_PROTEIN_BINDING", #HSP
  
  #ADHESION
  "GOBP_LEUKOCYTE_CELL_CELL_ADHESION", "GOMF_STRUCTURAL_MOLECULE_ACTIVITY", "GOMF_CELL_ADHESION_MOLECULE_BINDING", "GOBP_CELL_CELL_ADHESION", "GOBP_CELL_ADHESION_MEDIATED_BY_INTEGRIN", "GOBP_CELL_MATRIX_ADHESION", "GOBP_CELL_SUBSTRATE_ADHESION", "GOBP_CELL_ADHESION", #adhesion
  "GOCC_CELL_CELL_JUNCTION", "GOCC_ANCHORING_JUNCTION","GOCC_TIGHT_JUNCTION", "GOBP_TIGHT_JUNCTION_ORGANIZATION", "GOCC_APICAL_JUNCTION_COMPLEX", "GOCC_CELL_SUBSTRATE_JUNCTION",#tight junctions
  "GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX",  "GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY", "GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT", "GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT_CONFERRING_TENSILE_STRENGTH", #ECM
  "GOCC_ACTIN_CYTOSKELETON", "GOBP_ACTIN_FILAMENT_BASED_PROCESS", "GOCC_ACTIN_FILAMENT_BUNDLE", "GOBP_ACTIN_FILAMENT_ORGANIZATION", "GOMF_ACTIN_BINDING", #actin
  "GOCC_POSTSYNAPTIC_DENSITY_MEMBRANE", "GOMF_CADHERIN_BINDING", #other adhesion
  
  #EXPORT/SECRETION/TRANSPORT
  "GOBP_EXPORT_FROM_CELL", "GOBP_VESICLE_MEDIATED_TRANSPORT", "GOBP_EXOCYTOSIS", "GOBP_POSITIVE_REGULATION_OF_SECRETION", "GOBP_SECRETION", "GOBP_NEUROTRANSMITTER_SECRETION", "GOBP_REGULATION_OF_PROTEIN_SECRETION", "GOBP_NEGATIVE_REGULATION_OF_HORMONE_SECRETION", "GOBP_EXTRACELLULAR_TRANSPORT",  #export/secretion
  "GOBP_REGULATION_OF_TRANSPORT", "GOMF_MOLECULAR_TRANSDUCER_ACTIVITY", "GOBP_TRANSMEMBRANE_TRANSPORT", "GOMF_PASSIVE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY", "GOMF_TRANSMEMBRANE_TRANSPORTER_BINDING",#general transport
  "GOCC_CATION_CHANNEL_COMPLEX", "GOMF_INORGANIC_MOLECULAR_ENTITY_TRANSMEMBRANE_TRANSPORTER_ACTIVITY", "GOCC_POTASSIUM_CHANNEL_COMPLEX", "GOBP_LIGAND_GATED_ION_CHANNEL_SIGNALING_PATHWAY", "GOMF_TRANSPORTER_ACTIVITY", #ion transport
  "GOMF_GATED_CHANNEL_ACTIVITY", "GOMF_VOLTAGE_GATED_CHANNEL_ACTIVITY", "GOMF_SIGNALING_RECEPTOR_BINDING", "GOBP_ACTIVATION_OF_GTPASE_ACTIVITY",#transport related to signaling
  
  #METABOLISM
  "GOBP_CELLULAR_RESPONSE_TO_LIPOPROTEIN_PARTICLE_STIMULUS", "GOBP_STEROL_TRANSPORT", "GOBP_REGULATION_OF_STEROL_TRANSPORT", "GOBP_INTRACELLULAR_STEROL_TRANSPORT", "GOBP_STEROL_HOMEOSTASIS", "GOBP_LIPID_LOCALIZATION", "GOBP_REGULATION_OF_LIPID_METABOLIC_PROCESS", "GOBP_CHOLESTEROL_STORAGE",#str
  "GOBP_REGULATION_OF_INSULIN_SECRETION", "GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS",  #glucose
  "GOBP_PEPTIDE_BIOSYNTHETIC_PROCESS", "GOBP_PHOSPHORYLATION", "GOBP_GLYCOPROTEIN_BIOSYNTHETIC_PROCESS", #metabolism
  
  #CILIA
  "GOBP_CEREBROSPINAL_FLUID_CIRCULATION", "GOCC_CILIUM", "GOCC_CILIARY_TIP", "GOBP_CILIUM_MOVEMENT", "GOCC_CYTOPLASMIC_MICROTUBULE", "GOBP_MICROTUBULE_BASED_MOVEMENT", "GOMF_TUBULIN_BINDING", "GOCC_DYNEIN_COMPLEX", "GOBP_INTRACILIARY_TRANSPORT", "GOCC_INTRACILIARY_TRANSPORT_PARTICLE_B",
  
  #OTHER
  "GOBP_RESPONSE_TO_GROWTH_FACTOR", "GOBP_REGULATION_OF_RNA_SPLICING", "GOMF_RRNA_BINDING")


#per celltype
GSEA_join <- full_join(GO_labelled, GSEA_MainGenotype_Epithelial, by="pathway", suffix=c("", "_Epithelial")) %>% full_join(GSEA_MainGenotype_Fibroblast, by="pathway", suffix=c("", "_Fibroblast")) %>% full_join(GSEA_MainGenotype_Immune, by="pathway", suffix=c("", "_Immune"))

GSEA_join <- GSEA_join %>% column_to_rownames("pathway")
GSEA_join <- GSEA_join[selection,]
GSEA_join <- GSEA_join[match(selection, rownames(GSEA_join)),]
mat <- GSEA_join[,c(grep(pattern = "NES", x = colnames(GSEA_join), value = T))]
colnames(mat) <- c("Epithelial", colnames(mat)[2:ncol(mat)])
colnames(mat) <- gsub("NES_", "", colnames(mat))
mat[is.na(mat)] <- 0
mat[] <- lapply(mat, as.numeric)

mat2 <- GSEA_join[,c(grep(pattern = "padj", x = colnames(GSEA_join), value = T))]
colnames(mat2) <- c("Epithelial", colnames(mat2)[2:ncol(mat2)])
colnames(mat2) <- gsub("padj_", "", colnames(mat2))
mat2[is.na(mat2)] <- 1
mat2[] <- lapply(mat2, as.numeric)

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-3, 0, 3)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = col_fun, cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.15) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/APOE heatmap.pdf", sep = ""), font="Arial", width=w, height=h)



###Epi GO

selection <- c("GOBP_CEREBROSPINAL_FLUID_CIRCULATION", "GOCC_CILIUM", "GOBP_CILIUM_MOVEMENT", "GOBP_MICROTUBULE_BASED_MOVEMENT", "GOMF_TUBULIN_BINDING", "GOCC_DYNEIN_COMPLEX", "GOBP_INTRACILIARY_TRANSPORT",
               "GOBP_LIPID_BIOSYNTHETIC_PROCESS", "GOBP_STEROID_METABOLIC_PROCESS", "GOBP_REGULATION_OF_CHOLESTEROL_METABOLIC_PROCESS", "GOCC_LIPID_DROPLET", "GOBP_LIPID_METABOLIC_PROCESS", "GOBP_CELLULAR_LIPID_METABOLIC_PROCESS", "GOBP_CHOLESTEROL_STORAGE", "KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS",
               "KEGG_INSULIN_SIGNALING_PATHWAY", "GOBP_INTRACELLULAR_GLUCOSE_HOMEOSTASIS",
               "GOBP_INFLAMMATORY_RESPONSE", "GOBP_RESPONSE_TO_WOUNDING", "GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
               "GOBP_SECRETION", "GOBP_REGULATION_OF_PROTEIN_SECRETION", "GOBP_STEROL_TRANSPORT",
               "GOMF_CELL_ADHESION_MOLECULE_BINDING", "KEGG_ADHERENS_JUNCTION", "KEGG_TIGHT_JUNCTION", "GOCC_ANCHORING_JUNCTION", "GOCC_CELL_SUBSTRATE_JUNCTION", "GOMF_CADHERIN_BINDING")

#Ex vivo first
selection <- c("GOBP_PHOSPHOLIPID_TRANSPORT", "GOBP_LIPID_LOCALIZATION", "GOBP_REGULATION_OF_LIPID_METABOLIC_PROCESS", "GOBP_STEROID_METABOLIC_PROCESS", "GOBP_REGULATION_OF_LIPID_CATABOLIC_PROCESS", "GOBP_LIPID_BIOSYNTHETIC_PROCESS", "GOBP_REGULATION_OF_CHOLESTEROL_METABOLIC_PROCESS", 
               "GOBP_INFLAMMATORY_RESPONSE", "GOBP_ACTIVATION_OF_IMMUNE_RESPONSE", "GOBP_RESPONSE_TO_WOUNDING", "GOBP_T_CELL_DIFFERENTIATION",
               "GOMF_CELL_ADHESION_MOLECULE_BINDING", "GOCC_ANCHORING_JUNCTION", "GOCC_CELL_SUBSTRATE_JUNCTION", "GOMF_CADHERIN_BINDING",
               "GOBP_SECRETION", "GOBP_REGULATION_OF_PROTEIN_SECRETION", "GOBP_STEROL_TRANSPORT",
               "GOBP_CEREBROSPINAL_FLUID_CIRCULATION", "GOCC_CILIUM", "GOBP_CILIUM_MOVEMENT", "GOBP_MICROTUBULE_BASED_MOVEMENT", "GOMF_TUBULIN_BINDING", "GOCC_DYNEIN_COMPLEX", "GOBP_INTRACILIARY_TRANSPORT")

#per celltype
GSEA_EpiType_join <- full_join(GO_labelled, GSEA_MainGenotype_Epithelial, by="pathway", suffix=c("", "_Epithelial")) %>% full_join(GSEA_MainGenotype_Epi_1a, by="pathway", suffix=c("", "_Epi_1a")) %>% full_join(GSEA_MainGenotype_Epi_1b, by="pathway", suffix=c("", "_Epi_1b")) %>% full_join(GSEA_MainGenotype_Epi_2a, by="pathway", suffix=c("", "_Epi_2a")) %>% full_join(GSEA_MainGenotype_Epi_2b, by="pathway", suffix=c("", "_Epi_2b"))

GSEA_EpiType_join <- GSEA_EpiType_join %>% column_to_rownames("pathway")
GSEA_EpiType_join <- GSEA_EpiType_join[selection,]
GSEA_EpiType_join <- GSEA_EpiType_join[match(selection, rownames(GSEA_EpiType_join)),]
mat <- GSEA_EpiType_join[,c("NES", "NES_Epi_1a", "NES_Epi_1b", "NES_Epi_2a", "NES_Epi_2b")]
colnames(mat) <- c("Epithelial", "Epi_1a", "Epi_1b", "Epi_2a", "Epi_2b")
mat[is.na(mat)] <- 0
mat[] <- lapply(mat, as.numeric)

mat2 <- GSEA_EpiType_join[,c("padj", "padj_Epi_1a", "padj_Epi_1b", "padj_Epi_2a", "padj_Epi_2b")]
colnames(mat2) <- c("Epithelial", "Epi_1a", "Epi_1b", "Epi_2a", "Epi_2b")
mat2[is.na(mat2)] <- 1
mat2[] <- lapply(mat2, as.numeric)

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 30, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-3, 0, 3)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = col_fun, cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.1) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_EpiTypes.pdf", sep = ""), font="Arial", width=w, height=h)


pdf(paste(curr_dir, "/FGSEA/Bar_Epithelial.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat[[1]]), col=col_fun(mat[[1]]), horiz=T, xlim=c(0,2), border=NA)
dev.off()



###Fib GO

selection <- c("GOBP_CELL_MATRIX_ADHESION", "GOBP_CELL_ADHESION", "GOCC_CELL_SURFACE",
               "GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX", "GOCC_COLLAGEN_TRIMER", "KEGG_ECM_RECEPTOR_INTERACTION",
               "GOBP_REGULATION_OF_TRANSMEMBRANE_TRANSPORT", "GOBP_REGULATION_OF_POTASSIUM_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY", 
               "GOMF_IMMUNE_RECEPTOR_ACTIVITY", "GOBP_RESPONSE_TO_WOUNDING", "GOBP_INTERLEUKIN_6_PRODUCTION", "GOBP_INTERLEUKIN_8_PRODUCTION")

#Ex vivo first

selection <- c("GOBP_CELL_MATRIX_ADHESION", "GOBP_CELL_SUBSTRATE_ADHESION", "GOCC_ACTIN_CYTOSKELETON",
               "GOMF_IMMUNE_RECEPTOR_ACTIVITY", "GOBP_RESPONSE_TO_WOUNDING", "GOBP_INTERLEUKIN_6_PRODUCTION", "GOBP_LEUKOCYTE_CELL_CELL_ADHESION")


#per celltype
GSEA_FibType_join <- full_join(GO_labelled, GSEA_MainGenotype_Fibroblast, by="pathway", suffix=c("", "_Fibroblast")) %>% full_join(GSEA_MainGenotype_Fib_1, by="pathway", suffix=c("", "_Fib_1")) %>% full_join(GSEA_MainGenotype_Fib_2, by="pathway", suffix=c("", "_Fib_2")) %>% full_join(GSEA_MainGenotype_Fib_3, by="pathway", suffix=c("", "_Fib_3"))

GSEA_FibType_join <- GSEA_FibType_join %>% column_to_rownames("pathway")
GSEA_FibType_join <- GSEA_FibType_join[selection,]
GSEA_FibType_join <- GSEA_FibType_join[match(selection, rownames(GSEA_FibType_join)),]
mat <- GSEA_FibType_join[,c("NES", "NES_Fib_1", "NES_Fib_2", "NES_Fib_3")]
colnames(mat) <- c("Fibroblast", "Fib_1", "Fib_2", "Fib_3")
mat[is.na(mat)] <- 0
mat[] <- lapply(mat, as.numeric)

mat2 <- GSEA_FibType_join[,c("padj", "padj_Fib_1", "padj_Fib_2", "padj_Fib_3")]
colnames(mat2) <- c("Fibroblast", "Fib_1", "Fib_2", "Fib_3")
mat2[is.na(mat2)] <- 1
mat2[] <- lapply(mat2, as.numeric)

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 30, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-3, 0, 3)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = col_fun, cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.1) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_FibTypes.pdf", sep = ""), font="Arial", width=w, height=h)


pdf(paste(curr_dir, "/FGSEA/Bar_Fibroblast.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat[[1]]), col=col_fun(mat[[1]]), horiz=T, xlim=c(0,2), border=NA)
dev.off()


###Immune

selection <- c("GOBP_MYELOID_LEUKOCYTE_ACTIVATION", "GOBP_INTERLEUKIN_6_PRODUCTION", #Bcell, 
               "GOBP_T_CELL_PROLIFERATION", "GOBP_T_CELL_ACTIVATION", "GOBP_RESPONSE_TO_CYTOKINE",
               "GOBP_SMALL_MOLECULE_METABOLIC_PROCESS", "GOBP_LIPID_LOCALIZATION", "GOBP_CHOLESTEROL_STORAGE")
##per contrast
GSEA_BAM_join <- full_join(GO_labelled, GSEA_MainGenotype_Immune, by="pathway", suffix=c("", "_APOE")) %>% full_join(GSEA_MainGenotype_BAM_1, by="pathway", suffix=c("", "_BAM_1"))

GSEA_BAM_join <- GSEA_BAM_join %>% column_to_rownames("pathway")
GSEA_BAM_join <- GSEA_BAM_join[selection,]
GSEA_BAM_join <- GSEA_BAM_join[match(selection, rownames(GSEA_BAM_join)),]

mat <- GSEA_BAM_join[,c("NES", "NES_BAM_1")]
colnames(mat) <- c("Immune", "BAM_1")
mat[is.na(mat)] <- 0
mat[] <- lapply(mat, as.numeric)

mat2 <- GSEA_BAM_join[,c("padj", "padj_BAM_1")]
colnames(mat) <- c("Immune", "BAM_1")
mat2[is.na(mat2)] <- 1
mat2[] <- lapply(mat2, as.numeric)

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 30, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-3, 0, 3)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = col_fun, cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.1) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_Immune.pdf", sep = ""), font="Arial", width=w, height=h)


pdf(paste(curr_dir, "/FGSEA/Bar_Immune.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat[[1]]), col=col_fun(mat[[1]]), horiz=T, xlim=c(0,2), border=NA)
dev.off()


#GSVA
###Make snRNAseq modules

#Go term matrices
Temp <- cbind(GSEA_MainGenotype_Epithelial$pathway, GSEA_MainGenotype_Epithelial$leadingEdge)
sn_Epi_paths <- list("Epi_Cilia"=unique(c(unlist(Temp[Temp[,1]=="GOCC_CILIUM",2], Temp[Temp[,1]=="GOBP_MICROTUBULE_BASED_MOVEMENT",2], Temp[Temp[,1]=="GOCC_DYNEIN_COMPLEX",2]))),
                     "Epi_Secretion"=unique(c(unlist(Temp[Temp[,1]=="GOBP_REGULATION_OF_PROTEIN_SECRETION",2], Temp[Temp[,1]=="GOBP_SECRETION",2], Temp[Temp[,1]=="GOBP_STEROL_TRANSPORT",2]))),
                     "Epi_Adhesion"=unique(c(unlist(Temp[Temp[,1]=="GOCC_ADHERENS_JUNCTION",2], Temp[Temp[,1]=="GOCC_TIGHT_JUNCTION",2], Temp[Temp[,1]=="GOMF_CADHERIN_BINDING",2]))),
                     "Epi_Lipid"=unique(c(unlist(Temp[Temp[,1]=="GOBP_LIPID_BIOSYNTHETIC_PROCESS",2], Temp[Temp[,1]=="GOCC_LIPID_DROPLET",2], Temp[Temp[,1]=="GOBP_LIPID_METABOLIC_PROCESS",2]))),
                     "Epi_Immune"=unique(c(unlist(Temp[Temp[,1]=="GOBP_REGULATION_OF_T_CELL_ACTIVATION",2], Temp[Temp[,1]=="GOBP_INFLAMMATORY_RESPONSE",2], Temp[Temp[,1]=="GOBP_RESPONSE_TO_WOUNDING",2])))
)
saveRDS(sn_Epi_paths, paste0(curr_dir, "/GSVA_modules/Tsai_Epi_modules.rds"))

Temp <- cbind(GSEA_MainGenotype_Fibroblast$pathway, GSEA_MainGenotype_Fibroblast$leadingEdge)
sn_Fib_paths <- list("Fib_Barrier"=unique(c(unlist(Temp[Temp[,1]=="GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX",2], Temp[Temp[,1]=="GOBP_CELL_MATRIX_ADHESION",2]))), 
                     "Fib_transport"=unique(c(unlist(Temp[Temp[,1]=="GOMF_POTASSIUM_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",2]))), 
                     "Fib_Immune"=unique(c(unlist(Temp[Temp[,1]=="GOMF_IMMUNE_RECEPTOR_ACTIVITY",2], Temp[Temp[,1]=="GOBP_RESPONSE_TO_WOUNDING",2], Temp[Temp[,1]=="GOBP_INTERLEUKIN_6_PRODUCTION",2])))
)
saveRDS(sn_Fib_paths, paste0(curr_dir, "/GSVA_modules/Tsai_Fib_modules.rds"))

Temp <- cbind(GSEA_MainGenotype_Immune$pathway, GSEA_MainGenotype_Immune$leadingEdge)
sn_Imm_paths <- list("Imm_Immune"=unique(c(unlist(Temp[Temp[,1]=="GOBP_MYELOID_LEUKOCYTE_ACTIVATION",2], Temp[Temp[,1]=="GOBP_T_CELL_PROLIFERATION",2], Temp[Temp[,1]=="GOBP_RESPONSE_TO_CYTOKINE",2]))),
                     "Imm_Lipid"=unique(c(unlist(Temp[Temp[,1]=="GOBP_SMALL_MOLECULE_METABOLIC_PROCESS",2], Temp[Temp[,1]=="GOBP_LIPID_LOCALIZATION",2], Temp[Temp[,1]=="GOBP_CHOLESTEROL_STORAGE",2])))
)
saveRDS(sn_Imm_paths, paste0(curr_dir, "/GSVA_modules/Tsai_Imm_modules.rds"))


###Merge modules from all modalities

Exvivo_Epi_paths <- readRDS(paste0(curr_dir, "/GSVA_modules/Exvivo_Epi_modules.rds"))
tsai_Epi_paths <- readRDS(paste0(curr_dir, "/GSVA_modules/Tsai_Epi_modules.rds"))
jackson_Epi_paths <- readRDS(paste0(curr_dir, "/GSVA_modules/Jackson_Epi_modules.rds"))
Epi_modules <- StratigrapheR::merge_list(Exvivo_Epi_paths, tsai_Epi_paths, jackson_Epi_paths)
saveRDS(Epi_modules, paste0(curr_dir, "/GSVA_modules/combo_Epi_modules.rds"))

Exvivo_Fib_paths <- readRDS(paste0(curr_dir, "/GSVA_modules/Exvivo_Fib_modules.rds"))
tsai_Fib_paths <- readRDS(paste0(curr_dir, "/GSVA_modules/Tsai_Fib_modules.rds"))
jackson_Fib_paths <- readRDS(paste0(curr_dir, "/GSVA_modules/Jackson_Fib_modules.rds"))
Fib_modules <- StratigrapheR::merge_list(Exvivo_Fib_paths, tsai_Fib_paths, jackson_Fib_paths)
saveRDS(Fib_modules, paste0(curr_dir, "/GSVA_modules/combo_Fib_modules.rds"))

Exvivo_Imm_paths <- readRDS(paste0(curr_dir, "/GSVA_modules/Exvivo_Imm_modules.rds"))
tsai_Imm_paths <- readRDS(paste0(curr_dir, "/GSVA_modules/Tsai_Imm_modules.rds"))
jackson_Imm_paths <- readRDS(paste0(curr_dir, "/GSVA_modules/Jackson_Imm_modules.rds"))
Imm_modules <- StratigrapheR::merge_list(Exvivo_Imm_paths, tsai_Imm_paths, jackson_Imm_paths)
saveRDS(Imm_modules, paste0(curr_dir, "/GSVA_modules/combo_Imm_modules.rds"))


###load modules

Epi_modules <- readRDS(paste0(curr_dir, "/GSVA_modules/combo_Epi_modules.rds"))
Fib_modules <- readRDS(paste0(curr_dir, "/GSVA_modules/combo_Fib_modules.rds"))
Imm_modules <- readRDS(paste0(curr_dir, "/GSVA_modules/combo_Imm_modules.rds"))

##Dataprep

GSVA_fun(limma_MainGenotype_Epithelial, Epi_modules)
GSVA_fun(limma_MainGenotype_Fibroblast, Fib_modules)
GSVA_fun(limma_MainGenotype_Immune, Imm_modules)

GSVA_limmafun_lane(GSVA_Epithelial)
GSVA_limmafun_lane(GSVA_Fibroblast)
GSVA_limmafun_lane(GSVA_Immune)