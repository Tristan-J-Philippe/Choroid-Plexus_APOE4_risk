###
### 035 differential_RNA_InVitro_CSF.R
###
# Purpose: differential CSF proteomics
# Dependencies:
source("/030 limma_pattern_fgsea_util.R", verbose = FALSE)
curr_dir <- "/ChP_RNA" #output results to specific folder
set.seed(42)


#Set and load dependent variables for limma/DESeq functions
Gene.search <- read.csv("/DEG to GO.csv")
FDRcutoff <- 0.05
LFcutoff <- 0.25
pct_cutoff <- 0.05

####iCSF

cts <- read.csv(paste0(curr_dir, "/limma/olink iCSF protein per sample corrected and renamed.csv"), row.names=1)

#metadata
mdat <- data.frame(samples = colnames(cts), row.names = colnames(cts))
mdat <- mdat %>% mutate(Genotype = ifelse(grepl("33", samples), "33", "44")) #Extract APOE33 and APOE44 information
mdat <- mdat %>% mutate(Replicate = str_sub(samples, start=-1)) #Extract replication letter
mdat <- mdat %>% mutate(Differentiation = t(data.frame(regmatches(samples, gregexpr("[[:digit:]]+", samples)))[2,])) #Extract differentiation number

#limma
D = factor(paste(mdat$Differentiation))
R = factor(paste(mdat$Replicate))
G = factor(paste(mdat$Genotype))

Maindm <- model.matrix(~G)
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maincm <- makeContrasts(levels=Maindm, MainGenotype = G44)

heatmap(cor(cts), Colv=NA, Rowv=NA, symm=T, main="Sample correlations")
vdata <- cts

corfit <- duplicateCorrelation(cts, block=R)

Mainf <- lmFit(cts, Maindm, weights = cts, block=R, correlation=corfit$consensus.correlation)
Mainf <- contrasts.fit(Mainf, Maincm)
Mainf <- eBayes(Mainf, trend = TRUE)

limma_resfun(coef="MainGenotype", celltype="iCSF")


#Prep CSF proteomics
Path_fun(limma_MainGenotype_iCSF)


###CSF-DEP

selection <- c("GOCC_CELL_CELL_JUNCTION", "GOMF_CADHERIN_BINDING", "GOCC_ANCHORING_JUNCTION", "GOBP_ACTIN_FILAMENT_BASED_PROCESS",
               "GOBP_VESICLE_MEDIATED_TRANSPORT", "GOBP_MITOCHONDRIAL_MEMBRANE_ORGANIZATION")


GSEA_celltype_join <- full_join(GO_labelled, GSEA_MainGenotype_iCSF, by="pathway", suffix=c("", "_iCSF"))

GSEA_celltype_join <- GSEA_celltype_join[!is.na(GSEA_celltype_join$pathway),]
GSEA_celltype_join <- GSEA_celltype_join %>% rownames_to_column("XX")
GSEA_celltype_join <- GSEA_celltype_join %>% column_to_rownames("pathway")
GSEA_celltype_join <- GSEA_celltype_join[selection,]
GSEA_celltype_join <- GSEA_celltype_join[match(selection, rownames(GSEA_celltype_join)),]

mat <- as.matrix(GSEA_celltype_join[,c("NES")])
rownames(mat) <- rownames(GSEA_celltype_join)
mat[is.na(mat)] <- 0

mat2 <- as.matrix(GSEA_celltype_join[,c("padj")])
rownames(mat2) <- rownames(GSEA_celltype_join)
mat2[is.na(mat2)] <- 1


mat3 <- as.matrix(GSEA_celltype_join[,c("pval")])
rownames(mat3) <- rownames(GSEA_celltype_join)
mat3[is.na(mat3)] <- 1

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 30, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7), at = c(-3, 0, 3)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = col_fun, cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.15) {
    grid.text("*", x, y, gp=gpar(fontsize=7), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_iCSF.pdf", sep = ""), font="Arial", width=w, height=h)


pdf(paste(curr_dir, "/FGSEA/Bar_iCSF.pdf", sep = ""), width=3, height=nrow(mat)*.5)
barplot(abs(mat[,1]), col=col_fun(mat[,1]), horiz=T, xlim=c(0,2), border=NA)
dev.off()

###DEP

cm <- make_comb_mat(list_to_matrix(list(Epithelial = c(MainGenotype_Epithelial_sig$gene), Fibrobalst = c(MainGenotype_Fibroblast_sig$gene), DEP_iCSF = c(iCSF_DEP_corrected$protein))))
plts <- UpSet(cm, set_order = rownames(cm), pt_size = unit(2, "mm"), top_annotation = upset_top_annotation(cm, height = unit(2.5, "inch"), annotation_name_rot = 90))
plts
pdf(file=paste(curr_dir, "/limma/upset_DEP_AD.pdf", sep=""), width=4, height=12)
plts
dev.off()

