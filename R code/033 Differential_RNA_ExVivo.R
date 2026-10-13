###
### 033 Differential_RNA_InVitro.R
###
# Purpose: Differential in single nuclei RNA dataset
# Dependencies:
source("/030 limma_pattern_fgsea_util.R", verbose = FALSE)
curr_dir <- "/ChP_RNA" #output results to specific folder
set.seed(42)
set.seed(42)

#Set and load dependent variables for limma/DESeq functions
Gene.search <- read.csv("/DEG to GO.csv")
FDRcutoff <- 0.05
LFcutoff <- 0.25
pct_cutoff <- 0.05

#metadata
mdat <- read_csv(paste0(curr_dir, "/metadata/ex vivo ChP metadata (analytical) newAPOE RNAonly.csv"))
mdat <- column_to_rownames(mdat, "projid")
##Subset
###Final cell type

Ex_ChP <- readRDS(paste0(curr_dir, "/CellType/Finalcelltype/Ex_ChP.RDS"))
#Subcelltype
datprepfun_lane(Ex_ChP, "Epi_1a")
datprepfun_lane(Ex_ChP, "Epi_1b")
datprepfun_lane(Ex_ChP, "Epi_2a")
datprepfun_lane(Ex_ChP, "Epi_2b")
datprepfun_lane(Ex_ChP, "Fib_1")
datprepfun_lane(Ex_ChP, "Fib_2")
datprepfun_lane(Ex_ChP, "Fib_3")
datprepfun_lane(Ex_ChP, "Endo_1")
datprepfun_lane(Ex_ChP, "Endo_3")
datprepfun_lane(Ex_ChP, "Endo_2")
datprepfun_lane(Ex_ChP, "Mural_1")
datprepfun_lane(Ex_ChP, "Mural_2")
datprepfun_lane(Ex_ChP, "Mural_3")
datprepfun_lane(Ex_ChP, "BAM_1")
datprepfun_lane(Ex_ChP, "BAM_2")
datprepfun_lane(Ex_ChP, "BAM_3")
datprepfun_lane(Ex_ChP, "T_Cell")
#datprepfun_lane(Ex_ChP, "BAM")

datprepfun_lane(Ex_ChP, "Neuro")
datprepfun_lane(Ex_ChP, "Astro")
datprepfun_lane(Ex_ChP, "Oligo")


###Major cell type

#!!!!Change to Majorcelltype in datprepfun_lane
datprepfun_lane(Ex_ChP, "Epithelial")
datprepfun_lane(Ex_ChP, "Fibroblast")
datprepfun_lane(Ex_ChP, "Endothelial")
datprepfun_lane(Ex_ChP, "Mural")
datprepfun_lane(Ex_ChP, "Immune")

datprepfun_lane(Ex_ChP, "Parenchyma")

pcts <- cbind(pct_Epi_1a, pct_Epi_1b, pct_Epi_2a, pct_Epi_2b, pct_Fib_1, pct_Fib_2, pct_Fib_3, pct_Endo_1, pct_Endo_2, pct_Endo_3, pct_Mural_1, pct_Mural_2, pct_Mural_3, pct_BAM_1, pct_BAM_2, pct_BAM_3, pct_T_Cell, pct_Epithelial, pct_Fibroblast, pct_Endothelial, pct_Mural, pct_Immune)
colnames(pcts) <- c("Epi_1a", "Epi_1b", "Epi_2a", "Epi_2b", "Fib_1", "Fib_2", "Fib_3", "Endo_1", "Endo_2", "Endo_3", "Mural_1", "Mural_2", "Mural_3", "BAM_1", "BAM_2", "BAM_3", "T_Cell", "Epithelial", "Fibroblast", "Endothelial", "Mural", "Immune")
write.csv(pcts, paste0(curr_dir, "/Percent_genes_in_cells.csv"))

cts_names <- ls(pattern = "^cts_")

for (name in cts_names){
  outfile <- sub("^cts_", "", name)
  write.csv(get(name), file = paste0(curr_dir, "/limma_nolane/cnt_mtx/", outfile, ".csv"))
}


##limma runs
###factors set up

#discrete var
msex = factor(paste(mdat$msex))
X10X_batch = factor(paste(mdat$X10X_batch))
library_batch = factor(paste(mdat$library_batch))
apoe = factor(paste(mdat$apoe_genotype))
apoe_sex = factor(paste(mdat$apoe_genotype, mdat$msex, sep="."))
cogdx = factor(paste(mdat$cogdx))
Cog_Path = factor(paste(mdat$Cog_Path))
cogdx_simple = factor(paste(mdat$cogdx_simple), levels=c("ADd","MCI","NCI"))
cogdx_sex = factor(paste(mdat$cogdx, mdat$msex, sep="."))
cogdx_ad_adnc = factor(paste(mdat$cogdx_simple, mdat$ad_adnc, sep="."))
apoe_cogdx <- factor(paste(mdat$apoe_genotype, mdat$cogdx_simple, sep="."))
apoe_adnc <- factor(paste(mdat$apoe_genotype, mdat$ad_adnc, sep="."))

ad_adnc <- factor(paste(mdat$ad_adnc))

#continuous
cogdx_c = as.numeric(mdat$cogdx)
projid = as.numeric(mdat$projid)
age = as.numeric(mdat$age_death)
pmi = as.numeric(mdat$pmi)
educ = as.numeric(mdat$educ)
gpath = as.numeric(mdat$gpath)
Braaksc = as.numeric(mdat$braaksc)
tdp = as.numeric(mdat$tdp_st4)
tangles = as.numeric(mdat$tangsqrt_est_8reg)
amyl = as.numeric(mdat$amylsqrt_est_8reg)

arteriol = as.numeric(mdat$arteriol_scler)
caa = as.numeric(mdat$caa_4gp)
cvda = as.numeric(mdat$cvda_4gp2)
cogn_global = as.numeric(mdat$cogn_global_last)
sXrho = as.numeric(mdat$`soupX Rho`)

apoe <- factor(mdat$apoe_genotype, levels=c("33","34","other"))
cogdx_simple <- factor(mdat$cogdx_simple, levels=c("NCI","MCI","ADd"))
options(na.action='na.pass')
Maindm <- model.matrix(~apoe*cogdx_simple + msex + age + X10X_batch + pmi)
colnames(Maindm) <- make.names(colnames(Maindm))
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, cogdxAPOE_34v33=apoe34, cogdxAPOE_ADvNCI=cogdx_simpleADd, cogdxAPOE_cogdxapoe=apoe34.cogdx_simpleADd, cogdxAPOE_ADvNCIin34=cogdx_simpleADd+apoe34.cogdx_simpleADd)

limmafun_lane(cts_Epithelial, pct_Epithelial)
limmafun_lane(cts_Fibroblast, pct_Fibroblast)
limmafun_lane(cts_Endothelial, pct_Endothelial)
limmafun_lane(cts_BAM, pct_BAM)

limmafun_lane(cts_Epi_1a, pct_Epi_1a)
limmafun_lane(cts_Epi_1b, pct_Epi_1b)
limmafun_lane(cts_Epi_2a, pct_Epi_2a)
limmafun_lane(cts_Epi_2b, pct_Epi_2b)
limmafun_lane(cts_Fib_1, pct_Fib_1)
limmafun_lane(cts_Fib_2, pct_Fib_2)
limmafun_lane(cts_Fib_3, pct_Fib_3)
limmafun_lane(cts_BAM_1, pct_BAM_1)

#!!!!No batch correction
options(na.action='na.pass')
Maindm <- model.matrix(~apoe*cogdx_simple + msex + age + pmi)
colnames(Maindm) <- make.names(colnames(Maindm))
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, cogdxAPOE_34v33=apoe34, cogdxAPOE_ADvNCI=cogdx_simpleADd, cogdxAPOE_cogdxapoe=apoe34.cogdx_simpleADd, cogdxAPOE_ADvNCIin34=cogdx_simpleADd+apoe34.cogdx_simpleADd)
limmafun_lane(cts_Mural, pct_Mural)
limmafun_lane(cts_BAM_2, pct_BAM_2)

#FGSEA
####APOE cogdx ctrl (+cogdx+cogdx*APOE)

Path_fun(limma_cogdxAPOE_34v33_Epi_1a)
Path_fun(limma_cogdxAPOE_34v33_Epi_1b)
Path_fun(limma_cogdxAPOE_34v33_Epi_2a)
Path_fun(limma_cogdxAPOE_34v33_Epi_2b)
Path_fun(limma_cogdxAPOE_34v33_Fib_1)
Path_fun(limma_cogdxAPOE_34v33_Fib_2)
Path_fun(limma_cogdxAPOE_34v33_Fib_3)
Path_fun(limma_cogdxAPOE_34v33_BAM_1)
Path_fun(limma_cogdxAPOE_34v33_BAM_2)

Path_fun(limma_cogdxAPOE_34v33_Epithelial)
Path_fun(limma_cogdxAPOE_34v33_Fibroblast)
Path_fun(limma_cogdxAPOE_34v33_Endothelial)
Path_fun(limma_cogdxAPOE_34v33_Mural)
Path_fun(limma_cogdxAPOE_34v33_BAM)

####APOE cogdx ctrl (+cogdx+cogdx*APOE)

Path_fun(limma_cogdxAPOE_ADvNCI_Epi_1a)
Path_fun(limma_cogdxAPOE_ADvNCI_Epi_1b)
Path_fun(limma_cogdxAPOE_ADvNCI_Epi_2a)
Path_fun(limma_cogdxAPOE_ADvNCI_Epi_2b)
Path_fun(limma_cogdxAPOE_ADvNCI_Fib_1)
Path_fun(limma_cogdxAPOE_ADvNCI_Fib_2)
Path_fun(limma_cogdxAPOE_ADvNCI_Fib_3)
Path_fun(limma_cogdxAPOE_ADvNCI_BAM_1)
Path_fun(limma_cogdxAPOE_ADvNCI_BAM_2)

Path_fun(limma_cogdxAPOE_ADvNCI_Epithelial)
Path_fun(limma_cogdxAPOE_ADvNCI_Fibroblast)
Path_fun(limma_cogdxAPOE_ADvNCI_Endothelial)
Path_fun(limma_cogdxAPOE_ADvNCI_Mural)
Path_fun(limma_cogdxAPOE_ADvNCI_BAM)


####APOE_by_cogdx ctrl (+cogdx+cogdx*APOE)

Path_fun(limma_cogdxAPOE_ADvNCIin34_Epi_1a)
Path_fun(limma_cogdxAPOE_ADvNCIin34_Epi_1b)
Path_fun(limma_cogdxAPOE_ADvNCIin34_Epi_2a)
Path_fun(limma_cogdxAPOE_ADvNCIin34_Epi_2b)
Path_fun(limma_cogdxAPOE_ADvNCIin34_Fib_1)
Path_fun(limma_cogdxAPOE_ADvNCIin34_Fib_2)
Path_fun(limma_cogdxAPOE_ADvNCIin34_Fib_3)
Path_fun(limma_cogdxAPOE_ADvNCIin34_BAM_1)
Path_fun(limma_cogdxAPOE_ADvNCIin34_BAM_2)

Path_fun(limma_cogdxAPOE_ADvNCIin34_Epithelial)
Path_fun(limma_cogdxAPOE_ADvNCIin34_Fibroblast)
Path_fun(limma_cogdxAPOE_ADvNCIin34_Endothelial)
Path_fun(limma_cogdxAPOE_ADvNCIin34_Mural)
Path_fun(limma_cogdxAPOE_ADvNCIin34_BAM)


####APOE_by_cogdx ctrl (+cogdx+cogdx*APOE)

Path_fun(limma_cogdxAPOE_cogdxapoe_Epi_1a)
Path_fun(limma_cogdxAPOE_cogdxapoe_Epi_1b)
Path_fun(limma_cogdxAPOE_cogdxapoe_Epi_2a)
Path_fun(limma_cogdxAPOE_cogdxapoe_Epi_2b)
Path_fun(limma_cogdxAPOE_cogdxapoe_Fib_1)
Path_fun(limma_cogdxAPOE_cogdxapoe_Fib_2)
Path_fun(limma_cogdxAPOE_cogdxapoe_Fib_3)
Path_fun(limma_cogdxAPOE_cogdxapoe_BAM_1)
Path_fun(limma_cogdxAPOE_cogdxapoe_BAM_2)

Path_fun(limma_cogdxAPOE_cogdxapoe_Epithelial)
Path_fun(limma_cogdxAPOE_cogdxapoe_Fibroblast)
Path_fun(limma_cogdxAPOE_cogdxapoe_Endothelial)
Path_fun(limma_cogdxAPOE_cogdxapoe_Mural)
Path_fun(limma_cogdxAPOE_cogdxapoe_BAM)



#Figures DEGs

library(ggplot2)
library(ggrepel)
FDRcutoff <- 0.01
LFcutoff <- 0.25
Epi_genes <- c("TUBB4B", "CFAP43", "DNAH12", "FOXJ1", "CFAP47", "DYNC1I1", "CFAP299", "DGAT1", "LRP10", "P2RY14", "KCNIP1", "APOLD1", "SERPINA3", "IL6R", "CHI3L1", "UBB", "PLTP", "NQO1", "KCNN2", "CPE", "LAMB1", "CLDN5", "COL9A3")
Fib_genes <- c("COL8A1", "COL15A1", "LAMA2", "TJP1", "COL8A1", "COL1A2", "ADAMTS9", "CD96", "APOE", "ACSL4", "APOLD1", "KCNMA1", "SLC38A2", "TFRC")
BAM_genes <- c("HSP90B1", "HSPA6", "HSPA1B", "HSP90AA1")

ChP_proteins <- c("CLN6", "CFAP58", "CEP43", "ARL2BP", "CEP164", "CFAP58",
                  "APOA4", "HSD11B1", "CHAC2", "MAN2B2", "FABP5", "ALDH9A1", "IP6K1", "KDSR", "SLC37A4", "PTGR1", "A2M", "DDO",
                  "SLC7A8", "UAP1L1", "STX2", "IPO4", "RGL1", "ITPA", "TMEM63A", "CHMP1A",
                  "IL1RAP", "HSD11B1", "PP1R1A", "GABRA3", "PARP14", "GRIA2", "HLA-B", "VRK2", "IFIH1", "SERPINF1", "ORAI3", "DNAJB1",
                  "APC", "MOB4", "TMEM98", "GATAD2A", "B3GALNT2", "EGFL8", "ICAM3")
ChP_proteins[ChP_proteins %in% limma_cogdx_ADvNCI_Fibroblast$gene]
ChP_proteins[ChP_proteins %in% limma_cogdx_ADvNCI_Epithelial_sig$gen]

CSF_proteins <- c("BCAN", "CSH1", "SCG3", "NGF", "TGM2", "RPE", "RELN", "CDH3", "CPE", "PTK7", "VASH1", "EPGN", "ARHGAP5", "AMOT",
                  "HLA-E", "MIF", "IL20RA", "B2M", "PTGES2", "IL13RA2")
CSF_proteins[CSF_proteins %in% limma_cogdx_ADvNCI_Fibroblast$gene]
CSF_proteins[CSF_proteins %in% limma_cogdx_ADvNCI_Epithelial_sig$gene]


##custom Volcano

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

Volcano_fun(limma_cogdxAPOE_34v33NCI_Epithelial, Epi_genes)
Volcano_fun(limma_cogdxAPOE_34v33NCI_Fibroblast, Fib_genes)
Volcano_fun(limma_cogdxAPOE_34v33NCI_BAM, BAM_genes)

Volcano_fun(limma_cogdxAPOE_ADvNCI_Epithelial, Epi_genes)
Volcano_fun(limma_cogdxAPOE_ADvNCI_Fibroblast, Fib_genes)
Volcano_fun(limma_cogdxAPOE_ADvNCI_BAM, BAM_genes)

Volcano_fun <- function(TT_res, genes){
  name <- deparse(substitute(TT_res))
  name <- str_sub(name, 7)
  df.1b <- TT_res
  df.1b <- df.1b[!duplicated(df.1b$gene),]
  rownames(df.1b) <- NULL
  df.1b <- df.1b %>% column_to_rownames("gene")
  df.1b$diffexpressed <- "NS"
  df.1b$diffexpressed[df.1b$logFC > LFcutoff & df.1b$P.Value < 0.05] <- "UP"
  df.1b$diffexpressed[df.1b$logFC < -LFcutoff & df.1b$P.Value < 0.05] <- "DOWN"
  mycolors <- c("blue", "red", "black")
  names(mycolors) <- c("DOWN", "UP", "NS")
  
  lbl_gene <- df.1b[rownames(df.1b) %in% genes,]
  
  plts <- ggplot(data=df.1b, aes(x=logFC, y=-log10(P.Value), col = diffexpressed)) +
    ggrastr::rasterise(geom_point(size=.1), dpi=600) + theme_bw() + xlim(-2,2) + ylim(0,8) + 
    scale_color_manual(values=mycolors) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.position = "none", axis.text = element_text(size = 7), axis.title=element_text(size=7)) +
    geom_vline(xintercept=c(-LFcutoff, LFcutoff), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red") + geom_text_repel(data=lbl_gene, mapping=aes(label=rownames(lbl_gene)), size=2, max.overlaps=1000)
  plot(plts)
  graph2pdf(plts, paste(curr_dir, "/limma/volcano_", name, ".pdf", sep=""), width=1.574, height=1.574)
}

Volcano_fun(limma_cogdx_ADvNCI_ChP, ChP_proteins)
Volcano_fun(limma_cogdx_ADvMCI_ChP, ChP_proteins)
Volcano_fun(limma_cogdx_MCIvNCI_ChP, ChP_proteins)

Volcano_fun(limma_cogdx_ADvNCI_CSF, CSF_proteins)
Volcano_fun(limma_cogdx_ADvMCI_CSF, CSF_proteins)
Volcano_fun(limma_cogdx_MCIvNCI_CSF, CSF_proteins)


##Custom scatterplot

tcutoff=3.2
#t scatterplot function requires 2 toptable dfs and genes to label
tScat_fun <- function(datA, datB, genes=NULL){
  nameA <- deparse(substitute(datA))
  nameA <- str_sub(nameA, 7)
  nameB <- deparse(substitute(datB))
  nameB <- str_sub(nameB, 7)
  
  cutoff <- tcutoff
  
  #Genes in one in other
  datA <- datA[!duplicated(datA$gene),]
  datB <- datB[!duplicated(datB$gene),]
  datA <- data.frame(datA$gene, datA$t)
  datB <- data.frame(datB$gene, datB$t)
  colnames(datA) <- c("gene", "datA")
  colnames(datB) <- c("gene", "datB")
  
  #Sorting and log transforming
  val <- merge(datA, datB, by="gene", all.x=T, all.y=T)
  val[is.na(val)] <- runif(sum(is.na(val)), min = -cutoff, max = cutoff)
  val <- column_to_rownames(val, "gene")
  
  #Determining which are significant in both or either
  val <- mutate(val, AandB = ifelse(val$datA >= cutoff & val$datB >= cutoff, 1, 0) + ifelse(val$datA <= -cutoff & val$datB <= -cutoff, 1, 0))
  val <- mutate(val, AorB = ifelse(val$datA >= cutoff | val$datB >= cutoff, 1, 0) + ifelse(val$datA <= -cutoff | val$datB <= -cutoff, 1, 0))
  val <- mutate(val, Sums = rowSums(val[, c(3,4)]))
  
  lbl_gene <- val[rownames(val) %in% genes,]
  
  #plot for use in illustrator
  plts <- ggplot(val, aes(x = datA, y = datB)) +
    xlim(-8,8) + ylim(-8,8) +
    xlab(nameA) + ylab(nameB) +
    geom_point(aes(colour = cut(Sums, c(-Inf, 0, 1, 1.9, Inf))), show.legend = FALSE, size = .1) +
    scale_colour_manual(name = "Sums",
                        values = c("(-Inf,0]" = "black",
                                   "(0,1]" = "darkred",
                                   "(1.9, Inf]" = "red")) +
    coord_fixed(ratio = 1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.position = "none", axis.text = element_text(size = 7), axis.title=element_text(size=7)) +
    geom_vline(xintercept=c(-tcutoff, tcutoff), col="red") +
    geom_hline(yintercept=c(-tcutoff, tcutoff), col="red") +
    geom_text_repel(data=lbl_gene, mapping=aes(label=rownames(lbl_gene)), size=2, seed = 42, box.padding = 0.5, max.overlaps = Inf, arrow = arrow(length = unit(0.010, "npc")))
  plot(plts)
  graph2pdf(plts, paste(curr_dir, "/limma/scatterplot_", nameA, "_", nameB, "tval.pdf", sep=""))
  
  assign(paste(nameA, nameB, sep="_"), val, envir = .GlobalEnv)
}
#Epithelial
tScat_fun(limma_cogdx_ADvNCI_Epithelial, limma_cogdx_ADvMCI_Epithelial, Epi_genes)
tScat_fun(limma_cogdx_ADvNCI_Epithelial, limma_cogdx_MCIvNCI_Epithelial, Epi_genes)
tScat_fun(limma_cogdx_MCIvNCI_Epithelial, limma_cogdx_ADvMCI_Epithelial, Epi_genes)

#Epi sub
tScat_fun(limma_cogdx_ADvNCI_Epi_1a, limma_cogdx_ADvMCI_Epi_1a, Epi_genes)
tScat_fun(limma_cogdx_ADvNCI_Epi_1a, limma_cogdx_MCIvNCI_Epi_1a, Epi_genes)
tScat_fun(limma_cogdx_MCIvNCI_Epi_1a, limma_cogdx_ADvMCI_Epi_1a, Epi_genes)

tScat_fun(limma_cogdx_ADvNCI_Epi_1b, limma_cogdx_ADvMCI_Epi_1b, Epi_genes)
tScat_fun(limma_cogdx_ADvNCI_Epi_1b, limma_cogdx_MCIvNCI_Epi_1b, Epi_genes)
tScat_fun(limma_cogdx_MCIvNCI_Epi_1b, limma_cogdx_ADvMCI_Epi_1b, Epi_genes)

tScat_fun(limma_cogdx_ADvNCI_Epi_2a, limma_cogdx_ADvMCI_Epi_2a, Epi_genes)
tScat_fun(limma_cogdx_ADvNCI_Epi_2a, limma_cogdx_MCIvNCI_Epi_2a, Epi_genes)
tScat_fun(limma_cogdx_MCIvNCI_Epi_2a, limma_cogdx_ADvMCI_Epi_2a, Epi_genes)

tScat_fun(limma_cogdx_ADvNCI_Epi_2b, limma_cogdx_ADvMCI_Epi_2b, Epi_genes)
tScat_fun(limma_cogdx_ADvNCI_Epi_2b, limma_cogdx_MCIvNCI_Epi_2b, Epi_genes)
tScat_fun(limma_cogdx_MCIvNCI_Epi_2b, limma_cogdx_ADvMCI_Epi_2b, Epi_genes)

#Fibroblast
tScat_fun(limma_cogdx_ADvNCI_Fibroblast, limma_cogdx_ADvMCI_Fibroblast, unique(Fib_genes))
tScat_fun(limma_cogdx_ADvNCI_Fibroblast, limma_cogdx_MCIvNCI_Fibroblast, unique(Fib_genes))
tScat_fun(limma_cogdx_MCIvNCI_Fibroblast, limma_cogdx_ADvMCI_Fibroblast, unique(Fib_genes))

#Fibro sub
tScat_fun(limma_cogdx_MCIvNCI_Fib_1, limma_cogdx_ADvMCI_Fib_1, unique(Fib_genes))
tScat_fun(limma_cogdx_MCIvNCI_Fib_2, limma_cogdx_ADvMCI_Fib_2, unique(Fib_genes))
tScat_fun(limma_cogdx_MCIvNCI_Fib_3, limma_cogdx_ADvMCI_Fib_3, unique(Fib_genes))

tScat_fun(limma_cogdx_ADvNCI_Fib_1, limma_cogdx_MCIvNCI_Fib_1, unique(Fib_genes))
tScat_fun(limma_cogdx_ADvNCI_Fib_2, limma_cogdx_MCIvNCI_Fib_2, unique(Fib_genes))
tScat_fun(limma_cogdx_ADvNCI_Fib_3, limma_cogdx_MCIvNCI_Fib_3, unique(Fib_genes))

tScat_fun(limma_cogdx_ADvNCI_Fib_1, limma_cogdx_ADvMCI_Fib_1, unique(Fib_genes))
tScat_fun(limma_cogdx_ADvNCI_Fib_2, limma_cogdx_ADvMCI_Fib_2, unique(Fib_genes))
tScat_fun(limma_cogdx_ADvNCI_Fib_3, limma_cogdx_ADvMCI_Fib_3, unique(Fib_genes))

tScat_fun(limma_cogdx_ADvNCI_Fib_1, limma_cogdx_ADvNCI_Fib_2, unique(Fib_genes))
tScat_fun(limma_cogdx_ADvNCI_Fib_3, limma_cogdx_ADvNCI_Fib_2, unique(Fib_genes))

#Proteomics CHP
tcutoff=2.064
tScat_fun(limma_cogdx_ADvNCI_Fibroblast, limma_cogdx_ADvNCI_ChP, unique(ChP_proteins))
tScat_fun(limma_cogdx_ADvMCI_Fibroblast, limma_cogdx_ADvMCI_ChP, unique(ChP_proteins))
tScat_fun(limma_cogdx_MCIvNCI_Fibroblast, limma_cogdx_MCIvNCI_ChP, unique(ChP_proteins))
tScat_fun(limma_cogdx_ADvNCI_Epithelial, limma_cogdx_ADvNCI_ChP, unique(ChP_proteins))
tScat_fun(limma_cogdx_ADvMCI_Epithelial, limma_cogdx_ADvMCI_ChP, unique(ChP_proteins))
tScat_fun(limma_cogdx_MCIvNCI_Epithelial, limma_cogdx_MCIvNCI_ChP, unique(ChP_proteins))

#Proteomics CSF
tScat_fun(limma_cogdx_ADvNCI_Fibroblast, limma_cogdx_ADvNCI_CSF, unique(CSF_proteins))
tScat_fun(limma_cogdx_ADvMCI_Fibroblast, limma_cogdx_ADvMCI_CSF, unique(CSF_proteins))
tScat_fun(limma_cogdx_MCIvNCI_Fibroblast, limma_cogdx_MCIvNCI_CSF, unique(CSF_proteins))
tScat_fun(limma_cogdx_ADvNCI_Epithelial, limma_cogdx_ADvNCI_CSF, unique(CSF_proteins))
tScat_fun(limma_cogdx_ADvMCI_Epithelial, limma_cogdx_ADvMCI_CSF, unique(CSF_proteins))
tScat_fun(limma_cogdx_MCIvNCI_Epithelial, limma_cogdx_MCIvNCI_CSF, unique(CSF_proteins))

#Proteomics CHPvCSF
tScat_fun(limma_cogdx_ADvNCI_ChP, limma_cogdx_ADvNCI_CSF, unique(CSF_proteins))
tScat_fun(limma_cogdx_ADvMCI_ChP, limma_cogdx_ADvMCI_CSF, unique(CSF_proteins))
tScat_fun(limma_cogdx_MCIvNCI_ChP, limma_cogdx_MCIvNCI_CSF, unique(CSF_proteins))

#Proteomics cogdx
tScat_fun(limma_cogdx_ADvNCI_ChP, limma_cogdx_MCIvNCI_ChP, unique(CSF_proteins))
tScat_fun(limma_cogdx_ADvNCI_ChP, limma_cogdx_ADvMCI_ChP, unique(CSF_proteins))
tScat_fun(limma_cogdx_ADvMCI_ChP, limma_cogdx_MCIvNCI_ChP, unique(CSF_proteins))

tScat_fun(limma_cogdx_ADvNCI_CSF, limma_cogdx_MCIvNCI_CSF, unique(CSF_proteins))
tScat_fun(limma_cogdx_ADvNCI_CSF, limma_cogdx_ADvMCI_CSF, unique(CSF_proteins))
tScat_fun(limma_cogdx_ADvMCI_CSF, limma_cogdx_MCIvNCI_CSF, unique(CSF_proteins))



##Boxplots

box_fun <- function(limma_res, metadat, colour){
  #Takes a column of gene expression values and merges it with metadata from mdat in the environment to produce plot and return a box + point plot.
  #limma result (will trim the first 8 and last 2 cols since those are usually stats)
  #metadat: Column from mdat to aggregate data
  #Set of colours to use for box plot
  #Select relevant rows from Design Matrix 
  df <- data.frame(t(limma_res[limma_res$gene==gene, c(8:(length(colnames(limma_res))-2))])) #!!!!Trims the first 8 and last 2 cols since those are stats
  colnames(df) <- "Expression"
  dm <- matrix(ncol=ncol(mdat))
  dm <- dm[-1,]
  colnames(dm) <- colnames(mdat)
  
  for(i in rownames(df)){
    base_name <- sapply(str_split(i,"_",),'[',1) #get projid names
    base_name <- str_remove(base_name, "X")
    if(base_name %in% rownames(mdat)){
      name_rows <- rownames(dm)
      row_add <- mdat[base_name,]
      dm <- rbind(dm, row_add)
      rownames(dm) <- c(name_rows, paste(i))
    }else{print(paste("Warning colname of data:", i, "not in rowname of Design Matrix"))}
  }
  
  df <- merge(df, dm, by='row.names')
  df <- data.frame(df$Expression, df[[metadat]])
  colnames(df) <- c("Expression", "Group")
  
  df <- df[order(df$Group, decreasing=T),]
  df$Group <- factor(df$Group, levels=unique(df$Group))
  plts <- ggplot(df, aes(x=Group, y=Expression, fill=Group)) +
    geom_boxplot(fatten=1, colour="black") +
    scale_fill_manual(values=colour) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), legend.position="none", axis.line = element_line(colour = "black"), text=element_text(size=7, color="black"),
          axis.title.x=element_blank(), axis.title.y=element_text(face="italic", size=7)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
    ylab(gene)
  plot(plts)
  return(plts)
}

#Set up to rename data by disease
#Epithelial
Epi_genes <- Epi_genes[Epi_genes %in% limma_cogdx_ADvNCI_Epithelial$gene]
plt <- NULL
for(gene in Epi_genes){
  plts <- box_fun(limma_cogdx_ADvNCI_Epithelial, "Cog_Path", c("skyblue", "#F86247", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Epithelial_All.pdf", sep=""), plt, width=1, height=.5*length(Epi_genes))

gene_list <- c("COL9A1", "EFNA5", "ITGAV")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Epithelial$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Epithelial, "Cog_Path", c("skyblue", "#F86247", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Epithelial_cellchat.pdf", sep=""), plt, width=1, height=.5*length(gene_list))

gene_list <- c("CEP164", "APOA4", "ALDH9A1", "ALDH2", "PTPRM", "SLC44A4", "SLC7A8", "SLC37A4", "CACNB1", "IL20RA", "IL13RA2", "IL1RAP", "HLA-B", "DNAJB1")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Epithelial$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Epithelial, "Cog_Path", c("skyblue", "#F86247", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Epithelial_proteomic.pdf", sep=""), plt, width=1, height=.5*length(gene_list))


#Fibroblast
gene_list <- Fib_genes
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Fibroblast$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Fibroblast, "Cog_Path", c("skyblue", "#F86247", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Fibroblast_All.pdf", sep=""), plt, width=1, height=.5*length(gene_list))

gene_list <- c("LAMA2", "COL4A2", "COL1A2", "COL4A4", "SEMA3C", "BMP5", "JAM3", "CLDN11", "NAMPT")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Fibroblast$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Fibroblast, "Cog_Path", c("skyblue", "#F86247", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Fibroblast_Cellchat.pdf", sep=""), plt, width=1, height=.5*length(gene_list))

gene_list <- c("COL15A1", "COL18A1", "COL27A1", "COL4A1", "COL4A2", "COL4A4", "COL5A3", "COL7A1", "COL8A1")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Fibroblast$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Fibroblast, "Cog_Path", c("skyblue","#F86247",  "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Fibroblast_Collagen.pdf", sep=""), plt, width=.7, height=.5*length(Fib_genes))

gene_list <- c("CEP164", "APOA4", "ALDH9A1", "ALDH2", "PTPRM", "SLC44A4", "SLC7A8", "SLC37A4", "CACNB1", "IL20RA", "IL13RA2", "IL1RAP", "HLA-B", "DNAJB1")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Fibroblast$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Fibroblast, "Cog_Path", c("skyblue", "#F86247", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Fibroblast_Proteomic.pdf", sep=""), plt, width=1, height=.5*length(gene_list))

#Immune
gene_list <- c("HSPA1A", "HSP90B1", "HSPA6", "HSPA1B", "DNAJA1", "HSP90AA1")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_BAM$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_BAM, "Cog_Path", c("skyblue", "#F86247", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Immune.pdf", sep=""), plt, width=1, height=.5*length(gene_list))

gene_list <- c("CD68", "HLA-DRA", "HLA-DRB", "HLA-DQA1", "HLA-E")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_BAM$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_BAM, "Cog_Path", c("skyblue", "#F86247", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Immune_inflammation.pdf", sep=""), plt, width=.7, height=.5*length(gene_list))

gene_list <- c("CEP164", "APOA4", "ALDH9A1", "ALDH2", "PTPRM", "SLC44A4", "SLC7A8", "SLC37A4", "CACNB1", "IL20RA", "IL13RA2", "IL1RAP", "HLA-B", "DNAJB1")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_BAM$gene]
plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_BAM, "Cog_Path", c("skyblue", "#F86247", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Immune_Proteomic.pdf", sep=""), plt, width=1, height=.5*length(gene_list))

#Endothelial
gene_list <- c("LAMA5", "EFNB2", "PECAM1", "ITGA1", "ITGAV")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Endo$gene]

plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Endo, "Cog_Path", c("skyblue", "#F86247", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Endothelial_cellchat.pdf", sep=""), plt, width=1, height=.5*length(gene_list))

#Mural
gene_list <- c("FN1", "COL4A2", "ITGA1")
gene_list <- gene_list[gene_list %in% limma_cogdx_ADvNCI_Mural$gene]

plt <- NULL
for(gene in gene_list){
  plts <- box_fun(limma_cogdx_ADvNCI_Mural, "Cog_Path", c("skyblue", "#F86247", "red3"))
  plt <- plt/plts
}
ggsave2(paste(curr_dir, "/limma/Mural_cellchat.pdf", sep=""), plt, width=1, height=.5*length(gene_list))



##Heatmap

library(ComplexHeatmap)
library(circlize)
library(DESeq2)
library(viridis)
gc()


datprepfun <- function(srt, celltype){
  #Function to select and aggregate data by cell type and projid*Lane
  sub <- subset(srt, Subcelltype==celltype) #Majorcelltype/Subcelltype
  #Aggregate
  gc()
  cts <- AggregateExpression(sub, group.by = "projid", assay = "RNA", slot = "counts", return.seurat = FALSE)
  cts <- cts$RNA
  
  cts <- cts[,order(colnames(cts))]
  
  #Exclude Lanes that have too few counts
  min_to_exclude <- 60000
  low_sum <- colSums(cts)<min_to_exclude #Which cols fall below or above min
  cts <- cts[,low_sum==FALSE] #keep cols above min
  
  assign(paste0("cts_", celltype),  cts, envir = .GlobalEnv)
  
  #Find percentage of cells expressing genes
  pct <- pctincell(sub)
  assign(paste0("pct_", celltype),  pct, envir = .GlobalEnv)
  #write.csv(cts, paste0(curr_dir, "/limma/Pseudobulk_per_projid/cst_", celltype, ".csv"))
  gc()
}

limma_cogdx_ADvNCI_Epithelial_sig <- limma_cogdx_ADvNCI_Epithelial_sig[order(limma_cogdx_ADvNCI_Epithelial_sig$logFC, decreasing=T),]
limma_cogdx_ADvMCI_Epithelial_sig <- limma_cogdx_ADvMCI_Epithelial_sig[order(limma_cogdx_ADvMCI_Epithelial_sig$logFC, decreasing=T),]
limma_cogdx_MCIvNCI_Epithelial_sig <- limma_cogdx_MCIvNCI_Epithelial_sig[order(limma_cogdx_MCIvNCI_Epithelial_sig$logFC, decreasing=T),]

my_genes <- c(limma_cogdx_ADvNCI_Epithelial_sig$gene, limma_cogdx_ADvMCI_Epithelial_sig$gene, limma_cogdx_MCIvNCI_Epithelial_sig$gene) %>% unique()

cts <- data.frame(scale(cts_Epithelial))
cts <- cts[rownames(cts) %in% rownames(cts_atac),]
cts <- cts[order(match(rownames(cts), my_genes)),]
dm <- mdat[order(mdat$Cog_Path, decreasing=T),]
cts <- cts[,order(match(colnames(cts), rownames(dm)))]
cts <- t(scale(t(cts)))
htg <- Heatmap(cts, heatmap_legend_param = list(title = "Expression"), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = TRUE, show_row_names = TRUE, col = col_fun)
draw(htg)
graph2png(htg, file=paste(curr_dir, "/limmma/Heatmap_projid_Epi.png", sep = ""), width=7, height=7)
graph2pdf(htg, file=paste(curr_dir, "/limma/Heatmap_projid_Epi.pdf", sep = ""), width=7, height=7)



##Word cloud

library(wordcloud)
library(RColorBrewer)
library(wordcloud2)
library(stringr)
library(stringi)
library(tm)
library(htmlwidgets)
webshot::install_phantomjs(force = TRUE)

wordcloudfun <- function(FGSEA_pathways, filename){
  #clean up
  terms <- gsub("GOBP_", "", FGSEA_pathways)
  terms <- gsub("GOCC_", "", terms)
  terms <- gsub("GOMF_", "", terms)
  terms <- gsub("POSITIVE_REGULATION_OF_", "", terms)
  terms <- gsub("NEGATIVE_REGULATION_OF_", "", terms)
  terms <- gsub("OR_", "", terms)
  terms <- gsub("VIA_", "", terms)
  terms <- gsub("CELL_", "", terms)
  terms <- c(sapply(str_split(terms,"_",),'[',1),
             sapply(str_split(terms,"_",),'[',2),
             sapply(str_split(terms,"_",),'[',3),
             sapply(str_split(terms,"_",),'[',4), sapply(str_split(terms,"_",),'[',5), sapply(str_split(terms,"_",),'[',6), sapply(str_split(terms,"_",),'[',7))
  terms <- gsub("PROCESS", "", terms)
  terms <- gsub("PATHWAY", "", terms)
  terms <- gsub("RESPONSE", "", terms)
  terms <- gsub("REGULATION", "", terms)
  terms <- gsub("ACTIVITY", "", terms)
  terms <- gsub("LARGE", "", terms)
  terms <- gsub("COMPLEX", "", terms)
  terms <- gsub("FUNCTION", "", terms)
  terms <- gsub("CONTAINING", "", terms)
  terms <- gsub("OF", "", terms)
  terms <- gsub("BY", "", terms)
  terms <- gsub("FROM", "", terms)
  terms <- gsub("BASED", "", terms)
  terms <- gsub("NON", "", terms)
  terms <- gsub("SIGNAL", "SIGNALING", terms)
  terms <- gsub("SIGNALINGING", "SIGNALING", terms)
  terms <- gsub("RIBOSOMAL", "RIBOSOME", terms)
  terms <- gsub("ACTIVATING", "ACTIVATION", terms)
  terms <- gsub("CYTOSKELETAL", "CYTOSKELETON", terms)
  terms <- gsub("DEVELOPMENTAL", "DEVELOPMENT", terms)
  terms <- gsub("MITOCHONDRION", "MITOCHONDRIAL", terms)
  terms <- gsub("ORGANIZING", "ORGANIZATION", terms)
  terms <- gsub("RIBOSOME", "RIBOSOMAL", terms)
  terms <- gsub("SECRETION", "SECRETORY", terms)
  terms <- gsub("STRUCTURE", "STRUCTURAL", terms)
  terms <- gsub("SYNAPTIC", "SYNAPTSE", terms)
  terms <- gsub("TRANSPORTER", "TRANSPORT", terms)
  terms <- gsub("CILIARY", "CILIUM", terms)
  terms <- gsub("CILIA", "CILIUM", terms)
  terms <- gsub("INTRACILIARY", "CILIUM", terms)
  terms <- gsub("METABOLISM", "METABOLIC", terms)
  terms <- as.character(stri_remove_empty(terms, na_empty = FALSE)) %>% stri_remove_na()
  
  #remformat
  dtm <- as.matrix(TermDocumentMatrix(terms))
  dtm <- sort(rowSums(dtm),decreasing=TRUE)
  dtm <- data.frame(word = names(dtm),freq=dtm)
  
  plt <- wordcloud2(dtm, size=1, color="random-dark", backgroundColor = "white", shape="circle")
  saveWidget(plt,"1.html",selfcontained = F)
  webshot::webshot("1.html", paste(curr_dir,"/FGSEA/", filename, ".pdf", sep=""), vwidth = 1992, vheight = 1744, delay =10)
}

#make clouds
wordcloudfun(c(GSEA_cogdx_ADvNCI_Epi_1b$pathway, GSEA_cogdx_ADvNCI_Epi_1a$pathway), "Epi_cloud")
wordcloudfun(c(GSEA_cogdx_ADvNCI_Fib_3$pathway, GSEA_cogdx_ADvNCI_Fib_1$pathway, GSEA_cogdx_ADvNCI_Fib_2$pathway), "Fib_cloud")


##select GO plots

library(ggplot2)
library(fgsea)
library(export)
selectPath_fun <- function(TT_res, selection){
  #TT_res: limma top table result (not cutoff) for fgsea to analyse
  #selection of pathways to plot
  name <- deparse(substitute(TT_res))
  name <- str_sub(name, 7)
  DEG <- TT_res %>% arrange(by_group=t) %>% dplyr::select(gene, t) %>%  na.omit() %>% distinct() %>% deframe()
  fgseaRes <- fgsea(pathways = DB_list, stats = DEG, minSize = 10)
  
  DB_name <- deparse(substitute(DB_list))
  path <- paste(curr_dir, "/FGSEA/", sep="") #specify dir
  if(dir.exists(path)==TRUE){print("Directory exists")}else{dir.create(file.path(path), recursive = TRUE)} #check dir exists
  
  DEG <- DEG[abs(DEG)>4]
  plts <- plotGseaTable(DB_list[selection], DEG, fgseaRes, pathwayLabelStyle=list(size=7), headerLabelStyle=list(size=7), valueStyle=list(size=7), axisLabelStyle=list(size=7))
  graph2pdf(plts, paste(path, "Complex_", name, "_table.pdf", sep=""), font="Arial", width=5, height=1+length(selection)*.7)
  plot(plts)
  remove(plts)
}

#Epi2
selectPath_fun(limma_cogdx_ADvNCI_Epi_1b, c("GOCC_CILIUM", "GOCC_CILIARY_TIP", "GOCC_DYNEIN_COMPLEX", "GOBP_MICROTUBULE_BUNDLE_FORMATION", "GOBP_REGULATION_OF_CILIUM_MOVEMENT", "GOBP_REGULATION_OF_MICROTUBULE_BASED_MOVEMENT", "GOBP_CEREBROSPINAL_FLUID_CIRCULATION", "GOBP_INTRACILIARY_TRANSPORT", "GOBP_EXPORT_FROM_CELL", "GOBP_REGULATION_OF_SECRETION", "GOCC_CATION_CHANNEL_COMPLEX", "GOBP_VESICLE_MEDIATED_TRANSPORT", "GOMF_GATED_CHANNEL_ACTIVITY", "GOBP_PEPTIDE_BIOSYNTHETIC_PROCESS", "GOBP_PHOSPHORYLATION", "GOBP_INFLAMMATORY_RESPONSE", "GOBP_CELL_CELL_SIGNALING", "GOBP_CELL_ADHESION"))

#Fib
selectPath_fun(limma_cogdx_ADvNCI_Fib_3, c("GOCC_CELL_SUBSTRATE_JUNCTION", "GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX", "GOBP_CELL_ADHESION", "GOBP_REGULATION_OF_BINDING", "GOBP_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS", "GOBP_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS", "GOBP_REGULATION_OF_PROTEIN_CATABOLIC_PROCESS", "GOBP_IMMUNE_RESPONSE", "GOBP_MACROAUTOPHAGY", "GOBP_PEPTIDE_BIOSYNTHETIC_PROCESS", "GOMF_CATION_CHANNEL_ACTIVITY", "GOMF_GATED_CHANNEL_ACTIVITY", "GOCC_TRANSPORTER_COMPLEX", "GOBP_REGULATION_OF_INTRACELLULAR_TRANSPORT", "GOCC_ION_CHANNEL_COMPLEX"))

#BAM
selectPath_fun(limma_cogdx_ADvNCI_BAM_1, c("GOBP_REGULATION_OF_INFLAMMATORY_RESPONSE", "GOBP_REGULATION_OF_STEROL_TRANSPORT", "GOBP_CELL_MOTILITY", "GOBP_CELLULAR_MACROMOLECULE_CATABOLIC_PROCESS"))


##GO heatmap of overlap

library(ComplexHeatmap)
library(circlize)
library(DESeq2)
library(viridis)
library(export)
col_fun <-  colorRamp2(c(-2, 0, 2), c("blue","white", "red"))

expression_fun <- function(cts){
  cts <- cts[,-c(2:7, (length(cts)-2):length(cts))]
  cts <- cts[!(duplicated(cts$gene)),]
  rownames(cts) <- cts$gene
  cts$gene <- NULL
  colnames(cts) <- substr(colnames(cts), 2, 20)
  
  #Select relevant rows from Design Matrix  
  dm <- matrix(ncol=ncol(mdat))
  dm <- dm[-1,]
  colnames(dm) <- colnames(mdat)
  
  for(i in colnames(cts)){
    base_name <- sapply(str_split(i,"_",),'[',1) #get projid names
    if(base_name %in% rownames(mdat)){
      name_rows <- rownames(dm)
      row_add <- mdat[base_name,]
      dm <- rbind(dm, row_add)
      rownames(dm) <- c(name_rows, paste(i))
    }else{print(paste("Warning colname of data:", i, "not in rowname of Design Matrix"))}
  }
  
  #Enforce Maindm and cts name match
  if(all(rownames(dm)==colnames(cts))==FALSE){
    dm <- dm[rownames(dm) %in% colnames(cts),]      
    cts <- cts[,colnames(cts) %in% rownames(dm)]
  }
  
  #order dm and cts based on first variable of interest
  dm <- dm[order(dm$Group),]
  cts <- cts[,order(match(colnames(cts), rownames(dm)))]
  unique(colnames(cts)==rownames(dm))
  
  #gene_means <- tapply(as.numeric(cts["ABCA1", ]), dm$Group, mean) #test
  
  Rmeans <- as.data.frame(sapply(unique(dm$Group), function(grp) {
    rowMeans(cts[, dm$Group == grp, drop=FALSE])
  }))
  
  select <- pathway[selection]
  selectmeans <- data.frame(sapply(select, function(x){
    colMeans(Rmeans[rownames(Rmeans) %in% unique(x[[1]]),])
  }))
  selectmeans <- data.frame(t(scale(selectmeans)))
  selectmeans <- as.matrix(selectmeans[,c("NCI_33", "NCI_34", "MCI_33", "MCI_34", "ADd_33", "ADd_34")])
  return(selectmeans)
}
mdat <- data.frame(read_csv("/Ex_ChP/Ex_ChP7/metadata/ex vivo ChP metadata (analytical) newAPOE RNAonly 3334.csv"))
rownames(mdat) <- mdat$projid
mdat$Group <- paste(mdat$cogdx_simple, mdat$apoe_genotype, sep = "_")


###Overall APOE

selection <- c(      #IMMUNE
  "GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS", "GOBP_ADAPTIVE_IMMUNE_RESPONSE", "GOBP_REGULATION_OF_IMMUNE_RESPONSE", #general immune
  "GOBP_REGULATION_OF_RESPONSE_TO_STRESS", "GOBP_RESPONSE_TO_WOUNDING", "GOMF_CYTOKINE_RECEPTOR_ACTIVITY", #response immune
  "GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE", "GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE","GOBP_T_CELL_ACTIVATION",  #Tcell, Bcell, 
  "GOBP_INTERLEUKIN_8_PRODUCTION", "GOBP_INTERLEUKIN_6_PRODUCTION", "GOBP_PHAGOCYTOSIS",
  
  #ADHESION
  "GOMF_CELL_ADHESION_MOLECULE_BINDING", "GOBP_CELL_CELL_ADHESION", "GOBP_CELL_MATRIX_ADHESION", "GOBP_CELL_SUBSTRATE_ADHESION", #adhesion
  "GOCC_ANCHORING_JUNCTION", "GOCC_APICAL_JUNCTION_COMPLEX", "GOCC_CELL_SUBSTRATE_JUNCTION",#tight junctions
  "GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX",  "GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY", #ECM
  "GOCC_ACTIN_CYTOSKELETON", "GOBP_ACTIN_FILAMENT_BASED_PROCESS", "GOCC_ACTIN_FILAMENT_BUNDLE", "GOBP_ACTIN_FILAMENT_ORGANIZATION", "GOMF_ACTIN_BINDING", #actin
  "GOCC_POSTSYNAPTIC_DENSITY_MEMBRANE", "GOMF_CADHERIN_BINDING", #other adhesion
  
  #EXPORT/SECRETION/TRANSPORT
  "GOBP_VESICLE_MEDIATED_TRANSPORT", "GOBP_SECRETION", "GOBP_EXTRACELLULAR_TRANSPORT",  #export/secretion
  "GOMF_MOLECULAR_TRANSDUCER_ACTIVITY", #general transport
  "GOCC_CATION_CHANNEL_COMPLEX", "GOMF_VOLTAGE_GATED_CHANNEL_ACTIVITY",#transport related to signaling
  
  #METABOLISM
  "GOBP_STEROL_TRANSPORT", "GOBP_STEROL_HOMEOSTASIS", "GOBP_LIPID_LOCALIZATION", "GOBP_REGULATION_OF_LIPID_METABOLIC_PROCESS", "GOBP_CHOLESTEROL_STORAGE",#str
  
  #CILIA
  "GOBP_CEREBROSPINAL_FLUID_CIRCULATION", "GOCC_CILIUM", "GOCC_CILIARY_TIP", "GOBP_CILIUM_MOVEMENT", "GOCC_CYTOPLASMIC_MICROTUBULE", "GOBP_MICROTUBULE_BASED_MOVEMENT", "GOMF_TUBULIN_BINDING", "GOCC_DYNEIN_COMPLEX", "GOBP_INTRACILIARY_TRANSPORT", "GOCC_INTRACILIARY_TRANSPORT_PARTICLE_B")


##per contrast
GSEA_join <- full_join(GO_labelled, GSEA_apoe_gpath_gpathxAPOE_34v33_Epithelial, by="pathway", suffix=c("", "_Epithelial")) %>% full_join(GSEA_apoe_gpath_gpathxAPOE_34v33_Fibroblast, by="pathway", suffix=c("", "_Fibroblast")) %>% full_join(GSEA_apoe_gpath_gpathxAPOE_34v33_Endothelial, by="pathway", suffix=c("", "_Endothelial")) %>% full_join(GSEA_apoe_gpath_gpathxAPOE_34v33_Mural, by="pathway", suffix=c("", "_Mural")) %>% full_join(GSEA_apoe_gpath_gpathxAPOE_34v33_Immune, by="pathway", suffix=c("", "_Immune"))

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

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = col_fun, cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.1) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Ex vivo APOE heatmap apoe.pdf", sep = ""), font="Arial", width=w, height=h)



###Epi

#Organoid first
selection <- c("GOBP_LIPID_BIOSYNTHETIC_PROCESS", "GOBP_STEROID_METABOLIC_PROCESS", "GOBP_REGULATION_OF_CHOLESTEROL_METABOLIC_PROCESS", "GOCC_LIPID_DROPLET", "GOBP_LIPID_METABOLIC_PROCESS", "GOBP_CELLULAR_LIPID_METABOLIC_PROCESS", "KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS",
               "GOBP_INFLAMMATORY_RESPONSE", "GOBP_RESPONSE_TO_WOUNDING", "GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
               "GOBP_SECRETION", "GOBP_REGULATION_OF_PROTEIN_SECRETION", "GOBP_STEROL_TRANSPORT",
               "GOBP_CEREBROSPINAL_FLUID_CIRCULATION", "GOCC_CILIUM", "GOBP_CILIUM_MOVEMENT", "GOBP_MICROTUBULE_BASED_MOVEMENT", "GOMF_TUBULIN_BINDING", "GOCC_DYNEIN_COMPLEX", "GOBP_INTRACILIARY_TRANSPORT",
               "GOMF_CELL_ADHESION_MOLECULE_BINDING", "KEGG_TIGHT_JUNCTION", "GOCC_ANCHORING_JUNCTION", "GOCC_CELL_SUBSTRATE_JUNCTION", "GOMF_CADHERIN_BINDING")

#Ex vivo first
selection <- c("GOBP_PHOSPHOLIPID_TRANSPORT", "GOBP_LIPID_LOCALIZATION", "GOBP_REGULATION_OF_LIPID_METABOLIC_PROCESS", "GOBP_STEROID_METABOLIC_PROCESS", "GOBP_REGULATION_OF_LIPID_CATABOLIC_PROCESS", "GOBP_LIPID_BIOSYNTHETIC_PROCESS", "GOBP_REGULATION_OF_CHOLESTEROL_METABOLIC_PROCESS", 
               "GOBP_INFLAMMATORY_RESPONSE", "GOBP_RESPONSE_TO_WOUNDING", "GOBP_T_CELL_DIFFERENTIATION",
               "GOMF_CELL_ADHESION_MOLECULE_BINDING", "GOCC_ANCHORING_JUNCTION", "GOCC_CELL_SUBSTRATE_JUNCTION", "GOMF_CADHERIN_BINDING",
               "GOBP_SECRETION", "GOBP_REGULATION_OF_PROTEIN_SECRETION", "GOBP_STEROL_TRANSPORT",
               "GOBP_CEREBROSPINAL_FLUID_CIRCULATION", "GOCC_CILIUM", "GOBP_CILIUM_MOVEMENT", "GOBP_MICROTUBULE_BASED_MOVEMENT", "GOMF_TUBULIN_BINDING", "GOCC_DYNEIN_COMPLEX", "GOBP_INTRACILIARY_TRANSPORT")


##per contrast
GSEA_join <- full_join(GO_labelled, GSEA_cogdxAPOE_34v33NCI_Epithelial, by="pathway", suffix=c("", "_apoe33v34 in NCI")) %>% full_join(GSEA_cogdxAPOE_ADvNCI_Epithelial, by="pathway", suffix=c("", "_ADvNCI in apoe33")) %>% full_join(GSEA_cogdxAPOE_ADvNCIin34_Epithelial, by="pathway", suffix=c("", "_ADvNCI in apoe34")) %>% full_join(GSEA_cogdxAPOE_cogdxapoe_Epithelial, by="pathway", suffix=c("", "_ADvNCI in apoe34vapoe33"))# %>% full_join(GSEA_cogdx_ADvNCI_Epithelial, by="pathway", suffix=c("", "_cogdx_ADvNCI"))

GSEA_join <- GSEA_join %>% column_to_rownames("pathway")
GSEA_join <- GSEA_join[selection,]
GSEA_join <- GSEA_join[match(selection, rownames(GSEA_join)),]
mat <- GSEA_join[,c(grep(pattern = "NES", x = colnames(GSEA_join), value = T))]
colnames(mat) <- c("apoe33v34 in NCI", colnames(mat)[2:ncol(mat)]) #!!!!apoe33v34 in NCI
colnames(mat) <- gsub("NES_", "", colnames(mat))
mat[is.na(mat)] <- 0
mat[] <- lapply(mat, as.numeric)

mat2 <- GSEA_join[,c(grep(pattern = "padj", x = colnames(GSEA_join), value = T))]
colnames(mat2) <- c("apoe33v34 in NCI", colnames(mat2)[2:ncol(mat2)]) #!!!!apoe33v34 in NCI
colnames(mat2) <- gsub("padj_", "", colnames(mat2))
mat2[is.na(mat2)] <- 1
mat2[] <- lapply(mat2, as.numeric)

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = col_fun, cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.1) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/APOE heatmap Epithelial.pdf", sep = ""), font="Arial", width=w, height=h)



####Pathway expression

pathway1 <- split(x=GSEA_cogdx_33_ADvNCI_Epithelial$leadingEdge, f=GSEA_cogdx_33_ADvNCI_Epithelial$pathway)
pathway2 <- split(x=GSEA_cogdx_33_ADvMCI_Epithelial$leadingEdge, f=GSEA_cogdx_33_ADvMCI_Epithelial$pathway)
pathway3 <- split(x=GSEA_cogdx_33_MCIvNCI_Epithelial$leadingEdge, f=GSEA_cogdx_33_MCIvNCI_Epithelial$pathway)
pathway4 <- split(x=GSEA_apoe_cogdxNCI_34v33_Epithelial$leadingEdge, f=GSEA_apoe_cogdxNCI_34v33_Epithelial$pathway)
pathway5 <- split(x=GSEA_apoe_cogdxMCI_34v33_Epithelial$leadingEdge, f=GSEA_apoe_cogdxMCI_34v33_Epithelial$pathway)
pathway6 <- split(x=GSEA_apoe_cogdxAD_34v33_Epithelial$leadingEdge, f=GSEA_apoe_cogdxAD_34v33_Epithelial$pathway)

pathway <- modifyList(pathway1, c(pathway2, pathway3, pathway4, pathway5, pathway6), keep.null=T)
pathway <- lapply(pathway, unique)

selectmeans <- expression_fun(limma_cogdx_ADvNCI_Epithelial)

plts <- Heatmap(selectmeans, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(selectmeans)*unit(2.73, "mm"), height = nrow(selectmeans)*unit(2.73, "mm"), col=col_fun)
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/APOE expression heatmap Epithelial.pdf", sep = ""), font="Arial", width=w, height=h)


####per celltype

GSEA_EpiType_join <- full_join(GO_labelled, GSEA_apoe_cogdxNCI_34v33_Epi_1a, by="pathway", suffix=c("", "_Epi_1a")) %>% full_join(GSEA_apoe_cogdxNCI_34v33_Epi_1b, by="pathway", suffix=c("", "_Epi_1b")) %>% full_join(GSEA_apoe_cogdxNCI_34v33_Epi_2a, by="pathway", suffix=c("", "_Epi_2a")) %>% full_join(GSEA_apoe_cogdxNCI_34v33_Epi_2b, by="pathway", suffix=c("", "_Epi_2b"))

GSEA_EpiType_join <- GSEA_EpiType_join %>% column_to_rownames("pathway")
GSEA_EpiType_join <- GSEA_EpiType_join[selection,]
GSEA_EpiType_join <- GSEA_EpiType_join[match(selection, rownames(GSEA_EpiType_join)),]
mat <- GSEA_EpiType_join[,c("NES", "NES_Epi_1b", "NES_Epi_2a", "NES_Epi_2b")]
colnames(mat) <- c("Epi_1a", "Epi_1b", "Epi_2a", "Epi_2b")
mat[is.na(mat)] <- 0
mat[] <- lapply(mat, as.numeric)

mat2 <- GSEA_EpiType_join[,c("padj", "padj_Epi_1b", "padj_Epi_2a", "padj_Epi_2b")]
colnames(mat2) <- c("Epi_1a", "Epi_1b", "Epi_2a", "Epi_2b")
mat2[is.na(mat2)] <- 1
mat2[] <- lapply(mat2, as.numeric)

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = col_fun, cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.1) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
plot(plts)
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_EpiTypes.pdf", sep = ""), font="Arial", width=w, height=h)



###Fib

selection <- c("GOBP_CELL_MATRIX_ADHESION", "GOCC_CELL_SURFACE",
               "GOMF_IMMUNE_RECEPTOR_ACTIVITY", "GOBP_RESPONSE_TO_WOUNDING", "GOBP_INTERLEUKIN_6_PRODUCTION", "GOBP_INTERLEUKIN_8_PRODUCTION")

#Ex vivo first
selection <- c("GOBP_CELL_MATRIX_ADHESION", "GOBP_CELL_SUBSTRATE_ADHESION", "GOCC_ACTIN_CYTOSKELETON",
               "GOMF_IMMUNE_RECEPTOR_ACTIVITY", "GOBP_RESPONSE_TO_WOUNDING", "GOBP_INTERLEUKIN_6_PRODUCTION", "GOBP_LEUKOCYTE_CELL_CELL_ADHESION")

##per contrast
GSEA_join <- full_join(GO_labelled, GSEA_cogdxAPOE_34v33NCI_Fibroblast, by="pathway", suffix=c("", "_apoe33v34 in NCI")) %>% full_join(GSEA_cogdxAPOE_ADvNCI_Fibroblast, by="pathway", suffix=c("", "_ADvNCI in apoe33")) %>% full_join(GSEA_cogdxAPOE_ADvNCIin34_Fibroblast, by="pathway", suffix=c("", "_ADvNCI in apoe34")) %>% full_join(GSEA_cogdxAPOE_cogdxapoe_Fibroblast, by="pathway", suffix=c("", "_ADvNCI in apoe34vapoe33"))# %>% full_join(GSEA_cogdx_ADvNCI_Fibroblast, by="pathway", suffix=c("", "_cogdx_ADvNCI"))

GSEA_join <- GSEA_join %>% column_to_rownames("pathway")
GSEA_join <- GSEA_join[selection,]
GSEA_join <- GSEA_join[match(selection, rownames(GSEA_join)),]
mat <- GSEA_join[,c(grep(pattern = "NES", x = colnames(GSEA_join), value = T))]
colnames(mat) <- c("apoe33v34 in NCI", colnames(mat)[2:ncol(mat)]) #!!!!apoe33v34 in NCI
colnames(mat) <- gsub("NES_", "", colnames(mat))
mat[is.na(mat)] <- 0
mat[] <- lapply(mat, as.numeric)

mat2 <- GSEA_join[,c(grep(pattern = "padj", x = colnames(GSEA_join), value = T))]
colnames(mat2) <- c("apoe33v34 in NCI", colnames(mat2)[2:ncol(mat2)]) #!!!!apoe33v34 in NCI
colnames(mat2) <- gsub("padj_", "", colnames(mat2))
mat2[is.na(mat2)] <- 1
mat2[] <- lapply(mat2, as.numeric)
plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = col_fun, cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.1) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/APOE heatmap Fibroblast.pdf", sep = ""), font="Arial", width=w, height=h)


####Pathway expression

pathway1 <- split(x=GSEA_cogdx_33_ADvNCI_Fibroblast$leadingEdge, f=GSEA_cogdx_33_ADvNCI_Fibroblast$pathway)
pathway2 <- split(x=GSEA_cogdx_33_ADvMCI_Fibroblast$leadingEdge, f=GSEA_cogdx_33_ADvMCI_Fibroblast$pathway)
pathway3 <- split(x=GSEA_cogdx_33_MCIvNCI_Fibroblast$leadingEdge, f=GSEA_cogdx_33_MCIvNCI_Fibroblast$pathway)
pathway4 <- split(x=GSEA_apoe_cogdxNCI_34v33_Fibroblast$leadingEdge, f=GSEA_apoe_cogdxNCI_34v33_Fibroblast$pathway)
pathway5 <- split(x=GSEA_apoe_cogdxMCI_34v33_Fibroblast$leadingEdge, f=GSEA_apoe_cogdxMCI_34v33_Fibroblast$pathway)
pathway6 <- split(x=GSEA_apoe_cogdxAD_34v33_Fibroblast$leadingEdge, f=GSEA_apoe_cogdxAD_34v33_Fibroblast$pathway)

pathway <- modifyList(pathway1, c(pathway2, pathway3, pathway4, pathway5, pathway6), keep.null=T)
pathway <- lapply(pathway, unique)

selectmeans <- expression_fun(limma_cogdx_ADvNCI_Fibroblast)

plts <- Heatmap(selectmeans, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(selectmeans)*unit(2.73, "mm"), height = nrow(selectmeans)*unit(2.73, "mm"), col=col_fun)
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/APOE expression heatmap Fibroblast.pdf", sep = ""), font="Arial", width=w, height=h)


####per celltype

GSEA_FibType_join <- full_join(GO_labelled, GSEA_apoe_cogdxNCI_34v33_Fib_1, by="pathway", suffix=c("", "_Fib_1")) %>% full_join(GSEA_apoe_cogdxNCI_34v33_Fib_2, by="pathway", suffix=c("", "_Fib_2")) %>% full_join(GSEA_apoe_cogdxNCI_34v33_Fib_3, by="pathway", suffix=c("", "_Fib_3"))

GSEA_FibType_join <- GSEA_FibType_join %>% column_to_rownames("pathway")
GSEA_FibType_join <- GSEA_FibType_join[selection,]
GSEA_FibType_join <- GSEA_FibType_join[match(selection, rownames(GSEA_FibType_join)),]
mat <- GSEA_FibType_join[,c("NES", "NES_Fib_2", "NES_Fib_3")]
colnames(mat) <- c("Fib_1", "Fib_2", "Fib_3")
mat[is.na(mat)] <- 0
mat[] <- lapply(mat, as.numeric)

mat2 <- GSEA_FibType_join[,c("padj", "padj_Fib_2", "padj_Fib_3")]
colnames(mat2) <- c("Fib_1", "Fib_2", "Fib_3")
mat2[is.na(mat2)] <- 1
mat2[] <- lapply(mat2, as.numeric)

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = col_fun, cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.1) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_FibTypes apoe34v33 within NCI.pdf", sep = ""), font="Arial", width=w, height=h)


###Immune

selection <- c("GOBP_MYELOID_LEUKOCYTE_ACTIVATION", "GOBP_RESPONSE_TO_CYTOKINE", "GOBP_RESPONSE_TO_WOUNDING",
               "GOBP_LIPID_LOCALIZATION", "GOBP_CHOLESTEROL_STORAGE", "GOBP_LIPID_IMPORT_INTO_CELL")

##per contrast
GSEA_join <- full_join(GO_labelled, GSEA_cogdxAPOE_34v33NCI_BAM, by="pathway", suffix=c("", "_apoe33v34 in NCI")) %>% full_join(GSEA_cogdxAPOE_ADvNCI_BAM, by="pathway", suffix=c("", "_ADvNCI in apoe33")) %>% full_join(GSEA_cogdxAPOE_ADvNCIin34_BAM, by="pathway", suffix=c("", "_ADvNCI in apoe34")) %>% full_join(GSEA_cogdxAPOE_cogdxapoe_BAM, by="pathway", suffix=c("", "_ADvNCI in apoe34vapoe33"))# %>% full_join(GSEA_cogdx_ADvNCI_BAM, by="pathway", suffix=c("", "_cogdx_ADvNCI"))

GSEA_join <- GSEA_join %>% column_to_rownames("pathway")
GSEA_join <- GSEA_join[selection,]
GSEA_join <- GSEA_join[match(selection, rownames(GSEA_join)),]
mat <- GSEA_join[,c(grep(pattern = "NES", x = colnames(GSEA_join), value = T))]
colnames(mat) <- c("apoe33v34 in NCI", colnames(mat)[2:ncol(mat)]) #!!!!apoe33v34 in NCI
colnames(mat) <- gsub("NES_", "", colnames(mat))
mat[is.na(mat)] <- 0
mat[] <- lapply(mat, as.numeric)

mat2 <- GSEA_join[,c(grep(pattern = "padj", x = colnames(GSEA_join), value = T))]
colnames(mat2) <- c("apoe33v34 in NCI", colnames(mat2)[2:ncol(mat2)]) #!!!!apoe33v34 in NCI
colnames(mat2) <- gsub("padj_", "", colnames(mat2))
mat2[is.na(mat2)] <- 1
mat2[] <- lapply(mat2, as.numeric)
plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = col_fun, cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.1) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/APOE heatmap BAM.pdf", sep = ""), font="Arial", width=w, height=h)


####Pathway expression

pathway1 <- split(x=GSEA_cogdx_33_ADvNCI_Immune$leadingEdge, f=GSEA_cogdx_33_ADvNCI_Immune$pathway)
pathway2 <- split(x=GSEA_cogdx_33_ADvMCI_Immune$leadingEdge, f=GSEA_cogdx_33_ADvMCI_Immune$pathway)
pathway3 <- split(x=GSEA_cogdx_33_MCIvNCI_Immune$leadingEdge, f=GSEA_cogdx_33_MCIvNCI_Immune$pathway)
pathway4 <- split(x=GSEA_apoe_cogdxNCI_34v33_Immune$leadingEdge, f=GSEA_apoe_cogdxNCI_34v33_Immune$pathway)
pathway5 <- split(x=GSEA_apoe_cogdxMCI_34v33_Immune$leadingEdge, f=GSEA_apoe_cogdxMCI_34v33_Immune$pathway)
pathway6 <- split(x=GSEA_apoe_cogdxAD_34v33_Immune$leadingEdge, f=GSEA_apoe_cogdxAD_34v33_Immune$pathway)

pathway <- modifyList(pathway1, c(pathway2, pathway3, pathway4, pathway5, pathway6), keep.null=T)
pathway <- lapply(pathway, unique)

selectmeans <- expression_fun(limma_cogdx_ADvNCI_Immune)

plts <- Heatmap(selectmeans, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(selectmeans)*unit(2.73, "mm"), height = nrow(selectmeans)*unit(2.73, "mm"), col=col_fun)
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/APOE expression heatmap Immune.pdf", sep = ""), font="Arial", width=w, height=h)


###Endo

selection <- c("GOBP_PEPTIDE_BIOSYNTHETIC_PROCESS", "GOBP_PEPTIDE_METABOLIC_PROCESS")

##per contrast
GSEA_join <- full_join(GO_labelled, GSEA_cogdxAPOE_ADvNCI_Endothelial, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_cogdxAPOE_34v33NCI_Endothelial, by="pathway", suffix=c("", "_apoe33v34 in NCI")) %>% full_join(GSEA_cogdxAPOE_34v33ADd_Endothelial, by="pathway", suffix=c("", "_apoe33v34 in ADd"))

GSEA_join <- GSEA_join %>% column_to_rownames("pathway")
GSEA_join <- GSEA_join[selection,]
GSEA_join <- GSEA_join[match(selection, rownames(GSEA_join)),]
mat <- GSEA_join[,c(grep(pattern = "NES", x = colnames(GSEA_join), value = T))]
colnames(mat) <- c("apoe33v34 in NCI", colnames(mat)[2:ncol(mat)])
colnames(mat) <- gsub("NES_", "", colnames(mat))
mat[is.na(mat)] <- 0
mat[] <- lapply(mat, as.numeric)

mat2 <- GSEA_join[,c(grep(pattern = "padj", x = colnames(GSEA_join), value = T))]
colnames(mat2) <- c("apoe33v34 in NCI", colnames(mat2)[2:ncol(mat2)])
colnames(mat2) <- gsub("padj_", "", colnames(mat2))
mat2[is.na(mat2)] <- 1
mat2[] <- lapply(mat2, as.numeric)

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = col_fun, cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.1) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/APOE heatmap Endothelial.pdf", sep = ""), font="Arial", width=w, height=h)


####Pathway expression

pathway1 <- split(x=GSEA_cogdx_33_ADvNCI_Endothelial$leadingEdge, f=GSEA_cogdx_33_ADvNCI_Endothelial$pathway)
pathway2 <- split(x=GSEA_cogdx_33_ADvMCI_Endothelial$leadingEdge, f=GSEA_cogdx_33_ADvMCI_Endothelial$pathway)
pathway3 <- split(x=GSEA_cogdx_33_MCIvNCI_Endothelial$leadingEdge, f=GSEA_cogdx_33_MCIvNCI_Endothelial$pathway)
pathway4 <- split(x=GSEA_apoe_cogdxNCI_34v33_Endothelial$leadingEdge, f=GSEA_apoe_cogdxNCI_34v33_Endothelial$pathway)
pathway6 <- split(x=GSEA_apoe_cogdxAD_34v33_Endothelial$leadingEdge, f=GSEA_apoe_cogdxAD_34v33_Endothelial$pathway)

pathway <- modifyList(pathway1, c(pathway2, pathway3, pathway4, pathway6), keep.null=T)
pathway <- lapply(pathway, unique)


#no MCI_34
expression_fun <- function(cts){
  cts <- cts[,-c(2:7, (length(cts)-2):length(cts))]
  cts <- cts[!(duplicated(cts$gene)),]
  rownames(cts) <- cts$gene
  cts$gene <- NULL
  colnames(cts) <- substr(colnames(cts), 2, 20)
  
  #Select relevant rows from Design Matrix  
  dm <- matrix(ncol=ncol(mdat))
  dm <- dm[-1,]
  colnames(dm) <- colnames(mdat)
  
  for(i in colnames(cts)){
    base_name <- sapply(str_split(i,"_",),'[',1) #get projid names
    if(base_name %in% rownames(mdat)){
      name_rows <- rownames(dm)
      row_add <- mdat[base_name,]
      dm <- rbind(dm, row_add)
      rownames(dm) <- c(name_rows, paste(i))
    }else{print(paste("Warning colname of data:", i, "not in rowname of Design Matrix"))}
  }
  
  #Enforce Maindm and cts name match
  if(all(rownames(dm)==colnames(cts))==FALSE){
    dm <- dm[rownames(dm) %in% colnames(cts),]      
    cts <- cts[,colnames(cts) %in% rownames(dm)]
  }
  
  #order dm and cts based on first variable of interest
  dm <- dm[order(dm$Group),]
  cts <- cts[,order(match(colnames(cts), rownames(dm)))]
  unique(colnames(cts)==rownames(dm))
  
  #gene_means <- tapply(as.numeric(cts["ABCA1", ]), dm$Group, mean) #test
  
  Rmeans <- as.data.frame(sapply(unique(dm$Group), function(grp) {
    rowMeans(cts[, dm$Group == grp, drop=FALSE])
  }))
  
  select <- pathway[selection]
  selectmeans <- data.frame(sapply(select, function(x){
    colMeans(Rmeans[rownames(Rmeans) %in% unique(x[[1]]),])
  }))
  selectmeans <- data.frame(t(scale(selectmeans)))
  selectmeans <- as.matrix(selectmeans[,c("NCI_33", "NCI_34", "MCI_33", "ADd_33", "ADd_34")])
  return(selectmeans)
}

selectmeans <- expression_fun(limma_cogdx_ADvNCI_Endothelial)

plts <- Heatmap(selectmeans, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(selectmeans)*unit(2.73, "mm"), height = nrow(selectmeans)*unit(2.73, "mm"), col=col_fun)
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/APOE expression heatmap Endothelial.pdf", sep = ""), font="Arial", width=w, height=h)



####per celltype

GSEA_EndoType_join <- full_join(GO_labelled, GSEA_cogdx_ADvNCI_Endo_1, by="pathway", suffix=c("", "_Endo_1")) %>% full_join(GSEA_cogdx_ADvNCI_Endo_2, by="pathway", suffix=c("", "_Endo_2")) %>% full_join(GSEA_cogdx_ADvNCI_Endo_3, by="pathway", suffix=c("", "_Endo_3"))

GSEA_EndoType_join <- GSEA_EndoType_join %>% column_to_rownames("pathway")
GSEA_EndoType_join <- GSEA_EndoType_join[selection,]
GSEA_EndoType_join <- GSEA_EndoType_join[match(selection, rownames(GSEA_EndoType_join)),]
mat <- GSEA_EndoType_join[,c("NES", "NES_Endo_2", "NES_Endo_3")]
colnames(mat) <- c("Endo_1", "Endo_2", "Endo_3")
mat[is.na(mat)] <- 0
mat[] <- lapply(mat, as.numeric)

mat2 <- GSEA_EndoType_join[,c("padj", "padj_Endo_2", "padj_Endo_3")]
colnames(mat2) <- c("Endo_1", "Endo_2", "Endo_3")
mat2[is.na(mat2)] <- 1
mat2[] <- lapply(mat2, as.numeric)

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = col_fun, cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.1) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Overlap_heatmap_EndoTypes.pdf", sep = ""), font="Arial", width=w, height=h)


###Mural

selection <- c("GOBP_EXTRACELLULAR_TRANSPORT", "GOMF_SECONDARY_ACTIVE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",
               "GOBP_CHAPERONE_MEDIATED_PROTEIN_FOLDING", "GOMF_UNFOLDED_PROTEIN_BINDING") #immune

##per contrast
GSEA_join <- full_join(GO_labelled, GSEA_cogdxAPOE_ADvNCI_Mural, by="pathway", suffix=c("", "_ADvNCI")) %>% full_join(GSEA_cogdxAPOE_34v33NCI_Mural, by="pathway", suffix=c("", "_apoe33v34 in NCI")) %>% full_join(GSEA_cogdxAPOE_34v33ADd_Mural, by="pathway", suffix=c("", "_apoe33v34 in ADd"))

GSEA_join <- GSEA_join %>% column_to_rownames("pathway")
GSEA_join <- GSEA_join[selection,]
GSEA_join <- GSEA_join[match(selection, rownames(GSEA_join)),]
mat <- GSEA_join[,c(grep(pattern = "NES", x = colnames(GSEA_join), value = T))]
colnames(mat) <- c("apoe33v34 in NCI", colnames(mat)[2:ncol(mat)])
colnames(mat) <- gsub("NES_", "", colnames(mat))
mat[is.na(mat)] <- 0
mat[] <- lapply(mat, as.numeric)

mat2 <- GSEA_join[,c(grep(pattern = "padj", x = colnames(GSEA_join), value = T))]
colnames(mat2) <- c("apoe33v34 in NCI", colnames(mat2)[2:ncol(mat2)])
colnames(mat2) <- gsub("padj_", "", colnames(mat2))
mat2[is.na(mat2)] <- 1
mat2[] <- lapply(mat2, as.numeric)

plts <- Heatmap(mat, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(mat)*unit(2.7, "mm"), height = nrow(mat)*unit(2.7, "mm"), col = col_fun, cell_fun = function(j, i, x, y, w, h, fill){
  if(mat2[i, j] < 0.1) {
    grid.text("*", x, y, gp=gpar(fontsize=14), vjust=unit(.81, "mm"))
  }
})
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/APOE heatmap Mural.pdf", sep = ""), font="Arial", width=w, height=h)


####Pathway expression

pathway1 <- split(x=GSEA_cogdx_33_ADvNCI_Mural$leadingEdge, f=GSEA_cogdx_33_ADvNCI_Mural$pathway)
pathway2 <- split(x=GSEA_cogdx_33_ADvMCI_Mural$leadingEdge, f=GSEA_cogdx_33_ADvMCI_Mural$pathway)
pathway3 <- split(x=GSEA_cogdx_33_MCIvNCI_Mural$leadingEdge, f=GSEA_cogdx_33_MCIvNCI_Mural$pathway)
pathway4 <- split(x=GSEA_apoe_cogdxNCI_34v33_Mural$leadingEdge, f=GSEA_apoe_cogdxNCI_34v33_Mural$pathway)
pathway6 <- split(x=GSEA_apoe_cogdxAD_34v33_Mural$leadingEdge, f=GSEA_apoe_cogdxAD_34v33_Mural$pathway)

pathway <- modifyList(pathway1, c(pathway2, pathway3, pathway4, pathway6), keep.null=T)
pathway <- lapply(pathway, unique)

selectmeans <- expression_fun(limma_cogdx_ADvNCI_Mural)

plts <- Heatmap(selectmeans, cluster_rows = F, cluster_columns = F, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7), heatmap_legend_param = list(title="NES", title_gp=gpar(fontsize = 7), labels_gp=gpar(fontsize=7)), width = ncol(selectmeans)*unit(2.73, "mm"), height = nrow(selectmeans)*unit(2.73, "mm"), col=col_fun)
w = ComplexHeatmap:::width(draw(plts))
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(draw(plts))
h = convertY(h, "inch", valueOnly = TRUE)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/APOE expression heatmap Mural.pdf", sep = ""), font="Arial", width=w, height=h)


###Pathway summary figure

temp <- data.frame("Epithelial"=c(nrow(GSEA_apoe_gpath_gpathxAPOE_34v33_Epithelial_sig),
                                  nrow(GSEA_apoe_gpath_gpathxAPOE_gpath_Epithelial_sig),
                                  nrow(GSEA_apoe_gpath_gpathxAPOE_apoexgpath_Epithelial_sig)),
                   "Fibroblast"=c(nrow(GSEA_apoe_gpath_gpathxAPOE_34v33_Fibroblast_sig),
                                  nrow(GSEA_apoe_gpath_gpathxAPOE_gpath_Fibroblast_sig),
                                  nrow(GSEA_apoe_gpath_gpathxAPOE_apoexgpath_Fibroblast_sig)),
                   "Endothelial"=c(nrow(GSEA_apoe_gpath_gpathxAPOE_34v33_Endothelial_sig),
                                   nrow(GSEA_apoe_gpath_gpathxAPOE_gpath_Endothelial_sig),
                                   nrow(GSEA_apoe_gpath_gpathxAPOE_apoexgpath_Endothelial_sig)),
                   "Mural"=c(nrow(GSEA_apoe_gpath_gpathxAPOE_34v33_Mural_sig),
                             nrow(GSEA_apoe_gpath_gpathxAPOE_gpath_Mural_sig),
                             nrow(GSEA_apoe_gpath_gpathxAPOE_apoexgpath_Mural_sig)),
                   "Immune"=c(nrow(GSEA_apoe_gpath_gpathxAPOE_34v33_Immune_sig),
                              nrow(GSEA_apoe_gpath_gpathxAPOE_gpath_Immune_sig),
                              nrow(GSEA_apoe_gpath_gpathxAPOE_apoexgpath_Immune_sig)))
temp$contrast <- c("APOE34v33", "gpath", "apoe*gpath")
temp <- pivot_longer(temp, cols=1:5)

plts <- ggplot(temp, aes(value, name, color=contrast)) + geom_point() + 
  scale_color_manual(values = c("#440154FF", "#22A884FF", "#FDE725FF")) + theme_minimal() +
  labs(x = "Number of significant pathways", y = "Major cell type", title = "Cell type responses to apoe and gpath")
plot(plts)
graph2pdf(plts, file=paste(curr_dir, "/FGSEA/Number of pathways to celltype apoe by gpath.pdf", sep = ""), font="Arial", width=w, height=h)


#GSVA
###Make snRNAseq modules

#Go term matrices
Temp <- cbind(GSEA_cogdxAPOE_34v33NCI_Epithelial$pathway, GSEA_cogdxAPOE_34v33NCI_Epithelial$leadingEdge)
sn_Epi_paths <- list("Epi_Cilia"=unique(c(unlist(Temp[Temp[,1]=="GOCC_CILIUM",2], Temp[Temp[,1]=="GOBP_MICROTUBULE_BASED_MOVEMENT",2], Temp[Temp[,1]=="GOCC_DYNEIN_COMPLEX",2]))),
                     "Epi_Secretion"=unique(c(unlist(Temp[Temp[,1]=="GOBP_REGULATION_OF_PROTEIN_SECRETION",2], Temp[Temp[,1]=="GOBP_SECRETION",2], Temp[Temp[,1]=="GOBP_STEROL_TRANSPORT",2]))),
                     "Epi_Adhesion"=unique(c(unlist(Temp[Temp[,1]=="GOCC_ADHERENS_JUNCTION",2], Temp[Temp[,1]=="GOCC_TIGHT_JUNCTION",2], Temp[Temp[,1]=="GOMF_CADHERIN_BINDING",2]))),
                     "Epi_Lipid"=unique(c(unlist(Temp[Temp[,1]=="GOBP_LIPID_BIOSYNTHETIC_PROCESS",2], Temp[Temp[,1]=="GOCC_LIPID_DROPLET",2], Temp[Temp[,1]=="GOBP_LIPID_METABOLIC_PROCESS",2]))),
                     "Epi_Immune"=unique(c(unlist(Temp[Temp[,1]=="GOBP_REGULATION_OF_T_CELL_ACTIVATION",2], Temp[Temp[,1]=="GOBP_INFLAMMATORY_RESPONSE",2], Temp[Temp[,1]=="GOBP_RESPONSE_TO_WOUNDING",2])))
)
saveRDS(sn_Epi_paths, paste0(curr_dir, "/GSVA_modules/Exvivo_Epi_modules.rds"))

Temp <- cbind(GSEA_cogdxAPOE_34v33NCI_Fibroblast$pathway, GSEA_cogdxAPOE_34v33NCI_Fibroblast$leadingEdge)
sn_Fib_paths <- list("Fib_Barrier"=unique(c(unlist(Temp[Temp[,1]=="GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX",2], Temp[Temp[,1]=="GOBP_CELL_MATRIX_ADHESION",2]))), 
                     "Fib_transport"=unique(c(unlist(Temp[Temp[,1]=="GOMF_POTASSIUM_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",2]))), 
                     "Fib_Immune"=unique(c(unlist(Temp[Temp[,1]=="GOMF_IMMUNE_RECEPTOR_ACTIVITY",2], Temp[Temp[,1]=="GOBP_RESPONSE_TO_WOUNDING",2], Temp[Temp[,1]=="GOBP_INTERLEUKIN_6_PRODUCTION",2])))
)
saveRDS(sn_Fib_paths, paste0(curr_dir, "/GSVA_modules/Exvivo_Fib_modules.rds"))

Temp <- cbind(GSEA_cogdxAPOE_34v33NCI_Immune$pathway, GSEA_cogdxAPOE_34v33NCI_Immune$leadingEdge)
sn_Imm_paths <- list("Imm_Immune"=unique(c(unlist(Temp[Temp[,1]=="GOBP_MYELOID_LEUKOCYTE_ACTIVATION",2], Temp[Temp[,1]=="GOBP_T_CELL_PROLIFERATION",2], Temp[Temp[,1]=="GOBP_RESPONSE_TO_CYTOKINE",2]))),
                     "Imm_Lipid"=unique(c(unlist(Temp[Temp[,1]=="GOBP_SMALL_MOLECULE_METABOLIC_PROCESS",2], Temp[Temp[,1]=="GOBP_LIPID_LOCALIZATION",2], Temp[Temp[,1]=="GOBP_CHOLESTEROL_STORAGE",2])))
)
saveRDS(sn_Imm_paths, paste0(curr_dir, "/GSVA_modules/Exvivo_Imm_modules.rds"))


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


###Setup

library(purrr)
library(GSVA)
library(msigdbr)
library(stringr)

FDRcutoff <- 0.01
LFcutoff <- 0.25


##Dataprep

GSVA_fun(limma_cogdx_ADvNCI_Epithelial, Epi_modules)
GSVA_fun(limma_cogdx_ADvNCI_Fibroblast, Fib_modules)
GSVA_fun(limma_cogdx_ADvNCI_BAM, Imm_modules)


###Factors set up

mdat <- read_csv(paste0(curr_dir, "/metadata/ex vivo ChP metadata (analytical) newAPOE RNAonly.csv"))
mdat <- column_to_rownames(mdat, "projid")

#discrete var
msex = factor(paste(mdat$msex))
X10X_batch = factor(paste(mdat$X10X_batch))
library_batch = factor(paste(mdat$library_batch))
apoe = factor(paste(mdat$apoe_genotype))
apoe_sex = factor(paste(mdat$apoe_genotype, mdat$msex, sep="."))
cogdx = factor(paste(mdat$cogdx))
cogdx_sex = factor(paste(mdat$cogdx, mdat$msex, sep="."))
apoe_adnc <- factor(paste(mdat$apoe_genotype, mdat$ad_adnc, sep="."))
ad_adnc <- factor(paste(mdat$ad_adnc))

#continuous
cogdx_c = as.numeric(mdat$cogdx)
projid = as.numeric(mdat$projid)
age = as.numeric(mdat$age_death)
pmi = as.numeric(mdat$pmi)
educ = as.numeric(mdat$educ)
gpath = as.numeric(mdat$gpath)
Braaksc = as.numeric(mdat$braaksc)
tdp = as.numeric(mdat$tdp_st4)
tangles = as.numeric(mdat$tangsqrt_est_8reg)
amyl = as.numeric(mdat$amylsqrt_est_8reg)

arteriol = as.numeric(mdat$arteriol_scler)
caa = as.numeric(mdat$caa_4gp)
cvda = as.numeric(mdat$cvda_4gp2)
cogn_global = as.numeric(mdat$cogn_global_last)


####APOE cogdx ctrl

apoe <- factor(mdat$apoe_genotype, levels=c("33","34","other"))
cogdx_simple <- factor(mdat$cogdx_simple, levels=c("NCI","MCI","ADd"))
options(na.action='na.pass')
Maindm <- model.matrix(~apoe*cogdx_simple + msex + age + X10X_batch + pmi)
colnames(Maindm) <- make.names(colnames(Maindm))
rownames(Maindm) <- rownames(mdat)
colnames(Maindm)
Maindm <- Maindm[complete.cases(Maindm),]
Maincm <- makeContrasts(levels=Maindm, cogdxAPOE_34v33=apoe34, cogdxAPOE_ADvNCI=cogdx_simpleADd, cogdxAPOE_cogdxapoe=apoe34.cogdx_simpleADd, cogdxAPOE_ADvNCIin34=cogdx_simpleADd+apoe34.cogdx_simpleADd)

GSVA_limmafun_lane(GSVA_Epithelial)
GSVA_limmafun_lane(GSVA_Fibroblast)
GSVA_limmafun_lane(GSVA_BAM)


###Load limma Orgs

limma_GSVA_Org_Tsai_Epithelial <- read.csv("/Ex_ChP/Org_Tsai/limma/MainGenotype/Spreadsheet_GSVA_Epithelial_MainGenotype_uncorrected.csv", row.names=1)
limma_GSVA_Org_Tsai_Fibroblast <- read.csv("/Ex_ChP/Org_Tsai/limma/MainGenotype/Spreadsheet_GSVA_Fibroblast_MainGenotype_uncorrected.csv", row.names=1)
limma_GSVA_Org_Tsai_BAM <- read.csv("/Ex_ChP/Org_Tsai/limma/MainGenotype/Spreadsheet_GSVA_Immune_MainGenotype_uncorrected.csv", row.names=1)

limma_GSVA_Org_Jackson_Epithelial <- read.csv("/Ex_ChP/Org_Jackson/limma/MainGenotype/Spreadsheet_GSVA_Epithelial_MainGenotype_uncorrected.csv", row.names=1)
limma_GSVA_Org_Jackson_Fibroblast <- read.csv("/Ex_ChP/Org_Jackson/limma/MainGenotype/Spreadsheet_GSVA_Fibroblast_MainGenotype_uncorrected.csv", row.names=1)
limma_GSVA_Org_Jackson_BAM <- read.csv("/Ex_ChP/Org_Jackson/limma/MainGenotype/Spreadsheet_GSVA_Immune_MainGenotype_uncorrected.csv", row.names=1)

limma_GSVA_Org_Cocult_Epithelial <- read.csv("/Ex_ChP/Organoid_BAM2/MainGenotype/Spreadsheet_GSVA_Epithelial_MainGenotype_uncorrected.csv", row.names=1)
limma_GSVA_Org_Cocult_Fibroblast <- read.csv("/Ex_ChP/Organoid_BAM2/MainGenotype/Spreadsheet_GSVA_Fibroblast_MainGenotype_uncorrected.csv", row.names=1)
limma_GSVA_Org_Cocult_BAM <- read.csv("/Ex_ChP/Organoid_BAM2/MainGenotype/Spreadsheet_GSVA_Immune_MainGenotype_uncorrected.csv", row.names=1)


###Z score Custom scatterplot

library(ggplot2)
library(ggrepel)
library(ggpubr)

z_cutoff=.25

#Epithelial
GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_Epithelial, limma_GSVA_Org_Tsai_Epithelial)
GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_Epithelial, limma_GSVA_Org_Jackson_Epithelial)
GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_Epithelial, limma_GSVA_Org_Cocult_Epithelial)
GSVA_zScat_fun(limma_GSVA_Org_Tsai_Epithelial, limma_GSVA_Org_Jackson_Epithelial)
GSVA_zScat_fun(limma_GSVA_Org_Tsai_Epithelial, limma_GSVA_Org_Cocult_Epithelial)
GSVA_zScat_fun(limma_GSVA_Org_Jackson_Epithelial, limma_GSVA_Org_Cocult_Epithelial)
GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_Epithelial, limma_cogdxAPOE_ADvNCI_GSVA_Epithelial)
GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_Epithelial, limma_cogdxAPOE_ADvNCIin34_GSVA_Epithelial)
GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_Epithelial, limma_cogdxAPOE_cogdxapoe_GSVA_Epithelial)

GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_Epithelial, limma_cogdx_ADvNCI_GSVA_Epithelial)
GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_Epithelial, limma_cogdxAPOE_34v33ADd_GSVA_Epithelial)

#Fibroblast
GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_Fibroblast, limma_GSVA_Org_Tsai_Fibroblast)
GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_Fibroblast, limma_GSVA_Org_Jackson_Fibroblast)
GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_Fibroblast, limma_GSVA_Org_Cocult_Fibroblast)
GSVA_zScat_fun(limma_GSVA_Org_Tsai_Fibroblast, limma_GSVA_Org_Jackson_Fibroblast)
GSVA_zScat_fun(limma_GSVA_Org_Tsai_Fibroblast, limma_GSVA_Org_Cocult_Fibroblast)
GSVA_zScat_fun(limma_GSVA_Org_Jackson_Fibroblast, limma_GSVA_Org_Cocult_Fibroblast)
GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_Fibroblast, limma_cogdxAPOE_ADvNCI_GSVA_Fibroblast)
GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_Fibroblast, limma_cogdxAPOE_ADvNCIin34_GSVA_Fibroblast)
GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_Fibroblast, limma_cogdxAPOE_cogdxapoe_GSVA_Fibroblast)

GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_Fibroblast, limma_cogdx_ADvNCI_GSVA_Fibroblast)
GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_Fibroblast, limma_cogdxAPOE_34v33ADd_GSVA_Fibroblast)

#Immune
GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_BAM, limma_GSVA_Org_Tsai_BAM)
GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_BAM, limma_GSVA_Org_Jackson_BAM)
GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_BAM, limma_GSVA_Org_Cocult_BAM)
GSVA_zScat_fun(limma_GSVA_Org_Tsai_BAM, limma_GSVA_Org_Jackson_BAM)
GSVA_zScat_fun(limma_GSVA_Org_Tsai_BAM, limma_GSVA_Org_Cocult_BAM)
GSVA_zScat_fun(limma_GSVA_Org_Jackson_BAM, limma_GSVA_Org_Cocult_BAM)
GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_BAM, limma_cogdxAPOE_ADvNCI_GSVA_BAM)
GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_BAM, limma_cogdxAPOE_ADvNCIin34_GSVA_BAM)
GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_BAM, limma_cogdxAPOE_cogdxapoe_GSVA_BAM)

GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_BAM, limma_cogdx_ADvNCI_GSVA_BAM)
GSVA_zScat_fun(limma_cogdxAPOE_34v33NCI_GSVA_BAM, limma_cogdxAPOE_34v33ADd_GSVA_BAM)


