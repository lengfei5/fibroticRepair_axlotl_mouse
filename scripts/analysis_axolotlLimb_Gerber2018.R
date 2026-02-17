##########################################################################
##########################################################################
# Project: axolotl limb scRNA-seq dada reanalysis
# Script purpose: Gerber2018 scRNA-seq from Fluidigm.C1
# Usage example: The analysis is using Seurat v4
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Feb 16 13:30:53 2026
##########################################################################
##########################################################################
rm(list = ls())

library(Seurat)
#library(decoupleR)
library(tictoc)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(data.table)
#library("viridis")

version.analysis = '_axolotl_20250829'
resDir = paste0("../results/scRNAseq_analysis_immune", version.analysis, '/')
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '/groups/tanaka/People/current/jiwang/projects/bone_healing_CSD/fromTobie/CSD_batch1_batch2/'

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts'
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_scRNAseq.R')
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_Visium.R')

annot = readRDS(paste0('/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                       'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))

convert_to_geneSymbol = function(gene.ids, annot)
{
  # gene.ids = rownames(mat)
  mm = match(gene.ids, annot$geneID)
  ggs = annot$gene.symbol.toUse[mm]
  
  kk = which(is.na(ggs))
  ggs[kk] = gene.ids[kk]
  #return(make.unique(ggs))
  return(ggs)
  
}

#### TFs and  gene example of signaling pathways 
tfs = readRDS(file = paste0('/groups/tanaka/People/current/jiwang/projects/RA_competence/data',
                            '/annotations/curated_human_TFs_Lambert.rds'))
tfs = unique(tfs$`HGNC symbol`)

sps = readRDS(file = paste0("/groups/tanaka/People/current/jiwang/projects/RA_competence/",
                            '/data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
sps = unique(sps$gene)
sps = toupper(sps)


##########################################
# import the scRNA-seq data
##########################################
aa = readRDS(file = paste0("/groups/tanaka/People/current/jiwang/projects/scRNA_regeneration/data/", 
                           '/Gerber2018_Fluidigm.C1_batches_seuratObj.rds'))


aa = subset(aa, cells = colnames(aa)[which(aa$condition != '1apa')])
aa$condition = droplevels(aa$condition)

table(aa$sample)

#ggs = rownames(aa)

mm = match(rownames(aa), annot$geneID)
#ggs = paste0(annot$gene.symbol.toUse[mm], '_',  annot$geneID[mm])
#ggs[is.na(mm)] = rownames(aa)[is.na(mm)]
#ggs = get_geneName(ggs)
ggs = annot$gene.symbol.toUse[mm]

which(ggs == "KAZALD1")

Idents(aa) = factor(aa$condition, levels = c('0dpa', '3dpa', '5dpa', '8dpa', '11dpa', '18dpa',
                                             'Stage40', 'Stage44'))

DotPlot(aa, features = rownames(aa)[which(ggs == "KAZALD1")[2]]) + RotatedAxis()


xx = aa@assays$RNA@data[which(ggs == "KAZALD1")[2], ]
xx = data.frame(xx, colnames(aa), aa$condition)

colnames(xx) = c('logNorm_expression', 'cellID', 'axolotl_blastema_timePoint')


write.csv2(xx, file = paste0(resDir, 'axolotlLimb_Gerber2018_geneExpression_Kazald1.csv'), 
           quote = FALSE, row.names = TRUE)


##########################################
# import the brain data from Lust 2022
# some codes from Tomas
##########################################
meta = read.csv(file = paste0('/groups/tanaka/People/current/jiwang/projects/',
                'scRNA_regeneration/data/Lust2022_brain/6390083/pallium_metadata_simp.csv'),  
                header = T, row.names = 1)


cols_cc = c(
  #epen
  "#12400c", "#2d6624","#1d4f15", "#174711", "#2d6624", "#3d7f33", "#3b7b30", 
  "#468b3b", "#4f9843","#5dae50", "#66bb58", "#72cd64", "#306a26", "#78d669", "#81e472",
  #gaba
  "#700209", "#75090e","#7a0f13", "#801517", "#851a1b", "#8a1f1f", "#902423", "#952927", 
  "#9a2d2c","#a03230", "#a53634", "#aa3a39", "#b03f3d","#b54342", "#ba4846", "#c04c4b", "#c5504f",
  "#ca5554", "#d05959", "#d55e5e","#73050c", "#780c11","#8d2221", "#982b2a","#a23432", "#a83837", 
  "#b2413f", "#b84544", "#bd4a49", "#c85352", #"#cd5756",
  #glut
  "#054674", "#134d7b","#1d5481", "#265a88", "#2e618e", "#73a4cb", "#366995", "#3e709c", "#4677a2",
  "#4d7ea9", "#5586b0", "#5c8db7", "#6495bd","#6b9cc4", "#7bacd2", "#8ebfe4", "#96c7eb", "#9ecff2",
  "#18507e", "#18507e","#2a5e8b", "#497ba6","#5889b3", "#6fa0c8","#7fafd6", "#6091ba", "#5182ac",
  "#3a6c98",
  "#a6d7f9",
  #npc
  "#ffb120", "#feb72a","#fdbc34", "#fcc13d", "#fbc745", "#facc4e", "#f9d156", "#f8d65f", 
  "#f8da68","#f7df70", "#f7e479", "#f7e882", "#f7ed8a", "#f7f193", "#eca319"
)

ccnames = unique(sort(meta$cellclusters))
names(cols_cc) = c(ccnames[grepl("epen", ccnames)], ccnames[grepl("GABA", ccnames)], 
                   ccnames[grepl("glut", ccnames)],ccnames[grepl("npc", ccnames)])

cols_cc = c(cols_cc, "microglia_8" = "#E6530D", 
            "oligodendrocyte_15" = "#E43D88", "oligodendrocyte_10" = "#F662A5",
            "endothelial_11" = "#712166", "endothelial_12" = "#B0279D", 
            "endothelial_14" = "#BE5AB0")

reg_cols = c("other/unknown_pred" = "#C7CCC7", 
             "medial" = "#52168D", "medial_pred" = "#661CB0", 
             "dorsal" = "#C56007", "dorsal_pred" = "#ED7307", 
             "lateral" = "#118392", "lateral_pred" = "#16A3B6")
reg_cols_simp = c("medial" = "#52168D", "dorsal" = "#C56007", "lateral" = "#118392")

ep_wpi_srat = readRDS(file = paste0('/groups/tanaka/People/current/jiwang/projects/',
                                    'scRNA_regeneration/data/Lust2022_brain/6390083/ep_wpi_srat.rds'))

ep_wpi_srat$ep_type = ep_wpi_srat$ctlabels
ep_wpi_srat$ep_type[ep_wpi_srat$ep_type %in% c("epen_clus_4", "epen_clus_3")] = "active ependymal"
ep_wpi_srat$ep_type[ep_wpi_srat$ep_type %in% c("epen_clus_1", "epen_clus_7",
                                               "epen_clus_14", "epen_clus_13")] = "pro-neuronal ependymal"
ep_wpi_srat$ep_type[ep_wpi_srat$ep_type %in% c("epen_clus_0", "epen_clus_2", "epen_clus_12",
                                               "epen_clus_8", "epen_clus_6", "epen_clus_5",
                                               "epen_clus_9", "epen_clus_10", "epen_clus_11")] = "quiescent ependymal"

col_list = list(ctlabels = cols_cc, 
                ep_type = unname(cols_cc[grepl("epen", names(cols_cc))][c(1,15,8)]),
                Phase = MetBrewer::met.brewer(name="Egypt", n=3, type = "discrete"),
                integEp_res.1.5 = MetBrewer::met.brewer(name="Klimt",
                                                        n=length(unique(ep_wpi_srat@meta.data$integEp_res.1.5)),
                                                        type = "continuous"),
                reglabels = reg_cols)
csfunc = colorRampPalette(c("#df16df", "#ffd9f9"))
col_list$sample_simp = c("gray80", csfunc(6))
names(col_list$sample_simp) = c("SS", paste0(c(1,2,4,6,8,12), "_wpi"))


plt_umap_l = list()
for(v in c("integEp_res.1.5", "ctlabels", "ep_type", "Phase", "sample_simp", "reglabels")){
  # v = "integEp_res.1.5"
  pal_plt = col_list[[v]]
  
  plot_df = Embeddings(ep_wpi_srat, reduction = "umap_harmony")
  plot_df = merge(plot_df, ep_wpi_srat@meta.data, by = 0)
  plot_df = plot_df[order(plot_df$sample_simp, decreasing = F),]
  
  lab_df = plot_df[,c("umap_harmony_1","umap_harmony_2",v)]
  c1 = tapply(lab_df$umap_harmony_1, lab_df[,v], median)
  c2 = tapply(lab_df$umap_harmony_2, lab_df[,v], median)
  lab_df = data.frame(c1, c2, names(c1))
  colnames(lab_df) = c("umap_harmony_1","umap_harmony_2",v)
  
  plt_umap_l[[v]] = ggplot() + 
    geom_point(data = plot_df, 
               mapping = aes_string(x = "umap_harmony_1", y = "umap_harmony_2", colour = v),
               size = 0.7)+
    geom_label(data = lab_df, 
               mapping = aes_string(x = "umap_harmony_1", y = "umap_harmony_2", label = v),
               size = 2.7, label.padding = unit(0.1, "lines"))+
    scale_colour_manual(values = pal_plt)+
    theme_void()+
    theme(plot.title = element_blank(),
          legend.position = "none",
          aspect.ratio = 1, 
          axis.title = element_blank(), 
          axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          axis.line = element_blank())
  
}

for(n in names(plt_umap_l)){
  pdf(paste0("results/Div-seq/ep_wpi_srat_", n, ".pdf"), height = 4, width = 4)
  print(plt_umap_l[[n]])
  dev.off()
}

pdf(paste0("results/Div-seq/ep_wpi_srat_", v, "_legend.pdf"), height = 1, width = 0.7)
cowplot::plot_grid(cowplot::get_legend(plt))
dev.off()

eps = ep_wpi_srat
DimPlot(eps, reduction = 'umap_harmony', group.by = 'clusters_used', label = TRUE, repel = TRUE)

kk = which(rownames(eps) == 'KAZALD1')
DotPlot(eps, features =  "KAZALD1", group.by = 'sample_simp') + RotatedAxis()

xx = eps@assays$RNA@data[which(rownames(eps) == "KAZALD1"), ]
xx = data.frame(xx, colnames(eps), eps$sample_simp,  eps$clusters_used)

colnames(xx) = c('logNorm_expression', 'cellID', 'axolotl_Epen_timePoint', 'axolotl_Epen_cluster')


write.csv2(xx, file = paste0(resDir, 'axolotlBrain_Lust2022_geneExpression_Kazald1.csv'), 
           quote = FALSE, row.names = TRUE)
