##########################################################################
##########################################################################
# Project: Collaboration with Sabine and Sebastian from Cologne
# Script purpose: analyze the 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Jan 13 11:01:15 2023
##########################################################################
##########################################################################
rm(list = ls())

library(Seurat)
library(decoupleR)
library(tictoc)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library("viridis")

version.analysis = '_axolotl_20250702'
resDir = paste0("../results/scRNAseq_analysis_immune", version.analysis)
RdataDir = paste0(resDir, '/Rdata')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../fromTobie/CSD_batch1_batch2/'

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

levels = c('CSD_0dpa', 
           'BL_3and5dpa', 'BL_5dpa', 'BL_6dpa', 'BL_7dpa',  'BL_8dpa', 'BL_11dpa', 
           'CSD_3dpa', 'CSD_5dpa', 'CSD_6dpa', 'CSD_7dpa', 'CSD_8dpa', 'CSD_11dpa')

cols = rep(NA, length = length(levels))
names(cols) = levels

cols[1] = 'gray60'
#cols[1:3] = viridis(3)
cols[2:7] = colorRampPalette((brewer.pal(n = 6, name ="Blues")))(6)
cols[8:13] = colorRampPalette((brewer.pal(n = 6, name ="OrRd")))(6)

## gene example of signaling pathways 
sps = readRDS(file = paste0("/groups/tanaka/People/current/jiwang/projects/RA_competence/",
                            '/data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
sps = unique(sps$gene)
sps = toupper(sps)

## TFs 
tfs = readRDS(file = paste0('/groups/tanaka/People/current/jiwang/projects/RA_competence/data',
                            '/annotations/curated_human_TFs_Lambert.rds'))
tfs = unique(tfs$`HGNC symbol`)


########################################################
########################################################
# Section 0: use the scRNA-seq samples from Natasya and Tobie  
# import the processed scRNA data from Natasia and Tobi 
########################################################
########################################################
aa = readRDS(file = paste0(dataDir, 'BL_SeuratObj.RDS'))
#xx = readRDS(file = paste0("../fromTobie/CSD_batch1/", 'BL_SeuratObj.RDS'))

aa$batch = aa$exp
aa$time = paste0('dpa', aa$time)
aa$time = factor(aa$time, levels = c('dpa3', 'dpa5', 'dpa6', 'dpa7', 'dpa8', 'dpa11'))

p1 = DimPlot(aa, group.by = 'time', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'batch', label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(resDir, '/Tobie_umap.harmony_BL.pdf'), width = 18, height = 6)

aa$condition = aa$orig.ident
aa$condition[which(aa$condition == 'BL_3_5dpa')] = 'BL_3and5dpa'

aa$sample = aa$condition
aa$days = aa$condition
aa$days = gsub('BL_', '', aa$days)

p1 = DimPlot(aa, group.by = 'time', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(resDir, '/Tobie_umap.harmony_axloltol_BL_celltypes.pdf'), width = 18, height = 6)

saveRDS(aa, file = paste0(RdataDir, '/axoltol_limb_Blatema_twoBacthes_harmonyMerged_fromTobi.rds'))


## subsetting cell types excluding blood cells
Idents(aa) = factor(aa$celltype)
xx = subset(aa, idents = c('Connective Tissue', 'Epidermis', 'Macrophages', 'Neutrophils',
                           "B cells", 'Endothelial cells', "Eosinophils/Killer cells", "T cells"))


aa = xx 
rm(xx)


# RenameGenesSeurat  ------------------------------------------------------------------------------------
# RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = HGNC.updated[[i]]$Suggested.Symbol) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
#   print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
#   RNA <- obj@assays$RNA
#   
#   if (nrow(RNA) == length(newnames)) {
#     if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
#     if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
#     if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
#   } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
#   obj@assays$RNA <- RNA
#   return(obj)
# }

ggs = convert_to_geneSymbol(rownames(aa), annot = annot)
ids = rownames(aa)
kk = which(ggs != ids)
for(i in kk) ggs[i] = paste0(ggs[i], '_', ids[i])

counts = aa@assays$RNA@counts
meta = aa@meta.data
rds = aa@reductions$harmony@cell.embeddings

rownames(counts) = ggs

xx = CreateSeuratObject(counts = counts, meta.data = meta)
xx[["harmony"]] <- CreateDimReducObject(embeddings = rds, key = "harmony_", assay = DefaultAssay(xx))
#xx = RenameGenesSeurat(obj = aa, newnames = ggs)

aa = xx
rm(xx)


#Cluster cells
aa <- RunUMAP(aa, reduction = "harmony", dims = 1:50, n.neighbors = 50,  min.dist = 0.3)

p1 = DimPlot(aa, group.by = 'time', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(resDir, '/Tobie_umap.harmony_axloltol_BL_celltypes_seletedCellTypes.pdf'), 
       width = 16, height = 6)

#aa <- FindNeighbors(aa, dims = 1:100, reduction = "harmony")
#CSD_CT <- FindClusters(CSD_CT, resolution = 1)

saveRDS(aa, file = paste0(RdataDir, 
                          '/axoltol_limb_Blatema_twoBacthes_harmonyMerged_fromTobi_',
                          'filterCelltypes_geneNames.rds'))


########################################################
########################################################
# Section II: check marker genes for 
# - don't eat me pathway
# - Pan macrophage 
# - Senescence markers
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, 
                           '/axoltol_limb_Blatema_twoBacthes_harmonyMerged_fromTobi_',
                           'filterCelltypes_geneNames.rds'))

aa = NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 8000) # find subset-specific HVGs

aa <- ScaleData(aa)


ggs = rownames(aa)
ggs = get_geneName(ggs)


##########################################
# "don't eat me" pathway
##########################################
markers = c('CD47', "MER6","IAP", 'SIRPA', 'CD24', 'QPTC', 'QPTCL', 'B2M', 'AXL', "HAVCR2", "TIMD4",
            'CD274', 'PDCD1')

mm = match(markers, ggs)
mm = mm[which(!is.na(mm))]

FeaturePlot(aa, features = rownames(aa)[mm]) 
  #&
  #scale_color_distiller(palette = "RdYlBu")
  #scale_color_viridis_c()
ggsave(filename = paste0(resDir, '/Tobie_umap.harmony_axloltol_limbBlatema_donotEatMePathway.pdf'), 
       width = 8, height = 12)


##########################################
# senescence markers 
##########################################
markers = c('TP53', 'CDKN1A', "CDKN2A", "LMNB1", 
            "RB1", 'TP53BP1', 'MKI67')
mm = match(markers, ggs)
mm = mm[which(!is.na(mm))]

pdf(paste0(resDir, '/Tobie_umap.harmony_axloltol_limbBlatema_senescenceMarkers_time.pdf'), 
    height = 6, width =20, useDingbats = FALSE)
for(i in mm)
{
  p1 = FeaturePlot(aa, features = rownames(aa)[i], split.by = 'time') 
  plot(p1)
  
}

dev.off()

ggsave(filename = paste0(resDir, '/Tobie_umap.harmony_axloltol_limbBlatema_senescenceMarkers.pdf'), 
       width = 8, height = 12)

##########################################
# macrophage markers
##########################################
Idents(aa) = factor(aa$celltype)
oupMarker <- FindAllMarkers(aa)
oupMarker <- data.table(oupMarker)
oupMarker$pct.diff = oupMarker$pct.1 - oupMarker$pct.2
oupMarker <- oupMarker[, c("cluster","gene","avg_log2FC","pct.1","pct.2",
                           "pct.diff","p_val","p_val_adj")]
#fwrite(oupMarker, sep = "\t", file = "images/clustMarkers.txt")
aa@misc$marker <- oupMarker      # Store markers into Seurat object

# Get top genes for each cluster and do dotplot / violin plot
#oupMarker$cluster = factor(oupMarker$cluster, levels = reorderCluster)
oupMarker = oupMarker[order(cluster, -avg_log2FC)]
#genes.to.plot <- knownGenes

oupMarker = oupMarker[which(oupMarker$cluster == 'Macrophages'), ]

xx = readxl::read_xlsx('../Cell_marker_Mouse.xlsx', sheet = 1)
kk = grep('macrophage|Macrophage', xx$cell_name)
xx = xx[kk, ]
xx = data.frame(xx)

markers = unique(c(xx$Symbol, xx$marker))
markers = markers[!is.na(markers)]
markers = toupper(markers)

mm = match(get_geneName(oupMarker$gene), markers)

oupMarker = oupMarker[which(!is.na(mm)), ]

write.table(oupMarker, file = paste0(resDir, '/macrophage_markers_detected_inAxolotl.txt'), 
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

DotPlot(aa, features = oupMarker$gene[1:50], 
        group.by = "celltype"
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  RotatedAxis() + 
  coord_flip() 

ggsave(filename = paste0(resDir, '/macrophage_markers_detected_inAxolotl_top50.pdf'), 
       width = 8, height = 16)

