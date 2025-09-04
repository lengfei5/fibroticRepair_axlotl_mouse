##########################################################################
##########################################################################
# Project: Collaboration with Sabine and Sebastian from Cologne
# Script purpose: analyze macrophage subtypes from Natasya's data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Jan 13 11:01:15 2023
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
library("viridis")

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


## color schema 
levels = c('CSD_0dpa', 
           'BL_3and5dpa', 'BL_5dpa', 'BL_6dpa', 'BL_7dpa',  'BL_8dpa', 'BL_11dpa', 
           'CSD_3dpa', 'CSD_5dpa', 'CSD_6dpa', 'CSD_7dpa', 'CSD_8dpa', 'CSD_11dpa')

cols = rep(NA, length = length(levels))
names(cols) = levels

cols[1] = 'gray60'
#cols[1:3] = viridis(3)
cols[2:7] = colorRampPalette((brewer.pal(n = 6, name ="Blues")))(6)
cols[8:13] = colorRampPalette((brewer.pal(n = 6, name ="OrRd")))(6)


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

##########################################
# subset batch 1 day3, 8, 11 
##########################################
aa = subset(aa, cells = colnames(aa)[which(aa$batch == 'batch1')])

aa = NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 8000) # find subset-specific HVGs

aa <- ScaleData(aa)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)

ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, reduction = "pca", dims = 1:50, n.neighbors = 50,  min.dist = 0.3)

p1 = DimPlot(aa, group.by = 'time', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(resDir, '/Tobie_batch1Data_umap_axloltol_BL_celltypes.pdf'), 
       width = 14, height = 6)


saveRDS(aa, file = paste0(RdataDir, 
                          '/axoltol_limbBlatema_batch1_fromTobi_filterCelltypes_geneNames.rds'))

########################################################
########################################################
# Section II: check marker genes for 
# - don't eat me pathway
# - Pan macrophage 
# - Senescence markers
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, 
                           '/axoltol_limbBlatema_batch1_fromTobi_filterCelltypes_geneNames.rds'))

ggs = rownames(aa)
ggs = get_geneName(ggs)

p1 = DimPlot(aa, group.by = 'time', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE)

p1 + p2


##########################################
# subset macrophage and subclustering
##########################################
xx = readxl::read_xlsx(paste0('../../bone_healing_CSD/Cell_marker_Mouse.xlsx'), 
                       sheet = 1)
kk = grep('macrophage|Macrophage', xx$cell_name)
xx = xx[kk, ]
xx = data.frame(xx)

write.table(xx, file = '../data/Cell_marker_Mouse_macrophages.txt', 
            sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)


subs = subset(aa, cells = colnames(aa)[which(aa$celltype == "Macrophages")])

subs = NormalizeData(subs, normalization.method = "LogNormalize", scale.factor = 10000)
subs <- FindVariableFeatures(subs, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs

subs <- ScaleData(subs, features = rownames(subs))
subs <- RunPCA(subs, features = VariableFeatures(object = subs), verbose = FALSE, weight.by.var = TRUE)

subs <- RunUMAP(subs, reduction = "pca", dims = 1:50, n.neighbors = 50,  min.dist = 0.3)

p1 = DimPlot(subs, group.by = 'time', label = TRUE, repel = TRUE)
p2 = DimPlot(subs, group.by = 'Phase', label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(resDir, '/Tobie_batch1Data_umap_axloltol_BL_celltypes_macrophage.pdf'), 
       width = 14, height = 6)

ElbowPlot(subs, ndims = 30)

subs <- FindNeighbors(subs, dims = 1:20)
subs <- FindClusters(subs, verbose = FALSE, algorithm = 3, resolution = 0.7)

DimPlot(subs, label = TRUE, repel = TRUE)

ggsave(filename = paste0(resDir, '/Tobie_batch1Data_umap_axloltol_BL_celltypes_macrophage_subclustering.pdf'), 
       width = 10, height = 6)

subs$clusters = subs$seurat_clusters

Idents(subs) = factor(subs$clusters)
oupMarker <- FindAllMarkers(subs)

oupMarker %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1.0) %>%
  slice_head(n = 15) %>%
  ungroup() -> top10
DoHeatmap(subs, features = top10$gene) + NoLegend()


ggsave(filename = paste0(resDir, 
                         '/Tobie_batch1Data_umap_axloltol_BL_celltypes_macrophage_subclustering_top10Markers.pdf'), 
       width = 8, height = 16)

oupMarker <- data.table(oupMarker)
oupMarker$pct.diff = oupMarker$pct.1 - oupMarker$pct.2
oupMarker <- oupMarker[, c("cluster","gene","avg_log2FC","pct.1","pct.2",
                           "pct.diff","p_val","p_val_adj")]
#fwrite(oupMarker, sep = "\t", file = "images/clustMarkers.txt")
#aa@misc$marker <- oupMarker      # Store markers into Seurat object

# Get top genes for each cluster and do dotplot / violin plot
#oupMarker$cluster = factor(oupMarker$cluster, levels = reorderCluster)
oupMarker = oupMarker[order(cluster, -avg_log2FC)]
#genes.to.plot <- knownGenes

#oupMarker = oupMarker[which(oupMarker$cluster == 'Macrophages'), ]

markers = unique(c(xx$Symbol, xx$marker))
markers = markers[!is.na(markers)]
markers = toupper(markers)

mm = match(get_geneName(oupMarker$gene), markers)

oupMarker = oupMarker[which(!is.na(mm)), ]

write.table(oupMarker, file = paste0(resDir, '/macrophage_markers_detected_inAxolotl.txt'), 
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)


oupMarker %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1.0) %>%
  slice_head(n = 15) %>%
  ungroup() -> top10
DoHeatmap(subs, features = top10$gene) + NoLegend()

ggsave(filename = paste0(resDir, 
                         '/Tobie_batch1Data_umap_axloltol_BL_celltypes_macrophage_subclustering_',
                         'top10_annotatedMouseMarkers.pdf'), 
       width = 8, height = 16)


########################################################
########################################################
# Section III: 
# subclustering macrophage and fibroblast for Elly and Sabine's grant application
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, 
                           '/axoltol_limbBlatema_batch1_fromTobi_filterCelltypes_geneNames.rds'))

ggs = rownames(aa)
ggs = get_geneName(ggs)

p1 = DimPlot(aa, group.by = 'time', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE)

p1 + p2

aa$celltype[which(aa$celltype == 'Connective Tissue')] = 'Fibroblasts'
aa <- subset(aa, subset = celltype == 'Fibroblasts'| celltype == 'Macrophages')

aa = NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs

aa <- ScaleData(aa)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)

ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, reduction = "pca", dims = 1:20, n.neighbors = 30,  min.dist = 0.3)

p1 = DimPlot(aa, group.by = 'time', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE)

p1 + p2


##########################################
# quickly check the macrophage subtype markers 
##########################################
subs = subset(aa, cells = colnames(aa)[which(aa$celltype == "Macrophages")])

subs = NormalizeData(subs, normalization.method = "LogNormalize", scale.factor = 10000)
subs <- FindVariableFeatures(subs, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs

subs <- ScaleData(subs, features = rownames(subs))
subs <- RunPCA(subs, features = VariableFeatures(object = subs), verbose = FALSE, weight.by.var = TRUE)

ElbowPlot(subs, ndims = 50)

subs <- FindNeighbors(subs, dims = 1:20)
subs <- FindClusters(subs, verbose = FALSE, algorithm = 3, resolution = 0.5)

subs <- RunUMAP(subs, reduction = "pca", dims = 1:30, n.neighbors = 30,  min.dist = 0.3)

p1 = DimPlot(subs, group.by = 'time', label = TRUE, repel = TRUE)
p2 = DimPlot(subs, group.by = 'seurat_clusters', label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(resDir, '/Tobie_batch1Data_axloltolBlastema_macrophage_subclusters.pdf'), 
       width = 14, height = 6)


subs$clusters = subs$RNA_snn_res.0.5
DimPlot(subs, group.by = 'Phase')

subs$subtypes = NA
subs$subtypes[which(subs$clusters == 7)] = 'cycling'

# https://www.biocompare.com/Editorial-Articles/566347-A-Guide-to-Macrophage-Markers/
markers = c('ADGRE1', 'CD14', 'CD68', 'CX3CR1', 'ITGAM', 'CSF1R', ' H2AB1', 'MERTK',# pan macrophage
            'TLR2', 'NOS2', 'CD80', 'CD86', 'IFNG', # M1
             'ARG1', 'CD163', 'IL4', 'IRF4' # M2
            )

mm = match(markers, ggs)
mm = mm[which(!is.na(mm))]

FeaturePlot(subs, features = rownames(aa)[mm]) 
  #&
  #scale_color_distiller(palette = "RdYlBu")
  #scale_color_viridis_c()
ggsave(filename = paste0(resDir, '/Tobie_umap.harmony_axloltol_limbBlatema_someGeneMarkers.pdf'), 
       width = 12, height = 8)


Idents(subs) = factor(subs$clusters)
oupMarker <- FindAllMarkers(subs)

oupMarker = oupMarker[grep('^AME', oupMarker$gene, invert = TRUE), ]

oupMarker %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1.0) %>%
  slice_head(n = 15) %>%
  ungroup() -> top10
DoHeatmap(subs, features = top10$gene) + NoLegend()

ggsave(filename = paste0(resDir, 
                         '/Tobie_batch1Data_umap_axloltol_BL_celltypes_macrophage_subclustering_',
                         'top15_annotatedMouseMarkers.pdf'), 
       width = 12, height = 20)

kk = which(subs$clusters != 7)
subs$subtypes[kk] = paste0('C', subs$clusters[kk], '.ax')

DimPlot(subs, group.by = 'subtypes', label = TRUE, repel = TRUE)

saveRDS(subs, file = paste0(RdataDir, 
                          '/axoltol_limbBlatema_batch1_macrophage_subtypes.rds'))

