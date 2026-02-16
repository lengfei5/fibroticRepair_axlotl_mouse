##########################################################################
##########################################################################
# Project: Collaboration with Sabine and Sebastian from Cologne
# Script purpose: analyze macrophage subtypes from Natasya's data
# Usage example: This script is running Seurat_5.0.2 in R v4.3.0
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

saveRDS(aa, file = paste0(RdataDir, '/axoltol_limb_Blatema_twoBacthes_harmonyMerged_fromTobi_seuratV4.rds'))


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
# double check some marker genes of RUNX for Elly 
##########################################
aa = readRDS(paste0(RdataDir, 
                    '/axoltol_limb_Blatema_twoBacthes_harmonyMerged_fromTobi_',
                    'filterCelltypes_geneNames.rds'))

aa = NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 8000) # find subset-specific HVGs

aa <- ScaleData(aa)

genes = c(rownames(aa)[grep('RUNX1-|CBFB|ZEB1|SNAI', rownames(aa))])

FeaturePlot(aa, features = genes)

p1 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE) + NoLegend()
p11 = DimPlot(aa, group.by = 'time', label = TRUE, repel = TRUE) 
p2 = FeaturePlot(aa, features = genes)

(p1/p11) + p2 

ggsave(filename = paste0(resDir, '/Tobie_umap.harmony_axloltol_BL_celltypes_RunxGenes.pdf'), 
       width = 16, height = 8)

VlnPlot(aa, features = genes, group.by = 'time', split.by = 'celltype')

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
                          '/axoltol_limbBlatema_batch1_fromTobi_filterCelltypes_geneNames_seuratV4.rds'))

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



##########################################
# doublet searching  
##########################################
aa = readRDS(file = paste0(RdataDir, 
                           '/axoltol_limbBlatema_batch1_fromTobi_filterCelltypes_geneNames_seuratV4.rds'))

aa = NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)

aa$DF_out = NA

library(DoubletFinder)

aa$condition = factor(aa$condition)
Idents(aa) = aa$condition

cc = unique(aa$condition)

for(n in 1:length(cc))
{
  # n = 1
  cat(n, '-----', as.character(cc[n]), '\n')
  subs <- subset(aa, condition == cc[n])
  
  subs <- FindVariableFeatures(subs, selection.method = "vst", nfeatures = 5000)
  subs <- ScaleData(subs)
  
  subs <- RunPCA(subs, verbose = TRUE)
  subs <- FindNeighbors(subs, dims = 1:30)
  subs <- FindClusters(subs, resolution = 0.5)
  
  subs <- RunUMAP(subs, dims = 1:30)
  
  sweep.res.list_nsclc <- paramSweep_v3(subs)
  sweep.stats_nsclc <- summarizeSweep(sweep.res.list_nsclc, GT = FALSE)
  bcmvn_nsclc <- find.pK(sweep.stats_nsclc)
  
  pK <- bcmvn_nsclc %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    filter(BCmetric == max(BCmetric)) %>%
    select(pK) 
  
  pK <- as.numeric(as.character(pK[[1]]))
  annotations <- subs@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  
  
  nExp_poi <- round(0.076*nrow(subs@meta.data))  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  subs <- DoubletFinder::doubletFinder_v3(subs, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj,  
                                          reuse.pANN = FALSE , sct = FALSE)
  
  df_out = subs@meta.data
  subs$DF_out = df_out[, grep('DF.classification', colnames(df_out))]
  
  aa$DF_out[match(colnames(subs), colnames(aa))] = subs$DF_out
  
  DimPlot(subs, label = TRUE, repel = TRUE, group.by = 'DF_out',
          raster=FALSE)
  ggsave(filename = paste0(resDir, '/subs_doubletFinder_out_', cc[n], version.analysis, '.pdf'), 
         width = 12, height = 8)
  
  saveRDS(subs, file = paste0(RdataDir, 'subs_doubletFinder_out_', cc[n], version.analysis,  
                              '.rds'))
  
}

DimPlot(aa, group.by = 'DF_out')
saveRDS(aa, file = paste0(RdataDir, 
                          '/axoltol_limbBlatema_batch1_fromTobi_filterCelltypes_geneNames_DFout_seuratV4.rds'))


########################################################
########################################################
# Section III: 
# subclustering macrophage and fibroblast for Elly and Sabine's grant application
########################################################
########################################################
#aa = readRDS(file = paste0(RdataDir, 
#                           '/axoltol_limbBlatema_batch1_fromTobi_filterCelltypes_geneNames_seuratV4.rds'))
#aa = readRDS(file = paste0(RdataDir, 
#                                   '/axoltol_limbBlatema_batch1_DFout_filterCelltypes_geneNames_',
#                                   'FB_Mphg_subtypes_v2.rds'))
aa = readRDS(file = paste0(RdataDir, 
                                   '/axoltol_limbBlatema_batch1_DFout_filterCelltypes_geneNames_',
                                   'FB_Mphg.cleaned_subtypes_v3.rds'))
ggs = rownames(aa)
ggs = get_geneName(ggs)

which(ggs == "KAZALD1")[2]
Idents(aa) = factor(aa$condition, levels = c("BL_3and5dpa", "BL_8dpa", "BL_11dpa"))
DotPlot(aa, features = rownames(aa)[which(ggs == "KAZALD1")[2]]
       ) + RotatedAxis()


xx = aa@assays$RNA$data[which(ggs == "KAZALD1")[2], ]
xx = data.frame(xx, colnames(aa), aa$condition)
colnames(xx) = c('logNorm_expression', 'cellID', 'axolotl_blastema_timePoint')


write.csv2(xx, file = paste0(resDir, 'axolotlLimb_blastema_geneExpression_Kazald1.csv'), 
           quote = FALSE, row.names = TRUE)

# ## change the annotation the M5 macrophage to epidermis  
# kk = which(aa$subtypes == 'FB7.PITX2+.SIX1+')
# 
# aa$celltype[kk] = 'Muscle'
# aa$cluster[kk] = 'Muscle'
# aa$subtypes[kk] = 'Muscle'

# saveRDS(aa, file = paste0(RdataDir,
#                           '/axoltol_limbBlatema_batch1_DFout_filterCelltypes_geneNames_',
#                           'FB.clean_Mphg.cleaned_subtypes_v3.rds'))

p1 = DimPlot(aa, group.by = 'time', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE)

p1 + p2

Running_UMAP = FALSE
if(Running_UMAP){
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 8000) # find subset-specific HVGs
  aa <- ScaleData(aa)
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
  
  ElbowPlot(aa, ndims = 50)
  
  aa <- RunUMAP(aa, reduction = "pca", dims = 1:50, n.neighbors = 100,  min.dist = 0.4)
  
  aa$celltype[which(aa$celltype == 'Connective Tissue')] = 'FB'
  
  p1 = DimPlot(aa, group.by = 'time', label = TRUE, repel = TRUE)
  p2 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE)
  
  p1 + p2
  
  ggsave(filename = paste0(resDir, '/batch1Data_umap_axloltol_BL_celltypes_.pdf'), 
         width = 16, height = 6)
  
  saveRDS(aa, file = paste0(RdataDir, 
                            '/axoltol_limbBlatema_batch1_fromTobi_filterCelltypes_geneNames_umapUsed.rds'))
  
}


##########################################
# double check the macrophage and FB cell types markers 
##########################################
markers =  c('Procr', 'Dpt', 'Pi16', 'Col1a2', 'Acta2', 'Lum', 'Col3a1', 'Col1a1', 'Mmp2', 'Pdgfra')
markers = toupper(markers)

mm = match(markers, ggs)
mm = mm[which(!is.na(mm))]

FeaturePlot(aa, features = rownames(aa)[mm]) 

ggsave(filename = paste0(resDir, '/batch1Data_umap_axloltol_BL_celltypes_FB_markers.pdf'), 
       width = 16, height = 12)

markers =  toupper(c('Adgre1', 'Cd68', 'Itgam', 'Csf1r', "H2-Ab1", 'Mertk',
                             'Cd14',  'Cx3cr1', # pan macrophage
                             'Tlr2', 'Nos2', 'Cd80', 'Cd86', 'Ifng', # M1 (pro-inflamatory)
                             'Arg1', 'Cd163', 'Il4', 'Irf4' # M2 (anti-)
))
mm = match(markers, ggs)
mm = mm[which(!is.na(mm))]
FeaturePlot(aa, features = rownames(aa)[mm]) 

ggsave(filename = paste0(resDir, '/batch1Data_umap_axloltol_BL_celltypes_macrophage_markers.pdf'), 
       width = 16, height = 12)


markers =  toupper(c('Wnt3a', 'Wnt5a'))
mm = match(markers, ggs)
mm = mm[which(!is.na(mm))]
FeaturePlot(aa, features = rownames(aa)[mm])

aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)

aa <- FindNeighbors(aa, dims = 1:20)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)

DimPlot(aa, label = TRUE, repel = TRUE)

ggsave(filename = paste0(resDir, '/batch1Data_umap_axloltol_BL_clusters.pdf'), 
       width = 12, height = 8)


##########################################
# quickly check the macrophage subtype markers 
##########################################
#aa = readRDS(file = paste0(RdataDir, 
#                   '/axoltol_limbBlatema_batch1_DFout_FB_Mphg_subtypes.rds'))
#aa = subset(aa, cells = colnames(aa)[which(aa$DF_out == 'Singlet')])
#ggs = rownames(aa)
#ggs = get_geneName(ggs)
subs = subset(aa, cells = colnames(aa)[which(aa$celltype == "Macrophages")])

subs = NormalizeData(subs, normalization.method = "LogNormalize", scale.factor = 10000)
subs <- FindVariableFeatures(subs, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs

subs <- ScaleData(subs, features = rownames(subs))
subs <- RunPCA(subs, features = VariableFeatures(object = subs), verbose = FALSE, weight.by.var = TRUE)

ElbowPlot(subs, ndims = 50)

subs <- FindNeighbors(subs, dims = 1:20)
subs <- FindClusters(subs, verbose = FALSE, algorithm = 3, resolution = 0.3)

subs <- RunUMAP(subs, reduction = "pca", dims = 1:50, n.neighbors = 30,  min.dist = 0.3)

DimPlot(subs, group.by = 'subtypes', label = TRUE, repel = TRUE, label.size = 6)

p1 = DimPlot(subs, group.by = 'time', label = TRUE, repel = TRUE, label.size = 6)
p2 = DimPlot(subs, group.by = 'seurat_clusters', label = TRUE, repel = TRUE, label.size = 6)

p1 + p2

p3 = DimPlot(subs, group.by = 'subtypes', label = TRUE, repel = TRUE, label.size = 4) + NoLegend()
p1 + p3

#DimPlot(subs, group.by = 'DF_out', label = TRUE, repel = TRUE)
ggsave(filename = paste0(resDir, '/Tobie_batch1Data_axloltolBlastema_macrophage_subclusters.pdf'), 
       width = 14, height = 6)

subs$clusters =  subs$seurat_clusters

saveRDS(subs, file = paste0(RdataDir, 
                          '/axoltol_limbBlatema_batch1_macrophage_time_subtypeAnnotations_rmEpidermis_v3.rds'))


markers =  toupper(c('Adgre1', 'Cd68', 'Itgam', 'Csf1r', "H2-Ab1", 'Mertk',
                     'Cd14',  'Cx3cr1', # pan macrophage
                     'Tlr2', 'Nos2', 'Cd80', 'Cd86', 'Ifng', # M1 (pro-inflamatory)
                     'Arg1', 'Cd163', 'Il4', 'Irf4' # M2 (anti-)
))

markers = c('FN1', 'DPT', 'LUM', 'TWIST1', 'VIM', 'MMP9', 'OTOG', "FREM2", 'WNT5A', 'WNT3A', 'SPARC', 'ARGN')
mm = match(markers, ggs)
mm = mm[which(!is.na(mm))]
FeaturePlot(subs, features = rownames(subs)[mm])

ggsave(filename = paste0(resDir, '/batch1Data_umap_axloltol_BL_celltypes_macrophage_markers.pdf'), 
       width = 16, height = 12)


#subs$clusters = subs$RNA_snn_res.0.5


DimPlot(subs, group.by = 'clusters', label = TRUE, repel = TRUE)
DimPlot(subs, group.by = 'Phase')


# https://www.biocompare.com/Editorial-Articles/566347-A-Guide-to-Macrophage-Markers/
markers = c('ADGRE1', 'CD14', 'CD68', 'CX3CR1', 'ITGAM', 'CSF1R', 'H2AB1', 'MERTK',# pan macrophage
            'TLR2', 'NOS2', 'CD80', 'CD86', 'IFNG', # M1
             'ARG1', 'CD163', 'IL4', 'IRF4', # M2
            'FABP1', 'SIGLEC1', 'MARCO', 'C1QB', 'C1QC' #https://www.sciencedirect.com/science/article/pii/S0014482720303967
            )

mm = unique(c(match(markers, ggs), match(c("CTSK-AMEX60DD014742", "CTSK-AMEX60DD014744"), rownames(subs))))
mm = mm[which(!is.na(mm))]

FeaturePlot(subs, features = rownames(aa)[mm]) 
  #&
  #scale_color_distiller(palette = "RdYlBu")
  #scale_color_viridis_c()
ggsave(filename = paste0(resDir, '/Tobie_umap.harmony_axloltol_limbBlatema_someGeneMarkers.pdf'), 
       width = 12, height = 8)

## cycling cluster
DimPlot(subs, group.by = 'Phase')

## axolotl-specific cluster
genes = c("CTSK-AMEX60DD014744", "TNF-AMEX60DD010588", 'MARCO-AMEX60DD056028')
FeaturePlot(subs, features = genes)
VlnPlot(subs, features = genes, split.by = "clusters")

## resident cluster
genes = c("TREM2-AMEX60DD009193", 'MARCO-AMEX60DD056028', "MPEG1-AMEX60DD007152")
FeaturePlot(subs, features = genes)
VlnPlot(subs, features = genes, split.by = "clusters")

## 
genes = c("MPEG1-AMEX60DD007152", "ARG1-AMEX60DD034655", "CSF1R-AMEX60DD030145",  
          "TREML1-AMEX60DD009191", "CHIT1-AMEX60DD009133",   
          "CD274-AMEX60DD043441", "IRF4-AMEX60DD027194", 
          "FABP1-AMEX60DD046133",  'SIGLEC1-AMEX60DD015921','C1QB-AMEX60DD052070',
          rownames(subs)[grep('IRF5|PPARG-', rownames(subs))])

FeaturePlot(subs, features = genes)
VlnPlot(subs, features = genes, split.by = "clusters")

ggsave(filename = paste0(resDir, '/someMarkerGenes_macrophages_subtypes.pdf'), 
       width = 16, height = 10)


FeaturePlot(subs, features = c("FABP1-AMEX60DD046133",  'SIGLEC1-AMEX60DD015921','MARCO-AMEX60DD056028',
                             'C1QB-AMEX60DD052070', rownames(subs)[grep('MPEG1', ggs)]))

FeaturePlot(subs, features = c('ARG1-AMEX60DD034655', 'TLR2-AMEX60DD015195', 'MARCO-AMEX60DD056028'))

FeaturePlot(subs, features = c('ARG1-AMEX60DD034655', 'MARCO-AMEX60DD056028', "MPEG1-AMEX60DD007152",
                               rownames(subs)[grep('CSF1R|TREML1|CHIT1|TREM2', rownames(subs))]))

subs$clusters[which(subs$clusters == 2)] = 4
subs$clusters = droplevels(subs$clusters)

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
                         'top15_annotatedMouseMarkers_mergingM2M4_v2.pdf'), 
       width = 12, height = 20)


oupMarker %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1.0) %>%
  slice_head(n = 30) %>%
  ungroup() -> top10

write.csv(top10, file = paste0(resDir, '/macrophage_markerGenes4clusters.csv'), row.names = FALSE, quote = FALSE)


#kk = which(subs$clusters != 7)
#subs$subtypes[kk] = paste0('C', subs$clusters[kk], '.ax')

genes = c("FN1-AMEX60DD055420")
FeaturePlot(subs, features = genes)
VlnPlot(subs, features = genes, split.by = "clusters")


subs$subtypes = NA
subs$subtypes[which(subs$clusters == 0)] = 'M0.CXCR4.TLR5.pro-inflammatory'
subs$subtypes[which(subs$clusters == 1)] = 'M1.APOE.PPARG.anti-inflammatory'
subs$subtypes[which(subs$clusters == 2)] = 'M4.ARG1-.CSF1R-.MPEG1-'
subs$subtypes[which(subs$clusters == 3)] = 'M3.CXCR1.DYSF.pro-inflammatory'

subs$subtypes[which(subs$clusters == 4)] = 'M4.ARG1-.CSF1R-.MPEG1-'

subs$subtypes[which(subs$clusters == 5)] = 'M5.FN1+.TREM2-.resident??'
subs$subtypes[which(subs$clusters == 6)] = 'M1.APOE.PPARG.anti-inflammatory'
subs$subtypes[which(subs$clusters == 7)] = 'M7.CTSK+.TREML1+.MARCO.low'


p0 = DimPlot(subs, group.by = 'time', label = TRUE, repel = TRUE, label.size = 5)
p1 = DimPlot(subs, group.by = 'clusters', label = TRUE, repel = TRUE, label.size = 5)
p2 = DimPlot(subs, group.by = 'subtypes', label = TRUE, repel = TRUE)

p0 + p1 + p2

ggsave(filename = paste0(resDir, 
                         '/Tobie_batch1Data_umap_axloltol_BL_celltypes_macrophage_subtypes_time_annotations.pdf'), 
       width = 20, height = 6)


saveRDS(subs, file = paste0(RdataDir, 
                          '/axoltol_limbBlatema_batch1_macrophage_time_subtypeAnnotations_v2.rds'))



subs = readRDS(file = paste0(RdataDir, 
                             '/axoltol_limbBlatema_batch1_macrophage_time_subtypeAnnotations_v2.rds'))


aa$cluster = NA
aa$subtypes = NA

jj = match(colnames(subs), colnames(aa))
aa$cluster[jj] = paste0('M', subs$clusters)
aa$subtypes[jj] = subs$subtypes

rm(subs)

saveRDS(aa, file = paste0(RdataDir, 
                          '/axoltol_limbBlatema_batch1_with_macrophage.subtypes.rds'))

##########################################
# make plots 
##########################################
subs$cluster = NA
## proportions of cluster at each time points
#subs$clusters[which(subs$clusters == 6)] = 1
#subs$clusters = droplevels(subs$clusters)

pct = table(subs$subtypes, subs$time)
for(n in 1:ncol(pct)) pct[,n] = pct[,n]/sum(pct[,n])
pct = data.frame(pct)
df = pct 
colnames(df) = c('cluster', 'condition', 'pct')

ggplot(df, aes(x = condition, y = pct)) +
  geom_bar(stat = 'identity', aes(fill = cluster)) +
  
  ylab("% of cluster ") + 
  xlab("") + 
  #ylim(0, 1.2) + 
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 45, size = 14, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 14)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 12)) +
  #scale_fill_brewer(palette = "Set1")
  scale_fill_manual(values=c("#054674", '#25aff5', "#4d7ea9", '#D4D915','#ff9a36','#B95FBB'))
                             #'#31C53F', "darkgreen", "darkorange", "red", 'magenta', 'gray', 'green', 'blue', 'black')) 

ggsave(filename = paste0(resDir, '/axloltol_Blastema_macrophage_subclusters_proportions.pdf'), 
       width = 10, height = 6)


subs = readRDS(file = paste0(RdataDir, 
                    '/axoltol_limbBlatema_batch1_macrophage_subtypes.rds'))

kk = which(subs$subtypes == 'cycling')
subs$subtypes[kk] = paste0(subs$subtypes[kk], '.ax')

"#12400c", "#2d6624","#1d4f15", "#174711", "#2d6624", "#3d7f33", "#3b7b30", "#468b3b", "#4f9843","#5dae50", "#66bb58", "#72cd64", "#306a26", "#78d669", "#81e472"

#gaba
"#700209", "#75090e","#7a0f13", "#801517", "#851a1b", "#8a1f1f", "#902423", "#952927", "#9a2d2c","#a03230", "#a53634", "#aa3a39", "#b03f3d","#b54342", "#ba4846", "#c04c4b", "#c5504f", "#ca5554", "#d05959", "#d55e5e","#73050c", "#780c11","#8d2221", "#982b2a","#a23432", "#a83837", "#b2413f", "#b84544", "#bd4a49", "#c85352", #"#cd5756",
#glut
"#054674", "#134d7b","#1d5481", "#265a88", "#2e618e", "#73a4cb", "#366995", "#3e709c", "#4677a2","#4d7ea9", "#5586b0", "#5c8db7", "#6495bd","#6b9cc4", "#7bacd2", "#8ebfe4", "#96c7eb", "#9ecff2", "#18507e", "#18507e","#2a5e8b", "#497ba6","#5889b3", "#6fa0c8","#7fafd6", "#6091ba", "#5182ac", "#3a6c98", "#a6d7f9",
#npc
"#ffb120", "#feb72a","#fdbc34", "#fcc13d", "#fbc745", "#facc4e", "#f9d156", "#f8d65f", "#f8da68","#f7df70", "#f7e479", "#f7e882", "#f7ed8a", "#f7f193", "#eca319"

my_cols <- c('3'='#F68282','15'='#31C53F','5'='#1FA195','1'='#B95FBB','13'='#D4D915',
             '14'='#28CECA','9'='#ff9a36','8'='#2FF18B','11'='#aeadb3','6'='#faf4cf',
             '2'='#CCB1F1','12'='#25aff5','7'='#A4DFF2','4'='#4B4BF7','16'='#AC8F14',
             '10'='#E6C122')
'11'='#aeadb3','6'='#faf4cf',
'2'='#CCB1F1','12'='#25aff5')

cols = c('#F68282','#31C53F','#1FA195','#B95FBB','#D4D915',
         '#28CECA', '#ff9a36', '#2FF18B')

DimPlot(subs, group.by = 'subtypes', label = FALSE, repel = TRUE) +
  #labs(x = 'mouse', y = 'axolotl', fill = "Spearman's\nrho", size = "-log10\nadj. p-value")+
  theme_classic()+
  #scale_discrete_manual(palette = "Set1") +
  theme(axis.title = element_text(colour = "black", face = "bold"),
        axis.text = element_text(colour = "black"),
        #theme(axis.text.x = element_text(angle = 90, size = 10)),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))

ggsave(filename = paste0(resDir, '/Tobie_batch1Data_axloltolBlastema_macrophage_subclustersV2.pdf'), 
       width = 8, height = 6)


##########################################
# check the Maf genes 
##########################################
genes = c('MAF', 'MAFA', 'MAFB', 'NRL')
mm = grep('^MAF|NRL|KAZALD1', rownames(subs))

FeaturePlot(subs, features = c("MAFB-AMEX60DD028613", "KAZALD1-AMEX60DD046379")) &
  scale_color_gradient(low = "grey", high = "brown")

ggsave(filename = paste0(resDir, '/axoltol_blastema_macrophage_MAFB_KAZALD1.pdf'), 
       width = 16, height = 6)


##########################################
# check the transcripts expressed in embyonic macrophage in mouse 
## refs: https://www.sciencedirect.com/science/article/pii/S0012160622000859?via%3Dihub#bib57
## https://www.science.org/doi/full/10.1126/science.aaf4238
##########################################
subs = readRDS(file = paste0(RdataDir, 
              '/axoltol_limbBlatema_batch1_macrophage_time_subtypeAnnotations_rmEpidermis_v3.rds'))

## Timd4, Csf1r, Maf, Fcgr1 (CD64), Emr (F4.80), and Mrc1 (CD206), Timd4, Lyve1, Folr2 and/or Ccr2
genes = rownames(subs)[grep('TIMD4|CSF1R|LYVE1', rownames(subs))]

FeaturePlot(subs, features = genes) &
  scale_color_gradient(low = "grey", high = "brown")

ggsave(filename = paste0(resDir, '/axoltol_blastema_macrophage_embryonicMarkers.pdf'), 
       width = 16, height = 12)


##########################################
# compare the axolotl macrophages with mouse in pro-fibrotic activation  
##########################################
subs = readRDS(file = paste0(RdataDir, 
              '/axoltol_limbBlatema_batch1_macrophage_time_subtypeAnnotations_rmEpidermis_v3.rds'))


genes = rownames(subs)[grep('MERTK|CD36|MRC1|CD163', rownames(subs))]
FeaturePlot(subs, features = genes)

genes = rownames(subs)[grep('TREM2|FABP5|GPNMB|CD63|SPP1', rownames(subs))]
FeaturePlot(subs, features = genes)

genes = rownames(subs)[grep('PDGFC|IGF1|RETNLA', rownames(subs))]
FeaturePlot(subs, features = genes)


########################################################
########################################################
# Section : annotate the FB subtypes 
# 
########################################################
########################################################
#aa = readRDS(file = paste0(RdataDir, 
#                   '/axoltol_limbBlatema_batch1_with_macrophage.subtypes.rds'))
# aa = readRDS(file = paste0(RdataDir, 
#                            '/axoltol_limbBlatema_batch1_DFout_filterCelltypes_geneNames_',
#                            'FB_Mphg.cleaned_subtypes_v3.rds'))
aa = readRDS(file = paste0(RdataDir,
                           '/axoltol_limbBlatema_batch1_DFout_filterCelltypes_geneNames_',
                           'FB.clean_Mphg.cleaned_subtypes_v3.rds'))
ggs = rownames(aa)
ggs = get_geneName(ggs)

subs = subset(aa, cells = colnames(aa)[which(aa$celltype == "FB")])

subs = NormalizeData(subs, normalization.method = "LogNormalize", scale.factor = 10000)
subs <- FindVariableFeatures(subs, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs

subs <- ScaleData(subs, features = rownames(subs))
subs <- RunPCA(subs, features = VariableFeatures(object = subs), verbose = FALSE, weight.by.var = TRUE)

ElbowPlot(subs, ndims = 50)

subs <- FindNeighbors(subs, dims = 1:20)
subs <- FindClusters(subs, verbose = FALSE, algorithm = 3, resolution = 0.4)

subs <- RunUMAP(subs, reduction = "pca", dims = 1:50, n.neighbors = 100,  min.dist = 0.3)

p1 = DimPlot(subs, group.by = 'time', label = TRUE, repel = TRUE)
p2 = DimPlot(subs, group.by = 'seurat_clusters', label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(resDir, '/Tobie_batch1Data_axloltolBlastema_FB_subclusters.pdf'), 
       width = 14, height = 6)

p3 = DimPlot(subs, group.by = 'subtypes', label = TRUE, repel = TRUE) + NoLegend()

p1 + p3

subs$clusters = subs$RNA_snn_res.0.4
subs$clusters[which(subs$clusters == 5|subs$clusters == 6)] = 4
subs$clusters = droplevels(subs$clusters)

DimPlot(subs, group.by = 'Phase', label = TRUE, repel = TRUE)

p1 = DimPlot(subs, group.by = 'time', label = TRUE, repel = TRUE)
p2 = DimPlot(subs, group.by = 'clusters', label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(resDir, '/Tobie_batch1Data_axloltolBlastema_FB_subclusters_mergingEarlyClusters.pdf'), 
       width = 14, height = 6)

DimPlot(subs, group.by = 'Phase')

markers =  toupper(c('Procr', 'Dpt', 'Pi16', 'Col1a2', 'Acta2', 'Lum', 'Col3a1', 'Col1a1', 'Mmp2', 
                                'Pdgfra',  'Twist1', 'Vim', 'Mmp9'))
mm = match(markers, ggs)
mm = mm[which(!is.na(mm))]
FeaturePlot(subs, features = rownames(subs)[mm]) 

ggsave(filename = paste0(resDir, '/batch1Data_umap_axloltol_BL_celltypes_FB_markers.pdf'), 
       width = 16, height = 12)


### calculate marker genes for defined clusters 
#DimPlot(subs, group.by = 'clusters', label = TRUE, repel = TRUE)
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
                         '/Tobie_batch1Data_umap_axloltol_BL_celltypes_FB_subclustering_',
                         'top15_annotatedMouseMarkers_mergingEarlyClusters.pdf'), 
       width = 12, height = 20)


oupMarker %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1.0) %>%
  slice_head(n = 30) %>%
  ungroup() -> top10

write.csv(top10, file = paste0(resDir, '/FBs_markerGenes_top30_for_definedClusters.csv'), row.names = FALSE, 
          quote = FALSE)


## annotate the FB subtypes by manually checking marker genes
genes = c(rownames(subs)[grep('POSTN|ACTA2|FBN1|TNXB|MATN2|LRRC15', rownames(subs))],
          "RUNX2-AMEX60DD032922")

FeaturePlot(subs, features = genes)
VlnPlot(subs, features = genes, split.by = "clusters")

genes = c(rownames(subs)[grep('SFRP5|SFRP2|FSTL1|CTHRC1', rownames(subs))])
FeaturePlot(subs, features = genes)
VlnPlot(subs, features = genes, split.by = "clusters")


genes = c(rownames(subs)[grep('MMP8|CXCL2|IL11-|PDPN|CCL2-|CCL8-|TNC', rownames(subs))])
FeaturePlot(subs, features = genes)
VlnPlot(subs, features = genes, split.by = "clusters")


genes = c(rownames(subs)[grep('SFRP2|ALDH1A3|FSTL1|COL1A1|MGP|NREP|DCN|LUM|ACTA2|LRRC15|POSTN', 
                              rownames(subs))])
genes = c("LUM-AMEX60DD007586", "DCN-AMEX60DD007589", 
          "ACTA2-AMEX60DD052517","POSTN-AMEX60DD049175", "LRRC15-AMEX60DD001730",     
    "COL1A1-AMEX60DD009937", "SFRP2-AMEX60DD045083", "SFRP2-AMEX60DD049867",  "FSTL1-AMEX60DD046943", "MGP-AMEX60DD030520",
    "RUNX2-AMEX60DD032922")

FeaturePlot(subs, features = genes)
VlnPlot(subs, features = genes, split.by = "cluster")

genes = c(rownames(subs)[grep('MYOG|PAX7|ACTN3', rownames(subs))])

FeaturePlot(subs, features = genes)
VlnPlot(subs, features = genes, split.by = "clusters")


genes = c(rownames(subs)[grep('^FSP1|VIM', rownames(subs))])
FeaturePlot(subs, features = genes)

FeaturePlot(subs, features = genes)
VlnPlot(subs, features = genes, split.by = "clusters")


FeaturePlot(subs, features = rownames(subs)[grep('DPT|PI16', rownames(subs))])

FeaturePlot(subs, features = rownames(subs)[grep('POSTN|ACTA2|LRRC15|RUNX2', rownames(subs))])

FeaturePlot(subs, features = rownames(subs)[grep('SFRP2|CTHRC1|FASTL1', rownames(subs))])

FeaturePlot(subs, features = rownames(subs)[grep('PDPN|VCAM1|CSF2', rownames(subs))])

FeaturePlot(subs, features = rownames(subs)[grep('TNC|SERPINE1|IL11', rownames(subs))])

FeaturePlot(subs, features = rownames(subs)[grep('CD68|MARCO', rownames(subs))])

FeaturePlot(subs, features = rownames(subs)[grep('IL11-|MMP1-|CXCL8|IL7R', rownames(subs))])

#kk = which(subs$clusters != 7)
#subs$subtypes[kk] = paste0('C', subs$clusters[kk], '.ax')

subs$subtypes = NA
subs$subtypes[which(subs$clusters == 0)] = 'FB0.POSTN+.ACTA2+.LRRC15+.SFRP2+'
subs$subtypes[which(subs$clusters == 1)] = 'FB1.RUNX2+.LRRC15+.SFRP2+'
subs$subtypes[which(subs$clusters == 2)] = 'FB2.proliferation'

subs$subtypes[which(subs$clusters == 3)] = 'FB3.early.POSTN+.ACTA2+.LRRC15+.SFRP2+.TIG1+.MyoFB'

subs$subtypes[which(subs$clusters == 4)] = 'FB4.early.TNC+.IL11+.proInflammatory'
subs$subtypes[which(subs$clusters == 7)] = 'FB7.PITX2+.SIX1+'


p0 = DimPlot(subs, group.by = 'time', label = TRUE, repel = TRUE)
p1 = DimPlot(subs, group.by = 'clusters', label = TRUE, repel = TRUE)
p2 = DimPlot(subs, group.by = 'subtypes', label = TRUE, repel = TRUE) + NoLegend()

p0 + p1 + p2

ggsave(filename = paste0(resDir, 
                         '/Tobie_batch1Data_umap_axloltol_BL_celltypes_FB_subtypes_time_annotations.pdf'), 
       width = 24, height = 6)



#subs$clusters = as.factor(gsub('FB', '', subs$cluster))
subs$subtypes[which(subs$cluster == "FB0")] = 'FB0.POSTN+.ACTA2+.LRRC15+.SFRP2+.MGP+'
subs$subtypes[which(subs$cluster == 'FB1')] = 'FB1.RUNX2+.FOXC1+.GDF6+.TIMP1+'
subs$subtypes[which(subs$cluster == 'FB2')] = 'FB2.STMN1+'

subs$subtypes[which(subs$cluster == 'FB3')] = 'FB3.early.MyoFB.like.POSTN+.ACTA2+.LRRC15+.SFRP2+.MGP+.TIG1+'
subs$subtypes[which(subs$cluster == "FB4")] = 'FB4.early.InflammatoryResp.TNC+.IL11+'


saveRDS(subs, file = paste0(RdataDir, 
                            '/axoltol_limbBlatema_batch1_FB_time_subtypeAnnotations_v3.rds'))

subs = readRDS(file = paste0(RdataDir, 
                             '/axoltol_limbBlatema_batch1_FB_time_subtypeAnnotations_v3.rds'))

jj = match(colnames(subs), colnames(aa))
#aa$cluster[jj] = paste0('FB', subs$clusters)
aa$subtypes[jj] = subs$subtypes

#jj = which(is.na(aa$subtypes))
#aa$subtypes[jj] = aa$celltype[jj]
#aa$cluster[jj] = aa$celltype[jj]


DimPlot(aa, group.by = 'cluster', label = TRUE, repel = TRUE, label.size = 5) + NoLegend() +
  scale_color_brewer(palette = "Set1")

ggsave(filename = paste0(resDir, '/Tobie_batch1Data_axloltolBlastema_FB_M_others_v2.pdf'), 
       width = 8, height = 6)

# aa$time = droplevels(aa$time)
# DimPlot(aa, group.by = 'cluster', label = TRUE, repel = TRUE, label.size = 5, split.by = 'time') + 
#   NoLegend()
# 
# ggsave(filename = paste0(resDir, '/Tobie_batch1Data_axloltolBlastema_FB_M_others_timePoints.pdf'), 
#        width = 24, height = 6)

saveRDS(aa, file = paste0(RdataDir, 
              '/axoltol_limbBlatema_batch1_DFout_filterCelltypes_FBcleaned._Mphg.cleaned_subtypes_v3.rds'))

##########################################
# make plot of temporal changes
##########################################
subs = readRDS(file = paste0(RdataDir, 
                             '/axoltol_limbBlatema_batch1_FB_time_subtypeAnnotations_v3.rds'))

pct = table(subs$cluster, subs$time)
for(n in 1:ncol(pct)) pct[,n] = pct[,n]/sum(pct[,n])
pct = data.frame(pct)
df = pct 
colnames(df) = c('cluster', 'condition', 'pct')

ggplot(df, aes(x = condition, y = pct)) +
  geom_bar(stat = 'identity', aes(fill = cluster)) +
  
  ylab("% of cluster ") + 
  xlab("") + 
  #ylim(0, 1.2) + 
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 45, size = 14, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 14)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 10)) +
  scale_fill_brewer(palette = "Set1")

  #scale_fill_manual(values=c("#054674", '#25aff5', "#4d7ea9", '#D4D915','#ff9a36','#B95FBB'))
#'#31C53F', "darkgreen", "darkorange", "red", 'magenta', 'gray', 'green', 'blue', 'black')) 

ggsave(filename = paste0(resDir, '/axloltol_Blastema_FB_subclusters_proportions.pdf'), 
       width = 12, height = 6)



genes = c("LUM-AMEX60DD007586", "DCN-AMEX60DD007589", 
          "ACTA2-AMEX60DD052517","POSTN-AMEX60DD049175", "LRRC15-AMEX60DD001730",     
     "SFRP2-AMEX60DD049867",  "FSTL1-AMEX60DD046943", "MGP-AMEX60DD030520",
    "RUNX2-AMEX60DD032922",     
    "COL1A1-AMEX60DD009937", "STMN3-AMEX60DD027088", "TIMP1-AMEX60DD020992", "FOXC1-AMEX60DD038386", 
    "STMN1-AMEX60DD005746", "IL11-AMEX60DD016524", "GDF6-AMEX60DD033076",   "TNC-AMEX60DD050822")

subs$cluster = factor(subs$cluster, levels = c('FB4', 'FB3', 'FB2', 'FB1', 'FB0'))

DotPlot(subs, features = genes, group.by = 'cluster') + RotatedAxis() +
  scale_color_gradientn(colors = c("#F0F0F0", "#EFFAB6", "#69C6BE", "#007BB7", "#121D60"))  +
  labs( x = '', y = '' )

ggsave(filename = paste0(resDir, '/axloltol_Blastema_FB_subclusters_markers.pdf'), 
       width = 12, height = 6)


########################################################
########################################################
# Section : # make plots for grant figures 
# 
########################################################
########################################################

##########################################

##########################################
ax = readRDS(file = paste0(RdataDir, 
                           '/axoltol_limbBlatema_batch1_macrophage_time_subtypeAnnotations_rmEpidermis_v3.rds'))

#ax[["RNA"]] <- as(object = ax[["RNA5"]], Class = "Assay")

ax$species = 'ax'
ax$time = droplevels(ax$time)

p1 = DimPlot(ax, group.by = 'cluster', label = TRUE, repel =  TRUE)
p2 = DimPlot(ax, group.by = 'subtypes', label = TRUE, repel =  TRUE) + NoLegend()
p1 + p2

ax$cluster[which(ax$subtypes == "M0.CXCR4.TLR5.pro-inflammatory")] = 'M1'

ax$cluster[which(ax$subtypes == "M1.APOE.PPARG.anti-inflammatory" & ax$cluster == 'M1')] = 'M2' 
ax$cluster[which(ax$subtypes == "M3.CXCR1.DYSF.pro-inflammatory")] = 'M3' 
ax$cluster[which(ax$subtypes == "M7.CTSK+.TREML1+.MARCO.low")] = 'M4'
ax$cluster[which(ax$subtypes == "M4.ARG1-.CSF1R-.MPEG1-")] = 'M5' 
ax$cluster[which(ax$subtypes == "M1.APOE.PPARG.anti-inflammatory" & ax$cluster == 'M6')] = 'M6.cyling' 

#ax$subtypes[which(ax$subtypes == 'cycling')] = 'cycling.ax'
jj = which(ax$cluster == 'M6.cyling')
ax$cluster[jj] = 'M2'
ax$subtypes[jj] = 'M1.APOE.PPARG.anti-inflammatory'

DimPlot(ax, group.by = 'cluster', label = TRUE, repel = TRUE)

cols_macrophage = c( "#265a88", "#8BBEDC",  '#A4DFF2', '#25aff5', "#007BB7")

DimPlot(ax, group.by = 'cluster', label = TRUE, repel =  TRUE) + 
  #theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 14, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 14)) +
  scale_color_manual(values=cols_macrophage) + NoLegend() +
  ggtitle('')

ggsave(filename = paste0(resDir, '/axolotl_macrophage_subtypes_v3.pdf'), 
       width = 5, height = 4)


aa = readRDS(file = paste0(RdataDir, 
              '/axoltol_limbBlatema_batch1_DFout_filterCelltypes_FBcleaned._Mphg.cleaned_subtypes_v3.rds'))


jj = match(colnames(ax), colnames(aa))
aa$cluster[jj] = ax$cluster

"#12400c", "#2d6624","#1d4f15", "#174711", "#2d6624", "#3d7f33", "#3b7b30", "#468b3b", "#4f9843","#5dae50", "#66bb58", "#72cd64", "#306a26", "#78d669", "#81e472"

#gaba
"#700209", "#75090e","#7a0f13", "#801517", "#851a1b", "#8a1f1f", "#902423", "#952927", "#9a2d2c","#a03230", "#a53634", "#aa3a39", "#b03f3d","#b54342", "#ba4846", "#c04c4b", "#c5504f", "#ca5554", "#d05959", "#d55e5e","#73050c", "#780c11","#8d2221", "#982b2a","#a23432", "#a83837", "#b2413f", "#b84544", "#bd4a49", "#c85352", #"#cd5756",
#glut
"#054674", "#134d7b","#1d5481", "#265a88", "#2e618e", "#73a4cb", "#366995", "#3e709c", "#4677a2","#4d7ea9", "#5586b0", "#5c8db7", "#6495bd","#6b9cc4", "#7bacd2", "#8ebfe4", "#96c7eb", "#9ecff2", "#18507e", "#18507e","#2a5e8b", "#497ba6","#5889b3", "#6fa0c8","#7fafd6", "#6091ba", "#5182ac", "#3a6c98", "#a6d7f9",
#npc
"#ffb120", "#feb72a","#fdbc34", "#fcc13d", "#fbc745", "#facc4e", "#f9d156", "#f8d65f", "#f8da68","#f7df70", "#f7e479", "#f7e882", "#f7ed8a", "#f7f193", "#eca319"

my_cols <- c('3'='#F68282','15'='#31C53F','5'='#1FA195','1'='#B95FBB','13'='#D4D915',
             '14'='#28CECA','9'='#ff9a36','8'='#2FF18B','11'='#aeadb3','6'='#faf4cf',
             '2'='#CCB1F1','12'='#25aff5','7'='#A4DFF2','4'='#4B4BF7','16'='#AC8F14',
             '10'='#E6C122')
'11'='#aeadb3','6'='#faf4cf',
'2'='#CCB1F1','12'='#25aff5')

"#EDF8E9" "#C4E8BD" "#90D08E" "#59B668" "#27984B" "#006D2C"
"#EDF8E9" "#BAE4B3" "#74C476" "#31A354" "#006D2C"


cols = c('#F68282','#31C53F','#1FA195','#B95FBB','#D4D915',
         '#28CECA', '#ff9a36', '#2FF18B')

cols_cluster = c( "#7F7F7F", "#FFC000", "#70AD47", "#337f01", "#265401", "#CD00CF", "#800080", 
                  "#EFFAB6", "#69C6BE", "#007BB7", "#121D60",
                  "#69C6BE", "#007BB7", "#121D60")




cols_macrophage = c( "#265a88", "#8BBEDC",  '#A4DFF2', '#25aff5', "#007BB7")
cols_FB = c("#EFFAB6","#BAE4B3", "#69C6BE", '#2FF18B', '#31C53F')

cols_cluster = c("#12400c", "#952927", "#800080",'#F68282', cols_FB, cols_macrophage, 
                 "#FFC000", '#CCB1F1', "#CD00CF") 

DimPlot(aa, group.by = 'cluster', label = TRUE, repel = TRUE) + NoLegend() + 
  theme(axis.text.x = element_text(angle = 0, size = 14, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 14)) + 
  ggtitle('') + 
  scale_color_manual(values=cols_cluster)

ggsave(filename = paste0(resDir, '/axolotl_all_celltypes_v4.pdf'), 
       width = 6, height = 5)
  


### FB temporal changes
subs = readRDS(file = paste0(RdataDir, 
                             '/axoltol_limbBlatema_batch1_FB_time_subtypeAnnotations_v3.rds'))

pct = table(subs$cluster, subs$time)
for(n in 1:ncol(pct)) pct[,n] = pct[,n]/sum(pct[,n])
pct = data.frame(pct)
df = pct 
colnames(df) = c('cluster', 'condition', 'pct')

ggplot(df, aes(x = condition, y = pct)) +
  geom_bar(stat = 'identity', aes(fill = cluster)) +
  
  ylab("% of cluster ") + 
  xlab("") + 
  #ylim(0, 1.2) + 
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 45, size = 14, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 14)) +
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.8,.9), legend.direction = "vertical") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 10)) +
  scale_fill_manual(values = cols_FB)

#scale_fill_manual(values=c("#054674", '#25aff5', "#4d7ea9", '#D4D915','#ff9a36','#B95FBB'))
#'#31C53F', "darkgreen", "darkorange", "red", 'magenta', 'gray', 'green', 'blue', 'black')) 
ggsave(filename = paste0(resDir, '/axloltol_FB_subclusters_proportions.pdf'), 
       width = 6, height = 4)


### plots new marker genes in FB
subs = readRDS(file = paste0(RdataDir, 
                             '/axoltol_limbBlatema_batch1_FB_time_subtypeAnnotations_v3.rds'))

p1 = DimPlot(subs, group.by = 'cluster', label = TRUE, repel = TRUE, label.size = 5) +
  NoLegend()

##  Crabp1, Twist2, Lrg1
genes = c(rownames(subs)[grep('CRABP1|^LRG1|ACTA2', rownames(subs))], "TWIST2-AMEX60DD029436")
p2 = FeaturePlot(subs, features = genes, ncol = 2, pt.size = 0.8) &
  scale_color_gradient(low = "grey", high = "brown") 

p1 + p2

ggsave(filename = paste0(resDir, '/axoltol_FB_newMarkers_genes.pdf'), 
       width = 12, height = 6)

### check the pro-fibrotic macrophage activation
genes = c(rownames(ax)[grep('CLEC10A|NINJ1|TREM2|FABP5|GPNMB|IGF1-|MS4A7|GAS6', 
                            rownames(ax))], "SPP1-AMEX60DD043905")

FeaturePlot(ax, features = genes, ncol = 3, pt.size = 0.8) &
  scale_color_gradient(low = "grey", high = "brown") 

ggsave(filename = paste0(resDir, '/axoltol_macrophage_proFibrotic_genes.pdf'), 
       width = 10, height = 6)


genes = c(rownames(ax)[grep('CLEC10A|NINJ1|TREM2|FABP5|GPNMB|IGF1-|MS4A7|GAS6', 
                            rownames(ax))], "SPP1-AMEX60DD043905")

FeaturePlot(ax, features = genes, ncol = 3, pt.size = 0.8) &
  scale_color_gradient(low = "grey", high = "brown") 

ggsave(filename = paste0(resDir, '/axoltol_macrophage_proFibrotic_genes.pdf'), 
       width = 10, height = 6)



FeaturePlot(ax, features = c("MAFB-AMEX60DD028613", "LYVE1-AMEX60DD004804")) &
  #scale_color_viridis_c()
  scale_color_gradient(low = "grey", high = "brown") 

ggsave(filename = paste0(resDir, '/axoltol_macrophage_MAF.pdf'), 
       width = 10, height = 4)


##########################################
# plot for cross-species with scVI 
##########################################
library(Seurat)
library(SeuratData)
library(SeuratDisk)
mm = readRDS(file = paste0(RdataDir, 'mm_scRNAseq_for_crossSpecies_v3.rds'))
ax = readRDS(file = paste0(RdataDir, 'ax_scRNAseq_for_crossSpecies_v3.rds'))

aa = merge(mm, y = ax)

saveRDS(aa, file = paste0(RdataDir, 'SeuratObj_mm_ax_scRNAseq_merged_v2_seuratV4.rds'))

# Convert(paste0("../results/scRNAseq_analysis_immune_crossSpecies_20250903/",
#                "mm_ax_scRNAseq_merged_v2_seuratV4_analysis_output.h5ad"), dest = "h5seurat", 
#         overwrite = TRUE)
# pbmc3k <- LoadH5Seurat("pbmc3k_final.h5seurat")
# pbmc3k

embedding = read.csv(file = paste0('../results/scRNAseq_analysis_immune_crossSpecies_20250903/', 
                                   'mm_ax_scRNAseq_merged_v2_seuratV4_analysis_output_umap_coordinates.csv'),
                     row.names = c(1))

embedding = as.matrix(embedding)
mm = match(colnames(aa), rownames(embedding))


aa[['umap']] = Seurat::CreateDimReducObject(embeddings=embedding, assay = 'RNA', key = 'UMAP_')

kk = which(aa$cluster == 'M4')
aa$species[kk] = 'ax_M4'


"#facc4e", "#f9d156", "#f8d65f", "#f8da68","#f7df70", "#f7e479", "#f7e882", "#f7ed8a", "#f7f193", "#eca319"

cols_macrophage = c("#265a88", "#8BBEDC",  '#A4DFF2', '#25aff5', "#007BB7")

DimPlot(aa, group.by = 'species') + 
  scale_color_manual(values=c("#265a88", '#25aff5',  '#ff9a36'))

ggsave(filename = paste0(resDir, '/axolotl_mm_crossSpecies_scVI_v2.pdf'), 
       width = 6, height = 4)




##########################################
# FB annotation comparison between axolotl and mouse 
##########################################
library(ggalluvial)




########################################################
########################################################
# Section IV: cell-cell communication analysis using LIANA
# 
########################################################
########################################################
library(pryr) # monitor the memory usage
require(ggplot2)
#library(nichenetr)
library(Seurat) # please update to Seurat V4
library(tidyverse)
library(circlize)
library(RColorBrewer)
require(scran)
require(scater)
#library(nichenetr)
library(tidyverse)
library(circlize)
library(randomcoloR)

source(paste0(functionDir, '/functions_scRNAseq.R'))
source(paste0(functionDir, '/functions_Visium.R'))


refs = readRDS(file = paste0(RdataDir, 
                            '/axoltol_limbBlatema_batch1_DFout_filterCelltypes_geneNames_',
                           'FB_Mphg_subtypes_v2.rds'))

subref = subset(refs, cells = colnames(refs)[which(refs$celltype == 'FB'| refs$celltype == 'Macrophages')])


# run LIANA day by day
timepoint_specific = TRUE

times_slice = unique(subref$time)
#times_slice = c('d7')

subtypes = unique(subref$subtypes)
species = 'ax6'

# set parameter for ligand-receptor analysis
outDir_version = paste0(resDir, 'Ligand_Receptor_analysis_liana')
if(!dir.exists(outDir_version)) dir.create(outDir_version)


##########################################
# run LIANA for all pairs
##########################################
sce <- as.SingleCellExperiment(subref)
colLabels(sce) = as.factor(sce$celltype)

if(species == 'ax6'){
  rownames(sce) = toupper(get_geneName(rownames(sce)))
}else{
  rownames(sce) = toupper(rownames(sce))
}

ave.counts <- calculateAverage(sce, assay.type = "counts")

num.cells <- nexprs(sce, byrow=TRUE)
smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells",
              xlab=expression(Log[10]~"average count"))

# detected in >= 5 cells, ave.counts >=5 but not too high
genes.to.keep <- num.cells > 20 & ave.counts >= 10^-4  & ave.counts <10^3  
summary(genes.to.keep)

sce <- sce[genes.to.keep, ]

## run the liana wrap function by specifying resource and methods
# Resource currently included in OmniPathR (and hence `liana`) include:
show_resources()
# Resource currently included in OmniPathR (and hence `liana`) include:
show_methods()

if(min(exec("logcounts", sce)) < 0){
  xx = logcounts(sce)
  xx[which(xx<0)] = 0
  xx = Matrix(xx, sparse = TRUE)
  logcounts(sce) = xx
  rm(xx)
}

liana_test <- liana_wrap(sce,  
                         # method = c("natmi", "connectome", "logfc", "sca", "cytotalk"),
                         method = c("natmi", "connectome", "logfc", "sca"),
                         #resource = c("Consensus", 'CellPhoneDB', "OmniPath", "LRdb", "CellChatDB",  
                         # "CellTalkDB"), 
                         resource = c("Consensus"),
                         assay.type = "logcounts", 
                         idents_col = 'subtypes')

# Liana returns a list of results, each element of which corresponds to a method
# liana_test %>% glimpse

# We can aggregate these results into a tibble with consensus ranks
liana_test <- liana_test %>%
  liana_aggregate(resource = 'Consensus')

saveRDS(liana_test, file = paste0(outDir_version, '/res_lianaTest_Consensus_subtypes_allpairs.rds'))


##########################################
# make ligand-receptor plots 
##########################################
## test celltalker
library(celltalker)
library(SeuratData)
library(Connectome)
library(cowplot)

liana_test = readRDS(file = paste0(outDir_version, '/res_lianaTest_Consensus_subtypes_allpairs.rds'))

pct_cutoff = 0.05

for(n in 1:length(times_slice))
{
  # n = 3
  
  time = times_slice[n]
  cat(' run LIANA for time -- ', as.character(time), '\n')
  
  outDir = paste(outDir_version, '/', time, collapse = '')
  outDir = gsub(' ', '', outDir)
  if(!dir.exists(outDir)) dir.create(outDir)
  
  #for(n in 1:ncol(pct)) pct[,n] = pct[,n]/sum(pct[,n])
  
  pct_FB = table(subref$subtypes[which(subref$time == time & subref$celltype == 'FB')])
  pct_FB = pct_FB/sum(pct_FB)
  subtypes_FB = names(pct_FB)[which(pct_FB >= pct_cutoff)]
  
  pct_M = table(subref$subtypes[which(subref$time == time & subref$celltype == 'Macrophages')])
  pct_M = pct_M/sum(pct_M)
  subtypes_M = names(pct_M)[which(pct_M >= pct_cutoff)]
  
  source(paste0(functionDir, "/functions_cccInference.R"))
  
  #res = aggregate_output_LIANA(liana_out = paste(outDir))
  res = data.frame(liana_test)
  
  ii = which(!is.na(match(res$source, subtypes_FB)) & !is.na(match(res$target, subtypes_M)))
  jj = which(!is.na(match(res$source, subtypes_M)) & !is.na(match(res$target, subtypes_FB)))
  
  res = res[unique(c(ii, jj)), ]
  
  
  #colnames(res)[1:2] = c('source', 'target')
  colnames(res)[3:4] = c('ligand', 'receptor')
  
  res$weight_norm = res$sca.LRscore
  res$pair = paste0(res$ligand, ' - ', res$receptor)
  res$vector = paste0(res$source, ' - ', res$target)
  res$edge = paste0(res$source, ' - ', res$ligand, ' - ', res$receptor, ' - ', res$target)
  res$source.ligand = paste0(res$source, ' - ', res$ligand)
  res$receptor.target = paste0(res$receptor, ' - ', res$target)
  
  source(paste0(functionDir, '/functions_cccInference.R'))
  
  pdfname = paste0(outDir, '/LR_interactions_allPairs_LIANA_tops.pdf')
  pdf(pdfname, width=12, height = 8)
  
  for(ntop in c(100, 200, 300))
  {
    # ntop = 100
    test = res[c(1:ntop), ]
    test = test[which(test$ligand != 'ACTR2'), ] ## for unknow reason this ligand making problem for plots
    
    #colnames(test)[1:4] = c('sender', 'receiver', 'ligand', 'receptor')
    
    cells.of.interest = unique(c(test$source, test$target))
    cell_color = randomcoloR::distinctColorPalette(length(cells.of.interest))
    names(cell_color) <- cells.of.interest
    
    my_CircosPlot(test, 
                  weight.attribute = 'weight_norm',
                  cols.use = cell_color,
                  sources.include = cells.of.interest,
                  targets.include = cells.of.interest,
                  lab.cex = 0.5,
                  title = paste('LR scores top :', ntop))
    
  }
  
  dev.off()
  
  
}
