##########################################################################
##########################################################################
# Project: Collaboration with Sabine and Sebastian from Cologne
# Script purpose: analyze macrophage subtypes in mouse skin scRNA-seq data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Aug 29 15:45:18 2025
##########################################################################
##########################################################################
rm(list = ls())

library(Seurat)
library(SeuratObject)
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
library(DropletUtils)

species = 'mouse'
version.analysis = '_mouse_20250829'
resDir = paste0("../results/scRNAseq_analysis_immune", version.analysis, '/')
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts'
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_scRNAseq.R')
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_Visium.R')

levels = c('ctl_dpi0', 'ctl_dpi4',  'ctl_dpi7', 'ctl_dpi14')


########################################################
########################################################
# Section I : import data and quick check QCs
# 
########################################################
########################################################
dataDir = '../cellranger_out/'
design = cbind(c(310697, 228615, 228613, 228611),
               c('ctl_dpi0', 'ctl_dpi4', 'ctl_dpi7', 'ctl_dpi14'))
design = data.frame(design)
colnames(design) = c('sampleID', 'condition')

#design = design[-1, ]

# design = design[which(design$sampleID != '196323'), ]

# import data from cellranger output
for(n in 1:nrow(design))
{
  # n = 1
  cat(n, ' : ', design$condition[n], '\n')
  
  topdir = paste0(dataDir, design$sampleID[n], '/outs/raw_feature_bc_matrix/')
  exp = Matrix::readMM(paste0(topdir, "matrix.mtx.gz")) #read matrix
  bc = read.csv(paste0(topdir, "/barcodes.tsv.gz"), header = F, stringsAsFactors = F)
  g = read.csv(paste0(topdir, "/features.tsv.gz"), header = F, stringsAsFactors = F, sep = '\t')
  
  ## make unique gene names
  g$name = g$V2
  gg.counts = table(g$V2)
  gg.dup = names(gg.counts)[which(gg.counts>1)]
  index.dup = which(!is.na(match(g$V2, gg.dup)))
  g$name[index.dup] = paste0(g$V2[index.dup], '_', g$V1[index.dup])
  
  colnames(exp) = bc$V1
  rownames(exp) = g$name
  
  count.data = exp
  rm(exp);
  
  # get emptyDrops and default cutoff cell estimates
  iscell_dd = defaultDrops(count.data, expected = 12000) # default cell estimate, similar to 10x cellranger
  cat(sum(iscell_dd, na.rm=TRUE), ' cells found in ', design$condition[n], '\n')
  
  ## not used the emptyDrops too slow 
  # eout = emptyDrops(count.data, lower = 200)
  # eout$FDR[is.na(eout$FDR)] = 1
  # iscell_ed = eout$FDR<=0.01
  # sum(iscell_ed, na.rm=TRUE)
  
  meta = data.frame(row.names = colnames(count.data), condition = rep(design$condition[n], ncol(count.data)),
                    iscell_dd = iscell_dd)
  
  # plot rankings for number of UMI
  br.out <- barcodeRanks(count.data)

  pdf(paste0(resDir, "/UMIrank_emptyDrop_", design$condition[n], "_", design$sampleID[n],  ".pdf"),
      height = 6, width =10, useDingbats = FALSE)

  plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")

  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
  abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
  abline(v = sum(iscell_dd), col = 'darkgreen', lwd = 2.0)
  abline(v = c(8000, 10000, 12000), col = 'gray')
  text(x = c(8000, 10000, 12000), y =10000, labels = c( 8000, 10000, 12000),
       col = 'red')
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
         legend=c("knee", "inflection"))

  dev.off()
  
  
  # use defaultDrop to select cells.
  aa = CreateSeuratObject(counts = count.data[, iscell_dd],
                          meta.data = meta[iscell_dd, ], 
                          min.cells = 20, min.features = 100)
  aa$cell.id = paste0(colnames(aa), '_', design$condition[n], '_', design$sampleID[n])
  #aa$cell.id = paste0(colnames(aa), '_', design$condition[n], '_', design$sampleID[n])
  
  if(n == 1) {
    mnt = aa
  }else{
    #mnt = merge(mnt, aa, add.cell.ids = c("", design$condition[n]), collapse = TRUE)
    #mnt = merge(mnt, aa, add.cell.ids = c("", design$condition[n]))
    mnt = merge(mnt, aa)
  }
  
}

#levels = c('ctl_dpi0', 'ctl_dpi4',  'ctl_dpi7', 'ctl_dpi14')
aa$condition = factor(aa$condition, levels =  levels)

mnt[["percent.mt"]] <- PercentageFeatureSet(mnt, pattern = "^mt-")


saveRDS(mnt, 
     file = paste0(RdataDir, 
                   'seuratObject_design_variableGenes_', species, version.analysis, '.rds'))

rm(count.data); rm(br.out); rm(bc);rm(meta)

##########################################
# QCs and cell filtering
##########################################
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_design_variableGenes_', species, version.analysis, '.rds'))


pdfname = paste0(resDir, '/QCs_nCounts_nFeatures_percentMT.pdf')
pdf(pdfname, width=16, height = 8)

table(aa$condition) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(condition = factor(Var1, levels=levels)) %>%
  #mutate(cellNbs = integer(Freq))
  ggplot(aes(x=condition, y=Freq)) +
  geom_bar(stat="identity", width=0.5) +
  theme_classic() +
  labs( x = '', y = 'detected cell # from cellRanger barcodes' )  +
  theme(axis.text.x = element_text(angle = 90, size = 10)) + 
  geom_hline(yintercept = c(3000, 5000, 7000), col = 'red')

Idents(aa) = factor(aa$condition, levels = levels)

VlnPlot(aa, features = 'nFeature_RNA', y.max = 10000)
VlnPlot(aa, features = 'nCount_RNA', y.max = 100000)
VlnPlot(aa, features = 'percent.mt', y.max =20)


FeatureScatter(aa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(aa, feature1 = "nCount_RNA", feature2 = "percent.mt")

dev.off()


##########################################
## filter cells here
##########################################
VlnPlot(aa, features = 'nFeature_RNA', y.max = 10000, pt.size = 0.001) +
  geom_hline(yintercept = c(1000, 8000), col = 'red')

VlnPlot(aa, features = 'nCount_RNA', y.max = 100000) +
  geom_hline(yintercept = c(500, 8000), col = 'red')

VlnPlot(aa, features = 'percent.mt', y.max =20) + 
  geom_hline(yintercept = c(5), col = 'red')

aa <- subset(aa, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & percent.mt < 5)


##########################################
# first umap and clustering  
##########################################
aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000)
#all.genes <- rownames(aa)

aa <- ScaleData(aa)

aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE)
ElbowPlot(aa, ndims = 30)


aa$condition = factor(aa$condition, levels = levels)
Idents(aa) = aa$condition

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 50, min.dist = 0.1)
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_umap_', species, version.analysis, 
                          '.rds'))


##########################################
# overview marker genes 
##########################################
aa = readRDS(file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_umap_', species, version.analysis, 
                           '.rds'))

markers = c('Adgre1', 'Cd14', 'Cd68', 'Cx3cr1', 'Itgam-ENSMUSG00000030786',"Itgam-ENSMUSG00000108596", # pan macrophage
            'Tlr2', 'Nos2', 'Cd80', 'Cd86', 'Ifng', # M1
            'Arg1', 'Cd163', 'Il4', 'Irf4', 'Mrc1', # M2
            'En1', 'Ddit4', 'Ldha', 'Eno1', 'Serpine1', 'Lgals1', 'Hif1a', # protomyfibro- and Myofibro
            'Aldh1a3', 'Rdh10', 
            'Sfrp2', 'Cthrc1', 'Fstl1', #SFRP+ fibroblasts
            'Tnc', 'Stat3', 'Pcsk5', 'Cyp26b1', #Proto-Myofibroblasts
            'Acta2', 'Postn', 'Lrrc15', 'Runx2', # Myofibroblasts 
            'Pdpn', 'Ccl2', 'Cxcl1', 'Ccl11', 'Ccl7', 'Ccl8', #Pro-inflammatory
            'Dpt', 'Pi16', #Universal fibroblast markers and fascia
            'Ly6a', 'Procr', 'Plac8', #Fascia
            'Mgp', 'Cygb', 'Cxcl12',  # Reticular
            'Sparc', 'Dcn', 'Lum' #Papillary
)

mm = match(markers, rownames(aa))
markers[which(is.na(mm))]

pdfname = paste0(resDir, '/FeaturePlots_markerGenes.pdf')
pdf(pdfname, width=8, height = 6)

for(n in 1:length(markers))
#for(n in 1:5)
{
  cat(n, '\n')
  p = FeaturePlot(aa, features = markers[n])
  plot(p)
  
}

dev.off()


########################################################
########################################################
# Section II: identifying doublet 
# 
########################################################
########################################################
library(DoubletFinder)

aa = readRDS(file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_umap_', species, version.analysis, 
                           '.rds'))

aa$condition = factor(aa$condition, levels = levels)
Idents(aa) = aa$condition

cc = unique(aa$condition)

for(n in 1:length(cc))
{
  # n = 1
  subs <- subset(aa, condition == cc[n])
  
  subs <- FindVariableFeatures(subs, selection.method = "vst", nfeatures = 3000)
  subs <- ScaleData(subs)
  
  subs <- RunPCA(subs, verbose = TRUE)
  subs <- FindNeighbors(subs, dims = 1:30)
  subs <- FindClusters(subs, resolution = 1)
  
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
  
  subs <- doubletFinder_v3(subs, PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi.adj,  
                           reuse.pANN = FALSE, sct = FALSE)
  
  df_out = subs@meta.data
  subs$DF_out = df_out[, grep('DF.classification', colnames(df_out))]
  
  DimPlot(subs, label = TRUE, repel = TRUE, group.by = 'DF_out',
          raster=FALSE)
  ggsave(filename = paste0(resDir, '/subs_doubletFinder_out_', cc[n], '_', species, version.analysis, '.pdf'), 
         width = 12, height = 8)
  
  saveRDS(subs, file = paste0(RdataDir, 'subs_doubletFinder_out_', cc[n], '_', species, version.analysis,  
                              '.rds'))
  
}

##########################################
# save the doubletFinder in the main table  
##########################################
cc = unique(aa$condition)
aa$DF_out = NA

for(n in 1:length(cc))
{
  # n = 1
  cat(n, '--', cc[n], '\n')
  subs = readRDS(file =  paste0(RdataDir, 'subs_doubletFinder_out_', cc[n], '_', species, version.analysis,  
                                '.rds'))
  aa$DF_out[match(colnames(subs), colnames(aa))] = subs$DF_out
  
}


DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'DF_out', raster=FALSE)
ggsave(filename = paste0(resDir, '/umap_doubletFinder_results.pdf'), width = 12, height = 8)

VlnPlot(aa, features = 'nCount_RNA', group.by = 'DF_out', pt.size = 0.1) +
  geom_hline(yintercept = c(25000, 12500), col = 'red')

ggsave(filename = paste0(resDir, '/nCounts_RNA_double.vs.singlet_doubletFinder_results.pdf'), 
       width = 12, height = 8)

VlnPlot(aa, features = 'nFeature_RNA', group.by = 'DF_out', pt.size = 0.1) +
  geom_hline(yintercept = c(2000, 5000), col = 'red')

ggsave(filename = paste0(resDir, '/nFeatures_RNA_double.vs.singlet_doubletFinder_results.pdf'), 
       width = 12, height = 8)


as_tibble(data.frame(condition = aa$condition, group= aa$DF_out)) %>%
  group_by(condition, group) %>% tally() 

pcts = c()
for(n in 1:length(cc))
{
  # n =1
  pcts = c(pcts, length(which(aa$DF_out== 'Doublet' & 
                                aa$condition == cc[n]))/length(which(aa$condition == cc[n])))
  
}

data.frame(condition = cc, pct = pcts) %>%
  ggplot(aes(x = condition, y = pct, fill = condition)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none")  + 
  ggtitle('pct of doublets by DF ') + 
  theme(axis.text.x = element_text(angle = 90)) 

ggsave(filename = paste0(resDir, '/Percentages_doublet.vs.total_doubletFinder_results.pdf'), 
       width = 8, height = 6)

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_doubletFinderOut_', 
                          species, version.analysis, '.rds'))


##########################################
# discard the doublet and redo umap and clustering 
##########################################
aa = subset(aa, cells = colnames(aa)[which(aa$DF_out == 'Singlet')])

aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000)
#all.genes <- rownames(aa)

aa <- ScaleData(aa)

aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE)
ElbowPlot(aa, ndims = 30)


aa$condition = factor(aa$condition, levels = levels)
Idents(aa) = aa$condition

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 50, min.dist = 0.1)
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)


## clustering 
ElbowPlot(aa, ndims = 30)
aa <- FindNeighbors(aa, dims = 1:20)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)
DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE)

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE)

p1 + p2

ggsave(filename = paste0(resDir, '/first_UMAP_clusters.pdf'), 
       width = 16, height = 6)


markers = c('Adgre1', 'Cd14', 'Cd68', 'Cx3cr1', 'Itgam-ENSMUSG00000030786',"Itgam-ENSMUSG00000108596", # pan macrophage
            'Tlr2', 'Nos2', 'Cd80', 'Cd86', 'Ifng', # M1
            'Arg1', 'Cd163', 'Il4', 'Irf4', 'Mrc1', # M2
            'En1', 'Ddit4', 'Ldha', 'Eno1', 'Serpine1', 'Lgals1', 'Hif1a', # protomyfibro- and Myofibro
            'Aldh1a3', 'Rdh10', 
            'Sfrp2', 'Cthrc1', 'Fstl1', #SFRP+ fibroblasts
            'Tnc', 'Stat3', 'Pcsk5', 'Cyp26b1', #Proto-Myofibroblasts
            'Acta2', 'Postn', 'Lrrc15', 'Runx2', # Myofibroblasts 
            'Pdpn', 'Ccl2', 'Cxcl1', 'Ccl11', 'Ccl7', 'Ccl8', #Pro-inflammatory
            'Dpt', 'Pi16', #Universal fibroblast markers and fascia
            'Ly6a', 'Procr', 'Plac8', #Fascia
            'Mgp', 'Cygb', 'Cxcl12',  # Reticular
            'Sparc', 'Dcn', 'Lum' #Papillary
)

mm = match(markers, rownames(aa))
markers[which(is.na(mm))]

pdfname = paste0(resDir, '/FeaturePlots_markerGenes.pdf')
pdf(pdfname, width=8, height = 6)

for(n in 1:length(markers))
  #for(n in 1:5)
{
  cat(n, '\n')
  p = FeaturePlot(aa, features = markers[n])
  plot(p)
  
}

dev.off()


saveRDS(aa, file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_DFout_umapClustering_', 
                          species, version.analysis, '.rds'))


##########################################
# calculate the cell cycle scoring  
##########################################
aa = readRDS(file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_DFout_umapClustering_', 
                           species, version.analysis, '.rds'))

# Assign Cell-Cycle Scores
s.genes <- firstup(cc.genes$s.genes)
g2m.genes <- firstup(cc.genes$g2m.genes)

aa <- CellCycleScoring(aa, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
#aa <- RunPCA(aa, features = c(s.genes, g2m.genes))

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)

ggsave(filename = paste0(resDir, '/UMAP_cellcyclePhase.pdf'), 
       width = 8, height = 6)


saveRDS(aa, file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_DFout_cellCycle_umapClustering_', 
                          species, version.analysis, '.rds'))



########################################################
########################################################
# Section III: annotate subtypes of macrophage and fibroblast
# 
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_DFout_cellCycle_umapClustering_', 
                           species, version.analysis, '.rds'))


markers = c('Adgre1', 'Cd14', 'Cd68', 'Cx3cr1', 'Itgam-ENSMUSG00000030786',"Itgam-ENSMUSG00000108596", # pan macrophage
            'Tlr2', 'Nos2', 'Cd80', 'Cd86', 'Ifng', # M1
            'Arg1', 'Cd163', 'Il4', 'Irf4', 'Mrc1', # M2
            
            'En1', 'Ddit4', 'Ldha', 'Eno1', 'Serpine1', 'Lgals1', 'Hif1a', # protomyfibro- and Myofibro
            'Aldh1a3', 'Rdh10', 
            'Sfrp2', 'Cthrc1', 'Fstl1', #SFRP+ fibroblasts
            'Tnc', 'Stat3', 'Pcsk5', 'Cyp26b1', #Proto-Myofibroblasts
            'Acta2', 'Postn', 'Lrrc15', 'Runx2', # Myofibroblasts 
            'Pdpn', 'Ccl2', 'Cxcl1', 'Ccl11', 'Ccl7', 'Ccl8', #Pro-inflammatory
            'Dpt', 'Pi16', #Universal fibroblast markers and fascia
            'Ly6a', 'Procr', 'Plac8', #Fascia
            'Mgp', 'Cygb', 'Cxcl12',  # Reticular
            'Sparc', 'Dcn', 'Lum' #Papillary
)

mm = match(markers, rownames(aa))
markers[which(is.na(mm))]

aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.3)
DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE)

ggsave(filename = paste0(resDir, '/coarseClusters_RNA_snn_res_16Clusters.pdf'), 
       width = 8, height = 6)

markers = FindAllMarkers(aa, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)

saveRDS(markers, file = paste0(RdataDir, 'seuratObject_', species, version.analysis, 
                               '_markers_RNA_snn_res_16Clusters.rds'))

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
  #slice_max(n = 10, order_by = avg_log2FC) -> top10

#saveRDS(top10, file = paste0(RdataDir, 'top10_markerGenes_coarseCluster.rds'))

all.genes <- rownames(aa)
aa <- ScaleData(aa, features = all.genes)

xx = subset(aa, downsample = 1000)

DoHeatmap(xx, features = top10$gene) + NoLegend()

ggsave(filename = paste0(resDir, '/first_test_RNA_snn_res_16Clusters_clusterMarkers.pdf'), 
       width = 20, height = 20)




##########################################
# using known marker genes to annotate clusters
##########################################
aa$celltypes = NA

FeaturePlot(aa, features = c('Procr', 'Dpt', 'Pi16', 'Col1a2', 'Acta2', 'Lum', 'Col3a1', 'Col1a1', 'Mmp2', 
                             'Pdgfra'))

ggsave(filename = paste0(resDir, '/FeaturePlots_MarkerGenes_Fibroblast.pdf'), 
       width = 12, height = 8)


aa$celltypes[which(!is.na(match(aa$seurat_clusters, c(10, 7, 5, 0, 1, 3, 4))))] = 'fibroblast'


FeaturePlot(aa, features = c('Adgre1', 'Cd68', 'Itgam-ENSMUSG00000030786', 
                             'Cd14',  'Cx3cr1', # pan macrophage
                             'Tlr2', 'Nos2', 'Cd80', 'Cd86', 'Ifng', # M1 (pro-inflamatory)
                             'Arg1', 'Cd163', 'Il4', 'Irf4', 'Mrc1' # M2 (anti-)
))

ggsave(filename = paste0(resDir, '/FeaturePlots_MarkerGenes_Macrophage.pdf'), 
       width = 16, height = 10)

FeaturePlot(aa, features = c('Adgre1', 'Cd68', 'Itgam-ENSMUSG00000030786'
))

aa$celltypes[which(aa$seurat_clusters == 2| aa$seurat_clusters == 15| 
                  aa$seurat_clusters == 6)] = 'macrophage'


FeaturePlot(aa, features = c('Ly6g')) # neutrophil, likely is cluster 8 

aa$celltypes[which(!is.na(match(aa$seurat_clusters, c(8))))] = 'neutrophil'


FeaturePlot(aa, features = c('Cd3g', 'Cd3d', 'Cd3e')) # T cells, correspond to cluster 12 and cluster 9

aa$celltypes[which(!is.na(match(aa$seurat_clusters, c(9, 12))))] = 'T'

FeaturePlot(aa, features = c('Itgax')) # Dendritic cells, no clear cluster found

FeaturePlot(aa, features = c('Ncam1', 'Cd3g', 'Cd3d', 'Cd3e')) # NK cells

aa$celltypes[which(!is.na(match(aa$seurat_clusters, c(13))))] = 'NK'

FeaturePlot(aa, features = c('Pecam1', 'Vwf')) # Endothelial cells
aa$celltypes[which(!is.na(match(aa$seurat_clusters, c(11))))] = 'Endothelial'

FeaturePlot(aa, features = c('Pax7', 'Myod1', 'Ckm')) # 

aa$celltypes[which(!is.na(match(aa$seurat_clusters, c(14))))] = 'skeletalMuscle_others'


DimPlot(aa, group.by = 'celltypes', label = TRUE, repel = TRUE, raster=FALSE)

ggsave(filename = paste0(resDir, '/CoarseCluster_annotation_v1.pdf'), 
       width = 8, height = 6)

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_DFout_cellCycle_umapClustering_',
                          'celltypeAnnot.v1_',  species, version.analysis, '.rds'))


##########################################
# cell type annotation from David Saini
##########################################
sd = readRDS(file = paste0("../data/cell_annotation_mouseSkin_scRNAseq/", 
                           "WT.Skin.rds"))

DimPlot(sd, group.by = 'cell_type_coarse', label = TRUE, repel = TRUE)

DimPlot(sd, group.by = 'seurat_clusters', label = TRUE, repel = TRUE)


DimPlot(sd, group.by = 'cell_type_fine', label = TRUE, repel = TRUE)

DimPlot(sd, group.by = 'nn_cell', label = TRUE, repel = TRUE)

ggsave(filename = paste0(resDir, '/mouseSkin_WT_subtypes_annotation_v1.pdf'), 
       width = 12, height = 8)

sd$subtypes = sd$nn_cell

annots = sd@meta.data

saveRDS(annots, file = paste0(RdataDir, 'mouseSkin_WT_subtypes_annotation_v1.rds'))

rm(sd)

p1 = DimPlot(sd, group.by = 'cell_type_coarse', label = TRUE, repel = TRUE)

p2 = DimPlot(sd, group.by = 'cell_type_fine', label = TRUE, repel = TRUE)

p1 + p2
   

## annotate processed data 
aa = readRDS(file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_DFout_cellCycle_umapClustering_',
                           'celltypeAnnot.v1_',  species, version.analysis, '.rds'))

annots = readRDS(file = paste0(RdataDir, 'mouseSkin_WT_subtypes_annotation_v1.rds'))


aa$sampleID = NA
aa$sampleID[which(aa$condition == 'ctl_dpi0')] = 'SID310697'
aa$sampleID[which(aa$condition == 'ctl_dpi4')] = 'SID228615'
aa$sampleID[which(aa$condition == 'ctl_dpi7')] = 'SID228613'
aa$sampleID[which(aa$condition == 'ctl_dpi14')] = 'SID228611'

aa$cellID = sapply(aa$cell.id, function(x){return(unlist(strsplit(as.character(x), '_'))[1])})
aa$cellID = paste0(aa$sampleID, '_', aa$cellID)

mm = match(aa$cellID, rownames(annots))
aa$cell_type_coarse = annots$cell_type_coarse[mm]
aa$cell_type_fine = annots$cell_type_fine[mm]
aa$subtypes = annots$subtypes[mm]


aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000)
aa <- ScaleData(aa)

aa <- RunPCA(aa, verbose = FALSE, weight.by.var = FALSE)
ElbowPlot(aa, ndims = 30)

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.3)
DimPlot(aa, group.by = 'subtypes', label = TRUE, repel = TRUE)

p1 = DimPlot(aa, group.by = 'celltypes', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'subtypes', label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(resDir, '/mouseSkin_WT_subtypes_borrowedFromDS.pdf'), 
       width = 16, height = 8)


########################################################
########################################################
# Section IV: analyze the macrophage and fibroblast
# 
########################################################
########################################################
##########################################
# import directly the macrophage subtypes from David and Sebastian 
##########################################
sd = readRDS(file = paste0("../data/cell_annotation_mouseSkin_scRNAseq/", 
                           "WT.Skin.Macrophages.rds"))

p1 = DimPlot(sd, group.by = 'stage.m', label = TRUE, repel = TRUE)
p2 = DimPlot(sd, group.by = 'condition', label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(resDir, '/mouseSkin_WT_subtypes_borrowedFromDS.pdf'), 
       width = 16, height = 8)



# DimPlot(sd, group.by = 'cell_type_fine', label = TRUE, repel = TRUE)
# 
# DimPlot(sd, group.by = 'nn_cell', label = TRUE, repel = TRUE)
# 
# sd = subset(sd, cells = colnames(sd)[which(sd$cell_type_fine == 'MacDC'| sd$cell_type_fine == 'Macrophage')])
# 
# sd <- FindVariableFeatures(sd, selection.method = "vst", nfeatures = 3000)
# sd <- ScaleData(sd)
# 
# sd <- RunPCA(sd, verbose = FALSE, weight.by.var = FALSE)
# ElbowPlot(sd, ndims = 30)
# 
# sd <- RunUMAP(sd, dims = 1:20, n.neighbors = 30, min.dist = 0.3)
# 
# DimPlot(sd, group.by = 'nn_cell', label = TRUE, repel = TRUE)


saveRDS(sd, file = paste0(RdataDir, '/mouse_skin_macrophage_subtypes_SD.rds'))


##########################################
# subset macrophages
##########################################
aa = readRDS(file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_DFout_cellCycle_umapClustering_',
                     'celltypeAnnot.v1_',  species, version.analysis, '.rds'))


FeaturePlot(aa, features = c('Adgre1', 'Cd68', 'Itgam-ENSMUSG00000030786', 'Csf1r', 'H2-Ab1', 'Mertk',
                             'Cd14',  'Cx3cr1', # pan macrophage
                             'Tlr2', 'Nos2', 'Cd80', 'Cd86', 'Ifng', # M1 (pro-inflamatory)
                             'Arg1', 'Cd163', 'Il4', 'Irf4', 'Mrc1' # M2 (anti-)
))

ggsave(filename = paste0(resDir, '/CoarseCluster_annotation_macrophages_markGenes.pdf'), 
       width = 12, height = 10)

aa$clusters = aa$seurat_clusters

VlnPlot(aa, group.by = 'clusters', 
        features = c('Adgre1', 'Cd68', 'Itgam-ENSMUSG00000030786', 'Csf1r', 'H2-Ab1', 'Mertk',
                                    'Cd14',  'Cx3cr1', # pan macrophage
                                    'Tlr2', 'Nos2', 'Cd80', 'Cd86', 'Ifng', # M1 (pro-inflamatory)
                                    'Arg1', 'Cd163', 'Il4', 'Irf4', 'Mrc1' # M2 (anti-)),
                         
))


subs = subset(aa, cells = colnames(aa)[which(aa$celltypes == "macrophage")])

subs = NormalizeData(subs, normalization.method = "LogNormalize", scale.factor = 10000)
subs <- FindVariableFeatures(subs, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs

subs <- ScaleData(subs, features = rownames(subs))
subs <- RunPCA(subs, features = VariableFeatures(object = subs), verbose = FALSE, weight.by.var = FALSE)

ElbowPlot(subs, ndims = 50)

subs <- FindNeighbors(subs, dims = 1:30)
subs <- FindClusters(subs, verbose = FALSE, algorithm = 3, resolution = 0.5)

subs <- RunUMAP(subs, reduction = "pca", dims = 1:30, n.neighbors = 30,  min.dist = 0.3)

p1 = DimPlot(subs, group.by = 'condition', label = TRUE, repel = TRUE)
p2 = DimPlot(subs, group.by = 'seurat_clusters', label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(resDir, '/Macrophages_subtypes_clustering.pdf'), 
       width = 12, height = 10)


VlnPlot(subs, group.by = 'seurat_clusters', 
        features = c('Adgre1', 'Cd68', 'Itgam-ENSMUSG00000030786', 'Csf1r', 'H2-Ab1', 'Mertk',
                     'Cd14',  'Cx3cr1' # pan macrophage
                     
        ))

## remove two small clusters
subs = subset(subs, cells= colnames(subs)[which(subs$seurat_clusters != 10 & subs$seurat_clusters != 11)])

subs = NormalizeData(subs, normalization.method = "LogNormalize", scale.factor = 10000)
subs <- FindVariableFeatures(subs, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs

subs <- ScaleData(subs, features = rownames(subs))
subs <- RunPCA(subs, features = VariableFeatures(object = subs), verbose = FALSE, weight.by.var = FALSE)

ElbowPlot(subs, ndims = 50)

subs <- FindNeighbors(subs, dims = 1:20)
subs <- FindClusters(subs, verbose = FALSE, algorithm = 3, resolution = 0.5)

subs <- RunUMAP(subs, reduction = "pca", dims = 1:30, n.neighbors = 50,  min.dist = 0.3)

p1 = DimPlot(subs, group.by = 'condition', label = TRUE, repel = TRUE)
p2 = DimPlot(subs, group.by = 'seurat_clusters', label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(resDir, '/Macrophages_subtypes_clustering.pdf'), 
       width = 12, height = 6)

DimPlot(subs, group.by = 'Phase')

DimPlot(subs, group.by = 'subtypes')

subs$subtypes = NA
subs$subtypes[which(subs$seurat_clusters == 10)] = 'cycling'

markers = FindAllMarkers(subs, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
#slice_max(n = 10, order_by = avg_log2FC) -> top10

#saveRDS(top10, file = paste0(RdataDir, 'top10_markerGenes_coarseCluster.rds'))

all.genes <- rownames(subs)
subs <- ScaleData(subs, features = all.genes)

DoHeatmap(subs, features = top10$gene) + NoLegend()

ggsave(filename = paste0(resDir, '/macrophages_subclusters_clusterMarkers.pdf'), 
       width = 16, height = 20)

# VlnPlot(subs, group.by = 'seurat_clusters', 
#         features = c('Adgre1', 'Cd68', 'Itgam-ENSMUSG00000030786', 'Csf1r', 'H2-Ab1', 'Mertk',
#                      'Cd14',  'Cx3cr1' # pan macrophage
#                      
#         ))
# 
# subs = subset(subs, cells = colnames(subs)[which(subs$seurat_clusters != 10)])

subs$clusters = subs$seurat_clusters
kk = which(subs$clusters != 9)
subs$subtypes[kk] = paste0('C', subs$clusters[kk], '.mm')

DimPlot(subs, group.by = 'subtypes', label = TRUE, repel = TRUE)

saveRDS(subs, file = paste0(RdataDir, '/mouse_skin_macrophage_subtypes.rds'))

