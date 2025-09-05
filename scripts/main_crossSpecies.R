##########################################################################
##########################################################################
# Project: Collaboration with Sabine and Sebastian from Cologne
# Script purpose: cross-species analysis for axolotl and mouse
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Sep  3 13:55:39 2025
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
library(SeuratDisk)

version.analysis = '_crossSpecies_20250903'
resDir = paste0("../results/scRNAseq_analysis_immune", version.analysis, '/')
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '/groups/tanaka/People/current/jiwang/projects/bone_healing_CSD/fromTobie/CSD_batch1_batch2/'

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/'
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


########################################################
########################################################
# Section II: test Seurat-based integration methods
# 
########################################################
########################################################

##########################################
# import ax, nm, mm data 
##########################################
ax = readRDS(file = paste0(RdataDir, 'axoltol_limbBlatema_batch1_macrophage_subtypes.rds'))
mm = readRDS(file = paste0(RdataDir, 'mouse_skin_macrophage_subtypes.rds'))

ax$species = 'ax'
ax$time = droplevels(ax$time)
ax$subtypes[which(ax$subtypes == 'cycling')] = 'cycling.ax'

mm$species = 'mm'
mm$subtypes[which(mm$subtypes == 'cycling')] = 'cycling.mm'


## define the one-on-one ortholog between axololt and mice
an_orthologs = data.frame(ref = rownames(ax), query = rownames(ax))
rownames(an_orthologs) = an_orthologs$ref
an_orthologs$query = sapply(an_orthologs$query, 
                            function(x){firstup(unlist(strsplit(as.character(x), '-'))[1])})

jj = which(!is.na(match(an_orthologs$query, rownames(mm))))
an_orthologs = an_orthologs[jj, ]

counts = table(an_orthologs$query)
gg_uniq = names(counts)[which(counts == 1)]
jj2 = which(!is.na(match(an_orthologs$query, gg_uniq)))

an_orthologs = an_orthologs[jj2, ]

ax = subset(ax, features = an_orthologs$ref)
aa = subset(mm, features = an_orthologs$query)

rm(mm)

counts = ax@assays$RNA@layers$counts
metadata = ax@meta.data
counts = counts[match(an_orthologs$ref, rownames(ax)), ]
rownames(counts) = an_orthologs$query

new_ax <- CreateSeuratObject(counts=counts, assay = 'RNA', meta.data = metadata)
new_ax<- NormalizeData(new_ax, normalization.method = "LogNormalize", scale.factor = 10000)
new_ax <- FindVariableFeatures(new_ax, selection.method = "vst", nfeatures = 5000)

new_ax <- ScaleData(new_ax, features = rownames(new_ax))

aa = merge(aa, y = new_ax, add.cell.ids = c("m", "ax"), project = "skinRepair")

rm(list = c('counts', 'metadata', 'ax'))

rm(new_ax)

saveRDS(aa, file = paste0(RdataDir, 'mm_ax_scRNAseq_merged_forSeurat_v1.rds'))

##########################################
# save files for scVI and scANVI 
##########################################


###### WORKAROUND ###### from https://github.com/mojaveazure/seurat-disk/issues/147
# # assigning the previous version of the `[[` function for the Assay class to the SeuratDisk package environment
# "[[.Assay" <- function(x, i, ..., drop = FALSE) {
#   if (missing(x = i)) {
#     i <- colnames(x = slot(object = x, name = 'meta.features'))
#   }
#   data.return <- slot(object = x, name = 'meta.features')[, i, drop = FALSE, ...]
#   if (drop) {
#     data.return <- unlist(x = data.return, use.names = FALSE)
#     names(x = data.return) <- rep.int(x = rownames(x = x), times = length(x = i))
#   }
#   return(data.return)
# }
# environment(`[[.Assay`) <- asNamespace("SeuratObject")
# rlang::env_unlock(asNamespace("SeuratDisk"))
# assign("[[.Assay", `[[.Assay`, asNamespace("SeuratDisk"))
# lockEnvironment(asNamespace("SeuratDisk"), bindings = TRUE)
# rm(`[[.Assay`)

###### WORKAROUND ######
aa = readRDS(file = paste0(RdataDir, 'mm_ax_scRNAseq_merged_forSeurat_v1.rds'))

ax = subset(aa, cells = colnames(aa)[which(aa$species == 'ax')])
counts_ax = ax@assays$RNA@layers$counts.2
metadata_ax = ax@meta.data

mm = subset(aa, cells = colnames(aa)[which(aa$species == 'mm')])
counts_mm = mm@assays$RNA@layers$counts.1
metadata_mm = mm@meta.data

genes = rownames(aa)

save(counts_ax, metadata_ax, counts_mm, metadata_mm, genes, 
     file = paste0(resDir, 'metadata_counts_ax_mm.Rdata'))

## back to Seurat v4
load(file = paste0(resDir, 'metadata_counts_ax_mm.Rdata'))

colnames(counts_ax) = rownames(metadata_ax)
rownames(counts_ax) = genes
ax <- CreateSeuratObject(counts=counts_ax, assay = 'RNA', meta.data = metadata_ax)

colnames(counts_mm) = rownames(metadata_mm)
rownames(counts_mm) = genes
mm <- CreateSeuratObject(counts=counts_mm, assay = 'RNA', meta.data = metadata_mm)

aa = merge(mm, ax)

saveFile = paste0(resDir, 'mm_ax_scRNAseq_merged_seuratV4.h5Seurat')
SaveH5Seurat(aa, filename = saveFile, overwrite = TRUE)
Convert(saveFile, dest = "h5ad", overwrite = TRUE)


##########################################
# test Seurat 
##########################################
source(paste0(functionDir, 'functions_dataIntegration.R'))
aa = readRDS(file = paste0(RdataDir, 'nm_mm_ax_scRNAseq_merged_forSeurat_v1.rds'))


aa = ScaleData(aa, features = rownames(aa))

method = 'Harmony'

if(method == 'Seurat_RPCA'){
  ref.combined = IntegrateData_Seurat_RPCA(aa, 
                                           group.by = 'species', 
                                           nfeatures = 5000,
                                           #merge.order = matrix(c(-2, 1, -3, -1), ncol = 2),
                                           redo.normalization.scaling = FALSE,
                                           correct.all = FALSE)
  
  p1 = DimPlot(ref.combined, group.by = 'subtypes', label = TRUE, repel = TRUE, raster=FALSE) + 
    ggtitle('Seurat_RPCA')
  p2 = DimPlot(ref.combined, group.by = 'species', label = TRUE, repel = TRUE) +
    ggtitle("Seurat_RPCA")
  
  p1 + p2
  
  ggsave(filename = paste0(outDir, '/cross_species_mapping_Seurat_RPCA.pdf'), 
         width = 16, height = 12)
  
  
  
}

if(method == 'Seurat_CCA'){
  ref.combined = IntegrateData_Seurat_CCA(aa, 
                                          group.by = 'species', 
                                          nfeatures = 5000,
                                          #merge.order = matrix(c(-2, 1, -3, -1), ncol = 2),
                                          redo.normalization.scaling = FALSE,
                                          correct.all = FALSE)
  
  p1 = DimPlot(ref.combined, group.by = 'subtypes', label = TRUE, repel = TRUE, raster=FALSE) + 
    ggtitle('Seurat_CCA')
  p2 = DimPlot(ref.combined, group.by = 'species', label = TRUE, repel = TRUE) +
    ggtitle("Seurat_CCA")
  
  p1 + p2
  
  ggsave(filename = paste0(resDir, '/cross_species_mapping_Seurat_CCA.pdf'), 
         width = 16, height = 6)
  
  saveRDS(ref.combined, file = paste0(resDir, 'crossSpecies_mm_ax_scRNAseq_SeuratCCA.rds'))
  
}


if(method == 'Harmony'){
  source(paste0(functionDir, 'functions_dataIntegration.R'))
  ref.combined = IntegrateData_runHarmony(aa, 
                                          group.by = 'species',
                                          nfeatures = 5000,
                                          dims.use = c(1:50),
                                          redo.normalization.hvg.scale.pca = TRUE,
                                          max.iter.harmony = 30,
                                          epsilon.harmony = -Inf
                                          #correct.all = FALSE
  )
  
  p1 = DimPlot(ref.combined, group.by = 'subtypes', label = TRUE, repel = TRUE, raster=FALSE) + 
    ggtitle(method)
  p2 = DimPlot(ref.combined, group.by = 'species', label = TRUE, repel = TRUE) +
    ggtitle(method)
  
  p1 + p2
  
  ggsave(filename = paste0(resDir, '/cross_species_mapping_Harmony.pdf'), 
         width = 16, height = 6)
  
  saveRDS(ref.combined, file = paste0(resDir, 'mm_ax_scRNAseq_crossSpecies_Harmony.rds'))
  
  
}

########################################################
########################################################
# Section II : test neighboorhood analysis
# 
########################################################
########################################################
# Load packages
suppressPackageStartupMessages(library(scrabbitr))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(miloR))
suppressPackageStartupMessages(library(DelayedArray))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggraph))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(ggalluvial))
suppressPackageStartupMessages(library(ggrepel))
library(scuttle)
library(Seurat)
source(paste0(functionDir, 'functions_scRNAseq.R'))
source(paste0(functionDir, 'functions_Visium.R'))
source(paste0(functionDir, 'functions_cccInference.R'))

##########################################
# load ax and mm data and process them 
##########################################
ax = readRDS(file = paste0(RdataDir, 'axoltol_limbBlatema_batch1_macrophage_subtypes.rds'))
mm = readRDS(file = paste0(RdataDir, 'mouse_skin_macrophage_subtypes.rds'))

mm = subset(mm, cells = colnames(mm)[which(!is.na(mm$subtypes))])

ax$species = 'ax'
ax$time = droplevels(ax$time)
ax$subtypes[which(ax$subtypes == 'cycling')] = 'cycling.ax'

ax = as.SingleCellExperiment(ax, assay = 'RNA')

mm$species = 'mm'
mm$subtypes[which(mm$subtypes == 'cycling')] = 'cycling.mm'
mm = as.SingleCellExperiment(mm, assay = 'RNA')

##########################################
# prepare the one-to-one orthologs 
##########################################
an_orthologs = data.frame(ref = rownames(ax), query = rownames(ax))
rownames(an_orthologs) = an_orthologs$ref
an_orthologs$query = sapply(an_orthologs$query, 
                            function(x){firstup(unlist(strsplit(as.character(x), '-'))[1])})

jj = which(!is.na(match(an_orthologs$query, rownames(mm))))
an_orthologs = an_orthologs[jj, ]

counts = table(an_orthologs$query)
gg_uniq = names(counts)[which(counts == 1)]
jj2 = which(!is.na(match(an_orthologs$query, gg_uniq)))

an_orthologs = an_orthologs[jj2, ]

##########################################
# compute the milo graph
##########################################
reLoad_computated_r_milo_m_miol = TRUE
if(reLoad_computated_r_milo_m_miol)
{
  r_milo = readRDS(paste0('../published_dataset/scripts/Rdata/', "r_milo.rds"))
  m_milo = readRDS(paste0('../published_dataset/scripts/Rdata/', "m_milo.rds"))
  
}else{
  # Compute rabbit neighbourhoods
  a_milo <- Milo(ax)
  
  a_milo <- buildGraph(a_milo, k=30, d=50, reduced.dim="PCA")
  
  a_milo <- makeNhoods(a_milo, prop=0.05, k=30, d=50, refined=T, reduced_dims="PCA")
  
  a_milo <- buildNhoodGraph(a_milo)
  
  a_milo
  
  # Export miloR object
  writeMM(a_milo@nhoods, paste0(resDir, "ax_nhoods.mtx"))
  saveRDS(a_milo, paste0(resDir, "ax_milo.rds"))
  
  options(repr.plot.width = 12, repr.plot.height = 5, repr.plot.res = 300)
  p1 <- scrabbitr::plotNhoodSizeHist(a_milo, colour="blue")
  
  p1
  
  p2 <- plotNhoodGraph(a_milo, size_range=c(0.1,3), node_stroke=0.1) + 
    scale_fill_viridis(name = "Nhood size", option = "viridis", direction = 1) 
  #grid.arrange(p1, p2, nrow=1)
  p1  + p2
  
  ggsave(paste0(resDir, "ax_milo_nhood_size_hist_size_graph.pdf"), width=8, height=4, dpi=300)
  
  # Compute mouse neighbourhoods
  n_milo <- Milo(mm)
  
  n_milo <- buildGraph(n_milo, k=30, d=50, reduced.dim="PCA")
  
  n_milo <- makeNhoods(n_milo, prop=0.05, k=30, d=50, refined=T, reduced_dims="PCA")
  
  n_milo <- buildNhoodGraph(n_milo)
  
  
  # Plot nhoods and size distribution
  options(repr.plot.width = 12, repr.plot.height = 5, repr.plot.res = 300)
  
  p1 <- scrabbitr::plotNhoodSizeHist(n_milo, colour="red")
  #ggsave("../plots/compare_nhoods/m_milo_nhood_size_hist.pdf", p1, width=4, height=4, dpi=300)
  
  p2 <- plotNhoodGraph(n_milo, size_range=c(0.1,3),node_stroke=0.1) + 
    scale_fill_viridis(name = "Nhood size", option = "viridis", direction=1)
  #ggsave("../plots/compare_nhoods/m_nhood_size_graph.pdf", p2, width=6, height=5, dpi=300)
  
  #grid.arrange(p1, p2, nrow = 1)
  p1  + p2
  ggsave(paste0(resDir, "nm_milo_nhood_size_hist_size_graph.pdf"), width=8, height=4, dpi=300)
  
  saveRDS(n_milo, paste0(resDir, "n_milo.rds"))
  
}

# check mouse rowData
head(rowData(a_milo))
head(rowData(n_milo))

# Add mouse colData
head(colData(a_milo))


head(colData(n_milo))
#table(r_milo$somite_count)
#table(r_milo$dissection)
head(n_milo$subtypes)
#a_milo$celltype = a_milo$celltypes
#n_milo$celltype = n_milo$subtype

## Run neighbourhood comparison pipeline
# Run pipeline
out <- scrabbitr::calcNhoodSim(a_milo, n_milo, an_orthologs, 
                               sim_preprocessing="gene_spec", 
                               sim_measure="pearson",
                               hvg_join_type="intersection", 
                               max_hvgs=3000, 
                               #r_exclude = r_exclude, 
                               #m_exclude = m_exclude,
                               export_dir = resDir, 
                               verbose = TRUE)


out$nhood_sim[1:5, 1:5]

# Extract neighbourhood graph
a_graph <- nhoodGraph(a_milo)
n_graph <- nhoodGraph(n_milo)


# Add nhood attributes to igraph
a_nhoodIDs <- as.numeric(vertex_attr(a_graph)$name) 
a_indCells <- colnames(a_milo)[a_nhoodIDs]

V(a_graph)$cell_name <- a_indCells
V(a_graph)$celltype <- colData(a_milo)[a_indCells, "subtypes"]

n_nhoodIDs <- as.numeric(vertex_attr(n_graph)$name) 
n_indCells <- colnames(n_milo)[n_nhoodIDs]

V(n_graph)$cell_name <- n_indCells
V(n_graph)$celltype <- colData(n_milo)[n_indCells, "subtypes"]


# Calculate maximum correlations  
a_maxNhoods <- getMaxMappings(out$nhood_sim, 1, long_format=FALSE) # rabbit-mouse
n_maxNhoods <- getMaxMappings(out$nhood_sim, 2, long_format=FALSE) # mouse-rabbit
df_simFilt <- rbind(a_maxNhoods, n_maxNhoods)                                                  

options(repr.plot.width = 18, repr.plot.height = 8, repr.plot.res = 300)
p1 <- plotNhoodMaxSim(a_milo, a_maxNhoods)
p2 <- plotNhoodMaxSim(n_milo, n_maxNhoods)
p1 + p2

ggsave(paste0(resDir, "n_milo_a_milo_max_corr.pdf"), width=10, height=4, dpi=300)


##########################################
# highlight the distribution of correlation of celltypes
##########################################
celltypes <- unique(c(unique(colData(a_milo)$subtypes), unique(colData(n_milo)$subtypes)))
cols =  randomcoloR::distinctColorPalette(k = length(celltypes), altCol = FALSE, runTsne = FALSE)
cols = cols[length(cols):1]
names(cols) = celltypes

# Highlight 
options(repr.plot.width = 5.2, repr.plot.height = 5, repr.plot.res = 300)
#a_milo$isCM <- ifelse(a_milo$celltypes %in% exe_celltypes,"Extra-embryonic", "Embryonic")

p <- plotNhoodSimGroups(a_milo, a_maxNhoods$sim, 
                        group_by = "subtypes", 
                        xlabel="Correlation", ylabel="subtypes", 
                        colour_by="subtypes", 
                        group_colours = cols,
                        decreasing = FALSE,
                        size=0.15, 
                        rel_min_height=0.001, show_rank = TRUE
)

p <- p + 
  theme(text = element_text(size=14),
        axis.text = element_text(size=12), 
        axis.text.y = element_text(size=6),
        axis.ticks = element_line(size = 1),
        panel.grid.minor = element_line(size = 0.1), 
        panel.grid.major = element_line(size = 0.2))  +
  NoLegend()
p

ggsave(paste0(resDir, "ax_nhoods_celltype_ranked.pdf"), p, 
       width=5, height=8, dpi=300)


##########################################
# hightlight the projection between species
##########################################
options(repr.plot.width = 18, repr.plot.height = 8, repr.plot.res = 300)

celltypes <- unique(c(unique(colData(a_milo)$subtypes), unique(colData(n_milo)$subtypes)))
#celltypes = celltypes[grep('FB|CM|EC', celltypes)]

p_all <- plotTrajMappings(a_milo, n_milo, df_simFilt, 
                          group="subtypes", 
                          groups=celltypes, 
                          dimred="UMAP", 
                          #colour_by = cols, 
                          #rotate=90,
                          offset=c(0, 50), reflect.X=FALSE, reflect.Y=FALSE,
                          line_alpha=0.02,
                          edge_alpha=0.01, 
                          legend_pos="right") +
  guides(fill = guide_legend(override.aes = list(size=6))) + ggtitle("all")

p_all

ggsave(paste0(resDir, "axolotl_nm_mapping_test_v10.pdf"), width=18, height=10, dpi=300)


########################################################
########################################################
# Section IV: Cross-species cell type correlations 
# the same idea from https://www.science.org/doi/10.1126/science.abp9262
# Tomas and Ashley's paper
########################################################
########################################################
ax = readRDS(file = paste0(RdataDir, 'axoltol_limbBlatema_batch1_macrophage_subtypes.rds'))
mm = readRDS(file = paste0(RdataDir, 'mouse_skin_macrophage_subtypes.rds'))

ax$species = 'ax'
ax$time = droplevels(ax$time)
ax$subtypes[which(ax$subtypes == 'cycling')] = 'cycling.ax'

mm$species = 'mm'
mm$subtypes[which(mm$subtypes == 'cycling')] = 'cycling.mm'

##########################################
# find all DE genes for axolotl and mouse
##########################################
Idents(ax) = ax$subtypes
Idents(mm) = mm$subtypes
markers.ax = FindAllMarkers(ax, logfc.threshold = 0.5, min.pct = 0.1)
markers.mm = FindAllMarkers(mm, logfc.threshold = 0.5, min.pct = 0.1)

## define the one-on-one ortholog between axololt and mice
an_orthologs = data.frame(ref = rownames(ax), query = rownames(ax))
rownames(an_orthologs) = an_orthologs$ref
an_orthologs$query = sapply(an_orthologs$query, 
                            function(x){firstup(unlist(strsplit(as.character(x), '-'))[1])})

## intersect gene names with mouse data
jj = which(!is.na(match(an_orthologs$query, rownames(mm))))
an_orthologs = an_orthologs[jj, ]

## select unique genes 
#counts = table(an_orthologs$query)
gg_uniq = unique(an_orthologs$query)
jj2 = match(gg_uniq, an_orthologs$query)

an_orthologs = an_orthologs[jj2, ]

rownames(an_orthologs) = an_orthologs$query

## intersect with DE genes of axolotl and mouse
mm1 = match(an_orthologs$ref, markers.ax$gene)
mm2 = match(an_orthologs$query, markers.mm$gene)
kk = intersect(which(!is.na(mm1)), which(!is.na(mm2)))

an_orthologs = an_orthologs[kk, ]

an_orthologs$tfs = FALSE
an_orthologs$tfs[which(!is.na(match(toupper(an_orthologs$query), c(tfs))))] = TRUE

avg.ax = AggregateExpression(ax, features = an_orthologs$ref[which(an_orthologs$tfs == TRUE)], 
                             normalization.method = "LogNormalize",
                             scale.factor = 10000)
avg.ax = data.frame(avg.ax$RNA)

avg.mm = AggregateExpression(mm, features = an_orthologs$query[which(an_orthologs$tfs == TRUE)])
avg.mm = data.frame(avg.mm$RNA)


cort = psych::corr.test(avg.ax, avg.mm, method = "spearman", 
                        adjust = "fdr", alpha = 0.05, ci = F)

range(cort$r)

cort$maxrow = apply(cort$r, 1, which.max)
cort$maxcol = apply(cort$r, 2, which.max)

# cluster and order labels
hcr = hclust(dist(cort$r), method = "ward.D2")
hcc = hclust(dist(t(cort$r)), method = "ward.D2")
hcr = hcr$labels[hcr$order]
hcc = hcc$labels[hcc$order]

# reshaping the correlations
plot_df = reshape2::melt(cort$r)
plot_df$Var1 = factor(plot_df$Var1, levels = rev(hcr))
plot_df$Var2 = factor(plot_df$Var2, levels = hcc)

# add pvalue and max cor infor
plot_df$padj = -log10(reshape2::melt(cort$p.adj+min(cort$p.adj[cort$p.adj>0])/10)$value)
plot_df$rowmax = apply(Reduce(cbind, lapply(names(cort$maxrow), 
                                            function(n) plot_df$Var1==n &
                                              plot_df$Var2==colnames(cort$r)[cort$maxrow[n]])), 
                       1, any)
plot_df$colmax = apply(Reduce(cbind, lapply(names(cort$maxcol), 
                                            function(n) plot_df$Var2==n &
                                              plot_df$Var1==rownames(cort$r)[cort$maxcol[n]])), 
                       1, any)
plot_df$markcol = plot_df$value>quantile(plot_df$value, 0.98)


# getting a colourscale where 0 is white in the middle, and intensity leveled by max(abs(value))
cols = colorRampPalette(c(rev(RColorBrewer::brewer.pal(9, "Blues")),
                          RColorBrewer::brewer.pal(9, "Reds")))(101)
#cols = colorRampPalette(c(RColorBrewer::brewer.pal(9, "Reds")))(101)
br = seq((min(cort$r)), max(abs(cort$r)), length.out = 101)
cols = cols[!(br>max(cort$r) | br<min(cort$r))]

ggplot()+
  geom_point(data = plot_df, mapping = aes(x = Var2, y = Var1, fill = value), 
             shape = 21, size = 10) +
  geom_point(data = plot_df[plot_df$rowmax,], mapping = aes(x = Var2, y = Var1, size = 30), 
             shape = "—", show.legend = F, colour = "grey10")+
  geom_point(data = plot_df[plot_df$colmax,], mapping = aes(x = Var2, y = Var1, size = 30), 
             shape = "|", show.legend = F, colour = "grey10")+
  scale_x_discrete(expand = c(0,0.7)) +
  scale_y_discrete(expand = c(0,0.7)) +
  scale_fill_gradientn(breaks = signif(c(min(cort$r)+0.005, 0, max(cort$r)-0.005),2), 
                       #values = scales::rescale(c(min(br), 0, max(br))),
                       colours = cols) +
  labs(x = 'mouse', y = 'axolotl', fill = "Spearman's\nrho", size = "-log10\nadj. p-value")+
  theme_classic()+
  theme(axis.title = element_text(colour = "black", face = "bold"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8))

ggsave(paste0(resDir, "axolotl_nm_subtype_correlationAnalysis.pdf"), width=8, height=6, dpi=300)


