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
#library("viridis")
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
# Section I: test Seurat-based integration methods
# 
########################################################
########################################################

##########################################
# import ax, nm, mm data 
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

p2 = DimPlot(ax, group.by = 'cluster', label = TRUE, repel =  TRUE)

mm = readRDS(file = paste0(RdataDir, '/mouse_skin_macrophage_subtypes_SD.rds'))

mm$species = 'mm'
#mm$subtypes[which(mm$subtypes == 'cycling')] = 'cycling.mm'
mm$time = droplevels(mm$condition)
mm$subtypes = mm$stage.m
mm$cluster = mm$stage.m

p1 = DimPlot(mm, group.by = 'cluster', label = TRUE, repel =  TRUE)

p1 + p2

ggsave(filename = paste0(resDir, '/cross_species_macrophageData_mm_ax.pdf'), 
       width = 16, height = 6)


# convert a v5 assay to a v4 assay
xx = scCustomize::Convert_Assay(ax, assay = 'RNA', convert_to = 'V3')

saveRDS(xx, file = paste0(RdataDir, 
                '/axoltol_limbBlatema_batch1_macrophage_time_subtypeAnnotations_rmEpidermis_v3_seuratv3.rds'))

#ax[["RNA4"]] <- as(object = ax[["RNA"]], Class = "Assay")
#DefaultAssay(ax) = 'RNA4'

DefaultAssay(mm) = 'RNA'
mm = NormalizeData(mm, normalization.method = "LogNormalize", scale.factor = 10000)
mm <- FindVariableFeatures(mm, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs

mm <- ScaleData(mm, features = rownames(mm))


xx = scCustomize::Convert_Assay(mm, assay = 'RNA', convert_to = 'V3')

save(xx, file = paste0(RdataDir, '/mouse_skin_macrophage_subtypes_SD_seuratv3.rds'))

##########################################
# check the pro-fibrotic macrophage activation 
##########################################
ax = readRDS(file = paste0(RdataDir, 
                           '/axoltol_limbBlatema_batch1_macrophage_time_subtypeAnnotations_rmEpidermis_v3_seuratv3.rds'))

jj = which(ax$cluster == 'M6.cyling')
ax$cluster[jj] = 'M2'
ax$subtypes[jj] = 'M1.APOE.PPARG.anti-inflammatory'

DimPlot(ax, group.by = 'cluster', label = TRUE, repel = TRUE)


genes = c(rownames(ax)[grep('CLEC10A|NINJ1|TREM2|FABP5|GPNMB|IGF1-|MS4A7|GAS6', 
                            rownames(ax))], "SPP1-AMEX60DD043905")

FeaturePlot(ax, features = genes, ncol = 3) &
  scale_color_gradient(low = "grey", high = "brown")

ggsave(filename = paste0(resDir, '/axoltol_macrophage_proFibrotic_genes.pdf'), 
       width = 16, height = 8)


##########################################
# ## define the one-on-one ortholog between axololt and mice
##########################################
ax = readRDS(file = paste0(RdataDir, 
            '/axoltol_limbBlatema_batch1_macrophage_time_subtypeAnnotations_rmEpidermis_v3_seuratv3.rds'))

DimPlot(ax, group.by = 'cluster', label = TRUE, repel = TRUE)

jj = which(ax$cluster == 'M6.cyling')
ax$cluster[jj] = 'M2'
ax$subtypes[jj] = 'M1.APOE.PPARG.anti-inflammatory'



mm = readRDS(file = paste0(RdataDir, '/mouse_skin_macrophage_subtypes_SD.rds'))
mm$species = 'mm'
#mm$subtypes[which(mm$subtypes == 'cycling')] = 'cycling.mm'
mm$time = droplevels(mm$condition)
mm$subtypes = mm$stage.m
mm$cluster = mm$stage.m

an_orthologs = data.frame(ref = rownames(ax), query = rownames(ax))
rownames(an_orthologs) = an_orthologs$ref
an_orthologs$query = sapply(an_orthologs$query, 
                            function(x){firstup(unlist(strsplit(as.character(x), '-'))[1])})

jj = which(!is.na(match(an_orthologs$query, rownames(mm))))
an_orthologs = an_orthologs[jj, ]

## select unique genes 
#counts = table(an_orthologs$query)
#counts = table(an_orthologs$query)
#gg_uniq = names(counts)[which(counts == 1)]
#jj2 = which(!is.na(match(an_orthologs$query, gg_uniq)))

gg_uniq = unique(an_orthologs$query)
jj2 = match(gg_uniq, an_orthologs$query)

an_orthologs = an_orthologs[jj2, ]
#an_orthologs = an_orthologs[jj2, ]

ax = subset(ax, features = an_orthologs$ref)
aa = subset(mm, features = an_orthologs$query)

rm(mm)

counts = ax@assays$RNA@counts
metadata = ax@meta.data
counts = counts[match(an_orthologs$ref, rownames(ax)), ]
rownames(counts) = an_orthologs$query

new_ax <- CreateSeuratObject(counts=counts, assay = 'RNA', meta.data = metadata)
new_ax<- NormalizeData(new_ax, normalization.method = "LogNormalize", scale.factor = 10000)
new_ax <- FindVariableFeatures(new_ax, selection.method = "vst", nfeatures = 5000)
new_ax <- ScaleData(new_ax, features = rownames(new_ax))

new_ax[['umap']] = ax[['umap']] 

counts = aa@assays$SCT@counts
metadata = aa@meta.data
#counts = counts[match(an_orthologs$ref, rownames(ax)), ]
#rownames(counts) = an_orthologs$query

new_aa <- CreateSeuratObject(counts=counts, assay = 'RNA', meta.data = metadata)
new_aa<- NormalizeData(new_aa, normalization.method = "LogNormalize", scale.factor = 10000)
new_aa <- FindVariableFeatures(new_aa, selection.method = "vst", nfeatures = 5000)
new_aa <- ScaleData(new_aa, features = rownames(new_aa))

new_aa[["umap"]] <- aa[['umap']]

DimPlot(new_aa, group.by = 'cluster', label = TRUE, repel =  TRUE)

aa = merge(new_aa, y = new_ax, add.cell.ids = c("m", "ax"), project = "skinRepair")

#rm(list = c('counts', 'metadata', 'ax'))
#rm(new_ax)

saveRDS(new_aa, file = paste0(RdataDir, 'mm_scRNAseq_for_crossSpecies_v3.rds'))
saveRDS(new_ax, file = paste0(RdataDir, 'ax_scRNAseq_for_crossSpecies_v3.rds'))
saveRDS(aa, file = paste0(RdataDir, 'mm_ax_scRNAseq_for_crossSpecies_v3.rds'))

########################################################
########################################################
# Section I: # test Seurat integration 
# 
########################################################
########################################################
source(paste0(functionDir, 'functions_dataIntegration.R'))

aa = readRDS(file = paste0(RdataDir, 'mm_ax_scRNAseq_for_crossSpecies_v3.rds'))

aa = ScaleData(aa, features = rownames(aa))
#method = 'Harmony'

if(method == 'Seurat_RPCA'){
  ref.combined = IntegrateData_Seurat_RPCA(aa, 
                                           group.by = 'species', 
                                           nfeatures = 5000,
                                           #merge.order = matrix(c(-2, 1, -3, -1), ncol = 2),
                                           reference = 2,
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
  source(paste0(functionDir, 'functions_dataIntegration.R'))
  ref.combined = IntegrateData_Seurat_CCA(aa, 
                                          group.by = 'species', 
                                          nfeatures = 5000,
                                          reference = c(2),
                                          ndims = c(1:50),
                                          k.anchor = 10,
                                          k.weight = 100,
                                          #merge.order = matrix(c(-2, 1, -3, -1), ncol = 2),
                                          redo.normalization.scaling = TRUE,
                                          correct.all = FALSE)
  
  p1 = DimPlot(ref.combined, group.by = 'subtypes', label = TRUE, repel = TRUE, raster=FALSE) + 
    ggtitle('Seurat_CCA')
  p2 = DimPlot(ref.combined, group.by = 'species', label = TRUE, repel = TRUE) +
    ggtitle("Seurat_CCA")
  
  p1 + p2
  
  ggsave(filename = paste0(resDir, '/cross_species_mapping_Seurat_CCA.pdf'), 
         width = 16, height = 6)
  
  saveRDS(ref.combined, file = paste0(resDir, 'crossSpecies_mm_ax_scRNAseq_SeuratCCA.rds'))
  
  
  xx <- IntegrateLayers(object = xx, method = CCAIntegration, orig.reduction = "pca", 
                        new.reduction = "integrated.cca",
                          verbose = FALSE)
  
  
}

##########################################
# use CCA integration but use reference 
##########################################
mm = readRDS(file = paste0(RdataDir, 'mm_scRNAseq_for_crossSpecies_v3.rds'))
ax = readRDS(file = paste0(RdataDir, 'ax_scRNAseq_for_crossSpecies_v3.rds'))

ax$subtypes = ax$cluster

aa = merge(mm, y = ax, add.cell.ids = c("m", "ax"), project = "skinRepair")

source(paste0(functionDir, 'functions_dataIntegration.R'))
ref.combined = IntegrateData_Seurat_CCA(aa, 
                                        group.by = 'species', 
                                        nfeatures = 3000,
                                        reference = c(2),
                                        ndims = c(1:30),
                                        k.anchor = 10,
                                        k.weight = 100,
                                        #merge.order = matrix(c(-2, 1, -3, -1), ncol = 2),
                                        redo.normalization.scaling = TRUE,
                                        correct.all = FALSE)

p1 = DimPlot(ref.combined, group.by = 'subtypes', label = TRUE, repel = TRUE, raster=FALSE) + 
  ggtitle('Seurat_CCA')
p2 = DimPlot(ref.combined, group.by = 'species', label = TRUE, repel = TRUE) +
  ggtitle("Seurat_CCA")

p1 + p2

ggsave(filename = paste0(resDir, '/cross_species_mapping_Seurat_CCA.pdf'), 
       width = 16, height = 6)


##########################################
# Harmony integration
##########################################
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


##########################################
# label transferring with CCA 
##########################################
mapping_method = "seurat_query_ref_mapping"

ref = mm

data_version = 'mapping_mNT.noRA.RA.d2_d5_Marioni2019_selectedCelltypes_test2'

#features.common = intersect(rownames(aa), rownames(ref))
#aa = subset(aa, features = features.common)
#ref = subset(ref, features = features.common)

# aa$dataset = 'mNT'
# aa$stage = aa$condition
# aa$sequencing.batch = 'mNT'
# ref$dataset = 'ref'
# 
# aa$celltype = paste0('mNT_', aa$condition)
# 
# outDir = paste0(resDir,  mapping_method, '/', data_version, '/')
# system(paste0('mkdir -p ', outDir))

ElbowPlot(ref, ndims = 50, reduction = 'pca')
ref = RunUMAP(ref, reduction = "pca", dims = 1:30, n.neighbors = 30, 
              min.dist = 0.1, return.model = TRUE) 


DimPlot(ref, reduction = "umap", 
        group.by = "subtypes", label = TRUE,
        repel = TRUE, raster=FALSE) 

ref$labels = factor(ref$subtypes)

# In data transfer, Seurat has an option (set by default) to project the PCA structure of a reference 
# onto the query, instead of learning a joint structure with CCA. 
# We generally suggest using this option when projecting data between scRNA-seq datasets.
anchors <- FindTransferAnchors(reference = ref, 
                               query = ax, 
                               dims = 1:50,
                               normalization.method = "LogNormalize",
                               #reference.reduction = "pca.corrected",
                               max.features = 200,
                               k.anchor = 10,
                               reduction = "cca" 
)

predictions <- TransferData(anchorset = anchors, refdata = ref$labels, 
                            reference = ref,
                            query = ax,
                            weight.reduction = 'cca',
                            dims = 1:50)

saveRDS(anchors, file = paste0(RdataDir, '/cca_anchors_transferData.rds'))
query = ax;
query$predicted.id = predictions$predicted.id 
query$predicted.score = predictions$predicted.id.score

saveRDS(query, file = paste0(RdataDir, '/axolotl_macrophage_cca_anchors_transferData.rds'))

jj = which(query$predicted.score < 0.5)
query$predicted.id[jj] = NA


p1 = DimPlot(query, reduction = "umap", 
             group.by = "predicted.id", label = TRUE,
             repel = TRUE, raster=FALSE) 

p2 = FeaturePlot(query, features = 'predicted.score')

p1 + p2

ggsave(paste0(resDir, '/transfer_learning_seurat_cca.ref_v3.pdf'), 
       width = 16, height = 8)





## visualize the query cells alongside our reference and didn't work well
## (MapQuery in Seurat didn't work well)
# query <- Seurat::MapQuery(anchorset = anchors, 
#                   reference = ref, 
#                   query = aa,
#                   refdata = list(labels = "labels"),
#                   #refdata = list(celltype = "celltype"), 
#                   transferdata.args = list(weight.reduction = 'rpca.ref'), 
#                   reference.reduction = "pca", 
#                   reduction.model = "umap")
# 
# p1 <- DimPlot(ref, reduction = "umap", group.by = "celltype", 
#               label = TRUE, label.size = 3,
#               repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
# p2 <- DimPlot(query, reduction = "ref.umap", group.by = "predicted.labels", 
#               label = TRUE,
#               label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
# p1 + p2


##########################################
# label transferring with   
##########################################
Test_reference_mapping_Symphony = FALSE
if(Test_reference_mapping_Symphony){
  library(symphony)
  library(singlecellmethods)
  source('/groups/tanaka/People/current/jiwang/projects/RA_competence/scripts/utils_symphony.R')
  
  fig.size <- function (height, width) {
    options(repr.plot.height = height, repr.plot.width = width)
  }
  
  # mapping_method = 'symphony_mapping'
  # data_version = "mapping_mNT.noRA.RA.d2_d5_Marioni2019_selectedCelltypes"
  # 
  # outDir = paste0(resDir,  mapping_method, '/', data_version, '/')
  # system(paste0('mkdir -p ', outDir))
  
  ## import and prepare the ref and query
  #aa = readRDS(file = paste0(RdataDir, 
  #                           'seuratObject_mNT_selectedCondition_downsampled.1k.perCondition_reclustered.rds'))
  #ref = readRDS(file = paste0(RdataDir,  
  #                            'seuratObject_EmbryoAtlasData_all36sample_RNAassay_keep.relevant.celltypes_v3.rds'))
  
  # features.common = intersect(rownames(aa), rownames(ref))
  # aa = subset(aa, features = features.common)
  # ref = subset(ref, features = features.common)
  # 
  # aa$dataset = 'mNT'
  # aa$stage = aa$condition
  # aa$sequencing.batch = 'mNT'
  # ref$dataset = 'ref'
  # aa$celltype = paste0('mNT_', aa$condition)
  # 
  # ref$labels = paste0(ref$celltype, '_', ref$stage)
  
  
  #idx_query = which(metadata$donor == "5'") # use 5' dataset as the query
  ref_exp_full = ref@assays$RNA@data
  ref_metadata = ref@meta.data
  query_exp = ax@assays$RNA@data
  query_metadata = ax@meta.data
  
  # Sparse matrix with the normalized genes x cells matrix
  ref_exp_full[1:5, 1:2]
  
  # Select variable genes and subset reference expression by variable genes
  var_genes = vargenes_vst(ref_exp_full, groups = as.character(ref_metadata[['condition']]), topn = 500)
  ref_exp = ref_exp_full[var_genes, ]
  dim(ref_exp)
  
  # Build reference
  reference = symphony::buildReference(
    ref_exp,                   # reference expression (genes by cells)
    ref_metadata,              # reference metadata (cells x attributes)
    vars = NULL,         # variable(s) to integrate over
    K = 100,                   # number of Harmony soft clusters
    verbose = TRUE,            # display verbose output
    do_umap = TRUE,            # run UMAP and save UMAP model to file
    do_normalize = FALSE,      # perform log(CP10k) normalization on reference expression
    vargenes_method = 'vst',   # variable gene selection method: 'vst' or 'mvp'
    vargenes_groups = 'condition', # metadata column specifying groups for variable gene selection within each group
    topn = 500,               # number of variable genes (per group)
    theta = 2,                 # Harmony parameter(s) for diversity term
    d = 20,                    # number of dimensions for PCA
    save_uwot_path = '/groups/tanaka/People/current/jiwang/projects/fibroticRepair_axlotl_mouse/model_test/', # file path to save uwot UMAP model
    additional_genes = NULL    # vector of any additional genes to force include
  )
  
  #Run Harmony integration
  Run_Harmony_integration = FALSE
  if(Run_Harmony_integration){
    # Calculate and save the mean and standard deviations for each gene
    vargenes_means_sds = tibble(symbol = var_genes, mean = Matrix::rowMeans(ref_exp))
    vargenes_means_sds$stddev = singlecellmethods::rowSDs(ref_exp,  row_means = vargenes_means_sds$mean)
    head(vargenes_means_sds)
    
    #Scale data using calculated gene means and standard deviations
    ref_exp_scaled = singlecellmethods::scaleDataWithStats(ref_exp, vargenes_means_sds$mean, 
                                                           vargenes_means_sds$stddev, 1)
    
    #Run SVD, save gene loadings (s$u)
    set.seed(0)
    s = irlba::irlba(ref_exp_scaled, nv = 20)
    Z_pca_ref = diag(s$d) %*% t(s$v) # [pcs by cells]
    loadings = s$u
    
    set.seed(0)
    ref_harmObj = harmony::HarmonyMatrix(
      data_mat = t(Z_pca_ref),  ## PCA embedding matrix of cells
      meta_data = ref_metadata, ## dataframe with cell labels
      theta = c(2),             ## cluster diversity enforcement
      vars_use = NULL,    ## variable to integrate out
      nclust = NULL,             ## number of clusters in Harmony model
      max.iter.harmony = 20,
      return_object = TRUE,     ## return the full Harmony model object
      do_pca = FALSE            ## don't recompute PCs
    )
    
    # To run the next function buildReferenceFromHarmonyObj(), 
    # you need to input the saved gene loadings (loadings) and vargenes_means_sds.
    # Compress a Harmony object into a Symphony reference
    reference = symphony::buildReferenceFromHarmonyObj(
      ref_harmObj,            # output object from HarmonyMatrix()
      ref_metadata,           # reference cell metadata
      vargenes_means_sds,     # gene names, means, and std devs for scaling
      loadings,               # genes x PCs matrix
      verbose = TRUE,         # verbose output
      do_umap = TRUE,         # Set to TRUE only when UMAP model was saved for reference
      save_uwot_path = '/groups/tanaka/People/current/jiwang/projects/RA_competence/results/scRNAseq_R13547_10x_mNT_20220813/mapping_to_MouseGastrulationData/symphony_mapping/mapping_mNT.noRA.RA.d2_d5_Marioni2019_selectedCelltypes/model_test')
    
  }
  
  # Optionally, you can specify which normalization method was
  # used to build the reference as a custom slot inside the Symphony object to 
  # help record this information for future query users
  #reference$normalization_method = 'log(CP10k+1)'
  saveRDS(reference, paste0(outDir, 'testing_reference_mouseGastrulation_1.rds'))
  
  str(reference)
  
  # The harmonized embedding is located in the Z_corr slot of the reference object.
  dim(reference$Z_corr)
  reference$Z_corr[1:5, 1:5]
  
  # reference = readRDS(paste0(outDir, 'testing_reference1.rds'))
  umap_labels = cbind(ref_metadata, reference$umap$embedding)
  
  fig.size(3, 5)
  plotBasic(umap_labels, title = 'Reference', color.by = 'subtypes')
  
  # In order to map a new query dataset onto the reference, 
  # you will need a reference object saved from the steps above, 
  # as well as query cell expression and metadata.
  # Map query
  query = mapQuery(query_exp,             # query gene expression (genes x cells)
                   query_metadata,        # query metadata (cells x attributes)
                   reference,             # Symphony reference object
                   do_normalize = FALSE,  # perform log(CP10k+1) normalization on query
                   do_umap = TRUE)        # project query cells into reference UMAP
  
  ## Symphony assumes that the query is normalized in the same manner as the reference. 
  # Our implementation currently uses log(CP10k) normalization.
  
  str(query)
  
  #Let's take a look at what the query object contains:
  
  # Z: query cells in reference Harmonized embedding
  #Zq_pca: query cells in pre-Harmony reference PC embedding (prior to correction)
  # R: query cell soft cluster assignments
  #  Xq: query cell design matrix for correction step
  #  umap: query cells projected into reference UMAP coordinates (using uwot)
  #  meta_data: metadata
  # Predict query cell types using k-NN
  query = knnPredict(query, reference, reference$meta_data$subtypes, k = 5)
  
  ## Query cell type predictions are now in the cell_type_pred_knn column. 
  ## The cell_type_pred_knn_prob column reports the proportion of nearest neighbors with the winning vote 
  # (can help identify query cells that fall "on the border" between 2 reference cell types).
  
  head(query$meta_data)
  
  # Add the UMAP coordinates to the metadata
  reference$meta_data$cell_type_pred_knn = NA
  reference$meta_data$cell_type_pred_knn_prob = NA
  reference$meta_data$ref_query = 'reference'
  query$meta_data$ref_query = 'query'
  
  # Add the UMAP coordinates to the metadata
  umap_combined = rbind(query$umap, reference$umap$embedding)
  xx = query$meta_data[, c(24, 22)]
  colnames(xx)[2] = 'celltype' 
  meta_data_combined = rbind(xx, reference$meta_data[, c(27, 18)])
  
  umap_combined_labels = cbind(meta_data_combined, umap_combined)
  
  fig.size(6, 14)
  plotBasic(umap_combined_labels, title = 'Reference and query cells', 
            color.by = 'celltype', facet.by = 'ref_query')
  
  ggsave(paste0(outDir, '/mapping_symphony_ref_v1_celltype_coembedding.pdf'), 
         width = 16, height = 8)
  
  ax$predicted.id = query$meta_data$cell_type_pred_knn
  ax$predicted.score = query$meta_data$cell_type_pred_knn_prob
  
  p2 = FeaturePlot(ax, features = 'predicted.score')
  

  p1 = DimPlot(ax, reduction = "umap", 
               group.by = "predicted.id", label = TRUE,
               repel = TRUE, raster=FALSE) 
  
  p1 + p2
  
  
  jj = which(ax$predicted.score < 0.5)
  ax$predicted.id[jj] = NA
  
  p1 = DimPlot(ax, reduction = "umap", 
               group.by = "predicted.id", label = TRUE,
               repel = TRUE, raster=FALSE) 
  
  p1 + p2
  
  
  
  plot(p1)
  ggsave(paste0(outDir, '/mapping_symphony_ref_celltype_v1.pdf'), 
         width = 16, height = 8)
  
  
}


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

saveFile = paste0(resDir, 'mm_ax_scRNAseq_merged_v2_seuratV4.h5Seurat')
SaveH5Seurat(aa, filename = saveFile, overwrite = TRUE)
Convert(saveFile, dest = "h5ad", overwrite = TRUE)



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
mm = readRDS(file = paste0(RdataDir, 'mm_scRNAseq_for_crossSpecies_v3.rds'))
ax = readRDS(file = paste0(RdataDir, 'ax_scRNAseq_for_crossSpecies_v3.rds'))
mm <- FindVariableFeatures(mm, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs
mm <- ScaleData(mm)
mm <- RunPCA(mm, verbose = FALSE, weight.by.var = FALSE)

ax <- FindVariableFeatures(ax, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs
ax <- ScaleData(ax)
ax <- RunPCA(ax, verbose = FALSE, weight.by.var = FALSE)
#ax = readRDS(file = paste0(RdataDir, 'axoltol_limbBlatema_batch1_macrophage_subtypes.rds'))
#mm = readRDS(file = paste0(RdataDir, 'mouse_skin_macrophage_subtypes.rds'))
#mm = subset(mm, cells = colnames(mm)[which(!is.na(mm$subtypes))])

ax$subtypes = ax$cluster

#ax$species = 'ax'
#ax$time = droplevels(ax$time)
#ax$subtypes[which(ax$subtypes == 'cycling')] = 'cycling.ax'

ax = as.SingleCellExperiment(ax, assay = 'RNA')

#mm$species = 'mm'
#mm$subtypes[which(mm$subtypes == 'cycling')] = 'cycling.mm'
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
  
  a_milo <- buildGraph(a_milo, k=20, d=30, reduced.dim="PCA")
  
  a_milo <- makeNhoods(a_milo, prop=0.1, k=21, d=30, refined=T, reduced_dims="PCA")
  
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
head(a_milo$subtypes)

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
                               max_hvgs=5000, 
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

ggsave(paste0(resDir, "axolotl_nm_mapping_test_v11.pdf"), width=18, height=10, dpi=300)


########################################################
########################################################
# Section IV: Cross-species cell type correlations 
# the same idea from https://www.science.org/doi/10.1126/science.abp9262
# Tomas and Ashley's paper
########################################################
########################################################
mm = readRDS(file = paste0(RdataDir, 'mm_scRNAseq_for_crossSpecies_v3.rds'))
ax = readRDS(file = paste0(RdataDir, 'ax_scRNAseq_for_crossSpecies_v3.rds'))
#saveRDS(aa, file = paste0(RdataDir, 'mm_ax_scRNAseq_for_crossSpecies_v3.rds'))

p1 = DimPlot(mm, group.by = 'subtypes', label = TRUE, repel =  TRUE)
p2 = DimPlot(ax, group.by = 'cluster', label = TRUE, repel =  TRUE)

p1 + p2

ggsave(filename = paste0(resDir, '/cross_species_macrophageData_mm_ax_v3.pdf'), 
       width = 16, height = 6)

cols_cluster = c( "#7F7F7F", "#FFC000", "#70AD47", "#337f01", "#265401", "#CD00CF", "#800080", 
                  "#EFFAB6", "#69C6BE", "#007BB7", "#121D60",
                  "#69C6BE", "#007BB7", "#121D60")

cols_macrophage = c("#EFFAB6", "#70AD47", "#69C6BE", "#007BB7", 'royalblue')

DimPlot(ax, group.by = 'cluster', label = TRUE, repel =  TRUE) + 
  #theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 14, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 14)) +
  scale_color_manual(values=cols_macrophage)

ggsave(filename = paste0(resDir, '/axolotl_macrophage_subtypes_v3.pdf'), 
       width = 6, height = 4)








##########################################
# find all DE genes for axolotl and mouse
##########################################
Idents(ax) = ax$cluster
Idents(mm) = mm$subtypes
markers.ax = FindAllMarkers(ax, logfc.threshold = 0.25, min.pct = 0.1)

cat(length(unique(markers.ax$gene)), ' markers found in ax\n')
#markers.mm = FindAllMarkers(mm, logfc.threshold = 0.5, min.pct = 0.1)
#saveRDS(markers.mm, file = paste0(RdataDir, 'mm_scRNAseq_for_crossSpecies_v3_markersAll.rds'))
markers.mm = readRDS(file = paste0(RdataDir, 'mm_scRNAseq_for_crossSpecies_v3_markersAll.rds'))

cat(length(unique(markers.mm$gene)), ' markers found in mm\n')


## define the one-on-one ortholog between axololt and mice
an_orthologs = data.frame(ref = rownames(ax), query = rownames(ax))
rownames(an_orthologs) = an_orthologs$ref
an_orthologs$query = sapply(an_orthologs$query, 
                            function(x){firstup(unlist(strsplit(as.character(x), '-'))[1])})

## intersect gene names with mouse data
jj = which(!is.na(match(an_orthologs$query, rownames(mm))))
an_orthologs = an_orthologs[jj, ]

## select unique genes 
gg_uniq = unique(an_orthologs$query)
jj2 = match(gg_uniq, an_orthologs$query)

#counts = table(an_orthologs$query)
#gg_uniq = names(counts)[which(counts == 1)]
#jj2 = which(!is.na(match(an_orthologs$query, gg_uniq)))

an_orthologs = an_orthologs[jj2, ]

rownames(an_orthologs) = an_orthologs$query
cat(nrow(an_orthologs), 'genes used for orthologs\n')

## intersect with DE genes of axolotl and mouse
mm1 = match(an_orthologs$ref, markers.ax$gene)
mm2 = match(an_orthologs$query, markers.mm$gene)
kk = intersect(which(!is.na(mm1)), which(!is.na(mm2)))

an_orthologs = an_orthologs[kk, ]

an_orthologs$tfs = FALSE
an_orthologs$tfs[which(!is.na(match(toupper(an_orthologs$query), c(tfs, sps))))] = TRUE

avg.ax = AggregateExpression(ax, features = an_orthologs$ref[which(an_orthologs$tfs == TRUE)])
avg.ax = data.frame(avg.ax$RNA)

cat(nrow(avg.ax), 'TFs used for correlation analysis\n')

avg.mm = AggregateExpression(mm, features = an_orthologs$query[which(an_orthologs$tfs == TRUE)])
avg.mm = data.frame(avg.mm$RNA)

cat(nrow(avg.mm), 'TFs used for correlation analysis\n')

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
  geom_point(data = plot_df, mapping = aes(x = Var2, 
                                           y = factor(Var1, levels = c('M3', 'M5', 'M1', 'M2', 'M4')), 
                                                                fill = value), 
             shape = 21, size = 10) +
  #geom_point(data = plot_df[plot_df$rowmax,], mapping = aes(x = Var2, y = Var1, size = 30), 
  #           shape = "—", show.legend = F, colour = "grey10")+
  #geom_point(data = plot_df[plot_df$colmax,], mapping = aes(x = Var2, y = Var1, size = 30), 
  #           shape = "|", show.legend = F, colour = "grey10")+
  scale_x_discrete(expand = c(0,0.7)) +
  scale_y_discrete(expand = c(0,0.7)) +
  scale_fill_gradientn(breaks = signif(c(min(cort$r)+0.005, 0, max(cort$r)-0.005),2), 
                       #values = scales::rescale(c(min(br), 0, max(br))),
                       colours = cols) +
  labs(x = 'mouse', y = 'axolotl', fill = "Spearman's\nrho", size = "-log10\nadj. p-value")+
  theme_classic()+
  theme(axis.title = element_text(colour = "black", face = "bold", size = 14),
        axis.text = element_text(colour = "black", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1.0, vjust = 1.0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))

ggsave(paste0(resDir, "axolotl_nm_subtype_correlationAnalysis_v2.pdf"), width=8, height=6, dpi=300)


