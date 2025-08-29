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

species = 'mouse'
version.analysis = '_20250829'
resDir = paste0("../results/scRNAseq_analysis_immune", version.analysis, '/')
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts'
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_scRNAseq.R')
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_Visium.R')

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

design = design[-1, ]

# design = design[which(design$sampleID != '196323'), ]

# import data from cellranger output
for(n in 1:nrow(design))
{
  # n = 1
  cat(n, ' : ', design$condition[n], '\n')
  
  topdir = paste0(dataDir, design$sampleID[n], '/outs/filtered_feature_bc_matrix/')
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
  #iscell_dd = defaultDrops(count.data, expected = 8000) # default cell estimate, similar to 10x cellranger
  #sum(iscell_dd, na.rm=TRUE)
  
  ## not used the emptyDrops too slow 
  # eout = emptyDrops(count.data, lower = 200)
  # eout$FDR[is.na(eout$FDR)] = 1
  # iscell_ed = eout$FDR<=0.01
  # sum(iscell_ed, na.rm=TRUE)
  
  meta = data.frame(row.names = colnames(count.data), condition = rep(design$condition[n], ncol(count.data)))
  
  # plot rankings for number of UMI
  # br.out <- barcodeRanks(count.data)
  # 
  # pdf(paste0(resDir, "/UMIrank_emptyDrop_", design$condition[n], "_", design$sampleID[n],  ".pdf"), 
  #     height = 6, width =10, useDingbats = FALSE)
  # 
  # plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
  # 
  # o <- order(br.out$rank)
  # lines(br.out$rank[o], br.out$fitted[o], col="red")
  # abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
  # abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
  # abline(v = sum(iscell_dd), col = 'darkgreen', lwd = 2.0)
  # abline(v = c(3000, 5000, 8000, 10000, 12000), col = 'gray')
  # text(x = c(3000, 5000, 8000, 10000, 12000), y =10000, labels = c(3000, 5000, 8000, 10000, 12000), 
  #      col = 'red')
  # legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
  #        legend=c("knee", "inflection"))
  # 
  # dev.off()
  
  # use defaultDrop to select cells.
  aa = CreateSeuratObject(counts = count.data,
                          meta.data = meta, 
                          min.cells = 20, min.features = 100)
  aa$cell.id = paste0(colnames(aa), '_', design$condition[n], '_', design$sampleID[n])
  
  if(n == 1) {
    mnt = aa
  }else{
    mnt = merge(mnt, aa)
  }
  
}

mnt[["percent.mt"]] <- PercentageFeatureSet(mnt, pattern = "^mt-")

saveRDS(mnt, 
     file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, version.analysis, '.rds'))


##########################################
# QCs and cell filtering  
##########################################
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, version.analysis, '.Rdata'))
aa = mnt
rm(mnt)

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

#ggs = sapply(rownames(aa), function(x) {x = unlist(strsplit(x, '-')); x = x[grep('ENSMUSG', x, invert = TRUE)];
#                                                                               paste0(x, collapse = '-')})

features = c('Pou5f1', 'Sox2', 'Lef1', 'Otx2', 'Zfp703', 'Pax6', 'Foxa2', 'Shh', 'Nkx6-1', 'Nkx2-2', 'Olig2', 
             'Sox1', 'Tubb3')
features = rownames(aa)[!is.na(match(rownames(aa), features))]

Idents(aa) = factor(aa$condition, levels = levels)

for(n in 1:length(features))
{
  p1 = VlnPlot(aa, features = features[n])
  plot(p1)
}

dev.off()

rm(count.data); rm(br.out); rm(bc);rm(meta)

##########################################
## filter cells here
##########################################
VlnPlot(aa, features = 'nFeature_RNA', y.max = 10000) +
  geom_hline(yintercept = c(1000, 8000), col = 'red')

VlnPlot(aa, features = 'percent.mt', y.max =20) + 
  geom_hline(yintercept = c(5), col = 'red')

aa <- subset(aa, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 5)

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_', species, version.analysis, '.rds'))



