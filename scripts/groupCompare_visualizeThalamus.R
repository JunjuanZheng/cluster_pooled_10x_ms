# Kristin Muench
# 2019.03.28
# Pipeline to run Seurat setup for mouse scRNA-Seq data - with CCA
# This takes in a matrix with CCAs calculated, adjusts them, makes a tSNE, 
# and visualizes/counts clusters.

# Load older version of Seurat which has the functions I need
library('Seurat', lib.loc='/scg/apps/software/r/3.5.0/scg/seurat_2.3')
packageVersion('Seurat')

## run this in case you accidentally load Seurat 3
#detach("package:Seurat", unload=TRUE)

# # Set up workspace - uncomment for running job
# print('Run CCA Track: groupCompare.R')
# print('Setting up workspace...')
# ##remove all environmental vars for clarity
# rm(list=ls())
# 
# ## import from command line
# args <- commandArgs(TRUE)
# inputFile <- args[1]
# outputDir <- args[2]
# chosen_cc <- args[3]
# chosen_res <- args[4]
# print(paste0('Importing files from: ', inputFile))
# print(paste0('Output location: ', outputDir))
# print(paste0('Chosen # of CCs: ', chosen_cc))
# print(paste0('Chosen resolution: ', chosen_res))

# Set up workspace -  Uncomment to work with files
inputFile <- "/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190329_troubleshootRunthrough_demux/clusters/data.combined_withClust_r1.2_CC40.RData"
outputDir <- "/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190329_troubleshootRunthrough_demux"
chosen_cc <- 40
chosen_res <- 1.2

## load needed files
print('~*~')
print('Loading needed files...')
load(inputFile)

## make subdirectory
setwd(outputDir)
subDir <- 'groupCompare'
print(paste0('Name of Experiment: ', subDir))

if (file.exists(subDir)){
  setwd(file.path(outputDir, subDir))
} else {
  dir.create(file.path(outputDir, subDir))
  setwd(file.path(outputDir, subDir))
}

print('Names of metadata')
print(names(data.combined@meta.data))


# Declare Stanford colors for prettier visualization
stanColors <- data.frame(cond = c('WT.SAL', 'WT.LPS', 'HET.SAL', 'HET.LPS'),
                         colors = as.character(c('#009b76', '#8c1515', '#b26f16', '#53284f')))

# # DE Analysis
# 
# ## Make container file to store output
# setwd(paste0(outputDir, subDir))
# subSubDir <- 'differentialExpression'
# print(subSubDir)
# print(file.path(outputDir, subDir, subSubDir))
# 
# if (file.exists(subSubDir)){
#   setwd(file.path(outputDir, subDir, subSubDir))
# } else {
#   dir.create(file.path(outputDir, subDir, subSubDir))
#   setwd(file.path(outputDir, subDir, subSubDir))
# }
# 
# ## create new group in metadata for both cell type AND cond info
# data.combined@meta.data$celltype.cond <- paste0(data.combined@ident, "_",
#                                                 data.combined@meta.data$cond)
# 
# ## stash those labels for later use
# data.combined <- StashIdent(data.combined, save.name = "celltype")
# data.combined <- SetAllIdent(data.combined, id = "celltype.cond") # tell Seurat we are focusing on this grouping method now
# 
# ## DE between two groups
# deGroups <- function(data.combined, xAxis, yAxis){
#   # create title
#   myTitle <- paste0(xAxis, '_vs_', yAxis)
#   print(paste0('DE genes for ', myTitle, '...'))
#   
#   # do DE
#   ident1_xAxis_yAxis <- FindMarkers(data.combined, ident.1 = xAxis, ident.2 = yAxis, print.bar = TRUE)
#   head(ident1_xAxis_yAxis, 15)
#   
#   ident1_xAxis_yAxis
#   
#   # save results
#   write.csv(ident1_xAxis_yAxis, file = paste0(myTitle, '_r', chosen_res, '_cc', chosen_cc, '.csv'), quote=FALSE)
#   
#   # return DE lists
#   return(ident1_xAxis_yAxis)
# }
# 
# ## effect of LPS
# g0_WT.SAL_WT.LPS <- deGroups(data.combined, "0_WT.SAL", "0_WT.LPS")
# g1_WT.SAL_WT.LPS <- deGroups(data.combined, "1_WT.SAL", "1_WT.LPS")
# g2_WT.SAL_WT.LPS <- deGroups(data.combined, "2_WT.SAL", "2_WT.LPS")
# g3_WT.SAL_WT.LPS <- deGroups(data.combined, "3_WT.SAL", "3_WT.LPS")
# g4_WT.SAL_WT.LPS <- deGroups(data.combined, "4_WT.SAL", "4_WT.LPS")
# ## on microglia
# g30_WT.SAL_WT.LPS <- deGroups(data.combined, "30_WT.SAL", "30_WT.LPS")
# 
# 
# ## effect of genotype
# g0_WT.SAL_HET.SAL <- deGroups(data.combined, "0_WT.SAL", "0_HET.SAL")
# g1_WT.SAL_HET.SAL <- deGroups(data.combined, "1_WT.SAL", "1_HET.SAL")
# g4_WT.SAL_HET.SAL <- deGroups(data.combined, "4_WT.SAL", "4_HET.SAL")
# g7_WT.SAL_HET.SAL <- deGroups(data.combined, "7_WT.SAL", "7_HET.SAL")
# ### on neurons
# g9_WT.SAL_HET.SAL <- deGroups(data.combined, "9_WT.SAL", "9_HET.SAL")
# g14_WT.SAL_HET.SAL <- deGroups(data.combined, "14_WT.SAL", "14_HET.SAL")
# ### on interneurons
# g13_WT.SAL_HET.SAL <- deGroups(data.combined, "13_WT.SAL", "13_HET.SAL")
# g17_WT.SAL_HET.SAL <- deGroups(data.combined, "17_WT.SAL", "17_HET.SAL")
# g20_WT.SAL_HET.SAL <- deGroups(data.combined, "20_WT.SAL", "20_HET.SAL")
# g21_WT.SAL_HET.SAL <- deGroups(data.combined, "21_WT.SAL", "21_HET.SAL")
# g22_WT.SAL_HET.SAL <- deGroups(data.combined, "22_WT.SAL", "22_HET.SAL")
# g23_WT.SAL_HET.SAL <- deGroups(data.combined, "23_WT.SAL", "23_HET.SAL")
# g25_WT.SAL_HET.SAL <- deGroups(data.combined, "25_WT.SAL", "25_HET.SAL")
# ## on microglia
# g30_WT.SAL_HET.SAL <- deGroups(data.combined, "30_WT.SAL", "30_HET.SAL")
# 
# ## effect of both
# g13_WT.SAL_HET.LPS <- deGroups(data.combined, "13_WT.SAL", "13_HET.LPS")
# g23_WT.SAL_HET.LPS <- deGroups(data.combined, "23_WT.SAL", "23_HET.LPS")
# 
# 
# ## comparison between groups
# g13_WT.SAL_WT.SAL <- deGroups(data.combined, "13_WT.SAL", "23_WT.SAL")
# 
# ## !!!! Make sure this isnt the effect of one outlier that's high in Hmgb1 or something
# 
# 
# ## Set back identity
# data.combined <- SetAllIdent(data.combined, id = "res.1.2") # tell Seurat we are focusing on the cluster names now
# 
# 
# 
# # Differential Expression Analyses
# 
# ## When comparing two samples, what genes have a large difference?
# 
# ## Make container file to store output
# setwd(paste0(outputDir, subDir))
# subSubDir <- 'plotDifferencesInMeans'
# print(subSubDir)
# print(file.path(outputDir, subDir, subSubDir))
# 
# if (file.exists(subSubDir)){
#   setwd(file.path(outputDir, subDir, subSubDir))
# } else {
#   dir.create(file.path(outputDir, subDir, subSubDir))
#   setwd(file.path(outputDir, subDir, subSubDir))
# }
# 
# ## Function
# plotOneGroup <- function(data.combined, ident1, xAxis, yAxis, genesToLabel){
#   # name for this run
#   myTitle <- paste0(xAxis, '_vs_', yAxis, '_grp', ident1)
#   print(paste0('Plotting genes for ', myTitle, '...'))
#   
#   # get plotting data for first
#   cells1 <- SubsetData(data.combined, ident.use = ident1, subset.raw = T)
#   cells1 <- SetAllIdent(cells1, id = "cond")
#   avg.cells1 <- log1p(AverageExpression(cells1, show.progress = FALSE))
#   avg.cells1$gene <- rownames(avg.cells1)
#   
#   # clarify what genes to plot
#   avg.cells1$plotGene <- "" # effectively hiding all labels
#   avg.cells1[genesToLabel, 'plotGene'] <- rownames(avg.cells1[genesToLabel,])
#   
#   # plot
#   library(ggrepel)
#   p1 <- ggplot(avg.cells1, aes_string(x = xAxis, y = yAxis, label = 'plotGene')) + geom_point() + ggtitle(ident1) + geom_text_repel()
#   
#   ggsave(file = paste0(myTitle, '_r', chosen_res, '_cc', chosen_cc, '.pdf'), plot = p1, device='pdf')
#   
#   # what genes have a large difference in expression?
#   bigDiff <- data.frame(abs(avg.cells1[, xAxis] - avg.cells1[, yAxis]))
#   row.names(bigDiff) <- row.names(avg.cells1)
#   # write data to file
#   write.csv(bigDiff, file = paste0(myTitle, '_r', chosen_res, '_cc', chosen_cc, '.csv'), quote=FALSE)
#   
#   return(bigDiff)
#   
# }
# 
# ## Effect of LPS
# bigDiff_0_wt.sal_wt.lps <- plotOneGroup(data.combined, "0", 'WT.SAL', 'WT.LPS', row.names(g0_WT.SAL_WT.LPS[g0_WT.SAL_WT.LPS$p_val_adj < 0.05,]) )
# bigDiff_1_wt.sal_wt.lps <- plotOneGroup(data.combined, "1", 'WT.SAL', 'WT.LPS', row.names(g1_WT.SAL_WT.LPS[g1_WT.SAL_WT.LPS$p_val_adj < 0.05,]) )
# bigDiff_2_wt.sal_wt.lps <- plotOneGroup(data.combined, "2", 'WT.SAL', 'WT.LPS', row.names(g2_WT.SAL_WT.LPS[g2_WT.SAL_WT.LPS$p_val_adj < 0.05,]) )
# bigDiff_3_wt.sal_wt.lps <- plotOneGroup(data.combined, "3", 'WT.SAL', 'WT.LPS', row.names(g3_WT.SAL_WT.LPS[g3_WT.SAL_WT.LPS$p_val_adj < 0.05,] ) )
# bigDiff_4_wt.sal_wt.lps <- plotOneGroup(data.combined, "4", 'WT.SAL', 'WT.LPS', row.names(g4_WT.SAL_WT.LPS[g4_WT.SAL_WT.LPS$p_val_adj < 0.05,] ) )
# ## on microglia
# #bigDiff_30_wt.sal_wt.lps <- plotOneGroup(data.combined, "30", 'WT.SAL', 'WT.LPS', row.names(g30_WT.SAL_WT.LPS[g30_WT.SAL_WT.LPS$p_val_adj < 0.05,] ) )
# bigDiff_30_wt.sal_wt.lps <- plotOneGroup(data.combined, "30", 'WT.SAL', 'WT.LPS', row.names(g30_WT.SAL_WT.LPS) )
# 
# 
# # effect of genotype
# ## progenitors
# bigDiff_0_wt.sal_het.sal <- plotOneGroup(data.combined, "0", 'WT.SAL', 'HET.SAL', row.names(g0_WT.SAL_HET.SAL[g0_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
# bigDiff_1_wt.sal_het.sal <- plotOneGroup(data.combined, "1", 'WT.SAL', 'HET.SAL', row.names(g1_WT.SAL_HET.SAL[g1_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
# bigDiff_4_wt.sal_het.sal <- plotOneGroup(data.combined, "4", 'WT.SAL', 'HET.SAL', row.names(g4_WT.SAL_HET.SAL[g4_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
# bigDiff_7_wt.sal_het.sal <- plotOneGroup(data.combined, "7", 'WT.SAL', 'HET.SAL', row.names(g7_WT.SAL_HET.SAL[g7_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
# ## on neurons
# bigDiff_9_wt.sal_het.sal <- plotOneGroup(data.combined, "9", 'WT.SAL', 'HET.SAL', row.names(g9_WT.SAL_HET.SAL[g9_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
# bigDiff_14_wt.sal_het.sal <- plotOneGroup(data.combined, "14", 'WT.SAL', 'HET.SAL', row.names(g14_WT.SAL_HET.SAL[g14_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
# ## on interneurons - not all of these signif
# bigDiff_13_wt.sal_het.sal <- plotOneGroup(data.combined, "13", 'WT.SAL', 'HET.SAL', row.names(g13_WT.SAL_HET.SAL[g13_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
# bigDiff_17_wt.sal_het.sal <- plotOneGroup(data.combined, "17", 'WT.SAL', 'HET.SAL', row.names(g17_WT.SAL_HET.SAL[g17_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
# bigDiff_20_wt.sal_het.sal <- plotOneGroup(data.combined, "20", 'WT.SAL', 'HET.SAL', row.names(g20_WT.SAL_HET.SAL[g20_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
# bigDiff_21_wt.sal_het.sal <- plotOneGroup(data.combined, "21", 'WT.SAL', 'HET.SAL', row.names(g21_WT.SAL_HET.SAL[g21_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
# bigDiff_22_wt.sal_het.sal <- plotOneGroup(data.combined, "22", 'WT.SAL', 'HET.SAL', row.names(g22_WT.SAL_HET.SAL[g22_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
# bigDiff_23_wt.sal_het.sal <- plotOneGroup(data.combined, "23", 'WT.SAL', 'HET.SAL', row.names(g23_WT.SAL_HET.SAL[g23_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
# bigDiff_25_wt.sal_het.sal <- plotOneGroup(data.combined, "25", 'WT.SAL', 'HET.SAL', row.names(g25_WT.SAL_HET.SAL[g25_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
# ## on microglia
# #bigDiff_30_wt.sal_het.sal <- plotOneGroup(data.combined, "30", 'WT.SAL', 'HET.SAL', row.names(g30_WT.SAL_HET.SAL[g30_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
# bigDiff_30_wt.sal_het.sal <- plotOneGroup(data.combined, "30", 'WT.SAL', 'HET.SAL', row.names(g30_WT.SAL_HET.SAL) )
# ## combo compared to default
# bigDiff_13_wt.sal_het.lps <- plotOneGroup(data.combined, "13", 'WT.SAL', 'HET.LPS', row.names(g13_WT.SAL_HET.LPS[g13_WT.SAL_HET.LPS$p_val_adj < 0.05,] ) )
# bigDiff_23_wt.sal_het.lps <- plotOneGroup(data.combined, "23", 'WT.SAL', 'HET.LPS', row.names(g23_WT.SAL_HET.LPS[g23_WT.SAL_HET.LPS$p_val_adj < 0.05,] ) )






# Genes of interest to plot
## Make container file to store output
setwd(paste0(outputDir, subDir))
subSubDir <- 'visualizeExpression'
print(subSubDir)
print(file.path(outputDir, subDir, subSubDir))

if (file.exists(subSubDir)){
  setwd(file.path(outputDir, subDir, subSubDir))
} else {
  dir.create(file.path(outputDir, subDir, subSubDir))
  setwd(file.path(outputDir, subDir, subSubDir))
}


# Plot gene expression heatmap

## Function to plot gene expression
plotFeatureHeatmap <- function(data.combined, plotMarkers, myGroup, myTitle, chosen_res, chosen_cc){
  # plot
  p <- FeatureHeatmap(data.combined, features.plot = plotMarkers, group.by = myGroup, pt.size = 0.25, key.position = "top", max.exp = 3)
  
  # save
  ggsave(file = paste0(myTitle, '_', myGroup, '_r', chosen_res, '_cc', chosen_cc, '.pdf'), plot = p, device='pdf')
}


## WHy I originally thought cluster labels were wha tthey were
c13 <- c('Ccnd2', 'Ssbp2', 'Dlx1', 'Dlx2', 'Gad1', 'Gad2')
c15 <- c('Calb2', 'Reln', 'Sst', 'Ccnd2', 'Ssbp2','Mef2c', 'Sox2')
c16 <- c('Calb2', 'Reln', 'Sst', 'Ccnd2', 'Sox2','Pax6')
c18 <- c('Sst', 'Npy', 'Ccnd2', 'Otx2', 'Foxg1')
c19 <- c('Otx2','Ccnd2', 'Foxg1')
c23 <- c('Otx2', 'Lmo2', 'Cntnap2', 'Foxg1', 'Sox2')
c25 <- c('Ccnd2', 'Mef2c', 'Dlx1', 'Dlx2', 'Gad2', 'Cntnap2')
c26 <- c('Otx2', 'Mef2c', 'Lmo2', 'Ssbp2', 'Foxg1', 'Pax6')

FeatureHeatmap(data.combined, features.plot = c13, group.by='cond', pt.size = 0.25, key.position = "top", max.exp = 3)
FeatureHeatmap(data.combined, features.plot = c15, group.by='cond', pt.size = 0.25, key.position = "top", max.exp = 3)
FeatureHeatmap(data.combined, features.plot = c16, group.by='cond', pt.size = 0.25, key.position = "top", max.exp = 3)
FeatureHeatmap(data.combined, features.plot = c18, group.by='cond', pt.size = 0.25, key.position = "top", max.exp = 3)
FeatureHeatmap(data.combined, features.plot = c19, group.by='cond', pt.size = 0.25, key.position = "top", max.exp = 3)
FeatureHeatmap(data.combined, features.plot = c23, group.by='cond', pt.size = 0.25, key.position = "top", max.exp = 3)
FeatureHeatmap(data.combined, features.plot = c25, group.by='cond', pt.size = 0.25, key.position = "top", max.exp = 3)
FeatureHeatmap(data.combined, features.plot = c26, group.by='cond', pt.size = 0.25, key.position = "top", max.exp = 3)

## deploy

# example: plotFeatureHeatmap(data.combined, 'Irx3', "cond", 'FeatureHeatmap_prog_markers_A_', chosen_res, chosen_cc )

# info about 19
#FeatureHeatmap(data.combined, features.plot = 'Tcf7l2', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
#FeatureHeatmap(data.combined, features.plot = 'Wnt3', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
FeatureHeatmap(data.combined, features.plot = c('Ebf1', 'Fgfr1op', 'Hes1', 'Lef1'), group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)

FeatureHeatmap(data.combined, features.plot = c('Irx1', 'Irx3', 'Lhx9', 'Nkx2-1'), group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
FeatureHeatmap(data.combined, features.plot = c('Neurog1', 'Neurog2', 'Otx1', 'Otx2'), group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
FeatureHeatmap(data.combined, features.plot = c('Olig3', 'Pax6', 'Rara', 'Tcf4', 'Wnt1'), group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)

# info about 8
FeatureHeatmap(data.combined, features.plot = c('Fgf8', 'Foxa1', 'Gbx2'), group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
FeatureHeatmap(data.combined, features.plot = c('Lhx5', 'Olig2', 'Six3', 'Tle4'), group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)

# other int
FeatureHeatmap(data.combined, features.plot = 'Wnt8b', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Wnt7a', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Wnt5a', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Wnt3a', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Wnt2b', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Lmx1b', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Bmp2', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Bmp4', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Bmp7', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Calb1', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Cyp1b1', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Dlx2', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Ebf2', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Emx2', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Ebf2', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Ebf3', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Fezf1', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Gad1', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Lhx6', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Lhx1', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Lhx8', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               
               # other clusters
               FeatureHeatmap(data.combined, features.plot = 'Ascl1', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Calb2', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Fezf2', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Fgf12', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Fgfbp3', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Gli1', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Gli2', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Gli3', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Lhx2', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Spry2', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Tal1', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               
               # basically nothing/barely detected
               FeatureHeatmap(data.combined, features.plot = 'Bmp6', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Fgf10', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Fgf1', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Fgf21', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Fgf20', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Nkx2-2', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Pitx2', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Sfrp5', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Shh', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               FeatureHeatmap(data.combined, features.plot = 'Sox14', group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
               
# Visualization of other things that will help me understand Otx2+ cultures 
## Clusters 23/26
### From hem
FeatureHeatmap(data.combined, features.plot = c('Bmp7', 'Rspo1', 'Wnt8b', 'Wnt3a'), group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
FeatureHeatmap(data.combined, features.plot = c('Dkk3', 'Sostdc1'), group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)

### Theory I: Hem + Choroid Plexus 
#### C23
FeatureHeatmap(data.combined, features.plot = c('Otx2', 'Msx1', 'Hes1'), group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
FeatureHeatmap(data.combined, features.plot = c('Hes5', 'Lhx2',  'Foxg1'), group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
#### C26
FeatureHeatmap(data.combined, features.plot = c('Otx2', 'Folr1', 'Foxj1','Kcne2'), group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
FeatureHeatmap(data.combined, features.plot = c('Lhx2', 'Emx1', 'Emx2'), group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)

### THeory II: Midbrain dopaminergic cells
FeatureHeatmap(data.combined, features.plot = c('Lmx1a', 'Calml4', 'Folr1', 'Plxdc2'), group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)

### Theory III: Some kind of interneuron
FeatureHeatmap(data.combined, features.plot = c('Sulf1', 'Clic6', 'Hes1'), group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)

### Important in cilia
FeatureHeatmap(data.combined, features.plot = c('Rsph1', 'Foxj1', 'Hes1'), group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)

## Cluster 8: radial glia markers
FeatureHeatmap(data.combined, features.plot = c('Fabp7', 'Ddah1', 'Slc1a3', 'Pea15a', 'Dbi'), group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)

## Cluster 16: neuron markers
FeatureHeatmap(data.combined, features.plot = c('Trp73', 'Calb2', 'Stmn2', 'Meg3', 'Gap43'), group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
FeatureHeatmap(data.combined, features.plot = c('Ebf1', 'Mapt', 'Islr2', 'Tbr1'), group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)



# Visualization of things relevant to the Cluster 13/GABA story
FeatureHeatmap(data.combined, features.plot = c('Gad1', 'Gad2', 'Nkx2-1'), group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
FeatureHeatmap(data.combined, features.plot = c('Ndn', 'Zfp503', 'Dlx2', 'Ebf1', 'Zeb2', 'Maz'), group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)

## accompanying violin Plot
data.combined@meta.data$cond_clust <- paste0(data.combined@meta.data$cond, '_', data.combined@meta.data$res.1.2)

data.combined <- SetAllIdent(data.combined, "res.1.2")
data.combined_c13 <- SubsetData(data.combined, ident.use = '13')
VlnPlot(object = data.combined_c13, features.plot = c('Ebf1', 'Dlx2'), group.by = 'cond')

data.combined_c13 <- SetAllIdent(data.combined_c13, "cond")
data.combined_c13_justHet <- SubsetData(data.combined_c13, ident.use = c('WT.SAL', 'HET.SAL') )
VlnPlot(object = data.combined_c13_justHet, features.plot = c('Ebf1', 'Egr1', 'Zfp503', 'Ndn'), group.by = 'cond', nCol=2)
VlnPlot(object = data.combined_c13_justHet, features.plot = c('Maz', 'Aldoa', 'Dlx2', 'Ypel3'), group.by = 'cond', nCol=2)

## Dlx2 looks weird. What if we could cells with expression above threshold...?
data.combined_c13_normExprs <- FetchData(data.combined_c13, vars.all = c('Ebf1', 'Egr1', 'Zfp503', 'Ndn',
                                                                         'Maz', 'Aldoa', 'Dlx2', 'Ypel3'), 
                                                 use.raw=FALSE)
data.combined_c13_normExprs <- merge(data.combined_c13_normExprs, data.combined_c13@meta.data, by = 'row.names')

library(dplyr)
data.combined_c13_normExprs$Dlx2AboveThresh <- data.combined_c13_normExprs$Dlx2 > 0
tally_dlx <- as.tbl(data.frame(data.combined_c13_normExprs)) %>% group_by(Dlx2AboveThresh, cond) %>% tally()
mean_dlx <- as.tbl(data.frame(data.combined_c13_normExprs)) %>% group_by(cond) %>% dplyr::summarize(Mean = mean(Dlx2, na.rm=TRUE), SD = sd(Dlx2, na.rm=TRUE))


# plot individual genes and save
## function to plot individual genes w/violin plot
plotVlnPlot <- function(data.combined, plotGene, myGroup, myTitle, chosen_res, chosen_cc){
  # plot
  p <- VlnPlot(object = data.combined, features.plot = plotGene, group.by = myGroup, 
               do.return = TRUE)
  # save
  ggsave(file = paste0(myTitle, '_', myGroup, '_', plotGene, '_r', chosen_res, '_cc', chosen_cc, '.pdf'), plot = p, device='pdf')
  }
              
               ## deploy violin plot
               plotVlnPlot(data.combined, 'Hbb-y', 'cond', 'plotVlnPlot', chosen_res, chosen_cc)
               plotVlnPlot(data.combined, 'Hbb-y', 'sample', 'plotVlnPlot', chosen_res, chosen_cc)
               plotVlnPlot(data.combined, 'Gm10076', 'sample', 'plotVlnPlot', chosen_res, chosen_cc)
               plotVlnPlot(data.combined, 'Gm8730', 'sample', 'plotVlnPlot', chosen_res, chosen_cc)
               plotVlnPlot(data.combined, 'Gm10709', 'sample', 'plotVlnPlot', chosen_res, chosen_cc)
               plotVlnPlot(data.combined, 'Gm10116', 'sample', 'plotVlnPlot', chosen_res, chosen_cc)
               plotVlnPlot(data.combined, 'Rpl35', 'sample', 'plotVlnPlot', chosen_res, chosen_cc)
               plotVlnPlot(data.combined, 'Rpl29', 'sample', 'plotVlnPlot', chosen_res, chosen_cc)
               plotVlnPlot(data.combined, 'Hist1h2ap', 'sample', 'plotVlnPlot', chosen_res, chosen_cc)
               plotVlnPlot(data.combined, 'Fabp7', 'sample', 'plotVlnPlot', chosen_res, chosen_cc)
               
               
               
               # # 3.17.19: Visualize interneuron clusters that might change
               # data.combined <- SetAllIdent(data.combined, id = 'res.1.2')
               # data.combined@meta.data$res.1.2 <- as.factor(data.combined@meta.data$res.1.2 )
               # #plotFeatureHeatmap(data.combined, c('13','16','25','18','19','23'), "cond", 'FeatureHeatmap_clusters_justInhib_', chosen_res, chosen_cc )
               # FeatureHeatmap(data.combined, features.plot = 'PC1', group.by = "cond", pt.size = 0.25, key.position = "top", max.exp = 3)
               TSNEPlot(data.combined, group.by = "res.1.2", cells.use = row.names(data.combined@meta.data[which(data.combined@meta.data$cond == 'WT.SAL'), ]), 
                        colors.use = as.character(stanColors[stanColors$cond == 'WT.SAL', 'colors']))
               TSNEPlot(data.combined, group.by = "cond", cells.use = row.names(data.combined@meta.data[which(data.combined@meta.data$cond == 'WT.LPS'), ]), 
                        colors.use = as.character(stanColors[stanColors$cond == 'WT.LPS', 'colors']))
               TSNEPlot(data.combined, group.by = "cond", cells.use = row.names(data.combined@meta.data[which(data.combined@meta.data$cond == 'HET.SAL'), ]),
                        colors.use = as.character(stanColors[stanColors$cond == 'HET.SAL', 'colors']))
               TSNEPlot(data.combined, group.by = "cond", cells.use = row.names(data.combined@meta.data[which(data.combined@meta.data$cond == 'HET.LPS'), ]),
                        colors.use = as.character(stanColors[stanColors$cond == 'HET.LPS', 'colors']))
               
               
               TSNEPlot(data.combined, group.by = "res.1.2", cells.use = row.names(data.combined@meta.data[which(data.combined@meta.data$cond == 'WT.SAL'), ]))
               TSNEPlot(data.combined, group.by = "res.1.2", cells.use = row.names(data.combined@meta.data[which(data.combined@meta.data$cond == 'WT.LPS'), ]))
               TSNEPlot(data.combined, group.by = "res.1.2", cells.use = row.names(data.combined@meta.data[which(data.combined@meta.data$cond == 'HET.SAL'), ]))
               TSNEPlot(data.combined, group.by = "res.1.2", cells.use = row.names(data.combined@meta.data[which(data.combined@meta.data$cond == 'HET.LPS'), ]))
               
               
               
               # save variables
               print('~*~')
               #print('Saving variables...')
               #print(paste0('System time: ', Sys.time()))
               
               #setwd(paste0(outputDir, subDir))
               #save.image(file = paste0("allVars.RData"))
               
               print('~*~ All done! ~*~')
               print(paste0('System time: ', Sys.time()))