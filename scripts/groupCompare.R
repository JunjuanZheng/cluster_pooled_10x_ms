# Kristin Muench
# 2019.03.28
# Pipeline to run Seurat setup for mouse scRNA-Seq data - with CCA
# This takes in a matrix with CCAs calculated, adjusts them, makes a tSNE, 
# and visualizes/counts clusters.

# Load needed library
#install.packages('Seurat')

# Set up workspace
print('Run CCA Track: groupCompare.R')
print('Setting up workspace...')
##remove all environmental vars for clarity
rm(list=ls())
## Load older version of Seurat which has the functions I need
library('Seurat', lib.loc='/scg/apps/software/r/3.5.0/scg/seurat_2.3')
packageVersion('Seurat')

## import from command line
args <- commandArgs(TRUE)
inputFile <- args[1]
outputDir <- args[2]
chosen_cc <- args[3]
chosen_res <- args[4]
print(paste0('Importing files from: ', inputFile))
print(paste0('Output location: ', outputDir))
print(paste0('Chosen # of CCs: ', chosen_cc))
print(paste0('Chosen resolution: ', chosen_res))


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

setwd(paste0(outputDir, subDir))

print('Names of metadata')
print(names(data.combined@meta.data))


# Declare Stanford colors for prettier visualization
stanColors <- data.frame(cond = c('WT.SAL', 'WT.LPS', 'HET.SAL', 'HET.LPS'),
                         colors = as.character(c('#009b76', '#8c1515', '#b26f16', '#53284f')))

# DE Analysis

## Make container file to store output
setwd(paste0(outputDir, subDir))
subSubDir <- 'differentialExpression'
print(subSubDir)
print(file.path(outputDir, subDir, subSubDir))

if (file.exists(subSubDir)){
  setwd(file.path(outputDir, subDir, subSubDir))
} else {
  dir.create(file.path(outputDir, subDir, subSubDir))
  setwd(file.path(outputDir, subDir, subSubDir))
}

## create new group in metadata for both cell type AND cond info
data.combined@meta.data$celltype.cond <- paste0(data.combined@ident, "_",
                                                data.combined@meta.data$cond)

## stash those labels for later use
data.combined <- StashIdent(data.combined, save.name = "celltype")
data.combined <- SetAllIdent(data.combined, id = "celltype.cond") # tell Seurat we are focusing on this grouping method now

## DE between two groups
deGroups <- function(data.combined, xAxis, yAxis){
  # create title
  myTitle <- paste0(xAxis, '_vs_', yAxis)
  print(paste0('DE genes for ', myTitle, '...'))
  
  # do DE
  ident1_xAxis_yAxis <- FindMarkers(data.combined, ident.1 = xAxis, ident.2 = yAxis, print.bar = TRUE)
  head(ident1_xAxis_yAxis, 15)
  
  ident1_xAxis_yAxis
  
  # save results
  write.csv(ident1_xAxis_yAxis, file = paste0(myTitle, '_r', chosen_res, '_cc', chosen_cc, '.csv'), quote=FALSE)
  
  # return DE lists
  return(ident1_xAxis_yAxis)
}

## effect of LPS
g0_WT.SAL_WT.LPS <- deGroups(data.combined, "0_WT.SAL", "0_WT.LPS")
g1_WT.SAL_WT.LPS <- deGroups(data.combined, "1_WT.SAL", "1_WT.LPS")
g2_WT.SAL_WT.LPS <- deGroups(data.combined, "2_WT.SAL", "2_WT.LPS")
g3_WT.SAL_WT.LPS <- deGroups(data.combined, "3_WT.SAL", "3_WT.LPS")
g4_WT.SAL_WT.LPS <- deGroups(data.combined, "4_WT.SAL", "4_WT.LPS")
## on microglia
g30_WT.SAL_WT.LPS <- deGroups(data.combined, "30_WT.SAL", "30_WT.LPS")


## effect of genotype
g0_WT.SAL_HET.SAL <- deGroups(data.combined, "0_WT.SAL", "0_HET.SAL")
g1_WT.SAL_HET.SAL <- deGroups(data.combined, "1_WT.SAL", "1_HET.SAL")
g4_WT.SAL_HET.SAL <- deGroups(data.combined, "4_WT.SAL", "4_HET.SAL")
g7_WT.SAL_HET.SAL <- deGroups(data.combined, "7_WT.SAL", "7_HET.SAL")
### on neurons
g9_WT.SAL_HET.SAL <- deGroups(data.combined, "9_WT.SAL", "9_HET.SAL")
g14_WT.SAL_HET.SAL <- deGroups(data.combined, "14_WT.SAL", "14_HET.SAL")
### on interneurons
g13_WT.SAL_HET.SAL <- deGroups(data.combined, "13_WT.SAL", "13_HET.SAL")
g17_WT.SAL_HET.SAL <- deGroups(data.combined, "17_WT.SAL", "17_HET.SAL")
g20_WT.SAL_HET.SAL <- deGroups(data.combined, "20_WT.SAL", "20_HET.SAL")
g21_WT.SAL_HET.SAL <- deGroups(data.combined, "21_WT.SAL", "21_HET.SAL")
g22_WT.SAL_HET.SAL <- deGroups(data.combined, "22_WT.SAL", "22_HET.SAL")
g23_WT.SAL_HET.SAL <- deGroups(data.combined, "23_WT.SAL", "23_HET.SAL")
g25_WT.SAL_HET.SAL <- deGroups(data.combined, "25_WT.SAL", "25_HET.SAL")
## on microglia
g30_WT.SAL_HET.SAL <- deGroups(data.combined, "30_WT.SAL", "30_HET.SAL")

## effect of both
g13_WT.SAL_HET.LPS <- deGroups(data.combined, "13_WT.SAL", "13_HET.LPS")
g23_WT.SAL_HET.LPS <- deGroups(data.combined, "23_WT.SAL", "23_HET.LPS")


## comparison between groups
g13_WT.SAL_WT.SAL <- deGroups(data.combined, "13_WT.SAL", "23_WT.SAL")

## !!!! Make sure this isnt the effect of one outlier that's high in Hmgb1 or something


## Set back identity
data.combined <- SetAllIdent(data.combined, id = "res.1.2") # tell Seurat we are focusing on the cluster names now



# Differential Expression Analyses

## When comparing two samples, what genes have a large difference?

## Make container file to store output
setwd(paste0(outputDir, subDir))
subSubDir <- 'plotDifferencesInMeans'
print(subSubDir)
print(file.path(outputDir, subDir, subSubDir))

if (file.exists(subSubDir)){
  setwd(file.path(outputDir, subDir, subSubDir))
} else {
  dir.create(file.path(outputDir, subDir, subSubDir))
  setwd(file.path(outputDir, subDir, subSubDir))
}

## Function
plotOneGroup <- function(data.combined, ident1, xAxis, yAxis, genesToLabel){
  # name for this run
  myTitle <- paste0(xAxis, '_vs_', yAxis, '_grp', ident1)
  print(paste0('Plotting genes for ', myTitle, '...'))
  
  # get plotting data for first
  cells1 <- SubsetData(data.combined, ident.use = ident1, subset.raw = T)
  cells1 <- SetAllIdent(cells1, id = "cond")
  avg.cells1 <- log1p(AverageExpression(cells1, show.progress = FALSE))
  avg.cells1$gene <- rownames(avg.cells1)
  
  # clarify what genes to plot
  avg.cells1$plotGene <- "" # effectively hiding all labels
  avg.cells1[genesToLabel, 'plotGene'] <- rownames(avg.cells1[genesToLabel,])
  
  # plot
  library(ggrepel)
  p1 <- ggplot(avg.cells1, aes_string(x = xAxis, y = yAxis, label = 'plotGene')) + geom_point() + ggtitle(ident1) + geom_text_repel()
  
  ggsave(file = paste0(myTitle, '_r', chosen_res, '_cc', chosen_cc, '.pdf'), plot = p1, device='pdf')
  
  # what genes have a large difference in expression?
  bigDiff <- data.frame(abs(avg.cells1[, xAxis] - avg.cells1[, yAxis]))
  row.names(bigDiff) <- row.names(avg.cells1)
  # write data to file
  write.csv(bigDiff, file = paste0(myTitle, '_r', chosen_res, '_cc', chosen_cc, '.csv'), quote=FALSE)
  
  return(bigDiff)
  
}

## Effect of LPS
bigDiff_0_wt.sal_wt.lps <- plotOneGroup(data.combined, "0", 'WT.SAL', 'WT.LPS', row.names(g0_WT.SAL_WT.LPS[g0_WT.SAL_WT.LPS$p_val_adj < 0.05,]) )
bigDiff_1_wt.sal_wt.lps <- plotOneGroup(data.combined, "1", 'WT.SAL', 'WT.LPS', row.names(g1_WT.SAL_WT.LPS[g1_WT.SAL_WT.LPS$p_val_adj < 0.05,]) )
bigDiff_2_wt.sal_wt.lps <- plotOneGroup(data.combined, "2", 'WT.SAL', 'WT.LPS', row.names(g2_WT.SAL_WT.LPS[g2_WT.SAL_WT.LPS$p_val_adj < 0.05,]) )
bigDiff_3_wt.sal_wt.lps <- plotOneGroup(data.combined, "3", 'WT.SAL', 'WT.LPS', row.names(g3_WT.SAL_WT.LPS[g3_WT.SAL_WT.LPS$p_val_adj < 0.05,] ) )
bigDiff_4_wt.sal_wt.lps <- plotOneGroup(data.combined, "4", 'WT.SAL', 'WT.LPS', row.names(g4_WT.SAL_WT.LPS[g4_WT.SAL_WT.LPS$p_val_adj < 0.05,] ) )
## on microglia
#bigDiff_30_wt.sal_wt.lps <- plotOneGroup(data.combined, "30", 'WT.SAL', 'WT.LPS', row.names(g30_WT.SAL_WT.LPS[g30_WT.SAL_WT.LPS$p_val_adj < 0.05,] ) )
bigDiff_30_wt.sal_wt.lps <- plotOneGroup(data.combined, "30", 'WT.SAL', 'WT.LPS', row.names(g30_WT.SAL_WT.LPS) )


# effect of genotype
## progenitors
bigDiff_0_wt.sal_het.sal <- plotOneGroup(data.combined, "0", 'WT.SAL', 'HET.SAL', row.names(g0_WT.SAL_HET.SAL[g0_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
bigDiff_1_wt.sal_het.sal <- plotOneGroup(data.combined, "1", 'WT.SAL', 'HET.SAL', row.names(g1_WT.SAL_HET.SAL[g1_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
bigDiff_4_wt.sal_het.sal <- plotOneGroup(data.combined, "4", 'WT.SAL', 'HET.SAL', row.names(g4_WT.SAL_HET.SAL[g4_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
bigDiff_7_wt.sal_het.sal <- plotOneGroup(data.combined, "7", 'WT.SAL', 'HET.SAL', row.names(g7_WT.SAL_HET.SAL[g7_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
## on neurons
bigDiff_9_wt.sal_het.sal <- plotOneGroup(data.combined, "9", 'WT.SAL', 'HET.SAL', row.names(g9_WT.SAL_HET.SAL[g9_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
bigDiff_14_wt.sal_het.sal <- plotOneGroup(data.combined, "14", 'WT.SAL', 'HET.SAL', row.names(g14_WT.SAL_HET.SAL[g14_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
## on interneurons - not all of these signif
bigDiff_13_wt.sal_het.sal <- plotOneGroup(data.combined, "13", 'WT.SAL', 'HET.SAL', row.names(g13_WT.SAL_HET.SAL[g13_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
bigDiff_17_wt.sal_het.sal <- plotOneGroup(data.combined, "17", 'WT.SAL', 'HET.SAL', row.names(g17_WT.SAL_HET.SAL[g17_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
bigDiff_20_wt.sal_het.sal <- plotOneGroup(data.combined, "20", 'WT.SAL', 'HET.SAL', row.names(g20_WT.SAL_HET.SAL[g20_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
bigDiff_21_wt.sal_het.sal <- plotOneGroup(data.combined, "21", 'WT.SAL', 'HET.SAL', row.names(g21_WT.SAL_HET.SAL[g21_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
bigDiff_22_wt.sal_het.sal <- plotOneGroup(data.combined, "22", 'WT.SAL', 'HET.SAL', row.names(g22_WT.SAL_HET.SAL[g22_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
bigDiff_23_wt.sal_het.sal <- plotOneGroup(data.combined, "23", 'WT.SAL', 'HET.SAL', row.names(g23_WT.SAL_HET.SAL[g23_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
bigDiff_25_wt.sal_het.sal <- plotOneGroup(data.combined, "25", 'WT.SAL', 'HET.SAL', row.names(g25_WT.SAL_HET.SAL[g25_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
## on microglia
#bigDiff_30_wt.sal_het.sal <- plotOneGroup(data.combined, "30", 'WT.SAL', 'HET.SAL', row.names(g30_WT.SAL_HET.SAL[g30_WT.SAL_HET.SAL$p_val_adj < 0.05,] ) )
bigDiff_30_wt.sal_het.sal <- plotOneGroup(data.combined, "30", 'WT.SAL', 'HET.SAL', row.names(g30_WT.SAL_HET.SAL) )
## combo compared to default
bigDiff_13_wt.sal_het.lps <- plotOneGroup(data.combined, "13", 'WT.SAL', 'HET.LPS', row.names(g13_WT.SAL_HET.LPS[g13_WT.SAL_HET.LPS$p_val_adj < 0.05,] ) )
bigDiff_23_wt.sal_het.lps <- plotOneGroup(data.combined, "23", 'WT.SAL', 'HET.LPS', row.names(g23_WT.SAL_HET.LPS[g23_WT.SAL_HET.LPS$p_val_adj < 0.05,] ) )






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

## 16p genes
genesToPlot_16p <- c("Bola2", "Qprt", "Maz", "Mvp", "Cdipt","Sez6l2", "Asphd1", "Kctd13", "Tmem219", "Taok2", "Hirip3", "Ino80e", 
                     "Doc2a","Fam57b","Aldoa", "Ppp4c", "Tbx6", "Ypel3", "Gdpd3", "Mapk3", "Coro1a")
genesToPlot_16p_best_A <- c("Sez6l2", "Kctd13", "Hirip3", "Fam57b", "Aldoa", "Mapk3")
genesToPlot_16p_best_B <- c("Bola2", "Maz", "Cdipt", "Ppp4c", "Ypel3", "Coro1a")

## Other genes of interest
general_markers_prog <- c('Sox2', 'Pax6', 'Aldoc', 'Eomes', 'Axin2')
general_markers_YN <- c('Dcx', 'Bcl11b', 'Satb2')
general_markers_inhib <- c('Sst', 'Otx2', 'Reln', 'Gad1', 'Gad2', 'Dlx1', 'Dlx2')
general_markers_inhib_2 <- c('Pvalb', 'Calb1', 'Calb2', 'Nkx2-1', 'Ascl1', 'Gsx2')

## other
thingsMightBeDifferentBasedEarlierAnalysis <- c("Gm8730", "Rplp0", "Glo1", "Ldha", 'Erdr1', 'Gm10709', 'Htr3a') # can't find 'Rpl129'
sexMarkers <-  c("Xist", "Ddx3y", "Eif2s3y")
regByLPS <- c('Rpl29', 'Rpl35', 'Ubc', 'Gm10709', 'Gm10076', 'Hist1h2ap', 'Hes5', 'Jund', 'Gm8730')
regByLPS_small <- c('Rpl29', 'Gm10709', 'Hist1h2ap', 'Rpl35')
regByLPS_RGCsubgroup <- c('Hmgn2', 'Gm8186', 'Fabp7', 'Rhob', 'Jund', 'Atf3', 'Fam32a')


# Plot gene expression heatmap

## Function to plot gene expression
plotFeatureHeatmap <- function(data.combined, plotMarkers, myGroup, myTitle, chosen_res, chosen_cc){
  # plot
  p <- FeatureHeatmap(data.combined, features.plot = plotMarkers, group.by = myGroup, pt.size = 0.25, key.position = "top", max.exp = 3)
  
  # save
  ggsave(file = paste0(myTitle, '_', myGroup, '_r', chosen_res, '_cc', chosen_cc, '.pdf'), plot = p, device='pdf')
}

## deploy
### different types of markers
plotFeatureHeatmap(data.combined, general_markers_prog, "cond", 'FeatureHeatmap_prog_markers_A_', chosen_res, chosen_cc )
plotFeatureHeatmap(data.combined, general_markers_YN, "cond", 'FeatureHeatmap_YN_markers_B_', chosen_res, chosen_cc )
plotFeatureHeatmap(data.combined, general_markers_inhib, "cond", 'FeatureHeatmap_inhib_markers_C_', chosen_res, chosen_cc )
plotFeatureHeatmap(data.combined, general_markers_inhib_2, "cond", 'FeatureHeatmap_inhib_markers_C_2_', chosen_res, chosen_cc )

### 16p genes
plotFeatureHeatmap(data.combined, genesToPlot_16p[1:5], "cond", 'FeatureHeatmap_16p_A_', chosen_res, chosen_cc )
plotFeatureHeatmap(data.combined, genesToPlot_16p[6:10], "cond", 'FeatureHeatmap_16p_B_', chosen_res, chosen_cc )
plotFeatureHeatmap(data.combined, genesToPlot_16p[11:15], "cond", 'FeatureHeatmap_16p_C_', chosen_res, chosen_cc )
plotFeatureHeatmap(data.combined, genesToPlot_16p[16:21], "cond", 'FeatureHeatmap_16p_D_', chosen_res, chosen_cc )

plotFeatureHeatmap(data.combined, genesToPlot_16p_best_A, "cond", 'FeatureHeatmap_16p_best_A', chosen_res, chosen_cc )
plotFeatureHeatmap(data.combined, genesToPlot_16p_best_B, "cond", 'FeatureHeatmap_16p_best_B', chosen_res, chosen_cc )

### other
plotFeatureHeatmap(data.combined, thingsMightBeDifferentBasedEarlierAnalysis, "cond", 'FeatureHeatmap_mightBeDE_', chosen_res, chosen_cc )
plotFeatureHeatmap(data.combined, sexMarkers, "cond", 'FeatureHeatmap_sexMarkers_', chosen_res, chosen_cc )
plotFeatureHeatmap(data.combined, regByLPS, "cond", 'FeatureHeatmap_regByLPS', chosen_res, chosen_cc )
plotFeatureHeatmap(data.combined, regByLPS_small, "cond", 'FeatureHeatmap_regByLPS_small', chosen_res, chosen_cc )
plotFeatureHeatmap(data.combined, regByLPS_RGCsubgroup, "cond", 'FeatureHeatmap_regByLPS_RGCsubgroup', chosen_res, chosen_cc )


### are DE results driven by outliers?
plotFeatureHeatmap(data.combined, general_markers_YN, "sample", 'FeatureHeatmap_YN_markers_sample_', chosen_res, chosen_cc )
plotFeatureHeatmap(data.combined, thingsMightBeDifferentBasedEarlierAnalysis, "sample", 'FeatureHeatmap_mightBeDE_sample_', chosen_res, chosen_cc )
plotFeatureHeatmap(data.combined, sexMarkers, "sample", 'FeatureHeatmap_sexMarkers_sample_', chosen_res, chosen_cc )

## Plot long sex markers
### plot
p_sex <- FeatureHeatmap(data.combined, features.plot = sexMarkers, group.by = "sample", pt.size = 0.25, key.position = "top", max.exp = 3)
### save
ggsave(file = paste0('FeatureHeatmap_sexMarkers_long_sample_r', chosen_res, '_cc', chosen_cc, '.pdf'), plot = p_sex, device='pdf', width=36, height=16, units="cm")


# plot individual genes
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