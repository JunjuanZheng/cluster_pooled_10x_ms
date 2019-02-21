# Kristin Muench
# 2019.02.12
# Pipeline to run Seurat setup for mouse scRNA-Seq data - with CCA
# This takes in a matrix with CCAs calculated, adjusts them, makes a tSNE, 
# and visualizes/counts clusters.

# Load needed library
#install.packages('Seurat')

# Set up workspace
print('Run CCA Track: groupCompare.R')
print('Setting up workspace...')
## Load older version of Seurat which has the functions I need
library('Seurat', lib.loc='/scg/apps/software/r/3.5.0/scg/seurat_2.3')
packageVersion('Seurat')

## import from command line
args <- commandArgs(TRUE)
#load(args[1])
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


# Genes of interest to plot
## 16p genes
genesToPlot_16p <- c("Bola2", "Qprt", "Maz", "Mvp", "Cdipt","Sez6l2", "Asphd1", "Kctd13", "Tmem219", "Taok2", "Hirip3", "Ino80e", "Doc2a","Fam57b","Aldoa", "Ppp4c", "Tbx6", "Ypel3", "Gdpd3", "Mapk3", "Coro1a")

## Other genes of interest
general_markers_prog <- c('Sox2', 'Pax6', 'Aldoc', 'Eomes', 'Axin2')
general_markers_YN <- c('Dcx', 'Bcl11b', 'Satb2')
general_markers_inhib <- c('Sst', 'Otx2', 'Reln', 'Gad2', 'Dlx1', 'Dlx2')
general_markers_inhib_2 <- c('Pvalb', 'Calb1', 'Nkx2-1', 'Ascl1', 'Gsx2')

# other
thingsMightBeDifferentBasedEarlierAnalysis <- c("Gm8730", "Rplp0", "Glo1", "Ldha", 'Erdr1', 'Gm10709', 'Htr3a') # can't find 'Rpl129'
sexMarkers <-  c("Xist", "Ddx3y", "Eif2s3y")

# # !!!! UNCOMMENT WHEN YOU CAN GET CLUSTERS
# # Visualize gene expression across groups - split dot plot
# 
# ## set ident to factor
# data.combined@ident <- as.factor(data.combined@ident)
# 

# 
# ## Make split plot
# pdf(file = paste0('SplitDotPlot_genesOfInterest_r', chosen_res, '_cc', chosen_cc, '.pdf'))
# SplitDotPlotGG(data.combined, genes.plot = rev(markers.to.plot), cols.use = c("blue", "red", "green", "black"), x.lab.rot = T, plot.legend = T, do.return = T, grouping.var = "stim")
# dev.off()
# 
# pdf(file = paste0('SplitDotPlot_genes16p_r', chosen_res, '_cc', chosen_cc, '.pdf'))
# SplitDotPlotGG(data.combined, genes.plot = rev(genesToPlot_16p), cols.use = c("blue", "red", "green", "black"), x.lab.rot = T, plot.legend = T, do.return = T, grouping.var = "stim")
# dev.off()


# Create a heatmap that 

plotFeatureHeatmap <- function(data.combined, plotMarkers, myGroup, myTitle, chosen_res, chosen_cc){
    # plot
    p <- FeatureHeatmap(data.combined, features.plot = plotMarkers, group.by = myGroup, pt.size = 0.25, key.position = "top", max.exp = 3)
    # save
    ggsave(file = paste0(myTitle, '_', myGroup, '_r', chosen_res, '_cc', chosen_cc, '.pdf'), plot = p, device='pdf')
}

# plot areas
plotFeatureHeatmap(data.combined, general_markers_prog, "stim", 'FeatureHeatmap_prog_markers_A_', chosen_res, chosen_cc )
plotFeatureHeatmap(data.combined, general_markers_YN, "stim", 'FeatureHeatmap_YN_markers_B_', chosen_res, chosen_cc )
plotFeatureHeatmap(data.combined, general_markers_inhib, "stim", 'FeatureHeatmap_inhib_markers_C_', chosen_res, chosen_cc )
plotFeatureHeatmap(data.combined, general_markers_inhib_2, "stim", 'FeatureHeatmap_inhib_markers_C_2_', chosen_res, chosen_cc )

# 16p
plotFeatureHeatmap(data.combined, genesToPlot_16p[1:5], "stim", 'FeatureHeatmap_16p_A_', chosen_res, chosen_cc )
plotFeatureHeatmap(data.combined, genesToPlot_16p[6:10], "stim", 'FeatureHeatmap_16p_B_', chosen_res, chosen_cc )
plotFeatureHeatmap(data.combined, genesToPlot_16p[11:15], "stim", 'FeatureHeatmap_16p_C_', chosen_res, chosen_cc )
plotFeatureHeatmap(data.combined, genesToPlot_16p[16:21], "stim", 'FeatureHeatmap_16p_D_', chosen_res, chosen_cc )

# other
plotFeatureHeatmap(data.combined, thingsMightBeDifferentBasedEarlierAnalysis, "stim", 'FeatureHeatmap_mightBeDE_', chosen_res, chosen_cc )
plotFeatureHeatmap(data.combined, sexMarkers, "stim", 'FeatureHeatmap_sexMarkers_', chosen_res, chosen_cc )

# driven by outliers?
plotFeatureHeatmap(data.combined, general_markers_YN, "sample", 'FeatureHeatmap_YN_markers_sample_', chosen_res, chosen_cc )
plotFeatureHeatmap(data.combined, thingsMightBeDifferentBasedEarlierAnalysis, "sample", 'FeatureHeatmap_mightBeDE_sample_', chosen_res, chosen_cc )
plotFeatureHeatmap(data.combined, sexMarkers, "sample", 'FeatureHeatmap_sexMarkers_sample_', chosen_res, chosen_cc )


# save variables
print('~*~')
#print('Saving variables...')
#print(paste0('System time: ', Sys.time()))

#setwd(paste0(outputDir, subDir))
#save.image(file = paste0("allVars.RData"))

print('~*~ All done! ~*~')
print(paste0('System time: ', Sys.time()))




