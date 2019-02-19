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
print(paste0('Importing files from: ', inputFilePath))
print(paste0('Output location: ', outputDir))
print(paste0('Chosen # of CCs: ', chosen_cc))


## load needed files
print('~*~')
print('Loading needed files...')
setwd(inputFilePath)
load('data.combined_multiCCA_metadata.RData')

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



# # !!!! UNCOMMENT WHEN YOU CAN GET CLUSTERS
# # Visualize gene expression across groups - split dot plot
# 
# ## set ident to factor
# data.combined@ident <- as.factor(data.combined@ident)
# 
# ##Declare genes to plot
# ### 16p genes
# genesToPlot_16p <- c("Bola2", "Qprt", "Maz", "Mvp", "Cdipt","Sez6l2", "Asphd1", "Kctd13", "Tmem219", "Taok2", "Hirip3", "Ino80e", "Doc2a","Fam57b","Aldoa", "Ppp4c", "Tbx6", "Ypel3", "Gdpd3", "Mapk3", "Coro1a")
# 
# ### Other genes of interest
# markers.to.plot <- c('Sox2', 'Pax6', 'Aldoc', 'Dcx', 'Bcl11b', 'Sst', 'Otx2', 'Reln')
# 
# ## Make split plot
# pdf(file = paste0('SplitDotPlot_genesOfInterest_r', chosen_res, '_cc', chosen_cc, '.pdf'))
# SplitDotPlotGG(data.combined, genes.plot = rev(markers.to.plot), cols.use = c("blue", "red", "green", "black"), x.lab.rot = T, plot.legend = T, do.return = T, grouping.var = "stim")
# dev.off()
# 
# pdf(file = paste0('SplitDotPlot_genes16p_r', chosen_res, '_cc', chosen_cc, '.pdf'))
# SplitDotPlotGG(data.combined, genes.plot = rev(genesToPlot_16p), cols.use = c("blue", "red", "green", "black"), x.lab.rot = T, plot.legend = T, do.return = T, grouping.var = "stim")
# dev.off()








# save variables
print('~*~')
print('Saving variables...')
print(paste0('System time: ', Sys.time()))

setwd(paste0(outputDir, subDir))
save.image(file = paste0("allVars.RData"))

print('~*~ All done! ~*~')
print(paste0('System time: ', Sys.time()))




