# Kristin Muench
# 2019.02.12
# Pipeline to run Seurat setup for mouse scRNA-Seq data - with CCA
# This takes in a matrix with CCAs calculated, adjusts them, makes a tSNE, 
# and visualizes/counts clusters.
# SIMPLER VERSION WITH NO PICTURES TO TRY AND TROUBLESHOOT THE FINDCLUSTERS ERROR

# Set up workspace
print('Run CCA Track: clusters_simp.R')
print('Setting up workspace...')
## Load older version of Seurat which has the functions I need
library('Seurat', lib.loc='/scg/apps/software/r/3.5.0/scg/seurat_2.3')
print('Package Version of Seurat:')
packageVersion('Seurat')

## import from command line
args <- commandArgs(TRUE)
#load(args[1])
inputFile <- args[1]
outputDir <- args[2]
chosen_cc <- args[3]
chosen_res <- args[4]
print(paste0('Importing files : ', inputFile))
print(paste0('Output location: ', outputDir))
print(paste0('Chosen # of CCs: ', chosen_cc))


## load needed files
print('~*~')
print('Loading needed files...')
load(inputFile)

## make subdirectory
setwd(outputDir)
subDir <- 'clusters_simp'
print(paste0('Name of Experiment: ', subDir))

if (file.exists(subDir)){
    setwd(file.path(outputDir, subDir))
} else {
    dir.create(file.path(outputDir, subDir))
    setwd(file.path(outputDir, subDir))
}

setwd(paste0(outputDir, subDir))


# Align CCA subspaces
print('~*~')
print('Align CCA subspaces...')
print(paste0('System time: ', Sys.time()))

myGroupingVar = "stim"
data.combined <- AlignSubspace(data.combined, reduction.type = "cca", grouping.var = myGroupingVar,  dims.align = 1:chosen_cc)

print('Saving aligned subspaces...')
#setwd(paste0(outputDir, subDir))
#save(data.combined, file = paste0('data.combined_aligned_groupby_', myGroupingVar, '_cc_', chosen_cc, '.RData') )

# Run integrated analysis on all cells
print('~*~')
print('Run integrated analysis...')
print(paste0('System time: ', Sys.time()))

## make tSNE
### t-SNE and Clustering - NOTE USE OF 'REDUCTION.USE'
print('Making tSNE...')
print(paste0('System time: ', Sys.time()))
data.combined <- RunTSNE(data.combined, reduction.use = "cca.aligned", dims.use = 1:chosen_cc, do.fast = T)

#save(data.combined, file= paste0('data.combined_withTSNE_r', chosen_res, '_CC', chosen_cc, '.RData') )

## make clusters
print('Finding clusters...')
print(paste0('System time: ', Sys.time()))

data.combined <- FindClusters(data.combined, reduction.type = "cca.aligned", resolution = chosen_res, dims.use = 1:chosen_cc, nn.eps = 0.5) # the problem line
save(data.combined, file= paste0('data.combined_withClust_r', chosen_res, '_CC', chosen_cc, '.RData') )


## Function to visualize tSNE
visTSNE <- function(data.combined, myGroup){
  
  p3 <- TSNEPlot(data.combined, do.return = T, pt.size = 0.5, group.by = myGroup)
  
  # save output
  ggsave(file = paste0( 'tSNE_group_', myGroup, '_res', chosen_res, '_CC', chosen_cc, '.pdf' ), plot = p3, device='pdf')
}

visTSNE(data.combined, "stim")
visTSNE(data.combined, "Litter")
visTSNE(data.combined, "CellsPerSample")
visTSNE(data.combined, "SurgeryDate")
visTSNE(data.combined, "Condition")
visTSNE(data.combined, "Genotype")

### Visualize tSNE but with clusters
p2 <- TSNEPlot(data.combined, do.label = T, do.return = T, pt.size = 0.5)
ggsave(file = paste0( 'tSNE_clusters_res', chosen_res, '_CC', chosen_cc, '.pdf' ), plot = p2, device='pdf')

### Find cell marker types conserved across all groups
print('Finding cell markers conserved across all groups...')
print(paste0('System time: ', Sys.time()))

#### Make container file to store output
setwd(paste0(outputDir, subDir))
subSubDir <- 'marker_csv'
print(subSubDir)
print(file.path(outputDir, subDir, subSubDir))

if (file.exists(subSubDir)){
    setwd(file.path(outputDir, subDir, subSubDir))
} else {
    dir.create(file.path(outputDir, subDir, subSubDir))
    setwd(file.path(outputDir, subDir, subSubDir))
}

#### Helpful function
getMarkerInfo <- function(data.combined, i){

  # get conserved markers
  consMarkers <- FindConservedMarkers(data.combined, ident.1 = i, grouping.var = "stim", 
    print.bar = FALSE)
  print(head(consMarkers))
  
  # write answer
  write.csv(consMarkers, file = paste0('consMarkers_',i,'.csv') )
  
  # visualize
  pdf(file = paste0('FeaturePlot_', i))
  FeaturePlot(object = data.combined, features.plot = row.names(consMarkers[c(1:9),]), min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
  dev.off()
  
  # return results
  return(consMarkers)
}

consMarkers_0 <- getMarkerInfo(data.combined, 0)
consMarkers_1 <- getMarkerInfo(data.combined, 1)
consMarkers_2 <- getMarkerInfo(data.combined, 2)
consMarkers_3 <- getMarkerInfo(data.combined, 3)
consMarkers_4 <- getMarkerInfo(data.combined, 4)
consMarkers_5 <- getMarkerInfo(data.combined, 5)
consMarkers_6 <- getMarkerInfo(data.combined, 6)
consMarkers_7 <- getMarkerInfo(data.combined, 7)
consMarkers_8 <- getMarkerInfo(data.combined, 8)
consMarkers_9 <- getMarkerInfo(data.combined, 9)
consMarkers_10 <- getMarkerInfo(data.combined, 10)
consMarkers_11 <- getMarkerInfo(data.combined, 11)
consMarkers_12 <- getMarkerInfo(data.combined, 12)
consMarkers_13 <- getMarkerInfo(data.combined, 13)


#### plot some of these markers
print('Plotting some of the common markers we look for...')
genesToPlot <- c("Sox2", "Pax6", "Aldoc", "Hes1", 
                 "Eomes", "Neurog1", 
                 "Reln",
                 "Dlx1", "Dlx2", "Otx2", "Nrxn3", 'Sst',
                 "Aif1","Cx3cr1","C1qa",
                 "Neurod1", "Neurod2", "Lmo1", "Mapt", 'Bcl11b', 'Satb2')
pdf(file = paste0('FeaturePlot_FamousGenes'))
FeaturePlot(object = immune.combined, features.plot = genesToPlot, min.cutoff = "q9", cols.use = c("lightgrey",
    "blue"), pt.size = 0.5)
dev.off()







# save variables
print('~*~')
print('Saving variables...')
print(paste0('System time: ', Sys.time()))

setwd(paste0(outputDir, subDir))
save.image(file = paste0("allVars.RData"))

print('~*~ All done! ~*~')
print(paste0('System time: ', Sys.time()))




