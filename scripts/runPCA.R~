# Kristin Muench
# 2019.02.06
# Pipeline to import Seurat files and perform PCA on them. Separated from earlier scripts for ease of troubleshooting.

# install older version of Seurat which has the functions I need
library('Seurat', lib.loc='/scg/apps/software/r/3.5.0/scg/seurat_2.3')
packageVersion('Seurat')

# load needed files
print('Loading needed files...')
args <- commandArgs(TRUE) # load in arguments that accompanied script
setwd(args[1])
load('data_combined_normAndScale.RData') # load combined data
load('metadataFiles.RData') # load combined data
load('hv_genes.RData') # load combined data

data.combined <- data.combined.normAndScale # rename this variable since it's the one we'll use from now on

# Run Principal Components Analysis
## adapted for large datasets - just use hv.genes. Computing many PCs because default (20) doesn't seem enough
print('Running PCAs...')
data.combined <- RunPCA(object = data.combined, pc.genes = hv.genes, do.print = TRUE, pcs.print = 1:10, genes.print = 5, pcs.compute = 40)


# Basic PCA examination
print('Plotting PCAs...')

## Version 1: print PCA
png('printPCA.png', width = 20, height = 20)
PrintPCA(object = data.combined, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
dev.off()

## Version 2: viz PCA
png('VizPCA.png', width = 20, height = 20)
VizPCA(object = data.combined, pcs.use = 1:2)
dev.off()

## Version 3: PCA plot
png('PCAPlot.png', width = 20, height = 20)
PCAPlot(object = data.combined, dim.1 = 1, dim.2 = 2)
PCAPlot(object = data.combined, dim.1 = 3, dim.2 = 4)
dev.off()

## Are there markers strongly associated with cellular heterogeneity that have not passed through variable gene selection?
data.combined <- ProjectPCA(object = data.combined, do.print = FALSE) #The results of the projected PCA can be explored by setting use.full=T in the functions above


# Batch Effects Check: PCA
makePCAPic <- function(seuratObject, identVar, dim1, dim2, myTitle){
  # save as png
  png(myTitle, width = 20, height = 20)
  
  # build tree, which makes a plot
  PCAPlot(object = seuratObject, group.by = identVar, dim.1 = dim1, dim.2 = dim2)

  # close file
  dev.off()
}

print('Generating PCA plots for various possible batch effects...')
setwd(args[1])
makePCAPic(data.combined, 'orig.ident', 1, 2, 'testDataPCA_pc1_pc2_origIdent.png')
makePCAPic(data.combined, 'Litter', 1, 2, 'testDataPCA_pc1_pc2_Litter.png')
makePCAPic(data.combined, 'Genotype', 1, 2, 'testDataPCA_pc1_pc2_Genotype.png')
makePCAPic(data.combined, 'Condition', 1, 2, 'testDataPCA_pc1_pc2_Condition.png')
makePCAPic(data.combined, 'SurgeryDate', 1, 2, 'testDataPCA_pc1_pc2_SurgeryDate.png')
makePCAPic(data.combined, 'CellsPerSample', 1, 2, 'testDataPCA_pc1_pc2_CellsPerSample.png')



# Batch effects check: cluster tree

## function
makeDendrogramPic <- function(seuratObject, identVar, pcs, myTitle){
  seuratObject <- SetAllIdent(object = seuratObject, id = identVar)
  
  # save as png
  png(myTitle, width = 20, height = 20)
  
  # build tree, which makes a plot
  seuratObject <- BuildClusterTree(seuratObject, pcs.use = pcs, do.reorder=TRUE)
  
  # close file
  dev.off()
}

## deploy
print('Generating Dendrograms for various possible batch effects...')
setwd(args[1])
makeDendrogramPic(data.combined, 'Condition', c(1:36), 'clusterTestData_Condition.png') # Condition
makeDendrogramPic(data.combined, 'Genotype', c(1:36), 'clusterTestData_Genotype.png') # Genotype
makeDendrogramPic(data.combined, 'Litter', c(1:36), 'clusterTestData_Litter.png') # Genotype
makeDendrogramPic(data.combined, 'SurgeryDate', c(1:36), 'clusterTestData_SurgeryDate.png') # SurgeryDate
makeDendrogramPic(data.combined, 'CellsPerSample', c(1:36), 'clusterTestData_CellsPerSample.png') # CellsPerSample

print('~*~ Complete ~*~')