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
setwd('RData')
load('data_combined_normAndScale.RData') # load combined data
load('metadataFiles.RData') # load combined data
load('hv_genes.RData') # load combined data
setwd('..')

data.combined <- data.combined.normAndScale # rename this variable since it's the one we'll use from now on


# Create container folder for output of	    this script
subdir <- 'PCA'

if (file.exists(subDir)){
    setwd(file.path(mainDir, subDir))
} else {
    dir.create(file.path(mainDir, subDir))
    setwd(file.path(mainDir, subDir))

}


# Run Principal Components Analysis
## adapted for large datasets - just use hv.genes. Computing many PCs because default (20) doesn't seem enough
setwd(subdir)

print('Running PCAs...')
data.combined <- RunPCA(object = data.combined, pc.genes = hv.genes, do.print = TRUE, pcs.print = 1:10, genes.print = 5, pcs.compute = 40)


# Basic PCA examination
print('Plotting PCAs...')
#par("mar")
#par(c(10,10,10,10))

## Version 1: print PCA
PrintPCA(object = data.combined, pcs.print = 1:5, genes.print = 5, use.full = FALSE) # print list of genes that most strongly define Principal Components

## Version 2: viz PCA
pdf('VizPCA.pdf', width = 11, height = 8)
VizPCA(object = data.combined, pcs.use = 1:2)
dev.off()

## Version 3: PCA plot
png('PCAPlot.pdf', width = 11, height = 8)
PCAPlot(object = data.combined, dim.1 = 1, dim.2 = 2)
PCAPlot(object = data.combined, dim.1 = 3, dim.2 = 4)
dev.off()

## Are there markers strongly associated with cellular heterogeneity that have not passed through variable gene selection?
data.combined <- ProjectPCA(object = data.combined, do.print = FALSE) #The results of the projected PCA can be explored by setting use.full=T in the functions above


# Batch Effects Check: PCA
dev.new(width=10, height=10)
makePCAPic <- function(seuratObject, identVar, dim1, dim2, myTitle){
  # save as pdf
  pdf(myTitle, width = 10, height = 10)
  
  # build tree, which makes a plot
  PCAPlot(object = seuratObject, group.by = identVar, dim.1 = dim1, dim.2 = dim2)

  # close file
  dev.off()
}

print('Generating PCA plots for various possible batch effects...')
setwd(args[1])
makePCAPic(data.combined, 'orig.ident', 1, 2, 'testDataPCA_pc1_pc2_origIdent.pdf')
makePCAPic(data.combined, 'Litter', 1, 2, 'testDataPCA_pc1_pc2_Litter.pdf')
makePCAPic(data.combined, 'Genotype', 1, 2, 'testDataPCA_pc1_pc2_Genotype.pdf')
makePCAPic(data.combined, 'Condition', 1, 2, 'testDataPCA_pc1_pc2_Condition.pdf')
makePCAPic(data.combined, 'SurgeryDate', 1, 2, 'testDataPCA_pc1_pc2_SurgeryDate.pdf')
makePCAPic(data.combined, 'CellsPerSample', 1, 2, 'testDataPCA_pc1_pc2_CellsPerSample.pdf')



# Batch effects check: cluster tree

## function
makeDendrogramPic <- function(seuratObject, identVar, pcs, myTitle){
  seuratObject <- SetAllIdent(object = seuratObject, id = identVar)
  
  # save as png
  pdf(myTitle, width = 10, height = 10)
  
  # build tree, which makes a plot
  seuratObject <- BuildClusterTree(seuratObject, pcs.use = pcs, do.reorder=TRUE)
  
  # close file
  dev.off()
}

## deploy
print('Generating Dendrograms for various possible batch effects...')
setwd(args[1])
makeDendrogramPic(data.combined, 'Condition', c(1:36), 'clusterTestData_Condition.pdf') # Condition
makeDendrogramPic(data.combined, 'Genotype', c(1:36), 'clusterTestData_Genotype.pdf') # Genotype
makeDendrogramPic(data.combined, 'Litter', c(1:36), 'clusterTestData_Litter.pdf') # Genotype
makeDendrogramPic(data.combined, 'SurgeryDate', c(1:36), 'clusterTestData_SurgeryDate.pdf') # SurgeryDate
makeDendrogramPic(data.combined, 'CellsPerSample', c(1:36), 'clusterTestData_CellsPerSample.pdf') # CellsPerSample



# Determine statistically significant PCs
print('Determining statistically significant PCs...')

## Method 1: PC heatmap - useful when deciding what gene sets to use in downstream analyses
### More supervised - exploring PCs to determine relevant sources of heterogeneity, and could be used in conjunction with GSEA for example.
print('Method 1: PC Heatmap...')

pdf('pcStat_pcHeatmap_pc1_pc12.pdf', width=10, height=10)
PCHeatmap(object = data.combined, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
dev.off()

pdf('pcStat_pcHeatmap_pc13_pc24.pdf', width=10, height=10)
PCHeatmap(object = data.combined, pc.use = 13:24, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
dev.off()

pdf('pcStat_pcHeatmap_pc25_pc39.pdf', width=10, height=10)
PCHeatmap(object = data.combined, pc.use = 25:39, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
dev.off()

## Method 3: Elbow plot heuristic
print('Method3: Elbow Plot...')
pdf('pcStat_elbowPlot_pc40.pdf', width=10, height=10)
PCElbowPlot(object = data.combined, num.pc = 40) # for big datasets. In tutorials, it appears they select #PCs where it begins to plateau
dev.off()


#Save variables
setwd('../RData')
data.combined.pca <- data.combined
save(data.combined.pca, file='data_combined_pca.RData')


# !! Doesn't work, would be nice if it did

# Method 2: Resampling procedure - what is the distribution of p-values compared to a predicted uniform distribution? (Above line = significant PC)
## Statistical test based on a random null model, but is time-consuming for large datasets, and may not return a clear PC cutoff
#print('Method 2: JackStraw Plot...')

#data.combined <- JackStraw(object = data.combined, num.pc = 39, num.replicate = 100, display.progress = FALSE) # Calculate
#traceback()

#print('Plotting JackStraw...')
#pdf('pcStat_jackStraw_pc40.pdf', width=10, height=10)
#JackStrawPlot(object = data.combined, PCs = c(1:39) )# View
#dev.off()



print('~*~ Complete ~*~')