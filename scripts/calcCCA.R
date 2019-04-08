# calcCCA.R
# Kristin Muench
# 2019.03.26
# Pipeline to run Seurat setup for mouse scRNA-Seq data - with CCA
# This takes the loaded and filtered variable data and adjusts them with CCA.

# Set up workspace
print('Run CCA Track: calcCCA.R')
print('Setting up workspace...')
##remove all environmental vars for clarity
rm(list=ls())
## Load older version of Seurat which has the functions I need
library('Seurat', lib.loc='/scg/apps/software/r/3.5.0/scg/seurat_2.3')
packageVersion('Seurat')

## run this in case you accidentally load Seurat 3
#detach("package:Seurat", unload=TRUE)

## import from command line
args <- commandArgs(TRUE)
fileLoc <- args[1]
load(fileLoc)
print(paste0('Importing files from: ', fileLoc))

args <- commandArgs(TRUE)
fileLoc <- args[1]
outputDir <- args[2]
print(paste0('Output location2: ', outputDir))
numCCs = 75 # number of CCs to plot and run
print(paste0('Num CCs: ', numCCs))

# # Uncomment for troubleshooting
# fileLoc <- '/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190228_runThroughDemultiplex/makeVars/'
# outputDir <- '/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190228_runThroughDemultiplex/'
# ## load needed files - if troubleshooting/working in RStudio
# setwd(fileLoc)
# load('allVars.RData')

## make subdirectory
setwd(outputDir)
subDir <- 'calcCCA'
print('!!!!!!!!')
print(file.path(outputDir, subDir))

if (file.exists(subDir)){
  setwd(file.path(outputDir, subDir))
} else {
  dir.create(file.path(outputDir, subDir))
  setwd(file.path(outputDir, subDir))
}

setwd(file.path(outputDir, subDir))

# Figure out what genes should be inputs to CCA
print('~*~')
print('Figuring out what genes should be inputs to CCA...')
print(paste0('System time: ', Sys.time()))
## Find variable genes
wt.sal <- FindVariableGenes(wt.sal, do.plot = F)
wt.lps <- FindVariableGenes(wt.lps, do.plot = F)
het.sal <- FindVariableGenes(het.sal, do.plot = F)
het.lps <- FindVariableGenes(het.lps, do.plot = F)

## Glance at what those genes are
print('Highly variable genes:')
g.1 <- head(rownames(wt.sal@hvg.info), 1000)
g.2 <- head(rownames(wt.lps@hvg.info), 1000)
g.3 <- head(rownames(het.sal@hvg.info), 1000)
g.4 <- head(rownames(het.lps@hvg.info), 1000)

## what genes are in common
genes.use <- unique(c(g.1, g.2, g.3, g.4))
genes.use <- intersect(genes.use, rownames(wt.sal@scale.data))
genes.use <- intersect(genes.use, rownames(wt.lps@scale.data))
genes.use <- intersect(genes.use, rownames(het.sal@scale.data))
genes.use <- intersect(genes.use, rownames(het.lps@scale.data))


# Run Canonical Correlation Analysis
print('~*~')
print('Running Canonical Correlation Analysis...')
print(paste0('System time: ', Sys.time()))
setwd(file.path(outputDir, subDir))

objects <- c(wt.sal, wt.lps, het.sal, het.lps)
cellIDs <- c('WT.SAL', 'WT.LPS', 'HET.SAL', 'HET.LPS')

head(wt.sal@meta.data)
head(x = wt.sal@cell.names)

# Run multi CCA to find correlation between the four
data.combined <- RunMultiCCA(objects, genes.use = genes.use, num.ccs = numCCs) #, add.cell.ids = cellIDs)

head(data.combined@meta.data)
head(x = data.combined@cell.names)
table(data.combined@meta.data$orig.ident) # vewrify that it works and you can detect sample types still

## save variable because it takes so long to make
setwd(file.path(outputDir, subDir))
print(paste0('Directory should be ',paste0(outputDir, subDir)) )
setwd(file.path(outputDir, subDir))
save(data.combined, file = paste0("data.combined_multiCCA.RData"))

# # uncomment to just load earlier file ~~~~~~
# print('Loading CCA file you made before...')
# setwd(file.path(outputDir, subDir))
# print(paste0('file path is: ', file.path(outputDir, subDir)))
# load("data.combined_multiCCA.RData")

# add my own metadata to Seurat object
print('~*~')
print('Adding my metadata to Seurat...')

## prepare metadata
## !!!!!!! is this the best way to do this?
# files$ids <- unlist(lapply(strsplit(list.filenames, "_"), '[[', 1))
data.combined@meta.data$origSample <- unlist(lapply(strsplit(as.character(data.combined@meta.data$sample), split="_"), '[[', 1)) # gather original 10x SampleIDs

myCells <- data.frame(cellNames = data.combined@cell.names)
myCells$Identity <- data.combined@meta.data$origSample
metadataCols <- c('SampleID','CellsPerSample','SurgeryDate','Condition', 'Genotype', 'Litter', 'sex')
myMetadata <- merge(myCells, metadata[,metadataCols], by.x='Identity', by.y='SampleID')
rownames(myMetadata) <- myMetadata$cellNames

## add metadata to Seurat object
data.combined <- AddMetaData(data.combined, myMetadata[,-c(1:2)], col.name = metadataCols[-c(1:2)])

setwd(file.path(outputDir, subDir))
save(data.combined, file = paste0("data.combined_multiCCA_metadata.RData")) # save that data


print('Metadata looks like:')
head(data.combined@meta.data)

## If you already ran the above - load the data and start from here
#setwd(file.path(outputDir, subDir))
#load('data.combined_multiCCA.RData')
#load('data.combined_multiCCA_metadata.RData')



# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
print('Making visualizations of CCA results...')
## Function
visResultsCCA <- function(data.combined, myGroup){
  p1 <- DimPlot(object = data.combined, reduction.use = "cca", group.by = myGroup, 
                pt.size = 0.5, do.return = TRUE)
  
  p2 <- VlnPlot(object = data.combined, features.plot = "CC1", group.by = myGroup, 
                do.return = TRUE)
  
  plot_grid(p1, p2)
  
  # save method 2
  ggsave(file = paste0('CCA_DimPlot_', myGroup, '.pdf'), plot = p1, device='pdf')
  ggsave(file = paste0('CCA_VlnPlot_', myGroup, '.pdf'), plot = p2, device='pdf')
}

## Deploy
visResultsCCA(data.combined, "cond")
visResultsCCA(data.combined, "Litter")
visResultsCCA(data.combined, "CellsPerSample")
visResultsCCA(data.combined, "SurgeryDate")
visResultsCCA(data.combined, "Condition")
visResultsCCA(data.combined, "Genotype")


# Choose which CCs to use
print('~*~')
print('Making graphs to choose which CCs to use...')
print(paste0('System time: ', Sys.time()))
## bicor plot - when do we see a dropoff in signal/when does change in x not account for large change in y?

makeBicorPlot <- function(data.combined, myGroup, numCCs){
  p3 <- MetageneBicorPlot(data.combined, grouping.var = myGroup, dims.eval = 1:numCCs, 
                          display.progress = FALSE)
  
  ggsave(file = paste0('bicorPlot_', myGroup, '.pdf'), plot = p3, device='pdf')
}

makeBicorPlot(data.combined, "cond", numCCs)
makeBicorPlot(data.combined, "Litter", numCCs)
makeBicorPlot(data.combined, "CellsPerSample", numCCs)
makeBicorPlot(data.combined, "SurgeryDate", numCCs)
makeBicorPlot(data.combined, "Condition", numCCs)
makeBicorPlot(data.combined, "Genotype", numCCs)

## heatmap to assess which components are actually important?
print('Heatmap to demonstrate which components are actually important...')
pdf(file = 'dimHeatmap_1_9.pdf')
DimHeatmap(object = data.combined, reduction.type = "cca", cells.use = 500, dim.use = 1:9, do.balanced = TRUE)
dev.off()

pdf(file = 'dimHeatmap_10_19.pdf')
DimHeatmap(object = data.combined, reduction.type = "cca", cells.use = 500, dim.use = 10:19, do.balanced = TRUE)
dev.off()

pdf(file = 'dimHeatmap_20_29.pdf')
DimHeatmap(object = data.combined, reduction.type = "cca", cells.use = 500, dim.use = 20:29, do.balanced = TRUE)
dev.off()

pdf(file = 'dimHeatmap_30_39.pdf')
DimHeatmap(object = data.combined, reduction.type = "cca", cells.use = 500, dim.use = 30:39, do.balanced = TRUE)
dev.off()

pdf(file = 'dimHeatmap_40_49.pdf')
DimHeatmap(object = data.combined, reduction.type = "cca", cells.use = 500, dim.use = 40:49, do.balanced = TRUE)
dev.off()



# save variables
print('~*~')
print('Saving variables...')
print(paste0('System time: ', Sys.time()))

setwd(file.path(outputDir, subDir))
save.image(file = paste0("allVars.RData"))

print('~*~ All done! ~*~')
print(paste0('System time: ', Sys.time()))
