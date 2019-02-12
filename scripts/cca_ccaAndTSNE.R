# Kristin Muench
# 2019.02.12
# Pipeline to run Seurat setup for mouse scRNA-Seq data - with CCA
# This takes the loaded and filtered variable data and adjusts them with CCA.

# Load needed library
#install.packages('Seurat')

# Set up workspace
print('Setting up workspace...')
## Load older version of Seurat which has the functions I need
library('Seurat', lib.loc='/scg/apps/software/r/3.5.0/scg/seurat_2.3')
packageVersion('Seurat')

## import from command line
args <- commandArgs(TRUE)
#load(args[1])
fileLoc <- paste0(args[1], 'makeVars')
outputDir <- args[1]
print(paste0('Importing files from: ', fileLoc))
print(paste0('Output location: ', outputDir))

## make subdirectory
setwd(outputDir)
subDir <- 'ccaAndTSNE'

if (file.exists(subDir)){
    setwd(file.path(outputDir, subDir))
} else {
    dir.create(file.path(outputDir, subDir))
    setwd(file.path(outputDir, subDir))
}

## load needed files
print('\nLoading needed files...')
setwd(fileLoc)
load('allVars.RData')

setwd(paste0(outputDir, subDir))


# Figure out what genes should be inputs to CCA
print('\nFiguring out what genes should be inputs to CCA...')
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

print('\nRunning Canonical Correlation Analysis...')
print(paste0('System time: ', Sys.time()))

objects <- c(wt.sal, wt.lps, het.sal, het.lps)
cellIDs <- c('WT.SAL', 'WT.LPS', 'HET.SAL', 'HET.LPS')

data.combined <- RunMultiCCA(objects, genes.use = genes.use, num.ccs = 30, add.cell.ids = cellIDs)

head(x = data.combined@cell.names)
table(data.combined@meta.data$orig.ident)

# add my own metadata to Seurat object !!!!! FIX THIS
print('\nAdding my metadata to Seurat...')
myCells <- data.frame(cellNames = data.combined@cell.names)
myCells$Identity <- toupper(data.combined@meta.data$orig.ident)
metadataCols <- c('SampleID','CellsPerSample','SurgeryDate','Condition', 'Genotype', 'Litter')
myMetadata <- merge(myCells, metadata[,metadataCols], by.x='Identity', by.y='SampleID')
rownames(myMetadata) <- myMetadata$cellNames

# add to Seurat object
data.combined <- AddMetaData(data.combined, myMetadata[,-c(1:2)], col.name = metadataCols[-c(1:2)])

## visualize results of CCA plot CC1 versus CC2 and look at a violin plot
print('Making visualizations of CCA results...')
visResultsCCA <- function(data.combined, myGroup){
    p1 <- DimPlot(object = data.combined, reduction.use = "cca", group.by = myGroup, 
    pt.size = 0.5, do.return = TRUE)
    
    p2 <- VlnPlot(object = data.combined, features.plot = "CC1", group.by = "stim", 
    do.return = TRUE)
    
    pdf( paste0('CCA_DimPlot_VlnPlot_grpby_', myGroup, '.pdf') , width=10, height=10)
    plot_grid(p1, p2)
    dev.off()
}

visResultsCCA(data.combined, "stim")











# save variables
setwd(paste0(outputDir, subDir))
save.image(file = paste0("allVars.RData"))

print('~*~ All done! ~*~')