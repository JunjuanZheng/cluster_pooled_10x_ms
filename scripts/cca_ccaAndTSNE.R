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
print('~*~')
print('Loading needed files...')
setwd(fileLoc)
load('allVars.RData')

setwd(paste0(outputDir, subDir))


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

objects <- c(wt.sal, wt.lps, het.sal, het.lps)
cellIDs <- c('WT.SAL', 'WT.LPS', 'HET.SAL', 'HET.LPS')

data.combined <- RunMultiCCA(objects, genes.use = genes.use, num.ccs = 30, add.cell.ids = cellIDs)

head(x = data.combined@cell.names)
table(data.combined@meta.data$orig.ident)

## save variable because it takes so long to make
setwd(paste0(outputDir, subDir))
save(data.combined, file = paste0("data.combined_multiCCA.RData"))

# add my own metadata to Seurat object !!!!! FIX THIS
print('~*~')
print('Adding my metadata to Seurat...')

## prepare metadata
myCells <- data.frame(cellNames = data.combined@cell.names)
myCells$Identity <- toupper(data.combined@meta.data$orig.ident)
metadataCols <- c('SampleID','CellsPerSample','SurgeryDate','Condition', 'Genotype', 'Litter')
myMetadata <- merge(myCells, metadata[,metadataCols], by.x='Identity', by.y='SampleID')
rownames(myMetadata) <- myMetadata$cellNames

## add metadata to Seurat object
data.combined <- AddMetaData(data.combined, myMetadata[,-c(1:2)], col.name = metadataCols[-c(1:2)])


# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
print('Making visualizations of CCA results...')
## Function
visResultsCCA <- function(data.combined, myGroup){
    p1 <- DimPlot(object = data.combined, reduction.use = "cca", group.by = myGroup, 
    pt.size = 0.5, do.return = TRUE)
    
    p2 <- VlnPlot(object = data.combined, features.plot = "CC1", group.by = "stim", 
    do.return = TRUE)
    
    pdf( paste0('CCA_DimPlot_VlnPlot_grpby_', myGroup, '.pdf') , width=10, height=10)
    plot_grid(p1, p2)
    dev.off()
}

## Deploy
visResultsCCA(data.combined, "stim")
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

makeBicorPlot <- function(data.combined, myGroup){
    p3 <- MetageneBicorPlot(data.combined, grouping.var = myGroup, dims.eval = 1:30, 
    display.progress = FALSE)
    
    ggsave(file = paste0('bicorPlot_', myGroup, '.pdf'), plot = p3, device='pdf')
}

makeBicorPlot(data.combined, "stim")
makeBicorPlot(data.combined, "Litter")
makeBicorPlot(data.combined, "CellsPerSample") #!
makeBicorPlot(data.combined, "SurgeryDate")
makeBicorPlot(data.combined, "Condition")
makeBicorPlot(data.combined, "Genotype")

## heatmap to assess which components are actually important?
pdf(file = 'dimHeatmap_1_9.pdf')
DimHeatmap(object = data.combined, reduction.type = "cca", cells.use = 500, dim.use = 1:9, do.balanced = TRUE)
dev.off()

pdf(file = 'dimHeatmap_10_19.pdf')
DimHeatmap(object = data.combined, reduction.type = "cca", cells.use = 500, dim.use = 10:19, do.balanced = TRUE)
dev.off()

pdf(file = 'dimHeatmap_20_29.pdf')
DimHeatmap(object = data.combined, reduction.type = "cca", cells.use = 500, dim.use = 20:29, do.balanced = TRUE)
dev.off()

## choose number of CCs to use, based on above
chosen_cc <- 20

# Align CCA subspaces
print('~*~')
print('Align CCA subspaces...')
print(paste0('System time: ', Sys.time()))


# save variables
print('~*~')
print('Saving variables...')
print(paste0('System time: ', Sys.time()))

setwd(paste0(outputDir, subDir))
save.image(file = paste0("allVars.RData"))

print('~*~ All done! ~*~')