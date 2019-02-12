# Kristin Muench
# 2019.02.10
# Pipeline to import Seurat files generated earlier and identify clusters based on the number of PCs determined using runPCA.R


# install older version of Seurat which has the functions I need
library('Seurat', lib.loc='/scg/apps/software/r/3.5.0/scg/seurat_2.3')
packageVersion('Seurat')


# load needed files
print('Loading needed files...')
args <- commandArgs(TRUE) # load in arguments that accompanied script
mainDir <- args[1]
setwd(mainDir)
setwd('RData')
load('data_combined_pca.RData') # load combined data
load('metadataFiles.RData') # load combined data
load('hv_genes.RData') # load combined data
setwd('..')

data.combined <- data.combined.pca # rename this variable since it's the one we'll use from now on


# Create container folder for output of this script
subDir <- 'runClustering_pc20_r1p2'

if (file.exists(subDir)){
    setwd(file.path(mainDir, subDir))
} else {
    dir.create(file.path(mainDir, subDir))
    setwd(file.path(mainDir, subDir))

}


# Declare number of PCs to choose
pc_chosen <- 20
res_chosen <- 1.2

# Find Clusters of Cells
## save.SNN = T saves the SNN so that the clustering algorithm can be rerun
## using the same graph but with a different resolution value (see docs for
## full details)
## If this ends up taking too much time, lower n.start parameter to 10 (default 100) or employ an approximate nearest neighbor search via the RANN package by increasing the nn.eps parameter. Defafult is 0= exact nearest neighbor search
print('Finding clusters of cells...')

data.combined <- FindClusters(object = data.combined, reduction.type = "pca", dims.use = 1:pc_chosen, resolution = res_chosen, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

PrintFindClustersParams(object = data.combined) #print what settings you used



# Create a tSNE
print('Creating tSNE...')

getwd()
setwd(mainDir)
setwd(subDir)
data.combined <- RunTSNE(object = data.combined, dims.use = 1:pc_chosen, do.fast = TRUE, check_duplicates = FALSE)

## note that you can set do.label=T to help label individual clusters
pdf(file='tSNE_noLabels.pdf', height=8, width=8)
TSNEPlot(object = data.combined)
dev.off()

## save tSNE for later
setwd('../RData')
saveRDS(data.combined, file = paste0("data.combined_pc", pc_chosen, "_r", res_chosen, "r.rds") ) 

## build a dendrogram to visualize how close things are
setwd(mainDir)
setwd(subDir)
#pdf(file='buildClusterTree.pdf', height=8, width=8)
#data.combined <- BuildClusterTree(data.combined, pcs.use = pc_chosen)
#dev.off()

# Find cluster biomarkers
print('Finding clusters...')

library(dplyr) # so you can use magrittr's ceci n'est pas une pipe 
data.combined.markers <- FindAllMarkers(object = data.combined, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
data.combined.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

## print out clusters 
print('Printing out markers in clusters...')
nclust <- length(levels(data.combined.markers$cluster)) # total number of clusters

for (i in c(0:(length(levels(data.combined.markers$cluster)) -1 ) ) ) {
  print(i)
  clustMarkers <- data.combined.markers[data.combined.markers$cluster == i,]
  
  write.csv(clustMarkers, file = paste0('clustMarkers_',i,'.csv'), quote = FALSE, row.names=TRUE, col.names=TRUE )
}


## Plot heatmap of top ten genes per cluster
top10 <- data.combined.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
pdf(file='clusterMarkerHeatmap.pdf, height=15, width=15)
DoHeatmap(object = data.combined, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()

# save files
save.image(file = paste0("allVars_pc", pc_chosen, "_r", res_chosen, ".RData") )


print('~*~ Complete ~*~')