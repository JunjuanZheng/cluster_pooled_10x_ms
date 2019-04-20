# monocle.R
# Kristin Muench
# 2019.03.26
# Pipeline to run monocle
#
# Tutorial used: http://cole-trapnell-lab.github.io/monocle-release/docs/#constructing-single-cell-trajectories
# # # # # # # # # # # # # # # # # #

# Load needed libraries
library(monocle)
library(reshape2)

# make output dir
outputDir <- "/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190408_monocle"

# load needed paths
load('/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190329_troubleshootRunthrough_demux/makeVars/seuratObj_wt.sal.RData')
load('/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190329_troubleshootRunthrough_demux/clusters/data.combined_withClust_r1.2_CC40.RData')
# load('/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190408_monocle/allVars_20190408.RData')

# make cell ID to cluster map
mapIDtoClust <- data.frame(ID = row.names(data.combined@meta.data), Cluster = data.combined@meta.data$res.1.2)
wt.sal@meta.data$res.1.2 <- mapIDtoClust[match(row.names(wt.sal@meta.data), mapIDtoClust$ID), 'Cluster']

# make sparse to ease computational burden?
### Note: 10x .mtx files are already sparse matricies - don't want to make it a dense matrix!
#wt.sal <- MakeSparse(object = wt.sal)
#data.combined <- MakeSparse(object = data.combined)

# create CDS files for Monocle
## !!! Already has a negativeBinomial.size Expression Family
minimal_cds.wt.sal <- importCDS(wt.sal, import_all = TRUE)
cds.wt.sal <- importCDS(wt.sal)

# estimate size factors and dispersions
cds.wt.sal <- estimateSizeFactors(cds.wt.sal)
cds.wt.sal<- estimateDispersions(cds.wt.sal)

# Filtering low-quality cells
## how many genes detected per cell?
cds.wt.sal <- detectGenes(cds.wt.sal, min_expr = 0.1)
print(head(fData(cds.wt.sal)))
expressed_genes <- row.names(subset(fData(cds.wt.sal), num_cells_expressed >= 10))


## filter cells out - use judgment when picking parameters
## made up parameters - not sure what's best in this case
## inspiration: https://github.com/satijalab/seurat/issues/259
valid_cells <- row.names(subset(pData(cds.wt.sal),
                                nGene > 500 &
                                  nGene < 8000))
cds.wt.sal <- cds.wt.sal[,valid_cells]

## what is the distribution of mRNA across all cells?
pData(cds.wt.sal)$Total_mRNAs <- Matrix::colSums(exprs(cds.wt.sal))

cds.wt.sal <- cds.wt.sal[,pData(cds.wt.sal)$Total_mRNAs < 1e6]

upper_bound <- 10^(mean(log10(pData(cds.wt.sal)$Total_mRNAs)) + 
                     2*sd(log10(pData(cds.wt.sal)$Total_mRNAs))) # more reads than typical per cell; attempt to exclude doublets
lower_bound <- 10^(mean(log10(pData(cds.wt.sal)$Total_mRNAs)) - 
                     2*sd(log10(pData(cds.wt.sal)$Total_mRNAs))) # fewer reads than typical per cell

qplot(Total_mRNAs, data = pData(cds.wt.sal), color = cond, geom = "density") +
  geom_vline(xintercept = lower_bound) + geom_vline(xintercept = upper_bound)

cds.wt.sal <- cds.wt.sal[,pData(cds.wt.sal)$Total_mRNAs > lower_bound & pData(cds.wt.sal)$Total_mRNAs < upper_bound]
cds.wt.sal <- detectGenes(cds.wt.sal, min_expr = 0.1)

## Verify they follow a distribution that is roughly log-normal
### Log-transform each value in the expression matrix.
L <- log(exprs(cds.wt.sal[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.
qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")


# # Classify and count cells using markers  
# 
# ##Imputing cell type
# 
# EOMES_id <- row.names(subset(fData(cds.wt.sal), gene_short_name == "EOMES"))
# PAX6_id <- row.names(subset(fData(cds.wt.sal), gene_short_name == "PAX6"))
# DCX_id <- row.names(subset(fData(cds.wt.sal), gene_short_name == "DCX"))
# 
# cth <- newCellTypeHierarchy()
# cth <- addCellType(cth, "IntermediateProgenitor", classify_func = function(x) { x[EOMES_id,] >= 1 })
# cth <- addCellType(cth, "RadialGlia", classify_func = function(x) 
# { x[EOMES_id,] < 1 & x[DCX_id,] < 1 & x[PAX6_id,] > 1 })
# cth <- addCellType(cth, "YoungNeuron", classify_func = function(x) 
# { x[DCX_id,] > 1 })
# 
# ## add back into cds
# cds.wt.sal <- classifyCells(cds.wt.sal, cth, 0.1)
# 
# ## make pie chart showing how many there are
# table(pData(cds.wt.sal)$CellType)
# 
# pie <- ggplot(pData(cds.wt.sal), aes(x = factor(1), fill = factor(CellType))) +
#   geom_bar(width = 1)
# pie + coord_polar(theta = "y") +
#   theme(axis.title.x = element_blank(), axis.title.y = element_blank())
# 
# ## see tutorial for clustering
# 
# 
# 
# 
# # Classify and count cells without marker cells
# 
# ## pick useful subset of genes to use for classification
# disp_table <- dispersionTable(cds.wt.sal)
# unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
# cds.wt.sal <- setOrderingFilter(cds.wt.sal, unsup_clustering_genes$gene_id) ## genes that will be used for clustering
# plot_ordering_genes(cds.wt.sal)
# 
# ## cluster - look for batches
# # HSMM@auxClusteringData[["tSNE"]]$variance_explained <- NULL ## do I need this? not sure
# plot_pc_variance_explained(cds.wt.sal, return_all = F) # PCS on which clustering will be performed
# cds.wt.sal <- reduceDimension(cds.wt.sal, max_components = 2, num_dim = 6,
#                         reduction_method = 'tSNE', verbose = T) # is there a way to use cond here?
# cds.wt.sal <- clusterCells(cds.wt.sal, num_clusters = 2)
# plot_cell_clusters(cds.wt.sal, 1, 2, color = "cond", markers = c("SOX2", "SST")) #!!! why only 2392 samples in assayData?? Not 17197?
# plot_cell_clusters(cds.wt.sal, 1, 2, color = "sample")
# plot_cell_clusters(cds.wt.sal, 1, 2, color = "nUMI")
# 
# ## uncomment this to remove unwanted sources of variation, e.g. Litter
# # HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 2,
# #                         reduction_method = 'tSNE',
# #                         residualModelFormulaStr = "~Media + num_genes_expressed",
# #                         verbose = T)
# # HSMM <- clusterCells(HSMM, num_clusters = 2)
# # plot_cell_clusters(HSMM, 1, 2, color = "CellType")
# #cds.wt.sal <- clusterCells(cds.wt.sal, num_clusters = 2)
# cds.wt.sal <- clusterCells(cds.wt.sal)
# plot_cell_clusters(cds.wt.sal, 1, 2, color = "cond") + facet_wrap(~orig.ident)
# 
# ##!!! Another step where they choose subset of cells to cluster
# 


# Create subsets of only the cells we thing have related lineages
## RGC series: 0, 2, 4, 5, 12
## IP series: 6, 7, 10, 14, 20
## Neuron series: 1, 3, 11, 9, 21
## Excitatory neuron series: 0, 2, 4, 5, 12, 6, 7, 10, 14, 20, 1, 3, 11, 9, 21
## ZLI-adjacent series: 8, 19, 23, 26
## GABA series: 13, 18, 25
## C-R series: 15, 16
# !!!!Ω
cells_RGC <- row.names(subset(pData(cds.wt.sal),
                              row.names(pData(cds.wt.sal)) %in% mapIDtoClust[mapIDtoClust$Cluster=='0', 'ID'] |
                                row.names(pData(cds.wt.sal)) %in% mapIDtoClust[mapIDtoClust$Cluster=='2', 'ID'] |
                                row.names(pData(cds.wt.sal)) %in% mapIDtoClust[mapIDtoClust$Cluster=='4', 'ID'] |
                                row.names(pData(cds.wt.sal)) %in% mapIDtoClust[mapIDtoClust$Cluster=='5', 'ID'] |
                                row.names(pData(cds.wt.sal)) %in% mapIDtoClust[mapIDtoClust$Cluster=='12', 'ID']
                              ) )

cells_IP <- row.names(subset(pData(cds.wt.sal),
                              row.names(pData(cds.wt.sal)) %in% mapIDtoClust[mapIDtoClust$Cluster=='6', 'ID'] |
                                row.names(pData(cds.wt.sal)) %in% mapIDtoClust[mapIDtoClust$Cluster=='7', 'ID'] |
                                row.names(pData(cds.wt.sal)) %in% mapIDtoClust[mapIDtoClust$Cluster=='10', 'ID'] |
                                row.names(pData(cds.wt.sal)) %in% mapIDtoClust[mapIDtoClust$Cluster=='14', 'ID'] |
                                row.names(pData(cds.wt.sal)) %in% mapIDtoClust[mapIDtoClust$Cluster=='20', 'ID']
) )

cells_Neuron <- row.names(subset(pData(cds.wt.sal),
                             row.names(pData(cds.wt.sal)) %in% mapIDtoClust[mapIDtoClust$Cluster=='1', 'ID'] |
                               row.names(pData(cds.wt.sal)) %in% mapIDtoClust[mapIDtoClust$Cluster=='3', 'ID'] |
                               row.names(pData(cds.wt.sal)) %in% mapIDtoClust[mapIDtoClust$Cluster=='11', 'ID'] |
                               row.names(pData(cds.wt.sal)) %in% mapIDtoClust[mapIDtoClust$Cluster=='9', 'ID'] |
                               row.names(pData(cds.wt.sal)) %in% mapIDtoClust[mapIDtoClust$Cluster=='21', 'ID']
) )

cells_excite <- c(cells_RGC, cells_IP, cells_Neuron)

cells_ZLI <- row.names(subset(pData(cds.wt.sal),
                                 row.names(pData(cds.wt.sal)) %in% mapIDtoClust[mapIDtoClust$Cluster=='8', 'ID'] |
                                   row.names(pData(cds.wt.sal)) %in% mapIDtoClust[mapIDtoClust$Cluster=='19', 'ID'] |
                                   row.names(pData(cds.wt.sal)) %in% mapIDtoClust[mapIDtoClust$Cluster=='23', 'ID'] |
                                   row.names(pData(cds.wt.sal)) %in% mapIDtoClust[mapIDtoClust$Cluster=='26', 'ID']
) )

cells_C_R <- row.names(subset(pData(cds.wt.sal),
                                row.names(pData(cds.wt.sal)) %in% mapIDtoClust[mapIDtoClust$Cluster=='15', 'ID'] |
                                  row.names(pData(cds.wt.sal)) %in% mapIDtoClust[mapIDtoClust$Cluster=='16', 'ID']
                                  ) )


cds.wt.sal.C_R <- cds.wt.sal[,cells_C_R]
cds.wt.sal.RGC <- cds.wt.sal[,cells_RGC]

# Pseudotime trajectories
## Step 1. identify "interesting" genes to use
###DE genes from those at beginning and those at end
diff_test_res <- differentialGeneTest(cds.wt.sal.C_R,
                                      fullModelFormulaStr = "~res.1.2")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))

### set filter now that you've picked genes to use for ordering
cds.wt.sal.C_R <- setOrderingFilter(cds.wt.sal.C_R, ordering_genes)
plot_ordering_genes(cds.wt.sal.C_R)

## Step 2. Reduce dimensionality !!!!!!!!!
cds.wt.sal.C_R <- reduceDimension(cds.wt.sal.C_R, max_components = 2, method = 'DDRTree')
#cds.wt.sal.C_R <- reduceDimension(cds.wt.sal.C_R, max_components = 2, method = 'tSNE')

## Step 3. Visualize trajectory in reduced dimensional space
cds.wt.sal.C_R <- orderCells(cds.wt.sal.C_R, color.by='cond')



# Save variables

# save variables
print('~*~')
#print('Saving variables...')
#print(paste0('System time: ', Sys.time()))

setwd(outputDir)
save.image(file = "allVars_20190409.RData")

print('~*~ All done! ~*~')
print(paste0('System time: ', Sys.time()))
