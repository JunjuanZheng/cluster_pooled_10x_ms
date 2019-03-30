# monocle.R
# Kristin Muench
# 2019.03.26
# Pipeline to run monocle
# # # # # # # # # # # # # # # # # #

# Load needed libraries
library(monocle)
library(reshape2)

# load needed paths
load('/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190329_troubleshootRunthrough_demux/makeVars/seuratObj_wt.sal.RData')
load('/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190228_runThroughDemultiplex/clusters/data.combined_withTSNE_r1.2_CC40.RData')

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


# Classify and count cells using markers  

##Imputing cell type

EOMES_id <- row.names(subset(fData(cds.wt.sal), gene_short_name == "EOMES"))
PAX6_id <- row.names(subset(fData(cds.wt.sal), gene_short_name == "PAX6"))
DCX_id <- row.names(subset(fData(cds.wt.sal), gene_short_name == "DCX"))

cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "IntermediateProgenitor", classify_func = function(x) { x[EOMES_id,] >= 1 })
cth <- addCellType(cth, "RadialGlia", classify_func = function(x) 
{ x[EOMES_id,] < 1 & x[DCX_id,] < 1 & x[PAX6_id,] > 1 })
cth <- addCellType(cth, "YoungNeuron", classify_func = function(x) 
{ x[DCX_id,] > 1 })

## add back into cds
cds.wt.sal <- classifyCells(cds.wt.sal, cth, 0.1)

## make pie chart showing how many there are
table(pData(cds.wt.sal)$CellType)

pie <- ggplot(pData(cds.wt.sal), aes(x = factor(1), fill = factor(CellType))) +
  geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

## see tutorial for clustering




# Classify and count cells without marker cells

## pick useful subset of genes to use for classification
disp_table <- dispersionTable(cds.wt.sal)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds.wt.sal <- setOrderingFilter(cds.wt.sal, unsup_clustering_genes$gene_id) ## genes that will be used for clustering
plot_ordering_genes(cds.wt.sal)

## cluster - look for batches
# HSMM@auxClusteringData[["tSNE"]]$variance_explained <- NULL ## do I need this? not sure
plot_pc_variance_explained(cds.wt.sal, return_all = F) # PCS on which clustering will be performed
cds.wt.sal <- reduceDimension(cds.wt.sal, max_components = 2, num_dim = 6,
                        reduction_method = 'tSNE', verbose = T) # is there a way to use cond here?
cds.wt.sal <- clusterCells(cds.wt.sal, num_clusters = 2)
plot_cell_clusters(cds.wt.sal, 1, 2, color = "cond", markers = c("SOX2", "SST")) #!!! why only 2392 samples in assayData?? Not 17197?
plot_cell_clusters(cds.wt.sal, 1, 2, color = "sample")
plot_cell_clusters(cds.wt.sal, 1, 2, color = "nUMI")

## uncomment this to remove unwanted sources of variation, e.g. Litter
# HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 2,
#                         reduction_method = 'tSNE',
#                         residualModelFormulaStr = "~Media + num_genes_expressed",
#                         verbose = T)
# HSMM <- clusterCells(HSMM, num_clusters = 2)
# plot_cell_clusters(HSMM, 1, 2, color = "CellType")
#cds.wt.sal <- clusterCells(cds.wt.sal, num_clusters = 2)
cds.wt.sal <- clusterCells(cds.wt.sal)
plot_cell_clusters(cds.wt.sal, 1, 2, color = "cond") + facet_wrap(~orig.ident)

##!!! Another step where they choose subset of cells to cluster



# Pseudotime trajectories


