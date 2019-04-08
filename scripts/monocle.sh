# monocle.R
# Kristin Muench
# 2019.03.26
# Pipeline to run monocle
# # # # # # # # # # # # # # # # # #

# Load needed library
library('monocle')

# paths to files of interest
path_data.combined <- '/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190228_runThroughDemultiplex/clusters/data.combined_withTSNE_r1.2_CC40.RData'
path_wt.sal <- '/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190228_runThroughDemultiplex/makeVars/'

load(path_wt.sal)
load(path_data.combined)

## make sparse to ease computational burden?
### Note: 10x .mtx files are already sparse matricies - don't want to make it a dense matrix!
#wt.sal <- MakeSparse(object = wt.sal)
#data.combined <- MakeSparse(object = data.combined)

# create CDS files for Monocle
## !!! does this already have an expression family?
minimal_cds.wt.sal <- importCDS(data_to_be_imported, import_all = TRUE)
cds.wt.sal <- importCDS(data_to_be_importedt)

## specify their distribution
#HSMM <- newCellDataSet(count_matrix, phenoData = pd, featureData = fd, expressionFamily=negbinomial.size())

# estimate size factors and dispersions
cds.wt.sal <- estimateSizeFactors(cds.wt.sal)
cds.wt.sal<- estimateDispersions(cds.wt.sal)

# Filtering low-quality cells
## how many genes detected per cell?
cds.wt.sal <- detectGenes(cds.wt.sal, min_expr = 0.1)
print(head(fData(cds.wt.sal)))

## filter cells out - use judgment when picking parameters
valid_cells <- row.names(subset(pData(cds.wt.sal),
                    Cells.in.Well == 1 &
                    Control == FALSE &
                    Clump == FALSE &
                    Debris == FALSE &
                    Mapped.Fragments > 1000000))
        cds.wt.sal <- cds.wt.sal[,valid_cells]
        
## what is the distribution of mRNA across all cells?
pData(cds.wt.sal)$Total_mRNAs <- Matrix::colSums(exprs(cds.wt.sal))

cds.wt.sal <- cds.wt.sal[,pData(cds.wt.sal)$Total_mRNAs < 1e6]

upper_bound <- 10^(mean(log10(pData(cds.wt.sal)$Total_mRNAs)) + 2*sd(log10(pData(cds.wt.sal)$Total_mRNAs))) # more reads than typical per cell; attempt to exclude doublets
lower_bound <- 10^(mean(log10(pData(cds.wt.sal)$Total_mRNAs)) - 2*sd(log10(pData(cds.wt.sal)$Total_mRNAs))) # fewer reads than typical per cell

qplot(Total_mRNAs, data = pData(cds.wt.sal), color = Hours, geom = "density") +
geom_vline(xintercept = lower_bound) + geom_vline(xintercept = upper_bound)

lower_bound <- lower_bound[,pData(lower_bound)$Total_mRNAs > lower_bound & pData(lower_bound)$Total_mRNAs < upper_bound]
lower_bound <- detectGenes(lower_bound, min_expr = 0.1)

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


# Classify and count cells






