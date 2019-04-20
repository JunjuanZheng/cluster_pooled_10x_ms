# Kristin Muench
# 2019.03.28
# Pipeline to run Seurat setup for mouse scRNA-Seq data - with CCA
# This takes in a matrix with CCAs calculated, adjusts them, makes a tSNE, 
# and visualizes/counts clusters.

# Load older version of Seurat which has the functions I need
library('Seurat', lib.loc='/scg/apps/software/r/3.5.0/scg/seurat_2.3')
packageVersion('Seurat')

## run this in case you accidentally load Seurat 3
#detach("package:Seurat", unload=TRUE)

# # Set up workspace - uncomment for running job
# print('Run CCA Track: groupCompare.R')
# print('Setting up workspace...')
# ##remove all environmental vars for clarity
# rm(list=ls())
# 
# ## import from command line
# args <- commandArgs(TRUE)
# inputFile <- args[1]
# outputDir <- args[2]
# chosen_cc <- args[3]
# chosen_res <- args[4]
# print(paste0('Importing files from: ', inputFile))
# print(paste0('Output location: ', outputDir))
# print(paste0('Chosen # of CCs: ', chosen_cc))
# print(paste0('Chosen resolution: ', chosen_res))

# Set up workspace -  Uncomment to work with files
inputFile <- "/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190329_troubleshootRunthrough_demux/clusters/data.combined_withClust_r1.2_CC40.RData"
outputDir <- "/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190329_troubleshootRunthrough_demux"
chosen_cc <- 40
chosen_res <- 1.2

## load needed files
print('~*~')
print('Loading needed files...')
load(inputFile)

## make subdirectory
setwd(outputDir)
subDir <- 'groupCompare'
print(paste0('Name of Experiment: ', subDir))

if (file.exists(subDir)){
  setwd(file.path(outputDir, subDir))
} else {
  dir.create(file.path(outputDir, subDir))
  setwd(file.path(outputDir, subDir))
}

print('Names of metadata')
print(names(data.combined@meta.data))


# Declare Stanford colors for prettier visualization
stanColors <- data.frame(cond = c('WT.SAL', 'WT.LPS', 'HET.SAL', 'HET.LPS'),
                         colors = as.character(c('#009b76', '#8c1515', '#b26f16', '#53284f')))



# Extract Cell IDs related to the clusters of interest
subclust_IDs <- data.combined@meta.data[which(data.combined@meta.data$res.1.2 %in% c('19', '13', '25')),]

## save 
setwd(file.path(outputDir, subDir))
save(subclust_IDs, file = 'subclust_IDs.RData')

# Visualization of things relevant to the Cluster 13/GABA story
FeatureHeatmap(data.combined, features.plot = c('Gad1', 'Gad2', 'Nkx2-1'), group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)
FeatureHeatmap(data.combined, features.plot = c('Ndn', 'Zfp503', 'Dlx2', 'Ebf1', 'Zeb2', 'Maz'), group.by = 'cond', pt.size = 0.25, key.position = "top", max.exp = 3)

## accompanying violin Plot
data.combined@meta.data$cond_clust <- paste0(data.combined@meta.data$cond, '_', data.combined@meta.data$res.1.2)

data.combined <- SetAllIdent(data.combined, "res.1.2")
data.combined_c13 <- SubsetData(data.combined, ident.use = '13')
VlnPlot(object = data.combined_c13, features.plot = c('Ebf1', 'Dlx2'), group.by = 'cond')

data.combined_c13 <- SetAllIdent(data.combined_c13, "cond")
data.combined_c13_justHet <- SubsetData(data.combined_c13, ident.use = c('WT.SAL', 'HET.SAL') )
VlnPlot(object = data.combined_c13_justHet, features.plot = c('Ebf1', 'Egr1', 'Zfp503', 'Ndn'), group.by = 'cond', nCol=2)
VlnPlot(object = data.combined_c13_justHet, features.plot = c('Maz', 'Aldoa', 'Dlx2', 'Ypel3'), group.by = 'cond', nCol=2)

## Dlx2 looks weird. What if we could cells with expression above threshold...?
data.combined_c13_normExprs <- FetchData(data.combined_c13, vars.all = c('Ebf1', 'Egr1', 'Zfp503', 'Ndn',
                                                                         'Maz', 'Aldoa', 'Dlx2', 'Ypel3'), 
                                                 use.raw=FALSE)
data.combined_c13_normExprs <- merge(data.combined_c13_normExprs, data.combined_c13@meta.data, by = 'row.names')

library(dplyr)
data.combined_c13_normExprs$Dlx2AboveThresh <- data.combined_c13_normExprs$Dlx2 > 0
tally_dlx <- as.tbl(data.frame(data.combined_c13_normExprs)) %>% group_by(Dlx2AboveThresh, cond) %>% tally()
mean_dlx <- as.tbl(data.frame(data.combined_c13_normExprs)) %>% group_by(cond) %>% dplyr::summarize(Mean = mean(Dlx2, na.rm=TRUE), SD = sd(Dlx2, na.rm=TRUE))


# How large are GABA interneuron populations relative to each other
data.combined_normExprs <- FetchData(data.combined)
data.combined_normExprs <- merge(data.combined_normExprs, data.combined@meta.data, by = 'row.names')
tally_gaba_popSize <- as.tbl(data.frame(data.combined_normExprs)) %>% group_by(res.1.2, sample, cond) %>% filter(res.1.2 %in% c(13,25,18)) %>% tally(.drop = FALSE)
totals <- tally_gaba_popSize %>% group_by(sample) %>% summarize(Total = sum(n), SD = sd(n))

tally_gaba_popSize$totals <- totals$Total[match(unlist(tally_gaba_popSize$sample), totals$sample)]
tally_gaba_popSize$percent_totals <- (tally_gaba_popSize$n/tally_gaba_popSize$totals) *100
tally_gaba_popSize$cond <- factor(tally_gaba_popSize$cond)
tally_gaba_popSize$res.1.2 <- factor(tally_gaba_popSize$res.1.2)

# prep for plotting
plotData <- tally_gaba_popSize %>% group_by(cond, res.1.2) %>% summarize(Mean = mean(percent_totals), SD = sd(percent_totals))
plotData$res.1.2 <- factor(plotData$res.1.2, levels = c('25', '13', '18'))
plotData$cond <- factor(plotData$cond, levels = c('WT.SAL', 'HET.SAL', 'WT.LPS', 'HET.LPS'))

p <- ggplot(plotData, aes(fill=res.1.2, y=Mean, x=cond)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), colour="black", width=.2, position=position_dodge(.9)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## can we do an anova?
### ref: https://www.statmethods.net/stats/anovaAssumptions.html
### ref: https://www.statisticssolutions.com/what-to-do-when-the-assumptions-of-your-analysis-are-violated/
### ref: https://www.graphpad.com/support/faq/what-to-do-when-data-fail-tests-for-homogeneity-of-variance/
### Q-Q plot:univariate normality? (significant departures = violate normality)
qqnorm(tally_gaba_popSize$percent_totals)
qqline(tally_gaba_popSize$percent_totals)
### homogeneity of variances
# Figner-Killeen Test of Homogeneity of Variances
## fails! so  the significance level will be underestimated, which can cause the null hypothesis to be falsely rejected
fligner.test(percent_totals ~ cond, data = tally_gaba_popSize)
fligner.test(percent_totals ~ res.1.2, data = tally_gaba_popSize)

## two-way anova with factors for condition and cluster
inhibAnova <- aov(percent_totals~cond*res.1.2,data=
                    tally_gaba_popSize)
summary(inhibAnova)
posthoc <- TukeyHSD(x=inhibAnova)

               # save variables
               print('~*~')
               #print('Saving variables...')
               #print(paste0('System time: ', Sys.time()))
               
               #setwd(paste0(outputDir, subDir))
               #save.image(file = paste0("allVars.RData"))
               
               print('~*~ All done! ~*~')
               print(paste0('System time: ', Sys.time()))