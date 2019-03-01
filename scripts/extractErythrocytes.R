# Kristin Muench
# 2019.02.26
# Trying to extract erythrocytes

library('Seurat', lib.loc='/scg/apps/software/r/3.5.0/scg/seurat_2.3')

# load data
load('/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190219_makeSparse/clusters_simp/data.combined_withClust_r1.2_CC40.RData')

# set output dir
setwd('/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190219_makeSparse/20190226_sexLabel')

# rerun pipeline to access raw reads
# data.combined@raw.data

# cell names that are erythrocytes
cells_eryth <- row.names(data.combined@meta.data[data.combined@meta.data$res.1.2 == '24',])
cells_mesen <- row.names(data.combined@meta.data[data.combined@meta.data$res.1.2 == c('16','18'),])
cells_vasc <- row.names(data.combined@meta.data[data.combined@meta.data$res.1.2 == c('31','26',
                                                                                     '27','28'),])
cells_int <- row.names(data.combined@meta.data[data.combined@meta.data$res.1.2 == c('13','17',
                                                                                    '20','22',
                                                                                    '23', '25'),])

cells_RGC <- row.names(data.combined@meta.data[data.combined@meta.data$res.1.2 == c('1','2',
                                                                                    '3','4'),])


# cells from samples with known Sex
data.combined <- SetAllIdent(data.combined, id = "sample")
maleSamps <- c('KM2', 'KM13', 'KM1', 'KM18', 'J')
femSamps <- c('I', 'K', 'L')
defM <- SubsetData(data.combined, ident.use = maleSamps,
                   subset.raw = T)
defF <- SubsetData(data.combined, ident.use = femSamps,
                   subset.raw = T)


# plot Xist and other sex genes for Erythrocytes vs. non-Erythrocytes

## create labels
### Erythrocytes
data.combined@meta.data$isEryth <- 'Not Erythrocytes'
data.combined@meta.data[row.names(data.combined@meta.data) %in% cells_eryth, 'isEryth'] <- 'Erythrocytes'

### Interneurons
data.combined@meta.data$isInt <- 'Not Interneurons'
data.combined@meta.data[row.names(data.combined@meta.data) %in% cells_int, 'isInt'] <- 'Interneurons'

### Interneurons
data.combined@meta.data$isRGC <- 'Not RGC'
data.combined@meta.data[row.names(data.combined@meta.data) %in% cells_RGC, 'isRGC'] <- 'RGC'


## focus on these labels

## function
plotVlnPlot <- function(data.combined, plotGene, myGroup, myTitle, chosen_res, chosen_cc){
  data.combined <- SetAllIdent(data.combined, id = myGroup)
  
  # plot
  p <- VlnPlot(object = data.combined, features.plot = plotGene, group.by = myGroup, 
               do.return = TRUE)
  #p_raw <- VlnPlot(object = data.combined, features.plot = plotGene, group.by = myGroup, do.return = TRUE, use.raw=TRUE)
  
  #p_scaled <- VlnPlot(object = data.combined, features.plot = plotGene, group.by = myGroup,  do.return = TRUE, use.scaled=TRUE)
  
  # save
  ggsave(file = paste0(myTitle, '_', myGroup, '_', plotGene, '_r', chosen_res, 
                       '_cc', chosen_cc, '.pdf'), plot = p, device='pdf')
  
  #ggsave(file = paste0(myTitle, '_raw_', myGroup, '_', plotGene, '_r', chosen_res, '_cc', chosen_cc, '.pdf'), plot = p_raw, device='pdf')
  
  #ggsave(file = paste0(myTitle, '_scaled_', myGroup, '_', plotGene, '_r', chosen_res, '_cc', chosen_cc, '.pdf'), plot = p_scaled, device='pdf')
  
}

## deploy
plotVlnPlot(data.combined, 'Xist', 'isEryth', 'violinPlot', 1.2, 40)
# plotVlnPlot(data.combined, 'Ddx3y', 'isEryth', 'violinPlot', 1.2, 40)
# plotVlnPlot(data.combined, 'Eif2s3y', 'isEryth', 'violinPlot', 1.2, 40)

plotVlnPlot(data.combined, 'Xist', 'isInt', 'violinPlot', 1.2, 40)
plotVlnPlot(data.combined, 'Xist', 'isRGC', 'violinPlot', 1.2, 40)

plotVlnPlot(defF, 'Xist', 'isRGC', 'violinPlot_F', 1.2, 40)
plotVlnPlot(defM, 'Xist', 'isRGC', 'violinPlot_M', 1.2, 40)
plotVlnPlot(defF, 'Xist', 'isEryth', 'violinPlot_F', 1.2, 40)
plotVlnPlot(defM, 'Xist', 'isEryth', 'violinPlot_M', 1.2, 40)
plotVlnPlot(defF, 'Xist', 'isInt', 'violinPlot_F', 1.2, 40)
plotVlnPlot(defM, 'Xist', 'isInt', 'violinPlot_M', 1.2, 40)


# What % of cells have Xist=0, 0<x<1, 1<x<1.5, 1.5<x<2, in M and F, and for diff cell types?
## function for quick readout
countBins <- function(seuObj){
  print(paste0('Number of cells = 0: ', length(which(seuObj == 0)) ))
  print(paste0('Number of cells bewteen 0 and 1: ', length(which(seuObj > 0 & seuObj < 1)) ))
  print(paste0('Number of cells bewteen 1 and 1.5: ', length(which(seuObj > 1 & seuObj < 1.5)) ))
  print(paste0('Number of cells bewteen 1.5 and 2: ', length(which(seuObj > 1.5 & seuObj < 2)) ))
  print(paste0('Number of cells above 2: ', length(which(seuObj > 2)) ))
}


## for log2/normalized data
### extract count data
defF_Xist_all <- FetchData(defF, vars.all = 'Xist')
defF_Xist_eryth <- FetchData(defF, vars.all = 'Xist', cells.use = cells_eryth[cells_eryth %in% row.names(defF@meta.data)] )
defF_Xist_vasc <- FetchData(defF, vars.all = 'Xist', cells.use = cells_vasc[cells_vasc %in% row.names(defF@meta.data)] )
defF_Xist_RGC <- FetchData(defF, vars.all = 'Xist', cells.use = cells_RGC[cells_RGC %in% row.names(defF@meta.data)] )
defF_Xist_int <- FetchData(defF, vars.all = 'Xist', cells.use = cells_int[cells_int %in% row.names(defF@meta.data)] )

defM_Xist_all <- FetchData(defM, vars.all = 'Xist')
defM_Xist_eryth <- FetchData(defM, vars.all = 'Xist', cells.use = cells_eryth[cells_eryth %in% row.names(defM@meta.data)] )
defM_Xist_vasc <- FetchData(defM, vars.all = 'Xist', cells.use = cells_vasc[cells_vasc %in% row.names(defM@meta.data)] )
defM_Xist_RGC <- FetchData(defM, vars.all = 'Xist', cells.use = cells_RGC[cells_RGC %in% row.names(defM@meta.data)] )
defM_Xist_int <- FetchData(defM, vars.all = 'Xist', cells.use = cells_int[cells_int %in% row.names(defM@meta.data)] )

### use function to examine cells in each bin
countBins(defF_Xist_all)
countBins(defF_Xist_eryth)
countBins(defF_Xist_vasc)
countBins(defF_Xist_RGC)
countBins(defF_Xist_int)

countBins(defM_Xist_all)
countBins(defM_Xist_eryth)
countBins(defM_Xist_vasc)
countBins(defM_Xist_RGC)
countBins(defM_Xist_int)

## for raw
defM_Xist_all_raw <- FetchData(defM, vars.all = 'Xist', use.raw=TRUE)
defF_Xist_all_raw <- FetchData(defF, vars.all = 'Xist', use.raw=TRUE)
countBins(defF_Xist_all_raw)
countBins(defM_Xist_all_raw)



# What are the cell identities of F cells with low Xist, and M cells with high Xist?
thresh <- 1

## pull out what cells fit that criteria
F_LowXist <- row.names(data.frame(defF_Xist_all[which(defF_Xist_all < thresh),]) )
F_LowXist_raw <- row.names(data.frame(defF_Xist_all_raw[which(defF_Xist_all_raw < thresh),]) )
M_highXist <- row.names(data.frame(defM_Xist_all[which(defM_Xist_all > thresh),]) )
M_lowXist <- row.names(data.frame(defM_Xist_all[which(defM_Xist_all < thresh),]) )


F_LowXistMeta <- data.combined@meta.data[which(row.names(data.combined@meta.data) %in% F_LowXist ) , ]
F_LowXistMeta_raw <- data.combined@meta.data[which(row.names(data.combined@meta.data) %in% F_LowXist_raw ) , ]
M_HighXistMeta <- data.combined@meta.data[which(row.names(data.combined@meta.data) %in% M_highXist ) , ]

## tally
library(dplyr)
count_F_LowXist_Clust <- data.frame(F_LowXistMeta %>% count(res.1.2))
ggplot(count_F_LowXist_Clust, aes(x="", y=n, fill=res.1.2)) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0)

count_M_HighXist_Clust <- data.frame(M_HighXistMeta %>% count(res.1.2))
ggplot(count_M_HighXist_Clust, aes(x="", y=n, fill=res.1.2)) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0)

count_F_LowXist_Clust_raw <- data.frame( F_LowXistMeta_raw %>% count(res.1.2) )
ggplot(count_F_LowXist_Clust_raw, aes(x="", y=n, fill=res.1.2)) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0)

# Of the cells that are F-lowXist and M-highXist, what's their male gene expression like?
exprs_F_LowXist_Ddx3y <- FetchData(defF, vars.all = 'Ddx3y', cells.use = F_LowXist)
hist(exprs_F_LowXist_Ddx3y)
exprs_F_LowXist_Eif2s3y <- FetchData(defF, vars.all = 'Eif2s3y', cells.use = F_LowXist)
hist(exprs_F_LowXist_Eif2s3y)

exprs_M_lowXist_Ddx3y <- FetchData(defM, vars.all = 'Ddx3y', cells.use = M_lowXist)
exprs_M_lowXist_Eif2s3y <- FetchData(defM, vars.all = 'Eif2s3y', cells.use = M_lowXist)
exprs_M_lowXist_Ddx3y_raw <- FetchData(defM, vars.all = 'Ddx3y', cells.use = M_lowXist, use.raw = TRUE)
exprs_M_lowXist_Eif2s3y_raw <- FetchData(defM, vars.all = 'Eif2s3y', cells.use = M_lowXist, use.raw = TRUE)
exprs_M_lowXist_manGenes_raw <- data.frame(FetchData(defM, vars.all = c('Eif2s3y', 'Ddx3y'), cells.use = M_lowXist, use.raw = TRUE))


## try it again
allLowXist <- c(F_LowXist, M_lowXist)
lowXistCells <- SubsetData(data.combined, cells.use=allLowXist, subset.raw = T)
lowXistCells@meta.data$sex <- 'F'
lowXistCells@meta.data[ which(row.names(lowXistCells@meta.data) %in% M_lowXist) , 'sex'] <- 'M'

plotVlnPlot(lowXistCells, 'Xist', 'sex', 'violinPlot_lowXistCells', 1.2, 40)
plotVlnPlot(lowXistCells, 'Eif2s3y', 'sex', 'violinPlot_lowXistCells', 1.2, 40)
plotVlnPlot(lowXistCells, 'Ddx3y', 'sex', 'violinPlot_lowXistCells', 1.2, 40)

## how many of the low Xist cells express > 0 Eif2s3y? Ddx3y?
exprs_M_lowXist_Eif2s3y_raw <- FetchData(defM, vars.all = 'Eif2s3y', cells.use = M_lowXist, use.raw = TRUE)

## how many have one OR the other?
nrow(exprs_M_lowXist_manGenes_raw[which(exprs_M_lowXist_manGenes_raw$Ddx3y > 0 | exprs_M_lowXist_manGenes_raw$Eif2s3y > 0),])
nrow(exprs_M_lowXist_manGenes_raw[which(exprs_M_lowXist_manGenes_raw$Ddx3y == 0 & exprs_M_lowXist_manGenes_raw$Eif2s3y == 0),])


# Provide sex labels to the cells
exprs_all_sexGenes_raw <- data.frame(FetchData(data.combined, vars.all = c('Xist', 'Eif2s3y', 'Ddx3y'), use.raw = TRUE))
cells_probF <- row.names(exprs_all_sexGenes_raw[which(exprs_all_sexGenes_raw$Xist > thresh),])

## adjust metadata
data.combined@meta.data$sex <- 'M'
data.combined@meta.data[ which(row.names(data.combined@meta.data) %in% cells_probF) ,'sex'] <- 'F'

### quick check - for known maleSamps and femSamps, how'd we do?
#### helpful: https://www.r-graph-gallery.com/48-grouped-barplot-with-ggplot2/
tally_check_accuracy_Xist <- data.frame(data.combined@meta.data %>% group_by(sample) %>% count(sex))
tally_check_accuracy_Xist <- rbind(tally_check_accuracy_Xist, c('J', 'F', 0), c('KM1', 'F', 0), c('KM13', 'F', 0), c('KM18', 'F', 0), c('KM2', 'F', 0))
tally_check_accuracy_Xist$n <- as.numeric(tally_check_accuracy_Xist$n)

ggplot(tally_check_accuracy_Xist, aes(fill=sex, y=n, x=sample)) + 
  geom_bar( stat="identity")
ggplot(tally_check_accuracy_Xist, aes(fill=sex, y=n, x=sample)) + 
  geom_bar( stat="identity", position = "fill")

## adjust indivdual samples
data.combined@meta.data[ which( data.combined@meta.data$sample %in% femSamps) ,'sex' ] <- 'F'
data.combined@meta.data[ which( data.combined@meta.data$sample %in% maleSamps) ,'sex'] <- 'M'

tally_check_accuracy_Xist_afterCorrection <- data.combined@meta.data %>% group_by(sample) %>% count(sex)
ggplot(tally_check_accuracy_Xist_afterCorrection, aes(fill=sex, y=n, x=sample)) + 
  geom_bar( stat="identity")
ggplot(tally_check_accuracy_Xist_afterCorrection, aes(fill=sex, y=n, x=sample)) + 
  geom_bar( stat="identity", position = "fill")

## For samples F and G, what proportion of those samples express Male Sex Genes?
### If half or more express Male Sex Genes = more likely they are all male?
justFandG <- SubsetData(data.combined, ident.use = c('F','G', 'KM1', 'H'),
                        subset.raw = T)
cells_sampF <- row.names(justFandG@meta.data[ which(justFandG@meta.data$sample == 'F') ,])
cells_sampG <-  row.names(justFandG@meta.data[ which(justFandG@meta.data$sample == 'G') ,])
cells_sampKM1 <-  row.names(justFandG@meta.data[ which(justFandG@meta.data$sample == 'KM1') ,])
cells_sampH <-  row.names(justFandG@meta.data[ which(justFandG@meta.data$sample == 'H') ,])

exprs_justFandG_sexGenes <- data.frame(FetchData(justFandG, vars.all = c('Eif2s3y', 'Ddx3y', 'Xist'), use.raw = TRUE))
exprs_sampF <- exprs_justFandG_sexGenes[which(row.names(exprs_justFandG_sexGenes) %in% cells_sampF),]
exprs_sampG <- exprs_justFandG_sexGenes[which(row.names(exprs_justFandG_sexGenes) %in% cells_sampG),]
exprs_sampKM1 <- exprs_justFandG_sexGenes[which(row.names(exprs_justFandG_sexGenes) %in% cells_sampKM1),]
exprs_sampH <- exprs_justFandG_sexGenes[which(row.names(exprs_justFandG_sexGenes) %in% cells_sampH),]

# function to print out counts
howManyInEachSexGeneCategory <- function(exprs_samp, thresh){
  
  totalCells <- nrow(exprs_samp)
  totalNotXist <- nrow(exprs_samp[ which(exprs_samp$Xist < thresh) ,])
  totalXist <- nrow(exprs_samp[ which(exprs_samp$Xist > thresh) ,])
  totalEif <- nrow(exprs_samp[ which(exprs_samp$Eif2s3y > 0) ,])
  totalEif_NotXist <- nrow(exprs_samp[ which(exprs_samp$Eif2s3y > 0 &
                                               exprs_samp$Xist < thresh) ,])
  totalEif_Xist <- nrow(exprs_samp[ which(exprs_samp$Eif2s3y > 0 &
                                            exprs_samp$Xist > thresh) ,])
  
  totalDdx <- nrow(exprs_samp[ which(exprs_samp$Ddx3y > 0) ,])
  totalDdx_NotXist <- nrow(exprs_samp[ which(exprs_samp$Ddx3y > 0 &
                                               exprs_samp$Xist < thresh) ,])
  totalDdx_Xist <- nrow(exprs_samp[ which(exprs_samp$Ddx3y > 0 &
                                            exprs_samp$Xist > thresh) ,])
  
  print(paste0( 'How many total cells: ', totalCells))
  print(paste0( 'How many Xist-neg cells: ', totalNotXist,
                'In percent: ', (totalNotXist/totalCells)*100, '%'))
  print(paste0( 'Xist Neg, Eif2s3y+: ', totalEif_NotXist, 
                '-- Percent of Xist Neg that are Eif+:', (totalEif_NotXist/totalNotXist)*100, '%',
                '-- Percent of Xist Pos that are Eif+:', (totalEif_Xist/totalXist)*100, '%' ))
  
  print(paste0( 'Xist Neg, Ddx3y+: ', totalDdx_NotXist, 
                '-- Percent of Xist Neg that are Ddx+:', (totalDdx_NotXist/totalNotXist)*100, '%',
                '-- Percent of Xist Pos that are Ddx+:', (totalDdx_Xist/totalXist)*100, '%'))
  
}

howManyInEachSexGeneCategory(exprs_sampF, thresh)
howManyInEachSexGeneCategory(exprs_sampG, thresh)
howManyInEachSexGeneCategory(exprs_sampKM1, thresh)
howManyInEachSexGeneCategory(exprs_sampH, thresh)



# calculating specificity and sensitivity
A_trueM_calledM <- sum(tally_check_accuracy_Xist[which(tally_check_accuracy_Xist$sample %in% maleSamps &
                                                         tally_check_accuracy_Xist$sex == 'M'), 'n'])
B_trueM_calledF <- sum(tally_check_accuracy_Xist[which(tally_check_accuracy_Xist$sample %in% maleSamps &
                                                         tally_check_accuracy_Xist$sex == 'F'), 'n'])
C_trueF_calledM <- sum(tally_check_accuracy_Xist[which(tally_check_accuracy_Xist$sample %in% femSamps &
                                                         tally_check_accuracy_Xist$sex == 'M' ), 'n'])
D_trueF_calledF <- sum(tally_check_accuracy_Xist[which(tally_check_accuracy_Xist$sample %in% femSamps &
                                                         tally_check_accuracy_Xist$sex == 'F' ), 'n'])

sensitivity <- A_trueM_calledM/(A_trueM_calledM+C_trueF_calledM) * 100
specificity <- D_trueF_calledF/(D_trueF_calledF+B_trueM_calledF) * 100


# PLOT ROC CURVE
## On ROC curves: https://www.medcalc.org/manual/roc-curves.php

allLabeledSamps <- SubsetData(data.combined, ident.use = c(maleSamps, femSamps), subset.raw = T)
allLabeledSamps_exprs <- data.frame(FetchData(allLabeledSamps, vars.all = c('Eif2s3y', 'Ddx3y', 'Xist'), use.raw = TRUE))

test_roc <- data.frame(allLabeledSamps@meta.data[,'sex'])
row.names(test_roc) <- row.names(allLabeledSamps@meta.data)
test_roc <- merge(test_roc, allLabeledSamps_exprs, by = 'row.names')
colnames(test_roc) <- c('cellName', 'sex', 'Eif2s3y', 'Ddx3y', 'Xist')
write.csv(test_roc, file='test_roc_ExprsAndTrueSexLabel.csv')



# Save cell labels
metadata_withSex <- data.frame(data.combined@meta.data)
metadata_withSex$newID <- paste0(metadata_withSex$sample, '_', metadata_withSex$sex)
metadata_withSex$newID <- as.factor(metadata_withSex$newID)

## save all the barcodes
cells <- data.frame(Barcode = sapply(strsplit(row.names(metadata_withSex ), split="_"),'[', 2) ,
                    sample = metadata_withSex$newID)
setwd('/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190219_makeSparse/20190226_sexLabel')
write.csv(cells, file='allNewID_cells.csv', quote=FALSE, row.names=FALSE)

## save specific sample barcodes
extractBarcodesAndSave <- function(metadata_withSex, sampID){
  # extract barcodes
  cells <- data.frame(Barcode = sapply(strsplit(row.names(metadata_withSex[ which(metadata_withSex$newID == sampID) , ] ), split="_"),'[', 2))
  cells$Genotype <- sampID
  colnames(cells) <- c('Barcode', 'Genotype')
  
  # save csv
  write.csv(cells, file=paste0(sampID,'_cells.csv'), quote=FALSE, row.names=FALSE)
  
  # return
  return(cells)
}

cells_E_F <- extractBarcodesAndSave(metadata_withSex,'E_F')
cells_E_M <- extractBarcodesAndSave(metadata_withSex,'E_M')

# save variables
setwd('/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190219_makeSparse/20190226_sexLabel')
save.image(file = paste0("allVars.RData"))

