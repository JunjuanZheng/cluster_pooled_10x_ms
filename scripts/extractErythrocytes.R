# Kristin Muench
# 2019.02.26
# Trying to extract erythrocytes



# load data
#load('/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190219_makeSparse/clusters_simp/data.combined_withClust_r1.2_CC40.RData')

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
defM <- SubsetData(data.combined, ident.use = c('KM2', 'KM13', 'KM1', 'KM18', 'J'),
                      subset.raw = T)
defF <- SubsetData(data.combined, ident.use = c('I', 'K', 'L'),
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
data.combined <- SetAllIdent(data.combined, id = "isEryth")

## function
plotVlnPlot <- function(data.combined, plotGene, myGroup, myTitle, chosen_res, chosen_cc){
  # plot
  p <- VlnPlot(object = data.combined, features.plot = plotGene, group.by = myGroup, 
               do.return = TRUE)
  p_raw <- VlnPlot(object = data.combined, features.plot = plotGene, group.by = myGroup, 
                   do.return = TRUE, use.raw=TRUE)
  
  p_scaled <- VlnPlot(object = data.combined, features.plot = plotGene, group.by = myGroup, 
                   do.return = TRUE, use.scaled=TRUE)
  
  # save
  ggsave(file = paste0(myTitle, '_', myGroup, '_', plotGene, '_r', chosen_res, 
                       '_cc', chosen_cc, '.pdf'), plot = p, device='pdf')
  
  ggsave(file = paste0(myTitle, '_raw_', myGroup, '_', plotGene, '_r', chosen_res, 
                       '_cc', chosen_cc, '.pdf'), plot = p_raw, device='pdf')
  
  ggsave(file = paste0(myTitle, '_scaled_', myGroup, '_', plotGene, '_r', chosen_res, 
                       '_cc', chosen_cc, '.pdf'), plot = p_scaled, device='pdf')
  
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



countBins <- function(seuObj){
  print(paste0('Number of cells = 0: ', length(which(seuObj == 0)) ))
  print(paste0('Number of cells bewteen 0 and 1: ', length(which(seuObj > 0 & seuObj < 1)) ))
  print(paste0('Number of cells bewteen 1 and 1.5: ', length(which(seuObj > 1 & seuObj < 1.5)) ))
  print(paste0('Number of cells bewteen 1.5 and 2: ', length(which(seuObj > 1.5 & seuObj < 2)) ))
  print(paste0('Number of cells above 2: ', length(which(seuObj > 2)) ))
}

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


# What are the cell identities of F cells with low Xist, and M cells with high Xist?

# Repeat the above with raw data and see if that filters better?

# Create hypothetical filter for the male genes

# Are cells with high Xist also high in the other two?

# ident.remove the erythrocytes?

# for everything else, threshold and assing?

## If Xist expression is higher than 1? F
## Of remaining cells, is Y expression higher than thresh? M
## For remainder - what cell types are they and randomly re-permute labels

