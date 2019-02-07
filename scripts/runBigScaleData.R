# Kristin Muench
# 2019.01.21
# Pipeline to run Seurat setup for mouse scRNA-Seq data

# Load needed library
#install.packages('Seurat')


# Load older version of Seurat which has the functions I need
library('Seurat', lib.loc='/scg/apps/software/r/3.5.0/scg/seurat_2.3')
packageVersion('Seurat')


# info about this session
#print(sessionInfo())

# import from command line
args <- commandArgs(TRUE)
#load(args[1])
print(paste0('Output location: ', args[2]))

# get file names
list.filenames <- list.files(args[1])
identifier <- 'E' # make arbitrary identifier to initialize

print('Example file location:')
print(paste0(args[1], '/', list.filenames[1],'/mm10/'))

# import Seurat Objects
print('Import Seurat Objects...')
for (i in 1:length(list.filenames)){
  identifier[i] <- strsplit(list.filenames[i], "_")[[1]][[1]]
  print(identifier[i])
  print(identifier)
  
  # import data
  data <- Read10X(data.dir = paste0(args[1], '/', list.filenames[i], '/mm10/') )
  seuObj <- CreateSeuratObject(raw.data = data, project = identifier[i], min.cells = 3, min.genes = 200)

  # replace name with identifier
  assign(paste0(identifier[i],'_data'), data)
  assign(paste0(identifier[i],'_obj'), seuObj)
}

# combine data
print('Combining data...')
data.combined <- MergeSeurat(object1 = eval(parse(text = paste0(identifier[1], '_obj' ) ) ), 
                             object2 = eval(parse(text = paste0(identifier[2], '_obj' ) ) ),
                             add.cell.id1 = paste0(identifier[1], '_obj' ), 
                             add.cell.id2 = paste0(identifier[2], '_obj' ), project = "all")

## notice the cell names now have an added identifier
head(x = data.combined@cell.names)
table(data.combined@meta.data$orig.ident)

# use mergeSeurat iteratively to build data.combined
print('Adding samples...')
for (i in 3:length(identifier)) {
  data.combined <- MergeSeurat(data.combined,
                       eval(parse(text = paste0(identifier[i], '_obj' ) ) ),
                       add.cell.id2 = paste0(identifier[i], '_obj' ),
                       project = "all")
}

## notice the cell names now have an added identifier
head(x = data.combined@cell.names)
table(data.combined@meta.data$orig.ident)

data.combined.noNorm.noScale <- data.combined # so you can save this later as an intermediate variable



# Import metadata using cell names from data.combined
print('Importing metadata...')
metadata <- read.csv(args[3])

myCells <- data.frame(cellNames = data.combined@cell.names)
myCells$Identity <- data.combined@meta.data$orig.ident
metadataCols <- c('SampleID','CellsPerSample','SurgeryDate','Condition', 'Genotype', 'Litter')

myMetadata <- merge(myCells, metadata[,metadataCols], by.x='Identity', by.y='SampleID')
rownames(myMetadata) <- myMetadata$cellNames

## add my metadata to Seurat object
data.combined <- AddMetaData(data.combined, myMetadata[,-c(1:2)], col.name = metadataCols[-c(1:2)])



# Normalize Data
print('Normalizing data...')
data.combined <- NormalizeData(object = data.combined, normalization.method = "LogNormalize", 
    scale.factor = 10000)

data.combined.normOnly <- data.combined # so you can save intermediate variable


# Detect variable genes across cells
print('Identifying genes with highest variability...')
## export plot
setwd(args[2])
svg( filename = "findVariableGenes.svg" )
data.combined <- FindVariableGenes(object = data.combined, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
dev.off()

hv.genes <- head(rownames(data.combined@hvg.info), 1000) # only take top 1000 variable genes


# Scale Data
print('Scaling data...')
data.combined <- ScaleData(object = data.combined, genes.use = hv.genes, vars.to.regress = c("nUMI"), do.par = TRUE)

data.combined.normAndScale <- data.combined


# export output to file
print('Saving data...')
setwd(args[2])
save(data.combined.noNorm.noScale, file= 'data_combined_noNorm_noScale.RData')
save(data.combined.normOnly, file='data_combined_normOnly.RData')
save(data.combined.normAndScale, file='data_combined_normAndScale.RData')
save(metadata, myMetadata, file='metadataFiles.RData')
save(hv.genes, file='hv_genes.RData')

# Notify user that script is finished
print('~*~All finished!~*~')