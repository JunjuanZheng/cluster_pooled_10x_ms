sampl# Kristin Muench
# 2019.02.12
# Pipeline to run Seurat setup for mouse scRNA-Seq data - with CCA

# Load needed library
#install.packages('Seurat')

# Load older version of Seurat which has the functions I need
library('Seurat', lib.loc='/scg/apps/software/r/3.5.0/scg/seurat_2.3')
packageVersion('Seurat')

# Load cell IDs for demultiplexing samples
allNewID_cells <- read.csv("/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190219_makeSparse/20190226_sexLabel/allNewID_cells.csv")

# import from command line
args <- commandArgs(TRUE)
#load(args[1])
mtxPath <- args[1]
metadataPath <- args[2]
outputDir <- args[3]

# # !!! FOR TROUBLESHOOTING - UNCOMMENT FOR PATHS
# mtxPath <- '/labs/tpalmer/projects/cnv16p/data/scRNASeq/mouse/mtxData'
# metadataPath <- '/labs/tpalmer/projects/cnv16p/data/scRNASeq/mouse/metadata/20190211_sampleTable.csv'
# outputDir <- '/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190228_runThroughDemultiplex/'

print('Run CCA Track: cca_makeVars.R')
print(paste0('Mtx File location: ', mtxPath))
print(paste0('Metadata location: ', metadataPath))
print(paste0('Output location: ', outputDir))

# make subdirectory
setwd(outputDir)
subDir <- 'makeVars'

if (file.exists(subDir)){
  setwd(file.path(outputDir, subDir))
} else {
  dir.create(file.path(outputDir, subDir))
  setwd(file.path(outputDir, subDir))
}

# get file names
list.filenames <- list.files(mtxPath)

print('Example file location:')
print(paste0(mtxPath, '/', list.filenames[1],'/mm10/'))

# read metadata
metadata<-read.csv(metadataPath)

# in list filenames, what are the corresponding sample IDs?
files <- data.frame(filenames = list.filenames)
files$ids <- unlist(lapply(strsplit(list.filenames, "_"), '[[', 1))

# pull those sample IDs into a new list.filenames
list.filenames.wt.sal <- files[files$ids %in% metadata[metadata$Group == 'WT.SAL','SampleID'], 'filenames']
list.filenames.wt.lps <- files[files$ids %in% metadata[metadata$Group == 'WT.LPS','SampleID'], 'filenames']
list.filenames.het.sal <- files[files$ids %in% metadata[metadata$Group == 'HET.SAL','SampleID'], 'filenames']
list.filenames.het.lps <- files[files$ids %in% metadata[metadata$Group == 'HET.LPS','SampleID'], 'filenames']



# loop, but where first entry is used to make mergeSeurat, and subsequent used to tack on Seurat
# condName, e.g. 'WT.SAL'
# files: variable created above that has both file IDs and sample IDs
# !!!!! WORK ON MAKING ORIG.IDENT ACCURATE. FIX FILES HERE WHICH REFERES TO EVERYTHING RN
myCondSeurat <- function(filenames, mtxPath, condName, files){
  print('~*~')
  print(paste0('Working on ', condName))
  
  # pull out only files of interest
  files_relevant <- files[which(files$filenames %in% filenames),]
  print(files_relevant)
  
  # make first seurat object
  print(paste0('Making initial Seurat object with ', filenames[1], '...'))
  initialData <- Read10X(data.dir = paste0(mtxPath, '/', filenames[1], '/mm10/') )
  data <- CreateSeuratObject(raw.data = initialData, project = condName, min.cells = 3, min.features=200, names.field = files_relevant[1,2])
  #data <- MakeSparse(object = data) # make it sparse !!!!! NEW ADDITION FOR THIS RUN 2.19
  
  # # Are there any cells that didn't get categorized? Default to calling them 'Unknown sample' (think filtered out later)
  # # (otherwise find and troubleshoot using which(!(row.names(data@meta.data) %in% allNewID_cells$Barcode)))
  # data@meta.data$tst <- 'Unknown'
  # data@meta.data[ which(row.names(data@meta.data) %in% allNewID_cells$Barcode), 'tst'] <- as.character(allNewID_cells[ which(row.names(data@meta.data) %in% allNewID_cells$Barcode) , 'sample'])
  # Add sample IDs
  data@meta.data$orig.ident <- files_relevant[1,'ids'] 
  
  # now merge with the others
  for (i in 2:length(filenames)){
    
    # make this iteration's seurat object
    addingData <- Read10X(data.dir = paste0(mtxPath, '/', filenames[i], '/mm10/') )
    print(paste0('Reading in data from ',mtxPath, '/', filenames[i],'...'))
    addingData_obj <- CreateSeuratObject(raw.data = addingData, project = condName, min.cells = 3, min.features=200) # create object
    addingData_obj <- MakeSparse(object = addingData_obj) # make it sparse !!! NEW ADDITION FOR THIS RUN 2.19
    addingData_obj@meta.data$orig.ident <- files_relevant[i,'ids']
    
    # merge this new seurat object with existing data
    print(paste0('Merging Seurat data with data from ', filenames[i], '...'))
    data <- MergeSeurat(data,addingData_obj, project = condName, add.cell.id2 = filenames[i])
  }
  
  # give everything in this Seurat Object a condition name
  data@meta.data$stim <- condName
  
  # prefilter now that all the samples are in
  data <- FilterCells(data, subset.names = "nGene", low.thresholds = 200, high.thresholds = Inf)
  data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
  data <- ScaleData(data, display.progress = F, vars.to.regress = "nUMI", do.par = TRUE, num.cores = 4) # took out genes.use = hv.genes - seems complicated with separate seurat obj
  
  # give everything an alternative ID
  data@meta.data$trueBarcode <- sapply(strsplit(row.names(data@meta.data ), split="_"),tail, n=1L)
  data@meta.data$demuxSample <- allNewID_cells[ match(data@meta.data$trueBarcode,allNewID_cells$Barcode ) , 'sample']
  
  return(data)
}


wt.sal <- myCondSeurat(list.filenames.wt.sal, mtxPath, 'WT.SAL', files)
wt.lps <- myCondSeurat(list.filenames.wt.lps, mtxPath, 'WT.LPS', files)
het.sal <- myCondSeurat(list.filenames.het.sal, mtxPath, 'HET.SAL', files)
het.lps <- myCondSeurat(list.filenames.het.lps, mtxPath, 'HET.LPS', files)

# Print what is in each one
table(wt.sal@meta.data$orig.ident)
table(wt.lps@meta.data$orig.ident)
table(het.sal@meta.data$orig.ident)
table(het.lps@meta.data$orig.ident)

# save variables
setwd(paste0(outputDir, subDir))
save.image(file = paste0("allVars.RData"))

print('~*~ All done! ~*~')