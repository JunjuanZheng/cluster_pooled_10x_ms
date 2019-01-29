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

# notice the cell names now have an added identifier
head(x = data.combined@cell.names)
table(data.combined@meta.data$orig.ident)

# Append that data to existing seurat object
# # ## M
# 
# data.combined <- AddSamples(object = data.combined, 
#                               new.data = I_obj, 
#                               add.cell.id = 'I' ) # works if it's data

# add more samples
print('Adding samples...')
for (i in 1:length(identifier)) {
  data.combined <- AddSamples(object = data.combined, 
                              new.data = eval(parse (text = paste0(identifier[i], '_data' ) ) ), 
                              add.cell.id = paste0(identifier[i], '_data' ) )
}

# notice the cell names now have an added identifier
head(x = data.combined@cell.names)
table(data.combined@meta.data$orig.ident)

data.combined.noNorm.noScale <- data.combined

# Normalize Data
print('Normalizing data...')
data.combined <- NormalizeData(object = data.combined, normalization.method = "LogNormalize", 
    scale.factor = 10000)

data.combined.normOnly <- data.combined # so you can save intermediate variable

# Scale Data
print('Scaling data...')
data.combined <- ScaleData(object = data.combined, vars.to.regress = c("nUMI"))

data.combined.normAndScale <- data.combined

# export output to file
print('Saving data...')
setwd(args[2])
save(data.combined.noNorm.noScale, file= 'data_combined_noNorm_noScale.RData')
save(data.combined.normOnly, file='data_combined_normOnly.RData')
save(data.combined.normAndScale, file='data_combined_normAndScale.RData')

# Notify user that script is finished
print('~*~All finished!~*~')