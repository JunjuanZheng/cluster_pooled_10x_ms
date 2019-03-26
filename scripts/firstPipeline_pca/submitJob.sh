#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --partition=interactive
#SBATCH --job-name=R_CMD
#SBATCH --time=7-00:00:00

module purge
module load r/3.5.0

# Run R Script
## Order of variables: (1)Folder containin .mtx data from 10x, (2)Output file, (3) Metadata .csv
#Rscript old/regular/runBigScaleData.R $MTX_16p_SC_MS $OUTPUT_16p/201901_cluster_pooled_10x_ms/standard/20190218 $METADAT_16p_SC_MS

## PCA script. Order of variables: (1) Output file
#Rscript old/regular/runPCA.R $OUTPUT_16p/201901_cluster_pooled_10x_ms/standard/20190218

## Clustering script. Order of variables: (1) Output file    
Rscript old/regular/runClustering.R $OUTPUT_16p/201901_cluster_pooled_10x_ms/standard/20190218
