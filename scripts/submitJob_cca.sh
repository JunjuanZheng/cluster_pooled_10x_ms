#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=interactive
#SBATCH --job-name=R_CMD
#SBATCH --time=7-00:00:00

module purge
module load r/3.5.0

# Run R Script - CCA-based pipeline
## Order of variables: (1)Folder containin .mtx data from 10x, (2) Metadata.csv, (3) output loc                                                                                           
Rscript cca_makeVars.R $MTX_16p_SC_MS $METADAT_16p_SC_MS $OUTPUT_16p/201901_cluster_pooled_10x_ms/20190214/

## Order of variables: (1)Folder containin .mtx data from 10x, (2) Metadata.csv, (3) output loc                                                                                                   
Rscript cca_ccaAndTSNE.R $OUTPUT_16p/201901_cluster_pooled_10x_ms/20190214/

