#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=124G
#SBATCH --partition=interactive
#SBATCH --job-name=R_CMD
#SBATCH --time=2-00:00:00

# Could do 64 GB but takes over a day

module purge
module load r/3.5.0

# Run R Script - CCA-based pipeline
## Order of variables: (1)Folder containin .mtx data from 10x, (2) Metadata.csv, (3) output loc                                                                                           
#Rscript cca_makeVars.R $MTX_16p_SC_MS $METADAT_16p_SC_MS $OUTPUT_16p/201901_cluster_pooled_10x_ms/20190214/

## Order of variables: (1)Folder containing input from cca_makeVars, (2) output loc                                                                                                   
#Rscript cca_calcCCA.R $OUTPUT_16p/201901_cluster_pooled_10x_ms/20190214/ $OUTPUT_16p/201901_cluster_pooled_10x_ms/20190218/

## Order of variables: 
### (1) CCA-bearing variable set up in cca_calcCCA.R, 
### (2) output container folder,
### (3) chosen number of CCs
### (4) the chosen resolution
#Rscript cca_clusters.R $OUTPUT_16p/201901_cluster_pooled_10x_ms/20190214/ccaAndTSNE_sampleCol_actually0216 $OUTPUT_16p/201901_cluster_pooled_10x_ms/20190218/ 26 1.2
Rscript cca_clusters_simp.R $OUTPUT_16p/201901_cluster_pooled_10x_ms/20190214/ccaAndTSNE_sampleCol_actually0216/data.combined_multiCCA_metadata.RData $OUTPUT_16p/201901_cluster_pooled_10x_ms/20190219/ 26 1.2

## Order of variables: 
### (1) File containing CCA-bearing variable set up in cca_calcCCA.R, 
### (2) output container folder,
### (3) chosen number of CCs
### (4) the chosen resolution
#Rscript cca_clusters.R $OUTPUT_16p/201901_cluster_pooled_10x_ms/20190218/clusters_longTroubleshoot/data.combined_withTSNE_r1.2_CC26.RData $OUTPUT_16p/201901_cluster_pooled_10x_ms/20190218/ 26 1.2

