#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=interactive
#SBATCH --job-name=R_CMD
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=kmuench@stanford.edu
#SBATCH --mail-type=ALL

# Could do 64 GB but takes over a day

module purge
module load r/3.5.0

# Run R Script - Make individual variables for CCA
## Order of variables: (1)Folder containin .mtx data from 10x, (2) Metadata.csv, (3) output loc                                                                                           
#Rscript cca_makeVars.R $MTX_16p_SC_MS $METADAT_16p_SC_MS $OUTPUT_16p/201901_cluster_pooled_10x_ms/20190221_completeRunThrough/

# Merge the individual samples and perform CCA, visualize output to decide on number of CCs to use
## Order of variables: (1)Folder containing input from cca_makeVars, (2) output loc                                                                                                   
#Rscript cca_calcCCA.R $OUTPUT_16p/201901_cluster_pooled_10x_ms/20190221_completeRunThrough/ $OUTPUT_16p/201901_cluster_pooled_10x_ms/20190221_completeRunThrough/

# Align based on CCs, do tSNE, find clusters
## Order of variables: 
### (1) CCA-bearing variable set up in cca_calcCCA.R, 
### (2) output container folder,
### (3) chosen number of CCs
### (4) the chosen resolution
#Rscript cca_clusters.R $OUTPUT_16p/201901_cluster_pooled_10x_ms/20190219_makeSparse/calcCCA/data.combined_multiCCA_metadata.RData $OUTPUT_16p/201901_cluster_pooled_10x_ms/20190219_makeSparse/ 40 1.2
#Rscript cca_clusters.R $OUTPUT_16p/201901_cluster_pooled_10x_ms/20190221_completeRunThrough/calcCCA/data.combined_multiCCA_metadata.RData $OUTPUT_16p/201901_cluster_pooled_10x_ms/20190221_completeRunThrough/ 40 1.2

# Visualize output
## Order of variables: 
### (1) File containing CCA-bearing variable set up in cca_calcCCA.R, 
### (2) output container folder,
### (3) chosen number of CCs
### (4) the chosen resolution
Rscript cca_groupCompare.R $OUTPUT_16p/201901_cluster_pooled_10x_ms/20190219_makeSparse/clusters_simp/data.combined_withClust_r1.2_CC40.RData $OUTPUT_16p/201901_cluster_pooled_10x_ms/20190219_makeSparse/ 40 1.2

