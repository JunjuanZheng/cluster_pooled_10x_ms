#!/bin/bash
#
# submitJob.sh: Job submission script for CCA-based scRNA-Seq pipeline
# By Kristin Muench
# 2019.03.26
#
# SLURM submission flags:
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=124G
#SBATCH --partition=interactive
#SBATCH --job-name=troubsh_demux
#SBATCH --time=6-00:00:00
#SBATCH --mail-user=kmuench@stanford.edu
#SBATCH --mail-type=ALL
#
# Notes on choosing parameters:
# Got out of memory error using 64G on calcCCA.
#
# Example run:
# qsub submitJob.sh
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# prepare modules
module purge
module load r/3.5.0

# load in paths
barcodeSampleLUT=/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190326_troubleshootRunthrough/determineCellSex/allNewID_cells.csv # could be demultiplexed IDs or not
##, e.g. /scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190219_makeSparse/20190226_sexLabel/allNewID_cells.csv ## original, flawed
## /labs/tpalmer/projects/cnv16p/data/scRNASeq/mouse/metadata/originalSampleBarcodeLUT.csv ## no attempt to demultiplex
outputDir="$OUTPUT_16p/201901_cluster_pooled_10x_ms/20190329_troubleshootRunthrough_demux" # directory where output stored

# declare for documentation purposes what these variables are
#echo 'Barcode Sample Lookup Table Location: ' $barcodeSampleLUT
#echo 'Output Directory: ' $outputDir

# Make individual Seurat objects
## Order of variables: (1)Folder containin .mtx data from 10x, (2) Metadata.csv, (3) output loc   (4) barcodeSampleLUT     
#echo 'Making variables...'
#Rscript makeVars.R $MTX_16p_SC_MS $METADAT_16p_SC_MS $outputDir $barcodeSampleLUT

# Merge the individual samples and perform CCA, visualize output to decide on number of CCs to use
## Order of variables: (1)Folder containing input from makeVars.R, (2) output loc   
#outputDir_makeVars="$outputDir/makeVars/allVars.RData"
#echo 'Calculating CCA...'
#echo 'Importing this makeVars.RData file:' $outputDir_makeVars
#Rscript calcCCA.R $outputDir_makeVars $outputDir

# Align based on CCs, do tSNE, find clusters
## Order of variables: 
### (1) CCA-bearing variable set up in cca_calcCCA.R, 
### (2) output container folder,
### (3) chosen number of CCs
### (4) the chosen resolution
#echo 'Calculating Clusters...'
#Rscript clusters.R "$outputDir/calcCCA/data.combined_multiCCA_metadata.RData" $outputDir 40 1.2

# Visualize output
## Order of variables: 
### (1) File containing CCA-bearing variable set up in cca_calcCCA.R, 
### (2) output container folder,
### (3) chosen number of CCs
### (4) the chosen resolution
echo 'Making useful visualizations...'
outputDir_clusters="$outputDir/clusters/data.combined_withClust_r1.2_CC40.RData"
echo 'Importing this RData file with data and clusters:'
Rscript groupCompare.R $outputDir_clusters $outputDir 40 1.2
