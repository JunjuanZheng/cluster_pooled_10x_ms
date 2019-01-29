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
Rscript runBigScaleData.R $MTX_16p_SC_MS $OUTPUT_16p/201901_cluster_pooled_10x_ms

#source("20190121_runBigScaleData.R", echo=TRUE, max.deparse.length=10000)
