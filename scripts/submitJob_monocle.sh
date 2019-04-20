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
#SBATCH --mem=64G
#SBATCH --partition=interactive
#SBATCH --job-name=monocl
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
module load miniconda/3

# Visualize output
echo 'Running monocle:'
Rscript monocle.R
