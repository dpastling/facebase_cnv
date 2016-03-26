#!/bin/bash
#BSUB -J R[1-55]
#BSUB -e logs/vanilla_ice_%J.log
#BSUB -o logs/vanilla_ice_%J.out
#BSUB -q normal
#BSUB -n 12
#BSUB -R "select[mem>40] rusage[mem=40] span[hosts=1]"
#BSUB -P Tamim

# 3739

DATA=data/Spritz_release_genotype_files
RESULTS=results/vanilla_ice
#SAMPLES=($DATA/*_FinalReport_*.csv)
source code/missing_files.sh
sample_file=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

Rscript code/vanilla_ice_single.R $sample_file $RESULTS

