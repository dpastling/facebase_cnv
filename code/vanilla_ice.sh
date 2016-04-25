#!/bin/bash
#BSUB -J R[3001-3739]%3
#BSUB -e logs/vanilla_ice_%J.log
#BSUB -o logs/vanilla_ice_%J.out
#BSUB -q normal
#BSUB -n 12
#BSUB -R "select[mem>40] rusage[mem=40] span[hosts=1]"
#BSUB -P Tamim

# 3739

DATA=data/Spritz_release_genotype_files
#RESULTS=results/vanilla_ice
RESULTS=data/corrected_vals
SAMPLES=($DATA/*_FinalReport_*.csv)
sample_file=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

Rscript code/vanilla_ice_single.R $sample_file $RESULTS

