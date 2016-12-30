#!/bin/bash
#BSUB -J R
#BSUB -e logs/vanilla_ice_%J_%I.log
#BSUB -o logs/vanilla_ice_%J_%I.out
#BSUB -q normal
#BSUB -n 12
#BSUB -R "select[mem>40] rusage[mem=40] span[hosts=1]"
#BSUB -P Tamim

# 3739

DATA=data/Spritz_release_genotype_files
RESULTS=results/vanilla_ice
#RESULTS=data/corrected_vals
#SAMPLES=($DATA/*_FinalReport_*.csv)
#sample_file=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

#SAMPLES=(`cat vice_missing.txt`)
#sample_file=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
#sample_file=$DATA/*_FinalReport_$sample_file.csv

sample_file=$DATA/*FinalReport_1141@0123840970.csv

$HOME/R-3.2.3/bin/Rscript code/vanilla_ice_single.R $sample_file $RESULTS

