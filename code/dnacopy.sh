#!/bin/bash
#BSUB -J dnacopy[2001-3000]%5
#BSUB -e logs/dnacopy_%J.log
#BSUB -o logs/dnacopy_%J.out
#BSUB -q normal
#BSUB -R "select[mem>10] rusage[mem=10]"
#BSUB -P Tamim

#3739

DATA=data/penncnv
RESULTS=results/dnacopy_alpha0_1
SAMPLES=($DATA/FinalReport_*)
sample_file=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

Rscript code/dnacopy_bySample.R $sample_file $RESULTS

