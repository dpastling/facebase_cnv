#!/bin/bash
#BSUB -J dnacopy[701-1000]
#BSUB -e logs/dnacopy_%J.log
#BSUB -o logs/dnacopy_%J.out
#BSUB -q test
#BSUB -R "select[mem>10] rusage[mem=10]"
#BSUB -P Tamim

DATA=data/penncnv
RESULTS=results/dnacopy
SAMPLES=($DATA/FinalReport_*)
sample_file=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

Rscript code/dnacopy_bySample.R $sample_file $RESULTS

