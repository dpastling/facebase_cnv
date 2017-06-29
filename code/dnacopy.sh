#!/bin/bash
#BSUB -J dnacopy[1-3]%5
#BSUB -e logs/dnacopy_%J.log
#BSUB -o logs/dnacopy_%J.out
#BSUB -q normal
#BSUB -R "select[mem>10] rusage[mem=10]"
#BSUB -P Facebase

# TODO: need to check if input to DNAcopy is ok
# DATA=data/penncnv

source code/config_example.sh

sample_file=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

Rscript code/dnacopy_bySample.R $sample_file $DNACOPY_RESULTS

