#!/usr/bin/env bash
#BSUB -J convert[1-1000]
#BSUB -e logs/convert_%J.log
#BSUB -o logs/convert_%J.out
#BSUB -q normal
#BSUB -P Tamim

IMAGE=plots/dnacopy
#SAMPLES=($IMAGE/by_chr_*.pdf)
SAMPLES=($IMAGE/genome_wide_*.pdf)
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

new=`echo $sample | sed 's/pdf$/png/'`
convert -density 150 -depth 8 -flatten $sample $new

