#!/bin/bash
#BSUB -J dnacopy[1-14]
#BSUB -e logs/dnacopy_%J.log
#BSUB -o logs/dnacopy_%J.out
#BSUB -q test
#BSUB -R "select[mem>10] rusage[mem=10]"
#BSUB -P Tamim

DATA=data/penncnv
RESULTS=results/dnacopy_test
SAMPLES=(
"18508@1064691934"
"18508@1064692018"
"18858@1064691581"
"18858@1064691639"
"18861@1064691574"
"18861@1064691968"
"06985@1054736171"
"06985@1054736250"
"06985@1064733235"
"18502@1054736219"
"18502@1054736256"
"18505@1054736194"
"18505@1054736220"
"18505@1064692236"
)
sample_file=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
sample_file=$DATA/FinalReport_*${sample_file}.txt

Rscript code/dnacopy_bySample2.R $sample_file $RESULTS

