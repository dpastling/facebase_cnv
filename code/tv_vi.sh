#!/bin/bash
#BSUB -J ice_ice_baby[1-14]
#BSUB -e logs/ice_%J.log
#BSUB -o logs/ice_%J.out
#BSUB -q normal
#BSUB -R "select[mem>10] rusage[mem=10]"
#BSUB -P Tamim

DATA=data/Spritz_release_genotype_files
RESULTS=results/vanilla_ice
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
sample_file=$DATA/*FinalReport*${sample_file}*.csv

code/tv_vi.R $sample_file $RESULTS

