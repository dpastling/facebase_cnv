#!/usr/bin/env bash
#BSUB -J R[1-14]
#BSUB -e logs/logR_%J.log
#BSUB -o logs/logR_%J.out
#BSUB -q normal
#BSUB -P Tamim

SAMPLES=(
"CIDR_06985@1054736171"
"CIDR_06985@1054736250"
"CIDR_06985@1064733235"
"CIDR_18502@1054736219"
"CIDR_18502@1054736256"
"CIDR_18505@1054736194"
"CIDR_18505@1054736220"
"CIDR_18505@1064692236"
"ADDL_18508@1064691934"
"ADDL_18508@1064692018"
"ADDL_18858@1064691581"
"ADDL_18858@1064691639"
"ADDL_18861@1064691574"
"ADDL_18861@1064691968"
)

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

Rscript code/get_logR.R $sample


