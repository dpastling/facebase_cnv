#!/bin/bash
#BSUB -J penncnv[1-3]
#BSUB -e logs/penncnv_%J.log
#BSUB -o logs/penncnv_%J.out
#BSUB -q test
#BSUB -R "select[mem>10] rusage[mem=10]"
#BSUB -P Facebase

#3739

source code/config_example.sh

sample_file=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
sample=`echo $sample_file | perl -pe 's/^.+?\/[^\/]+$/$1/'`
sample=`echo $sample | perl -pe 's/FinalReport_/$1/'`
sample=`echo $sample | perl -pe 's/\.txt$/$1/'`

detect_cnv.pl \
--test \
--minsnp 5 \
--confidence \
--hmm $hmm_file \
--pfb $pfb_file \
--gcmodel $gc_model \
--log $RESULTS/sample_${sample}.log \
--out $RESULTS/sample_${sample}.rawcnv \
$sample_file

convert_cnv.pl --intype penncnv --outtype tab --output $RESULTS/sample_${sample}.tabcnv $RESULTS/sample_${sample}.rawcnv


