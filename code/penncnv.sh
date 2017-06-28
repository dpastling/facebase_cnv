#!/bin/bash
#BSUB -J penncnv[1-1000]
#BSUB -e logs/penncnv_%J.log
#BSUB -o logs/penncnv_%J.out
#BSUB -q test
#BSUB -R "select[mem>10] rusage[mem=10]"
#BSUB -P Facebase

#3739

DATA=data/penncnv
RESULTS=results/penncnv
SAMPLES=($DATA/FinalReport_*)

sample_file=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
sample=`echo $sample_file | perl -pe 's/^.+?\/[^\/]+$/$1/'`
sample=`echo $sample | perl -pe 's/FinalReport_/$1/'`
sample=`echo $sample | perl -pe 's/\.txt$/$1/'`

detect_cnv.pl \
--test \
--minsnp 5 \
--confidence \
--hmm data/IlluminaHumanCoreExome_v12-A.hmm \
--pfb data/Human_Omni25exome.pfb \
--gcmodel data/Human_Omni25exome_GCModel.txt \
--log $RESULTS/sample_${sample}.log \
--out $RESULTS/sample_${sample}.rawcnv \
$sample_file

convert_cnv.pl --intype penncnv --outtype tab --output $RESULTS/sample_${sample}.tabcnv $RESULTS/sample_${sample}.rawcnv


