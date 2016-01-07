#!/bin/bash
#BSUB -J penncnv[3001-3739]
#BSUB -e logs/penncnv_%J.log
#BSUB -o logs/penncnv_%J.out
#BSUB -q test
#BSUB -R "select[mem>10] rusage[mem=10]"
#BSUB -P Tamim

#3739

DATA=data/penncnv
RESULTS=results/penncnv
#SAMPLES=($DATA/FinalReport_*)
SAMPLES=(
NA18508
NA18858
NA18861
NA06985
NA18502
NA18505
)

sample_file=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
sample=`echo $sample_file | perl -pe 's/^.+?FinalReport_(.+?)\.txt$/$1/'`

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

#detect_cnv.pl \
#--test \
#--minsnp 5 \
#--minlength 50000 \
#--confidence \
#--hmm data/IlluminaHumanCoreExome_v12-A.hmm \
#--pfb data/Human_Omni25exome.pfb \
#--gcmodel data/Human_Omni25exome_GCModel.txt \
#--log $RESULTS/sample_${sample}.log \
#--out $RESULTS/sample_${sample}.rawcnv \
#$sample_file
#filter_cnv.pl --numsnp 5 --length 50k --output sampleall.filtered.cnv sampleall.rawcnv

convert_cnv.pl --intype penncnv --outtype tab --output $RESULTS/sample_${sample}.tabcnv $RESULTS/sample_${sample}.rawcnv


