#!/usr/bin/env bash

SAMPLES=(
06985
18502
18505
18508
18858
18861
)

results=results/vanilla_ice_qc
valid=data/validation

for sample in ${SAMPLES[@]}
do
cp $valid/estd20_NA${sample}_5plus_probes.bed $results/${sample}_valid.bed

sample=$results/$sample

bedtools intersect -u -a ${sample}_all_vice_del.bed -b ${sample}_all_penn_del.bed > ${sample}_vice.penn_del.bed
bedtools intersect -u -a ${sample}_all_vice_dup.bed -b ${sample}_all_penn_dup.bed > ${sample}_vice.penn_dup.bed

bedtools intersect -u -a ${sample}_all_penn_del.bed -b ${sample}_all_vice_del.bed > ${sample}_penn.vice_del.bed
bedtools intersect -u -a ${sample}_all_penn_dup.bed -b ${sample}_all_vice_dup.bed > ${sample}_penn.vice_dup.bed

bedtools intersect -u -a ${sample}_all_vice_del.bed -b ${sample}_valid.bed > ${sample}_vice.valid_del.bed
bedtools intersect -u -a ${sample}_all_vice_dup.bed -b ${sample}_valid.bed > ${sample}_vice.valid_dup.bed

bedtools intersect -u -a ${sample}_all_penn_del.bed -b ${sample}_valid.bed > ${sample}_penn.valid_del.bed
bedtools intersect -u -a ${sample}_all_penn_dup.bed -b ${sample}_valid.bed > ${sample}_penn.valid_dup.bed

bedtools intersect -u -a ${sample}_all_dna_del.bed -b ${sample}_all_penn_del.bed > ${sample}_dna.penn_del.bed
bedtools intersect -u -a ${sample}_all_dna_dup.bed -b ${sample}_all_penn_dup.bed > ${sample}_dna.penn_dup.bed

bedtools intersect -u -a ${sample}_all_dna_del.bed -b ${sample}_all_vice_del.bed > ${sample}_dna.vice_del.bed
bedtools intersect -u -a ${sample}_all_dna_dup.bed -b ${sample}_all_vice_dup.bed > ${sample}_dna.vice_dup.bed

bedtools intersect -u -a ${sample}_all_dna_del.bed -b ${sample}_all_valid_del.bed > ${sample}_dna.valid_del.bed
bedtools intersect -u -a ${sample}_all_dna_dup.bed -b ${sample}_all_valid_dup.bed > ${sample}_dna.valid_dup.bed
done

