#!/usr/bin/env bash

SAMPLES=(
18508.1064691934
18508.1064692018
18858.1064691581
18858.1064691639
18861.1064691574
18861.1064691968
06985.1054736171
06985.1054736250
06985.1064733235
18502.1054736219
18502.1054736256
18505.1054736194
18505.1054736220
18505.1064692236
)

valid_folder=data/validation
penn_folder=$valid_folder
dna_folder=results/dnacopy_test
results_folder=results/dnacopy_test
#dna_folder=$valid_folder
#results_folder=$valid_folder

for sample in ${SAMPLES[@]}
do
	id=`echo $sample | perl -pe 's/\.\d+$//'`
	validation_bed=$valid_folder/estd20_NA${id}_5plus_probes.bed
	penn=$penn_folder/sample_*${sample}*.tabcnv

	penn_new=`echo $penn | sed 's/.tabcnv/_overlap.bed/'`

	dna=$dna_folder/sdundo_*${sample}*.txt
	dna_new=`echo $dna | sed 's/.txt/_overlap.bed/'`
	bedtools intersect -wo -a $validation_bed -b $dna > $dna_new
	bedtools intersect -wo -a $penn -b $dna > $results_folder/penn_sdundo_pvalue_${sample}_overlap.bed

	#Compare validated penn and validated DNAcopy
	cut -f 1-4 $penn_new | sort | uniq > foo.bed
	awk -F "\t" '$11 < -0.1 || $11 > 0.1' $dna_new | cut -f 1-4 | sort | uniq > bar.bed
	bedtools intersect -wo -a foo.bed -b bar.bed > $results_folder/validated_penn_vs_dna_sdundo_${sample}_overlap.bed
	rm foo.bed
	rm bar.bed

	dna=$dna_folder/segment_smoothed_*${sample}*.txt
	new=`echo $dna | sed 's/.txt/_overlap.bed/'`
	bedtools intersect -wo -a $validation_bed -b $dna > $new
	bedtools intersect -wo -a $penn -b $dna > $results_folder/penn_segment_smoothed_${sample}_overlap.bed

	cut -f 1-4 $penn_new | sort | uniq > foo.bed
	temp=$dna_folder/segment_*${sample}*_overlap.bed
	awk -F "\t" '$11 < -0.1 || $11 > 0.1' $temp | cut -f 1-4 | sort | uniq > bar.bed
	bedtools intersect -wo -a foo.bed -b bar.bed > $results_folder/validated_penn_vs_dna_segment_${sample}_overlap.bed
	rm foo.bed
	rm bar.bed

done


