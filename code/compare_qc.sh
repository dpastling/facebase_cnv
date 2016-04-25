#!/usr/bin/env bash
#BSUB -J compare
#BSUB -o logs/compare_qc_%J.out
#BSUB -e logs/compare_qc_%J.err

SAMPLES=(
"18508.1064691934"
"18508.1064692018"
"18858.1064691581"
"18858.1064691639"
"18861.1064691574"
"18861.1064691968"
"06985.1054736171"
"06985.1054736250"
"06985.1064733235"
"18502.1054736219"
"18502.1054736256"
"18505.1054736194"
"18505.1054736220"
"18505.1064692236"
)

valid_folder=data/validation
data_folder=results/logR
penn_folder=$data_folder/penncnv
dna_folder=$data_folder/dnacopy
vice_folder=$data_folder/vice
results_folder=results/vanilla_ice_qc

for sample in ${SAMPLES[@]}
do
    id=`echo $sample | perl -pe 's/\.\d+$//'`
    orig_id=`echo $sample | sed 's/\./\@/'`
    
    validation_bed=$valid_folder/estd20_NA${id}_5plus_probes.bed
    vice_file=${vice_folder}/*$orig_id.txt
    penn_file=${penn_folder}/*${orig_id}.txt
    dna_file=${dna_folder}/*${orig_id}.txt
    
    # for DNAcopy
    #results/dnacopy/segment_pvalue_X10.1064639586.txt | cut -f 2- | tail -n +2
    
    #n=`tail -n +2 $vice_file | sed 's/chr//' | awk '$7 < 3 || $7 > 3' | bedtools intersect -b $validation_bed -a - | wc -l`
    #n=`tail -n +2 $vice_file | sed 's/chr//' | awk '$7 < 3 || $7 > 3' | wc -l`
    #echo $sample $n
    
    vice_bed=$results_folder/${sample}_vice.bed
    dna_bed=$results_folder/${sample}_dna.bed
    penn_bed=$results_folder/${sample}_penn.bed
    tail -n +2 $vice_file | sed 's/chr//' | awk '$7 < 3 || $7 > 3' | grep -v "^[XYM]" > $vice_bed
    tail -n +2 $dna_file | cut -f 2- | grep -v "^[XYM]" > $dna_bed
    tail -n +2 $penn_file > $penn_bed

    # filter copy number 2 regions from the DNAcopy data
    awk '$11 < -0.1 || $11 > 0.1' $dna_bed > temp.bed
    mv temp.bed $dna_bed

    prefix=$results_folder/${sample}
    bedtools intersect -u -a $validation_bed -b $vice_bed > ${prefix}_valid_vice_overlap.txt
    bedtools intersect -u -a $validation_bed -b $penn_bed > ${prefix}_valid_penn_overlap.txt
    bedtools intersect -u -a $validation_bed -b $dna_bed > ${prefix}_valid_dna_overlap.txt
    bedtools intersect -u -a $vice_bed -b $penn_bed > ${prefix}_vice_penn_overlap.txt
    bedtools intersect -u -a $vice_bed -b $dna_bed > ${prefix}_vice_dna_overlap.txt
    bedtools intersect -u -a $vice_bed -b $validation_bed > ${prefix}_vice_valid_overlap.txt
    bedtools intersect -u -a $penn_bed -b $dna_bed > ${prefix}_penn_dna_overlap.txt
    bedtools intersect -u -a $penn_bed -b $vice_bed > ${prefix}_penn_vice_overlap.txt
    bedtools intersect -u -a $penn_bed -b $validation_bed > ${prefix}_penn_valid_overlap.txt
    bedtools intersect -u -a $dna_bed -b $penn_bed > ${prefix}_dna_penn_overlap.txt
    bedtools intersect -u -a $dna_bed -b $vice_bed > ${prefix}_dna_vice_overlap.txt
    bedtools intersect -u -a $dna_bed -b $validation_bed > ${prefix}_dna_valid_overlap.txt

    cp $validation_bed ${prefix}_valid.bed
    
done

## Get inverse overlaps, regions not overlapping with the other file
for bed in $results_folder/*.bed
do
    prefix=`echo $bed | sed 's/.bed//'`
    for x in ${prefix}*_overlap.txt
    do
        output=`echo $x | sed 's/_overlap/_inverse/'`
        bedtools intersect -v -a $bed -b $x > $output
    done
done


# print results to STDOUT
for sample in ${SAMPLES[@]}
do
    prefix=$results_folder/${sample}

    for x in ${prefix}_valid*
    do
        awk -F "\t" 'OFS="\t" {print $1, $2, $3, $6, "NA", "NA", "NA", "NA"}' $x | sed "s|^|$x\t|"
        #echo $x `awk -F "\t" '{l += $3 - $2; n += $6} END { print NR, l/NR, n/NR, "NA"}' $x`
    done
    
    for x in ${prefix}_vice*
    do
        awk -F "\t" 'OFS="\t" {print $1, $2, $3, $6, $8, $10, $11, $12}' $x | sed "s|^|$x\t|"
        #echo $x `awk -F "\t" '{l += $3 - $2; n += $6; s += $8} END { print NR, l/NR, n/NR, s/NR}' $x`
    done
    
    for x in ${prefix}_penn*
    do
        awk -F "\t" 'OFS="\t" {print $1, $2, $3, $9, $8, $10, $11, $12}' $x | sed "s|^|$x\t|"
        #echo $x `awk -F "\t" '{l += $3 - $2; n += $9; s += $8} END { print NR, l/NR, n/NR, s/NR}' $x`
    done
    
    for x in ${prefix}_dna*
    do
        awk -F "\t" 'OFS="\t" {print $1, $2, $3, $4, $7, $10, $11, $12}' $x | sed "s|^|$x\t|"
        #echo $x `awk '{l += $3 - $2; n += $4; s += $7} END { print NR, l/NR, n/NR, s/NR}' $x`
    done
done

