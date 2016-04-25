#!/usr/bin/env bash

results="results/vanilla_ice_qc"
for bed in $results/*.bed
do
    prefix=`echo $bed | sed 's/.bed//'`
    for x in ${prefix}*_overlap.txt
    do
        output=`echo $x | sed 's/_overlap/_inverse/'`
        bedtools intersect -v -a $bed -b $x > $output
    done
done


