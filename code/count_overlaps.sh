#!/usr/bin/env bash

folder=data/validation
#folder=results/dnacopy_test
high_cutoff=0.1
low_cutoff=-0.1

for x in $folder/sample*_overlap.bed
do
n=`cut -f 7,8 $x | sort | uniq | wc -l`
v=`cut -f 4 $x | sort | uniq | wc -l`
echo $x $n $v
done

for x in $folder/segment_smoothed_*_overlap.bed
do
n=`awk -F "\t" '$11 < -0.1 || $11 > 0.1' $x | cut -f 7,8 | sort | uniq | wc -l`
v=`awk -F "\t" '$11 < -0.1 || $11 > 0.1' $x | cut -f 4 | sort | uniq | wc -l`
echo $x $n $v
done

for x in $folder/penn_*_overlap.bed
do
n=`awk -F "\t" '$14 < -0.1 || $14 > 0.1' $x | cut -f 10,11 | sort | uniq | wc -l`
v=`awk -F "\t" '$14 < -0.1 || $14 > 0.1' $x | cut -f 1,2 | sort | uniq | wc -l`
echo $x $n $v
done

for x in $folder/sdundo_pvalue_*_overlap.bed
do
n=`awk -F "\t" '$11 < -0.1 || $11 > 0.1' $x | cut -f 7,8 | sort | uniq | wc -l`
v=`awk -F "\t" '$11 < -0.1 || $11 > 0.1' $x | cut -f 4 | sort | uniq | wc -l`
echo $x $n $v
done

for x in $folder/validated_penn*_overlap.bed
do
v=`cut -f 4 $x | sort | uniq | wc -l`
echo $x $v
done


