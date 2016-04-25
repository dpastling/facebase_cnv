#!/usr/bin/env bash

18508.1064691934
18508.1064692018

18858.1064691581
18858.1064691639

18861.1064691574
18861.1064691968

18502.1054736219
18502.1054736256

18505.1054736194
18505.1054736220
18505.1064692236

sample=06985
replicate1=1054736171
replicate2=1054736250
replicate3=1064733235

sample=18505
replicate1=1054736194
replicate2=1054736220
replicate3=1064692236

method=vice
sample1=$sample.${replicate1}_${method}.bed
sample2=$sample.${replicate2}_${method}.bed
sample3=$sample.${replicate3}_${method}.bed

bedtools intersect -wo -a $sample1 -b $sample2 | \
  awk '($11 > 0 && $23 > 0) || ($11 <0 && $23 < 0)' | \
  cut -f 1-12 > ${sample}_1_2.bed

bedtools intersect -wo -a $sample2 -b $sample3 | \
  awk '($11 > 0 && $23 > 0) || ($11 <0 && $23 < 0)' | \
  cut -f 1-12 > ${sample}_2_3.bed

bedtools intersect -wo -a $sample1 -b $sample3 | \
  awk '($11 > 0 && $23 > 0) || ($11 <0 && $23 < 0)' | \
  cut -f 1-12 > ${sample}_1_3.bed

bedtools intersect -wo -a ${sample}_1_2.bed -b $sample3 | \
  awk '($11 > 0 && $23 > 0) || ($11 <0 && $23 < 0)' | \
  cut -f 1-12 > ${sample}_all_${method}.bed


