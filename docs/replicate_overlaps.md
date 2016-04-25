[astlingd@compute12 vanilla_ice_qc]$ ll 18505.*_penn.bed
	-rw-r----- 1 astlingd analysiscore 102245 Apr 11 13:50 18505.1054736194_penn.bed
	-rw-r----- 1 astlingd analysiscore  46734 Apr 11 13:50 18505.1054736220_penn.bed
	-rw-r----- 1 astlingd analysiscore 167551 Apr 11 13:50 18505.1064692236_penn.bed
[astlingd@compute12 vanilla_ice_qc]$ file1=18505.1054736194_penn.bed
[astlingd@compute12 vanilla_ice_qc]$ file2=18505.1054736220_penn.bed
[astlingd@compute12 vanilla_ice_qc]$ file3=18505.1064692236_penn.bed
[astlingd@compute12 vanilla_ice_qc]$ 
[astlingd@compute12 vanilla_ice_qc]$ bedtools intersect -a $file1 -b $file2 > 1_2_overlap.bed
[astlingd@compute12 vanilla_ice_qc]$ bedtools intersect -a $file2 -b $file3 > 2_3_overlap.bed
[astlingd@compute12 vanilla_ice_qc]$ bedtools intersect -a $file1 -b $file3 > 1_3_overlap.bed
[astlingd@compute12 vanilla_ice_qc]$ bedtools intersect -a 1_2_overlap.bed -b $file3 > all_overlap.bed
[astlingd@compute12 vanilla_ice_qc]$ wc -l 18505.*_penn.bed
	   724 18505.1054736194_penn.bed
	   332 18505.1054736220_penn.bed
	  1188 18505.1064692236_penn.bed
	  2244 total
[astlingd@compute12 vanilla_ice_qc]$ wc -l 1_2_overlap.bed 
	269 1_2_overlap.bed
[astlingd@compute12 vanilla_ice_qc]$ wc -l 2_3_overlap.bed 
	268 2_3_overlap.bed
[astlingd@compute12 vanilla_ice_qc]$ wc -l 1_3_overlap.bed 
	620 1_3_overlap.bed
[astlingd@compute12 vanilla_ice_qc]$ wc -l all_overlap.bed 
	249 all_overlap.bed
[astlingd@compute12 vanilla_ice_qc]$ ls 18505*_valid.bed
	18505.1054736194_valid.bed  18505.1054736220_valid.bed  18505.1064692236_valid.bed
[astlingd@compute12 vanilla_ice_qc]$ valid=18505.1054736194_valid.bed
[astlingd@compute12 vanilla_ice_qc]$ 
[astlingd@compute12 vanilla_ice_qc]$ bedtools intersect -a $file1 -b $valid | wc -l
	65
[astlingd@compute12 vanilla_ice_qc]$ bedtools intersect -a $file2 -b $valid | wc -l
	65
[astlingd@compute12 vanilla_ice_qc]$ bedtools intersect -a $file3 -b $valid | wc -l
	63
[astlingd@compute12 vanilla_ice_qc]$ bedtools intersect -a 1_2_overlap.bed -b $valid | wc -l
	58
[astlingd@compute12 vanilla_ice_qc]$ bedtools intersect -a 2_3_overlap.bed -b $valid | wc -l
	56
[astlingd@compute12 vanilla_ice_qc]$ bedtools intersect -a 1_3_overlap.bed -b $valid | wc -l
	56
[astlingd@compute12 vanilla_ice_qc]$ bedtools intersect -a all_overlap.bed -b $valid | wc -l
	53
[astlingd@compute12 vanilla_ice_qc]$ wc -l $valid
	230 18505.1054736194_valid.bed

