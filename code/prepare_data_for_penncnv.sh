#!/usr/bin/env bash
#BSUB -J split[3001-3739]%5
#BSUB -e logs/split_illumina_report_%J.err
#BSUB -o logs/split_illumina_report_%J.out
#BSUB -R "span[hosts=1]"
#BSUB -P Tamim
#BSUB -q short

SAMPLES=(data/Spritz_release_genotype_files/Spritz_Omni2-5PlusExome_release_FinalReport_*.csv)
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

split_illumina_report.pl --comma --prefix data/penncnv/FinalReport_ --suffix .txt ${sample}

