#!/usr/bin/env bash
#BSUB -J build
#BSUB -e logs/gcmodel_%J.log
#BSUB -o logs/crlmm_%J.out
#BSUB -q test
#BSUB -R "span[hosts=1]"
#BSUB -R "select[mem>51] rusage[mem=51]"
#BSUB -P Tamim

# catch unset variables, non-zero exits in pipes and calls, enable x-trace.
set -o nounset -o pipefail -o errexit -x

code/make_gcmodel.pl data/Human_Omni25exome.bed > data/Human_Omni25exome_GCModel.txt

