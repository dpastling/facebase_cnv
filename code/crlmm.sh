#!/usr/bin/env bash
#BSUB -J R
#BSUB -e logs/crlmm_%J.log
#BSUB -o logs/crlmm_%J.out
#BSUB -q normal
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -R "select[mem>51] rusage[mem=51]"
#BSUB -P Tamim

# catch unset variables, non-zero exits in pipes and calls, enable x-trace.
set -o nounset -o pipefail -o errexit -x


Rscript code/crlmm.R

