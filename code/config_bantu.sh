#!/usr/bin/env bash

DATA=data/penncnv
PENNCNV_RESULTS=results/penncnv
DNACOPY_RESULTS=results/dnacopy
VICE_RESULTS=results/vanilla_ice
SAMPLES=($DATA/FinalReport_*)

hmm_file=data/IlluminaHumanCoreExome_v12-A.hmm
pfb_file=data/Human_Omni25exome.pfb
gc_model=data/Human_Omni25exome_GCModel.txt
marker_info=data/Marker_Info_Files/HumanOmni25Exome-8v1_A.csv

