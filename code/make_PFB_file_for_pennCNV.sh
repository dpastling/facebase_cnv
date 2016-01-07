#!/usr/bin/env bash

echo "SNP,Chr,Position,PFB" | perl -pe 's/,/\t/g' > data/header.pfb
cut -f 1-4 data/Marker_Info_Files/Allele_Frequencies/Omni25exome_MAF.txt | grep -v "^Name" > data/temp.pfb
cat data/header.pfb data/temp.pfb > data/Human_Omni25exome.pfb
rm data/header.pfb
rm data/temp.pfb

