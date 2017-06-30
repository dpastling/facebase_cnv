
# Copy Number Variation (CNV) Analysis from SNP Genotype Array Data

TODO: summarize the project here

- assumes data have been processed by commercial software and results stored in FinalReport
- This code was used to generate CNV calls for the [Facebase project](https://www.facebase.org)
- calls using three tools: VanillaICE, pennCNV, and DNACopy. 
- Rather than relying on one algorithm, we seek consensus among the three for more robust results.

## Getting Started

### Prerequisites

- PennCNV
[](http://penncnv.openbioinformatics.org/en/latest/)

The following R/Bioconductor packages
- DNAcopy
- VanillaICE
- ArrayTV

These can be installed within an R session by the following:

```r
source("https://bioconductor.org/biocLite.R")
biocLite(c("DNAcopy", "VanillaICE", "ArrayTV"))
```

- bedtools

[](http://bedtools.readthedocs.io/)


[https://support.illumina.com/array/array_kits/humanomni2-5exome-8-beadchip-kit.html](Illumina Omni2.5Exome support files)

ftp://webdata2:webdata2@ussd-ftp.illumina.com/MyIllumina/a1b2a21a-38fe-4f5f-9b03-466cf066d85c/Omni25exome_MAF.zip
ftp://webdata:webdata@ussd-ftp.illumina.com/Downloads/ProductFiles/HumanOmni2-5Exome-8/v1-0/HumanOmni2-5Exome-8-v1-0-B.bpm
ftp://webdata:webdata@ussd-ftp.illumina.com/Downloads/ProductFiles/HumanOmni2-5Exome-8/v1-0/HumanOmni2-5Exome-8-v1-0-B.csv

### Installing

```
git clone ????
```

### Tesing

TODO: find a small example file online
TODO: also need support files 

## Running the code

- modify the config.sh file 

## License

This project is licensed under the MIT License - see the LICENSE.md file for details

## Acknowledgments


---------------------------------------

## Overview of the analysis code

## Notes on the overlaps

## TODO items from March 17th, 2016 meeting

+ append LogR to the output of PennCNV and VanillaICE
- compare quality score metrics of validated CNVs versus un-validated
- compare other metrics such as CNV length, number of probes, for validated and non-validated CNVs
- calculate probe density: (number of probes) / (CNV length). Maybe probes per KB? or probes per 50kb?
- how many total probes in the genome.
- summarise probe density for the original validation set. Perhaps the CNVs that don't overlap have low probe density
- make a few plots of LogR for validated and non-validated CNVs

### Total number of probes

    2559461

    1000 * (2559461 / 3e9) = 0.8531537

There are roughly 0.85 probes per kb. Or 8.5 probes per 10 kb

