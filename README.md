
# Copy Number Variation (CNV) Analysis from SNP Array Data

- This code was used to generate CNV calls for the [Facebase project](https://www.facebase.org)
- calls using three tools: VanillaICE, pennCNV, and DNACopy. 
- Rather than relying on one algorithm, we seek consensus among the three for more robust results.

## Getting Started

### Prerequisites

- assumes data have been processed by Illumina GenomeStudio and results stored in FinalReport

**PennCNV**
[](http://penncnv.openbioinformatics.org/en/latest/)

The following R/Bioconductor packages can be installed within an R session by the following:

- DNAcopy
- VanillaICE
- ArrayTV


```r
source("https://bioconductor.org/biocLite.R")
biocLite(c("DNAcopy", "VanillaICE", "ArrayTV"))
```

- bedtools

[](http://bedtools.readthedocs.io/)


[https://support.illumina.com/array/array_kits/humanomni2-5exome-8-beadchip-kit.html](Illumina Omni2.5Exome support files)

[ftp://webdata2:webdata2@ussd-ftp.illumina.com/MyIllumina/a1b2a21a-38fe-4f5f-9b03-466cf066d85c/Omni25exome_MAF.zip]()
[ftp://webdata:webdata@ussd-ftp.illumina.com/Downloads/ProductFiles/HumanOmni2-5Exome-8/v1-0/HumanOmni2-5Exome-8-v1-0-B.bpm]()
[ftp://webdata:webdata@ussd-ftp.illumina.com/Downloads/ProductFiles/HumanOmni2-5Exome-8/v1-0/HumanOmni2-5Exome-8-v1-0-B.csv]()

### Installing

```
git clone https://github.com/dpastling/facebase_cnv
```

### Tesing

TODO: find a small example file online
TODO: also need support files 

## Running the code

- modify the config.sh file 

## License

This project is licensed under the MIT License - see the LICENSE.md file for details

## Acknowledgments

