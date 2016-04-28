
# Facebase Analysis

TODO: summarize the project here

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

