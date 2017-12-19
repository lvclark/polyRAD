# polyRAD: Genotype Calling with Uncertainty from Sequencing Data in Polyploids and Diploids

This R package is currently under development.

## Purpose

Genotypes derived from genotyping-by-sequencing (GBS) and restriction site-associated DNA sequencing (RAD-seq) have inherent uncertainty associated with them due to sampling error, i.e. some alleles might not get sequenced at all, or might not be sequenced in exact proportion to their copy number in the genome.  This package imports read depth in a variety of formats output by various bioinformatics pipelines and estimates the probability of each possible genotype for each taxon and locus.  Unlike similar pipelines, polyRAD can account for population structure and variable inheritance modes (autopolyploid, allopolyploid, intermediate).  Gentoypes and/or probability distributions can then be exported for downstream analysis such as genome-wide association, genomic selection, QTL mapping, or population structure analysis.

## Funding

This material is based upon work supported by the National Science Foundation under Grant No. 
[1661490](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1661490&HistoricalAwards=false).  
Any opinions, findings, and conclusions or recommendations expressed in this material are those 
of the author(s) and do not necessarily reflect the views of the National Science Foundation.
