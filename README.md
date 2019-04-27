# polyRAD: Genotype Calling with Uncertainty from Sequencing Data in Polyploids and Diploids

This R package is ready for use, although new features are still being developed.  See the [list of future features](https://github.com/lvclark/polyRAD/wiki/todo) to see where it is headed.

polyRAD is part of an upcoming suite of interoperable software called the
[ploidyverse](https://github.com/ploidyverse).

I'm always interested in new collaboration!  If you have new feature requests I would love to hear them, especially if you have datasets that you would be willing to share that could be used for testing these new features.  I expect that there will be more publications as polyRAD is developed further, and I am open to adding new coauthors.  Contact: Lindsay Clark, University of Illinois, Urbana-Champaign.

## Purpose

Genotypes derived from genotyping-by-sequencing (GBS) and restriction site-associated DNA sequencing (RAD-seq) have inherent uncertainty associated with them due to sampling error, i.e. some alleles might not get sequenced at all, or might not be sequenced in exact proportion to their copy number in the genome.  This package imports read depth in a variety of formats output by various bioinformatics pipelines and estimates the probability of each possible genotype for each taxon and locus.  Unlike similar pipelines, polyRAD can account for population structure and variable inheritance modes (autopolyploid, allopolyploid, intermediate).  Genotypes and/or probability distributions can then be exported for downstream analysis such as genome-wide association, genomic selection, QTL mapping, or population structure analysis.

## Why polyRAD?

If you're like me, you don't want to waste a lot of money sequencing your DNA samples at a higher depth than is necessary.  You would rather spend that money adding more samples to the project, or using a different restriction enzyme to get more markers!  You may have also noticed that some loci get sequenced at a much higher depth than others, which means that if you sequence the same library a second time, you aren't likely to get a lot of reads for the loci that need it most.  So how can we get the maximum amount of information out of sequencing data where many loci are low depth?  And, for example, if we only have five reads, how can we estimate allele dosage in a heterozygous octoploid?

The answer that polyRAD provides is a Bayesian genotype caller with many options for specifying genotype prior probabilties.  When read depth is low, **accurate priors make a big difference in the accuracy of genotype calls.**  And because some genotype calls are going to be uncertain no matter how sophisticated our algorithm is, polyRAD can export genotypes as continuous numeric variables reflecting the probabilities of all possible allele copy numbers.  This includes genotypes with zero reads, where the priors themselves are used for imputation.

**Genotype priors in diversity panels and natural populations:**

* Either assume no population structure (HWE), or let polyRAD infer population structure and model allele frequency gradients.
* The user specifies a rate of self-fertilization ranging anywhere from zero to one.
* If loci have known positions in a reference genome, polyRAD can search for loci in linkage disequilibrium and use those loci to update priors.

**Genotype priors in biparental mapping populations:**

* The user specifies the number of generations of backcrossing, intermating, and/or selfing.  (For an F1 population, all three would be zero.)
* Based on likely parental genotypes and allele frequencies in the progeny, polyRAD determines the segregation pattern of each marker.
* If the loci have known positions in a reference genome or on a map, linked markers can be used for updating priors.

In particular, by using population structure and linkage to inform genotype priors on a per-individual basis, high depth markers are used by polyRAD to improve the accuracy of genotyping at low depth markers.  All pipelines allow autopolyploidy, allopolyploidy, or some mixture of the two.  And because non-model organisms need some love, reference genomes are optional.

## Formats supported

To hopefully answer the question, "Can I use polyRAD?":

polyRAD requires as input the sequence read depth at each allele for each sample.  Alleles must also be grouped into loci.  The bioinformatics pipeline that you used for SNP discovery did not have to assume polyploidy, as long as it faithfully reported allelic read depth.  Genomic alignment information is optional.  Right now there are data import functions for the following formats:

* [Variant Call Format (VCF)](https://samtools.github.io/hts-specs/).  The allele depth (AD) genotype field must be present.  I have tested the import function on files produced by the [TASSEL GBSv2](https://bitbucket.org/tasseladmin/tassel-5-source/wiki/Tassel5GBSv2Pipeline) pipeline.  It should also work for [GATK](https://software.broadinstitute.org/gatk/) and [Tassel4-Poly](https://github.com/guilherme-pereira/tassel4-poly).
* [TagDigger](https://github.com/lvclark/tagdigger).  This is another piece of software I created, which reads FASTQ files and expected tag sequences and outputs a CSV file with read counts at samples x tags.
* [UNEAK](https://tassel.bitbucket.io/TasselArchived.html).  The UNEAK pipeline outputs read depth in a file called HapMap.hmc.txt, which can be read by polyRAD.  (Beware that read depth is capped at 127 by UNEAK; TagDigger can help you if you expect high depth to be common in your dataset.)
* [Stacks](http://catchenlab.life.illinois.edu/stacks/).  If you have catalog files (catalog.alleles.tsv etc.) and matches files (matches.tsv) generated by Stacks, they can be imported by polyRAD.
* [TASSEL-GBSv2](https://bitbucket.org/tasseladmin/tassel-5-source/wiki/Tassel5GBSv2Pipeline).  Rather than running the entire pipeline, you can run GBSSeqToTagDBPlugin, TagExportToFastqPlugin, Bowtie2 or BWA as described in "Run Alignment Program(s)", then GetTagTaxaDistFromDBPlugin to get the TagTaxaDist and SAM files needed for import by polyRAD.

Currently there are export functions for the following software.  Genotypes are exported as continuous variables for these three formats.  There are also functions to generate matrices of continuous or discrete genotypes, which can be used in custom export functions.

* [GAPIT](http://www.zzlab.net/GAPIT/) and [FarmCPU](http://www.zzlab.net/FarmCPU/)
* [rrBLUP](https://cran.r-project.org/web/packages/rrBLUP/)
* [TASSEL](http://www.maizegenetics.net/tassel)

There are export functions for discrete genotypes for the following software:

* [polymapR](https://cran.r-project.org/package=polymapR)
* [GWASpoly](https://potatobreeding.cals.wisc.edu/software/)

Genotype probabilities can be exported to:

* [MAPpoly](https://github.com/mmollina/MAPpoly)

## Installation

polyRAD depends on some Bioconductor packages.  Before attempting to install polyRAD, run

```
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

BiocManager::install("pcaMethods")
```

If you plan to import from VCF, also run

```
BiocManager::install("VariantAnnotation")
```

polyRAD can then be installed from CRAN with

```
install.packages("polyRAD")
```

Alternatively, if there are new features not yet on the CRAN version that you
want to use, you can install the development version here on GitHub at your own
risk.  There are R packages such as `devtools` and `githubinstall` that 
facilitate installing directly from GitHub.

## Tutorial

The tutorial document for the package is available [on Github](https://github.com/lvclark/polyRAD/blob/master/vignettes/polyRADtutorial.md).

## Citation

If you use polyRAD, please cite this manuscript:

Clark LV, Lipka AE, and Sacks EJ (2019) polyRAD: Genotype calling with uncertainty from sequencing data in
polyploids and diploids.  *G3* 9(3) 663--673. doi:[10.1534/g3.118.200913](https://doi.org/10.1534/g3.118.200913)

### Additional citations

Citable Zenodo DOI for the software:
[![DOI](https://zenodo.org/badge/99379777.svg)](https://zenodo.org/badge/latestdoi/99379777)

polyRAD has also been presented in these posters:

Clark LV, Lipka AE, and Sacks EJ (2019) Improvements to Genotype Calling in 
Polyploids Using the polyRAD R Package.  Plant and Animal Genome Conference 
XXVII, January 12-16, San Diego, California, USA.
doi:[10.13140/RG.2.2.18358.75847](https://doi.org/10.13140/RG.2.2.18358.75847)

Clark LV, Lipka AE, and Sacks EJ (2018) polyRAD: Genotype Calling with Uncertainty from Sequencing Data
in Polyploids and Diploids.  Plant and Animal Genome Conference XXVI, January 13-17, San Diego, California, USA.
doi:[10.13140/RG.2.2.27134.08001](https://doi.org/10.13140/RG.2.2.27134.08001)

## Funding

This material is based upon work supported by the National Science Foundation under Grant No. 
[1661490](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1661490&HistoricalAwards=false).  
Any opinions, findings, and conclusions or recommendations expressed in this material are those 
of the author(s) and do not necessarily reflect the views of the National Science Foundation.
