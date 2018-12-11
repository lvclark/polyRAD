# polyRAD tutorial

## Table of Contents
* [Introduction](#introduction)
* [Summary of available functions](#functions)
* [Estimating genotype probabilities in a mapping population](#mapping)
* [Estimating genotype probabilities in a diversity panel](#diversity)
* [Considerations for RAM and processing time](#considerations)
* [Citing polyRAD](#citation)

Introduction <a name="introduction"></a>
----------------------------------------

polyRAD is an R package that assists with genotype calling from DNA
sequence datasets such as genotyping-by-sequencing (GBS) or restriction
site-associated DNA sequencing (RAD) in polyploids and diploids.
Genotype likelihoods are estimated from allelic read depth, genotype
prior probabilities are estimated from population parameters, and then
genotype posterior probabilities are estimated from likelihoods and
prior probabilities. Posterior probabilities can be used directly in
downstream analysis, converted to posterior mean genotypes for analyses
of additive genetic effects, or used for export of the most probable
genotypes for analyses that require discrete genotypic data.

Analyses in polyRAD center around objects of an S3 class called
`RADdata`. A single `RADdata` object contains the entire dataset of read
depth and locus information, as well as parameters that are estimated
during the course of analysis.

Summary of available functions <a name="functions"></a>
-------------------------------------------------------

For any function named in this section, see its help page for more
information. (For example by typing `?VCF2RADdata` into the R console.)

Several functions are available for import of read depth data and
(optionally) alignment information into a RADdata object:

-   `VCF2RADdata`
-   `readTagDigger`
-   `readStacks`
-   `readHMC`
-   `readTASSELGBSv2`

More generally, the `RADdata` function is used for constructing RADdata
objects; see the help page for that function for more information on
what data are needed.

Several pipelines are available for genotype estimation, depending on
how the population is structured (i.e. what the genotype prior
probabilities should be.):

-   `PipelineMapping2Parents`
-   `IterateHWE`
-   `IterateHWE_LD`
-   `IteratePopStruct`
-   `IteratePopStructLD`

For exporting the estimated genotypes to other software:

-   `ExportGAPIT`
-   `Export_rrBLUP_Amat`
-   `Export_rrBLUP_GWAS`
-   `Export_TASSEL_Numeric`
-   `Export_polymapR`

If you need continuous numerical genotypes exported in some other
format, see `GetWeightedMeanGenotypes`. If you need discrete numerical
genotypes, see `GetProbableGenotypes`. Also, `GetLikelyGen` returns the
most likely genotypes (based on read depth only) for a single sample.

There are also various utilities for manipulating RADdata objects:

-   `SubsetByTaxon`
-   `SubsetByLocus`
-   `SplitByChromosome`
-   `MergeRareHaplotypes`
-   `RemoveMonomorphicLoci`
-   `RemoveHighDepthLoci`
-   `EstimateContaminationRate`
-   `StripDown`
-   `LocusInfo`

See `?GetTaxa` for a list of accessor functions as well.

Estimating genotype probabilities in a mapping population <a name="mapping"></a>
--------------------------------------------------------------------------------

In this example, we’ll import some data from an F1 mapping population of
*Miscanthus sinensis* that were output by the
[UNEAK](https://doi.org/10.1371/journal.pgen.1003215) pipeline. These
data are from a study by Liu *et al.* (2015;
[doi:10.1111/gcbb.12275](https://doi.org/10.1111/gcbb.12275); data
available at <http://hdl.handle.net/2142/79522>), and can be found in
the “extdata” folder of the polyRAD installation. *Miscanthus* is an
ancient tetraploid that has undergone diploidization. Given the ability
of the UNEAK pipeline to filter paralogs, we expect most loci to behave
in a diploid fashion, but some may behave in an allotetraploid fashion.

We’ll start by loading polyRAD and importing the data into a `RADdata`
object. The `possiblePloidies` argument indicates the expected
inheritance modes: diploid (2) and allotetraploid (2 2).

With your own dataset, you will not need to use `system.file`. Instead,
directly create a text string indicating the name of your file (and its
location if it is not in the working directory.)

``` r
library(polyRAD)
maphmcfile <- system.file("extdata", "ClareMap_HapMap.hmc.txt", 
                          package = "polyRAD")
maphmcfile
```

    ## [1] "C:/Users/lvclark/Documents/R/win-library/3.5/polyRAD/extdata/ClareMap_HapMap.hmc.txt"

``` r
mydata <- readHMC(maphmcfile,
                  possiblePloidies = list(2, c(2, 2)))
mydata
```

    ## ## RADdata object ##
    ## 299 taxa and 50 loci
    ## 766014 total reads
    ## Assumed sample cross-contamination rate of 0.001
    ## 
    ## Possible ploidies:
    ## Autodiploid (2)
    ## Allotetraploid (2 2)

We can view the imported taxa names (subsetted here for space).

``` r
GetTaxa(mydata)[c(1:10,293:299)]
```

    ##  [1] "IGR-2011-001"    "Kaskade-Justin"  "Map1-001"       
    ##  [4] "Map1-002"        "Map1-003"        "Map1-005"       
    ##  [7] "Map1-008"        "Map1-011"        "Map1-016"       
    ## [10] "Map1-018"        "Map1-488"        "Map1-489"       
    ## [13] "Map1-490"        "Map1-491"        "Zebrinus-Justin"
    ## [16] "p196-150A-c"     "p877-348-b"

All names starting with “Map” are progeny. “Kaskade-Justin” and
“Zebrinus-Justin” are the parents. “IGR-2011-001”, “p196-150A-c”“,
and”p877-348-b" aren’t part of the population, but were doubled haploid
lines that were used to screen for paralogous markers. We can tell
polyRAD which taxa are the parents; since this is an F1 population it
doesn’t matter which is “donor” and which is “recurrent”.

``` r
mydata <- SetDonorParent(mydata, "Kaskade-Justin")
mydata <- SetRecurrentParent(mydata, "Zebrinus-Justin")
```

The next thing we’ll want to do is add our genomic alignment data. For
this dataset, we have alignment data stored in a CSV file, also in the
“extdata” directory of the polyRAD installation. We’ll add it to the
`locTable` slot of our `RADdata` object. Be sure to name the new columns
“Chr” and “Pos”.

``` r
alignfile <- system.file("extdata", "ClareMap_alignments.csv", 
                         package = "polyRAD")

aligndata <- read.csv(alignfile, row.names = 1)
head(aligndata)
```

    ##         Sorghum.LG Position.on.Sorghum.LG..bp.
    ## TP5489           1                     4560204
    ## TP13305          1                     4584260
    ## TP18261          1                     2911329
    ## TP18674          1                      387849
    ## TP19030          1                     7576879
    ## TP26698          1                     6972841

``` r
mydata$locTable$Chr <- aligndata[GetLoci(mydata), 1]
mydata$locTable$Pos <- aligndata[GetLoci(mydata), 2]
head(mydata$locTable)
```

    ##         Chr     Pos
    ## TP5489    1 4560204
    ## TP13305   1 4584260
    ## TP18261   1 2911329
    ## TP18674   1  387849
    ## TP19030   1 7576879
    ## TP26698   1 6972841

If you don’t have alignment data in your own dataset, you can still use
the pipeline described here. Just set `useLinkage = FALSE` in the code
below. The advantage of including alignment data is that gentoypes at
linked markers are used for imputing missing or correcting erroneous
genotypes.

It is important that the only individuals included in the analysis are
those that are truly part of the population. Allele frequencies are used
for inferring segregation pattern, and could be incorrect if many
individuals are included that are not part of the population.
Additionally, the genotype priors will be incorrect for individuals that
are not part of the population, leading to incorrect genotypes.

At this point we would normally do

``` r
mydata <- AddPCA(mydata)
```

However, because a very small number of markers was used in this example
dataset, the PCA does not accurately reflect the relatedness of
individuals. Here I will load a PCA that was done with the full set of
markers.

``` r
load(system.file("extdata", "examplePCA.RData", package = "polyRAD"))
mydata$PCA <- examplePCA
```

Now a plot can be used for visualizing the relationship among taxa.

``` r
plot(mydata)
```

![](polyRADtutorial_files/figure-markdown_github/unnamed-chunk-7-1.png)

Now we’ll extract a subset of taxa that we actually want to analyze. We
can see from the plot that a fair number of them were the product of
self-fertilization of “Zebrinus-Justin” and should be eliminated.

``` r
realprogeny <- GetTaxa(mydata)[mydata$PCA[,"PC1"] > -10 &
                                 mydata$PCA[,"PC1"] < 10]
# eliminate the one doubled haploid line in this group
realprogeny <- realprogeny[!realprogeny %in% c("IGR-2011-001", "p196-150A-c",
                                               "p877-348-b")]
# also retain parents
keeptaxa <- c(realprogeny, GetDonorParent(mydata), GetRecurrentParent(mydata))

mydata <- SubsetByTaxon(mydata, taxa = keeptaxa)
plot(mydata)
```

![](polyRADtutorial_files/figure-markdown_github/unnamed-chunk-8-1.png)

Now we can perform a preliminary run of the pipeline. The
`allowedDeviation` argument indicates how different the apparent allele
frequency (based on read depth ratios) can be from an expected allele
frequency (determined based on ploidy and mapping population type) and
still be classified as that allele frequency. The default settings
assume an F1 population, but the population type can be adjusted using
the `n.gen.backcrossing`, `n.gen.intermating`, and `n.gen.selfing`
arguments. We’ll also lower `minLikelihoodRatio` from the default
because one of the parents has many uncertain genotypes under the
tetraploid model (which was determined by exploration of the dataset
outside of this tutorial; many NA values were observed in `priorProb`
under the default). Since this first run is for a rough estimate of
genotypes, we’ll set `useLinkage = FALSE` to save a little computational
time.

``` r
mydata2 <- PipelineMapping2Parents(mydata, 
                                   freqAllowedDeviation = 0.06,
                                   useLinkage = FALSE,
                                   minLikelihoodRatio = 2)
```

    ## Making initial parameter estimates...

    ## Done.

We can use these preliminary estimates to determine whether we need to
adjust the overdispersion parameter. Exactly how much does read depth
distribution deviate from what would be expected under binomial
distibution? The `TestOverdispersion` function will help us here. We
will use the `qqman` package to visualize the results.

``` r
library(qqman)
overdispersionP <- TestOverdispersion(mydata2, to_test = 6:10)
qq(overdispersionP[["6"]])
```

![](polyRADtutorial_files/figure-markdown_github/unnamed-chunk-10-1.png)

``` r
qq(overdispersionP[["7"]])
```

![](polyRADtutorial_files/figure-markdown_github/unnamed-chunk-10-2.png)

``` r
qq(overdispersionP[["8"]])
```

![](polyRADtutorial_files/figure-markdown_github/unnamed-chunk-10-3.png)

``` r
qq(overdispersionP[["9"]])
```

![](polyRADtutorial_files/figure-markdown_github/unnamed-chunk-10-4.png)

``` r
qq(overdispersionP[["10"]])
```

![](polyRADtutorial_files/figure-markdown_github/unnamed-chunk-10-5.png)

It looks like `9` follows the red line most closely, so we’ll use that
for the overdispersion parameter. Now we can re-run the pipeline to
properly call the genotypes.

``` r
mydata <- PipelineMapping2Parents(mydata, 
                                  freqAllowedDeviation = 0.06,
                                  useLinkage = TRUE, overdispersion = 9,
                                  minLikelihoodRatio = 2)
```

    ## Making initial parameter estimates...

    ## Updating priors using linkage...

    ## Done.

We can examine the allele frequencies. Allele frequencies that fall
outside of the expected ranges will be recorded as they were estimated
from read depth. In this case all are within the expected ranges.

``` r
table(mydata$alleleFreq)
```

    ## 
    ## 0.125  0.25 0.375   0.5 0.625  0.75 0.875 
    ##     1    42     1    12     1    42     1

Genotype likelihood is also stored in the object for each possible
genotype at each locus, taxon, and ploidy. This is the probability of
seeing the observed distribution of reads.

``` r
mydata$alleleDepth[8,19:26]
```

    ## TP31810_0 TP31810_1 TP34632_0 TP34632_1 TP34939_0 TP34939_1 TP35570_0 
    ##        26        33        12        18         9         5        18 
    ## TP35570_1 
    ##         4

``` r
mydata$genotypeLikelihood[[1]][,8,19:26]
```

    ##      TP31810_0    TP31810_1    TP34632_0    TP34632_1    TP34939_0
    ## 0 1.274105e-06 5.796763e-07 1.819889e-05 3.238446e-07 3.064777e-06
    ## 1 3.520993e-02 3.520993e-02 6.090279e-02 6.090279e-02 1.084014e-01
    ## 2 5.796763e-07 1.274105e-06 3.238446e-07 1.819889e-05 3.431537e-05
    ##      TP34939_1    TP35570_0    TP35570_1
    ## 0 3.431537e-05 2.688637e-08 2.257439e-04
    ## 1 1.084014e-01 2.751453e-02 2.751453e-02
    ## 2 3.064777e-06 2.257439e-04 2.688637e-08

``` r
mydata$genotypeLikelihood[[2]][,8,19:26]
```

    ##      TP31810_0    TP31810_1    TP34632_0    TP34632_1    TP34939_0
    ## 0 1.274105e-06 5.796763e-07 1.819889e-05 3.238446e-07 3.064777e-06
    ## 1 1.733971e-02 6.779369e-03 4.353570e-02 1.028542e-02 2.003410e-02
    ## 2 3.520993e-02 3.520993e-02 6.090279e-02 6.090279e-02 1.084014e-01
    ## 3 6.779369e-03 1.733971e-02 1.028542e-02 4.353570e-02 1.067533e-01
    ## 4 5.796763e-07 1.274105e-06 3.238446e-07 1.819889e-05 3.431537e-05
    ##      TP34939_1    TP35570_0    TP35570_1
    ## 0 3.431537e-05 2.688637e-08 2.257439e-04
    ## 1 1.067533e-01 1.128698e-03 1.118037e-01
    ## 2 1.084014e-01 2.751453e-02 2.751453e-02
    ## 3 2.003410e-02 1.118037e-01 1.128698e-03
    ## 4 3.064777e-06 2.257439e-04 2.688637e-08

Above, for one individal (Map1-018), we see its read depth at eight
alleles (four loci), followed by the genotype likelihoods under diploid
and tetraploid models. For example, at locus TP35570, heterozygosity is
the most likely state, although there is a chance that this individual
is homozygous for allele 0 and the four reads of allele 1 were due to
contamination. If this locus is allotetraploid, it is most likely that
there is one copy of allele 1 and three copies of allele 0. Other loci
have higher depth and as a result there is less uncertainty in the
genotype, particularly for the diploid model.

The prior genotype probabilities (expected genotype distributions) are
also stored in the object for each possible ploidy. These distributions
are estimated based on the most likely parent genotypes. Low confidence
parent genotypes can be ignored by increasing the `minLikelihoodRatio`
argument to `PipelineMapping2Parents`.

``` r
mydata$priorProb[[1]][,19:26]
```

    ##   TP31810_0 TP31810_1 TP34632_0 TP34632_1 TP34939_0 TP34939_1 TP35570_0
    ## 0       0.5       0.0       0.0       0.5       0.0       0.5      0.25
    ## 1       0.5       0.5       0.5       0.5       0.5       0.5      0.50
    ## 2       0.0       0.5       0.5       0.0       0.5       0.0      0.25
    ##   TP35570_1
    ## 0      0.25
    ## 1      0.50
    ## 2      0.25

``` r
mydata$priorProb[[2]][,19:26]
```

    ##   TP31810_0 TP31810_1 TP34632_0 TP34632_1 TP34939_0 TP34939_1 TP35570_0
    ## 0         0         0         0         0        NA        NA       0.0
    ## 1         1         0         0         1        NA        NA       0.0
    ## 2         0         0         0         0        NA        NA       0.5
    ## 3         0         1         1         0        NA        NA       0.5
    ## 4         0         0         0         0        NA        NA       0.0
    ##   TP35570_1
    ## 0       0.0
    ## 1       0.5
    ## 2       0.5
    ## 3       0.0
    ## 4       0.0

Here we see some pretty big differences under the diploid and
allotetraploid models. For example, if TP35570 is behaving in a diploid
fashion we expect F2-like segregation since both parents were
heterozygous. However, if TP35570 is behaving in an allotetraploid
fashion, a 1:1 segregation ratio is expected due to one parent being
heterozygous at one isolocus and the other being homozygous at both
isoloci.

Now we want to determine which ploidy is the best fit for each locus.
This is done by comparing genotype prior probabilities to genotype
likelihoods and estimating a *χ*<sup>2</sup> statistic. Lower values
indicate a better fit.

``` r
mydata$ploidyChiSq[,19:26]
```

    ##        TP31810_0   TP31810_1  TP34632_0  TP34632_1 TP34939_0 TP34939_1
    ## [1,]   0.1262056   0.1262056   2.654915   2.654915 0.8418978 0.8418978
    ## [2,] 197.3757006 197.3757006 189.289082 189.289082        NA        NA
    ##      TP35570_0 TP35570_1
    ## [1,]  5.818905  5.818905
    ## [2,] 88.956011 88.956011

We can make a plot to get an overall sense of how well the markers fit
the diploid versus tetraploid model.

``` r
plot(mydata$ploidyChiSq[1,], mydata$ploidyChiSq[2,], 
     xlab = "Chi-squared for diploid model",
     ylab = "Chi-squared for tetraploid model")
```

![](polyRADtutorial_files/figure-markdown_github/unnamed-chunk-16-1.png)

For each allele, whichever model gives the lower Chi-squared value is
the one with the best fit. In this case it looks like everything is
diploid with fairly high confidence.

Now we’ll examine the posterior genotype probabilities. These are still
estimated separately for each ploidy.

``` r
mydata$posteriorProb[[1]][,10,19:26]
```

    ##      TP31810_0    TP31810_1    TP34632_0    TP34632_1    TP34939_0
    ## 0 2.166937e-06 0.000000e+00 0.000000e+00 1.000000e+00 0.000000e+00
    ## 1 9.999978e-01 9.999978e-01 3.567191e-19 3.567191e-19 8.693566e-14
    ## 2 0.000000e+00 2.166937e-06 1.000000e+00 0.000000e+00 1.000000e+00
    ##      TP34939_1    TP35570_0    TP35570_1
    ## 0 1.000000e+00 1.930766e-22 9.999999e-01
    ## 1 8.693566e-14 1.418476e-07 1.418476e-07
    ## 2 0.000000e+00 9.999999e-01 1.930766e-22

``` r
mydata$posteriorProb[[2]][,10,19:26]
```

    ##   TP31810_0 TP31810_1 TP34632_0 TP34632_1 TP34939_0 TP34939_1    TP35570_0
    ## 0         0         0         0         0       NaN       NaN 0.0000000000
    ## 1         1         0         0         1       NaN       NaN 0.0000000000
    ## 2         0         0         0         0       NaN       NaN 0.0006169573
    ## 3         0         1         1         0       NaN       NaN 0.9993830427
    ## 4         0         0         0         0       NaN       NaN 0.0000000000
    ##      TP35570_1
    ## 0 0.0000000000
    ## 1 0.9993830427
    ## 2 0.0006169573
    ## 3 0.0000000000
    ## 4 0.0000000000

We can export the results for use in downstream analysis. The function
below weights possible ploidies for each allele based on the results in
`mydata$ploidyChiSq`, and for each taxon outputs a continuous, numerical
genotype that is the mean of all possible genotypes weighted by genotype
posterior probabilities (*i.e.* the posterior mean genotype). By
default, one allele per locus is discarded in order to avoid
mathematical singularities in downstream analysis. The continuous
genotypes also range from zero to one by default, which can be changed
with the `minval` and `maxval` arguments.

``` r
mywm <- GetWeightedMeanGenotypes(mydata)
mywm[c(276, 277, 1:5), 10:13]
```

    ##                    TP31810_0    TP34632_1    TP34939_1  TP35570_1
    ## Kaskade-Justin  7.155005e-07 5.005078e-01 5.001090e-01 0.49322780
    ## Zebrinus-Justin 5.000697e-01 5.513391e-06 4.831070e-07 0.49893123
    ## Map1-001        4.998047e-01 3.457929e-03 8.646485e-18 0.49290029
    ## Map1-002        1.625360e-04 3.457929e-03 1.606743e-16 0.01542878
    ## Map1-003        6.890262e-03 4.965421e-01 5.000000e-01 0.02755674
    ## Map1-005        4.998177e-01 3.457929e-03 1.314400e-15 0.49444951
    ## Map1-008        1.616524e-04 4.965420e-01 1.514873e-08 0.01535917

Note that the parent posterior mean genotypes were estimated using
gentoype likelihood only, ignoring the priors set for the progeny. In
some places they may not match the progeny genotypes, indicating a
likely error in parental genotype calling. We can see the parental
genotypes that were used for estimating progeny priors using
`$likelyGeno_donor` and `$likelyGeno_recurrent`.

``` r
mydata$likelyGeno_donor[,19:26]
```

    ##   TP31810_0 TP31810_1 TP34632_0 TP34632_1 TP34939_0 TP34939_1 TP35570_0
    ## 2         0         2         1         1         1         1         1
    ## 4         0         4         2         2        NA        NA         3
    ##   TP35570_1
    ## 2         1
    ## 4         1

``` r
mydata$likelyGeno_recurrent[,19:26]
```

    ##   TP31810_0 TP31810_1 TP34632_0 TP34632_1 TP34939_0 TP34939_1 TP35570_0
    ## 2         1         1         2         0         2         0         1
    ## 4         2         2         4         0         4         0         2
    ##   TP35570_1
    ## 2         1
    ## 4         2

Estimating genotype probabilities in a diversity panel <a name="diversity"></a>
-------------------------------------------------------------------------------

Pipelines in polyRAD for processing a diversity panel (i.e. a germplasm
collection, a set of samples collected in the wild, or a panel for
genome-wide association analysis or genomic prediction) use iterative
algorithms. Essentially, allele frequencies are re-estimated with each
iteration until convergence is reached.

Here we’ll import a RAD-seq dataset from a large collection of wild and
ornamental *Miscanthus* from Clark *et al.* (2014;
[doi:10.1093/aob/mcu084](http://hdl.handle.net/10.1093/aob/mcu084)).

Since the data are in VCF format, we will need the Bioconductor package
VariantAnnotation to load them. See
<https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html>
for installation instructions.

Again, with your own dataset you will not need to use `system.file` (see
section on mapping populations).

``` r
library(VariantAnnotation)

myVCF <- system.file("extdata", "Msi01genes.vcf", package = "polyRAD")
```

For your own VCF files, you will want to compress and index them before
reading them. This has already been done for the file supplied with
polyRAD, but here is how you would do it:

``` r
mybg <- bgzip(myVCF)
indexTabix(mybg, format = "vcf")
```

Now we can make our `RADdata` object. Because this is a small example
dataset, we are setting `expectedLoci` and `expectedAlleles` to very low
values; in a real dataset they should reflect how much data you are
actually expecting. It is best to slightly overestimate the number of
expected alleles and loci.

``` r
mydata <- VCF2RADdata(myVCF, possiblePloidies = list(2, c(2,2)),
                      expectedLoci = 100, expectedAlleles = 500)
```

    ## Reading file...

    ## Unpacking data from VCF...

    ## Filtering markers...

    ## Phasing 55 SNPs on chromosome 01

    ## Reading file...

    ## 24 loci imported.

    ## Building RADdata object...

    ## Merging rare haplotypes...

    ## 24 markers retained out of 24 originally.

``` r
mydata
```

    ## ## RADdata object ##
    ## 585 taxa and 24 loci
    ## 422433 total reads
    ## Assumed sample cross-contamination rate of 0.001
    ## 
    ## Possible ploidies:
    ## Autodiploid (2)
    ## Allotetraploid (2 2)

For natural populations and diversity panels, we can run
`TestOverdispersion` before performing any genotype calling.

``` r
overdispersionP <- TestOverdispersion(mydata, to_test = 8:10)
```

    ## Genotype estimates not found in object. Performing rough genotype estimation under HWE.

``` r
qq(overdispersionP[["8"]])
```

![](polyRADtutorial_files/figure-markdown_github/unnamed-chunk-23-1.png)

``` r
qq(overdispersionP[["9"]])
```

![](polyRADtutorial_files/figure-markdown_github/unnamed-chunk-23-2.png)

``` r
qq(overdispersionP[["10"]])
```

![](polyRADtutorial_files/figure-markdown_github/unnamed-chunk-23-3.png)

Again, nine looks like a good value.

We can iteratively estimate genotype probabilities assuming
Hardy-Weinberg equilibrium. The argument `tol` is set to a higher value
than the default here in order to help the tutorial run more quickly.
Since *Miscanthus* is highly outcrossing, we will leave the
`selfing.rate` argument at its default of zero.

``` r
mydataHWE <- IterateHWE(mydata, tol = 1e-3, overdispersion = 9)
```

Let’s take a look at allele frequencies:

``` r
hist(mydataHWE$alleleFreq, breaks = 20)
```

![](polyRADtutorial_files/figure-markdown_github/unnamed-chunk-25-1.png)

We can do a different genotype probability estimation that models
population structure and variation in allele frequencies among
populations. We don’t need to specify populations, since principal
components analysis is used to assess population structure assuming an
isolation-by-distance model, with gradients of gene flow across many
groups of individuals. This dataset includes a very broad sampling of
*Miscanthus* across Asia, so it is very appropriate to model population
structure in this case.

For this example, since random number generation is used internally by
`IteratePopStruct` for probabalistic principal components analysis, I am
setting a seed so that the vignette always renders in the same way.

``` r
set.seed(3908)
mydataPopStruct <- IteratePopStruct(mydata, nPcsInit = 8, tol = 5e-03,
                                    overdispersion = 9)
```

Allele frequency estimates have changed slightly:

``` r
hist(mydataPopStruct$alleleFreq, breaks = 20)
```

![](polyRADtutorial_files/figure-markdown_github/unnamed-chunk-27-1.png)

Here’s some of the population structure that was used for modeling
allele frequencies (fairly weak in this case because so few markers were
used):

``` r
plot(mydataPopStruct)
```

![](polyRADtutorial_files/figure-markdown_github/unnamed-chunk-28-1.png)

And here’s an example of allele frequency varying across the
environment. Allele frequencies were estimated for each taxon, and are
stored in the `$alleleFreqByTaxa` slot. In the plot below, color
indicates estimated local allele frequency.

``` r
myallele <- 1
freqcol <- heat.colors(101)[round(mydataPopStruct$alleleFreqByTaxa[,myallele] * 100) + 1]
plot(mydataPopStruct, pch = 21, bg = freqcol)
```

![](polyRADtutorial_files/figure-markdown_github/unnamed-chunk-29-1.png)

As before, we can export the posterior mean genotypes for downstream
analysis.

``` r
wmgenoPopStruct <- GetWeightedMeanGenotypes(mydataPopStruct)
wmgenoPopStruct[1:10,1:5]
```

    ##                       S01_138045_G S01_139820_TT S01_139820_CT
    ## KMS207-8              2.121035e-06    0.07826646    0.18363045
    ## JM0051.003            1.960683e-01    0.31724781    0.20377208
    ## JM0034.001            1.000000e-04    0.00010000    0.12903134
    ## JM0220.001            1.000000e-04    0.50485295    0.20647566
    ## NC-2010-003-001       4.238161e-01    0.00010000    0.00010000
    ## JM0026.001            1.267297e-01    0.00010000    0.17797905
    ## JM0026.002            8.283588e-04    0.07539674    0.20510083
    ## PI294605-US64-0007-01 9.597573e-07    0.01082588    0.07149441
    ## JM0058.001            6.569028e-01    0.00010000    0.23944415
    ## UI10-00086-Silberfeil 9.307387e-01    0.00010000    0.12717438
    ##                       S01_139883_G S01_150928_GG
    ## KMS207-8              7.169509e-02  3.800123e-05
    ## JM0051.003            5.730691e-03  3.169649e-06
    ## JM0034.001            5.005010e-03  2.022024e-03
    ## JM0220.001            1.616635e-01  1.219048e-06
    ## NC-2010-003-001       2.764331e-05  1.351630e-01
    ## JM0026.001            1.000000e-04  1.182203e-04
    ## JM0026.002            7.560629e-03  1.128420e-01
    ## PI294605-US64-0007-01 8.961110e-02  2.685860e-05
    ## JM0058.001            1.000000e-04  3.471470e-03
    ## UI10-00086-Silberfeil 3.114721e-03  5.711384e-06

If you expect that your species has high linkage disequilibrium, the
functions `IterateHWE_LD` and `IteratePopStructLD` behave like
`IterateHWE` and `IteratePopStruct`, respectively, but also update
priors based on genotypes at linked loci.

Considerations for RAM and processing time <a name="considerations"></a>
------------------------------------------------------------------------

`RADdata` objects contain large matrices and arrays for storing read
depth and the parameters that are estimated by the pipeline functions,
and as a result require a lot of RAM (computer memory) in comparison to
the posterior mean genotypes that are exported. A `RADdata` object that
has just been imported will take up less RAM than one that has been
processed by a pipeline function. `RADdata` objects will also take up
more RAM (and take longer for pipeline functions to process) if they
have more possible ploidies and/or higher ploidies.

If you have hundreds of thousands, or possibly even tens of thousands,
of markers in your dataset, it may be too large to process as one object
on a typical computer. In that case, I recommend using the
`SplitByChromosome` function immediately after import. This function
will create separate `RADdata` objects by chromosomes or groups of
chromosomes, and will save those objects to separate R workspace
(.RData) files on your hard drive. You can then run a loop to re-import
those objects one at a time, process each one with a pipeline function,
and export posterior mean geneotypes (or any other parameters you wish
to keep) to a file or a smaller R object. If you have access to a high
performance computing cluster, you may instead wish to process
individual chromosomes as parallel jobs.

If you don’t have alignment positions for your markers, or if you want
to divide them up some other way than by chromosome, see
`SubsetByLocus`. If you are importing from VCF but don’t want to import
the whole genome at once, see the examples on the help page for
`VCF2RADdata` for how to import just a particular genomic region.

You might use `SubsetByLocus` and select a random subset of ~1000 loci
to use with `TestOverdispersion` for estimating the overdispersion
parameter.

If you are using one of the iterative pipelines, it is possible to set
the `tol` argument higher in order to reduce processing time at the
expense of accuracy.

After you have run a pipeline, if you want to keep the `RADdata` object
but discard any components that are not needed for genotype export, you
can use the `StripDown` function.

Below is an example script showing how I processed a real dataset with
hundreds of thousands of SNPs. Note that the (very large) VCF files are
not included with the polyRAD installation.

``` r
library(polyRAD)
library(VariantAnnotation)

# Two files produced by the TASSEL-GBSv2 pipeline using two different
# enzyme systems.
NsiI_file <- "170705Msi_NsiI_genotypes.vcf.bgz"
PstI_file <- "170608Msi_PstI_genotypes.vcf.bgz"

# The vector allSam was defined outside of this script, and contains the 
# names of all samples that I wanted to import.  Below I find sample names
# within the VCF files that match those samples.
NsiI_sam <- allSam[allSam %in% samples(scanVcfHeader(NsiI_file))]
PstI_sam <- allSam[allSam %in% samples(scanVcfHeader(PstI_file))]

# Import two RADdata objects, assuming diploidy.  A large yield size was
# used due to the computer having 64 Gb RAM; on a typical laptop you
# would probably want to keep the default of 5000.
PstI_RAD <- VCF2RADdata(PstI_file, samples = PstI_sam, yieldSize = 5e4,
                        expectedAlleles = 1e6, expectedLoci = 2e5)
NsiI_RAD <- VCF2RADdata(NsiI_file, samples = NsiI_sam, yieldSize = 5e4,
                        expectedAlleles = 1e6, expectedLoci = 2e5)

# remove any loci duplicated across the two sets
nLoci(PstI_RAD)    # 116757
nLoci(NsiI_RAD)    # 187434
nAlleles(PstI_RAD) # 478210
nAlleles(NsiI_RAD) # 952511
NsiI_keeploci <- which(!GetLoci(NsiI_RAD) %in% GetLoci(PstI_RAD))
cat(nLoci(NsiI_RAD) - length(NsiI_keeploci), 
    file = "180522Num_duplicate_loci.txt") #992 duplicate
NsiI_RAD <- SubsetByLocus(NsiI_RAD, NsiI_keeploci)

# combine allele depth into one matrix
PstI_depth <- PstI_RAD$alleleDepth
NsiI_depth <- NsiI_RAD$alleleDepth
total_depth <- matrix(0L, nrow = length(allSam), 
                      ncol = ncol(PstI_depth) + ncol(NsiI_depth),
                      dimnames = list(allSam, 
                                      c(colnames(PstI_depth), 
                                        colnames(NsiI_depth))))
total_depth[,colnames(PstI_depth)] <- PstI_depth[allSam,]
total_depth[rownames(NsiI_depth),colnames(NsiI_depth)] <- NsiI_depth

# combine other slots
total_alleles2loc <- c(PstI_RAD$alleles2loc, 
                       NsiI_RAD$alleles2loc + nLoci(PstI_RAD))
total_locTable <- rbind(PstI_RAD$locTable, NsiI_RAD$locTable)
total_alleleNucleotides <- c(PstI_RAD$alleleNucleotides, 
                             NsiI_RAD$alleleNucleotides)

# build new RADdata object and save
total_RAD <- RADdata(total_depth, total_alleles2loc, total_locTable,
                     list(2L), 0.001, total_alleleNucleotides)
#save(total_RAD, file = "180524_RADdata_NsiIPstI.RData")

# Make groups representing pairs of chromosomes, and one group for all 
# non-assembled scaffolds.
splitlist <- list(c("^01$", "^02$"),
                  c("^03$", "^04$"),
                  c("^05$", "^06$"),
                  c("^07$", "^08$"),
                  c("^09$", "^10$"),
                  c("^11$", "^12$"),
                  c("^13$", "^14$", "^15$"),
                  c("^16$", "^17$"),
                  c("^18$", "^194"), "^SCAFFOLD")
# split by chromosome and save seperate objects
SplitByChromosome(total_RAD, chromlist = splitlist, 
                  chromlist.use.regex = TRUE, fileprefix = "180524splitRAD")

# files with RADdata objects
splitfiles <- grep("^180524splitRAD", list.files("."), value = TRUE)

# list to hold markers formatted for GAPIT/FarmCPU
GAPITlist <- list()
length(GAPITlist) <- length(splitfiles)

# loop through RADdata objects
for(i in 1:length(splitfiles)){
  load(splitfiles[i])
  splitRADdata <- IteratePopStructLD(splitRADdata)
  GAPITlist[[i]] <- ExportGAPIT(splitRADdata)
}
#save(GAPITlist, file = "180524GAPITlist.RData")

# put together into one dataset for FarmCPU
GM.all <- rbind(GAPITlist[[1]]$GM, GAPITlist[[2]]$GM, GAPITlist[[3]]$GM,
                GAPITlist[[4]]$GM, GAPITlist[[5]]$GM, GAPITlist[[6]]$GM, 
                GAPITlist[[7]]$GM, GAPITlist[[8]]$GM,
                GAPITlist[[9]]$GM, GAPITlist[[10]]$GM)
GD.all <- cbind(GAPITlist[[1]]$GD, GAPITlist[[2]]$GD[,-1],
                GAPITlist[[3]]$GD[,-1], GAPITlist[[4]]$GD[,-1],
                GAPITlist[[5]]$GD[,-1], GAPITlist[[6]]$GD[,-1],
                GAPITlist[[7]]$GD[,-1], GAPITlist[[8]]$GD[,-1], 
                GAPITlist[[9]]$GD[,-1], GAPITlist[[10]]$GD[,-1])
#save(GD.all, GM.all, file = "180525GM_GD_all_polyRAD.RData") # 1076888 markers
```

Citing polyRAD <a name="citation"></a>
--------------------------------------

A preprint is available describing polyRAD:

Clark LV, Lipka AE, and Sacks EJ (2018) polyRAD: Genotype calling with
uncertainty from sequencing data in polyploids and diploids. bioRxiv,
[doi:10.1101/380899](https://doi.org/10.1101/380899).
