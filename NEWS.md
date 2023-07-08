# polyRAD 2.0.1

* Bug fix in `TestOverdispersion` when the optimal value is the minimum tested.
* `VCF2RADdata` now can make use of .csi indices on VCFs that are already bgzipped.
* Bug fix in `Export_MAPpoly` relating to how ploidy is coded.

# polyRAD 2.0

* Multiploid populations are now supported. This includes natural populations with
variable ploidy, as well as better support for mapping populations such as
triploid populations derived from diploid x tetraploid crosses. Internally,
slots containing genotype prior probabilities, likelihoods, and posterior
probabilities are now formatted as two-dimensional lists, with possible
inheritance modes across loci in the first dimension and varying ploidy across
taxa in the second dimension, whereas previously these were one-dimensional
lists reflecting inheritance modes across loci only. The `priorProbPloidies` slot
has been removed as it is no longer needed.

* The `ExamineGenotype` function has been added.

* The `AddPriorTimesLikelihood` function has been removed.

* A warning is now returned if the user attempts to simulate self-fertilization
in odd-ploidy individuals (i.e. using the selfing.rate argument in
`IterateHWE` and related functions). The selfing rate is set to zero in these
individuals to prevent an error.

* `IteratePopStruct` and `IteratePopStructLD` now have an argument called
`maxR2changeratio` to allow more fine-tuning of the number of principal components
used.

* The `readDArTag` function has been updated to support a newer format from
Excellence in Breeding, as well as their original format.

* `MergeIdenticalHaplotypes` now takes IUPAC ambiguity codes into account. It is
now used internally by `VCF2RADdata`, `readStacks`, `readTASSELGBSv2`, and
`readProcessIsoloci`. Some additional loci will now be filtered out by these
functions.

* `GetProbableGenotypes` with `multiallelic = "correct"`, and by extension
`RADdata2VCF`, now run much faster by searching a much narrower range of
possible multiallelic genotypes.

* Optimizations to improve speed of data import.

* Bug fix in `MergeRareHaplotypes` when some alleles have zero reads.

# polyRAD 1.6

* `TestOverdispersion` now prints a helpful suggestion for which overdispersion
parameter to use, as well as returning that value in its output.
The use of `qqman` has been eliminated from the vignettes.

* `ExpectedHindHe`, `ExpectedHindHeMapping`, and `SimAlleleDepth` now have a
`contamRate` parameter to simulate sample cross contamination when sampling
allelic read depth, and an `errorRate` parameter to simulate sequencing error.

* A bug has been fixed in `Export_GWASpoly` so that special characters in sample
names are not converted to periods.

# polyRAD 1.5

* `SimGenotypesMapping`, `ExpectedHindHeMapping`, `Export_polymapR_probs`,
`reverseComplement`, and `readDArTag` functions added.

* Bug fixed in `GetProbableGenotypes` when there are multiple `possiblePloidies`
and `multiallelic` is set to `"na"` or `"correct"`.

* Bug fixed in `process_sam_multi.py` for cases where sample names contain spaces.

* Bug fixed in `HindHeMapping` when loci are discarded due to difficulty estimating
parental genotypes.

* Vignettes have been updated with a figure to assist in estimating inbreeding.

# polyRAD 1.4

* `Export_GWASpoly` can now export posterior mean genotypes.

* `Export_adegenet_genind` function added.

# polyRAD 1.3

* Added the `ExpectedHindHe`, `SimGenotypes`, and `SimAlleleDepth` functions.
The first in particular is intended to help the user to select a Hind/He
threshold for filtering markers.

* Added the `Export_Structure` function.

* In the variant calling pipeline, `process_isoloci.py` has been modified so
that the user can set an expected and maximum Hind/He. The previous
functionality, where these thresholds were calculated automatically from
inbreeding, remains available.

* Fixed a bug in `VCF2RADdata` that would occur if `expectedAlleles` or
`expectedLoci` were set too low.  `VCF2RADdata` now also sets the
`"Variable_sites_only"` attribute to be `FALSE` if `phaseSNPs = FALSE`.

* Fixed a bug in `SubsetByLocus`, `MergeRareHaplotypes`, and
`MergeIdenticalHaplotypes` that would delete the `"Variable_sites_only"`
attribute, preventing VCF export.

# polyRAD 1.2

* The functions `HindHe` and `HindHeMapping` have been added, to assist with filtering
loci that are behaving in a non-Mendelian fashion, as well as individuals that
do not behave like the expected ploidy.  The function `InbreedingFromHIndHe` has
been added to facilitate estimating the inbreeding statistic F from the results
of `HindHe`.

* The Python scripts `process_sam_multi.py` and `process_isoloci.py` have been added,
to serve as a pipeline for variant calling in highly duplicated genomes.  These
are described in a new vignette, "Variant and Genotype Calling in Highly
Duplicated Genomes" (`isolocus_sorting.Rmd`).  The functions `readProcessSamMulti`
and `readProcessIsoloci` have been added for import of data from the Python
scripts into polyRAD.

* The function `RADdata2VCF` has been added for export to VCF.

* A Python script has been added to help identify full tag sequences for alleles
imported from TASSEL-GBSv2 using `VCF2RADdata`.

* The `GetProbableGenotypes` function has been updated to optionally recognize
genotypes where copy number does not sum to the ploidy across all alleles for a
locus.  It can now either set these genotypes to NA or correct them to sum to
the ploidy.  Some internal `Rcpp` functions were added for this purpose.

* A bug has been fixed in `SubsetByPloidy`.

* `MergeTaxaDepth` now generates a more useful error if taxa are specified that
do not exist in the object.

* If duplicate locus names are passed to `SubsetByLocus`, it now generates a
warning and deals with the duplicated loci more appropriately.

* The `AddDepthSamplingPermutations` function has been added.  Permutations for
estimating genotype likelihoods are no longer calculated at the time of data
import, but instead are calculated the first time that genotype likelihoods
are calculated.  This is to avoid unnessasary calculation and recalculation
when alleles and taxa are filtered or merged before running genotype calling.
From the user perspective, there isn't any change except faster import.

* The `EstimateParentalGenotypes` function has been added, using code that was
previously in the `AddGenotypePriorProb_Mapping2Parents`.  This does not
impact anything from the user's perspective, but makes the code slightly
less unwieldy and allows parental genotype estimations to be used for other
purposes.

* The `MergeIdenticalHaplotypes` function has been added.

* Internal functions for simulating gametes and calculating their frequencies have
been corrected so that they still work when only one allele is being simulated.

# polyRAD 1.1

* A plot method has been added for `RADdata` objects.

* The `MergeTaxaDepth`, `RemoveUngenotypedLoci`, and `SubsetByPloidy` utility
functions have been added.

* The `Export_MAPpoly` and `Export_GWASpoly` functions have been added.

* Bug fixes have been made in `TestOverdispersion`, `VCF2RADdata`,
`SubsetByLocus`, and `SubsetByTaxon`.

* Some code in `AddAlleleFreqByTaxa` has been translated to `Rcpp` to speed
computation time for `IteratePopStruct` and `IteratePopStructLD`.

# polyRAD 1.0

* Genotype likelihoods are now estimated under a beta-binomial distribution 
rather than the binomial distribution.  This change was made so that real
sequencing data would be accurately modeled; even in diploid heterozygotes,
read depth of two alleles is often very different from a 1:1 ratio, due to
many underlying issues with sequencing data that would be difficult to model.
Under the beta-binomial with respect to the binomial, there is an increased 
probability of read depth ratios that differ from the true allele copy 
ratio.  In a practical sense, this means reduced certainty in the estimation of
allele copy number from read depth alone, and an increased importance of 
genotype prior probabilities.  The exact shape of the beta-binomial 
distribution is determined by an overdispersion parameter, which the user can
optimize using the `TestOverdispersion` function.

* When using linkage disequilibrium to update genotype priors, the square of
Pearson's correlation coefficient is now used for weighting markers, where
Pearson's correlation coefficient was used previously without being squared.
This applies to both mapping populations and diversity panels, and results
in improved genotyping accuracy.

* The functions `Export_polymapR`, `readTASSELGBSv2`, `RemoveHighDepthLoci`,
`AddGenotypePriorProb_Even`, and `TestOverdispersion` have been added.

* This version of polyRAD is incompatible with `RADdata` objects generated by
previous versions of polyRAD due to a change in format of the 
`depthSamplingPermutations` slot.  This slot was changed to simplify the
estimation of genotype likelihood.
