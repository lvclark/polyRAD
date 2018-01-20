## Functions for seeing if read counts match inheretance modes and sorting
## alleles into isoloci.

# Function to test if there are more alleles than expected for many samples
# for a given ploidy.
# depth is a read depth matrix for one locus, with alleles in columns and taxa
# in rows.
# ploidy is a single integer indicating the max number of alleles expected per
# genotype.
# contamRate is a number close to zero indicating the sample cross-contamination
# rate.
# freq is the subset of the allele frequency vector corresponding to this locus.
# genotypeLikelihoods is a 3d array, if present, from the RADdata object.
CheckAllelesPerLocus <- function(depth, ploidy, contamRate, freq = NULL,
                                 genotypeLikelihoods = NULL){
  
  # Internal function to sum likelihoods for alleles having zero copies.
  # Returns the likelihood that extra alleles are all due to contamination.
  # likelihoods is a vector containing likelihood for the allele being
  # present in zero copies, for all alleles with reads.
  sumLikelihoods <- function(likelihoods){
    # number of additional alleles beyond what is expected
    nExtra <- length(likelihoods) - ploidy
    if(nExtra == 0){
      return(1)
    }
    if(nExtra == 1){
      return(sum(likelihoods))
    }
    if(nExtra > 1){
      outlik <- 0
      for(i in 1:length(likelihoods)){
        outlik <- outlik + likelihoods[i] * sumLikelihoods(likelihoods[-i])
      }
      return(outlik/nExtra)
    }
  }
  
  # get total depth for each individual
  locDepth <- rowSums(depth)
  # get allele frequencies
  if(is.null(freq)){
    freq <- colMeans(depth/locDepth, na.rm = TRUE)
  }
  # get probability of sampling each allele due to contamination
  contamProb <- freq * contamRate
  # get identities of individuals with too many alleles
  toomany <- which(rowSums(depth > 0) > ploidy)
  # vector to store likelihood of seeing extra alleles due to contamination
  contamlik <- numeric(length(toomany))
  names(contamlik) <- names(toomany)
  
  if(length(toomany) > 0){
    for(i in 1:length(toomany)){
      taxon <- toomany[i]
      # For each individual, and all alleles that individual has above zero reads,
      # get the probability that all reads of that allele were sampled due to
      # contamination.
      if(is.null(genotypeLikelihoods)){
        thislik <- sapply(which(depth[taxon,] > 0),
                          function(x) dbinom(depth[taxon, x], locDepth[taxon],
                                             contamProb[x]))
      } else {
        thislik <- genotypeLikelihoods[1, taxon, depth[taxon,] > 0]
      }
      # Get probability that in reality there are ploidy or fewer alleles
      contamlik[i] <- sumLikelihoods(thislik)
    }
#    return(prod(contamlik))
  }# else {
#    return(1)
#  }
  return(contamlik) # return one likelihood for each ind. with too many alleles
}
