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

# Function to check whether all individuals would have alleles from all isoloci
# if a locus were split a certain way.  For individuals missing alleles for
# some isoloci, returns the likelihood of those alleles being missing due to
# sampling error.

# Function to cluster alleles into putative isoloci based on sequence similarities.
# Use negative associations as well?
ClusterAlleles <- function(depth, nucleotides, nclust){
  
  nAl <- length(nucleotides) # number of alleles
  if(ncol(depth) != nAl){
    stop("Number of alleles does not match between depth and nucleotides")
  }
  # get distances between sequences
  if(requireNamespace("Biostrings", quietly = TRUE)){
    nucdist <- -Biostrings::stringDist(nucleotides, method = "substitutionMatrix",
                                      substitutionMatrix = -polyRADsubmat)
  } else {
    nucsplit <- strsplit(nucleotides, "")
    nucmat <- matrix(0L, nrow = nAl, ncol = nAl)
    for(i in 1:(nAl-1)){
      for(j in (i+1):nAl){
        nucmat[i,j] <- nucmat[j,i] <- sum(sapply(1:length(splitnucA),
                                                 function(i) polyRADsubmat[nucsplit[[i]],
                                                                           nucsplit[[j]]]))
      }
    }
    nucdist <- as.dist(nucmat)
  }
  
  # make clusters based on DNA sequences
  nucclust <- cutree(hclust(nucdist), k = nclust)
}
