# wrapper function to estimate Hind/He for each locus in a RADdata object.
HindHe <- function(object, ...){
  UseMethod("HindHe", object)
}
HindHe.RADdata <- function(object, omitTaxa = GetBlankTaxa(object), ...){
  taxa <- GetTaxa(object)[!GetTaxa(object) %in% omitTaxa]
  
  hindhe <- HindHeMat(object$alleleDepth[taxa,, drop = FALSE],
                      object$depthRatio[taxa,, drop = FALSE],
                      object$alleles2loc, nLoci(object), numeric(0))
  colnames(hindhe) <- GetLoci(object)
  rownames(hindhe) <- taxa
  return(hindhe)
}

# Function to get Hind/He matrix in a mapping population
HindHeMapping <- function(object, ...){
  UseMethod("HindHeMapping", object)
}
HindHeMapping.RADdata <- function(object, n.gen.backcrossing = 0,
                                  n.gen.intermating = 0, n.gen.selfing = 0,
                                  ploidy = object$possiblePloidies[[1]],
                                  minLikelihoodRatio = 10,
                                  omitTaxa = c(GetDonorParent(object), 
                                               GetRecurrentParent(object), 
                                               GetBlankTaxa(object)), ...){
  if(length(ploidy) != 1){
    stop("Current implementation for autopolyploids only")
  }
  if(n.gen.intermating > 0){
    stop("If the most recent generation was produced by random mating among progeny, use HindHe instead.")
  }
  donorParent <- GetDonorParent(object)
  recurrentParent <- GetRecurrentParent(object)
  progeny <- GetTaxa(object)[!GetTaxa(object) %in% omitTaxa]
  object <- EstimateParentalGenotypes(object, donorParent = donorParent,
                                      recurrentParent = recurrentParent,
                                      n.gen.backcrossing = n.gen.backcrossing,
                                      n.gen.selfing = n.gen.selfing,
                                      n.gen.intermating = n.gen.intermating,
                                      minLikelihoodRatio = minLikelihoodRatio,
                                      donorParentPloidies = list(ploidy),
                                      recurrentParentPloidies = list(ploidy))
  likelyGenDon <- object$likelyGeno_donor[as.character(ploidy),]
  likelyGenRec <- object$likelyGeno_recurrent[as.character(ploidy),]
  # Get probabilities of pairs of alleles from a random progeny coming from
  # different locus copies one parent or the other, or from different parents.
  progAlProbs <- .progAlProbs(ploidy, n.gen.backcrossing, n.gen.selfing)
  
  # Identify loci where multiallelic genotypes can be determined
  goodLocDon <- tapply(likelyGenDon, object$alleles2loc,
                       function(x) !any(is.na(x)) && sum(x) == ploidy)
  goodLocRec <- tapply(likelyGenRec, object$alleles2loc,
                       function(x) !any(is.na(x)) && sum(x) == ploidy)
  keeploc <- which(goodLocDon & goodLocRec)
  object <- SubsetByLocus(object, keeploc)
  
  # Get within- and across- parent probabilties of sampling two different alleles.
  parentHo <- matrix(c(HoOneParent(likelyGenRec, object$alleles2loc, keeploc, ploidy),
                       HoOneParent(likelyGenDon, object$alleles2loc, keeploc, ploidy),
                       HoTwoParents(likelyGenRec, likelyGenDon, object$alleles2loc,
                                    keeploc, ploidy)),
                     byrow = TRUE, ncol = length(keeploc), nrow = 3)
  # Get per-locus 'He' values, which indicate the probability of two alleles sampled
  # without replacement from one progeny being different.
  heByLoc <- progAlProbs %*% parentHo
  
  # Get Hind/He matrix
  outmat <- HindHeMat(object$alleleDepth[progeny,, drop = FALSE],
                      object$depthRatio[progeny,, drop = FALSE],
                      object$alleles2loc, nLoci(object), heByLoc)
  outmat[is.infinite(outmat)] <- NaN
  colnames(outmat) <- GetLoci(object)
  rownames(outmat) <- progeny
  
  return(outmat)
}
