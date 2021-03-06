# Function to help user identify the best overdispersion parameter
TestOverdispersion <- function(object, ...){
  UseMethod("TestOverdispersion", object)
}
TestOverdispersion.RADdata <- function(object, to_test = seq(6, 20, by = 2),
                                       ...){
  # rough genotype estimation if not done already
  if(is.null(object$posteriorProb)){
    message("Genotype estimates not found in object. Performing rough genotype estimation under HWE.")
    object <- AddAlleleFreqHWE(object)
    object <- AddGenotypePriorProb_HWE(object)
    object <- AddGenotypeLikelihood(object)
    object <- AddPloidyChiSq(object)
    object <- AddGenotypePosteriorProb(object)
  }
  # get the best ploidy for each marker
  bp <- BestPloidies(object$ploidyChiSq)
  # find which was most commonly the best ploidy
  bpt <- table(bp)
  pldindex <- as.integer(names(bpt)[which.max(bpt)])[1]
  theseal <- which(bp == pldindex)
  # identify genotypes that probably have one allele copy
  one_copy_gen <- which(object$posteriorProb[[pldindex]][2,,theseal] >= 0.95)
  if(length(one_copy_gen) == 0){
    stop("Genotype quality too low for this protocol.")
  }
  # get rough probability of sampling allele from this gen (1/ploidy)
  samprob <- 1 / sum(object$priorProbPloidies[[pldindex]])
  
  # get allele counts for these genotypes
  alcnt <- object$alleleDepth[,theseal][one_copy_gen]
  # get total read counts for thise genotypes
  totcnt <- alcnt + object$antiAlleleDepth[,theseal][one_copy_gen]
  # set up values where we need distribution fn
  testcnt <- mapply(function(n, k){
    expected <- n * samprob
    if(k == expected){
      return(0:n)
    } else {
      if(k > expected){
        hi <- k
        lo <- floor(2 * expected - k)
        out <- hi:n
        if(lo >= 0){
          out <- c(0:lo, out)
        }
      } else {
        lo <- k
        hi <- ceiling(2 * expected - k)
        out <- 0:lo
        if(hi >= lo + 1){
          out <- c(out, hi:n)
        }
      }
      return(out)
    }
  }, totcnt, alcnt)
  # set up sampling permutations (except for alcnt == totcnt/2)
  perm <- mapply(lchoose, totcnt, testcnt)

  # set up list for p-value output
  outlist <- list()
  length(outlist) <- length(to_test)
  names(outlist) <- as.character(to_test)
  
  # loop through possible overdispersion parameters
  for(i in 1:length(to_test)){
    op <- to_test[i] * samprob # overdispersion parameter times sampling prob.
    np <- to_test[i] * (1 - samprob)
    outlist[[i]] <- rep(NA_real_, length(totcnt)) # vector for p-values
    for(g in 1:length(totcnt)){ # loop through genotypes
      outlist[[i]][g] <- sum(exp(perm[[g]] +
                                   lbeta(testcnt[[g]] + op, 
                                         totcnt[g] - testcnt[[g]] + np) - 
                                   lbeta(op, np)))
    }
  }
  
  return(outlist)
}
