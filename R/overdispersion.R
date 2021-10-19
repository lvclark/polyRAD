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
  one_copy_gen <-
    sapply(object$posteriorProb[pldindex,],
           function(x){
             which(x[2,,theseal] >= 0.95)
           }, simplify = FALSE)
  if(all(lengths(one_copy_gen) == 0)){
    stop("Genotype quality too low for this protocol.")
  }
  # get rough probability of sampling allele from this gen (1/ploidy)
  samprobs <- 1 / sum(object$priorProbPloidies[[pldindex]]) /
    as.integer(colnames(object$priorProb)) * 2
  
  # set up list for p-value output
  outlist <- list(numeric(0))[rep(1, length(to_test))]
  names(outlist) <- as.character(to_test)
  
  for(h in seq_len(ncol(object$priorProb))){
    # samples with this ploidy
    thesesam <- dimnames(object$posteriorProb[[pldindex,h]])[[2]]
    # get allele counts for these genotypes
    alcnt <- object$alleleDepth[thesesam,theseal][one_copy_gen[[h]]]
    # get total read counts for these genotypes
    totcnt <- alcnt + object$antiAlleleDepth[thesesam,theseal][one_copy_gen[[h]]]
    # set up values where we need distribution fn
    testcnt <- mapply(function(n, k){
      expected <- n * samprobs[h]
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
    
    # loop through possible overdispersion parameters
    for(i in 1:length(to_test)){
      op <- to_test[i] * samprobs[h] # overdispersion parameter times sampling prob.
      np <- to_test[i] * (1 - samprobs[h])
      temp <- rep(NA_real_, length(totcnt)) # vector for p-values
      for(g in 1:length(totcnt)){ # loop through genotypes
        temp[g] <- sum(exp(perm[[g]] +
                             lbeta(testcnt[[g]] + op, 
                                   totcnt[g] - testcnt[[g]] + np) - 
                             lbeta(op, np)))
      }
      outlist[[i]] <- c(outlist[[i]], temp)
    }
  }
  
  # print out suggestion for best value to use; determine which follows
  # expected values most closely.
  expP <- ppoints(length(outlist[[1]]))
  odscores <- sapply(outlist,
                     function(x){
                       sqrt(mean((sort(x) - expP) ^ 2))
                     })
  best <- as.numeric(names(odscores)[which.min(odscores)])
  cat(paste0("Optimal value is ", best, "."), sep = "\n")
  if(best == min(to_test)){
    cat("Consider testing lower values.", sep = "\n")
  }
  if(best == max(to_test)){
    cat("Consider testing higher values.", sep = "\n")
  }
  
  return(outlist)
}
