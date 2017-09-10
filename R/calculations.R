# additional polyRAD functions that perform calculations

# internal function to take allele frequencies and get prior probs under HWE
# freqs is a vector of allele frequencies
# ploidy is a vector indicating the ploidy
.HWEpriors <- function(freqs, ploidy){
  if(length(unique(ploidy)) != 1){
    stop("All subgenomes must be same ploidy")
  }
  nsubgen <- length(ploidy)
  if(nsubgen == 1){ # for diploid/autopolyploid, or single subgenome with recursion
    priors <- matrix(NA, nrow = ploidy+1, ncol = length(freqs),
                     dimnames = list(as.character(0:ploidy), names(freqs)))
    antifreqs <- 1 - freqs
    for(i in 0:ploidy){
      priors[i+1,] <- choose(ploidy, i) * freqs ^ i * antifreqs ^ (ploidy - i)
    }
  } else {
    remainingfreqs <- freqs
    priors <- matrix(1, nrow = 1, ncol = length(freqs))
    for(pld in ploidy){
      thesefreqs <- remainingfreqs
      # allele frequencies partitioned for this subgenome
      thesefreqs[thesefreqs > 1/nsubgen] <- 1/nsubgen
      # allele frequencies reserved for remaining subgenomes
      remainingfreqs <- remainingfreqs - thesefreqs
      # priors just for this subgenome
      thesepriors <- .HWEpriors(nsubgen * thesefreqs, pld)
      
      # multiply by priors already calculated to get overall priors
      oldpriors <- priors
      newpld <- dim(thesepriors)[1] + dim(oldpriors)[1] - 2
      priors <- matrix(0, nrow = newpld + 1, ncol = length(freqs), 
                       dimnames = list(as.character(0:newpld),
                                       names(freqs)))
      for(i in 1:dim(oldpriors)[1]){
        for(j in 1:dim(thesepriors)[1]){
          thisalcopy <- i + j - 2
          priors[thisalcopy + 1,] <- priors[thisalcopy + 1,] + 
            oldpriors[i,] * thesepriors[j,]
        }
      }
    }
  }
  return(priors)
}
