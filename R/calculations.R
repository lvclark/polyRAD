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

# internal function to multiply genotype priors by genotype likelihood
.priorTimesLikelihood <- function(object){
  if(!"RADdata" %in% class(object)){
    stop("RADdata object needed for .priorTimesLikelihood")
  }
  if(is.null(object$priorProb)){
    stop("Genotype prior probabilities must be added first.")
  }
  if(is.null(object$genotypeLikelihood)){
    stop("Genotype likelihoods must be added first.")
  }
  
  ploidytotpriors <- sapply(object$priorProb, function(x) dim(x)[1] - 1)
  ploidytotlikeli <- sapply(object$genotypeLikelihood, function(x) dim(x)[1] - 1)
  
  results <- list()
  length(results) <- length(object$priorProb)
  
  for(i in 1:length(object$priorProb)){
    j <- which(ploidytotlikeli == ploidytotpriors[i])
    stopifnot(length(j) == 1)
    if(attr(object, "priorType") == "population"){
      # expand priors out by individuals
      thispriorarr <- array(object$priorProb[[i]], 
                            dim = c(dim(object$priorProb[[i]])[1], 1, 
                                    dim(object$priorProb[[i]])[2]))[,rep(1, length(GetTaxa(object))),]
      dimnames(thispriorarr) <- dimnames(object$genotypeLikelihoods)[[j]]
    } else {
      thispriorarr <- object$priorProb[[i]]
    }
    stopifnot(identical(dim(thispriorarr), dim(object$genotypeLikelihood[[j]])))
    results[[i]] <- thispriorarr * object$genotypeLikelihood[[j]]
  }
  
  return(results)
}
