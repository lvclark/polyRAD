# Wrapper functions that run multiple steps in a pipeline for polyRAD

IterateHWE <- function(object, tol = 1e-8, excludeTaxa = GetBlankTaxa(object)){
  if(!"RADdata" %in% class(object)){
    stop("RADdata object needed.")
  }
  if(tol > 1e-2){
    warning("tol unexpectedly high.")
  }
  if(tol < 0){
    stop("tol must be above zero")
  }
  nIter <- 1 # which round of iteration are we on
  meanDiff <- 1 # mean difference between allele frequencies from round to round
                # (1 is a dummy value for while loop)
  
  object <- AddAlleleFreqHWE(object, excludeTaxa = excludeTaxa)
  
  while(meanDiff > tol){
    cat(paste("Starting iteration", nIter), sep = "\n")
    oldAlFreq <- object$alleleFreq
    object <- AddGenotypePriorProb_HWE(object)
    object <- AddGenotypeLikelihood(object)
    object <- AddPloidyChiSq(object, excludeTaxa = excludeTaxa)
    object <- AddGenotypePosteriorProb(object)
    object <- AddAlleleFreqHWE(object, excludeTaxa = excludeTaxa)
    nIter <- nIter + 1
    meanDiff <- mean(abs(oldAlFreq - object$alleleFreq), na.rm = TRUE)
    cat(paste("Mean difference in allele frequencies of", meanDiff), sep = "\n")
  }
  
  return(object)
}

IteratePopStruct <- function(object, tol = 1e-2, 
                             excludeTaxa = GetBlankTaxa(object),
                             nPcsInit = 50){
  if(!"RADdata" %in% class(object)){
    stop("RADdata object needed.")
  }
  if(tol > 1e-2){
    warning("tol unexpectedly high.")
  }
  if(tol < 0){
    stop("tol must be above zero")
  }
  
  nIter <- 1 # which round of iteration are we on
  meanDiff <- 1 # mean difference between allele frequencies from round to round
                # (1 is a dummy value for while loop)
  
  object <- AddPCA(object, nPcsInit = nPcsInit)
  object <- AddAlleleFreqByTaxa(object)
  
  while(meanDiff > tol){
    cat(paste("Starting iteration", nIter), sep = "\n")
    oldAlFreq <- object$alleleFreqByTaxa
    object <- AddAlleleFreqHWE(object, excludeTaxa = excludeTaxa)
    object <- AddGenotypePriorProb_ByTaxa(object)
    object <- AddGenotypeLikelihood(object)
    object <- AddGenotypePosteriorProb(object)
    # later insert object <- AddPloidyChiSq(object, excludeTaxa = excludeTaxa)
    object <- AddPCA(object, nPcsInit = nPcsInit)
    object <- AddAlleleFreqByTaxa(object)
    nIter <- nIter + 1
    meanDiff <- mean(abs(oldAlFreq - object$alleleFreqByTaxa), na.rm = TRUE)
    cat(paste("Mean difference in allele frequencies of", meanDiff), sep = "\n")
  }
  
  return(object)
}
