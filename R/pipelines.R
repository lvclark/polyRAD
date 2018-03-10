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
    message(paste("Starting iteration", nIter))
    oldAlFreq <- object$alleleFreq
    object <- AddGenotypePriorProb_HWE(object)
    object <- AddGenotypeLikelihood(object)
    object <- AddPloidyChiSq(object, excludeTaxa = excludeTaxa)
    object <- AddGenotypePosteriorProb(object)
    object <- AddAlleleFreqHWE(object, excludeTaxa = excludeTaxa)
    nIter <- nIter + 1
    meanDiff <- mean(abs(oldAlFreq - object$alleleFreq), na.rm = TRUE)
    message(paste("Mean difference in allele frequencies of", meanDiff))
  }
  
  return(object)
}

IteratePopStruct <- function(object, tol = 1e-3, 
                             excludeTaxa = GetBlankTaxa(object),
                             nPcsInit = 50, minfreq = 0.0001){
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
  
  message("Performing initial PCA and allele frequency estimation.")
  object <- AddPCA(object, nPcsInit = nPcsInit)
  object <- AddAlleleFreqByTaxa(object)
  
  while(meanDiff > tol){
    message(paste("Starting iteration", nIter))
    message(paste("PCs used:", dim(object$PCA)[2]))
    oldAlFreq <- object$alleleFreqByTaxa
    object <- AddAlleleFreqHWE(object, excludeTaxa = excludeTaxa)
    object <- AddGenotypePriorProb_ByTaxa(object)
    object <- AddGenotypeLikelihood(object)
    object <- AddPloidyChiSq(object, excludeTaxa = excludeTaxa)
    object <- AddGenotypePosteriorProb(object)
    object <- AddPCA(object, nPcsInit = dim(object$PCA)[2] + 1,
                     minPcsOut = dim(object$PCA)[2])
      # -> reasoning for PC number constraints: 
      #     tend to get more accuracy with more PCs,
      #     number of PCs going up and down prevents convergence of algorithm
    object <- AddAlleleFreqByTaxa(object, minfreq = minfreq)
    nIter <- nIter + 1
    meanDiff <- mean(abs(oldAlFreq - object$alleleFreqByTaxa), na.rm = TRUE)
    message(paste("Mean difference in allele frequencies of", meanDiff))
  }
  
  return(object)
}

# Function to run genotype estimation using both population structure and
# linkage disequilibrium.
# LDdist is the distance in nucleotides within which to search for alleles in LD.
# minLDcorr is the minimum correlation coefficient between the residuals of 
# allelic values after regression on PCA axes, with the values of a nearby 
# allele, for that allele to be used in genotype prediction.
IteratePopStructLD <- function(object, tol = 1e-3, 
                             excludeTaxa = GetBlankTaxa(object),
                             nPcsInit = 50, minfreq = 0.0001,
                             LDdist = 1e4, minLDcor = 0.01){
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
  
  # Initialization before testing for LD
  message("Performing initial PCA and allele frequency estimation.")
  object <- AddPCA(object, nPcsInit = nPcsInit)
  object <- AddAlleleFreqByTaxa(object)
  message("Performing preliminary genotype estimation")
  object <- AddAlleleFreqHWE(object, excludeTaxa = excludeTaxa)
  object <- AddGenotypePriorProb_ByTaxa(object)
  object <- AddGenotypeLikelihood(object)
  object <- AddPloidyChiSq(object, excludeTaxa = excludeTaxa)
  object <- AddGenotypePosteriorProb(object)
  object <- AddPCA(object, nPcsInit = dim(object$PCA)[2] + 1,
                   minPcsOut = dim(object$PCA)[2])
  object <- AddAlleleFreqByTaxa(object, minfreq = minfreq)
  
  # Test for LD
  message("Finding alleles in LD")
  wmgeno <- GetWeightedMeanGenotypes(object, omit1allelePerLocus = FALSE)
  object$alleleLinkages <- list() # to hold LD information for each allele
  length(object$alleleLinkages) <- nAlleles(object)
  for(L in 1:attr(object, "nLoci")){
    theseAlleles <- which(object$alleles2loc == L) # alleles for this locus
    nearbyAlleles <- FindNearbyAlleles(object, L, LDdist) # alleles to test for LD
    for(a in theseAlleles){
      # get residuals of allelic values not predicted by population structure
      thisLM <- lm(wmgeno[,a] ~ object$PCA)
      thisResid <- residuals(thisLM)
      # set up vectors to store correlated alleles and correlation stats
      linkedalleles <- integer(0)
      correlations <- numeric(0)
      # see how each nearby allele correlates with residuals
      for(a2 in nearbyAlleles){
        thiscor <- cor(wmgeno[,a2], thisResid)
        if(thiscor >= minLDcor){
          linkedalleles <- c(linkedalleles, a2)
          correlations <- c(correlations, thiscor)
        }
      }
      # store linkages for this allele
      object$alleleLinkages[[a]] <- data.frame(allele = linkedalleles,
                                        corr = correlations)
    }
  }
  
  # Iteratively estimate genotypes
  while(meanDiff > tol){
    message(paste("Starting iteration", nIter))
    message(paste("PCs used:", dim(object$PCA)[2]))
    oldAlFreq <- object$alleleFreqByTaxa
    object <- AddAlleleFreqHWE(object, excludeTaxa = excludeTaxa)
    object <- AddGenotypePriorProb_ByTaxa(object)
    ## Add a call here to get prior probs with LD
    object <- AddGenotypeLikelihood(object)
    object <- AddPloidyChiSq(object, excludeTaxa = excludeTaxa)
    object <- AddGenotypePosteriorProb(object)
    object <- AddPCA(object, nPcsInit = dim(object$PCA)[2] + 1,
                     minPcsOut = dim(object$PCA)[2])
    # -> reasoning for PC number constraints: 
    #     tend to get more accuracy with more PCs,
    #     number of PCs going up and down prevents convergence of algorithm
    object <- AddAlleleFreqByTaxa(object, minfreq = minfreq)
    nIter <- nIter + 1
    meanDiff <- mean(abs(oldAlFreq - object$alleleFreqByTaxa), na.rm = TRUE)
    message(paste("Mean difference in allele frequencies of", meanDiff))
  }
  
  return(object)
}

PipelineMapping2Parents <- function(object, donorParent = GetDonorParent(object),
                                    recurrentParent = GetRecurrentParent(object),
                                    n.gen.backcrossing = 0,
                                    n.gen.intermating = 0,
                                    n.gen.selfing = 0, 
                                    donorParentPloidies = object$possiblePloidies,
                                    recurrentParentPloidies = object$possiblePloidies,
                                    minLikelihoodRatio = 10,
                                    freqAllowedDeviation = 0.05,
                                    freqExcludeTaxa = c(GetDonorParent(object),
                                                    GetRecurrentParent(object),
                                                    GetBlankTaxa(object))){
  # estimate possible allele frequencies
  allelesin <- max(sapply(donorParentPloidies, sum)) + 
    max(sapply(recurrentParentPloidies, sum))
  possfreq <- seq(0, 1, length.out = (n.gen.backcrossing + 1) * allelesin + 1)
  
  object <- AddAlleleFreqMapping(object, expectedFreqs = possfreq,
                                 allowedDeviation = freqAllowedDeviation,
                                 excludeTaxa = freqExcludeTaxa)
  # calculations for rest of pipeline
  object <- AddGenotypeLikelihood(object)
  object <- AddGenotypePriorProb_Mapping2Parents(object, donorParent = donorParent,
                                                 recurrentParent = recurrentParent,
                                                 n.gen.backcrossing = n.gen.backcrossing,
                                                 n.gen.intermating = n.gen.intermating,
                                                 n.gen.selfing = n.gen.selfing,
                                                 donorParentPloidies = donorParentPloidies,
                                                 recurrentParentPloidies = recurrentParentPloidies,
                                                 minLikelihoodRatio = minLikelihoodRatio)
  object <- AddPloidyChiSq(object, excludeTaxa = freqExcludeTaxa)
  object <- AddGenotypePosteriorProb(object)
  
  return(object)
}
