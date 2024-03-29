# Wrapper functions that run multiple steps in a pipeline for polyRAD

IterateHWE <- function(object, selfing.rate = 0, tol = 1e-5, 
                       excludeTaxa = GetBlankTaxa(object),
                       overdispersion = 9){
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
    object <- AddGenotypePriorProb_HWE(object, selfing.rate)
    object <- AddGenotypeLikelihood(object, overdispersion = overdispersion)
    object <- AddPloidyChiSq(object, excludeTaxa = excludeTaxa)
    object <- AddGenotypePosteriorProb(object)
    object <- AddAlleleFreqHWE(object, excludeTaxa = excludeTaxa)
    nIter <- nIter + 1
    meanDiff <- mean(abs(oldAlFreq - object$alleleFreq), na.rm = TRUE)
    message(paste("Mean difference in allele frequencies of", meanDiff))
  }
  
  return(object)
}

IterateHWE_LD <- function(object, selfing.rate = 0, tol = 1e-5, 
                          excludeTaxa = GetBlankTaxa(object),
                          LDdist = 1e4, minLDcorr = 0.2,
                          overdispersion = 9){
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
  
  message("Performing preliminary genotype estimation.")
  object <- AddAlleleFreqHWE(object, excludeTaxa = excludeTaxa)
  object <- AddGenotypePriorProb_HWE(object, selfing.rate)
  object <- AddGenotypeLikelihood(object, overdispersion = overdispersion)
  object <- AddPloidyChiSq(object, excludeTaxa = excludeTaxa)
  object <- AddGenotypePosteriorProb(object)
  object <- AddAlleleFreqHWE(object, excludeTaxa = excludeTaxa)
  message("Finding alleles in LD.")
  object <- AddAlleleLinkages(object, type = "hwe",
                              linkageDist = LDdist, minCorr = minLDcorr,
                              excludeTaxa = GetBlankTaxa(object))
  
  while(meanDiff > tol){
    message(paste("Starting iteration", nIter))
    oldAlFreq <- object$alleleFreq
    object <- AddGenotypePriorProb_HWE(object, selfing.rate)
    object <- AddGenotypePriorProb_LD(object, type = "hwe")
    object <- AddGenotypeLikelihood(object, overdispersion = overdispersion)
    object <- AddPloidyChiSq(object, excludeTaxa = excludeTaxa)
    object <- AddGenotypePosteriorProb(object)
    object <- AddAlleleFreqHWE(object, excludeTaxa = excludeTaxa)
    nIter <- nIter + 1
    meanDiff <- mean(abs(oldAlFreq - object$alleleFreq), na.rm = TRUE)
    message(paste("Mean difference in allele frequencies of", meanDiff))
  }
  
  return(object)
}

IteratePopStruct <- function(object, selfing.rate = 0, tol = 1e-3, 
                             excludeTaxa = GetBlankTaxa(object),
                             nPcsInit = 10, minfreq = 0.0001,
                             overdispersion = 9, maxR2changeratio = 0.05){
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
  object <- AddPCA(object, nPcsInit = nPcsInit, maxR2changeratio = maxR2changeratio)
  object <- AddAlleleFreqByTaxa(object)
  
  while(meanDiff > tol){
    message(paste("Starting iteration", nIter))
    message(paste("PCs used:", dim(object$PCA)[2]))
    oldAlFreq <- object$alleleFreqByTaxa
    object <- AddAlleleFreqHWE(object, excludeTaxa = excludeTaxa)
    object <- AddGenotypePriorProb_ByTaxa(object, selfing.rate)
    object <- AddGenotypeLikelihood(object, overdispersion = overdispersion)
    object <- AddPloidyChiSq(object, excludeTaxa = excludeTaxa)
    object <- AddGenotypePosteriorProb(object)
    object <- AddPCA(object, nPcsInit = dim(object$PCA)[2] + 1,
                     minPcsOut = dim(object$PCA)[2],
                     maxR2changeratio = maxR2changeratio)
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
IteratePopStructLD <- function(object, selfing.rate = 0, tol = 1e-3, 
                             excludeTaxa = GetBlankTaxa(object),
                             nPcsInit = 10, minfreq = 0.0001,
                             LDdist = 1e4, minLDcorr = 0.2,
                             overdispersion = 9, maxR2changeratio = 0.05){
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
  object <- AddPCA(object, nPcsInit = nPcsInit, maxR2changeratio = maxR2changeratio)
  object <- AddAlleleFreqByTaxa(object)
  message("Performing preliminary genotype estimation.")
  object <- AddAlleleFreqHWE(object, excludeTaxa = excludeTaxa)
  object <- AddGenotypePriorProb_ByTaxa(object, selfing.rate)
  object <- AddGenotypeLikelihood(object, overdispersion = overdispersion)
  object <- AddPloidyChiSq(object, excludeTaxa = excludeTaxa)
  object <- AddGenotypePosteriorProb(object)
  object <- AddPCA(object, nPcsInit = dim(object$PCA)[2] + 1,
                   minPcsOut = dim(object$PCA)[2],
                   maxR2changeratio = maxR2changeratio)
  object <- AddAlleleFreqByTaxa(object, minfreq = minfreq)
  
  # Test for LD
  message("Finding alleles in LD.")
  object <- AddAlleleLinkages(object, type = "popstruct",
                              linkageDist = LDdist, minCorr = minLDcorr,
                              excludeTaxa = GetBlankTaxa(object))
  
  # Iteratively estimate genotypes
  while(meanDiff > tol){
    message(paste("Starting iteration", nIter))
    message(paste("PCs used:", dim(object$PCA)[2]))
    oldAlFreq <- object$alleleFreqByTaxa
    object <- AddAlleleFreqHWE(object, excludeTaxa = excludeTaxa)
    object <- AddGenotypePriorProb_ByTaxa(object, selfing.rate)
    object <- AddGenotypePriorProb_LD(object, type = "popstruct")
    object <- AddGenotypeLikelihood(object, overdispersion = overdispersion)
    object <- AddPloidyChiSq(object, excludeTaxa = excludeTaxa)
    object <- AddGenotypePosteriorProb(object)
    object <- AddPCA(object, nPcsInit = dim(object$PCA)[2] + 1,
                     minPcsOut = dim(object$PCA)[2],
                     maxR2changeratio = maxR2changeratio)
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

PipelineMapping2Parents <- function(object,
                                    n.gen.backcrossing = 0,
                                    n.gen.intermating = 0,
                                    n.gen.selfing = 0,
                                    minLikelihoodRatio = 10,
                                    freqAllowedDeviation = 0.05,
                                    freqExcludeTaxa = c(GetDonorParent(object),
                                                    GetRecurrentParent(object),
                                                    GetBlankTaxa(object)),
                                    useLinkage = TRUE, linkageDist = 1e7,
                                    minLinkageCorr = 0.5,
                                    overdispersion = 9){
  if(useLinkage && (is.null(object$locTable$Chr) ||
                    is.null(object$locTable$Pos))){
    stop("Set useLinkage = FALSE if alignment data unavailable.")
  }
  donorParent <- GetDonorParent(object)
  recurrentParent <- GetRecurrentParent(object)
  pld.don <- GetTaxaPloidy(object)[donorParent]
  pld.rec <- GetTaxaPloidy(object)[recurrentParent]
  # estimate possible allele frequencies
  message("Making initial parameter estimates...")
  pld.max <- max(sapply(object$possiblePloidies, sum))
  allelesin <- (pld.don + pld.rec) * pld.max / 2
  possfreq <- seq(0, 1, length.out = (n.gen.backcrossing + 1) * allelesin + 1)
  
  object <- AddAlleleFreqMapping(object, expectedFreqs = possfreq,
                                 allowedDeviation = freqAllowedDeviation,
                                 excludeTaxa = freqExcludeTaxa)
  # calculations for rest of pipeline
  object <- AddGenotypeLikelihood(object, overdispersion = overdispersion)
  object <- AddGenotypePriorProb_Mapping2Parents(object, donorParent = donorParent,
                                                 recurrentParent = recurrentParent,
                                                 n.gen.backcrossing = n.gen.backcrossing,
                                                 n.gen.intermating = n.gen.intermating,
                                                 n.gen.selfing = n.gen.selfing,
                                                 minLikelihoodRatio = minLikelihoodRatio)
  object <- AddPloidyChiSq(object, excludeTaxa = freqExcludeTaxa)
  object <- AddGenotypePosteriorProb(object)
  
  # add in linkage data if available
  if(useLinkage){
    message("Updating priors using linkage...")
    # find linkages
    object <- AddAlleleLinkages(object, type = "mapping", 
                                linkageDist = linkageDist, 
                                minCorr = minLinkageCorr,
                                excludeTaxa = freqExcludeTaxa)
    # update genotype probabilities
    object <- AddGenotypePriorProb_LD(object, type = "mapping")
    object <- AddGenotypePosteriorProb(object)
  } # end IF statement for using linked alleles
  message("Done.")
  
  return(object)
} # end PipelineMapping2Parents
