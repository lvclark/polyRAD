# class definition for storing data in polyRAD

# RADdata class constructor ####
RADdata <- function(alleleDepth, alleles2loc, locTable, possiblePloidies, 
                    contamRate, alleleNucleotides,
                    taxaPloidy = 2L){
  if(!is.integer(alleleDepth)){
    stop("alleleDepth must be in integer format.")
  }
  if(!is.matrix(alleleDepth)){
    stop("alleleDepth must be in matrix format.")
  }
  if(any(is.na(alleleDepth))){
    stop("There should be no missing data in alleleDepth; put 0 for zero depth.")
  }
  if(is.null(rownames(alleleDepth))){
    stop("alleleDepth must have taxa names as row names.")
  }
  if(length(alleles2loc) != dim(alleleDepth)[[2]]){
    stop("There must be one value of alleles2loc for each column of alleleDepth.")
  }
  if(!is.integer(alleles2loc)){
    stop("alleles2loc must be an integer.")
  }
  if(any(is.na(alleles2loc))){
    stop("No NA values allowed in alleles2loc.")
  }
  if(max(alleles2loc) > dim(locTable)[[1]]){
    stop("Each locus number in alleles2loc must correspond to a row in loctable.")
  }
  if(!is.data.frame(locTable)){
    stop("locTable must be a data frame.")
  }
  if(!is.null(locTable$Chr) && is.factor(locTable$Chr)){
    warning("Chromosomes should be coded as character or integer rather than factor in locTable.")
  }
  if(!is.list(possiblePloidies)){
    stop("possiblePloidies must be list")
  }
  possiblePloidies <- lapply(possiblePloidies, as.integer)
  if(!all(sapply(possiblePloidies, function(x) all(!is.na(x) & x > 0)))){
    stop("Each element of possiblePloidies should be a vector of integers greater than zero.")
  }
  if(contamRate < 0 || contamRate > 1){
    stop("contamRate can only range from zero to one.")
  }
  if(contamRate > 0.01){
    warning("contamRate higher than expected.")
  }
  if(!is.character(alleleNucleotides)){
    stop("alleleNucleotides must be a character vector.")
  }
  if(length(alleleNucleotides) != length(alleles2loc)){
    stop("Length of alleleNucleotides must be same as length of alleles2loc.")
  }
  taxaPloidy <- as.integer(taxaPloidy)
  if(length(taxaPloidy) == 1){
    # repeat ploidy for all samples if only one is provided
    taxaPloidy <- rep(taxaPloidy, times = nrow(alleleDepth))
  }
  if(length(taxaPloidy) != nrow(alleleDepth)){
    stop("Need one ploidy for each taxon.")
  }
  if(any(is.na(taxaPloidy))){
    stop("taxaPloidy must be integer.")
  }
  if(any(taxaPloidy < 1)){
    stop("taxaPloidy must be positive integer.")
  }
  if(is.null(names(taxaPloidy))){
    names(taxaPloidy) <- rownames(alleleDepth)
  } else {
    if(!setequal(names(taxaPloidy), rownames(alleleDepth))){
      stop("Sample names must match between taxaPloidy and alleleDepth.")
    }
    taxaPloidy <- taxaPloidy[rownames(alleleDepth)]
  }
  
  taxa <- dimnames(alleleDepth)[[1]]
  names(taxaPloidy) <- taxa
  nTaxa <- dim(alleleDepth)[1]
  nLoci <- dim(locTable)[1]
  locDepth <- t(apply(alleleDepth, 1, function(x) tapply(x, alleles2loc, sum)))
     # dimnames(locDepth)[[2]] is integer from alleles2loc converted to character
  if(length(unique(alleles2loc)) == 1){
    locDepth <- matrix(locDepth, nrow = nTaxa, ncol = nLoci,
                       dimnames = list(taxa, as.character(unique(alleles2loc))))
  }
  
  expandedLocDepth <- locDepth[,as.character(alleles2loc), drop = FALSE]
  
  # for each allele and taxon, get proportion of reads for that locus
  depthRatio <- alleleDepth/expandedLocDepth
  # depth of reads for each locus that do NOT belong to a given allele
  antiAlleleDepth <- expandedLocDepth - alleleDepth
  dimnames(antiAlleleDepth)[[2]] <- dimnames(alleleDepth)[[2]]
  
  return(structure(list(alleleDepth = alleleDepth, alleles2loc = alleles2loc,
                        locTable = locTable, possiblePloidies = possiblePloidies,
                        locDepth = locDepth,
                        depthRatio = depthRatio, antiAlleleDepth = antiAlleleDepth,
                        alleleNucleotides = alleleNucleotides,
                        taxaPloidy = taxaPloidy), 
                   class = "RADdata", taxa = taxa, nTaxa = nTaxa, nLoci = nLoci,
                   contamRate = contamRate))
}

# print method for RADdata (just some summary statistics) ####
print.RADdata <- function(x, ...){
  cat("## RADdata object ##", 
      paste(nTaxa(x), "taxa and", nLoci(x), "loci"),
      paste(sum(x$locDepth), "total reads"), 
      paste("Assumed sample cross-contamination rate of", 
            attr(x, "contamRate")),
      "\nPossible ploidies:", sep = "\n")
  printPloidies <- function(y){
    if(length(y) == 1){
      prefix <- "Auto"
    } else {
      prefix <- "Allo"
    }
    tot <- sum(y)
    if(tot == 1){
      suffix <- "ha"
    }
    if(tot == 2){
      suffix <- "di"
    }
    if(tot == 3){
      suffix <- "tri"
    }
    if(tot == 4){
      suffix <- "tetra"
    }
    if(tot == 5){
      suffix <- "penta"
    }
    if(tot == 6){
      suffix <- "hexa"
    }
    if(tot == 7){
      suffix <- "septa"
    }
    if(tot == 8){
      suffix <- "octo"
    }
    if(tot > 8){
      suffix <- "poly"
    }
    return(paste(prefix, suffix, "ploid (", paste(y, collapse = " "), ")", sep = ""))
  }
  txp <- unique(x$taxaPloidy)
  for(pl in x$possiblePloidies){
    for(tp in txp){
      cat(printPloidies(pl * tp / 2L), sep = "\n")
    }
  }
  if(!is.null(attr(x, "alleleFreqType"))){
    cat("", paste("Allele frequencies estimated for", attr(x, "alleleFreqType")),
        sep = "\n")
  }
}

#### parameter estimation generic functions and methods ####
AddDepthSamplingPermutations <- function(object, ...){
  UseMethod("AddDepthSamplingPermutations", object)
}
AddDepthSamplingPermutations.RADdata <- function(object, ...){
  # get log of number of permutations of order in which each allele could have
  # been sampled from total depth from that locus.
  expandedLocDepth <- object$alleleDepth + object$antiAlleleDepth
  depthSamplingPermutations <- lchoose(expandedLocDepth, object$alleleDepth)
  dimnames(depthSamplingPermutations)[[2]] <- dimnames(object$alleleDepth)[[2]]
  object$depthSamplingPermutations <- depthSamplingPermutations
  return(object)
}
# estimate allele frequencies assuming mapping population.
# Simple version using depth ratios.
AddAlleleFreqMapping <- function(object, ...){
  UseMethod("AddAlleleFreqMapping", object)
}
AddAlleleFreqMapping.RADdata <- function(object, 
                                         expectedFreqs = seq(0, 1, 0.25),
                                         allowedDeviation = 0.05,
                                         excludeTaxa = c(GetDonorParent(object),
                                                         GetRecurrentParent(object),
                                                         GetBlankTaxa(object)), ...){#,
#                                         deleteLociOutsideFreqRange = FALSE){
  maxdev <- min(dist(expectedFreqs, method = "manhattan"))/2
  if(maxdev < allowedDeviation && !isTRUE(all.equal(maxdev, allowedDeviation))){
    stop("allowedDeviation is too large given intervals within expectedFreqs")
  }
  if(!is.character(excludeTaxa)){
    stop("excludeTaxa must be a character vector (taxa names).")
  }
  for(ext in excludeTaxa[!excludeTaxa %in% GetTaxa(object)]){
    warning(paste(ext, "not found in taxa list."))
  }
  taxaToKeep <- GetTaxa(object)[!GetTaxa(object) %in% excludeTaxa]
  meanRat <- .alleleFreq(object, type = "choose", taxaToKeep = taxaToKeep)
  outFreq <- rep(NA, length(meanRat))
  names(outFreq) <- names(meanRat)
  attr(outFreq, "type") <- attr(meanRat, "type")
  for(f in expectedFreqs){
    outFreq[which(meanRat >= f - allowedDeviation & meanRat < f + allowedDeviation)] <- f
  }
  
  # fill in any missing allele frequencies for partially completed loci
  incompleteLoci <- unique(object$alleles2loc[is.na(outFreq)])
  for(L in incompleteLoci){
    theseal <- which(object$alleles2loc == L)
    if(sum(is.na(outFreq[theseal])) == 1){
      subtot <- sum(outFreq[theseal], na.rm = TRUE)
      outFreq[theseal[is.na(outFreq[theseal])]] <- 1 - subtot
    }
  }
  
  # for unexpected allele frequencies, keep them as they are
  outFreq[is.na(outFreq)] <- meanRat[is.na(outFreq)]
  
  # add allele frequencies to the object
  object$alleleFreq <- outFreq
  attr(object, "alleleFreqType") <- "mapping"
  
  # delete loci that don't match expected allele frequencies
#  if(deleteLociOutsideFreqRange){
#    lociToDelete <- unique(object$allele2loc[is.na(outFreq)])
#    object <- SubsetLoci(object, loci = lociToDelete, delete = TRUE)
#  }
  
  return(object)
}

# estimate allele frequencies assuming a diversity panel in Hardy-Weinberg
# Equilibrium.  Simple version using depth ratios.
AddAlleleFreqHWE <- function(object, ...){
  UseMethod("AddAlleleFreqHWE", object)
}
AddAlleleFreqHWE.RADdata <- function(object, excludeTaxa = GetBlankTaxa(object),
                                     ...){
  taxaToKeep <- GetTaxa(object)[!GetTaxa(object) %in% excludeTaxa]
  object$alleleFreq <- .alleleFreq(object, type = "choose", 
                                   taxaToKeep = taxaToKeep)
  attr(object, "alleleFreqType") <- "HWE"
  return(object)
}

# estimate likelihood of genotypes given the read depth
# (i.e. the probability of that read count distribution given each 
# possible genotype)
AddGenotypeLikelihood <- function(object, ...){
  UseMethod("AddGenotypeLikelihood", object)
}
AddGenotypeLikelihood.RADdata <- function(object, overdispersion = 9, ...){
  if(is.null(object$alleleFreq)){
    message("Allele frequencies not found; estimating under HWE from depth ratios.")
    object <- AddAlleleFreqHWE(object)
  }
  alFreq <- object$alleleFreq
  if(is.null(object$depthSamplingPermutations)){
    message("Generating sampling permutations for allele depth.")
    object <- AddDepthSamplingPermutations(object)
  }
  
  # get ploidies, ignoring inheritance pattern
  ploidies <- sort(unique(sapply(object$possiblePloidies, sum)))
  # fix any allele freq that are zero, to prevent NaN likelihood
  minfreq <- 1/nTaxa(object)/max(ploidies)
  alFreq[alFreq == 0] <- minfreq
  alFreq[alFreq == 1] <- 1 - minfreq
  # get ploidies by taxa
  tx_pld_unique <- sort(unique(GetTaxaPloidy(object)))
  
  # set up list for genotype likelihoods and loop through
  object$genotypeLikelihood <- array(list(),
                                     dim = c(length(ploidies),
                                             length(tx_pld_unique)),
                                     dimnames = list(NULL, as.character(tx_pld_unique)))
  # probability of getting each allele from contamination
  sampleContam <- attr(object, "contamRate") * alFreq
  for(i in seq_along(ploidies)){
    for(h in seq_along(tx_pld_unique)){
      pldtot <- sum(object$possiblePloidies[[i]]) * tx_pld_unique[h] / 2L
      # get probability of sampling each allele from each possible genotype
      sampleReal <- (0:pldtot)/pldtot * (1 - attr(object, "contamRate"))
      alleleProb <- matrix(0, nrow = length(sampleReal), 
                           ncol = length(sampleContam))
      for(j in seq_along(sampleReal)){
        alleleProb[j,] <- sampleReal[j] + sampleContam
      }
      antiAlleleProb <- 1 - alleleProb
      # multiply probabilities by overdispersion factor for betabinomial
      alleleProb <- alleleProb * overdispersion
      antiAlleleProb <- antiAlleleProb * overdispersion
      
      thesetaxa <- GetTaxaByPloidy(object, tx_pld_unique[h])
      
      # get likelihoods
      object$genotypeLikelihood[[i,h]] <- array(0, dim = c(pldtot+1, 
                                                         length(thesetaxa), nAlleles(object)),
                                              dimnames = list(as.character(0:pldtot),
                                                              thesetaxa,
                                                              GetAlleleNames(object)))
      for(j in 1:(pldtot+1)){
        # likelihoods under beta-binomial distribution
        object$genotypeLikelihood[[i,h]][j,,] <- 
          exp(object$depthSamplingPermutations[thesetaxa,] +
                sweep(lbeta(sweep(object$alleleDepth[thesetaxa,], 2, alleleProb[j,], "+"),
                            sweep(object$antiAlleleDepth[thesetaxa,], 2, antiAlleleProb[j,], "+")),
                      2, lbeta(alleleProb[j,], antiAlleleProb[j,]), "-"))
      }
      # fix likelihoods where all are zero
      totlik <- colSums(object$genotypeLikelihood[[i,h]])
      toRecalculate <- which(totlik == 0, arr.ind = TRUE)
      if(dim(toRecalculate)[1] > 0){
        for(k in 1:dim(toRecalculate)[1]){
          taxon <- toRecalculate[k,1]
          allele <- toRecalculate[k,2]
          # for rare cases where likelihood still not estimated, set to one
          if(sum(object$genotypeLikelihood[[i,h]][, taxon, allele]) == 0){
            object$genotypeLikelihood[[i,h]][, taxon, allele] <- 1
          }
        }
      }
    }
  }
  
  return(object)
}

# for a given taxon, return the most likely genotypes (based on likelihoods only)
GetLikelyGen <- function(object, taxon, minLikelihoodRatio = 10){
  UseMethod("GetLikelyGen", object)
}
GetLikelyGen.RADdata <- function(object, taxon, minLikelihoodRatio = 10){
  if(length(taxon) != 1){
    stop("Only one taxon can be passed to GetLikelyGen.")
  }
  taxind <- match(taxon, GetTaxa(object))
  if(is.na(taxind)){
    stop("taxon not found in object.")
  }
  if(length(minLikelihoodRatio) != 1 || is.na(minLikelihoodRatio) ||
     !is.numeric(minLikelihoodRatio)){
    stop("A single numeric value is needed for minLikelihoodRatio.")
  }
  if(minLikelihoodRatio < 1){
    warning("minimumLikelihoodRatio less than 1 does not make sense.")
  }
  if(is.null(object$genotypeLikelihood)){
    cat("Genotype likelihoods not found.  Estimating.", sep = "\n")
    object <- AddGenotypeLikelihood(object)
  }
  npld <- length(object$genotypeLikelihood)
  ploidies <- sapply(object$genotypeLikelihood, function(x) dim(x)[1] - 1)
  nAllele <- nAlleles(object)
  
  outmat <- matrix(NA_integer_, nrow = npld, ncol = nAllele,
                   dimnames = list(as.character(ploidies), 
                                   GetAlleleNames(object)))
  for(i in 1:npld){
    nonNaAlleles <- which(!is.na(object$genotypeLikelihood[[i]][1,taxind,]))
    # get the most likely genotype
#    outmat[i,nonNaAlleles] <- apply(object$genotypeLikelihood[[i]][,taxind,nonNaAlleles], 2, which.max) - 1 # R version
    outmat[i,nonNaAlleles] <- BestGenos(object$genotypeLikelihood[[i]][,taxind,nonNaAlleles],
                                        ploidies[i], 1, length(nonNaAlleles)) # Rcpp function
    # remove genotypes that don't meet the likelihood ratio threshold
    if(minLikelihoodRatio > 1){
      # get likelihood ratios
      myrat <- apply(object$genotypeLikelihood[[i]][,taxind,nonNaAlleles], 2, 
                     function(x){
                       xmxind <- which.max(x)
                       xsecond <- max(x[-xmxind])
                       return(x[xmxind]/xsecond)
                     })
      outmat[i,nonNaAlleles[myrat < minLikelihoodRatio]] <- NA_integer_
    }
  }
  return(outmat)
}

# Estimate parental genotypes.  Used for getting genotype priors in mapping
# population, and for Hind/He statistic.
EstimateParentalGenotypes <- function(object, ...){
  UseMethod("EstimateParentalGenotypes", object)
}
EstimateParentalGenotypes.RADdata <- 
  function(object,
           donorParent = GetDonorParent(object), 
           recurrentParent = GetRecurrentParent(object), n.gen.backcrossing = 0,
           n.gen.intermating = 0,
           n.gen.selfing = 0, donorParentPloidies = object$possiblePloidies,
           recurrentParentPloidies = object$possiblePloidies,
           minLikelihoodRatio = 10, ...){
    if(any(!donorParentPloidies %in% object$possiblePloidies) ||
       any(!recurrentParentPloidies %in% object$possiblePloidies)){
      # make sure we have all parental ploidies so we can get likelihoods 
      ### change this? ploidy by individual in future?
      stop("All parent ploidies must be in the possible ploidies for the object")
    }
    if(is.null(object$alleleFreq) || attr(object,"alleleFreqType") != "mapping"){
      message("Allele frequencies for mapping population not found.  Estimating.")
      allelesin <- max(sapply(donorParentPloidies, sum)) + 
        max(sapply(recurrentParentPloidies, sum))
      possfreq <- seq(0, 1, length.out = (n.gen.backcrossing + 1) * allelesin + 1)
      alldev <- (possfreq[2] - possfreq[1])/2
      object <- AddAlleleFreqMapping(object, allowedDeviation = alldev, 
                                     expectedFreqs = possfreq)
    }
    if(is.null(object$genotypeLikelihood)){
      message("Genotype likelihoods not found.  Estimating.")
      object <- AddGenotypeLikelihood(object)
    }
    
    pldtot <- sapply(object$genotypeLikelihood, function(x) dim(x)[1] - 1)
    pldtot.don <- pldtot[pldtot %in% sapply(donorParentPloidies, sum)]
    pldtot.rec <- pldtot[pldtot %in% sapply(recurrentParentPloidies, sum)]
    likelyGen.don <- GetLikelyGen(object, donorParent, 
                                  minLikelihoodRatio = minLikelihoodRatio)[as.character(pldtot.don),,drop = FALSE]
    likelyGen.rec <- GetLikelyGen(object, recurrentParent,
                                  minLikelihoodRatio = minLikelihoodRatio)[as.character(pldtot.rec),,drop = FALSE]
    # combinations of parent ploidies that match list of progeny ploidies
    pldcombos <- matrix(NA, nrow = 0, ncol = 2, 
                        dimnames = list(NULL,c("donor","recurrent")))
    # matrix of expected allele frequencies for all possible genotypes and ploidies
    expfreq_byPloidy <- list()
    # find possible combinations of parent ploidies, and possible expected allele frequencies
    for(pl.d in pldtot.don){
      for(pl.r in pldtot.rec){
        if((pl.d/2 + pl.r/2) %in% pldtot && (n.gen.backcrossing == 0 || pl.d == pl.r)){
          pldcombos <- rbind(pldcombos, matrix(c(pl.d, pl.r), nrow = 1, ncol = 2))
          
          expfreq_byPloidy[[length(expfreq_byPloidy) + 1]] <- matrix(nrow = pl.d+1, ncol = pl.r+1)
          for(gen.d in 0:pl.d){
            for(gen.r in 0:pl.r){
              expfreq_byPloidy[[length(expfreq_byPloidy)]][gen.d + 1, gen.r + 1] <- 
                (gen.d * 0.5^n.gen.backcrossing + gen.r * (2 - 0.5^n.gen.backcrossing))/
                (pl.d + pl.r)
            }
          }
        }
      }
    }
    
    # do allele frequencies match parent genotypes?
    freqMatchGen <- matrix(FALSE, nrow = dim(pldcombos)[1], ncol = nAlleles(object))
    for(i in 1:dim(pldcombos)[1]){
      thisgen.don <- likelyGen.don[as.character(pldcombos[i,"donor"]),]
      thisgen.rec <- likelyGen.rec[as.character(pldcombos[i,"recurrent"]),]
      expfreq <- (thisgen.don * 0.5^n.gen.backcrossing + 
                    thisgen.rec * (2 - 0.5^n.gen.backcrossing))/(pldcombos[i,"donor"] + 
                                                                   pldcombos[i,"recurrent"])
      freqMatchGen[i,] <- expfreq == object$alleleFreq
    }
    freqMatchGen[is.na(freqMatchGen)] <- FALSE
    allelesToFix <- which(colSums(freqMatchGen) == 0)
    # correct parental genotypes where appropriate, using rounded allele frequencies
    for(a in allelesToFix){
      thisfreq <- object$alleleFreq[a]
      for(i in 1:dim(pldcombos)[1]){
        poss_matches <- which(expfreq_byPloidy[[i]] == thisfreq, arr.ind = TRUE) - 1
        if(nrow(poss_matches) == 0) next
        if(nrow(poss_matches) == 1){
          # only one possible match (i.e. when there is backcrossing)
          likelyGen.don[as.character(pldcombos[i,"donor"]),a] <- 
            unname(poss_matches[,1])
          likelyGen.rec[as.character(pldcombos[i,"recurrent"]),a] <- 
            unname(poss_matches[,2])
        } else { # multiple possible matches
          # vector to contain genotype combo likelihoods
          thislikeli <- numeric(nrow(poss_matches)) 
          # indices for current ploidies
          plind.d <- which(pldtot == pldcombos[i,"donor"])
          plind.r <- which(pldtot == pldcombos[i,"recurrent"])
          for(m in 1:nrow(poss_matches)){
            gen.d <- poss_matches[m,1]
            gen.r <- poss_matches[m,2]
            
            thislikeli[m] <- 
              object$genotypeLikelihood[[plind.d]][gen.d + 1, donorParent, a] *
              object$genotypeLikelihood[[plind.r]][gen.r + 1, recurrentParent, a]
          }
          bestcombo <- which(thislikeli == max(thislikeli))
          if(length(bestcombo) == 1){
            likelyGen.don[as.character(pldcombos[i,"donor"]),a] <- 
              unname(poss_matches[bestcombo, 1])
            likelyGen.rec[as.character(pldcombos[i,"recurrent"]),a] <- 
              unname(poss_matches[bestcombo, 2])
          }
        }
      }
    }
    
    # save corrected parental genotypes to object
    object$likelyGeno_donor <- likelyGen.don
    object$likelyGeno_recurrent <- likelyGen.rec
    object$pldcombos <- pldcombos # this is needed to add genotype prior prob
    return(object)
}
# for a mapping population with two parents, get prior genotype probabilities
# based on parent genotypes and progeny allele frequencies
AddGenotypePriorProb_Mapping2Parents <- function(object, ...){
  UseMethod("AddGenotypePriorProb_Mapping2Parents", object)
}
AddGenotypePriorProb_Mapping2Parents.RADdata <- function(object,
    donorParent = GetDonorParent(object), 
    recurrentParent = GetRecurrentParent(object), n.gen.backcrossing = 0,
    n.gen.intermating = 0,
    n.gen.selfing = 0, donorParentPloidies = object$possiblePloidies,
    recurrentParentPloidies = object$possiblePloidies,
    minLikelihoodRatio = 10, ...){
  # get most likely genotype for the parents
  object <- EstimateParentalGenotypes(object, donorParent = donorParent, 
                                      recurrentParent = recurrentParent, n.gen.backcrossing = n.gen.backcrossing,
                                      n.gen.intermating = n.gen.intermating,
                                      n.gen.selfing = n.gen.selfing, donorParentPloidies = donorParentPloidies,
                                      recurrentParentPloidies = recurrentParentPloidies,
                                      minLikelihoodRatio = minLikelihoodRatio)
  likelyGen.don <- object$likelyGeno_donor
  likelyGen.rec <- object$likelyGeno_recurrent
  pldcombos <- object$pldcombos
  
  # function for deteriming ploidy of offspring (by index in possiblePloidies)
  offspringPloidy <- function(pld1, pld2){ 
    if(length(pld1) == length(pld2)){
      newPld <- (pld1 + pld2)/2
    } else {
      newPld <- (sum(pld1) + sum(pld2))/2
    }
    newPld <- as.integer(newPld)
    out <- which(sapply(object$possiblePloidies, function(x) identical(x, newPld)))
    if(length(out) == 0){
      out <- which(sapply(object$possiblePloidies, sum) == sum(newPld))
    }
    return(out)
  }
  
  # expand ploidy combinations across allo and auto types, using index in possiblePloidies
  pldtot2 <- sapply(object$possiblePloidies, sum)
  if(n.gen.backcrossing == 0){
    bcnames <- character(0)
    bcloop <- integer(0)
  } else {
    bcloop <- 1:n.gen.backcrossing
    bcnames <- paste("BC", bcloop, sep = "")
  }
  pldcombosExpand <- array(0L, dim = c(0, 4 + n.gen.backcrossing), 
                           dimnames = list(NULL, c("donor","recurrent","F1",
                                                   bcnames, "final")))
  for(i in 1:dim(pldcombos)[1]){
    thesepl.don <- which(pldtot2 == pldcombos[i,"donor"])
    thesepl.rec <- which(pldtot2 == pldcombos[i,"recurrent"])
    for(pl.d in thesepl.don){
      for(pl.r in thesepl.rec){
        pldcombosExpand <- rbind(pldcombosExpand, 
                                 array(c(pl.d, pl.r, rep(NA, 2+n.gen.backcrossing)),
                                       dim = c(1,4 + n.gen.backcrossing)))
        newrow <- dim(pldcombosExpand)[1]
        pldcombosExpand[newrow, "F1"] <- offspringPloidy(object$possiblePloidies[[pl.d]],
                                                         object$possiblePloidies[[pl.r]])
        thiscol <- 3
        for(b in bcloop){
          thiscol <- match(paste("BC", b, sep = ""), dimnames(pldcombosExpand)[[2]])
          pldcombosExpand[newrow, thiscol] <-
            offspringPloidy(object$possiblePloidies[[pl.r]], 
                            object$possiblePloidies[[pldcombosExpand[newrow, thiscol - 1]]])
        }
        pldcombosExpand[newrow, "final"] <- pldcombosExpand[newrow, thiscol]
      }
    }
  }

  # get prior genotype probabilities for F1
  OutPriors <- list()
  length(OutPriors) <- dim(pldcombosExpand)[1]
  for(i in 1:length(OutPriors)){
    donorPld <- object$possiblePloidies[[pldcombosExpand[i,"donor"]]]
    recurPld <- object$possiblePloidies[[pldcombosExpand[i,"recurrent"]]]
    theseDonorGen <- likelyGen.don[as.character(sum(donorPld)),]
    theseRecurGen <- likelyGen.rec[as.character(sum(recurPld)),]
    donorGamProb <- .gameteProb(.makeGametes(theseDonorGen, donorPld), donorPld)
    recurGamProb <- .gameteProb(.makeGametes(theseRecurGen, recurPld), recurPld)
    OutPriors[[i]] <- .progenyProb(donorGamProb, recurGamProb)
  }
  # backcross
  for(gen in bcloop){
    # reestimate prior probs
    OutPriorsLastGen <- OutPriors
    OutPriors <- list()
    length(OutPriors) <- dim(pldcombosExpand)[1]
    for(i in 1:length(OutPriors)){
      ### consider just estimating recurrent gametes once since this is repetitive
      recurPld <- object$possiblePloidies[[pldcombosExpand[i,"recurrent"]]]
      theseRecurGen <- likelyGen.rec[as.character(sum(recurPld)),]
      recurGamProb <- .gameteProb(.makeGametes(theseRecurGen, recurPld), recurPld)
      # get gamete prob for current population
      currPld <- object$possiblePloidies[[pldcombosExpand[i,paste("BC", gen, sep = "")]]]
      currGamProb <- .gameteProbPop(OutPriorsLastGen[[i]], currPld)
      # update genotype priors
      OutPriors[[i]] <- .progenyProb(recurGamProb, currGamProb)
    }
  }
  # intermate (random mating within the population)
  if(n.gen.intermating == 0){
    mateloop <- integer(0)
  } else {
    mateloop <- 1:n.gen.intermating
  }
  for(gen in mateloop){
    OutPriorsLastGen <- OutPriors
    OutPriors <- list()
    length(OutPriors) <- dim(pldcombosExpand)[1]
    for(i in 1:length(OutPriors)){
      currPld <- object$possiblePloidies[[pldcombosExpand[i,"final"]]]
      currGamProb <- .gameteProbPop(OutPriorsLastGen[[i]], currPld)
      OutPriors[[i]] <- .progenyProb(currGamProb, currGamProb)
    }
  }
  # self (everything in population is self-fertilized)
  if(n.gen.selfing == 0){
    selfloop <- integer(0)
  } else {
    selfloop <- 1:n.gen.selfing
  }
  for(gen in selfloop){
    # reestimate prior probs
    OutPriorsLastGen <- OutPriors
    OutPriors <- list()
    length(OutPriors) <- dim(pldcombosExpand)[1]
    for(i in 1:length(OutPriors)){
      currPld <- object$possiblePloidies[[pldcombosExpand[i,"final"]]]
      OutPriors[[i]] <- .selfPop(OutPriorsLastGen[[i]], currPld)
    }
  }
  
  for(i in 1:length(OutPriors)){
    dimnames(OutPriors[[i]]) <- list(as.character(0:(dim(OutPriors[[i]])[1] - 1)),
                                     GetAlleleNames(object))
  }

  object$priorProb <- OutPriors
  object$priorProbPloidies <- object$possiblePloidies[pldcombosExpand[,"final"]]
  object$donorPloidies <- object$possiblePloidies[pldcombosExpand[,"donor"]]
  object$recurrentPloidies <- object$possiblePloidies[pldcombosExpand[,"recurrent"]]
  attr(object, "priorType") <- "population" 
  # --> indicates prior probs are estimated for whole pop, not by taxa
  return(object)
}

AddGenotypePriorProb_HWE <- function(object, ...){
  UseMethod("AddGenotypePriorProb_HWE", object)
}
AddGenotypePriorProb_HWE.RADdata <- function(object, selfing.rate = 0, ...){
  if(is.null(object$alleleFreq) || attr(object, "alleleFreqType") != "HWE"){
    stop("Allele frequencies not estimated under HWE.")
  }

  tx_pld_unique <- sort(unique(GetTaxaPloidy(object)))
  # array of priors by marker, with marker inheritance patterns in rows and
  # individual ploidies in columns.
  priors <- array(list(),
                  dim = c(length(object$possiblePloidies),
                          length(tx_pld_unique)),
                  dimnames = list(NULL, as.character(tx_pld_unique)))
  
  for(i in seq_along(object$possiblePloidies)){
    for(j in seq_along(tx_pld_unique)){
      priors[[i, j]] <- .HWEpriors(object$alleleFreq,
                                   object$possiblePloidies[[i]] * tx_pld_unique[j] / 2L,
                                   selfing.rate)
    }
  }
  
  object$priorProb <- priors
  object$priorProbPloidies <- object$possiblePloidies
  attr(object, "priorType") <- "population"
  return(object)
}

AddGenotypePriorProb_Even <- function(object, ...){
  UseMethod("AddGenotypePriorProb_Even", object)
}
AddGenotypePriorProb_Even.RADdata <- function(object, ...){
  tx_pld_unique <- sort(unique(GetTaxaPloidy(object)))
  priors <- array(list(),
                  dim = c(length(object$possiblePloidies),
                          length(tx_pld_unique)),
                  dimnames = list(NULL, as.character(tx_pld_unique)))
  
  for(i in seq_along(object$possiblePloidies)){
    for(j in seq_along(tx_pld_unique)){
      thispld <- sum(object$possiblePloidies[[i]] * tx_pld_unique[j] / 2L)
      priors[[i, j]] <- matrix(1/(thispld + 1), nrow = thispld + 1, 
                            ncol = nAlleles(object),
                            dimnames = list(as.character(0:thispld),
                                            GetAlleleNames(object)))
    }
  }
  object$priorProb <- priors
  object$priorProbPloidies <- object$possiblePloidies
  attr(object, "priorType") <- "population"
  return(object)
}

AddPloidyLikelihood <- function(object, ...){
  UseMethod("AddPloidyLikelihood", object)
}
AddPloidyLikelihood.RADdata <- function(object, excludeTaxa = GetBlankTaxa(object), 
                                        minLikelihoodRatio = 50, ...){
  if(attr(object, "priorType") != "population"){
    stop("AddPloidyLikelihood not yet defined for priors estimated on a per-taxon basis.")
  }
  if(ncol(object$priorProb) > 1){
    stop("AddPloidyLikelihood not yet defined for multiploid populations.")
  }
  taxa <- GetTaxa(object)
  if(!is.null(attr(object, "donorParent"))){
    taxa <- taxa[taxa != GetDonorParent(object)]
  }
  if(!is.null(attr(object, "recurrentParent"))){
    taxa <- taxa[taxa != GetRecurrentParent(object)]
  }
  taxa <- taxa[!taxa %in% GetBlankTaxa(object)]
  taxa <- taxa[!taxa %in% excludeTaxa]
  
  nAlleles <- nAlleles(object)
  
  likgen <- lapply(taxa, function(x) GetLikelyGen(object, x,
                                                  minLikelihoodRatio = minLikelihoodRatio))
  object$ploidyLikelihood <- matrix(nrow = length(object$priorProb),
                                    ncol = nAlleles,
                                    dimnames = list(NULL, GetAlleleNames(object)))
  for(i in 1:nrow(object$priorProb)){
    thisploidy <- dim(object$priorProb[[i,1]])[1] - 1
    thesegen <- sapply(likgen, function(x) x[as.character(thisploidy),])
    countstable <- matrix(0L, nrow = thisploidy + 1, ncol = dim(thesegen)[1],
                          dimnames = list(as.character(0:thisploidy),
                                          dimnames(thesegen)[[1]]))
    for(j in 0:thisploidy){
      countstable[j + 1, ] <- rowSums(thesegen == j, na.rm = TRUE)
    }
    # get ploidy likelihood
    thislikehd <- sapply(1:nAlleles, function(x){
      if(any(is.na(object$priorProb[[i,1]][,x]))){
        NA
      } else {
        dmultinom(countstable[,x], 
                  prob = object$priorProb[[i,1]][,x])
      }
    })
    object$ploidyLikelihood[i,] <- thislikehd
    # then do chi-squared test?
  }
  return(object)
}

AddPloidyChiSq <- function(object, ...){
  UseMethod("AddPloidyChiSq", object)
}
AddPloidyChiSq.RADdata <- function(object, excludeTaxa = GetBlankTaxa(object),
                                   ...){
  taxa <- GetTaxa(object)
  if(!is.null(attr(object, "donorParent"))){
    taxa <- taxa[taxa != GetDonorParent(object)]
  }
  if(!is.null(attr(object, "recurrentParent"))){
    taxa <- taxa[taxa != GetRecurrentParent(object)]
  }
  taxa <- taxa[!taxa %in% GetBlankTaxa(object)]
  taxa <- taxa[!taxa %in% excludeTaxa]
  
  if(is.null(object$priorProb)){
    stop("Prior genotype probabilities must be estimated first.")
  }
  if(is.null(object$genotypeLikelihood)){
    object <- AddGenotypeLikelihood(object)
  }
  
  nAllele <- nAlleles(object)
  object$ploidyChiSq <- matrix(0, nrow = nrow(object$priorProb),
                               ncol = nAllele,
                               dimnames = list(NULL, GetAlleleNames(object)))
  # object$ploidyChiSqP <- matrix(NA, nrow = length(object$priorProb),
  #                               ncol = nAllele,
  #                               dimnames = list(NULL, GetAlleleNames(object)))
  
  # get weighted genotype tallies from genotype likelihoods
  gental <- array(list(),
                  dim = dim(object$genotypeLikelihood),
                  dimnames = dimnames(object$genotypeLikelihood))
  for(i in seq_len(nrow(gental))){
    for(h in seq_len(ncol(gental))){
      thesetaxa <- intersect(taxa, dimnames(object$genotypeLikelihood[[i,h]])[[2]])
      # likelihood total for each individual and locus at this ploidy
      totlik <- colSums(object$genotypeLikelihood[[i,h]][,thesetaxa,])
      # normalize likelihoods by total for each individual and locus
      normlik <- sweep(object$genotypeLikelihood[[i,h]][,thesetaxa,],
                       2:3, totlik, FUN = "/")
      # get population proportion of genotype likelihoods for allele and copy number
      gental[[i,h]] <- apply(normlik, c(1,3), mean)
    }
  }
  
  # loop through ploidies
  for(i in seq_len(nrow(object$priorProb))){
    thisploidy <- sum(object$priorProbPloidies[[i]])
    whichlik <-
      which(sapply(object$genotypeLikelihood[,1], 
                   function(x){
                     dim(x)[1] - 1L == thisploidy *
                       as.integer(colnames(object$genotypeLikelihood)[1]) / 2L
                   }))
    stopifnot(length(whichlik) == 1L)
    stopifnot(all(sapply(colnames(object$genotypeLikelihood),
                         function(x){
                           dim(object$genotypeLikelihood[[whichlik,x]])[1] - 1L ==
                             thisploidy * as.integer(x) / 2L
                         })))
    for(h in seq_len(ncol(object$priorProb))){
      # get priors
      if(attr(object, "priorType") == "population"){
        thesepriors <- object$priorProb[[i,h]]
      } else {
        # convert priors by taxon to population priors
        thesepriors <- rowMeans(aperm(object$priorProb[[i,h]], c(1,3,2)), dims = 2)
      }
      # estimate the components that are summed to make chi square
      chisqcomp <- (gental[[whichlik,h]] - thesepriors)^2/
        thesepriors * length(taxa)
      # chi-squared statistic
      thesechisq <- apply(chisqcomp, 2, function(x) sum(x[x != Inf]))
      object$ploidyChiSq[i,] <- object$ploidyChiSq[i,] + thesechisq
    }
    # degrees of freedom
    # theseDF <- colSums(thesepriors != 0) - 1
    # p-values
    # object$ploidyChiSqP[i,] <- pchisq(thesechisq, theseDF, lower.tail = FALSE)
  }
  
  return(object)
}

AddPriorTimesLikelihood <- function(object, ...){
  UseMethod("AddPriorTimesLikelihood", object)
}
AddPriorTimesLikelihood.RADdata <- function(object, ...){
  results <- .priorTimesLikelihood(object)
  
  object$priorTimesLikelihood <- results
  return(object)
}

AddGenotypePosteriorProb <- function(object, ...){
  UseMethod("AddGenotypePosteriorProb", object)
}
AddGenotypePosteriorProb.RADdata <- function(object, ...){
  if(is.null(object$priorTimesLikelihood)){
    PTL <- .priorTimesLikelihood(object)
  } else {
    PTL <- object$priorTimesLikelihood
  }
  object$posteriorProb <- array(list(),
                                dim = dim(PTL),
                                dimnames = dimnames(PTL))

  for(i in seq_len(nrow(object$posteriorProb))){
    for(h in seq_len(ncol(object$posteriorProb))){
      totPriorTimesLikeli <- colSums(PTL[[i,h]])
      object$posteriorProb[[i,h]] <- sweep(PTL[[i,h]], c(2,3),
                                         totPriorTimesLikeli, FUN = "/")
    }
  }
  return(object)
}


GetWeightedMeanGenotypes <- function(object, ...){
  UseMethod("GetWeightedMeanGenotypes", object)
}
GetWeightedMeanGenotypes.RADdata <- function(object, minval = 0, maxval = 1,
                                             omit1allelePerLocus = TRUE, 
                                             omitCommonAllele = TRUE,
                                             naIfZeroReads = FALSE, 
                                             onePloidyPerAllele = FALSE, ...){
  # maybe include an argument for selecting a specific ploidy rather than
  # letting the function pick what seems to be best?
  if(is.null(object$posteriorProb)){
    stop("Need to estimate genotype posterior probabilities first.")
  }
  if(!CanDoGetWeightedMeanGeno(object)){
    stop("Need to estimate ploidy chi-squared first.")
  }

  altokeep <- 1:nAlleles(object)
  if(omit1allelePerLocus){
    # make allele subset, to remove mathematical redundancy in dataset
    mymatch <- OneAllelePerMarker(object, commonAllele = omitCommonAllele)
    altokeep <- altokeep[-mymatch]
  }  
  
  nPloidies <- nrow(object$priorProb)
  
  # get weights for ploidies to use for each allele
  if(is.null(object$ploidyChiSq)){
    ploidyweights <- matrix(1, nrow = 1, ncol = length(altokeep))
  } else {
    nPloidies <- dim(object$ploidyChiSq)[1]
    if(onePloidyPerAllele){
      ploidyweights <- matrix(0, nrow = nPloidies,
                              ncol = length(altokeep))
      bestploidies <- BestPloidies(object$ploidyChiSq[,altokeep, drop = FALSE])
      for(i in 1:nPloidies){
        ploidyweights[i,bestploidies == i] <- 1
      }
    } else {
      chisqInverse <- 1/object$ploidyChiSq[,altokeep, drop = FALSE]
      chisqInverse[is.na(chisqInverse)] <- 0
      ploidyweights <- sweep(chisqInverse, 2, colSums(chisqInverse), "/")
    }
    ploidyweights[is.na(ploidyweights)] <- 0
  }

  # set up emtpy matrix to contain results
  wmgeno <- matrix(NA, nrow = nTaxa(object),
                   ncol = length(altokeep),
                   dimnames = list(GetTaxa(object),
                                   GetAlleleNames(object)[altokeep]))
  # loop through locus ploidies
  for(i in 1:nPloidies){
    # loop through taxa ploidies
    for(h in seq_len(ncol(object$priorProb))){
      # taxa with this ploidy
      thesetaxa <- dimnames(object$posteriorProb[[i,h]])[[2]]
      # values to represent each allele copy number
      thesegenval <- seq(minval, maxval, length.out = dim(object$posteriorProb[[i,h]])[1])
      # weighted mean genotypes for this ploidy
      thesewm <- rowSums(aperm(sweep(object$posteriorProb[[i,h]][,,altokeep],
                                     1, thesegenval, "*"), 
                               c(2,3,1)), 
                         dims = 2)
      # multiply by weight for this ploidy and add to total
      thesewm <- sweep(thesewm, 2, ploidyweights[i,], "*")
      nasofar <- is.na(wmgeno[thesetaxa,])
      nanew <- is.na(thesewm)
      add <- !nasofar & !nanew
      wmgeno[thesetaxa,][nasofar] <- thesewm[nasofar]
      wmgeno[thesetaxa,][add] <- wmgeno[thesetaxa,][add] + thesewm[add]
    }
  }
  # where there was no good ploidy, put in missing data
  wmgeno[,colSums(ploidyweights) == 0] <- NA
  
  if(naIfZeroReads){ # if there were zero reads for that locus, replace with NA
    wmgeno[object$locDepth[,as.character(object$alleles2loc)[altokeep]] == 0] <- NA
  }
  
  return(wmgeno)
}

# function to get the most probable genotypes at the most probable ploidy
GetProbableGenotypes <- function(object, ...){
  UseMethod("GetProbableGenotypes", object)
}
GetProbableGenotypes.RADdata <- function(object, omit1allelePerLocus = TRUE,
                                         omitCommonAllele = TRUE,
                                         naIfZeroReads = FALSE, 
                                         correctParentalGenos = TRUE, 
                                         multiallelic = "correct", ...){
  if(!CanDoGetWeightedMeanGeno(object)){
    stop("Need posteriorProb and ploidyChiSq.")
  }
  if(!multiallelic %in% c("ignore", "na", "correct")){
    stop("multiallelic should equal 'ignore', 'na', or 'correct'.")
  }
  # determine which alleles should be processed
  allelesToExport <- 1:nAlleles(object)
  if(omit1allelePerLocus && multiallelic == "ignore"){
    allelesToExport <- allelesToExport[-OneAllelePerMarker(object,
                                                           commonAllele = omitCommonAllele)]
  }
  
  # matrix for output
  outmat <- matrix(NA_integer_, nrow = nTaxa(object),
                   ncol = length(allelesToExport),
                   dimnames = list(GetTaxa(object), 
                                   GetAlleleNames(object)[allelesToExport]))
  # index of which ploidy was exported for each allele
  if(nrow(object$posteriorProb) == 1){
    pldindex <- rep(1L, length(allelesToExport))
    ploidyindices <- 1L
  } else {
    if(multiallelic == "ignore"){
      pldindex <- BestPloidies(object$ploidyChiSq[,allelesToExport,drop=FALSE]) # Rcpp function
    } else {
      # If we are going to look at multiallelic genotypes, get the best ploidy
      # for each locus instead of each allele.
      locChiSq <- t(rowsum(t(object$ploidyChiSq), object$alleles2loc))
      pldindex_loc <- BestPloidies(locChiSq)
      names(pldindex_loc) <- colnames(locChiSq)
      pldindex <- unname(pldindex_loc[as.character(object$alleles2loc)])
    }
    ploidyindices <- sort(unique(pldindex))
  }
  names(pldindex) <- GetAlleleNames(object)[allelesToExport]
  
  # loop through locus ploidies and fill in matrix
  for(p in ploidyindices){
    thesealleles <- which(pldindex == p)
    # loop through taxa ploidies
    for(h in seq_len(ncol(object$posteriorProb))){
      thispld <- dim(object$posteriorProb[[p,h]])[1] - 1L
      thesetaxa <- dimnames(object$posteriorProb[[p,h]])[[2]]
      outmat[thesetaxa,thesealleles] <-
        BestGenos(object$posteriorProb[[p,h]][,thesetaxa,allelesToExport[thesealleles]],
                  thispld, length(thesetaxa), length(thesealleles)) # Rcpp function
      # correct or delete genotypes that don't sum to ploidy
      if(multiallelic %in% c("na", "correct")){
        # fix up alleles2loc and determine number of loci
        thisA2L <- object$alleles2loc[thesealleles]
        theseLoci <- unique(thisA2L)
        thisA2L <- match(thisA2L, theseLoci)
        # run the Rcpp function do to the correction
        outmat[thesetaxa,thesealleles] <- 
          CorrectGenos(outmat[thesetaxa,thesealleles], 
                       object$posteriorProb[[p,h]][,thesetaxa,allelesToExport[thesealleles]],
                       thisA2L, length(thesetaxa), thispld, length(thesealleles),
                       length(theseLoci), multiallelic == "correct")
      }
    }
    
    # correct parent genotypes if this is a mapping population
    if(correctParentalGenos && !is.null(object$likelyGeno_donor) &&
       !is.null(object$likelyGeno_recurrent)){
      ### Fix this section when working through mapping pipeline
      pld.d <- sum(object$donorPloidies[[p]])
      pld.r <- sum(object$recurrentPloidies[[p]])
      outmat[GetDonorParent(object), thesealleles] <- 
        object$likelyGeno_donor[as.character(pld.d), allelesToExport[thesealleles]]
      outmat[GetRecurrentParent(object), thesealleles] <- 
        object$likelyGeno_recurrent[as.character(pld.r), allelesToExport[thesealleles]]
    }
  }
  
  # put in NA for zero reads if necessary
  if(naIfZeroReads){
    isZero <- (object$locDepth == 0)[,as.character(object$alleles2loc[allelesToExport])]
    outmat[isZero] <- NA_integer_
  }
  
  # subset if desired and if correction was done
  if(omit1allelePerLocus && multiallelic != "ignore"){
    allelesToExport <- allelesToExport[-OneAllelePerMarker(object,
                                                           commonAllele = omitCommonAllele)]
    outmat <- outmat[,allelesToExport]
    pldindex <- pldindex[allelesToExport]
  }
  
  return(list(genotypes = outmat, ploidy_index = pldindex))
}

AddPCA <- function(object, ...){
  UseMethod("AddPCA", object)
}
# some key additional arguments: nPcs is the number of PC axes to return
AddPCA.RADdata <- function(object, nPcsInit = 10, maxR2changeratio = 0.05, 
                           minPcsOut = 1, ...){
  if(minPcsOut > nPcsInit){
    stop("minPcsOut can not be greater than nPcsInit.")
  }
  # matrix for input to PCA; depth ratios or posterior probs
  if(!CanDoGetWeightedMeanGeno(object)){
    genmat <- object$depthRatio[,-OneAllelePerMarker(object)]
  } else {
    genmat <- GetWeightedMeanGenotypes(object, omit1allelePerLocus = TRUE,
                                       naIfZeroReads = FALSE)
  }
  
  # replace NaN with NA
  genmat[is.na(genmat)] <- NA
  # check that no columns are completely NA
  if(any(colMeans(is.na(genmat)) == 1)){
    message("Alleles with completely missing data:")
    cat(colnames(genmat)[colMeans(is.na(genmat)) == 1], sep = "\n")
    stop("Alleles found with completely missing data.")
  }
  
  genfreq <- colMeans(genmat, na.rm = TRUE)
  # if any individuals are completely missing, fill them in with mean.
  # This can happen with blanks or very low quality samples when looping
  # through the genome in chunks.
  missind <- which(rowMeans(is.na(genmat)) == 1)
  for(i in missind){
    genmat[i,] <- genfreq
  }
  
  # remove non-variable sites
  genmat <- genmat[, which(genfreq > 0 & genfreq < 1)]
  
  # adjust number of PC axes if necessary
  if(nPcsInit > dim(genmat)[2]){
    nPcsInit <- dim(genmat)[2]
  }
  
  # run principal components analysis
  pc <- pcaMethods::pca(genmat, method = "ppca", nPcs = nPcsInit, ...)
  # get rate of change in R2 values
  roc <- pc@R2[1:(nPcsInit - 1)] - pc@R2[2:nPcsInit]
  cutoff <- which(roc < roc[1] * maxR2changeratio)
  if(length(cutoff) == 0){
    cutoff <- nPcsInit
  }
  # make sure number of PCs meets minimum threshold
  if(cutoff[1] < minPcsOut){
    cutoff <- minPcsOut
  }
  
  object$PCA <- pc@scores[,1:cutoff[1]]

  return(object)
}

AddAlleleFreqByTaxa <- function(object, ...){
  UseMethod("AddAlleleFreqByTaxa")
}
AddAlleleFreqByTaxa.RADdata <- function(object, minfreq = 0.0001, ...){
  if(is.null(object$PCA)){
    stop("Need to run AddPCA first.")
  }
  if(minfreq <= 0){
    stop("minfreq must be greater than zero.")
  }
  if(minfreq >= 0.5){
    stop("minfreq must be less than 0.5 (typically much less).")
  }
  if(!CanDoGetWeightedMeanGeno(object)){
    genmat <- object$depthRatio
    genmat[is.na(genmat)] <- NA
  } else {
    genmat <- GetWeightedMeanGenotypes(object, omit1allelePerLocus = FALSE,
                                       minval = 0, maxval = 1,
                                       naIfZeroReads = FALSE)
  }
  if(sum(is.na(genmat)) == 0){
    # regress estimated genotypes on PC axes
    PCcoef <- lm(genmat ~ object$PCA)$coefficients
  } else {
    PCcoef <- matrix(NA, nrow = dim(object$PCA)[2] + 1, 
                     ncol = dim(genmat)[2])
    for(i in 1:dim(genmat)[2]){
      PCcoef[,i] <- lm(genmat[,i] ~ object$PCA)$coefficients
    }
  }
  PCcoef[is.na(PCcoef)] <- 0 # for non-variable sites
  
  # predict allele frequencies from PC axes
  predAl <- object$PCA %*% PCcoef[-1,] + 
    matrix(PCcoef[1,], byrow = TRUE, nrow = nTaxa(object), 
           ncol = nAlleles(object))
  
  # adjust allele frequencies to possible values
  predAl <- AdjustAlleleFreq(predAl, object$alleles2loc, minfreq)
  
  dimnames(predAl) <- list(GetTaxa(object), GetAlleleNames(object))
  object$alleleFreqByTaxa <- predAl
  return(object)
}

AddGenotypePriorProb_ByTaxa <- function(object, ...){
  UseMethod("AddGenotypePriorProb_ByTaxa", object)
}
AddGenotypePriorProb_ByTaxa.RADdata <- function(object, selfing.rate = 0, ...){
  if(is.null(object$alleleFreqByTaxa)){
    stop("Need to run AddAlleleFreqByTaxa first.")
  }
  tx_pld_unique <- sort(unique(GetTaxaPloidy(object)))
  priors <- array(list(),
                  dim = c(length(object$possiblePloidies),
                          length(tx_pld_unique)),
                  dimnames = list(NULL, as.character(tx_pld_unique)))
  
  for(i in seq_along(object$possiblePloidies)){
    for(j in seq_along(tx_pld_unique)){
      pldtot <- sum(object$possiblePloidies[[i]]) * tx_pld_unique[j] / 2L
      thesetaxa <- GetTaxaByPloidy(object, tx_pld_unique[j])
      priors[[i, j]] <- array(NA, dim = c(pldtot + 1,
                                       length(thesetaxa), nAlleles(object)),
                           dimnames = list(as.character(0:pldtot),
                                           thesetaxa,
                                           GetAlleleNames(object)))
      for(k in thesetaxa){
        priors[[i, j]][,k,] <- .HWEpriors(object$alleleFreqByTaxa[k,], 
                                       object$possiblePloidies[[i]] * tx_pld_unique[j] / 2L,
                                       selfing.rate)
      }
    }
  }
  
  object$priorProb <- priors
  object$priorProbPloidies <- object$possiblePloidies
  attr(object, "priorType") <- "taxon"
  return(object)
}

AddGenotypePriorProb_LD <- function(object, ...){
  UseMethod("AddGenotypePriorProb_LD", object)
}
AddGenotypePriorProb_LD.RADdata <- function(object, type, ...){
  if(is.null(object$posteriorProb) || is.null(object$alleleLinkages)){
    stop("posteriorProb and alleleLinkages slots needed.")
  }
  if(!type %in% c("mapping", "hwe", "popstruct")){
    stop("type must be 'mapping', 'hwe', or 'popstruct'.")
  }
  # Set up list of arrays to contain prior probabilities based on linked
  # loci alone.
  nPld <- length(object$posteriorProb)
  object$priorProbLD <- list()
  length(object$priorProbLD) <- nPld
  # Find parents, to exclude from correlations in mapping populations
  if(type == "mapping"){
    parents <- c(GetDonorParent(object), GetRecurrentParent(object))
    progeny <- which(!GetTaxa(object) %in% parents)
  }
  
  # Loop through the possible ploidies
  for(pldIndex in 1:nPld){
    # set up array to be put into priorProbLD slot.
    # (Setting up the array *in* the slot caused a mysterious problem
    # where object would get copied to a new memory address on each
    # iteration through the alleles, but only when the Rcpp function
    # ThirdDimProd was used instead of apply.)
    thisarr <- 
      array(dim = dim(object$posteriorProb[[pldIndex]]),
            dimnames = dimnames(object$posteriorProb[[pldIndex]]))
    # number of possible genotypes
    ngen <- dim(object$posteriorProb[[pldIndex]])[1]
    # Loop through alleles
    for(a in 1:nAlleles(object)){
      atab <- object$alleleLinkages[[a]]
      if(length(atab$allele) == 0){ # no linked alleles
        thisarr[,,a] <- 1/ngen
      } else {             # linked alleles exist
        # get posterior probabilities for linked alleles
        thispost <- object$posteriorProb[[pldIndex]][,, atab$allele, drop = FALSE]
        
        # in a mapping population, make sure we are considering genotypes that are possible
        if(type == "mapping"){
          # genotypes possible for this allele
          possibleThisAllele <- which(object$priorProb[[pldIndex]][,a] > 0)
          # array to hold new probabilities for getting priors
          newpost <- array(0, dim = dim(thispost))
          for(a2 in atab$allele){
            # genotypes possible for this linked allele
            possibleLinked <- which(object$priorProb[[pldIndex]][,a2] > 0)
            i <- match(a2, atab$allele)
            if(length(possibleThisAllele) == 2 && length(possibleLinked) == 2){
              # If there are only two possible genotypes for each, we know that all
              # copies of alleles are linked and can treat them that way.
              newpost[possibleThisAllele,progeny,i] <- 
                thispost[possibleLinked,progeny,i] * atab$corr[i]^2 +
                0.5 * (1 - atab$corr[i]^2)
            } else {
              # If any allele has more than two genotypes, we aren't sure of
              # complete phasing.
              # regress posterior probs for each allele copy number on all posterior probs
              for(j in possibleThisAllele){
                thisX <- thispost[possibleLinked[-1],, i]
                if(is.vector(thisX)){
                  thisX <- matrix(thisX, nrow = length(thisX), ncol = 1)
                } else {
                  thisX <- t(thisX)
                }
                thislm <- lm.fit(x = cbind(rep(1, length(progeny)), thisX[progeny,, drop=FALSE]),
                                 y = object$posteriorProb[[pldIndex]][j, progeny, a])
                newpost[j,progeny,i] <- thislm$fitted.values
              }
            }
          }
          newpost[newpost < 0] <- 0
          newpost[newpost > 1] <- 1
          thispost <- newpost
        } else { # for hwe or pop structure situations
          # multiply by correlation coefficient
          thispost <- sweep(thispost, 3, atab$corr^2, "*")
          # add even priors for the remainder of the coefficient
          thispost <- sweep(thispost, 3, (1 - atab$corr^2)/ngen, "+")
        }
        
        # multiply across alleles to get priors
        if(length(atab$allele) == 1){
          thisarr[,,a] <- thispost[,, 1]
        } else {
#          thisLDprior <- apply(thispost, c(1, 2), prod) # non-compiled version
          thisLDprior <- ThirdDimProd(thispost, ngen, nTaxa(object)) # Rcpp function
          thisLDprior <- sweep(thisLDprior, 2, colSums(thisLDprior), "/")
          thisarr[,,a] <- thisLDprior
        }
      } # end of chunk for if there are linked alleles
    } # end of loop through alleles
    object$priorProbLD[[pldIndex]] <- thisarr
  } # end of loop through ploidies
  
  return(object)
} # end of AddGenotypePriorProb_LD.RADdata

AddAlleleLinkages <- function(object, ...){
  UseMethod("AddAlleleLinkages", object)
}
AddAlleleLinkages.RADdata <- function(object, type, linkageDist, minCorr,
                                      excludeTaxa = character(0), ...){
  if(!type %in% c("mapping", "hwe", "popstruct")){
    stop("type must be 'mapping', 'hwe', or 'popstruct'.")
  }
  if(minCorr < 0){
    stop("minCorr for linkage disequilibrium cannot be negative.")
  }
  # get weighted mean genotypes for doing correlations
  wmgeno <- GetWeightedMeanGenotypes(object, omit1allelePerLocus = FALSE)
  wmgeno <- wmgeno[!rownames(wmgeno) %in% excludeTaxa, ]
  # set up new slot in RADdata object
  object$alleleLinkages <- list()
  length(object$alleleLinkages) <- nAlleles(object)
  # set up PCA values if they will be used
  if(type == "popstruct"){
    if(is.null(object$PCA)) stop("Run AddPCA first.")
    PCA <- object$PCA[!rownames(object$PCA) %in% excludeTaxa, , drop = FALSE]
  }
  
  # loop through loci
  for(L in 1:nLoci(object)){
    theseAlleles <- which(object$alleles2loc == L) # alleles for this locus
    nearbyAlleles <- FindNearbyAlleles(object, L, linkageDist)
    
    for(a in theseAlleles){
      # empty vectors to store alleles that are linked, and correlations
      linkedalleles <- integer(0)
      correlations <- numeric(0)
      # get weighted mean genotype values for this allele
      thisGenVal <- wmgeno[,a]
      if(type == "popstruct"){
        # get residuals after population structure accounted for
        thislm <- lm.fit(y = thisGenVal,
                         x = cbind(rep(1, nrow(PCA)), PCA))
        thisGenVal <- thislm$residuals
      }
      
      # see how each allele correlates with residuals
      for(a2 in nearbyAlleles){
        if(all(is.na(thisGenVal)) || all(thisGenVal == thisGenVal[1])) break
        if(all(is.na(wmgeno[,a2])) || all(wmgeno[,a2] == wmgeno[1,a2])) next
        thiscor <- cor(wmgeno[,a2], thisGenVal)
        if(thiscor >= minCorr){
          linkedalleles <- c(linkedalleles, a2)
          correlations <- c(correlations, thiscor)
        }
      }
      
      # store linkages for this allele
      object$alleleLinkages[[a]] <- list(allele = linkedalleles,
                                         corr = correlations)
    } # end of loop through alleles for one locus
  } # end of loop through loci
  
  return(object)
} # end of AddAlleleLinkages function

## The functions AddNormalizedDepthProp and AddAlleleBias are not currently
## used in polyRAD, but may be incorporated into future pipelines.
# Function to get a normalized proportion of reads belonging to each allele
# at each locus.
AddNormalizedDepthProp <- function(object, ...){
  UseMethod("AddNormalizedDepthProp", object)
}
AddNormalizedDepthProp.RADdata <- function(object, ...){
  # get the total number of reads per individual
  depthPerInd <- rowSums(object$locDepth)
  # get proportion of total reads for a taxon belonging to each allele
  propDepthAl <- sweep(object$alleleDepth, 1, depthPerInd, "/")
  totAl <- colSums(propDepthAl, na.rm = TRUE)
  # get proportion of total reads for a taxon belonging to other alleles at locus
  propDepthAnti <- sweep(object$antiAlleleDepth, 1, depthPerInd, "/")
  totAnti <- colSums(propDepthAnti, na.rm = TRUE)
  # get normalized proportion of reads for a locus belonging to an allele
  object$normalizedDepthProp <- totAl/(totAl + totAnti)
  
  return(object)
}
# function to estimate bias of each allele at each locus
# equivalent to h in Gerard et al. (2018)
AddAlleleBias <- function(object, ...){
  UseMethod("AddAlleleBias", object)
}
AddAlleleBias.RADdata <- function(object, maxbias = 4, ...){
  if(is.null(object$normalizedDepthProp)){
    object <- AddNormalizedDepthProp(object)
  }
  
  # estimate bias directly from data
  bias <- ((1 - object$normalizedDepthProp)/(1 - object$alleleFreq))/
    (object$normalizedDepthProp / object$alleleFreq)
  
  object$alleleBias <- bias
  
  return(object)
}

#### Accessors ####
GetTaxa <- function(object, ...){
  UseMethod("GetTaxa", object)
}
GetTaxa.RADdata <- function(object, ...){
  return(attr(object, "taxa"))
}
nTaxa <- function(object, ...){
  UseMethod("nTaxa", object)
}
nTaxa.RADdata <- function(object, ...){
  return(attr(object, "nTaxa"))
}
GetLoci <- function(object, ...){
  UseMethod("GetLoci", object)
}
GetLoci.RADdata <- function(object, ...){
  return(row.names(object$locTable))
}
nLoci <- function(object, ...){
  UseMethod("nLoci", object)
}
nLoci.RADdata <- function(object, ...){
  return(attr(object, "nLoci"))
}
nAlleles <- function(object, ...){
  UseMethod("nAlleles", object)
}
nAlleles.RADdata <- function(object, ...){
  return(dim(object$alleleDepth)[2])
}
GetAlleleNames <- function(object, ...){
  UseMethod("GetAlleleNames", object)
}
GetAlleleNames.RADdata <- function(object, ...){
  return(dimnames(object$alleleDepth)[[2]])
}
GetLocDepth <- function(object, ...){
  UseMethod("GetLocDepth", object)
}
GetLocDepth.RADdata <- function(object, ...){
  locnames <- row.names(object$locTable)[
    as.integer(dimnames(object$locDepth)[[2]])]
  outdepth <- object$locDepth
  dimnames(outdepth)[[2]] <- locnames
  return(outdepth)
}
GetContamRate <- function(object, ...){
  UseMethod("GetContamRate", object)
}
GetContamRate.RADdata <- function(object, ...){
  return(attr(object, "contamRate"))
}
SetContamRate <- function(object, value, ...){
  UseMethod("SetContamRate", object)
}
SetContamRate.RADdata <- function(object, value, ...){
  if(value < 0 || value > 1){
    stop("contamRate must range from zero to one.")
  }
  if(value > 0.01){
    warning("contamRate higher than expected.")
  }
  attr(object, "contamRate") <- value
  return(object)
}

SetTaxaPloidy <- function(object, value, ...){
  UseMethod("SetTaxaPloidy", object)
}
SetTaxaPloidy.RADdata <- function(object, value, ...){
  storage.mode(value) <- "integer"
  if(length(value) != nTaxa(object)){
    stop("Need one ploidy for each taxon.")
  }
  if(any(is.na(value))){
    stop("taxaPloidy must be integer.")
  }
  if(any(value < 1)){
    stop("taxaPloidy must be positive integer.")
  }
  if(!is.null(names(value))){
    if(!setequal(names(value), GetTaxa(object))){
      stop("Vector names do not match taxa names in object")
    }
    value <- value[GetTaxa(object)]
  } else {
    names(value) <- GetTaxa(object)
  }
  object$taxaPloidy <- value
  return(object)
}
GetTaxaPloidy <- function(object, ...){
  UseMethod("GetTaxaPloidy", object)
}
GetTaxaPloidy.RADdata <- function(object, ...){
  return(object$taxaPloidy)
}
GetTaxaByPloidy <- function(object, ...){
  UseMethod("GetTaxaByPloidy")
}
GetTaxaByPloidy.RADdata <- function(object, ploidy, ...){
  return(GetTaxa(object)[GetTaxaPloidy(object) == ploidy])
}

# Functions for assigning taxa to specific roles ####
SetDonorParent <- function(object, value){
  UseMethod("SetDonorParent", object)
}
SetDonorParent.RADdata <- function(object, value){
  if(!is.character(value) || length(value) != 1){
    stop("value must be one character string indicating the donor parent.")
  }
  if(!value %in% GetTaxa(object)){
    stop("value must be one of the taxa listed in the object.")
  }
  attr(object, "donorParent") <- value
  return(object)
}
GetDonorParent <- function(object, ...){
  UseMethod("GetDonorParent", object)
}
GetDonorParent.RADdata <- function(object, ...){
  dp <- attr(object, "donorParent")
  if(is.null(dp)){
    stop("Need to assign a donor parent with SetDonorParent.")
  }
  return(dp)
}
SetRecurrentParent <- function(object, value){
  UseMethod("SetRecurrentParent", object)
}
SetRecurrentParent.RADdata <- function(object, value){
  if(!is.character(value) || length(value) != 1){
    stop("value must be one character string indicating the recurrent parent.")
  }
  if(!value %in% GetTaxa(object)){
    stop("value must be one of the taxa listed in the object.")
  }
  attr(object, "recurrentParent") <- value
  return(object)
}
GetRecurrentParent <- function(object, ...){
  UseMethod("GetRecurrentParent", object)
}
GetRecurrentParent.RADdata <- function(object, ...){
  rp <- attr(object, "recurrentParent")
  if(is.null(rp)){
    stop("Need to assign a recurrent parent with SetRecurrentParent.")
  }
  return(rp)
}
SetBlankTaxa <- function(object, value){
  UseMethod("SetBlankTaxa", object)
}
SetBlankTaxa.RADdata <- function(object, value){
  if(!is.character(value)){
    stop("value must be a character vector")
  }
  if(!all(value %in% GetTaxa(object))){
    stop("Every element in value must be a taxon listed in the object.")
  }
  attr(object, "blankTaxa") <- value
  return(object)
}
GetBlankTaxa <- function(object, ...){
  UseMethod("GetBlankTaxa", object)
}
GetBlankTaxa.RADdata <- function(object, ...){
  bt <- attr(object, "blankTaxa")
  if(is.null(bt)){
    bt <- character(0)
  }
  return(bt)
}

# some basic utilities ####

# plot method that shows the PCA results
plot.RADdata <- function(x, ...){
  if(is.null(x$PCA)){
    message("Performing principal components analysis...")
    x <- AddPCA(x)
  }
  
  plot(x$PCA[,1], x$PCA[,2], xlab = "PC1", ylab = "PC2",
       main = "PCA of RADdata object", ...)
  if(!is.null(attr(x, "donorParent"))){
    don <- GetDonorParent(x)
    points(x$PCA[don, 1], x$PCA[don, 2], col = "darkgreen")
    text(x$PCA[don, 1], x$PCA[don, 2], "donor", col = "darkgreen", pos = 1)
  }
  if(!is.null(attr(x, "recurrentParent"))){
    rec <- GetRecurrentParent(x)
    points(x$PCA[rec, 1], x$PCA[rec, 2], col = "blue")
    text(x$PCA[rec, 1], x$PCA[rec, 2], "recurrent", col = "blue", pos = 1)
  }
}

OneAllelePerMarker <- function(object, ...){
  UseMethod("OneAllelePerMarker", object)
}
OneAllelePerMarker.RADdata <- function(object, commonAllele = FALSE, ...){
  if(commonAllele){
    # get index of most common allele for each marker
    if(is.null(object$alleleFreq)){
      # estimate allele frequencies if not done already
      object <- AddAlleleFreqHWE(object)
    }
    # sort by locus numbers, then by allele freq within locus
    myorder <- order(object$alleles2loc, object$alleleFreq, 
                     decreasing = c(FALSE, TRUE), method = "radix")
    # get position of most common allele within sorted vector
    mymatch1 <- fastmatch::fmatch(1:max(object$alleles2loc), 
                                  object$alleles2loc[myorder])
    mymatch1 <- na.omit(mymatch1)
    # translate that to position in unsorted vector
    mymatch <- myorder[mymatch1]
  } else {
    # get the index of the first allele for each marker
    mymatch <- fastmatch::fmatch(1:max(object$alleles2loc), object$alleles2loc)
    mymatch <- na.omit(mymatch)
  }
  return(mymatch)
}
CanDoGetWeightedMeanGeno <- function(object, ...){
  UseMethod("CanDoGetWeightedMeanGeno", object)
}
CanDoGetWeightedMeanGeno.RADdata <- function(object, ...){
  return(!is.null(object$posteriorProb) && 
           (!is.null(object$ploidyChiSq) || nrow(object$posteriorProb) == 1))
}

SubsetByTaxon <- function(object, ...){
  UseMethod("SubsetByTaxon", object)
}
SubsetByTaxon.RADdata <- function(object, taxa, ...){
  # check and convert taxa
  if(length(taxa) == 0){
    stop("At least one taxon must be provided for subsetting.")
  }
  if(is.character(taxa)){
    taxa <- fastmatch::fmatch(taxa, GetTaxa(object))
    if(any(is.na(taxa))) stop("Some taxa don't match RADdata object.")
  } else {
    if(any(is.na(taxa))) stop("No missing data allowed in taxa.")
  }
  if(!is.numeric(taxa)){
    stop("taxa must be a numeric or character vector")
  }
  
  # set up object and transfer attributes (including class)
  splitRADdata <- object
  attr(splitRADdata, "nTaxa") <- length(taxa)
  attr(splitRADdata, "taxa") <- GetTaxa(object)[taxa]
  
  # mandatory slots
  splitRADdata$alleleDepth <- object$alleleDepth[taxa, , drop = FALSE]
  splitRADdata$alleles2loc <- object$alleles2loc
  splitRADdata$locTable <- object$locTable
  splitRADdata$possiblePloidies <- object$possiblePloidies
  splitRADdata$locDepth <- object$locDepth[taxa, , drop = FALSE]
  splitRADdata$depthRatio <- object$depthRatio[taxa, , drop = FALSE]
  splitRADdata$antiAlleleDepth <- object$antiAlleleDepth[taxa, , drop = FALSE]
  splitRADdata$alleleNucleotides <- object$alleleNucleotides
  
  # slots that may have been added by other functions
  if(!is.null(object$depthSamplingPermutations)){
    splitRADdata$depthSamplingPermutations <- 
      object$depthSamplingPermutations[taxa, , drop = FALSE]
  }
  if(!is.null(object$alleleFreq)){
    splitRADdata$alleleFreq <- object$alleleFreq
  }
  if(!is.null(object$genotypeLikelihood)){
    splitRADdata$genotypeLikelihood <- 
      lapply(object$genotypeLikelihood, function(x) x[, taxa,, drop = FALSE])
  }
  if(!is.null(object$priorProb)){
    if(length(dim(object$priorProb[[1]])) == 3){
      splitRADdata$priorProb <- 
        lapply(object$priorProb, function(x) x[, taxa,, drop = FALSE])
    } else {
      splitRADdata$priorProb <- object$priorProb
    }
  }
  if(!is.null(object$priorProbPloidies)){
    splitRADdata$priorProbPloidies <- object$priorProbPloidies
  }
  if(!is.null(object$ploidyChiSq)){
    splitRADdata$ploidyChiSq <- object$ploidyChiSq
  }
  if(!is.null(object$ploidyChiSqP)){
    splitRADdata$ploidyChiSqP <- object$ploidyChiSqP
  }
  if(!is.null(object$posteriorProb)){
    splitRADdata$posteriorProb <- 
      lapply(object$posteriorProb, function(x) x[, taxa,, drop = FALSE])
  }
  if(!is.null(object$alleleFreqByTaxa)){
    splitRADdata$alleleFreqByTaxa <- object$alleleFreqByTaxa[taxa,, drop = FALSE]
  }
  if(!is.null(object$PCA)){
    splitRADdata$PCA <- object$PCA[taxa,, drop = FALSE]
  }
  if(!is.null(object$priorProbLD)){
    splitRADdata$priorProbLD <- 
      lapply(object$priorProbLD, function(x) x[, taxa,, drop = FALSE])
  }
  
  return(splitRADdata)
}

SubsetByLocus <- function(object, ...){
  UseMethod("SubsetByLocus", object)
}
SubsetByLocus.RADdata <- function(object, loci, ...){
  # check and convert loci
  if(length(loci) == 0){
    stop("At least one locus must be provided for subsetting.")
  }
  if(is.character(loci)){
    loci <- fastmatch::fmatch(loci, rownames(object$locTable))
    if(any(is.na(loci))) stop("Some loci don't match RADdata object.")
  } else {
    if(any(is.na(loci))) stop("No missing data allowed in loci.")
  }
  if(!is.numeric(loci)){
    stop("loci must be a numeric or character vector")
  }
  if(anyDuplicated(loci)){
    loci <- unique(loci)
    warning("Duplicate loci ignored.")
  }
  
  # set up object and transfer attributes (including class)
  thesealleles <- object$alleles2loc %fin% loci
  splitRADdata <- object
  attr(splitRADdata, "nLoci") <- length(loci)
  
  # mandatory slots
  splitRADdata$alleleDepth <- object$alleleDepth[, thesealleles, drop = FALSE]
  oldAl2loc <- object$alleles2loc[thesealleles]
  splitRADdata$alleles2loc <- fastmatch::fmatch(oldAl2loc, loci)
  splitRADdata$locTable <- object$locTable[loci,, drop = FALSE]
  splitRADdata$possiblePloidies <- object$possiblePloidies
  splitRADdata$locDepth <- object$locDepth[, as.character(loci), drop = FALSE]
  dimnames(splitRADdata$locDepth)[[2]] <- as.character(1:length(loci))
  splitRADdata$depthRatio <- object$depthRatio[, thesealleles, drop = FALSE]
  splitRADdata$antiAlleleDepth <- object$antiAlleleDepth[, thesealleles, drop = FALSE]
  splitRADdata$alleleNucleotides <- object$alleleNucleotides[thesealleles]
  attr(splitRADdata$alleleNucleotides, "Variable_sites_only") <- attr(object$alleleNucleotides, "Variable_sites_only")
  
  # additional components that may exist if some processing has already been done
  if(!is.null(object$depthSamplingPermutations)){
    splitRADdata$depthSamplingPermutations <- 
      object$depthSamplingPermutations[, thesealleles, drop = FALSE]
  }
  if(!is.null(object$alleleFreq)){
    splitRADdata$alleleFreq <- object$alleleFreq[thesealleles]
    attr(splitRADdata$alleleFreq, "type") <- attr(object$alleleFreq, "type")
  }
  if(!is.null(object$genotypeLikelihood)){
    splitRADdata$genotypeLikelihood <- 
      lapply(object$genotypeLikelihood, function(x) x[,, thesealleles, drop = FALSE])
  }
  if(!is.null(object$priorProb)){
    if(length(dim(object$priorProb[[1]])) == 3){
      splitRADdata$priorProb <- 
        lapply(object$priorProb, function(x) x[,, thesealleles, drop = FALSE])
    }
    if(length(dim(object$priorProb[[1]])) == 2){
      splitRADdata$priorProb <- 
        lapply(object$priorProb, function(x) x[, thesealleles, drop = FALSE])
    }
  }
  if(!is.null(object$priorProbPloidies)){
    splitRADdata$priorProbPloidies <- object$priorProbPloidies
  }
  if(!is.null(object$ploidyChiSq)){
    splitRADdata$ploidyChiSq <- object$ploidyChiSq[, thesealleles, drop = FALSE]
  }
  if(!is.null(object$ploidyChiSqP)){
    splitRADdata$ploidyChiSqP <- object$ploidyChiSqP[, thesealleles, drop = FALSE]
  }
  if(!is.null(object$posteriorProb)){
    splitRADdata$posteriorProb <- 
      lapply(object$posteriorProb, function(x) x[,, thesealleles, drop = FALSE])
  }
  if(!is.null(object$alleleFreqByTaxa)){
    splitRADdata$alleleFreqByTaxa <- object$alleleFreqByTaxa[, thesealleles, drop = FALSE]
  }
  if(!is.null(object$PCA)){
    splitRADdata$PCA <- object$PCA
  }
  if(!is.null(object$likelyGeno_donor)){
    splitRADdata$likelyGeno_donor <- 
      object$likelyGeno_donor[, thesealleles, drop = FALSE]
  }
  if(!is.null(object$likelyGeno_recurrent)){
    splitRADdata$likelyGeno_recurrent <- 
      object$likelyGeno_recurrent[, thesealleles, drop = FALSE]
  }
  if(!is.null(object$priorProbLD)){
    splitRADdata$priorProbLD <- 
      lapply(object$priorProbLD, function(x) x[,, thesealleles, drop = FALSE])
  }
  if(!is.null(object$alleleLinkages)){
    splitRADdata$alleleLinkages <- object$alleleLinkages[thesealleles]
  }
  
  return(splitRADdata)
}
SubsetByPloidy <- function(object, ...){
  UseMethod("SubsetByPloidy", object)
}
SubsetByPloidy.RADdata <- function(object, ploidies, ...){
  ploidies <- lapply(ploidies, as.integer)
  object$possiblePloidies <- ploidies
  # don't need to subset much if genotype calling not done
  if(is.null(object$priorProbPloidies)){
    return(object)
  }
  
  npld <- length(ploidies)
  # get indices of the desired ploidies in the object
  pldindex <- integer(npld)
  for(i in 1:npld){
    for(j in 1:length(object$priorProbPloidies)){
      if(identical(ploidies[[i]], object$priorProbPloidies[[j]])){
        pldindex[i] <- j
        break
      }
    }
    if(pldindex[i] == 0){
      stop(paste("Ploidy", paste(ploidies[[i]], collapse = " "), 
                 "not found in dataset."))
    }
  }
  object$priorProbPloidies <- ploidies
  
  if(!is.null(object$priorProb)){
    object$priorProb <- object$priorProb[pldindex]
  }
  if(!is.null(object$posteriorProb)){
    object$posteriorProb <- object$posteriorProb[pldindex]
  }
  if(!is.null(object$ploidyChiSq)){
    object$ploidyChiSq <- object$ploidyChiSq[pldindex,, drop = FALSE]
  }
  if(!is.null(object$ploidyChiSqP)){
    object$ploidyChiSqP <- object$ploidyChiSqP[pldindex,, drop = FALSE]
  }
  if(!is.null(object$ploidyLikelihood)){
    object$ploidyLikelihood <- object$ploidyLikelihood[pldindex,, drop = FALSE]
  }
  if(!is.null(object$priorTimesLikelihood)){
    object$priorTimesLikelihood <- object$priorTimesLikelihood[pldindex]
  }
  if(!is.null(object$priorProbLD)){
    object$priorProbLD <- object$priorProbLD[pldindex]
  }
  
  # subset items relating to likelihood, which ignore auto/allo differences
  pldsums <- sapply(ploidies, sum)
  if(!is.null(object$genotypeLikelihood)){
    likpld <- sapply(object$genotypeLikelihood, function(x) dim(x)[1] - 1L)
    object$genotypeLikelihood <- object$genotypeLikelihood[likpld %in% pldsums]
  }
  if(!is.null(object$likelyGeno_donor)){
    keeprow <- as.integer(rownames(object$likelyGeno_donor)) %in% pldsums
    object$likelyGeno_donor <- object$likelyGeno_donor[keeprow,, drop = FALSE]
  }
  if(!is.null(object$likelyGeno_recurrent)){
    keeprow <- as.integer(rownames(object$likelyGeno_recurrent)) %in% pldsums
    object$likelyGeno_recurrent <- 
      object$likelyGeno_recurrent[keeprow,, drop = FALSE]
  }
  
  # subset parental ploidies if mapping population
  ploidyfound <- function(x){
    for(pld in ploidies){
      if(identical(pld, x)) return(TRUE)
    }
    return(FALSE)
  }
  if(!is.null(object$donorPloidies)){
    keeppld <- sapply(object$donorPloidies, ploidyfound)
    object$donorPloidies <- object$donorPloidies[keeppld]
  }
  if(!is.null(object$recurrentPloidies)){
    keeppld <- sapply(object$recurrentPloidies, ploidyfound)
    object$recurrentPloidies <- object$recurrentPloidies[keeppld]
  }
  
  return(object)
}
SplitByChromosome <- function(object, ...){
  UseMethod("SplitByChromosome", object)
}
SplitByChromosome.RADdata <- function(object, chromlist = NULL, 
                                      chromlist.use.regex = FALSE, 
                                      fileprefix = "splitRADdata", ...){
  # set up the list of chromosomes if necessary
  if(is.null(chromlist)){
    chromlist <- unique(object$locTable$Chr)
    chromlist.use.regex <- FALSE
  }
  ngroups <- length(chromlist)
  
  # get vectors of locus numbers to keep for each item in chromlist
  locgroups <- list()
  length(locgroups) <- ngroups
  for(i in 1:ngroups){
    thesechr <- chromlist[[i]]
    if(chromlist.use.regex){ 
      # matching with regular expresions
      theseindices <- integer(0)
      for(chr in thesechr){
        theseindices <- c(theseindices, grep(chr, object$locTable$Chr))
      }
      theseindices <- sort(theseindices)
    } else { 
      # matching exactly
      theseindices <- which(object$locTable$Chr %fin% thesechr)
    }
    locgroups[[i]] <- theseindices
    message(paste(length(theseindices), "loci for chromosome(s)", 
                  paste(thesechr, collapse = " ")))
  }
  # check on total number of loci
  if(sum(sapply(locgroups, length)) > nLoci(object)){
    warning("Some loci are in multiple groups.")
  }
  
  # generate file names
  filenames <- paste(fileprefix, "_", 
                     sapply(chromlist, function(x) paste(x, collapse = "")),
                     ".RData", sep = "")
  
  # make new RADdata objects and save as .RData files.
  for(i in 1:ngroups){
    message(paste("Making new RADdata object for", 
                  paste(chromlist[[i]], collapse = " "), "..."))
    # build new object 
    splitRADdata <- SubsetByLocus(object, locgroups[[i]])

    # save the object to a file
    save(splitRADdata, file = filenames[i])
  }
  
  return(filenames)
}

# Function to discard slots that are no longer needed after pipelines have been
# run.
StripDown <- function(object, ...){
  UseMethod("StripDown", object)
}
StripDown.RADdata <- function(object, 
                              remove.slots = c("depthSamplingPermutations",
                                               "depthRatio", 
                                               "antiAlleleDepth",
                                               "genotypeLikelihood",
                                               "priorProb",
                                               "priorProbLD"),
                              ...){
  always.keep <- c("alleles2loc", "alleleNucleotides", "locTable", 
                   "priorProbPloidies", "possiblePloidies", "ploidyChiSq", 
                   "posteriorProb")
  if(any(always.keep %in% remove.slots)){
    ak <- always.keep[always.keep %in% remove.slots]
    stop(paste(c("Removal of the following would interfere with data export:", 
                 ak), 
               collapse = " "))
  }
  
  for(slot in remove.slots){
    object[[slot]] <- NULL
  }
  return(object)
}

# Function to find allele indices for nearby loci.
# locus can be the number or name of the locus.
# distance is the distance in basepairs within which to search.
# allele indices (not locus indices) are returned
FindNearbyAlleles <- function(object, ...){
  UseMethod("FindNearbyAlleles", object)
}
FindNearbyAlleles.RADdata <- function(object, locus, distance){
  if(!all(c("Chr", "Pos") %in% names(object$locTable))){
    stop("Alignment data not present in RADdata object.")
  }
  
  if(is.character(locus)){
    locus <- fastmatch::fmatch(locus, rownames(object$locTable))
  }
  thischr <- object$locTable[locus, "Chr"]
  thispos <- object$locTable[locus, "Pos"]
  
  # numbers for loci on this chromosome
  chrLocNum <- which(object$locTable$Chr == thischr)
  # positions on this chromosome
  chrPos <- object$locTable$Pos[chrLocNum]
  # numbers for loci near this locus
  matchingLocNum <- chrLocNum[chrPos >= thispos - distance & chrPos <= thispos + distance]
  ## (could redo this with GRanges)
  
  # remove this locus from matches
  matchingLocNum <- matchingLocNum[matchingLocNum != locus]
  
  # get matching alleles
  allelesout <- which(object$alleles2loc %fin% matchingLocNum)
  
  return(allelesout)
}

# function to find rare tags and merge them in with other tags from the locus.
MergeRareHaplotypes <- function(object, ...){
  UseMethod("MergeRareHaplotypes", object)
}
MergeRareHaplotypes.RADdata <- function(object, min.ind.with.haplotype = 10,
                                        ...){
  if(!is.null(object$alleleFreq) || !is.null(object$depthSamplingPermutations)){
    stop("Run MergeRareHaplotypes before running any pipeline functions.")
  }

  # find alleles that are rare
  Nindwithal <- colSums(object$alleleDepth > 0)
  rarealleles <- which(Nindwithal < min.ind.with.haplotype &
                         Nindwithal > 0)
  # find alleles that are completely absent and remove them
  zeroalleles <- which(Nindwithal == 0)
  remIndex <- length(zeroalleles) + 1L

  # numbers for corresponding loci
  lociToFix <- unique(object$alleles2loc[rarealleles])
  # preallocate for other alleles to remove
  allelesToRemove <- c(zeroalleles, integer(length(rarealleles)))

  # loop through loci that need fixing
  for(L in lociToFix){
    thesealleles <- which(object$alleles2loc == L)
    theserare <- thesealleles[thesealleles %fin% rarealleles]
    while(length(theserare) > 0 && length(thesealleles) > 1){
      thisAl <- theserare[1]
      # get nucleotide distance between this allele and others
      thesenuc <- object$alleleNucleotides[thesealleles]
      splitnuc <- strsplit(thesenuc, split = "")
      splitnucA <- splitnuc[[which(thesealleles == thisAl)]]
      nucdist <- integer(length(thesealleles))
      for(a in 1:length(thesealleles)){
        nucdist[a] <- sum(sapply(1:length(splitnucA),
                                 function(i) polyRADsubmat[splitnucA[i],
                                                           splitnuc[[a]][i]]))
      }
      # find the closest allele
      alToMerge <- thesealleles[nucdist == min(nucdist[-match(thisAl, thesealleles)])]
      if(length(alToMerge) > 1){
        # if multiple are closest, merge to the most common one
        alToMerge <- alToMerge[which.max(Nindwithal[alToMerge])]
      }
      if(length(alToMerge) > 1){
        # if still don't have one, pick at random
        alToMerge <- sample(alToMerge, 1)
      }
      # merge read counts
      object$alleleDepth[,alToMerge] <- object$alleleDepth[,alToMerge] +
        object$alleleDepth[,thisAl]
      object$antiAlleleDepth[,alToMerge] <-
        object$antiAlleleDepth[,alToMerge] - object$alleleDepth[,thisAl]
      object$depthRatio[,alToMerge] <- 
        object$depthRatio[,alToMerge] + object$depthRatio[,thisAl]
      Nindwithal[alToMerge] <- Nindwithal[alToMerge] + Nindwithal[thisAl]
      # merge nucleotides
      newNt <- .mergeNucleotides(object$alleleNucleotides[[alToMerge]],
                                 object$alleleNucleotides[[thisAl]])
      object$alleleNucleotides[[alToMerge]] <- newNt
      
      # remove this allele
      thesealleles <- thesealleles[thesealleles != thisAl]
      theserare <- thesealleles[Nindwithal[thesealleles] < 
                                  min.ind.with.haplotype]
      allelesToRemove[remIndex] <- thisAl
      remIndex <- remIndex + 1L
    } # end of while loop for merging rare alleles
    # if an insertion was discarded, get rid of the placeholder
    thesenuc <- object$alleleNucleotides[thesealleles]
    splitnuc <- strsplit(thesenuc, split = "")
    notplaceholder <- sapply(1:nchar(thesenuc[1]), 
                             function(i) any(sapply(splitnuc,
                                                    function(x) x[i] != '.')))
    if(any(!notplaceholder)){
      splitnuc <- lapply(splitnuc, function(x) x[notplaceholder])
      thesenuc <- sapply(splitnuc, function(x) paste(x, collapse = ""))
      object$alleleNucleotides[thesealleles] <- thesenuc
    }
  } # end of loop through loci
  
  # quit if there is nothing to remove
  if(remIndex == 1L){
    return(object)
  }
  
  allelesToRemove <- allelesToRemove[1:(remIndex - 1L)]
  
  # update and return RADdata object
  object$alleleDepth <- object$alleleDepth[,-allelesToRemove]
  object$antiAlleleDepth <- object$antiAlleleDepth[,-allelesToRemove]
  object$depthRatio <- object$depthRatio[,-allelesToRemove]
  varsiteonly <- attr(object$alleleNucleotides, "Variable_sites_only")
  object$alleleNucleotides <- object$alleleNucleotides[-allelesToRemove]
  attr(object$alleleNucleotides, "Variable_sites_only") <- varsiteonly
  object$alleles2loc <- object$alleles2loc[-allelesToRemove]
  return(object)
} # end of MergeRareHaplotypes

# Function to remove monomorphic markers
RemoveMonomorphicLoci <- function(object, ...){
  UseMethod("RemoveMonomorphicLoci", object)
}
RemoveMonomorphicLoci.RADdata <- function(object, verbose = TRUE, ...){
  alleleTally <- table(object$alleles2loc)
  locToKeep <- as.integer(names(alleleTally)[alleleTally > 1])
  
  if(verbose){
    message(paste(length(locToKeep), "markers retained out of", nLoci(object),
                  "originally."))
  }
  
  object <- SubsetByLocus(object, locToKeep)
  
  return(object)
}

# Function to remove high-depth markers (likely paralogs)
RemoveHighDepthLoci <- function(object, ...){
  UseMethod("RemoveHighDepthLoci", object)
}
RemoveHighDepthLoci.RADdata <- function(object, max.SD.above.mean = 2,
                                        verbose = TRUE, ...){
  # get total depth for each locus
  totdepth <- colSums(object$locDepth)
  stopifnot(all(!is.na(totdepth))) # there should never be NA values
  # mean and SD for depth
  meandepth <- mean(totdepth)
  sddepth <- sd(totdepth)
  # identify markers to keep
  cutoff <- meandepth + max.SD.above.mean * sddepth
  tokeep <- as.integer(names(totdepth)[totdepth <= cutoff])
  # subset object
  object <- SubsetByLocus(object, tokeep)
  
  if(verbose){
    message(paste(length(tokeep), "markers retained out of", length(totdepth),
                  "originally."))
    graphics::hist(totdepth/nTaxa(object), col = "lightgrey",
                   main = "Histogram of mean read depth", xlab = "Depth")
    graphics::abline(v = cutoff/nTaxa(object), col = "blue")
    graphics::text("Cutoff", x = cutoff/nTaxa(object), 
                   y = mean(graphics::par("yaxp")[1:2]),
                   col = "blue", pos = 2, srt = 90)
  }
  
  return(object)
}

# function to estimate contamination rate
EstimateContaminationRate <- function(object, ...){
  UseMethod("EstimateContaminationRate", object)
}
EstimateContaminationRate.RADdata <- function(object, multiplier = 1, ...){
  blanks <- GetBlankTaxa(object)
  if(length(blanks) == 0){
    stop("Run SetBlankTaxa before running EstimateContaminationRate.")
  }
  if(length(multiplier) > 1 && length(multiplier) != length(blanks)){
    stop("Need one value for multiplier, or one value per blank taxa.")
  }
  if(length(multiplier) > 1 && 
     (is.null(names(multiplier)) || !all(blanks %in% names(multiplier)))){
    stop("Names for multiplier vector should match names of blank taxa.")
  }
  
  # get depths
  nonblanks <- GetTaxa(object)[!GetTaxa(object) %in% blanks]
  meandepth <- mean(rowMeans(object$locDepth[nonblanks,, drop = FALSE]))
  
  blankdepth <- rowMeans(object$locDepth[blanks,, drop = FALSE])
  if(length(multiplier) == 1){
    blankdepth <- blankdepth * multiplier
  } else {
    blankdepth <- blankdepth * multiplier[blanks]
  }
  meanblankdepth <- mean(blankdepth)
  
  # get contamination rate and assign to object
  newcontam <- meanblankdepth/meandepth
  object <- SetContamRate(object, newcontam)
  
  cat(paste("Contamination rate:", newcontam), sep = "\n")

  return(object)
}

# function to return some useful info on the locus
LocusInfo <- function(object, ...){
  UseMethod("LocusInfo", object)
}
LocusInfo.RADdata <- function(object, locus, genome = NULL, annotation = NULL, verbose = TRUE, ...){
  if(length(locus) != 1){
    stop("LocusInfo function is designed for just one locus.")
  }
  # identify the number for the locus
  locnum <- fastmatch::fmatch(locus, GetLoci(object))
  # check if it is an allele name instead
  if(is.na(locnum)){
    alnum <- fastmatch::fmatch(locus, GetAlleleNames(object))
    if(is.na(alnum)){
      stop("'locus' does not match any locus or allele in 'object'.")
    }
    locnum <- object$alleles2loc[alnum]
  }
  
  # basic info on loci and alleles
  out <- list()
  out$Locus <- GetLoci(object)[locnum]
  if(verbose) cat(out$Locus, sep = "\n")
  aligned <- "Chr" %in% names(object$locTable)
  if(aligned){
    out$Chromosome <- object$locTable$Chr[locnum]
    out$Position <- object$locTable$Pos[locnum]
    if(verbose) cat(c(paste("Chromosome:", out$Chromosome),
                      paste("Position:", out$Position)), sep = "\n")
  }
  alnums <- which(object$alleles2loc == locnum)
  out$Alleles <- GetAlleleNames(object)[alnums]
  out$Haplotypes <- object$alleleNucleotides[alnums]
  if(verbose){
    cat(paste(length(out$Alleles), "alleles"), sep = "\n")
  }
  
  # allele frequencies if existing
  if(!is.null(object$alleleFreq)){
    out$Frequencies <- object$alleleFreq[alnums]
  }
  
  # functional annotation
  if(!is.null(annotation) && aligned){
    if(!requireNamespace("VariantAnnotation", quietly = TRUE)){
      stop("VariantAnnotation and other BioConductor packages needed if genome and annotation are provided.")
    }
    getfunc <- TRUE
    # get correct chromosome name
    chr <- out$Chr # chromosome name
    chrToMatch <- GenomeInfoDb::seqlevels(annotation) # chrom. in TxDb
    if(!chr %in% chrToMatch){
      chr2 <- grep(paste("^(Chr|chr|CHR)(omosome|OMOSOME)?0?", chr, "$", sep = ""),
                   chrToMatch, value = TRUE)
      if(length(chr2) == 1){
        chr <- chr2
      } else {
        warning("Could not match chromosome between RADdata and annotation.")
        getfunc <- FALSE # stop looking for functional annotation
      }
    }
    if(getfunc){
      nal <- length(out$Alleles)
      # set up GRanges for where this marker is
      gr <- GenomicRanges::GRanges(rep(chr, nal), 
                                   IRanges::IRanges(rep(out$Pos, nal), 
                                                    rep(out$Pos + nchar(out$Haplotypes[1]) - 1, nal)),
                                   strand = rep("+", nal))
      
      varsite <- attr(object$alleleNucleotides, "Variable_sites_only")
      if(((isTRUE(varsite) || is.null(varsite)) &&
          nchar(out$Haplotypes[1]) > 1) || 
         is.null(genome)){
        # warning that we can't get functional annotation of alleles
        if(is.null(genome)){
          warning("No reference genome provided; will not search for protein coding changes.")
        } else {
          warning("Haplotypes were imported as variable sites only; will not search for protein coding changes.")
        }
        # add genes overlapping the locus
        mygenes <- unique(GenomicFeatures::transcriptsByOverlaps(annotation, gr)$tx_name)
        out$Transcripts <- mygenes
        # print genes to console
        if(verbose){
          cat(c("Transcripts:", out$Transcripts), sep = "\n")
        }
      } else { # get functional annotation of alleles
        # alleles as DNAStrings
        ds <- Biostrings::DNAStringSet(out$Haplotypes)
        # predict coding changes
        message("Predicting protein coding changes...")
        out$PredictCoding <- VariantAnnotation::predictCoding(gr, annotation, genome, ds)
        # print results to console
        if(verbose && length(out$PredictCoding) > 0){
          cat(paste("Gene:", out$PredictCoding$GENEID[1]), sep = "\n")
        }
      }
    }
  } # end of search for functional annotation
  
  # print haplotypes
  if(verbose){
    stringsprint <- out$Haplotypes
    if(!is.null(out$PredictCoding) && nal == length(out$PredictCoding)){
      stringsprint <- paste(stringsprint, out$PredictCoding$CONSEQUENCE)
    }
    if(!is.null(out$Frequencies)){
      stringsprint <- paste(stringsprint, round(out$Frequencies, 2))
    }
    cat(stringsprint, sep = "\n")
  }
  
  return(out)
}

# function to merge multiple taxa into one before doing genotyping
MergeTaxaDepth <- function(object, ...){
  UseMethod("MergeTaxaDepth", object)
}
MergeTaxaDepth.RADdata <- function(object, taxa, ...){
  if(!is.null(object$alleleFreq) || !is.null(object$depthSamplingPermutations)){
    stop("Run MergeTaxaCounts before running any pipeline functions.")
  }
  if(length(taxa) < 2){
    stop("At least two taxa must be merged.")
  }
  message(paste(length(taxa), " taxa will be merged into the taxon ", taxa[1],
                ".", sep = ""))
  
  taxanum <- match(taxa, GetTaxa(object))
  if(any(is.na(taxanum))){
    stop(paste("Taxa not found in object:",
               paste(taxa[is.na(taxanum)], collapse = " ")))
  }
  
  # sum read depths and remove taxa from matrices
  object$alleleDepth[taxanum[1],] <- 
    as.integer(colSums(object$alleleDepth[taxanum,]))
  object$alleleDepth <- object$alleleDepth[-taxanum[-1],]
  
  object$locDepth[taxanum[1],] <- 
    as.integer(colSums(object$locDepth[taxanum,]))
  object$locDepth <- object$locDepth[-taxanum[-1],]
  
  object$antiAlleleDepth[taxanum[1],] <- 
    as.integer(colSums(object$antiAlleleDepth[taxanum,]))
  object$antiAlleleDepth <- object$antiAlleleDepth[-taxanum[-1],]
  
  # remove taxa from attributes
  attr(object, "taxa") <- attr(object, "taxa")[-taxanum[-1]]
  attr(object, "nTaxa") <- length(attr(object, "taxa"))
  
  # recalculated depth ratio and sampling permutations
  object$depthRatio <- object$depthRatio[-taxanum[-1],]
  object$depthRatio[taxa[1],] <- object$alleleDepth[taxa[1],]/
    (object$alleleDepth[taxa[1],] + object$antiAlleleDepth[taxa[1],])
  
  return(object)
}
# function to filter out loci that couldn't be genotyped, generally because
# parent genotypes did not match allele frequency.
RemoveUngenotypedLoci <- function(object, ...){
  UseMethod("RemoveUngenotypedLoci", object)
}
RemoveUngenotypedLoci.RADdata <- function(object, removeNonvariant = TRUE,
                                          ...){
  if(is.null(object$posteriorProb)){
    stop("Perform genotype calling before running RemoveUngenotypedLoci.")
  }
  
  # if this is a mapping population, we will exclude the parents when
  # looking for missing or non-variable genotypes
  have_parents <- !is.null(attr(object, "donorParent")) && 
    !is.null(attr(object, "recurrentParent"))
  if(have_parents){
    parents <- c(match(GetDonorParent(object), GetTaxa(object)),
                 match(GetRecurrentParent(object), GetTaxa(object)))
  }
  
  # vector of alleles to discard; becomes FALSE if non-missing for any ploidy
  alleles_discard <- rep(TRUE, nAlleles(object))
  
  # loop through ploidies
  for(prob in object$posteriorProb){
    if(have_parents){
      prob <- prob[, -parents, , drop = FALSE]
    }
    all_missing <- apply(prob, 3, function(x) all(is.na(x)))
    
    if(removeNonvariant){
      nmprob <- prob[,, !all_missing, drop = FALSE]
      nonvar <- apply(nmprob, 3, 
            function(x) all(apply(x, 1, 
                                  function(y) sd(y, na.rm = TRUE) == 0),
                            na.rm = TRUE))
      all_missing[!all_missing] <- nonvar
    }
    
    alleles_discard <- alleles_discard & all_missing
  }
  
  # discard loci that have any alleles to discard
  loci_discard <- unique(object$alleles2loc[alleles_discard])
  if(length(loci_discard) > 0){
    loci_keep <- (1:nLoci(object))[-loci_discard]
    object <- SubsetByLocus(object, loci_keep)
  }
  
  return(object)
}

# Function to merge alleles with identical nucleotides, generally because
# one tag was truncated with respect to the variable region.
MergeIdenticalHaplotypes <- function(object, ...){
  UseMethod("MergeIdenticalHaplotypes", object)
}
MergeIdenticalHaplotypes.RADdata <- function(object, ...){
  if(!is.null(object$alleleFreq) || !is.null(object$depthSamplingPermutations)){
    stop("Run MergeIdenticalHaplotypes before running any pipeline functions.")
  }
  
  remal <- integer(0) # indices of alleles to remove
  for(L in 1:nLoci(object)){
    thesecol <- which(object$alleles2loc == L)
    dup <- duplicated(object$alleleNucleotides[thesecol])
    for(al in thesecol[dup]){
      # find allele to merge this one into
      alM <- min(thesecol[object$alleleNucleotides[thesecol] == 
                            object$alleleNucleotides[al]])
      stopifnot(al != alM)
      # consolidate read depth
      object$alleleDepth[,alM] <- object$alleleDepth[,alM] +
        object$alleleDepth[,al]
      object$antiAlleleDepth[,alM] <- object$antiAlleleDepth[,alM] -
        object$alleleDepth[,al]
    }
    remal <- c(remal, thesecol[dup])
  }
  
  # remove duplicated alleles from all slots
  object$alleleDepth <- object$alleleDepth[,-remal]
  object$antiAlleleDepth <- object$antiAlleleDepth[,-remal]
  object$alleles2loc <- object$alleles2loc[-remal]
  varsiteonly <- attr(object$alleleNucleotides, "Variable_sites_only")
  object$alleleNucleotides <- object$alleleNucleotides[-remal]
  attr(object$alleleNucleotides, "Variable_sites_only") <- varsiteonly
  object$depthRatio <- object$depthRatio[,-remal]
  
  return(object)
}
