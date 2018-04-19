# class definition for storing data in polyRAD

# RADdata class constructor ####
RADdata <- function(alleleDepth, alleles2loc, locTable, possiblePloidies, 
                    contamRate, alleleNucleotides){
  if(!is.integer(alleleDepth)){
    stop("alleleDepth must be in integer format.")
  }
  if(!is.matrix(alleleDepth)){
    stop("alleleDepth must be in matrix format.")
  }
  if(any(is.na(alleleDepth))){
    stop("There should be no missing data in alleleDepth; put 0 for zero depth.")
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
    stop("loctable must be a data frame.")
  }
  if(!is.list(possiblePloidies)){
    stop("possiblePloidies must be list")
  }
  possiblePloidies <- lapply(possiblePloidies, as.integer)
  if(!all(sapply(possiblePloidies, function(x) all(!is.na(x) && x > 0)))){
    stop("Each element of possiblePloidies should be a vector of integers greater than zero.")
  }
  if(contamRate < 0 || contamRate > 1){
    stop("contamRate can only range from zero to one.")
  }
  if(contamRate > 0.01){
    warning("contamRate higher than expected.")
  }
  if(!(is.character(alleleNucleotides) | 
       "DNAStringSet" %in% class(alleleNucleotides))){
    stop("alleleNucleotides must be a character vector or DNAStringSet.")
  }
  if(length(alleleNucleotides) != length(alleles2loc)){
    stop("Length of alleleNucleotides must be same as length of alleles2loc.")
  }
  
  taxa <- dimnames(alleleDepth)[[1]]
  nTaxa <- dim(alleleDepth)[1]
  nLoci <- dim(locTable)[1]
  locDepth <- t(apply(alleleDepth, 1, function(x) tapply(x, alleles2loc, sum)))
     # dimnames(locDepth)[[2]] is integer from alleles2loc converted to character
  if(length(unique(alleles2loc)) == 1){
    locDepth <- matrix(locDepth, nrow = nTaxa, ncol = nLoci,
                       dimnames = list(taxa, as.character(unique(alleles2loc))))
  }
  
  expandedLocDepth <- locDepth[,as.character(alleles2loc), drop = FALSE]
  
  # get number of permutations of order in which each allele could have been sampled from total depth from that locus
  depthSamplingPermutations <- choose(expandedLocDepth, alleleDepth)
  dimnames(depthSamplingPermutations)[[2]] <- dimnames(alleleDepth)[[2]]
  # for each allele and taxon, get proportion of reads for that locus
  depthRatio <- alleleDepth/expandedLocDepth
  # depth of reads for each locus that do NOT belong to a given allele
  antiAlleleDepth <- expandedLocDepth - alleleDepth
  dimnames(antiAlleleDepth)[[2]] <- dimnames(alleleDepth)[[2]]
  
  # convert alleleNucleotides to DNAStringSet if Bioconductor installed
  if(requireNamespace("Biostrings", quietly = TRUE) && is.character(alleleNucleotides)){
    alleleNucleotides <- Biostrings::DNAStringSet(alleleNucleotides)
  }
  
  return(structure(list(alleleDepth = alleleDepth, alleles2loc = alleles2loc,
                        locTable = locTable, possiblePloidies = possiblePloidies,
                        locDepth = locDepth, 
                        depthSamplingPermutations = depthSamplingPermutations,
                        depthRatio = depthRatio, antiAlleleDepth = antiAlleleDepth,
                        alleleNucleotides = alleleNucleotides), 
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
  for(pl in x$possiblePloidies){
    cat(printPloidies(pl), sep = "\n")
  }
  if(!is.null(attr(x, "alleleFreqType"))){
    cat("", paste("Allele frequencies estimated for", attr(x, "alleleFreqType")),
        sep = "\n")
  }
}

#### parameter estimation generic functions and methods ####
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
  if(min(dist(expectedFreqs, method = "manhattan"))/2 < allowedDeviation){
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
AddGenotypeLikelihood.RADdata <- function(object, ...){
  if(is.null(object$alleleFreq)){
    cat("Allele frequencies not found; estimating under HWE from depth ratios.",
        sep = "\n")
    object <- AddAlleleFreqHWE(object)
  }
  alFreq <- object$alleleFreq
  
  # get ploidies, ignoring inheritance pattern
  ploidies <- sort(unique(sapply(object$possiblePloidies, sum)))
  # set up list for genotype likelihoods and loop through
  object$genotypeLikelihood <- list()
  length(object$genotypeLikelihood) <- length(ploidies)
  # probability of getting each allele from contamination
  sampleContam <- attr(object, "contamRate") * alFreq
  for(i in 1:length(ploidies)){
    # get probability of sampling each allele from each possible genotype
    sampleReal <- (0:ploidies[i])/ploidies[i] * (1 - attr(object, "contamRate"))
    alleleProb <- matrix(0, nrow = length(sampleReal), 
                         ncol = length(sampleContam))
    for(j in 1:length(sampleReal)){
      alleleProb[j,] <- sampleReal[j] + sampleContam
    }
    antiAlleleProb <- 1 - alleleProb
    
    # get likelihoods
    object$genotypeLikelihood[[i]] <- array(0, dim = c(ploidies[i]+1, 
                                                  dim(object$alleleDepth)),
                                       dimnames = list(as.character(0:ploidies[i]),
                                                       GetTaxa(object),
                                                       GetAlleleNames(object)))
    for(j in 1:(ploidies[i]+1)){
      object$genotypeLikelihood[[i]][j,,] <- 
        object$depthSamplingPermutations * 
          AlleleProbExp(object$alleleDepth, alleleProb[j,]) * 
            AlleleProbExp(object$antiAlleleDepth, antiAlleleProb[j,]) # Rcpp function
#          t(apply(object$alleleDepth, 1, function(x) alleleProb[j,] ^ x) * 
#            apply(object$antiAlleleDepth, 1, function(x) antiAlleleProb[j,] ^ x)) # non-C version of above two lines
      # when depth is too high, use dbinom instead.
      toRecalculate <- which(is.na(object$genotypeLikelihood[[i]][j,,]) |
                               object$genotypeLikelihood[[i]][j,,] == Inf,
                             arr.ind = TRUE)
      if(dim(toRecalculate)[1] == 0) next
      for(k in 1:dim(toRecalculate)[1]){
        # note repetitive code below
        taxon <- toRecalculate[k,1]
        allele <- toRecalculate[k,2]
        object$genotypeLikelihood[[i]][j,taxon,allele] <-
          dbinom(object$alleleDepth[taxon,allele],
                 object$locDepth[taxon, as.character(object$alleles2loc[allele])],
                 alleleProb[j,allele])
      }
    }
    # fix likelihoods where all are zero
    totlik <- colSums(object$genotypeLikelihood[[i]])
    toRecalculate <- which(totlik == 0, arr.ind = TRUE)
    if(dim(toRecalculate)[1] > 0){
      for(k in 1:dim(toRecalculate)[1]){
        for(j in 1:(ploidies[i] + 1)){
          # note repetitive code from just above
          taxon <- toRecalculate[k,1]
          allele <- toRecalculate[k,2]
          object$genotypeLikelihood[[i]][j,taxon,allele] <-
            dbinom(object$alleleDepth[taxon,allele],
                   object$locDepth[taxon, as.character(object$alleles2loc[allele])],
                   alleleProb[j,allele])
        }
        # for rare cases where likelihood still not estimated, set to one
        if(sum(object$genotypeLikelihood[[i]][, taxon, allele]) == 0){
          object$genotypeLikelihood[[i]][, taxon, allele] <- 1
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
    outmat[i,nonNaAlleles] <- apply(object$genotypeLikelihood[[i]][,taxind,nonNaAlleles], 2, which.max) - 1
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
                                   expectedFreq = possfreq)
  }
  if(is.null(object$genotypeLikelihood)){
    message("Genotype likelihoods not found.  Estimating.")
    object <- AddGenotypeLikelihood(object)
  }
  # get most likely genotype for the parents
  pldtot <- sapply(object$genotypeLikelihood, function(x) dim(x)[1] - 1)
  pldtot.don <- pldtot[pldtot %in% sapply(donorParentPloidies, sum)]
  pldtot.rec <- pldtot[pldtot %in% sapply(recurrentParentPloidies, sum)]
  nAlleles <- nAlleles(object)
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
  freqMatchGen <- matrix(FALSE, nrow = dim(pldcombos)[1], ncol = nAlleles)
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
      if(nrow(poss_matches) == 1){ # only one possible match (i.e. when there is backcrossing)
        likelyGen.don[as.character(pldcombos[i,"donor"]),a] <- unname(poss_matches[,1])
        likelyGen.rec[as.character(pldcombos[i,"recurrent"]),a] <- unname(poss_matches[,2])
      } else { # multiple possible matches
        current.don <- likelyGen.don[as.character(pldcombos[i,"donor"]),a]
        current.rec <- likelyGen.rec[as.character(pldcombos[i,"recurrent"]),a]
        # if one genotype is already a match, fix the other one
        if(!is.na(current.don) && any(poss_matches[,1] == current.don)){
          likelyGen.rec[as.character(pldcombos[i,"recurrent"]),a] <- unname(poss_matches[poss_matches[,1] == current.don,2])
        } else if(!is.na(current.rec) && any(poss_matches[,2] == current.rec)){
          likelyGen.don[as.character(pldcombos[i,"donor"]),a] <- unname(poss_matches[poss_matches[,2] == current.rec,1])
        }
      }
    }
  }
  
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

  # function to generate all gamete genotypes for a set of genotypes.
  # alCopy is a vector of values ranging from zero to ploidy indicating 
  # allele copy number.
  # ploidy is the ploidy
  # rnd indicates which round of the recursive algorithm we are on
  # Output is a matrix.  Alleles are in columns, which should be treated
  # independently.  Rows indicate gametes, with values indicating how many
  # copies of the allele that gamete has.
  makeGametes <- function(alCopy, ploidy, rnd = ploidy[1]/2){
    if(rnd %% 1 != 0 || ploidy[1] < 2){
      stop("Even numbered ploidy needed to simulate gametes.")
    }
    if(length(unique(ploidy)) != 1){
      stop("Currently all subgenome ploidies must be equal.")
    }
    if(any(alCopy > sum(ploidy), na.rm = TRUE)){
      stop("Cannot have alCopy greater than ploidy.")
    }
    if(length(ploidy) > 1){ # allopolyploids
      thisAl <- matrix(0L, nrow = 1, ncol = length(alCopy))
      for(pl in ploidy){
        # get allele copies for this isolocus.  Currently simplified; 
        # minimizes the number of isoloci to which an allele can belong.
        ### update in the future ###
        thisCopy <- alCopy
        thisCopy[thisCopy > pl] <- pl
        alCopy <- alCopy - thisCopy
        # make gametes for this isolocus (recursive, goes to autopoly version)
        thisIsoGametes <- makeGametes(thisCopy, pl, rnd)
        # add to current gamete set
        nGamCurr <- dim(thisAl)[1]
        nGamNew <- dim(thisIsoGametes)[1]
        thisAl <- thisAl[rep(1:nGamCurr, each = nGamNew),] + 
          thisIsoGametes[rep(1:nGamNew, times = nGamCurr),]
      }
    } else { # autopolyploids
      thisAl <- sapply(alCopy, function(x){
        if(is.na(x)){
          rep(NA, ploidy)
        } else {
          c(rep(0L, ploidy - x), rep(1L, x))
        }})
      if(rnd > 1){
        # recursively add alleles to gametes for polyploid
        nReps <- factorial(ploidy-1)/factorial(ploidy - rnd)
        thisAl <- thisAl[rep(1:ploidy, each = nReps),]
        for(i in 1:ploidy){
          thisAl[((i-1)*nReps+1):(i*nReps),] <- 
            thisAl[((i-1)*nReps+1):(i*nReps),] + 
            makeGametes(alCopy - thisAl[i*nReps,], ploidy - 1, rnd - 1)
        }
      }
    }
    return(thisAl)
  }
  # function to get probability of a gamete with a given allele copy number,
  # given output from makeGametes
  gameteProb <- function(makeGamOutput, ploidy){
    return(t(sapply(0:(sum(ploidy)/2), 
                    function(x) colMeans(makeGamOutput == x))))
  }
  # function to take two sets of gamete probabilities (from two parents)
  # and output genotype probabilities
  progenyProb <- function(gamProb1, gamProb2){
    outmat <- matrix(0L, nrow = dim(gamProb1)[1] + dim(gamProb2)[1] - 1,
                     ncol = dim(gamProb1)[2])
    for(i in 1:dim(gamProb1)[1]){
      copy1 <- i - 1
      for(j in 1:dim(gamProb2)[1]){
        copy2 <- j - 1
        thisrow <- copy1 + copy2 + 1
        outmat[thisrow,] <- outmat[thisrow,] + gamProb1[i,] * gamProb2[j,]
      }
    }
    return(outmat)
  }
  # function to get gamete probabilities for a population, given genotype priors
  gameteProbPop <- function(priors, ploidy){
    # get gamete probs for all possible genotypes
    possGamProb <- gameteProb(makeGametes(1:dim(priors)[1] - 1, ploidy),ploidy)
    # output matrix
    outmat <- possGamProb %*% priors
    return(outmat)
  }
  # function to adjust genotype probabilities from one generation of selfing
  selfPop <- function(priors, ploidy){
    # get gamete probs for all possible genotypes
    possGamProb <- gameteProb(makeGametes(1:dim(priors)[1] - 1, ploidy),ploidy)
    # progeny probs for all possible genotypes, selfed
    possProgenyProb <- matrix(0, nrow = dim(priors)[1], 
                              ncol = dim(possGamProb)[2])
    for(i in 1:dim(possGamProb)[2]){
      possProgenyProb[,i] <- progenyProb(possGamProb[,i,drop = FALSE],
                                         possGamProb[,i,drop = FALSE])
    }
    # multiple progeny probs by prior probabilities of those genotypes
    outmat <- possProgenyProb %*% priors
    return(outmat)
  }
  # get prior genotype probabilities for F1
  OutPriors <- list()
  length(OutPriors) <- dim(pldcombosExpand)[1]
  for(i in 1:length(OutPriors)){
    donorPld <- object$possiblePloidies[[pldcombosExpand[i,"donor"]]]
    recurPld <- object$possiblePloidies[[pldcombosExpand[i,"recurrent"]]]
    theseDonorGen <- likelyGen.don[as.character(sum(donorPld)),]
    theseRecurGen <- likelyGen.rec[as.character(sum(recurPld)),]
    donorGamProb <- gameteProb(makeGametes(theseDonorGen, donorPld), donorPld)
    recurGamProb <- gameteProb(makeGametes(theseRecurGen, recurPld), recurPld)
    OutPriors[[i]] <- progenyProb(donorGamProb, recurGamProb)
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
      recurGamProb <- gameteProb(makeGametes(theseRecurGen, recurPld), recurPld)
      # get gamete prob for current population
      currPld <- object$possiblePloidies[[pldcombosExpand[i,paste("BC", gen, sep = "")]]]
      currGamProb <- gameteProbPop(OutPriorsLastGen[[i]], currPld)
      # update genotype priors
      OutPriors[[i]] <- progenyProb(recurGamProb, currGamProb)
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
      currGamProb <- gameteProbPop(OutPriorsLastGen[[i]], currPld)
      OutPriors[[i]] <- progenyProb(currGamProb, currGamProb)
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
      OutPriors[[i]] <- selfPop(OutPriorsLastGen[[i]], currPld)
    }
  }
  
  for(i in 1:length(OutPriors)){
    dimnames(OutPriors[[i]]) <- list(as.character(0:(dim(OutPriors[[i]])[1] - 1)),
                                     GetAlleleNames(object))
  }

  object$priorProb <- OutPriors
  object$priorProbPloidies <- object$possiblePloidies[pldcombosExpand[,"final"]]
  attr(object, "priorType") <- "population" 
  # --> indicates prior probs are estimated for whole pop, not by taxa
  return(object)
}

AddGenotypePriorProb_HWE <- function(object, ...){
  UseMethod("AddGenotypePriorProb_HWE", object)
}
AddGenotypePriorProb_HWE.RADdata <- function(object, ...){
  if(is.null(object$alleleFreq) || attr(object, "alleleFreqType") != "HWE"){
    stop("Allele frequencies not estimated under HWE.")
  }
  priors <- list()
  length(priors) <- length(object$possiblePloidies)
  
  for(i in 1:length(priors)){
    priors[[i]] <- .HWEpriors(object$alleleFreq, object$possiblePloidies[[i]])
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
  for(i in 1:length(object$priorProb)){
    thisploidy <- dim(object$priorProb[[i]])[1] - 1
    thesegen <- sapply(likgen, function(x) x[as.character(thisploidy),])
    countstable <- matrix(0L, nrow = thisploidy + 1, ncol = dim(thesegen)[1],
                          dimnames = list(as.character(0:thisploidy),
                                          dimnames(thesegen)[[1]]))
    for(j in 0:thisploidy){
      countstable[j + 1, ] <- rowSums(thesegen == j, na.rm = TRUE)
    }
    # get ploidy likelihood
    thislikehd <- sapply(1:nAlleles, function(x){
      if(any(is.na(object$priorProb[[i]][,x]))){
        NA
      } else {
        dmultinom(countstable[,x], 
                  prob = object$priorProb[[i]][,x])
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
  object$ploidyChiSq <- matrix(NA, nrow = length(object$priorProb),
                               ncol = nAllele,
                               dimnames = list(NULL, GetAlleleNames(object)))
  object$ploidyChiSqP <- matrix(NA, nrow = length(object$priorProb),
                                ncol = nAllele,
                                dimnames = list(NULL, GetAlleleNames(object)))
  
  # get weighted genotype tallies from genotype likelihoods
  gental <- list()
  length(gental) <- length(object$genotypeLikelihood)
  for(i in 1:length(gental)){
    # likelihood total for each individual and locus at this ploidy
    totlik <- colSums(object$genotypeLikelihood[[i]][,taxa,])
    # normalize likelihoods by total for each individual and locus
    normlik <- sweep(object$genotypeLikelihood[[i]][,taxa,],
                     2:3, totlik, FUN = "/")
    # get population proportion of genotype likelihoods for allele and copy number
    gental[[i]] <- apply(normlik, c(1,3), mean)
  }
  
  # loop through ploidies
  for(i in 1:length(object$priorProb)){
    thisploidy <- dim(object$priorProb[[i]])[1] - 1
    whichlik <- which(sapply(object$genotypeLikelihood, 
                             function(x) dim(x)[1] - 1) == thisploidy)
    stopifnot(length(whichlik) == 1)
    # get priors
    if(attr(object, "priorType") == "population"){
      thesepriors <- object$priorProb[[i]]
    } else {
      # convert priors by taxon to population priors
      thesepriors <- rowMeans(aperm(object$priorProb[[i]], c(1,3,2)), dims = 2)
    }
    # estimate the components that are summed to make chi square
    chisqcomp <- (gental[[whichlik]] - thesepriors)^2/
      thesepriors * length(taxa)
    # degrees of freedom
    theseDF <- colSums(thesepriors != 0) - 1
    # chi-squared statistic
    thesechisq <- apply(chisqcomp, 2, function(x) sum(x[x != Inf]))
    object$ploidyChiSq[i,] <- thesechisq
    # p-values
    object$ploidyChiSqP[i,] <- pchisq(thesechisq, theseDF, lower.tail = FALSE)
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
  object$posteriorProb <- list()
  length(object$posteriorProb) <- length(PTL)
  for(i in 1:length(object$posteriorProb)){
    totPriorTimesLikeli <- colSums(PTL[[i]])
    object$posteriorProb[[i]] <- sweep(PTL[[i]], c(2,3),
                                       totPriorTimesLikeli, FUN = "/")
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
  
  nPloidies <- length(object$priorProb)
  
  # get weights for ploidies to use for each allele
  if(is.null(object$ploidyChiSq)){
    ploidyweights <- matrix(1, nrow = 1, ncol = length(altokeep))
  } else {
    nPloidies <- dim(object$ploidyChiSq)[1]
    if(onePloidyPerAllele){
      ploidyweights <- matrix(0, nrow = nPloidies,
                              ncol = length(altokeep))
      bestploidies <- apply(object$ploidyChiSq[,altokeep, drop = FALSE], 2, 
                            function(x){
                              if(all(is.na(x))){
                                   return(0)
                                 } else {
                                   return(which.min(x))
                                 }})
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
  wmgeno <- matrix(0, nrow = nTaxa(object),
                   ncol = length(altokeep),
                   dimnames = list(GetTaxa(object),
                                   GetAlleleNames(object)[altokeep]))
  # loop through ploidies
  for(i in 1:nPloidies){
    # values to represent each allele copy number
    thesegenval <- seq(minval, maxval, length.out = dim(object$posteriorProb[[i]])[1])
    # weighted mean genotypes for this ploidy
    thesewm <- rowSums(aperm(sweep(object$posteriorProb[[i]][,,altokeep],
                                   1, thesegenval, "*"), 
                             c(2,3,1)), 
                       dims = 2)
    thesewm[is.na(thesewm)] <- 0
    # multiply by weight for this ploidy and add to total
    wmgeno <- wmgeno + sweep(thesewm, 2, ploidyweights[i,], "*")
  }
  
  if(naIfZeroReads){ # if there were zero reads for that locus, replace with NA
    wmgeno[object$locDepth[,as.character(object$alleles2loc)[altokeep]] == 0] <- NA
  }
  
  return(wmgeno)
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
  # remove non-variable sites
  genfreq <- colMeans(genmat, na.rm = TRUE)
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
  for(loc in sort(unique(object$alleles2loc))){
    thesecol <- which(object$alleles2loc == loc) # columns for this locus
    for(taxa in 1:nTaxa(object)){
      thesefreq <- predAl[taxa,thesecol]
      thesefreq <- thesefreq/sum(thesefreq) # must sum to 1
      while(any(thesefreq < minfreq)){
        oldfreq <- thesefreq
        toolow <- thesefreq < minfreq
        thesefreq[toolow] <- minfreq
        canadjust <- thesefreq > minfreq
        adjust <- sum(thesefreq - oldfreq)/sum(canadjust)
        thesefreq[canadjust] <- thesefreq[canadjust] - adjust
      }
      while(any(thesefreq > (1 - minfreq))){
        oldfreq <- thesefreq
        toohigh <- thesefreq > 1 - minfreq
        thesefreq[toohigh] <- 1 - minfreq
        canadjust <- thesefreq < 1 - minfreq
        adjust <- sum(oldfreq - thesefreq)/sum(canadjust)
        thesefreq[canadjust] <- thesefreq[canadjust] + adjust
      }
      predAl[taxa,thesecol] <- thesefreq
    }
  }
  
  dimnames(predAl) <- list(GetTaxa(object), GetAlleleNames(object))
  object$alleleFreqByTaxa <- predAl
  return(object)
}

AddGenotypePriorProb_ByTaxa <- function(object, ...){
  UseMethod("AddGenotypePriorProb_ByTaxa", object)
}
AddGenotypePriorProb_ByTaxa.RADdata <- function(object, ...){
  if(is.null(object$alleleFreqByTaxa)){
    stop("Need to run AddAlleleFreqByTaxa first.")
  }
  priors <- list()
  length(priors) <- length(object$possiblePloidies)
  
  for(i in 1:length(priors)){
    pldtot <- sum(object$possiblePloidies[[i]])
    priors[[i]] <- array(NA, dim = c(pldtot + 1,
                                     nTaxa(object), nAlleles(object)),
                         dimnames = list(as.character(0:pldtot),
                                         GetTaxa(object),
                                         GetAlleleNames(object)))
    for(j in 1:nTaxa(object)){
      priors[[i]][,j,] <- .HWEpriors(object$alleleFreqByTaxa[j,], 
                                     object$possiblePloidies[[i]])
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
    # set up array
    object$priorProbLD[[pldIndex]] <- 
      array(dim = dim(object$posteriorProb[[pldIndex]]),
            dimnames = dimnames(object$posteriorProb[[pldIndex]]))
    # number of possible genotypes
    ngen <- dim(object$posteriorProb[[pldIndex]])[1]
    # Loop through alleles
    for(a in 1:nAlleles(object)){
      atab <- object$alleleLinkages[[a]]
      if(length(atab$allele) == 0){ # no linked alleles
        object$priorProbLD[[pldIndex]][,,a] <- 1/ngen
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
                thispost[possibleLinked,progeny,i] * atab$corr[i] +
                0.5 * (1 - atab$corr[i])
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
          thispost <- sweep(thispost, 3, atab$corr, "*")
          # add even priors for the remainder of the coefficient
          thispost <- sweep(thispost, 3, (1 - atab$corr)/ngen, "+")
        }
        
        # multiply across alleles to get priors
        if(length(atab$allele) == 1){
          object$priorProbLD[[pldIndex]][,,a] <- thispost[,, 1]
        } else {
          thisLDprior <- apply(thispost, c(1, 2), prod)
          thisLDprior <- sweep(thisLDprior, 2, colSums(thisLDprior), "/")
          object$priorProbLD[[pldIndex]][,,a] <- thisLDprior
        }
      } # end of chunk for if there are linked alleles
    } # end of loop through alleles
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
           (!is.null(object$ploidyChiSq) || length(object$posteriorProb) == 1))
}

SubsetByTaxon <- function(object, ...){
  UseMethod("SubsetByTaxon", object)
}
SubsetByTaxon.RADdata <- function(object, taxa, ...){
  # check and convert taxa
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
  splitRADdata <- list()
  oldAttributes <- attributes(object)
  oldAttributes <- oldAttributes[-match("names", names(oldAttributes))]
  attributes(splitRADdata) <- oldAttributes
  attr(splitRADdata, "nTaxa") <- length(taxa)
  attr(splitRADdata, "taxa") <- GetTaxa(object)[taxa]
  
  # mandatory slots
  splitRADdata$alleleDepth <- object$alleleDepth[taxa, , drop = FALSE]
  splitRADdata$alleles2loc <- object$alleles2loc
  splitRADdata$locTable <- object$locTable
  splitRADdata$possiblePloidies <- object$possiblePloidies
  splitRADdata$locDepth <- object$locDepth[taxa, , drop = FALSE]
  splitRADdata$depthSamplingPermutations <- 
    object$depthSamplingPermutations[taxa, , drop = FALSE]
  splitRADdata$depthRatio <- object$depthRatio[taxa, , drop = FALSE]
  splitRADdata$antiAlleleDepth <- object$antiAlleleDepth[taxa, , drop = FALSE]
  splitRADdata$alleleNucleotides <- object$alleleNucleotides
  
  # slots that may have been added by other functions
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
  
  return(splitRADdata)
}

SubsetByLocus <- function(object, ...){
  UseMethod("SubsetByLocus", object)
}
SubsetByLocus.RADdata <- function(object, loci, ...){
  # check and convert loci
  if(is.character(loci)){
    loci <- fastmatch::fmatch(loci, rownames(object$locTable))
    if(any(is.na(loci))) stop("Some loci don't match RADdata object.")
  } else {
    if(any(is.na(loci))) stop("No missing data allowed in loci.")
  }
  if(!is.numeric(loci)){
    stop("loci must be a numeric or character vector")
  }
  
  # set up object and transfer attributes (including class)
  thesealleles <- object$alleles2loc %fin% loci
  splitRADdata <- list()
  oldAttributes <- attributes(object)
  oldAttributes <- oldAttributes[-match("names", names(oldAttributes))]
  attributes(splitRADdata) <- oldAttributes
  attr(splitRADdata, "nLoci") <- length(loci)
  
  # mandatory slots
  splitRADdata$alleleDepth <- object$alleleDepth[, thesealleles, drop = FALSE]
  oldAl2loc <- object$alleles2loc[thesealleles]
  splitRADdata$alleles2loc <- fastmatch::fmatch(oldAl2loc, loci)
  splitRADdata$locTable <- object$locTable[loci,, drop = FALSE]
  splitRADdata$possiblePloidies <- object$possiblePloidies
  splitRADdata$locDepth <- object$locDepth[, as.character(loci), drop = FALSE]
  dimnames(splitRADdata$locDepth)[[2]] <- as.character(1:length(loci))
  splitRADdata$depthSamplingPermutations <- 
    object$depthSamplingPermutations[, thesealleles, drop = FALSE]
  splitRADdata$depthRatio <- object$depthRatio[, thesealleles, drop = FALSE]
  splitRADdata$antiAlleleDepth <- object$antiAlleleDepth[, thesealleles, drop = FALSE]
  splitRADdata$alleleNucleotides <- object$alleleNucleotides[thesealleles]
  
  # additional components that may exist if some processing has already been done
  if(!is.null(object$alleleFreq)){
    splitRADdata$alleleFreq <- object$alleleFreq[thesealleles]
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
  
  return(splitRADdata)
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
