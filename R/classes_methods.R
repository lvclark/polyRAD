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
  if(!is.character(alleleNucleotides)){
    stop("alleleNucleotides must be a character vector.")
  }
  if(length(alleleNucleotides) != length(alleles2loc)){
    stop("Length of alleleNucleotides must be same as length of alleles2loc.")
  }
  
  taxa <- dimnames(alleleDepth)[[1]]
  nTaxa <- dim(alleleDepth)[1]
  nLoci <- dim(locTable)[1]
  locDepth <- t(apply(alleleDepth, 1, function(x) tapply(x, alleles2loc, sum)))
     # dimnames(locDepth)[[2]] is integer from alleles2loc converted to character
  
  # get number of permutations of order in which each allele could have been sampled from total depth from that locus
  expandedLocDepth <- locDepth[,as.character(alleles2loc)]
  depthSamplingPermutations <- choose(expandedLocDepth, alleleDepth)
  # for each allele and taxon, get proportion of reads for that locus
  depthRatio <- alleleDepth/expandedLocDepth
  # depth of reads for each locus that do NOT belong to a given allele
  antiAlleleDepth <- expandedLocDepth - alleleDepth
  
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
      paste(attr(x, "nTaxa"), "taxa and", attr(x, "nLoci"), "loci"),
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
                                                         GetBlankTaxa(object))){#,
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
  taxaToKeep <- !GetTaxa(object) %in% excludeTaxa
  meanRat <- colMeans(object$depthRatio[taxaToKeep,], na.rm = TRUE)
  outFreq <- rep(NA, length(meanRat))
  names(outFreq) <- names(meanRat)
  for(f in expectedFreqs){
    outFreq[which(meanRat >= f - allowedDeviation & meanRat < f + allowedDeviation)] <- f
  }
  
  # fill in any missing allele frequencies for partially completed loci
  incompleteLoci <- unique(object$allele2loc[is.na(outFreq)])
  for(L in incompleteLoci){
    theseal <- which(allele2loc == L)
    if(sum(is.na(outFreq[theseal])) == 1){
      subtot <- sum(outFreq[theseal], na.rm = TRUE)
      outFreq[theseal[is.na(outFreq[theseal])]] <- 1 - subtot
    }
  }
  
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
AddAlleleFreqHWE.RADdata <- function(object, ...){
  meanRat <- colMeans(object$depthRatio, na.rm = TRUE)
  object$alleleFreq <- meanRat
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
  if(is.matrix(object$alleleFreq)){
    # If allele frequency was estimated seperately for each taxon, get mean
    # across the whole sample, since contamination could come from any taxon.
    alFreq <- colMeans(object$alleleFreq, na.rm = TRUE)
  } else {
    alFreq <- object$alleleFreq
  }
  
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
                                                       attr(object, "taxa"),
                                                       dimnames(object$alleleDepth)[[2]]))
    for(j in 1:(ploidies[i]+1)){
      object$genotypeLikelihood[[i]][j,,] <- 
        object$depthSamplingPermutations * 
          t(apply(object$alleleDepth, 1, function(x) alleleProb[j,] ^ x) * 
            apply(object$antiAlleleDepth, 1, function(x) antiAlleleProb[j,] ^ x))
    }
  }
  
  return(object)
}

# for a given taxon, return the most likely genotypes (based on likelihoods only)
GetLikelyGen <- function(object, taxon){
  UseMethod("GetLikelyGen", object)
}
GetLikelyGen.RADdata <- function(object, taxon){
  if(length(taxon) != 1){
    stop("Only one taxon can be passed to GetLikelyGen.")
  }
  taxind <- match(taxon, GetTaxa(object))
  if(is.na(taxind)){
    stop("taxon not found in object.")
  }
  if(is.null(object$genotypeLikelihood)){
    cat("Genotype likelihoods not found.  Estimating.", sep = "\n")
    object <- AddGenotypeLikelihood(object)
  }
  npld <- length(object$genotypeLikelihood)
  ploidies <- sapply(object$genotypeLikelihood, function(x) dim(x)[1] - 1)
  nAllele <- dim(object$alleleDepth)[2]
  
  outmat <- matrix(NA_integer_, nrow = npld, ncol = nAllele,
                   dimnames = list(as.character(ploidies), 
                                   dimnames(object$alleleDepth)[[2]]))
  naAlleles <- which(is.na(object$genotypeLikelihood[[i]][1,taxind,]))
  for(i in 1:npld){
    outmat[i,-naAlleles] <- apply(object$genotypeLikelihood[[i]][,taxind,-naAlleles], 2, which.max) - 1
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
    n.gen.selfing = 0, donorParentPloidies = object$possiblePloidies,
    recurrentParentPloidies = object$possiblePloidies){
  if(any(!donorParentPloidies %in% object$possiblePloidies) ||
     any(!recurrentParentPloidies %in% object$possiblePloidies)){
    # make sure we have all parental ploidies so we can get likelihoods 
    ### change this? ploidy by individual in future?
    stop("All parent ploidies must be in the possible ploidies for the object")
  }
  if(is.null(object$alleleFreq) || attr(object,"alleleFreqType") != "mapping"){
    cat("Allele frequencies for mapping population not found.  Estimating.", 
        sep = "\n")
    allelesin <- max(sapply(donorParentPloidies, sum)) + 
      max(sapply(recurrentParentPloidies, sum))
    possfreq <- seq(0, 1, length.out = (n.gen.backcrossing + 1) * allelesin + 1)
    alldev <- (possfreq[2] - possfreq[1])/2
    object <- AddAlleleFreqMapping(object, allowedDeviation = alldev, 
                                   expectedFreq = possfreq)
  }
  if(is.null(object$genotypeLikelihood)){
    cat("Genotype likelihoods not found.  Estimating.", sep = "\n")
    object <- AddGenotypeLikelihood(object)
  }
  # get most likely genotype for the parents
  pldtot <- sapply(object$genotypeLikelihood, function(x) dim(x)[1] - 1)
  pldtot.don <- pldtot[pldtot %in% sapply(donorParentPloidies, sum)]
  pldtot.rec <- pldtot[pldtot %in% sapply(recurrentParentPloidies, sum)]
  nAlleles <- dim(object$alleleDepth)[2]
  likelyGen.don <- GetLikelyGen(object, donorParent)[as.character(pldtot.don),]
  likelyGen.rec <- GetLikelyGen(object, recurrentParent)[as.character(pldtot.rec),]
  # combinations of parent ploidies that match list of progeny ploidies
  pldcombos <- matrix(NA, nrow = 0, ncol = 2, 
                      dimnames = list(NULL,c("donor","recurrent")))
  for(pl.d in pldtot.don){
    for(pl.r in pldtot.rec){
      if((pl.d/2 + pl.r/2) %in% pldtot && (n.gen.backcrossing == 0 || pl.d == pl.r)){
        pldcombos <- rbind(pldcombos, matrix(c(pl.d, pl.r), nrow = 1, ncol = 2))
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
  
  # expand ploidy combinations across allo and auto types
  ploidies.don <- sapply(donorParentPloidies, sum)
  ploidies.rec <- sapply(recurrentParentPloidies, sum)
  pldcombosExpand <- array(list(), dim = c(0,2), 
                           dimnames = list(NULL, c("donor","recurrent")))
  for(i in 1:dim(pldcombos)[1]){
    thesepl.don <- donorParentPloidies[ploidies.don == pldcombos[i,"donor"]]
    thesepl.rec <- recurrentParentPloidies[ploidies.rec == pldcombos[i,"recurrent"]]
    for(pl.d in thesepl.don){
      for(pl.r in thesepl.rec){
        pldcombosExpand <- rbind(pldcombosExpand, array(list(pl.d, pl.r), dim = c(1,2)))
      }
    }
  }
  # function to generate all gamete genotypes for a set of genotypes.
  # alCopy is a vector of values ranging from zero to ploidy indicating 
  # allele copy number.
  # ploidy is the ploidy
  # rnd indicates which round of the recursive algorithm we are on
  # Output is a matrix.  Alleles are in columns, which should be treated
  # independently.  Rows indicate how many copies of the allele that gamete
  # has.
  makeGametes <- function(alCopy, ploidy, rnd = ploidy[1]/2){
    if(rnd %% 1 != 0 || ploidy < 2){
      stop("Even numbered ploidy needed to simulate gametes.")
    }
    if(length(unique(ploidy)) != 1){
      stop("Currently all subgenome ploidies must be equal.")
    }
    if(any(alCopy > ploidy)){
      stop("Cannot have alCopy greater than ploidy.")
    }
    if(length(ploidy) > 1){ # allopolyploids
      for(pl in ploidy){
        # get allele copies for this isolocus.  Currently simplified; 
        # minimizes the numbe of isoloci to which an allele can belong.
        thisCopy <- alCopy
        thisCopy[thisCopy > pl] <- pl
        alCopy <- alCopy - thisCopy
        # make gametes for this isolocus (recursive, goes to autopoly version)
        thisIsoGametes <- makeGametes(thisCopy, pl, rnd)
      }
      ### continue from here; need to look at all combinations across isoloci ###
      ### maybe use expand.grid? ###
    } else { # autopolyploids
      thisAl <- sapply(alCopy, function(x) c(rep(0L, ploidy - x), rep(1L, x)))
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
  # get prior genotype probabilities for F1
  # backcross
  if(n.gen.backcrossing == 0){
    bcloop <- integer(0)
  } else {
    bcloop <- 1:n.gen.backcrossing
  }
  for(gen in bcloop){
    # reestimate prior probs
  }
  # self
  
  return(object)
}

#### Accessors ####
GetTaxa <- function(object, ...){
  UseMethod("GetTaxa", object)
}
GetTaxa.RADdata <- function(object, ...){
  return(attr(object, "taxa"))
}
GetLoci <- function(object, ...){
  UseMethod("GetLoci", object)
}
GetLoci.RADdata <- function(object, ...){
  return(row.names(object$locTable))
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

# Functions for assigning taxa to specific roles
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
