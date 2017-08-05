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
                        depthRatio = depthRatio, antiAlleleDepth = antiAlleleDepth), 
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
                                         deleteLociOutsideFreqRange = FALSE){
  if(min(dist(expectedFreqs, method = "manhattan"))/2 < allowedDeviation){
    stop("allowedDeviation is too large given intervals within expectedFreqs")
  }
  meanRat <- colMeans(object$depthRatio, na.rm = TRUE)
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
  if(deleteLociOutsideFreqRange){
    lociToDelete <- unique(object$allele2loc[is.na(outFreq)])
    object <- SubsetLoci(object, loci = lociToDelete, delete = TRUE)
  }
  
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

#### Accessors ####
GetTaxa <- function(object, ...){
  UseMethod("GetTaxa", object)
}
GetTaxa.RADdata <- function(object, ...){
  return(attr(object, "taxa"))
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