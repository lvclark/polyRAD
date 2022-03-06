# additional polyRAD functions that perform calculations.
# various internal functions.

# internal function to take allele frequencies and get prior probs under HWE
# freqs is a vector of allele frequencies
# ploidy is a vector indicating the ploidy
.HWEpriors <- function(freqs, ploidy, selfing.rate){
  if(length(unique(ploidy)) != 1){
    stop("All subgenomes must be same ploidy")
  }
  if(selfing.rate < 0 || selfing.rate > 1){
    stop("selfing.rate must not be less than zero or more than one.")
  }
  nsubgen <- length(ploidy)
  if(nsubgen == 1){ # for diploid/autopolyploid, or single subgenome with recursion
    priors <- matrix(NA, nrow = ploidy+1, ncol = length(freqs),
                     dimnames = list(as.character(0:ploidy), names(freqs)))
    antifreqs <- 1 - freqs
    # genotype probabilities under random mating
    for(i in 0:ploidy){
      priors[i+1,] <- choose(ploidy, i) * freqs ^ i * antifreqs ^ (ploidy - i)
    }
    if(selfing.rate > 0 && selfing.rate < 1 && ploidy %% 2 != 0){
      warning("Not simulating self-fertilization for individuals of odd ploidy.")
    }
    # adjust for self fertilization if applicable
    if(selfing.rate > 0 && selfing.rate < 1 && ploidy %% 2 == 0){
      sm <- .selfmat(ploidy)
      # Equation 6 from de Silva et al. 2005 (doi:10.1038/sj.hdy.6800728)
      priors <- (1 - selfing.rate) * 
        solve(diag(ploidy + 1) - selfing.rate * sm, priors)
      rownames(priors) <- as.character(0:ploidy)
    }
    if(selfing.rate == 1){
      priors <- matrix(0, nrow = ploidy+1, ncol = length(freqs),
                       dimnames = list(as.character(0:ploidy), names(freqs)))
      priors[1,] <- antifreqs
      priors[ploidy + 1, ] <- freqs
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
      thesepriors <- .HWEpriors(nsubgen * thesefreqs, pld, selfing.rate)
      
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
  
  ploidytotpriors <- sapply(object$priorProb[,1], function(x) dim(x)[1] - 1)
  ploidytotlikeli <- sapply(object$genotypeLikelihood[,1], function(x) dim(x)[1] - 1)
  
  results <- array(list(),
                   dim = dim(object$priorProb),
                   dimnames = dimnames(object$priorProb))
  
  for(i in seq_len(nrow(object$priorProb))){
    j <- which(ploidytotlikeli == ploidytotpriors[i])
    stopifnot(length(j) == 1)
    for(h in seq_len(ncol(object$priorProb))){
      thesetaxa <- dimnames(object$genotypeLikelihood[[j,h]])[[2]]
      if(attr(object, "priorType") == "population"){
        # expand priors out by individuals
        thispriorarr <- array(object$priorProb[[i,h]], 
                              dim = c(dim(object$priorProb[[i,h]])[1], 1, 
                                      dim(object$priorProb[[i,h]])[2]))[,rep(1, length(thesetaxa)),, drop = FALSE]
        dimnames(thispriorarr) <- dimnames(object$genotypeLikelihoods[[j,h]])
      } else {
        thispriorarr <- object$priorProb[[i,h]]
      }
      stopifnot(identical(dim(thispriorarr), dim(object$genotypeLikelihood[[j,h]])))
      results[[i,h]] <- thispriorarr * object$genotypeLikelihood[[j,h]]
      # factor in LD if present
      if(!is.null(object$priorProbLD)){
        results[[i,h]] <- results[[i,h]] * object$priorProbLD[[i,h]]
      }
      # find any that total to zero (within taxon x allele) and replace with priors
      totzero <- which(colSums(results[[i,h]]) == 0)
      if(length(totzero) > 0){
        for(a in 1:dim(thispriorarr)[1]){
          results[[i,h]][a,,][totzero] <- thispriorarr[a,,][totzero]
        }
      }
      # in a mapping population, don't use priors for parents
      if(!is.null(attr(object, "donorParent")) &&
         !is.null(attr(object, "recurrentParent"))){
        parents <- c(GetDonorParent(object), GetRecurrentParent(object))
        parents <- intersect(parents, thesetaxa)
        if(length(parents) > 0){ # only fix if this is the ploidy for at least one parent
          results[[i,h]][, parents, ] <- object$genotypeLikelihood[[j,h]][, parents, ]
        }
      }
    }
  }
  
  return(results)
}

# internal function to get best estimate of allele frequencies depending
# on what parameters are available
.alleleFreq <- function(object, type = "choose", taxaToKeep = GetTaxa(object)){
  if(!"RADdata" %in% class(object)){
    stop("RADdata object needed for .priorTimesLikelihood")
  }
  if(!type %in% c("choose", "individual frequency", "posterior prob",
                  "depth ratio")){
    stop("Type must be 'choose', 'individual frequency', 'posterior prob', or 'depth ratio'.")
  }
  if(type == "individual frequency" && is.null(object$alleleFreqByTaxa)){
    stop("Need alleleFreqByTaxa if type = 'individual frequency'.")
  }
  if(type == "posterior prob" && 
     (is.null(object$posteriorProb) || is.null(object$ploidyChiSq))){
    stop("Need posteriorProb and ploidyChiSq if type = 'posterior prob'.")
  }
  
  if(type %in% c("choose", "individual frequency") &&
     !is.null(object$alleleFreqByInd)){
    outFreq <- colMeans(object$alleleFreqByTaxa[taxaToKeep,,drop = FALSE], 
                        na.rm = TRUE)
    attr(outFreq, "type") <- "individual frequency"
  } else {
    if(type %in% c("choose", "posterior prob") &&
       CanDoGetWeightedMeanGeno(object)){
      wmgeno <- GetWeightedMeanGenotypes(object, minval = 0, maxval = 1,
                                         omit1allelePerLocus = FALSE)
      outFreq <- colMeans(wmgeno[taxaToKeep,,drop = FALSE], na.rm = TRUE)
      attr(outFreq, "type") <- "posterior prob"
    } else {
      if(type %in% c("choose", "depth ratio")){
        outFreq <- colMeans(object$depthRatio[taxaToKeep,,drop = FALSE],
                            na.rm = TRUE)
        attr(outFreq, "type") <- "depth ratio"
      }
    }
  }
  return(outFreq)
}

# Internal function to get a consensus between two DNA sequences, using
# IUPAC ambiguity codes.  Should work for either character strings
# or DNAStrings.  Puts in ambiguity codes any time there is not 100%
# consensus.
.mergeNucleotides <- function(seq1, seq2){
  split1 <- strsplit(seq1, split = "")[[1]]
  split2 <- strsplit(seq2, split = "")[[1]]
  splitNew <- character(length(split1))
  splitNew[split1 == split2] <- split1[split1 == split2]
  
  IUPAC_key <- list(M = list(c('A', 'C'), c('A', 'M'), c('C', 'M')), 
                    R = list(c('A', 'G'), c('A', 'R'), c('G', 'R')), 
                    W = list(c('A', 'T'), c('A', 'W'), c('T', 'W')), 
                    S = list(c('C', 'G'), c('C', 'S'), c('G', 'S')),
                    Y = list(c('C', 'T'), c('C', 'Y'), c('T', 'Y')),
                    K = list(c('G', 'T'), c('G', 'K'), c('T', 'K')),
                    V = list(c('A', 'S'), c('C', 'R'), c('G', 'M'),
                             c('A', 'V'), c('C', 'V'), c('G', 'V')),
                    H = list(c('A', 'Y'), c('C', 'W'), c('T', 'M'),
                             c('A', 'H'), c('C', 'H'), c('T', 'H')),
                    D = list(c('A', 'K'), c('G', 'W'), c('T', 'R'),
                             c('A', 'D'), c('G', 'D'), c('T', 'D')),
                    B = list(c('C', 'K'), c('G', 'Y'), c('T', 'S'),
                             c('C', 'B'), c('G', 'B'), c('T', 'B')),
                    # throw out deletions with respect to reference
                    A = list(c('A', '-')), C = list(c('C', '-')),
                    G = list(c('G', '-')), T = list(c('T', '-')),
                    # throw out insertions with respect to reference
                    . = list(c('A', '.'), c('C', '.'), c('G', '.'),
                             c('T', '.'), c('M', '.'), c('R', '.'),
                             c('W', '.'), c('S', '.'), c('Y', '.'),
                             c('K', '.'), c('V', '.'), c('H', '.'),
                             c('D', '.'), c('B', '.'), c('N', '.')))
  for(i in which(split1 != split2)){
    nucset <- c(split1[i], split2[i])
    splitNew[i] <- "N"
    for(nt in names(IUPAC_key)){
      for(matchset in IUPAC_key[[nt]]){
        if(setequal(nucset, matchset)){
          splitNew[i] <- nt
          break
        }
      }
    }
  }
  out <- paste(splitNew, collapse = "")

  return(out)
}

# substitution matrix for distances between alleles in polyRAD.
# Any partial match based on ambiguity counts as a complete match.
polyRADsubmat <- matrix(c(0,1,1,1, 0,0,0,1,1,1, 0,0,0,1, 0,1,1, # A
                          1,0,1,1, 0,1,1,0,0,1, 0,0,1,0, 0,1,1, # C
                          1,1,0,1, 1,0,1,0,1,0, 0,1,0,0, 0,1,1, # G
                          1,1,1,0, 1,1,0,1,0,0, 1,0,0,0, 0,1,1, # T
                          0,0,1,1, 0,0,0,0,0,1, 0,0,0,0, 0,1,1, # M
                          0,1,0,1, 0,0,0,0,1,0, 0,0,0,0, 0,1,1, # R
                          0,1,1,0, 0,0,0,1,0,0, 0,0,0,0, 0,1,1, # W
                          1,0,0,1, 0,0,1,0,0,0, 0,0,0,0, 0,1,1, # S
                          1,0,1,0, 0,1,0,0,0,0, 0,0,0,0, 0,1,1, # Y
                          1,1,0,0, 1,0,0,0,0,0, 0,0,0,0, 0,1,1, # K
                          0,0,0,1, 0,0,0,0,0,0, 0,0,0,0, 0,1,1, # V
                          0,0,1,0, 0,0,0,0,0,0, 0,0,0,0, 0,1,1, # H
                          0,1,0,0, 0,0,0,0,0,0, 0,0,0,0, 0,1,1, # D
                          1,0,0,0, 0,0,0,0,0,0, 0,0,0,0, 0,1,1, # B
                          0,0,0,0, 0,0,0,0,0,0, 0,0,0,0, 0,0,0, # N
                          1,1,1,1, 1,1,1,1,1,1, 1,1,1,1, 0,0,1, # -
                          1,1,1,1, 1,1,1,1,1,1, 1,1,1,1, 0,1,0  # .
                          ),
                         nrow = 17, ncol = 17,
                         dimnames = list(c('A', 'C', 'G', 'T', 'M', 'R', 'W', 
                                           'S', 'Y', 'K', 'V', 'H', 'D', 'B', 
                                           'N', '-', '.'),
                                         c('A', 'C', 'G', 'T', 'M', 'R', 'W', 
                                           'S', 'Y', 'K', 'V', 'H', 'D', 'B', 
                                           'N', '-', '.')))

# Reverse complement; set up as S4 method to be compatible with Biostrings
setGeneric("reverseComplement", signature="x",
           function(x, ...) standardGeneric("reverseComplement")
)
setMethod("reverseComplement", "character",
          function(x, ...){
            stri_reverse(stri_trans_char(x,
                                         "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
                                         "TGCAYRSWMKVHDBNtgcayrswmkvhdbn"))
          } 
)

#### Getting genotype priors in mapping pops and simulating selfing ####

# function to generate all gamete genotypes for a set of genotypes.
# alCopy is a vector of values ranging from zero to ploidy indicating 
# allele copy number.
# ploidy is the ploidy
# rnd indicates which round of the recursive algorithm we are on
# Output is a matrix.  Alleles are in columns, which should be treated
# independently.  Rows indicate gametes, with values indicating how many
# copies of the allele that gamete has.
.makeGametes <- function(alCopy, ploidy, rnd = ploidy[1]/2){
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
      thisIsoGametes <- .makeGametes(thisCopy, pl, rnd)
      # add to current gamete set
      nGamCurr <- dim(thisAl)[1]
      nGamNew <- dim(thisIsoGametes)[1]
      thisAl <- thisAl[rep(1:nGamCurr, each = nGamNew),, drop = FALSE] + 
        thisIsoGametes[rep(1:nGamNew, times = nGamCurr),, drop = FALSE]
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
      thisAl <- thisAl[rep(1:ploidy, each = nReps),, drop = FALSE]
      for(i in 1:ploidy){
        thisAl[((i-1)*nReps+1):(i*nReps),] <- 
          thisAl[((i-1)*nReps+1):(i*nReps),, drop=FALSE] + 
          .makeGametes(alCopy - thisAl[i*nReps,, drop = FALSE], ploidy - 1, rnd - 1)
      }
    }
  }
  return(thisAl)
}
# function to get probability of a gamete with a given allele copy number,
# given output from makeGametes
.gameteProb <- function(makeGamOutput, ploidy){
  outmat <- sapply(0:(sum(ploidy)/2), 
                   function(x) colMeans(makeGamOutput == x))
  if(ncol(makeGamOutput) == 1){
    outmat <- matrix(outmat, nrow = ploidy/2 + 1, ncol = 1)
  } else {
    outmat <- t(outmat)
  }
  
  return(outmat)
}
# function to take two sets of gamete probabilities (from two parents)
# and output genotype probabilities
.progenyProb <- function(gamProb1, gamProb2){
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
.gameteProbPop <- function(priors, ploidy){
  # get gamete probs for all possible genotypes
  possGamProb <- .gameteProb(.makeGametes(1:dim(priors)[1] - 1, ploidy),ploidy)
  # output matrix
  outmat <- possGamProb %*% priors
  return(outmat)
}
# function just to make selfing matrix.
# Every parent genotype is a column.  The frequency of its offspring after
# selfing is in rows.
.selfmat <- function(ploidy){
  # get gamete probs for all possible genotypes
  possGamProb <- .gameteProb(.makeGametes(0:ploidy, ploidy), ploidy)
  # genotype probabilities after selfing each possible genotype
  outmat <- matrix(0, nrow = ploidy + 1, ncol = ploidy + 1)
  for(i in 1:(ploidy + 1)){
    outmat[,i] <- .progenyProb(possGamProb[,i,drop = FALSE],
                               possGamProb[,i,drop = FALSE])
  }
  return(outmat)
}
# function to adjust genotype probabilities from one generation of selfing
.selfPop <- function(priors, ploidy){
  # progeny probs for all possible genotypes, selfed
  possProgenyProb <- .selfmat(ploidy)
  # multiple progeny probs by prior probabilities of those genotypes
  outmat <- possProgenyProb %*% priors
  return(outmat)
}

# Functions used by HindHeMapping ####

# Function for making multiallelic gametes.  Autopolyploid segretation only.
# geno is a vector of length ploidy with values indicating allele identity.
# rnd is how many alleles will be in the gamete.
.makeGametes2 <- function(geno, rnd = length(geno)/2){
  nal <- length(geno) # number of remaining alleles in genotype
  if(rnd == 1){
    return(matrix(geno, nrow = nal, ncol = 1))
  } else {
    outmat <- matrix(0L, nrow = choose(nal, rnd), ncol = rnd)
    currrow <- 1
    for(i in 1:(nal - rnd + 1)){
      submat <- .makeGametes2(geno[-(1:i)], rnd - 1)
      theserows <- currrow + 1:nrow(submat) - 1
      currrow <- currrow + nrow(submat)
      outmat[theserows,1] <- geno[i]
      outmat[theserows,2:ncol(outmat)] <- submat
    }
    return(outmat)
  }
}
# take output of .makeGametes2 and get progeny probabilities.
# Return a list, where the first item is a matrix with progeny genotypes in
# rows, and the second item is a vector of corresponding probabilities.
.progenyProb2 <- function(gam1, gam2){
  allprog <- cbind(gam1[rep(1:nrow(gam1), each = nrow(gam2)),],
                   gam2[rep(1:nrow(gam2), times = nrow(gam1)),])
  allprog <- t(apply(allprog, 1, sort)) # could be sped up with Rcpp
  # set up output
  outprog <- matrix(allprog[1,], nrow = 1, ncol = ncol(allprog))
  outprob <- 1/nrow(allprog)
  # tally up unique genotypes
  for(i in 2:nrow(allprog)){
    thisgen <- allprog[i,]
    existsYet <- FALSE
    for(j in 1:nrow(outprog)){
      if(identical(thisgen, outprog[j,])){
        outprob[j] <- outprob[j] + 1/nrow(allprog)
        existsYet <- TRUE
        break
      }
    }
    if(!existsYet){
      outprog <- rbind(outprog, thisgen, deparse.level = 0)
      outprob <- c(outprob, 1/nrow(allprog))
    }
  }
  return(list(outprog, outprob))
}
# Consolidate a list of outputs from .progenyProb2 into a single output.
# It is assumed the all probabilities in pplist have been multiplied by
# some factor so they sum to one.
.consolProgProb <- function(pplist){
  outprog <- pplist[[1]][[1]]
  outprob <- pplist[[1]][[2]]
  if(length(pplist) > 1){
    for(i in 2:length(pplist)){
      for(pr1 in 1:nrow(pplist[[i]][[1]])){
        thisgen <- pplist[[i]][[1]][pr1,]
        existsYet <- FALSE
        for(pr2 in 1:nrow(outprog)){
          if(identical(thisgen, outprog[pr2,])){
            outprob[pr2] <- outprob[pr2] + pplist[[i]][[2]][pr1]
            existsYet <- TRUE
            break
          }
        }
        if(!existsYet){
          outprog <- rbind(outprog, thisgen, deparse.level = 0)
          outprob <- c(outprob, pplist[[i]][[2]][pr1])
        }
      }
    }
  }
  
  return(list(outprog, outprob))
}

# Build progeny probability object for given mapping population properties
.buildProgProb <- function(ploidy1, ploidy2, gen_backcrossing, gen_selfing){
  # set up parents; number indexes locus copy
  p1 <- 1:ploidy1
  p2 <- (1:ploidy2) + ploidy1
  # create F1 progeny probabilities
  progprob <- .progenyProb2(.makeGametes2(p1), .makeGametes2(p2))
  # backcross
  if(gen_backcrossing > 0){
    gam1 <- .makeGametes2(p1)
    for(g in 1:gen_backcrossing){
      allprogprobs <- lapply(1:nrow(progprob[[1]]),
                             function(x){
                               pp <- .progenyProb2(gam1, .makeGametes2(progprob[[1]][x,]))
                               pp[[2]] <- pp[[2]] * progprob[[2]][x]
                               return(pp)
                             })
      progprob <- .consolProgProb(allprogprobs)
    }
  }
  # self-fertilize
  if(gen_selfing > 0){
    for(g in 1:gen_selfing){
      allprogprobs <- lapply(1:nrow(progprob[[1]]),
                             function(x){
                               gam <- .makeGametes2(progprob[[1]][x,])
                               pp <- .progenyProb2(gam, gam)
                               pp[[2]] <- pp[[2]] * progprob[[2]][x]
                               return(pp)
                             })
      progprob <- .consolProgProb(allprogprobs)
    }
  }
  return(progprob)
}

# Probability that, depending on the generation, two alleles sampled in a progeny
# are different locus copies from the same parent, or from different parents.
# If there is backcrossing, parent 1 is the recurrent parent.
# (No double reduction.)
# Can take a few seconds to run if there are many generations, but it is only
# intended to be run once for the whole dataset.
.progAlProbs <- function(ploidy1, ploidy2, gen_backcrossing, gen_selfing){
  # set up parents; number indexes locus copy
  p1 <- 1:ploidy1
  p2 <- (1:ploidy2) + ploidy1
  # probabilities of progeny genotypes
  progprob <- .buildProgProb(ploidy1, ploidy2, gen_backcrossing, gen_selfing)
  ploidy.prog <- ncol(progprob[[1]])
  
  # total probability that (without replacement, from individual progeny):
  diffp1 <- 0 # two different locus copies, both from parent 1, are sampled
  diffp2 <- 0 # two different locus copies, both from parent 2, are sampled
  diff12 <- 0 # locus copies from two different parents are sampled
  
  ncombo <- choose(ploidy.prog, 2) # number of ways to choose 2 alleles from a genotype
  # examine each progeny genotype
  for(p in 1:nrow(progprob[[1]])){
    thisgen <- progprob[[1]][p,]
    thisprob <- progprob[[2]][p]
    for(m in 1:(ploidy.prog - 1)){
      al1 <- thisgen[m]
      for(n in (m + 1):ploidy.prog){
        al2 <- thisgen[n]
        if(al1 == al2) next
        if((al1 %in% p1) && (al2 %in% p1)){
          diffp1 <- diffp1 + (thisprob / ncombo)
          next
        }
        if((al1 %in% p2) && (al2 %in% p2)){
          diffp2 <- diffp2 + (thisprob / ncombo)
          next
        }
        if(((al1 %in% p1) && (al2 %in% p2)) ||
           ((al2 %in% p1) && (al1 %in% p2))){
          diff12 <- diff12 + (thisprob / ncombo)
          next
        }
        stop("Allele indexing messed up.")
      }
    }
  }
  
  return(c(diffp1, diffp2, diff12))
}
