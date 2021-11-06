# testing
# sampleReads(c(1, 2, 1), 30, overdispersion = 100)
# 
# testprob <- c(0.5, 0.25, 0.25)
# 
# test1 <- t(sapply(1:1000, function(x) sampleReads(testprob * 4, 30, overdispersion = 20)))
# test2 <- MultiRNG::draw.multinomial(1000, 3, testprob, 30)
# 
# median(apply(test1, 1, dmultinom, prob = testprob))
# median(apply(test2, 1, dmultinom, prob = testprob))

# Wrapper function to simulate an allele depth matrix, given locus depth and genotypes
SimAlleleDepth <- function(locDepth, genotypes, alleles2loc, overdispersion = 20){
  nsam <- nrow(genotypes)
  if(nrow(locDepth) != nsam) stop("genotypes and locDepth should have the same number of rows")
  if(length(alleles2loc) !=  ncol(genotypes)) stop("length of alleles2loc should match number of columns in genotypes.")
  
  locnames <- as.character(seq_len(max(alleles2loc)))
  if(!all(locnames %in% colnames(locDepth))){
    stop("Loci missing from locDepth or named incorrectly.")
  }
  locDepth <- locDepth[,locnames] # ensure correct order for Rcpp fn
  
  alleleDepth <- simAD(locDepth, genotypes, alleles2loc, overdispersion) # Rcpp fn
  dimnames(alleleDepth) <- dimnames(genotypes)
  
  return(alleleDepth)
}

# data(exampleRAD)
# exampleRAD <- IterateHWE(exampleRAD)
# mygeno <- GetProbableGenotypes(exampleRAD, omit1allelePerLocus = FALSE)[[1]]
# 
# testdepth <- simAlleleDepth(exampleRAD$locDepth, mygeno, exampleRAD$alleles2loc)

# Wrapper function to simulate genotype matrix
SimGenotypes <- function(alleleFreq, alleles2loc, nsam, inbreeding, ploidy){
  if(length(alleleFreq) != length(alleles2loc)) stop("Need same number of alleles in alleleFreq and alleles2loc.")
  
  geno <- simGeno(alleleFreq, alleles2loc, nsam, inbreeding, ploidy) # Rcpp fn
  colnames(geno) <- names(alleleFreq)
  
  return(geno)
}

# testgeno <- simGenotypes(exampleRAD$alleleFreq, exampleRAD$alleles2loc, 100, 0.2, 2)
# mean(testgeno[,1] == 1) # HO of 0.26
# 1 - sum(exampleRAD$alleleFreq[1:2]^2) # HE of 0.316
# 1 - (0.26/0.316) # 0.18, about right with sampling error

# Get expected Hind/He distribution based on depths and allele freqs in a RADdata object
ExpectedHindHe <- function(object, ploidy = object$possiblePloidies[[1]],
                           inbreeding = 0, overdispersion = 20,
                           reps = ceiling(5000 / nLoci(object)),
                           quiet = FALSE, plot = TRUE){
  if(length(ploidy) != 1){
    stop("Please give a single value for ploidy; function assumes diploid or autopolyploid inheritance")
  }
  if(is.null(object$alleleFreq)){
    object <- AddAlleleFreqHWE(object)
  }
  ### Add something here to filter based on allele frequency or preliminary Hind/He
  
  # matrix of hind/he values for simulated loci
  out <- matrix(0, nrow = nLoci(object), ncol = reps,
                dimnames = list(GetLoci(object), NULL))
  
  for(i in seq_len(reps)){
    if(!quiet && i %% 10 == 1) message(paste("Simulating rep", i))
    geno <- matrix(NA_integer_, nrow = nTaxa(object), ncol = nAlleles(object),
                   dimnames = list(GetTaxa(object), GetAlleleNames(object)))
    for(h in unique(GetTaxaPloidy(object))){
      if((ploidy * h) %% 2L != 0){
        stop("Either ploidy or all taxa ploidies must be even")
      }
      thesesam <- which(GetTaxaPloidy(object) == h)
      geno[thesesam,] <- SimGenotypes(object$alleleFreq, object$alleles2loc,
                                      length(thesesam),
                                      inbreeding, (ploidy * h) %/% 2L)
    }
    depths <- SimAlleleDepth(object$locDepth, geno, object$alleles2loc,
                             overdispersion)
    rownames(depths) <- GetTaxa(object)
    simrad <- RADdata(depths, object$alleles2loc, object$locTable,
                      object$possiblePloidies, GetContamRate(object),
                      object$alleleNucleotides, GetTaxaPloidy(object))
    out[,i] <- colMeans(HindHe(simrad), na.rm = TRUE)
  }
  
  if(!quiet){
    message(paste("Completed", reps, "simulation reps."))
    q <- quantile(out, probs = c(0.025, 0.975), na.rm = TRUE)
    cat(c(paste("Mean Hind/He:", formatC(mean(out, na.rm = TRUE), digits = 3)),
          paste("Standard deviation:", formatC(sd(out, na.rm = TRUE), digits = 3)),
          paste("95% of observations are between", formatC(q[1], digits = 3),
                "and", formatC(q[2], digits = 3))),
        sep = "\n")
  }
  
  if(plot){
    hist(out, xlab = "Hind/He", main = "Expected distribution of Hind/He",
         breaks = 30)
  }
  
  invisible(out)
}

# data(exampleRAD)
# expectedHindHe(exampleRAD)
# expectedHindHe(mydata, reps = 10) # get from VCF in tutorial

# Simulated genotypes in a mapping population
SimGenotypesMapping <- function(donorGen, recurGen, alleles2loc, nsam, ploidy,
                                n.gen.backcrossing, n.gen.selfing){
  if(length(donorGen) != length(recurGen) ||
     length(donorGen) != length(alleles2loc)){
    stop("Need same number of alleles in donorGen, recurGen, and alleles2loc.")
  }
  # get possible progeny and their probabilities, treating each allele
  # from each parent as unique.  (multiallelic genotypes)
  progprob <- .buildProgProb(ploidy, n.gen.backcrossing, n.gen.selfing)
  
  geno <- simGenoMapping(donorGen, recurGen, progprob[[1]], progprob[[2]],
                         alleles2loc, nsam, ploidy) # Rcpp fn
  colnames(geno) <- names(donorGen)
  
  return(geno)
}

ExpectedHindHeMapping <- function(object, ploidy = object$possiblePloidies[[1]],
                           n.gen.backcrossing = 0, n.gen.selfing = 0, overdispersion = 20,
                           freqAllowedDeviation = 0.05, minLikelihoodRatio = 10,
                           reps = ceiling(5000 / nLoci(object)),
                           quiet = FALSE, plot = TRUE){
  if(length(ploidy) != 1){
    stop("Please give a single value for ploidy; function assumes diploid or autopolyploid inheritance")
  }
  possfreq <- seq(0, 1, length.out = (n.gen.backcrossing + 1) * ploidy * 2 + 1)
  if(is.null(object$likelyGeno_donor) || is.null(object$likelyGeno_recurrent)){
    object <- AddAlleleFreqMapping(object, expectedFreqs = possfreq,
                                   allowedDeviation = freqAllowedDeviation)
    object <- AddGenotypeLikelihood(object, overdispersion = overdispersion)
    object <- EstimateParentalGenotypes(object, n.gen.backcrossing = n.gen.backcrossing,
                                        n.gen.intermating = 0, n.gen.selfing = n.gen.selfing,
                                        donorParentPloidies = list(ploidy),
                                        recurrentParentPloidies = list(ploidy),
                                        minLikelihoodRatio = minLikelihoodRatio)
  }
  
  # matrix of hind/he values for simulated loci
  out <- matrix(0, nrow = nLoci(object), ncol = reps,
                dimnames = list(GetLoci(object), NULL))
  donParent <- GetDonorParent(object)
  recParent <- GetRecurrentParent(object)
  progeny <- setdiff(GetTaxa(object), c(donParent, recParent, GetBlankTaxa(object)))
  
  for(i in seq_len(reps)){
    if(!quiet && i %% 10 == 1) message(paste("Simulating rep", i))
    geno <- SimGenotypesMapping(object$likelyGeno_donor[as.character(ploidy),],
                                object$likelyGeno_recurrent[as.character(ploidy),],
                                object$alleles2loc, length(progeny), ploidy,
                                n.gen.backcrossing, n.gen.selfing)
    geno <- rbind(object$likelyGeno_donor[as.character(ploidy),],
                  object$likelyGeno_recurrent[as.character(ploidy),],
                  geno)
    depths <- SimAlleleDepth(object$locDepth[c(donParent, recParent, progeny),],
                             geno, object$alleles2loc,
                             overdispersion)
    rownames(depths) <- c(donParent, recParent, progeny)
    simrad <- RADdata(depths, object$alleles2loc, object$locTable,
                      object$possiblePloidies, GetContamRate(object),
                      object$alleleNucleotides)
    simrad <- SetDonorParent(simrad, donParent)
    simrad <- SetRecurrentParent(simrad, recParent)
    simrad <- AddAlleleFreqMapping(simrad, expectedFreqs = possfreq,
                                   allowedDeviation = freqAllowedDeviation)
    simrad <- AddDepthSamplingPermutations(simrad)
    simrad <- AddGenotypeLikelihood(simrad, overdispersion = overdispersion)
    simrad <- EstimateParentalGenotypes(simrad, n.gen.backcrossing = n.gen.backcrossing,
                                        n.gen.intermating = 0, n.gen.selfing = n.gen.selfing,
                                        donorParentPloidies = list(ploidy),
                                        recurrentParentPloidies = list(ploidy),
                                        minLikelihoodRatio = minLikelihoodRatio)
    hh <- HindHeMapping(simrad, n.gen.backcrossing, 0, n.gen.selfing, ploidy,
                        minLikelihoodRatio)
    out[colnames(hh),i] <- colMeans(hh, na.rm = TRUE)
  }
  
  if(!quiet){
    message(paste("Completed", reps, "simulation reps."))
    q <- quantile(out, probs = c(0.025, 0.975), na.rm = TRUE)
    cat(c(paste("Mean Hind/He:", formatC(mean(out, na.rm = TRUE), digits = 3)),
          paste("Standard deviation:", formatC(sd(out, na.rm = TRUE), digits = 3)),
          paste("95% of observations are between", formatC(q[1], digits = 3),
                "and", formatC(q[2], digits = 3))),
        sep = "\n")
  }
  
  if(plot){
    hist(out, xlab = "Hind/He", main = "Expected distribution of Hind/He",
         breaks = 30)
  }
  
  invisible(out)
}
