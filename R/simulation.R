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
simAlleleDepth <- function(locDepth, genotypes, alleles2loc, overdispersion = 20){
  nsam <- nrow(genotypes)
  if(nrow(locDepth) != nsam) stop("genotypes and locDepth should have the same number of rows")
  if(length(alleles2loc) !=  ncol(genotypes)) stop("length of alleles2loc should match number of columns in genotypes.")
  
  alleleDepth <- matrix(0L, nrow = nsam, ncol = ncol(genotypes),
                        dimnames = dimnames(genotypes))
  for(L in unique(alleles2loc)){
    thesecol <- which(alleles2loc == L)
    thisLocDepth <- locDepth[, as.character(L)]
    for(s in 1:nsam){
      alleleDepth[s, thesecol] <-
        sampleReads(genotypes[s, thesecol], thisLocDepth[s], overdispersion)
    }
  }
  
  return(alleleDepth)
}

# data(exampleRAD)
# exampleRAD <- IterateHWE(exampleRAD)
# mygeno <- GetProbableGenotypes(exampleRAD, omit1allelePerLocus = FALSE)[[1]]
# 
# testdepth <- simAlleleDepth(exampleRAD$locDepth, mygeno, exampleRAD$alleles2loc)

# Wrapper function to simulate genotype matrix
simGenotypes <- function(alleleFreq, alleles2loc, nsam, inbreeding, ploidy){
  if(length(alleleFreq) != length(alleles2loc)) stop("Need same number of alleles in alleleFreq and alleles2loc.")
  
  geno <- matrix(0L, nrow = nsam, ncol = length(alleleFreq),
                 dimnames = list(NULL, names(alleleFreq)))
  
  for(L in unique(alleles2loc)){
    thesecol <- which(alleles2loc == L)
    thesefreq <- alleleFreq[thesecol]
    for(s in 1:nsam){
      geno[s, thesecol] <- sampleGenotype(thesefreq, inbreeding, ploidy)
    }
  }
  
  return(geno)
}

# testgeno <- simGenotypes(exampleRAD$alleleFreq, exampleRAD$alleles2loc, 100, 0.2, 2)
# mean(testgeno[,1] == 1) # HO of 0.26
# 1 - sum(exampleRAD$alleleFreq[1:2]^2) # HE of 0.316
# 1 - (0.26/0.316) # 0.18, about right with sampling error

# Get expected Hind/He distribution based on depths and allele freqs in a RADdata object
expectedHindHe <- function(object, ploidy = object$possiblePloidies[[1]],
                           inbreeding = 0, overdispersion = 20, reps = 10,
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
    if(!quiet) message(paste("Simulating rep", i))
    geno <- simGenotypes(object$alleleFreq, object$alleles2loc, nTaxa(object),
                         inbreeding, ploidy)
    depths <- simAlleleDepth(object$locDepth, geno, object$alleles2loc,
                             overdispersion)
    rownames(depths) <- GetTaxa(object)
    simrad <- RADdata(depths, object$alleles2loc, object$locTable,
                      object$possiblePloidies, GetContamRate(object),
                      object$alleleNucleotides)
    out[,i] <- colMeans(HindHe(simrad), na.rm = TRUE)
  }
  
  if(!quiet){
    q <- quantile(out, probs = c(0.025, 0.975))
    cat(c(paste("Mean Hind/He:", formatC(mean(out), digits = 3)),
          paste("Standard deviation:", formatC(sd(out), digits = 3)),
          paste("95% of observations between", formatC(q[1], digits = 3),
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
# expectedHindHe(mydata) # get from VCF in tutorial
