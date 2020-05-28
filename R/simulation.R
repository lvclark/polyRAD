# Randomly generate a genotype based on allele frequencies and inbreeding
sampleGenotype <- function(freq, inbreeding, ploidy){
  nal <- length(freq)
  geno <- integer(nal)
  for(k in 1:ploidy){
    # get probability that this allele is the same as a previous one
    repprob <- 1 - (1 - inbreeding) ^ (k - 1)
    if(sample(c(TRUE, FALSE), size = 1, prob = c(repprob, 1 - repprob))){
      a <- sample(1:nal, size = 1, prob = geno)
    } else {
      a <- sample(1:nal, size = 1, prob = freq)
    }
    geno[a] <- geno[a] + 1L
  }
  
  return(geno)
}
