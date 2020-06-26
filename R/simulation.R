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

# Randomly generate reads based on a genotype and overdispersion
sampleReads <- function(geno, nreads, overdispersion = 20){
  initprobs <- geno / sum(geno)
  alpha <- initprobs * overdispersion
  newprobs <- rgamma(length(geno), alpha, 1)
  
  reads <- sample.int(length(geno), size = nreads, replace = TRUE, prob = newprobs)
  out <- sapply(1:length(geno), function(x) sum(reads == x))
  return(out)
}

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
