#source('D:/GitHub/polyRAD/R/classes_methods.R')

# make a simple dataset to use in documentation examples

alleles2loc <- c(1,1,2,2,2,3,3,4,4)
alNuc <- c("A", "T", "AC", "GC", "GT", "T", "C", "A", "G")
taxa <- paste("sample", formatC(1:100, format='d', width = 3, flag = '0'), sep = "")
loci <- paste("loc", 1:4, sep = "")
realFreq <- c(0.2, 0.8, 0.1, 0.3, 0.6, 0.95, 0.05, 0.3, 0.7)
locPloidies <- c(2,2,4,2)
taxaPloidy <- rep(c(2L, 3L), times = c(80,20))
realGen <- matrix(0L, nrow = 100, ncol = length(alleles2loc))
depth <- matrix(NA_integer_, nrow = 100, ncol = length(alleles2loc),
                dimnames = list(taxa, paste(loci[alleles2loc], alNuc, sep = "_")))
for(L in 1:max(alleles2loc)){
  thesecol <- which(alleles2loc == L)
  thesefreq <- realFreq[thesecol]
  for(s in 1:length(taxa)){
    thispld <- locPloidies[L] * taxaPloidy[s] / 2L
    for(i in 1:thispld){
      thisAl <- sample(thesecol, size = 1, prob = thesefreq)
      realGen[s,thisAl] <- realGen[s,thisAl] + 1
    }
    totcounts <- round(rnorm(1, mean = 20, sd = 10))
    if(totcounts < 0) totcounts <- 0
    reads <- sample(thesecol, size = totcounts, replace = TRUE, prob = realGen[s, thesecol])
    for(i in thesecol){
      depth[s,i] <- sum(reads == i)
    }
  }
}

realGen[1:10,]
tail(realGen)
depth[1:10,]

exampleRAD <- RADdata(depth, as.integer(alleles2loc), 
                      data.frame(row.names = loci, Chr = c(1,4,6,9), 
                                 Pos = c(1111020, 33069, 2637920, 5549872)),
                      list(2L, 4L), 0.001, alNuc, taxaPloidy)
exampleRAD
exampleRAD$depthRatio[1:10,]

save(exampleRAD, file = "data/exampleRAD.RData")

# make a simple mapping dataset to use in documentation examples
alleles2loc <- rep(1:4, each = 2)
alNuc <- c("A", "G", "A", "C", "G", "T", "T", "A")
taxa <- c("parent1", "parent2",
          paste("progeny", formatC(1:100, format='d', width = 3, flag = '0'), sep = ""))
loci <- paste("loc", 1:4, sep = "")
depth <- matrix(NA_integer_, nrow = length(taxa), ncol = length(alleles2loc),
                dimnames = list(taxa, paste(loci[alleles2loc], alNuc, sep = "_")))
# simulate a diploid BC1 pop
for(L in 1:max(alleles2loc)){
  thesecol <- which(alleles2loc == L)
  for(s in 1:length(taxa)){
    totcounts <- round(rnorm(1, mean = 20, sd = 10))
    if(totcounts < 0) totcounts <- 0
    if(s %in% 1:2){ # make sure parents have data
      while(totcounts == 0){
        totcounts <- as.integer(round(rnorm(1, mean = 20, sd = 10)))
      }
    }
    if(taxa[s] == "parent1"){
      depth[s,thesecol] <- as.integer(c(totcounts, 0L))
    }
    if(taxa[s] == "parent2"){
      depth[s,thesecol] <- as.integer(c(0L, totcounts))
    }
    if(s %in% 3:length(taxa)){
      whichgen <- sample(2, 1)
      if(whichgen == 1){
        depth[s,thesecol] <- as.integer(c(0L, totcounts))
      } else {
        reads <- sample(2, size = totcounts, replace = TRUE)
        depth[s,thesecol] <- as.integer(c(sum(reads == 1), sum(reads == 2)))
      }
    }
  }
}

depth[1:20,]

exampleRAD_mapping <- RADdata(depth, alleles2loc, 
                              data.frame(row.names = loci, Chr = c(1,4,6,9), 
                                         Pos = c(1111020, 33069, 2637920, 5549872)),
                              list(2L), 0.001, alNuc, 2L)
exampleRAD_mapping

save(exampleRAD_mapping, file = "data/exampleRAD_mapping.RData")

# Simulate a triploid F1 pop

alleles2loc <- rep(1:4, each = 2)
alNuc <- c("A", "G", "A", "C", "G", "T", "T", "A")
taxa <- c("parent1", "parent2",
          paste("progeny", formatC(1:100, format='d', width = 3, flag = '0'), sep = ""))
loci <- paste("loc", 1:4, sep = "")
depth <- matrix(NA_integer_, nrow = length(taxa), ncol = length(alleles2loc),
                dimnames = list(taxa, paste(loci[alleles2loc], alNuc, sep = "_")))

gen1 <- c(2, 0, 1, 1, 0, 2, 2, 0) # diploid parent
gen2 <- c(2, 2, 0, 4, 3, 1, 3, 1) # tetraploid parent

for(L in 1:max(alleles2loc)){
  thesecol <- which(alleles2loc == L)
  gen1t <- gen1[thesecol]
  gen2t <- gen2[thesecol]
  gamsam1 <- rep(c(0,1), times = gen1t)
  gamsam2 <- rep(c(0,1), times = gen2t)
  for(s in 1:length(taxa)){
    totcounts <- as.integer(round(rnorm(1, mean = 20, sd = 10)))
    if(totcounts < 0) totcounts <- 0L
    if(s %in% 1:2){ # make sure parents have data
      while(totcounts == 0){
        totcounts <- as.integer(round(rnorm(1, mean = 20, sd = 10)))
      }
    }
    if(taxa[s] == "parent1"){
      geno <- gen1t
    }
    if(taxa[s] == "parent2"){
      geno <- gen2t
    }
    if(s %in% 3:length(taxa)){
      s1 <- sample(gamsam1, 1)
      s2 <- sample(gamsam2, 2)
      geno <- c(sum(c(s1, s2) == 0),
                sum(c(s1, s2) == 1))
    }
    alcounts <- integer(length(thesecol))
    if(sum(geno > 0) == 1){
      alcounts[geno > 0] <- totcounts
    } else {
      thesereads <- sample(seq_len(sum(geno > 0)), size = totcounts,
                           replace = TRUE, prob = geno[geno > 0])
      alcounts[geno > 0] <- sapply(seq_len(sum(geno > 0)),
                                   function(x) sum(thesereads == x))
    }
    depth[s,thesecol] <- alcounts
  }
}

exampleRAD_mapping3x <- RADdata(depth, alleles2loc, 
                              data.frame(row.names = loci, Chr = c(1,4,6,9), 
                                         Pos = c(1111020, 33069, 2637920, 5549872)),
                              list(2L, c(2L, 2L)), 0.001, alNuc, c(2L, 4L, rep(3L, 100)))
# Diploid genome but allotetraploid in there as red herring

# Testing 
exampleRAD_mapping3x <- SetDonorParent(exampleRAD_mapping3x, "parent1")
exampleRAD_mapping3x <- SetRecurrentParent(exampleRAD_mapping3x, "parent2")
exampleRAD_mapping3x <- AddAlleleFreqMapping(exampleRAD_mapping3x,
                                             expectedFreqs = seq(0, 1, length.out = 7),
                                             allowedDeviation = 0.05)
exampleRAD_mapping3x <- AddGenotypeLikelihood(exampleRAD_mapping3x)
exampleRAD_mapping3x <- EstimateParentalGenotypes(exampleRAD_mapping3x)

exampleRAD_mapping3x$likelyGeno_donor
exampleRAD_mapping3x$likelyGeno_recurrent
