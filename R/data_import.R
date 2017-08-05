# functions for import of data into RADgen object

# Function to import read counts from UNEAK output (HapMap.hmc.txt).
# includeloci should be a character vector of loci names to keep in output.
# If shortIndNames is TRUE, just keep the part of the individual name before
# the first underscore.
readHMC <- function(file, includeLoci=NULL, shortIndNames=TRUE,
                    possiblePloidies = list(2), contamRate = 0.001,
                    fastafile = sub("hmc.txt", "fas.txt", file, fixed = TRUE)){
  # read in file; not using read.table, to preserve indiv. names
  filelines <- strsplit(readLines(file), split="\t")
  # extract names of individuals
  indstrings <- filelines[[1]][2 :
                                 (match("HetCount_allele1", filelines[[1]]) - 1)]
  if(shortIndNames){
    indstrings <- sapply(strsplit(indstrings, split="_", fixed=TRUE),
                         function(x) x[1])
  }
  # number of individuals
  nInd <- length(indstrings)
  # number of loci
  nLoc <- ifelse(is.null(includeLoci), length(filelines) - 1,
                 length(includeLoci))
  # set up matrix to contain results
  alDepth <- matrix(as.integer(0), nrow = nInd, ncol = 2*nLoc,
                    dimnames = list(indstrings, NULL))
  # loop to fill in data
  locnames <- character(nLoc) # names of loci
  if(is.null(includeLoci)){  # if we are using all loci in the file
    for(i in 1:nLoc){
      locnames[i] <- filelines[[i+1]][1]
      thesecounts <- strsplit(filelines[[i+1]][2:(nInd+1)],
                              split="|", fixed=TRUE)
      alDepth[,2*i - 1] <- sapply(thesecounts, function(x) as.integer(x[1]))
      alDepth[,2*i] <- sapply(thesecounts, function(x) as.integer(x[2]))
    }
  } else { # if we are just using a subset of loci
    locnames <- includeLoci
    for(i in 2:length(filelines)){
      # get locus name for this line of the file and see if it is in
      # the list of loci that we are keeping. If so, what position?
      Lnum <- fastmatch::fmatch(filelines[[i]][1], locnames)
      # go to next line of file if we don't want this locus
      if(is.na(Lnum)) next
      thesecounts <- strsplit(filelines[[i]][2:(nInd+1)],
                              split="|", fixed=TRUE)
      alDepth[,Lnum*2 - 1] <- sapply(thesecounts, function(x) as.integer(x[1]))
      alDepth[,Lnum*2] <- sapply(thesecounts, function(x) as.integer(x[2]))
    }
  }
  dimnames(alDepth)[[2]] <- paste(rep(locnames, each = 2), c(0,1), sep = "_")
  
  # get nucleotides for the alleles
  fastalines <- readLines(fastafile)
  fastaseq <- fastalines[seq(2,length(fastalines), by = 2)] # just sequences
  fastaloc <- sapply(strsplit(fastalines[seq(1,length(fastalines)-1, by = 2)],
                              split = "_"), function(x) x[1])
  alleleNucleotides <- character(nLoc * 2)
  for(i in 1:nLoc){
    rowid <- fastmatch::fmatch(locnames[i], fastaloc)
    theseseq <- strsplit(fastaseq[c(rowid, rowid + 1)], split = "")
    varpos <- theseseq[[1]] != theseseq[[2]]
    alleleNucleotides[2*i - 1] <- theseseq[[1]][varpos]
    alleleNucleotides[2*i] <- theseseq[[2]][varpos]
  }
  
  return(RADdata(alleleDepth = alDepth, alleles2loc = rep(1:nLoc, each = 2),
                 locTable = data.frame(row.names = locnames), 
                 possiblePloidies = possiblePloidies,
                 contamRate = contamRate,
                 alleleNucleotides = alleleNucleotides))
}
