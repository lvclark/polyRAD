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
                              split = "_"), function(x) substring(x[1], 2, nchar(x[1])))
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

# function to import read counts from TagDigger
readTagDigger <- function(countfile, includeLoci = NULL, 
                          possiblePloidies = list(2), contamRate = 0.001, 
                          dbfile = NULL, dbColumnsToKeep = NULL,
                          dbChrCol = "Chr", dbPosCol = "Pos",
                          dbNameCol = "Marker name"){
  # read marker database, if applicable
  if(!is.null(dbfile)){
    mydb <- read.csv(dbfile, header = TRUE, stringsAsFactors = FALSE,
                     row.names = make.names(dbNameCol))
    if(!is.null(dbColumnsToKeep)){ # only retain subset of markers
      mydb <- mydb[, make.names(dbColumnsToKeep)]
    }
    if(!all(make.names(c(dbChrCol, dbPosCol)) %in% names(mydb))){
      stop(paste(dbChrCol, "and", dbPosCol, "not found in column names of", 
                 dbfile))
    }
    names(mydb)[match(make.names(dbChrCol), names(mydb))] <- "Chr"
    names(mydb)[match(make.names(dbPosCol), names(mydb))] <- "Pos"
    # subset by loci
    if(!is.null(includeLoci)){
      mydb <- mydb[row.names(mydb) %fin% includeLoci,]
    }
    if(dim(mydb)[1] == 0) stop("includeLoci and mydb don't match")
  }
  # read the counts
#  mycounts <- as.matrix(read.csv(countfile, row.names = 1, header = TRUE))
  # use scan; more code but a lot less processing time
  mycon <- file(countfile, open = 'r')
  myheader <- scan(mycon, what = "", nlines = 1, sep = ",")
  alleles <- myheader[-1]
  nalleles <- length(alleles)
  whatlist <- list(0L)[rep(1, nalleles)]
  names(whatlist) <- alleles
  whatlist <- c(list(Taxa = ""), whatlist)
  mydata <- scan(mycon, what = whatlist, sep = ",")
  close(mycon)
  taxa <- mydata[[1]]
  mydata <- mydata[-1]
  mycounts <- matrix(unlist(mydata), nrow = length(taxa), ncol = nalleles,
                  dimnames = list(taxa, alleles))
  
  # extract marker names 
  mrkrNamesByAl <- sapply(strsplit(alleles, split = "_"),
                          function(x) x[1])
  # subset by loci
  if(!is.null(includeLoci)){
    alToKeep <- mrkrNamesByAl %fin% includeLoci
    if(sum(alToKeep) == 0) stop("countsfile and includeLoci don't match")
    mrkrNamesByAl <- mrkrNamesByAl[alToKeep]
    mycounts <- mycounts[, alToKeep]
  }
  # make locTable if still necessary
  if(is.null(dbfile)){
    mydb <- data.frame(row.names = unique(mrkrNamesByAl))
  }
  # get locTable indices for counts matrix columns
  alleles2loc <- fastmatch::fmatch(mrkrNamesByAl, row.names(mydb))
  # error if there isn't a match in marker names between the two files
  if(any(is.na(alleles2loc))){
    cat(paste("Some markers in", file, "not found in", dbfile), sep = "\n")
    sadMarkers <- unique(mrkrNamesByAl[is.na(alleles2loc)])
    maxMrkrToPrint <- 10
    cat(sadMarkers[1:max(c(maxMrkrToPrint, length(sadMarkers)))], sep = "\n")
    if(length(sadMarkers) > maxMrkrToPrint){
      cat("...", sep = "\n")
    }
    stop(paste("Some markers in", file, "not found in", dbfile))
  }
  
  # get variable nucleotides for tags
  myNT <- sapply(strsplit(dimnames(mycounts)[[2]], split = "_"),
                 function(x) x[2])
  
  
  # make RADdata object
  return(RADdata(mycounts, alleles2loc, mydb, possiblePloidies, contamRate,
                 myNT))
}

