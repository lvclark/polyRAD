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

# Function to consolidate loci imported from VCF back into tags.
# Essentially do phasing on the SNPs in the VCF to make haplotypes, using read depth.
# Assume locTable is already sorted by chromosome and position.
# A reference genome can be provided as FASTA or compressed, indexed FASTA
# from Bioconductor in order to get nucleotides in between variable sites.
# tol indicates how dissimilar read depth can be and still have two alleles
# combined into a tag.
consolidateSNPs <- function(alleleDepth, alleles2loc, locTable, alleleNucleotides,
                            tagsize = 80, refgenome = NULL, tol = 0.01){
  if(!all(c("Chr", "Pos") %in% names(locTable))){
    stop("locTable needs Chr for chromosome and Pos for position")
  }
  
  nLoc <- dim(locTable)[1] # number of loci to start
  rowsByChr <- tapply(1:nLoc, locTable$Chr) # rows in locTable for each chromosome
  nAl <- dim(alleleDepth)[2] # number of alleles to start
  nInd <- dim(alleleDepth)[1] # number of individuals
  
  # preallocate objects to output
  alleleDepthOut <- matrix(0L, nrow = nInd, ncol = nAl,
                           dimnames = list(dimnames(alleleDepth)[[1]],
                                           as.character(1:nAl)))
  alleles2locOut <- integer(nAl)
  alleleNucleotidesOut <- character(nAl)
  locTableOut <- data.frame(row.names = as.character(1:nLoc),
                            Chr = character(nLoc),
                            Pos = integer(nLoc))
  
  # variables to keep track of which allele and locus we are on
  currAlOut <- 1
  currLocOut <- 1
  
  ### Internal functions for consolidateSNP ####
  
  # function for finding allele matches by looking for individuals with only one allele.
  # only searches for homozygotes in marker 1; run this function in both directions.
  findMatchesHomoz <- function(depth1, depth2){
    # which individuals just have one allele for marker 1, 
    # and equal depth for both markers.
    homoz1 <- which(rowSums(depth1 > 0) == 1 & rowSums(depth1) == rowSums(depth2))
    match1 <- integer(0) # to hold allele indices for marker 1
    match2 <- integer(0) # to hold allele indices for marker 2
    for(alCol in 1:dim(depth1)[2]){
      # homozygotes for this particular allele in marker 1
      homoz1Al <- homoz1[depth1[homoz1,alCol] > 0]
      # alleles in marker 2 that have reads for those homozygotes
      thisAlMatch <- which(colSums(depth2[homoz1Al,,drop=FALSE]) > 0)
      # add them to match output
      nNewMatches <- length(thisAlMatch)
      match1 <- c(match1, rep(alCol, nNewMatches))
      match2 <- c(match2, thisAlMatch)
    }
    return(matrix(c(match1, match2), nrow = length(match1), ncol = 2))
  }
  
  # function to identify remaining matches not found using homozygotes
  findRemainingMatches <- function(depth1, depth2, alMatch){
    # find individuals with more than one allele for both markers, 
    # and equal depth for both markers.
    notHomoz <- which(rowSums(depth1 > 0) > 1 & rowSums(depth2 > 0) > 1 &
                        rowSums(depth1) == rowSums(depth2))
    notdone <- TRUE
    while(notdone){ # If there are unmatched alleles and progress is being 
                    # made, it will keep going.
      notdone <- FALSE
      for(ind in notHomoz){
        counts1 <- depth1[ind,] # counts for this individual
        counts2 <- depth2[ind,]
        # try to find allele combos that match, remove them from counts
        for(m in 1:dim(alMatch)[1]){
          al1 <- alMatch[m,1] # allele numbers for this pair
          al2 <- alMatch[m,2]
          if(counts1[al1] > 0 && counts2[al2] > 0){
            toSubtract <- min(c(counts1[al1], counts2[al2]))
            counts1[al1] <- counts1[al1] - toSubtract
            counts2[al2] <- counts2[al2] - toSubtract
          }
        }
        # if we have gotten it resolvable, add to match list
        rem1 <- which(counts1 > 0)
        rem2 <- which(counts2 > 0)
        if(length(rem1) == 1 && length(rem2) > 0){
          alMatch <- rbind(alMatch,
                           matrix(c(rep(rem1, length(rem2)), rem2),
                                  nrow = length(rem2), ncol = 2))
          notdone <- TRUE
        } else {
          if(length(rem1) > 1 && length(rem2) == 1){
            alMatch <- rbind(alMatch,
                             matrix(c(rem1, rep(rem2, length(rem1))),
                                    nrow = length(rem1), ncol = 2))
            notdone <- TRUE
          }
        }
      }
      # check if all alleles in both markers have a match
      if(all(1:dim(depth1)[2] %in% alMatch[,1]) &&
         all(1:dim(depth2)[2] %in% alMatch[,2])) notdone <- FALSE
    }
    alMatch <- unique(alMatch, MARGIN = 1) # remove duplicate matches
    return(alMatch)
  }
  
  # Function to consolidate read depth across two markers being merged into one.
  # depth1 and depth2 are read depth matrices for the two markers, respectively,
  # with individuals in rows and alleles in columns.
  # alMatch is a matrix with two columns, with one row for each new allele.  The
  # values indicate the corresponding alleles for the first and second markers.
  consolidateDepth <- function(depth1, depth2, alMatch){
    # depth to output for the new marker
    newdepth <- matrix(0L, nrows = dim(depth1)[1],
                       nrol = dim(alMatch)[1])
    # boolean to figure out if an allele has been processed yet
    allelesDone <- rep(FALSE, dim(alMatch)[1])
    # whether progress has been made
    progress <- TRUE
    
    while(!all(allelesDone) && progress){
      progress <- FALSE
      # loop through alleles for the first marker
      for(al1 in 1:dim(depth1)[2]){
        # which new alleles match this original allele
        newAl <- which(alMatch[,1] == al1 & !allelesDone)
        # skip for now if this allele doesn't correspond to exactly one new, unprocessed allele
        if(length(newAl) != 1) next
        # find individuals with this allele
        indWithAl <- which(depth1[,al1] > 0)
        # add depth to output matrix
        newdepth[indWithAl, newAl] <- depth1[indWithAl, al1]
        # remove these counts from original matrices
        al2 <- alMatch[newAl, 2]
        depth2[indWithAl, al2] <- depth2[indWithAl, al2] - depth1[indWithAl, al1]
        depth1[indWithAl, al1] <- 0L
        # adjust new depth if depth2 was lower than depth1
        d2neg <- depth2[indWithAl, al2] < 0
        newdepth[indWithAl[d2neg], newAl] <- newdepth[indWithAl[d2neg], newAl] +
          depth2[indWithAl[d2neg], al2]
        depth2[indWithAl[d2neg], al2] <- 0L
        
        # mark allele as done
        allelesDone[newAl] <- TRUE
        progress <- TRUE
      }
      # do the same for the second marker
      for(al2 in 1:dim(depth2)[2]){
        newAl <- which(alMatch[,2] == al2 & !allelesDone)
        if(length(newAl) != 1) next
        indWithAl <- which(depth2[,al1] > 0)
        newdepth[indWithAl, newAl] <- depth2[indWithAl, al2]
        al1 <- alMatch[newAl, 1]
        depth1[indWithAl, al1] <- depth1[indWithAl, al1] - depth2[indWithAl, al2]
        depth2[indWithAl, al2] <- 0L
        d1neg <- depth1[indWithAl, al1] < 0
        newdepth[indWithAl[d1neg], newAl] <- newdepth[indWithAl[d1neg], newAl] +
          depth1[indWithAl[d1neg], al1]
        depth1[indWithAl[d1neg], al1] <- 0L
        allelesDone[newAl] <- TRUE
        progress <- TRUE
      }
    }
    
    # Find any individuals with reads remaining to assign.
    # This should only come up if there is homoplasy or recombination within 
    # tags, e.g. alleles AC, AD, BC, and BD all exist.
    indRemaining <- which(rowSums(depth1) > 0 & rowSums(depth2) > 0)
    alRemaining <- which(!allelesDone) # alleles that still need to be assigned
    for(ind in indRemaining){
      als1 <- which(depth1[ind,] > 0) # alleles for marker 1
      als2 <- which(depth2[ind,] > 0) # alleles for marker 2
      # new alleles that are possible matches
      possAl <- alRemaining[alMatch[alRemaining, 1] %in% als1 & 
                              alMatch[alRemaining, 2] %in% als2]
      
      progress <- TRUE
      while(any(depth1[ind,] > 0) && any(depth2[ind,] > 0) && progress){
        progress <- FALSE
        for(al1 in als1){
          if(depth1[ind, al1] <= 0) next
          newAl <- possAl[alMatch[possAl, 1] == al1] # number for new allele(s)
          if(length(newAl) == 1){
            # Resolve cases where only one allele is possible based on 
            # presence/absence.
            al2 <- alMatch[newAl, 2]
          } else {
            # Try to find a match based on depth
            al2 <- which(depth2[ind,] == depth1[ind, al1])
            if(length(al2) != 1) next
            newAl <- newAl[alMatch[newAl, 2] == al2]
            if(length(newAl) != 1) next
          }
          newCount <- min(c(depth1[ind, al1], 
                            depth2[ind, al2]))
          newdepth[ind, newAl] <- newCount
          depth1[ind, al1] <- depth1[ind, al1] - newCount
          depth2[ind, al2] <- depth2[ind, al2] - newCount
          possAl <- possAl[possAl != newAl]
          progress <- TRUE
        }
        for(al2 in als2){ # same procedure with the second marker
          if(depth2[ind, al2] <= 0) next
          newAl <- possAl[alMatch[possAl, 2] == al2]
          if(length(newAl) == 1){
            al1 <- alMatch[newAl, 1]
          } else {
            al1 <- which(depth1[ind,] == depth2[ind, al2])
            if(length(al1) != 1) next
            newAl <- newAl[alMatch[newAl, 1] == al1]
            if(length(newAl) != 1) next
          }
          newCount <- min(c(depth1[ind, al1], 
                            depth2[ind, al2]))
          newdepth[ind, newAl] <- newCount
          depth1[ind, al1] <- depth1[ind, al1] - newCount
          depth2[ind, al2] <- depth2[ind, al2] - newCount
          possAl <- possAl[possAl != newAl]
          progress <- TRUE
        }
      } # end of while loop for unresolved reads
    } # end of loop through unresolved individuals
    return(newdepth)
  } # end of consolidateDepth internal function
  
  ### End internal functions for consolidateSNP ####
  
  # loop through chromosomes
  for(chrset in rowsByChr){
    thisChrom <- locTable$Chr[chrset[1]] # current chromosome name
    cat(paste("Phasing SNPs:", thisChrom), sep = "\n")
    if(!identical(locTable$Pos[chrset], sort(locTable$Pos[chrset]))){
      stop(paste(thisChrom, ": Loci must be sorted by position.", sep = ""))
    }
    currLocIn <- 1 # current locus number index WITHIN chrset (this chromosome)
    # data for this locus
    lastDepth <- alleleDepth[, alleles2loc == chrset[1], drop = FALSE]
    lastName <- row.names(locTable)[chrset[1]]
    lastPos <- locTable$Pos[chrset[1]]
    lastSeq <- alleleNucleotides[alleles2loc == chrset[1]]
    
    # loop through loci on this chromosome
    for(currLocIn in 2:(length(chrset)+1)){
      if(currLocIn <= length(chrset)){
        # data for next locus
        thisDepth <- alleleDepth[, alleles2loc == chrset[currLocIn], drop = FALSE]
        thisName <- row.names(locTable)[chrset[currLocIn]]
        thisPos <- locTable$Pos[chrset[currLocIn]]
        thisSeq <- alleleNucleotides[alleles2loc == chrset[currLocIn]]
      }

      # get proportion difference in depth between these two loci
      diff <- abs(rowSums(thisDepth) - rowSums(lastDepth)) / 
        ((sum(thisDepth) + sum(lastDepth))/2)
      
      if(diff > tol || thisPos - lastPos + 1 > tagsize || currLocIn > length(chrset)){
        ## If these are different loci (either due to counts or distance)
        ## put the "last" data into the output, and make the new data the last.
        
        thisNAl <- dim(lastDepth)[2] # number of alleles for locus to output
        thisAlOut <- (1:thisNAl) + currAlOut - 1 # indices for all alleles to output
        
        # add data to output objects
        alleleDepthOut[,thisAlOut] <- lastDepth
        dimnames(alleleDepthOut)[[2]][thisAlOut] <- dimnames(lastDepth)[[2]]
        alleles2locOut[thisAlOut] <- currLocOut
        alleleNucleotidesOut[thisAlOut] <- lastSeq
        row.names(locTable)[currLocOut] <- lastName
        locTable$Chr[currLocOut] <- thisChrom
        locTable$Pos[currLocOut] <- lastPos
        
        # shift "this" locus to "last" locus
        lastDepth <- thisDepth
        lastName <- thisName
        lastPos <- thisPos
        lastSeq <- thisSeq
        
        # increment current allele and locus
        currAlOut <- currAlOut + thisNAl
        currLocOut <- currLocOut + 1
      } else {
        ## If these are the same locus, merge them.
        
        # matrix to indicate how alleles match up
        alMatch <- matrix(NA_integer_, nrow = 0, ncol = 2,
                          dimnames = list(NULL, c("last", "this")))
        # find individuals in "last" matrix that only have one allele
        alMatch <- rbind(alMatch, findMatchesHomoz(lastDepth, thisDepth))
        # find individuals in "this" matrix that only have one allele
        alMatch <- rbind(alMatch, findMatchesHomoz(thisDepth, lastDepth)[,2:1])
        # reduce to unique set of matches
        alMatch <- unique(alMatch, MARGIN = 1)
        # match any remaining alleles that don't have "homozygotes"
        alMatch <- findRemainingMatches(lastDepth, thisDepth, alMatch)
        
        # make new set of haplotype sequences
        startPosFromReference <- lastPos + nchar(lastSeq[1])
        endPosFromReference <- thisPos - 1
        if(endPosFromReference >= startPosFromReference){ ## add check that we have a reference
          # retrieve reference sequence and past onto end of lastSeq
        }
        newSeq <- paste(lastSeq[alMatch[,1]], thisSeq[alMatch[,2]], sep = "")
        
        # make new depth matrix
        lastDepth <- consolidateDepth(lastDepth, thisDepth, alMatch)
      }
    } # end of loop through loci
  } # end of loop through chromosomes
  # trim output to remove columns not used.
}