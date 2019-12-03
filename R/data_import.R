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
# refgenome should either be an FaFile object or a path to a FASTA file for the genome
# tol indicates how dissimilar read depth can be and still have two alleles
# combined into a tag.
consolidateSNPs <- function(alleleDepth, alleles2loc, locTable, alleleNucleotides,
                            tagsize = 80, refgenome = NULL, tol = 0.01){
  if(!all(c("Chr", "Pos") %in% names(locTable))){
    stop("locTable needs Chr for chromosome and Pos for position")
  }
  
  # Set up reference genome if provided.  Make FaFile object.
  if(!is.null(refgenome) && !is(refgenome, "FaFile")){
    if(!is.character(refgenome)){
      stop("refgenome must be character string or FaFile object.")
    }
    if(!file.exists(refgenome)){
      stop(paste("File", refgenome, "not found."))
    }
    if(!file.exists(sprintf("%s.fai", refgenome))){
      # index the fasta file if necessary
      message("Creating FASTA file index...")
      Rsamtools::indexFa(refgenome)
    }
    refgenome <- Rsamtools::FaFile(refgenome)
  }
  # get chromosome names in FASTA
  if(!is.null(refgenome)){
    FAchrnames <- GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(refgenome))
  }
  
  nLoc <- dim(locTable)[1] # number of loci to start
  # rows in locTable for each chromosome
  rowsByChr <- tapply(1:nLoc, locTable$Chr, function(x) x) 
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
                            Pos = integer(nLoc),
                            stringsAsFactors = FALSE)
  
  # variables to keep track of which allele and locus we are on
  currAlOut <- 1L
  currLocOut <- 1L
  
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
    # try consolidating read depth with existing new markers
    cons <- consolidateDepth(depth1[notHomoz,, drop = FALSE], 
                             depth2[notHomoz,, drop = FALSE], alMatch)
    # find individuals that weren't fully matched
    nm <- rowSums(cons) < rowSums(depth1[notHomoz,, drop = FALSE])
    nmInd <- notHomoz[nm]
    cons <- cons[nm,, drop = FALSE]
    # add new alleles until we can match everything
    progress <- TRUE
    while(length(nmInd) > 0 && progress){
      progress <- FALSE
      # Temporary copies of depth matrices that we can subtract from
      depth1T <- depth1[nmInd,, drop = FALSE]
      depth2T <- depth2[nmInd,, drop = FALSE]
      # loop through non-matched individuals
      for(i in 1:length(nmInd)){
        # loop through known new alleles for this ind.
        for(newAl in which(cons[i,] > 0)){
          # corresponding alleles in old markers
          al1 <- alMatch[newAl, 1]
          al2 <- alMatch[newAl, 2]
          # remove reads for this allele
          depth1T[i, al1] <- depth1T[i, al1] - cons[i, newAl]
          depth2T[i, al2] <- depth2T[i, al2] - cons[i, newAl]
        }
      }
      
      # look for matches with what's left
      newAlMatch <- unique(rbind(findMatchesHomoz(depth1T, depth2T),
                                 findMatchesHomoz(depth2T, depth1T)[,2:1]))
      # add these new matches to the list of new alleles
      if(dim(unique(rbind(alMatch, newAlMatch)))[1] > dim(alMatch)[1]){
        progress <- TRUE
        alMatch <- unique(rbind(alMatch, newAlMatch))
        
      } else { # if that didn't work, try looking for alleles with equal depth
        
        # get rid of rows where multiple alleles have same depth
        keepInd <- apply(depth1T, 1, function(x) length(unique(x)) == length(x)) &
          apply(depth2T, 1, function(x) length(unique(x)) == length(x))
        depth1T <- depth1T[keepInd,, drop = FALSE]
        depth2T <- depth2T[keepInd,, drop = FALSE]
        # alleles to potentially match
        als1 <- which(colSums(depth1T) > 0)
        als2 <- which(colSums(depth2T) > 0)
        for(al1 in als1){
          for(al2 in als2){
            # individuals with equal read depth for these two alleles
            if(sum(depth1T[, al1] == depth2T[, al2] & depth1T[, al1] > 0) > 0){
              newAlMatch <- c(al1, al2)
              if(dim(unique(rbind(alMatch, newAlMatch)))[1] > dim(alMatch)[1]){
                alMatch <- rbind(alMatch, newAlMatch)
                progress <- TRUE
              }
            }
          }
        }
        
      }
      
      if(progress){
        # retry consolidation
        cons <- consolidateDepth(depth1[nmInd,, drop = FALSE], 
                                 depth2[nmInd,, drop = FALSE], alMatch)
        nm <- rowSums(cons) < rowSums(depth1[nmInd,, drop = FALSE])
        nmInd <- nmInd[nm]
        cons <- cons[nm,, drop = FALSE]
      }
    }

    return(alMatch)
  }
  
  # Function to consolidate read depth across two markers being merged into one.
  # depth1 and depth2 are read depth matrices for the two markers, respectively,
  # with individuals in rows and alleles in columns.
  # alMatch is a matrix with two columns, with one row for each new allele.  The
  # values indicate the corresponding alleles for the first and second markers.
  consolidateDepth <- function(depth1, depth2, alMatch){
    # depth to output for the new marker
    newdepth <- matrix(0L, nrow = dim(depth1)[1],
                       ncol = dim(alMatch)[1])
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
        indWithAl <- which(depth2[,al2] > 0)
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
    if(!is.null(refgenome)){ # matching chromosome name in ref. genome
      if(thisChrom %in% FAchrnames){
        thisFAchr <- thisChrom
      } else {
        thisFAchr <- grep(paste("(Chr|chr|CHR)(omosome|OMOSOME)?", thisChrom, sep = ""),
                          FAchrnames, value = TRUE)
      }
      if(length(thisFAchr) == 0){
        warning(paste("Couldn't match", thisChrom, "to reference genome."))
      }
      if(length(thisFAchr) > 1){
        stop(paste(thisChrom, "matches multiple sequence names in reference genome."))
      }
    }
    message(paste("Phasing", length(chrset), "SNPs on chromosome", thisChrom))
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
        
        # get proportion difference in depth between these two loci
        diff <- sum(abs(rowSums(thisDepth) - rowSums(lastDepth))) / 
          ((sum(thisDepth) + sum(lastDepth))/2)
      }
      
      if(length(chrset) == 1 || diff > tol || thisPos - lastPos + 1 > tagsize 
         || currLocIn > length(chrset)){
        ## If these are different loci (either due to counts or distance)
        ## put the "last" data into the output, and make the new data the last.
        
        thisNAl <- dim(lastDepth)[2] # number of alleles for locus to output
        thisAlOut <- (1:thisNAl) + currAlOut - 1L # indices for all alleles to output
        
        # add data to output objects
        alleleDepthOut[,thisAlOut] <- lastDepth
        dimnames(alleleDepthOut)[[2]][thisAlOut] <- dimnames(lastDepth)[[2]]
        alleles2locOut[thisAlOut] <- currLocOut
        alleleNucleotidesOut[thisAlOut] <- lastSeq
        row.names(locTableOut)[currLocOut] <- lastName
        locTableOut$Chr[currLocOut] <- thisChrom
        locTableOut$Pos[currLocOut] <- lastPos
        
        if(length(chrset) > 1){
          # shift "this" locus to "last" locus
          lastDepth <- thisDepth
          lastName <- thisName
          lastPos <- thisPos
          lastSeq <- thisSeq
        }
        
        # increment current allele and locus
        currAlOut <- currAlOut + thisNAl
        currLocOut <- currLocOut + 1L
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
        
        # if matching didn't work, skip adding this SNP
        if(dim(alMatch)[1] < dim(lastDepth)[2]) next
        
        # make new set of haplotype sequences
        startPosFromReference <- lastPos + SummarizedExperiment::width(lastSeq[1])
        endPosFromReference <- thisPos - 1
        if(endPosFromReference >= startPosFromReference &&
           !is.null(refgenome) && length(thisFAchr) == 1){
          # retrieve reference sequence and past onto end of lastSeq
          nonvarSeq <- Biostrings::getSeq(refgenome,
                              GenomicRanges::GRanges(seqnames = thisFAchr,
                                IRanges::IRanges(startPosFromReference,
                                                 endPosFromReference),
                                        strand = "+"))[[1]]
          lastSeq <- paste(lastSeq, nonvarSeq, sep = "")
        }
        lastSeq <- paste(lastSeq[alMatch[,1]], thisSeq[alMatch[,2]], sep = "")
        
        # make new depth matrix
#        print(c(lastName, thisName)) # debug
        lastDepth <- consolidateDepth(lastDepth, thisDepth, alMatch)
        dimnames(lastDepth)[[2]] <- paste(lastName, lastSeq, sep = "_")
      }
    } # end of loop through loci
  } # end of loop through chromosomes
  
  # trim output to remove columns not used.
  alleleDepthOut <- alleleDepthOut[, 1:(currAlOut - 1)]
  alleles2locOut <- alleles2locOut[1:(currAlOut - 1)]
  alleleNucleotidesOut <- alleleNucleotidesOut[1:(currAlOut - 1)]
  locTableOut <- locTableOut[1:(currLocOut - 1),]
  
  return(list(alleleDepth = alleleDepthOut, alleles2loc = alleles2locOut,
              alleleNucleotides = alleleNucleotidesOut,
              locTable = locTableOut))
}

# Function to make a function for filtering a VCF file.
# Assumes TASSEL GBSv2 format. (Diploid GT is listed first for each genotype.)
# The output function is used in the prefilters argument for 
# VariantAnnotation::filterVcf.
MakeTasselVcfFilter <- function(min.ind.with.reads = 200,
                                min.ind.with.minor.allele = 10){
  function(lines){
    # vector to indicate the number of individuals with reads for each line
    ind.with.reads <- sapply(gregexpr("[[:blank:]][[:digit:]]/[[:digit:]]:",
                                      lines), 
                             function(x){
                               if(x[1] == -1){
                                 0
                               } else {
                                 length(x)
                               }
                             })
    
    # set up vector to indicate whether each line should be kept
    result <- ind.with.reads >= min.ind.with.reads
    # check on ind with minor allele
    if(min.ind.with.minor.allele > 0){
      ind.with.ref <- sapply(gregexpr("[[:blank:]](0/[[:digit:]]|[[:digit:]]/0):",
                                      lines[result]), 
                             function(x){
                               if(x[1] == -1){
                                 0
                               } else {
                                 length(x)
                               }
                             })
      ind.with.alt <- sapply(gregexpr("[[:blank:]]([123]/[[:digit:]]|[[:digit:]]/[123]):",
                                      lines[result]), 
                             function(x){
                               if(x[1] == -1){
                                 0
                               } else {
                                 length(x)
                               }
                             })
      result[result] <- ind.with.ref >= min.ind.with.minor.allele &
        ind.with.alt >= min.ind.with.minor.allele
    }
    return(result)
  }
}

# Function to import data from a VCF file.
# Reads in VCF using BioConductor, then phases SNPs into haplotypes using 
# consolidateSNPs.
# file can be a character string or a TabixFile.
# samples is a character vector of the names of samples to retain.
# svparam can generally be left as-is unless you have additional columns to 
# import to locTable (from fixed or info) or if you have specific genomic
# ranges to import with "which".
# expectedAlleles and expectedLoci are for after SNP phasing and
# consolidation to haplotypes.
VCF2RADdata <- function(file, phaseSNPs = TRUE, tagsize = 80, refgenome = NULL, 
                        tol = 0.01, al.depth.field = "AD", 
                        min.ind.with.reads = 200, 
                        min.ind.with.minor.allele = 10,
                        possiblePloidies = list(2),
                        contamRate = 0.001,
                        samples = VariantAnnotation::samples(VariantAnnotation::scanVcfHeader(file)),
                        svparam = VariantAnnotation::ScanVcfParam(fixed = "ALT",
                                                                  info = NA, 
                                                                  geno = al.depth.field,
                                                                  samples = samples),
                        yieldSize = 5000, expectedAlleles = 5e5, 
                        expectedLoci = 1e5, maxLoci = NA){
  
  # clean up parameters for import
  if(is.na(VariantAnnotation::vcfFixed(svparam))){
    VariantAnnotation::vcfFixed(svparam) <- "ALT"
  }
  if(length(VariantAnnotation::vcfFixed(svparam)) > 0 &&
     !"ALT" %in% VariantAnnotation::vcfFixed(svparam)){
    VariantAnnotation::vcfFixed(svparam) <- 
      c(VariantAnnotation::vcfFixed(svparam), "ALT")
  }
  if(length(VariantAnnotation::vcfGeno(svparam)) > 0 && 
     is.na(VariantAnnotation::vcfGeno(svparam)[1])){
    stop("geno field must be provided in svparam.")
  }
  if(!identical(VariantAnnotation::vcfGeno(svparam), "AD")){
    warning("Allele depth field not set to AD.")
  }
  if(!is.na(VariantAnnotation::vcfSamples(svparam)[1])){
    samples <- VariantAnnotation::vcfSamples(svparam)
  }
  if(!all(samples %in% VariantAnnotation::samples(VariantAnnotation::scanVcfHeader(file)))){
    stop("Not all samples given in argument found in file.")
  }
  if(length(samples) < min.ind.with.reads){
    stop("Need to adjust min.ind.with.reads for number of samples in dataset.")
  }
  if(length(samples)/2 < min.ind.with.minor.allele){
    stop("Need to adjust min.ind.with.minor.allele for number of samples in dataset.")
  }
  
  # determine what extra columns to add to locTable
  extracols <- c(VariantAnnotation::vcfFixed(svparam),
                 VariantAnnotation::vcfInfo(svparam))
  hdr <- VariantAnnotation::scanVcfHeader(file)
  if(length(VariantAnnotation::vcfFixed(svparam)) == 0){
    extracols <- c("QUAL", "FILTER", extracols)
  }
  if(length(VariantAnnotation::vcfInfo(svparam)) == 0){
    extracols <- c(extracols, row.names(VariantAnnotation::info(hdr)))
  }
  extracols <- extracols[!extracols %in% c("CHROM", "POS", "ID", "REF", "ALT")]
  extracols <- extracols[!is.na(extracols)]
  # preallocate objects for constructing RADdata object
  alleleDepth <- matrix(0L, nrow = length(samples), ncol = expectedAlleles,
                        dimnames = list(samples, 1:expectedAlleles))
  locTable <- data.frame(row.names = as.character(1:expectedLoci),
                         Chr = character(expectedLoci), 
                         Pos = integer(expectedLoci),
                         matrix(nrow = expectedLoci, ncol = length(extracols),
                                dimnames = list(NULL, extracols)),
                         stringsAsFactors = FALSE)
  alleles2loc <- integer(expectedAlleles)
  alleleNucleotides <- character(expectedAlleles)
  # to track which allele we're on
  currAl <- 0L
  currLoc <- 0L
  
  # create Tabix file
  if(is.character(file) && !endsWith(file, ".bgz")){
    if(file.exists(paste(file, ".bgz", sep = ""))){
      file <- paste(file, ".bgz", sep = "")
    } else {
      message("Compressing file with bgzip.")
      tempbgz <- tempfile(fileext = ".bgz") # temporary file for example
      file <- Rsamtools::bgzip(file, dest = tempbgz)
      message(paste("Compressed VCF sent to", file))
    }
  }
  if(is.character(file) && !file.exists(paste(file, ".tbi", sep = ""))){
    message("Indexing VCF.")
    Rsamtools::indexTabix(file, format = "vcf")
  }
  tfile <- Rsamtools::TabixFile(file, yieldSize = yieldSize)
  
  # set up genome argument if needed
  genome <- VariantAnnotation::seqinfo(hdr)
  if(length(genome) == 0 && length(VariantAnnotation::vcfWhich(svparam)) > 0){
    genome <- unique(names(VariantAnnotation::vcfWhich(svparam)))
    genome <- GenomeInfoDb::Seqinfo(seqnames = genome)
  }

  # Read data one chunk at a time
  open(tfile)
  message("Reading file...")
  while(nrow(vcf <- VariantAnnotation::readVcf(tfile, genome = genome, 
                                               param = svparam))){
    thisNloc <- nrow(vcf) # number of loci in this chunk
    message("Unpacking data from VCF...")
    # reference alleles
    thisRef <- as.character(VariantAnnotation::ref(vcf))
    # alternative alleles; clean out ones w/o alt allele
    thisAltList <- lapply(VariantAnnotation::alt(vcf), function(x){
      char <- as.character(x)
      char <- char[char != ""]
      return(char)}) 
    nAlt <- sapply(thisAltList, length) # n alt alleles per locus
    thisNallele <- thisNloc + sum(nAlt) # n alleles in this chunk
    thisAlt <- unlist(thisAltList)
    
    # put reference and alternative alleles together into alleleNucleotides
    thisAlleleNucleotides <- character(thisNallele)
    alsums <- cumsum(nAlt + 1)
    refpos <- c(1, alsums[-thisNloc] + 1)
    thisAlleleNucleotides[refpos] <- thisRef
    thisAlleleNucleotides[-refpos] <- thisAlt
    # set up alleles2loc
    thisAlleles2loc <- rep(1:thisNloc, times = nAlt + 1)
    # set up locTable
    thisLocTable <- data.frame(row.names = make.unique(row.names(vcf)),
                               Chr = as.character(SummarizedExperiment::seqnames(vcf)),
                               Pos = Biostrings::start(vcf),
                               S4Vectors::mcols(vcf)[,extracols],
                               stringsAsFactors = FALSE)
    # set up depth matrix
    thisDepthVal <- unlist(VariantAnnotation::geno(vcf)[[al.depth.field]])
    thisAlNames <- paste(row.names(vcf)[thisAlleles2loc],
                         thisAlleleNucleotides, sep = "_")
    if(length(thisDepthVal) == length(samples) * thisNallele){
      thisAlDepth <- matrix(thisDepthVal, 
                            nrow = length(samples), ncol = thisNallele,
                            dimnames = list(samples, thisAlNames),
                            byrow = TRUE)
    } else {
      if(length(thisDepthVal) %% length(samples) != 0){
        stop("Unexpected number of allele depth fields.")
      }
      # For issue seen with GATK where there may be more AD values than alleles
      thisAlDepth <- matrix(thisDepthVal, nrow = length(samples),
                            ncol = length(thisDepthVal) %/% length(samples),
                            dimnames = list(samples, NULL),
                            byrow = TRUE)
      thisNADfield <- sapply(VariantAnnotation::geno(vcf)[[al.depth.field]][,1],
                             length)
      ADmismatchLoc <- which(thisNADfield != nAlt + 1)
      # Trim down to right number of alleles, and put in zero reads so locus
      # will be discarded.
      remal <- integer(0)
      for(L in ADmismatchLoc){
        firstal <- ifelse(L == 1, 1L, sum(thisNADfield[1:(L - 1)]) + 1)
        lastal <- sum(thisNADfield[1:L])
        thisAlDepth[, firstal:lastal] <- 0L
        remal <- c(remal, (firstal:lastal)[-(1:(nAlt[L] + 1))])
      }
      thisAlDepth <- thisAlDepth[,-remal]
      colnames(thisAlDepth) <- thisAlNames
    }
    
    # replace NA with zero
    thisAlDepth[is.na(thisAlDepth)] <- 0L
    # how many individuals have each allele
    indperal <- colSums(thisAlDepth > 0)
    # loop to filter markers
    message("Filtering markers...")
    keepLoc <- logical(thisNloc) # should loci be retained?
    remAl <- integer(0) # list of alleles to be removed
    iAl <- 0L # keep track of allele position for last locus
    for(i in 1:thisNloc){
      thiscol <- 1:(nAlt[i] + 1L) + iAl # allele columns for this locus
      # check if it passes filtering
      iDepth <- thisAlDepth[,thiscol, drop = FALSE]
      if(sum(rowSums(iDepth) > 0) < min.ind.with.reads){
        keepLoc[i] <- FALSE
      } else {
        keepLoc[i] <- sum(indperal[thiscol] >= min.ind.with.minor.allele) >= 2
      }
      # get rid of funny stuff from TASSEL-GBSv2
      if(keepLoc[i] && any(grepl("N", thisAlleleNucleotides[thiscol]))){
        keepLoc[i] <- FALSE
      }
      # pad out indels, using VCF convention of listing nucleotide before indel
      if(keepLoc[i] && 
         length(unique(nchar(thisAlleleNucleotides[thiscol]))) > 1){
        thisAN <- thisAlleleNucleotides[thiscol]
        maxWidth <- max(nchar(thisAN))
        for(a in 1:length(thisAN)){
          while(nchar(thisAN[a]) < maxWidth){
            thisAN[a] <- paste(thisAN[a], "-", sep = "")
          }
        }
        thisAlleleNucleotides[thiscol] <- thisAN
      }
      # cut the locus if it does not pass filtering
      if(!keepLoc[i]){
        remAl <- c(remAl, thiscol)
      }
      iAl <- iAl + nAlt[i] + 1L # update last allele index
    }
    # remove cut alleles from allele objects
    if(length(remAl) > 0){
      thisAlleleNucleotides <- thisAlleleNucleotides[-remAl]
      thisAlDepth <- thisAlDepth[, -remAl]
      thisAlleles2loc <- thisAlleles2loc[-remAl]
    }
    # update locTable to reflect cut loci
    thisLocTable <- thisLocTable[keepLoc,]
    # update locus numbers
    thisNloc <- sum(keepLoc)
    if(thisNloc > 0){
      thisAlleles2loc <- rep(1:thisNloc, times = table(thisAlleles2loc))
    } else {
      thisAlleles2loc <- integer(0)
    }
    thisNallele <- length(thisAlleles2loc)
    
    # group SNPs into tags
    if(phaseSNPs && thisNloc > 0){
      # add the last marker to the current set if appropriate
      # (i.e. if the break between file chunks may have been within a tag)
      if(currLoc > 0 && thisLocTable$Chr[1] == locTable$Chr[currLoc] &&
         thisLocTable$Pos[1] - locTable$Pos[currLoc] < tagsize){
        
        thisLocTable <- rbind(locTable[currLoc,], thisLocTable)
        oldAlCol <- which(alleles2loc == currLoc)
        thisAlleleNucleotides <- c(alleleNucleotides[oldAlCol], 
                                   thisAlleleNucleotides)
        thisAlleles2loc <- c(rep(1, length(oldAlCol)),
                             thisAlleles2loc + 1)
        thisAlDepth <- cbind(alleleDepth[,oldAlCol],
                             thisAlDepth)
        currLoc <- currLoc - 1L
        currAl <- currAl - length(oldAlCol)
      }
      
      # perform grouping + phasing
      consTags <- consolidateSNPs(thisAlDepth, thisAlleles2loc, thisLocTable,
                                  thisAlleleNucleotides, tagsize = tagsize,
                                  refgenome = refgenome, tol = tol)
      thisAlDepth <- consTags$alleleDepth
      thisAlleles2loc <- consTags$alleles2loc
      thisAlleleNucleotides <- consTags$alleleNucleotides
      thisLocTable <- consTags$locTable
      thisNloc <- dim(thisLocTable)[1]
      thisNallele <- length(thisAlleles2loc)
    }
    
    if(thisNloc > 0){
      # add data from this chunk to objects for whole dataset
      thisAlCol <- (1:thisNallele) + currAl
      if(currAl + thisNallele > ncol(alleleDepth)){
        message("Exceeded expected number of alleles; processing may slow.")
        alleleDepth <- cbind(alleleDepth[,1:currAl], thisAlDepth)
      } else {
        alleleDepth[,thisAlCol] <- thisAlDepth
      }
      dimnames(alleleDepth)[[2]][thisAlCol] <- dimnames(thisAlDepth)[[2]]
      alleles2loc[thisAlCol] <- thisAlleles2loc + currLoc
      alleleNucleotides[thisAlCol] <- thisAlleleNucleotides
      if(currLoc + thisNloc > nrow(locTable)){
        message("Exceeded expected number of loci; processing may slow.")
        locTable <- rbind(locTable[1:currLoc,], thisLocTable)
      } else {
        locTable[(1:thisNloc) + currLoc, ] <- thisLocTable
      }
      row.names(locTable)[(1:thisNloc) + currLoc] <- row.names(thisLocTable)
      # update position in the output
      currAl <- currAl + thisNallele
      currLoc <- currLoc + thisNloc
    }
    
    # don't loop if we are using "which" in svparam
    if(is.na(yieldSize)) break
    # quit if we have a limit on how many loci to import.
    if(!is.na(maxLoci) && currLoc >= maxLoci) break
    message("Reading file...")
  }
  close(tfile)
  
  if(currAl == 0 || currLoc == 0){
    stop("No loci passed the missing data and allele frequency thresholds.")
  }
  message(paste(currLoc, "loci imported."))
  
  # trim output
  alleleDepth <- alleleDepth[, 1:currAl]
  alleles2loc <- alleles2loc[1:currAl]
  alleleNucleotides <- alleleNucleotides[1:currAl]
  locTable <- locTable[1:currLoc, ]
  # indicate whether non-variable sites included in alleleNucleotides
  attr(alleleNucleotides, "Variable_sites_only") <- is.null(refgenome)

  # build RADdata object
  message("Building RADdata object...")
  radout <- RADdata(alleleDepth, alleles2loc, locTable, possiblePloidies,
                    contamRate, alleleNucleotides)
  message("Merging rare haplotypes...")
  radout <- MergeRareHaplotypes(radout, 
                                min.ind.with.haplotype = min.ind.with.minor.allele)
  radout <- RemoveMonomorphicLoci(radout)
  
  return(radout)
}

# Function to import from Stacks 1.48
readStacks <- function(allelesFile, matchesFolder, version = 2, 
                       min.ind.with.reads = 200,
                       min.ind.with.minor.allele = 10,
                       readAlignmentData = FALSE,
                       sumstatsFile = "populations.sumstats.tsv",
                       possiblePloidies = list(2),
                       contamRate = 0.001){
  # get columns depending on version number
  if(!version %in% c(1,2)){
    stop("Version must be equal to 1 or 2.")
  }
  if(version == 2){
    # for catalog alleles file
    hapcol <- 3
    loccol <- 2
    afcols <- list(NULL, integer(0), character(0), NULL, NULL)
    # for matches file
    mloccol <- 1
    msamcol <- 2
    mhapcol <- 4
    mdepcol <- 5
    mfcols <- list(integer(0), integer(0), NULL, character(0), integer(0),
                   NULL)
    # for sumstats file
    sloccol <- 1
    schrcol <- 2
    sposcol <- 3
    sfcols <- list(integer(0), character(0), integer(0), NULL, NULL, NULL,
                   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                   NULL, NULL, NULL, NULL, NULL, NULL)
  } else {
    # for catalog alleles file
    hapcol <- 4
    loccol <- 3
    afcols <- list(NULL, NULL, integer(0), character(0), NULL, NULL)
    # for matches file
    mloccol <- 3
    msamcol <- 4
    mhapcol <- 6
    mdepcol <- 7
    mfcols <- list(NULL, NULL, integer(0), integer(0), NULL, character(0),
                   integer(0), NULL)
    # for catalog tags file
    tloccol <- 3
    tchrcol <- 4
    tposcol <- 5
    tstrcol <- 6
    tfcols <- list(NULL, NULL, integer(0), character(0), integer(0),
                   character(0), NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                   NULL)
  }
  
  # read in catalog.alleles.tsv file
  af <- scan(allelesFile, what = afcols,
             sep = "\t", comment.char = "#", na.strings = character(0))
  # get locus names (numbers) and haplotypes for variable sites
  keep <- af[[hapcol]] != ""
  locNames <- af[[loccol]][keep]
  alleleNucleotides <- af[[hapcol]][keep]
  attr(alleleNucleotides, "Variable_sites_only") <- TRUE
  alleleNames <- paste(locNames, alleleNucleotides, sep = "_")
  
  # get all of the sstacks files to read
  sstacksFiles <- list.files(matchesFolder, "\\.matches\\.tsv")
  if(length(sstacksFiles) == 0){
    stop("No .matches.tsv files found.")
  }
  sampleNames <- sub("\\.matches\\.tsv(\\.gz)?$", "", sstacksFiles)
  sstacksFiles <- file.path(matchesFolder, sstacksFiles)
  if(length(sampleNames) < min.ind.with.reads){
    stop("Need to adjust min.ind.with.reads for number of samples in dataset.")
  }
  if(length(sampleNames)/2 < min.ind.with.minor.allele){
    stop("Need to adjust min.ind.with.minor.allele for number of samples in dataset.")
  }
  
  # set up allele depth matrix
  alleleDepth <- matrix(0L, nrow = length(sampleNames), 
                        ncol = length(alleleNames),
                        dimnames = list(sampleNames, alleleNames))
  # vector for reordering samples according to numbers in matches files
  reorder <- integer(length(sampleNames))
  
  # read sstacks files
  for(i in 1:length(sampleNames)){
    mf <- scan(sstacksFiles[i], 
               what = mfcols, sep = "\t", comment.char = "#", 
               na.strings = character(0))
    keep <- mf[[mhapcol]] != "consensus"
    theseLocNames <- mf[[mloccol]][keep]
    theseAlNuc <- mf[[mhapcol]][keep]
    theseDepth <- mf[[mdepcol]][keep]
    theseAlNames <- paste(theseLocNames, theseAlNuc, sep = "_")
    alleleDepth[i, theseAlNames] <- theseDepth
    reorder[mf[[msamcol]][1]] <- i
  }
  # eliminate any  sample numbers that were skipped, and reorder matrix
  reorder <- reorder[!is.na(reorder) & reorder > 0]
  alleleDepth <- alleleDepth[reorder,]
  
  # filter loci
  alRemove <- integer(0)
  alIndex <- 1L
  indperal <- colSums(alleleDepth > 0)
  while(alIndex <= length(alleleNames)){
    # find allele columns for this locus
    thisLoc <- locNames[alIndex]
    thesecol <- alIndex
    alIndex <- alIndex + 1L
    while(alIndex <= length(alleleNames) && 
          locNames[alIndex] == thisLoc){
      thesecol <- c(thesecol, alIndex)
      alIndex <- alIndex + 1L
    }
    
    # determine whether to keep the locus
    keepLoc <- sum(rowSums(alleleDepth[,thesecol, drop = FALSE]) > 0) >= min.ind.with.reads
    if(keepLoc){
      commonAllele <- which.max(indperal[thesecol])
      keepLoc <- sum(indperal[thesecol[-commonAllele]]) >= min.ind.with.minor.allele
    }
    # update list of alleles to remove
    if(!keepLoc) alRemove <- c(alRemove, thesecol)
  }
  # subset objects
  alleleDepth <- alleleDepth[, -alRemove]
  locNames <- locNames[-alRemove]
  alleleNucleotides <- alleleNucleotides[-alRemove]
  alleleNames <- alleleNames[-alRemove]
  
  # Get chromosome and position
  uniqueLocNames <- unique(locNames)
  if(readAlignmentData){
    if(version == 1){
      tagFile <- sub("alleles", "tags", allelesFile)
      tf <- scan(tagFile, what = tfcols,
                 sep = "\t", comment.char = "#", na.strings = character(0))
      keeprows <- fastmatch::fmatch(uniqueLocNames, tf[[tloccol]])
      locTable <- data.frame(row.names = as.character(uniqueLocNames),
                             Chr = tf[[tchrcol]][keeprows],
                             Pos = tf[[tposcol]][keeprows],
                             Strand = tf[[tstrcol]][keeprows],
                             stringsAsFactors = FALSE)
    }
    if(version == 2){
      sf <- scan(sumstatsFile, what = sfcols,
                 sep = "\t", comment.char = "#", na.strings = character(0))
      keeprows <- fastmatch::fmatch(uniqueLocNames, sf[[sloccol]])
      locTable <- data.frame(row.names = as.character(uniqueLocNames),
                             Chr = sf[[schrcol]][keeprows],
                             Pos = sf[[sposcol]][keeprows],
                             stringsAsFactors = FALSE)
    }
  } else {
    # no alignment data
    locTable <- data.frame(row.names = as.character(uniqueLocNames))
  }
  
  # match depth matrix columns to locus table
  alleles2loc <- fastmatch::fmatch(locNames, uniqueLocNames)
  
  # build RADdata object
  radout <- RADdata(alleleDepth, alleles2loc, locTable, possiblePloidies,
                    contamRate, alleleNucleotides)
  message("Merging rare haplotypes...")
  radout <- MergeRareHaplotypes(radout, 
                                min.ind.with.haplotype = min.ind.with.minor.allele)
  radout <- RemoveMonomorphicLoci(radout)
  return(radout)
}

# Read in from the TASSEL GBSv2 database with minimal processing.
# Use counts matrix output by GetTagTaxaDistFromDBPlugin, plus SAM file.
readTASSELGBSv2 <- function(tagtaxadistFile, samFile, min.ind.with.reads = 200,
                            min.ind.with.minor.allele = 10,
                            possiblePloidies = list(2), contamRate = 0.001,
                            chromosomes = NULL){
  # read the SAM file
  message("Reading SAM file...")
  samwhat <- list(NULL, 0L, "", 0L, NULL, NULL, NULL, NULL, NULL, "", NULL,
                  NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
  samchunk <- scan(samFile, what = samwhat, comment.char = "@", sep = "\t",
                   fill = TRUE, quiet = TRUE)
  samflag <- samchunk[[2]]
  samchr <- samchunk[[3]]
  sampos <- samchunk[[4]]
  samseq <- samchunk[[10]]
  rm(samchunk)
  
  samind <- seq_along(samflag) # index for the tags
  
  # eliminate unaligned tags
  aligned <- bitwAnd(samflag, 4) == 0
  samflag <- samflag[aligned]
  samchr <- samchr[aligned]
  sampos <- sampos[aligned]
  samseq <- samseq[aligned]
  samind <- samind[aligned]
  
  # subset to desired chromosomes
  if(!is.null(chromosomes)){
    chrnotfound <- chromosomes[!chromosomes %fin% samchr]
    if(length(chrnotfound) > 0){
      warning(paste("Chromosomes not found:", paste(chrnotfound, collapse = " ")))
    }
    chrmatch <- samchr %fin% chromosomes
    if(sum(chrmatch) == 0){
      stop("No chromosome names from SAM file match those provided.")
    }
    samflag <- samflag[chrmatch]
    samchr <- samchr[chrmatch]
    sampos <- sampos[chrmatch]
    samseq <- samseq[chrmatch]
    samind <- samind[chrmatch]
  }
  
  # get strand
  samstrand <- ifelse(bitwAnd(samflag, 16) == 0, "top", "bot")
  # get tag sequence length
  samtlen <- nchar(samseq)
  
  # combine to make marker names
  sammrkr <- paste(samchr, sampos, samstrand, samtlen, sep = "-")
  
  # find all non-monomorphic markers
  mrkrTable <- table(sammrkr)
  polymrkr <- names(mrkrTable)[mrkrTable > 1]
  # subset to get only tags belonging to polymorphic markers
  polytags <- sammrkr %fin% polymrkr
  samflag <- samflag[polytags]
  samchr <- samchr[polytags]
  sampos <- sampos[polytags]
  samseq <- samseq[polytags]
  samind <- samind[polytags]
  samstrand <- samstrand[polytags]
  samtlen <- samtlen[polytags]
  sammrkr <- sammrkr[polytags]
  
  # build locTable
  lookup <- fastmatch::fmatch(polymrkr, sammrkr)
  locTable <- data.frame(row.names = polymrkr, Chr = samchr[lookup],
                         Pos = sampos[lookup], Strand = samstrand[lookup],
                         SequenceLength = samtlen[lookup],
                         stringsAsFactors = FALSE)
  locTable <- locTable[order(locTable$Chr, locTable$Pos),]
  
  # open tagtaxadist file
  message("Reading TagTaxaDist file...")
  ttdcon <- file(tagtaxadistFile, open = 'rt')
  # get taxa names
  taxa <- scan(ttdcon, what = character(), sep = "\t", nlines = 1, 
               quiet = TRUE)[-1]
  # set up matrix to contain read counts
  alleleDepth <- matrix(0L, nrow = length(taxa), ncol = length(sammrkr),
                        dimnames = list(taxa, paste(sammrkr, samseq, sep = "_")))
  # setup for reading file in a loop
  whatlist <- list(0L)[rep(1, length(taxa))]
  whatlist <- c(list(""), whatlist) ## might change "" to NULL if we don't care about seq
  chunksize <- 1e5
  # read first chunk
  dataIn <- scan(ttdcon, what = whatlist, sep = "\t", nlines = chunksize,
                 quiet = TRUE)
  nread <- length(dataIn[[1]]) # number of tags read
  lastTag <- 0
  # loop through file
  while(nread){
    # convert imported data to matrix
    datamat <- matrix(unlist(dataIn[-1]), nrow = length(taxa), ncol = nread,
                      byrow = TRUE)
    # find the tags we will keep
    thesetags <- samind > lastTag & samind <= lastTag + nread
    # add to matrix
    alleleDepth[, thesetags] <- datamat[, samind[thesetags] - lastTag]
    # read next chunk
    dataIn <- scan(ttdcon, what = whatlist, sep = "\t", nlines = chunksize,
                   quiet = TRUE)
    nread <- length(dataIn[[1]]) # number of tags read
    lastTag <- lastTag + nread
  }
  # close tagtaxadist file
  close(ttdcon)
  
  # filter loci
  message("Filtering loci...")
  lociToRemove <- integer(0)
  tagsToRemove <- integer(0)
  for(i in 1:nrow(locTable)){
    thesealleles <- which(sammrkr == row.names(locTable)[i])
    thesepres <- alleleDepth[,thesealleles] > 0 # presence or absence of alleles in individuals
    if(sum(rowSums(thesepres) > 0) < min.ind.with.reads ||
       sum(colSums(thesepres) >= min.ind.with.minor.allele) < 2){
      lociToRemove <- c(lociToRemove, i)
      tagsToRemove <- c(tagsToRemove, thesealleles)
    }
  }
  if(length(lociToRemove) > 0){
    locTable <- locTable[-lociToRemove,]
    alleleDepth <- alleleDepth[,-tagsToRemove]
    sammrkr <- sammrkr[-tagsToRemove]
    samseq <- samseq[-tagsToRemove]
  }
  if(nrow(locTable) == 0) stop("No loci passed filtering threshold.")
  
  # build vector to match alleles to loci
  alleles2loc <- fastmatch::fmatch(sammrkr, rownames(locTable))
  # sort alleles to keep things tidy
  tagorder <- order(alleles2loc)
  alleles2loc <- sort(alleles2loc)
  alleleDepth <- alleleDepth[,tagorder]
  samseq <- samseq[tagorder]
  
  # build RADdata object
  message("Building RADdata object...")
  attr(samseq, "Variable_sites_only") <- FALSE
  radout <- RADdata(alleleDepth, alleles2loc, locTable, possiblePloidies,
                    contamRate, samseq)
  radout <- MergeRareHaplotypes(radout, 
                                min.ind.with.haplotype = min.ind.with.minor.allele)
  radout <- RemoveMonomorphicLoci(radout)
  return(radout)
}

# Function to read output of process_sam_multi.py, so the user can get a
# preliminary look at Hind/He distribution before going forward with
# the isolocus_splitter script.
readProcessSamMulti <- function(alignfile, depthfile = sub("align", "depth", alignfile),
                                expectedLoci = 1000, min.ind.with.reads = 200,
                                min.ind.with.minor.allele = 10, possiblePloidies = list(2),
                                contamRate = 0.001, expectedAlleles = expectedLoci * 15,
                                maxLoci = expectedLoci){
  # read file headers
  aligncon <- file(alignfile, open = 'r')
  depthcon <- file(depthfile, open = 'r')
  alignheader <- scan(aligncon, sep = ",", nlines = 1, what = character(),
                      quiet = TRUE)
  depthheader <- scan(depthcon, sep = ",", nlines = 1, what = character(),
                      quiet = TRUE)
  nalign <- (length(alignheader) - 1) / 3 # number of alignment positions reported
  samples <- depthheader[-1]
  nsam <- length(samples)
  
  # set up objects to go into RADdata object
  locTable <- data.frame(Chr = rep(NA_character_, expectedLoci),
                         Pos = rep(NA_integer_, expectedLoci),
                         stringsAsFactors = FALSE)
  locnames <- rep(NA_character_, expectedLoci)
  alleles2loc <- rep(NA_integer_, expectedAlleles)
  alleleNucleotides <- rep(NA_character_, expectedAlleles)
  alleleDepth <- matrix(0L, nrow = nsam, ncol = expectedAlleles,
                        dimnames = list(samples, NULL))
  allelenames <- rep(NA_character_, expectedAlleles)
  
  # count loci and alleles read in (passing any filtering)
  loccount <- 0L
  alcount <- 0L
  
  # read the files
  nscan <- 1e4 # number of lines to read at once
  whatlistalign <- list(character(), integer(), NULL)
  whatlistalign <- whatlistalign[c(rep(1, nalign + 1), rep(2, nalign), rep(3, nalign))]
  whatlistdepth <- list(character(), integer())
  whatlistdepth <- whatlistdepth[c(1, rep(2, nsam))]
  # dummy objects; will hold the last partial marker read
  lastdepthmat <- matrix(integer(0), nrow = 0, ncol = nsam)
  lastaligndata <- whatlistalign
  lastchunk <- FALSE # indicates if we have read the last chunk of the file
  while(loccount < maxLoci && !lastchunk){
    aligndata <- scan(aligncon, sep = ",", what = whatlistalign, nlines = nscan,
                      quiet = TRUE)
    depthdata <- scan(depthcon, sep = ",", what = whatlistdepth, nlines = nscan,
                      quiet = TRUE)
    nread <- length(aligndata[[1]])
    lastchunk <- nread == 0 # is this the last chunk of the file?
    stopifnot(nread == length(depthdata[[1]]))
    if(lastchunk){ # end of file
      depthmat <- lastdepthmat
      aligndata <- lastaligndata
      if(nrow(depthmat) == 0) break
    } else {
      depthmat <- matrix(unlist(depthdata[-1]), nrow = nread, ncol = nsam)
      # merge in any potentially unfinished markers from last chunk
      depthmat <- rbind(lastdepthmat, depthmat)
      for(i in 1:(length(whatlistalign) - nalign)){
        aligndata[[i]] <- c(lastaligndata[[i]], aligndata[[i]])
      }
    }
    nread <- nrow(depthmat) # updated if any loci from last chunk added in
    if(nread < 2) break # quit if there's nothing left to process
    # convert number of mutations to matrix for easier processing
    NM <- matrix(unlist(aligndata[(1:nalign) + (1 + nalign)]),
                 nrow = nread, ncol = nalign)
    
    # find groups of alleles from the same set of alignment positions and process them
    theseal <- 1L
    for(i in 2:nread){
      samegroup <- all(sapply(1:nalign, function(x) aligndata[[x]][i] == aligndata[[x]][i-1]))
      if(samegroup){
        theseal <- c(theseal, i)
      }
      if(i == nread){
        lastaligndata <- lapply(aligndata, function(x) x[theseal])
        lastdepthmat <- depthmat[theseal,]
      }
      if(!samegroup || (lastchunk && i == nread)) { # we have one complete group of alleles
        # divide up into loci based on number of mutations from reference
        thisNM <- NM[theseal,, drop = FALSE]
        thisNM <- thisNM[,!is.na(thisNM[1,]), drop = FALSE]
        hapAssign <- InitHapAssign(thisNM)
        hapAssign <- lapply(1:nalign, function(x) which(hapAssign == x))
        for(L in 1:nalign){
          # skip monomorphic loci or those with no alleles
          if(length(hapAssign[[L]]) < 2) next
          # filter by missing data rate
          thesealSub <- theseal[hapAssign[[L]]]
          if(sum(colSums(depthmat[thesealSub,]) > 0) < min.ind.with.reads) next
          # filter by minor allele frequency
          if(sum(rowSums(depthmat[thesealSub,] > 0) >= min.ind.with.minor.allele) < 2) next
          # add the locus in if it passed everything
          thislocname <-aligndata[[L]][thesealSub[1]]
          firstal <- alcount + 1L
          alcount <- alcount + length(thesealSub)
          if(thislocname %in% locnames){ # same locus may show up in multiple groups
            alleles2loc[firstal:alcount] <- match(thislocname, locnames)
          } else {
            loccount <- loccount + 1L
            locnames[loccount] <- thislocname
            splitlocname <- strsplit(locnames[loccount], "-")[[1]]
            locTable$Chr[loccount] <- splitlocname[1]
            locTable$Pos[loccount] <- as.integer(splitlocname[2])
            alleles2loc[firstal:alcount] <- loccount
          }
          alleleNucleotides[firstal:alcount] <- aligndata[[nalign+1]][thesealSub]
          allelenames[firstal:alcount] <-
            paste(locnames[loccount], alleleNucleotides[firstal:alcount], sep = "_")
          if(alcount > ncol(alleleDepth)){
            alleleDepth <- alleleDepth[,1:(firstal-1)]
            alleleDepth <- cbind(alleleDepth, t(depthmat[thesealSub,]))
          } else {
            alleleDepth[,firstal:alcount] <- t(depthmat[thesealSub,])
          }
          if(loccount >= maxLoci) break
        }
      }
      if(!samegroup){
        theseal <- i
      }
      if(loccount >= maxLoci) break
    }
  }
  
  # wrap-up and build RADdata object
  close(aligncon)
  close(depthcon)
  
  locnames <- locnames[1:loccount]
  locTable <- locTable[1:loccount,]
  allelenames <- allelenames[1:alcount]
  alleleNucleotides <- alleleNucleotides[1:alcount]
  alleleDepth <- alleleDepth[,1:alcount]
  alleles2loc <- alleles2loc[1:alcount]
  
  alleleorder <- order(alleles2loc)
  allelenames <- allelenames[alleleorder]
  alleleNucleotides <- alleleNucleotides[alleleorder]
  alleleDepth <- alleleDepth[,alleleorder]
  alleles2loc <- alleles2loc[alleleorder]
  
  rownames(locTable) <- locnames
  colnames(alleleDepth) <- allelenames
  out <- RADdata(alleleDepth, alleles2loc, locTable, possiblePloidies,
                 contamRate, alleleNucleotides)
#  out <- MergeRareHaplotypes(out, min.ind.with.haplotype = min.ind.with.minor.allele)
#  out <- RemoveMonomorphicLoci(out)
  return(out)
}

readProcessIsoloci <- function(sortedfile, min.ind.with.reads = 200,
                               min.ind.with.minor.allele = 10,
                               min.median.read.depth = 10,
                               possiblePloidies = list(2), contamRate = 0.001,
                               nameFromTagStart = TRUE,
                               mergeRareHap = TRUE){
  message("Reading file...")
  incon <- file(sortedfile, open = "r")
  # read header
  header <- scan(incon, sep = ",", nlines = 1, what = character(), quiet = TRUE)
  samples <- header[-(1:4)]
  nSam <- length(samples)
  scanwhat <- list(character(), integer(), character(), NULL, integer())
  scanwhat <- scanwhat[c(1:4, rep(5, nSam))]
  
  # read file
  mydata <- scan(incon, sep = ",", what = scanwhat, quiet = TRUE)
  close(incon)
  nAl <- length(mydata[[1]])
  
  # get depth matrix
  alleleDepth <- matrix(unlist(mydata[(1:nSam) + 4]), nrow = nSam, ncol = nAl,
                         byrow = TRUE, dimnames = list(samples, NULL))
  if(any(is.na(alleleDepth))){
    stop("Missing data in depth matrix.")
  }
  mydata <- mydata[1:3] # free up space
  # factor by locus, sorting locus names
  alleles2loc_factor <- as.factor(mydata[[1]])
  loci <- levels(alleles2loc_factor)
  nLoc <- length(loci)
  alleles2loc <- as.integer(alleles2loc_factor)
  
  # perform filtering
  message("Filtering and sorting loci...")
  keeploc <- integer(0)
  for(L in 1:nLoc){
    submat <- alleleDepth[,alleles2loc == L]
    depthperind <- rowSums(submat)
    if(sum(depthperind > 0) >= min.ind.with.reads &&
       sum(colSums(submat > 0) >= min.ind.with.minor.allele) > 1 &&
       median(depthperind) >= min.median.read.depth){
      keeploc <- c(keeploc, L)
    }
  }
  keepal <- which(alleles2loc %fin% keeploc)
  
  alleles2loc_factor <- droplevels(alleles2loc_factor[keepal])
  loci <- levels(alleles2loc_factor)
  alleleDepth <- alleleDepth[,keepal, drop = FALSE]
  alleles2loc <- as.integer(alleles2loc_factor)
  alleleNucleotides <- mydata[[3]][keepal]
  
  # sort by locus name (i.e. position and chromosome)
  alorder <- order(alleles2loc)
  alleleDepth <- alleleDepth[, alorder, drop = FALSE]
  alleleNucleotides <- alleleNucleotides[alorder]
  alleles2loc <- alleles2loc[alorder]
  
  # build locTable
  chrom <- sub("\\-.*$", "", loci)
  pos <- mydata[[2]][fastmatch::fmatch(loci, mydata[[1]])]
  if(!nameFromTagStart){
    loci <- paste(chrom, pos, sep = "-")
  }
  locTable <- data.frame(row.names = loci,
                         Chr = chrom, Pos = pos,
                         stringsAsFactors = FALSE)
  
  # build RADdata object
  message("Building RADdata object...")
  attr(alleleNucleotides, "Variable_sites_only") <- FALSE
  radout <- RADdata(alleleDepth, alleles2loc, locTable, possiblePloidies,
                    contamRate, alleleNucleotides)
  if(mergeRareHap){
    radout <- MergeRareHaplotypes(radout, 
                                  min.ind.with.haplotype = min.ind.with.minor.allele)
    radout <- RemoveMonomorphicLoci(radout)
  }
  return(radout)
}
