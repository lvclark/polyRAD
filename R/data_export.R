# functions to export data into formats for use by other software

.ordered_gen_and_loctable <- function(object, numgen, 
                                      omit1allelePerLocus = TRUE,
                                      chromAsInteger = TRUE){
  ## Internal function to make a data frame of loci with one row per     ##
  ## allele.  The data frame, with alignment data, and the matrix of     ##
  ## genotypes are output with alleles in the same order.                ##
  ## object = RADdata object                                             ##
  ## numgen = output of GetWeightedMeanGenotypes                         ##
  ## omit1allelePerLocus = TRUE if the same argument was used for GWMG   ##
  ## chromAsInteger = should chromosomes be forced to be integers?       ##
  
  # look up alleles in loc table
  if(omit1allelePerLocus){
    alindex <- object$alleles2loc[-OneAllelePerMarker(object)]
  } else {
    alindex <- object$alleles2loc
  }
  
  if(all(c("Chr", "Pos") %in% names(object$Loctable))){
    loctable <- object$locTable[alindex,c("Chr", "Pos")]
    # sort the genotypes by chromosome and position
    alOrder <- order(loctable$Chr, loctable$Pos)
    numgen <- numgen[,alOrder]
    loctable <- loctable[alOrder,]
    
    if(chromAsInteger && !is.integer(loctable$Chr)){
      if(!class(loctable$Chr) %in% c("numeric", "character", "factor")){
        stop("Unexpected values in chromosome column.")
      }
      # conversion of chromosomes to integer format if needed
      if(is.numeric(loctable$Chr)){
        if(any(loctable$Chr %% 1 != 0)){
          stop("Unexpected numeric non-integer values in chromosome column.")
        }
        loctable$Chr <- as.integer(loctable$Chr)
      }
      if(is.factor(loctable$Chr)){
        loctable$Chr <- as.character(loctable$Chr)
      }
      if(is.character(loctable$Chr)){
        # check if there are multiple sections with digits in chromosome strings
        weirdstrings <- grep("[[:digit:]+][^[:digit:]]+[[:digit:]+]",
                             loctable$Chr)
        if(length(weirdstrings) > 0){
          warning("Strings naming chromosomes have multiple sections with numbers.")
          cat(loctable$Chr[weirdstrings][1:max(c(10, length(weirdstrings)))],
              sep = "\n")
        }
        # remove everything that is not a number from the character strings
        chrnum <- as.integer(gsub("[^[:digit:]]", "", loctable$Chr))
        # get indices for chromosome names starting with "chr"
        chromInd <- grep("^chr", loctable$Chr, ignore.case = TRUE)
        # get highest chromosome number, and add this to non-"chr" sequences
        if(length(chromInd) > 0 && length(chromInd) < length(chrnum)){
          nChr <- max(chrnum[chromInd], na.rm = TRUE)
          chrnum[-chromInd] <- chrnum[-chromInd] + nChr
        }
        # check that strings with same num are the same
        for(cn in sort(unique(chrnum))){
          thisnumind <- which(chrnum == cn)
          if(length(unique(loctable$Chr[thisnumind])) > 1){
            stop("Chromosome names not resolvable to chromosome numbers; number manually.")
          }
        }
        loctable$Chr <- chrnum
      }
    }
  } else {
    # insert fake alignment data if it doesn't exist
    loctable <- data.frame(Chr = rep(1, length(alindex)),
                           Pos = alindex)
    cat("Alignment data not found in object.  Creating dummy data", sep = "\n")
  }
  
  return(list(numgen = numgen, loctable = loctable))
}
 
ExportGAPIT <- function(object, onePloidyPerAllele = FALSE){
  ## From a RADdata object, export weighted mean genotypes and alignment ##
  ## information formatted for GAPIT.                                    ##
  
  # get numeric genotypes
  numgen <- GetWeightedMeanGenotypes(object, minval = 0, maxval = 2,
                                     omit1allelePerLocus = TRUE,
                                     onePloidyPerAllele = onePloidyPerAllele)
  
  # look up alleles in loc table
  ord <- .ordered_gen_and_loctable(object, numgen, omit1allelePerLocus = TRUE,
                                   chromAsInteger = TRUE)
  rm(numgen)
  
  taxa <- GetTaxa(object)
  
  # finish converting to GAPIT format
  GD <- data.frame(taxa = taxa, ord$numgen, stringsAsFactors = FALSE)
  GM <- data.frame(Name = dimnames(ord$numgen)[[2]],
                   Chromosome = ord$loctable$Chr,
                   Position = ord$loctable$Pos,
                   stringsAsFactors = FALSE)
  
  return(list(GD = GD, GM = GM))
}

Export_rrBLUP_Amat <- function(object, naIfZeroReads = FALSE,
                               onePloidyPerAllele = FALSE){
  ## From a RADdata object, export weighted mean genotypes formatted ##
  ## for the A.mat function in rrBLUP.                               ##
  
  numgen <- GetWeightedMeanGenotypes(object, minval = -1, maxval = 1,
                                    omit1allelePerLocus = TRUE,
                                    naIfZeroReads = naIfZeroReads,
                                    onePloidyPerAllele = onePloidyPerAllele)
  return(numgen)
}

Export_rrBLUP_GWAS <- function(object, naIfZeroReads = FALSE,
                               onePloidyPerAllele = FALSE){
  ## From a RADdata object, export weighted mean gentoypes formatted ##
  ## for the GWAS function in rrBLUP.                                ##
  
  numgen <- Export_rrBLUP_Amat(object, naIfZeroReads = naIfZeroReads,
                               onePloidyPerAllele = onePloidyPerAllele)
  
  # look up alleles in loc table
  ord <- .ordered_gen_and_loctable(object, numgen, omit1allelePerLocus = TRUE,
                                   chromAsInteger = FALSE)
  rm(numgen)
  
  # set up output data frame
  geno <- data.frame(Marker = dimnames(ord$numgen)[[2]],
                     Chr = ord$loctable$Chr,
                     Pos = ord$loctable$Pos,
                     t(ord$numgen),
                     stringsAsFactors = FALSE)
  return(geno)
}