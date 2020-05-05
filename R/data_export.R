# functions to export data into formats for use by other software

.chromosome_to_integer <- function(chrvect){
  ## Convert chromosomes as strings to chromosomes as integers ##
  
  if(is.integer(chrvect)) return (chrvect)
  
  if(!class(chrvect) %in% c("numeric", "character", "factor")){
    stop("Unexpected values in chromosome column.")
  }
  
  if(is.numeric(chrvect)){
    if(any(chrvect %% 1 != 0)){
      stop("Unexpected numeric non-integer values in chromosome column.")
    }
    chrnum <- as.integer(chrvect)
  }
  if(is.factor(chrvect)){
    chrvect <- as.character(chrvect)
  }
  if(is.character(chrvect)){
    # check if there are multiple sections with digits in chromosome strings
    weirdstrings <- grep("[[:digit:]+][^[:digit:]]+[[:digit:]+]",
                         chrvect)
    if(length(weirdstrings) > 0){
      warning("Strings naming chromosomes have multiple sections with numbers.")
      cat(chrvect[weirdstrings][1:max(c(10, length(weirdstrings)))],
          sep = "\n")
    }
    # remove everything that is not a number from the character strings
    chrnum <- as.integer(gsub("[^[:digit:]]", "", chrvect))
    # get indices for chromosome names starting with "chr"
    chromInd <- grep("^chr", chrvect, ignore.case = TRUE)
    # get highest chromosome number, and add this to non-"chr" sequences
    if(length(chromInd) > 0 && length(chromInd) < length(chrnum)){
      nChr <- max(chrnum[chromInd], na.rm = TRUE)
      chrnum[-chromInd] <- chrnum[-chromInd] + nChr
    }
    # check that strings with same num are the same
    for(cn in sort(unique(chrnum))){
      thisnumind <- which(chrnum == cn)
      if(length(unique(chrvect[thisnumind])) > 1){
        stop("Chromosome names not resolvable to chromosome numbers; number manually.")
      }
    }
  }
  
  return(chrnum)
}

.ordered_gen_and_loctable <- function(object, numgen, 
                                      omit1allelePerLocus = TRUE,
                                      omitCommonAllele = TRUE,
                                      chromAsInteger = TRUE){
  ## Internal function to make a data frame of loci with one row per     ##
  ## allele.  The data frame, with alignment data, and the matrix of     ##
  ## genotypes are output with alleles in the same order.                ##
  ## object = RADdata object                                             ##
  ## numgen = output of GetWeightedMeanGenotypes                         ##
  ## omit1allelePerLocus = TRUE if the same argument was used for GWMG   ##
  ## omitCommonAllele = TRUE if the same argument was used for GwMG      ##
  ## chromAsInteger = should chromosomes be forced to be integers?       ##
  
  # look up alleles in loc table
  if(omit1allelePerLocus){
    alindex <- object$alleles2loc[-OneAllelePerMarker(object,
                                                      commonAllele = omitCommonAllele)]
  } else {
    alindex <- object$alleles2loc
  }
  
  if(all(c("Chr", "Pos") %in% names(object$locTable))){
    loctable <- object$locTable[alindex,c("Chr", "Pos")]
    # sort the genotypes by chromosome and position
    alOrder <- order(loctable$Chr, loctable$Pos)
    numgen <- numgen[,alOrder]
    loctable <- loctable[alOrder,]
    
    if(chromAsInteger && !is.integer(loctable$Chr)){
      # conversion of chromosomes to integer format if needed
      loctable$Chr <- .chromosome_to_integer(loctable$Chr)
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
                                     omitCommonAllele = TRUE,
                                     onePloidyPerAllele = onePloidyPerAllele)
  
  # look up alleles in loc table
  ord <- .ordered_gen_and_loctable(object, numgen, omit1allelePerLocus = TRUE,
                                   chromAsInteger = TRUE, omitCommonAllele = TRUE)
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
                                    omitCommonAllele = TRUE,
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
                                   omitCommonAllele = TRUE,
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

Export_TASSEL_Numeric <- function(object, file, naIfZeroReads = FALSE,
                                  onePloidyPerAllele = FALSE){
  ## From a RADdata object, export numeric genotypes for TASSEL ##
  
  numgen <- GetWeightedMeanGenotypes(object, minval = 0, maxval = 1,
                                     omit1allelePerLocus = TRUE,
                                     omitCommonAllele = TRUE,
                                     naIfZeroReads = naIfZeroReads,
                                     onePloidyPerAllele = onePloidyPerAllele)
  # file header
  cat(c("<Numeric>",
        paste(c("<Marker>", colnames(numgen)), collapse = "\t")),
      sep = "\n", file = file)
  # genotypes table
  write.table(round(numgen, 3), file = file, row.names = TRUE, col.names = FALSE,
              append = TRUE, sep = "\t", quote = FALSE)
}

Export_polymapR <- function(object, naIfZeroReads = TRUE,
                            progeny = GetTaxa(object)[!GetTaxa(object) %in% 
                                                        c(GetDonorParent(object),
                                                          GetRecurrentParent(object),
                                                          GetBlankTaxa(object))]){
  if(!is(object, "RADdata")){
    stop("RADdata object needed")
  }
  if(length(object$posteriorProb) > 1){
    stop("Only one ploidy allowed for Export_polymapR.")
  }
  
  out <- t(GetProbableGenotypes(object, naIfZeroReads = naIfZeroReads,
                                correctParentalGenos = TRUE)[[1]])
  
  # make sure parents are first
  don <- GetDonorParent(object)
  rec <- GetRecurrentParent(object)
  neworder <- c(don, rec, progeny)
  
  return(out[, neworder])
}

Export_MAPpoly <- function(object, file, pheno = NULL, ploidyIndex = 1,
                          progeny = GetTaxa(object)[!GetTaxa(object) %in% c(GetDonorParent(object),
                                                                            GetRecurrentParent(object),
                                                                            GetBlankTaxa(object))],
                          digits = 3){
  # Confirm that genotype calling has been performed
  if(is.null(object$likelyGeno_donor) || is.null(object$posteriorProb)){
    stop("PipelineMapping2Parents needs to be run before using Export_MAPpoly.")
  }
  # Check phenotypes
  if(!is.null(pheno) && is.null(colnames(pheno))){
    stop("pheno should be a matrix or data frame with column names")
  }
  if(!is.null(pheno) && nrow(pheno) != length(progeny)){
    stop("Need one row of pheno for every progeny.")
  }
  if(!is.null(pheno) && !is.null(row.names(pheno)) && !identical(progeny, row.names(pheno))){
    warning("Please check that progeny vector and rows of pheno are in same order.")
  }
  if(!is.null(pheno) && any(grepl(" ", colnames(pheno)))){
    stop("Phenotype names should not have spaces.")
  }
  # Check progeny names
  if(!all(progeny %in% GetTaxa(object))){
    stop("Not all progeny names found in object.")
  }
  if(any(grepl(" ", progeny))){
    stop("Taxa names should not have spaces.")
  }
  # Determine the ploidy
  if(ploidyIndex > length(object$priorProbPloidies)){
    stop("ploidyIndex should be the index of the desired ploidy within object$priorProbPloidies (not the ploidy itself).")
  }
  ploidy <- object$priorProbPloidies[[ploidyIndex]]
  if(length(ploidy) != 1){
    stop("Export is for autopolyploids only.")
  }
  
  # Get parent genotypes
  donorGen <- object$likelyGeno_donor[as.character(ploidy),]
  recurGen <- object$likelyGeno_recurrent[as.character(ploidy),]
  
  # Identify markers to use
  keepal <- which(!is.na(donorGen) & !is.na(recurGen) & 
    !(donorGen == 0 & recurGen == ploidy) & !(donorGen == ploidy & recurGen == 0))
  keepal <- keepal[!keepal %in% OneAllelePerMarker(object, commonAllele = TRUE)]
  if(any(grepl(" ", GetAlleleNames(object)[keepal]))){
    stop("Allele and locus names should not have spaces.")
  }
  
  # Get chromosome and position
  if(is.null(object$locTable$Chr) || all(is.na(object$locTable$Chr))){
    chrnum <- NA_integer_
  } else {
    chrnum <- .chromosome_to_integer(object$locTable$Chr[object$alleles2loc[keepal]])
  }
  if(is.null(object$locTable$Pos) || all(is.na(object$locTable$Pos))){
    position <- NA_integer_
  } else {
    position <- object$locTable$Pos[object$alleles2loc[keepal]]
  }
  
  # Write file header
  cat(c(paste("ploidy", ploidy),
        paste("nind", length(progeny)),
        paste("nmrk", length(keepal)),
        paste("mrknames", paste(GetAlleleNames(object)[keepal], collapse = " ")),
        paste("indnames", paste(progeny, collapse = " ")),
        paste("dosageP", paste(donorGen[keepal], collapse = " ")),
        paste("dosageQ", paste(recurGen[keepal], collapse = " ")),
        paste("seq", paste(chrnum, collapse = " ")),
        paste("seqpos", paste(position, collapse = " ")),
        paste("nphen", ifelse(is.null(pheno), 0, ncol(pheno))),
        "pheno---------------------------------------"), file = file, sep = "\n")
  # Write phenotypes
  if(!is.null(pheno)){
    for(i in 1:ncol(pheno)){
      cat(paste(colnames(pheno)[i], paste(pheno[,i], collapse = " ")),
          file = file, sep = "\n", append = TRUE)
    }
  }
  cat("geno----------------------------------------", 
      file = file, sep = "\n", append = TRUE) # line 12 + nphen
  
  # Write genotype posterior probabilities
  genotab <- data.frame(rep(GetAlleleNames(object)[keepal], each = length(progeny)),
                        rep(progeny, times = length(keepal)),
                        matrix(round(object$posteriorProb[[ploidyIndex]][, progeny, keepal], digits),
                               byrow = TRUE, nrow = length(progeny) * length(keepal),
                               ncol = ploidy + 1))
  write.table(genotab, file = file, append = TRUE, quote = FALSE,
              col.names = FALSE, row.names = FALSE)
}

Export_GWASpoly <- function(object, file, naIfZeroReads = TRUE){
  # matrix of discrete genotypes
  mygeno <- t(GetProbableGenotypes(object,
                                   omit1allelePerLocus = TRUE,
                                   omitCommonAllele = TRUE,
                                   naIfZeroReads = naIfZeroReads)$genotypes)
  
  # get loci to correspond to these alleles
  locindex <- object$alleles2loc[-OneAllelePerMarker(object,
                                                     commonAllele = TRUE)]
  loctable <- object$locTable
  if(is.null(loctable$Chr)){
    loctable$Chr <- rep(0L, nrow(loctable))
  }
  if(is.null(loctable$Pos)){
    loctable$Pos <- 1:nrow(loctable)
  }
  
  # data frame for export
  outdata <- data.frame(Marker = rownames(mygeno),
                        Chrom = .chromosome_to_integer(loctable$Chr[locindex]),
                        Position = loctable$Pos[locindex],
                        mygeno)
  
  # export
  write.csv(outdata, file = file, row.names = FALSE, quote = FALSE)
}

RADdata2VCF <- function(object, file = NULL, asSNPs = TRUE){
  varok <- attr(object$alleleNucleotides, "Variable_sites_only")
  if(is.null(varok) || varok){
    stop("Complete haplotype information not provided; unable to determine SNP positions.  Use refgenome argument in VCF2RADdata.")
  }
  
  # Determine most probable genotypes, and their ploidies
  temp <- GetProbableGenotypes(object, omit1allelePerLocus = FALSE)
  geno <- temp$genotypes
  pld_index <- temp$ploidy_index
  pld_ind_per_loc <- 
    tapply(pld_index, object$alleles2loc,
           function(x){
             u <- unique(x)
             if(length(u) == 1){
               return(u)
             } else {
               return(as.integer(names(which.max(table(x)))))
             }
           })
  pld_per_loc <- sapply(object$priorProbPloidies, sum)[pld_ind_per_loc]
  
  # Process data with internal RCpp function
  temp <- PrepVCFexport(genotypes, object$alleles2loc, object$alleleDepth,
                        object$alleleNucleotides, object$locTable, pld_per_loc,
                        asSNPs)
  REF <- Biostrings::DNAStringSet(temp$REF)
  ALT <- Biostrings::DNAStringSetList(temp$ALT)
  CHROM <- object$locTable$Chr[temmp$Lookup]
  rr <- GenomicRanges::GRanges(CHROM,
                IRanges::IRanges(start = temp$POS, width = nchar(REF)))
  fixed <- S4Vectors::DataFrame(REF = REF, ALT = ALT)
  cd <- S4Vectors::DataFrame(row.names = GetTaxa(object))
  
  vcf <- VariantAnnotation::VCF(rowRanges = rr, fixed = fixed, colData = cd,
                                collapsed = TRUE) ## need to add header
  VariantAnnotation::geno(vcf)$GT <- temp$GT
  VariantAnnotation::geno(vcf)$AD <- temp$AD
  
}
