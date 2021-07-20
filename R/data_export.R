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

Export_GWASpoly <- function(object, file, naIfZeroReads = TRUE, postmean = TRUE, digits = 3){
  if(postmean){
    mygeno <- t(GetWeightedMeanGenotypes(object, maxval = max(sapply(object$possiblePloidies, sum)),
                                         omit1allelePerLocus = TRUE,
                                         omitCommonAllele = TRUE,
                                         naIfZeroReads = naIfZeroReads))
    mygeno <- round(mygeno, digits = digits)
  } else {
    # matrix of discrete genotypes
    mygeno <- t(GetProbableGenotypes(object,
                                     omit1allelePerLocus = TRUE,
                                     omitCommonAllele = TRUE,
                                     naIfZeroReads = naIfZeroReads)$genotypes)
  }
  
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

RADdata2VCF <- function(object, file = NULL, asSNPs = TRUE, hindhe = TRUE,
                        sampleinfo = data.frame(row.names = GetTaxa(object)),
                        contigs = data.frame(row.names = unique(object$locTable$Chr))){
  # shortcuts to functions to use
  DataFrame <- S4Vectors::DataFrame
  
  varok <- attr(object$alleleNucleotides, "Variable_sites_only")
  if(is.null(varok) || varok){
    stop("Complete haplotype information not provided; unable to determine SNP positions.  Use refgenome argument in VCF2RADdata.")
  }
  if(is.null(object$locTable$Chr) || is.null(object$locTable$Pos)){
    stop("Need chromosome and position information (Chr and Pos in locTable).")
  }
  if(is.null(object$locTable$Ref)){
    warning("Reference allele not indicated.  Using the major allele for each locus.")
    object$locTable$Ref <- object$alleleNucleotides[OneAllelePerMarker(object, commonAllele = TRUE)]
  }
  if(nrow(sampleinfo) != nTaxa(object) || !all(rownames(sampleinfo) %in% GetTaxa(object))){
    "sampleinfo doesn't match taxa in RADdata object."
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
  temp <- PrepVCFexport(geno, object$alleles2loc, object$alleleDepth,
                        object$alleleNucleotides, object$locTable, pld_per_loc,
                        asSNPs)
  REF <- Biostrings::DNAStringSet(temp$REF)
  ALT <- Biostrings::DNAStringSetList(temp$ALT)
  CHROM <- as.character(object$locTable$Chr)[temp$Lookup]
  rr <- GenomicRanges::GRanges(CHROM,
                IRanges::IRanges(start = temp$POS, 
                                 width = BiocGenerics::width(REF)))
  fixed <- DataFrame(REF = REF, ALT = ALT)
  cd <- DataFrame(row.names = GetTaxa(object))
  if(ncol(sampleinfo) > 0){
    cd <- cbind(cd, sampleinfo[GetTaxa(object),])
    colnames(cd) <- colnames(sampleinfo)
  }
  DP <- t(object$locDepth[,as.character(temp$Lookup)])
  rownames(DP) <- NULL
  info <- DataFrame(NS = rowSums(DP > 0), DP = rowSums(DP), LU = temp$Lookup)
  infohdr <- DataFrame(row.names = c("NS", "DP", "LU"), Number = c("1", "1", "1"),
                       Type = c("Integer", "Integer", "Integer"),
                       Description = c("Number of samples with data", "Combined depth across samples",
                                       "Lookup index of marker in RADdata object"))
  metahdr <- DataFrame()
  if(ncol(contigs) > 0){
    ctg <- DataFrame(row.names = rownames(contigs), contigs)
    colnames(ctg) <- colnames(contigs)
  } else {
    ctg <- DataFrame(row.names = rownames(contigs))
  }
  
  # Add Hind/He if desired
  if(hindhe){
    infohdr <- rbind(infohdr,
                     DataFrame(row.names = "HH", Number = "1", Type = "Float",
                               Description = "Hind/He for the locus in the RADdata object"))
    metahdr <- DataFrame(row.names = "HH", Number = "1", Type = "Float",
                         Description = "Hind/He for the sample, averaged across loci in the RADdata object")
    hh <- HindHe(object)
    info$HH <- colMeans(hh, na.rm = TRUE)[temp$Lookup]
    cd$HH <- rowMeans(hh, na.rm = TRUE)
  }
  
  # Build VCF object
  hdr <- VariantAnnotation::VCFHeader(reference = rownames(ctg),
    samples = GetTaxa(object),
    IRanges::DataFrameList(fileformat = DataFrame(row.names = "fileformat", Value = "VCFv4.3"),
                           fileDate = DataFrame(row.names = "fileDate", Value = gsub("-", "", Sys.Date())),
                           source = DataFrame(row.names = "source", Value = paste0("polyRADv", packageVersion("polyRAD"))),
                           FORMAT = DataFrame(row.names = c("GT", "AD", "DP"),
                                              Number = c("1", "R", "1"), Type = c("String", "Integer", "Integer"),
                                              Description = c("Genotype", "Read depth for each allele", "Read depth")),
                           INFO = infohdr, META = metahdr, contig = ctg))
  if(ncol(cd) > 0){
    VariantAnnotation::meta(hdr)$SAMPLE <- cd
  }
  vcf <- VariantAnnotation::VCF(rowRanges = rr, fixed = fixed, info = info, colData = cd,
                                geno = S4Vectors::SimpleList(GT = temp$GT, AD = temp$AD, DP = DP),
                                exptData = list(header = hdr), collapsed = TRUE)
  
  # output
  if(!is.null(file)){
    VariantAnnotation::writeVcf(vcf, file)
  }
  return(vcf)
}

Export_Structure <- function(object, file, includeDistances = FALSE,
                             extraCols = NULL, missingIfZeroReads = TRUE){
  # sort data if necessary
  if(includeDistances){
    if(is.null(object$locTable$Chr) || is.null(object$locTable$Pos)){
      stop("Chromosome and position needed in locTable if distances are to be output.")
    }
    ord <- order(object$locTable$Chr, object$locTable$Pos)
    if(!identical(ord, seq_len(nLoci(object)))){
      object <- SubsetByLocus(object, ord)
    }
    # get inter marker distances
    distrow <- rep(-1, nLoci(object))
    for(chr in unique(object$locTable$Chr)){
      theserow <- which(object$locTable$Chr == chr)
      distrow[theserow[-1]] <-
        object$locTable$Pos[theserow[-1]] - object$locTable$Pos[theserow[-length(theserow)]]
    }
  }
  # get genotype matrix
  geno <- GetProbableGenotypes(object, omit1allelePerLocus = FALSE,
                               naIfZeroReads = missingIfZeroReads)
  # get ploidy
  pind <- unique(geno$ploidy_index)
  ploidies <- sapply(object$priorProbPloidies[pind], sum)
  ploidy <- max(ploidies)
  stopifnot(all(geno$genotypes <= ploidy, na.rm = TRUE))
  # put data into Structure format (Rcpp function)
  strdata <- FormatStructure(geno$genotypes, object$alleles2loc, ploidy)
  #colnames(strdata) <- GetLoci(object)
  
  # write file
  if(is.null(extraCols)){
    ncolx <- 1
    outdf <- data.frame(rep(GetTaxa(object), each = ploidy),
                        strdata)
  } else {
    ncolx <- 1 + ncol(extraCols)
    if(!is.null(rownames(extraCols)) &&
       !identical(rownames(extraCols), as.character(seq_len(nTaxa(object)))) &&
       !identical(rownames(extraCols), GetTaxa(object))){
      if(!all(GetTaxa(object) %in% rownames(extraCols))){
        stop("Row names of extraCols don't match sample names.")
      }
      extraCols <- extraCols[GetTaxa(object),]
    }
    if(nrow(extraCols) != nTaxa(object)){
      stop("Number of rows in extraCols is different from number of taxa in object.")
    }
    outdf <- data.frame(rep(GetTaxa(object), each = ploidy),
                        extraCols[rep(seq_len(nTaxa(object)), each = ploidy),],
                        strdata)
  }
  cat(paste(c(rep("", ncolx), GetLoci(object)), collapse = "\t"),
      file = file, sep = "\n")
  if(includeDistances){
    cat(paste(c(rep("", ncolx), distrow), collapse = "\t"), file = file,
        sep = "\n", append = TRUE)
  }
  write.table(outdf, file = file, append = TRUE, sep = "\t", col.names = FALSE,
              row.names = FALSE, quote = FALSE)
  cat(c(paste("Number of individuals:", nTaxa(object)),
        paste("Number of loci:", nLoci(object)),
        paste("Ploidy of data:", ploidy),
        "Missing data value: -9", "", "File contains:",
        "Row of marker names"),
      sep = "\n")
  if(includeDistances){
    cat("Map distances between loci", sep = "\n")
  }
  cat("Individual ID for each individual", sep = "\n")
  if(ncolx > 1){
    cat(paste(ncolx - 1, "extra columns that you should define for Structure"))
  }
}

Export_adegenet_genind <- function(object, ploidyIndex = 1){
  object <- SubsetByPloidy(object, ploidies = object$priorProbPloidies[ploidyIndex])
  
  tab <- GetProbableGenotypes(object, omit1allelePerLocus = FALSE,
                              multiallelic = "correct")$genotypes
  colnames(tab) <- paste(sub("\\.", "_", GetLoci(object)[object$alleles2loc]),
                         object$alleleNucleotides, sep = ".")
  
  out <- methods::new("genind", tab = tab,
             ploidy = sum(object$priorProbPloidies[[1]]),
             type = "codom")

  return(out)
}

Export_polymapR_probs <- function(object, maxPcutoff = 0.9,
                                  correctParentalGenos = TRUE,
                                  multiallelic = "correct"){
  if(!is(object, "RADdata")){
    stop("RADdata object needed")
  }
  if(length(object$posteriorProb) > 1){
    stop("Only one ploidy allowed for Export_polymapR_probs.")
  }
  object <- RemoveUngenotypedLoci(object, removeNonvariant = TRUE)
  p1 <- dim(object$posteriorProb[[1]])[1] # ploidy plus one
  omitals <- OneAllelePerMarker(object, commonAllele = TRUE)
  keepals <- GetAlleleNames(object)[-omitals]
  
  probmat <- matrix(object$posteriorProb[[1]][,,keepals],
                    nrow = nTaxa(object) * length(keepals),
                    ncol = p1,
                    dimnames = list(NULL, paste0("P", seq_len(p1) - 1L)),
                    byrow = TRUE)
  out <- data.frame(SampleName = rep(GetTaxa(object), times = length(keepals)),
                    MarkerName = rep(keepals, each = nTaxa(object)),
                    probmat)
  genomat <- GetProbableGenotypes(object, omit1allelePerLocus = TRUE,
                                  omitCommonAllele = TRUE,
                                  correctParentalGenos = correctParentalGenos,
                                  multiallelic = multiallelic)$genotypes
  genovect <- as.vector(genomat)
  out$maxP <- numeric(nTaxa(object) * length(keepals))
  for(i in seq_len(p1) - 1L){
    theserows <- which(genovect == i)
    out$maxP[theserows] <- out[[paste0("P", i)]][theserows]
  }
  out$maxgeno <- genovect
  out$geno <- genovect
  out$geno[out$maxP < maxPcutoff] <- NA_real_
  
  return(out)
}
