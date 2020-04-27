#include <Rcpp.h>
#include <string>
#include "Hap2SNP.h"
using namespace Rcpp;

// Prepare GT strings for the VCF from a matrix of allele copy number.

// [[Rcpp::export]]
StringVector MakeGTstrings(IntegerMatrix genotypes, int ploidy) {
  int nal = genotypes.ncol();
  int nsam = genotypes.nrow();
  IntegerVector thisgen(nal);
  StringVector out(nsam);
  int cp;
  String thisstring;
  
  for(int s = 0; s < nsam; s++){
    thisstring = "";
    thisgen = genotypes(s, _);
    if(any(is_na(thisgen))){
      for(int k = 0; k < ploidy; k++){
        thisstring += "./";
      }
    } else {
      for(int a = 0; a < nal; a++){
        cp = thisgen[a];
        for(int c = 0; c < cp; c++){
          thisstring += a;
          thisstring += "/";
        }
      }
    }

    thisstring.replace_last("/", "");
    out[s] = thisstring;
  }
  return out;
}

// Make a faster lookup for allele indices, for when loci will be looped through
// multiple times.  For internal use; indices start at zero in output.
// alleles2loc is assumed to be from RADdata and have indices starting at one.
List AlleleIndex(IntegerVector alleles2loc, int nloc){
  IntegerVector alleles = seq(0, alleles2loc.size() - 1);
  IntegerVector thesecol;
  List out(nloc);
  
  for(int L = 0; L < nloc; L++){
    thesecol = alleles[alleles2loc == L + 1];
    out[L] = thesecol;
  }
  return out;
}

// Subset a matrix using a vector of column indices
IntegerMatrix SubsetMatrixCol(IntegerMatrix mat, IntegerVector cols){
  int nr = mat.nrow();
  int co = cols.size();
  IntegerMatrix out(nr, co);
  
  for(int r = 0; r < nr; r++){
    for(int c = 0; c < co; c++){
      out(r, c) = mat(r, cols[c]);
    }
  }
  return out;
}

// Do the matrix multiplication with the conversion matrices, to avoid needing
// RcppArmadillo.  Assumes all values in the conversion matrix are 0 and 1, and
// the sum of each row is 1.
IntegerMatrix ConvMatMult(IntegerMatrix origmat, IntegerMatrix convmat){
  int nsam = origmat.nrow();
  int nhap = origmat.ncol();
  int nal = convmat.ncol();
  IntegerMatrix out(nsam, nal);
  
  for(int i = 0; i < nhap; i++){
    for(int j = 0; j < nal; j++){
      if(convmat(i, j) == 1){
        for(int s = 0; s < nsam; s++){
          out(s, j) += origmat(s, i);
        }
        break;
      }
    }
  }
  return out;
}

// Format allelic read depth for the AD slot in a VCF object.
List FormatAD(IntegerMatrix depthmat){
  int nsam = depthmat.nrow();
  IntegerVector thisdep;
  List out(nsam);
  
  for(int i = 0; i < nsam; i++){
    thisdep = depthmat(i, _ );
    out[i] = clone(thisdep);
  }
  return out;
}

// Function to take genotype calls and slots from a RADdata object and prepare
// data for export to VCF.

// [[Rcpp::export]]
List PrepVCFexport(IntegerMatrix genotypes, IntegerVector alleles2loc,
                   IntegerMatrix alleleDepth, StringVector alleleNucleotides,
                   DataFrame locTable, int ploidy, bool asSNPs) {
  int nloc = locTable.nrows();
  //int nsam = genotypes.nrow();
  List alleleLookup = AlleleIndex(alleles2loc, nloc);
  IntegerVector thesecol;
  StringVector thesehap;
  IntegerVector posvect = locTable["Pos"];
  StringVector refvect = locTable["Ref"];
  IntegerMatrix thesegeno;
  IntegerMatrix thesedepths;
  IntegerMatrix thesegeno1;
  IntegerMatrix thesedepths1;
  StringVector theseGT;
  List theseAD;
  int pos;
  int nsubloc;
  std::string ref;
  List hapconv;
  List outpos(nloc);
  List outal(nloc);
  List outgen(nloc);
  List outdepth(nloc);
  List thesemats;
  
  for(int L = 0; L < nloc; L++){
    // Subset data for this locus
    thesecol = alleleLookup[L];
    thesehap = alleleNucleotides[thesecol];
    
    thesegeno = SubsetMatrixCol(genotypes, thesecol);
    thesedepths = SubsetMatrixCol(alleleDepth, thesecol);
    
    pos = posvect[L];
    ref = refvect[L];
    
    // Convert haplotypes to SNPs
    if(asSNPs){
      hapconv = Hap2SNP(thesehap, ref, pos);
    } else {
      hapconv = Hap2Hap(thesehap, ref, pos);
    }
    outpos[L] = hapconv[0];
    outal[L] = hapconv[1];
    // Consider separating alleles into reference and alts.
    thesemats = hapconv[2];
    
    // Format genotypes and depths
    nsubloc = thesemats.size();
    for(int i = 0; i < nsubloc; i++){
      IntegerMatrix thismat = thesemats(i);
      thesegeno1 = ConvMatMult(thesegeno, thismat);
      thesedepths1 = ConvMatMult(thesedepths, thismat);
      theseGT = MakeGTstrings(thesegeno1, ploidy);
      theseAD = FormatAD(thesedepths1);
      // Need to fill these into a matrix or list.
    }
    
  }
  List out = List::create(outpos, outal, theseGT, theseAD);
  return out;
}

// Testing
/*** R
alleles2loc <- c(1,1,1,2,2)
genotypes <- matrix(c(2,0,2,0,4,
                      0,4,0,1,3,
                      1,1,2,2,2), nrow = 3, ncol = 5, byrow = TRUE)
depth <- matrix(sample(100, 15), nrow = 3, ncol = 5)
depth[genotypes == 0] <- 0
alnuc <- c("AA", "AG", "CA", "G", "T")
locTable <- data.frame(Chr = c("Chr01", "Chr03"),
                       Pos = c(101, 501),
                       Ref = c("AG", "G"),
                       stringsAsFactors = FALSE)
PrepVCFexport(genotypes, alleles2loc, depth, alnuc, locTable, 4, TRUE)
*/

