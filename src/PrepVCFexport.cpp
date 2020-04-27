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

// Function to take genotype calls and slots from a RADdata object and prepare
// data for export to VCF.

// [[Rcpp::export]]
List PrepVCFexport(IntegerMatrix genotypes, IntegerVector alleles2loc,
                   IntegerMatrix alleleDepth, StringVector alleleNucleotides,
                   DataFrame locTable, int ploidy, bool asSNPs) {
  int nloc = locTable.nrows();
  int nsam = genotypes.nrow();
  List alleleLookup = AlleleIndex(alleles2loc, nloc);
  IntegerVector thesecol;
  int thisnal;
  StringVector thesehap;
  IntegerVector posvect = locTable["Pos"];
  StringVector refvect = locTable["Ref"];
  IntegerMatrix thesegeno;
  IntegerMatrix thesedepths;
  int pos;
  std::string ref;
  List hapconv;
  List outpos(nloc);
  List outal(nloc);
  List outgen(nloc);
  List outdepth(nloc);
  
  for(int L = 0; L < nloc; L++){
    // Subset data for this locus
    thesecol = alleleLookup[L];
    thisnal = thesecol.size();
    thesehap = alleleNucleotides[thesecol];
    
    thesegeno = SubsetMatrixCol(genotypes, thesecol);
    thesedepths = SubsetMatrixCol(alleleDepth, thesecol);
    
    pos = posvect[L];
    ref = refvect[L];
    
    // Convert haplotypes to SNPs
    if(asSNPs){
      hapconv = Hap2SNP(thesehap, ref, pos);
      
      // Insert code to convert genotypes and depths
      // Make separate function to put depths in to AD format for VariantAnnotation
    } else {
      hapconv = Hap2Hap(thesehap, ref, pos);

      // Insert code to add genotypes and depths

    }
    
    outpos[L] = hapconv[0];
    outal[L] = hapconv[1];
    
  }
  List out;
  return out;
}



