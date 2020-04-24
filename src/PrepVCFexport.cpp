#include <Rcpp.h>
#include <string>
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

// Function to take genotype calls and slots from a RADdata object and prepare
// data for export to VCF.

// [[Rcpp::export]]
List PrepVCFexport(IntegerMatrix genotypes, IntegerVector alleles2loc,
                   IntegerMatrix alleleDepth, StringVector alleleNucleotides,
                   DataFrame locTable, int ploidy) {
  List out;
  return out;
}



