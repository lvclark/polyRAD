#include <Rcpp.h>
using namespace Rcpp;

// Convert a set of haplotypes to a set of SNPs, with a position for each SNP
// and a matrix to indicate conversion of haplotypes to SNP alleles.
// [[Rcpp::export]]
List Hap2SNP(std::vector< std::string > haps, std::string refhap, int pos) {
  int nhap = haps.size();
  int npos = haps[0].size();
  std::vector< std::string > thesenuc;
  std::string thisnuc;
  bool isvar;
  
  // Need to deal with insertions making alleleNucleotides shorter than
  // the reference haplotype.  Possibly will have to deal with that before
  // this function, in order to get the right size for the reference haplotype.
  // Search for * to locate insertion points.  Then there is the issue I knew
  // I would run into, with all haplotypes having the same insertion with
  // respect to the reference.  Maybe something in locTable could indicate these.
  // Should MakeAlleleStrings in isoloci_fun.py also return the reference?
  
  // check that all haplotypes are same length
  for(int i = 1; i < nhap; i++){
    if(haps[i].size() != npos){
      stop("All haplotypes must be same length.");
    }
  }
  
  // loop to search for variable sites
  for(int j = 0; j < npos; j++){
    thesenuc = (refhap[j]);
    isvar = false;
    for(int i = 0; i < nhap; i++){
      thisnuc = haps[i][j];
      thesenuc.emplace_back(thisnuc);
      if(thisnuc != refhap[j]){
        isvar = true;
      }
    }
    if(isvar){
      
    }
  }
}
