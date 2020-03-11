#include <Rcpp.h>
using namespace Rcpp;

// Convert a set of haplotypes to a set of SNPs, with a position for each SNP
// and a matrix to indicate conversion of haplotypes to SNP alleles.
// [[Rcpp::export]]
List Hap2SNP(StringVector haps, std::string refhap, int pos) {
  int nhap = haps.size();
  unsigned int npos = haps[0].size();
  StringVector thesenuc(nhap);
  std::string thishap;
  std::string thisnuc;
  std::string thisref;
  bool isvar;
  List out(3);
  
  // check that all haplotypes are same length
  for(int i = 1; i < nhap; i++){
    if(haps[i].size() != npos){
      stop("All haplotypes must be same length.");
    }
  }
  
  // loop to search for variable sites
  for(unsigned int j = 0; j < npos; j++){
    isvar = false;
    for(int i = 0; i < nhap; i++){
      thishap = haps[i];
      thisnuc = thishap[j];
      thisref = refhap[j];
      thesenuc[i] = thisnuc;
      if(thisnuc != thisref){
        isvar = true;
      }
    }
    if(isvar){
      for(int i = 0; i < nhap; i++){
        Rcout << thesenuc[i]; //testing
      }
      // will need padding for insertions and deletions
    }
    if(refhap[j] != '.'){
      pos += 1;
      // possibly adjust for indels bigger than 1 nt
    }
  }
  
  return out;
}
