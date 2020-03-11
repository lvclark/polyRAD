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
  
  int nvar = 0; // number of variable sites found
  IntegerVector outpos(npos);
  List outnuc(npos);
  List outmat(npos);
  IntegerMatrix thismat;
  
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
      outpos[nvar] = pos;
      // pad insertions or deletions
      if(any(thesenuc == "." | thesenuc == "-") && !any(outpos == pos - 1)){
        outpos[nvar] -= 1;
        thisref = refhap[j-1] + thisref;
        for(int i = 0; i < nhap; i++){
          thesenuc[i] = refhap[j-1] + thesenuc[i]
        }
        // add something here to add downstream sites for longer indels
      }
      
      nvar += 1;
    }
    if(refhap[j] != '.'){
      pos += 1;
      // possibly adjust for indels bigger than 1 nt
    }
  }
  
  List out = List::create(outpos, outnuc, outmat);
  return out;
}
