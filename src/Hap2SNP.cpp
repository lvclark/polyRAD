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
  LogicalVector isvar(npos);
  LogicalVector isindel(npos);
  IntegerVector allpos(npos);
  
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
    thisref = refhap[j];
    if(thisref == "."){
      isindel[j] = true;
    }
    for(int i = 0; i < nhap; i++){
      thishap = haps[i];
      thisnuc = thishap[j];
      if(thisnuc != thisref){
        isvar[j] = true;
      }
      if(thisnuc == "." || thisnuc == "-"){
        isindel[j] = true;
      }
    }
    if(j == 0){
      allpos[j] = pos;
    } else {
      if(refhap[j] != '.'){
        allpos[j] = allpos[j - 1] + 1;
      } else {
        allpos[j] = allpos[j - 1];
      }
    }
  }
  
  //List out = List::create(outpos, outnuc, outmat);
  List out = List::create(allpos, isvar, isindel);
  return out;
}
