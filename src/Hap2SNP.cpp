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
  
  IntegerVector starts(npos);
  IntegerVector ends(npos);
  int nsites = 0;
  
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
  
  // loop to determine sites to report
  unsigned int j = 0;
  while(j < npos){
    if(isvar[j]){
      if(isindel[j]){
        // add assert that j > 0?
        starts[nsites] = j - 1;
      } else {
        starts[nsites] = j;
      }
      while(j < npos - 1 && isindel[j + 1]){
        j++;
      }
      ends[nsites] = j;
      nsites++;
    }
    j++;
  }
  
  IntegerVector outpos(nsites);
  List outnuc(nsites);
  List outmat(nsites);
  IntegerMatrix thismat;
  int start;
  int end;
  
  // loop to determine which haplotypes have which alleles
  for(int n = 0; n < nsites; n++){
    start = starts[n];
    end = ends[n];
  }
  
  //List out = List::create(outpos, outnuc, outmat);
  List out = List::create(allpos, isvar, isindel, starts, ends);
  return out;
}

// Testing
/*** R
myref <- "ACGT.AAGCGCTT.AC"
myhaps <- c("ACGTTAAGCGCTT.AA", "ACGTTCAGCGCTT.AC", "ACGTTCAGCGCTG.AC",
            "ACGTTA--CGCTT.AC", "ACGTTAAGCGCTTCAC", "TCGTTAAGCGCTTCAC")
Hap2SNP(myhaps, myref, 201)
*/
