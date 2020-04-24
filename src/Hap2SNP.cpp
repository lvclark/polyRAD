#include <Rcpp.h>
#include <string>
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
      if(thisnuc != thisref && 
         (thisnuc == "A" || thisnuc == "C" || thisnuc == "G" || thisnuc == "T" ||
         thisnuc == "." || thisnuc == "-")){
        // Nucleotide differs from reference AND
        // is not ambiguous due to haplotype merging.
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
  
  // loop to determine sites to report.
  // Use VCF format for indels (includes previous nucleotide)
  unsigned int j = 0;
  while(j < npos){
    if(isvar[j]){
      if(isindel[j]){
        // add assert that j > 0? (does give error on its own)
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
  int start;
  int vlen;
  StringVector uniquenuc;
  bool nucfound;
  int nal;
  String thisal;
  
  // loop to determine which haplotypes have which alleles
  for(int n = 0; n < nsites; n++){
    start = starts[n];
    vlen = ends[n] - start + 1;
    outpos[n] = allpos[start];
    thisref = refhap.substr(start, vlen);
    uniquenuc = StringVector::create(thisref);
    for(int i = 0; i < nhap; i++){
      thishap = haps[i];
      thisnuc = thishap.substr(start, vlen);
      thesenuc[i] = thisnuc;
      // See if this allele has been found already.
      // Add to unique list, maintaining order, esp. reference first.
      nucfound = false;
      for(unsigned int u = 0; u < uniquenuc.size(); u++){
        if(uniquenuc[u] == thisnuc){
          nucfound = true;
          break;
        }
      }
      if(!nucfound){
        uniquenuc.push_back(thisnuc);
      }
    }
    
    // Conversion matrix
    nal = uniquenuc.size();
    IntegerMatrix thismat(nhap, nal);
    for(int i = 0; i < nhap; i++){
      for(int u = 0; u < nal; u++){
        if(thesenuc[i] == uniquenuc[u]){
          thismat(i, u) = 1;
          break;
        }
      }
    }
    
    // Remove insertion and deletion characters
    for(int u = 0; u < nal; u++){
      thisal = uniquenuc[u];
      thisal.replace_all(".", "");
      thisal.replace_all("-", "");
      uniquenuc[u] = thisal;
    }

    outnuc[n] = uniquenuc;
    outmat[n] = thismat;
  }
  
  List out = List::create(outpos, outnuc, outmat);
  //List out = List::create(allpos, isvar, isindel, starts, ends, outpos, outnuc);
  return out;
}

// Testing
/*** R
myref <- "ACGT.AAGCGCTT.AC"
myhaps <- c("ACGTTAAGCGCTT.AA", "ACGTTCAGCGCTT.AC", "ACGTTCAGCGCTG.AC",
            "ACGTTA--CGCYT.AC", "ACGTTAAGCGCTTCAC", "TCGTTAAGCGCTTCAC")
Hap2SNP(myhaps, myref, 201)
*/
