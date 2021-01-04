#include <Rcpp.h>
using namespace Rcpp;

// Convert genotypes to Structure format

// [[Rcpp::export]]
IntegerMatrix FormatStructure(IntegerMatrix genotypes, IntegerVector alleles2loc, int ploidy) {
  int nloc = max(alleles2loc);
  int nsam = genotypes.rows();
  int nal = alleles2loc.size();
  IntegerVector alleles = seq(0, nal - 1);
  IntegerVector thesecol;
  int thisnal;
  int curr_row = 0;
  int thiscopynum;
  IntegerMatrix out(nsam * ploidy, nloc);
  std::fill(out.begin(), out.end(), -9);
  
  for(int L = 0; L < nloc; L++){
    thesecol = alleles[alleles2loc == L + 1];
    thisnal = thesecol.size();
    curr_row = 0;
    for(int s = 0; s < nsam; s++){
      //if(curr_row > s * ploidy){
      //  Rcout << "Locus " << L << " sample " << s << "\n";
      //  stop("Ploidy doesn't match genotype matrix");
      //}
      curr_row = s * ploidy;
      for(int a = 0; a < thisnal; a ++){
        if(IntegerVector::is_na(genotypes(s, thesecol[a]))){
          continue;
        }
        thiscopynum = genotypes(s, thesecol[a]);
        for(int c = 0; c < thiscopynum; c++){
          out(curr_row, L) = a + 1;
          curr_row++;
        }
      }
    }
  }
  return out;
}
