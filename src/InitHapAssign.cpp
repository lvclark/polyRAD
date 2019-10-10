#include <Rcpp.h>
using namespace Rcpp;

// Rcpp version of the InitHapAssign python function to take a matrix of NM
// values and assign tags to alignment locations.
// Returns a vector with locus index, rather than a list with allele index.

// [[Rcpp::export]]
IntegerVector InitHapAssign(IntegerMatrix NMmat){
  int nloc = NMmat.ncol();
  int nhap = NMmat.nrow();
  IntegerVector out(nhap);
  IntegerVector thisrow;
  int minNM;
  IntegerVector bestLoc;
  IntegerVector loci = seq(1, nloc);
  
  for(int h = 0; h < nhap; h++){
    thisrow = NMmat(h, _ );
    minNM = min(thisrow);
    bestLoc = loci[thisrow == minNM];
    if(bestLoc.size() == 1){
      out[h] = bestLoc[0];
    } else {
      out[h] = sample(bestLoc, 1)[0];
    }
  }
  
  return out;
}
