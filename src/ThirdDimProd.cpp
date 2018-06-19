#include <Rcpp.h>
using namespace Rcpp;

// Take a three-dimensional array (as a vector) and return a matrix
// corresponding to the first two dimensions, containing the product across
// the third dimension.  Used internally by AddGenotypePriorProb_LD.

// [[Rcpp::export]]
NumericMatrix ThirdDimProd(NumericVector probs, int ngen, int ntaxa) {
  NumericMatrix out(ngen, ntaxa);
  int copynum;
  int taxon;
  int probsize = probs.size();
  
  // Replace zeros in matrix with values from first linked alleles
  for(int i = 0; i < ngen * ntaxa; i++){
    copynum = i % ngen;
    taxon = i / ngen % ntaxa;
    
    out(copynum, taxon) = probs[i];
  }
  // Multiply by remaining alleles
  for(int i = ngen * ntaxa; i < probsize; i++){
    copynum = i % ngen;
    taxon = i / ngen % ntaxa;
    
    out(copynum, taxon) *= probs[i];
  }
  
  return out;
}
