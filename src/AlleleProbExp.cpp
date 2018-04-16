#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

// Function to take a vector of probabilities of sampling a read of a given 
// allele for a given genotype (vector of all alleles in the dataset, for just
// one allele copy number), and a matrix of read depth individual x allele,
// and raise the sampling probabilities to the power of the read depth, as 
// part of calculating binomial probability.
NumericMatrix AlleleProbExp(IntegerMatrix depth, NumericVector alleleProb) {
  int alleles = alleleProb.size();
  int samples = depth.nrow();
  NumericMatrix out(samples, alleles);
  
  for(int s = 0; s < samples; ++s){
    for(int a = 0; a < alleles; ++a){
      out(s, a) = 1;
      for(int i = 0; i < depth(s, a); ++i){
        out(s, a) *= alleleProb[a];
      }
    }
  }
  
  return out;
}
