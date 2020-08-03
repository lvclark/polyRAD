#include <Rcpp.h>
using namespace Rcpp;

// Randomly genearte a genotype, given allele frequencies and inbreeding
// [[Rcpp::export]]
NumericVector sampleGenotype(NumericVector freq, double inbreeding, int ploidy) {
  int nal = freq.size();
  NumericVector geno(nal);
  double repprob;
  double rnd;
  int a;
  
  for(double k = 0; k < ploidy; k++){
    // get probability that this allele is the same as a previous one
    repprob = 1 - pow(1 - inbreeding, k);
    rnd = runif(1)[0];
    if(rnd <= repprob){
      a = sample(nal, 1, true, geno, false)[0];
    }
    else{
      a = sample(nal, 1, true, freq, false)[0];
    }
    geno[a] += 1;
  }
  
  return geno;
}
