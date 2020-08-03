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

// Randomly generate reads based on a genotype and overdispersion
// [[Rcpp::export]]
IntegerVector sampleReads(NumericVector geno, int nreads, double overdispersion){
  int nal = geno.size();
  NumericVector initprobs = geno / sum(geno);
  NumericVector alpha = initprobs * overdispersion;
  NumericVector newprobs(nal);
  IntegerVector out(nal);
  
  // Get read probabilities from gamma distribution
  for(int a = 0; a < nal; a++){
    newprobs[a] = Rcpp::rgamma(1, alpha[a], 1.0)[0];
  }
  
  // Sample the reads
  IntegerVector readassign = sample(nal, nreads, true, newprobs, false);
  for(int i = 0; i < nreads; i++){
    out[readassign[i]] += 1;
  }
  
  return out;
}
