#include <Rcpp.h>
using namespace Rcpp;

// Function to take a 3D matrix of genotype likelihoods or probabilities
// (as a vector) and output a matrix of the best genotypes.
// [[Rcpp::export]]
IntegerMatrix BestGenos(NumericVector probs, int ploidy, int ntaxa, int nalleles) {
  IntegerMatrix bestgenos(ntaxa, nalleles);
  int nprobs = probs.size();
  int ngen = ploidy + 1;
  int ngentimestaxa = ngen * ntaxa;
  int bestgen = NA_INTEGER;
  float bestprob = 0;
  int copynum;
  int taxon;
  int allele;
  
  for(int i = 0; i < nprobs; i++){
    copynum = i % ngen;
    taxon = i / ngen % ntaxa;
    allele = i / ngentimestaxa;
    
    if(copynum == 0 || probs[i] > bestprob){
      bestgen = copynum;
      bestprob = probs[i];
      if(NumericVector::is_na(probs[i])){
        bestgen = NA_INTEGER;
        bestprob = 0;
      }
    }
    if(copynum == ploidy){
      bestgenos(taxon, allele) = bestgen;
    }
  }
  
  return bestgenos;
}

// Function to find best ploidies from ploidyChiSq slot
// [[Rcpp::export]]
IntegerVector BestPloidies(NumericMatrix chisq) {
  int nalleles = chisq.ncol();
  int npld = chisq.nrow();
  IntegerVector bestploidies(nalleles);
  int bestpld;
  float bestchisq;
  
  for(int a = 0; a < nalleles; a++){
    bestpld = 0;
    bestchisq = chisq(0,a);
    for(int pld = 0; pld < npld; pld++){
      if(chisq(pld, a) < bestchisq || (NumericVector::is_na(bestchisq) && !NumericVector::is_na(chisq(pld, a)))){
        bestpld = pld;
        bestchisq = chisq(pld, a);
      }
    }
    if(NumericVector::is_na(bestchisq)){
      bestpld = -1; // Zero in the output will indicate they were all missing.
    }
    bestploidies[a] = bestpld + 1;
  }
  
  return bestploidies;
}
