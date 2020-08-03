#include <Rcpp.h>
using namespace Rcpp;

// Randomly genearte a genotype, given allele frequencies and inbreeding
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

// Simulate a genotype matrix
// [[Rcpp::export]]
NumericMatrix simGeno(NumericVector alleleFreq, IntegerVector alleles2loc, int nsam, double inbreeding, int ploidy){
  int nal = alleles2loc.size();
  IntegerVector alleles = seq(0, nal - 1);
  IntegerVector thesecol;
  int thisnal;
  NumericVector thesefreq;
  NumericVector thisgeno;
  NumericMatrix out(nsam, nal);
  int maxL = max(alleles2loc);
  
  for(int L = 1; L <= maxL; L++){
    thesecol = alleles[alleles2loc == L];
    thesefreq = alleleFreq[thesecol];
    thisnal = thesecol.size();
    for(int s = 0; s < nsam; s++){
      thisgeno = sampleGenotype(thesefreq, inbreeding, ploidy);
      for(int a = 0; a < thisnal; a++){
        out(s, thesecol[a]) = thisgeno[a];
      }
    }
  }
  
  return out;
}

// simulate an allele depth matrix, given locus depth and genotypes
// [[Rcpp::export]]
IntegerMatrix simAD(IntegerMatrix locDepth, NumericMatrix genotypes, IntegerVector alleles2loc, double overdispersion){
  int nsam = genotypes.rows();
  int nal = alleles2loc.size();
  int nloc = locDepth.cols();
  IntegerVector alleles = seq(0, nal - 1);
  IntegerVector thesecol;
  int thisnal;
  NumericVector thisgeno;
  IntegerVector thesedepths;
  IntegerMatrix out(nsam, nal);
  
  for(int L = 1; L <= nloc; L++){
    thesecol = alleles[alleles2loc == L];
    thisnal = thesecol.size();
    thisgeno = NumericVector(thisnal);
    for(int s = 0; s < nsam; s++){
      // Retrieve genotype
      for(int a = 0; a < thisnal; a++){
        thisgeno[a] = genotypes(s, thesecol[a]);
      }
      // Simulate depths and add to matrix
      thesedepths = sampleReads(thisgeno, locDepth(s, L - 1), overdispersion);
      for(int a = 0; a < thisnal; a++){
        out(s, thesecol[a]) = thesedepths[a];
      }
    }
  }
  
  return out;
}
