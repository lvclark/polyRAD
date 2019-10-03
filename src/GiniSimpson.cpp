#include <Rcpp.h>
using namespace Rcpp;

// Function to get the Gini-Simpson index, for a group of sequencing reads or
// anything else.  Calculated WITH replacement.
// [[Rcpp::export]]
double GiniSimpson(NumericVector counts) {
  double N = sum(counts);
  double out = 1.0;
  
  for(int i = 0; i < counts.size(); i++){
    out -= pow(counts[i] / N, 2);
  }
  
  return out;
}

// Function to get Hind/He on allelic read depth as encoded in polyRAD.
// Returns a vector of values, individual x locus.
// He is calculated on the fly for natural populations and diversity panels,
// or input for mapping populations.
// [[Rcpp::export]]
NumericMatrix HindHeMat(IntegerMatrix alleleDepth, NumericMatrix depthRatio,
                        IntegerVector alleles2loc, int nLoci, NumericVector He){
  int nTaxa = alleleDepth.nrow();
  IntegerVector alleles = seq(0, alleles2loc.size() - 1);
  IntegerVector thesecol;
  NumericVector thesecounts;
  int thesenal;
  int thisdepth;
  NumericVector thesefreq;
  double thisHind;
  double thisHe;
  bool calcHe = He.size() == 0;
  NumericMatrix HindHe(nTaxa, nLoci);
  
  for(int L = 0; L < nLoci; L++){
    thesecol = alleles[alleles2loc == L + 1];
    thesenal = thesecol.size();
    thesecounts = NumericVector(thesenal);
    if(calcHe){
      thesefreq = NumericVector(thesenal);
      // fill in allele frequencies from depth ratios
      for(int a = 0; a < thesenal; a++){
        thesefreq[a] = mean(na_omit(depthRatio( _ , thesecol[a])));
      }
      thisHe = GiniSimpson(thesefreq);
    } else {
      thisHe = He[L];
    }
    
    for(int t = 0; t < nTaxa; t++){
      // For this taxon and locus, copy to temporary vector
      // and add to vector for total across taxa.
      for(int a = 0; a < thesenal; a++){
        thesecounts[a] = alleleDepth(t, thesecol[a]);
      }
      thisdepth = sum(thesecounts);
      if(thisdepth < 2){ // can't calculate for depth below 2
        HindHe(t, L) = R_NaN;
      } else {
        thisHind = GiniSimpson(thesecounts); // H for t x L
        HindHe(t, L) = thisHind / thisHe * thisdepth / (thisdepth - 1);
      }
    }
  }
  
  return HindHe;
}

// Function to get a per-locus estimate of probability of sampling two different
// alleles (without replacement) from a parental genotype.
// [[Rcpp::export]]
NumericVector HoOneParent(IntegerVector genotypes, IntegerVector alleles2loc,
                          IntegerVector keeploc, double ploidy){
  int nloc = keeploc.size();
  int L;
  IntegerVector thisgen;
  NumericVector out(nloc, 1.0);
  
  for(int i = 0; i < nloc; i++){
    L = keeploc[i];
    thisgen = genotypes[alleles2loc == L];
    for(int a = 0; a < thisgen.size(); a++){
      out[i] -= (thisgen[a] / ploidy) * ((thisgen[a] - 1)/(ploidy - 1));
    }
  }
  
  return out;
}

// Function to get a per-locus estimate of the probability of sampling two
// different alleles if one is selected from each of two parental genotypes.
// [[Rcpp::export]]
NumericVector HoTwoParents(IntegerVector genotypes1, IntegerVector genotypes2,
                           IntegerVector alleles2loc, IntegerVector keeploc,
                           double ploidy){
  int nloc = keeploc.size();
  int L;
  IntegerVector thisgen1;
  IntegerVector thisgen2;
  NumericVector out(nloc, 1.0);
  
  for(int i = 0; i < nloc; i++){
    L = keeploc[i];
    thisgen1 = genotypes1[alleles2loc == L];
    thisgen2 = genotypes2[alleles2loc == L];
    for(int a = 0; a < thisgen1.size(); a++){
      out[i] -= (thisgen1[a] / ploidy) * (thisgen2[a]/ploidy);
    }
  }
  
  return out;
}
