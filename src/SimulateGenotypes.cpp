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

// Simulate a genotype matrix in a mapping population
// progGeno and genoProbs come from .buildProgProb
// [[Rcpp::export]]
NumericMatrix simGenoMapping(NumericVector donorGeno, NumericVector recurGeno, NumericMatrix progGeno,
                             NumericVector genoProbs, IntegerVector alleles2loc, int nsam, int ploidy){
  int nal = alleles2loc.size();
  IntegerVector alleles = seq(0, nal - 1);
  IntegerVector thesecol;
  int thisnal;
  NumericVector thisgeno;
  NumericVector thisDonorGeno;
  NumericVector thisRecurGeno;
  NumericMatrix out(nsam, nal);
  int maxL = max(alleles2loc);
  int ngeno = genoProbs.size();
  int genoindex;
  int thisal1;
  int thisal2 = 0;
  
  for(int L = 1; L <= maxL; L++){
    thesecol = alleles[alleles2loc == L];
    thisnal = thesecol.size();
    thisDonorGeno = donorGeno[thesecol];
    thisRecurGeno = recurGeno[thesecol];
    // Cumulative sums of genotypes for allele lookup
    // For some reason I get an error if declaring in header
    NumericVector thisDonorCum = cumsum(thisDonorGeno);
    NumericVector thisRecurCum = cumsum(thisRecurGeno);
    for(int s = 0; s < nsam; s++){
      // Determine which genotype this sample has
      // NEED TO FIX; always gives 0
      genoindex = sample(ngeno, 1, false, genoProbs, false)[0];
      // For each allele in the genotype, determine its index in the output matrix
      for(int a = 0; a < ploidy; a++){
        // Number from 1 to 2*ploidy, indicating parent and copy
        thisal1 = progGeno(genoindex, a);
        // Higher numbers are for donor parent
        if(thisal1 > ploidy){
          thisal1 -= ploidy;
          for(int b = 0; b < thisnal; b++){
            if(thisDonorCum[b] >= thisal1){
              thisal2 = thesecol[b];
              break;
            }
          }
        }
        // Lower numbers are for recurrent parent
        else {
          for(int b = 0; b < thisnal; b++){
            if(thisRecurCum[b] >= thisal1){
              thisal2 = thesecol[b];
              break;
            }
          }
        }
        out(s, thisal2) += 1;
      }
    }
  }
  
  return out;
}

// simulate an allele depth matrix, given locus depth and genotypes
// [[Rcpp::export]]
IntegerMatrix simAD(IntegerMatrix locDepth, NumericMatrix genotypes,
                    IntegerVector alleles2loc, double overdispersion,
                    double contamRate, NumericVector alleleFreq,
                    double errorRate){
  int nsam = genotypes.rows();
  int nal = alleles2loc.size();
  int nloc = locDepth.cols();
  IntegerVector alleles = seq(0, nal - 1);
  IntegerVector thesecol;
  int thisnal;
  NumericVector thesefreq;
  NumericVector thesecontam;
  NumericVector thisgeno;
  NumericVector theseprobs;
  NumericVector newprobs;
  double evendist;
  double toadd;
  IntegerVector thesedepths;
  IntegerMatrix out(nsam, nal);
  
  for(int L = 1; L <= nloc; L++){
    thesecol = alleles[alleles2loc == L];
    thisnal = thesecol.size();
    thesefreq = alleleFreq[thesecol];
    thesecontam = thesefreq * contamRate;
    thisgeno = NumericVector(thisnal);
    for(int s = 0; s < nsam; s++){
      // Retrieve genotype
      for(int a = 0; a < thisnal; a++){
        thisgeno[a] = genotypes(s, thesecol[a]);
      }
      // Get allele sampling probabilities
      theseprobs = thisgeno / sum(thisgeno) * (1 - contamRate) + thesecontam;
      newprobs = theseprobs * (1 - errorRate);
      if(errorRate > 0){
        // Probability of sampling some other allele with error
        evendist = 1.0 / (thisnal - 1.0) * errorRate;
        for(int a = 0; a < thisnal; a++){
          // Prob. that a was the true allele but b was detected due to error
          toadd = evendist * theseprobs[a];
          for(int b = 0; b < thisnal; b++){
            if(b != a){
              newprobs[b] += toadd;
            }
          }
        }
      }
      // Simulate depths and add to matrix
      thesedepths = sampleReads(newprobs, locDepth(s, L - 1), overdispersion);
      for(int a = 0; a < thisnal; a++){
        out(s, thesecol[a]) = thesedepths[a];
      }
    }
  }
  
  return out;
}
