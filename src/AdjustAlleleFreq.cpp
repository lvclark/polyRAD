#include <Rcpp.h>
using namespace Rcpp;

// Internal function for adjusting allele frequency at one taxon and locus
NumericVector AdjustOneFreq(NumericVector thesefreq,
                            double minfreq, double maxfreq, int thesenal){
  IntegerVector alseq = seq(0, thesenal - 1);
  NumericVector oldfreq(thesenal);
  LogicalVector toolow(thesenal);
  LogicalVector toohigh(thesenal);
  LogicalVector canadjustL(thesenal);
  IntegerVector canadjust;
  double adjust;
  
  thesefreq = thesefreq/sum(thesefreq); // make frequencies sum to one
  
  while(is_true(any(thesefreq < minfreq))){
    oldfreq = clone(thesefreq);
    toolow = thesefreq < minfreq;
    thesefreq[toolow] = minfreq;
    canadjustL = thesefreq > minfreq;
    canadjust = alseq[canadjustL];
    adjust = sum(thesefreq - oldfreq)/canadjust.size();
    for(int a = 0; a < canadjust.size(); a++){
      thesefreq[canadjust[a]] -= adjust;
    }
  }
  while(is_true(any(thesefreq > maxfreq))){
    oldfreq = clone(thesefreq);
    toohigh = thesefreq > maxfreq;
    thesefreq[toohigh] = maxfreq;
    canadjustL = thesefreq < maxfreq;
    canadjust = alseq[canadjustL];
    adjust = sum(oldfreq - thesefreq)/canadjust.size();
    for(int a = 0; a < canadjust.size(); a++){
      thesefreq[canadjust[a]] += adjust;
    }
  }
  
  return thesefreq;
}

// Internal function for AddAlleleFreqByTaxa
// [[Rcpp::export]]
NumericMatrix AdjustAlleleFreq(NumericMatrix predAl, IntegerVector alleles2loc,
                               double minfreq) {
  IntegerVector loci = unique(alleles2loc);
  int nloc = loci.size();   // number of loci
  int nal = predAl.ncol();  // number of alleles
  int ntax = predAl.nrow(); // number of taxa
  double maxfreq = 1.0 - minfreq;
  IntegerVector alleles = seq(0, nal - 1);
  IntegerVector thesecol;   // columns of matrix for a given locus
  NumericVector thesefreq;  // allele frequencies for a given locus and taxon
  int thesenal;             // number of alleles for a given locus
  
  for(int L = 0; L < nloc; L++){
    thesecol = alleles[alleles2loc == L + 1];
    thesenal = thesecol.size();
    thesefreq = NumericVector(thesenal);
    for(int taxa = 0; taxa < ntax; taxa++){
      // copy this row of matrix into temporary vector
      for(int a = 0; a < thesenal; a++){
        thesefreq[a] = predAl(taxa, thesecol[a]);
      }
      
      thesefreq = AdjustOneFreq(thesefreq, minfreq, maxfreq, thesenal);
      
      // copy vector back to matrix
      for(int a = 0; a < thesenal; a++){
        predAl(taxa, thesecol[a]) = thesefreq[a];
      }
    }
  }
  
  return predAl;
}
