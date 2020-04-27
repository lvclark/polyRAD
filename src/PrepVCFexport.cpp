#include <Rcpp.h>
#include <string>
using namespace Rcpp;

// Prepare GT strings for the VCF from a matrix of allele copy number.

// [[Rcpp::export]]
StringVector MakeGTstrings(IntegerMatrix genotypes, int ploidy) {
  int nal = genotypes.ncol();
  int nsam = genotypes.nrow();
  IntegerVector thisgen(nal);
  StringVector out(nsam);
  int cp;
  String thisstring;
  
  for(int s = 0; s < nsam; s++){
    thisstring = "";
    thisgen = genotypes(s, _);
    if(any(is_na(thisgen))){
      for(int k = 0; k < ploidy; k++){
        thisstring += "./";
      }
    } else {
      for(int a = 0; a < nal; a++){
        cp = thisgen[a];
        for(int c = 0; c < cp; c++){
          thisstring += a;
          thisstring += "/";
        }
      }
    }

    thisstring.replace_last("/", "");
    out[s] = thisstring;
  }
  return out;
}

// Function to take genotype calls and slots from a RADdata object and prepare
// data for export to VCF.

// [[Rcpp::export]]
List PrepVCFexport(IntegerMatrix genotypes, IntegerVector alleles2loc,
                   IntegerMatrix alleleDepth, StringVector alleleNucleotides,
                   DataFrame locTable, int ploidy, bool asSNPs) {
  int nloc = locTable.nrows();
  int nsam = genotypes.nrow();
  IntegerVector alleles = seq(0, alleles2loc.size() - 1);
  IntegerVector thesecol;
  int thisnal;
  StringVector thesehap;
  int pos;
  std::string ref;
  List outpos(nloc);
  List outal(nloc);
  List outgen(nloc);
  List outdepth(nloc);
  
  for(int L = 0; L < nloc; L++){
    // Subset data for this locus
    thesecol = alleles[alleles2loc == L + 1];
    thisnal = thesecol.size();
    IntegerMatrix thesegeno(nsam, thisnal);
    IntegerMatrix thesedepths(nsam, thisnal);
    thesehap = alleleNucleotides[thesecol];
    
    for(int a = 0; a < thisnal; a++){
      for(int s = 0; s < nsam; s++){
        thesegeno(s, a) = genotypes(s, thesecol[a]);
        thesedepths(s, a) = alleleDepth(s, thesecol[a]);
      }
    }
    
    pos = locTable["Pos"][L];
    ref = locTable["Ref"][L]
    
    // Convert haplotypes to SNPs
    if(asSNPs){
      hapconv = Hap2SNP(thesehap, ref, pos);
      outpos[L] = hapconv[0];
      outal[L] = hapconv[1];
      // Insert code to convert genotypes and depths
      // Make separate function to put depths in to AD format for VariantAnnotation
    } else {
      outpos[L] = IntegerVector::create(pos);
      outal[L] = List::create(alleleNucleotides);
      // Insert code to add genotypes and depths
    }
    
  }
  List out;
  return out;
}



