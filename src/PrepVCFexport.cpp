#include <Rcpp.h>
#include <string>
using namespace Rcpp;

// Remove indel characters from allele strings
StringVector RemoveIndelChar(StringVector uniquenuc, int nal){
  String thisal;
  
  for(int u = 0; u < nal; u++){
    thisal = uniquenuc[u];
    thisal.replace_all(".", "");
    thisal.replace_all("-", "");
    uniquenuc[u] = thisal;
  }
  return uniquenuc;
}

// Convert a set of haplotypes to a set of SNPs, with a position for each SNP
// and a matrix to indicate conversion of haplotypes to SNP alleles.
// [[Rcpp::export]]
List Hap2SNP(StringVector haps, std::string refhap, int pos) {
  int nhap = haps.size();
  unsigned int npos = haps[0].size();
  StringVector thesenuc(nhap);
  std::string thishap;
  std::string thisnuc;
  std::string thisref;
  LogicalVector isvar(npos);
  LogicalVector isindel(npos);
  IntegerVector allpos(npos);
  
  IntegerVector starts(npos);
  IntegerVector ends(npos);
  int nsites = 0;
  
  // check that all haplotypes are same length
  for(int i = 1; i < nhap; i++){
    if(haps[i].size() != npos){
      stop("All haplotypes must be same length.");
    }
  }
  
  // loop to search for variable sites
  for(unsigned int j = 0; j < npos; j++){
    thisref = refhap[j];
    if(thisref == "."){
      isindel[j] = true;
    }
    for(int i = 0; i < nhap; i++){
      thishap = haps[i];
      thisnuc = thishap[j];
      if(thisnuc != thisref && 
         (thisnuc == "A" || thisnuc == "C" || thisnuc == "G" || thisnuc == "T" ||
         thisnuc == "." || thisnuc == "-")){
        // Nucleotide differs from reference AND
        // is not ambiguous due to haplotype merging.
        isvar[j] = true;
      }
      if(thisnuc == "." || thisnuc == "-"){
        isindel[j] = true;
      }
    }
    if(j == 0){
      allpos[j] = pos;
    } else {
      if(refhap[j] != '.'){
        allpos[j] = allpos[j - 1] + 1;
      } else {
        allpos[j] = allpos[j - 1];
      }
    }
  }
  
  // loop to determine sites to report.
  // Use VCF format for indels (includes previous nucleotide)
  unsigned int j = 0;
  while(j < npos){
    if(isvar[j]){
      if(isindel[j]){
        // add assert that j > 0? (does give error on its own)
        starts[nsites] = j - 1;
      } else {
        starts[nsites] = j;
      }
      while(j < npos - 1 && isindel[j + 1]){
        j++;
      }
      ends[nsites] = j;
      nsites++;
    }
    j++;
  }
  
  IntegerVector outpos(nsites);
  List outnuc(nsites);
  List outmat(nsites);
  int start;
  int vlen;
  StringVector uniquenuc;
  bool nucfound;
  int nal;
  
  // loop to determine which haplotypes have which alleles
  for(int n = 0; n < nsites; n++){
    start = starts[n];
    vlen = ends[n] - start + 1;
    outpos[n] = allpos[start];
    thisref = refhap.substr(start, vlen);
    uniquenuc = StringVector::create(thisref);
    for(int i = 0; i < nhap; i++){
      thishap = haps[i];
      thisnuc = thishap.substr(start, vlen);
      thesenuc[i] = thisnuc;
      // See if this allele has been found already.
      // Add to unique list, maintaining order, esp. reference first.
      nucfound = false;
      for(unsigned int u = 0; u < uniquenuc.size(); u++){
        if(uniquenuc[u] == thisnuc){
          nucfound = true;
          break;
        }
      }
      if(!nucfound){
        uniquenuc.push_back(thisnuc);
      }
    }
    
    // Conversion matrix
    nal = uniquenuc.size();
    IntegerMatrix thismat(nhap, nal);
    for(int i = 0; i < nhap; i++){
      for(int u = 0; u < nal; u++){
        if(thesenuc[i] == uniquenuc[u]){
          thismat(i, u) = 1;
          break;
        }
      }
    }
    
    // Remove insertion and deletion characters
    uniquenuc = RemoveIndelChar(uniquenuc, nal);
    
    outnuc[n] = uniquenuc;
    outmat[n] = thismat;
  }
  
  List out = List::create(outpos, outnuc, outmat);
  //List out = List::create(allpos, isvar, isindel, starts, ends, outpos, outnuc);
  return out;
}

// Dummy fn for Hap2SNP that returns the same info but for haplotypes
// [[Rcpp::export]]
List Hap2Hap(StringVector haps, std::string refhap, int pos){
  int nhap = haps.size();
  int nhapOut = nhap;
  IntegerVector outpos = IntegerVector::create(pos);
  List outnuc(1);
  List outmat(1);
  List out;
  int refind = -1;
  
  // Is the reference among the alleles?
  for(int h = 0; h < nhap; h++){
    if(haps[h] == refhap){
      refind = h;
      break;
    }
  }
  
  // Arrange the sequences and conversion matrix, putting reference first
  int outind = 1;
  if(refind == -1){
    nhapOut++;
  }
  StringVector uniquenuc(nhapOut);
  uniquenuc[0] = refhap;
  IntegerMatrix thismat(nhap, nhapOut);
  for(int h = 0; h < nhap; h++){
    if(h == refind){
      thismat(h, 0) = 1;
    } else {
      thismat(h, outind) = 1;
      uniquenuc[outind] = haps[h];
      outind++;
    }
  }
  
  // Remove insertion and deletion characters
  uniquenuc = RemoveIndelChar(uniquenuc, nhapOut);
  
  outnuc[0] = uniquenuc;
  outmat[0] = thismat;
  out = List::create(outpos, outnuc, outmat);
  return out;
}

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

// Make a faster lookup for allele indices, for when loci will be looped through
// multiple times.  For internal use; indices start at zero in output.
// alleles2loc is assumed to be from RADdata and have indices starting at one.
List AlleleIndex(IntegerVector alleles2loc, int nloc){
  IntegerVector alleles = seq(0, alleles2loc.size() - 1);
  IntegerVector thesecol;
  List out(nloc);
  
  for(int L = 0; L < nloc; L++){
    thesecol = alleles[alleles2loc == L + 1];
    out[L] = thesecol;
  }
  return out;
}

// Subset a matrix using a vector of column indices
IntegerMatrix SubsetMatrixCol(IntegerMatrix mat, IntegerVector cols){
  int nr = mat.nrow();
  int co = cols.size();
  IntegerMatrix out(nr, co);
  
  for(int r = 0; r < nr; r++){
    for(int c = 0; c < co; c++){
      out(r, c) = mat(r, cols[c]);
    }
  }
  return out;
}

// Do the matrix multiplication with the conversion matrices, to avoid needing
// RcppArmadillo.  Assumes all values in the conversion matrix are 0 and 1, and
// the sum of each row is 1.
IntegerMatrix ConvMatMult(IntegerMatrix origmat, IntegerMatrix convmat){
  int nsam = origmat.nrow();
  int nhap = origmat.ncol();
  int nal = convmat.ncol();
  IntegerMatrix out(nsam, nal);
  
  for(int i = 0; i < nhap; i++){
    for(int j = 0; j < nal; j++){
      if(convmat(i, j) == 1){
        for(int s = 0; s < nsam; s++){
          out(s, j) += origmat(s, i);
        }
        break;
      }
    }
  }
  return out;
}

// Format allelic read depth for the AD slot in a VCF object.
List FormatAD(IntegerMatrix depthmat){
  int nsam = depthmat.nrow();
  IntegerVector thisdep;
  List out(nsam);
  
  for(int i = 0; i < nsam; i++){
    thisdep = depthmat(i, _ );
    out[i] = clone(thisdep);
  }
  return out;
}

// Function to take genotype calls and slots from a RADdata object and prepare
// data for export to VCF.

// [[Rcpp::export]]
List PrepVCFexport(IntegerMatrix genotypes, IntegerVector alleles2loc,
                   IntegerMatrix alleleDepth, StringVector alleleNucleotides,
                   DataFrame locTable, IntegerVector ploidy, bool asSNPs) {
  int nloc = locTable.nrows();
  int nsam = genotypes.nrow();
  List alleleLookup = AlleleIndex(alleles2loc, nloc);
  IntegerVector thesecol;
  StringVector thesehap;
  IntegerVector posvect = locTable["Pos"];
  StringVector refvect = locTable["Ref"];
  IntegerMatrix thesegeno;
  IntegerMatrix thesedepths;
  IntegerMatrix thesegeno1;
  IntegerMatrix thesedepths1;
  StringVector theseGT;
  List theseAD;
  int pos;
  int nsubloc;
  std::string ref;
  List hapconv(nloc);
  List thishapconv;
  IntegerVector sitesperloc(nloc);
  List thesemats;
  IntegerVector thesepos;
  List theseal;
  CharacterVector samplenames = rownames(genotypes);
  int pld;
  
  // Convert haplotypes to SNPs and tally variants
  for(int L = 0; L < nloc; L++){
    thesecol = alleleLookup[L];
    thesehap = alleleNucleotides[thesecol];
    pos = posvect[L];
    ref = refvect[L];
    
    if(asSNPs){
      thishapconv = Hap2SNP(thesehap, ref, pos);
    } else {
      thishapconv = Hap2Hap(thesehap, ref, pos);
    }
    hapconv[L] = clone(thishapconv);
    thesepos = thishapconv[0];
    sitesperloc[L] = thesepos.size();
  }
  
  // Convert data for each SNP/indel/haplotype
  int nsites = sum(sitesperloc);
  int currsite = 0;
  StringVector outREF(nsites);
  List outALT(nsites);
  IntegerVector outPos(nsites);
  StringMatrix outGT(nsam, nsites);
  rownames(outGT) = samplenames;
  List outAD(nsam * nsites);
  outAD.attr("dim") = Dimension(nsites, nsam);
  colnames(outAD) = samplenames;
  IntegerVector outLookup(nsites);
  
  for(int L = 0; L < nloc; L++){
    // Subset data for this locus
    thesecol = alleleLookup[L];
    thesegeno = SubsetMatrixCol(genotypes, thesecol);
    thesedepths = SubsetMatrixCol(alleleDepth, thesecol);
    thishapconv = hapconv[L];
    thesemats = thishapconv[2];
    thesepos = thishapconv[0];
    theseal = thishapconv[1];
    pld = ploidy[L];
    
    nsubloc = sitesperloc[L];
    for(int i = 0; i < nsubloc; i++){
      // Format genotypes and depths
      IntegerMatrix thismat = thesemats(i);
      thesegeno1 = ConvMatMult(thesegeno, thismat);
      thesedepths1 = ConvMatMult(thesedepths, thismat);
      theseGT = MakeGTstrings(thesegeno1, pld);
      theseAD = FormatAD(thesedepths1);
      outGT( _ , currsite) = theseGT;
      for(int s = 0; s < nsam; s++){
        IntegerVector theseAD1 = theseAD[s];
        outAD(currsite, s) = clone(theseAD1);
      }
      
      // Fill other data
      outPos[currsite] = thesepos[i];
      StringVector theseSNPal = theseal(i);
      outREF[currsite] = theseSNPal[0];
      theseSNPal.erase(0);
      outALT[currsite] = theseSNPal;
      outLookup[currsite] = L + 1;
      
      currsite++;
    }
  }
  List out = List::create(_["POS"] = outPos, _["REF"] = outREF, _["ALT"] = outALT,
                          _["GT"] = transpose(outGT), _["AD"] = outAD,
                          _["Lookup"] = outLookup);
  return out;
}

// Testing
/*** R
myref <- "ACGT.AAGCGCTT.AC"
myhaps <- c("ACGTTAAGCGCTT.AA", "ACGTTCAGCGCTT.AC", "ACGTTCAGCGCTG.AC",
            "ACGTTA--CGCYT.AC", "ACGTTAAGCGCTTCAC", "TCGTTAAGCGCTTCAC")
Hap2SNP(myhaps, myref, 201)
Hap2Hap(myhaps, myref, 201)

alleles2loc <- c(1,1,1,2,2)
genotypes <- matrix(c(2,0,2,0,4,
                      0,4,0,1,3,
                      1,1,2,2,2,
                      4,0,0,4,0), nrow = 4, ncol = 5, byrow = TRUE)
depth <- matrix(sample(100, 20), nrow = 4, ncol = 5)
depth[genotypes == 0] <- 0
alnuc <- c("AA", "AG", "CA", "G", "T")
rownames(depth) <- rownames(genotypes) <- paste0("Sam", 1:4)
colnames(depth) <- colnames(genotypes) <- 
  paste(c("Loc1", "Loc1", "Loc1", "Loc2", "Loc2"), alnuc, sep = "_")
locTable <- data.frame(Chr = c("Chr01", "Chr03"),
                       Pos = c(101, 501),
                       Ref = c("AG", "G"),
                       stringsAsFactors = FALSE)
out <- PrepVCFexport(genotypes, alleles2loc, depth, alnuc, locTable, c(4, 4), TRUE)
out
*/

