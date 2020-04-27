#ifndef HAP2SNP_H
#define HAP2SNP_H

#include <Rcpp.h>
Rcpp::List Hap2SNP(Rcpp::StringVector haps, std::string refhap, int pos);

#endif
