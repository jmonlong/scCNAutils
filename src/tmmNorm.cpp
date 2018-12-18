#include <Rcpp.h>
using namespace Rcpp;

//' Compute the normalization factor
//' 
//' @param ge gene expression amtrix
//' @param cont index of the control sample. R index so starting at 1.
//' @export
// [[Rcpp::export]]
NumericMatrix tmmNormC(NumericMatrix ge, int cont) {
  // R index -> Cpp index
  cont = cont - 1;
  // Init
  int nsamps = ge.ncol();
  int ngenes = ge.nrow();
  NumericMatrix gen(ngenes, nsamps);
  // Iterate over sample (columns)
  for(int samp = 0; samp < nsamps; samp++){
    // Compute normalization factor
    std::vector<double> ms;
    for(int i = 0; i < ngenes; i++){
      if(ge(i, samp) > 0 & ge(i, cont) > 0){
	double mi = log(ge(i, samp) / ge(i, cont));
	ms.insert(std::upper_bound(ms.begin(), ms.end(), mi), mi);
      }
    }
    double nf = 0;
    int nonzero = 0;
    for(int i = ms.size() * 0.3; i < ms.size() * 0.7; i++){
      nf += ms[i];
      nonzero++;
    }
    nf = exp(nf/nonzero);
    // Normalize column
    for(int i = 0; i < ngenes; i++){
      gen(i, samp) = ge(i, samp) / nf;
    }
  }
  return gen;
}
