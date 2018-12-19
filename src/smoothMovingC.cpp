#include <Rcpp.h>
using namespace Rcpp;

//' Smooth signal
//' 
//' @param ge an ordered gene expression matrix
//' @param winsize the size of the sliding window (assumed to be odd).
// [[Rcpp::export]]
NumericMatrix smoothMovingC(NumericMatrix ge, int winsize) {
  int nsamps = ge.ncol();
  int ngenes = ge.nrow();
  NumericMatrix ges(ngenes, nsamps);
  for(int gene = 0; gene < ngenes; gene++){
    int winS = gene - (winsize - 1) / 2;
    if(winS < 0){
      winS = 0;
    }
    int winE = gene + (winsize - 1) / 2 + 1;
    if(winE > ngenes){
      winE = ngenes;
    }
    // Replace value at window start by median over window
    int mid = (winE - winS - 1) / 2;
    for(int samp = 0; samp < nsamps; samp++){
      // Compute median
      std::vector<double> ms;
      for(int winI = winS; winI < winE; winI++){
	double mi = ge(winI, samp);
	ms.insert(std::upper_bound(ms.begin(), ms.end(), mi), mi);
      }
      ges(gene, samp) = ms[mid];
    }
  }
  return ges;
}
