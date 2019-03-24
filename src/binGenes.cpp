#include <Rcpp.h>
using namespace Rcpp;

//' Bin genes
//' 
//' @param ge an ordered gene expression matrix
//' @param bins a vector with bin ids
// [[Rcpp::export]]
NumericMatrix binGenesC(NumericMatrix ge, CharacterVector bins) {
  int nsamps = ge.ncol();
  int ngenes = ge.nrow();
  // Get number of unique bins
  int nbins = 0;
  String curBin = "";
  for(int gene = 0; gene < ngenes; gene++){
    if(curBin != bins[gene]){
      nbins++;
      curBin = bins[gene];
    }
  }
  // Init results
  NumericMatrix gem(nbins, nsamps);
  // Iterate over genes
  int gemIdx = 0;
  curBin = bins[0];
  NumericVector curExp(nsamps);
  for(int samp = 0; samp < nsamps; samp++){
    curExp[samp] = ge(0, samp);
  }
  for(int gene = 1; gene < ngenes; gene++){
    if(bins[gene] == curBin){
      // Add expression
      for(int samp = 0; samp < nsamps; samp++){
	curExp[samp] = curExp[samp] + ge(gene, samp);
      }
    } else {
      // Update results
      for(int samp = 0; samp < nsamps; samp++){
	gem(gemIdx, samp) = curExp[samp];
      }
      gemIdx++;
      curBin = bins[gene];
      // Prepare vector 
      for(int samp = 0; samp < nsamps; samp++){
	curExp[samp] = ge(gene, samp);
      }
    }
  }
  // Update last bin
  for(int samp = 0; samp < nsamps; samp++){
    gem(gemIdx, samp) = curExp[samp];
  }
  return gem;
}
