##' @title Run tSNE
##' @param pca_o the output of \code{\link{run_pca}} 
##' @param nb_pcs the number of PCs to use. Default 10.
##' @param nb_it the number of iterations. Default 1000.
##' @return a data.frame with columns: cell, tsne1, tsne2
##' @author Jean Monlong
run_tsne <- function(pca_o, nb_pcs=10, nb_it=1000){
  nb_pcs = min(nb_pcs, ncol(pca_o$x))
  tsne.o = Rtsne::Rtsne(pca_o$x[,1:nb_pcs], pca=FALSE, max_iter=nb_it)
  tsne.df = data.frame(cell=rownames(pca_o$x), tsne1=tsne.o$Y[,1], tsne2=tsne.o$Y[,2],
                       stringsAsFactors=FALSE)
  return(tsne.df)
}
