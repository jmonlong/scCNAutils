##' @title Community detection
##' @param pca_o the output of \code{\link{run_pca}} 
##' @param nb_pcs the number of PCs to use. Default 10.
##' @param k the number of nearest neighbor for the KNN graph. Default 100.
##' @return a data.frame with two columns: 'cell' and 'community'.
##' @author Jean Monlong
##' @export
find_communities <- function(pca_o, nb_pcs=10, k=100){
  nb_pcs = min(nb_pcs, ncol(pca_o$x))
  knn.o = FNN::get.knn(as.matrix(pca_o$x[,1:nb_pcs]), k = k)
  knn.df = data.frame(from=rep(1:nrow(knn.o$nn.index), k),
                      to=as.vector(knn.o$nn.index),
                      weight=1/(1+as.vector(knn.o$nn.dist)))
  nw = igraph::graph_from_data_frame(knn.df, directed=FALSE)
  nw = igraph::simplify(nw)
  wc = igraph::cluster_louvain(nw)
  comm.df = data.frame(cell=rownames(pca_o$x),
                       community=as.factor(igraph::membership(wc)),
                       stringsAsFactors=FALSE)  
  return(comm.df)
}
