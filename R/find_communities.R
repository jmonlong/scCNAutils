##' Build a KNN graph and run Louvain algorithm for community detection.
##'
##' @title Community detection
##' @param pca_o the output of \code{\link{run_pca}} 
##' @param nb_pcs the number of PCs to use. Default 10.
##' @param k the number of nearest neighbor for the KNN graph. Default 100.
##' @param gamma a vector of gamma. Default 1.
##' @param nreps the number of repetition for each gamma, Default 1.
##' @param nb_cores the number of processors to use. Default is 1.
##' @return a list with:
##' \item{comm}{a data.frame with two columns: 'cell' and 'community'.}
##' \item{comm.all}{a matrix with communities for each gamma}
##' \item{gamma}{the list of input gamma corresponding to each comm.all column.}
##' \item{best.gamma}{the gamma resulting on the highest ARI mean}
##' \item{ari.df}{data.frame with ARI stats for each gamma}
##' @author Jean Monlong
##' @export
find_communities <- function(pca_o, nb_pcs=10, k=100, gamma=1, nreps=1,
                             nb_cores=1){
  nb_pcs = min(nb_pcs, ncol(pca_o$x))
  all.comm = best.gamma = ari.df = NULL

  knn.o = FNN::get.knn(as.matrix(pca_o$x[,1:nb_pcs]), k = k)
  knn.df = data.frame(from=rep(1:nrow(knn.o$nn.index), k),
                      to=as.vector(knn.o$nn.index),
                      weight=1/(1+as.vector(knn.o$nn.dist)))
  nw = igraph::graph_from_data_frame(knn.df, directed=FALSE)
  nw = igraph::simplify(nw)

  if(length(gamma)==1 && gamma==1 && nreps==1){
    ## Default: gamma=1, only once -> Use igraph
    wc = igraph::cluster_louvain(nw)
    comm.df = data.frame(cell=rownames(pca_o$x),
                         community=as.factor(igraph::membership(wc)),
                         stringsAsFactors=FALSE)
  } else {
    ## Exploration: multiple gammas and/or permutations -> Use pyhton wrapper
    wc = run_louvain(nw, gamma=gamma, nreps=nreps, nb_cores=nb_cores)

    ## Adjusted Rand Index
    adjRandInd <- function(x,y){
      cont = table(x,y)
      xsums.ch2 = sum(sapply(rowSums(cont), choose, k=2))
      ysums.ch2 = sum(sapply(colSums(cont), choose, k=2))
      n = length(x)
      index = sum(sapply(as.vector(cont), choose, k=2))
      expindex = xsums.ch2*ysums.ch2/choose(n,2)
      maxindex = (xsums.ch2+ysums.ch2)/2
      (index-expindex)/(maxindex-expindex)
    }

    pw.ss = 100 ## only use up to pw.ss random pairs
    exp.l = parallel::mclapply(unique(wc$gamma), function(gamma){
      cl.mat = wc$comm[,which(wc$gamma==gamma), drop=FALSE]
      nbcls = apply(cl.mat, 2, function(x)length(unique(x)))
      if(ncol(cl.mat) == 1){
        comm.rep = cl.mat[,1]
        ari.df = data.frame(gamma=gamma, ari.mean=NA,
                            ari.sd=NA,
                            nb.cluster=mean(nbcls))
      } else {
        ## Stats: number of communities and adusted Rand index
        pw.ij = utils::combn(ncol(cl.mat),2)
        if(ncol(pw.ij) > pw.ss){
          pw.ij = pw.ij[, sample.int(ncol(pw.ij), pw.ss)]
        }
        ari = sapply(1:ncol(pw.ij), function(ij){
          adjRandInd(cl.mat[,pw.ij[1,ij]], cl.mat[,pw.ij[2,ij]])
        }) 
        ari.df = data.frame(gamma=gamma, ari.mean=mean(ari),
                            ari.sd=stats::sd(ari),
                            nb.cluster=mean(nbcls))
        ## Pick representative run as the one with highest mean ARI
        sim.mat = matrix(NA, ncol(cl.mat), ncol(cl.mat))
        sim.mat[pw.ij[1,] + (pw.ij[2,]-1)*ncol(cl.mat)] = ari
        sim.mat[pw.ij[2,] + (pw.ij[1,]-1)*ncol(cl.mat)] = ari
        sim.mean = apply(sim.mat, 1, mean, na.rm=TRUE)
        comm.rep = cl.mat[,which.max(sim.mean)]
      }
      return(list(comm=comm.rep, ari.df=ari.df))
    }, mc.cores=nb_cores)
    ari.df = do.call(rbind, lapply(exp.l, function(e) e$ari.df))
    best.gamma.idx = order(ari.df$ari.mean, decreasing=TRUE)[1]
    all.comm = matrix(unlist(lapply(exp.l, function(e) e$comm)),
                      ncol=length(exp.l))
    comm.df = data.frame(cell=rownames(pca_o$x),
                         community=as.factor(all.comm[,best.gamma.idx]),
                         stringsAsFactors=FALSE)
    best.gamma = gamma[best.gamma.idx]
  }
  
  return(list(comm=comm.df, comm.all=all.comm, gamma=gamma,
              best.gamma=best.gamma, ari.df=ari.df))
}
