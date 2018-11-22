##' UMAP on the PCA results.
##'
##' This functions depends on Python and UMAP being installed. Make sure umap-mearn,
##' sklearn, numpy, scipy and pandas are installed. For example with something like:
##' 'pip install sklearn numpy scipy pandas umap-learn'.
##' @title Run UMAP
##' @param pca_o the output of \code{\link{run_pca}} 
##' @param nb_pcs the number of PCs to use. Default 10.
##' @param nb_neighbors the number of neighbors. Default 5.
##' @return a data.frame with columns: cell, umap1, umap2
##' @author Jean Monlong
##' @export
run_umap <- function(pca_o, nb_pcs=10, nb_neighbors=5){
  nb_pcs = min(nb_pcs, ncol(pca_o$x))
  temp.pref = paste0('tempforumap', stats::runif(1, 0, 1e4))
  ## Write Python script
  write(c('import umap',
          'import pandas as pd',
          paste0('data = pd.read_csv("', temp.pref, '.csv", sep=",", header=None)'),
          paste0('reducer = umap.UMAP(n_neighbors=', nb_neighbors, ')'),
          'embedding = reducer.fit_transform(data)',
          paste0('pd.DataFrame(embedding).to_csv("', temp.pref, '.csv", header=False, index=False)')),
        file=paste0(temp.pref, '.py'))
  ## Write input data
  utils::write.table(pca_o$x[,1:nb_pcs], file=paste0(temp.pref, '.csv'),
                     row.names=FALSE, col.names=FALSE, sep=',')
  ## Run command
  runcmd = tryCatch({
    system2('python', paste0(temp.pref, '.py'))
  }, error = function(err){
    stop('Error when running UMAP in Python. Make sure umap-mearn, sklearn, numpy, scipy and pandas are installed. For example something like:\n',
         'pip install sklearn numpy scipy pandas umap-learn\n\n',
         'The python error was:\n',
         err)
  })
  ## Load results and cleanup
  umap.df = utils::read.csv(paste0(temp.pref, '.csv'), header=FALSE, as.is=TRUE)
  file.remove(paste0(temp.pref, c('.csv', '.py')))
  colnames(umap.df) = paste0('umap', 1:2)
  umap.df$cell = rownames(pca_o$x)
  return(umap.df)
}
