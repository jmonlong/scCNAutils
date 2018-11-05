##' @title tSNE graphs
##' @param tsne_df the output data.frame from \code{\link{run_tsne}} (columns: cell, tsne1, tsne2)
##' @param qc_df a data.frame with QC metrics (output from \code{\link{qc_cells}}). Default is NULL
##' (i.e. not used)
##' @return a list of ggplot objects
##' @author Jean Monlong
##' @import ggplot2
##' @export
plot_tsne <- function(tsne_df, qc_df=NULL){
  tsne1 = tsne2 = tot = mito = NULL
  ptalpha = .5
  if(nrow(tsne_df) > 1e4){
    ptalpha = .1
  }
  ggp.l = list()
  ggp.l$nocolor = ggplot(tsne_df, aes(tsne1, tsne2)) + geom_point(alpha=ptalpha) +
    theme_bw()
  if(!is.null(qc_df)){
    nrows = nrow(tsne_df)
    tsne_df = merge(tsne_df, qc_df)
    if(nrow(tsne_df) < nrows){
      warning('Some cells in tsne_df are missing from qc_df.')
    }
    ggp.l$depth = ggplot(tsne_df, aes(tsne1, tsne2, colour=tot)) +
      geom_point(alpha=ptalpha) + theme_bw() +
      scale_colour_gradientn(name='depth', colors=grDevices::terrain.colors(10))
    if(any(qc_df$mito > 0)){
      ggp.l$mito = ggplot(tsne_df, aes(tsne1, tsne2, colour=mito/tot)) +
        geom_point(alpha=ptalpha) + theme_bw() +
        scale_colour_gradientn(name='proportion\nof\nmitochondrial\nRNA',
                               colors=grDevices::terrain.colors(10))
    }
  }
  return(ggp.l)
}
