##' @title QC graphs
##' @param qc_df the output data.frame from qc_cells
##' @return a list of ggplots
##' @author Jean Monlong
##' @import ggplot2
plot_qc_cells <- function(qc_df){
  tot = zeros = mito = NULL
  ggp.l = list()
  ggp.l$depth = ggplot(qc_df, aes(x=tot)) + geom_histogram() + theme_bw() +
    xlab('depth') + ylab('cell')
  ggp.l$zeros = ggplot(qc_df, aes(x=zeros)) + geom_histogram() + theme_bw() +
    xlab('number of 0s') + ylab('cell')
  if(any(qc_df$mito > 0)){
    ggp.l$mito = ggplot(qc_df, aes(x=mito/tot)) + geom_histogram() + theme_bw() +
      xlab('proportion of mitochondrial RNA') + ylab('cell')
  }
  return(ggp.l)
}
