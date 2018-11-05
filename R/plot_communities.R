##' @title Community graphs
##' @param comm_df the output data.frame from \code{\link{find_communities}}
##' @param qc_df a data.frame with QC metrics (output from \code{\link{qc_cells}}).
##' Default is NULL (i.e. not used)
##' @return a list of ggplot2 graphs.
##' @author Jean Monlong
##' @import ggplot2
##' @export
plot_communities <- function(comm_df, qc_df=NULL){
  community = mito = tot = G1.S = G2.M = NULL
  ggp.l = list()
  ggp.l$n = ggplot(comm_df, aes(community)) + geom_bar() + theme_bw() + ylab('cell')
  if(!is.null(qc_df)){
    nrows = nrow(comm_df)
    comm_df = merge(comm_df, qc_df)
    if(nrow(comm_df) < nrows){
      warning('Some cells in comm_df are missing from qc_df.')
    }
    ggp.l$depth = ggplot(comm_df, aes(community, tot)) + geom_boxplot() + theme_bw() +
      ylab('depth')
    if(any(qc_df$mito > 0)){
      ggp.l$mito = ggplot(comm_df, aes(community, mito/tot)) + geom_boxplot() + theme_bw() +
        ylab('proportion of mitochondrial RNA')
    }
    if(all(c('G1.S', 'G2.M') %in% colnames(comm_df))){
      ggp.l$ccellcycle = ggplot(comm_df, aes(community, G1.S+G2.M)) + geom_boxplot() +
        theme_bw()
    }
  }
  return(ggp.l)
}
