##' @title Define cycling cells
##' @param qc_df the output data.frame from \code{\link{qc_cells}} (ran with a
##' non-null \emph{cell_cycle} parameter)
##' @param sd_th the number of SD used for the thresholds.
##' @return a list with
##' \item{cells.noc}{a vector with the names of non-cycling cells}
##' \item{graphs}{a list of ggplot2 graphs}
##' @author Jean Monlong
##' @import ggplot2
##' @export
define_cycling_cells <- function(qc_df, sd_th=3){
  ## Compute thresholds
  g1s.th = stats::median(qc_df$G1.S) + sd_th*stats::mad(qc_df$G1.S)
  g2m.th = stats::median(qc_df$G2.M) + sd_th*stats::mad(qc_df$G2.M)
  ## List cells below thresholds
  cells.noc = qc_df$cell[which(qc_df$G1.S<g1s.th & qc_df$G2.M<g2m.th)]
  ## Graphs
  ptalpha = ifelse(nrow(qc_df)<1e4, .5, .1)
  G1.S = G2.M = NULL
  ggp.l = list()
  ggp.l$g1.s = ggplot(qc_df, aes(G1.S)) + geom_histogram() + theme_bw() +
    geom_vline(xintercept=g1s.th, linetype=2) + ylab('cell')
  ggp.l$g2.m = ggplot(qc_df, aes(G2.M)) + geom_histogram() + theme_bw() +
    geom_vline(xintercept=g2m.th, linetype=2) + ylab('cell')
  ggp.l$both = ggplot(qc_df, aes(G1.S, G2.M)) + geom_point(alpha=ptalpha) + theme_bw() +
    geom_hline(yintercept=g2m.th, linetype=2) +
    geom_vline(xintercept=g1s.th, linetype=2)
  return(list(cells.noc=cells.noc, graphs=ggp.l))
}
