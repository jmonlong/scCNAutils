##' @title Heatmap of CNA
##' @param cna_df CNA from \code{\link{call_cna}}.
##' @param chrs_order order of the chroosomes in the graph.
##' @return a ggplot2 graph
##' @author Jean Monlong
##' @export
##' @import ggplot2
plot_cna <- function(cna_df, chrs_order=c(1:22, 'X', 'Y')){
  CN = end = prop = start = NULL
  cna_df$CN = factor(cna_df$CN, levels=c('loss','neutral','gain'))
  cna_df$chr = factor(cna_df$chr, levels=chrs_order)

  ggp = ggplot(cna_df, aes(xmin=start, xmax=end, ymin=-.5, ymax=.5, fill=CN, alpha=prop)) +
    geom_rect() + facet_grid(community~chr, scales='free', space='free') + theme_bw() +
    scale_fill_manual(values=c('steelblue','white','indianred')) + xlab('position') +
    ylab('community') + scale_alpha_continuous(name='prop of\nmetacells') + 
    theme(strip.text.y=element_text(angle=0), axis.text.y=element_blank(),
          axis.text.x=element_blank())
  
  return(ggp)
}
