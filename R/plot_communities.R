##' @title Community graphs
##' @param comm_df the output data.frame from \code{\link{find_communities}}
##' @return a list of ggplot2 graphs.
##' @author Jean Monlong
##' @import ggplot2
##' @export
plot_communities <- function(comm_df){
  community = NULL
  ggp.l = list()
  ggp.l$n = ggplot(comm_df, aes(community)) + geom_bar() + theme_bw() + ylab('cell')
  return(ggp.l)
}
