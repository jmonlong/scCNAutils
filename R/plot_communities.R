##' Graphs about the communities found by \code{\link{find_communities}}.
##' For example the size of the communities or the distribution of QC
##' metrics in each community.
##'
##' If the QC data.frame is provided, the distribution of QC metrics is
##' shown to investigate if some communities are batch effects.
##' 
##' If multiple samples were merged (\code{\link{merge_samples}}), the
##' proportion of cells from each sample of origin can be shown if the
##' info_df data.frame is provided.
##' 
##' If qc_df and/or info_df are null but their columns present in comm_df,
##' their corresponding graphs will be generated. Hence a merged version of
##' comm_df, qc_df and info_df works (e.g. output of
##' \code{\link{auto_cna_signal}}.
##'
##' @title Community graphs
##' @param comm_df the output data.frame from \code{\link{find_communities}}
##' @param qc_df a data.frame with QC metrics (output from \code{\link{qc_cells}}).
##' Default is NULL (i.e. not used)
##' @param info_df a data.frame with sample merge info (output from
##' \code{\link{merge_samples}}).
##' Default is NULL (i.e. not used)
##' @return a list of ggplot2 graphs.
##' @author Jean Monlong
##' @import ggplot2
##' @importFrom magrittr %>%
##' @export
##' @examples
##' \dontrun{
##' ggp.l = plot_communities(comm_df, qc_df)
##'
##' ## Print first graph
##' ggpl.l[[1]]
##'
##' ## Customize ggplot
##' ggpl.l[[1]] + ggtitle('First graph about communities')
##' }
plot_communities <- function(comm_df, qc_df=NULL, info_df=NULL){
  community = mito = tot = G1.S = G2.M = prop = NULL
  ggp.l = list()
  ggp.l$n = ggplot(comm_df, aes(community)) + geom_bar() + theme_bw() + ylab('cell')

  ## QC metrics
  if(!is.null(qc_df) | all(c('tot', 'mito') %in% colnames(comm_df))){
    if(!all(c('tot', 'mito') %in% colnames(comm_df))){
      nrows = nrow(comm_df)
      comm_df = merge(comm_df, qc_df)
      if(nrow(comm_df) < nrows){
        warning('Some cells in comm_df are missing from qc_df.')
      }
    }
    ggp.l$depth = ggplot(comm_df, aes(community, tot)) + geom_boxplot() + theme_bw() +
      ylab('depth')
    if(any(comm_df$mito > 0)){
      ggp.l$mito = ggplot(comm_df, aes(community, mito/tot)) + geom_boxplot() + theme_bw() +
        ylab('proportion of mitochondrial RNA')
    }
    if(all(c('G1.S', 'G2.M') %in% colnames(comm_df))){
      ggp.l$ccellcycle = ggplot(comm_df, aes(community, G1.S+G2.M)) + geom_boxplot() +
        theme_bw()
    }
  }

  ## Sample info
  if(!is.null(info_df) | 'sample' %in% colnames(comm_df)){
    if(!('sample' %in% colnames(comm_df))){
      nrows = nrow(comm_df)
      comm_df = merge(comm_df, info_df)
      if(nrow(comm_df) < nrows){
        warning('Some cells in comm_df are missing from info_df.')
      }
    }
    samp.df = comm_df %>% dplyr::group_by(.data$sample) %>%
      dplyr::mutate(ntot=dplyr::n()) %>% dplyr::ungroup(.data) %>% 
      dplyr::group_by(.data$sample, .data$community) %>%
      dplyr::summarize(prop=dplyr::n()/.data$ntot[1])
    ggp.l$sample = ggplot(samp.df, aes(community, prop, fill=sample)) +
      geom_bar(stat='identity') + theme_bw() +
      ylab('proportion of cells in a sample')
  }
  return(ggp.l)
}
