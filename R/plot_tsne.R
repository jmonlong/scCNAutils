##' tSNE graphs colored according to QC metrics or sample labels.
##'
##' If the QC data.frame is provided, the distribution of QC metrics is
##' shown to investigate if some communities are batch effects.
##' 
##' If multiple samples were merged (\code{\link{merge_samples}}), the
##' points can be colored by sample of origin by providing the info_df
##' data.frame.
##' 
##' If any qc_df/comm_df/info_df are null but their columns present in tsne_df,
##' their corresponding graphs will be generated. Hence a merged version of
##' tsne_df, comm_df, qc_df and info_df works (e.g. output of
##' \code{\link{auto_cna_signal}}.
##'
##' @title tSNE graphs
##' @param tsne_df the output data.frame from \code{\link{run_tsne}} (columns: cell, tsne1, tsne2)
##' @param qc_df a data.frame with QC metrics (output from \code{\link{qc_cells}}).
##' Default is NULL (i.e. not used)
##' @param comm_df a data.frame with communities (output from \code{\link{find_communities}}).
##' Default is NULL (i.e. not used)
##' @param info_df a data.frame with sample merge info (output from
##' \code{\link{merge_samples}}).
##' @return a list of ggplot objects
##' @author Jean Monlong
##' @import ggplot2
##' @importFrom magrittr %>%
##' @export
##' @examples
##' \dontrun{
##' ggp.l = plot_tsne(tsne_df, qc_df, comm_df)
##'
##' ## Print first graph
##' ggpl.l[[1]]
##'
##' ## Customize ggplot
##' ggpl.l[[1]] + ggtitle('First tSNE graph')
##' }
plot_tsne <- function(tsne_df, qc_df=NULL, comm_df=NULL, info_df=NULL){
  tsne1 = tsne2 = tot = mito = G1.S = G2.M = community = NULL
  ptalpha = .5
  if(nrow(tsne_df) > 1e4){
    ptalpha = .1
  }
  ggp.l = list()
  ggp.l$nocolor = ggplot(tsne_df, aes(tsne1, tsne2)) + geom_point(alpha=ptalpha) +
    theme_bw()

  ## Sample info
  if(!is.null(info_df) | 'sample' %in% colnames(tsne_df)){
    if(!('sample' %in% colnames(tsne_df))){
      nrows = nrow(tsne_df)
      tsne_df = merge(tsne_df, info_df)
      if(nrow(tsne_df) < nrows){
        warning('Some cells in tsne_df are missing from info_df.')
      }
    }
    if(!is.factor(tsne_df$sample)){
      tsne_df$sample = factor(tsne_df$sample)
    }
    cent =  tsne_df %>% dplyr::select(.data$sample, .data$tsne1, .data$tsne2) %>%
      dplyr::group_by(.data$sample) %>% dplyr::summarize_all(stats::median)
    ggp.l$sample = ggplot(tsne_df, aes(tsne1, tsne2, colour=sample)) +
      geom_point(alpha=ptalpha) + theme_bw() +
      ggrepel::geom_label_repel(aes(label=sample), data=cent) + 
      guides(colour=FALSE)
  }

  ## Community info
  if(!is.null(comm_df) | 'community' %in% colnames(tsne_df)){
    if(!('community' %in% colnames(tsne_df))){
      nrows = nrow(tsne_df)
      tsne_df = merge(tsne_df, comm_df)
      if(nrow(tsne_df) < nrows){
        warning('Some cells in tsne_df are missing from comm_df.')
      }
    }
    if(!is.factor(tsne_df$community)){
      tsne_df$community = factor(tsne_df$community)
    }
    cent =  tsne_df %>% dplyr::select(.data$community, .data$tsne1, .data$tsne2) %>%
      dplyr::group_by(.data$community) %>% dplyr::summarize_all(stats::median)
    ggp.l$comm = ggplot(tsne_df, aes(tsne1, tsne2, colour=community)) +
      geom_point(alpha=ptalpha) + theme_bw() +
      ggrepel::geom_label_repel(aes(label=community), data=cent) + 
      guides(colour=FALSE)
  }

  ## QC metrics
  if(!is.null(qc_df) | all(c('tot', 'mito') %in% colnames(tsne_df))){
    if(!all(c('tot', 'mito') %in% colnames(tsne_df))){
      nrows = nrow(tsne_df)
      tsne_df = merge(tsne_df, qc_df)
      if(nrow(tsne_df) < nrows){
        warning('Some cells in tsne_df are missing from qc_df.')
      }
    }
    ggp.l$depth = ggplot(tsne_df, aes(tsne1, tsne2, colour=tot)) +
      geom_point(alpha=ptalpha) + theme_bw() +
      scale_colour_gradientn(name='depth', colors=grDevices::terrain.colors(10))
    if(any(tsne_df$mito > 0)){
      ggp.l$mito = ggplot(tsne_df, aes(tsne1, tsne2, colour=mito/tot)) +
        geom_point(alpha=ptalpha) + theme_bw() +
        scale_colour_gradientn(name='proportion\nof\nmitochondrial\nRNA',
                               colors=grDevices::terrain.colors(10))
    }
    if(all(c('G1.S', 'G2.M') %in% colnames(tsne_df))){
      ggp.l$cellcycle = ggplot(tsne_df, aes(tsne1, tsne2, colour=G1.S+G2.M)) +
        geom_point(alpha=ptalpha) + theme_bw() +
        scale_colour_gradientn(name='cycling score',
                               colors=grDevices::terrain.colors(10))
    }
  }

  return(ggp.l)
}
