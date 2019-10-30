##' UMAP graphs colored according to QC metrics or sample labels.
##'
##' If the QC data.frame is provided, the distribution of QC metrics is
##' shown to investigate if some communities are batch effects.
##' 
##' If multiple samples were merged (\code{\link{merge_samples}}), the
##' points can be colored by sample of origin by providing the info_df
##' data.frame.
##' 
##' If any qc_df/comm_df/info_df are null but their columns present in umap_df,
##' their corresponding graphs will be generated. Hence a merged version of
##' umap_df, comm_df, qc_df and info_df works (e.g. output of
##' \code{\link{auto_cna_signal}}.
##'
##' @title UMAP graphs
##' @param umap_df the output data.frame from \code{\link{run_umap}} (columns: cell, umap1, umap2)
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
##' ggp.l = plot_umap(umap_df, qc_df, comm_df)
##'
##' ## Print first graph
##' ggpl.l[[1]]
##'
##' ## Customize ggplot
##' ggpl.l[[1]] + ggtitle('First umap graph')
##' }
plot_umap <- function(umap_df, qc_df=NULL, comm_df=NULL, info_df=NULL){
  umap1 = umap2 = tot = mito = G1.S = G2.M = community = NULL
  ptalpha = .5
  if(nrow(umap_df) > 1e4){
    ptalpha = .1
  }
  ggp.l = list()
  ggp.l$nocolor = ggplot(umap_df, aes(umap1, umap2)) + geom_point(alpha=ptalpha) +
    theme_bw()

  ## Sample info
  if(!is.null(info_df) | 'sample' %in% colnames(umap_df)){
    if(!('sample' %in% colnames(umap_df))){
      nrows = nrow(umap_df)
      umap_df = merge(umap_df, info_df)
      if(nrow(umap_df) < nrows){
        warning('Some cells in umap_df are missing from info_df.')
      }
    }
    if(!is.factor(umap_df$sample)){
      umap_df$sample = factor(umap_df$sample)
    }
    cent =  umap_df %>% dplyr::select(.data$sample, .data$umap1, .data$umap2) %>%
      dplyr::group_by(.data$sample) %>% dplyr::summarize_all(stats::median)
    ggp.l$sample = ggplot(umap_df, aes(umap1, umap2, colour=sample)) +
      geom_point(alpha=ptalpha) + theme_bw() +
      ggrepel::geom_label_repel(aes(label=sample), data=cent) + 
      guides(colour=FALSE)
  }

  ## Community info
  if(!is.null(comm_df) | 'community' %in% colnames(umap_df)){
    if(!('community' %in% colnames(umap_df))){
      nrows = nrow(umap_df)
      umap_df = merge(umap_df, comm_df)
      if(nrow(umap_df) < nrows){
        warning('Some cells in umap_df are missing from comm_df.')
      }
    }
    if(!is.factor(umap_df$community)){
      umap_df$community = factor(umap_df$community)
    }
    cent =  umap_df %>% dplyr::select(.data$community, .data$umap1, .data$umap2) %>%
      dplyr::group_by(.data$community) %>% dplyr::summarize_all(stats::median)
    ggp.l$comm = ggplot(umap_df, aes(umap1, umap2, colour=community)) +
      geom_point(alpha=ptalpha) + theme_bw() +
      ggrepel::geom_label_repel(aes(label=community), data=cent) + 
      guides(colour=FALSE)
  }

  ## QC metrics
  if(!is.null(qc_df) | all(c('tot', 'mito') %in% colnames(umap_df))){
    if(!all(c('tot', 'mito') %in% colnames(umap_df))){
      nrows = nrow(umap_df)
      umap_df = merge(umap_df, qc_df)
      if(nrow(umap_df) < nrows){
        warning('Some cells in umap_df are missing from qc_df.')
      }
    }
    ggp.l$depth = ggplot(umap_df, aes(umap1, umap2, colour=log10(tot+1))) +
      geom_point(alpha=ptalpha) + theme_bw() +
      scale_colour_gradientn(name='depth (log10[x+1])', colors=grDevices::terrain.colors(10))
    if(any(umap_df$mito > 0)){
      ggp.l$mito = ggplot(umap_df, aes(umap1, umap2, colour=mito/tot)) +
        geom_point(alpha=ptalpha) + theme_bw() +
        scale_colour_gradientn(name='proportion\nof\nmitochondrial\nRNA',
                               colors=grDevices::terrain.colors(10))
    }
    if(all(c('G1.S', 'G2.M') %in% colnames(umap_df))){
      ggp.l$cellcycle = ggplot(umap_df, aes(umap1, umap2, colour=G1.S+G2.M)) +
        geom_point(alpha=ptalpha) + theme_bw() +
        scale_colour_gradientn(name='cycling score',
                               colors=grDevices::terrain.colors(10))
    }
  }

  return(ggp.l)
}
