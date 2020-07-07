##' @title QC graphs
##' @param qc_df the output data.frame from \code{\link{qc_cells}}
##' @param info_df a data.frame with sample merge info (output from
##' \code{\link{merge_samples}}). Default is NULL (i.e. not used)
##' @return a list of ggplots
##' @author Jean Monlong
##' @import ggplot2
##' @export
##' @examples
##' \dontrun{
##' ggp.l = plot_qc_cells(qc_df)
##'
##' ## Print first graph
##' ggpl.l[[1]]
##'
##' ## Customize ggplot
##' ggpl.l[[1]] + ggtitle('First QC graph')
##' }
plot_qc_cells <- function(qc_df, info_df=NULL){
  tot = zeros = mito = NULL
  ggp.l = list()
  ggp.l$depth = ggplot(qc_df, aes(x=tot)) + geom_histogram(bins=30) + theme_bw() +
    xlab('depth') + ylab('cell')
  ggp.l$depth.log = ggplot(qc_df, aes(x=log10(tot+1))) + geom_histogram(bins=30) + theme_bw() +
    xlab('depth (log10[x+1])') + ylab('cell')
  ggp.l$depth.w90 = ggplot(qc_df, aes(x=winsor(tot, uq=.9))) + geom_histogram(bins=30) + theme_bw() +
    xlab('depth (winsorized at 90th percentile)') + ylab('cell')
  ggp.l$zeros = ggplot(qc_df, aes(x=zeros)) + geom_histogram(bins=30) + theme_bw() +
    xlab('number of 0s') + ylab('cell')
  if(any(qc_df$mito > 0)){
    ggp.l$mito = ggplot(qc_df, aes(x=mito/tot)) + geom_histogram(bins=30) + theme_bw() +
      xlab('proportion of mitochondrial RNA') + ylab('cell')
  }

  ## Sample info
  if(!is.null(info_df) | 'sample' %in% colnames(qc_df)){
    if(!('sample' %in% colnames(qc_df))){
      nrows = nrow(qc_df)
      qc_df = merge(qc_df, info_df)
      if(nrow(qc_df) < nrows){
        warning('Some cells in qc_df are missing from info_df.')
      }
    }
    ggp.l$depth.sample = ggplot(qc_df, aes(x=tot)) + geom_histogram(bins=30) +
      theme_bw() + xlab('depth') + ylab('cell') +
      facet_grid(sample~., scales='free') +
      theme(strip.text.y=element_text(angle=0))
    ggp.l$depth.log.sample = ggplot(qc_df, aes(x=log10(tot+1))) + geom_histogram(bins=30) + theme_bw() +
      xlab('depth (log10[x+1])') + ylab('cell') +
      facet_grid(sample~., scales='free') +
      theme(strip.text.y=element_text(angle=0)) +
      scale_x_continuous(breaks=0:10)
    ggp.l$depth.w90.sample = ggplot(qc_df, aes(x=winsor(tot, uq=.9))) + geom_histogram(bins=30) + theme_bw() +
      xlab('depth (winsorized at 90th percentile)') + ylab('cell') +
      facet_grid(sample~., scales='free') +
      theme(strip.text.y=element_text(angle=0))
    ggp.l$zeros.sample = ggplot(qc_df, aes(x=zeros)) + geom_histogram(bins=30) + theme_bw() +
      xlab('number of 0s') + ylab('cell') +
      facet_grid(sample~., scales='free') +
      theme(strip.text.y=element_text(angle=0))
    if(any(qc_df$mito > 0)){
      ggp.l$mito.sample = ggplot(qc_df, aes(x=mito/tot)) + geom_histogram(bins=30) +
        theme_bw() + ylab('cell') +
        xlab('proportion of mitochondrial RNA') +
        facet_grid(sample~., scales='free') +
        theme(strip.text.y=element_text(angle=0))
    }
  }

  return(ggp.l)
}
