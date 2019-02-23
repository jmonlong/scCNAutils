##' The new bins are overlapped with the regions in cov_df to compute a weight.
##' The coverage in the new bin is the weighted sum of the coverage in overlapping
##' regions of cov_df.
##'
##' Coordinates in data.frame are expected to be defined by columns names 'chr',
##' 'start' and 'end'.
##' @title Re-bin coverage data
##' @param cov_df a data.frame with coverage information.
##' @param bins a data.frame or GRanges object with the new bins.
##' @return a data.frame with new bins
##' @author Jean Monlong
##' @export
rebin_cov <- function(cov_df, bins){
  ## GRanges
  if(is.data.frame(bins)){
    bins = GenomicRanges::makeGRangesFromDataFrame(bins)
  }  
  bins.cov = GenomicRanges::makeGRangesFromDataFrame(cov_df)

  ## Overlap bins and compute proportion overlapped by new bins
  ol = GenomicRanges::findOverlaps(bins, bins.cov)
  ol = GenomicRanges::as.data.frame(ol)
  ol$qsw = GenomicRanges::width(GenomicRanges::pintersect(bins[ol$queryHits], bins.cov[ol$subjectHits]))
  ol = dplyr::mutate(dplyr::group_by(ol, .data$queryHits), ol=.data$qsw/sum(.data$qsw))

  ## Create a weight matrix from the overlap proportion
  merge.mat = matrix(0, length(bins), nrow(cov_df))
  merge.mat[cbind(ol$queryHits, ol$subjectHits)] = ol$ol

  ## Weighted sum of the coverage
  cells = setdiff(colnames(cov_df), c('chr', 'start', 'end'))
  cov.m = merge.mat %*% as.matrix(cov_df[,cells])
  cov.m = as.data.frame(cov.m)
  colnames(cov.m) = cells
  bins.df = data.frame(chr=as.character(GenomicRanges::seqnames(bins)),
                       start=GenomicRanges::start(bins),
                       end=GenomicRanges::end(bins),
                       stringsAsFactors=FALSE)
  cov.m = cbind(bins.df, cov.m)

  return(cov.m)
}
