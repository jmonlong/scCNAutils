##' Annotate the CN segment predicted by the HMM. The signal in the segment is compared
##' to the signal in "neutral" segments nearby using a Wilcoxon test.
##' @title Annotate CN segments
##' @param seg.df a data.frame with segment information
##' @param hmm.df a data.frame with bin information
##' @return an annotated version of seg.df with a column wt.pv wit the pvalue of the
##' Wilcoxon test.
##' @author Jean Monlong
##' @importFrom magrittr %>%
annotate_cna_seg <- function(seg.df, hmm.df){
  seg.gr = GenomicRanges::makeGRangesFromDataFrame(seg.df, keep.extra.columns=TRUE)
  hmm.gr = GenomicRanges::makeGRangesFromDataFrame(hmm.df, keep.extra.columns=TRUE)
  neutral.gr = hmm.gr[which(hmm.gr$CN=='neutral')]
  pvs = sapply(1:nrow(seg.df), function(ii){
    seg.mean = hmm.gr$mean[which(IRanges::overlapsAny(hmm.gr, seg.gr[ii]))]
    neut.d = GenomicRanges::distance(seg.gr[ii], neutral.gr)
    cont.mean = neutral.gr$mean[utils::head(order(neut.d), length(seg.mean))]
    stats::wilcox.test(seg.mean, cont.mean)$p.value
  })
  seg.df$wt.pv = pvs
  return(seg.df)
}
