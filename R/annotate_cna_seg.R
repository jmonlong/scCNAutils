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
  if(length(neutral.gr)==0){
    warning('No copy neutral segment found. This is fishy, community(ies) ', unique(seg.df$community), ' may be mostly made of low quality/depth cells.\n Suggestions: 1) filter them out in auto_cna_signal with min_total_exp; 2) increase metacell sizes with metacell_size or bin_mean_exp.')
    pvs = rep(NA, nrow(seg.df))
  } else {
    pvs = sapply(1:nrow(seg.df), function(ii){
      seg.mean = hmm.gr$mean[which(IRanges::overlapsAny(hmm.gr, seg.gr[ii]))]
      neut.d = GenomicRanges::distance(seg.gr[ii], neutral.gr)
      ## handle cases where there are no "neutral CN" segments in this chromomosome
      if(all(is.na(neut.d))){
        neut.d = rep(1, length(neutral.gr))
      }
      ## otherwise put them at the end of the preference list
      neut.d = ifelse(is.na(neut.d), max(neut.d, na.rm=TRUE) + 1, neut.d)
      ## put the segments right next to the one of interest even farther
      ## down the list (bc not independent enough)
      neut.d = ifelse(neut.d==0, max(neut.d, na.rm=TRUE) + 1, neut.d)
      ## retrieve mean in neutral segments nearby, or in other chromosomes,
      ## or in worst case right next to the segment of interest
      cont.mean = neutral.gr$mean[utils::head(order(neut.d), length(seg.mean))]
      return(suppressWarnings(stats::wilcox.test(seg.mean, cont.mean)$p.value))
    })
  }
  seg.df$wt.pv = pvs
  return(seg.df)
}
