##' Calls CNA using a HMM approach considering multiple samples at
##' the same time.
##'
##' @title Call CNA
##' @param z_df the Z-scores, from \code{\link{zscore}}.
##' @param mc_info the information about the metacells, if relevant. Default is NULL.
##' @param trans_prob the transition probability for the HMM.
##' @param nb_cores the number of processor to use.
##' @return a data.frame with the CNA calls.
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @export
call_cna_multisamps <- function(z_df, mc_info, trans_prob=1e-4, nb_cores=1){
  . = V1 = end = queryHits = subjectHits = NULL
  options('dplyr.show_progress'=FALSE)
  
  cells = setdiff(colnames(z_df), c("chr","start","end", 'symbol'))
  mc_info = mc_info[which(mc_info$cell %in% cells),]
  comms = unique(mc_info$community)
  ## Call CNV in each community across multiple metacells
  hmm.o = parallel::mclapply(comms, function(comm){
    cell = mc_info$cell[which(mc_info$community==comm)]
    cn.df = cnaHMM(z_df[,c('chr','start','end',cell)], trans.prob=trans_prob)
    rle.o = rle(paste(cn.df$chr, cn.df$CN))
    ## Merge into segments
    cn.df$seg = rep(1:length(rle.o$lengths), rle.o$lengths)
    seg.df = dplyr::ungroup(cn.df) %>% dplyr::group_by(.data$chr, .data$CN, .data$seg) %>%
      dplyr::summarize(start=min(.data$start), end=max(.data$end),
                       mean=mean(.data$mean), length=dplyr::n())
    seg.df$seg = NULL
    seg.df = seg.df[order(seg.df$chr, seg.df$start),]
    seg.df$community = as.character(comm)
    cn.df$community = as.character(comm)
    list(seg=seg.df, cn=cn.df)
  }, mc.cores=nb_cores)
  seg.df = do.call(rbind, lapply(hmm.o, function(e)e$seg))
  seg.df$pass.filter = seg.df$length>5
  cn.df = do.call(rbind, lapply(hmm.o, function(e)e$cn))
  return(list(seg.df=seg.df, hmm.df=cn.df))
}
