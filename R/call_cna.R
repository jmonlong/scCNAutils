##' @title Call CNA
##' @param z_df the Z-scores, from \code{\link{zscore}}.
##' @param trans_prob the transition probability for the HMM.
##' @param nb_cores the number of processor to use.
##' @param mc_info the information about the metacells, if relevant. Default is NULL.
##' @return a data.frame with the CNA calls.
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @export
call_cna <- function(z_df, trans_prob=1e-4, nb_cores=1, mc_info=NULL){
  . = V1 = end = queryHits = subjectHits = NULL
  freq.range <- function(df, nb.samp){
    gr.all = GenomicRanges::GRanges(df$chr,
                                    IRanges::IRanges(df$start, df$end),
                                    cell=df$cell)
    gr.d = GenomicRanges::disjoin(gr.all)
    fr.df = GenomicRanges::as.data.frame(gr.d)[, 1:3]
    colnames(fr.df)[1] = "chr"
    fr.df$chr = as.character(fr.df$chr)
    ol = as.data.frame(GenomicRanges::findOverlaps(gr.d,gr.all))
    ol = ol %>% dplyr::group_by(queryHits) %>% dplyr::summarize(n=length(unique(gr.all$cell[subjectHits])))
    fr.df$nb = 0
    fr.df$nb[ol$queryHits] = ol$n
    fr.df$prop = fr.df$nb/nb.samp
    fr.df
  }
  
  cells = setdiff(colnames(z_df), c("chr","start","end", 'symbol'))
  ## Call CNV in each cell
  hmm.o = parallel::mclapply(cells, function(cell){
    cn.df = cnaHMM(z_df[,c('chr','start','end',cell)], trans.prob=trans_prob)
    rle.o = rle(paste(cn.df$chr, cn.df$CN))
    ## Merge into segments
    cn.df$seg = rep(1:length(rle.o$lengths), rle.o$lengths)
    seg.df = dplyr::ungroup(cn.df) %>% dplyr::group_by(.data$chr, .data$CN, .data$seg) %>%
      dplyr::summarize(start=min(.data$start), end=max(.data$end),
                       mean=mean(.data$mean), length=dplyr::n())
    seg.df$seg = NULL
    seg.df = seg.df[order(seg.df$chr, seg.df$start),]
    ## Filter short segments and merge again
    seg.c = seg.df[which(seg.df$length>5),]
    seg.c = seg.c[order(seg.c$chr, seg.c$start),]
    rle.o = rle(paste(seg.c$chr, seg.c$CN))
    seg.c$seg = rep(1:length(rle.o$lengths), rle.o$lengths)
    seg.c = dplyr::ungroup(seg.c) %>% dplyr::group_by(.data$chr, .data$CN, .data$seg) %>%
      dplyr::summarize(start=min(.data$start), end=max(end),
                       mean=sum(.data$mean*.data$length)/sum(.data$length),
                       length=sum(.data$length))
    seg.c$seg = NULL
    seg.c = seg.c[order(seg.c$chr, seg.c$start),]
    seg.c$cell = as.character(cell)
    cn.df$cell = as.character(cell)
    list(seg=seg.c, cn=cn.df)
  }, mc.cores=nb_cores)
  seg.df = do.call(rbind, lapply(hmm.o, function(e)e$seg))
  if(!is.null(mc_info)){
    seg.df = merge(seg.df, mc_info)
    nb_metacells = unique(table(mc_info$community))
    if(length(nb_metacells)>1){
      warning('Different number of metacells per community. Normal?')
    }
    seg.df = seg.df %>% dplyr::group_by(.data$community, .data$CN) %>%
      dplyr::do(freq.range(.data, nb_metacells))
  }
  return(seg.df)
}
