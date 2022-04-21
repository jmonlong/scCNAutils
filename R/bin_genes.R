##' @title Merge consecutive genes into expressed bins
##' @param ge_df the input gene expression with coordinate columns (chr, start, end)
##' and then one column per cell.
##' @param mean_exp the desired minimum mean expression in the bin.
##' @param nb_cores the number of processors to use.
##' @param rcpp use Rcpp function. Default is FALSE. More memory-efficient and
##' faster when running on one core.
##' @return a data.frame with bin expression.
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @export
bin_genes <- function(ge_df, mean_exp=3, nb_cores=1, rcpp=FALSE){
  options('dplyr.show_progress'=FALSE)
  cells = setdiff(colnames(ge_df), c("chr","start","end"))
  ge_df = ge_df[order(ge_df$chr, ge_df$start),]
  m.df = ge_df[, c('chr', 'start', 'end')]
  ## mean gene expression
  m.df$m = apply(as.matrix(ge_df[,cells]),1,mean)
  ## cumulative expression in each chromosome
  m.df = m.df %>% dplyr::group_by(.data$chr) %>%
    dplyr::mutate(cumm=cumsum(.data$m),
                  bins=as.character(cut(.data$cumm,
                                        breaks=seq(0,
                                                   max(.data$cumm) + mean_exp,
                                                   min(max(.data$cumm),
                                                       mean_exp)),
                                        include.lowest=TRUE)))
  ## create bin name and merge with main data
  m.df$bins = paste0(as.character(m.df$chr), '_', m.df$bins)
  m.df$cumm = NULL
  message(length(unique(m.df$bins)), ' expressed bins created.')
  if(rcpp){
    m.df$bins = factor(m.df$bins, levels=unique(m.df$bins))
    coord.df = m.df %>% dplyr::ungroup() %>% 
      dplyr::group_by(.data$chr, .data$bins) %>%
      dplyr::summarize(start=min(.data$start), end=max(.data$end))
    geb = binGenesC(as.matrix(ge_df[,cells]), as.character(m.df$bins))
    colnames(geb) = cells
    coord.df$bins = NULL
    res = cbind(as.data.frame(coord.df), geb)
  } else {
    ## sum expression in each bin per cell
    mergeGenes.f <- function(ge_df){
      vv = colSums(as.matrix(ge_df[, cells]))
      as.data.frame(t(vv))
    }
    ge_df.m = parallel::mclapply(unique(ge_df$chr), function(chri){
      ge_df %>% dplyr::filter(.data$chr == chri) %>%
        merge(m.df[which(m.df$chr == chri),]) %>% 
        dplyr::group_by(.data$bins) %>%
          dplyr::mutate(start=min(.data$start), end=max(.data$end)) %>%
          dplyr::ungroup() %>% dplyr::mutate(bins=NULL) %>%
          dplyr::group_by(.data$chr, .data$start, .data$end) %>%
          dplyr::do(mergeGenes.f(.data))
    }, mc.cores=nb_cores)
    res = as.data.frame(data.table::rbindlist(ge_df.m))
  }
  return(res)
}
