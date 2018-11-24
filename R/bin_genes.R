##' @title Merge consecutive genes into expressed bins
##' @param ge_df the input gene expression with coordinate columns (chr, start, end)
##' and then one column per cell.
##' @param mean_exp the desired minimum mean expression in the bin.
##' @param nb_cores the number of processors to use.
##' @return a data.frame with bin expression.
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @export
bin_genes <- function(ge_df, mean_exp=3, nb_cores=1){
  options('dplyr.show_progress'=FALSE)
  cells = setdiff(colnames(ge_df), c("chr","start","end"))
  m.df = ge_df[, c('chr', 'start', 'end')]
  ## mean gene expression
  m.df$m = apply(as.matrix(ge_df[,cells]),1,mean)
  ## cumulative expression in each chromosome
  m.df = m.df %>% dplyr::group_by(.data$chr) %>%
    dplyr::arrange(.data$start) %>%
    dplyr::mutate(cumm=cumsum(.data$m))
  ## bin cumulative expression
  m.df = dplyr::ungroup(m.df) %>% dplyr::group_by(.data$chr) %>%
    dplyr::mutate(bins=as.character(cut(.data$cumm, breaks=seq(0,max(.data$cumm),
                                                  min(max(.data$cumm), mean_exp)),
                                        include.lowest=TRUE)))
  ## create bin name and merge with main data
  m.df = dplyr::ungroup(m.df) %>%
    dplyr::mutate(bins=paste0(as.character(.data$chr), '_', .data$bins))
  m.df$cumm = NULL
  message(length(unique(m.df$bins)), ' expressed bins created.')
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
        dplyr::ungroup(.data) %>% dplyr::mutate(bins=NULL) %>%
        dplyr::group_by(.data$chr, .data$start, .data$end) %>%
        dplyr::do(mergeGenes.f(.data))
  }, mc.cores=nb_cores)
  return(as.data.frame(data.table::rbindlist(ge_df.m)))
}
