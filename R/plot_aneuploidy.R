##' Graphs showing the median expression in each chromosome for each community.
##' @title Aneuploidy graph
##' @param ge_df a data.frame with gene expression (better if binned and
##' normalized).
##' @param comm_df a data.frame with community information for each cell.
##' @param baseline_cells cells to use as baseline.
##' @param baseline_communities the communities to use as baseline.
##' @param max_cells the maximum number of cells to consider in the boxplot of
##' each community. Default: 100.
##' @param chrs_order order of the chromosomes in the graph.
##' @return a list of ggplot2 object, one for each chromosome.
##' @author Jean Monlong
##' @export
##' @import ggplot2
##' @importFrom magrittr %>%
plot_aneuploidy <- function(ge_df, comm_df=NULL, baseline_cells=NULL, baseline_communities=NULL, max_cells=100, chrs_order=c(1:22, 'X', 'Y')){
  comm_df = comm_df[which(comm_df$cell %in% colnames(ge_df)),]
  if(is.null(baseline_cells) & !is.null(baseline_communities)){
    baseline_cells = comm_df$cell[which(comm_df$community %in% baseline_communities)]
  }
  if(is.null(baseline_cells)){
    baseline_cells = comm_df$cell
  }
  chrs_order = intersect(chrs_order, unique(ge_df$chr))
  if(length(chrs_order)==0){
    stop('Inconsistent chromosome names in "chrs_order".')
  }
  cells.ss = lapply(unique(comm_df$community), function(comm){
    cells = comm_df$cell[which(comm_df$community == comm)]
    if(length(cells)>max_cells){
      cells = sample(cells, max_cells)
    }
    cells
  })
  cells.ss = do.call(c, cells.ss)
  ge_df = ge_df[,c('chr', 'start', cells.ss)]
  ge.chr = tidyr::gather(ge_df, 'cell', 'exp', 3:ncol(ge_df))
  ge.chr = merge(ge.chr, comm_df)
  ge.chr = ge.chr %>% dplyr::group_by(.data$community, .data$chr, .data$cell) %>%
    dplyr::summarize(exp=stats::median(.data$exp))
  norm.ge = ge.chr %>% dplyr::filter(.data$cell %in% baseline_cells) %>%
    dplyr::ungroup(.data) %>% 
    dplyr::group_by(.data$chr) %>%
    dplyr::summarize(exp.norm=stats::median(.data$exp))
  ge.chr = merge(ge.chr, norm.ge)
  ge.chr$exp = ge.chr$exp / ge.chr$exp.norm

  ## Reorder community if it looks like numeric
  if(all(!is.na(suppressWarnings(as.numeric(unique(ge.chr$community)))))){
    ge.chr$community = as.numeric(ge.chr$community)
  }
  
  community = exp = NULL
  ggp = lapply(chrs_order, function(ch){
    ge.chr %>% dplyr::filter(.data$chr==ch) %>%
      ggplot(aes(x=factor(community), y=winsor(exp, 2))) +
      geom_hline(yintercept=1) + 
      geom_hline(yintercept=c(.5,1.5), linetype=2) + 
      geom_boxplot() + theme_bw() + xlab('community') +
      ylab('median normalized chromosomal expression') +
      scale_y_continuous(breaks=seq(0, 2, .25),
                         labels=c(seq(0, 1.75, .25), '2+'),
                         limits=c(0,2)) +
      ggtitle(ch)
  })
  return(ggp)
}
