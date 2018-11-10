##' From raw gene expression, a few QC metrics are computed.
##'
##' If cell_cycle is provided it should be a data.frame (or a tsv file) with
##' two columns: 'symbol' with gene names, and 'phase' with the cell cycle phase
##' (e.g. either 'G1.S' or 'G2.M').
##' @title Compute quality control metrics for each cell
##' @param ge_df the input gene expression with a 'symbol' column and then one
##' column per cell.
##' @param cell_cycle if non-null, either a file or data.frame to compute cell
##' cycle scores. See details.
##' @return a data.frame with qc metrics per cell.
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @export
qc_cells <- function(ge_df, cell_cycle=NULL){
  options('dplyr.show_progress'=FALSE)
  cells = setdiff(colnames(ge_df), 'symbol')
  tot = colSums(ge_df[,cells])
  zeros = apply(ge_df[,cells], 2, function(x)sum(x==0))
  mito.idx = grep('MT-', ge_df$symbol)
  mito = rep(0, length(cells))
  if(length(mito.idx) > 0){
    mito = colSums(ge_df[mito.idx,cells])
  }
  qc.df = data.frame(tot=as.numeric(tot),
                     zeros=as.numeric(zeros),
                     mito=as.numeric(mito),
                     cell=cells,
                     stringsAsFactors=FALSE)
  ## Cell cycle score
  if(!is.null(cell_cycle)){
    if(!is.data.frame(cell_cycle)){
      if(is.character(cell_cycle) & length(cell_cycle) == 1){
        cell_cycle = utils::read.table(cell_cycle, as.is=TRUE, header=TRUE)
      } else {
        stop('cell_cycle must be either a file name or a data.frame.')
      }
    }
    ## For each phase compute the normalized expression of genes in each cell
    med.tot = stats::median(qc.df$tot)
    ge.cc = ge_df %>% dplyr::filter(.data$symbol %in% cell_cycle$symbol) %>%
      tidyr::gather('cell', 'ge', -'symbol') %>%
      merge(cell_cycle) %>%
      merge(qc.df) %>% dplyr::group_by(.data$symbol) %>%
      dplyr::mutate(ge.mean=mean(.data$ge*med.tot/.data$tot)) %>%
      dplyr::filter(.data$ge.mean>0) %>% dplyr::group_by(.data$phase, .data$cell) %>%
      dplyr::summarize(ge=mean(log(.data$ge+1)))
    ## Scale and add columns for each phase to the main df
    ge.cc = ge.cc %>% dplyr::ungroup(.data) %>% dplyr::group_by(.data$phase) %>%
      dplyr::mutate(ge=scale(.data$ge)) %>% tidyr::spread('phase', 'ge')
    if(nrow(ge.cc) > 0){
      qc.df = merge(qc.df, ge.cc)
    } else {
      warning('No cell-cycle gene expressed, suspicious...')
      qc.df$G1.S = qc.df$G2.M = NA
    }
  }
  return(qc.df)
}
