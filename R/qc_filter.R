##' @title Filter cells based on QC results
##' @param ge_df the input gene expression with a 'symbol' column and then one
##' column per cell.
##' @param qc_df the output data.frame from qc_cells
##' @param max_mito_prop the maximum proportion of mitochondrial RNA.
##' @param min_total_exp the minimum total cell expression
##' @return \emph{ge_df} with only the cells that passed the filters
##' @author Jean Monlong
qc_filter <- function(ge_df, qc_df, max_mito_prop=.2, min_total_exp=0){
  ## Find cells to keep
  mito.filt = qc_df$mito / qc_df$tot <= max_mito_prop
  exp.filt = qc_df$tot >= min_total_exp
  cells.tokeep = qc_df$cell[which(mito.filt & exp.filt)]
  ## Keep cells
  info.cols = intersect(colnames(ge_df), c("chr","start","end", 'symbol'))
  ge_df[,c(info.cols, intersect(colnames(ge_df), cells.tokeep))]
}
