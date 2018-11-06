##' @title Filter cells based on QC results
##' @param ge_df the input gene expression with a 'symbol' column and then one
##' column per cell.
##' @param qc_df the output data.frame from qc_cells
##' @param max_mito_prop the maximum proportion of mitochondrial RNA.
##' @param min_total_exp the minimum total cell expression
##' @param cells_sel consider only these cells. Other cells filtered no matter what.
##' @return \emph{ge_df} with only the cells that passed the filters
##' @author Jean Monlong
##' @export
qc_filter <- function(ge_df, qc_df, max_mito_prop=.2, min_total_exp=0, cells_sel=NULL){
  ## Find cells to keep
  mito.filt = qc_df$mito / qc_df$tot <= max_mito_prop
  exp.filt = qc_df$tot >= min_total_exp
  cells.tokeep = qc_df$cell[which(mito.filt & exp.filt)]
  if(!is.null(cells_sel)){
    cells.tokeep = intersect(cells.tokeep, cells_sel)
  }
  if(length(cells.tokeep)==0){
    stop("No cells match the filtering criteria.")
  }
  ## Keep cells
  info.cols = intersect(colnames(ge_df), c("chr","start","end", 'symbol'))
  ge_df[,c(info.cols, intersect(colnames(ge_df), cells.tokeep))]
}
