##' @title Normalize gene expression using total expression
##' @param ge_df the input gene expression
##' @param nb.cores the number of processors to use.
##' @return a data.frame with the normalized expression.
##' @author Jean Monlong
##' @keywords internal
norm_ge_total <- function(ge_df, nb.cores=1){
  cells = setdiff(colnames(ge_df), c("chr","start","end", 'symbol'))
  total.reads = colSums(ge_df[,cells])
  total.reads.med = stats::median(total.reads)
  ge_df[, cells] = t(t(ge_df[, cells]) / total.reads * total.reads.med)
  return(ge_df)
}
