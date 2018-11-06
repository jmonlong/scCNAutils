##' The expression of each cell is normalized to account for depth differences.
##'
##' @title Normalize gene expression
##' @param ge_df the input gene expression
##' @param method the normalization method
##' @param nb_cores the number of processors to use.
##' @return a data.frame with the normalized expression.
##' @author Jean Monlong
##' @export
norm_ge <- function(ge_df, method=c('tmm', 'total'), nb_cores=1){
  if(method[1] == 'tmm'){
    return(norm_ge_tmm(ge_df=ge_df, nb_cores=nb_cores))
  } else if(method[1] == 'total') {
    return(norm_ge_total(ge_df=ge_df))
  } else {
    stop('method must be one of: tmm, total')
  }
}
    
