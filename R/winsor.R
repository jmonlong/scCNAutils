##' Convenience function to winsorize a vector.
##' @title Winsorize
##' @param x input vector
##' @param u upper limit
##' @param l lower limit
##' @param uq the quantile for the upper limit. Used is u is NULL.
##' @return winsorized vector
##' @author Jean Monlong
##' @export
winsor <- function(x, u=NULL, l=NULL, uq=NULL){
  if(is.null(u) & !is.null(uq)) u = stats::quantile(x, probs=uq, na.rm=TRUE)
  if(!is.null(u) && any(x>u, na.rm=TRUE)) x[which(x>u)] = u
  if(!is.null(l) && any(x<l, na.rm=TRUE)) x[which(x<l)] = l
  x
}
