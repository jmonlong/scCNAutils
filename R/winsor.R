##' Convenience function to winsorize a vector.
##' @title Winsorize
##' @param x input vector
##' @param u upper limit
##' @param l lower limit
##' @return winsorized vector
##' @author Jean Monlong
##' @export
winsor <- function(x, u=NULL, l=NULL){
  if(!is.null(u) && any(x>u, na.rm=TRUE)) x[which(x>u)] = u
  if(!is.null(l) && any(x<l, na.rm=TRUE)) x[which(x<l)] = l
  x
}
