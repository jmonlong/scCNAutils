##' @title Winsorize
##' @param x input vector
##' @param u upper limit
##' @param l lower limit
##' @param qu upper quartile, if 'u' is null
##' @return winsorized vector
##' @author Jean Monlong
##' @keywords internal
winsor <- function(x, u=NULL, l=NULL, qu=NULL){
  if(is.null(u) && !is.null(qu)) u = stats::quantile(x, probs=qu)
  if(!is.null(u) && any(x>u)) x[x>u] = u
  if(!is.null(l) && any(x<l)) x[x<l] = l
  x
}
