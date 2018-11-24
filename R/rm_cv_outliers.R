##' Compute the coefficient of variation for each gene/bin and remove the ones
##' with the highest values, either based on a quantile or SD-based threshold.
##' The genes/bins that satisfy both quantile and SD-based thresholds are
##' removed.
##' 
##' @title Remove outliers based on the coefficient of variation
##' @param ge_df the expression data.frame.
##' @param ol_quant_th the quantile threshold. Default is 0.99 (removes the
##' top 1\% with highest values).
##' @param ol_sd_th the SD-based treshold. Default is 5.
##' @return a subset of the ge_df data.frame
##' @author Jean Monlong
##' @export
rm_cv_outliers <- function(ge_df, ol_quant_th=.99, ol_sd_th=5){
  cells = setdiff(colnames(ge_df), c("chr","start","end", "symbol"))
  ## Compute coefficient of variation
  bin.m = apply(ge_df[,cells], 1, mean, na.rm=TRUE)
  bin.sd = apply(ge_df[,cells], 1, stats::sd, na.rm=TRUE)
  bin.cv = bin.sd / bin.m

  ## Quantile-based threshold
  quant.ol = bin.cv > stats::quantile(bin.cv, probs=ol_quant_th)

  ## MAD-based threshold
  cv.m = mean(bin.cv)
  cv.mad = stats::mad(bin.cv)
  mad.ol = bin.cv > cv.m + ol_sd_th * cv.mad

  return(ge_df[which(!(quant.ol & mad.ol)),])
}
