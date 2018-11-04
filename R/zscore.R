##' @title Compute Z-score
##' @param ge_df the inpue expression data.frame
##' @param wins.th the threshold to winsorize Z-score. Default is 3
##' @param method the normalization method. Either norm or z.
##' @param normals the cells to use as normals. If NULL (default) all cells are used as normals
##' @return a data.frame with Z-scores.
##' @author Jean Monlong
zscore <- function(ge_df, wins.th=3, method=c('norm', 'z'), normals=NULL){
  samps = setdiff(colnames(ge_df), c("chr","start","end", "symbol"))
  if(is.null(normals)){
    normals = 1:length(samps)
  } else {
    normals = which(samps %in% normals)
  }
  if(method[1]=='norm'){
    ge_df[,samps] = t(apply(ge_df[,samps], 1, function(x) winsor(log(x/max(mean(x[normals], na.rm=TRUE),1)), u=wins.th, l=-wins.th)))
  } else if(method[1]=='z'){
    ge_df[,samps] = t(apply(ge_df[,samps], 1, function(x) winsor((x-mean(x[normals], na.rm=TRUE))/max(stats::sd(x, na.rm=TRUE),1), u=wins.th, l=-wins.th)))
  }
  ge_df
}
