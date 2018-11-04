##' @title Normalize gene expression using TMM
##' @param ge_df the input gene expression.
##' @param nb_cores the number of processors to use.
##' @return a data.frame with the normalized expression.
##' @author Jean Monlong
##' @keywords internal
norm_ge_tmm <- function(ge_df, nb_cores=1){
  tmmNorm <- function(ge.samp, ge.cont, trim.level=.3){
    no0 = which(ge.samp!=0 & ge.cont!=0)
    m = exp(mean(log(ge.samp[no0]/ge.cont[no0]), trim=trim.level))
    ge.samp / m
  }
  cells = setdiff(colnames(ge_df), c("chr","start","end", 'symbol'))
  ## Pick control sample with total expression close to median.
  total.reads = colSums(as.matrix(ge_df[,cells]))
  cont.sample = cells[order(total.reads)[length(cells)/2]]
  cont.ge = as.numeric(ge_df[,cont.sample])
  info.cols = intersect(colnames(ge_df), c("chr","start","end", 'symbol'))
  norm.df = ge_df[, info.cols, drop=FALSE]
  ## Chunks for parallel computation
  if(nb_cores == 1){
    cells.chunks = list(x=cells)
  } else {
    cells.chunks = tapply(cells, cut(1:length(cells), nb_cores), identity)
  }
  ge.mat.l = parallel::mclapply(cells.chunks, function(cells){
    apply(ge_df[, cells], 2, tmmNorm, cont.ge)
  }, mc.cores=nb_cores)
  ## Merge chunks
  norm.mat = do.call(cbind, ge.mat.l)
  return(cbind(norm.df, norm.mat))
}
