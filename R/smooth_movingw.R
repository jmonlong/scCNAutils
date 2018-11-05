##' @title Moving-window smoothing
##' @param df the input data.frame with coordinate columns (chr, start, end) and then
##' one column per cell
##' @param wsize the window size. Default is 3.
##' @param nb_cores the number of processors to use.
##' @param FUN the function to apply to each window. Default is median.
##' @return a data.frame with smoothed signal.
##' @author Jean Monlong
##' @export
smooth_movingw <- function(df, wsize=3, nb_cores=1, FUN=stats::median){
  ## Function to smooth on cell (vector)
  smooth.cell <- function(vec){
    if(length(vec) < wsize){
      return(rep(FUN(vec), length(vec)))
    }
    mmm = lapply(1:wsize, function(ii){
      c(rep(NA, ii-1), vec, rep(NA, wsize-ii))
    })
    mmm = matrix(unlist(mmm), ncol=wsize)
    mmm = mmm[round(wsize/2):(round(wsize/2)+length(vec)-1),]
    apply(mmm, 1, FUN, na.rm=TRUE)
  }
  ## function to smooth all cells
  smooth.f <- function(df){
    ## reorder just in case
    df = df[order(df$chr, df$start),]
    ## analyze each cell
    df[, cells] = apply(df[, cells], 2, smooth.cell)
    df
  }
  cells = setdiff(colnames(df), c("chr","start","end"))
  ## Smooth each chromosome
  res.l = parallel::mclapply(unique(df$chr), function(chri){
    smooth.f(df[which(df$chr == chri),])
  }, mc.cores=nb_cores)
  return(as.data.frame(data.table::rbindlist(res.l)))
}
