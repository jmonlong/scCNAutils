##' PCA analysis, eventually using a subset of core cells for the PC construction.
##'
##' Cells in core_cells are used to build the principal components to which
##' all cells are then projected to. Usually used to reduce the effect of
##' cell cycle in the PCA, by using only cells that don't cycle (see \code{\link{qc_cells}})
##' as \emph{core_cells}.
##'
##' The graph (\emph{sdev.graph}) shows the standard deviation for the top 50 PCs. To
##' show more/less PCs, add \code{xlim(1,N)} to the \emph{sdev.graph}. See examples.
##' @title Run PCA
##' @param z_df a data.frame with z-scores for each cell
##' @param core_cells if non-NULL, a vector with the names of the cells to use as core
##' cells. See details. Default is NULL.
##' @param out_pcs the number of top PCs to report. Default is 100.
##' @return a list with
##' \item{x}{the PC matrix}
##' \item{sdev}{the standard deviations of the PCs}
##' \item{sdev.graph}{a ggplot graph of the sdev}
##' @author Jean Monlong
##' @import ggplot2
##' @examples
##' \dontrun{
##' pca.o = run_pca(z)
##' 
##' ## Zoom in to the top 20 PCs
##' pca.o$sdev.graph + xlim(1,20)
##' }
##' @export
run_pca <- function(z_df, core_cells=NULL, out_pcs=100){
  cells = setdiff(colnames(z_df), c("chr","start","end", "symbol"))
  z.mat = t(as.matrix(z_df[, cells]))
  if(!is.null(core_cells)){
    ## build with core cells and project all cells
    cells.ol = intersect(core_cells, rownames(z.mat))
    if(length(cells.ol) == 0){
      stop('None of core_cells are in z_df.')
    } else if(length(cells.ol) < length(core_cells)){
      warning('Some core_cells are not in z_df.')
    } 
    pca.noc = stats::prcomp(z.mat[cells.ol,], scale=TRUE)
    pca.o = list(x=scale(z.mat) %*% pca.noc$rotation, sdev=pca.noc$sdev)
  } else {
    # build with all cells
    pca.o = stats::prcomp(z.mat, scale=TRUE)
    pca.o = list(x=pca.o$x, sdev=pca.o$sdev)
  }
  ## keep only top PCs
  out_pcs = min(out_pcs, ncol(pca.o$x))
  pca.o$x = pca.o$x[, 1:out_pcs]
  ## graph
  PC = sdev = type = NULL
  sd.df = rbind(data.frame(PC=1:length(pca.o$sdev), sdev=pca.o$sdev, type='sdev'),
                data.frame(PC=1:length(pca.o$sdev), sdev=cumsum(pca.o$sdev),
                           type='cumulative sdev'))
  pca.o$sdev.graph = ggplot(sd.df[which(sd.df$PC<=50),], aes(x=PC, y=sdev, colour=type)) + geom_point(na.rm=TRUE) +
    geom_line(na.rm=TRUE) + theme_bw() + ylab('standard deviation') +
    theme(legend.title=element_blank(), legend.position=c(.01,.99),
          legend.justification=c(0,1)) 
  return(pca.o)
}
