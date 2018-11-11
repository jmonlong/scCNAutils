##' @title Metacells by clustering within communities.
##' @param groups.df the filtered data.frame with community information
##' @param ge_df normalized gene expression of all cells (e.g. output from
##' \code{\link{norm_ge}}.
##' @param nb_metacells the number of metacells per comunity.
##' @param metacell_size the number of cells in a metacell.
##' @param nb_cores the number of processor to use.
##' @param hclust_method the 'method' used when clustering (see hclust).
##' @return a list with
##' \item{cells.tm}{a list with cell names for each metacell.}
##' \item{cells.info}{a data.frame recording the community of origina for each metacell.}
##' @author Jean Monlong
##' @keywords internal
metacells_cluster <- function(groups.df, ge_df, nb_metacells, metacell_size, nb_cores=1, 
                              hclust_method='ward.D'){
  cells.tm = parallel::mclapply(unique(groups.df$community), function(comm){
    group.comm = groups.df[which(groups.df$community==comm),]
    mat = t(as.matrix(ge_df[, group.comm$cell]))
    dmat = as.matrix(stats::dist(mat))
    clsize.rle = rle(as.numeric(cut(1:nrow(mat), ceiling(nrow(mat)/metacell_size))))
    clsizes = clsize.rle$lengths
    cpt = 1
    lab = rep(NA, nrow(mat))
    for(clss in clsizes[-1]){
      lab.ii = which(is.na(lab))
      hc.o = stats::hclust(stats::as.dist(dmat[lab.ii, lab.ii]), method=hclust_method)
      clt = 0
      ct = length(lab.ii)-clss
      while(max(clt)<clss){
        cls = stats::cutree(hc.o, ct)
        clt = table(cls)
        ct = ct - 1
      }
      cl.sel = which(cls == as.numeric(names(clt)[which.max(clt)]))
      lab[lab.ii[utils::head(cl.sel, clss)]] = cpt
      cpt = cpt + 1
    }
    lab[which(is.na(lab))] = cpt
    split(group.comm$cell, lab)
  }, mc.cores=nb_cores)
  cells.info = data.frame(community=rep(unique(groups.df$community),
                                        unlist(lapply(cells.tm, length))),
                          stringsAsFactors=FALSE)
  if(is.factor(groups.df$community)){
    cells.info$community = factor(cells.info$community,
                                  levels=levels(groups.df$community))
  }
  cells.info$cell = paste0("mc",1:nrow(cells.info))
  return(list(cells.tm=do.call(c, cells.tm), cells.info=cells.info))
}
