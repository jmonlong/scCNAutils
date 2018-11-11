##' @title Metacells by subsampling within communities.
##' @param groups.df the filtered data.frame with community information
##' @param nb_metacells the number of metacells per comunity.
##' @param metacell_size the number of cells in a metacell.
##' @return a list with
##' \item{cells.tm}{a list with cell names for each metacell.}
##' \item{cells.info}{a data.frame recording the community of origina for each metacell.}
##' @author Jean Monlong
##' @keywords internal
metacells_subsample <- function(groups.df, nb_metacells, metacell_size){
  ## Select cells for metacells
  cells.tm = lapply(unique(groups.df$community), function(comm){
    group.comm = groups.df[which(groups.df$community==comm),]
    recycle = FALSE
    if(nrow(group.comm) < metacell_size*nb_metacells) {
      warning("Recycling cells because ", comm, " is a small group.")
      recycle = TRUE
    }
    cells.s = sample(group.comm$cell, metacell_size*nb_metacells, recycle)
    split(cells.s, rep(1:nb_metacells, metacell_size))
  })
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
