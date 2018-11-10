##' Randomly select cells in each community and merge them to create metacells
##' with higher resolution.
##'
##' @title Make metacells
##' @param ge_df normalized gene expression of all cells (e.g. output from
##' \code{\link{norm_ge}}.
##' @param comm_df a data.frame with community information, output from
##' \code{\link{find_communities}}.
##' @param nb_metacells the number of metacells per comunity.
##' @param metacell_size the number of cells in a metacell.
##' @param baseline_cells the cells to use for baseline communities.
##' @param nb_cores the number of processor to use.
##' @return a list with
##' \item{ge}{a data.frame with coordinates and gene expression for each metacell.}
##' \item{info}{information about which metacell correspond to which community.}
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @export
make_metacells <- function(ge_df, comm_df, nb_metacells=10, metacell_size=3,
                           baseline_cells=NULL, nb_cores=1){
  cells = setdiff(colnames(ge_df), c("chr","start","end"))
  options('dplyr.show_progress'=FALSE)
  if(is.null(baseline_cells)){
    baseline_cells = cells
  }

  ## Add "baseline" communities
  groups.df = data.frame(cell=baseline_cells, stringsAsFactors=FALSE)
  baseline.comm = sample.int(ceiling(length(baseline_cells)/
                                     (2*metacell_size*nb_metacells)),
                             length(baseline_cells), replace=TRUE)
  groups.df$community = paste0('baseline_', baseline.comm)
  groups.df = rbind(groups.df, comm_df[, c('cell', 'community')])

  ## Filter small communities
  comm.size = table(groups.df$community)
  small.groups = names(comm.size)[which(comm.size < metacell_size)]
  if(length(small.groups)>0){
    warning('Communities smaller than the number of cells to merge have been removed')
    groups.df = groups.df[which(!(groups.df$community %in% small.groups)), ]
  }

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
  cells.tm = do.call(c, cells.tm)

  ## Merge gene expression
  ge.mc = ge_df[, c('chr', 'start', 'end')]
  ge.mc.l = parallel::mclapply(1:length(cells.tm), function(ii){
    ge_df[,cells.tm[[ii]]] %>% as.matrix %>% apply(1,sum)
  }, mc.cores=nb_cores)
  names(ge.mc.l) = cells.info$cell
  ge.mc = cbind(ge.mc, as.data.frame(ge.mc.l))
  return(list(ge=ge.mc, info=cells.info))
}
