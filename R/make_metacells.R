##' Randomly select cells in each community and merge them to create metacells
##' with higher resolution.
##'
##' @title Make metacells
##' @param ge_df gene expression of all cells
##' @param comm_df a data.frame with community information, output from
##' \code{\link{find_communities}}.
##' @param nb_metacells the number of metacells per comunity.
##' @param metacell_size the number of cells in a metacell.
##' @param baseline_cells the cells to use for baseline communities.
##' @param nb_cores the number of processor to use.
##' @param max_baseline_comm the maximum number of baseline communities to generate.
##' @return a list with
##' \item{ge}{a data.frame with coordinates and gene expression for each metacell.}
##' \item{info}{information about which metacell correspond to which community.}
##' \item{mc_cells}{information about which cells were used for each metacell.}
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @export
make_metacells <- function(ge_df, comm_df, nb_metacells=10, metacell_size=3,
                           baseline_cells=NULL, nb_cores=1, max_baseline_comm=3){
  coords.symbol = intersect(colnames(ge_df), c("chr","start","end", 'symbol'))
  cells = setdiff(colnames(ge_df), c("chr","start","end", 'symbol'))
  options('dplyr.show_progress'=FALSE)
  if(is.null(baseline_cells)){
    baseline_cells = cells
  } else if(is.logical(baseline_cells) && !baseline_cells){
    baseline_cells = NULL
  }

  ## Add "baseline" communities
  groups.df = NULL
  if(is.null(baseline_cells) | is.null(nb_metacells)){
    groups.df = comm_df[, c('cell', 'community')]
  } else {
    groups.df = data.frame(cell=baseline_cells, stringsAsFactors=FALSE)
    nbm = ifelse(is.null(nb_metacells), 1, nb_metacells)
    baseline.comm = sample.int(min(max_baseline_comm,
                                   ceiling(length(baseline_cells)/
                                           (2*metacell_size*nbm))),
                               length(baseline_cells), replace=TRUE)
    groups.df$community = paste0('baseline_', baseline.comm)
    groups.df = rbind(groups.df, comm_df[, c('cell', 'community')])
  }

  ## Filter small communities
  comm.size = table(groups.df$community)
  small.groups = names(comm.size)[which(comm.size < metacell_size)]
  if(length(small.groups)>0){
    warning('Communities smaller than the number of cells to merge have been removed')
    groups.df = groups.df[which(!(groups.df$community %in% small.groups)), ]
  }

  if(!is.null(nb_metacells)){    
    mc.o = metacells_subsample(groups.df, nb_metacells, metacell_size)
  } else {
    mc.o = metacells_cluster(groups.df, ge_df, metacell_size, nb_cores)
  }
  
  ## Merge gene expression
  ge.mc = ge_df[, coords.symbol, drop=FALSE]
  ge.mc.l = parallel::mclapply(1:length(mc.o$cells.tm), function(ii){
    ##message(ii)
    ge_df[,mc.o$cells.tm[[ii]]] %>% as.matrix %>% rowSums
  }, mc.cores=nb_cores)
  names(ge.mc.l) = mc.o$cells.info$cell
  ge.mc = cbind(ge.mc, as.data.frame(ge.mc.l))

  ## Format cells info into data.frame
  mc_cells = data.frame(cell=unlist(mc.o$cells.tm),
                        metacell=rep(mc.o$cells.info$cell, unlist(lapply(mc.o$cells.tm, length))),
                        stringsAsFactors=FALSE)
  
  return(list(ge=ge.mc, info=mc.o$cells.info, mc_cells=mc_cells))
}
