##' @title Automated pipeline to call CNA
##' @param ge_df normalized gene expression of all cells (e.g. output from
##' \code{\link{norm_ge}}.
##' @param comm_df a data.frame with community information, output from
##' \code{\link{find_communities}}.
##' @param nb_metacells the number of metacells per comunity.
##' @param metacell_size the number of cells in a metacell.
##' @param trans_prob the transition probability for the HMM.
##' @param baseline_cells cells to use as baseline.
##' @param baseline_communities communities to use as baseline. Used if
##' baseline.cells is NULL.
##' @param prefix the prefix to use for the files created by this function (e.g. graphs).
##' @param nb_cores the number of processors to use.
##' @param chrs the chromosome names to keep. NULL to include all the chromosomes.
##' @param bin_mean_exp the desired minimum mean expression in the bin.
##' @param z_wins_th the threshold to winsorize Z-score. Default is 3
##' @param smooth_wsize the window size for smoothing. Default is 3.
##' @return a data.frame with CNAs
##' @author Jean Monlong
##' @export
auto_cna_call <- function(ge_df, comm_df, nb_metacells=10, metacell_size=3,
                          trans_prob=1e-4,
                          baseline_cells=NULL, baseline_communities=NULL,
                          prefix='scCNAutils_out',
                          nb_cores=1, chrs=c(1:22,"X","Y"),
                          bin_mean_exp=3, z_wins_th=3, smooth_wsize=3){

  ## Make metacells
  message('Creating metacells...')
  if(is.null(baseline_cells) & !is.null(baseline_communities)){
    baseline_cells = comm_df$cell[which(comm_df$community %in% baseline_communities)]
  }
  mc.o = make_metacells(ge_df, comm_df, nb_metacells, metacell_size, baseline_cells,
                        nb_cores)

  ## Normalize and bin genes
  message('Normalizing and binning...')
  ge_df = norm_ge(mc.o$ge, nb_cores=nb_cores)
  ge_df = bin_genes(ge_df, bin_mean_exp, nb_cores=nb_cores)

  ## Normalize using baseline metacells
  message('Scaling...')
  baseline.metacells = mc.o$info$cell[grep('baseline_', mc.o$info$community)]
  z_df = zscore(ge_df, z_wins_th, method='norm', normals=baseline.metacells)

  ## Remove baseline metacells and smooth
  message('Smoothing...')
  comm.metacells = mc.o$info$cell[grep('baseline_', mc.o$info$community, invert=TRUE)]
  z_df = z_df[, c('chr','start','end', comm.metacells)]
  z_df = smooth_movingw(z_df, smooth_wsize, nb_cores)

  ## Call CNAs and generate graph
  message('Calling CNAs...')
  cna.df = call_cna(z_df, trans_prob, nb_cores, mc.o$info)
  ggp = plot_cna(cna.df, chrs)
  grDevices::pdf(paste0(prefix, '-CNA.pdf', 9, 7))
  print(ggp)
  grDevices::dev.off()
  
  return(cna.df)
}
