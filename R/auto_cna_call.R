##' Automated pipeline to call CNA using metacells.
##' 
##' Once the metacells are created there are two ways to call CNA. First, if
##' \code{multisamps=FALSE}, to call CNA on each metacell and merge the result
##' per community, keeping the information about how many metacell support the
##' CNA. Second, if \code{multisamps=TRUE} (default), to run the HMM on all the
##' metacells for a community. The multi-sample approach should be more robust.
##'
##' The transition probability (\code{trans_prob}) is going to affect the HMM
##' segmentation. Smaller values will create longer segments. One approach,
##' often advocated by HMM aficionados, is to try different values and use the
##' ones that gives the best results, for example based on the QC graphs (TODO).
##' Another approach is to use a loose transition probability and then filter
##' short segments ('length' column or 'pass.filter' column).
##' @title Automated pipeline to call CNA
##' @param ge_df normalized gene expression of all cells (e.g. output from
##' \code{\link{norm_ge}}.
##' @param comm_df a data.frame with community information, output from
##' \code{\link{find_communities}}.
##' @param nb_metacells the number of metacells per comunity.
##' @param metacell_size the number of cells in a metacell.
##' @param multisamps use the multi-sample version of the HMM segmentation? Default is TRUE.
##' See details.
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
                          multisamps=TRUE, trans_prob=.1,
                          baseline_cells=NULL, baseline_communities=NULL,
                          prefix='scCNAutils_out',
                          nb_cores=1, chrs=c(1:22,"X","Y"),
                          bin_mean_exp=3, z_wins_th=3, smooth_wsize=3){
  comm_df = comm_df[which(comm_df$cell %in% colnames(ge_df)),]
  if(is.null(baseline_cells) & !is.null(baseline_communities)){
    baseline_cells = comm_df$cell[which(comm_df$community %in% baseline_communities)]
  }

  ## Aneuploidy graph from cells
  ggp = plot_aneuploidy(ge_df, comm_df, baseline_cells=baseline_cells)
  grDevices::pdf(paste0(prefix, '-aneuploidy.pdf'), 9, 7)
  tmp = lapply(ggp, print)
  grDevices::dev.off()

  ## Make metacells
  message('Creating metacells...')
  mc.o = make_metacells(ge_df, comm_df, nb_metacells, metacell_size, baseline_cells,
                        nb_cores)
  save(mc.o, file=paste0(prefix, '-metacells.RData'))

  ## Normalize and bin genes
  message('Normalizing and binning...')
  ge_df = norm_ge(mc.o$ge, nb_cores=nb_cores)
  ge_df = bin_genes(ge_df, bin_mean_exp, nb_cores=nb_cores)
  ge_df = norm_ge(ge_df, nb_cores=nb_cores)

  ## Aneuploidy graph from metacells
  baseline.metacells = mc.o$info$cell[grep('baseline_', mc.o$info$community)]
  info.nobaseline = mc.o$info
  if(!is.null(baseline_communities)){
    info.nobaseline = mc.o$info[grep('baseline_', mc.o$info$community, invert=TRUE), ]
    baseline.metacells = mc.o$info$cell[which(mc.o$info$community %in% baseline_communities)]
  }
  ggp = plot_aneuploidy(ge_df, info.nobaseline, baseline_cells=baseline.metacells)
  grDevices::pdf(paste0(prefix, '-aneuploidy-metacells.pdf'), 9, 7)
  tmp = lapply(ggp, print)
  grDevices::dev.off()
    
  ## Normalize using baseline metacells
  message('Scaling...')
  baseline.metacells = mc.o$info$cell[grep('baseline_', mc.o$info$community)]
  z_df = zscore(ge_df, z_wins_th, method='norm', normals=baseline.metacells)

  ## Remove baseline metacells and smooth
  message('Smoothing...')
  comm.metacells = mc.o$info$cell[grep('baseline_', mc.o$info$community, invert=TRUE)]
  z_df = z_df[, c('chr','start','end', comm.metacells)]
  z_df = smooth_movingw(z_df, smooth_wsize, nb_cores)
  save(z_df, file=paste0(prefix, '-cnaz.RData'))

  ## Call CNAs and generate graph
  message('Calling CNAs...')
  if(multisamps){
    cna.o = call_cna_multisamps(z_df, mc.o$info, trans_prob, nb_cores)
  } else {
    cna.o = call_cna(z_df, trans_prob, nb_cores, mc.o$info)
  }
  ggp = plot_cna(cna.o, chrs)
  grDevices::pdf(paste0(prefix, '-CNA.pdf'), 9, 7)
  tmp = lapply(ggp, print)
  grDevices::dev.off()

  return(cna.o)
}
