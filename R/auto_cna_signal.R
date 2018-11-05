##' Goes from reading raw gene counts to CNA-level signal, tSNE and community detection.
##' @title Automated pipeline to compute CNA signal from scRNA expression
##' @param data a data.frame with gene expression of the path to the folder with the
##' 'matrix.mtx', 'genes.tsv' and 'barcodes.tsv' files.
##' @param genes_coord either a file name or a data.frame with coordinates
##' and gene names.
##' @param prefix the prefix to use for the files created by this function (e.g. graphs).
##' @param nb_cores the number of processors to use.
##' @param pause_after_qc pause after the QC to pick custom QC thresholds.
##' @param max_mito_prop the maximum proportion of mitochondrial RNA.
##' @param min_total_exp the minimum total cell expression
##' @param chrs the chromosome names to keep. NULL to include all the chromosomes.
##' @param cell_cycle if non-null, either a file or data.frame to compute cell
##' cycle scores. See details.
##' @param bin_mean_exp the desired minimum mean expression in the bin.
##' @param z_wins_th the threshold to winsorize Z-score. Default is 3
##' @param smooth_wsize the window size for smoothing. Default is 3.
##' @param cc_sd_th the number of SD used for the thresholds when defining cycling cells.
##' @param nb_pcs the number of PCs used in the community detection or tSNE. 
##' @param comm_k the number of nearest neighbor for the KNN graph. Default 100.
##' @return a data.frame with QC, community and tSNE for each cell.
##' @author Jean Monlong
##' @export
auto_cna_signal <- function(data, genes_coord, prefix='scCNAutils_out', nb_cores=1,
                            pause_after_qc=FALSE,
                            max_mito_prop=.2, min_total_exp=0, chrs=c(1:22,"X","Y"),
                            cell_cycle=NULL, bin_mean_exp=3, z_wins_th=3, smooth_wsize=3,
                            cc_sd_th=3, nb_pcs=10, comm_k=100){
  if(is.character(data) & length(data)==1){
    data = read_mtx(path=data)
  }

  ## QC
  qc.df = qc_cells(data, cell_cycle)
  qc.ggp = plot_qc_cells(qc.df)
  grDevices::pdf(paste0(prefix, '-qc.pdf'), 9, 7)
  lapply(qc.ggp, print)
  grDevices::dev.off()
  if(pause_after_qc){
    message('Inspect ', paste0(prefix, '-qc.pdf'), ' and rerun with "pause_after_qc=FALSE", eventually adjusting "max_mito_prop" and "min_total_exp".')
    return(list(qc=qc.df))
  }
  data = qc_filter(data, qc.df, max_mito_prop, min_total_exp)
  
  ## Expression to CNA signal
  data = convert_to_coord(data, genes_coord, chrs)
  data = norm_ge(data, nb_cores=nb_cores)
  save(data, file=paste0(prefix, '-coord-norm.RData'))
  data = bin_genes(data, bin_mean_exp, nb_cores)
  data = zscore(data, z_wins_th, method='z')
  data = smooth_movingw(data, smooth_wsize, nb_cores)
  save(data, file=paste0(prefix, '-coord-norm-bin', bin_mean_exp, '-z', z_wins_th,
                         '-smooth', smooth_wsize, '.RData'))

  ## Potentially find cycling cells
  core_cells = NULL
  if(!is.null(cell_cycle)){
    cc.l = define_cycling_cells(qc.df, cc_sd_th)
    core_cells = cc.l$cells.noc
    grDevices::pdf(paste0(prefix, '-qc-cellcycle.pdf'), 9, 7)
    lapply(cc.l$graphs, print)
    grDevices::dev.off()
  }
  
  ## PCA
  pca.o = run_pca(data, core_cells)
  save(pca.o, file=paste0(prefix, '-coord-norm-bin', bin_mean_exp, '-z', z_wins_th,
                         '-smooth', smooth_wsize, '-pca.RData'))
  grDevices::pdf(paste0(prefix, '-coord-norm-bin', bin_mean_exp, '-z', z_wins_th,
                        '-smooth', smooth_wsize, '-pcaSD.pdf'), 9, 7)
  print(pca.o$sdev.graph)
  grDevices::dev.off()

  ## Community
  comm.df = find_communities(pca.o, nb_pcs, comm_k)
  comm.ggp = plot_communities(comm.df)
  grDevices::pdf(paste0(prefix, '-coord-norm-bin', bin_mean_exp, '-z', z_wins_th,
                        '-smooth', smooth_wsize, '-comm', nb_pcs, 'PCs', comm_k,
                        '.pdf'), 9, 7)
  lapply(comm.ggp, print)
  grDevices::dev.off()

  ## tSNE
  tsne.df = run_tsne(pca.o, nb_pcs)
  tsne.ggp = plot_tsne(tsne.df, qc.df)
  grDevices::pdf(paste0(prefix, '-coord-norm-bin', bin_mean_exp, '-z', z_wins_th,
                        '-smooth', smooth_wsize, '-tsne', nb_pcs, 'PCs.pdf'), 9, 7)
  lapply(tsne.ggp, print)
  grDevices::dev.off()

  ## Master data.frame with QC, community and tSNE info
  res.df = merge(qc.df, comm.df)
  res.df = merge(res.df, tsne.df)
  return(res.df)
}
