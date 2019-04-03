##' @title Heatmap of the CNA scores
##' @param z_df the data.frame with Z-scores.
##' @param cells_df a data.frame with cell information.
##' @param nb_subsamp the number of cells to sub-sample. Default is 1000.
##' @param hc.method the hierarchical clustering method to order cells.
##' @param z_win threshold to winsorized Z scores in the color scale.
##' @return a list of ggplot2 objects.
##' @author Jean Monlong
##' @export
##' @import ggplot2
plot_heatmap <- function(z_df, cells_df=NULL, nb_subsamp=1000,
                         hc.method='ward.D', z_win=3){
  chr = start = end = cell = cello = z = NULL
  cells = setdiff(colnames(z_df), c('chr', 'start', 'end'))

  ## Subsample Z scores
  cells = sample(cells, min(c(nb_subsamp, length(cells))))
  z_df = z_df[, c('chr', 'start', 'end', cells)]

  ## Hierarchical clustering to order cells.
  z.d = stats::dist(t(z_df[,cells]))
  hc.o = stats::hclust(z.d, method=hc.method)

  ## Melt score matrix
  z_df = tidyr::gather(z_df, cell, z, -chr, -start, -end)
  
  ## Merge cells info
  cells_df = cells_df[, intersect(colnames(cells_df), c('cell', 'community', 'sample'))]
  z_df = merge(z_df, cells_df, all.x=TRUE)

  ## Order for cells and chromosomes
  z_df$cell = factor(as.character(z_df$cell), levels=hc.o$labels[hc.o$order])
  if(all(grepl('chr', z_df$chr))){
    z_df$chr = factor(z_df$chr, levels=paste0('chr', c(1:22, 'X','Y')))
  } else {
    z_df$chr = factor(z_df$chr, levels=c(1:22, 'X','Y'))
  }

  ## ggplot2
  ggp.l = list()
  ggp.l$all = ggplot(z_df, aes(xmin=start, xmax=end, ymin=as.numeric(cell)-.5,
                               ymax=as.numeric(cell)+.5,
                               fill=winsor(z, z_win, -z_win))) +
    geom_rect() + theme_bw() +
    facet_grid(.~chr, scales='free', space='free') +
    ylab('cell') + xlab('position') +
    scale_fill_gradient2(name="z", low=scales::muted("blue"),
                         high=scales::muted("red")) + 
    theme(axis.text.y=element_blank(), axis.text.x=element_blank(), legend.position='bottom')

  ## Group by communities
  if('community' %in% colnames(z_df)){
    cell.order = unique(z_df$cell[order(z_df$community, z_df$cell)])
    z_df$cello = factor(z_df$cell, levels=cell.order)
    ggp.l$community = ggplot(z_df, aes(xmin=start, xmax=end,
                                       ymin=as.numeric(cello)-.5,
                                       ymax=as.numeric(cello)+.5,
                                       fill=winsor(z, z_win, -z_win))) +
    geom_rect() + theme_bw() +
    facet_grid(community~chr, scales='free', space='free') +
    ylab('cell') + xlab('position') + 
    scale_fill_gradient2(name="z", low=scales::muted("blue"),
                         high=scales::muted("red")) + 
    theme(axis.text.y=element_blank(), axis.text.x=element_blank(), legend.position='bottom')
  }

  ## Group by sample
  if('sample' %in% colnames(z_df)){
    cell.order = unique(z_df$cell[order(z_df$sample, z_df$cell)])
    z_df$cello = factor(z_df$cell, levels=cell.order)
    ggp.l$sample = ggplot(z_df, aes(xmin=start, xmax=end,
                                       ymin=as.numeric(cello)-.5,
                                       ymax=as.numeric(cello)+.5,
                                       fill=winsor(z, z_win, -z_win))) +
    geom_rect() + theme_bw() +
    facet_grid(sample~chr, scales='free', space='free') +
    ylab('cell') + xlab('position') + 
    scale_fill_gradient2(name="z", low=scales::muted("blue"),
                         high=scales::muted("red")) + 
    theme(axis.text.y=element_blank(), axis.text.x=element_blank(), legend.position='bottom')
  }

  return(ggp.l)
}
