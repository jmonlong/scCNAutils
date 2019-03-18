##' @title Heatmap of the CNA scores
##' @param z_df the data.frame with Z-scores.
##' @param cells_df a data.frame with cell information.
##' @param nb_subsamp the number of cells to sub-sample. Default is 1000.
##' @param hc.method the hierarchical clustering method to order cells.
##' @return a list of ggplot2 objects.
##' @author Jean Monlong
##' @external
##' @import ggplot2
plot_heatmap <- function(z_df, cells_df=NULL, nb_subsamp=1000, hc.method='ward.D'){
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
  z_df$cell = factor(z_df$cell, levels=hc.o$labels[hc.o$order])
  
  ## Merge cells info
  cells_df = cells_df[, intersect(colnames(cells_df), c('cell', 'community', 'sample'))]
  z_df = merge(z_df, cells_df)

  ## Order for chromosomes
  if(all(grepl('chr', z_df$chr))){
    z_df$chr = factor(z_df$chr, levels=paste0('chr', c(1:22, 'X','Y')))
  } else {
    z_df$chr = factor(z_df$chr, levels=c(1:22, 'X','Y'))
  }

  ## ggplot2
  ggp.l = list()
  ggp.l$all = ggplot(z_df, aes(x=start, xend=end, y=as.numeric(cell),
                               yend=as.numeric(cell)+1)) +
    geom_rect() + theme_bw() +
    facet_grid(.~chr, scales='free') +
    ylab('cell') + xlab('position') + 
    theme(axis.text.y=element_blank(), axis.text.x=element_blank())

  ## Group by communities
  if('community' %in% colnames(z_df)){
    cell.order = unique(z_df$cell[order(z_df$community, z_df$cell)])
    z_df$cello = factor(z_df$cell, levels=cell.order)
    ggp.l$community = ggplot(z_df, aes(x=start, xend=end, y=as.numeric(cello),
                                       yend=as.numeric(cello)+1)) +
    geom_rect() + theme_bw() +
    facet_grid(community~chr, scales='free') +
    ylab('cell') + xlab('position') + 
    theme(axis.text.y=element_blank(), axis.text.x=element_blank())
  }

  ## Group by sample
  if('sample' %in% colnames(z_df)){
    cell.order = unique(z_df$cell[order(z_df$sample, z_df$cell)])
    z_df$cello = factor(z_df$cell, levels=cell.order)
    ggp.l$community = ggplot(z_df, aes(x=start, xend=end, y=as.numeric(cello),
                                       yend=as.numeric(cello)+1)) +
    geom_rect() + theme_bw() +
    facet_grid(sample~chr, scales='free') +
    ylab('cell') + xlab('position') + 
    theme(axis.text.y=element_blank(), axis.text.x=element_blank())
  }

  return(ggp.l)
}
