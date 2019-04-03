##' Add columns with the names of the genes in the region, total
##' and expressed genes only.
##'
##' The subset of "expressed" genes is made out of genes with
##' non-zero expression in at least 10% of cells or with a mean
##' expression higher than 0.5.
##' @title Annotate CNAs with gene information
##' @param cna_df a data.frame with CNA calls. E.g. data.frames created
##' by \code{call_cna_*} functions.
##' @param gene_info a data.frame with gene information created by the
##' \code{gene_info} function
##' @param cancer_genes a vector with the names of cancer genes. For an
##' additional column. Use if non-NULL.
##' @return the input data.frame with two new columns with all/expressed genes.
##' @author Jean Monlong
##' @export
##' @importFrom magrittr %>%
##' @importFrom rlang .data
annotate_cna <- function(cna_df, gene_info, cancer_genes=NULL){
  ## Define expressed genes (see gene_info function for info columns)
  gene_info$expressed = gene_info$prop.non0 > .1 | gene_info$exp.mean > .5

  ## Create GRanges objects
  cna.gr = GenomicRanges::makeGRangesFromDataFrame(cna_df)
  gene_info = gene_info[order(gene_info$chr, gene_info$start, gene_info$end),]
  gene.gr = GenomicRanges::makeGRangesFromDataFrame(gene_info)

  ## Overlap CNAs and genes
  ol = GenomicRanges::findOverlaps(cna.gr, gene.gr)
  ol = GenomicRanges::as.data.frame(ol)
  ol$gene = gene_info$symbol[ol$subjectHits]
  ol$expressed = gene_info$expressed[ol$subjectHits]

  ## Compute value for new columns
  ol.exp = ol %>% dplyr::filter(.data$expressed) %>%
    dplyr::group_by(.data$queryHits) %>%
    dplyr::summarize(genes=paste(.data$gene, collapse=';'))
  ol.all = ol %>% 
    dplyr::group_by(.data$queryHits) %>%
    dplyr::summarize(genes=paste(.data$gene, collapse=';'))

  ## Add columns to CNA data.frame
  cna_df$exp.genes = NA
  cna_df$exp.genes[ol.exp$queryHits] = ol.exp$genes
  cna_df$all.genes = NA
  cna_df$all.genes[ol.all$queryHits] = ol.all$genes

  ## Cancer genes if specified
  if(!is.null(cancer_genes)){
    ol$cancer = ol$gene %in% cancer_genes
    ol.cancer = ol %>% dplyr::filter(.data$cancer) %>%
      dplyr::group_by(.data$queryHits) %>%
      dplyr::summarize(genes=paste(.data$gene, collapse=';'))
    cna_df$cancer.genes = NA
    cna_df$cancer.genes[ol.cancer$queryHits] = ol.cancer$genes
  }
  
  return(cna_df)
}
