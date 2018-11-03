##' If genes_coord is a filename, the file is expected to be a tab-delimited
##' file with four columns: 'chr', 'start', 'end', 'symbol'. The order of the
##' columns is not important.
##'
##' The gene names in column 'symbol' should match the gene names in the input
##' ge_df. 
##' @title Convert gene symbols to coordinates
##' @param ge_df the data.frame with gene expression and one column 'symbol'
##' with gene names.
##' @param genes_coord either a file name or a data.frame with coordinates
##' and gene names.
##' @param chrs the chromosome names to keep. NULL to include all the
##' chromosomes.
##' @param rm_dup remove duplicated coordinates? Default is TRUE.
##' @return a data.frame with columns 'chr', 'start', 'end' columns with genes
##' coordinates (and still one column per barcode).
##' @author Jean Monlong
convert_to_coord <- function(ge_df, genes_coord, chrs=c(1:22,"X","Y"),
                             rm_dup=TRUE){
  if(!is.data.frame(genes_coord)){
    if(is.character(genes_coord) & length(genes_coord) == 1){
      genes_coord = utils::read.table(genes_coord, as.is=TRUE, header=TRUE)
    } else {
      stop('genes_coord must be either a file name or a data.frame.')
    }
  }
  if(any(!(ge_df$symbol %in% genes_coord$symbol))){
    missing.genes = setdiff(ge_df$symbol, genes_coord$symbol)
    warning('Some genes are missing and will be removed: ', missing.genes)
  }
  ## Merge coordinates and remove gene names
  ge_df = merge(genes_coord, ge_df)
  ge_df$symbol = NULL
  ## Filter chromosomes
  if(!is.null(chrs)){
    if(!any(ge_df$chr %in% chrs)){
      stop("Chromosome names in chrs don't match the input genes coord.")
    }
    ge_df = ge_df[which(ge_df$chr %in% chrs),]
  }
  ## Remove duplicates
  if(rm_dup & any(duplicated(ge_df[, c('chr','start','end')]))){
    ge_df = ge_df[which(!duplicated(ge_df[, c('chr','start','end')])),]
  }
  ## Order by coordinates
  ge_df = ge_df[order(ge_df$chr, ge_df$start),]
  return(ge_df)
}
