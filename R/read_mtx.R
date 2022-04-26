##' @title Read a trio of genes, barcodes and mtx files.
##' @param mtx_file the path to the mtx file
##' @param genes_file the path to the genes file. 
##' @param barcodes_file the path to the barcodes file
##' @param path the path to the folder containing the files
##' @param rm_dup remove duplicated gene names? Default is TRUE.
##' @param genes_col the column to use in genes_file. Default is 2.
##' @param min_barcode_exp minimum total expression in barcodes/cells. Default is 1.
##' @param filter_zero_genes filter genes with no reads in any cells? Default is TRUE.
##' @return a data.frame with a 'symbol' column with gene names and one
##' column per barcode.
##' @author Jean Monlong
##' @export
read_mtx <- function(mtx_file='matrix.mtx', genes_file='genes.tsv',
                     barcodes_file='barcodes.tsv', path='.', rm_dup=TRUE,
                     genes_col=2, min_barcode_exp=1, filter_zero_genes=TRUE){
  if(!file.exists(paste0(path, '/', mtx_file))){
    if(file.exists(paste0(path, '/', mtx_file, '.gz'))){
      mtx_file = paste0(mtx_file, '.gz')
    } else {
      stop(paste0(path, '/', mtx_file), ' not found')
    }
  }
  mat = Matrix::readMM(paste0(path, '/', mtx_file))

  ## Find genes with non-zero expression
  if(filter_zero_genes){
    chunk.size = ifelse(nrow(mat)<20, nrow(mat), 20)
    nonzero.g = tapply(1:nrow(mat), cut(1:nrow(mat), chunk.size), function(iis){
      mat = as.matrix(mat[iis,])
      iis[which(rowSums(mat)>0)]
    })
    nonzero.g = unlist(nonzero.g)
  } else {
    nonzero.g = 1:nrow(mat)
  }

  ## Find barcodes with enough total expression
  if(min_barcode_exp > 0){
    chunk.size = ifelse(ncol(mat)<20, ncol(mat), 20)
    exp.b = tapply(1:ncol(mat), cut(1:ncol(mat), chunk.size), function(iis){
      mat = as.matrix(mat[,iis])
      colSums(mat)
    })
    exp.b = unlist(exp.b)
    nonzero.b = which(exp.b>=min_barcode_exp)
  } else {
    nonzero.b = 1:ncol(mat)
  }

  ## Subset the matrix and convert to data.frame
  mat = as.data.frame(as.matrix(mat[nonzero.g,nonzero.b]))

  ## Get gene and barcode names
  genes.df = utils::read.table(paste0(path, "/", genes_file),
                               as.is = TRUE)
  genes = genes.df[nonzero.g, 2]
  genes = gsub('hg19_', '', genes)
  mat = cbind(genes, mat)
  barcodes = scan(paste0(path, "/", barcodes_file), "", quiet = TRUE)
  colnames(mat) = c("symbol", barcodes[nonzero.b])

  ## Remove duplicated gene names
  if(rm_dup & any(duplicated(mat$symbol))){
    mat = mat[which(!duplicated(mat$symbol)),]
  }
  return(mat)
}
