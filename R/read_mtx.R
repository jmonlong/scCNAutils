##' @title Read a trio of genes, barcodes and mtx files.
##' @param mtx_file the path to the mtx file
##' @param genes_file the path to the genes file. 
##' @param barcodes_file the path to the barcodes file
##' @param path the path to the folder containing the files
##' @param rm_dup remove duplicated gene names? Default is TRUE.
##' @param genes_col the column to use in genes_file. Default is 2.
##' @return a data.frame with a 'symbol' column with gene names and one
##' column per barcode.
##' @author Jean Monlong
##' @export
read_mtx <- function(mtx_file='matrix.mtx', genes_file='genes.tsv',
                     barcodes_file='barcodes.tsv', path='.', rm_dup=TRUE,
                     genes_col=2){
  mat = Matrix::readMM(paste0(path, '/', mtx_file))
  mat = as.data.frame(as.matrix(mat))
  genes.df = utils::read.table(paste0(path, '/', genes_file), as.is=TRUE)
  mat = cbind(genes.df[,genes_col], mat)
  colnames(mat) = c('symbol', scan(paste0(path, '/', barcodes_file),'',
                                   quiet=TRUE))
  ## Remove duplicated gene names
  if(rm_dup & any(duplicated(mat$symbol))){
    mat = mat[which(!duplicated(mat$symbol)),]
  }
  return(mat)
}
