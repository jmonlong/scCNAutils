context("Reading files")

## Generate files
nb.genes = 100
nb.cells = 10

## Expression
mat = matrix(rpois(nb.cells*nb.genes, 1), nb.genes, nb.cells)
mat = Matrix::Matrix(mat, sparse=TRUE)
Matrix::writeMM(mat, 'matrix.mtx')

## Genes
genes.df = data.frame(geneid=1:nb.genes, symbol=paste0('gene', 1:nb.genes),
                      stringsAsFactors=FALSE)
## Add one duplicate
genes.df[1,] = genes.df[2,]
write.table(genes.df, file='genes.tsv', col.names=FALSE,
            row.names=FALSE, sep='\t')

## Barcodes
barcodes = paste('barcode', 1:nb.cells)
write.table(barcodes, file='barcodes.tsv', col.names=FALSE,
            row.names=FALSE, sep='\t')

test_that("reads correctly", {
  df = read_mtx(rm_dup=FALSE)
  expect_true(all(df$symbol == genes.df$symbol))
  expect_true(all(colnames(df) == c('symbol', barcodes)))
  df.exp = as.matrix(df[,-1])
  expect_true(all(df.exp == as.matrix(mat)))
})

test_that("removes duplicate", {
  df = read_mtx()
  expect_equal(nrow(df), nrow(genes.df) - 1)
})

## Remove files
file.remove('matrix.mtx', 'genes.tsv', 'barcodes.tsv')
