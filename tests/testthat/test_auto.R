context('Automated pipelines')

## Generate files
nb.genes = 100
nb.cells = 200

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

## Coordinates
genes.coord = data.frame(symbol=genes.df$symbol,
                         chr=sample(c(1:22,"X","Y", 'MT'), nb.genes, TRUE),
                         start=sample.int(10000, nb.genes),
                         stringsAsFactors=FALSE)
genes.coord$end = genes.coord$start + 10

## Cell cycle genes
cell_cycle = data.frame(symbol=sample(genes.df$symbol, 4),
                        phase=rep(c('G1.S', 'G2.M'), 2),
                        stringsAsFactors=FALSE)


test_that("CNA signal runs", {
  res.df = auto_cna_signal('.', genes.coord, prefix='tempfortest', cell_cycle=cell_cycle)
  outfiles = list.files('.', 'tempfortest')
  expect_gt(length(outfiles), 0)
  expect_true(all(file.remove(outfiles)))
})

## Remove files
file.remove('matrix.mtx', 'genes.tsv', 'barcodes.tsv')
