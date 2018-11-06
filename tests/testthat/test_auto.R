context('Automated pipelines')

## Generate files
nb.genes = 100
nb.cells = 200

## Sample A
dir.create('sampleA')
## Expression
mat = matrix(rpois(nb.cells*nb.genes, 1), nb.genes, nb.cells)
mat = Matrix::Matrix(mat, sparse=TRUE)
Matrix::writeMM(mat, 'sampleA/matrix.mtx')

## Genes
genes.df = data.frame(geneid=1:nb.genes, symbol=paste0('gene', 1:nb.genes),
                      stringsAsFactors=FALSE)
## Add one duplicate
genes.df[1,] = genes.df[2,]
write.table(genes.df, file='sampleA/genes.tsv', col.names=FALSE,
            row.names=FALSE, sep='\t')

## Barcodes
barcodes = paste0('barcode', 1:nb.cells)
write.table(barcodes, file='sampleA/barcodes.tsv', col.names=FALSE,
            row.names=FALSE, sep='\t')

## Sample B
dir.create('sampleB')
## Expression
mat = matrix(rpois(nb.cells*nb.genes, 1), nb.genes, nb.cells)
mat = Matrix::Matrix(mat, sparse=TRUE)
Matrix::writeMM(mat, 'sampleB/matrix.mtx')

## Genes
genes.df = data.frame(geneid=1:nb.genes, symbol=paste0('gene', 1:nb.genes),
                      stringsAsFactors=FALSE)
## Add one duplicate
genes.df[1,] = genes.df[2,]
write.table(genes.df, file='sampleB/genes.tsv', col.names=FALSE,
            row.names=FALSE, sep='\t')

## Barcodes
barcodes = paste0('barcode', 1:nb.cells)
write.table(barcodes, file='sampleB/barcodes.tsv', col.names=FALSE,
            row.names=FALSE, sep='\t')

## Coordinates
genes.coord = data.frame(symbol=genes.df$symbol,
                         chr=sample(1:2, nb.genes, TRUE),
                         start=sample.int(10000, nb.genes),
                         stringsAsFactors=FALSE)
genes.coord$end = genes.coord$start + 10

## Cell cycle genes
cell_cycle = data.frame(symbol=sample(genes.df$symbol, 4),
                        phase=rep(c('G1.S', 'G2.M'), 2),
                        stringsAsFactors=FALSE)


test_that("CNA signal runs for one sample", {
  res.df = auto_cna_signal('sampleA', genes.coord, prefix='tempfortest', cell_cycle=cell_cycle)
  load('tempfortest-ge-coord-norm.RData')
  res.df$community = factor(sample.int(2, nrow(res.df), TRUE))
  cna.df = auto_cna_call(data, res.df, prefix='tempfortest',
                         baseline_communities=1)
  outfiles = list.files('.', 'tempfortest')
  expect_gt(length(outfiles), 0)
  expect_true(all(file.remove(outfiles)))
})

test_that("CNA signal runs for two samples", {
  expect_warning({res.df = auto_cna_signal(c('sampleA','sampleB'), genes.coord, prefix='tempfortest', cell_cycle=cell_cycle)}, 'path')
  outfiles = list.files('.', 'tempfortest')
  expect_gt(length(outfiles), 0)
  expect_true(all(file.remove(outfiles)))
})

## Remove files
unlink(c('sampleA', 'sampleB'), recursive = TRUE)
