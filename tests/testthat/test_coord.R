context("Convert to coordinates")

## Generate gene expression
nb.genes = 100
nb.cells = 10
mat = matrix(rpois(nb.cells*nb.genes, 1), nb.genes, nb.cells)
colnames(mat) = paste0('barcode', 1:nb.cells)
ge.df = data.frame(symbol=paste0('gene', 1:nb.genes), stringsAsFactors=FALSE)
ge.df = cbind(ge.df, mat)

## Generate coordinates
genes.coord = data.frame(symbol=ge.df$symbol,
                         chr=sample(c(1:22,"X","Y", 'MT'), nb.genes, TRUE),
                         start=sample.int(10000, nb.genes),
                         stringsAsFactors=FALSE)
genes.coord$end = genes.coord$start + 10

test_that("converts from data.frame and keep everything", {
  coord.df = convert_to_coord(ge.df, genes.coord, chrs=NULL)
  expect_equal(nrow(ge.df), nrow(coord.df))
})

test_that("warning when missing gene", {
  genes.coord = genes.coord[-1,]
  expect_warning(convert_to_coord(ge.df, genes.coord, chrs=NULL), 'missing')
})

test_that("errors when chr names don't match", {
  expect_error(convert_to_coord(ge.df, genes.coord, chrs='XX'), 'match')
})

test_that("removes duplicates", {
  genes.coord[1,2:4] = genes.coord[2,2:4]
  coord.df = convert_to_coord(ge.df, genes.coord, chrs=NULL)
  expect_equal(nrow(coord.df) + 1, nrow(ge.df))
})

test_that("filter chromosomes", {
  coord.df = convert_to_coord(ge.df, genes.coord, chrs=unique(genes.coord$chr)[-1])
  expect_true(nrow(coord.df) < nrow(ge.df))
})

test_that("converts from file and keep everything", {
  write.table(genes.coord, file='temp.tsv', row.names=FALSE, sep='\t')
  coord.df = convert_to_coord(ge.df, 'temp.tsv', chrs=NULL)
  expect_equal(nrow(ge.df), nrow(coord.df))
  file.remove('temp.tsv')
})

