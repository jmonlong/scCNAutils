context('PCA and tSNE')

## Generate data
nb.genes = 100
nb.cells = 100
mat = matrix(rnorm(nb.cells*nb.genes),
             nb.genes, nb.cells, byrow=TRUE)
colnames(mat) = paste('barcode', 1:nb.cells)
df = data.frame(chr=sample(c(1:5), nb.genes, TRUE),
                start=sample.int(10000, nb.genes),
                stringsAsFactors=FALSE)
df$end = df$start + 10
df = cbind(df, mat)
## Gene expression for qc
mat = matrix(rpois(nb.cells*nb.genes, 1), nb.genes, nb.cells)
colnames(mat) = paste('barcode', 1:nb.cells)
ge.df = data.frame(symbol=paste0('gene', 1:nb.genes), stringsAsFactors=FALSE)
ge.df = cbind(ge.df, mat)


test_that("PCA runs with default", {
  pca.o = run_pca(df)
  expect_equal(length(pca.o), 3)
  expect_gt(ncol(pca.o$x), 0)
  expect_equal(length(pca.o$sdev), min(nb.cells, nb.genes))
  expect_equal(nrow(pca.o$x), nb.cells)
})

test_that("PCA runs with core cells", {
  nb.core.cells = nb.cells / 2
  pca.o = run_pca(df, core_cells=sample(colnames(mat), nb.core.cells))
  expect_equal(length(pca.o), 3)
  expect_gt(ncol(pca.o$x), 0)
  expect_equal(length(pca.o$sdev), min(nb.core.cells, nb.genes))
  expect_equal(nrow(pca.o$x), nb.cells)
})

test_that("PCA graph", {
  pca.o = run_pca(df)
  pdf('temp.pdf')
  print(pca.o$sdev.graph)
  dev.off()
  expect_true(file.remove('temp.pdf'))
})

test_that("tSNE runs with default", {
  pca.o = run_pca(df)
  tsne.df = run_tsne(pca.o, nb_it=10)
  expect_equal(nrow(tsne.df), nb.cells)
})

test_that("tSNE graphs", {
  qc.df = qc_cells(ge.df)
  pca.o = run_pca(df)
  tsne.df = run_tsne(pca.o, nb_it=10)
  graphs = plot_tsne(tsne.df, qc.df)
  pdf('temp.pdf')
  lapply(graphs, print)
  dev.off()
  expect_true(file.remove('temp.pdf'))
})

