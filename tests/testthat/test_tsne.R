context('tSNE')

## Dummy PCA
nb.cells = 200
cells = paste('barcode', 1:nb.cells)
nb.pcs = 100
pca.o = list(x=matrix(rnorm(nb.cells*nb.pcs), nb.cells))
rownames(pca.o$x) = cells
## Dummy QC data.frame
qc.df = data.frame(cell=cells,
                   tot=rpois(nb.cells, 100),
                   mito=rpois(nb.cells, 1),
                   stringsAsFactors=FALSE)


test_that("tSNE runs with default", {
  tsne.df = run_tsne(pca.o, nb_it=10)
  expect_equal(nrow(tsne.df), nb.cells)
})

test_that("tSNE graphs", {
  tsne.df = run_tsne(pca.o, nb_it=10)
  graphs = plot_tsne(tsne.df, qc.df)
  pdf('temp.pdf')
  lapply(graphs, print)
  dev.off()
  expect_true(file.remove('temp.pdf'))
})
