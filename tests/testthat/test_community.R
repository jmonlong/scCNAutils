context('Community detection')

## Dummy PCA
nb.cells = 200
cells = paste('barcode', 1:nb.cells)
nb.pcs = 100
pca.o = list(x=matrix(rnorm(nb.cells, nb.pcs), nb.cells))
rownames(pca.o$x) = cells

test_that("communities found with default", {
  comm.df = find_communities(pca.o)
  expect_equal(nrow(comm.df), nb.cells)
  expect_gt(nlevels(comm.df$community), 0)
})

test_that("communities graphs", {
  comm.df = find_communities(pca.o)
  ggp = plot_communities(comm.df)
  pdf('temp.pdf')
  lapply(ggp, print)
  dev.off()
  expect_true(file.remove('temp.pdf'))
})

