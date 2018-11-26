context('Community detection')

## Dummy PCA
nb.cells = 200
cells = paste0('barcode', 1:nb.cells)
nb.pcs = 100
pca.o = list(x=matrix(rnorm(nb.cells, nb.pcs), nb.cells))
rownames(pca.o$x) = cells

## Dummy QC data.frame
qc.df = data.frame(cell=cells,
                   tot=rpois(nb.cells, 100),
                   mito=rpois(nb.cells, 1),
                   G1.S=rnorm(nb.cells),
                   G2.M=rnorm(nb.cells),
                   stringsAsFactors=FALSE)

## Dummy merge info
info.df = data.frame(cell=cells,
                     sample=paste0('s', sample.int(2, nb.cells, TRUE)),
                     stringsAsFactors=FALSE)

test_that("communities found with default", {
  comm.l = find_communities(pca.o)
  comm.df = comm.l$comm
  expect_equal(nrow(comm.df), nb.cells)
  expect_gt(nlevels(comm.df$community), 0)
})

test_that("communities graphs", {
  comm.l = find_communities(pca.o)
  comm.df = comm.l$comm
  ggp = plot_communities(comm.df, qc.df, info.df)
  pdf('temp.pdf')
  lapply(ggp, print)
  dev.off()
  expect_true(file.remove('temp.pdf'))
})

