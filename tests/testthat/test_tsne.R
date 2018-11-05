context('tSNE')

## Dummy PCA
nb.cells = 200
cells = paste0('barcode', 1:nb.cells)
nb.pcs = 100
pca.o = list(x=matrix(rnorm(nb.cells*nb.pcs), nb.cells))
rownames(pca.o$x) = cells

## Dummy QC data.frame
qc.df = data.frame(cell=cells,
                   tot=rpois(nb.cells, 100),
                   mito=rpois(nb.cells, 1),
                   G1.S=rnorm(nb.cells),
                   G2.M=rnorm(nb.cells),
                   stringsAsFactors=FALSE)

## Dummy communities
comm.df = data.frame(cell=cells,
                     community=factor(sample.int(3, nb.cells, TRUE)),
                     stringsAsFactors=FALSE)

## Dummy merge info
info.df = data.frame(cell=cells,
                     sample=paste0('s', sample.int(2, nb.cells, TRUE)),
                     stringsAsFactors=FALSE)

test_that("tSNE runs with default", {
  tsne.df = run_tsne(pca.o, nb_it=10)
  expect_equal(nrow(tsne.df), nb.cells)
})

test_that("tSNE graphs", {
  tsne.df = run_tsne(pca.o, nb_it=10)
  graphs = plot_tsne(tsne.df, qc.df, comm.df, info.df)
  pdf('temp.pdf')
  lapply(graphs, print)
  dev.off()
  expect_true(file.remove('temp.pdf'))
})
