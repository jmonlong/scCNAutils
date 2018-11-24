context('Metacells')

## Gene expression
nb.genes = 100
nb.cells = 100
mat = matrix(rpois(nb.cells*nb.genes, sample.int(3, nb.cells, TRUE)),
             nb.genes, nb.cells, byrow=TRUE)
cells = paste0('barcode', 1:nb.cells)
colnames(mat) = cells
tot.raw = colSums(mat)
ge.df = data.frame(chr=sample(c(1:22,"X","Y", 'MT'), nb.genes, TRUE),
                   start=sample.int(10000, nb.genes),
                   stringsAsFactors=FALSE)
ge.df$end = ge.df$start + 10
ge.df = cbind(ge.df, mat)

## Dummy communities
comm.df = data.frame(cell=cells,
                     community=factor(sample.int(3, nb.cells, TRUE)),
                     stringsAsFactors=FALSE)

test_that("makes metacells", {
  mc.l = make_metacells(ge.df, comm.df, nb_metacells=2)
  expect_gt(nrow(mc.l$ge), 0)
})

test_that("makes metacells by clustering", {
  mc.l = make_metacells(ge.df, comm.df, nb_metacells=NULL)
  expect_gt(nrow(mc.l$ge), 0)
})

test_that("makes metacells without baseline cells", {
  mc.l = make_metacells(ge.df, comm.df, baseline_cells=FALSE, nb_metacells=2)
  expect_true(!any(grepl('baseline_', mc.l$info$cell)))
})
