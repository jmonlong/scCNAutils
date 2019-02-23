context("Rebin coverage")

## Generate scRNA data
nb.genes = 50
scrna.df = data.frame(chr=1,
                      start=sample.int(10000, nb.genes),
                      stringsAsFactors=FALSE)
scrna.df$end = scrna.df$start + 100

## Generate scDNA data
nb.genes = 500
nb.cells = 10
mat = matrix(rnorm(nb.cells*nb.genes),
             nb.genes, nb.cells, byrow=TRUE)
colnames(mat) = paste0('barcode', 1:nb.cells)
cov.df = data.frame(chr=1,
                    start=sample.int(10000, nb.genes),
                    stringsAsFactors=FALSE)
cov.df$end = cov.df$start + 10
cov.df = cbind(cov.df, mat)


test_that("rebin without errors", {
  cov.m = rebin_cov(cov.df, scrna.df)
  expect_equal(nrow(cov.m), nrow(scrna.df))
})
