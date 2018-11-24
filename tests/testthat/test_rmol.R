context("Remove outliers")

## Generate gene expression
nb.genes = 100
nb.cells = 10
mat = matrix(rpois(nb.cells*nb.genes, 3),
             nb.genes, nb.cells, byrow=TRUE)
## Add some outliers
mat[1:10, 1:5] = matrix(rpois(50, 100), 10, 5)
colnames(mat) = paste0('barcode', 1:nb.cells)
df = data.frame(symbol=paste0('gene', 1:nb.genes), stringsAsFactors=FALSE)
df = cbind(df, mat)

test_that("some genes are removed", {
  rm.df = rm_cv_outliers(df, ol_quant_th=.6, ol_sd_th=1)
  expect_lt(nrow(rm.df), nrow(df))
})

