context("Gene expression normalization")

## Generate gene expression
nb.genes = 100
nb.cells = 10
mat = matrix(rpois(nb.cells*nb.genes, sample.int(3, nb.cells, TRUE)),
             nb.genes, nb.cells, byrow=TRUE)
colnames(mat) = paste0('barcode', 1:nb.cells)
tot.raw = colSums(mat)
df = data.frame(symbol=paste0('gene', 1:nb.genes), stringsAsFactors=FALSE)
df = cbind(df, mat)

test_that("normalize using the total method", {
  norm.df = norm_ge(df, method='total')
  tot.norm = colSums(norm.df[, colnames(mat)])
  expect_gt(sd(tot.raw), sd(tot.norm))
  expect_true(all(abs(tot.norm-tot.norm[1])<.0000001))
})

test_that("normalize using the tmm method", {
  norm.df = norm_ge(df, method='tmm')
  tot.norm = colSums(norm.df[, colnames(mat)])
  expect_gt(sd(tot.raw), sd(tot.norm))
})

test_that("normalize using the tmm method in parallel", {
  norm.df = norm_ge(df, nb_cores=2)
  tot.norm = colSums(norm.df[, colnames(mat)])
  expect_gt(sd(tot.raw), sd(tot.norm))
})
