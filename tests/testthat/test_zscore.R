context("Compute z-score")

## Generate gene expression
nb.genes = 100
nb.cells = 10
mat = matrix(rpois(nb.cells*nb.genes, 1),
             nb.genes, nb.cells, byrow=TRUE)
colnames(mat) = paste('barcode', 1:nb.cells)
df = data.frame(chr=sample(c(1:22,"X","Y", 'MT'), nb.genes, TRUE),
                start=sample.int(10000, nb.genes),
                stringsAsFactors=FALSE)
df$end = df$start + 10
df = cbind(df, mat)

test_that("zscore norm", {
  z.df = zscore(df, method='norm')
  expect_gt(nrow(z.df), 0)
})

test_that("zscore z", {
  z.df = zscore(df, method='z')
  expect_gt(nrow(z.df), 0)
})


