context("Bin genes")

## Generate gene expression
nb.genes = 100
nb.cells = 10
mat = matrix(rpois(nb.cells*nb.genes, 1),
             nb.genes, nb.cells, byrow=TRUE)
colnames(mat) = paste0('barcode', 1:nb.cells)
df = data.frame(chr=sample(c(1:22,"X","Y", 'MT'), nb.genes, TRUE),
                start=sample.int(10000, nb.genes),
                stringsAsFactors=FALSE)
df$end = df$start + 10
df = cbind(df, mat)

test_that("bin genes", {
  bin.df = bin_genes(df)
  expect_lt(nrow(bin.df), nrow(df))
})

test_that("bin genes in parallel", {
  bin.df = bin_genes(df, nb_cores=2)
  expect_lt(nrow(bin.df), nrow(df))
})
