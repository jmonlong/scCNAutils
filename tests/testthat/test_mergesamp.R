context('Merge multiple samples')

## Gene expression for sample A
nb.genes = 100
nb.cells = 10
mat = matrix(rpois(nb.cells*nb.genes, 1), nb.genes, nb.cells)
colnames(mat) = paste0('barcode', 1:nb.cells)
df = data.frame(symbol=paste0('gene', 1:nb.genes), stringsAsFactors=FALSE)
df = cbind(df, mat)

## Gene expression for sample B
nb.genes = 100
nb.cells2 = 20
mat = matrix(rpois(nb.cells2*nb.genes, 1), nb.genes, nb.cells2)
colnames(mat) = paste0('barcode', 1:nb.cells2)
df2 = data.frame(symbol=paste0('gene', 1:nb.genes), stringsAsFactors=FALSE)
df2 = cbind(df2, mat)

## Make some genes unique to each sample
idx = sample.int(nb.genes, 3)
df$symbol[idx] = paste0('A', df$symbol[idx])
df2$symbol[idx] = paste0('B', df2$symbol[idx])

test_that("merge with different genes", {
  ge_df = merge_samples(list(sampleA=df, sampleB=df2))
  expect_equal(ncol(ge_df$ge), 1 + nb.cells + nb.cells2)
  expect_equal(nrow(ge_df$info), nb.cells + nb.cells2)
  expect_equal(nrow(ge_df$ge), nb.genes + 3)
})
