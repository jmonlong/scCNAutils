context("Smoothing")

## Generate data
nb.genes = 500
nb.cells = 10
mat = matrix(rnorm(nb.cells*nb.genes),
             nb.genes, nb.cells, byrow=TRUE)
colnames(mat) = paste0('barcode', 1:nb.cells)
df = data.frame(chr=sample(c(1:5), nb.genes, TRUE),
                start=sample.int(10000, nb.genes),
                stringsAsFactors=FALSE)
df$end = df$start + 10
df = cbind(df, mat)

sd.raw = sd(as.matrix(df[, colnames(mat)]))

test_that("smooth on single core", {
  smooth.df = smooth_movingw(df)
  sd.smooth = sd(as.matrix(smooth.df[, colnames(mat)]))
  expect_lt(sd.smooth, sd.raw)  
})

test_that("smooth in parallel", {
  smooth.df = smooth_movingw(df, nb_cores=2)
  sd.smooth = sd(as.matrix(smooth.df[, colnames(mat)]))
  expect_lt(sd.smooth, sd.raw)  
})

test_that("smooth with mean", {
  smooth.df = smooth_movingw(df, FUN=mean)
  sd.smooth = sd(as.matrix(smooth.df[, colnames(mat)]))
  expect_lt(sd.smooth, sd.raw)  
})

