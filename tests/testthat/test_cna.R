context('CNA calling')

## Zscores
nb.genes = 500
nb.cells = 10
mat = matrix(rnorm(nb.cells*nb.genes),
             nb.genes, nb.cells, byrow=TRUE)
cells = paste0('barcode', 1:nb.cells)
colnames(mat) = cells
df = data.frame(chr=sample(c(1:5), nb.genes, TRUE),
                start=sample.int(10000, nb.genes),
                stringsAsFactors=FALSE)
df$end = df$start + 10
df = cbind(df, mat)

## Metacell info
mc.info = data.frame(cell=cells, community=factor(rep(1:2,5)), stringsAsFactors=FALSE)

test_that("CNAs are called without errors", {
  cna.df =call_cna(df, mc_info=mc.info)
  expect_gt(nrow(cna.df), 0)  
})

test_that("graphs don't throw errors", {
  cna.df =call_cna(df, mc_info=mc.info)
  ggp = plot_cna(cna.df, chrs_order=5:1)
  pdf('temp.pdf')
  print(ggp)
  dev.off()
  expect_true(file.remove('temp.pdf'))  
})
