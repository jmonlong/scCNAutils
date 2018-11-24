context('Aneuploidy functions')

## Generate gene expression
nb.genes = 100
nb.cells = 50
mat = matrix(rpois(nb.cells*nb.genes, 3),
             nb.genes, nb.cells, byrow=TRUE)
colnames(mat) = paste0('barcode', 1:nb.cells)
df = data.frame(chr=sample(1:3, nb.genes, TRUE),
                start=sample.int(10000, nb.genes),
                stringsAsFactors=FALSE)
df$end = df$start + 10
df = cbind(df, mat)

## Community information
comm_df = data.frame(cell=colnames(mat),
                     community=sample(1:2,nb.cells, TRUE),
                     stringsAsFactors=FALSE)

test_that("graphs don't raise errors", {
  ggp = plot_aneuploidy(df, comm_df, max_cells=10)
  pdf('temp.pdf')
  tmp = lapply(ggp, print)
  dev.off()
  expect_true(file.remove('temp.pdf'))  
})
