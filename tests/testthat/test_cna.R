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

## Genes
genes.df = data.frame(geneid=1:nb.genes, symbol=paste0('gene', 1:nb.genes),
                      stringsAsFactors=FALSE)
## Add one duplicate
genes.df[1,] = genes.df[2,]
## Coordinates
genes.coord = data.frame(symbol=genes.df$symbol,
                         chr=sample(1:2, nb.genes, TRUE),
                         start=sample.int(10000, nb.genes),
                         stringsAsFactors=FALSE)
genes.coord$end = genes.coord$start + 10

## Gene expression
mat = matrix(rpois(nb.cells*nb.genes, 1), nb.genes, nb.cells)
cells = paste0('barcode', 1:nb.cells)
colnames(mat) = cells
ge = data.frame(symbol=paste0('gene', 1:nb.genes), stringsAsFactors=FALSE)
ge = cbind(ge, mat)




test_that("CNAs are called without errors", {
  cna.o = call_cna(df, mc_info=mc.info)
  expect_gt(nrow(cna.o$seg.df), 0)  
  expect_gt(nrow(cna.o$hmm.df), 0)  
  ggp = plot_cna(cna.o, chrs_order=5:1)
  pdf('temp.pdf')
  tmp = lapply(ggp, print)
  dev.off()
  expect_true(file.remove('temp.pdf'))  
})

test_that("CNAs are called without errors with multisamp", {
  cna.o = call_cna_multisamps(df, mc.info)
  expect_gt(nrow(cna.o$seg.df), 0)  
  expect_gt(nrow(cna.o$hmm.df), 0)  
  ggp = plot_cna(cna.o, chrs_order=5:1)
  pdf('temp.pdf')
  tmp = lapply(ggp, print)
  dev.off()
  expect_true(file.remove('temp.pdf'))  
})

test_that("Get gene info and annotate CNAs", {
  gene.info = gene_info(ge, genes.coord)
  expect_true(any(!is.na(gene.info$exp.mean)))
  expect_true(any(!is.na(gene.info$exp.sd)))
  expect_true(any(!is.na(gene.info$prop.non0)))
  cna.o = call_cna_multisamps(df, mc.info)
  seg.df = annotate_cna(cna.o$seg.df, gene.info)
  expect_gt(nrow(seg.df), 0)
  expect_true(any(!is.na(seg.df$exp.genes)))
  expect_true(any(!is.na(seg.df$all.genes)))
})
