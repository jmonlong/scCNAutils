context("QC per cell")

## Generate gene expression
nb.genes = 100
nb.cells = 10
mat = matrix(rpois(nb.cells*nb.genes, 1), nb.genes, nb.cells)
colnames(mat) = paste('barcode', 1:nb.cells)
df = data.frame(symbol=paste0('gene', 1:nb.genes), stringsAsFactors=FALSE)
df = cbind(df, mat)

test_that("computes simple QC, no mito genes", {
  qc.df = qc_cells(df)
  expect_equal(nrow(qc.df), ncol(df) - 1)
  expect_gt(mean(qc.df$tot), 0)
  expect_gt(mean(qc.df$zeros), 0)
  expect_true(all(qc.df$mito==0))
})

test_that("computes simple QC with mito genes", {
  df$symbol[sample.int(nrow(df), 3)] = paste0('MT-', 1:3)
  qc.df = qc_cells(df)
  expect_true(any(qc.df$mito>0))
})

test_that("computes simple QC with one mito gene", {
  df$symbol[sample.int(nrow(df), 1)] = 'MT-1'
  qc.df = qc_cells(df)
  expect_true(any(qc.df$mito>0))
})

test_that("computes cell cycle score", {
  cell_cycle = data.frame(symbol=sample(df$symbol, 4),
                          phase=rep(c('G1.S', 'G2.M'), 2),
                          stringsAsFactors=FALSE)
  qc.df = qc_cells(df, cell_cycle=cell_cycle)
  expect_true(all(c('G1.S', 'G2.M') %in% colnames(qc.df)))
  expect_true(all(!is.na(qc.df$G1.S)))
  expect_true(all(!is.na(qc.df$G2.M)))
})

test_that("apply mito filters", {
  df$symbol[sample.int(nrow(df), 3)] = paste0('MT-', 1:3)
  qc.df = qc_cells(df)
  mito.prop = sort(unique(qc.df$mito/qc.df$tot))
  df.filt = qc_filter(df, qc.df, max_mito_prop=mito.prop[2])
  expect_gt(ncol(df), ncol(df.filt))
})

test_that("apply total filters", {
  qc.df = qc_cells(df)
  df.filt = qc_filter(df, qc.df, min_total_exp=max(qc.df$tot)-10)
  expect_gt(ncol(df), ncol(df.filt))
})

test_that("apply no filters", {
  qc.df = qc_cells(df)
  df.filt = qc_filter(df, qc.df)
  expect_equal(ncol(df), ncol(df.filt))
})

test_that("qc graphs don't thow an error", {
  qc.df = qc_cells(df)
  ggp = plot_qc_cells(qc.df)
  pdf('temp.pdf')
  lapply(ggp, print)
  dev.off()
  expect_gt(length(ggp), 0)
  file.remove('temp.pdf')
})

