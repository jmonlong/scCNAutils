## Data from Tirosh et al. 2016 https://portals.broadinstitute.org/single_cell/study/oligodendroglioma-intra-tumor-heterogeneity
## "log2(TPM/10+1) of all genes in all cells"
## https://portals.broadinstitute.org/single_cell/data/public/oligodendroglioma-intra-tumor-heterogeneity?filename=OG_processed_data_portal.txt
ge = read.table('OG_processed_data_portal.txt', as.is=TRUE, header=TRUE)
## "cell type assignment"
## https://portals.broadinstitute.org/single_cell/data/public/oligodendroglioma-intra-tumor-heterogeneity?filename=cell_type_assignment_portal.txt
info.df = read.table('cell_type_assignment_portal.txt', as.is=TRUE, header=TRUE, skip=1, sep='\t')

## Prepare input gene expression data.frame
colnames(ge)[1] = 'symbol'
any(is.na(ge))
## NAs !?
ge.mat = as.matrix(ge[,-1])
row.with.na = which(apply(ge.mat, 1, function(x)any(is.na(x))))
ge.na.median = apply(ge.mat[row.with.na,], 1, median, na.rm=TRUE)
summary(ge.na.median)
## Most of those genes have a median expression of 0 so I'll change the NAs to 0
ge.mat[which(is.na(ge.mat))] = 0
## Counts
ge.mat = (exp(ge.mat*log(2))-1)*10
ge[,-1] = ge.mat

## Prepare cell labels
## Pretend that assigned cell types are samples (to color graphs)
info.df = info.df[,1:2]
colnames(info.df) = c('cell', 'sample')

## Genes coordinates from gencode
if(!file.exists('gencode.v29.basic.annotation.gtf.gz')){
  download.file('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.basic.annotation.gtf.gz', 'gencode.v29.basic.annotation.gtf.gz')
}
genc.df = read.table('gencode.v29.basic.annotation.gtf.gz', as.is=TRUE, sep='\t')
genc.df = subset(genc.df, V3 == 'gene')
genc.df$symbol = gsub('.* gene_name ([^;]*);.*', '\\1', genc.df$V9)
genc.df$type = gsub('.* gene_type ([^;]*);.*', '\\1', genc.df$V9)
genc.df = subset(genc.df, type == 'protein_coding')
genc.df = genc.df[,c(1,4,5,10)]
colnames(genc.df) = c('chr','start','end', 'symbol')
genc.df = subset(genc.df, !duplicated(symbol))
genc.df$chr = gsub('chr', '', genc.df$chr)
save(genc.df, file='genes.v29.RData')

## Cell cycle
cc.df = read.csv('cell-cycle-genes-tirosh2016.csv', header=TRUE, as.is=TRUE)

library(scCNAutils)

load('genes.v29.RData')
ge = NULL
cells.df = auto_cna_signal(ge, genc.df, 'og', cell_cycle=cc.df, nb_cores=3)
