- [Reference manual](scCNAutils-manual.pdf).

## Install on a HPC

### Install the package locally

To install locally I found that the easiest if to create a directory in your home, e.g. `~/R/library`, and add this path to your library path using `.libPaths`.

To install locally:

```r
.libPaths('~/R/library/')
source('http://bioconductor.org/biocLite.R')
biocLite('jmonlong/scCNAutils')
```

From now on, to use the package:

```r
.libPaths('~/R/library/')
library(scCNAutils)
```

### Compilation problems

Some dependencies, like *RcppHMM*,  require C++ compilation and might raise installation errors if the HPC is not well configured (or uses modules).
I usually try different *gcc* modules until the installation/compilation works.

On Abacus, the modules that worked for me were `mugqic/R_Bioconductor/3.5.0_3.7` and `mugqic_dev/gcc/4.7.2`.


## Getting genes coordinates

One option if to use the [Gencode annotation](https://www.gencodegenes.org/human/).
For example to retrieve the coordinates for each gene name (removing duplicated names):

```r
download.file('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.basic.annotation.gtf.gz',
              'gencode.v29.basic.annotation.gtf.gz')
genc.df = read.table('gencode.v29.basic.annotation.gtf.gz', as.is=TRUE, sep='\t')
## Keep protein-coding genes
genc.df = subset(genc.df, V3 == 'gene')
genc.df$symbol = gsub('.* gene_name ([^;]*);.*', '\\1', genc.df$V9)
genc.df$type = gsub('.* gene_type ([^;]*);.*', '\\1', genc.df$V9)
genc.df = subset(genc.df, type == 'protein_coding')
## Format and remove duplicates
genc.df = genc.df[,c(1,4,5,10)]
colnames(genc.df) = c('chr','start','end', 'symbol')
genc.df$chr = gsub('chr', '', genc.df$chr)
genc.df = subset(genc.df, !duplicated(symbol))
## Save
write.table(genc.df, file='genes.v29.tsv', quote=FALSE, row.names=FALSE, sep='\t')
```

## Genes associated with cell cycle

We use genes from [Tirosh et al 2016](https://www.nature.com/articles/nature20123). 
`cell-cycle-genes-tirosh2016.tsv` contains the genes for phase G1/S and G2/M.

