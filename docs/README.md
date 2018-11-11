## Reference manual

- Access within R. E.g. `?plot_tsne`.
- PDF version: [scCNAutils-manual.pdf](scCNAutils-manual.pdf).

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

## Metacells

Because each cell has a limited depth, some analysis is performed at the metacell level. 
A metacell is just several cells merged to increase depth and resolution. 
In general we aim at creating metacells from around the same number of cells, to make sure the metacells have comparable depth in later analysis.
There are currently two ways implemented in the `make_metacells` functions to create metacells:

1. Random selection of cells in pre-defined groups of cells (e.g. community).
1. Clustering of cells into similar-size clusters, then merged into metacells.

With the first approach, we can investigate CNAs at the community level. 
In practice, we create a few metacells (e.g. 10) for each community that are then used to call CNA.
CNA consistently present in metacells from a community are most likely shared by the majority of the cells in this community.

```r
mc.o = make_metacells(ge_df, comm_df, nb_metacells=10, metacell_size=5)
```

where 

- *ge_df* contains coordinates columns (chr/start/end) and expression (or coverage) for each cell.
- *comm_df* contains a column *cell* and and column *community*, defining the grouping for each cell.

The second approach is used to investigate all the cells, creating metacells that capture as much as possible rare signal shared by just a few cells.
The clustering is still performed separately in each community but the cells are not randomly selected, they are clustered to form metacells of similar cells within each community.
The communities are used to speed up the clustering and simplify summary at the community level.
Hence the community definition has much less an effect on the final metacells than in the first strategy.
To use this approach, specify `nb_metacells=NULL`.

```r
mc.o = make_metacells(ge_df, comm_df, nb_metacells=NULL, metacell_size=5)
```

The output (*mc.o*) is a list with 

- *ge* the gene expression (or coverage) in each metacell.
- *info* a data.frame with information about each metacell and its corresponding community.
- *mc_cells* a data.frame with information about which cells were used in each metacell.
