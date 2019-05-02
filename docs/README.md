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

### Download problems with `biocLite`/`install_github`

I sometimes encounter problems when `biocLite` or *devtools* tries to download the package on the HPC.
Usually it's a problem with Curl or OpenSSL and it's difficult to fix on a HPC where I don't have permissions.
In that case it's easier to **clone the repo locally and install with devtools** (`devtools::install()` from within the cloned repo). 

Of note, the Bioconductor dependencies must be installed separately.
All the dependencies are listed in the [DESCRIPTION](../DESCRIPTION) file or the [manual](scCNAutils-manual.pdf).

So in bash:

```sh
cd FOLDER/WHERE/TO/CLONE/REPO
git clone https://github.com/jmonlong/scCNAutils.git
cd scCNAutils
```

Then in R:

```r
.libPaths('~/R/library/')
source('http://bioconductor.org/biocLite.R')
biocLite('GenomicRanges')
devtools::install()
```

### Compilation problems

Some dependencies, like *RcppHMM*,  require C++ compilation and might raise installation errors if the HPC is not well configured (or uses modules).
I usually try different *gcc* modules until the installation/compilation works.

On Abacus, the modules that worked for me were `mugqic/R_Bioconductor/3.5.0_3.7` and `mugqic_dev/gcc/4.7.2`.

## Optimized version

Some of the slow R functions were converted into Cpp functions. 
It requires some compilation when installing the package and is not as multi-threaded but it is much more memory efficient and faster.
When the data is very big (e.g. >20K cells), we recommend using `rcpp=TRUE` in relevant functions (e.g. `auto_cna_signal` functions).

Note: The Rcpp version is likely faster than with multi-threaded R functions even in small dataset so might be useful in general.

## Getting genes coordinates

[Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build) from 10X Genomics seems to use Ensembl annotations. 
For example:

- `ftp://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz` for hg19.
- `ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz` for GRCh38

Depending on the Cellranger version it might be another version from with similar path.

For example to retrieve the coordinates for each gene name (removing duplicated names) for hg19:

```r
download.file('ftp://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz',
              'Homo_sapiens.GRCh37.87.gtf.gz')
gene.df = read.table('Homo_sapiens.GRCh37.87.gtf.gz', as.is=TRUE, sep='\t')
## Keep genes
gene.df = subset(gene.df, V3 == 'gene')
gene.df$symbol = gsub('.* gene_name ([^;]*);.*', '\\1', gene.df$V9)
## Format and remove duplicates
gene.df = gene.df[,c(1,4,5,10)]
colnames(gene.df) = c('chr','start','end', 'symbol')
gene.df = subset(gene.df, !duplicated(symbol))
## Save
write.table(gene.df, file='genes.coord.tsv', quote=FALSE, row.names=FALSE, sep='\t')
## Check that gene names match with the scRNA-seq data     
ge = read_mtx(path='path/to/some/data')
mean(ge$symbol %in% gene.df$symbol) ## Poportion of genes in annotation
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

## Tuning the community detection

The Louvain algorithm can take a *resolution* parameter gamma as input. 
By default gamma is 1.
Larger gammas result in more communities and vice versa.
One strategy is then to pick the gamma that leads to the most stable communities.
Briefly Louvain is ran several time for each gamma and each resulting is compared to the others with the Rand index. 
The higher the Rand index and the less variable, the more similar/stable are the communities for a particular gamma.

The command to explore gammas, comparing 20 runs for gammas between 0.1 and 1.5:

```r
## pca.o from run_pca
comm.exp = find_communities(pca.o, gamma=seq(.1,1.5,.05), nreps=20, nb_cores=4)
comm.exp$ari.df ## Adjusted Rand Index for each gamma
head(comm.exp$comm) ## communities for best gamma (highest average ARI).
```

## Rebinning coverage data

To compare scRNA-seq data (e.g. CNA Z-scores) with coverage data (e.g. WGS or scDNA-seq) we might want to work on the same regions/bins.
The coverage data will most likely have a better resolution so we could merge these bins into larger bins matching a scRNA-seq analysis.
The function *rebin_cov* does that by overlapping both sets of bins and computing a coverage on the new/larger bins as the weighted sum of the coverage in overlapping bins, with the proportion of overlap as weight.

```r
load('scrnaseq-coord-norm-bin5-z2-smooth7.RData') ## CNA Z-scores from a scRNA-seq analysis
scrna.bins = data[, c('chr','start','end')] ## Optional: keep only coords columns

## cov.df: data.frame with columns chr, start, end and then sample/cell names
cov.new = rebin_cov(cov.df, scrna.bins)
```
