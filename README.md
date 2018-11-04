# scCNAutils
Functions to analyze copy number aberrations in single-cell data


## Workflow

### From gene expression to CNA signal

```
read_mtx
  |-----
  |    V
  |  qc_cells
  |    |----------------> plot_qc_cells
  V    V            |
qc_filter           |
  |                 |
  V                 |
convert_to_coord    |
  |                 |
  V                 | 
norm_ge             |
  |                 |
  V                 |
bin_genes           |
  |                 |
  V                 |
zscore              |
  |                 |
  V                 |
smooth_movingm      |
  |                 |
  V                 |
run_pca     <--------
  |-------------------------------
  V                              V
comm_detection X              run_tsne X
  |       |_ _ _ _ _ _ _ _ _ _    |
  V                          V    V
plot_communities X          plot_tsne X
```

### From communities to CNA

