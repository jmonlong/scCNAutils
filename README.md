# scCNAutils
Functions to analyze copy number aberrations in single-cell data


## Workflow

### From gene expression to CNA signal

```
read_mtx
  |-----
  |    V
  |  qc_cells
  |    |----------------> graph_qc_cells X
  V    V            |
qc_filter X         |
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
zscore X            |
  |                 |
  V                 |
smooth_movingm      |
  |                 |
  V                 |
run_pca X   <--------
  |-------------------------------
  V                              V
comm_detection X              run_tsne X
  |       |_ _ _ _ _ _ _ _ _ _    |
  V                          V    V
graph_communities X        graph_tsne X
```

### From communities to CNA

