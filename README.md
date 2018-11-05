# scCNAutils
Functions to analyze copy number aberrations in single-cell data


## Workflow

### From gene expression to CNA signal

The `auto_cna_signal` function calls the appropriate functions to go from raw expression to communities and tSNE based on CNA signal.
The internal workflow is as follow:

![](docs/flowchart-cnasignal.png)

### From communities to CNA


## Todo

- Example on public data.
- Integrate QC info in `plot_communities`.
- Integrate communities in `plot_tsne`.
- Add multi-samples capability and graphs.
- Add meta-cell based CNA calling.
