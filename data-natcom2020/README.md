The code from this package was written and used primarily for the following study. 

> Couturier, C.P., Ayyadhury, S., Le, P.U. et al. *Single-cell RNA-seq reveals that glioblastoma recapitulates a normal neurodevelopmental hierarchy.* Nat Commun 11, 3406 (2020). https://doi.org/10.1038/s41467-020-17186-5

The list of cancer/non-cancer cells are available in the [`cells-community-normal-natcom2020.tsv.gz`](cells-community-normal-natcom2020.tsv.gz).
This tsv has three columns: the sample ID, the barcode of the cell, and the community name (in the form of `<SAMPLE>.<CLONE>`, or `normal` for non-cancer cells).

```
sample	barcode	community.patient
OPK322	AAACCTGAGAAGGTTT-1	OPK322.2
OPK363	CCCAGTTAGCTGTTCA-1	OPK363.1
OPK402B	TCGCGTTTCCTCAACC-1	normal
OPK322	GCTGCAGCAGCTGCTG-1	OPK322.1
...
```

Of note, the sample naming in the article's figures uses `BT` instead of `OPK`. 
Also when both **B**ulk (**W**hole sample in the figures) and **S**tem-cells (**GSC**: glioma stem cell in the figures) were sequenced for a sample, a `B`/`S` suffix is added to the sample name here. 
For example, `OPK402B` corresponds to `BT402-W`
