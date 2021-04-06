The code from this package was written and used primarily for the following study. 

> Couturier, C.P., Ayyadhury, S., Le, P.U. et al. *Single-cell RNA-seq reveals that glioblastoma recapitulates a normal neurodevelopmental hierarchy.* Nat Commun 11, 3406 (2020). https://doi.org/10.1038/s41467-020-17186-5

The list of cancer/non-cancer cells are available in the [`cells-community-normal-natcom2020.tsv.gz`](cells-community-normal-natcom2020.tsv.gz).
This tsv has three columns: the sample ID, the barcode of the cell, and the community name (in the form of `<SAMPLE>.<CLONE>`, or `normal` for non-cancer cells).

```
sample	barcode	community.patient
BT402_W	ATCATCTTCTTCCTTC-1	normal
BT363_W	CATCGGGAGTCATCCA-1	BT363.1
BT326_GSC	ACGCCGAAGCGATAGC-1	BT326.2
...
```
