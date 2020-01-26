# MARVEL: Output

This document describes the output produced by the MARVEL. 

## Pipeline overview
MARVEL is built using [Nextflow](https://www.nextflow.io/).
The overview of MARVEL can be summarized into the following figure.

![Overview of MARVEL](images/marvel_overview.png)

## Motif scanning & test

**Output directory: `results/test/chunks`**

* `*_profiles.npz`
  * QQ-plot of the corresponding test
* `*_real_results.npz`
  * Test statistics, motif selections and corresponding motif coefficients of tested regioins/genes.
* `*_perm_results.npz`
  * Test statistics, motif selection counts of permutated regions/genes.


## Summary of results

**Output directory: `results/summary`**

* `*.qqplot.pdf`
  * QQ-plot of the corresponding test
* `*.table.xlsx`
  * Excel spreadsheet of top association enhancers/promoters/genes.
