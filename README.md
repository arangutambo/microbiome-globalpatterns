# Do microbial communities separate by habitat? Global Patterns 16S re-analysis

Alpha and beta diversity, PERMANOVA and phylum-level differential abundance on the public `GlobalPatterns` dataset
(Caporaso et al. 2011, [doi:10.1073/pnas.1000080107](https://doi.org/10.1073/pnas.1000080107)), shipped with the
Bioconductor package `phyloseq`.

**Report:** [`report/report.pdf`](report/report.pdf) (source: `report/report.qmd`).

## Result in brief
- 23 non-control samples; habitat (host-associated vs free-living) explains 12% of Bray-Curtis dissimilarity and 22% of weighted-UniFrac distance (PERMANOVA, p = 0.001, the minimum for 999 permutations); sample type explains 69% and 77%.
- Same result without rarefying (12%). Dispersion differs between habitats at p = 0.065, which the report treats as a caveat.
- 21 of 34 phyla differ by habitat (DESeq2, `poscounts`).
- Limits: 2 to 4 samples per sample type; samples are not independent replicates of "habitat".

## Reproduce
```r
# R >= 4.4; Bioconductor: phyloseq, DESeq2; CRAN: vegan, ggplot2, dplyr, gridExtra, knitr
quarto render report/report.qmd
```

## Provenance
Written and run with the assistance of Claude Code (Anthropic). I am responsible for checking and understanding the results. No part derives from university coursework or assessed material.
`report/apa.csl` is the APA style from the [CSL styles repository](https://github.com/citation-style-language/styles) (CC BY-SA 3.0).
