# Do microbial communities separate by habitat? Global Patterns 16S re-analysis

Alpha and beta diversity, PERMANOVA and phylum-level differential abundance on the public `GlobalPatterns` dataset
(Caporaso et al. 2011, [doi:10.1073/pnas.1000080107](https://doi.org/10.1073/pnas.1000080107)), shipped with the
Bioconductor package `phyloseq`.

**Report:** [`report/report.pdf`](report/report.pdf) (source: `report/report.qmd`).

## Result in brief
- 23 non-control samples; habitat (host-associated vs free-living) explains 12% of Bray-Curtis dissimilarity and 22% of weighted-UniFrac distance (PERMANOVA, p = 0.001, the minimum for 999 permutations); sample type explains 69% and 77%.
- Same result without rarefying (12%). Dispersion differs between habitats at p = 0.065, which the report treats as a caveat.
- 21 of 34 phyla differ by habitat (DESeq2, `poscounts`).

## What this adds to the original analysis
Rarefying and DESeq2 are both contested for this kind of data, so the report re-runs each question the other way
(`R/compositional.R`, written to `results/compositional.json`):
- **Rarefying is not load-bearing.** On an Aitchison distance (CLR, no rarefying, every read used) habitat explains
  14% and sample type 78%, stable at 13-14% across OTU prevalence filters. But `betadisper` on that distance gives
  p = 0.0022 against 0.065 on Bray-Curtis, so the unequal-spread caveat gets *stronger*, not weaker.
- **The taxon count is a property of the method.** ANCOM-BC2, which models each sample's sampling fraction, calls
  5 of 34 phyla and 75 of 292 genera where DESeq2-poscounts calls 21 and 208 - while agreeing on effect size
  (r = 0.91 and 0.94). Firmicutes, Tenericutes and Fusobacteria survive both, so the headline stands; "21 of 34"
  does not.
- Limits: 2 to 4 samples per sample type; samples are not independent replicates of "habitat".

## Reproduce
```r
# R >= 4.4; Bioconductor: phyloseq, DESeq2, ANCOMBC, microbiome; CRAN: vegan, ggplot2, dplyr, gridExtra, knitr, jsonlite
Rscript R/compositional.R      # ~1 min; writes results/compositional.json
quarto render report/report.qmd
```

## Provenance
Written and run with the assistance of Claude Code (Anthropic). I am responsible for checking and understanding the results. No part derives from university coursework or assessed material.
`report/apa.csl` is the APA style from the [CSL styles repository](https://github.com/citation-style-language/styles) (CC BY-SA 3.0).
