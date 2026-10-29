# Compositional re-analysis, written to results/compositional.json for the report.
#
# The main report rarefies for diversity and uses DESeq2 with poscounts size
# factors for differential abundance. Both choices are contested: McMurdie &
# Holmes (2014) argue rarefying discards data, and compositional methods
# (Gloor et al. 2017; Lin & Peddada 2020) argue that count tables carry only
# relative information, so a method that does not model the sampling fraction
# will call too many taxa. This script re-runs both questions the other way:
#   1. Aitchison distance (CLR + Euclidean) with no rarefying at all.
#   2. ANCOM-BC2 against DESeq2-poscounts, at phylum and at genus level.
# Run:  Rscript R/compositional.R   (a few minutes; ANCOM-BC2 is the slow part)

suppressPackageStartupMessages({library(phyloseq); library(vegan); library(DESeq2); library(ANCOMBC); library(jsonlite)})
data(GlobalPatterns)
gp <- subset_samples(GlobalPatterns, SampleType != "Mock"); gp <- prune_taxa(taxa_sums(gp) > 0, gp)
host <- c("Feces","Skin","Tongue")
sample_data(gp)$habitat <- factor(ifelse(sample_data(gp)$SampleType %in% host, "host-associated","free-living"),
                                  levels = c("free-living","host-associated"))
md <- data.frame(sample_data(gp)); out <- list()

otu <- as(otu_table(gp), "matrix"); if (!taxa_are_rows(gp)) otu <- t(otu)
clr_dist <- function(minp) {
  keep <- rowSums(otu > 0) >= minp
  o <- otu[keep, , drop = FALSE] + 1
  list(d = dist(t(apply(o, 2, function(x) log(x) - mean(log(x))))),
       n_otu = sum(keep), read_share = sum(otu[keep, ]) / sum(otu))
}
out$aitchison <- lapply(c("2","3","5"), function(k) {
  cd <- clr_dist(as.integer(k))
  set.seed(1); h <- adonis2(cd$d ~ habitat, data = md, permutations = 999)
  set.seed(1); t <- adonis2(cd$d ~ SampleType, data = md, permutations = 999)
  list(min_prevalence = as.integer(k), n_otu = cd$n_otu, read_share = round(cd$read_share, 4),
       habitat_r2 = round(h$R2[1], 3), habitat_p = h$`Pr(>F)`[1],
       sampletype_r2 = round(t$R2[1], 3), sampletype_p = t$`Pr(>F)`[1])
})
names(out$aitchison) <- paste0("prev", c(2,3,5))
cd <- clr_dist(3)
out$aitchison_betadisper_p <- round(anova(betadisper(cd$d, md$habitat))$`Pr(>F)`[1], 4)
out$n_otu_total <- nrow(otu)

da_at <- function(rank, min_reads) {
  g <- tax_glom(gp, rank, NArm = TRUE); g <- prune_taxa(taxa_sums(g) >= min_reads, g)
  dds <- suppressMessages(phyloseq_to_deseq2(g, ~ habitat))
  dds <- DESeq(dds, sfType = "poscounts", quiet = TRUE)
  rd <- as.data.frame(results(dds, alpha = 0.05, contrast = c("habitat","host-associated","free-living")))
  rd$taxon <- rownames(rd); rd$name <- as.character(tax_table(g)[rownames(rd), rank])
  rd$sig_d <- !is.na(rd$padj) & rd$padj < 0.05
  set.seed(1)
  ab <- suppressWarnings(ancombc2(data = g, fix_formula = "habitat", p_adj_method = "BH",
          prv_cut = 0, lib_cut = 0, alpha = 0.05, verbose = FALSE, n_cl = 1))
  r <- ab$res
  a <- data.frame(taxon = r$taxon, lfc = r[[grep("^lfc_habitat", names(r), value=TRUE)]],
                  q = r[[grep("^q_habitat", names(r), value=TRUE)]],
                  sig_a = r[[grep("^diff_habitat", names(r), value=TRUE)]])
  m <- merge(rd[, c("taxon","name","log2FoldChange","padj","sig_d")], a, by = "taxon")
  up <- m$sig_d & m$sig_a & m$log2FoldChange > 0
  list(rank = rank, n_taxa = nrow(m),
       deseq_sig = sum(m$sig_d), ancombc_sig = sum(m$sig_a),
       both = sum(m$sig_d & m$sig_a), deseq_only = sum(m$sig_d & !m$sig_a),
       ancombc_only = sum(!m$sig_d & m$sig_a),
       lfc_correlation = round(cor(m$log2FoldChange, m$lfc, use = "complete.obs"), 3),
       both_host_enriched = m$name[up][order(-m$log2FoldChange[up])])
}
out$da_phylum <- da_at("Phylum", 1000)
out$da_genus  <- da_at("Genus", 1000)
write_json(out, "results/compositional.json", auto_unbox = TRUE, pretty = TRUE, digits = 6)
cat("\n== Aitchison ==\n"); for (a in out$aitchison) cat(sprintf("  OTUs in >=%d samples: %d (%.1f%% reads)  habitat R2 %.3f p %.3f  sampletype R2 %.3f\n", a$min_prevalence, a$n_otu, 100*a$read_share, a$habitat_r2, a$habitat_p, a$sampletype_r2))
cat(sprintf("  betadisper p = %.4f\n", out$aitchison_betadisper_p))
for (d in list(out$da_phylum, out$da_genus)) cat(sprintf("\n== %s (%d taxa) ==\n  DESeq2 %d | ANCOM-BC %d | both %d | DESeq2-only %d | ANCOM-BC-only %d | lfc cor %.3f\n  both, host-enriched: %s\n",
  d$rank, d$n_taxa, d$deseq_sig, d$ancombc_sig, d$both, d$deseq_only, d$ancombc_only, d$lfc_correlation, paste(d$both_host_enriched, collapse=", ")))
