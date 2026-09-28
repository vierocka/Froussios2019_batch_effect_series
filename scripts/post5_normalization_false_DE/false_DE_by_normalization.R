# Post 5: can normalization choice change false-DE counts? Run from the repository root.
suppressMessages({ library(DESeq2); library(ggplot2) })
work <- Sys.getenv("WORK_DIR", file.path(getwd(), "work"))
set.seed(1)

# coding-gene counts (featureCounts -t gene has no rRNA rows) for the 14 ExpA/ExpB samples
map <- read.delim("data/replicate_map.tsv"); map <- map[map$paper_experiment != "ExpC", ]
pct <- read.csv("data/paper_TableS2C_rRNA.csv")
map$p <- pct$pct_rRNA[match(map$paper_replicate, pct$replicate)] / 100
read_count <- function(f) { d <- read.delim(f, comment.char = "#"); setNames(d[[ncol(d)]], d$Geneid) }
coding <- sapply(map$run, function(r) read_count(file.path(work, "STAR_both/fC", paste0(r, ".count"))))
coding <- coding[rowSums(coding) > 0, ]

# SYNTHETIC rRNA row (see post 4): coding genes only, plus one calibrated aggregate rRNA feature
# added at the paper's published fraction. Never part of the tested DE matrix below -- used only
# as an alternative source of size factors.
rrna <- round(colSums(coding) * map$p / (1 - map$p))
with_rrna <- rbind(coding, rRNA = rrna)

# Same batch (ExpA), same genotype -- there is no DESIGNED group difference between the two
# labels, so a padj < 0.05 call is false relative to the randomized labels (this does not mean
# the replicates have no real biological or technical heterogeneity among themselves).
#
# Replicate 11 (31% rRNA, post 4) was excluded from the original study's own analysis. To test
# consequences among RETAINED samples from one experiment, replicate 6 is used here instead --
# the highest-rRNA replicate that was actually kept in ExpA (23.7% reported rRNA per Table S2C).
# It is deliberately placed in group 2 in every split below.
#
# DE is always tested on the SAME coding-only matrix; only the size factors are estimated from
# with_rrna vs. coding-only (or from their column sums, for total-count scaling), then assigned
# onto the coding-only DESeqDataSet. This isolates the normalisation effect: any difference
# between "included" and "ignored" can only come from the size factors, not from DESeq2's
# dispersion-trend fitting or independent filtering seeing a different set of rows.
n_de <- function(coding_m, norm_m, run_ids, total_norm) {
  cd <- data.frame(group = factor(rep(c("g1", "g2"), each = 3)), row.names = run_ids)
  test_dds <- DESeqDataSetFromMatrix(coding_m[, run_ids], cd, ~group)
  if (total_norm) {
    sf <- colSums(norm_m[, run_ids]); sf <- sf / mean(sf)
  } else {
    norm_dds <- estimateSizeFactors(DESeqDataSetFromMatrix(norm_m[, run_ids], cd, ~1))
    sf <- sizeFactors(norm_dds)
  }
  sizeFactors(test_dds) <- sf
  res <- results(DESeq(test_dds, quiet = TRUE))
  sum(res$padj < 0.05, na.rm = TRUE)
}
expA <- map$run[map$paper_experiment == "ExpA"]; rep6 <- map$run[map$paper_replicate == 6]
others <- setdiff(expA, rep6)  # 6 of the 7 ExpA replicates besides rep 6
# 8 sampled 3-vs-3 comparisons: eight of the 15 possible choices for rep 6's two group-2 partners,
# then randomly selecting three of the remaining four samples for group 1 (not an exhaustive
# enumeration of the 60 possible comparisons: choose(6,2)*choose(4,3) = 15*4 = 60).
splits <- lapply(sample(combn(length(others), 2, simplify = FALSE), 8), function(i) {
  g2 <- c(rep6, others[i]); c(sample(setdiff(others, g2), 3), g2) })
dB <- do.call(rbind, lapply(seq_along(splits), function(k) {
  ids <- splits[[k]]
  data.frame(split = k, normalisation = rep(c("DESeq2 median-of-ratios", "total-count scaling"), each = 2),
             rRNA = rep(c("included for size factors", "ignored for size factors"), 2),
             n = c(n_de(coding, with_rrna, ids, FALSE), n_de(coding, coding, ids, FALSE),
                   n_de(coding, with_rrna, ids, TRUE),  n_de(coding, coding, ids, TRUE)))
}))
print(dB)
write.csv(dB, "scripts/post5_normalization_false_DE/false_DE_by_normalization_results.csv", row.names = FALSE)
dB$rRNA <- factor(dB$rRNA, levels = c("included for size factors", "ignored for size factors"))
pB <- ggplot(dB, aes(rRNA, n, group = split)) + geom_line(colour = "grey60") + geom_point(size = 2.5, colour = "#2C7FB8") +
  facet_wrap(~normalisation) + scale_y_sqrt() +
  scale_x_discrete(labels = c("included for size factors" = "included", "ignored for size factors" = "ignored")) +
  labs(x = "synthetic rRNA feature: source of size factors",
       y = "false-relative-to-label DE genes\n(padj < 0.05, sqrt scale)",
       title = "Same batch and genotype, no designed group difference:\n8 sampled 3-vs-3 splits from the 7 ExpA replicates (one line per split)",
       subtitle = "Replicate 6 (23.7% rRNA, the highest-rRNA replicate retained in ExpA) is in every split") +
  theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold", size = 12))

ggsave("scripts/post5_normalization_false_DE/false_DE_by_normalization.png", pB, width = 9, height = 6.5, dpi = 300, bg = "white")
