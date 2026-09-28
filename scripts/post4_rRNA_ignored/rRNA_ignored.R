# Post 4: what changes if rRNA loci are ignored? Run from the repository root.
suppressMessages({ library(DESeq2); library(ggplot2); library(patchwork) })
work <- Sys.getenv("WORK_DIR", file.path(getwd(), "work"))
set.seed(1)

# coding-gene counts (featureCounts -t gene has no rRNA rows) for the 14 ExpA/ExpB samples
map <- read.delim("data/replicate_map.tsv"); map <- map[map$paper_experiment != "ExpC", ]
pct <- read.csv("data/paper_TableS2C_rRNA.csv")
map$p <- pct$pct_rRNA[match(map$paper_replicate, pct$replicate)] / 100
read_count <- function(f) { d <- read.delim(f, comment.char = "#"); setNames(d[[ncol(d)]], d$Geneid) }
coding <- sapply(map$run, function(r) read_count(file.path(work, "STAR_both/fC", paste0(r, ".count"))))
coding <- coding[rowSums(coding) > 0, ]

# SYNTHETIC rRNA row, not a recovered read count: how much rRNA would have to be added to the
# counted coding-gene total for it to represent the paper's published fraction of coding+rRNA?
# rrna = coding_total * p/(1-p)  =>  rrna/(coding_total+rrna) = p (solve to check). Our own current-
# annotation recount is known to badly undercount real rRNA reads (see NOTES), so this sensitivity
# analysis uses the paper's published Table S2C fraction instead of a real per-locus count.
rrna <- round(colSums(coding) * map$p / (1 - map$p))
with_rrna <- rbind(coding, rRNA = rrna)

# Panel A: one high-rRNA sample (rep 11); per-gene SHARE of the library with and without the
# synthetic rRNA row in the total. Raw coding-gene counts never change -- only their computed
# share of the (smaller or larger) total does.
s <- map$run[map$paper_replicate == 11]; p11 <- map$p[map$paper_replicate == 11]
x <- coding[, s]; x <- x[x > 0]
dA <- data.frame(share = c(x / (sum(x) + rrna[s]), x / sum(x)),
                 total = rep(c(sprintf("rRNA in library total (rRNA = %.0f%%)", 100 * p11), "rRNA ignored (coding genes only)"), each = length(x)))
pA <- ggplot(dA, aes(log2(share), colour = total)) + geom_density(linewidth = 1) +
  scale_colour_manual(values = c("#D95F0E", "#2C7FB8"), name = NULL) +
  labs(x = "log2 share of library, per coding gene", y = "density",
       title = sprintf("A. One library (rep 11, %.0f%% rRNA, synthetic):\nall coding-gene library shares shift by %.2f log2", 100 * p11, log2(1 / (1 - p11)))) +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 11)) +
  guides(colour = guide_legend(nrow = 2))

# Panel B: same batch (ExpA), same genotype -- no DESIGNED group difference, so any padj < 0.05
# call is false relative to the randomized labels (this does not rule out real biological/technical
# heterogeneity among the replicates themselves -- rep 6, 23.7% rRNA, is in every split below).
#
# Isolates the normalisation effect from the fitted-gene-set effect: DE is always tested on the
# SAME coding-only matrix; only the size factors are estimated from with_rrna vs. coding-only (or
# from their column sums, for total-count scaling), then assigned onto the coding-only DESeqDataSet.
# This way any difference between "in the matrix" and "ignored" can only come from the size factors,
# not from DESeq2's dispersion-trend fitting or independent filtering seeing a different gene set.
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
splits <- lapply(sample(combn(length(others), 2, simplify = FALSE), 8), function(i) {
  g2 <- c(rep6, others[i]); c(sample(setdiff(others, g2), 3), g2) })  # g1 = 3 of the remaining 4; 1 of 7 unused
dB <- do.call(rbind, lapply(seq_along(splits), function(k) {
  ids <- splits[[k]]
  data.frame(split = k, normalisation = rep(c("DESeq2 median-of-ratios", "total-count scaling"), each = 2),
             rRNA = rep(c("in the matrix", "ignored"), 2),
             n = c(n_de(coding, with_rrna, ids, FALSE), n_de(coding, coding, ids, FALSE),
                   n_de(coding, with_rrna, ids, TRUE),  n_de(coding, coding, ids, TRUE)))
}))
print(dB)
write.csv(dB, "scripts/post4_rRNA_ignored/rRNA_ignored_panelB_results.csv", row.names = FALSE)
dB$rRNA <- factor(dB$rRNA, levels = c("in the matrix", "ignored"))
pB <- ggplot(dB, aes(rRNA, n, group = split)) + geom_line(colour = "grey60") + geom_point(size = 2, colour = "#2C7FB8") +
  facet_wrap(~normalisation) + scale_y_sqrt() +
  labs(x = "rRNA (source of size factors; DE always tested on the same coding-gene matrix)",
       y = "false-relative-to-label DE genes (padj < 0.05, sqrt scale)",
       title = "B. Same batch and genotype, no designed group difference:\n8 sampled 3-vs-3 splits from the 7 ExpA replicates (one line per split)") +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold", size = 11))

ggsave("scripts/post4_rRNA_ignored/rRNA_ignored.png", pA + pB, width = 12, height = 5, dpi = 300, bg = "white")
