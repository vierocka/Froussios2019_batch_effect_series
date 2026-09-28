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

# add rRNA back as one row, at the paper's published fraction of each library (Table S2C)
rrna <- round(colSums(coding) * map$p / (1 - map$p))
with_rrna <- rbind(coding, rRNA = rrna)

# Panel A: one high-rRNA sample (rep 11); per-gene share of the library with and without rRNA in the total
s <- map$run[map$paper_replicate == 11]; p11 <- map$p[map$paper_replicate == 11]
x <- coding[, s]; x <- x[x > 0]
dA <- data.frame(share = c(x / (sum(x) + rrna[s]), x / sum(x)),
                 total = rep(c(sprintf("rRNA in library total (rRNA = %.0f%%)", 100 * p11), "rRNA ignored (coding genes only)"), each = length(x)))
pA <- ggplot(dA, aes(log10(share), colour = total)) + geom_density(linewidth = 1) +
  scale_colour_manual(values = c("#D95F0E", "#2C7FB8"), name = NULL) +
  labs(x = "log10 share of library, per coding gene", y = "density",
       title = sprintf("A. One library (rep 11, %.0f%% rRNA):\nall coding genes shift by %.2f log2", 100 * p11, log2(1 / (1 - p11)))) +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 11)) +
  guides(colour = guide_legend(nrow = 2))

# Panel B: same batch (ExpA), 3 vs 3 replicates, rep 6 (23.7% rRNA) always in group 2; every DE call is a false positive
n_de <- function(m, run_ids, total_norm) {
  cd <- data.frame(group = factor(rep(c("g1", "g2"), each = 3)), row.names = run_ids)
  dds <- DESeqDataSetFromMatrix(m[, run_ids], cd, ~group)
  if (total_norm) sizeFactors(dds) <- colSums(counts(dds)) / mean(colSums(counts(dds))) else dds <- estimateSizeFactors(dds)
  res <- results(DESeq(dds, quiet = TRUE))
  sum(res$padj < 0.05 & rownames(res) != "rRNA", na.rm = TRUE)
}
expA <- map$run[map$paper_experiment == "ExpA"]; rep6 <- map$run[map$paper_replicate == 6]
others <- setdiff(expA, rep6)
splits <- lapply(sample(combn(length(others), 2, simplify = FALSE), 8), function(i) {
  g2 <- c(rep6, others[i]); c(sample(setdiff(others, g2), 3), g2) })
dB <- do.call(rbind, lapply(seq_along(splits), function(k) {
  ids <- splits[[k]]
  data.frame(split = k, normalisation = rep(c("DESeq2 median-of-ratios", "total-count scaling"), each = 2),
             rRNA = rep(c("in the matrix", "ignored"), 2),
             n = c(n_de(with_rrna, ids, FALSE), n_de(coding, ids, FALSE), n_de(with_rrna, ids, TRUE), n_de(coding, ids, TRUE)))
}))
print(dB)
dB$rRNA <- factor(dB$rRNA, levels = c("in the matrix", "ignored"))
pB <- ggplot(dB, aes(rRNA, n, group = split)) + geom_line(colour = "grey60") + geom_point(size = 2, colour = "#2C7FB8") +
  facet_wrap(~normalisation) + scale_y_sqrt() +
  labs(x = "rRNA", y = "false DE genes (padj < 0.05, sqrt scale)",
       title = "B. Same batch and genotype: 8 random 3 vs 3 splits\n(one line per split)") +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold", size = 11))

ggsave("scripts/post4_rRNA_ignored/rRNA_ignored.png", pA + pB, width = 12, height = 5, dpi = 300, bg = "white")
