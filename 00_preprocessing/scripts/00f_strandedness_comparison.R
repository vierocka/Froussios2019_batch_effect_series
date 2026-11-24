# Strandedness check: gene counts with featureCounts -s 0 (fC/), -s 1 and -s 2, and -s 2 with
# --countReadPairs. Run from the repository root after 00e_count_strandedness.sh (pipeline_data/STAR_both/fC_*/).
# 1. assigned fraction per sample and variant (which strand setting fits the library?)
# 2. sanity check: fC_s0 must equal fC/
# 3. gene level: how many genes change between -s 0 and -s 2, and which ones
# 4. post 2's result (mean rlog correlation with the other 13 samples) on -s 0 vs. -s 2 counts
suppressMessages({ library(DESeq2); library(ggplot2); library(patchwork) })
work <- Sys.getenv("PIPELINE_DIR", file.path(getwd(), "pipeline_data"))
out  <- "00_preprocessing/results"
cnt  <- function(v) file.path(work, "STAR_both", v)
variants <- c("fC_s0", "fC_s1", "fC_s2", "fC_s2_pairs")
variants <- variants[dir.exists(cnt(variants))]
stopifnot(c("fC_s0", "fC_s2") %in% variants)

map <- read.delim("data/replicate_map.tsv")
pct <- read.csv("data/paper_TableS2C_rRNA.csv")
map$pct_rRNA <- pct$pct_rRNA[match(map$paper_replicate, pct$replicate)]
read_count <- function(f) { d <- read.delim(f, comment.char = "#"); setNames(d[[ncol(d)]], d$Geneid) }
mat <- function(v, runs) sapply(runs, function(r) read_count(file.path(cnt(v), paste0(r, ".count"))))

# 1. featureCounts summaries
summ <- do.call(rbind, lapply(variants, function(v) do.call(rbind, lapply(map$run, function(r) {
  s <- read.delim(file.path(cnt(v), paste0(r, ".count.summary")))
  n <- setNames(s[[2]], s$Status)
  data.frame(variant = v, run = r, replicate = map$paper_replicate[map$run == r],
             experiment = map$paper_experiment[map$run == r], total = sum(n), assigned = n[["Assigned"]],
             no_features = n[["Unassigned_NoFeatures"]], ambiguity = n[["Unassigned_Ambiguity"]])
}))))
summ$pct_assigned <- 100 * summ$assigned / summ$total
summ$pct_ambiguity <- 100 * summ$ambiguity / summ$total
write.csv(summ, file.path(out, "table00a_strandedness_assignment.csv"), row.names = FALSE)
cat("assigned % per variant (median and range over the 17 samples):\n")
print(aggregate(pct_assigned ~ variant, summ, function(x) round(c(median = median(x), min = min(x), max = max(x)), 1)))
cat("ambiguity % per variant (median):\n"); print(aggregate(pct_ambiguity ~ variant, summ, function(x) round(median(x), 2)))
if ("fC_s1" %in% variants) {
  w <- reshape(summ[summ$variant %in% c("fC_s1", "fC_s2"), c("variant", "run", "assigned")],
               idvar = "run", timevar = "variant", direction = "wide")
  cat(sprintf("assigned -s 2 / -s 1: median %.1f (range %.1f-%.1f); > 1 means reverse-stranded (read 1 antisense)\n",
              median(w$assigned.fC_s2 / w$assigned.fC_s1), min(w$assigned.fC_s2 / w$assigned.fC_s1),
              max(w$assigned.fC_s2 / w$assigned.fC_s1)))
}
if ("fC_s2_pairs" %in% variants) {
  r <- summ$assigned[summ$variant == "fC_s2_pairs"] / summ$assigned[summ$variant == "fC_s2"]
  cat(sprintf("assigned fragments (--countReadPairs) / assigned reads, -s 2: median %.3f (range %.3f-%.3f)\n",
              median(r), min(r), max(r)))
}

# 2. fC_s0 must reproduce fC/, the counts used by topics 02-07
pub <- mat("fC", map$run); s0 <- mat("fC_s0", map$run); s2 <- mat("fC_s2", map$run)
same <- identical(unname(pub), unname(s0)) && identical(rownames(pub), rownames(s0))
cat(sprintf("\nfC_s0 identical to fC/: %s\n", same))
if (!same) warning("fC_s0 differs from fC/: check featureCounts version and BAMs before comparing")

# 3. gene level, 14 ExpA/ExpB samples (as in post 2)
ab <- map$run[map$paper_experiment != "ExpC"]
keep <- rowSums(s0[, ab]) > 0 | rowSums(s2[, ab]) > 0
a0 <- s0[keep, ab]; a2 <- s2[keep, ab]
lr <- log2((a2 + 1) / (a0 + 1))
per_sample <- data.frame(run = ab, replicate = map$paper_replicate[match(ab, map$run)],
  total_s0 = colSums(a0), total_s2 = colSums(a2),
  genes_expressed_s0 = colSums(a0 >= 10),
  genes_change_10pct = colSums(a0 >= 10 & abs(lr) > log2(1.1)),
  genes_change_2fold = colSums(a0 >= 10 & abs(lr) > 1),
  genes_gain_2fold   = colSums(a0 >= 10 & lr > 1))
cat("\nper sample, genes with >= 10 reads at -s 0 that change between -s 0 and -s 2:\n")
print(per_sample, row.names = FALSE)
m0 <- rowMeans(a0); m2 <- rowMeans(a2)
gene_ratio <- data.frame(gene = rownames(a0), mean_s0 = m0, mean_s2 = m2,
                         log2_s2_vs_s0 = log2((m2 + 1) / (m0 + 1)))
gene_ratio <- gene_ratio[order(gene_ratio$log2_s2_vs_s0), ]
write.csv(gene_ratio, file.path(out, "table00b_strandedness_gene_ratio.csv"), row.names = FALSE)
g <- log2((m2 + 1) / (m0 + 1)); e <- m0 >= 10
cat(sprintf("\ngenes (mean >= 10 reads at -s 0): %d; lose > 10%%: %d; lose > 2-fold: %d; gain > 10%%: %d; gain > 2-fold: %d\n",
            sum(e), sum(e & g < -log2(1.1)), sum(e & g < -1), sum(e & g > log2(1.1)), sum(e & g > 1)))
cat("genes losing most (antisense reads counted on them when unstranded):\n")
print(head(gene_ratio[gene_ratio$mean_s0 >= 10, ], 10), row.names = FALSE)
cat("genes gaining most (unstranded counting discarded their reads as ambiguous with an antisense neighbour):\n")
print(head(gene_ratio[gene_ratio$mean_s2 >= 10 & gene_ratio$log2_s2_vs_s0 > 0, ][
  order(-gene_ratio$log2_s2_vs_s0[gene_ratio$mean_s2 >= 10 & gene_ratio$log2_s2_vs_s0 > 0]), ], 10), row.names = FALSE)

# 4. post 2's panel B (genes only, 4 rRNA rows not added) on -s 0 vs. -s 2
smp <- map[map$paper_experiment != "ExpC", ]; smp <- smp[order(smp$pct_rRNA), ]
mean_cor <- function(m) {
  m <- m[rowSums(m) > 0, smp$run]; storage.mode(m) <- "integer"
  rl <- assay(rlog(DESeqDataSetFromMatrix(m, data.frame(row.names = colnames(m), x = rep(1, ncol(m))), ~1), blind = TRUE))
  cm <- cor(rl); diag(cm) <- NA; rowMeans(cm, na.rm = TRUE)
}
cors <- data.frame(run = smp$run, replicate = smp$paper_replicate, experiment = smp$paper_experiment,
                   pct_rRNA = smp$pct_rRNA, s0 = mean_cor(s0), s2 = mean_cor(s2))
write.csv(cors, file.path(out, "table00c_strandedness_topic02_correlation.csv"), row.names = FALSE)
cat(sprintf("\npost 2 mean rlog correlation: -s 0 %.4f-%.4f, -s 2 %.4f-%.4f; Spearman (rank of samples) %.2f\n",
            min(cors$s0), max(cors$s0), min(cors$s2), max(cors$s2), cor(cors$s0, cors$s2, method = "spearman")))
cat(sprintf("Spearman with rRNA %%: -s 0 %.2f, -s 2 %.2f\n",
            cor(cors$s0, cors$pct_rRNA, method = "spearman"), cor(cors$s2, cors$pct_rRNA, method = "spearman")))

# Figure
lab <- c(fC_s0 = "-s 0 (unstranded, fC/)", fC_s1 = "-s 1", fC_s2 = "-s 2", fC_s2_pairs = "-s 2, fragments")
pal <- c(fC_s0 = "#333333", fC_s1 = "#66A61E", fC_s2 = "#D95F0E", fC_s2_pairs = "#2C7FB8")
dA <- summ[summ$variant %in% c("fC_s0", "fC_s1", "fC_s2"), ]
dA$rep <- factor(dA$replicate, levels = sort(unique(dA$replicate)))
pA <- ggplot(dA, aes(rep, pct_assigned, colour = variant, group = variant)) +
  geom_line() + geom_point(size = 2) +
  scale_colour_manual(values = pal, labels = lab, name = NULL) +
  labs(x = "replicate (paper numbering)", y = "% of reads assigned to a gene",
       title = "A. Which strand setting fits the library?") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

dB <- gene_ratio[gene_ratio$mean_s0 >= 1 | gene_ratio$mean_s2 >= 1, ]
pB <- ggplot(dB, aes(log2(mean_s0 + 1), log2_s2_vs_s0)) +
  geom_point(size = 0.4, alpha = 0.3, colour = "grey30") +
  geom_hline(yintercept = c(-1, 0, 1), linetype = c("dotted", "solid", "dotted"), colour = "grey50") +
  labs(x = "log2 mean count at -s 0 (+1), 14 samples", y = "log2 (-s 2 / -s 0)",
       title = "B. Per gene: stranded vs. unstranded counts",
       subtitle = sprintf("%d genes with >= 10 reads; %d change more than 2-fold", sum(m0 >= 10),
                          sum(m0 >= 10 & abs(log2((m2 + 1) / (m0 + 1))) > 1))) +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))

dC <- rbind(data.frame(cors[, c("run", "replicate", "pct_rRNA")], counts = "-s 0 (fC/)", mean_cor = cors$s0),
            data.frame(cors[, c("run", "replicate", "pct_rRNA")], counts = "-s 2", mean_cor = cors$s2))
dC$rep <- factor(dC$replicate, levels = cors$replicate)
pC <- ggplot(dC, aes(rep, mean_cor, colour = counts, group = counts)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  scale_colour_manual(values = c("-s 0 (fC/)" = "#333333", "-s 2" = "#D95F0E"), name = NULL) +
  labs(x = "replicate, ordered by rRNA % (as in post 2)", y = "mean rlog correlation\nwith the other 13 samples",
       title = "C. Post 2 panel B on unstranded vs. stranded counts") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

p <- (pA | pB) / pC + plot_annotation(
  title = "Froussios et al. 2019, TruSeq Stranded Total RNA: featureCounts -s 0 vs. -s 1 vs. -s 2",
  theme = theme(plot.title = element_text(size = 12, face = "bold")))
ggsave(file.path(out, "figure00a_strandedness.png"), p, width = 13, height = 10, dpi = 300, bg = "white")
cat("written: figure00a, table00a-c in", out, "\n")
