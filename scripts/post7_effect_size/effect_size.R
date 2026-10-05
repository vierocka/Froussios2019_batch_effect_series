# Post 7: does normalization only move genes across padj < 0.05, or also shift effect sizes?
# Run from the repository root. Same 8 splits and 4 size-factor treatments as post 5; only the
# rule that calls a gene "DE" changes.
suppressMessages({ library(DESeq2); library(ggplot2); library(patchwork) })
work <- Sys.getenv("WORK_DIR", file.path(getwd(), "work"))
out <- "scripts/post7_effect_size"
set.seed(1)

# identical to post 5: coding-gene counts, synthetic rRNA row, splits
map <- read.delim("data/replicate_map.tsv"); map <- map[map$paper_experiment != "ExpC", ]
pct <- read.csv("data/paper_TableS2C_rRNA.csv")
map$p <- pct$pct_rRNA[match(map$paper_replicate, pct$replicate)] / 100
read_count <- function(f) { d <- read.delim(f, comment.char = "#"); setNames(d[[ncol(d)]], d$Geneid) }
coding <- sapply(map$run, function(r) read_count(file.path(work, "STAR_both/fC", paste0(r, ".count"))))
coding <- coding[rowSums(coding) > 0, ]
rrna <- round(colSums(coding) * map$p / (1 - map$p))
with_rrna <- rbind(coding, rRNA = rrna)
expA <- map$run[map$paper_experiment == "ExpA"]; rep6 <- map$run[map$paper_replicate == 6]
others <- setdiff(expA, rep6)
splits <- lapply(sample(combn(length(others), 2, simplify = FALSE), 8), function(i) {
  g2 <- c(rep6, others[i]); c(sample(setdiff(others, g2), 3), g2) })

# fit once per split x treatment; DE is always tested on the coding-only matrix
fit <- function(norm_m, ids, total_norm) {
  cd <- data.frame(group = factor(rep(c("g1", "g2"), each = 3)), row.names = ids)
  dds <- DESeqDataSetFromMatrix(coding[, ids], cd, ~group)
  if (total_norm) { sf <- colSums(norm_m[, ids]); sf <- sf / mean(sf) }
  else sf <- sizeFactors(estimateSizeFactors(DESeqDataSetFromMatrix(norm_m[, ids], cd, ~1)))
  sizeFactors(dds) <- sf
  DESeq(dds, quiet = TRUE)
}
# call rules: significance, raw effect size, shrunken effect size, and a proper threshold test
calls <- function(dds) {
  r <- results(dds)                                     # MLE log2FC, H0: log2FC = 0
  s <- lfcShrink(dds, coef = "group_g2_vs_g1", type = "apeglm", quiet = TRUE)
  t1 <- results(dds, lfcThreshold = 1)                  # H0: |log2FC| <= 1
  lfc <- r$log2FoldChange; sh <- s$log2FoldChange
  c(`padj < 0.05`                 = sum(r$padj < 0.05, na.rm = TRUE),
    `|log2FC| > 1 (raw)`          = sum(abs(lfc) > 1, na.rm = TRUE),
    `|log2FC| > 2 (raw)`          = sum(abs(lfc) > 2, na.rm = TRUE),
    `|log2FC| > 1 (shrunken)`     = sum(abs(sh) > 1, na.rm = TRUE),
    `|log2FC| > 2 (shrunken)`     = sum(abs(sh) > 2, na.rm = TRUE),
    `padj < 0.05 & |log2FC| > 1`  = sum(r$padj < 0.05 & abs(lfc) > 1, na.rm = TRUE),
    `test |log2FC| > 1, padj < 0.05` = sum(t1$padj < 0.05, na.rm = TRUE))
}
treat <- data.frame(normalisation = rep(c("DESeq2 median-of-ratios", "total-count scaling"), each = 2),
                    rRNA = rep(c("included for size factors", "ignored for size factors"), 2),
                    total = rep(c(FALSE, TRUE), each = 2), with = rep(c(TRUE, FALSE), 2))
gsf <- function(dds) { l <- log2(sizeFactors(dds)); mean(l[4:6]) - mean(l[1:3]) }
res <- list(); shift <- list(); lowc <- list()
for (k in seq_along(splits)) {
  ids <- splits[[k]]
  fits <- lapply(seq_len(nrow(treat)), function(j) fit(if (treat$with[j]) with_rrna else coding, ids, treat$total[j]))
  # who passes a raw log2FC cutoff in a split with no designed difference? (median-of-ratios, rRNA ignored)
  r0 <- results(fits[[2]]); ok <- !is.na(r0$log2FoldChange)
  for (t in 1:2) { x <- ok & abs(r0$log2FoldChange) > t
    lowc[[length(lowc) + 1]] <- data.frame(split = k, cutoff = t, n = sum(x), median_baseMean = median(r0$baseMean[x]),
      pct_baseMean_below_10 = 100 * mean(r0$baseMean[x] < 10), median_baseMean_all = median(r0$baseMean[ok])) }
  for (j in seq_len(nrow(treat))) {
    n <- calls(fits[[j]])
    res[[length(res) + 1]] <- data.frame(split = k, treat[j, 1:2], rule = names(n), n = unname(n), row.names = NULL)
  }
  # how far does ignoring the rRNA feature move every gene's log2FC? (per normalisation)
  for (m in c(1, 3)) {
    d <- results(fits[[m + 1]])$log2FoldChange - results(fits[[m]])$log2FoldChange
    shift[[length(shift) + 1]] <- data.frame(split = k, normalisation = treat$normalisation[m],
      median_shift = median(d, na.rm = TRUE), iqr_shift = IQR(d, na.rm = TRUE),
      # expected shift: change in the between-group difference of mean log2 size factors
      sf_shift = gsf(fits[[m]]) - gsf(fits[[m + 1]]))
  }
}
res <- do.call(rbind, res); shift <- do.call(rbind, shift); lowc <- do.call(rbind, lowc)
write.csv(res, file.path(out, "effect_size_results.csv"), row.names = FALSE)
write.csv(shift, file.path(out, "log2FC_shift.csv"), row.names = FALSE)

# summary: per rule and normalisation, in how many splits does ignoring rRNA change the count?
w <- reshape(res, idvar = c("split", "normalisation", "rule"), timevar = "rRNA", direction = "wide")
w$changed <- w$`n.included for size factors` != w$`n.ignored for size factors`
cat("splits (of 8) whose count changes when the rRNA feature is ignored:\n")
print(xtabs(changed ~ rule + normalisation, w))
cat("\nrange of counts per rule and normalisation (min-max over splits and both rRNA treatments):\n")
print(aggregate(n ~ rule + normalisation, res, function(x) sprintf("%d-%d", min(x), max(x))))
cat("\nlog2FC shift when rRNA is ignored (median over genes; IQR over genes):\n")
print(transform(shift, median_shift = round(median_shift, 3), iqr_shift = signif(iqr_shift, 2), sf_shift = round(sf_shift, 3)))

cat("\ngenes passing a raw |log2FC| cutoff (median-of-ratios, rRNA ignored): how well expressed are they?\n")
print(transform(lowc, median_baseMean = round(median_baseMean, 1), pct_baseMean_below_10 = round(pct_baseMean_below_10),
                median_baseMean_all = round(median_baseMean_all)), row.names = FALSE)
write.csv(lowc, file.path(out, "raw_log2FC_expression.csv"), row.names = FALSE)

# Figure
rules <- c("padj < 0.05", "|log2FC| > 1 (raw)", "|log2FC| > 2 (raw)", "|log2FC| > 1 (shrunken)", "|log2FC| > 2 (shrunken)")
dA <- res[res$rule %in% rules, ]
dA$rule <- factor(dA$rule, levels = rules)
dA$rRNA <- factor(dA$rRNA, levels = c("included for size factors", "ignored for size factors"))
pA <- ggplot(dA, aes(rRNA, n, group = split)) + geom_line(colour = "grey60") +
  geom_point(size = 2, colour = "#2C7FB8") +
  facet_grid(normalisation ~ rule) + scale_y_sqrt() +
  scale_x_discrete(labels = c("included for size factors" = "incl.", "ignored for size factors" = "ign.")) +
  labs(x = "synthetic rRNA feature: source of size factors", y = "genes called (sqrt scale)",
       title = "A. Genes called under five rules (one line per split; same 8 splits as post 5)") +
  theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"), strip.text = element_text(size = 9))
dB <- shift[shift$normalisation == "total-count scaling", ]; dB$split <- factor(dB$split)
pB <- ggplot(dB, aes(split)) +
  geom_col(aes(y = median_shift), fill = "#D95F0E", width = 0.6) +
  geom_point(aes(y = sf_shift), size = 2.5, colour = "grey20") +
  geom_hline(yintercept = 0, colour = "grey40") +
  labs(x = "split", y = "change in log2FC\n(rRNA ignored minus included)",
       title = "B. Total-count scaling: ignoring the rRNA feature shifts every gene's log2FC by about the same amount",
       subtitle = "bar = median over genes (IQR over genes about 0.02); dot = shift expected from the size factors alone. Median-of-ratios: 0 in every split") +
  theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"))
p <- pA / pB + plot_layout(heights = c(1.6, 1))
ggsave(file.path(out, "effect_size.png"), p, width = 12, height = 9, dpi = 300, bg = "white")
