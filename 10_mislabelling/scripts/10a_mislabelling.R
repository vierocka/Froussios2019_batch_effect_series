# Mislabelling: what happens to DE when one sample per group is swapped? Run from the repository root after
# 08a_tin_cluster.sh (TIN for ExpA and ExpB in pipeline_data/TIN/).
# Samples: the 3 ExpA and 3 ExpB replicates with the highest TIN, ranked ExpA1-3 and ExpB1-3 by TIN.
#   mislabelled: group 1 = ExpA1, ExpA2 + ExpB3; group 2 = ExpB1, ExpB2 + ExpA3
#   correct:     group 1 = ExpA1-3;            group 2 = ExpB1-3
# The batch is ignored on purpose (~ group): we pretend that ExpA and ExpB are two conditions.
# Sample-level QC on raw counts of genes with > 2 reads in every sample (log2): replicate-replicate (R-R)
# plots of all 15 pairs, Pearson correlation matrix, and a dendrogram of Euclidean distances.
suppressMessages({ library(DESeq2); library(ggplot2); library(patchwork) })
work <- Sys.getenv("PIPELINE_DIR", file.path(getwd(), "pipeline_data"))
out <- "10_mislabelling/results"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

map <- read.delim("data/replicate_map.tsv"); map <- map[map$paper_experiment %in% c("ExpA", "ExpB"), ]
map$TIN <- sapply(map$run, function(r) read.delim(file.path(work, "TIN", paste0(r, "_Aligned.sortedByCoord.out.summary.txt")))[1, 3])
map <- map[order(map$paper_experiment, -map$TIN), ]
sel <- do.call(rbind, lapply(split(map, map$paper_experiment), head, 3))
sel$label <- paste0(sel$paper_experiment, ave(seq_len(nrow(sel)), sel$paper_experiment, FUN = seq_along), " (r", sel$paper_replicate, ")")
sel$id <- sub(" .*", "", sel$label)
sel$mislabelled_group <- ifelse(sel$id %in% c("ExpA1", "ExpA2", "ExpB3"), "group 1", "group 2")
sel$correct_group <- ifelse(sel$paper_experiment == "ExpA", "group 1", "group 2")
write.csv(sel[, c("id", "label", "run", "paper_replicate", "paper_experiment", "TIN", "mislabelled_group", "correct_group")],
          file.path(out, "table10a_samples.csv"), row.names = FALSE)
print(sel[, c("id", "run", "paper_replicate", "TIN", "mislabelled_group", "correct_group")], row.names = FALSE)

rd <- function(r) { d <- read.delim(file.path(work, "STAR_both/fC", paste0(r, ".count")), comment.char = "#")
  setNames(d[[ncol(d)]], d$Geneid) }
counts <- sapply(sel$run, rd); colnames(counts) <- sel$id

# DE: mislabelled vs. correct grouping ----------------------------------------------------------------------
de <- function(groups) {
  m <- counts[rowSums(counts) > 0, ]
  cd <- data.frame(group = factor(ifelse(groups == "group 1", "g1", "g2")), row.names = colnames(m))
  res <- results(DESeq(DESeqDataSetFromMatrix(m, cd, ~ group), quiet = TRUE), contrast = c("group", "g2", "g1"))
  data.frame(gene = rownames(res), baseMean = res$baseMean, log2FC = res$log2FoldChange, padj = res$padj)
}
d_mis <- de(sel$mislabelled_group); d_cor <- de(sel$correct_group)
sig <- function(d, lfc = 0) d$gene[which(d$padj < 0.05 & abs(d$log2FC) > lfc)]
cnt <- data.frame(grouping = c("mislabelled (one sample swapped per group)", "correct (ExpA vs. ExpB)"),
                  group1 = c(paste(sel$id[sel$mislabelled_group == "group 1"], collapse = ", "), paste(sel$id[sel$correct_group == "group 1"], collapse = ", ")),
                  group2 = c(paste(sel$id[sel$mislabelled_group == "group 2"], collapse = ", "), paste(sel$id[sel$correct_group == "group 2"], collapse = ", ")),
                  DE_padj05 = c(length(sig(d_mis)), length(sig(d_cor))),
                  DE_padj05_lfc1 = c(length(sig(d_mis, 1)), length(sig(d_cor, 1))),
                  tested = c(sum(!is.na(d_mis$padj)), sum(!is.na(d_cor$padj))))
cnt$mislabelled_DE_also_in_correct <- c(length(intersect(sig(d_mis), sig(d_cor))), NA)
write.csv(cnt, file.path(out, "table10b_DE_counts.csv"), row.names = FALSE)
print(cnt, row.names = FALSE)
write.csv(d_mis[which(d_mis$padj < 0.05), ], file.path(out, "table10c_DE_genes_mislabelled.csv"), row.names = FALSE)
write.csv(d_cor[which(d_cor$padj < 0.05), ], file.path(out, "table10d_DE_genes_correct.csv"), row.names = FALSE)

# sample QC on raw counts -----------------------------------------------------------------------------------
keep <- rowSums(counts > 2) == ncol(counts)
lg <- log2(counts[keep, ])
cat(sprintf("\n%d genes with > 2 reads in all 6 samples\n", sum(keep)))
ids <- sel$id
pairs <- t(combn(ids, 2))
ptype <- function(a, b) ifelse(substr(a, 1, 4) == substr(b, 1, 4), paste("within", substr(a, 1, 4)), "ExpA vs. ExpB")
pr <- data.frame(x = pairs[, 1], y = pairs[, 2], type = ptype(pairs[, 1], pairs[, 2]))
pr$r <- mapply(function(a, b) cor(lg[, a], lg[, b]), pr$x, pr$y)
pr$type <- factor(pr$type, c("within ExpA", "within ExpB", "ExpA vs. ExpB"))
pr <- pr[order(pr$type), ]
write.csv(pr, file.path(out, "table10e_pairwise_correlation.csv"), row.names = FALSE)

cols <- c("within ExpA" = "#2C7FB8", "within ExpB" = "#D95F0E", "ExpA vs. ExpB" = "#636363")
lim <- range(lg)
rr <- lapply(seq_len(nrow(pr)), function(i) {
  d <- data.frame(x = lg[, pr$x[i]], y = lg[, pr$y[i]])
  ggplot(d, aes(x, y)) + geom_point(size = 0.15, alpha = 0.15, colour = cols[as.character(pr$type[i])]) +
    geom_abline(colour = "black", linewidth = 0.3, linetype = "dashed") +
    annotate("text", x = lim[1], y = lim[2], hjust = 0, vjust = 1, size = 3, label = sprintf("r = %.4f", pr$r[i])) +
    coord_equal(xlim = lim, ylim = lim) +
    labs(x = sel$label[sel$id == pr$x[i]], y = sel$label[sel$id == pr$y[i]], title = as.character(pr$type[i])) +
    theme_minimal(base_size = 8) + theme(plot.title = element_text(size = 8, colour = cols[as.character(pr$type[i])], face = "bold"))
})
pRR <- wrap_plots(rr, ncol = 5) + plot_annotation(
  title = "Replicate-replicate plots, raw counts (log2), genes with > 2 reads in all 6 samples",
  subtitle = sprintf("%d genes; row 1: within ExpA and within ExpB; rows 2-3: ExpA vs. ExpB", sum(keep)),
  theme = theme(plot.title = element_text(face = "bold")))
ggsave(file.path(out, "figure10a_RR_plots.png"), pRR, width = 15, height = 9.5, dpi = 300, bg = "white")

cm <- cor(lg)
hc <- hclust(dist(t(lg), method = "euclidean"), method = "average")
o <- hc$labels[hc$order]
cdf <- expand.grid(a = colnames(cm), b = colnames(cm)); cdf$r <- as.vector(cm)
cdf$a <- factor(cdf$a, o); cdf$b <- factor(cdf$b, rev(o))
pC <- ggplot(cdf, aes(a, b, fill = r)) + geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.4f", r)), size = 3) +
  scale_fill_gradient(low = "#FFF7EC", high = "#B30000", name = "Pearson r") + coord_equal() +
  labs(x = NULL, y = NULL, title = "Correlation of raw log2 counts", subtitle = "order of the dendrogram") +
  theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(out, "figure10b_correlation_matrix.png"), pC, width = 7, height = 6, dpi = 300, bg = "white")
write.csv(data.frame(sample = rownames(cm), round(cm, 5)), file.path(out, "table10f_correlation_matrix.csv"), row.names = FALSE)

png(file.path(out, "figure10c_dendrogram.png"), width = 7, height = 5.5, units = "in", res = 300)
par(mar = c(6, 4, 4, 1))
dg <- as.dendrogram(hc)
lab_col <- ifelse(substr(labels(dg), 1, 4) == "ExpA", "#2C7FB8", "#D95F0E")
dendextend_free <- function(n) { if (is.leaf(n)) attr(n, "nodePar") <- list(lab.col = lab_col[match(attr(n, "label"), labels(dg))], pch = NA); n }
plot(dendrapply(dg, dendextend_free), ylab = "Euclidean distance (raw log2 counts)",
     main = "Sample dendrogram (average linkage)\nshortest join = most similar samples")
dev.off()

# DE counts figure
cnt$grouping <- factor(cnt$grouping, cnt$grouping)
pD <- ggplot(cnt, aes(grouping, DE_padj05, fill = grouping)) + geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%d DE genes\n(%d with |log2FC| > 1)", DE_padj05, DE_padj05_lfc1)), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("#D95F0E", "#2C7FB8"), guide = "none") + scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = NULL, y = "DE genes (padj < 0.05)", title = "One swapped sample per group",
       subtitle = sprintf("mislabelled: %s vs. %s\ncorrect: %s vs. %s", cnt$group1[1], cnt$group2[1], cnt$group1[2], cnt$group2[2])) +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(out, "figure10d_DE_counts.png"), pD, width = 7, height = 5.5, dpi = 300, bg = "white")
