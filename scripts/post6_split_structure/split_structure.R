# Post 6: why does one same-batch 3-vs-3 split give 689 DE genes? Run from the repository root.
suppressMessages({ library(DESeq2); library(ggplot2); library(patchwork) })
work <- Sys.getenv("WORK_DIR", file.path(getwd(), "work"))
cores <- as.integer(Sys.getenv("THREADS", "4"))
out <- "scripts/post6_split_structure"

# coding-gene counts (featureCounts -t gene) for the 14 ExpA/ExpB samples, filtered as in post 5
map <- read.delim("data/replicate_map.tsv"); map <- map[map$paper_experiment != "ExpC", ]
pct <- read.csv("data/paper_TableS2C_rRNA.csv")
map$pct <- pct$pct_rRNA[match(map$paper_replicate, pct$replicate)]
rd <- function(r) read.delim(file.path(work, "STAR_both/fC", paste0(r, ".count")), comment.char = "#")
chr <- with(rd(map$run[1]), setNames(Chr, Geneid))
coding <- sapply(map$run, function(r) { d <- rd(r); setNames(d[[ncol(d)]], d$Geneid) })
coding <- coding[rowSums(coding) > 0, ]
colnames(coding) <- paste0("r", map$paper_replicate)
expA <- paste0("r", map$paper_replicate[map$paper_experiment == "ExpA"])
pr <- setNames(map$pct, colnames(coding))

# DESeq2 with its default median-of-ratios size factors (rRNA never enters: post 5 showed it
# does not change this split under median-of-ratios)
de <- function(g1, g2) {
  cd <- data.frame(group = factor(rep(c("g1", "g2"), c(length(g1), length(g2)))), row.names = c(g1, g2))
  DESeq(DESeqDataSetFromMatrix(coding[, c(g1, g2)], cd, ~group), quiet = TRUE)
}
n_de <- function(dds) sum(results(dds)$padj < 0.05, na.rm = TRUE)

# Post 5's split 6: replicates 1, 2, 3 (group 1) vs. 4, 5, 6 (group 2); replicate 7 unused
g1 <- c("r1", "r2", "r3"); g2 <- c("r4", "r5", "r6")
dds <- de(g1, g2); res <- results(dds)
sig <- res[which(res$padj < 0.05), ]
cat(sprintf("split 6: %d DE genes (padj < 0.05), %d higher in group 2; %d nuclear, %d Pt, %d Mt\n",
            nrow(sig), sum(sig$log2FoldChange > 0), sum(!chr[rownames(sig)] %in% c("Pt", "Mt")),
            sum(chr[rownames(sig)] == "Pt"), sum(chr[rownames(sig)] == "Mt")))
cat(sprintf("median |log2FC| of DE genes %.2f (max %.2f); median baseMean DE %.0f vs. all tested %.0f\n",
            median(abs(sig$log2FoldChange)), max(abs(sig$log2FoldChange)),
            median(sig$baseMean), median(res$baseMean[!is.na(res$padj)])))
sig_df <- data.frame(gene = rownames(sig), chr = chr[rownames(sig)], as.data.frame(sig))
write.csv(sig_df[order(sig_df$padj), ], file.path(out, "split6_DE_genes.csv"), row.names = FALSE)

# Is it rRNA, replicate 6, or something shared by the replicates? (numbers for the README)
# 1. PCA of the 7 ExpA replicates
v <- assay(vst(DESeqDataSetFromMatrix(coding[, expA], data.frame(x = rep(1, 7), row.names = expA), ~1), blind = TRUE))
pc <- prcomp(t(v[order(-apply(v, 1, var))[1:500], ]))
cat(sprintf("\nPCA, 500 most variable genes: PC1 %.0f%%, PC2 %.0f%% of variance\nPC1 per replicate: %s\n",
            100 * pc$sdev[1]^2 / sum(pc$sdev^2), 100 * pc$sdev[2]^2 / sum(pc$sdev^2),
            paste(sprintf("%s %.1f", names(sort(pc$x[, 1])), sort(pc$x[, 1])), collapse = ", ")))
# 2. all 70 distinct 3-vs-3 splits of the 7 ExpA replicates (one left out)
splits <- list()
for (u in expA) { rest <- setdiff(expA, u); for (a in combn(rest, 3, simplify = FALSE)) {
  b <- setdiff(rest, a); k <- paste(sort(c(paste(a, collapse = "+"), paste(b, collapse = "+"))), collapse = " vs ")
  splits[[k]] <- list(a, b) } }
all70 <- do.call(rbind, parallel::mclapply(names(splits), function(k) {
  s <- splits[[k]]
  data.frame(split = k, n_de = n_de(de(s[[1]], s[[2]])),
             log_rRNA_separation = abs(mean(log(pr[s[[1]]])) - mean(log(pr[s[[2]]]))))
}, mc.cores = cores))
all70 <- all70[order(-all70$n_de), ]
write.csv(all70, file.path(out, "all_70_splits.csv"), row.names = FALSE)
cat(sprintf("\nall %d splits: median %g DE genes, 90th percentile %g, max %d (%s); split 6 ranks %d\n",
            nrow(all70), median(all70$n_de), quantile(all70$n_de, 0.9), all70$n_de[1], all70$split[1],
            which(all70$split == "r1+r2+r3 vs r4+r5+r6")))
cat(sprintf("Spearman rho, DE count vs. rRNA separation of the two groups: %.2f\n",
            cor(all70$n_de, all70$log_rRNA_separation, method = "spearman")))
# 3. leave one out, and swap replicate 6 or 3 for replicate 7
cat("\nvariants of split 6:\n")
for (v2 in list(list(g1, c("r4", "r5")), list(g1, c("r4", "r6")), list(g1, c("r5", "r6")),
                list(c("r2", "r3"), g2), list(c("r1", "r3"), g2), list(c("r1", "r2"), g2),
                list(g1, c("r4", "r5", "r7")), list(c("r1", "r2", "r7"), g2)))
  cat(sprintf("  %-10s vs %-10s %4d DE genes\n", paste(v2[[1]], collapse = ","), paste(v2[[2]], collapse = ","),
              n_de(de(v2[[1]], v2[[2]]))))
# 4. per-sample: how strongly does each replicate follow the group-2 direction on the DE genes?
nc <- counts(estimateSizeFactors(DESeqDataSetFromMatrix(coding[, expA], data.frame(x = rep(1, 7), row.names = expA), ~1)), normalized = TRUE)
l <- log2(nc[rownames(sig), ] + 1); l <- (l - rowMeans(l[, g1])) * sign(sig$log2FoldChange)
cat("\nmean signed log2 deviation from the group-1 mean on the DE genes (positive = group-2 direction):\n")
print(round(colMeans(l), 2))

# Figure
nc6 <- counts(dds, normalized = TRUE)[!is.na(res$padj), ]
med <- data.frame(gene = rownames(nc6),
                  g1 = apply(nc6[, g1], 1, median), g2 = apply(nc6[, g2], 1, median))
all_lab <- sprintf("all %s tested genes", format(nrow(med), big.mark = ","))
up <- rownames(sig)[sig$log2FoldChange > 0]; dn <- rownames(sig)[sig$log2FoldChange < 0]
dA <- rbind(transform(med, set = all_lab),
            transform(med[med$gene %in% up, ], set = sprintf("%d DE genes higher in replicates 4, 5, 6", length(up))),
            transform(med[med$gene %in% dn, ], set = sprintf("%d DE genes lower in replicates 4, 5, 6", length(dn))))
dA <- rbind(data.frame(set = dA$set, group = "replicates 1, 2, 3", median = dA$g1),
            data.frame(set = dA$set, group = "replicates 4, 5, 6", median = dA$g2))
dA$set <- factor(dA$set, levels = unique(dA$set))
pal <- c("replicates 1, 2, 3" = "#2C7FB8", "replicates 4, 5, 6" = "#D95F0E")
pA <- ggplot(dA, aes(log2(median + 1), colour = group)) + geom_density(linewidth = 1) +
  facet_wrap(~set, ncol = 1, scales = "free_y") + scale_colour_manual(values = pal, name = NULL) +
  coord_cartesian(xlim = c(2, 18)) +
  labs(x = "per-gene median of normalized counts in the group (log2 + 1)", y = "density",
       title = "A. Gene medians per group") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

dB <- data.frame(rep = c(g1, g2), pct = pr[c(g1, g2)],
                 group = rep(c("replicates 1, 2, 3", "replicates 4, 5, 6"), each = 3))
bars <- aggregate(pct ~ group, dB, median)
pB <- ggplot(dB, aes(group, pct, colour = group)) +
  geom_col(data = bars, aes(fill = group), colour = NA, alpha = 0.25, width = 0.6) +
  geom_point(size = 3.5) +
  ggrepel::geom_text_repel(aes(label = sub("r", "rep ", rep)), size = 3.5, direction = "y", nudge_x = 0.35,
                           colour = "grey20", min.segment.length = 0) +
  scale_colour_manual(values = pal, guide = "none") + scale_fill_manual(values = pal, guide = "none") +
  labs(x = NULL, y = "% rRNA (paper's Table S2C)", title = "B. rRNA per replicate",
       subtitle = "bar = group median") +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))

p <- pA + pB + plot_layout(widths = c(1.6, 1)) + plot_annotation(
  title = sprintf("Same genotype, same batch (ExpA), no designed difference: replicates 1, 2, 3 vs. 4, 5, 6 -> %d DE genes (padj < 0.05)", nrow(sig)),
  theme = theme(plot.title = element_text(size = 11, face = "bold")))
ggsave(file.path(out, "split_structure.png"), p, width = 11, height = 7.5, dpi = 300, bg = "white")
