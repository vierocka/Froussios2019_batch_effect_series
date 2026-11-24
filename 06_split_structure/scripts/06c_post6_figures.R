# One single-panel figure per short post 6A-6D (volcano, rRNA swap, PCA, GO enrichment).
# Run from the repository root after 06a and 06b (6D reads table06c). Same counts and DESeq2
# settings as 06a_split_structure.R.
suppressMessages({ library(DESeq2); library(ggplot2) })
work <- Sys.getenv("PIPELINE_DIR", file.path(getwd(), "pipeline_data"))
out <- "06_split_structure/results"

map <- read.delim("data/replicate_map.tsv"); map <- map[map$paper_experiment == "ExpA", ]
pct <- read.csv("data/paper_TableS2C_rRNA.csv")
rd <- function(r) read.delim(file.path(work, "STAR_both/fC", paste0(r, ".count")), comment.char = "#")
coding <- sapply(map$run, function(r) { d <- rd(r); setNames(d[[ncol(d)]], d$Geneid) })
colnames(coding) <- paste0("r", map$paper_replicate)
coding <- coding[rowSums(coding) > 0, ]
pr <- setNames(pct$pct_rRNA[match(map$paper_replicate, pct$replicate)], colnames(coding))

fit <- function(g1, g2) {
  cd <- data.frame(group = factor(rep(c("g1", "g2"), each = 3)), row.names = c(g1, g2))
  m <- coding[, c(g1, g2)]; m <- m[rowSums(m) > 0, ]
  results(DESeq(DESeqDataSetFromMatrix(m, cd, ~group), quiet = TRUE))
}
n_de <- function(res) sum(res$padj < 0.05, na.rm = TRUE)

blue <- "#2C7FB8"; orange <- "#D95F0E"
th <- theme_minimal(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 17), plot.subtitle = element_text(colour = "grey30"),
        plot.caption = element_text(colour = "grey40", size = 10, hjust = 0), panel.grid.minor = element_blank(),
        plot.title.position = "plot", plot.caption.position = "plot", plot.margin = margin(16, 20, 12, 16))
save <- function(p, f) ggsave(file.path(out, f), p, width = 7, height = 7, dpi = 200, bg = "white")

# 6A: volcano of split 6, DE = padj < 0.05 only (no fold-change cutoff)
res6 <- fit(c("r1", "r2", "r3"), c("r4", "r5", "r6"))
saved <- read.csv(file.path(out, "table06a_split6_DE_genes.csv"))$gene
stopifnot(setequal(rownames(res6)[which(res6$padj < 0.05)], saved))
d <- as.data.frame(res6)[!is.na(res6$padj), ]
d$cls <- ifelse(d$padj >= 0.05, "not DE", ifelse(d$log2FoldChange > 0, "up", "down"))
de <- d[d$cls != "not DE", ]
small <- sum(abs(de$log2FoldChange) < 1)
pA <- ggplot(d, aes(log2FoldChange, -log10(padj))) +
  geom_point(data = d[d$cls == "not DE", ], colour = "grey80", size = 0.7) +
  geom_point(data = de, aes(colour = cls), size = 1.1, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey35") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", colour = "grey55") +
  annotate("text", x = -3.9, y = -log10(0.05), label = "padj = 0.05", vjust = -0.5, hjust = 0, size = 4, colour = "grey35") +
  annotate("text", x = 0, y = max(-log10(d$padj)) * 0.97, size = 4.2, colour = "grey20",
           label = sprintf("%d of %d DE genes\nchange less than 2-fold", small, nrow(de))) +
  scale_colour_manual(values = c(up = orange, down = blue), guide = "none") +
  coord_cartesian(xlim = c(-4, 4)) +
  labs(x = "log2 fold change (replicates 4, 5, 6 vs. 1, 2, 3)", y = "-log10 padj",
       title = sprintf("Identical plants, %d DE genes", nrow(de)),
       subtitle = "Same genotype, same batch, no designed difference.\nDE = padj < 0.05 only, no fold-change cutoff.",
       caption = "Froussios et al. 2019, ExpA; DESeq2 defaults. Dotted lines: 2-fold.") + th
save(pA, "figure06c_post6A_volcano.png")

# 6B: swap the high-rRNA replicate 6 for low-rRNA replicate 7
res7 <- fit(c("r1", "r2", "r3"), c("r4", "r5", "r7"))
dB <- data.frame(split = c(sprintf("1, 2, 3 vs. 4, 5, 6\nreplicate 6: %.1f%% rRNA", pr["r6"]),
                           sprintf("1, 2, 3 vs. 4, 5, 7\nreplicate 7: %.2f%% rRNA", pr["r7"])),
                 n = c(n_de(res6), n_de(res7)))
dB$split <- factor(dB$split, levels = dB$split)
pB <- ggplot(dB, aes(split, n)) +
  geom_col(width = 0.55, fill = c(orange, "grey55")) +
  geom_text(aes(label = format(n, big.mark = ",")), vjust = -0.4, size = 7, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "DE genes (padj < 0.05)", title = "Is it the rRNA?",
       subtitle = "Swap the highest-rRNA replicate for a low-rRNA one:\nDE genes go up, not down.",
       caption = "Froussios et al. 2019, ExpA; rRNA % from the paper's Table S2C; DESeq2 defaults.") +
  th + theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(size = 14, colour = "grey10"))
save(pB, "figure06d_post6B_rRNA_swap.png")

# 6C: PCA of the 7 replicates (vst, 500 most variable genes), as in 06a_split_structure.R
v <- assay(vst(DESeqDataSetFromMatrix(coding, data.frame(x = rep(1, 7), row.names = colnames(coding)), ~1), blind = TRUE))
pc <- prcomp(t(v[order(-apply(v, 1, var))[1:500], ]))
ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2))
dC <- data.frame(rep = sub("r", "", colnames(coding)), PC1 = pc$x[, 1], PC2 = pc$x[, 2])
if (dC$PC1[dC$rep == "1"] > 0) dC$PC1 <- -dC$PC1   # orient 1 -> 7 left to right
dC$grp <- ifelse(dC$rep %in% c("1", "2", "3"), "g1", ifelse(dC$rep %in% c("4", "5", "6"), "g2", "out"))
res_mix <- fit(c("r1", "r2", "r7"), c("r4", "r5", "r6"))
pC <- ggplot(dC, aes(PC1, PC2)) +
  geom_point(aes(fill = grp), shape = 21, size = 11, colour = "white", stroke = 1) +
  geom_text(aes(label = rep), colour = "white", fontface = "bold", size = 5.5) +
  annotate("segment", x = min(dC$PC1), xend = max(dC$PC1), y = min(dC$PC2) - 4, yend = min(dC$PC2) - 4,
           arrow = arrow(length = unit(0.25, "cm")), colour = "grey45") +
  annotate("text", x = mean(range(dC$PC1)), y = min(dC$PC2) - 4, vjust = 1.6, size = 4.2, colour = "grey30",
           label = sprintf("1, 2, 3 vs. 4, 5, 7 (across the order): %s DE genes\n1, 2, 7 vs. 4, 5, 6 (mixed): %d DE genes",
                           format(n_de(res7), big.mark = ","), n_de(res_mix))) +
  scale_fill_manual(values = c(g1 = blue, g2 = orange, out = "grey55"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.25, 0.12))) +
  labs(x = sprintf("PC1 (%d%% of variance)", ve[1]), y = sprintf("PC2 (%d%%)", ve[2]),
       title = "A hidden order",
       subtitle = "The replicates line up almost by their number.\nBlue / orange: the two groups of split 6.",
       caption = "Froussios et al. 2019, ExpA, 7 replicates; vst, 500 most variable genes.") + th
save(pC, "figure06e_post6C_PCA.png")

# 6D: GO enrichment, top terms with distinct gene sets; transcription factor term highlighted
goea <- read.csv(file.path(out, "table06c_GO_enrichment_split6.csv"))
goea <- goea[goea$FDR < 0.05 & !duplicated(goea$DE_genes), ]
dD <- head(goea, 8)
dD$label <- sprintf("%s  (%d)", dD$name, dD$n_de)
dD$label <- factor(dD$label, levels = rev(dD$label))
pD <- ggplot(dD, aes(-log10(FDR), label, fill = go == "GO:0003700")) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(`TRUE` = orange, `FALSE` = "grey65"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "-log10 FDR", y = NULL, title = "What are these genes?",
       subtitle = "Transcription factors are the top GO term,\nin a comparison with no designed difference.",
       caption = sprintf("GO enrichment of the %d DE genes (hypergeometric, BH).\nIn brackets: DE genes in the term.",
                         length(saved))) +
  th + theme(axis.text.y = element_text(size = 12, colour = "grey10"), panel.grid.major.y = element_blank())
save(pD, "figure06f_post6D_GO.png")

cat(sprintf("6A: %d DE (%d < 2-fold); 6B: %d vs %d; 6C: PC1 %d%%, mix %d; written figure06c-f\n",
            nrow(de), small, n_de(res6), n_de(res7), ve[1], n_de(res_mix)))
