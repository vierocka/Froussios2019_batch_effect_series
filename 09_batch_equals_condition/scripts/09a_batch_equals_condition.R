# When batch = condition: the 4 ExpA replicates with the highest TIN ("condition A") vs. the 4 ExpB
# replicates with the highest TIN ("condition B"). All 8 are Col-0 wild type, so every DE gene is a batch effect.
#   model 1, batch ignored:     DESeq2 ~ condition on raw counts
#   model 2, batch removed:     DESeq2 ~ batch + condition is not estimable (batch and condition are the same
#                               column); ComBat-seq removes the batch from the counts, then DESeq2 ~ condition
# Per model: DE genes and their GO enrichment, PCA (PC1 vs. PC2, PC2 vs. PC3), correlation of PC1-3 with batch,
# gene contributions to PC1-3 (top 20 of PC1, GO enrichment of these 20).
# Run from the repository root after 08a_tin_cluster.sh has been run for ExpA and ExpB (pipeline_data/TIN/).
# Also writes the tables read by the Shiny apps (results/app_data/).
suppressMessages({ library(DESeq2); library(sva); library(ggplot2); library(patchwork) })
source("lib/go_enrichment.R")
work <- Sys.getenv("PIPELINE_DIR", file.path(getwd(), "pipeline_data"))
out <- "09_batch_equals_condition/results"; app <- file.path(out, "app_data")
dir.create(app, showWarnings = FALSE, recursive = TRUE)

# samples: 4 best TIN (median transcript TIN) per experiment
map <- read.delim("data/replicate_map.tsv"); map <- map[map$paper_experiment %in% c("ExpA", "ExpB"), ]
map$TIN <- sapply(map$run, function(r) read.delim(file.path(work, "TIN", paste0(r, "_Aligned.sortedByCoord.out.summary.txt")))[1, 3])
map <- map[order(map$paper_experiment, -map$TIN), ]
sel <- do.call(rbind, lapply(split(map, map$paper_experiment), head, 4))
sel$sample <- paste0(sel$paper_experiment, "_r", sel$paper_replicate)
sel$batch <- sel$paper_experiment
sel$condition <- ifelse(sel$batch == "ExpA", "A", "B")
rownames(sel) <- sel$sample
write.csv(sel[, c("sample", "run", "paper_replicate", "batch", "condition", "TIN")],
          file.path(out, "table09a_samples.csv"), row.names = FALSE)
print(sel[, c("run", "paper_replicate", "batch", "condition", "TIN")])

rd <- function(r) { d <- read.delim(file.path(work, "STAR_both/fC", paste0(r, ".count")), comment.char = "#")
  setNames(d[[ncol(d)]], d$Geneid) }
counts <- sapply(sel$run, rd); colnames(counts) <- sel$sample
counts <- counts[rowSums(counts >= 10) >= 4, ]        # >= 10 reads in at least 4 of the 8 samples
cat(sprintf("%d genes with >= 10 reads in at least 4 samples\n", nrow(counts)))
cd <- sel[, c("batch", "condition")]; cd[] <- lapply(cd, factor)

# model 2: the batch + condition design is not estimable
msg <- tryCatch({ DESeqDataSetFromMatrix(counts, cd, ~ batch + condition); "estimable" }, error = function(e) conditionMessage(e))
cat("DESeq2 ~ batch + condition:", gsub("\\s+", " ", msg), "\n")
msg_cs <- tryCatch({ ComBat_seq(counts, batch = cd$batch, group = cd$condition); "ran" }, error = function(e) conditionMessage(e))
cat("ComBat-seq with condition as covariate:", gsub("\\s+", " ", msg_cs), "\n")
combat <- suppressMessages(ComBat_seq(counts, batch = cd$batch))         # batch only, as one must when batch = condition
storage.mode(combat) <- "integer"
write.csv(data.frame(design = c("DESeq2 ~ batch + condition", "ComBat-seq, batch with condition as covariate"),
                     result = gsub("\\s+", " ", c(msg, msg_cs))), file.path(out, "table09b_confounded_designs.csv"), row.names = FALSE)

go_universe <- rownames(counts)
info <- gene_info()
analyse <- function(m, label, tag) {
  dds <- DESeq(DESeqDataSetFromMatrix(m, cd, ~ condition), quiet = TRUE)
  res <- results(dds, contrast = c("condition", "B", "A"))
  de <- data.frame(gene = rownames(res), baseMean = res$baseMean, log2FC = res$log2FoldChange, lfcSE = res$lfcSE,
                   pvalue = res$pvalue, padj = res$padj)
  de$symbol <- info$symbol[match(de$gene, info$gene)]
  sig <- de$gene[which(de$padj < 0.05)]
  sig_strong <- de$gene[which(de$padj < 0.05 & abs(de$log2FC) > 1)]
  write.csv(de[which(de$padj < 0.05), ][order(de$padj[which(de$padj < 0.05)]), ],
            file.path(out, sprintf("table09c_DE_genes_%s.csv", tag)), row.names = FALSE)
  go_de <- if (length(sig_strong) >= 5) go_enrich(sig_strong, go_universe) else go_enrich(character(0), go_universe)[0, ]
  write.csv(go_de, file.path(out, sprintf("table09d_GO_DE_genes_%s.csv", tag)), row.names = FALSE)

  vs <- assay(vst(dds, blind = TRUE))
  pca <- prcomp(t(vs), center = TRUE, scale. = FALSE)
  ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  sc <- data.frame(sample = rownames(pca$x), pca$x[, 1:3], batch = cd$batch, TIN = sel$TIN)
  bnum <- as.numeric(cd$batch == "ExpB")
  pcb <- do.call(rbind, lapply(1:3, function(k) {
    t <- cor.test(pca$x[, k], bnum)
    data.frame(model = label, PC = paste0("PC", k), pct_variance = ve[k], r_with_batch = unname(t$estimate),
               R2 = unname(t$estimate)^2, p = t$p.value)
  }))
  # contribution of a gene to a PC = squared loading, in % (the loadings of a PC have unit length)
  contrib <- data.frame(gene = rownames(pca$rotation), 100 * pca$rotation[, 1:3]^2)
  names(contrib)[2:4] <- paste0("contrib_PC", 1:3)
  for (k in 1:3) contrib[[paste0("rank_PC", k)]] <- rank(-contrib[[paste0("contrib_PC", k)]], ties.method = "first")
  for (k in 1:3) contrib[[paste0("loading_PC", k)]] <- pca$rotation[, k]
  contrib$symbol <- info$symbol[match(contrib$gene, info$gene)]
  top20 <- contrib[order(contrib$rank_PC1), ][1:20, ]
  write.csv(top20[, c("gene", "symbol", "contrib_PC1", "loading_PC1")], file.path(out, sprintf("table09e_top20_PC1_genes_%s.csv", tag)), row.names = FALSE)
  go_pc <- go_enrich(top20$gene, go_universe, min_size = 5)
  write.csv(go_pc, file.path(out, sprintf("table09f_GO_top20_PC1_%s.csv", tag)), row.names = FALSE)

  # app tables
  nc <- counts(dds, normalized = TRUE)
  write.csv(data.frame(gene = rownames(nc), round(nc, 2), check.names = FALSE), file.path(app, sprintf("normalized_counts_%s.csv", tag)), row.names = FALSE)
  write.csv(de, file.path(app, sprintf("DE_all_genes_%s.csv", tag)), row.names = FALSE)
  write.csv(sc, file.path(app, sprintf("PCA_scores_%s.csv", tag)), row.names = FALSE)
  write.csv(data.frame(PC = paste0("PC", seq_along(ve)), pct_variance = ve), file.path(app, sprintf("PCA_variance_%s.csv", tag)), row.names = FALSE)
  write.csv(contrib, file.path(app, sprintf("PCA_contributions_%s.csv", tag)), row.names = FALSE)

  cat(sprintf("\n== %s\nDE genes padj < 0.05: %d (%.1f%% of %d tested); with |log2FC| > 1: %d\n", label, length(sig),
              100 * length(sig) / sum(!is.na(de$padj)), sum(!is.na(de$padj)), length(sig_strong)))
  cat(sprintf("GO terms FDR < 0.05 for DE genes: %d; for top 20 PC1 genes: %d\n", sum(go_de$FDR < 0.05), sum(go_pc$FDR < 0.05)))
  print(pcb, row.names = FALSE, digits = 3)
  list(label = label, n_de = length(sig), n_strong = length(sig_strong), sc = sc, ve = ve, pcb = pcb, top20 = top20,
       go_de = go_de, go_pc = go_pc)
}
r1 <- analyse(counts, "batch ignored (~ condition)", "batch_ignored")
r2 <- analyse(combat, "batch removed (ComBat-seq, then ~ condition)", "batch_removed")
write.csv(rbind(r1$pcb, r2$pcb), file.path(out, "table09g_PC_batch_correlation.csv"), row.names = FALSE)
write.csv(data.frame(model = c(r1$label, r2$label), DE_padj05 = c(r1$n_de, r2$n_de), DE_padj05_lfc1 = c(r1$n_strong, r2$n_strong)),
          file.path(out, "table09h_DE_counts.csv"), row.names = FALSE)

# figures ---------------------------------------------------------------------------------------------------
pal <- c(ExpA = "#2C7FB8", ExpB = "#D95F0E")
th <- theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
pca_plot <- function(r, x, y) {
  ggplot(r$sc, aes(.data[[x]], .data[[y]], colour = batch)) + geom_point(size = 4) +
    geom_text(aes(label = sub("Exp[AB]_", "", sample)), colour = "black", size = 2.8, vjust = -1.2) +
    scale_colour_manual(values = pal, labels = c(ExpA = "ExpA = condition A", ExpB = "ExpB = condition B"), name = NULL) +
    labs(x = sprintf("%s (%.0f%%)", x, r$ve[as.integer(sub("PC", "", x))]), y = sprintf("%s (%.0f%%)", y, r$ve[as.integer(sub("PC", "", y))])) +
    scale_y_continuous(expand = expansion(mult = 0.15)) + th
}
row <- function(r) (pca_plot(r, "PC1", "PC2") + labs(title = r$label, subtitle = sprintf("%d DE genes (padj < 0.05)", r$n_de))) + pca_plot(r, "PC2", "PC3")
pA <- (row(r1) / row(r2)) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(out, "figure09a_PCA.png"), pA, width = 10, height = 9, dpi = 300, bg = "white")

pcb <- rbind(r1$pcb, r2$pcb)
pB <- ggplot(pcb, aes(PC, abs(r_with_batch), fill = model)) + geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%% var.\np = %.2g", pct_variance, p)), position = position_dodge(0.8), vjust = -0.3, size = 2.8) +
  scale_fill_manual(values = c("#636363", "#74C476"), name = NULL) + scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "|Pearson r| between PC score and batch", title = "How much of each PC is batch?") + th
ggsave(file.path(out, "figure09b_PC_batch_correlation.png"), pB, width = 7.5, height = 5, dpi = 300, bg = "white")

bars <- function(r) {
  t <- r$top20; t$lab <- ifelse(is.na(t$symbol) | t$symbol == t$gene, t$gene, paste0(t$gene, " (", t$symbol, ")"))
  ggplot(t, aes(reorder(lab, contrib_PC1), contrib_PC1, fill = loading_PC1 > 0)) + geom_col() + coord_flip() +
    scale_fill_manual(values = c(`TRUE` = "#D95F0E", `FALSE` = "#2C7FB8"), labels = c(`TRUE` = "positive loading", `FALSE` = "negative loading"), name = NULL) +
    labs(x = NULL, y = "% of PC1 variance explained by the gene", title = r$label,
         subtitle = sprintf("top 20 genes together: %.1f%% of PC1; positive PC1 side: %s", sum(t$contrib_PC1),
                           paste(sub("Exp[AB]_", "", r$sc$sample[r$sc$PC1 > 0]), collapse = ", "))) + th
}
pC <- bars(r1) + bars(r2)
ggsave(file.path(out, "figure09c_top20_PC1_contributions.png"), pC, width = 13, height = 6.5, dpi = 300, bg = "white")

gop <- function(g, title) {
  g <- head(g[g$FDR < 0.05, ], 12)
  if (!nrow(g)) return(ggplot() + annotate("text", 0, 0, label = "no GO term with FDR < 0.05") + theme_void() + labs(title = title))
  g$name <- substr(g$name, 1, 55)
  ggplot(g, aes(fold, reorder(name, fold), size = n_hit, colour = -log10(FDR))) + geom_point() +
    scale_colour_gradient(low = "#9ECAE1", high = "#08519C") + labs(x = "fold enrichment", y = NULL, title = title, size = "genes") + th +
    theme(legend.position = "right")
}
pD <- (gop(r1$go_de, "GO: DE genes, |log2FC| > 1 (batch ignored)") / gop(r1$go_pc, "GO: top 20 PC1 genes (batch ignored)") /
       gop(r2$go_pc, "GO: top 20 PC1 genes (batch removed)"))
ggsave(file.path(out, "figure09d_GO.png"), pD, width = 10, height = 12, dpi = 300, bg = "white")
write.csv(info[info$gene %in% rownames(counts), ], file.path(app, "gene_info.csv"), row.names = FALSE)
write.csv(sel[, c("sample", "run", "paper_replicate", "batch", "condition", "TIN")], file.path(app, "samples.csv"), row.names = FALSE)
