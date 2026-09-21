# Post figure: rRNA fraction per sample (paper Table S2C) vs. mean rlog correlation with the other samples.
# Run from the repository root. Needs outputs of 02_align_count and 03_rRNA_qc.
suppressMessages({ library(DESeq2); library(ggplot2); library(patchwork) })

work <- Sys.getenv("WORK_DIR", file.path(getwd(), "work"))

# 14 ExpA + ExpB samples; replicate numbers follow the paper, not ENA sample titles
map  <- read.delim("data/replicate_map.tsv")
map  <- map[map$paper_experiment %in% c("ExpA", "ExpB"), ]
pct  <- read.csv("data/paper_TableS2C_rRNA.csv")
samples <- merge(map, pct[, c("replicate", "pct_rRNA")], by.x = "paper_replicate", by.y = "replicate")
samples <- samples[order(samples$pct_rRNA), ]
samples$run <- factor(samples$run, levels = samples$run)
samples$flag <- ifelse(samples$paper_replicate == 11, "rep 11\n(excluded)",
                ifelse(samples$paper_replicate == 6,  "rep 6\n(kept)", ""))

# gene counts (featureCounts -t gene has no rRNA rows), plus the 4 rRNA genes from the rRNA count files
rrna_genes <- c("AT2G01010", "AT2G01020", "AT3G41768", "AT3G41979")
read_count <- function(f) { d <- read.delim(f, comment.char = "#"); setNames(d[[ncol(d)]], d$Geneid) }
gene_mat <- sapply(as.character(samples$run), function(r) read_count(file.path(work, "STAR_both/fC", paste0(r, ".count"))))
rrna_mat <- sapply(as.character(samples$run), function(r) read_count(file.path(work, "rRNA_qc/counts", paste0(r, "_rRNA.count")))[rrna_genes])
stopifnot(!any(rrna_genes %in% rownames(gene_mat)))
counts <- rbind(gene_mat, rrna_mat)
storage.mode(counts) <- "integer"

dds <- DESeqDataSetFromMatrix(counts, colData = samples, design = ~1)
dds <- dds[rowSums(counts(dds)) > 0, ]
rl  <- assay(rlog(dds, blind = TRUE))

mean_cor <- function(genes) { cm <- cor(rl[genes, ]); diag(cm) <- NA; rowMeans(cm, na.rm = TRUE) }
set.seed(42)
lab <- c(sprintf("All %s genes (incl. 4 rRNA)", format(nrow(rl), big.mark = ",")),
         sprintf("%s genes (4 rRNA removed)", format(nrow(rl) - 4, big.mark = ",")),
         "Random 5,000 genes (avg. of 50 draws)")
cor_df <- data.frame(
  run = rep(samples$run, 3),
  condition = factor(rep(lab, each = nrow(samples)), levels = lab),
  mean_cor = c(mean_cor(rownames(rl))[as.character(samples$run)],
               mean_cor(setdiff(rownames(rl), rrna_genes))[as.character(samples$run)],
               rowMeans(replicate(50, mean_cor(sample(rownames(rl), 5000))))[as.character(samples$run)]))

pal <- c(ExpA = "#2C7FB8", ExpB = "#D95F0E")
panelA <- ggplot(samples, aes(run, pct_rRNA, fill = paper_experiment)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = flag), vjust = -0.6, size = 3, fontface = "bold", lineheight = 0.9) +
  scale_fill_manual(values = pal, name = "Experiment") +
  labs(x = NULL, y = "% of reads mapping to rRNA\n(paper's Table S2C)",
       title = sprintf("A. rRNA content varies ~%d-fold across 14 'same-condition' replicates",
                       round(max(samples$pct_rRNA) / min(samples$pct_rRNA)))) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold", size = 12),
        legend.position = "top", plot.margin = margin(5, 95, 5, 5)) +
  expand_limits(y = max(samples$pct_rRNA) * 1.25)

panelB <- ggplot(cor_df, aes(run, mean_cor, color = condition, group = condition)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  scale_color_manual(values = c("#333333", "#E7298A", "#66A61E"), name = NULL) +
  labs(x = "sample (ordered by rRNA%, same order as panel A)",
       y = "mean rlog correlation\nwith the other 13 samples",
       title = "B. ...but genome-wide correlation barely moves either way") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold", size = 12),
        legend.position = "top", legend.text = element_text(size = 10)) +
  coord_cartesian(ylim = c(min(cor_df$mean_cor) - 0.002, 1))

combined <- panelA / panelB + plot_annotation(
  title = "Froussios et al. 2019 (Bioinformatics), Arabidopsis Col-0, ExpA+ExpB (n=14)",
  theme = theme(plot.title = element_text(size = 11, face = "italic", hjust = 0)))
ggsave("scripts/post2_rRNA_correlation/rRNA_correlation_blindspot.png", combined, width = 8.5, height = 9.5, dpi = 300, bg = "white")
