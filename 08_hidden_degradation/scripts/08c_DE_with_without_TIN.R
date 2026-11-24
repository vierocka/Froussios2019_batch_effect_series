# False DE in ExpA 3-vs-4 splits, without and with TIN as a covariate. Run from the repository root
# after 08b_degradation_metrics.R (reads table08a).
# All 7 ExpA replicates are used in every split, so replicate 6 (lowest TIN, 23.7% rRNA) is always in one
# group; with 3 vs. 4 samples it weighs more in the smaller group.
# EXCLUDE_REPS (e.g. "6") leaves replicates out (run 08b with the same setting first); with 6 replicates the
# splits are 3 vs. 3 (10 unique splits) and the output goes to results/without_rep<N>/.
suppressMessages({ library(DESeq2); library(ggplot2) })
work <- Sys.getenv("PIPELINE_DIR", file.path(getwd(), "pipeline_data"))
cores <- as.integer(Sys.getenv("THREADS", "4"))
excl <- as.integer(strsplit(Sys.getenv("EXCLUDE_REPS", ""), ",")[[1]])
out <- file.path("08_hidden_degradation/results", if (length(excl)) paste0("without_rep", paste(excl, collapse = "_")) else "")

met <- read.csv(file.path(out, "table08a_degradation_metrics.csv"))
read_count <- function(r) { d <- read.delim(file.path(work, "STAR_both/fC", paste0(r, ".count")), comment.char = "#")
  setNames(d[[ncol(d)]], d$Geneid) }
counts <- sapply(met$run, read_count)
colnames(counts) <- paste0("r", met$replicate)
counts <- counts[rowSums(counts) > 0, ]
tin <- setNames(met$TIN_median, colnames(counts))

# DE between the two groups; with_tin adds the centred sample TIN as a covariate (~ TIN + group)
n_de <- function(g1, g2, with_tin) {
  ids <- c(g1, g2)
  cd <- data.frame(group = factor(rep(c("g1", "g2"), c(length(g1), length(g2)))), TIN = tin[ids] - mean(tin[ids]),
                   row.names = ids)
  dds <- DESeqDataSetFromMatrix(counts[, ids], cd, if (with_tin) ~ TIN + group else ~ group)
  res <- results(DESeq(dds, quiet = TRUE), name = "group_g2_vs_g1")
  sum(res$padj < 0.05, na.rm = TRUE)
}

# all splits into floor(n/2) vs. the rest (7 replicates: 35 3-vs-4 splits; 6: 10 unique 3-vs-3 splits)
reps <- colnames(counts); k <- floor(length(reps) / 2)
g1s <- combn(reps, k, simplify = FALSE)
if (2 * k == length(reps)) g1s <- Filter(function(a) reps[1] %in% a, g1s)      # drop mirror images
splits <- lapply(g1s, function(a) list(g1 = a, g2 = setdiff(reps, a)))
design_lab <- sprintf("%d-vs-%d", k, length(reps) - k)
res <- do.call(rbind, parallel::mclapply(seq_along(splits), function(k) {
  s <- splits[[k]]
  data.frame(split = k, group1 = paste(sub("r", "", s$g1), collapse = ","), group2 = paste(sub("r", "", s$g2), collapse = ","),
             rep6_in = if (!"r6" %in% reps) "not included" else if ("r6" %in% s$g1) "group of 3" else "group of 4",
             TIN_group1 = mean(tin[s$g1]), TIN_group2 = mean(tin[s$g2]),
             abs_dTIN = abs(mean(tin[s$g1]) - mean(tin[s$g2])),
             DE_without_TIN = n_de(s$g1, s$g2, FALSE), DE_with_TIN = n_de(s$g1, s$g2, TRUE))
}, mc.cores = cores))
write.csv(res, file.path(out, "table08c_DE_all_splits.csv"), row.names = FALSE)

# 5 random splits for the figures
set.seed(8)
pick <- sort(sample(nrow(res), 5))
five <- res[pick, ]
five$label <- sprintf("%s vs. %s", five$group1, five$group2)
write.csv(five, file.path(out, "table08d_DE_5_random_splits.csv"), row.names = FALSE)
print(five[, c("label", "rep6_in", "abs_dTIN", "DE_without_TIN", "DE_with_TIN")], row.names = FALSE)

cat("\nall splits, median (mean) DE genes:\n")
print(aggregate(cbind(DE_without_TIN, DE_with_TIN) ~ rep6_in, res, function(x) sprintf("%g (%.0f)", median(x), mean(x))))
cat(sprintf("Spearman rho, DE vs. |dTIN|: without TIN %.2f, with TIN %.2f\n",
            cor(res$DE_without_TIN, res$abs_dTIN, method = "spearman"), cor(res$DE_with_TIN, res$abs_dTIN, method = "spearman")))

# figures: one panel per split, samples placed by group at their TIN, replicate 6 in orange
pts <- do.call(rbind, lapply(seq_len(nrow(five)), function(i) {
  g1 <- strsplit(five$group1[i], ",")[[1]]; g2 <- strsplit(five$group2[i], ",")[[1]]
  data.frame(panel = i, replicate = c(g1, g2), group = rep(c("group 1", "group 2"), c(length(g1), length(g2))),
             TIN = tin[paste0("r", c(g1, g2))])
}))
plot_splits <- function(col, title, file) {
  lab <- setNames(sprintf("%s\n%d DE genes", five$label, five[[col]]), seq_len(nrow(five)))
  p <- ggplot(pts, aes(group, TIN)) +
    geom_point(aes(colour = replicate == "6"), size = 7, position = position_jitter(width = 0.18, height = 0, seed = 1)) +
    geom_text(aes(label = replicate), colour = "white", size = 3.2, fontface = "bold",
              position = position_jitter(width = 0.18, height = 0, seed = 1)) +
    facet_wrap(~panel, nrow = 1, labeller = as_labeller(lab)) +
    scale_colour_manual(values = c(`TRUE` = "#D95F0E", `FALSE` = "#2C7FB8"), labels = c(`TRUE` = "replicate 6", `FALSE` = "other replicates"), name = NULL) +
    labs(x = NULL, y = "sample TIN (median)", title = title,
         subtitle = sprintf("5 random %s splits of %d ExpA replicates (same genotype and batch)%s; DE = padj < 0.05",
                           design_lab, length(reps), if (length(excl)) paste0(", replicate ", paste(excl, collapse = ", "), " left out") else "")) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom", strip.text = element_text(face = "bold"))
  ggsave(file.path(out, file), p, width = 13, height = 5, dpi = 300, bg = "white")
}
plot_splits("DE_without_TIN", "DE without TIN in the model (~ group)", "figure08d_splits_without_TIN.png")
plot_splits("DE_with_TIN", "DE with TIN as a covariate (~ TIN + group)", "figure08e_splits_with_TIN.png")

# all splits: DE vs. TIN imbalance between the groups
long <- rbind(data.frame(res[, c("abs_dTIN", "rep6_in")], model = "~ group", DE = res$DE_without_TIN),
              data.frame(res[, c("abs_dTIN", "rep6_in")], model = "~ TIN + group", DE = res$DE_with_TIN))
pF <- ggplot(long, aes(abs_dTIN, DE + 1, colour = rep6_in)) +
  geom_point(size = 2.5, alpha = 0.8) + facet_wrap(~model) + scale_y_log10() +
  scale_colour_manual(values = c("group of 3" = "#D95F0E", "group of 4" = "#2C7FB8", "not included" = "#2C7FB8"), name = "replicate 6 in the") +
  labs(x = "|difference in mean TIN between the groups|", y = "DE genes + 1 (log scale)",
       title = sprintf("All %d %s splits of ExpA%s", nrow(res), design_lab, if (length(excl)) paste0(" without replicate ", paste(excl, collapse = ", ")) else "")) +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"), legend.position = if (length(excl)) "none" else "bottom")
ggsave(file.path(out, "figure08f_DE_vs_TIN_imbalance.png"), pF, width = 9, height = 5, dpi = 300, bg = "white")
