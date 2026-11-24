# Degradation metrics of the 7 ExpA replicates: TIN per sample, TIN vs. rRNA %, TIN vs. singleton reads.
# Run from the repository root after 08a_tin_cluster.sh (pipeline_data/TIN/) and 00c_trim_star_featurecounts.sh
# (STAR logs in pipeline_data/STAR_both, STAR_Fonly, STAR_Ronly).
# EXCLUDE_REPS (e.g. "6") leaves replicates out; the output then goes to results/without_rep<N>/.
suppressMessages({ library(ggplot2); library(patchwork) })
work <- Sys.getenv("PIPELINE_DIR", file.path(getwd(), "pipeline_data"))
excl <- as.integer(strsplit(Sys.getenv("EXCLUDE_REPS", ""), ",")[[1]])
out <- file.path("08_hidden_degradation/results", if (length(excl)) paste0("without_rep", paste(excl, collapse = "_")) else "")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

map <- read.delim("data/replicate_map.tsv"); map <- map[map$paper_experiment == "ExpA" & !map$paper_replicate %in% excl, ]
pct <- read.csv("data/paper_TableS2C_rRNA.csv")
map$pct_rRNA <- pct$pct_rRNA[match(map$paper_replicate, pct$replicate)]

# TIN: median / mean / SD over transcripts (tin.py summary) and the per-transcript values (TIN > 0)
tin_sum <- do.call(rbind, lapply(map$run, function(r) {
  s <- read.delim(file.path(work, "TIN", paste0(r, "_Aligned.sortedByCoord.out.summary.txt")))
  data.frame(run = r, TIN_mean = s[1, 2], TIN_median = s[1, 3], TIN_sd = s[1, 4])
}))
tin_tx <- do.call(rbind, lapply(map$run, function(r) {
  x <- read.delim(file.path(work, "TIN", paste0(r, "_Aligned.sortedByCoord.out.tin.xls")))
  data.frame(run = r, TIN = x$TIN[x$TIN > 0])
}))

# singletons: reads whose mate was dropped by trimming, mapped as forward-only or reverse-only (STAR logs)
star_input <- function(dir, r) {
  l <- readLines(file.path(work, dir, paste0(r, "_Log.final.out")))
  as.numeric(sub(".*\\|\\s*", "", grep("Number of input reads", l, value = TRUE)))
}
reads <- data.frame(run = map$run,
                    pairs = sapply(map$run, star_input, dir = "STAR_both"),
                    fwd_only = sapply(map$run, star_input, dir = "STAR_Fonly"),
                    rev_only = sapply(map$run, star_input, dir = "STAR_Ronly"))
reads$pct_singletons <- 100 * (reads$fwd_only + reads$rev_only) / (reads$pairs + reads$fwd_only + reads$rev_only)
reads$pct_fwd_only <- 100 * reads$fwd_only / (reads$pairs + reads$fwd_only + reads$rev_only)

d <- merge(merge(map[, c("run", "paper_replicate", "pct_rRNA")], tin_sum), reads)
d <- d[order(d$paper_replicate), ]
names(d)[names(d) == "paper_replicate"] <- "replicate"
write.csv(d, file.path(out, "table08a_degradation_metrics.csv"), row.names = FALSE)
print(d, row.names = FALSE)

cors <- function(x, y) {
  s <- cor.test(x, y, method = "spearman", exact = TRUE); p <- cor.test(x, y, method = "pearson")
  data.frame(spearman_rho = unname(s$estimate), spearman_p = s$p.value, pearson_r = unname(p$estimate), pearson_p = p$p.value)
}
ct <- rbind(cbind(comparison = "TIN median vs. rRNA %", cors(d$TIN_median, d$pct_rRNA)),
            if (6 %in% d$replicate) cbind(comparison = "TIN median vs. rRNA % without replicate 6",
                                          cors(d$TIN_median[d$replicate != 6], d$pct_rRNA[d$replicate != 6])),
            cbind(comparison = "TIN median vs. singleton %", cors(d$TIN_median, d$pct_singletons)),
            cbind(comparison = "rRNA % vs. singleton %", cors(d$pct_rRNA, d$pct_singletons)))
write.csv(ct, file.path(out, "table08b_correlations.csv"), row.names = FALSE)
print(ct, row.names = FALSE, digits = 3)

col6 <- function(r) ifelse(r == 6, "replicate 6", "other ExpA replicates")
pal <- c("replicate 6" = "#D95F0E", "other ExpA replicates" = "#2C7FB8")
th <- theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

tin_tx$replicate <- factor(map$paper_replicate[match(tin_tx$run, map$run)], levels = sort(map$paper_replicate))
pA <- ggplot(tin_tx, aes(replicate, TIN, fill = col6(as.integer(as.character(replicate))))) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.2, width = 0.6) +
  geom_point(data = d, aes(factor(replicate), TIN_median), inherit.aes = FALSE, shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = pal, name = NULL) +
  labs(x = "ExpA replicate", y = "transcript integrity number (TIN)",
       title = "TIN per sample", subtitle = "boxes: all transcripts with TIN > 0; diamond: median (tin.py summary)") + th
ggsave(file.path(out, "figure08a_TIN_per_sample.png"), pA, width = 7, height = 5.5, dpi = 300, bg = "white")

scatter <- function(x, xlab, title, ct_row) {
  ggplot(d, aes(.data[[x]], TIN_median, colour = col6(replicate))) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey60", linetype = "dashed", linewidth = 0.6) +
    geom_point(size = 3.5) +
    geom_text(aes(label = replicate), colour = "white", size = 2.6, fontface = "bold") +
    scale_colour_manual(values = pal, name = NULL) +
    labs(x = xlab, y = "TIN (median over transcripts)", title = title,
         subtitle = sprintf("Spearman rho = %.2f (p = %.2g); Pearson r = %.2f (p = %.2g); n = %d",
                            ct_row$spearman_rho, ct_row$spearman_p, ct_row$pearson_r, ct_row$pearson_p, nrow(d))) + th
}
pB <- scatter("pct_rRNA", "rRNA % (paper's Table S2C)", "TIN vs. rRNA", ct[1, ])
pC <- scatter("pct_singletons", "singleton reads, % of mapped fragments (STAR logs)", "TIN vs. singletons",
              ct[ct$comparison == "TIN median vs. singleton %", ])
ggsave(file.path(out, "figure08b_TIN_vs_rRNA.png"), pB, width = 6.5, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(out, "figure08c_TIN_vs_singletons.png"), pC, width = 6.5, height = 5.5, dpi = 300, bg = "white")
