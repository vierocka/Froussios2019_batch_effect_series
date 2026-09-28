# Post 4: what does an omitted feature do to normalization? Run from the repository root.
suppressMessages(library(ggplot2))

work <- Sys.getenv("WORK_DIR", file.path(getwd(), "work"))

# coding-gene counts (featureCounts -t gene has no rRNA rows) for the 14 ExpA/ExpB samples
map <- read.delim("data/replicate_map.tsv"); map <- map[map$paper_experiment != "ExpC", ]
pct <- read.csv("data/paper_TableS2C_rRNA.csv")
map$p <- pct$pct_rRNA[match(map$paper_replicate, pct$replicate)] / 100
read_count <- function(f) { d <- read.delim(f, comment.char = "#"); setNames(d[[ncol(d)]], d$Geneid) }
coding <- sapply(map$run, function(r) read_count(file.path(work, "STAR_both/fC", paste0(r, ".count"))))
coding <- coding[rowSums(coding) > 0, ]

# SYNTHETIC rRNA total, not a recovered read count: under a simplified two-component model
# (counted coding genes plus rRNA), how much rRNA would have to exist for the counted coding
# total to represent the paper's published fraction p of coding+rRNA? rrna = coding_total*p/(1-p)
# solves R/(C+R) = p by construction. Our own current-annotation rRNA recount is known to badly
# undercount (see the wider investigation's notes), so this uses the paper's published fraction
# (Table S2C) instead of a real per-locus count.
s <- map$run[map$paper_replicate == 11]; p11 <- map$p[map$paper_replicate == 11]
rrna_s <- round(sum(coding[, s]) * p11 / (1 - p11))

# One library (replicate 11, 31% rRNA per the paper). log2 share of the library per coding gene,
# computed two ways: with the synthetic rRNA total in the denominator, and with rRNA ignored
# (coding genes only). Raw coding-gene counts never change -- only their computed share of the
# total does, by a constant log2(1/(1-p)) for this sample.
x <- coding[, s]; x <- x[x > 0]
dA <- data.frame(share = c(x / (sum(x) + rrna_s), x / sum(x)),
                 total = rep(c(sprintf("rRNA in library total (rRNA = %.0f%%)", 100 * p11), "rRNA ignored (coding genes only)"), each = length(x)))
pA <- ggplot(dA, aes(log2(share), colour = total)) + geom_density(linewidth = 1) +
  scale_colour_manual(values = c("#D95F0E", "#2C7FB8"), name = NULL) +
  labs(x = "log2 share of library, per coding gene", y = "density",
       title = sprintf("One library (rep 11, %.0f%% rRNA, synthetic):\nall coding-gene library shares shift by %.2f log2", 100 * p11, log2(1 / (1 - p11))),
       subtitle = sprintf("A gene worth 1%% of the complete library becomes %.2f%% of the coding-only total\n(its share rises %.0f%%; it did not gain a single read)",
                           1 / (1 - p11), 100 * (1 / (1 - p11) - 1))) +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 12)) +
  guides(colour = guide_legend(nrow = 2))

ggsave("scripts/post4_rRNA_ignored/rRNA_ignored.png", pA, width = 7.5, height = 6.5, dpi = 300, bg = "white")
