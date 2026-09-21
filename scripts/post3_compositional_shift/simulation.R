# Teaching example: 10 extreme genes (like rRNA) barely move the correlation of two log2-scale expression vectors.
# Run from the repository root.
set.seed(1)
n <- 10000; n_changed <- 10; n_loops <- 1000; rho <- 0.9
lo <- log2(0.01); hi <- 12                      # log2 scale, range of typical rlog / log2 values

# two vectors with correlation ~0.9 that both span [lo, hi] (Gaussian copula)
z <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, rho, rho, 1), 2))
x <- lo + (hi - lo) * pnorm(z[, 1])
y <- lo + (hi - lo) * pnorm(z[, 2])
r0 <- cor(x, y)

# linear-scale share of the total taken by a set of genes
share <- function(v, idx) sum(2^v[idx]) / sum(2^v)

r_new <- numeric(n_loops); share_new <- numeric(n_loops)
for (i in seq_len(n_loops)) {
  idx <- sample(n, n_changed)
  y2 <- y
  y2[idx] <- runif(n_changed, hi - 1, hi)        # 10 random genes jump to the top of the range
  r_new[i] <- cor(x, y2)
  share_new[i] <- share(y2, idx)
}

# same test, but the 10 genes are set to hold 30% of the linear-scale total (like the 31% rRNA library)
r_30 <- replicate(n_loops, {
  idx <- sample(n, n_changed); y2 <- y
  y2[idx] <- log2(0.3 / 0.7 * sum(2^y[-idx]) / n_changed) + rnorm(n_changed, 0, 0.1)
  cor(x, y2)
})

cat(sprintf("baseline r = %.4f\n", r0))
cat(sprintf("after changing %d of %d genes, %d loops: r range %.4f - %.4f, mean shift %.4f\n",
            n_changed, n, n_loops, min(r_new), max(r_new), mean(r_new) - r0))
cat(sprintf("30%% share scenario: r range %.4f - %.4f\n", min(r_30), max(r_30)))
cat(sprintf("the %d genes hold %.1f%% of the linear-scale total (mean)\n", n_changed, 100 * mean(share_new)))

library(ggplot2); library(patchwork)
idx <- sample(n, n_changed); y2 <- y; y2[idx] <- runif(n_changed, hi - 1, hi)
d <- data.frame(x = x, y = y2, changed = seq_len(n) %in% idx)
pA <- ggplot(d[!d$changed, ], aes(x, y)) +
  geom_point(alpha = 0.15, size = 0.5, colour = "grey40") +
  geom_point(data = d[d$changed, ], colour = "#2C7FB8", size = 2.5) +
  labs(x = "sample 1 (log2 scale)", y = "sample 2 (log2 scale)",
       title = sprintf("A. r = %.3f; 10 genes (blue) pushed to the top of the range", r0)) +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold", size = 11))
pB <- ggplot(data.frame(r = c(r_new, r_30), scenario = rep(c("~1% of total signal", "30% of total signal"), each = n_loops)),
             aes(r, fill = scenario)) +
  geom_histogram(binwidth = 0.0005, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("#2C7FB8", "#D95F0E"), name = "10 changed genes hold", guide = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = r0, linetype = 2) +
  coord_cartesian(xlim = c(0.86, 0.91)) +
  labs(x = "correlation after the change", y = "loops (of 1000)",
       title = sprintf("B. %d random draws of 10 genes; dashed = original r", n_loops)) +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold", size = 11), legend.position = "bottom")
ggsave("scripts/post3_compositional_shift/compositional_shift.png", pA + pB, width = 11, height = 4.5, dpi = 300, bg = "white")
