# Volcano plot and GO enrichment of split 6's DE genes, combined with panels A and B of
# 06a_split_structure.R into one 4-panel figure. Run from the repository root after 06a.
suppressMessages({ library(DESeq2); library(ggplot2); library(patchwork) })
work <- Sys.getenv("PIPELINE_DIR", file.path(getwd(), "pipeline_data"))
out <- "06_split_structure/results"
go_dir <- file.path(work, "GO_annotation")

# GO annotation: TAIR GAF and go-basic.obo from the GO Consortium (downloaded once if missing)
dir.create(go_dir, showWarnings = FALSE, recursive = TRUE)
gaf_file <- file.path(go_dir, "tair.gaf.gz"); obo_file <- file.path(go_dir, "go-basic.obo")
if (!file.exists(gaf_file)) download.file("http://current.geneontology.org/annotations/tair.gaf.gz", gaf_file, mode = "wb")
if (!file.exists(obo_file)) download.file("http://current.geneontology.org/ontology/go-basic.obo", obo_file, mode = "wb")

# coding-gene counts and split 6, exactly as in 06a_split_structure.R
map <- read.delim("data/replicate_map.tsv"); map <- map[map$paper_experiment != "ExpC", ]
pct <- read.csv("data/paper_TableS2C_rRNA.csv")
map$pct <- pct$pct_rRNA[match(map$paper_replicate, pct$replicate)]
rd <- function(r) read.delim(file.path(work, "STAR_both/fC", paste0(r, ".count")), comment.char = "#")
coding <- sapply(map$run, function(r) { d <- rd(r); setNames(d[[ncol(d)]], d$Geneid) })
coding <- coding[rowSums(coding) > 0, ]
colnames(coding) <- paste0("r", map$paper_replicate)
pr <- setNames(map$pct, colnames(coding))
g1 <- c("r1", "r2", "r3"); g2 <- c("r4", "r5", "r6")
cd <- data.frame(group = factor(rep(c("g1", "g2"), each = 3)), row.names = c(g1, g2))
dds <- DESeq(DESeqDataSetFromMatrix(coding[, c(g1, g2)], cd, ~group), quiet = TRUE)
res <- results(dds)
tested <- rownames(res)[!is.na(res$padj)]
sig <- rownames(res)[which(res$padj < 0.05)]
# check against the saved list of 06a_split_structure.R
saved <- read.csv(file.path(out, "table06a_split6_DE_genes.csv"))$gene
stopifnot(setequal(sig, saved))
cat(sprintf("split 6: %d DE genes (padj < 0.05) out of %d tested; matches table06a\n",
            length(sig), length(tested)))

# GO: direct annotations (NOT qualifiers removed), propagated to all ancestors (is_a, part_of)
gaf <- read.delim(gaf_file, comment.char = "!", header = FALSE, quote = "", colClasses = "character")
gaf <- gaf[!grepl("NOT", gaf$V4), c("V2", "V5")]
names(gaf) <- c("gene", "go")
obo <- readLines(obo_file)
starts <- grep("^\\[Term\\]$", obo); ends <- c(starts[-1] - 1, length(obo))
terms <- lapply(seq_along(starts), function(i) {
  b <- obo[starts[i]:ends[i]]
  if (any(b == "is_obsolete: true")) return(NULL)
  id <- sub("^id: ", "", grep("^id: ", b, value = TRUE)[1])
  par <- c(sub("^is_a: (GO:[0-9]+).*", "\\1", grep("^is_a: ", b, value = TRUE)),
           sub("^relationship: part_of (GO:[0-9]+).*", "\\1", grep("^relationship: part_of ", b, value = TRUE)))
  data.frame(id = id, name = sub("^name: ", "", grep("^name: ", b, value = TRUE)[1]),
             ns = sub("^namespace: ", "", grep("^namespace: ", b, value = TRUE)[1]),
             parents = paste(par, collapse = ";"))
})
terms <- do.call(rbind, terms)
parents <- setNames(strsplit(terms$parents, ";"), terms$id)
anc_cache <- new.env()
ancestors <- function(t) {
  if (!is.null(anc_cache[[t]])) return(anc_cache[[t]])
  p <- parents[[t]]; p <- p[p %in% names(parents)]
  a <- unique(c(t, unlist(lapply(p, ancestors))))
  assign(t, a, envir = anc_cache); a
}
gaf <- unique(gaf[gaf$go %in% terms$id & gaf$gene %in% tested, ])
anc <- lapply(unique(gaf$go), ancestors); names(anc) <- unique(gaf$go)
g2go <- unique(data.frame(gene = rep(gaf$gene, lengths(anc[gaf$go])), go = unlist(anc[gaf$go])))
universe <- unique(g2go$gene); sig_go <- intersect(sig, universe)
cat(sprintf("GO-annotated: %d of %d tested genes, %d of %d DE genes\n",
            length(universe), length(tested), length(sig_go), length(sig)))

# hypergeometric test per term (terms with at least 5 annotated tested genes, no upper limit so that
# broad terms such as GO:0003700 DNA-binding transcription factor activity are tested), BH within each ontology
size <- table(g2go$go); hits <- table(g2go$go[g2go$gene %in% sig_go])
keep <- names(size)[size >= 5]
N <- length(universe); n <- length(sig_go)
goea <- data.frame(go = keep, n_term = as.integer(size[keep]),
                   n_de = as.integer(ifelse(is.na(hits[keep]), 0, hits[keep])))
goea$expected <- goea$n_term * n / N
goea$fold <- goea$n_de / goea$expected
goea$p <- phyper(goea$n_de - 1, goea$n_term, N - goea$n_term, n, lower.tail = FALSE)
goea$name <- terms$name[match(goea$go, terms$id)]; goea$ontology <- terms$ns[match(goea$go, terms$id)]
goea$FDR <- ave(goea$p, goea$ontology, FUN = function(p) p.adjust(p, "BH"))
de_genes <- split(g2go$gene[g2go$gene %in% sig_go], g2go$go[g2go$gene %in% sig_go])
goea$DE_genes <- vapply(goea$go, function(t) paste(sort(de_genes[[t]]), collapse = ";"), "")
goea <- goea[order(goea$FDR, goea$p), c("ontology", "go", "name", "n_term", "n_de", "expected", "fold", "p", "FDR", "DE_genes")]
write.csv(goea, file.path(out, "table06c_GO_enrichment_split6.csv"), row.names = FALSE)
cat(sprintf("GO terms tested: %d; FDR < 0.05: %d\n", nrow(goea), sum(goea$FDR < 0.05)))
print(head(goea[goea$FDR < 0.05, c("ontology", "go", "name", "n_de", "n_term", "fold", "FDR")], 20), row.names = FALSE)
tf_terms <- goea[grepl("transcription", goea$name), ]
cat("\ntranscription-related terms (smallest FDR):\n")
print(head(tf_terms[, c("ontology", "go", "name", "n_de", "n_term", "fold", "FDR")], 8), row.names = FALSE)
tf_genes <- unique(g2go$gene[g2go$go == "GO:0003700"])   # DNA-binding transcription factor activity

# Figure
pal <- c("replicates 1, 2, 3" = "#2C7FB8", "replicates 4, 5, 6" = "#D95F0E")
# A: per-gene medians per group (as in 06a_split_structure.R)
nc6 <- counts(dds, normalized = TRUE)[tested, ]
med <- data.frame(gene = rownames(nc6), g1 = apply(nc6[, g1], 1, median), g2 = apply(nc6[, g2], 1, median))
up <- sig[res[sig, "log2FoldChange"] > 0]; dn <- sig[res[sig, "log2FoldChange"] < 0]
dA <- rbind(transform(med, set = sprintf("all %s tested genes", format(nrow(med), big.mark = ","))),
            transform(med[med$gene %in% up, ], set = sprintf("%d DE genes higher in replicates 4, 5, 6", length(up))),
            transform(med[med$gene %in% dn, ], set = sprintf("%d DE genes lower in replicates 4, 5, 6", length(dn))))
dA <- rbind(data.frame(set = dA$set, group = "replicates 1, 2, 3", median = dA$g1),
            data.frame(set = dA$set, group = "replicates 4, 5, 6", median = dA$g2))
dA$set <- factor(dA$set, levels = unique(dA$set))
pA <- ggplot(dA, aes(log2(median + 1), colour = group)) + geom_density(linewidth = 1) +
  facet_wrap(~set, ncol = 1, scales = "free_y") + scale_colour_manual(values = pal, name = NULL) +
  coord_cartesian(xlim = c(2, 18)) +
  labs(x = "per-gene median of normalized counts in the group (log2 + 1)", y = "density",
       title = "A. Gene medians per group") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# B: volcano, DESeq2 log2FC (unshrunken, as used for the DE call) vs. -log10 padj
dB <- data.frame(gene = tested, lfc = res[tested, "log2FoldChange"], padj = res[tested, "padj"])
dB$class <- ifelse(dB$padj >= 0.05, "not DE",
                   ifelse(dB$lfc > 0, "DE, higher in replicates 4, 5, 6", "DE, lower in replicates 4, 5, 6"))
dB$class <- factor(dB$class, levels = c("not DE", "DE, higher in replicates 4, 5, 6", "DE, lower in replicates 4, 5, 6"))
dB$tf <- dB$gene %in% tf_genes & dB$padj < 0.05
lab <- data.frame(gene = c("AT4G16780", "AT4G17460", "AT2G44910"), name = c("ATHB-2/HAT4", "HAT1", "ATHB4"))
lab <- merge(lab, dB)
pB <- ggplot(dB, aes(lfc, -log10(padj))) +
  geom_point(aes(colour = class), size = 0.8, alpha = 0.6) +
  geom_point(data = dB[dB$tf, ], shape = 21, size = 1.9, stroke = 0.5, colour = "black", fill = NA) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  ggrepel::geom_text_repel(data = lab, aes(label = name), size = 3.3, min.segment.length = 0,
                           nudge_x = 0.6, colour = "grey15") +
  scale_colour_manual(values = c("not DE" = "grey75", "DE, higher in replicates 4, 5, 6" = "#D95F0E",
                                 "DE, lower in replicates 4, 5, 6" = "#2C7FB8"), name = NULL) +
  labs(x = "log2 fold change (replicates 4, 5, 6 vs. 1, 2, 3)", y = "-log10 padj",
       title = "B. Volcano plot of split 6",
       subtitle = sprintf("dashed line: padj = 0.05; black ring: DE gene annotated as a transcription factor (GO:0003700, %d genes)",
                          sum(dB$tf))) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 9)) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1), nrow = 2))

# C: rRNA per replicate (panel B of 06a_split_structure.R)
dC <- data.frame(rep = c(g1, g2), pct = pr[c(g1, g2)], group = rep(names(pal), each = 3))
bars <- aggregate(pct ~ group, dC, median)
pC <- ggplot(dC, aes(group, pct, colour = group)) +
  geom_col(data = bars, aes(fill = group), colour = NA, alpha = 0.25, width = 0.6) +
  geom_point(size = 3.5) +
  ggrepel::geom_text_repel(aes(label = sub("r", "rep ", rep)), size = 3.5, direction = "y", nudge_x = 0.35,
                           colour = "grey20", min.segment.length = 0) +
  scale_colour_manual(values = pal, guide = "none") + scale_fill_manual(values = pal, guide = "none") +
  labs(x = NULL, y = "% rRNA (paper's Table S2C)", title = "C. rRNA per replicate",
       subtitle = "bar = group median") +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))

# D: GO enrichment, 12 terms with the smallest FDR; GO:0003700 (DNA-binding transcription factor activity) in orange
dD <- head(goea[goea$FDR < 0.05, ], 12)
if (nrow(dD) == 0) dD <- head(goea, 12)
dD$label <- sprintf("%s (%s, %d/%d)", dD$name, c(biological_process = "BP", molecular_function = "MF",
                                                 cellular_component = "CC")[dD$ontology], dD$n_de, dD$n_term)
dD$label <- factor(dD$label, levels = rev(dD$label))
dD$tf <- dD$go == "GO:0003700"
pD <- ggplot(dD, aes(-log10(FDR), label, fill = tf)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c(`TRUE` = "#D95F0E", `FALSE` = "grey60"), guide = "none") +
  labs(x = "-log10 FDR (hypergeometric, BH within ontology)", y = NULL,
       title = sprintf("D. GO enrichment of the %d DE genes (%d GO-annotated)", length(sig), length(sig_go)),
       subtitle = sprintf("background: %s GO-annotated tested genes; brackets: DE genes / genes in term",
                          format(length(universe), big.mark = ","))) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 9),
        axis.text.y = element_text(size = 9))

p <- (pA | pB) / (pC | pD) + plot_layout(widths = c(1, 1.3), heights = c(1.15, 1)) + plot_annotation(
  title = sprintf("Same genotype, same batch (ExpA), no designed difference: replicates 1, 2, 3 vs. 4, 5, 6 -> %d DE genes (padj < 0.05)", length(sig)),
  theme = theme(plot.title = element_text(size = 12, face = "bold")))
ggsave(file.path(out, "figure06b_split_structure_4panels.png"), p, width = 16, height = 12, dpi = 300, bg = "white")
cat("written: figure06b, table06c\n")
