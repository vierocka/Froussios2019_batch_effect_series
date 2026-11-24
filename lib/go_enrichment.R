# GO enrichment for Arabidopsis gene lists (source()d by topic scripts; run from the repository root).
# TAIR GAF (NOT qualifiers removed) propagated to all ancestors via is_a and part_of (go-basic.obo);
# hypergeometric test per term with >= min_size annotated universe genes; BH within each ontology.
# The GO files are downloaded once into PIPELINE_DIR/GO_annotation.
go_files <- function() {
  d <- file.path(Sys.getenv("PIPELINE_DIR", file.path(getwd(), "pipeline_data")), "GO_annotation")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  f <- c(gaf = file.path(d, "tair.gaf.gz"), obo = file.path(d, "go-basic.obo"))
  if (!file.exists(f["gaf"])) download.file("http://current.geneontology.org/annotations/tair.gaf.gz", f["gaf"], mode = "wb")
  if (!file.exists(f["obo"])) download.file("http://current.geneontology.org/ontology/go-basic.obo", f["obo"], mode = "wb")
  f
}

# gene -> symbol and description (for labels and the Shiny apps)
gene_info <- function() {
  gaf <- read.delim(go_files()["gaf"], comment.char = "!", header = FALSE, quote = "", colClasses = "character")
  x <- unique(gaf[, c("V2", "V3", "V10", "V11")]); x <- x[!duplicated(x$V2), ]
  data.frame(gene = x$V2, symbol = x$V3, description = x$V10, synonyms = x$V11)
}

.go_cache <- new.env()
go_annotation <- function() {
  if (!is.null(.go_cache$g2go)) return(.go_cache)
  f <- go_files()
  gaf <- read.delim(f["gaf"], comment.char = "!", header = FALSE, quote = "", colClasses = "character")
  gaf <- unique(gaf[!grepl("NOT", gaf$V4), c("V2", "V5")]); names(gaf) <- c("gene", "go")
  obo <- readLines(f["obo"])
  starts <- grep("^\\[Term\\]$", obo); ends <- c(starts[-1] - 1, length(obo))
  terms <- do.call(rbind, lapply(seq_along(starts), function(i) {
    b <- obo[starts[i]:ends[i]]
    if (any(b == "is_obsolete: true")) return(NULL)
    par <- c(sub("^is_a: (GO:[0-9]+).*", "\\1", grep("^is_a: ", b, value = TRUE)),
             sub("^relationship: part_of (GO:[0-9]+).*", "\\1", grep("^relationship: part_of ", b, value = TRUE)))
    data.frame(id = sub("^id: ", "", grep("^id: ", b, value = TRUE)[1]),
               name = sub("^name: ", "", grep("^name: ", b, value = TRUE)[1]),
               ns = sub("^namespace: ", "", grep("^namespace: ", b, value = TRUE)[1]), parents = paste(par, collapse = ";"))
  }))
  parents <- setNames(strsplit(terms$parents, ";"), terms$id)
  anc_cache <- new.env()
  ancestors <- function(t) {
    if (!is.null(anc_cache[[t]])) return(anc_cache[[t]])
    p <- parents[[t]]; p <- p[p %in% names(parents)]
    a <- unique(c(t, unlist(lapply(p, ancestors)))); assign(t, a, envir = anc_cache); a
  }
  gaf <- gaf[gaf$go %in% terms$id, ]
  ug <- unique(gaf$go); anc <- setNames(lapply(ug, ancestors), ug)
  .go_cache$g2go <- unique(data.frame(gene = rep(gaf$gene, lengths(anc[gaf$go])), go = unlist(anc[gaf$go])))
  .go_cache$terms <- terms
  .go_cache
}

go_enrich <- function(genes, universe, min_size = 5) {
  a <- go_annotation(); terms <- a$terms
  g2go <- a$g2go[a$g2go$gene %in% universe, ]
  u <- unique(g2go$gene); s <- intersect(genes, u)
  size <- table(g2go$go); hits <- table(g2go$go[g2go$gene %in% s])
  keep <- names(size)[size >= min_size]
  N <- length(u); n <- length(s)
  r <- data.frame(go = keep, n_term = as.integer(size[keep]), n_hit = as.integer(ifelse(is.na(hits[keep]), 0, hits[keep])))
  r$expected <- r$n_term * n / N; r$fold <- r$n_hit / r$expected
  r$p <- phyper(r$n_hit - 1, r$n_term, N - r$n_term, n, lower.tail = FALSE)
  r$name <- terms$name[match(r$go, terms$id)]; r$ontology <- terms$ns[match(r$go, terms$id)]
  r$FDR <- ave(r$p, r$ontology, FUN = function(p) p.adjust(p, "BH"))
  hg <- split(g2go$gene[g2go$gene %in% s], g2go$go[g2go$gene %in% s])
  r$genes <- vapply(r$go, function(t) paste(sort(hg[[t]]), collapse = ";"), "")
  r <- r[order(r$FDR, r$p), c("ontology", "go", "name", "n_term", "n_hit", "expected", "fold", "p", "FDR", "genes")]
  attr(r, "n_annotated") <- c(list = n, universe = N)
  r
}
