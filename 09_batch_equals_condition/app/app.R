# Shiny app for topic 09 (batch = condition): 4 best-TIN ExpA vs. 4 best-TIN ExpB replicates.
# Reads the tables written by 09a_batch_equals_condition.R (results/app_data/). Start from the repository root:
#   Rscript -e 'shiny::runApp("09_batch_equals_condition/app", launch.browser = TRUE)'
suppressMessages({ library(shiny); library(bslib); library(DT); library(ggplot2) })

dd <- file.path("..", "results", "app_data")
if (!dir.exists(dd)) dd <- "09_batch_equals_condition/results/app_data"
rd <- function(f) read.csv(file.path(dd, f), check.names = FALSE)
samples <- rd("samples.csv")
info <- rd("gene_info.csv")
models <- c("batch ignored (~ condition)" = "batch_ignored", "batch removed (ComBat-seq, then ~ condition)" = "batch_removed")
D <- lapply(models, function(m) list(
  nc = rd(sprintf("normalized_counts_%s.csv", m)), de = rd(sprintf("DE_all_genes_%s.csv", m)),
  sc = rd(sprintf("PCA_scores_%s.csv", m)), ve = rd(sprintf("PCA_variance_%s.csv", m)),
  ct = rd(sprintf("PCA_contributions_%s.csv", m))))
names(D) <- models
genes <- D[[1]]$nc$gene
gi <- info[match(genes, info$gene), ]
search_lab <- ifelse(is.na(gi$symbol), genes,
                     sprintf("%s | %s | %s", genes, gi$symbol, ifelse(is.na(gi$description), "", substr(gi$description, 1, 70))))
choices <- setNames(genes, search_lab)
pal <- c(ExpA = "#2C7FB8", ExpB = "#D95F0E")
th <- theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
tair <- function(g) sprintf("https://www.arabidopsis.org/locus?name=%s", g)
string_net <- function(g) sprintf("https://string-db.org/cgi/network?species=3702&identifiers=%s", paste(g, collapse = "%0d"))
STRING_INPUT <- "https://string-db.org/cgi/input?input_page_active_form=multiple_identifiers"
model_pick <- function(id) radioButtons(id, NULL, choices = models, inline = TRUE)
lg <- function(nc) log2(as.matrix(nc[, samples$sample]) + 1)

ui <- page_navbar(
  title = "Batch = condition: ExpA (4) vs. ExpB (4)", theme = bs_theme(bootswatch = "flatly"),
  nav_panel(title = "Data structure", icon = icon("table-cells"),
    navset_card_tab(
      nav_panel("Sample-wise",
        p("Eight Col-0 wild-type replicates: the 4 replicates with the highest transcript integrity number (TIN) from ExpA",
          "(labelled condition A) and from ExpB (labelled condition B). Batch and condition are the same column, so",
          "every difference between A and B is a batch effect."),
        layout_columns(col_widths = c(6, 6), DTOutput("samples"), plotOutput("design", height = "330px"))),
      nav_panel("Data-wise",
        uiOutput("datasentence"), model_pick("model_data"),
        h5("Replicate-replicate plots: each ExpA replicate vs. each ExpB replicate (log2 normalized counts + 1)"),
        plotOutput("rr", height = "820px"),
        h5("Distribution of normalized values per sample"),
        plotOutput("dist", height = "430px")))),
  nav_panel(title = "Gene search", icon = icon("magnifying-glass"),
    layout_columns(col_widths = c(8, 4),
      selectizeInput("gene", NULL, choices = NULL, width = "100%",
                     options = list(placeholder = "Type an AGI code, symbol or word of the description (e.g. CCA1, AT2G46830, light)")),
      model_pick("model_gene")),
    uiOutput("geneinfo"),
    layout_columns(col_widths = c(7, 5), plotOutput("geneplot", height = "420px"),
                   card(card_header("DE result and PCA rank"), tableOutput("genede"), tableOutput("generank")))),
  nav_panel(title = "Analysis", icon = icon("chart-line"),
    model_pick("model_an"),
    layout_columns(col_widths = c(6, 6), plotOutput("pca12", height = "420px"), plotOutput("pca23", height = "420px")),
    h5("The 20 genes contributing most to PC1, PC2 and PC3"),
    p("Contribution = squared loading in % (the contributions of all genes to one PC sum to 100%).",
      "Paste a list into STRING (organism Arabidopsis thaliana): ", a("STRING multiple-protein input", href = STRING_INPUT, target = "_blank"), "."),
    layout_columns(col_widths = c(4, 4, 4),
      card(card_header("PC1"), DTOutput("top_PC1"), card_footer(downloadButton("dl_PC1", "download"), uiOutput("str_PC1", inline = TRUE))),
      card(card_header("PC2"), DTOutput("top_PC2"), card_footer(downloadButton("dl_PC2", "download"), uiOutput("str_PC2", inline = TRUE))),
      card(card_header("PC3"), DTOutput("top_PC3"), card_footer(downloadButton("dl_PC3", "download"), uiOutput("str_PC3", inline = TRUE)))),
    h5("All DE genes (padj < 0.05)"), uiOutput("desentence"), DTOutput("detable"), downloadButton("dl_de", "download DE genes")))

server <- function(input, output, session) {
  updateSelectizeInput(session, "gene", choices = choices, server = TRUE, selected = "AT2G46830")

  output$samples <- renderDT(datatable(samples, rownames = FALSE, options = list(dom = "t")) |> formatRound("TIN", 2))
  output$design <- renderPlot({
    s <- samples; s$sample <- factor(s$sample, rev(s$sample))
    m <- rbind(data.frame(s, col = "batch", val = s$batch), data.frame(s, col = "condition", val = s$condition))
    ggplot(m, aes(col, sample, fill = val)) + geom_tile(colour = "white", linewidth = 1) + geom_text(aes(label = val), colour = "white", fontface = "bold") +
      scale_fill_manual(values = c(ExpA = "#2C7FB8", A = "#2C7FB8", ExpB = "#D95F0E", B = "#D95F0E"), guide = "none") +
      geom_text(aes(x = 2.75, label = sprintf("TIN %.1f", TIN)), size = 3.5) + scale_x_discrete(expand = expansion(add = c(0.6, 0.9))) +
      labs(x = NULL, y = NULL, title = "Design: batch and condition are identical") + th
  })
  output$datasentence <- renderUI(p(sprintf("%d genes with at least 10 reads in at least 4 of the 8 samples are analysed (normalized counts: DESeq2 median-of-ratios).",
                                            length(genes))))
  output$rr <- renderPlot({
    x <- lg(D[[input$model_data]]$nc); a <- samples$sample[samples$batch == "ExpA"]; b <- samples$sample[samples$batch == "ExpB"]
    d <- do.call(rbind, lapply(a, function(i) do.call(rbind, lapply(b, function(j)
      data.frame(A = i, B = j, x = x[, i], y = x[, j], r = cor(x[, i], x[, j]))))))
    lab <- unique(d[, c("A", "B", "r")])
    ggplot(d, aes(x, y)) + geom_point(size = 0.2, alpha = 0.15, colour = "#636363") + geom_abline(linetype = "dashed") +
      geom_text(data = lab, aes(x = 0, y = max(d$y), label = sprintf("r = %.4f", r)), hjust = 0, vjust = 1, size = 3.5) +
      facet_grid(B ~ A) + coord_equal() + labs(x = "ExpA replicate", y = "ExpB replicate") + th
  })
  output$dist <- renderPlot({
    x <- lg(D[[input$model_data]]$nc)
    d <- data.frame(sample = rep(colnames(x), each = nrow(x)), v = as.vector(x))
    d$batch <- samples$batch[match(d$sample, samples$sample)]; d$sample <- factor(d$sample, samples$sample)
    ggplot(d, aes(v, fill = batch)) + geom_histogram(bins = 60) + facet_wrap(~sample, nrow = 2) +
      scale_fill_manual(values = pal) + labs(x = "log2(normalized count + 1)", y = "genes") + th
  })

  output$geneinfo <- renderUI({
    req(input$gene); g <- input$gene; i <- info[info$gene == g, ]
    tagList(h4(g, if (nrow(i)) paste0(" (", i$symbol, ")")), if (nrow(i)) p(i$description),
            a("open in TAIR", href = tair(g), target = "_blank", class = "btn btn-outline-primary btn-sm"))
  })
  output$geneplot <- renderPlot({
    req(input$gene); nc <- D[[input$model_gene]]$nc; v <- unlist(nc[nc$gene == input$gene, samples$sample])
    d <- data.frame(sample = factor(samples$sample, samples$sample), v = v, batch = samples$batch)
    ggplot(d, aes(sample, v, fill = batch)) + geom_col(width = 0.7) + scale_fill_manual(values = pal) +
      labs(x = NULL, y = "normalized count", title = input$gene) + th + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  })
  output$genede <- renderTable({
    req(input$gene); d <- D[[input$model_gene]]$de; d <- d[d$gene == input$gene, ]
    data.frame(statistic = c("baseMean", "log2 fold change (B vs. A)", "p-value", "adjusted p-value (FDR)"),
               value = c(sprintf("%.1f", d$baseMean), sprintf("%.3f", d$log2FC), formatC(d$pvalue, format = "e", digits = 2), formatC(d$padj, format = "e", digits = 2)))
  })
  output$generank <- renderTable({
    req(input$gene); ct <- D[[input$model_gene]]$ct; d <- ct[ct$gene == input$gene, ]
    data.frame(PC = paste0("PC", 1:3), `rank of contribution` = sprintf("%d of %d", unlist(d[, paste0("rank_PC", 1:3)]), nrow(ct)),
               `contribution (%)` = sprintf("%.4f", unlist(d[, paste0("contrib_PC", 1:3)])), check.names = FALSE)
  })

  pca <- function(x, y) renderPlot({
    r <- D[[input$model_an]]; s <- r$sc; ve <- setNames(r$ve$pct_variance, r$ve$PC)
    ggplot(s, aes(.data[[x]], .data[[y]], colour = batch)) + geom_point(size = 5) +
      geom_text(aes(label = sample), colour = "black", vjust = -1.3, size = 3.5) + scale_colour_manual(values = pal) +
      scale_y_continuous(expand = expansion(mult = 0.15)) +
      labs(x = sprintf("%s (%.0f%%)", x, ve[x]), y = sprintf("%s (%.0f%%)", y, ve[y])) + th
  })
  output$pca12 <- pca("PC1", "PC2"); output$pca23 <- pca("PC2", "PC3")
  top <- function(k) { ct <- D[[input$model_an]]$ct; ct <- ct[order(ct[[paste0("rank_PC", k)]]), ][1:20, ]
    data.frame(gene = ct$gene, symbol = ct$symbol, contribution = round(ct[[paste0("contrib_PC", k)]], 3)) }
  for (k in 1:3) local({
    kk <- k
    output[[paste0("top_PC", kk)]] <- renderDT(datatable(top(kk), rownames = FALSE, fillContainer = FALSE, options = list(dom = "t", pageLength = 20)))
    output[[paste0("dl_PC", kk)]] <- downloadHandler(function() sprintf("top20_PC%d_%s.txt", kk, input$model_an),
                                                    function(f) writeLines(top(kk)$gene, f))
    output[[paste0("str_PC", kk)]] <- renderUI(a("open these 20 in STRING", href = string_net(top(kk)$gene), target = "_blank"))
  })
  sigde <- reactive({ d <- D[[input$model_an]]$de; d <- d[which(d$padj < 0.05), ]; d[order(d$padj), ] })
  output$desentence <- renderUI(p(sprintf("%d DE genes (padj < 0.05) of %d tested.", nrow(sigde()), sum(!is.na(D[[input$model_an]]$de$padj)))))
  output$detable <- renderDT(datatable(sigde(), rownames = FALSE, filter = "top") |> formatSignif(c("baseMean", "log2FC", "lfcSE", "pvalue", "padj"), 3))
  output$dl_de <- downloadHandler(function() sprintf("DE_genes_%s.csv", input$model_an), function(f) write.csv(sigde(), f, row.names = FALSE))
}

shinyApp(ui, server)
