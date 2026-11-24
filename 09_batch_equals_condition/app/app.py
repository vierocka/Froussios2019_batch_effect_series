"""Shiny for Python app for topic 09 (batch = condition): 4 best-TIN ExpA vs. 4 best-TIN ExpB replicates.

Reads the tables written by 09a_batch_equals_condition.R (results/app_data/). Start from the repository root:
    shiny run 09_batch_equals_condition/app/app.py --launch-browser
Needs: pip install shiny faicons pandas matplotlib
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from faicons import icon_svg
from shiny import App, reactive, render, ui

DD = Path(__file__).resolve().parent.parent / "results" / "app_data"
MODELS = {"batch_ignored": "batch ignored (~ condition)", "batch_removed": "batch removed (ComBat-seq, then ~ condition)"}
PAL = {"ExpA": "#2C7FB8", "ExpB": "#D95F0E"}
STRING_INPUT = "https://string-db.org/cgi/input?input_page_active_form=multiple_identifiers"

samples = pd.read_csv(DD / "samples.csv")
S = samples["sample"].tolist()
info = pd.read_csv(DD / "gene_info.csv").drop_duplicates("gene").set_index("gene")
D = {m: {"nc": pd.read_csv(DD / f"normalized_counts_{m}.csv").set_index("gene"),
         "de": pd.read_csv(DD / f"DE_all_genes_{m}.csv").set_index("gene", drop=False),
         "sc": pd.read_csv(DD / f"PCA_scores_{m}.csv"),
         "ve": pd.read_csv(DD / f"PCA_variance_{m}.csv").set_index("PC").pct_variance,
         "ct": pd.read_csv(DD / f"PCA_contributions_{m}.csv").set_index("gene", drop=False)} for m in MODELS}
GENES = D["batch_ignored"]["nc"].index.tolist()


def search_label(g):
    if g not in info.index or pd.isna(info.symbol.get(g)):
        return g
    desc = info.description.get(g)
    return f"{g} | {info.symbol[g]} | {'' if pd.isna(desc) else str(desc)[:70]}"


CHOICES = {g: search_label(g) for g in GENES}


def tair(g):
    return f"https://www.arabidopsis.org/locus?name={g}"


def string_net(genes):
    return "https://string-db.org/cgi/network?species=3702&identifiers=" + "%0d".join(genes)


def lg(m):
    return np.log2(D[m]["nc"][S] + 1)


def model_pick(id_):
    return ui.input_radio_buttons(id_, None, MODELS, inline=True)


def top_card(k):
    return ui.card(ui.card_header(f"PC{k}"), ui.output_data_frame(f"top_PC{k}"),
                   ui.download_button(f"dl_PC{k}", "download"), ui.output_ui(f"str_PC{k}"))


app_ui = ui.page_navbar(
    ui.nav_panel(
        "Data structure",
        ui.navset_card_tab(
            ui.nav_panel(
                "Sample-wise",
                ui.p("Eight Col-0 wild-type replicates: the 4 replicates with the highest transcript integrity number "
                     "(TIN) from ExpA (labelled condition A) and from ExpB (labelled condition B). Batch and condition "
                     "are the same column, so every difference between A and B is a batch effect."),
                ui.layout_columns(ui.output_data_frame("samples_tab"), ui.output_plot("design", height="330px"),
                                  col_widths=(6, 6))),
            ui.nav_panel(
                "Data-wise",
                ui.output_ui("datasentence"), model_pick("model_data"),
                ui.h5("Replicate-replicate plots: each ExpA replicate vs. each ExpB replicate (log2 normalized counts + 1)"),
                ui.output_plot("rr", height="820px"),
                ui.h5("Distribution of normalized values per sample"),
                ui.output_plot("dist", height="430px"))),
        icon=icon_svg("table-cells")),
    ui.nav_panel(
        "Gene search",
        ui.layout_columns(
            ui.input_selectize("gene", None, choices=[], width="100%",
                               options={"placeholder": "Type an AGI code, symbol or word of the description "
                                                       "(e.g. CCA1, AT2G46830, light)"}),
            model_pick("model_gene"), col_widths=(8, 4)),
        ui.output_ui("geneinfo"),
        ui.layout_columns(ui.output_plot("geneplot", height="420px"),
                          ui.card(ui.card_header("DE result and PCA rank"), ui.output_table("genede"),
                                  ui.output_table("generank")), col_widths=(7, 5)),
        icon=icon_svg("magnifying-glass")),
    ui.nav_panel(
        "Analysis",
        model_pick("model_an"),
        ui.layout_columns(ui.output_plot("pca12", height="420px"), ui.output_plot("pca23", height="420px"),
                          col_widths=(6, 6)),
        ui.h5("The 20 genes contributing most to PC1, PC2 and PC3"),
        ui.p("Contribution = squared loading in % (the contributions of all genes to one PC sum to 100%). Paste a list into STRING (organism Arabidopsis thaliana): ",
             ui.a("STRING multiple-protein input", href=STRING_INPUT, target="_blank"), "."),
        ui.layout_columns(top_card(1), top_card(2), top_card(3), col_widths=(4, 4, 4)),
        ui.h5("All DE genes (padj < 0.05)"), ui.output_ui("desentence"), ui.output_data_frame("detable"),
        ui.download_button("dl_de", "download DE genes"),
        icon=icon_svg("chart-line")),
    title="Batch = condition: ExpA (4) vs. ExpB (4)",
)


def server(input, output, session):
    ui.update_selectize("gene", choices=CHOICES, selected="AT2G46830", server=True)

    @render.data_frame
    def samples_tab():
        return render.DataGrid(samples.round({"TIN": 2}))

    @render.plot
    def design():
        fig, ax = plt.subplots()
        for i, r in samples.iterrows():
            y = len(samples) - i
            for x, val in ((0, r.batch), (1, r.condition)):
                ax.add_patch(plt.Rectangle((x, y - 0.45), 0.95, 0.9, color=PAL[r.batch]))
                ax.text(x + 0.475, y, val, ha="center", va="center", color="white", fontweight="bold")
            ax.text(2.1, y, f"TIN {r.TIN:.1f}", va="center")
            ax.text(-0.1, y, r["sample"], ha="right", va="center")
        ax.set_xlim(-1.2, 2.9); ax.set_ylim(0.3, len(samples) + 0.9); ax.axis("off")
        ax.text(0.475, len(samples) + 0.7, "batch", ha="center"); ax.text(1.475, len(samples) + 0.7, "condition", ha="center")
        ax.set_title("Design: batch and condition are identical", loc="left", fontweight="bold")
        return fig

    @render.ui
    def datasentence():
        return ui.p(f"{len(GENES)} genes with at least 10 reads in at least 4 of the 8 samples are analysed "
                    "(normalized counts: DESeq2 median-of-ratios).")

    @render.plot
    def rr():
        x = lg(input.model_data())
        a = samples["sample"][samples.batch == "ExpA"].tolist(); b = samples["sample"][samples.batch == "ExpB"].tolist()
        fig, axes = plt.subplots(4, 4, figsize=(10, 10), sharex=True, sharey=True)
        hi = x.values.max()
        for i, bb in enumerate(b):
            for j, aa in enumerate(a):
                ax = axes[i, j]
                ax.scatter(x[aa], x[bb], s=0.3, alpha=0.15, color="#636363", rasterized=True)
                ax.plot([0, hi], [0, hi], "k--", lw=0.5)
                ax.text(0.5, hi, f"r = {np.corrcoef(x[aa], x[bb])[0, 1]:.4f}", va="top", fontsize=8)
                if i == 3:
                    ax.set_xlabel(aa)
                if j == 0:
                    ax.set_ylabel(bb)
        fig.tight_layout()
        return fig

    @render.plot
    def dist():
        x = lg(input.model_data())
        fig, axes = plt.subplots(2, 4, figsize=(12, 5), sharex=True, sharey=True)
        for ax, s, b in zip(axes.flat, S, samples.batch):
            ax.hist(x[s], bins=60, color=PAL[b]); ax.set_title(s, fontsize=9)
        fig.supxlabel("log2(normalized count + 1)"); fig.supylabel("genes")
        fig.tight_layout()
        return fig

    @render.ui
    def geneinfo():
        g = input.gene()
        if not g:
            return ui.p("Type a gene above.")
        sym = f" ({info.symbol[g]})" if g in info.index else ""
        desc = info.description[g] if g in info.index else ""
        return ui.TagList(ui.h4(g + sym), ui.p("" if pd.isna(desc) else desc),
                          ui.a("open in TAIR", href=tair(g), target="_blank", class_="btn btn-outline-primary btn-sm"))

    @render.plot
    def geneplot():
        g = input.gene()
        if not g:
            return None
        v = D[input.model_gene()]["nc"].loc[g, S]
        fig, ax = plt.subplots()
        ax.bar(S, v, color=[PAL[b] for b in samples.batch])
        ax.set_ylabel("normalized count"); ax.set_title(g, loc="left", fontweight="bold")
        plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
        fig.tight_layout()
        return fig

    @render.table
    def genede():
        g = input.gene()
        if not g:
            return None
        d = D[input.model_gene()]["de"].loc[g]
        return pd.DataFrame({"statistic": ["baseMean", "log2 fold change (B vs. A)", "p-value", "adjusted p-value (FDR)"],
                             "value": [f"{d.baseMean:.1f}", f"{d.log2FC:.3f}", f"{d.pvalue:.2e}", f"{d.padj:.2e}"]})

    @render.table
    def generank():
        g = input.gene()
        if not g:
            return None
        ct = D[input.model_gene()]["ct"]; d = ct.loc[g]
        return pd.DataFrame({"PC": ["PC1", "PC2", "PC3"],
                             "rank of contribution": [f"{int(d[f'rank_PC{k}'])} of {len(ct)}" for k in (1, 2, 3)],
                             "contribution (%)": [f"{d[f'contrib_PC{k}']:.4f}" for k in (1, 2, 3)]})

    def pca(x, y):
        r = D[input.model_an()]
        fig, ax = plt.subplots()
        for _, p in r["sc"].iterrows():
            ax.scatter(p[x], p[y], s=90, color=PAL[p.batch])
            ax.annotate(p["sample"], (p[x], p[y]), textcoords="offset points", xytext=(0, 8), ha="center", fontsize=8)
        ax.set_xlabel(f"{x} ({r['ve'][x]:.0f}%)"); ax.set_ylabel(f"{y} ({r['ve'][y]:.0f}%)")
        ax.margins(0.15)
        return fig

    @render.plot
    def pca12():
        return pca("PC1", "PC2")

    @render.plot
    def pca23():
        return pca("PC2", "PC3")

    def top(k):
        ct = D[input.model_an()]["ct"].sort_values(f"rank_PC{k}").head(20)
        return pd.DataFrame({"gene": ct.gene, "symbol": ct.symbol, "contribution": ct[f"contrib_PC{k}"].round(3)})

    def make_top(k):
        @render.data_frame
        def _tab():
            return render.DataGrid(top(k), height="560px")

        @render.download(filename=lambda: f"top20_PC{k}_{input.model_an()}.txt")
        def _dl():
            yield "\n".join(top(k).gene) + "\n"

        @render.ui
        def _str():
            return ui.a("open these 20 in STRING", href=string_net(top(k).gene), target="_blank")
        return _tab, _dl, _str

    for k in (1, 2, 3):
        t, dl, s = make_top(k)
        output(id=f"top_PC{k}")(t); output(id=f"dl_PC{k}")(dl); output(id=f"str_PC{k}")(s)

    @reactive.calc
    def sigde():
        d = D[input.model_an()]["de"]
        d = d[d.padj < 0.05].sort_values("padj").reset_index(drop=True)
        return d.round({"baseMean": 1, "log2FC": 3, "lfcSE": 3})

    @render.ui
    def desentence():
        return ui.p(f"{len(sigde())} DE genes (padj < 0.05) of {D[input.model_an()]['de'].padj.notna().sum()} tested.")

    @render.data_frame
    def detable():
        return render.DataGrid(sigde(), filters=True, height="500px")

    @render.download(filename=lambda: f"DE_genes_{input.model_an()}.csv")
    def dl_de():
        yield sigde().to_csv(index=False)


app = App(app_ui, server)
