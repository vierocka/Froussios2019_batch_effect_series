"""GO enrichment for Arabidopsis gene lists (imported by topic scripts; run from the repository root).

TAIR GAF (NOT qualifiers removed) propagated to all ancestors via is_a and part_of (go-basic.obo);
hypergeometric test per term with >= min_size annotated universe genes; BH within each ontology.
The GO files are downloaded once into PIPELINE_DIR/GO_annotation.
"""
import os
import urllib.request
from functools import lru_cache
from pathlib import Path

import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


def go_files():
    d = Path(os.environ.get("PIPELINE_DIR", "pipeline_data")) / "GO_annotation"
    d.mkdir(parents=True, exist_ok=True)
    gaf, obo = d / "tair.gaf.gz", d / "go-basic.obo"
    for f, url in ((gaf, "http://current.geneontology.org/annotations/tair.gaf.gz"),
                   (obo, "http://current.geneontology.org/ontology/go-basic.obo")):
        if not f.exists():
            urllib.request.urlretrieve(url, f)
    return gaf, obo


def _read_gaf():
    return pd.read_csv(go_files()[0], sep="\t", comment="!", header=None, dtype=str, quoting=3, compression="gzip")


def gene_info():
    """gene -> symbol and description (for labels and the Shiny apps)."""
    x = _read_gaf()[[1, 2, 9, 10]].drop_duplicates(subset=1)
    x.columns = ["gene", "symbol", "description", "synonyms"]
    return x.reset_index(drop=True)


@lru_cache(maxsize=1)
def go_annotation():
    gaf = _read_gaf()
    gaf = gaf[~gaf[3].fillna("").str.contains("NOT")][[1, 4]].drop_duplicates()
    gaf.columns = ["gene", "go"]
    terms, cur = {}, None
    with open(go_files()[1]) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("["):
                cur = {"parents": []} if line == "[Term]" else None
                continue
            if cur is None or ": " not in line:
                continue
            k, val = line.split(": ", 1)
            if k == "id":
                cur["id"] = val
                terms[val] = cur
            elif k in ("name", "namespace"):
                cur[k] = val
            elif k == "is_a":
                cur["parents"].append(val.split()[0])
            elif k == "relationship" and val.startswith("part_of "):
                cur["parents"].append(val.split()[1])
            elif k == "is_obsolete" and val == "true":
                terms.pop(cur["id"], None)
                cur = None

    @lru_cache(maxsize=None)
    def ancestors(t):
        a = {t}
        for p in terms[t]["parents"]:
            if p in terms:
                a |= ancestors(p)
        return frozenset(a)

    gaf = gaf[gaf.go.isin(terms.keys())]
    g2go = pd.DataFrame([(g, a) for g, t in zip(gaf.gene, gaf.go) for a in ancestors(t)],
                        columns=["gene", "go"]).drop_duplicates()
    return g2go, terms


def go_enrich(genes, universe, min_size=5):
    g2go, terms = go_annotation()
    g2go = g2go[g2go.gene.isin(set(universe))]
    u = set(g2go.gene)
    s = set(genes) & u
    size = g2go.go.value_counts()
    hits = g2go[g2go.gene.isin(s)].go.value_counts()
    keep = size.index[size >= min_size]
    N, n = len(u), len(s)
    r = pd.DataFrame({"go": keep, "n_term": size[keep].values, "n_hit": hits.reindex(keep, fill_value=0).values})
    r["expected"] = r.n_term * n / N
    r["fold"] = r.n_hit / r.expected
    r["p"] = hypergeom.sf(r.n_hit - 1, N, r.n_term, n)
    r["name"] = [terms[t]["name"] for t in r.go]
    r["ontology"] = [terms[t]["namespace"] for t in r.go]
    r["FDR"] = r.groupby("ontology").p.transform(lambda p: multipletests(p, method="fdr_bh")[1])
    hg = g2go[g2go.gene.isin(s)].groupby("go").gene.apply(lambda g: ";".join(sorted(g)))
    r["genes"] = r.go.map(hg).fillna("")
    r = r.sort_values(["FDR", "p"], kind="stable")[
        ["ontology", "go", "name", "n_term", "n_hit", "expected", "fold", "p", "FDR", "genes"]].reset_index(drop=True)
    r.attrs["n_annotated"] = {"list": n, "universe": N}
    return r
