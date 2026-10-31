"""
Cut Discover 465 (human sympathetic ganglia multiome, Yang & Gereau; CC BY 4.0) down to a local
test set: two aggregated libraries (barcode suffixes -1 and -2, as samples sg1 and sg2), up to
245 cells per author cell type in each (about 2,000 cells per sample), and write

  <out>/fragments/sg<N>/atac_fragments.tsv.gz   the cells' fragments (one pass over the 17.7 GB file)
  <out>/cells/sg<N>.h5ad                         their RNA counts (Gene Expression only), obs_names sg<N>:<barcode>
  <out>/authors_celltypes.csv                    sample_id, obs_id, cell_type (the authors' Celltypes)

  python scripts/cut_sympathetic.py <data/sympathetic-multiome/primary> <data/work/sympathetic/cut>

Barcodes keep their aggregation suffix. meta.data.csv's per-cell counts match these barcodes
(nCount_RNA r = 0.9997), but its "sample id" column is filled in blocks the size of each suffix's
cell count, in file order, and its Celltypes show no marker-gene specificity (C1QA, CD163, MPZ,
PECAM1, DCN ...): the annotation columns do not line up with cell_name in the published file. So
samples are named by suffix, and authors_celltypes.csv is only kept for reference; atac-regulome's
demo uses clusters of the ATAC data instead (atac-regulome/scripts/demo_labels.py).
"""
from __future__ import annotations

import csv
import gzip
import subprocess
import sys
from pathlib import Path

import numpy as np

SUFFIXES = ("1", "2")
PER_TYPE = 245
SEED = 0


def main(src: Path, out: Path) -> None:
    rng = np.random.default_rng(SEED)
    meta = list(csv.DictReader(open(src / "meta.data.csv", newline="")))
    chosen: dict[str, tuple[str, str]] = {}          # barcode -> (sample_id, cell_type)
    for sfx in SUFFIXES:
        rows = [r for r in meta if r["cell_name"].rsplit("-", 1)[1] == sfx]
        for t in sorted({r["Celltypes"] for r in rows}):
            cells = sorted(r["cell_name"] for r in rows if r["Celltypes"] == t)
            for b in rng.choice(cells, size=min(PER_TYPE, len(cells)), replace=False):
                chosen[str(b)] = (f"sg{sfx}", t)
    (out / "cells").mkdir(parents=True, exist_ok=True)
    with open(out / "authors_celltypes.csv", "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["sample_id", "obs_id", "cell_type"])
        for b, (sid, t) in sorted(chosen.items(), key=lambda x: (x[1][0], x[0])):
            w.writerow([sid, f"{sid}:{b}", t])
    print(f"{len(chosen)} cells chosen", flush=True)
    write_rna(src / "filtered_feature_bc_matrix", out, chosen)
    write_fragments(src / "atac_fragments.tsv.gz", out, chosen)


def write_rna(mtx: Path, out: Path, chosen: dict[str, tuple[str, str]]) -> None:
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp

    barcodes = [l.strip() for l in gzip.open(mtx / "barcodes.tsv.gz", "rt")]
    feats = [l.rstrip("\n").split("\t") for l in gzip.open(mtx / "features.tsv.gz", "rt")]
    col = {i + 1: b for i, b in enumerate(barcodes) if b in chosen}
    gex = [i + 1 for i, f in enumerate(feats) if f[2] == "Gene Expression"]
    max_gex = max(gex)
    # one awk pass over the 585 M entries keeps only the chosen columns' Gene Expression rows
    keep = out / "cells" / "columns.txt"
    keep.write_text("\n".join(map(str, col)) + "\n")
    awk = f'NR==FNR{{k[$1]=1;next}} /^%/{{next}} !h{{h=1;next}} ($2 in k) && $1<={max_gex}'
    unzip = subprocess.Popen(["gzip", "-dc", str(mtx / "matrix.mtx.gz")], stdout=subprocess.PIPE)
    proc = subprocess.Popen(["awk", awk, str(keep), "-"], stdin=unzip.stdout, stdout=subprocess.PIPE, text=True)
    unzip.stdout.close()
    r, c, v = [], [], []
    for line in proc.stdout:
        a, b, x = line.split()
        r.append(int(a) - 1)
        c.append(int(b))
        v.append(int(x))
    proc.wait()
    unzip.wait()
    keep.unlink()
    cols = sorted(col)
    ci = {k: i for i, k in enumerate(cols)}
    M = sp.csr_matrix((v, ([ci[x] for x in c], r)), shape=(len(cols), len(feats)))[:, [g - 1 for g in gex]]
    var = pd.DataFrame({"gene_ids": [feats[g - 1][0] for g in gex], "feature_types": "Gene Expression"},
                       index=[feats[g - 1][1] for g in gex])
    for sid in sorted({s for s, _ in chosen.values()}):
        idx = [i for i, k in enumerate(cols) if chosen[col[k]][0] == sid]
        a = ad.AnnData(X=M[idx].astype(np.float32), obs=pd.DataFrame(index=[f"{sid}:{col[cols[i]]}" for i in idx]),
                       var=var.copy())
        a.var_names_make_unique()
        a.obs["sample_id"] = pd.Categorical([sid] * a.n_obs)
        a.uns["pennsieve"] = {"organism": "human", "chemistry": "multiome", "source": "Discover 465 (cut)"}
        a.write_h5ad(out / "cells" / f"{sid}.h5ad", compression="gzip")
        print(f"{sid}: RNA {a.n_obs} cells x {a.n_vars} genes", flush=True)


def write_fragments(frag: Path, out: Path, chosen: dict[str, tuple[str, str]]) -> None:
    """One awk pass: header lines to every sample, a chosen barcode's fragments to its sample; then gzip."""
    plain = {s: out / "fragments" / s / "atac_fragments.tsv" for s in sorted({s for s, _ in chosen.values()})}
    for p in plain.values():
        p.parent.mkdir(parents=True, exist_ok=True)
    bc = out / "fragments" / "barcodes.tsv"
    bc.write_text("".join(f"{b}\t{plain[s]}\n" for b, (s, _) in chosen.items()))
    heads = " ".join(f'print > "{p}";' for p in plain.values())
    awk = f'NR==FNR{{f[$1]=$2;next}} /^#/{{{heads} next}} ($4 in f){{print > f[$4]}}'
    unzip = subprocess.Popen(["gzip", "-dc", str(frag)], stdout=subprocess.PIPE)
    subprocess.run(["awk", "-F", "\t", awk, str(bc), "-"], stdin=unzip.stdout, check=True)
    unzip.stdout.close()
    unzip.wait()
    bc.unlink()
    for p in plain.values():
        subprocess.run(["gzip", "-f", str(p)], check=True)
    print("fragments written", flush=True)


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
