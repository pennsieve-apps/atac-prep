"""
One sample: fragments -> SnapATAC2 AnnData with QC in obs and a 500 bp tile matrix.

1. Which barcodes are cells: the RNA side's barcodes when `cells` is linked (they passed
   sc-prep/sc-qc, so both modalities describe the same nuclei), else Cell Ranger's cell calls
   (per_barcode_metrics.csv is_cell / singlecell.csv is__cell_barcode), else every barcode with
   at least MIN_FRAGMENTS fragments.
2. snap.pp.import_fragments (backed, straight into the output file; unsorted input is sorted
   on disk), snap.metrics.tsse against the GENCODE annotation.
3. ATAC QC: n_fragment in [MIN_FRAGMENTS, MAX_FRAGMENTS], tsse >= MIN_TSSE. Failing cells go.
4. Fragment-size distribution of the kept cells, snap.pp.add_tile_matrix (500 bp).
5. Without `cells`: Scrublet doublets (snap.pp.select_features + scrublet + filter_doublets);
   with `cells` the RNA side already removed doublets.
6. obs_names <sample_id>:<barcode>, obs sample_id + sample metadata (categorical), uns["pennsieve"].
"""
from __future__ import annotations

import csv
import logging
from pathlib import Path

import numpy as np
import polars as pl

from processor.genome import Genome

log = logging.getLogger("atac-prep")
BIN_SIZE = 500
DOUBLET_PROBABILITY = 0.5
MIN_CELLS_FOR_SCRUBLET = 200


def cell_calls(path: Path) -> list[str]:
    """Barcodes Cell Ranger called as cells."""
    column = next(c for s, c in (("per_barcode_metrics.csv", "is_cell"), ("singlecell.csv", "is__cell_barcode"))
                  if path.name.endswith(s))
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        if column not in (reader.fieldnames or []) or "barcode" not in reader.fieldnames:
            raise ValueError(f"{path.name} has no barcode / {column} column")
        return [r["barcode"] for r in reader if r[column].strip() in ("1", "True", "true")]


def rna_barcodes(path: Path, sample_id: str) -> tuple[list[str], str | None]:
    """Barcodes of the RNA side (obs_names without the '<sample_id>:' prefix) and its organism."""
    import anndata as ad

    a = ad.read_h5ad(path, backed="r")
    try:
        pre = f"{sample_id}:"
        names = [n[len(pre):] if n.startswith(pre) else n for n in a.obs_names.astype(str)]
        organism = (a.uns.get("pennsieve") or {}).get("organism")
    finally:
        a.file.close()
    return names, (str(organism) if organism else None)


def _categorical(values: list[str]) -> pl.Series:
    return pl.Series(values, dtype=pl.Categorical)


def process(sample, genome: Genome, whitelist: list[str] | None, calls_from: str, params: dict, out_h5ad: Path,
            threads: int, prep_version: str, source: str) -> dict:
    """Write out_h5ad; returns the QC rows, figure data and numbers for the summary."""
    import snapatac2 as snap

    sid = sample.sample_id
    adata = snap.pp.import_fragments(
        sample.fragments, chrom_sizes=genome.chrom_sizes, file=out_h5ad, sorted_by_barcode=False,
        min_num_fragments=0 if whitelist is not None else params["MIN_FRAGMENTS"], whitelist=whitelist,
        tempdir=Path("/tmp"), n_jobs=threads)
    try:
        n_imported = adata.n_obs
        if n_imported == 0:
            raise ValueError("no cell has fragments on the genome's chromosomes ("
                             + ", ".join(list(genome.chrom_sizes)[:3]) + ", …); check the genome and the barcodes")
        snap.metrics.tsse(adata, genome.annotation, n_jobs=threads)
        obs = adata.obs[:].to_pandas()
        barcodes = list(adata.obs_names)
        nfrag, tsse = obs["n_fragment"].to_numpy(), obs["tsse"].to_numpy()
        reasons = [[] for _ in barcodes]
        for i in np.flatnonzero(nfrag < params["MIN_FRAGMENTS"]):
            reasons[i].append("low_fragments")
        for i in np.flatnonzero(nfrag > params["MAX_FRAGMENTS"]):
            reasons[i].append("high_fragments")
        for i in np.flatnonzero(~(tsse >= params["MIN_TSSE"])):
            reasons[i].append("low_tsse")
        keep = np.array([not r for r in reasons])
        figure = {"n_fragment": nfrag.copy(), "tsse": tsse.copy(), "keep": keep.copy(),
                  "tss_profile": np.asarray(adata.uns["TSS_profile"], dtype=float),
                  "library_tsse": float(adata.uns["library_tsse"])}
        if not keep.any():
            raise ValueError(f"no cell passed ATAC QC ({n_imported} evaluated)")
        adata.subset(obs_indices=np.flatnonzero(keep))
        snap.metrics.frag_size_distr(adata, n_jobs=threads)
        figure["frag_size"] = np.asarray(adata.uns["frag_size_distr"], dtype=float)
        snap.pp.add_tile_matrix(adata, bin_size=BIN_SIZE, n_jobs=threads)

        scrublet = None
        if calls_from != "rna":
            if adata.n_obs >= MIN_CELLS_FOR_SCRUBLET:
                snap.pp.select_features(adata, n_jobs=threads, verbose=False)
                snap.pp.scrublet(adata, random_state=params["SEED"], n_jobs=threads, verbose=False)
                prob = adata.obs["doublet_probability"].to_numpy()
                kept_bcs = list(adata.obs_names)
                dbl = prob > DOUBLET_PROBABILITY
                scrublet = {"ran": True, "n_doublets": int(dbl.sum())}
                if dbl.all():
                    raise ValueError("Scrublet called every cell a doublet")
                if dbl.any():
                    adata.subset(obs_indices=np.flatnonzero(~dbl))
                    idx = {b: i for i, b in enumerate(barcodes)}
                    for b, d in zip(kept_bcs, dbl):
                        if d:
                            reasons[idx[b]].append("doublet")
            else:
                scrublet = {"ran": False, "reason": f"fewer than {MIN_CELLS_FOR_SCRUBLET} cells"}

        kept = list(adata.obs_names)
        adata.obs_names = [f"{sid}:{b}" for b in kept]
        n = adata.n_obs
        adata.obs["sample_id"] = _categorical([sid] * n)
        for col, val in sample.metadata.items():
            if col in ("sample_id", "n_fragment", "frac_dup", "frac_mito", "tsse"):
                continue
            adata.obs[col] = _categorical([val] * n)
        adata.uns["pennsieve"] = {
            "source": source, "genome": genome.name, "genome_from": genome.how,
            "organism": "human" if genome.name == "GRCh38" else "mouse", "cells_from": calls_from, "tile_bin_size": BIN_SIZE, "prep_version": prep_version, "sample_id": sid,
            "qc": {"min_fragments": params["MIN_FRAGMENTS"], "max_fragments": params["MAX_FRAGMENTS"],
                   "min_tsse": params["MIN_TSSE"], "scrublet": bool(scrublet and scrublet["ran"])},
        }
        kept_obs = adata.obs[:].to_pandas()
    finally:
        adata.close()

    rows = [{"sample_id": sid, "obs_id": f"{sid}:{b}", "n_fragment": int(nfrag[i]),
             "frac_dup": float(obs["frac_dup"].iat[i]), "frac_mito": float(obs["frac_mito"].iat[i]),
             "tsse": float(tsse[i]), "passed": not reasons[i], "reason": ";".join(reasons[i])}
            for i, b in enumerate(barcodes)]
    removed: dict[str, int] = {}
    for r in reasons:
        for x in r:
            removed[x] = removed.get(x, 0) + 1
    stats = {
        "n_evaluated": n_imported, "n_kept": len(kept), "removed": removed, "cells_from": calls_from,
        "n_whitelist": len(whitelist) if whitelist is not None else None,
        "genome": genome.name, "genome_from": genome.how, "library_tsse": figure["library_tsse"],
        "median_n_fragment": float(np.median(kept_obs["n_fragment"])), "median_tsse": float(np.median(kept_obs["tsse"])),
        "median_frac_dup": float(np.median(kept_obs["frac_dup"])), "median_frac_mito": float(np.median(kept_obs["frac_mito"])),
        "scrublet": scrublet,
    }
    return {"rows": rows, "figure": figure, "stats": stats}
