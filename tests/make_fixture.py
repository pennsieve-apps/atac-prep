"""
A tiny synthetic input in the shapes atac-prep reads, plus a tiny genome for it:

  python -m tests.make_fixture data/fixture [--layout linked|flat|both]

writes data/fixture/genomes/GRCh38/ (annotation.gff3.gz, chrom_sizes.tsv: two short chromosomes
and chrM with 60 genes; the processor reads it through GENOMES_DIR) and data/fixture/<layout>/:

  knee1/outs/          Cell Ranger ARC per-sample layout: atac_fragments.tsv.gz (with ARC's header
  knee2/outs/          lines, reference GRCh38), per_barcode_metrics.csv (is_cell), atac_peaks.bed and
                       filtered_feature_bc_matrix.h5 (the RNA half, not read by atac-prep)
  plain/S3_fragments.tsv.gz (+ .tbi)   a bare 10x-style fragment file: no header, no cell calls,
                       220 cells and 20 near-empty barcodes
  samples.csv          sample sheet: sample_id, path, condition
  cells/<id>.h5ad      (linked) the RNA side for knee1 and knee2 as sc-qc writes it: obs_names
                       <sample_id>:<barcode>, a subset of ARC's cells (RNA QC dropped the rest)
                       (flat: the same files at the root)

Each ARC sample has 70 called cells: 60 good ones (many fragments, enriched at TSSs), 5 with
few fragments and 5 without TSS enrichment (PLANTED), plus 40 uncalled background barcodes.
"""
from __future__ import annotations

import argparse
import gzip
import struct
from pathlib import Path

import numpy as np
import pandas as pd

CHROMS = {"chr1": 3_000_000, "chr2": 2_000_000, "chrM": 16_569}
N_GENES = 60
N_GOOD, N_LOW_FRAG, N_LOW_TSSE, N_BACKGROUND = 60, 5, 5, 40
N_PLAIN = 220          # enough cells for Scrublet (atac-prep runs it from 200 cells)
SAMPLES = ("knee1", "knee2")
HUMAN_CONTIGS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]   # as a GRCh38 .tbi lists them
REFERENCE = "/refdata/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"


def genes(rng: np.random.Generator) -> list[tuple[str, int, str, str]]:
    """(chrom, tss, strand, name), spread over chr1 and chr2."""
    out = []
    for i in range(N_GENES):
        chrom = "chr1" if i % 3 else "chr2"
        tss = int(rng.integers(20_000, CHROMS[chrom] - 20_000))
        out.append((chrom, tss, "+" if i % 2 else "-", f"GENE{i:02d}"))
    return sorted(out)


def write_genome(root: Path, gene_list) -> Path:
    d = root / "genomes" / "GRCh38"
    d.mkdir(parents=True, exist_ok=True)
    (d / "chrom_sizes.tsv").write_text("".join(f"{c}\t{n}\n" for c, n in CHROMS.items()))
    with gzip.open(d / "annotation.gff3.gz", "wt") as fh:
        fh.write("##gff-version 3\n")
        for i, (chrom, tss, strand, name) in enumerate(gene_list):
            start, end = (tss, tss + 5000) if strand == "+" else (tss - 5000, tss)
            gid, tid = f"ENSG{i:011d}.1", f"ENST{i:011d}.1"
            common = f"gene_id={gid};gene_type=protein_coding;gene_name={name}"
            fh.write(f"{chrom}\tHAVANA\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gid};{common}\n")
            fh.write(f"{chrom}\tHAVANA\ttranscript\t{start}\t{end}\t.\t{strand}\t.\tID={tid};Parent={gid};{common};"
                     f"transcript_id={tid};transcript_type=protein_coding;transcript_name={name}-201\n")
            fh.write(f"{chrom}\tHAVANA\texon\t{start}\t{end}\t.\t{strand}\t.\tID=exon:{tid}:1;Parent={tid};{common};"
                     f"transcript_id={tid};transcript_type=protein_coding\n")
    return d.parent


def barcodes(n: int, rng: np.random.Generator) -> list[str]:
    seen: set[str] = set()
    while len(seen) < n:
        seen.add("".join(rng.choice(list("ACGT"), 16)) + "-1")
    return sorted(seen)


def cell_fragments(bc: str, n: int, tss_frac: float, gene_list, rng) -> pd.DataFrame:
    """n fragments of one barcode: a share centred on TSSs, 3 % on chrM, the rest uniform."""
    u = rng.random(n)
    length = np.where(u < 0.65, rng.integers(40, 150, n), np.where(u < 0.93, rng.integers(180, 260, n),
                                                                   rng.integers(350, 450, n)))
    kind = rng.random(n)
    g = rng.integers(len(gene_list), size=n)
    gchrom = np.array([x[0] for x in gene_list])[g]
    gtss = np.array([x[1] for x in gene_list])[g]
    uchrom = np.where(rng.random(n) < 0.6, "chr1", "chr2")
    ustart = (rng.random(n) * (np.where(uchrom == "chr1", CHROMS["chr1"], CHROMS["chr2"]) - 500)).astype(int)
    chrom = np.where(kind < 0.03, "chrM", np.where(kind < 0.03 + tss_frac, gchrom, uchrom))
    start = np.where(kind < 0.03, rng.integers(0, CHROMS["chrM"] - 500, n),
                     np.where(kind < 0.03 + tss_frac,
                              np.maximum(0, gtss + rng.normal(0, 150, n).astype(int) - length // 2), ustart))
    return pd.DataFrame({"chrom": chrom, "start": start, "end": start + length, "barcode": bc,
                         "count": np.where(rng.random(n) < 0.1, 2, 1)})


def write_fragments(path: Path, frags: list[pd.DataFrame], header: bool) -> None:
    df = pd.concat(frags, ignore_index=True)
    df["order"] = df["chrom"].map({c: i for i, c in enumerate(CHROMS)})
    df = df.sort_values(["order", "start", "end", "barcode"]).drop(columns="order")
    with gzip.open(path, "wt") as fh:
        if header:
            fh.write("# id=knee\n# description=synthetic\n#\n# pipeline_name=cellranger-arc\n"
                     "# pipeline_version=cellranger-arc-2.0.2\n#\n"
                     f"# reference_path={REFERENCE}\n# reference_fasta_hash=0\n"
                     "# reference_gtf_hash=0\n# reference_version=2020-A\n# mkref_version=cellranger-arc-2.0.0\n#\n"
                     + "".join(f"# primary_contig={c}\n" for c in CHROMS))
        df.to_csv(fh, sep="\t", header=False, index=False)


def write_tbi(path: Path, names: list[str]) -> None:
    """A tabix index header (the sequence names; no bins): what the processor reads from a .tbi."""
    nm = b"".join(n.encode() + b"\0" for n in names)
    head = b"TBI\1" + struct.pack("<8i", len(names), 0, 1, 2, 3, ord("#"), 0, len(nm)) + nm
    with gzip.open(path, "wb") as fh:
        fh.write(head + struct.pack("<i", 0) * len(names))


def write_rna(path: Path, sample_id: str, bcs: list[str], gene_list, rng) -> None:
    """An sc-qc-like RNA AnnData: counts in layers['counts'], log-normalized X, obs_names <sample_id>:<barcode>."""
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp

    X = sp.csr_matrix(rng.poisson(2.0, size=(len(bcs), len(gene_list))).astype(np.float32))
    a = ad.AnnData(X=X.copy(), obs=pd.DataFrame(index=[f"{sample_id}:{b}" for b in bcs]),
                   var=pd.DataFrame(index=[g[3] for g in gene_list]))
    a.layers["counts"] = X
    a.obs["sample_id"] = pd.Categorical([sample_id] * len(bcs))
    a.uns["pennsieve"] = {"organism": "human", "chemistry": "multiome"}
    a.write_h5ad(path)


def write_rna_h5(path: Path) -> None:
    """A stand-in for ARC's filtered_feature_bc_matrix.h5 (atac-prep must ignore it)."""
    import h5py

    with h5py.File(path, "w") as h:
        h.attrs["filetype"] = "matrix"


PLANTED: dict[str, dict[str, list[str]]] = {}


def make_fixture(out: Path, layout: str = "both", seed: int = 0) -> dict[str, Path]:
    rng = np.random.default_rng(seed)
    gene_list = genes(rng)
    write_genome(out, gene_list)
    trees: dict[str, Path] = {"genomes": out / "genomes"}
    per_sample = {}
    for sid in SAMPLES:
        bcs = barcodes(N_GOOD + N_LOW_FRAG + N_LOW_TSSE + N_BACKGROUND, rng)
        good, low_frag = bcs[:N_GOOD], bcs[N_GOOD:N_GOOD + N_LOW_FRAG]
        low_tsse = bcs[N_GOOD + N_LOW_FRAG:N_GOOD + N_LOW_FRAG + N_LOW_TSSE]
        background = bcs[N_GOOD + N_LOW_FRAG + N_LOW_TSSE:]
        frags = []
        for b in good:
            frags.append(cell_fragments(b, int(rng.integers(2200, 3200)), 0.35, gene_list, rng))
        for b in low_frag:
            frags.append(cell_fragments(b, 400, 0.35, gene_list, rng))
        for b in low_tsse:
            frags.append(cell_fragments(b, 2500, 0.0, gene_list, rng))
        for b in background:
            frags.append(cell_fragments(b, int(rng.integers(5, 60)), 0.05, gene_list, rng))
        PLANTED[sid] = {"good": good, "low_fragments": low_frag, "low_tsse": low_tsse, "background": background}
        per_sample[sid] = frags
    plain_bcs = barcodes(N_PLAIN + 20, rng)
    plain = []
    for b in plain_bcs[:N_PLAIN]:
        plain.append(cell_fragments(b, int(rng.integers(1300, 1800)), 0.35, gene_list, rng))
    for b in plain_bcs[N_PLAIN:]:
        plain.append(cell_fragments(b, 50, 0.35, gene_list, rng))
    PLANTED["S3"] = {"good": plain_bcs[:N_PLAIN], "background": plain_bcs[N_PLAIN:]}

    for lay in (("linked", "flat") if layout == "both" else (layout,)):
        root = out / lay
        for sid in SAMPLES:
            outs = root / sid / "outs"
            outs.mkdir(parents=True, exist_ok=True)
            write_fragments(outs / "atac_fragments.tsv.gz", list(per_sample[sid]), header=True)
            p = PLANTED[sid]
            with open(outs / "per_barcode_metrics.csv", "w") as fh:
                fh.write("barcode,gex_barcode,atac_barcode,is_cell,excluded_reason,atac_fragments\n")
                for b in p["good"] + p["low_fragments"] + p["low_tsse"]:
                    fh.write(f"{b},X,X,1,0,1000\n")
                for b in p["background"]:
                    fh.write(f"{b},X,X,0,0,10\n")
            (outs / "atac_peaks.bed").write_text("chr1\t100\t600\nchr2\t100\t600\n")
            write_rna_h5(outs / "filtered_feature_bc_matrix.h5")
            cells_dir = root / "cells" if lay == "linked" else root
            cells_dir.mkdir(parents=True, exist_ok=True)
            # RNA QC kept 50 of the 60 good cells and every planted bad one
            kept = p["good"][:50] + p["low_fragments"] + p["low_tsse"]
            write_rna(cells_dir / f"{sid}.h5ad", sid, kept, gene_list, rng)
        (root / "plain").mkdir(parents=True, exist_ok=True)
        write_fragments(root / "plain" / "S3_fragments.tsv.gz", list(plain), header=False)
        write_tbi(root / "plain" / "S3_fragments.tsv.gz.tbi", HUMAN_CONTIGS)
        (root / "samples.csv").write_text("sample_id,path,condition\nknee1,knee1,OA\nknee2,knee2/outs,control\n"
                                          "S3,plain/S3_fragments.tsv.gz,OA\n")
        trees[lay] = root
    return trees


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path)
    ap.add_argument("--layout", choices=("linked", "flat", "both"), default="both")
    args = ap.parse_args()
    for name, path in make_fixture(args.out, args.layout).items():
        print(f"{name}: {path}")
