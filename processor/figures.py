"""Per-sample QC figures (PNG, 150 dpi, sample id in the title)."""
from __future__ import annotations

from pathlib import Path

import numpy as np

DPI = 150


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def tsse(sample_id: str, fig: dict, params: dict, path: Path) -> Path:
    """Left: TSS enrichment against fragments per cell, kept and removed, with the thresholds.
    Right: the sample's aggregate insertion profile around TSSs (±2 kb, 21 bp running mean)."""
    plt = _plt()
    f, (a, b) = plt.subplots(1, 2, figsize=(10, 4.2))
    x, y, keep = np.maximum(fig["n_fragment"], 1), fig["tsse"], fig["keep"]
    a.scatter(x[~keep], y[~keep], s=4, c="#c0392b", alpha=0.5, lw=0, label=f"removed ({(~keep).sum():,})")
    a.scatter(x[keep], y[keep], s=4, c="#34495e", alpha=0.5, lw=0, label=f"kept ({keep.sum():,})")
    a.set_xscale("log")
    a.axvline(params["MIN_FRAGMENTS"], c="grey", ls="--", lw=1)
    a.axvline(params["MAX_FRAGMENTS"], c="grey", ls="--", lw=1)
    a.axhline(params["MIN_TSSE"], c="grey", ls="--", lw=1)
    a.set_xlabel("unique fragments per cell")
    a.set_ylabel("TSS enrichment")
    a.legend(loc="upper right", fontsize=8, markerscale=3)
    prof = np.convolve(fig["tss_profile"], np.ones(21) / 21, mode="same")   # as sample.library_tsse
    flank = np.mean(np.r_[prof[:100], prof[-100:]]) or 1.0
    pos = np.arange(len(prof)) - len(prof) // 2
    b.plot(pos, prof / flank, c="#34495e", lw=1)
    b.set_xlabel("distance to TSS (bp)")
    b.set_ylabel("insertions / flank mean")
    b.set_title(f"library TSS enrichment {fig['library_tsse']:.1f}", fontsize=9)
    f.suptitle(f"{sample_id}: TSS enrichment")
    f.tight_layout()
    f.savefig(path, dpi=DPI)
    plt.close(f)
    return path


def fragment_sizes(sample_id: str, fig: dict, path: Path) -> Path:
    """Fragment-length distribution of the kept cells: the nucleosome-free peak below ~150 bp, then
    mono- and di-nucleosome peaks near 200 and 400 bp in a good library."""
    plt = _plt()
    f, a = plt.subplots(figsize=(6, 4))
    d = fig["frag_size"]
    a.plot(np.arange(len(d)), d / max(d.sum(), 1), c="#34495e", lw=1)
    a.set_yscale("log")
    a.set_xlabel("fragment length (bp)")
    a.set_ylabel("fraction of fragments")
    a.set_title(f"{sample_id}: fragment sizes (kept cells)")
    f.tight_layout()
    f.savefig(path, dpi=DPI)
    plt.close(f)
    return path
