"""
atac-prep: the first ATAC step of the Multiome Gene Regulation workflow.

Finds the fragment files under INPUT_DIR (Cell Ranger ARC / ATAC outputs or any 10x-style
fragments.tsv.gz), optionally keeps only the nuclei that passed RNA QC (port `cells`), and writes
per sample

  atac/<sample_id>.h5ad   SnapATAC2 AnnData: fragments (obsm["fragment_paired"]), 500 bp tile matrix
                          in X, QC in obs, obs_names <sample_id>:<barcode>, uns["pennsieve"]
  qc/atac_qc.csv          every evaluated barcode with its QC metrics and why it was removed
  qc/<sample_id>_tsse.png, qc/<sample_id>_fragment_sizes.png
  qc/summary.json, provenance.json (in both output folders)

One sample at a time; see processor/sample.py for the method.
"""
from __future__ import annotations

import csv
import json
import logging
import math
import os
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("atac-prep")
for noisy in ("matplotlib", "PIL", "fontTools", "numba", "h5py", "anndata"):
    logging.getLogger(noisy).setLevel(logging.WARNING)

APP_NAME = "atac-prep"
REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_METADATA = ("subject_id,subject_experimental_group,sample_experimental_group,sex,age,species,"
                    "sample_anatomical_location,laterality")
DEFAULTS = {
    "SAMPLE_SHEET": "samples.csv",
    "METADATA_COLUMNS": DEFAULT_METADATA,
    "CONDITION_COLUMN": "auto",
    "GENOME": "auto",
    "MIN_FRAGMENTS": "1000",
    "MAX_FRAGMENTS": "100000",
    "MIN_TSSE": "5",
    "SEED": "42",
    "THREADS": "0",
}
PARAMETERS = tuple(DEFAULTS)
QC_COLUMNS = ["sample_id", "obs_id", "n_fragment", "frac_dup", "frac_mito", "tsse", "passed", "reason"]
REASONS = ("low_fragments", "high_fragments", "low_tsse", "doublet")
CITATIONS = {   # each checked on Crossref and PubMed (README "Citations")
    "snapatac2": {"doi": "10.1038/s41592-023-02139-9", "pmid": "38191932", "text": "Zhang et al. Nat Methods 2024 (SnapATAC2)"},
    "gencode": {"doi": "10.1093/nar/gkac1071", "pmid": "36420896", "text": "Frankish et al. Nucleic Acids Res 2023 (GENCODE)"},
    "best_practices": {"doi": "10.1038/s41576-023-00586-w", "pmid": "37002403", "text": "Heumos et al. Nat Rev Genet 2023"},
    "scrublet": {"doi": "10.1016/j.cels.2018.11.005", "pmid": "30954476", "text": "Wolock et al. Cell Syst 2019 (Scrublet)"},
}


class ParameterError(ValueError):
    pass


def read_parameters(env=os.environ) -> dict:
    from processor.genome import GENOMES

    raw = {k: (env.get(k) if env.get(k) not in (None, "") else v) for k, v in DEFAULTS.items()}
    p: dict = {"SAMPLE_SHEET": raw["SAMPLE_SHEET"].strip(), "CONDITION_COLUMN": raw["CONDITION_COLUMN"].strip()}
    p["METADATA_COLUMNS"] = [c.strip() for c in raw["METADATA_COLUMNS"].split(",") if c.strip()]
    genome = raw["GENOME"].strip()
    match = {g.lower(): g for g in GENOMES + ("auto",)}.get(genome.lower())
    if match is None:
        raise ParameterError(f"GENOME must be auto, {' or '.join(GENOMES)}; got {genome!r}")
    p["GENOME"] = match
    try:
        p["MIN_FRAGMENTS"] = int(raw["MIN_FRAGMENTS"])
        p["MAX_FRAGMENTS"] = int(raw["MAX_FRAGMENTS"])
        p["MIN_TSSE"] = float(raw["MIN_TSSE"])
        p["SEED"] = int(raw["SEED"])
        p["THREADS"] = int(raw["THREADS"])
    except ValueError as exc:
        raise ParameterError(f"bad numeric parameter: {exc}") from exc
    if p["MIN_FRAGMENTS"] < 1 or p["MAX_FRAGMENTS"] <= p["MIN_FRAGMENTS"]:
        raise ParameterError("need 1 <= MIN_FRAGMENTS < MAX_FRAGMENTS")
    if not math.isfinite(p["MIN_TSSE"]) or p["MIN_TSSE"] < 0 or p["THREADS"] < 0:
        raise ParameterError("MIN_TSSE and THREADS must be >= 0")
    return p


def write_json(path: Path, data) -> None:
    def clean(o):
        if isinstance(o, dict):
            return {str(k): clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(v) for v in o]
        if isinstance(o, float) and not math.isfinite(o):
            return None
        if hasattr(o, "item"):
            return clean(o.item())
        if isinstance(o, Path):
            return str(o)
        return o

    with open(path, "w") as fh:
        json.dump(clean(data), fh, indent=2)
        fh.write("\n")


def tool_versions() -> dict:
    from importlib.metadata import PackageNotFoundError, version

    out = {"python": platform.python_version()}
    for name in ("snapatac2", "anndata", "numpy", "polars", "pandas", "scipy", "matplotlib", "openpyxl"):
        try:
            out[name] = version(name)
        except PackageNotFoundError:
            out[name] = None
    return out


def peak_rss_mb() -> float | None:
    try:
        import resource

        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return round(peak / (1024 * 1024 if sys.platform == "darwin" else 1024), 1)
    except (ImportError, OSError):
        return None


def app_version() -> str:
    try:
        return (REPO_ROOT / "VERSION").read_text().strip()
    except OSError:
        return "unknown"


def main() -> int:
    started = datetime.now(timezone.utc)
    t0 = time.monotonic()
    input_dir, output_dir = os.environ.get("INPUT_DIR"), os.environ.get("OUTPUT_DIR")
    if not input_dir or not output_dir:
        log.error("INPUT_DIR and OUTPUT_DIR must be set")
        return 2
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    if not input_dir.is_dir():
        log.error("INPUT_DIR %s is not a folder", input_dir)
        return 2
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    try:
        params = read_parameters()
    except ParameterError as exc:
        log.error("bad parameter: %s", exc)
        return 2
    threads = params["THREADS"] or int(os.environ.get("OMP_NUM_THREADS") or 2)

    from processor import figures as F
    from processor import genome as G
    from processor import inputs as I
    from processor import sample as S

    try:
        found = I.discover(input_dir, params["SAMPLE_SHEET"], params["METADATA_COLUMNS"], params["CONDITION_COLUMN"])
    except (I.InputError, OSError, ValueError) as exc:
        log.error("cannot read the inputs: %s", exc)
        return 1
    rel = lambda p: str(p.relative_to(input_dir)) if p else None  # noqa: E731
    log.info("input: %d fragment file(s), %d RNA file(s) (cells %s), sample sheet %s, SPARC sheets %s, threads %d",
             len(found.samples), len(found.rna), "linked" if found.rna else "not linked", rel(found.sheet) or "none",
             ", ".join(rel(p) for p in found.sparc.values()) or "none", threads)
    for s in found.skipped:
        log.error("sample sheet row %s (%s): %s", s["sample_id"] or "?", s["path"], s["reason"])
    if not found.samples:
        log.error("no fragment file (*fragments.tsv.gz) found under INPUT_DIR")
        return 1

    atac_dir, qc_dir = output_dir / "atac", output_dir / "qc"
    atac_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)
    version = app_version()
    qc_rows: list[dict] = []
    per_sample: dict[str, dict] = {}
    statuses: list[dict] = list(found.skipped)
    used_scrublet = False
    for s in found.samples:
        t = time.monotonic()
        sid, src = s.sample_id, s.rel(input_dir)
        out = atac_dir / f"{sid}.h5ad"
        for note in s.notes:
            log.warning("%s: %s", sid, note)
        try:
            organism = G.organism_of(s.metadata.get("species") or s.metadata.get("organism"))
            if found.rna:
                if sid not in found.rna:
                    raise ValueError(f"cells is linked but holds no {sid}.h5ad (RNA files: "
                                     + ", ".join(sorted(found.rna)) + ")")
                whitelist, rna_org = S.rna_barcodes(found.rna[sid], sid)
                calls_from, organism = "rna", organism or G.organism_of(rna_org)
            elif s.cell_calls is not None:
                whitelist, calls_from = S.cell_calls(s.cell_calls), "cellranger"
            else:
                whitelist, calls_from = None, "min_fragments"
            if whitelist is not None and not whitelist:
                raise ValueError(f"no cell barcodes in {calls_from} cell calls")
            genome = G.resolve(params["GENOME"], s.fragments, organism, s.tbi, s.peaks_bed)
            result = S.process(s, genome, whitelist, calls_from, params, out, threads, version, src)
            st = result["stats"]
            F.tsse(sid, result["figure"], params, qc_dir / f"{sid}_tsse.png")
            F.fragment_sizes(sid, result["figure"], qc_dir / f"{sid}_fragment_sizes.png")
        except Exception as exc:  # noqa: BLE001 - one sample's failure must not sink the run
            out.unlink(missing_ok=True)
            if isinstance(exc, ValueError):     # an input problem we describe (GenomeError is one): no traceback
                log.error("%s (%s): FAILED: %s", sid, src, exc)
            else:
                log.exception("%s (%s): FAILED: %s", sid, src, exc)
            statuses.append({"sample_id": sid, "path": src, "status": "failed", "reason": f"{type(exc).__name__}: {exc}",
                             "seconds": round(time.monotonic() - t, 1)})
            continue
        qc_rows.extend(result["rows"])
        used_scrublet |= bool(st["scrublet"] and st["scrublet"]["ran"])
        st["seconds"] = round(time.monotonic() - t, 1)
        per_sample[sid] = st
        statuses.append({"sample_id": sid, "path": src, "status": "ok", "output": f"atac/{sid}.h5ad",
                         "n_kept": st["n_kept"], "seconds": st["seconds"], "peak_rss_mb": peak_rss_mb()})
        log.info("%s: %s, cells from %s%s, genome %s (%s): %d evaluated -> %d kept (removed: %s); median %.0f fragments, "
                 "TSSE %.1f; %.1fs", sid, src, calls_from, f" ({len(whitelist)})" if whitelist is not None else "",
                 genome.name, genome.how, st["n_evaluated"], st["n_kept"],
                 ", ".join(f"{k} {v}" for k, v in st["removed"].items()) or "none", st["median_n_fragment"],
                 st["median_tsse"], st["seconds"])

    with open(qc_dir / "atac_qc.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=QC_COLUMNS, lineterminator="\n")
        w.writeheader()
        w.writerows(qc_rows)
    with open(qc_dir / "qc_summary.csv", "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["sample_id", "n_cells_before", "n_cells_after", *(f"n_removed_{r}" for r in REASONS),
                    "median_tsse", "median_fragments", "cells_from", "genome"])
        for sid, st in per_sample.items():
            w.writerow([sid, st["n_evaluated"], st["n_kept"], *(st["removed"].get(r, 0) for r in REASONS),
                        round(st["median_tsse"], 3), st["median_n_fragment"], st["cells_from"], st["genome"]])
    finished = datetime.now(timezone.utc)
    citations = [CITATIONS[k] for k in ("snapatac2", "gencode", "best_practices")]
    if used_scrublet:
        citations.append(CITATIONS["scrublet"])
    summary = {
        "app": APP_NAME, "version": version, "n_samples": len(per_sample),
        "n_failed": sum(1 for s in statuses if s["status"] != "ok"),
        "n_cells_kept": sum(s["n_kept"] for s in per_sample.values()),
        "n_cells_evaluated": sum(s["n_evaluated"] for s in per_sample.values()),
        "thresholds": {k: params[k] for k in ("MIN_FRAGMENTS", "MAX_FRAGMENTS", "MIN_TSSE")},
        "samples": per_sample,
    }
    write_json(qc_dir / "summary.json", summary)
    provenance = {
        "app": APP_NAME, "version": version, "tools": tool_versions(), "parameters": params,
        "inputs": {"fragments": [s.rel(input_dir) for s in found.samples],
                   "cells": sorted(rel(p) for p in found.rna.values()),
                   "sample_sheet": rel(found.sheet), "sparc": {k: rel(v) for k, v in found.sparc.items()},
                   "ignored": found.ignored},
        "started": started.isoformat(timespec="seconds"), "finished": finished.isoformat(timespec="seconds"),
        "seconds": round(time.monotonic() - t0, 1), "threads": threads,
        "samples": statuses, "peak_rss_mb": peak_rss_mb(), "citations": citations,
    }
    for d in (atac_dir, qc_dir):
        write_json(d / "provenance.json", provenance)
    log.info("%d sample(s) written, %d failed, %d cells kept, peak RSS %s MB, %.1fs", len(per_sample),
             summary["n_failed"], summary["n_cells_kept"], provenance["peak_rss_mb"], time.monotonic() - t0)
    if not per_sample:
        log.error("every sample failed")
        return 1
    return 0
