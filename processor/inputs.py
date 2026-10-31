"""
Finding the inputs under INPUT_DIR (symlinks followed): fragment files and the Cell Ranger cell
calls next to them (port `fragments`), the RNA side's <sample_id>.h5ad (port `cells`, in
INPUT_DIR/cells/ when linked, else anywhere), the optional sample sheet and SPARC
samples/subjects sheets. Sample ids and metadata columns follow docs/conventions.md of
workflow-single-cell, so the ATAC files get the same sample_id as sc-prep gives the RNA half.
"""
from __future__ import annotations

import csv
import logging
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

log = logging.getLogger("atac-prep")

SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9_-]+$")
FRAGMENTS_SUFFIX = "fragments.tsv.gz"
# Cell Ranger cell calls: ARC per_barcode_metrics.csv (is_cell), ATAC singlecell.csv (is__cell_barcode)
CELL_CALLS = {"per_barcode_metrics.csv": "is_cell", "singlecell.csv": "is__cell_barcode"}
OUTS_NAMES = {"outs", "atac", "output", "outputs"}
MAX_SHEET_BYTES = 20 * 1024 * 1024


class InputError(ValueError):
    pass


def sanitise(name: str) -> str:
    s = re.sub(r"[^A-Za-z0-9_-]+", "_", name).strip("_")
    return re.sub(r"_+", "_", s) or "sample"


def snake(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(name).strip().lower()).strip("_")


def allowed_roots(root: Path) -> list[Path]:
    """INPUT_DIR itself and the targets of its top-level symlinks (how the platform stages inputs)."""
    return [Path(os.path.realpath(root))] + [Path(os.path.realpath(e)) for e in root.iterdir() if e.is_symlink()]


def walk(root: Path) -> list[Path]:
    """Every regular file under root, following symlinks, sorted; hidden names and loops skipped, and
    anything whose resolved path escapes INPUT_DIR (and its staged entries) left out."""
    found, seen = [], set()
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        real = os.path.realpath(dirpath)
        if real in seen:
            dirnames[:] = []
            continue
        seen.add(real)
        dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
        found.extend(Path(dirpath) / f for f in sorted(filenames)
                     if not f.startswith(".") and (Path(dirpath) / f).is_file())
    allowed = allowed_roots(root)
    return [p for p in found if any(Path(os.path.realpath(p)).is_relative_to(a) for a in allowed)]


def _prefix(name: str, suffix: str) -> str:
    """The sample prefix of a file name: S1_atac_fragments.tsv.gz -> S1, atac_fragments.tsv.gz -> ''."""
    p = name[: -len(suffix)]
    p = re.sub(r"(^|[_.-])atac[_.-]?$", "", p)
    return p.strip("_-.")


@dataclass
class Sample:
    fragments: Path
    folder: Path
    default_id: str
    cell_calls: Path | None = None
    tbi: Path | None = None
    peaks_bed: Path | None = None
    sample_id: str = ""
    metadata: dict[str, str] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    def rel(self, root: Path) -> str:
        return str(self.fragments.relative_to(root))


@dataclass
class Found:
    samples: list[Sample]
    rna: dict[str, Path]                 # sample_id -> RNA .h5ad (the `cells` port)
    sheet: Path | None
    sparc: dict[str, Path]               # "samples" / "subjects" -> file
    ignored: list[dict]
    skipped: list[dict]


def default_id(folder: Path, prefix: str, root: Path) -> str:
    if prefix:
        return sanitise(prefix)
    p = folder
    while p.name in OUTS_NAMES and p != root:
        p = p.parent
    return sanitise(p.name) if p != root else "sample"


def find_samples(files: list[Path], root: Path) -> list[Sample]:
    by_folder: dict[Path, list[Path]] = {}
    for f in files:
        by_folder.setdefault(f.parent, []).append(f)
    out = []
    for folder, fs in sorted(by_folder.items()):
        frags = [f for f in fs if f.name.endswith(FRAGMENTS_SUFFIX)]
        calls = {}
        for f in fs:
            for suffix in CELL_CALLS:
                if f.name.endswith(suffix):
                    calls[_prefix(f.name, suffix)] = f
        for frag in frags:
            prefix = _prefix(frag.name, FRAGMENTS_SUFFIX)
            s = Sample(fragments=frag, folder=folder, default_id=default_id(folder, prefix, root))
            s.cell_calls = calls.get(prefix) or (next(iter(calls.values())) if len(frags) == 1 and len(calls) == 1 else None)
            tbi = frag.with_name(frag.name + ".tbi")
            s.tbi = tbi if tbi.is_file() else None
            bed = [f for f in fs if f.name.endswith("peaks.bed") and _prefix(f.name, "peaks.bed") == prefix]
            s.peaks_bed = bed[0] if bed else None
            out.append(s)
    return out


# --- sample sheet --------------------------------------------------------------------------

def read_table(path: Path) -> list[dict[str, str]]:
    """Rows of a .csv / .tsv / .xlsx sheet as stripped strings (header row first)."""
    if path.stat().st_size > MAX_SHEET_BYTES:
        raise InputError(f"{path.name} is larger than {MAX_SHEET_BYTES // 2**20} MB")
    if path.suffix.lower() == ".xlsx":
        from openpyxl import load_workbook

        wb = load_workbook(path, read_only=True, data_only=True)
        rows = [["" if v is None else str(v).strip() for v in r] for r in wb.worksheets[0].iter_rows(values_only=True)]
        wb.close()
    else:
        with open(path, newline="", encoding="utf-8-sig") as fh:
            rows = [[v.strip() for v in r] for r in csv.reader(fh, delimiter="\t" if path.suffix.lower() == ".tsv" else ",")]
    rows = [r for r in rows if any(r)]
    if not rows:
        return []
    header = rows[0]
    return [dict(zip(header, r + [""] * (len(header) - len(r)))) for r in rows[1:]]


def path_matches(sheet_path: str, s: Sample, root: Path) -> bool:
    """A sheet `path` names a sample when it equals or ends its fragment file's path, its folder,
    or the folder above outs/."""
    want = sheet_path.replace("\\", "/").strip().strip("/")
    if not want:
        return False
    forms = {str(s.fragments.relative_to(root)), str(s.folder.relative_to(root))}
    p = s.folder
    while p.name in OUTS_NAMES and p.parent != root and p != root:
        p = p.parent
        forms.add(str(p.relative_to(root)))
    return any(f == want or f.endswith("/" + want) for f in forms)


def apply_sheet(samples: list[Sample], sheet: Path, root: Path, skipped: list[dict], ignored: list[dict]) -> list[Sample]:
    rows = read_table(sheet)
    if rows and not {"sample_id", "path"} <= set(rows[0]):
        raise InputError(f"{sheet.name} needs sample_id and path columns (has: {', '.join(rows[0])})")
    chosen: dict[int, Sample] = {}
    for row in rows:
        hits = [s for s in samples if path_matches(row["path"], s, root)]
        if len(hits) != 1:
            skipped.append({"sample_id": row.get("sample_id", ""), "path": row["path"],
                            "reason": f"sample sheet path matches {len(hits)} fragment files"})
            continue
        s = hits[0]
        if id(s) in chosen:
            skipped.append({"sample_id": row.get("sample_id", ""), "path": row["path"],
                            "reason": "sample sheet row names a sample an earlier row already named"})
            continue
        sid = row.get("sample_id") or s.default_id
        if not SAMPLE_ID_RE.match(sid):
            s.notes.append(f"sample_id {sid!r} sanitised to {sanitise(sid)!r}")
            sid = sanitise(sid)
        s.sample_id = sid
        s.metadata.update({k: v for k, v in row.items() if k not in ("sample_id", "path") and k})
        chosen[id(s)] = s
    for s in samples:
        if id(s) not in chosen:
            ignored.append({"path": s.rel(root), "reason": "fragment file not in the sample sheet"})
    return [s for s in samples if id(s) in chosen]


# --- SPARC dataset metadata ----------------------------------------------------------------

def sparc_metadata(samples: list[Sample], files: dict[str, Path], root: Path, columns: list[str]) -> None:
    """Join SPARC samples/subjects rows onto each sample under a sam-… folder (columns snake_case,
    subject columns prefixed only on a clash, empty columns dropped, only `columns` kept)."""
    sam_rows = {r.get("sample id", ""): r for r in read_table(files["samples"])} if "samples" in files else {}
    sub_rows = {r.get("subject id", ""): r for r in read_table(files["subjects"])} if "subjects" in files else {}

    def tidy(rows: dict) -> dict[str, dict[str, str]]:
        keep = {k for r in rows.values() for k, v in r.items() if v}
        return {i: {snake(k): v for k, v in r.items() if k in keep} for i, r in rows.items()}

    sam_rows, sub_rows = tidy(sam_rows), tidy(sub_rows)
    for s in samples:
        parts = s.fragments.relative_to(root).parts[:-1]
        sam = next((p for p in reversed(parts) if p.startswith("sam-")), None)
        row = dict(sam_rows.get(sam, {})) if sam else {}
        sub_id = row.get("subject_id") or next((p for p in reversed(parts) if p.startswith("sub-")), None)
        sub = sub_rows.get(sub_id, {}) if sub_id else {}
        for k, v in sub.items():
            row[f"subject_{k}" if k in row and k != "subject_id" else k] = v
        for col in columns:
            if col in row and col not in s.metadata:
                s.metadata[col] = row[col]


def add_condition(samples: list[Sample], condition_column: str) -> None:
    for s in samples:
        if "condition" in s.metadata:
            continue
        cols = ([condition_column] if condition_column != "auto"
                else ["subject_experimental_group", "sample_experimental_group"])
        for c in cols:
            if s.metadata.get(c):
                s.metadata["condition"] = s.metadata[c]
                break


# --- everything ----------------------------------------------------------------------------

def _is_sparc(path: Path, kind: str) -> bool:
    try:
        rows = read_table(path)
    except (InputError, OSError, ValueError, UnicodeDecodeError):
        return False
    return bool(rows) and f"{kind} id" in rows[0]


def discover(root: Path, sheet_name: str, metadata_columns: list[str], condition_column: str) -> Found:
    files = walk(root)
    cells_dir = root / "cells"
    rna_files = [f for f in files if f.suffix.lower() == ".h5ad"
                 and (not cells_dir.is_dir() or f.is_relative_to(cells_dir))]
    rna: dict[str, Path] = {}
    for f in rna_files:
        sid = f.stem
        if sid in rna:
            raise InputError(f"two RNA files for sample {sid}: {rna[sid].relative_to(root)}, {f.relative_to(root)}")
        rna[sid] = f
    samples = find_samples([f for f in files if not f.is_relative_to(cells_dir)], root)
    skipped: list[dict] = []
    ignored: list[dict] = []

    shallow = lambda ps: sorted(ps, key=lambda p: (len(p.relative_to(root).parts), str(p)))  # noqa: E731
    sparc = {}
    for kind in ("samples", "subjects"):
        hits = [f for f in files if f.stem.lower() == kind and f.suffix.lower() in (".xlsx", ".csv", ".tsv")
                and _is_sparc(f, kind[:-1])]
        if hits:
            sparc[kind] = shallow(hits)[0]
    sheets = [f for f in files if f.name == sheet_name and f not in sparc.values()]
    sheet = shallow(sheets)[0] if sheets else None

    for s in samples:
        s.sample_id = s.default_id
    if sheet is not None:
        samples = apply_sheet(samples, sheet, root, skipped, ignored)
    if sparc:
        sparc_metadata(samples, sparc, root, metadata_columns)
    add_condition(samples, condition_column)
    seen: dict[str, Sample] = {}
    for s in samples:
        if s.sample_id in seen:
            raise InputError(f"two fragment files would be sample {s.sample_id}: {seen[s.sample_id].rel(root)} and "
                             f"{s.rel(root)}; name them in a sample sheet")
        seen[s.sample_id] = s
    return Found(samples=samples, rna=rna, sheet=sheet, sparc=sparc, ignored=ignored, skipped=skipped)
