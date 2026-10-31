"""
The genome a fragment file is aligned to, and the annotation baked into the image for it.

GENOMES_DIR (default /opt/genomes, filled at build time) holds one folder per genome with
annotation.gff3.gz (GENCODE basic, as SnapATAC2 ships it) and chrom_sizes.tsv (primary
chromosomes). GENOME=auto decides per sample from, in order: the reference named in the fragment
file's header (Cell Ranger writes `# reference_path=…`), the organism (sample sheet / SPARC
`species`, or the RNA file's uns["pennsieve"]["organism"]), the chromosome names in the tabix
index or Cell Ranger's peaks.bed (chr20-chr22 only in human).
"""
from __future__ import annotations

import gzip
import os
import re
import struct
from dataclasses import dataclass
from pathlib import Path

GENOMES = ("GRCh38", "GRCm39")
ORGANISM_GENOME = {"human": "GRCh38", "mouse": "GRCm39"}
SUPPORTED = {"grch38": "GRCh38", "hg38": "GRCh38", "grcm39": "GRCm39", "mm39": "GRCm39"}
UNSUPPORTED = ("grch37", "hg19", "grcm38", "mm10")
MAX_HEADER_LINES = 200
MAX_LINE = 65536


class GenomeError(ValueError):
    pass


@dataclass
class Genome:
    name: str
    annotation: Path
    chrom_sizes: dict[str, int]
    how: str
    inferred: bool = False      # from organism or chromosome names: cannot tell mm10 from GRCm39, hg19 from GRCh38


def genomes_dir() -> Path:
    return Path(os.environ.get("GENOMES_DIR", "/opt/genomes"))


def load(name: str, how: str, inferred: bool = False) -> Genome:
    d = genomes_dir() / name
    sizes = {}
    for line in (d / "chrom_sizes.tsv").read_text().splitlines():
        chrom, n = line.split("\t")
        sizes[chrom] = int(n)
    return Genome(name=name, annotation=d / "annotation.gff3.gz", chrom_sizes=sizes, how=how, inferred=inferred)


def organism_of(text: str | None) -> str | None:
    t = (text or "").strip().lower()
    if t in ("human", "homo sapiens", "ncbitaxon:9606"):
        return "human"
    if t in ("mouse", "mus musculus", "ncbitaxon:10090"):
        return "mouse"
    return None


def header_reference(fragments: Path) -> str | None:
    """The reference named in the fragment file's '#' header lines, or None."""
    with gzip.open(fragments, "rt", errors="replace") as fh:
        for _ in range(MAX_HEADER_LINES):
            line = fh.readline(MAX_LINE)
            if not line.startswith("#"):
                return None
            m = re.match(r"#\s*reference_(?:path|genome)\s*=\s*(\S+)", line)
            if m:
                return m.group(1)
    return None


def tbi_names(tbi: Path) -> list[str]:
    """Sequence names from a tabix index header."""
    with gzip.open(tbi, "rb") as fh:
        head = fh.read(36)
        if len(head) < 36 or head[:4] != b"TBI\1":
            raise GenomeError(f"{tbi.name} is not a tabix index")
        n_ref, *_rest, l_nm = struct.unpack("<8i", head[4:36])
        if not 0 < l_nm < 1_000_000:
            raise GenomeError(f"{tbi.name}: implausible name block")
        return [n.decode() for n in fh.read(l_nm).split(b"\0") if n][:n_ref]


def bed_names(bed: Path, limit: int = 2_000_000) -> set[str]:
    names = set()
    with open(bed, errors="replace") as fh:
        for _ in range(limit):
            line = fh.readline(MAX_LINE)
            if not line:
                break
            if not line.startswith(("#", "track", "browser")):
                names.add(line.split("\t", 1)[0])
    return names


def resolve(param: str, fragments: Path, organism: str | None, tbi: Path | None, peaks_bed: Path | None) -> Genome:
    if param != "auto":
        return load(param, "GENOME parameter")
    ref = header_reference(fragments)
    if ref:
        low = Path(ref).name.lower()
        for tag, name in SUPPORTED.items():
            if tag in low:
                return load(name, f"fragment header reference {Path(ref).name}")
        if any(tag in low for tag in UNSUPPORTED):
            raise GenomeError(f"aligned to {Path(ref).name}; only GRCh38 and GRCm39 annotations are in the image")
    if organism in ORGANISM_GENOME:
        return load(ORGANISM_GENOME[organism], f"organism {organism}", inferred=True)
    names: set[str] = set()
    if tbi is not None:
        names = set(tbi_names(tbi))
    elif peaks_bed is not None:
        names = bed_names(peaks_bed)
    if names & {"chr20", "chr21", "chr22"}:
        return load("GRCh38", "chromosome names (chr20-22: human)", inferred=True)
    if "chr19" in names:
        return load("GRCm39", "chromosome names (chr1-19, no chr20: mouse)", inferred=True)
    raise GenomeError("cannot tell the genome (no reference in the fragment header, no organism, no .tbi or "
                      "peaks.bed); set GENOME")
