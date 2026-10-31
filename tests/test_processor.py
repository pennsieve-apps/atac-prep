"""The synthetic fixture through the processor: both input layouts, with and without the RNA side,
SPARC metadata, genome detection and parameter checks."""
from __future__ import annotations

import csv
import gzip
import json
import shutil
from pathlib import Path

import pytest
import snapatac2 as snap

from processor import genome as G
from processor import main as app
from tests.make_fixture import CHROMS, HUMAN_CONTIGS, PLANTED, make_fixture, write_tbi


@pytest.fixture(scope="session")
def fx(tmp_path_factory) -> dict[str, Path]:
    return make_fixture(tmp_path_factory.mktemp("fixture"), "both")


def run(input_dir: Path, output_dir: Path, fx, monkeypatch, **params) -> int:
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    monkeypatch.setenv("GENOMES_DIR", str(fx["genomes"]))
    monkeypatch.setenv("MPLCONFIGDIR", "/tmp/matplotlib")
    for name in app.PARAMETERS:
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv("THREADS", "2")
    for name, value in params.items():
        monkeypatch.setenv(name, str(value))
    return app.main()


def rows(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def kept_barcodes(out: Path, sid: str) -> list[str]:
    a = snap.read(out / "atac" / f"{sid}.h5ad", backed="r")
    try:
        return [n.split(":", 1)[1] for n in a.obs_names]
    finally:
        a.close()


@pytest.mark.parametrize("layout", ["linked", "flat"])
def test_rna_linked(fx, tmp_path, monkeypatch, layout):
    out = tmp_path / "out"
    assert run(fx[layout], out, fx, monkeypatch) == 0
    assert sorted(p.name for p in (out / "atac").iterdir()) == ["knee1.h5ad", "knee2.h5ad", "provenance.json"]
    for sid in ("knee1", "knee2"):
        # RNA QC kept 50 good cells + the 10 planted bad ones; ATAC QC drops the bad ones
        assert kept_barcodes(out, sid) == sorted(PLANTED[sid]["good"][:50])
        a = snap.read(out / "atac" / f"{sid}.h5ad", backed="r")
        try:
            assert all(n.startswith(f"{sid}:") for n in a.obs_names)
            assert {"n_fragment", "frac_dup", "frac_mito", "tsse", "sample_id", "condition"} <= set(a.obs[:].columns)
            assert "fragment_paired" in a.obsm
            n_tiles = sum(-(-n // 500) for c, n in CHROMS.items() if c != "chrM")
            assert a.shape == (50, n_tiles)
            pen = a.uns["pennsieve"]
            assert pen["genome"] == "GRCh38" and pen["cells_from"] == "rna" and pen["sample_id"] == sid
        finally:
            a.close()
        assert (out / "qc" / f"{sid}_tsse.png").exists() and (out / "qc" / f"{sid}_fragment_sizes.png").exists()
    qc = rows(out / "qc" / "atac_qc.csv")
    k1 = {r["obs_id"].split(":")[1]: r for r in qc if r["sample_id"] == "knee1"}
    assert len(k1) == 60
    assert all("low_fragments" in k1[b]["reason"] for b in PLANTED["knee1"]["low_fragments"])
    assert all("low_tsse" in k1[b]["reason"] for b in PLANTED["knee1"]["low_tsse"])
    # S3 has no RNA file while cells is linked: skipped, not processed ATAC-only
    prov = json.loads((out / "qc" / "provenance.json").read_text())
    s3 = next(s for s in prov["samples"] if s["sample_id"] == "S3")
    assert s3["status"] == "failed" and "S3.h5ad" in s3["reason"]
    assert prov["parameters"]["MIN_TSSE"] == 5 and "10.1038/s41592-023-02139-9" in [c["doi"] for c in prov["citations"]]
    summary = json.loads((out / "qc" / "summary.json").read_text())
    assert summary["n_samples"] == 2 and summary["n_cells_kept"] == 100


def test_atac_only(fx, tmp_path, monkeypatch):
    """Without `cells`: Cell Ranger's cell calls for the ARC sample, a fragment threshold and
    Scrublet for the bare file; the sample sheet is optional."""
    src = tmp_path / "in"
    shutil.copytree(fx["linked"] / "knee1", src / "knee1")
    shutil.copytree(fx["linked"] / "plain", src / "plain")
    out = tmp_path / "out"
    assert run(src, out, fx, monkeypatch) == 0
    assert kept_barcodes(out, "knee1") == sorted(PLANTED["knee1"]["good"])
    s3 = set(kept_barcodes(out, "S3"))
    assert s3 <= set(PLANTED["S3"]["good"]) and len(s3) > 0.9 * len(PLANTED["S3"]["good"])
    qc = [r for r in rows(out / "qc" / "atac_qc.csv") if r["sample_id"] == "S3"]
    assert not {r["obs_id"].split(":")[1] for r in qc} & set(PLANTED["S3"]["background"])
    summary_rows = {r["sample_id"]: r for r in rows(out / "qc" / "qc_summary.csv")}
    assert int(summary_rows["S3"]["n_cells_after"]) == len(s3)
    assert int(summary_rows["knee1"]["n_removed_low_fragments"]) == 5
    summary = json.loads((out / "qc" / "summary.json").read_text())
    assert summary["samples"]["knee1"]["cells_from"] == "cellranger"
    assert summary["samples"]["S3"]["cells_from"] == "min_fragments" and summary["samples"]["S3"]["scrublet"]["ran"]
    assert "chromosome names" in summary["samples"]["S3"]["genome_from"]


def test_sparc_metadata(fx, tmp_path, monkeypatch):
    from openpyxl import Workbook

    src = tmp_path / "in"
    sam = src / "primary" / "sub-1" / "sam-1"
    sam.mkdir(parents=True)
    shutil.copy(fx["linked"] / "knee1" / "outs" / "atac_fragments.tsv.gz", sam)
    for name, header, values in (("samples", ["sample id", "subject id", "sample experimental group", "notes"],
                                  ["sam-1", "sub-1", "OA", ""]),
                                 ("subjects", ["subject id", "sex", "species", "age", "BMI"],
                                  ["sub-1", "Female", "Homo sapiens", "54 years", "31"])):
        wb = Workbook()
        wb.active.append(header)
        wb.active.append(values)
        wb.save(src / f"{name}.xlsx")
    out = tmp_path / "out"
    assert run(src, out, fx, monkeypatch, GENOME="auto") == 0
    a = snap.read(out / "atac" / "sam-1.h5ad", backed="r")
    try:
        obs = a.obs[:].to_pandas()
    finally:
        a.close()
    assert obs["sex"].iloc[0] == "Female" and obs["sample_experimental_group"].iloc[0] == "OA"
    assert obs["condition"].iloc[0] == "OA" and obs["age"].iloc[0] == "54 years"
    assert "bmi" not in obs.columns and "notes" not in obs.columns   # only METADATA_COLUMNS reach the cells


def test_symlink_out_of_input_ignored(fx, tmp_path, monkeypatch):
    # links directly in INPUT_DIR or a port folder are the platform's; one deeper that leaves is not
    src = tmp_path / "in"
    (src / "a" / "b").mkdir(parents=True)
    (src / "a" / "b" / "atac_fragments.tsv.gz").symlink_to(fx["linked"] / "knee1" / "outs" / "atac_fragments.tsv.gz")
    assert run(src, tmp_path / "out", fx, monkeypatch) == 1


def test_genome_detection(tmp_path, fx):
    names = HUMAN_CONTIGS
    write_tbi(tmp_path / "x.tbi", names)
    assert G.tbi_names(tmp_path / "x.tbi") == names
    frag = tmp_path / "f.tsv.gz"
    with gzip.open(frag, "wt") as fh:
        fh.write("# reference_path=/ref/refdata-cellranger-arc-mm10-2020-A-2.0.0\nchr1\t1\t100\tA-1\t1\n")
    with pytest.raises(G.GenomeError, match="mm10"):
        G.resolve("auto", frag, None, None, None)
    with gzip.open(frag, "wt") as fh:
        fh.write("chr1\t1\t100\tA-1\t1\n")
    with pytest.raises(G.GenomeError, match="set GENOME"):
        G.resolve("auto", frag, None, None, None)


@pytest.mark.parametrize("params", [{"GENOME": "hg19"}, {"MIN_FRAGMENTS": "5000", "MAX_FRAGMENTS": "100"},
                                    {"MIN_TSSE": "x"}, {"MIN_TSSE": "-1"}])
def test_bad_parameters(fx, tmp_path, monkeypatch, params):
    assert run(fx["linked"], tmp_path / "out", fx, monkeypatch, **params) == 2


def test_no_input(fx, tmp_path, monkeypatch):
    (tmp_path / "in").mkdir()
    assert run(tmp_path / "in", tmp_path / "out", fx, monkeypatch) == 1
