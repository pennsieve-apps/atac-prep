# atac-prep

The first ATAC step of the **Multiome Gene Regulation** workflow:

- **B. Multiome Gene Regulation:** `atac-prep → atac-regulome (+ labels from sc-annotate) → sc-report`

A 10x multiome run measures two things in the same nucleus: its RNA and which parts of its
chromatin are open (ATAC). The RNA half goes through sc-prep and sc-qc; this app prepares the ATAC
half so the two describe **the same nuclei under the same names**. It reads the fragment files
(every Tn5-cut DNA fragment with the barcode of its nucleus), keeps the nuclei that passed RNA QC
(or Cell Ranger's cell calls), and measures how good each nucleus' chromatin signal is: how many
unique fragments it has and how strongly they pile up at gene promoters (TSS enrichment; a nucleus
with intact, accessible chromatin has fragments concentrated at transcription start sites, a broken
or empty droplet's fragments are spread evenly). Nuclei that fail are removed and listed with the
reason. What is left is one file per sample that the next step uses for peak calling per cell type.

For the computational analyst: per sample, the called nuclei's fragments sorted by barcode on disk,
`snapatac2.pp.import_fragments` (backed, streaming, written straight to the output `.h5ad`),
`snapatac2.metrics.tsse` (on 250-nucleus subsets, to bound memory) against the GENCODE annotation
baked into the image, cut-offs on
`n_fragment` and `tsse`, `snapatac2.metrics.frag_size_distr`, `snapatac2.pp.add_tile_matrix`
(500 bp), and for ATAC-only runs Scrublet (`select_features` + `scrublet`, nuclei with doublet
probability > 0.5 removed, the rule of `filter_doublets`).
`obs_names` become `<sample_id>:<barcode>`, exactly the RNA side's names (Cell Ranger ARC gives a
nucleus the same barcode in both modalities).

## Method

1. **Samples.** Every `*fragments.tsv.gz` under the input is a sample (`atac_fragments.tsv.gz` of
   Cell Ranger ARC, `fragments.tsv.gz` of Cell Ranger ATAC, `<prefix>_fragments.tsv.gz`). sample_id:
   a sample sheet's `sample_id` (`SAMPLE_SHEET`, columns `sample_id,path,…`; `path` equal to or
   ending the fragment file, its folder or the folder above `outs/`), else the file-name prefix,
   else the folder name above `outs/`, which is how sc-prep names the RNA half of the same ARC run.
   SPARC `samples`/`subjects` sheets (`.xlsx`/`.csv`/`.tsv`) are joined for samples under a `sam-…`
   folder as sc-prep does it; only `METADATA_COLUMNS` reach the cells.
2. **Which barcodes are nuclei** (`uns["pennsieve"]["cells_from"]`):
   - `rna` — the `cells` input is linked: the barcodes of `cells/<sample_id>.h5ad` (the nuclei that
     passed sc-prep/sc-qc, doublets already removed). A sample with no RNA file is skipped.
   - `cellranger` — Cell Ranger's calls next to the fragments: `per_barcode_metrics.csv`
     (`is_cell`, ARC) or `singlecell.csv` (`is__cell_barcode`, ATAC).
   - `min_fragments` — neither: every barcode with at least `MIN_FRAGMENTS` fragments.
3. **Genome** (`GENOME=auto`): the reference in the fragment file's header (Cell Ranger writes
   `# reference_path=…/refdata-cellranger-arc-GRCh38-2020-A-2.0.0`), else the organism (sample
   sheet / SPARC `species`, or the RNA file's `uns["pennsieve"]["organism"]`), else the chromosome
   names in the `.tbi` index or `atac_peaks.bed` (chr20–22 exist only in human). GRCh38 and GRCm39
   are built in (GENCODE v41 / vM30 basic, the annotations SnapATAC2 uses). An mm10 or hg19
   reference named in the header is refused. Organism and chromosome names cannot tell GRCm39 from
   mm10 (most mouse Cell Ranger ARC data use `mm10-2020-A`) or GRCh38 from hg19, so a genome guessed
   that way is logged as a warning and the sample fails ("wrong genome?") when its library TSS
   enrichment is below 2; set `GENOME` when the files carry no header. Malformed fragment lines
   (fewer than five columns, non-numeric start or end) are skipped and counted
   (`summary.json` → `n_malformed_lines`).
4. **ATAC QC.** A nucleus is removed when `n_fragment < MIN_FRAGMENTS` (`low_fragments`),
   `n_fragment > MAX_FRAGMENTS` (`high_fragments`) or `tsse < MIN_TSSE` (`low_tsse`). Defaults
   1,000 / 100,000 / 5: 1,000 fragments and TSSE 5 are SnapATAC2's `filter_cells` defaults (ArchR's
   are 1,000 and 4, Granja et al. 2021); 100,000 is the upper cut of SnapATAC2's standard tutorial.
   The tutorial's 5,000 fragments / TSSE 10 suit fresh PBMCs; nuclei from frozen DRG, ganglia,
   cartilage or synovium are sparser, and with RNA linked the nuclei are already vetted by RNA QC,
   so the permissive end is the default. The QC figures show where each sample sits.
5. **Doublets** (only when `cells` is not linked, and from 200 cells): Scrublet on the tile matrix,
   probability > 0.5 removed (`doublet`). With RNA linked, sc-qc removed doublets on the RNA side.
6. **Output per sample**: the kept nuclei, input order.

Deviations from the best-practices guide (Heumos et al. 2023): thresholds are fixed values, not
MAD-based, because TSS enrichment is a ratio with a known scale (values below ~4–5 mean little
promoter signal whatever the sample) and SnapATAC2/ArchR/Signac all use fixed cuts; nucleosome
signal and fraction of reads in peaks are not used (no peaks exist yet; FRiP is a per-label
property computed after atac-regulome's peak calling).

## Ports

| Port | Kind | Files | Notes |
|---|---|---|---|
| `fragments` | input (root) | `*fragments.tsv.gz` (+ `.tbi`), `per_barcode_metrics.csv` / `singlecell.csv`, `atac_peaks.bed`, `samples.csv`, SPARC sheets | Searched recursively; RNA files of an ARC `outs/` are ignored. |
| `cells` | input (subfolder, optional) | `<sample_id>.h5ad` from sc-prep / sc-qc | Linked: `INPUT_DIR/cells/`; unlinked: any `.h5ad` under the input. |
| `atac` | output | `atac/<sample_id>.h5ad`, `atac/provenance.json` | Below. |
| `qc` | output | `qc/atac_qc.csv`, `qc/qc_summary.csv`, `qc/<sample_id>_tsse.png`, `qc/<sample_id>_fragment_sizes.png`, `qc/summary.json`, `qc/provenance.json` | Below. |

### `atac/<sample_id>.h5ad`

A SnapATAC2 AnnData (HDF5, `.h5ad` layout). Read it with `snapatac2.read(path, backed="r")`;
plain `anndata.read_h5ad` also works after `import hdf5plugin` (SnapATAC2 compresses with Blosc).

| Slot | Content |
|---|---|
| `X` | tile matrix: Tn5 insertions per 500 bp bin, CSR (chrM, chrY excluded) |
| `obsm["fragment_paired"]` | every fragment of every kept nucleus (SnapATAC2's layout; atac-regulome calls peaks from it) |
| `obs_names` | `<sample_id>:<barcode>` |
| `obs` | `n_fragment`, `frac_dup`, `frac_mito`, `tsse`, (`doublet_probability`, `doublet_score` ATAC-only), `sample_id` and the sample's metadata columns (categorical) |
| `var["selected"]` | ATAC-only runs: the bins Scrublet used |
| `uns` | `reference_sequences`, `TSS_profile` (insertions ±2 kb around TSSs, summed over the sample's nuclei), `library_tsse` (peak of its 21 bp running mean over the mean of the outer 100 bp), `frag_size_distr`; `pennsieve`: `source`, `genome`, `genome_from`, `organism`, `cells_from`, `tile_bin_size`, `qc` (thresholds, scrublet), `prep_version`, `sample_id` |

### `qc/atac_qc.csv`

Every evaluated barcode: `sample_id, obs_id, n_fragment, frac_dup, frac_mito, tsse, passed,
reason` (`low_fragments`, `high_fragments`, `low_tsse`, `doublet`, joined with `;`).

### `qc/qc_summary.csv`

One row per sample: `sample_id, n_cells_before, n_cells_after, n_removed_low_fragments,
n_removed_high_fragments, n_removed_low_tsse, n_removed_doublet, median_tsse, median_fragments,
cells_from, genome`.

### Figures, `summary.json`, `provenance.json`

- `<sample_id>_tsse.png`: TSS enrichment against unique fragments for every evaluated nucleus,
  kept and removed, with the cut-offs; right, the sample's insertion profile around TSSs.
- `<sample_id>_fragment_sizes.png`: fragment-length distribution of the kept nuclei: a good
  library shows the nucleosome-free peak below ~150 bp and mono-/di-nucleosome peaks near 200
  and 400 bp.
- `summary.json`: per sample counts, removals, medians, library TSSE, genome and where the cell
  calls came from. `provenance.json` (both folders): version, tools, parameters, inputs, per-sample
  status, timings, peak memory, citations.

## Parameters

| Parameter | Default | What it does |
|---|---|---|
| `SAMPLE_SHEET` | `samples.csv` | Optional sample sheet (`sample_id,path,…`); its other columns go onto the cells. With a sheet only its rows are samples. |
| `METADATA_COLUMNS` | `subject_id,subject_experimental_group,sample_experimental_group,sex,age,species,sample_anatomical_location,laterality` | SPARC sheet columns copied onto the cells (never other clinical fields). |
| `CONDITION_COLUMN` | `auto` | Column copied to `condition` (`auto`: subject, then sample experimental group). |
| `GENOME` | `auto` | `auto`, `GRCh38` or `GRCm39`. |
| `MIN_FRAGMENTS` | `1000` | Minimum unique fragments per nucleus. |
| `MAX_FRAGMENTS` | `100000` | Maximum unique fragments per nucleus. |
| `MIN_TSSE` | `5` | Minimum TSS enrichment. |
| `SEED` | `42` | Scrublet's random state. |
| `THREADS` | `0` | SnapATAC2 / BLAS threads; 0 = the container's CPU quota. |

Bad values exit with code 2 before anything is read.

## Results on the demo data

Discover 465 (human sympathetic ganglia, single-nucleus multiome) is one aggregated fragment file;
`scripts/cut_sympathetic.py` cut libraries -1 and -2 to ~2,000 nuclei each (up to 245 per author
cell type, so all eleven types are represented) and the matching RNA counts. Run with the RNA
side linked (default parameters, `GENOME=auto` resolved from the RNA file's organism, since the
Discover file has no header); outputs in `workflow-single-cell/data/work/sympathetic/atac-prep/`
for atac-regulome.

| Sample | Fragments in the cut | Nuclei (RNA) | Kept | Removed: low fragments / high / low TSSE | Median fragments | Median TSSE | Library TSSE | Median dup. fraction |
|---|---|---|---|---|---|---|---|---|
| sg1 | 26.5 M | 2,011 | 1,996 | 8 / 1 / 6 | 11,428 | 16.8 | 20.2 | 0.33 |
| sg2 | 15.4 M | 2,019 | 2,005 | 6 / 0 / 8 | 6,516 | 15.6 | 20.1 | 0.30 |

These nuclei had already passed the authors' joint RNA+ATAC QC, so almost all pass here; the
fragment-size plots show the nucleosome-free peak and the mono- and di-nucleosome shoulders of a
good library, and `frac_mito` is 0 (nuclei). The same sg2 file **without** the RNA link
(`GENOME=GRCh38`, every barcode with ≥ 1,000 fragments): 2,013 evaluated → 1,929 kept, Scrublet
called 76 doublets (3.8 %), 8 low TSSE.

The RE-JOIN knee layout (Cell Ranger ARC per sample: `atac_fragments.tsv.gz` + `.tbi`,
`atac_peaks.bed`, `filtered_feature_bc_matrix.h5`, `per_barcode_metrics.csv`) is the fixture's
`knee1/outs`, `knee2/outs` shape; the tests run it with and without the RNA link.

## Runtime

`computeTypes: [standard]`, 4 vCPU / 8 GiB, 6 h timeout. Measured in the linux/amd64 image under
emulation on an M-series Mac (Docker `--memory 4g --cpus 4`, `--user 1000`); peak memory = the
process's max RSS and the cgroup's `memory.peak` (which also counts page cache and, with a tmpfs
`/tmp`, the sorted fragments):

| Run | Wall time | Peak RSS | cgroup peak |
|---|---|---|---|
| sg1 + sg2 with RNA (42 M fragments, 4,030 nuclei), `/tmp` on disk | 127 s | 1.47 GB | 2.23 GB |
| the same, read-only root, `/tmp` as tmpfs | 187 s | 1.32 GB | 1.82 GB |
| sg2 ATAC-only, with Scrublet (2,013 nuclei), tmpfs `/tmp` | 94 s | 1.93 GB | 2.51 GB |

**Full-size samples.** Memory no longer follows the sample's size: the fragments are filtered to
the called nuclei and sorted by barcode on disk (GNU sort), imported 100 nuclei at a time, and TSS
enrichment runs on 250-nucleus subsets, because SnapATAC2's own unsorted import and `tsse` hold
every fragment in memory (~130–140 bytes each, measured: 2.1 GB and 3.4 GB for 15 M and 26 M
fragments; a 10,000-nucleus ARC sample with ~20,000 fragments per nucleus would need ~30 GB that
way). What still grows: Scrublet (ATAC-only runs) holds the cells × selected-bins matrix, about
+0.6 GB for these 2,000 nuclei, so roughly 3–5 GB for a 10,000-nucleus sample (an extrapolation);
with the RNA link it does not run. 8 GiB covers the measured 2.2 GB cgroup peak and that estimate.
Disk in `/tmp`: about twice the compressed size of the sample's fragment file (GNU sort's compressed
temporaries, then the sorted copy, gzip -1) plus one 250-nucleus subset. Output: ~16 bytes per kept fragment (412 MB for sg1), so ~3–4 GB per
full ARC sample. Time: ~4.5 µs per fragment under emulation, i.e. ~15–20 min for a 200 M-fragment
sample, less on native amd64. Samples are processed one after another.

Image: 1.54 GB (`python:3.12-slim`, SnapATAC2 and its dependencies, the two GENCODE annotations).

## Running it locally

```sh
python3.12 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
make test PYTHON=.venv/bin/python       # pytest on the synthetic fixture (no Docker)
make run PYTHON=.venv/bin/python        # build the image, run it over data/fixture/linked -> data/output
make validate PYTHON=.venv/bin/python   # app.yml against the manifest schema, CITATION.cff
```

The fixture (`tests/make_fixture.py`, see `sample-input/README.md`) has two Cell Ranger ARC-shaped
samples (the RE-JOIN layout) with planted low-fragment and no-TSS-enrichment nuclei, a bare fragment
file, and the RNA side; the tests check the exact nuclei kept with and without RNA, Cell Ranger's
calls, Scrublet, SPARC metadata (clinical columns stay off), genome detection and parameter errors.

## Limitations

- **Two genomes.** GRCh38 and GRCm39 only. Mouse data aligned to mm10 (Cell Ranger ARC's
  `mm10-2020-A` reference) must be realigned or lifted over; rat is not supported.
- **One fragment file = one sample.** An aggregated (`cellranger-arc aggr`) file is treated as one
  sample whose barcodes carry `-1`, `-2`, … suffixes; split it per library (as
  `scripts/cut_sympathetic.py` does for the demo) so samples match the RNA side.
- **Fixed thresholds.** They do not adapt to a sample; a sample with uniformly low TSSE (poor
  nuclei prep) loses most nuclei, which `qc_summary.csv` shows. Lower `MIN_TSSE` deliberately
  rather than silently.
- **No ambient-DNA or barcode-multiplet correction** beyond Cell Ranger's calling and (ATAC-only)
  Scrublet. AMULET (Thibodeau et al. 2021), which counts loci with more than two reads per nucleus,
  is the ATAC-specific multiplet detector the best-practices guide mentions; it is not run here
  (multiome data get their doublet calls from the RNA side).
- **Disk.** Fragments are sorted on disk under `/tmp`: allow about twice the compressed size of the
  largest sample's fragment file.

## Security

Inputs are read through `processor/safe_inputs.py`, vendored verbatim from
`workflow-single-cell/shared/` (do not edit the copy): only files that resolve inside `INPUT_DIR` or
the platform's links are used, sheets are size-capped (`.xlsx` also by unpacked size, parsed with
defusedxml), RNA `.h5ad` files with external HDF5 links or storage are refused. Fragment files are
streamed through `gzip`/`awk`/`sort` subprocesses without a shell. The container runs as uid 1000
with a read-only root; it writes only to `OUTPUT_DIR` and `/tmp`.

## Licences

Apache-2.0 for this app (`LICENSE`). Bundled software and data and their licences:
`THIRD_PARTY.md` (no GPL Python packages; the Debian base layer contains GPL system tools).

## Credits

- **SnapATAC2** (Kai Zhang, Bing Ren and colleagues; scverse): fragment import, TSS enrichment,
  fragment sizes, tile matrix and Scrublet; the default thresholds come from its documentation.
- **GENCODE** (the GENCODE consortium, EMBL-EBI and partners): the human and mouse gene annotations
  whose TSSs define TSS enrichment.
- **Scrublet** (Samuel Wolock, Romain Lopez, Allon Klein): doublet scores for ATAC-only runs, as
  reimplemented in SnapATAC2.
- **ArchR** (Jeffrey Granja, William Greenleaf, Howard Chang and colleagues): the comparison for the
  default thresholds.
- **anndata** (Isaac Virshup and the scverse team): reading the RNA side.
- **Single-cell best practices** (Heumos et al. and the Single-cell Best Practices Consortium): the
  overall QC approach followed here.
- **Demo data:** Lite Yang and Robert W. Gereau (Washington University in St. Louis), *A single-cell
  multi-omic atlas of human sympathetic ganglia*, SPARC / Pennsieve Discover 465.

## Citations

> Zhang K, Zemke NR, Armand EJ, Ren B. **A fast, scalable and versatile tool for analysis of
> single-cell omics data.** Nat Methods 21(2):217-227 (2024). https://doi.org/10.1038/s41592-023-02139-9,
> PMID 38191932.

> Frankish A, Carbonell-Sala S, Diekhans M, Jungreis I, Loveland JE, Mudge JM, Sisu C, Wright JC,
> et al. (62 authors). **GENCODE: reference annotation for the human and mouse genomes in 2023.**
> Nucleic Acids Res 51(D1):D942-D949 (2023). https://doi.org/10.1093/nar/gkac1071, PMID 36420896.

> Wolock SL, Lopez R, Klein AM. **Scrublet: Computational Identification of Cell Doublets in
> Single-Cell Transcriptomic Data.** Cell Syst 8(4):281-291.e9 (2019).
> https://doi.org/10.1016/j.cels.2018.11.005, PMID 30954476.

> Granja JM, Corces MR, Pierce SE, Bagdatli ST, Choudhry H, Chang HY, Greenleaf WJ. **ArchR is a
> scalable software package for integrative single-cell chromatin accessibility analysis.** Nat
> Genet 53(3):403-411 (2021). https://doi.org/10.1038/s41588-021-00790-6, PMID 33633365.

> Thibodeau A, Eroglu A, McGinnis CS, et al. (14 authors). **AMULET: a novel read count-based method
> for effective multiplet detection from single nucleus ATAC-seq data.** Genome Biol 22(1):252
> (2021). https://doi.org/10.1186/s13059-021-02469-x, PMID 34465366 (mentioned, not used).

> Heumos L, Schaar AC, Lance C, et al. (43 authors incl. the Single-cell Best Practices
> Consortium). **Best practices for single-cell analysis across modalities.** Nat Rev Genet
> 24(8):550-572 (2023). https://doi.org/10.1038/s41576-023-00586-w, PMID 37002403.

> Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. **anndata: Access and store annotated data
> matrices.** J Open Source Softw 9(101):4371 (2024). https://doi.org/10.21105/joss.04371 (not in
> PubMed; checked on Crossref).

Demo data:

> Yang L, Gereau R. **A single-cell multi-omic atlas of human sympathetic ganglia** (Version 1)
> [Data set]. SPARC Portal (2025). https://doi.org/10.26275/oqxx-si9h, CC BY 4.0 (checked on
> DataCite and the Pennsieve Discover API, dataset 465). Probably described by: Yang L, Dourson AJ,
> Tao R, et al. (18 authors), Gereau RW. **A multi-omic atlas of human autonomic and sensory ganglia
> implicates cell types in peripheral neuropathies.** bioRxiv (2025).
> https://doi.org/10.1101/2025.09.18.677119, PMID 41000833 (same authors and tissue; the dataset
> record does not link it).

> Yang L, Meriau P, Cavalli V, Gereau RW. **Nuclear extraction using opti-prep gradient
> centrifugation for 10X Genomics Multiome Assay v1.** protocols.io (2024).
> https://doi.org/10.17504/protocols.io.5jyl8q4y8l2w/v1 (the dataset's protocol; checked on Crossref).

Every DOI was checked on Crossref (DataCite for the dataset) and every PMID on PubMed (E-utilities
`esummary`: title, authors, journal, volume, pages); full author lists are in `CITATION.cff`,
generated from the Crossref records.
