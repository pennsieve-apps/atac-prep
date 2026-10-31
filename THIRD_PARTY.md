# Third-party software and data

This app is Apache-2.0 (see `LICENSE`). The image contains the following third-party components,
installed from PyPI at the versions pinned in `requirements.txt`, and reference data downloaded at
build time; their licences apply to the image. No GPL Python packages (the Debian base layer has GPL system
tools): SnapATAC2's optional dependencies igraph (GPL), plotly/kaleido and MACS3 are deliberately not installed (atac-prep does
not cluster, plot interactively or call peaks).

| Component | What | Licence | Source |
|---|---|---|---|
| SnapATAC2 2.10.0 | `import_fragments`, `metrics.tsse`, `metrics.frag_size_distr`, `add_tile_matrix`, `select_features`, `scrublet`; its Rust core (anndata-rs, precellar) | MIT | https://github.com/scverse/SnapATAC2 |
| anndata 0.13.4, h5py 3.16.0 | Reading the RNA `.h5ad` (`cells`) | BSD-3-Clause | https://github.com/scverse/anndata |
| polars 1.44.2, polars-runtime-32, pyarrow 26.0.0 | SnapATAC2's data frames | MIT; Apache-2.0 | https://pola.rs |
| NumPy 2.5.4, SciPy 1.18.1, pandas 2.3.3, scikit-learn 1.9.1 | Arrays, QC tables; Scrublet's PCA / neighbours | BSD-3-Clause | numpy.org, scipy.org, pandas.pydata.org, scikit-learn.org |
| Matplotlib 3.11.2, Pillow, fonttools, contourpy, kiwisolver, cycler, pyparsing | PNG figures | Matplotlib licence (PSF-based); MIT-CMU; MIT; BSD | https://matplotlib.org |
| openpyxl 3.1.5, et_xmlfile, defusedxml 0.7.1 | SPARC `samples.xlsx` / `subjects.xlsx` (parsed through defusedxml) | MIT; MIT; PSF | https://openpyxl.readthedocs.io |
| rustworkx 0.18.1, pooch 1.9.0, scverse-misc 0.1.8, typeguard 4.6.0, natsort 8.4.0, tqdm 4.70.1, pyfaidx 0.8.2 | Imported by SnapATAC2 | Apache-2.0; BSD-3-Clause; BSD-3-Clause; MIT; MIT; MPL-2.0 + MIT; BSD-3-Clause | PyPI |
| zarr, numcodecs, donfig, google-crc32c, legacy-api-wrap, array-api-compat, pydantic(-core, -settings), requests, urllib3, certifi, idna, charset-normalizer, platformdirs, PyYAML, python-dotenv, msgspec, narwhals, cloudpickle, joblib, threadpoolctl, packaging, python-dateutil, pytz, tzdata, six, typing_extensions, typing-inspection, annotated-types, session-info2 | Dependencies of the above (versions in `requirements.txt`) | MIT / BSD / Apache-2.0 / MPL-2.0 (certifi) / PSF | PyPI |
| python:3.12.16-slim-bookworm | Base image (Debian 12 + CPython 3.12) | PSF; Debian packages under their own free licences | https://hub.docker.com/_/python |

## Data in the image

| File | What | Licence / terms | Source, checksum |
|---|---|---|---|
| `/opt/genomes/GRCh38/annotation.gff3.gz` | GENCODE release 41 basic annotation (GRCh38), TSSs for TSS enrichment | GENCODE data are freely available without restriction (EMBL-EBI terms of use) | ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41, sha256 b82a655b… (the hash SnapATAC2 2.10.0 pins) |
| `/opt/genomes/GRCm39/annotation.gff3.gz` | GENCODE release M30 basic annotation (GRCm39) | as above | …/Gencode_mouse/release_M30, sha256 6f433e26… |
| `/opt/genomes/*/chrom_sizes.tsv` | Primary chromosome lengths | from SnapATAC2's genome tables (MIT) | written at build time |

## Not bundled

Demo data (not in the image or the repository): Yang L, Gereau RW, *A single-cell multi-omic atlas of
human sympathetic ganglia*, SPARC / Pennsieve Discover dataset 465 v1, doi:10.26275/oqxx-si9h
(CC BY 4.0).
