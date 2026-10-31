# sample-input

atac-prep reads fragment files. The contract is in `workflow-single-cell/docs/conventions.md` ("ATAC").

## Synthetic fixture (no download)

```sh
make fixture PYTHON=.venv/bin/python   # data/fixture/{genomes,linked,flat}/
make run PYTHON=.venv/bin/python       # Docker over data/fixture/linked/ -> data/output/{atac,qc}/
```

`knee1/outs/` and `knee2/outs/` are Cell Ranger ARC per-sample layouts (`atac_fragments.tsv.gz`
with ARC's header naming the GRCh38 reference, `per_barcode_metrics.csv`, `atac_peaks.bed`,
`filtered_feature_bc_matrix.h5`), the shape of the RE-JOIN knee multiome data. Each has 70 called
cells: 60 good, 5 with ~400 fragments and 5 without TSS enrichment, plus 40 uncalled barcodes.
`plain/S3_fragments.tsv.gz` (+ `.tbi`) is a bare fragment file with 220 cells. `cells/` holds the
RNA side for knee1/knee2 (50 good cells + the 10 bad ones, as if RNA QC kept them). The fixture's
genome (`genomes/GRCh38`: two short chromosomes, 60 genes) replaces the baked annotation.

## Real data

**Cell Ranger ARC, one folder per sample** (RE-JOIN, PRECISION multiome): give the folders as they
are, e.g. `<sample>/outs/{atac_fragments.tsv.gz,atac_fragments.tsv.gz.tbi,per_barcode_metrics.csv,
atac_peaks.bed,filtered_feature_bc_matrix.h5}`. sample_id is the folder above `outs/` (as sc-prep
names the RNA half), or a `samples.csv` (`sample_id,path,...`) names them. Link sc-qc's `cells`
output to the `cells` input to keep only the nuclei that passed RNA QC.

**Demo** (Pennsieve Discover 465, human sympathetic ganglia, CC BY 4.0): one aggregated 17.7 GB
fragment file for 8 libraries. `scripts/cut_sympathetic.py` cuts it to libraries -1 and -2 (as
samples `sg1`, `sg2`), ~2,000 cells each (up to 245 per author cell type), and writes the matching
RNA counts and the authors' cell types as `labels.csv`:

```sh
python scripts/cut_sympathetic.py ../workflow-single-cell/data/sympathetic-multiome/primary \
    ../workflow-single-cell/data/work/sympathetic/cut
docker run --rm --platform linux/amd64 --memory 4g --network none --read-only --tmpfs /tmp \
  -e INPUT_DIR=/data/input -e OUTPUT_DIR=/data/output \
  -v .../cut/fragments:/data/input:ro -v .../cut/cells:/data/input/cells:ro \
  -v .../work/sympathetic/atac-prep:/data/output atac-prep:dev
```
