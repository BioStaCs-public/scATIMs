# scATIMs

> Private release candidate. Do not redistribute or make this repository public
> until authorship, licensing and reference-database release decisions are approved.

scATIMs is a two-step workflow for assigning microbial signals in single-cell
RNA-sequencing data and identifying cell- or empty-droplet-enriched taxa:

1. **scMic** corrects cell barcodes, collapses UMIs directionally and performs
   anchor-weighted allocation of ambiguous microbial alignments separately in
   cell and empty-droplet data streams.
2. **ATIMs** compares cell and empty droplets with Fisher exact,
   Mann–Whitney U and cell-side dispersion tests, applies Benjamini–Hochberg
   correction and reports cell- and empty-enriched taxa.

The package consumes PathSeq-tagged BAM files plus STARsolo barcode artifacts.
PathSeq and host-removal are prerequisites, not additional scATIMs algorithm steps.

## Release status

- Candidate version: `0.1.0rc1`
- Paper runtime source is frozen separately and is not edited by this repository.
- Server-specific paths have been replaced with CLI/configuration parameters.
- The paper parameter configuration is in `configs/paper_parameters.yaml`.
- No real sample data, reference database, GATK binary or patient information is included.
- Git has not yet been initialized for this candidate.

## Installation

Create the tested Python 3.11 environment:

```bash
conda env create -f environment/environment.yml
conda activate scatims
python -m pip install -e .
scatims --version
```

For exact internal provenance, the task archive retains the original
`scMic_v2` explicit conda specification. The curated environment here contains
only dependencies needed to run the release candidate.

## Run both steps

Copy the paper configuration and replace its three placeholder paths:

```bash
cp configs/paper_parameters.yaml config.local.yaml
scatims run --config config.local.yaml
```

`config.local.yaml` is ignored by Git. The paper configuration locks all six
taxonomic ranks and `min_evidence_score: 0.0`.

Individual steps remain available:

```bash
scatims scmic --help
scatims atims --help
```

## Required input layout

By default, each sample under `project.root_dir` is expected to contain:

```text
<sample_id>/
├── 05.pathseq_v10/<sample_id>_rna_output.pathseq.bam
├── 01.starSolo/Solo.out/Gene/filtered_gz/barcodes.tsv.gz
└── 01.starSolo_bam/CR_CB_mapping.tsv
```

All three relative paths are configurable. Detailed schemas and output files
are documented in `docs/input_output.md`.

## Reference resources

The taxonomy hierarchy is supplied with `scmic.taxonomy_hierarchy`. The custom
V7 database used for the paper is not redistributed. See `resources/README.md`
for the required PathSeq and taxonomy components.

## Tests

The tests use synthetic data only and include a complete scMic-to-ATIMs run:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

The integration test constructs a PathSeq-tagged BAM, barcode whitelist,
CR-to-CB mapping and taxonomy table, then verifies Statistics, CellEnriched,
EmptyEnriched and Parameters outputs.

## Documentation

- `docs/algorithm.md`: algorithm stages and terminology.
- `docs/input_output.md`: file contracts and outputs.
- `docs/parameters.md`: paper-locked and resource parameters.
- `docs/release_provenance.md`: mapping from frozen runtime scripts to this candidate.
