# scATIMs

> Private release candidate. Do not redistribute or make this repository public
> until authorship, licensing and reference-database release decisions are approved.

scATIMs is an end-to-end workflow from raw single-cell RNA-sequencing FASTQ
files to statistically decontaminated microbial signals. It comprises stages
00–07:

| Stage | Role |
|---|---|
| 00 | FASTQ quality control |
| 01 | STARsolo host-expression quantification and CR-to-CB identity recovery |
| 02 | Barcode/UMI identity embedding in candidate reads |
| 03 | GRCh38 host filtering |
| 04 | T2T-CHM13 host filtering |
| 04b | Optional UniVec/vector filtering |
| 05 | GATK PathSeq microbial classification |
| 06 | scMic barcode correction, directional UMI collapse and anchor-weighted assignment |
| 07 | ATIMs cell-versus-empty-droplet statistical decontamination |

Stages 00–05 orchestrate established third-party tools and preserve microbial
candidate-read identity. The original computational method is concentrated in
stage 06 (scMic) and stage 07 (ATIMs); this distinction must be retained in
method descriptions and software attribution.

## Release status

- Candidate version: `0.1.0rc1`
- Paper runtime source is frozen separately and is not edited by this repository.
- Server-specific paths have been replaced with CLI/configuration parameters.
- The paper parameter configuration is in `configs/paper_parameters.yaml`.
- No real sample data, reference database, GATK binary or patient information is included.
- The current private repository contains the validated 06–07 core and the
  parameterized stage-05 PathSeq wrapper. Integration of stages 00–04b into
  the unified runner is in progress and must not yet be described as complete.

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

## Run the validated 06–07 core

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
- `docs/workflow_00_07.md`: complete stage boundary, contracts and implementation status.
