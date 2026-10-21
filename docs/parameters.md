# Parameters

## Paper-locked scATIMs configuration

| Parameter | Paper value |
|---|---:|
| Taxonomic levels | species, genus, family, order, class, phylum |
| Fisher FDR threshold | 0.05 |
| Mann–Whitney FDR threshold | 0.05 |
| Cell dispersion FDR threshold | 0.05 |
| Enrichment ratio threshold | 1.0, strict `>` comparison |
| log2FC threshold | 1.0, strict `>` comparison |

Use `configs/paper_parameters.yaml` for paper reproduction. There is no
configurable ATIMs minimum-evidence-score threshold in this release. Before
statistical testing, ATIMs retains assignments with `score > 0`, excludes
`unknown` taxa, and then applies the tabled rank-level criteria.

## Resource parameters

The paper batch runner used 32 sample workers. This is a resource choice, not a
scientific threshold; reduce it for smaller machines. PathSeq Spark and JVM
memory defaults in `workflows/run_pathseq.sh` are also resource settings and
can be overridden by CLI options.

## PathSeq parameters recorded for the paper workflow

- `--min-clipped-read-length 60`
- `--is-host-aligned false`
- `--filter-duplicates false`
- `--min-score-identity 0.7`

PathSeq is stage 05 of the complete scATIMs workflow. These values are workflow
parameters, not parameters of the scMic or ATIMs statistical algorithms.
