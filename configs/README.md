# Configuration files

`paper_parameters.yaml` contains the scientific thresholds and six-rank setup
used for the paper. Copy it to a Git-ignored local file and replace the three
placeholder paths; do not edit the locked template for each dataset.

`full_workflow.example.yaml` is the complete stage 00–07 template. Copy it to
`config.local.yaml`, replace every `/path/to/...` value and select the explicit
10x chemistry before running `scatims workflow`. The example contains no
server path and is not itself a paper execution record.

Resource choices such as `workers` may be reduced without changing the
scientific decision rules. Any change to Fisher, Mann–Whitney, dispersion,
enrichment-ratio or log2FC thresholds should be recorded as a new analysis
configuration. ATIMs has no configurable minimum-evidence-score threshold in
this release.
