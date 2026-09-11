# Configuration files

`paper_parameters.yaml` contains the scientific thresholds and six-rank setup
used for the paper. Copy it to a Git-ignored local file and replace the three
placeholder paths; do not edit the locked template for each dataset.

Resource choices such as `workers` may be reduced without changing the
scientific decision rules. Any change to Fisher, Mann–Whitney, dispersion,
enrichment-ratio, log2FC or minimum-evidence thresholds should be recorded as a
new analysis configuration.
