# Release provenance

## Frozen runtime sources

| Runtime file | Candidate destination | Frozen SHA-256 |
|---|---|---|
| `scMic_zyx6.py` | `src/scatims/scmic.py` | `78e0825f29c4e7273e943d35135422fde37f58f1c3aa74b8565ec467d58808a6` |
| `ATIMs_zyx9.py` | `src/scatims/atims.py` | `4fac7bcae4d27bb18df20bf0fad95729f4286c4e90a082bca15bf04793315ba0` |
| `run_pathseq.sh` | `workflows/run_pathseq.sh` | `39ac9e99cee8002e7690590d46d80668892e831901302892c77831ed75b32aa3` |
| `run_scMic_ATIMs_pipeline_zyx.sh` | replaced by `scatims run --config ...` | `bbc78e54017be745e6ebfccd033a0be8917a42faae150ff6a018628fca537c29` |

## Candidate interface changes

- Unified the public project name as scATIMs with scMic and ATIMs subcommands.
- Converted fixed sample-relative input locations to configurable templates.
- Removed server absolute paths and local sample-list defaults.
- Added a YAML-driven two-step runner and paper parameter file.
- Changed the candidate standalone ATIMs minimum-evidence default from the
  historical 1.0 to the actual paper value 0.0.
- Required Parameters JSON to exist before ATIMs treats a rank as complete.
- Normalized the candidate scMic source from CRLF to LF.
- Updated public terminology from global anchoring to compartment-specific
  anchor construction.

The Fisher, Mann–Whitney, fold-change, enrichment-ratio, dispersion,
Benjamini–Hochberg and anchor-weighted assignment calculations were not
intentionally changed. Equivalence beyond the included synthetic tests still
requires comparison against frozen real-run outputs before public release.
