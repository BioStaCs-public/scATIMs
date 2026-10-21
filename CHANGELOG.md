# Changelog

## 0.2.0rc2 — 2026-10-10

- Migrated the finalized 05–07 release contract: `05.pathseq_v10`,
  `06.scMic_6_V10_EM_Dis`, and `07.ATIMs_9_V10_EM_Dis` are fixed,
  documented output locations.
- Made YAML the sole control point for the inclusive stage range through
  `pipeline.start_step` and `pipeline.end_step`.
- Removed the retired ATIMs per-droplet minimum-evidence-score option from
  code, configuration, CLI and output parameter JSON. ATIMs now uses the
  pre-existing positive-score input filter (`score > 0`) as its only
  per-droplet inclusion rule.
- Corrected CR-to-CB conflict handling to reproduce the historical
  `sort -u` table: every unique valid `(CR_Raw, CB_Corrected)` pair is emitted
  and conflicts are audited without stopping the normal workflow. The earlier
  candidate incorrectly stopped by default and, in its opt-in non-strict path,
  selected only one corrected barcode.
- Updated defaults to the finalized output locations and bumped the private
  release candidate to `0.2.0rc2`.
- Completed the release hygiene audit documented in
  `docs/release_audit_20261010.md`.

## 0.2.0rc1 — 2026-09-12

- Defined scATIMs as the complete human scRNA-seq stage 00–07 workflow.
- Added parameterized fastp, STARsolo, identity embedding, GRCh38, T2T,
  optional UniVec and PathSeq execution stages.
- Added deterministic audited CR-to-CB mapping infrastructure. In this
  candidate it still stopped by default on conflicting raw barcodes and did
  not yet reproduce the historical all-pairs `sort -u` behavior; that
  compatibility correction is part of `0.2.0rc2` above.
- Added stage ranges, dry-run, force-rerun controls and a complete configuration template.
- Retained scMic and ATIMs as the original stage-06/07 method core.
- Recorded exact checksums and observed composition of the paper V7 resource snapshot.

## 0.1.0rc1 — 2026-09-11

- Established the private scATIMs release candidate.
- Defined scMic as step 1 and ATIMs as step 2 under one package name.
- Preserved the frozen paper algorithm while parameterizing sample-relative input paths.
- Locked the paper configuration used by the original private candidate.
- Added environment, documentation and test scaffolding.
