# Changelog

## 0.2.0rc1 — 2026-09-12

- Defined scATIMs as the complete human scRNA-seq stage 00–07 workflow.
- Added parameterized fastp, STARsolo, identity embedding, GRCh38, T2T,
  optional UniVec and PathSeq execution stages.
- Added deterministic audited CR-to-CB mapping with conflict gating.
- Added stage ranges, dry-run, force-rerun controls and a complete configuration template.
- Retained scMic and ATIMs as the original stage-06/07 method core.
- Recorded exact checksums and observed composition of the paper V7 resource snapshot.

## 0.1.0rc1 — 2026-09-11

- Established the private scATIMs release candidate.
- Defined scMic as step 1 and ATIMs as step 2 under one package name.
- Preserved the frozen paper algorithm while parameterizing sample-relative input paths.
- Locked the paper configuration, including `min_evidence_score: 0.0`.
- Added environment, documentation and test scaffolding.
