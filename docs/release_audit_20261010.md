# Release hygiene audit — 2026-10-21

Scope: private release candidate `0.2.0rc2`, before any commit or remote push.

## Sensitive-information check

The tracked worktree was scanned for credential-like assignments, known local
credential strings, private IP addresses, and machine-specific data roots.
No credentials, tokens, passwords, email authorization codes, private IPs, or
laboratory data-root paths were found in tracked source, configuration, or
documentation.  Example configuration values use `/path/to/...` placeholders.

The repository intentionally retains its private Git remote URL and its
private-release notice; neither is a credential.

## Reproducibility check

- Python package metadata pins the supported interpreter to Python 3.11 and
  declares all runtime Python dependencies.
- `configs/full_workflow.example.yaml` externalizes tool and reference paths.
- The inclusive workflow range is declared in YAML as `pipeline.start_step`
  and `pipeline.end_step`; no command-line range can silently supersede it.
- The release fixes the producer/consumer directory contract for stages 05–07:
  `05.pathseq_v10`, `06.scMic_6_V10_EM_Dis`, and
  `07.ATIMs_9_V10_EM_Dis`.
- Reference databases and real sequencing data are intentionally excluded;
  their expected components and provenance are documented under `resources/`
  and `docs/`.
- The unit suite uses synthetic inputs only and covers CR-to-CB mapping,
  preprocessing command construction, scMic-to-ATIMs integration, fixed output
  contracts, the all-pairs CR-to-CB conflict behavior, and the absence of the
  retired min-score CLI/configuration field.

## Verification evidence — 2026-10-21

The candidate source was copied to an isolated temporary directory and tested
with the declared paper-runtime Python environment:

- Python 3.11.5; numpy 1.26.4; pandas 2.2.3; scipy 1.14.1;
  statsmodels 0.14.3; pysam 0.24.0; PyYAML 6.0.2.
- `PYTHONPATH=src python -m unittest discover -s tests -v`: 15/15 passed.
- `python -m pip wheel --no-deps --no-build-isolation .`: succeeded; generated
  `scATIMs-0.2.0rc2-py3-none-any.whl` (SHA-256
  `d73d048d7012a6e0b48c518592a7c992494452179934d5ca1347967946f13cf9`).
- `git diff --check` passed. A tracked-file scan found no credentials, private
  IPs, laboratory data roots, email addresses or authorization strings. The
  only lexical matches for credential terms occur in this audit's negative
  statement.

The historical N08 controller directory is not used as release-validation
evidence: it did not record an immutable candidate commit or source digest, and
the current server-side trial copy differs from this candidate. It must not be
described as an end-to-end validation of the public release without rerunning
or comparing outputs with an identified release commit.

## Follow-up before publication

Run one controlled, real-FASTQ equivalence test against the frozen paper
runtime before describing any stage as production-validated. Do not add
reference databases, raw reads, or credentials to this repository.

Before public release, obtain author-team approval for citation metadata,
authorship, licensing and reference-database distribution; see
`docs/CITATION_REQUIRED_BEFORE_PUBLIC.md`.
