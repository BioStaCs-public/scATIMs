# Real-FASTQ release validation

## Controlled equivalence test

Commit `c446950598d9b6a61f4d5388c52f951ff7db7d31` was run from stages 00
through 07 on one representative real 10x 3' v2 FASTQ library, using the same
tool and reference-resource configuration as an archived legacy workflow
result. The candidate output was written to a new isolated location; neither
the raw reads nor the legacy result was modified.

The comparison evaluated semantic output contracts rather than compressed-file
bytes: fastp counts; STARsolo barcode content; CR-to-CB mapping table; sorted
read-name sets after GRCh38, T2T-CHM13 and UniVec filtering; PathSeq score and
`(read_name, YP)` assignments; and all stage-06/07 TSV outputs.

## Result

All compared computational and biological outputs were identical. The only
intentional schema difference was the removal of the retired
`min_evidence_score: 0.0` key from the stage-07 parameter JSON files; all
shared parameter values and all stage-07 statistical outputs were identical.

This test clears the controlled real-FASTQ equivalence gate for the identified
commit. It is a release reproducibility check, not a claim of validation across
all chemistries, organisms, sequencing platforms or reference-database
versions.
