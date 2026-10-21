# scATIMs 00–07 workflow contract

## Stage map

| Stage | Required input | Primary output | Implementation status |
|---|---|---|---|
| 00 QC | paired R1/R2 FASTQ | `00.QC/<sample>/*.QC.R[12].fq.gz` and fastp reports | implemented; real-run validation pending |
| 01 STARsolo | stage-00 FASTQ, STAR index, chemistry whitelist | `01.starSolo/Solo.out/Gene/{raw,filtered}_gz` | implemented; real-run validation pending |
| 01 identity | same-run FASTQ, STAR index and whitelist | audited `01.starSolo_bam/CR_CB_mapping.tsv` | implemented; real-run validation pending |
| 02 identity embedding | paired FASTQ and chemistry definition | `02.umiTools/<sample>.raw.reads_barcoded.fq.gz` | implemented; real-run validation pending |
| 03 GRCh38 filter | stage-02 candidate reads | `03.hg38/*_hg38_unmap.{bam,fq}` | implemented; real-run validation pending |
| 04 T2T filter | stage-03 unmapped reads | `04.t2t/*_t2t_unmap.bam` | implemented; real-run validation pending |
| 04b vector filter | stage-04 unmapped BAM and UniVec | `04b.vector/02.bam/*.vector_unmapped.bam` | implemented, optional; real-run validation pending |
| 05 PathSeq | stage-04 or 04b BAM and PathSeq resources | `05.pathseq_v10/*.pathseq.{bam,txt}` | implemented; real-run validation pending |
| 06 scMic | PathSeq BAM, whitelist, CR/CB map, taxonomy hierarchy | estimated-UMI assignments and correction metrics | implemented and tested |
| 07 ATIMs | stage-06 assignments, whitelist and CR/CB map | Statistics/CellEnriched/EmptyEnriched/Parameters | implemented and tested |

“Real-run validation pending” is deliberate and must remain visible until the
parameterized 00–05 implementation completes a controlled raw-FASTQ run. Unit
and dry-run tests do not by themselves establish full production validation.

## Human workflow scope

The inspected paper workflow is human scRNA-seq and uses sequential GRCh38 and
T2T-CHM13 host filtering. Mouse code exists in historical directories but is
not equivalent or validated and is outside the current release scope.

The actual runtime versions observed for these stages were fastp 0.23.2,
STAR 2.7.10a, BWA 0.7.17-r1188, samtools 1.5 and GATK 4.6.0.0. The complete
runtime snapshot remains the provenance authority.

Library chemistry controls barcode and UMI coordinates and the STARsolo
whitelist. Those settings must be explicit per run; they must not be inferred
from sample names.

## Identity continuity

STARsolo must emit both `CR` and `CB` tags. `scatims cr-cb` streams the tagged
BAM and writes deterministic mapping, audit and conflict files. To reproduce
the original `run_CR2CB_multi_samples_del.sh` behavior, every unique valid
`(CR, CB)` pair is retained in sorted order, including a raw barcode mapping to
multiple corrected barcodes. The conflicts file records those cases; strict
failure is available only as an explicit validation option.
Stage 02 embeds the raw barcode and UMI in the candidate-read identity so stage
06 can recover the corrected cell identity after host and vector filtering.

## Attribution

Stages 00–05 wrap or coordinate fastp, STARsolo, BWA, samtools and GATK
PathSeq. Their inclusion defines the scATIMs software workflow, not ownership
of those third-party algorithms. scMic and ATIMs are the original method core.
# Execution control

The complete workflow is dispatched by `scatims workflow --config CONFIG.yaml`.
Set the inclusive range in the configuration itself:

```yaml
pipeline:
  start_step: "00"
  end_step: "07"
```

Valid step identifiers are `00`, `01`, `02`, `03`, `04`, `04b`, `05`, `06`,
and `07`. This release deliberately does not offer command-line range
overrides, so the configuration archived with a run is the authoritative record
of the requested interval. Stages 05–07 use fixed producer/consumer folders:
`05.pathseq_v10`, `06.scMic_6_V10_EM_Dis`, and
`07.ATIMs_9_V10_EM_Dis`.
