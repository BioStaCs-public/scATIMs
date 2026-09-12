# scATIMs 00–07 workflow contract

## Stage map

| Stage | Required input | Primary output | Implementation status |
|---|---|---|---|
| 00 QC | paired R1/R2 FASTQ | `00.QC/<sample>/*.QC.R[12].fq.gz` and fastp reports | integration pending |
| 01 STARsolo | stage-00 FASTQ, STAR index, chemistry whitelist | `01.starSolo/Solo.out/Gene/{raw,filtered}_gz` | integration pending |
| 01 identity | same-run FASTQ, STAR index and whitelist | audited `01.starSolo_bam/CR_CB_mapping.tsv` | extractor implemented; STAR invocation pending |
| 02 identity embedding | paired FASTQ and chemistry definition | `02.umiTools/<sample>.raw.reads_barcoded.fq.gz` | integration pending |
| 03 GRCh38 filter | stage-02 candidate reads | `03.hg38/*_hg38_unmap.{bam,fq}` | integration pending |
| 04 T2T filter | stage-03 unmapped reads | `04.t2t/*_t2t_unmap.bam` | integration pending |
| 04b vector filter | stage-04 unmapped BAM and UniVec | `04b.vector/02.bam/*.vector_unmapped.bam` | integration pending; optional by design |
| 05 PathSeq | stage-04 or 04b BAM and PathSeq resources | `05.pathseq_v10/*.pathseq.{bam,txt}` | parameterized wrapper present |
| 06 scMic | PathSeq BAM, whitelist, CR/CB map, taxonomy hierarchy | estimated-UMI assignments and correction metrics | implemented and tested |
| 07 ATIMs | stage-06 assignments, whitelist and CR/CB map | Statistics/CellEnriched/EmptyEnriched/Parameters | implemented and tested |

“Integration pending” is deliberate and must remain visible until the
parameterized implementation is tested. The repository must not claim a
complete raw-FASTQ run merely because the stage boundary is documented.

## Human workflow scope

The inspected paper workflow is human scRNA-seq and uses sequential GRCh38 and
T2T-CHM13 host filtering. Mouse code exists in historical directories but is
not equivalent or validated and is outside the current release scope.

Library chemistry controls barcode and UMI coordinates and the STARsolo
whitelist. Those settings must be explicit per run; they must not be inferred
from sample names.

## Identity continuity

STARsolo must emit both `CR` and `CB` tags. `scatims cr-cb` streams the tagged
BAM and writes deterministic mapping, audit and conflict files. A raw barcode
mapping to multiple corrected barcodes is a contract violation by default.
Stage 02 embeds the raw barcode and UMI in the candidate-read identity so stage
06 can recover the corrected cell identity after host and vector filtering.

## Attribution

Stages 00–05 wrap or coordinate fastp, STARsolo, BWA, samtools and GATK
PathSeq. Their inclusion defines the scATIMs software workflow, not ownership
of those third-party algorithms. scMic and ATIMs are the original method core.
