# Algorithm and workflow overview

## Workflow boundary

The scATIMs software workflow spans stages 00–07: FASTQ QC (00), STARsolo
expression/barcode processing (01), barcode/UMI identity embedding (02),
GRCh38 filtering (03), T2T-CHM13 filtering (04), optional UniVec filtering
(04b), PathSeq microbial classification (05), scMic assignment (06), and
ATIMs decontamination (07).

Stages 00–05 are reproducible workflow components built around established
third-party tools. Stages 06–07 constitute the scATIMs core method. Calling the
complete workflow scATIMs does not imply authorship of fastp, STARsolo, BWA,
samtools, GATK PathSeq or their underlying algorithms.

## Stage 06: scMic

scMic reads microbial alignments carrying the PathSeq `YP` taxonomy tag and
uses the STARsolo-derived CR-to-CB mapping to correct raw read barcodes. Reads
are separated into cell and empty-droplet streams using the filtered barcode
whitelist.

Within each corrected barcode and physical alignment pattern, UMIs are
collapsed with the directional abundance rule and Hamming distance one. For
each taxonomic rank, unambiguous UMI-pattern evidence is used to construct
anchors. Ambiguous candidates are allocated in proportion to those anchors; if
no candidate has anchor support, the signal is assigned to `unknown`.

Anchor construction is performed independently for the cell and empty-droplet
streams. The historical internal function name `calculate_global_anchors` is
retained in release candidate 0.1.0rc1 for provenance, but “global anchoring” is
not the public method description.

## Stage 07: ATIMs

ATIMs aggregates scMic scores by barcode and taxon. For each taxon and rank it
computes:

- a one-sided Fisher exact test on positive-droplet counts;
- a one-sided Mann–Whitney U test on full droplet score distributions,
  including zeros;
- log2 fold change and enrichment ratio;
- a cell-side dispersion test based on the positive cell score distribution.

P values are corrected within each tested taxonomic-rank result table using
the Benjamini–Hochberg method. CellEnriched requires all five cell-side
conditions; EmptyEnriched requires Fisher, Mann–Whitney, enrichment-ratio and
log2FC conditions and does not use an empty-side dispersion condition.

The implementation reports each taxonomic rank separately. A genus-level
classification does not imply the same result at species level because each
rank has its own candidate set and multiple-testing correction.
