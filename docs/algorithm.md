# Algorithm overview

## Step 1: scMic

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

## Step 2: ATIMs

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
