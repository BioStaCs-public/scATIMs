# Input and output contract

## Project and sample list

`project.root_dir` contains one directory per sample. `project.sample_list` is
a UTF-8 text file with one sample identifier per non-comment line. Blank lines
and lines beginning with `#` are ignored.

## scMic inputs

### PathSeq BAM

Default relative template:

`05.pathseq_v10/{sample_id}_rna_output.pathseq.bam`

Requirements:

- readable BAM with reference names and alignment coordinates;
- microbial reads carry a PathSeq `YP` tag containing one or more comma-separated tax IDs;
- query names end in `_<raw_barcode>_<raw_UMI>`.

### Barcode whitelist

Default relative path:

`01.starSolo/Solo.out/Gene/filtered_gz/barcodes.tsv.gz`

One corrected cell barcode per line, plain text or gzip.

### CR-to-CB mapping

Default relative path:

`01.starSolo_bam/CR_CB_mapping.tsv`

Tab-separated table with `CR_Raw` and `CB_Corrected`. If the header names are
absent, scMic/ATIMs fall back to the first and second columns respectively.

### Taxonomy hierarchy

Tab-separated table containing `tax_id` and:

- `species_scientific_name`
- `genus_scientific_name`
- `family_scientific_name`
- `order_scientific_name`
- `class_scientific_name`
- `phylum_scientific_name`

## scMic outputs

For each sample and cell/empty stream:

- `<sample>_<cell|empty>_read_table.tsv`
- `<sample>_<cell|empty>_UMI_assignment.tsv`

Additional sample outputs:

- `<sample>_cell_resolution_summary.tsv`
- `<sample>_correction_metrics.json`

The assignment score is a fractional estimated-UMI quantity and must not be
described as an integer UMI count.

## ATIMs outputs

For each sample and requested taxonomic rank:

- `<sample>_ATIMs_<rank>_Statistics.tsv`
- `<sample>_ATIMs_<rank>_CellEnriched.tsv`
- `<sample>_ATIMs_<rank>_EmptyEnriched.tsv`
- `<sample>_ATIMs_<rank>_Parameters.json`

The Parameters JSON records the thresholds used for that result. Existing
outputs are skipped only when all four files are present and Statistics is non-empty.
