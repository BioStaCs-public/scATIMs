# Input and output contract

## Project and sample list

`project.root_dir` contains one directory per sample. `project.sample_list` is
a UTF-8 text file with one sample identifier per non-comment line. Blank lines
and lines beginning with `#` are ignored.

For a complete run, `project.input_dir` contains one directory per sample with
paired FASTQ files. The default templates are
`<sample>/<sample>.R1.fq.gz` and `<sample>/<sample>.R2.fq.gz`; both templates
are configurable. Sample identifiers are restricted to letters, digits,
period, underscore and hyphen, and duplicates are rejected.

## Stages 00–05 outputs

- `00.QC/<sample>/`: QC FASTQ pair and fastp JSON/HTML reports.
- `01.starSolo/`: expression matrices plus gzip copies under
  `Solo.out/Gene/{raw,filtered}_gz`.
- `01.starSolo_bam/`: CR/CB-tagged STARsolo BAM,
  `CR_CB_mapping.tsv`, `CR_CB_mapping_audit.tsv` and conflict table.
- `02.umiTools/`: barcode/UMI-embedded candidate-read FASTQ, unique raw
  barcode table and per-barcode read counts. The directory name is retained
  for paper provenance; the release implementation performs streaming
  identity embedding internally and does not call the `umi_tools` executable.
- `03.hg38/`: GRCh38 aligned, unmapped and flagstat artifacts.
- `04.t2t/`: T2T-CHM13 aligned, unmapped and flagstat artifacts.
- `04b.vector/`: optional UniVec aligned and vector-unmapped artifacts.
- `05.pathseq_v10/`: PathSeq BAM, score table and score/filter metrics.

Stage completion requires the documented terminal outputs to exist and be
non-empty. `--force` reruns selected stages; otherwise completed stages resume
without overwriting their results.

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
