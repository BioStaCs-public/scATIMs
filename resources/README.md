# Reference resources

scATIMs does not redistribute human references, the GATK JAR, PathSeq indices
or the custom V7 microbial database in this repository.

The scMic step requires a taxonomy hierarchy table with at least `tax_id` and
the scientific-name columns for species, genus, family, order, class and
phylum. The paper run used the custom V7 hierarchy recorded in the internal
provenance archive.

If the optional PathSeq workflow is used, the user must separately provide a
GATK 4.6.0.0 JAR, host k-mer/BWA images and the microbial BWA image, sequence
dictionary and taxonomy database. Public download/build instructions remain a
release blocker until database redistribution and provenance are approved.
