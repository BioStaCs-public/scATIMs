# Reference resources

scATIMs does not redistribute human references, the GATK JAR, PathSeq indices
or the custom V7 microbial database in this repository.

The scMic step requires a taxonomy hierarchy table with at least `tax_id` and
the scientific-name columns for species, genus, family, order, class and
phylum. The paper run used the custom V7 hierarchy recorded in the internal
provenance archive.

The frozen V7 snapshot contains 130,359 unique genome accessions mapped to
4,458 tax_ids. Its hierarchy labels 3,196 tax_ids as Bacteria and 1,262 as
Viruses; it contains no Fungi-labelled tax_id. The runtime-defining SHA-256
values are:

| External artifact | SHA-256 |
|---|---|
| `pathseq_custom_V7_fixed.fa` | `33273dd36928df542bbea84486dd13715a0925eb858d9cccf87e1d0dbeb8890a` |
| `pathseq_custom_V7_fixed.fa.img` | `da200fddbe045b22e217d67c784cc89ea92bc53922c34ef0405879d28829da23` |
| `pathseq_custom_V7_fixed.dict` | `707b691a22853d47d6bbbee132d3091487323153153f92c6266890b5c00baa9a` |
| `pathseq_custom_V7_fixed.db` | `21cc709016f8d47a0465feb813e0662781bd0e10e20e3f35f4e23e155b3cadee` |
| `V7_taxid_hierarchy.txt` | `6d734afeabc0ee992fd25f9f44203cc3aad84a95741a68dc8ef0a2222db3e168` |

These checksums make the paper artifacts identifiable but do not by themselves
recover the historical taxon/accession selection rule or grant redistribution
permission.

For stage 05, the user must separately provide a
GATK 4.6.0.0 JAR, host k-mer/BWA images and the microbial BWA image, sequence
dictionary and taxonomy database. Public download/build instructions remain a
release blocker until database redistribution and provenance are approved.
