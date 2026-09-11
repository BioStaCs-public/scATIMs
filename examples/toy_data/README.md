# Toy data

The automated integration test creates a minimal synthetic PathSeq-tagged BAM,
barcode whitelist, CR-to-CB mapping and taxonomy table at runtime. No real
sample or patient data are included.

Run it from the repository root with:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

After the release interface is approved, a persistent command-line toy example
can be exported here from the same synthetic generator.
