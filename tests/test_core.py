import gzip
import io
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd
import pysam

from scatims import atims, cli, cr_cb, preprocess, scmic, workflow


class ScMicCoreTests(unittest.TestCase):
    def test_directional_umi_collapse(self):
        observed = ["AAAA", "AAAA", "AAAA", "AAAT"]
        mapping = scmic.umi_tools_directional_collapse(observed)
        self.assertEqual(mapping["AAAA"], "AAAA")
        self.assertEqual(mapping["AAAT"], "AAAA")

    def test_configurable_sample_path(self):
        path = scmic._resolve_sample_path(Path("/project/S1"), "inputs/{sample_id}.bam", "S1")
        self.assertEqual(path, Path("/project/S1/inputs/S1.bam"))


class AtimsCoreTests(unittest.TestCase):
    def test_paper_threshold_filters_cell_enrichment(self):
        cell_rows = []
        for index in range(20):
            count = 20 if index == 0 else 1
            for umi_index in range(count):
                cell_rows.append(
                    {"cell_barcode": f"C{index:02d}", "UMI": f"U{umi_index}", "feature_name": "Taxon A", "score": 1.0}
                )
        empty_rows = [{"cell_barcode": "E00", "UMI": "U0", "feature_name": "Taxon A", "score": 1.0}]
        result = atims.analyze_features(
            pd.DataFrame(cell_rows), pd.DataFrame(empty_rows), 20, 20, "species", 0.0
        )
        thresholds = {
            "fisher_pval": 0.05,
            "mw_pval": 0.05,
            "dispersion_pval": 0.05,
            "enrichment_ratio": 1.0,
            "log2fc": 1.0,
        }
        cell_enriched, empty_enriched, _ = atims.filter_microbes(result, thresholds)
        self.assertEqual(cell_enriched["feature_name"].tolist(), ["Taxon A"])
        self.assertTrue(empty_enriched.empty)


class EndToEndToyTest(unittest.TestCase):
    @staticmethod
    def _write_read(bam, sample_id, raw_barcode, umi, position):
        read = pysam.AlignedSegment()
        read.query_name = f"{sample_id}_read{position}_{raw_barcode}_{umi}"
        read.query_sequence = "A" * 50
        read.flag = 0
        read.reference_id = 0
        read.reference_start = position
        read.mapping_quality = 60
        read.cigar = ((0, 50),)
        read.query_qualities = pysam.qualitystring_to_array("I" * 50)
        read.set_tag("YP", "1", value_type="Z")
        bam.write(read)

    def test_scmic_to_atims(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sample_id = "S1"
            sample_dir = root / sample_id
            (sample_dir / "05.pathseq_v10").mkdir(parents=True)
            (sample_dir / "01.starSolo/Solo.out/Gene/filtered_gz").mkdir(parents=True)
            (sample_dir / "01.starSolo_bam").mkdir(parents=True)

            taxonomy = root / "taxonomy.tsv"
            taxonomy.write_text(
                "tax_id\tspecies_scientific_name\tgenus_scientific_name\t"
                "family_scientific_name\torder_scientific_name\t"
                "class_scientific_name\tphylum_scientific_name\n"
                "1\tTaxon alpha\tTaxon\tTaxonaceae\tTaxonales\tTaxonia\tTaxonota\n",
                encoding="utf-8",
            )

            cell_barcodes = [f"C{index:02d}" for index in range(20)]
            empty_barcodes = [f"E{index:02d}" for index in range(20)]
            mapping_lines = ["CR_Raw\tCB_Corrected"]
            mapping_lines.extend(f"R{barcode}\t{barcode}" for barcode in cell_barcodes + empty_barcodes)
            (sample_dir / "01.starSolo_bam/CR_CB_mapping.tsv").write_text(
                "\n".join(mapping_lines) + "\n", encoding="utf-8"
            )
            with gzip.open(
                sample_dir / "01.starSolo/Solo.out/Gene/filtered_gz/barcodes.tsv.gz", "wt"
            ) as handle:
                handle.write("\n".join(cell_barcodes) + "\n")

            bam_path = sample_dir / "05.pathseq_v10/S1_rna_output.pathseq.bam"
            header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "microbe", "LN": 10000}]}
            with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
                position = 100
                for index, barcode in enumerate(cell_barcodes):
                    copies = 20 if index == 0 else 1
                    for copy_index in range(copies):
                        self._write_read(
                            bam,
                            sample_id,
                            f"R{barcode}",
                            f"U{copy_index:03d}",
                            position,
                        )
                        position += 1
                self._write_read(bam, sample_id, "RE00", "U000", position)

            ok, message = scmic.process_sample(root, sample_id, taxonomy, "06.scMic_test")
            self.assertTrue(ok, message)

            thresholds = {
                "fisher_pval": 0.05,
                "mw_pval": 0.05,
                "dispersion_pval": 0.05,
                "enrichment_ratio": 1.0,
                "log2fc": 1.0,
            }
            ok, message = atims.run_one_sample(
                root,
                sample_id,
                ["species"],
                thresholds,
                0.0,
                "06.scMic_test",
                "07.ATIMs_test",
            )
            self.assertTrue(ok, message)

            output_dir = sample_dir / "07.ATIMs_test"
            expected = [
                output_dir / "S1_ATIMs_species_Statistics.tsv",
                output_dir / "S1_ATIMs_species_CellEnriched.tsv",
                output_dir / "S1_ATIMs_species_EmptyEnriched.tsv",
                output_dir / "S1_ATIMs_species_Parameters.json",
            ]
            self.assertTrue(all(path.exists() for path in expected))
            cell_enriched = pd.read_csv(expected[1], sep="\t")
            self.assertEqual(cell_enriched["feature_name"].tolist(), ["Taxon alpha"])
            parameters = json.loads(expected[3].read_text(encoding="utf-8"))
            self.assertEqual(parameters["min_evidence_score"], 0.0)


class UnifiedCliTests(unittest.TestCase):
    def test_paper_config_forwards_locked_parameter(self):
        config = {
            "project": {"root_dir": "/project", "sample_list": "/samples.txt"},
            "scmic": {"taxonomy_hierarchy": "/taxonomy.tsv", "workers": 1},
            "atims": {"tax_levels": "all", "min_evidence_score": 0.0, "workers": 1},
        }
        fake_yaml = types.SimpleNamespace(safe_load=lambda handle: config)
        with tempfile.TemporaryDirectory() as tmp:
            config_path = Path(tmp) / "config.yaml"
            config_path.write_text("test fixture", encoding="utf-8")
            with mock.patch.dict(sys.modules, {"yaml": fake_yaml}), mock.patch.object(
                cli.scmic, "main"
            ) as scmic_main, mock.patch.object(cli.atims, "main") as atims_main:
                cli.run_from_config(config_path)
        scmic_main.assert_called_once()
        atims_args = atims_main.call_args.args[0]
        index = atims_args.index("--min_evidence_score")
        self.assertEqual(atims_args[index + 1], "0.0")


class CrCbMappingTests(unittest.TestCase):
    def test_deterministic_mapping(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bam_path = root / "tagged.bam"
            header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 100}]}
            with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
                for index, (cr, cb) in enumerate((("RAW2", "CB2"), ("RAW1", "CB1"))):
                    read = pysam.AlignedSegment()
                    read.query_name = f"r{index}"
                    read.query_sequence = "A" * 20
                    read.flag = 4
                    read.query_qualities = pysam.qualitystring_to_array("I" * 20)
                    read.set_tag("CR", cr)
                    read.set_tag("CB", cb)
                    bam.write(read)
            output = root / "CR_CB_mapping.tsv"
            metrics = cr_cb.extract_mapping(bam_path, output)
            self.assertEqual(metrics["mapping_rows"], 2)
            self.assertEqual(output.read_text().splitlines(), [
                "CR_Raw\tCB_Corrected", "RAW1\tCB1", "RAW2\tCB2"
            ])
            self.assertIn("status\tsuccess", (root / "CR_CB_mapping_audit.tsv").read_text())

    def test_conflicting_raw_barcode_fails_with_audit(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bam_path = root / "tagged.bam"
            header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 100}]}
            with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
                for index, cb in enumerate(("CB1", "CB2")):
                    read = pysam.AlignedSegment()
                    read.query_name = f"r{index}"
                    read.query_sequence = "A" * 20
                    read.flag = 4
                    read.query_qualities = pysam.qualitystring_to_array("I" * 20)
                    read.set_tag("CR", "RAW1")
                    read.set_tag("CB", cb)
                    bam.write(read)
            output = root / "CR_CB_mapping.tsv"
            with self.assertRaisesRegex(RuntimeError, "conflicting raw barcodes"):
                cr_cb.extract_mapping(bam_path, output)
            self.assertFalse(output.exists())
            self.assertIn(
                "status\tfailed_conflicting_cr",
                (root / "CR_CB_mapping_audit.tsv").read_text(encoding="utf-8"),
            )
            conflicts = (root / "CR_CB_mapping_conflicts.tsv").read_text(encoding="utf-8")
            self.assertIn("RAW1\tCB1\t1", conflicts)
            self.assertIn("RAW1\tCB2\t1", conflicts)


class PreprocessingTests(unittest.TestCase):
    def test_paper_chemistry_coordinates(self):
        self.assertEqual(preprocess.chemistry_parameters("3", "v2")[:2], (16, 10))
        self.assertEqual(preprocess.chemistry_parameters("3", "v3")[:2], (16, 12))
        self.assertEqual(preprocess.chemistry_parameters("5", "v3")[:2], (16, 12))

    def test_star_identity_command_contains_cr_cb(self):
        command = preprocess.star_solo_command(
            "STAR", "/ref", "R1.gz", "R2.gz", "/out", "/wl", 8, "3", "v3", True
        )
        self.assertIn("--outSAMattributes CR CB", command)
        self.assertIn("--soloCBlen 16", command)
        self.assertIn("--soloUMIlen 12", command)
        five_prime = preprocess.star_solo_command(
            "STAR", "/ref", "R1.gz", "R2.gz", "/out", "/wl", 8, "5", "v3"
        )
        self.assertIn("--soloBarcodeMate 1", five_prime)
        self.assertNotIn("--soloBarcodeReadLength", five_prime)
        c4 = preprocess.star_solo_command(
            "STAR", "/ref", "R1.gz", "R2.gz", "/out", None, 8, "3", "C4"
        )
        self.assertIn("--soloCBwhitelist None", c4)
        self.assertIn("--clipAdapterType CellRanger4", c4)

    def test_streaming_identity_embedding(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            r1, r2, output = root / "R1.gz", root / "R2.gz", root / "out.gz"
            with gzip.open(r1, "wt") as handle:
                handle.write("@read1 1:N:0:1\nAAAACCCC\n+\nIIIIIIII\n")
            with gzip.open(r2, "wt") as handle:
                handle.write("@read1 2:N:0:1\nTTTT\n+\nIIII\n")
            barcodes, stats = root / "barcodes.tsv", root / "stats.tsv"
            count = preprocess.embed_barcode_umi(r1, r2, output, 4, 4, barcodes, stats)
            self.assertEqual(count, 1)
            with gzip.open(output, "rt") as handle:
                self.assertEqual(handle.readline().strip(), "@read1_AAAA_CCCC 2:N:0:1")
            self.assertEqual(barcodes.read_text(), "AAAA\n")
            self.assertEqual(stats.read_text(), "barcode\treads\nAAAA\t1\n")

    def test_stage_range(self):
        self.assertEqual(workflow.selected_stages("03", "05"), ("03", "04", "04b", "05"))
        with self.assertRaises(ValueError):
            workflow.selected_stages("07", "00")

    def test_sample_identifiers_are_safe_and_unique(self):
        with tempfile.TemporaryDirectory() as temporary:
            sample_list = Path(temporary) / "samples.txt"
            sample_list.write_text("S1\nS1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate"):
                workflow.read_samples(sample_list)
            sample_list.write_text("../outside\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Unsafe"):
                workflow.read_samples(sample_list)

    def test_resume_requires_all_nonempty_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            nonempty, empty = root / "nonempty", root / "empty"
            nonempty.write_text("ok\n", encoding="utf-8")
            empty.touch()
            self.assertTrue(preprocess._complete(nonempty))
            self.assertFalse(preprocess._complete(nonempty, empty))
            with self.assertRaisesRegex(RuntimeError, "non-empty"):
                preprocess._require(empty)

    def test_full_preprocessing_dry_run_does_not_create_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample_list = root / "samples.txt"
            sample_list.write_text("S1\n")
            config_path = root / "config.yaml"
            config_path.write_text("test fixture\n")
            output = root / "output"
            config = {
                "project": {"input_dir": str(root / "input"), "root_dir": str(output),
                            "sample_list": str(sample_list)},
                "library": {"species": "human", "three_or_five_prime": "3", "tenx_version": "v3"},
                "tools": {"fastp": "fastp", "star": "STAR", "bwa": "bwa",
                          "samtools": "samtools", "java": "java"},
                "resources": {"human_star_index": "/ref/star", "whitelist_2018": "/ref/wl",
                              "host_fasta_grch38": "/ref/hg38", "host_fasta_t2t": "/ref/t2t",
                              "vector_fasta": "/ref/vector", "gatk_jar": "/ref/gatk.jar",
                              "pathseq_host_db": "/ref/host", "pathseq_microbe_db": "/ref/microbe"},
                "compute": {"threads": 2, "pathseq_sample_workers": 1},
                "workflow": {"vector_filter": True},
            }
            fake_yaml = types.SimpleNamespace(safe_load=lambda handle: config)
            stream = io.StringIO()
            with mock.patch.dict(sys.modules, {"yaml": fake_yaml}), mock.patch("sys.stdout", stream):
                workflow.run(config_path, "00", "05", dry_run=True)
            rendered = stream.getvalue()
            self.assertIn("fastp", rendered)
            self.assertIn("--outSAMattributes CR CB", rendered)
            self.assertIn("scatims embed", rendered)
            self.assertIn("/04b.vector/02.bam/S1.t2t_unmap.vector_unmapped.bam", rendered)
            self.assertIn("PathSeqPipelineSpark", rendered)
            self.assertIn("--min-clipped-read-length 60", rendered)
            self.assertIn("--is-host-aligned false", rendered)
            self.assertIn("--filter-duplicates false", rendered)
            self.assertIn("--min-score-identity 0.7", rendered)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
