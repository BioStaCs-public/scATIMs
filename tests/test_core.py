import gzip
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd
import pysam

from scatims import atims, cli, cr_cb, scmic


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


if __name__ == "__main__":
    unittest.main()
