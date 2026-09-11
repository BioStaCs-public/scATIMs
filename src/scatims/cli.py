"""Unified command-line interface for the two scATIMs steps."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__, atims, scmic


def _value(mapping: dict, key: str, section: str):
    value = mapping.get(key)
    if value is None or value == "":
        raise ValueError(f"Missing required configuration value: {section}.{key}")
    return value


def _append(args: list[str], flag: str, value) -> None:
    if value is not None:
        args.extend([flag, str(value)])


def run_from_config(config_path: str | Path) -> None:
    import yaml

    with Path(config_path).open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}

    project = config.get("project", {})
    scmic_config = config.get("scmic", {})
    atims_config = config.get("atims", {})

    root_dir = _value(project, "root_dir", "project")
    sample_list = _value(project, "sample_list", "project")
    taxonomy = _value(scmic_config, "taxonomy_hierarchy", "scmic")

    scmic_args: list[str] = [
        "--root-dir", str(root_dir),
        "--sample-list", str(sample_list),
        "--taxonomy-hierarchy", str(taxonomy),
    ]
    for key, flag in (
        ("output_subdir", "--output-subdir"),
        ("workers", "--workers"),
        ("pathseq_bam_relative", "--pathseq-bam-relative"),
        ("barcode_whitelist_relative", "--barcode-whitelist-relative"),
        ("barcode_mapping_relative", "--barcode-mapping-relative"),
        ("error_log_name", "--error-log-name"),
    ):
        _append(scmic_args, flag, scmic_config.get(key))

    atims_args: list[str] = [
        "--root-dir", str(root_dir),
        "--sample-list", str(sample_list),
    ]
    for key, flag in (
        ("tax_levels", "--tax-levels"),
        ("scmic_dir_name", "--scmic-dir-name"),
        ("output_dir_name", "--output-dir-name"),
        ("workers", "--workers"),
        ("barcode_whitelist_relative", "--barcode-whitelist-relative"),
        ("barcode_mapping_relative", "--barcode-mapping-relative"),
        ("error_log_name", "--error-log-name"),
        ("fisher_pval", "--fisher_pval"),
        ("mw_pval", "--mw_pval"),
        ("dispersion_pval", "--dispersion_pval"),
        ("enrichment_ratio", "--enrichment_ratio"),
        ("log2fc", "--log2fc"),
        ("min_evidence_score", "--min_evidence_score"),
    ):
        _append(atims_args, flag, atims_config.get(key))

    print("[scATIMs] Step 1/2: scMic")
    scmic.main(scmic_args)
    print("[scATIMs] Step 2/2: ATIMs")
    atims.main(atims_args)


def _print_help() -> None:
    print(
        "scATIMs: scMic microbial-signal assignment followed by ATIMs decontamination\n\n"
        "Usage:\n"
        "  scatims run --config CONFIG.yaml\n"
        "  scatims scmic [scMic options]\n"
        "  scatims atims [ATIMs options]\n"
        "  scatims --version"
    )


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in {"-h", "--help"}:
        _print_help()
        return 0
    if args[0] in {"-V", "--version"}:
        print(f"scATIMs {__version__}")
        return 0
    if args[0] == "scmic":
        return scmic.main(args[1:])
    if args[0] == "atims":
        return atims.main(args[1:])
    if args[0] == "run":
        parser = argparse.ArgumentParser(prog="scatims run")
        parser.add_argument("--config", required=True)
        parsed = parser.parse_args(args[1:])
        run_from_config(parsed.config)
        return 0
    raise SystemExit(f"Unknown scATIMs command: {args[0]!r}. Use --help.")


if __name__ == "__main__":
    main()
