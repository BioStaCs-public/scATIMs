"""Generate a deterministic raw-to-corrected barcode map from STARsolo tags."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import os
from pathlib import Path

import pysam


def _atomic_write(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        handle.writelines(lines)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def extract_mapping(bam_path, mapping_path, audit_path=None, conflict_path=None,
                    decompress_threads=2, fail_on_conflict=False):
    """Extract the legacy CR/CB mapping table with deterministic audit outputs.

    The original ``run_CR2CB_multi_samples_del.sh`` implementation emitted every
    unique valid ``(CR_Raw, CB_Corrected)`` pair via ``sort -u``.  It did not
    reject a raw barcode that had more than one corrected-barcode value.  Keep
    that table-level behavior here: downstream legacy scMic code reads the
    rows sequentially into a dictionary, so for a conflicting CR the final
    lexicographically sorted row is the effective value.  Conflicts remain in
    a separate audit file rather than being silently hidden.

    ``fail_on_conflict`` is retained only for an explicit validation mode; the
    workflow itself uses the historical non-failing behavior.
    """
    bam_path, mapping_path = Path(bam_path), Path(mapping_path)
    audit_path = Path(audit_path or mapping_path.with_name("CR_CB_mapping_audit.tsv"))
    conflict_path = Path(conflict_path or mapping_path.with_name("CR_CB_mapping_conflicts.tsv"))
    if not bam_path.is_file() or bam_path.stat().st_size == 0:
        raise FileNotFoundError(f"STARsolo BAM is missing or empty: {bam_path}")
    counts = defaultdict(Counter)
    metrics = {key: 0 for key in ("records_total", "records_missing_cr", "records_missing_cb",
                                  "records_invalid_cb", "records_valid_pair")}
    metrics.update({"bam_path": str(bam_path), "bam_size_bytes": bam_path.stat().st_size})
    with pysam.AlignmentFile(str(bam_path), "rb", threads=max(1, int(decompress_threads))) as bam:
        for read in bam.fetch(until_eof=True):
            metrics["records_total"] += 1
            if not read.has_tag("CR"):
                metrics["records_missing_cr"] += 1
                continue
            if not read.has_tag("CB"):
                metrics["records_missing_cb"] += 1
                continue
            cr, cb = str(read.get_tag("CR")).strip(), str(read.get_tag("CB")).strip()
            if not cr:
                metrics["records_missing_cr"] += 1
                continue
            if not cb or cb == "-":
                metrics["records_invalid_cb"] += 1
                continue
            counts[cr][cb] += 1
            metrics["records_valid_pair"] += 1
    conflicts = {cr: value for cr, value in counts.items() if len(value) > 1}
    metrics.update({"unique_cr": len(counts),
                    "unique_cr_cb_pairs": sum(len(value) for value in counts.values()),
                    "conflicting_cr": len(conflicts),
                    "fail_on_conflict": str(bool(fail_on_conflict)).lower()})
    conflict_lines = ["CR_Raw\tCB_Corrected\tn_records\n"]
    for cr in sorted(conflicts):
        for cb, count in sorted(conflicts[cr].items()):
            conflict_lines.append(f"{cr}\t{cb}\t{count}\n")
    _atomic_write(conflict_path, conflict_lines)

    def write_audit(status):
        metrics["status"] = status
        _atomic_write(audit_path, ["metric\tvalue\n"] +
                      [f"{key}\t{value}\n" for key, value in metrics.items()])

    if not counts:
        write_audit("failed_no_valid_mapping")
        raise RuntimeError(f"No valid CR/CB pairs found in {bam_path}")
    if conflicts and fail_on_conflict:
        write_audit("failed_conflicting_cr")
        raise RuntimeError(f"Found {len(conflicts)} conflicting raw barcodes; see {conflict_path}")
    pairs = sorted((cr, cb) for cr, values in counts.items() for cb in values)
    _atomic_write(mapping_path, ["CR_Raw\tCB_Corrected\n"] +
                  [f"{cr}\t{cb}\n" for cr, cb in pairs])
    metrics.update({"mapping_rows": len(pairs),
                    "mapping_policy": "legacy_all_unique_cr_cb_pairs",
                    "mapping_path": str(mapping_path),
                    "conflict_path": str(conflict_path)})
    write_audit("success")
    return metrics


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--audit")
    parser.add_argument("--conflicts")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument(
        "--fail-on-conflict", action="store_true",
        help="Validation only: stop instead of reproducing the legacy all-pairs mapping behavior.",
    )
    args = parser.parse_args(argv)
    extract_mapping(args.bam, args.output, args.audit, args.conflicts,
                    args.threads, args.fail_on_conflict)
    return 0


if __name__ == "__main__":
    main()
