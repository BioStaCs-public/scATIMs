"""Unified scATIMs stage 00-07 workflow runner."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
import re
import shlex
import shutil

from . import atims, scmic
from .preprocess import HumanPreprocessor, quote, run_shell


STAGES = ("00", "01", "02", "03", "04", "04b", "05", "06", "07")


def read_samples(path):
    with Path(path).open(encoding="utf-8") as handle:
        samples = [line.strip() for line in handle if line.strip() and not line.lstrip().startswith("#")]
    for sample in samples:
        if sample in {".", ".."} or not re.fullmatch(r"[A-Za-z0-9._-]+", sample):
            raise ValueError(f"Unsafe sample identifier: {sample!r}")
    if len(samples) != len(set(samples)):
        raise ValueError("Sample list contains duplicate identifiers")
    return samples


def selected_stages(first="00", last="07"):
    if first not in STAGES or last not in STAGES:
        raise ValueError(f"Stages must be one of: {', '.join(STAGES)}")
    begin, end = STAGES.index(first), STAGES.index(last)
    if begin > end:
        raise ValueError("pipeline.start_step must precede pipeline.end_step")
    return STAGES[begin:end + 1]


def preflight(config, samples, stages, dry_run=False):
    """Validate the workflow contract before starting any non-dry-run stage."""
    for section in ("project", "library", "tools", "resources"):
        if section not in config:
            raise ValueError(f"Missing configuration section: {section}")
    chemistry = (str(config["library"].get("three_or_five_prime")),
                 str(config["library"].get("tenx_version")))
    from .preprocess import chemistry_parameters
    chemistry_parameters(*chemistry)
    if dry_run:
        return
    project, resources, tools = config["project"], config["resources"], config["tools"]
    Path(project["root_dir"]).mkdir(parents=True, exist_ok=True)
    required_tools = set()
    required_resources = []
    if set(stages) & {"00"}:
        required_tools.add("fastp")
    if set(stages) & {"01"}:
        required_tools.add("star")
        required_resources.append(("human_star_index", "dir"))
        _, _, whitelist_key = chemistry_parameters(*chemistry)
        if whitelist_key:
            required_resources.append((whitelist_key, "file"))
    if set(stages) & {"03", "04", "04b"}:
        required_tools.update({"bwa", "samtools"})
    if "03" in stages:
        required_resources.append(("host_fasta_grch38", "file"))
    if "04" in stages:
        required_resources.append(("host_fasta_t2t", "file"))
    if "04b" in stages and config.get("workflow", {}).get("vector_filter", True):
        required_resources.append(("vector_fasta", "file"))
    if "05" in stages:
        required_tools.add("java")
        required_resources += [("gatk_jar", "file"), ("pathseq_host_db", "dir"),
                               ("pathseq_microbe_db", "dir")]
    for key in sorted(required_tools):
        value = str(tools.get(key, key))
        resolved = Path(value) if "/" in value else shutil.which(value)
        if not resolved or not Path(resolved).is_file():
            raise FileNotFoundError(f"Executable not found for tools.{key}: {value}")
    for key, kind in required_resources:
        value = resources.get(key)
        path = Path(value) if value else None
        valid = path and (path.is_dir() if kind == "dir" else path.is_file())
        if not valid:
            raise FileNotFoundError(f"Missing {kind} for resources.{key}: {value}")
    if set(stages) & {"00", "01", "02"}:
        input_root = Path(project["input_dir"])
        patterns = config.get("inputs", {})
        for sample in samples:
            for key, default in (("r1", "{sample_id}.R1.fq.gz"),
                                 ("r2", "{sample_id}.R2.fq.gz")):
                path = input_root / sample / patterns.get(key, default).format(sample_id=sample)
                if not path.is_file():
                    raise FileNotFoundError(f"Missing {key.upper()} FASTQ for {sample}: {path}")


def _append(arguments, flag, value):
    if value is not None:
        arguments += [flag, str(value)]


def _pathseq_command(config, sample):
    project, resources = config["project"], config["resources"]
    tools, compute = config.get("tools", {}), config.get("compute", {})
    settings = config.get("pathseq", {})
    vector = config.get("workflow", {}).get("vector_filter", True)
    base = Path(project["root_dir"]) / sample
    input_bam = base / ((f"04b.vector/02.bam/{sample}.t2t_unmap.vector_unmapped.bam")
                        if vector else f"04.t2t/{sample}_rna_hg38_unmap_t2t_unmap.bam")
    # These are release contracts, not user-tunable output names: downstream
    # stages and archival tooling consume these exact directories.
    output = base / "05.pathseq_v10"
    host = Path(resources["pathseq_host_db"])
    microbe = Path(resources["pathseq_microbe_db"])
    prefix = output / f"{sample}_rna_output.pathseq"
    args = [tools.get("java", "java"), f"-Xmx{compute.get('pathseq_java_xmx', '64g')}",
            f"-Xms{compute.get('pathseq_java_xms', '32g')}", "-jar", resources["gatk_jar"],
            "PathSeqPipelineSpark", "--spark-master", compute.get("pathseq_spark_master", "local[10]"),
            "--conf", f"spark.driver.memory={compute.get('pathseq_driver_memory', '80g')}",
            "--conf", f"spark.executor.memory={compute.get('pathseq_executor_memory', '20g')}",
            "--conf", f"spark.executor.cores={compute.get('pathseq_executor_cores', 5)}",
            "--input", input_bam, "--min-clipped-read-length", settings.get("min_clipped_read_length", 60),
            "--kmer-file", host / settings.get("host_kmer", "GCA_000001405.15_GRCh38_no_alt_analysis_set.hss"),
            "--filter-bwa-image", host / settings.get("host_image", "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.img"),
            "--microbe-bwa-image", microbe / settings.get("microbe_image", "pathseq_custom_V7_fixed.fa.img"),
            "--microbe-dict", microbe / settings.get("microbe_dict", "pathseq_custom_V7_fixed.dict"),
            "--taxonomy-file", microbe / settings.get("taxonomy_file", "pathseq_custom_V7_fixed.db"),
            "--output", f"{prefix}.bam", "--scores-output", f"{prefix}.txt",
            "--score-metrics", f"{prefix}.score_metrics.txt",
            "--filter-metrics", f"{prefix}.filter_metrics.txt",
            "--is-host-aligned", "false", "--filter-duplicates", "false",
            "--min-score-identity", settings.get("min_score_identity", 0.7)]
    return " ".join(quote(item) for item in args), output, Path(f"{prefix}.bam"), Path(f"{prefix}.txt")


def _run_pathseq_sample(config, sample, log_dir, dry_run, force):
    command, output, bam, scores = _pathseq_command(config, sample)
    if bam.is_file() and bam.stat().st_size > 0 and scores.is_file() and scores.stat().st_size > 0 and not force:
        return
    if not dry_run:
        output.mkdir(parents=True, exist_ok=True)
    run_shell(command, log_dir / f"{sample}.05.log", dry_run)
    if not dry_run and (not bam.is_file() or bam.stat().st_size == 0 or
                        not scores.is_file() or scores.stat().st_size == 0):
        raise RuntimeError(f"PathSeq completed without required outputs for {sample}")


def run(config_path, dry_run=False, force=False):
    import yaml
    config_path = Path(config_path).resolve()
    with config_path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    if config.get("library", {}).get("species", "human") != "human":
        raise ValueError("This release candidate validates only the human workflow")
    samples = read_samples(config["project"]["sample_list"])
    if not samples:
        raise ValueError("Sample list contains no samples")
    pipeline = config.get("pipeline", {})
    stages = selected_stages(str(pipeline.get("start_step", "00")),
                             str(pipeline.get("end_step", "07")))
    preflight(config, samples, stages, dry_run)
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    log_dir = Path(config["project"]["root_dir"]) / ".scatims_runs" / run_id
    if not dry_run:
        log_dir.mkdir(parents=True, exist_ok=True)
    print(f"scATIMs stages: {','.join(stages)}; samples: {len(samples)}; dry_run={dry_run}")

    for stage in stages:
        print(f"[scATIMs] stage {stage}")
        if stage in {"00", "01", "02", "03", "04", "04b"}:
            for sample in samples:
                processor = HumanPreprocessor(config, sample, dry_run, force, log_dir)
                getattr(processor, f"stage{stage}")()
        elif stage == "05":
            workers = max(1, int(config.get("compute", {}).get("pathseq_sample_workers", 4)))
            with ThreadPoolExecutor(max_workers=workers) as executor:
                futures = [executor.submit(_run_pathseq_sample, config, sample, log_dir,
                                           dry_run, force) for sample in samples]
                for future in futures:
                    future.result()
        elif stage == "06":
            section, project = config["scmic"], config["project"]
            args = ["--root-dir", str(project["root_dir"]), "--sample-list", str(project["sample_list"]),
                    "--taxonomy-hierarchy", str(section["taxonomy_hierarchy"])]
            args += ["--output-subdir", "06.scMic_6_V10_EM_Dis",
                     "--pathseq-bam-relative", "05.pathseq_v10/{sample_id}_rna_output.pathseq.bam"]
            for key, flag in (("workers", "--workers"),
                              ("barcode_whitelist_relative", "--barcode-whitelist-relative"),
                              ("barcode_mapping_relative", "--barcode-mapping-relative"),
                              ("error_log_name", "--error-log-name")):
                _append(args, flag, section.get(key))
            print("scatims scmic", shlex.join(args)) if dry_run else scmic.main(args)
        elif stage == "07":
            section, project = config["atims"], config["project"]
            args = ["--root-dir", str(project["root_dir"]), "--sample-list", str(project["sample_list"])]
            args += ["--scmic-dir-name", "06.scMic_6_V10_EM_Dis",
                     "--output-dir-name", "07.ATIMs_9_V10_EM_Dis"]
            for key, flag in (("tax_levels", "--tax-levels"), ("workers", "--workers"),
                              ("barcode_whitelist_relative", "--barcode-whitelist-relative"),
                              ("barcode_mapping_relative", "--barcode-mapping-relative"),
                              ("error_log_name", "--error-log-name"), ("fisher_pval", "--fisher_pval"),
                              ("mw_pval", "--mw_pval"), ("dispersion_pval", "--dispersion_pval"),
                              ("enrichment_ratio", "--enrichment_ratio"), ("log2fc", "--log2fc")):
                _append(args, flag, section.get(key))
            print("scatims atims", shlex.join(args)) if dry_run else atims.main(args)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    run(args.config, args.dry_run, args.force)
    return 0


if __name__ == "__main__":
    main()
