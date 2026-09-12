"""Human scRNA-seq preprocessing stages 00-04b for scATIMs."""

from __future__ import annotations

import gzip
from collections import Counter
from pathlib import Path
import shlex
import shutil
import subprocess

from .cr_cb import extract_mapping


CHEMISTRIES = {
    ("3", "v1"): (14, 10, "whitelist_2014"),
    ("3", "v2"): (16, 10, "whitelist_2016"),
    ("3", "v3"): (16, 12, "whitelist_2018"),
    ("3", "v3.1"): (16, 12, "whitelist_2018"),
    ("3", "C4"): (25, 10, None),
    ("5", "v1"): (16, 10, "whitelist_2016"),
    ("5", "v1.1"): (16, 10, "whitelist_2016"),
    ("5", "v2"): (16, 10, "whitelist_2016"),
    ("5", "v3"): (16, 12, "whitelist_5p_2023"),
}


def chemistry_parameters(prime: str, version: str):
    try:
        return CHEMISTRIES[(str(prime), str(version))]
    except KeyError as error:
        raise ValueError(f"Unsupported 10x chemistry: {prime!r} prime, {version!r}") from error


def quote(value) -> str:
    return shlex.quote(str(value))


def run_shell(command: str, log_path: Path, dry_run: bool = False) -> None:
    if dry_run:
        print(command)
        return
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as log:
        log.write(f"$ {command}\n")
        completed = subprocess.run(
            ["bash", "-o", "pipefail", "-c", command],
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )
    if completed.returncode:
        raise RuntimeError(f"Command failed ({completed.returncode}); see {log_path}")


def star_solo_command(star, index, r1, r2, output, whitelist, threads,
                      prime, version, tagged_bam=False, bam_sort_ram=60000000000):
    cb_len, umi_len, _ = chemistry_parameters(prime, version)
    reads = [r2, r1] if str(prime) == "3" else [r1, r2]
    command = [
        star, "--genomeDir", index, "--runThreadN", threads,
        "--outFileNamePrefix", f"{output}/", "--readFilesIn", *reads,
        "--readFilesCommand", "zcat", "--soloType", "CB_UMI_Simple",
        "--soloCBstart", 1, "--soloCBlen", cb_len,
        "--soloUMIstart", cb_len + 1, "--soloUMIlen", umi_len,
        "--soloCBwhitelist", whitelist or "None",
        "--soloUMIdedup", "1MM_CR",
        "--soloCBmatchWLtype", "1MM_multi_Nbase_pseudocounts",
        "--soloUMIfiltering", "MultiGeneUMI_CR",
        "--outFilterScoreMin", 30,
    ]
    if str(prime) == "5":
        command += ["--soloBarcodeMate", 1, "--clip5pNbases", 39, 0]
    else:
        command += ["--soloBarcodeReadLength", 0, "--clipAdapterType", "CellRanger4"]
    if tagged_bam:
        command += [
            "--genomeLoad", "NoSharedMemory", "--outSAMtype", "BAM", "SortedByCoordinate",
            "--outSAMunmapped", "Within", "--outSAMattributes", "CR", "CB",
            "--limitBAMsortRAM", bam_sort_ram, "--soloCellFilter", "None",
        ]
    else:
        command += [
            "--genomeLoad", "NoSharedMemory", "--outSAMtype", "None",
            "--soloCellFilter", "EmptyDrops_CR", "--soloFeatures", "Gene",
            "--soloMultiMappers", "EM",
        ]
    return " ".join(quote(item) for item in command)


def gzip_solo_tables(solo_gene_dir: Path) -> None:
    for source_name in ("filtered", "raw"):
        source = solo_gene_dir / source_name
        destination = solo_gene_dir / f"{source_name}_gz"
        destination.mkdir(parents=True, exist_ok=True)
        for path in source.iterdir():
            if not path.is_file():
                continue
            target = destination / f"{path.name}.gz"
            with path.open("rb") as reader, gzip.open(target, "wb") as writer:
                shutil.copyfileobj(reader, writer)


def _fastq_record(handle):
    record = [handle.readline() for _ in range(4)]
    if not record[0]:
        return None
    if any(not line for line in record[1:]):
        raise ValueError("Truncated FASTQ record")
    return record


def embed_barcode_umi(r1_path, r2_path, output_path, barcode_length, umi_length,
                      barcode_table=None, barcode_stats=None):
    """Stream paired FASTQ and append R1 raw barcode/UMI to each R2 read name."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(f".{output_path.name}.tmp")
    records, barcode_counts = 0, Counter()
    with gzip.open(r1_path, "rt") as r1, gzip.open(r2_path, "rt") as r2, gzip.open(temporary, "wt") as output:
        while True:
            left, right = _fastq_record(r1), _fastq_record(r2)
            if left is None and right is None:
                break
            if left is None or right is None:
                raise ValueError("R1 and R2 contain different numbers of FASTQ records")
            left_name = left[0].split()[0].removeprefix("@")
            right_name = right[0].split()[0].removeprefix("@")
            if left_name.split("/")[0] != right_name.split("/")[0]:
                raise ValueError(f"FASTQ read-name mismatch: {left_name} != {right_name}")
            sequence = left[1].strip()
            if len(sequence) < barcode_length + umi_length:
                raise ValueError(f"R1 sequence shorter than barcode+UMI for {left_name}")
            barcode = sequence[:barcode_length]
            umi = sequence[barcode_length:barcode_length + umi_length]
            barcode_counts[barcode] += 1
            header_tail = right[0].strip().split(maxsplit=1)
            header = f"@{right_name}_{barcode}_{umi}"
            if len(header_tail) == 2:
                header += f" {header_tail[1]}"
            output.write(header + "\n" + "".join(right[1:]))
            records += 1
    temporary.replace(output_path)
    if barcode_table:
        Path(barcode_table).write_text("".join(f"{barcode}\n" for barcode in sorted(barcode_counts)), encoding="utf-8")
    if barcode_stats:
        Path(barcode_stats).write_text(
            "barcode\treads\n" + "".join(f"{barcode}\t{barcode_counts[barcode]}\n" for barcode in sorted(barcode_counts)),
            encoding="utf-8",
        )
    return records


def _complete(*paths):
    return all(Path(path).is_file() and Path(path).stat().st_size > 0 for path in paths)


def _require(*paths):
    missing = [str(path) for path in paths if not _complete(path)]
    if missing:
        raise RuntimeError("Stage completed without required non-empty output(s): " + ", ".join(missing))


class HumanPreprocessor:
    """Execute or render stages 00-04b for one human sample."""

    def __init__(self, config, sample, dry_run=False, force=False, log_dir=None):
        self.cfg, self.sample = config, sample
        self.dry_run, self.force = dry_run, force
        self.input = Path(config["project"]["input_dir"]) / sample
        self.base = Path(config["project"]["root_dir"]) / sample
        self.tools = config["tools"]
        self.resources = config["resources"]
        self.compute = config.get("compute", {})
        self.library = config["library"]
        self.log_dir = Path(log_dir or Path(config["project"]["root_dir"]) / ".scatims_logs")

    def _paths(self):
        r1 = self.input / self.cfg.get("inputs", {}).get("r1", "{sample_id}.R1.fq.gz").format(sample_id=self.sample)
        r2 = self.input / self.cfg.get("inputs", {}).get("r2", "{sample_id}.R2.fq.gz").format(sample_id=self.sample)
        qc = self.base / "00.QC" / self.sample
        return r1, r2, qc / f"{self.sample}.QC.R1.fq.gz", qc / f"{self.sample}.QC.R2.fq.gz"

    def _run(self, stage, command):
        run_shell(command, self.log_dir / f"{self.sample}.{stage}.log", self.dry_run)

    def stage00(self):
        r1, r2, q1, q2 = self._paths()
        report_json, report_html = q1.parent / f"{self.sample}.fastp.json", q1.parent / f"{self.sample}.fastp.html"
        if _complete(q1, q2, report_json, report_html) and not self.force:
            return
        if not self.dry_run:
            q1.parent.mkdir(parents=True, exist_ok=True)
        command = [self.tools.get("fastp", "fastp"), "-w", self.compute.get("threads", 32),
                   "-i", r1, "-I", r2, "-o", q1, "-O", q2,
                   "-j", report_json, "-h", report_html]
        self._run("00", " ".join(quote(item) for item in command))
        if not self.dry_run:
            _require(q1, q2, report_json, report_html)

    def stage01(self):
        _, _, r1, r2 = self._paths()
        prime, version = str(self.library["three_or_five_prime"]), str(self.library["tenx_version"])
        _, _, whitelist_key = chemistry_parameters(prime, version)
        whitelist = self.resources.get(whitelist_key) if whitelist_key else None
        expression = self.base / "01.starSolo"
        filtered = expression / "Solo.out/Gene/filtered_gz/barcodes.tsv.gz"
        if not _complete(filtered) or self.force:
            if not self.dry_run:
                expression.mkdir(parents=True, exist_ok=True)
            self._run("01.expression", star_solo_command(
                self.tools.get("star", "STAR"), self.resources["human_star_index"], r1, r2,
                expression, whitelist, self.compute.get("threads", 32), prime, version))
            if not self.dry_run:
                gzip_solo_tables(expression / "Solo.out/Gene")
                _require(filtered)
        identity = self.base / "01.starSolo_bam"
        mapping = identity / "CR_CB_mapping.tsv"
        if not _complete(mapping) or self.force:
            if not self.dry_run:
                identity.mkdir(parents=True, exist_ok=True)
            bam = identity / "Aligned.sortedByCoord.out.bam"
            self._run("01.identity", star_solo_command(
                self.tools.get("star", "STAR"), self.resources["human_star_index"], r1, r2,
                identity, whitelist, self.compute.get("threads", 32), prime, version, True,
                self.compute.get("star_bam_sort_ram", 60000000000)))
            if not self.dry_run:
                extract_mapping(bam, mapping,
                                decompress_threads=self.compute.get("bam_threads", 2),
                                fail_on_conflict=True)
                _require(mapping, identity / "CR_CB_mapping_audit.tsv")

    def stage02(self):
        _, _, r1, r2 = self._paths()
        output = self.base / "02.umiTools" / f"{self.sample}.raw.reads_barcoded.fq.gz"
        barcode_table = output.parent / f"{self.sample}.reads.barcodes.tsv"
        barcode_stats = output.parent / f"{self.sample}.barcodes.stat.tsv"
        if _complete(output, barcode_table, barcode_stats) and not self.force:
            return
        cb_len, umi_len, _ = chemistry_parameters(
            str(self.library["three_or_five_prime"]), str(self.library["tenx_version"]))
        if self.dry_run:
            print(f"scatims embed {quote(r1)} {quote(r2)} -> {quote(output)} (CB={cb_len}, UMI={umi_len})")
        else:
            embed_barcode_umi(r1, r2, output, cb_len, umi_len, barcode_table, barcode_stats)
            _require(output, barcode_table, barcode_stats)

    def stage03(self):
        source = self.base / "02.umiTools" / f"{self.sample}.raw.reads_barcoded.fq.gz"
        out = self.base / "03.hg38"
        bam, unmapped_bam = out / f"{self.sample}_rna_hg38_align_sort.bam", out / f"{self.sample}_rna_hg38_unmap.bam"
        unmapped_fq = out / f"{self.sample}_rna_hg38_unmap.fq"
        flagstat = Path(str(bam) + ".flagstat.txt")
        if _complete(bam, unmapped_bam, unmapped_fq, flagstat) and not self.force:
            return
        if not self.dry_run:
            out.mkdir(parents=True, exist_ok=True)
        bwa, samtools, threads = self.tools.get("bwa", "bwa"), self.tools.get("samtools", "samtools"), self.compute.get("threads", 32)
        rg = f"@RG\\tID:{self.sample}\\tSM:{self.sample}"
        self._run("03.align", f"{quote(bwa)} mem -q -t {threads} -R {quote(rg)} {quote(self.resources['host_fasta_grch38'])} {quote(source)} | {quote(samtools)} view -bS -@ {threads} | {quote(samtools)} sort -@ {threads} -m {quote(self.compute.get('samtools_mem_per_thread', '1G'))} -o {quote(bam)}")
        self._run("03.extract", f"{quote(samtools)} view -b -@ {threads} -f 4 -o {quote(unmapped_bam)} {quote(bam)} && {quote(samtools)} fastq -@ {threads} {quote(unmapped_bam)} > {quote(unmapped_fq)} && {quote(samtools)} flagstat -@ {threads} {quote(bam)} > {quote(flagstat)}")
        if not self.dry_run:
            _require(bam, unmapped_bam, unmapped_fq, flagstat)

    def stage04(self):
        source = self.base / "03.hg38" / f"{self.sample}_rna_hg38_unmap.fq"
        out = self.base / "04.t2t"
        bam = out / f"{self.sample}_rna_hg38_unmap_t2t_align_sort.bam"
        unmapped = out / f"{self.sample}_rna_hg38_unmap_t2t_unmap.bam"
        flagstat = Path(str(unmapped) + ".flagstat.txt")
        if _complete(bam, unmapped, flagstat) and not self.force:
            return
        if not self.dry_run:
            out.mkdir(parents=True, exist_ok=True)
        bwa, samtools, threads = self.tools.get("bwa", "bwa"), self.tools.get("samtools", "samtools"), self.compute.get("threads", 32)
        rg = f"@RG\\tID:{self.sample}\\tSM:{self.sample}"
        self._run("04.align", f"{quote(bwa)} mem -q -t {threads} -R {quote(rg)} {quote(self.resources['host_fasta_t2t'])} {quote(source)} | {quote(samtools)} view -bS -@ {threads} | {quote(samtools)} sort -@ {threads} -m {quote(self.compute.get('samtools_mem_per_thread', '1G'))} -o {quote(bam)}")
        self._run("04.extract", f"{quote(samtools)} view -b -@ {threads} -f 4 -o {quote(unmapped)} {quote(bam)} && {quote(samtools)} flagstat -@ {threads} {quote(unmapped)} > {quote(flagstat)}")
        if not self.dry_run:
            _require(bam, unmapped, flagstat)

    def stage04b(self):
        if not self.cfg.get("workflow", {}).get("vector_filter", True):
            return
        source = self.base / "04.t2t" / f"{self.sample}_rna_hg38_unmap_t2t_unmap.bam"
        out = self.base / "04b.vector"
        fastq, bam_dir = out / "01.fastq", out / "02.bam"
        se = fastq / f"{self.sample}.t2t_unmap.SE.fq.gz"
        aligned = bam_dir / f"{self.sample}.t2t_unmap.vector.bam"
        unmapped = bam_dir / f"{self.sample}.t2t_unmap.vector_unmapped.bam"
        if _complete(se, aligned, unmapped) and not self.force:
            return
        if not self.dry_run:
            fastq.mkdir(parents=True, exist_ok=True)
            bam_dir.mkdir(parents=True, exist_ok=True)
        bwa, samtools, threads = self.tools.get("bwa", "bwa"), self.tools.get("samtools", "samtools"), self.compute.get("vector_threads", self.compute.get("threads", 32))
        self._run("04b.fastq", f"{quote(samtools)} sort -n -@ {threads} {quote(source)} | {quote(samtools)} fastq -@ {threads} - -0 {quote(se)} -1 {quote(fastq / (self.sample + '.R1.fq.gz'))} -2 {quote(fastq / (self.sample + '.R2.fq.gz'))} -s {quote(fastq / (self.sample + '.singles.fq.gz'))}")
        rg = f"@RG\\tID:{self.sample}\\tSM:{self.sample}"
        self._run("04b.align", f"{quote(bwa)} mem -q -t {threads} -R {quote(rg)} {quote(self.resources['vector_fasta'])} {quote(se)} | {quote(samtools)} view -b -@ {threads} -o {quote(aligned)} -")
        self._run("04b.extract", f"{quote(samtools)} view -b -@ {threads} -f 4 -o {quote(unmapped)} {quote(aligned)}")
        if not self.dry_run:
            _require(se, aligned, unmapped)
