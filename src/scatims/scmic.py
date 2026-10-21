#!/usr/bin/env python3
"""scMic step of scATIMs.

Single-cell microbial signal assignment with barcode correction, directional
UMI collapse and anchor-weighted allocation across taxonomic ranks.

【核心升级说明 - 终极版】：
1. [新增] 宿主级 Barcode 纠错：引入 CR->CB 映射表，实时纠正测序错误，保证细胞物理身份一致性。
2. [新增] 病原级 UMI 纠错：引入 UMI-tools 的 Directional 网络算法 (大鱼吃小鱼)。在同一物理断点
   (Pattern) 内，利用丰度和汉明距离吸收测序错误的 UMI 子代，防止丰度虚高。
3. 分区锚点加权分配：分别在 cell 与 empty 数据流中扫描无争议 UMI
   建立锚点，再分配多重比对。
4. 全链路质控统计 (Metrics Tracker)：记录纠错前后的变化情况，生成评估报告。
"""

import argparse
import gzip
import csv
import traceback
import json
from pathlib import Path
from collections import defaultdict, Counter
from typing import Dict, List, Set, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
import pysam

# ================= 全局配置 =================
TAX_LEVELS = ["species", "genus", "family", "order", "class", "phylum"]

# ================= 工具与预处理函数 =================

def open_maybe_gzip(path: Path, mode: str = "rt"):
    """智能判断是否为 gzip 文件并打开"""
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)

def parse_qname(qname: str) -> Tuple[str, str]:
    """从 BAM query_name 提取原始 Barcode (CR) 和 原始 UMI"""
    parts = qname.split("_")
    return (parts[-2], parts[-1]) if len(parts) >= 3 else ("", "")

def load_barcode_mapping(path: Path) -> Dict[str, str]:
    """
    [新增] 加载 CR 到 CB 的映射表
    文件格式：两列，CR_Raw 和 CB_Corrected
    """
    mapping = {}
    if not path.exists():
        return mapping

    # 逐行读取，避免大文件一次性 load 内存爆炸
    with open_maybe_gzip(path) as f:
        header = f.readline().strip().split("\t")
        cr_idx = header.index("CR_Raw") if "CR_Raw" in header else 0
        cb_idx = header.index("CB_Corrected") if "CB_Corrected" in header else 1

        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                mapping[parts[cr_idx]] = parts[cb_idx]
    return mapping

def normalize_tax_name(x: object) -> str:
    s = str(x).strip() if x is not None else ""
    if not s or s.lower() in {"nan", "none", "unknown", "unclassified"}:
        return "unknown"
    return " ".join(s.replace("_", " ").split())

def canonicalize_name(name: str, level: str) -> str:
    name = normalize_tax_name(name)
    if name == "unknown": return "unknown"
    parts = name.split()
    if level == "species":
        if len(parts) >= 2:
            if parts[0] == "Candidatus": return " ".join(parts[:3])
            return " ".join(parts[:2])
    elif level == "genus":
        if parts[0] == "Candidatus" and len(parts) >= 2: return " ".join(parts[:2])
        return parts[0]
    return name

def load_taxonomy_hierarchy(path: Path) -> Dict[str, Dict[str, str]]:
    tax_map = {}
    df = pd.read_csv(path, sep="\t", dtype=str)
    col_mapping = {
        "species": "species_scientific_name", "genus": "genus_scientific_name",
        "family": "family_scientific_name", "order": "order_scientific_name",
        "class": "class_scientific_name", "phylum": "phylum_scientific_name"
    }
    for _, row in df.iterrows():
        tid = str(row["tax_id"])
        level_dict = {}
        for lvl in TAX_LEVELS:
            col_name = col_mapping.get(lvl)
            raw_name = row[col_name] if col_name in row else "unknown"
            level_dict[lvl] = canonicalize_name(raw_name, lvl)
        tax_map[tid] = level_dict
    return tax_map

# ================= UMI 物理纠错算法 (Directional) =================

def get_hamming_distance(s1: str, s2: str) -> int:
    """计算汉明距离"""
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def umi_tools_directional_collapse(umis: list) -> dict:
    """
    [新增核心] 等效实现 UMI-tools 的 Directional (有向图网络) 纠错算法
    输入: 一个 Pattern 物理断点下的所有 UMI 列表 (含重复)
    输出: {原始UMI: 纠正后的UMI} 映射字典
    """
    umi_counts = Counter(umis)

    # 1. 按照频数降序排列 (确定吞噬优先级)
    sorted_umis = sorted(umi_counts.keys(), key=lambda x: umi_counts[x], reverse=True)

    # 初始化：每个 UMI 默认映射给自己
    correction_map = {umi: umi for umi in sorted_umis}

    # 2. 贪心判定逻辑
    for i, target_umi in enumerate(sorted_umis):
        target_count = umi_counts[target_umi]

        for j in range(i + 1, len(sorted_umis)):
            query_umi = sorted_umis[j]
            query_count = umi_counts[query_umi]

            # 【保险栓】UMI-tools 定向网络的丰度判定数学界限：
            # 母代丰度必须大于等于子代丰度的两倍减一
            if target_count >= (2 * query_count) - 1:
                # 只吸收汉明距离为 1 的 1-bp 突变
                if get_hamming_distance(target_umi, query_umi) == 1:
                    # 如果这只“小鱼”还没有被吃掉，就被当前的高频“大鱼”吃掉
                    if correction_map[query_umi] == query_umi:
                        correction_map[query_umi] = target_umi

    # 3. 路径压缩：防止偶尔发生的链式吸收 (A吸收B，B吸收C，确保C直接指向A)
    for k in correction_map:
        root = k
        while correction_map[root] != root:
            root = correction_map[root]
        correction_map[k] = root

    return correction_map

# ================= 核心 EM 算法逻辑 =================

def calculate_global_anchors(read_df: pd.DataFrame) -> Dict[str, Dict[str, float]]:
    """【第一轮扫描】：只寻找无争议的 Unique 物理片段建立权威分"""
    anchors = {lvl: defaultdict(float) for lvl in TAX_LEVELS}
    grouped = read_df.groupby(["barcode", "UMI"], sort=False)

    for (bc, umi), g in grouped:
        tmp = g.copy()
        # 注意：这里继续使用 pattern 去重，因为前面 UMI 虽然纠正了，但可能还有完全重复的 PCR 冗余行
        tmp["pattern_key"] = (tmp["ref_name"].astype(str) + "_" + tmp["ref_start"].astype(str) + "_" + tmp["strand"].astype(str) + "_" + tmp["cigar"].astype(str))
        pattern_df = tmp.drop_duplicates("pattern_key")
        n_patterns = len(pattern_df)
        if n_patterns == 0: continue

        base_weight = 1.0 / n_patterns
        for _, row in pattern_df.iterrows():
            for lvl in TAX_LEVELS:
                cands = [x.strip() for x in str(row[f"cand_{lvl}"]).split(",") if x.strip() and x.strip() != "unknown"]
                if len(cands) == 1:
                    anchors[lvl][cands[0]] += base_weight

    return anchors

def compute_multilevel_support_pattern(sub_df: pd.DataFrame, global_anchors: Dict[str, Dict[str, float]]):
    """【第二轮扫描】：智能加权分配"""
    tmp = sub_df.copy()
    tmp["pattern_key"] = (tmp["ref_name"].astype(str) + "_" + tmp["ref_start"].astype(str) + "_" + tmp["strand"].astype(str) + "_" + tmp["cigar"].astype(str))
    pattern_df = tmp.drop_duplicates("pattern_key")
    n_patterns = len(pattern_df)

    supports = {lvl: defaultdict(float) for lvl in TAX_LEVELS}
    if n_patterns == 0:
        return supports, len(sub_df), 0, 0.0

    base_weight = 1.0 / n_patterns

    for _, row in pattern_df.iterrows():
        for lvl in TAX_LEVELS:
            cands = [x.strip() for x in str(row[f"cand_{lvl}"]).split(",") if x.strip() and x.strip() != "unknown"]

            if len(cands) == 1:
                supports[lvl][cands[0]] += base_weight
            elif len(cands) > 1:
                anchor_scores = [global_anchors[lvl].get(c, 0.0) for c in cands]
                total_anchor = sum(anchor_scores)

                if total_anchor > 0:
                    for c, a_score in zip(cands, anchor_scores):
                        supports[lvl][c] += base_weight * (a_score / total_anchor)
                else:
                    supports[lvl]["unknown"] += base_weight
            else:
                supports[lvl]["unknown"] += base_weight

    raw_reads = len(sub_df)
    ratio = raw_reads / n_patterns if n_patterns > 0 else 0.0

    return supports, raw_reads, n_patterns, ratio

# ================= 数据流处理与探针记录 =================

def build_read_table(sample_id, source_bam, whitelist, hits, tax_map, out_file, mode, barcode_mapping):
    base_cols = ["sample_id", "barcode_raw", "barcode", "UMI", "read_name", "entity_type",
                 "ref_name", "ref_start", "strand", "cigar"]
    cand_cols = [f"cand_{lvl}" for lvl in TAX_LEVELS]
    columns = base_cols + cand_cols

    # 统计探针
    stats = {
        "total_reads_parsed": 0,
        "reads_mapped_to_CB": 0,
        "reads_discarded_no_CB": 0,
        "reads_kept_for_entity": 0
    }

    with pysam.AlignmentFile(source_bam, "rb") as bam, open(out_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(columns)
        for r in bam:
            bc_raw, umi = parse_qname(r.query_name)
            if not bc_raw: continue

            stats["total_reads_parsed"] += 1

            # [升级逻辑] 实时判定并转换 Barcode
            if bc_raw not in barcode_mapping:
                stats["reads_discarded_no_CB"] += 1
                continue

            bc_corrected = barcode_mapping[bc_raw]
            stats["reads_mapped_to_CB"] += 1

            # 使用校正后的 CB 判断归属
            if mode == "cell" and bc_corrected not in whitelist: continue
            if mode == "empty" and bc_corrected in whitelist: continue

            if r.query_name not in hits: continue

            stats["reads_kept_for_entity"] += 1

            tids = hits[r.query_name]["ids"].split(",")
            cand_strings = []
            for lvl in TAX_LEVELS:
                names = {tax_map.get(t, {}).get(lvl, "unknown") for t in tids}
                cand_strings.append(",".join(sorted({n for n in names if n != "unknown"})))

            ref_start = r.reference_start if r.reference_start is not None else -1
            strand = "-" if r.is_reverse else "+"
            cigar = r.cigarstring if r.cigarstring is not None else "*"

            # 记录时保留 bc_raw 以便溯源，但下游核心聚类全部使用 bc_corrected (保存为 barcode)
            writer.writerow([sample_id, bc_raw, bc_corrected, umi, r.query_name, mode,
                             r.reference_name or "*", ref_start, strand, cigar] + cand_strings)

    return pd.read_csv(out_file, sep="\t"), stats

def build_umi_table(read_df, out_file):
    if read_df.empty:
        return pd.DataFrame(), {}

    # [升级逻辑] 物理位置聚类 + Directional UMI 纠错
    read_df["pattern_key"] = (
        read_df["ref_name"].astype(str) + "_" +
        read_df["ref_start"].astype(str) + "_" +
        read_df["strand"].astype(str) + "_" +
        read_df["cigar"].astype(str)
    )

    corrected_umis_list = []
    correction_count = 0  # 有多少条 Read 的 UMI 被强行改变了
    unique_umis_before = 0
    unique_umis_after = 0

    grouped = read_df.groupby(["barcode", "pattern_key"], sort=False)
    for (bc, pat), g in grouped:
        raw_umis = g["UMI"].tolist()
        unique_before = len(set(raw_umis))
        unique_umis_before += unique_before

        if unique_before == 1:
            corrected_umis_list.extend(raw_umis)
            unique_umis_after += 1
        else:
            correction_map = umi_tools_directional_collapse(raw_umis)
            for u in raw_umis:
                c_u = correction_map[u]
                corrected_umis_list.append(c_u)
                if u != c_u:
                    correction_count += 1
            unique_umis_after += len(set(correction_map.values()))

    # 替换 Dataframe 中的 UMI 为纠错后的 UMI，后续所有逻辑均基于这个
    read_df["UMI_raw"] = read_df["UMI"]
    read_df["UMI"] = corrected_umis_list

    umi_stats = {
        "total_reads_in_umi_collapse": len(read_df),
        "reads_umi_corrected_count": correction_count,
        "unique_umis_before_collapse": unique_umis_before,
        "unique_umis_after_collapse": unique_umis_after,
        "collapse_reduction_rate(%)": round((1 - unique_umis_after/unique_umis_before)*100, 2) if unique_umis_before > 0 else 0
    }

    rows = []
    global_anchors = calculate_global_anchors(read_df)

    # 注意：此时的 grouped 需要基于最新的纠错后 UMI 进行重组
    grouped_final = read_df.groupby(["barcode", "UMI"], sort=False)
    for (bc, umi), g in grouped_final:
        supports, r_cnt, p_cnt, red = compute_multilevel_support_pattern(g, global_anchors)

        for lvl in TAX_LEVELS:
            if not supports[lvl]:
                rows.append({
                    "barcode": bc, "UMI": umi, "taxonomic_level": lvl,
                    "taxon_name": "unknown", "score": 1.0,
                    "reads": r_cnt, "patterns": p_cnt, "redundancy": round(red, 4)
                })
                continue

            for taxon_name, score in supports[lvl].items():
                if score > 0:
                    rows.append({
                        "barcode": bc, "UMI": umi, "taxonomic_level": lvl,
                        "taxon_name": taxon_name, "score": round(score, 6),
                        "reads": r_cnt, "patterns": p_cnt, "redundancy": round(red, 4)
                    })

    if not rows: return pd.DataFrame(), umi_stats
    cols = ["barcode", "UMI", "taxonomic_level", "taxon_name", "score", "reads", "patterns", "redundancy"]
    df = pd.DataFrame(rows, columns=cols)
    df.to_csv(out_file, sep="\t", index=False)

    return df, umi_stats

def generate_resolution_summary(umi_df: pd.DataFrame, sample_id: str, out_file: Path):
    if umi_df.empty: return
    stats = []
    for lvl in TAX_LEVELS:
        lvl_df = umi_df[umi_df["taxonomic_level"] == lvl]
        if lvl_df.empty: continue

        total_score = lvl_df["score"].sum()
        total_umis = int(round(total_score))

        unknown_mask = lvl_df["taxon_name"].str.lower() == "unknown"
        unknown_score = lvl_df.loc[unknown_mask, "score"].sum()

        assigned_score = total_score - unknown_score
        assignment_rate = (assigned_score / total_score) * 100 if total_score > 0 else 0

        stats.append({
            "sample_id": sample_id, "taxonomic_level": lvl, "total_umis": total_umis,
            "assigned_score": round(assigned_score, 2), "unknown_score": round(unknown_score, 2),
            "assignment_rate(%)": round(assignment_rate, 2)
        })
    pd.DataFrame(stats).to_csv(out_file, sep="\t", index=False)

def _resolve_sample_path(sample_dir: Path, template: str, sample_id: str) -> Path:
    """Resolve a configurable sample-relative path template."""
    rendered = Path(template.format(sample_id=sample_id))
    return rendered if rendered.is_absolute() else sample_dir / rendered


def process_sample(
    root_dir,
    sample_id,
    taxonomy_hierarchy,
    output_subdir,
    pathseq_bam_relative="05.pathseq_v10/{sample_id}_rna_output.pathseq.bam",
    barcode_whitelist_relative="01.starSolo/Solo.out/Gene/filtered_gz/barcodes.tsv.gz",
    barcode_mapping_relative="01.starSolo_bam/CR_CB_mapping.tsv",
):
    try:
        s_dir = root_dir / sample_id
        out_dir = s_dir / output_subdir

        cell_umi_file = out_dir / f"{sample_id}_cell_UMI_assignment.tsv"
        empty_umi_file = out_dir / f"{sample_id}_empty_UMI_assignment.tsv"
        summary_file = out_dir / f"{sample_id}_cell_resolution_summary.tsv"
        metrics_file = out_dir / f"{sample_id}_correction_metrics.json"

        if (cell_umi_file.exists() and empty_umi_file.exists() and summary_file.exists() and metrics_file.exists() and
            cell_umi_file.stat().st_size > 0 and empty_umi_file.stat().st_size > 0):
            print(f"[SKIP] {sample_id} - 结果已存在，跳过处理")
            return True, "Skipped - already processed"

        bam_path = _resolve_sample_path(s_dir, pathseq_bam_relative, sample_id)
        whitelist_path = _resolve_sample_path(s_dir, barcode_whitelist_relative, sample_id)
        mapping_path = _resolve_sample_path(s_dir, barcode_mapping_relative, sample_id)

        if not bam_path.exists(): raise FileNotFoundError(f"Missing BAM: {bam_path}")
        if not whitelist_path.exists(): raise FileNotFoundError(f"Missing Whitelist: {whitelist_path}")
        if not mapping_path.exists(): raise FileNotFoundError(f"Missing Barcode Mapping: {mapping_path}")

        # 加载所有参照表
        tax_map = load_taxonomy_hierarchy(taxonomy_hierarchy)
        barcode_mapping = load_barcode_mapping(mapping_path)
        whitelist = {line.strip() for line in open_maybe_gzip(whitelist_path)}

        # 预先获取命中微生物的 query_names 字典
        hits = {}
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for r in bam:
                if r.has_tag("YP"): hits[r.query_name] = {"ids": r.get_tag("YP")}

        out_dir.mkdir(parents=True, exist_ok=True)

        cell_umi_df = None
        global_metrics = {"sample_id": sample_id}

        for mode in ["cell", "empty"]:
            # 构建 Read 表并返回 Barcode 矫正统计
            r_df, bc_stats = build_read_table(
                sample_id, bam_path, whitelist, hits, tax_map,
                out_dir/f"{sample_id}_{mode}_read_table.tsv", mode, barcode_mapping
            )
            # 构建 UMI 表并返回 UMI 去重纠错统计
            umi_df, umi_stats = build_umi_table(r_df, out_dir/f"{sample_id}_{mode}_UMI_assignment.tsv")

            global_metrics[f"{mode}_barcode_stats"] = bc_stats
            global_metrics[f"{mode}_umi_stats"] = umi_stats

            if mode == "cell": cell_umi_df = umi_df

        if cell_umi_df is not None:
            generate_resolution_summary(cell_umi_df, sample_id, out_dir/f"{sample_id}_cell_resolution_summary.tsv")

        # 落盘验证报告 (便于后续发文章和评估)
        with open(metrics_file, "w", encoding="utf-8") as mf:
            json.dump(global_metrics, mf, indent=4, ensure_ascii=False)

        return True, "Success"

    except Exception as e:
        error_msg = traceback.format_exc()
        return False, f"{str(e)}\n{error_msg}"

# ================= 主入口 =================

def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--root-dir", required=True)
    parser.add_argument("--sample-list", required=True)
    parser.add_argument("--taxonomy-hierarchy", required=True)
    parser.add_argument("--output-subdir", default="06.scMic_6_V10_EM_Dis", help="输出文件夹名称")
    parser.add_argument("--workers", type=int, default=1, help="并行样本数")
    parser.add_argument(
        "--pathseq-bam-relative",
        default="05.pathseq_v10/{sample_id}_rna_output.pathseq.bam",
        help="PathSeq BAM path relative to each sample directory; supports {sample_id}",
    )
    parser.add_argument(
        "--barcode-whitelist-relative",
        default="01.starSolo/Solo.out/Gene/filtered_gz/barcodes.tsv.gz",
        help="Filtered barcode file relative to each sample directory",
    )
    parser.add_argument(
        "--barcode-mapping-relative",
        default="01.starSolo_bam/CR_CB_mapping.tsv",
        help="CR-to-CB mapping file relative to each sample directory",
    )
    parser.add_argument("--error-log-name", default="scMic_multilevel_error.log")
    args = parser.parse_args(argv)

    with open(args.sample_list, encoding="utf-8") as sample_handle:
        samples = [line.strip() for line in sample_handle if line.strip() and not line.lstrip().startswith("#")]
    root_path, tax_path = Path(args.root_dir), Path(args.taxonomy_hierarchy)

    log_file = root_path / args.error_log_name
    with open(log_file, "w") as f:
        f.write(f"--- Log started at {pd.Timestamp.now()} ---\n")

    print(f"[INFO] 启动高级定管线 (CR->CB + UMI-Directional + EM): {args.workers} 进程，目标样本: {len(samples)}")

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        future_to_sample = {
            executor.submit(
                process_sample,
                root_path,
                s,
                tax_path,
                args.output_subdir,
                args.pathseq_bam_relative,
                args.barcode_whitelist_relative,
                args.barcode_mapping_relative,
            ): s
            for s in samples
        }

        for future in as_completed(future_to_sample):
            s = future_to_sample[future]
            try:
                ok, msg = future.result()
                if ok:
                    print(f"[ OK ] {s}")
                else:
                    print(f"[FAIL] {s} - 报错详情见 {log_file.name}")
                    with open(log_file, "a") as lf:
                        lf.write(f"--- Sample: {s} ---\n{msg}\n\n")
            except Exception as exc:
                print(f"[CRASH] {s}: {exc}")


if __name__ == "__main__":
    main()
