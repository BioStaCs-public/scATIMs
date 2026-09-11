#!/usr/bin/env python3
"""ATIMs step of scATIMs.

Statistical decontamination of scMic assignments across taxonomic ranks.

【核心说明】：
1. 精简防线：由于上游引入 UMI 物理纠错，移除下游 Redundancy 冗余度质控。
2. 引入“双向平行法庭”：同时筛查 Cell-Enriched (真实感染) 与 Empty-Enriched (环境污染)。
3. 统计检验：全面采用 Fisher 检验、Mann-Whitney 检验与 Dispersion 检验进行严格双向判决。
4. P值校正：加入 FDR(Benjamini-Hochberg) 检验校正，增强显著性检验可靠性。移除了 Empty 侧离散度检验。
"""

import argparse
import gzip
import json
import traceback
from collections import defaultdict
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, chi2
from statsmodels.stats.multitest import multipletests # 新增：用于P值校正

# ================= 数据加载与清洗 =================

def open_maybe_gzip(path, mode='rt'):
    path = str(path)
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)


def _resolve_sample_path(sample_dir: Path, template: str, sample_id: str) -> Path:
    """Resolve a configurable sample-relative path template."""
    rendered = Path(template.format(sample_id=sample_id))
    return rendered if rendered.is_absolute() else sample_dir / rendered

def load_droplet_statistics(mapping_file: Path, whitelist_file: Path):
    all_valid_cbs = set()
    with open_maybe_gzip(mapping_file) as f:
        header = f.readline().strip().split("\t")
        cb_idx = header.index("CB_Corrected") if "CB_Corrected" in header else 1
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) > cb_idx:
                all_valid_cbs.add(parts[cb_idx])

    whitelist_barcodes = {line.strip() for line in open_maybe_gzip(whitelist_file) if line.strip()}
    true_cells = whitelist_barcodes.intersection(all_valid_cbs)

    total_cell = len(true_cells)
    total_empty = len(all_valid_cbs) - total_cell

    return total_cell, total_empty

def prepare_umi_data(umi_df, level):
    df = umi_df[umi_df['taxonomic_level'] == level].copy()
    df = df[df['taxon_name'].notna()]
    df = df[df['taxon_name'].str.lower() != 'unknown']
    df = df[df['score'] > 0]
    df = df.rename(columns={'barcode': 'cell_barcode', 'taxon_name': 'feature_name'})
    return df[['cell_barcode', 'UMI', 'feature_name', 'score']].copy()

# ================= 核心统计引擎 =================

def calculate_barcode_statistics(df_umi):
    if df_umi.empty: return {}
    stats_df = df_umi.groupby(['cell_barcode', 'feature_name'])['score'].sum().reset_index()
    bc_feat_score = defaultdict(dict)
    for _, row in stats_df.iterrows():
        bc_feat_score[row['cell_barcode']][row['feature_name']] = float(row['score'])
    return bc_feat_score

def analyze_features(df_cell_umi, df_empty_umi, total_cell, total_empty, level, min_evidence_score):
    all_features = sorted(set(df_cell_umi['feature_name'].unique()) | set(df_empty_umi['feature_name'].unique()))
    cell_score = calculate_barcode_statistics(df_cell_umi)
    empty_score = calculate_barcode_statistics(df_empty_umi)

    results = []
    pseudocount = 1e-4

    for feature in all_features:
        c_scores = [d[feature] for d in cell_score.values() if feature in d]
        e_scores = [d[feature] for d in empty_score.values() if feature in d]

        cells_with_mic = sum(1 for s in c_scores if s >= min_evidence_score)
        empty_with_mic = sum(1 for s in e_scores if s >= min_evidence_score)

        _, p_fisher_cell = fisher_exact(
            [[cells_with_mic, total_cell - cells_with_mic],
             [empty_with_mic, total_empty - empty_with_mic]], alternative='greater'
        )
        _, p_fisher_empty = fisher_exact(
            [[cells_with_mic, total_cell - cells_with_mic],
             [empty_with_mic, total_empty - empty_with_mic]], alternative='less'
        )

        cp = (cells_with_mic + pseudocount) / (total_cell + pseudocount)
        ep = (empty_with_mic + pseudocount) / (total_empty + pseudocount)
        enrichment_ratio_cell = cp / ep
        enrichment_ratio_empty = ep / cp

        full_c = c_scores + [0.0] * max(0, total_cell - len(c_scores))
        full_e = e_scores + [0.0] * max(0, total_empty - len(e_scores))

        try: _, p_mw_cell = mannwhitneyu(full_c, full_e, alternative='greater')
        except: p_mw_cell = np.nan

        try: _, p_mw_empty = mannwhitneyu(full_c, full_e, alternative='less')
        except: p_mw_empty = np.nan

        c_score_sum = sum(c_scores)
        e_score_sum = sum(e_scores)
        cm = c_score_sum / total_cell if total_cell else 0
        em = e_score_sum / total_empty if total_empty else 0
        log2fc_cell = np.log2((cm + pseudocount) / (em + pseudocount))
        log2fc_empty = np.log2((em + pseudocount) / (cm + pseudocount))

        pos_c_scores = [s for s in c_scores if s >= min_evidence_score]
        n_pos_c = len(pos_c_scores)
        if n_pos_c > 1:
            chi2_c = (n_pos_c - 1) * (np.var(pos_c_scores, ddof=1) / (np.mean(pos_c_scores) + 1e-8))
            p_dispersion_cell = 1.0 - chi2.cdf(chi2_c, df=n_pos_c - 1)
        else: p_dispersion_cell = 1.0

        # 修改点：移除了 p_dispersion_empty 的计算

        results.append({
            'taxonomic_level': level,
            'feature_name': feature,
            'total_cell': total_cell,
            'total_empty': total_empty,
            'cells_with_mic': cells_with_mic,
            'empty_with_mic': empty_with_mic,
            'cell_total_score': round(c_score_sum, 4),
            'empty_total_score': round(e_score_sum, 4),

            'log2fc_cell': round(log2fc_cell, 4),
            'enrichment_ratio_cell': round(enrichment_ratio_cell, 4),
            'p_value_fisher_cell': p_fisher_cell,
            'p_value_mw_cell': p_mw_cell,
            'p_value_dispersion_cell': p_dispersion_cell,

            'log2fc_empty': round(log2fc_empty, 4),
            'enrichment_ratio_empty': round(enrichment_ratio_empty, 4),
            'p_value_fisher_empty': p_fisher_empty,
            'p_value_mw_empty': p_mw_empty
            # 修改点：移除了 p_value_dispersion_empty 输出
        })

    df_results = pd.DataFrame(results)

    # 修改点：统一进行 FDR (Benjamini-Hochberg) 多重检验校正
    if not df_results.empty:
        # Cell 侧检验 P 值校正
        df_results['p_adj_fisher_cell'] = multipletests(df_results['p_value_fisher_cell'].fillna(1.0), method='fdr_bh')[1]
        df_results['p_adj_mw_cell'] = multipletests(df_results['p_value_mw_cell'].fillna(1.0), method='fdr_bh')[1]
        df_results['p_adj_dispersion_cell'] = multipletests(df_results['p_value_dispersion_cell'].fillna(1.0), method='fdr_bh')[1]

        # Empty 侧检验 P 值校正
        df_results['p_adj_fisher_empty'] = multipletests(df_results['p_value_fisher_empty'].fillna(1.0), method='fdr_bh')[1]
        df_results['p_adj_mw_empty'] = multipletests(df_results['p_value_mw_empty'].fillna(1.0), method='fdr_bh')[1]

    return df_results

# ================= 阈值宣判 (产生 2 份榜单) =================

def filter_microbes(df_results, thresholds):
    """法官宣判：仅输出符合 MW 检验在内的严苛条件结果 (使用 FDR 校正后 P 值)"""
    if df_results.empty:
        return df_results, df_results, thresholds

    t_fisher = thresholds['fisher_pval']
    t_mw = thresholds['mw_pval']
    t_enrich = thresholds['enrichment_ratio']
    t_log2fc = thresholds['log2fc']
    t_dispersion = thresholds['dispersion_pval']

    # --- 宣判基础防线准备 (使用校正后的 P 值 p_adj_*) ---
    cond_enrich_c = df_results['enrichment_ratio_cell'] > t_enrich
    cond_fc_c = df_results['log2fc_cell'] > t_log2fc
    cond_fisher_c = df_results['p_adj_fisher_cell'] < t_fisher
    cond_mw_c = df_results['p_adj_mw_cell'] < t_mw
    cond_disp_c = df_results['p_adj_dispersion_cell'] < t_dispersion

    cond_enrich_e = df_results['enrichment_ratio_empty'] > t_enrich
    cond_fc_e = df_results['log2fc_empty'] > t_log2fc
    cond_fisher_e = df_results['p_adj_fisher_empty'] < t_fisher
    cond_mw_e = df_results['p_adj_mw_empty'] < t_mw
    # 修改点：Empty 侧不再需要 cond_disp_e

    # --- 构建掩码 ---
    cell_mask = cond_fisher_c & cond_enrich_c & cond_mw_c & cond_fc_c & cond_disp_c
    empty_mask = cond_fisher_e & cond_enrich_e & cond_mw_e & cond_fc_e # 修改点：移除了 cond_disp_e

    df_cell = df_results[cell_mask].copy()
    df_empty = df_results[empty_mask].copy()

    used = thresholds.copy()
    # 标注现在使用的是校正后P值作为阈值
    used['p_value_type'] = 'FDR_adjusted (Benjamini-Hochberg)'

    return df_cell, df_empty, used

# ================= 调度与输出 =================

def process_sample_level(sample_id, scmic_dir, out_dir, total_cell, total_empty, level, thresholds, min_ev_score):
    # 定义输出文件路径
    cell_file = out_dir / f'{sample_id}_ATIMs_{level}_CellEnriched.tsv'
    empty_file = out_dir / f'{sample_id}_ATIMs_{level}_EmptyEnriched.tsv'
    stats_file = out_dir / f'{sample_id}_ATIMs_{level}_Statistics.tsv'
    param_file = out_dir / f'{sample_id}_ATIMs_{level}_Parameters.json'

    if (cell_file.exists() and empty_file.exists() and stats_file.exists() and param_file.exists() and
        stats_file.stat().st_size > 0):
        df_stats = pd.read_csv(stats_file, sep='\t')
        df_c = pd.read_csv(cell_file, sep='\t')
        df_e = pd.read_csv(empty_file, sep='\t')
        return len(df_stats), len(df_c), len(df_e)

    cell_umi_file = scmic_dir / f"{sample_id}_cell_UMI_assignment.tsv"
    empty_umi_file = scmic_dir / f"{sample_id}_empty_UMI_assignment.tsv"

    df_cell = pd.read_csv(cell_umi_file, sep='\t')
    df_empty = pd.read_csv(empty_umi_file, sep='\t')

    df_cell_clean = prepare_umi_data(df_cell, level)
    df_empty_clean = prepare_umi_data(df_empty, level)

    df_results = analyze_features(df_cell_clean, df_empty_clean, total_cell, total_empty, level, min_ev_score)
    df_c, df_e, t_used = filter_microbes(df_results, thresholds)

    t_used['min_evidence_score'] = min_ev_score
    # 结果落盘
    df_results.to_csv(stats_file, sep='\t', index=False)
    df_c.to_csv(cell_file, sep='\t', index=False)
    df_e.to_csv(empty_file, sep='\t', index=False)

    with open(param_file, 'w', encoding='utf-8') as f:
        json.dump(t_used, f, indent=2, ensure_ascii=False)

    return len(df_results), len(df_c), len(df_e)

def run_one_sample(
    root_dir,
    sample_id,
    levels,
    thresholds,
    min_ev_score,
    scmic_dir_name,
    out_dir_name,
    barcode_mapping_relative="01.starSolo_bam/CR_CB_mapping.tsv",
    barcode_whitelist_relative="01.starSolo/Solo.out/Gene/filtered_gz/barcodes.tsv.gz",
):
    try:
        sample_dir = root_dir / sample_id
        scmic_dir = sample_dir / scmic_dir_name
        out_dir = sample_dir / out_dir_name
        out_dir.mkdir(parents=True, exist_ok=True)

        mapping_file = _resolve_sample_path(sample_dir, barcode_mapping_relative, sample_id)
        filt_bc = _resolve_sample_path(sample_dir, barcode_whitelist_relative, sample_id)

        if not mapping_file.exists() or not filt_bc.exists():
            raise FileNotFoundError(f"Missing required barcode files in {sample_id}.")

        total_cell, total_empty = load_droplet_statistics(mapping_file, filt_bc)

        print(f"\n{'='*60}\n样本: {sample_id} | 真实细胞: {total_cell} | 空液滴: {total_empty}\n{'='*60}")

        for lvl in levels:
            c_all, c_enriched, e_enriched = process_sample_level(sample_id, scmic_dir, out_dir, total_cell, total_empty, lvl, thresholds, min_ev_score)
            print(f"  > 层级 [{lvl}]: 发现 {c_all} 种实体 | Cell富集: {c_enriched} | Empty富集: {e_enriched}")

        return True, "Success"
    except Exception as e:
        return False, f"{str(e)}\n{traceback.format_exc()}"

# ================= 命令行接口 =================

def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--root-dir', required=True)
    parser.add_argument('--sample-list', required=True)
    parser.add_argument('--tax-levels', default='species,genus,family')
    parser.add_argument('--scmic-dir-name', default='06.scMic_V10_EM')
    parser.add_argument('--output-dir-name', default='07.ATIMs_V10_Dual')
    parser.add_argument('--workers', type=int, default=12)
    parser.add_argument(
        '--barcode-whitelist-relative',
        default='01.starSolo/Solo.out/Gene/filtered_gz/barcodes.tsv.gz',
        help='Filtered barcode file relative to each sample directory',
    )
    parser.add_argument(
        '--barcode-mapping-relative',
        default='01.starSolo_bam/CR_CB_mapping.tsv',
        help='CR-to-CB mapping file relative to each sample directory',
    )
    parser.add_argument('--error-log-name', default='ATIMs_v10_Dual_error.log')

    parser.add_argument('--fisher_pval', type=float, default=0.05)
    parser.add_argument('--mw_pval', type=float, default=0.05)
    parser.add_argument('--enrichment_ratio', type=float, default=1.0)
    parser.add_argument('--log2fc', type=float, default=1.0)

    parser.add_argument(
        '--min_evidence_score',
        type=float,
        default=0.0,
        help='Paper-locked minimum per-droplet evidence score (default: 0.0)',
    )
    parser.add_argument('--dispersion_pval', type=float, default=0.05)

    args = parser.parse_args(argv)

    levels = ['species', 'genus', 'family', 'order', 'class', 'phylum'] if args.tax_levels == 'all' else [x.strip() for x in args.tax_levels.split(',')]

    thresholds = {
        'fisher_pval': args.fisher_pval,
        'mw_pval': args.mw_pval,
        'enrichment_ratio': args.enrichment_ratio,
        'log2fc': args.log2fc,
        'dispersion_pval': args.dispersion_pval
    }

    with open(args.sample_list, encoding='utf-8') as sample_handle:
        samples = [line.strip() for line in sample_handle if line.strip() and not line.lstrip().startswith('#')]
    root_path = Path(args.root_dir)

    log_file = root_path / args.error_log_name
    with open(log_file, "w") as f:
        f.write(f"--- Log started at {pd.Timestamp.now()} ---\n")

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        future_to_sample = {
            executor.submit(
                run_one_sample,
                root_path,
                s,
                levels,
                thresholds,
                args.min_evidence_score,
                args.scmic_dir_name,
                args.output_dir_name,
                args.barcode_mapping_relative,
                args.barcode_whitelist_relative,
            ): s
            for s in samples
        }

        for future in as_completed(future_to_sample):
            s = future_to_sample[future]
            ok, msg = future.result()
            if not ok:
                with open(log_file, "a") as lf:
                    lf.write(f"\n--- Sample: {s} ---\n{msg}\n")

if __name__ == '__main__':
    main()
