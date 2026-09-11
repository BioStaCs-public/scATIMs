#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Optional PathSeq prerequisite workflow for scATIMs
# ==============================================================================

# -----------------------------
# 默认参数
# -----------------------------

ROOT_DIR=""
SAMPLE_LIST=""

# PathSeq
JAVA_BIN="java"
GATK_JAR=""

HOST_DB=""
MICROBE_DB=""
INPUT_BAM_RELATIVE="04b.vector/02.bam/{sample_id}.t2t_unmap.vector_unmapped.bam"

PATHSEQ_DIR_NAME="05.pathseq_v10"

PATHSEQ_SPARK_MASTER="local[10]"
PATHSEQ_DRIVER_MEMORY="80g"
PATHSEQ_EXECUTOR_MEMORY="20g"
PATHSEQ_EXECUTOR_CORES="5"

PATHSEQ_XMX="64g"
PATHSEQ_XMS="32g"
PATHSEQ_MIN_CLIPPED_READ_LENGTH="60"

FORCE_RERUN="false"

# 并行控制参数
PARALLEL_JOBS=4
PARALLEL_METHOD="auto"

# -----------------------------
# 工具函数
# -----------------------------

usage() {
    cat <<EOF
用法:
  bash $(basename "$0") --root-dir ROOT_DIR --sample-list samples.txt [选项]

必需参数:
  --root-dir PATH                 样本根目录
  --sample-list FILE              每行一个 sample_id

PathSeq 参数:
  --java-bin PATH                 Java 可执行文件
  --gatk-jar FILE                 GATK jar 路径
  --host-db DIR                   host db 目录
  --microbe-db DIR                microbe db 目录
  --input-bam-relative TEMPLATE   样本目录下的输入BAM模板，可用 {sample_id}
  --pathseq-dir-name NAME         默认: 05.pathseq_v10
  --min-clipped-read-length INT   默认: 60
  --spark-master VALUE            默认: local[10]
  --driver-memory VALUE           默认: 80g
  --executor-memory VALUE         默认: 20g
  --executor-cores N              默认: 5
  --java-xmx VALUE                默认: 64g
  --java-xms VALUE                默认: 32g

并行控制:
  --parallel-jobs N               并行任务数，默认: 4
  --parallel-method METHOD        并行方式: auto/parallel/xargs/background

其他:
  --force-rerun                   强制重跑所有样本
  -h, --help                      显示帮助
EOF
}

log_info() {
    echo "[$(date '+%F %T')] [INFO] $*"
}

log_warning() {
    echo "[$(date '+%F %T')] [WARN] $*"
}

log_error() {
    echo "[$(date '+%F %T')] [ERROR] $*" >&2
}

log_success() {
    echo "[$(date '+%F %T')] [OK] $*"
}

require_file() {
    local f="$1"
    local msg="${2:-文件不存在}"
    [[ -f "$f" ]] || { log_error "$msg: $f"; exit 1; }
}

require_dir() {
    local d="$1"
    local msg="${2:-目录不存在}"
    [[ -d "$d" ]] || { log_error "$msg: $d"; exit 1; }
}

require_value() {
    local value="$1"
    local option="$2"
    [[ -n "$value" ]] || { log_error "缺少必需参数: $option"; usage; exit 2; }
}

require_executable() {
    local executable="$1"
    if [[ "$executable" == */* ]]; then
        [[ -x "$executable" ]] || { log_error "Java不可执行: $executable"; exit 1; }
    else
        command -v "$executable" >/dev/null 2>&1 || { log_error "找不到可执行程序: $executable"; exit 1; }
    fi
}

# -----------------------------
# 参数解析
# -----------------------------

while [[ $# -gt 0 ]]; do
    case "$1" in
        --root-dir) ROOT_DIR="$2"; shift 2 ;;
        --sample-list) SAMPLE_LIST="$2"; shift 2 ;;
        --java-bin) JAVA_BIN="$2"; shift 2 ;;
        --gatk-jar) GATK_JAR="$2"; shift 2 ;;
        --host-db) HOST_DB="$2"; shift 2 ;;
        --microbe-db) MICROBE_DB="$2"; shift 2 ;;
        --input-bam-relative) INPUT_BAM_RELATIVE="$2"; shift 2 ;;
        --pathseq-dir-name) PATHSEQ_DIR_NAME="$2"; shift 2 ;;
        --min-clipped-read-length) PATHSEQ_MIN_CLIPPED_READ_LENGTH="$2"; shift 2 ;;
        --spark-master) PATHSEQ_SPARK_MASTER="$2"; shift 2 ;;
        --driver-memory) PATHSEQ_DRIVER_MEMORY="$2"; shift 2 ;;
        --executor-memory) PATHSEQ_EXECUTOR_MEMORY="$2"; shift 2 ;;
        --executor-cores) PATHSEQ_EXECUTOR_CORES="$2"; shift 2 ;;
        --java-xmx) PATHSEQ_XMX="$2"; shift 2 ;;
        --java-xms) PATHSEQ_XMS="$2"; shift 2 ;;
        --force-rerun) FORCE_RERUN="true"; shift ;;
        --parallel-jobs) PARALLEL_JOBS="$2"; shift 2 ;;
        --parallel-method) PARALLEL_METHOD="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "未知参数: $1"; usage; exit 1 ;;
    esac
done

# -----------------------------
# 基础检查
# -----------------------------

require_value "$ROOT_DIR" --root-dir
require_value "$SAMPLE_LIST" --sample-list
require_value "$GATK_JAR" --gatk-jar
require_value "$HOST_DB" --host-db
require_value "$MICROBE_DB" --microbe-db
require_dir "$ROOT_DIR"
require_file "$SAMPLE_LIST"
require_executable "$JAVA_BIN"
require_file "$GATK_JAR"
require_dir "$HOST_DB"
require_dir "$MICROBE_DB"

# -----------------------------
# 读取样本列表
# -----------------------------

samples=()
while IFS= read -r line || [[ -n "$line" ]]; do
    sample=$(echo "$line" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    [[ -z "$sample" || "$sample" =~ ^# ]] && continue
    # 允许任何有效的样本名（不只是 SAMN 或数字开头）
    samples+=("$sample")
done < "$SAMPLE_LIST"

total_samples=${#samples[@]}
if [[ $total_samples -eq 0 ]]; then
    log_error "未找到任何有效样本"
    exit 1
fi

log_info "共 $total_samples 个样本待处理"

# -----------------------------
# 进度跟踪变量
# -----------------------------

# 创建临时目录存放进度文件
PROGRESS_DIR=$(mktemp -d)
trap "rm -rf $PROGRESS_DIR" EXIT

# 创建状态文件
SUCCESS_FILE="$PROGRESS_DIR/success"
FAILED_FILE="$PROGRESS_DIR/failed"
SKIPPED_FILE="$PROGRESS_DIR/skipped"
RUNNING_FILE="$PROGRESS_DIR/running"
> "$SUCCESS_FILE"
> "$FAILED_FILE"
> "$SKIPPED_FILE"
> "$RUNNING_FILE"

# 进度显示函数（带正在运行的样本）
print_progress() {
    local success=$(wc -l < "$SUCCESS_FILE" 2>/dev/null || echo 0)
    local failed=$(wc -l < "$FAILED_FILE" 2>/dev/null || echo 0)
    local skipped=$(wc -l < "$SKIPPED_FILE" 2>/dev/null || echo 0)
    local completed=$((success + failed + skipped))
    local percent=$((completed * 100 / total_samples))

    # 读取正在运行的样本
    local running_samples=""
    if [[ -f "$RUNNING_FILE" ]] && [[ -s "$RUNNING_FILE" ]]; then
        running_samples=$(cat "$RUNNING_FILE" | tr '\n' ',' | sed 's/,$//' | sed 's/,/, /g')
    fi

    # 清空当前行
    printf "\r\033[K"

    # 显示进度条
    printf "[%-50s] %d/%d (%d%%) | ✓:%d ✗:%d ⊘:%d" \
        "$(printf '#%.0s' $(seq 1 $((percent / 2))))" \
        "$completed" "$total_samples" "$percent" \
        "$success" "$failed" "$skipped"

    # 显示正在运行的样本
    if [[ -n "$running_samples" ]]; then
        printf " | ▶:%s" "$running_samples"
    fi
}

# 添加正在运行的样本
add_running() {
    local sample="$1"
    echo "$sample" >> "$RUNNING_FILE"
    print_progress
}

# 移除正在运行的样本
remove_running() {
    local sample="$1"
    grep -v "^$sample$" "$RUNNING_FILE" > "${RUNNING_FILE}.tmp" 2>/dev/null || true
    mv "${RUNNING_FILE}.tmp" "$RUNNING_FILE" 2>/dev/null || true
    print_progress
}

# 打印最终统计
print_final_stats() {
    echo ""  # 换行
    local success=$(wc -l < "$SUCCESS_FILE")
    local failed=$(wc -l < "$FAILED_FILE")
    local skipped=$(wc -l < "$SKIPPED_FILE")

    echo ""
    log_info "========== 执行完成 =========="
    log_info "✓ 成功: $success"
    log_info "✗ 失败: $failed"
    log_info "⊘ 跳过: $skipped"
    log_info "=============================="

    if [[ $failed -gt 0 ]]; then
        echo ""
        log_warning "失败的样本:"
        cat "$FAILED_FILE" | while read s; do echo "  - $s"; done
    fi
}

# -----------------------------
# PathSeq 运行函数（带进度记录和实时显示）
# -----------------------------

run_pathseq_one_sample() {
    local sample="$1"

    local sample_dir="${ROOT_DIR}/${sample}"
    local input_relative="${INPUT_BAM_RELATIVE//\{sample_id\}/$sample}"
    local input_bam="${sample_dir}/${input_relative}"
    local pathseq_dir="${sample_dir}/${PATHSEQ_DIR_NAME}"
    local pathseq_bam="${pathseq_dir}/${sample}_rna_output.pathseq.bam"
    local pathseq_txt="${pathseq_dir}/${sample}_rna_output.pathseq.txt"
    local sample_log="${pathseq_dir}/${sample}.pathseq.log"

    # 添加到运行列表
    add_running "$sample"

    # 检查输入文件
    if [[ ! -f "$input_bam" ]]; then
        remove_running "$sample"
        echo "$sample" >> "$FAILED_FILE"
        print_progress
        return 1
    fi

    mkdir -p "$pathseq_dir"

    # 检查是否需要跳过
    if [[ "$FORCE_RERUN" != "true" && -f "$pathseq_bam" && -f "$pathseq_txt" ]]; then
        remove_running "$sample"
        echo "$sample" >> "$SKIPPED_FILE"
        print_progress
        return 0
    fi

    # 运行 PathSeq（重定向所有输出到日志文件）
    set +e
    "$JAVA_BIN" \
      -Xmx"${PATHSEQ_XMX}" \
      -Xms"${PATHSEQ_XMS}" \
      -jar "$GATK_JAR" \
      PathSeqPipelineSpark \
      --spark-master "$PATHSEQ_SPARK_MASTER" \
      --conf spark.driver.memory="$PATHSEQ_DRIVER_MEMORY" \
      --conf spark.executor.memory="$PATHSEQ_EXECUTOR_MEMORY" \
      --conf spark.executor.cores="$PATHSEQ_EXECUTOR_CORES" \
      --input "$input_bam" \
      --min-clipped-read-length "$PATHSEQ_MIN_CLIPPED_READ_LENGTH" \
      --kmer-file "${HOST_DB}/GCA_000001405.15_GRCh38_no_alt_analysis_set.hss" \
      --filter-bwa-image "${HOST_DB}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.img" \
      --microbe-bwa-image "${MICROBE_DB}/pathseq_custom_V7_fixed.fa.img" \
      --microbe-dict "${MICROBE_DB}/pathseq_custom_V7_fixed.dict" \
      --taxonomy-file "${MICROBE_DB}/pathseq_custom_V7_fixed.db" \
      --output "$pathseq_bam" \
      --scores-output "$pathseq_txt" \
      --score-metrics "${pathseq_dir}/${sample}_rna_output.pathseq.score_metrics.txt" \
      --filter-metrics "${pathseq_dir}/${sample}_rna_output.pathseq.filter_metrics.txt" \
      --is-host-aligned false \
      --filter-duplicates false \
      --min-score-identity 0.7 \
      > "$sample_log" 2>&1

    local rc=$?
    set -e

    # 从运行列表中移除
    remove_running "$sample"

    # 记录结果
    if [[ $rc -eq 0 ]] && [[ -f "$pathseq_bam" ]] && [[ -f "$pathseq_txt" ]]; then
        echo "$sample" >> "$SUCCESS_FILE"
    else
        echo "$sample" >> "$FAILED_FILE"
    fi

    print_progress
    return $rc
}

# 导出函数和变量
export -f run_pathseq_one_sample add_running remove_running print_progress
export ROOT_DIR PATHSEQ_DIR_NAME PATHSEQ_SPARK_MASTER PATHSEQ_DRIVER_MEMORY
export PATHSEQ_EXECUTOR_MEMORY PATHSEQ_EXECUTOR_CORES PATHSEQ_XMX PATHSEQ_XMS
export PATHSEQ_MIN_CLIPPED_READ_LENGTH HOST_DB MICROBE_DB JAVA_BIN GATK_JAR
export INPUT_BAM_RELATIVE FORCE_RERUN PROGRESS_DIR SUCCESS_FILE FAILED_FILE SKIPPED_FILE RUNNING_FILE total_samples

# -----------------------------
# 选择并行方式
# -----------------------------

if [[ "$PARALLEL_METHOD" == "auto" ]]; then
    if command -v parallel &> /dev/null; then
        PARALLEL_METHOD="parallel"
    elif command -v xargs &> /dev/null; then
        PARALLEL_METHOD="xargs"
    else
        PARALLEL_METHOD="background"
    fi
fi

log_info "并行方式: $PARALLEL_METHOD, 任务数: $PARALLEL_JOBS"
echo ""
echo "开始处理... (Ctrl+C 可中断，已完成样本将自动跳过)"
echo ""

# 初始进度显示
print_progress

# -----------------------------
# 执行并行任务
# -----------------------------

case "$PARALLEL_METHOD" in
    parallel)
        printf "%s\n" "${samples[@]}" | \
            parallel -j "$PARALLEL_JOBS" --line-buffer 'run_pathseq_one_sample {}' 2>/dev/null
        ;;
    xargs)
        printf "%s\n" "${samples[@]}" | \
            xargs -I {} -P "$PARALLEL_JOBS" bash -c 'run_pathseq_one_sample "$@"' _ {} 2>/dev/null
        ;;
    background)
        declare -a pids=()
        declare -a running_samples_arr=()

        for sample in "${samples[@]}"; do
            while [[ ${#pids[@]} -ge $PARALLEL_JOBS ]]; do
                for i in "${!pids[@]}"; do
                    if ! kill -0 "${pids[$i]}" 2>/dev/null; then
                        wait "${pids[$i]}" 2>/dev/null
                        unset 'pids[i]'
                        unset 'running_samples_arr[i]'
                    fi
                done
                pids=("${pids[@]}")
                running_samples_arr=("${running_samples_arr[@]}")
                sleep 0.5
            done

            run_pathseq_one_sample "$sample" &
            pids+=($!)
            running_samples_arr+=("$sample")
        done

        for pid in "${pids[@]}"; do
            wait "$pid" 2>/dev/null
        done
        ;;
esac

# 打印最终统计
print_final_stats

# 如果有失败的任务，退出码为1
if [[ $(wc -l < "$FAILED_FILE") -gt 0 ]]; then
    exit 1
fi

exit 0
