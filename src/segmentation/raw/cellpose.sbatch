#!/bin/bash
#SBATCH -p gpu
#SBATCH --time=1-00:00:00
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=512G
#SBATCH --gres=gpu:1
#SBATCH -o logs/cellpose_%j.out
#SBATCH -e logs/cellpose_%j.err

set -euo pipefail

echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME}"
echo "Start: $(date)"

module load micromamba

INDIR="/projects/shared/ST_VISIUM_HD_Sandri/H&E"
OUTDIR="/home/rceccaroni/projects/HE-segmentation/cellpose"

BENCHDIR="${OUTDIR}/benchmarks/job_${SLURM_JOB_ID}"
TIME_LOG_DIR="${BENCHDIR}/time_logs"
GPU_LOG_DIR="${BENCHDIR}/gpu_logs"

mkdir -p "$OUTDIR" "$BENCHDIR" "$TIME_LOG_DIR" "$GPU_LOG_DIR" "$OUTDIR/segmentations"

SUMMARY_TSV="${BENCHDIR}/cellpose_benchmark_summary.tsv"

echo -e "job_id\timage\tstart_time\tend_time\telapsed_seconds\tstatus\ttime_log\tgpu_log" > "$SUMMARY_TSV"

{
  echo "Job ID: ${SLURM_JOB_ID}"
  echo "Node: ${SLURMD_NODENAME}"
  echo "CPUs: ${SLURM_CPUS_PER_TASK:-NA}"
  echo "Memory requested: ${SLURM_MEM_PER_NODE:-NA} MB"
  echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-NA}"
  echo "Start: $(date)"
  echo
  nvidia-smi
} > "${BENCHDIR}/job_info.txt"

for IMG in "$INDIR"/*.tif; do
  [ -e "$IMG" ] || continue

  BASENAME="$(basename "$IMG" .tif)"
  SAFE_NAME="$(echo "$BASENAME" | tr ' /' '__')"

  IMAGE_TIME_LOG="${TIME_LOG_DIR}/${SAFE_NAME}.time.txt"
  IMAGE_GPU_LOG="${GPU_LOG_DIR}/${SAFE_NAME}.gpu.csv"

  echo "Processing: $IMG"

  START_EPOCH=$(date +%s)
  START_ISO=$(date --iso-8601=seconds)

  # Monitor GPU
  nvidia-smi \
    --query-gpu=timestamp,index,name,utilization.gpu,utilization.memory,memory.used,memory.total,power.draw,temperature.gpu \
    --format=csv \
    -l 2 \
    > "$IMAGE_GPU_LOG" &

  GPU_MONITOR_PID=$!

  STATUS="OK"

  if micromamba run -n env-cellpose time -v -o "$IMAGE_TIME_LOG" \
    micromamba run -n env-cellpose python -m cellpose \
      --image_path "$IMG" \
      --pretrained_model cpsam \
      --use_gpu \
      --channel_axis -1 \
      --diameter 60 \
      --batch_size 128 \
      --flow_threshold 0.6 \
      --niter 2000 \
      --min_size 100 \
      --save_tif \
      --savedir "$OUTDIR/segmentations"
  then
    STATUS="OK"
  else
    STATUS="FAILED"
  fi

  # Stop monitor GPU
  if ps -p "$GPU_MONITOR_PID" > /dev/null 2>&1; then
    kill "$GPU_MONITOR_PID" || true
    wait "$GPU_MONITOR_PID" 2>/dev/null || true
  fi

  END_EPOCH=$(date +%s)
  END_ISO=$(date --iso-8601=seconds)
  ELAPSED=$((END_EPOCH - START_EPOCH))

  echo -e "${SLURM_JOB_ID}\t${IMG}\t${START_ISO}\t${END_ISO}\t${ELAPSED}\t${STATUS}\t${IMAGE_TIME_LOG}\t${IMAGE_GPU_LOG}" >> "$SUMMARY_TSV"

  echo "Done: $IMG"
  echo "Elapsed seconds: $ELAPSED"
  echo "Status: $STATUS"

  if [ "$STATUS" = "FAILED" ]; then
    echo "Error while processing $IMG"
    exit 1
  fi
done

if command -v sacct >/dev/null 2>&1; then
  sacct -j "$SLURM_JOB_ID" \
    --format=JobID,JobName,Partition,AllocCPUS,Elapsed,TotalCPU,MaxRSS,ReqMem,State,ExitCode \
    > "${BENCHDIR}/slurm_sacct.txt" || true
fi

DETAILED_TSV="${BENCHDIR}/cellpose_benchmark_detailed.tsv"

echo -e "image\telapsed_seconds\tstatus\tmax_ram_kb\tcpu_percent\tmax_gpu_util_percent\tmax_gpu_mem_used_mib\tmax_gpu_temp_c\tmax_gpu_power_w" > "$DETAILED_TSV"

tail -n +2 "$SUMMARY_TSV" | while IFS=$'\t' read -r job_id image start_time end_time elapsed_seconds status time_log gpu_log; do

  max_ram_kb="NA"
  cpu_percent="NA"
  max_gpu_util="NA"
  max_gpu_mem="NA"
  max_gpu_temp="NA"
  max_gpu_power="NA"

  if [ -f "$time_log" ]; then
    max_ram_kb=$(grep "Maximum resident set size" "$time_log" | awk '{print $NF}' || echo "NA")
    cpu_percent=$(grep "Percent of CPU this job got" "$time_log" | awk '{print $NF}' | tr -d '%' || echo "NA")
  fi

  if [ -f "$gpu_log" ]; then
    max_gpu_util=$(awk -F',' 'NR > 1 {
      gsub(/ %/, "", $4);
      if ($4+0 > max) max=$4+0
    } END {if (max=="") print "NA"; else print max}' "$gpu_log")

    max_gpu_mem=$(awk -F',' 'NR > 1 {
      gsub(/ MiB/, "", $6);
      if ($6+0 > max) max=$6+0
    } END {if (max=="") print "NA"; else print max}' "$gpu_log")

    max_gpu_power=$(awk -F',' 'NR > 1 {
      gsub(/ W/, "", $8);
      if ($8+0 > max) max=$8+0
    } END {if (max=="") print "NA"; else print max}' "$gpu_log")

    max_gpu_temp=$(awk -F',' 'NR > 1 {
      gsub(/ C/, "", $9);
      if ($9+0 > max) max=$9+0
    } END {if (max=="") print "NA"; else print max}' "$gpu_log")
  fi

  echo -e "${image}\t${elapsed_seconds}\t${status}\t${max_ram_kb}\t${cpu_percent}\t${max_gpu_util}\t${max_gpu_mem}\t${max_gpu_temp}\t${max_gpu_power}" >> "$DETAILED_TSV"

done

echo "Detailed benchmark table:"
echo "$DETAILED_TSV"

echo "End: $(date)"
echo "Benchmark files written to: $BENCHDIR"
