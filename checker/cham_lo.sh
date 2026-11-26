#!/bin/bash
set -u
DS=$(readlink -f "$1"); N=${2:-4}; GB=${3:-6}
K=$(dirname "$(readlink -f "$0")"); P=~/miniforge3/envs/b2demo/bin/python; GPU=~/Vision_team_data_189/sonhh20/chay_gpu.sh
O=$DS.ket_qua; mkdir -p $O/khoa $O/xong $O/loi $O/log
GPU_CHAM=${GPU_CHAM:-auto}; LOI_DAU=${LOI_DAU:-10}
if [ "${CHO_BO_DUNG:-1}" = 1 ]; then
  co_bo_dung() { for p in /proc/[0-9]*; do [ "$p" = /proc/$$ ] && continue; grep -qa "botri_chung" $p/cmdline 2>/dev/null && return 0; done; return 1; }
  while co_bo_dung; do echo "$(date +%T) cho bo dung xong..." >> $O/nhat_ky.txt; sleep 60; done
fi
cong_nhan() {
  local cn=$1 loi=$2
  while true; do
    local con=$(df -BG --output=avail $K 2>/dev/null | tail -1 | tr -dc 0-9); con=${con:-0}
    if [ "$con" -lt "${DIA_MIN:-8}" ]; then echo "$(date +%T) $cn DUNG: dia con ${con} GB < ${DIA_MIN:-8} GB" >> $O/nhat_ky.txt; touch $O/DUNG_DIA; break; fi
    [ -e $O/DUNG ] && break
    local lay=""
    while read -r d v l; do
      [ -z "${d:-}" ] && continue; [[ "$d" == \#* ]] && continue
      local k=$(echo "$d" | md5sum | cut -c1-12)
      [ -e $O/xong/$k ] && continue
      [ -e $O/xong/$k.loi ] && [ "${CHAY_LAI_LOI:-0}" != 1 ] && continue
      if mkdir $O/khoa/$k 2>/dev/null; then lay="$d|${v:-dong_vat}|${l:-}|$k"; break; fi
    done < $DS
    [ -z "$lay" ] && break
    IFS='|' read -r d v l k <<< "$lay"
    echo "$(date +%T) $cn BAT_DAU $d" >> $O/nhat_ky.txt; t0=$(date +%s)
    ( cd $K && OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 NUMBA_NUM_THREADS=4 IGL_NUM_THREADS=4 PYTHONNOUSERSITE=1 \
      CHUNG=1 $GPU $GPU_CHAM $GB /usr/bin/time -v taskset -c $loi $P cham_mot.py "$d" --vat $v ${l:+--lop $l} ) > $O/log/$k.log 2>&1
    rc=$?
    echo "$(date +%T) $cn XONG $d rc=$rc $(( $(date +%s) - t0 ))s $(grep -o "Maximum resident set size (kbytes): [0-9]*" $O/log/$k.log | grep -o "[0-9]*$")KB" >> $O/nhat_ky.txt
    if [ $rc = 0 ]; then touch $O/xong/$k; else touch $O/loi/$k; rmdir $O/khoa/$k; touch $O/xong/$k.loi; fi
  done
}
( while ls $O/khoa >/dev/null 2>&1 && [ ! -e $O/HET ]; do
    echo "$(date +%T) gpu[$(nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv,noheader,nounits | tr '\n' ' ')] tai $(cut -d' ' -f1 /proc/loadavg) xong $(ls $O/xong | wc -l)" >> $O/giam_sat.txt; sleep 60
  done ) &
GS=$!
for i in $(seq 0 $((N - 1))); do cong_nhan cn$i "$((LOI_DAU + i * 4))-$((LOI_DAU + i * 4 + 3))" & done
wait $(jobs -p | grep -v "^$GS$")
touch $O/HET
first=1; : > $O/tong.csv
while read -r d v l; do
  [ -z "${d:-}" ] && continue; [[ "$d" == \#* ]] && continue
  f="$d/kiem3.csv"; [ -f "$f" ] || f="$K/cham_ra/$(basename "$(dirname "$d")")__$(basename "$d")/kiem3.csv"; [ -f "$f" ] || continue
  if [ $first = 1 ]; then cat "$f" >> $O/tong.csv; first=0; else tail -n +2 "$f" >> $O/tong.csv; fi
done < $DS
echo "$(date +%T) XONG_TAT_CA $(ls $O/xong | grep -vc loi) te bao, loi: $(ls $O/loi | wc -l)" >> $O/nhat_ky.txt
