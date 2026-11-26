T_DG=$P/bo_dung_chung/tang_dg; N1G=/tmp/n1gpu; mkdir -p $N1G 2>/dev/null
tdg_song() { [ -e $T_DG/DUNG_LOI ] || [ -e $T_DG/XONG_TDG ] && return 1; [ -f $T_DG/dieu_phoi.pid ] && kill -0 $(cat $T_DG/dieu_phoi.pid) 2>/dev/null; }
u5_chua_khoi() { local id x; [ "$DRY" = 1 ] && return 1; tdg_song || return 1; [ -e $T_DG/DUNG_U5 ] && return 1
  while read -r id x; do case $id in u5_c_*|u5_x_*|u5_r_*) [ -e $T_DG/xong/$id ] || [ -e $T_DG/loi/$id ] || [ -d $T_DG/khoa/$id ] || return 0;; esac; done < $T_DG/viec.txt; return 1; }
tdg_con_viec() { local id x; tdg_song || return 1
  while read -r id x; do [ -e $T_DG/xong/$id ] || [ -e $T_DG/loi/$id ] || [ -d $T_DG/khoa/$id ] || return 0; done < $T_DG/viec.txt; return 1; }
dat_cho_g() { local dat=0 p_ gb_; [ -f /tmp/gpu$1.datcho ] && while read -r p_ gb_; do kill -0 $p_ 2>/dev/null && dat=$((dat + gb_)); done < /tmp/gpu$1.datcho; echo $dat; }
n_ta() { local f n=0; for f in $N1G/$1_*; do [ -e "$f" ] || continue; kill -0 ${f##*_} 2>/dev/null && n=$((n + 1)) || rm -f "$f"; done; echo $n; }
chon_gpu() {
  local gb=$1 g tr thu="1 0" n0 n1 cv=0
  u5_chua_khoi && return 1
  n0=$(n_ta 0); n1=$(n_ta 1); [ $n0 -gt 0 ] && [ $n1 = 0 ] && thu="0 1"; tdg_con_viec && cv=1
  for g in $thu; do
    if [ $cv = 1 ]; then if [ $g = 0 ]; then [ $n1 -gt 0 ] && continue; else [ $n0 -gt 0 ] && continue; fi; fi
    flock -n -s /tmp/gpu$g.khoa -c true 2>/dev/null || continue
    tr=$(( $(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits -i $g) / 1024 ))
    [ $(( $(dat_cho_g $g) + gb )) -le 22 ] && [ $tr -ge $((gb + 1)) ] && { echo $g; return 0; }
  done; return 1; }
gpu_san() {
  [ "$DRY" = 1 ] && return 0
  local f=$HN/.gpu_san_$1 t_; t_=$(stat -c %Y $f 2>/dev/null || echo 0)
  if [ $(( $(date +%s) - t_ )) -ge 20 ]; then chon_gpu $1 > $f.$$ 2>/dev/null; mv -f $f.$$ $f; fi; [ -s $f ]; }
do_dinh_gpu() {
  local p=$1 mx=0 m; mkdir -p $HN/bo_nho
  while kill -0 $p 2>/dev/null; do
    m=$(nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader,nounits 2>/dev/null | while IFS=', ' read -r q u; do
          [ "$(ps -o pgid= -p $q 2>/dev/null | tr -d ' ')" = "$p" ] && echo $u; done | awk '{s+=$1} END{print s+0}')
    [ "${m:-0}" -gt "$mx" ] && { mx=$m; echo $mx > $HN/bo_nho/$2; }; sleep 3; done; }
gpu_chay() {
  local gb=$1 lg=$2 g p t=0 rc w; shift 2
  [ "$DRY" = 1 ] && { chay $lg "$@"; return; }
  while true; do
    [ -e $DUNG_F ] && { echo "DUNG" > $HN/khong_dem_$ID_HT; return 1; }
    exec 6>/tmp/n1_gpu_vao.khoa; flock 6
    g=$(chon_gpu $gb) && break
    flock -u 6; exec 6>&-; sleep 30
  done
  CHUNG=1 GB_TOI_DA=$gb setsid $CG $g $gb /usr/bin/time -f "THOI_GIAN %e s %M KB" "$@" > $lg 2>&1 < /dev/null & p=$!
  while kill -0 $p 2>/dev/null && ! grep -aq "BAT_DAU_GPU$g" $lg; do sleep 3; t=$((t + 3))
    if [ $t -ge 45 ]; then kill -- -$p 2>/dev/null; wait $p 2>/dev/null; flock -u 6; exec 6>&-
      echo "khong vao duoc GPU$g trong 45 s (khong cho trong chay_gpu.sh)" > $HN/khong_dem_$ID_HT; log "GPU_KHONG_VAO $ID_HT g$g ($gb GB)"; return 1; fi; done
  touch $N1G/${g}_$p; flock -u 6; exec 6>&-
  log "GPU$g VAO $ID_HT ($gb GB, so truoc $(( $(dat_cho_g $g) - gb )) GB)"
  do_dinh_gpu $p $ID_HT & w=$!
  wait $p; rc=$?; kill $w 2>/dev/null; wait $w 2>/dev/null; rm -f $N1G/${g}_$p; return $rc; }
