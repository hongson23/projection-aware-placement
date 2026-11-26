#!/bin/bash
set -u
GOC=$HOME/Vision_team_data_189/sonhh20; P=$GOC/ch3_tebao_dep; PK=$P/ngoai_packmol; W=$PK/n1; CG=$GOC/chay_gpu.sh
B=$HOME/miniforge3/envs/b2demo/bin/python; PKB=$PK/cai/packmol-100.0.0-linux-x86_64/bin/packmol
NHAN_CN=${NHAN_CN:-"29 30 31"}; NHAN_N=${NHAN_N:-24-27}; DIA_MIN_GB=11
mkdir -p $W/khoa $W/xong $W/hong $W/log $W/thu $W/kq $W/bo_nho
log() { echo "$(date '+%m-%d %T') $*" >> $W/tien_do.txt; }
ENV1="OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONNOUSERSITE=1"
bd_args() {
  local TV=T SO="--so ti_the=600 --so lysosome=200"
  case $1 in D4) SO="--so ti_the=300 --so lysosome=100";; D4A2) TV=A; SO="--so ti_the=300 --so lysosome=100";; esac
  echo "--dac_ta dac_ta/v12/dong_vat_v12.json --seed $2 --N1 128 --K 4 --thu_vien $TV $SO --tang2 --khong_tu_do --eval_dung_tu_dau --out $3"; }
ds_ai() { [ $1 = D4A2 ] && echo "DS_AI=ai/khuon_chon/danh_sach_v2.json"; }
dir_cf() { echo $W/thu/$1_$2_k$3; }
xong() { [ -e $W/xong/$1 ]; }
hong2() { [ $(cat $W/hong/$1 2>/dev/null || echo 0) -ge 2 ]; }
chay() { local lg=$1; shift; /usr/bin/time -f "THOI_GIAN %e s %M KB" "$@" > $lg 2>&1; }
dia_du() { [ $(df -BG --output=avail $P | tail -1 | tr -dc 0-9) -ge $DIA_MIN_GB ]; }
ram_du() { [ $(awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo) -ge 40 ]; }
DRY=0; HN=$W; DUNG_F=$W/DUNG_PK; source $GOC/cellpack_ngoai/v122e3/gpu_n1.sh
gb_n() { [ -s $W/GB_N_$1 ] && { cat $W/GB_N_$1; return; }; case $1 in D5) echo 20;; *) echo 12;; esac; }
sinh_viec() { [ -s $W/viec.txt ] && return
  { for o in D4A2 D4 D5; do for s in $(seq 6201 6210); do for k in 0 1 2 3; do echo "K $o $s $k"; done; done; done
    for o in D4A2 D4 D5; do for s in $(seq 6201 6210); do echo "C $o $s"; echo "N $o $s"; done; done; } > $W/viec.txt; }
id_cua() { echo "$*" | tr ' ' '_'; }
san_sang() {
  set -- $@; [ -e $W/DUYET_PK ] || return 1
  case $1 in
    K) dia_du && ram_du || return 1; return 0;;
    C) local k; for k in 0 1 2 3; do xong K_$2_$3_$k || hong2 K_$2_$3_$k || [ -e $(dir_cf $2 $3 $k)/XUAT_HONG ] || return 1; done; return 0;;
    N) xong C_$2_$3 || { hong2 C_$2_$3 && return 2; return 1; }; [ -e $W/thu/$2_$3_CHON_HONG ] && return 2
       gpu_san $(gb_n $2) || return 1; dia_du || return 1; return 0;;
  esac; return 1; }
viec_K() { local o=$1 s=$2 k=$3 cpu=$4 D=$(dir_cf $1 $2 $3) h=0.25 t0 t1 rc gp
  [ -e $D/THU_H0125 ] && h=0.125
  if [ -s $D/pk/MD5_PK_OUT ] && [ -s $D/tu_the_pk.json ]; then log "K $o $s k$k: packmol da chay (MD5_PK_OUT) => chi chon lai"; viec_K_chon $o $s $k $cpu; return; fi
  if [ -s $D/pk/pk_out.xyz ] && grep -q '^RC ' $D/pk/pk.log 2>/dev/null; then log "K $o $s k$k: packmol da chay => chi doc + chon lai"; gp=$(grep '^RC ' $D/pk/pk.log | tail -1 | awk '{print $4}'); viec_K_doc $o $s $k $cpu $gp; return; fi
  rm -rf $D/pk $D/bd $D/*.obj $D/dau_vao_cellpack.json $D/tu_the_pk.json $D/chon_k$k.json; mkdir -p $D
  ( cd $P; chay $D/xuat_h$h.log env $ENV1 CUDA_VISIBLE_DEVICES="" $(ds_ai $o) taskset -c $cpu nice -n 10 $B $PK/cau_noi_pk_e3.py xuat --thu_muc $D --k $k --thiet_bi cpu --h_luoi_ngan $h -- $(bd_args $o $s $D/bd) )
  if ! grep -aq XUAT_XONG $D/xuat_h$h.log; then
    if grep -aq "không kín" $D/xuat_h$h.log; then
      if [ $h = 0.25 ]; then touch $D/THU_H0125; echo "khong kin h 0.25 => thu h 0.125" > $W/khong_dem_$ID_HT; log "XUAT_KHONG_KIN $o $s k$k => h 0.125"; return 1
      else grep -ao "AssertionError: [^\"]*" $D/xuat_h$h.log | head -1 > $D/XUAT_HONG; echo 2 > $W/hong/K_${o}_${s}_$k; log "XUAT_HONG $o $s k$k"; return 1; fi; fi
    return 1; fi
  rm -rf $D/bd
  ( cd $P; chay $D/dich.log env $ENV1 taskset -c $cpu nice -n 10 $B $PK/pk_dau_vao.py --thu_muc $D --seed_pk $((10 * s + k)) ) || return 1
  grep -aq PK_DAU_VAO_XONG $D/dich.log || return 1
  t0=$(date +%s.%N)
  ( cd $D/pk && env $ENV1 taskset -c $cpu nice -n 10 timeout -k 60 36000 $PKB -t 1 pk.inp > pk.log 2>&1 ); rc=$?
  t1=$(date +%s.%N); gp=$(awk "BEGIN{printf \"%.2f\", $t1 - $t0}"); echo "RC $rc GIAY $gp" >> $D/pk/pk.log
  [ -s $D/pk/pk_out.xyz ] || { log "PK_KHONG_RA $o $s k$k rc=$rc (ha tang)"; return 1; }
  [ $rc = 124 ] && log "PK_HET_GIO $o $s k$k (dung pk_out.xyz tot nhat da ghi)"
  viec_K_doc $o $s $k $cpu $gp; }
viec_K_doc() { local o=$1 s=$2 k=$3 cpu=$4 gp=$5 D=$(dir_cf $1 $2 $3)
  ( cd $P; chay $D/doc.log env $ENV1 taskset -c $cpu $B $PK/pk_doc.py --thu_muc $D --giay_pk $gp ) || return 1
  grep -aq PK_DOC_XONG $D/doc.log || return 1
  md5sum $D/pk/pk_out.xyz > $D/pk/MD5_PK_OUT; rm -f $D/pk/*.dat $D/pk/bp_*.xyz; gzip -f $D/pk/pk.inp
  viec_K_chon $o $s $k $cpu; }
viec_K_chon() { local o=$1 s=$2 k=$3 cpu=$4 D=$(dir_cf $1 $2 $3)
  ( cd $P; chay $D/chon.log env $ENV1 CUDA_VISIBLE_DEVICES="" $(ds_ai $o) taskset -c $cpu nice -n 10 $B $PK/cau_noi_pk_e3.py chon --thu_muc $D --k $k --ket_qua_cp $D/tu_the_pk.json --thiet_bi cpu -- $(bd_args $o $s $D/bd) ) || return 1
  rm -rf $D/bd; [ -s $D/chon_k$k.json ]; }
viec_C() { local o=$1 s=$2 ok=0 k
  for k in 0 1 2 3; do [ -s $(dir_cf $o $s $k)/chon_k$k.json ] && ok=$((ok + 1)); done
  if [ $ok = 0 ]; then touch $W/thu/${o}_${s}_CHON_HONG; log "CHON_HONG $o $s"; return 0; fi
  $B - $W/thu $o $s > $W/thu/${o}_${s}_chon.json <<'PY' || return 1
import json, sys, os
thu, o, s = sys.argv[1:]; ds = []
for k in range(4):
    f = f"{thu}/{o}_{s}_k{k}/chon_k{k}.json"
    if os.path.exists(f): ds.append(json.load(open(f)))
best = min(ds, key=lambda d: (d["so_khong_dat"], d["phat_cung_day_du"]))
tong = sum(d["giay_dat_cp"] + (d["giay_khoi_L"] or 0) + d["giay_cham_chon"] for d in ds)
print(json.dumps({"k_chon": best["k"], "K_hieu_luc": len(ds), "cac_k": [d["k"] for d in ds], "giay_dat_tong": tong, "chi_tiet": ds}))
PY
  log "CHON $o $s => $(cut -c1-80 $W/thu/${o}_${s}_chon.json)"; }
viec_N() { local o=$1 s=$2 k gd D GB=$(gb_n $1) O
  k=$($B -c "import json;print(json.load(open('$W/thu/${o}_${s}_chon.json'))['k_chon'])") || return 1
  gd=$($B -c "import json;print(json.load(open('$W/thu/${o}_${s}_chon.json'))['giay_dat_tong'])") || return 1
  D=$(dir_cf $o $s $k); O=$W/kq/${o}_pkchinh_${s}; rm -rf $D/bd $O
  ( cd $P; gpu_chay $GB $D/nhap.log env OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 PYTHONNOUSERSITE=1 PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True $(ds_ai $o) taskset -c $NHAN_N $B $PK/cau_noi_pk_e3.py nhap \
       --thu_muc $D --k $k --ket_qua_cp $D/tu_the_pk.json --thiet_bi cuda --giay_dat $gd -- $(bd_args $o $s $D/bd) ) || return 1
  grep -aq NHAP_XONG $D/nhap.log || return 1
  mv $D/bd $O && echo "k$k" > $O/NHAN_PACKMOL
  $B - $O/thong_tin.json $D/tu_the_pk.json > $O/KHOP_TU_THE.json <<'PY'
import json, sys, numpy as np
tt = json.load(open(sys.argv[1]))["tu_the"]; kq = json.load(open(sys.argv[2]))["dat"]; lech = 0.0; sai = []
for ten, v in kq.items():
    u = tt.get(ten)
    if u is None or u.get("t") is None: sai.append(ten); continue
    lech = max(lech, float(np.abs(np.array(u["t"]) - np.array(v["t"])).max()), float(np.abs(np.array(u["R"]) - np.array(v["R"])).max()))
print(json.dumps({"dat": lech <= 1e-6 and not sai, "lech_max": lech, "so_sai": len(sai), "sai": sai[:20], "so_bo_phan_packmol": len(kq)}))
PY
  log "N $o $s k$k khop_tu_the $(cut -c1-60 $O/KHOP_TU_THE.json)"; }
cong_nhan() { local w=$1 cpu=$2
  while true; do
    [ -e $W/DUNG_PK ] && { log "W$w dung (DUNG_PK)"; exit 0; }
    local lam=0 con=0 L id
    while read -r L; do [ -z "$L" ] && continue; id=$(id_cua $L)
      xong $id && continue; hong2 $id && continue
      san_sang $L; local r=$?; [ $r = 2 ] && continue; [ $r = 1 ] && { con=1; continue; }
      mkdir $W/khoa/$id 2>/dev/null || { local op=$(cat $W/khoa/$id/pid 2>/dev/null); if [ -n "$op" ] && ! kill -0 $op 2>/dev/null; then log "W$w nhan lai khoa $id"; else con=1; continue; fi; }
      echo $$ > $W/khoa/$id/pid; log "W$w BAT_DAU $id (nhan $cpu)"; local t0=$(date +%s) ok=1
      set -- $L; ID_HT=$id
      case $1 in K) viec_K $2 $3 $4 $cpu || ok=0;; C) viec_C $2 $3 || ok=0;; N) viec_N $2 $3 || ok=0;; esac
      if [ $ok = 1 ]; then touch $W/xong/$id; log "W$w XONG $id $(( $(date +%s) - t0 )) s"
      elif [ -e $W/khong_dem_$id ]; then log "W$w THU_LAI $id ($(head -c 100 $W/khong_dem_$id))"; rm -f $W/khong_dem_$id
      else local n=$(( $(cat $W/hong/$id 2>/dev/null || echo 0) + 1 )); echo $n > $W/hong/$id; log "W$w HONG $id (lan $n$( [ $n -ge 2 ] && echo ', bo' ))"; fi
      rm -rf $W/khoa/$id; lam=1; break
    done < $W/viec.txt
    if [ $lam = 0 ]; then [ $con = 0 ] && { log "W$w HET_VIEC"; touch $W/W_${w}_HET; exit 0; }; sleep 60; fi
  done; }
giam_sat() {
  local n_cn=$(wc -w <<< "$NHAN_CN") c
  while [ ! -e $W/DUNG_GS ]; do
    echo "$(date '+%m-%d %T') $(nvidia-smi --query-gpu=index,utilization.gpu,memory.used --format=csv,noheader | tr '\n' ' ') so0=$(dat_cho_g 0) so1=$(dat_cho_g 1) tai=$(cut -d' ' -f1 /proc/loadavg) RAM=$(awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo)G dia=$(df -BG --output=avail $P | tail -1 | tr -dc 0-9)G khoa=$(ls $W/khoa | wc -l) xong=$(ls $W/xong | wc -l)/$(wc -l < $W/viec.txt)" >> $W/giam_sat.log
    for p in $(pgrep -f "[p]ackmol -t 1 pk.inp"); do r=$(awk '/VmRSS/{print $2}' /proc/$p/status 2>/dev/null); [ -n "$r" ] && [ $r -gt $((8*1048576)) ] && { kill $p; log "CANH_RSS kill packmol $p"; }; done
    if [ -s $W/MO_RONG ]; then for c in $(cat $W/MO_RONG); do [ -e $W/cn_mo_$c.pid ] && kill -0 $(cat $W/cn_mo_$c.pid) 2>/dev/null && continue
        setsid nohup taskset -c $c bash $0 cong_nhan m$c $c > $W/log/cn_m$c.out 2>&1 < /dev/null & echo $! > $W/cn_mo_$c.pid; n_cn=$((n_cn + 1)); log "MO_RONG: them cong nhan nhan $c"; done; fi
    [ $(ls $W/W_*_HET 2>/dev/null | wc -l) -ge $n_cn ] && [ -z "$(ls $W/khoa)" ] && { log "GIAM_SAT: moi cong nhan het viec"; touch $W/XONG_HANG_PK; exit 0; }
    sleep 60; done; }
case ${1:-} in
  khoi_dong)
    exec 9>$W/khoi_dong.khoa; flock -n 9 || { echo "da co ban dang chay"; exit 1; }
    sinh_viec; log "KHOI_DONG cong nhan tren nhan $NHAN_CN (cho DUYET_PK)"
    while [ ! -e $W/DUYET_PK ] && [ ! -e $W/DUNG_PK ]; do sleep 60; done
    { md5sum $PK/cau_noi_pk_e3.py $PK/pk_dau_vao.py $PK/pk_doc.py $0 $P/bo_dung_chung/ma3/botri_chung_v122e3.py $PKB $GOC/cellpack_ngoai/v122e3/gpu_n1.sh; } > $W/moi_truong.txt 2>&1
    grep -q "^4131f86d.*cau_noi_pk_e3.py" $W/moi_truong.txt && grep -q "^a98211a1.*pk_dau_vao.py" $W/moi_truong.txt && grep -q "^289f6588.*pk_doc.py" $W/moi_truong.txt \
      && grep -q "^6e3beba8.*botri_chung_v122e3.py" $W/moi_truong.txt && grep -q "^55c2b577.*bin/packmol" $W/moi_truong.txt || { log "MD5_SAI (xem moi_truong.txt) => khong chay"; exit 1; }
    i=0; for c in $NHAN_CN; do setsid nohup taskset -c $c bash $0 cong_nhan $i $c > $W/log/cn_$i.out 2>&1 < /dev/null & echo $! > $W/cn_$i.pid; i=$((i + 1)); done
    setsid nohup taskset -c 31 nice -n 19 bash $0 giam_sat > $W/log/gs.out 2>&1 < /dev/null & echo $! > $W/gs.pid;;
  cong_nhan) cong_nhan $2 $3;;
  giam_sat) giam_sat;;
  *) echo "dung: $0 khoi_dong | cong_nhan <w> <nhan> | giam_sat"; exit 1;;
esac
