#!/bin/bash
set -u
DRY=${DRY:-0}
GOC=${GOC:-$HOME/Vision_team_data_189/sonhh20}
P=$GOC/ch3_tebao_dep; C=$GOC/cellpack_ngoai; G=$C/cellpack_gpu; V=$C/v122e3; H=$V/hang
XONG=$V/DUYET_N1
E=$GOC/envs/cellpack/bin/python; B=$HOME/miniforge3/envs/b2demo/bin/python; CG=$GOC/chay_gpu.sh
TS=$HOME/miniforge3/envs/b2demo/lib/python3.11/site-packages
KQ=$P/bo_dung_chung/tang_dg/kq
[ "$DRY" = 1 ] && E=/usr/bin/python3
SO_CN=${SO_CN:-5}; NGAN_SACH_GB=${NGAN_SACH_GB:-24}; NGAN_SACH_GOC_GB=${NGAN_SACH_GOC_GB:-60}
DIA_MIN_GB=${DIA_MIN_GB:-11}
NHAN_N=${NHAN_N:-24-27}
mkdir -p $H/khoa $H/log $H/hong $H/ss $H/xong $V/thu $V/kq
log() { echo "$(date '+%m-%d %T') $*" >> $H/tien_do.txt; }
ENV1="OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONNOUSERSITE=1"

bd_args() {
  local TV=T SO="--so ti_the=600 --so lysosome=200"
  case $1 in D1) SO="--so ti_the=15 --so lysosome=5";; D4) SO="--so ti_the=300 --so lysosome=100";; D4A2) TV=A; SO="--so ti_the=300 --so lysosome=100";; esac
  echo "--dac_ta dac_ta/v12/dong_vat_v12.json --seed $2 --N1 128 --K 4 --thu_vien $TV $SO --tang2 --khong_tu_do --eval_dung_tu_dau --out $3"
}
ds_ai() { [ $1 = D4A2 ] && echo "DS_AI=ai/khuon_chon/danh_sach_v2.json"; }
rss_hang_gb() {
  local s=0 p r
  for p in $(pgrep -f "[c]au_noi_cp_e3.py .*--thu_muc $V|[c]hay_cp_gpu.py --vao $V|[k]iem_cau_noi_e3.py --thu_muc $V"); do
    r=$(awk '/VmRSS/{print $2}' /proc/$p/status 2>/dev/null); [ -n "$r" ] && s=$((s + r)); done
  echo $((s / 1048576))
}
cho_ram() {
  [ "$DRY" = 1 ] && return 0
  while true; do
    local da=$(rss_hang_gb) trong=$(awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo)
    if [ $2 = goc ]; then local can_=$(($1 + 30)); [ $can_ -lt 80 ] && can_=80
      [ $((da + $1)) -le $((NGAN_SACH_GB + NGAN_SACH_GOC_GB)) ] && [ $trong -ge $can_ ] && return 0
    else [ $((da + $1)) -le $NGAN_SACH_GB ] && [ $trong -ge $(($1 + 30)) ] && return 0; fi
    sleep ${NGU:-30}
  done
}
gpu_v() {
  local T=$P/bo_dung_chung/tang_dg id x
  { [ -e $T/DUNG_U5 ] || [ -e $T/XONG_PT_U5 ]; } && { echo auto; return; }
  while read -r id x; do case $id in u5_*) [ -e $T/xong/$id ] || [ -e $T/loi/$id ] || { echo 1; return; };; esac; done < $T/viec.txt
  echo auto; }
HN=$H; DUNG_F=$V/DUNG; source $C/v122e3/gpu_n1.sh
gb_v() { cat $H/GB_V 2>/dev/null || echo 2; }
gb_n() { [ -s $H/GB_N_$1 ] && { cat $H/GB_N_$1; return; }; case $1 in D5) echo 20;; D4|D4A2) echo 12;; *) echo 2;; esac; }
m_phay_san() {
  [ "$DRY" = 1 ] && return 0
  local o=$1 s=$2 nh=$3 k=0 f fc mf=$V/cham/M_phay/MANIFEST.txt w=$H/cho_M_phay_$1_$2_$3
  [ $nh = chinh ] && k=$(grep -o '"k_chon": *[0-9]*' $V/thu/${o}_${s}_chon.json 2>/dev/null | tr -dc 0-9)
  [ -n "$k" ] || return 1
  f=$V/cham/M_phay/${o}_${s}_k${k}_${nh}.json; fc=$V/thu/${o}_${s}_k${k}/cp_$nh/tu_the_cp.json
  if [ -s $f ] && [ -s $fc ] && grep -q " $(md5sum $f | cut -c1-32) $(md5sum $fc | cut -c1-32) ${o}_${s}_k${k}_${nh}$" $mf 2>/dev/null; then rm -f $w $w.canh; return 0; fi
  [ -e $w ] || touch $w
  local t_; t_=$(stat -c %Y $w 2>/dev/null) || return 1
  if [ ! -e $w.canh ] && [ $(( $(date +%s) - t_ )) -gt 7200 ]; then touch $w.canh; log "!!!!! CANH_BAO: N_${o}_${s}_${nh} CHỜ M′ ĐÓNG BĂNG > 2 GIỜ (${o}_${s}_k${k}_${nh}) — kiểm m_phay_tat_ca.sh; KHÔNG nhập khi thiếu M′ !!!!!"; fi
  return 1
}
dia_du() {
  [ "$DRY" = 1 ] && [ -z "${DRY_DIA_MB:-}" ] && return 0
  local tr=${DRY_DIA_MB:-$(df -BM --output=avail $P | tail -1 | tr -dc 0-9)} can=$((DIA_MIN_GB * 1024)) k_
  for k_ in $(ls $H/khoa 2>/dev/null); do case $k_ in V_*_h1) can=$((can + 1000));; V_*|G_*) can=$((can + 400));; esac; done
  case $1 in V|G) can=$((can + 400));; esac
  [ $tr -ge $can ] && { rm -f $H/CHO_DIA; return 0; }
  [ -e $H/CHO_DIA ] || { touch $H/CHO_DIA; log "CHO_DIA: ổ trống ${tr} MB < ${can} MB => tạm ngừng nhận việc (không giết việc đang chạy)"; }
  return 1
}
chay() {
  local lg=$1; shift
  if [ "$DRY" = 1 ]; then echo "DRY $*" > $lg; sleep ${DRY_GIAY:-1}
    if [ -n "${DRY_KHONG_KIN:-}" ] && [[ "$*" == *xuat* ]] && [[ "$*" =~ $DRY_KHONG_KIN ]]; then echo "AssertionError: vung_cam: mesh ngăn không kín" >> $lg; return 1; fi
    [ -n "${DRY_HONG:-}" ] && [[ "$*" =~ $DRY_HONG ]] && { echo "DRY_HONG" >> $lg; return 1; }
    return 0; fi
  /usr/bin/time -f "THOI_GIAN %e s %M KB" "$@" > $lg 2>&1
}

sinh_viec() {
  [ -s $H/viec.txt ] && return
  {
    echo "X0 D1 6201 0"
    for o in D5 D4A2 D4 D1; do for s in $(seq 6201 6210); do for k in 0 1 2 3; do [ $o$s$k = D162010 ] || echo "X $o $s $k"; done; done; done
    echo "G D1 6201 0 co_dinh"
    echo "V D1 6201 0 chinh"
    for o in D5 D4A2 D4 D1; do for s in $(seq 6201 6210); do for k in 0 1 2 3; do [ $o$s$k = D162010 ] || echo "V $o $s $k chinh"; done; done; done
    echo "T D1"
    echo "P chon_cap"
    for o in D5 D4A2 D4 D1; do for s in $(seq 6201 6210); do echo "C $o $s"; done; done
    for o in D5 D4A2 D4 D1; do for s in $(seq 6201 6210); do echo "N $o $s chinh"; done; done
  } > $H/viec.txt
}
id_cua() { echo "$*" | tr ' ,' '__'; }
xuat_xong() { if [ $1 = D1 ] && [ $2 = 6201 ]; then xong X0_D1_6201_0; else xong X_$1_$2_$3; fi; }
seed_cap() {
  local s; for s in ${2//,/ }; do
    xuat_xong $1 $s $3 && { echo $s; return; }
    { [ -e $(dir_cf $1 $s $3)/XUAT_HONG ] || hong2 X_${1}_${s}_${3}; } && continue
    return; done; echo HONG; }
dir_cf() { echo $V/thu/$1_$2_k$3; }
xong() { [ -e $H/xong/$1 ]; }
hong2() { [ $(cat $H/hong/$1 2>/dev/null || echo 0) -ge 2 ]; }
kiem_dat() { [ "$DRY" = 1 ] && return 0; grep -q "^KIEM_DAT" $V/KIEM_CAU_NOI_D1.txt 2>/dev/null; }
v_bo() {
  [ -e $(dir_cf $1 $2 $3)/XUAT_HONG ] || hong2 X_$1_$2_$3 || hong2 V_$1_$2_$3_$4; }
dung_nhap() { [ -e $H/DUNG_NHAP ] || { touch $H/DUNG_NHAP; log "DUNG_NHAP: $* => không nhập gì, PI quyết"; }; }

san_sang() {
  set -- $@
  [ -e $XONG ] || return 1
  [ $1 = X0 ] && return 0
  kiem_dat || { [ -e $H/CHO_KIEM ] || { touch $H/CHO_KIEM; log "CHO_KIEM: pha 0 chưa KIEM_DAT => mọi việc khác chờ (xem KIEM_CAU_NOI_D1.txt)"; }; return 1; }
  case $1 in
    X) [ -e $V/DUYET_CP ] || return 1; return 0;;
    G) [ -e $V/DUYET_CP ] && [ -e $V/CHO_PHEP_RAM_GOC ] || return 1
       [ -n "$(seed_cap $2 $3 $4)" ] || return 1
       [ "$(seed_cap $2 $3 $4)" = HONG ] && return 2
       local g_ n_=0; for g_ in $(ls $H/khoa | grep '^G_'); do n_=$((n_ + 1)); set -- ${g_//_/ }
         [ "$DRY" = 1 ] || grep -aq LUOI_XONG $(dir_cf $2 $3 $4)/cp_goc/cp.log 2>/dev/null || return 1; done
       [ $n_ -lt ${SO_GOC:-3} ] || return 1; return 0;;
    T) [ -e $V/DUYET_CP ] && [ -e $V/MO_GPU_CP ] || return 1; gpu_san $(gb_n D1) || return 1; [ -e $H/KHAC ] && return 2
       xong V_D1_6201_0_chinh || { hong2 V_D1_6201_0_chinh && return 2; return 1; }; return 0;;
    V) [ -e $V/DUYET_CP ] && [ -e $V/MO_GPU_CP ] || return 1; gpu_san $(gb_v) || return 1; [ -e $H/KHAC ] && return 2
       { [ $2 = D1 ] && [ $3 = 6201 ] && xong X0_D1_6201_0; } || xong X_$2_$3_$4 || { [ -e $(dir_cf $2 $3 $4)/XUAT_HONG ] && return 2; hong2 X_$2_$3_$4 && return 2; return 1; }
       return 0;;
    P) [ -e $V/DUYET_CP ] || return 1
       local o_ s_ k_; for o_ in D5 D4A2; do for s_ in $(seq 6201 6210); do for k_ in 0 1 2 3; do
         xong V_${o_}_${s_}_${k_}_chinh || [ -e $(dir_cf $o_ $s_ $k_)/XUAT_HONG ] || hong2 V_${o_}_${s_}_${k_}_chinh || hong2 X_${o_}_${s_}_${k_} || return 1; done; done; done; return 0;;
    C) [ -e $H/KHAC ] && return 2; [ -e $V/DUYET_CP ] || return 1
       for k in 0 1 2 3; do local d=$(dir_cf $2 $3 $k); xong V_$2_$3_${k}_chinh || [ -e $d/XUAT_HONG ] || hong2 V_$2_$3_${k}_chinh || hong2 X_$2_$3_$k || return 1; done; return 0;;
    N) [ -e $V/MO_GPU_CP ] || return 1; gpu_san $(gb_n $2) || return 1; [ -e $H/KHAC ] && return 2; [ -e $H/DUNG_NHAP ] && return 2; [ -e $H/TUONG_DUONG_DAT ] || return 1
       xong T_D1 || { hong2 T_D1 && dung_nhap "lượt khói T_D1 hỏng hai lần"; return 1; }
       if [ $4 = chinh ]; then xong C_$2_$3 || { hong2 C_$2_$3 && { [ -e $V/thu/$2_$3_C_HONG ] || { touch $V/thu/$2_$3_C_HONG; log "C_HONG $2 $3: chọn k hỏng hai lần (cặp tính PLACE THUA, mục 3.7)"; }; return 2; }; return 1; }
         [ -e $V/thu/$2_$3_CHON_HONG ] && return 2
       else xong V_$2_$3_0_$4 || { v_bo $2 $3 0 $4 && return 2; return 1; }; fi
       m_phay_san $2 $3 $4 || return 1
       [ -n "$(ls $H/khoa | grep "^N_$2_$3_" 2>/dev/null)" ] && return 1; return 0;;
  esac
  return 1
}

viec_X() {
  local o=$1 s=$2 k=$3 cpu=$4 D=$(dir_cf $1 $2 $3) h=0.25
  [ -e $D/THU_H0125 ] && h=0.125
  rm -rf $D/dau_vao_cellpack.json $D/*.obj; mkdir -p $D
  cho_ram 3 thuong
  ( cd $P; chay $D/xuat_h$h.log env $ENV1 CUDA_VISIBLE_DEVICES="" $(ds_ai $o) taskset -c $cpu nice -n 5 $B $V/cau_noi_cp_e3.py xuat --thu_muc $D --k $k --thiet_bi cpu --h_luoi_ngan $h -- $(bd_args $o $s $D/bd) )
  local rc=$?
  if [ $rc = 0 ] && { [ "$DRY" = 1 ] || grep -aq XUAT_XONG $D/xuat_h$h.log; }; then rm -rf $D/bd; return 0; fi
  if grep -aq "không kín" $D/xuat_h$h.log 2>/dev/null; then
    if [ $h = 0.25 ]; then touch $D/THU_H0125 $H/khong_dem_$(id_cua X $o $s $k); log "XUAT_KHONG_KIN $o $s k$k (h 0.25) => thử lại h 0.125 (không tính là lần hỏng)"; return 1
    else grep -ao "AssertionError: [^\"]*" $D/xuat_h$h.log | head -1 > $D/XUAT_HONG; echo 1 > $H/hong/$(id_cua X $o $s $k); log "XUAT_HONG $o $s k$k: $(cat $D/XUAT_HONG) (không thử thêm)"; return 1; fi
  fi
  return 1
}
viec_G() {
  local o=$1 s=$(seed_cap $1 $2 $3) k=$3 cpu=$4 R; local D=$(dir_cf $o $s $k); echo "$o $s $k" > $H/G_${o}_k${k}.cap
  R=$D/cp_goc; rm -rf $R; mkdir -p $R
  cho_ram 40 goc
  chay $R/cp.log env $ENV1 CUDA_VISIBLE_DEVICES="" taskset -c $cpu $E $G/chay_cp_gpu.py --vao $D --ra_dir $R --seed $((10 * s + k)) --thu_tu ban_kinh_giam --bo tac_gia || { find $R/cp_out -name '*_grid.dat' -delete 2>/dev/null; return 1; }
  rm -rf $R/cp_out $R/sph; [ "$DRY" = 1 ] || grep -aq CP_XONG $R/cp.log
}
viec_V() {
  local o=$1 s=$2 k=$3 nh=$4 cpu=$5 D=$(dir_cf $1 $2 $3) R TT=ban_kinh_giam BO=tac_gia SP=""
  case $nh in thutu_md) TT=mac_dinh;; h1) SP="--spacing 0.1796875";; ja5) BO=mac_dinh;; esac
  R=$D/cp_$nh; rm -rf $R; mkdir -p $R
  cho_ram 4 thuong
  gpu_chay $(gb_v) $R/cp.log env $ENV1 PYTHONPATH=$G/ban/v4 CELLPACK_TORCH_SITE=$TS taskset -c $cpu \
      $E $G/chay_cp_gpu.py --vao $D --ra_dir $R --seed $((10 * s + k)) --thu_tu $TT --bo $BO $SP || { find $R/cp_out -name '*_grid.dat' -delete 2>/dev/null; return 1; }
  rm -rf $R/cp_out $R/sph; [ "$DRY" = 1 ] || grep -aq CP_XONG $R/cp.log
}
viec_P() {
  $E $C/v122/rs_dieu_kien_m.py $V/thu/*_k*/dau_vao_cellpack.json > $H/dieu_kien_m.txt 2>&1
  grep -q "CO THE ROI" $H/dieu_kien_m.txt && log "DIEU_KIEN_M: có đầu vào có thể rơi về đường dự phòng (xem dieu_kien_m.txt; bộ đếm sẽ cho biết)"
  if [ "$DRY" = 1 ]; then printf "G D5 6203 2 i_nhieu_bo_phan\nG D5 6207 1 ii_nhieu_loai\nG D5 6205 0 iii_thuong\nG D4A2 6206 0 i_nhieu_bo_phan+ii_nhieu_loai\nG D4A2 6205 0 iii_thuong\n" > $H/cap.txt
  else $E $V/chon_cap_e3.py $V/thu > $H/cap.txt 2> $H/cap_bang.txt || return 1; fi
  grep -q "^DUNG_RA_LAI" $H/cap.txt && { touch $H/DUNG_RA_LAI; log "DUNG_RA_LAI ghi nhan (quyet dinh 1 (B), PI 20/11): KHONG dung nhap; nhap chi mo khi moi cap (gom cap bo dem khac 0) TRUNG va khong KHAC; luot bo dem khac 0 mang nhan khong chung nhan tung bit"; }
  grep "^G " $H/cap.txt >> $H/viec.txt; log "CAP $(grep -c '^G ' $H/cap.txt) cặp: $(grep '^G ' $H/cap.txt | tr '\n' ';')"
}
kiem_cap() {
  local tong=0 dat=0 L
  while read -r L; do set -- $L; [ $1 = G ] || continue
    local sc=$(seed_cap $2 $3 $4); [ "$sc" = HONG ] && { [ -e $H/G_$2_k$4.khong_cap ] || { touch $H/G_$2_k$4.khong_cap; log "CAP_KHONG_LAP_DUOC $2 k$4: mọi seed ứng viên hỏng xuất (ghi vào báo cáo; cổng xét các cặp còn lại)"; }; continue; }
    tong=$((tong + 1)); [ -z "$sc" ] && continue
    local D=$(dir_cf $2 $sc $4) ten=$2_${sc}_k$4
    if hong2 $(id_cua $L) || v_bo $2 $sc $4 chinh; then dung_nhap "cặp $2 ${sc} k$4 không chạy xong (gốc hoặc v4 hỏng hai lần)"; continue; fi
    xong $(id_cua $L) && xong V_$2_${sc}_$4_chinh || continue
    if [ "$DRY" != 1 ] && { ! grep -aq "CELLPACK_FILE .*cellpack_ngoai/src/" $D/cp_goc/cp.log || ! grep -aq "CELLPACK_FILE .*cellpack_gpu/ban/v4/" $D/cp_chinh/cp.log; }; then
      [ -e $H/KHAC ] || { echo "$ten CELLPACK_FILE" > $H/KHAC; log "KHAC $ten: bản cellPACK nạp sai (CELLPACK_FILE)"; }; continue; fi
    [ -e $H/ss/$ten.json ] || { if [ "$DRY" = 1 ]; then echo '{"ket_luan": "TRUNG"}' > $H/ss/$ten.json; else $E $G/so_sanh.py $D/cp_goc $D/cp_chinh --ra $H/ss/$ten.json > /dev/null 2>&1; fi; log "SO_SANH $ten $(grep -o '"ket_luan": "[A-Z]*"' $H/ss/$ten.json)"; }
    if grep -q '"ket_luan": "TRUNG"' $H/ss/$ten.json; then dat=$((dat + 1)); else [ -e $H/KHAC ] || { echo $ten > $H/KHAC; log "KHAC $ten => DỪNG mọi việc nhập (N), PI quyết"; }; fi
  done < $H/viec.txt
  hong2 P_chon_cap && dung_nhap "chọn cặp hỏng hai lần"
  xong P_chon_cap && [ $tong -gt 0 ] && [ $dat = $tong ] && [ ! -e $H/TUONG_DUONG_DAT ] && [ ! -e $H/KHAC ] && { touch $H/TUONG_DUONG_DAT; log "TUONG_DUONG_DAT ($dat/$tong cặp TRUNG) => mở nhập"; }
}
viec_C() {
  local o=$1 s=$2 cpu=$3 k D tong=0 ok=0
  cho_ram 4 thuong
  for k in 0 1 2 3; do D=$(dir_cf $o $s $k); xong V_${o}_${s}_${k}_chinh || continue
    ( cd $P; chay $D/chon.log env $ENV1 CUDA_VISIBLE_DEVICES="" $(ds_ai $o) taskset -c $cpu $B $V/cau_noi_cp_e3.py chon --thu_muc $D --k $k --ket_qua_cp $D/cp_chinh/tu_the_cp.json --thiet_bi cpu -- $(bd_args $o $s $D/bd) ) || return 1
    [ "$DRY" = 1 ] && echo "{\"k\": $k, \"so_khong_dat\": $((RANDOM % 9)), \"phat_cung_day_du\": 0.0$RANDOM, \"giay_dat_cp\": 1, \"giay_khoi_L\": 1, \"giay_cham_chon\": 1}" > $D/chon_k$k.json
    ok=$((ok + 1)); done
  if [ $ok = 0 ]; then touch $V/thu/${o}_${s}_CHON_HONG; log "CHON_HONG $o $s: không cấu hình nào xuất/đặt được (mục 3.7: cặp tính PLACE THUA)"; return 0; fi
  $E - $V/thu $o $s > $V/thu/${o}_${s}_chon.json <<'PY' || return 1
import json, sys, os
thu, o, s = sys.argv[1:]; ds = []
for k in range(4):
    f = f"{thu}/{o}_{s}_k{k}/chon_k{k}.json"
    if os.path.exists(f): ds.append(json.load(open(f)))
best = min(ds, key=lambda d: (d["so_khong_dat"], d["phat_cung_day_du"]))
tong = sum(d["giay_dat_cp"] + (d["giay_khoi_L"] or 0) + d["giay_cham_chon"] for d in ds)
print(json.dumps({"k_chon": best["k"], "K_hieu_luc": len(ds), "cac_k": [d["k"] for d in ds], "giay_dat_tong": tong, "chi_tiet": ds}))
PY
  log "CHON $o $s => $(cat $V/thu/${o}_${s}_chon.json | cut -c1-80)"
}
viec_T() {
  local cpu=$1 D=$(dir_cf D1 6201 0) O=$V/kq/D1_cpthu_6201
  ( cd $P; chay $D/chon.log env $ENV1 CUDA_VISIBLE_DEVICES="" taskset -c $cpu $B $V/cau_noi_cp_e3.py chon --thu_muc $D --k 0 --ket_qua_cp $D/cp_chinh/tu_the_cp.json --thiet_bi cpu -- $(bd_args D1 6201 $D/bd) ) || return 1
  local gd=0; [ "$DRY" = 1 ] || gd=$($E -c "import json;d=json.load(open('$D/chon_k0.json'));print(d['giay_dat_cp']+(d['giay_khoi_L'] or 0)+d['giay_cham_chon'])") || return 1
  rm -rf $D/bd $O
  ( cd $P; gpu_chay $(gb_n D1) $D/nhap_thu.log env OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 PYTHONNOUSERSITE=1 taskset -c $NHAN_N $B $V/cau_noi_cp_e3.py nhap --thu_muc $D --k 0 \
       --ket_qua_cp $D/cp_chinh/tu_the_cp.json --thiet_bi cuda --giay_dat $gd -- $(bd_args D1 6201 $D/bd) ) || return 1
  [ "$DRY" = 1 ] || grep -aq NHAP_XONG $D/nhap_thu.log || return 1
  [ "$DRY" = 1 ] && mkdir -p $D/bd; mv $D/bd $O && log "T_D1: chon + nhap D1 6201 k0 XONG ($O)"
}
viec_N() {
  local o=$1 s=$2 nh=$3 k gd D GB=17 O PY=$E
  GB=$(gb_n $o)
  if [ $nh = chinh ]; then k=$($PY -c "import json;print(json.load(open('$V/thu/${o}_${s}_chon.json'))['k_chon'])" 2>/dev/null || echo 0)
       gd=$($PY -c "import json;print(json.load(open('$V/thu/${o}_${s}_chon.json'))['giay_dat_tong'])" 2>/dev/null || echo 0)
  else k=0; gd=$($PY -c "import json;d=json.load(open('$(dir_cf $o $s 0)/cp_$nh/tu_the_cp.json'));print(d['giay_dat'])" 2>/dev/null || echo 0); fi
  D=$(dir_cf $o $s $k); O=$V/kq/${o}_cp${nh}_${s}
  rm -rf $D/bd $O
  ( cd $P; gpu_chay $GB $D/nhap_$nh.log env OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 PYTHONNOUSERSITE=1 PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True $(ds_ai $o) taskset -c $NHAN_N $B $V/cau_noi_cp_e3.py nhap --thu_muc $D --k $k \
       --ket_qua_cp $D/cp_$nh/tu_the_cp.json --thiet_bi cuda --giay_dat $gd -- $(bd_args $o $s $D/bd) ) || return 1
  [ "$DRY" = 1 ] || grep -aq NHAP_XONG $D/nhap_$nh.log || return 1
  [ "$DRY" = 1 ] && mkdir -p $D/bd; mv $D/bd $O && echo "$nh k$k" > $O/NHAN_CELLPACK
}

cong_nhan() {
  local w=$1 cpu=$((24 + $1))
  while true; do
    [ -e $V/DUNG ] && { log "W$w dừng (DUNG)"; exit 0; }
    local lam=0 con=0 L id
    while read -r L; do
      [ -z "$L" ] && continue; id=$(id_cua $L)
      xong $id && continue; hong2 $id && continue
      san_sang $L; local r=$?; [ $r = 2 ] && continue; [ $r = 1 ] && { con=1; continue; }
      dia_du ${L%% *} || { con=1; continue; }
      mkdir $H/khoa/$id 2>/dev/null || { local op=$(cat $H/khoa/$id/pid 2>/dev/null); if [ -n "$op" ] && ! kill -0 $op 2>/dev/null; then log "W$w nhận lại khoá $id"; else con=1; continue; fi; }
      echo $$ > $H/khoa/$id/pid; log "W$w BAT_DAU $id (nhân $cpu)"; local t0=$(date +%s) ok=1
      set -- $L; ID_HT=$id
      case $1 in
        X0) viec_X D1 6201 0 $cpu && { ( cd $P; chay $(dir_cf D1 6201 0)/kiem.log env $ENV1 CUDA_VISIBLE_DEVICES="" taskset -c $cpu $B $V/kiem_cau_noi_e3.py --thu_muc $(dir_cf D1 6201 0) \
              --tt_nn1x $KQ/e3b_B_D1_nn1x_6201/thong_tin.json --tt_td $KQ/e3b_B_D1_td_6201/thong_tin.json -- $(bd_args D1 6201 $(dir_cf D1 6201 0)/bd) ); \
              grep -aho "KIEM_[A-Z]*.*" $(dir_cf D1 6201 0)/kiem.log | head -1 > $V/KIEM_CAU_NOI_D1.txt; log "KIEM_CAU_NOI D1 6201 k0: $(cut -c1-200 $V/KIEM_CAU_NOI_D1.txt)"; true; } || ok=0;;
        X) viec_X $2 $3 $4 $cpu || ok=0;;
        G) viec_G $2 $3 $4 $cpu || ok=0;;
        V) viec_V $2 $3 $4 $5 $cpu || ok=0;;
        C) viec_C $2 $3 $cpu || ok=0;;
        P) viec_P || ok=0;;
        T) viec_T $cpu || ok=0;;
        N) viec_N $2 $3 $4 || ok=0;;
      esac
      if [ $ok = 1 ]; then touch $H/xong/$id; log "W$w XONG $id $(( $(date +%s) - t0 )) s"
      elif [ -e $H/khong_dem_$id ]; then log "W$w THU_LAI $id (khong tinh lan hong: $(head -c 120 $H/khong_dem_$id 2>/dev/null))"; rm -f $H/khong_dem_$id
      else local n=$(( $(cat $H/hong/$id 2>/dev/null || echo 0) + 1 )); echo $n > $H/hong/$id; log "W$w HONG $id (lần $n$( [ $n -ge 2 ] && echo ', bỏ' ))"; fi
      rm -rf $H/khoa/$id
      { [ $1 = G ] || [ $1 = V ] || [ $1 = P ]; } && kiem_cap
      lam=1; break
    done < $H/viec.txt
    if [ $lam = 0 ]; then
      if [ $con = 0 ]; then log "W$w HET_VIEC"; touch $H/W_${w}_HET; exit 0; fi
      sleep ${NGU:-60}
    fi
  done
}

giam_sat() {
  while [ ! -e $H/DUNG_GS ]; do
    if [ -e $XONG ] && [ "$DRY" != 1 ]; then
      for p in $(pgrep -f "[c]au_noi_cp_e3.py .*--thu_muc $V|[c]hay_cp_gpu.py --vao $V|[k]iem_cau_noi_e3.py --thu_muc $V"); do
        r=$(awk '/VmRSS/{print $2}' /proc/$p/status 2>/dev/null); [ -n "$r" ] && [ $r -gt $((50*1048576)) ] && { kill $p; log "CANH_RSS kill $p rss_kb=$r"; }
      done
      echo "$(date '+%m-%d %T') $(nvidia-smi --query-gpu=index,utilization.gpu,memory.used --format=csv,noheader | tr '\n' ' ') tai=$(cut -d' ' -f1 /proc/loadavg) RAM_trong=$(awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo)G rss_hang=$(rss_hang_gb)G khoa=$(ls $H/khoa | wc -l) xong=$(ls $H/xong | wc -l)/$(wc -l < $H/viec.txt)" >> $H/giam_sat.log
    fi
    [ $(ls $H/W_*_HET 2>/dev/null | wc -l) -ge $SO_CN ] && { log "GIAM_SAT: mọi công nhân hết việc"; touch $V/XONG_HANG; exit 0; }
    sleep ${NGU:-60}
  done
}

case ${1:-} in
  khoi_dong)
    if [ "$DRY" != 1 ]; then exec 9>$H/khoi_dong.khoa; flock -n 9 || { echo "đã có bản đang chạy"; exit 1; }; fi
    sinh_viec; echo $SO_CN > $H/so_cong_nhan; log "KHOI_DONG $SO_CN công nhân (chờ $XONG)"
    while [ ! -e $XONG ] && [ ! -e $V/DUNG ]; do sleep ${NGU:-60}; done
    [ -e $V/DUNG ] && { log "DUNG trước khi bắt đầu"; exit 0; }
    log "THAY XONG_DUNG_CHINH => khởi công nhân"
    if [ "$DRY" != 1 ]; then
      taskset -c 24 $E -c "import numpy,scipy,trimesh,rtree,sys;print('numpy',numpy.__version__,'scipy',scipy.__version__,'trimesh',trimesh.__version__,'rtree',rtree.__version__,'py',sys.version.split()[0])" > $H/moi_truong.txt 2>&1
      ( cd $G && md5sum -c ban/v4.md5 ) >> $H/moi_truong.txt 2>&1 || { touch $H/KHAC; log "MD5_V4_SAI => KHAC (không chạy v4/nhập)"; }
      md5sum $V/cau_noi_cp_e3.py $V/kiem_cau_noi_e3.py $V/chon_cap_e3.py $C/v122/rs_dieu_kien_m.py $0 $G/chay_cp_gpu.py $G/so_sanh.py >> $H/moi_truong.txt
      cat $C/src/.git/refs/heads/main >> $H/moi_truong.txt
      bn=$(taskset -c 24 git -C $C/src status --porcelain --untracked-files=no | wc -l); echo "src_thay_doi $bn" >> $H/moi_truong.txt
      [ "$bn" = 0 ] || { touch $H/KHAC; log "SRC_GOC_KHONG_SACH ($bn tệp) => KHAC"; }; fi
    if [ "$DRY" = 1 ]; then for w in $(seq 0 $((SO_CN - 1))); do bash $0 cong_nhan $w & done; bash $0 giam_sat & wait
    else for w in $(seq 0 $((SO_CN - 1))); do setsid nohup taskset -c $((24 + w)) bash $0 cong_nhan $w > $H/log/cn_$w.out 2>&1 < /dev/null & echo $! > $H/cn_$w.pid; done
         setsid nohup taskset -c 24-28 nice -n 19 bash $0 giam_sat > $H/log/gs.out 2>&1 < /dev/null & echo $! > $H/gs.pid; fi;;
  cong_nhan) cong_nhan $2;;
  giam_sat) giam_sat;;
  *) echo "dùng: $0 khoi_dong | cong_nhan <w> | giam_sat"; exit 1;;
esac
