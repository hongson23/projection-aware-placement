#!/bin/bash
R=~/Vision_team_data_189/sonhh20/ch3_tebao_dep/bo_kiem_v3; Q=$R/cong_v3; P=~/miniforge3/envs/b2demo/bin/python
mkdir -p $Q/chay $Q/xong $Q/loi $Q/log
cd $R
mk(){ d=$R/$1; mkdir -p $d; for f in cay_loi.py cay_loi_tv.py dung_tv_thu.py kiem2.py kiem3.py kiem_thucvat.py luat_thucvat_v14_van_ban.json luat_v12 luat_v12.py luat_v41.json so_sanh.py tang_toc.py thu_roi.py tv_thu v3_cay_loi_chung.py; do ln -sfn ../$f $d/$f; done; }
: > $Q/viec.txt
for i in 0 1 2 3; do mk hn3_D4A2_td_1001_p$i; echo "D4A2_p$i|3|cd $R/hn3_D4A2_td_1001_p$i && CHIA=$i/4 $P cay_loi.py D4A2_td_1001" >> $Q/viec.txt; done
for i in 0 1; do mk hn3_TV1_p$i; echo "TV1_p$i|3|cd $R/hn3_TV1_p$i && CHIA=$i/2 $P -c \"import cay_loi_tv as C; C.chay(ds=('TV1',))\"" >> $Q/viec.txt; done
for i in 0 1; do mk hn3_TV6_p$i; echo "TV6_p$i|3|cd $R/hn3_TV6_p$i && CHIA=$i/2 $P -c \"import cay_loi_tv as C; C.chay(ds=('TV6',))\"" >> $Q/viec.txt; done
for i in 0 1; do mk hn3_D1_td_1001_p$i; echo "D1_p$i|3|cd $R/hn3_D1_td_1001_p$i && CHIA=$i/2 $P cay_loi.py D1_td_1001" >> $Q/viec.txt; done
mk hn3_TV2; echo "TV2|3|cd $R/hn3_TV2 && $P -c \"import cay_loi_tv as C; C.chay(ds=('TV2',))\"" >> $Q/viec.txt

cat > $Q/cong_nhan.sh <<"EOF"
#!/bin/bash
Q=~/Vision_team_data_189/sonhh20/ch3_tebao_dep/bo_kiem_v3/cong_v3; CN=$1
export OMP_NUM_THREADS=3 MKL_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 NUMBA_NUM_THREADS=3 IGL_NUM_THREADS=3 PYTHONNOUSERSITE=1
while true; do
  [ -f $Q/DUNG ] && break
  lay=""
  while IFS='|' read -r ten gb lenh; do
    [ -z "$ten" ] && continue
    if mkdir $Q/chay/$ten 2>/dev/null; then lay="$ten|$gb|$lenh"; break; fi
  done < $Q/viec.txt
  [ -z "$lay" ] && break
  IFS='|' read -r ten gb lenh <<< "$lay"
  echo "$(date +%F_%T) $CN BAT_DAU $ten" >> $Q/nhat_ky.txt; t0=$(date +%s)
  CHUNG=1 ~/Vision_team_data_189/sonhh20/chay_gpu.sh auto $gb /usr/bin/time -v taskset -c 12-31 bash -c "$lenh" > $Q/log/$ten.log 2>&1
  rc=$?
  echo "$(date +%F_%T) $CN XONG $ten rc=$rc $(( $(date +%s) - t0 ))s" >> $Q/nhat_ky.txt
  [ $rc = 0 ] && touch $Q/xong/$ten || touch $Q/loi/$ten
done
echo "$(date +%F_%T) $CN THOAT" >> $Q/nhat_ky.txt
EOF
cat > $Q/giam_sat.sh <<"EOF"
#!/bin/bash
Q=~/Vision_team_data_189/sonhh20/ch3_tebao_dep/bo_kiem_v3/cong_v3; R=$Q/..
while [ ! -e $R/XONG_CONG_V3 ]; do
  g=$(nvidia-smi --query-gpu=index,utilization.gpu,memory.used --format=csv,noheader | paste -sd' ')
  echo "$(date +%F_%T) $g | load $(cut -d' ' -f1 /proc/loadavg) | xong $(ls $Q/xong | wc -l) loi $(ls $Q/loi | wc -l) chay $(ls $Q/chay | wc -l)/$(grep -c . $Q/viec.txt)" >> $Q/giam_sat.txt
  sleep 60
done
EOF
cat > $Q/buoc_cuoi.sh <<"EOF"
#!/bin/bash
Q=~/Vision_team_data_189/sonhh20/ch3_tebao_dep/bo_kiem_v3/cong_v3; R=$Q/..; P=~/miniforge3/envs/b2demo/bin/python
N=$(grep -c . $Q/viec.txt)
while [ $(( $(ls $Q/xong | wc -l) + $(ls $Q/loi | wc -l) )) -lt $N ]; do sleep 60; done
cd $R && taskset -c 12-31 $P tong_hop_cong_v3.py > $Q/tong_hop.log 2>&1
touch $R/XONG_CONG_V3
EOF
chmod +x $Q/*.sh
for k in 1 2 3 4 5 6; do setsid nohup $Q/cong_nhan.sh cn$k > /dev/null 2>&1 & echo $! >> $Q/pid.txt; sleep 2; done
setsid nohup $Q/giam_sat.sh > /dev/null 2>&1 & echo $! >> $Q/pid.txt
setsid nohup $Q/buoc_cuoi.sh > /dev/null 2>&1 & echo $! >> $Q/pid.txt
cat $Q/viec.txt | cut -d'|' -f1 | paste -sd' '; cat $Q/pid.txt | paste -sd' '
