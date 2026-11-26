#!/bin/bash
set -e
cd ~/Vision_team_data_189/sonhh20/ch3_tebao_dep/bo_dung_chung/hinh_cad
K=../tang_dg/kq; VMAX=0.2062; O=$PWD/r3; mkdir -p $O
B="taskset -c 30,31 $HOME/blender-4.2.1-linux-x64/blender -b -t 2 --python-exit-code 1 -P scripts/ve_hd_blender.py"
FIX='{"mang": "#D9D9D9", "nhan": "#B8B8D0", "luoi_noi_chat": "#8FA3B5"}'
for a in U5b td; do
  J=$K/u5_n14d5_d5_707681_${a}_6543
  for v in 3d cat; do
    [ $v = 3d ] && C="Y0=100 ELEV=35" || C="Y0=0 ELEV=55"
    env MODE=cat CLIP_AX=z CLIP_SIGN=gt $C MESH=$J/mesh_vo_bao.npz HD=du_lieu/hd_u5_n14d5_d5_707681_${a}_6543.json LUT=du_lieu/lut_magma_r.json VMAX=$VMAX \
        FIXED_COL="$FIX" ALPHA_mang=0.12 ALPHA_ER=0.30 OUT=$O/${v}_$a.png BLEND=$O/${v}_$a.blend RES=1600 RESY=1200 SPP=256 ORTHO=27 $B
  done
done
touch $O/XONG_VE
