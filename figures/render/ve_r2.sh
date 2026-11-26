#!/bin/bash
set -e
cd ~/Vision_team_data_189/sonhh20/ch3_tebao_dep/bo_dung_chung/hinh_cad
K=../tang_dg/kq; VMAX=${VMAX:-1.2479}; O=$PWD/r2; mkdir -p $O
B="taskset -c 30,31 $HOME/blender-4.2.1-linux-x64/blender -b -t 2 --python-exit-code 1 -P scripts/ve_hd_blender.py"
FIX='{"vach": "#C8C8C8", "mang": "#B0C4DE", "khong_bao": "#9FC5E8", "luoi_noi_chat": "#8FA3B5"}'
for a in ${NHANH:-u5ktt td nn1x}; do
  J=$K/u5_tv_c_real_${a}_7311
  env MODE=cat MESH=$J/mesh_vo_bao.npz HD=du_lieu/hd_u5_tv_c_real_${a}_7311.json LUT=du_lieu/lut_magma_r.json VMAX=$VMAX \
      FIXED_COL="$FIX" LOAI_LON=mang,vach,khong_bao,luoi_noi_chat,nhan LOAI_NHO=luc_lap,ti_the,peroxisome,golgi \
      ALPHA_mang=0.08 ALPHA_vach=0.08 ALPHA_khong_bao=0.15 BO_LOAI=luoi_noi_chat \
      OUT=$O/cat_$a.png BLEND=${BLEND_ON:+$O/cat_$a.blend} RES=3000 RESY=1500 SPP=256 ELEV=15 ORTHO=84 $B
done
touch $O/XONG_VE
