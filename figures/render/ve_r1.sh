#!/bin/bash
set -e
cd ~/Vision_team_data_189/sonhh20/ch3_tebao_dep/bo_dung_chung/hinh_cad
K=../tang_dg/kq; S=6512; VMAX=0.2771; O=$PWD/r1; mkdir -p $O
B="taskset -c 30,31 $HOME/blender-4.2.1-linux-x64/blender -b -t 2 --python-exit-code 1 -P scripts/ve_hd_blender.py"
for a in P8 U5b; do
  J=$K/u5_r_${a}_k08_$S
  env MODE=cat MESH=$J/mesh_vo_bao.npz HD=du_lieu/hd_u5_r_${a}_k08_$S.json LUT=du_lieu/lut_magma_r.json VMAX=$VMAX \
      OUT=$O/cat_$a.png BLEND=$O/cat_$a.blend RES=1600 SPP=256 ELEV=15 ORTHO=22 $B
  i=0
  for zc in "0.57,5.83,6.67" "-4.84,5.72,4.59" "5.56,2.75,-6.47"; do
    i=$((i+1))
    env MODE=zoom ZC=$zc ZR=1.2 MESH=$J/mesh_vo_bao.npz TRUOC=$J/mesh_vo_bao_truoc.npz HD=du_lieu/hd_u5_r_${a}_k08_$S.json \
        LUT=du_lieu/lut_magma_r.json VMAX=$VMAX OUT=$O/zoom${i}_$a.png RES=900 SPP=256 ELEV=15 $B
  done
done
touch $O/XONG_VE
