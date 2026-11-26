#!/bin/bash
set -u
DS=$(readlink -f "$1"); VS=$(readlink -m "$2"); N=${3:-4}; GB=${4:-6}
K=$(dirname "$(readlink -f "$0")"); D=$(dirname "$K"); P=~/miniforge3/envs/b2demo/bin/python
TT=${TIEN_TO:-vs_}; NEN=${VS_NEN:-lzma}
mkdir -p "$VS"; NK=$VS/nhat_ky_ve_sinh.txt
echo "$(date '+%F %T') BAT_DAU ds=$DS vs=$VS tien_to=$TT nen=$NEN" >> $NK
( cd $D && awk '/^## A\./{a=1;next} /^## /{a=0} a && /^[0-9a-f]{32}  /' $K/DONG_BANG_BO_KIEM_V3.txt | md5sum -c --quiet ) >> $NK 2>&1 \
  || { echo "MD5 muc A KHAC dong bang — dung" >> $NK; touch $VS/LOI_VE_SINH; exit 1; }
[ "$(awk '/^## A\./{a=1;next} /^## /{a=0} a && /^[0-9a-f]{32}  /' $K/DONG_BANG_BO_KIEM_V3.txt | wc -l)" = 16 ] || { echo "muc A khong du 16 dong — dung" >> $NK; touch $VS/LOI_VE_SINH; exit 1; }
( cd $K && md5sum -c --quiet VE_SINH_LUOI_MD5.txt ) >> $NK 2>&1 || { echo "MD5 ve_sinh_luoi.py KHAC — dung" >> $NK; touch $VS/LOI_VE_SINH; exit 1; }
md5sum $K/ve_sinh_luoi.py $K/cham_lo_vs.sh $K/cham_lo.sh >> $NK
DS2=$VS/$(basename "$DS"); TAM=$DS2.tam; : > $TAM
while read -r d v l; do
  if [ -z "${d:-}" ] || [[ "$d" == \#* ]]; then echo "${d:-} ${v:-} ${l:-}" >> $TAM; continue; fi
  [ -f "$d/mesh_vo_bao.npz" ] || { echo "thieu $d/mesh_vo_bao.npz" >> $NK; touch $VS/LOI_VE_SINH; exit 1; }
  cha=$TT$(basename "$(dirname "$d")"); id=$(basename "$d"); o=$VS/$cha/$id
  cr=$K/cham_ra/${cha}__$id
  if [ -e "$cr/mesh.npz" ] && [ "$(readlink -f "$cr/mesh.npz")" != "$(readlink -f "$o/mesh_vo_bao.npz")" ]; then
    echo "cham_ra/${cha}__$id da co va tro noi khac — dung" >> $NK; touch $VS/LOI_VE_SINH; exit 1; fi
  if [ ! -e "$o/XONG_VE_SINH" ]; then
    rm -rf "$o"; mkdir -p "$o"
    ${VS_LOI:+taskset -c $VS_LOI} env OMP_NUM_THREADS=1 PYTHONNOUSERSITE=1 $P $K/ve_sinh_luoi.py "$d/mesh_vo_bao.npz" "$o/mesh_vo_bao.npz" \
       --nhat_ky "$o/ve_sinh.json" --nen $NEN >> $NK 2>&1 || { echo "LOI ve sinh $d (rc=$?)" >> $NK; touch $VS/LOI_VE_SINH; exit 1; }
    [ -f "$d/thong_tin.json" ] && cp -p "$d/thong_tin.json" "$o/thong_tin.json"
    touch "$o/XONG_VE_SINH"
  fi
  echo "$o ${v:-dong_vat} ${l:-}" >> $TAM
done < "$DS"
mv $TAM $DS2; echo "$(date '+%F %T') XONG_VE_SINH $(grep -vc '^#' $DS2) dong -> $DS2" >> $NK; touch $VS/XONG_VE_SINH
[ "${CHI_VE_SINH:-0}" = 1 ] && exit 0
exec bash $K/cham_lo.sh "$DS2" "$N" "$GB"
