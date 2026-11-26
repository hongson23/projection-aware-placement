#!/bin/bash
set -u
DS=$(readlink -f "$1"); VSB=$(readlink -m "$2"); N=${3:-4}; GB=${4:-6}
K=$(dirname "$(readlink -f "$0")"); P=~/miniforge3/envs/b2demo/bin/python
TT=${TIEN_TO:-vsb_}; NEN=${VS_NEN:-lzma}
mkdir -p "$VSB"; NK=$VSB/nhat_ky_sua_mesh.txt
echo "$(date '+%F %T') BAT_DAU ds=$DS vsb=$VSB tien_to=$TT" >> $NK
( cd $K && md5sum -c --quiet SUA_MESH_MD5.txt ) >> $NK 2>&1 || { echo "MD5 sua_mesh KHAC — dung" >> $NK; touch $VSB/LOI_SUA_MESH; exit 1; }
md5sum $K/sua_mesh_dn.py $K/cham_lo_vsb.sh >> $NK
if [ ! -e "$VSB/vs/XONG_VE_SINH" ]; then
  CHI_VE_SINH=1 TIEN_TO=$TT bash $K/cham_lo_vs.sh "$DS" "$VSB/vs" "$N" "$GB" || { echo "cham_lo_vs.sh loi (rc=$?)" >> $NK; touch $VSB/LOI_SUA_MESH; exit 1; }
fi
[ -e "$VSB/vs/XONG_VE_SINH" ] && [ ! -e "$VSB/vs/LOI_VE_SINH" ] || { echo "ve sinh chua xong / loi" >> $NK; touch $VSB/LOI_SUA_MESH; exit 1; }
DSV=$VSB/vs/$(basename "$DS")
D=$(dirname "$K")
( cd $D && awk '/^## A\./{a=1;next} /^## /{a=0} a && /^[0-9a-f]{32}  /' $K/DONG_BANG_BO_KIEM_V3.txt | md5sum -c --quiet ) >> $NK 2>&1 \
  || { echo "MD5 muc A KHAC dong bang — dung" >> $NK; touch $VSB/LOI_SUA_MESH; exit 1; }
[ "$(awk '/^## A\./{a=1;next} /^## /{a=0} a && /^[0-9a-f]{32}  /' $K/DONG_BANG_BO_KIEM_V3.txt | wc -l)" = 16 ] || { echo "muc A khong du 16 dong — dung" >> $NK; touch $VSB/LOI_SUA_MESH; exit 1; }
( cd $K && md5sum -c --quiet VE_SINH_LUOI_MD5.txt && md5sum -c --quiet SUA_MESH_MD5.txt ) >> $NK 2>&1 || { echo "MD5 ve_sinh/sua_mesh KHAC — dung" >> $NK; touch $VSB/LOI_SUA_MESH; exit 1; }
DS2=$VSB/sb/$(basename "$DS"); mkdir -p "$VSB/sb"; TAM=$DS2.tam; : > $TAM
while read -r d v l; do
  if [ -z "${d:-}" ] || [[ "$d" == \#* ]]; then echo "${d:-} ${v:-} ${l:-}" >> $TAM; continue; fi
  [ -f "$d/mesh_vo_bao.npz" ] || { echo "thieu $d/mesh_vo_bao.npz" >> $NK; touch $VSB/LOI_SUA_MESH; exit 1; }
  cha=$(basename "$(dirname "$d")"); id=$(basename "$d"); o=$VSB/sb/$cha/$id
  cr=$K/cham_ra/${cha}__$id
  if [ -e "$cr/mesh.npz" ] && [ "$(readlink -f "$cr/mesh.npz")" != "$(readlink -f "$o/mesh_vo_bao.npz")" ]; then
    echo "cham_ra/${cha}__$id da co va tro noi khac — dung" >> $NK; touch $VSB/LOI_SUA_MESH; exit 1; fi
  if [ ! -e "$o/XONG_SUA_MESH" ]; then
    rm -rf "$o"; mkdir -p "$o"
    ${VS_LOI:+taskset -c $VS_LOI} env OMP_NUM_THREADS=1 PYTHONNOUSERSITE=1 $P $K/sua_mesh_dn.py "$d/mesh_vo_bao.npz" "$o/mesh_vo_bao.npz" \
       --vat "${v:-dong_vat}" --nhat_ky "$o/sua_mesh.json" --nen $NEN --bang "$VSB/sb/trang_thai_sua.tsv" >> $NK 2>&1 || { echo "LOI sua_mesh $d (rc=$?)" >> $NK; touch $VSB/LOI_SUA_MESH; exit 1; }
    [ -f "$d/thong_tin.json" ] && cp -p "$d/thong_tin.json" "$o/thong_tin.json"
    touch "$o/XONG_SUA_MESH"
  fi
  echo "$o ${v:-dong_vat} ${l:-}" >> $TAM
done < "$DSV"
mv $TAM $DS2; echo "$(date '+%F %T') XONG_SUA_MESH $(grep -vc '^#' $DS2) dong -> $DS2" >> $NK; touch $VSB/XONG_SUA_MESH
[ "${CHI_SUA:-0}" = 1 ] && exit 0
exec bash $K/cham_lo.sh "$DS2" "$N" "$GB"
