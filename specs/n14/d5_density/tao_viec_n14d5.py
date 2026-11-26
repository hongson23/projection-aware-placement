# -*- coding: utf-8 -*-
"""tao_viec_n14d5.py (22/11/2026) — danh sách việc N14 MẬT ĐỘ CAO (W3, Bài 1 ch3). CHỈ SINH TỆP n14d5/viec_<giai_đoạn>.txt; KHÔNG chèn hàng
(chèn = chen_n14d5.py, chỉ sau khi PI duyệt). Bộ dựng: bo_dung_chung/ma3/botri_chung_v122u5n14.py (sinh bởi tao_ma_u5n14.py, md5 trong MD5_N14D5.txt).
Số bào quan: n14d5/so_luong_n14d5.json (đóng băng từ kha_thi.json, luật C: cùng φ_eff TỔNG như ô tổng hợp D5 / k08, ti : ly = 3 : 1).
Nhánh (cờ = ô D5 lượt chính 1×, như bảng thống nhất N1; mọi nhánh cùng bộ dựng, cùng đặc tả, cùng số, cùng seed):
  U5b  = PLACE K 4 + --u5 --u5_lambda 1.0 --u5_kappa 0 --u5_chon gc --u5_er --u5_tam   (u5_c_U5b_D5_*)
  td   = PLACE K 4                                                               (e3b_B_D5_td_*)
  nn1x = thử–loại --so_lo 3 --K 6                                                (e3b_B_D5_nn1x_*)
  unc8 = SA cấu hình c08 (b2/SA_D5.txt, token @SA_D5@ do chay_tdg.sh thay) 1102/3222 vòng (e1_D5_unt_*)
Chung: --tang2 --khong_tu_do --eval_dung_tu_dau --hd_truoc_sau; E3_ATLAS_XA = atlas N14; thu_vien T.
Giai đoạn: id (4 việc đồng nhất cờ-tắt ≡ v122n14: lệnh CŨ của n14_798836_{td,nn1x,un1x}_6501, + td có --hd_truoc_sau) | thu (2 tế bào × 2 biến thể × 4 nhánh, seed 6549) |
           d5_1, k8_1 (seed 6541; mốc 1/3) | d5_23, k8_23 (seed 6542–6543).
Dùng: python tao_viec_n14d5.py   (đọc so_luong_n14d5.json và gb_n14d5.json nếu có — gb do kiem_n14d5.py thu ghi sau bước thử)."""
import os, json
R = os.path.expanduser("~/Vision_team_data_189/sonhh20/ch3_tebao_dep"); B = "bo_dung_chung"; D = f"{R}/{B}/n14d5"; KQ = f"{B}/tang_dg/kq"
PY = "/home/thainq97/miniforge3/envs/b2demo/bin/python"; MA = f"{B}/ma3/botri_chung_v122u5n14.py"
ENV = f"GO_LOI=1 E3_ATLAS_XA={B}/n14/atlas_xa_cache env PYTHONNOUSERSITE=1 PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True"
CELL = ["798836", "989830", "593409", "707681", "997092"]; THU = ["798836", "997092"]
NH = {"U5b": "--cach tu_dong --u5 --u5_lambda 1.0 --u5_kappa 0 --u5_chon gc --u5_er --u5_tam", "td": "--cach tu_dong",
      "nn1x": "--cach ngau_nhien --so_lo 3 --K 6", "unc8": "--cach u_nhiet --un_vong_lon 1102 --un_vong_nho 3222 @SA_D5@"}
DT = {"d5": lambda c: f"{B}/n14/dac_ta/dong_vat_n14_{c}.json", "k8": lambda c: f"{B}/n14d5/dac_ta/dong_vat_n14k08_{c}.json"}
SO = json.load(open(f"{D}/so_luong_n14d5.json"))["so_luong"]
GB = json.load(open(f"{D}/gb_n14d5.json")) if os.path.exists(f"{D}/gb_n14d5.json") else {"798836": 4, "989830": 7, "593409": 6, "707681": 7, "997092": 16}
def dong(i, gb, cmd): return f"{i} C {gb} u5 1 {ENV} {PY} {cmd} --out {KQ}/{i}"
def viec(v, c, nh, s, pre):
    i = f"u5_n14d5_{pre}{v}_{c}_{nh}_{s}"; n = SO[v][c]
    return dong(i, GB[c], f"{MA} --dac_ta {DT[v](c)} --seed {s} --N1 128 --K 4 --thu_vien T --so ti_the={n['ti_the']} --so lysosome={n['lysosome']} "
                         f"{NH[nh]} --tang2 --khong_tu_do --eval_dung_tu_dau --hd_truoc_sau")
G = {}
CU = {"td": "--cach tu_dong", "nn1x": "--cach ngau_nhien --so_lo 3 --K 17",
      "un1x": "--cach u_nhiet --un_vong_lon 976 --un_vong_nho 2982 --un_T1 1e-6 --un_p0 0.5 --un_buoc_t 0.3 0.02 --un_buoc_r 0.5 0.03 --un_ti_le_lon 0.4 --un_p_doi 0.5 --un_lam_moi 50"}
G["id"] = [dong(f"u5_n14d5_id_798836_{m}_6501", 4, f"{MA} --dac_ta {B}/n14/dac_ta/dong_vat_n14_798836.json --seed 6501 --N1 128 --K 4 --thu_vien T --so ti_the=23 --so lysosome=8 "
                f"{CU[m]} --tang2 --khong_tu_do --eval_dung_tu_dau") for m in ("td", "nn1x", "un1x")]
G["id"].append(G["id"][0].replace("_td_6501", "_tdhd_6501").replace(" --out ", " --hd_truoc_sau --out "))
G["thu"] = [viec(v, c, nh, 6549, "thu_") for v in ("d5", "k8") for nh in NH for c in THU]
for v in ("d5", "k8"):
    G[f"{v}_1"] = [viec(v, c, nh, 6541, "") for nh in NH for c in CELL]
    G[f"{v}_23"] = [viec(v, c, nh, s, "") for s in (6542, 6543) for nh in NH for c in CELL]
for k, L in G.items():
    open(f"{D}/viec_{k}.txt", "w").write("\n".join(L) + "\n"); print(k, len(L))
print("tong", sum(len(L) for L in G.values()))
