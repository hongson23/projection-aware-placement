# -*- coding: utf-8 -*-
"""chon_dev.py (22/11/2026, bản 3: trần ngân sách theo SỐ ĐẾM bước tương đương — PI duyệt 22/11 09:15 "ok"; giây chỉ mô tả; bản trước: chon_dev.py.truoc_dem)
bản 2 (21/11/2026, sau rà soát) — LUẬT CHỌN GHI TRƯỚC cho CX (thiết kế §6), seed dev 6597, 6598, chỉ số phía bộ dựng.
Thước mỗi việc (kq/<id>/thong_tin.json):
  G_HD = Σ_{mọi trường kể cả ER} (HD_q − h)₊ (HD hai chiều mesh trước↔sau chiếu) [thước chọn; rà soát: Spearman 0,99 với G_c^{K+ER}, 90/90 dấu cặp]
  HDmax; Vmat = 100·Σ(V_trước − V_sau)₊/ΣV_trước (kể cả ER); mat = |mat_bo_phan|; d50 = |con_lai_duoi_50pt|; ti_gio = giay_tong / giay_tong(U5b cùng lô, cùng seed)
Mỗi nhóm g ∈ {best, eq} (luoi_dev.json):
  (1) đủ 2 seed; (2) giờ: trung vị ti_gio ≤ trần_g (best: 1,05·r_P4·1,05 ≈ giờ P4 + 10 %; eq: 1,10);
  (3) an toàn: Σ mat ≤ Σ mat(U5b) + 1, Σ d50 ≤ Σ d50(U5b) + 2;
  (4) nhỏ nhất trung vị G_HD; trong 5 % / 0,005 là hoà ⇒ HDmax (2 %) ⇒ Vmat (2 %) ⇒ Bh, B, A ⇒ L1 trước L2 ⇒ cx_l0 nhỏ.
  Không cấu hình nào qua (2) ở nhóm g ⇒ g = None (BẤT KHẢ, ghi rõ; KHÔNG lấy cấu hình yếu nhất thay).
Cổng 'CX có hiệu lực' (báo, không chặn): trung vị G_HD(CX_best) < trung vị G_HD(P4 cùng lô). Báo kèm: 'tran' (bỏ qua giờ), tất định GPU (U5b cùng lô so u5_v2E_K4_6597).
Ra: u5x/chon_dev.json + chon_dev.md."""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import dem_buoc as DB
R = os.environ.get("CX_R", os.path.expanduser("~/Vision_team_data_189/sonhh20/ch3_tebao_dep")); O = f"{R}/bo_dung_chung/tang_dg"; U = f"{O}/u5x"; H = 0.1796875
SEEDS = (6597, 6598)
def thuoc(d):
    hd = d["ban_u5b"]["hd_truoc_sau"] or {}; v = [max(x) for k, x in hd.items() if not k.startswith("_") and x]
    tt = d["the_tich_truoc_sau"]; vb = sum(a for a, b in tt.values()); vm = sum(max(a - b, 0.0) for a, b in tt.values())
    cx = (d.get("bao_cao_bo_tri") or {}).get("cx") or {}
    return {"G_HD": float(sum(max(x - H, 0.0) for x in v)), "HDmax": float(max(v)) if v else None, "Vmat": 100 * vm / max(vb, 1e-12),
            "mat": len(d.get("mat_bo_phan") or []), "d50": len(d.get("con_lai_duoi_50pt") or []), "giay": d["giay_tong"], "k_chon": d["bao_cao_bo_tri"].get("cau_hinh_chon"),
            "n_su_kien": cx.get("n_su_kien"), "lam_cuoi": [(c_.get("lam_lon"), c_.get("lam_nho")) for c_ in cx.get("cau_hinh", [])]}
def doc(id_):
    f = f"{O}/kq/{id_}/thong_tin.json"
    if not (os.path.exists(f) and os.path.exists(f"{O}/xong/{id_}")): return None
    d = json.load(open(f)); r = thuoc(d); r["dem"] = DB.chi_phi(d, HS); r["sach_loi"] = DB.sach(id_); return r
LU = json.load(open(f"{U}/luoi_dev.json")); MM = LU["m_theo_muc"]; HS = MM["hang_so_buoc"]
S_P8 = {s: DB.so_buoc(json.load(open(f"{O}/kq/u5_v2P8_K8_{s}/thong_tin.json"))) for s in SEEDS}
ref = {s: doc(f"u5_x_dev_U5b_{s}") for s in SEEDS}; p4 = {s: doc(f"u5_x_dev_P4_{s}") for s in SEEDS}
assert all(ref.values()), "thiếu U5b tham chiếu cùng lô"
TRAN = {"best": 1.10, "eq": 1.10}
bang = []
for c, p in LU["cfg"].items():
    r = {s: doc(f"u5_x_dev_{c}_{s}") for s in SEEDS}; row = {"cfg": c, **p, "du_seed": all(r.values()), "theo_seed": r}
    if row["du_seed"]:
        for k in ("G_HD", "HDmax", "Vmat"): row[k] = float(np.median([r[s][k] for s in SEEDS]))
        row["ti_gio"] = float(np.median([r[s]["giay"] / ref[s]["giay"] for s in SEEDS]))
        row["ti_gio_p4"] = float(np.median([r[s]["giay"] / p4[s]["giay"] for s in SEEDS])) if all(p4.values()) else None
        row["ti_dem"] = float(np.median([r[s]["dem"] / ref[s]["dem"] for s in SEEDS])); row["ti_dem_p8"] = float(np.median([r[s]["dem"] / S_P8[s] for s in SEEDS]))
        row["mat"] = sum(r[s]["mat"] for s in SEEDS); row["d50"] = sum(r[s]["d50"] for s in SEEDS)
        row["an_toan"] = row["mat"] <= sum(ref[s]["mat"] for s in SEEDS) + 1 and row["d50"] <= sum(ref[s]["d50"] for s in SEEDS) + 2
    bang.append(row)
UU = {"Bh": 0, "B": 1, "A": 2}
def chon(nhom):
    gio = (lambda x: x["ti_dem_p8"] <= TRAN["best"]) if nhom == "best" else (lambda x: x["ti_dem"] <= TRAN["eq"])
    ds = [x for x in bang if nhom in x["nhom"] and x["du_seed"] and x["an_toan"] and gio(x)]
    if not ds: return None
    g0 = min(x["G_HD"] for x in ds); ds = [x for x in ds if x["G_HD"] <= max(g0 * 1.05, g0 + 0.005)]
    h0 = min(x["HDmax"] for x in ds); ds = [x for x in ds if x["HDmax"] <= h0 * 1.02]
    v0 = min(x["Vmat"] for x in ds); ds = [x for x in ds if x["Vmat"] <= v0 * 1.02 + 1e-9]
    x = sorted(ds, key=lambda x: (UU[x["neo"]], x["muc"] != "L1", x["l0"]))[0]
    return {k: x[k] for k in ("cfg", "neo", "muc", "m", "l0", "G_HD", "HDmax", "Vmat", "ti_dem", "ti_dem_p8", "ti_gio")}
best, eq = chon("best"), chon("eq")
assert best is not None or eq is not None, "không cấu hình nào (best, eq) qua luật — dừng, không chèn lượt chính"
p4m = float(np.median([p4[s]["G_HD"] for s in SEEDS])) if all(p4.values()) else None
tran = [x for x in bang if "tran" in x["nhom"] and x["du_seed"]]
tat_dinh = None
f0 = f"{O}/kq/u5_v2E_K4_6597/thong_tin.json"
if os.path.exists(f0):
    a = thuoc(json.load(open(f0))); tat_dinh = {"G_HD": [a["G_HD"], ref[6597]["G_HD"]], "HDmax": [a["HDmax"], ref[6597]["HDmax"]], "giay": [a["giay"], ref[6597]["giay"]]}
out = {"best": best, "eq": eq, "CUNG_GIO_BAT_KHA": eq is None, "BEST_BAT_KHA": best is None, "CX_HIEU_LUC_DEV": (p4m is not None and (best or eq)["G_HD"] < p4m), "P4_G_HD_tv": p4m,
       "tran": [{k: x.get(k) for k in ("cfg", "G_HD", "HDmax", "Vmat", "ti_gio")} for x in tran], "tat_dinh_GPU_U5b_6597": tat_dinh, "U5b_ref": ref, "P4_ref": p4, "tran_gio": TRAN, "bang": bang}
json.dump(out, open(f"{U}/chon_dev.json", "w"), indent=1, ensure_ascii=False)
L = ["# Chọn cấu hình CX trên dev (luật ghi trước, chon_dev.py) — CHƯA rà soát chéo", "", f"best: **{best['cfg'] if best else 'BẤT KHẢ'}**; eq: **{eq['cfg'] if eq else 'BẤT KHẢ'}**; CX có hiệu lực trên dev (G_HD best < P4 {p4m}): {out['CX_HIEU_LUC_DEV']}", "",
     "| cfg | nhóm | G_HD tv | HDmax tv | Vmat % tv | ti_dem tv (vs U5b) | ti_dem tv (vs P8) | ti_gio tv (mô tả) | mat | d50 | an toàn |", "|---|---|---|---|---|---|---|---|---|---|---|"]
for x in bang:
    if x["du_seed"]: L.append(f"| {x['cfg']} | {','.join(x['nhom'])} | {x['G_HD']:.4f} | {x['HDmax']:.4f} | {x['Vmat']:.3f} | {x['ti_dem']:.3f} | {x['ti_dem_p8']:.3f} | {x['ti_gio']:.3f} | {x['mat']} | {x['d50']} | {x['an_toan']} |")
    else: L.append(f"| {x['cfg']} | {','.join(x['nhom'])} | thiếu seed | | | | | | |")
L += ["", "U5b cùng lô: " + json.dumps(ref, ensure_ascii=False), "P4 cùng lô: " + json.dumps(p4, ensure_ascii=False), "Tất định GPU (u5_v2E_K4_6597 vs U5b cùng lô): " + json.dumps(tat_dinh)]
open(f"{U}/chon_dev.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))
