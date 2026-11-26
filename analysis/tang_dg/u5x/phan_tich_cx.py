# -*- coding: utf-8 -*-
"""phan_tich_cx.py (21/11/2026, bản 2 sau rà soát) — U5b so với đối chứng kiểu ComPLx ở k08, seed 6511–6520, từ u5v2/doc_lap/*.json
(phan_tich_doc_lap_u5.py: bộ kiểm v3 đóng băng qua khoá mù + số bộ dựng). LUẬT GHI TRƯỚC (thiết kế §8), viết trước khi có số CX nào.
Mỗi bản CX ∈ {CXb (ngân sách ≤ giờ P4), CXe (cùng giờ, nếu có)}, thước chính G_c^{K+ER} = Gc_K + (HD_ER − h)₊.
Seed thiếu/lỗi ở BẤT KỲ bên nào ⇒ giá trị XẤU NHẤT cho U5b: không thắng; giảm tương đối = −1; ΔV = +∞; số mất U5b = 99.
  Cổng hiệu lực: CX < P4 ở ≥ 6/10 seed VÀ trung vị giảm tương đối (P4 → CX) ≥ 10 %; không qua ⇒ CX_KHONG_HIEU_LUC
                 (không được viết "U5 beats a ComPLx-style baseline"; chỉ được viết "a ComPLx-style anchoring baseline did not improve over PLACE here").
  U5_HON_CX ⇔ (a) U5b < CX ở ≥ 8/10; (b) trung vị giảm ≥ 30 %; (c) trung vị ΔV (U5b − CX, % kể cả ER) ≤ +0,5 điểm %; (d) Σ mất U5b ≤ CX + 1.
  CX_NGANG_HOAC_HON ⇔ U5b thắng ≤ 5/10 HOẶC trung vị giảm < 10 %; còn lại CHUA_RO.
  [22/11, PI duyệt: r tính theo SỐ ĐẾM bước tương đương (dem_buoc.py, hằng số u5x/m_theo_muc.json); tỉ số giây chỉ mô tả]
  Nhãn giờ r = trung vị bước_tương_đương(CX)/bước(U5b): r < 0,90 ⇒ CX_THIEU_GIO (không được kết luận); 0,90 ≤ r ≤ 1,10 ⇒ CUNG_GIO;
  r > 1,10 ⇒ CX_NHIEU_GIO_HON (U5_HON_CX ⇒ "U5b wins even with CX given r× the wall-clock"; CX_NGANG_HOAC_HON ⇒ CX_NGANG_NHUNG_TON_GIO).
Báo kèm: HD max/p95 kể cả ER, theo loại, % V, giờ, số sự kiện / λ cuối / Π mỗi việc CX, k_chon trùng P4, trôi giờ U5b (u5_x_U5bcal vs u5_r_U5b).
Ra: u5x/ket_qua_cx.json + .md (CHƯA rà soát chéo)."""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import dem_buoc as DB
O = os.environ.get("CX_O", os.path.expanduser("~/Vision_team_data_189/sonhh20/ch3_tebao_dep/bo_dung_chung/tang_dg")); D = os.environ.get("DOC_LAP", f"{O}/u5v2/doc_lap"); U = f"{O}/u5x"
H = 0.1797; SEEDS = list(range(6511, 6521))
try:
    from scipy.stats import wilcoxon
except Exception: wilcoxon = None
def g(r, *k):
    for x in k:
        if r is None: return None
        r = r.get(x) if isinstance(r, dict) else None
    return r
def gcer(r):
    if not r: return None
    gk, er = r.get("Gc_K"), g(r, "hd_ts_theo_loai", "luoi_noi_chat", "max")
    return None if gk is None or er is None or gk != gk else gk + max(er - H, 0.0)
THUOC = {"Gc_K_ER": gcer, "Gc_K": lambda r: r.get("Gc_K"), "HD_max_tat_ca": lambda r: g(r, "hd_ts_tat_ca", "max"), "HD_p95_tat_ca": lambda r: g(r, "hd_ts_tat_ca", "p95"),
         "HD_ER": lambda r: g(r, "hd_ts_theo_loai", "luoi_noi_chat", "max"), "HD_nhan": lambda r: g(r, "hd_ts_theo_loai", "nhan", "max"), "HD_golgi": lambda r: g(r, "hd_ts_theo_loai", "golgi", "max"),
         "HD_ti_the_max": lambda r: g(r, "hd_ts_theo_loai", "ti_the", "max"), "HD_lysosome_max": lambda r: g(r, "hd_ts_theo_loai", "lysosome", "max"),
         "Vmat_tat_ca_pt": lambda r: 100 * r["V_mat_ti_le_tat_ca"], "Vmat_tu_the_pt": lambda r: 100 * r["V_mat_ti_le_tu_the"],
         "Vmat_ER_pt": lambda r: 100 * (g(r, "theo_loai", "luoi_noi_chat", "V_mat_ti_le") or 0), "n_mat_kiem": lambda r: r.get("n_mat_kiem"),
         "con_lai_duoi_50pt": lambda r: len(r.get("con_lai_duoi_50pt") or []), "V_test_tru_43b_duoi_luoi": lambda r: r.get("V_test_tru_43b_duoi_luoi"),
         "giay_tong": lambda r: r.get("giay_tong"), "giay_dat": lambda r: r.get("giay_dat"), "giay_chieu": lambda r: r.get("giay_chieu")}
def doc(ten, s):
    f = f"{D}/{ten}_{s}.json"; return json.load(open(f)) if os.path.exists(f) else None
NH = {"U5b": "u5_r_U5b_k08", "P4": "u5_r_P4_k08", "P8": "u5_r_P8_k08", "CXb": "u5_x_CXb_k08", "CXe": "u5_x_CXe_k08"}
B = {n: {s: doc(p, s) for s in SEEDS} for n, p in NH.items()}
CX_CO = [n for n in ("CXb", "CXe") if any(B[n].values())]
def cap(a, b):
    res = {}
    for m, fn in THUOC.items():
        pr = [(s, fn(B[a][s]), fn(B[b][s])) for s in SEEDS if B[a][s] and B[b][s]]; pr = [(s, x, y) for s, x, y in pr if x is not None and y is not None]
        if not pr: continue
        dd = np.array([x - y for _, x, y in pr]); rel = [(x - y) / max(abs(y), 1e-9) for _, x, y in pr]; p = None
        if wilcoxon is not None and np.any(np.abs(dd) > 1e-12):
            try: p = float(wilcoxon(dd).pvalue)
            except Exception: p = None
        res[m] = {"n": len(pr), "thang": int((dd < -1e-9).sum()), "thua": int((dd > 1e-9).sum()), "hoa": int((np.abs(dd) <= 1e-9).sum()),
                  "tv_a": float(np.median([x for _, x, _ in pr])), "tv_b": float(np.median([y for _, _, y in pr])), "tv_hieu": float(np.median(dd)), "tv_tuong_doi": float(np.median(rel)), "p": p,
                  "theo_seed": {s: [x, y] for s, x, y in pr}}
    return res
def so_chinh(a, b, x_xau=True):
    """thắng / giảm tương đối trên G_c^{K+ER}, seed thiếu = a thua (rel −1)."""
    th, rel = 0, []
    for s in SEEDS:
        x, y = gcer(B[a][s]), gcer(B[b][s])
        if x is None or y is None: rel.append(-1.0); continue
        if x < y - 1e-9: th += 1
        rel.append((y - x) / max(abs(y), 1e-9))
    return th, float(np.median(rel))
def thong_ke_cx(n):
    out = {}
    for s in SEEDS:
        f = f"{O}/kq/{NH[n]}_{s}/thong_tin.json"
        if not os.path.exists(f): continue
        b = json.load(open(f))["bao_cao_bo_tri"]; c = b.get("cx") or {}; ev = c.get("su_kien") or []
        out[s] = {"n_su_kien": c.get("n_su_kien"), "giay_cx": c.get("giay"), "lam_cuoi": [(x.get("lam_lon"), x.get("lam_nho")) for x in c.get("cau_hinh", [])],
                  "Pi_dau_cuoi_nho": [ev_["Pi"] for ev_ in ev if ev_["pha"] == "nho"][:1] + [ev_["Pi"] for ev_ in ev if ev_["pha"] == "nho"][-1:], "k_chon": b.get("cau_hinh_chon"),
                  "k_chon_trung_P4": (B["P4"][s] or {}).get("k_chon") == b.get("cau_hinh_chon")}
    return out
out = {"cap": {}, "quyet_dinh": {}, "cx_moi_viec": {}}
for n in CX_CO:
    for a, b in (("U5b", n), (n, "P4"), (n, "P8")): out["cap"][f"{a}-{b}"] = cap(a, b)
    th_p4, giam_p4 = so_chinh(n, "P4"); hieu_luc = th_p4 >= 6 and giam_p4 >= 0.10
    th, giam = so_chinh("U5b", n)
    dV = [(THUOC["Vmat_tat_ca_pt"](B["U5b"][s]) - THUOC["Vmat_tat_ca_pt"](B[n][s])) if (B["U5b"][s] and B[n][s]) else float("inf") for s in SEEDS]
    mat = [sum(((B["U5b"][s] or {}).get("n_mat_kiem") if B["U5b"][s] and B[n][s] else 99) or 0 for s in SEEDS), sum((B[n][s] or {}).get("n_mat_kiem") or 0 for s in SEEDS)]
    tg = [B[n][s]["giay_tong"] / B["U5b"][s]["giay_tong"] for s in SEEDS if B["U5b"][s] and B[n][s]]; r_giay = float(np.median(tg)) if tg else None
    HS = json.load(open(f"{U}/m_theo_muc.json"))["hang_so_buoc"]; kq_ = lambda i: json.load(open(f"{O}/kq/{i}/thong_tin.json"))
    td = [DB.chi_phi(kq_(f"{NH[n]}_{s}"), HS) / DB.chi_phi_u5(kq_(f"{NH['U5b']}_{s}"), HS) for s in SEEDS if B["U5b"][s] and B[n][s]]; r = float(np.median(td)) if td else None
    dk = {"a_thang_ge8": th >= 8, "b_giam_ge30": giam >= 0.30, "c_dV_le_0.5": float(np.median(dV)) <= 0.5, "d_mat": mat[0] <= mat[1] + 1}
    kl = "U5_HON_CX" if all(dk.values()) else ("CX_NGANG_HOAC_HON" if th <= 5 or giam < 0.10 else "CHUA_RO")
    gio = None if r is None else ("CX_THIEU_GIO" if r < 0.90 else ("CUNG_GIO" if r <= 1.10 else "CX_NHIEU_GIO_HON"))
    if kl == "CX_NGANG_HOAC_HON" and gio == "CX_NHIEU_GIO_HON": kl = "CX_NGANG_NHUNG_TON_GIO"
    if not hieu_luc: kl = "CX_KHONG_HIEU_LUC (" + kl + ")"
    out["quyet_dinh"][n] = {"ket_luan": kl, "nhan_gio": gio, "ti_dem_CX_U5b_tv": r, "ti_giay_CX_U5b_tv_mo_ta": r_giay, "sach_loi": {s: DB.sach(f"{NH[n]}_{s}") for s in SEEDS}, "dieu_kien": dk, "thang": th, "giam_tv": giam, "dV_tv": float(np.median(dV)), "mat": mat,
                            "hieu_luc": {"thang_vs_P4": th_p4, "giam_vs_P4": giam_p4, "dat": hieu_luc}, "seed_thieu": [s for s in SEEDS if not (B["U5b"][s] and B[n][s])]}
    out["cx_moi_viec"][n] = thong_ke_cx(n)
cal = {}
for s in (6511, 6512):
    f1, f2 = f"{O}/kq/u5_x_U5bcal_k08_{s}/thong_tin.json", f"{O}/kq/u5_r_U5b_k08_{s}/thong_tin.json"
    if os.path.exists(f1) and os.path.exists(f2): cal[s] = [json.load(open(f1))["giay_tong"], json.load(open(f2))["giay_tong"]]
out["troi_gio_U5b"] = {"theo_seed": cal, "ti_tv": float(np.median([a / b for a, b in cal.values()])) if cal else None}
json.dump(out, open(f"{U}/ket_qua_cx.json", "w"), indent=1, ensure_ascii=False)
L = ["# U5b so với đối chứng kiểu ComPLx (CX), k08 — tổng hợp tự động (CHƯA rà soát chéo)", "", "```", json.dumps({"quyet_dinh": out["quyet_dinh"], "troi_gio_U5b": out["troi_gio_U5b"]}, ensure_ascii=False, indent=1), "```", ""]
for k, v in out["cap"].items():
    L += [f"## {k}", "| thước | thắng/thua/hoà (trái nhỏ hơn) | tv trái | tv phải | tv hiệu | tv tương đối | p |", "|---|---|---|---|---|---|---|"]
    for m, x in v.items(): L.append(f"| {m} | {x['thang']}/{x['thua']}/{x['hoa']} | {x['tv_a']:.4g} | {x['tv_b']:.4g} | {x['tv_hieu']:.4g} | {100 * x['tv_tuong_doi']:.1f} % | {x['p'] if x['p'] is None else round(x['p'], 4)} |")
    L.append("")
L += ["## Mỗi việc CX", "```", json.dumps(out["cx_moi_viec"], ensure_ascii=False, indent=1), "```"]
open(f"{U}/ket_qua_cx.md", "w").write("\n".join(L) + "\n"); print("\n".join(L[:40]))
