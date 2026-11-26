# -*- coding: utf-8 -*-
"""phan_tich_rong.py (20/11/2026) — tổng hợp lượt rộng W1 từ u5v2/doc_lap/u5_r_*.json (phan_tich_doc_lap_u5.py, bộ kiểm v3 đóng băng + số bộ dựng).
Mỗi vùng × cặp nhánh × thước (đều: nhỏ hơn = tốt hơn): số seed thắng/thua/hoà, trung vị hiệu (A − B), trung vị hiệu tương đối, p Wilcoxon hai phía (cặp theo seed).
Ra: u5v2/ket_qua_rong.json + .md (CHƯA rà soát chéo)."""
import os, json, glob, collections, math
import numpy as np
U = os.environ.get("RA_RONG", os.path.expanduser("~/Vision_team_data_189/sonhh20/ch3_tebao_dep/bo_dung_chung/tang_dg/u5v2")); D = os.environ.get("DOC_LAP", U + "/doc_lap")
try:
    from scipy.stats import wilcoxon
except Exception: wilcoxon = None
def g(r, *k):
    for x in k:
        if r is None: return None
        r = r.get(x) if isinstance(r, dict) else None
    return r
H = 0.1797
def _gcer(r):
    gk, er = r.get("Gc_K"), g(r, "hd_ts_theo_loai", "luoi_noi_chat", "max")
    if gk is None or er is None or gk != gk: return None
    return gk + max(er - H, 0.0)
THUOC = {
    "HD_max_tat_ca": lambda r: g(r, "hd_ts_tat_ca", "max"),
    "Gc_K_ER": _gcer,
    "Gc_K": lambda r: r.get("Gc_K"),
    "V_test": lambda r: r.get("V_test"),
    "V_test_tru_43b_duoi_luoi": lambda r: r.get("V_test_tru_43b_duoi_luoi"),
    "n_mat_kiem": lambda r: r.get("n_mat_kiem"),
    "HD_tu_max": lambda r: g(r, "hd_ts_tu_the", "max"), "HD_tu_p95": lambda r: g(r, "hd_ts_tu_the", "p95"), "HD_tu_p50": lambda r: g(r, "hd_ts_tu_the", "p50"),
    "HD_tu_tong": lambda r: g(r, "hd_ts_tu_the", "tong"),
    "HD_ER": lambda r: g(r, "hd_ts_theo_loai", "luoi_noi_chat", "max"),
    "HD_nhan": lambda r: g(r, "hd_ts_theo_loai", "nhan", "max"), "HD_golgi": lambda r: g(r, "hd_ts_theo_loai", "golgi", "max"),
    **{f"HD_{l}_{q}": (lambda r, l=l, q=q: g(r, "hd_ts_theo_loai", l, q)) for l in ("luoi_noi_chat", "nhan", "golgi", "trung_the", "nhan_con", "ti_the", "lysosome") for q in ("max", "p95")},
    "n_43b": lambda r: (r.get("luat_truot") or {}).get("4.3b", {}).get("so_ca_truot", 0), "n_43b_duoi_luoi": lambda r: r.get("n_43b_duoi_luoi"),
    "Vmat_tat_ca_pt": lambda r: 100 * r["V_mat_ti_le_tat_ca"], "Vmat_tu_the_pt": lambda r: 100 * r["V_mat_ti_le_tu_the"],
    "Vmat_ER_pt": lambda r: 100 * g(r, "theo_loai", "luoi_noi_chat", "V_mat_ti_le"),
    "vo_manh": lambda r: len(r.get("vo_manh") or []), "con_lai_duoi_50pt": lambda r: len(r.get("con_lai_duoi_50pt") or []),
    "giay_tong": lambda r: r.get("giay_tong"), "giay_dat": lambda r: r.get("giay_dat")}
MO_TA_THEM = {"V_nhan_truoc": lambda r: r.get("V_nhan_truoc"), "k_chon": lambda r: r.get("k_chon")}
CAP = [("U5b", "P8"), ("U5b", "P4"), ("U5b", "a1"), ("U5b", "l0"), ("U5b8", "P8"), ("U5b8", "U5b")]
B = collections.defaultdict(dict)
for f in glob.glob(D + "/u5_r_*.json"):
    if f.endswith("_hd.json"): continue
    r = json.load(open(f)); p = os.path.basename(f)[:-5].split("_")
    B[(p[2], p[3])][int(p[4])] = r
out = {"cap": {}, "mo_ta": {}}
for (nh, vg), d in sorted(B.items()):
    out["mo_ta"][f"{nh}|{vg}"] = {m: (lambda xs: {"n": len(xs), "tv": float(np.median(xs)) if xs else None, "min": min(xs) if xs else None, "max": max(xs) if xs else None})(
        [v for v in (fn(r) for r in d.values()) if v is not None]) for m, fn in THUOC.items()}
for vg in sorted({v for _, v in B}):
    for a, b in CAP:
        A, Bb = B.get((a, vg), {}), B.get((b, vg), {}); ss = sorted(set(A) & set(Bb))
        if not ss: continue
        thieu = sorted((set(A) | set(Bb) | set(range(6511, 6521))) - set(ss))
        cung_k = [s_ for s_ in ss if A[s_].get("k_chon") == Bb[s_].get("k_chon") and abs((A[s_].get("V_nhan_truoc") or 0) - (Bb[s_].get("V_nhan_truoc") or 0)) < 1e-3]
        res = {}
        for m, fn in THUOC.items():
            pr = [(fn(A[s]), fn(Bb[s])) for s in ss]; pr = [(x, y) for x, y in pr if x is not None and y is not None]
            if not pr: continue
            dd = np.array([x - y for x, y in pr], float); rel = [(x - y) / max(abs(y), 1e-9) for x, y in pr]
            p = None
            if wilcoxon is not None and np.any(np.abs(dd) > 1e-12):
                try: p = float(wilcoxon(dd, zero_method="wilcox", alternative="two-sided").pvalue)
                except Exception: p = None
            res[m] = {"n": len(pr), "thang": int((dd < -1e-9).sum()), "thua": int((dd > 1e-9).sum()), "hoa": int((np.abs(dd) <= 1e-9).sum()),
                      "tv_hieu": float(np.median(dd)), "tv_hieu_tuong_doi": float(np.median(rel)), "p_wilcoxon": p}
        out["cap"][f"{vg}|{a}-{b}"] = {"seed": ss, "seed_thieu": thieu, "seed_cung_cau_hinh": cung_k, "thuoc": res}
def thang_vung(a, b, vg):
    A, Bb = B.get((a, vg), {}), B.get((b, vg), {})
    if not A or not Bb: return None
    def tv(D, m):
        xs = [THUOC[m](r) for r in D.values()]; xs = [x for x in xs if x is not None and x == x]
        return float(np.median(xs)) if xs else None
    hA, hB, vA, vB = tv(A, "HD_max_tat_ca"), tv(Bb, "HD_max_tat_ca"), tv(A, "Vmat_tat_ca_pt"), tv(Bb, "Vmat_tat_ca_pt")
    mA = sum((r.get("n_mat_kiem") or 0) for r in A.values()); mB = sum((r.get("n_mat_kiem") or 0) for r in Bb.values())
    ok = None if None in (hA, hB, vA, vB) else bool(hA < hB and vA < vB and mA <= mB)
    return {"HD_max_tv": [hA, hB], "Vmat_pt_tv": [vA, vB], "mat_tong": [mA, mB], "seed": [sorted(A), sorted(Bb)], "THANG_VUNG": ok}
out["quyet_dinh"] = {f"{vg}|{a}-{b}": thang_vung(a, b, vg) for vg in ("k08", "D5", "D4A2") for a, b in CAP}
json.dump(out, open(U + "/ket_qua_rong.json", "w"), indent=1, ensure_ascii=False)
L = ["# Lượt rộng W1 — tổng hợp tự động (CHƯA rà soát chéo)", "", "thắng = nhánh trái nhỏ hơn (tốt hơn) ở seed đó.", "", "Luật: A thắng vùng ⇔ trung vị HD_max (kể cả ER) nhỏ hơn VÀ trung vị % thể tích mất (kể cả ER) nhỏ hơn VÀ tổng bộ phận mất không nhiều hơn.", "", "\n".join(f"- {k}: {v}" for k, v in out["quyet_dinh"].items()), ""]
for k, v in out["cap"].items():
    L.append(f"## {k} (n = {len(v['seed'])}; thiếu {v['seed_thieu']}; cùng cấu hình {v['seed_cung_cau_hinh']})"); L.append("| thước | thắng/thua/hoà | trung vị hiệu | tv tương đối | p |"); L.append("|---|---|---|---|---|")
    for m, x in v["thuoc"].items():
        L.append(f"| {m} | {x['thang']}/{x['thua']}/{x['hoa']} | {x['tv_hieu']:.4g} | {100 * x['tv_hieu_tuong_doi']:.1f} % | {x['p_wilcoxon'] if x['p_wilcoxon'] is None else round(x['p_wilcoxon'], 4)} |")
    L.append("")
open(U + "/ket_qua_rong.md", "w").write("\n".join(L)); print("\n".join(L[:60]))
