# -*- coding: utf-8 -*-
"""phan_tich_u5c.py (21/11/2026) — bảng kiểu CG cho lượt u5_c (CHƯA rà soát chéo). Thước lấy NGUYÊN từ phan_tich_rong.py (md5 b4f56e61,
đoạn trước "CAP = [", exec), bản ghi u5v2/doc_lap/<id>.json (phan_tich_doc_lap_u5.py; e3b_B_* thêm hd_ts_* bởi hd_e3b.py).
Nhánh: u5_r_<nh>_<vùng>_<seed>, u5_c_<nh>_<vùng>_<seed>, e3b_B_<ô>_{td→P4, nn1x→TL, un1x→SA}_<seed>.
Cặp (trái nhỏ hơn = tốt hơn, cặp theo seed): k08: U5b so a1l03/a1/a1l3/ah/wv, mỗi đối chứng so P8, P4; ô N1: U5b so P4/TL/SA.
Luật vùng như phan_tich_rong (HD_max, %V mất, Σ mất). Ra: u5v2/u5c/ket_qua_u5c.json + .md"""
import os, json, glob, collections, hashlib
import numpy as np
U = os.path.expanduser("~/Vision_team_data_189/sonhh20/ch3_tebao_dep/bo_dung_chung/tang_dg/u5v2"); D = U + "/doc_lap"
src = open(U + "/phan_tich_rong.py").read(); exec(src.split("\nCAP = [")[0].split("import numpy as np", 1)[1])
D = U + "/doc_lap"
B = collections.defaultdict(dict); MAP = {"td": "P4", "nn1x": "TL", "un1x": "SA"}
for f in glob.glob(D + "/*.json"):
    if f.endswith("_hd.json"): continue
    n = os.path.basename(f)[:-5]; p = n.split("_")
    if n.startswith("u5_r_") or n.startswith("u5_c_"): k = (p[2], p[3]); s = int(p[4])
    elif n.startswith("e3b_B_") and p[3] in MAP: k = (MAP[p[3]], p[2]); s = int(p[4])
    else: continue
    B[k][s] = json.load(open(f))
CHINH = ["HD_max_tat_ca", "Gc_K_ER", "Gc_K", "Vmat_tat_ca_pt", "Vmat_tu_the_pt", "n_mat_kiem", "con_lai_duoi_50pt", "V_test", "HD_ER", "HD_nhan", "HD_golgi", "HD_ti_the_max", "HD_lysosome_max", "giay_dat", "giay_tong"]
CAP = [("k08", a, b) for a, b in [("U5b", "a1l03"), ("U5b", "a1"), ("U5b", "a1l3"), ("U5b", "ah"), ("U5b", "wv"),
                                  ("a1l03", "P8"), ("a1", "P8"), ("a1l3", "P8"), ("ah", "P8"), ("wv", "P8"), ("U5b", "P8"),
                                  ("a1l03", "P4"), ("a1l3", "P4"), ("ah", "P4"), ("wv", "P4"), ("U5b", "P4")]]
CAP += [(o, "U5b", b) for o in ("D1", "D4", "D4A2", "D5") for b in ("P4", "TL", "SA")]
def tv(xs): xs = [x for x in xs if x is not None and x == x]; return float(np.median(xs)) if xs else None
out = {"cap": {}, "dong_nhat": None, "md5_phan_tich_rong": hashlib.md5(src.encode()).hexdigest()}
for vg, a, b in CAP:
    A, Bb = B.get((a, vg), {}), B.get((b, vg), {}); ss = sorted(set(A) & set(Bb))
    if not ss: continue
    res = {}
    for m in CHINH:
        fn = THUOC[m]; pr = [(fn(A[s]), fn(Bb[s])) for s in ss]; pr = [(x, y) for x, y in pr if x is not None and y is not None and x == x and y == y]
        if not pr: continue
        dd = np.array([x - y for x, y in pr], float); rel = [(x - y) / max(abs(y), 1e-9) for x, y in pr]; p = None
        if wilcoxon is not None and np.any(np.abs(dd) > 1e-12):
            try: p = float(wilcoxon(dd, zero_method="wilcox", alternative="two-sided").pvalue)
            except Exception: p = None
        res[m] = {"n": len(pr), "tv_A": tv([x for x, _ in pr]), "tv_B": tv([y for _, y in pr]), "thang": int((dd < -1e-9).sum()), "thua": int((dd > 1e-9).sum()),
                  "hoa": int((np.abs(dd) <= 1e-9).sum()), "tv_hieu_tuong_doi": float(np.median(rel)), "p_wilcoxon": p}
    hA, hB = tv([THUOC["HD_max_tat_ca"](r) for r in A.values()]), tv([THUOC["HD_max_tat_ca"](r) for r in Bb.values()])
    vA, vB = tv([THUOC["Vmat_tat_ca_pt"](r) for r in A.values()]), tv([THUOC["Vmat_tat_ca_pt"](r) for r in Bb.values()])
    mA, mB = sum((A[s].get("n_mat_kiem") or 0) for s in ss), sum((Bb[s].get("n_mat_kiem") or 0) for s in ss)
    ok = None if None in (hA, hB, vA, vB) else bool(hA < hB and vA < vB and mA <= mB)
    out["cap"][f"{vg}|{a}-{b}"] = {"seed": ss, "seed_chi_A": sorted(set(A) - set(Bb)), "seed_chi_B": sorted(set(Bb) - set(A)), "thuoc": res, "THANG_VUNG": ok, "mat_tong": [mA, mB]}
a_, b_ = B.get(("dong", "D4A2"), {}).get(6511), B.get(("U5b", "D4A2"), {}).get(6511)
if a_ and b_: out["dong_nhat"] = {m: [THUOC[m](a_), THUOC[m](b_)] for m in CHINH}
json.dump(out, open(U + "/u5c/ket_qua_u5c.json", "w"), indent=1, ensure_ascii=False)
f4 = lambda x: "–" if x is None else f"{x:.4g}"
L = ["# Lượt u5_c — bảng tự động kiểu CG (CHƯA rà soát chéo)", "", "trái/phải = nhánh A/B; thắng = A nhỏ hơn (tốt hơn) ở seed đó; tv = trung vị theo seed; Δ% = trung vị hiệu tương đối (A − B)/B.", "",
     f"Đồng nhất v122u5f (cờ tắt) ≡ v122u5e, D4A2 6511: {out['dong_nhat']}", ""]
for k, v in out["cap"].items():
    L += [f"## {k} — n = {len(v['seed'])}; THẮNG_VÙNG = {v['THANG_VUNG']}; Σ mất {v['mat_tong'][0]}/{v['mat_tong'][1]}; seed chỉ A {v['seed_chi_A']}, chỉ B {v['seed_chi_B']}",
          "| thước | tv A | tv B | thắng/thua/hoà | Δ% | p |", "|---|---|---|---|---|---|"]
    for m, x in v["thuoc"].items():
        L.append(f"| {m} | {f4(x['tv_A'])} | {f4(x['tv_B'])} | {x['thang']}/{x['thua']}/{x['hoa']} | {100 * x['tv_hieu_tuong_doi']:.1f} | {f4(x['p_wilcoxon'])} |")
    L.append("")
open(U + "/u5c/ket_qua_u5c.md", "w").write("\n".join(L)); print("\n".join(L[:40]))
