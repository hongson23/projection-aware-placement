# -*- coding: utf-8 -*-
"""phan_tich_n14d5.py (22/11/2026) — bảng N14 MẬT ĐỘ CAO (W3, vật chứa thật Allen hiPSC), lối đồ hoạ máy tính, cùng thước bảng thống nhất N1.
Nguồn mỗi việc u5_n14d5_<v>_<CellId>_<nhánh>_<seed> (v ∈ d5, k8; nhánh U5b | td (PLACE) | nn1x (thử–loại 1×) | unc8 (SA c08 1×)):
  u5v2/doc_lap/<id>.json (phan_tich_doc_lap_u5.py --lo <lô chấm>, không sửa: HD trước/sau của bộ dựng, %V mất, G_c^K, bộ phận mất), kiem3.csv của
  bộ kiểm v3 (tách luật sự thật theo TỪNG CA thành ba lớp như phan_tich_n14.py bản 2), thong_tin.json (giờ, khả thi PROJECT, φ_eff đo).
Thước CHÍNH (nhỏ hơn = tốt hơn; ER bị màng thật cắt ~81 % ở mọi nhánh — tính chất dữ liệu — nên thước chính KHÔNG gồm ER):
  HD_max_dat, HD_p95_dat (Hausdorff hai chiều mesh trước/sau chiếu, mọi bộ phận đặt; nhân cố định = 0), Vmat_dat_pt (% thể tích bộ phận đặt mất),
  n_mat (bộ phận mất mesh, bộ kiểm), n_bep (bộ phận đặt còn < 50 %), Gc_K, luat_dat (ca luật sự thật lớp phép đặt, trừ luật mềm 4.5b/4.6b/4.8b,
  trừ 4.3b trên ống dưới lưới, + phạt LUAT_MAT mỗi bộ phận mất như V_test).
Thước PHỤ (cùng cột bảng N1): HD_max_ca_ER, HD_ER, Vmat_ca_ER_pt, Vmat_ER_pt, luat_du_lieu, luat_phu_thuoc, kha_thi (PROJECT hội tụ), con_vp_max_um,
  φ_eff đo, danh_gia_phat, giay_dat, giay_tong (máy dùng chung ⇒ chỉ mô tả; số giờ đưa vào bài lấy từ lượt đo độc quyền N3′).
Báo: trung vị [Q1; Q3] trên 15 cặp (tế bào, seed) mỗi biến thể × nhánh; trung vị theo seed từng tế bào; đếm A tốt/hoà/kém theo cặp và theo tế bào
  cho U5b − {td, nn1x, unc8} và td − {nn1x, unc8}; hoà ⇔ |hiệu| ≤ TOL. Lượt hỏng (loi/ hoặc thiếu số): báo theo nhánh, trong cặp = THUA cho nhánh hỏng
  (luật điều phối 22/11), trung vị trên lượt có số kèm n. Không kiểm định. CHƯA rà soát chéo khi in.
Dùng: python phan_tich_n14d5.py [--seed 6541] ⇒ n14d5/ket_qua/N14D5[_s<seed>].json + .md"""
import os, sys, json, csv, glob, math, re, collections, hashlib, argparse
import numpy as np
R = os.path.expanduser("~/Vision_team_data_189/sonhh20/ch3_tebao_dep"); B = R + "/bo_dung_chung"; T = B + "/tang_dg"; D = B + "/n14d5"
DL = T + "/u5v2/doc_lap"; CR = R + "/bo_kiem_v3/cham_ra"; sys.path.insert(0, D)
CELL = ["798836", "989830", "593409", "707681", "997092"]; NH = ["U5b", "td", "nn1x", "unc8"]; SEED = [6541, 6542, 6543]
TEN = {"U5b": "U5b", "td": "PLACE", "nn1x": "thử–loại 1×", "unc8": "SA c08 1×"}; BT = {"d5": "D5-like (φ_eff ≈ 0,45)", "k8": "k08-like (khe 0,8, φ_eff ≈ 0,51)"}
MEM = {"4.5b", "4.6b", "4.8b"}; LUAT_MAT = {"ti_the": ["4.3a", "4.3b", "M.a", "M.c₀", "M.h"], "lysosome": ["4.6a", "M.a", "M.c₀", "M.h"]}; LUAT_MAT_MD = ["M.a", "M.c₀", "M.h"]
CHO_SAN = {"mang", "nhan"}; PHU = {"luoi_noi_chat"}
CHINH = ["HD_max_dat", "HD_p95_dat", "Vmat_dat_pt", "n_mat", "n_bep", "Gc_K", "luat_dat"]
PHU_T = ["HD_max_ca_ER", "HD_ER", "Vmat_ca_ER_pt", "Vmat_ER_pt", "ER_ngoai_mang_pt", "Vmat_ER_trong_mang_pt", "n_khong_dat", "luat_du_lieu", "luat_phu_thuoc", "kha_thi", "con_vp_max_um", "phi_eff", "danh_gia_phat", "giay_dat", "giay_tong"]
TOL = {"HD_max_dat": 0.01, "HD_p95_dat": 0.01, "Vmat_dat_pt": 0.1, "n_mat": 0.5, "n_bep": 0.5, "Gc_K": 0.01, "luat_dat": 0.5}
CAP = [("U5b", "td"), ("U5b", "nn1x"), ("U5b", "unc8"), ("td", "nn1x"), ("td", "unc8")]
def jl(p):
    try: return json.load(open(p))
    except (OSError, ValueError): return None
def tl(b): return b.rstrip("0123456789").rstrip("_")
def _er_ngoai(t, mode):
    try:
        a, b = t["the_tich_truoc_sau"]["luoi_noi_chat"]; ng = t["n14"]["truoc_chieu"]["ngoai_cha_phu_thuoc_um3"]
        return 100 * ng / a if mode == 0 else 100 * max(a - ng - b, 0) / max(a - ng, 1e-12)
    except (KeyError, TypeError, ValueError): return None
def mot(i):
    r = jl(f"{DL}/{i}.json"); t = jl(f"{T}/kq/{i}/thong_tin.json")
    if r is None or t is None: return None
    f = glob.glob(f"{CR}/vstdg*_mu_sau__{r['mu']}/kiem3.csv"); rs = list(csv.DictReader(open(f[0], encoding="utf-8"))) if len(f) == 1 else []
    lo = dict(t.get("loai") or {}); lo.setdefault("luoi_noi_chat", "luoi_noi_chat"); L = lambda k: lo.get(k, tl(k))
    dl = set(t.get("duoi_luoi") or []); lop = collections.Counter()
    for x in rs:
        if x["phan_quyet"] != "su_that" or x["loai_luat"] != "su_that_v1.2" or x["luat"] in MEM or not int(x["so_ca_truot"] or 0): continue
        cas = [c for c in (x["bo_phan_truot"] or "").split(";") if c.strip()]
        if len(cas) != int(x["so_ca_truot"]): cas = [";".join(cas)] * int(x["so_ca_truot"])
        for c in cas:
            ten = re.findall(r"[A-Za-z_][A-Za-z_0-9]*", c); bp = {tl(b) for b in ten}
            if x["luat"] == "4.3b" and any(b in dl for b in ten): lop["duoi_luoi_43b"] += 1; continue
            lop["du_lieu" if bp and bp <= CHO_SAN else ("phu_thuoc" if bp and bp <= CHO_SAN | PHU else "dat")] += 1
    mat = r.get("mat_kiem") or []
    hl = r.get("hd_ts_theo_loai") or {}; th = r.get("theo_loai") or {}
    vd = [(v["V_truoc"], v["V_mat"]) for l, v in th.items() if l not in ("luoi_noi_chat", "nhan")]
    er = th.get("luoi_noi_chat") or {}
    tu = r.get("hd_ts_tu_the") or {}
    x = {"id": i, "HD_max_dat": tu.get("max"), "HD_p95_dat": tu.get("p95"),
         "Vmat_dat_pt": 100 * sum(b for a, b in vd) / max(sum(a for a, b in vd), 1e-12) if vd else None,
         "n_mat": r.get("n_mat_kiem"), "n_bep": len([k for k in (r.get("con_lai_duoi_50pt") or []) if L(k) not in ("luoi_noi_chat",)]), "Gc_K": r.get("Gc_K"),
         "luat_dat": lop["dat"] + sum(len(LUAT_MAT.get(L(k), LUAT_MAT_MD)) for k in mat) if rs else None,
         "HD_max_ca_ER": (r.get("hd_ts_tat_ca") or {}).get("max"), "HD_ER": (hl.get("luoi_noi_chat") or {}).get("max"),
         "Vmat_ca_ER_pt": 100 * r["V_mat_ti_le_tat_ca"] if r.get("V_mat_ti_le_tat_ca") is not None else None,
         "Vmat_ER_pt": 100 * er["V_mat_ti_le"] if er.get("V_mat_ti_le") is not None else None,
         "ER_ngoai_mang_pt": _er_ngoai(t, 0), "Vmat_ER_trong_mang_pt": _er_ngoai(t, 1), "n_khong_dat": (t.get("bao_cao_bo_tri") or {}).get("khong_dat"),
         "luat_du_lieu": lop["du_lieu"] if rs else None, "luat_phu_thuoc": lop["phu_thuoc"] if rs else None, "n_43b_duoi_luoi": lop["duoi_luoi_43b"],
         "kha_thi": None if t.get("kha_thi") is None else (0 if t["kha_thi"] else 1), "con_vp_max_um": t.get("con_vp_max_um"),
         "danh_gia_phat": (t.get("bao_cao_bo_tri") or {}).get("danh_gia_phat"), "giay_dat": (t.get("bao_cao_bo_tri") or {}).get("giay_dat"), "giay_tong": t.get("giay_tong"),
         "hd_loai": {l: (v or {}).get("max") for l, v in hl.items()}}
    try:
        import kha_thi_n14d5 as KT; x["phi_eff"] = KT.mot(i)["phi_eff"]
    except Exception: x["phi_eff"] = None
    return x
def tk(v):
    v = [float(a) for a in v if a is not None and a == a]
    if not v: return None
    q = np.quantile(v, [0.25, 0.5, 0.75]); return {"tv": float(q[1]), "q1": float(q[0]), "q3": float(q[2]), "n": len(v)}
def f_(s):
    if s is None: return "–"
    g = lambda a: f"{a:.3g}" if abs(a) < 1000 else f"{a:.0f}"
    return f"{g(s['tv'])} [{g(s['q1'])}; {g(s['q3'])}]" + ("" if s["n"] == 15 else f" (n={s['n']})")
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--seed", type=int, default=None); a = ap.parse_args()
    sd = [a.seed] if a.seed else SEED; KQ = {"seed": sd, "bien_the": {}}
    md = ["# N14 mật độ cao — bảng tự động (CHƯA rà soát chéo)\n", f"Seed {sd}; 5 tế bào Allen hiPSC (màng + nhân thật cố định); bộ dựng v122u5n14; bộ kiểm v3; trung vị [Q1; Q3] trên cặp (tế bào, seed).\n"]
    for v in ("d5", "k8"):
        X = {}; hong = collections.defaultdict(list)
        for c in CELL:
            for nh in NH:
                for s in sd:
                    i = f"u5_n14d5_{v}_{c}_{nh}_{s}"
                    if not os.path.exists(f"{T}/xong/{i}") and not os.path.exists(f"{T}/loi/{i}"): continue
                    x = mot(i) if os.path.exists(f"{T}/xong/{i}") else None
                    if x is None: hong[nh].append(f"{c}/{s}")
                    X[(c, nh, s)] = x
        if not X: continue
        Kv = {"bang": {}, "theo_te_bao": {}, "cap": {}, "hong": dict(hong)}
        md.append(f"\n## {BT[v]}\n\n| nhánh | " + " | ".join(CHINH + PHU_T) + " |\n|" + "---|" * (1 + len(CHINH) + len(PHU_T)) + "\n")
        for nh in NH:
            row = {m: tk([X[k][m] for k in X if k[1] == nh and X[k]]) for m in CHINH + PHU_T}; Kv["bang"][nh] = row
            md.append(f"| {TEN[nh]} | " + " | ".join(f_(row[m]) for m in CHINH + PHU_T) + " |\n")
            Kv["theo_te_bao"][nh] = {c: {m: tk([X[k][m] for k in X if k[0] == c and k[1] == nh and X[k]]) for m in CHINH} for c in CELL}
        md.append("\nLượt hỏng: " + (json.dumps(dict(hong), ensure_ascii=False) if hong else "0") + "\n")
        md.append("\n### Đếm cặp (A tốt / hoà / kém; theo cặp (tế bào, seed) | theo tế bào, trung vị seed)\n\n| A − B | " + " | ".join(CHINH) + " |\n|" + "---|" * (1 + len(CHINH)) + "\n")
        for A, Bn in CAP:
            cells = []
            for m in CHINH:
                w = t_ = l_ = 0; wc = tc = lc = 0
                for c in CELL:
                    da, db = [], []
                    for s in sd:
                        ka, kb = (c, A, s), (c, Bn, s)
                        if ka not in X or kb not in X: continue
                        xa, xb = (X[ka] or {}).get(m), (X[kb] or {}).get(m)
                        if xa is None and xb is None: continue
                        da.append(math.inf if xa is None else xa); db.append(math.inf if xb is None else xb)
                        if xa is None: l_ += 1; continue
                        if xb is None: w += 1; continue
                        d = xa - xb
                        if abs(d) <= TOL[m]: t_ += 1
                        elif d < 0: w += 1
                        else: l_ += 1
                    if da:
                        ma, mb = float(np.median(da)), float(np.median(db))
                        d = 0.0 if ma == mb else (ma - mb if math.isfinite(ma) and math.isfinite(mb) else (-1e9 if ma < mb else 1e9))
                        if abs(d) <= TOL[m]: tc += 1
                        elif d < 0: wc += 1
                        else: lc += 1
                Kv["cap"].setdefault(f"{A}-{Bn}", {})[m] = {"cap": [w, t_, l_], "te_bao": [wc, tc, lc]}; cells.append(f"{w}/{t_}/{l_} \\| {wc}/{tc}/{lc}")
            md.append(f"| {TEN[A]} − {TEN[Bn]} | " + " | ".join(cells) + " |\n")
        md.append("\n### Hausdorff max theo loại (trung vị, µm)\n\n| nhánh | " + " | ".join(["golgi", "trung_the", "nhan_con", "ti_the", "lysosome", "luoi_noi_chat"]) + " |\n|---|---|---|---|---|---|---|\n")
        for nh in NH:
            md.append(f"| {TEN[nh]} | " + " | ".join(f_(tk([(X[k] or {}).get('hd_loai', {}).get(l) for k in X if k[1] == nh and X[k]])) for l in ["golgi", "trung_the", "nhan_con", "ti_the", "lysosome", "luoi_noi_chat"]) + " |\n")
        Kv["viec"] = {"_".join(map(str, k)): x for k, x in X.items()}; KQ["bien_the"][v] = Kv
    md.append(f"\nmd5 phan_tich_n14d5.py {hashlib.md5(open(__file__, 'rb').read()).hexdigest()}\n")
    os.makedirs(f"{D}/ket_qua", exist_ok=True); tg = "" if not a.seed else f"_s{a.seed}"
    json.dump(KQ, open(f"{D}/ket_qua/N14D5{tg}.json", "w"), ensure_ascii=False, indent=1, default=float); open(f"{D}/ket_qua/N14D5{tg}.md", "w").write("".join(md))
    print("".join(md)[:3000])
if __name__ == "__main__":
    main()
