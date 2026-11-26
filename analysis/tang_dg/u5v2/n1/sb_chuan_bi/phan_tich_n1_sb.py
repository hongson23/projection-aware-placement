# -*- coding: utf-8 -*-
"""(BẢN SB 22/11: = phan_tich_n1.py 69037361, chỉ đổi DL → u5v2/doc_lap_sb và tệp ra → n1/sb_chuan_bi/*_sb: bảng N1 CÓ SỬA MESH ĐỒNG NHẤT; bảng không sửa = độ nhạy.) phan_tich_n1.py (22/11/2026) — BẢNG THỐNG NHẤT N1 (PI 22/11 06:50 "ok, bảng thống nhất"): U5b (u5_c_U5b_*), PLACE td / thử–loại nn1x / SA un1x
(e3b_B_*), cellPACK (<ô>_cpchinh_*), Packmol (<ô>_pkchinh_*), ô D1/D4/D4A2/D5 × seed 6201–6210, bộ dựng v122e3 --eval_dung_tu_dau, bộ kiểm v3.
Mọi số lấy từ u5v2/doc_lap/<id>.json (phan_tich_doc_lap_u5.py + hd_e3b.py, cùng thước cho mọi nhánh) và kiem3.csv / thong_tin.json gốc.
Bảng kiểu C&G: trung vị [Q1; Q3] theo seed; cặp theo seed (thắng/thua/hoà, trung vị hiệu, p Wilcoxon hai phía — chỉ mô tả).
Thêm (không đổi thước gốc): %V mất trừ ER; bộ phận KHÔNG ĐẶT (loại (a), chỉ đối chứng ngoài có) và G_c^K / bộ phận mất / luật trượt khi bỏ loại (a);
G_c^K của cellPACK trên M′ (bỏ danh sách M′ đóng băng, quyết định 3) — mô tả; phát biểu H5 chỉ theo phan_tich_cp.py chinh (KQ_CP.md).
Dùng: python phan_tich_n1.py [--cuoi]  ⇒ u5v2/n1/ket_qua_n1.json + BangThongNhat_N1_20261122.md (TẠM khi thiếu nhánh). CHƯA rà soát chéo."""
import os, sys, json, glob, csv, hashlib, collections, datetime
import numpy as np
R = os.path.expanduser("~/Vision_team_data_189/sonhh20/ch3_tebao_dep"); O = R + "/bo_dung_chung/tang_dg"; U = O + "/u5v2"; N1 = U + "/n1"; DL = U + "/doc_lap_sb"
CPV = os.path.expanduser("~/Vision_team_data_189/sonhh20/cellpack_ngoai/v122e3"); PKW = R + "/ngoai_packmol/n1"; CHAM_RA = R + "/bo_kiem_v3/cham_ra"
H = 0.1797; SEED = list(range(6201, 6211)); O_ = ["D1", "D4", "D4A2", "D5"]
try:
    from scipy.stats import wilcoxon
except Exception: wilcoxon = None
sys.path.insert(0, R + "/bo_kiem_v3/cham_chinh_dv/phan_tich"); import phan_tich_chinh_dv as PT
NHANH = {"U5b": ("U5b", lambda o, s: f"u5_c_U5b_{o}_{s}"), "P4": ("PLACE", lambda o, s: f"e3b_B_{o}_td_{s}"), "TL": ("thử–loại", lambda o, s: f"e3b_B_{o}_nn1x_{s}"),
         "SA": ("SA", lambda o, s: f"e3b_B_{o}_un1x_{s}"), "CP": ("cellPACK", lambda o, s: f"{o}_cpchinh_{s}"), "PK": ("Packmol", lambda o, s: f"{o}_pkchinh_{s}")}
LOAI = ["luoi_noi_chat", "nhan", "golgi", "trung_the", "nhan_con", "ti_the", "lysosome"]
TEN_LOAI = {"luoi_noi_chat": "ER", "nhan": "nhân", "golgi": "Golgi", "trung_the": "trung thể", "nhan_con": "nhân con", "ti_the": "ti thể", "lysosome": "lysosome"}
def g(r, *k):
    for x in k:
        r = r.get(x) if isinstance(r, dict) else None
        if r is None: return None
    return r
def jl(p):
    try: return json.load(open(p))
    except (OSError, ValueError): return None
def md5(p):
    try: return hashlib.md5(open(p, "rb").read()).hexdigest()
    except OSError: return None
_K3 = {}
def kiem3(mu):
    if mu not in _K3:
        f = glob.glob(f"{CHAM_RA}/*_mu_sau__{mu}/kiem3.csv"); _K3[mu] = list(csv.DictReader(open(f[0], encoding="utf-8"))) if len(f) == 1 else None
    return _K3[mu]
LUAT_MAT = {"ti_the": ["4.3a", "4.3b", "M.a", "M.c₀", "M.h"], "lysosome": ["4.6a", "M.a", "M.c₀", "M.h"]}; LUAT_MAT_MD = ["M.a", "M.c₀", "M.h"]
def m_phay(o, s, tt):
    k = (tt.get("cau_noi") or {}).get("k"); j = jl(f"{CPV}/cham/M_phay/{o}_{s}_k{k}_chinh.json")
    return set(j["M_phay"]["bo_khoi_M_phay"]) if j and "M_phay" in j else None
def mot(nh, o, s):
    id_ = NHANH[nh][1](o, s); r = jl(f"{DL}/{id_}.json")
    if r is None: return None
    tt = jl(f"{O}/kq/{id_}/thong_tin.json") or {}; tu = tt.get("tu_the") or {}; lo = dict(tt.get("loai") or {}); lo.setdefault("luoi_noi_chat", "luoi_noi_chat")
    loai_cua = lambda k: lo.get(k, k.rstrip("0123456789").rstrip("_"))
    kd = {k for k, v in tu.items() if v.get("khong_dat") or v.get("t") is None}
    x = {"id": id_, "HD_max": g(r, "hd_ts_tat_ca", "max"), "HD_p95": g(r, "hd_ts_tat_ca", "p95"),
         "Vmat_pt": 100 * r["V_mat_ti_le_tat_ca"] if r.get("V_mat_ti_le_tat_ca") is not None else None,
         "n_mat": r.get("n_mat_kiem"), "n_bep": len(r.get("con_lai_duoi_50pt") or []), "n_vo_manh": len(r.get("vo_manh") or []), "n_khong_dat": len(kd),
         "Gc_K": r.get("Gc_K"), "V_test_tru_43b": r.get("V_test_tru_43b_duoi_luoi"), "giay_dat": r.get("giay_dat"), "giay_tong": r.get("giay_tong")}
    for l in LOAI:
        for q in ("max", "p95"): x[f"HD_{l}_{q}"] = g(r, "hd_ts_theo_loai", l, q)
    tl = r.get("theo_loai") or {}; vt = sum(v["V_truoc"] for l, v in tl.items() if l != "luoi_noi_chat"); vm = sum(v["V_mat"] for l, v in tl.items() if l != "luoi_noi_chat")
    x["Vmat_tru_ER_pt"] = 100 * vm / vt if vt > 0 else None
    rs = kiem3(r["mu"]) if r.get("mu") else None
    if rs is not None:
        def rbs(k):
            try: return PT.ban_kinh_bao(tu[k]["khuon"], tu[k]["s"])
            except (ValueError, KeyError): return float("nan")
        Di = {rr["phan_quyet"]: rr for rr in rs if rr["luat"] == "D_i"}; gc = float([rr for rr in rs if rr["luat"] == "G_c"][0]["trang_thai"])
        mat = [k for k, d in Di.items() if d["trang_thai"] == "KHONG_CO_MESH"]; mat_a = [k for k in mat if k in kd]; mat_b = [k for k in mat if k not in kd]
        x["Gc_K_tru_khong_dat"] = gc + sum(max(rbs(k) - H, 0) if rbs(k) == rbs(k) else float("nan") for k in mat_b)
        x["n_mat_tru_khong_dat"] = len(mat_b)
        x["V_test_tru_43b_tru_khong_dat"] = None if x["V_test_tru_43b"] is None else x["V_test_tru_43b"] - sum(len(LUAT_MAT.get(loai_cua(k), LUAT_MAT_MD)) for k in mat_a)
        if nh == "CP":
            mp = m_phay(o, s, tt); x["n_M_phay_bo"] = None if mp is None else len(mp)
            if mp is not None:
                bo = 0.0
                for k in mp - kd:
                    d = Di.get(k)
                    if d is None: continue
                    if d["trang_thai"] == "KHONG_CO_MESH": bo += max(rbs(k) - H, 0) if rbs(k) == rbs(k) else float("nan")
                    elif d["trang_thai"] not in ("KHONG_DUNG_LAI_DUOC", ""): bo += max(float(d["trang_thai"]) - H, 0)
                x["Gc_K_M_phay"] = x["Gc_K_tru_khong_dat"] - bo
    if nh == "PK": x["khop_tu_the"] = (jl(f"{PKW}/kq/{id_}/KHOP_TU_THE.json") or {}).get("dat")
    return x
def _h2(d, i):
    try: return int(open(f"{d}/hong/{i}").read().strip() or 0) >= 2
    except (OSError, ValueError): return False
def trang_thai(nh, o, s, cuoi, co):
    """Luật điều phối 22/11 (ghi TRƯỚC khi có bảng): ok = có số; hong = lượt hỏng (không bao giờ bỏ im lặng); cho = chưa xong (chỉ ở bản tạm)."""
    if co: return "ok"
    id_ = NHANH[nh][1](o, s)
    if nh in ("U5b", "P4", "TL", "SA"):
        if os.path.exists(f"{O}/loi/{id_}"): return "hong"
        if os.path.exists(f"{O}/xong/{id_}") and os.path.exists(U + "/u5c/XONG_U5C"): return "hong"
    elif nh == "CP":
        H_ = CPV + "/hang"
        if _h2(H_, f"N_{o}_{s}_chinh") or _h2(H_, f"C_{o}_{s}") or os.path.exists(f"{CPV}/thu/{o}_{s}_C_HONG") or os.path.exists(f"{CPV}/thu/{o}_{s}_CHON_HONG"): return "hong"
        if any(os.path.exists(f"{N1}/{m}_cpn1") for m in ("XONG_CHAM", "LOI_CHAM", "DUNG_TU_THE")): return "hong"
    elif nh == "PK":
        lo = "pkn1b" if o == "D5" else "pkn1a"
        if _h2(PKW, f"N_{o}_{s}") or _h2(PKW, f"C_{o}_{s}") or os.path.exists(f"{PKW}/thu/{o}_{s}_CHON_HONG"): return "hong"
        if any(os.path.exists(f"{N1}/{m}_{lo}") for m in ("XONG_CHAM", "LOI_CHAM", "DUNG_TU_THE")): return "hong"
    return "hong" if cuoi else "cho"
def cau_hinh(nh, o):
    """số cấu hình (seed × k) hỏng và (Packmol) không hội tụ, ở mọi k và ở k được chọn — báo, không loại."""
    r = {"cau_hinh_hong": 0, "khong_hoi_tu_moi_k": 0, "khong_hoi_tu_k_chon": 0, "so_cau_hinh": 0}
    if nh == "CP":
        for s in SEED:
            for k in range(4):
                d = f"{CPV}/thu/{o}_{s}_k{k}"; r["so_cau_hinh"] += 1
                if os.path.exists(d + "/XUAT_HONG") or _h2(CPV + "/hang", f"X_{o}_{s}_{k}") or _h2(CPV + "/hang", f"V_{o}_{s}_{k}_chinh"): r["cau_hinh_hong"] += 1
    if nh == "PK":
        for s in SEED:
            ch = jl(f"{PKW}/thu/{o}_{s}_chon.json") or {}
            for k in range(4):
                d = f"{PKW}/thu/{o}_{s}_k{k}"; r["so_cau_hinh"] += 1
                if os.path.exists(d + "/XUAT_HONG") or _h2(PKW, f"K_{o}_{s}_{k}"): r["cau_hinh_hong"] += 1
                j = jl(d + "/tu_the_pk.json")
                if j and not j["packmol"].get("hoi_tu"):
                    r["khong_hoi_tu_moi_k"] += 1
                    if ch.get("k_chon") == k: r["khong_hoi_tu_k_chon"] += 1
    return r
def tk(v):
    v = [float(a) for a in v if a is not None and a == a]
    if not v: return None
    q = np.quantile(v, [0.25, 0.5, 0.75]); return {"tv": float(q[1]), "q1": float(q[0]), "q3": float(q[2]), "n": len(v)}
CHINH = ["HD_max", "HD_p95", "HD_luoi_noi_chat_max", "Vmat_pt", "Vmat_tru_ER_pt", "n_mat", "n_bep", "n_khong_dat", "Gc_K", "Gc_K_tru_khong_dat", "V_test_tru_43b", "giay_dat", "giay_tong"]
TEN = {"HD_max": "HD max (mọi bộ phận, kể ER)", "HD_p95": "HD p95", "HD_luoi_noi_chat_max": "HD max ER", "Vmat_pt": "%V mất (kể ER)", "Vmat_tru_ER_pt": "%V mất (trừ ER)",
       "n_mat": "bộ phận mất mesh", "n_bep": "bộ phận bẹp (< 50 % thể tích)", "n_khong_dat": "không đặt (loại a)", "Gc_K": "G_c^K", "Gc_K_tru_khong_dat": "G_c^K bỏ loại (a)",
       "V_test_tru_43b": "luật sự thật trượt (trừ 4.3b dưới lưới)", "giay_dat": "giây đặt", "giay_tong": "giây tổng (bộ dựng/nhập)"}
CAP_A = ["CP", "PK", "U5b"]; CAP_B = ["U5b", "P4", "TL", "SA"]; THUOC_CAP = ["HD_max", "HD_luoi_noi_chat_max", "Vmat_pt", "Vmat_tru_ER_pt", "n_mat", "Gc_K", "Gc_K_tru_khong_dat", "V_test_tru_43b"]
def chinh(cuoi):
    B = {(nh, o): {s: mot(nh, o, s) for s in SEED} for nh in NHANH for o in O_}
    B = {k: {s: v for s, v in d.items() if v is not None} for k, d in B.items()}
    TT = {(nh, o): {s: trang_thai(nh, o, s, cuoi, s in B[(nh, o)]) for s in SEED} for nh in NHANH for o in O_ if not (nh == "PK" and o == "D1")}
    out = {"tao_luc": datetime.datetime.now().isoformat(timespec="seconds"), "cuoi": cuoi, "so_te_bao": {f"{nh}|{o}": len(B[(nh, o)]) for nh, o in B},
           "seed_thieu": {f"{nh}|{o}": [s for s in SEED if s not in B[(nh, o)]] for nh, o in B if not (nh == "PK" and o == "D1")},
           "hong": {f"{nh}|{o}": {"ok": sum(v == "ok" for v in d.values()), "hong": [s for s, v in d.items() if v == "hong"], "cho": [s for s, v in d.items() if v == "cho"], **cau_hinh(nh, o)} for (nh, o), d in TT.items()}, "bang": {}, "cap": {}}
    for (nh, o), d in B.items():
        if not d: continue
        out["bang"][f"{nh}|{o}"] = {m: tk([x.get(m) for x in d.values()]) for m in CHINH + [f"HD_{l}_{q}" for l in LOAI for q in ("max", "p95")] + ["Gc_K_M_phay", "n_M_phay_bo", "n_vo_manh"]}
    for o in O_:
        for a in CAP_A:
            for b in CAP_B:
                if a == b: continue
                A, Bb = B[(a, o)], B[(b, o)]; TA, TB = TT.get((a, o)), TT.get((b, o))
                if TA is None or TB is None: continue
                ss = sorted(s for s in SEED if TA[s] != "cho" and TB[s] != "cho")
                if not ss: continue
                res = {}
                for m in THUOC_CAP:
                    th_h = sum(TA[s] == "ok" and TB[s] == "hong" for s in ss); tu_h = sum(TA[s] == "hong" and TB[s] == "ok" for s in ss); hh = sum(TA[s] == "hong" and TB[s] == "hong" for s in ss)
                    pr = [(A[s].get(m), Bb[s].get(m)) for s in ss if TA[s] == "ok" and TB[s] == "ok"]; pr = [(x, y) for x, y in pr if x is not None and y is not None and x == x and y == y]
                    dd = np.array([x - y for x, y in pr], float); p = None
                    if len(dd) and wilcoxon is not None and np.any(np.abs(dd) > 1e-12):
                        try: p = float(wilcoxon(dd, zero_method="wilcox", alternative="two-sided").pvalue)
                        except Exception: p = None
                    res[m] = {"n": len(pr) + th_h + tu_h + hh, "n_ca_hai_co_so": len(pr), "A_tot_hon": int((dd < -1e-9).sum()) + th_h, "A_kem_hon": int((dd > 1e-9).sum()) + tu_h,
                              "hoa": int((np.abs(dd) <= 1e-9).sum()) + hh, "do_hong": [th_h, tu_h, hh], "tv_hieu": float(np.median(dd)) if len(dd) else None, "p": p}
                out["cap"][f"{o}|{a}-{b}"] = {"seed": ss, "thuoc": res}
    ma = {p: md5(p) for p in [os.path.abspath(__file__), U + "/phan_tich_doc_lap_u5.py", U + "/hd_e3b.py", N1 + "/cham_n1.sh", N1 + "/cho_bang_n1.sh", R + "/bo_dung_chung/ma3/botri_chung_v122e3.py",
                              R + "/bo_dung_chung/ma3/botri_chung_v122u5f.py", CPV + "/cau_noi_cp_e3.py", CPV + "/hang_n1.sh", R + "/ngoai_packmol/cau_noi_pk_e3.py", R + "/ngoai_packmol/n1/hang_pk_n1.sh",
                              os.path.expanduser("~/Vision_team_data_189/sonhh20/cellpack_ngoai/v122/cham/phan_tich_cp.py")]}
    out["md5"] = ma
    json.dump(out, open(N1 + "/sb_chuan_bi/ket_qua_n1_sb.json", "w"), ensure_ascii=False, indent=1)
    f = lambda t, n=3: "–" if not t else f"{t['tv']:.{n}g} [{t['q1']:.{n}g}; {t['q3']:.{n}g}]"
    L = [f"# Bảng thống nhất N1 — {'BẢN CUỐI' if cuoi else 'TẠM (còn thiếu nhánh)'} ({out['tao_luc']})", "",
         'DUYET: PI 22/11/2026 06:50 — "ok, bảng thống nhất".', "",
         "Trạng thái: CHƯA rà soát chéo (số tự động; không dùng để kết luận trước khi rà).", "",
         "Seed 6201–6210; ô D1, D4, D4A2, D5 (Packmol không có D1); bộ dựng v122e3 `--eval_dung_tu_dau` (U5b: v122u5f); bộ kiểm v3; trung vị [Q1; Q3] theo seed.", "",
         "Lệnh: `python u5v2/n1/phan_tich_n1.py" + (" --cuoi`" if cuoi else "`") + " (bộ chờ `u5v2/n1/cho_bang_n1.sh`).", "", "md5:", ""]
    L += [f"- `{os.path.relpath(p, os.path.expanduser('~/Vision_team_data_189/sonhh20'))}` {m}" for p, m in ma.items()]
    L += ["", "Số tế bào có số liệu (nhánh|ô): " + ", ".join(f"{k} {v}" for k, v in out["so_te_bao"].items() if not k.startswith("PK|D1")), ""]
    L += ["## Lượt hỏng / không hội tụ (luật điều phối 22/11: không bỏ im lặng; trong cặp, lượt hỏng = thua cho nhánh hỏng; trung vị trên lượt thành công, kèm n)", "",
          "| nhánh | ô | thành công | hỏng (seed) | chưa xong (seed) | cấu hình hỏng / tổng | Packmol không hội tụ: mọi k / k chọn |", "|---|---|---|---|---|---|---|"]
    for k, v in out["hong"].items():
        nh, o = k.split("|"); L.append(f"| {NHANH[nh][0]} | {o} | {v['ok']} | {v['hong'] or '–'} | {v['cho'] or '–'} | {v['cau_hinh_hong']}/{v['so_cau_hinh']} | {str(v['khong_hoi_tu_moi_k']) + ' / ' + str(v['khong_hoi_tu_k_chon']) if nh == 'PK' else '–'} |")
    L += ["", "Packmol không hội tụ vẫn có pk_out.xyz ⇒ chấm nguyên như kết quả của Packmol (đánh dấu ở trên).", ""]
    for o in O_:
        L += [f"## {o}", "", "| nhánh | n | " + " | ".join(TEN[m] for m in CHINH) + " |", "|" + "---|" * (len(CHINH) + 2)]
        for nh in NHANH:
            t = out["bang"].get(f"{nh}|{o}")
            if not t: continue
            n = max((v["n"] for v in t.values() if v), default=0)
            L.append(f"| {NHANH[nh][0]} | {n} | " + " | ".join(f(t[m]) for m in CHINH) + " |")
        L += ["", f"Hausdorff theo loại ({o}; max / p95, µm):", "", "| nhánh | " + " | ".join(TEN_LOAI[l] for l in LOAI) + " |", "|" + "---|" * (len(LOAI) + 1)]
        for nh in NHANH:
            t = out["bang"].get(f"{nh}|{o}")
            if t: L.append(f"| {NHANH[nh][0]} | " + " | ".join(f"{f(t[f'HD_{l}_max'])} / {f(t[f'HD_{l}_p95'])}" for l in LOAI) + " |")
        t = out["bang"].get(f"CP|{o}")
        if t: L += ["", f"cellPACK {o}: G_c^K trên M′ (bỏ M′ đóng băng, mô tả) {f(t['Gc_K_M_phay'])}; số bộ phận bỏ khỏi M′ {f(t['n_M_phay_bo'])}."]
        L += ["", f"Cặp theo seed ({o}; A tốt hơn = A nhỏ hơn; trung vị hiệu A − B; p Wilcoxon hai phía, mô tả):", "", "| A − B | thước | n (cả hai có số) | A tốt/kém/hoà | trong đó do hỏng (A thắng/A thua/cả hai) | tv hiệu (cả hai có số) | p |", "|---|---|---|---|---|---|---|"]
        for k, v in out["cap"].items():
            if not k.startswith(o + "|"): continue
            a, b = k.split("|")[1].split("-")
            for m, x in v["thuoc"].items():
                ps = "–" if x["p"] is None else f"{x['p']:.3g}"
                L.append(f"| {NHANH[a][0]} − {NHANH[b][0]} | {TEN.get(m, m)} | {x['n']} ({x['n_ca_hai_co_so']}) | {x['A_tot_hon']}/{x['A_kem_hon']}/{x['hoa']} | {'/'.join(map(str, x['do_hong']))} | {'–' if x['tv_hieu'] is None else format(x['tv_hieu'], '.4g')} | {ps} |")
        L.append("")
    pk_sai = [x["id"] for (nh, o), d in B.items() if nh == "PK" for x in d.values() if x.get("khop_tu_the") is not True]
    L += ["## Ghi chú cố định", "",
          f"- Packmol khớp tư thế (tu_the ↔ tu_the_pk.json ≤ 1e-6) KHÔNG đạt / thiếu: {pk_sai if pk_sai else 'không có'}.",
          f"- cellPACK DUNG_RA_LAI (quyết định 1 (B), lượt bộ đếm ≠ 0: không chứng nhận từng bit): {'CÓ' if os.path.exists(CPV + '/hang/DUNG_RA_LAI') else 'không'}.",
          "- Luật lượt hỏng (điều phối 22/11, ghi trước khi có bảng): trung vị trên lượt thành công (n in kèm); trong cặp theo seed lượt hỏng tính THUA cho nhánh hỏng (cả hai hỏng = hoà); bản TẠM bỏ qua seed chưa xong.",
          "- Hausdorff U5b lấy từ bộ dựng (--hd_truoc_sau, trong bộ nhớ); các nhánh khác từ hd_e3b.py trên npz float32 (cùng định nghĩa đỉnh → mặt hai chiều, bỏ màng + tầng 2).",
          "- Chấm cellPACK/Packmol thử lại tối đa 3 lần khi thiếu (cham_n1.sh); U5b/e3b_B một lần (cham_lo_u5.sh). Giờ cellPACK/Packmol = CPU, không so tốc độ.", ""]
    kq = os.path.expanduser("~/Vision_team_data_189/sonhh20/cellpack_ngoai/v122e3/cham/ket_qua")
    L += ["## Tệp kèm", "", f"- cổng trước chiếu cellPACK: `{kq}/KIEM_TRUOC.json` (md5 {md5(kq + '/KIEM_TRUOC.json')})",
          f"- H5 cellPACK trên M và M′ (phan_tich_cp.py chinh): `{kq}/KQ_CP.md` (md5 {md5(kq + '/KQ_CP.md')})", f"- số máy: `{N1}/ket_qua_n1.json`", ""]
    ten = "BangThongNhat_N1_SB_20261122.md"; open(f"{N1}/sb_chuan_bi/{ten}", "w").write("\n".join(L)); print("XONG", ten, "cuoi" if cuoi else "tam")
if __name__ == "__main__":
    chinh("--cuoi" in sys.argv)
