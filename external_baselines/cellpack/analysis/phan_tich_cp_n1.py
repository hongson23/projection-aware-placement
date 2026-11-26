# -*- coding: utf-8 -*-
"""phan_tich_cp_n1.py (24/11/2026, BẢN SAO của phan_tich_cp.py 3d109673; chỉ đổi quy_uoc_quay: 0 = 0 ⇒ không xác định, không trượt) — phan_tich_cp.py (20/11/2026) — phân tích GHI TRƯỚC nhánh cellPACK (phụ lục Bài 1 ch3, DangKyTruoc_BoSung_KiemToan_cellPACK_20261112.md
BẢN 3: §1 Phần B, §3.2, §3.4, §3.6–3.8, §4.B, §5.B, §6, §7, §9.1 N1). Không thêm kiểm định, không đổi ngưỡng.
Nhập (không chép) các hàm của phan_tich_chinh_dv.py a48ce789 (H, TOL0, ban_kinh_bao, wilcoxon_chinh_xac, spearman, LUAT_*).

Dùng:  python phan_tich_cp.py kiem_truoc   → cham/ket_qua/KIEM_TRUOC.json   (§6: D_i^K trên bản TRƯỚC chiếu, 2 seed đầu mỗi ô; chỉ báo,
                                             như tiền lệ v12.1 PI 12/11 + 14/11 — không chặn; KHÔNG tính G_c nào)
       python phan_tich_cp.py chinh        → cham/ket_qua/KQ_CP.json, KQ_CP.md, bang_cp.tex (bảng booktabs kiểu C&G)
Biến: CP_V (gốc v122, chỉ để thử), CP_C (thư mục cham, mặc định cellpack_ngoai/v122/cham), CP_O (danh sách ô, mặc định "D5,D4A2"; thử: "D1"),
      CP_SEED (mặc định 4001-4010), CP_CHAM_RA (cham_ra của bộ kiểm).
Đầu vào: lô chấm cham/lo_<lo>/ (lo = chinh, do_nhay) do cham_cp.sh tạo: khoa/KHOA_MU.json, vs_<ban>/ds_<ban>.txt.ket_qua/tong.csv.
Tế bào nguồn: cellPACK = v122/kq/<ô>_cp<nhãn>_<seed>/ (nhãn chinh | thutu_md | h1 | ja5); PLACE và đối chứng = v121/kq_chinh_dv/<ô>_<nhánh>_<seed>/.

Định nghĩa (đăng ký §3.2, §4.B, §5.B; ĐKC 4.1):
  D_i = D_i^K của bộ kiểm v3 trên bản ĐÃ VỆ SINH; bộ phận mất mesh sau PROJECT (KHONG_CO_MESH, loại (b)) ⇒ D_i = rb·s (bán kính bao).
  M_cp = bộ phận NHỎ (ti_the, lysosome, nhan_con) cellPACK đặt được (tu_the.t ≠ null) ở cấu hình chọn; loại (a) (khong_dat) ngoài M.
  M_pl: cùng k (cau_noi.k = bao_cao_bo_tri.cau_hinh_chon) ⇒ đúng các tên của M_cp; khác k ⇒ với mỗi loại τ bỏ khỏi PLACE đúng số
        bộ phận cellPACK thiếu ở τ, chọn các bộ phận rb·s LỚN NHẤT (bản 2 §3.2); so trên phần còn lại.
  G_c(M) = Σ_{i∈M} max(D_i − h, 0), h = 0,1797 µm.   d = G_c(cellPACK) − G_c(PLACE td)  (H1: d > 0, PLACE ít méo hơn).
  PLACE THUA (CHON_HONG / C_HONG / nhập hỏng / không chấm được ở nhánh nào): d = −(max|d| + 1) (hạng lớn nhất, dấu âm).
  H5: Wilcoxon chính xác một phía (liệt kê 2^{n′}, bỏ |d| < 1e-6), chuỗi D5 → D4A2, α = 0,025 mỗi bước; ô bất đồng (ĐKC 4.1 trên
      20 tế bào của ô: Spearman(G_c^K, G_c^B) < 0,8 hoặc dấu d^K ≠ d^B ở > 2/10 seed; G_c^B từ doi_mat_um trên cùng M) ⇒ chỉ mô tả, bỏ
      qua không chặn; ô không đạt ⇒ ô sau không kiểm. Báo thêm α = 0,05/3. Câu chữ: đạt ∧ trung vị d ≥ h ⇒ "PLACE ít méo hơn" (còn
      cần độ nhạy (ii) cùng chiều); đạt ∧ trung vị d < h ⇒ "có ý nghĩa thống kê nhưng dưới một bước lưới".
  H5a: tỉ lệ đặt được theo loại và ô, gộp 10 seed, KTC 95 % Clopper–Pearson (từ tu_the_cp.json của cấu hình chọn).
  H6: V_test (ĐKC 4.2 trừ 4.5b, 4.6b, 4.8b) — luật đếm trên cả tế bào; phạt bộ phận mất mesh chỉ cho bộ phận trong M (loại (b));
      Wilcoxon p chỉ là số mô tả.
"""
import os, sys, json, csv, math, re, ast, hashlib, collections, statistics, glob
import numpy as np
GOC = os.path.expanduser("~/Vision_team_data_189/sonhh20/ch3_tebao_dep")
V = os.environ.get("CP_V", os.path.expanduser("~/Vision_team_data_189/sonhh20/cellpack_ngoai/v122")); H_ = os.path.join(V, "hang")
C = os.environ.get("CP_C", os.path.join(V, "cham")); RA = os.path.join(C, "ket_qua")
KQP = os.environ.get("CP_KQP", os.path.join(GOC, "bo_dung_chung", "v121", "kq_chinh_dv"))
PL_MAU = os.environ.get("CP_PL_MAU", "{o}_{nh}_{s}")
def ten_pl(o, nh, s): return PL_MAU.format(o=o, nh=nh, s=s)
M_PHAY = os.path.join(C, "M_phay"); DEN = 0.01
CHAM_RA = os.environ.get("CP_CHAM_RA", os.path.join(GOC, "bo_kiem_v3", "cham_ra"))
CC_CHINH = os.path.join(GOC, "bo_kiem_v3", "cham_chinh_dv")
PTF = os.path.join(CC_CHINH, "phan_tich", "phan_tich_chinh_dv.py")
def md5(p): return hashlib.md5(open(p, "rb").read()).hexdigest()
MD5_DO_LAN = os.environ.get("CP_MD5_DO_LAN", "")
assert md5(PTF) == "a48ce7894bef9ca860bf574312ac2206", "phan_tich_chinh_dv.py khác bản đóng băng"
sys.path.insert(0, os.path.dirname(PTF)); import phan_tich_chinh_dv as PT
H = PT.H; TOL = PT.TOL0
O_CP = os.environ.get("CP_O", "D5,D4A2").split(","); a_, b_ = os.environ.get("CP_SEED", "4001-4010").split("-"); SEED = list(range(int(a_), int(b_) + 1))
NHO = ("ti_the", "lysosome", "nhan_con"); NHAN_DN = ("thutu_md", "h1", "ja5"); DOI = ("nn1x", "un1x")
BO_DEM = ("va_cham_cpu_du_phong", "chua_mo_ho", "chua_broken", "v3_va_cham_du_phong", "v3_chua_mo_ho")
ALPHA = 0.025; ALPHA3 = 0.05 / 3
def loai(ten): return ten.rstrip("0123456789").rstrip("_")
def tv(x):
    x = [float(v) for v in x if v is not None and not (isinstance(v, float) and math.isnan(v))]
    if not x: return None
    q = np.quantile(x, [0.25, 0.5, 0.75])
    return {"tv": float(q[1]), "q1": float(q[0]), "q3": float(q[2]), "min": float(min(x)), "max": float(max(x)), "n": len(x)}
def cp_ci(k, n, a=0.05):
    from scipy.stats import beta
    if n == 0: return [float("nan"), float("nan")]
    lo = 0.0 if k == 0 else float(beta.ppf(a / 2, k, n - k + 1)); hi = 1.0 if k == n else float(beta.ppf(1 - a / 2, k + 1, n - k))
    return [lo, hi]
def jl(p):
    try: return json.load(open(p))
    except (OSError, ValueError): return None

def doc_lo(lo):
    """→ {(nguon, ban): {"rows": [...], "id": idm}} cho lô cham/lo_<lo>; kèm danh sách thiếu."""
    B = os.path.join(C, "lo_" + lo); kf = os.path.join(B, "khoa", "KHOA_MU.json")
    if not os.path.exists(kf): return {}, {"khong_co_lo": lo}
    khoa = json.load(open(kf))["id"]; out = {}; thieu = []
    for ban in ("sau", "truoc"):
        tong = os.path.join(B, "vs_" + ban, "ds_%s.txt.ket_qua" % ban, "tong.csv")
        theo = collections.defaultdict(list)
        if os.path.exists(tong):
            for r in csv.DictReader(open(tong, encoding="utf-8")): theo[r["te_bao"]].append(r)
        for tb, rs in theo.items():
            idm = tb.split("__", 1)[1]; k = khoa[idm]; assert k["ban"] == ban, tb
            out[(k["nguon"], ban)] = {"rows": rs, "id": idm}
        for idm, k in khoa.items():
            if k["ban"] == ban and (k["nguon"], ban) not in out: thieu.append([k["nguon"], ban])
    return out, {"thieu_cham": thieu}
def bang_D(rows):
    D = {}
    for r in rows:
        if r["luat"] != "D_i": continue
        st = r["trang_thai"]
        D[r["phan_quyet"]] = ("MAT", None) if st == "KHONG_CO_MESH" else ("KDL", None) if st in ("KHONG_DUNG_LAI_DUOC", "") else ("D", float(st))
    return D
def rbs(tu): return PT.ban_kinh_bao(tu["khuon"], tu["s"])
def gc_tap(D, tu_the, ten_ds):
    """G_c trên tập tên; MAT ⇒ rb·s; KDL hoặc thiếu dòng ⇒ NaN. Trả (G_c, số mất mesh, số tên thiếu)."""
    g = 0.0; n_mat = 0; n_thieu = 0
    for t in ten_ds:
        v = D.get(t)
        if v is None: n_thieu += 1; continue
        if v[0] == "D": g += max(v[1] - H, 0.0)
        elif v[0] == "MAT": n_mat += 1; g += max(rbs(tu_the[t]) - H, 0.0)
        else: return float("nan"), n_mat, n_thieu
    return (float("nan") if n_thieu else g), n_mat, n_thieu
def gc_B(tt, ten_ds):
    doi = tt.get("doi_mat_um") or {}; g = 0.0; n0 = 0
    for t in ten_ds:
        v = doi.get(t)
        if not v or v[1] is None: n0 += 1; continue
        g += max(v[1] - H, 0.0)
    return g, n0
def v_test(rows, tu_the, ten_mat_tinh):
    st = {r["luat"]: r for r in rows if r["phan_quyet"] == "su_that" and r["loai_luat"] == "su_that_v1.2"}
    tr = {k: int(v["so_ca_truot"]) for k, v in st.items()}
    pen = sum(len(PT.LUAT_MAT.get(loai(t), PT.LUAT_MAT_MAC_DINH)) for t in ten_mat_tinh)
    return sum(v for k, v in tr.items() if k not in PT.LUAT_MEM) + pen, sum(tr.get(k, 0) for k in PT.LUAT_CUNG_KP)
def gc_tap_lon(D, tu_the, M):
    return gc_tap(D, tu_the, list(M) + [t for t in tu_the if loai(t) not in NHO])[0]
def gc_ca_te_bao(rows, tu_the, bu_a):
    """G_c cả tế bào (dòng G_c của bộ kiểm + mất mesh rb·s); bu_a: cách tính bộ phận loại (a) — 'ngoai' (bỏ), 'rbs', '0'."""
    g = [r for r in rows if r["luat"] == "G_c"]; assert len(g) == 1; G = float(g[0]["trang_thai"])
    for t, v in bang_D(rows).items():
        if v[0] != "MAT": continue
        a = tu_the[t].get("t") is None
        if a and bu_a in ("ngoai", "0"): continue
        G += max(rbs(tu_the[t]) - H, 0.0)
    return G

def tt_nguon(ten):
    d = os.path.join(V, "kq", ten) if ten.split("_")[1].startswith("cp") else os.path.join(KQP, ten)
    return jl(os.path.join(d, "thong_tin.json"))
def k_cp(tt): return int(tt["cau_noi"]["k"])
def k_pl(tt): return int(tt["bao_cao_bo_tri"]["cau_hinh_chon"])
def ghep(tt_cp, tt_pl, bo_them=()):
    """→ (M_cp, M_pl, thông tin). bo_them: tên cellPACK coi như không đặt (độ nhạy mép nón)."""
    tc, tp = tt_cp["tu_the"], tt_pl["tu_the"]
    for lo_ in NHO:
        assert sum(loai(t) == lo_ for t in tc) == sum(loai(t) == lo_ for t in tp), ("so luong khac", lo_)
    M_cp = sorted(t for t, v in tc.items() if loai(t) in NHO and v.get("t") is not None and t not in bo_them)
    thieu = {lo_: sum(1 for t in tc if loai(t) == lo_) - sum(1 for t in M_cp if loai(t) == lo_) for lo_ in NHO}
    kc, kp = k_cp(tt_cp), k_pl(tt_pl)
    if kc == kp:
        for t in M_cp: assert t in tp and tp[t]["khuon"] == tc[t]["khuon"] and abs(tp[t]["s"] - tc[t]["s"]) < 1e-6, ("cung k nhung lech", t)
        M_pl = list(M_cp); cach = "cung_k"
    else:
        M_pl = []
        for lo_ in NHO:
            ds = sorted((t for t in tp if loai(t) == lo_), key=lambda t: (-rbs(tp[t]), t))
            M_pl += ds[thieu[lo_]:]
        M_pl.sort(); cach = "khac_k"
    return M_cp, M_pl, {"k_cp": kc, "k_pl": kp, "cach": cach, "thieu_theo_loai": thieu}

def mot_cap(KQ, o, nh_cp, s, nh_pl="td", ban="sau", bo_them=()):
    """một cặp (ô, seed): → dict hoặc {"place_thua": lý do}."""
    cp, pl = f"{o}_cp{nh_cp}_{s}", ten_pl(o, nh_pl, s)
    tc, tp = tt_nguon(cp), tt_nguon(pl)
    if tc is None: return {"place_thua": "cellPACK khong co ket qua nhap (" + ly_do_hong(o, s, nh_cp) + ")"}
    if tp is None: return {"place_thua": f"{pl} khong co"}
    rc, rp = KQ.get((cp, ban)), KQ.get((pl, ban))
    if rc is None or rp is None: raise RuntimeError(f"THIEU_CHAM {cp if rc is None else pl} ({ban}) — cham lai, khong phan tich (SAI-3)")
    M_cp, M_pl, info = ghep(tc, tp, bo_them)
    Dc, Dp = bang_D(rc["rows"]), bang_D(rp["rows"])
    gc, mc, _ = gc_tap(Dc, tc["tu_the"], M_cp); gp, mp, _ = gc_tap(Dp, tp["tu_the"], M_pl)
    if math.isnan(gc) or math.isnan(gp): return {"place_thua": "G_c NaN (khuon khong dung lai duoc / thieu dong D_i)", **info}
    bc, _ = gc_B(tc, M_cp); bp, _ = gc_B(tp, M_pl)
    mat_c = [t for t in M_cp if Dc.get(t, ("",))[0] == "MAT"]; mat_p = [t for t in M_pl if Dp.get(t, ("",))[0] == "MAT"]
    vc, kpc = v_test(rc["rows"], tc["tu_the"], mat_c); vp, kpp = v_test(rp["rows"], tp["tu_the"], mat_p)
    a_c = [t for t, v in tc["tu_the"].items() if v.get("t") is None]
    vc_full = v_test(rc["rows"], tc["tu_the"], mat_c + a_c)[0]
    gl_c, gl_p = gc_tap_lon(Dc, tc["tu_the"], M_cp), gc_tap_lon(Dp, tp["tu_the"], M_pl)
    return {"Gc_cp": gc, "Gc_pl": gp, "d": gc - gp, "GcB_cp": bc, "GcB_pl": bp, "dB": bc - bp, "V_cp": vc, "V_pl": vp, "dV": vc - vp,
            "dV_toan_tap": vc_full - vp, "Gc_lon_cp": gl_c, "Gc_lon_pl": gl_p, "d_ca_te_bao_cung_tap": gl_c - gl_p,
            "KP_cp": kpc, "KP_pl": kpp, "n_M": len(M_cp), "n_mat_cp": mc, "n_mat_pl": mp, **info}
def ly_do_hong(o, s, nh="chinh"):
    for f, l in ((f"{o}_{s}_CHON_HONG", "CHON_HONG"), (f"{o}_{s}_C_HONG", "C_HONG")):
        if os.path.exists(os.path.join(V, "thu", f)): return l
    n = os.path.join(H_, "hong", f"N_{o}_{s}_{nh}")
    if os.path.exists(n): return "nhap hong %s lan" % open(n).read().strip()
    return "chua co / khong ro"

def kiem_dinh(ds, seeds):
    """ds: {seed: cặp}. PLACE THUA ⇒ −(max|d| + 1). Trả số liệu Wilcoxon + mô tả."""
    ok = {s: c["d"] for s, c in ds.items() if "d" in c}; thua = [s for s in seeds if s not in ok]
    M = max([abs(x) for x in ok.values()] or [0.0]) + 1.0
    d = [ok[s] for s in seeds if s in ok] + [-M] * len(thua)
    if not d: return None
    p, W, n_, hoa = PT.wilcoxon_chinh_xac(d)
    return {"n": len(d), "n_khac0": n_, "bo_duoi_tol": hoa, "W_cong": W, "p": p, "place_thua": thua,
            "thang": sum(x >= TOL for x in d), "hoa": sum(abs(x) < TOL for x in d), "thua": sum(x <= -TOL for x in d),
            "d": tv(d), "d_chi_cap_du": tv(list(ok.values())), "d_theo_seed": {str(s): (ok[s] if s in ok else "PLACE_THUA") for s in seeds}}
def bat_dong(ds):
    r = [c for c in ds.values() if "d" in c]
    a = [c["Gc_cp"] for c in r] + [c["Gc_pl"] for c in r]; b = [c["GcB_cp"] for c in r] + [c["GcB_pl"] for c in r]
    rho = PT.spearman(a, b) if len(a) >= 3 else float("nan")
    lech = [s for s, c in ds.items() if "d" in c and (0 if abs(c["d"]) < TOL else np.sign(c["d"])) != (0 if abs(c["dB"]) < TOL else np.sign(c["dB"]))]
    return {"spearman": rho, "so_seed_lech_dau": len(lech), "seed_lech": [int(s) for s in lech], "bat_dong": bool((not (rho >= 0.8)) or len(lech) > 2)}

def thu_dir(o, s, k): return os.path.join(V, "thu", f"{o}_{s}_k{k}")
def chon_k(o, s):
    j = jl(os.path.join(V, "thu", f"{o}_{s}_chon.json")); return None if j is None else j
def cp_log(d):
    t = open(d, errors="replace").read() if os.path.exists(d) else ""
    m = re.search(r"GPU_THONG_KE (\{.*?\})", t); bd = ast.literal_eval(m.group(1)) if m else None
    x = re.search(r"CP_XONG dat (\d+) / (\d+) giay_luoi ([\d.]+) giay_dat ([\d.]+)", t)
    r = re.findall(r"THOI_GIAN ([\d.]+) s (\d+) KB", t)
    return {"bo_dem": bd, "dat": int(x.group(1)) if x else None, "yeu_cau": int(x.group(2)) if x else None,
            "giay_luoi": float(x.group(3)) if x else None, "giay_dat": float(x.group(4)) if x else None,
            "giay_tuong": float(r[-1][0]) if r else None, "rss_gb": int(r[-1][1]) / 1048576 if r else None}
def cap_da_chay():
    out = set()
    for f in glob.glob(os.path.join(H_, "G_*_k*.cap")):
        o_, s_, k_ = open(f).read().split(); out.add((o_, int(s_), int(k_)))
    return out
def mo_ta_hang(o):
    """K hiệu lực, xuất lại 0,125, XUAT_HONG, bộ đếm v4, giờ từng pha, H5a (từ tu_the_cp.json cấu hình chọn)."""
    out = {"seed": {}, "dat_theo_loai": {lo_: [0, 0] for lo_ in NHO}, "bo_dem_khac0": [], "gio": collections.defaultdict(list)}
    CAP = cap_da_chay()
    for s in SEED:
        r = {"K_hieu_luc": 0, "K_v4_xong": 0, "xuat_lai_0125": 0, "xuat_hong": 0, "v4_hong": 0}
        for k in range(4):
            d = thu_dir(o, s, k)
            if os.path.exists(os.path.join(d, "THU_H0125")): r["xuat_lai_0125"] += 1
            if os.path.exists(os.path.join(d, "XUAT_HONG")): r["xuat_hong"] += 1; continue
            r["K_hieu_luc"] += 1
            dv = jl(os.path.join(d, "dau_vao_cellpack.json"))
            if dv: out["gio"]["giay_xuat"].append(dv.get("giay_xuat"))
            L = cp_log(os.path.join(d, "cp_chinh", "cp.log"))
            if os.path.exists(os.path.join(H_, "xong", f"V_{o}_{s}_{k}_chinh")):
                r["K_v4_xong"] += 1
                for g_ in ("giay_luoi", "giay_dat", "giay_tuong", "rss_gb"): out["gio"]["v4_" + g_].append(L[g_])
                if L["bo_dem"] is None or any(L["bo_dem"].get(b, 0) for b in BO_DEM):
                    out["bo_dem_khac0"].append({"seed": s, "k": k, "bo_dem": {b: (L["bo_dem"] or {}).get(b) for b in BO_DEM}, "trong_cap": (o, s, k) in CAP,
                                                "nhan": "co cap tuong duong" if (o, s, k) in CAP else "khong chung nhan"})
            else: r["v4_hong"] += 1
            cj = jl(os.path.join(d, f"chon_k{k}.json"))
            if cj: out["gio"]["giay_cham_chon"].append(cj.get("giay_cham_chon"))
        ch = chon_k(o, s); r["k_chon"] = ch["k_chon"] if ch else None
        if ch:
            d = thu_dir(o, s, ch["k_chon"]); dv = jl(os.path.join(d, "dau_vao_cellpack.json")); tc = jl(os.path.join(d, "cp_chinh", "tu_the_cp.json"))
            if dv and tc:
                for b in dv["bo_phan"]:
                    lo_ = b["loai"]; out["dat_theo_loai"][lo_][1] += 1; out["dat_theo_loai"][lo_][0] += b["ten"] in tc["dat"]
                r["dat_k_chon"] = len(tc["dat"]); r["yeu_cau"] = len(dv["bo_phan"])
            d0 = jl(os.path.join(thu_dir(o, s, 0), "cp_chinh", "tu_the_cp.json"))
            r["dat_k0_K1"] = len(d0["dat"]) if d0 else None
        tn = tt_nguon(f"{o}_cpchinh_{s}")
        if tn:
            out["gio"]["giay_nhap_tong"].append(tn.get("giay_tong")); r["lech_khoi_L_xuat_nhap"] = tn.get("bao_cao_bo_tri", {}).get("cellpack", {}).get("lech_khoi_L")
            r["phat_cung_truoc_chieu"] = tn.get("bao_cao_bo_tri", {}).get("phat_cung_day_du")
        nn = tt_nguon(ten_pl(o, "nn1x", s))
        if nn and k_pl(nn) in range(4):
            dv = jl(os.path.join(thu_dir(o, s, k_pl(nn)), "dau_vao_cellpack.json"))
            if dv:
                r["lech_khoi_L_vs_nn1x"] = max(max(float(np.abs(np.array(e["t"]) - np.array(nn["tu_the"][t]["t"])).max()), float(np.abs(np.array(e["R"]) - np.array(nn["tu_the"][t]["R"])).max()))
                                              for t, e in dv["lon_co_dinh"].items())
        out["seed"][str(s)] = r
    goc = []
    for f in glob.glob(os.path.join(H_, "G_*_k*.cap")):
        o_, s_, k_ = open(f).read().split()
        if o_ != o: continue
        Lg = cp_log(os.path.join(thu_dir(o_, s_, k_), "cp_goc", "cp.log")); Lv = cp_log(os.path.join(thu_dir(o_, s_, k_), "cp_chinh", "cp.log"))
        goc.append({"seed": int(s_), "k": int(k_), "goc_giay_tuong": Lg["giay_tuong"], "goc_rss_gb": Lg["rss_gb"], "v4_giay_tuong": Lv["giay_tuong"],
                    "ti_so_v4_goc": (Lv["giay_tuong"] / Lg["giay_tuong"]) if Lg["giay_tuong"] and Lv["giay_tuong"] else None,
                    "so_sanh": (jl(os.path.join(H_, "ss", f"{o_}_{s_}_k{k_}.json")) or {}).get("ket_luan")})
    out["cap_tuong_duong"] = goc
    out["H5a"] = {lo_: {"dat": v[0], "yeu_cau": v[1], "ti_le": v[0] / v[1] if v[1] else None, "ktc95_CP": cp_ci(v[0], v[1])} for lo_, v in out["dat_theo_loai"].items()}
    ks = [min(r["K_hieu_luc"], r["K_v4_xong"]) for r in out["seed"].values()]
    out["tom_tat"] = {"seed_K4": sum(k == 4 for k in ks), "xuat_lai_0125": sum(r["xuat_lai_0125"] for r in out["seed"].values()),
                      "xuat_hong": sum(r["xuat_hong"] for r in out["seed"].values()), "so_luot_bo_dem_khac0": len(out["bo_dem_khac0"])}
    out["gio"] = {k: tv(v) for k, v in out["gio"].items()}
    return out
def n1_cau_hinh(o):
    rs = [jl(os.path.join(C, "lan_cau_hinh", f"{o}_{s}_k{k}.json")) for s in SEED for k in range(4)]; rs = [r for r in rs if r]
    return {"so_cau_hinh_do": len(rs), "so_vuot_h2": sum(r["cau_hinh"]["vuot_h2"] for r in rs), "lan_max": tv([r["cau_hinh"]["lan_max_um"] for r in rs]),
            "ti_le_F_mat": tv([r["cau_hinh"]["ti_le_F_mat"] for r in rs]), "lech_khoi_L_max": max([r["lech_khoi_L_vs_xuat"] for r in rs]) if rs else None,
            "bien_hop_cham_F": sum(r["bien_hop_cham_F"] for r in rs)}
def mep_non(o, s, nh="chinh"):
    """Quyết định 3 (20/11): danh sách bỏ khỏi M′ = bộ phận chạm vùng mất do phép mở (đệm 0,052 µm), ĐÓNG BĂNG từ tu_the_cp.json của cấu hình
    được nhập, TRƯỚC khi nhập (mtime tệp M′ < mtime thong_tin.json của kết quả nhập và md5 tu_the_cp khớp). Không đạt ⇒ None (không có M′)."""
    tt = tt_nguon(f"{o}_cp{nh}_{s}")
    if tt is None: return None
    k = int(tt["cau_noi"]["k"]); f = os.path.join(M_PHAY, f"{o}_{s}_k{k}_{nh}.json"); j = jl(f)
    if j is None or "M_phay" not in j: return None
    fc = os.path.join(V, "thu", f"{o}_{s}_k{k}", "cp_" + nh, "tu_the_cp.json"); ft = os.path.join(V, "kq", f"{o}_cp{nh}_{s}", "thong_tin.json")
    if md5(fc) != j["M_phay"]["md5_tu_the_cp"] or os.path.getmtime(f) >= os.path.getmtime(ft): return None
    if j.get("md5_ma") != MD5_DO_LAN: return None
    man = {}; mf = os.path.join(M_PHAY, "MANIFEST.txt")
    for l in (open(mf) if os.path.exists(mf) else []):
        p_ = l.split()
        if len(p_) == 5: man[p_[4]] = (p_[0] + " " + p_[1], p_[2], p_[3])
    m_ = man.get(f"{o}_{s}_k{k}_{nh}")
    if m_ is None or m_[1] != md5(f) or m_[2] != md5(fc): return None
    bd = gio_bat_dau_nhap(o, s, nh)
    if bd is None or not (m_[0] < bd): return None
    return set(j["M_phay"]["bo_khoi_M_phay"])
def gio_bat_dau_nhap(o, s, nh):
    """giờ BAT_DAU sớm nhất của việc N_<ô>_<seed>_<nhãn> trong hang/tien_do.txt (dạng 'MM-DD HH:MM:SS'), đổi sang 'YYYY-MM-DD HH:MM:SS'."""
    f = os.path.join(H_, "tien_do.txt"); g = None
    for l in (open(f, errors="replace") if os.path.exists(f) else []):
        if f" BAT_DAU N_{o}_{s}_{nh} " in l + " " and re.match(r"\d\d-\d\d \d\d:\d\d:\d\d", l):
            t = "2026-" + l[:14]; g = t if g is None or t < g else g
    return g

def bao_theo_khuon(o):
    """Quyết định 2: bao sai số dựng mesh = max D^K TRƯỚC chiếu theo LOẠI KHUÔN trên 30 tế bào lượt chính cùng ô (td, nn1x, un1x × 4001–4010;
    chấm thô cham_chinh_dv, khoá đã mở 15/11). Tham chiếu cố định (cùng bộ dựng mesh phần C) cho mọi biến thể của mục 4."""
    k = json.load(open(os.path.join(CC_CHINH, "khoa", "KHOA_MU.json")))["id"]; ma = {(v["nguon"], v["ban"]): i for i, v in k.items()}
    env = collections.defaultdict(float); n = 0
    for nh in ("td",) + DOI:
        for s in range(4001, 4011):
            i = ma.get((f"{o}_{nh}_{s}", "truoc")); f = os.path.join(CHAM_RA, f"mu_truoc__{i}", "kiem3.csv") if i else None
            if not f or not os.path.exists(f): continue
            n += 1
            for r in csv.DictReader(open(f, encoding="utf-8")):
                if r["luat"] == "D_i" and r["trang_thai"] not in ("KHONG_CO_MESH", "KHONG_DUNG_LAI_DUOC", ""):
                    env[r["bo_phan_truot"]] = max(env[r["bo_phan_truot"]], float(r["trang_thai"]))
    return dict(env), n
def quy_uoc_quay(dv, kc):
    """Cổng tất định (rà soát v3 SAI-2): số cặp cầu (cây cầu dau_vao, bán kính − a = thân) của hai bộ phận KHÁC nhau chồng nhau dưới
    x = R·c + t (quy ước bộ dựng) phải ít hơn dưới x = Rᵀ·c + t. Chỉ đọc tu_the_cp + dau_vao."""
    from scipy.spatial import cKDTree
    out = {}
    for ten_qu, tr in (("Rc_t", False), ("RTc_t", True)):
        C_ = []; r_ = []; id_ = []
        for i_, b in enumerate(dv["bo_phan"]):
            v = kc["dat"].get(b["ten"])
            if v is None or b["ngan"] != "te_bao_chat": continue
            R = np.asarray(v["R"], float); R = R.T if tr else R
            tam = np.asarray(b["tam"], float) @ R.T + np.asarray(v["t"], float); rr = np.asarray(b["r"], float) - dv["nguong"]["a_nho"]
            C_.append(tam); r_.append(rr); id_ += [i_] * len(rr)
        C_ = np.concatenate(C_); r_ = np.concatenate(r_); id_ = np.asarray(id_); t_ = cKDTree(C_)
        n = 0
        for a, b in t_.query_pairs(2 * float(r_.max())):
            if id_[a] != id_[b] and np.linalg.norm(C_[a] - C_[b]) < r_[a] + r_[b] - 1e-6: n += 1
        out[ten_qu] = n
    out["dat"] = out["Rc_t"] < out["RTc_t"] or (out["Rc_t"] == 0 and out["RTc_t"] == 0)
    out["khong_xac_dinh"] = out["Rc_t"] == 0 and out["RTc_t"] == 0
    return out
def khop_tu_the(o, s, nh="chinh", quy_uoc=False):
    """Quyết định 2 (i), CỔNG TẤT ĐỊNH CỨNG (không có dấu PI nào vượt qua): tu_the của kết quả nhập khớp tu_the_cp.json (t, R ≤ 1e-6), khuôn + s khớp
    dau_vao (≤ 1e-6), tập không đặt trùng, không thiếu/thừa tên, bộ phận lớn khớp tệp xuất; k khớp chon.json; (quy_uoc) quy ước quay đúng."""
    tt = tt_nguon(f"{o}_cp{nh}_{s}")
    if tt is None: return {"dat": False, "ly_do": "khong co ket qua nhap"}
    k = int(tt["cau_noi"]["k"]); d = os.path.join(V, "thu", f"{o}_{s}_k{k}"); kc = jl(os.path.join(d, "cp_" + nh, "tu_the_cp.json")); dv = jl(os.path.join(d, "dau_vao_cellpack.json"))
    if kc is None or dv is None: return {"dat": False, "ly_do": "thieu tu_the_cp/dau_vao"}
    lech = 0.0; sai = []
    if nh == "chinh":
        ch = chon_k(o, s)
        if ch is None or int(ch["k_chon"]) != k: sai.append("k_khac_chon_json")
    ten_dv = {b["ten"]: b for b in dv["bo_phan"]}
    for ten, b in ten_dv.items():
        tu = tt["tu_the"].get(ten)
        if tu is None: sai.append(ten); continue
        if tu.get("khuon") != b["khuon"] or abs(float(tu["s"]) - float(b["s"])) > 1e-6: sai.append(ten + ":khuon/s"); continue
        if ten in kc["dat"]:
            if tu.get("t") is None: sai.append(ten); continue
            lech = max(lech, float(np.abs(np.array(tu["t"]) - np.array(kc["dat"][ten]["t"])).max()), float(np.abs(np.array(tu["R"]) - np.array(kc["dat"][ten]["R"])).max()))
        elif tu.get("t") is not None: sai.append(ten)
    for ten, e in dv["lon_co_dinh"].items():
        tu = tt["tu_the"].get(ten)
        if tu is None or tu.get("t") is None or tu.get("khuon") != e["khuon"] or abs(float(tu["s"]) - float(e["s"])) > 1e-6: sai.append(ten + ":lon"); continue
        lech = max(lech, float(np.abs(np.array(tu["t"]) - np.array(e["t"])).max()), float(np.abs(np.array(tu["R"]) - np.array(e["R"])).max()))
    thua = sorted(set(tt["tu_the"]) - set(ten_dv) - set(dv["lon_co_dinh"]))
    if thua: sai += ["thua:" + t for t in thua]
    out = {"dat": lech <= 1e-6 and not sai, "lech_max": lech, "sai_ten": sai[:20], "so_sai": len(sai)}
    if quy_uoc: out["quy_uoc_quay"] = quy_uoc_quay(dv, kc); out["dat"] = out["dat"] and out["quy_uoc_quay"]["dat"]
    return out
def buoc_kiem_truoc():
    """Quyết định 2 (20/11): cổng = khớp tư thế (tất định) ∧ mọi bộ phận đã đặt có D^K trước chiếu ≤ bao theo loại khuôn của lượt chính + 0,01 µm;
    xét 2 seed đầu mỗi ô (§6); báo mọi seed. Chỉ đọc bản TRƯỚC chiếu; không G_c."""
    KQ, meta = doc_lo("chinh"); out = {"meta": meta, "dung_sai_um": DEN, "o": {}}; dat_all = True; tu_the_all = True
    for o in O_CP:
        env, n_ref = bao_theo_khuon(o); out["o"][o] = {"bao_theo_khuon": env, "so_te_bao_tham_chieu": n_ref, "seed": {}}
        if n_ref != 30: dat_all = False; out["o"][o]["loi"] = "so te bao tham chieu %d != 30" % n_ref
        for s in SEED:
            ten = f"{o}_cpchinh_{s}"; x = KQ.get((ten, "truoc")); tt = tt_nguon(ten); r = {"cong": s in SEED[:2]}
            r["khop_tu_the"] = khop_tu_the(o, s, quy_uoc=s in SEED[:2])
            if tt is not None and not r["khop_tu_the"]["dat"]: tu_the_all = False
            if x is None or tt is None: r["bao"] = {"dat": False, "ly_do": "khong co ban truoc chieu da cham"}
            else:
                D = bang_D(x["rows"]); vuot = []; thieu_khuon = set()
                mat_truoc = [t for t, v in tt["tu_the"].items() if v.get("t") is not None and D.get(t, ("",))[0] != "D"]
                for t, v in tt["tu_the"].items():
                    if v.get("t") is None or D.get(t, ("",))[0] != "D": continue
                    kh = v["khuon"]
                    if kh not in env: thieu_khuon.add(kh); continue
                    if D[t][1] > env[kh] + DEN: vuot.append([t, kh, D[t][1], env[kh]])
                r["bao"] = {"dat": not vuot and not thieu_khuon and not mat_truoc, "mat_hoac_kdl_truoc_chieu": mat_truoc[:20], "so_vuot": len(vuot), "vuot": vuot[:30], "khuon_khong_co_tham_chieu": sorted(thieu_khuon)}
            r["dat"] = r["khop_tu_the"]["dat"] and r["bao"]["dat"]
            if r["cong"] and not r["dat"]: dat_all = False
            out["o"][o]["seed"][str(s)] = r
    k_ = json.load(open(os.path.join(CC_CHINH, "khoa", "KHOA_MU.json")))["id"]; ma_ = {(v["nguon"], v["ban"]): i for i, v in k_.items()}
    vs = {}
    for o in O_CP:
        for s in SEED:
            pl = ten_pl(o, "td", s); x = KQ.get((pl, "truoc")); i = ma_.get((pl, "truoc"))
            f = os.path.join(CHAM_RA, f"mu_truoc__{i}", "kiem3.csv") if i else None
            if i is None: vs[pl] = "khong_ap_dung (khong co cham tho luot chinh)"; continue
            if x is None or not os.path.exists(f): vs[pl] = "THIEU"; dat_all = False; continue
            Dv, Dt = bang_D(x["rows"]), bang_D(list(csv.DictReader(open(f, encoding="utf-8"))))
            lech = [t for t in set(Dv) | set(Dt) if (Dv.get(t, ("?",))[0] != Dt.get(t, ("?",))[0]) or (Dv.get(t, ("",))[0] == "D" and not (abs(Dv[t][1] - Dt[t][1]) <= 1e-6))]
            vs[pl] = {"so_bo_phan": len(Dt), "so_lech": len(lech), "lech": sorted(lech)[:20],
                      "max_abs_dD": max([abs(Dv[t][1] - Dt[t][1]) for t in Dt if Dt[t][0] == "D" and Dv.get(t, ("",))[0] == "D"] or [0.0])}
            if lech: dat_all = False
    out["ve_sinh_vs_tho_truoc_chieu"] = vs
    out["dat"] = dat_all and tu_the_all; out["khop_tu_the_moi_te_bao"] = tu_the_all
    os.makedirs(RA, exist_ok=True); json.dump(out, open(os.path.join(RA, "KIEM_TRUOC.json"), "w"), ensure_ascii=False, indent=1, default=_jd)
    print("KIEM_TRUOC dat" if dat_all else "KIEM_TRUOC CHUA_DAT", json.dumps({o: {s: (r["khop_tu_the"]["dat"], r["bao"].get("so_vuot")) for s, r in d["seed"].items()} for o, d in out["o"].items()}))
    sys.exit(0 if (dat_all and tu_the_all) else (4 if not tu_the_all else 3))

def buoc_chinh():
    if os.environ.get("CP_THU") != "1":
        for f in (os.path.join(H_, "TUONG_DUONG_DAT"), os.path.join(C, "DUYET_MA_CHAM")): assert os.path.exists(f), "thieu " + f
        kt = jl(os.path.join(RA, "KIEM_TRUOC.json")); assert (kt and kt.get("dat")) or os.path.exists(os.path.join(C, "DUYET_TU_THE_CP")), "kiem truoc chieu chua dat / chua co DUYET_TU_THE_CP"
    KQ, meta = doc_lo("chinh"); KQd, meta_d = doc_lo("do_nhay")
    if meta.get("thieu_cham") or meta.get("khong_co_lo"): sys.exit("THIEU_CHAM lo chinh: %s — cham lai, khong phan tich (SAI-3)" % meta)
    KQall = {**KQ, **KQd}
    KQS = {"ban": "phan_tich_cp", "md5_ma": md5(os.path.abspath(__file__)), "meta_chinh": meta, "meta_do_nhay": meta_d,
           "co_TUONG_DUONG_DAT": os.path.exists(os.path.join(H_, "TUONG_DUONG_DAT")), "KHAC": os.path.exists(os.path.join(H_, "KHAC")),
           "kiem_truoc": jl(os.path.join(RA, "KIEM_TRUOC.json")) is not None}
    H5 = {}; CHUOI = []; tiep = True
    for o in O_CP:
        ds = {s: mot_cap(KQ, o, "chinh", s) for s in SEED}
        kd = kiem_dinh(ds, SEED); bd = bat_dong(ds)
        e = {"cap": {str(s): c for s, c in ds.items()}, "kiem_dinh": kd, "bat_dong": bd}
        if bd["bat_dong"]: e["trang_thai_chuoi"] = "BAT_DONG_chi_mo_ta"
        elif not tiep: e["trang_thai_chuoi"] = "KHONG_KIEM (o truoc khong dat)"
        else:
            dat = kd is not None and kd["p"] <= ALPHA; e["trang_thai_chuoi"] = "DAT" if dat else "KHONG_DAT"; tiep = dat
        e["dat_alpha_0_0167"] = bool(kd and kd["p"] <= ALPHA3)
        e["H6_V_test"] = kiem_dinh({s: {"d": c["dV"]} for s, c in ds.items() if "dV" in c}, SEED)
        e["H6_V_test_toan_tap_yeu_cau"] = kiem_dinh({s: {"d": c["dV_toan_tap"]} for s, c in ds.items() if "dV" in c}, SEED)
        e["ca_te_bao_cung_tap"] = kiem_dinh({s: {"d": c["d_ca_te_bao_cung_tap"]} for s, c in ds.items() if "d" in c}, SEED)
        med = kd["d"]["tv"] if kd and kd["d"] else None
        mn = {s: mep_non(o, s) for s in SEED}
        kdM = kiem_dinh({s: (mot_cap(KQ, o, "chinh", s, bo_them=mn[s]) if mn[s] is not None else {"place_thua": "khong co M' dong bang truoc nhap"}) for s in SEED}, SEED)
        _c = [(len(mn[s]), ds[s]["n_M"]) for s in SEED if mn[s] is not None and "n_M" in ds[s]]
        ti_bo = [sum(a for a, b in _c) / max(sum(b for a, b in _c), 1)] if _c else []
        e["M_phay"] = {"kiem_dinh": kdM, "so_bo": {str(s): (len(v) if v is not None else None) for s, v in mn.items()}, "ti_le_bo_max": max(ti_bo) if ti_bo else None,
                       "seed_khong_co_M_phay": [s for s in SEED if mn[s] is None]}
        dat_Mp = bool(kdM and kdM["p"] <= ALPHA and kdM["d"] and kdM["d"]["tv"] >= H and ti_bo and max(ti_bo) <= 0.20 and not e["M_phay"]["seed_khong_co_M_phay"])
        e["cau_chu"] = ("H5 dat tren M va M', trung vi d >= h — CHUA DU DIEU KIEN PHAT BIEU 'PLACE it meo hon' (thieu do nhay (ii), §1)" if e["trang_thai_chuoi"] == "DAT" and med is not None and med >= H and dat_Mp else
                        "H5 dat tren M nhung KHONG dat tren M' (hoac bo > 20 %) — cam cau 'it meo hon' (quyet dinh 3)" if e["trang_thai_chuoi"] == "DAT" and med is not None and med >= H else
                        "co y nghia thong ke nhung duoi mot buoc luoi" if e["trang_thai_chuoi"] == "DAT" else "khong phat bieu 'it meo hon'")
        e["do_nhay"] = {}
        e["do_nhay"]["chi_seed_K4"] = (lambda ss: kiem_dinh({s: ds[s] for s in ss}, ss) if ss else None)(
            [s for s in SEED if sum(os.path.exists(os.path.join(H_, "xong", f"V_{o}_{s}_{k}_chinh")) for k in range(4)) == 4])
        mnb = {s: (lambda j: None if j is None else set(j["bo_phan"]["mep_non_tren_h2"]))(jl(os.path.join(C, "lan_bo_phan", f"{o}_cpchinh_{s}.json"))) for s in SEED}
        e["do_nhay"]["bo_lan_than_tren_h2_N1b"] = kiem_dinh({s: (mot_cap(KQ, o, "chinh", s, bo_them=mnb[s]) if mnb[s] is not None else {"place_thua": "khong co do lan"}) for s in SEED}, SEED)
        e["do_nhay"]["truoc_chieu"] = kiem_dinh({s: mot_cap(KQ, o, "chinh", s, ban="truoc") for s in SEED}, SEED)
        for nh in DOI:
            try: e["do_nhay"]["so_" + nh] = kiem_dinh({s: mot_cap(KQ_THO_DOI(KQ), o, "chinh", s, nh_pl=nh) for s in SEED}, SEED)
            except RuntimeError as ex: e["do_nhay"]["so_" + nh] = {"khong_co": str(ex)}
            if e["do_nhay"]["so_" + nh] and "khong_co" not in e["do_nhay"]["so_" + nh]: e["do_nhay"]["so_" + nh]["nhan"] = "doi chung cham THO luot chinh, cellPACK ve sinh"
        for nh in NHAN_DN:
            ss = SEED[:5] if nh != "h1" else SEED[:3]
            e["do_nhay"]["cp_" + nh] = kiem_dinh({s: mot_cap(KQall, o, nh, s) for s in ss}, ss) if KQd else None
            if e["do_nhay"]["cp_" + nh]: e["do_nhay"]["cp_" + nh]["nhan"] = "v4, chua chung nhan tuong duong o cau hinh nay"
        bu = {"bu_0": {}, "bu_rbs": {}}
        tho = {}
        for s in SEED:
            tc, tp = tt_nguon(f"{o}_cpchinh_{s}"), tt_nguon(ten_pl(o, "td", s))
            rc, rp = KQ.get((f"{o}_cpchinh_{s}", "sau")), KQ.get((ten_pl(o, "td", s), "sau"))
            if tc and tp and rc and rp:
                gp = gc_ca_te_bao(rp["rows"], tp["tu_the"], "ngoai")
                bu["bu_0"][s] = {"d": gc_ca_te_bao(rc["rows"], tc["tu_the"], "0") - gp}
                bu["bu_rbs"][s] = {"d": gc_ca_te_bao(rc["rows"], tc["tu_the"], "rbs") - gp}
            c = ds[s]
            if "d" in c and tp:
                tho[str(s)] = {"Gc_pl_ve_sinh": c["Gc_pl"], "Gc_pl_tho": gc_tho(ten_pl(o, "td", s), tc, tp)}
        for k_, v in bu.items(): e["do_nhay"][k_] = kiem_dinh(v, SEED)
        e["place_tho_vs_ve_sinh"] = tho
        e["hang"] = mo_ta_hang(o); e["N1_cau_hinh"] = n1_cau_hinh(o)
        H5[o] = e
    KQS["H5"] = H5
    os.makedirs(RA, exist_ok=True)
    json.dump(KQS, open(os.path.join(RA, "KQ_CP.json"), "w"), ensure_ascii=False, indent=1, default=_jd)
    viet_md(KQS); viet_tex(KQS)
    print("XONG phan_tich_cp", {o: (v["trang_thai_chuoi"], v["kiem_dinh"] and round(v["kiem_dinh"]["p"], 4)) for o, v in H5.items()})
_KHOA_THO = None
def KQ_THO_DOI(KQ):
    global _KHOA_THO
    if _KHOA_THO is None:
        k = json.load(open(os.path.join(CC_CHINH, "khoa", "KHOA_MU.json")))["id"]; _KHOA_THO = {(v["nguon"], v["ban"]): i for i, v in k.items()}
    out = dict(KQ)
    for (ng, ban), i in _KHOA_THO.items():
        if any(f"_{nh}_" in ng for nh in DOI):
            f = os.path.join(CHAM_RA, f"mu_{ban}__{i}", "kiem3.csv")
            if os.path.exists(f): out[(ng, ban)] = {"rows": list(csv.DictReader(open(f, encoding="utf-8"))), "id": i}
    return out
def _jd(x):
    if isinstance(x, np.bool_): return bool(x)
    if isinstance(x, np.floating): return float(x)
    if isinstance(x, np.integer): return int(x)
    if isinstance(x, set): return sorted(x)
    return str(x)
def gc_tho(pl, tc, tp):
    """G_c của PLACE trên cùng M_pl, đọc chấm THÔ lượt chính (cham_chinh_dv, khoá đã mở từ 15/11)."""
    global _KHOA_THO
    if _KHOA_THO is None:
        k = json.load(open(os.path.join(CC_CHINH, "khoa", "KHOA_MU.json")))["id"]; _KHOA_THO = {(v["nguon"], v["ban"]): i for i, v in k.items()}
    i = _KHOA_THO.get((pl, "sau"))
    f = os.path.join(CHAM_RA, f"mu_sau__{i}", "kiem3.csv") if i else None
    if not f or not os.path.exists(f): return None
    _, M_pl, _ = ghep(tc, tp)
    return gc_tap(bang_D(list(csv.DictReader(open(f, encoding="utf-8")))), tp["tu_the"], M_pl)[0]

def f3(x, n=3): return "—" if x is None or (isinstance(x, float) and math.isnan(x)) else (("%." + str(n) + "f") % x)
def tvs(t, n=3): return "—" if not t else "%s [%s; %s]" % (f3(t["tv"], n), f3(t["q1"], n), f3(t["q3"], n))
def viet_md(K):
    L = ["# Kết quả nhánh cellPACK (phan_tich_cp.py, md5 %s)" % K["md5_ma"][:8], "",
         "TUONG_DUONG_DAT: %s; KHAC: %s; kiểm trước chiếu đã chạy: %s" % (K["co_TUONG_DUONG_DAT"], K["KHAC"], K["kiem_truoc"]), "",
         "## H5 (chuỗi D5 → D4A2, α = 0,025; d = G_c(cellPACK) − G_c(PLACE td) trên tập M, µm)", "",
         "| Ô | n | n′ | thắng/hoà/thua | PLACE THUA | trung vị d [Q1; Q3] | p (một phía) | bất đồng (ρ, lệch dấu) | chuỗi | α=0,0167 | câu chữ |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for o, e in K["H5"].items():
        k = e["kiem_dinh"] or {}; b = e["bat_dong"]
        L.append("| %s | %s | %s | %s/%s/%s | %d | %s | %s | %s, %d | %s | %s | %s |" % (o, k.get("n"), k.get("n_khac0"), k.get("thang"), k.get("hoa"), k.get("thua"), len(k.get("place_thua", [])),
                 tvs(k.get("d")), f3(k.get("p"), 4), f3(b["spearman"]), b["so_seed_lech_dau"], e["trang_thai_chuoi"], e["dat_alpha_0_0167"], e["cau_chu"]))
    L += ["", "## H5a — tỉ lệ cellPACK đặt được (gộp 10 seed, KTC 95 % Clopper–Pearson)", "", "| Ô | loại | đặt / yêu cầu | tỉ lệ | KTC 95 % |", "|---|---|---|---|---|"]
    for o, e in K["H5"].items():
        for lo_, v in e["hang"]["H5a"].items(): L.append("| %s | %s | %d / %d | %s | [%s; %s] |" % (o, lo_, v["dat"], v["yeu_cau"], f3(v["ti_le"]), f3(v["ktc95_CP"][0]), f3(v["ktc95_CP"][1])))
    L += ["", "## Mô tả và độ nhạy (không p cho chuỗi)", "", "| Ô | thước | n | thắng/hoà/thua | trung vị d [Q1; Q3] | p mô tả |", "|---|---|---|---|---|---|"]
    for o, e in K["H5"].items():
        for k_ in ("H6_V_test_toan_tap_yeu_cau", "ca_te_bao_cung_tap"):
            v = e.get(k_)
            if v: L.append("| %s | %s | %d | %d/%d/%d | %s | %s |" % (o, k_, v["n"], v["thang"], v["hoa"], v["thua"], tvs(v["d"]), f3(v["p"], 4)))
        v = e["H6_V_test"]
        if v: L.append("| %s | H6 V_test (dV; PLACE THUA cho cặp thiếu) | %d | %d/%d/%d | %s | %s |" % (o, v["n"], v["thang"], v["hoa"], v["thua"], tvs(v["d"], 1), f3(v["p"], 4)))
        for k_, v in e["do_nhay"].items():
            if isinstance(v, dict) and "p" in v: L.append("| %s | %s | %d | %d/%d/%d | %s | %s |" % (o, k_, v["n"], v["thang"], v["hoa"], v["thua"], tvs(v["d"]), f3(v["p"], 4)))
    L += ["", "## Hàng, giờ, bộ đếm, N1", ""]
    for o, e in K["H5"].items():
        h = e["hang"]; L.append("- %s: seed K hiệu lực = 4: %d/10; xuất lại 0,125: %d; XUAT_HONG: %d; lượt v4 có bộ đếm ≠ 0 (không chứng nhận): %d" % (
            o, h["tom_tat"]["seed_K4"], h["tom_tat"]["xuat_lai_0125"], h["tom_tat"]["xuat_hong"], h["tom_tat"]["so_luot_bo_dem_khac0"]))
        L.append("  - giờ (trung vị [Q1; Q3], s): " + "; ".join("%s %s" % (k_, tvs(v, 1)) for k_, v in h["gio"].items()))
        n1 = e["N1_cau_hinh"]; L.append("  - N1: %d cấu hình đo; vượt h/2: %d; độ lấn max %s µm; tỉ lệ F mất %s" % (n1["so_cau_hinh_do"], n1["so_vuot_h2"], tvs(n1["lan_max"]), tvs(n1["ti_le_F_mat"], 5)))
    open(os.path.join(RA, "KQ_CP.md"), "w").write("\n".join(L) + "\n")
def viet_tex(K):
    T = [r"% bang_cp.tex — sinh bởi phan_tich_cp.py md5 " + K["md5_ma"][:8] + r" (booktabs, kiểu Computers \& Graphics)",
         r"\begin{table}[t]\centering\small", r"\caption{External baseline: cellPACK (a576c5c, accelerated port v4) vs.\ PLACE on the same set of placed parts. "
         r"$d = G_c(\mathrm{cellPACK}) - G_c(\mathrm{PLACE})$ per seed ($\mu$m), exact one-sided Wilcoxon signed-rank, fixed chain D5$\rightarrow$D4A2, $\alpha=0.025$.}",
         r"\label{tab:cellpack}", r"\begin{tabular}{lrrcccl}", r"\toprule", r"Cell & $n$ & $n'$ & W/T/L & median $d$ [IQR] & $p$ & Decision \\", r"\midrule"]
    for o, e in K["H5"].items():
        k = e["kiem_dinh"] or {}; d = k.get("d")
        T.append(r"%s & %s & %s & %s/%s/%s & %s & %s & %s \\" % (o, k.get("n"), k.get("n_khac0"), k.get("thang"), k.get("hoa"), k.get("thua"),
                 ("%s [%s, %s]" % (f3(d["tv"]), f3(d["q1"]), f3(d["q3"]))) if d else "--", f3(k.get("p"), 4), e["trang_thai_chuoi"].replace("_", r"\_")))
    T += [r"\bottomrule", r"\end{tabular}", r"\end{table}", "",
          r"\begin{table}[t]\centering\small", r"\caption{cellPACK placement rate (pooled over 10 seeds; 95\% Clopper--Pearson CI).}", r"\label{tab:cellpack_place}",
          r"\begin{tabular}{llrrc}", r"\toprule", r"Cell & Part type & placed & requested & rate [95\% CI] \\", r"\midrule"]
    for o, e in K["H5"].items():
        for lo_, v in e["hang"]["H5a"].items():
            T.append(r"%s & %s & %d & %d & %s [%s, %s] \\" % (o, lo_.replace("_", r"\_"), v["dat"], v["yeu_cau"], f3(v["ti_le"]), f3(v["ktc95_CP"][0]), f3(v["ktc95_CP"][1])))
    T += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    open(os.path.join(RA, "bang_cp.tex"), "w").write("\n".join(T) + "\n")

if __name__ == "__main__":
    {"kiem_truoc": buoc_kiem_truoc, "chinh": buoc_chinh}[sys.argv[1]]()
