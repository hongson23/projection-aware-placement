# -*- coding: utf-8 -*-
"""phan_tich_chinh_dv.py — phân tích XÁC NHẬN (ghi trước) lượt chính ĐỘNG VẬT v12.2, Bài 1 chương 3.

Theo đúng: DangKyTruoc_LuotChinh_Bai1_20261111.md (bản 3.1) + DangKyTruoc_BoSung_v121_DongVat_20261114.md (bản 3).
Không thêm kiểm định, không đổi ngưỡng.

BA BƯỚC TÁCH RỜI (chạy tuần tự, mỗi bước ghi tệp riêng):
  1) mu   --ban sau|truoc : CHỈ đọc bảng MÙ (tong.csv của bộ kiểm v3 + mu_<ban>/<id>/thong_tin.json["tu_the"]),
                            tính thước theo TẾ BÀO MÙ → phan_tich/bang_mu_<ban>.csv. Không đọc khoá, không đọc kết quả bộ dựng.
  2) mo                   : bước MỞ MÙ duy nhất: đọc khoa/KHOA_MU.json, gắn cột (ô, nhánh, seed) cho mỗi id; rồi đọc các số
                            BỘ DỰNG tự đo (thong_tin.json gốc) → phan_tich/bang_mo.csv.
  3) thong_ke             : thống kê ghi trước trên bang_mo.csv → phan_tich/ket_qua_so.json + phan_tich/bang_ket_qua.md.

Định nghĩa (trích đăng ký):
  G_c  = Σ_i max(D_i^K − h, 0), h = 0,1797 µm (bộ kiểm do_khuon). Bộ phận MẤT MESH: D_i^K = bán kính bao của khuôn
         (bản 3.1 §4.1) — do_khuon bỏ qua bộ phận mất (KHONG_CO_MESH) nên bước 1 cộng thêm phần này.
  V_test = Σ số ca trượt (phán quyết su_that) của các luật sự thật v1.2, TRỪ luật trùng số hạng mềm 4.5b, 4.6b, 4.8b
         (1.3: vế lệch tâm chỉ nằm ở phán quyết tham_so, phán quyết su_that không đo lệch tâm → giữ 1.3).
         Bộ phận mất mesh: +1 cho mỗi luật áp lên nó (ti thể: 4.3a, 4.3b, M.a, M.c₀, M.h; lysosome: 4.6a, M.a, M.c₀, M.h;
         loại khác: M.a, M.c₀, M.h). V_mem = 4.5b + 4.6b + 4.8b (+1 cho 4.6b mỗi lysosome mất). KKD báo riêng.
  d    = thước(đối chứng) − thước(PLACE) theo seed; H2/H3 một phía: PLACE nhỏ hơn ⇔ d > 0.
  Wilcoxon hạng có dấu CHÍNH XÁC bằng liệt kê 2^{n′} cách gán dấu; |d| < 1e-6 bỏ; hạng trung bình khi trùng.
  Chuỗi R (nn1x) và S (un1x), α = 0,025 mỗi chuỗi, thứ tự D5 → D4A2 → D4; ô "chỉ mô tả" (giờ lệch / bất đồng) bỏ qua, không chặn.
  Hỏng (thiếu việc dựng hoặc thiếu chấm) ở nhánh nào của cặp → cặp tính PLACE THUA (d = −(max|d| + 1)).
  H3: Wilcoxon chính xác trên V_test, Holm 6 kiểm định (2 đối chứng × 3 ô), α = 0,05.
  Bất đồng G_c^K / G_c^B (theo ô, 30 tế bào): Spearman < 0,8 hoặc dấu hiệu theo seed khác nhau ở > 2/10 seed (từng đối chứng).
  G_c^B = Σ_i max(doi_mat_um[i][1] − h, 0) trên cùng tập bộ phận tu_the.
"""
import os, sys, json, csv, argparse, hashlib, itertools, collections, math
import numpy as np

GOC = os.path.expanduser("~/Vision_team_data_189/sonhh20/ch3_tebao_dep")
CC = os.path.join(GOC, "bo_kiem_v3", "cham_chinh_dv")
PT = os.path.join(CC, "phan_tich")
KQ = os.path.join(GOC, "bo_dung_chung", "v121", "kq_chinh_dv")
MAY = os.path.join(GOC, "bo_dung_chung", "v121", "chinh_dv", "hang", "may")
H = 0.1797
TOL0 = 1e-6
LUAT_MEM = {"4.5b", "4.6b", "4.8b"}
LUAT_CUNG_KP = ["1.4a", "4.3a", "M.a", "M.c₀", "M.h"]
LUAT_MAT = {"ti_the": ["4.3a", "4.3b", "M.a", "M.c₀", "M.h"], "lysosome": ["4.6a", "M.a", "M.c₀", "M.h"]}
LUAT_MAT_MAC_DINH = ["M.a", "M.c₀", "M.h"]
LUAT_MAT_MEM = {"lysosome": ["4.6b"]}
O_CHUOI = ["D5", "D4A2", "D4"]
O_TAT = ["D1", "D4", "D4A2", "D5"]
DOI_CHUNG = ["nn1x", "un1x"]
SEED10 = list(range(4001, 4011)); SEED5 = list(range(4001, 4006))


def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def ban_kinh_bao(khuon, s):
    """bán kính bao của khuôn quanh tâm t (x = t + s·R·x_khuôn): s · max_{x ∈ khuôn} |x|."""
    sys.path.insert(0, os.path.join(GOC, "bo_kiem_v3"))
    import luat_v12 as LV
    if khuon == "cau":
        return float(s)
    if khuon in ("nhan", "golgi"):
        a = LV.P("dv.khuon_nhan_truc" if khuon == "nhan" else "dv.khuon_golgi_truc"); return float(s) * float(max(a))
    if khuon.startswith("ong") and khuon[3:].isdigit():
        r, c = LV.P("dv.khuon_ong")[int(khuon[3:])]
        tau = np.linspace(-1, 1, 20001)
        Q = np.stack([tau * (1 - r), c * (1 - tau ** 2), 0 * tau], 1)
        return float(s) * float(np.linalg.norm(Q, axis=1).max() + r)
    raise ValueError("chua co ban kinh bao cho khuon %s" % khuon)


def buoc_mu(ban):
    R = os.path.join(CC, "ds_%s.txt.ket_qua" % ban)
    tong = os.path.join(R, "tong.csv")
    rows = list(csv.DictReader(open(tong, encoding="utf-8")))
    theo = collections.defaultdict(list)
    for r in rows:
        theo[r["te_bao"]].append(r)
    out = []
    for tb, rs in sorted(theo.items()):
        idm = tb.split("__", 1)[1]
        tu_the = json.load(open(os.path.join(CC, "mu_%s" % ban, idm, "thong_tin.json")))["tu_the"]
        g = [r for r in rs if r["luat"] == "G_c"]
        assert len(g) == 1, tb
        g = g[0]
        Gc_goc = float(g["trang_thai"]); n_do = int(g["so_ca"]); G_dem_goc = int(g["so_ca_truot"])
        Di = [r for r in rs if r["luat"] == "D_i"]
        mat = [(r["phan_quyet"], r["bo_phan_truot"]) for r in Di if r["trang_thai"] == "KHONG_CO_MESH"]
        kdl = [r["phan_quyet"] for r in Di if r["trang_thai"] == "KHONG_DUNG_LAI_DUOC"]
        Dv = [float(r["trang_thai"]) for r in Di if r["trang_thai"] not in ("KHONG_CO_MESH", "KHONG_DUNG_LAI_DUOC", "")]
        assert len(Dv) == n_do, (tb, len(Dv), n_do)
        assert abs(sum(max(x - H, 0) for x in Dv) - Gc_goc) < 1e-3, tb
        R_mat = []
        for ten, kh in mat:
            assert tu_the[ten]["khuon"] == kh
            R_mat.append(ban_kinh_bao(kh, tu_the[ten]["s"]))
        Gc_K = Gc_goc + sum(max(x - H, 0) for x in R_mat)
        G_dem = G_dem_goc + sum(1 for x in R_mat if x > H)
        st = {r["luat"]: r for r in rs if r["phan_quyet"] == "su_that" and r["loai_luat"] == "su_that_v1.2"}
        tr = {k: int(v["so_ca_truot"]) for k, v in st.items()}
        V_goc = sum(v for k, v in tr.items() if k not in LUAT_MEM)
        Vm_goc = sum(v for k, v in tr.items() if k in LUAT_MEM)
        pen = 0; pen_m = 0
        for ten, kh in mat:
            lo = ten.rstrip("0123456789").rstrip("_")
            pen += len(LUAT_MAT.get(lo, LUAT_MAT_MAC_DINH)); pen_m += len(LUAT_MAT_MEM.get(lo, []))
        kk = [r for r in rs if r["loai_luat"] == "dem_kkd_nhan"]
        kkd_nhan = int(kk[0]["so_ca"]) if kk else 0
        kkd_luat = sorted(k for k, v in st.items() if v["trang_thai"] == "KHONG_KIEM_DUOC")
        mc = st["M.c₀"]
        out.append({
            "id": idm, "ban": ban, "Gc_K": round(Gc_K, 5), "Gc_kiem_goc": Gc_goc, "n_bo_phan_do": n_do, "G_dem": G_dem,
            "D_max": round(max(Dv), 5) if Dv else "", "n_D_tren_h2": sum(1 for x in Dv if x > H / 2),
            "n_mat": len(mat), "mat": ";".join(t for t, _ in mat), "R_mat_tong": round(sum(R_mat), 5), "n_khong_dung_lai": len(kdl),
            "V_test": V_goc + pen, "V_test_chua_phat_mat": V_goc, "V_mem": Vm_goc + pen_m,
            "V_test_kkd_la_vp": V_goc + pen + kkd_nhan,
            "t_4.3b": tr.get("4.3b", 0), "t_4.5a": tr.get("4.5a", 0), "t_Mc0": tr.get("M.c₀", 0),
            "t_4.5b": tr.get("4.5b", 0), "t_4.6b": tr.get("4.6b", 0), "t_4.8b": tr.get("4.8b", 0),
            "Mc0_trang_thai": mc["trang_thai"], "Mc0_cap": mc["bo_phan_truot"], "kkd_nhan_Mc0": kkd_nhan,
            "n_luat_kkd": len(kkd_luat), "luat_kkd": ";".join(kkd_luat),
            "trang_thai_4.4a": st["4.4a"]["trang_thai"],
            "kp_c_truot_bo_kiem": sum(tr.get(k, 0) for k in LUAT_CUNG_KP),
        })
    os.makedirs(PT, exist_ok=True)
    f = os.path.join(PT, "bang_mu_%s.csv" % ban)
    with open(f, "w", newline="", encoding="utf-8") as fo:
        w = csv.DictWriter(fo, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
    json.dump({"tong_csv": tong, "md5_tong_csv": md5(tong), "so_te_bao": len(out), "het": os.path.exists(os.path.join(R, "HET")),
               "so_loi": len(os.listdir(os.path.join(R, "loi"))), "md5_bang_mu": md5(f), "md5_ma": md5(os.path.abspath(__file__))},
              open(os.path.join(PT, "bang_mu_%s.meta.json" % ban), "w"), indent=1)
    print("bang mu %s: %d te bao -> %s" % (ban, len(out), f))


def tach_ten(job):
    p = job.split("_")
    return p[0], "_".join(p[1:-1]), int(p[-1])


def buoc_mo():
    KH = os.path.join(CC, "khoa", "KHOA_MU.json")
    khoa = json.load(open(KH))
    hang = []
    for ban in ("sau", "truoc"):
        f = os.path.join(PT, "bang_mu_%s.csv" % ban)
        if not os.path.exists(f):
            continue
        for r in csv.DictReader(open(f, encoding="utf-8")):
            k = khoa["id"][r["id"]]
            assert k["ban"] == ban, (r["id"], k["ban"], ban)
            job = os.path.basename(k["nguon"].rstrip("/"))
            o, nh, sd = tach_ten(job)
            tt = json.load(open(os.path.join(KQ, job, "thong_tin.json")))
            tu_mu = json.load(open(os.path.join(CC, "mu_%s" % ban, r["id"], "thong_tin.json")))["tu_the"]
            r.update({"job": job, "o": o, "nhanh": nh, "seed": sd, "tu_the_khop": int(tu_mu == tt["tu_the"])})
            bc = tt["bao_cao_bo_tri"]; doi = tt.get("doi_mat_um") or {}
            gB = 0.0; nB_none = 0
            for ten in tt["tu_the"]:
                v = doi.get(ten)
                if not v or v[1] is None:
                    nB_none += 1; continue
                gB += max(v[1] - H, 0.0)
            may = os.path.join(MAY, job + ".txt"); tai_cao = ""
            if os.path.exists(may):
                s_ = open(may).read(); tai_cao = int("tai_cao 1" in s_)
            r.update({"giay_dat": bc["giay_dat"], "danh_gia_phat": bc.get("danh_gia_phat"), "J": bc["J_truoc_chieu"]["J"],
                      "Gc_B": round(gB, 5), "Gc_B_khong_co": nB_none,
                      "goT_B": sum(1 for kk, v in doi.items() if v and v[1] is not None and v[1] > tt["h1"] and not kk.startswith("luoi_noi_chat")),
                      "mat_B": ";".join(tt.get("mat_bo_phan") or []), "n_mat_B": len(tt.get("mat_bo_phan") or []),
                      "nut_vo_nghiem": tt.get("diem_khong_kha_thi_truoc_chieu"), "con_vp_max_um": tt.get("con_vp_max_um"),
                      "diem_vp_con_lai": tt.get("diem_vi_pham_con_lai"), "khong_hoi_tu_nut": tt["thong_ke_chieu"].get("khong_hoi_tu_nut"),
                      "so_nut_tong": sum(x.get("so_nut", 0) for x in (tt.get("diem_chieu_theo_luot") or []) if isinstance(x, dict)),
                      "ban_dung": tt.get("ban"), "sdf_md5": (tt.get("ban_v121") or {}).get("sdf_md5"), "tai_cao": tai_cao})
            hang.append(r)
    f = os.path.join(PT, "bang_mo.csv")
    keys = list(hang[0].keys())
    with open(f, "w", newline="", encoding="utf-8") as fo:
        w = csv.DictWriter(fo, fieldnames=keys); w.writeheader(); w.writerows(hang)
    jobs = set(os.listdir(KQ)); thieu = sorted(j for j in jobs if j not in set(h["job"] for h in hang) and j not in set(khoa["bo_qua"]))
    json.dump({"md5_khoa": md5(KH), "bo_qua_khoa": khoa["bo_qua"], "thieu_khoa": khoa["thieu"], "job_dung_khong_co_cham": thieu,
               "md5_bang_mo": md5(f), "md5_ma": md5(os.path.abspath(__file__))}, open(os.path.join(PT, "bang_mo.meta.json"), "w"), indent=1)
    print("bang mo: %d hang -> %s; job thieu cham: %s" % (len(hang), f, thieu))


def wilcoxon_chinh_xac(d, tol=TOL0):
    """một phía H1: d > 0. Liệt kê đủ 2^{n'} cách gán dấu. Trả (p, W+, n', số hoà)."""
    d = np.asarray(d, float)
    nz = d[np.abs(d) >= tol]; n = len(nz)
    if n == 0:
        return 1.0, 0.0, 0, len(d)
    a = np.round(np.abs(nz), 9)
    order = np.argsort(a, kind="mergesort"); rk = np.empty(n)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and a[order[j + 1]] == a[order[i]]:
            j += 1
        rk[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    Wp = float(rk[nz > 0].sum())
    dem = 0
    for dau in itertools.product((0, 1), repeat=n):
        if float(np.dot(dau, rk)) >= Wp - 1e-9:
            dem += 1
    return dem / 2.0 ** n, Wp, n, len(d) - n


def dau_chinh_xac(thang, thua):
    n = thang + thua
    if n == 0:
        return 1.0
    return sum(math.comb(n, k) for k in range(thang, n + 1)) / 2.0 ** n


def tv(x):
    x = [v for v in x if v is not None and not (isinstance(v, float) and math.isnan(v))]
    if not x:
        return None
    q = np.percentile(x, [25, 50, 75])
    return [float(q[1]), float(q[0]), float(q[2])]


def spearman(a, b):
    from scipy.stats import spearmanr
    a = np.asarray(a, float); b = np.asarray(b, float)
    if np.ptp(a) == 0 or np.ptp(b) == 0:
        return float("nan")
    return float(spearmanr(a, b).correlation)


def buoc_thong_ke():
    rows = list(csv.DictReader(open(os.path.join(PT, "bang_mo.csv"), encoding="utf-8")))
    num = ["Gc_K", "Gc_kiem_goc", "G_dem", "D_max", "n_D_tren_h2", "n_mat", "V_test", "V_test_chua_phat_mat", "V_mem", "V_test_kkd_la_vp",
           "t_4.3b", "t_4.5a", "t_Mc0", "kkd_nhan_Mc0", "n_luat_kkd", "kp_c_truot_bo_kiem", "giay_dat", "J", "Gc_B", "goT_B", "n_mat_B",
           "nut_vo_nghiem", "con_vp_max_um", "diem_vp_con_lai", "seed", "tu_the_khop", "n_khong_dung_lai"]
    for r in rows:
        for k in num:
            if k in r and r[k] not in ("", None, "None"):
                r[k] = float(r[k])
            elif k in r:
                r[k] = None
    B = {(r["ban"], r["o"], r["nhanh"], int(r["seed"])): r for r in rows}
    KQS = {"ghi_chu": "d = doi_chung - PLACE; H1: d > 0 (PLACE nho hon)"}
    md = []

    def lay(ban, o, nh, sd):
        return B.get((ban, o, nh, sd))

    def cap(ban, o, nh, thuoc, seeds, goc="td"):
        """→ d (đã xử lý hỏng), số cặp hỏng, thắng/hoà/thua, trung vị d."""
        ds = []; hong = 0
        for sd in seeds:
            a = lay(ban, o, goc, sd); b = lay(ban, o, nh, sd)
            if a is None or b is None or a.get(thuoc) is None or b.get(thuoc) is None:
                ds.append(None); hong += 1
            else:
                ds.append(b[thuoc] - a[thuoc])
        ok = [x for x in ds if x is not None]
        m = max([abs(x) for x in ok] + [0.0]) + 1.0
        ds = [x if x is not None else -m for x in ds]
        th = sum(1 for x in ds if x >= TOL0); hoa = sum(1 for x in ds if abs(x) < TOL0); thua = sum(1 for x in ds if x <= -TOL0)
        return ds, hong, th, hoa, thua, float(np.median(ds))

    def ti_gio(o, nh, seeds, goc="td"):
        t = [lay("sau", o, nh, s)["giay_dat"] / lay("sau", o, goc, s)["giay_dat"] for s in seeds if lay("sau", o, nh, s) and lay("sau", o, goc, s)]
        return float(np.median(t)) if t else None, t

    A = {}
    for ban in ("sau", "truoc"):
        for (b_, o, nh, sd), r in B.items():
            if b_ != ban:
                continue
            A.setdefault((ban, o, nh), []).append(r)
    bangA = []
    for (ban, o, nh), rs in sorted(A.items()):
        bangA.append({"ban": ban, "o": o, "nhanh": nh, "n": len(rs),
                      "Gc_K": tv([r["Gc_K"] for r in rs]), "V_test": tv([r["V_test"] for r in rs]), "V_mem": tv([r["V_mem"] for r in rs]),
                      "G_dem": tv([r["G_dem"] for r in rs]),
                      "Mc0_te_bao_truot": sum(1 for r in rs if r["t_Mc0"] and r["t_Mc0"] > 0), "Mc0_cap_truot": int(sum(r["t_Mc0"] or 0 for r in rs)),
                      "Mc0_cap_loai": sorted(set(r["Mc0_cap"] for r in rs if r["Mc0_cap"])),
                      "kkd_nhan_te_bao": sum(1 for r in rs if (r["kkd_nhan_Mc0"] or 0) > 0), "n_luat_kkd": sorted(set(r["luat_kkd"] for r in rs)),
                      "t_4.3b": tv([r["t_4.3b"] for r in rs]), "t_4.5a_tong": int(sum(r["t_4.5a"] or 0 for r in rs)),
                      "n_mat_tong": int(sum(r["n_mat"] for r in rs)), "te_bao_mat": sum(1 for r in rs if r["n_mat"] > 0),
                      "kp_c_dat_bo_kiem": "%d/%d" % (sum(1 for r in rs if r["kp_c_truot_bo_kiem"] == 0), len(rs)),
                      "kp_c_dat_kem_mat": "%d/%d" % (sum(1 for r in rs if r["kp_c_truot_bo_kiem"] == 0 and r["n_mat"] == 0), len(rs)),
                      "kp_a_0_nut_vo_nghiem": "%d/%d" % (sum(1 for r in rs if (r["nut_vo_nghiem"] or 0) == 0), len(rs)),
                      "nut_vo_nghiem": tv([r["nut_vo_nghiem"] for r in rs]), "con_vp_max_um": max((r["con_vp_max_um"] or 0) for r in rs),
                      "T_giay_dat": tv([r["giay_dat"] for r in rs]), "J": tv([r["J"] for r in rs]),
                      "Gc_B": tv([r["Gc_B"] for r in rs]), "goT_B": tv([r["goT_B"] for r in rs]),
                      "n_D_tren_h2": tv([r["n_D_tren_h2"] for r in rs]), "D_max": tv([r["D_max"] for r in rs]),
                      "tu_the_khop": "%d/%d" % (sum(1 for r in rs if r["tu_the_khop"] == 1), len(rs)),
                      "tai_cao": sum(1 for r in rs if str(r.get("tai_cao")) in ("1", "1.0"))})
    KQS["A_mo_ta"] = bangA

    BD = {}
    for o in O_TAT:
        rs = [lay("sau", o, nh, s) for nh in ["td"] + DOI_CHUNG for s in SEED10]
        rs = [r for r in rs if r]
        rho = spearman([r["Gc_K"] for r in rs], [r["Gc_B"] for r in rs])
        lech = {}
        for nh in DOI_CHUNG:
            n_ = 0; ds_ = []
            for s in SEED10:
                a, b = lay("sau", o, "td", s), lay("sau", o, nh, s)
                if not (a and b):
                    continue
                dK = b["Gc_K"] - a["Gc_K"]; dB = b["Gc_B"] - a["Gc_B"]
                sK = 0 if abs(dK) < TOL0 else int(np.sign(dK)); sB = 0 if abs(dB) < TOL0 else int(np.sign(dB))
                if sK != sB:
                    n_ += 1; ds_.append([s, round(dK, 5), round(dB, 5)])
            lech[nh] = {"so_seed_lech_dau": n_, "seed": ds_}
        bat = (not (rho >= 0.8)) or any(v["so_seed_lech_dau"] > 2 for v in lech.values())
        BD[o] = {"n": len(rs), "spearman": rho, "lech_dau": lech, "bat_dong": bool(bat)}
    KQS["B_bat_dong"] = BD

    C = {}
    for nh in DOI_CHUNG:
        dung = False
        for o in O_CHUOI + ["D1"]:
            ds, hong, th, hoa, thua, mdd = cap("sau", o, nh, "Gc_K", SEED10)
            p, Wp, n1, z = wilcoxon_chinh_xac(ds)
            tg, tl = ti_gio(o, nh, SEED10)
            chi_mo_ta = []
            if o == "D1":
                chi_mo_ta.append("D1 chi mo ta (ban 3.1)")
            if tg is None or not (0.8 <= tg <= 1.25):
                chi_mo_ta.append("ti so gio ngoai [0,8; 1,25]")
            if BD[o]["bat_dong"]:
                chi_mo_ta.append("bat dong G_c^K/G_c^B")
            if o == "D1":
                trang = "mo_ta"
            elif chi_mo_ta:
                trang = "bo_qua_chi_mo_ta"
            elif dung:
                trang = "khong_kiem (chuoi da dung)"
            else:
                trang = "BAC_BO_H0" if p <= 0.025 else "KHONG_BAC_BO"
                if p > 0.025:
                    dung = True
            C[(nh, o)] = {"doi_chung": nh, "o": o, "d": [round(x, 5) for x in ds], "hong": hong, "thang": th, "hoa": hoa, "thua": thua,
                          "trung_vi_d": mdd, "p": p if o != "D1" else None, "p_tham_khao": p, "W+": Wp, "n_khac_0": n1, "ti_gio_tv": tg,
                          "chi_mo_ta": chi_mo_ta, "ket_qua": trang}
    KQS["C_chuoi_H2"] = list(C.values())

    H3 = []
    for nh in DOI_CHUNG:
        for o in O_CHUOI:
            ds, hong, th, hoa, thua, mdd = cap("sau", o, nh, "V_test", SEED10)
            p, Wp, n1, z = wilcoxon_chinh_xac(ds)
            ds2 = cap("sau", o, nh, "V_test_kkd_la_vp", SEED10)
            p2 = wilcoxon_chinh_xac(ds2[0])[0]
            H3.append({"doi_chung": nh, "o": o, "d": [round(x, 5) for x in ds], "thang": th, "hoa": hoa, "thua": thua, "trung_vi_d": mdd,
                       "p": p, "n_khac_0": n1, "chi_mo_ta": C[(nh, o)]["chi_mo_ta"] if "ti so gio" in " ".join(C[(nh, o)]["chi_mo_ta"]) else [],
                       "p_do_nhay_kkd_la_vp": p2})
    ks = sorted([h for h in H3 if not h["chi_mo_ta"]], key=lambda h: h["p"]); m = len(ks); dung = False
    for i, h in enumerate(ks):
        nguong = 0.05 / (m - i); h["holm_nguong"] = nguong
        if not dung and h["p"] <= nguong:
            h["holm"] = "BAC_BO_H0"
        else:
            dung = True; h["holm"] = "KHONG_BAC_BO"
    KQS["D_H3"] = H3

    GD = []
    for nh in DOI_CHUNG:
        for o in O_TAT:
            ds, hong, th, hoa, thua, mdd = cap("sau", o, nh, "G_dem", SEED10)
            GD.append({"doi_chung": nh, "o": o, "thang": th, "hoa": hoa, "thua": thua, "trung_vi_d": mdd, "p_dau_mo_ta": dau_chinh_xac(th, thua)})
    KQS["E_G_dem"] = GD

    KQS["F_J_PLACE_vs_SA"] = [dict(zip(["o", "thang", "hoa", "thua", "trung_vi_d"], [o] + list(cap("sau", o, "un1x", "J", SEED10)[2:]))) for o in O_TAT]

    QT = []
    for o in O_TAT:
        for nh in ["td", "nn1x", "nn4x", "un1x", "un4x"]:
            rs = [lay("sau", o, nh, s) for s in SEED5]; rs = [r for r in rs if r]
            tg = ti_gio(o, nh, SEED5)[0] if nh != "td" else 1.0
            e = {"o": o, "nhanh": nh, "n": len(rs), "T_tv": tv([r["giay_dat"] for r in rs]), "ti_gio_tv": tg,
                 "Gc_K": tv([r["Gc_K"] for r in rs]), "V_test": tv([r["V_test"] for r in rs])}
            if nh != "td":
                for th_ in ("Gc_K", "V_test"):
                    x = cap("sau", o, nh, th_, SEED5)
                    e["vs_PLACE_" + th_] = {"thang_PLACE": x[2], "hoa": x[3], "thua_PLACE": x[4], "trung_vi_d": x[5]}
            QT.append(e)
    KQS["G_chat_luong_thoi_gian"] = QT

    PH = []
    for o in ("D4A2", "D5"):
        for nh in ("deu", "km"):
            for ban in ("sau", "truoc"):
                rs = [lay(ban, o, nh, s) for s in SEED10]; rs = [r for r in rs if r]
                if not rs:
                    continue
                e = {"ban": ban, "o": o, "nhanh": nh, "n": len(rs), "Gc_K": tv([r["Gc_K"] for r in rs]), "V_test": tv([r["V_test"] for r in rs]),
                     "J": tv([r["J"] for r in rs]), "Mc0_te_bao_truot": sum(1 for r in rs if (r["t_Mc0"] or 0) > 0),
                     "ti_gio_tv": ti_gio(o, nh, SEED10)[0]}
                for th_ in ("Gc_K", "V_test", "J"):
                    x = cap(ban, o, nh, th_, SEED10)
                    e["vs_PLACE_" + th_] = {"thang_PLACE": x[2], "hoa": x[3], "thua_PLACE": x[4], "trung_vi_d": x[5]}
                PH.append(e)
    KQS["H_doi_chung_phu"] = PH
    KHE = []
    for o in ("D4A2", "D5"):
        for khe in ("khe0.1_", "", "khe0.5_"):
            for nh in ("td", "nn1x", "un1x"):
                rs = [lay("sau", o, khe + nh, s) for s in SEED5]; rs = [r for r in rs if r]
                e = {"o": o, "khe": khe.rstrip("_") or "khe0.3 (chinh)", "nhanh": nh, "n": len(rs), "Gc_K": tv([r["Gc_K"] for r in rs]),
                     "V_test": tv([r["V_test"] for r in rs]), "Mc0_te_bao_truot": sum(1 for r in rs if (r["t_Mc0"] or 0) > 0)}
                if nh != "td":
                    e["ti_gio_tv"] = ti_gio(o, khe + nh, SEED5, goc=khe + "td")[0]
                    for th_ in ("Gc_K", "V_test"):
                        x = cap("sau", o, khe + nh, th_, SEED5, goc=khe + "td")
                        e["vs_PLACE_" + th_] = {"thang_PLACE": x[2], "hoa": x[3], "thua_PLACE": x[4], "trung_vi_d": x[5]}
                KHE.append(e)
    KQS["H_khe"] = KHE

    BDS = []
    for r in rows:
        if r["ban"] != "sau":
            continue
        ly = []
        if (r["t_Mc0"] or 0) > 0 and (r["diem_vp_con_lai"] or 0) == 0:
            ly.append("M.c0 truot (%s, %d cap) du bo dung bao 0 vi pham con lai" % (r["Mc0_cap"], r["t_Mc0"]))
        if set(filter(None, r["mat"].split(";"))) != set(filter(None, (r["mat_B"] or "").split(";"))):
            ly.append("mat bo phan: bo kiem %d vs bo dung %d" % (r["n_mat"], r["n_mat_B"]))
        if (r["kp_c_truot_bo_kiem"] or 0) > 0 and (r["diem_vp_con_lai"] or 0) == 0:
            ly.append("K-PROJECT(c) truot o bo kiem, bo dung bao 0 vi pham con lai")
        if ly:
            BDS.append({"job": r["job"], "ly_do": ly})
    KQS["I_bat_dong_bo_dung"] = BDS
    KQS["I_got"] = [{"o": e["o"], "nhanh": e["nhanh"], "goT_B": e["goT_B"], "G_dem_kiem": e["G_dem"]} for e in bangA if e["ban"] == "sau"]

    KQS["J_tai_cao_so_viec"] = sum(1 for r in rows if str(r.get("tai_cao")) in ("1", "1.0"))
    KQS["kiem_toan_ven"] = {"so_hang": len(rows), "tu_the_khong_khop": [r["job"] for r in rows if r["tu_the_khop"] != 1],
                            "ban_dung_khac": sorted(set(r["ban_dung"] for r in rows)), "sdf_md5": sorted(set(r["sdf_md5"] for r in rows)),
                            "n_khong_dung_lai": int(sum(r["n_khong_dung_lai"] or 0 for r in rows))}
    json.dump(KQS, open(os.path.join(PT, "ket_qua_so.json"), "w"), ensure_ascii=False, indent=1, default=str)

    def f3(x, k=3):
        if x is None:
            return "—"
        return "%.*g [%.*g; %.*g]" % (k, x[0], k, x[1], k, x[2])
    L = ["# Bảng số (sinh tự động bởi phan_tich_chinh_dv.py, md5 %s)\n" % md5(os.path.abspath(__file__))]
    L.append("## A. Mô tả theo ô × nhánh (trung vị [Q1; Q3])\n")
    L.append("| bản | ô | nhánh | n | G_c^K (µm) | V_test | V_mem | G đếm | M.c₀ tế bào trượt (cặp) | KKD nhân M.c₀ | luật KKD | 4.3b trượt | mất bộ phận (tế bào) | K-PROJECT(c) bộ kiểm | (c) kể mất | (a) 0 nút vô nghiệm | T đặt (s) | G_c^B | gọt bộ dựng |")
    L.append("|" + "---|" * 19)
    for e in bangA:
        L.append("| %s | %s | %s | %d | %s | %s | %s | %s | %d (%d) | %d | %s | %s | %d (%d) | %s | %s | %s | %s | %s | %s |" % (
            e["ban"], e["o"], e["nhanh"], e["n"], f3(e["Gc_K"]), f3(e["V_test"]), f3(e["V_mem"]), f3(e["G_dem"]), e["Mc0_te_bao_truot"], e["Mc0_cap_truot"],
            e["kkd_nhan_te_bao"], "/".join(e["n_luat_kkd"]), f3(e["t_4.3b"]), e["n_mat_tong"], e["te_bao_mat"], e["kp_c_dat_bo_kiem"], e["kp_c_dat_kem_mat"],
            e["kp_a_0_nut_vo_nghiem"], f3(e["T_giay_dat"], 4), f3(e["Gc_B"]), f3(e["goT_B"])))
    L.append("\n## B. Bất đồng G_c^K / G_c^B\n")
    for o, v in BD.items():
        L.append("- %s: Spearman %.3f; lệch dấu nn1x %d, un1x %d → %s" % (o, v["spearman"], v["lech_dau"]["nn1x"]["so_seed_lech_dau"], v["lech_dau"]["un1x"]["so_seed_lech_dau"], "BẤT ĐỒNG" if v["bat_dong"] else "đồng thuận"))
    L.append("\n## C. Chuỗi H2 (G_c), α = 0,025 mỗi chuỗi\n")
    L.append("| chuỗi | ô | thắng/hoà/thua PLACE | n′ | W+ | p chính xác | trung vị d (µm) | tỉ số giờ | hỏng | ghi | kết quả |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for e in KQS["C_chuoi_H2"]:
        L.append("| %s | %s | %d/%d/%d | %d | %g | %s | %.4g | %.3f | %d | %s | %s |" % (e["doi_chung"], e["o"], e["thang"], e["hoa"], e["thua"], e["n_khac_0"], e["W+"],
                 ("%.5f" % e["p"]) if e["p"] is not None else "(%.5f, mô tả)" % e["p_tham_khao"], e["trung_vi_d"], e["ti_gio_tv"] or float("nan"), e["hong"], "; ".join(e["chi_mo_ta"]), e["ket_qua"]))
    L.append("\n## D. H3 (V_test), Holm 6 kiểm định α = 0,05\n")
    L.append("| đối chứng | ô | thắng/hoà/thua | n′ | p | ngưỡng Holm | Holm | trung vị d | p (KKD là vi phạm, mô tả) |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for e in H3:
        L.append("| %s | %s | %d/%d/%d | %d | %.5f | %s | %s | %g | %.5f |" % (e["doi_chung"], e["o"], e["thang"], e["hoa"], e["thua"], e["n_khac_0"], e["p"],
                 ("%.4f" % e["holm_nguong"]) if "holm_nguong" in e else "—", e.get("holm", "mô tả"), e["trung_vi_d"], e["p_do_nhay_kkd_la_vp"]))
    L.append("\n## E. G (đếm), phép thử dấu chính xác (mô tả)\n")
    for e in GD:
        L.append("- %s %s: %d/%d/%d, p dấu %.4f, trung vị d %g" % (e["doi_chung"], e["o"], e["thang"], e["hoa"], e["thua"], e["p_dau_mo_ta"], e["trung_vi_d"]))
    L.append("\n## F. J trước chiếu PLACE vs SA (mô tả)\n")
    for e in KQS["F_J_PLACE_vs_SA"]:
        L.append("- %s: PLACE thắng/hoà/thua %d/%d/%d, trung vị d %.4g" % (e["o"], e["thang"], e["hoa"], e["thua"], e["trung_vi_d"]))
    L.append("\n## G. Chất lượng–thời gian (seed 4001–4005)\n")
    L.append("| ô | nhánh | T đặt (s) | tỉ số giờ | G_c^K | V_test | PLACE vs nhánh G_c (thắng/hoà/thua; tv d) | V_test |")
    L.append("|---|---|---|---|---|---|---|---|")
    for e in QT:
        a = e.get("vs_PLACE_Gc_K"); b = e.get("vs_PLACE_V_test")
        L.append("| %s | %s | %s | %s | %s | %s | %s | %s |" % (e["o"], e["nhanh"], f3(e["T_tv"], 4), "%.3f" % e["ti_gio_tv"] if e["ti_gio_tv"] else "—", f3(e["Gc_K"]), f3(e["V_test"]),
                 ("%d/%d/%d; %.4g" % (a["thang_PLACE"], a["hoa"], a["thua_PLACE"], a["trung_vi_d"])) if a else "—",
                 ("%d/%d/%d; %.4g" % (b["thang_PLACE"], b["hoa"], b["thua_PLACE"], b["trung_vi_d"])) if b else "—"))
    L.append("\n## H. Đối chứng phụ (deu, km; 10 seed) và độ nhạy khe (5 seed) — mô tả\n")
    L.append("| bản | ô | nhánh | G_c^K | V_test | J | M.c₀ tb trượt | tỉ số giờ | PLACE vs: G_c | V_test | J |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for e in PH:
        z = [e["vs_PLACE_" + k] for k in ("Gc_K", "V_test", "J")]
        L.append("| %s | %s | %s | %s | %s | %s | %d | %.3f | %s |" % (e["ban"], e["o"], e["nhanh"], f3(e["Gc_K"]), f3(e["V_test"]), f3(e["J"]), e["Mc0_te_bao_truot"], e["ti_gio_tv"] or float("nan"),
                 " | ".join("%d/%d/%d; %.4g" % (x["thang_PLACE"], x["hoa"], x["thua_PLACE"], x["trung_vi_d"]) for x in z)))
    L.append("\n| ô | khe | nhánh | G_c^K | V_test | M.c₀ tb trượt | tỉ số giờ | PLACE vs: G_c | V_test |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for e in KHE:
        a = e.get("vs_PLACE_Gc_K"); b = e.get("vs_PLACE_V_test")
        L.append("| %s | %s | %s | %s | %s | %d | %s | %s | %s |" % (e["o"], e["khe"], e["nhanh"], f3(e["Gc_K"]), f3(e["V_test"]), e["Mc0_te_bao_truot"], ("%.3f" % e["ti_gio_tv"]) if e.get("ti_gio_tv") else "—",
                 ("%d/%d/%d; %.4g" % (a["thang_PLACE"], a["hoa"], a["thua_PLACE"], a["trung_vi_d"])) if a else "—",
                 ("%d/%d/%d; %.4g" % (b["thang_PLACE"], b["hoa"], b["thua_PLACE"], b["trung_vi_d"])) if b else "—"))
    L.append("\n## I. Bất đồng bộ kiểm / bộ dựng tự đo (sau chiếu)\n")
    dem = collections.Counter()
    for e in BDS:
        for l_ in e["ly_do"]:
            dem[(e["job"].rsplit("_", 1)[0], l_.split(" (")[0].split(":")[0])] += 1
    for (k, l_), n in sorted(dem.items()):
        L.append("- %s: %s × %d" % (k, l_, n))
    L.append("\nToàn vẹn: %s; tải cao: %d việc" % (json.dumps(KQS["kiem_toan_ven"], ensure_ascii=False), KQS["J_tai_cao_so_viec"]))
    open(os.path.join(PT, "bang_ket_qua.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("buoc", choices=["mu", "mo", "thong_ke"]); ap.add_argument("--ban", default="sau")
    a = ap.parse_args()
    {"mu": lambda: buoc_mu(a.ban), "mo": buoc_mo, "thong_ke": buoc_thong_ke}[a.buoc]()
