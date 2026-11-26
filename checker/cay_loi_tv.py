# -*- coding: utf-8 -*-
"""
cay_loi_tv.py — kiểm thử đột biến cho bộ kiểm tế bào THỰC VẬT (kiem_thucvat.py).

Gốc: tế bào thử đúng-theo-cách-dựng TV1 (Sinh10), TV2 (Sinh10, có ribosome), TV6 (KHTN6) do
dung_tv_thu.py dựng. Mỗi đột biến cấy MỘT lỗi; "phát hiện" ⇔ luật ĐÍCH có ca TRƯỢT MỚI (khoá không
chứa giá trị) nhắc đúng bộ phận đích, và bộ phận đó ĐẠT luật đích ở gốc (nếu không: bỏ, ghi không hợp lệ).
Đối chứng: 'sat_nguong_dat' (luật đích không được trượt mới tại bộ phận) và 'trung_tinh' (không luật nào
được trượt mới).  Chạy: python cay_loi_tv.py   (ghi dot_bien_tv/, cay_loi_tv_ket_qua.json)
"""
import os, sys, json, copy, time, shutil, collections
import numpy as np
import trimesh
import manifold3d as mf
import kiem2, kiem_thucvat
import kiem3 as _K3M, v3_cay_loi_chung as V3, luat_v12 as LV


def _K3TV(d, im_lang=True):
    return _K3M.kiem3(d, "thuc_vat", im_lang=True)
import dung_tv_thu as D
import cay_loi as CA

HERE = os.path.dirname(os.path.abspath(__file__))
GOC = os.path.join(HERE, "tv_thu")
OUT = os.path.join(HERE, "dot_bien_tv")


def doc(tb):
    z = np.load(os.path.join(GOC, tb, "mesh.npz")); L = json.load(open(os.path.join(GOC, tb, "loai.json")))
    M = {}
    for k in z.files:
        n, s = k.rsplit("__", 1); M.setdefault(n, {})[s] = np.array(z[k])
    return M, L


def ghi(M, L, d):
    os.makedirs(d, exist_ok=True)
    np.savez_compressed(os.path.join(d, "mesh.npz"), **{n + "__" + s: a for n, x in M.items() for s, a in x.items()})
    json.dump(L, open(os.path.join(d, "loai.json"), "w"), ensure_ascii=False)


def tl(L, lo):
    return sorted([n for n, t in L["loai"].items() if t == lo], key=lambda s: (len(s), s))


def con_chau(L, ten):
    out = [ten]; i = 0
    while i < len(out):
        out += [c for c, p in L["cha"].items() if p == out[i]]; i += 1
    return out


def man(M, n):
    return mf.Manifold(mf.Mesh64(np.ascontiguousarray(M[n]["v"], float), np.ascontiguousarray(M[n]["f"], np.uint64)))


def dat_m(M, n, m):
    V, F = D.mesh_of(m.simplify(1e-6)); M[n]["v"] = V; M[n]["f"] = F


def dich_cum(M, L, ten, t):
    for n in con_chau(L, ten):
        if n in M: M[n]["v"] = M[n]["v"] + np.asarray(t)[None]


def bien_doi_cum(M, L, ten, fn):
    for n in con_chau(L, ten):
        if n in M: M[n]["v"] = fn(M[n]["v"])


def xoa(M, L, ten, du_kien=True, cum=True):
    for n in (con_chau(L, ten) if cum else [ten]):
        M.pop(n, None); L["loai"].pop(n, None); L["cha"].pop(n, None)
        if du_kien: L["bo_phan_du_kien"].pop(n, None)


def them(M, L, ten, m, lo, cha=None):
    V, F = D.mesh_of(m); M[ten] = {"v": V, "f": F}; L["loai"][ten] = lo; L["bo_phan_du_kien"][ten] = lo
    if cha: L["cha"][ten] = cha


def huong_mang(M, n):
    """Pháp tuyến ra ngoài của mặt hộp gần tâm bộ phận nhất (±x, ±y, ±z)."""
    c = M[n]["v"].mean(0); g = np.array([30 - abs(c[0]), 15 - abs(c[1]), 15 - abs(c[2])]); k = int(np.argmin(g))
    u = np.zeros(3); u[k] = np.sign(c[k]); return u


def L_(M, n):
    return kiem2.Luoi(M[n]["v"], M[n]["f"])


def kc(M, a, b):
    return float(man(M, a).min_gap(man(M, b), 20.0))


def khe_toi(M, L, a, b, khe, cum=True):
    """Dịch cụm a theo hướng tâm a → tâm b cho tới khi khe chính xác a–b = khe."""
    for _ in range(12):
        d = kc(M, a, b)
        if abs(d - khe) < 2e-4: break
        u = M[b]["v"].mean(0) - M[a]["v"].mean(0); u /= np.linalg.norm(u)
        dich_cum(M, L, a, u * (d - khe)) if cum else None
    return "khe %s–%s = %.3f (do lai %.4f)" % (a, b, khe, kc(M, a, b))


def khe_mang(M, L, a, khe):
    """Dịch cụm a dọc pháp tuyến mặt hộp gần nhất để khe tới màng = khe."""
    u = huong_mang(M, a)
    for _ in range(12):
        d = kc(M, a, "mang")
        if abs(d - khe) < 2e-4: break
        dich_cum(M, L, a, u * (d - khe))
    return "khe %s–mang = %.3f (do lai %.4f)" % (a, khe, kc(M, a, "mang"))


def tam(M, n):
    return trimesh.Trimesh(M[n]["v"], M[n]["f"], process=False).center_mass


def phong_cum(M, L, ten, k, c=None):
    c = tam(M, ten) if c is None else c
    bien_doi_cum(M, L, ten, lambda V: c + (V - c) * k)


def phong_phang(M, L, ten, k):
    """Phóng các cạnh trong mặt phẳng của lục lạp/granum (giữ bề dày theo trục ngắn)."""
    V = M[ten]["v"]; c = V.mean(0); _, U = np.linalg.eigh(np.cov((V - c).T)); n = U[:, 0]
    def f(X):
        Y = X - c; h = (Y @ n)[:, None] * n; return c + h + (Y - h) * k
    bien_doi_cum(M, L, ten, f)


def xoay_cum(M, L, ten, goc, truc=2):
    """Xoay cụm quanh trục chính số 'truc' (0 = ngắn nhất/pháp tuyến, 2 = dài nhất) qua tâm."""
    V = M[ten]["v"]; c = V.mean(0); _, U = np.linalg.eigh(np.cov((V - c).T)); ax = U[:, truc]
    Rm = trimesh.transformations.rotation_matrix(np.radians(goc), ax)[:3, :3]
    bien_doi_cum(M, L, ten, lambda X: c + (X - c) @ Rm.T)


def ke_hoach(tbn, M, L, base):
    P = []
    def A(loai, fn, dich, bp, kieu="dot_bien", ngoai_le=(), bp2=None):
        P.append({"loai": loai, "fn": fn, "dich": dich, "bp": bp, "kieu": kieu, "ngoai_le": list(ngoai_le), "bp2": bp2})
    LL = tl(L, "luc_lap"); TT = tl(L, "ti_the"); GO = tl(L, "golgi"); PE = tl(L, "peroxisome")
    if tbn == "TV6":
        A("lop6_bat_3_ti_the", lambda M, L: [them(M, L, "ti_the_x%d" % i, D.nang(np.array([x, 12.85, -1.0]) - [0.4, 0, 0], np.array([x, 12.85, -1.0]) + [0.4, 0, 0], 0.35), "ti_the") for i, x in enumerate((-18, 0, 18))] and "them 3 ti the", ["2.8"], "ti_the")
        A("lop6_ve_golgi", lambda M, L: them(M, L, "golgi_x", D.tru(0.3, 0.4, np.array([0, 13.6 - 0.55, -5.0]), [0, 1, 0]), "golgi") or "them golgi o lop 6", ["2.9"], "golgi")
        A("dc_dich_ca_te_bao", lambda M, L: [M[n].__setitem__("v", M[n]["v"] + [3.0, -2.0, 1.0]) for n in M] and "dich ca te bao", [], "*", "trung_tinh")
        return P
    nh = "nhan"; mn = "mang_nhan"; er = "luoi_noi_chat"
    l0, l1, l2, l3, l4 = LL[0], LL[1], LL[2], LL[12], LL[5]
    A("bo_khong_bao", lambda M, L: xoa(M, L, "khong_bao") or "xoa khong bao", ["2.1"], "khong_bao")
    A("keo_dai_ca_te_bao_x1_03", lambda M, L: [M[n].__setitem__("v", M[n]["v"] * [1.03, 1, 1]) for n in M] and "keo dai x 1,03", ["2.2a"], "vach")
    A("keo_dai_ca_te_bao_x1_015_doi_chung", lambda M, L: [M[n].__setitem__("v", M[n]["v"] * [1.015, 1, 1]) for n in M] and "keo dai x 1,015", ["2.2a"], "vach", "sat_nguong_dat")
    P[-1]["toan_te_bao"] = P[-2]["toan_te_bao"] = True
    A("thanh_goc_nhon", lambda M, L: dat_m(M, "vach", mf.Manifold.cube((60, 30, 30), True) - D.hop(1.2)) or "thanh hop goc nhon", ["2.2b"], "vach")
    A("mang_lech_ra_ngoai", lambda M, L: M["mang"].__setitem__("v", M["mang"]["v"] + [1.5, 0, 0]) or "dich mang 1,5", ["2.3a"], "mang")
    A("thanh_day_0_9", lambda M, L: dat_m(M, "vach", D.hop(0) - D.hop(0.9)) or "thanh day 0,9 (mang van o 1,2: co khe)", ["2.3b"], "vach")
    A("thanh_day_1_05_doi_chung", lambda M, L: dat_m(M, "vach", D.hop(0.15) - D.hop(1.2)) or "thanh 1,05 (ngoai lui 0,15)", ["2.3b"], "vach", "sat_nguong_dat")
    A("mang_ho_thanh_0_05", lambda M, L: dat_m(M, "mang", D.hop(1.25) - D.hop(1.4)) or "mang cach thanh 0,05", ["2.3c"], "mang")
    A("mang_day_0_14", lambda M, L: dat_m(M, "mang", D.hop(1.2) - D.hop(1.34)) or "mang day 0,14", ["2.3d"], "mang")
    A("mang_day_0_16_doi_chung", lambda M, L: dat_m(M, "mang", D.hop(1.2) - D.hop(1.36)) or "mang day 0,16", ["2.3d"], "mang", "sat_nguong_dat")
    def tach_kb(M, L):
        k = man(M, "khong_bao"); a, b = k.split_by_plane((1, 0, 0), 0.25); b = b ^ mf.Manifold.cube((200, 200, 200), True).translate((-100 - 0.25, 0, 0))
        dat_m(M, "khong_bao", a); them(M, L, "khong_bao_2", k.split_by_plane((1, 0, 0), -0.25)[1], "khong_bao"); return "chia khong bao lam 2 (khe 0,5)"
    A("hai_khong_bao", tach_kb, ["2.4a"], "khong_bao")
    A("khong_bao_thu_0_9", lambda M, L: phong_cum(M, L, "khong_bao", 0.9) or "khong bao x0,9", ["2.4b"], "khong_bao")
    A("khong_bao_khoet_lo", lambda M, L: dat_m(M, "khong_bao", man(M, "khong_bao") - D.sphere([0, 0, 2.0], 4.0)) or "khoet cau r4 trong khong bao", ["2.4c"], "khong_bao")
    def u_kb(depth):
        def f(M, L):
            c = np.array([9.0, 0.0, 15 - 1.4 - depth - 0.8])
            dat_m(M, "khong_bao", man(M, "khong_bao") + D.sphere(c, 0.8, 64)); return "u khong bao toi do sau %.2f" % depth
        return f
    A("khong_bao_sat_mang_0_9", u_kb(0.9), ["2.5b"], "khong_bao")
    A("khong_bao_sat_mang_0_98", u_kb(0.98), ["2.5b"], "khong_bao")
    A("khong_bao_sat_mang_1_03_doi_chung", u_kb(1.03), ["2.5b"], "khong_bao", "sat_nguong_dat")
    A("hai_nhan", lambda M, L: them(M, L, "nhan_2", D.sphere([-9.0, 0, -11.5], 0.5), "nhan") or "them nhan thu hai", ["2.6a"], "nhan")
    A("nhan_xa_mang_0_55", lambda M, L: khe_mang(M, L, nh, 0.55), ["2.6b"], nh)
    A("nhan_xa_mang_0_45_doi_chung", lambda M, L: khe_mang(M, L, nh, 0.45), ["2.6b"], nh, "sat_nguong_dat")
    A("nhan_nho_ti_le", lambda M, L: phong_cum(M, L, nh, 0.62) or "nhan x0,62 (D 3,7)", ["2.6e"], nh)
    A("nhan_d_7_2", lambda M, L: phong_cum(M, L, nh, 1.2) or "nhan D 7,2", ["2.6f"], nh)
    A("nhan_d_6_95_doi_chung", lambda M, L: phong_cum(M, L, nh, 6.95 / 6.0) or "nhan D 6,95", ["2.6f"], nh, "sat_nguong_dat")
    A("nhan_det", lambda M, L: bien_doi_cum(M, L, nh, lambda V: (lambda c: c + (V - c) * [1, 1, 0.85])(M[nh]["v"].mean(0))) or "nhan nen z 0,85", ["2.6f"], nh)
    A("luc_lap_xa_mang_0_6", lambda M, L: khe_mang(M, L, l0, 0.6), ["2.7b"], l0)
    A("luc_lap_xa_mang_0_52", lambda M, L: khe_mang(M, L, l1, 0.52), ["2.7b"], l1)
    A("luc_lap_xa_mang_0_45_doi_chung", lambda M, L: khe_mang(M, L, l2, 0.45), ["2.7b"], l2, "sat_nguong_dat")
    def elip(M, L):
        V = M[l3]["v"]; c = V.mean(0); _, U = np.linalg.eigh(np.cov((V - c).T)); n = U[:, 0]
        u, v, nn = D.khung(n); e = mf.Manifold.sphere(1.0, 64).scale((3.0, 3.0, 1.1))
        dat_m(M, l3, D.dat(e, c, u, v, nn)); return "luc lap thanh elip 6 x 2,2"
    A("luc_lap_elip", elip, ["2.7c"], l3)
    A("luc_lap_d_7_3", lambda M, L: phong_phang(M, L, LL[13], 7.3 / 6.0) or "D 7,3", ["2.7d"], LL[13], "st_dat_ts_truot")
    A("luc_lap_d_7_05", lambda M, L: phong_phang(M, L, LL[14], 7.05 / 6.0) or "D 7,05", ["2.7d"], LL[14], "st_dat_ts_truot")
    A("luc_lap_d_6_95_doi_chung", lambda M, L: phong_phang(M, L, LL[15], 6.95 / 6.0) or "D 6,95", ["2.7d"], LL[15], "sat_nguong_dat")
    A("luc_lap_con_29", lambda M, L: xoa(M, L, LL[-1]) or "con 29 luc lap", ["2.7e"], "luc_lap")
    A("luc_lap_xoay_40", lambda M, L: xoay_cum(M, L, LL[16], 40) or "xoay 40°", ["2.7f"], LL[16])
    A("luc_lap_xoay_20_doi_chung", lambda M, L: (dich_cum(M, L, LL[17], -huong_mang(M, LL[17]) * 1.0), xoay_cum(M, L, LL[17], 20), "xoay 20° (lui vao 1 µm)")[2], ["2.7f"], LL[17], "sat_nguong_dat")
    A("luc_lap_chi_2_mat", lambda M, L: [xoa(M, L, n) for n in LL if huong_mang(M, n)[2] == 0] and "chi giu luc lap mat tren/duoi", ["2.7b"], "luc_lap")
    A("ti_the_khe_luc_lap_0_29", lambda M, L: khe_toi(M, L, TT[0], tl(L, "luc_lap")[18], 0.29), ["K.1"], TT[0], bp2=tl(L, "luc_lap")[18])
    A("ti_the_khe_luc_lap_0_31_doi_chung", lambda M, L: khe_toi(M, L, TT[1], tl(L, "luc_lap")[19], 0.31), ["K.1"], TT[1], "sat_nguong_dat", bp2=tl(L, "luc_lap")[19])
    gr = [c for c, p in L["cha"].items() if p == LL[3] and L["loai"][c] == "granum"]
    A("granum_khe_0_09", lambda M, L: khe_toi(M, L, gr[0], gr[1], 0.09), ["K.1"], gr[0], bp2=gr[1])
    A("granum_khe_0_11_doi_chung", lambda M, L: khe_toi(M, L, gr[3], gr[4], 0.11), ["K.1"], gr[3], "sat_nguong_dat", bp2=gr[4])
    A("nhan_con_sat_mang_nhan", lambda M, L: dich_cum(M, L, "nhan_con_1", [0, 0, 1.85]) or "dich nhan con 1,85 (khe ~0,23)", ["K.1"], "nhan_con_1")
    A("peroxisome_ra_ngoai_mang", lambda M, L: dich_cum(M, L, PE[0], huong_mang(M, PE[0]) * 1.2) or "day ra 1,2", ["M.a"], PE[0])
    A("granum_ra_ngoai_luc_lap", lambda M, L: dich_cum(M, L, gr[5], (lambda V: (np.linalg.eigh(np.cov((V - V.mean(0)).T))[1][:, 0]))(M[LL[3]]["v"]) * 1.5) or "day granum ra", ["M.a"], gr[5])
    A("golgi_xa_mang_1_2", lambda M, L: khe_mang(M, L, GO[0], 1.2), ["M.d"], GO[0])
    A("golgi_xa_mang_1_05", lambda M, L: khe_mang(M, L, GO[1], 1.05), ["M.d"], GO[1])
    A("golgi_xa_mang_0_95_doi_chung", lambda M, L: khe_mang(M, L, GO[2], 0.95), ["M.d"], GO[2], "sat_nguong_dat")
    def tach_pe(M, L):
        m = man(M, PE[1]); a, b = m.split_by_plane((1, 0, 0), float(M[PE[1]]["v"].mean(0)[0]))
        dat_m(M, PE[1], a + b.translate((-0.05, 0, 0))); return "cat peroxisome lam 2 manh"
    A("peroxisome_hai_manh", tach_pe, ["M.e"], PE[1])
    def gop_pe(M, L):
        a, b = PE[2], PE[3]; V = np.vstack([M[a]["v"], M[b]["v"]]); F = np.vstack([M[a]["f"], M[b]["f"] + len(M[a]["v"])])
        M[a] = {"v": V, "f": F}; xoa(M, L, b); return "gop %s vao %s (mot vat)" % (b, a)
    A("gop_hai_peroxisome", gop_pe, ["M.g"], PE[2])
    th = [c for c, p in L["cha"].items() if L["loai"][c] == "thylakoid" and p == gr[2]]
    th = sorted(th)
    A("thylakoid_chong", lambda M, L: dich_cum(M, L, th[0], (M[th[1]]["v"].mean(0) - M[th[0]]["v"].mean(0)) * 0.5) or "thylakoid chong nhau", ["M.h"], th[0])
    def pe_chong(M, L):
        a, b = PE[4], PE[5]; ca = M[a]["v"].mean(0); cb = M[b]["v"].mean(0); u = (ca - cb) / np.linalg.norm(ca - cb)
        dich_cum(M, L, a, cb + u * 0.7 * 0.7 - ca); return "tam %s cach tam %s 0,49 µm (chong)" % (a, b)
    A("peroxisome_chong", pe_chong, ["M.h"], PE[4], bp2=PE[5])
    if "nhan_con_1" in M:
        A("nhan_con_3_cai", lambda M, L: [them(M, L, "nhan_con_x%d" % i, D.sphere(M["nhan_con_1"]["v"].mean(0) + o, 0.5), "nhan_con", "nhan") for i, o in enumerate(([1.4, 0, 0], [-1.4, 0, 0]))] and "them 2 nhan con", ["5.1a"], "nhan_con", "st_dat_ts_truot")
        A("nhan_con_d_2_2", lambda M, L: phong_cum(M, L, "nhan_con_1", 2.2 / 1.5) or "nhan con D 2,2", ["5.1b"], "nhan_con_1")
        A("nhan_con_d_1_95_doi_chung", lambda M, L: phong_cum(M, L, "nhan_con_1", 1.95 / 1.5) or "D 1,95", ["5.1b"], "nhan_con_1", "sat_nguong_dat")
    def env_kh(KH, L1=0.05):
        def f(M, L):
            c = tam(M, nh); Rn = float(np.linalg.norm(M[nh]["v"] - c, axis=1).max())
            r1, r2, r3, r4 = Rn, Rn - L1, Rn - L1 - KH, Rn - 2 * L1 - KH
            lop = D.sphere(c, r1, 96) - D.sphere(c, r4, 96)
            dirs = [(M[p]["v"].mean(0) - c) / np.linalg.norm(M[p]["v"].mean(0) - c) for p in tl(L, "lo_nhan")]
            k = (D.sphere(c, r2, 96) - D.sphere(c, r3, 96)) - mf.Manifold.batch_boolean([D.tru(0.6, 0.175, c + d * (Rn - 0.1), d, 24) for d in dirs], mf.OpType.Add) if KH > 0 else None
            e = lop if k is None else lop - k
            e = e - mf.Manifold.batch_boolean([D.tru(0.6, 0.125, c + d * (Rn - 0.1), d, 24) for d in dirs], mf.OpType.Add)
            dat_m(M, mn, e)
            for p, d in zip(tl(L, "lo_nhan"), dirs):
                dat_m(M, p, D.tru(0.6, 0.125, c + d * (Rn - 0.1), d, 24) ^ lop)
            return "mang nhan: lop %.2f, khoang %.2f" % (L1, KH)
        return f
    A("mang_nhan_mot_lop", env_kh(0.0, 0.085), ["5.2a"], mn)
    A("mang_nhan_khoang_0_12", env_kh(0.12), ["5.2c"], mn)
    A("lo_nhan_con_29", lambda M, L: [xoa(M, L, p) for p in tl(L, "lo_nhan")[:3]] and "xoa 3 lo nhan", ["5.2b"], "lo_nhan")
    lo0 = tl(L, "lo_nhan")[5]
    A("lo_nhan_rong_0_35", lambda M, L: phong_cum(M, L, lo0, 0.35 / 0.25) or "lo rong 0,35", ["5.2b"], lo0)
    A("ti_the_dai_2_2", lambda M, L: CA.cl_dai(M, L, TT[2], 2.2), ["5.3b"], TT[2])
    A("ti_the_rong_0_45", lambda M, L: CA.cl_rong(M, L, TT[3], 0.45), ["5.3b"], TT[3])
    A("ti_the_rong_0_52_doi_chung", lambda M, L: CA.cl_rong(M, L, TT[4], 0.52), ["5.3b"], TT[4], "sat_nguong_dat")
    A("ti_the_con_9", lambda M, L: xoa(M, L, TT[-1]) or "con 9 ti the", ["5.3c"], "ti_the")
    mao = sorted(c for c, p in L["cha"].items() if p == TT[5])
    A("ti_the_3_mao", lambda M, L: xoa(M, L, mao[0]) or "xoa 1 mao (con 3)", ["5.3e"], TT[5], "st_dat_ts_truot")
    mao6 = sorted(c for c, p in L["cha"].items() if p == TT[6])
    A("mao_day_0_06", lambda M, L: bien_doi_cum(M, L, mao6[1], lambda V: (lambda c, ax: c + (V - c) + (((V - c) @ ax)[:, None] * ax) * 0.5)(V.mean(0), np.array([1.0, 0, 0]))) or "mao day 0,06", ["5.3e"], mao6[1], "st_dat_ts_truot")
    A("er_tach_nhan_0_02", lambda M, L: CA.cl_er_tach(M, L, 0.02), ["5.4a"], er)
    def er_vanh(day, d_ong=0.15):
        def f(M, L):
            c = tam(M, nh); Rn = float(np.linalg.norm(M[nh]["v"] - c, axis=1).max())
            vanh = (D.sphere(c, Rn + day, 96) - D.sphere(c, Rn, 96)) ^ D.hop(1.4 + 0.31)
            zt = -15 + 1.4 + 0.45
            ong = [D.nang(np.array([0, 0, zt]) + np.array(d) * 1.5, np.array([0, 0, zt]) + np.array(d) * 8.0, d_ong / 2, 16) for d in ([1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0])]
            dat_m(M, er, vanh + mf.Manifold.batch_boolean(ong, mf.OpType.Add)); return "vanh ER %.2f, ong %.2f" % (day, d_ong)
        return f
    A("er_vanh_1_3", er_vanh(1.3), ["5.4b"], er)
    A("er_ong_0_3", er_vanh(0.7, 0.3), ["5.4b"], er)
    A("golgi_gan_nhau_1_5", lambda M, L: khe_toi(M, L, GO[3], GO[4], 1.5), ["5.6a"], GO[3], "st_dat_ts_truot", bp2=GO[4])
    A("golgi_gan_nhau_1_95", lambda M, L: khe_toi(M, L, GO[4], GO[3], 1.95), ["5.6a"], GO[4], "st_dat_ts_truot", bp2=GO[3])
    A("golgi_con_4", lambda M, L: xoa(M, L, GO[-1]) or "con 4 golgi", ["5.6b"], "golgi")
    tui = sorted(c for c, p in L["cha"].items() if p == GO[2])
    A("golgi_3_tui", lambda M, L: [xoa(M, L, t) for t in tui[:2]] and "xoa 2 tui", ["5.6c"], GO[2])
    tui1 = sorted(c for c, p in L["cha"].items() if p == GO[1])
    A("tui_khe_0_05", lambda M, L: dich_cum(M, L, tui1[-1], huong_mang(M, GO[1]) * -0.02) or "tui cuoi lech 0,02", ["5.6c"], tui1[-1])
    A("peroxisome_d_1_1", lambda M, L: phong_cum(M, L, PE[6], 1.1 / 0.7) or "D 1,1", ["5.7a"], PE[6])
    A("peroxisome_d_0_98_doi_chung", lambda M, L: phong_cum(M, L, PE[7], 0.98 / 0.7) or "D 0,98", ["5.7a"], PE[7], "sat_nguong_dat")
    A("peroxisome_con_19", lambda M, L: xoa(M, L, PE[-1]) or "con 19", ["5.7c"], "peroxisome")
    def pe_xa(khe):
        def f(M, L):
            p = PE[8]; ll = min(tl(L, "luc_lap"), key=lambda l: kc(M, p, l))
            for _ in range(12):
                d = kc(M, p, ll)
                if abs(d - khe) < 2e-4: break
                u = M[p]["v"].mean(0) - M[ll]["v"].mean(0); u -= (u @ huong_mang(M, p)) * huong_mang(M, p); u /= np.linalg.norm(u)
                dich_cum(M, L, p, u * (khe - d))
            return "peroxisome xa luc lap %.2f (do %.4f)" % (khe, kc(M, p, ll))
        return f
    A("peroxisome_xa_luc_lap_1_2", pe_xa(1.2), ["5.7d"], PE[8])
    for rid, lo_ in (("5.8", "lysosome"), ("5.9", "trung_the"), ("5.10", "cau_sinh_chat"), ("5.11", "khung_xuong")):
        A("them_%s" % lo_, lambda M, L, lo_=lo_: them(M, L, lo_ + "_1", D.sphere([-9.0, 12.5, -9.0], 0.3), lo_) or "them " + lo_, [rid], lo_)
    lp = LL[4]; mt = [c for c, p in L["cha"].items() if p == lp and L["loai"][c] == "mang_trong_lap"][0]
    def mt_lui(lui):
        def f(M, L):
            V = M[lp]["v"]; c = V.mean(0); _, U = np.linalg.eigh(np.cov((V - c).T)); n = U[:, 0]
            dat_m(M, mt, D.thau_kinh(c, n, 6.0, 2.2, lui=lui)); return "mang trong lui %.3f" % lui
        return f
    A("mang_trong_khe_0_08", mt_lui(0.08), ["10.1a"], mt, "st_dat_ts_truot")
    A("mang_trong_khe_0_055_doi_chung", mt_lui(0.055), ["10.1a"], mt, "sat_nguong_dat")
    lp6 = LL[6]; mt6 = [c for c, p in L["cha"].items() if p == lp6 and L["loai"][c] == "mang_trong_lap"][0]
    def gap_nep(M, L):
        V = M[lp6]["v"]; c = V.mean(0); _, U = np.linalg.eigh(np.cov((V - c).T)); n = U[:, 0]; u = U[:, 2]
        ran = D.tru(0.3, 0.4, c + u * 2.8 + n * 0.0, U[:, 1], 24).scale((1, 1, 1))
        dat_m(M, mt6, man(M, mt6) - D.dat(mf.Manifold.cube((0.8, 0.2, 0.9), True), c + u * 2.6, U[:, 1], n, u)); return "khoet ranh (gap nep) vao mang trong"
    A("mang_trong_gap_nep", gap_nep, ["10.1b"], mt6)
    grl = [c for c, p in L["cha"].items() if p == LL[7] and L["loai"][c] == "granum"]
    A("granum_con_9", lambda M, L: xoa(M, L, sorted(grl)[0]) or "xoa 1 granum", ["10.2a"], LL[7], "st_dat_ts_truot")
    gr8 = sorted(c for c, p in L["cha"].items() if p == LL[8] and L["loai"][c] == "granum")
    A("granum_d_0_65", lambda M, L: phong_phang(M, L, gr8[0], 0.65 / 0.45) or "granum D 0,65", ["10.2b"], gr8[0])
    th8 = sorted(c for c, p in L["cha"].items() if p == gr8[1])
    A("granum_2_thylakoid", lambda M, L: [xoa(M, L, t) for t in th8[:2]] and "con 2 thylakoid", ["10.2c"], gr8[1])
    th8b = sorted(c for c, p in L["cha"].items() if p == gr8[2])
    A("thylakoid_day_0_07", lambda M, L: bien_doi_cum(M, L, th8b[0], lambda V: (lambda c, n: c + (V - c) + (((V - c) @ n)[:, None] * n) * 0.4)(V.mean(0), np.linalg.eigh(np.cov((V - V.mean(0)).T))[1][:, 0])) or "thylakoid day 0,07", ["10.2d"], th8b[0])
    tn9 = sorted(c for c, p in L["cha"].items() if p == LL[9] and L["loai"][c] == "tam_noi")
    A("bo_tam_noi_giua_hang", lambda M, L: xoa(M, L, tn9[1]) or "xoa 1 tam noi trong hang", ["10.3a"], LL[9])
    tn10 = sorted(c for c, p in L["cha"].items() if p == LL[10] and L["loai"][c] == "tam_noi")
    gr10 = sorted((c for c, p in L["cha"].items() if p == LL[10] and L["loai"][c] == "granum"), key=lambda s: int(s.rsplit("_", 1)[1]))
    A("granum_dau_khong_tam_noi", lambda M, L: xoa(M, L, tn10[0]) or "xoa tam noi dau hang", ["10.3b"], gr10[0])
    def tam_hep(M, L):
        n = tn10[2]; V = M[n]["v"]; c = V.mean(0); _, U = np.linalg.eigh(np.cov((V - c).T)); w = U[:, 1]
        M[n]["v"] = c + (V - c) - (((V - c) @ w)[:, None] * w) * (1 - 0.05 / 0.15); return "tam noi hep 0,05 (cham 0,0025 µm²)"
    tn12 = sorted(c for c, p in L["cha"].items() if p == LL[12] and L["loai"][c] == "tam_noi")
    gr12 = sorted((c for c, p in L["cha"].items() if p == LL[12] and L["loai"][c] == "granum"), key=lambda s: int(s.rsplit("_", 1)[1]))
    def tam_hep2(M, L):
        n = [t for t in tn12 if t.endswith("_1")][0]; V = M[n]["v"]; c = V.mean(0); _, U = np.linalg.eigh(np.cov((V - c).T)); w = U[:, 1]
        M[n]["v"] = c + (V - c) - (((V - c) @ w)[:, None] * w) * (1 - 0.05 / 0.15); return "tam noi dau hang hep 0,05 (granum dau hang chi con vung cham 0,0025 µm²)"
    A("tam_noi_hep", tam_hep2, ["10.3b"], gr12[0], "st_dat_ts_truot")
    A("tam_noi_nghieng_20", lambda M, L: xoay_cum(M, L, tn10[5], 20, truc=0) or "tam noi xoay 20° trong mat phang (quanh phap tuyen)", ["10.3c"], tn10[5])
    tn11 = sorted(c for c, p in L["cha"].items() if p == LL[11] and L["loai"][c] == "tam_noi")
    A("tam_noi_it_hon_n_tru_1", lambda M, L: xoa(M, L, tn11[0]) or "con 8 tam noi / 10 granum", ["10.5"], LL[11])
    def golgi_gom(M, L):
        c = tam(M, nh)
        for i, g in enumerate(GO[:3]):
            u = huong_mang(M, g); p = c + np.array([(-2.5, 0, 0), (2.5, 0, 0), (0, 2.5, 0)][i]) + np.array([0, 0, 2.6])
            dich_cum(M, L, g, p - M[g]["v"].mean(0))
        return "dua 3/5 the Golgi ve sat nhan"
    A("golgi_gom_canh_nhan", golgi_gom, ["5.6a"], "golgi")
    def golgi_mot_goc(M, L):
        pos = [(6.0, 12.6, -5.0), (10.0, 12.6, -5.0), (14.0, 12.6, -5.0), (18.0, 12.6, -5.0), (22.0, 12.6, -5.0)]
        for g, p in zip(GO, pos):
            dich_cum(M, L, g, np.array(p) - M[g]["v"].mean(0))
        return "dua ca 5 the Golgi vao mot phan tam (x>0, y>0, z<0)"
    A("golgi_mot_phan_tam", golgi_mot_goc, ["5.6a"], "golgi")
    A("luc_lap_chi_4_mat_doi_chung", lambda M, L: [xoa(M, L, n) for n in tl(L, "luc_lap") if abs(huong_mang(M, n)[0]) > 0] and "xoa luc lap hai mat dau (con 4 mat)", ["2.7b"], "luc_lap", "sat_nguong_dat")
    P[-1]["ten_cu"] = "luc_lap_chi_4_mat"
    A("dan_nhan_sai_ti_the_la_peroxisome", lambda M, L: (L["loai"].__setitem__(TT[7], "peroxisome"), L["bo_phan_du_kien"].__setitem__(TT[7], "peroxisome"), "dan nhan %s la peroxisome" % TT[7])[2], ["nghi_sai_loai"], TT[7])
    A("dan_nhan_sai_luc_lap_la_lysosome", lambda M, L: (L["loai"].__setitem__(LL[25], "lysosome"), L["bo_phan_du_kien"].__setitem__(LL[25], "lysosome"), "dan nhan %s la lysosome" % LL[25])[2], ["nghi_sai_loai"], LL[25])
    A("thieu_tang2_mao", lambda M, L: [xoa(M, L, m_) for m_ in tl(L, "mao")] and "khong dung mao nao (tang 2)", ["5.3e"], TT[8])
    A("dc_dich_ca_te_bao", lambda M, L: [M[n].__setitem__("v", M[n]["v"] + [3.0, -2.0, 1.0]) for n in M] and "dich ca te bao (3,-2,1)", [], "*", "trung_tinh")
    A("dc_chia_tam_giac_luc_lap", lambda M, L: M[LL[20]].update(dict(zip(["v", "f"], (lambda m: (np.asarray(m.vertices), np.asarray(m.faces)))(trimesh.Trimesh(M[LL[20]]["v"], M[LL[20]]["f"], process=False).subdivide())))) or "chia tam giac luc lap", [], LL[20], "trung_tinh")
    return P


def ke_hoach_v12(tbn, M, L, base):
    P = []
    def A(loai, fn, dich, bp, kieu, bp2=None):
        P.append({"loai": "v12_" + loai, "fn": fn, "dich": dich, "bp": bp, "kieu": kieu, "ngoai_le": [], "bp2": bp2})
    if tbn in ("TV6", "TV2"):
        return P
    LL = tl(L, "luc_lap"); GO = tl(L, "golgi"); PE = tl(L, "peroxisome"); er = "luoi_noi_chat"
    d1 = LV.KHOANG("tv", "2.7d", "đường kính")[1]; dv1 = LV.P("tv.kich_thuoc_luc_lap_ve")["duong_kinh"][1]
    dk = min(d1 - 0.5, dv1 + 1.0)
    def ll_dk(M, L, a=LL[3]):
        e = kiem2.BoPhan("x", "x", M[a]["v"], M[a]["f"]).kich_thuoc()[0]; phong_phang(M, L, a, dk / e); return "phong mat phang %s: duong kinh %.2f -> %.2f" % (a, e, dk)
    A("luc_lap_duong_kinh_%.1f" % dk, ll_dk, ["2.7d"], LL[3], "st_dat_ts_truot")
    k0 = LV.P("tv.khe_nhan_mang")[0]
    A("nhan_khe_mang_%.2f" % (k0 - 0.1), lambda M, L: khe_mang(M, L, "nhan", k0 - 0.1), ["2.6b"], "nhan", "st_dat_ts_truot")
    k0 = LV.P("tv.khe_luc_lap_mang")[0]
    A("luc_lap_khe_mang_%.2f" % (k0 - 0.1), lambda M, L, a=LL[6]: khe_mang(M, L, a, k0 - 0.1), ["2.7b"], LL[6], "st_dat_ts_truot")
    if len(GO) >= 2:
        kg = LV.P("tv.khe_golgi_golgi")
        A("golgi_khe_%.1f" % (kg / 2), lambda M, L: khe_toi(M, L, GO[0], GO[1], kg / 2), ["5.6a"], GO[0], "st_dat_ts_truot")
    if PE:
        p0 = LV.P("tv.duong_kinh_peroxisome_ve")[0]; f0 = LV.KHOANG("tv", "5.7a", "đường")[0]; canh = (p0 + f0) / 2 + 0.05
        def pe_nho(M, L, a=PE[0]):
            e = kiem2.BoPhan("x", "x", M[a]["v"], M[a]["f"]).kich_thuoc()[0]; phong_cum(M, L, a, canh / e); return "phong %s: canh %.3f -> %.3f" % (a, e, canh)
        A("peroxisome_canh_%.2f" % canh, pe_nho, ["5.7a"], PE[0], "st_dat_ts_truot")
    if er in M and "nhan" in M:
        def er_vao(M, L):
            u = tam(M, "nhan") - M[er]["v"].mean(0); u /= max(np.linalg.norm(u), 1e-12)
            dich_cum(M, L, er, u * 0.15); return "dich ER 0,15 µm vao phia tam nhan (cap ER–nhan/mang nhan truoc day mien khe)"
        A("er_lan_vao_nhan_M.c0", er_vao, ["M.c₀@su_that"], er, "dot_bien")
    if PE and LL:
        a1, b1, t1, g1 = _pe_ll_chon(M, L, PE, LL, 0.1, bat_dau=1, cho_con=True)
        A("peroxisome_chong_luc_lap_M.c0", lambda M, L, a=a1, b=b1, t=t1, g=g1: (dich_cum(M, L, a, t), "dich %s (huong tranh khong bao) toi khe 0,05 voi %s roi theo phap tuyen tiep xuc lan 0,1 µm; giao truoc kiem %s" % (a, b, g))[1],
          ["M.c₀@su_that"], a1, "dot_bien", bp2=b1)
        _n2 = []; a2, b2, t2, g2 = _pe_ll_chon(M, L, PE, LL, 0.02, bat_dau=2, cho_con=False, ra_phap_tuyen=_n2)
        for dl, kieu in ((0.02, "dot_bien"), (0.005, "sat_nguong_dat")):
            t_ = t2 + (dl - 0.02) * _n2[0]
            A("peroxisome_chong_luc_lap_%g_M.c0" % dl, lambda M, L, a=a2, b=b2, t=t_, dl=dl, g=g2: (dich_cum(M, L, a, t), "dich %s (huong tranh khong bao) toi khe 0,05 voi %s roi theo phap tuyen tiep xuc lan %g µm (cap chon o 0,02; giao truoc kiem o 0,02: %s)" % (a, b, dl, g))[1],
              ["M.c₀@su_that"], a2, kieu, bp2=(b2 if kieu == "dot_bien" else None))
    return P


def _pe_ll_chon(M, L, PE, LL, dl, bat_dau=0, cho_con=False, ra_phap_tuyen=None):
    """Chọn (peroxisome a, lục lạp b, phép dịch t): a trượt theo hướng cặp điểm gần nhất a→b (bỏ thành phần hướng vào không bào)
    tới khi khe còn 0,05 µm, rồi theo hướng cặp điểm gần nhất đi hết khe và lấn dl; sau khi dịch, khối a CHỈ giao b (cho_con: và con cháu của b); trả thêm thể tích giao đã kiểm."""
    def hop(n):
        return M[n]["v"].min(0), M[n]["v"].max(0)
    ten_khac = [n for n in M if L["loai"].get(n) not in (None,)]
    for k in range(len(PE)):
        a = PE[(bat_dau + k) % len(PE)]
        ho_a = set(con_chau(L, a)); c = M[a]["v"].mean(0)
        for b in sorted(LL, key=lambda n: np.linalg.norm(M[n]["v"].mean(0) - c))[:3]:
            p, q, _ = CA.kc_gan(M, a, b); d = kc(M, a, b)
            if d > 1.0:
                continue
            u0 = (q - p) / max(np.linalg.norm(q - p), 1e-12); hs = [u0]
            if "mang" in M:
                pm, qm, _ = CA.kc_gan(M, a, "mang"); nm = (qm - pm) / max(np.linalg.norm(qm - pm), 1e-12)
                w = u0 - (u0 @ nm) * nm
                if np.linalg.norm(w) > 0.2:
                    w /= np.linalg.norm(w)
                    hs += [np.cos(np.radians(g_)) * w + np.sin(np.radians(g_)) * nm for g_ in np.arange(0, 32.5, 2.5)]
            m0, mb = man(M, a), man(M, b)
            for u in hs:
                gap = lambda s_: float(m0.translate(tuple((u * s_).tolist())).min_gap(mb, 5.0))
                s_lo, s_hi = 0.0, None
                for s_ in np.arange(d, d + 2.5, 0.05):
                    if gap(s_) <= 0: s_hi = s_; break
                    s_lo = s_
                if s_hi is None:
                    continue
                for _ in range(40):
                    sm = 0.5 * (s_lo + s_hi)
                    if gap(sm) > 0.05: s_lo = sm
                    else: s_hi = sm
                t0_ = u * s_lo; M_ = {a: {"v": M[a]["v"] + t0_, "f": M[a]["f"]}, b: M[b]}
                p2, q2, _ = CA.kc_gan(M_, a, b); n2 = (q2 - p2) / max(np.linalg.norm(q2 - p2), 1e-12); g2 = gap(s_lo)
                t = t0_ + n2 * (g2 + dl)
                ma = man(M, a).translate(tuple(t.tolist())); lo = M[a]["v"].min(0) + t - 1e-6; hi = M[a]["v"].max(0) + t + 1e-6
                cho = {b} | (set(con_chau(L, b)) if cho_con else set())
                giao = {}; hong = False
                for n in ten_khac:
                    if n in ho_a: continue
                    l2, h2 = hop(n)
                    if np.any(h2 < lo) or np.any(l2 > hi): continue
                    v = float((ma ^ man(M, n)).volume())
                    if v > 1e-12:
                        giao[n] = round(v, 8)
                        if n not in cho: hong = True; break
                if os.environ.get("GO_PE"): print("  thu", a, b, u.round(3), "giao", giao, flush=True)
                if not hong and b in giao:
                    if ra_phap_tuyen is not None: ra_phap_tuyen.append(n2)
                    return a, b, t, giao
    raise RuntimeError("khong tim duoc peroxisome chi giao luc lap (dl = %g)" % dl)


def nhac(c, bp, bp2=None):
    b = c.get("bo_phan")
    if bp2 is not None:
        return isinstance(b, list) and bp in b and bp2 in b
    return bp == "*" or b == bp or (isinstance(b, list) and bp in b)


def k_nv(c): return json.dumps(c.get("bo_phan"), ensure_ascii=False) + "|" + str(c.get("nguong"))


def moi_truot(base, r):
    B = {x["id"]: x for x in base["luat"]}; out = {}
    for x in r["luat"]:
        cu = set(k_nv(c) for c in B[x["id"]]["ca_truot"])
        m = [c for c in x["ca_truot"] if k_nv(c) not in cu]
        if m: out[x["id"]] = m
    return out


def chay(ds=("TV1", "TV2", "TV6"), chi=None):
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    kq_tep = os.path.join(HERE, "cay_loi_tv_ket_qua.json")
    KQ = {"ket_qua": []}
    if chi and os.path.exists(kq_tep):
        KQ = json.load(open(kq_tep)); KQ["ket_qua"] = [x for x in KQ["ket_qua"] if x["loai"] not in chi]
    for tbn in ds:
        M0, L0 = doc(tbn)
        dg = os.path.join(OUT, tbn + "__goc"); ghi(M0, L0, dg)
        base = _K3TV(dg)
        P = ke_hoach(tbn, M0, L0, base) + ke_hoach_v12(tbn, M0, L0, base)
        if tbn == "TV2":
            rb = "ribosome"
            def ri(M, L):
                V = M[rb]["v"].copy(); V[:42] += [0, 0, 0.2]; M[rb]["v"] = V; return "nang 1 cham ribosome 0,2 µm khoi ER"
            P = [{"loai": "ribosome_roi_er", "fn": ri, "dich": ["5.5"], "bp": rb, "kieu": "dot_bien", "ngoai_le": []},
                 {"loai": "dc_dich_ca_te_bao", "fn": lambda M, L: [M[n].__setitem__("v", M[n]["v"] + [1.0, 1.0, -1.0]) for n in M] and "dich ca te bao", "dich": [], "bp": "*", "kieu": "trung_tinh", "ngoai_le": []}]
        co_san = {q["id"] for q in base["luat"]}
        CHIA = os.environ.get("CHIA", "")
        ci, cn = (int(CHIA.split("/")[0]), int(CHIA.split("/")[1])) if CHIA else (0, 1)
        for k, x in enumerate(P):
            if chi and x["loai"] not in chi: continue
            if k % cn != ci: continue
            x = dict(x, dich=V3.map_dich(x["dich"], "tv", co_san), ngoai_le=V3.map_dich(x.get("ngoai_le", []), "tv", co_san))
            M, L = copy.deepcopy(M0), copy.deepcopy(L0)
            rec = {"te_bao": tbn, "loai": x["loai"], "kieu": x["kieu"], "luat_dich": x["dich"], "bo_phan": x["bp"] if not x.get("bp2") else [x["bp"], x["bp2"]]}
            Bmap = {q["id"]: q for q in base["luat"]}
            if x.get("ten_cu"): rec["ten_cu"] = x["ten_cu"]
            truot_bp = lambda rid: rid in Bmap and (any(nhac(c, x["bp"], x.get("bp2")) for c in Bmap[rid]["ca_truot"]) or Bmap[rid]["trang_thai"] == "KHONG_AP_DUNG")
            khl, con, bo = V3.loc_hop_le(x["dich"], x["kieu"], Bmap, truot_bp, x.get("toan_te_bao", False))
            if khl:
                rec["khong_hop_le"] = khl; KQ["ket_qua"].append(rec); print(tbn, x["loai"], "KHONG_HOP_LE", khl, flush=True); continue
            if bo:
                rec["muc_dich_bo_vi_truot_o_goc"] = bo; x = dict(x, dich=con); rec["luat_dich"] = con
            try:
                rec["mo_ta"] = x["fn"](M, L)
            except Exception as ex:
                rec["loi_cay"] = repr(ex); KQ["ket_qua"].append(rec); print(tbn, x["loai"], "LOI CAY", repr(ex), flush=True); continue
            d = os.path.join(OUT, "%s__%02d_%s" % (tbn, k, x["loai"])); ghi(M, L, d)
            r = _K3TV(d)
            mv = moi_truot(base, r); rec["luat_truot_moi"] = sorted(mv)
            if x["kieu"] == "dot_bien":
                if x["dich"] == ["nghi_sai_loai"]:
                    cu = set(y["bo_phan"] for y in base["nghi_sai_loai"]["ds"])
                    hit = [y for y in r["nghi_sai_loai"]["ds"] if y["bo_phan"] == x["bp"] and y["bo_phan"] not in cu]
                else:
                    hit = [c for rid in x["dich"] for c in mv.get(rid, []) if nhac(c, x["bp"], x.get("bp2"))]
                if x["loai"].startswith("thieu_tang2"):
                    lt = [q.get("loai_truot") for q in r["luat"] if q["id"] in x["dich"]]
                    rec["loai_truot"] = lt; hit = hit if all(v == "THIEU_BO_PHAN_TANG2" for v in lt) else []
                rec["phat_hien"] = bool(hit); rec["do_duoc"] = hit[:2]
                rec["phat_hien_theo_phan_quyet"] = V3.theo_pq(x["dich"], mv, lambda c: nhac(c, x["bp"], x.get("bp2")))
            elif x["kieu"] == "st_dat_ts_truot":
                pq = V3.theo_pq(x["dich"], mv, lambda c: nhac(c, x["bp"], x.get("bp2")))
                rec["phat_hien_theo_phan_quyet"] = pq; rec["dung_ky_vong"] = pq["tham_so"] and not pq["su_that"]
            elif x["kieu"] == "sat_nguong_dat":
                g = [c for rid in x["dich"] for c in mv.get(rid, []) if nhac(c, x["bp"], x.get("bp2"))]
                rec["bao_dong_gia"] = bool(g); rec["chi_tiet"] = g[:2]
            else:
                g = sorted(q for q in mv if q not in x["ngoai_le"]); rec["bao_dong_gia"] = bool(g); rec["luat_bao_gia"] = g
                rec["chi_tiet"] = {q: mv[q][:2] for q in g}
            rec["thoi_gian_s"] = r["thoi_gian_s"]
            KQ["ket_qua"].append(rec)
            tag = ("PHAT_HIEN" if rec.get("phat_hien") else "BO_SOT") if x["kieu"] == "dot_bien" else (("DUNG" if rec["dung_ky_vong"] else "SAI") if x["kieu"] == "st_dat_ts_truot" else ("BAO_GIA" if rec["bao_dong_gia"] else "SACH"))
            print("%-4s %-34s %-22s %-9s moi=%s" % (tbn, x["loai"], x["bp"], tag, rec["luat_truot_moi"]), flush=True)
            json.dump(r, open(os.path.join(OUT, "kq__" + os.path.basename(d) + ".json"), "w"), ensure_ascii=False, indent=1, default=kiem2._json_default)
            shutil.rmtree(d, ignore_errors=True)
    db = [x for x in KQ["ket_qua"] if x["kieu"] == "dot_bien" and "phat_hien" in x]
    dc = [x for x in KQ["ket_qua"] if x["kieu"] != "dot_bien" and "bao_dong_gia" in x]
    phu = collections.defaultdict(lambda: [0, 0])
    for x in db:
        for rid in x["luat_dich"]: phu[rid][1] += 1; phu[rid][0] += int(x["phat_hien"])
    ids = sorted({q["id"] for q in base["luat"]}) + ["nghi_sai_loai"]
    KQ["do_phu_theo_luat"] = {i: ("%d/%d" % tuple(phu[i]) if i in phu else "khong co dot bien") for i in ids}
    KQ["tong"] = {"so_dot_bien": len(db), "phat_hien": sum(x["phat_hien"] for x in db), "khong_hop_le": sum(1 for x in KQ["ket_qua"] if "khong_hop_le" in x),
                  "loi_cay": sum(1 for x in KQ["ket_qua"] if "loi_cay" in x), "so_doi_chung": len(dc), "doi_chung_bao_gia": sum(x["bao_dong_gia"] for x in dc),
                  "so_luat_co_dot_bien": sum(1 for v in KQ["do_phu_theo_luat"].values() if v != "khong co dot bien"), "thoi_gian_s": round(time.time() - t0, 1)}
    KQ["v12"] = V3.tong_ket_v12(KQ["ket_qua"])
    json.dump(KQ, open(kq_tep, "w"), ensure_ascii=False, indent=1, default=kiem2._json_default)
    print(json.dumps(KQ["v12"], ensure_ascii=False))
    print(json.dumps(KQ["do_phu_theo_luat"], ensure_ascii=False)); print(KQ["tong"])


if __name__ == "__main__":
    chi = [x for x in os.environ.get("CHI_LOAI", "").split(",") if x] or None
    chay(chi=chi)
