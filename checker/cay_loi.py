# -*- coding: utf-8 -*-
"""
cay_loi.py (bản 2) — kiểm thử đột biến cho kiem2.py.

Trên tế bào gốc thật (D1_td_1001, D4A2_td_1001), cấy TỪNG lỗi một lần, chạy kiem2, so với bản gốc.
  * Đột biến "phát hiện" CHỈ KHI luật ĐÍCH (luật bị nhắm) có ca TRƯỢT MỚI nhắc tới bộ phận bị cấy lỗi,
    hoặc luật đích chuyển sang TRƯỢT. Các luật khác trượt theo chỉ ghi thông tin.
  * Đối chứng (không lỗi), hai kiểu:
      - "trung_tinh": phép biến đổi không được làm thay đổi gì → MỌI ca trượt mới ở BẤT KỲ luật nào là báo
        động giả (trừ ngoại lệ ghi rõ, vd M.d khi nhân thành vỏ một lớp vì M.d trở nên kiểm được);
      - "sat_nguong_dat": đặt sát ngưỡng phía ĐẠT (khe 0,31; rộng 0,51; dài 1,99) → báo động giả nếu luật
        đích có ca trượt mới nhắc tới bộ phận đó.
  * Mọi đối chứng dựng từ tế bào thật (nhân thật khoan hầm, nhân thật khoét rỗng, ti thể thật…).
Ghi: dot_bien/<te_bao>__<ten>/, cay_loi_ket_qua.json; in tỉ lệ theo loại và độ phủ theo luật.
"""
import os, sys, json, copy, time, shutil, collections
import numpy as np
import trimesh
import manifold3d as mf
import kiem2
import kiem3 as _K3M, v3_cay_loi_chung as V3, luat_v12 as LV


def _K3(d, im_lang=True):
    out = _K3M.kiem3(d, "dong_vat", im_lang=True)
    out["dem_trang_thai"] = {"su_that": out["dem_su_that"], "tham_so": out["dem_tham_so"]}
    return out

HERE = os.path.dirname(os.path.abspath(__file__))
DL = os.path.join(HERE, "..", "bo_kiem_v2_dulieu")
OUT = os.path.join(HERE, "dot_bien")


def doc(te_bao):
    z = np.load(os.path.join(DL, te_bao, "mesh.npz"))
    L = json.load(open(os.path.join(DL, te_bao, "loai.json")))
    M = {}
    for k in z.files:
        n, s = k.rsplit("__", 1)
        M.setdefault(n, {})[s] = np.array(z[k])
    return M, L


def ghi(M, L, thu_muc):
    os.makedirs(thu_muc, exist_ok=True)
    arr = {}
    for n, d in M.items():
        arr[n + "__v"] = np.asarray(d["v"], float); arr[n + "__f"] = np.asarray(d["f"], np.int64)
    np.savez(os.path.join(thu_muc, "mesh.npz"), **arr)
    json.dump(L, open(os.path.join(thu_muc, "loai.json"), "w"), ensure_ascii=False, indent=1)


def theo_loai(L, loai):
    return sorted([n for n, t in L["loai"].items() if t == loai], key=lambda s: (len(s), s))


def mot(L, loai):
    return theo_loai(L, loai)[0]


def luoi(M, n):
    return kiem2.Luoi(M[n]["v"], M[n]["f"])


def man(M, n):
    return mf.Manifold(mf.Mesh64(np.ascontiguousarray(M[n]["v"], float), np.ascontiguousarray(M[n]["f"], np.uint64)))


def dat_man(M, n, m):
    mm = m.to_mesh64()
    M[n]["v"] = np.asarray(mm.vert_properties)[:, :3].copy(); M[n]["f"] = np.asarray(mm.tri_verts).astype(np.int64)


def kc_gan(M, a, b):
    A = luoi(M, a); B = luoi(M, b)
    d1, _, c1 = B.aabb().squared_distance(B.V, B.F, A.V)
    d2, _, c2 = A.aabb().squared_distance(A.V, A.F, B.V)
    i, j = int(d1.argmin()), int(d2.argmin())
    if d1[i] <= d2[j]:
        return A.V[i], c1[i], float(np.sqrt(d1[i]))
    return c2[j], B.V[j], float(np.sqrt(d2[j]))


def kc_dung(M, a, b):
    """Khe CHÍNH XÁC (dùng chính hàm của kiem2) giữa a và b (≤ 5 µm)."""
    return kiem2.kc_mat_mat(luoi(M, a), luoi(M, b), 5.0)[0]


def truc(v):
    c = v.mean(0); w, U = np.linalg.eigh(np.cov((v - c).T))
    return c, U[:, ::-1]


def dich(M, n, t):
    M[n]["v"] = M[n]["v"] + np.asarray(t)[None]


def tam(M, n):
    return trimesh.Trimesh(M[n]["v"], M[n]["f"], process=False).center_mass


def cl_khe(M, L, a, b, khe):
    """Dịch a theo hướng cặp điểm gần nhất để khe tới b = khe (âm = chồng). Lặp 3 lần cho đúng số."""
    for _ in range(4):
        p, q, d = kc_gan(M, a, b); u = (q - p) / max(d, 1e-12)
        dd = kc_dung(M, a, b) if khe >= 0 else d
        if khe >= 0 and abs(dd - khe) < 2e-4: break
        dich(M, a, u * ((dd if khe >= 0 else d) - khe))
        if khe < 0: break
    return "dich %s de khe toi %s = %.3f µm (do lai: %s)" % (a, b, khe, round(kc_dung(M, a, b), 4) if khe >= 0 else "chong")


def tia_toi_mat(M, ten_mat, p, u):
    Lm = luoi(M, ten_mat)
    _, t, _ = Lm.aabb().intersect_ray_first(Lm.V, Lm.F, p[None], u[None])
    return float(t[0])


def cl_day_ra(M, L, a, ten_mat, vuot, tam_tu=None):
    """Dịch a theo hướng từ tâm 'tam_tu' ra ngoài cho tới khi điểm xa nhất vượt mặt 'ten_mat' một đoạn vuot
    (vuot < 0: còn cách mặt |vuot| theo tia)."""
    c0 = tam(M, tam_tu or ten_mat); ca = M[a]["v"].mean(0); u = (ca - c0) / np.linalg.norm(ca - c0)
    v = M[a]["v"]; p = v[np.argmax(v @ u)]
    s = tia_toi_mat(M, ten_mat, p, u)
    dich(M, a, u * (s + vuot))
    return "day %s theo huong ra ngoai, diem xa nhat %s mat %s %.3f µm" % (a, "vuot" if vuot >= 0 else "con cach", ten_mat, abs(vuot))


def cl_khe_chua(M, L, a, ten_mat, khe):
    """Đặt a bên trong ten_mat sao cho khe chính xác tới mặt = khe (lặp)."""
    c0 = tam(M, ten_mat); ca = M[a]["v"].mean(0); u = (ca - c0) / np.linalg.norm(ca - c0)
    for _ in range(6):
        d = kc_dung(M, a, ten_mat)
        if abs(d - khe) < 2e-4: break
        dich(M, a, u * (d - khe))
    return "dat %s trong %s, khe toi mat = %.3f (do lai %.4f)" % (a, ten_mat, khe, kc_dung(M, a, ten_mat))


def cl_co_gian_truc(M, a, k_truc=1.0, k_ngang=1.0):
    v = M[a]["v"]; c, U = truc(v); p = (v - c) @ U
    p = p * np.array([k_truc, k_ngang, k_ngang])
    M[a]["v"] = c + p @ U.T


def cl_dai(M, L, a, dai):
    for _ in range(8):
        e = kiem2.BoPhan(a, "x", M[a]["v"], M[a]["f"]).kich_thuoc()[0]
        if abs(e - dai) < 5e-4: break
        cl_co_gian_truc(M, a, k_truc=dai / e)
    return "co gian %s doc truc dai toi %.3f µm (do lai %.4f)" % (a, dai, kiem2.BoPhan(a, "x", M[a]["v"], M[a]["f"]).kich_thuoc()[0])


def cl_rong(M, L, a, rong):
    for _ in range(8):
        w = kiem2.BoPhan(a, "x", M[a]["v"], M[a]["f"]).do_day()[0]
        if abs(w - rong) < 1e-3: break
        cl_co_gian_truc(M, a, k_ngang=rong / w)
    return "co gian ngang %s toi do day %.3f µm (do lai %.4f)" % (a, rong, kiem2.BoPhan(a, "x", M[a]["v"], M[a]["f"]).do_day()[0])


def cl_phong(M, L, a, k, tam_=None):
    c = tam(M, a) if tam_ is None else tam_
    M[a]["v"] = c + (M[a]["v"] - c) * k
    return "phong %s he so %s" % (a, k)


def cl_xoa(M, L, a, ca_du_kien=False):
    del M[a]; L["loai"].pop(a, None)
    if ca_du_kien: L["bo_phan_du_kien"].pop(a, None)
    return "xoa %s%s" % (a, " (ca khoi du kien)" if ca_du_kien else " (con trong du kien → mat)")


def cl_cat(M, L, a, so_manh=2, ho=0.05):
    m = man(M, a); c, U = truc(M[a]["v"]); nrm = U[:, 0]
    p = (M[a]["v"] - c) @ nrm
    cuts = np.linspace(p.min(), p.max(), so_manh + 1)[1:-1][::-1]
    parts = []; rest = m
    for k, t in enumerate(cuts):
        tren, duoi = rest.split_by_plane(tuple(nrm.tolist()), float(c @ nrm + t))
        parts.append(tren); rest = duoi
    parts.append(rest)
    out = parts[0]
    for k, q in enumerate(parts[1:]):
        out = out + q.translate(tuple((-(k + 1) * ho * nrm).tolist()))
    dat_man(M, a, out)
    return "cat %s thanh %d manh, tach %.2f µm (cung ten)" % (a, so_manh, ho)


def diem_tu_do(M, L, trong, tranh, ban_kinh):
    """Điểm trong mặt 'trong' xa mặt đó và xa các bộ phận 'tranh' nhất."""
    Lt = trong if isinstance(trong, kiem2.Luoi) else luoi(M, trong)
    xs = [np.linspace(Lt.lo[i], Lt.hi[i], 26) for i in range(3)]
    G = np.stack(np.meshgrid(*xs, indexing='ij'), -1).reshape(-1, 3)
    G = G[Lt.quan(G) > 0.5]
    d = np.sqrt(Lt.kc2_diem(G))
    for n in tranh:
        d = np.minimum(d, np.sqrt(luoi(M, n).kc2_diem(G)))
    i = int(d.argmax()); return G[i], float(d[i])


def cl_vao_nhan(M, L, a):
    nh = mot(L, "nhan"); p, dd = diem_tu_do(M, L, nh, theo_loai(L, "nhan_con"), 0)
    dich(M, a, p - M[a]["v"].mean(0))
    return "dat %s vao giua nhan (cach mat/nhan con %.2f µm truoc khi tru ban kinh)" % (a, dd)


def lam_nhan_rong(M, L, day=0.4, r_lo=0.3):
    nh = mot(L, "nhan"); v, f = M[nh]["v"], M[nh]["f"]
    N = man(M, nh); c = tam(M, nh)
    R = float(np.linalg.norm(v - c, axis=1).min()); k = (R - day) / R
    inner = N.translate(tuple((-c).tolist())).scale((k, k, k)).translate(tuple(c.tolist()))
    vo = N - inner
    H = float(np.linalg.norm(v - c, axis=1).max()) + 1
    for rot in [(0, 0, 0), (180, 0, 0), (0, 90, 0), (0, -90, 0), (-90, 0, 0), (90, 0, 0)]:
        vo = vo - mf.Manifold.cylinder(H, r_lo, r_lo, 24).rotate(rot).translate(tuple(c.tolist()))
    dat_man(M, nh, vo)
    im = inner.to_mesh64()
    return kiem2.Luoi(np.asarray(im.vert_properties)[:, :3], np.asarray(im.tri_verts).astype(np.int64)), vo.genus()


def cl_khoang(M, L, a):
    Lk, g = lam_nhan_rong(M, L)
    p, dd = diem_tu_do(M, L, Lk, theo_loai(L, "nhan_con"), 0)
    dich(M, a, p - M[a]["v"].mean(0))
    return "nhan thanh vo rong 0,4 µm, 6 lo (giong %d); dat %s vao khoang" % (g, a)


def dc_nhan_rong(M, L):
    Lk, g = lam_nhan_rong(M, L)
    return "nhan THAT khoet rong (vo 0,4 µm, 6 lo, giong %d), khoang chi chua nhan con" % g


def dc_khoan_ham(M, L, r=0.4):
    """Nhân thật khoan một đường hầm (giống 1, KHÔNG khoang kín), tránh nhân con ≥ r + 0.45 µm."""
    nh = mot(L, "nhan"); N = man(M, nh); c = tam(M, nh)
    NC = [luoi(M, n) for n in theo_loai(L, "nhan_con")]
    H = float(np.linalg.norm(M[nh]["v"] - c, axis=1).max()) * 2 + 2
    best = None
    for rot, ax in [((0, 0, 0), np.array([0, 0, 1.])), ((0, 90, 0), np.array([1., 0, 0])), ((-90, 0, 0), np.array([0, 1., 0]))]:
        b1 = np.cross(ax, [1, 0, 0] if abs(ax[0]) < 0.9 else [0, 1, 0]); b1 /= np.linalg.norm(b1); b2 = np.cross(ax, b1)
        for x in np.linspace(-2.5, 2.5, 11):
            for y in np.linspace(-2.5, 2.5, 11):
                o = c + x * b1 + y * b2
                cl = min(float(np.min(np.linalg.norm(np.cross(L_.V - o, ax), axis=1))) for L_ in NC) if NC else 9
                if best is None or cl > best[0]: best = (cl, rot, o, ax)
    cl, rot, o, ax = best
    cy = mf.Manifold.cylinder(H, r, r, 32, True).rotate(rot).translate(tuple(o.tolist()))
    out = N - cy
    dat_man(M, nh, out)
    return "nhan THAT khoan ham r=%.2f (giong %d, khong khoang), ham cach nhan con %.2f µm" % (r, out.genus(), cl - r)


def cl_er_tach(M, L, khe=0.02):
    """Bỏ lớp ER sát nhân: ER − (nhân nở theo pháp tuyến đỉnh khe+0,005 µm)."""
    nh = mot(L, "nhan"); er = mot(L, "luoi_noi_chat")
    tmn = trimesh.Trimesh(M[nh]["v"], M[nh]["f"], process=False)
    V2 = tmn.vertices + (khe + 0.005) * tmn.vertex_normals
    No = mf.Manifold(mf.Mesh64(np.ascontiguousarray(V2, float), np.ascontiguousarray(tmn.faces, np.uint64)))
    dat_man(M, er, man(M, er) - No)
    return "ER tach khoi nhan ~%.3f µm (do lai %.4f)" % (khe, kc_dung(M, er, nh))


def cl_vo_trung_nguoc(M, L, a):
    v, f = M[a]["v"], M[a]["f"]
    M[a]["v"] = np.vstack([v, v]); M[a]["f"] = np.vstack([f, f[:, ::-1] + len(v)])
    return "them vo trung lat nguoc vao %s" % a


def cl_tu_cat(M, L, a):
    v, f = M[a]["v"], M[a]["f"]; c, U = truc(v)
    M[a]["v"] = np.vstack([v, v + 0.15 * U[:, 0]]); M[a]["f"] = np.vstack([f, f + len(v)])
    return "chong hai ban sao lech 0,15 µm trong cung %s (tu cat)" % a


def cl_them_nhan(M, L):
    nc = theo_loai(L, "nhan_con")[0]
    M["nhan_2"] = copy.deepcopy(M[nc]); L["loai"]["nhan_2"] = "nhan"; L["bo_phan_du_kien"]["nhan_2"] = "nhan"
    dich(M, "nhan_2", np.zeros(3))
    return "them bo phan nhan thu hai (ban sao %s)" % nc


def cl_bo_nhieu(M, L, loai, con_lai):
    xs = theo_loai(L, loai)
    for n in xs[con_lai:]:
        cl_xoa(M, L, n, ca_du_kien=True)
    return "bo bot %s, con %d (ca khoi du kien)" % (loai, con_lai)


def cl_dep_mang(M, L, kz):
    mg = mot(L, "mang"); c = tam(M, mg); v = M[mg]["v"] - c; v[:, 2] *= kz; M[mg]["v"] = v + c
    return "nen mang theo z he so %s" % kz


def cl_trung_the_gan_mang(M, L):
    t = mot(L, "trung_the"); mg = mot(L, "mang"); nh = mot(L, "nhan")
    cn = tam(M, nh); cm = tam(M, mg); u = (cn - cm); u = -u / (np.linalg.norm(u) + 1e-9) if np.linalg.norm(u) > 1e-3 else np.array([0, 0, 1.])
    dich(M, t, cm + u * 8.0 - M[t]["v"].mean(0))
    return cl_khe_chua(M, L, t, mg, 0.7) + " (phia xa nhan)"


def ten_dat(base, rid):
    r = [x for x in base["luat"] if x["id"] == rid][0]
    hong = set(json.dumps(c["bo_phan"], ensure_ascii=False) for c in r["ca_truot"])
    return hong


def cap_gan(M, L, la, lb, n, bo=()):
    cs = []
    for a in theo_loai(L, la):
        if a in bo: continue
        ca = M[a]["v"].mean(0)
        for b in theo_loai(L, lb):
            if b == a: continue
            cs.append((np.linalg.norm(M[b]["v"].mean(0) - ca), a, b))
    cs.sort(); out = []; dung = set()
    for _, a, b in cs:
        if a in dung or b in dung: continue
        out.append((a, b)); dung |= {a, b}
        if len(out) >= n: break
    return out


def dat_het(base, rids, n):
    return not any(bp_truot_o_goc(base, rid, n) for rid in rids)


def ong_cong(r, R_, goc_do, n=24):
    th = np.radians(goc_do); ts = np.linspace(-th / 2, th / 2, n)
    sp = [mf.Manifold.sphere(r, 32).translate((R_ * np.sin(t), R_ * (1 - np.cos(t)), 0.0)) for t in ts]
    out = None
    for a, b in zip(sp[:-1], sp[1:]):
        hs = mf.Manifold.batch_hull([a, b]); out = hs if out is None else out + hs
    return out


def cl_ong_cong(M, L, a, dai_that, r=0.3, goc=150):
    """Thay ti thể a bằng ống cong (cung goc°) có dài THẬT theo trục = dai_that (gồm hai chỏm)."""
    arc = dai_that - 2 * r; m = ong_cong(r, arc / np.radians(goc), goc)
    c0 = M[a]["v"].mean(0); dat_man(M, a, m); dich(M, a, c0 - M[a]["v"].mean(0))
    return "thay %s bang ong cong %d°, r=%.2f, dai that theo truc %.2f µm" % (a, goc, r, dai_that)


def cl_chong_cung_loai(M, L, a, b, ti=0.7):
    """Đặt tâm a cách tâm b đúng ti·(ra + rb) (r = nửa cạnh OBB lớn nhất) → chắc chắn chồng nhau."""
    ca = M[a]["v"].mean(0); cb = M[b]["v"].mean(0)
    ra = kiem2.BoPhan("x", "x", M[a]["v"], M[a]["f"]).kich_thuoc()[0] / 2
    rb = kiem2.BoPhan("x", "x", M[b]["v"], M[b]["f"]).kich_thuoc()[0] / 2
    u = (ca - cb) / (np.linalg.norm(ca - cb) + 1e-12)
    dich(M, a, cb + u * ti * (ra + rb) - ca)
    return "dat tam %s cach tam %s %.2f·(ra+rb) de chong cung loai" % (a, b, ti)


def ke_hoach(te_bao, M, L, base):
    tt = theo_loai(L, "ti_the"); ly0 = theo_loai(L, "lysosome"); nc0 = theo_loai(L, "nhan_con")
    ok_kt = [n for n in tt if dat_het(base, ["4.3b", "4.3d", "M.f"], n)]
    ly_ok = [n for n in ly0 if dat_het(base, ["4.6a", "4.6d", "M.c", "M.a", "M.b", "M.f"], n)]
    ly = ly_ok + [n for n in ly0 if n not in ly_ok]
    nc = [n for n in nc0 if dat_het(base, ["4.2a", "4.2b"], n)] + [n for n in nc0 if not dat_het(base, ["4.2a", "4.2b"], n)]
    tt = [n for n in tt if dat_het(base, ["4.3a", "4.3d", "M.a", "M.b", "M.f"], n)] + [n for n in tt if not dat_het(base, ["4.3a", "4.3d", "M.a", "M.b", "M.f"], n)]
    nh = mot(L, "nhan"); mg = mot(L, "mang"); er = mot(L, "luoi_noi_chat"); go = mot(L, "golgi"); tr = mot(L, "trung_the")
    P = []
    def A(loai, fn, dich_, bp, kieu="dot_bien", ngoai_le=()):
        P.append({"loai": loai, "fn": fn, "dich": dich_, "bp": bp, "kieu": kieu, "ngoai_le": list(ngoai_le)})
    caps = cap_gan(M, L, "ti_the", "ti_the", 9)
    for a, b in caps[:3]:
        A("chong_lan_0_1", lambda M, L, a=a, b=b: cl_khe(M, L, a, b, -0.1), ["4.3d"], a)
    for a, b in caps[3:5]:
        A("khe_0_2", lambda M, L, a=a, b=b: cl_khe(M, L, a, b, 0.2), ["4.3d"], a)
    for a, b in caps[5:7]:
        A("khe_0_29_sat_nguong", lambda M, L, a=a, b=b: cl_khe(M, L, a, b, 0.29), ["4.3d"], a)
    for a, b in (caps[7:9] if len(caps) >= 9 else caps[:2]):
        A("khe_0_31_doi_chung", lambda M, L, a=a, b=b: cl_khe(M, L, a, b, 0.31), ["4.3d"], a, "sat_nguong_dat")
    for a, b in cap_gan(M, L, "lysosome", "ti_the", 1):
        A("khe_0_2_lysosome_ti_the", lambda M, L, a=a, b=b: cl_khe(M, L, a, b, 0.2), ["M.c"], a)
    A("ra_ngoai_mang_0_2", lambda M, L, a=tt[1]: cl_day_ra(M, L, a, mg, 0.2), ["M.a"], tt[1])
    A("ra_ngoai_mang_0_2", lambda M, L, a=ly[1]: cl_day_ra(M, L, a, mg, 0.2), ["M.a"], ly[1])
    A("khe_mang_0_4", lambda M, L, a=ly[2]: cl_khe_chua(M, L, a, mg, 0.4), ["M.b"], ly[2])
    ok_ts = [n for n in ok_kt if not bp_truot_o_goc(base, "4.3b@tham_so", n)] or ok_kt
    ts_ = lambda i: ok_kt[i % len(ok_kt)] if ok_kt[i % len(ok_kt)] in ok_ts else ok_ts[i % len(ok_ts)]
    A("dai_2_3", lambda M, L, a=ts_(0): cl_dai(M, L, a, 2.3), ["4.3b"], ts_(0), "st_dat_ts_truot")
    A("dai_2_01_sat_nguong", lambda M, L, a=ts_(1): cl_dai(M, L, a, 2.01), ["4.3b"], ts_(1), "st_dat_ts_truot")
    A("dai_1_99_doi_chung", lambda M, L, a=ok_kt[2]: cl_dai(M, L, a, 1.99), ["4.3b"], ok_kt[2], "sat_nguong_dat")
    A("rong_0_4_doi_chung", lambda M, L, a=ok_kt[3]: cl_rong(M, L, a, 0.4), ["4.3b"], ok_kt[3], "sat_nguong_dat"); P[-1]["ten_cu"] = "rong_0_4"
    A("rong_0_49_doi_chung", lambda M, L, a=ok_kt[4]: cl_rong(M, L, a, 0.49), ["4.3b"], ok_kt[4], "sat_nguong_dat"); P[-1]["ten_cu"] = "rong_0_49_sat_nguong"
    A("rong_0_51_doi_chung", lambda M, L, a=ok_kt[5 % len(ok_kt)]: cl_rong(M, L, a, 0.51), ["4.3b"], ok_kt[5 % len(ok_kt)], "sat_nguong_dat")
    A("xoa_ti_the", lambda M, L, a=tt[2]: cl_xoa(M, L, a), ["M.f"], tt[2])
    A("xoa_lysosome", lambda M, L, a=ly[3]: cl_xoa(M, L, a), ["M.f"], ly[3])
    A("xoa_trung_the", lambda M, L: cl_xoa(M, L, tr), ["4.8c"], tr)
    A("tach_hai_manh", lambda M, L, a=tt[3]: cl_cat(M, L, a), ["M.f"], tt[3])
    A("tach_hai_manh", lambda M, L, a=ly[4]: cl_cat(M, L, a), ["M.f"], ly[4])
    A("vo_trung_lat_nguoc", lambda M, L, a=ly[0]: cl_vo_trung_nguoc(M, L, a), ["M.f"], ly[0])
    A("tu_cat", lambda M, L, a=tt[6]: cl_tu_cat(M, L, a), ["M.f"], tt[6])
    A("lysosome_vao_nhan", lambda M, L, a=ly[0]: cl_vao_nhan(M, L, a), ["M.c"], ly[0])
    A("ti_the_vao_nhan", lambda M, L, a=tt[7]: cl_vao_nhan(M, L, a), ["4.3a"], tt[7])
    A("vao_khoang_nhan_rong", lambda M, L, a=ly[1]: cl_khoang(M, L, a), ["M.c"], ly[1])
    A("vao_khoang_nhan_rong", lambda M, L, a=tt[8]: cl_khoang(M, L, a), ["4.3a"], tt[8])
    A("nhan_con_bang_0", lambda M, L: cl_bo_nhieu(M, L, "nhan_con", 0), ["4.2a"], "nhan_con")
    A("nhan_con_lo_0_01", lambda M, L: cl_day_ra(M, L, nc[0], nh, 0.01, tam_tu=nh), ["4.2a"], nc[0])
    A("nhan_con_khe_0_2", lambda M, L: cl_khe_chua(M, L, nc[1], nh, 0.2), ["4.2b"], nc[1])
    A("nhan_con_rong_0_4", lambda M, L: cl_phong(M, L, nc[0], 0.4 / kiem2.BoPhan("x", "x", M[nc[0]]["v"], M[nc[0]]["f"]).kich_thuoc()[0]), ["4.2a"], nc[0])
    A("bo_nhan", lambda M, L: cl_xoa(M, L, nh, ca_du_kien=True), ["1.1"], "nhan")
    A("mang_phong_1_6", lambda M, L: cl_phong(M, L, mg, 1.6), ["1.2a"], mg, "st_dat_ts_truot")
    A("mang_dep_0_7", lambda M, L: cl_dep_mang(M, L, 0.7), ["1.2b"], mg)
    A("hai_nhan", lambda M, L: cl_them_nhan(M, L), ["1.3"], "nhan")
    A("nhan_ra_ngoai_mang", lambda M, L: cl_day_ra(M, L, nh, mg, 0.2), ["1.4a"], nh)
    A("nhan_khe_mang_0_8", lambda M, L: cl_khe_chua(M, L, nh, mg, 0.8), ["1.4b"], nh)
    A("nhan_thu_0_5", lambda M, L: cl_phong(M, L, nh, 0.5), ["1.5"], nh)
    A("er_tach_nhan_0_02", lambda M, L: cl_er_tach(M, L, 0.02), ["4.4a"], er)
    A("er_phong_1_25", lambda M, L: cl_phong(M, L, er, 1.25, tam_=tam(M, nh)), ["4.4b"], er)
    A("er_lot_vao_nhan", lambda M, L: cl_phong(M, L, er, 0.9, tam_=tam(M, nh)), ["4.4b"], er)
    for rl in (0.9, 1.2, 1.5):
        A("nhan_rong_lo_%.1f_lysosome" % rl, lambda M, L, rl=rl: (lambda Lk_g: (dich(M, ly[2], diem_tu_do(M, L, Lk_g[0], theo_loai(L, "nhan_con"), 0)[0] - M[ly[2]]["v"].mean(0)), "vo 0,4 µm, 6 lo r=%.1f (giong %d), dat %s vao khoang" % (rl, Lk_g[1], ly[2]))[1])(lam_nhan_rong(M, L, 0.4, rl)), ["M.c"], ly[2])
    cl2 = cap_gan(M, L, "lysosome", "lysosome", 1, bo=[n for n in ly0 if not dat_het(base, ["chong_lan_cung_loai"], n)])
    for a, b in cl2:
        A("chong_cung_loai_lysosome", lambda M, L, a=a, b=b: cl_chong_cung_loai(M, L, a, b, 0.75), ["M.h"], a)
    ncf = [n for n in nc0 if dat_het(base, ["chong_lan_cung_loai"], n)]
    if ncf and len(nc0) >= 2:
        b_ = [n for n in nc0 if n != ncf[0]][0]
        A("chong_cung_loai_nhan_con", lambda M, L, a=ncf[0], b=b_: cl_chong_cung_loai(M, L, a, b, 0.6), ["M.h"], ncf[0])
    A("ti_the_cong_dai_2_4", lambda M, L, a=ts_(6): cl_ong_cong(M, L, a, 2.4), ["4.3b"], ts_(6), "st_dat_ts_truot")
    A("ti_the_cong_dai_1_8_doi_chung", lambda M, L, a=ok_kt[7 % len(ok_kt)]: cl_ong_cong(M, L, a, 1.8), ["4.3b"], ok_kt[7 % len(ok_kt)], "sat_nguong_dat")
    A("golgi_khe_nhan_0_2", lambda M, L: cl_khe(M, L, go, nh, 0.2), ["4.5c"], go)
    A("lysosome_1_5um", lambda M, L, a=ly_ok[0]: cl_phong(M, L, a, 1.5 / kiem2.BoPhan("x", "x", M[a]["v"], M[a]["f"]).kich_thuoc()[0]), ["4.6a"], ly_ok[0])
    A("lysosome_keo_dai", lambda M, L, a=ly_ok[1]: (cl_co_gian_truc(M, a, k_truc=1.5), "keo dai %s 1,5 lan" % a)[1], ["4.6d"], ly_ok[1])
    A("trung_the_3_manh", lambda M, L: cl_cat(M, L, tr, 3), ["4.8a"], tr)
    A("trung_the_gan_mang", lambda M, L: cl_trung_the_gan_mang(M, L), ["4.8b"], tr)
    A("trung_the_khe_nhan_0_2", lambda M, L: cl_khe(M, L, tr, nh, 0.2), ["4.8c"], tr)
    A("nhan_vo_mot_lop", lambda M, L: dc_nhan_rong(M, L), ["M.d"], nh)
    if len(tt) <= 20:
        A("ti_the_con_9", lambda M, L: cl_bo_nhieu(M, L, "ti_the", 9), ["4.3c"], "ti_the")
    if len(ly) <= 8:
        A("lysosome_con_2", lambda M, L: cl_bo_nhieu(M, L, "lysosome", 2), ["4.6c"], "lysosome")
    A("dc_nhan_khoet_rong", lambda M, L: dc_nhan_rong(M, L), [], nh, "trung_tinh", ngoai_le=["M.d"])
    A("dc_nhan_khoan_ham", lambda M, L: dc_khoan_ham(M, L), [], nh, "trung_tinh")
    A("dc_dich_0_01", lambda M, L, a=tt[5]: (dich(M, a, 0.01 * (M[a]["v"].mean(0) - tam(M, nh)) / np.linalg.norm(M[a]["v"].mean(0) - tam(M, nh))), "dich %s 0,01 µm ra xa nhan" % a)[1], [], tt[5], "trung_tinh")
    A("dc_chia_tam_giac", lambda M, L, a=tt[6]: (M[a].update(dict(zip(["v", "f"], (lambda m: (np.asarray(m.vertices), np.asarray(m.faces)))(trimesh.Trimesh(M[a]["v"], M[a]["f"], process=False).subdivide())))), "chia tam giac %s" % a)[1], [], tt[6], "trung_tinh")
    return P


def k_nv(c): return json.dumps(c.get("bo_phan"), ensure_ascii=False) + "|" + str(c.get("nguong"))
def k_v(c): return k_nv(c) + "|" + json.dumps(c.get("gia_tri"), ensure_ascii=False, default=str)


def truot_moi(base, mut, co_gia_tri):
    B = {r["id"]: r for r in base["luat"]}; out = {}
    k = k_v if co_gia_tri else k_nv
    for r in mut["luat"]:
        b = B[r["id"]]; cu = set(k(c) for c in b["ca_truot"])
        moi = [c for c in r["ca_truot"] if k(c) not in cu]
        doi = (b["trang_thai"] != "TRUOT" and r["trang_thai"] == "TRUOT")
        if moi or doi:
            out[r["id"]] = {"moi": moi, "chuyen_truot": doi}
    return out


def nhac(c, bp):
    """Ca nhắc ĐÚNG bộ phận bp (so khớp nguyên tên, không khớp chuỗi con)."""
    b = c.get("bo_phan")
    return b == bp or (isinstance(b, list) and bp in b)


def bp_truot_o_goc(base, rid, bp):
    if rid == "chong_lan_cung_loai":
        return any(bp in x["bo_phan"] for x in base["chong_lan_cung_loai"]["cap"])
    ids_ = {q["id"] for q in base["luat"]}
    if rid not in ids_:
        m_ = V3.map_dich([rid], "dv", ids_)
        m_ = [x for x in m_ if x.endswith("@su_that")] or m_
        return any(bp_truot_o_goc(base, x, bp) for x in m_)
    r = [q for q in base["luat"] if q["id"] == rid][0]
    return any(nhac(c, bp) for c in r["ca_truot"])


def cung_loai_moi(base, mut, bp):
    cu = set(json.dumps(x["bo_phan"]) for x in base["chong_lan_cung_loai"]["cap"])
    return [x for x in mut["chong_lan_cung_loai"]["cap"] if json.dumps(x["bo_phan"]) not in cu and bp in x["bo_phan"]]


def _them_nhan_con(M, L, so_them):
    nh = mot(L, "nhan"); nc = theo_loai(L, "nhan_con"); goc = nc[0]
    for i in range(so_them):
        p, d = diem_tu_do(M, L, nh, theo_loai(L, "nhan_con"), 0)
        V = M[goc]["v"]; c = V.mean(0); r0 = float(np.linalg.norm(V - c, axis=1).max())
        r = min(r0, 0.6 * d)
        if r < 0.3: raise RuntimeError("khong du cho cho nhan con moi (r = %.3f)" % r)
        ten = "nhan_con_v12_%d" % i
        M[ten] = {"v": p + (V - c) * (r / r0), "f": M[goc]["f"].copy()}; L["loai"][ten] = "nhan_con"; L["bo_phan_du_kien"][ten] = "nhan_con"
    return "them %d nhan con (tong %d)" % (so_them, len(theo_loai(L, "nhan_con")))


def _phong_den_canh(M, L, a, canh):
    e = kiem2.BoPhan("x", "x", M[a]["v"], M[a]["f"]).kich_thuoc()[0]
    return cl_phong(M, L, a, canh / e) + " (canh OBB lon %.3f -> %.3f)" % (e, canh)


def ke_hoach_v12(te_bao, M, L, base):
    P = []
    def A(loai, fn, dich_, bp, kieu):
        P.append({"loai": "v12_" + loai, "fn": fn, "dich": dich_, "bp": bp, "kieu": kieu, "ngoai_le": []})
    nc = theo_loai(L, "nhan_con"); ly = theo_loai(L, "lysosome"); nh = mot(L, "nhan"); go = mot(L, "golgi")
    a0, a1, _ = LV.KHOANG("dv", "4.2a", "nhân con", 0); p1 = LV.P("dv.so_nhan_con")[1]
    if len(nc) <= p1 and p1 + 1 <= a1:
        A("nhan_con_%d" % (p1 + 1), lambda M, L: _them_nhan_con(M, L, p1 + 1 - len(nc)), ["4.2a"], "nhan_con", "st_dat_ts_truot")
    if ly:
        rmax = LV.P("dv.dai_lysosome")[1]; l1 = LV.KHOANG("dv", "4.6a", "đường kính")[1]
        canh = min(l1 - 0.1, 2 * rmax + 0.2)
        A("lysosome_canh_%.2f" % canh, lambda M, L, a=ly[-1]: _phong_den_canh(M, L, a, canh), ["4.6a"], ly[-1], "st_dat_ts_truot")
    for a, b in cap_gan(M, L, "ti_the", "ti_the", 2)[:1]:
        A("ti_the_khe_0_2_M.h", lambda M, L, a=a, b=b: cl_khe(M, L, a, b, 0.2), ["M.h"], a, "st_dat_ts_truot")
    for a, b in cap_gan(M, L, "lysosome", "ti_the", 2)[1:2]:
        A("lysosome_ti_the_khe_0_01_M.c0", lambda M, L, a=a, b=b: cl_khe(M, L, a, b, 0.01), ["M.c₀@su_that", "M.c@tham_so"], a, "st_dat_ts_truot")
    for a, b in cap_gan(M, L, "lysosome", "ti_the", 3)[2:3]:
        A("lysosome_chong_ti_the_M.c0", lambda M, L, a=a, b=b: cl_khe(M, L, a, b, -0.1), ["M.c₀@su_that"], a, "dot_bien")
    if go and nh:
        A("golgi_chong_nhan_M.c0", lambda M, L: cl_khe(M, L, go, nh, -0.1), ["M.c₀@su_that"], go, "dot_bien")
    er = mot(L, "luoi_noi_chat")
    if er and nh:
        def er_vao(M, L):
            c = tam(M, nh); u = c - M[er]["v"].mean(0); u /= max(np.linalg.norm(u), 1e-12)
            dich(M, er, u * 0.15); return "dich ER 0,15 µm vao phia tam nhan (cap ER–nhan truoc day mien khe)"
        A("er_lan_vao_nhan_M.c0", er_vao, ["M.c₀@su_that"], er, "dot_bien")
    if er and nh:
        for dl, kieu in ((0.05, "dot_bien"), (0.02, "dot_bien"), (0.005, "sat_nguong_dat")):
            A("er_loi_vao_nhan_%g_M.c0" % dl, lambda M, L, dl=dl: _er_loi(M, er, nh, dl), ["M.c₀@su_that"], er, kieu)
    for a, b in cap_gan(M, L, "lysosome", "ti_the", 4)[3:4]:
        A("lysosome_chong_ti_the_0.02_M.c0", lambda M, L, a=a, b=b: cl_khe(M, L, a, b, -0.02), ["M.c₀@su_that"], a, "dot_bien")
        A("lysosome_chong_ti_the_0.005_M.c0", lambda M, L, a=a, b=b: cl_khe(M, L, a, b, -0.005), ["M.c₀@su_that"], a, "sat_nguong_dat")
    mg = mot(L, "mang"); tts = theo_loai(L, "trung_the"); tt_ = tts[0] if tts else None
    if go and nh and mg and tt_ and er:
        S45 = {}
        def g45(M, L):
            return _g45_goc(M, L, go, tt_, nh, mg, er, S45)
        def pha45_khe(M, L):
            return _g45_dich_ra(M, go, tt_, nh, mg, S45["u"], S45["v_muc"])
        def pha45_tam(M, L):
            dich(M, tt_, S45["tt_ngang"] - tam(M, tt_))
            return "doi trung the NGANG ra khoi bao loi Golgi (khe Golgi–trung the %.3f µm); Golgi giu nguyen" % kc_dung(M, go, tt_)
        P.append({"loai": "v12_golgi_lech_mang_4.5b", "ten_cu": "v12_golgi_sat_mang_4.5b", "fn": pha45_khe, "dich": ["4.5b@su_that"], "bp": go, "kieu": "dot_bien", "ngoai_le": [], "goc_rieng": ("g45", g45)})
        P.append({"loai": "v12_trung_the_ngang_khoi_bao_loi_4.5b", "ten_cu": "v12_trung_the_ra_khoi_bao_loi_4.5b", "fn": pha45_tam, "dich": ["4.5b@su_that"], "bp": [go, tt_], "kieu": "dot_bien", "ngoai_le": [], "goc_rieng": ("g45", g45)})
    if len(ly) >= 3 and nh and mg:
        def xa_nhan(M, n):
            return kc_dung(M, n, nh)
        def g46(M, L):
            ds = sorted(ly, key=lambda n: xa_nhan(M, n)); giu = ds[-1]
            for n in ds[:-1]: cl_khe(M, L, n, nh, 2.1)
            return "goc rieng 4.6b: moi lysosome tru %s dat khe nhan 2,1 µm" % giu
        def pha46_xa(M, L):
            for n in ly: cl_khe_chua(M, L, n, mg, 0.4)
            return "moi lysosome sat mang (khe 0,4 µm) → trung vi u > 0,5"
        def pha46_gan(M, L):
            for n in ly: cl_khe(M, L, n, nh, 2.1)
            return _doi_lysosome_u_nho(M, L, ly, nh, mg, 2.1, 0.45)
        P.append({"loai": "v12_lysosome_sat_mang_4.6b", "fn": pha46_xa, "dich": ["4.6b@su_that"], "bp": "lysosome", "kieu": "dot_bien", "ngoai_le": [], "goc_rieng": ("g46", g46)})
        P.append({"loai": "v12_lysosome_tat_ca_gan_nhan_4.6b", "fn": pha46_gan, "dich": ["4.6b@su_that"], "bp": "lysosome", "kieu": "dot_bien", "ngoai_le": [], "goc_rieng": ("g46", g46)})
    return P


def _g45_dich_ra(M, go, tt, nh, mg, u, v_muc):
    """Vế 1 của 4.5b: dời Golgi + trung thể theo u (tâm nhân → Golgi ở gốc) tới khi khe nhân − khe màng = v_muc (khe CHÍNH XÁC)."""
    s_tong = 0.0
    for _ in range(10):
        gn, gm = kc_dung(M, go, nh), kc_dung(M, go, mg); v = gn - gm
        if abs(v - v_muc) < 0.01:
            break
        s = (v_muc - v) / 2; dich(M, go, u * s); dich(M, tt, u * s); s_tong += s
    gn, gm = kc_dung(M, go, nh), kc_dung(M, go, mg)
    if not gn > gm:
        raise RuntimeError("khong tao duoc vi pham ve 1 (khe nhan %.4f, khe mang %.4f)" % (gn, gm))
    return "doi Golgi + trung the %.3f µm theo ban kinh ra phia mang: khe nhan %.4f > khe mang %.4f (vi pham ve 1 %.3f µm)" % (s_tong, gn, gm, gn - gm)


def _g45_goc(M, L, go, tt, nh, mg, er, S, n_thu=4000, seed=45, v_muc=0.3):
    """Gốc riêng 4.5b (15/11): Golgi ĐẠT cả hai vế, không chạm bộ phận nào; ghi vào S hướng u, vị trí ngang của trung thể.
    1) Thử n_thu phép xoay cứng Golgi quanh tâm nhân (hạt cố định), dời theo bán kính để khe đỉnh tới nhân = (xa nhất của ER) + k,
       k = 0,3 → 0,2 → 0,1 (vỏ ER đặc 0–2 µm quanh nhân); giữ chỗ khe màng − khe nhân ≥ 0,15 và Golgi trong màng, không chạm ER.
    2) Trung thể đặt TRONG Golgi, lệch về phía nhân (tâm cách mặt Golgi ≥ r + 0,05) → tâm trong bao lồi, không chạm ER/bộ phận khác.
       (Golgi là khối lồi liền: không có khoảng trống trong bao lồi nên cặp Golgi–trung thể ở gốc tất yếu chồng — M.c₀ trượt sẵn ở gốc.)
    3) Chọn chỗ ít ti thể/lysosome vướng nhất (khe < 0,35) với Golgi ở gốc, ở vị trí đột biến vế 1 (ước lượng) và với trung thể ở chỗ
       ngang của đột biến vế 2; dời các bào quan vướng tới chỗ trống gần chỗ cũ nhất (khe ≥ 0,35 tới mọi bộ phận, giữ số lượng).
    4) Kiểm lại bằng khe CHÍNH XÁC: khe nhân < khe màng và tâm trung thể trong bao lồi (như kiem3), không thì lỗi cấy."""
    from scipy.spatial.transform import Rotation as Rot
    from scipy.spatial import ConvexHull, Delaunay
    rng = np.random.default_rng(seed)
    cN = tam(M, nh); Ln = luoi(M, nh); Lm = luoi(M, mg); Le = luoi(M, er)
    dEmax = float(np.sqrt(Ln.kc2_diem(M[er]["v"])).max())
    G0 = M[go]["v"]; GF = M[go]["f"]; sub = rng.choice(len(G0), min(300, len(G0)), replace=False)
    T0 = M[tt]["v"]; cT = tam(M, tt); rT = float(np.linalg.norm(T0 - cT, axis=1).max())
    nho = [n for n, t in L["loai"].items() if t in ("ti_the", "lysosome") and n in M]
    BB = {n: (M[n]["v"].min(0), M[n]["v"].max(0)) for n in nho}; LL = {}
    def Lx(n):
        if n not in LL: LL[n] = luoi(M, n)
        return LL[n]
    def dmin(Lp, P):
        return float(np.sqrt(Lp.kc2_diem(P)).min())
    def dat(Rm, gn_muc):
        Gv = (G0 - cN) @ Rm.T + cN; u = Gv.mean(0) - cN; u /= np.linalg.norm(u)
        for _ in range(3):
            Gv = Gv + u * (gn_muc - dmin(Ln, Gv[sub]))
        return Gv, u
    def uoc_dich(Gv, u):
        s = 0.0
        for _ in range(5):
            Gs = Gv[sub] + u * s; s += (v_muc - (dmin(Ln, Gs) - dmin(Lm, Gs))) / 2
        return s
    def tt_trong(Gv, u):
        Lg = kiem2.Luoi(Gv, GF); c = Gv.mean(0); p = c
        for t in np.linspace(0, 2.0, 41):
            q = c - u * t
            if Lg.quan(q[None])[0] > 0.5 and np.sqrt(Lg.kc2_diem(q[None]))[0] >= rT + 0.05: p = q
            else: break
        return p
    def tt_ngang(Gv, u, p):
        c = Gv.mean(0); w_, U = np.linalg.eigh(np.cov((Gv - c).T)); ra = []
        for e in U.T[::-1][:2]:
            w = e - (e @ u) * u; w /= np.linalg.norm(w)
            for sg in (1, -1):
                ww = sg * w; q = p + ww * (float(((Gv - p) @ ww).max()) + rT + 0.35)
                ok = Lm.quan(q[None])[0] > 0.5 and np.sqrt(Lm.kc2_diem(q[None]))[0] >= rT + 0.35 and \
                    Le.quan(q[None])[0] < 0.5 and np.sqrt(Le.kc2_diem(q[None]))[0] >= rT + 0.35
                if ok: ra.append(q)
        return ra
    def vuong(Gvs, qs):
        out = set()
        for Gv in Gvs:
            Lg = kiem2.Luoi(Gv, GF); lo, hi = Gv.min(0) - 0.4, Gv.max(0) + 0.4
            for n in nho:
                a, b = BB[n]
                if np.any(b < lo) or np.any(a > hi) or n in out: continue
                if dmin(Lg, M[n]["v"]) < 0.35 or dmin(Lx(n), Gv) < 0.35 or (Lg.quan(M[n]["v"]) > 0.5).any() or (Lx(n).quan(Gv) > 0.5).any():
                    out.add(n)
        for q in qs:
            for n in nho:
                if n not in out and float(np.linalg.norm(M[n]["v"] - q, axis=1).min()) < rT + 0.35 + 0.05:
                    out.add(n)
        return out
    ung = []
    for k_ in (0.3, 0.2, 0.1):
        gn_muc = dEmax + k_
        for i in range(n_thu):
            Rm = Rot.random(random_state=int(rng.integers(1 << 30))).as_matrix(); Gv, u = dat(Rm, gn_muc)
            gm = dmin(Lm, Gv[sub])
            if gm - gn_muc >= 0.15 and (Lm.quan(Gv[sub]) > 0.5).all():
                ung.append((gm - gn_muc, i, Gv, u))
        if ung: break
    if not ung:
        raise RuntimeError("khong co cho dat Golgi voi khe nhan < khe mang ngoai vo ER")
    tot = None
    for dd, i, Gv, u in sorted(ung, key=lambda r: -r[0])[:40]:
        if (Le.quan(Gv) > 0.5).any() or (kiem2.Luoi(Gv, GF).quan(M[er]["v"]) > 0.5).any() or dmin(Le, Gv) < 0.05:
            continue
        s = uoc_dich(Gv, u); p = tt_trong(Gv, u); qs = tt_ngang(Gv, u, p)
        if not qs: continue
        bien = {i_: float(np.sqrt(Lm.kc2_diem(q[None]))[0] - np.sqrt(Ln.kc2_diem(q[None]))[0]) for i_, q in enumerate(qs)}
        qi = [i_ for i_ in bien if bien[i_] >= 0.2] or list(bien)
        best_q = qs[min(qi, key=lambda i_: (len(vuong([], [qs[i_]])), -bien[i_]))]
        vg = vuong([Gv, Gv + u * (0.5 * s), Gv + u * s, Gv + u * (s + 0.15)], [best_q])
        key = (len(vg), -dd)
        if tot is None or key < tot[0]:
            tot = (key, Gv, u, s, p, best_q, vg, dd, k_)
    if tot is None:
        raise RuntimeError("moi cho ung vien deu cham ER hoac khong co cho ngang cho trung the")
    _, Gv, u, s, p, q, vg, dd, k_ = tot
    M[go]["v"] = Gv; dich(M, tt, p - tam(M, tt))
    doi = []
    if vg:
        cam = [(Gv + u * t, GF) for t in (0.5 * s, s, s + 0.15)] + [(M[tt]["v"] - tam(M, tt) + q, M[tt]["f"])]
        parts = [(M[n]["v"], M[n]["f"]) for n in M if n != mg and n not in vg] + cam
        V_, F_, o_ = [], [], 0
        for Vn, Fn in parts:
            V_.append(np.asarray(Vn, float)); F_.append(np.asarray(Fn, np.int64) + o_); o_ += len(Vn)
        LU = kiem2.Luoi(np.concatenate(V_), np.concatenate(F_))
        xs = [np.arange(Lm.lo[j], Lm.hi[j], 0.3) for j in range(3)]
        Gd = np.stack(np.meshgrid(*xs, indexing="ij"), -1).reshape(-1, 3); Gd = Gd[Lm.quan(Gd) > 0.5]
        Dm = np.sqrt(Lm.kc2_diem(Gd)); keep = LU.quan(Gd) < 0.5; Gd, Dm = Gd[keep], Dm[keep]; D = np.sqrt(LU.kc2_diem(Gd))
        for n in sorted(vg):
            c = M[n]["v"].mean(0); r = float(np.linalg.norm(M[n]["v"] - c, axis=1).max())
            ok = np.where((D >= r + 0.35) & (Dm >= r + 0.35))[0]
            if len(ok) == 0:
                raise RuntimeError("khong con cho trong de doi %s" % n)
            j = ok[np.argmin(np.linalg.norm(Gd[ok] - c, axis=1))]; dich(M, n, Gd[j] - c); doi.append(n)
            Ln_ = luoi(M, n); D = np.minimum(D, np.sqrt(Ln_.kc2_diem(Gd))); inn = Ln_.quan(Gd) > 0.5; D[inn] = 0
    gn, gm = kc_dung(M, go, nh), kc_dung(M, go, mg)
    inh = bool(Delaunay(M[go]["v"][ConvexHull(M[go]["v"]).vertices]).find_simplex(np.array(trimesh.Trimesh(M[tt]["v"], M[tt]["f"], process=False).center_mass)[None])[0] >= 0)
    if not (gn < gm and inh):
        raise RuntimeError("goc g45 khong dat 4.5b: khe nhan %.4f, khe mang %.4f, tam trung the trong bao loi %s" % (gn, gm, inh))
    S.update({"u": u, "tt_ngang": q, "v_muc": v_muc, "s_uoc": s})
    return ("goc rieng 4.5b (15/11): xoay cung Golgi quanh tam nhan (hat %d, chon 1/%d ung vien, k = %.1f): khe nhan %.4f < khe mang %.4f, "
            "khe ER %.3f; trung the trong Golgi lech ve phia nhan (tam trong bao loi); doi %d bao quan vuong (khe >= 0,35): %s"
            % (seed, len(ung), k_, gn, gm, kc_dung(M, go, er), len(doi), ",".join(doi)))


def _u_lyso(M, nh, mg, P):
    Ln = luoi(M, nh); Lm = luoi(M, mg); P = np.atleast_2d(P)
    dn = np.sqrt(Ln.kc2_diem(P)); dm = np.sqrt(Lm.kc2_diem(P))
    return dn / (dn + dm + 1e-12)


def _doi_lysosome_u_nho(M, L, ly, nh, mg, khe, u_max):
    tm = trimesh.Trimesh(M[nh]["v"], M[nh]["f"], process=False); V = np.asarray(tm.vertices); N = np.asarray(tm.vertex_normals)
    Lm = luoi(M, mg); dung = [tam(M, n) for n in ly if _u_lyso(M, nh, mg, tam(M, n))[0] < u_max]; doi = 0
    for n in ly:
        c = tam(M, n)
        if _u_lyso(M, nh, mg, c)[0] < u_max:
            continue
        r = float(np.linalg.norm(M[n]["v"] - c, axis=1).max())
        Pc = V + N * (khe + r); uu = _u_lyso(M, nh, mg, Pc); trg = Lm.quan(Pc) > 0.5
        ok = False
        for i in np.argsort(uu):
            if uu[i] >= u_max: break
            if not trg[i] or any(np.linalg.norm(Pc[i] - q) < 2 * r + 0.2 for q in dung): continue
            dich(M, n, Pc[i] - tam(M, n)); cl_khe(M, L, n, nh, khe)
            c2 = tam(M, n)
            if _u_lyso(M, nh, mg, c2)[0] < u_max and Lm.quan(c2[None])[0] > 0.5:
                dung.append(c2); doi += 1; ok = True; break
        if not ok:
            raise RuntimeError("khong tim duoc cho u < %g cho %s" % (u_max, n))
    um = float(_u_lyso(M, nh, mg, np.array([tam(M, n) for n in ly])).max())
    if um >= 0.5:
        raise RuntimeError("sau khi doi van con u max = %.3f >= 0,5" % um)
    return "moi lysosome khe nhan %.1f µm; doi %d lysosome toi cho u < %g → u max = %.3f < 0,5 (khong co lysosome rai ra ngoai)" % (khe, doi, u_max, um)


def _er_loi(M, er, nh, dl, r=0.15):
    """Lồi cục bộ: đỉnh ER trong bán kính r quanh một điểm tiếp xúc ER–nhân dời vào phía tâm nhân dl µm."""
    Ln = luoi(M, nh); d2, _, _ = Ln.aabb().squared_distance(Ln.V, Ln.F, M[er]["v"])
    tx = np.where(d2 < 1e-8)[0]
    if len(tx) == 0:
        tx = np.argsort(d2)[:1]
    c0 = M[er]["v"][tx[len(tx) // 2]]; gan = np.linalg.norm(M[er]["v"] - c0, axis=1) < r
    n = c0 - M[nh]["v"].mean(0); n /= np.linalg.norm(n)
    V = M[er]["v"].copy(); V[gan] -= dl * n; M[er]["v"] = V
    return "loi ER %d dinh (r = %g µm) vao nhan %g µm" % (int(gan.sum()), r, dl)


def _don_thu_muc(d):
    """Xoá mesh.npz, kiem2_dang_chay.json, kiem2_bo_dem.pkl; nén gzip kiem2.json, kiem3.json (giữ loai.json, kiem3.csv)."""
    import gzip
    for f in ("mesh.npz", "kiem2_dang_chay.json", "kiem2_bo_dem.pkl"):
        try: os.remove(os.path.join(d, f))
        except FileNotFoundError: pass
    for f in ("kiem2.json", "kiem3.json"):
        p = os.path.join(d, f)
        if os.path.exists(p):
            with open(p, "rb") as a, gzip.open(p + ".gz", "wb") as b: shutil.copyfileobj(a, b)
            os.remove(p)


def chay(te_bao_list=("D1_td_1001", "D4A2_td_1001")):
    t0 = time.time()
    chi = [x for x in os.environ.get("CHI_LOAI", "").split(",") if x]
    if chi and os.path.exists(os.path.join(HERE, "cay_loi_ket_qua.json")):
        cu = json.load(open(os.path.join(HERE, "cay_loi_ket_qua.json")))
        KQ = {"ket_qua": [x for x in cu["ket_qua"] if x["loai"] not in chi], "goc": cu["goc"]}
    else:
        if os.path.isdir(OUT): shutil.rmtree(OUT)
        KQ = {"ket_qua": [], "goc": {}}
    os.makedirs(OUT, exist_ok=True)
    for te_bao in te_bao_list:
        M0, L0 = doc(te_bao)
        dg = os.path.join(OUT, te_bao + "__goc"); ghi(M0, L0, dg)
        base = _K3(dg)
        KQ["goc"][te_bao] = {"dem": base["dem_trang_thai"], "truot": [r["id"] for r in base["luat"] if r["trang_thai"] == "TRUOT"]}
        co_san = {q["id"] for q in base["luat"]}
        CHIA = os.environ.get("CHIA", "")
        ci, cn = (int(CHIA.split("/")[0]), int(CHIA.split("/")[1])) if CHIA else (0, 1)
        goc_r = {}; base_toan = base
        for k, x in enumerate(ke_hoach(te_bao, M0, L0, base_toan) + ke_hoach_v12(te_bao, M0, L0, base_toan)):
            if chi and x["loai"] not in chi:
                continue
            if k % cn != ci:
                continue
            x = dict(x, dich=V3.map_dich(x["dich"], "dv", co_san), ngoai_le=V3.map_dich(x.get("ngoai_le", []), "dv", co_san))
            base, M_src, L_src = base_toan, M0, L0
            if x.get("goc_rieng"):
                gten, gfn = x["goc_rieng"]
                if gten not in goc_r:
                    M1, L1 = copy.deepcopy(M0), copy.deepcopy(L0)
                    try:
                        mt = gfn(M1, L1); dgr = os.path.join(OUT, "%s__goc_rieng_%s" % (te_bao, gten)); ghi(M1, L1, dgr)
                        goc_r[gten] = (_K3(dgr), M1, L1, mt)
                    except Exception as ex:
                        goc_r[gten] = (None, None, None, repr(ex))
                if goc_r[gten][0] is None:
                    KQ["ket_qua"].append({"te_bao": te_bao, "loai": x["loai"], "kieu": x["kieu"], "luat_dich": x["dich"], "bo_phan": x["bp"],
                                          "loi_cay": "goc rieng %s loi: %s" % (gten, goc_r[gten][3])}); continue
                base, M_src, L_src = goc_r[gten][0], goc_r[gten][1], goc_r[gten][2]
            M, L = copy.deepcopy(M_src), copy.deepcopy(L_src)
            rec = {"te_bao": te_bao, "loai": x["loai"], "kieu": x["kieu"], "luat_dich": x["dich"], "bo_phan": x["bp"]}
            if x.get("goc_rieng"):
                rec["goc_rieng"] = {"ten": x["goc_rieng"][0], "mo_ta": goc_r[x["goc_rieng"][0]][3]}
            if x.get("ten_cu"): rec["ten_cu"] = x["ten_cu"]
            Bmap = {q["id"]: q for q in base["luat"]}
            khl, con, bo = V3.loc_hop_le(x["dich"], x["kieu"], Bmap, lambda rid: bp_truot_o_goc(base, rid, x["bp"]), x.get("toan_te_bao", False))
            if khl:
                rec["khong_hop_le"] = khl
                KQ["ket_qua"].append(rec); print(te_bao, x["loai"], x["bp"], "KHONG_HOP_LE", khl, flush=True); continue
            if bo:
                rec["muc_dich_bo_vi_truot_o_goc"] = bo; x = dict(x, dich=con); rec["luat_dich"] = con
            try:
                rec["mo_ta"] = x["fn"](M, L)
            except Exception as ex:
                rec["loi_cay"] = repr(ex); KQ["ket_qua"].append(rec); print(te_bao, x["loai"], "LOI CAY", ex, flush=True); continue
            d = os.path.join(OUT, "%s__%02d_%s" % (te_bao, k, x["loai"])); ghi(M, L, d)
            r = _K3(d)
            mv = truot_moi(base, r, True); mnv = truot_moi(base, r, False)
            rec["luat_truot_moi"] = sorted(mnv.keys())
            rec["du_phong"] = [n for n, e in r["loi_luoi"].items() if e.get("CHAN")]
            if x["kieu"] == "dot_bien":
                thay = False; do = []
                for rid in x["dich"]:
                    if rid == "chong_lan_cung_loai":
                        moi_ = cung_loai_moi(base, r, x["bp"]); thay |= bool(moi_); do += [{"luat": rid, "bo_phan": c["bo_phan"], "gia_tri": c["giao_um3"]} for c in moi_[:2]]
                    elif rid in mnv:
                        moi_ = [c for c in mnv[rid]["moi"] if nhac(c, x["bp"])]
                        thay |= bool(moi_); do += [{"luat": rid, "bo_phan": c.get("bo_phan"), "gia_tri": c.get("gia_tri")} for c in moi_[:2]]
                rec["phat_hien"] = thay; rec["do_duoc"] = do
                rec["phat_hien_theo_phan_quyet"] = V3.theo_pq(x["dich"], {q: v["moi"] for q, v in mnv.items()}, lambda c: nhac(c, x["bp"]))
            elif x["kieu"] == "st_dat_ts_truot":
                pq = V3.theo_pq(x["dich"], {q: v["moi"] for q, v in mnv.items()}, lambda c: nhac(c, x["bp"]))
                rec["phat_hien_theo_phan_quyet"] = pq; rec["dung_ky_vong"] = pq["tham_so"] and not pq["su_that"]
            elif x["kieu"] == "sat_nguong_dat":
                gia = [c for rid in x["dich"] if rid in mnv for c in mnv[rid]["moi"] if nhac(c, x["bp"])]
                rec["bao_dong_gia"] = bool(gia); rec["chi_tiet"] = gia[:2]
                rr = [c for q in r["luat"] if q["id"] in x["dich"] for c in q["ca_dat_mau"] + q["ca_truot"] if nhac(c, x["bp"])]
                rec["gia_tri_luat_dich"] = [c.get("gia_tri") for c in rr][:2]
            else:
                gia = sorted(q for q in mnv if q not in x["ngoai_le"])
                rec["bao_dong_gia"] = bool(gia); rec["luat_bao_gia"] = gia
                rec["chi_tiet"] = {q: mnv[q]["moi"][:2] for q in gia}
                rec["khoi_dac_nhan"] = {k_: v for k_, v in (r["khoi_dac_nhan"] or {}).items() if k_ != "nhat_ky_tim_r"}
                rec["M.d"] = [q["trang_thai"] for q in r["luat"] if q["id"] == "M.d@su_that"][0]
            rec["thoi_gian_s"] = r["thoi_gian_s"]
            rec["M.c0_kkd_nhan"] = (r.get("M.c0_khong_kiem_duoc_nhan") or {}).get("so_cap")
            KQ["ket_qua"].append(rec)
            if os.environ.get("XOA_MESH_DOT_BIEN", "1") == "1":
                _don_thu_muc(d)
            tag = ("PHAT_HIEN" if rec.get("phat_hien") else "BO_SOT") if x["kieu"] == "dot_bien" else (("DUNG" if rec["dung_ky_vong"] else "SAI") if x["kieu"] == "st_dat_ts_truot" else ("BAO_GIA" if rec["bao_dong_gia"] else "SACH"))
            print("%-13s %-26s %-14s %-9s moi=%s %s" % (te_bao, x["loai"], x["bp"], tag, rec["luat_truot_moi"], ("du_phong=%s" % rec["du_phong"]) if rec["du_phong"] else ""), flush=True)
    if os.environ.get("XOA_MESH_DOT_BIEN", "1") == "1":
        for g_ in os.listdir(OUT):
            if "__goc" in g_: _don_thu_muc(os.path.join(OUT, g_))
    db = [x for x in KQ["ket_qua"] if x["kieu"] == "dot_bien" and "phat_hien" in x]
    dc = [x for x in KQ["ket_qua"] if x["kieu"] != "dot_bien" and "bao_dong_gia" in x]
    ti = collections.defaultdict(lambda: [0, 0])
    for x in db: ti[x["loai"]][1] += 1; ti[x["loai"]][0] += int(x["phat_hien"])
    phu = collections.defaultdict(lambda: [0, 0])
    for x in db:
        for rid in x["luat_dich"]: phu[rid][1] += 1; phu[rid][0] += int(x["phat_hien"])
    luat_all = sorted({q["id"] for q in base_toan["luat"]})
    KQ["khong_hop_le"] = [(x["te_bao"], x["loai"], x["bo_phan"], x["khong_hop_le"]) for x in KQ["ket_qua"] if "khong_hop_le" in x]
    KQ["ti_le_theo_loai"] = {k: "%d/%d" % tuple(v) for k, v in ti.items()}
    KQ["do_phu_theo_luat"] = {rid: ("%d/%d" % tuple(phu[rid]) if rid in phu else "khong co dot bien") for rid in luat_all}
    KQ["doi_chung"] = [{"te_bao": x["te_bao"], "loai": x["loai"], "kieu": x["kieu"], "bao_dong_gia": x["bao_dong_gia"]} for x in dc]
    KQ["loi_cay"] = [x for x in KQ["ket_qua"] if "loi_cay" in x]
    KQ["du_phong_da_chay"] = [(x["te_bao"], x["loai"], x["du_phong"]) for x in KQ["ket_qua"] if x.get("du_phong")]
    KQ["tong"] = {"so_dot_bien": len(db), "phat_hien": sum(x["phat_hien"] for x in db),
                  "so_khong_hop_le_bo_qua": sum(1 for x in KQ["ket_qua"] if "khong_hop_le" in x),
                  "so_doi_chung": len(dc), "doi_chung_bao_gia": sum(x["bao_dong_gia"] for x in dc),
                  "so_luat_co_dot_bien": sum(1 for v in KQ["do_phu_theo_luat"].values() if v != "khong co dot bien"),
                  "thoi_gian_s": round(time.time() - t0, 1)}
    KQ["v12"] = V3.tong_ket_v12(KQ["ket_qua"])
    json.dump(KQ, open(os.path.join(HERE, "cay_loi_ket_qua.json"), "w"), ensure_ascii=False, indent=1, default=kiem2._json_default)
    print(json.dumps(KQ["v12"], ensure_ascii=False))
    print(json.dumps(KQ["ti_le_theo_loai"], ensure_ascii=False)); print(json.dumps(KQ["do_phu_theo_luat"], ensure_ascii=False)); print(KQ["tong"])


if __name__ == "__main__":
    chay(tuple(sys.argv[1:]) if len(sys.argv) > 1 else ("D1_td_1001", "D4A2_td_1001"))
