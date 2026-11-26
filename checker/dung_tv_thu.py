# -*- coding: utf-8 -*-
"""
dung_tv_thu.py — dựng TẾ BÀO THỰC VẬT THỬ, đúng theo cách dựng, cho bộ kiểm kiem_thucvat.py.

Các tế bào này do người viết bộ kiểm tự thiết kế từ LỜI LUẬT v1.4 (và bảng số ở Phụ lục 1b):
không dùng mã dựng hay trường máy của bộ dựng. Đây chỉ là dữ liệu kiểm thử, không phải bộ dựng.

Quy ước định dạng (do bộ kiểm định nghĩa, bộ dựng thật nên xuất theo):
  mesh.npz : <ten>__v, <ten>__f  (µm)
  loai.json: {"lop": "Sinh10"|"KHTN6", "loai": {ten: loai}, "cha": {ten_con: ten_cha},
              "bo_phan_du_kien": {ten: loai}}
  loại: vach, mang, khong_bao, nhan, mang_nhan, lo_nhan, nhan_con, luc_lap, mang_trong_lap,
        granum, thylakoid, tam_noi, ti_the, mao, golgi, tui_golgi, peroxisome, luoi_noi_chat,
        ribosome (một bộ phận gồm nhiều chấm)
Dùng:  python dung_tv_thu.py <thu_muc_ra> <TV1|TV2|TV6>
"""
import os, sys, json
import numpy as np
import manifold3d as mf

SEG = 48


def mesh_of(m):
    mm = m.to_mesh64()
    return np.asarray(mm.vert_properties)[:, :3].copy(), np.asarray(mm.tri_verts).astype(np.int64)


def from_mesh(V, F):
    return mf.Manifold(mf.Mesh64(np.ascontiguousarray(V, float), np.ascontiguousarray(F, np.uint64)))


def sphere(c, r, seg=SEG):
    return mf.Manifold.sphere(r, seg).translate(tuple(map(float, c)))


def khung(n):
    n = np.asarray(n, float); n = n / np.linalg.norm(n)
    a = np.array([1.0, 0, 0]) if abs(n[0]) < 0.9 else np.array([0, 1.0, 0])
    u = np.cross(n, a); u /= np.linalg.norm(u); v = np.cross(n, u)
    return u, v, n


def dat(m, o, u, v, n):
    M = np.zeros((3, 4)); M[:, 0] = u; M[:, 1] = v; M[:, 2] = n; M[:, 3] = o
    return m.transform(M)


def tru(h, r, o, n, seg=32):
    u, v, nn = khung(n)
    return dat(mf.Manifold.cylinder(h, r, r, seg, True), o, u, v, nn)


def nang(p, q, r, seg=24):
    return mf.Manifold.batch_hull([sphere(p, r, seg), sphere(q, r, seg)])


A_, B_, C_ = 60.0, 30.0, 30.0; RB = 3.0
TAM_GOC = [(sx * (A_ / 2 - RB), sy * (B_ / 2 - RB), sz * (C_ / 2 - RB)) for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]


def hop(d, seg=64):
    return mf.Manifold.batch_hull([sphere(c, RB - d, seg) for c in TAM_GOC])


def thau_kinh(o, n, D, t, seg=96, lui=0.0):
    a = D / 2; h = t / 2; Rs = (a * a + h * h) / (2 * h)
    n = np.asarray(n, float); n /= np.linalg.norm(n)
    return sphere(o + n * (Rs - h), Rs - lui, seg) ^ sphere(o - n * (Rs - h), Rs - lui, seg)


_ICO = None


def nở(m, r):
    """Nở một khối LỒI thêm ≥ r (bao lồi của đỉnh ⊕ điểm cầu đã phóng để ngoại tiếp)."""
    global _ICO
    if _ICO is None:
        Vs, _ = mesh_of(mf.Manifold.sphere(1.0, 16)); _ICO = Vs / 0.97
    V, _ = mesh_of(m)
    P = (V[:, None, :] + r * _ICO[None]).reshape(-1, 3)
    import trimesh
    P = np.unique(np.round(P, 9), axis=0)
    for lam in (0.0, 1e-7, 1e-6):
        Q = P + (np.random.RandomState(1).randn(*P.shape) * lam if lam else 0)
        h = trimesh.convex.convex_hull(Q)
        h = trimesh.Trimesh(h.vertices, h.faces, process=True)
        h.update_faces(h.area_faces > 1e-12); h.remove_unreferenced_vertices(); h.fix_normals()
        m = from_mesh(h.vertices, h.faces)
        if m.status() == mf.Error.NoError and not m.is_empty():
            return m
    raise RuntimeError("khong dung duoc bao loi no")


def fibo(n):
    i = np.arange(n) + 0.5; phi = np.arccos(1 - 2 * i / n); th = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], 1)


def dung(kieu):
    P = {}; loai = {}; cha = {}
    def them(ten, m, lo, ch=None):
        P[ten] = m; loai[ten] = lo
        if ch: cha[ten] = ch
    big = kieu == "TV2"; lop6 = kieu == "TV6"
    WALL, MEM = 1.2, 0.2
    IN = WALL + MEM
    them("vach", hop(0) - hop(WALL), "vach")
    them("mang", hop(WALL) - hop(IN), "mang")
    zin = -C_ / 2 + IN
    Rn = 3.4 if big else 3.0
    cn = np.array([0.0, 0.0, zin + 0.4 + Rn])
    them("nhan", sphere(cn, Rn, 96), "nhan")
    if not lop6:
        L1, KH = 0.05, 0.07
        r1, r2, r3, r4 = Rn, Rn - L1, Rn - L1 - KH, Rn - 2 * L1 - KH
        dirs = fibo(40 if big else 32)
        lop = sphere(cn, r1, 96) - sphere(cn, r4, 96)
        khoang = sphere(cn, r2, 96) - sphere(cn, r3, 96)
        vanh = [tru(0.6, 0.175, cn + d * (Rn - 0.1), d, 24) for d in dirs]
        lo = [tru(0.6, 0.125, cn + d * (Rn - 0.1), d, 24) for d in dirs]
        khoang = khoang - mf.Manifold.batch_boolean(vanh, mf.OpType.Add)
        env = lop - khoang - mf.Manifold.batch_boolean(lo, mf.OpType.Add)
        them("mang_nhan", env, "mang_nhan", "nhan")
        for k, (d, c_) in enumerate(zip(dirs, lo)):
            them("lo_nhan_%d" % (k + 1), c_ ^ lop, "lo_nhan", "mang_nhan")
        if big:
            them("nhan_con_1", sphere(cn + [0.9, 0, 0.3], 0.8, 48), "nhan_con", "nhan")
            them("nhan_con_2", sphere(cn + [-1.0, 0.2, -0.3], 0.6, 48), "nhan_con", "nhan")
        else:
            them("nhan_con_1", sphere(cn, 0.75, 48), "nhan_con", "nhan")
    D, T = (6.8, 2.4) if big else (6.0, 2.2)
    depth_c = IN + 0.4 + T / 2
    vt = []
    xs = [-22.5, -13.5, -4.5, 4.5, 13.5, 22.5]
    for x in xs:
        for y in (-6.5, 6.5):
            vt.append((np.array([x, y, C_ / 2 - depth_c]), np.array([0, 0, 1.0]), "tren"))
    for x in ([-22.5, -13.5, 13.5, 22.5] if (big or lop6) else [-22.5, 22.5]):
        for y in (-6.5, 6.5):
            vt.append((np.array([x, y, -C_ / 2 + depth_c]), np.array([0, 0, -1.0]), "duoi"))
    for s in (1, -1):
        for x in xs:
            vt.append((np.array([x, s * (B_ / 2 - depth_c), 5.0]), np.array([0, s, 0.0]), "ben"))
    for s in (1, -1):
        vt.append((np.array([s * (A_ / 2 - depth_c), 0.0, 0.0]), np.array([s, 0, 0.0]), "dau"))
    if not (big or lop6):
        vt = vt[:30]
    ll = []
    ngr = 12 if big else 10; nth = 6 if big else 4
    for i, (o, n, mat) in enumerate(vt):
        ten = "luc_lap_%d" % (i + 1)
        them(ten, thau_kinh(o, n, D, T), "luc_lap"); ll.append((ten, o, n, mat))
        if lop6:
            continue
        them("mang_trong_lap_%d" % (i + 1), thau_kinh(o, n, D, T, lui=0.05), "mang_trong_lap", ten)
        u, v, nn = khung(n)
        gx = [-1.875, -1.125, -0.375, 0.375, 1.125, 1.875] if big else [-1.5, -0.75, 0.0, 0.75, 1.5]
        Dg = 0.5 if big else 0.45; tt, gap = 0.05, 0.015
        hg = nth * tt + (nth - 1) * gap
        tam = []
        k = 0
        for yy in (-0.5, 0.5):
            for xx in gx:
                k += 1
                og = o + u * xx + v * yy
                gname = "granum_%d_%d" % (i + 1, k)
                them(gname, tru(hg, Dg / 2, og, nn, 32), "granum", ten)
                tam.append((gname, og))
                for j in range(nth):
                    z0 = -hg / 2 + tt / 2 + j * (tt + gap)
                    them("thylakoid_%d_%d_%d" % (i + 1, k, j + 1), tru(tt, Dg / 2, og + nn * z0, nn, 32), "thylakoid", gname)
        g_all = mf.Manifold.batch_boolean([P[g] for g, _ in tam], mf.OpType.Add)
        cap = []
        m_ = len(gx)
        for row in (0, 1):
            for c in range(m_ - 1):
                cap.append((tam[row * m_ + c], tam[row * m_ + c + 1]))
        mid = m_ // 2
        cap.append((tam[mid], tam[m_ + mid]))
        for j, ((ga, pa), (gb, pb)) in enumerate(cap):
            d = pb - pa; L = np.linalg.norm(d); e = d / L
            w = np.cross(nn, e)
            box = mf.Manifold.cube((L, 0.15, 0.05), True)
            M = np.zeros((3, 4)); M[:, 0] = e; M[:, 1] = w; M[:, 2] = nn; M[:, 3] = (pa + pb) / 2
            them("tam_noi_%d_%d" % (i + 1, j + 1), box.transform(M) - g_all, "tam_noi", ten)
    tis = []
    if not lop6:
        for s in (1, -1):
            for x in (-18, -9, 0, 9, 18):
                tis.append((np.array([x, s * (B_ / 2 - IN - 0.4 - 0.35), -1.0]), np.array([1.0, 0, 0])))
        if big:
            for x in (-18, -9, 0, 9, 18):
                tis.append((np.array([x, 0.0, C_ / 2 - IN - 0.4 - 0.35]), np.array([1.0, 0, 0])))
    for i, (o, ax) in enumerate(tis):
        ten = "ti_the_%d" % (i + 1)
        them(ten, nang(o - ax * 0.4, o + ax * 0.4, 0.35, 32), "ti_the")
        u, v, nn = khung(ax)
        for j, s in enumerate([-0.45, -0.15, 0.15, 0.45]):
            them("mao_%d_%d" % (i + 1, j + 1), tru(0.04, 0.25, o + ax * s, ax, 24), "mao", ten)
    gol = []
    if not lop6:
        ga = [(1, -20), (1, -5), (1, 10), (-1, -12), (-1, 3)]
        if big:
            ga = [(1, -24), (1, -12), (1, 0), (1, 12), (1, 24), (-1, -18), (-1, -6), (-1, 6)]
        nt, Dt, tt, gp = (6, 0.9, 0.05, 0.035) if big else (5, 0.8, 0.04, 0.03)
        hs = nt * tt + (nt - 1) * gp
        for i, (s, x) in enumerate(ga):
            nrm = np.array([0, s, 0.0])
            o = np.array([x, s * (B_ / 2 - IN - 0.4 - hs / 2), -5.0])
            ten = "golgi_%d" % (i + 1)
            them(ten, tru(hs, Dt / 2, o, nrm, 32), "golgi")
            for j in range(nt):
                z0 = -hs / 2 + tt / 2 + j * (tt + gp)
                them("tui_golgi_%d_%d" % (i + 1, j + 1), tru(tt, Dt / 2, o + nrm * z0, nrm, 32), "tui_golgi", ten)
    import kiem2
    per = []
    if not lop6:
        chon = [x for x in ll if x[3] in ("tren", "duoi")]
        ben = [x for x in ll if x[3] == "ben"]
        chon += (ben[:5] + ben[6:11]) if big else [ben[0], ben[2], ben[9], ben[11]]
        for (ten, o, n, mat) in chon:
            if mat in ("tren", "duoi"):
                huong = np.array([0, -np.sign(o[1]), 0.0])
            else:
                huong = np.array([0, 0, -1.0])
            dep = n * (T / 2 - 0.4 - 0.35)
            Vl, Fl = mesh_of(P[ten]); Ll = kiem2.Luoi(Vl, Fl)
            lo_, hi_ = D / 2, D / 2 + 1.5
            for _ in range(30):
                rho = (lo_ + hi_) / 2
                ps = sphere(o + dep + huong * rho, 0.35, 32); Vp, Fp = mesh_of(ps)
                g = kiem2.kc_mat_mat(kiem2.Luoi(Vp, Fp), Ll, 5.0)[0]
                if g < 0.5: lo_ = rho
                else: hi_ = rho
            per.append(sphere(o + dep + huong * hi_, 0.35, 32))
        for i, m in enumerate(per):
            them("peroxisome_%d" % (i + 1), m, "peroxisome")
    ong = []
    if not lop6:
        RE = Rn + 0.7
        vanh = (sphere(cn, RE, 96) - sphere(cn, Rn, 96)) ^ hop(IN + 0.31)
        zt = zin + 0.45
        x0 = 1.5 if not big else 1.8
        for d in ([1, 0], [-1, 0], [0, 1], [0, -1]):
            d = np.array([d[0], d[1], 0.0])
            ong.append(nang(np.array([0, 0, zt]) + d * x0, np.array([0, 0, zt]) + d * 8.0, 0.075, 16))
        er = vanh + mf.Manifold.batch_boolean(ong, mf.OpType.Add)
        them("luoi_noi_chat", er, "luoi_noi_chat")
        if big:
            ds = fibo(9000); ds = ds[(cn[2] + RE * ds[:, 2] >= zin + 0.45)]
            ok = []
            for dd in ds:
                p = cn + dd * (RE + 0.0175)
                if abs(p[2] - zt) < 0.6 and (abs(p[1]) < 0.6 or abs(p[0]) < 0.6):
                    continue
                ok.append(p)
                if len(ok) >= 3500: break
            base = mf.Manifold.sphere(0.0175, 8)
            Vb, Fb = mesh_of(base)
            V = np.concatenate([Vb + p for p in ok]); F = np.concatenate([Fb + k * len(Vb) for k in range(len(ok))])
            P["ribosome"] = ("raw", V, F); loai["ribosome"] = "ribosome"
    core = hop(IN + 1.01)
    no = [nở(P[t], 0.31) for t, l in loai.items() if l in ("luc_lap", "ti_the", "golgi", "peroxisome")]
    rn = (Rn + 0.7 + 0.31) if not lop6 else Rn + 0.31
    no.append(sphere(cn, rn / 0.995, 96))
    if big:
        Vr = P["ribosome"][1].reshape(-1, len(mesh_of(mf.Manifold.sphere(0.0175, 8))[0]), 3).mean(1)
        no += [sphere(p_, (0.0175 + 0.315) / 0.95, 16) for p_ in Vr]
    no += [nở(o_, 0.31) for o_ in ong]
    kb = core - mf.Manifold.batch_boolean(no, mf.OpType.Add)
    parts = sorted(kb.decompose(), key=lambda m: -m.volume())
    them("khong_bao", parts[0], "khong_bao")
    return P, loai, cha


def ghi(thu_muc, kieu):
    P, loai, cha = dung(kieu)
    os.makedirs(thu_muc, exist_ok=True)
    arr = {}
    for t, m in P.items():
        if isinstance(m, tuple):
            V, F = m[1], m[2]
        else:
            V, F = mesh_of(m.simplify(1e-6))
        arr[t + "__v"] = V.astype(np.float64); arr[t + "__f"] = F.astype(np.int64)
    np.savez_compressed(os.path.join(thu_muc, "mesh.npz"), **arr)
    json.dump({"lop": "KHTN6" if kieu == "TV6" else "Sinh10", "loai": loai, "cha": cha, "bo_phan_du_kien": dict(loai),
               "ghi_chu": "te bao thuc vat THU %s do bo kiem tu dung (dung_tv_thu.py)" % kieu},
              open(os.path.join(thu_muc, "loai.json"), "w"), ensure_ascii=False, indent=1)
    print(thu_muc, len(P), "bo phan")


if __name__ == "__main__":
    ghi(sys.argv[1], sys.argv[2])
