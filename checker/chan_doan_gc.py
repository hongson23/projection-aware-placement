# -*- coding: utf-8 -*-
"""
chan_doan_gc.py — chẩn đoán thử vệ sinh G_c trên mesh TRƯỚC chiếu (bo_kiem_v3, độc lập bộ dựng).

Tách ba nguyên nhân khi D_i > h/2 trên mesh_vo_bao_truoc:
  (1) lệch ĐỊNH NGHĨA khuôn   — so D theo định nghĩa tài liệu (đường gấp khúc Q(τ), quy_uoc_ma.ong) với cách hiểu cũ
                                 (cung tròn qua ba điểm (±(1−r),0,0), (0,c,0));
  (2) độ PHÂN GIẢI lưới        — mô phỏng: lấy mẫu SDF khuôn CHÍNH XÁC tại đúng lưới của mesh bộ dựng (bước a và gốc o
                                 SUY TỪ ĐỈNH mesh: đỉnh marching cubes có 2/3 toạ độ nằm trên nút lưới), chạy marching
                                 cubes (skimage) → mesh_mo_phong; đo D_mp bằng cùng thước; H = Hausdorff(mesh bộ dựng,
                                 mesh_mo_phong). Tư thế đúng + chỉ do lưới ⇒ D_bd ≈ D_mp và H ≪ a.
  (3) lệch TƯ THẾ thật         — H lớn (≈ độ dịch) dù định nghĩa và lưới giống nhau.
Đột biến: dịch tư thế δ ∈ {h/4, h/2, h} theo +x thế giới → D_bd(δ), H(δ) = Hausdorff(mesh bộ dựng, mô phỏng ở tư thế dịch).

Chỉ đọc: thong_tin.json['tu_the'], mesh_vo_bao_truoc.npz; định nghĩa khuôn qua do_khuon.dinh_nghia_khuon (tham_so v1.2 +
khuon_thuc_vat_bo_sung v1.2). Không GPU (KIEM_TANG_TOC=0, CUDA_VISIBLE_DEVICES=).
Dùng: python chan_doan_gc.py <thu_muc> [--vat dong_vat|thuc_vat] [--ra tep.json] [--mesh mesh_vo_bao_truoc.npz]
"""
import os, sys, json, argparse, time
os.environ.setdefault("KIEM_TANG_TOC", "0")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import do_khuon as DK
import kiem2 as K


def suy_luoi(V, a, tinh_chinh=True):
    """gốc lưới o (3,) và bước a; a ban đầu = h, rồi TINH CHỈNH bằng hồi quy x = o + a·k trên các toạ độ nằm trên nút
    (rà soát 12/11: bước thật 23/128 = 0,1796875 ≠ h = 0,1797; không tinh chỉnh thì gốc lệch ≤ 5e-4 µm ở vật lớn).
    Trả (o, tỉ lệ đỉnh trên nút theo từng trục, a, tỉ lệ đỉnh có ≥ 2 toạ độ trên nút — phải ≈ 1 với mesh marching cubes)."""
    def _goc(a):
        o = np.zeros(3)
        for ax in range(3):
            r = np.mod(V[:, ax], a)
            hst, ed = np.histogram(r, bins=2000, range=(0, a))
            m = 0.5 * (ed[hst.argmax()] + ed[hst.argmax() + 1])
            rr = np.mod(r - m + a / 2, a) - a / 2
            k = np.abs(rr) < a * 1e-3
            o[ax] = np.mod(m + (np.median(rr[k]) if k.any() else 0.0), a)
        return o
    o = _goc(a)
    if tinh_chinh:
        for tol in (2e-4, 1e-4, 2e-5):
            ks, xs = [], []
            for ax in range(3):
                kf = (V[:, ax] - o[ax]) / a; k = np.round(kf); m = np.abs(kf - k) * a < tol
                ks.append(k[m]); xs.append(V[m, ax] - 0 * ax)
            k_all = np.concatenate(ks); x_all = np.concatenate(xs)
            if len(k_all) < 10 or np.ptp(k_all) < 2: break
            A_ = np.zeros((len(k_all), 4)); A_[:, 0] = k_all; c = 0
            for ax in range(3):
                A_[c:c + len(ks[ax]), 1 + ax] = 1; c += len(ks[ax])
            sol = np.linalg.lstsq(A_, x_all, rcond=None)[0]; a = float(sol[0]); o = np.mod(sol[1:], a)
    on = np.stack([np.abs(np.mod(V[:, ax] - o[ax] + a / 2, a) - a / 2) < max(a * 1e-4, 2e-6) for ax in range(3)], 1)
    return o, on.mean(0), a, float((on.sum(1) >= 2).mean())


def mo_phong(kh, tt, a, o, dich=None):
    """marching cubes của SDF khuôn CHÍNH XÁC lấy mẫu tại lưới (a, o) — đúng phép rời rạc hoá giả định của mesh trước chiếu."""
    from skimage.measure import marching_cubes
    import trimesh
    sdf_t, nua, _ = kh
    R = np.asarray(tt["R"], float); s = float(tt["s"]); t = np.asarray(tt["t"], float) + (0 if dich is None else np.asarray(dich, float))
    he = np.abs(R) @ (np.asarray(nua, float) * s) + 2.5 * a
    k0 = np.floor((t - he - o) / a).astype(int); k1 = np.ceil((t + he - o) / a).astype(int)
    ax = [o[i] + a * np.arange(k0[i], k1[i] + 1) for i in range(3)]
    G = np.stack(np.meshgrid(*ax, indexing="ij"), -1).reshape(-1, 3)
    psi = s * sdf_t(((G - t) @ R) / s)
    vol = psi.reshape(len(ax[0]), len(ax[1]), len(ax[2]))
    if vol.min() >= 0:
        return None, None
    v, f, _, _ = marching_cubes(vol, 0.0, spacing=(a, a, a))
    v = v + np.array([ax[0][0], ax[1][0], ax[2][0]])
    m = trimesh.Trimesh(v, f, process=True)
    if m.volume < 0:
        m.invert()
    return np.asarray(m.vertices, float), np.asarray(m.faces, np.int64)


def _mau(V, F):
    return np.vstack([V, V[F].mean(1)])


def hausdorff(V1, F1, V2, F2):
    L1 = K.Luoi(V1, F1); L2 = K.Luoi(V2, F2)
    d12 = np.sqrt(L2.kc2_diem(_mau(V1, F1))); d21 = np.sqrt(L1.kc2_diem(_mau(V2, F2)))
    return float(max(d12.max(), d21.max())), float(np.median(np.concatenate([d12, d21])))


def khuon_cung_tron(r, c):
    """cách hiểu CŨ (để so): đường tâm là cung tròn qua (±(1−r),0,0) và (0,c,0)."""
    L = 1 - r
    if c == 0:
        Q = np.stack([np.linspace(-L, L, 257), np.zeros(257), np.zeros(257)], 1)
    else:
        Rr = (L ** 2 + c ** 2) / (2 * c); cy = c - Rr; th = np.arcsin(L / Rr)
        g = np.linspace(-th, th, 257); Q = np.stack([Rr * np.sin(g), cy + Rr * np.cos(g), 0 * g], 1)
    return (lambda P, Q=Q, r=r: DK._sdf_polyline(P, Q) - r), np.array([1.0, r + c, r]), "cung tron r=%g c=%g (cach hieu cu)" % (r, c)


def ban_kinh_dac_trung(kh_ten, tt, V):
    """bán kính nét nhỏ nhất (µm) — để so với bước lưới."""
    s = float(tt["s"])
    if kh_ten.startswith("ong") and kh_ten[3:].isdigit():
        return DK.LV.P("dv.khuon_ong")[int(kh_ten[3:])][0] * s
    if kh_ten in ("cau", "cau_nhan", "cau_per"):
        return s
    if kh_ten == "nhan":
        return min(DK.LV.P("dv.khuon_nhan_truc")) * s
    if kh_ten == "golgi":
        return min(DK.LV.P("dv.khuon_golgi_truc")) * s
    return None


def chan_doan(thu_muc, vat, tep_mesh="mesh_vo_bao_truoc.npz", chi=None):
    h = DK.H_VAT[vat]; DK.BUOC_CHIA = 4.0 if vat == "dong_vat" else 2.0
    a = h
    tu_the = json.load(open(os.path.join(thu_muc, "thong_tin.json")))["tu_the"]
    z = np.load(os.path.join(thu_muc, tep_mesh))
    out = {}
    for ten, tt in tu_the.items():
        if chi and not any(ten.startswith(c) for c in chi):
            continue
        if ten + "__v" not in z.files:
            continue
        kh = DK.dinh_nghia_khuon(tt["khuon"])
        if kh is None:
            out[ten] = {"khuon": tt["khuon"], "trang_thai": "KHONG_DUNG_LAI_DUOC"}; continue
        V = z[ten + "__v"].astype(np.float64); F = z[ten + "__f"].astype(np.int64)
        o, fr, a_tc, f2 = suy_luoi(V, a)
        r_bd = DK.do_mot(V, F, tt, h)
        rec = {"khuon": tt["khuon"], "s": tt["s"], "so_dinh": len(V), "rho_um": ban_kinh_dac_trung(tt["khuon"], tt, V),
               "ti_le_dinh_tren_nut": fr.round(3).tolist(), "a_tinh_chinh": a_tc, "ti_le_dinh_2_toa_do_tren_nut": round(f2, 4), "D_bd": r_bd["D_um"], "trung_vi_bd": r_bd["trung_vi_lech_um"]}
        Ve, Fe = mo_phong(kh, tt, a_tc, o)
        if Ve is not None:
            r_mp = DK.do_mot(Ve, Fe, tt, h)
            rec["D_mp"] = r_mp["D_um"]; rec["so_dinh_mp"] = len(Ve)
            rec["H_bd_mp"], rec["H_trung_vi"] = [round(x, 5) for x in hausdorff(V, F, Ve, Fe)]
        k = tt["khuon"]
        if k.startswith("ong") and k[3:].isdigit():
            r_, c_ = DK.LV.P("dv.khuon_ong")[int(k[3:])]
            kh_cu = khuon_cung_tron(r_, c_)
            _orig = DK.dinh_nghia_khuon
            DK.dinh_nghia_khuon = lambda n, kh_cu=kh_cu: kh_cu
            try:
                rec["D_bd_cung_tron"] = DK.do_mot(V, F, tt, h)["D_um"]
            finally:
                DK.dinh_nghia_khuon = _orig
        db = {}
        for nhan, d in (("h/4", h / 4), ("h/2", h / 2), ("h", h)):
            dich = np.array([d, 0, 0])
            Dd = DK.do_mot(V, F, tt, h, dich)["D_um"]
            Vd, Fd = mo_phong(kh, tt, a_tc, o, dich)
            Hd = hausdorff(V, F, Vd, Fd)[0] if Vd is not None else None
            db[nhan] = {"D_bd": Dd, "H_bd_mp": round(Hd, 5) if Hd is not None else None}
        rec["dot_bien_dich"] = db
        out[ten] = rec
    return {"thu_muc": thu_muc, "mesh": tep_mesh, "vat": vat, "h_um": h, "a_luoi_um": a, "bo_phan": out}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("thu_muc"); ap.add_argument("--vat", default=None); ap.add_argument("--ra", default=None)
    ap.add_argument("--mesh", default="mesh_vo_bao_truoc.npz"); ap.add_argument("--chi", nargs="*", default=None)
    A = ap.parse_args(); t0 = time.time()
    ten_tm = os.path.basename(os.path.normpath(A.thu_muc))
    vat = A.vat or ("thuc_vat" if ten_tm.startswith("TV") else "dong_vat")
    kq = chan_doan(A.thu_muc, vat, A.mesh, A.chi); kq["thoi_gian_s"] = round(time.time() - t0, 1)
    ra = A.ra or os.path.join(HERE, "chan_doan_gc", ten_tm + ".json")
    os.makedirs(os.path.dirname(ra), exist_ok=True)
    json.dump(kq, open(ra, "w"), ensure_ascii=False, indent=1)
    h = kq["h_um"]; bp = [v for v in kq["bo_phan"].values() if "D_bd" in v]
    tr = [v for v in bp if v["D_bd"] > h / 2]
    print("%s n=%d truot_h/2=%d D_bd_max=%.4f H_max=%.4f t=%.0fs" % (ten_tm, len(bp), len(tr), max(v["D_bd"] for v in bp),
          max(v.get("H_bd_mp", 0) for v in bp), kq["thoi_gian_s"]))
