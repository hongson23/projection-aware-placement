# -*- coding: utf-8 -*-
"""
do_khuon.py — đo độ lệch của mesh bộ phận cuối so với KHUÔN dựng lại từ tư thế đã lưu (điều kiện khởi đầu lượt chính).

  D_i = max_{p : |SDF_khuôn(p)| ≤ h} | sd_mesh_i(p) − SDF_khuôn_i(p) |        (cả hai có dấu: âm ở trong)
  G_c = Σ_i max(D_i − h, 0)

Khuôn được dựng lại CHỈ từ: tư thế đã lưu của bộ dựng (thong_tin.json → "tu_the": {bộ phận: {khuon, t, s, R}}; quy ước
x = t + s·R·x_khuôn, đã kiểm trên mesh trước chiếu) và định nghĩa khuôn trong tham_so_quy_uoc_v1_2:
  * nhan  : elip trục dv.khuon_nhan_truc (nửa trục, đơn vị s)          * golgi : elip trục dv.khuon_golgi_truc
  * cau, cau_nhan, cau_per : cầu bán kính 1
  * ong<k>: ống ti thể dv.khuon_ong[k] = (bán kính r, độ cong c), ĐÚNG quy_uoc_ma.ong của khuon_thuc_vat_bo_sung_v1_2:
            trục là đường gấp khúc 33 điểm Q(τ) = (τ(1−r), c(1−τ²), 0), τ ∈ [−1, 1]; φ = khoảng cách tới đường gấp khúc − r;
            nửa dài kể cả chỏm = 1, s = nửa trục dài. (Bản trước 11/11 dùng cung tròn — đã thay; so sánh trong chan_doan_gc.py:
            |ΔD| ≤ 0,001 µm, c = 0 trùng hẳn.)
  * thực vật (thau_kinh*, nang*, tru*, cau_nhan, cau_per) và khuôn AI (ai:*, ai_*): theo công thức / công thức chuẩn hoá
    của khuon_thuc_vat_bo_sung_v1_2 (mục thuc_vat, dong_vat_doi_chieu.khuon_ai_D4A2); khuôn không có định nghĩa → KHONG_DUNG_LAI_DUOC.
Điểm đo: lưới đều bước h/4 trong hộp bao khuôn (+h), giữ điểm |SDF| ≤ h. sd_mesh: khoảng cách chính xác (igl AABB) + dấu
theo số vòng quấn. Thử vệ sinh: trên mesh_vo_bao_truoc mọi D_i ≤ h/2; đột biến dịch tư thế (+h theo x) phải làm D_i > h/2.

Dùng: python do_khuon.py <thu_muc_ket_qua_bo_dung> [--vat dong_vat|thuc_vat] [--mesh mesh_vo_bao.npz] [--ra tep.json] [--ve_sinh]
"""
import os, sys, json, argparse, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import luat_v12 as LV
import kiem2 as K

H_VAT = {"dong_vat": 0.1797, "thuc_vat": 0.05}
BUOC_CHIA = None


def _sdf_elip(P, a):
    """khoảng cách có dấu CHÍNH XÁC tới mặt elip trục a (vector hoá; nhân tử Lagrange bằng chia đôi)."""
    a = np.asarray(a, float); a2 = a ** 2
    f0 = ((P / a) ** 2).sum(1) - 1.0
    ra = f0 > 0
    lo = np.where(ra, 0.0, -a2.min() + 1e-15); hi = np.where(ra, np.maximum(1.0, (np.abs(P) * a).max(1) * 2 + a2.max()), 0.0)
    def g(l):
        return ((a * P / (a2 + l[:, None])) ** 2).sum(1) - 1.0
    for _ in range(200):
        m = 0.5 * (lo + hi); gm = g(m)
        lo = np.where(gm > 0, m, lo); hi = np.where(gm > 0, hi, m)
    l = 0.5 * (lo + hi)
    Q = a2 * P / (a2 + l[:, None])
    d = np.linalg.norm(P - Q, axis=1)
    return np.where(ra, d, -d)


def _sdf_polyline(P, Q):
    """khoảng cách chính xác điểm → đường gấp khúc Q (vector hoá theo đoạn)."""
    best = np.full(len(P), np.inf)
    for i in range(len(Q) - 1):
        a, b = Q[i], Q[i + 1]; ab = b - a; L2 = float(ab @ ab)
        u = np.clip(((P - a) @ ab) / L2, 0.0, 1.0)
        best = np.minimum(best, np.linalg.norm(P - (a + u[:, None] * ab), axis=1))
    return best


KB = json.load(open(os.path.join(HERE, "luat_v12", "khuon_thuc_vat_bo_sung_v1_2_20261111.json")))
GOC_SV = os.path.expanduser("~/Vision_team_data_189/sonhh20/")
_AI = {}


def _ai_luoi(muc):
    """khuôn AI: chuẩn hoá đúng công thức của tệp bổ sung; trả Luoi (khoảng cách có dấu tới mesh)."""
    import hashlib, trimesh
    if muc["ten"] in _AI:
        return _AI[muc["ten"]]
    p = GOC_SV + muc["tep_server"]
    h_ = hashlib.sha256(open(p, "rb").read()).hexdigest()
    if h_ != muc["sha256"]:
        raise ValueError("sha256 khuon AI %s khong khop" % muc["ten"])
    z = np.load(p); V = np.asarray(z["v"], float).copy()
    V -= (V.max(0) + V.min(0)) / 2; V = V[:, np.argsort(-np.ptp(V, 0))]; V *= 2.0 / np.ptp(V, 0).max()
    tm = trimesh.Trimesh(V, np.asarray(z["f"]), process=True)
    if tm.volume < 0: tm.invert()
    L = K.Luoi(tm.vertices, tm.faces); _AI[muc["ten"]] = L
    return L


def dinh_nghia_khuon(ten):
    """→ (φ trong khung khuôn, nửa cạnh hộp bao khuôn, mô tả) theo tham_so v1.2 + khuon_*_bo_sung (chỉ hai tệp này)."""
    for m in KB["thuc_vat"]:
        if m.get("ten") != ten: continue
        hinh = m.get("hinh", ""); tp = m.get("tham_so", {})
        if hinh == "cầu":
            return (lambda P: np.linalg.norm(P, axis=1) - tp["ban_kinh"]), np.ones(3) * tp["ban_kinh"], "cau (bo sung)"
        if hinh.startswith("thau_kinh"):
            q = tp["nua_day_z"]; R_ = (q + 1 / q) / 2; c = (1 / q - q) / 2
            def f(P, q=q, R_=R_, c=c):
                x = np.abs(P[:, 2]); y = np.linalg.norm(P[:, :2], axis=1)
                return np.where((y - 1) * c > x, np.sqrt(x ** 2 + (y - 1) ** 2), np.sqrt((x + c) ** 2 + y ** 2) - R_)
            return f, np.array([1.0, 1.0, q]), "thau kinh nua day %g (bo sung)" % q
        if hinh.startswith("viên nang"):
            L_, r = tp["nua_dai_ke_ca_chom"], tp["ban_kinh"]
            f = lambda P, L_=L_, r=r: np.sqrt(np.maximum(np.abs(P[:, 0]) - (L_ - r), 0) ** 2 + P[:, 1] ** 2 + P[:, 2] ** 2) - r
            return f, np.array([L_, r, r]), "vien nang r=%g (bo sung)" % r
        if hinh.startswith("trụ bo tròn"):
            rb, hc, e = tp["ban_kinh"], tp["nua_cao"], tp["ban_kinh_bo"]
            def f(P, rb=rb, hc=hc, e=e):
                d1 = np.linalg.norm(P[:, :2], axis=1) - (rb - e); d2 = np.abs(P[:, 2]) - (hc - e)
                return np.minimum(np.maximum(d1, d2), 0) + np.sqrt(np.maximum(d1, 0) ** 2 + np.maximum(d2, 0) ** 2) - e
            return f, np.array([rb, rb, hc]), "tru bo tron (bo sung)"
        if hinh.startswith("mesh AI"):
            L = _ai_luoi(m)
            def f(P, L=L):
                d = np.sqrt(L.kc2_diem(P)); return np.where(L.quan(P) > 0.5, -d, d)
            return f, np.abs(L.V).max(0), "khuon AI %s" % m["tep_server"]
    if ten in ("cau",):
        return (lambda P: np.linalg.norm(P, axis=1) - 1.0), np.ones(3), "cau ban kinh 1"
    if ten == "nhan":
        a = LV.P("dv.khuon_nhan_truc"); return (lambda P: _sdf_elip(P, a)), np.asarray(a, float), "elip %s (dv.khuon_nhan_truc)" % (a,)
    if ten == "golgi":
        a = LV.P("dv.khuon_golgi_truc"); return (lambda P: _sdf_elip(P, a)), np.asarray(a, float), "elip %s (dv.khuon_golgi_truc)" % (a,)
    if ten.startswith("ong") and ten[3:].isdigit():
        ds = LV.P("dv.khuon_ong"); i = int(ten[3:])
        if i < len(ds):
            r, c = ds[i]; tau = np.linspace(-1, 1, 33)
            Q = np.stack([tau * (1 - r), c * (1 - tau ** 2), 0 * tau], 1)
            return (lambda P, Q=Q, r=r: _sdf_polyline(P, Q) - r), np.array([1.0, r + c, r]), "ong r=%g c=%g (dv.khuon_ong[%d] + quy_uoc_ma)" % (r, c, i)
    for m in KB["dong_vat_doi_chieu"]["khuon_ai_D4A2"]:
        if m["ten"] == ten:
            L = _ai_luoi(m)
            def f(P, L=L):
                d = np.sqrt(L.kc2_diem(P)); return np.where(L.quan(P) > 0.5, -d, d)
            return f, np.abs(L.V).max(0), "khuon AI %s" % m["tep_server"]
    return None


def do_mot(V, F, tt, h, dich=None):
    kh = dinh_nghia_khuon(tt["khuon"])
    if kh is None:
        return {"khuon": tt["khuon"], "trang_thai": "KHONG_DUNG_LAI_DUOC", "ly_do": "khuon khong co dinh nghia trong tep tham so/luat"}
    sdf_t, nua, mo_ta = kh
    R = np.asarray(tt["R"], float); t = np.asarray(tt["t"], float) + (0 if dich is None else np.asarray(dich)); s = float(tt["s"])
    b = h / (4.0 if BUOC_CHIA is None else BUOC_CHIA); m = nua * s + 1.5 * h
    c = max(8 * b, float(m.max()) / 40)
    xs = [np.arange(-m[i], m[i] + 1e-12, c) for i in range(3)]
    G0 = np.stack(np.meshgrid(*xs, indexing="ij"), -1).reshape(-1, 3)
    st0 = s * sdf_t(G0 / s)
    G0 = G0[np.abs(st0) <= h + c * np.sqrt(3) / 2]
    o = np.arange(-c / 2, c / 2, b) + b / 2
    O = np.stack(np.meshgrid(o, o, o, indexing="ij"), -1).reshape(-1, 3)
    Gs = []; Ss = []
    for k0 in range(0, len(G0), max(1, 2_000_000 // len(O))):
        G = (G0[k0:k0 + max(1, 2_000_000 // len(O)), None, :] + O[None]).reshape(-1, 3)
        st = s * sdf_t(G / s); k = np.abs(st) <= h
        Gs.append(G[k]); Ss.append(st[k])
    G = np.concatenate(Gs) if Gs else np.zeros((0, 3)); st = np.concatenate(Ss) if Ss else np.zeros(0)
    Pw = t + G @ R.T
    L = K.Luoi(V, F)
    d = np.sqrt(L.kc2_diem(Pw)); w = L.quan(Pw)
    sd = np.where(w > 0.5, -d, d)
    e = np.abs(sd - st); i = int(e.argmax()) if len(e) else -1
    return {"khuon": tt["khuon"], "mo_ta": mo_ta, "trang_thai": "DO_DUOC", "D_um": round(float(e.max()), 5) if len(e) else None,
            "so_diem": int(len(e)), "trung_vi_lech_um": round(float(np.median(e)), 5) if len(e) else None,
            "diem_lech_max": Pw[i].round(4).tolist() if i >= 0 else None}


def do_te_bao(thu_muc, vat, tep_mesh="mesh_vo_bao.npz", dich_bo_phan=None, h=None):
    global BUOC_CHIA
    h = H_VAT[vat] if h is None else h
    BUOC_CHIA = 4.0 if vat == "dong_vat" else 2.0
    tu_the = json.load(open(os.path.join(thu_muc, "thong_tin.json")))["tu_the"]
    z = np.load(os.path.join(thu_muc, tep_mesh))
    out = {}
    for ten, tt in tu_the.items():
        if ten + "__v" not in z.files:
            out[ten] = {"khuon": tt.get("khuon"), "trang_thai": "KHONG_CO_MESH"}; continue
        dich = None
        if dich_bo_phan and ten in dich_bo_phan:
            dich = dich_bo_phan[ten]
        out[ten] = do_mot(z[ten + "__v"].astype(np.float64), z[ten + "__f"].astype(np.int64), tt, h, dich)
    Ds = {k: v["D_um"] for k, v in out.items() if v.get("trang_thai") == "DO_DUOC"}
    Gc = float(sum(max(D - h, 0.0) for D in Ds.values()))
    return {"thu_muc": thu_muc, "mesh": tep_mesh, "vat": vat, "h_um": h, "buoc_luoi_mau_um": h / BUOC_CHIA, "G_c_um": round(Gc, 5),
            "so_bo_phan_do": len(Ds), "so_khong_dung_lai_duoc": sum(1 for v in out.values() if v.get("trang_thai") == "KHONG_DUNG_LAI_DUOC"),
            "D_max_um": round(max(Ds.values()), 5) if Ds else None, "bo_phan": out}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("thu_muc"); ap.add_argument("--vat", default="dong_vat"); ap.add_argument("--mesh", default="mesh_vo_bao.npz")
    ap.add_argument("--ra", default=None); ap.add_argument("--ve_sinh", action="store_true")
    a = ap.parse_args(); t0 = time.time()
    kq = {"cuoi": do_te_bao(a.thu_muc, a.vat, a.mesh)}
    if a.ve_sinh:
        h = H_VAT[a.vat]
        tr = do_te_bao(a.thu_muc, a.vat, "mesh_vo_bao_truoc.npz")
        tr["ve_sinh_moi_D_le_h_2"] = all(v["D_um"] <= h / 2 for v in tr["bo_phan"].values() if v.get("trang_thai") == "DO_DUOC")
        kq["truoc_chieu"] = tr
        db = do_te_bao(a.thu_muc, a.vat, "mesh_vo_bao_truoc.npz", dich_bo_phan={k: [h, 0, 0] for k in tr["bo_phan"]})
        db["bat_duoc_moi_bo_phan"] = all(v["D_um"] > h / 2 for v in db["bo_phan"].values() if v.get("trang_thai") == "DO_DUOC")
        kq["dot_bien_dich_tu_the_h"] = {k: db[k] for k in ("G_c_um", "D_max_um", "bat_duoc_moi_bo_phan")} | {"D_min_um": min(v["D_um"] for v in db["bo_phan"].values() if v.get("trang_thai") == "DO_DUOC")}
    kq["thoi_gian_s"] = round(time.time() - t0, 1)
    ra = a.ra or os.path.join(HERE, "do_khuon", os.path.basename(os.path.normpath(a.thu_muc)) + ".json")
    os.makedirs(os.path.dirname(ra), exist_ok=True)
    json.dump(kq, open(ra, "w"), ensure_ascii=False, indent=1)
    c = kq["cuoi"]
    print("%s G_c=%.4f D_max=%s do=%d khong_dung_lai=%d" % (os.path.basename(a.thu_muc), c["G_c_um"], c["D_max_um"], c["so_bo_phan_do"], c["so_khong_dung_lai_duoc"]),
          ("| ve_sinh truoc: %s, D_max %s | dot bien bat %s (D_min %s)" % (kq["truoc_chieu"]["ve_sinh_moi_D_le_h_2"], kq["truoc_chieu"]["D_max_um"],
           kq["dot_bien_dich_tu_the_h"]["bat_duoc_moi_bo_phan"], kq["dot_bien_dich_tu_the_h"]["D_min_um"]) if a.ve_sinh else ""))
