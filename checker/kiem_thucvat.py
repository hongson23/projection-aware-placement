# -*- coding: utf-8 -*-
"""
kiem_thucvat.py — BỘ KIỂM ĐỘC LẬP cho mô hình TẾ BÀO THỰC VẬT (luật v1.4, 59 luật).

Chỉ dựa trên LỜI LUẬT (tài liệu giáo viên v1.4 + các trường id/lop/phat_bieu/nhan/luat_cung/
pham_vi/so_that/thuc_te_cho_giao_vien và nguong_tu_dat của tệp JSON; KHÔNG dùng trường máy).
Dùng lại hình học chính xác của kiem2.py (bool manifold3d, khoảng cách mặt–mặt chính xác).

Gọi:  python kiem2.py <thu_muc> --vat thuc_vat [--lop KHTN6|Sinh10]
Định dạng đầu vào: xem dung_tv_thu.py (loai.json có "lop", "loai", "cha", "bo_phan_du_kien").
"""
import os, json, time, collections
import numpy as np
import trimesh
import igl
from scipy.spatial import cKDTree
import kiem2
from kiem2 import (BoPhan, Luoi, kc_mat_mat, tao_manifold, ket_qua, tong_hop, trong_khoang, EPS_VOL, _json_default)

HERE = os.path.dirname(os.path.abspath(__file__))

TV_TIEP_XUC = 0.01
TV_TOL_HOP = 0.02
TV_TOL_THE_TICH_HOP = 0.10
TV_HOP_BO = (0.97, 0.9995)
TV_TOL_KHE_LAP = 0.01
TV_NG_CAU = 0.9
TV_NG_LOI = 0.98
TV_NEP_SAU = 0.02
TV_THAU_KINH = (0.48, 0.60)
TV_NG_TRON_DIA = 0.9
TV_SO_MAT_RAI = 5
TV_KB_THIEU_TOI_DA = 1.0
TV_KB_DELTA = 0.05
TV_KB_BUOC = 0.2
TV_TOL_ER = 0.05
TV_GOLGI_GAN_NHAN = 2.0
TV_GOLGI_TI_LE_GAN = 0.5
TV_GOLGI_SO_VUNG = 3
TV_SAI_LOAI_DU = 0.3
TV_PHANG = 0.5
TV_GOC_THANG = 10.0
TV_KE_NHAU = 3
TV_LECH_CHONG = 0.05
TV_MD_TIA = 0.5
import luat_v12 as _LV
_P = _LV.P
for _k in ("TV_TOL_HOP", "TV_TOL_THE_TICH_HOP", "TV_TOL_KHE_LAP", "TV_NG_CAU", "TV_NG_LOI", "TV_NEP_SAU", "TV_NG_TRON_DIA",
           "TV_KB_THIEU_TOI_DA", "TV_TOL_ER", "TV_GOLGI_GAN_NHAN", "TV_PHANG", "TV_GOC_THANG", "TV_KE_NHAU", "TV_LECH_CHONG", "TV_MD_TIA"):
    globals()[_k] = _P("kiem.tv." + _k)
TV_HOP_BO = tuple(_P("kiem.tv.TV_HOP_BO")); TV_THAU_KINH = tuple(_P("kiem.tv.TV_THAU_KINH"))
TV_SO_MAT_RAI = int(_P("tv.so_mat_rai"))
TV_TIEP_XUC = float(_P("chung.tiep_xuc_toi_da_um"))
KHE_ROI = float(_P("chung.khe_roi")); KHE_TRONG_LL = float(_P("tv.khe_trong_luc_lap"))
KCL_TV = dict(_P("chung.khe_cung_loai")["tv"])
MIEN_K1 = {frozenset(x) for x in _P("tv.mien_khe_chung")}
KT_VACH = _P("tv.kich_thuoc_vach"); DAY_VACH = tuple(_P("tv.day_vach")); DAY_MANG = tuple(_P("tv.day_mang"))
KHE_MANG_VACH = float(_P("tv.khe_mang_vach")); LOP_CTB = float(_P("tv.lop_chat_te_bao"))
CACH_KB = _P("tv.cach_dung_khong_bao")
KHE_NHAN_MANG = tuple(_P("tv.khe_nhan_mang")); DK_NHAN = tuple(_P("tv.duong_kinh_nhan"))
KHE_LL_MANG = tuple(_P("tv.khe_luc_lap_mang")); KT_LL_VE = _P("tv.kich_thuoc_luc_lap_ve"); SO_LL = tuple(_P("tv.so_luc_lap"))
LECH_LL = float(_P("tv.lech_huong_luc_lap")); TT_LOP6 = tuple(_P("tv.ti_the_lop6")["so_luong"])
GAN_MANG_NHAN = float(_P("tv.gan_mang_hoac_nhan"))
SO_NHAN_CON = tuple(_P("tv.so_nhan_con")); DK_NHAN_CON = tuple(_P("tv.duong_kinh_nhan_con"))
LO_NHAN_VE = _P("tv.lo_nhan_ve"); MANG_NHAN_VE = _P("tv.mang_nhan_ve"); KT_TT_VE = _P("tv.kich_thuoc_ti_the_ve"); SO_TT = tuple(_P("tv.so_ti_the"))
MAO_VE = _P("tv.mao_ve"); VUNG_ER = float(_P("tv.vung_noi_er_mang_nhan")); ER_VE = _P("tv.luoi_noi_chat_ve"); RB_VE = _P("tv.ribosome_ve")
KHE_GOLGI = float(_P("tv.khe_golgi_golgi")); SO_GOLGI = tuple(_P("tv.so_golgi")); GOLGI_VE = _P("tv.golgi_ve")
DK_PE_VE = tuple(_P("tv.duong_kinh_peroxisome_ve")); SO_PE = tuple(_P("tv.so_peroxisome")); GAN_PE_LL = float(_P("tv.gan_peroxisome_luc_lap"))
KHE_MT_LAP = float(_P("tv.khe_mang_trong_lap")); SO_GRANUM = tuple(_P("tv.so_granum")); SO_THYL_VE = tuple(_P("tv.so_thylakoid_ve"))
THYL_VE = _P("tv.thylakoid_ve"); VUNG_TAM = float(_P("tv.vung_noi_tam_granum")); SO_TAM_NOI_MAX = int(_P("tv.so_tam_noi")["toi_da"])
LIEN_KHOI_MIEN_TV = set(_P("tv.lien_khoi_mien"))
CHO_DUYET_TV = {k: v for k, v in globals().items() if k.startswith("TV_")}

LOAI_CHA = {"mang_nhan": "nhan", "lo_nhan": "mang_nhan", "nhan_con": "nhan", "mang_trong_lap": "luc_lap",
            "granum": "luc_lap", "tam_noi": "luc_lap", "thylakoid": "granum", "mao": "ti_the", "tui_golgi": "golgi"}
TREN = ["nhan", "luc_lap", "ti_the", "golgi", "peroxisome", "luoi_noi_chat", "ribosome", "khong_bao"]
NHIEU_KHOI = {"ribosome", "luoi_noi_chat", "golgi"}
LOP10_LOAI = {"mang_nhan", "lo_nhan", "nhan_con", "mang_trong_lap", "granum", "thylakoid", "tam_noi", "mao",
              "golgi", "tui_golgi", "peroxisome", "luoi_noi_chat", "ribosome"}


def _flip_man(c):
    c2 = c.copy(); c2.invert(); m, _ = tao_manifold(c2.vertices, c2.faces); return m, Luoi(c2.vertices, c2.faces)


class Vo:
    """Vỏ (thành, màng, màng nhân): mặt ngoài O (thể tích dương lớn nhất) và mặt trong I (âm lớn nhất)."""
    def __init__(self, bp):
        pos = [c for c in bp.comps if c.volume > 0]; neg = [c for c in bp.comps if c.volume < 0]
        self.O = max(pos, key=lambda c: c.volume) if pos else None
        self.I = min(neg, key=lambda c: c.volume) if neg else None
        self.fill, _ = tao_manifold(self.O.vertices, self.O.faces) if self.O is not None else (None, None)
        self.LO = Luoi(self.O.vertices, self.O.faces) if self.O is not None else None
        if self.I is not None:
            self.trong, self.LI = _flip_man(self.I)
        else:
            self.trong, self.LI = None, None

    def do_day(self):
        if self.LO is None or self.LI is None:
            return None, None
        dmin = min(float(np.sqrt(self.LO.kc2_diem(self.LI.V)).min()), float(np.sqrt(self.LI.kc2_diem(self.LO.V)).min()))
        dmax = max(float(np.sqrt(self.LO.kc2_diem(self.LI.V)).max()), float(np.sqrt(self.LI.kc2_diem(self.LO.V)).max()))
        return dmin, dmax


def giao(A, B):
    if A.man is not None and B.man is not None:
        return float((A.man ^ B.man).volume())
    return kiem2.du_phong_the_tich(A, B, "giao")


def hop_cach(A, B):
    return kiem2.hop_cach(A, B)


def kc(A, B, R):
    """Khoảng cách chính xác mặt–mặt (< R) giữa hai bộ phận RỜI nhau: manifold3d min_gap (BVH, chính xác;
    đã đối chiếu = kc_mat_mat trên 1 114 cặp); dự phòng kc_mat_mat khi thiếu manifold."""
    if A.man is not None and B.man is not None:
        P, L_ = (A.luoi.V, B.luoi) if len(A.luoi.V) <= len(B.luoi.V) else (B.luoi.V, A.luoi)
        U = float(np.sqrt(L_.kc2_diem(P)).min())
        s_ = min(float(R), U * (1 + 1e-9) + 1e-9)
        d = float(A.man.min_gap(B.man, s_))
        if d >= s_:
            d = min(float(R), U)
        return d, d < R
    if max(A.luoi.Lmax, B.luoi.Lmax) > 1.0:
        d = min(float(np.sqrt(B.luoi.kc2_diem(A.luoi.V)).min()), float(np.sqrt(A.luoi.kc2_diem(B.luoi.V)).min()))
        return (d if d < R else R), d < R
    d, cx = kc_mat_mat(A.luoi, B.luoi, R)
    return (d if cx else R), cx


def khe(A, B, R):
    """(khe chính xác nếu < R, giao µm³, chính xác?). Chồng lấn → khe 0."""
    if hop_cach(A, B) >= R:
        return R, 0.0, False
    g = giao(A, B)
    if g > EPS_VOL:
        return 0.0, g, True
    d, cx = kc(A, B, R)
    return d, g, cx


def trong(A, Cman, CL):
    """A nằm trong khối Cman (mặt CL): bool + phép thử đỉnh khi chạm."""
    if A.man is not None and Cman is not None:
        hv = float((A.man - Cman).volume())
    else:
        w = CL.quan(A.luoi.V); hv = float((w < 0.5).mean()) * 1.0
    if hv > EPS_VOL:
        return False, hv
    w = CL.quan(A.luoi.V); dd = np.sqrt(CL.kc2_diem(A.luoi.V))
    if ((w < 0.5) & (dd > 1e-9)).any():
        return False, hv
    gan = (CL.V >= A.lo - 1e-6).all(1) & (CL.V <= A.hi + 1e-6).all(1)
    if gan.any():
        P = CL.V[gan]; wa = A.luoi.quan(P); da = np.sqrt(A.luoi.kc2_diem(P))
        if ((wa > 0.5) & (da > 1e-9)).any():
            return False, hv
    return True, hv


def dien_tich_cham(A, LB, tol=TV_TIEP_XUC):
    """Diện tích mặt A có cả 3 đỉnh cách mặt B ≤ tol (µm²)."""
    d = np.sqrt(LB.kc2_diem(A.luoi.V))
    ok = (d[A.luoi.F] <= tol).all(1)
    V = A.luoi.V; F = A.luoi.F[ok]
    if len(F) == 0:
        return 0.0
    return float(0.5 * np.linalg.norm(np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]]), axis=1).sum())


def do_lom(bp):
    """Độ lõm sâu nhất: khoảng cách lớn nhất từ đỉnh lưới tới mặt bao lồi của nó (lồi ⇒ ≈ 0)."""
    h = bp.tm.convex_hull
    return float(np.sqrt(Luoi(h.vertices, h.faces).kc2_diem(bp.luoi.V)).max())


def truc_pca(V):
    c = V.mean(0); w, U = np.linalg.eigh(np.cov((V - c).T)); return c, U


def cap_ung_vien(As, Bs, R, cung=False):
    """Cặp (a, b) có khoảng cách hộp bao < R (lọc bằng KD-tree tâm hộp)."""
    if not As or not Bs:
        return []
    ca = np.array([(a.lo + a.hi) / 2 for a in As]); ra = np.array([np.linalg.norm(a.hi - a.lo) / 2 for a in As])
    cb = np.array([(b.lo + b.hi) / 2 for b in Bs]); rb = np.array([np.linalg.norm(b.hi - b.lo) / 2 for b in Bs])
    t = cKDTree(cb); out = []
    for i, a in enumerate(As):
        for j in t.query_ball_point(ca[i], r=ra[i] + rb.max() + R):
            b = Bs[j]
            if b is a or (cung and id(b) <= id(a)):
                continue
            if hop_cach(a, b) < R:
                out.append((a, b))
    return out


DAI_LOAI = {
    "nhan": ((5.2, 7), (5.2, 7), (5.2, 7)), "nhan_con": ((1, 2), (1, 2), (1, 2)),
    "luc_lap": ((5, 7), (5, 7), (2, 2.5)), "mang_trong_lap": ((4.8, 7), (4.8, 7), (1.8, 2.5)),
    "granum": ((0.3, 0.6), (0.3, 0.6), (0.14, 0.62)), "thylakoid": ((0.3, 0.6), (0.3, 0.6), (0.04, 0.06)),
    "ti_the": ((1, 2), (0.5, 1), (0.5, 1)), "mao": (None, None, (0.03, 0.05)),
    "golgi": ((0.5, 1), (0.5, 1), (0.18, 0.56)), "tui_golgi": ((0.5, 1), (0.5, 1), (0.03, 0.06)),
    "peroxisome": ((0.5, 1), (0.5, 1), (0.5, 1)), "lo_nhan": ((0.2, 0.3), None, None),
    "tam_noi": (None, None, (0.0, 0.1)),
}
PHU_THUOC_T2 = {"5.1a": ["nhan_con"], "5.1b": ["nhan_con"], "5.2a": ["mang_nhan", "lo_nhan"], "5.2b": ["lo_nhan"],
                "5.2c": ["mang_nhan"], "5.3e": ["mao"], "5.6c": ["tui_golgi"], "10.1a": ["mang_trong_lap"],
                "10.1b": ["mang_trong_lap"], "10.2a": ["granum", "thylakoid"], "10.2b": ["granum"],
                "10.2c": ["thylakoid"], "10.2d": ["thylakoid"], "10.3a": ["tam_noi"], "10.3b": ["tam_noi"],
                "10.3c": ["tam_noi"], "10.5": ["tam_noi"]}


def hop_loai(e, lo_):
    d = DAI_LOAI.get(lo_)
    if d is None: return None
    m = TV_SAI_LOAI_DU
    return all(r is None or (r[0] / (1 + m) <= x <= r[1] * (1 + m)) for x, r in zip(e, d))


def kiem_loai(tb):
    """Mỗi bộ phận phải có hình/cỡ hợp loại khai báo (biên ± TV_SAI_LOAI_DU). Không hợp → 'nghi_sai_loai',
    kèm các loại khác mà cỡ của nó khớp. Loại chưa có dải (vach, mang, khong_bao, ER, ribosome…) chỉ xét
    các loại 'không được có' (lysosome, trung_the, cau_sinh_chat, khung_xuong) bằng cách gợi ý loại khớp."""
    out = []
    for b in tb.bp.values():
        if b.loai in ("vach", "mang", "khong_bao", "luoi_noi_chat", "ribosome", "mang_nhan"):
            continue
        e = b.kich_thuoc()
        ok = hop_loai(e, b.loai)
        khop = [lo_ for lo_ in DAI_LOAI if lo_ != b.loai and hop_loai(e, lo_)]
        if ok is False or (ok is None and khop):
            out.append({"bo_phan": b.ten, "loai_khai": b.loai, "kich_thuoc": [round(x, 3) for x in e], "khop_loai": khop})
    return out


class TeBaoTV:
    def __init__(self, thu_muc):
        L = json.load(open(os.path.join(thu_muc, "loai.json")))
        z = np.load(os.path.join(thu_muc, "mesh.npz"))
        self.lop_khai = L.get("lop")
        self.loai = dict(L.get("loai", {})); self.du_kien = dict(L.get("bo_phan_du_kien", self.loai))
        self.cha = dict(L.get("cha", {}))
        ten = sorted(set(k.rsplit("__", 1)[0] for k in z.files))
        self.bp = {}
        for t in ten:
            lo = self.loai.get(t) or self.du_kien.get(t) or "khong_ro"
            nhieu = lo in ("ribosome",)
            self.bp[t] = BoPhan(t, lo, z[t + "__v"], z[t + "__f"], kiem_long_nhau=not nhieu, kiem_tu_cat=not nhieu)
        self.mat = sorted(set(self.du_kien) - set(self.bp))
        self.tl = collections.defaultdict(list)
        for b in self.bp.values():
            self.tl[b.loai].append(b)
        self.mat_tl = collections.defaultdict(list)
        for t in self.mat:
            self.mat_tl[self.du_kien[t]].append(t)
        self.cha_suy = {}
        for b in self.bp.values():
            if b.loai in LOAI_CHA and b.ten not in self.cha:
                cand = [p for p in self.tl.get(LOAI_CHA[b.loai], []) if np.all(p.lo <= b.hi) and np.all(b.lo <= p.hi)]
                c = b.luoi.V.mean(0)[None]
                cand = [p for p in cand if p.luoi.quan(c)[0] > 0.5]
                if cand:
                    self.cha[b.ten] = min(cand, key=lambda p: p.the_tich).ten; self.cha_suy[b.ten] = self.cha[b.ten]
        self.con = collections.defaultdict(list)
        for c_, p_ in self.cha.items():
            if c_ in self.bp:
                self.con[p_].append(self.bp[c_])

    def mot(self, lo):
        x = self.tl.get(lo, []); return x[0] if x else None

    def con_loai(self, p, lo):
        return [c for c in self.con.get(p.ten, []) if c.loai == lo]

    def to_tien(self, ten):
        out = []
        while ten in self.cha:
            ten = self.cha[ten]; out.append(ten)
        return out


def ca_mat(tb, loai):
    return [{"bo_phan": t, "dat": False, "gia_tri": "MAT (co trong bo_phan_du_kien, khong co trong mesh.npz)", "nguong": "phai ton tai"}
            for t in tb.mat_tl.get(loai, [])]


def mat_nhieu(tb, loais):
    out = []
    for l in loais: out += ca_mat(tb, l)
    return out


def kiem_tv(thu_muc, lop=None, duong_luat=None, im_lang=False):
    t0 = time.time()
    LU = json.load(open(duong_luat or os.path.join(HERE, "luat_thucvat_v14_van_ban.json")))
    luat = {l["id"]: l for l in LU["luat"]}
    tb = TeBaoTV(thu_muc)
    if lop is None:
        lop = tb.lop_khai or ("Sinh10" if any(b.loai in LOP10_LOAI for b in tb.bp.values()) else "KHTN6")
    vach = tb.mot("vach"); mang = tb.mot("mang"); kb = tb.mot("khong_bao"); nhan = tb.mot("nhan")
    VV = Vo(vach) if vach else None; VM = Vo(mang) if mang else None
    IN_man = VM.trong if VM else None; IN_L = VM.LI if VM else None
    V_in = float(IN_man.volume()) if IN_man is not None else float('nan')
    mn = tb.mot("mang_nhan"); VN = Vo(mn) if mn else None
    LL = tb.tl.get("luc_lap", []); TT = tb.tl.get("ti_the", []); GO = tb.tl.get("golgi", []); PE = tb.tl.get("peroxisome", [])
    ER = tb.tl.get("luoi_noi_chat", []); RB = tb.tl.get("ribosome", [])
    tren = [b for b in tb.bp.values() if b.loai in TREN]
    KQ = []
    def L_(i): return luat[i]
    _dbg = os.environ.get("KTV_DEBUG")
    def _t(n):
        if _dbg:
            import resource, sys
            print("[%.1fs %.1fGB] %s" % (time.time() - t0, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, n), file=sys.stderr, flush=True)
    def add(i, st, ch, ca, **kw):
        KQ.append(ket_qua(L_(i), st, ch, ca, **kw))

    _t("bat dau")
    khe_mang = {}
    for b in tren + tb.tl.get("vach", []):
        if mang is None or b is mang: continue
        if b.loai == "vach": continue
        d, g, cx = khe(b, mang, 50.0)
        khe_mang[b.ten] = (d, g)

    _t('2.1')
    ca = []
    for lo_ in ("vach", "mang", "khong_bao", "nhan", "luc_lap"):
        ca.append({"bo_phan": lo_, "dat": bool(tb.tl.get(lo_)), "gia_tri": len(tb.tl.get(lo_, [])), "nguong": ">= 1"})
    if vach and mang and VV.fill is not None and IN_man is not None:
        ok, hv = trong(mang, VV.fill, VV.LO); ca.append({"bo_phan": "mang", "dat": ok, "gia_tri": round(hv, 6), "nguong": "mang trong khoi thanh (thanh o ngoai cung)"})
        if kb and nhan:
            v_ct = float((IN_man - kb.man - nhan.man).volume()) if kb.man is not None and nhan.man is not None else float('nan')
            ca.append({"bo_phan": "chat_te_bao", "dat": v_ct > EPS_VOL, "gia_tri": round(v_ct, 2), "nguong": "the tich (trong mang − khong bao − nhan) > 0"})
        for b in (kb, nhan) + tuple(LL):
            if b is None: continue
            ok, _ = trong(b, IN_man, IN_L)
            if not ok: ca.append({"bo_phan": b.ten, "dat": False, "gia_tri": "ngoai mang", "nguong": "ben trong mang"})
    add("2.1", tong_hop(ca, False), "Có thành, màng, chất tế bào (thể tích trong màng trừ không bào, nhân > 0), nhân, không bào trung tâm, lục lạp; màng nằm trong khối thành; nhân, không bào, lục lạp nằm trong màng.", ca)

    _t('2.2a / 2.2b')
    ca = []; ca2 = []
    if vach and VV.fill is not None:
        e = sorted(trimesh.Trimesh(VV.O.vertices, VV.O.faces).bounding_box_oriented.primitive.extents)[::-1]
        for k, (x, m) in enumerate(zip(e, (KT_VACH["dai"], KT_VACH["rong"], KT_VACH["cao"]))):
            ca.append({"bo_phan": "vach", "dat": abs(x - m) <= TV_TOL_HOP * m, "gia_tri": round(float(x), 3), "nguong": "%g µm ± %g%%" % (m, 100 * TV_TOL_HOP)})
        v = float(VV.fill.volume())
        ca.append({"bo_phan": "vach", "dat": True, "gia_tri": round(v, 0), "nguong": "thong tin: the tich khoi thanh (v1.2 khong co tham so the tich ngoai; 2.2a su that kiem the tich TRONG mang)"})
        r = v / float(np.prod(e)); hull = float(trimesh.Trimesh(VV.O.vertices, VV.O.faces).convex_hull.volume)
        ca2.append({"bo_phan": "vach", "dat": trong_khoang(r, *TV_HOP_BO) and v / hull >= 0.999, "gia_tri": {"V_tren_hop": round(r, 4), "V_tren_bao_loi": round(v / hull, 4)}, "nguong": "hop bo goc: V/hop trong %s, loi" % (TV_HOP_BO,)})
    else:
        ca.append({"bo_phan": "vach", "dat": False, "gia_tri": "khong co", "nguong": ""}); ca2 = list(ca)
    add("2.2a", tong_hop(ca, False), "Kích thước ngoài của thành (hộp bao định hướng) 60 × 30 × 30 µm (± TV_TOL_HOP) và thể tích khối thành khoảng 54 pL (± TV_TOL_THE_TICH_HOP).", ca)
    add("2.2b", tong_hop(ca2, False), "Hộp chữ nhật bo góc: khối thành lồi, thể tích / thể tích hộp bao trong TV_HOP_BO (nhỏ hơn 1 vì bo góc, gần 1 vì là hộp).", ca2)

    _t('2.3a / 2.3b / 2.3c / 2.3d')
    ca = []
    if vach and mang and VV.fill is not None:
        ok, hv = trong(mang, VV.fill, VV.LO)
        g = giao(mang, vach)
        ca.append({"bo_phan": "mang", "dat": ok and g <= EPS_VOL, "gia_tri": {"ngoai_thanh_um3": round(hv, 6), "giao_vat_lieu_thanh_um3": round(g, 6)}, "nguong": "mang ben trong thanh, khong giao vat lieu thanh"})
        ca.append({"bo_phan": "vach", "dat": VV.I is not None and vach.so_khoi == 1, "gia_tri": vach.so_khoi, "nguong": "thanh la mot vo kin bao quanh toan bo"})
    else:
        ca.append({"bo_phan": "vach/mang", "dat": False, "gia_tri": "thieu", "nguong": ""})
    add("2.3a", tong_hop(ca, False), "Màng nằm trong khối bao bởi mặt ngoài thành, không giao vật liệu thành; thành là một vỏ kín (có mặt trong) bao toàn bộ.", ca)
    def do_vo(V_, b, lo_, hi_, rid, ten_):
        ca = mat_nhieu(tb, [b.loai if b else ten_])
        if b is None or V_ is None:
            ca.append({"bo_phan": ten_, "dat": False, "gia_tri": "khong co", "nguong": ""}); return ca
        dmin, dmax = V_.do_day()
        ca.append({"bo_phan": b.ten, "dat": dmin is not None and dmin >= lo_ - 1e-6 and dmax <= hi_ + 1e-6, "gia_tri": {"day_min": round(dmin, 4) if dmin is not None else None, "day_max": round(dmax, 4) if dmax is not None else None}, "nguong": "[%g, %g] µm" % (lo_, hi_)})
        return ca
    ca = do_vo(VV, vach, DAY_VACH[0], DAY_VACH[1], "2.3b", "vach")
    add("2.3b", tong_hop(ca, False), "Bề dày thành = khoảng cách giữa mặt trong và mặt ngoài: nhỏ nhất chính xác (mặt–mặt), lớn nhất theo đỉnh (hai chiều); cả hai trong [1; 1,5] µm.", ca)
    ca = []
    if vach and mang and VM.LO is not None and VV.LI is not None:
        d = np.sqrt(VV.LI.kc2_diem(VM.LO.V)); f = float((d <= KHE_MANG_VACH + TV_TIEP_XUC).mean())
        ca.append({"bo_phan": "mang", "dat": f >= 1.0 and giao(mang, vach) <= EPS_VOL, "gia_tri": {"ti_le_dinh_ap_sat": round(f, 4), "khe_max_um": round(float(d.max()), 4)}, "nguong": "moi dinh mat ngoai mang cach mat trong thanh <= %g µm, khong giao" % TV_TIEP_XUC})
    else:
        ca.append({"bo_phan": "mang", "dat": False, "gia_tri": "thieu", "nguong": ""})
    add("2.3c", tong_hop(ca, False), "Áp sát, không khe: mọi đỉnh mặt ngoài màng cách mặt trong thành ≤ TV_TIEP_XUC, và màng không giao thành.", ca)
    ca = do_vo(VM, mang, DAY_MANG[0], DAY_MANG[1], "2.3d", "mang")
    add("2.3d", tong_hop(ca, False), "Bề dày lớp màng (mặt ngoài–mặt trong) trong [0,15; 0,25] µm.", ca)

    _t('2.4a / 2.4b / 2.4c / 2.5b')
    KB = tb.tl.get("khong_bao", [])
    ca = mat_nhieu(tb, ["khong_bao"]) + [{"bo_phan": "khong_bao", "dat": len(KB) == 1, "gia_tri": len(KB), "nguong": "dung 1"}]
    for b in KB:
        ca.append({"bo_phan": b.ten, "dat": b.so_khoi == 1, "gia_tri": b.so_khoi, "nguong": "mot khoi"})
    add("2.4a", tong_hop(ca, True), "Đúng một không bào, một khối. 'Có màng riêng bao quanh' (màng không bào) không kiểm được nếu mô hình không có bộ phận màng không bào.", ca, ly_do="Không có loại bộ phận 'màng không bào' trong định dạng; không kiểm được lớp màng riêng.")
    ca = mat_nhieu(tb, ["khong_bao"])
    if kb and IN_man is not None:
        r = float(kb.man.volume()) / V_in if kb.man is not None else float('nan')
        _a, _b, _ = _LV.KHOANG("tv", "2.4b", "68", 1)
        ca.append({"bo_phan": kb.ten, "dat": trong_khoang(r, _a, _b), "gia_tri": round(r, 4), "nguong": "[%g; %g] (luat su that 2.4b v1.2)" % (_a, _b), "them": {"V_trong_mang_um3": round(V_in, 1)}})
    else:
        ca.append({"bo_phan": "khong_bao", "dat": False, "gia_tri": "thieu", "nguong": ""})
    add("2.4b", tong_hop(ca, False), "Thể tích không bào / thể tích khối phía trong mặt trong của màng ∈ [0,70; 0,85].", ca)
    ca5 = mat_nhieu(tb, ["khong_bao"])
    if kb and IN_man is not None:
        ok, hv = trong(kb, IN_man, IN_L)
        d = kc(kb, mang, LOP_CTB)
        ca5.append({"bo_phan": kb.ten, "dat": ok and (not d[1] or d[0] >= LOP_CTB - 1e-6), "gia_tri": round(d[0], 5) if d[1] else ">=%g" % LOP_CTB, "nguong": ">= %g µm toi mat trong mang (tv.lop_chat_te_bao), khong bao ben trong mang" % LOP_CTB})
    add("2.5b", tong_hop(ca5, False), "Không bào nằm trong màng và khoảng cách chính xác từ mặt không bào tới mặt trong màng ≥ 1 µm.", ca5)
    ca = mat_nhieu(tb, ["khong_bao"])
    if kb and IN_man is not None:
        ca.append({"bo_phan": kb.ten, "dat": kb.so_khoi == 1, "gia_tri": kb.so_khoi, "nguong": "lien khoi"})
        ca += [dict(c) for c in KQ[-1]["ca_truot"]] if KQ[-1]["ca_truot"] else []
        for b in tren:
            if b is kb: continue
            d, g, cx = khe(kb, b, CACH_KB["cach_bo_phan_um"])
            if d < CACH_KB["cach_bo_phan_um"] - 1e-6:
                ca.append({"bo_phan": [kb.ten, b.ten], "dat": False, "gia_tri": round(d, 4), "giao_um3": round(g, 6), "nguong": ">= %g µm toi moi bo phan (tv.cach_dung_khong_bao)" % CACH_KB["cach_bo_phan_um"]})
        h = TV_KB_BUOC; lo_, hi_ = mang.lo, mang.hi
        xs = [np.arange(lo_[i] + h / 2, hi_[i], h) for i in range(3)]
        G = np.stack(np.meshgrid(*xs, indexing='ij'), -1).reshape(-1, 3)
        G = G[IN_L.quan(G) > 0.5]
        G = G[np.sqrt(IN_L.kc2_diem(G)) > CACH_KB["cach_mang_um"] + TV_KB_DELTA]
        khac = [b for b in tb.bp.values() if b.loai not in ("vach", "mang", "khong_bao")]
        if khac:
            Vs = np.concatenate([b.luoi.V for b in khac]); off = np.cumsum([0] + [len(b.luoi.V) for b in khac])[:-1]
            Fs = np.concatenate([b.luoi.F + o for b, o in zip(khac, off)])
            U = Luoi(Vs, Fs)
            G = G[(np.sqrt(U.kc2_diem(G)) > CACH_KB["cach_bo_phan_um"] + TV_KB_DELTA) & (U.quan(G) < 0.5)]
        thieu = float((kb.luoi.quan(G) <= 0.5).sum()) * h ** 3 if len(G) else 0.0
        ca.append({"bo_phan": kb.ten, "dat": thieu <= TV_KB_THIEU_TOI_DA, "gia_tri": round(thieu, 3), "nguong": "the tich cho trong hop le ngoai khong bao <= %g µm³ (luoi %g µm, bien %g)" % (TV_KB_THIEU_TOI_DA, h, TV_KB_DELTA), "so_diem": int(len(G))})
        r = float(kb.man.volume()) / V_in if kb.man is not None else float('nan')
        _a, _b, _ = _LV.KHOANG("tv", "2.4b", "68", 1)
        ca.append({"bo_phan": kb.ten, "dat": trong_khoang(r, _a, _b), "gia_tri": round(r, 4), "nguong": "[%g; %g] (luat su that 2.4b v1.2)" % (_a, _b)})
    add("2.4c", tong_hop(ca, False), "Không bào là phần còn lại: (i) liền một khối; (ii) cách màng ≥ 1 µm (như 2.5b); (iii) cách mọi bộ phận ≥ 0,3 µm; (iv) TỐI ĐA: thể tích các điểm lưới nằm trong màng, cách màng > 1 + δ, cách mọi bộ phận > 0,3 + δ mà KHÔNG thuộc không bào ≤ TV_KB_THIEU_TOI_DA; (v) tỉ lệ thể tích 70–85 %.", ca)

    _t('2.6a / 2.6b / 2.6e / 2.6f')
    NH = tb.tl.get("nhan", [])
    ca = mat_nhieu(tb, ["nhan"]) + [{"bo_phan": "nhan", "dat": len(NH) == 1, "gia_tri": len(NH), "nguong": "dung 1"}]
    add("2.6a", tong_hop(ca, False), "Đúng một bộ phận nhân.", ca)
    ca = mat_nhieu(tb, ["nhan"])
    for b in NH:
        d, g = khe_mang.get(b.ten, (None, None))
        ok_kb = (kb is None) or giao(b, kb) <= EPS_VOL
        mat_ = None
        if IN_L is not None:
            _, fi, _ = IN_L.aabb().squared_distance(IN_L.V, IN_L.F, b.luoi.V.mean(0)[None])
            nrm = trimesh.Trimesh(IN_L.V, IN_L.F, process=False).face_normals[fi[0]]
            k = int(np.argmax(np.abs(nrm))); mat_ = ("-" if nrm[k] > 0 else "+") + "xyz"[k]
        ca.append({"bo_phan": b.ten, "dat": d is not None and trong_khoang(d, *KHE_NHAN_MANG) and ok_kb, "gia_tri": {"khe_mang_um": round(d, 4) if d is not None else None, "mat_gan_nhat": mat_, "khong_trong_khong_bao": ok_kb}, "nguong": "khe toi mang %s µm (tv.khe_nhan_mang); khong giao khong bao ('thuong o mat day' chi ghi)" % (KHE_NHAN_MANG,)})
    add("2.6b", tong_hop(ca, False), "Nhân sát màng: khe nhân–mặt trong màng ∈ [0,3; 0,5] µm (số Phụ lục 1b); không giao không bào (không bào lõm quanh nhân). 'Thường ở mặt đáy' chỉ ghi mặt gần nhất.", ca)
    ca = mat_nhieu(tb, ["nhan"])
    for b in NH:
        r = b.the_tich / V_in
        ca.append({"bo_phan": b.ten, "dat": trong_khoang(r, 0.0016, 0.023), "gia_tri": round(r, 5), "nguong": "[0,0016; 0,023]"})
    add("2.6e", tong_hop(ca, False), "Thể tích nhân / thể tích phía trong màng ∈ [0,16 %; 2,3 %].", ca)
    ca = mat_nhieu(tb, ["nhan"])
    for b in NH:
        e = b.kich_thuoc()
        ca.append({"bo_phan": b.ten, "dat": all(trong_khoang(x, *DK_NHAN) for x in e) and e[2] / e[0] >= TV_NG_CAU, "gia_tri": [round(x, 3) for x in e], "nguong": "moi canh %s µm (tv.duong_kinh_nhan), canh nho/lon >= %g (hinh cau)" % (DK_NHAN, TV_NG_CAU)})
    add("2.6f", tong_hop(ca, False), "Hình cầu (cạnh nhỏ/lớn của OBB ≥ TV_NG_CAU) đường kính (mọi cạnh OBB) ∈ [5,2; 7] µm.", ca)

    _t('2.7b / 2.7c / 2.7d / 2.7e / 2.7f')
    ca = mat_nhieu(tb, ["luc_lap"]); mats = set()
    for b in LL:
        d, g = khe_mang.get(b.ten, (None, None))
        _, fi, _ = IN_L.aabb().squared_distance(IN_L.V, IN_L.F, b.luoi.V.mean(0)[None])
        nrm = trimesh.Trimesh(IN_L.V, IN_L.F, process=False).face_normals[fi[0]]; k = int(np.argmax(np.abs(nrm))); mats.add((k, nrm[k] > 0))
        ca.append({"bo_phan": b.ten, "dat": d is not None and trong_khoang(d, *KHE_LL_MANG), "gia_tri": round(d, 4) if d is not None else None, "nguong": "khe toi mang %s µm (tv.khe_luc_lap_mang)" % (KHE_LL_MANG,)})
    phu = sum(np.pi * (b.kich_thuoc()[0] / 2) * (b.kich_thuoc()[1] / 2) for b in LL) / max(1e-9, float(trimesh.Trimesh(IN_L.V, IN_L.F, process=False).area))
    ca.append({"bo_phan": "luc_lap", "dat": len(mats) >= TV_SO_MAT_RAI, "gia_tri": {"so_mat": len(mats), "thong_tin_ti_le_phu_mat_trong": round(phu, 3)}, "nguong": "rai tren >= %d / 6 mat hop (ti le phu chi la thong tin; luat ghi mo hinh phu 9–28 %%)" % TV_SO_MAT_RAI})
    add("2.7b", tong_hop(ca, False), "Mỗi lục lạp sát màng: khe tới mặt trong màng ∈ [0,3; 0,5] µm (Phụ lục 1b); 'rải trên mặt trong' = có lục lạp trên ≥ TV_SO_MAT_RAI trong 6 mặt hộp (tỉ lệ phủ mặt trong ghi thông tin); 'mặt dẹt hướng ra màng' kiểm ở 2.7f.", ca)
    ca = mat_nhieu(tb, ["luc_lap"]); ca_d = mat_nhieu(tb, ["luc_lap"])
    for b in LL:
        e = b.kich_thuoc(); hull = float(b.tm.convex_hull.volume)
        prof = b.the_tich / (np.pi * (e[0] / 2) * (e[1] / 2) * e[2])
        lom = do_lom(b)
        ca.append({"bo_phan": b.ten, "dat": b.the_tich / hull >= TV_NG_LOI and lom <= TV_NEP_SAU and e[1] / e[0] >= TV_NG_TRON_DIA and trong_khoang(prof, *TV_THAU_KINH),
                   "gia_tri": {"loi": round(b.the_tich / hull, 4), "lom_sau_um": round(lom, 4), "tron": round(e[1] / e[0], 3), "V_tren_dia_tru": round(prof, 3)}, "nguong": "loi >= %g va lom <= %g µm, e2/e1 >= %g, V/(pi a b t) trong %s" % (TV_NG_LOI, TV_NEP_SAU, TV_NG_TRON_DIA, TV_THAU_KINH)})
        ca_d.append({"bo_phan": b.ten, "dat": trong_khoang(e[0], *KT_LL_VE["duong_kinh"]) and trong_khoang(e[1], *KT_LL_VE["duong_kinh"]) and trong_khoang(e[2], *KT_LL_VE["day"]), "gia_tri": [round(x, 3) for x in e], "nguong": "duong kinh %s, day %s µm (tv.kich_thuoc_luc_lap_ve)" % (KT_LL_VE["duong_kinh"], KT_LL_VE["day"])})
    add("2.7c", tong_hop(ca, False), "Thấu kính hai mặt lồi: khối lồi (V/bao lồi ≥ TV_NG_LOI và độ lõm ≤ TV_NEP_SAU), đĩa tròn (e2/e1 ≥ TV_NG_TRON_DIA) và hồ sơ thấu kính V/(π·(e1/2)(e2/2)·e3) ∈ TV_THAU_KINH (thấu kính ≈ 0,5; elip 0,667; đĩa trụ 1).", ca)
    add("2.7d", tong_hop(ca_d, False), "Hai cạnh lớn OBB (đường kính) ∈ [5; 7] µm, cạnh nhỏ (bề dày) ∈ [2; 2,5] µm.", ca_d)
    ca = mat_nhieu(tb, []) + [{"bo_phan": "luc_lap", "dat": trong_khoang(len(LL), *SO_LL), "gia_tri": len(LL), "nguong": "%s (tv.so_luc_lap)" % (SO_LL,)}]
    add("2.7e", tong_hop(ca, False), "Số lục lạp có mặt ∈ [30; 50].", ca)
    ca = mat_nhieu(tb, ["luc_lap"])
    for b in LL:
        _, U = truc_pca(b.luoi.V); ax = U[:, 0]
        _, fi, _ = IN_L.aabb().squared_distance(IN_L.V, IN_L.F, b.luoi.V.mean(0)[None])
        nrm = trimesh.Trimesh(IN_L.V, IN_L.F, process=False).face_normals[fi[0]]
        g = float(np.degrees(np.arccos(min(1.0, abs(float(ax @ nrm))))))
        ca.append({"bo_phan": b.ten, "dat": g <= LECH_LL, "gia_tri": round(g, 2), "nguong": "lech truc ngan / phap tuyen mang <= %g° (tv.lech_huong_luc_lap)" % LECH_LL})
    add("2.7f", tong_hop(ca, False), "Trục ngắn nhất (PCA) của lục lạp lệch pháp tuyến mặt màng gần nhất ≤ 30° (Phụ lục 1b). Luật mềm.", ca)

    _t('2.8 / 2.9 (lớp 6)')
    ca = [{"bo_phan": "ti_the", "dat": len(TT) == 0 or trong_khoang(len(TT), *TT_LOP6), "gia_tri": len(TT), "nguong": "0 (tat) hoac %s (tv.ti_the_lop6)" % (TT_LOP6,)}]
    add("2.8", tong_hop(ca, False), "Lớp 6: ti thể tắt (0) hoặc 5–10.", ca)
    ca = []
    for lo_ in ("mang_trong_lap", "granum", "thylakoid", "tam_noi", "luoi_noi_chat", "golgi", "tui_golgi", "peroxisome", "ribosome"):
        ca.append({"bo_phan": lo_, "dat": len(tb.tl.get(lo_, [])) == 0, "gia_tri": len(tb.tl.get(lo_, [])), "nguong": "0 o lop 6"})
    add("2.9", tong_hop(ca, True), "Lớp 6: lục lạp không có cấu tạo bên trong (không màng trong, granum, thylakoid, tấm nối); không có lưới nội chất, Golgi, peroxisome (và ribosome). 'Màu xanh đặc' không kiểm được (không có dữ liệu màu).", ca, ly_do="Màu không có trong mesh.")

    _t('K.1')
    ca = []; MIEN = set(MIEN_K1) | {frozenset(["ribosome", "nhan"])}
    def _ng(a, b): return KCL_TV.get(a.loai, KHE_ROI) if a.loai == b.loai else KHE_ROI
    for b in tren:
        if b.loai in ("khong_bao",): continue
        d, g = khe_mang.get(b.ten, (None, None))
        if d is not None:
            ca.append({"bo_phan": [b.ten, "mang"], "dat": d >= KHE_ROI - 1e-9 and g <= EPS_VOL, "gia_tri": round(d, 4), "nguong": ">= %g µm toi mang" % KHE_ROI})
    for a, b in cap_ung_vien(tren, tren, KHE_ROI, cung=True):
        if frozenset([a.loai, b.loai]) in MIEN or a.loai == b.loai == "ribosome":
            continue
        ng_ = _ng(a, b)
        d, g, cx = khe(a, b, max(ng_, 1e-6))
        ca.append({"bo_phan": [a.ten, b.ten], "dat": d >= ng_ - 1e-9 and g <= EPS_VOL, "gia_tri": round(d, 4), "giao_um3": round(g, 6), "nguong": ">= %g µm" % ng_})
    for n_ in NH:
        ncs = tb.con_loai(n_, "nhan_con"); env = tb.con_loai(n_, "mang_nhan")
        for c in ncs:
            for e_ in env:
                d, g, cx = khe(c, e_, KHE_ROI)
                ca.append({"bo_phan": [c.ten, e_.ten], "dat": d >= KHE_ROI - 1e-9, "gia_tri": round(d, 4), "nguong": ">= %g µm (nhan con – mang nhan)" % KHE_ROI})
        _n = KCL_TV.get("nhan_con", KHE_ROI)
        for a, b in cap_ung_vien(ncs, ncs, max(_n, 1e-9), cung=True):
            d, g, cx = khe(a, b, max(_n, 1e-6))
            ca.append({"bo_phan": [a.ten, b.ten], "dat": d >= _n - 1e-9 and g <= EPS_VOL, "gia_tri": round(d, 4), "nguong": ">= %g µm" % _n})
    for l_ in LL:
        gr = tb.con_loai(l_, "granum"); tn = tb.con_loai(l_, "tam_noi"); mt = tb.con_loai(l_, "mang_trong_lap")
        for a, b in cap_ung_vien(gr, gr, KCL_TV.get("granum", KHE_TRONG_LL), cung=True) + cap_ung_vien(tn, tn, KCL_TV.get("tam_noi", KHE_TRONG_LL), cung=True):
            ng_ = KCL_TV.get(a.loai, KHE_TRONG_LL)
            d, g, cx = khe(a, b, ng_)
            ca.append({"bo_phan": [a.ten, b.ten], "dat": d >= ng_ - 1e-9, "gia_tri": round(d, 4), "nguong": ">= %g µm (trong luc lap)" % ng_})
        for m_ in mt:
            Vm = Vo(m_)
            for x in gr + tn:
                ok, _ = trong(x, m_.man, m_.luoi)
                d = float(np.sqrt(m_.luoi.kc2_diem(x.luoi.V)).min()); cx = d < KHE_TRONG_LL
                ca.append({"bo_phan": [x.ten, m_.ten], "dat": ok and (not cx or d >= KHE_TRONG_LL - 1e-9), "gia_tri": round(d, 4) if cx else ">=%g" % KHE_TRONG_LL, "nguong": ">= %g µm toi mang trong, ben trong (tv.khe_trong_luc_lap)" % KHE_TRONG_LL})
    add("K.1", tong_hop(ca, False),
        "Khe ≥ 0,3 µm: mọi bộ phận trên cùng (nhân, lục lạp, ti thể, Golgi, peroxisome, lưới nội chất, ribosome, không bào) với mặt trong màng và với nhau (kể cả cùng loại); nhân con với màng nhân và với nhau. Trong lục lạp: granum–granum, tấm nối–tấm nối, granum/tấm nối với màng trong ≥ 0,1 µm. MIỄN (theo lời luật): lưới nội chất–nhân (nối liền màng nhân), ribosome–lưới nội chất, ribosome–nhân (bám màng ngoài nhân), tấm nối–granum. Không xét cặp tổ tiên–con cháu và bên trong một granum/một chồng Golgi/một ti thể (khe riêng ở 10.2d, 5.6c).",
        ca, bat_dong="Không bào–màng theo 2.5b (≥ 1 µm, chặt hơn).")

    _t('M.a')
    ca = []
    for loai_ in sorted(set(tb.du_kien.values()) - {"vach", "mang"}):
        ca += ca_mat(tb, loai_)
    for b in tb.bp.values():
        if b.loai in ("vach", "mang") or IN_man is None: continue
        if b.loai in ("thylakoid", "lo_nhan", "mao", "tui_golgi", "tam_noi", "granum", "mang_trong_lap", "nhan_con", "mang_nhan"):
            continue
        ok, hv = trong(b, IN_man, IN_L)
        ca.append({"bo_phan": b.ten, "dat": ok, "gia_tri": round(hv, 6), "nguong": "ben trong mat trong mang"})
    vo_cache = {}
    for c_, p_ in tb.cha.items():
        if c_ in tb.bp and p_ in tb.bp:
            C, Pp = tb.bp[c_], tb.bp[p_]
            if Pp.loai == "mang_nhan" and p_ in tb.cha and tb.cha[p_] in tb.bp:
                p_ = tb.cha[p_]; Pp = tb.bp[p_]
            if p_ not in vo_cache: vo_cache[p_] = Vo(Pp)
            Vp = vo_cache[p_]; M_ = Vp.fill if Vp.fill is not None else Pp.man; Lp = Vp.LO if Vp.LO is not None else Pp.luoi
            ok, hv = trong(C, M_, Lp)
            if not ok:
                ca.append({"bo_phan": [c_, p_], "dat": False, "gia_tri": round(hv, 6), "nguong": "con nam trong cha"})
    add("M.a", tong_hop(ca, False), "Mọi bộ phận trên cùng nằm trong khối phía trong màng; mọi bộ phận con nằm trong khối (bao ngoài) của cha. Bộ phận bị mất = trượt.", ca)

    _t('M.d')
    ca = mat_nhieu(tb, ["ti_the", "golgi", "peroxisome"])
    for b in TT + GO + PE:
        dm = khe_mang.get(b.ten, (99, 0))[0]
        dn = min([khe(b, n_, GAN_MANG_NHAN + 1e-6)[0] for n_ in NH] or [99])
        ca.append({"bo_phan": b.ten, "dat": min(dm, dn) <= GAN_MANG_NHAN + 1e-9, "gia_tri": {"khe_mang": round(dm, 4), "khe_nhan": round(dn, 4)}, "nguong": "cach mang hoac nhan <= %g µm (tv.gan_mang_hoac_nhan)" % GAN_MANG_NHAN})
    add("M.d", tong_hop(ca, False), "Ti thể, Golgi, peroxisome cách mặt trong màng hoặc nhân ≤ 1 µm (Phụ lục 1b).", ca)

    _t('M.e / M.g / M.h')
    ca = []; cag = []
    for loai_ in sorted(set(tb.du_kien.values()) - LIEN_KHOI_MIEN_TV):
        ca += ca_mat(tb, loai_)
    for b in tb.bp.values():
        if b.loai in LIEN_KHOI_MIEN_TV: continue
        ca.append({"bo_phan": b.ten, "dat": b.so_khoi == 1 and not b.loi_chan, "gia_tri": b.so_khoi if not b.loi_chan else {"so_khoi": b.so_khoi, "loi_chan": b.loi_chan}, "nguong": "mot khoi lien hop le"})
    add("M.e", tong_hop(ca, False), "Mỗi bộ phận (trừ Golgi, lưới nội chất, ribosome) là đúng một khối đặc hợp lệ (vỏ có khoang tính là 1).", ca)
    for b in tb.bp.values():
        if b.loai == "khong_ro":
            cag.append({"bo_phan": b.ten, "dat": False, "gia_tri": "loai khong ro", "nguong": "moi vat co loai"})
        if b.loai not in ("ribosome", "luoi_noi_chat") and b.so_khoi > 1:
            cag.append({"bo_phan": b.ten, "dat": False, "gia_tri": b.so_khoi, "nguong": "mot vat = mot bao quan (khong gop nhieu ban sao vao mot vat)"})
    cag.append({"bo_phan": "*", "dat": True, "gia_tri": len(tb.bp), "nguong": "so vat rieng"})
    r = ket_qua(L_("M.g"), tong_hop(cag, False), "Mỗi bào quan là một vật riêng trong mesh.npz (một khoá riêng, có loại); không vật nào gộp nhiều bản sao (> 1 khối), trừ ribosome (tập chấm) và lưới nội chất (vùng). Báo riêng: không tính vào lỗi sinh học.", cag)
    r["kiem_rieng_day_hoc"] = True; KQ.append(r)
    ca = []; nhom = collections.defaultdict(list)
    for b in tb.bp.values():
        if b.loai not in ("vach", "mang"): nhom[b.loai].append(b)
    for lo_, xs in nhom.items():
        if lo_ == "ribosome": continue
        for a, b in cap_ung_vien(xs, xs, 1e-9, cung=True):
            g = giao(a, b)
            if g > EPS_VOL:
                ca.append({"bo_phan": [a.ten, b.ten], "dat": False, "gia_tri": round(g, 7), "nguong": "giao <= %g µm³" % EPS_VOL})
    add("M.h", tong_hop(ca, False), "Mọi cặp cùng loại (mọi loại, kể cả thylakoid trong một granum, túi trong một chồng): thể tích giao chính xác ≤ EPS_VOL; chạm được phép.", ca)

    _t('5.1a / 5.1b / 5.2a / 5.2b / 5.2c')
    NC = tb.tl.get("nhan_con", [])
    ca = mat_nhieu(tb, ["nhan_con"]) + [{"bo_phan": "nhan_con", "dat": trong_khoang(len(NC), *SO_NHAN_CON), "gia_tri": len(NC), "nguong": "%s (tv.so_nhan_con)" % (SO_NHAN_CON,)}]
    for c in NC:
        p = tb.bp.get(tb.cha.get(c.ten, ""), None)
        ok = p is not None and p.loai == "nhan" and trong(c, p.man, p.luoi)[0]
        ca.append({"bo_phan": c.ten, "dat": ok, "gia_tri": "trong nhan" if ok else "ngoai nhan / khong cha", "nguong": "nam trong nhan"})
    add("5.1a", tong_hop(ca, True), "1–2 nhân con, mỗi cái nằm trong nhân. 'Không có màng bao' không kiểm được từ mesh.", ca, ly_do="Không phân biệt được có/không màng bao nhân con.")
    ca = mat_nhieu(tb, ["nhan_con"])
    for c in NC:
        e = c.kich_thuoc(); ca.append({"bo_phan": c.ten, "dat": all(trong_khoang(x, *DK_NHAN_CON) for x in e), "gia_tri": [round(x, 3) for x in e], "nguong": "%s µm (tv.duong_kinh_nhan_con)" % (DK_NHAN_CON,)})
    add("5.1b", tong_hop(ca, False), "Mọi cạnh OBB của nhân con ∈ [1; 2] µm.", ca)
    tia = None
    if mn is not None and NH:
        c0 = NH[0].luoi.V.mean(0); dirs = kiem2._fibo(400)
        hits = mn.luoi.aabb().intersect_ray(mn.luoi.V, mn.luoi.F, np.repeat(c0[None], len(dirs), 0), dirs)
        tia = []
        lo_dir = []
        for p in tb.tl.get("lo_nhan", []):
            v = p.luoi.V.mean(0) - c0; R_ = np.linalg.norm(v); w = np.ptp(p.luoi.V, 0).max() / 2
            lo_dir.append((v / R_, (w + 0.1) / R_))
        for d_, h_ in zip(dirs, hits):
            if any(np.arccos(np.clip(d_ @ u, -1, 1)) < a for u, a in lo_dir):
                continue
            ts = sorted(set(round(x[1], 7) for x in h_))
            tia.append(ts)
    ca = mat_nhieu(tb, ["mang_nhan", "lo_nhan"])
    LO = tb.tl.get("lo_nhan", [])
    if tia is not None:
        bon = [t for t in tia if len(t) == 4]
        f = len(bon) / max(1, len(tia))
        ca.append({"bo_phan": mn.ten, "dat": f >= TV_MD_TIA, "gia_tri": round(f, 3), "nguong": ">= %g tia thay hai lop (4 mat)" % TV_MD_TIA})
        kh = [t[2] - t[1] for t in bon]
        ca.append({"bo_phan": mn.ten, "dat": bool(kh) and min(kh) > 1e-4, "gia_tri": round(min(kh), 4) if kh else None, "nguong": "co khoang giua hai lop (> 0)"})
        g = int(round((2 - max(mn.comps, key=lambda c: c.volume).euler_number) / 2))
        ca.append({"bo_phan": mn.ten, "dat": g >= 1 and len(LO) >= 1, "gia_tri": {"giong": g, "so_lo_nhan": len(LO)}, "nguong": "co lo nhan"})
        er_ok = any(dien_tich_cham(e_, VN.LO) > 0 for e_ in ER) if VN is not None else False
        ca.append({"bo_phan": mn.ten, "dat": er_ok, "gia_tri": er_ok, "nguong": "mang ngoai noi lien luoi noi chat (cham)"})
    else:
        ca.append({"bo_phan": "mang_nhan", "dat": False, "gia_tri": "khong co", "nguong": "phai co"})
    add("5.2a", tong_hop(ca, False), "Màng nhân: theo 400 tia từ tâm nhân, ≥ TV_MD_TIA tia gặp đúng 4 mặt (hai lớp) và có khoang > 0 giữa hai lớp; có lỗ (giống ≥ 1 và có bộ phận lỗ nhân); màng ngoài chạm lưới nội chất. Số thật 20–50 nm chỉ là thông tin (mô hình vẽ theo 5.2c).", ca)
    ca = mat_nhieu(tb, ["lo_nhan"]) + [{"bo_phan": "lo_nhan", "dat": trong_khoang(len(LO), *LO_NHAN_VE["so"]), "gia_tri": len(LO), "nguong": "%s (tv.lo_nhan_ve)" % (LO_NHAN_VE["so"],)}]
    if NH:
        c0 = NH[0].luoi.V.mean(0)
        for p in LO:
            pc = p.luoi.V.mean(0); rr = (pc - c0) / np.linalg.norm(pc - c0)
            a1 = np.cross(rr, [1, 0, 0] if abs(rr[0]) < 0.9 else [0, 1, 0]); a1 /= np.linalg.norm(a1); a2 = np.cross(rr, a1)
            w1 = np.ptp(p.luoi.V @ a1); w2 = np.ptp(p.luoi.V @ a2)
            ok = trong_khoang(w1, *LO_NHAN_VE["duong_kinh"]) and trong_khoang(w2, *LO_NHAN_VE["duong_kinh"]) and (mn is None or giao(p, mn) <= EPS_VOL)
            ca.append({"bo_phan": p.ten, "dat": ok, "gia_tri": [round(w1, 3), round(w2, 3)], "nguong": "rong [0,2; 0,3] µm (vuong goc huong ban kinh), nam trong lo (khong giao mang nhan)"})
    add("5.2b", tong_hop(ca, False), "30–50 lỗ nhân; bề rộng mỗi lỗ (đo vuông góc hướng bán kính) ∈ [0,2; 0,3] µm; lỗ không giao vật liệu màng nhân.", ca)
    ca = mat_nhieu(tb, ["mang_nhan"])
    if tia is not None:
        bon = [t for t in tia if len(t) == 4]
        if bon:
            tong = [t[3] - t[0] for t in bon]; l1 = [t[1] - t[0] for t in bon]; l2 = [t[3] - t[2] for t in bon]; kh = [t[2] - t[1] for t in bon]
            ca.append({"bo_phan": mn.ten, "dat": min(tong) >= MANG_NHAN_VE["day"][0] - 1e-4 and max(tong) <= MANG_NHAN_VE["day"][1] + 1e-4, "gia_tri": [round(min(tong), 4), round(max(tong), 4)], "nguong": "day [0,15; 0,2]"})
            ok = all(abs(x - MANG_NHAN_VE["day_lop"]) <= TV_TOL_KHE_LAP for x in l1 + l2)
            ca.append({"bo_phan": mn.ten, "dat": ok, "gia_tri": [round(min(l1 + l2), 4), round(max(l1 + l2), 4)], "nguong": "moi lop 0,05 ± %g" % TV_TOL_KHE_LAP})
            ca.append({"bo_phan": mn.ten, "dat": min(kh) >= MANG_NHAN_VE["khoang"][0] - 1e-4 and max(kh) <= MANG_NHAN_VE["khoang"][1] + 1e-4, "gia_tri": [round(min(kh), 4), round(max(kh), 4)], "nguong": "khoang [0,05; 0,1]"})
        else:
            ca.append({"bo_phan": mn.ten, "dat": False, "gia_tri": "khong tia nao thay 2 lop", "nguong": ""})
    else:
        ca.append({"bo_phan": "mang_nhan", "dat": False, "gia_tri": "khong co", "nguong": ""})
    add("5.2c", tong_hop(ca, False), "Theo tia (không qua lỗ): tổng bề dày màng nhân ∈ [0,15; 0,2] µm, mỗi lớp 0,05 ± TV_TOL_KHE_LAP µm, khoang ∈ [0,05; 0,1] µm.", ca)

    _t('5.3b / 5.3c / 5.3e')
    ca = mat_nhieu(tb, ["ti_the"])
    for b in TT:
        e = b.kich_thuoc(); w, sw = b.do_day(); Lg, _ = kiem2.dai_truc_giua(b)
        dai = Lg if (Lg is not None and Lg > e[0] * (1 + kiem2.NG_CONG_43B)) else e[0]
        _r0, _r1 = KT_TT_VE["rong"]
        ok = trong_khoang(dai, *KT_TT_VE["dai"]) and (trong_khoang(w, _r0, _r1) or (w < _r0 and w + sw >= _r0))
        ca.append({"bo_phan": b.ten, "dat": ok, "gia_tri": {"dai": round(dai, 4), "rong_day_cuc_bo": round(w, 4), "sai_so": round(sw, 4), "dat_trong_sai_so": bool(ok and w < 0.5)}, "nguong": "dai [1;2], rong (duong kinh mat cat = do day cuc bo) [0,5;1]"})
    add("5.3b", tong_hop(ca, False), "Như 4.3b tế bào động vật (PI đã quyết): dài = OBB hoặc trục giữa nếu cong; rộng = độ dày cục bộ (2 × bán kính nội tiếp), sai số đo nghiêng về ĐẠT.", ca)
    ca = [{"bo_phan": "ti_the", "dat": trong_khoang(len(TT), *SO_TT), "gia_tri": len(TT), "nguong": "%s (tv.so_ti_the)" % (SO_TT,)}]
    add("5.3c", tong_hop(ca, False), "Số ti thể có mặt ∈ [10; 20] (lớp 10).", ca)
    ca = mat_nhieu(tb, ["mao"])
    for b in TT:
        ms = tb.con_loai(b, "mao")
        ca.append({"bo_phan": b.ten, "dat": trong_khoang(len(ms), *MAO_VE["so"]), "gia_tri": len(ms), "nguong": "%s mao (tv.mao_ve)" % (MAO_VE["so"],)})
        for m_ in ms:
            e = m_.kich_thuoc()
            ca.append({"bo_phan": m_.ten, "dat": trong_khoang(e[2], *MAO_VE["day"]) and trong(m_, b.man, b.luoi)[0], "gia_tri": round(e[2], 4), "nguong": "day [0,03; 0,05], ben trong ti the"})
    add("5.3e", tong_hop(ca, False), "Mỗi ti thể có 4–8 mào (bộ phận con 'mao'), mỗi mào dày (cạnh OBB nhỏ nhất) ∈ [0,03; 0,05] µm và nằm trong ti thể.", ca)

    _t('5.4a / 5.4b / 5.5')
    ca = mat_nhieu(tb, ["luoi_noi_chat"])
    for e_ in ER:
        S = dien_tich_cham(e_, VN.LO) if VN is not None else 0.0
        ca.append({"bo_phan": e_.ten, "dat": S >= VUNG_ER, "gia_tri": round(S, 3), "nguong": "vung noi voi mang ngoai nhan >= %g µm² (tv.vung_noi_er_mang_nhan)" % VUNG_ER})
    if not ER: ca.append({"bo_phan": "luoi_noi_chat", "dat": False, "gia_tri": "khong co", "nguong": "phai co"})
    add("5.4a", tong_hop(ca, True), "Lưới nội chất chạm mặt ngoài màng nhân (diện tích mặt ER cách ≤ TV_TIEP_XUC) ≥ 1 µm². Hạt/trơn (có ribosome hay không) không phân biệt được trên một vùng ER.", ca, ly_do="Không phân biệt ER hạt và ER trơn trong một bộ phận.")
    ca = mat_nhieu(tb, ["luoi_noi_chat"])
    for e_ in ER:
        ca.append({"bo_phan": e_.ten, "dat": e_.so_khoi == 1, "gia_tri": e_.so_khoi, "nguong": "lien khoi"})
        if not NH or IN_L is None: continue
        nh_ = NH[0]; c0 = nh_.luoi.V.mean(0); dirs = kiem2._fibo(300)
        hn = nh_.luoi.aabb().intersect_ray(nh_.luoi.V, nh_.luoi.F, np.repeat(c0[None], len(dirs), 0), dirs)
        he = e_.luoi.aabb().intersect_ray(e_.luoi.V, e_.luoi.F, np.repeat(c0[None], len(dirs), 0), dirs)
        khac = [b for b in tren if b.loai not in ("nhan", "luoi_noi_chat", "khong_bao", "ribosome")]
        def bi_chan(p):
            if np.sqrt(IN_L.kc2_diem(p[None]))[0] <= ER_VE["cach_mang_um"] + TV_TOL_ER: return True
            return any(np.sqrt(b.luoi.kc2_diem(p[None]))[0] <= ER_VE["cach_bo_phan_um"] + TV_TOL_ER for b in khac if np.all(p >= b.lo - 0.5) and np.all(p <= b.hi + 0.5))
        sai = collections.Counter(); day = []
        for d, h1, h2 in zip(dirs, hn, he):
            if not h1: continue
            tn = max(x[1] for x in h1)
            te = sorted(x[1] for x in h2 if x[1] >= tn - TV_TIEP_XUC)
            if te and abs(te[0] - tn) <= TV_TIEP_XUC:
                t_ra = te[1] if len(te) > 1 else te[0]
                th = t_ra - tn; day.append(th)
                p = c0 + d * t_ra
                if th > ER_VE["day_quanh_nhan_um"][1] + TV_TOL_ER:
                    q = c0 + d * (tn + ER_VE["day_quanh_nhan_um"][1] + TV_TOL_ER)
                    if np.sqrt(e_.luoi.kc2_diem(q[None]))[0] > 0.1 + TV_TOL_ER: sai["vanh_day_qua_1um"] += 1
                elif th < ER_VE["day_quanh_nhan_um"][0] - TV_TOL_ER and not bi_chan(p):
                    sai["vanh_mong_khong_ly_do"] += 1
            else:
                if not bi_chan(c0 + d * (tn + 0.02)): sai["thieu_vanh"] += 1
        ca.append({"bo_phan": e_.ten, "dat": sum(sai.values()) == 0, "gia_tri": {"loi_tia": dict(sai), "day_vanh_min_max": [round(min(day), 3), round(max(day), 3)] if day else None}, "nguong": "vanh [0,5;1] µm quanh nhan (mong/thieu chi khi bi cat boi vung 0,3 µm quanh mang/bao quan), dung sai %g" % TV_TOL_ER})
        dn = np.sqrt(nh_.luoi.kc2_diem(e_.luoi.V)); idx = np.nonzero(dn > ER_VE["day_quanh_nhan_um"][1] + TV_TOL_ER)[0]
        if len(idx):
            idx = idx[:: max(1, len(idx) // 1500)]
            nr = trimesh.Trimesh(e_.luoi.V, e_.luoi.F, process=False).vertex_normals[idx]
            O = e_.luoi.V[idx] - nr * 1e-5
            _, tt_, _ = e_.luoi.aabb().intersect_ray_first(e_.luoi.V, e_.luoi.F, O, -nr)
            dk = tt_[np.isfinite(tt_)]
            dm = np.sqrt(IN_L.kc2_diem(e_.luoi.V[idx]))
            med = float(np.median(dk)) if len(dk) else float('nan')
            ca.append({"bo_phan": e_.ten, "dat": bool(len(dk)) and dk.min() >= ER_VE["ong_duong_kinh_um"][0] - TV_TOL_ER / 5 and med <= ER_VE["ong_duong_kinh_um"][1] + TV_TOL_ER / 5,
                       "gia_tri": {"duong_kinh_ong_min_trungvi": [round(float(dk.min()), 4), round(med, 4)] if len(dk) else None}, "nguong": "ong duong kinh [0,1; 0,2] µm: min va trung vi cua tia phap tuyen vao trong (tia o dau ong/cho noi dai hon nen dung trung vi)"})
            ca.append({"bo_phan": e_.ten, "dat": dm.min() >= ER_VE["ong_do_sau_um"][0] - 1e-6 and dm.max() <= ER_VE["ong_do_sau_um"][1] + ER_VE["ong_duong_kinh_um"][1] / 2, "gia_tri": [round(float(dm.min()), 4), round(float(dm.max()), 4)], "nguong": "mat ong cach mang [0,3; 0,7] (truc ong o do sau 0,3–0,6)"})
    add("5.4b", tong_hop(ca, False), "Lưới nội chất liền khối; theo 300 tia từ tâm nhân: vành dày [0,5; 1] µm áp nhân (mỏng/thiếu chỉ được khi bị cắt bởi vùng 0,3 µm quanh màng hay bào quan khác; phần dài quá 1 µm chỉ được nếu là ống); phần cách nhân > 1 µm là ống đường kính [0,1; 0,2] µm, trục ở độ sâu 0,3–0,6 µm.", ca)
    if not RB:
        add("5.5", "KHONG_AP_DUNG", "Tuỳ chọn giáo viên bật/tắt: mô hình không có ribosome → coi là TẮT.", [], ly_do="Không có bộ phận ribosome (tuỳ chọn tắt).")
    else:
        ca = []
        for r_ in RB:
            comps = r_.comps
            ca.append({"bo_phan": r_.ten, "dat": trong_khoang(len(comps), *RB_VE["so_cham"]), "gia_tri": len(comps), "nguong": "3 000–5 000 cham (Phu luc 1b)"})
            Ls = [VN.LO] if VN is not None else []
            Ls += [e_.luoi for e_ in ER]
            xa = 0; sai_co = 0
            for c in comps:
                ext = np.ptp(c.vertices, 0)
                if not (RB_VE["duong_kinh_um"][0] - 1e-3 <= ext.max() <= RB_VE["duong_kinh_um"][1] + 1e-3): sai_co += 1
                cc = c.vertices.mean(0)[None]; rr = ext.max() / 2
                dmin = min(float(np.sqrt(L.kc2_diem(cc))[0]) for L in Ls) if Ls else 9
                if dmin > rr + TV_TIEP_XUC: xa += 1
            ca.append({"bo_phan": r_.ten, "dat": sai_co == 0, "gia_tri": sai_co, "nguong": "so cham sai co (rong [0,03; 0,04])"})
            ca.append({"bo_phan": r_.ten, "dat": xa == 0, "gia_tri": xa, "nguong": "so cham khong bam ER / mang ngoai nhan"})
        add("5.5", tong_hop(ca, False), "Nếu bật: 3 000–5 000 chấm rộng 0,03–0,04 µm (Phụ lục 1b), mỗi chấm chạm (≤ TV_TIEP_XUC) lưới nội chất hoặc mặt ngoài màng nhân.", ca)

    _t('5.6a / 5.6b / 5.6c')
    ca = mat_nhieu(tb, ["golgi"])
    for a, b in cap_ung_vien(GO, GO, KHE_GOLGI, cung=True):
        d, g, cx = khe(a, b, KHE_GOLGI)
        ca.append({"bo_phan": [a.ten, b.ten], "dat": d >= KHE_GOLGI - 1e-9, "gia_tri": round(d, 4), "nguong": ">= %g µm (tv.khe_golgi_golgi)" % KHE_GOLGI})
    if GO and NH:
        dn = [khe(g_, NH[0], 50.0)[0] for g_ in GO]
        f_gan = float(np.mean([d <= TV_GOLGI_GAN_NHAN for d in dn]))
        ca.append({"bo_phan": "golgi", "dat": f_gan <= TV_GOLGI_TI_LE_GAN, "gia_tri": round(f_gan, 3), "nguong": "khong gom canh nhan: ti le the Golgi cach nhan <= %g µm phai <= %g" % (TV_GOLGI_GAN_NHAN, TV_GOLGI_TI_LE_GAN)})
        c_tb = (mang.lo + mang.hi) / 2 if mang else np.zeros(3)
        vung = set(tuple((g_.luoi.V.mean(0) > c_tb).astype(int)) for g_ in GO)
        ca.append({"bo_phan": "golgi", "dat": len(vung) >= TV_GOLGI_SO_VUNG, "gia_tri": len(vung), "nguong": "rai khap: co mat o >= %d / 8 phan tam te bao" % TV_GOLGI_SO_VUNG})
    ca.append({"bo_phan": "golgi", "dat": len(GO) >= 2, "gia_tri": len(GO), "nguong": "nhieu the roi"})
    add("5.6a", tong_hop(ca, False), "Nhiều thể Golgi rời nhau, cách nhau ≥ 2 µm (Phụ lục 1b); 'không gom cạnh nhân' = tỉ lệ thể Golgi cách nhân ≤ TV_GOLGI_GAN_NHAN không quá TV_GOLGI_TI_LE_GAN; 'rải khắp chất tế bào' = có mặt ở ≥ TV_GOLGI_SO_VUNG trong 8 phần tám của tế bào.", ca)
    ca = [{"bo_phan": "golgi", "dat": trong_khoang(len(GO), *SO_GOLGI), "gia_tri": len(GO), "nguong": "%s (tv.so_golgi)" % (SO_GOLGI,)}]
    add("5.6b", tong_hop(ca, False), "Số thể Golgi ∈ [5; 10].", ca)
    ca = mat_nhieu(tb, ["tui_golgi"])
    for g_ in GO:
        tu = tb.con_loai(g_, "tui_golgi"); e = g_.kich_thuoc()
        ca.append({"bo_phan": g_.ten, "dat": trong_khoang(len(tu), *GOLGI_VE["so_tui"]), "gia_tri": len(tu), "nguong": "%s tui (tv.golgi_ve)" % (GOLGI_VE["so_tui"],)})
        ca.append({"bo_phan": g_.ten, "dat": trong_khoang(e[0], *GOLGI_VE["rong"]) and trong_khoang(e[1], *GOLGI_VE["rong"]) and trong_khoang(e[2], *GOLGI_VE["day_chong"]), "gia_tri": [round(x, 3) for x in e], "nguong": "rong [0,5;1], cao chong [0,18; 0,56]"})
        if tu:
            _, U = truc_pca(np.concatenate([t.luoi.V for t in tu])); ax = U[:, 0]
            tu = sorted(tu, key=lambda t: float(t.luoi.V.mean(0) @ ax))
            for t in tu:
                ca.append({"bo_phan": t.ten, "dat": trong_khoang(t.kich_thuoc()[2], *GOLGI_VE["day_tui"]), "gia_tri": round(t.kich_thuoc()[2], 4), "nguong": "tui day [0,03; 0,06]"})
            for a, b in zip(tu[:-1], tu[1:]):
                d, cx = kc(a, b, 1.0)
                ca.append({"bo_phan": [a.ten, b.ten], "dat": trong_khoang(d, GOLGI_VE["khe_tui"][0] - 1e-6, GOLGI_VE["khe_tui"][1] + 1e-6), "gia_tri": round(d, 4), "nguong": "khe giua tui %s (tv.golgi_ve)" % (GOLGI_VE["khe_tui"],)})
    add("5.6c", tong_hop(ca, False), "Mỗi thể Golgi: 4–6 túi con; bao chồng rộng (hai cạnh lớn) ∈ [0,5; 1], cao (cạnh nhỏ) ∈ [0,18; 0,56] µm; mỗi túi dày ∈ [0,03; 0,06]; khe chính xác giữa hai túi liên tiếp (theo trục chồng) ∈ [0,02; 0,04] µm.", ca)

    _t('5.7a / 5.7c / 5.7d')
    ca = mat_nhieu(tb, ["peroxisome"])
    for b in PE:
        e = b.kich_thuoc(); ca.append({"bo_phan": b.ten, "dat": all(trong_khoang(x, *DK_PE_VE) for x in e) and e[2] / e[0] >= TV_NG_CAU, "gia_tri": [round(x, 3) for x in e], "nguong": "%s µm (tv.duong_kinh_peroxisome_ve), tron (>= %g)" % (DK_PE_VE, TV_NG_CAU)})
    add("5.7a", tong_hop(ca, False), "Mọi cạnh OBB ∈ [0,5; 1] µm và tròn (cạnh nhỏ/lớn ≥ TV_NG_CAU).", ca)
    add("5.7c", tong_hop([{"bo_phan": "peroxisome", "dat": trong_khoang(len(PE), *SO_PE), "gia_tri": len(PE), "nguong": "20–40"}], False), "Số peroxisome ∈ [20; 40].", [{"bo_phan": "peroxisome", "dat": trong_khoang(len(PE), *SO_PE), "gia_tri": len(PE), "nguong": "20–40"}])
    ca = mat_nhieu(tb, ["peroxisome"])
    for b in PE:
        dl = min([khe(b, l_, GAN_PE_LL + 1e-6)[0] for l_ in LL] or [99])
        ca.append({"bo_phan": b.ten, "dat": dl <= GAN_PE_LL + 1e-9, "gia_tri": round(dl, 4), "nguong": "cach luc lap gan nhat <= %g µm (tv.gan_peroxisome_luc_lap)" % GAN_PE_LL})
    add("5.7d", tong_hop(ca, False), "Mỗi peroxisome cách lục lạp gần nhất ≤ 1 µm (Phụ lục 1b). Luật mềm.", ca)

    _t('5.8–5.11')
    for rid, lo_ in (("5.8", "lysosome"), ("5.9", "trung_the"), ("5.10", "cau_sinh_chat"), ("5.11", "khung_xuong")):
        ca = [{"bo_phan": lo_, "dat": len(tb.tl.get(lo_, [])) == 0, "gia_tri": len(tb.tl.get(lo_, [])), "nguong": "0"}]
        add(rid, tong_hop(ca, False), "Không có bộ phận loại '%s'." % lo_, ca)

    _t('10.x')
    ca1 = mat_nhieu(tb, ["mang_trong_lap"]); cb1 = []; c2a = mat_nhieu(tb, ["granum", "thylakoid"]); c2b = []; c2c = []; c2d = []
    c3a = []; c3b = []; c3c = []; c5 = []
    for l_ in LL:
        mt = tb.con_loai(l_, "mang_trong_lap")
        ca1.append({"bo_phan": l_.ten, "dat": len(mt) == 1, "gia_tri": len(mt), "nguong": "1 mang trong"})
        for m_ in mt:
            ok, _ = trong(m_, l_.man, l_.luoi)
            dmin, cx = kc_mat_mat(m_.luoi, l_.luoi, 1.0)
            dmax = float(np.sqrt(l_.luoi.kc2_diem(m_.luoi.V)).max())
            ca1.append({"bo_phan": m_.ten, "dat": ok and abs(dmin - KHE_MT_LAP) <= TV_TOL_KHE_LAP and abs(dmax - KHE_MT_LAP) <= TV_TOL_KHE_LAP, "gia_tri": [round(dmin, 4), round(dmax, 4)], "nguong": "khe hai mang 0,05 ± %g (Phu luc 1b), ben trong" % TV_TOL_KHE_LAP})
            hull = float(m_.tm.convex_hull.volume); lom = do_lom(m_)
            cb1.append({"bo_phan": m_.ten, "dat": m_.the_tich / hull >= TV_NG_LOI and lom <= TV_NEP_SAU, "gia_tri": {"V_tren_bao_loi": round(m_.the_tich / hull, 4), "lom_sau_um": round(lom, 4)}, "nguong": "khong gap nep: V/bao loi >= %g va do lom <= %g µm" % (TV_NG_LOI, TV_NEP_SAU)})
        gr = tb.con_loai(l_, "granum"); tn = tb.con_loai(l_, "tam_noi")
        c2a.append({"bo_phan": l_.ten, "dat": trong_khoang(len(gr), *SO_GRANUM), "gia_tri": len(gr), "nguong": "%s grana (tv.so_granum)" % (SO_GRANUM,)})
        for g_ in gr:
            th = tb.con_loai(g_, "thylakoid"); e = g_.kich_thuoc()
            c2b.append({"bo_phan": g_.ten, "dat": trong_khoang(e[0], 0.3, 0.6) and trong_khoang(e[1], 0.3, 0.6), "gia_tri": [round(e[0], 3), round(e[1], 3)], "nguong": "duong kinh [0,3; 0,6]"})
            c2c.append({"bo_phan": g_.ten, "dat": trong_khoang(len(th), *SO_THYL_VE), "gia_tri": len(th), "nguong": "%s thylakoid (tv.so_thylakoid_ve)" % (SO_THYL_VE,)})
            c2d.append({"bo_phan": g_.ten, "dat": trong_khoang(e[2], *THYL_VE["granum_cao"]), "gia_tri": round(e[2], 4), "nguong": "cao granum %s (tv.thylakoid_ve)" % (THYL_VE["granum_cao"],)})
            if th:
                V = np.concatenate([t.luoi.V for t in th]); cg, U = truc_pca(V); ax = U[:, 0]
                if np.ptp(V @ ax) < np.ptp(V @ U[:, 2]) * 0.2:
                    pass
                th = sorted(th, key=lambda t: float(t.luoi.V.mean(0) @ ax))
                ok_dep = True; lech = 0.0; goc = 0.0
                for t in th:
                    e_t = t.kich_thuoc(); _, Ut = truc_pca(t.luoi.V)
                    ok_dep &= e_t[2] / e_t[1] <= TV_PHANG
                    goc = max(goc, float(np.degrees(np.arccos(min(1.0, abs(float(Ut[:, 0] @ ax)))))))
                    cc = t.luoi.V.mean(0) - cg; lech = max(lech, float(np.linalg.norm(cc - (cc @ ax) * ax)))
                    c2d.append({"bo_phan": t.ten, "dat": trong_khoang(e_t[2], *THYL_VE["day"]), "gia_tri": round(e_t[2], 4), "nguong": "thylakoid day %s (tv.thylakoid_ve)" % (THYL_VE["day"],)})
                c2a.append({"bo_phan": g_.ten, "dat": ok_dep and goc <= TV_GOC_THANG and lech <= TV_LECH_CHONG, "gia_tri": {"det": ok_dep, "goc_truc_max": round(goc, 2), "lech_tam_max": round(lech, 4)}, "nguong": "chong thylakoid det, song song, dong truc"})
                for a, b in zip(th[:-1], th[1:]):
                    d, cx = kc(a, b, 1.0)
                    c2d.append({"bo_phan": [a.ten, b.ten], "dat": trong_khoang(d, THYL_VE["khe"][0] - 1e-6, THYL_VE["khe"][1] + 1e-6), "gia_tri": round(d, 4), "nguong": "khe thylakoid %s (tv.thylakoid_ve)" % (THYL_VE["khe"],)})
        cham = {}
        for t in tn:
            ds = []
            for g_ in gr:
                if hop_cach(t, g_) > TV_TIEP_XUC: continue
                S = dien_tich_cham(t, g_.luoi)
                if S > 0: ds.append((g_, S))
            cham[t.ten] = ds
        c5.append({"bo_phan": l_.ten, "dat": trong_khoang(len(tn), max(0, len(gr) - 1), SO_TAM_NOI_MAX) and len(tn) >= len(gr) - 1, "gia_tri": {"tam_noi": len(tn), "granum": len(gr)}, "nguong": "[so granum − 1; 40]"})
        par = {g_.ten: g_.ten for g_ in gr}
        def f_(x):
            while par[x] != x: par[x] = par[par[x]]; x = par[x]
            return x
        for t in tn:
            ds = cham[t.ten]
            for (ga, _), (gb, _) in zip(ds[:-1], ds[1:]):
                par[f_(ga.ten)] = f_(gb.ten)
            if not ds:
                c3a.append({"bo_phan": t.ten, "dat": False, "gia_tri": "khong cham granum nao", "nguong": "tam noi cham granum"})
        nthanh = len(set(f_(g_.ten) for g_ in gr)) if gr else 0
        c3a.append({"bo_phan": l_.ten, "dat": nthanh <= 1, "gia_tri": nthanh, "nguong": "he thylakoid lien 1 khoi (so thanh phan do thi granum–tam noi)"})
        for g_ in gr:
            Smax = max([S for t in tn for (gg, S) in cham[t.ten] if gg is g_] or [0.0])
            c3b.append({"bo_phan": g_.ten, "dat": Smax >= VUNG_TAM, "gia_tri": round(Smax, 5), "nguong": "co tam noi cham, vung noi >= %g µm² (tv.vung_noi_tam_granum)" % VUNG_TAM, "S_cham": Smax})
        cent = {g_.ten: g_.luoi.V.mean(0) for g_ in gr}
        for t in tn:
            ds = cham[t.ten]
            if len(ds) != 2:
                c3c.append({"bo_phan": t.ten, "dat": False, "gia_tri": len(ds), "nguong": "noi dung 2 granum"}); continue
            (ga, _), (gb, _) = ds; seg = cent[gb.ten] - cent[ga.ten]; seg /= np.linalg.norm(seg)
            e_t = t.kich_thuoc(); _, Ut = truc_pca(t.luoi.V)
            goc = float(np.degrees(np.arccos(min(1.0, abs(float(Ut[:, 2] @ seg))))))
            gan = sorted(gr, key=lambda g_: np.linalg.norm(cent[g_.ten] - cent[ga.ten]))[1:1 + TV_KE_NHAU]
            ke = gb in gan
            c3c.append({"bo_phan": t.ten, "dat": e_t[2] / e_t[1] <= TV_PHANG and goc <= TV_GOC_THANG and ke, "gia_tri": {"phang": round(e_t[2] / e_t[1], 3), "goc": round(goc, 2), "ke_nhau": ke}, "nguong": "phang, thang (goc <= %g°), noi hai granum ke nhau (trong %d gan nhat)" % (TV_GOC_THANG, TV_KE_NHAU)})
    add("10.1a", tong_hop(ca1, True), "Mỗi lục lạp có đúng một màng trong (bộ phận con) nằm trong, khe hai màng 0,05 ± TV_TOL_KHE_LAP µm (min chính xác, max theo đỉnh). DNA và ribosome chất nền không vẽ → không kiểm được.", ca1, ly_do="DNA/ribosome chất nền chưa vẽ (lời luật).")
    add("10.1b", tong_hop(cb1, False), "Màng trong không gấp nếp: khối màng trong lồi — V / bao lồi ≥ TV_NG_LOI VÀ độ lõm sâu nhất (đỉnh → mặt bao lồi) ≤ TV_NEP_SAU.", cb1)
    add("10.2a", tong_hop(c2a, False), "10–20 grana mỗi lục lạp; mỗi granum là chồng thylakoid dẹt (cạnh nhỏ/giữa ≤ TV_PHANG), song song (≤ TV_GOC_THANG°) và đồng trục (lệch tâm ≤ TV_LECH_CHONG).", c2a)
    add("10.3a", tong_hop(c3a, False), "Trong mỗi lục lạp, đồ thị granum–tấm nối (cạnh = tấm nối chạm hai granum) liền một thành phần; mọi tấm nối chạm ít nhất một granum.", c3a)
    add("10.3c", tong_hop(c3c, False), "Mỗi tấm nối chạm đúng 2 granum, phẳng (cạnh nhỏ/giữa ≤ TV_PHANG), trục dài lệch đoạn nối hai tâm ≤ TV_GOC_THANG°, granum kia trong TV_KE_NHAU granum gần nhất.", c3c)
    add("10.2b", tong_hop(c2b, False), "Đường kính granum (hai cạnh lớn OBB) ∈ [0,3; 0,6] µm.", c2b)
    add("10.2c", tong_hop(c2c, False), "Mỗi granum có 3–8 thylakoid con.", c2c)
    add("10.2d", tong_hop(c2d, False), "Thylakoid dày ∈ [0,04; 0,06], khe chính xác giữa hai thylakoid liên tiếp ∈ [0,01; 0,02], granum cao (cạnh nhỏ OBB) ∈ [0,14; 0,62] µm.", c2d)
    add("10.5", tong_hop(c5, False), "Số tấm nối mỗi lục lạp ∈ [số granum − 1; 40].", c5)
    add("10.3b", tong_hop(c3b, False), "Mỗi granum có ít nhất một tấm nối chạm vào, diện tích chạm ≥ 0,005 µm² (Phụ lục 1b).", c3b)

    _t("nghi_sai_loai")
    nghi = kiem_loai(tb)

    thu_tu = [l["id"] for l in LU["luat"]]
    by = {r["id"]: r for r in KQ}
    assert set(by) == set(thu_tu), (set(thu_tu) - set(by), set(by) - set(thu_tu))
    KQ = [by[i] for i in thu_tu]
    for r in KQ:
        if lop not in (luat[r["id"]].get("lop") or []):
            r["neu_ap_dung"] = r["trang_thai"]; r["trang_thai"] = "KHONG_AP_DUNG"
            r["ly_do"] = "Luật dành cho %s; mô hình lớp %s. " % ("/".join(luat[r["id"]]["lop"]), lop) + (r.get("ly_do") or "")
        mat_ = [c for c in r["ca_truot"] if str(c.get("gia_tri", "")).startswith("MAT")]
        r["so_ca_truot_mat"] = len(mat_); r["so_ca_truot_hinh"] = r["so_ca_truot"] - len(mat_)
        thieu_ = [lo_ for lo_ in PHU_THUOC_T2.get(r["id"], []) if not tb.tl.get(lo_)]
        r["thieu_bo_phan_tang2"] = thieu_
        if r["trang_thai"] == "TRUOT":
            r["loai_truot"] = ("THIEU_BO_PHAN_TANG2" if thieu_ else ("HINH_HOC" if r["so_ca_truot_hinh"] > 0 else "CHI_DO_MAT_BO_PHAN"))
    dem = collections.Counter(r["trang_thai"] for r in KQ)
    sinh_hoc_truot = [r["id"] for r in KQ if r["trang_thai"] == "TRUOT" and not r.get("kiem_rieng_day_hoc")]
    truot_tach = dict(collections.Counter(r.get("loai_truot") for r in KQ if r["trang_thai"] == "TRUOT"))
    out = {"te_bao": os.path.basename(os.path.normpath(thu_muc)), "vat": "thuc_vat", "lop": lop, "so_luat": len(KQ),
           "dem_trang_thai": dict(dem), "truot_sinh_hoc": sinh_hoc_truot, "truot_tach": truot_tach,
           "nghi_sai_loai": {"ghi": "khong tinh vao 59 luat; bo phan co hinh/co ngoai dai cua loai khai bao", "so": len(nghi), "ds": nghi},
           "truot_day_hoc_bao_rieng": [r["id"] for r in KQ if r["trang_thai"] == "TRUOT" and r.get("kiem_rieng_day_hoc")],
           "bo_phan_mat": tb.mat, "so_bo_phan": len(tb.bp), "cha_suy_ra": tb.cha_suy,
           "loi_luoi": {b.ten: {"da_xu_ly": b.loi, "CHAN": b.loi_chan} for b in tb.bp.values() if b.loi or b.loi_chan},
           "tham_so_cho_duyet": CHO_DUYET_TV, "thoi_gian_s": round(time.time() - t0, 1), "luat": KQ}
    global _NGU_CANH
    _NGU_CANH = dict(locals())
    json.dump(out, open(os.path.join(thu_muc, "kiem2.json"), "w"), ensure_ascii=False, indent=1, default=_json_default)
    if not im_lang:
        print("== %s (thuc vat, lop %s, %d bo phan, mat %d) %.1f s" % (out["te_bao"], lop, len(tb.bp), len(tb.mat), out["thoi_gian_s"]))
        print("   dem:", out["dem_trang_thai"], "/ %d" % len(KQ), " truot sinh hoc:", sinh_hoc_truot, " day hoc:", out["truot_day_hoc_bao_rieng"], " tach:", truot_tach)
        print("   nghi_sai_loai: %d" % len(nghi), nghi[:3])
        for r in KQ:
            if r["trang_thai"] == "TRUOT":
                c = r["ca_truot"][0]; print("   %-6s TRUOT %d ca; vd %s: %s (%s)" % (r["id"], r["so_ca_truot"], c.get("bo_phan"), str(c.get("gia_tri"))[:90], c.get("nguong")))
    return out


if os.environ.get("KIEM_TANG_TOC", "1") != "0":
    import sys as _sys, tang_toc as _tang_toc
    _tang_toc.cai_tv(_sys.modules[__name__])
