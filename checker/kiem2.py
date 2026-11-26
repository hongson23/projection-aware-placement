# -*- coding: utf-8 -*-
"""
kiem2.py — BỘ KIỂM ĐỘC LẬP (bản 2) cho mô hình tế bào động vật dạy học.

Viết độc lập: chỉ dựa trên LỜI luật (luat_v41.json, trường "phat_bieu") và các mesh bao
ngoài bào quan (mesh.npz + loai.json). Không đọc mã bộ dựng hay bộ kiểm cũ.

Cách dùng:  python kiem2.py <thu_muc_te_bao> [--lop Sinh10|KHTN6] [--luat duong_dan]
Ghi:        <thu_muc_te_bao>/kiem2.json, in tóm tắt.

Nguyên tắc hình học (bảo thủ):
  * Chồng lấn / chứa: phép bool chính xác của manifold3d (Mesh64, số thực kép).
      - chồng lấn  ⇔ thể tích giao  > EPS_VOL
      - A nằm trong B ⇔ thể tích (A − B) ≤ EPS_VOL
  * Khoảng cách (khe) giữa hai mặt: CHÍNH XÁC, không lấy mẫu:
      với hai lưới tam giác không cắt nhau, khoảng cách nhỏ nhất đạt ở cặp (đỉnh, tam giác)
      hoặc (cạnh, cạnh). Ta tính đỉnh→tam giác bằng igl (chính xác) theo cả hai chiều, và
      cạnh↔cạnh bằng công thức đoạn–đoạn đóng (Ericson 2005) cho mọi cặp cạnh có thể cho
      khoảng cách < ngưỡng (lọc bằng KD-tree với bán kính đã cộng nửa độ dài cạnh lớn nhất,
      nên lọc không bỏ sót). Sai số = sai số làm tròn số thực (~1e-12 µm).
      Nếu khoảng cách ≥ ngưỡng R thì chỉ báo cận dưới "≥ R" (không cần số chính xác).
  * Nhân rỗng / có lỗ: "trong nhân" = khối ĐẶC của nhân (khoang nhân kể cả phần rỗng bên
      trong vỏ). Cách dựng khối đặc:
        - nếu thành phần mặt lớn nhất của nhân có giống 0 (mặt cầu tôpô) → khối đặc chính
          xác = khối bao bởi riêng mặt ngoài đó (bỏ mọi khoang trong);
        - nếu mặt ngoài có giống > 0 (vỏ thủng lỗ) → ĐÓNG HÌNH THÁI trên lưới voxel:
          giãn bán kính r → lấp khoang kín → co lại r (r nhỏ nhất làm kín khoang, + 1 bước dư);
          rồi marching cubes → manifold. Giới hạn: lõm hẹp hơn 2r ở mặt ngoài cũng bị lấp;
          sai số mặt ≈ một cạnh voxel (H_VOXEL). Bao lồi bị loại vì lấp mọi chỗ lõm của nhân
          hình hạt đậu.
  * Lưới lỗi (không kín, không đa tạp, mặt suy biến): báo rõ; dọn thành phần diện tích 0;
      nếu vẫn không dựng được manifold thì dùng phương pháp dự phòng (số vòng quấn trên lưới
      điểm trong) và ghi "phuong_phap": "du_phong_mau_diem".
"""
import sys, os, json, time, argparse, collections
import numpy as np
import trimesh
import igl
import manifold3d as mf
from scipy.spatial import cKDTree, ConvexHull, Delaunay
from scipy import ndimage

EPS_VOL = 1e-6
H_VOXEL = 0.05
R_DONG_MAX = 1.0
R_DONG_MAX2 = 3.0
V_KHOANG_MIN = 0.05
TU_CAT_TOL = 1e-6
H_DAY = 0.02
NG_12B_TRON = 0.8
NG_46D_CAU = 0.85
DAY_ER_UM = 2.0
DUNG_SAI_44B = 0.05
DUNG_SAI_KHOANG = 0.0
TIEP_XUC_UM = 0.01
NG_45D_TI_LE = 0.5
MD_KHOANG_MAX = 0.3
MD_TI_LE_TIA = 0.5
NG_48A_GOC = 75.0
NG_DAC_VO = 0.6
NG_DET_43B = 1.5
NG_CONG_43B = 0.05
MIEN_MC = [("luoi_noi_chat", "nhan"), ("nhan_con", "nhan"), ("ribosome", "luoi_noi_chat")]
DINH_NGHIA_MEM = {
    "4.5b": "cạnh nhân ⇔ khe Golgi–nhân < khe Golgi–màng; ôm trung thể ⇔ tâm trung thể trong bao lồi Golgi",
    "4.6b": "u = d(nhân)/(d(nhân)+d(màng)) tại tâm lysosome; tập trung ⇔ trung vị u < 0,5; rải ra ⇔ max u ≥ 0,5",
    "4.8b": "gần nhân ⇔ khe trung thể–nhân < khe trung thể–màng",
}
CHO_DUYET = {"NG_12B_TRON": NG_12B_TRON, "NG_46D_CAU": NG_46D_CAU, "DAY_ER_UM": DAY_ER_UM, "DUNG_SAI_44B": DUNG_SAI_44B,
             "DUNG_SAI_KHOANG": DUNG_SAI_KHOANG, "TIEP_XUC_UM": TIEP_XUC_UM, "NG_45D_TI_LE": NG_45D_TI_LE,
             "MD_KHOANG_MAX": MD_KHOANG_MAX, "MD_TI_LE_TIA": MD_TI_LE_TIA, "NG_48A_GOC": NG_48A_GOC,
             "NG_DAC_VO": NG_DAC_VO, "R_DONG_MAX2": R_DONG_MAX2, "NG_DET_43B": NG_DET_43B, "NG_CONG_43B": NG_CONG_43B,
             "4.3b_rong": "rong = do day cuc bo (2 x ban kinh cau noi tiep); dung sai do (can sai so) nghieng ve DAT",
             "M.h": "chong lan cung loai: giao > EPS_VOL (1e-6 µm³) la truot; cham duoc phep",
             "MIEN_MC": MIEN_MC, "DINH_NGHIA_MEM": DINH_NGHIA_MEM}

HERE = os.path.dirname(os.path.abspath(__file__))
import luat_v12 as _LV
_P = _LV.P
EPS_VOL = float(_P("chung.eps_giao_um3"))
NG_12B_TRON = float(_P("kiem.dv.NG_12B_TRON")); NG_46D_CAU = float(_P("kiem.dv.NG_46D_CAU"))
DAY_ER_UM = float(_P("dv.D_er")); DUNG_SAI_44B = float(_P("kiem.dv.DUNG_SAI_44B"))
DUNG_SAI_KHOANG = float(_P("chung.dung_sai_khoang")); TIEP_XUC_UM = float(_P("chung.tiep_xuc_toi_da_um"))
MD_KHOANG_MAX = float(_P("kiem.dv.MD_KHOANG_MAX")); MD_TI_LE_TIA = float(_P("kiem.dv.MD_TI_LE_TIA"))
NG_48A_GOC = float(_P("dv.goc_vuong_toi_thieu")); NG_CONG_43B = float(_P("kiem.dv.NG_CONG_43B")); NG_DET_43B = float(_P("kiem.dv.NG_DET_43B"))
KHE_ROI = float(_P("chung.khe_roi")); KHE_CUNG_LOAI_DV = dict(_P("chung.khe_cung_loai")["dv"])
T_NHAN_MANG = float(_P("dv.T_nhan_mang")); T_MANG = float(_P("dv.T_mang")); T_NHANCON = float(_P("dv.T_nhancon"))
SO_TI_THE = tuple(_P("dv.so_ti_the")); SO_LYSOSOME = tuple(_P("dv.so_lysosome")); TI_THE_LOP6 = tuple(_P("dv.ti_the_lop6")["so_luong"])
LIEN_KHOI_MIEN = set(_P("dv.lien_khoi_mien"))
MIEN_MC = [tuple(x) for x in _P("dv.mien_cap")] + [("nhan_con", "nhan")]
LOAI_LOP10 = {"nhan_con", "luoi_noi_chat", "golgi", "trung_the", "lysosome"}


def kc_doan_doan(P1, Q1, P2, Q2):
    """Khoảng cách CHÍNH XÁC giữa các cặp đoạn thẳng [P1,Q1]–[P2,Q2] (vector hoá, Ericson 2005)."""
    d1 = Q1 - P1; d2 = Q2 - P2; r = P1 - P2
    a = np.einsum('ij,ij->i', d1, d1); e = np.einsum('ij,ij->i', d2, d2)
    f = np.einsum('ij,ij->i', d2, r); c = np.einsum('ij,ij->i', d1, r)
    b = np.einsum('ij,ij->i', d1, d2)
    tiny = 1e-20
    a_ = np.where(a > tiny, a, 1.0); e_ = np.where(e > tiny, e, 1.0)
    den = a * e - b * b
    s = np.where(den > tiny, np.clip((b * f - c * e) / np.where(den > tiny, den, 1.0), 0, 1), 0.0)
    t = (b * s + f) / e_
    s = np.where(t < 0, np.clip(-c / a_, 0, 1), np.where(t > 1, np.clip((b - c) / a_, 0, 1), s))
    t = np.clip(t, 0, 1)
    s = np.where(a > tiny, s, 0.0)
    t = np.where(e > tiny, t, np.clip(f / e_, 0, 1) * 0.0)
    t = np.where((a <= tiny) & (e > tiny), np.clip(f / e_, 0, 1), t)
    C1 = P1 + d1 * s[:, None]; C2 = P2 + d2 * t[:, None]
    return np.linalg.norm(C1 - C2, axis=1)


def tao_manifold(V, F):
    try:
        m = mf.Manifold(mf.Mesh64(np.ascontiguousarray(V, dtype=np.float64),
                                  np.ascontiguousarray(F, dtype=np.uint64)))
        st = m.status()
        if st != mf.Error.NoError or m.is_empty():
            return None, str(st)
        return m, "NoError"
    except Exception as ex:
        return None, "ngoai_le:" + str(ex)


class Luoi:
    """Một lưới tam giác (đã dọn) + các cấu trúc tra cứu dựng lười."""
    def __init__(self, V, F):
        self.V = np.ascontiguousarray(V, dtype=np.float64)
        self.F = np.ascontiguousarray(F, dtype=np.int64)
        E = np.vstack([self.F[:, [0, 1]], self.F[:, [1, 2]], self.F[:, [2, 0]]])
        E.sort(axis=1)
        self.E = np.unique(E, axis=0)
        L = np.linalg.norm(self.V[self.E[:, 0]] - self.V[self.E[:, 1]], axis=1)
        self.Lmax = float(L.max()) if len(L) else 0.0
        self.lo = self.V.min(0); self.hi = self.V.max(0)
        self._aabb = None; self._kdt = None; self._kde = None

    def aabb(self):
        if self._aabb is None:
            self._aabb = igl.AABB(); self._aabb.init(self.V, self.F)
        return self._aabb

    def kdt(self):
        if self._kdt is None:
            self._kdt = cKDTree(self.V)
        return self._kdt

    def kc2_diem(self, P):
        """Bình phương khoảng cách chính xác từ điểm tới mặt (tam giác)."""
        P = np.ascontiguousarray(P, dtype=np.float64)
        if len(P) == 0:
            return np.zeros(0)
        d, _, _ = self.aabb().squared_distance(self.V, self.F, P)
        return d

    def quan(self, P):
        """Số vòng quấn (igl.fast_winding_number, Barnes–Hut; sai lệch chỉ đáng kể sát mặt) — >0,5 là bên trong."""
        return igl.fast_winding_number(self.V, self.F, np.ascontiguousarray(P, dtype=np.float64))


def kc_vf(A, B):
    """Khoảng cách đỉnh→tam giác theo hai chiều (cận trên của khoảng cách mặt–mặt; sai lệch
    ≤ nửa cạnh dài nhất). Dùng cho các luật vị trí tương đối (không có ngưỡng tuyệt đối)."""
    return float(min(np.sqrt(B.kc2_diem(A.V)).min(), np.sqrt(A.kc2_diem(B.V)).min()))


def kc_mat_mat(A, B, R):
    """Khoảng cách nhỏ nhất CHÍNH XÁC giữa mặt A và mặt B nếu < R; ngược lại trả (R, False)
    nghĩa là đã chứng minh ≥ R.  Trả (giá_trị, chinh_xac:bool)."""
    dA = np.sqrt(B.kc2_diem(A.V))
    U = float(dA.min())
    rr = R + A.Lmax + B.Lmax
    idxB = A.kdt().query_ball_point(B.V, r=rr, return_length=True) if len(B.V) < len(A.V) else None
    if idxB is not None:
        candB = np.nonzero(idxB > 0)[0]
    else:
        lst = B.kdt().query_ball_point(A.V, r=rr)
        candB = np.unique(np.concatenate([np.asarray(x, dtype=np.int64) for x in lst])) if len(lst) else np.zeros(0, np.int64)
    dB = np.full(len(B.V), np.inf)
    if len(candB):
        dB[candB] = np.sqrt(A.kc2_diem(B.V[candB]))
        U = min(U, float(dB[candB].min()))
    lim = min(U, R)
    eA = A.E[np.minimum(dA[A.E[:, 0]], dA[A.E[:, 1]]) < lim + A.Lmax / 2 + 1e-12]
    eB = B.E[np.minimum(dB[B.E[:, 0]], dB[B.E[:, 1]]) < lim + B.Lmax / 2 + 1e-12]
    best = U
    if len(eA) and len(eB):
        mA = 0.5 * (A.V[eA[:, 0]] + A.V[eA[:, 1]]); mB = 0.5 * (B.V[eB[:, 0]] + B.V[eB[:, 1]])
        tA = cKDTree(mA); tB = cKDTree(mB)
        pairs = tA.query_ball_tree(tB, r=lim + (A.Lmax + B.Lmax) / 2 + 1e-12)
        ii = np.repeat(np.arange(len(pairs)), [len(p) for p in pairs])
        if len(ii):
            jj = np.concatenate([np.asarray(p, dtype=np.int64) for p in pairs if len(p)])
            for s in range(0, len(ii), 2_000_000):
                a_ = eA[ii[s:s + 2_000_000]]; b_ = eB[jj[s:s + 2_000_000]]
                dd = kc_doan_doan(A.V[a_[:, 0]], A.V[a_[:, 1]], B.V[b_[:, 0]], B.V[b_[:, 1]])
                best = min(best, float(dd.min()))
    if best < R:
        return best, True
    return R, False


def tu_cat(V, F, toi_da=3000):
    """Số mặt TỰ CẮT đã xác minh. pymeshlab chỉ dùng để sinh ứng viên (nó có dung sai, báo nhầm trên lưới
    sinh bởi bool); mỗi ứng viên được kiểm lại CHÍNH XÁC (igl.simplex_simplex_squared_distance):
      - hai tam giác không chung đỉnh: cắt ⇔ khoảng cách = 0;
      - chung một đỉnh v: cắt ⇔ cạnh đối diện v của tam giác này chạm tam giác kia (hoặc ngược lại);
      - chung một cạnh: chồng gấp ⇔ đồng phẳng ngược hướng và hai đỉnh còn lại cùng phía cạnh chung.
    Trả -1 nếu không chạy được."""
    try:
        import pymeshlab as ml
        V = np.asarray(V, float); F = np.asarray(F, np.int64)
        ms = ml.MeshSet(); ms.add_mesh(ml.Mesh(V, F.astype(np.int32)))
        ms.compute_selection_by_self_intersections_per_face()
        ung = np.nonzero(ms.current_mesh().face_selection_array())[0][:toi_da]
    except Exception:
        return -1
    if len(ung) == 0:
        return 0
    T = V[F]
    blo = T.min(1) - 1e-9; bhi = T.max(1) + 1e-9
    xac = set()
    for i in ung:
        nb = np.nonzero(np.all(blo <= bhi[i], 1) & np.all(bhi >= blo[i], 1))[0]
        for j in nb:
            if j == i or (i, j) in xac: continue
            chung = set(F[i]) & set(F[j])
            if len(chung) == 0:
                d, _, _ = igl.simplex_simplex_squared_distance(T[i], T[j]); cat = d <= 1e-24
                if cat:
                    def hai_phia(P, Q):
                        nq = np.cross(Q[1] - Q[0], Q[2] - Q[0]); ln = np.linalg.norm(nq)
                        if ln < 1e-18: return False
                        sd = (P - Q[0]) @ (nq / ln)
                        return sd.max() > TU_CAT_TOL and sd.min() < -TU_CAT_TOL
                    cat = hai_phia(T[i], T[j]) and hai_phia(T[j], T[i])
            elif len(chung) == 1:
                v = list(chung)[0]
                ei = V[[x for x in F[i] if x != v]]; ej = V[[x for x in F[j] if x != v]]
                d1, _, _ = igl.simplex_simplex_squared_distance(ei, T[j]); d2, _, _ = igl.simplex_simplex_squared_distance(ej, T[i])
                cat = min(d1, d2) <= 1e-24
                if cat:
                    def hai_phia2(P, Q):
                        nq = np.cross(Q[1] - Q[0], Q[2] - Q[0]); ln = np.linalg.norm(nq)
                        if ln < 1e-18: return False
                        sd = (P - Q[0]) @ (nq / ln)
                        return sd.max() > TU_CAT_TOL and sd.min() < -TU_CAT_TOL
                    cat = hai_phia2(ej, T[i]) and hai_phia2(ei, T[j])
            elif len(chung) == 2:
                a, b = [V[x] for x in chung]
                pi = V[[x for x in F[i] if x not in chung][0]]; pj = V[[x for x in F[j] if x not in chung][0]]
                ni = np.cross(T[i][1] - T[i][0], T[i][2] - T[i][0]); nj = np.cross(T[j][1] - T[j][0], T[j][2] - T[j][0])
                ni /= np.linalg.norm(ni) + 1e-30; nj /= np.linalg.norm(nj) + 1e-30
                e = (b - a); e /= np.linalg.norm(e) + 1e-30
                si = np.cross(e, pi - a); sj = np.cross(e, pj - a)
                cat = (ni @ nj < -1 + 1e-9) and (si @ sj > 0)
            else:
                cat = True
            if cat:
                xac.add((i, j)); xac.add((j, i))
    return len({i for i, _ in xac})


class BoPhan:
    def __init__(self, ten, loai, V, F, kiem_long_nhau=True, kiem_tu_cat=True):
        self.ten = ten; self.loai = loai
        self.loi = []
        self.loi_chan = []
        self.phuong_phap = "bool_chinh_xac"
        tm = trimesh.Trimesh(np.asarray(V, float), np.asarray(F, np.int64), process=False)
        ar = tm.area_faces
        if (ar <= 1e-14).any():
            self.loi.append("co %d mat dien tich 0 (da bo)" % int((ar <= 1e-14).sum()))
            tm.update_faces(ar > 1e-14); tm.remove_unreferenced_vertices()
        comps = tm.split(only_watertight=False)
        giu = []
        for c in comps:
            if c.area <= 1e-9 or abs(c.volume) <= 1e-9:
                self.loi.append("thanh phan suy bien (%d dinh, dien tich %.2e) da bo" % (len(c.vertices), c.area))
            else:
                giu.append(c)
        if len(giu) != len(comps):
            tm = trimesh.util.concatenate(giu) if giu else tm
        self.kin = bool(tm.is_watertight); self.cung_chieu = bool(tm.is_winding_consistent)
        if not self.kin: self.loi_chan.append("luoi khong kin")
        if not self.cung_chieu: self.loi_chan.append("huong mat khong nhat quan")
        vol = float(tm.volume) if self.kin else float('nan')
        if self.kin and vol < 0:
            self.loi.append("huong mat nguoc (the tich am) — da lat")
            tm.invert()
        self.tm = tm
        self.comps = [c for c in tm.split(only_watertight=False)] if len(giu) else []
        nt = tu_cat(tm.vertices, tm.faces) if kiem_tu_cat else 0
        self.so_mat_tu_cat = nt
        if nt > 0: self.loi_chan.append("tu cat: %d mat" % nt)
        if nt < 0: self.loi.append("khong chay duoc phep thu tu cat")
        vols = [float(c.volume) for c in self.comps]
        self.vols_tp = vols
        if kiem_long_nhau:
            self._kiem_long_nhau()
        self.luoi = Luoi(tm.vertices, tm.faces)
        self.man, st = tao_manifold(tm.vertices, tm.faces)
        if self.man is None:
            self.loi_chan.append("khong dung duoc manifold (%s)" % st)
        if self.loi_chan:
            self.man = None
            self.phuong_phap = "du_phong_mau_diem"
        self.so_khoi = int(sum(1 for v in vols if v > 1e-9))
        self.the_tich = float(self.man.volume()) if self.man is not None else (vol if self.kin else float('nan'))
        self.lo, self.hi = self.luoi.lo, self.luoi.hi
        self.dac = None; self.dac_luoi = None; self.dac_cach = None; self.dac_loi = None; self.dac_mat = None
        self.co_khoang = False

    def _kiem_long_nhau(self):
        """Mỗi thành phần thể tích âm (mặt khoang) phải nằm HẲN trong một thành phần dương và không trùng
        nó; các thành phần dương không được giao nhau. Vi phạm → lỗi chặn (vd vỏ trùng lật ngược)."""
        pos = []; neg = []
        for c in self.comps:
            m, st = tao_manifold(*(lambda cc: (cc.vertices, cc.faces))(c if c.volume > 0 else _lat(c)))
            (pos if c.volume > 0 else neg).append((c, m))
        for i in range(len(pos)):
            for j in range(i + 1, len(pos)):
                a, b = pos[i], pos[j]
                if a[1] is None or b[1] is None: continue
                if np.any(a[0].bounds[0] > b[0].bounds[1]) or np.any(b[0].bounds[0] > a[0].bounds[1]): continue
                g = (a[1] ^ b[1]).volume()
                if g > EPS_VOL:
                    self.loi_chan.append("hai thanh phan duong giao nhau (%.2e µm³)" % g)
        for c, m in neg:
            ok = False
            for p, mp in pos:
                if m is None or mp is None: continue
                if (m - mp).volume() <= EPS_VOL and mp.volume() - m.volume() > EPS_VOL:
                    ok = True; break
            if not ok:
                self.loi_chan.append("thanh phan the tich am (%.3g µm³) khong nam han trong thanh phan duong nao (vo trung lat nguoc / khoang ho)" % c.volume)

    def kich_thuoc(self):
        if getattr(self, "_kt", None) is not None:
            return self._kt
        self._kt = self._kich_thuoc()
        return self._kt

    def _kich_thuoc(self):
        try:
            e = np.sort(self.tm.bounding_box_oriented.primitive.extents)[::-1]
        except Exception:
            v = self.tm.vertices - self.tm.vertices.mean(0)
            _, U = np.linalg.eigh(np.cov(v.T)); p = v @ U
            e = np.sort(p.max(0) - p.min(0))[::-1]
        return [float(x) for x in e]

    def do_day(self):
        """Độ dày cục bộ lớn nhất = 2 × bán kính cầu nội tiếp lớn nhất (trục giữa). Trả (độ dày, cận sai số).
        Cách đo: lưới thô bước h trong hộp bao → điểm trong (số vòng quấn) → khoảng cách chính xác tới mặt;
        mọi điểm thô có d ≥ d_max − h√3/2 (có thể chứa cực đại toàn cục) được tinh chỉnh bằng lưới 5³ bước
        h/4 phủ ô của nó; rồi đánh bóng Nelder–Mead. Cận sai số trên độ dày (Lipschitz) = 2·(h/4)·√3/2."""
        if getattr(self, "_day", None) is not None:
            return self._day
        h = min(H_DAY, max(1e-3, (self.hi - self.lo).min() / 8))
        xs = [np.arange(self.lo[i] + h / 2, self.hi[i], h) for i in range(3)]
        G = np.stack(np.meshgrid(*xs, indexing='ij'), -1).reshape(-1, 3)
        G = G[self.luoi.quan(G) > 0.5]
        if len(G) == 0:
            self._day = (0.0, h); return self._day
        d = np.sqrt(self.luoi.kc2_diem(G))
        C = G[d >= d.max() - h * np.sqrt(3) / 2]
        C = C[np.argsort(-np.sqrt(self.luoi.kc2_diem(C)))[:400]]
        o = (np.arange(5) - 2) * (h / 4)
        O = np.stack(np.meshgrid(o, o, o, indexing='ij'), -1).reshape(-1, 3)
        G2 = (C[:, None, :] + O[None]).reshape(-1, 3)
        G2 = G2[self.luoi.quan(G2) > 0.5]
        d2 = np.sqrt(self.luoi.kc2_diem(G2)) if len(G2) else np.array([d.max()])
        p0 = G2[int(d2.argmax())] if len(G2) else G[int(d.argmax())]
        r = float(max(d.max(), d2.max()))
        try:
            from scipy.optimize import minimize
            def f(x):
                x = x[None]
                return -float(np.sqrt(self.luoi.kc2_diem(x))[0]) if self.luoi.quan(x)[0] > 0.5 else 1.0
            res = minimize(f, p0, method="Nelder-Mead", options={"xatol": 1e-6, "fatol": 1e-8, "maxiter": 400})
            if -res.fun > r: r = float(-res.fun)
        except Exception:
            pass
        self._day = (2 * r, float(2 * (h / 4) * np.sqrt(3) / 2))
        return self._day


def dai_truc_giua(bp, h=0.02, k_muot=3):
    """Chiều dài theo TRỤC GIỮA (cho hạt cong): voxel hoá bước h (≤ 3 triệu điểm) → khung xương
    (skimage skeletonize) → đường đi dài nhất trên đồ thị 26-láng giềng (hai lượt Dijkstra) → làm trơn
    trung bình trượt (±k_muot điểm, giữ hai đầu) để bỏ bậc thang → độ dài đường gấp khúc + khoảng cách
    chính xác từ hai đầu khung tới mặt. Trả (độ dài, cận sai số ≈ 2h) hoặc (None, None)."""
    try:
        from skimage.morphology import skeletonize
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import dijkstra
    except Exception:
        return None, None
    ext = bp.hi - bp.lo
    h = max(h, float(np.prod(ext + 4 * h) / 3e6) ** (1 / 3))
    lo = bp.lo - 2 * h
    n = np.ceil((ext + 4 * h) / h).astype(int) + 1
    xs = [lo[i] + h * np.arange(n[i]) for i in range(3)]
    G = np.stack(np.meshgrid(*xs, indexing='ij'), -1).reshape(-1, 3)
    occ = (bp.luoi.quan(G) > 0.5).reshape(n)
    if occ.sum() < 8:
        return None, None
    sk = skeletonize(occ)
    idx = np.argwhere(sk)
    if len(idx) < 2:
        return None, None
    P = lo + idx * h
    kd = cKDTree(idx)
    pr = kd.query_pairs(r=np.sqrt(3) + 1e-6, output_type='ndarray')
    if len(pr) == 0:
        return None, None
    wgt = np.linalg.norm(idx[pr[:, 0]] - idx[pr[:, 1]], axis=1) * h
    A = coo_matrix((np.r_[wgt, wgt], (np.r_[pr[:, 0], pr[:, 1]], np.r_[pr[:, 1], pr[:, 0]])), shape=(len(idx), len(idx))).tocsr()
    d0 = dijkstra(A, indices=0); d0[~np.isfinite(d0)] = -1
    a = int(np.argmax(d0))
    da, pred = dijkstra(A, indices=a, return_predecessors=True); da[~np.isfinite(da)] = -1
    b = int(np.argmax(da))
    path = [b]
    while path[-1] != a and pred[path[-1]] >= 0:
        path.append(int(pred[path[-1]]))
    Q = P[path]
    if len(Q) > 2 * k_muot + 1:
        Qs = Q.copy()
        for i in range(1, len(Q) - 1):
            k = min(k_muot, i, len(Q) - 1 - i)
            Qs[i] = Q[i - k:i + k + 1].mean(0)
        Q = Qs
    L = float(np.linalg.norm(np.diff(Q, axis=0), axis=1).sum())
    m_ = min(6, len(Q) - 1); them = 0.0
    for o_, q_ in [(Q[0], Q[m_]), (Q[-1], Q[-1 - m_])]:
        u = o_ - q_; nu = np.linalg.norm(u)
        if nu < 1e-12: continue
        u = u / nu
        _, t, _ = bp.luoi.aabb().intersect_ray_first(bp.luoi.V, bp.luoi.F, o_[None], u[None])
        tt = float(t[0])
        them += tt if np.isfinite(tt) and tt < 10 else float(np.sqrt(bp.luoi.kc2_diem(o_[None]))[0])
    return L + them, 2 * h


def _lat(c):
    c2 = c.copy(); c2.invert(); return c2


def _mc_manifold(solid, lo, h):
    from skimage import measure
    pad = np.pad(solid.astype(np.float32), 1)
    vv, ff, _, _ = measure.marching_cubes(pad, 0.5, spacing=(h, h, h))
    vv = vv + lo - h
    mc = trimesh.Trimesh(vv, ff[:, ::-1], process=True)
    if mc.volume < 0:
        mc.invert()
    m, st = tao_manifold(mc.vertices, mc.faces)
    if m is None:
        m2 = trimesh.Trimesh(vv, ff, process=True); m, st = tao_manifold(m2.vertices, m2.faces); mc = m2
    return m, mc


def lam_day(bp):
    """Dựng khối ĐẶC (khoang) của bộ phận vỏ (màng, nhân).
    - mặt ngoài giống 0: khối bao bởi mặt ngoài (chính xác);
    - mặt ngoài giống > 0: tìm KHOANG KÍN bằng đóng hình thái (r = 0,1…R_DONG_MAX); nếu có khoang
      (lấp thêm > V_KHOANG_MIN) → khối đặc = đóng hình thái; nếu KHÔNG có khoang (vd nhân đặc có
      đường hầm) → khối đặc = chính vật liệu (chính xác), không lấp gì.
    Không bao giờ dùng bao lồi. Lỗi dựng → lỗi chặn trong loi_chan."""
    if bp.dac is not None or bp.dac_cach is not None:
        return
    if not bp.comps:
        bp.loi_chan.append("khong co thanh phan hop le de dung khoi dac"); bp.dac_cach = {"cach": "loi"}; return
    lon = max(bp.comps, key=lambda c: abs(c.volume))
    g = int(round((2 - lon.euler_number) / 2))
    if lon.is_watertight and g == 0:
        m, st = tao_manifold(lon.vertices, lon.faces)
        if m is not None:
            bp.dac = m; bp.dac_luoi = None; bp.dac_mat = Luoi(lon.vertices, lon.faces)
            bp.co_khoang = any(v < 0 for v in bp.vols_tp)
            bp.dac_cach = {"cach": "mat_ngoai_giong0_chinh_xac", "giong_mat_ngoai": 0,
                           "so_thanh_phan_mat": len(bp.comps), "co_khoang_kin": bp.co_khoang}
            return
        bp.loi_chan.append("mat ngoai giong 0 nhung khong dung duoc manifold (%s)" % st)
    def tim(h, r_tu, r_den):
        lo = bp.lo - r_den - 4 * h; hi = bp.hi + r_den + 4 * h
        n = np.ceil((hi - lo) / h).astype(int) + 1
        xs = [lo[i] + h * np.arange(n[i]) for i in range(3)]
        G = np.stack(np.meshgrid(*xs, indexing='ij'), -1).reshape(-1, 3)
        w = igl.fast_winding_number(bp.luoi.V, bp.luoi.F, G)
        occ = (w > 0.5).reshape(n); del G, w
        dout = ndimage.distance_transform_edt(~occ) * h
        nhat = []
        for r in np.arange(r_tu, r_den + 1e-9, 0.1):
            dil = dout <= r
            gain = (ndimage.binary_fill_holes(dil).sum() - dil.sum()) * h ** 3
            nhat.append((round(float(r), 2), round(float(gain), 3)))
            if gain > V_KHOANG_MIN:
                return float(r) + 0.1, lo, occ, dout, h, nhat
        return None, lo, occ, dout, h, nhat
    chon, lo, occ, dout, h, nhat = tim(H_VOXEL, 0.1, R_DONG_MAX)
    vat = occ.sum() * h ** 3
    try:
        v_bao = float(trimesh.convex.convex_hull(bp.luoi.V).volume)
    except Exception:
        v_bao = float('nan')
    do_dac = float(vat / v_bao) if v_bao > 0 else float('nan')
    canh_bao = None
    if chon is None and do_dac < NG_DAC_VO:
        chon, lo, occ, dout, h, nhat2 = tim(2 * H_VOXEL, R_DONG_MAX + 0.1, R_DONG_MAX2)
        nhat += nhat2
        if chon is not None:
            canh_bao = "khoang chi kin o r = %.1f µm (lo vo lon); lom hep hon 2r bi lap" % chon
    if chon is None:
        bp.dac = bp.man; bp.dac_luoi = None; bp.dac_mat = bp.luoi; bp.co_khoang = False
        bp.dac_cach = {"cach": "vat_lieu_giong%d_khong_khoang_kin_chinh_xac" % g, "giong_mat_ngoai": g,
                       "do_dac_vat_lieu_tren_bao_loi": round(do_dac, 3),
                       "nhat_ky_tim_r": nhat, "ghi": "khong co khoang kin → dung chinh vat lieu"}
        if do_dac < NG_DAC_VO:
            bp.loi_chan.append("vo rong (do dac %.2f < %.2f) giong %d nhung khong dong duoc khoang den r = %.1f µm — "
                               "KHONG XAC DINH duoc 'trong nhan'; vat trong khoang co the bi bo sot" % (do_dac, NG_DAC_VO, g, R_DONG_MAX2))
        if bp.man is None:
            bp.loi_chan.append("khong co manifold vat lieu cho khoi dac")
        return
    dil = dout <= chon
    fil = ndimage.binary_fill_holes(dil)
    ero = ndimage.distance_transform_edt(fil) * h > chon
    solid = ndimage.binary_fill_holes(ero | occ)
    m, mc = _mc_manifold(solid, lo, h)
    if m is None:
        bp.loi_chan.append("dong hinh thai khong dung duoc manifold — CHAN")
        bp.dac = bp.man; bp.dac_mat = bp.luoi; bp.dac_cach = {"cach": "loi_dong_hinh_thai"}; return
    bp.dac = m; bp.dac_luoi = Luoi(mc.vertices, mc.faces); bp.dac_mat = bp.dac_luoi; bp.co_khoang = True
    core = ndimage.distance_transform_edt(solid) * h > 2 * h
    bp.dac_loi, _ = _mc_manifold(core, lo, h)
    if canh_bao: bp.loi.append(canh_bao)
    bp.dac_cach = {"cach": "dong_hinh_thai_voxel", "giong_mat_ngoai": g, "h_voxel_um": h,
                   "ban_kinh_dong_um": round(chon, 2), "the_tich_vat_lieu_um3": round(float(vat), 3),
                   "do_dac_vat_lieu_tren_bao_loi": round(do_dac, 3),
                   "the_tich_khoi_dac_um3": round(float(m.volume()), 3), "nhat_ky_tim_r": nhat,
                   "gioi_han": "lom hep hon 2r o mat ngoai bi lap; sai so mat ~ 1 voxel"}


def hop_cach(a, b):
    """Khoảng cách giữa hai hộp bao trục (cận dưới của khoảng cách mặt/khối)."""
    d = np.maximum(0, np.maximum(a.lo - b.hi, b.lo - a.hi))
    return float(np.linalg.norm(d))


def khoi_cua(bp):
    """Khối dùng để xét chồng lấn: nhân/màng dùng khối đặc; còn lại dùng chính nó."""
    if bp.loai in ("nhan", "mang"):
        lam_day(bp); return bp.dac
    return bp.man


def mat_chua(bp):
    """Mặt dùng cho phép thử điểm trong/ngoài (số vòng quấn) của khối đặc."""
    if bp.loai in ("nhan", "mang"):
        lam_day(bp)
        if bp.dac_mat is not None: return bp.dac_mat
    return bp.luoi


def the_tich_giao(A, B):
    ma, mb = khoi_cua(A), khoi_cua(B)
    if ma is not None and mb is not None:
        return float((ma ^ mb).volume()), "bool_chinh_xac"
    return du_phong_the_tich(A, B, "giao"), "du_phong_mau_diem"


def the_tich_hieu(A, C):
    """Thể tích phần A nằm NGOÀI khối đặc của C."""
    ma, mc = khoi_cua(A), khoi_cua(C)
    if ma is not None and mc is not None:
        return float((ma - mc).volume()), "bool_chinh_xac"
    return du_phong_the_tich(A, C, "hieu"), "du_phong_mau_diem"


def du_phong_the_tich(A, B, kieu, h=0.02):
    """DỰ PHÒNG khi không có manifold: lưới điểm bước h trong hộp A; điểm trong A và trong/ngoài B theo
    số vòng quấn. Sai số ~ diện tích mặt × h/2."""
    lo, hi = A.lo.copy(), A.hi.copy()
    if kieu == "giao":
        lo = np.maximum(lo, B.lo); hi = np.minimum(hi, B.hi)
        if np.any(hi <= lo): return 0.0
    h = max(h, float(np.prod(hi - lo) / 2e6) ** (1 / 3))
    xs = [np.arange(lo[i] + h / 2, hi[i], h) for i in range(3)]
    G = np.stack(np.meshgrid(*xs, indexing='ij'), -1).reshape(-1, 3)
    inA = mat_chua(A).quan(G) > 0.5
    inB = mat_chua(B).quan(G[inA]) > 0.5
    return float((inB.sum() if kieu == "giao" else (~inB).sum()) * h ** 3)


def mat_kc(bp):
    """Các mặt dùng để đo khe: mặt vật liệu thật; thêm mặt khối đặc nếu khối đặc là xấp xỉ voxel."""
    out = [bp.luoi]
    if bp.dac_luoi is not None:
        out.append(bp.dac_luoi)
    return out


def dinh_ngoai(A, C):
    """Số đỉnh của A nằm NGOÀI khối đặc C (số vòng quấn < 0,5 và cách mặt > 1e-9 µm)."""
    L = mat_chua(C)
    w = L.quan(A.luoi.V)
    d = np.sqrt(L.kc2_diem(A.luoi.V))
    return int(((w < 0.5) & (d > 1e-9)).sum())


def khe(A, B, R):
    """Khe CHÍNH XÁC giữa hai bộ phận (0 nếu chồng lấn).  Trả dict."""
    hb = hop_cach(A, B)
    if hb >= R:
        return {"kc_um": round(hb, 4), "chinh_xac": False, "ghi": "cach hop bao >= nguong (can duoi)",
                "giao_um3": 0.0}
    gv, pp = the_tich_giao(A, B)
    if gv > EPS_VOL:
        return {"kc_um": 0.0, "chinh_xac": True, "giao_um3": round(gv, 7), "ghi": "chong lan", "pp": pp}
    best, cx = R, False
    for LA in mat_kc(A):
        for LB in mat_kc(B):
            d, c = kc_mat_mat(LA, LB, R)
            if c and d < best:
                best, cx = d, True
    return {"kc_um": round(best, 5), "chinh_xac": cx, "giao_um3": round(gv, 8), "pp": pp,
            "ghi": "" if cx else "da chung minh >= nguong"}


def chua(A, C, t):
    """A nằm trong khối đặc của C, cách mặt C ít nhất t.
    Trong ⇔ thể tích (A − C) ≤ EPS_VOL VÀ (nếu hai mặt chạm/cắt) không đỉnh nào của A ra ngoài C."""
    hv, pp = the_tich_hieu(A, C)
    out = {"ngoai_um3": round(hv, 8), "pp": pp}
    if hv > EPS_VOL:
        out.update({"trong": False, "kc_um": 0.0, "chinh_xac": True})
        return out
    R = t if t > 0 else 1e-6
    best, cx = R, False
    for LC in mat_kc(C) + ([C.dac_mat] if C.dac_mat is not None and C.dac_mat is not C.luoi and C.dac_mat is not C.dac_luoi else []):
        d, c = kc_mat_mat(A.luoi, LC, R)
        if c and d < best:
            best, cx = d, True
    if cx and best < 1e-9:
        nd = dinh_ngoai(A, C)
        out["so_dinh_ngoai"] = nd
        if nd > 0:
            out.update({"trong": False, "kc_um": 0.0, "chinh_xac": True}); return out
    out["trong"] = True
    out.update({"kc_um": round(best, 5) if t > 0 else None, "chinh_xac": cx})
    return out


def ket_qua(luat, trang_thai, cach_hieu, ca=None, ly_do="", bat_dong="", them=None):
    d = {"id": luat["id"], "phat_bieu": luat["phat_bieu"], "luat_cung": luat.get("luat_cung"),
         "lop": luat.get("lop"), "trang_thai": trang_thai, "cach_hieu": cach_hieu,
         "ly_do": ly_do, "bat_dong_voi_truong_may": bat_dong}
    ca = ca or []
    hong = [c for c in ca if not c.get("dat", True)]
    d["so_ca"] = len(ca); d["so_ca_truot"] = len(hong)
    d["ca_truot"] = hong
    d["ca_dat_mau"] = [c for c in ca if c.get("dat", True)][:30]
    d["ca_tat_ca"] = list(ca)
    if them: d.update(them)
    return d


def tong_hop(ca, co_phan_khong_kiem, co_kiem=True):
    if any(not c.get("dat", True) for c in ca):
        return "TRUOT"
    if not co_kiem:
        return "KHONG_KIEM_DUOC"
    return "MOT_PHAN" if co_phan_khong_kiem else "DAT"


def trong_khoang(x, lo, hi):
    return (x >= lo) and (x <= hi)


class TeBao:
    def __init__(self, thu_muc):
        self.thu_muc = thu_muc
        L = json.load(open(os.path.join(thu_muc, "loai.json")))
        z = np.load(os.path.join(thu_muc, "mesh.npz"))
        ten = sorted(set(k.rsplit("__", 1)[0] for k in z.files))
        self.loai = dict(L.get("loai", {}))
        self.du_kien = dict(L.get("bo_phan_du_kien", self.loai))
        self.bp = {}
        for t in ten:
            lo = self.loai.get(t) or self.du_kien.get(t)
            if lo is None:
                lo = "khong_ro"
            self.bp[t] = BoPhan(t, lo, z[t + "__v"], z[t + "__f"])
        self.mat = sorted(set(self.du_kien) - set(self.bp))
        self.thua = sorted(set(self.bp) - set(self.du_kien))
        self.theo_loai = collections.defaultdict(list)
        for t, b in self.bp.items():
            self.theo_loai[b.loai].append(b)
        self.mat_theo_loai = collections.defaultdict(list)
        for t in self.mat:
            self.mat_theo_loai[self.du_kien[t]].append(t)

    def mot(self, loai):
        x = self.theo_loai.get(loai, [])
        return x[0] if len(x) == 1 else (x[0] if x else None)


def ca_mat(tb, loai, luat_id):
    return [{"bo_phan": t, "dat": False, "gia_tri": "MAT (co trong bo_phan_du_kien, khong co trong mesh.npz)",
             "nguong": "phai ton tai"} for t in tb.mat_theo_loai.get(loai, [])]


def kiem(thu_muc, lop=None, duong_luat=None, im_lang=False):
    t0 = time.time()
    duong_luat = duong_luat or os.path.join(HERE, "luat_v41.json")
    LU = json.load(open(duong_luat))
    luat = {l["id"]: l for l in LU["luat"]}
    tb = TeBao(thu_muc)
    if lop is None:
        lop = "Sinh10" if (set(tb.du_kien.values()) & LOAI_LOP10) else "KHTN6"
    mang = tb.mot("mang"); nhan = tb.mot("nhan")
    nhan_list = tb.theo_loai.get("nhan", [])
    if mang is not None: lam_day(mang)
    if nhan is not None: lam_day(nhan)
    bao_quan = [b for b in tb.bp.values() if b.loai != "mang"]

    mien = {frozenset(x) for x in MIEN_MC}
    S = KHE_ROI
    cap = {}
    for i in range(len(bao_quan)):
        A = bao_quan[i]
        for j in range(i + 1, len(bao_quan)):
            B = bao_quan[j]
            khac = A.loai != B.loai
            if not khac and not (A.loai == "ti_the"):
                continue
            if khac and frozenset([A.loai, B.loai]) in mien:
                continue
            if hop_cach(A, B) >= S:
                cap[(A.ten, B.ten)] = {"kc_um": round(hop_cach(A, B), 4), "chinh_xac": False, "giao_um3": 0.0,
                                       "ghi": "cach hop bao >= nguong (can duoi)"}
                continue
            cap[(A.ten, B.ten)] = khe(A, B, S)
    t_cap = time.time() - t0
    cl_cung = []
    for i in range(len(bao_quan)):
        A = bao_quan[i]
        for j in range(i + 1, len(bao_quan)):
            B = bao_quan[j]
            if A.loai != B.loai or hop_cach(A, B) > 0:
                continue
            k = cap.get((A.ten, B.ten))
            gv = k["giao_um3"] if k is not None else the_tich_giao(A, B)[0]
            if gv > EPS_VOL:
                cl_cung.append({"bo_phan": [A.ten, B.ten], "loai": A.loai, "giao_um3": round(float(gv), 6)})

    trong_mang = {}
    if mang is not None:
        for b in bao_quan:
            tmin = T_NHAN_MANG if b.loai == "nhan" else T_MANG
            trong_mang[b.ten] = chua(b, mang, tmin)

    KQ = []

    def L_(i): return luat[i]

    l = L_("1.1"); ca = []
    ca.append({"bo_phan": "mang", "dat": mang is not None, "gia_tri": "co" if mang else "khong co / mat", "nguong": "phai co"})
    ca.append({"bo_phan": "nhan", "dat": nhan is not None, "gia_tri": "co" if nhan else "khong co / mat", "nguong": "phai co"})
    if mang is not None and nhan is not None:
        vtb = float((mang.dac - nhan.dac).volume()) if (mang.dac is not None and nhan.dac is not None) else \
            float(mang.the_tich - nhan.the_tich)
        ca.append({"bo_phan": "chat_te_bao", "dat": vtb > EPS_VOL, "gia_tri": round(vtb, 3),
                   "nguong": "the tich (mang - nhan) > 0 µm³"})
        ca.append({"bo_phan": "nhan", "dat": trong_mang[nhan.ten]["trong"],
                   "gia_tri": "the tich nhan ngoai mang = %s µm³" % trong_mang[nhan.ten]["ngoai_um3"],
                   "nguong": "nhan trong mang (mang bao ngoai)"})
    KQ.append(ket_qua(l, tong_hop(ca, False),
                      "Có đúng một bộ phận màng và một nhân; chất tế bào = khối màng trừ khối nhân có thể tích > 0; nhân nằm trong khối màng (màng bao ngoài).",
                      ca))

    l = L_("1.2a"); ca = []
    if mang is not None:
        e = mang.kich_thuoc()
        ok = all(trong_khoang(x, 10, 30) for x in e) and mang.the_tich > 0 and min(e) > 0
        ca.append({"bo_phan": mang.ten, "dat": ok, "gia_tri": {"kich_thuoc_hop_dinh_huong_um": [round(x, 3) for x in e],
                   "the_tich_um3": round(mang.the_tich, 2)}, "nguong": "khoi 3 chieu (the tich > 0); moi kich thuoc trong [10, 30] µm"})
    else:
        ca.append({"bo_phan": "mang", "dat": False, "gia_tri": "khong co", "nguong": "phai co"})
    KQ.append(ket_qua(l, tong_hop(ca, False),
                      "Đường kính tế bào = ba cạnh hộp bao định hướng nhỏ nhất của màng; cả ba phải trong 10–30 µm; khối 3 chiều = thể tích > 0 và cạnh nhỏ nhất > 0.",
                      ca, bat_dong="Trường máy ban_kinh_um=10 là cách dựng (mô hình), không phải tiêu chí; không dùng."))

    l = L_("1.2b"); ca = []
    if mang is not None:
        e = mang.kich_thuoc()
        c = mang.luoi.V.mean(0) if mang.man is None else np.array(mang.tm.center_mass)
        r = np.linalg.norm(mang.luoi.V - c, axis=1)
        ti_truc = e[2] / e[0]; ti_r = float(r.min() / r.max())
        ok = ti_truc >= NG_12B_TRON and ti_r >= NG_12B_TRON
        ca.append({"bo_phan": mang.ten, "dat": ok, "gia_tri": {"canh_nho_tren_canh_lon": round(ti_truc, 3),
                   "r_min_tren_r_max": round(ti_r, 3)}, "nguong": ">= NG_12B_TRON = %g ca hai (cho duyet)" % NG_12B_TRON})
    else:
        ca.append({"bo_phan": "mang", "dat": False, "gia_tri": "khong co", "nguong": "phai co"})
    KQ.append(ket_qua(l, tong_hop(ca, False),
                      "'Khối tròn': tỉ số cạnh hộp định hướng nhỏ nhất/lớn nhất ≥ NG_12B_TRON và r_min/r_max từ tâm khối ≥ NG_12B_TRON (= %g, chờ duyệt; lời luật không cho số)." % NG_12B_TRON,
                      ca))

    l = L_("1.3"); ca = []
    ca.append({"bo_phan": "nhan", "dat": len(nhan_list) == 1 and not tb.mat_theo_loai.get("nhan"),
               "gia_tri": len(nhan_list), "nguong": "dung 1 bo phan nhan"})
    if nhan is not None:
        nk = len(nhan.dac.decompose()) if nhan.dac is not None else nhan.so_khoi
        ca.append({"bo_phan": nhan.ten, "dat": nk == 1, "gia_tri": nk, "nguong": "khoi dac cua nhan lien 1 khoi"})
        if mang is not None:
            ca.append({"bo_phan": nhan.ten, "dat": trong_mang[nhan.ten]["trong"],
                       "gia_tri": "the tich ngoai mang %s µm³" % trong_mang[nhan.ten]["ngoai_um3"], "nguong": "nam trong khoi te bao"})
            cm = np.array(mang.tm.center_mass); cn = np.array(trimesh.Trimesh(*_mesh(nhan.dac)).center_mass) if nhan.dac is not None else nhan.luoi.V.mean(0)
            if cn is not None:
                lech = float(np.linalg.norm(cn - cm) / (mang.kich_thuoc()[0] / 2))
                ca.append({"bo_phan": nhan.ten, "dat": True, "gia_tri": round(lech, 3),
                           "nguong": "thong tin: do lech tam / ban kinh (khong bat buoc o giua; khong kiem)"})
    KQ.append(ket_qua(l, tong_hop(ca, False),
                      "Đếm bộ phận nhân = 1 và khối đặc của nhân liền một khối; nhân nằm trong khối tế bào. 'Không bắt buộc ở giữa' → không đặt ràng buộc vị trí, chỉ ghi độ lệch tâm.",
                      ca, bat_dong="Trường máy lech_tam_ti_le_ban_kinh [0,11; 0,25] (mềm) không có trong lời luật — chỉ báo thông tin."))

    l = L_("1.4a"); ca = ca_mat(tb, "nhan", "1.4a")
    for b in nhan_list:
        if mang is None: ca.append({"bo_phan": b.ten, "dat": False, "gia_tri": "khong co mang", "nguong": ""}); continue
        r = trong_mang[b.ten]
        ca.append({"bo_phan": b.ten, "dat": r["trong"], "gia_tri": "ngoai mang %s µm³" % r["ngoai_um3"], "nguong": "the tich ngoai mang <= %g" % EPS_VOL})
    KQ.append(ket_qua(l, tong_hop(ca, False), "Khối nhân ⊂ khối bên trong màng (bool: nhân − màng rỗng).", ca))

    l = L_("1.4b"); ca = ca_mat(tb, "nhan", "1.4b")
    for b in nhan_list:
        if mang is None: continue
        r = trong_mang[b.ten]
        ok = r["trong"] and r["kc_um"] is not None and r["kc_um"] >= T_NHAN_MANG
        ca.append({"bo_phan": b.ten, "dat": ok, "gia_tri": r["kc_um"], "chinh_xac": r.get("chinh_xac"), "nguong": ">= 1,0 µm (khe mat nhan – mat mang, chinh xac)"})
    KQ.append(ket_qua(l, tong_hop(ca, False), "Nhân nằm trong màng và khoảng cách nhỏ nhất giữa mặt nhân và mặt màng ≥ 1 µm (đỉnh–tam giác + cạnh–cạnh, chính xác).", ca))

    l = L_("1.5"); ca = []
    if nhan is not None and mang is not None:
        vn = float(nhan.dac.volume()) if nhan.dac is not None else nhan.the_tich
        vt = float(mang.dac.volume()) if mang.dac is not None else mang.the_tich
        dn = (6 * vn / np.pi) ** (1 / 3); dt = (6 * vt / np.pi) ** (1 / 3)
        ti = dn / dt
        ca.append({"bo_phan": nhan.ten, "dat": trong_khoang(ti, 0.35, 0.6), "gia_tri": round(ti, 4),
                   "nguong": "[0,35; 0,6] (duong kinh cau tuong duong khoi dac nhan / te bao)",
                   "them": {"d_nhan_um": round(dn, 3), "d_te_bao_um": round(dt, 3),
                            "ti_canh_lon": round(nhan.kich_thuoc()[0] / mang.kich_thuoc()[0], 4)}})
        e = nhan.kich_thuoc()
        ca.append({"bo_phan": nhan.ten, "dat": all(trong_khoang(x, 5, 10) for x in e),
                   "gia_tri": [round(x, 3) for x in e], "nguong": "nhan rong: ca ba canh hop dinh huong trong [5, 10] µm"})
    else:
        ca.append({"bo_phan": "nhan/mang", "dat": False, "gia_tri": "thieu", "nguong": "phai co"})
    KQ.append(ket_qua(l, tong_hop(ca, False),
                      "Hai mệnh đề, cả hai phải đạt: (i) tỉ số đường kính = đường kính cầu cùng thể tích của khối đặc nhân / của khối tế bào ∈ [0,35; 0,6]; (ii) 'rộng 5–10 µm' = mọi cạnh hộp bao định hướng nhỏ nhất của nhân ∈ [5, 10] µm (đo trực tiếp, không suy từ tỉ số).",
                      ca, bat_dong="Ghi chú máy 'với tế bào 20 µm, dải dùng được 0,35–0,5' hẹp hơn lời luật (0,6) — theo lời luật."))

    l = L_("1.6"); ti = tb.theo_loai.get("ti_the", [])
    ca = []
    n_ti = len(ti)
    ca.append({"bo_phan": "ti_the", "dat": (n_ti == 0) or trong_khoang(n_ti, *TI_THE_LOP6), "gia_tri": n_ti, "nguong": "0 (tat) hoac %s (dv.ti_the_lop6)" % (TI_THE_LOP6,)})
    for b in ti:
        if nhan is None: break
        k = cap.get((b.ten, nhan.ten)) or cap.get((nhan.ten, b.ten))
        if k is not None:
            ca.append({"bo_phan": b.ten, "dat": bool(k["kc_um"] > TIEP_XUC_UM and k.get("giao_um3", 0) <= EPS_VOL),
                       "gia_tri": k["kc_um"], "nguong": "khong cham nhan (kc > 0,01 µm, khong giao)"})
    KQ.append(ket_qua(l, tong_hop(ca, False), "Lớp 6: ti thể tắt (0) hoặc 5–10 cái, không cái nào chạm nhân (kc > %g µm, không giao)." % TIEP_XUC_UM, ca))

    l = L_("4.2a"); nc = tb.theo_loai.get("nhan_con", [])
    ca = ca_mat(tb, "nhan_con", "4.2a")
    ca.append({"bo_phan": "nhan_con", "dat": trong_khoang(len(nc), 1, 3), "gia_tri": len(nc), "nguong": "mo hinh ve 1–3 (that 1–5)"})
    for b in nc:
        e = b.kich_thuoc()
        ca.append({"bo_phan": b.ten, "dat": all(trong_khoang(x, 0.5, 5) for x in e), "gia_tri": [round(x, 3) for x in e], "nguong": "moi canh hop dinh huong trong [0,5; 5] µm"})
        if nhan is not None:
            r = chua(b, nhan, 0.0)
            ca.append({"bo_phan": b.ten, "dat": r["trong"], "gia_tri": "ngoai khoi nhan %s µm³" % r["ngoai_um3"], "nguong": "nam trong khoi dac nhan"})
        else:
            ca.append({"bo_phan": b.ten, "dat": False, "gia_tri": "khong co nhan", "nguong": "nam trong nhan"})
    KQ.append(ket_qua(l, tong_hop(ca, True),
                      "Số nhân con 1–3 (con số của mô hình trong ngoặc); mỗi cái 'rộng 0,5–5 µm' = mọi cạnh hộp định hướng trong [0,5; 5]; nằm trong khối đặc của nhân (kể cả khoang của nhân rỗng). 'Không có màng bao' KHÔNG kiểm được từ mặt bao (là chi tiết bên trong/vẽ viền).",
                      ca, ly_do="Phần 'không có màng bao' không kiểm được từ mesh bao ngoài."))

    l = L_("4.2b"); ca = ca_mat(tb, "nhan_con", "4.2b")
    for b in nc:
        if nhan is None: ca.append({"bo_phan": b.ten, "dat": False, "gia_tri": "khong co nhan", "nguong": ""}); continue
        r = chua(b, nhan, T_NHANCON)
        ok = r["trong"] and r["kc_um"] is not None and r["kc_um"] >= T_NHANCON
        ca.append({"bo_phan": b.ten, "dat": ok, "gia_tri": r["kc_um"], "trong_nhan": r["trong"], "nguong": ">= 0,3 µm toi mat nhan (mat vat lieu va mat khoi dac)"})
    KQ.append(ket_qua(l, tong_hop(ca, False), "Nhân con nằm trong khối đặc nhân và khoảng cách chính xác tới mặt màng nhân (mặt lưới nhân thật; với nhân thủng lỗ thêm cả mặt khối đặc phủ qua lỗ, lấy nhỏ nhất) ≥ 0,3 µm.", ca))

    l = L_("4.3a"); ca = ca_mat(tb, "ti_the", "4.3a")
    for b in ti:
        ok = True; gt = {}
        if mang is not None:
            r = trong_mang[b.ten]; ok &= r["trong"]; gt["ngoai_mang_um3"] = r["ngoai_um3"]
        else:
            ok = False
        if nhan is not None:
            k = cap.get((b.ten, nhan.ten)) or cap.get((nhan.ten, b.ten))
            g = k.get("giao_um3", 0.0) if k else 0.0
            ok &= g <= EPS_VOL; gt["giao_nhan_um3"] = g
        ca.append({"bo_phan": b.ten, "dat": bool(ok), "gia_tri": gt, "nguong": "trong mang va khong giao khoi dac nhan"})
    KQ.append(ket_qua(l, tong_hop(ca, False), "Chất tế bào = trong màng, ngoài nhân (GO). Ti thể ⊂ khối màng và thể tích giao với khối đặc nhân (kể cả khoang nhân rỗng) ≤ ε.", ca))

    l = L_("4.3b"); ca = ca_mat(tb, "ti_the", "4.3b")
    tk = DUNG_SAI_KHOANG
    for b in ti:
        e = b.kich_thuoc(); w, sw = b.do_day()
        Lg, sL = dai_truc_giua(b)
        cong = Lg is not None and Lg > e[0] * (1 + NG_CONG_43B)
        dai = Lg if cong else e[0]
        ok_dai = trong_khoang(dai, 1.0 - tk, 2.0 + tk)
        ok_w_chac = trong_khoang(w, 0.5 - tk, 1.0 + tk)
        ok_w = ok_w_chac or (w < 0.5 - tk and w + sw >= 0.5 - tk)
        det = (e[1] / w) if w > 0 else float('inf')
        ca.append({"bo_phan": b.ten, "dat": bool(ok_dai and ok_w),
                   "gia_tri": {"dai": round(dai, 4), "dai_cach": "truc_giua (hat cong)" if cong else "obb",
                               "dai_obb": round(e[0], 4), "dai_truc_giua": round(Lg, 4) if Lg is not None else None,
                               "sai_so_dai_truc_giua": sL,
                               "rong_day_cuc_bo": round(w, 4), "sai_so_day": round(sw, 4),
                               "rong_truc_giua_obb": round(e[1], 4), "canh_nho_obb": round(e[2], 4),
                               "ti_obb_giua_tren_day": round(det, 3), "co_det": bool(det > NG_DET_43B),
                               "dat_trong_sai_so": bool(ok_w and not ok_w_chac)},
                   "nguong": "dai [1, 2] µm (OBB; truc giua neu cong > %g%%); rong = do day cuc bo [0,5; 1] µm, sai so do nghieng ve DAT; co_det neu canh giua OBB / do day > %g" % (100 * NG_CONG_43B, NG_DET_43B)})
    KQ.append(ket_qua(l, tong_hop(ca, False),
                      "Dài = cạnh lớn nhất OBB, hoặc dài theo TRỤC GIỮA (khung xương) nếu dài hơn OBB quá NG_CONG_43B (hạt cong) ∈ [1, 2] µm. Rộng = ĐỘ DÀY CỤC BỘ = 2 × bán kính cầu nội tiếp lớn nhất ∈ [0,5; 1] µm; độ dày đo là cận dưới nên nếu thiếu không quá cận sai số thì tính ĐẠT (cờ dat_trong_sai_so). Ghi kèm cạnh giữa OBB và cờ co_det (cạnh giữa/độ dày > NG_DET_43B) — nghĩa của 'rộng' (độ dày hay cạnh giữa) chờ PI quyết. Dung sai 'khoảng' = DUNG_SAI_KHOANG.",
                      ca, bat_dong="Trường máy ti_le_truc_rong_tren_dai [0,25; 0,6] không có trong lời luật — không dùng."))

    l = L_("4.3c"); ca = []
    ca.append({"bo_phan": "ti_the", "dat": trong_khoang(len(ti), *SO_TI_THE), "gia_tri": len(ti), "nguong": "%s ti the co mat trong mesh (dv.so_ti_the)" % (SO_TI_THE,),
               "them": {"so_mat": len(tb.mat_theo_loai.get("ti_the", []))}})
    for b in ti:
        if b.so_khoi != 1:
            ca.append({"bo_phan": b.ten, "dat": False, "gia_tri": b.so_khoi, "nguong": "moi ti the la 1 hat"})
    n_chong = 0
    for (a, c), k in cap.items():
        if tb.bp[a].loai == "ti_the" and tb.bp[c].loai == "ti_the" and k.get("giao_um3", 0) > EPS_VOL:
            n_chong += 1
            ca.append({"bo_phan": [a, c], "dat": False, "gia_tri": "chong lan %s µm³" % k["giao_um3"], "nguong": "hat roi (khong giao)"})
    KQ.append(ket_qua(l, tong_hop(ca, False), "Số ti thể có mặt ∈ [10, 20]; 'hạt rời' = mỗi ti thể một khối và không cặp ti thể nào giao nhau (bool).", ca))

    l = L_("4.3d"); ca = ca_mat(tb, "ti_the", "4.3d")
    for (a, c), k in cap.items():
        la, lc = tb.bp[a].loai, tb.bp[c].loai
        if la == "ti_the" and lc in ("ti_the", "nhan") or lc == "ti_the" and la == "nhan":
            ng_ = KHE_CUNG_LOAI_DV.get("ti_the", KHE_ROI) if la == lc else KHE_ROI
            ca.append({"bo_phan": [a, c], "dat": k["kc_um"] >= ng_ and k.get("giao_um3", 0) <= EPS_VOL,
                       "gia_tri": k["kc_um"], "chinh_xac": k.get("chinh_xac"), "giao_um3": k.get("giao_um3"), "nguong": ">= %g µm (chung.khe_cung_loai / chung.khe_roi)" % ng_})
    KQ.append(ket_qua(l, tong_hop(ca, False), "Mọi cặp ti thể–ti thể và ti thể–nhân: không giao (bool, với nhân dùng khối đặc) và khe chính xác ≥ 0,3 µm.", ca))

    l = L_("4.4a"); er = tb.theo_loai.get("luoi_noi_chat", [])
    ca = ca_mat(tb, "luoi_noi_chat", "4.4a")
    if not er:
        ca.append({"bo_phan": "luoi_noi_chat", "dat": False, "gia_tri": "khong co", "nguong": "phai co"})
    for b in er:
        if nhan is None: ca.append({"bo_phan": b.ten, "dat": False, "gia_tri": "khong co nhan", "nguong": ""}); continue
        d = min(kc_mat_mat(b.luoi, L2, 1.0)[0] for L2 in [nhan.luoi])
        ca.append({"bo_phan": b.ten, "dat": d <= TIEP_XUC_UM, "gia_tri": round(d, 6), "nguong": "lien voi mat ngoai nhan: kc <= 0,01 µm"})
    KQ.append(ket_qua(l, tong_hop(ca, True),
                      "Phần kiểm được: có lưới nội chất và nó liền (chạm, khoảng cách ≤ 0,01 µm) mặt ngoài nhân. KHÔNG kiểm được từ mặt bao: phân biệt ER hạt/trơn, ribosome, dạng tấm/ống. 'Trải khắp chất tế bào' bị luật 4.4b (mô hình chỉ vẽ quanh nhân) thay thế nên không kiểm.",
                      ca, ly_do="ER hạt/trơn, ribosome, dạng tấm–ống là chi tiết không có trong mesh bao.",
                      bat_dong="Trường máy ghi 'bộ dựng chưa hỗ trợ — chỉ miễn khe'; lời luật đòi liền màng nhân → bộ kiểm vẫn kiểm sự tiếp xúc."))

    l = L_("4.4b"); ca = ca_mat(tb, "luoi_noi_chat", "4.4b")
    for b in er:
        if nhan is None or mang is None: ca.append({"bo_phan": b.ten, "dat": False, "gia_tri": "thieu nhan/mang", "nguong": ""}); continue
        khoi_n = nhan.dac_loi if nhan.dac_loi is not None else nhan.dac
        g = float((b.man ^ khoi_n).volume()) if b.man is not None and khoi_n is not None else du_phong_the_tich(b, nhan, "giao")
        ok1 = g <= 1e-3
        ca.append({"bo_phan": b.ten, "dat": ok1, "gia_tri": round(g, 5), "nguong": "ER khong lot vao khoi dac nhan: giao <= 0,001 µm³ (dung sai mat trung nhau; nhan dong hinh thai dung loi co 2 voxel)"})
        dn = np.sqrt(nhan.luoi.kc2_diem(b.luoi.V)); dm = np.sqrt(mang.luoi.kc2_diem(b.luoi.V))
        xa = int((dn >= dm).sum())
        dmax = float(dn.max()); dtren = dmax + b.luoi.Lmax
        ok2 = dmax <= DAY_ER_UM * (1 + DUNG_SAI_44B)
        ca.append({"bo_phan": b.ten, "dat": ok2, "gia_tri": {"kc_max_dinh_toi_nhan_um": round(dmax, 4), "can_tren_um": round(dtren, 4),
                   "thong_tin_so_dinh_gan_mang_hon_nhan": xa, "thong_tin_tong_dinh": int(len(dn))},
                   "nguong": "vanh ER quanh nhan: kc max toi mat nhan <= DAY_ER_UM x (1 + DUNG_SAI_44B) = %g µm (cho duyet)" % (DAY_ER_UM * (1 + DUNG_SAI_44B))})
    KQ.append(ket_qua(l, tong_hop(ca, False),
                      "'Chỉ vẽ phần ER quanh nhân' = (i) ER nằm trong chất tế bào, không lọt vào khối đặc nhân (giao ≤ 0,001 µm³, dung sai cho mặt chung ER–nhân); (ii) ER là một vành quanh nhân: khoảng cách lớn nhất từ đỉnh ER tới mặt nhân ≤ 2,0 µm (+5% dung sai). Lời luật KHÔNG cho số; 2,0 µm lấy từ trường máy day_um (có ghi bất đồng). Số đỉnh ER gần màng hơn gần nhân chỉ ghi thông tin.",
                      ca, bat_dong="Lời luật không có số độ dày; dùng day_um = 2,0 của trường máy làm cách hiểu 'quanh nhân', thêm dung sai 5% vì đây là thông số vẽ chứ không phải số khoa học."))

    l = L_("4.5a"); go = tb.theo_loai.get("golgi", [])
    ca = ca_mat(tb, "golgi", "4.5a")
    co_kiem = bool(ca)
    if not go and not ca:
        ca.append({"bo_phan": "golgi", "dat": False, "gia_tri": "khong co", "nguong": "phai co"}); co_kiem = True
    thong_tin = {}; ly_do_45a = []
    for b in go:
        lop_ = dem_lop(b)
        thong_tin[b.ten] = {"so_khoi": b.so_khoi, "kich_thuoc_bao_um": [round(x, 3) for x in b.kich_thuoc()], "dem_lop_theo_tia": lop_}
        lop_pb = max(x["lop_pho_bien"] for x in lop_.values())
        if b.so_khoi == 1 and lop_pb >= 2:
            co_kiem = True
            ca.append({"bo_phan": b.ten, "dat": trong_khoang(lop_pb, 4, 8), "gia_tri": lop_pb, "nguong": "4–8 tui (dem lop pho bien theo tia)"})
            ly_do_45a.append("%s: đếm được %d lớp theo tia nhưng túi dính nhau nên không đo được bề rộng từng túi" % (b.ten, lop_pb))
        elif b.so_khoi == 1:
            ly_do_45a.append("%s: 1 khối; trên %s tia song song 3 trục chính, số lớp lớn nhất = %s → túi không tách được" % (
                b.ten, "/".join(str(x["so_tia_gap"]) for x in lop_.values()), "/".join(str(x["lop_max"]) for x in lop_.values())))
        if b.so_khoi >= 2:
            co_kiem = True
            ca.append({"bo_phan": b.ten, "dat": trong_khoang(b.so_khoi, 4, 8), "gia_tri": b.so_khoi, "nguong": "4–8 tui (mo hinh)"})
            for k, c in enumerate(b.comps):
                if c.volume <= 0: continue
                e = np.sort(c.bounding_box_oriented.primitive.extents)[::-1]
                ca.append({"bo_phan": "%s#%d" % (b.ten, k), "dat": trong_khoang(e[0], 0.7, 1.1), "gia_tri": round(float(e[0]), 3), "nguong": "rong tui [0,7; 1,1] µm"})
    st = tong_hop(ca, any("không đo được" in x for x in ly_do_45a), co_kiem)
    KQ.append(ket_qua(l, st,
                      "Nếu Golgi có các túi là khối rời: đếm túi (4–8) và đo bề rộng (cạnh lớn nhất) từng túi ∈ [0,7; 1,1] µm. Nếu Golgi là MỘT mặt bao liền (các túi bị gộp) thì không đếm/đo được túi → KHÔNG KIỂM ĐƯỢC.",
                      ca, ly_do="; ".join(ly_do_45a),
                      them={"thong_tin_bao": thong_tin,
                            "canh_bao": "Bao Golgi rộng %s µm > 1,1 µm: nếu các túi xếp chồng đồng trục thì bao không thể rộng hơn túi; chỉ hợp lý nếu túi xếp lệch/dải." % (
                                ", ".join(str(round(b.kich_thuoc()[0], 2)) for b in go)) if go else ""}))

    l = L_("4.5b"); tt = tb.theo_loai.get("trung_the", [])
    ca = ca_mat(tb, "golgi", "4.5b") + ca_mat(tb, "trung_the", "4.5b")
    for b in go:
        if nhan is None or mang is None: break
        dn = kc_vf(b.luoi, nhan.luoi); dm = kc_vf(b.luoi, mang.luoi)
        ca.append({"bo_phan": b.ten, "dat": dn < dm, "gia_tri": {"khe_nhan_um": round(dn, 4), "khe_mang_um": round(dm, 4)}, "nguong": "canh nhan: gan nhan hon gan mang"})
        for c in tt:
            cc = np.array(c.tm.center_mass)
            try:
                inh = bool(Delaunay(b.luoi.V[ConvexHull(b.luoi.V).vertices]).find_simplex(cc[None])[0] >= 0)
            except Exception:
                inh = False
            ca.append({"bo_phan": [b.ten, c.ten], "dat": inh, "gia_tri": "tam trung the %s bao loi Golgi" % ("TRONG" if inh else "NGOAI"),
                       "nguong": "om quanh: tam trung the nam trong bao loi Golgi"})
    KQ.append(ket_qua(l, tong_hop(ca, False),
                      "'Cạnh nhân' = khe Golgi–nhân nhỏ hơn khe Golgi–màng; 'ôm quanh trung thể' = tâm khối trung thể nằm trong bao lồi của Golgi. Lời luật không cho số; định nghĩa tương đối do bộ kiểm đặt. Luật mềm.",
                      ca))

    l = L_("4.5c"); ca = ca_mat(tb, "golgi", "4.5c")
    for b in go:
        if nhan is None: break
        k = cap.get((b.ten, nhan.ten)) or cap.get((nhan.ten, b.ten))
        ca.append({"bo_phan": [b.ten, nhan.ten], "dat": k["kc_um"] >= KHE_ROI and k.get("giao_um3", 0) <= EPS_VOL, "gia_tri": k["kc_um"], "giao_um3": k.get("giao_um3"), "nguong": ">= %g µm (chung.khe_roi)" % KHE_ROI})
    KQ.append(ket_qua(l, tong_hop(ca, False), "Golgi không giao khối đặc nhân (bool) và khe chính xác tới nhân ≥ 0,3 µm.", ca))

    l = L_("4.5d"); ca = []
    for b in go:
        if nhan is None or not er: break
        P = b.luoi.V[:: max(1, len(b.luoi.V) // 300)]
        _, _, C = nhan.luoi.aabb().squared_distance(nhan.luoi.V, nhan.luoi.F, P)
        pts = []; idx = []
        for k_, (p, q) in enumerate(zip(P, C)):
            n_ = max(3, int(np.ceil(np.linalg.norm(q - p) / 0.02)) + 1)
            S_ = p + np.linspace(0, 1, n_)[1:-1, None] * (q - p)
            pts.append(S_); idx.append(np.full(len(S_), k_))
        pts = np.vstack(pts); idx = np.concatenate(idx)
        trong = np.zeros(len(pts), bool)
        for e in er:
            trong |= e.luoi.quan(pts) > 0.5
        cat = len(np.unique(idx[trong]))
        f = cat / len(P)
        ca.append({"bo_phan": b.ten, "dat": f >= NG_45D_TI_LE, "gia_tri": round(f, 3), "nguong": ">= NG_45D_TI_LE = %g doan ngan nhat Golgi→nhan cat qua ER" % NG_45D_TI_LE})
    KQ.append(ket_qua(l, tong_hop(ca, True, bool(ca)),
                      "Phần kiểm được: ER xen giữa = ít nhất một nửa số đoạn ngắn nhất từ điểm trên mặt Golgi tới nhân đi xuyên khối ER (bước 0,02 µm, số vòng quấn chính xác). 'Túi vận chuyển' không có loại bộ phận nào → không kiểm. Luật mềm, chưa có nguồn.",
                      ca, ly_do="Không có bộ phận 'túi vận chuyển' trong mesh."))

    l = L_("4.6a"); ly = tb.theo_loai.get("lysosome", [])
    ca = ca_mat(tb, "lysosome", "4.6a")
    for b in ly:
        e = b.kich_thuoc()
        ca.append({"bo_phan": b.ten, "dat": all(trong_khoang(x, 0.1, 1.2) for x in e), "gia_tri": [round(x, 3) for x in e], "nguong": "[0,1; 1,2] µm moi canh"})
    KQ.append(ket_qua(l, tong_hop(ca, False), "Đường kính lysosome = các cạnh hộp định hướng nhỏ nhất, mọi cạnh ∈ [0,1; 1,2] µm ('thường 0,1–0,5' chỉ là thông tin).", ca))

    l = L_("4.6d"); ca = ca_mat(tb, "lysosome", "4.6d")
    for b in ly:
        e = b.kich_thuoc(); ti_ = e[2] / e[0]
        ca.append({"bo_phan": b.ten, "dat": ti_ >= NG_46D_CAU and b.so_khoi == 1, "gia_tri": round(ti_, 3), "nguong": "canh nho / canh lon >= NG_46D_CAU = %g (cho duyet)" % NG_46D_CAU})
    KQ.append(ket_qua(l, tong_hop(ca, False), "'Hình cầu' = một khối, tỉ số cạnh hộp định hướng nhỏ nhất/lớn nhất ≥ 0,85. Lời luật không cho số; 0,85 do bộ kiểm đặt (cho phép lưới ít đỉnh).", ca))

    l = L_("4.6b"); ca = []
    if ly and nhan is not None and mang is not None:
        u = []
        for b in ly:
            c = np.array(b.tm.center_mass)[None]
            dn = float(np.sqrt(nhan.luoi.kc2_diem(c))[0]); dm = float(np.sqrt(mang.luoi.kc2_diem(c))[0])
            u.append(dn / (dn + dm + 1e-12))
        u = np.array(u)
        ca.append({"bo_phan": "lysosome", "dat": bool(np.median(u) < 0.5), "gia_tri": round(float(np.median(u)), 3), "nguong": "trung vi vi tri tuong doi u=d_nhan/(d_nhan+d_mang) < 0,5 (tap trung quanh nhan)"})
        ca.append({"bo_phan": "lysosome", "dat": bool(u.max() >= 0.5), "gia_tri": round(float(u.max()), 3), "nguong": "co it nhat 1 lysosome u >= 0,5 (rai ra phia ngoai)"})
    else:
        ca.append({"bo_phan": "lysosome", "dat": False, "gia_tri": "thieu lysosome/nhan/mang", "nguong": ""})
    KQ.append(ket_qua(l, tong_hop(ca, False),
                      "Vị trí tương đối của tâm mỗi lysosome u = d(nhân)/(d(nhân)+d(màng)): 'tập trung quanh nhân' = trung vị u < 0,5; 'rải ra phía ngoài' = có lysosome u ≥ 0,5. Định nghĩa do bộ kiểm đặt; luật mềm.",
                      ca))

    l = L_("4.6c"); ca = [{"bo_phan": "lysosome", "dat": trong_khoang(len(ly), *SO_LYSOSOME), "gia_tri": len(ly), "nguong": "%s (dv.so_lysosome)" % (SO_LYSOSOME,),
                          "them": {"so_mat": len(tb.mat_theo_loai.get("lysosome", []))}}]
    KQ.append(ket_qua(l, tong_hop(ca, False), "Số lysosome có mặt trong mesh ∈ [3, 8].", ca))

    l = L_("4.7"); rib = tb.theo_loai.get("ribosome", [])
    KQ.append(ket_qua(l, "KHONG_KIEM_DUOC",
                      "Chấm ribosome trên ER hạt: cần bộ phận ribosome hoặc chi tiết mặt ER.",
                      [], ly_do="Mesh là mặt bao, không có loại 'ribosome'; chấm 0,04 µm trên ER (nếu có) bị gộp vào bao ER; luật lại là tuỳ chọn bật/tắt của giáo viên."))

    l = L_("4.8a"); ca = ca_mat(tb, "trung_the", "4.8a"); co_kiem = bool(ca)
    if not tt and not ca:
        ca.append({"bo_phan": "trung_the", "dat": False, "gia_tri": "khong co", "nguong": "phai co"}); co_kiem = True
    tt_info = {}; ly_do_48a = []
    for b in tt:
        e_ = b.kich_thuoc()
        tt_info[b.ten] = {"so_khoi": b.so_khoi, "kich_thuoc_bao_um": [round(x, 3) for x in e_]}
        if b.so_khoi == 1:
            ly_do_48a.append("%s: 1 khối, hộp bao %s µm (cạnh nhỏ/lớn = %.2f), không có hai khối trung tử để đo" % (
                b.ten, [round(x, 3) for x in e_], e_[2] / e_[0]))
        if b.so_khoi >= 3 or b.so_khoi == 0:
            co_kiem = True
            ca.append({"bo_phan": b.ten, "dat": False, "gia_tri": b.so_khoi, "nguong": "dung 2 trung tu (1 khoi bao thi khong kiem duoc)"})
        if b.so_khoi == 2:
            co_kiem = True; ax = []
            for k, c in enumerate([c for c in b.comps if c.volume > 0]):
                e = np.sort(c.bounding_box_oriented.primitive.extents)[::-1]
                ok = trong_khoang(e[0], 0.3, 0.5) and trong_khoang(e[1], 0.15, 0.25) and trong_khoang(e[2], 0.15, 0.25)
                ca.append({"bo_phan": "%s#%d" % (b.ten, k), "dat": ok, "gia_tri": [round(float(x), 3) for x in e], "nguong": "dai [0,3;0,5], rong [0,15;0,25]"})
                v = c.vertices - c.vertices.mean(0); w, U = np.linalg.eigh(np.cov(v.T)); ax.append(U[:, -1])
            g = float(np.degrees(np.arccos(abs(np.clip(ax[0] @ ax[1], -1, 1)))))
            ca.append({"bo_phan": b.ten, "dat": g >= NG_48A_GOC, "gia_tri": round(g, 1), "nguong": "goc hai truc >= NG_48A_GOC = %g°" % NG_48A_GOC})
    KQ.append(ket_qua(l, tong_hop(ca, True, co_kiem),
                      "Trung thể 2 khối: đo từng trung tử (dài, rộng) và góc hai trục. 1 khối (mặt bao gộp): không kiểm được. ≥ 3 khối: TRƯỢT (không phải hai trung tử).",
                      ca, ly_do="; ".join(ly_do_48a),
                      them={"thong_tin_bao": tt_info}))

    l = L_("4.8b"); ca = ca_mat(tb, "trung_the", "4.8b")
    for b in tt:
        if nhan is None or mang is None: break
        dn = kc_vf(b.luoi, nhan.luoi); dm = kc_vf(b.luoi, mang.luoi)
        ca.append({"bo_phan": b.ten, "dat": dn < dm, "gia_tri": {"khe_nhan_um": round(dn, 4), "khe_mang_um": round(dm, 4), "goi_y_1_5_um": dn <= 1.5},
                   "nguong": "gan nhan hon gan mang (1,5 µm chi la thong tin)"})
    KQ.append(ket_qua(l, tong_hop(ca, False), "'Gần nhân' = khe trung thể–nhân nhỏ hơn khe trung thể–màng. Số 1,5 µm là trung bình 2D, luật ghi không dùng làm giới hạn cứng → chỉ báo. Luật mềm.", ca,
                      bat_dong="Trường máy goi_y_um=1,5 — lời phạm vi nói không dùng làm giới hạn cứng; chỉ báo."))

    l = L_("4.8c"); ca = ca_mat(tb, "trung_the", "4.8c")
    for b in tt:
        if nhan is None: break
        k = cap.get((b.ten, nhan.ten)) or cap.get((nhan.ten, b.ten))
        ca.append({"bo_phan": [b.ten, nhan.ten], "dat": k["kc_um"] >= KHE_ROI and k.get("giao_um3", 0) <= EPS_VOL, "gia_tri": k["kc_um"], "giao_um3": k.get("giao_um3"), "nguong": ">= %g µm (chung.khe_roi)" % KHE_ROI})
    KQ.append(ket_qua(l, tong_hop(ca, False), "Trung thể không giao khối đặc nhân và khe chính xác ≥ 0,3 µm.", ca))

    l = L_("M.a"); ca = []
    for loai_ in sorted(set(tb.du_kien.values()) - {"mang"}):
        ca += ca_mat(tb, loai_, "M.a")
    for b in bao_quan:
        if mang is None: ca.append({"bo_phan": b.ten, "dat": False, "gia_tri": "khong co mang", "nguong": ""}); continue
        r = trong_mang[b.ten]
        ca.append({"bo_phan": b.ten, "dat": r["trong"], "gia_tri": r["ngoai_um3"], "nguong": "the tich ngoai mang <= %g µm³" % EPS_VOL})
    KQ.append(ket_qua(l, tong_hop(ca, False), "Mọi bộ phận khác màng: thể tích (bộ phận − khối màng) ≤ ε (bool chính xác). Bộ phận bị mất = trượt.", ca))

    l = L_("M.b"); ca = []
    for loai_ in sorted(set(tb.du_kien.values()) - {"mang"}):
        ca += ca_mat(tb, loai_, "M.b")
    for b in bao_quan:
        if mang is None: continue
        r = trong_mang[b.ten]
        kc = r["kc_um"]
        ca.append({"bo_phan": b.ten, "dat": bool(r["trong"] and kc is not None and kc >= T_MANG), "gia_tri": kc, "chinh_xac": r.get("chinh_xac"), "nguong": ">= %g µm toi mat mang (dv.T_mang)" % T_MANG})
    KQ.append(ket_qua(l, tong_hop(ca, False), "Mọi bào quan nằm trong màng và khe chính xác tới mặt màng ≥ 0,5 µm (kể cả nhân, ER, nhân con).", ca))

    l = L_("M.c"); ca = []
    for loai_ in sorted(set(tb.du_kien.values()) - {"mang"}):
        ca += ca_mat(tb, loai_, "M.c")
    for (a, c), k in cap.items():
        if tb.bp[a].loai == tb.bp[c].loai:
            continue
        ca.append({"bo_phan": [a, c], "dat": k["kc_um"] >= KHE_ROI and k.get("giao_um3", 0) <= EPS_VOL, "gia_tri": k["kc_um"], "giao_um3": k.get("giao_um3"), "chinh_xac": k.get("chinh_xac"), "nguong": ">= %g µm (chung.khe_roi)" % KHE_ROI})
    KQ.append(ket_qua(l, tong_hop(ca, False),
                      "Mọi cặp bộ phận KHÁC loại (trừ màng): không giao (bool; nhân dùng khối đặc nên vật nằm trong khoang nhân rỗng là giao) và khe chính xác ≥ 0,3 µm. Miễn: ER–nhân (4.4a buộc liền), nhân con–nhân (4.2a/4.2b quy định riêng), ribosome–ER (4.7). Cặp cùng loại không thuộc luật này.",
                      ca, bat_dong="Lời luật không ghi miễn; miễn suy từ 4.4a, 4.2a–b, 4.7 (trùng trường máy 'mien'). Miễn trung_tu–trung_tu không cần vì trung tử không là bộ phận riêng."))

    l = L_("M.d"); ca = []
    if nhan is not None:
        lon = max(nhan.comps, key=lambda c: abs(c.volume)); g = int(round((2 - lon.euler_number) / 2))
        if not nhan.co_khoang:
            KQ.append(ket_qua(l, "KHONG_KIEM_DUOC", "Màng nhân hai lớp, khoang, lỗ nhân: chỉ kiểm được khi mesh nhân là VỎ bao quanh một khoang kín (nhân rỗng). Nhân đặc (kể cả có đường hầm) là mặt bao, không thể hiện các lớp màng.", [],
                              ly_do="Mesh nhân: giống mặt ngoài = %d, có khoang kín = không (cách dựng khối đặc: %s) — hai lớp/khoang/lỗ không thể hiện trong mặt bao." % (g, (nhan.dac_cach or {}).get("cach")),
                              them={"giong_mat_ngoai": g}))
        else:
            c0 = nhan.luoi.V.mean(0); dirs = _fibo(400); R0 = float(np.linalg.norm(nhan.hi - nhan.lo))
            A = nhan.luoi.aabb(); n_lop2 = 0; n_xuyen = 0; n_tia = 0
            O = c0 + dirs * R0; D = -dirs
            hits = A.intersect_ray(nhan.luoi.V, nhan.luoi.F, O, D)
            for h_ in hits:
                ts = sorted(set(round(x[1] * 1.0, 9) for x in h_))
                ts = [t * 1.0 for t in ts]
                if not ts: continue
                n_tia += 1
                if len(ts) >= 4 and (ts[2] - ts[1]) <= MD_KHOANG_MAX:
                    n_lop2 += 1
            ca.append({"bo_phan": nhan.ten, "dat": g >= 1, "gia_tri": g, "nguong": "co lo nhan: mat vo co giong >= 1"})
            f2 = n_lop2 / max(1, n_tia)
            ca.append({"bo_phan": nhan.ten, "dat": f2 >= MD_TI_LE_TIA, "gia_tri": round(f2, 3), "nguong": ">= MD_TI_LE_TIA = %g tia gap hai lop cach nhau <= MD_KHOANG_MAX = %g µm" % (MD_TI_LE_TIA, MD_KHOANG_MAX)})
            KQ.append(ket_qua(l, tong_hop(ca, False), "Mesh nhân có cấu trúc vỏ → kiểm: (i) có lỗ = giống mặt vỏ ≥ 1; (ii) hai lớp = theo tia từ ngoài vào, hai lớp liên tiếp cách nhau ≤ 0,3 µm ở ≥ 50% tia. Độ rộng 20–40 nm không bắt buộc vì 'mo_hinh' của luật cho vẽ to hơn thật.", ca,
                              them={"giong_mat_ngoai": g, "uoc_so_lo": g + 1}))
    else:
        KQ.append(ket_qua(l, "TRUOT", "Không có nhân.", [{"bo_phan": "nhan", "dat": False, "gia_tri": "mat", "nguong": "phai co"}]))

    l = L_("M.e")
    KQ.append(ket_qua(l, "KHONG_KIEM_DUOC", "Mào ti thể là màng trong — chi tiết bên trong.", [],
                      ly_do="Mesh ti thể là mặt bao ngoài; màng trong và mào không có trong mesh."))

    l = L_("M.f"); ca = []
    for loai_ in sorted(set(tb.du_kien.values()) - LIEN_KHOI_MIEN):
        ca += ca_mat(tb, loai_, "M.f")
    for b in tb.bp.values():
        if b.loai in LIEN_KHOI_MIEN:
            continue
        if b.loai == "trung_the" and b.so_khoi == 2:
            ca.append({"bo_phan": b.ten, "dat": True, "gia_tri": 2, "nguong": "trung the = 2 trung tu (4.8a) — mien"}); continue
        ok = b.so_khoi == 1 and not b.loi_chan
        ca.append({"bo_phan": b.ten, "dat": ok, "gia_tri": b.so_khoi if not b.loi_chan else {"so_khoi": b.so_khoi, "loi_chan": b.loi_chan},
                   "nguong": "1 khoi dac hop le (vo co khoang trong tinh la 1; khong tu cat, khong vo trung lat nguoc)", "loi_luoi": b.loi})
    KQ.append(ket_qua(l, tong_hop(ca, False), "Mỗi bộ phận (trừ Golgi) = đúng một khối đặc HỢP LỆ: đếm thành phần mặt có thể tích dương (mặt khoang có thể tích âm, nằm hẳn trong thành phần dương, không tính thêm khối); lưới tự cắt / vỏ trùng lật ngược / khoang hở = trượt. Bộ phận bị mất = trượt.", ca))

    l = {"id": "M.h", "phat_bieu": "Hai bộ phận cùng loại không chồng lên nhau (được chạm, không được giao).",
         "nhan": "sự thật", "luat_cung": True, "lop": ["KHTN6", "Sinh10"]}
    luat["M.h"] = l
    ca = [{"bo_phan": x["bo_phan"], "dat": False, "gia_tri": x["giao_um3"], "nguong": "the tich giao <= %g µm³" % EPS_VOL} for x in cl_cung]
    KQ.append(ket_qua(l, tong_hop(ca, False),
                      "Mọi cặp bộ phận cùng loại (kể cả ti thể–ti thể, nhân con–nhân con, lysosome–lysosome; trừ màng): thể tích giao chính xác (bool manifold3d) ≤ EPS_VOL = %g µm³. Chạm (giao = 0) được phép. Bộ phận bị mất không tạo ca (không thể chồng)." % EPS_VOL,
                      ca, them={"so_cap_cung_loai_da_xet": sum(1 for i in range(len(bao_quan)) for j in range(i + 1, len(bao_quan)) if bao_quan[i].loai == bao_quan[j].loai and hop_cach(bao_quan[i], bao_quan[j]) == 0)}))

    assert len(KQ) == 35, len(KQ)
    for r in KQ:
        lop_luat = luat[r["id"]].get("lop") or []
        if lop not in lop_luat:
            r["neu_ap_dung"] = r["trang_thai"]; r["trang_thai"] = "KHONG_AP_DUNG"
            r["ly_do"] = ("Luật dành cho %s; mô hình kiểm theo lớp %s. " % ("/".join(lop_luat), lop)) + r.get("ly_do", "")
        mat_ = [c for c in r["ca_truot"] if str(c.get("gia_tri", "")).startswith("MAT")]
        r["so_ca_truot_mat"] = len(mat_); r["so_ca_truot_hinh"] = r["so_ca_truot"] - len(mat_)
        if r["trang_thai"] == "TRUOT":
            r["loai_truot"] = "HINH_HOC" if r["so_ca_truot_hinh"] > 0 else "CHI_DO_MAT_BO_PHAN"
    dem = collections.Counter(r["trang_thai"] for r in KQ)
    dem_tach = collections.Counter(r.get("loai_truot") for r in KQ if r["trang_thai"] == "TRUOT")
    loi_luoi = {b.ten: {"da_xu_ly": b.loi, "CHAN": b.loi_chan} for b in tb.bp.values() if b.loi or b.loi_chan}
    out = {"te_bao": os.path.basename(os.path.normpath(thu_muc)), "lop": lop, "so_luat": len(KQ),
           "dem_trang_thai": dict(dem), "truot_tach": dict(dem_tach),
           "bo_phan_loi_chan": sorted(b.ten for b in tb.bp.values() if b.loi_chan),
           "chong_lan_cung_loai": {"luat": "M.h", "ghi": "da thanh luat M.h (PI duyet 05/11)",
                                   "so_cap": len(cl_cung), "cap": cl_cung},
           "tham_so_cho_duyet": CHO_DUYET,
           "luat_cung_truot": [r["id"] for r in KQ if r["trang_thai"] == "TRUOT" and r["luat_cung"]],
           "luat_mem_truot": [r["id"] for r in KQ if r["trang_thai"] == "TRUOT" and not r["luat_cung"]],
           "bo_phan_mat": tb.mat, "bo_phan_thua": tb.thua, "so_bo_phan": len(tb.bp),
           "loi_luoi": loi_luoi,
           "khoi_dac_nhan": nhan.dac_cach if nhan is not None else None,
           "khoi_dac_mang": mang.dac_cach if mang is not None else None,
           "so_cap_khe_da_xet": len(cap), "thoi_gian_s": round(time.time() - t0, 2), "thoi_gian_cap_s": round(t_cap, 2),
           "tham_so_ky_thuat": {"EPS_VOL_um3": EPS_VOL, "H_VOXEL_um": H_VOXEL, "R_DONG_MAX_um": R_DONG_MAX,
                                "V_KHOANG_MIN_um3": V_KHOANG_MIN, "H_DAY_um": H_DAY},
           "luat": KQ}
    global _NGU_CANH
    _NGU_CANH = {k_: v_ for k_, v_ in locals().items() if k_ not in ("fo",)}
    with open(os.path.join(thu_muc, "kiem2.json"), "w") as fo:
        json.dump(out, fo, ensure_ascii=False, indent=1, default=_json_default)
    if not im_lang:
        in_tom_tat(out)
    return out


def dem_lop(bp, n_tia=15, buoc=0.01):
    """Đếm số LỚP (đoạn 'trong' rời nhau) dọc tia song song ba trục chính của bộ phận.
    Trả {truc: {"so_tia_gap": k, "lop_max": m, "lop_pho_bien": mode}} — dùng cho 4.5a (túi Golgi)."""
    v = bp.luoi.V; c = v.mean(0); _, U = np.linalg.eigh(np.cov((v - c).T))
    out = {}
    for ax in range(3):
        a = U[:, ax]; b1 = U[:, (ax + 1) % 3]; b2 = U[:, (ax + 2) % 3]
        pa = (v - c) @ a; p1 = (v - c) @ b1; p2 = (v - c) @ b2
        ts = np.arange(pa.min() - 2 * buoc, pa.max() + 2 * buoc, buoc)
        g1 = np.linspace(p1.min(), p1.max(), n_tia + 2)[1:-1]; g2 = np.linspace(p2.min(), p2.max(), n_tia + 2)[1:-1]
        S1, S2 = np.meshgrid(g1, g2, indexing='ij')
        Q = (c[None, None, None] + ts[None, None, :, None] * a + S1[..., None, None] * b1 + S2[..., None, None] * b2).reshape(-1, 3)
        ins = (bp.luoi.quan(Q) > 0.5).reshape(n_tia * n_tia, len(ts))
        runs = (np.diff(ins.astype(np.int8), axis=1) == 1).sum(1) + ins[:, 0]
        gap = runs[runs > 0]
        out["truc_%d" % ax] = {"so_tia_gap": int(len(gap)), "lop_max": int(gap.max()) if len(gap) else 0,
                               "lop_pho_bien": int(np.bincount(gap).argmax()) if len(gap) else 0}
    return out


def _json_default(o):
    if isinstance(o, (np.integer,)): return int(o)
    if isinstance(o, (np.floating,)): return float(o)
    if isinstance(o, (np.bool_,)): return bool(o)
    if isinstance(o, np.ndarray): return o.tolist()
    return str(o)


def _mesh(m):
    mm = m.to_mesh64()
    return np.asarray(mm.vert_properties)[:, :3], np.asarray(mm.tri_verts)


def _fibo(n):
    i = np.arange(n) + 0.5; phi = np.arccos(1 - 2 * i / n); th = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], 1)


def in_tom_tat(out):
    print("== %s  (lop %s, %d bo phan, mat %d)  %.1f s" % (out["te_bao"], out["lop"], out["so_bo_phan"], len(out["bo_phan_mat"]), out["thoi_gian_s"]))
    print("   dem:", out["dem_trang_thai"], "/ %d" % out["so_luat"], " truot:", out["truot_tach"])
    if out["bo_phan_loi_chan"]: print("   LOI CHAN:", out["bo_phan_loi_chan"][:10])
    print("   M.h chong lan cung loai: %d cap" % out["chong_lan_cung_loai"]["so_cap"], out["chong_lan_cung_loai"]["cap"][:3])
    for r in out["luat"]:
        x = ""
        if r["trang_thai"] == "TRUOT":
            c = r["ca_truot"][0] if r["ca_truot"] else {}
            gt = c.get("gia_tri")
            if c.get("chinh_xac") is False and isinstance(gt, (int, float)):
                gt = ">=%s" % gt
            x = " — %d ca truot; vd %s: %s (nguong %s)" % (r["so_ca_truot"], c.get("bo_phan"), str(gt)[:80], c.get("nguong"))
        print("   %-5s %-15s%s%s" % (r["id"], r["trang_thai"], "" if r["luat_cung"] else "(mem)", x))
    if out["loi_luoi"]:
        print("   loi luoi:", out["loi_luoi"])
    if out["khoi_dac_nhan"]:
        print("   khoi dac nhan:", {k: v for k, v in out["khoi_dac_nhan"].items() if k != "nhat_ky_tim_r"})


if os.environ.get("KIEM_TANG_TOC", "1") != "0":
    import tang_toc as _tang_toc
    _tang_toc.cai(sys.modules[__name__])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("thu_muc"); ap.add_argument("--lop", default=None); ap.add_argument("--luat", default=None)
    ap.add_argument("--vat", default="dong_vat", choices=["dong_vat", "thuc_vat"])
    a = ap.parse_args()
    if a.vat == "thuc_vat":
        import kiem_thucvat
        kiem_thucvat.kiem_tv(a.thu_muc, a.lop, a.luat)
    else:
        kiem(a.thu_muc, a.lop, a.luat)
