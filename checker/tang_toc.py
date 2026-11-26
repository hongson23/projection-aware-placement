# -*- coding: utf-8 -*-
"""
tang_toc.py — TĂNG TỐC bộ kiểm độc lập (kiem2.py / kiem_thucvat.py) KHÔNG đổi ngữ nghĩa.

Cài bằng một dòng ở cuối kiem2.py / kiem_thucvat.py (tắt bằng KIEM_TANG_TOC=0 → mã gốc nguyên vẹn).
Mọi thay đổi đều là CÙNG phép tính, chỉ đổi cách tìm:

 1) Luoi.__init__: np.unique(E, axis=0) (chậm) → np.unique trên khoá 1 chiều a·n + b. Cùng tập cạnh,
    cùng thứ tự từ điển (a < n, b < n) → mảng E giống hệt từng phần tử.

 2) kc_mat_mat, bước cạnh–cạnh: bản gốc liệt kê mọi cặp cạnh có trung điểm cách < lim + (Lmax_A + Lmax_B)/2
    bằng KD-tree (hàng tỉ cặp khi hai mặt lớn song song cách ~lim) rồi tính khoảng cách đoạn–đoạn.
    Kết quả gốc chỉ phụ thuộc vào các cặp có khoảng cách d < lim (cặp d ≥ lim không đổi được giá trị trả về:
    lim = min(U, R); nếu lim = U thì best ≤ U; nếu lim = R < U thì best ≥ R ⇒ trả (R, False) như nhau),
    và mọi cặp d < lim đều nằm trong tập lọc gốc (lọc gốc bảo thủ). Bản mới: trên GPU (float64) lọc cặp
    bằng hộp bao theo khối (sắp Morton) → cận dưới trung điểm theo TỪNG cạnh → khoảng cách đoạn–đoạn fp64;
    giữ cặp d_gpu < lim + 1e-7 rồi TÍNH LẠI chúng bằng chính kiem2.kc_doan_doan (numpy, CPU) ⇒ giá trị
    trả về giống hệt bản gốc. Không có GPU → chạy mã gốc.

 3) Ghi dần + chạy tiếp: mỗi luật xong ghi vào <thu_muc>/kiem2_dang_chay.json; các phép tính đắt
    (khe, kc, giao, chứa, khoảng cách mặt–mặt, thể tích bool/dự phòng) được ghi đệm vào
    <thu_muc>/kiem2_bo_dem.pkl theo khoá = băm nội dung lưới + tham số; chạy lại sau khi đứt thì lấy từ đệm
    (cùng giá trị đã tính). Đệm vô hiệu khi mã, mesh.npz hay loai.json đổi.
"""
import os, sys, json, time, pickle, hashlib, functools
import numpy as np

PHIEN_BAN = "tang_toc_v2"
_JSON_DEFAULT = str
_GPU = None
_THONG_KE = {"canh_canh_gpu": 0, "canh_canh_cpu": 0, "cap_tinh_lai": 0, "dem_trung": 0, "dem_ghi": 0,
             "t_canh_canh_gpu_s": 0.0}


def _gpu():
    global _GPU
    if _GPU is None:
        if os.environ.get("KIEM_GPU", "1") == "0":
            _GPU = False
        else:
            try:
                import torch
                _GPU = torch.device("cuda") if torch.cuda.is_available() else False
            except Exception:
                _GPU = False
    return _GPU


def _luoi_init_nhanh(self, V, F):
    self.V = np.ascontiguousarray(V, dtype=np.float64)
    self.F = np.ascontiguousarray(F, dtype=np.int64)
    E = np.vstack([self.F[:, [0, 1]], self.F[:, [1, 2]], self.F[:, [2, 0]]])
    E.sort(axis=1)
    if len(E):
        n = int(E.max()) + 1
        u = np.unique(E[:, 0] * n + E[:, 1])
        self.E = np.stack([u // n, u % n], 1).astype(E.dtype, copy=False)
    else:
        self.E = np.unique(E, axis=0)
    L = np.linalg.norm(self.V[self.E[:, 0]] - self.V[self.E[:, 1]], axis=1)
    self.Lmax = float(L.max()) if len(L) else 0.0
    self.lo = self.V.min(0); self.hi = self.V.max(0)
    self._aabb = None; self._kdt = None; self._kde = None


def _quan_nhanh(self, P):
    """= igl.fast_winding_number(V, F, P) nhưng dựng cây BVH MỘT lần cho mỗi lưới (igl.FastWindingNumberBVH,
    accuracy_scale = 2 như mặc định). Đã kiểm: GIỐNG HỆT TỪNG BIT với lời gọi gốc trên không bào thật
    (500 000 điểm, kể cả sát mặt). Tránh dựng lại cây (và bùng luồng lồng nhau) mỗi lần gọi."""
    import igl
    P = np.ascontiguousarray(P, dtype=np.float64)
    B = getattr(self, "_fwn", None)
    if B is None:
        try:
            B = igl.FastWindingNumberBVH(); B.init(self.V, self.F); self._fwn = B
        except Exception:
            self._fwn = False; B = False
    if B is False or len(P) == 0:
        return igl.fast_winding_number(self.V, self.F, P)
    return B.winding_number(P, 2.0)


def _morton(P, lo, hi):
    import torch
    q = ((P - lo) / torch.clamp(hi - lo, min=1e-12) * 1023).clamp(0, 1023).to(torch.int64)
    def tach(x):
        x = (x | (x << 16)) & 0x030000FF
        x = (x | (x << 8)) & 0x0300F00F
        x = (x | (x << 4)) & 0x030C30C3
        x = (x | (x << 2)) & 0x09249249
        return x
    return tach(q[:, 0]) | (tach(q[:, 1]) << 1) | (tach(q[:, 2]) << 2)


def _doan_doan_t(P1, Q1, P2, Q2):
    """Bản torch (fp64) của kiem2.kc_doan_doan — cùng công thức, chỉ để chọn ứng viên."""
    import torch
    d1 = Q1 - P1; d2 = Q2 - P2; r = P1 - P2
    a = (d1 * d1).sum(-1); e = (d2 * d2).sum(-1); f = (d2 * r).sum(-1); c = (d1 * r).sum(-1); b = (d1 * d2).sum(-1)
    tiny = 1e-20
    one = torch.ones_like(a)
    a_ = torch.where(a > tiny, a, one); e_ = torch.where(e > tiny, e, one)
    den = a * e - b * b
    s = torch.where(den > tiny, torch.clamp((b * f - c * e) / torch.where(den > tiny, den, one), 0, 1), torch.zeros_like(a))
    t = (b * s + f) / e_
    s = torch.where(t < 0, torch.clamp(-c / a_, 0, 1), torch.where(t > 1, torch.clamp((b - c) / a_, 0, 1), s))
    t = torch.clamp(t, 0, 1)
    s = torch.where(a > tiny, s, torch.zeros_like(s))
    t = torch.where(e > tiny, t, torch.clamp(f / e_, 0, 1) * 0.0)
    t = torch.where((a <= tiny) & (e > tiny), torch.clamp(f / e_, 0, 1), t)
    C1 = P1 + d1 * s[..., None]; C2 = P2 + d2 * t[..., None]
    return torch.linalg.norm(C1 - C2, dim=-1)


def _cap_ung_vien_gpu(VA, eA, VB, eB, lim, dev, S=64, CT=128, CP=2048):
    """Trả (ia, jb): chỉ số (vào eA, eB) của các cặp cạnh có khoảng cách đoạn–đoạn (fp64 GPU) ≤ min_gpu + 1e-7
    và ≤ lim + 1e-7. Lọc BẢO THỦ hai tầng, đều là CẬN DƯỚI của khoảng cách đoạn–đoạn:
      (1) nhóm 64 cạnh liên tiếp theo đường Morton → hộp bao nhóm; cặp nhóm có khoảng cách hộp > ngưỡng bị loại;
      (2) trong cặp nhóm: |m_a − m_b| − l_a/2 − l_b/2 (m trung điểm, l độ dài) > ngưỡng bị loại.
    Hai tầng lọc tính fp32 trên toạ độ đã trừ tâm (sai số làm tròn ≤ ~1e-5 µm) với biên M32 = 1e-4 µm ⇒ không loại
    nhầm cặp nào có khoảng cách < ngưỡng. Cặp còn lại: khoảng cách đoạn–đoạn fp64 (cùng công thức kc_doan_doan)."""
    import torch
    m = 1e-7; M32 = 1e-4; _t0 = time.time()
    DO = os.environ.get("KIEM_DO_GPU") == "1"
    A0 = torch.as_tensor(VA[eA[:, 0]], device=dev); A1 = torch.as_tensor(VA[eA[:, 1]], device=dev)
    B0 = torch.as_tensor(VB[eB[:, 0]], device=dev); B1 = torch.as_tensor(VB[eB[:, 1]], device=dev)
    mA = 0.5 * (A0 + A1); mB = 0.5 * (B0 + B1)
    lo = torch.minimum(mA.min(0).values, mB.min(0).values); hi = torch.maximum(mA.max(0).values, mB.max(0).values)
    c0 = 0.5 * (lo + hi)
    oA = torch.argsort(_morton(mA, lo, hi)); oB = torch.argsort(_morton(mB, lo, hi))
    A0, A1, mA = A0[oA], A1[oA], mA[oA]
    B0, B1, mB = B0[oB], B1[oB], mB[oB]
    def nhom(X0, X1, M):
        n = X0.shape[0]; k = (n + S - 1) // S; pad = k * S - n
        l = (0.5 * torch.linalg.norm(X1 - X0, dim=1)).float()
        Ml = (M - c0).float()
        blo = (torch.minimum(X0, X1) - c0).float(); bhi = (torch.maximum(X0, X1) - c0).float()
        if pad:
            Ml = torch.cat([Ml, Ml[-1:].expand(pad, 3)]); l = torch.cat([l, torch.full((pad,), -float("inf"), device=dev)])
            blo = torch.cat([blo, blo[-1:].expand(pad, 3)]); bhi = torch.cat([bhi, bhi[-1:].expand(pad, 3)])
        return Ml.view(k, S, 3), l.view(k, S), blo.view(k, S, 3).min(1).values, bhi.view(k, S, 3).max(1).values
    MA, LA, tlo, thi = nhom(A0, A1, mA); MB, LB, blo, bhi = nhom(B0, B1, mB)
    nA = A0.shape[0]; nB = B0.shape[0]
    if DO: torch.cuda.synchronize(); _THONG_KE["do_t_chuan_bi"] = _THONG_KE.get("do_t_chuan_bi", 0) + time.time() - _t0
    cur = float(lim)
    out_i = []; out_j = []; out_d = []
    for t0 in range(0, MA.shape[0], CT):
        if DO: torch.cuda.synchronize(); _a = time.time()
        g = torch.clamp(torch.maximum(blo[None] - thi[t0:t0 + CT, None], tlo[t0:t0 + CT, None] - bhi[None]), min=0)
        ti, bj = torch.nonzero(torch.linalg.norm(g, dim=2) <= cur + M32, as_tuple=True)
        del g
        ti = ti + t0
        if DO: torch.cuda.synchronize(); _THONG_KE["do_t_hop"] = _THONG_KE.get("do_t_hop", 0) + time.time() - _a; _THONG_KE["do_n_cap_nhom"] = _THONG_KE.get("do_n_cap_nhom", 0) + ti.numel()
        for p0 in range(0, ti.numel(), CP):
            if DO: torch.cuda.synchronize(); _a = time.time()
            tt = ti[p0:p0 + CP]; bb = bj[p0:p0 + CP]
            dm = torch.linalg.norm(MA[tt][:, :, None, :] - MB[bb][:, None, :, :], dim=-1)
            ok = (dm - LA[tt][:, :, None] - LB[bb][:, None, :]) <= cur + M32
            del dm
            p, i, j = torch.nonzero(ok, as_tuple=True)
            del ok
            gi = tt[p] * S + i; gj = bb[p] * S + j
            if DO: torch.cuda.synchronize(); _THONG_KE["do_t_loc"] = _THONG_KE.get("do_t_loc", 0) + time.time() - _a; _THONG_KE["do_n_pass"] = _THONG_KE.get("do_n_pass", 0) + gi.numel(); _a = time.time()
            for v in range(0, gi.numel(), 4_000_000):
                i_ = gi[v:v + 4_000_000]; j_ = gj[v:v + 4_000_000]
                d = _doan_doan_t(A0[i_], A1[i_], B0[j_], B1[j_])
                sel = d <= cur + m
                if sel.any():
                    ds = d[sel]
                    out_i.append(oA[i_[sel]]); out_j.append(oB[j_[sel]]); out_d.append(ds)
                    cur = min(cur, float(ds.min()))
            if DO: torch.cuda.synchronize(); _THONG_KE["do_t_seg"] = _THONG_KE.get("do_t_seg", 0) + time.time() - _a
    if not out_i:
        z = np.zeros(0, np.int64); return z, z
    I = torch.cat(out_i); J = torch.cat(out_j); D = torch.cat(out_d)
    k = D <= cur + m
    return I[k].cpu().numpy(), J[k].cpu().numpy()


def _lam_kc_mat_mat(K):
    goc = K.kc_mat_mat
    cKDTree = K.cKDTree

    def kc_mat_mat(A, B, R):
        dev = _gpu()
        if not dev:
            _THONG_KE["canh_canh_cpu"] += 1
            return goc(A, B, R)
        gA = np.linalg.norm(np.maximum(0.0, np.maximum(B.lo - A.V, A.V - B.hi)), axis=1)
        can = gA < R + A.Lmax / 2 + 1e-9
        dA = np.full(len(A.V), np.inf)
        if can.any():
            dA[can] = np.sqrt(B.kc2_diem(A.V[can]))
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
        if len(eA) and len(eB) and lim > 0:
            t0 = time.time()
            ia, jb = _cap_ung_vien_gpu(A.V, eA, B.V, eB, lim, dev)
            _THONG_KE["canh_canh_gpu"] += 1; _THONG_KE["t_canh_canh_gpu_s"] += time.time() - t0
            if len(ia):
                _THONG_KE["cap_tinh_lai"] += int(len(ia))
                a_ = eA[ia]; b_ = eB[jb]
                dd = K.kc_doan_doan(A.V[a_[:, 0]], A.V[a_[:, 1]], B.V[b_[:, 0]], B.V[b_[:, 1]])
                best = min(best, float(dd.min()))
        if best < R:
            return best, True
        return R, False

    kc_mat_mat.__doc__ = (goc.__doc__ or "") + "\n[tang_toc] buoc canh–canh tren GPU, tinh lai ung vien bang kc_doan_doan goc."
    kc_mat_mat._goc = goc
    return kc_mat_mat


class _Dem:
    def __init__(self):
        self.tep = None; self.d = {}; self.fo = None; self.dau = None

    def mo(self, thu_muc, dau):
        self.dong()
        self.tep = os.path.join(thu_muc, "kiem2_bo_dem.pkl"); self.d = {}; self.dau = dau
        hop_le = False
        if os.path.exists(self.tep):
            try:
                with open(self.tep, "rb") as f:
                    if pickle.load(f) == dau:
                        hop_le = True
                        while True:
                            try:
                                k, v = pickle.load(f); self.d[k] = v
                            except Exception:
                                break
            except Exception:
                hop_le = False
        if hop_le:
            tam = self.tep + ".tam"
            with open(tam, "wb") as f:
                pickle.dump(dau, f)
                for k, v in self.d.items(): pickle.dump((k, v), f)
            os.replace(tam, self.tep)
            self.fo = open(self.tep, "ab")
        else:
            self.d = {}; self.fo = open(self.tep, "wb"); pickle.dump(dau, self.fo); self.fo.flush()
        return len(self.d)

    def ghi(self, k, v):
        self.d[k] = v
        if self.fo is not None:
            pickle.dump((k, v), self.fo); self.fo.flush(); _THONG_KE["dem_ghi"] += 1

    def dong(self):
        if self.fo is not None:
            self.fo.close()
        self.fo = None; self.tep = None; self.d = {}


DEM = _Dem()


def _bam_luoi(L):
    h = getattr(L, "_bam", None)
    if h is None:
        b = hashlib.blake2b(digest_size=16)
        b.update(np.ascontiguousarray(L.V).tobytes()); b.update(np.ascontiguousarray(L.F).tobytes())
        h = b.hexdigest(); L._bam = h
    return h


def _khoa_bp(b):
    if b is None:
        return None
    if hasattr(b, "luoi") and hasattr(b, "ten"):
        return ("bp", b.ten, b.loai, _bam_luoi(b.luoi), b.man is None)
    if hasattr(b, "V") and hasattr(b, "F"):
        return ("luoi", _bam_luoi(b))
    if b.__class__.__name__ == "Manifold":
        return ("man", float(b.volume()), int(b.num_tri()), int(b.num_vert()))
    return ("x", repr(b))


def _ghi_nho(ten, fn, nguong_s=0.2):
    @functools.wraps(fn)
    def g(*a, **kw):
        if DEM.fo is None:
            return fn(*a, **kw)
        try:
            k = (ten,) + tuple(_khoa_bp(x) if not isinstance(x, (int, float, str)) else x for x in a) + tuple(sorted(kw.items()))
        except Exception:
            return fn(*a, **kw)
        if k in DEM.d:
            _THONG_KE["dem_trung"] += 1
            return DEM.d[k]
        t0 = time.time(); v = fn(*a, **kw)
        if time.time() - t0 >= nguong_s:
            try:
                pickle.dumps(v); DEM.ghi(k, v)
            except Exception:
                pass
        return v
    g._goc = fn
    return g


class _DangChay:
    def __init__(self):
        self.tep = None; self.ds = []; self.t0 = None

    def mo(self, thu_muc, so_dem):
        self.tep = os.path.join(thu_muc, "kiem2_dang_chay.json"); self.ds = []; self.t0 = time.time()
        self.so_dem = so_dem; self._ghi(False)

    def them(self, r):
        if self.tep is None:
            return
        x = dict(r); x["_t_s"] = round(time.time() - self.t0, 2)
        self.ds.append(x); self._ghi(False)

    def _ghi(self, xong):
        if self.tep is None:
            return
        tam = self.tep + ".tam"
        with open(tam, "w") as f:
            json.dump({"xong": xong, "ghi": "ket qua TUNG LUAT theo thu tu chay, TRUOC hau xu ly lop (KHONG_AP_DUNG) — ban cuoi o kiem2.json",
                       "so_ban_ghi_dem_nap_lai": self.so_dem, "thong_ke_tang_toc": dict(_THONG_KE), "luat": self.ds},
                      f, ensure_ascii=False, indent=1, default=_JSON_DEFAULT)
        os.replace(tam, self.tep)

    def dong(self, xong=True):
        self._ghi(xong); self.tep = None


DANG = _DangChay()


def _dau_van_tay(thu_muc):
    b = hashlib.blake2b(digest_size=16)
    for f in ("kiem2.py", "kiem_thucvat.py", "tang_toc.py"):
        p = os.path.join(os.path.dirname(os.path.abspath(__file__)), f)
        if os.path.exists(p):
            b.update(open(p, "rb").read())
    for f in ("mesh.npz", "loai.json"):
        p = os.path.join(thu_muc, f)
        if os.path.exists(p):
            st = os.stat(p); b.update(("%s|%d|%d" % (f, st.st_size, st.st_mtime_ns)).encode())
    b.update(PHIEN_BAN.encode()); b.update(os.environ.get("KIEM_GPU", "1").encode())
    return b.hexdigest()


def _boc_kiem(fn):
    @functools.wraps(fn)
    def g(thu_muc, *a, **kw):
        n = DEM.mo(thu_muc, _dau_van_tay(thu_muc))
        DANG.mo(thu_muc, n)
        try:
            out = fn(thu_muc, *a, **kw)
            DANG.dong(True)
            return out
        finally:
            DEM.dong()
    g._goc = fn
    return g


def _boc_ket_qua(fn):
    @functools.wraps(fn)
    def g(*a, **kw):
        r = fn(*a, **kw)
        try:
            DANG.them(r)
        except Exception:
            pass
        return r
    g._goc = fn
    return g


def _boc_the_tich(K, kieu):
    """the_tich_giao / the_tich_hieu: giữ NGUYÊN tác dụng phụ (khoi_cua → lam_day) rồi mới đệm phép bool."""
    def f(A, B):
        ma, mb = K.khoi_cua(A), K.khoi_cua(B)
        if ma is not None and mb is not None:
            k = ("bool_" + kieu, _khoa_bp(A), _khoa_bp(B))
            if DEM.fo is not None and k in DEM.d:
                _THONG_KE["dem_trung"] += 1; v = DEM.d[k]
            else:
                t0 = time.time(); v = float((ma ^ mb).volume()) if kieu == "giao" else float((ma - mb).volume())
                if DEM.fo is not None and time.time() - t0 >= 0.2: DEM.ghi(k, v)
            return v, "bool_chinh_xac"
        return K.du_phong_the_tich(A, B, kieu), "du_phong_mau_diem"
    return f


def _boc_du_phong(K):
    goc = K.du_phong_the_tich
    nho = _ghi_nho("du_phong_the_tich", goc)
    def f(A, B, kieu, h=0.02):
        lo, hi = A.lo.copy(), A.hi.copy()
        if kieu == "giao":
            lo = np.maximum(lo, B.lo); hi = np.minimum(hi, B.hi)
            if np.any(hi <= lo): return 0.0
        K.mat_chua(A); K.mat_chua(B)
        return nho(A, B, kieu, h)
    f._goc = goc
    return f


def cai(K):
    """K = module kiem2 (hoặc __main__ khi chạy kiem2.py trực tiếp)."""
    if getattr(K, "_tang_toc_da_cai", False):
        return
    global _JSON_DEFAULT
    _JSON_DEFAULT = K._json_default
    K.Luoi.__init__ = _luoi_init_nhanh
    if os.environ.get("KIEM_FWN_CACHE", "1") != "0":
        K.Luoi.quan = _quan_nhanh
    K.kc_mat_mat = _ghi_nho("kc_mat_mat", _lam_kc_mat_mat(K))
    K.du_phong_the_tich = _boc_du_phong(K)
    K.the_tich_giao = _boc_the_tich(K, "giao")
    K.the_tich_hieu = _boc_the_tich(K, "hieu")
    K.ket_qua = _boc_ket_qua(K.ket_qua)
    K.kiem = _boc_kiem(K.kiem)
    K._tang_toc_da_cai = True


def cai_tv(T):
    """T = module kiem_thucvat. Chỉ đệm hàm THUẦN (không tác dụng phụ lên bộ phận)."""
    if getattr(T, "_tang_toc_da_cai", False):
        return
    import kiem2 as K
    T.kc_mat_mat = K.kc_mat_mat; T.ket_qua = K.ket_qua; T.Luoi = K.Luoi
    goc_giao = T.giao
    def giao(A, B):
        if A.man is not None and B.man is not None:
            k = ("tv_giao", _khoa_bp(A), _khoa_bp(B))
            if DEM.fo is not None and k in DEM.d:
                _THONG_KE["dem_trung"] += 1; return DEM.d[k]
            t0 = time.time(); v = float((A.man ^ B.man).volume())
            if DEM.fo is not None and time.time() - t0 >= 0.2: DEM.ghi(k, v)
            return v
        return K.du_phong_the_tich(A, B, "giao")
    giao._goc = goc_giao
    T.giao = giao
    for ten in ("kc", "trong", "dien_tich_cham"):
        setattr(T, ten, _ghi_nho("tv_" + ten, getattr(T, ten)))
    T.kiem_tv = _boc_kiem(T.kiem_tv)
    T._tang_toc_da_cai = True
