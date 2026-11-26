"""gpu_tang_toc.py — tăng tốc cellPACK (bản sao cellpack_gpu) mà KHÔNG đổi kết quả.

Nguyên tắc (quyết định PI 15/11/2026): kết quả đóng gói phải TRÙNG HỆT bản gốc a576c5c (cùng bộ phận đặt, cùng vị trí,
cùng xoay, cùng seed). Vì vậy:
  * mọi vị từ quyết định (va chạm, trong/ngoài, lọc điểm tự do) giữ ĐÚNG công thức số thực float64 của bản gốc
    (cùng phép toán, cùng thứ tự cộng: ((dx²+dy²)+dz²), sqrt làm tròn đúng IEEE; không dùng FMA/rút gọn);
  * GPU chỉ được (a) tính cùng công thức trên tập điểm lớn hơn rồi chứng minh được phần thêm không đổi kết quả,
    hoặc (b) bỏ bớt ứng viên mà ta chứng minh được không thể làm vị từ đúng;
  * không đổi thứ tự/số lần gọi bộ sinh số ngẫu nhiên; thứ tự chèn của dict (ảnh hưởng thứ tự xoá điểm tự do) giữ nguyên.
Biến môi trường:
  CELLPACK_GPU=0        tắt mọi đường tăng tốc (chạy đúng mã gốc trong bản sao);
  CELLPACK_GPU_THIET_BI=cpu|cuda  đường va chạm dùng numpy (cpu) hay torch (cuda, mặc định nếu có);
  CELLPACK_GPU_KIEM=1   mỗi lần gọi chạy THÊM mã gốc và assert trùng (chậm; dùng khi thử tương đương);
  CELLPACK_TORCH_SITE   thư mục site-packages chứa torch (nối vào CUỐI sys.path, không che gói của env cellpack).
"""
import os, sys
import numpy

BAT = os.environ.get("CELLPACK_GPU", "1") != "0"
KIEM = os.environ.get("CELLPACK_GPU_KIEM", "0") == "1"
BAT_CHUA = BAT and os.environ.get("CELLPACK_GPU_CHUA", "1") != "0"
_TB = os.environ.get("CELLPACK_GPU_THIET_BI", "cuda")
THONG_KE = {"va_cham_goi": 0, "va_cham_co": 0, "va_cham_cpu_du_phong": 0, "kiem_ok": 0, "loc_tu_do": 0}
_torch = None


def torch_mod():
    global _torch
    if _torch is None:
        p = os.environ.get("CELLPACK_TORCH_SITE")
        if p and p not in sys.path:
            sys.path.append(p)
        import torch
        _torch = torch
    return _torch


def dung_gpu():
    if not BAT or _TB != "cuda":
        return False
    try:
        return torch_mod().cuda.is_available()
    except Exception:
        return False


class _LuoiGPU:
    """Bản sao trên GPU của toạ độ lưới (đổi float32 -> float64 là phép đổi chính xác) + bảng khoảng cách (cập nhật cùng lúc
    với mảng numpy). masterGridPositions = N0 = nx·ny·nz điểm lưới đều (thứ tự (j·nx + i)·nz + k, kiểm từng bit) rồi các
    điểm mặt ngăn nối thêm (extendGridArrays) — phần nối thêm được kiểm toàn bộ (vét cạn) trên GPU."""

    def __init__(self, grid, distances):
        t = torch_mod()
        self.dev = t.device("cuda")
        P = numpy.asarray(grid.masterGridPositions)
        assert P.dtype in (numpy.float32, numpy.float64) and P.ndim == 2 and P.shape[1] == 3
        self.grid = grid
        self.P_np = grid.masterGridPositions
        x, y, z = (numpy.asarray(v) for v in (grid._x, grid._y, grid._z))
        self.nx, self.ny, self.nz = len(x), len(y), len(z)
        self.N0 = N0 = self.nx * self.ny * self.nz
        self.mang_luoi_dung = False
        if N0 <= len(P) and min(self.nx, self.ny, self.nz) >= 2:
            J, I, K = numpy.meshgrid(numpy.arange(self.ny), numpy.arange(self.nx), numpy.arange(self.nz), indexing="ij")
            idx = ((J * self.nx + I) * self.nz + K).ravel()
            Q = numpy.empty((N0, 3), P.dtype)
            Q[idx, 0] = x[I.ravel()]; Q[idx, 1] = y[J.ravel()]; Q[idx, 2] = z[K.ravel()]
            self.mang_luoi_dung = bool(numpy.array_equal(Q, P[:N0]))
        self.x0, self.y0, self.z0 = float(x[0]), float(y[0]), float(z[0])
        self.s = (float(x[1]) - float(x[0]), float(y[1]) - float(y[0]), float(z[1]) - float(z[0]))
        self.P = t.from_numpy(numpy.ascontiguousarray(P, dtype=numpy.float64)).to(self.dev)
        self.nE = len(P) - N0
        self.E_idx = t.arange(N0, len(P), device=self.dev)
        self.dist_np = distances
        assert isinstance(distances, numpy.ndarray) and distances.dtype == numpy.float64 and len(distances) == len(P)
        self.D = t.from_numpy(numpy.ascontiguousarray(distances)).to(self.dev)

    def cap_nhat(self, idx, vals):
        """Ghi cùng giá trị đã ghi vào mảng numpy (copy từng bit)."""
        if len(idx) == 0:
            return
        t = torch_mod()
        self.D[t.as_tensor(numpy.asarray(idx, dtype=numpy.int64), device=self.dev)] = t.as_tensor(numpy.asarray(vals, dtype=numpy.float64), device=self.dev)

    def kiem_khop(self):
        t = torch_mod()
        a = self.D.cpu().numpy()
        return numpy.array_equal(a.view(numpy.uint64), numpy.ascontiguousarray(self.dist_np, dtype=numpy.float64).view(numpy.uint64))


_LUOI = {"obj": None}


def lay_luoi(grid, distances):
    """Trả bản GPU của lưới cho đúng mảng khoảng cách `distances` đang dùng trong pack_grid (đối tượng numpy, so bằng `is`)."""
    L = _LUOI["obj"]
    if L is None or L.dist_np is not distances or L.P_np is not grid.masterGridPositions:
        if not isinstance(distances, numpy.ndarray) or distances.dtype != numpy.float64:
            return None
        L = _LUOI["obj"] = _LuoiGPU(grid, distances)
        if KIEM:
            assert L.kiem_khop()
    return L


def bao_cap_nhat(distances, idx, vals):
    """Gọi từ BaseGrid.updateDistances sau khi ghi distance[pt] = dist (giữ bản GPU khớp từng bit)."""
    L = _LUOI["obj"]
    if L is not None and L.dist_np is distances:
        L.cap_nhat(idx, vals)
    bao_min(distances, vals)


def co_va_cham_gpu(L, centers, radii, dpad):
    """Vị từ va chạm của MultiSphereIngr.collision_jitter ở tầng sâu nhất:
        ∃ cầu s, ∃ điểm lưới g trong cKDTree.query_ball_point(c_s, r_s + dpad):  dist[g] + |P[g] − c_s| ≤ r_s
    Tính trên GPU trên hộp lưới bao quả cầu bán kính ρ_s = r_s − m, m = min(dist).
    Chứng minh trùng: vị từ đúng ⇒ |P[g] − c_s| ≤ r_s − dist[g] ≤ ρ_s; nếu ρ_s < r_s + dpad − biên (tức −m < dpad − biên)
    thì g chắc chắn nằm trong quả cầu của cKDTree, nên điểm thêm ngoài cầu cKDTree không bao giờ làm vị từ đúng; điểm trong
    cầu cKDTree mà ngoài hộp có |P[g] − c_s| > ρ_s nên vị từ sai. Vị từ trên từng điểm tính bằng cùng công thức float64.
    Trả None nếu điều kiện chứng minh không thoả (khi đó gọi đường CPU gốc)."""
    t = torch_mod()
    m = float(t.min(L.D).item())
    bien = 1e-6 * (1.0 + abs(dpad))
    if not L.mang_luoi_dung or not (-m < dpad - bien):
        return None
    if not all(isinstance(c, numpy.ndarray) and c.dtype == numpy.float64 for c in centers):
        return None
    C = numpy.asarray(centers, dtype=numpy.float64).reshape(-1, 3)
    R = numpy.asarray(radii, dtype=numpy.float64).reshape(-1)
    if len(R) == 0:
        return False
    rho = R - m
    o = numpy.array([L.x0, L.y0, L.z0]); s = numpy.array(L.s); n = numpy.array([L.nx, L.ny, L.nz])
    lo = numpy.floor((C - rho[:, None] - o) / s).astype(numpy.int64) - 1
    hi = numpy.ceil((C + rho[:, None] - o) / s).astype(numpy.int64) + 1
    lo = numpy.clip(lo, 0, n - 1); hi = numpy.clip(hi, 0, n - 1)
    w = hi - lo + 1
    if (w <= 0).any():
        w = numpy.maximum(w, 0)
    B = int(w.max()) if len(w) else 0
    if B == 0:
        return False
    dev = L.dev
    ar = t.arange(B, device=dev)
    lo_t = t.as_tensor(lo, device=dev); w_t = t.as_tensor(w, device=dev)
    ii = lo_t[:, 0:1] + ar[None, :]; jj = lo_t[:, 1:2] + ar[None, :]; kk = lo_t[:, 2:3] + ar[None, :]
    vi = ar[None, :] < w_t[:, 0:1]; vj = ar[None, :] < w_t[:, 1:2]; vk = ar[None, :] < w_t[:, 2:3]
    ii = t.where(vi, ii, lo_t[:, 0:1]); jj = t.where(vj, jj, lo_t[:, 1:2]); kk = t.where(vk, kk, lo_t[:, 2:3])
    idx = (jj[:, None, :, None] * L.nx + ii[:, :, None, None]) * L.nz + kk[:, None, None, :]
    hop = vi[:, :, None, None] & vj[:, None, :, None] & vk[:, None, None, :]
    idx = idx.reshape(len(R), -1); hop = hop.reshape(len(R), -1)
    Ct = t.as_tensor(C, device=dev); Rt = t.as_tensor(R, device=dev)
    Pg = L.P[idx]
    dx = Pg[..., 0] - Ct[:, 0:1]; dy = Pg[..., 1] - Ct[:, 1:2]; dz = Pg[..., 2] - Ct[:, 2:3]
    d2 = dx * dx; d2 = d2 + dy * dy; d2 = d2 + dz * dz
    dist = t.sqrt(d2)
    vc = ((L.D[idx] + dist) <= Rt[:, None]) & hop
    co = vc.any()
    if L.nE > 0:
        Pe = L.P[L.E_idx]
        ex = Pe[None, :, 0] - Ct[:, 0:1]; ey = Pe[None, :, 1] - Ct[:, 1:2]; ez = Pe[None, :, 2] - Ct[:, 2:3]
        e2 = ex * ex; e2 = e2 + ey * ey; e2 = e2 + ez * ez
        co = co | ((L.D[L.E_idx][None, :] + t.sqrt(e2)) <= Rt[:, None]).any()
    return bool(co.item())


def co_va_cham_cpu(grid, centers, radii, dpad, distances, gridPointsCoords):
    """Cùng vị từ, numpy trên CPU, đúng tập điểm cKDTree như bản gốc (không cần chứng minh tỉa)."""
    for r, c in zip(radii, centers):
        pts = grid.getPointsInSphere(c, r + dpad)
        delta = numpy.take(gridPointsCoords, pts, 0) - c
        delta *= delta
        distA = numpy.sqrt(delta.sum(1))
        if len(pts) and bool(((numpy.take(distances, pts) + distA) <= r).any()):
            return True
    return False


def diem_moi_va_trong(grid, centers, radii, dpad, distances, gridPointsCoords):
    """Khi KHÔNG va chạm: dựng (insidePoints, newDistPoints) đúng như vòng lặp gốc (giá trị + thứ tự chèn dict).
    Gốc, theo từng cầu s rồi từng điểm g (thứ tự cKDTree): d = |P[g]−c_s| − r_s;
      d ≤ 0: inside[g] = d nếu g chưa có hoặc |d| < |inside[g]|  (mọi d ≤ 0 ⇒ giữ d lớn nhất; hoà thì giá trị bằng nhau);
      ngược lại nếu d < distances[g]: newdist[g] = min(d, newdist[g]).
    Thứ tự chèn dict = lần xuất hiện đầu tiên của g trong dãy thoả điều kiện."""
    G = []; Dd = []
    for r, c in zip(radii, centers):
        pts = grid.getPointsInSphere(c, r + dpad)
        delta = numpy.take(gridPointsCoords, pts, 0) - c
        delta *= delta
        distA = numpy.sqrt(delta.sum(1))
        G.append(numpy.asarray(pts, dtype=numpy.int64)); Dd.append(distA - r)
    if not G:
        return {}, {}
    G = numpy.concatenate(G); Dd = numpy.concatenate(Dd)
    pos = numpy.arange(len(G))
    m_in = Dd <= 0
    m_nd = (~m_in) & (Dd < numpy.take(distances, G))
    return _gom(G[m_in], Dd[m_in], pos[m_in], lay_max=True), _gom(G[m_nd], Dd[m_nd], pos[m_nd], lay_max=False)


def _gom(g, d, p, lay_max):
    if len(g) == 0:
        return {}
    o = numpy.lexsort((p, g))
    g = g[o]; d = d[o]; p = p[o]
    dau = numpy.ones(len(g), bool); dau[1:] = g[1:] != g[:-1]
    st = numpy.nonzero(dau)[0]
    v = (numpy.maximum if lay_max else numpy.minimum).reduceat(d, st)
    thu_tu = numpy.argsort(p[st], kind="stable")
    keys = g[st][thu_tu].tolist(); vals = v[thu_tu]
    return {k: vals[i] for i, k in enumerate(keys)}


def loc_diem_tu_do(starting_array, array_length, distances, comp_ids, current_comp_id, cuttoff):
    """Thay vòng lặp Python của Ingredient.get_list_of_free_indices (nhánh update): giữ pt theo đúng thứ tự,
    điều kiện comp_ids[pt] == current_comp_id and distances[pt] >= cuttoff (cùng kiểu so sánh numpy)."""
    THONG_KE["loc_tu_do"] += 1
    if array_length <= 0:
        return []
    fp = numpy.asarray(starting_array[:array_length])
    m = (numpy.asarray(comp_ids)[fp] == current_comp_id) & (numpy.asarray(distances)[fp] >= cuttoff)
    return [starting_array[i] for i in numpy.flatnonzero(m).tolist()]


THONG_KE.update({"chua_goi": 0, "chua_diem": 0, "chua_mo_ho": 0, "chua_broken": 0, "chua_kiem_ok": 0})
_D0 = numpy.array([0.4395064455, 0.617598629942, 0.652231566745])
_LUOI_MESH = {}


def _mesh_gpu(mesh):
    t = torch_mod(); dev = t.device("cuda")
    k = id(mesh)
    M = _LUOI_MESH.get(k)
    if M is not None and M["mesh"] is mesh and M["nf"] == len(mesh.faces):
        return M
    tri = numpy.asarray(mesh.triangles, dtype=numpy.float64)
    N = numpy.asarray(mesh.face_normals, dtype=numpy.float64)
    d = _D0
    a = numpy.cross(d, [1.0, 0.0, 0.0]); a /= numpy.linalg.norm(a); b = numpy.cross(d, a)
    A2 = tri @ a; B2 = tri @ b
    e = numpy.linalg.norm(tri[:, [1, 2, 0]] - tri, axis=2)
    ar3 = 0.5 * numpy.linalg.norm(numpy.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
    emax = e.max(1); h3 = numpy.where(emax > 0, 2 * ar3 / numpy.maximum(emax, 1e-300), 0.0)
    den = N @ d
    nv0 = numpy.einsum("ij,ij->i", N, tri[:, 0])
    X = float(numpy.abs(tri).max())
    c2 = (A2[:, 1] - A2[:, 0]) * (B2[:, 2] - B2[:, 0]) - (B2[:, 1] - B2[:, 0]) * (A2[:, 2] - A2[:, 0])
    T = lambda x: t.as_tensor(numpy.ascontiguousarray(x), dtype=t.float64, device=dev)
    M = {"mesh": mesh, "nf": len(mesh.faces), "a": a, "b": b, "X": X,
         "A2": T(A2), "B2": T(B2), "den": T(den), "nv0": T(nv0), "N": T(N), "h3": T(h3), "emax": T(emax), "c2": T(c2)}
    if len(_LUOI_MESH) > 16:
        _LUOI_MESH.clear()
    _LUOI_MESH[k] = M
    return M


def _phan_loai_gpu(M, P):
    """P (n,3) float64 điểm đã trong AABB. Trả (chua, mo_ho) bool (n,)."""
    t = torch_mod(); dev = t.device("cuda")
    EPS = 2.220446049250313e-16
    n = len(P); F = M["nf"]
    q1 = P @ M["a"]; q2 = P @ M["b"]
    Xp = max(M["X"], float(numpy.abs(P).max()) if n else 0.0)
    den = M["den"]; aden = den.abs()
    con_hieu = aden > 1e-5 * (1 + 1e-9) + 1e-15
    chac_loai = aden < 1e-5 * (1 - 1e-9) - 1e-15
    out_c = numpy.zeros(n, bool); out_m = numpy.zeros(n, bool)
    cs = max(1, int(4_000_000 // max(F, 1)))
    A2, B2, c2 = M["A2"], M["B2"], M["c2"]
    sdc2 = t.where(c2 == 0, t.full_like(c2, float("nan")), c2)
    for i0 in range(0, n, cs):
        Q1 = t.as_tensor(q1[i0:i0 + cs], device=dev)[:, None]; Q2 = t.as_tensor(q2[i0:i0 + cs], device=dev)[:, None]
        Pt = t.as_tensor(P[i0:i0 + cs], device=dev)
        l0 = ((A2[:, 1] - Q1) * (B2[:, 2] - Q2) - (B2[:, 1] - Q2) * (A2[:, 2] - Q1)) / sdc2
        l1 = ((A2[:, 2] - Q1) * (B2[:, 0] - Q2) - (B2[:, 2] - Q2) * (A2[:, 0] - Q1)) / sdc2
        l2 = 1.0 - l0 - l1
        npd = Pt @ M["N"].T
        tt = (M["nv0"][None, :] - npd) / t.where(den == 0, t.ones_like(den), den)[None, :]
        Xt = Xp + tt.abs()
        mb = 1e3 * EPS * (Xt / (aden * M["h3"])[None, :] + (M["emax"] / M["h3"])[None, :] ** 2 + (Xt * Xt) / c2.abs()[None, :]) + 1e-9
        mt = 1e3 * EPS * Xt / aden[None, :] + 1e-12
        lmin = t.minimum(t.minimum(l0, l1), l2)
        chac_trong = lmin > mb
        chac_ngoai = lmin < -mb
        xuoi_c = tt > -1e-6 + mt; xuoi_k = tt < -1e-6 - mt
        nguoc_c = tt < 1e-6 - mt; nguoc_k = tt > 1e-6 + mt
        loai_hinh = chac_loai[None, :] | chac_ngoai
        trung_hinh = con_hieu[None, :] & chac_trong
        mo_hinh = ~(loai_hinh | trung_hinh)
        hf = trung_hinh & xuoi_c; hb = trung_hinh & nguoc_c
        mo = (mo_hinh & ~(xuoi_k & nguoc_k)) | (trung_hinh & ~(xuoi_c | xuoi_k)) | (trung_hinh & ~(nguoc_c | nguoc_k))
        mo_p = mo.any(1)
        nf = hf.sum(1); nb = hb.sum(1)
        for H, nH in ((hf, nf), (hb, nb)):
            k = int(nH.max().item()) if nH.numel() else 0
            if k >= 2:
                ts = t.where(H, tt, t.full_like(tt, float("inf"))).topk(k, dim=1, largest=False).values
                g = ts[:, 1:] - ts[:, :-1]
                mo_p = mo_p | ((g < 1e-6) & t.isfinite(g)).any(1)
        cf = (nf % 2) == 1; cb = (nb % 2) == 1
        dong_y = cf == cb
        broken = (~dong_y) & (nf > 0) & (nb > 0)
        mo_p = mo_p | broken
        THONG_KE["chua_broken"] += int(broken.sum().item())
        chua = dong_y & cf
        out_c[i0:i0 + cs] = chua.cpu().numpy(); out_m[i0:i0 + cs] = mo_p.cpu().numpy()
    return out_c, out_m


def chua_diem(mesh, points):
    """Thay mesh.contains(points) (trimesh, không embree): cùng kết quả từng điểm (xem chú thích trên)."""
    from trimesh import bounds as _tb
    THONG_KE["chua_goi"] += 1
    P = numpy.asanyarray(points, dtype=numpy.float64)
    if P.ndim != 2 or P.shape[1] != 3:
        return mesh.contains(points)
    kq = numpy.zeros(len(P), dtype=bool)
    trong_hop = _tb.contains(mesh.bounds, P)
    idx = numpy.flatnonzero(trong_hop)
    THONG_KE["chua_diem"] += len(idx)
    if len(idx) == 0:
        return kq
    c, m = _phan_loai_gpu(_mesh_gpu(mesh), P[idx])
    kq[idx[~m]] = c[~m]
    if m.any():
        THONG_KE["chua_mo_ho"] += int(m.sum())
        kq[idx[m]] = mesh.contains(P[idx[m]])
    if KIEM:
        g = mesh.contains(points)
        assert numpy.array_equal(g, kq), ("KIEM contains", int((g != kq).sum()), len(kq))
        THONG_KE["chua_kiem_ok"] += 1
    return kq


V3 = BAT and os.environ.get("CELLPACK_V3", "1") != "0"
THONG_KE.update({"v3_va_cham": 0, "v3_va_cham_du_phong": 0, "v3_chua_1": 0, "v3_chua_mo_ho": 0, "v3_kiem_ok": 0})
_MIN = {"arr": None, "m": None}


def min_khoang_cach(distances):
    """Cận dưới của min(distances), giữ bảo thủ: chỉ hạ xuống khi updateDistances ghi giá trị nhỏ hơn (không bao giờ nâng)."""
    if _MIN["arr"] is not distances:
        _MIN["arr"] = distances; _MIN["m"] = float(numpy.min(distances))
    return _MIN["m"]


def bao_min(distances, vals):
    if _MIN["arr"] is distances and len(vals):
        _MIN["m"] = min(_MIN["m"], float(numpy.min(vals)))


def co_va_cham_cat_tia(grid, centers, radii, dpad, distances, gridPointsCoords):
    """Cùng vị từ va chạm (xem co_va_cham_gpu); ứng viên = cKDTree.query_ball_point(c_s, ρ_s + biên), ρ_s = r_s − m.
    Vị từ đúng ⇒ |P[g] − c_s| ≤ r_s − dist[g] ≤ ρ_s; với −m < dpad − biên mọi điểm như vậy vừa ở trong cầu gốc (r_s + dpad)
    vừa ở trong cầu tỉa (ρ_s + biên); điểm trong cầu gốc mà ngoài cầu tỉa có khoảng cách > ρ_s nên vị từ sai. Vị từ trên từng
    điểm tính bằng đúng công thức numpy của bản gốc. Trả None nếu không thoả điều kiện chứng minh."""
    m = min_khoang_cach(distances)
    bien = 1e-6 * (1.0 + abs(dpad))
    if not (-m < dpad - bien):
        return None
    THONG_KE["v3_va_cham"] += 1
    if grid.tree is None:
        grid.getPointsInSphere(centers[0], 0.0)
    for r, c in zip(radii, centers):
        rho = r - m
        pts = grid.tree.query_ball_point(c, rho + 1e-6 * (1.0 + abs(rho)))
        if not pts:
            continue
        delta = numpy.take(gridPointsCoords, pts, 0) - c
        delta *= delta
        distA = numpy.sqrt(delta.sum(1))
        if bool(((numpy.take(distances, pts) + distA) <= r).any()):
            return True
    return False


_MESH_CPU = {}
_O_BIN = 0.25
_DELTA = 0.01


def _mesh_cpu(mesh):
    k = id(mesh)
    M = _MESH_CPU.get(k)
    if M is not None and M["mesh"] is mesh and M["nf"] == len(mesh.faces):
        return M
    EPS = 2.220446049250313e-16
    tri = numpy.asarray(mesh.triangles, dtype=numpy.float64)
    N = numpy.asarray(mesh.face_normals, dtype=numpy.float64)
    d = _D0
    a = numpy.cross(d, [1.0, 0.0, 0.0]); a /= numpy.linalg.norm(a); b = numpy.cross(d, a)
    A2 = tri @ a; B2 = tri @ b
    e = numpy.linalg.norm(tri[:, [1, 2, 0]] - tri, axis=2); emax = e.max(1)
    ar3 = 0.5 * numpy.linalg.norm(numpy.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
    h3 = numpy.where(emax > 0, 2 * ar3 / numpy.maximum(emax, 1e-300), 0.0)
    den = N @ d; aden = numpy.abs(den)
    nv0 = numpy.einsum("ij,ij->i", N, tri[:, 0])
    c2 = (A2[:, 1] - A2[:, 0]) * (B2[:, 2] - B2[:, 0]) - (B2[:, 1] - B2[:, 0]) * (A2[:, 2] - A2[:, 0])
    bmin = numpy.asarray(mesh.bounds[0], float); bmax = numpy.asarray(mesh.bounds[1], float)
    X = float(max(numpy.abs(tri).max(), numpy.abs(bmin).max(), numpy.abs(bmax).max()))
    Dg = float(numpy.linalg.norm(bmax - bmin)) + 1.0
    chac_loai = aden < 1e-5 * (1 - 1e-9) - 1e-15
    with numpy.errstate(divide="ignore", invalid="ignore"):
        Xt = X + Dg / aden
        m_max = 1e3 * EPS * (Xt / (aden * h3) + (emax / h3) ** 2) + 1e-9
        W = numpy.maximum(A2.max(1) - A2.min(1), B2.max(1) - B2.min(1))
        tia_duoc = (~chac_loai) & numpy.isfinite(m_max) & (aden > 1e-2) & (m_max < 0.5 * _DELTA / (3 * numpy.maximum(W, 1e-300)))
    luon = numpy.flatnonzero(~chac_loai & ~tia_duoc)
    tb = numpy.flatnonzero(tia_duoc)
    lo_a = A2[tb].min(1) - _DELTA; hi_a = A2[tb].max(1) + _DELTA; lo_b = B2[tb].min(1) - _DELTA; hi_b = B2[tb].max(1) + _DELTA
    oa = (lo_a.min() if len(tb) else 0.0) - _O_BIN; ob = (lo_b.min() if len(tb) else 0.0) - _O_BIN
    na = int(((hi_a.max() if len(tb) else 0.0) - oa) // _O_BIN) + 2; nb = int(((hi_b.max() if len(tb) else 0.0) - ob) // _O_BIN) + 2
    i0 = numpy.floor((lo_a - oa) / _O_BIN).astype(numpy.int64); i1 = numpy.floor((hi_a - oa) / _O_BIN).astype(numpy.int64)
    j0 = numpy.floor((lo_b - ob) / _O_BIN).astype(numpy.int64); j1 = numpy.floor((hi_b - ob) / _O_BIN).astype(numpy.int64)
    cnt = (i1 - i0 + 1) * (j1 - j0 + 1)
    owner = numpy.repeat(numpy.arange(len(tb)), cnt)
    off = numpy.arange(cnt.sum()) - numpy.repeat(numpy.cumsum(cnt) - cnt, cnt)
    wj = numpy.repeat(j1 - j0 + 1, cnt)
    bi = numpy.repeat(i0, cnt) + off // wj; bj = numpy.repeat(j0, cnt) + off % wj
    key = bi * nb + bj
    o = numpy.argsort(key, kind="stable"); key = key[o]; owner = tb[owner[o]]
    start = numpy.searchsorted(key, numpy.arange(na * nb)); stop = numpy.searchsorted(key, numpy.arange(na * nb), side="right")
    M = {"mesh": mesh, "nf": len(mesh.faces), "a": a, "b": b, "A2": A2, "B2": B2, "den": den, "aden": aden, "nv0": nv0, "N": N,
         "h3": h3, "emax": emax, "c2": c2, "X": X, "chac_loai": chac_loai, "luon": luon, "owner": owner, "start": start, "stop": stop,
         "oa": oa, "ob": ob, "na": na, "nb": nb, "so_tia": int(len(tb)), "so_luon": int(len(luon))}
    if len(_MESH_CPU) > 16:
        _MESH_CPU.clear()
    _MESH_CPU[k] = M
    return M


def _phan_loai_1(M, p):
    """Một điểm (đã trong AABB), numpy, cùng phép phân loại có biên như _phan_loai_gpu, trên ứng viên = ô 2D ∪ 'luôn xét'.
    Tam giác bị tỉa: hộp bao 2D (lề δ) không chứa hình chiếu q ⇒ toạ độ trọng tâm thật < −δ/(3W) và sai số của trimesh
    < m_max < δ/(6W) ⇒ trimesh chắc chắn trượt."""
    EPS = 2.220446049250313e-16
    q1 = float(p @ M["a"]); q2 = float(p @ M["b"])
    ia = int((q1 - M["oa"]) // _O_BIN); ib = int((q2 - M["ob"]) // _O_BIN)
    if 0 <= ia < M["na"] and 0 <= ib < M["nb"]:
        kk = ia * M["nb"] + ib
        cand = numpy.concatenate([M["owner"][M["start"][kk]:M["stop"][kk]], M["luon"]])
    else:
        cand = M["luon"]
    if len(cand) == 0:
        return False, False
    A2 = M["A2"][cand]; B2 = M["B2"][cand]; c2 = M["c2"][cand]; den = M["den"][cand]; aden = M["aden"][cand]
    with numpy.errstate(divide="ignore", invalid="ignore"):
        sdc2 = numpy.where(c2 == 0, numpy.nan, c2)
        l0 = ((A2[:, 1] - q1) * (B2[:, 2] - q2) - (B2[:, 1] - q2) * (A2[:, 2] - q1)) / sdc2
        l1 = ((A2[:, 2] - q1) * (B2[:, 0] - q2) - (B2[:, 2] - q2) * (A2[:, 0] - q1)) / sdc2
        l2 = 1.0 - l0 - l1
        npd = M["N"][cand] @ p
        tt = (M["nv0"][cand] - npd) / numpy.where(den == 0, 1.0, den)
        Xp = max(M["X"], float(numpy.abs(p).max()))
        Xt = Xp + numpy.abs(tt)
        mb = 1e3 * EPS * (Xt / (aden * M["h3"][cand]) + (M["emax"][cand] / M["h3"][cand]) ** 2 + (Xt * Xt) / numpy.abs(c2)) + 1e-9
        mt = 1e3 * EPS * Xt / aden + 1e-12
        lmin = numpy.minimum(numpy.minimum(l0, l1), l2)
        con_hieu = aden > 1e-5 * (1 + 1e-9) + 1e-15
        chac_trong = lmin > mb; chac_ngoai = lmin < -mb
        xuoi_c = tt > -1e-6 + mt; xuoi_k = tt < -1e-6 - mt
        nguoc_c = tt < 1e-6 - mt; nguoc_k = tt > 1e-6 + mt
    loai_hinh = M["chac_loai"][cand] | chac_ngoai
    trung_hinh = con_hieu & chac_trong
    mo_hinh = ~(loai_hinh | trung_hinh)
    hf = trung_hinh & xuoi_c; hb = trung_hinh & nguoc_c
    mo = (mo_hinh & ~(xuoi_k & nguoc_k)) | (trung_hinh & ~(xuoi_c | xuoi_k)) | (trung_hinh & ~(nguoc_c | nguoc_k))
    if mo.any():
        return False, True
    for H in (hf, hb):
        if H.sum() >= 2:
            ts = numpy.sort(tt[H])
            if (numpy.diff(ts) < 1e-6).any():
                return False, True
    nf = int(hf.sum()); nb_ = int(hb.sum())
    cf = nf % 2 == 1; cb = nb_ % 2 == 1
    if cf != cb and nf > 0 and nb_ > 0:
        return False, True
    return (cf and cf == cb), False


def chua_1_diem(mesh, point):
    """Thay mesh.contains([point])[0]."""
    from trimesh import bounds as _tb
    THONG_KE["v3_chua_1"] += 1
    P = numpy.asanyarray([point], dtype=numpy.float64)
    if P.shape != (1, 3):
        return mesh.contains([point])[0]
    trong = _tb.contains(mesh.bounds, P)
    if not trong[0]:
        kq = numpy.zeros(1, dtype=bool)
    else:
        c, mh = _phan_loai_1(_mesh_cpu(mesh), P[0])
        if mh:
            THONG_KE["v3_chua_mo_ho"] += 1
            kq = mesh.contains(P)
        else:
            kq = numpy.array([c], dtype=bool)
    if KIEM:
        g = mesh.contains([point])
        assert numpy.array_equal(g, kq), ("KIEM chua_1", g, kq)
        THONG_KE["v3_kiem_ok"] += 1
    return kq[0]
