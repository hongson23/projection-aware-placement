# -*- coding: utf-8 -*-
"""sua_mesh_dn.py (22/11/2026) — SỬA MESH ĐỒNG NHẤT, bước tiền xử lý thứ hai trước khi chấm (sau ve_sinh_luoi.py, trước cham_lo.sh).
Thiết kế: Inside3D/ThuNghiem/ThietKe_SuaMesh_DongNhat_20261122.md (lớp lỗi L2 của LoiMesh_va_TrichXuatChungNhan_20261121.md;
cùng phép với phân tích (b) đã rà hai vòng, bo_dung_chung/v122/hau_xu_ly_b.py bản 3).

Dùng:  python sua_mesh_dn.py <vao.npz> <ra.npz> --vat dong_vat|thuc_vat [--nhat_ky <ra.json>] [--nen lzma|deflate]

Phép (MỘT bảng cặp cố định, áp như nhau cho mọi nhánh / mọi phương pháp, không đọc nhãn nhánh):
  với mỗi cặp (chính A, phụ B) trong CAP[vat] mà cả hai có mặt:  B ← B ∖ A  (hiệu boolean manifold3d, double), A giữ nguyên.
    dong_vat: (nhan, luoi_noi_chat)  — cặp được chạm duy nhất (dv.mien_cap; M.c₀ + 4.4a)
    thuc_vat: (vach, mang)           — màng áp mặt trong thành (2.3a, 2.3c: màng không giao vật liệu thành); ở đây mặt ngoài màng và
                                       mặt trong thành TRÙNG nhau (phẳng), phần "giao" là phim suy biến do nhiễu ký hiệu (phụ thuộc thứ tự toán hạng)
  Chỉ sửa khi phần B lấn vào A ở MỨC HIỆN VẬT LƯỚI (giao khác rỗng VÀ độ sâu xuyên ≤ D_MAX = 0,5·h VÀ thể tích bỏ ≤ R_MAX·V(B)).
  Lấn lớn hơn là vi phạm thật của phương pháp ⇒ KHÔNG sửa (giữ nguyên để bộ kiểm báo), trạng thái "vuot_nguong_hien_vat".
Đỉnh mới là giao điểm double ⇒ B′ lưu đỉnh float64 (làm tròn float32 đưa lại ≈1e-6 µm xuyên); mặt giữ kiểu số nguyên gốc.
Mọi thành phần zip khác chép NGUYÊN BYTE, giữ thứ tự. Không đổi gì ⇒ <ra.npz> là liên kết cứng tới <vao.npz>.
Kiểm sau (trên đúng mảng sẽ ghi): Manifold NoError, kín + cùng chiều, số khối (V > 1e-9) và giống không đổi, cổng kiem2.BoPhan
không loi_chan / không tự cắt, |V(B′ ∩ A)| ≤ 1e-9 µm³ theo CẢ HAI thứ tự toán hạng, độ sâu xuyên sau (đo trên mặt khối giao,
cùng đại lượng kiem3.do_sau_xuyen) ≤ SAU_SAU_MAX[vat] (1e-5 µm động vật, 1e-4 µm thực vật), cục bộ: trọng tâm mọi mặt
B′ cách ∂A ∪ ∂B ≤ 1e-6 µm. Nếu đầu ra boolean có mặt diện tích ≤ 1e-14 µm² (bộ kiểm sẽ xoá ⇒ hở) thì TRƯỚC khi kiểm, gộp các
đỉnh nối bởi cạnh ≤ 2e-6 µm của chính các mặt đó (han_vi_mo) rồi bỏ mặt suy biến và các cặp mặt trùng ngược chiều sinh ra; các mặt
kề có đỉnh bị gộp dời theo (đo ≤ 2,4e-7 µm). Hàn rồi mà kiểm trượt ⇒ "loi_han:<ly_do>" (giữ nguyên B).
Trạng thái mỗi cặp: da_sua | khong_giao | giu_nguyen:<ly_do> | vuot_nguong_hien_vat | loi_han:<ly_do> | thieu_bo_phan; --bang ghi TSV. Trượt bất kỳ ⇒ giữ nguyên B (không sửa tay), trạng thái "giu_nguyen:<ly_do>".
Không sửa tệp vào; ghi nguyên tử. Nhật ký json: mỗi cặp mọi số đo trước/sau + giây.
Mã chỉ ĐỌC (import) bo_kiem_v3/kiem2.py đóng băng (BoPhan = cổng hợp lệ M.f của bộ kiểm), không gọi luật nào."""
import os, sys, io, json, time, zipfile, argparse, hashlib, shutil
import numpy as np
import manifold3d as mf
import trimesh
import igl

HERE = os.path.dirname(os.path.abspath(__file__))
KDIR = os.environ.get("BO_KIEM_V3", HERE)
if KDIR not in sys.path: sys.path.insert(0, KDIR)
CAP = {"dong_vat": [("nhan", "luoi_noi_chat")], "thuc_vat": [("vach", "mang")]}
H_VAT = {"dong_vat": 23.0 / 128.0, "thuc_vat": 0.05}
D_MAX_H = 0.5
R_MAX = 1e-4
EPS_GIAO_SAU = 1e-9
SAU_SAU_MAX = {"dong_vat": 1e-5,
               "thuc_vat": 1e-4}
CUC_BO_MAX = 1e-6


def man(V, F):
    m = mf.Manifold(mf.Mesh64(np.ascontiguousarray(V, np.float64), np.ascontiguousarray(F, np.uint64)))
    return m, str(m.status())


def ra_mesh(m):
    me = m.to_mesh64()
    return np.asarray(me.vert_properties, np.float64)[:, :3].copy(), np.asarray(me.tri_verts, np.int64).copy()


def mau(V, F, k=8):
    u, v = np.meshgrid(np.arange(k + 1), np.arange(k + 1), indexing="ij"); m = (u + v) <= k
    b = np.stack([u[m], v[m]], 1) / k; w = np.stack([1 - b.sum(1), b[:, 0], b[:, 1]], 1)
    return np.einsum("bj,fjd->fbd", w, V[F]).reshape(-1, 3)


def do_sau_giao(I, VA, FA, VB, FB):
    """độ sâu xuyên hai chiều đo trên MẶT của khối giao I = A ∩ B: ∂I = (∂A ∩ B) ∪ (∂B ∩ A), nên max trên điểm mẫu của ∂I
    (k = 8 mỗi tam giác, như kiem3.do_sau_xuyen) của max(d(·, ∂A), d(·, ∂B)) = cùng đại lượng với kiem3.do_sau_xuyen, nhưng chỉ
    lấy mẫu phần giao (rẻ với cặp vách–màng phủ cả tế bào). Cận dưới như hàm của bộ kiểm."""
    if I.is_empty():
        return 0.0
    VI, FI = ra_mesh(I)
    if len(FI) == 0:
        return 0.0
    P = mau(VI, FI)
    return float(np.maximum(kc(P, VA, FA), kc(P, VB, FB)).max())


def so_khoi(V, F):
    t = trimesh.Trimesh(np.asarray(V, float), np.asarray(F, np.int64), process=False)
    c = t.split(only_watertight=False)
    return int(sum(1 for x in c if abs(x.volume) > 1e-9)), bool(t.is_watertight), bool(t.is_winding_consistent)


def kc(P, V, F):
    out = []
    for i in range(0, len(P), 400000):
        d2, _, _ = igl.point_mesh_squared_distance(np.ascontiguousarray(P[i:i + 400000], np.float64), np.ascontiguousarray(V, np.float64), np.ascontiguousarray(F, np.int64))
        out.append(np.sqrt(np.maximum(d2, 0)))
    return np.concatenate(out) if out else np.zeros(0)


A_VI_MO = 1e-14
E_HAN = 2e-6


def bo_cap_trung_nguoc(F):
    """Gộp đỉnh có thể biến hai mặt kề thành CÙNG tập đỉnh, ngược chiều (một nếp gấp thể tích 0: cạnh dùng 4 lần ⇒ lưới không kín
    theo trimesh; đo ở cầu nối cellPACK thử D1/D1b/D4_9001, rà soát 2 K1). Bỏ cả hai mặt của mỗi cặp như vậy (thể tích, mặt hình
    học không đổi). Cặp trùng CÙNG chiều không bỏ (để kiểm sau bắt)."""
    S = np.sort(F, axis=1); key = S[:, 0] * (2 ** 42) + S[:, 1] * (2 ** 21) + S[:, 2] if F.max() < 2 ** 21 else None
    if key is None:
        _, inv, cnt = np.unique(S, axis=0, return_inverse=True, return_counts=True); inv = inv.ravel()
    else:
        _, inv, cnt = np.unique(key, return_inverse=True, return_counts=True)
    bo = np.zeros(len(F), bool)
    for g in np.nonzero(cnt == 2)[0]:
        i, j = np.nonzero(inv == g)[0]
        a_, b_ = list(F[i]), list(F[j]); k = b_.index(a_[0]); cung = b_[k:] + b_[:k] == a_
        if not cung: bo[i] = bo[j] = True
    return F[~bo]


def han_vi_mo(V, F):
    """Gộp các cụm đỉnh nối bởi cạnh ≤ E_HAN thuộc mặt diện tích ≤ A_VI_MO (đỉnh giữ lại = chỉ số nhỏ nhất của cụm, giữ toạ độ
    của nó), bỏ mặt suy biến + đỉnh thừa. Tối đa 3 vòng. Trả (V, F, số mặt vi mô trước, số đỉnh gộp). Kiểm hợp lệ làm ở bên gọi."""
    n0 = None; gop = 0
    for _ in range(3):
        T = V[F]; ar = 0.5 * np.linalg.norm(np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]), axis=1)
        vm = np.nonzero(ar <= A_VI_MO)[0]
        if n0 is None: n0 = int(len(vm))
        if not len(vm): break
        cha = np.arange(len(V))
        def goc(i):
            while cha[i] != i: cha[i] = cha[cha[i]]; i = cha[i]
            return i
        for f in vm:
            for i, j in ((0, 1), (1, 2), (2, 0)):
                a_, b_ = F[f, i], F[f, j]
                if np.linalg.norm(V[a_] - V[b_]) <= E_HAN:
                    ra_, rb = goc(a_), goc(b_)
                    if ra_ != rb: cha[max(ra_, rb)] = min(ra_, rb); gop += 1
        m = np.array([goc(i) for i in range(len(V))])
        F = m[F]; F = F[(F[:, 0] != F[:, 1]) & (F[:, 1] != F[:, 2]) & (F[:, 0] != F[:, 2])]
        F = bo_cap_trung_nguoc(F)
        dung = np.zeros(len(V), bool); dung[F.ravel()] = True; so = np.cumsum(dung) - 1
        V = V[dung]; F = so[F]
    return V, F, n0, gop


def sua_cap(VA, FA, VB, FB, h, sau_max):
    """trả (V', F' hoặc None, bản ghi)."""
    t0 = time.time(); r = {}
    mA, sA = man(VA, FA); mB, sB = man(VB, FB)
    r["manifold_vao"] = [sA, sB]
    if sA != "Error.NoError" or sB != "Error.NoError" or mA.is_empty() or mB.is_empty():
        r["trang_thai"] = "giu_nguyen:dau_vao_khong_dung_duoc_manifold"; r["giay"] = round(time.time() - t0, 2); return None, r
    I1, I2_ = mA ^ mB, mB ^ mA
    r["giao_truoc_um3"] = [float(I1.volume()), float(I2_.volume())]; r["V_B_truoc_um3"] = float(mB.volume()); r["giong_B"] = int(mB.genus())
    if I1.is_empty() and I2_.is_empty():
        r["trang_thai"] = "khong_giao"; r["giay"] = round(time.time() - t0, 2); return None, r
    r["sau_truoc_um"] = max(do_sau_giao(I1, VA, FA, VB, FB), do_sau_giao(I2_, VA, FA, VB, FB))
    D = mB - mA
    if str(D.status()) != "Error.NoError" or D.is_empty():
        r["trang_thai"] = "giu_nguyen:bool_loi_" + str(D.status()); r["giay"] = round(time.time() - t0, 2); return None, r
    r["V_bo_um3"] = r["V_B_truoc_um3"] - float(D.volume()); r["V_bo_tuong_doi"] = r["V_bo_um3"] / max(r["V_B_truoc_um3"], 1e-300)
    if r["sau_truoc_um"] > D_MAX_H * h or r["V_bo_tuong_doi"] > R_MAX:
        r["trang_thai"] = "vuot_nguong_hien_vat"; r["nguong"] = {"D_MAX_um": D_MAX_H * h, "R_MAX": R_MAX}
        r["giay"] = round(time.time() - t0, 2); return None, r
    V2, F2 = ra_mesh(D)
    V2, F2, nvm, ngop = han_vi_mo(V2, F2)
    r["mat_vi_mo_bool"] = nvm; r["dinh_han"] = ngop
    m2, s2 = man(V2, F2)
    kb, kin, cc = so_khoi(V2, F2); kb0, _, _ = so_khoi(VB, FB)
    r.update({"so_dinh": [int(len(VB)), int(len(V2))], "so_mat": [int(len(FB)), int(len(F2))], "manifold_sau": s2,
              "kin_cung_chieu": [kin, cc], "so_khoi": [kb0, kb]})
    ly = []
    if s2 != "Error.NoError": ly.append("manifold_sau")
    else:
        J1, J2 = m2 ^ mA, mA ^ m2
        r["giong_B_sau"] = int(m2.genus()); r["giao_sau_um3"] = [float(J1.volume()), float(J2.volume())]
        r["sau_sau_um"] = max(do_sau_giao(J1, VA, FA, V2, F2), do_sau_giao(J2, VA, FA, V2, F2))
        if r["giong_B_sau"] != r["giong_B"]: ly.append("giong_doi")
        if max(abs(x) for x in r["giao_sau_um3"]) > EPS_GIAO_SAU: ly.append("giao_sau")
        if not r["sau_sau_um"] <= sau_max: ly.append("sau_sau")
    if not (kin and cc): ly.append("khong_kin_hoac_nguoc_chieu")
    if kb != kb0: ly.append("so_khoi_doi")
    import kiem2 as K2
    b = K2.BoPhan("phu", "phu", V2, F2)
    r["bophan_sau"] = {"loi_chan": list(b.loi_chan), "loi": list(b.loi), "so_mat_tu_cat": int(b.so_mat_tu_cat)}
    if b.loi_chan or b.so_mat_tu_cat: ly.append("cong_BoPhan")
    C = V2[F2].mean(1); dmin = np.minimum(kc(C, VB, FB), kc(C, VA, FA))
    r["cuc_bo_max_um"] = float(dmin.max()) if len(dmin) else 0.0
    if not r["cuc_bo_max_um"] <= CUC_BO_MAX: ly.append("cuc_bo")
    r["giay"] = round(time.time() - t0, 2)
    if ly:
        r["trang_thai"] = ("loi_han:" if nvm else "giu_nguyen:") + ",".join(ly); return None, r
    r["trang_thai"] = "da_sua"
    return (V2, F2), r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vao"); ap.add_argument("ra"); ap.add_argument("--vat", required=True, choices=sorted(CAP))
    ap.add_argument("--nhat_ky", default=None); ap.add_argument("--nen", default="lzma", choices=["lzma", "deflate"])
    ap.add_argument("--bang", default=None, help="tệp TSV: nối MỘT dòng trạng thái mỗi cặp (luôn ghi, kể cả khong_giao / giu_nguyen)")
    a = ap.parse_args()
    vao, ra = os.path.abspath(a.vao), os.path.abspath(a.ra)
    if os.path.lexists(ra): sys.exit("ra da co: %s (khong ghi de)" % ra)
    t0 = time.time()
    zf = zipfile.ZipFile(vao); tho = [(i, zf.read(i.filename)) for i in zf.infolist()]
    ten_tp = [i.filename for i, _ in tho]; assert all(t.endswith(".npy") for t in ten_tp)
    khoa = [t[:-4] for t in ten_tp]
    z = np.load(vao, allow_pickle=False)
    nk = {"vao": vao, "ra": ra, "vat": a.vat, "md5_ma": hashlib.md5(open(os.path.abspath(__file__), "rb").read()).hexdigest(),
          "md5_vao": hashlib.md5(open(vao, "rb").read()).hexdigest(), "manifold3d": getattr(mf, "__version__", None),
          "tham_so": {"CAP": CAP[a.vat], "h": H_VAT[a.vat], "D_MAX_um": D_MAX_H * H_VAT[a.vat], "R_MAX": R_MAX,
                      "EPS_GIAO_SAU": EPS_GIAO_SAU, "SAU_SAU_MAX": SAU_SAU_MAX[a.vat], "CUC_BO_MAX": CUC_BO_MAX}, "cap": {}}
    moi = {}
    for A_, B_ in CAP[a.vat]:
        ten = A_ + "|" + B_
        if not all(k in khoa for k in (A_ + "__v", A_ + "__f", B_ + "__v", B_ + "__f")):
            nk["cap"][ten] = {"trang_thai": "thieu_bo_phan"}; continue
        FBo = z[B_ + "__f"]
        try:
            kq, r = sua_cap(z[A_ + "__v"].astype(np.float64), z[A_ + "__f"].astype(np.int64),
                            z[B_ + "__v"].astype(np.float64), FBo.astype(np.int64), H_VAT[a.vat], SAU_SAU_MAX[a.vat])
        except Exception as ex:
            import traceback
            kq, r = None, {"trang_thai": "giu_nguyen:ngoai_le", "loi": traceback.format_exc()[-1500:]}
        nk["cap"][ten] = r
        if kq is not None:
            V2, F2 = kq
            if F2.max() > np.iinfo(FBo.dtype).max: nk["cap"][ten]["trang_thai"] = "giu_nguyen:tran_kieu_mat"; continue
            moi[B_ + "__v"] = np.ascontiguousarray(V2, np.float64); moi[B_ + "__f"] = np.ascontiguousarray(F2.astype(FBo.dtype))
    z.close()
    tam = os.path.join(os.path.dirname(ra), ".tam_%d_%s" % (os.getpid(), os.path.basename(ra)))
    if not moi:
        try: os.link(vao, tam); nk["cach_ghi"] = "lien_ket_cung (khong doi gi)"
        except OSError: shutil.copy2(vao, tam); nk["cach_ghi"] = "chep (khong doi gi)"
    else:
        kieu = zipfile.ZIP_LZMA if a.nen == "lzma" else zipfile.ZIP_DEFLATED
        with zipfile.ZipFile(tam, "w", kieu, allowZip64=True) as o:
            for (i, data), k in zip(tho, khoa):
                if k in moi:
                    buf = io.BytesIO(); np.lib.format.write_array(buf, moi[k], allow_pickle=False); o.writestr(i.filename, buf.getvalue())
                else:
                    o.writestr(i.filename, data)
        nk["cach_ghi"] = "npz moi (%s)" % a.nen
        z2 = zipfile.ZipFile(tam); assert [i.filename for i in z2.infolist()] == ten_tp
        for (i, data), k in zip(tho, khoa):
            d2 = z2.read(i.filename)
            if k in moi:
                x = np.load(io.BytesIO(d2), allow_pickle=False); assert x.dtype == moi[k].dtype and x.tobytes() == moi[k].tobytes(), k
            else:
                assert d2 == data, k
    os.replace(tam, ra)
    nk["md5_ra"] = hashlib.md5(open(ra, "rb").read()).hexdigest(); nk["giay"] = round(time.time() - t0, 2)
    if a.nhat_ky: json.dump(nk, open(a.nhat_ky, "w"), ensure_ascii=False, indent=1, default=float)
    if a.bang:
        with open(a.bang, "a") as f_:
            for ten, v in nk["cap"].items():
                tt = v.get("trang_thai", "")
                lop = tt.split(":")[0]
                g = v.get("giao_truoc_um3") or [None, None]
                f_.write("\t".join(str(x) for x in (os.path.dirname(ra), a.vat, ten, lop, tt, v.get("sau_truoc_um"), v.get("V_bo_um3"),
                                                       v.get("mat_vi_mo_bool"), v.get("dinh_han"), v.get("giay"), g[0], g[1])) + "\n")
    print("%s: %s, %s, %.1fs" % (os.path.basename(os.path.dirname(ra)), {k: v.get("trang_thai") for k, v in nk["cap"].items()}, nk["cach_ghi"], nk["giay"]))


if __name__ == "__main__":
    main()
