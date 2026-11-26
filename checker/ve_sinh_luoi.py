# -*- coding: utf-8 -*-
"""ve_sinh_luoi.py — VỆ SINH LƯỚI, bước tiền xử lý trước khi chấm (PI 16/11 §16 điểm 1; dựa trên chan_doan_cd/a8_sua_luoi.py).

Dùng:  python ve_sinh_luoi.py <vao.npz> <ra.npz> [--nhat_ky <ra.json>] [--nen lzma|deflate]

Với MỖI bộ phận <ten> (cặp mảng <ten>__v, <ten>__f trong npz), đúng ba phép, theo thứ tự:
  1) gộp các đỉnh TRÙNG HỆT toạ độ, so trên giá trị đã lưu (float32; dung sai 0; chỉ đồng nhất -0.0 với +0.0);
     đỉnh giữ lại là đỉnh xuất hiện ĐẦU TIÊN, giữ nguyên thứ tự đỉnh;
  2) bỏ các mặt có chỉ số đỉnh lặp (sau phép 1), giữ nguyên thứ tự các mặt còn lại và thứ tự đỉnh trong mỗi mặt;
  3) bỏ các đỉnh không còn mặt nào dùng, giữ nguyên thứ tự.
KHÔNG làm gì khác: không bỏ mặt trùng / mặt lưng-đối-lưng, không lật mặt, không tách thành phần, không đổi dtype.
Kiểm: thể tích có dấu (float64) của từng bộ phận đổi tương đối ≤ 1e-9, ngược lại thoát mã 2 và KHÔNG ghi đầu ra.
Mọi thay đổi ghi vào nhật ký json (mỗi bộ phận: số đỉnh gộp, số mặt bỏ, số đỉnh thừa bỏ, thể tích trước/sau,
số cạnh biên / cạnh dùng > 2 mặt / mặt trùng trước và sau).
Bộ phận không đổi: các mảng của nó được chép NGUYÊN BYTE (đọc thô từ zip). Nếu KHÔNG bộ phận nào đổi: <ra.npz> là liên kết
cứng tới <vao.npz> (khác ổ đĩa thì chép), tức trùng từng byte. Có thay đổi: ghi npz mới (zip, nén lzma mặc định — np.load đọc
được như thường), các thành phần zip theo đúng thứ tự gốc.
Không sửa tệp vào. Ghi nguyên tử (tệp tạm + os.replace)."""
import os, sys, json, io, zipfile, argparse, hashlib, shutil
import numpy as np

NGUONG_DV = 1e-9


def the_tich(V, F):
    if len(F) == 0:
        return 0.0
    P = np.asarray(V, np.float64)[np.asarray(F, np.int64)]
    return float(np.einsum("ij,ij->i", P[:, 0], np.cross(P[:, 1], P[:, 2])).sum() / 6.0)


def dem_canh(F):
    """(số cạnh biên, số cạnh dùng > 2 mặt, số mặt trùng tập đỉnh) — cạnh vô hướng."""
    F = np.asarray(F, np.int64)
    if len(F) == 0:
        return 0, 0, 0
    n = int(F.max()) + 1
    E = np.sort(np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]), axis=1)
    _, c = np.unique(E[:, 0] * n + E[:, 1], return_counts=True)
    S = np.sort(F, axis=1); S = S[np.lexsort((S[:, 2], S[:, 1], S[:, 0]))]
    trung = int((S[1:] == S[:-1]).all(axis=1).sum())
    return int((c == 1).sum()), int((c > 2).sum()), trung


def ve_sinh(V, F):
    """Trả (V2, F2, thong_ke) — V2/F2 là CHÍNH V/F nếu không có gì đổi."""
    if not np.isfinite(V).all():
        raise ValueError("toa do khong huu han")
    n = len(V)
    khoa = np.ascontiguousarray(V + V.dtype.type(0))
    khoa = khoa.view(np.dtype((np.void, khoa.dtype.itemsize * khoa.shape[1]))).ravel()
    _, idx_dau, inv = np.unique(khoa, return_index=True, return_inverse=True)
    inv = inv.ravel()
    thu_tu = np.argsort(idx_dau, kind="stable"); hang = np.empty_like(thu_tu); hang[thu_tu] = np.arange(len(thu_tu))
    map1 = hang[inv]
    giu_dinh1 = np.sort(idx_dau)
    so_gop = n - len(giu_dinh1)
    F1 = map1[np.asarray(F, np.int64)]
    lap = (F1[:, 0] == F1[:, 1]) | (F1[:, 1] == F1[:, 2]) | (F1[:, 0] == F1[:, 2])
    F1 = F1[~lap]
    dung = np.zeros(len(giu_dinh1), bool); dung[F1.ravel()] = True
    map2 = np.cumsum(dung) - 1
    F2 = map2[F1]
    giu_dinh = giu_dinh1[dung]
    tk = {"dinh_gop": int(so_gop), "mat_chi_so_lap_bo": int(lap.sum()), "dinh_thua_bo": int((~dung).sum())}
    if so_gop == 0 and not lap.any() and dung.all():
        return V, F, tk
    return V[giu_dinh], F2.astype(F.dtype), tk


def doc_tho(zf):
    return [(i, zf.read(i.filename)) for i in zf.infolist()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vao"); ap.add_argument("ra"); ap.add_argument("--nhat_ky", default=None)
    ap.add_argument("--nen", default="lzma", choices=["lzma", "deflate"])
    a = ap.parse_args()
    vao, ra = os.path.abspath(a.vao), os.path.abspath(a.ra)
    if os.path.lexists(ra):
        sys.exit("ra da co: %s (khong ghi de)" % ra)
    zf = zipfile.ZipFile(vao)
    tho = doc_tho(zf)
    ten_tp = [i.filename for i, _ in tho]
    assert all(t.endswith(".npy") for t in ten_tp), "thanh phan khong phai .npy"
    khoa = [t[:-4] for t in ten_tp]
    bo_phan = sorted(set(k.rsplit("__", 1)[0] for k in khoa if k.endswith("__v") or k.endswith("__f")))
    z = np.load(vao, allow_pickle=False)
    nk = {"vao": vao, "ra": ra, "md5_ma": hashlib.md5(open(os.path.abspath(__file__), "rb").read()).hexdigest(),
          "md5_vao": hashlib.md5(open(vao, "rb").read()).hexdigest(), "nguong_dV_tuong_doi": NGUONG_DV,
          "so_bo_phan": len(bo_phan), "bo_phan_doi": {}, "khong_kin_sau": {}}
    moi = {}
    for b in bo_phan:
        if b + "__v" not in khoa or b + "__f" not in khoa:
            raise SystemExit("bo phan thieu __v hoac __f: %s" % b)
        V = z[b + "__v"]; F = z[b + "__f"]
        V2, F2, tk = ve_sinh(V, F)
        bien0, qua0, trung0 = dem_canh(F)
        if V2 is V and F2 is F:
            if bien0 or qua0:
                nk["khong_kin_sau"][b] = {"canh_bien": bien0, "canh_hon_2_mat": qua0, "mat_trung": trung0, "doi": False}
            continue
        v0, v1 = the_tich(V, F), the_tich(V2, F2)
        dv = abs(v1 - v0) / max(abs(v0), 1e-300) if v0 != 0 else abs(v1 - v0)
        bien1, qua1, trung1 = dem_canh(F2)
        tk.update({"so_dinh": [int(len(V)), int(len(V2))], "so_mat": [int(len(F)), int(len(F2))], "V_truoc": v0, "V_sau": v1,
                   "dV_tuong_doi": dv, "canh_bien": [bien0, bien1], "canh_hon_2_mat": [qua0, qua1], "mat_trung": [trung0, trung1]})
        nk["bo_phan_doi"][b] = tk
        if bien1 or qua1:
            nk["khong_kin_sau"][b] = {"canh_bien": bien1, "canh_hon_2_mat": qua1, "mat_trung": trung1, "doi": True}
        if not dv <= NGUONG_DV:
            print("LOI: %s dV tuong doi %.3e > %g — khong ghi" % (b, dv, NGUONG_DV)); sys.exit(2)
        moi[b + "__v"] = V2; moi[b + "__f"] = F2
    nk["so_bo_phan_doi"] = len(nk["bo_phan_doi"])
    nk["dV_lon_nhat"] = max([x["dV_tuong_doi"] for x in nk["bo_phan_doi"].values()] or [0.0])
    tam = os.path.join(os.path.dirname(ra), ".tam_%d_%s" % (os.getpid(), os.path.basename(ra)))
    if not moi:
        try:
            os.link(vao, tam); nk["cach_ghi"] = "lien_ket_cung (khong doi gi)"
        except OSError:
            shutil.copy2(vao, tam); nk["cach_ghi"] = "chep (khong doi gi)"
    else:
        kieu = zipfile.ZIP_LZMA if a.nen == "lzma" else zipfile.ZIP_DEFLATED
        with zipfile.ZipFile(tam, "w", kieu, allowZip64=True) as o:
            for (i, data), k in zip(tho, khoa):
                if k in moi:
                    buf = io.BytesIO(); np.lib.format.write_array(buf, np.ascontiguousarray(moi[k]), allow_pickle=False)
                    o.writestr(i.filename, buf.getvalue())
                else:
                    o.writestr(i.filename, data)
        nk["cach_ghi"] = "npz moi (%s)" % a.nen
        z2 = zipfile.ZipFile(tam)
        assert [i.filename for i in z2.infolist()] == ten_tp
        for (i, data), k in zip(tho, khoa):
            d2 = z2.read(i.filename)
            if k in moi:
                x = np.load(io.BytesIO(d2), allow_pickle=False)
                assert x.dtype == moi[k].dtype and x.shape == moi[k].shape and x.tobytes() == np.ascontiguousarray(moi[k]).tobytes(), k
            else:
                assert d2 == data, k
    os.replace(tam, ra)
    nk["md5_ra"] = hashlib.md5(open(ra, "rb").read()).hexdigest()
    if a.nhat_ky:
        json.dump(nk, open(a.nhat_ky, "w"), ensure_ascii=False, indent=1)
    print("%s: bo phan doi %d/%d, dV max %.2e, khong kin sau: %s, %s" % (os.path.basename(os.path.dirname(ra)) or ra, nk["so_bo_phan_doi"], len(bo_phan),
          nk["dV_lon_nhat"], {k: (v["canh_bien"], v["canh_hon_2_mat"]) for k, v in nk["khong_kin_sau"].items()}, nk["cach_ghi"]))


if __name__ == "__main__":
    main()
