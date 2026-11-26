# -*- coding: utf-8 -*-
"""cham_mot.py — chấm MỘT tế bào từ thư mục kết quả bộ dựng (mesh_vo_bao.npz + thong_tin.json) hoặc thư mục định dạng
bộ kiểm (mesh.npz + loai.json):
  1) (nếu là thư mục bộ dựng) dựng thư mục bộ kiểm bo_kiem_v3/cham_ra/<ten>/: mesh.npz = liên kết tới mesh_vo_bao.npz,
     loai.json suy từ TÊN bộ phận (bỏ hậu tố số, như tv_conv.py) — không đọc trường 'loai'/'cha' của bộ dựng;
  2) kiem3 (luật sự thật v1.2, hai phán quyết) → kiem3.json / kiem3.csv;
  3) (nếu có thong_tin.json["tu_the"]) do_khuon → D_i từng bộ phận, G_c; thêm dòng 'G_c' và 'D_i' vào kiem3.csv.
  4) (NEN=1, mặc định, 15/11) xoá kiem2_dang_chay.json, kiem2_bo_dem.pkl; nén gzip kiem2.json, kiem3.json, do_khuon.json (giữ kiem3.csv).
Dùng: python cham_mot.py <thu_muc> [--vat dong_vat|thuc_vat] [--lop Sinh10]"""
import os, sys, json, csv, argparse
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)


def chuan_bi(d, vat, lop):
    d = os.path.abspath(d)
    if os.path.exists(os.path.join(d, "mesh.npz")) and os.path.exists(os.path.join(d, "loai.json")):
        return d, None
    if not os.path.exists(os.path.join(d, "mesh_vo_bao.npz")):
        raise FileNotFoundError("khong co mesh.npz/loai.json hay mesh_vo_bao.npz trong %s" % d)
    ten = os.path.basename(os.path.dirname(d)) + "__" + os.path.basename(d)
    o = os.path.join(HERE, "cham_ra", ten); os.makedirs(o, exist_ok=True)
    if not os.path.lexists(os.path.join(o, "mesh.npz")):
        os.symlink(os.path.join(d, "mesh_vo_bao.npz"), os.path.join(o, "mesh.npz"))
    z = np.load(os.path.join(d, "mesh_vo_bao.npz"))
    ns = sorted(set(k.rsplit("__", 1)[0] for k in z.files))
    loai = {n: n.rstrip("0123456789").rstrip("_") for n in ns}
    json.dump({"lop": lop or "Sinh10", "loai": loai, "bo_phan_du_kien": loai, "nguon": d}, open(os.path.join(o, "loai.json"), "w"), ensure_ascii=False, indent=1)
    return o, d


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("thu_muc"); ap.add_argument("--vat", default="dong_vat"); ap.add_argument("--lop", default=None)
    a = ap.parse_args()
    o, goc = chuan_bi(a.thu_muc, a.vat, a.lop)
    import kiem3
    out = kiem3.kiem3(o, a.vat, a.lop)
    if goc and os.path.exists(os.path.join(goc, "thong_tin.json")) and "tu_the" in json.load(open(os.path.join(goc, "thong_tin.json"))):
        import do_khuon
        g = do_khuon.do_te_bao(goc, a.vat)
        json.dump(g, open(os.path.join(o, "do_khuon.json"), "w"), ensure_ascii=False, indent=1)
        with open(os.path.join(o, "kiem3.csv"), "a", newline="") as f:
            w = csv.writer(f)
            w.writerow([out["te_bao"], a.vat, out["lop"], "G_c", "khuon", "gia_tri_um", g["G_c_um"], g["so_bo_phan_do"], sum(1 for v in g["bo_phan"].values() if (v.get("D_um") or 0) > g["h_um"]), "h=%g" % g["h_um"], ""])
            for n, v in g["bo_phan"].items():
                w.writerow([out["te_bao"], a.vat, out["lop"], "D_i", "khuon", n, v.get("D_um", v.get("trang_thai")), v.get("so_diem", 0), int((v.get("D_um") or 0) > g["h_um"]), v.get("khuon"), ""])
        print("   G_c = %.4f µm (h = %g), D_max = %s" % (g["G_c_um"], g["h_um"], g["D_max_um"]))
    if os.environ.get("NEN", "1") == "1":
        import gzip, shutil
        for f in ("kiem2_dang_chay.json", "kiem2_bo_dem.pkl"):
            try: os.remove(os.path.join(o, f))
            except FileNotFoundError: pass
        for f in ("kiem2.json", "kiem3.json", "do_khuon.json"):
            p_ = os.path.join(o, f)
            if os.path.exists(p_):
                with open(p_, "rb") as x, gzip.open(p_ + ".gz", "wb") as y: shutil.copyfileobj(x, y)
                os.remove(p_)
