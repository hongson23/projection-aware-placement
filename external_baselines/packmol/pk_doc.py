"""pk_doc.py (20/11/2026, bản thử dev) — đọc kết quả Packmol (pk/pk_out.xyz) ⇒ tư thế cứng từng bộ phận (Kabsch, x = t + R·y)
⇒ tu_the_pk.json CÙNG ĐỊNH DẠNG tu_the_cp.json của nhánh cellPACK ({dat: {ten: {t, R}}, seed_cp, md5_dau_vao, giay_dat}),
để cau_noi_*.py nhap (PROJECT + mesh + thong_tin.json) chạy y như nhánh cellPACK. Đồng thời tự đo (chỉ chẩn đoán, KHÔNG phải
thước của bài): chồng cầu–cầu lớn nhất giữa hai bộ phận, và RMSD Kabsch.
Dùng: python pk_doc.py --thu_muc D --giay_pk <giây tường Packmol> [--ra pk]
"""
import argparse, json, os, re
import numpy as np
from scipy.spatial import cKDTree
ap = argparse.ArgumentParser(); ap.add_argument("--thu_muc", required=True); ap.add_argument("--giay_pk", type=float, required=True)
ap.add_argument("--ra", default="pk"); ap.add_argument("--log", default="pk.log")
A = ap.parse_args(); D = A.thu_muc; R = os.path.join(D, A.ra)
meta = json.load(open(os.path.join(R, "pk_meta.json"))); dv = json.load(open(os.path.join(D, "dau_vao_cellpack.json"))); U = meta["U"]
lines = open(os.path.join(R, "pk_out.xyz")).read().split("\n"); n = int(lines[0]); X = np.array([[float(v) for v in l.split()[1:4]] for l in lines[2:2 + n]]) / U
assert n == meta["so_nguyen_tu"], (n, meta["so_nguyen_tu"])
bp = {b["ten"]: b for b in dv["bo_phan"]}; dat = {}; o = 0; rmsd_max = 0.0; C_all = []; R_all = []; I_all = []
for j, e in enumerate(meta["thu_tu"]):
    b = bp[e["ten"]]; Y = np.asarray(b["tam"], float); Xi = X[o:o + len(Y)]; o += len(Y)
    yc, xc = Y.mean(0), Xi.mean(0); H = (Y - yc).T @ (Xi - xc); Uu, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ Uu.T)) or 1.0; Rm = Vt.T @ np.diag([1, 1, d]) @ Uu.T; t = xc - Rm @ yc
    rmsd = float(np.sqrt(((Y @ Rm.T + t - Xi) ** 2).sum(1).mean())); rmsd_max = max(rmsd_max, rmsd)
    assert abs(np.linalg.det(Rm) - 1) < 1e-9
    dat[e["ten"]] = {"t": t.tolist(), "R": Rm.tolist()}
    C_all.append(Xi); R_all.append(np.asarray(b["r"], float)); I_all.append(np.full(len(Y), j))
C = np.vstack(C_all); Rr = np.concatenate(R_all); I = np.concatenate(I_all)
tr = cKDTree(C); pr = tr.query_pairs(2 * Rr.max(), output_type="ndarray"); pr = pr[I[pr[:, 0]] != I[pr[:, 1]]]
cho = Rr[pr[:, 0]] + Rr[pr[:, 1]] - np.linalg.norm(C[pr[:, 0]] - C[pr[:, 1]], axis=1)
log = open(os.path.join(R, A.log)).read() if os.path.exists(os.path.join(R, A.log)) else ""
def tim(p):
    m = re.findall(p, log); return m[-1] if m else None
pk = {"phien_ban": "Packmol.jl 100.0.0 (bản thực thi linux-x86_64, sha256 d6994083…)", "seed_pk": meta["seed_pk"], "U": U, "buoc_diem": meta["buoc_diem"],
      "nloop": meta["nloop"], "so_nguyen_tu": n, "rmsd_kabsch_max_um": rmsd_max,
      "chong_cau_max_um": float(cho.max()) if len(cho) else 0.0, "so_cap_cau_chong_tren_0.01um": int((cho > 0.01).sum()),
      "so_bo_phan_co_cau_chong_tren_0.01um": int(len(np.unique(I[pr[cho > 0.01].ravel()]))) if len(cho) else 0,
      "hoi_tu": ("Success!" in log) or ("Packing converged" in log), "khoang_cach_min_log": tim(r"Minimum distance:\s+([0-9.eE+-]+)"),
      "vi_pham_rang_buoc_max_log": tim(r"Maximum violation of the constraints:\s+([0-9.eE+-]+)"), "giay_pk": A.giay_pk}
json.dump({"dat": dat, "seed_cp": meta["seed_pk"], "md5_dau_vao": meta["md5_dau_vao"], "giay_dat": A.giay_pk, "packmol": pk},
          open(os.path.join(D, "tu_the_pk.json"), "w"), ensure_ascii=False, indent=1)
print("PK_DOC_XONG", json.dumps(pk, ensure_ascii=False), flush=True)
