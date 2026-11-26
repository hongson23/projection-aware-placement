# -*- coding: utf-8 -*-
"""kiem_kin.py (24/11/2026, R4 hình bài CAD) — đếm kín/biên/không đa tạp theo bộ phận trên mesh CUỐI.
Mỗi việc <d> (thư mục kq/ có mesh_vo_bao.npz), làm y như đường chấm cham_lo_vsb.sh (bước 1–2), trên BẢN SAO ở /tmp/hinh_cad_$USER:
  (1) ve_sinh_luoi.py <d>/mesh_vo_bao.npz vs.npz                       (md5 kiểm theo bo_kiem_v3/VE_SINH_LUOI_MD5.txt)
  (2) sua_mesh_dn.py vs.npz sb.npz --vat <dong_vat|thuc_vat> --nhat_ky  (md5 kiểm theo bo_kiem_v3/SUA_MESH_MD5.txt)
      (nén --nen deflate thay lzma mặc định: chỉ đổi nén zip, không đổi mảng)
Rồi với MỖI bộ phận (cặp <ten>__v/__f) ở BA mức: sau chiếu (gốc), sau vệ sinh, sau sửa (CUỐI):
  canh_bien = số cạnh không hướng thuộc đúng 1 mặt; canh_kdt = số cạnh thuộc > 2 mặt;
  wt = trimesh.Trimesh(V, F, process=False).is_watertight; wc = is_winding_consistent (mức cuối).
Chép npz nguồn sang /tmp trước (nguồn có thể bị bộ điều phối xoá giữa chừng; nguồn mất ⇒ MAT_NGUON, bỏ qua). Không sửa tệp vào; xoá bản sao sau khi đếm. Dùng: python kiem_kin.py <ra_thu_muc> <vat> <d>...  ⇒ <ra>/<viec>.json"""
import os, sys, json, subprocess, shutil, hashlib, time
import numpy as np, trimesh
GOC = os.path.expanduser("~/Vision_team_data_189/sonhh20/ch3_tebao_dep"); K = GOC + "/bo_kiem_v3"
P = os.path.expanduser("~/miniforge3/envs/b2demo/bin/python")
TMP = "/tmp/hinh_cad_" + os.environ.get("USER", "u") + "/kin"
def md5(f): return hashlib.md5(open(f, "rb").read()).hexdigest()
for lst in ("VE_SINH_LUOI_MD5.txt", "SUA_MESH_MD5.txt"):
    for l in open(f"{K}/{lst}"):
        h, f = l.split()
        if f.endswith(".py"): assert md5(f"{K}/{f}") == h, (f, "md5 khác")
MD5 = {"ve_sinh_luoi.py": md5(f"{K}/ve_sinh_luoi.py"), "sua_mesh_dn.py": md5(f"{K}/sua_mesh_dn.py")}

def dem(npz, wc=False):
    z = np.load(npz); r = {}
    for k in sorted(f[:-3] for f in z.files if f.endswith("__v")):
        V = np.asarray(z[k + "__v"], np.float64); F = np.asarray(z[k + "__f"], np.int64)
        if len(F) == 0: r[k] = {"nF": 0, "canh_bien": 0, "canh_kdt": 0, "wt": False}; continue
        E = np.sort(np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]), axis=1)
        _, c = np.unique(E[:, 0] * (int(F.max()) + 1) + E[:, 1], return_counts=True)
        m = trimesh.Trimesh(V, F, process=False)
        e = {"nF": int(len(F)), "canh_bien": int((c == 1).sum()), "canh_kdt": int((c > 2).sum()), "wt": bool(m.is_watertight)}
        if wc: e["wc"] = bool(m.is_winding_consistent)
        r[k] = e
    return r

ra, vat = sys.argv[1], sys.argv[2]; os.makedirs(ra, exist_ok=True)
for d in sys.argv[3:]:
    j = os.path.basename(os.path.normpath(d)); fo = f"{ra}/{j}.json"
    if os.path.exists(fo): continue
    t0 = time.time(); w = f"{TMP}/{j}"; shutil.rmtree(w, ignore_errors=True); os.makedirs(w)
    env = dict(os.environ, OMP_NUM_THREADS="1", PYTHONNOUSERSITE="1")
    try: shutil.copy2(f"{d}/mesh_vo_bao.npz", f"{w}/goc.npz")
    except FileNotFoundError: print("MAT_NGUON", j, flush=True); shutil.rmtree(w, ignore_errors=True); continue
    out_md5_goc = md5(f"{w}/goc.npz")
    a = subprocess.run([P, f"{K}/ve_sinh_luoi.py", f"{w}/goc.npz", f"{w}/vs.npz", "--nhat_ky", f"{w}/vs.json", "--nen", "deflate"],
                       env=env, capture_output=True, text=True)
    out = {"viec": j, "vat": vat, "md5": MD5, "md5_goc": out_md5_goc, "ve_sinh_rc": a.returncode}
    if a.returncode != 0: out["loi"] = a.stderr[-2000:]; json.dump(out, open(fo, "w")); print("LOI_VS", j, flush=True); continue
    b = subprocess.run([P, f"{K}/sua_mesh_dn.py", f"{w}/vs.npz", f"{w}/sb.npz", "--vat", vat, "--nhat_ky", f"{w}/sb.json", "--nen", "deflate"],
                       env=env, capture_output=True, text=True)
    out["sua_rc"] = b.returncode
    if b.returncode != 0: out["loi"] = b.stderr[-2000:]; json.dump(out, open(fo, "w")); print("LOI_SB", j, flush=True); continue
    try: out["sua_mesh"] = json.load(open(f"{w}/sb.json"))
    except Exception: out["sua_mesh"] = None
    out["goc"] = dem(f"{w}/goc.npz"); out["ve_sinh"] = dem(f"{w}/vs.npz"); out["cuoi"] = dem(f"{w}/sb.npz", wc=True)
    out["giay"] = round(time.time() - t0, 1); shutil.rmtree(w, ignore_errors=True)
    json.dump(out, open(fo + ".tam", "w")); os.replace(fo + ".tam", fo)
    c = out["cuoi"]; print("XONG", j, "kin", sum(e["wt"] for e in c.values()), "/", len(c), out["giay"], "s", flush=True)
