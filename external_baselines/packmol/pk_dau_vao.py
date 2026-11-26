"""pk_dau_vao.py (20/11/2026, bản 2 sau rà soát: tách vung_cam theo thành phần) — dịch ĐẦU VÀO CHUNG của nhánh cellPACK (dau_vao_cellpack.json + 3 mesh ngăn do
cau_noi_v122.py xuat sinh) sang tệp nhập Packmol (Packmol.jl v100.0.0, bản thực thi, không cần Julia).

Không đổi gì của bài toán đã xuất: cùng mesh ngăn (đã dời ±a theo SDF thật), cùng cây cầu đã nở a.
  * mỗi bộ phận (ti thể, lysosome, nhân con) = một 'structure' number 1; nguyên tử = tâm cầu (khung bộ phận, µm·U),
    bán kính nguyên tử = bán kính cầu đã nở (r·s + a) ⇒ hai cầu khác bộ phận cách nhau ≥ r_i + r_j (khe 2a = S_chung + β như cellPACK).
  * ngăn = mặt cho bằng điểm có pháp tuyến (lấy mẫu đều trên mesh ngăn, bước --buoc_diem µm, pháp tuyến mặt hướng ra ngoài).
    Ràng buộc trên TÂM nguyên tử, lùi đúng bán kính của nó: 'inside surface vo_chua.dat −r', 'outside surface vung_cam_<i>.dat +r' (mỗi thành phần)
    (bộ phận tế bào chất); 'inside surface long_nhan.dat −r' (nhân con). ⇒ cả cầu nằm trong ngăn (cellPACK chỉ kiểm tâm + lưới).
  * đơn vị: 1 đơn vị Packmol = 1/U µm (U = 100 ⇒ 0,01 µm); mọi tham số khác của Packmol để MẶC ĐỊNH của tác giả.
Dùng: python pk_dau_vao.py --thu_muc D --seed_pk S [--U 100] [--buoc_diem 0.05] [--nloop 0 (= mặc định)]
"""
import argparse, hashlib, json, os, sys, time
import numpy as np, trimesh
ap = argparse.ArgumentParser()
ap.add_argument("--thu_muc", required=True); ap.add_argument("--seed_pk", type=int, required=True)
ap.add_argument("--U", type=float, default=100.0); ap.add_argument("--buoc_diem", type=float, default=0.05)
ap.add_argument("--nloop", type=int, default=0); ap.add_argument("--ra", default="pk")
A = ap.parse_args(); t0 = time.time()
D = A.thu_muc; R = os.path.join(D, A.ra); os.makedirs(R, exist_ok=True); U = A.U
dv_path = os.path.join(D, "dau_vao_cellpack.json"); dv = json.load(open(dv_path))
md5_dv = hashlib.md5(open(dv_path, "rb").read()).hexdigest()
mat = {}; TEP = {}
for ten in ("vo_chua", "vung_cam", "long_nhan"):
    m0 = trimesh.load(os.path.join(D, dv["ngan"][ten]["tep"]), process=True, force="mesh")
    assert m0.is_watertight, ten
    cs = m0.split(only_watertight=False); assert all(c.is_watertight for c in cs), ten
    if ten != "vung_cam": assert len(cs) == 1, (ten, len(cs))
    TEP[ten] = []
    for ic, m in enumerate(cs):
        if m.volume < 0: m.invert()
        tep = ten + (f"_{ic}" if ten == "vung_cam" else "") + ".dat"
        n = int(np.ceil(m.area / A.buoc_diem ** 2))
        P, fi = trimesh.sample.sample_surface_even(m, n, seed=12345)
        N = m.face_normals[fi]
        np.savetxt(os.path.join(R, tep), np.hstack([P * U, N]), fmt="%.4f %.4f %.4f %.6f %.6f %.6f",
                   header=f"{ten} thanh phan {ic}: {len(P)} diem, buoc {A.buoc_diem} um, U={U}")
        TEP[ten].append(tep)
        mat[tep] = {"so_diem": int(len(P)), "dien_tich_um2": float(m.area), "the_tich_um3": float(m.volume)}
rmin = min(min(b["r"]) for b in dv["bo_phan"])
L = [f"# sinh boi pk_dau_vao.py tu {dv_path} md5 {md5_dv}", f"tolerance {2 * rmin * U:.4f}", "filetype xyz",
     "output pk_out.xyz", f"seed {A.seed_pk}"]
if A.nloop > 0: L.append(f"nloop {A.nloop}")
thu_tu = []
for b in dv["bo_phan"]:
    tam = np.asarray(b["tam"], float); r = np.asarray(b["r"], float); f = f"bp_{len(thu_tu):04d}.xyz"
    with open(os.path.join(R, f), "w") as fh:
        fh.write(f"{len(r)}\n{b['ten']}\n"); [fh.write(f"C {x * U:.5f} {y * U:.5f} {z * U:.5f}\n") for x, y, z in tam]
    L += [f"structure {f}", "  number 1"]
    for i, ri in enumerate(r):
        L += [f"  atoms {i + 1}", f"    radius {ri * U:.5f}"]
        if b["ngan"] == "te_bao_chat":
            L += [f"    inside surface vo_chua.dat {-ri * U:.5f}"] + [f"    outside surface {t_} {ri * U:.5f}" for t_ in TEP["vung_cam"]]
        elif b["ngan"] == "long_nhan":
            L += [f"    inside surface long_nhan.dat {-ri * U:.5f}"]
        else: sys.exit("ngan la " + b["ngan"])
        L += ["  end atoms"]
    L += ["end structure"]
    thu_tu.append({"ten": b["ten"], "tep": f, "so_nguyen_tu": len(r)})
open(os.path.join(R, "pk.inp"), "w").write("\n".join(L) + "\n")
json.dump({"md5_dau_vao": md5_dv, "U": U, "buoc_diem": A.buoc_diem, "seed_pk": A.seed_pk, "nloop": A.nloop or "mac_dinh",
           "mat": mat, "thu_tu": thu_tu, "so_nguyen_tu": int(sum(t["so_nguyen_tu"] for t in thu_tu)),
           "giay_dich": round(time.time() - t0, 2)}, open(os.path.join(R, "pk_meta.json"), "w"), ensure_ascii=False, indent=1)
print("PK_DAU_VAO_XONG", len(thu_tu), "bo_phan", sum(t["so_nguyen_tu"] for t in thu_tu), "nguyen_tu", mat, round(time.time() - t0, 1), "s", flush=True)
