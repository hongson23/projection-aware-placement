# -*- coding: utf-8 -*-
"""chuan_bi_n14.py (19/11/2026) — N14 (ChinhChuan_DoHoa_Bai1_20261119 §6.3): vật chứa THẬT = màng + nhân tế bào Allen hiPSC đã phân đoạn.
Chạy MỘT lần trên server (CPU, lõi 0–7), TRƯỚC mọi lượt đặt. Ghi bo_dung_chung/n14/:
  chon_te_bao.csv     — bảng đủ điều kiện của cả 200 tế bào (luật chọn ghi trước, mục LUAT_CHON) + thứ tự rút
  te_bao/<CellId>/    — mesh.npz (mang, nhan; µm, tâm = tâm hộp bao màng), khuon_nhan.npz (nhân trong khung đơn vị),
                        sdf.npz (SDF ĐÚNG tới mesh màng + mesh nhân trên lưới h1/2 thẳng hàng với lưới chiếu), meta.json
  dac_ta/dong_vat_n14_<CellId>.json — đặc tả bộ dựng (sinh từ dac_ta/v12/dong_vat_v12.json, các thay đổi ghi trong 'n14_thay_doi')
  chon.json           — 5 tế bào chọn + số lượng bào quan + số liệu dữ liệu cho sẵn (khe nhân–màng ...)
LUAT_CHON (ghi trước, không nhìn kết quả đặt):
  (0) tập gốc: 200 tế bào của cổng 26/10 (mau.json: M0, outlier No, edge_flag 0; 8 tế bào × 25 dòng tế bào, rút phân tầng theo đĩa);
  (1) gate cham_bien = False (màng không chạm biên ảnh cắt);
  (2) vệ sinh mặt nạ: thành phần liên thông lớn nhất của màng ≥ 99 % số voxel màng, của nhân ≥ 99 % số voxel nhân; nhân ⊂ màng ≥ 99,9 % voxel;
  (3) lưới chiếu cùng bước h1 = 23/128 µm như ô động vật D4: số nút Nx·Ny·Nz ≤ 3,0 triệu (bộ nhớ GPU; 128³ = 2,1 triệu);
  (4) rút PHÂN TẦNG THEO CỠ (bản 2, sau rà soát NÊN-1, trước mọi kết quả lượt chính): các tế bào đạt (1)–(3) sắp theo thể tích voxel màng
      (thành phần lớn nhất), chia 5 tầng bằng nhau theo thứ hạng; một bộ sinh numpy default_rng(20261119), mỗi tầng (nhỏ → lớn) một hoán vị
      các CellId của tầng (sắp theo số); lấy tế bào ĐẦU TIÊN của tầng dựng được mesh hợp lệ (kín nước, nhất quán chiều); hỏng thì lấy tế bào kế tiếp.
      (Bản 1 rút đều không phân tầng cho 4/5 tế bào trong 8 % nhỏ nhất — rà soát 19/11 — nên đổi; bản 1 chỉ dùng cho việc thử.)
  KHÔNG loại tế bào có nhân chạm màng: khi cả màng và nhân là dữ liệu cho sẵn, khe nhân–màng (T_nhan_mang = 1 µm, 1.4b) và lớp lưới
  nội chất bao nhân (day_nguoc, 4.4b) là tính chất của DỮ LIỆU, không phải biến quyết định; bộ dựng không đặt hàng nào chỉ gồm hình
  cho sẵn nên chạm màng không tự tạo nút vô nghiệm. Số liệu khe được báo theo từng tế bào (min, % diện tích mặt nhân < 0,5 / 1 / 2,5 µm).
SỐ LƯỢNG (ghi trước): ti_thể = round(300·ρ), lysosome = round(100·ρ) (ô D4: 300 + 100), ρ = V_khả_dụng(thật) / V_khả_dụng(D4);
  V_khả_dụng = thể tích (nút lưới chiếu) cách mặt màng vào trong ≥ T_mang (0,5 µm) và cách mặt nhân ra ngoài > D_er + S_roi (2,3 µm)
  (cùng công thức cho cả hai; bỏ nón khuyết lưới nội chất ở cả hai). D4: màng cầu nhiễu (SDF v12.1 đã lưu), nhân elip trục dv.khuon_nhan_truc
  với s = E[3,9 + 1,1·σ(z)], z ~ N(0,1), đặt ở tâm. ⇒ mật độ bào quan trong vùng khả dụng ≈ D4.
"""
import os, sys, json, csv, hashlib, math, time
import numpy as np, tifffile, trimesh, igl
from scipy import ndimage as ndi
from skimage.measure import marching_cubes

GOC = os.path.expanduser("~/Vision_team_data_189/sonhh20/")
ALLEN = GOC + "b2demo/ch3_tebao/"; R = GOC + "ch3_tebao_dep/"; OUT = R + "bo_dung_chung/n14/"
VX = 0.108333; H1 = 23 / 128; LE = 2.5; NMAX = 3_000_000; HAT = 20261119; N_CHON = 5
SIGMA = (2.0, 1.0, 1.0)
T_MANG, D_ER, S_ROI, T_NHAN_MANG = 0.5, 2.0, 0.3, 1.0
DT_GOC = R + "dac_ta/v12/dong_vat_v12.json"; SDF_D4 = R + "bo_dung_chung/v121/sdf_mang/sdf_75005fbc7d88018f.npz"
SIGN = igl.SignedDistanceType.SIGNED_DISTANCE_TYPE_PSEUDONORMAL


def lon_nhat(b):
    l, n = ndi.label(b)
    if n == 0: return b, 0.0
    c = np.bincount(l.ravel())[1:]; k = int(c.argmax()) + 1
    return l == k, float(c.max() / c.sum())


def nap(m):
    a = tifffile.imread(ALLEN + "crops/" + os.path.basename(m["crop_seg"]))
    return a[:, 0] > 0, a[:, 1] > 0


def mesh_tu_mat_na(b):
    """mặt nạ (z,y,x) → mesh µm (x,y,z): lọc Gauss σ = SIGMA voxel, marching cubes mức 0,5 (đệm 6 voxel ⇒ kín), một khối lớn nhất."""
    f = ndi.gaussian_filter(np.pad(b.astype(np.float32), 6), SIGMA)
    v, fa, _, _ = marching_cubes(f, 0.5, spacing=(VX, VX, VX))
    v = (v - 6 * VX)[:, ::-1].copy(); fa = fa[:, ::-1].copy()
    m = trimesh.Trimesh(v, fa, process=True)
    ps = m.split(only_watertight=False)
    if len(ps) > 1: m = max(ps, key=lambda x: abs(x.volume))
    if m.volume < 0: m.invert()
    return m, len(ps)


def kiem_mesh(m):
    return {"kin_nuoc": bool(m.is_watertight), "nhat_quan_chieu": bool(m.is_winding_consistent), "the_tich_um3": float(m.volume),
            "so_mat": int(len(m.faces)), "euler": int(m.euler_number), "dien_tich_um2": float(m.area)}


def sdf_luoi(m, axes, lat=16):
    """SDF có dấu ĐÚNG tới mesh (igl: khoảng cách điểm–tam giác chính xác, dấu theo giả pháp tuyến) tại mọi nút lưới axes."""
    V = np.asarray(m.vertices, float); F = np.asarray(m.faces, np.int64); n = [len(a) for a in axes]
    out = np.empty(n, np.float32)
    for i in range(0, n[0], lat):
        i1 = min(n[0], i + lat); Q = np.stack(np.meshgrid(axes[0][i:i1], axes[1], axes[2], indexing="ij"), -1).reshape(-1, 3)
        out[i:i1] = igl.signed_distance(Q, V, F, sign_type=SIGN)[0].reshape(i1 - i, n[1], n[2])
    return out


def sdf_elip(P, a):
    a = np.asarray(a, float); a2 = a ** 2; f0 = ((P / a) ** 2).sum(1) - 1.0; ra = f0 > 0
    lo = np.where(ra, 0.0, -a2.min() + 1e-15); hi = np.where(ra, np.maximum(1.0, (np.abs(P) * a).max(1) * 2 + a2.max()), 0.0)
    for _ in range(200):
        mm = 0.5 * (lo + hi); gm = ((a * P / (a2 + mm[:, None])) ** 2).sum(1) - 1.0; lo = np.where(gm > 0, mm, lo); hi = np.where(gm > 0, hi, mm)
    l = 0.5 * (lo + hi); Q = a2 * P / (a2 + l[:, None]); d = np.linalg.norm(P - Q, axis=1)
    return np.where(ra, d, -d)


def kd_D4():
    z = np.load(SDF_D4); S = z["sdf"]; N1 = 128; A = 11.5
    s_nut = S[2:2 * N1 + 2:2, 2:2 * N1 + 2:2, 2:2 * N1 + 2:2].reshape(-1)
    ax = np.linspace(-A + H1 / 2, A - H1 / 2, N1); X = np.stack(np.meshgrid(ax, ax, ax, indexing="ij"), -1).reshape(-1, 3)
    zz = np.linspace(-8, 8, 160001); w = np.exp(-zz ** 2 / 2); w /= w.sum(); s_tb = float(((3.9 + 1.1 / (1 + np.exp(-zz))) ** 3 * w).sum() ** (1 / 3))
    trong = s_nut < -T_MANG; d = np.full(len(X), 1e9); d[trong] = sdf_elip(X[trong], np.array([1.0, 0.943, 0.896]) * s_tb)
    kd = float(((s_nut < -T_MANG) & (d > D_ER + S_ROI)).sum()) * H1 ** 3
    return kd, s_tb


def main():
    os.makedirs(OUT + "te_bao", exist_ok=True); os.makedirs(OUT + "dac_ta", exist_ok=True)
    mau = json.load(open(ALLEN + "mau.json")); gate = {r["CellId"]: r for r in csv.DictReader(open(ALLEN + "ket_qua/te_bao.csv"))}
    bang = []
    for m in sorted(mau, key=lambda x: int(x["CellId"])):
        nuc, mem = nap(m); g = gate[m["CellId"]]
        mem1, fm = lon_nhat(mem); nuc1, fn = lon_nhat(nuc)
        trong = float((nuc & mem).sum() / max(nuc.sum(), 1))
        p = np.argwhere(mem1); ext = (p.max(0) - p.min(0) + 1)[::-1] * VX
        NV = 2 * np.ceil((ext / 2 + LE) / H1).astype(int); nut = int(np.prod(NV))
        r = {"CellId": m["CellId"], "structure_name": m["structure_name"], "PlateId": m["PlateId"], "cham_bien": g["cham_bien"],
             "min_nhan_mang_cong": float(g["min_tat_ca"]), "lon_nhat_mang": round(fm, 5), "lon_nhat_nhan": round(fn, 5), "nhan_trong_mang": round(trong, 6),
             "ext_x": round(ext[0], 2), "ext_y": round(ext[1], 2), "ext_z": round(ext[2], 2), "so_nut": nut, "V_mang_vox_um3": round(float(mem1.sum()) * VX ** 3, 2)}
        r["dat"] = (g["cham_bien"] == "False") and fm >= 0.99 and fn >= 0.99 and trong >= 0.999 and nut <= NMAX
        bang.append(r)
    dat = [r for r in bang if r["dat"]]; dat_sx = sorted(dat, key=lambda r: (r["V_mang_vox_um3"], int(r["CellId"])))
    tang = [dat_sx[i * len(dat_sx) // N_CHON:(i + 1) * len(dat_sx) // N_CHON] for i in range(N_CHON)]; rg = np.random.default_rng(HAT); thu_tu_tang = []
    for i, t in enumerate(tang):
        ids = sorted((r["CellId"] for r in t), key=int); thu_tu_tang.append([ids[j] for j in rg.permutation(len(ids))])
        for r in t: r["tang"] = i
    for r in bang:
        r.setdefault("tang", ""); r["thu_tu_rut"] = thu_tu_tang[r["tang"]].index(r["CellId"]) if r["tang"] != "" else ""
    dat = [r["CellId"] for r in dat]
    with open(OUT + "chon_te_bao.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(bang[0])); w.writeheader(); w.writerows(bang)
    n_cham = sum(r["min_nhan_mang_cong"] == 0 for r in bang)
    print("du dieu kien", len(dat), "/", len(bang), "; nhan cham mang (gate min_tat_ca = 0):", n_cham, "/", len(bang), flush=True)

    KD4, s_tb = kd_D4(); print("V_kha_dung D4", round(KD4, 1), "s_tb", round(s_tb, 4), flush=True)
    DT0 = json.load(open(DT_GOC)); md5_dt = hashlib.md5(open(DT_GOC, "rb").read()).hexdigest()
    chon, hong = [], []
    for cid in [c for tt_ in thu_tu_tang for c in tt_]:
        if len(chon) >= N_CHON: break
        if any(m_["tang"] == [r for r in bang if r["CellId"] == cid][0]["tang"] for m_ in chon): continue
        t0 = time.time(); m = [x for x in mau if x["CellId"] == cid][0]; nuc, mem = nap(m)
        mem, _ = lon_nhat(mem); mem = ndi.binary_fill_holes(mem); nuc, _ = lon_nhat(nuc & mem); nuc = ndi.binary_fill_holes(nuc) & mem
        Mm, km = mesh_tu_mat_na(mem); Mn, kn = mesh_tu_mat_na(nuc)
        km_, kn_ = kiem_mesh(Mm), kiem_mesh(Mn)
        ly = [k for k, v in (("mang_khong_kin", not km_["kin_nuoc"]), ("nhan_khong_kin", not kn_["kin_nuoc"]), ("mang_chieu", not km_["nhat_quan_chieu"]),
                              ("nhan_chieu", not kn_["nhat_quan_chieu"])) if v]
        if ly: hong.append({"CellId": cid, "ly_do": ly}); print("HONG", cid, ly, flush=True); continue
        c0 = (Mm.vertices.max(0) + Mm.vertices.min(0)) / 2; Mm.apply_translation(-c0); Mn.apply_translation(-c0)
        ext = Mm.vertices.max(0) - Mm.vertices.min(0)
        NV = 2 * np.ceil((ext / 2 + LE) / H1).astype(int); AV = NV * H1 / 2
        NF = 2 * NV + 3; LOF = -AV - H1 / 2; HIF = AV + H1 / 2
        axes = [np.linspace(LOF[i], HIF[i], NF[i]) for i in range(3)]
        assert all(abs((axes[i][1] - axes[i][0]) - H1 / 2) < 1e-12 for i in range(3))
        sm = sdf_luoi(Mm, axes); sn = sdf_luoi(Mn, axes)
        doc = [float(np.abs(np.diff(s, axis=a)).max()) / (H1 / 2) for s in (sm, sn) for a in range(3)]
        ax1 = [np.linspace(-AV[i] + H1 / 2, AV[i] - H1 / 2, NV[i]) for i in range(3)]
        assert all(np.allclose(axes[i][2:2 * NV[i] + 2:2], ax1[i], atol=1e-9) for i in range(3))
        sm_n = sm[2:2 * NV[0] + 2:2, 2:2 * NV[1] + 2:2, 2:2 * NV[2] + 2:2]; sn_n = sn[2:2 * NV[0] + 2:2, 2:2 * NV[1] + 2:2, 2:2 * NV[2] + 2:2]
        kd = float(((sm_n < -T_MANG) & (sn_n > D_ER + S_ROI)).sum()) * H1 ** 3; rho = kd / KD4
        n_ti, n_ly = max(1, int(round(300 * rho))), max(1, int(round(100 * rho)))
        io = np.unravel_index(int(np.argmin(sn)), sn.shape); o = np.array([axes[i][io[i]] for i in range(3)]); r_in = float(-sn.min())
        S = float(np.linalg.norm(Mn.vertices - o, axis=1).max()); Vu = (np.asarray(Mn.vertices) - o) / S
        dv = -igl.signed_distance(np.asarray(Mn.vertices, float), np.asarray(Mm.vertices, float), np.asarray(Mm.faces, np.int64), sign_type=SIGN)[0]
        wa = np.zeros(len(Mn.vertices)); np.add.at(wa, Mn.faces.ravel(), np.repeat(Mn.area_faces / 3, 3)); wa /= wa.sum()
        khe = {"min_um": float(dv.min()), "p1_um": float(np.quantile(dv, 0.01)), "trung_vi_um": float(np.median(dv)),
               **{f"dien_tich_duoi_{t}": float(wa[dv < t].sum()) for t in (H1, 0.5, T_NHAN_MANG, D_ER + T_MANG)},
               "nhan_ngoai_mang_dinh": int((dv < 0).sum()), "cham": bool(dv.min() < H1)}
        d = OUT + f"te_bao/{cid}/"; os.makedirs(d, exist_ok=True)
        np.savez_compressed(d + "mesh.npz", mang__v=np.asarray(Mm.vertices, np.float64), mang__f=Mm.faces.astype(np.int32),
                            nhan__v=np.asarray(Mn.vertices, np.float64), nhan__f=Mn.faces.astype(np.int32))
        np.savez_compressed(d + "khuon_nhan.npz", v=Vu.astype(np.float64), f=np.asarray(Mn.faces, np.int64))
        np.savez(d + "sdf.npz", sdf=sm, sdf_nhan=sn, lo=LOF, hi=HIF, NV=NV, AV=AV, h1=H1)
        meta = {"tang": [r for r in bang if r["CellId"] == cid][0]["tang"], "CellId": cid, "structure_name": m["structure_name"], "PlateId": m["PlateId"], "FOVId": m["FOVId"], "crop_seg": m["crop_seg"],
                "nguon": "Allen Cell Explorer hiPSC single-cell image dataset (allencell.org; trích dẫn theo Citation Policy)",
                "vx_um": VX, "tam_hop_bao_goc_um_xyz": c0.tolist(), "ext_um": ext.tolist(), "NV": NV.tolist(), "AV": AV.tolist(), "so_nut": int(np.prod(NV)),
                "mesh_mang": km_, "mesh_nhan": kn_, "so_khoi_mc": [km, kn], "doc_roi_rac_max": doc, "V_kha_dung_um3": kd, "V_kha_dung_D4_um3": KD4,
                "rho": rho, "so_ti_the": n_ti, "so_lysosome": n_ly, "nhan_goc_um": o.tolist(), "nhan_S_um": S, "nhan_r_trong_um": r_in, "khe_nhan_mang": khe,
                "V_te_bao_chat_um3": km_["the_tich_um3"] - kn_["the_tich_um3"], "giay": round(time.time() - t0, 1),
                "md5_sdf": hashlib.md5(open(d + "sdf.npz", "rb").read()).hexdigest()}
        json.dump(meta, open(d + "meta.json", "w"), ensure_ascii=False, indent=1)
        json.dump(dac_ta(DT0, md5_dt, cid, d, meta), open(OUT + f"dac_ta/dong_vat_n14_{cid}.json", "w"), ensure_ascii=False, indent=1)
        chon.append(meta); print("CHON", cid, m["structure_name"], "ext", np.round(ext, 1), "nut", int(np.prod(NV)), "rho", round(rho, 3), "ti/ly", n_ti, n_ly,
                                 "khe min", round(khe["min_um"], 3), "giay", meta["giay"], flush=True)
    json.dump({"luat": __doc__, "V_kha_dung_D4_um3": KD4, "s_tb_D4": s_tb, "so_du_dieu_kien": len(dat), "so_goc": len(bang),
               "nhan_cham_mang_cong_2610": [n_cham, len(bang)], "thu_tu_rut_theo_tang": [t[:6] for t in thu_tu_tang],
               "bien_tang_V_mang_um3": [[t[0]["V_mang_vox_um3"], t[-1]["V_mang_vox_um3"]] for t in tang], "chon": chon, "hong": hong},
              open(OUT + "chon.json", "w"), ensure_ascii=False, indent=1)


def dac_ta(DT0, md5_dt, cid, d, meta):
    """đặc tả N14 từ dong_vat_v12.json: chỉ đổi phần hình học vật chứa + khởi tạo/đặt lại theo hình thật; mọi luật/tham số luật giữ nguyên."""
    D = json.loads(json.dumps(DT0)); TS = D["tham_so_khong_thuoc_luat"]; TH = D["tham_so_so_hoc"]
    kn = f"nhan_that_{cid}"; S = meta["nhan_S_um"]; o = meta["nhan_goc_um"]
    D["phien_ban"] = f"n14v3c-dong_vat-{cid}-20261120"
    D["vo"] = {"ten": "mang", "hinh": "luoi", "tep_sdf": d + "sdf.npz", "md5_sdf": meta["md5_sdf"], "nua_hop_xyz": meta["AV"], "nua_bao_xyz": (np.asarray(meta["ext_um"]) / 2).tolist(), "h": H1,
               "nguon": f"1.2a/1.2b; hình màng THẬT Allen hiPSC CellId {cid}"}
    D["khuon"].append({"ten": kn, "hinh": "luoi_tep", "tep": d + "khuon_nhan.npz", "nguon": f"nhân THẬT Allen CellId {cid} (khung đơn vị: gốc = điểm sâu nhất, cỡ S = {S:.4f} µm)"})
    D["thu_vien"] = {"T": {"co_san": {"nhan": [kn], "golgi": ["golgi"], "lysosome": ["@khuon_lysosome"], "ti_the": ["ong0", "ong1", "ong2", "ong3", "ong4"]}}}
    for b in D["bo_phan"]:
        if b["ten"] == "nhan": b.update(co_dinh=True, co_dinh_tep_sdf=d + "sdf.npz", co_dinh_t=o, co_dinh_r=[1, 0, 0, 0, 1, 0])
    TS["dai_nhan"] = {"gia_tri": [S, S], "nguon": "N14: cỡ nhân thật (cố định)"}
    R_nc = TS["kt_nhancon_trong"]["gia_tri"] * meta["nhan_r_trong_um"] / 0.896
    TS["n14_nhancon_ban_kinh"] = {"gia_tri": R_nc, "nguon": "N14: kt_nhancon_trong × s_tương_đương, s_tương_đương = r_nội_tiếp / 0,896 (elip dv.khuon_nhan_truc có cùng bán kính nội tiếp)"}
    TS["n14_tam_nhan"] = {"gia_tri": o, "nguon": "N14: điểm sâu nhất của nhân thật"}
    TH["quanh_mat_xa"] = {"gia_tri": 15.0, "nguon": "N14: tầm tia tìm khe quanh mặt nhân (µm)"}
    D["quan_he"]["lech_tam"] = []
    D["ngau_nhien"] = {"goc_lam_lai": "nhan"}
    tt_trong = lambda sau: ["vung_mau", ["+", "@T_mang", "$DEL"], "nhan", ["+", "@S_roi", "$DEL"], 1e9]
    tt_tt = lambda s_: ["vung_mau", ["+", "@T_mang", "$DEL", s_], "nhan", ["+", "@S_roi", "$DEL", s_], ["+", "@gan_tt_nhan", s_]]
    a_g = lambda s_: ["*", "@kt_golgi_nua_day", s_]
    tt_g = lambda s_: ["vung_mau", ["+", "@T_mang", "$DEL", a_g(s_)], "nhan", ["+", "@S_roi", "$DEL", a_g(s_)], "@kt_lai_golgi#1"]
    tt_cam = ["vung_mau_cam", ["+", "@T_mang", "$DEL"], "nhan", ["+", "@S_roi", "$DEL"], 1e9, ["+", "@S_roi", "$DEL"]]
    nc = ["vung_mau", 0.0, "nhan", -1e9, ["-", 0.0, "@T_nhancon", "$DEL", "@dai_nhan_con#0"]]
    D["khoi_tao"]["luat"] = [
        {"gan": "cN", "bt": ["vec", "@n14_tam_nhan"]},
        {"dat": "nhan", "t": "$cN", "r": ["vec", "@tay_quay_don_vi"]},
        {"dat": "trung_the", "t": tt_tt("$s0.trung_the"), "r": ["chuan", 6]},
        {"dat": "golgi", "t": tt_g("$s0.golgi"), "r": ["chuan", 6]},
        {"dat": "nhan_con", "t": nc, "r": ["chuan", 6]},
        {"dat": "ti_the", "t": tt_trong("@kt_sau_ti_the"), "r": ["chuan", 6]},
        {"dat": "lysosome", "t": tt_trong("@kt_sau_lyso"), "r": ["chuan", 6]}]
    D["khoi_tao"]["deu"] = [
        {"gan": "cN", "bt": ["vec", "@n14_tam_nhan"]}, {"gan": "tT", "bt": ["trong_deu", "@kt_deu_sau_tt"]}, {"gan": "tG", "bt": ["trong_deu", "@kt_deu_sau_golgi"]},
        {"gan": "rG", "bt": ["chuan", 6]}, {"dat": "nhan", "t": "$cN", "r": ["vec", "@tay_quay_don_vi"]},
        {"dat": "trung_the", "t": "$tT", "r": ["chuan", 6]}, {"dat": "golgi", "t": "$tG", "r": "$rG"},
        {"dat": "nhan_con", "t": nc, "r": ["chuan", 6]},
        {"dat": "ti_the", "t": ["trong_deu", "@kt_sau_ti_the"], "r": ["chuan", 6]}, {"dat": "lysosome", "t": ["trong_deu", "@kt_sau_lyso"], "r": ["chuan", 6]}]
    ncl = nc
    D["dat_lai"]["luat"] = {"nhan": D["dat_lai"]["luat"]["nhan"],
                            "trung_the": {"t": tt_tt("$s.trung_the"), "r": ["chuan", 6]},
                            "golgi": {"t": tt_g("$s.golgi"), "r": ["chuan", 6]},
                            "nhan_con": {"t": ncl, "r": ["chuan", 6]},
                            "ti_the": {"t": tt_cam, "r": ["chuan", 6]}, "lysosome": {"t": tt_cam, "r": ["chuan", 6]}}
    D["dat_lai"]["deu"]["nhan_con"] = {"t": ncl, "r": ["chuan", 6]}
    D.pop("tu_the_co_dinh", None)
    D["n14_thay_doi"] = {"tu": "dac_ta/v12/dong_vat_v12.json", "md5_tu": md5_dt, "te_bao": cid,
                         "doi": ["vo: cau_nhieu → luoi (SDF đúng tới mesh màng thật)", "khuon nhan_that (mesh nhân thật), nhan co_dinh (tư thế + cỡ cố định)",
                                 "lech_tam: bỏ (nhân cố định)", "bản 3 (20/11): cận dưới cộng β; nhân con đề xuất đều trong vùng hợp lệ của nhân thật", "khoi_tao/dat_lai (bản 2, ĐK-1): tâm nhân = điểm sâu nhất nhân thật; mọi bộ phận đề xuất ĐỀU trên vùng hợp lệ của vật chứa thật (vung_mau: độ sâu trong màng + dải khe tới mặt nhân, xem chú thích mã), cùng phép cho mọi nhánh; Golgi xoay ngẫu nhiên",
                                 "goc_lam_lai: nhan như D4 (bản 3b; bộ dựng: gốc cố định ⇒ làm lại cả nhóm)",
                                 "thu_vien chỉ T", "mọi luật + tham số luật (T_*, S_roi, D_er, dải cỡ bào quan, khe) GIỮ NGUYÊN"]}
    return D


def viet_lai_dac_ta():
    DT0 = json.load(open(DT_GOC)); md5_dt = hashlib.md5(open(DT_GOC, "rb").read()).hexdigest()
    for c in json.load(open(OUT + "chon.json"))["chon"]:
        d = OUT + f"te_bao/{c['CellId']}/"; meta = json.load(open(d + "meta.json"))
        assert meta["md5_sdf"] == hashlib.md5(open(d + "sdf.npz", "rb").read()).hexdigest()
        json.dump(dac_ta(DT0, md5_dt, c["CellId"], d, meta), open(OUT + f"dac_ta/dong_vat_n14_{c['CellId']}.json", "w"), ensure_ascii=False, indent=1); print("dac_ta", c["CellId"])


if __name__ == "__main__":
    viet_lai_dac_ta() if sys.argv[1:] == ["dac_ta"] else main()
