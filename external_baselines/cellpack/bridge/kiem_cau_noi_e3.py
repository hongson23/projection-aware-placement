"""kiem_cau_noi_v122.py — kiểm MỘT lần xuất của cau_noi_v122.py (18/11/2026) có khớp bộ dựng v12.2 của lượt chính không.

Chạy SAU `cau_noi_v122.py xuat` (cwd = ch3_tebao_dep, cùng đối số bộ dựng sau '--'):
  python kiem_cau_noi_v122.py --thu_muc D --tt_nn1x <thong_tin.json nhánh nn1x cùng ô/seed> [--tt_td <thong_tin.json nhánh td>] -- <đối số bộ dựng>
Chạy lại phần A + khoi_tao(k) + khối L của botri_chung_v122.py trên CPU (như cầu nối), rồi so:
  (1) dấu: md5 bộ dựng, SDF màng (tệp, md5), L_M = 1, β, β_M = số của lượt chính (thong_tin.json ban_v12);
  (2) mức luật: muc_vo_luat = −(T_nho + β_M), muc_cam_luat = S_CN + β, muc_long_nhan_luat = −(T_nhancon + β), a = (S_nho_chung + β)/2;
  (3) TRƯỜNG: trên đỉnh mesh ngăn đã xuất, giá trị trường của bộ dựng v12.2 so với mức (mức luật ± a):
        vo_chua : φ_vo(v) − (muc_vo_luat + a_nho)      (φ_vo là SDF đúng ⇒ phải ≈ 0)
        long_nhan: ψ_nhan(v) − (muc_long_nhan + a_nc)   (ψ_nhan = s·SDF khuôn ⇒ ≈ 0)
        vung_cam : f_cam(v) (hợp các vùng cấm ở mức luật; không phải SDF đúng ⇒ chỉ đòi f_cam(v) < 0 và báo phân vị)
  (4) nguyên liệu: tên, loại, khuôn, cỡ s của mọi bộ phận nhỏ + nhân con = bộ dựng (cấu hình k) = thong_tin.json của nhánh nn1x/td
      (khuôn, s chỉ so được khi nhánh đó chọn đúng cấu hình k); r cây cầu = r_khuôn·s + a;
  (5) bộ phận lớn cố định (nhân, trung thể, Golgi): tư thế trong tệp xuất so với tu_the của nn1x khi nn1x chọn đúng k (khối L CPU so với GPU).
Ngưỡng (ghi trước khi chạy): (1)(2)(4) bằng nhau tuyệt đối (số thực: |Δ| ≤ 1e-12); (3) max |sai| ≤ 0,05 µm ở vo_chua và long_nhan,
f_cam(v) < 0 ở mọi đỉnh vung_cam; (5) chỉ báo (không chặn): lech_khoi_L_vs_nn1x (tu_the của thong_tin.json làm tròn 7–9 chữ số ⇒ lệch ~1e-7 là khớp). Kết quả: <D>/kiem_cau_noi.json, dòng KIEM_DAT / KIEM_TRUOT.
"""
import argparse, hashlib, json, os, sys
ap = argparse.ArgumentParser()
ap.add_argument("--thu_muc", required=True); ap.add_argument("--tt_nn1x", default=None); ap.add_argument("--tt_td", default=None)
ap.add_argument("--bo_dung", default="bo_dung_chung/ma3/botri_chung_v122e3.py"); ap.add_argument("--nguong_truong", type=float, default=0.05)
if "--" not in sys.argv: sys.exit("cần '--' trước đối số bộ dựng")
_i = sys.argv.index("--"); A_ = ap.parse_args(sys.argv[1:_i]); ARGV_BD = sys.argv[_i + 1:]
dv = json.load(open(os.path.join(A_.thu_muc, "dau_vao_cellpack.json")))
MD5_BD = "6e3beba8706763433586bfbda1a9151f"; SDF_TEP = "bo_dung_chung/v121/sdf_mang/sdf_75005fbc7d88018f.npz"; SDF_MD5 = "ee14b688e6225311735c3ebfbeecc902"
src = open(A_.bo_dung, encoding="utf-8").read(); m5 = hashlib.md5(src.encode("utf-8")).hexdigest()
if m5 != MD5_BD: sys.exit(f"md5 bộ dựng {m5} ≠ {MD5_BD}")
if not os.path.isfile(SDF_TEP) or hashlib.md5(open(SDF_TEP, "rb").read()).hexdigest() != SDF_MD5: sys.exit("SDF màng thiếu/lệch: dừng (không để bộ dựng dựng lại)")
assert dv["argv_bo_dung"][:len(ARGV_BD)] == ARGV_BD or dv["argv_bo_dung"] == ARGV_BD + ["--cach", "ngau_nhien"], "đối số bộ dựng khác lúc xuất"
MOC_A = 't_dat0 = time.time(); t_day_du = 0.0\n'
MOC_L0 = '        if ARG.cach == "ngau_nhien":   # bộ phận lớn thử–loại tuần tự; rồi kiểm chung cả nhóm, hỏng thì làm lại cả nhóm\n'
MOC_L1 = '        thu_loai_nho(nh, lon, g, bao)\n'
pA = src[:src.index(MOC_A)]; pL = src[src.index(MOC_L0) + len(MOC_L0):src.index(MOC_L1)]
pL = "\n".join(l[12:] if l.startswith(" " * 12) else l for l in pL.splitlines()) + "\n"
import torch, numpy as np, trimesh
pA = pA.replace('DEV = "cuda"', 'DEV = "cpu"', 1)
torch.cuda.synchronize = lambda *a, **k: None; torch.cuda.max_memory_allocated = lambda *a, **k: 0
torch.cuda.get_rng_state = lambda *a, **k: torch.get_rng_state(); torch.cuda.set_rng_state = lambda *a, **k: None
sys.argv = [A_.bo_dung] + dv["argv_bo_dung"]
NS = {"__name__": "__bo_dung__", "__file__": os.path.abspath(A_.bo_dung)}
if os.environ.get("SDF_MANG_THU_MUC") not in (None, "bo_dung_chung/v121/sdf_mang"): sys.exit("SDF_MANG_THU_MUC lệch lượt chính: dừng")
exec(compile(pA, A_.bo_dung + "#A", "exec"), NS)
assert NS["V121"].get("sdf_tep") == SDF_TEP and NS["V121"].get("sdf_md5") == SDF_MD5, "SDF màng của phần A (kiểm) lệch lượt chính"
K = int(dv["k"])
exec(compile(f"lon, nh, g = khoi_tao({K}); HT = {{b.loai: b for b in lon if b.loai in DON}}\n", "kiem#k", "exec"), NS)
exec(compile(pL, A_.bo_dung + "#L", "exec"), NS)
T = NS["T"]; lon, nh = NS["lon"], NS["nh"]; TH = dv["nguong"]; kq = {"thu_muc": A_.thu_muc, "k": K, "loi": [], "bao": {}}
def bang(ten, a, b, tol=1e-12):
    ok = (a == b) if not isinstance(a, float) else abs(a - b) <= tol
    kq["bao"][ten] = {"xuat": a, "bo_dung": b, "khop": bool(ok)}
    if not ok: kq["loi"].append(ten)
bang("md5_bo_dung", dv["md5_bo_dung"], MD5_BD); bang("ban", dv["ban"], "cau_noi_cp_e3")
bang("sdf_tep", dv["ban_v121"]["sdf_tep"], SDF_TEP); bang("sdf_md5", dv["ban_v121"]["sdf_md5"], SDF_MD5)
bang("L_M", float(TH["L_M"]), 1.0); bang("beta", float(TH["beta"]), float(NS["DEL"])); bang("beta_M", float(TH["beta_M"]), float(NS["DEL_M"]))
for nhan_, tep_ in (("nn1x", A_.tt_nn1x), ("td", A_.tt_td)):
    if tep_:
        j = json.load(open(tep_)); b12 = j["ban_v12"]
        bang(f"beta_vs_{nhan_}", float(TH["beta"]), float(b12["beta"])); bang(f"beta_M_vs_{nhan_}", float(TH["beta_M"]), float(b12["beta_M"]))
        bang(f"sdf_md5_vs_{nhan_}", dv["ban_v121"]["sdf_md5"], j["ban_v121"]["sdf_md5"]); bang(f"ban_bo_dung_{nhan_}", j["ban"], "botri_chung_v122e3")
CHUA = NS["CHUA"]; DEL = float(NS["DEL"]); DEL_M = float(NS["DEL_M"])
S_NHO = max([NS["S_NN"]] + list(NS["KCL_NN"].values())); a_nho = (S_NHO + DEL) / 2; a_nc = (NS["KCL_LON"].get("nhan_con", 0.0) + DEL) / 2
bang("muc_vo_luat", float(TH["muc_vo_luat"]), -(CHUA["ti_the"]["t"] + DEL_M)); bang("muc_cam_luat", float(TH["muc_cam_luat"]), NS["S_CN"] + DEL)
bang("muc_long_nhan_luat", float(TH["muc_long_nhan_luat"]), -(CHUA["nhan_con"]["t"] + DEL)); bang("a_nho", float(TH["a_nho"]), a_nho); bang("a_nhan_con", float(TH["a_nhan_con"]), a_nc)
def dinh(ten): return np.asarray(trimesh.load(os.path.join(A_.thu_muc, dv["ngan"][ten]["tep"]), process=False).vertices, float)
def thong_ke(x): x = np.asarray(x, float); return {"n": int(len(x)), "max_abs": float(np.abs(x).max()), "p99_abs": float(np.quantile(np.abs(x), 0.99)), "tb": float(x.mean()), "min": float(x.min()), "max": float(x.max())}
with torch.no_grad():
    v = dinh("vo_chua"); s_vo = NS["phi_vo"](T(v)).cpu().numpy() - (TH["muc_vo_luat"] + TH["doi_vo_um"]); kq["bao"]["truong_vo_chua"] = thong_ke(s_vo)
    q = NS["HT"]["nhan"]; v = dinh("long_nhan"); s_ln = q.phi(T(v)).cpu().numpy() - (TH["muc_long_nhan_luat"] + TH["doi_long_nhan_um"]); kq["bao"]["truong_long_nhan"] = thong_ke(s_ln)
    dt = [b for b in lon if b.loai in NS["RNL"]["doi_tac"]]; v = T(dinh("vung_cam"))
    fc = NS["phi_cam"](v, NS["VANH"][NS["CAM_NHO"]]) - TH["muc_cam_luat"]
    for b in dt: fc = torch.minimum(fc, b.phi(v) - TH["muc_doi_tac_lon_luat"])
    fc = fc.cpu().numpy(); kq["bao"]["truong_vung_cam"] = thong_ke(fc); kq["bao"]["truong_vung_cam"]["so_dinh_f_cam_khong_am"] = int((fc >= 0).sum())
for ten_, x_ in (("truong_vo_chua", s_vo), ("truong_long_nhan", s_ln)):
    if float(np.abs(x_).max()) > A_.nguong_truong: kq["loi"].append(ten_)
if (fc >= 0).any(): kq["loi"].append("truong_vung_cam")
TEN_K = NS["TEN_K"]; S_ = nh.sc(); bd = {}
for i in range(nh.n): bd[nh.ten[i]] = (nh.loai[i], TEN_K[int(nh.tid[i])], float(S_[i]))
for b in lon:
    if b.loai == "nhan_con": bd[b.ten] = (b.loai, b.khuon, float(b.sc()))
xu = {b["ten"]: (b["loai"], b["khuon"], float(b["s"])) for b in dv["bo_phan"]}
bang("tap_nguyen_lieu", sorted(xu), sorted(bd)); kq["bao"]["so_nguyen_lieu"] = len(xu)
lech_nl = [t for t in xu if t in bd and (xu[t][:2] != bd[t][:2] or abs(xu[t][2] - bd[t][2]) > 1e-12)]
kq["bao"]["nguyen_lieu_lech_bo_dung"] = lech_nl[:20]; kq["bao"]["so_nguyen_lieu_lech_bo_dung"] = len(lech_nl)
if lech_nl: kq["loi"].append("nguyen_lieu_vs_bo_dung")
ct = {t["ten"]: t for t in dv["bo_phan"]}
sai_r = 0
for t, (loai_, kh_, s_) in xu.items():
    a_ = a_nc if loai_ == "nhan_con" else a_nho; r = np.asarray(ct[t]["r"]); sai_r += int(not (r.min() > a_ - 1e-12))
bang("r_cay_cau_lon_hon_a", sai_r, 0)
for nhan_, tep_ in (("nn1x", A_.tt_nn1x), ("td", A_.tt_td)):
    if not tep_: continue
    j = json.load(open(tep_)); kc = j.get("bao_cao_bo_tri", {}).get("cau_hinh_chon"); tt = j["tu_the"]; kq["bao"][f"k_chon_{nhan_}"] = kc
    if kc != K: kq["bao"][f"so_voi_{nhan_}"] = f"không so được: {nhan_} chọn k = {kc} ≠ {K}"; continue
    lech = [t for t, (l_, kh_, s_) in xu.items() if t not in tt or tt[t].get("khuon") != kh_ or abs(float(tt[t].get("s", -1)) - s_) > 1e-6]
    kq["bao"][f"nguyen_lieu_lech_{nhan_}"] = lech[:20]; kq["bao"][f"so_nguyen_lieu_lech_{nhan_}"] = len(lech)
    if lech: kq["loi"].append(f"nguyen_lieu_vs_{nhan_}")
    if nhan_ == "nn1x":
        lm = 0.0
        for ten_l, e in dv["lon_co_dinh"].items():
            if ten_l in tt and tt[ten_l].get("t") is not None:
                lm = max(lm, float(np.abs(np.asarray(e["t"]) - np.asarray(tt[ten_l]["t"])).max()))
                if tt[ten_l].get("R") is not None: lm = max(lm, float(np.abs(np.asarray(e["R"]) - np.asarray(tt[ten_l]["R"])).max()))
        kq["bao"]["lech_khoi_L_vs_nn1x"] = lm
kq["ket_luan"] = "KIEM_DAT" if not kq["loi"] else "KIEM_TRUOT"
json.dump(kq, open(os.path.join(A_.thu_muc, "kiem_cau_noi.json"), "w"), ensure_ascii=False, indent=1, default=str)
print(kq["ket_luan"], kq["loi"], json.dumps({k: kq["bao"][k] for k in kq["bao"] if k.startswith("truong") or k.startswith("lech") or k.startswith("k_chon")}, ensure_ascii=False, default=str), flush=True)
sys.exit(0 if not kq["loi"] else 3)
