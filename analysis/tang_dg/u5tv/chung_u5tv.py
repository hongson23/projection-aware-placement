# -*- coding: utf-8 -*-
"""chung_u5tv.py (22/11/2026, W4 thực vật U5) — định nghĩa chung: mức, nhánh, seed, dựng lệnh, chèn/bỏ việc trong hàng tdg.
Bộ dựng: ma3/botri_chung_v123u5b.py (sinh bởi u5tv/tao_ma_tv_u5b.py). Seed: chính 7311–7320; dev 7321–7340 (M2 của rà soát trước chạy).
TV-real: (lục lạp, ti thể) rút bằng np.random.default_rng(20261120), integers(80,121), integers(200,451) lần lượt theo seed 7311, 7312, …
(7311–7320 trùng bảng thiết kế ThietKe_ChonBuocDat §4; 7321–7340 là phần rút tiếp theo CÙNG luật)."""
import os, time, hashlib, shutil, json
import numpy as np
R = os.path.expanduser("~/Vision_team_data_189/sonhh20/ch3_tebao_dep")
O = "bo_dung_chung/tang_dg"; U = f"{O}/u5tv"; KQ = f"{O}/kq"
PYB = "/home/thainq97/miniforge3/envs/b2demo/bin/python"
ENV = "env PYTHONNOUSERSITE=1 PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True"
ENV_IGL = ENV + " IGL_PHU=@bo_dung_chung/tang_dg/u5tv/IGL_PHU.txt IGL_PHU_LUONG=4"
NEW = "bo_dung_chung/ma3/botri_chung_v123u5c.py"
TV = "--dac_ta dac_ta/v12/thuc_vat_v12.json"
_rng = np.random.default_rng(20261120)
REAL = {s: (int(_rng.integers(80, 121)), int(_rng.integers(200, 451))) for s in range(7311, 7341)}
MUC = ("D1", "D5", "real")
def so_muc(muc, s):
    if muc == "D1": return ""
    nl, nm = {"D4": (71, 15), "D5": (113, 15)}.get(muc) or REAL[s]
    return f"--so luc_lap={nl} --so ti_the={nm} --so golgi=8 --so peroxisome=30"
SA_C00 = "--un_T1 1e-4 --un_p0 0.5 --un_buoc_t 0.15 0.01 --un_buoc_r 0.25 0.015 --un_ti_le_lon 0.4 --un_p_doi 0.5 --un_lam_moi 50"
U5KTT = "--cach tu_dong --u5 --u5_lambda 1.0 --u5_kappa 0 --u5_chon gc --u5_tam --ktt --s_co_hat --khoa_kt"
def nhanh(nh, ns=None):
    """ns: ngân sách {'K_nn', 'K_nnm', 'un': [a, b], 'sa': chuỗi, 'saL': bool} (None ⇒ ngân sách v12.3 TV-D1 cũ, chỉ cho M1)."""
    ns = ns or {"K_nn": 109, "K_nnm": 109, "un": [12777, 1866], "sa": SA_C00, "saL": False}
    if nh == "u5ktt": return U5KTT
    if nh == "td": return "--cach tu_dong"
    if nh == "nn1x": return f"--cach ngau_nhien --so_lo 3 --K {ns['K_nn']}"
    if nh == "nn1xm": return f"--cach ngau_nhien --so_lo 3 --K {ns['K_nnm']} --chon_mem --s_co_hat --khoa_kt"
    if nh == "un1x": return f"--cach u_nhiet --un_vong_lon {ns['un'][0]} --un_vong_nho {ns['un'][1]} {ns['sa']}" + (" --sa_lich_p1" if ns.get("saL") else "")
    if nh == "drej": return "--cach ngau_nhien --so_lo 3 --K 4"
    if nh == "dun": return f"--cach u_nhiet --un_vong_lon 100 --un_vong_nho 300 {SA_C00}"
    raise KeyError(nh)
def lenh(muc, s, nh_chuoi, id_, day_du=True, vung=True, env=None):
    base = f"{TV} --seed {s} --N1 128 --K 4 --thu_vien T {so_muc(muc, s)}".replace("  ", " ")
    duoi = "--eval_dung_tu_dau --tang2 --khong_tu_do"
    if day_du: duoi += " --hd_truoc_sau --giu_truoc" + (" --hd_vung_truoc" if vung else "")
    else: duoi += " --chi_dat"
    return f"{env or (ENV_IGL if day_du else ENV)} {PYB} {NEW} {base} {nh_chuoi} {duoi} --out {KQ}/{id_}"
def dong(id_, gd, gb, lo, gi, cmd): return f"{id_} {gd} {gb} {lo} {gi} {cmd}"

V = f"{O}/viec.txt"; M = f"{O}/MD5_TDG.txt"
MA_MD5 = ("bo_dung_chung/ma3/botri_chung_v123u5c.py", f"{U}/tao_ma_tv_u5c.py", "bo_dung_chung/ma3/botri_chung_v123u5b.py", "bo_dung_chung/ma3/botri_cucbo_v123u5b.py", f"{U}/igl_tho.py", f"{U}/tao_ma_tv_u5b.py",
          "bo_dung_chung/ma3/botri_chung_v123u5.py", "bo_dung_chung/ma3/botri_cucbo_v123u5.py", f"{U}/tao_ma_tv_u5.py")
def md5(f): return hashlib.md5(open(f, "rb").read()).hexdigest()
def _chua_khoi(i): return not (os.path.exists(f"{O}/xong/{i}") or os.path.exists(f"{O}/loi/{i}") or os.path.isdir(f"{O}/khoa/{i}"))
def sua_hang(them=(), bo=(), ghi=""):
    """Chèn các dòng 'them' NGAY TRƯỚC việc u5_* chưa khởi đầu tiên; bỏ các id 'bo' (chỉ khi chưa khởi). Sao lưu, ghi .moi rồi os.replace,
    cập nhật md5 viec.txt + mã u5tv trong MD5_TDG.txt. Trả (số chèn, số bỏ)."""
    os.chdir(R)
    them = [l.rstrip("\n") for l in them if l.strip()]
    ids = [l.split()[0] for l in them]
    assert all(i.startswith("u5_tv_") for i in ids) and len(set(ids)) == len(ids), ids
    for l in them:
        p = l.split(None, 5); assert len(p) == 6 and p[1] in ("C", "T1") and p[3] in ("u5", "kt") and p[4] in ("0", "1") and int(p[2]) >= 1, l
        assert p[5].rstrip().endswith(f"--out {KQ}/{p[0]}"), l
    cu = open(V).read().splitlines()
    bo = set(bo); assert all(i.startswith("u5_tv_") and _chua_khoi(i) for i in bo), [i for i in bo if not _chua_khoi(i)]
    cu = [l for l in cu if l.split()[0] not in bo]; co = {l.split()[0] for l in cu if l.strip()}
    assert not (co & set(ids)), co & set(ids)
    w1 = max([n for n, l in enumerate(cu) if l.split()[0].startswith(("u5_c_", "u5_x_", "u5_r_"))], default=-1)
    tv = [n for n, l in enumerate(cu) if n > w1 and l.split()[0].startswith("u5_tv_")]
    vt = next((n for n in tv if _chua_khoi(cu[n].split()[0])), (tv[-1] + 1) if tv else w1 + 1)
    ts = time.strftime("%Y%m%d_%H%M%S"); shutil.copy2(V, f"{U}/sao_luu/viec.txt.{ts}"); shutil.copy2(M, f"{U}/sao_luu/MD5_TDG.txt.{ts}")
    ra = cu[:vt] + them + cu[vt:]
    open(V + ".moi", "w").write("\n".join(ra) + "\n"); os.replace(V + ".moi", V)
    ml = open(M).read().splitlines(); ml = [f"{md5(V)}  {V}" if l.endswith(f"  {V}") else l for l in ml]
    for f in MA_MD5:
        cu_ = [l for l in ml if l.endswith(f"  {f}")]
        if cu_: assert cu_[0].split()[0] == md5(f), f"md5 {f} doi so voi MD5_TDG — khong chen"
        elif os.path.exists(f): ml.append(f"{md5(f)}  {f}")
    open(M + ".moi", "w").write("\n".join(ml) + "\n"); os.replace(M + ".moi", M)
    open(f"{U}/nhat_ky_hang.txt", "a").write(f"{time.strftime('%F %T')} {ghi} chen {len(them)} (truoc dong {vt + 1}) bo {len(bo)} sao_luu {ts} viec.txt {md5(V)}\n")
    return len(them), len(bo)
def xong(i): return os.path.exists(f"{R}/{O}/xong/{i}")
def loi(i): return os.path.exists(f"{R}/{O}/loi/{i}")
def ket(i): return xong(i) or loi(i)
def doc(i, dat=False):
    f = f"{R}/{KQ}/{i}/" + ("thong_tin_dat.json" if dat else "thong_tin.json")
    return json.load(open(f)) if os.path.exists(f) else None
