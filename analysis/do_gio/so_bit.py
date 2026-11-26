# -*- coding: utf-8 -*-
"""so_bit.py (22/11/2026) — so TỪNG BIT kết quả của một lượt đo giờ với kết quả gốc (kiểm "cùng lệnh ⇒ cùng chất lượng").
Dùng: python so_bit.py --loai dung|cpn|pkn|cpv|cpg|pkk --moi <thư mục mới> --goc <thư mục gốc> [--ra f.json]
  dung/cpn/pkn : thong_tin.json (so sâu, số thực so bằng ==) + mọi *.npz có ở GỐC (cùng tên khoá, cùng dtype, cùng giá trị)
  cpv/cpg      : <moi>/cp_chinh|cp_goc/tu_the_cp.json  vs  <goc>/tu_the_cp.json
  pkk          : md5 <moi>/pk/pk_out.xyz vs <goc>/pk/MD5_PK_OUT ; <moi>/tu_the_pk.json vs <goc>/tu_the_pk.json
Bỏ qua (ghi rõ trong kết quả): khoá GIỜ/BỘ NHỚ/GPU (regex BO_QUA) và khối sổ sách tiến trình phụ igl (dict có đủ khoá dung, ly_do, giay, so_lan)
— khối này mô tả MÔI TRƯỜNG (lõi, khoá chỗ, số lõi bận), được báo riêng ở "moi_truong", không tính vào phán quyết.
Khoá CHỈ có ở bản mới (bộ dựng đời sau thêm mục báo cáo) được liệt kê ở "chi_co_o_moi", không tính; khoá có ở gốc mà thiếu ở mới là KHÁC.
Phán quyết: TRUNG (0 khác ngoài phần bỏ qua) | KHAC | KHONG_CO_GOC. Dòng in cuối cùng (stdout) bắt đầu bằng phán quyết.
"""
import os, re, sys, json, hashlib, argparse, math
import numpy as np
ap = argparse.ArgumentParser(); ap.add_argument("--loai", required=True); ap.add_argument("--moi", required=True); ap.add_argument("--goc", required=True); ap.add_argument("--ra")
A = ap.parse_args()
BO_QUA = re.compile(r"(^|_)giay($|_)|^_giay$|^bo_nho|^gpu$|thoi_gian|^cellpack_file$")
CHI_MOI_OK = re.compile(r"^(ban_|u5|du_bao|D_hat)")
def la_igl(d): return isinstance(d, dict) and {"dung", "ly_do", "giay", "so_lan"} <= set(d)
khac, bo, mt, chi_moi = [], [], [], []
def so(a, b, p):
    if isinstance(a, dict) and isinstance(b, dict):
        if la_igl(a) or la_igl(b):
            mt.append({"duong": p, "goc": a, "moi": b}); return
        for k in sorted(set(a) | set(b)):
            q = f"{p}.{k}"
            if BO_QUA.search(str(k)): bo.append(q); continue
            if k not in a and CHI_MOI_OK.match(str(k)): chi_moi.append(q); continue
            if k not in a: khac.append((q, "chi_co_o_moi_khong_phai_khoa_bao_cao")); continue
            if k not in b: khac.append((q, "thieu_khoa_o_moi")); continue
            so(a[k], b[k], q)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b): khac.append((p, "do_dai", len(a), len(b))); return
        for i, (x, y) in enumerate(zip(a, b)): so(x, y, f"{p}[{i}]")
    else:
        if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b): return
        if a != b or type(a) != type(b): khac.append((p, str(a)[:60], str(b)[:60]))
def so_json(fa, fb):
    if not os.path.exists(fa): return "KHONG_CO_GOC"
    if not os.path.exists(fb): khac.append((fb, "khong_co_tep_moi")); return "KHAC"
    so(json.load(open(fa)), json.load(open(fb)), os.path.basename(fa)); return None
mang = {}
def so_npz(dg, dm):
    for f in sorted(os.listdir(dg)):
        if not f.endswith(".npz"): continue
        fm = os.path.join(dm, f)
        if not os.path.exists(fm): khac.append((f, "khong_co_npz_moi")); mang[f] = "THIEU"; continue
        a, b = np.load(os.path.join(dg, f)), np.load(fm); ka, kb = set(a.files), set(b.files)
        if ka != kb: khac.append((f, "khoa_khac", sorted(ka ^ kb)[:5])); mang[f] = "KHOA_KHAC"; continue
        sai = [k for k in sorted(ka) if a[k].dtype != b[k].dtype or a[k].shape != b[k].shape or not np.array_equal(a[k], b[k])]
        for k in sai[:10]: khac.append((f"{f}:{k}", "mang_khac"))
        mang[f] = {"so_mang": len(ka), "khac": len(sai)}
ket = None
if A.loai in ("dung", "cpn", "pkn"):
    ket = so_json(os.path.join(A.goc, "thong_tin.json"), os.path.join(A.moi, "thong_tin.json"))
    if ket is None and os.path.isdir(A.goc): so_npz(A.goc, A.moi)
elif A.loai in ("cpv", "cpg"):
    sub = "cp_chinh" if A.loai == "cpv" else "cp_goc"
    ket = so_json(os.path.join(A.goc, "tu_the_cp.json"), os.path.join(A.moi, sub, "tu_the_cp.json"))
elif A.loai == "pkk":
    fg = os.path.join(A.goc, "pk", "MD5_PK_OUT"); fm = os.path.join(A.moi, "pk", "pk_out.xyz")
    if not os.path.exists(fg): ket = "KHONG_CO_GOC"
    else:
        m0 = open(fg).read().split()[0]; m1 = hashlib.md5(open(fm, "rb").read()).hexdigest() if os.path.exists(fm) else None
        mang["pk_out.xyz"] = {"goc": m0, "moi": m1}
        if m0 != m1: khac.append(("pk_out.xyz", m0, m1))
        ket = so_json(os.path.join(A.goc, "tu_the_pk.json"), os.path.join(A.moi, "tu_the_pk.json"))
else: sys.exit("loai sai")
if ket is None: ket = "TRUNG" if not khac else "KHAC"
r = {"ket": ket, "loai": A.loai, "moi": A.moi, "goc": A.goc, "so_khac": len(khac), "khac": [list(map(str, x)) for x in khac[:30]],
     "so_bo_qua": len(bo), "bo_qua_mau": bo[:40], "chi_co_o_moi": chi_moi, "moi_truong": mt, "mang": mang}
s = json.dumps(r, ensure_ascii=False, indent=1)
if A.ra: open(A.ra, "w").write(s)
print(ket, len(khac), "khac;", len(bo), "bo qua;", len(mt), "khoi moi truong")
