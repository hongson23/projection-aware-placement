# -*- coding: utf-8 -*-
"""tong_kin.py (24/11/2026) — gộp kiem_kin/<viec>.json (kiem_kin.py) thành KIEM_KIN.md: theo nhóm (k08 P8/U5b, N14 d5 theo nhánh, thực vật theo
mức × nhánh), mỗi mức (sau chiếu / sau vệ sinh / CUỐI = sau sửa mesh đồng nhất): số bộ phận kín (trimesh is_watertight, process=False) / tổng,
tổng cạnh biên, tổng cạnh không đa tạp, số bộ phận có cạnh biên / có cạnh không đa tạp; mức cuối thêm số bộ phận không cùng chiều (is_winding_consistent).
Liệt kê mọi bộ phận không kín ở mức cuối (việc, bộ phận, cạnh biên, cạnh không đa tạp). Dùng: python tong_kin.py <kiem_kin_dir> <ra.md>"""
import sys, os, json, glob, re, collections
D, RA = sys.argv[1], sys.argv[2]
def nhom(j):
    m = re.match(r"u5_r_(P8|U5b)_k08_\d+$", j)
    if m: return f"động vật k08 | {'P8' if m.group(1) == 'P8' else 'PA-PLACE (U5b)'}"
    m = re.match(r"u5_n14d5_d5_\d+_(\w+?)_\d+$", j)
    if m: return f"N14 d5 (Allen) | {dict(U5b='PA-PLACE (U5b)', td='PLACE (td)', nn1x='thử–loại 1× (nn1x)', unc8='SA c08 (unc8)').get(m.group(1), m.group(1))}"
    m = re.match(r"u5_tv_c_(\w+?)_(u5ktt|td|nn1xm|nn1x|un1x)_\d+$", j)
    if m: return f"thực vật TV-{m.group(1)} | {dict(u5ktt='PA-PLACE+F (u5ktt)', td='PLACE (td)', nn1x='thử–loại (nn1x)', nn1xm='thử–loại J mềm (nn1xm)', un1x='SA (un1x)')[m.group(2)]}"
    return "khác | " + j
R = [json.load(open(f)) for f in sorted(glob.glob(f"{D}/*.json"))]
G = collections.defaultdict(list); loi = []; bq = [r for r in R if "bo_qua" in r]; R = [r for r in R if "bo_qua" not in r]
for r in R:
    if "cuoi" not in r: loi.append(r); continue
    G[nhom(r["viec"])].append(r)
L = ["# KIEM_KIN — kiểm kín (watertight) theo bộ phận trên mesh cuối", "",
     "Trạng thái: **CHƯA rà** soát chéo (số tự động từ kiem_kin.py; chưa có agent độc lập tính lại).", "",
     "Mesh cuối = `mesh_vo_bao.npz` (sau chiếu) → `bo_kiem_v3/ve_sinh_luoi.py` → `bo_kiem_v3/sua_mesh_dn.py --vat <dong_vat|thuc_vat>`, chạy y như bước 1–2 "
     "của `cham_lo_vsb.sh` trên bản sao ở `/tmp/hinh_cad_$USER` (md5 hai script kiểm theo VE_SINH_LUOI_MD5.txt / SUA_MESH_MD5.txt; "
     f"ve_sinh_luoi {R[0]['md5']['ve_sinh_luoi.py'] if R else '?'}, sua_mesh_dn {R[0]['md5']['sua_mesh_dn.py'] if R else '?'}). "
     "Nén zip deflate thay lzma (không đổi mảng).",
     "Mỗi bộ phận = một cặp `<tên>__v/__f` (gồm màng, vách, nhân, ER, không bào, mọi bào quan; tầng 2 kể riêng từng túi/trung tử).",
     "Kín = `trimesh.Trimesh(V, F, process=False).is_watertight`; cạnh biên = cạnh không hướng thuộc đúng 1 mặt; cạnh không đa tạp = cạnh thuộc > 2 mặt.",
     "Lưu ý: kín (watertight) chỉ nói về liên thông cạnh; KHÔNG kiểm tự cắt (self-intersection) — số tự cắt có ở n2_thuoc_chuan.py (hop_le.tu_cat) nếu cần.", ""]
for lv, ten in (("cuoi", "CUỐI (sau chiếu + vệ sinh + sửa mesh đồng nhất)"), ("ve_sinh", "sau vệ sinh (chưa sửa)"), ("goc", "sau chiếu (npz gốc của bộ dựng)")):
    L += [f"## Mức {ten}", "", "| nhóm | nhánh | số việc | bộ phận kín / tổng | % kín | bộ phận có cạnh biên | Σ cạnh biên | bộ phận có cạnh không đa tạp | Σ cạnh không đa tạp"
          + (" | bộ phận không cùng chiều |" if lv == "cuoi" else " |"), "|---|---|---|---|---|---|---|---|---" + ("|---|" if lv == "cuoi" else "|")]
    for g in sorted(G):
        rs = G[g]; E = [e for r in rs for e in r[lv].values()]
        kin = sum(e["wt"] for e in E); n = len(E)
        row = (f"| {g.split(' | ')[0]} | {g.split(' | ')[1]} | {len(rs)} | {kin} / {n} | {100 * kin / n:.2f} | {sum(e['canh_bien'] > 0 for e in E)} | "
               f"{sum(e['canh_bien'] for e in E)} | {sum(e['canh_kdt'] > 0 for e in E)} | {sum(e['canh_kdt'] for e in E)}")
        if lv == "cuoi": row += f" | {sum(not e.get('wc', True) for e in E)}"
        L.append(row + " |")
    L.append("")
L += ["## Bộ phận không kín ở mức cuối", ""]
bad = [(r["viec"], k, e) for g in sorted(G) for r in G[g] for k, e in r["cuoi"].items() if not e["wt"]]
if not bad: L.append("Không có.")
else:
    L += ["| việc | bộ phận | số mặt | cạnh biên | cạnh không đa tạp |", "|---|---|---|---|---|"]
    L += [f"| {j} | {k} | {e['nF']} | {e['canh_bien']} | {e['canh_kdt']} |" for j, k, e in bad]
L += ["", "## Bộ phận không kín ở mức sau chiếu (npz gốc), theo nhóm × loại bộ phận", "",
      "Mọi trường hợp đều là cạnh không đa tạp (không có cạnh biên); vệ sinh lưới (gộp đỉnh trùng hệt, bỏ mặt suy biến) đã làm chúng kín.", "",
      "| nhóm | nhánh | loại | số bộ phận không kín | Σ cạnh không đa tạp |", "|---|---|---|---|---|"]
for g in sorted(G):
    c = collections.Counter(); e2 = collections.Counter()
    for r in G[g]:
        for k, e in r["goc"].items():
            if not e["wt"]: l = re.sub(r"_?\d+$", "", k); c[l] += 1; e2[l] += e["canh_kdt"]
    L += [f"| {g.split(' | ')[0]} | {g.split(' | ')[1]} | {l} | {c[l]} | {e2[l]} |" for l in sorted(c)]
L += ["", "## Trạng thái sửa mesh đồng nhất (sua_mesh_dn, theo cặp)", ""]
tt = collections.Counter()
for g in sorted(G):
    for r in G[g]:
        for c, v in ((r.get("sua_mesh") or {}).get("cap") or {}).items():
            tt[(g.split(" | ")[0], c, str(v.get("trang_thai")).split(":")[0])] += 1
L += ["| nhóm | cặp | trạng thái | số việc |", "|---|---|---|---|"] + [f"| {a} | {b} | {c} | {n} |" for (a, b, c), n in sorted(tt.items())]
if bq:
    L += ["", "## Việc không kiểm (npz tạm thời)", "", "Thực vật: chỉ kiểm 6 việc giữ npz lâu dài (TV-real và TV-D5 × u5ktt/td/nn1x, seed 7311). Các việc dưới đây có npz lúc "
          "lập danh sách nhưng bộ điều phối xoá npz sau khi chấm (đã mất u5_tv_c_D1_nn1x(m)_7313 giữa chừng lượt đầu), nên không đưa vào:", ""] + [f"- {r['viec']}" for r in bq]
if loi:
    L += ["", "## Việc lỗi", ""] + [f"- {r['viec']}: ve_sinh rc {r.get('ve_sinh_rc')}, sua rc {r.get('sua_rc')}: `{(r.get('loi') or '')[-300:]}`" for r in loi]
open(RA, "w").write("\n".join(L) + "\n"); print("\n".join(L[:60]))
