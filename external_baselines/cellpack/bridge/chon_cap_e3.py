"""chon_cap_v122.py <thu_v122> — chọn cặp thử tương đương gốc ↔ v4 trên đầu vào v12.2 theo luật ĐĂNG KÝ TRƯỚC của rà soát độc lập
RA_SOAT_V4.md §7.2 (luật của rs_chon_cap.py, chỉ đổi cách đọc tên thư mục): mỗi ô D4A2, D5, trên các lượt v4 cấu hình chính (cp_chinh):
 (i) lượt nhiều bộ phận nhất (hoà: nhiều lần loại hơn); (ii) lượt nhiều lần loại nhất; (iii) lượt 'thường': seed nhỏ nhất ≥ 6205, k nhỏ nhất;
 (iv) MỌI lượt có bộ đếm ≠ 0 (v3_va_cham_du_phong, va_cham_cpu_du_phong, chua_mo_ho, chua_broken, v3_chua_mo_ho, gpu_kiem, cellpack_file ngoài ban/v4),
      tối đa 2 mỗi ô; hơn 2 ⇒ in DUNG_RA_LAI (hàng dừng nhập, PI rà).
In ra stdout các dòng 'G <ô> <seed> <k> <lý do>' (trùng thì gộp lý do); stderr: bảng mọi lượt. Chỉ đọc JSON."""
import json, glob, os, re, sys
DEM = ["v3_va_cham_du_phong", "va_cham_cpu_du_phong", "chua_mo_ho", "chua_broken", "v3_chua_mo_ho"]
rows = []
for d in sorted(glob.glob(f"{sys.argv[1]}/*_*_k*/cp_chinh")):
    m = re.fullmatch(r"(D4A2|D5)_(\d+)_k(\d+)", os.path.basename(os.path.dirname(d)))
    if not m or not (os.path.exists(f"{d}/van_tay.json") and os.path.exists(f"{d}/tu_the_cp.json")): continue
    v = json.load(open(f"{d}/van_tay.json")); t = json.load(open(f"{d}/tu_the_cp.json")); st = v.get("gpu_thong_ke") or {}
    bad = {c: st.get(c) for c in DEM if st.get(c)}
    if v.get("gpu_kiem"): bad["gpu_kiem"] = True
    if "ban/v4/" not in str(v.get("cellpack_file")): bad["cellpack_file"] = v.get("cellpack_file")
    rows.append(dict(o=m.group(1), s=int(m.group(2)), k=int(m.group(3)), n=t.get("so_yeu_cau") or 0, loai=t.get("so_lan_loai") or 0, bad=bad))
    print(f"{m.group(0):16s} n={rows[-1]['n']} loai={rows[-1]['loai']} {'BO_DEM ' + str(bad) if bad else 'bo_dem=0'}", file=sys.stderr)
out = {}
for o in ("D5", "D4A2"):
    R_ = [r for r in rows if r["o"] == o]
    if not R_: continue
    def them(r, ly): out.setdefault((o, r["s"], r["k"]), []).append(ly)
    them(max(R_, key=lambda r: (r["n"], r["loai"])), "i_nhieu_bo_phan"); them(max(R_, key=lambda r: r["loai"]), "ii_nhieu_loai")
    th = sorted([r for r in R_ if r["s"] >= 6205], key=lambda r: (r["s"], r["k"]))
    if th: them(th[0], "iii_thuong")
    xau = [r for r in R_ if r["bad"]]
    if len(xau) > 2: print(f"DUNG_RA_LAI {o} {len(xau)}")
    for r in xau[:2]: them(r, "iv_bo_dem_khac_0")
for (o, s, k), ly in out.items(): print(f"G {o} {s} {k} {'+'.join(ly)}")
