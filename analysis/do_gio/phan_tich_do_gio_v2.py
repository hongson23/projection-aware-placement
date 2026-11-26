# -*- coding: utf-8 -*-
"""phan_tich_do_gio.py BẢN 3 (23/11/2026) — bảng giờ chạy của bài CAD, CHỈ từ lượt đo giờ độc quyền (bo_dung_chung/do_gio).
Dùng: python phan_tich_do_gio.py  ⇒ BANG_GIO_DOC_QUYEN_v2.md, bang_k08_v2.csv, bang_tv_v2.csv, doc_quyen_v2.csv
Pha (giây, đồng hồ tường; mã bộ dựng không đổi):
  dat = bao_cao_bo_tri.giay_dat | khoi = giay_bo_tri − dat | chieu = giay_chieu | hd = Σ ban_u5b.hd_*._giay (đánh giá, ngoài phương pháp)
  cb_chieu = t(dòng SO_HANG) − (t00 + giay_bo_tri), t00 ≈ t(dòng V12/THOI_GIAN cuối) − giay_tong   (raster trường + dựng hàng ràng buộc; chỉ bộ dựng động vật)
  sau = giay_tong − giay_bo_tri − cb_chieu − chieu − hd (mesh + đo) | pp = giay_tong − hd | khoi_dong = tường /usr/bin/time − giay_tong
  Thực vật: thêm thoi_gian_buoc_giay (bộ dựng tự ghi từng bước).
Chỉ dùng lượt DOC_QUYEN (giam_doc_quyen.py bản 3: GPU không có tiến trình lạ, không có tính toán của người dùng khác trên vùng lõi, mọi tiến trình ta ghim).
"""
import os, sys, json, csv, statistics as st
R = os.path.expanduser("~/Vision_team_data_189/sonhh20/ch3_tebao_dep"); os.chdir(R)
DG = os.environ.get("DG_THU_MUC", "bo_dung_chung/do_gio"); VL = f"{DG}/viec_do_gio.txt"
rows = {}
for l in open(VL):
    p = l.rstrip("\n").split("\t")
    if len(p) >= 9: rows[p[0]] = dict(zip(["id", "pha", "lan", "loai", "gb", "nguon", "nhom", "uoc", "lenh"], p))
nh = lambda i: rows[i]["nhom"].split("|")
def xong(i):
    f = f"{DG}/xong/{i}"
    if not os.path.exists(f): return None
    t = open(f).read().split(); d = {"run": t[0]}
    for x in t[1:]: k, _, v = x.partition("="); d[k] = v
    return d
def jl(f):
    try: return json.load(open(f))
    except Exception: return None
def tgian(run):
    o = {}
    try:
        for l in open(f"{DG}/log/{run}.time"):
            if "Maximum resident set size" in l: o["rss_gb"] = int(l.split(":")[-1]) / 1048576
            if "Elapsed (wall clock)" in l:
                v = l.rsplit(" ", 1)[-1].strip().split(":"); o["wall"] = sum(float(x) * 60 ** i for i, x in enumerate(reversed(v)))
    except OSError: pass
    return o
def gio_dong(run):
    o = {}
    try:
        for l in open(f"{DG}/log/{run}.log", errors="replace"):
            p = l.split(" ", 2)
            if len(p) < 2: continue
            try: t = float(p[0])
            except ValueError: continue
            if p[1].startswith("SO_HANG") and "SO_HANG" not in o: o["SO_HANG"] = t
            if p[1].startswith("V12") or p[1].startswith("THOI_GIAN"): o["KET"] = t
    except OSError: pass
    return o
def hd_giay(t):
    return sum(v["_giay"] for k, v in (t.get("ban_u5b") or {}).items() if k.startswith("hd") and isinstance(v, dict) and isinstance(v.get("_giay"), (int, float)))
def igl_dung(t):
    r = []
    def tim(d):
        if isinstance(d, dict):
            if {"dung", "ly_do", "so_lan"} <= set(d): r.append(d.get("dung"))
            for v in d.values(): tim(v)
        elif isinstance(d, list):
            for v in d: tim(v)
    tim(t); return r
def doc(i):
    x = xong(i)
    if not x: return None
    run = x["run"]; t = jl(f"{DG}/kq/{run}/thong_tin.json"); dq = jl(f"{DG}/doc_quyen/{run}.json") or {}
    if not t: return None
    dat = (t.get("bao_cao_bo_tri") or {}).get("giay_dat"); hd = hd_giay(t)
    r = {"id": i, "run": run, "doc_quyen": dq.get("ket_luan"), "vuot_nen": [k for k, v in (dq.get("vuot_nen") or {}).items() if v], "so_bit": x.get("so_bit"),
         "dat": dat, "khoi": t["giay_bo_tri"] - dat if dat is not None else None, "chieu": t["giay_chieu"], "hd": hd, "cb_chieu": None,
         "sau": t["giay_tong"] - t["giay_bo_tri"] - t["giay_chieu"] - hd, "pp": t["giay_tong"] - hd, "tong": t["giay_tong"],
         "gpu_torch_gb": t.get("bo_nho_gb"), "gpu_smi_gb": (dq.get("gpu_ta_dinh_mib") or 0) / 1024, "tan_mhz": dq.get("tan_so_L_mhz_tb"), "nhiet": dq.get("nhiet_max"),
         "he_thong_loi": (dq.get("la_theo_chu_loi") or {}).get("he_thong"), "rss_gb": None, "wall": None, "buoc": t.get("thoi_gian_buoc_giay"), "igl": igl_dung(t)}
    r.update(tgian(run))
    gd = gio_dong(run)
    if "SO_HANG" in gd and "KET" in gd:
        cb = gd["SO_HANG"] - (gd["KET"] - t["giay_tong"] + t["giay_bo_tri"])
        if -2 < cb and r["sau"] - cb > -1: r["cb_chieu"] = max(0.0, cb); r["sau"] -= r["cb_chieu"]
    r["khoi_dong"] = r["wall"] - t["giay_tong"] if r["wall"] else None
    return r
med = lambda v: st.median(v) if v else None
def fm(v, d=1): return "—" if v is None else (f"{v:.{d}f}" if isinstance(v, (int, float)) else str(v))
def tom(v, d=1):
    v = [x for x in v if isinstance(x, (int, float))]
    return "—" if not v else (fm(v[0], d) if len(v) == 1 else f"{fm(med(v), d)} [{fm(min(v), d)}–{fm(max(v), d)}]")
D = {i: doc(i) for i in rows}
import glob as _g
SAI_SO = {}
for _f in _g.glob(f"{DG}/so_bit/*.json"):
    _d = json.load(open(_f)); _m = 0.0
    for _k, _a, _b in _d.get("khac", []):
        try:
            _a = float(_a); _b = float(_b); _m = max(_m, abs(_a - _b) / max(abs(_a), abs(_b), 1e-300))
        except Exception:
            _m = float("inf")
    SAI_SO[os.path.basename(_f)[:-5]] = _m
sach = lambda r: r is not None and r["doc_quyen"] == "DOC_QUYEN"
out = []; W = out.append
mt = open(f"{DG}/moi_truong_do.txt").read() if os.path.exists(f"{DG}/moi_truong_do.txt") else ""
W("# Bảng giờ chạy — lượt đo giờ ĐỘC QUYỀN (Bài 1 ch3, bản 3)\n")
W(f"Bản v2 (25/11): mọi lượt DOC_QUYEN vào trung vị. So bit: lệch chỉ ở sai số làm tròn (≤ 1e-12): {sorted(i for i, m in SAI_SO.items() if m <= 1e-12)}; PHÂN KỲ (lệch lớn): {sorted((i, round(m, 4)) for i, m in SAI_SO.items() if m > 1e-12)}.")
W(f"(Dòng gốc bản 3:) Chỉ lượt DOC_QUYEN và không KHÁC bit. KHÁC bit (loại, phải giải thích): {[i for i, r in D.items() if r and r['so_bit'] == 'KHAC']}. Thiếu/lỗi: {[i for i, r in D.items() if r is None]}; không độc quyền (giữ, loại khỏi bảng): {[i for i, r in D.items() if r and not sach(r)]}.\n")
COT = ["khoi_dong", "khoi", "dat", "cb_chieu", "chieu", "sau", "pp", "hd", "wall", "gpu_torch_gb", "gpu_smi_gb", "rss_gb"]
TEN = {"khoi_dong": "khởi động", "khoi": "khởi tạo", "dat": "đặt", "cb_chieu": "chuẩn bị PROJECT", "chieu": "PROJECT", "sau": "mesh + đo", "pp": "tổng phương pháp",
       "hd": "HD (đánh giá)", "wall": "tường", "gpu_torch_gb": "đỉnh GPU torch (GB)", "gpu_smi_gb": "đỉnh GPU nvidia-smi (GB)", "rss_gb": "đỉnh RSS (GB)"}
K = [r for i, r in D.items() if sach(r) and rows[i]["pha"] == "P1" and nh(i)[0] == "k08"]
W("## k08 (D5 khe 0,8), 3 seed × 1 lượt, hai làn song song — trung vị [min–max] trên seed\n")
W("| nhánh | n | " + " | ".join(TEN[c] for c in COT) + " | pp / P8 (theo seed) | so bit |\n|" + "---|" * (len(COT) + 4))
cs = []
for a in ("U5b", "P4", "P8", "CXe", "CX2e", "CX2b"):
    rs = [r for r in K if nh(r["id"])[1] == a]
    if not rs: continue
    ts = []
    for r in rs:
        c = [x for x in K if nh(x["id"])[1] == "P8" and nh(x["id"])[3] == nh(r["id"])[3]]
        if c: ts.append(r["pp"] / c[0]["pp"])
    W(f"| {a} | {len(rs)} | " + " | ".join(tom([r.get(c) for r in rs], 2 if "gb" in c else 1) for c in COT) + f" | {tom(ts, 3)} | {sum(r['so_bit'] == 'TRUNG' for r in rs)}/{len(rs)} TRUNG |")
    cs += [{"nhanh": a, "seed": nh(r["id"])[3], "lan": rows[r["id"]]["lan"], **{c: r.get(c) for c in COT + ["tong", "so_bit", "tan_mhz", "he_thong_loi"]}, "vuot_nen": "+".join(r["vuot_nen"])} for r in rs]
W("\nGhi chú: CX* chạy bộ dựng v122u5x, U5b/P4/P8 chạy v122u5e (cùng phép tính khi cờ tắt); ở lượt cũ U5b trên u5x chậm hơn ≈ 3,5 % (u5_x_U5bcal 661/773 s so với 637/747 s) — tỉ số CX/P8 mang độ lệch bộ dựng này. Tầng 1 (U5b/P4/P8) mỗi seed trọn một làn; CX seed 6513 ở làn khác P8 6513.")
if cs:
    with open(f"{DG}/bang_k08_v2.csv", "w", newline="") as f: w = csv.DictWriter(f, list(cs[0])); w.writeheader(); w.writerows(cs)
T = [r for i, r in D.items() if r and rows[i]["pha"] == "P1" and nh(i)[0] == "tv"]
if T:
    W("\n## Thực vật TV-real seed 7311, 1 lượt — giờ từng bước do bộ dựng tự ghi\n")
    W("| nhánh | độc quyền | khởi tạo | đặt | PROJECT | tổng pp | HD | tường | đỉnh GPU torch (GB) | đỉnh RSS (GB) | igl phụ dùng | so bit | các bước (s) |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    cs = []
    for r in sorted(T, key=lambda r: r["id"]):
        b = r["buoc"] or {}
        W(f"| {nh(r['id'])[1]} | {r['doc_quyen']} | {fm(r['khoi'])} | {fm(r['dat'])} | {fm(r['chieu'])} | {fm(r['pp'])} | {fm(r['hd'])} | {fm(r['wall'])} | "
          f"{fm(r['gpu_torch_gb'], 2)} | {fm(r['rss_gb'], 2)} | {r['igl']} | {r['so_bit']} | " + ", ".join(f"{k} {fm(v)}" for k, v in b.items() if isinstance(v, (int, float))) + " |")
        cs.append({"nhanh": nh(r["id"])[1], "doc_quyen": r["doc_quyen"], **{c: r.get(c) for c in COT + ["tong", "so_bit"]}, "igl": r["igl"], **{f"buoc_{k}": v for k, v in b.items()}})
    keys = sorted({k for c in cs for k in c}); [c.setdefault(k, None) for c in cs for k in keys]
    with open(f"{DG}/bang_tv_v2.csv", "w", newline="") as f: w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(cs)
    W("\nCột igl phụ: cờ `dung` của tiến trình igl phụ; phải True ở mọi nhánh để giờ so được với nhau (lượt chính thực vật có khi chạy KHÔNG có igl phụ ⇒ bước D_igl chậm ≈ 2–3×; lượt đo luôn dùng igl phụ).")
    W("Tổng pp của thực vật = giay_tong − Σ hd_*: CÒN GỒM bước `thong_ke_va_D_igl` = D tới mesh sau chiếu (phương pháp: dùng cho vùng phụ thuộc) + D0 tới mesh TRƯỚC chiếu (chỉ đánh giá, khi có --hd_vung_truoc) + thể tích/IoU (đánh giá) — không tách được nếu không sửa mã; báo bước này riêng, và tổng pp là CẬN TRÊN của giờ phương pháp. Lượt --tv_thuan bỏ phần đánh giá.")
a = D.get("don_U5b_6511_g0"); b = D.get("k08_U5b_6511")
W("\n## Nhiễu chéo hai làn\n")
W(f"- U5b k08 6511 chạy cặp (GPU0) / một mình (GPU0): {b['pp'] / a['pp']:.4f}; theo pha: " + ", ".join(f"{k} {b[k] / a[k]:.4f}" for k in ("dat", "chieu", "sau") if a.get(k) and b.get(k))
  if a and b and sach(a) and sach(b) else "- Không đo (không bật --don hoặc thiếu lượt độc quyền): hai việc trên hai GPU dùng chung L3/bộ nhớ của máy; ghi là giới hạn.")
W("\n## Độc quyền (giam_doc_quyen.py)\n")
dqr = []
for i in rows:
    x = xong(i)
    if not x: continue
    j = jl(f"{DG}/doc_quyen/{x['run']}.json") or {}
    dqr.append({"id": i, "run": x["run"], "ket": j.get("ket_luan"), "nguoi_dung_loi": (j.get("la_theo_chu_loi") or {}).get("nguoi_dung"),
                "he_thong_loi": (j.get("la_theo_chu_loi") or {}).get("he_thong"), **(j.get("ti_le") or {}), "vuot_nen": "+".join(k for k, v in (j.get("vuot_nen") or {}).items() if v),
                "gpu_la_max": j.get("gpu_la_max"), "gpu_kia": j.get("gpu_kia_so_tien_trinh_max"), "tan": j.get("tan_so_L_mhz_tb"), "nhiet": j.get("nhiet_max"), "so_bit": x.get("so_bit")})
if dqr:
    with open(f"{DG}/doc_quyen_v2.csv", "w", newline="") as f: w = csv.DictWriter(f, sorted({k for d in dqr for k in d})); w.writeheader(); w.writerows(dqr)
    W(f"- {sum(d['ket'] == 'DOC_QUYEN' for d in dqr)}/{len(dqr)} lượt DOC_QUYEN; chạy lại (_l2): {sum(d['run'].endswith('_l2') for d in dqr)}; GPU lạ max {max((d['gpu_la_max'] or 0) for d in dqr)}; "
      f"tính toán người dùng khác max {max((d['nguoi_dung_loi'] or 0) for d in dqr):.4f} lõi.")
    W(f"- Tải hệ thống (falco, tailscaled, kworker…; được phép) trên vùng lõi: trung vị {tom([d['he_thong_loi'] for d in dqr], 3)} lõi; lượt vượt nền p95: {[d['id'] for d in dqr if d['vuot_nen']]}.")
    W(f"- Tần số lõi việc TB {tom([d['tan'] for d in dqr], 0)} MHz; nhiệt GPU max {max((d['nhiet'] or 0) for d in dqr)} °C.")
sbs = [(i, xong(i).get("so_bit")) for i in rows if rows[i]["nguon"] != "-" and xong(i)]
W(f"\n## So từng bit với lượt gốc\n\n{sum(s == 'TRUNG' for _, s in sbs)}/{len(sbs)} TRUNG; KHÁC: {[i for i, s in sbs if s == 'KHAC']}; không có gốc: {[i for i, s in sbs if s == 'KHONG_CO_GOC']}. "
  f"Thực vật không có gốc để so: {[i for i in rows if rows[i]['nguon'] == '-']}.")
if mt: W("\n## Môi trường đo\n\n```\n" + "\n".join(l for l in mt.splitlines() if any(k in l for k in ("Model name", "CPU(s):", "NVIDIA", "python", "driver", "Thread", "ondemand", "performance")))[:3000] + "\n```")
open(f"{DG}/BANG_GIO_DOC_QUYEN_v2.md", "w").write("\n".join(out) + "\n"); print("\n".join(out))
