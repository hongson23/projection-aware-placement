"""so_sanh.py A B [--ra out.json] — so kết quả cellPACK bản gốc (A) và bản GPU (B) cùng đầu vào, cùng seed.
TRUNG khi: cùng tập + thứ tự bộ phận đặt, vị trí t và ma trận xoay R trùng từng bit (số thực JSON ghi bằng repr => khôi phục đúng bit),
cùng số lần loại; và (nếu có van_tay.json) cùng dấu vân tay lưới sau dựng, sau đặt, và trạng thái bộ sinh ngẫu nhiên sau đặt."""
import json, sys, struct
A, B = sys.argv[1], sys.argv[2]
ra = sys.argv[sys.argv.index("--ra") + 1] if "--ra" in sys.argv else None
ta = json.load(open(f"{A}/tu_the_cp.json")); tb = json.load(open(f"{B}/tu_the_cp.json"))
loi = []
if ta["md5_dau_vao"] != tb["md5_dau_vao"]: loi.append("md5 dau vao khac")
if ta["seed_cp"] != tb["seed_cp"]: loi.append("seed khac")
ka, kb = list(ta["dat"]), list(tb["dat"])
if ka != kb: loi.append(f"tap/thu tu bo phan dat khac: {len(ka)} vs {len(kb)}; chi A {sorted(set(ka)-set(kb))[:5]} chi B {sorted(set(kb)-set(ka))[:5]}")
bit = lambda v: [struct.pack("<d", float(x)) for x in (v if not isinstance(v[0], list) else sum(v, []))]
lech = 0; lech_max = 0.0
for k in ka:
    if k not in tb["dat"]: continue
    a, b = ta["dat"][k], tb["dat"][k]
    if a["ten_cp"] != b["ten_cp"]: loi.append(f"ten_cp {k}")
    if bit(a["t"]) != bit(b["t"]) or bit(a["R"]) != bit(b["R"]):
        lech += 1; lech_max = max(lech_max, max(abs(x - y) for x, y in zip(a["t"], b["t"])))
if lech: loi.append(f"{lech} bo phan lech tu the (max |dt| {lech_max:.3e})")
for f in ("so_dat_cp", "so_lan_loai", "so_diem_luoi", "spacing_luoi"):
    if ta.get(f) != tb.get(f): loi.append(f"{f}: {ta.get(f)} vs {tb.get(f)}")
vt = {}
try:
    va = json.load(open(f"{A}/van_tay.json")); vb = json.load(open(f"{B}/van_tay.json"))
    for nhom in ("luoi", "sau_dat"):
        for k in sorted(set(va.get(nhom, {})) | set(vb.get(nhom, {}))):
            if va.get(nhom, {}).get(k) != vb.get(nhom, {}).get(k): loi.append(f"van tay {nhom}.{k} khac")
    for k in ("dat_sha", "np_rng_sha", "py_rng_sha"):
        if va.get(k) != vb.get(k): loi.append(f"van tay {k} khac")
    vt = {"giay_luoi": [va.get("giay_luoi"), vb.get("giay_luoi")], "giay_dat": [va.get("giay_dat"), vb.get("giay_dat")],
          "file": [va.get("cellpack_file"), vb.get("cellpack_file")]}
except FileNotFoundError:
    loi.append("thieu van_tay.json (chi so tu_the)")
kq = {"A": A, "B": B, "ket_luan": "TRUNG" if not loi else "KHAC", "loi": loi, "so_dat": [len(ka), len(kb)], "so_yeu_cau": ta.get("so_yeu_cau"),
      "giay": {"A": [ta.get("giay_luoi"), ta.get("giay_dat")], "B": [tb.get("giay_luoi"), tb.get("giay_dat")]}, **vt}
print(json.dumps(kq, ensure_ascii=False))
if ra: json.dump(kq, open(ra, "w"), indent=1, ensure_ascii=False)
sys.exit(0 if not loi else 3)
