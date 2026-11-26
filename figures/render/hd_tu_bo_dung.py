# -*- coding: utf-8 -*-
"""hd_tu_bo_dung.py (24/11/2026) — chuyển HD_i của BỘ DỰNG (thong_tin.json → ban_u5b.hd_truoc_sau, [đỉnh sau→mặt trước, đỉnh trước→mặt sau],
cùng định nghĩa với n2_tv.py; khớp hd_n2_dv.py tới 5e-5 µm trên 20 việc k08) sang dạng json của ve_hd_blender.py. Dùng cho việc KHÔNG giữ mesh trước
chiếu (N14). Dùng: python hd_tu_bo_dung.py <ra_thu_muc> <thu_muc_viec>..."""
import sys, os, json
ra = sys.argv[1]
for d in sys.argv[2:]:
    j = os.path.basename(os.path.normpath(d)); t = json.load(open(f"{d}/thong_tin.json")); h = t["ban_u5b"]["hd_truoc_sau"]
    bp = {k: {"hd": v, "hd_max": max(v)} for k, v in h.items() if not k.startswith("_") and v}
    dat = [e["hd_max"] for k, e in bp.items() if k in (t.get("tu_the") or {})]
    json.dump({"viec": j, "nguon": "thong_tin.json ban_u5b.hd_truoc_sau", "hd_max_dat": max(dat), "bo_phan": bp}, open(f"{ra}/hd_{j}.json", "w"))
    print(j, "hd_max_dat", max(dat), "n", len(bp))
