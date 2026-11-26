# -*- coding: utf-8 -*-
"""hd_tu_n2.py (24/11/2026) — chuyển n2.json của thước N2 thực vật (n2_tv.py md5 14b46b12…: bo_phan[k].hd = [đỉnh sau→mặt trước, đỉnh trước→mặt sau];
vung = ER, không bào từ bộ dựng) sang json của ve_hd_blender.py. Dùng: python hd_tu_n2.py <ra_thu_muc> <thu_muc_viec>..."""
import sys, os, json
ra = sys.argv[1]
for d in sys.argv[2:]:
    j = os.path.basename(os.path.normpath(d)); n = json.load(open(f"{d}/n2.json"))
    bp = {k: {"hd": e["hd"], "hd_max": max(e["hd"])} for k, e in n["bo_phan"].items() if e.get("hd")}
    for k, e in (n.get("vung") or {}).items():
        if e.get("hd"): bp[k] = {"hd": e["hd"], "hd_max": max(e["hd"]), "vung": True}
    json.dump({"viec": j, "nguon": "n2.json", "hd_max_ran": n["hd_max_ran"], "bo_phan": bp}, open(f"{ra}/hd_{j}.json", "w"))
    print(j, "hd_max_ran", n["hd_max_ran"], "n", len(bp))
