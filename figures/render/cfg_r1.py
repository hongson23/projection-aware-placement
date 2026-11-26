# -*- coding: utf-8 -*-
"""cfg_r1.py — sinh cấu hình ghép Hình R1 (k08, seed 6512) từ du_lieu/hd_*.json; nhãn ô phóng = bào quan nhỏ có HD_i lớn nhất trong hộp."""
import json, numpy as np
S = 6512; ZR = 1.2; VMAX = 0.2771
Z = [("ti_the_390", [0.57, 5.83, 6.67]), ("ti_the_295", [-4.84, 5.72, 4.59]), ("ti_the_11", [5.56, 2.75, -6.47])]
cot = []
for a, ten in (("P8", "P8"), ("U5b", "PA-PLACE")):
    H = json.load(open(f"du_lieu/hd_u5_r_{a}_k08_{S}.json")); bp = H["bo_phan"]; nh = []
    for i, (n, c) in enumerate(Z):
        trong = [(e["hd_max"], k) for k, e in bp.items() if e.get("hd") is not None and k.startswith(("ti_the", "lysosome"))
                 and np.all(np.abs(np.array(e["c"]) - c) <= ZR)]
        m = max(trong)
        nh.append(f"{i + 1}  box max HDᵢ = {m[0]:.3f} µm")
    cot.append({"tieu_de": ten, "phu_de": f"cell HD max = {H['hd_max_tru_mang_tang2']:.3f} µm", "cat": f"r1/cat_{a}.png",
                "zoom": [f"r1/zoom{i}_{a}.png" for i in (1, 2, 3)], "nhan_zoom": nh, "diem": [c for _, c in Z]})
json.dump({"ra": "R1_k08_P8_vs_PAPLACE.png", "vmax": VMAX, "lut": "du_lieu/lut_magma_r.json", "elev": 15, "ortho": 22, "res": 1600, "cw": 1500,
           "nhan_thang": "HDᵢ (µm)", "hang": ["Half section", "Zoom"], "cot": cot,
           "ghi_chu": f"k08, seed {S}. Zoom: solid = final mesh, blue wireframe = mesh before projection. Same camera and colour scale in both columns."},
          open("du_lieu/cfg_r1.json", "w"), ensure_ascii=False, indent=1)
print(json.dumps(cot, ensure_ascii=False)[:800])
