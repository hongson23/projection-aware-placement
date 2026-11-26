# -*- coding: utf-8 -*-
"""hd_n2_dv.py (24/11/2026, phiên thực nghiệm, cho hình bài CAD) — HD_i theo bộ phận giữa mesh TRƯỚC chiếu và SAU chiếu, dùng ĐÚNG
hàm `mesh` và `kc` của thước N2 bo_dung_chung/tang_dg/u5tv/n2_tv.py (md5 14b46b125905736cb4217fd7d7b27a51; hai hàm được exec từ chính
tệp đó sau khi kiểm md5): kc(P, m) = sqrt(igl.point_mesh_squared_distance(P, V, F)), HD_i = [max kc(đỉnh sau → mặt trước),
max kc(đỉnh trước → mặt sau)]; hai mesh trùng hệt ⇒ [0, 0] (như n2_tv.py). Khác n2_tv.py duy nhất: tên tệp trước chiếu
(động vật: mesh_vo_bao_truoc.npz; thực vật: mesh_truoc_giu.npz) và duyệt MỌI bộ phận có ở cả hai npz (kể cả màng, nhân, ER).
Kiểm chéo: chênh lệch tối đa với thong_tin.json ban_u5b.hd_truoc_sau (bộ dựng) ⇒ 'lech_bo_dung'.
Dùng: python hd_n2_dv.py <ra_thu_muc> <thu_muc_viec>...  ⇒ <ra_thu_muc>/hd_<viec>.json"""
import sys, os, json, hashlib, numpy as np, igl
N2 = os.path.expanduser("~/Vision_team_data_189/sonhh20/ch3_tebao_dep/bo_dung_chung/tang_dg/u5tv/n2_tv.py")
src = open(N2, "rb").read(); assert hashlib.md5(src).hexdigest() == "14b46b125905736cb4217fd7d7b27a51", "n2_tv.py đổi md5"
L = src.decode("utf-8").splitlines(); ns = {"np": np, "igl": igl}
for i, l in enumerate(L):
    if l.startswith("def mesh(") or l.startswith("def kc("):
        blk = [l]; j = i + 1
        while j < len(L) and L[j].startswith("    "): blk.append(L[j]); j += 1
        exec("\n".join(blk), ns)
mesh, kc = ns["mesh"], ns["kc"]
ra = sys.argv[1]; os.makedirs(ra, exist_ok=True)
for d in sys.argv[2:]:
    j = os.path.basename(os.path.normpath(d)); T = json.load(open(f"{d}/thong_tin.json"))
    ft = f"{d}/mesh_vo_bao_truoc.npz" if os.path.exists(f"{d}/mesh_vo_bao_truoc.npz") else f"{d}/mesh_truoc_giu.npz"
    A = np.load(ft); B = np.load(f"{d}/mesh_vo_bao.npz")
    ks = sorted(set(k[:-3] for k in A.files if k.endswith("__v")) | set(k[:-3] for k in B.files if k.endswith("__v")))
    hd_bd = (T.get("ban_u5b") or {}).get("hd_truoc_sau") or {}; bp = {}; lech = 0.0
    for k in ks:
        m0, m1 = mesh(A, k), mesh(B, k)
        if m0 is None or m1 is None: bp[k] = {"hd": None, "mat": m1 is None}; continue
        if m0[0].shape == m1[0].shape and m0[1].shape == m1[1].shape and np.array_equal(m0[0], m1[0]) and np.array_equal(m0[1], m1[1]): h = [0.0, 0.0]
        else: h = [float(kc(m1[0], m0).max()), float(kc(m0[0], m1).max())]
        bp[k] = {"hd": h, "hd_max": max(h), "c": m1[0].mean(0).tolist(), "lo": m1[0].min(0).tolist(), "hi": m1[0].max(0).tolist()}
        if k in hd_bd and hd_bd[k]: lech = max(lech, abs(max(h) - max(hd_bd[k])))
    hs = [e["hd_max"] for k, e in bp.items() if e.get("hd") is not None and k not in (T.get("cha") or {}) and k != "mang"]
    out = {"viec": j, "nguon_truoc": ft, "md5_n2_tv": "14b46b125905736cb4217fd7d7b27a51", "hd_max_tru_mang_tang2": max(hs), "n": len(hs),
           "hd_max_tat_ca": max(e["hd_max"] for e in bp.values() if e.get("hd") is not None), "lech_bo_dung": lech, "bo_phan": bp}
    json.dump(out, open(f"{ra}/hd_{j}.json", "w")); print("HD", j, round(out["hd_max_tru_mang_tang2"], 4), "n", len(hs), "lech", round(lech, 5), flush=True)
