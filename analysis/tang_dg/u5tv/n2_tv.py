# -*- coding: utf-8 -*-
"""n2_tv.py (22/11/2026, M7 của rà soát trước chạy) — thước N2 cho tế bào THỰC VẬT, cùng lưới, từ mesh TRƯỚC chiếu (mesh_truoc_giu.npz) và
mesh SAU chiếu (mesh_vo_bao.npz) của chính việc; không đọc nhánh. Mỗi bộ phận rắn trong tu_the:
  HD hai chiều (đỉnh → mặt, igl.point_mesh_squared_distance, như --hd_truoc_sau của v122u5b); hai mesh trùng hệt ⇒ 0;
  thể tích mesh trước/sau (công thức phân kỳ, |V|); mất ⇔ không có mesh sau; < 50 % ⇔ V_sau < 0,5 V_trước.
Theo loại: Σ V trước, Σ V sau (tuyệt đối), % mất = Σ max(V_trước − V_sau, 0) / Σ V_trước. Toàn tế bào: HD max, HD p95 theo bộ phận (bộ phận
không đổi = 0). ER, không bào: lấy từ thong_tin.json ban_u5b.hd_vung_truoc (dựng lại từ mesh trước chiếu trong bộ dựng) nếu có.
Kiểm chéo: chênh lệch tối đa với HD của bộ dựng (ban_u5b.hd_truoc_sau) ⇒ ghi 'lech_bo_dung'.
Dùng: python n2_tv.py <thu_muc_viec>  ⇒ <thu_muc_viec>/n2.json"""
import sys, os, json, numpy as np, igl
d = sys.argv[1]; T = json.load(open(f"{d}/thong_tin.json"))
A = np.load(f"{d}/mesh_truoc_giu.npz"); B = np.load(f"{d}/mesh_vo_bao.npz")
def mesh(Z, k):
    return (np.asarray(Z[f"{k}__v"], np.float64), np.asarray(Z[f"{k}__f"], np.int64)) if f"{k}__v" in Z.files else None
def vol(m):
    V, F = m; a, b, c = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]; return abs(float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum()) / 6.0)
def kc(P, m): return np.sqrt(np.maximum(igl.point_mesh_squared_distance(P, m[0], m[1])[0], 0))
lo = T["loai"]; bp = {}; hd_bd = (T.get("ban_u5b") or {}).get("hd_truoc_sau") or {}; lech = 0.0
for k in T["tu_the"]:
    m0, m1 = mesh(A, k), mesh(B, k)
    e = {"loai": lo.get(k, k), "V0": vol(m0) if m0 else None, "V1": vol(m1) if m1 else 0.0, "mat": m1 is None}
    if m0 is None or m1 is None: e["hd"] = None
    elif m0[0].shape == m1[0].shape and m0[1].shape == m1[1].shape and np.array_equal(m0[0], m1[0]) and np.array_equal(m0[1], m1[1]): e["hd"] = [0.0, 0.0]
    else: e["hd"] = [float(kc(m1[0], m0).max()), float(kc(m0[0], m1).max())]
    if e["hd"] is not None and k in hd_bd and hd_bd[k] is not None: lech = max(lech, abs(max(e["hd"]) - max(hd_bd[k])))
    e["duoi50"] = bool(e["V0"] and e["V1"] < 0.5 * e["V0"]); bp[k] = e
hds = np.array([max(e["hd"]) if e["hd"] is not None else np.nan for e in bp.values()])
theo = {}
for e in bp.values():
    t = theo.setdefault(e["loai"], {"n": 0, "V0": 0.0, "V1": 0.0, "mat_V": 0.0, "mat": 0, "duoi50": 0, "hd_max": 0.0})
    t["n"] += 1; t["V0"] += e["V0"] or 0.0; t["V1"] += e["V1"]; t["mat_V"] += max((e["V0"] or 0.0) - e["V1"], 0.0); t["mat"] += e["mat"]; t["duoi50"] += e["duoi50"]
    t["hd_max"] = max(t["hd_max"], max(e["hd"]) if e["hd"] is not None else np.inf)
for t in theo.values(): t["pt_mat"] = 100 * t["mat_V"] / t["V0"] if t["V0"] else None
V0 = sum(t["V0"] for t in theo.values()); ML = sum(t["mat_V"] for t in theo.values())
vung = {}
for k, e in ((T.get("ban_u5b") or {}).get("hd_vung_truoc") or {}).items():
    if k.startswith("_") or not isinstance(e, dict): continue
    vung[k] = {"hd": e.get("hd"), "V0": e.get("the_tich_truoc"), "V1": e.get("the_tich_sau"),
               "pt_mat": 100 * max(e["the_tich_truoc"] - e["the_tich_sau"], 0) / e["the_tich_truoc"] if e.get("the_tich_truoc") else None}
hd_vung = [max(v["hd"]) for v in vung.values() if v.get("hd")]
out = {"hd_max_ran": float(np.nanmax(hds)) if len(hds) else None, "hd_p95_ran": float(np.nanpercentile(hds, 95)) if len(hds) else None,
       "hd_max_tat_ca": max([float(np.nanmax(hds))] + hd_vung) if len(hds) else None, "n_nan": int(np.isnan(hds).sum()),
       "pt_mat_ran": 100 * ML / V0 if V0 else None, "mat": sum(t["mat"] for t in theo.values()), "duoi50": sum(t["duoi50"] for t in theo.values()),
       "theo_loai": theo, "vung": vung, "lech_bo_dung": lech, "bo_phan": bp}
json.dump(out, open(f"{d}/n2.json.moi", "w"), default=float); os.replace(f"{d}/n2.json.moi", f"{d}/n2.json")
print("N2", d, {k: out[k] for k in ("hd_max_ran", "hd_p95_ran", "hd_max_tat_ca", "pt_mat_ran", "mat", "duoi50", "lech_bo_dung")})
