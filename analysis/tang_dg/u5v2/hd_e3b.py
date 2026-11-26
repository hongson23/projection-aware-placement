# -*- coding: utf-8 -*-
"""hd_e3b.py (21/11/2026) — Hausdorff hai chiều trước/sau chiếu cho việc e3b_B_* (v122e3 không có --hd_truoc_sau), tính từ mesh đã giữ
kq/<id>/mesh_vo_bao_truoc.npz ↔ mesh_vo_bao.npz đúng như --hd_truoc_sau của bộ dựng (igl.point_mesh_squared_distance, đỉnh → mặt, hai chiều;
bỏ màng và bộ phận tầng 2 như bộ dựng). Thêm hd_ts_* vào u5v2/doc_lap/<id>.json (cùng hàm thong_ke của phan_tich_doc_lap_u5.py),
đánh dấu "hd_ts_nguon". Dùng: python hd_e3b.py id1 id2 ..."""
import os, sys, json, time, collections
import numpy as np, igl
U = os.path.expanduser("~/Vision_team_data_189/sonhh20/ch3_tebao_dep/bo_dung_chung/tang_dg/u5v2"); O = os.path.dirname(U)
DL = os.environ.get("DOC_LAP", U + "/doc_lap")
sys.path.insert(0, U); from phan_tich_doc_lap_u5 import thong_ke
def kc(P, V, F): return np.sqrt(np.maximum(igl.point_mesh_squared_distance(np.asarray(P, float), np.asarray(V, float), np.asarray(F, np.int64))[0], 0))
for id_ in sys.argv[1:]:
    f = f"{DL}/{id_}.json"; d = f"{O}/kq/{id_}"
    if not os.path.exists(f) or not os.path.exists(f"{d}/mesh_vo_bao_truoc.npz"): print(id_, "thieu doc_lap/mesh truoc"); continue
    t0 = time.time(); z1, z0 = np.load(f"{d}/mesh_vo_bao.npz"), np.load(f"{d}/mesh_vo_bao_truoc.npz"); tt = json.load(open(f"{d}/thong_tin.json"))
    BO_QUA = {"mang"} | {c["ten"] for x in (tt.get("tang2") or []) for c in x.get("con", [])}
    ten = sorted({k[:-3] for k in list(z0.files) + list(z1.files) if k.endswith("__v")}); hts = {}
    for k in ten:
        if k + "__v" not in z1.files or k + "__v" not in z0.files: hts[k] = None; continue
        v1, f1, v0, f0 = z1[k + "__v"], z1[k + "__f"], z0[k + "__v"], z0[k + "__f"]
        if k in BO_QUA: continue
        hts[k] = [round(float(kc(v1, v0, f0).max()), 4), round(float(kc(v0, v1, f1).max()), 4)]
    lo = tt["loai"]; lo.setdefault("luoi_noi_chat", "luoi_noi_chat"); loai_cua = lambda k: lo.get(k, k.rstrip("0123456789").rstrip("_")); tu = tt["tu_the"]
    byt = collections.defaultdict(list)
    for k, v in hts.items():
        if v: byt[loai_cua(k)].append(max(v))
    r = json.load(open(f))
    r["hd_ts_theo_loai"] = {l: thong_ke(v) for l, v in byt.items()}; r["hd_ts_tat_ca"] = thong_ke([x for v in byt.values() for x in v])
    r["hd_ts_tu_the"] = thong_ke([max(v) for k, v in hts.items() if v and k in tu]); r["hd_ts_giay"] = round(time.time() - t0, 1)
    r["hd_ts_nguon"] = "hd_e3b.py tu mesh_vo_bao_truoc.npz/mesh_vo_bao.npz (float32)"; r["hd_ts_tung_bo_phan"] = hts
    json.dump(r, open(f + ".tam", "w"), ensure_ascii=False, indent=1); os.replace(f + ".tam", f); print("xong", id_, r["hd_ts_giay"], "s", flush=True)
