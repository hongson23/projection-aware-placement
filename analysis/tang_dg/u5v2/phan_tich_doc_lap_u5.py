# -*- coding: utf-8 -*-
"""phan_tich_doc_lap_u5.py (20/11/2026) — kiểm độc lập U5 bậc 0 từ dữ liệu GỐC (không đọc U5.json).
Nguồn: kq/<id>/thong_tin.json (bộ dựng), bo_kiem_v3/cham_ra/vstdgu5_mu_sau__<mu>/kiem3.csv + do_khuon.json.gz (bộ kiểm v3),
khoá mù cham_u5/khoa/KHOA_MU.json, mesh sau chiếu cham_u5/vs/vstdgu5_mu_sau/<mu>/mesh_vo_bao.npz.
Ra: u5v2/doc_lap/<id>.json (mỗi việc) + tong.json.  Dùng: python phan_tich_doc_lap_u5.py [--hd] [--j N]
--hd: thêm Hausdorff HAI CHIỀU từng bộ phận tu_the (khuôn tại tư thế lưu ↔ mesh sau chiếu), CPU.
"""
import os, sys, json, csv, gzip, collections, argparse, time
import numpy as np
R = os.path.expanduser("~/Vision_team_data_189/sonhh20/ch3_tebao_dep"); O = R + "/bo_dung_chung/tang_dg"
KV3 = R + "/bo_kiem_v3"; sys.path.insert(0, KV3)
H = 0.1797; MEM = {"4.5b", "4.6b", "4.8b"}
LUAT_MAT = {"ti_the": ["4.3a", "4.3b", "M.a", "M.c₀", "M.h"], "lysosome": ["4.6a", "M.a", "M.c₀", "M.h"]}; LUAT_MAT_MD = ["M.a", "M.c₀", "M.h"]
OUT = O + "/u5v2/doc_lap"

def thong_ke(xs):
    xs = [x for x in xs if x is not None]
    if not xs: return None
    a = np.asarray(xs, float)
    return {"n": len(a), "max": float(a.max()), "p95": float(np.percentile(a, 95)), "p50": float(np.median(a)), "tong": float(a.sum()),
            "n_tren_h": int((a > H).sum()), "n_tren_h2": int((a > H / 2).sum()), "vuot_h": float(np.maximum(a - H, 0).sum()),
            "vuot_h2": float(np.maximum(a - H / 2, 0).sum())}

def _fib(n):
    k = np.arange(n) + 0.5; ph = np.arccos(1 - 2 * k / n); th = np.pi * (1 + 5 ** 0.5) * k
    return np.stack([np.cos(th) * np.sin(ph), np.sin(th) * np.sin(ph), np.cos(ph)], 1)
def mat_khuon(ten, sdf_t, n=6000):
    """điểm trên MẶT khuôn đơn vị (chính xác cho elip/cầu; ống: vòng quanh trục + chỏm, lọc |φ| ≤ 1e-4)."""
    import do_khuon as DK, luat_v12 as LV
    u = _fib(n)
    if ten == "cau": return u
    if ten in ("nhan", "golgi"):
        a = np.asarray(LV.P("dv.khuon_nhan_truc" if ten == "nhan" else "dv.khuon_golgi_truc"), float)
        return u / np.sqrt(((u / a) ** 2).sum(1, keepdims=True))
    if ten.startswith("ong"):
        r, c = LV.P("dv.khuon_ong")[int(ten[3:])]; tau = np.linspace(-1, 1, 241)
        Q = np.stack([tau * (1 - r), c * (1 - tau ** 2), 0 * tau], 1); T_ = np.stack([np.full_like(tau, 1 - r), -2 * c * tau, 0 * tau], 1)
        T_ /= np.linalg.norm(T_, axis=1, keepdims=True); N1 = np.stack([-T_[:, 1], T_[:, 0], 0 * tau], 1); N2 = np.cross(T_, N1)
        an = np.linspace(0, 2 * np.pi, 48, endpoint=False)
        P = (Q[:, None] + r * (np.cos(an)[None, :, None] * N1[:, None] + np.sin(an)[None, :, None] * N2[:, None])).reshape(-1, 3)
        P = np.concatenate([P, Q[0] + r * u, Q[-1] + r * u])
        return P[np.abs(sdf_t(P)) <= 1e-4]
    return None
def hd_hai_chieu(V, F, tt):
    import do_khuon as DK, kiem2 as K
    kh = DK.dinh_nghia_khuon(tt["khuon"])
    if kh is None: return None
    sdf_t, nua, _ = kh; Rm = np.asarray(tt["R"], float); t = np.asarray(tt["t"], float); s = float(tt["s"])
    d_ab = np.abs(s * sdf_t(((V - t) @ Rm) / s)); sau_truoc = float(d_ab.max())
    Pk = mat_khuon(tt["khuon"], sdf_t)
    if Pk is None or not len(Pk): return {"sau_den_truoc": sau_truoc, "truoc_den_sau": None, "hd": sau_truoc}
    Pw = t + s * Pk @ Rm.T; L = K.Luoi(V, F); truoc_sau = float(np.sqrt(L.kc2_diem(Pw)).max())
    return {"sau_den_truoc": sau_truoc, "truoc_den_sau": truoc_sau, "hd": max(sau_truoc, truoc_sau), "n_diem_khuon": int(len(Pk))}

LO = "u5"
def mot(id_, mu, hd):
    tt = json.load(open(f"{O}/kq/{id_}/thong_tin.json")); C = f"{KV3}/cham_ra/vstdg{LO}_mu_sau__{mu}"
    rs = list(csv.DictReader(open(C + "/kiem3.csv", encoding="utf-8")))
    tu = tt["tu_the"]; lo = tt["loai"]; lo.setdefault("luoi_noi_chat", "luoi_noi_chat")
    loai_cua = lambda k: lo.get(k, k.rstrip("0123456789").rstrip("_"))
    st = [r for r in rs if r["phan_quyet"] == "su_that" and r["loai_luat"] == "su_that_v1.2"]
    truot = {r["luat"]: {"trang_thai": r["trang_thai"], "so_ca": int(r["so_ca"]), "so_ca_truot": int(r["so_ca_truot"]), "bo_phan": r["bo_phan_truot"]} for r in st if int(r["so_ca_truot"]) > 0}
    Di = {r["phan_quyet"]: r for r in rs if r["luat"] == "D_i"}
    mat = [k for k, r in Di.items() if r["trang_thai"] == "KHONG_CO_MESH"]
    pen = sum(len(LUAT_MAT.get(loai_cua(k), LUAT_MAT_MD)) for k in mat)
    Vt = sum(v["so_ca_truot"] for k, v in truot.items() if k not in MEM) + pen
    dl = set(tt.get("duoi_luoi") or []); b43 = [x.strip('"') for x in (truot.get("4.3b", {}).get("bo_phan") or "").split(";") if x.strip()]
    n43_dl = sum(1 for x in b43 if x in dl)
    gc = float([r for r in rs if r["luat"] == "G_c"][0]["trang_thai"])
    try:
        sys.path.insert(0, O); import phan_tich_tdg as PTD; gk = PTD.gc_K_kiem(rs, tu)
    except Exception as e: gk = {"loi": repr(e)}
    D_k = {k: float(r["trang_thai"]) for k, r in Di.items() if r["trang_thai"] not in ("KHONG_CO_MESH", "KHONG_DUNG_LAI_DUOC", "")}
    by = collections.defaultdict(lambda: collections.defaultdict(list))
    for k, v in (tt.get("doi_mat_um") or {}).items():
        if v and v[1] is not None: by[loai_cua(k)]["doi_max"].append(v[1]); by[loai_cua(k)]["doi_tb"].append(v[0])
    for k, d in D_k.items(): by[loai_cua(k)]["D_kiem"].append(d)
    vol = tt.get("the_tich_truoc_sau") or {}; V_l = collections.defaultdict(lambda: [0.0, 0.0, 0.0])
    for k, (a, b) in vol.items(): x = V_l[loai_cua(k)]; x[0] += a; x[1] += b; x[2] += max(a - b, 0)
    iou = collections.defaultdict(list)
    for k, v in (tt.get("iou_truoc_sau") or {}).items(): iou[loai_cua(k)].append(v)
    theo_loai = {l: {"doi_max": thong_ke(by[l]["doi_max"]), "doi_tb_tb": float(np.mean(by[l]["doi_tb"])) if by[l]["doi_tb"] else None,
                     "D_kiem": thong_ke(by[l]["D_kiem"]), "V_truoc": V_l[l][0], "V_sau": V_l[l][1], "V_mat": V_l[l][2],
                     "V_mat_ti_le": V_l[l][2] / max(V_l[l][0], 1e-12), "iou_min": min(iou[l]) if iou[l] else None,
                     "iou_tb": float(np.mean(iou[l])) if iou[l] else None} for l in set(by) | set(V_l)}
    tatca_doi = [v[1] for k, v in (tt.get("doi_mat_um") or {}).items() if v and v[1] is not None]
    tu_doi = [v[1] for k, v in (tt.get("doi_mat_um") or {}).items() if v and v[1] is not None and k in tu]
    Vt_all = sum(a for a, b in vol.values()); Vm_all = sum(max(a - b, 0) for a, b in vol.values())
    Vt_tu = sum(a for k, (a, b) in vol.items() if k in tu); Vm_tu = sum(max(a - b, 0) for k, (a, b) in vol.items() if k in tu)
    b = tt["bao_cao_bo_tri"]
    out = {"id": id_, "mu": mu, "Gc_K": gk.get("Gc_K"), "V_test_cong": gk.get("V_test"), "Gc_K_csv": gc, "Gc_K_tinh_lai": float(sum(max(d - H, 0) for d in D_k.values())), "n_mat_kiem": len(mat), "mat_kiem": mat,
           "V_test": Vt, "V_test_tru_43b_duoi_luoi": Vt - n43_dl, "n_43b_duoi_luoi": n43_dl, "luat_truot": truot, "theo_loai": theo_loai,
           "doi_tat_ca": thong_ke(tatca_doi), "doi_tu_the": thong_ke(tu_doi),
           "V_mat_ti_le_tat_ca": Vm_all / max(Vt_all, 1e-12), "V_mat_ti_le_tu_the": Vm_tu / max(Vt_tu, 1e-12), "V_mat_um3_tat_ca": Vm_all,
           "mat_bo_phan_bd": tt.get("mat_bo_phan"), "con_lai_duoi_50pt": tt.get("con_lai_duoi_50pt"), "vo_manh": tt.get("vo_manh"), "vo_do_chieu": tt.get("vo_do_chieu"),
           "n_duoi_luoi": len(tt.get("duoi_luoi") or []), "con_vp_max_um": tt.get("con_vp_max_um"), "diem_khong_kha_thi": tt.get("diem_khong_kha_thi_truoc_chieu"),
           "diem_vi_pham_ban_dau": tt.get("diem_vi_pham_ban_dau"), "P": b.get("phat_cung_day_du"), "k_chon": b.get("cau_hinh_chon"), "V_nhan_truoc": (vol.get("nhan") or [None])[0],
           "giay_dat": b.get("giay_dat"), "giay_chieu": tt.get("giay_chieu"), "giay_tong": tt.get("giay_tong"),
           "phat_truoc_chieu_tong_nho": float(sum(v for k, v in (tt.get("phat_cung_truoc_chieu") or {}).items() if loai_cua(k) in ("ti_the", "lysosome"))),
           "J_truoc_chieu": b.get("J_truoc_chieu"), "u5_cau_hinh": [{kk: c.get(kk) for kk in ("k", "P", "G_du_bao", "n_tren_h", "V_vanh", "phu_thuoc_ngoai_Gc")} for c in (b.get("u5") or {}).get("cau_hinh", [])]}
    hts = (tt.get("ban_u5b") or {}).get("hd_truoc_sau")
    if hts:
        byt = collections.defaultdict(list)
        for k, v in hts.items():
            if k.startswith("_") or not v: continue
            byt[loai_cua(k)].append(max(v))
        out["hd_ts_theo_loai"] = {l: thong_ke(v) for l, v in byt.items()}; out["hd_ts_tat_ca"] = thong_ke([x for v in byt.values() for x in v])
        out["hd_ts_tu_the"] = thong_ke([max(v) for k, v in hts.items() if not k.startswith("_") and v and k in tu]); out["hd_ts_giay"] = hts.get("_giay")
    if hd:
        z = np.load(f"{O}/cham_{LO}/vs/vstdg{LO}_mu_sau/{mu}/mesh_vo_bao.npz"); HD = {}
        for k, t_ in tu.items():
            if k + "__v" not in z.files: HD[k] = None; continue
            HD[k] = hd_hai_chieu(z[k + "__v"].astype(np.float64), z[k + "__f"].astype(np.int64), t_)
        out["hd"] = HD; byh = collections.defaultdict(lambda: collections.defaultdict(list))
        for k, v in HD.items():
            if v: byh[loai_cua(k)]["hd"].append(v["hd"]); byh[loai_cua(k)]["sau_den_truoc"].append(v["sau_den_truoc"]); byh[loai_cua(k)]["truoc_den_sau"].append(v["truoc_den_sau"])
        out["hd_theo_loai"] = {l: {m: thong_ke(byh[l][m]) for m in byh[l]} for l in byh}
        out["hd_tat_ca"] = thong_ke([v["hd"] for v in HD.values() if v])
    return out

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--hd", action="store_true"); ap.add_argument("--j", type=int, default=1); ap.add_argument("--chi", default=None)
    ap.add_argument("--lo", default="u5")
    a = ap.parse_args(); os.makedirs(OUT, exist_ok=True); LO = a.lo
    khoa = json.load(open(f"{O}/cham_{LO}/khoa/KHOA_MU.json"))["id"]
    viec = sorted((v["nguon"], mu) for mu, v in khoa.items() if v["ban"] == "sau")
    if a.chi: viec = [x for x in viec if x[0] == a.chi]
    def chay(x):
        t0 = time.time(); r = mot(x[0], x[1], a.hd); r["giay_phan_tich"] = round(time.time() - t0, 1)
        json.dump(r, open(f"{OUT}/{x[0]}{'_hd' if a.hd else ''}.json", "w"), ensure_ascii=False, indent=1); return x[0]
    if a.j > 1:
        from multiprocessing import Pool
        with Pool(a.j) as p:
            for i in p.imap_unordered(chay, viec): print("xong", i, flush=True)
    else:
        for x in viec: print("xong", chay(x), flush=True)
