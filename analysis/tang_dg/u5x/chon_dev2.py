# -*- coding: utf-8 -*-
"""chon_dev2.py (23/11/2026) — LUẬT CHỌN GHI TRƯỚC cho lưới dev CX GỘP (14 cấu hình cũ u5_x_dev_* + 6 mới u5_x_dev2_*), seed dev 6597, 6598.
Thiết kế: Inside3D/ThuNghiem/doi_chung_cx_20261121/ThietKe_CX_MoRongDev_20261122.md bản 2 §3 (PI duyệt 23/11 07:45). Luật = chon_dev.py, chỉ thêm:
  - gộp lưới (nhóm theo nhãn; cấu hình eq KHÔNG xét vào best); HS = hằng số cũ (luoi_dev.json) + hang_so_bo_sung (luoi_dev2.json: lớn 48:10 := lớn 64:20);
  - phá hoà: Bh, B, A ⇒ mức L1, L0, L2 ⇒ λ1 nhỏ ⇒ dem nhỏ ⇒ tên cfg (tất định);
  - báo kèm: số seed k_chon trùng P4, G_HD từng seed so với P4 cùng seed; cổng dev "CX có hiệu lực" = G_HD(thắng) ≤ 0,90 · G_HD(P4) (≥ 10 %, như lượt chính);
  - EQ_MOI / BEST_MOI = cấu hình thắng thuộc lưới mới (⇒ tao_viec_cx2.py chinh2 chèn lượt chính); kiem_cu = chọn lại chỉ trên lưới cũ (phải ra eq BL1l1m800, best None).
Trần (không đổi): eq trung vị ti_dem ≤ 1,10 (vs U5b dev cùng seed); best trung vị ti_dem_p8 ≤ 1,10. An toàn (không nới): Σ mất ≤ Σ mất(U5b) + 1; Σ d50 ≤ Σ d50(U5b) + 2.
Nhóm không có cấu hình nào qua ⇒ None = BẤT KHẢ (ghi lý do từng cấu hình), không lấy cấu hình yếu nhất thay. Nhóm eq luôn còn BL1l1m800 (đã qua) nên không assert dừng.
Ra: u5x/chon_dev2.json + chon_dev2.md. CX_R (môi trường) chỉ để chạy thử trên bản sao."""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import dem_buoc as DB
R = os.environ.get("CX_R", os.path.expanduser("~/Vision_team_data_189/sonhh20/ch3_tebao_dep")); O = f"{R}/bo_dung_chung/tang_dg"; U = f"{O}/u5x"; H = 0.1796875
SEEDS = (6597, 6598); TRAN = {"best": 1.10, "eq": 1.10}; UU = {"Bh": 0, "B": 1, "A": 2}; MUCO = {"L1": 0, "L0": 1, "L2": 2}
EQ_CU = "BL1l1m800"
def thuoc(d):
    hd = d["ban_u5b"]["hd_truoc_sau"] or {}; v = [max(x) for k, x in hd.items() if not k.startswith("_") and x]
    tt = d["the_tich_truoc_sau"]; vb = sum(a for a, b in tt.values()); vm = sum(max(a - b, 0.0) for a, b in tt.values())
    cx = (d.get("bao_cao_bo_tri") or {}).get("cx") or {}
    return {"G_HD": float(sum(max(x - H, 0.0) for x in v)), "HDmax": float(max(v)) if v else None, "Vmat": 100 * vm / max(vb, 1e-12),
            "mat": len(d.get("mat_bo_phan") or []), "d50": len(d.get("con_lai_duoi_50pt") or []), "giay": d["giay_tong"], "k_chon": d["bao_cao_bo_tri"].get("cau_hinh_chon"),
            "n_su_kien": cx.get("n_su_kien"), "lam_cuoi": [(c_.get("lam_lon"), c_.get("lam_nho")) for c_ in cx.get("cau_hinh", [])]}
LU1 = json.load(open(f"{U}/luoi_dev.json")); LU2 = json.load(open(f"{U}/luoi_dev2.json"))
HS = json.loads(json.dumps(LU1["m_theo_muc"]["hang_so_buoc"]))
for k, v in LU2["hang_so_bo_sung"].items(): HS[k].update(v)
def doc(id_):
    f = f"{O}/kq/{id_}/thong_tin.json"
    if not (os.path.exists(f) and os.path.exists(f"{O}/xong/{id_}")): return None
    d = json.load(open(f)); r = thuoc(d); r["dem"] = DB.chi_phi(d, HS); r["sach_loi"] = DB.sach(id_); return r
S_P8 = {s: DB.so_buoc(json.load(open(f"{O}/kq/u5_v2P8_K8_{s}/thong_tin.json"))) for s in SEEDS}
ref = {s: doc(f"u5_x_dev_U5b_{s}") for s in SEEDS}; p4 = {s: doc(f"u5_x_dev_P4_{s}") for s in SEEDS}
assert all(ref.values()) and all(p4.values()), "thiếu U5b/P4 tham chiếu dev"
CFG = [(c, p, "u5_x_dev_", "cu") for c, p in LU1["cfg"].items()] + [(c, p, LU2["tien_to"], "moi") for c, p in LU2["cfg"].items()]
assert len({c for c, *_ in CFG}) == len(CFG) == 20, "lưới gộp phải có đúng 20 cấu hình, không trùng tên"
bang = []
for c, p, tt_, vong in CFG:
    r = {s: doc(f"{tt_}{c}_{s}") for s in SEEDS}; row = {"cfg": c, **p, "vong": vong, "du_seed": all(r.values()), "theo_seed": r}
    if row["du_seed"]:
        for k in ("G_HD", "HDmax", "Vmat"): row[k] = float(np.median([r[s][k] for s in SEEDS]))
        row["ti_gio"] = float(np.median([r[s]["giay"] / ref[s]["giay"] for s in SEEDS]))
        row["dem_tv"] = float(np.median([r[s]["dem"] for s in SEEDS]))
        row["ti_dem"] = float(np.median([r[s]["dem"] / ref[s]["dem"] for s in SEEDS])); row["ti_dem_p8"] = float(np.median([r[s]["dem"] / S_P8[s] for s in SEEDS]))
        row["mat"] = sum(r[s]["mat"] for s in SEEDS); row["d50"] = sum(r[s]["d50"] for s in SEEDS)
        row["an_toan"] = row["mat"] <= sum(ref[s]["mat"] for s in SEEDS) + 1 and row["d50"] <= sum(ref[s]["d50"] for s in SEEDS) + 2
        row["k_chon_trung_P4"] = sum(1 for s in SEEDS if r[s]["k_chon"] == p4[s]["k_chon"])
        row["G_HD_vs_P4_theo_seed"] = {s: [r[s]["G_HD"], p4[s]["G_HD"]] for s in SEEDS}
    bang.append(row)
def qua_gio(nhom, x): return (x["ti_dem_p8"] <= TRAN["best"]) if nhom == "best" else (x["ti_dem"] <= TRAN["eq"])
def chon(nhom, chi_cu=False):
    ds = [x for x in bang if nhom in x["nhom"] and x["du_seed"] and x["an_toan"] and qua_gio(nhom, x) and (x["vong"] == "cu" or not chi_cu)]
    if not ds: return None
    g0 = min(x["G_HD"] for x in ds); ds = [x for x in ds if x["G_HD"] <= max(g0 * 1.05, g0 + 0.005)]
    h0 = min(x["HDmax"] for x in ds); ds = [x for x in ds if x["HDmax"] <= h0 * 1.02]
    v0 = min(x["Vmat"] for x in ds); ds = [x for x in ds if x["Vmat"] <= v0 * 1.02 + 1e-9]
    x = sorted(ds, key=lambda x: (UU[x["neo"]], MUCO[x["muc"]], x["l0"], x["dem_tv"], x["cfg"]))[0]
    return {k: x[k] for k in ("cfg", "vong", "neo", "muc", "m", "l0", "G_HD", "HDmax", "Vmat", "ti_dem", "ti_dem_p8", "ti_gio", "k_chon_trung_P4")}
def ly_do(nhom):
    out = []
    for x in bang:
        if nhom not in x["nhom"]: continue
        if not x["du_seed"]: out.append({"cfg": x["cfg"], "ly_do": "thiếu seed/lỗi"}); continue
        l = ([] if qua_gio(nhom, x) else ["trần đếm"]) + ([] if x["an_toan"] else [f"an toàn (mất {x['mat']}, d50 {x['d50']})"])
        out.append({"cfg": x["cfg"], "vong": x["vong"], "G_HD": x["G_HD"], "ly_do": ", ".join(l) or "qua"})
    return out
best, eq = chon("best"), chon("eq"); kiem_cu = {"eq": (chon("eq", True) or {}).get("cfg"), "best": (chon("best", True) or {}).get("cfg")}
assert kiem_cu == {"eq": EQ_CU, "best": None}, f"chọn lại trên lưới cũ không ra như 22/11: {kiem_cu}"
assert eq is not None, "nhóm eq rỗng — không thể (BL1l1m800 đã qua)"
p4m = float(np.median([p4[s]["G_HD"] for s in SEEDS]))
hl = {g: (x is not None and x["G_HD"] <= 0.90 * p4m) for g, x in (("eq", eq), ("best", best))}
out = {"eq": eq, "best": best, "EQ_MOI": eq["vong"] == "moi", "BEST_MOI": bool(best and best["vong"] == "moi"), "BEST_BAT_KHA": best is None,
       "CX_HIEU_LUC_DEV": hl, "nguong_hieu_luc": "G_HD(thắng) ≤ 0,90 · G_HD(P4) trung vị", "P4_G_HD_tv": p4m, "kiem_cu": kiem_cu,
       "ly_do_best": ly_do("best"), "ly_do_eq": ly_do("eq"), "hang_so": HS, "tran": TRAN, "U5b_ref": ref, "P4_ref": p4, "bang": bang}
json.dump(out, open(f"{U}/chon_dev2.json", "w"), indent=1, ensure_ascii=False)
L = ["# Chọn cấu hình CX trên dev — lưới GỘP 20 cấu hình (luật ghi trước, chon_dev2.py) — CHƯA rà soát chéo", "",
     f"eq: **{eq['cfg']}** ({eq['vong']}); best: **{best['cfg'] + ' (' + best['vong'] + ')' if best else 'BẤT KHẢ'}**; EQ_MOI {out['EQ_MOI']}; BEST_MOI {out['BEST_MOI']}",
     f"CX có hiệu lực trên dev (≥ 10 % dưới P4 {p4m:.4f}): {hl}; kiểm lưới cũ: {kiem_cu}", "",
     "| cfg | vòng | nhóm | G_HD tv | HDmax tv | Vmat % tv | ti_dem (U5b) | ti_dem (P8) | ti_gio (mô tả) | mất | d50 | an toàn | k_chon = P4 | G_HD seed vs P4 |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for x in bang:
    if x["du_seed"]:
        gs = "; ".join(f"{a:.3f}/{b:.3f}" for a, b in x["G_HD_vs_P4_theo_seed"].values())
        L.append(f"| {x['cfg']} | {x['vong']} | {','.join(x['nhom'])} | {x['G_HD']:.4f} | {x['HDmax']:.4f} | {x['Vmat']:.3f} | {x['ti_dem']:.4f} | {x['ti_dem_p8']:.3f} | {x['ti_gio']:.3f} | {x['mat']} | {x['d50']} | {x['an_toan']} | {x['k_chon_trung_P4']}/2 | {gs} |")
    else: L.append(f"| {x['cfg']} | {x['vong']} | {','.join(x['nhom'])} | thiếu seed | | | | | | | | | | |")
L += ["", "Lý do nhóm best: " + json.dumps(out["ly_do_best"], ensure_ascii=False)]
open(f"{U}/chon_dev2.md", "w").write("\n".join(L) + "\n"); print("\n".join(L[:4]))
