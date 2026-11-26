# -*- coding: utf-8 -*-
"""kha_thi_n14d5.py (22/11/2026) — N14 mật độ cao: số học khả thi (CPU, một tiến trình, chỉ đọc đầu ra cũ + SDF đã có).
Đọc: thong_tin.json của 45 việc N14 cũ (D4-like), việc D5 / k08 / D4 tổng hợp đã xong, sdf.npz 5 tế bào, he_so_steiner_khuon.json (khám phá 18/11).
Định nghĩa φ_eff = Σ V(K_i ⊕ δ_i B) / V_kd như KhamPha_NguongMatDo_20261118 §2 (δ_ti = (kcl_ti + β)/2, δ_ly = (kcl_ly + β)/2,
chướng ngại Golgi/trung thể δ = S_roi + β − δ_ti), V_kd = ti_le_chiem.kha_dung_um3 của bộ dựng (cùng định nghĩa mọi ô).
Ghi: n14d5/kha_thi.json (+ in bảng)."""
import os, json, glob, math, statistics as st
import numpy as np
R = os.path.expanduser("~/Vision_team_data_189/sonhh20/ch3_tebao_dep"); B = R + "/bo_dung_chung"; KQ = B + "/tang_dg/kq"; OUT = B + "/n14d5"
K = json.load(open(OUT + "/he_so_steiner_khuon.json")); PI = math.pi
CELL = ["798836", "989830", "593409", "707681", "997092"]
def dil(V, S, M, d): return V + S * d + M * d * d + 4 / 3 * PI * d ** 3
def P_(dt): d = json.load(open(R + "/" + dt)); return {k: v["gia_tri"] for k, v in d["tham_so_khong_thuoc_luat"].items()}
def mot(job, S_roi=None):
    d = json.load(open(f"{KQ}/{job}/thong_tin.json")); P = P_(d["dac_ta"]); beta = d["ban_v12"]["beta"]; Sr = P["S_roi"] if S_roi is None else S_roi
    dl = {"ti_the": (P["kcl_ti_the"] + beta) / 2, "lysosome": (P["kcl_lysosome"] + beta) / 2}; dch = Sr + beta - dl["ti_the"]
    r = {"V_kd": d["ti_le_chiem"]["kha_dung_um3"], "n": {}, "Veff": {}, "V": {}}
    for k, t in d["tu_the"].items():
        L = d["loai"].get(k)
        if L not in ("ti_the", "lysosome", "golgi", "trung_the"): continue
        V0, S0, M0 = K[t["khuon"]]; s = t["s"]; dd = dl.get(L, dch)
        r["n"][L] = r["n"].get(L, 0) + 1; r["V"][L] = r["V"].get(L, 0) + s ** 3 * V0; r["Veff"][L] = r["Veff"].get(L, 0) + dil(s ** 3 * V0, s * s * S0, s * M0, dd)
    r["phi_nho"] = (r["Veff"].get("ti_the", 0) + r["Veff"].get("lysosome", 0)) / r["V_kd"]
    r["phi_cn"] = (r["Veff"].get("golgi", 0) + r["Veff"].get("trung_the", 0)) / r["V_kd"]; r["phi_eff"] = r["phi_nho"] + r["phi_cn"]
    r["v1_ti"] = r["Veff"].get("ti_the", 0) / max(r["n"].get("ti_the", 1), 1); r["v1_ly"] = r["Veff"].get("lysosome", 0) / max(r["n"].get("lysosome", 1), 1)
    r["cn_um3"] = r["Veff"].get("golgi", 0) + r["Veff"].get("trung_the", 0)
    return r
def co(m): return [os.path.basename(x) for x in sorted(glob.glob(f"{KQ}/{m}")) if os.path.exists(x + "/thong_tin.json")]
def tv(a): a = [x for x in a if x is not None]; return float(np.median(a)) if a else None
def main():
    kq = {"tong_hop": {}, "te_bao": {}}
    for ten, m in (("D4", "e3b_B_D4_td_62*"), ("D5", "e3b_B_D5_td_62*"), ("D5_U5b", "u5_c_U5b_D5_62*"), ("k08_P4", "u5_r_P4_k08_65*"), ("k08_U5b", "u5_r_U5b_k08_65*")):
        rr = [mot(j) for j in co(m)]
        kq["tong_hop"][ten] = {"n_viec": len(rr), **{k: tv([r[k] for r in rr]) for k in ("V_kd", "phi_nho", "phi_cn", "phi_eff", "v1_ti", "v1_ly", "cn_um3")}}
    D5 = kq["tong_hop"]["D5"]; K8 = kq["tong_hop"]["k08_P4"]
    PHI_D5 = tv([mot(j)["phi_eff"] for j in co("e3b_B_D5_td_62*") + co("u5_c_U5b_D5_62*")])
    PHI_K8 = tv([mot(j)["phi_eff"] for j in co("u5_r_P4_k08_65*") + co("u5_r_U5b_k08_65*")])
    kq["muc_tieu"] = {"phi_eff_D5": PHI_D5, "phi_eff_k08": PHI_K8, "nguon": "trung vị gộp PLACE + U5b tổng hợp (e3b_B_D5_td, u5_c_U5b_D5 | u5_r_P4_k08, u5_r_U5b_k08)"}
    DV_syn = D5["V_kd"] - K8["V_kd"]
    kn = K["nhan"]; z = np.linspace(-8, 8, 160001); w = np.exp(-z ** 2 / 2); w /= w.sum(); s_tb = float(((3.9 + 1.1 / (1 + np.exp(-z))) ** 3 * w).sum() ** (1 / 3))
    Vn = [kn[0] * s_tb ** 3, kn[1] * s_tb ** 2, kn[2] * s_tb]
    DVc_syn = dil(*Vn, 2.0 + 0.8) - dil(*Vn, 2.0 + 0.3)
    c_cal = DV_syn / DVc_syn
    kq["hieu_chinh_k08"] = {"DV_kd_tong_hop_do": DV_syn, "DV_vo_steiner_tong_hop": DVc_syn, "he_so": c_cal, "s_tb": s_tb}
    zs = np.load(R + "/bo_dung_chung/v121/sdf_mang/sdf_75005fbc7d88018f.npz")["sdf"]; Ns = 128; ss = zs[2:2 * Ns + 2:2, 2:2 * Ns + 2:2, 2:2 * Ns + 2:2]
    Vsyn = {r_: float((ss < -(0.5 + r_)).sum()) * (23 / 128) ** 3 - dil(*Vn, 2.3 + r_) for r_ in (0.0, 0.25, 0.5, 1.0)}
    kq["co_vung_tong_hop_D4"] = {str(r_): Vsyn[r_] / Vsyn[0.0] for r_ in (0.25, 0.5, 1.0)}; kq["V_chon_tong_hop_D4_steiner"] = Vsyn[0.0]
    chon = {c["CellId"]: c for c in json.load(open(B + "/n14/chon.json"))["chon"]}
    H1 = 23 / 128
    for cid in CELL:
        rr = [mot(j) for j in co(f"n14_{cid}_*_650[123]")]; cn8 = tv([mot(j, S_roi=0.8)["cn_um3"] for j in co(f"n14_{cid}_*_650[123]")])
        z = np.load(f"{B}/n14/te_bao/{cid}/sdf.npz"); NV = z["NV"]; sm = z["sdf"][2:2 * NV[0] + 2:2, 2:2 * NV[1] + 2:2, 2:2 * NV[2] + 2:2]; sn = z["sdf_nhan"][2:2 * NV[0] + 2:2, 2:2 * NV[1] + 2:2, 2:2 * NV[2] + 2:2]
        tr = sm < -0.5; Vc = {o: float((tr & (sn > o)).sum()) * H1 ** 3 for o in (2.3, 2.8)}
        Vtb = float((sm < 0).sum()) * H1 ** 3; Vnh = float((sn < 0).sum()) * H1 ** 3
        x = {k: tv([r[k] for r in rr]) for k in ("V_kd", "phi_nho", "phi_cn", "phi_eff", "v1_ti", "v1_ly", "cn_um3")}
        x["V_kd_min_max"] = [min(r["V_kd"] for r in rr), max(r["V_kd"] for r in rr)]; x["n_viec"] = len(rr); x["n_cu"] = rr[0]["n"]
        x["V_te_bao_luoi"] = Vtb; x["V_nhan_luoi"] = Vnh; x["V_chon_2.3"] = Vc[2.3]; x["V_chon_2.8"] = Vc[2.8]
        x["rho_chon"] = chon[cid]["rho"]; x["rho_kd"] = x["V_kd"] / D5["V_kd"]
        nA = (max(1, round(600 * x["rho_chon"])), max(1, round(200 * x["rho_chon"])))
        nB = (max(1, round(600 * x["rho_kd"])), max(1, round(200 * x["rho_kd"])))
        def phi_n(n, Vkd, cn): return (n[0] * x["v1_ti"] + n[1] * x["v1_ly"]) / Vkd, cn / Vkd
        for ten, n in (("A_600rho_chon", nA), ("B_600rho_kd", nB)):
            a, b = phi_n(n, x["V_kd"], x["cn_um3"]); x[ten] = {"ti_the": n[0], "lysosome": n[1], "phi_nho": a, "phi_cn": b, "phi_eff": a + b}
            Vk8 = x["V_kd"] - c_cal * (Vc[2.3] - Vc[2.8]); a8, b8 = phi_n(n, Vk8, x["cn_um3"] * 1.0)
            x[ten]["k08"] = {"V_kd_uoc": Vk8, "phi_nho": a8, "phi_cn_duoi": b8, "phi_eff": a8 + b8}
        def luat_C(phi_t, Vkd, cn):
            k = max(phi_t - cn / Vkd, 0) * Vkd / (3 * x["v1_ti"] + x["v1_ly"]); n = (max(1, round(3 * k)), max(1, round(k)))
            a, b = (n[0] * x["v1_ti"] + n[1] * x["v1_ly"]) / Vkd, cn / Vkd; return {"ti_the": n[0], "lysosome": n[1], "phi_nho": a, "phi_cn": b, "phi_eff": a + b, "V_kd": Vkd}
        Vk8 = x["V_kd"] - c_cal * (Vc[2.3] - Vc[2.8]); x["cn_um3_k08"] = cn8
        x["C_D5"] = luat_C(PHI_D5, x["V_kd"], x["cn_um3"]); x["C_k08"] = luat_C(PHI_K8, Vk8, cn8)
        nC = (x["C_D5"]["ti_the"], x["C_D5"]["lysosome"]); x["C_D5_dem_o_k08"] = {"phi_eff": (nC[0] * x["v1_ti"] + nC[1] * x["v1_ly"] + cn8) / Vk8}
        x["co_vung"] = {str(r_): float((sm < -(0.5 + r_)).__and__(sn > 2.3 + r_).sum()) * H1 ** 3 / Vc[2.3] for r_ in (0.25, 0.5, 1.0)}
        kq["te_bao"][cid] = x
    json.dump(kq, open(OUT + "/kha_thi.json", "w"), indent=1, ensure_ascii=False)
    print(json.dumps(kq["tong_hop"], indent=0)); print(kq["hieu_chinh_k08"])
    print("muc tieu", kq["muc_tieu"], "co vung tong hop", kq["co_vung_tong_hop_D4"], "Vchon_syn", round(Vsyn[0.0], 1))
    for cid, x in kq["te_bao"].items():
        print(cid, "C_D5", x["C_D5"]["ti_the"], x["C_D5"]["lysosome"], round(x["C_D5"]["phi_eff"], 3), "(cn", round(x["C_D5"]["phi_cn"], 3), ") C_k08", x["C_k08"]["ti_the"], x["C_k08"]["lysosome"], round(x["C_k08"]["phi_eff"], 3), "(cn", round(x["C_k08"]["phi_cn"], 3), "Vkd8", round(x["C_k08"]["V_kd"], 1), ") demD5@k08", round(x["C_D5_dem_o_k08"]["phi_eff"], 3), "co_vung", {k: round(v, 3) for k, v in x["co_vung"].items()}, "v1", round(x["v1_ti"], 3), round(x["v1_ly"], 3))
    for cid, x in kq["te_bao"].items():
        print(cid, "Vkd", round(x["V_kd"], 1), x["V_kd_min_max"], "rho_chon", round(x["rho_chon"], 4), "rho_kd", round(x["rho_kd"], 4), "phi_cu", round(x["phi_eff"], 3), "(nho", round(x["phi_nho"], 3), "cn", round(x["phi_cn"], 3), ")",
              "| A", x["A_600rho_chon"]["ti_the"], x["A_600rho_chon"]["lysosome"], round(x["A_600rho_chon"]["phi_eff"], 3), "k08", round(x["A_600rho_chon"]["k08"]["phi_eff"], 3),
              "| B", x["B_600rho_kd"]["ti_the"], x["B_600rho_kd"]["lysosome"], round(x["B_600rho_kd"]["phi_eff"], 3), "k08", round(x["B_600rho_kd"]["k08"]["phi_eff"], 3), round(x["B_600rho_kd"]["k08"]["V_kd_uoc"], 1),
              "Vtb", round(x["V_te_bao_luoi"]), "Vnh", round(x["V_nhan_luoi"]), "V2.3", round(x["V_chon_2.3"], 1), "V2.8", round(x["V_chon_2.8"], 1))

if __name__ == "__main__":
    main()
