# -*- coding: utf-8 -*-
"""phan_tich_tdg.py BẢN 2 (19/11/2026; thêm B2, SA đã chỉnh trong E1, E3a 260 tế bào, E3b D4/D4A2/D5, S-P) — phân tích GHI TRƯỚC của lượt 'tăng đóng góp' (E1, E2, E3a, E3b), Bài 1 ch3.
Đăng ký trước: Inside3D/ThuNghiem/DangKyTruoc_ThucNghiemTangDongGop_20261119.md (§8). md5 tệp này ghi ở §11 TRƯỚC việc đầu tiên.
Dùng:  python phan_tich_tdg.py e12   → $O/ket_qua/E12.json + E12.md  (E1 đường cong ngân sách, E2 cắt bỏ PLACE theo mật độ)
       python phan_tich_tdg.py e3    → $O/ket_qua/E3.json + E3.md    (E3a phát lại PROJECT, E3b dựng lại; cần chấm xong)
Biến DATA_TDG = thư mục gốc dữ liệu (mặc định bo_dung_chung/tang_dg) — để chạy thử trên dữ liệu giả.
Định nghĩa:
  P   = phạt cứng đầy đủ trước chiếu của cấu hình được chọn (bao_cao_bo_tri.phat_cung_day_du; với E3: ban_e3.P_goc / P_dung);
  y   = P / N, N = tổng số bộ phận (so_luong);  τ = 1,6e-3, τ_lo = 4e-4 µm²/bộ phận (đăng ký quét mật độ §7).
  Ngân sách b (đơn vị = ngân sách 1× theo SỐ ĐẾM của lượt chính): PLACE b = k/4 (tốt nhất trong k cấu hình đầu của việc K = 16);
  thử–loại b = k/K1 (K1 = 6 ở D5, 16 ở D4A2; tốt nhất trong k cấu hình đầu của việc K = 4·K1); SA trục vòng b = m (việc un_m, K = 4,
  tốt nhất trong 4); SA trục khởi lại (phụ) b = m·k/4.  Phạt từng cấu hình đọc từ dòng 'k <i> phat <L> <giây>' của log (GO_LOI=1).
  G_c^K = như phan_tich_chinh_dv.py (nhập hàm, không chép): Σ max(D_i − h, 0) của bộ kiểm + bán kính bao của bộ phận mất mesh;
  G_c^B = Σ max(doi_mat_um[i][1] − h, 0) (bộ dựng tự đo); h = 0,1797 µm.
"""
import os, sys, json, csv, re, math, collections, hashlib
import numpy as np
GOC = os.path.expanduser("~/Vision_team_data_189/sonhh20/ch3_tebao_dep")
O = os.environ.get("DATA_TDG", os.path.join(GOC, "bo_dung_chung", "tang_dg"))
KQ = os.path.join(O, "kq"); LOG = os.path.join(O, "log"); RA = os.path.join(O, "ket_qua")
KQ_CHINH = os.path.join(GOC, "bo_dung_chung", "v121", "kq_chinh_dv")
BANG_MO = os.path.join(GOC, "bo_kiem_v3", "cham_chinh_dv", "phan_tich", "bang_mo.csv")
sys.path.insert(0, os.path.join(GOC, "bo_kiem_v3", "cham_chinh_dv", "phan_tich"))
import phan_tich_chinh_dv as PT
H = PT.H; TAU = 1.6e-3; TAU_LO = 4e-4; SAN = 1e-9
def md5(p): return hashlib.md5(open(p, "rb").read()).hexdigest()
def tv(x):
    x = [v for v in x if v is not None and not (isinstance(v, float) and math.isnan(v))]
    return None if not x else {"tv": float(np.median(x)), "min": float(min(x)), "max": float(max(x)), "n": len(x)}
def doc_k(id_):
    f = os.path.join(LOG, id_ + ".log"); out = {}
    if not os.path.exists(f): return []
    for l in open(f, errors="replace"):
        m = re.match(r"^k (\d+) phat (\S+) (\d+)", l)
        if m: out[int(m.group(1))] = float(m.group(2))
    return [out[k] for k in sorted(out)] if out and sorted(out) == list(range(len(out))) else []
def doc_k_goc(id_):
    f = os.path.join(LOG, id_ + ".log"); out = {}
    if not os.path.exists(f): return []
    for l in open(f, errors="replace"):
        m = re.match(r"^k (\d+) phat_goc (\S+)", l)
        if m: out[int(m.group(1))] = float(m.group(2))
    return [out[k] for k in sorted(out)] if out and sorted(out) == list(range(len(out))) else []
def tot_dau(L, k): return min(L[:k]) if len(L) >= k else None
def tot_dau_goc(L, Lg, k):
    if len(L) < k or len(Lg) < k: return None
    i = min(range(k), key=lambda j: (L[j], j)); return Lg[i]

def e12():
    viec = [l.split()[0] for l in open(os.path.join(O, "viec.txt")) if l.strip()]
    J = {}
    for id_ in viec:
        if not id_.startswith(("e1_", "e2_", "b2_")): continue
        f = os.path.join(KQ, id_, "thong_tin_dat.json")
        if not os.path.exists(f): J[id_] = None; continue
        d = json.load(open(f)); b = d["bao_cao_bo_tri"]; N = sum(d["so_luong"].values()); Lk = doc_k(id_)
        J[id_] = {"P": b["phat_cung_day_du"], "N": N, "Lk": Lk, "Lg": doc_k_goc(id_), "khong_dat": b.get("khong_dat"), "giay_dat": b.get("giay_dat"),
                  "S_mem": (b.get("J_truoc_chieu") or {}).get("mem_lon", 0) + (b.get("J_truoc_chieu") or {}).get("mem_nho", 0),
                  "danh_gia_phat": b.get("danh_gia_phat"), "K": b["ngan_sach_nhanh"]["K"], "seed": d["seed"]}
        if Lk and not id_.startswith("b2_"): assert abs(min(Lk) - b["phat_cung_day_du"]) <= 1e-9 * max(1.0, abs(b["phat_cung_day_du"])), (id_, min(Lk), b["phat_cung_day_du"])
    KQ1 = {"thieu": sorted(k for k, v in J.items() if v is None), "o": {}}
    f_ = os.path.join(O, "b2", "chon.json"); KQ1["B2"] = json.load(open(f_)) if os.path.exists(f_) else None
    for o in ("D5", "D4A2"):
        if KQ1["B2"] and o in KQ1["B2"]:
            KQ1["B2"][o]["danh_gia_phat"] = {i: tv([J.get(f"b2_{o}_c{i:02d}_{sd}", {}) and J[f"b2_{o}_c{i:02d}_{sd}"]["danh_gia_phat"] for sd in (6401, 6402, 6403)]) for i in range(16)}
            ch_ = KQ1["B2"][o]["chon"]; g0, gc_ = (KQ1["B2"][o]["danh_gia_phat"].get(i) for i in (0, ch_))
            KQ1["B2"][o]["co_chi_phi_tren_1p2x"] = bool(g0 and gc_ and gc_["tv"] > 1.2 * g0["tv"])
            KQ1["B2"][o]["giay_dat"] = {i: tv([J.get(f"b2_{o}_c{i:02d}_{sd}", {}) and J[f"b2_{o}_c{i:02d}_{sd}"]["giay_dat"] for sd in (6401, 6402, 6403)]) for i in range(16)}
    for o, K1, seeds, ms, mt in (("D5", 6, range(6101, 6106), (1, 4), (1, 4, 16)), ("D4A2", 16, range(6101, 6104), (1, 4), (0.5, 1, 4))):
        diem = collections.defaultdict(dict)
        for sd in seeds:
            j = J.get(f"e1_{o}_td_K16_{sd}")
            if j and len(j["Lk"]) == 16:
                for k in (1, 2, 4, 8, 16):
                    diem[("PLACE", k / 4)][sd] = tot_dau(j["Lk"], k) / j["N"]
                    if len(j["Lg"]) == 16: diem[("goc:PLACE", k / 4)][sd] = tot_dau_goc(j["Lk"], j["Lg"], k) / j["N"]
            j = J.get(f"e1_{o}_nn_K{4 * K1}_{sd}")
            if j and len(j["Lk"]) == 4 * K1:
                for k in sorted({1, max(1, K1 // 4), K1 // 2, K1, 2 * K1, 4 * K1}): diem[("thu_loai", k / K1)][sd] = tot_dau(j["Lk"], k) / j["N"]
            for pp, tien, mm in (("SA", "un_m", ms), ("SA_chinh", "unt_m", mt)):
                for m in mm:
                    j = J.get(f"e1_{o}_{tien}" + f"{m:g}".replace(".", "p") + f"_{sd}")
                    if not j: continue
                    diem[(pp, m)][sd] = j["P"] / j["N"]
                    if len(j["Lg"]) == 4: diem[("goc:" + pp, m)][sd] = tot_dau_goc(j["Lk"], j["Lg"], 4) / j["N"]
                    if len(j["Lk"]) == 4:
                        for k in (1, 2, 4): diem[(f"{pp}_khoi_lai_m{m:g}", m * k / 4)][sd] = tot_dau(j["Lk"], k) / j["N"]
        bang = {}
        for (pp, b), v in sorted(diem.items()):
            ys = [v[s] for s in sorted(v)]
            bang[f"{pp}@{b:g}"] = {"pp": pp, "b": b, "y": tv(ys), "theo_seed": {str(s): v[s] for s in sorted(v)},
                                    "n_duoi_tau": sum(y <= TAU for y in ys), "n_duoi_tau_lo": sum(y <= TAU_LO for y in ys)}
        def tvy(pp, b): x = bang.get(f"{pp}@{b:g}"); return None if not x or not x["y"] else x["y"]["tv"]
        def b_sao(pp, muc):
            if muc is None: return None
            bs = sorted(x["b"] for x in bang.values() if x["pp"] == pp and x["y"])
            for b in bs:
                if tvy(pp, b) <= muc + SAN: return b
            return "> %g" % bs[-1] if bs else None
        r = {"bang": bang}
        for ref in (1, 0.25):
            mu = tvy("PLACE", ref)
            r[f"b_sao_khop_PLACE_{ref:g}x"] = {pp: b_sao(pp, mu) for pp in ("SA", "SA_chinh", "thu_loai")}
        def cap(ppA, bA, ppB, bB):
            a = bang.get(f"{ppA}@{bA:g}", {}).get("theo_seed", {}); b = bang.get(f"{ppB}@{bB:g}", {}).get("theo_seed", {})
            ss = sorted(set(a) & set(b)); th = sum(a[s] > b[s] for s in ss); tu = sum(a[s] < b[s] for s in ss)
            return {"n": len(ss), "A_lon_hon": th, "A_nho_hon": tu, "p_dau_mot_phia": PT.dau_chinh_xac(th, tu),
                    "log10_ti_tv": float(np.median([math.log10((a[s] + SAN) / (b[s] + SAN)) for s in ss])) if ss else None}
        r["cap"] = {"SAc4_vs_PLACE1": cap("SA_chinh", 4, "PLACE", 1), "SAc16_vs_PLACE1": cap("SA_chinh", 16, "PLACE", 1), "SA4_vs_PLACE1": cap("SA", 4, "PLACE", 1),
                    "SAc1_vs_PLACE025": cap("SA_chinh", 1, "PLACE", 0.25), "SAc1_vs_SA1": cap("SA_chinh", 1, "SA", 1), "SAc4_vs_SA4": cap("SA_chinh", 4, "SA", 4),
                    "TL4_vs_TL1": cap("thu_loai", 4, "thu_loai", 1), "SAc4_vs_PLACE025": cap("SA_chinh", 4, "PLACE", 0.25)}
        if o == "D5":
            c1 = r["cap"]["SAc4_vs_PLACE1"]
            r["C1"] = {"gia_thuyet": "y_SA_da_chinh(4x) > y_PLACE(1x) theo seed (một phía), D5, dấu chính xác, alpha 0,05; thước = phần dư P (không phải G_c)",
                       "n": c1["n"], "thang": c1["A_lon_hon"], "p": c1["p_dau_mot_phia"], "dat": bool(c1["n"] == 5 and c1["p_dau_mot_phia"] <= 0.05)}
            g1 = cap("goc:SA_chinh", 4, "goc:PLACE", 1)
            r["C1_thuoc_goc"] = {"n": g1["n"], "thang": g1["A_lon_hon"], "p": g1["p_dau_mot_phia"], "ghi_chu": "CÙNG cấu hình (chọn theo thước sửa) đo bằng thước gốc — không phải điều bản cài đặt cũ sẽ tạo ra"}
        tl1, tl4 = tvy("thu_loai", 1), tvy("thu_loai", 4)
        r["du_doan"] = {"P2_SAc4_chua_khop_PLACE1": None if None in (tvy("SA_chinh", 4), tvy("PLACE", 1)) else tvy("SA_chinh", 4) > tvy("PLACE", 1) + SAN,
                        "P3_thu_loai_bao_hoa": (tl4 >= 0.5 * tl1) if (tl1 and tl4) else None,
                        "P4_PLACE025_nho_hon_SAc1": None if None in (tvy("PLACE", 0.25), tvy("SA_chinh", 1)) else tvy("PLACE", 0.25) < tvy("SA_chinh", 1),
                        "P5_SAc1_khong_kem_SA1": None if None in (tvy("SA_chinh", 1), tvy("SA", 1)) else tvy("SA_chinh", 1) <= tvy("SA", 1) + SAN}
        r["giay_dat_dinh_huong"] = {k: tv([J[k2]["giay_dat"] for k2 in J if J[k2] and k2.startswith(f"e1_{o}_{k}_")]) for k in
                                    ["td_K16", f"nn_K{4 * K1}"] + [("un_m%g" % m).replace(".", "p") for m in ms] + [("unt_m%g" % m).replace(".", "p") for m in mt]}
        r["khong_dat_thu_loai"] = tv([J[k2]["khong_dat"] for k2 in J if J[k2] and k2.startswith(f"e1_{o}_nn_")])
        KQ1["o"][o] = r
    e2 = {}
    for N in (440, 576, 724, 900):
        dd = {}
        for nh in ("td", "km", "deu"):
            v = {sd: J[f"e2_A{N:03d}_{nh}_{sd}"]["P"] / J[f"e2_A{N:03d}_{nh}_{sd}"]["N"] for sd in range(6301, 6304) if J.get(f"e2_A{N:03d}_{nh}_{sd}")}
            dd[nh] = {"y": tv(list(v.values())), "theo_seed": {str(k): x for k, x in v.items()}, "n_duoi_tau": sum(x <= TAU for x in v.values()), "n_duoi_tau_lo": sum(x <= TAU_LO for x in v.values()),
                      "S_mem": tv([J[f"e2_A{N:03d}_{nh}_{sd}"]["S_mem"] for sd in range(6301, 6304) if J.get(f"e2_A{N:03d}_{nh}_{sd}")])}
        for nh in ("km", "deu"):
            a, b = dd[nh]["theo_seed"], dd["td"]["theo_seed"]; ss = sorted(set(a) & set(b))
            dd[nh]["log10_ti_voi_td"] = [math.log10((a[s] + SAN) / (b[s] + SAN)) for s in ss]
            dd[nh]["kem_hon_td_so_seed"] = sum(a[s] > b[s] for s in ss)
        e2[f"A{N:03d}"] = dd
    KQ1["E2"] = e2
    KQ1["du_doan_E2"] = {"deu_kem_hon_td_moi_muc": all(e2[k]["deu"]["kem_hon_td_so_seed"] == len(e2[k]["deu"]["log10_ti_voi_td"]) > 0 for k in e2),
                         "km_trong_2x_td_moi_muc": all(all(abs(x) <= math.log10(2) for x in e2[k]["km"]["log10_ti_voi_td"]) for k in e2)}
    KQ1["md5_ma"] = md5(os.path.abspath(__file__))
    os.makedirs(RA, exist_ok=True); json.dump(KQ1, open(os.path.join(RA, "E12.json"), "w"), indent=1, ensure_ascii=False)
    L = ["# E1/E2 — kết quả tự động (phan_tich_tdg.py e12; CHƯA rà soát chéo)\n", f"thiếu: {KQ1['thieu']}\n"]
    for o, r in KQ1["o"].items():
        L.append(f"\n## E1 {o}\n\n| phương pháp | b | y trung vị | min | max | n | ≤τ | ≤τ_lo |\n|---|---|---|---|---|---|---|---|")
        for k, x in r["bang"].items():
            if x["y"]: L.append(f"| {x['pp']} | {x['b']:g} | {x['y']['tv']:.3g} | {x['y']['min']:.3g} | {x['y']['max']:.3g} | {x['y']['n']} | {x['n_duoi_tau']} | {x['n_duoi_tau_lo']} |")
        L.append(f"\nb* khớp PLACE 1×: {r['b_sao_khop_PLACE_1x']}; khớp PLACE 0,25×: {r['b_sao_khop_PLACE_0.25x']}\n\ncặp: {json.dumps(r['cap'], ensure_ascii=False)}\n")
        if "C1" in r: L.append(f"\n**C1**: {r['C1']}\n")
        L.append(f"\ndự đoán: {r['du_doan']}\n")
    L.append("\n## E2\n")
    for k, dd in e2.items():
        L.append(f"- {k}: " + "; ".join(f"{nh} y_tv {dd[nh]['y']['tv']:.3g} (≤τ {dd[nh]['n_duoi_tau']})" for nh in dd if dd[nh]["y"]))
    L.append(f"\ndự đoán E2: {KQ1['du_doan_E2']}\n")
    open(os.path.join(RA, "E12.md"), "w").write("\n".join(L)); print("\n".join(L))

def gc_B(tt):
    doi = tt.get("doi_mat_um") or {}; g = 0.0; n_none = 0
    for ten in tt["tu_the"]:
        v = doi.get(ten)
        if not v or v[1] is None: n_none += 1; continue
        g += max(v[1] - H, 0.0)
    return g, n_none
def so_bd(tt, P=None):
    g, nn = gc_B(tt); doi = tt.get("doi_mat_um") or {}
    return {"P": P, "n_vp": tt.get("diem_vi_pham_ban_dau"), "n_kkd": tt.get("diem_khong_kha_thi_truoc_chieu"), "Gc_B": g, "Gc_B_khong_co": nn,
            "n_mat": len(tt.get("mat_bo_phan") or []), "n_duoi50": len(tt.get("con_lai_duoi_50pt") or []),
            "n_doi_tren_h": sum(1 for v in doi.values() if v and v[1] is not None and v[1] > H),
            "tong_doi_max": sum(v[1] for v in doi.values() if v and v[1] is not None),
            "V_mat_ti_le": (lambda vt: sum(max(a - b, 0) for a, b in vt.values()) / max(sum(a for a, b in vt.values()), 1e-12))(tt.get("the_tich_truoc_sau") or {})}
def gc_K_kiem(rs, tu_the):
    g = [r for r in rs if r["luat"] == "G_c"]; assert len(g) == 1; g = g[0]
    Gc = float(g["trang_thai"]); Di = [r for r in rs if r["luat"] == "D_i"]
    mat = [(r["phan_quyet"], r["bo_phan_truot"]) for r in Di if r["trang_thai"] == "KHONG_CO_MESH"]; Rm = []; nai = 0
    for ten, kh in mat:
        assert tu_the[ten]["khuon"] == kh
        try: Rm.append(PT.ban_kinh_bao(kh, tu_the[ten]["s"]))
        except ValueError: nai += 1
    st = {r["luat"]: r for r in rs if r["phan_quyet"] == "su_that" and r["loai_luat"] == "su_that_v1.2"}
    tr = {k: int(v["so_ca_truot"]) for k, v in st.items()}; pen = 0
    for ten, kh in mat:
        lo = ten.rstrip("0123456789").rstrip("_"); pen += len(PT.LUAT_MAT.get(lo, PT.LUAT_MAT_MAC_DINH))
    return {"Gc_K": float("nan") if nai else Gc + sum(max(x - H, 0) for x in Rm), "n_mat_K": len(mat), "n_mat_ai_khong_R": nai,
            "V_test": sum(v for k, v in tr.items() if k not in PT.LUAT_MEM) + pen}
def doc_kiem(lo):
    C = os.path.join(O, "cham_" + lo); KIEM = {}
    if not os.path.exists(os.path.join(C, "khoa", "KHOA_MU.json")): return KIEM
    khoa = json.load(open(os.path.join(C, "khoa", "KHOA_MU.json")))
    for l in open(os.path.join(C, "vs", "ds_sau.txt")):
        if not l.strip() or l.startswith("#"): continue
        d = l.split()[0].rstrip("/"); idm = os.path.basename(d); cha = os.path.basename(os.path.dirname(d))
        f = os.path.join(os.environ.get("CHAM_RA_TDG", os.path.join(GOC, "bo_kiem_v3", "cham_ra")), f"{cha}__{idm}", "kiem3.csv")
        if not os.path.exists(f): continue
        tu = json.load(open(os.path.join(C, "mu_sau", idm, "thong_tin.json")))["tu_the"]
        KIEM[os.path.basename(khoa["id"][idm]["nguon"].rstrip("/"))] = gc_K_kiem(list(csv.DictReader(open(f, encoding="utf-8"))), tu)
    return KIEM
def tach(job):
    p = job.split("_"); o = p[0]; sd = int(p[-1]); mid = p[1:-1]
    return o, ("_".join(mid[:-1]) + "_") if len(mid) > 1 else "", mid[-1], sd
def dem(a, b, seeds=None):
    if seeds is None:
        ss = sorted(set(a) & set(b)); d = [a[s] - b[s] for s in ss if not (math.isnan(a[s]) or math.isnan(b[s]))]
    else:
        ok = [s for s in seeds if s in a and s in b and not (math.isnan(a[s]) or math.isnan(b[s]))]; d0 = [a[s] - b[s] for s in ok]
        M = max([abs(x) for x in d0] or [0.0]) + 1.0; d = d0 + [-M] * (len(seeds) - len(ok))
    if not d: return None
    p, W, n_, hoa = PT.wilcoxon_chinh_xac(d)
    return {"n": len(d), "thang": sum(x > 1e-6 for x in d), "thua": sum(x < -1e-6 for x in d), "hoa": hoa, "p_wilcoxon": p, "d_tv": float(np.median(d))}
def e3():
    KIEM = {**doc_kiem("e3"), **doc_kiem("sp")}
    MO = {r["job"]: r for r in csv.DictReader(open(BANG_MO, encoding="utf-8")) if r["ban"] == "sau"}
    out = {"E3a": {}, "E3b": {}, "SP": {}, "thieu_cham": [], "thieu_dung": []}
    viec = [l.split()[0] for l in open(os.path.join(O, "viec.txt")) if l.strip()]
    def doc(id_):
        f = os.path.join(KQ, id_, "thong_tin.json")
        return json.load(open(f)) if os.path.exists(f) else None
    def gan(r, v, id_):
        if r.get(v) is None: return
        if id_ in KIEM: r[v].update(KIEM[id_])
        else: out["thieu_cham"].append(id_)
    hang = []
    for id_ in viec:
        if not id_.startswith("e3a_B_"): continue
        job = id_[6:]; o, tien, nh, sd = tach(job)
        tB = doc(id_); tA = doc("e3a_A_" + job) if sd == 4001 else None; tO = json.load(open(os.path.join(KQ_CHINH, job, "thong_tin.json")))
        if tB is None: out["thieu_dung"].append(id_)
        r = {"job": job, "o": o, "tien": tien, "nhanh": nh, "seed": sd}
        r["goc"] = so_bd(tO, tO["bao_cao_bo_tri"]["phat_cung_day_du"])
        r["goc"].update({"Gc_K": float(MO[job]["Gc_K"] or "nan"), "V_test": int(MO[job]["V_test"]), "n_mat_K": int(MO[job]["n_mat"])})
        r["B"] = so_bd(tB, tB["ban_e3"]["P_dung"]) if tB else None
        r["A"] = so_bd(tA, tA["ban_e3"]["P_goc"]) if tA else None
        if tB: r["B"]["P_goc_cung_tu_the"] = tB["ban_e3"]["P_goc"]
        gan(r, "B", id_)
        if sd == 4001: gan(r, "A", "e3a_A_" + job)
        hang.append(r)
    tA = [r for r in hang if r["A"]]
    out["E3a"]["phat_lai"] = {"n": len(tA), "n_vp_bang": sum(r["A"]["n_vp"] == r["goc"]["n_vp"] for r in tA), "n_mat_bang": sum(r["A"]["n_mat"] == r["goc"]["n_mat"] for r in tA),
        "max_dGcB": max([abs(r["A"]["Gc_B"] - r["goc"]["Gc_B"]) for r in tA] or [None]),
        "max_dP_tuong_doi": max([abs(r["A"]["P"] - r["goc"]["P"]) / max(abs(r["goc"]["P"]), 1e-9) for r in tA] or [None]),
        "dGcK_A_goc": {r["job"]: r["A"]["Gc_K"] - r["goc"]["Gc_K"] for r in tA if "Gc_K" in r["A"]}}
    pl = out["E3a"]["phat_lai"]
    pl["H0_dat"] = bool(pl["n"] == 36 and len([x for x in pl["dGcK_A_goc"].values() if not math.isnan(x)]) == 36 and all(abs(r["A"]["n_vp"] - r["goc"]["n_vp"]) <= 1e-4 * r["goc"]["n_vp"] for r in tA) and pl["n_mat_bang"] == pl["n"]
                        and (pl["max_dGcB"] or 0) <= 1e-3 and all(abs(r["A"]["P"] - r["goc"]["P"]) <= 1e-5 * r["goc"]["P"] + 1e-9 for r in tA)
                        and all(abs(x) <= 0.01 for x in pl["dGcK_A_goc"].values() if not math.isnan(x)))
    nhom = collections.defaultdict(list)
    for r in hang:
        if r["B"]: nhom[(r["o"], r["tien"], r["nhanh"])].append(r)
    lech = {}
    for (o, tien, nh), rs in sorted(nhom.items()):
        e = {}
        for m in ("Gc_K", "Gc_B", "P", "n_mat", "n_mat_K", "n_vp", "n_kkd", "n_doi_tren_h", "V_test", "V_mat_ti_le"):
            dd = [r["B"][m] - r["goc"][m] for r in rs if r["B"].get(m) is not None and r["goc"].get(m) is not None
                  and not (isinstance(r["B"][m], float) and math.isnan(r["B"][m])) and not (isinstance(r["goc"][m], float) and math.isnan(r["goc"][m]))]
            e[m] = {"goc": tv([r["goc"][m] for r in rs if r["goc"].get(m) is not None]), "B": tv([r["B"][m] for r in rs if r["B"].get(m) is not None]), "hieu_B_tru_goc": tv(dd)}
        e["P_gia_ti_le"] = tv([(r["B"]["P_goc_cung_tu_the"] - r["B"]["P"]) / max(r["B"]["P_goc_cung_tu_the"], 1e-12) for r in rs if r["B"]["P_goc_cung_tu_the"] > 0])
        lech[f"{o}_{tien}{nh}"] = e
    out["E3a"]["lech"] = lech
    ss = {}; tl = {}
    for (o, tien) in sorted({(r["o"], r["tien"]) for r in hang}):
        for dc in ("nn1x", "un1x", "nn4x", "un4x"):
            for m in ("Gc_K", "Gc_B"):
                for v in ("goc", "B"):
                    a = {r["seed"]: r[v][m] for r in hang if (r["o"], r["tien"], r["nhanh"]) == (o, tien, dc) and r.get(v) and m in r[v]}
                    b = {r["seed"]: r[v][m] for r in hang if (r["o"], r["tien"], r["nhanh"]) == (o, tien, "td") and r.get(v) and m in r[v]}
                    nseed = 10 if (tien == "" and dc in ("nn1x", "un1x")) else 5
                    x = dem(a, b, list(range(4001, 4001 + nseed)))
                    if x: ss[f"{o}_{tien}{dc}_{m}_{v}"] = x
                R_ = {(r["nhanh"], r["seed"]): r for r in hang if (r["o"], r["tien"]) == (o, tien) and r["B"]}
                xx = []
                for s_ in range(4001, 4011):
                    a, b = R_.get((dc, s_)), R_.get(("td", s_))
                    if not a or not b or m not in a["B"] or m not in b["B"]: continue
                    d0 = a["goc"][m] - b["goc"][m]; d1 = a["B"][m] - b["B"][m]
                    if not (math.isnan(d0) or math.isnan(d1)): xx.append((d0 - d1, d1 / d0 if abs(d0) > 1e-9 else None))
                if xx: tl[f"{o}_{tien}{dc}_{m}"] = {"d_goc_tru_d_B": tv([u for u, _ in xx]), "ti_d_B_tren_d_goc": tv([w for _, w in xx if w is not None])}
    out["E3a"]["so_sanh"] = ss; out["E3a"]["thien_lech_PLACE"] = tl
    def giu(k):
        a, b = ss.get(k + "_Gc_K_B"), ss.get(k + "_Gc_K_goc")
        return None if not a or not b else bool(a["p_wilcoxon"] <= 0.025 and a["thang"] >= b["thang"] - 1)
    out["E3a"]["ket_luan_giu"] = {k: (giu(k) if pl["H0_dat"] else "CHI_MO_TA (H0 phát lại trượt)") for k in ("D5_nn1x", "D5_un1x", "D4A2_nn1x")}
    out["E3a"]["loi_phinh_uu_the_PLACE"] = {k: (v["d_goc_tru_d_B"]["tv"] > 0 and v["ti_d_B_tren_d_goc"] is not None and v["ti_d_B_tren_d_goc"]["tv"] < 0.9)
                                           for k, v in tl.items() if k.endswith("_Gc_K") and v["d_goc_tru_d_B"]}
    def doi_dem(k):
        a, b = ss.get(k + "_Gc_K_B"), ss.get(k + "_Gc_K_goc")
        return None if not a or not b else {"thang": [b["thang"], a["thang"]], "thua": [b["thua"], a["thua"]], "doi_hon_1": abs(a["thang"] - b["thang"]) > 1 or abs(a["thua"] - b["thua"]) > 1}
    out["E3a"]["D4_D1_khe_dem"] = {k: doi_dem(k) for k in ["D4_nn1x", "D4_un1x", "D1_nn1x", "D1_un1x", "D4A2_un1x"] + [f"{o}_khe{x}_{dc}" for o in ("D5", "D4A2") for x in ("0.1", "0.5") for dc in ("nn1x", "un1x")]}
    hb = []
    for o in ("D5", "D4A2", "D4"):
        for sd in range(6201, 6206):
            for nh in ("td", "nn1x", "un1x"):
                r = {"o": o, "seed": sd, "nhanh": nh}
                for v in ("A", "B", "R"):
                    id_ = f"e3b_{v}_{o}_{nh}_{sd}"; t_ = doc(id_)
                    r[v] = so_bd(t_, (t_["ban_e3"]["P_dung"] if v == "R" else t_["bao_cao_bo_tri"]["phat_cung_day_du"])) if t_ else None
                    if t_ is None: out["thieu_dung"].append(id_)
                    gan(r, v, id_)
                hb.append(r)
    eb = {}
    for o in ("D5", "D4A2", "D4"):
        for nh in ("td", "nn1x", "un1x"):
            rs = [r for r in hb if r["o"] == o and r["nhanh"] == nh and r["A"] and r["B"]]
            eb[f"{o}_{nh}"] = {m: {"A": tv([r["A"].get(m) for r in rs]), "B": tv([r["B"].get(m) for r in rs]), "R": tv([r["R"].get(m) for r in rs if r.get("R")]),
                                   "hieu_B_tru_A": tv([r["B"][m] - r["A"][m] for r in rs if r["A"].get(m) is not None and r["B"].get(m) is not None]),
                                   "do_dat_B_tru_R": tv([r["B"][m] - r["R"][m] for r in rs if r.get("R") and r["R"].get(m) is not None and r["B"].get(m) is not None]),
                                   "do_chieu_R_tru_A": tv([r["R"][m] - r["A"][m] for r in rs if r.get("R") and r["R"].get(m) is not None and r["A"].get(m) is not None])}
                               for m in ("Gc_K", "Gc_B", "P", "n_mat", "n_mat_K", "n_kkd", "V_test")}
        for v in ("A", "B", "R"):
            for dc in ("nn1x", "un1x"):
                for m in ("Gc_K", "Gc_B", "P"):
                    a = {r["seed"]: r[v][m] for r in hb if r["o"] == o and r["nhanh"] == dc and r[v] and m in r[v]}
                    b = {r["seed"]: r[v][m] for r in hb if r["o"] == o and r["nhanh"] == "td" and r[v] and m in r[v]}
                    x = dem(a, b)
                    if x: eb[f"{o}_PLACE_vs_{dc}_{m}_{v}"] = {k: x[k] for k in ("n", "thang", "thua", "d_tv")}
    for o in ("D5", "D4A2", "D4"):
        for nh in ("td", "nn1x", "un1x"):
            xs = []
            for sd in range(6201, 6206):
                a_, b_ = doc(f"e3b_A_{o}_{nh}_{sd}"), doc(f"e3b_B_{o}_{nh}_{sd}")
                if a_ and b_: xs.append(b_["bao_cao_bo_tri"]["giay_dat"] / max(a_["bao_cao_bo_tri"]["giay_dat"], 1e-9))
            eb[f"{o}_{nh}_ti_gio_dat_B_tren_A"] = tv(xs)
    out["E3b"] = eb
    hs = []
    for N in (440, 576, 724, 900):
        for nh in ("td", "nn1x", "un1x"):
            for sd in (7101, 7102, 7103):
                q = f"A{N:03d}_{nh}_{sd}"; r = {"muc": N, "nhanh": nh, "seed": sd}
                fq = os.path.join(GOC, "bo_dung_chung", "quet_mat_do", "kq", q, "thong_tin.json")
                r["y_quet"] = (lambda d: d["bao_cao_bo_tri"]["phat_cung_day_du"] / sum(d["so_luong"].values()))(json.load(open(fq))) if os.path.exists(fq) else None
                for v in ("A", "B"):
                    id_ = f"sp_{v}_{q}"; t_ = doc(id_)
                    r[v] = so_bd(t_, t_["ban_e3"]["P_dung" if v == "B" else "P_goc"]) if t_ else None
                    if t_ is None: out["thieu_dung"].append(id_)
                    gan(r, v, id_)
                hs.append(r)
    sp = {}
    for N in (440, 576, 724, 900):
        for nh in ("td", "nn1x", "un1x"):
            rs = [r for r in hs if r["muc"] == N and r["nhanh"] == nh]
            sp[f"A{N:03d}_{nh}"] = {"y_quet": tv([r["y_quet"] for r in rs]), **{f"{m}_{v}": tv([r[v].get(m) for r in rs if r[v]]) for m in ("Gc_K", "Gc_B", "n_mat_K", "n_kkd") for v in ("A", "B")}}
        for v in ("A", "B"):
            for dc in ("nn1x", "un1x"):
                a = {r["seed"]: r[v]["Gc_K"] for r in hs if r["muc"] == N and r["nhanh"] == dc and r[v] and "Gc_K" in r[v]}
                b = {r["seed"]: r[v]["Gc_K"] for r in hs if r["muc"] == N and r["nhanh"] == "td" and r[v] and "Gc_K" in r[v]}
                x = dem(a, b)
                if x: sp[f"A{N:03d}_PLACE_vs_{dc}_Gc_K_{v}"] = {k: x[k] for k in ("n", "thang", "thua", "d_tv")}
    for v in ("A", "B"):
        cc = [(r["y_quet"] > TAU, r[v]["Gc_K"] >= 1.0) for r in hs if r["y_quet"] is not None and r[v] and "Gc_K" in r[v] and not math.isnan(r[v]["Gc_K"])]
        sp[f"tau_ngoai_mau_{v}"] = {"ghi_chu": "y_quét đo bằng thước gốc; G_c dưới " + ("eval gốc" if v == "A" else "phép chiếu đã sửa"), "n": len(cc), "dung": sum(a == b for a, b in cc), "duong_that": sum(b for _, b in cc), "duong_du_doan": sum(a for a, _ in cc)}
        pr = [(r["y_quet"], r[v]["Gc_K"]) for r in hs if r["y_quet"] is not None and r[v] and "Gc_K" in r[v] and not math.isnan(r[v]["Gc_K"])]
        sp[f"spearman_y_Gc_{v}"] = PT.spearman([a for a, _ in pr], [b for _, b in pr]) if len(pr) > 2 else None
    out["SP"] = sp
    out["hang_E3a"] = hang; out["hang_E3b"] = hb; out["hang_SP"] = hs; out["md5_ma"] = md5(os.path.abspath(__file__))
    os.makedirs(RA, exist_ok=True); json.dump(out, open(os.path.join(RA, "E3.json"), "w"), indent=1, ensure_ascii=False, default=str)
    L = ["# E3 + S-P — kết quả tự động (phan_tich_tdg.py e3; CHƯA rà soát chéo)\n", f"thiếu dựng: {out['thieu_dung']}\n\nthiếu chấm: {out['thieu_cham']}\n",
         f"phát lại: {json.dumps(out['E3a']['phat_lai'], ensure_ascii=False)}\n", "\n| nhóm | thước | gốc tv | B tv | B−gốc tv [min, max] |\n|---|---|---|---|---|"]
    for k, e in lech.items():
        for m in ("Gc_K", "Gc_B", "P", "n_mat_K", "n_kkd", "V_test"):
            x = e[m]
            if x["hieu_B_tru_goc"] and x["goc"] and x["B"]: L.append(f"| {k} | {m} | {x['goc']['tv']:.4g} | {x['B']['tv']:.4g} | {x['hieu_B_tru_goc']['tv']:.4g} [{x['hieu_B_tru_goc']['min']:.4g}, {x['hieu_B_tru_goc']['max']:.4g}] |")
    L.append("\n## So sánh PLACE vs đối chứng (gốc / B)\n"); L += [f"- {k}: {v}" for k, v in ss.items()]
    L.append(f"\n## Thiên lệch có lợi cho PLACE\n{json.dumps(tl, ensure_ascii=False)}\n\n## Kết luận giữ: {out['E3a']['ket_luan_giu']}\n\nD4/D1/khe: {json.dumps(out['E3a']['D4_D1_khe_dem'], ensure_ascii=False)}\n")
    L.append(f"\n## E3b\n{json.dumps(eb, ensure_ascii=False)}\n\n## S-P\n{json.dumps(sp, ensure_ascii=False)}\n")
    open(os.path.join(RA, "E3.md"), "w").write("\n".join(L)); print("\n".join(L))

if __name__ == "__main__":
    {"e12": e12, "e3": e3}[sys.argv[1]]()
