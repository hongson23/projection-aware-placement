# -*- coding: utf-8 -*-
"""phan_tich_tv_u5.py BẢN 3 (22/11/2026; sửa theo RaSoat_LuotThucVat_U5_20261122 B2, B3, n3, n7 — TRƯỚC lượt chính, trước mọi kết quả chính).
LUẬT ĐỌC GHI TRƯỚC (đóng băng md5 trong nhat_ky_hang.txt) — bảng kiểu CG cho W4 thực vật (TV-D1, TV-D5, TV-real × 5 nhánh × seed).
Danh sách id KỲ VỌNG = mức × nhánh × seed (chính: u5_tv_c_<mức>_<nhánh>_<seed>, seed 7311–7320 hoặc --seeds; --dev: các id m1 trong viec_s0.txt).
Trạng thái mỗi id: ok | hong_dung (loi/<id>, hết giờ rc=124 ghi riêng, hoặc thiếu thong_tin/n2) | hong_cham (có n2 nhưng thiếu kiem3.csv).
Thước CHÍNH (N2 cùng lưới trước ↔ sau PROJECT): HDmax = max HD hai chiều trên MỌI cấu trúc (bộ phận rắn ∪ ER ∪ không bào);
  việc có bộ phận MẤT ⇒ HDmax = +∞ (B3); HDp95 = phân vị 95 qua bộ phận rắn (kèm p95 chỉ trên bộ phận bị chạm, mô tả);
  %V mất = Σ max(V0 − V1, 0)/Σ V0 (rắn) và theo loại; MẤT, BẸP (< 50 %) đếm bộ phận.
Thước PHỤ (bộ kiểm v3 qua bản mù): G_c^K tại h = 0,05 và 2h; luật 2.7b, 2.7f, M.d theo TỈ LỆ ca trượt; luật sự thật trượt Σ (trừ 2.7b/2.7f/M.d;
  2.7c dưới độ phân giải và 2.3a mặt chung tách cột riêng); khả thi sau PROJECT; ngân sách đếm (K′ / vòng SA / Adam), số lần đánh giá; giờ đặt (MÔ TẢ).
Việc hỏng (B2): vào trung vị với giá trị TỆ NHẤT (HD = +∞, %V = 100, MẤT = số bộ phận yêu cầu, tỉ lệ trượt = 1, G_c = +∞, luật = —);
  thắng/hoà/thua của U5+KTT: U5+KTT hỏng ⇒ L; đối chứng hỏng mà U5+KTT không ⇒ W; cả hai hỏng ⇒ T và đếm riêng "cả hai hỏng". Thước của bộ kiểm:
  hong_cham cũng là hỏng. Hoà: |Δ| < 1e-4 µm (HD), < 1e-3 điểm % (thể tích), bằng nhau (đếm).
Lời nhận "U5+KTT tốt hơn X ở mức L trên thước m" ⇔ thắng > thua VÀ trung vị cùng chiều; "rõ" ⇔ thắng ở mọi seed (6/6; kiểm định dấu 6/6 cho p = 0,03). Mức mà mọi nhánh có HDmax < h ⇒ chỉ nói
  "không gây hại". Không gộp các mức.
Dùng: python phan_tich_tv_u5.py [--dev] [--seeds 7311,7312] [--ra <tệp.md>]   (chạy từ ch3_tebao_dep)"""
import os, sys, json, glob, csv, math, statistics as st, collections, hashlib, re
R = os.path.expanduser("~/Vision_team_data_189/sonhh20/ch3_tebao_dep"); os.chdir(R)
O = "bo_dung_chung/tang_dg"; U = f"{O}/u5tv"; KV3 = "bo_kiem_v3"; INF = float("inf")
DEV = "--dev" in sys.argv
SEEDS = [int(x) for x in sys.argv[sys.argv.index("--seeds") + 1].split(",")] if "--seeds" in sys.argv else list(range(7311, 7317))
RA = sys.argv[sys.argv.index("--ra") + 1] if "--ra" in sys.argv else f"{U}/{'KET_QUA_DEV' if DEV else 'KET_QUA_CHINH'}.md"
MUC = ("D1", "D5", "real"); NH = ("u5ktt", "td", "nn1x", "nn1xm", "un1x")
TEN = {"u5ktt": "U5+KTT", "td": "PLACE (P1)", "nn1x": "thử–loại", "nn1xm": "thử–loại chọn J mềm", "un1x": "SA (chỉnh B2)"}
MEM = ("2.7b", "2.7f", "M.d"); RIENG = ("2.7c", "2.3a"); LOAI = ("luc_lap", "ti_the", "nhan", "golgi", "peroxisome"); H = 0.05
if DEV: KY = [l.split()[0] for l in open(f"{U}/viec_s0.txt") if l.startswith("u5_tv_m1_")]
else: KY = [f"u5_tv_c_{m}_{nh}_{s}" for s in SEEDS for m in MUC for nh in NH
           if (m == "D1" and s in (7311, 7312, 7313)) or (m != "D1" and nh != "nn1xm")]
nk = open(f"{O}/nhat_ky.txt", errors="replace").read()
mu = {}
for f in glob.glob(f"{U}/cham/lo_*/mu/khoa/KHOA_MU.json"):
    lo = f.split("/cham/lo_")[1].split("/")[0]
    for m, v in json.load(open(f))["id"].items():
        if v["ban"] == "sau": mu[os.path.basename(v["nguon"].rstrip("/"))] = (lo, m)
def k3(i):
    if i not in mu: return None
    lo, m = mu[i]; f = f"{KV3}/cham_ra/vsu5tv{lo}_mu_sau__{m}/kiem3.csv"
    return list(csv.DictReader(open(f, encoding="utf-8"))) if os.path.exists(f) else None
def mot(i):
    t = f"{O}/kq/{i}/thong_tin.json"; n = f"{O}/kq/{i}/n2.json"
    het_gio = bool(re.search(rf"LOI(_LAN1)? {re.escape(i)} rc=124", nk))
    if os.path.exists(f"{O}/loi/{i}") or not (os.path.exists(t) and os.path.exists(n)):
        return {"tt": "hong_dung", "het_gio": het_gio}
    T = json.load(open(t)); N = json.load(open(n)); b = T["bao_cao_bo_tri"]; r = {"tt": "ok", "het_gio": het_gio}
    nyc = sum((T.get("so_luong") or {}).values()); r["n_yeu_cau"] = nyc
    r.update({k: N[k] for k in ("hd_max_tat_ca", "hd_p95_ran", "pt_mat_ran", "mat", "duoi50", "lech_bo_dung")})
    if N["mat"] > 0: r["hd_max_tat_ca"] = INF
    cham = [max(e["hd"]) for e in N["bo_phan"].values() if e.get("hd") and max(e["hd"]) > 0]
    r["hd_p95_cham"] = sorted(cham)[int(0.95 * (len(cham) - 1))] if cham else 0.0
    r["pt_loai"] = {l: (N["theo_loai"].get(l) or {}).get("pt_mat") for l in LOAI}
    r.update({"kha_thi": T.get("kha_thi"), "giay_dat": b.get("giay_dat"), "danh_gia": b.get("danh_gia_phat"), "so_luong": T.get("so_luong")})
    ns = b.get("ngan_sach_nhanh") or {}
    r["ngan_sach"] = f"K {ns.get('K')}" if T.get("cach") == "ngau_nhien" else (f"vòng {ns.get('un_vong_lon')}/{ns.get('un_vong_nho')}" if T.get("cach") == "u_nhiet" else f"K {ns.get('K')}×Adam {ns.get('buoc')}")
    kt = (b.get("ktt") or {}); ck = [x for x in kt.get("chon_khoa", []) if x[0] == b.get("cau_hinh_chon")]
    r["ktt_nhan"] = bool(ck and ck[0][-1]) if kt.get("bat") else None
    u5 = b.get("u5")
    if u5: r["u5_khac_P"] = min(u5["cau_hinh"], key=lambda x: (x["P"], x["k"]))["k"] != b.get("cau_hinh_chon")
    rs = k3(i); cham_bp = set((T.get("ban_u5b") or {}).get("bo_phan_cham") or [])
    if rs is None: r["tt_k"] = "hong_cham"; return r
    r["tt_k"] = "ok"
    Di = [(x["phan_quyet"], x["trang_thai"]) for x in rs if x["luat"] == "D_i"]
    D = [float(v) for p, v in Di if v not in ("KHONG_CO_MESH", "", None)]
    r["mat_K"] = sum(1 for p, v in Di if v == "KHONG_CO_MESH")
    r["Gc_h"] = INF if r["mat_K"] else sum(max(x - H, 0) for x in D); r["Gc_2h"] = INF if r["mat_K"] else sum(max(x - 2 * H, 0) for x in D)
    r["san"] = [float(v) for p, v in Di if p not in cham_bp and v not in ("KHONG_CO_MESH", "", None)]
    st_ = [x for x in rs if x["loai_luat"] == "su_that_v1.2" and x["phan_quyet"] == "su_that"]
    for m in MEM + RIENG:
        x = [y for y in st_ if y["luat"] == m]; r[f"tl_{m}"] = (int(x[0]["so_ca_truot"]) / max(int(x[0]["so_ca"]), 1)) if x else None
    r["that_truot"] = sum(int(y["so_ca_truot"]) for y in st_ if y["luat"] not in MEM + RIENG)
    return r
XAU = {"hd_max_tat_ca": INF, "hd_p95_ran": INF, "pt_mat_ran": 100.0, "Gc_h": INF, "Gc_2h": INF}
def gt(r, m):
    K_ = m in ("Gc_h", "Gc_2h") or m.startswith("tl_")
    if r["tt"] != "ok" or (K_ and r.get("tt_k") != "ok"):
        return 1.0 if m.startswith("tl_") else (r.get("n_yeu_cau") or INF if m in ("mat", "duoi50") else XAU.get(m, INF))
    return r.get(m)
def kv(i): p = i.split("_"); return p[3], p[4], int(p[5])
D = {i: mot(i) for i in KY}
sys.path.insert(0, U); from chung_u5tv import REAL
def yc(i):
    m, nh, s = kv(i)
    if m == "D5": return {"luc_lap": 113, "ti_the": 15, "golgi": 8, "peroxisome": 30}
    if m == "real": return {"luc_lap": REAL[s][0], "ti_the": REAL[s][1], "golgi": 8, "peroxisome": 30}
    return None
def fm(x, n=3):
    if x is None: return "—"
    return "∞" if x == INF else f"{x:.{n}f}"
def mi(xs, n=3):
    xs = sorted(x for x in xs if x is not None)
    if not xs: return "—"
    q = lambda p: xs[min(len(xs) - 1, max(0, int(round(p * (len(xs) - 1)))))]
    if len(xs) < 5: return f"{fm(st.median(xs), n)} [{fm(xs[0], n)}–{fm(xs[-1], n)}]"
    return f"{fm(st.median(xs), n)} [{fm(q(0.25), n)}–{fm(q(0.75), n)}]"
TOL = {"hd_max_tat_ca": 1e-4, "hd_p95_ran": 1e-4, "pt_mat_ran": 1e-3, "mat": 0.5, "duoi50": 0.5, "Gc_h": 1e-4, "tl_2.7b": 1e-9, "tl_2.7f": 1e-9, "tl_M.d": 1e-9}
def hong(r, m): return r["tt"] != "ok" or ((m in ("Gc_h", "Gc_2h") or m.startswith("tl_")) and r.get("tt_k") != "ok")
def wtl(u, c, m):
    hu, hc = hong(u, m), hong(c, m)
    if hu and hc: return "B"
    if hu: return "L"
    if hc: return "W"
    a, b = u.get(m), c.get(m)
    if a is None or b is None: return "T"
    if a == b or abs(a - b) < TOL[m]: return "T"
    return "W" if a < b else "L"
out = [f"# {'DEV (M1)' if DEV else 'Lượt chính'} W4 thực vật U5 — bảng tự động (phan_tich_tv_u5.py md5 {hashlib.md5(open(__file__, 'rb').read()).hexdigest()[:8]}; seed {SEEDS[0]}–{SEEDS[-1]})", "", "Ô số: trung vị [IQR] khi n ≥ 5; [min–max] khi n < 5.", ""]
for muc in MUC:
    ids_m = [i for i in KY if kv(i)[0] == muc]
    if not ids_m: continue
    seeds = sorted({kv(i)[2] for i in ids_m})
    out += [f"## TV-{muc} (seed {', '.join(map(str, seeds))}; n = {len(seeds)})", "",
            "| nhánh | hỏng dựng/chấm | HDmax µm (mọi cấu trúc) | HDp95 µm | p95 bộ phận bị chạm | %V mất | MẤT Σ | BẸP <50 % Σ | khả thi | G_c^K h | G_c^K 2h | 2.7b | 2.7f | M.d | luật thật trượt Σ | 2.7c | 2.3a | ngân sách đếm | đánh giá | giờ đặt s (mô tả) |",
            "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for nh in NH:
        R_ = [D[i] for i in ids_m if kv(i)[1] == nh]
        if not R_: continue
        okk = [r for r in R_ if r["tt"] == "ok"]; okc = [r for r in okk if r.get("tt_k") == "ok"]
        out.append(f"| {TEN[nh]} | {sum(r['tt'] != 'ok' for r in R_)} ({sum(r['het_gio'] for r in R_ if r['tt'] != 'ok')} hết giờ)/{sum(r.get('tt_k') == 'hong_cham' for r in okk)} | "
                   f"{mi([gt(r, 'hd_max_tat_ca') for r in R_])} | {mi([gt(r, 'hd_p95_ran') for r in R_])} | {mi([r['hd_p95_cham'] for r in okk])} | {mi([gt(r, 'pt_mat_ran') for r in R_], 2)} | "
                   f"{sum(gt(r, 'mat') for r in R_) if all(r['tt'] == 'ok' for r in R_) else str(sum(r['mat'] for r in okk)) + ' + hỏng'} | {sum(r['duoi50'] for r in okk)} | {sum(1 for r in okk if r['kha_thi'])}/{len(R_)} | "
                   f"{mi([gt(r, 'Gc_h') for r in R_])} | {mi([gt(r, 'Gc_2h') for r in R_])} | " + " | ".join(mi([gt(r, f'tl_{m}') for r in R_], 3) for m in MEM) +
                   f" | {sum(r['that_truot'] for r in okc) if len(okc) == len(R_) else (str(sum(r['that_truot'] for r in okc)) + ' (thiếu ' + str(len(R_) - len(okc)) + ')')} | "
                   + " | ".join(mi([r.get(f'tl_{m}') for r in okc], 3) for m in RIENG) +
                   f" | {okk[0]['ngan_sach'] if okk else '—'} | {mi([r['danh_gia'] for r in okk], 0)} | {mi([r['giay_dat'] for r in okk], 1)} |")
    out += ["", "%V mất theo loại (trung vị [IQR]; chỉ việc dựng được):", "", "| nhánh | " + " | ".join(LOAI) + " |", "|---|" + "---|" * len(LOAI)]
    for nh in NH:
        okk = [D[i] for i in ids_m if kv(i)[1] == nh and D[i]["tt"] == "ok"]
        if okk: out.append(f"| {TEN[nh]} | " + " | ".join(mi([r['pt_loai'].get(l) for r in okk], 2) for l in LOAI) + " |")
    out += ["", "Thắng/hoà/thua của U5+KTT theo seed (cả hai hỏng đếm riêng ở ô cuối):", "",
            "| so với | HDmax | HDp95 | %V mất | MẤT | BẸP | G_c^K h | 2.7b | 2.7f | M.d | cả hai hỏng |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for nh in NH[1:]:
        if not any(kv(i)[1] == nh for i in ids_m): continue
        c = []; bb = 0
        for m in ("hd_max_tat_ca", "hd_p95_ran", "pt_mat_ran", "mat", "duoi50", "Gc_h", "tl_2.7b", "tl_2.7f", "tl_M.d"):
            k = collections.Counter(wtl(D[f"u5_tv_{'m1' if DEV else 'c'}_{muc}_u5ktt_{s}"], D[f"u5_tv_{'m1' if DEV else 'c'}_{muc}_{nh}_{s}"], m)
                                    for s in seeds if f"u5_tv_{'m1' if DEV else 'c'}_{muc}_{nh}_{s}" in D and f"u5_tv_{'m1' if DEV else 'c'}_{muc}_u5ktt_{s}" in D)
            c.append(f"{k['W']}/{k['T'] + k['B']}/{k['L']}"); bb = max(bb, k["B"])
        out.append(f"| {TEN[nh]} | " + " | ".join(c) + f" | {bb} |")
    okall = [D[i] for i in ids_m if D[i]["tt"] == "ok"]
    san = sorted(x for r in okall for x in r.get("san", []))
    u = [D[i] for i in ids_m if kv(i)[1] == "u5ktt" and D[i]["tt"] == "ok"]
    duoi_h = all((st.median([gt(D[i], "hd_max_tat_ca") for i in ids_m if kv(i)[1] == nh]) or 0) < H for nh in NH if any(kv(i)[1] == nh for i in ids_m))
    out += ["", f"- Sàn D_i bộ kiểm (bộ phận PROJECT không chạm): p99 = {fm(san[int(0.99 * (len(san) - 1))] if san else None, 4)} µm (n = {len(san)}).",
            f"- U5+KTT: KTT nhận {sum(1 for r in u if r['ktt_nhan'])}/{len(u)}; U5 chọn khác argmin P ở {sum(1 for r in u if r.get('u5_khac_P'))}/{len(u)}.",
            f"- Lệch N2 so với HD bộ dựng (max): {fm(max((r['lech_bo_dung'] for r in okall), default=None), 6)} µm.",
            f"- so_luong khác yêu cầu (S7): {[i for i in ids_m if D[i]['tt'] == 'ok' and yc(i) and any(D[i]['so_luong'].get(k_) != v_ for k_, v_ in yc(i).items())]}.",
            f"- Mọi nhánh HDmax < h: {'CÓ ⇒ chỉ nói «không gây hại»' if duoi_h else 'không'}.", ""]
out += ["Việc hỏng: " + (", ".join(f"{i} ({D[i]['tt']}{', hết giờ' if D[i]['het_gio'] else ''})" for i in KY if D[i]["tt"] != "ok") or "không"),
        "Việc chưa chấm được: " + (", ".join(i for i in KY if D[i]["tt"] == "ok" and D[i].get("tt_k") != "ok") or "không"), ""]
open(RA, "w").write("\n".join(out) + "\n"); print("\n".join(out))
json.dump(D, open(RA.replace(".md", ".json"), "w"), default=lambda x: None if x == INF else float(x), ensure_ascii=False)
