#!/usr/bin/env python3
"""tao_so_lieu.py — sinh so_lieu.tex (macro LaTeX) từ các tệp kết quả trong nguon/.

Bài chỉ dùng macro; không chép số tay vào main.tex / main_vi.tex.
Nguồn và trạng thái rà soát của từng tệp: NGUON.md. Số chỉ có trong tệp rà soát dạng văn bản nằm ở nguon/tay.json,
mỗi mục kèm tệp nguồn và vị trí. Kết quả chưa có ⇒ macro \\TBD{...} (in đậm [TBD]).

Chạy:  python3 tao_so_lieu.py   (từ thư mục so_lieu/)
"""
import json, os, re, hashlib, datetime

D = os.path.dirname(os.path.abspath(__file__)); N = os.path.join(D, "nguon")
J = lambda f: json.load(open(os.path.join(N, f)))
OUT = []


def ten(s):
    """Tên macro: chỉ chữ cái (LaTeX không cho chữ số trong tên macro)."""
    tbl = {"0": "Zero", "1": "One", "2": "Two", "3": "Three", "4": "Four", "5": "Five", "6": "Six", "7": "Seven", "8": "Eight", "9": "Nine"}
    s = "".join(tbl.get(c, c) for c in s)
    s = re.sub(r"[^A-Za-z]", "", s)
    assert s and s[0].isalpha(), s
    return s


def so(x, nd=3):
    """Số có nd chữ số có nghĩa, dấu chấm thập phân (bản tiếng Việt cũng dùng macro này; dấu phẩy đổi trong main_vi bằng \\VNnum nếu cần)."""
    if x is None: return r"\TBD{n/a}"
    if isinstance(x, int) or (float(x).is_integer() and abs(x) >= 10): return f"{int(round(x))}"
    if x == 0: return "0"
    from math import log10, floor
    d = max(0, nd - 1 - int(floor(log10(abs(x)))))
    return f"{x:.{d}f}"


def m(name, val, chu=""):
    OUT.append(f"\\newcommand{{\\{ten(name)}}}{{{val}}}" + (f"  % {chu}" if chu else ""))


def pct(x):
    return f"{abs(100 * x):.0f}"


R = J("ket_qua_rong.json"); MT = R["mo_ta"]; CP = R["cap"]
THUOC = {"HDmax": "HD_max_tat_ca", "HDp": "HD_tu_p95", "HDER": "HD_ER", "Vlost": "Vmat_tat_ca_pt", "VlostER": "Vmat_ER_pt",
         "Gc": "Gc_K_ER", "Lost": "n_mat_kiem", "Half": "con_lai_duoi_50pt", "Ttot": "giay_tong", "Tplace": "giay_dat", "Frag": "vo_manh"}
NH = {"U5b": "U", "P4": "Pfour", "P8": "Peight", "U5b8": "Ueight", "a1": "Aone", "l0": "Lzero"}
VUNG = {"k08": "K", "D5": "Dfive", "D4A2": "DfourAtwo"}
for v, vv in VUNG.items():
    for nh, nn in NH.items():
        k = f"{nh}|{v}"
        if k not in MT: continue
        for t, tk in THUOC.items():
            if tk in MT[k]:
                m(f"r{vv}{nn}{t}", so(MT[k][tk]["tv"]), f"ket_qua_rong.json mo_ta[{k}][{tk}].tv (trung vị, n={MT[k][tk]['n']})")
for kc, tg in [("k08|U5b-P8", "rKUvsPeight"), ("k08|U5b-P4", "rKUvsPfour"), ("D5|U5b-P4", "rDfiveUvsPfour"), ("D5|U5b-P8", "rDfiveUvsPeight"),
               ("k08|U5b-a1", "rKUvsAone"), ("k08|U5b-l0", "rKUvsLzero")]:
    th = CP[kc]["thuoc"]
    for t, tk in THUOC.items():
        if tk not in th: continue
        x = th[tk]
        m(f"{tg}{t}Win", f"{x['thang']}/{x['thua']}", f"cap[{kc}][{tk}] thắng/thua (hoà {x['hoa']})")
        if x.get("tv_hieu_tuong_doi") is not None: m(f"{tg}{t}Rel", pct(x["tv_hieu_tuong_doi"]), f"|tv_hieu_tuong_doi| %, cap[{kc}][{tk}]")
        if x.get("p_wilcoxon") is not None: m(f"{tg}{t}P", so(x["p_wilcoxon"], 2), f"Wilcoxon ký hạng hai phía, cap[{kc}][{tk}]")
m("rNseeds", "10", "seed 6511–6520 mỗi vùng × nhánh (ket_qua_rong.json mo_ta n)")
m("rDfiveHDratio", so(MT["U5b|D5"]["HD_max_tat_ca"]["tv"] / MT["P4|D5"]["HD_max_tat_ca"]["tv"], 2), "tv U5b/P4 HD max ở D5")
m("rDfiveHDfactor", f'{MT["P4|D5"]["HD_max_tat_ca"]["tv"] / MT["U5b|D5"]["HD_max_tat_ca"]["tv"]:.0f}', "tỉ số P4/U5b HD max ở D5 (≈ 9×)")

UC = J("ket_qua_u5c.json")["cap"]
CTL = {"a1l03": "AoneLa", "a1": "AoneLb", "a1l3": "AoneLc", "ah": "Ahalf", "wv": "WV"}
for c, cn in CTL.items():
    x = UC[f"k08|U5b-{c}"]; th = x["thuoc"]
    for t, tk in [("HDmax", "HD_max_tat_ca"), ("Vlost", "Vmat_tat_ca_pt"), ("Half", "con_lai_duoi_50pt"), ("Gc", "Gc_K_ER")]:
        if tk not in th: continue
        m(f"c{cn}{t}", so(th[tk]["tv_B"]), f"ket_qua_u5c.json cap[k08|U5b-{c}][{tk}].tv_B (n={th[tk]['n']})")
        m(f"c{cn}{t}Win", f"{th[tk]['thang']}/{th[tk]['thua']}", f"U5b thắng/thua so với {c}, {tk}")
    m(f"c{cn}LostTot", str(x["mat_tong"][1]), f"tổng bộ phận mất (bộ kiểm) qua mọi seed, {c}; U5b = {x['mat_tong'][0]}")
m("cUHDmax", so(UC["k08|U5b-wv"]["thuoc"]["HD_max_tat_ca"]["tv_A"]), "U5b cùng lô u5c (≡ ket_qua_rong)")
m("cULostTot", str(UC["k08|U5b-wv"]["mat_tong"][0]), "tổng bộ phận mất của U5b, k08 (mat_tong[0])")

X = J("ket_qua_cx.json"); XC = X["cap"]
for kc, tg in [("U5b-CXe", "xUvsCX"), ("CXe-P4", "xCXvsPfour"), ("CXe-P8", "xCXvsPeight")]:
    for t, tk in [("HDmax", "HD_max_tat_ca"), ("HDp", "HD_p95_tat_ca"), ("Vlost", "Vmat_tat_ca_pt"), ("Gc", "Gc_K_ER"), ("Half", "con_lai_duoi_50pt"),
                  ("Ttot", "giay_tong"), ("Tplace", "giay_dat"), ("Tproj", "giay_chieu")]:
        if tk not in XC[kc]: continue
        x = XC[kc][tk]
        m(f"{tg}{t}A", so(x["tv_a"]), f"{kc} {tk} tv trái"); m(f"{tg}{t}B", so(x["tv_b"]), f"{kc} {tk} tv phải")
        m(f"{tg}{t}Win", f"{x['thang']}/{x['thua']}", f"{kc} {tk} thắng/thua (hoà {x['hoa']})")
        if x.get("tv_tuong_doi") is not None: m(f"{tg}{t}Rel", pct(x["tv_tuong_doi"]), f"{kc} {tk} |tv tương đối| %")
        if x.get("p") is not None: m(f"{tg}{t}P", so(x["p"], 2), f"{kc} {tk} p")
q = X["quyet_dinh"]["CXe"]
m("xStepRatio", so(q["ti_dem_CX_U5b_tv"], 3), "ket_qua_cx.json quyet_dinh.CXe.ti_dem_CX_U5b_tv (bước tương đương CX/U5b, trung vị)")
m("xWinsVsPfour", str(q["hieu_luc"]["thang_vs_P4"]), "CXe thắng P4 (Gc_K_ER) trên 10 seed")
cd = J("chon_dev.json")
m("xEqCfg", cd["eq"]["cfg"], "chon_dev.json eq (cấu hình CX chạy ở lượt chính vòng 1)")
for k_, t_ in [("neo", "xEqAnchor"), ("muc", "xEqLevel"), ("m", "xEqM"), ("l0", "xEqLambda")]:
    m(t_, str(cd["eq"][k_]), f"chon_dev.json eq.{k_}")

X2 = J("ket_qua_cx2.json"); C2 = X2["cap"]; VB = {"CXe": "Ce", "CX2e": "Cte", "CX2b": "Ctb"}
TH2 = [("HDmax", "HD_max_tat_ca"), ("HDp", "HD_p95_tat_ca"), ("Vlost", "Vmat_tat_ca_pt"), ("Gc", "Gc_K_ER"), ("Ttot", "giay_tong")]
for v, vn in VB.items():
    for t, tk in TH2:
        x = C2[f"U5b-{v}"][tk]
        m(f"y{vn}{t}", so(x["tv_b"]), f"ket_qua_cx2.json cap[U5b-{v}][{tk}].tv_b")
        m(f"yU{vn}{t}Win", f"{x['thang']}/{x['thua']}", f"U5b vs {v}, {tk}")
        m(f"yU{vn}{t}Rel", pct(x["tv_tuong_doi"]), f"|tv tương đối| U5b vs {v}, {tk}")
        if v == "CXe": m(f"yU{t}", so(x["tv_a"]), f"U5b tv, {tk}")
    for t, tk in [("Gc", "Gc_K_ER"), ("HDmax", "HD_max_tat_ca"), ("Vlost", "Vmat_tat_ca_pt")]:
        x = C2[f"{v}-P4"][tk]; m(f"y{vn}vPfour{t}Win", f"{x['thang']}/{x['thua']}", f"{v} vs P4, {tk}")
        if t == "Gc": m(f"y{vn}vPfourGcRel", ("+" if x["tv_tuong_doi"] > 0 else "-") + pct(x["tv_tuong_doi"]), f"tv tương đối {v} vs P4 (dấu: + là CX tệ hơn)")
    x = C2[f"{v}-P8"]["Gc_K_ER"]; m(f"y{vn}vPeightGcWin", f"{x['thang']}/{x['thua']}", f"{v} vs P8, Gc_K_ER")
for t, tk in TH2:
    m(f"yPfour{t}", so(C2["CX2e-P4"][tk]["tv_b"]), f"P4 tv, {tk}"); m(f"yPeight{t}", so(C2["CX2e-P8"][tk]["tv_b"]), f"P8 tv, {tk}")
cd2 = J("chon_dev2.json")
for g, gn in [("eq", "Cte"), ("best", "Ctb")]:
    for k_, t_ in [("neo", "Anchor"), ("muc", "Level"), ("m", "M"), ("l0", "Lambda")]:
        m(f"y{gn}{t_}", str(cd2[g][k_]), f"chon_dev2.json {g}.{k_}")

T = J("tay.json")
for k, v in T.items():
    if k.startswith("_"): continue
    m(k, v["gia_tri"], f"{v['nguon']}")

TBD = {
        "nfourteen": "n14d5 phan_tich_n14d5.py", "timeExcl": "đo giờ độc quyền rút gọn",
}
for k, v in TBD.items(): m(k, r"\TBD{" + k + "}", "CHƯA CÓ: " + v)

import statistics as _st
TV_MUC = {"D1": "Done", "D5": "Dfive", "real": "Real"}
TV_NH = {"u5ktt": "U", "td": "Place", "nn1x": "Rej", "un1x": "SA", "nn1xm": "RejSoft"}
TV_SEED = {"D1": [7311, 7312, 7313], "D5": list(range(7311, 7317)), "real": list(range(7311, 7317))}
TV_TH = {"HDmax": ("hd_max_tat_ca", 1, 1e-4), "HDp": ("hd_p95_ran", 1, 1e-4), "Vlost": ("pt_mat_ran", 1, 1e-3),
         "Lost": ("mat", 1, 0), "Half": ("duoi50", 1, 0), "SevenB": ("tl_2.7b", 100, 0), "SevenF": ("tl_2.7f", 100, 0),
         "Md": ("tl_M.d", 100, 0), "Eval": ("danh_gia", 1e-6, 0)}
TVJ = J("ket_qua_tv.json")
for mu, mm in TV_MUC.items():
    for nh, nn in TV_NH.items():
        ids = [f"u5_tv_c_{mu}_{nh}_{sd}" for sd in TV_SEED[mu]]
        if nh == "nn1xm" and mu != "D1":
            for th in TV_TH: m(f"tv{mm}{nn}{th}", "--", "không chạy ở mức này")
            continue
        R_ = [TVJ[i] for i in ids]; assert all(r["tt"] == "ok" for r in R_), (mu, nh)
        for th, (k, f, _) in TV_TH.items():
            v = _st.median([r[k] * f for r in R_])
            nd = 2 if th in ("SevenB", "SevenF", "Md") else 3
            m(f"tv{mm}{nn}{th}", so(v, nd) if th not in ("Lost", "Half") else f"{int(round(v))}", f"ket_qua_tv.json trung vị {k}×{f}, {mu} {nh}, n={len(R_)}")
        if nh == "u5ktt":
            m(f"tv{mm}Probe", f"{sum(1 for r in R_ if r.get('ktt_nhan'))}/{len(R_)}", f"KTT nhận (ktt_nhan), {mu}")
    _sl = TVJ[f"u5_tv_c_{mu}_u5ktt_7311"]["so_luong"]
    for k_, n_ in (("luc_lap", "Chl"), ("ti_the", "Mito"), ("peroxisome", "Perox"), ("golgi", "Golgi")):
        m(f"tv{mm}N{n_}", str(_sl.get(k_, 0)), f"ket_qua_tv.json so_luong[{k_}], {mu} (giống mọi seed và nhánh)")
    for nh, nn in TV_NH.items():
        if nh == "u5ktt" or (nh == "nn1xm" and mu != "D1"): continue
        for th in ("HDmax", "HDp", "Vlost"):
            k, f, tol = TV_TH[th]; w = t = l = 0
            for sd in TV_SEED[mu]:
                a = TVJ[f"u5_tv_c_{mu}_u5ktt_{sd}"][k]; b = TVJ[f"u5_tv_c_{mu}_{nh}_{sd}"][k]
                if abs(a - b) < tol: t += 1
                elif a < b: w += 1
                else: l += 1
            m(f"tv{mm}Uvs{nn}{th}Win", f"{w}/{l}", f"thắng/thua U5+KTT theo seed (hoà {t}), {mu} vs {nh}")

NF_NH = {"U5b": "U", "td": "Place", "nn1x": "Rej", "unc8": "SA"}
NF_TH = {"HD_max_dat": "HDmax", "HD_p95_dat": "HDp", "Vmat_dat_pt": "Vlost", "n_bep": "Half", "n_mat": "Lost", "Gc_K": "Gc",
         "luat_dat": "Rules", "phi_eff": "Phi", "ER_ngoai_mang_pt": "ERout", "Vmat_ER_pt": "ERlost"}
NF_CAP = {"U5b-td": "UvsPlace", "U5b-nn1x": "UvsRej", "U5b-unc8": "UvsSA"}
NF = J("N14D5.json"); B_ = NF["bien_the"]["d5"]
for nh, nn in NF_NH.items():
    for th, tk in NF_TH.items():
        v = B_["bang"][nh][th]; m(f"nf{nn}{tk}", so(v["tv"]), f"N14D5.json d5 bang[{nh}][{th}].tv (n={v['n']}; trung vị trên cặp tế bào×seed)")
for cp, cn in NF_CAP.items():
    for th, tk in NF_TH.items():
        if th not in B_["cap"][cp]: continue
        t, h, k = B_["cap"][cp][th]["cap"]; tt, ht, kt = B_["cap"][cp][th]["te_bao"]
        m(f"nf{cn}{tk}Win", f"{t}/{k}", f"N14D5.json cap[{cp}][{th}].cap tốt/hoà/kém = {t}/{h}/{k}")
        m(f"nf{cn}{tk}Tie", f"{h}", "số cặp hoà")
        m(f"nf{cn}{tk}CellWin", f"{tt}/{kt}", f"theo tế bào (trung vị seed) tốt/hoà/kém = {tt}/{ht}/{kt}")
m("nfHong", str(len(B_["hong"])), "N14D5.json d5 hong")

def _n1(fn):
    T = open(os.path.join(N, fn)).read(); out = {}; cap = {}
    for o in ("D1", "D4", "D4A2", "D5"):
        blk = T[T.index(f"\n## {o}\n") + 1:]
        if "\n## " in blk: blk = blk[:blk.index("\n## ")]
        hd = None
        for l in blk.split("\n"):
            if l.startswith("| nhánh | n |"): hd = [c.strip() for c in l.strip("|").split("|")]; continue
            if hd and l.startswith("| ") and not l.startswith("|---"):
                c = [x.strip() for x in l.strip("|").split("|")]
                if len(c) == len(hd) and c[0] in ("U5b", "PLACE", "thử–loại", "SA", "cellPACK", "Packmol"):
                    out[(o, c[0])] = dict(zip(hd, c))
            mm_ = re.match(r"\| (.+?) − (.+?) \| (HD max \(mọi bộ phận, kể ER\)|%V mất \(kể ER\)) \| (\d+) \((\d+)\) \| (\d+)/(\d+)/(\d+) \|", l)
            if mm_: cap[(o, mm_.group(1), mm_.group(2), mm_.group(3))] = tuple(int(mm_.group(i)) for i in (4, 6, 7, 8))
    return out, cap
def _tv(x):
    v = x.split(" [")[0].strip()
    return None if v in ("–", "-", "") else float(v)
N1_NH = {"U5b": "U", "PLACE": "Place", "thử–loại": "Rej", "SA": "SA", "cellPACK": "CP", "Packmol": "PK"}
N1_O = {"D1": "Done", "D4": "Dfour", "D4A2": "DfourAtwo", "D5": "Dfive"}
N1_TH = {"HD max (mọi bộ phận, kể ER)": "HDmax", "HD p95": "HDp", "%V mất (kể ER)": "Vlost", "bộ phận mất mesh": "Lost",
         "không đặt (loại a)": "NotPlaced", "G_c^K": "Gc", "luật sự thật trượt (trừ 4.3b dưới lưới)": "Rules"}
for fn, pre in (("bang_n1.md", "nl"), ("bang_n1_sb.md", "nlsb")):
    B1, C1 = _n1(fn); assert len(B1) == 23, (fn, len(B1))
    for (o, nh), row in B1.items():
        for col, tk in N1_TH.items():
            v = _tv(row[col]); val = "--" if v is None else (f"{int(round(v))}" if tk in ("Lost", "NotPlaced", "Rules") or abs(v) >= 100 else so(v))
            m(f"{pre}{N1_O[o]}{N1_NH[nh]}{tk}", val, f"{fn} ô {o} nhánh {nh} cột {col} (trung vị)")
    for (o, A, Bn, th), (n_, t, k, h) in C1.items():
        tk = "HDmax" if th.startswith("HD") else "Vlost"
        if A == "U5b": w, l, other = t, k, Bn
        elif Bn == "U5b": w, l, other = k, t, A
        else: continue
        m(f"{pre}{N1_O[o]}Uvs{N1_NH[other]}{tk}Win", f"{w}/{l}", f"{fn} cặp {A} − {Bn} {th}: tốt/kém/hoà = {t}/{k}/{h}, n={n_}")

import csv as _csv
_K = list(_csv.DictReader(open(os.path.join(N, "bang_k08_v2.csv"))))
_TMN = {"U5b": "U", "P4": "Pfour", "P8": "Peight", "CXe": "Ce", "CX2e": "Cte", "CX2b": "Ctb"}
_pp = {}
for r in _K: _pp.setdefault(r["nhanh"], {})[r["seed"]] = float(r["pp"])
for nh, nn in _TMN.items():
    m(f"tm{nn}", f"{int(round(_st.median(_pp[nh].values())))}", f"bang_k08_v2.csv trung vị pp (s), {nh}, n={len(_pp[nh])}")
    if nh != "P8":
        rr = [_pp[nh][sd] / _pp["P8"][sd] for sd in _pp[nh]]
        m(f"tm{nn}RatioPeight", f"{_st.median(rr):.2f}", f"trung vị tỉ số pp/pp_P8 theo seed, {nh}")
_TV = list(_csv.DictReader(open(os.path.join(N, "bang_tv_v2.csv"))))
_ppv = [float(r["pp"]) for r in _TV]
m("tmTvMin", f"{min(_ppv)/1000:.1f}", "bang_tv_v2.csv min pp (ks), TV-real 4 nhánh"); m("tmTvMax", f"{max(_ppv)/1000:.1f}", "bang_tv_v2.csv max pp (ks)")

hdr = ["% so_lieu.tex — SINH TỰ ĐỘNG bởi tao_so_lieu.py, đừng sửa tay. " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
       "% md5 nguồn:"]
for f in sorted(os.listdir(N)):
    p = os.path.join(N, f)
    if os.path.isfile(p): hdr.append(f"%   {hashlib.md5(open(p, 'rb').read()).hexdigest()[:8]}  {f}")
hdr.append(r"\providecommand{\TBD}[1]{\textbf{[TBD]}}")
open(os.path.join(D, "so_lieu.tex"), "w").write("\n".join(hdr + OUT) + "\n")
def vi(l):
    mm = re.match(r"(\\newcommand\{\\[A-Za-z]+\}\{)([^}]*)(\}.*)", l)
    if not mm or not re.fullmatch(r"-?\d+\.\d+", mm.group(2)): return l
    return mm.group(1) + mm.group(2).replace(".", "{,}") + mm.group(3)
open(os.path.join(D, "so_lieu_vi.tex"), "w").write("\n".join(hdr + [vi(l) for l in OUT]) + "\n")
print(len(OUT), "macro")
