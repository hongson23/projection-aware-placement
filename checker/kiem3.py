# -*- coding: utf-8 -*-
"""
kiem3.py — BỘ KIỂM ĐỘC LẬP bản 3: đánh giá theo bộ luật SỰ THẬT v1.2 (đóng băng 11/11/2026) với HAI phán quyết.

  * su_that : giá trị đo so với DẢI SỰ THẬT tách từ chữ (so_trong_van_ban) của luật v1.2; nhóm B (lời luật nói rõ
              "≤ tham số X") thì dùng tham số X làm lượng hoá vị từ sự thật.
  * tham_so : giá trị đo so với GIÁ TRỊ VẼ/CHỌN trong tham_so_quy_uoc_v1_2 (nhóm A của DoiChieu: số vẽ, dải vẽ, khe quy ước).
  * Mỗi luật quy ước cũ (v4.1 / v1.4) không còn là luật → "kiem_tham_so" (chỉ phán quyết tham_so), ngưỡng đọc từ tệp tham số.
  * Phần lời luật không đo được trên mesh → KHONG_KIEM_DUOC (cả luật) hoặc MOT_PHAN (một vế); KHÔNG BAO GIỜ tính là ĐẠT.

Độc lập: chỉ đọc phat_bieu / so_trong_van_ban / khong_the_hien_co_chu_dich (+ id, lop, luat_cung, nhan) của luật và tệp tham số
(luat_v12.py); không đọc trường `may`, mã hay đặc tả bộ dựng. Hình học: kiem2.py / kiem_thucvat.py (bản v2 đã qua cổng tương đương
GPU = CPU), ngưỡng quy ước đã chuyển sang đọc tệp tham số.

Dùng:  python kiem3.py <thu_muc> [--vat dong_vat|thuc_vat] [--lop KHTN6|Sinh10]
Ghi:   <thu_muc>/kiem3.json, <thu_muc>/kiem3.csv (một dòng mỗi tế bào × luật × phán quyết).
"""
import os, sys, json, csv, time, argparse, collections
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import luat_v12 as LV
import kiem2 as K
import kiem_thucvat as T

EPS = float(LV.P("chung.eps_giao_um3"))
TX = float(LV.P("chung.tiep_xuc_toi_da_um"))
DS = float(LV.P("chung.dung_sai_khoang"))
DAT, TRUOT, KKD, MP, KAD, KCTS = "DAT", "TRUOT", "KHONG_KIEM_DUOC", "MOT_PHAN", "KHONG_AP_DUNG", "KHONG_CO_THAM_SO"


def trong(x, lo, hi, tol=0.0):
    return x is not None and not (isinstance(x, float) and np.isnan(x)) and lo * (1 - tol) - 1e-12 <= x <= hi * (1 + tol) + 1e-12


def ca(bo_phan, gia_tri, *, st=None, ts=None):
    """st = (dat, nguong_chu) cho phán quyết sự thật; ts = (dat, nguong_chu, [id tham số]) cho phán quyết tham số."""
    d = {"bo_phan": bo_phan, "gia_tri": gia_tri}
    if st is not None:
        d["dat_su_that"] = bool(st[0]); d["nguong_su_that"] = st[1]
    if ts is not None:
        d["dat_tham_so"] = bool(ts[0]); d["nguong_tham_so"] = ts[1]; d["tham_so"] = list(ts[2])
    return d


def tong(cas, khoa, khong_do=None, co_do=True):
    """khong_do: mô tả phần lời luật KHÔNG đo được (→ MOT_PHAN nếu phần còn lại đạt; KKD nếu không còn gì đo)."""
    xs = [c for c in cas if ("dat_" + khoa) in c]
    if any(not c["dat_" + khoa] for c in xs):
        return TRUOT
    if not xs or not co_do:
        return KKD if (khoa == "su_that" or khong_do) else KCTS
    return MP if khong_do else DAT


class BoLuat:
    def __init__(self, vat, lop):
        self.vat = vat; self.lop = lop; self.out = []

    def them(self, rid, cas, *, khong_do_st=None, khong_do_ts=None, kkd_st=None, cach_hieu="", nhom="", ghi=""):
        L = LV.luat(self.vat, rid)
        st = KKD if kkd_st else tong(cas, "su_that", khong_do_st)
        ts = tong(cas, "tham_so", khong_do_ts)
        r = {"id": rid, "phat_bieu": L["phat_bieu"], "so_trong_van_ban": L["so_trong_van_ban"],
             "khong_the_hien_co_chu_dich": L["khong_the_hien_co_chu_dich"], "lop": L["lop"], "luat_cung": L["luat_cung"],
             "nhom_doi_chieu": nhom, "cach_hieu": cach_hieu, "ghi": ghi,
             "phan_khong_do_duoc_su_that": kkd_st or khong_do_st, "phan_khong_do_duoc_tham_so": khong_do_ts,
             "trang_thai_su_that": st, "trang_thai_tham_so": ts, "so_ca": len(cas), "ca": cas}
        if self.lop not in (L["lop"] or []):
            r["neu_ap_dung"] = {"su_that": st, "tham_so": ts}
            r["trang_thai_su_that"] = r["trang_thai_tham_so"] = KAD
        self.out.append(r)
        return r


def dien_tich_cham(A, LB, tol):
    d = np.sqrt(LB.kc2_diem(A.luoi.V))
    ok = (d[A.luoi.F] <= tol).all(1)
    V = A.luoi.V; F = A.luoi.F[ok]
    if len(F) == 0:
        return 0.0
    return float(0.5 * np.linalg.norm(np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]]), axis=1).sum())


def khe_mat(A, B):
    """khoảng cách CHÍNH XÁC mặt A – mặt B: U0 = min khoảng cách đỉnh→tam giác hai chiều là cận trên đạt được;
    kc_mat_mat với bán kính R = U0(1+1e-9) cho giá trị chính xác (≤ U0) mà không liệt kê cặp xa."""
    U0 = min(float(np.sqrt(B.luoi.kc2_diem(A.luoi.V)).min()), float(np.sqrt(A.luoi.kc2_diem(B.luoi.V)).min()))
    if U0 <= 0:
        return 0.0
    d, cx = K.kc_mat_mat(A.luoi, B.luoi, U0 * (1 + 1e-9) + 1e-9)
    return float(min(d, U0))


def chac_roi(SA, SB):
    """True ⇔ CHỨNG MINH được hai khối (mặt SA, SB) rời nhau ⇒ thể tích giao = 0 đúng như phép bool sẽ trả:
    (i) mọi đỉnh của lưới ít đỉnh hơn cách mặt kia > cạnh dài nhất của nó ⇒ hai mặt không chạm;
    (ii) một đỉnh mỗi bên nằm ngoài khối kia (số vòng quấn < 0,5, xa mặt) ⇒ không khối nào chứa khối kia.
    Không chứng minh được → False (khi đó tính bool đầy đủ)."""
    X, Y = (SA, SB) if len(SA.V) <= len(SB.V) else (SB, SA)
    U = float(np.sqrt(Y.kc2_diem(X.V)).min())
    if not (U > X.Lmax > 0):
        return False
    return bool(Y.quan(X.V[:1])[0] < 0.5 and X.quan(Y.V[:1])[0] < 0.5)


def do_sau_xuyen(SA, SB, k=8, le=0.05):
    """ĐỘ SÂU XUYÊN lớn nhất (µm) của phần giao hai khối có mặt SA, SB (kiem2.Luoi): lấy mẫu dày (k+1)(k+2)/2 điểm
    barycentric trên MỌI tam giác của mỗi mặt nằm trong hộp chung (+le), giữ điểm nằm trong khối kia (số vòng quấn > 0,5),
    độ sâu = khoảng cách tới mặt kia; max của hai chiều. (Rà soát 12/11: không đo ở đỉnh — đỉnh mù với nêm dây cung.)"""
    lo = np.maximum(SA.lo, SB.lo) - le; hi = np.minimum(SA.hi, SB.hi) + le
    if np.any(hi < lo):
        return 0.0
    u, v = np.meshgrid(np.arange(k + 1), np.arange(k + 1), indexing="ij"); m = (u + v) <= k
    b = np.stack([u[m], v[m]], 1) / k; w = np.stack([1 - b.sum(1), b[:, 0], b[:, 1]], 1)
    best = 0.0
    for X, Y in ((SA, SB), (SB, SA)):
        C = X.V[X.F].mean(1); R = np.linalg.norm(X.V[X.F] - C[:, None], axis=2).max(1)
        sel = np.where(np.all((C + R[:, None] >= lo) & (C - R[:, None] <= hi), 1))[0]
        for i in range(0, len(sel), 20000):
            P = np.einsum("bj,fjd->fbd", w, X.V[X.F[sel[i:i + 20000]]]).reshape(-1, 3)
            kk = Y.quan(P) > 0.5
            if kk.any():
                best = max(best, float(np.sqrt(Y.kc2_diem(P[kk])).max()))
    return best


def nhan_khong_giong0(C):
    """Đăng ký trước v12.1 §9: trả lý do (chuỗi) nếu KHỐI ĐẶC NHÂN (kiem2 khoi_dac_nhan) không phải trường hợp genus 0 chính xác,
    None nếu là genus 0. 'Genus 0' ⇔ khối đặc dựng từ mặt ngoài giống 0 (cach = mat_ngoai_giong0_chinh_xac). Mặt ngoài giống ≥ 1
    (khối đặc = đóng hình thái voxel, hoặc = chính vật liệu có đường hầm) hay dựng lỗi → không genus 0."""
    nhan = C.get("nhan")
    if nhan is None:
        return None
    K.lam_day(nhan)
    dc = nhan.dac_cach or {}
    g = dc.get("giong_mat_ngoai")
    if dc.get("cach") == "mat_ngoai_giong0_chinh_xac" and g == 0:
        return None
    try:
        gd = int(nhan.dac.genus()) if nhan.dac is not None else None
    except Exception:
        gd = None
    return "khoi dac nhan khong genus 0 (cach = %s, giong mat ngoai = %s, giong khoi dac dung = %s)" % (dc.get("cach"), g, gd)


NG_MC0 = "độ sâu xuyên ≤ chung.tiep_xuc_toi_da_um = %g µm (PI duyệt 12/11/2026; khe ≥ 0, được chạm)" % TX


TO_TIEN_DV = {("nhan_con", "nhan"), ("trung_tu", "trung_the"), ("tui_golgi", "golgi"), ("mao", "ti_the"),
              ("lo_nhan", "nhan"), ("mang_nhan", "nhan")}


def _thay_quy_uoc(vat):
    """id luật quy ước cũ → các id tham số thay nó (trường thay_luat_quy_uoc của tệp tham số)."""
    m = collections.defaultdict(list)
    for p in LV._TS["tham_so"]:
        for x in (p.get("thay_luat_quy_uoc") or []):
            v, _, rid = str(x).partition(":")
            if v == vat:
                m[rid.split(" ")[0].strip()].append(p["id"])
    return m


def kiem_tham_so_cu(KQ_cu, vat, bo_qua_ids):
    """Luật quy ước cũ (không còn là luật sự thật) → phép kiểm THAM SỐ (phán quyết tham_so; ca giữ nguyên bản v2,
    ngưỡng đã đọc từ tệp tham số)."""
    thay = _thay_quy_uoc(vat); out = []
    for r in KQ_cu:
        if r["id"] in bo_qua_ids:
            continue
        cas = r.get("ca_tat_ca", r["ca_truot"] + r["ca_dat_mau"])
        out.append({"id": r["id"], "phat_bieu_cu": r["phat_bieu"], "tham_so": thay.get(r["id"], []),
                    "trang_thai_tham_so": r["trang_thai"], "neu_ap_dung": r.get("neu_ap_dung"), "so_ca": len(cas),
                    "ca": [dict(c, dat_tham_so=bool(c.get("dat", True)), nguong_tham_so=c.get("nguong")) for c in cas],
                    "cach_hieu": r.get("cach_hieu"), "ly_do": r.get("ly_do")})
    return out


def _cu(KQ_cu):
    return {r["id"]: r for r in KQ_cu}


def _cas(r):
    return r.get("ca_tat_ca", r["ca_truot"] + r["ca_dat_mau"]) if r else []


def danh_gia_dv(C, lop):
    tb, mang, nhan = C["tb"], C["mang"], C["nhan"]; cu = _cu(C["KQ"]); B = BoLuat("dv", lop); P = LV.P
    theo = tb.theo_loai
    n_m, s_m = LV.SO("dv", "1.1", "màng"); n_n, s_n = LV.SO("dv", "1.1", "nhân")
    cs = [ca("mang", len(theo.get("mang", [])) - 0 if not tb.mat_theo_loai.get("mang") else "MAT", st=(len(theo.get("mang", [])) == n_m and not tb.mat_theo_loai.get("mang"), s_m)),
          ca("nhan", len(theo.get("nhan", [])), st=(len(theo.get("nhan", [])) == n_n and not tb.mat_theo_loai.get("nhan"), s_n))]
    for c in _cas(cu.get("1.1"))[2:]:
        cs.append(ca(c["bo_phan"], c["gia_tri"], st=(c["dat"], c["nguong"])))
    B.them("1.1", cs, cach_hieu="Đếm màng = 1, nhân = 1 (so_trong_van_ban); chất tế bào = khối màng − khối nhân > ε; nhân trong màng.")
    cs = []
    if mang is not None:
        e = mang.kich_thuoc(); lo, s_lo = LV.CAN_DUOI("dv", "1.2a", "đường kính"); pid, hi, s_hi = LV.THAM_SO_CHUOI("dv", "1.2a", "cận trên")
        cs.append(ca(mang.ten, {"kich_thuoc_obb_um": [round(x, 3) for x in e], "the_tich_um3": round(mang.the_tich, 2)},
                     st=(min(e) >= lo * (1 - DS) and mang.the_tich > 0, "khối 3 chiều, mọi cạnh OBB " + s_lo),
                     ts=(max(e) <= hi * (1 + DS), "mọi cạnh OBB ≤ %g µm" % hi, [pid])))
    else:
        cs.append(ca("mang", "khong co", st=(False, "phai co mang")))
    B.them("1.2a", cs, nhom="B", cach_hieu="Đường kính = ba cạnh hộp bao định hướng; sự thật: ≥ 10 µm (lời); cận trên 30 µm là THAM SỐ chọn (phán quyết tham_so).")
    cs = [ca(c["bo_phan"], c["gia_tri"], st=(c["dat"], c["nguong"])) for c in _cas(cu.get("1.3")) if "thong tin" not in str(c.get("nguong"))]
    for c in _cas(cu.get("1.3")):
        if "thong tin" in str(c.get("nguong")) and isinstance(c.get("gia_tri"), (int, float)):
            lt = tuple(P("dv.lech_tam"))
            cs.append(ca(c["bo_phan"], c["gia_tri"], ts=(trong(c["gia_tri"], *lt), "lệch tâm / bán kính trong %s" % (lt,), ["dv.lech_tam"])))
    B.them("1.3", cs, cach_hieu="1 nhân, khối đặc liền, trong tế bào (sự thật); lệch tâm chỉ là tham số vẽ dv.lech_tam.")
    B.them("1.4a", [ca(c["bo_phan"], c["gia_tri"], st=(c["dat"], c["nguong"])) for c in _cas(cu.get("1.4a"))], cach_hieu="Khối nhân − khối màng ≤ ε.")
    cs = []; r0, r1, s_r = LV.KHOANG("dv", "1.5", "tỉ số"); w0, w1, s_w = LV.KHOANG("dv", "1.5", "rộng")
    dn = tuple(P("dv.dai_nhan"))
    for c in _cas(cu.get("1.5")):
        g = c["gia_tri"]
        if isinstance(g, (int, float)):
            cs.append(ca(c["bo_phan"], g, st=(trong(g, r0, r1, DS), s_r)))
        elif isinstance(g, list):
            cs.append(ca(c["bo_phan"], g, st=(all(trong(x, w0, w1, DS) for x in g), s_w + " (mọi cạnh OBB)"),
                         ts=(trong(g[0] / 2, *dn), "nửa trục dài (cạnh OBB lớn / 2) trong %s µm" % (dn,), ["dv.dai_nhan"])))
        else:
            cs.append(ca(c["bo_phan"], g, st=(c["dat"], c["nguong"])))
    B.them("1.5", cs, cach_hieu="(i) tỉ số đường kính cầu tương đương khối đặc nhân / khối tế bào; (ii) mọi cạnh OBB nhân. Tham số: nửa trục dài dv.dai_nhan.")
    nc = theo.get("nhan_con", []); n0, n1, s_c = LV.KHOANG("dv", "4.2a", "nhân con", 0); k0, k1, s_k = LV.KHOANG("dv", "4.2a", "rộng")
    snc = tuple(P("dv.so_nhan_con")); dnc = tuple(P("dv.dai_nhan_con"))
    cs = [ca("nhan_con", len(nc), st=(trong(len(nc), n0, n1) and not tb.mat_theo_loai.get("nhan_con"), s_c),
             ts=(trong(len(nc), *snc), "số vẽ %s" % (snc,), ["dv.so_nhan_con"]))]
    for b in nc:
        e = b.kich_thuoc()
        cs.append(ca(b.ten, [round(x, 3) for x in e], st=(all(trong(x, k0, k1, DS) for x in e), s_k + " (mọi cạnh OBB)"),
                     ts=(trong(e[0] / 2, *dnc), "bán kính (cạnh OBB lớn / 2) trong %s" % (dnc,), ["dv.dai_nhan_con"])))
        r = K.chua(b, nhan, 0.0) if nhan is not None else {"trong": False, "ngoai_um3": None}
        cs.append(ca(b.ten, "ngoai khoi nhan %s µm³" % r["ngoai_um3"], st=(r["trong"], "nằm trong khối đặc nhân")))
    B.them("4.2a", cs, khong_do_st="'không có màng bao' — chi tiết bên trong, không có trên mesh bao", nhom="A",
           cach_hieu="Sự thật: số 1–5, mỗi cái rộng 0,5–5 µm, trong nhân. Tham số: số vẽ dv.so_nhan_con, bán kính vẽ dv.dai_nhan_con.")
    B.them("4.3a", [ca(c["bo_phan"], c["gia_tri"], st=(c["dat"], c["nguong"])) for c in _cas(cu.get("4.3a"))], cach_hieu="Ti thể ⊂ khối màng, không giao khối đặc nhân.")
    t0, t1, s_t = LV.KHOANG("dv", "4.3b", "tiết diện"); dt = tuple(P("dv.dai_ti_the")); tl = tuple(P("dv.ti_le_rong_dai_ti_the"))
    cs = []
    for c in _cas(cu.get("4.3b")):
        g = c["gia_tri"]
        if not isinstance(g, dict):
            cs.append(ca(c["bo_phan"], g, st=(c["dat"], c["nguong"]))); continue
        w, sw, dai = g["rong_day_cuc_bo"], g["sai_so_day"], g["dai"]
        ok = trong(w, t0, t1, DS) or (w < t0 and w + sw >= t0)
        cs.append(ca(c["bo_phan"], {"tiet_dien_um": w, "sai_so": sw, "dai_um": dai},
                     st=(ok, s_t + " (độ dày cục bộ = 2 × bán kính cầu nội tiếp lớn nhất; sai số đo nghiêng về ĐẠT)"),
                     ts=(trong(dai / 2, *dt) and trong(w / dai if dai else float('nan'), *tl),
                         "nửa dài %s µm và rộng/dài %s" % (dt, tl), ["dv.dai_ti_the", "dv.ti_le_rong_dai_ti_the", "kiem.dv.NG_CONG_43B"])))
    B.them("4.3b", cs, nhom="A", cach_hieu="Tiết diện = độ dày cục bộ ∈ 0,2–1 µm (sự thật v1.2). Dài/tỉ lệ rộng–dài là tham số vẽ.")
    er = theo.get("luoi_noi_chat", []); pid, tx, s_tx = LV.THAM_SO_CHUOI("dv", "4.4a", "ER")
    cs = [ca("luoi_noi_chat", t, st=(False, "phai co")) for t in tb.mat_theo_loai.get("luoi_noi_chat", [])]
    if not er:
        cs.append(ca("luoi_noi_chat", "khong co", st=(False, "phai co")))
    for b in er:
        if nhan is None:
            cs.append(ca(b.ten, "khong co nhan", st=(False, s_tx))); continue
        d = float(K.kc_mat_mat(b.luoi, nhan.luoi, 1.0)[0]); S = dien_tich_cham(b, nhan.luoi, tx)
        cs.append(ca(b.ten, {"khe_um": round(d, 6), "dien_tich_cham_um2": round(S, 5)}, st=(d <= tx and S > 0, s_tx)))
    B.them("4.4a", cs, nhom="B", khong_do_st="ER hạt/trơn, ribosome, 'trải khắp chất tế bào' (khong_the_hien_co_chu_dich)",
           cach_hieu="ER chạm mặt ngoài nhân: khe chính xác ≤ chung.tiep_xuc_toi_da_um và diện tích mặt ER có cả 3 đỉnh cách nhân ≤ ngưỡng > 0.")
    go = theo.get("golgi", []); g0, g1, s_g = LV.KHOANG("dv", "4.5a", "túi"); r0_, r1_, s_rg = LV.KHOANG("dv", "4.5a", "rộng")
    stv = tuple(P("dv.so_tui_golgi_ve")); cs = []; kkd = None
    for b in go:
        if b.so_khoi >= 2:
            cs.append(ca(b.ten, b.so_khoi, st=(trong(b.so_khoi, g0, g1), s_g), ts=(trong(b.so_khoi, *stv), "số túi vẽ %s" % (stv,), ["dv.so_tui_golgi_ve"])))
            for k, c in enumerate([c for c in b.comps if c.volume > 0]):
                e = np.sort(c.bounding_box_oriented.primitive.extents)[::-1]
                cs.append(ca("%s#%d" % (b.ten, k), round(float(e[0]), 3), st=(trong(float(e[0]), r0_, r1_, DS), s_rg + " (cạnh OBB lớn)")))
        else:
            kkd = "%s là MỘT khối liền: không đếm/đo được túi" % b.ten
    if not go:
        cs.append(ca("golgi", "khong co", st=(False, "phai co")))
    B.them("4.5a", cs, kkd_st=kkd if (kkd and not cs) else None, khong_do_st=kkd if (kkd and cs) else None,
           cach_hieu="Golgi có túi là khối rời: đếm túi 4–11, rộng mỗi túi 0,7–1,1 µm; một khối liền → KHÔNG KIỂM ĐƯỢC.")
    tt = theo.get("trung_the", []); cs = []
    for b in go:
        if nhan is None or mang is None: break
        dn_, dm_ = khe_mat(b, nhan), khe_mat(b, mang)
        cs.append(ca(b.ten, {"khe_nhan_um": round(dn_, 4), "khe_mang_um": round(dm_, 4)}, st=(dn_ < dm_, LV._chuoi("dv", "4.5b", "khe"))))
        for c in tt:
            cc = np.array(c.tm.center_mass)
            try:
                from scipy.spatial import ConvexHull, Delaunay
                inh = bool(Delaunay(b.luoi.V[ConvexHull(b.luoi.V).vertices]).find_simplex(cc[None])[0] >= 0)
            except Exception:
                inh = False
            gt = khe_mat(b, c); g_ = float(P("dv.gan_tt_golgi"))
            cs.append(ca([b.ten, c.ten], {"tam_trung_the_trong_bao_loi": inh, "khe_golgi_trung_the_um": round(gt, 4)},
                         st=(inh, LV._chuoi("dv", "4.5b", "tâm")), ts=(gt <= g_, "khe Golgi–trung thể ≤ %g µm" % g_, ["dv.gan_tt_golgi"])))
    if not go or not tt:
        cs.append(ca("golgi/trung_the", "thieu", st=(False, "phai co Golgi va trung the")))
    B.them("4.5b", cs, cach_hieu="Vị từ tương đối của lời luật (khe chính xác mặt–mặt). Tham số: ngưỡng tuyệt đối dv.gan_tt_golgi của bộ dựng (cách hiểu của bộ kiểm: khe Golgi–trung thể).")
    cs = []
    for b in tt:
        if nhan is None or mang is None: break
        dn_, dm_ = khe_mat(b, nhan), khe_mat(b, mang); g_ = float(P("dv.gan_tt_nhan"))
        cs.append(ca(b.ten, {"khe_nhan_um": round(dn_, 4), "khe_mang_um": round(dm_, 4)}, st=(dn_ < dm_, LV._chuoi("dv", "4.8b", "khe")),
                     ts=(dn_ <= g_, "khe trung thể–nhân ≤ %g µm" % g_, ["dv.gan_tt_nhan"])))
    if not tt:
        cs.append(ca("trung_the", "khong co", st=(False, "phai co")))
    B.them("4.8b", cs, cach_hieu="Khe chính xác trung thể–nhân < trung thể–màng. Tham số: dv.gan_tt_nhan.")
    ly = theo.get("lysosome", []); l0, l1, s_l = LV.KHOANG("dv", "4.6a", "đường kính"); dl = tuple(P("dv.dai_lysosome"))
    cs = [ca(t, "MAT", st=(False, "phai ton tai")) for t in tb.mat_theo_loai.get("lysosome", [])]
    for b in ly:
        e = b.kich_thuoc()
        cs.append(ca(b.ten, [round(x, 3) for x in e], st=(all(trong(x, l0, l1, DS) for x in e), s_l + " (mọi cạnh OBB)"),
                     ts=(trong(e[0] / 2, *dl), "bán kính vẽ (cạnh OBB lớn / 2) %s" % (dl,), ["dv.dai_lysosome"])))
    B.them("4.6a", cs, nhom="", cach_hieu="Đường kính = cạnh OBB. Tham số: bán kính vẽ dv.dai_lysosome.")
    cs = [ca(c["bo_phan"], c["gia_tri"], st=(c["dat"], c["nguong"])) for c in _cas(cu.get("4.6b"))]
    if ly and nhan is not None:
        dn_ = [khe_mat(b, nhan) for b in ly]; g_ = float(P("dv.gan_lyso_nhan"))
        cs.append(ca("lysosome", round(float(np.median(dn_)), 4), ts=(float(np.median(dn_)) <= g_, "trung vị khe lysosome–nhân ≤ %g µm" % g_, ["dv.gan_lyso_nhan"])))
    B.them("4.6b", cs, cach_hieu="u = d(nhân)/(d(nhân)+d(màng)) tại tâm lysosome; trung vị u < 0,5 và có u ≥ 0,5 (lời luật v1.2). Tham số: dv.gan_lyso_nhan (cách hiểu: trung vị khe tới nhân).")
    n_tt, s_ntt = LV.SO("dv", "4.8a", "trung tử"); a0, a1, s_a = LV.KHOANG("dv", "4.8a", "dài"); b0, b1, s_b = LV.KHOANG("dv", "4.8a", "rộng")
    goc = float(P("dv.goc_vuong_toi_thieu")); dtt = tuple(P("dv.dai_trung_the")); cs = []; kkd = None
    for b in tt:
        e_ = b.kich_thuoc()
        if b.so_khoi == 1:
            kkd = "%s: 1 khối (mặt bao gộp) — không có hai trung tử để đo" % b.ten
            cs.append(ca(b.ten, round(e_[0] / 2, 3), ts=(trong(e_[0] / 2, *dtt), "bán kính vẽ %s" % (dtt,), ["dv.dai_trung_the"])))
        elif b.so_khoi != n_tt:
            cs.append(ca(b.ten, b.so_khoi, st=(False, s_ntt)))
        else:
            ax = []
            for k, c in enumerate([c for c in b.comps if c.volume > 0]):
                e = np.sort(c.bounding_box_oriented.primitive.extents)[::-1]
                cs.append(ca("%s#%d" % (b.ten, k), [round(float(x), 3) for x in e], st=(trong(e[0], a0, a1, DS) and trong(e[1], b0, b1, DS) and trong(e[2], b0, b1, DS), s_a + "; " + s_b)))
                v = c.vertices - c.vertices.mean(0); w, U = np.linalg.eigh(np.cov(v.T)); ax.append(U[:, -1])
            g = float(np.degrees(np.arccos(abs(np.clip(ax[0] @ ax[1], -1, 1)))))
            cs.append(ca(b.ten, round(g, 1), st=(g >= goc, "vuông góc: góc hai trục ≥ dv.goc_vuong_toi_thieu = %g°" % goc)))
    if not tt:
        cs.append(ca("trung_the", "khong co", st=(False, "phai co")))
    stc = [c for c in cs if "dat_su_that" in c]
    B.them("4.8a", cs, kkd_st=kkd if not stc else None, khong_do_st=kkd if stc else None, nhom="B",
           cach_hieu="2 khối trung tử: dài/rộng theo OBB, góc hai trục ≥ tham số (lượng hoá 'vuông góc'); 1 khối → KHÔNG KIỂM ĐƯỢC (chỉ còn tham số bán kính vẽ).")
    B.them("M.a", [ca(c["bo_phan"], c["gia_tri"], st=(c["dat"], c["nguong"])) for c in _cas(cu.get("M.a"))], cach_hieu="Mọi bộ phận − khối màng ≤ ε.")
    r = cu.get("M.d")
    B.them("M.d", [ca(c["bo_phan"], c["gia_tri"], st=(c["dat"], c["nguong"])) for c in _cas(r)],
           kkd_st=(r.get("ly_do") or "nhan dac") if r["trang_thai"] == "KHONG_KIEM_DUOC" else None,
           khong_do_st="khoang 20–50 nm: lời luật ghi 'chỉ ghi, không kiểm trên mô hình'",
           cach_hieu="Chỉ kiểm được khi nhân là vỏ có khoang: có lỗ (giống ≥ 1) và hai lớp theo tia (kiem.dv.MD_*).")
    mao = theo.get("mao", []); ti = theo.get("ti_the", [])
    if mao:
        cs = [ca(b.ten, sum(1 for m_ in mao if np.all(m_.lo >= b.lo - 1e-6) and np.all(m_.hi <= b.hi + 1e-6)), st=None) for b in ti]
        for c_ in cs:
            c_["dat_su_that"] = c_["gia_tri"] >= LV.CAN_DUOI("dv", "M.e", "mào")[0]; c_["nguong_su_that"] = LV._chuoi("dv", "M.e", "mào")
        B.them("M.e", cs, cach_hieu="Đếm bộ phận 'mao' nằm trong hộp từng ti thể.")
    else:
        B.them("M.e", [], kkd_st="mesh ti thể là mặt bao; không có bộ phận mào", cach_hieu="Không có bộ phận mào → KHÔNG KIỂM ĐƯỢC.")
    kdn = nhan_khong_giong0(C)
    bq = [b for b in tb.bp.values() if b.loai != "mang"]; cs = []; n_xet = 0; kkd_nhan = []
    for i in range(len(bq)):
        for j in range(i + 1, len(bq)):
            A_, B_ = bq[i], bq[j]
            if A_.loai == B_.loai or (A_.loai, B_.loai) in TO_TIEN_DV or (B_.loai, A_.loai) in TO_TIEN_DV:
                continue
            if K.hop_cach(A_, B_) > 0:
                continue
            n_xet += 1
            if kdn and "nhan" in (A_.loai, B_.loai):
                kkd_nhan.append([A_.ten, B_.ten])
                cs.append({"bo_phan": [A_.ten, B_.ten], "gia_tri": {"trang_thai": KKD, "ly_do": kdn}})
                continue
            K.khoi_cua(A_); K.khoi_cua(B_)
            if chac_roi(K.mat_chua(A_), K.mat_chua(B_)):
                cs.append(ca([A_.ten, B_.ten], {"do_sau_um": 0.0, "giao_um3": 0.0}, st=(True, NG_MC0))); continue
            g, pp = K.the_tich_giao(A_, B_)
            dsx = do_sau_xuyen(K.mat_chua(A_), K.mat_chua(B_)) if g > 0 else 0.0
            cs.append(ca([A_.ten, B_.ten], {"do_sau_um": round(dsx, 6), "giao_um3": round(g, 8)}, st=(dsx <= TX, NG_MC0)))
    cs.append(ca("*", n_xet, st=(True, "so cap khac loai hop cham da xet (thong tin)")))
    r_ = B.them("M.c₀", cs, khong_do_st=("%d cặp có nhân KHÔNG KIỂM ĐƯỢC: %s" % (len(kkd_nhan), kdn)) if kkd_nhan else None,
                cach_hieu="Mọi cặp bộ phận khác loại có hộp bao chạm nhau (kể cả ER–nhân trước đây miễn khe): độ sâu xuyên của phần giao (bool > 0 thì lấy mẫu dày; nhân dùng khối đặc) ≤ chung.tiep_xuc_toi_da_um. Cặp hộp rời: không thể giao. Khối đặc nhân không genus 0 → mọi cặp có nhân KHÔNG KIỂM ĐƯỢC (đăng ký trước v12.1 §9).",
                ghi="so cap hop cham da xet: %d; cap co nhan khong kiem duoc: %d" % (n_xet, len(kkd_nhan)))
    r_["khong_kiem_duoc_nhan"] = {"so_cap": len(kkd_nhan), "ly_do": kdn, "cap": kkd_nhan}
    C["_mc0_kkd_nhan"] = r_["khong_kiem_duoc_nhan"]
    kcl = dict(P("chung.khe_cung_loai")["dv"]); cs = []
    r = cu.get("M.h")
    for c in _cas(r):
        cs.append(ca(c["bo_phan"], c["gia_tri"], st=(False, c["nguong"])))
    for lo_, ds_ in theo.items():
        if lo_ == "mang": continue
        for i in range(len(ds_)):
            for j in range(i + 1, len(ds_)):
                a_, b_ = ds_[i], ds_[j]; ng = float(kcl.get(lo_, 0.0))
                if ng <= 0 or K.hop_cach(a_, b_) >= ng:
                    continue
                k = C["cap"].get((a_.ten, b_.ten)) or C["cap"].get((b_.ten, a_.ten)) or K.khe(a_, b_, ng)
                cs.append(ca([a_.ten, b_.ten], k["kc_um"], ts=(k["kc_um"] >= ng and k.get("giao_um3", 0) <= EPS, "khe cùng loại ≥ %g µm" % ng, ["chung.khe_cung_loai"])))
    cs.append(ca("*", r.get("so_cap_cung_loai_da_xet") if r else None, st=(True, "so cap cung loai hop cham da xet (thong tin)")))
    B.them("M.h", cs, cach_hieu="Mọi cặp cùng loại có hộp chạm: giao chính xác ≤ ε (chạm được). Tham số: khe cùng loại chung.khe_cung_loai (ti thể 0,3).")
    return B.out


def danh_gia_tv(C, lop):
    tb = C["tb"]; cu = _cu(C["KQ"]); B = BoLuat("tv", lop); P = LV.P; tl = tb.tl
    mang, vach, kb, nhan = C["mang"], C["vach"], C["kb"], C["nhan"]; IN_man, IN_L, V_in = C["IN_man"], C["IN_L"], C["V_in"]
    def lay(rid, khoa="su_that"):
        return [ca(c["bo_phan"], c["gia_tri"], st=(c["dat"], c["nguong"])) if khoa == "su_that" else
                ca(c["bo_phan"], c["gia_tri"], ts=(c["dat"], c["nguong"], [])) for c in _cas(cu.get(rid))]
    cs = []
    for lo_, khoa in (("vach", "thành"), ("mang", "màng")):
        n, s = LV.SO("tv", "2.1", khoa); cs.append(ca(lo_, len(tl.get(lo_, [])), st=(len(tl.get(lo_, [])) == n, s)))
    for lo_, khoa in (("nhan", "nhân"), ("khong_bao", "không bào"), ("luc_lap", "lục lạp")):
        n, s = LV.CAN_DUOI("tv", "2.1", khoa); cs.append(ca(lo_, len(tl.get(lo_, [])), st=(len(tl.get(lo_, [])) >= n, s)))
    cs.append(ca("nhan", len(tl.get("nhan", [])), ts=(len(tl.get("nhan", [])) == int(P("tv.so_nhan")), "số nhân vẽ = %s" % P("tv.so_nhan"), ["tv.so_nhan"])))
    cs += [ca(c["bo_phan"], c["gia_tri"], st=(c["dat"], c["nguong"])) for c in _cas(cu.get("2.1")) if c["nguong"] != ">= 1"]
    B.them("2.1", cs, cach_hieu="Đếm theo so_trong_van_ban (1 thành, 1 màng, ≥ 1 nhân/không bào/lục lạp) + màng trong thành, chất tế bào > 0, nhân/không bào/lục lạp trong màng. Tham số: tv.so_nhan.")
    v0, v1, s_v = LV.KHOANG("tv", "2.2a", "thể tích"); cs = []
    if IN_man is not None:
        cs.append(ca("mang", {"V_trong_mang_pL": round(V_in / 1000.0, 3)}, st=(trong(V_in / 1000.0, v0, v1, DS), s_v)))
    else:
        cs.append(ca("mang", "khong co mat trong mang", st=(False, s_v)))
    for c in _cas(cu.get("2.2a")):
        if "thong tin" in str(c["nguong"]): continue
        cs.append(ca(c["bo_phan"], c["gia_tri"], ts=(c["dat"], c["nguong"] + " (cạnh OBB thành)", ["tv.kich_thuoc_vach", "kiem.tv.TV_TOL_HOP"])))
    B.them("2.2a", cs, nhom="A", cach_hieu="Sự thật: thể tích khối phía TRONG mặt trong của màng ∈ 60–100 pL (PI 11/11). Tham số: hộp ngoài của thành = tv.kich_thuoc_vach ± TV_TOL_HOP.")
    B.them("2.3a", lay("2.3a"), cach_hieu="Màng trong khối thành, không giao vật liệu thành; thành là vỏ kín.")
    n, s = LV.SO("tv", "2.4a", "không bào"); KB = tl.get("khong_bao", [])
    cs = [ca("khong_bao", len(KB), st=(len(KB) == n and not tb.mat_tl.get("khong_bao"), s))] + [ca(b.ten, b.so_khoi, st=(b.so_khoi == 1, "một khối")) for b in KB]
    B.them("2.4a", cs, khong_do_st="'có màng riêng bao quanh' — không có bộ phận màng không bào", cach_hieu="Đúng 1 không bào, 1 khối.")
    B.them("2.4b", lay("2.4b"), cach_hieu="V(không bào) / V(trong màng) ∈ 0,68–0,85 (dải đọc từ lời luật v1.2).")
    pid, g_, s_g = LV.THAM_SO_CHUOI("tv", "2.6b", "khe"); knm = tuple(P("tv.khe_nhan_mang")); cs = []
    for c in _cas(cu.get("2.6b")):
        gv = c["gia_tri"]
        if not isinstance(gv, dict): cs.append(ca(c["bo_phan"], gv, st=(False, c["nguong"]))); continue
        d = gv.get("khe_mang_um")
        cs.append(ca(c["bo_phan"], gv, st=(d is not None and d <= g_ + 1e-9 and gv.get("khong_trong_khong_bao"), s_g + "; không giao không bào"),
                     ts=(d is not None and trong(d, *knm), "khe nhân–màng trong %s (khe quy ước)" % (knm,), ["tv.khe_nhan_mang"])))
    B.them("2.6b", cs, nhom="A", cach_hieu="Sự thật: nhân sát màng (khe ≤ tv.gan_nhan_mang_toi_da, lời luật nêu tham số), không ở trong không bào. Tham số: dải khe quy ước tv.khe_nhan_mang.")
    pid, gl, s_gl = LV.THAM_SO_CHUOI("tv", "2.7b", "khe"); pid2, smr, s_smr = LV.THAM_SO_CHUOI("tv", "2.7b", "mặt"); klm = tuple(P("tv.khe_luc_lap_mang")); cs = []
    for c in _cas(cu.get("2.7b")):
        gv = c["gia_tri"]
        if isinstance(gv, dict) and "so_mat" in gv:
            cs.append(ca(c["bo_phan"], gv, st=(gv["so_mat"] >= smr, s_smr))); continue
        if not isinstance(gv, (int, float)): cs.append(ca(c["bo_phan"], gv, st=(False, c["nguong"]))); continue
        cs.append(ca(c["bo_phan"], gv, st=(gv <= gl + 1e-9, s_gl), ts=(trong(gv, *klm), "khe lục lạp–màng trong %s (khe quy ước)" % (klm,), ["tv.khe_luc_lap_mang"])))
    B.them("2.7b", cs, nhom="A", cach_hieu="Sự thật: mỗi lục lạp khe tới mặt trong màng ≤ tv.gan_luc_lap_mang_toi_da; rải trên ≥ tv.so_mat_rai mặt hộp. Tham số: dải khe quy ước.")
    B.them("2.7c", lay("2.7c"), cach_hieu="Khối lồi, đĩa tròn, hồ sơ thấu kính (kiem.tv.TV_NG_LOI, TV_NEP_SAU, TV_NG_TRON_DIA, TV_THAU_KINH).")
    d0, d1, s_d = LV.KHOANG("tv", "2.7d", "đường kính"); t0, t1, s_t = LV.KHOANG("tv", "2.7d", "dày"); cs = []
    for c in _cas(cu.get("2.7d")):
        e = c["gia_tri"]
        if not isinstance(e, list): cs.append(ca(c["bo_phan"], e, st=(False, c["nguong"]))); continue
        cs.append(ca(c["bo_phan"], e, st=(trong(e[0], d0, d1, DS) and trong(e[1], d0, d1, DS) and trong(e[2], t0, t1, DS), s_d + "; " + s_t),
                     ts=(c["dat"], c["nguong"], ["tv.kich_thuoc_luc_lap_ve"])))
    B.them("2.7d", cs, nhom="A", cach_hieu="Hai cạnh OBB lớn = đường kính, cạnh nhỏ = dày. Sự thật 4–10 × 1,11–2,5; tham số vẽ 5–7 × 2–2,5.")
    B.them("M.d", lay("M.d"), nhom="B", cach_hieu="Ti thể, Golgi, peroxisome cách mặt trong màng hoặc nhân ≤ tv.gan_mang_hoac_nhan (lời luật nêu tham số).")
    cs = lay("M.h")
    nh_ = collections.defaultdict(list)
    for b in tb.bp.values():
        if b.loai not in ("vach", "mang", "ribosome"): nh_[b.loai].append(b)
    cs.append(ca("*", sum(len(T.cap_ung_vien(xs, xs, 1e-9, cung=True)) for xs in nh_.values()), st=(True, "so cap cung loai hop cham da xet (thong tin)")))
    for c in _cas(cu.get("K.1")):
        bp = c["bo_phan"]
        if isinstance(bp, list) and len(bp) == 2 and bp[0] in tb.bp and bp[1] in tb.bp and tb.bp[bp[0]].loai == tb.bp[bp[1]].loai:
            cs.append(ca(bp, c["gia_tri"], ts=(c["dat"], c["nguong"], ["chung.khe_cung_loai", "tv.khe_trong_luc_lap"])))
    B.them("M.h", cs, cach_hieu="Mọi cặp cùng loại: giao ≤ ε (chạm được). Tham số: khe cùng loại (chung.khe_cung_loai).")
    NC = tl.get("nhan_con", []); n_, s_ = LV.CAN_DUOI("tv", "5.1a", "nhân con")
    cs = [ca("nhan_con", len(NC), st=(len(NC) >= n_, s_))]
    for c in _cas(cu.get("5.1a")):
        if c["bo_phan"] == "nhan_con": cs.append(ca(c["bo_phan"], c["gia_tri"], ts=(c["dat"], c["nguong"], ["tv.so_nhan_con"])))
        else: cs.append(ca(c["bo_phan"], c["gia_tri"], st=(c["dat"], c["nguong"])))
    cs += [ca(c["bo_phan"], c["gia_tri"], ts=(c["dat"], c["nguong"], ["tv.duong_kinh_nhan_con"])) for c in _cas(cu.get("5.1b"))]
    B.them("5.1a", cs, nhom="A", khong_do_st="'không có màng bao' không phân biệt được trên mesh",
           cach_hieu="Sự thật: ≥ 1 nhân con, mỗi cái trong nhân. Tham số: số vẽ tv.so_nhan_con, đường kính vẽ tv.duong_kinh_nhan_con.")
    B.tham_5_2 = None
    B.them("5.2a", lay("5.2a"), khong_do_st="khoang 20–50 nm: lời luật ghi 'chỉ ghi, không kiểm trên mô hình'",
           cach_hieu="Theo 400 tia: hai lớp (kiem.tv.TV_MD_TIA), có khoang, có lỗ nhân, màng ngoài chạm ER.")
    n_, s_ = LV.CAN_DUOI("tv", "5.3e", "mào"); cs = []
    for c in _cas(cu.get("5.3e")):
        if isinstance(c["gia_tri"], int) and "mao" in str(c["nguong"]) and str(c["bo_phan"]) in tb.bp and tb.bp[str(c["bo_phan"])].loai == "ti_the":
            cs.append(ca(c["bo_phan"], c["gia_tri"], st=(c["gia_tri"] >= n_, s_), ts=(c["dat"], c["nguong"], ["tv.mao_ve"])))
        else:
            cs.append(ca(c["bo_phan"], c["gia_tri"], ts=(c["dat"], c["nguong"], ["tv.mao_ve"])))
    if not tl.get("ti_the"):
        cs.append(ca("ti_the", 0, st=(True, "khong co ti the (thong tin)")))
    B.them("5.3e", cs, nhom="A", kkd_st=None if any("dat_su_that" in c and c["bo_phan"] != "ti_the" for c in cs) or not tl.get("ti_the") else "khong co bo phan mao",
           cach_hieu="Sự thật: ≥ 1 mào mỗi ti thể. Tham số: số mào/dày vẽ tv.mao_ve.")
    ER = tl.get("luoi_noi_chat", []); pid, tx, s_tx = LV.THAM_SO_CHUOI("tv", "5.4a", "ER"); vung = float(P("tv.vung_noi_er_mang_nhan"))
    mn = C.get("mn"); VN = C.get("VN"); mat = VN.LO if VN is not None and VN.LO is not None else (nhan.luoi if nhan is not None else None)
    cs = [ca(t, "MAT", st=(False, "phai ton tai")) for t in tb.mat_tl.get("luoi_noi_chat", [])]
    if not ER: cs.append(ca("luoi_noi_chat", "khong co", st=(False, "phai co")))
    for e_ in ER:
        S = dien_tich_cham(e_, mat, tx) if mat is not None else 0.0
        cs.append(ca(e_.ten, {"dien_tich_cham_um2": round(S, 4), "mat": "mang_nhan (mat ngoai)" if (VN is not None and VN.LO is not None) else "nhan"},
                     st=(S > 0, s_tx), ts=(S >= vung, "diện tích nối ≥ %g µm²" % vung, ["tv.vung_noi_er_mang_nhan"])))
    B.them("5.4a", cs, nhom="A", khong_do_st="ER hạt/trơn không phân biệt được trên một vùng ER",
           cach_hieu="ER chạm mặt ngoài màng nhân (hoặc mặt nhân nếu không có bộ phận màng nhân): diện tích mặt ER cách ≤ chung.tiep_xuc_toi_da_um > 0. Tham số: ≥ tv.vung_noi_er_mang_nhan.")
    GO = tl.get("golgi", []); n_, s_n = LV.CAN_DUOI("tv", "5.6a", "thể Golgi"); s_k = LV._chuoi("tv", "5.6a", "khe")
    kg = float(P("tv.khe_golgi_golgi")); cs = [ca("golgi", len(GO), st=(len(GO) >= n_, s_n))]
    for a, b in T.cap_ung_vien(GO, GO, kg, cung=True):
        d, g, cx = T.khe(a, b, kg)
        cs.append(ca([a.ten, b.ten], round(d, 4), st=(d > TX and g <= EPS, s_k + " (chạm ⇔ khe ≤ chung.tiep_xuc_toi_da_um)"), ts=(d >= kg - 1e-9, "khe ≥ %g µm" % kg, ["tv.khe_golgi_golgi"])))
    for c in _cas(cu.get("5.6a")):
        if c["bo_phan"] == "golgi" and "khe" not in str(c["nguong"]) and "nhieu" not in str(c["nguong"]):
            cs.append(ca(c["bo_phan"], c["gia_tri"], ts=(c["dat"], c["nguong"] + " (định nghĩa bộ kiểm cho 'rải khắp / không gom')", ["kiem.tv.TV_GOLGI_GAN_NHAN"])))
    B.them("5.6a", cs, nhom="A", khong_do_st="'rải khắp chất tế bào; không gom thành dải' — lời luật không cho số; bộ kiểm chỉ báo ở phán quyết tham số",
           cach_hieu="Sự thật: ≥ 2 thể Golgi, không cặp nào chạm. Tham số: khe ≥ tv.khe_golgi_golgi; định nghĩa 'không gom' của bộ kiểm.")
    p0, p1, s_p = LV.KHOANG("tv", "5.7a", "đường"); cs = []
    for c in _cas(cu.get("5.7a")):
        e = c["gia_tri"]
        if not isinstance(e, list): cs.append(ca(c["bo_phan"], e, st=(False, c["nguong"]))); continue
        cs.append(ca(c["bo_phan"], e, st=(all(trong(x, p0, p1, DS) for x in e), s_p + " (mọi cạnh OBB)"), ts=(c["dat"], c["nguong"], ["tv.duong_kinh_peroxisome_ve", "kiem.tv.TV_NG_CAU"])))
    B.them("5.7a", cs, nhom="A", khong_do_st="'một màng' không kiểm được trên mesh bao", cach_hieu="Sự thật: mọi cạnh OBB 0,1–1 µm ('cầu hoặc bầu dục' → không đòi độ tròn). Tham số: dải vẽ + độ tròn TV_NG_CAU.")
    B.them("5.9", lay("5.9"), cach_hieu="Không có bộ phận loại trung_the.")
    LL = tl.get("luc_lap", []); n_, s_n = LV.SO("tv", "10.1a", "màng trong"); cs = []
    for l_ in LL:
        mt = tb.con_loai(l_, "mang_trong_lap")
        cs.append(ca(l_.ten, len(mt), st=(len(mt) == n_, s_n)))
        for m_ in mt:
            ok, hv = T.trong._goc(m_, l_.man, l_.luoi) if hasattr(T.trong, "_goc") else T.trong(m_, l_.man, l_.luoi)
            cs.append(ca([m_.ten, l_.ten], round(hv, 6), st=(ok, "màng trong nằm trong lục lạp (chứa, t = 0)")))
    for c in _cas(cu.get("10.1a")):
        if isinstance(c["gia_tri"], list):
            cs.append(ca(c["bo_phan"], c["gia_tri"], ts=(c["dat"], c["nguong"], ["tv.khe_mang_trong_lap", "kiem.tv.TV_TOL_KHE_LAP"])))
    B.them("10.1a", cs, nhom="A", khong_do_st="DNA và ribosome chất nền không vẽ", cach_hieu="Sự thật: đúng 1 màng trong mỗi lục lạp, nằm trong. Tham số: khe hai màng vẽ.")
    n_, s_n = LV.CAN_DUOI("tv", "10.2a", "granum"); cs = []
    for c in _cas(cu.get("10.2a")):
        bp = c["bo_phan"]
        if bp in tb.bp and tb.bp[bp].loai == "luc_lap" and isinstance(c["gia_tri"], int):
            cs.append(ca(bp, c["gia_tri"], st=(c["gia_tri"] >= n_, s_n), ts=(c["dat"], c["nguong"], ["tv.so_granum"])))
        else:
            cs.append(ca(bp, c["gia_tri"], st=(c["dat"], c["nguong"] + " (kiem.tv.TV_PHANG/TV_GOC_THANG/TV_LECH_CHONG)")))
    B.them("10.2a", cs, nhom="A", cach_hieu="Sự thật: ≥ 1 granum mỗi lục lạp, mỗi granum là chồng thylakoid dẹt. Tham số: số grana vẽ tv.so_granum.")
    B.them("10.3a", lay("10.3a"), cach_hieu="Đồ thị granum–tấm nối liền một thành phần mỗi lục lạp.")
    B.them("M.a", lay("M.a"), cach_hieu="Mọi bộ phận trên cùng trong khối trong màng; con trong cha.")
    B.them("2.7f", lay("2.7f"), nhom="B", cach_hieu="Trục ngắn PCA lệch pháp tuyến màng ≤ tv.lech_huong_luc_lap (lời luật nêu tham số).")
    e0, e1, s_e = LV.KHOANG("tv", "2.6e", "%", 1)
    cs = [ca(c["bo_phan"], c["gia_tri"], st=(trong(c["gia_tri"], e0, e1, DS) if isinstance(c["gia_tri"], (int, float)) else False, s_e)) for c in _cas(cu.get("2.6e"))]
    B.tham_26f = None
    B.them("2.6e", cs, cach_hieu="V(nhân) / V(trong màng).")
    g0, g1, s_g = LV.KHOANG("tv", "10.2b", "đường")
    B.them("10.2b", [ca(c["bo_phan"], c["gia_tri"], st=(all(trong(x, g0, g1, DS) for x in c["gia_tri"]) if isinstance(c["gia_tri"], list) else False, s_g + " (hai cạnh OBB lớn)")) for c in _cas(cu.get("10.2b"))],
           cach_hieu="Hai cạnh OBB lớn của granum.")
    h0, h1, s_h = LV.KHOANG("tv", "10.2c", "thylakoid")
    B.them("10.2c", [ca(c["bo_phan"], c["gia_tri"], st=(trong(c["gia_tri"], h0, h1) if isinstance(c["gia_tri"], (int, float)) else False, s_h), ts=(c["dat"], c["nguong"], ["tv.so_thylakoid_ve"])) for c in _cas(cu.get("10.2c"))],
           nhom="A", cach_hieu="Số thylakoid con mỗi granum. Tham số: số vẽ tv.so_thylakoid_ve.")
    cs = []
    for c in _cas(cu.get("10.3b")):
        S = c.get("S_cham", c["gia_tri"] if isinstance(c["gia_tri"], (int, float)) else 0.0)
        cs.append(ca(c["bo_phan"], c["gia_tri"], st=(S > 0, LV._chuoi("tv", "10.3b", "granum")), ts=(c["dat"], c["nguong"], ["tv.vung_noi_tam_granum"])))
    B.them("10.3b", cs, nhom="A", cach_hieu="Mỗi granum có tấm nối chạm (≤ chung.tiep_xuc_toi_da_um), diện tích > 0. Tham số: ≥ tv.vung_noi_tam_granum.")
    bs = [b for b in tb.bp.values() if b.loai not in ("vach", "mang")]; cs = []; n_xet = 0
    def to_tien(a, b):
        if a.ten in tb.to_tien(b.ten) or b.ten in tb.to_tien(a.ten):
            return True
        for x, y in ((a, b), (b, a)):
            if x.loai == "mang_trong_lap" and y.loai in ("granum", "thylakoid", "tam_noi") and set(tb.to_tien(x.ten)) & set(tb.to_tien(y.ten)):
                return True
        return False
    for a, b in T.cap_ung_vien(bs, bs, 1e-9, cung=True):
        if a.loai == b.loai or to_tien(a, b):
            continue
        n_xet += 1
        if chac_roi(a.luoi, b.luoi):
            cs.append(ca([a.ten, b.ten], {"do_sau_um": 0.0, "giao_um3": 0.0}, st=(True, NG_MC0))); continue
        g = T.giao(a, b)
        dsx = do_sau_xuyen(a.luoi, b.luoi) if g > 0 else 0.0
        cs.append(ca([a.ten, b.ten], {"do_sau_um": round(dsx, 6), "giao_um3": round(g, 8)}, st=(dsx <= TX, NG_MC0)))
    cs.append(ca("*", n_xet, st=(True, "so cap khac loai hop cham da xet (thong tin)")))
    B.them("M.c₀", cs, cach_hieu="Mọi cặp bộ phận khác loại có hộp chạm (kể cả ER–nhân, ER–màng nhân, ribosome–ER, tấm nối–granum… trước đây miễn khe): độ sâu xuyên (lấy mẫu dày khi bool > 0) ≤ chung.tiep_xuc_toi_da_um. Trừ thành, màng và cặp tổ tiên–con cháu ('*').")
    return B.out


def _phang(KQ3, KTS, cu=None):
    """Dạng phẳng cho bộ cấy lỗi: id '<luật>@su_that' / '<luật>@tham_so'; mỗi ca có 'dat', 'nguong'."""
    out = []
    for r in KQ3:
        for khoa in ("su_that", "tham_so"):
            xs = [dict(bo_phan=c["bo_phan"], gia_tri=c["gia_tri"], dat=c["dat_" + khoa], nguong=c["nguong_" + khoa]) for c in r["ca"] if ("dat_" + khoa) in c]
            out.append({"id": r["id"] + "@" + khoa, "trang_thai": r["trang_thai_" + khoa], "luat_cung": r["luat_cung"],
                        "loai_truot": (cu or {}).get(r["id"], {}).get("loai_truot"),
                        "ca_truot": [c for c in xs if not c["dat"]], "ca_dat_mau": [c for c in xs if c["dat"]][:30]})
    for r in KTS:
        xs = [dict(bo_phan=c.get("bo_phan"), gia_tri=c.get("gia_tri"), dat=c["dat_tham_so"], nguong=c["nguong_tham_so"]) for c in r["ca"]]
        out.append({"id": r["id"] + "@tham_so", "trang_thai": r["trang_thai_tham_so"], "luat_cung": False, "loai_truot": (cu or {}).get(r["id"], {}).get("loai_truot"),
                    "ca_truot": [c for c in xs if not c["dat"]], "ca_dat_mau": [c for c in xs if c["dat"]][:30]})
    return out


def kiem3(thu_muc, vat="dong_vat", lop=None, im_lang=False):
    t0 = time.time()
    if vat == "dong_vat":
        cu = K.kiem(thu_muc, lop, os.path.join(HERE, "luat_v41.json"), im_lang=True)
        C = K._NGU_CANH; lop = cu["lop"]; v = "dv"
        KQ3 = danh_gia_dv(C, lop)
    else:
        cu = T.kiem_tv(thu_muc, lop, os.path.join(HERE, "luat_thucvat_v14_van_ban.json"), im_lang=True)
        C = T._NGU_CANH; lop = cu["lop"]; v = "tv"
        KQ3 = danh_gia_tv(C, lop)
    ids_st = {l["id"] for l in LV.LUAT[v]["luat"]}
    assert {r["id"] for r in KQ3} == ids_st, ({r["id"] for r in KQ3} ^ ids_st)
    bo = set(ids_st) | {"4.5d", "M.h"}
    KTS = kiem_tham_so_cu(C["KQ"], v, bo)
    for r in KTS:
        if r["neu_ap_dung"] is None and r["trang_thai_tham_so"] == "KHONG_AP_DUNG":
            pass
    dem = {k: dict(collections.Counter(r["trang_thai_" + k] for r in KQ3)) for k in ("su_that", "tham_so")}
    out = {"te_bao": os.path.basename(os.path.normpath(thu_muc)), "vat": vat, "lop": lop, "phien_ban_luat": LV.PHIEN_BAN,
           "so_luat_su_that": len(KQ3), "dem_su_that": dem["su_that"], "dem_tham_so": dem["tham_so"],
           "truot_su_that": [r["id"] for r in KQ3 if r["trang_thai_su_that"] == TRUOT],
           "truot_tham_so": [r["id"] for r in KQ3 if r["trang_thai_tham_so"] == TRUOT],
           "kiem_tham_so_cu_truot": [r["id"] for r in KTS if r["trang_thai_tham_so"] == TRUOT],
           "loi_luoi": cu.get("loi_luoi"), "bo_phan_mat": cu.get("bo_phan_mat"),
           "thoi_gian_s": round(time.time() - t0, 1), "luat_su_that": KQ3, "kiem_tham_so": KTS, "luat": _phang(KQ3, KTS, _cu(C["KQ"]))}
    for k_ in ("chong_lan_cung_loai", "khoi_dac_nhan", "nghi_sai_loai"):
        if k_ in cu: out[k_] = cu[k_]
    out["M.c0_khong_kiem_duoc_nhan"] = C.get("_mc0_kkd_nhan") or {"so_cap": 0, "ly_do": None, "cap": []}
    json.dump(out, open(os.path.join(thu_muc, "kiem3.json"), "w"), ensure_ascii=False, indent=1, default=K._json_default)
    with open(os.path.join(thu_muc, "kiem3.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["te_bao", "vat", "lop", "luat", "loai_luat", "phan_quyet", "trang_thai", "so_ca", "so_ca_truot", "bo_phan_truot", "tham_so"])
        for r in KQ3:
            for k in ("su_that", "tham_so"):
                tr = [c for c in r["ca"] if c.get("dat_" + k) is False]
                ts_ = sorted({p for c in r["ca"] for p in c.get("tham_so", [])}) if k == "tham_so" else []
                w.writerow([out["te_bao"], vat, lop, r["id"], "su_that_v1.2", k, r["trang_thai_" + k], sum(1 for c in r["ca"] if ("dat_" + k) in c), len(tr),
                            ";".join(json.dumps(c["bo_phan"], ensure_ascii=False) for c in tr[:50]), ";".join(ts_)])
        kk = out["M.c0_khong_kiem_duoc_nhan"]
        w.writerow([out["te_bao"], vat, lop, "M.c₀", "dem_kkd_nhan", "su_that", KKD if kk["so_cap"] else "", kk["so_cap"], 0,
                    ";".join(json.dumps(c, ensure_ascii=False) for c in kk["cap"][:50]), kk["ly_do"] or ""])
        for r in KTS:
            tr = [c for c in r["ca"] if not c["dat_tham_so"]]
            w.writerow([out["te_bao"], vat, lop, r["id"], "quy_uoc_cu", "tham_so", r["trang_thai_tham_so"], r["so_ca"], len(tr),
                        ";".join(json.dumps(c.get("bo_phan"), ensure_ascii=False) for c in tr[:50]), ";".join(r["tham_so"])])
    if not im_lang:
        print("== %s (%s, lop %s) %.1f s  su_that %s  tham_so %s" % (out["te_bao"], vat, lop, out["thoi_gian_s"], dem["su_that"], dem["tham_so"]))
        print("   truot su_that:", out["truot_su_that"], " truot tham_so:", out["truot_tham_so"], " kiem tham so cu truot:", out["kiem_tham_so_cu_truot"])
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("thu_muc"); ap.add_argument("--vat", default="dong_vat", choices=["dong_vat", "thuc_vat"]); ap.add_argument("--lop", default=None)
    a = ap.parse_args()
    kiem3(a.thu_muc, a.vat, a.lop)
