# -*- coding: utf-8 -*-
"""luat_v12.py — nạp bộ luật SỰ THẬT v1.2 (đóng băng 11/11/2026) và tệp THAM SỐ QUY ƯỚC v1.2 cho bo_kiem_v3.

ĐỘC LẬP với bộ dựng: chỉ đọc các trường id, lop, luat_cung, nhan, phat_bieu, so_trong_van_ban,
khong_the_hien_co_chu_dich của tệp luật (KHÔNG đọc `may`, `may_goc`, `quy_uoc_may`), và tệp tham số.
Mọi ngưỡng sự thật được TÁCH TỪ CHỮ trong so_trong_van_ban (KHOANG/CAN_DUOI/SO) — đổi tệp luật thì đổi ngưỡng;
không tìm thấy chuỗi thì báo lỗi (không im lặng dùng số cũ).
"""
import os, re, json

HERE = os.path.dirname(os.path.abspath(__file__))
THU_MUC = os.path.join(HERE, "luat_v12")
TEP = {"dv": "luat_su_that_dong_vat_v1_2_20261111.json", "tv": "luat_su_that_thuc_vat_v1_2_20261111.json",
       "tham_so": "tham_so_quy_uoc_v1_2_20261111.json"}
TRUONG_DUOC_DOC = ("id", "lop", "luat_cung", "nhan", "phat_bieu", "so_trong_van_ban", "khong_the_hien_co_chu_dich")


def _doc_luat(vat):
    d = json.load(open(os.path.join(THU_MUC, TEP[vat])))
    luat = [{k: r.get(k) for k in TRUONG_DUOC_DOC} for r in d["luat"]]
    return {"phien_ban": d.get("phien_ban"), "nghia_dau_sao": d.get("nghia_dau_sao"), "luat": luat}


LUAT = {v: _doc_luat(v) for v in ("dv", "tv")}
_TS = json.load(open(os.path.join(THU_MUC, TEP["tham_so"])))
THAM_SO = {p["id"]: p["gia_tri"] for p in _TS["tham_so"]}
PHIEN_BAN = {"dv": LUAT["dv"]["phien_ban"], "tv": LUAT["tv"]["phien_ban"], "tham_so": _TS.get("phien_ban")}


def P(pid):
    """giá trị tham số (lỗi nếu thiếu)"""
    if pid not in THAM_SO:
        raise KeyError("thieu tham so %s trong %s" % (pid, TEP["tham_so"]))
    return THAM_SO[pid]


def luat(vat, rid):
    for r in LUAT[vat]["luat"]:
        if r["id"] == rid:
            return r
    raise KeyError("khong co luat %s/%s" % (vat, rid))


_SO = r"(\d+(?:[.,]\d+)?)"


def _f(s):
    return float(s.replace(",", "."))


def _chuoi(vat, rid, khoa):
    for s in luat(vat, rid)["so_trong_van_ban"]:
        if khoa in s:
            return s
    raise KeyError("luat %s/%s: khong thay so_trong_van_ban chua '%s'" % (vat, rid, khoa))


def KHOANG(vat, rid, khoa, chi_so=0):
    """dải 'a–b' thứ chi_so trong chuỗi so_trong_van_ban có chứa `khoa` → (a, b, chuỗi)."""
    s = _chuoi(vat, rid, khoa)
    m = re.findall(_SO + r"\s*[–-]\s*" + _SO, s)
    if len(m) <= chi_so:
        raise ValueError("luat %s/%s: chuoi '%s' khong co dai thu %d" % (vat, rid, s, chi_so))
    return _f(m[chi_so][0]), _f(m[chi_so][1]), s


def CAN_DUOI(vat, rid, khoa):
    """'≥ a' trong chuỗi có `khoa` → (a, chuỗi)."""
    s = _chuoi(vat, rid, khoa)
    m = re.search(r"≥\s*" + _SO, s)
    if not m:
        raise ValueError("luat %s/%s: chuoi '%s' khong co '≥ so'" % (vat, rid, s))
    return _f(m.group(1)), s


def SO(vat, rid, khoa):
    """số nguyên đứng đầu chuỗi (vd '1 màng tế bào', '2 trung tử', '0 trung thể') → (n, chuỗi)."""
    s = _chuoi(vat, rid, khoa)
    m = re.match(r"\s*" + _SO, s)
    if not m:
        raise ValueError("luat %s/%s: chuoi '%s' khong bat dau bang so" % (vat, rid, s))
    return _f(m.group(1)), s


def THAM_SO_CHUOI(vat, rid, khoa):
    """id tham số được nhắc trong chuỗi so_trong_van_ban (vd 'khe nhân–màng ≤ tham số tv.gan_nhan_mang_toi_da')."""
    s = _chuoi(vat, rid, khoa)
    m = re.search(r"((?:chung|dv|tv|kiem)\.[A-Za-z0-9_.]+)", s)
    if not m:
        raise ValueError("luat %s/%s: chuoi '%s' khong nhac tham so" % (vat, rid, s))
    pid = m.group(1).rstrip(".")
    return pid, P(pid), s
