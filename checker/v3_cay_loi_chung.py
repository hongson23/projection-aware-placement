# -*- coding: utf-8 -*-
"""v3_cay_loi_chung.py — phần dùng chung cho bộ cấy lỗi bản 3 (cay_loi.py, cay_loi_tv.py trong bo_kiem_v3).

* Ánh xạ luật đích cũ (v4.1 / v1.4) → mục của kiem3: luật sự thật v1.2 → '<id>@su_that' và '<id>@tham_so';
  luật quy ước cũ → '<id>@tham_so'; 'chong_lan_cung_loai' (động vật) giữ nguyên + 'M.h@su_that'.
* Phát hiện tách theo phán quyết: phat_hien_theo_phan_quyet = {su_that, tham_so}.
* Kiểu mới 'st_dat_ts_truot' (đúng sự thật, sai tham số): ĐÚNG kỳ vọng ⇔ có ca trượt MỚI nhắc bộ phận ở mục @tham_so
  và KHÔNG có ca trượt mới nhắc bộ phận ở mục @su_that của luật đích.
"""
import luat_v12 as LV

ID_ST = {v: {l["id"] for l in LV.LUAT[v]["luat"]} for v in ("dv", "tv")}


def map_dich(ds, vat, co_san):
    out = []
    for d in ds:
        if "@" in d or d in ("chong_lan_cung_loai", "nghi_sai_loai"):
            out.append(d)
            if d == "chong_lan_cung_loai":
                out.append("M.h@su_that")
            continue
        if d in ID_ST[vat]:
            out += [d + "@su_that", d + "@tham_so"]
        else:
            out.append(d + "@tham_so")
    seen = []
    for d in out:
        if (d in co_san or d in ("chong_lan_cung_loai", "nghi_sai_loai")) and d not in seen:
            seen.append(d)
    return seen


def loc_hop_le(dich, kieu, base_luat, truot_bp, toan_te_bao=False):
    """Bộ lọc tính hợp lệ (15/11/2026, sau rà soát cổng v3 14/11).
    dich: các mục đích đã ánh xạ; base_luat: {id: mục luật phẳng của tế bào GỐC}; truot_bp(rid) → True nếu ở gốc mục rid có
    ca trượt nhắc đúng bộ phận đích (hoặc luật không áp dụng); toan_te_bao: phép cấy biến đổi CẢ tế bào (vd kéo dài mọi bộ phận).
    Trả (khong_hop_le: str | None, dich_con_lai, muc_bo).
      * Mục @su_that: gốc đã trượt PHẦN SỰ THẬT (ở đúng bộ phận; hoặc ở bất kỳ ca nào khi phép cấy là toàn tế bào) → KHÔNG HỢP LỆ.
      * Mục @tham_so đã trượt ở gốc: kiểu 'st_dat_ts_truot' → KHÔNG HỢP LỆ (không thử được 'tham số chuyển sang trượt');
        kiểu khác → bỏ riêng mục đó (như bản 13/11), không hợp lệ chỉ khi không còn mục nào.
      * Mục khác (chong_lan_cung_loai, nghi_sai_loai, M.h cũ…): đã trượt ở gốc → bỏ riêng; hết mục → không hợp lệ."""
    st = [d for d in dich if d.endswith("@su_that") and (truot_bp(d) or (toan_te_bao and d in base_luat and base_luat[d]["trang_thai"] == "TRUOT"))]
    if st:
        return "goc da truot phan SU THAT cua luat dich%s: %s" % (" (phep cay toan te bao)" if toan_te_bao else "", st), dich, st
    khac = [d for d in dich if not d.endswith("@su_that") and truot_bp(d)]
    if khac and kieu == "st_dat_ts_truot" and any(d.endswith("@tham_so") for d in khac):
        return "st_dat_ts_truot nhung goc da truot phan THAM SO: %s" % khac, dich, khac
    con = [d for d in dich if d not in khac]
    if khac and not con:
        return "bo phan da truot luat dich o goc: %s" % khac, dich, khac
    return None, con, khac


def theo_pq(dich, moi_theo_luat, nhac_fn):
    """moi_theo_luat: {id: [ca mới]}; trả {su_that: bool, tham_so: bool, ca: [...]}."""
    st = [c for d in dich if d.endswith("@su_that") for c in moi_theo_luat.get(d, []) if nhac_fn(c)]
    ts = [c for d in dich if d.endswith("@tham_so") for c in moi_theo_luat.get(d, []) if nhac_fn(c)]
    return {"su_that": bool(st), "tham_so": bool(ts), "vd_su_that": st[:2], "vd_tham_so": ts[:2]}


def tong_ket_v12(ket_qua):
    db = [x for x in ket_qua if x.get("kieu") == "dot_bien" and "phat_hien" in x]
    sv = [x for x in ket_qua if x.get("kieu") == "st_dat_ts_truot" and "dung_ky_vong" in x]
    def dem(f, loc=lambda x: True):
        xs = [x for x in db if loc(x)]
        return "%d/%d" % (sum(1 for x in xs if f(x)), len(xs))
    co = lambda k: (lambda x: any(d.endswith("@" + k) for d in x.get("luat_dich", [])))
    return {"dot_bien": len(db),
            "phat_hien_bat_ky": dem(lambda x: x["phat_hien"]),
            "phat_hien_su_that (tren dot bien co muc dich @su_that)": dem(lambda x: x.get("phat_hien_theo_phan_quyet", {}).get("su_that"), co("su_that")),
            "phat_hien_tham_so (tren dot bien co muc dich @tham_so)": dem(lambda x: x.get("phat_hien_theo_phan_quyet", {}).get("tham_so"), co("tham_so")),
            "bo_sot_co_muc_su_that (xem lai: dai v1.2 co the da hop le)": [(x["te_bao"], x["loai"], x["bo_phan"], x["luat_dich"]) for x in db if not x["phat_hien"] and co("su_that")(x)],
            "st_dat_ts_truot": {"so": len(sv), "dung_ky_vong": sum(1 for x in sv if x["dung_ky_vong"]),
                                "ca": [(x["te_bao"], x["loai"], x["bo_phan"], x["dung_ky_vong"]) for x in sv]},
            "M.c0": [(x["te_bao"], x["loai"], x["bo_phan"], x.get("phat_hien"), x.get("kieu")) for x in ket_qua if "M.c₀@su_that" in (x.get("luat_dich") or [])]}
