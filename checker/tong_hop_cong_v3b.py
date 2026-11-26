# -*- coding: utf-8 -*-
"""tong_hop_cong_v3b.py — bảng CỔNG ĐỘT BIẾN CUỐI của bo_kiem_v3 (15/11/2026).

Gộp kết quả cổng v3 (hn3_*/cay_loi*_ket_qua.json, chạy 13/11) với lượt chạy lại các đột biến bị ảnh hưởng bởi bản sửa 15/11
(hn4_*/cay_loi*_ket_qua.json) và lượt chạy lại B1 (hn5_*: hai đột biến 4.5b trên gốc g45 mới, 15/11 chiều): mọi bản ghi hn3 có cùng
(tế bào, tên loại) với một bản ghi hn4/hn5 (tên cũ nếu đã đổi tên, trường 'ten_cu') bị THAY; hn5 thắng hn4 nếu trùng.
Ghi TOM_TAT_CONG_V3_CUOI.md, cong_v3_cuoi.json.

Thước của bảng (theo kiểu bản ghi):
  * dot_bien: 'phát hiện' (luật đích có ca trượt MỚI nhắc đúng bộ phận, bất kỳ phán quyết); 'phát hiện mức sự thật' trên các đột biến
    có mục @su_that;
  * st_dat_ts_truot ('đúng sự thật / sai tham số'): đúng kỳ vọng ⇔ có ca trượt mới ở @tham_so và KHÔNG ở @su_that;
  * đối chứng (trung_tinh, sat_nguong_dat): báo động giả;
  * không hợp lệ (gốc đã trượt phần sự thật của luật đích, hoặc tham số với st_dat_ts_truot); lỗi cấy.
"""
import os, sys, json, glob, time, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)


def doc(mau):
    kq = []
    for f in sorted(glob.glob(os.path.join(HERE, mau))):
        for x in json.load(open(f))["ket_qua"]:
            x = dict(x); x["_nguon"] = os.path.relpath(f, HERE); kq.append(x)
    return kq


def ket_luan(x):
    if "khong_hop_le" in x: return "KHONG_HOP_LE"
    if "loi_cay" in x: return "LOI_CAY"
    if x["kieu"] == "dot_bien": return "PHAT_HIEN" if x.get("phat_hien") else "BO_SOT"
    if x["kieu"] == "st_dat_ts_truot": return "DUNG" if x.get("dung_ky_vong") else "SAI"
    return "BAO_GIA" if x.get("bao_dong_gia") else "SACH"


DOI_TV = os.path.join(HERE, "sua_1611", "doi_nhan_tv_1611.json")


def ap_doi_nhan_tv(kq):
    """16/11 (PI): đổi nhãn kỳ vọng 10 đột biến (2 Golgi nhóm C giữ dot_bien sau rà soát chéo) thực vật dot_bien -> st_dat_ts_truot trên bản ghi CÓ SẴN (không chạy lại).
    Bản ghi dot_bien đã có phat_hien_theo_phan_quyet tính bằng đúng V3.theo_pq như nhánh st_dat_ts_truot của cay_loi_tv.py, và
    luat_dich giữ cả @su_that lẫn @tham_so (không mục nào bị bỏ vì trượt ở gốc) ⇒ dung_ky_vong = tham_so and not su_that là chính xác;
    tính hợp lệ theo kiểu mới (V3.loc_hop_le với 'st_dat_ts_truot') và khớp bản ghi đã kiểm lại trong sua_1611/doi_nhan_tv_1611.py."""
    D = {(e["te_bao"], e["loai"]): e for e in json.load(open(DOI_TV))["doi_nhan"]}
    out = []; n = 0
    for x in kq:
        e = D.get((x["te_bao"], x["loai"]))
        if e and x["kieu"] == "dot_bien" and "khong_hop_le" not in x and "loi_cay" not in x:
            assert e["khong_hop_le_neu_st_dat_ts_truot"] is None and e["khop_ban_ghi_cu"] and e["luat_dich"] == x["luat_dich"], e["loai"]
            x = dict(x); pq = x["phat_hien_theo_phan_quyet"]
            x["phat_hien_khi_con_la_dot_bien"] = x.pop("phat_hien", None)
            x.update(kieu="st_dat_ts_truot", kieu_cu="dot_bien", nhom_doi_nhan=e["nhom"], ly_do_doi_nhan=e["ly_do"],
                     dung_ky_vong=bool(pq["tham_so"] and not pq["su_that"]))
            n += 1
        out.append(x)
    assert n == len(D) == 10, n
    return out


def gia_tri(x):
    return x.get("do_duoc") or x.get("chi_tiet") or x.get("gia_tri_luat_dich") or \
        ({k: x["phat_hien_theo_phan_quyet"].get(k) for k in ("vd_su_that", "vd_tham_so")} if x.get("phat_hien_theo_phan_quyet") else "")


o = []; P = o.append
P("# Cổng đột biến bo_kiem_v3 — bảng CUỐI (%s)" % time.strftime("%Y-%m-%d %H:%M"))
P("")
P("Sinh tự động bởi tong_hop_cong_v3b.py từ hn3_* (cổng 13/11) + hn4_* (chạy lại 15/11 các đột biến bị ảnh hưởng bởi bản sửa theo rà soát 14/11) + hn5_* (chạy lại 15/11 chiều hai đột biến 4.5b trên gốc g45 mới — mục chặn B1 của rà soát 15/11).")
P("16/11 (PI): trong 12 đột biến thực vật 'chỉ bắt ở tham số', 10 đã ĐỔI NHÃN dot_bien → st_dat_ts_truot trên kết quả CÓ SẴN (không chạy lại); 2 đột biến Golgi (nhóm C) GIỮ dot_bien vì vi phạm LỜI sự thật 5.6a mà vị từ máy khai báo khong_hinh_hoc. Bảng lý do: mục 'Đổi nhãn kỳ vọng thực vật 16/11', sua_1611/doi_nhan_tv_1611.json.")
P("Rà soát chéo độc lập lần đổi nhãn (16/11, sua_1611/ra_soat_doi_nhan/PHAN_QUYET.txt): NOT READY với bản đổi 12 — đồng ý 10/12, bác 2 ca Golgi nhóm C (CHẶN 1) và số ca tác dụng phụ bị cắt [:6] (CHẶN 2). Bản này sửa theo cả hai mục chặn (phương án (a)); vòng 1b cùng rà soát (16/11, mục 'Vong 1b' của PHAN_QUYET.txt): READY — CHẶN 1, 2 đã đóng, số tính lại độc lập khớp: 71/71, 29/31, 15/15, 0/20, không hợp lệ 5; động vật không đổi; md5 mục A 16/16.")
P("Trạng thái: rà soát chéo vòng 1 (15/11 sáng) + vòng 2 (15/11 chiều, mục B1 và B2; chứng cứ sua_1511/ra_soat_vong2/PHAN_QUYET.txt): B1 READY; B2 văn bản sửa theo vòng 2 và 2b (bản này đã sửa 2 chỗ vòng 2b chỉ ra). Phán quyết và phản biện: DONG_BANG_BO_KIEM_V3.txt, sua_1511/ra_soat_1511.txt.")
P("Bản ghi hn3 cùng (tế bào, loại) với bản ghi hn4/hn5 bị thay (hn5 thắng hn4). Mã: kiem3.py bản 15/11 (thêm luật đăng ký trước §9: khối đặc nhân không genus 0 → M.c₀ cặp có nhân KHÔNG KIỂM ĐƯỢC); cay_loi.py, v3_cay_loi_chung.py bản 15/11; cay_loi_tv.py bản 16/11 (so với 15/11 chỉ đổi kiểu kỳ vọng 10 đột biến, bộ cấy không đổi).")
P("")
tong = {}
for vat, tep in (("dong_vat", "cay_loi_ket_qua.json"), ("thuc_vat", "cay_loi_tv_ket_qua.json")):
    cu = doc("hn3_*/" + tep); m5 = doc("hn5_*/" + tep); k5 = {(x["te_bao"], x.get("ten_cu") or x["loai"]) for x in m5}
    moi = [x for x in doc("hn4_*/" + tep) if (x["te_bao"], x.get("ten_cu") or x["loai"]) not in k5] + m5
    thay = {(x["te_bao"], x.get("ten_cu") or x["loai"]) for x in moi}
    bi_thay = [x for x in cu if (x["te_bao"], x["loai"]) in thay]
    kq = [x for x in cu if (x["te_bao"], x["loai"]) not in thay] + moi
    if vat == "thuc_vat":
        kq = ap_doi_nhan_tv(kq)
    db = [x for x in kq if x["kieu"] == "dot_bien" and "phat_hien" in x]
    st_db = [x for x in db if any(d.endswith("@su_that") for d in x.get("luat_dich", []))]
    sv = [x for x in kq if x["kieu"] == "st_dat_ts_truot" and "dung_ky_vong" in x]
    dc = [x for x in kq if x["kieu"] in ("trung_tinh", "sat_nguong_dat") and "bao_dong_gia" in x]
    khl = [x for x in kq if "khong_hop_le" in x]; lc = [x for x in kq if "loi_cay" in x]
    t = {"so_ban_ghi": len(kq), "so_ban_ghi_hn3_bi_thay": len(bi_thay), "so_ban_ghi_hn4": len(moi),
         "dot_bien": len(db), "phat_hien": sum(bool(x["phat_hien"]) for x in db),
         "dot_bien_muc_su_that": len(st_db), "phat_hien_muc_su_that": sum(bool(x.get("phat_hien_theo_phan_quyet", {}).get("su_that")) for x in st_db),
         "st_dat_ts_truot": len(sv), "st_dat_ts_truot_dung": sum(bool(x["dung_ky_vong"]) for x in sv),
         "doi_chung": len(dc), "doi_chung_bao_gia": sum(bool(x["bao_dong_gia"]) for x in dc),
         "khong_hop_le": len(khl), "loi_cay": len(lc)}
    tong[vat] = t
    if vat == "thuc_vat":
        kq_tv = kq
    P("## %s" % vat)
    P("| mục | số |")
    P("|---|---|")
    P("| đột biến: phát hiện (bất kỳ phán quyết) | %d / %d |" % (t["phat_hien"], t["dot_bien"]))
    P("| đột biến có mục @su_that: phát hiện ở mức SỰ THẬT | %d / %d |" % (t["phat_hien_muc_su_that"], t["dot_bien_muc_su_that"]))
    P("| đúng sự thật / sai tham số: đúng kỳ vọng | %d / %d |" % (t["st_dat_ts_truot_dung"], t["st_dat_ts_truot"]))
    P("| đối chứng (trung tính + sát ngưỡng đạt + sạch): báo động giả | %d / %d |" % (t["doi_chung_bao_gia"], t["doi_chung"]))
    P("| không hợp lệ (loại khỏi mẫu số) | %d |" % t["khong_hop_le"])
    P("| lỗi cấy | %d |" % t["loi_cay"])
    P("| bản ghi: tổng / hn3 bị thay / hn4+hn5 mới | %d / %d / %d |" % (len(kq), len(bi_thay), len(moi)))
    kkd_dc = [x for x in dc if not x["bao_dong_gia"] and (x.get("M.c0_kkd_nhan") or 0) > 0]
    t["doi_chung_sach_nho_KKD_nhan"] = [(x["te_bao"], x["loai"], x.get("M.c0_kkd_nhan")) for x in kkd_dc]
    if kkd_dc:
        P("| đối chứng: báo động giả NẾU coi M.c₀ 'không kiểm được' (nhân không genus 0) là vi phạm (quy ước độ nhạy §9) | %d / %d |" % (t["doi_chung_bao_gia"] + len(kkd_dc), t["doi_chung"]))
    P("")
    if kkd_dc:
        P("Lưu ý: %d đối chứng sạch CHỈ vì các cặp M.c₀ có nhân thành KHÔNG KIỂM ĐƯỢC (khối đặc nhân không genus 0; dc_nhan_khoet_rong: bản 13/11 báo giả ER–nhân do khối đặc đóng hình thái voxel; dc_nhan_khoan_ham: giống 1, sạch cả hai bản nhưng nay M.c₀ cặp có nhân cũng không kiểm): %s. Theo quy ước độ nhạy (KKD = vi phạm) đây là báo động giả." % (len(kkd_dc), t["doi_chung_sach_nho_KKD_nhan"]))
        P("")
    for ten, xs in (("Bỏ sót (đột biến)", [x for x in db if not x["phat_hien"]]),
                    ("Đột biến có mục sự thật nhưng chỉ bắt ở tham số", [x for x in st_db if x["phat_hien"] and not x.get("phat_hien_theo_phan_quyet", {}).get("su_that")]),
                    ("Đúng sự thật / sai tham số SAI kỳ vọng", [x for x in sv if not x["dung_ky_vong"]]),
                    ("Báo động giả (đối chứng)", [x for x in dc if x["bao_dong_gia"]]),
                    ("Không hợp lệ", khl), ("Lỗi cấy", lc)):
        P("%s: %d" % (ten, len(xs)))
        for x in xs:
            P("- %s %s [%s] %s đích %s — %s" % (x["te_bao"], x["loai"], x["kieu"], x["bo_phan"], x.get("luat_dich"),
                                               (x.get("khong_hop_le") or x.get("loi_cay") or json.dumps(x.get("luat_bao_gia") or x.get("phat_hien_theo_phan_quyet"), ensure_ascii=False, default=str))[:260]))
        P("")
    P("### Bản ghi chạy lại 15/11 (hn4; hn5 = B1 4.5b) và kết quả 13/11 tương ứng")
    P("| tế bào | loại (tên cũ) | nguồn | kiểu 13/11 → 15/11 | kết quả 13/11 → 15/11 | bộ phận 15/11 | luật trượt mới 15/11 | giá trị đo / mô tả |")
    P("|---|---|---|---|---|---|---|---|")
    for x in moi:
        ten_cu = x.get("ten_cu") or x["loai"]
        c_ = [y for y in bi_thay if y["te_bao"] == x["te_bao"] and y["loai"] == ten_cu]
        c0 = c_[0] if c_ else None
        P("| %s | %s%s | %s | %s → %s | %s → %s | %s | %s | %s |" % (x["te_bao"], x["loai"], (" (%s)" % ten_cu) if ten_cu != x["loai"] else "", x["_nguon"].split("/")[0],
                                                       c0["kieu"] if c0 else "—", x["kieu"], ket_luan(c0) if c0 else "—", ket_luan(x), x["bo_phan"],
                                                       ", ".join(x.get("luat_truot_moi") or []),
                                                       json.dumps(gia_tri(x) or x.get("mo_ta") or x.get("khong_hop_le"), ensure_ascii=False, default=str)[:230].replace("|", "/")))
    P("")
    P("### M.c₀ (độ sâu xuyên), 4.5b, 4.6b — mọi bản ghi")
    P("| tế bào | đột biến | kiểu | bộ phận | kết quả | giá trị đo |")
    P("|---|---|---|---|---|---|")
    for x in kq:
        if not any(r in (x.get("luat_dich") or []) for r in ("M.c₀@su_that", "4.5b@su_that", "4.6b@su_that")) and "M.c0" not in x["loai"]:
            continue
        P("| %s | %s | %s | %s | %s | %s |" % (x["te_bao"], x["loai"], x["kieu"], x["bo_phan"], ket_luan(x), json.dumps(gia_tri(x), ensure_ascii=False, default=str)[:200].replace("|", "/")))
    P("")
    if vat == "dong_vat":
        kk = [(x["te_bao"], x["loai"], x.get("M.c0_kkd_nhan")) for x in moi if x.get("M.c0_kkd_nhan")]
        P("### M.c₀ KHÔNG KIỂM ĐƯỢC do khối đặc nhân không genus 0 (bản ghi hn4): %s" % (kk or "không có"))
        P("")
    tong[vat]["ban_ghi_hn4"] = [{k: x.get(k) for k in ("te_bao", "loai", "ten_cu", "kieu", "bo_phan", "luat_dich", "phat_hien", "phat_hien_theo_phan_quyet",
                                                      "dung_ky_vong", "bao_dong_gia", "khong_hop_le", "loi_cay", "luat_truot_moi", "mo_ta", "do_duoc", "chi_tiet",
                                                      "gia_tri_luat_dich", "M.c0_kkd_nhan", "thoi_gian_s")} for x in moi]
P("## Đổi nhãn kỳ vọng thực vật 16/11 (PI duyệt hướng; cùng cách động vật 15/11)")
P("- Trước đổi nhãn (bảng 15/11): phát hiện 81/81; mức sự thật 29/41; st_dat_ts_truot 5/5; đối chứng 0/20; không hợp lệ 5. Sau: xem bảng thuc_vat trên.")
P("- Cách làm: chỉ đổi KIỂU kỳ vọng (cay_loi_tv.py: 10 lời gọi A(...) thêm 'st_dat_ts_truot', bản cũ cay_loi_tv.py.bak_20261116_truoc_doi_nhan_tv) và áp lên bản ghi hn3 có sẵn (ap_doi_nhan_tv). KHÔNG chạy lại, KHÔNG đổi mã chấm mục A. Kiểm lại từ kq__*.json + TV1__goc/kiem3.json bằng V3.map_dich/loc_hop_le/theo_pq (sua_1611/doi_nhan_tv_1611.py): 10/10 hợp lệ theo kiểu mới (ở gốc bộ phận đích không có ca trượt ở @su_that lẫn @tham_so của luật đích; @tham_so gốc = DAT; @su_that gốc DAT hoặc MOT_PHAN do vế khong_hinh_hoc), 10/10 khớp bản ghi cũ, 10/10 DUNG. Gốc hn3_TV1_p0 và p1 tương đương theo từng luật.")
P("- Nhóm A = giá trị sau cấy nằm trong dải số của vị từ sự thật; B = đại lượng bị đổi không có trong vị từ sự thật. Nhóm C (KHÔNG đổi nhãn) = vi phạm vế lời 'rải khắp / không gom cạnh nhân' của 5.6a mà chính luật khai báo khong_hinh_hoc và tách sang tham số kiểm kiem.tv.TV_GOLGI_GAN_NHAN.")
P("| đột biến | nhóm | luật | sự thật / tham số sau cấy | lý do (theo văn bản luật + tham số) | sự thật MỚI ngoài luật đích (tác dụng phụ) |")
P("|---|---|---|---|---|---|")
_J = json.load(open(DOI_TV)); _D = _J["doi_nhan"] + _J["giu_dot_bien"]
for e in _D:
    P("| %s%s | %s | %s | %s / %s | %s | %s |" % (e["loai"], "" if e["kieu_moi"] == "st_dat_ts_truot" else " (GIỮ dot_bien)", e["nhom"], e["luat"], e["trang_thai_su_that_sau_cay"], e["trang_thai_tham_so_sau_cay"], e["ly_do"].replace("|", "/"),
      "; ".join("%s: %d ca, sâu tới %.4g µm, vd %s" % (k, e["so_ca_tac_dung_phu"][k], e["do_sau_max_tac_dung_phu_um"][k], v[:2]) for k, v in e["su_that_moi_ngoai_luat_dich_tac_dung_phu"].items()) or "không"))
P("")
_C = _J["giu_dot_bien"]; _t = tong["thuc_vat"]
P("- Nhóm C GIỮ dot_bien (rà soát chéo 16/11, CHẶN 1): %s vi phạm lời sự thật 5.6a ('rải khắp chất tế bào; không gom thành một dải cạnh nhân'); vị từ sự thật hình học chỉ có '≥ 2 thể' + 'khe > 0' nên @su_that trả MOT_PHAN, không bắt được — đây là KHOẢNG TRỐNG PHỦ của vị từ 5.6a theo thiết kế luật, không phải 'tế bào đúng sự thật'. Hai ca này là 2 bỏ sót mức sự thật trong bảng (29 / %d); cả hai bắt được ở @tham_so (2/2). Nếu PI/GV nâng vế này thành vị từ sự thật thì đây là giới hạn của bộ kiểm v3, cần ghi trong bài." % ([e["loai"] for e in _C], _t["dot_bien_muc_su_that"]))
P("- Cách đếm khác (độ nhạy): (b) tách 2 ca nhóm C thành dòng riêng 'vi phạm lời sự thật, vế khong_hinh_hoc → chỉ bắt ở @tham_so: 2/2', loại khỏi mẫu số: mức sự thật %d / %d, st_dat_ts_truot %d / %d; (c) đổi nhãn cả 12 (bản đầu 16/11, bị rà soát bác): 29 / 29, st_dat_ts_truot %d / %d — chỉ là phương án độ nhạy, không dùng." % (
  _t["phat_hien_muc_su_that"], _t["dot_bien_muc_su_that"] - len(_C), _t["st_dat_ts_truot_dung"], _t["st_dat_ts_truot"], _t["st_dat_ts_truot_dung"] + len(_C), _t["st_dat_ts_truot"] + len(_C)))
P("- Tác dụng phụ: 3/10 đột biến đã đổi nhãn tạo thêm ca trượt M.c₀@su_that ở cặp khác luật đích (lục lạp phóng to lấn không bào 0,21 / 0,34 µm; tấm nối hẹp xuyên 0,0112 µm ở CHÍNH cặp tấm nối–granum đang thử, sát sàn phát hiện 0,0100–0,0105). Thước st_dat_ts_truot chỉ xét luật đích (như động vật: mang_phong_1_6 có 1.5/4.6b; v12_luc_lap_duong_kinh_8.0 có M.c₀) nên kết quả không đổi; các ô này KHÔNG phải tế bào 'sạch sự thật toàn phần' — nếu trích làm ví dụ 'đúng sự thật / sai tham số' trong bài phải cấy lại bản không chạm (lượt sau). Hai ca nhóm C cũng có M.c₀ phụ (17 ca lấn không bào ≤ 0,15 µm; 18 ca lấn ER ≤ 0,35 µm).")
_phu = collections.OrderedDict()
for x in kq_tv:
    if x["kieu"] == "dot_bien" and "phat_hien" in x:
        for d in x.get("luat_dich", []):
            if d.endswith("@su_that"):
                a = _phu.setdefault(d, [0, 0]); a[1] += 1; a[0] += int(bool(x.get("phat_hien_theo_phan_quyet", {}).get("su_that")))
_ids = sorted({q["id"] for q in json.load(open(os.path.join(HERE, "luat_v12", "luat_su_that_thuc_vat_v1_2_20261111.json")))["luat"]})
_khong = [i for i in _ids if i + "@su_that" not in _phu]
P("- Phủ mức sự thật thực vật sau đổi nhãn (đột biến dot_bien có mục @su_that, bắt ở sự thật / tổng): %s. Luật sự thật thực vật KHÔNG có đột biến dot_bien nào ở mức sự thật: %s." % (
  ", ".join("%s %d/%d" % (k.split("@")[0], v[0], v[1]) for k, v in sorted(_phu.items())), ", ".join(_khong)))
P("")
P("## B1 (rà soát 15/11): 4.5b trên gốc riêng g45 mới (hn5)")
P("- Gốc g45 13/11 trượt vế 1 (khe Golgi–nhân 2,1 > khe Golgi–màng) nên v12_golgi_sat_mang_4.5b KHÔNG HỢP LỆ ở cả hai tế bào. Gốc mới (_g45_goc, cay_loi.py): xoay cứng Golgi quanh tâm nhân tới chỗ khe nhân < khe màng ngoài vỏ ER, dời ti thể/lysosome vướng (giữ số lượng), trung thể trong Golgi lệch về phía nhân. Vế 1: dời Golgi + trung thể theo bán kính tới khi khe nhân − khe màng = 0,3 µm; vế 2: dời trung thể ngang ra khỏi bao lồi.")
P("- Tác dụng phụ của gốc g45 (rà soát vòng 2): trung thể nằm TRỌN trong khối Golgi (thể tích giao / thể tích trung thể = 1,000) → M.c₀ Golgi–trung thể trượt sẵn ở gốc (độ sâu 0,77 µm D1 / 0,83 µm D4A2; không thành ca mới ở đột biến); M.c và 4.8b tham số trượt ở gốc; D4A2 dời cứng 17 bào quan (13 ti thể 0,83–6,09 µm, 4 lysosome 0,40–0,79 µm; số lượng giữ nguyên), kéo theo M.b tham số ở ti_the_132, ti_the_97. Bốn đột biến: luật trượt MỚI so với g45 chỉ có 4.5b@su_that.")
P("- Phạm vi: Golgi ở cả hai tế bào là MỘT khối gần lồi tuyệt đối (V / V_bao lồi = 0,9989) nên tâm trung thể trong bao lồi ⇔ trong khối Golgi: vế 2 của 4.5b và M.c₀ KHÔNG thể cùng đạt với Golgi dạng này (vấn đề lượng hoá luật ôm quanh trung thể — cần PI biết), và cổng không có ca tâm trong bao lồi nhưng ngoài khối (cần Golgi lõm/nhiều túi). 4.5b không có đối chứng sát ngưỡng (biên các đột biến 0,22–0,86 µm).")
for x in doc("hn5_*/cay_loi_ket_qua.json"):
    P("- %s %s (tên cũ %s): %s; gốc: %s; đột biến: %s; luật trượt mới: %s; đo: %s" % (x["te_bao"], x["loai"], x.get("ten_cu"), ket_luan(x), (x.get("goc_rieng") or {}).get("mo_ta"),
      x.get("mo_ta") or x.get("khong_hop_le") or x.get("loi_cay"), x.get("luat_truot_moi"), json.dumps(x.get("do_duoc"), ensure_ascii=False, default=str)))
P("")
P("## Giới hạn độ nhạy đã biết (M.c₀)")
P("- Độ sâu xuyên đo bằng lấy mẫu barycentric (k = 8) trên tam giác: cận DƯỚI của độ sâu thật; sai số lấy mẫu ≤ 0,73 nm trên mọi ca đo (rà soát vòng 2, 15/11: dựng lại đột biến từ tế bào gốc bằng chính hàm cấy; độ sâu thật = max(khối giao manifold3d lấy mẫu k = 24 trên ∂I; tự lấy mẫu k = 96 quanh điểm sâu nhất), hai cách lệch nhau ≤ 0,11 nm; chứng cứ sua_1511/ra_soat_vong2/q6_*.json, q6_tom.out). D4A2 lồi ER 0,02 → k=8 0,0155 (0,015478; = JSON cổng) / k=32 0,0162 / chỉ đỉnh 0,0136 / thật 0,0162; D4A2 lồi ER 0,05 → 0,0377 (k=8 = k=32 = thật); D1 lồi ER 0,02 → 0,0199 (= thật); D1 lysosome 0,02 → 0,0162 (= thật); D4A2 lysosome 0,02 → 0,0197 (thật 0,0199).")
P("- Phần hụt so với danh nghĩa ở D4A2 lồi ER (0,02 → 0,0155: 22,6 %; 0,05 → 0,0377: 24,6 %) chủ yếu do BỘ CẤY: ở D4A2 _er_loi chỉ dời 2 đỉnh, theo hướng tâm khối nhân lệch pháp tuyến mặt nhân 52°, nên độ lấn thật chỉ 0,0162 / 0,0377 (với 0,02: bộ cấy 19,0 % + lấy mẫu 3,6 điểm %; với 0,05: toàn bộ do bộ cấy). Độ sâu thật lớn hơn dl·cos (0,0123 / 0,0307) nên không phải phép chiếu đơn giản: hướng dời xiên cộng hình học cục bộ của lưới ER thô. D1 lysosome 0,02 → 0,0162: hụt 19 % cũng của bộ cấy (cl_khe), k=8 đúng bằng thật. Tỉ số k=8 / thật 0,955–1,000 trên mọi ca có độ sâu > 0,005 µm.")
P("- Ngưỡng 0,01 µm → sàn phát hiện ≈ 0,0100–0,0105 µm độ sâu thật trên hai lưới này (quét lồi theo pháp tuyến ở D4A2: thật 0,01038 → k=8 0,01034, vượt; thật 0,00989 < 0,01 → k=8 0,00985, không vượt — đúng, không phải số đo sàn; D1: thật 0,0105 → k=8 0,0105, vượt; tỉ số k=8/thật xấu nhất 0,955 ⇒ sàn ≤ 0,01047); xấu nhất ở chỗ lồi 2 đỉnh trên lưới ER thô; phụ thuộc cỡ tam giác, KHÔNG tổng quát hoá cho lưới khác.")
P("- D4A2 gốc (chưa cấy) đã có xuyên ER–nhân ≈ 0,005 µm (thật 0,00497; k=8 0,00457) ở vùng tiếp xúc. Đối chứng sat_nguong_dat lồi ER 0,005 ở D4A2 thực chất đo chỗ tiếp xúc gốc (0,00457), chỗ lồi chỉ ≈ 0,0034 ⇒ đối chứng này không thử được gì ở D4A2.")
P("- Luật đăng ký trước §9: tế bào có khối đặc nhân không genus 0 → mọi cặp M.c₀ có nhân KHÔNG KIỂM ĐƯỢC (đếm trong kiem3.json 'M.c0_khong_kiem_duoc_nhan' và dòng CSV 'dem_kkd_nhan'); phân tích độ nhạy coi các ca này là vi phạm.")
P("")
P("## Tương đương GPU–CPU, bộ cấy v2 (chạy 11/11, so lại 13/11) — không đổi")
f = os.path.join(HERE, "..", "bo_kiem_v2_gpu", "so_sanh_dot_bien.json")
if os.path.exists(f):
    for x in json.load(open(f)):
        if "bo" in x:
            P("- %s: %d/%d tệp giống hệt, lệch số max %g; khác: %d tệp" % (x["bo"], x["giong_het"], x["so_te_bao_dot_bien"], x["lech_so_max"], len(x["khac"])))
    P("- Mọi khác biệt thực vật nằm ở THỨ TỰ danh sách ca trượt K.1; so theo tập thì 0 tệp khác.")
json.dump(tong, open(os.path.join(HERE, "cong_v3_cuoi.json"), "w"), ensure_ascii=False, indent=1, default=str)
open(os.path.join(HERE, "TOM_TAT_CONG_V3_CUOI.md"), "w").write("\n".join(o) + "\n")
print(json.dumps({v: {k: t[k] for k in t if k != "ban_ghi_hn4"} for v, t in tong.items()}, ensure_ascii=False))
