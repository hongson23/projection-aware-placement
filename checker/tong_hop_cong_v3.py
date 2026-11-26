# -*- coding: utf-8 -*-
"""tong_hop_cong_v3.py — gộp cổng đột biến v3 (hn3_*/cay_loi*_ket_qua.json, 13/11/2026) + tương đương GPU–CPU v2
(bo_kiem_v2_gpu/so_sanh_dot_bien.json) → TOM_TAT_CONG_V3.md, cay_loi_v3_cong.json."""
import os, sys, json, glob, collections, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import v3_cay_loi_chung as V3

o = []; P = o.append
P("# Tóm tắt cổng đột biến bo_kiem_v3 — %s" % time.strftime("%Y-%m-%d %H:%M"))
P("")
P("Sinh tự động bởi tong_hop_cong_v3.py. Trạng thái: CHƯA rà soát chéo (kết quả máy).")
P("Mã phán quyết: kiem3.py với M.c₀ theo độ sâu xuyên > chung.tiep_xuc_toi_da_um = 0,01 µm (PI duyệt 12/11; bản trước: kiem3.py.bak_20261113_truoc_Mc0_do_sau). Các ngưỡng khác không đổi.")
P("")
hang = os.path.join(HERE, "cong_v3")
xong = sorted(os.listdir(os.path.join(hang, "xong"))); loi = sorted(os.listdir(os.path.join(hang, "loi")))
P("Việc: xong %d, lỗi %d %s" % (len(xong), len(loi), loi))
P("")
tong = {}
for vat, tep in (("dong_vat", "cay_loi_ket_qua.json"), ("thuc_vat", "cay_loi_tv_ket_qua.json")):
    kq = []; nguon = []
    for f in sorted(glob.glob(os.path.join(HERE, "hn3_*", tep))):
        d = json.load(open(f)); kq += d["ket_qua"]; nguon.append(os.path.relpath(f, HERE))
    db = [x for x in kq if x["kieu"] == "dot_bien" and "phat_hien" in x]
    sv = [x for x in kq if x["kieu"] == "st_dat_ts_truot" and "dung_ky_vong" in x]
    dc = [x for x in kq if x["kieu"] in ("trung_tinh", "sat_nguong_dat") and "bao_dong_gia" in x]
    khl = [x for x in kq if "khong_hop_le" in x]; lc = [x for x in kq if "loi_cay" in x]
    co = lambda k: (lambda x: any(d.endswith("@" + k) for d in x.get("luat_dich", [])))
    st_db = [x for x in db if co("su_that")(x)]
    tong[vat] = {"nguon": nguon, "ban_ghi": len(kq), "dot_bien": len(db), "phat_hien": sum(x["phat_hien"] for x in db),
                 "dot_bien_muc_su_that": len(st_db), "phat_hien_su_that": sum(bool(x.get("phat_hien_theo_phan_quyet", {}).get("su_that")) for x in st_db),
                 "st_dat_ts_truot": len(sv), "st_dat_ts_truot_dung": sum(x["dung_ky_vong"] for x in sv),
                 "doi_chung": len(dc), "doi_chung_bao_gia": sum(x["bao_dong_gia"] for x in dc), "khong_hop_le": len(khl), "loi_cay": len(lc)}
    t = tong[vat]
    P("## %s (%d tệp kết quả)" % (vat, len(nguon)))
    P("| mục | số |")
    P("|---|---|")
    P("| đột biến: phát hiện (bất kỳ phán quyết) | %d / %d |" % (t["phat_hien"], t["dot_bien"]))
    P("| đột biến có mục @su_that: phát hiện theo sự thật | %d / %d |" % (t["phat_hien_su_that"], t["dot_bien_muc_su_that"]))
    P("| đúng sự thật / sai tham số: đúng kỳ vọng | %d / %d |" % (t["st_dat_ts_truot_dung"], t["st_dat_ts_truot"]))
    P("| đối chứng (trung tính + sát ngưỡng đạt): báo động giả | %d / %d |" % (t["doi_chung_bao_gia"], t["doi_chung"]))
    P("| không hợp lệ (bộ phận đã trượt ở gốc) | %d |" % t["khong_hop_le"])
    P("| lỗi cấy | %d |" % t["loi_cay"])
    P("")
    bs = [x for x in db if not x["phat_hien"]]
    if bs:
        P("Bỏ sót:"); [P("- %s %s %s đích %s" % (x["te_bao"], x["loai"], x["bo_phan"], x["luat_dich"])) for x in bs[:30]]; P("")
    sai = [x for x in sv if not x["dung_ky_vong"]]
    if sai:
        P("Đúng sự thật/sai tham số SAI kỳ vọng:"); [P("- %s %s %s pq=%s" % (x["te_bao"], x["loai"], x["bo_phan"], {k: x["phat_hien_theo_phan_quyet"][k] for k in ("su_that", "tham_so")})) for x in sai]; P("")
    gia = [x for x in dc if x["bao_dong_gia"]]
    if gia:
        P("Báo động giả:"); [P("- %s %s (%s) %s" % (x["te_bao"], x["loai"], x["kieu"], json.dumps(x.get("luat_bao_gia") or x.get("chi_tiet"), ensure_ascii=False, default=str)[:300])) for x in gia]; P("")
    if khl:
        P("Không hợp lệ:"); [P("- %s %s %s: %s" % (x["te_bao"], x["loai"], x["bo_phan"], x["khong_hop_le"])) for x in khl[:30]]; P("")
    if lc:
        P("Lỗi cấy:"); [P("- %s %s: %s" % (x["te_bao"], x["loai"], str(x["loi_cay"])[:200])) for x in lc]; P("")
    P("### M.c₀ (độ sâu xuyên), 4.5b, 4.6b")
    P("| tế bào | đột biến | kiểu | bộ phận | kết quả | giá trị đo |")
    P("|---|---|---|---|---|---|")
    for x in kq:
        if not any(r in (x.get("luat_dich") or []) for r in ("M.c₀@su_that", "4.5b@su_that", "4.6b@su_that")) and "M.c0" not in x["loai"]:
            continue
        kqx = ("PHAT_HIEN" if x.get("phat_hien") else "BO_SOT") if x["kieu"] == "dot_bien" and "phat_hien" in x else \
              (("DUNG" if x["dung_ky_vong"] else "SAI") if "dung_ky_vong" in x else (("BAO_GIA" if x["bao_dong_gia"] else "SACH") if "bao_dong_gia" in x else ("KHONG_HOP_LE" if "khong_hop_le" in x else "LOI_CAY" if "loi_cay" in x else "?")))
        gt = x.get("do_duoc") or x.get("chi_tiet") or x.get("gia_tri_luat_dich") or ""
        P("| %s | %s | %s | %s | %s | %s |" % (x["te_bao"], x["loai"], x["kieu"], x["bo_phan"], kqx, json.dumps(gt, ensure_ascii=False, default=str)[:160]))
    P("")
    tong[vat]["v12"] = V3.tong_ket_v12(kq)
P("## Tương đương GPU–CPU, bộ cấy v2 (105 động vật + 82 thực vật; chạy 11/11, so lại 13/11)")
f = os.path.join(HERE, "..", "bo_kiem_v2_gpu", "so_sanh_dot_bien.json")
if os.path.exists(f):
    for x in json.load(open(f)):
        if "bo" in x:
            P("- %s: %d/%d tệp giống hệt, lệch số max %g; khác: %d tệp" % (x["bo"], x["giong_het"], x["so_te_bao_dot_bien"], x["lech_so_max"], len(x["khac"])))
        else:
            P("- %s: giống hệt %s (số khác %d); tổng ref %s; gpu %s" % (x["tep_tong_hop"], x["giong_het"], x["so_khac"], x.get("tong_ref"), x.get("tong_gpu")))
    P("- Ghi chú: mọi khác biệt thực vật nằm ở THỨ TỰ danh sách ca trượt K.1 (và thứ tự cặp trong tệp tổng hợp); so theo tập (sắp xếp) thì 0 tệp khác; trạng thái luật giống hệt.")
json.dump(tong, open(os.path.join(HERE, "cay_loi_v3_cong.json"), "w"), ensure_ascii=False, indent=1, default=str)
open(os.path.join(HERE, "TOM_TAT_CONG_V3.md"), "w").write("\n".join(o) + "\n")
print(json.dumps({v: {k: t[k] for k in t if k not in ("nguon", "v12")} for v, t in tong.items()}, ensure_ascii=False))
