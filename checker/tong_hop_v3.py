# -*- coding: utf-8 -*-
"""Gộp kết quả bộ cấy lỗi bản 3 chạy tách theo tế bào (hn_<tb>/cay_loi*_ket_qua.json) → cay_loi_v3_tong.json."""
import os, sys, json, glob, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import v3_cay_loi_chung as V3

out = {}
for vat, tep in (("dong_vat", "cay_loi_ket_qua.json"), ("thuc_vat", "cay_loi_tv_ket_qua.json")):
    kq = []; nguon = []
    for f in sorted(glob.glob(os.path.join(HERE, "hn_*", tep))):
        d = json.load(open(f)); kq += d["ket_qua"]; nguon.append(os.path.relpath(f, HERE))
    db = [x for x in kq if x["kieu"] == "dot_bien" and "phat_hien" in x]
    dc = [x for x in kq if x["kieu"] in ("trung_tinh", "sat_nguong_dat") and "bao_dong_gia" in x]
    phu = collections.defaultdict(lambda: [0, 0])
    for x in db:
        for rid in x["luat_dich"]:
            pq = rid.rsplit("@", 1)[-1] if "@" in rid else None
            ok = x.get("phat_hien_theo_phan_quyet", {}).get(pq) if pq else x["phat_hien"]
            phu[rid][1] += 1; phu[rid][0] += int(bool(ok))
    out[vat] = {"nguon": nguon, "so_ban_ghi": len(kq),
                "khong_hop_le": sum(1 for x in kq if "khong_hop_le" in x), "loi_cay": sum(1 for x in kq if "loi_cay" in x),
                "doi_chung": len(dc), "doi_chung_bao_gia": sum(1 for x in dc if x["bao_dong_gia"]),
                "doi_chung_bao_gia_ds": [(x["te_bao"], x["loai"], x.get("luat_bao_gia") or x.get("chi_tiet")) for x in dc if x["bao_dong_gia"]][:20],
                "v12": V3.tong_ket_v12(kq), "do_phu_theo_muc": {k: "%d/%d" % tuple(v) for k, v in sorted(phu.items())}}
json.dump(out, open(os.path.join(HERE, "cay_loi_v3_tong.json"), "w"), ensure_ascii=False, indent=1, default=str)
print(json.dumps({v: {k: o[k] for k in ("so_ban_ghi", "khong_hop_le", "loi_cay", "doi_chung", "doi_chung_bao_gia")} | {"v12": {k: o["v12"][k] for k in o["v12"] if not k.startswith("bo_sot") and k != "M.c0"}} for v, o in out.items()}, ensure_ascii=False, indent=1, default=str))
