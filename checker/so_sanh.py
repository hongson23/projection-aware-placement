# -*- coding: utf-8 -*-
"""so_sanh.py — đối chiếu hai kiem2.json (bản gốc CPU vs bản tăng tốc).
Cổng: (1) trạng thái từng luật giống hệt; (2) danh sách bộ phận trượt (ca_truot: bo_phan + nguong) giống hệt;
(3) mọi giá trị số khác nhau ≤ 1e-6 (tuyệt đối); (4) mọi chuỗi/bool giống hệt. Bỏ qua trường thời gian.
Dùng: python so_sanh.py goc.json moi.json [--json]   → in một dòng tóm tắt; mã thoát 0 nếu khớp."""
import sys, json, math

BO_QUA = {"thoi_gian_s", "thoi_gian_cap_s", "te_bao"}
TOL = 1e-6


def duyet(a, b, duong, khac, so_max):
    if isinstance(a, dict) and isinstance(b, dict):
        ka = set(a) - BO_QUA; kb = set(b) - BO_QUA
        if ka != kb:
            khac.append((duong, "khoa", sorted(ka ^ kb)[:5]))
        for k in sorted(ka & kb):
            duyet(a[k], b[k], duong + "/" + str(k), khac, so_max)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            khac.append((duong, "do_dai", (len(a), len(b))))
        for i, (x, y) in enumerate(zip(a, b)):
            duyet(x, y, duong + "[%d]" % i, khac, so_max)
    elif isinstance(a, bool) or isinstance(b, bool):
        if a != b:
            khac.append((duong, "bool", (a, b)))
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if (isinstance(a, float) and math.isnan(a)) and (isinstance(b, float) and math.isnan(b)):
            return
        d = abs(a - b) if not (math.isinf(a) or math.isinf(b)) else (0.0 if a == b else float("inf"))
        so_max[0] = max(so_max[0], d)
        if d > TOL:
            khac.append((duong, "so", (a, b)))
    else:
        if a != b:
            khac.append((duong, "gia_tri", (str(a)[:80], str(b)[:80])))


def so(goc, moi):
    A = json.load(open(goc)); B = json.load(open(moi))
    LA = {r["id"]: r for r in A["luat"]}; LB = {r["id"]: r for r in B["luat"]}
    tt = [i for i in LA if LA[i]["trang_thai"] != LB.get(i, {}).get("trang_thai")]
    def hong(r):
        return sorted(json.dumps(c.get("bo_phan"), ensure_ascii=False) + "|" + str(c.get("nguong")) for c in r["ca_truot"])
    ds = [i for i in LA if i in LB and hong(LA[i]) != hong(LB[i])]
    khac = []; so_max = [0.0]
    duyet(A, B, "", khac, so_max)
    return {"luat": len(LA), "khac_trang_thai": tt, "khac_ds_truot": ds, "so_khac_bat_ky": len(khac),
            "lech_so_max": so_max[0], "vd_khac": khac[:5], "khop": not tt and not ds and not khac,
            "dem_goc": A.get("dem_trang_thai"), "dem_moi": B.get("dem_trang_thai"),
            "t_goc": A.get("thoi_gian_s"), "t_moi": B.get("thoi_gian_s")}


if __name__ == "__main__":
    r = so(sys.argv[1], sys.argv[2])
    print(json.dumps(r, ensure_ascii=False))
    sys.exit(0 if r["khop"] else 1)
