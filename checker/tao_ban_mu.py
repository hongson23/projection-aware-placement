# -*- coding: utf-8 -*-
"""tao_ban_mu.py — dựng BẢN SAO MÙ để chấm lượt chính (đăng ký trước v12.1 §9: "chấm mù trên bản sao đã xoá siêu dữ liệu nhánh;
chấm cả trước và sau chiếu"; bỏ chấm các việc kc).

Dùng:  python tao_ban_mu.py <thu_muc_kq> <thu_muc_ra> [--bo _kc_] [--vat dong_vat]
  <thu_muc_kq>  vd bo_dung_chung/v121/kq_chinh_dv (mỗi thư mục con: mesh_vo_bao.npz [sau chiếu], mesh_vo_bao_truoc.npz [trước chiếu],
                thong_tin.json)
  <thu_muc_ra>  vd bo_kiem_v3/cham_chinh_dv → tạo
     mu_sau/<ID>/   mesh_vo_bao.npz  = LIÊN KẾT CỨNG tới mesh_vo_bao.npz (không tốn đĩa, không lộ đường dẫn)
                    thong_tin.json   = CHỈ trường "tu_the" (cho G_c); mọi trường khác (nhánh, seed, ô, giờ, ngân sách…) bị bỏ
     mu_truoc/<ID>/ như trên với mesh_vo_bao_truoc.npz
     ds_sau.txt, ds_truoc.txt  — danh sách cho cham_lo.sh, thứ tự XÁO NGẪU NHIÊN
     khoa/KHOA_MU.json (chmod 600) — ID → {nguồn, bản}; CHỈ mở khi phân tích sau khi chấm xong.
ID = 12 ký tự hex ngẫu nhiên (secrets), riêng cho từng (tế bào, bản). Không ghi đè thư mục ra đã có."""
import os, sys, json, secrets, random, argparse, shutil

ap = argparse.ArgumentParser()
ap.add_argument("kq"); ap.add_argument("ra"); ap.add_argument("--bo", default="_kc_"); ap.add_argument("--vat", default="dong_vat")
a = ap.parse_args()
kq = os.path.abspath(a.kq); ra = os.path.abspath(a.ra)
if os.path.exists(ra):
    sys.exit("thu muc ra da co: %s (khong ghi de)" % ra)
te_bao = sorted(d for d in os.listdir(kq) if os.path.isdir(os.path.join(kq, d)))
bo = [d for d in te_bao if a.bo and a.bo in d]; te_bao = [d for d in te_bao if d not in bo]
khoa = {}; ds = {"sau": [], "truoc": []}; thieu = []
for ban, tep in (("sau", "mesh_vo_bao.npz"), ("truoc", "mesh_vo_bao_truoc.npz")):
    os.makedirs(os.path.join(ra, "mu_" + ban))
    for d in te_bao:
        src = os.path.join(kq, d, tep); tt = os.path.join(kq, d, "thong_tin.json")
        if not os.path.exists(src):
            thieu.append([d, ban]); continue
        i = secrets.token_hex(6)
        while i in khoa: i = secrets.token_hex(6)
        o = os.path.join(ra, "mu_" + ban, i); os.makedirs(o)
        try:
            os.link(src, os.path.join(o, "mesh_vo_bao.npz"))
        except OSError:
            shutil.copy2(src, os.path.join(o, "mesh_vo_bao.npz"))
        if os.path.exists(tt):
            t = json.load(open(tt))
            if "tu_the" in t:
                json.dump({"tu_the": t["tu_the"]}, open(os.path.join(o, "thong_tin.json"), "w"), ensure_ascii=False)
        khoa[i] = {"nguon": d, "ban": ban}
        ds[ban].append(o)
os.makedirs(os.path.join(ra, "khoa"))
kf = os.path.join(ra, "khoa", "KHOA_MU.json")
json.dump({"kq": kq, "bo_qua": bo, "thieu": thieu, "id": khoa}, open(kf, "w"), ensure_ascii=False, indent=1)
os.chmod(kf, 0o600)
for ban in ("sau", "truoc"):
    xs = ds[ban][:]; random.SystemRandom().shuffle(xs)
    with open(os.path.join(ra, "ds_%s.txt" % ban), "w") as f:
        f.write("# ban sao mu (%s chieu), %d te bao, tao boi tao_ban_mu.py\n" % (ban, len(xs)))
        for o in xs:
            f.write("%s %s\n" % (o, a.vat))
print("te bao: %d (bo %d: %s); sau %d, truoc %d; thieu %s; khoa: %s" % (len(te_bao), len(bo), a.bo, len(ds["sau"]), len(ds["truoc"]), thieu, kf))
