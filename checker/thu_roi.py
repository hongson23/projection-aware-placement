import sys, time, numpy as np
sys.path.insert(0, "/home/thainq97/Vision_team_data_189/sonhh20/ch3_tebao_dep/bo_kiem_v3")
import kiem3, kiem2 as K, kiem_thucvat as T
tb = T.TeBaoTV(sys.argv[1])
bs = [b for b in tb.bp.values() if b.loai not in ("vach", "mang")]
n = n_roi = sai = 0; t_roi = t_bool = 0.0; gmax = 0.0
for a, b in T.cap_ung_vien(bs, bs, 1e-9, cung=True):
    if a.loai == b.loai: continue
    n += 1; t = time.time(); r = kiem3.chac_roi(a.luoi, b.luoi); t_roi += time.time() - t
    if r:
        n_roi += 1; t = time.time(); g = T.giao._goc(a, b); t_bool += time.time() - t; gmax = max(gmax, g); sai += g > 0
print({"cap": n, "chung_minh_roi": n_roi, "trong_do_bool_khac_0": sai, "giao_max": gmax, "t_loc_s": round(t_roi, 1), "t_bool_tren_cap_roi_s": round(t_bool, 1)})
