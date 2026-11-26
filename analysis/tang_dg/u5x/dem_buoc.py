# -*- coding: utf-8 -*-
"""dem_buoc.py (22/11/2026, PI duyệt 09:15 "ok") — ngân sách CX theo SỐ ĐẾM, không theo giây.
so_buoc(d): số bước Adam thật của một việc PLACE (tu_dong), suy ĐÚNG từ bộ đếm danh_gia_phat của bộ dựng:
  danh_gia = K·buoc·L + R·(buoc·n + n) + K·(L + n)   (L = số bộ phận lớn, n = số bộ phận nhỏ, R = tổng số vòng nhỏ đã chạy)
  ⇒ bước = K·buoc + R·buoc.
chi_phi(d, hs): bước tương đương = bước Adam + Σ sự kiện CX × hằng số mức (nhỏ: hs['nho'][N:lap], lớn: hs['lon'][N:lap]) + K lần chiếu chọn cấu hình × hs['nho'];
  lần chiếu chỉ để đo giờ (cx_do_gio) KHÔNG tính. Hằng số ghi trong u5x/m_theo_muc.json TRƯỚC khi chọn (chon_m.py).
  U5 (đối xứng, điều phối 22/11): việc có bao['u5'] cộng thêm (giay_du_bao + giay_vanh)/hs['t_buoc'] bước (chi phí Ĝ/vành mỗi cấu hình).
sach(id) (điều phối 22/11, rà CX_dem): chỉ tính chồng với tiến trình TÍNH TOÁN (python, botri*, cham_mot, julia/packmol, blender, cellpack, …),
  bỏ tiến trình hệ thống (tailscaled, falco, Xorg, gdm, systemd*, containerd, dockerd, calico, …); sạch ⇔ ≥ 3 mẫu và tỉ lệ mẫu có chồng tính toán ≤ 5 %."""
import re
TINH = re.compile(r"python|botri|cham_mot|julia|packmol|blender|cellpack|igl|torch|ncnn|ffmpeg|matlab|\bR\b|Rscript", re.I)
HE = re.compile(r"tailscaled|falco|Xorg|gdm|systemd|containerd|dockerd|calico|snapd|sshd|kworker|journald|rsyslog|gnome|pulseaudio|nvidia-persistenced|irqbalance|cron", re.I)
def la_tinh(cmd): return bool(TINH.search(cmd)) and not HE.search(cmd)
import os, json
R0 = os.environ.get("CX_R", os.path.expanduser("~/Vision_team_data_189/sonhh20/ch3_tebao_dep")); O = f"{R0}/bo_dung_chung/tang_dg"
_VAI = {}
def _vai(dac_ta):
    if dac_ta not in _VAI: _VAI[dac_ta] = {b["ten"]: b["vai"] for b in json.load(open(f"{R0}/{dac_ta}"))["bo_phan"]}
    return _VAI[dac_ta]
def so_buoc(d):
    b = d["bao_cao_bo_tri"]; ns = b["ngan_sach_nhanh"]; K, bu = ns["K"], ns["buoc"]; v = _vai(d["dac_ta"])
    L = sum(c for l, c in d["so_luong"].items() if v[l] == "lon"); n = sum(c for l, c in d["so_luong"].items() if v[l] == "nho")
    r = (b["danh_gia_phat"] - K * bu * L - K * (L + n)) / (bu * n + n); assert abs(r - round(r)) < 1e-6 and 0 <= round(r) <= K * 4, ("không suy được số vòng", r)
    return K * bu + int(round(r)) * bu
def muc(c): return f"{c['N']}:{c['lap']}"
def chi_phi(d, hs):
    S = so_buoc(d); c = d["bao_cao_bo_tri"].get("cx")
    if not c: return chi_phi_u5(d, hs)
    m_ = muc(c); ev = c.get("su_kien") or []
    nl = sum(1 for e in ev if e["pha"] == "lon"); nn = len(ev) - nl; nsel = len(c.get("cau_hinh") or []) if c.get("chon") == "gc" else 0
    return S + nl * hs["lon"][m_] + (nn + nsel) * hs["nho"][m_]
def chi_phi_u5(d, hs):
    u = d["bao_cao_bo_tri"].get("u5")
    return so_buoc(d) + ((u.get("giay_du_bao", 0) + u.get("giay_vanh", 0)) / hs["t_buoc"] if u else 0.0)
def sach(id_):
    f = f"{O}/u5x/loi_cpu/{id_}.jsonl"
    if not os.path.exists(f): return None
    xs = [json.loads(l) for l in open(f) if l.strip()]
    nt = sum(1 for x in xs if any(la_tinh(c[-1]) for c in x.get("chong", [])))
    return {"mau": len(xs), "chong_moi_loai": sum(1 for x in xs if x["n_chong"] > 0), "chong_tinh_toan": nt, "ti_le": nt / max(len(xs), 1), "sach": len(xs) >= 3 and nt <= 0.05 * len(xs)}
