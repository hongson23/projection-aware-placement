# -*- coding: utf-8 -*-
"""ghep_v2.py (24/11/2026, bản in; v3: thêm "co_chu" trong cấu hình = hệ số cỡ chữ) — ghép lại R1/R2/R3 từ các PNG Blender ĐÃ CÓ (không render lại) theo yêu cầu phiên viết:
chữ ≥ 7 pt ở bề rộng in 16 cm (ảnh W = 3150 px ⇒ 7 pt ≈ 49 px; mọi cỡ chữ ở đây ≥ 52 px), chính tả Mỹ, nhãn ô (a)(b)…, xuất PNG + TIFF LZW 500 dpi.
Dùng: python ghep_v2.py <cau_hinh.json>
cau_hinh: {"ra": <tên không đuôi>, "W": 3150, "tieu_de": [..]|null, "phu_de": [..]|null,
  "hang": [{"kieu": "anh", "nhan": str, "anh": [png/cột], "o": ["(a)", ...]}
           | {"kieu": "phong", "nhan": str, "anh": [[png ô...]/cột], "nhan_o": [[str...]/cột], "o": ["(c)", ...]}],
  "diem": [[[x,y,z]...]/cột] (đánh số trên hàng "anh" đầu tiên), "elev", "ortho", "res", "tam",
  "lut", "vmax", "buoc_vach", "nhan_thang", "ghi_chu": [dòng...]}
Hàng "anh": cắt theo HỘP BAO CHUNG của các cột (cùng camera ⇒ cùng tỉ lệ)."""
import sys, json, math, os
import numpy as np
from PIL import Image, ImageDraw, ImageFont
C = json.load(open(sys.argv[1]))
FD = os.path.join(os.path.dirname(__import__("matplotlib").__file__), "mpl-data/fonts/ttf/")
def font(sz, b=False): return ImageFont.truetype(FD + ("DejaVuSans-Bold.ttf" if b else "DejaVuSans.ttf"), sz)
SC = C.get("co_chu", 1.0)
F_TD, F_PD, F_NH, F_VA, F_O, F_GC = (int(round(x * SC)) for x in (84, 60, 62, 56, 66, 54))
def trang(im):
    bg = Image.new("RGBA", im.size, (255, 255, 255, 255)); bg.alpha_composite(im.convert("RGBA")); return bg.convert("RGB")
W = C.get("W", 3150); M = 40; GAP = 50; TG = 20; RG = 50
td = C.get("tieu_de"); nc = len(td) if td else len(C["hang"][0]["anh"])
LBL = (F_NH + 40) if any(h.get("nhan") for h in C["hang"]) else 0
CW = (W - 2 * M - LBL - (nc - 1) * GAP) // nc
H_HEAD = (F_TD + 20 + (F_PD + 30 if C.get("phu_de") else 0) + 20) if td else 0
rows = []
for h in C["hang"]:
    if h["kieu"] == "anh":
        ims = [Image.open(f).convert("RGBA") for f in h["anh"]]; bb = None
        for im in ims:
            b = im.getchannel("A").point(lambda a: 255 if a > 8 else 0).getbbox()
            bb = b if bb is None else (min(bb[0], b[0]), min(bb[1], b[1]), max(bb[2], b[2]), max(bb[3], b[3]))
        p = 16; bb = (max(bb[0] - p, 0), max(bb[1] - p, 0), min(bb[2] + p, ims[0].width), min(bb[3] + p, ims[0].height))
        s = CW / (bb[2] - bb[0]); rows.append((h, ims, bb, int(round((bb[3] - bb[1]) * s)), s))
    else:
        nt = len(h["anh"][0]); TW = (CW - (nt - 1) * TG) // nt; rows.append((h, None, None, TW, TW))
gc = C.get("ghi_chu") or []
H_CB = 60 + 70 + 20 + F_VA + 20
H = M + H_HEAD + sum(r[3] for r in rows) + RG * len(rows) + H_CB + (F_GC + 16) * len(gc) + M
out = Image.new("RGB", (W, H), "white"); dr = ImageDraw.Draw(out)
def x_cot(i): return M + LBL + i * (CW + GAP)
def nhan_o(x, y, t):
    f = font(F_O, True); w = dr.textlength(t, font=f); dr.rectangle((x, y, x + w + 16, y + F_O + 14), fill="white"); dr.text((x + 8, y + 4), t, fill="black", font=f)
if td:
    for i, t in enumerate(td):
        dr.text((x_cot(i) + CW / 2, M), t, fill="black", font=font(F_TD, True), anchor="mt")
        if C.get("phu_de"): dr.text((x_cot(i) + CW / 2, M + F_TD + 20), C["phu_de"][i], fill="#222222", font=font(F_PD), anchor="mt")
y = M + H_HEAD; da_diem = False
EL = math.radians(C.get("elev", 15)); R = C.get("res", 1600); OR = C.get("ortho", 22)
right = np.array([1.0, 0, 0]); up = np.array([0, math.sin(EL), math.cos(EL)]); tam = np.array(C.get("tam", [0, 0, 0]), float)
for h, ims, bb, rh, s in rows:
    if h["kieu"] == "anh":
        for i, im in enumerate(ims):
            out.paste(trang(im.crop(bb)).resize((CW, rh), Image.LANCZOS), (x_cot(i), y))
            if C.get("diem") and not da_diem:
                for k, pt in enumerate(C["diem"][i]):
                    q = np.asarray(pt, float) - tam; u = R / 2 + q @ right / OR * R; v = R / 2 - q @ up / OR * R
                    u = x_cot(i) + (u - bb[0]) * s; v = y + (v - bb[1]) * s; r = int(30 * SC)
                    dr.ellipse((u - r, v - r, u + r, v + r), outline="#0072B2", width=7)
                    dr.text((u + r + 2, v - r + 4), str(k + 1), fill="#0072B2", font=font(F_O, True), anchor="lb")
            if h.get("o"): nhan_o(x_cot(i), y, h["o"][i])
        if C.get("diem"): da_diem = True
    else:
        TW = rh
        for i, tiles in enumerate(h["anh"]):
            for k, zf in enumerate(tiles):
                tx = x_cot(i) + k * (TW + TG)
                out.paste(trang(Image.open(zf)).resize((TW, TW), Image.LANCZOS), (tx, y))
                dr.rectangle((tx, y, tx + TW - 1, y + TW - 1), outline="#888888", width=3)
                if h.get("nhan_o"):
                    t = h["nhan_o"][i][k]; f = font(F_VA, True); w = dr.textlength(t, font=f)
                    yy = y + TW - F_VA - 20 if C.get("nhan_o_duoi") else y + 4
                    dr.rectangle((tx + TW - w - 22, yy, tx + TW - 4, yy + F_VA + 12), fill="white")
                    dr.text((tx + TW - 12, yy + 4), t, fill="black", font=f, anchor="rt")
            if h.get("o"): nhan_o(x_cot(i), y, h["o"][i])
    if h.get("nhan"):
        f = font(F_NH, True); tw = int(dr.textlength(h["nhan"], font=f)); tim = Image.new("RGB", (tw + 10, F_NH + 20), "white")
        ImageDraw.Draw(tim).text((5, 5), h["nhan"], fill="black", font=f); tim = tim.rotate(90, expand=True)
        out.paste(tim, (M, int(y + rh / 2 - tim.height / 2)))
    y += rh + RG
lut = json.load(open(C["lut"]))["lut"]; vmax = C["vmax"]
f_lab = font(F_NH, True); lab = C.get("nhan_thang", "HDᵢ (µm)"); lw = dr.textlength(lab, font=f_lab)
bw = min(1700, W - 2 * M - int(lw) - 260); bx = int((W - bw + lw + 30) / 2) - 60; yb = y + 10
for k in range(bw): dr.line((bx + k, yb, bx + k, yb + 60), fill=lut[int(round(k / (bw - 1) * (len(lut) - 1)))])
dr.rectangle((bx, yb, bx + bw, yb + 60), outline="black", width=2)
st = C.get("buoc_vach", 0.05); t = 0.0; fv = font(F_VA)
while t <= vmax + 1e-9:
    u = bx + t / vmax * bw; dr.line((u, yb + 60, u, yb + 80), fill="black", width=3)
    dr.text((u, yb + 86), ("%.2f" % t), fill="black", font=fv, anchor="mt"); t += st
if abs(round(vmax / st) * st - vmax) > 1e-6:
    dr.line((bx + bw, yb + 60, bx + bw, yb + 80), fill="black", width=3)
    dr.text((bx + bw + 14, yb + 30), "%.3f" % vmax, fill="black", font=fv, anchor="lm")
dr.text((bx - 24, yb + 30), lab, fill="black", font=f_lab, anchor="rm")
yg = yb + H_CB - 10
for k, l in enumerate(gc): dr.text((W / 2, yg + k * (F_GC + 16)), l, fill="#222222", font=font(F_GC), anchor="mt")
out.save(C["ra"] + ".png", dpi=(500, 500)); out.save(C["ra"] + ".tif", compression="tiff_lzw", dpi=(500, 500))
print("GHEP_V2", C["ra"], out.size, "cm in o 500 dpi:", round(out.size[0] / 500 * 2.54, 2), "x", round(out.size[1] / 500 * 2.54, 2))
