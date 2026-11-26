"""Hinh cap seed K08: moi diem la mot seed; truc x = doi chung, truc y = PA-PLACE (U5b); duoi duong cheo = PA-PLACE tot hon.
Nguon duy nhat: so_lieu/nguon/ket_qua_cx2.json (DA ra cheo: RA_SOAT_CX2_20261123.md, khop tung o). Chay: python3 hinh/fig4_paired_seeds.py"""
import json, os, hashlib, sys
VI = "--vi" in sys.argv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
D = os.path.dirname(os.path.abspath(__file__)); F = os.path.join(D, "..", "so_lieu", "nguon", "ket_qua_cx2.json")
R = json.load(open(F)); C = R["cap"]
def seed(cap, th, side):
    d = C[cap][th]["theo_seed"]; return [d[s][side] for s in sorted(d)]
TH = [("HD_max_tat_ca", "HD max (µm)"), ("Vmat_tat_ca_pt", "thể tích mất (%)" if VI else "lost volume (%)")]
DC = [("P8", "CXe-P8", 1, "#2a78d6", "o_rong"), ("P4", "CXe-P4", 1, "#eb6834", "o"),
      ("CXe", "U5b-CXe", 1, "#1baf7a", "s"), ("CX2e", "U5b-CX2e", 1, "#1baf7a", "^"), ("CX2b", "U5b-CX2b", 1, "#1baf7a", "D")]
fams = [f.name for f in font_manager.fontManager.ttflist]
plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "serif", "font.serif": ["TeX Gyre Termes", "Times New Roman", "Times", "DejaVu Serif"],
                     "font.size": 10, "axes.edgecolor": "#52514e", "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e", "ytick.color": "#52514e"})
fig, ax = plt.subplots(1, 2, figsize=(6.6, 3.0))
for a, (th, nhan) in zip(ax, TH):
    u = seed("U5b-CXe", th, 0)
    for cap in ("U5b-CX2e", "U5b-CX2b"): assert seed(cap, th, 0) == u
    lo, hi = 0, 0
    for ten, cap, side, col, mk in DC:
        x = seed(cap, th, side); hi = max(hi, max(x), max(u))
        if mk == "o_rong":
            a.scatter(x, u, s=70, facecolors="none", edgecolors=col, linewidths=1.6, marker="o", label=ten, zorder=4)
        else:
            a.scatter(x, u, s=30, c=col, marker=mk, edgecolors="#fcfcfb", linewidths=0.8, label=ten, zorder=3)
    hi *= 1.05
    a.plot([0, hi], [0, hi], ls="--", lw=1, c="#9a9893", zorder=1)
    a.set_xlim(0, hi); a.set_ylim(0, hi); a.set_aspect("equal")
    a.set_xlabel(f"{'Đối chứng' if VI else 'Baseline'}, {nhan}"); a.set_ylabel(f"PA-PLACE, {nhan}")
    a.grid(True, lw=0.4, c="#e4e3de", zorder=0)
    for s in ("top", "right"): a.spines[s].set_visible(False)
    a.text(0.62 * hi, 0.78 * hi, "y = x", ha="right", va="bottom", fontsize=9, color="#52514e", rotation=45, rotation_mode="anchor")
ax[1].legend(frameon=False, loc="upper left", fontsize=9, handletextpad=0.3, borderaxespad=0.2)
for a, nh in zip(ax, ("(a)", "(b)")): a.set_title(nh, loc="left", fontsize=10, fontweight="bold")
if VI:
    plt.rcParams["axes.formatter.use_locale"] = False
    from matplotlib.ticker import FuncFormatter
    for a in ax:
        for axis in (a.xaxis, a.yaxis): axis.set_major_formatter(FuncFormatter(lambda v, _: (f"{v:g}").replace(".", ",")))
fig.tight_layout()
out = os.path.join(D, "fig4_paired_seeds_vi.pdf" if VI else "fig4_paired_seeds.pdf"); fig.savefig(out); fig.savefig(out.replace(".pdf", ".png"), dpi=200)
print("ghi", out, "| md5 nguon", hashlib.md5(open(F, "rb").read()).hexdigest()[:8])
