"""botri_chung_v122e3 (19/11/2026) = botri_chung_v122 (md5 8abb7b1c…) + E3: nạp tư thế đã lưu (--tu_the_tu), sửa phần nối
  ngoài hộp atlas của eval_at bằng atlas xa hai mức (--eval_dung / --eval_dung_tu_dau), --chi_P. Mọi cờ tắt ⇒ cùng phép tính như v12.2.
  Sinh bởi bo_dung_chung/tang_dg/tao_ma_e3.py.
botri_chung_v122 (14/11/2026) = botri_chung_v121 (md5 06d37cd6…) + ĐÚNG MỘT thay đổi (ứng viên, chưa dùng cho lượt chính):
  mọi cặp KHÁC loại trong roi_chieu.mien_cap (động vật: lưới nội chất–nhân) — M.c₀ vẫn áp dụng (s = 0, được chạm) — nhận hàng không giao
  u_a + u_b ≥ τ (τ = tol_con = 1e-6; không có s, không có β) thay vì bị bỏ hẳn; giữ được chạm 4.4a, cấm xuyên. Ngoài ra chỉ đổi nhãn 'ban', thêm số hàng
  mới vào dòng SO_HANG và mục báo cáo 'ban_v122'. Sinh bởi bo_dung_chung/v122/tao_ma_v122.py.
botri_chung_v121 (13/11/2026) = botri_chung_v12 (md5 c10b69cd…) + ĐÚNG HAI thay đổi (ứng viên, chưa dùng cho lượt chính):
  (a) vỏ 'cau_nhieu' (động vật): trường màng = khoảng cách có dấu ĐÚNG tới CÙNG mặt màng (mặt mức 0 của trường cũ |P| − R + nhiễu,
      trích một lần bằng marching cubes bước 0,04 µm, khoảng cách tới mesh bằng igl, dấu theo trường cũ), lấy mẫu trên lưới h1/2 thẳng hàng
      với lưới chiếu (nút chiếu = nút lưới ⇒ ψ_vo tại nút đúng), nội suy ba tuyến tính cho PLACE; L_M = 1 ⇒ β_M = β. Lưu sẵn theo băm.
  (b) bước Adam của PLACE: lr = adam_lr · min(1, h1 / 0.1796875) (= 0,05 với động vật h1 = 23/128; 0,0139 với thực vật h = 0,05).
  Ngoài ra chỉ đổi nhãn 'ban' và thêm mục báo cáo 'ban_v121'.
Bản v12:
botri_chung_v12 (11/11/2026) = botri_chung 802cf6f0 + sửa cho lượt chính Bài 1 (luật sự thật v1.2):
  (1) PROJECT dừng theo KKT (khả thi nguyên thuỷ + bù trừ), lượt 2–3 chiếu lại TỪ ψ, kiểm khả thi đóng Mệnh đề 4 mỗi nút, biên hàng màng là tham số (--chieu kkt; --chieu cu = bản cũ);
  (2) ba nhánh đối xứng: cùng số mẫu kiểm nhận, cùng số hạng mềm (SA thêm luật gần đo bằng mẫu), cùng chọn tốt nhất trong K, cùng cách đo giờ đặt
      (đồng hồ tường trên GPU, có làm nóng + đồng bộ CUDA, không gồm dựng khuôn/đọc tệp), SA theo SỐ VÒNG (--un_vong_lon/--un_vong_nho), thử–loại theo số lô (--so_lo);
  (3) tầng 2 có luật (--tang2): hai trung tử vuông góc (4.8a), chồng túi Golgi (4.5a), dựng giải tích trong vỏ bao đã chiếu;
  (4) --kiem_chieu N: chiếu tham chiếu (mọi hàng, KKT 1e-9) trên N nút để đo khoảng cách tới phép chiếu; --khong_tu_do: bỏ tự đo trên mesh.
Bản gốc:
botri_chung — bộ sắp đặt + chiếu MÃ CHUNG (05/11/2026), tách từ botri2.py (md5 df883002…).
Mã KHÔNG chứa tên loại bộ phận nào: loại, cỡ, số lượng, khuôn, quan hệ (chứa / rời / vùng cấm / gần / lệch tâm),
miễn trừ, trường phụ thuộc (vành + khuyết nón), luật khởi tạo, đặt lại, tư thế cố định và mọi hằng số đều đọc từ
tệp đặc tả (--dac_ta). Thuật toán giữ nguyên botri2: atlas khuôn, bộ phận nhỏ theo lô, Adam trên (t, R), cỡ cố định,
đặt lại bộ phận còn phạt, chọn tốt nhất trong K cấu hình, Dykstra cục bộ 12/24/48 hàng, giữ chỗ; cùng thứ tự gọi số ngẫu nhiên.
Biểu thức khởi tạo trong đặc tả: danh sách tiền tố [phép, đối...]; '$x' là biến lúc chạy, '@x' là hằng số đặc tả.
"""
import argparse, json, math, os, sys, time
import numpy as np, torch, trimesh, igl
import torch.nn.functional as Fn
from scipy import ndimage as ndi
from skimage.measure import marching_cubes

ap = argparse.ArgumentParser()
ap.add_argument("--dac_ta", default="dac_ta/dong_vat_botri2.json")
ap.add_argument("--seed", type=int, default=101); ap.add_argument("--N1", type=int, default=128)
ap.add_argument("--cach", default="tu_dong", choices=["tu_dong", "ngau_nhien", "tay", "u_nhiet"])
ap.add_argument("--thu_vien", default="T")
ap.add_argument("--khong_chieu", action="store_true"); ap.add_argument("--trong_so", type=int, default=1)
ap.add_argument("--K", type=int, default=4); ap.add_argument("--buoc", type=int, default=400)
ap.add_argument("--so", action="append", default=[], help="LOAI=N: ghi đè số lượng của một loại (0 = giữ số rút ngẫu nhiên)")
ap.add_argument("--khoi_tao", default="luat", choices=["luat", "deu"])
ap.add_argument("--out", default="do_botri/chung"); ap.add_argument("--luu_truong", action="store_true")
ap.add_argument("--h", type=float, default=0.0, help="> 0: dùng lưới cục bộ bước h (ghi đè đặc tả)")
ap.add_argument("--luoi_day_du", action="store_true", help="phép thử tương đương: mọi bộ phận lấy hộp = toàn miền")
ap.add_argument("--ngan_sach_giay", type=float, default=0.0, help="u_nhiet: ngân sách thời gian (giây) của TOÀN bước đặt (mọi K cấu hình)")
ap.add_argument("--un_p0", type=float, default=0.5, help="u_nhiet: xác suất nhận bước lên dốc trung bình lúc đầu (định T0)")
ap.add_argument("--un_T1", type=float, default=1e-4, help="u_nhiet: nhiệt độ cuối (tuyệt đối, cùng đơn vị J)")
ap.add_argument("--un_buoc_t", type=float, nargs=2, default=[0.3, 0.02], help="u_nhiet: độ lệch dời (× bán kính bao) đầu, cuối")
ap.add_argument("--un_buoc_r", type=float, nargs=2, default=[0.5, 0.03], help="u_nhiet: góc xoay (rad) đầu, cuối")
ap.add_argument("--un_p_doi", type=float, default=0.5, help="u_nhiet: xác suất bước dời (còn lại là xoay)")
ap.add_argument("--un_ti_le_lon", type=float, default=0.2, help="u_nhiet: phần ngân sách mỗi cấu hình cho bộ phận lớn")
ap.add_argument("--un_lam_moi", type=int, default=50, help="u_nhiet: số vòng giữa hai lần rút lại bộ mẫu cố định")
ap.add_argument("--un_vong_lon", type=int, default=0, help="u_nhiet: số vòng chặng bộ phận lớn MỖI cấu hình (> 0 cùng --un_vong_nho: ngân sách theo vòng, không theo giây)")
ap.add_argument("--un_vong_nho", type=int, default=0, help="u_nhiet: số vòng chặng bộ phận nhỏ MỖI cấu hình")
ap.add_argument("--so_lo", type=int, default=0, help="ngau_nhien: ghi đè thu_loai_so_lo (núm ngân sách thử–loại)")
ap.add_argument("--chieu", default="kkt", choices=["kkt", "cu"], help="kkt: PROJECT v12 (dừng KKT, chiếu lại từ ψ); cu: như botri_chung 802cf6f0")
ap.add_argument("--kkt_tol", type=float, default=1e-6, help="ngưỡng phần dư bù trừ KKT (µm, metric W); khả thi nguyên thuỷ dùng tol_chieu")
ap.add_argument("--lap_kkt", type=int, default=5000, help="số vòng Dykstra tối đa mỗi lượt (chế độ kkt)")
ap.add_argument("--bien_mang", default="", help="hệ số biên hàng màng β_M/β: số, hoặc 'tu_dong' = (1+L_M)/2; rỗng = theo đặc tả (bien_mang_he_so, mặc định 1)")
ap.add_argument("--kiem_chieu", type=int, default=0, help="> 0: chiếu tham chiếu trên chừng ấy nút vi phạm (mọi hàng, KKT 1e-9) và báo khoảng cách tới phép chiếu")
ap.add_argument("--tang2", action="store_true", help="dựng đối tượng con tầng 2 theo mục 'tang2' của đặc tả")
ap.add_argument("--khong_tu_do", action="store_true", help="bỏ các phép tự đo trên mesh (bộ kiểm độc lập đo)")
ap.add_argument("--khong_mem", action="store_true", help="đối chứng: bỏ mọi số hạng luật mềm S(g) khỏi J (PLACE và SA)")
ap.add_argument("--chi_dat", action="store_true", help="chỉ chạy bước sắp đặt rồi dừng (đo giờ khi hiệu chỉnh ngân sách); ghi thong_tin_dat.json")
ap.add_argument("--tu_the_tu", default="", help="E3: nạp tư thế đã lưu (thong_tin.json), bỏ vòng K (chỉ đọc tệp)")
ap.add_argument("--eval_dung", action="store_true", help="E3: sửa phần nối ngoài hộp atlas, bật SAU khi nạp tư thế")
ap.add_argument("--eval_dung_tu_dau", action="store_true", help="E3b: sửa phần nối ngoài hộp atlas từ trước vòng K (đặt + chiếu)")
ap.add_argument("--chi_P", action="store_true", help="E3: chỉ tính phạt cứng đầy đủ dưới eval gốc và eval sửa rồi dừng (cần --tu_the_tu)")
ARG = ap.parse_args()
assert 1 <= ARG.K < 10 ** 6
BUOC_UN = ARG.un_vong_lon > 0 and ARG.un_vong_nho > 0
if ARG.cach == "u_nhiet" and ARG.ngan_sach_giay <= 0 and not BUOC_UN: ap.error("--cach u_nhiet cần --ngan_sach_giay > 0 hoặc --un_vong_lon/--un_vong_nho > 0")

DT0 = json.load(open(ARG.dac_ta)); TS = {**DT0["tham_so_khong_thuoc_luat"], **DT0.get("tham_so_so_hoc", {})}
def giai(x):
    if isinstance(x, str) and x.startswith("@"):
        ten, _, chi = x[1:].partition("#"); v = giai(TS[ten]["gia_tri"]); return v[int(chi)] if chi else v
    if isinstance(x, list): return [giai(v) for v in x]
    if isinstance(x, dict): return {k: giai(v) for k, v in x.items()}
    return x
DT = giai({k: v for k, v in DT0.items() if not k.startswith("tham_so_")}); SH = {k: giai(v["gia_tri"]) for k, v in TS.items()}
if ARG.so_lo > 0: SH["thu_loai_so_lo"] = ARG.so_lo
VO = DT["vo"]; VO_TEN = VO["ten"]; QH = DT["quan_he"]
BP = DT["bo_phan"]; LOAI = {b["ten"]: b for b in BP}
LON = [b["ten"] for b in BP if b["vai"] == "lon"]; NHO = [b["ten"] for b in BP if b["vai"] == "nho"]; PHU = [b["ten"] for b in BP if b["vai"] == "phu_thuoc"]
LID = {l: i for i, l in enumerate(NHO)}
DAI = {l: tuple(LOAI[l]["dai"]) for l in LON + NHO}
LUOI = dict(DT.get("luoi", {"kieu": "toan_phan"}))
if ARG.h > 0: LUOI.update(kieu="cuc_bo", h=ARG.h)
CUCBO = LUOI["kieu"] == "cuc_bo"
DAI_K = {kd["ten"]: tuple(kd["dai_co"]) for kd in DT["khuon"] if "dai_co" in kd}

DEV = "cuda"; F64 = torch.float64; torch.manual_seed(ARG.seed)
if os.environ.get("GB_TOI_DA"): torch.cuda.set_per_process_memory_fraction(min(1.0, float(os.environ["GB_TOI_DA"]) / 23.5))
OUT = ARG.out; os.makedirs(OUT, exist_ok=True); t00 = time.time()
def T(v, dt=F64): return torch.as_tensor(np.asarray(v, float), device=DEV, dtype=dt)
def relu2(x): return torch.clamp(x, min=0) ** 2
def cong(x, ds):
    for d in ds: x = x + d
    return x + DEL

if VO["hinh"] == "cau_nhieu": R_VO = VO["ban_kinh"]; A = VO["nua_hop"]
else: R_VO = None; A = max(VO["nua_canh"]) + 1.0
N1 = ARG.N1; h1 = LUOI["h"] if CUCBO else 2 * A / N1; bu = math.sqrt(SH["bu_can"]) * h1; DEL = bu

G = SH["atlas_G"]; GB = SH["atlas_GB"]; NM = SH["khuon_so_mat"]; NTR = SH["khuon_so_trong"]
gax = torch.linspace(-GB, GB, G, device=DEV, dtype=F64); GP = torch.stack(torch.meshgrid(gax, gax, gax, indexing="ij"), -1).reshape(-1, 3)
def elip(P, a):
    Y = P / T(a); k0 = torch.linalg.norm(Y, dim=1); k1 = torch.linalg.norm(Y / T(a), dim=1).clamp_min(1e-12); return k0 * (k0 - 1) / k1
def ong(P, nua, r, cong_=0.0):
    nd = SH["ong_so_diem"]; t = torch.linspace(-1, 1, nd, device=DEV, dtype=F64)[:, None]
    Q = torch.cat([t * nua, cong_ * (1 - t ** 2), 0 * t], 1); d = torch.full((P.shape[0],), 1e9, device=DEV, dtype=F64)
    for i in range(nd - 1):
        p0, p1 = Q[i], Q[i + 1]; ab = p1 - p0; s = ((P - p0) @ ab / (ab @ ab)).clamp(0, 1); d = torch.minimum(d, torch.linalg.norm(P - (p0 + s[:, None] * ab), dim=1))
    return d - r
def sdf_hinh(P, kd):
    hh = kd["hinh"]
    if hh == "thau_kinh":
        q = kd["nua_day"]; Rr = (q + 1 / q) / 2; cc = (1 / q - q) / 2; x = P[..., 2].abs(); y = torch.linalg.norm(P[..., :2], dim=-1)
        return torch.where((y - 1.0) * cc > x * 1.0, torch.sqrt(x ** 2 + (y - 1.0) ** 2), torch.sqrt((x + cc) ** 2 + y ** 2) - Rr)
    if hh == "nang":
        rr = kd["r"]; px = (P[..., 0].abs() - (kd["nua_dai"] - rr)).clamp_min(0); return torch.sqrt(px ** 2 + P[..., 1] ** 2 + P[..., 2] ** 2) - rr
    if hh == "tru_bo":
        e = kd["bo"]; d1 = torch.linalg.norm(P[..., :2], dim=-1) - (kd["ban_kinh"] - e); d2 = P[..., 2].abs() - (kd["nua_cao"] - e)
        return torch.clamp(torch.maximum(d1, d2), max=0) + torch.sqrt(d1.clamp_min(0) ** 2 + d2.clamp_min(0) ** 2) - e
    if hh == "cau": return torch.linalg.norm(P, dim=-1) - kd["ban_kinh"]
    raise NotImplementedError(hh)
HINH_GT = ("cau", "thau_kinh", "nang", "tru_bo")
KHUON = {}
def them_khuon(ten, f_grid=None, V=None, Fa=None):
    if V is None:
        f = f_grid.reshape(G, G, G); v, fa, _, _ = marching_cubes(f.cpu().numpy().astype(np.float32), 0.0); v = v / (G - 1) * 2 * GB - GB
        m = trimesh.Trimesh(v, fa, process=True)
    else: m = trimesh.Trimesh(V, Fa, process=True)
    if m.volume < 0: m.invert()
    ck = None
    if "atlas_cache" in SH:
        import hashlib; ck = f'{SH["atlas_cache"]}/{hashlib.md5(np.asarray(m.vertices, np.float64).tobytes() + np.asarray(m.faces, np.int64).tobytes() + f"{G}_{GB}".encode()).hexdigest()}.npy'
    if ck and os.path.exists(ck): sdf = T(np.load(ck)).reshape(G, G, G)
    else:
        sdf = T(igl.signed_distance(GP.cpu().numpy(), np.asarray(m.vertices, float), np.asarray(m.faces, np.int64))[0]).reshape(G, G, G)
        if ck: os.makedirs(SH["atlas_cache"], exist_ok=True); np.save(ck, sdf.cpu().numpy())
    mat, _ = trimesh.sample.sample_surface_even(m, NM, seed=SH["khuon_hat_mat"]); mat = np.asarray(mat)
    mat = mat[np.arange(NM) % len(mat)]
    trong = GP[(sdf.reshape(-1) < -SH["khuon_nguong_trong"] * float(-sdf.min()))]; tr = trong[torch.randperm(trong.shape[0], device=DEV)[:NTR]]; tr = tr[torch.arange(NTR, device=DEV) % tr.shape[0]]
    KHUON[ten] = dict(grid=sdf.permute(2, 1, 0)[None, None].float(), mat=T(mat), trong=tr, dien_tich=float(m.area), r_trong=float(-sdf.min()),
                      rb=float(np.linalg.norm(m.vertices, axis=1).max()), sau=GP[int(torch.argmin(sdf.reshape(-1)))], nz=float(np.abs(m.vertices[:, 2]).max()))
    KHUON[ten]["_V"] = np.asarray(m.vertices, float); KHUON[ten]["_F"] = np.asarray(m.faces, np.int64)
for kd in DT["khuon"]:
    if kd["hinh"] == "cau": them_khuon(kd["ten"], torch.linalg.norm(GP, dim=1) - kd["ban_kinh"])
    elif kd["hinh"] == "elip": them_khuon(kd["ten"], elip(GP, kd["truc"]))
    elif kd["hinh"] == "ong": them_khuon(kd["ten"], ong(GP, kd["dai"] - kd["rc"][0], kd["rc"][0], kd["rc"][1]))
    else: them_khuon(kd["ten"], sdf_hinh(GP, kd))
KD = {kd["ten"]: kd for kd in DT["khuon"]}
TV = DT["thu_vien"][ARG.thu_vien]; LIB = {l: list(v) for l, v in TV["co_san"].items()}
if "ai" in TV:
    ai = TV["ai"]
    for d_ in json.load(open(os.environ.get(ai["bien_moi_truong"], ai["tep_mac_dinh"])))["khuon"]:
        lo_ = ai["anh_xa_loai"].get(d_["loai"], d_["loai"])
        if lo_ in ai["bo_qua_loai"]: continue
        z = np.load(f"{ai['thu_muc']}/{d_['ten']}.npz", allow_pickle=True); V = z["v"].astype(float); V = V - (V.max(0) + V.min(0)) / 2
        ext = np.ptp(V, axis=0); V = V[:, np.argsort(-ext)]; ext = np.ptp(V, axis=0); V = V * SH["ai_dai_truc"] / ext.max()
        them_khuon(ai["tien_to"] + d_["ten"], V=V, Fa=np.asarray(z["f"], np.int64))
        LIB[lo_].append(ai["tien_to"] + d_["ten"])
        if d_["ten"] in ai.get("dai_co", {}): DAI_K[ai["tien_to"] + d_["ten"]] = tuple(ai["dai_co"][d_["ten"]])
KHL = {l: (LIB[l] if LOAI[l]["khuon"] == "thu_vien" else [LOAI[l]["khuon"]]) for l in LON + NHO}
print("KHUON", KHL, flush=True)
TEN_K = list(KHUON); TID = {k: i for i, k in enumerate(TEN_K)}; nT = len(TEN_K)
AT = torch.cat([Fn.pad(KHUON[k]["grid"], (0, 0, 0, 0, 1, 1), mode="replicate") for k in TEN_K], dim=2); DT_ = AT.shape[2]
MATS = torch.stack([KHUON[k]["mat"] for k in TEN_K]); TRONGS = torch.stack([KHUON[k]["trong"] for k in TEN_K])
RB = T([KHUON[k]["rb"] for k in TEN_K]); RTR = T([KHUON[k]["r_trong"] for k in TEN_K]); SAU = torch.stack([KHUON[k]["sau"] for k in TEN_K])
def _eval_at_goc(Pl, tid):
    sh = Pl.shape[:-1]; Pl = Pl.reshape(-1, 3); tid = tid.reshape(-1)
    q = Pl.clamp(-GB, GB); zi = tid.to(F64) * (G + 2) + 1 + (q[:, 2] + GB) / (2 * GB) * (G - 1)
    g = torch.stack([q[:, 0] / GB, q[:, 1] / GB, zi / (DT_ - 1) * 2 - 1], 1).float()[None, None, None]
    v = Fn.grid_sample(AT, g, align_corners=True, padding_mode="border").reshape(-1).to(F64)
    a = torch.linalg.norm(Pl - q, dim=1)
    return torch.where(a > 0, torch.sqrt((v.clamp_min(0) ** 2 + a ** 2).clamp_min(1e-12)), v).reshape(sh)
E3 = {"bat": False, "G2": 97, "GB2": [8.0, 64.0], "AT": None, "kiem_xa": None, "giay_xay": None}
def _eval_luoi(AT_, GB_, G_, Pl, tid):
    q = Pl.clamp(-GB_, GB_); DT2 = AT_.shape[2]; zi = tid.to(F64) * (G_ + 2) + 1 + (q[:, 2] + GB_) / (2 * GB_) * (G_ - 1)
    g = torch.stack([q[:, 0] / GB_, q[:, 1] / GB_, zi / (DT2 - 1) * 2 - 1], 1).float()[None, None, None]
    v = Fn.grid_sample(AT_, g, align_corners=True, padding_mode="border").reshape(-1).to(F64)
    a = torch.linalg.norm(Pl - q, dim=1)
    return torch.where(a > 0, torch.sqrt((v.clamp_min(0) ** 2 + a ** 2).clamp_min(1e-12)), v)
def _xay_xa():
    import hashlib as _hl
    if E3["AT"] is not None: return
    t_ = time.time(); G2 = E3["G2"]; ck_dir = os.environ.get("E3_ATLAS_XA", "bo_dung_chung/tang_dg/atlas_xa_cache"); os.makedirs(ck_dir, exist_ok=True); AT2 = []
    for GB2 in E3["GB2"]:
        ax2 = np.linspace(-GB2, GB2, G2); GP2 = np.stack(np.meshgrid(ax2, ax2, ax2, indexing="ij"), -1).reshape(-1, 3); gr = []
        for k in TEN_K:
            V_, F_ = KHUON[k]["_V"], KHUON[k]["_F"]
            ck = f"{ck_dir}/{_hl.md5(V_.tobytes() + F_.tobytes() + f'{G2}_{GB2}'.encode()).hexdigest()}.npy"
            if os.path.exists(ck): sd_ = np.load(ck)
            else: sd_ = igl.signed_distance(GP2, V_, F_)[0].astype(np.float32); _tm = f"{ck}.{os.getpid()}.tam.npy"; np.save(_tm, sd_); os.replace(_tm, ck)
            gr.append(Fn.pad(T(sd_).reshape(G2, G2, G2).permute(2, 1, 0)[None, None].float(), (0, 0, 0, 0, 1, 1), mode="replicate"))
        AT2.append(torch.cat(gr, dim=2))
    E3["AT"] = AT2; E3["giay_xay"] = round(time.time() - t_, 1)
    rg = np.random.default_rng(20261119); kx = {}
    for k in TEN_K:
        V_, F_ = KHUON[k]["_V"], KHUON[k]["_F"]; o = {}
        for ten_, lo_, hi_ in (("vo1", GB, E3["GB2"][0]), ("vo2", E3["GB2"][0], E3["GB2"][1])):
            d_ = rg.normal(size=(3000, 3)); d_ /= np.abs(d_).max(1, keepdims=True); Q = d_ * rg.uniform(lo_, hi_, size=(3000, 1))
            ex = igl.signed_distance(Q, V_, F_)[0]
            Qt = T(Q); tt = torch.full((3000,), TID[k], device=DEV)
            E3["bat"] = True; v1 = eval_at(Qt, tt).cpu().numpy(); E3["bat"] = False; v0 = eval_at(Qt, tt).cpu().numpy()
            o[ten_] = {"sai_max_dung": float(np.abs(v1 - ex).max()), "sai_tb_dung": float(np.abs(v1 - ex).mean()), "sai_max_goc": float(np.abs(v0 - ex).max()),
                       "thieu_max_goc": float((ex - v0).max()), "khoang_min": float(ex.min())}
        Qi = T((rg.uniform(-GB, GB, size=(3000, 3)))); ti = torch.full((3000,), TID[k], device=DEV)
        E3["bat"] = True; a1 = eval_at(Qi, ti); E3["bat"] = False; a0 = eval_at(Qi, ti); o["trong_hop_bang_nhau"] = bool(torch.equal(a0, a1))
        kx[k] = o
    E3["kiem_xa"] = kx
def eval_at(Pl, tid):
    v0 = _eval_at_goc(Pl, tid)
    if not E3["bat"]: return v0
    sh = Pl.shape[:-1]; P_ = Pl.reshape(-1, 3); t_ = tid.reshape(-1); v0 = v0.reshape(-1); m_ = P_.abs().amax(1)
    v1 = _eval_luoi(E3["AT"][0], E3["GB2"][0], E3["G2"], P_, t_); v2 = _eval_luoi(E3["AT"][1], E3["GB2"][1], E3["G2"], P_, t_)
    return torch.where(m_ <= GB, v0, torch.where(m_ <= E3["GB2"][0], v1, v2)).reshape(sh)
if os.environ.get("TU_KIEM"):
    for k in TEN_K:
        q = (torch.rand(5000, 3, device=DEV, dtype=F64) * 2 - 1) * 2.3; g = (q.clamp(-GB, GB) / GB).float()[None, None, None]
        v0 = Fn.grid_sample(KHUON[k]["grid"], g, align_corners=True, padding_mode="border").reshape(-1).to(F64)
        a = torch.linalg.norm(q - q.clamp(-GB, GB), dim=1); v0 = torch.where(a > 0, torch.sqrt((v0.clamp_min(0) ** 2 + a ** 2).clamp_min(1e-12)), v0)
        print("TU_KIEM atlas", k, float((eval_at(q, torch.full((5000,), TID[k], device=DEV)) - v0).abs().max()), flush=True)

def rot6(r):
    a, b = r[..., :3], r[..., 3:]; e1 = a / a.norm(dim=-1, keepdim=True); b = b - (e1 * b).sum(-1, keepdim=True) * e1; e2 = b / b.norm(dim=-1, keepdim=True)
    return torch.stack([e1, e2, torch.cross(e1, e2, dim=-1)], -1)
def phi_khung(P, t, R, S, tid):
    loc = ((P - t)[..., :, None] * R).sum(-2) / S[..., None]; return S * eval_at(loc, tid)

class BoPhan:
    def __init__(s, ten, loai, khuon, t, r, z):
        s.ten, s.loai, s.khuon = ten, loai, khuon; s.smin, s.smax = DAI_K.get(khuon, DAI[loai]); s.tid = TID[khuon]
        s.t = T(t).clone(); s.r = T(r).clone(); s.z = T([z]).clone()
    def sc(s): return s.smin + (s.smax - s.smin) * torch.sigmoid(s.z)[0]
    def R(s): return rot6(s.r[None])[0]
    def phi(s, P): S = s.sc(); return S * eval_at(((P - s.t) @ s.R()) / S, torch.full(P.shape[:-1], s.tid, device=DEV))
    def mau(s, day_du=False):
        Sv = float(s.sc()); n = min(NM, max(SH["mau_toi_thieu"], int(KHUON[s.khuon]["dien_tich"] * Sv * Sv / DEL ** 2)))
        if day_du: idx = torch.arange(NM, device=DEV)
        elif getattr(s, "idx_co_dinh", None) is not None: idx = s.idx_co_dinh
        else: idx = torch.randint(0, NM, (n,), device=DEV)
        return (torch.cat([MATS[s.tid][idx], TRONGS[s.tid]], 0) * s.sc()) @ s.R().T + s.t
    def tam(s): return (SAU[s.tid] * s.sc()) @ s.R().T + s.t
    def r_trong(s): return float(RTR[s.tid] * s.sc())
    def rb(s): return float(RB[s.tid] * s.sc())

MIEN_NHO = set(QH["roi_nho_nho"]["mien_cung_loai"])
class Nhom:
    def __init__(s, ds):
        s.ten = [d[0] for d in ds]; s.loai = [d[1] for d in ds]; n = len(ds); s.n = n
        s.tid = torch.tensor([TID[d[2]] for d in ds], device=DEV); s.lid = torch.tensor([LID[d[1]] for d in ds], device=DEV)
        s.mien = torch.tensor([d[1] in MIEN_NHO for d in ds], device=DEV, dtype=torch.bool)
        s.smin = T([DAI_K.get(d[2], DAI[d[1]])[0] for d in ds]); s.smax = T([DAI_K.get(d[2], DAI[d[1]])[1] for d in ds])
        s.t = T(np.array([d[3] for d in ds])).reshape(n, 3); s.r = T(np.array([d[4] for d in ds])).reshape(n, 6); s.z = T(np.array([d[5] for d in ds])).reshape(n)
    def sc(s, z=None): return s.smin + (s.smax - s.smin) * torch.sigmoid(s.z if z is None else z)
    def R(s): return rot6(s.r)
    def mask(s, loai): return torch.tensor([l == loai for l in s.loai], device=DEV, dtype=torch.bool)
    def diem(s, M, idx=None, t=None, r=None, z=None, tid=None):
        t = s.t if t is None else t; r = s.r if r is None else r; z = s.z if z is None else z; tid = s.tid if tid is None else tid
        n = t.shape[0]; S = DAI_T(tid, z); Rm = rot6(r)
        ii = torch.randint(0, NM, (n, M), device=DEV) if idx is None else idx
        L = torch.cat([MATS[tid[:, None], ii], TRONGS[tid][:, :SH["mau_so_trong_lo"]]], 1)
        return torch.einsum("nij,nmj->nmi", Rm, L * S[:, None, None]) + t[:, None, :]
def DAI_T(tid, z): return LO_T[tid] + (HI_T[tid] - LO_T[tid]) * torch.sigmoid(z)
LO_T = torch.zeros(nT, device=DEV, dtype=F64); HI_T = torch.ones(nT, device=DEV, dtype=F64)
for l in NHO:
    for k in KHL[l]: LO_T[TID[k]], HI_T[TID[k]] = DAI_K.get(k, DAI[l])

def nhieu(P, nh_):
    g = np.random.default_rng(nh_["hat"]); out = torch.zeros(P.shape[:-1], dtype=F64, device=DEV); k0, k1 = nh_["he_so_k"]
    for _ in range(nh_["so_song"]):
        k = g.normal(size=3); k = k / np.linalg.norm(k) * (2 * math.pi / nh_["quy_mo"]) * g.uniform(k0, k1); out = out + torch.sin(P @ T(k) + g.uniform(0, 2 * math.pi))
    return nh_["bien_do"] * out / nh_["chia"]
V121 = {}
if VO["hinh"] == "cau_nhieu":
    def phi_vo_goc(P): return torch.linalg.norm(P, dim=-1) - R_VO + nhieu(P, VO["nhieu"])
    NF = 2 * N1 + 3; LOF = -A - h1 / 2; HIF = A + h1 / 2; HF = (HIF - LOF) / (NF - 1)
    HMC = 0.04; _ext = R_VO + VO["nhieu"]["bien_do"] * VO["nhieu"]["so_song"] / VO["nhieu"]["chia"] + 3 * HMC
    _key = __import__("hashlib").md5(json.dumps({"vo": VO, "A": A, "N1": N1, "HMC": HMC, "ext": _ext, "ma": "v121a"}, sort_keys=True).encode()).hexdigest()[:16]
    _thu = os.environ.get("SDF_MANG_THU_MUC", "bo_dung_chung/v121/sdf_mang"); os.makedirs(_thu, exist_ok=True); _tep = f"{_thu}/sdf_{_key}.npz"
    if not os.path.exists(_tep):
        import fcntl; _kh_ = open(f"{_tep}.khoa", "w"); fcntl.flock(_kh_, fcntl.LOCK_EX)
        if not os.path.exists(_tep):
            _t0 = time.time(); am = torch.arange(-_ext, _ext + HMC / 2, HMC, device=DEV, dtype=F64); nm_ = am.numel(); fm = np.empty((nm_, nm_, nm_), np.float32)
            with torch.no_grad():
                YZ = torch.stack(torch.meshgrid(am, am, indexing="ij"), -1).reshape(-1, 2)
                for i_ in range(nm_): fm[i_] = phi_vo_goc(torch.cat([am[i_].expand(YZ.shape[0], 1), YZ], 1)).reshape(nm_, nm_).cpu().numpy()
            Vm, Fm, _, _ = marching_cubes(fm, 0.0, spacing=(HMC, HMC, HMC)); Vm = Vm - _ext; del fm
            Fm = np.ascontiguousarray(Fm.astype(np.int64)); Vm = np.ascontiguousarray(Vm.astype(np.float64))
            af = np.linspace(LOF, HIF, NF); sdf_f = np.empty((NF, NF, NF), np.float32)
            with torch.no_grad():
                for i_ in range(0, NF, 32):
                    i1 = min(NF, i_ + 32); Q = np.stack(np.meshgrid(af[i_:i1], af, af, indexing="ij"), -1).reshape(-1, 3)
                    d2 = igl.point_mesh_squared_distance(Q, Vm, Fm)[0]; sg = np.sign(phi_vo_goc(T(Q)).cpu().numpy())
                    sdf_f[i_:i1] = (sg * np.sqrt(d2)).reshape(i1 - i_, NF, NF); print("SDF_MANG lat", i1, "/", NF, round(time.time() - _t0, 1), flush=True)
                eV = float(phi_vo_goc(T(Vm)).abs().max())
            np.savez(_tep + ".tmp.npz", sdf=sdf_f, V=Vm.astype(np.float32), F=Fm.astype(np.int32), sai_mat_mc=eV, giay=time.time() - _t0, HMC=HMC)
            os.replace(_tep + ".tmp.npz", _tep)
        fcntl.flock(_kh_, fcntl.LOCK_UN); _kh_.close()
    _z = np.load(_tep); SDF_F = _z["sdf"]
    _dl = max(float(np.abs(np.diff(SDF_F, axis=a_)).max()) for a_ in range(3)) / HF
    V121.update(sdf_tep=_tep, sdf_md5=__import__("hashlib").md5(open(_tep, "rb").read()).hexdigest(), sdf_luoi=[NF, LOF, HIF, HF], mc_buoc=HMC, mc_so_mat=int(_z["F"].shape[0]),
                sai_mat_mc_max=float(_z["sai_mat_mc"]), doc_roi_rac_max=_dl, giay_dung=float(_z["giay"]))
    SDF_G = torch.as_tensor(SDF_F, device=DEV).permute(2, 1, 0)[None, None].contiguous()
    PSI_VO_NUT = torch.as_tensor(np.ascontiguousarray(SDF_F[2:2 * N1 + 2:2, 2:2 * N1 + 2:2, 2:2 * N1 + 2:2]).reshape(-1), device=DEV)
    del _z
    def phi_vo(P):
        sh = P.shape[:-1]; Pl = P.reshape(-1, 3); q = Pl.clamp(LOF, HIF)
        g = ((q - LOF) / (HIF - LOF) * 2 - 1).float()[None, None, None]
        v = Fn.grid_sample(SDF_G, g, align_corners=True, padding_mode="border").reshape(-1).to(F64)
        a = torch.linalg.norm(Pl - q, dim=1)
        return torch.where(a > 0, torch.sqrt((v.clamp_min(0) ** 2 + a ** 2).clamp_min(1e-12)), v).reshape(sh)
elif VO["hinh"] == "hop_bo_tron":
    _hv = SH["hat_vo"]; _rv = np.random.default_rng(ARG.seed * _hv[0] + _hv[1]); LOP_VO = []; _d = 0.0
    for lp in VO["lop"]: _x = float(_rv.uniform(*lp["day"])); LOP_VO.append({"ten": lp["ten"], "ngoai": _d, "day": _x}); _d += _x
    NC = np.array(VO["nua_canh"], float); RBO = VO["ban_kinh_bo"]; D_VO = _d
    def sdf_hop(P, lui):
        b = T(NC - lui); rr = RBO - lui; q = P.abs() - (b - rr)
        return torch.linalg.norm(q.clamp_min(0), dim=-1) + q.max(-1).values.clamp(max=0) - rr
    def phi_vo(P): return sdf_hop(P, D_VO)
else: raise NotImplementedError(VO["hinh"])
def _L_mang():
    return 1.0
    nh_ = VO["nhieu"]; g = np.random.default_rng(nh_["hat"]); k0, k1 = nh_["he_so_k"]; s_ = 0.0
    for _ in range(nh_["so_song"]): g.normal(size=3); s_ += (2 * math.pi / nh_["quy_mo"]) * g.uniform(k0, k1); g.uniform(0, 2 * math.pi)
    return 1.0 + nh_["bien_do"] / nh_["chia"] * s_
L_M = _L_mang(); _bm = ARG.bien_mang or str(SH.get("bien_mang_he_so", 1.0))
HS_M = (1.0 + L_M) / 2 if _bm == "tu_dong" else float(_bm); DEL_M = HS_M * bu
print("BIEN_MANG L_M", round(L_M, 4), "he_so", round(HS_M, 4), "beta", round(bu, 4), "beta_M", round(DEL_M, 4), flush=True)
def cong_b(x, ds, b): return cong(x, ds) + ((DEL_M - DEL) if b == VO_TEN else 0.0)
HT = {}
VANH = {v["ten"]: v for v in DT.get("vanh", [])}
def u_kh(v):
    q = HT[v["quanh"]]; u = HT[v["khuyet"]["huong_toi"]].t - q.t; return u / u.norm().clamp_min(1e-9)
def _non(P, v):
    q = HT[v["quanh"]]; u = u_kh(v); rr = torch.linalg.norm(P - q.t, dim=-1); cosv = ((P - q.t) @ u) / rr.clamp_min(1e-9)
    return (cosv - v["khuyet"]["cos"]) * rr
def phi_vanh(P, v):
    fN = HT[v["quanh"]].phi(P); f = torch.maximum(fN - v["day"], -(fN - v["trong"]))
    return torch.maximum(f, _non(P, v)) if v.get("khuyet") else f
def phi_cam(P, v):
    fN = HT[v["quanh"]].phi(P); f = torch.maximum(fN - v["day"], _non(P, v)) if v.get("khuyet") else fN - v["day"]
    return torch.minimum(f, fN)
def phi_ten(ten, P): return phi_vo(P) if ten == VO_TEN else HT[ten].phi(P)

CHUA = {r["a"]: r for r in QH["chua"]}
DAYNG = {}
for r in QH.get("day_nguoc", []): DAYNG.setdefault(r["a"], []).append(r)
CAMR = {}
for r in QH.get("roi_vung_cam", []):
    for a in r["a"]: CAMR.setdefault(a, []).append(r)
ROI_LON = QH["roi_lon"]; RL = set(ROI_LON["loai"]); S_RL = ROI_LON["s"]
RNL = QH["roi_nho_lon"]; S_NL = RNL["s"]; S_NN = QH["roi_nho_nho"]["s"]
KCL_NN = QH["roi_nho_nho"].get("khe_cung_loai", {}); KCL_LON = QH["roi_lon"].get("khe_cung_loai", {}); KCAP_NN = QH["roi_nho_nho"].get("khe_cap", [])
assert not (set(KCL_NN) & MIEN_NHO) and not (set(KCL_LON) & RL)
S_MAT = None
if KCL_NN or KCAP_NN:
    S_MAT = torch.full((len(NHO), len(NHO)), S_NN, dtype=F64, device=DEV)
    for l, s_ in KCL_NN.items(): S_MAT[LID[l], LID[l]] = s_
    for a_, b_, s_ in KCAP_NN: S_MAT[LID[a_], LID[b_]] = s_; S_MAT[LID[b_], LID[a_]] = s_
S_LOC = S_NN if S_MAT is None else max(S_NN, float(S_MAT.max()))
def khe_nn(lA, lB):
    if S_MAT is None: return S_NN
    return S_MAT[lA, lB][:, None]
_tn = {CHUA[l]["b"] for l in NHO}; _kn = {CHUA[l]["t"] for l in NHO}
if _tn != {VO_TEN} or len(_kn) != 1: raise NotImplementedError("bộ phận nhỏ phải cùng chứa trong vỏ với một khe")
T_NHO = _kn.pop(); _cn = {(r["vung"], r["s"]) for l in NHO for r in CAMR.get(l, [])}
if any(CAMR.get(l) for l in NHO):
    if len(_cn) != 1 or any(len(CAMR.get(l, [])) != 1 for l in NHO): raise NotImplementedError("bộ phận nhỏ phải cùng một vùng cấm")
    CAM_NHO, S_CN = _cn.pop()
else: CAM_NHO = None

def phat_lon(bp, cac, day_du=False):
    """luật cứng cho bộ phận lớn bp với các bộ phận lớn khác trong cac."""
    P = bp.mau(day_du); L = 0.0; r = CHUA[bp.loai]
    L = L + relu2(cong_b(phi_ten(r["b"], P), [r["t"]], r["b"])).mean()
    for d in DAYNG.get(bp.loai, []):
        v = VANH[d["vanh"]]; q = HT[v["quanh"]]
        if v.get("khuyet"):
            u = u_kh(v); cosv = ((P - q.t) @ u) / torch.linalg.norm(P - q.t, dim=1).clamp_min(1e-9)
            L = L + (relu2(cong_b(phi_ten(d["b"], P), d["khe"], d["b"])) * (cosv < v["khuyet"]["cos"])).mean()
        else: L = L + relu2(cong_b(phi_ten(d["b"], P), d["khe"], d["b"])).mean()
    for c in CAMR.get(bp.loai, []): L = L + relu2(c["s"] + DEL - phi_cam(P, VANH[c["vung"]])).mean()
    if bp.loai in RL:
        for q in cac:
            if q is bp or q.loai not in RL: continue
            L = L + relu2(S_RL + DEL - q.phi(P)).mean() + relu2(S_RL + DEL - bp.phi(q.mau(day_du))).mean()
    s_cl = KCL_LON.get(bp.loai)
    if s_cl is not None:
        for q in cac:
            if q is bp or q.loai != bp.loai: continue
            L = L + relu2(s_cl + DEL - q.phi(P)).mean() + relu2(s_cl + DEL - bp.phi(q.mau(day_du))).mean()
    return L
def cap_lan_can(tA, SA, tidA, tB, SB, tidB, du=SH["du_cap_toi_uu"], cung=False, mA=None, mB=None, lA=None, lB=None):
    """cặp (i,j) có cầu bao cách nhau < S+DEL+du; cung=True: A≡B (bỏ i=j). Bỏ cặp cùng loại thuộc 'mien_cung_loai'."""
    with torch.no_grad():
        d = torch.cdist(tA, tB); thr = (RB[tidA] * SA)[:, None] + (RB[tidB] * SB)[None, :] + S_LOC + DEL + du; m = d < thr
        if cung: m.fill_diagonal_(False)
        if mA is not None: m &= ~(mA[:, None] & mB[None, :] & (lA[:, None] == lB[None, :]))
        return m.nonzero(as_tuple=True)
def cap_nhom(nh, du): return cap_lan_can(nh.t, nh.sc(), nh.tid, nh.t, nh.sc(), nh.tid, du=du, cung=True, mA=nh.mien, mB=nh.mien, lA=nh.lid, lB=nh.lid)
def phat_nho(nh, lon, M=SH["M_nho"], cap=None):
    """phạt luật cứng của từng bộ phận nhỏ (n,)."""
    P = nh.diem(M); n, m, _ = P.shape; S = nh.sc(); Rm = nh.R()
    L = relu2(cong_b(phi_vo(P), [T_NHO], VO_TEN)).mean(1)
    if CAM_NHO: L = L + relu2(S_CN + DEL - phi_cam(P, VANH[CAM_NHO])).mean(1)
    for b in lon:
        if b.loai not in RNL["doi_tac"]: continue
        L = L + relu2(S_NL + DEL - b.phi(P)).mean(1)
        Pb = b.mau(); vb = phi_khung(Pb[None].expand(n, -1, -1), nh.t[:, None, :], Rm[:, None], S[:, None], nh.tid[:, None].expand(n, Pb.shape[0]))
        L = L + relu2(S_NL + DEL - vb).mean(1)
    if cap is None: cap = cap_nhom(nh, SH["du_cap_toi_uu"])
    i, j = cap
    if i.numel():
        v = phi_khung(P[i], nh.t[j][:, None, :], Rm[j][:, None], S[j][:, None], nh.tid[j][:, None].expand(-1, m))
        c = relu2(khe_nn(nh.lid[i], nh.lid[j]) + DEL - v).mean(1); L = L.index_add(0, i, c).index_add(0, j, c)
    return L
def phat_mem_lon():
    L = 0.0
    if ARG.khong_mem: return L
    for r in QH.get("lech_tam", []):
        q = HT[r["a"]]; L = L + relu2(q.t.norm() - r["khoang"][1] * r["ban_kinh"]) + relu2(r["khoang"][0] * r["ban_kinh"] - q.t.norm())
    for r in QH.get("gan", []):
        a = HT[r["a"]]
        if r.get("do") == "mau":
            P = a.mau(); d = (-phi_vo(P)).min() if r["b"] == VO_TEN else HT[r["b"]].phi(P).min(); L = L + relu2(d - r["toi_da"])
        elif r["b"] == VO_TEN: L = L + relu2(-phi_vo(a.t[None])[0] - a.sc() - r["toi_da"])
        else: L = L + relu2(HT[r["b"]].phi(a.t[None])[0] - a.sc() - r["toi_da"])
    return L
def kc_gan_mau(nh, r):
    """khoảng mặt–mặt (đo bằng điểm mẫu) từ mỗi bộ phận nhỏ tới đích gần nhất trong danh sách r['b'] (vỏ, bộ phận lớn, hoặc loại nhỏ)."""
    P = nh.diem(SH["M_gan"]); n, m, _ = P.shape; out = torch.full((n,), 1e9, dtype=F64, device=DEV)
    for b in (r["b"] if isinstance(r["b"], list) else [r["b"]]):
        if b == VO_TEN: d = (-phi_vo(P)).min(1).values
        elif b in HT: d = HT[b].phi(P).min(1).values
        else:
            jd = nh.mask(b).nonzero().flatten()
            if not jd.numel(): continue
            kk = min(SH["gan_so_lang_gieng"], jd.numel())
            with torch.no_grad(): nb = jd[torch.cdist(nh.t, nh.t[jd]).topk(kk, largest=False).indices]
            S = nh.sc(); Rm = nh.R()
            v = phi_khung(P[:, None].expand(n, kk, m, 3), nh.t[nb][:, :, None], Rm[nb][:, :, None], S[nb][:, :, None], nh.tid[nb][:, :, None].expand(n, kk, m))
            d = v.min(2).values.min(1).values
        out = torch.minimum(out, d)
    return out
def phap_vo(P):
    with torch.no_grad():
        e = SH["eps_phap"]; I3 = torch.eye(3, dtype=F64, device=DEV)
        g_ = torch.stack([phi_vo(P + e * I3[k]) - phi_vo(P - e * I3[k]) for k in range(3)], -1); return g_ / g_.norm(dim=-1, keepdim=True).clamp_min(1e-12)
def phat_mem_nho(nh, L):
    if ARG.khong_mem: return L
    for r in QH.get("gan_nho", []):
        if r.get("do") == "mau": L = L + r["trong_so"] * (relu2(kc_gan_mau(nh, r) - r["toi_da"]) * nh.mask(r["a"])).sum()
        else: L = L + r["trong_so"] * (relu2(HT[r["b"]].phi(nh.t) - r["toi_da"]) * nh.mask(r["a"])).sum()
    for r in QH.get("huong", []):
        cs = (nh.R()[:, :, 2] * phap_vo(nh.t)).sum(-1).abs()
        L = L + r["trong_so"] * (relu2(math.cos(math.radians(r["lech_toi_da_do"])) - cs) * nh.mask(r["a"])).sum()
    return L

_h = SH["hat_so_luong"]; _rd = np.random.default_rng(ARG.seed * _h[0] + _h[1])
SO = {l: 1 for l in LON + NHO}
for l in DT["so_luong_thu_tu"]: SO[l] = int(_rd.integers(*LOAI[l]["so_luong"]["ngau_nhien"]))
for l in LON + NHO:
    if isinstance(LOAI[l]["so_luong"], int): SO[l] = LOAI[l]["so_luong"]
for x in ARG.so:
    l, n_ = x.split("="); n_ = int(n_)
    if l not in SO: raise KeyError(l)
    if n_: SO[l] = n_
DON = [l for l in LON if LOAI[l]["dat_ten"] == "don"]
for l in DON: assert SO[l] == 1, l
def ten_ban(l, i): return l if LOAI[l]["dat_ten"] == "don" else f"{l}_{i + 1}"
def sig(lo, hi, z): return lo + (hi - lo) / (1 + math.exp(-z))
def trong_deu(g, sau):
    while True:
        p = g.uniform(-R_VO, R_VO, 3)
        if np.linalg.norm(p) < R_VO - sau: return p
def trong_vung(g, sau, c, nguong):
    for _ in range(SH["trong_vung_so_lan"]):
        p = g.uniform(-R_VO, R_VO, 3)
        if np.linalg.norm(p) < R_VO - sau and np.linalg.norm(p - c) > nguong: return p
    return trong_deu(g, sau)
class MoiTruong(dict):
    def __missing__(s, k):
        if k.startswith("t."): return HT[k[2:]].t.cpu().numpy()
        if k.startswith("s."): return float(HT[k[2:]].sc())
        raise KeyError(k)
def tinh(e, env, g):
    if isinstance(e, (int, float)): return e
    if isinstance(e, str):
        if e.startswith("$"): return env[e[1:]]
        raise ValueError(e)
    op, a = e[0], e[1:]
    if op == "vec": return np.array(a[0] if len(a) == 1 else a, float)
    if op == "bang": return np.array(a[0][int(tinh(a[1], env, g))], float)
    if op == "huong": v = g.normal(size=3); return v / np.linalg.norm(v)
    if op == "deu": lo = tinh(a[0], env, g); hi = tinh(a[1], env, g); return g.uniform(lo, hi)
    if op == "chuan": return g.normal(size=int(a[0]))
    if op == "trong_deu": return trong_deu(g, tinh(a[0], env, g))
    if op == "trong_vung": s_ = tinh(a[0], env, g); c = tinh(a[1], env, g); return trong_vung(g, s_, c, tinh(a[2], env, g))
    if op in ("+", "*", "-", "/"):
        x = tinh(a[0], env, g)
        for b in a[1:]:
            y = tinh(b, env, g); x = x + y if op == "+" else x * y if op == "*" else x - y if op == "-" else x / y
        return x
    if op == "cheo": x = tinh(a[0], env, g); return np.cross(x, tinh(a[1], env, g))
    if op == "chuan_hoa": v = tinh(a[0], env, g); return v / np.linalg.norm(v)
    if op == "noi": x = tinh(a[0], env, g); return np.concatenate([x, tinh(a[1], env, g)])
    if op == "tia":
        o = tinh(a[1], env, g); u = tinh(a[2], env, g); lo, hi, n = a[3]
        with torch.no_grad():
            ds_ = torch.linspace(lo, hi, n, device=DEV, dtype=F64); f = HT[a[0]].phi(T(o)[None] + ds_[:, None] * T(u)[None])
            return float(ds_[(f >= tinh(a[4], env, g)).nonzero()[0, 0]])
    raise ValueError(op)
KT = DT["khoi_tao"]; DL = DT.get("dat_lai", {}).get(ARG.khoi_tao, {})
DO_THI = KT.get("kieu") == "do_thi"
def gan_vo(l):
    for r in QH.get("gan", []) + QH.get("gan_nho", []):
        bs = r["b"] if isinstance(r["b"], list) else [r["b"]]
        if r["a"] == l and VO_TEN in bs: return r["toi_da"]
    return None
HUONG_L = {r["a"] for r in QH.get("huong", [])}
def diem_tren_vo(g, lui, le):
    b = NC - lui; core = b - (RBO - lui); ar = np.array([core[1] * core[2], core[0] * core[2], core[0] * core[1]]); ar = ar / ar.sum()
    for _ in range(SH["tren_vo_so_lan"]):
        ax = int(g.choice(3, p=ar)); sg = 1.0 if g.uniform() < 0.5 else -1.0; u = g.uniform(-core, core)
        if any(abs(u[j]) > core[j] - le for j in range(3) if j != ax): continue
        p = u.copy(); p[ax] = sg * b[ax]; n_ = np.zeros(3); n_[ax] = -sg; return p, n_
    raise RuntimeError("không tìm được điểm trên vỏ")
def trong_hop(g, sau):
    for _ in range(SH["trong_vung_so_lan"]):
        p = g.uniform(-(NC - D_VO), NC - D_VO)
        if float(phi_vo(T(p)[None])[0]) < -sau: return p
    raise RuntimeError("không tìm được điểm trong vỏ")
def dat_mot(l, kh, s, g):
    """tư thế ban đầu cho một bộ phận loại l (khuôn kh, cỡ s) theo đồ thị luật: gần vỏ thì đặt trên mặt trong vỏ ở độ sâu ngẫu nhiên
    trong [khe chứa + β, khoảng tối đa], xoay theo luật hướng; còn lại rải đều trong vỏ."""
    gv = gan_vo(l); KH_ = KHUON[kh]; tc = CHUA[l]["t"]
    if ARG.khoi_tao == "deu" or gv is None:
        p = trong_hop(g, KH_["rb"] * s + tc + DEL_M); return p, g.normal(size=6)
    nua = (KH_["nz"] if l in HUONG_L else KH_["rb"]) * s
    p, n_ = diem_tren_vo(g, D_VO, KH_["rb"] * s); gap = g.uniform(tc + DEL_M, max(gv, tc + DEL_M)); p = p + n_ * (nua + gap)
    if l in HUONG_L:
        e1 = np.cross(n_, g.normal(size=3)); e1 /= np.linalg.norm(e1); e2 = np.cross(n_, e1); return p, np.concatenate([e1, e2])
    return p, g.normal(size=6)
def khoi_tao_do_thi(k):
    gk = np.random.default_rng([ARG.seed, SH["hat_khuon_co"], k]); KC, Z = {}, {}
    for l in LON + NHO:
        KC[l] = [KHL[l][int(gk.integers(len(KHL[l])))] for _ in range(SO[l])]; Z[l] = [gk.normal() for _ in range(SO[l])]
    g = np.random.default_rng([ARG.seed, SH["hat_khoi_tao"], k] + ([SH["hat_deu_cong"]] if ARG.khoi_tao == "deu" else [])); TT = {}
    for l in LON + NHO:
        TT[l] = [dat_mot(l, KC[l][i], sig(*DAI_K.get(KC[l][i], DAI[l]), Z[l][i]), g) for i in range(SO[l])]
    lon = [BoPhan(ten_ban(l, i), l, KC[l][i], TT[l][i][0], TT[l][i][1], Z[l][i]) for l in LON for i in range(SO[l])]
    ds = [(ten_ban(l, i), l, KC[l][i], TT[l][i][0], TT[l][i][1], Z[l][i]) for l in NHO for i in range(SO[l])]
    return lon, Nhom(ds), g
def khoi_tao(k):
    if DO_THI: return khoi_tao_do_thi(k)
    gk = np.random.default_rng([ARG.seed, SH["hat_khuon_co"], k])
    KC = {l: [KHL[l][0]] * SO[l] for l in LON + NHO if len(KHL[l]) == 1 and LOAI[l]["khuon"] != "thu_vien"}; Z = {}
    for st in KT["gk"]:
        if "khuon" in st: l = st["khuon"]; KC[l] = [KHL[l][int(gk.integers(len(KHL[l])))] for _ in range(SO[l])]
        elif "z" in st: l = st["z"]; Z[l] = [gk.normal() for _ in range(SO[l])]
        else:
            ds = st["z_gop"]; zz = gk.normal(size=sum(SO[l] for l in ds)); o = 0
            for l in ds: Z[l] = list(zz[o:o + SO[l]]); o += SO[l]
    g = np.random.default_rng([ARG.seed, SH["hat_khoi_tao"], k] + ([SH["hat_deu_cong"]] if ARG.khoi_tao == "deu" else []))
    env = MoiTruong(DEL=DEL); TT = {}
    for l in DON: env["s0." + l] = sig(*DAI[l], Z[l][0])
    for st in KT[ARG.khoi_tao]:
        if "gan" in st: env[st["gan"]] = tinh(st["bt"], env, g); continue
        l = st["dat"]; TT[l] = []
        for i in range(SO[l]): env["i"] = i; t = tinh(st["t"], env, g); TT[l].append((t, tinh(st["r"], env, g)))
    lon = [BoPhan(ten_ban(l, i), l, KC[l][i], TT[l][i][0], TT[l][i][1], Z[l][i]) for l in LON for i in range(SO[l])]
    ds = [(ten_ban(l, i), l, KC[l][i], TT[l][i][0], TT[l][i][1], Z[l][i]) for l in NHO for i in range(SO[l])]
    return lon, Nhom(ds), g
def mau_tu_the(nh, idx, g):
    """tư thế ngẫu nhiên mới cho các bộ phận nhỏ idx (cùng phân bố khởi tạo của nhánh); cỡ GIỮ NGUYÊN."""
    env = MoiTruong(DEL=DEL); t, r = [], []
    if DO_THI: S_ = nh.sc()
    for i in idx:
        if DO_THI: p_, r_ = dat_mot(nh.loai[i], TEN_K[int(nh.tid[i])], float(S_[i]), g); t.append(p_); r.append(r_); continue
        e = DL[nh.loai[i]]; t.append(tinh(e["t"], env, g)); r.append(tinh(e["r"], env, g))
    return T(np.array(t)).reshape(-1, 3), T(np.array(r)).reshape(-1, 6), nh.z[torch.tensor(idx, device=DEV)].clone()
def lay_mau_lon(b, g):
    if DO_THI: p_, r_ = dat_mot(b.loai, b.khuon, float(b.sc()), g); b.t = T(p_); b.r = T(r_); return
    env = MoiTruong(DEL=DEL); e = DL[b.loai]; t = tinh(e["t"], env, g)
    b.t = T(t); b.r = T(tinh(e["r"], env, g))
def dat_co_dinh(lon):
    env = MoiTruong(DEL=DEL)
    for st in DT["tu_the_co_dinh"]:
        if "gan" in st: env[st["gan"]] = tinh(st["bt"], env, None); continue
        for i, b in enumerate([b for b in lon if b.loai == st["dat"]]):
            env["i"] = i
            if "t" in st: b.t = T(tinh(st["t"], env, None))
            if "r" in st: b.r = T(tinh(st["r"], env, None))

LR_PLACE = SH["adam_lr"] * min(1.0, h1 / 0.1796875)
V121["adam_lr"] = LR_PLACE
def toi_uu(ps, ham, buoc, lr=LR_PLACE):
    for p in ps: p.requires_grad_(True)
    opt = torch.optim.Adam(ps, lr=lr)
    for _ in range(buoc): opt.zero_grad(); L = ham(); L.backward(); opt.step()
    for p in ps: p.requires_grad_(False)
def phat_day_du(lon, nh):
    with torch.no_grad():
        pl = {b.ten: float(phat_lon(b, lon, day_du=True)) for b in lon}
        torch.manual_seed(SH["hat_day_du"]); pn = phat_nho(nh, lon, M=SH["M_day_du"], cap=cap_nhom(nh, SH["du_cap_kiem"]))
    return pl, {nh.ten[i]: float(pn[i]) for i in range(nh.n)}
NG = SH["nguong_phat"]
def thu_loai_nho(nh, lon, g, bao):
    """đặt tuần tự từng bộ phận nhỏ: ≤ lo·C ứng viên, chỉ kiểm với bộ phận đã đặt; không được thì giữ ứng viên phạt nhỏ nhất."""
    da = torch.zeros(nh.n, dtype=torch.bool, device=DEV); DTL = [b for l in RNL["doi_tac_thu_loai"] for b in lon if b.loai == l]
    C = SH["thu_loai_C"]; MM = SH["thu_loai_M"]; NL = SH["thu_loai_so_lo"]
    for i in g.permutation(nh.n).tolist():
        tot = (1e18, None)
        for lo in range(NL):
            tc, rc, zc = mau_tu_the(nh, [i] * C, g); tid = nh.tid[i].expand(C); S = DAI_T(tid, zc); Rc = rot6(rc)
            with torch.no_grad():
                P = nh.diem(MM, t=tc, r=rc, z=zc, tid=tid); m = P.shape[1]
                L = relu2(cong_b(phi_vo(P), [T_NHO], VO_TEN)).mean(1)
                if CAM_NHO: L = L + relu2(S_CN + DEL - phi_cam(P, VANH[CAM_NHO])).mean(1)
                for b in DTL:
                    L = L + relu2(S_NL + DEL - b.phi(P)).mean(1)
                    Pb = b.mau(True); L = L + relu2(S_NL + DEL - phi_khung(Pb[None].expand(C, -1, -1), tc[:, None], Rc[:, None], S[:, None], tid[:, None].expand(C, Pb.shape[0]))).mean(1)
                pj = da.nonzero().flatten()
                if pj.numel():
                    Sj = nh.sc()[pj]; ci, cj = cap_lan_can(tc, S, tid, nh.t[pj], Sj, nh.tid[pj], du=SH["du_cap_kiem"], mA=nh.mien[i].expand(C), mB=nh.mien[pj], lA=nh.lid[i].expand(C), lB=nh.lid[pj])
                    if ci.numel():
                        q = pj[cj]; Rj = nh.R()[q]
                        v = phi_khung(P[ci], nh.t[q][:, None], Rj[:, None], nh.sc()[q][:, None], nh.tid[q][:, None].expand(-1, m))
                        kc = khe_nn(nh.lid[i].expand(ci.numel()), nh.lid[q]); L = L.index_add(0, ci, relu2(kc + DEL - v).mean(1))
                        Pq = nh.diem(MM, t=nh.t[q], r=nh.r[q], z=nh.z[q], tid=nh.tid[q])
                        v2 = phi_khung(Pq, tc[ci][:, None], Rc[ci][:, None], S[ci][:, None], tid[ci][:, None].expand(-1, Pq.shape[1]))
                        L = L.index_add(0, ci, relu2(kc + DEL - v2).mean(1))
            ok = (L <= NG).nonzero().flatten(); DEM["danh_gia"] += C
            if ok.numel(): c = int(ok[0]); bao["dat_lai"] += lo * C + c + 1; tot = (0.0, c); break
            c = int(torch.argmin(L))
            if float(L[c]) < tot[0]: tot = (float(L[c]), c, tc[c], rc[c], zc[c])
        else: bao["dat_lai"] += NL * C; bao["khong_dat"] += 1
        if tot[0] == 0.0: nh.t[i], nh.r[i], nh.z[i] = tc[tot[1]], rc[tot[1]], zc[tot[1]]
        else: nh.t[i], nh.r[i], nh.z[i] = tot[2], tot[3], tot[4]
        da[i] = True

DEM = {"danh_gia": 0}
def _xoay(Rm, goc, gen):
    ax = torch.randn(Rm.shape[:-2] + (3,), generator=gen, device=DEV, dtype=F64); ax = ax / ax.norm(dim=-1, keepdim=True)
    Kx = torch.zeros(Rm.shape[:-2] + (3, 3), device=DEV, dtype=F64)
    Kx[..., 0, 1], Kx[..., 0, 2], Kx[..., 1, 2] = -ax[..., 2], ax[..., 1], -ax[..., 0]; Kx = Kx - Kx.transpose(-1, -2)
    th = goc[..., None, None]; Q = torch.eye(3, device=DEV, dtype=F64) + torch.sin(th) * Kx + (1 - torch.cos(th)) * (Kx @ Kx)
    return Q @ Rm
def _r6(Rm): return torch.cat([Rm[..., :, 0], Rm[..., :, 1]], -1)
def _T0(dE, p0):
    up = dE[dE > 0]
    return float(up.mean()) / math.log(1 / p0) if up.numel() else ARG.un_T1
def u_nhiet(lon, nh, k, ngan, bc):
    """Ủ nhiệt trên J(g) = Λ·Σ phạt cứng + S(g): cùng hàm phạt, cùng β, cùng khuôn/cỡ/khởi tạo; chỉ dời + xoay (cỡ cố định).
    Hai chặng như tu_dong: (1) bộ phận lớn (tính lại toàn phần, ít bộ phận), (2) bộ phận nhỏ giữ cố định bộ phận lớn,
    mỗi vòng đề xuất SONG SONG cho một tập độc lập (không cặp nào có thể chạm nhau trước/sau bước) ⇒ tương đương tuần tự.
    Bộ mẫu cố định trong mỗi đoạn un_lam_moi vòng (năng lượng tất định ⇒ Metropolis đúng), rút lại định kỳ."""
    gen = torch.Generator(device=DEV); gen.manual_seed(ARG.seed * 7919 + k)
    t_bd = time.time(); t_lon = t_bd + ARG.un_ti_le_lon * ngan; t_het = t_bd + ngan
    def lich(u, ab): return ab[0] * (ab[1] / ab[0]) ** min(max(u, 0.0), 1.0)
    def gan_mau():
        for b in lon:
            Sv = float(b.sc()); n = min(NM, max(SH["mau_toi_thieu"], int(KHUON[b.khuon]["dien_tich"] * Sv * Sv / DEL ** 2)))
            b.idx_co_dinh = torch.randint(0, NM, (n,), device=DEV, generator=gen)
    def E_lon():
        DEM["danh_gia"] += len(lon); return float(sum(phat_lon(b, lon) for b in lon) * W + phat_mem_lon())
    def de_xuat(b, u):
        if float(torch.rand((), generator=gen, device=DEV)) < ARG.un_p_doi:
            d = torch.randn(3, generator=gen, device=DEV, dtype=F64) * lich(u, ARG.un_buoc_t) * b.rb(); return b.t + d, b.r
        g_ = torch.randn((), generator=gen, device=DEV, dtype=F64) * lich(u, ARG.un_buoc_r); return b.t, _r6(_xoay(b.R(), g_, gen))
    s1 = {"de_xuat": 0, "nhan": 0}
    with torch.no_grad():
        gan_mau(); E = E_lon(); E_dau1 = E; dE = []
        for _ in range(30):
            b = lon[int(torch.randint(len(lon), (1,), generator=gen, device=DEV))]; t0_, r0_ = b.t, b.r; b.t, b.r = de_xuat(b, 0.0); dE.append(E_lon() - E); b.t, b.r = t0_, r0_
        T0 = max(_T0(torch.tensor(dE), ARG.un_p0), ARG.un_T1); s1["T0"] = T0; vong = 0
        while (vong < ARG.un_vong_lon) if BUOC_UN else (time.time() < t_lon):
            u = vong / ARG.un_vong_lon if BUOC_UN else (time.time() - t_bd) / max(t_lon - t_bd, 1e-9); Tn = T0 * (ARG.un_T1 / T0) ** u
            if vong and vong % ARG.un_lam_moi == 0: gan_mau(); E = E_lon()
            b = lon[int(torch.randint(len(lon), (1,), generator=gen, device=DEV))]; t0_, r0_ = b.t, b.r; b.t, b.r = de_xuat(b, u)
            En = E_lon(); s1["de_xuat"] += 1; vong += 1
            if En <= E or float(torch.rand((), generator=gen, device=DEV, dtype=F64)) < math.exp(-(En - E) / Tn): E = En; s1["nhan"] += 1
            else: b.t, b.r = t0_, r0_
        s1["J_dau"] = E_dau1; s1["J_cuoi_mau_co_dinh"] = E; s1["vong"] = vong; s1["giay"] = round(time.time() - t_bd, 3)
        for b in lon: b.idx_co_dinh = None; b.r = _r6(b.R())
    n = nh.n; S = nh.sc(); tid = nh.tid; rbS = RB[tid] * S; lid = nh.lid; mien = nh.mien; MN = SH["M_nho"]; ntl = SH["mau_so_trong_lo"]
    DTN = [b for b in lon if b.loai in RNL["doi_tac"]]; s2 = {"de_xuat": 0, "nhan": 0, "vong": 0, "co_tap": 0}
    for r in QH.get("gan_nho", []):
        if r.get("do") == "mau" and any(b not in HT and b != VO_TEN for b in (r["b"] if isinstance(r["b"], list) else [r["b"]])): raise NotImplementedError("u_nhiet: gan_nho đo bằng mẫu tới loại nhỏ")
    MK = {l: nh.mask(l) for l in set(r["a"] for r in QH.get("gan_nho", []) + QH.get("huong", []))}
    du = SH["du_cap_toi_uu"]
    with torch.no_grad():
        tc = nh.t.clone(); Rc = nh.R().clone()
        def mau():
            ii = torch.randint(0, NM, (n, MN), device=DEV, generator=gen)
            Ls = torch.cat([MATS[tid[:, None], ii], TRONGS[tid][:, :ntl]], 1) * S[:, None, None]
            gan_mau(); Pb = [b.mau() for b in DTN]
            for b in lon: b.idx_co_dinh = None
            return Ls, Pb
        def diem(t, Rm, L): return torch.einsum("nij,nmj->nmi", Rm, L) + t[:, None, :]
        def E_mot(a, t, Rm, P):
            """phần J do các bộ phận nhỏ a (tư thế t, Rm, điểm P) gánh: Λ·(phạt riêng + 2·Σ_j (c_aj + c_ja)) + phạt mềm riêng.
            Tổng J nhóm nhỏ = Λ·(Σ riêng + ½ Σ_a cặp_a) + Σ mềm — khớp phat_nho (index_add hai chiều trên cặp có thứ tự)."""
            DEM["danh_gia"] += a.numel(); na, m, _ = P.shape; Sa = S[a]; ta_ = tid[a]
            L = relu2(cong_b(phi_vo(P), [T_NHO], VO_TEN)).mean(1)
            if CAM_NHO: L = L + relu2(S_CN + DEL - phi_cam(P, VANH[CAM_NHO])).mean(1)
            for b, Pbb in zip(DTN, Pb):
                L = L + relu2(S_NL + DEL - b.phi(P)).mean(1)
                vb = phi_khung(Pbb[None].expand(na, -1, -1), t[:, None, :], Rm[:, None], Sa[:, None], ta_[:, None].expand(na, Pbb.shape[0]))
                L = L + relu2(S_NL + DEL - vb).mean(1)
            ci, cj = cap_lan_can(t, Sa, ta_, tc, S, tid, du=du, mA=mien[a], mB=mien, lA=lid[a], lB=lid)
            k_ = cj != a[ci]; ci, cj = ci[k_], cj[k_]
            cp = torch.zeros(na, dtype=F64, device=DEV)
            if ci.numel():
                kc = khe_nn(lid[a][ci], lid[cj])
                v1 = phi_khung(P[ci], tc[cj][:, None, :], Rc[cj][:, None], S[cj][:, None], tid[cj][:, None].expand(-1, m))
                v2 = phi_khung(Pc[cj], t[ci][:, None, :], Rm[ci][:, None], Sa[ci][:, None], ta_[ci][:, None].expand(-1, Pc.shape[1]))
                cp = cp.index_add(0, ci, 2 * (relu2(kc + DEL - v1).mean(1) + relu2(kc + DEL - v2).mean(1)))
            mem = torch.zeros(na, dtype=F64, device=DEV)
            for r in QH.get("gan_nho", []):
                if r.get("do") == "mau":
                    Pm = torch.cat([P[:, :SH["M_gan"]], P[:, MN:]], 1); dmin = torch.full((na,), 1e9, dtype=F64, device=DEV)
                    for b in (r["b"] if isinstance(r["b"], list) else [r["b"]]):
                        dmin = torch.minimum(dmin, (-phi_vo(Pm)).min(1).values if b == VO_TEN else HT[b].phi(Pm).min(1).values)
                    mem = mem + r["trong_so"] * relu2(dmin - r["toi_da"]) * MK[r["a"]][a]
                else: mem = mem + r["trong_so"] * relu2(HT[r["b"]].phi(t) - r["toi_da"]) * MK[r["a"]][a]
            for r in QH.get("huong", []):
                cs = (Rm[:, :, 2] * phap_vo(t)).sum(-1).abs(); mem = mem + r["trong_so"] * relu2(math.cos(math.radians(r["lech_toi_da_do"])) - cs) * MK[r["a"]][a]
            if ARG.khong_mem: mem = mem * 0
            return W * (L + cp) + mem, W * (L + 0.5 * cp) + mem
        def tap_doc_lap(mt):
            """tập độc lập (Luby, 3 lượt): i, j trong tập ⇒ cầu bao cách nhau ≥ ngưỡng cặp + bước dời tối đa của cả hai."""
            thr = rbS[:, None] + rbS[None, :] + S_LOC + DEL + du + mt[:, None] + mt[None, :]
            C = torch.cdist(tc, tc) < thr; C.fill_diagonal_(False)
            pr = torch.rand(n, generator=gen, device=DEV, dtype=F64); con = torch.ones(n, dtype=torch.bool, device=DEV); sel = torch.zeros_like(con)
            for _ in range(3):
                pm = torch.where(con, pr, torch.full_like(pr, -1.0)); nbm = torch.where(C, pm[None, :], torch.full_like(C, -1.0, dtype=F64)).max(1).values
                win = con & (pm > nbm); sel |= win; con &= ~win & ~(C & win[None, :]).any(1)
            return sel.nonzero().flatten()
        def de_xuat_lo(a, u):
            na = a.numel(); doi = torch.rand(na, generator=gen, device=DEV, dtype=F64) < ARG.un_p_doi
            sg = lich(u, ARG.un_buoc_t) * rbS[a]; d = torch.randn(na, 3, generator=gen, device=DEV, dtype=F64) * sg[:, None]
            nd = d.norm(dim=1, keepdim=True); d = d * torch.clamp(3 * sg[:, None] / nd.clamp_min(1e-12), max=1.0)
            g_ = torch.randn(na, generator=gen, device=DEV, dtype=F64) * lich(u, ARG.un_buoc_r)
            tn = tc[a] + d * doi[:, None]; Rn = torch.where(doi[:, None, None], Rc[a], _xoay(Rc[a], g_, gen))
            return tn, Rn
        Ls, Pb = mau(); Pc = diem(tc, Rc, Ls)
        ia = torch.arange(n, device=DEV); _, eh = E_mot(ia, tc, Rc, Pc); s2["J_dau"] = float(eh.sum())
        dEs = []; u0 = 0.0
        while sum(x.numel() for x in dEs) < 50 and len(dEs) < 50:
            a = tap_doc_lap(3 * lich(u0, ARG.un_buoc_t) * rbS); tn, Rn = de_xuat_lo(a, u0)
            Eo, _ = E_mot(a, tc[a], Rc[a], Pc[a]); En, _ = E_mot(a, tn, Rn, diem(tn, Rn, Ls[a])); dEs.append(En - Eo)
        T0 = max(_T0(torch.cat(dEs), ARG.un_p0), ARG.un_T1); s2["T0"] = T0; t2 = time.time(); vong = 0
        while (vong < ARG.un_vong_nho) if BUOC_UN else (time.time() < t_het):
            u = vong / ARG.un_vong_nho if BUOC_UN else (time.time() - t2) / max(t_het - t2, 1e-9); Tn = T0 * (ARG.un_T1 / T0) ** min(u, 1.0)
            if vong and vong % ARG.un_lam_moi == 0: Ls, Pb = mau(); Pc = diem(tc, Rc, Ls)
            a = tap_doc_lap(3 * lich(u, ARG.un_buoc_t) * rbS); tn, Rn = de_xuat_lo(a, u); Pn = diem(tn, Rn, Ls[a])
            Eo, _ = E_mot(a, tc[a], Rc[a], Pc[a]); En, _ = E_mot(a, tn, Rn, Pn); d_ = En - Eo
            ok = (d_ <= 0) | (torch.rand(a.numel(), generator=gen, device=DEV, dtype=F64) < torch.exp(-d_.clamp_min(0) / Tn))
            ac = a[ok]; tc[ac] = tn[ok]; Rc[ac] = Rn[ok]; Pc[ac] = Pn[ok]
            s2["de_xuat"] += a.numel(); s2["nhan"] += int(ok.sum()); s2["co_tap"] += a.numel(); vong += 1
        s2["vong"] = vong; s2["giay"] = round(time.time() - t2, 3); s2["co_tap_tb"] = s2["co_tap"] / max(vong, 1); del s2["co_tap"]
        _, eh = E_mot(ia, tc, Rc, Pc); s2["J_cuoi_mau_co_dinh"] = float(eh.sum())
        nh.t = tc.clone(); nh.r = _r6(Rc)
    bc.setdefault("u_nhiet", []).append({"k": k, "ngan_sach": round(ngan, 2), "giay": round(time.time() - t_bd, 2), "lon": s1, "nho": s2})

NS = None
if DT.get("ngan_sach"):
    NS_ = DT["ngan_sach"]
    def V_hop(lui):
        a_ = NC - lui - (RBO - lui); rr = RBO - lui
        return 8 * a_.prod() + 8 * (a_[0] * a_[1] + a_[1] * a_[2] + a_[0] * a_[2]) * rr + 2 * math.pi * a_.sum() * rr ** 2 + 4 / 3 * math.pi * rr ** 3
    def V_loi(l):
        out = []
        for kh in KHL[l]:
            s = sum(DAI_K.get(kh, DAI[l])) / 2; KH_ = KHUON[kh]; gv = gan_vo(l); tc = CHUA[l]["t"]; nua = (KH_["nz"] if l in HUONG_L else KH_["rb"]) * s
            zc = nua + ((tc + DEL_M + max(gv, tc + DEL_M)) / 2 if gv is not None else NS_["cach_vo"] + KH_["rb"] * s)
            nn = NS_["no"] + NS_.get("no_them", {}).get(l, 0.0); Lb = KH_["rb"] * s + nn; hs = NS_["h"]
            ax_ = torch.arange(-Lb, Lb + hs, hs, device=DEV, dtype=F64); Q = torch.stack(torch.meshgrid(ax_, ax_, ax_, indexing="ij"), -1).reshape(-1, 3)
            ps = s * eval_at(Q / s, torch.full((Q.shape[0],), TID[kh], device=DEV))
            out.append(float(((ps < nn) & (Q[:, 2] + zc > NS_["cach_vo"])).sum()) * hs ** 3)
        return sum(out) / len(out)
    Vin = V_hop(D_VO); Vc = V_hop(D_VO + NS_["cach_vo"]); VL = {l: V_loi(l) for l in LON + NHO}
    def uoc(): return (Vc - sum(VL[l] * SO[l] for l in LON + NHO)) / Vin
    NS = {"uoc_dau": uoc(), "giam": {}}
    for l in NS_["giam"]:
        while uoc() < NS_["nguong"] and SO[l] > LOAI[l]["so_luong"]["ngau_nhien"][0]: SO[l] -= 1; NS["giam"][l] = NS["giam"].get(l, 0) + 1
    NS["uoc_sau"] = uoc(); NS["V_trong_vo"] = Vin; NS["V_loi"] = VL; print("NGAN_SACH", NS, flush=True)
bao = {"so_luong": [SO[l] for l in DT["so_luong_thu_tu"]], "dat_lai": 0, "khong_dat": 0}; tot = None; W = SH["trong_so_cung"]
if NS: bao["ngan_sach"] = NS
def _lam_nong():
    with torch.enable_grad():
        q_ = torch.linspace(-1, 1, 768, device=DEV, dtype=F64).reshape(256, 3).requires_grad_(True)
        eval_at(q_, torch.zeros(256, dtype=torch.long, device=DEV)).sum().backward()
    torch.cdist(q_.detach(), q_.detach()); torch.cuda.synchronize()
_lam_nong()
def _lam_nong_nhanh():
    global HT
    l_, n_, _ = khoi_tao(0); HT = {b.loai: b for b in l_ if b.loai in DON}
    with torch.no_grad(): phat_nho(n_, l_, M=SH["M_kiem_vong"], cap=cap_nhom(n_, SH["du_cap_kiem"])); [phat_lon(b, l_, day_du=True) for b in l_]
    if ARG.cach == "tu_dong":
        toi_uu([n_.t, n_.r], lambda: phat_mem_nho(n_, phat_nho(n_, l_).sum() * W), 2); toi_uu([p for b in l_ for p in (b.t, b.r)], lambda: sum(phat_lon(b, l_) for b in l_) * W + phat_mem_lon(), 2)
    if ARG.cach == "u_nhiet":
        g_ = torch.Generator(device=DEV); g_.manual_seed(1); _xoay(torch.eye(3, dtype=F64, device=DEV)[None].expand(4, 3, 3), torch.zeros(4, dtype=F64, device=DEV), g_)
    torch.cuda.synchronize(); torch.manual_seed(ARG.seed)
_lam_nong_nhanh(); torch.cuda.synchronize()
if ARG.eval_dung_tu_dau:
    assert not ARG.tu_the_tu; _xay_xa(); E3["bat"] = True; torch.cuda.synchronize()
t_dat0 = time.time(); t_day_du = 0.0
for k in range(0 if ARG.tu_the_tu else ARG.K):
    lon, nh, g = khoi_tao(k); HT = {b.loai: b for b in lon if b.loai in DON}
    if ARG.cach == "tay": dat_co_dinh(lon)
    if ARG.cach == "tu_dong":
        toi_uu([p for b in lon for p in (b.t, b.r)], lambda: sum(phat_lon(b, lon) for b in lon) * W + phat_mem_lon(), ARG.buoc)
        DEM["danh_gia"] += ARG.buoc * len(lon)
        for vong in range(SH["so_vong_nho"]):
            st = {"cap": None, "i": 0}
            def ham():
                if st["i"] % SH["cap_lam_moi"] == 0: st["cap"] = cap_nhom(nh, SH["du_cap_toi_uu"])
                st["i"] += 1; return phat_mem_nho(nh, phat_nho(nh, lon, cap=st["cap"]).sum() * W)
            toi_uu([nh.t, nh.r], ham, ARG.buoc)
            DEM["danh_gia"] += ARG.buoc * nh.n + nh.n
            with torch.no_grad():
                torch.manual_seed(SH["hat_kiem_vong"] + vong); pn = phat_nho(nh, lon, M=SH["M_kiem_vong"], cap=cap_nhom(nh, SH["du_cap_kiem"]))
            xau = (pn > NG).nonzero().flatten().tolist()
            if os.environ.get("GO_LOI"): print("vong", vong, "xau", len(xau), round(time.time() - t00), flush=True)
            if not xau or vong == SH["so_vong_nho"] - 1: break
            t_, r_, z_ = mau_tu_the(nh, xau, g); nh.t[xau], nh.r[xau], nh.z[xau] = t_, r_, z_; bao["dat_lai"] += len(xau)
    elif ARG.cach == "u_nhiet":
        if BUOC_UN: u_nhiet(lon, nh, k, 0.0, bao)
        else:
            con_lai = ARG.ngan_sach_giay - (time.time() - t_dat0) - t_day_du * (ARG.K - k)
            u_nhiet(lon, nh, k, max(con_lai / (ARG.K - k), 0.5), bao)
    else:
        if ARG.cach == "ngau_nhien":
            for lam in range(SH["lam_lai_nhom"]):
                da = []
                for b in lon:
                    for _ in range(SH["thu_lon_so_lan"]):
                        with torch.no_grad():
                            DEM["danh_gia"] += 1
                            if float(phat_lon(b, da + [b], day_du=True)) <= NG: break
                        lay_mau_lon(b, g); bao["dat_lai"] += 1
                    da.append(b)
                with torch.no_grad():
                    DEM["danh_gia"] += len(lon)
                    if sum(float(phat_lon(b, lon, day_du=True)) for b in lon) <= NG: break
                lay_mau_lon(HT[DT["ngau_nhien"]["goc_lam_lai"]], g)
        thu_loai_nho(nh, lon, g, bao)
    _t = time.time(); pl, pn = phat_day_du(lon, nh); L = sum(pl.values()) + sum(pn.values()); t_day_du = max(t_day_du, time.time() - _t)
    DEM["danh_gia"] += len(lon) + nh.n
    if os.environ.get("GO_LOI"): print("k", k, "phat", L, round(time.time() - t00), flush=True)
    if os.environ.get("GO_LOI") and os.environ.get("E3_PHAT_GOC") and ARG.eval_dung_tu_dau:
        E3["bat"] = False; _pg, _ng = phat_day_du(lon, nh); E3["bat"] = True; print("k", k, "phat_goc", sum(_pg.values()) + sum(_ng.values()), flush=True)
    if tot is None or L < tot[0]: tot = (L, lon, nh, k, pl, pn)
E3_BC = {"tu_the_tu": ARG.tu_the_tu, "eval_dung": ARG.eval_dung, "eval_dung_tu_dau": ARG.eval_dung_tu_dau, "G2": E3["G2"], "GB2": E3["GB2"]}
if ARG.tu_the_tu:
    import hashlib as _hl2
    _TT = json.load(open(ARG.tu_the_tu)); k_n = int(_TT["bao_cao_bo_tri"]["cau_hinh_chon"])
    assert _TT["seed"] == ARG.seed and _TT["cach"] == ARG.cach and _TT["thu_vien"] == ARG.thu_vien and _TT.get("khoi_tao", "luat") == ARG.khoi_tao, "lệnh khác lượt gốc"
    assert _TT["so_luong"] == {l: SO[l] for l in LON + NHO}, "số lượng khác lượt gốc"
    lon, nh, g = khoi_tao(k_n); HT = {b.loai: b for b in lon if b.loai in DON}
    assert set(_TT["tu_the"]) == set(b.ten for b in lon) | set(nh.ten), "tập tên khác"
    with torch.no_grad():
        for b in lon:
            r_ = _TT["tu_the"][b.ten]; assert r_["khuon"] == b.khuon, b.ten; assert abs(float(b.sc()) - r_["s"]) < 2e-6, b.ten
            b.t = T(r_["t"]); b.r = _r6(T(r_["R"]))
        S0_ = nh.sc()
        for i in range(nh.n):
            r_ = _TT["tu_the"][nh.ten[i]]; assert r_["khuon"] == TEN_K[int(nh.tid[i])], nh.ten[i]; assert abs(float(S0_[i]) - r_["s"]) < 2e-6, nh.ten[i]
            nh.t[i] = T(r_["t"]); nh.r[i] = _r6(T(r_["R"]))
        lech_R = max([float((b.R() - T(_TT["tu_the"][b.ten]["R"])).abs().max()) for b in lon] + [float((nh.R()[i] - T(_TT["tu_the"][nh.ten[i]]["R"])).abs().max()) for i in range(nh.n)])
    pl0, pn0 = phat_day_du(lon, nh); _xay_xa(); E3["bat"] = True; pl1, pn1 = phat_day_du(lon, nh); E3["bat"] = False
    E3_BC.update({"md5_tu_the_tu": _hl2.md5(open(ARG.tu_the_tu, "rb").read()).hexdigest(), "cau_hinh_chon_luu": k_n, "lech_R_max": lech_R,
                  "P_luu": _TT["bao_cao_bo_tri"]["phat_cung_day_du"], "P_goc": sum(pl0.values()) + sum(pn0.values()), "P_dung": sum(pl1.values()) + sum(pn1.values()),
                  "phat_goc": {**pl0, **pn0}, "phat_dung": {**pl1, **pn1}, "kiem_xa": E3["kiem_xa"], "giay_xay_xa": E3["giay_xay"]})
    if ARG.chi_P:
        json.dump({"ban": "botri_chung_v122e3", "chi_P": True, "seed": ARG.seed, "ban_e3": E3_BC}, open(f"{OUT}/thong_tin_P.json", "w"), ensure_ascii=False, indent=1)
        print("XONG_BOTRI chi_P", E3_BC["P_goc"], E3_BC["P_dung"], flush=True); sys.exit(0)
    E3["bat"] = bool(ARG.eval_dung); pl, pn = (pl1, pn1) if ARG.eval_dung else (pl0, pn0)
    tot = (sum(pl.values()) + sum(pn.values()), lon, nh, k_n, pl, pn)
if ARG.eval_dung_tu_dau: E3_BC.update({"kiem_xa": E3["kiem_xa"], "giay_xay_xa": E3["giay_xay"]})
assert not (ARG.eval_dung and not ARG.tu_the_tu) and not (ARG.chi_P and not ARG.tu_the_tu), "cờ E3 sai"
torch.cuda.synchronize(); GIAY_DAT = time.time() - t_dat0
_, lon, nh, k_chon, pl, pn = tot; HT = {b.loai: b for b in lon if b.loai in DON}
bao["cau_hinh_chon"] = k_chon; bao["phat_cung_day_du"] = tot[0]; t_botri = time.time() - t00
if True:
    _rc, _rg = torch.get_rng_state(), torch.cuda.get_rng_state()
    _km = ARG.khong_mem; ARG.khong_mem = False
    with torch.no_grad(): _ml = float(phat_mem_lon()); _mn = float(phat_mem_nho(nh, torch.zeros((), dtype=F64, device=DEV)))
    ARG.khong_mem = _km; bao["khong_mem"] = _km; bao["khoi_tao"] = ARG.khoi_tao; bao["un_tham_so"] = {"T1": ARG.un_T1, "p0": ARG.un_p0, "buoc_t": ARG.un_buoc_t, "buoc_r": ARG.un_buoc_r}
    torch.set_rng_state(_rc); torch.cuda.set_rng_state(_rg)
    bao["J_truoc_chieu"] = {"J": W * tot[0] + _ml + _mn, "Lambda_x_phat_cung": W * tot[0], "mem_lon": _ml, "mem_nho": _mn}
    bao["giay_dat"] = round(GIAY_DAT, 3); bao["danh_gia_phat"] = DEM["danh_gia"]; bao["so_lo"] = SH["thu_loai_so_lo"]
    bao["ngan_sach_nhanh"] = {"so_lo": SH["thu_loai_so_lo"], "un_vong_lon": ARG.un_vong_lon, "un_vong_nho": ARG.un_vong_nho, "ngan_sach_giay": ARG.ngan_sach_giay, "K": ARG.K, "buoc": ARG.buoc}
if ARG.chi_dat:
    json.dump({"ban": "botri_chung_v122e3", "ban_e3": E3_BC, "ban_v121": V121, "chi_dat": True, "dac_ta": ARG.dac_ta, "seed": ARG.seed, "cach": ARG.cach, "thu_vien": ARG.thu_vien, "so_luong": {l: SO[l] for l in LON + NHO},
               "bao_cao_bo_tri": bao, "giay_dat": bao["giay_dat"], "gpu": os.environ.get("CUDA_VISIBLE_DEVICES")}, open(f"{OUT}/thong_tin_dat.json", "w"), ensure_ascii=False, indent=1)
    print("XONG_BOTRI chi_dat", ARG.cach, "giay_dat", bao["giay_dat"], "J", bao["J_truoc_chieu"]["J"], flush=True); sys.exit(0)
phat_truoc = {**pl, **pn}
LAPK_KKT = 10
TKC = {"luot": [], "khong_hoi_tu": 0, "kkt_max": 0.0, "lap_max": 0}
def dykstra_kkt(x, sl, cf, Bs, hop, wsr, nrm, step_c, lap_max, tol_kkt, ten_luot="?", cap=None):
    """Dykstra trên kk hàng mỗi nút. Dừng TỪNG nút (nút xong bị loại khỏi tập hoạt động) khi: khả thi nguyên thuỷ ≤ tol_chieu VÀ
    cận khoảng tới phép chiếu √(2·Σ_r λ_r·khe_r) ≤ tol_kkt (metric W; λ_r = bước cuối của hàng r, bất biến x = x0 + Σ λ_r W⁻¹a_r ⇒
    dừng đúng; ½‖x − x*‖²_W ≤ khe đối ngẫu Σλ·khe khi x khả thi). cap: số vòng tối đa riêng từng nút (nút vô nghiệm theo Mệnh đề 4).
    Trả x, số vòng mỗi nút, cận √(2·khe) mỗi nút (∞ nếu chưa khả thi)."""
    Nv, kk = Bs.shape; X0 = x; out = x.clone(); lapn = torch.full((Nv,), lap_max, dtype=torch.long, device=DEV); rk = torch.full((Nv,), float("inf"), dtype=F64, device=DEV)
    ok_all = torch.zeros(Nv, dtype=torch.bool, device=DEV); cap = torch.full((Nv,), lap_max, dtype=torch.long, device=DEV) if cap is None else cap
    act = torch.arange(Nv, device=DEV); hopf = hop.to(F64)
    Pd = torch.zeros(kk, Nv, x.shape[1], dtype=F64, device=DEV); mu = torch.zeros(kk, Nv, dtype=F64, device=DEV)
    lap = -1
    for lap in range(lap_max):
        if not act.numel(): break
        for r in range(kk):
            y = x + Pd[r]; s_ = sl[:, r]; ay = (cf[:, r] * y.gather(1, s_)).sum(1)
            st_ = torch.clamp(Bs[:, r] - ay, min=0) / nrm[:, r] * hopf[:, r]
            xn = y.scatter_add(1, s_, st_[:, None] * step_c[:, r]); Pd[r] = y - xn; x = xn; mu[r] = st_
        if lap % LAPK_KKT == LAPK_KKT - 1 or lap == lap_max - 1:
            ay_ = (cf * x.gather(1, sl.reshape(x.shape[0], -1)).reshape(x.shape[0], kk, 2)).sum(2); kh = (ay_ - Bs) * hopf
            prim = torch.clamp(-kh, min=0).max(1).values; gap = (mu.T * torch.clamp(kh, min=0)).sum(1); bd = torch.sqrt(2 * gap)
            ok = (prim < TOLC) & (bd <= tol_kkt); xong_ = ok | (lap + 1 >= cap[act])
            rk[act] = torch.where(prim < TOLC, bd, torch.full_like(bd, float("inf"))); ok_all[act[ok]] = True; lapn[act[ok]] = lap + 1
            if bool(xong_.any()):
                out[act[xong_]] = x[xong_]; giu = ~xong_
                act = act[giu]; x = x[giu]; Pd = Pd[:, giu]; mu = mu[:, giu]; sl = sl[giu]; cf = cf[giu]; Bs = Bs[giu]; hopf = hopf[giu]; nrm = nrm[giu]; step_c = step_c[giu]
    if act.numel(): out[act] = x
    nk = int((~ok_all).sum())
    TKC["luot"].append((ten_luot, Nv, kk, lap + 1, lapn.cpu(), float(rk[ok_all].max()) if bool(ok_all.any()) else 0.0)); TKC["khong_hoi_tu"] += nk
    TKC["kkt_max"] = max(TKC["kkt_max"], float(rk[ok_all].max()) if bool(ok_all.any()) else 0.0); TKC["lap_max"] = max(TKC["lap_max"], lap + 1)
    return out, lapn, rk
def chuan_bi_hang(Uv, Bv, sub, RMAX, hang=None):
    """chọn ≤ RMAX hàng gần nhất (sub < dia_phuong_h·h1) mỗi nút; trả các tensor cho dykstra_kkt. hang: chỉ số hàng toàn cục ứng với cột của sub."""
    Nv = sub.shape[0]; kk = int(min(RMAX, int((sub < DPH * h1).sum(1).max())))
    vals, rs_ = torch.topk(-sub, kk, dim=1); hop = (-vals) < DPH * h1; rs = rs_ if hang is None else hang[rs_]
    vid = RV[rs]; cf = RC[rs] * hop[..., None]
    flat = vid.reshape(Nv, -1); sv, perm = torch.sort(flat, 1); moi = torch.ones_like(sv, dtype=torch.bool); moi[:, 1:] = sv[:, 1:] != sv[:, :-1]
    slot_s = torch.cumsum(moi.long(), 1) - 1; slot = torch.empty_like(slot_s).scatter_(1, perm, slot_s); Ks = int(slot_s.max()) + 1
    x = torch.zeros(Nv, Ks, dtype=F64, device=DEV).scatter_(1, slot, Uv.gather(1, flat))
    ws = torch.ones(Nv, Ks, dtype=F64, device=DEV).scatter_(1, slot, w[flat]); var_slot = torch.full((Nv, Ks), DUM, device=DEV).scatter_(1, slot, flat)
    sl = slot.reshape(Nv, kk, 2); Bs = Bv.gather(1, rs_); wsr = ws.gather(1, sl.reshape(Nv, -1)).reshape(Nv, kk, 2)
    nrm = (cf ** 2 / wsr).sum(2); nrm = torch.where(nrm > 0, nrm, torch.ones_like(nrm)); step_c = cf / wsr
    return x, sl, cf, Bs, hop, wsr, nrm, step_c, var_slot
def kha_thi_dong(B, RVh, RCh, RKh, nvar, tol=None):
    """Mệnh đề 4 (kiểm khả thi đóng) tại các nút (hàng của B): biến chỉ có hệ số +1 đặt bằng trần giữ chỗ; biến cha (hệ số −1 ở hàng chứa)
    đặt bằng cận dưới lớn nhất (không có cận dưới: −∞); rồi kiểm mọi hàng. Trả mặt nạ nút VÔ NGHIỆM."""
    tol = TOLC if tol is None else tol; N = B.shape[0]; BIG = 1e6
    ub = torch.full((N, nvar), BIG, dtype=F64, device=DEV); m2 = (RKh == 2).nonzero().flatten()
    if m2.numel():
        Bm = B[:, m2]; ub.scatter_reduce_(1, RVh[m2, 0][None].expand(N, -1), torch.where(Bm > -1e8, -Bm, torch.full_like(Bm, BIG)), "amin")
    x = ub.clone(); x[:, nvar - 1] = 0.0; bad = torch.zeros(N, dtype=torch.bool, device=DEV)
    m0 = RKh == 0; cha = sorted(set(RVh[m0 & (RCh[:, 1] < 0), 1].tolist()))
    for p_ in cha:
        lb = torch.full((N,), -BIG, dtype=F64, device=DEV); upb = ub[:, p_].clone()
        r1 = ((RKh == 1) & (RVh[:, 0] == p_)).nonzero().flatten()
        if r1.numel(): lb = torch.maximum(lb, B[:, r1].max(1).values)
        for c_ in (0, 1):
            rr = (m0 & (RVh[:, c_] == p_) & (RCh[:, 1] > 0)).nonzero().flatten()
            if rr.numel(): lb = torch.maximum(lb, (B[:, rr] - x[:, RVh[rr, 1 - c_]]).max(1).values)
        rr = (m0 & (RVh[:, 1] == p_) & (RCh[:, 1] < 0)).nonzero().flatten()
        if rr.numel(): upb = torch.minimum(upb, (x[:, RVh[rr, 0]] - B[:, rr]).min(1).values)
        bad |= lb > upb + tol; x[:, p_] = torch.where(lb > -BIG / 2, lb, torch.full_like(lb, -BIG))
    lhs = RCh[:, 0][None] * x[:, RVh[:, 0]] + RCh[:, 1][None] * x[:, RVh[:, 1]]
    act = (B > -1e8) | (RKh[None] != 2)
    return bad | ((lhs < B - tol) & act).any(1)
def tom_tat_chieu():
    out = []
    for ten_, Nv, kk, lap, lapn, rkm in TKC["luot"]:
        if Nv: out.append({"luot": ten_, "so_nut": Nv, "so_hang": kk, "vong_khoi": lap, "vong_tb": round(float(lapn.double().mean()), 1), "vong_p99": int(lapn.double().quantile(0.99)), "kkt_max": rkm})
    gom_ = {}
    for x_ in out:
        g_ = gom_.setdefault(x_["luot"], {"so_nut": 0, "so_khoi": 0, "vong_khoi_max": 0, "vong_tb_trong_so": 0.0, "kkt_max": 0.0})
        g_["so_nut"] += x_["so_nut"]; g_["so_khoi"] += 1; g_["vong_khoi_max"] = max(g_["vong_khoi_max"], x_["vong_khoi"]); g_["vong_tb_trong_so"] += x_["vong_tb"] * x_["so_nut"]; g_["kkt_max"] = max(g_["kkt_max"], x_["kkt_max"])
    for g_ in gom_.values(): g_["vong_tb_trong_so"] = round(g_["vong_tb_trong_so"] / max(g_["so_nut"], 1), 1)
    return {"theo_luot": gom_, "khong_hoi_tu_nut": TKC["khong_hoi_tu"], "kkt_max": TKC["kkt_max"], "vong_toi_da_mot_khoi": TKC["lap_max"], "lap_kkt": ARG.lap_kkt, "kkt_tol": ARG.kkt_tol, "tol_nguyen_thuy": SH["tol_chieu"], "tieu_chi": "khả thi ≤ tol_chieu và √(2·Σλ·khe) ≤ kkt_tol (cận ‖u−u*‖_W)"}
import trimesh.creation as _tc
def dung_tang2(cha_ds, VBm):
    """Đối tượng con tầng 2 (mục 'tang2' của đặc tả). cha_ds: [(tên, loại, khuôn, t, R, s, tâm sâu)] (numpy). Trả {tên: mesh}, {tên: loại}, {tên: cha}, báo cáo."""
    ME, LO_, CH_, BC = {}, {}, {}, []
    for e in DT.get("tang2", []):
        for (ten, lo_, kh, t_, R_, s_, c_) in [d for d in cha_ds if d[1] == e["cha"]]:
            g = np.random.default_rng([ARG.seed, sum(map(ord, e["loai"])), sum(map(ord, ten))]); parts = []
            if e["kieu"] == "hai_tru_vuong_goc":
                n_ = int(e["so"][0]); assert n_ == 2 == int(e["so"][1])
                L_ = g.uniform(*e["dai"], size=2); r_ = g.uniform(*e["rong"], size=2) / 2; d_ = r_[0] + r_[1]
                for i_, (ax_, z_) in enumerate(((0, -d_ / 2), (1, d_ / 2))):
                    m = _tc.capsule(height=max(L_[i_] - 2 * r_[i_], 1e-3), radius=r_[i_], count=[32, 32])
                    rot = np.eye(4); rot[:3, :3] = (np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]]) if ax_ == 0 else np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])).T
                    m.apply_transform(rot); m.apply_translation([0, 0, z_]); parts.append((m, {"dai": float(L_[i_]), "rong": float(2 * r_[i_]), "truc_cuc_bo": int(ax_)}))
                info = {"goc_do": float(e["goc_do"])}
            elif e["kieu"] == "chong_dia":
                kd_ = KD.get(kh, {})
                if kd_.get("hinh") != "elip": BC.append({"cha": ten, "bo_qua": f"khuôn cha {kh} không phải elip"}); continue
                a_ = np.array(kd_["truc"]) * s_ - h1; ix_ = int(np.argmin(a_)); o_ = [j for j in range(3) if j != ix_]; amin = float(a_[o_].min())
                lo, hi = e["so"]; n_ = int(g.integers(lo, hi + 1)); W_ = g.uniform(*e["rong"], size=n_)
                if W_.max() / 2 >= amin: BC.append({"cha": ten, "bo_qua": "túi rộng hơn khối Golgi"}); continue
                H_ = 2 * a_[ix_] * math.sqrt(1 - (W_.max() / 2 / amin) ** 2); tau = H_ / (2 * n_ - 1)
                for k_ in range(n_):
                    m = _tc.cylinder(radius=W_[k_] / 2, height=tau, sections=48); xk = -H_ / 2 + tau / 2 + 2 * k_ * tau
                    rot = np.eye(4)
                    if ix_ == 0: rot[:3, :3] = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]]).T
                    elif ix_ == 1: rot[:3, :3] = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]]).T
                    m.apply_transform(rot); off = np.zeros(3); off[ix_] = xk; m.apply_translation(off); parts.append((m, {"rong": float(W_[k_]), "day": float(tau)}))
                info = {"so_tui": n_, "day_tui": float(tau), "H": float(H_)}
            else: raise NotImplementedError(e["kieu"])
            ten_con = []
            for i_, (m, inf) in enumerate(parts):
                V_ = (np.asarray(m.vertices) @ R_.T) + t_; m2 = trimesh.Trimesh(V_, np.asarray(m.faces), process=True)
                nm = f"{e['loai']}_{len([k for k in LO_ if LO_[k] == e['loai']]) + 1}"; ME[nm] = m2; LO_[nm] = e["loai"]; CH_[nm] = ten; ten_con.append(nm)
                inf["ngoai_cha_dinh"] = int((igl.signed_distance(V_, np.asarray(VBm[ten].vertices, float), np.asarray(VBm[ten].faces, np.int64))[0] > 0).sum()) if VBm.get(ten) is not None else None
                info.setdefault("con", []).append(dict(inf, ten=nm))
            if e["kieu"] == "hai_tru_vuong_goc":
                info["goc_do_theo_cau_tao"] = float(e["goc_do"]); assert abs(e["goc_do"] - 90) < 1e-9
            BC.append(dict(info, cha=ten, loai=e["loai"], nguon=e["nguon"]))
    return ME, LO_, CH_, BC
if CUCBO:
    _f = os.path.join(os.path.dirname(os.path.abspath(__file__)), "botri_cucbo_v121.py"); exec(compile(open(_f).read(), _f, "exec")); sys.exit(0)
with torch.no_grad():
    mem = {}
    for r in QH.get("lech_tam", []): mem[r["bao_cao"]] = float(HT[r["a"]].t.norm()) / r["ban_kinh"]
    for r in QH.get("gan", []): a = HT[r["a"]]; mem[r["bao_cao"]] = float(HT[r["b"]].phi(a.t[None])[0]) - float(a.sc())
    mem["do_mat_mat"] = True
    for r in QH.get("gan_nho", []):
        mk = nh.mask(r["a"]); mem[r["bao_cao"]] = float(HT[r["b"]].phi(nh.t[mk]).mean()) if bool(mk.any()) else 0.0

ax1 = torch.linspace(-A + h1 / 2, A - h1 / 2, N1, device=DEV, dtype=F64); X = torch.stack(torch.meshgrid(ax1, ax1, ax1, indexing="ij"), -1).reshape(-1, 3); NX = X.shape[0]
TK = SH["truong_khoi"]
def truong(f):
    out = torch.empty(NX, dtype=torch.float32, device=DEV)
    for a0 in range(0, NX, TK): out[a0:a0 + TK] = f(X[a0:a0 + TK]).float()
    return out
psi, tam, co, rbk, loai = {}, {}, {}, {}, {}; cam_g = {}
with torch.no_grad():
    psi_vo = PSI_VO_NUT.clone() if VO["hinh"] == "cau_nhieu" else truong(phi_vo)
    for b in lon: psi[b.ten] = truong(b.phi); tam[b.ten] = b.tam().cpu().numpy(); co[b.ten] = b.r_trong(); rbk[b.ten] = b.rb(); loai[b.ten] = b.loai
    S_ = nh.sc(); R_ = nh.R()
    for i in range(nh.n):
        ti, Ri, Si, di = nh.t[i], R_[i], S_[i], int(nh.tid[i])
        psi[nh.ten[i]] = truong(lambda P: Si * eval_at(((P - ti) @ Ri) / Si, torch.full(P.shape[:-1], di, device=DEV)))
        tam[nh.ten[i]] = ((SAU[di] * Si) @ Ri.T + ti).cpu().numpy(); co[nh.ten[i]] = float(RTR[di] * Si); rbk[nh.ten[i]] = float(RB[di] * Si); loai[nh.ten[i]] = nh.loai[i]
    for vt in PHU:
        v = VANH[vt]; q = HT[v["quanh"]]
        psi[vt] = truong(lambda P, v=v: phi_vanh(P, v))
        tam[vt] = (q.t - u_kh(v) * (q.sc() + v["tam_lui"])).cpu().numpy() if v.get("khuyet") else q.t.cpu().numpy()
        co[vt] = v["co"]; rbk[vt] = v["rb"]; loai[vt] = vt
        cam_g[vt] = truong(lambda P, v=v: phi_cam(P, v))
bien = list(psi); ix = {k: i for i, k in enumerate(bien)}; nv = len(bien); DUM = nv
ctr = {b.ten: b.t.cpu().numpy() for b in lon}; ctr.update({nh.ten[i]: nh.t[i].cpu().numpy() for i in range(nh.n)})
rv, rc, rk, rval, rtam, mota = [], [], [], [], [], []
rrb = []
def them_hang(a, ca, b, cb, kind, val, tm=(0, 0, 0), rb_=0.0): rv.append([a, b]); rc.append([ca, cb]); rk.append(kind); rval.append(val); rtam.append(tm); rrb.append(rb_)
def ten_mot(l): return ten_ban(l, 0)
for k in bien:
    r = CHUA[loai[k]]
    if r["b"] == VO_TEN: them_hang(ix[k], 1, DUM, 0, 1, r["t"] + DEL_M); mota.append(("chứa", k, VO_TEN, r["t"], r["nguon"]))
    else: p_ = ten_mot(r["b"]); them_hang(ix[k], 1, ix[p_], -1, 0, r["t"] + bu); mota.append(("chứa", k, p_, r["t"], r["nguon"]))
RCH = QH["roi_chieu"]; S_C = RCH["s"]; MCAP = {frozenset(p) for p in RCH["mien_cap"]}; MCUNG = set(RCH["mien_cung_loai"])
KCL_C = RCH.get("khe_cung_loai", {}); assert not (set(KCL_C) & MCUNG); NG_CL = RCH.get("nguon_cung_loai", RCH["nguon"]); n_cl = 0
ngoai = [k for k in bien if loai[k] in RCH["loai"]]
LH = SH["loc_cap_h"]
TAU_C0 = SH["tol_con"]; N_C0 = []
with torch.no_grad():
    for i_, a in enumerate(ngoai):
        for b in ngoai[i_ + 1:]:
            la, lb = loai[a], loai[b]
            if la != lb and frozenset((la, lb)) in MCAP:
                if float((psi[a] + psi[b]).min()) > TAU_C0 + LH * h1: continue
                them_hang(ix[a], 1, ix[b], 1, 0, TAU_C0); mota.append(("rời", a, b, 0.0, "M.c₀ (v12.2: miễn khe, không miễn giao)")); N_C0.append((a, b)); continue
            if la == lb and la in MCUNG: continue
            cl_ = la == lb and la in KCL_C; s = KCL_C[la] if cl_ else S_C
            if la not in PHU and lb not in PHU and np.linalg.norm(ctr[a] - ctr[b]) > rbk[a] + rbk[b] + s + bu + LH * h1: continue
            if float((psi[a] + psi[b]).min()) > s + bu + LH * h1: continue
            them_hang(ix[a], 1, ix[b], 1, 0, s + bu); mota.append(("rời", a, b, s, NG_CL if cl_ else RCH["nguon"])); n_cl += cl_
    for lt, s in KCL_C.items():
        if lt in RCH["loai"]: continue
        cl = [k for k in bien if loai[k] == lt]
        for i_, a in enumerate(cl):
            for b in cl[i_ + 1:]:
                if np.linalg.norm(ctr[a] - ctr[b]) > rbk[a] + rbk[b] + s + bu + LH * h1: continue
                if float((psi[a] + psi[b]).min()) > s + bu + LH * h1: continue
                them_hang(ix[a], 1, ix[b], 1, 0, s + bu); mota.append(("rời", a, b, s, NG_CL)); n_cl += 1
gc_ = SH["giu_cho"]
for k in bien:
    rr = min(max(gc_[0] * bu / gc_[1], gc_[2] * co[k]), gc_[3] * co[k]); them_hang(ix[k], -1, DUM, 0, 2, rr, tam[k], rbk[k] + gc_[4] * h1); mota.append(("giữ chỗ", k, None, 0, ""))
print("SO_HANG", len(rv), "roi", sum(1 for x in mota if x[0] == "rời"), "cung_loai", n_cl, "v122_khong_giao", len(N_C0), flush=True)
RRB = T(rrb); RV = torch.tensor(rv, device=DEV); RC = T(rc); RK = torch.tensor(rk, device=DEV); RVAL = T(rval); RTAM = T(np.array(rtam)); nr = RV.shape[0]
the_tich = {k: float((psi[k] < 0).sum()) * h1 ** 3 for k in bien}
w = torch.ones(nv + 1, dtype=F64, device=DEV)
if ARG.trong_so:
    for k in bien: w[ix[k]] = max(the_tich[k], SH["the_tich_san"]) ** (-1 / SH["trong_so_mu_mau"])
    w[:nv] = w[:nv] / w[:nv].min()
PHI = {k: psi[k].clone() for k in bien}
def res_khoi(j, nguon):
    Uj = torch.stack([nguon[k][j].double() for k in bien] + [torch.zeros(len(j), dtype=F64, device=DEV)], 1)
    B = RVAL[None].expand(len(j), nr).clone(); m1 = RK == 1; B[:, m1] = B[:, m1] + psi_vo[j].double()[:, None]
    m2 = RK == 2; dd_ = torch.linalg.norm(X[j][:, None, :] - RTAM[m2][None], dim=2)
    B[:, m2] = torch.where(dd_ < RRB[m2][None], -(dd_ - RVAL[m2][None]), torch.full_like(dd_, -1e9))
    res = (Uj[:, RV[:, 0]] * RC[:, 0] + Uj[:, RV[:, 1]] * RC[:, 1]) - B
    return Uj, B, res
_kh = SH["KH"]; n_vp = 0; t0 = time.time(); KH = max(_kh[0], min(_kh[1], int(_kh[2] / (_kh[3] * max(nr, 1)))))
DPH = SH["dia_phuong_h"]; TOLC = SH["tol_chieu"]; LAPK = SH["lap_kiem"]
def chieu_diem(jj, Uv, Bv, sub, RMAX):
    """Dykstra cục bộ cho các điểm jj (giá trị xuất phát Uv): giữ ≤ RMAX hàng gần nhất (res < 4h1)."""
    Nv = sub.shape[0]; kk = int(min(RMAX, int((sub < DPH * h1).sum(1).max())))
    vals, rs = torch.topk(-sub, kk, dim=1); hop = (-vals) < DPH * h1
    vid = RV[rs]; cf = RC[rs] * hop[..., None]
    flat = vid.reshape(Nv, -1); sv, perm = torch.sort(flat, 1); moi = torch.ones_like(sv, dtype=torch.bool); moi[:, 1:] = sv[:, 1:] != sv[:, :-1]
    slot_s = torch.cumsum(moi.long(), 1) - 1; slot = torch.empty_like(slot_s).scatter_(1, perm, slot_s); Ks = int(slot_s.max()) + 1
    x = torch.zeros(Nv, Ks, dtype=F64, device=DEV).scatter_(1, slot, Uv.gather(1, flat))
    ws = torch.ones(Nv, Ks, dtype=F64, device=DEV).scatter_(1, slot, w[flat]); var_slot = torch.full((Nv, Ks), DUM, device=DEV).scatter_(1, slot, flat)
    sl = slot.reshape(Nv, kk, 2); Bs = Bv.gather(1, rs); wsr = ws.gather(1, sl.reshape(Nv, -1)).reshape(Nv, kk, 2)
    nrm = (cf ** 2 / wsr).sum(2); nrm = torch.where(nrm > 0, nrm, torch.ones_like(nrm)); step_c = cf / wsr
    Pd = torch.zeros(kk, Nv, Ks, dtype=F64, device=DEV)
    for lap in range(SH["lap_toi_da"]):
        if lap % LAPK == LAPK - 1:
            ay_ = (cf * x.gather(1, sl.reshape(Nv, -1)).reshape(Nv, kk, 2)).sum(2)
            if float((torch.clamp(Bs - ay_, min=0) * hop).max()) < TOLC: break
        for r in range(kk):
            y = x + Pd[r]; s_ = sl[:, r]; ay = (cf[:, r] * y.gather(1, s_)).sum(1)
            st_ = torch.clamp(Bs[:, r] - ay, min=0) / nrm[:, r] * hop[:, r]
            xn = y.scatter_add(1, s_, st_[:, None] * step_c[:, r]); Pd[r] = y - xn; x = xn
    for vk in var_slot.unique().tolist():
        if vk == DUM: continue
        m = var_slot == vk; pi_, si_ = m.nonzero(as_tuple=True); PHI[bien[vk]][jj[pi_]] = x[pi_, si_].float()
RM_ = SH["RMAX"]; so_luot = [0] * len(RM_); khong_kha_thi = 0; CHI_VP = []
if ARG.chieu == "cu":
    with torch.no_grad():
        for a0 in range(0, NX, KH):
            j = torch.arange(a0, min(NX, a0 + KH), device=DEV); Uj, B, res = res_khoi(j, psi); vi = (res < -TOLC).any(1)
            if not vi.any(): continue
            n_vp += int(vi.sum())
            if ARG.khong_chieu: continue
            chieu_diem(j[vi], Uj[vi], B[vi], res[vi], RM_[0]); so_luot[0] += int(vi.sum())
            for luot in range(1, len(RM_)):
                Uj2, B2, res2 = res_khoi(j, PHI); v2 = (res2 < -TOLC).any(1)
                if not v2.any(): break
                chieu_diem(j[v2], Uj2[v2], B2[v2], res2[v2], RM_[luot]); so_luot[luot] += int(v2.sum())
        torch.cuda.synchronize(); t_chieu = time.time() - t0
        con_vp = 0
        for a0 in range(0, NX, KH):
            j = torch.arange(a0, min(NX, a0 + KH), device=DEV); _, _, res = res_khoi(j, PHI); bad = (res < -SH["tol_con"]); con_vp += int(bad.any(1).sum())
            if os.environ.get("GO_LOI") and bad.any():
                for kd_ in (0, 1, 2):
                    mk = bad[:, RK == kd_]
                    if mk.any(): print("con_vp loai", kd_, int(mk.any(1).sum()), "min", float(res[:, RK == kd_].min()), flush=True)


else:
  with torch.no_grad():
    for a0 in range(0, NX, KH):
        j = torch.arange(a0, min(NX, a0 + KH), device=DEV); Uj, B, res = res_khoi(j, psi); vi = (res < -TOLC).any(1)
        if not vi.any(): continue
        n_vp += int(vi.sum()); vn = torch.zeros_like(vi); vn[vi] = kha_thi_dong(B[vi], RV, RC, RK, nv + 1); khong_kha_thi += int(vn.sum())
        capn = torch.where(vn, torch.full_like(vn, SH["lap_toi_da"], dtype=torch.long), torch.full_like(vn, ARG.lap_kkt, dtype=torch.long))
        ep = torch.zeros_like(res, dtype=torch.bool)
        if ARG.kiem_chieu: CHI_VP.append(j[vi].cpu())
        if ARG.khong_chieu: continue
        cham = set()
        def _chieu(idx, sub_, RMAX_, ten_):
            x, sl, cf, Bs, hop, wsr, nrm, step_c, var_slot = chuan_bi_hang(Uj[idx], B[idx], sub_, RMAX_)
            x, _, _ = dykstra_kkt(x, sl, cf, Bs, hop, wsr, nrm, step_c, ARG.lap_kkt, ARG.kkt_tol, ten_, cap=capn[idx])
            jj = j[idx]
            for vk in cham:
                PHI[bien[vk]][jj] = psi[bien[vk]][jj]
            for vk in var_slot.unique().tolist():
                if vk == DUM: continue
                cham.add(vk); m = var_slot == vk; pi_, si_ = m.nonzero(as_tuple=True); PHI[bien[vk]][jj[pi_]] = x[pi_, si_].float()
        iv = vi.nonzero().flatten(); _chieu(iv, res[iv], RM_[0], "1"); so_luot[0] += int(iv.numel())
        for luot in range(1, len(RM_)):
            _, _, res2 = res_khoi(j, PHI); v2 = (res2 < -SH["tol_con"]).any(1)
            if not v2.any(): break
            ep |= res2 < -SH["tol_con"]; iv2 = v2.nonzero().flatten(); sub2 = torch.where(ep[iv2], torch.full_like(res2[iv2], -1e9), res[iv2])
            _chieu(iv2, sub2, RM_[luot], str(luot + 1)); so_luot[luot] += int(iv2.numel())
    torch.cuda.synchronize(); t_chieu = time.time() - t0
    con_vp = 0; con_theo_loai = {0: 0, 1: 0, 2: 0}; con_max = 0.0
    for a0 in range(0, NX, KH):
        j = torch.arange(a0, min(NX, a0 + KH), device=DEV); _, _, res = res_khoi(j, PHI); bad = (res < -SH["tol_con"]); con_vp += int(bad.any(1).sum())
        if bad.any():
            con_max = max(con_max, float(-res.min()))
            for kd_ in (0, 1, 2): con_theo_loai[kd_] += int(bad[:, RK == kd_].any(1).sum())
KIEM_CHIEU = None
if ARG.kiem_chieu and CHI_VP and not ARG.khong_chieu:
    with torch.no_grad():
        allv = torch.cat(CHI_VP); gk_ = torch.Generator().manual_seed(ARG.seed); pick = allv[torch.randperm(allv.numel(), generator=gk_)[:ARG.kiem_chieu]].to(DEV)
        _, Bk_, _ = res_khoi(pick, psi); pick = pick[~kha_thi_dong(Bk_, RV, RC, RK, nv + 1)]
        dmax = 0.0; dW = 0.0; kk_max = 0; vp_ref = 0; luot_cu = len(TKC["luot"]); nk0 = TKC["khong_hoi_tu"]
        for c0 in range(0, pick.numel(), 2000):
            jj = pick[c0:c0 + 2000]; Uj, B, res = res_khoi(jj, psi); _, _, resP = res_khoi(jj, PHI)
            sub = torch.where(resP < -TOLC, torch.full_like(resP, -1e9), res / 5)
            x, sl, cf, Bs, hop, wsr, nrm, step_c, var_slot = chuan_bi_hang(Uj, B, sub, 10 ** 6)
            kk_max = max(kk_max, Bs.shape[1]); x, _, _ = dykstra_kkt(x, sl, cf, Bs, hop, wsr, nrm, step_c, 20000, 1e-9, "tham_chieu")
            cur = torch.zeros_like(x); cur0 = torch.zeros_like(x)
            for vk in var_slot.unique().tolist():
                if vk == DUM: continue
                m = var_slot == vk; pi_, si_ = m.nonzero(as_tuple=True); cur[pi_, si_] = PHI[bien[vk]][jj[pi_]].double(); cur0[pi_, si_] = psi[bien[vk]][jj[pi_]].double()
            real = var_slot != DUM; dd = (x - cur).abs() * real; dmax = max(dmax, float(dd.max()))
            Ur = Uj.clone(); Ur.scatter_(1, var_slot, x); Ur[:, DUM] = 0.0; rr_ = (Ur[:, RV[:, 0]] * RC[:, 0] + Ur[:, RV[:, 1]] * RC[:, 1]) - B
            vp_ref = vp_ref + int((rr_ < -SH["tol_con"]).any(1).sum())
            dW = max(dW, float(((x - cur) ** 2 * w[var_slot] * real).sum(1).sqrt().max()))
        KIEM_CHIEU = {"so_nut": int(pick.numel()), "dmax_um": dmax, "dmax_W": dW, "so_hang_toi_da": kk_max, "khong_hoi_tu_tham_chieu": TKC["khong_hoi_tu"] - nk0, "nut_tham_chieu_vi_pham_hang_nao_do": vp_ref,
                      "ghi_chu": "dmax = max |u − u*| trên mọi biến của nút mẫu; u* = Dykstra trên mọi hàng có khe ψ < 20h ∪ hàng đang vi phạm, dừng √(2·khe) ≤ 1e-9, rồi kiểm u* với MỌI hàng; PHI lưu float32 (sàn ~1e-6)"}
        TKC["luot"] = TKC["luot"][:luot_cu]; TKC["khong_hoi_tu"] = nk0
VBL = SH["vobao_le"]
def vobao(f_flat, k):
    f = f_flat.reshape(N1, N1, N1).cpu().numpy().astype(np.float32)
    if k in ctr:
        c = ctr[k]; r_ = rbk[k] + VBL[0] * h1 + VBL[1]; lo = np.clip(((c - r_ + A) / h1).astype(int), 0, N1); hi = np.clip(((c + r_ + A) / h1).astype(int) + SH["vobao_cong_o"], 0, N1)
    else: lo, hi = np.zeros(3, int), np.full(3, N1)
    sub = f[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]
    if (sub < 0).sum() < SH["vobao_o_toi_thieu"]: return None, 0
    v, fa, _, _ = marching_cubes(np.pad(sub, 1, constant_values=1.0), 0.0)
    return trimesh.Trimesh((v - 1 + lo) * h1 - A + h1 / 2, fa, process=True), int(ndi.label(sub < 0)[1])
VB, b0, VB0, b0t = {}, {}, {}, {}
for k in [VO_TEN] + bien:
    VB[k], b0[k] = vobao(psi_vo if k == VO_TEN else PHI[k], k)
    VB0[k], b0t[k] = vobao(psi_vo if k == VO_TEN else psi[k], k)
def luu(M_, ten):
    np.savez_compressed(f"{OUT}/{ten}", **{f"{k}__v": np.asarray(M_[k].vertices, np.float32) for k in M_ if M_[k] is not None},
                        **{f"{k}__f": np.asarray(M_[k].faces, np.int32) for k in M_ if M_[k] is not None})
T2 = ({}, {}, {}, [])
if ARG.tang2:
    with torch.no_grad():
        S2_ = nh.sc(); R2_ = nh.R()
        _cd = [(b.ten, b.loai, b.khuon, b.t.cpu().numpy(), b.R().cpu().numpy(), float(b.sc()), b.tam().cpu().numpy()) for b in lon]
    T2 = dung_tang2(_cd, VB)
    for k_, m_ in T2[0].items(): VB[k_] = m_; VB0[k_] = m_
luu(VB, "mesh_vo_bao.npz"); luu(VB0, "mesh_vo_bao_truoc.npz")
def sd(Pt, m): return igl.signed_distance(np.asarray(Pt, float), np.asarray(m.vertices, float), np.asarray(m.faces, np.int64))[0]
do = []
for lo_, a, b, can, src in ([] if ARG.khong_tu_do else mota):
    if lo_ == "giữ chỗ": continue
    if VB[a] is None or VB[b] is None: do.append({"loai": lo_, "a": a, "b": b, "can": can, "do": None, "dat": False, "dat_khong_chong": False, "nguon": src}); continue
    d = float((-sd(VB[a].vertices, VB[b])).min()) if lo_ == "chứa" else float(min(sd(VB[a].vertices, VB[b]).min(), sd(VB[b].vertices, VB[a]).min()))
    do.append({"loai": lo_, "a": a, "b": b, "can": can, "do": round(d, 3), "dat": bool(d >= can - SH["tol_do"]), "dat_khong_chong": bool(d >= -SH["tol_do"]), "nguon": src})
CH = DT["chiem"]
with torch.no_grad():
    vol = {k: (the_tich[k], float((PHI[k] < 0).sum()) * h1 ** 3) for k in bien}
    iou = {k: float(((psi[k] < 0) & (PHI[k] < 0)).sum()) / max(float(((psi[k] < 0) | (PHI[k] < 0)).sum()), 1.0) for k in bien}
    doi = {}
    for k in bien:
        band = psi[k].abs() < SH["doi_mat_dai"] * h1
        if band.any(): dd = (PHI[k][band] - psi[k][band]).abs(); doi[k] = [round(float(dd.mean()), 4), round(float(dd.max()), 4)]
        else: doi[k] = [None, None]
    kd = (psi_vo < -CH["khe_vo"]) & (cam_g[CH["vung_cam"]] > CH["khe_cam"]); nho_ = [k for k in bien if loai[k] in CH["loai"]]
    mn = torch.full((NX,), 1e9, device=DEV)
    for k in nho_: mn = torch.minimum(mn, psi[k])
    chiem = {"kha_dung_um3": float(kd.sum()) * h1 ** 3, "chiem": float(((mn < 0) & kd).sum()) / max(float(kd.sum()), 1.0),
             "chiem_co_khe": float(((mn < CH["khe_cam"] / SH["chiem_khe_chia"]) & kd).sum()) / max(float(kd.sum()), 1.0)}
KHE_ER_NHAN = None
if not ARG.khong_tu_do:
    for vt in PHU:
        q_ = ten_mot(VANH[vt]["quanh"])
        if VB.get(vt) is not None and VB.get(q_) is not None:
            d1_ = sd(VB[vt].vertices, VB[q_]); KHE_ER_NHAN = {"vung": vt, "quanh": q_, "khe_min_um": float(np.abs(d1_).min()), "dinh_trong_quanh": int((d1_ < -SH["tol_do"] - 1e-3).sum()), "so_dinh": int(len(d1_))}
tu_the = {b.ten: {"khuon": b.khuon, "t": b.t.cpu().numpy().round(7).tolist(), "s": round(float(b.sc()), 7), "R": b.R().cpu().numpy().round(9).tolist()} for b in lon}
for i in range(nh.n): tu_the[nh.ten[i]] = {"khuon": TEN_K[int(nh.tid[i])], "t": nh.t[i].cpu().numpy().round(7).tolist(), "s": round(float(S_[i]), 7), "R": R_[i].cpu().numpy().round(9).tolist()}
LKM = set(QH.get("lien_khoi_mien", []))
kq = {"ban": "botri_chung_v122e3", "ban_v121": V121, "dac_ta": ARG.dac_ta, "phien_ban_dac_ta": DT.get("phien_ban"), "seed": ARG.seed, "cach": ARG.cach, "thu_vien": ARG.thu_vien, "khoi_tao": ARG.khoi_tao,
      "so_luong": {l: SO[l] for l in LON + NHO}, "khong_chieu": ARG.khong_chieu,
      "trong_so": ARG.trong_so, "h1": h1, "N1": N1, "K": ARG.K, "so_luat": len(do), "dat": sum(x["dat"] for x in do), "dat_khong_chong": sum(x["dat_khong_chong"] for x in do),
      "diem_vi_pham_ban_dau": n_vp, "diem_chieu_theo_luot": so_luot, "diem_vi_pham_con_lai": con_vp, "giay_bo_tri": round(t_botri, 1), "giay_chieu": round(t_chieu, 1), "giay_tong": round(time.time() - t00, 1),
      "bo_nho_gb": round(torch.cuda.max_memory_allocated() / 2 ** 30, 2), "phat_cung_truoc_chieu": phat_truoc, "luat_mem": mem, "bao_cao_bo_tri": bao, "do": do,
      "the_tich_truoc_sau": vol, "b0": b0, "doi_mat_um": doi, "ti_le_chiem": chiem,
      "duoi_luoi": [k for k in bien if 2 * co[k] < SH["duoi_luoi_h"] * h1], "mat_bo_phan": [k for k in bien if vol[k][1] <= SH["mat_bo_phan_o"] * h1 ** 3],
      "con_lai_duoi_50pt": [k for k in bien if vol[k][1] < SH["con_lai_ti_le"] * vol[k][0]],
      "vo_manh": [k for k in bien if b0[k] > 1 and loai[k] not in LKM], "b0_truoc": b0t,
      "vo_do_chieu": [k for k in bien if b0[k] > b0t[k] and loai[k] not in LKM], "vo_do_luoi": [k for k in bien if b0t[k] > 1 and loai[k] not in LKM],
      "iou_truoc_sau": iou, "loai": {**{k: loai[k] for k in bien if loai[k] not in PHU}, **T2[1]}, "tu_the": tu_the,
      "ban_v12": {"chieu": ARG.chieu, "bien_mang_he_so": HS_M, "L_M": L_M, "beta": bu, "beta_M": DEL_M, "tang2": ARG.tang2, "khong_tu_do": ARG.khong_tu_do},
      "diem_khong_kha_thi_truoc_chieu": khong_kha_thi if ARG.chieu == "kkt" else None, "kha_thi": (khong_kha_thi == 0) if ARG.chieu == "kkt" else None, "thong_ke_chieu": tom_tat_chieu() if ARG.chieu == "kkt" else None,
      "con_vp_theo_loai_hang": con_theo_loai if ARG.chieu == "kkt" else None, "con_vp_max_um": con_max if ARG.chieu == "kkt" else None, "kiem_chieu": KIEM_CHIEU,
      "cha": T2[2], "tang2": T2[3], "luoi_noi_chat_nhan_khe_um": KHE_ER_NHAN, "ban_v122": {"hang_khong_giao_cap_mien_khe": [list(p) for p in N_C0], "tau": TAU_C0}, "ban_e3": E3_BC}
json.dump(kq, open(f"{OUT}/thong_tin.json", "w"), ensure_ascii=False, indent=1)
print("V12 kha_thi", khong_kha_thi == 0, "khong_kha_thi", khong_kha_thi, "kiem_chieu", KIEM_CHIEU, "thong_ke", json.dumps(kq["thong_ke_chieu"], ensure_ascii=False), "giay_dat", bao["giay_dat"], flush=True)
print(f"XONG_BOTRI {ARG.cach}/{ARG.khoi_tao} dat {kq['dat']}/{kq['so_luat']} con_vp {con_vp} <50% {kq['con_lai_duoi_50pt'][:5]} vo {kq['vo_manh'][:5]} "
      f"dat_lai {bao['dat_lai']} khong_dat {bao['khong_dat']} chiem {chiem['chiem']:.3f}/{chiem['chiem_co_khe']:.3f} t {kq['giay_bo_tri']}+{kq['giay_chieu']}s gb {kq['bo_nho_gb']}", flush=True)
