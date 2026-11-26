"""cau_noi_v122.py — cầu nối cellPACK ↔ bộ dựng v12.2 (18/11/2026; sinh từ cau_noi_v12b.py, md5 cau_noi_v12b ghi trong CAU_NOI_GOC).
KHÔNG sửa botri_chung_v122.py (md5 8abb7b1c… kiểm trước khi chạy). So với cau_noi_v12b.py CHỈ đổi:
  (1) bộ dựng v12.2 (md5, đường dẫn mặc định, nhãn 'ban');
  (2) chốt tệp SDF màng của lượt chính động vật (sdf_75005fbc7d88018f.npz, md5 ee14b688…): kiểm khoá băm + tệp + md5 TRƯỚC khi chạy
      phần A (để bộ dựng không bao giờ dựng lại SDF), kiểm lại sau phần A (V121 của bộ dựng);
  (3) ghi thêm vào đầu ra: V121 (sdf_tep, sdf_md5, adam_lr…), beta, beta_M, L_M, md5 cầu nối, md5 đặc tả, h_luoi_ngan.
Cách cắt A/L/C, khối L, cây cầu, dựng mesh ngăn, nhập, chọn: giữ nguyên văn cau_noi_v12b.py.
Chạy với cwd = ch3_tebao_dep (đường dẫn tương đối của bộ dựng, đặc tả, khuôn, SDF như lượt chính).

Cách làm: đọc nguyên văn mã bộ dựng đã đóng băng, cắt theo ba mốc văn bản cố định, thực thi:
  phần A  = đầu tệp → trước 't_dat0 = …' (đặc tả, khuôn, atlas, số lượng, khởi tạo, làm nóng)   — nguyên văn
  phần L  = khối 'bộ phận lớn thử–loại' của nhánh ngau_nhien (vòng lam_lai_nhom)                — nguyên văn
  phần C  = từ '_, lon, nh, k_chon, pl, pn = tot' → hết tệp (PROJECT, mesh, thong_tin.json)        — nguyên văn
Chỉ phần giữa (vòng sắp đặt K cấu hình) được thay bằng mã ở đây.

Hai chế độ (cùng --seed, cùng đối số bộ dựng sau '--'):
  xuat : A → khoi_tao(0) (khuôn, cỡ, khởi tạo CHUNG mọi nhánh ở k = 0) → L (bộ phận lớn đặt bằng thử–loại, như nhánh thử–loại)
         → ghi đầu vào cellPACK: lưới vỏ chứa (màng lùi T_mang + β_M − a), vùng cấm (ER ∪ trung thể ∪ Golgi nở S + β − a),
         lòng nhân cho nhân con (nhân lùi T_nhancon + β − a_nc), cây cầu bao từng bộ phận nhỏ / nhân con (cỡ đã rút, nở a).
  nhap : A → khoi_tao(0) → L (cùng số ngẫu nhiên ⇒ cùng kết quả) → gán tư thế cellPACK cho bộ phận nhỏ + nhân con, bỏ bộ phận
         cellPACK không đặt được → phat_day_du → C (PROJECT + mesh + thong_tin.json) → thêm khối 'cellpack' và mục tu_the
         {khuon, s, t: null, khong_dat: true} cho bộ phận không đặt được.
--thiet_bi cpu: chỉ để chạy thử trên CPU (thay DEV = "cuda" trong bản sao trên bộ nhớ, giả các hàm torch.cuda); lượt thật dùng cuda.

Dùng: python cau_noi_v122.py xuat|chon|nhap --thu_muc D [--ket_qua_cp tu_the_cp.json] [--thiet_bi cpu|cuda] -- <đối số botri_chung_v122.py>
"""
import argparse, hashlib, json, math, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("che_do", choices=["xuat", "chon", "nhap"]); ap.add_argument("--thu_muc", required=True)
ap.add_argument("--ket_qua_cp", default=None); ap.add_argument("--k", type=int, default=0, help="cấu hình k (khuôn, cỡ, khởi tạo chung với cấu hình k của mọi nhánh)"); ap.add_argument("--thiet_bi", default="cuda", choices=["cpu", "cuda"])
ap.add_argument("--bo_dung", default="bo_dung_chung/ma3/botri_chung_v122e3.py"); ap.add_argument("--giay_dat", type=float, default=None, help="nhap: giờ đặt tổng (Σ_k) do bộ điều phối cộng")
ap.add_argument("--h_luoi_ngan", type=float, default=0.25, help="bước lưới (µm) dựng mesh ngăn cho cellPACK")
ap.add_argument("--dien_tich_cay", type=int, default=4000, help="số điểm mặt cho bước tham lam (bước vá + kiểm dùng MỌI điểm mặt và điểm trong)")
ap.add_argument("--delta_cay", type=float, default=0.05, help="cầu phủ lồi khỏi mặt khuôn ≤ delta (đơn vị khuôn)")
if "--" not in sys.argv: sys.exit("cần '--' trước đối số bộ dựng")
_i = sys.argv.index("--"); A_ = ap.parse_args(sys.argv[1:_i]); ARGV_BD = sys.argv[_i + 1:]
MD5_DONG_BANG = "6e3beba8706763433586bfbda1a9151f"
CAU_NOI_GOC = "cau_noi_v12b.py"
SDF_KHOA = "75005fbc7d88018f"; SDF_TEP = "bo_dung_chung/v121/sdf_mang/sdf_75005fbc7d88018f.npz"; SDF_MD5 = "ee14b688e6225311735c3ebfbeecc902"
NHAN = "cau_noi_cp_e3"
MD5_CAU_NOI = hashlib.md5(open(os.path.abspath(__file__), "rb").read()).hexdigest()
src = open(A_.bo_dung, encoding="utf-8").read()
m5 = hashlib.md5(src.encode("utf-8")).hexdigest()
if m5 != MD5_DONG_BANG: sys.exit(f"md5 bộ dựng {m5} ≠ {MD5_DONG_BANG}: dừng")
MOC_A = 't_dat0 = time.time(); t_day_du = 0.0\n'
MOC_C = '_, lon, nh, k_chon, pl, pn = tot; HT = {b.loai: b for b in lon if b.loai in DON}\n'
MOC_L0 = '        if ARG.cach == "ngau_nhien":   # bộ phận lớn thử–loại tuần tự; rồi kiểm chung cả nhóm, hỏng thì làm lại cả nhóm\n'
MOC_L1 = '        thu_loai_nho(nh, lon, g, bao)\n'
for m_ in (MOC_A, MOC_C, MOC_L0, MOC_L1): assert src.count(m_) == 1, m_[:40]


def _kiem_sdf_truoc():
    """Chốt SDF màng TRƯỚC phần A: tính lại khoá băm như bộ dựng v12.2 (dòng _key) từ đặc tả + N1; đòi khoá = SDF_KHOA, tệp có sẵn,
    md5 tệp = SDF_MD5. Nếu không, bộ dựng sẽ dựng SDF mới (nặng) => dừng ngay."""
    ap0 = argparse.ArgumentParser(add_help=False); ap0.add_argument("--dac_ta", default="dac_ta/dong_vat_botri2.json"); ap0.add_argument("--N1", type=int, default=128)
    a0, _ = ap0.parse_known_args(ARGV_BD)
    DT0 = json.load(open(a0.dac_ta)); TS = {**DT0["tham_so_khong_thuoc_luat"], **DT0.get("tham_so_so_hoc", {})}
    def giai(x):
        if isinstance(x, str) and x.startswith("@"):
            ten, _, chi = x[1:].partition("#"); v = giai(TS[ten]["gia_tri"]); return v[int(chi)] if chi else v
        if isinstance(x, list): return [giai(v) for v in x]
        if isinstance(x, dict): return {k: giai(v) for k, v in x.items()}
        return x
    VO = giai(DT0["vo"]); assert VO["hinh"] == "cau_nhieu", "cầu nối v12.2 chỉ cho tế bào động vật (vỏ cau_nhieu)"
    HMC = 0.04; ext = VO["ban_kinh"] + VO["nhieu"]["bien_do"] * VO["nhieu"]["so_song"] / VO["nhieu"]["chia"] + 3 * HMC
    key = hashlib.md5(json.dumps({"vo": VO, "A": VO["nua_hop"], "N1": a0.N1, "HMC": HMC, "ext": ext, "ma": "v121a"}, sort_keys=True).encode()).hexdigest()[:16]
    if key != SDF_KHOA: sys.exit(f"khoá SDF màng {key} ≠ {SDF_KHOA} (đặc tả/N1 khác lượt chính): dừng")
    thu = os.environ.get("SDF_MANG_THU_MUC", "bo_dung_chung/v121/sdf_mang"); tep = f"{thu}/sdf_{key}.npz"
    if os.path.normpath(tep) != os.path.normpath(SDF_TEP): sys.exit(f"tệp SDF {tep} ≠ {SDF_TEP}: chạy với cwd = ch3_tebao_dep, không đặt SDF_MANG_THU_MUC")
    if not os.path.isfile(tep): sys.exit(f"thiếu {tep} (cwd sai?): dừng, KHÔNG để bộ dựng dựng lại SDF")
    m = hashlib.md5(open(tep, "rb").read()).hexdigest()
    if m != SDF_MD5: sys.exit(f"md5 SDF {m} ≠ {SDF_MD5}: dừng")
    return {"dac_ta": a0.dac_ta, "md5_dac_ta": hashlib.md5(open(a0.dac_ta, "rb").read()).hexdigest(), "N1": a0.N1, "sdf_khoa": key}


KIEM_TRUOC = _kiem_sdf_truoc()
pA = src[:src.index(MOC_A)]; pC = src[src.index(MOC_C):]
pL = src[src.index(MOC_L0) + len(MOC_L0):src.index(MOC_L1)]
pL = "\n".join(l[12:] if l.startswith(" " * 12) else l for l in pL.splitlines()) + "\n"
import torch
if A_.thiet_bi == "cpu":
    assert pA.count('DEV = "cuda"') == 1; pA = pA.replace('DEV = "cuda"', 'DEV = "cpu"', 1)
    torch.cuda.synchronize = lambda *a, **k: None; torch.cuda.max_memory_allocated = lambda *a, **k: 0
    torch.cuda.get_rng_state = lambda *a, **k: torch.get_rng_state(); torch.cuda.set_rng_state = lambda *a, **k: None
os.makedirs(A_.thu_muc, exist_ok=True)
sys.argv = [A_.bo_dung] + ARGV_BD
if "--cach" not in ARGV_BD: sys.argv += ["--cach", "ngau_nhien"]
NS = {"__name__": "__bo_dung__", "__file__": os.path.abspath(A_.bo_dung)}
exec(compile(pA, A_.bo_dung + "#A", "exec"), NS)
assert NS["ARG"].cach == "ngau_nhien", "cầu nối chỉ dùng --cach ngau_nhien (khối bộ phận lớn của thử–loại)"
assert NS["ARG"].eval_dung_tu_dau and not NS["ARG"].tu_the_tu, "cau_noi_cp_e3: bắt buộc --eval_dung_tu_dau (bộ dựng đã sửa, như N1)"
assert not NS["CUCBO"], "cầu nối chỉ hỗ trợ lưới toàn phần (tế bào động vật)"
V121 = NS["V121"]; assert V121.get("sdf_tep") == SDF_TEP and V121.get("sdf_md5") == SDF_MD5, ("SDF màng của phần A lệch lượt chính", V121.get("sdf_tep"), V121.get("sdf_md5"))
assert NS["L_M"] == 1.0, "v12.2: L_M phải = 1 (SDF đúng)"
exec(compile(f"lon, nh, g = khoi_tao({int(A_.k)}); HT = {{b.loai: b for b in lon if b.loai in DON}}\n", "cau_noi#k", "exec"), NS)
_tL = time.time(); exec(compile(pL, A_.bo_dung + "#L", "exec"), NS); GIAY_L = time.time() - _tL
import numpy as np
F64 = torch.float64; DEV = NS["DEV"]; T = NS["T"]
lon, nh = NS["lon"], NS["nh"]; SH = NS["SH"]; DEL = NS["DEL"]; DEL_M = NS["DEL_M"]
CHUA = NS["CHUA"]; T_NHO = NS["T_NHO"]; S_NN = NS["S_NN"]; S_NL = NS["S_NL"]; S_CN = NS["S_CN"]; CAM = NS["CAM_NHO"]
KCL_NN = NS["KCL_NN"]; KCL_LON = NS["KCL_LON"]
S_NHO = max([S_NN] + list(KCL_NN.values())); a_nho = (S_NHO + DEL) / 2
L_NC = [b for b in lon if b.loai == "nhan_con"]
a_nc = (KCL_LON.get("nhan_con", 0.0) + DEL) / 2
TH = {"a_nho": a_nho, "a_nhan_con": a_nc, "beta": DEL, "beta_M": DEL_M, "L_M": NS["L_M"], "S_nho_chung": S_NHO,
      "muc_vo_luat": -(CHUA["ti_the"]["t"] + DEL_M), "muc_cam_luat": S_CN + DEL, "muc_doi_tac_lon_luat": S_NL + DEL,
      "muc_long_nhan_luat": -(CHUA["nhan_con"]["t"] + DEL), "doi_vo_um": +a_nho, "doi_cam_um": -a_nho, "doi_long_nhan_um": +a_nc}
assert CHUA["ti_the"]["t"] == CHUA["lysosome"]["t"] == T_NHO


def cay_cau(tid, delta, n_mat):
    """cây cầu 2 tầng phủ khuôn tid (khung đơn vị): tầng 1 = cầu bao tâm gốc; tầng 2 = cầu tâm tại nút atlas trong khuôn, bán kính
    |φ(c)| + delta, chọn tham lam trên ứng viên, rồi VÁ: điểm mặt/điểm trong chưa phủ lấy cầu tại nút trong gần nhất (luôn phủ vì
    khoảng tới nút gần nhất ≤ nửa đường chéo ô = 0,043 < delta). Mọi cầu ⊂ khuôn ⊕ B_delta. Kiểm phủ trên MỌI điểm mặt (MATS)
    và mọi điểm trong (TRONGS) mà bộ dựng dùng."""
    from scipy.spatial import cKDTree
    grid = NS["KHUON"][NS["TEN_K"][tid]]["grid"][0, 0].permute(2, 1, 0).double().cpu().numpy()
    G = NS["G"]; GB = NS["GB"]; ax = np.linspace(-GB, GB, G); hg = ax[1] - ax[0]
    assert delta > hg * math.sqrt(3) / 2 + 1e-9, "delta phải > nửa đường chéo ô atlas"
    P = np.stack(np.meshgrid(ax, ax, ax, indexing="ij"), -1).reshape(-1, 3); f = grid.reshape(-1)
    ung_all = np.nonzero(f < 0)[0]; rng = np.random.default_rng(12345)
    ung = ung_all if len(ung_all) <= 3000 else np.unique(np.concatenate([ung_all[np.argsort(f[ung_all])[:200]], rng.choice(ung_all, 2800, replace=False)]))
    mat_all = np.concatenate([NS["MATS"][tid].cpu().numpy(), NS["TRONGS"][tid].cpu().numpy()])
    mat = mat_all[rng.permutation(len(mat_all))[:n_mat]] if n_mat < len(mat_all) else mat_all
    C = P[ung]; rC = -f[ung] + delta
    tr = cKDTree(mat); phu = [np.array(tr.query_ball_point(C[j], rC[j]), int) for j in range(len(C))]
    con = np.ones(len(mat), bool); chon = []
    while con.any():
        gain = np.array([con[p].sum() if len(p) else 0 for p in phu]); j = int(gain.argmax())
        if gain[j] == 0: break
        chon.append(j); con[phu[j]] = False
    tam = list(C[chon]); r = list(rC[chon])
    def chua_phu(Q):
        d = np.full(len(Q), np.inf)
        for c_, r_ in zip(tam, r): d = np.minimum(d, np.linalg.norm(Q - c_, axis=1) - r_)
        return d > 0
    tri = cKDTree(P[ung_all]); n_va = 0
    for _ in range(50):
        cp_ = chua_phu(mat_all)
        if not cp_.any(): break
        _, ii = tri.query(mat_all[cp_]); ii = np.unique(ii)
        for q_ in ii: tam.append(P[ung_all[q_]]); r.append(-f[ung_all[q_]] + delta); n_va += 1
    tam = np.array(tam); r = np.array(r); con_lai = int(chua_phu(mat_all).sum())
    assert con_lai == 0, f"cây cầu {NS['TEN_K'][tid]} còn {con_lai} điểm chưa phủ"
    R0 = float((np.linalg.norm(tam, axis=1) + r).max())
    return {"tam": tam.tolist(), "r": r.tolist(), "R0": R0, "so_cau": len(r), "so_cau_va": n_va, "so_diem_kiem_phu": int(len(mat_all)), "phu_ti_le": 1.0}


def luoi_tu_truong(f, lo, hi, h, muc, doi, ten):
    """(1) mesh M mặt mức trường {f = muc} (mức luật của bộ dựng); (2) SDF THẬT tới M (igl, dấu theo số vòng quấn) trên cùng lưới;
    (3) mesh cuối = {sdf_M = doi} (doi > 0: nở a µm; doi < 0: co a µm). Giữ mọi thành phần; đòi mesh kín."""
    from skimage.measure import marching_cubes
    import trimesh, igl
    ax = [np.arange(lo[i], hi[i] + h / 2, h) for i in range(3)]; sh = [len(a) for a in ax]
    Q = np.stack(np.meshgrid(*ax, indexing="ij"), -1).reshape(-1, 3); vals = np.empty(len(Q))
    with torch.no_grad():
        for a0 in range(0, len(Q), 500000): vals[a0:a0 + 500000] = f(T(Q[a0:a0 + 500000])).cpu().numpy()
    def mc(V):
        v, fa, _, _ = marching_cubes(np.pad(V.reshape(sh), 1, constant_values=1.0).astype(np.float32), 0.0)
        m = trimesh.Trimesh((v - 1) * h + np.array(lo), fa, process=True)
        if m.volume < 0: m.invert()
        return m
    M = mc(vals - muc); assert M.is_watertight, f"{ten}: mesh mức luật không kín"
    sd = igl.signed_distance(Q, np.asarray(M.vertices, float), np.asarray(M.faces, np.int64), sign_type=igl.SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER)[0]
    m = mc(sd - doi); kin = bool(m.is_watertight); assert kin, f"{ten}: mesh ngăn không kín"
    m.export(os.path.join(A_.thu_muc, ten + ".obj"))
    return {"tep": ten + ".obj", "so_dinh": int(len(m.vertices)), "kin": kin, "the_tich_um3": float(m.volume), "so_thanh_phan": len(m.split(only_watertight=False)),
            "the_tich_muc_luat_um3": float(M.volume), "doi_um": doi}


if A_.che_do == "xuat":
    t0 = time.time(); HT = NS["HT"]
    Av = NS["A"]; lo = [-Av] * 3; hi = [Av] * 3; h = A_.h_luoi_ngan
    out = {"ban": NHAN, "k": A_.k, "md5_bo_dung": m5, "md5_cau_noi": MD5_CAU_NOI, "argv_bo_dung": sys.argv[1:], "seed": NS["ARG"].seed, "nguong": TH, "ngan": {}, "bo_phan": [], "lon_co_dinh": {},
           "ban_v121": {k_: V121[k_] for k_ in V121}, "kiem_truoc": KIEM_TRUOC, "thiet_bi_xuat": A_.thiet_bi, "h_luoi_ngan": h}
    out["ngan"]["vo_chua"] = luoi_tu_truong(NS["phi_vo"], lo, hi, h, TH["muc_vo_luat"], TH["doi_vo_um"], "vo_chua")
    dt = [b for b in lon if b.loai in NS["RNL"]["doi_tac"]]
    def f_cam(P):
        v = NS["phi_cam"](P, NS["VANH"][CAM]) - TH["muc_cam_luat"]
        for b in dt: v = torch.minimum(v, b.phi(P) - TH["muc_doi_tac_lon_luat"])
        return v
    out["ngan"]["vung_cam"] = luoi_tu_truong(f_cam, lo, hi, h, 0.0, TH["doi_cam_um"], "vung_cam")
    q = HT["nhan"]; rq = q.rb() + 0.5; c = q.t.cpu().numpy()
    out["ngan"]["long_nhan"] = luoi_tu_truong(q.phi, (c - rq).tolist(), (c + rq).tolist(), h / 2, TH["muc_long_nhan_luat"], TH["doi_long_nhan_um"], "long_nhan")
    for b in lon:
        if b.loai != "nhan_con": out["lon_co_dinh"][b.ten] = {"loai": b.loai, "khuon": b.khuon, "t": b.t.cpu().numpy().tolist(), "R": b.R().cpu().numpy().tolist(), "s": float(b.sc())}
    cay = {}
    def them(ten, loai, tid, s, a, ngan):
        if tid not in cay: cay[tid] = cay_cau(tid, A_.delta_cay, A_.dien_tich_cay)
        cc = cay[tid]; tam = (np.array(cc["tam"]) * s).tolist(); r = (np.array(cc["r"]) * s + a).tolist()
        out["bo_phan"].append({"ten": ten, "loai": loai, "khuon": NS["TEN_K"][tid], "s": s, "ngan": ngan, "tam": tam, "r": r, "R0": cc["R0"] * s + a,
                               "lo_khoi_khuon_um": A_.delta_cay * s + a})
    S_ = nh.sc()
    for i in range(nh.n): them(nh.ten[i], nh.loai[i], int(nh.tid[i]), float(S_[i]), a_nho, "te_bao_chat")
    for b in L_NC: them(b.ten, b.loai, b.tid, float(b.sc()), a_nc, "long_nhan")
    out["cay_khuon"] = {NS["TEN_K"][k]: {kk: v[kk] for kk in ("so_cau", "so_cau_va", "so_diem_kiem_phu", "phu_ti_le", "R0")} for k, v in cay.items()}
    out["giay_khoi_L"] = round(GIAY_L, 3)
    out["hop"] = [lo, hi]; out["giay_xuat"] = round(time.time() - t0, 2)
    json.dump(out, open(os.path.join(A_.thu_muc, "dau_vao_cellpack.json"), "w"), ensure_ascii=False, indent=1)
    print("XUAT_XONG", len(out["bo_phan"]), "bo_phan", json.dumps(out["ngan"], ensure_ascii=False), json.dumps(out["cay_khuon"]), flush=True)
    sys.exit(0)

kq = json.load(open(A_.ket_qua_cp)); dat = kq["dat"]
dv = json.load(open(os.path.join(A_.thu_muc, "dau_vao_cellpack.json")))
Nhom = NS["Nhom"]; TEN_K = NS["TEN_K"]
def r6(R): R = np.asarray(R, float); return np.concatenate([R[:, 0], R[:, 1]])
assert dv["k"] == A_.k and dv["seed"] == NS["ARG"].seed and dv["argv_bo_dung"] == sys.argv[1:], "không khớp lúc xuất (k/seed/đối số)"
assert kq.get("seed_cp") is not None
import hashlib as _hl
assert kq.get("md5_dau_vao") == _hl.md5(open(os.path.join(A_.thu_muc, "dau_vao_cellpack.json"), "rb").read()).hexdigest(), "tu_the_cp.json không sinh từ dau_vao_cellpack.json này"
if A_.che_do == "nhap" and A_.thiet_bi == "cuda": assert A_.giay_dat is not None, "lượt thật phải truyền --giay_dat = Σ_k (pack_grid + khối L + chấm chọn)"
_bp = {b["ten"]: b for b in dv["bo_phan"]}; S0 = nh.sc()
for i in range(nh.n):
    e = _bp[nh.ten[i]]; assert e["khuon"] == TEN_K[int(nh.tid[i])] and abs(e["s"] - float(S0[i])) < 1e-9, ("lệch khuôn/cỡ", nh.ten[i])
lech_L = 0.0
for b in lon:
    if b.loai == "nhan_con":
        e = _bp[b.ten]; assert e["khuon"] == b.khuon and abs(e["s"] - float(b.sc())) < 1e-9; continue
    e = dv["lon_co_dinh"][b.ten]; assert e["khuon"] == b.khuon and abs(e["s"] - float(b.sc())) < 1e-9
    lech_L = max(lech_L, float(np.abs(b.t.cpu().numpy() - np.array(e["t"])).max()), float(np.abs(b.R().cpu().numpy() - np.array(e["R"])).max()))
    b.t = T(e["t"]); b.r = T(r6(e["R"]))
print("LECH_KHOI_L_so_voi_luc_xuat", lech_L, flush=True)
khong_dat = {}
ds = []
Rn = nh.R().cpu().numpy(); Sn = nh.sc()
for i in range(nh.n):
    ten = nh.ten[i]
    if ten not in dat: khong_dat[ten] = {"khuon": TEN_K[int(nh.tid[i])], "s": round(float(Sn[i]), 7), "t": None, "R": None, "khong_dat": True, "loai": nh.loai[i]}; continue
    R = np.asarray(dat[ten]["R"], float); assert abs(np.linalg.det(R) - 1) < 1e-6, (ten, np.linalg.det(R))
    ds.append((ten, nh.loai[i], TEN_K[int(nh.tid[i])], np.asarray(dat[ten]["t"], float), r6(R), float(nh.z[i])))
nh2 = Nhom(ds) if ds else None
assert nh2 is not None, "cellPACK không đặt được bộ phận nhỏ nào"
lon2 = []
for b in lon:
    if b.loai == "nhan_con":
        if b.ten not in dat: khong_dat[b.ten] = {"khuon": b.khuon, "s": round(float(b.sc()), 7), "t": None, "R": None, "khong_dat": True, "loai": b.loai}; continue
        R = np.asarray(dat[b.ten]["R"], float); assert abs(np.linalg.det(R) - 1) < 1e-6
        b.t = T(dat[b.ten]["t"]); b.r = T(r6(R))
    lon2.append(b)
NS["lon"], NS["nh"] = lon2, nh2; NS["HT"] = {b.loai: b for b in lon2 if b.loai in NS["DON"]}
NS["_GIAY_CP"] = float(A_.giay_dat) if A_.giay_dat is not None else kq["giay_dat"] + dv.get("giay_khoi_L", 0.0)
_tp = time.time()
exec(compile(f"pl, pn = phat_day_du(lon, nh); tot = (sum(pl.values()) + sum(pn.values()), lon, nh, {int(A_.k)}, pl, pn)\n"
             "GIAY_DAT = float(_GIAY_CP)\n", "cau_noi#nhap", "exec"), NS)
GIAY_CHAM = time.time() - _tp
NS["bao"]["khong_dat_cellpack_vang"] = len(khong_dat)
_pt = {**NS["pl"], **NS["pn"]}; NG_ = SH["nguong_phat"]
print("PHAT_TRUOC_CHIEU cellpack: bo_phan phat > nguong", sum(v > NG_ for v in _pt.values()), "/", len(_pt), "tong", round(sum(_pt.values()), 6),
      "max", max(_pt, key=_pt.get), round(max(_pt.values()), 6), flush=True)
if A_.che_do == "chon":
    json.dump({"k": A_.k, "so_khong_dat": len(khong_dat), "phat_cung_day_du": float(NS["tot"][0]), "so_bo_phan_phat_tren_nguong": sum(v > NG_ for v in _pt.values()),
               "giay_dat_cp": kq["giay_dat"], "giay_khoi_L": dv.get("giay_khoi_L"), "giay_cham_chon": round(GIAY_CHAM, 3), "lech_khoi_L": lech_L}, open(os.path.join(A_.thu_muc, f"chon_k{A_.k}.json"), "w"), indent=1)
    print("CHON_XONG", A_.k, len(khong_dat), float(NS["tot"][0]), flush=True); sys.exit(0)
NS["bao"]["cellpack"] = {k: kq[k] for k in kq if k != "dat"}
NS["bao"]["cellpack"]["bo_phan_phat_tren_nguong"] = sorted(k for k, v in _pt.items() if v > NG_); NS["bao"]["cellpack"]["khong_dat"] = sorted(khong_dat)
NS["bao"]["cellpack"]["so_dat"] = len(dat); NS["bao"]["cellpack"]["so_yeu_cau"] = len(dat) + len(khong_dat); NS["bao"]["cellpack"]["lech_khoi_L"] = lech_L
exec(compile('E3_BC = {"tu_the_tu": ARG.tu_the_tu, "eval_dung": ARG.eval_dung, "eval_dung_tu_dau": ARG.eval_dung_tu_dau, "G2": E3["G2"], "GB2": E3["GB2"]}\n'
             'if ARG.eval_dung_tu_dau: E3_BC.update({"kiem_xa": E3["kiem_xa"], "giay_xay_xa": E3["giay_xay"]})\n', "cau_noi#E3_BC", "exec"), NS)
exec(compile(pC, A_.bo_dung + "#C", "exec"), NS)
tt_ = os.path.join(NS["OUT"], "thong_tin.json"); j = json.load(open(tt_))
j["tu_the"].update(khong_dat); j["cach"] = "cellpack"; j["cau_noi"] = {"ban": NHAN, "md5_bo_dung": m5, "md5_cau_noi": MD5_CAU_NOI, "thiet_bi": A_.thiet_bi, "k": A_.k,
                                                                       "sdf_md5": V121["sdf_md5"], "thiet_bi_xuat": dv.get("thiet_bi_xuat"), "md5_dau_vao": kq.get("md5_dau_vao")}
json.dump(j, open(tt_, "w"), ensure_ascii=False, indent=1)
print("NHAP_XONG", NS["OUT"], "khong_dat", len(khong_dat), flush=True)
