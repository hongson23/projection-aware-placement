import scipy.sparse as _sp
from scipy.sparse.csgraph import connected_components as _cc
from scipy.spatial import cKDTree
TG = {"bo_tri": round(time.time() - t00, 1)}; _tg0 = [time.time()]
def tg(ten):
    t_ = time.time(); TG[ten] = round(TG.get(ten, 0.0) + t_ - _tg0[0], 1); _tg0[0] = t_
t_cb = time.time(); h = h1; BETA = bu; LV = SH["loc_cap_h"]
RCH = QH["roi_chieu"]; S_C = RCH["s"]; MCAP = {frozenset(p) for p in RCH["mien_cap"]}; MCUNG = set(RCH["mien_cung_loai"])
KCL_C = RCH.get("khe_cung_loai", {}); KCAP_C = {frozenset((a_, b_)): s_ for a_, b_, s_ in RCH.get("khe_cap", [])}
def khe_c(la, lb):
    if la != lb and frozenset((la, lb)) in MCAP: return None
    if la == lb and la in MCUNG: return None
    if frozenset((la, lb)) in KCAP_C: return KCAP_C[frozenset((la, lb))]
    if la == lb and la in KCL_C: return KCL_C[la]
    return S_C
VPT = DT.get("vung_phu_thuoc")
if VPT is None:
    VPT = [{"ten": v["ten"], "quanh": v["quanh"], "day": v["day"], "trong": v["trong"], "khuyet": v.get("khuyet"), "cach_vo": CHUA[v["ten"]]["t"],
            "tru": "*", "cach_bo_phan": S_C, "mien_cach": [v["quanh"]], "nguon": v.get("nguon", "")} for v in VANH.values()]
QUANH = {v["quanh"]: v for v in VPT if v.get("quanh")}
BV = SH["beta_V"]

DS = []
with torch.no_grad():
    for b in lon: DS.append(dict(ten=b.ten, loai=b.loai, khuon=b.khuon, t=b.t.clone(), R=b.R().clone(), S=float(b.sc()), tid=b.tid))
    S_ = nh.sc(); R_ = nh.R()
    for i in range(nh.n): DS.append(dict(ten=nh.ten[i], loai=nh.loai[i], khuon=TEN_K[int(nh.tid[i])], t=nh.t[i].clone(), R=R_[i].clone(), S=float(S_[i]), tid=int(nh.tid[i])))
for d in DS:
    d["tam"] = ((SAU[d["tid"]] * d["S"]) @ d["R"].T + d["t"]).cpu().numpy(); d["co"] = float(RTR[d["tid"]] * d["S"]); d["rb"] = float(RB[d["tid"]] * d["S"])
    d["kd"] = KD.get(d["khuon"]) if KD.get(d["khuon"], {}).get("hinh") in HINH_GT else None
def psi_bp(d, P):
    loc = ((P - d["t"]) @ d["R"]) / d["S"]
    if d["kd"] is not None: return d["S"] * sdf_hinh(loc, d["kd"])
    return d["S"] * eval_at(loc, torch.full(loc.shape[:-1], d["tid"], device=DEV))
NP = len(DS); bien = [d["ten"] for d in DS]; ix = {k: i for i, k in enumerate(bien)}; loai = {d["ten"]: d["loai"] for d in DS}; DUM = NP
LCO = sorted({d["loai"] for d in DS})
Wl = {}
for la in LCO:
    ks = [khe_c(la, lb) for lb in LCO]; Wl[la] = max([x for x in ks if x is not None] + [S_C]) + BETA + SH["dia_phuong_h"] * h
    if la in QUANH: v = QUANH[la]; Wl[la] = max(Wl[la], (max(v["day"]) if isinstance(v["day"], list) else v["day"]) + v["cach_bo_phan"] + SH["le_quanh"])
if VO["hinh"] == "hop_bo_tron": EXT = NC - D_VO + 2 * h
else: EXT = np.full(3, A)
GLO = np.floor(-EXT / h).astype(int); GHI = np.ceil(EXT / h).astype(int)
for d in DS:
    m_ = d["rb"] + Wl[d["loai"]] + float(np.linalg.norm(d["tam"] - d["t"].cpu().numpy())); c_ = d["t"].cpu().numpy()
    lo_ = np.floor((c_ - m_) / h).astype(int); hi_ = np.ceil((c_ + m_) / h).astype(int)
    if ARG.luoi_day_du: lo_, hi_ = GLO.copy(), GHI.copy()
    d["lo"] = np.maximum(lo_, GLO); d["hi"] = np.minimum(hi_, GHI); d["shp"] = d["hi"] - d["lo"] + 1; d["n"] = int(d["shp"].prod())
TK = SH["truong_khoi"]
def diem_hop(d, a0, a1):
    j = torch.arange(a0, a1, device=DEV); sh = d["shp"]
    gi = torch.stack([j // (sh[1] * sh[2]), (j // sh[2]) % sh[1], j % sh[2]], 1) + T(d["lo"], torch.long)
    return gi, gi.to(F64) * h
with torch.no_grad():
    for d in DS:
        d["psi"] = torch.empty(d["n"], dtype=torch.float32, device=DEV)
        for a0 in range(0, d["n"], TK): _, X_ = diem_hop(d, a0, min(d["n"], a0 + TK)); d["psi"][a0:a0 + X_.shape[0]] = psi_bp(d, X_).float()
        d["PHI"] = d["psi"].clone()
tg("hop_con_psi"); print("HOP_CON", NP, "diem", sum(d["n"] for d in DS), "GB", round(sum(d["n"] for d in DS) * 8 / 2 ** 30, 2), round(time.time() - t_cb, 1), flush=True)
LO = np.array([d["lo"] for d in DS]); HI = np.array([d["hi"] for d in DS]); LOt = T(LO, torch.long); HIt = T(HI, torch.long)
def giao(a, b): lo_ = np.maximum(LO[a], LO[b]); hi_ = np.minimum(HI[a], HI[b]); return (lo_, hi_) if (hi_ >= lo_).all() else None
def lat(d, lo_, hi_, ten="psi"):
    a = lo_ - d["lo"]; b = hi_ - d["lo"] + 1; return d[ten].reshape(*d["shp"].tolist())[a[0]:b[0], a[1]:b[1], a[2]:b[2]]

rv, rc, rk, rval, rtam, rrb, mota = [], [], [], [], [], [], []
def them_hang(a, ca, b, cb, kind, val, tm=(0, 0, 0), rb_=0.0): rv.append([a, b]); rc.append([ca, cb]); rk.append(kind); rval.append(val); rtam.append(tm); rrb.append(rb_)
def ten_mot(l): return ten_ban(l, 0)
for k in bien:
    r = CHUA[loai[k]]
    if r["b"] == VO_TEN: them_hang(ix[k], 1, DUM, 0, 1, r["t"] + DEL_M); mota.append(("chứa", k, VO_TEN, r["t"], r["nguon"]))
    else: p_ = ten_mot(r["b"]); them_hang(ix[k], 1, ix[p_], -1, 0, r["t"] + bu); mota.append(("chứa", k, p_, r["t"], r["nguon"]))
ngoai = [k for k in bien if loai[k] in RCH["loai"] or loai[k] in KCL_C]
ctr = {d["ten"]: d["t"].cpu().numpy() for d in DS}; rbk = {d["ten"]: d["rb"] for d in DS}; co = {d["ten"]: d["co"] for d in DS}; tam = {d["ten"]: d["tam"] for d in DS}
TR = cKDTree(np.array([ctr[k] for k in ngoai])) if ngoai else None; RBmax = max(rbk.values())
with torch.no_grad():
    for i_, a in enumerate(ngoai):
        for j_ in sorted(TR.query_ball_point(ctr[a], 2 * RBmax + max(Wl.values()) + 3.0)):
            if j_ <= i_: continue
            b = ngoai[j_]; la, lb = loai[a], loai[b]
            if la == lb and la not in RCH["loai"] and la not in KCL_C: continue
            if (la not in RCH["loai"] or lb not in RCH["loai"]) and not (la == lb and la in KCL_C): continue
            s = khe_c(la, lb)
            if s is None: continue
            if np.linalg.norm(ctr[a] - ctr[b]) > rbk[a] + rbk[b] + s + bu + LV * h: continue
            gi = giao(ix[a], ix[b])
            if gi is None or float((lat(DS[ix[a]], *gi) + lat(DS[ix[b]], *gi)).min()) > s + bu + LV * h: continue
            them_hang(ix[a], 1, ix[b], 1, 0, s + bu); mota.append(("rời", a, b, s, RCH.get("nguon_cung_loai", RCH["nguon"]) if la == lb else RCH["nguon"]))
gc_ = SH["giu_cho"]
for k in bien:
    rr = min(max(gc_[0] * bu / gc_[1], gc_[2] * co[k]), gc_[3] * co[k]); them_hang(ix[k], -1, DUM, 0, 2, rr, tam[k], rbk[k] + gc_[4] * h); mota.append(("giữ chỗ", k, None, 0, ""))
print("SO_HANG", len(rv), "roi", sum(1 for x in mota if x[0] == "rời"), flush=True); tg("hang")
RRB = T(rrb); RV = torch.tensor(rv, device=DEV); RC = T(rc); RK = torch.tensor(rk, device=DEV); RVAL = T(rval); RTAM = T(np.array(rtam)); nr = RV.shape[0]
the_tich = {d["ten"]: float((d["psi"] < 0).sum()) * h ** 3 for d in DS}
w = torch.ones(NP + 1, dtype=F64, device=DEV)
if ARG.trong_so:
    for k in bien: w[ix[k]] = max(the_tich[k], SH["the_tich_san"]) ** (-1 / SH["trong_so_mu_mau"])
    w[:NP] = w[:NP] / w[:NP].min()
LANG = []; HANG = []
rv_np = np.array(rv)
for p in range(NP):
    nb = np.nonzero(((LO <= HI[p]) & (HI >= LO[p])).all(1))[0]; LANG.append(nb)
    ok = np.isin(rv_np[:, 0], nb) & (np.isin(rv_np[:, 1], nb) | (rv_np[:, 1] == DUM)); HANG.append(torch.tensor(np.nonzero(ok)[0], device=DEV))
DPH = SH["dia_phuong_h"]; TOLC = SH["tol_chieu"]; LAPK = SH["lap_kiem"]; RM_ = SH["RMAX"]; NGOAI_HOP = SH["ngoai_hop"]
def gom(gi, p, ten):
    N = gi.shape[0]; U = torch.full((N, NP + 1), NGOAI_HOP, dtype=F64, device=DEV); U[:, DUM] = 0.0; li = {}
    for q in LANG[p]:
        dq = DS[q]; m = ((gi >= LOt[q]) & (gi <= HIt[q])).all(1)
        if not m.any(): continue
        l3 = gi[m] - LOt[q]; sh = dq["shp"]; lf = (l3[:, 0] * sh[1] + l3[:, 1]) * sh[2] + l3[:, 2]
        U[m, q] = dq[ten][lf].double(); li[int(q)] = (m, lf)
    return U, li
def B_khoi(X_, hang):
    rkh = RK[hang]; B = RVAL[hang][None].expand(X_.shape[0], hang.numel()).clone(); m1 = rkh == 1
    if m1.any(): B[:, m1] = B[:, m1] + phi_vo(X_)[:, None]
    m2 = rkh == 2
    if m2.any():
        dd_ = torch.linalg.norm(X_[:, None, :] - RTAM[hang][m2][None], dim=2)
        B[:, m2] = torch.where(dd_ < RRB[hang][m2][None], -(dd_ - RVAL[hang][m2][None]), torch.full_like(dd_, -1e9))
    return B
def res_of(U, B, hang): return (U[:, RV[hang, 0]] * RC[hang, 0] + U[:, RV[hang, 1]] * RC[hang, 1]) - B
def chieu_diem_cb(Uv, Bv, sub, hang, RMAX):
    """Dykstra cục bộ (như botri2) cho các điểm có hàng hang; trả x và var_slot."""
    Nv = sub.shape[0]; kk = int(min(RMAX, int((sub < DPH * h).sum(1).max())))
    vals, rs_ = torch.topk(-sub, kk, dim=1); hop = (-vals) < DPH * h; rs = hang[rs_]
    vid = RV[rs]; cf = RC[rs] * hop[..., None]
    flat = vid.reshape(Nv, -1); sv, perm = torch.sort(flat, 1); moi = torch.ones_like(sv, dtype=torch.bool); moi[:, 1:] = sv[:, 1:] != sv[:, :-1]
    slot_s = torch.cumsum(moi.long(), 1) - 1; slot = torch.empty_like(slot_s).scatter_(1, perm, slot_s); Ks = int(slot_s.max()) + 1
    x = torch.zeros(Nv, Ks, dtype=F64, device=DEV).scatter_(1, slot, Uv.gather(1, flat))
    ws = torch.ones(Nv, Ks, dtype=F64, device=DEV).scatter_(1, slot, w[flat]); var_slot = torch.full((Nv, Ks), DUM, device=DEV).scatter_(1, slot, flat)
    sl = slot.reshape(Nv, kk, 2); Bs = Bv.gather(1, rs_); wsr = ws.gather(1, sl.reshape(Nv, -1)).reshape(Nv, kk, 2)
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
    return x, var_slot
def ghi_lai(x, var_slot, sel, li):
    for vk in var_slot.unique().tolist():
        if vk == DUM or vk not in li: continue
        m = var_slot == vk; pi_, si_ = m.nonzero(as_tuple=True); mq, lf = li[vk]
        pos = torch.full((mq.shape[0],), -1, dtype=torch.long, device=DEV); pos[mq.nonzero().flatten()] = torch.arange(int(mq.sum()), device=DEV)
        pp = pos[sel[pi_]]; ok = pp >= 0
        if (x[pi_[~ok], si_[~ok]] != NGOAI_HOP).any(): raise RuntimeError("biến không phủ điểm bị đổi giá trị")
        DS[vk]["PHI"][lf[pp[ok]]] = x[pi_[ok], si_[ok]].float()
n_vp = 0; so_luot = [0] * len(RM_); khong_kha_thi = 0; t0 = time.time()
KHc = max(SH["KH"][0], min(SH["KH"][1], int(SH["KH"][2] / (SH["KH"][3] * max(nr, 1)))))
CHI_VP = []
def dat_lai_psi(sel, li, cham):
    for q in cham:
        if q == DUM or q not in li: continue
        mq, lf = li[q]; pos = torch.full((mq.shape[0],), -1, dtype=torch.long, device=DEV); pos[mq.nonzero().flatten()] = torch.arange(int(mq.sum()), device=DEV)
        pp = pos[sel]; ok = pp >= 0; DS[q]["PHI"][lf[pp[ok]]] = DS[q]["psi"][lf[pp[ok]]]
if ARG.chieu == "cu":
    with torch.no_grad():
        for p in range(NP):
            d = DS[p]; hang = HANG[p]
            if not hang.numel(): continue
            for a0 in range(0, d["n"], KHc):
                gi, X_ = diem_hop(d, a0, min(d["n"], a0 + KHc))
                own = ((gi[:, None, :] >= LOt[None, :p]) & (gi[:, None, :] <= HIt[None, :p])).all(2).any(1) if p else torch.zeros(gi.shape[0], dtype=torch.bool, device=DEV)
                keep = (~own).nonzero().flatten()
                if not keep.numel(): continue
                gi, X_ = gi[keep], X_[keep]
                U, li = gom(gi, p, "psi"); B = B_khoi(X_, hang); res = res_of(U, B, hang); vi = (res < -TOLC).any(1)
                if not vi.any(): continue
                ub = torch.full_like(U, float("inf")); lb = torch.full_like(U, -float("inf")); rkh = RK[hang]
                for kind_, sgn in ((2, -1), (1, 1)):
                    mm = (rkh == kind_).nonzero().flatten()
                    if not mm.numel(): continue
                    vv = RV[hang[mm], 0]; Bm = B[:, mm]
                    if kind_ == 2: ub.scatter_reduce_(1, vv[None].expand(Bm.shape[0], -1), torch.where(Bm > -1e8, -Bm, torch.full_like(Bm, float("inf"))), "amin")
                    else: lb.scatter_reduce_(1, vv[None].expand(Bm.shape[0], -1), Bm, "amax")
                bad = (lb > ub + TOLC).any(1)
                m0 = ((rkh == 0) & (RC[hang, 1] > 0)).nonzero().flatten()
                if m0.numel(): bad |= ((ub[:, RV[hang[m0], 0]] + ub[:, RV[hang[m0], 1]]) < B[:, m0] - TOLC).any(1)
                khong_kha_thi += int(bad.sum())
                n_vp += int(vi.sum())
                if ARG.khong_chieu: continue
                sel = vi.nonzero().flatten(); x, vs = chieu_diem_cb(U[sel], B[sel], res[sel], hang, RM_[0]); ghi_lai(x, vs, sel, li); so_luot[0] += int(sel.numel())
                for luot in range(1, len(RM_)):
                    U2, li2 = gom(gi, p, "PHI"); res2 = res_of(U2, B, hang); v2 = (res2 < -TOLC).any(1)
                    if not v2.any(): break
                    sel = v2.nonzero().flatten(); x, vs = chieu_diem_cb(U2[sel], B[sel], res2[sel], hang, RM_[luot]); ghi_lai(x, vs, sel, li2); so_luot[luot] += int(sel.numel())

else:
  with torch.no_grad():
    for p in range(NP):
        d = DS[p]; hang = HANG[p]
        if not hang.numel(): continue
        for a0 in range(0, d["n"], KHc):
            gi, X_ = diem_hop(d, a0, min(d["n"], a0 + KHc))
            own = ((gi[:, None, :] >= LOt[None, :p]) & (gi[:, None, :] <= HIt[None, :p])).all(2).any(1) if p else torch.zeros(gi.shape[0], dtype=torch.bool, device=DEV)
            keep = (~own).nonzero().flatten()
            if not keep.numel(): continue
            gi, X_ = gi[keep], X_[keep]
            U, li = gom(gi, p, "psi"); B = B_khoi(X_, hang); res = res_of(U, B, hang); vi = (res < -TOLC).any(1)
            if not vi.any(): continue
            n_vp += int(vi.sum()); vn = torch.zeros_like(vi); vn[vi] = kha_thi_dong(B[vi], RV[hang], RC[hang], RK[hang], NP + 1); khong_kha_thi += int(vn.sum())
            capn = torch.where(vn, torch.full_like(vn, SH["lap_toi_da"], dtype=torch.long), torch.full_like(vn, ARG.lap_kkt, dtype=torch.long)); ep = torch.zeros_like(res, dtype=torch.bool)
            if ARG.kiem_chieu: CHI_VP.append((p, gi[vi].cpu()))
            if ARG.khong_chieu: continue
            cham = set()
            def _chieu(sel, sub_, RMAX_, ten_):
                x, sl, cf, Bs, hop, wsr, nrm, step_c, var_slot = chuan_bi_hang(U[sel], B[sel], sub_, RMAX_, hang)
                x, _, _ = dykstra_kkt(x, sl, cf, Bs, hop, wsr, nrm, step_c, ARG.lap_kkt, ARG.kkt_tol, ten_, cap=capn[sel])
                dat_lai_psi(sel, li, cham); ghi_lai(x, var_slot, sel, li); cham.update(var_slot.unique().tolist())
            sel = vi.nonzero().flatten(); _chieu(sel, res[sel], RM_[0], "1"); so_luot[0] += int(sel.numel())
            for luot in range(1, len(RM_)):
                U2, _ = gom(gi, p, "PHI"); res2 = res_of(U2, B, hang); v2 = (res2 < -SH["tol_con"]).any(1)
                if not v2.any(): break
                ep |= res2 < -SH["tol_con"]; sel = v2.nonzero().flatten(); sub2 = torch.where(ep[sel], torch.full_like(res2[sel], -1e9), res[sel])
                _chieu(sel, sub2, RM_[luot], str(luot + 1)); so_luot[luot] += int(sel.numel())
with torch.no_grad():
    torch.cuda.synchronize(); t_chieu = time.time() - t0
    con_vp = 0; con_max = 0.0
    for p in range(NP):
        d = DS[p]; hang = HANG[p]
        if not hang.numel(): continue
        for a0 in range(0, d["n"], KHc):
            gi, X_ = diem_hop(d, a0, min(d["n"], a0 + KHc))
            own = ((gi[:, None, :] >= LOt[None, :p]) & (gi[:, None, :] <= HIt[None, :p])).all(2).any(1) if p else torch.zeros(gi.shape[0], dtype=torch.bool, device=DEV)
            keep = (~own).nonzero().flatten()
            if not keep.numel(): continue
            U, _ = gom(gi[keep], p, "PHI"); rr_ = res_of(U, B_khoi(X_[keep], hang), hang); con_vp += int((rr_ < -SH["tol_con"]).any(1).sum()); con_max = max(con_max, float(-rr_.min()) if rr_.numel() else 0.0)
KIEM_CHIEU = None
if ARG.kiem_chieu and CHI_VP and not ARG.khong_chieu:
    with torch.no_grad():
        nn_ = [g_.shape[0] for _, g_ in CHI_VP]; tot_ = sum(nn_); gk_ = np.random.default_rng(ARG.seed); pick = set(gk_.choice(tot_, min(ARG.kiem_chieu, tot_), replace=False).tolist())
        dmax = 0.0; dW = 0.0; kk_max = 0; so_ = 0; vp_ref = 0; luot_cu = len(TKC["luot"]); nk0 = TKC["khong_hoi_tu"]; o_ = 0
        for (p, g_), n_ in zip(CHI_VP, nn_):
            ch_ = [i_ - o_ for i_ in range(o_, o_ + n_) if i_ in pick]; o_ += n_
            if not ch_: continue
            gi = g_[ch_].to(DEV); X_ = gi.to(F64) * h; hang = HANG[p]
            U, _ = gom(gi, p, "psi"); UP, _ = gom(gi, p, "PHI"); B = B_khoi(X_, hang); res = res_of(U, B, hang); resP = res_of(UP, B, hang)
            okn = ~kha_thi_dong(B, RV[hang], RC[hang], RK[hang], NP + 1)
            if not okn.any(): continue
            U, UP, B, res, resP = U[okn], UP[okn], B[okn], res[okn], resP[okn]
            sub = torch.where(resP < -TOLC, torch.full_like(resP, -1e9), res / 5)
            x, sl, cf, Bs, hop, wsr, nrm, step_c, var_slot = chuan_bi_hang(U, B, sub, 10 ** 6, hang); kk_max = max(kk_max, Bs.shape[1])
            x, _, _ = dykstra_kkt(x, sl, cf, Bs, hop, wsr, nrm, step_c, 20000, 1e-9, "tham_chieu")
            real = var_slot != DUM; cur = UP.gather(1, var_slot); dd = (x - cur).abs() * real; dmax = max(dmax, float(dd.max())); so_ += int(okn.sum())
            Ur = U.clone(); Ur.scatter_(1, var_slot, x); Ur[:, DUM] = 0.0; vp_ref += int((res_of(Ur, B, hang) < -SH["tol_con"]).any(1).sum())
            dW = max(dW, float(((x - cur) ** 2 * w[var_slot] * real).sum(1).sqrt().max()))
        KIEM_CHIEU = {"so_nut": so_, "dmax_um": dmax, "dmax_W": dW, "so_hang_toi_da": kk_max, "khong_hoi_tu_tham_chieu": TKC["khong_hoi_tu"] - nk0, "nut_tham_chieu_vi_pham_hang_nao_do": vp_ref,
                      "ghi_chu": "dmax = max |u − u*|; u* = Dykstra trên hàng có khe ψ < 20h ∪ đang vi phạm, √(2·khe) ≤ 1e-9, kiểm với mọi hàng của hộp; PHI float32"}
        TKC["luot"] = TKC["luot"][:luot_cu]; TKC["khong_hoi_tu"] = nk0
tg("chieu_dykstra"); print("CHIEU", n_vp, so_luot, "con", con_vp, "khong_kha_thi", khong_kha_thi, round(t_chieu, 1), flush=True)

def mc_hop(d, ten):
    f = d[ten].reshape(*d["shp"].tolist()).cpu().numpy().astype(np.float32)
    if (f < 0).sum() < SH["vobao_o_toi_thieu"]: return None, 0
    v, fa, _, _ = marching_cubes(np.pad(f, 1, constant_values=1.0), 0.0)
    return trimesh.Trimesh((v - 1 + d["lo"]) * h, fa, process=True), int(ndi.label(f < 0)[1])
VB, b0, VB0, b0t = {}, {}, {}, {}
for d in DS: VB[d["ten"]], b0[d["ten"]] = mc_hop(d, "PHI"); VB0[d["ten"]], b0t[d["ten"]] = mc_hop(d, "psi")
tg("mc_bo_phan_ran")
with torch.no_grad():
    vol = {d["ten"]: (the_tich[d["ten"]], float((d["PHI"] < 0).sum()) * h ** 3) for d in DS}
    iou = {d["ten"]: float(((d["psi"] < 0) & (d["PHI"] < 0)).sum()) / max(float(((d["psi"] < 0) | (d["PHI"] < 0)).sum()), 1.0) for d in DS}
    doi = {}
    for d in DS:
        band = d["psi"].abs() < SH["doi_mat_dai"] * h
        doi[d["ten"]] = [round(float((d["PHI"][band] - d["psi"][band]).abs().mean()), 4), round(float((d["PHI"][band] - d["psi"][band]).abs().max()), 4)] if band.any() else [None, None]
if os.environ.get("LUU_PHI"):
    np.savez_compressed(f"{OUT}/phi_hop.npz", **{f"{d['ten']}__phi": d["PHI"].cpu().numpy() for d in DS}, **{f"{d['ten']}__lo": d["lo"] for d in DS},
                        **{f"{d['ten']}__shp": d["shp"] for d in DS}, **{f"{d['ten']}__psi": d["psi"].cpu().numpy() for d in DS})

t1 = time.time()
for d in DS:
    W_ = Wl[d["loai"]]; D = torch.where(d["PHI"] <= -SH["D_trong"], torch.full_like(d["PHI"], -SH["D_trong"]), torch.full_like(d["PHI"], W_))
    m = ((d["psi"] < W_) & (d["PHI"] > -SH["D_trong"])).nonzero().flatten()
    if VB[d["ten"]] is not None and m.numel():
        for a0 in range(0, m.numel(), 4_000_000):
            mm = m[a0:a0 + 4_000_000]; sh = d["shp"]
            gi = torch.stack([mm // (sh[1] * sh[2]), (mm // sh[2]) % sh[1], mm % sh[2]], 1) + T(d["lo"], torch.long)
            dv = igl.signed_distance((gi.to(F64) * h).cpu().numpy(), np.asarray(VB[d["ten"]].vertices, float), np.asarray(VB[d["ten"]].faces, np.int64))[0]
            D[mm] = torch.as_tensor(np.minimum(dv, W_), dtype=torch.float32, device=DEV)
    elif VB[d["ten"]] is None: D = torch.full_like(d["PHI"], W_)
    d["D"] = D; del d["psi"]
torch.cuda.empty_cache(); print("D_MESH", round(time.time() - t1, 1), flush=True); tg("thong_ke_va_D_igl")
def tra_D(q, X_):
    d = DS[q]; u = X_ / h - T(d["lo"]); i0 = torch.floor(u).long(); f = u - i0; sh = T(d["shp"], torch.long); W_ = Wl[d["loai"]]
    ok = ((i0 >= 0) & (i0 + 1 < sh)).all(1); out = torch.full((X_.shape[0],), W_, dtype=F64, device=DEV)
    if not ok.any(): return out
    i0 = i0[ok]; f = f[ok]; Dv = d["D"]; acc = 0.0
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                ii = i0 + T([dx, dy, dz], torch.long); lf = (ii[:, 0] * sh[1] + ii[:, 1]) * sh[2] + ii[:, 2]
                wt = (f[:, 0] if dx else 1 - f[:, 0]) * (f[:, 1] if dy else 1 - f[:, 1]) * (f[:, 2] if dz else 1 - f[:, 2])
                acc = acc + wt * Dv[lf].double()
    out[ok] = acc; return out
def tra_D_luoi(q, gi):
    d = DS[q]; l3 = gi - LOt[q]; sh = d["shp"]; ok = ((l3 >= 0) & (l3 < T(sh, torch.long))).all(1)
    out = torch.full((gi.shape[0],), Wl[d["loai"]], dtype=F64, device=DEV)
    if ok.any(): l3 = l3[ok]; out[ok] = d["D"][(l3[:, 0] * sh[1] + l3[:, 1]) * sh[2] + l3[:, 2]].double()
    return out

OT = SH["o_vung"]
def tru_cua(v, da_dung):
    rigid = [q for q in range(NP)] if v.get("tru", "*") == "*" else [q for q in range(NP) if DS[q]["loai"] in v["tru"]]
    mien = set(v.get("mien_cach", [])); return [q for q in rigid if DS[q]["loai"] not in mien], [q for q in rigid if DS[q]["loai"] in mien], list(da_dung)
def xay_ong(v, day_q):
    """ống ngoại vi: đồ thị nút trên mặt lệch của vỏ; cạnh chỉ giữ khi mọi điểm mẫu cách bộ phận ≥ r + khe + β_V; giữ thành phần chạm vành."""
    o = v["ong_ngoai_vi"]; r_ = o["ban_kinh"]; lui = D_VO + o["do_sau_tam"]; bb = NC - lui; rr = RBO - lui; core = bb - rr
    g = np.random.default_rng(ARG.seed * SH["hat_ong"][0] + SH["hat_ong"][1])
    ar = np.array([core[1] * core[2], core[0] * core[2], core[0] * core[1]]); Atot = 8 * ar.sum(); nc_ = int(Atot / (o["khoang_nut"] ** 2) * 4)
    ax_ = g.choice(3, size=nc_, p=ar / ar.sum()); sg = np.where(g.uniform(size=nc_) < 0.5, -1.0, 1.0); U_ = g.uniform(-1, 1, (nc_, 3)) * bb
    U_[np.arange(nc_), ax_] = sg * bb[ax_]; cc = np.clip(U_, -core, core); dv_ = U_ - cc; nn_ = np.linalg.norm(dv_, axis=1, keepdims=True)
    P0 = cc + rr * dv_ / np.maximum(nn_, 1e-12)
    keep = []; tr_ = None
    for i in range(nc_):
        if keep and np.min(np.linalg.norm(P0[keep] - P0[i], axis=1)) < o["khoang_nut"]: continue
        keep.append(i)
    Pn = P0[keep]
    def chieu_mat(Q):
        Qt = T(Q)
        for _ in range(2):
            Qt.requires_grad_(True); f = sdf_hop(Qt, lui); gr, = torch.autograd.grad(f.sum(), Qt); Qt = (Qt - f[:, None] * gr / gr.norm(dim=1, keepdim=True).clamp_min(1e-12)).detach()
        return Qt
    tru_, mien_, _ = tru_cua(v, [])
    def thoang(Qt):
        ok = torch.ones(Qt.shape[0], dtype=torch.bool, device=DEV)
        for q in tru_: ok &= tra_D(q, Qt) >= r_ + v["cach_bo_phan"] + BV
        for q in mien_: ok &= tra_D(q, Qt) >= r_
        return ok
    Qn = chieu_mat(Pn)
    with torch.no_grad(): okn = thoang(Qn).cpu().numpy()
    Pn = Qn.cpu().numpy()[okn]; tree = cKDTree(Pn); cap_ = np.array(sorted(tree.query_pairs(o["canh_toi_da"])), dtype=int).reshape(-1, 2)
    ns = int(o["canh_toi_da"] / h) + 2; tt = T(np.linspace(0, 1, ns))[None, :, None]; canh = []
    for a0 in range(0, len(cap_), 2000):
        cp = cap_[a0:a0 + 2000]; Q = T(Pn[cp[:, 0]])[:, None] * (1 - tt) + T(Pn[cp[:, 1]])[:, None] * tt
        Qp = chieu_mat(Q.reshape(-1, 3).cpu().numpy())
        with torch.no_grad(): okc = thoang(Qp).reshape(len(cp), ns).all(1).cpu().numpy()
        Qp = Qp.reshape(len(cp), ns, 3).cpu().numpy(); st = max(1, int(round(o["buoc_mau"] / (o["canh_toi_da"] / (ns - 1)))))
        for k_ in np.nonzero(okc)[0]:
            L_ = np.linalg.norm(Pn[cp[k_, 1]] - Pn[cp[k_, 0]]); m_ = max(2, int(np.ceil(L_ / o["buoc_mau"])) + 1)
            canh.append((int(cp[k_, 0]), int(cp[k_, 1]), chieu_mat(np.linspace(Pn[cp[k_, 0]], Pn[cp[k_, 1]], m_)).cpu().numpy(), None))
    n_ = len(Pn); A_ = _sp.coo_matrix((np.ones(len(canh)), ([c_[0] for c_ in canh], [c_[1] for c_ in canh])), shape=(n_, n_))
    ncomp, lab = _cc(A_, directed=False); qq = ix[ten_mot(v["quanh"])] if v.get("quanh") else None
    cham = set()
    if qq is not None:
        with torch.no_grad():
            for i, j, Q, _ in canh:
                if bool((tra_D(qq, T(Q)) <= day_q + r_).any()): cham.add(lab[i])
    giu = cham if cham else ({int(np.bincount(lab[[c_[0] for c_ in canh]]).argmax())} if canh else set())
    caps = []
    for i, j, Q, _ in canh:
        if lab[i] not in giu: continue
        for a_ in range(len(Q) - 1): caps.append(np.concatenate([Q[a_], Q[a_ + 1]]))
    info = {"so_nut": int(n_), "so_canh": len(canh), "so_thanh_phan": int(ncomp), "giu_thanh_phan": len(giu), "cham_vanh": bool(cham), "so_vien_nang": len(caps), "ban_kinh": r_}
    return (T(np.array(caps)) if caps else None), r_, info
def truong_vung(v, X_, gi, ong, r_ong, day_q, truoc):
    """trường vùng phụ thuộc tại điểm X_ (điểm mạng gi)."""
    f = phi_vo(X_) + v["cach_vo"]
    core = None
    if v.get("quanh"):
        qq = ix[ten_mot(v["quanh"])]; Dq = tra_D_luoi(qq, gi); core = torch.maximum(-Dq + v.get("trong", 0.0), Dq - day_q)
        if v.get("khuyet"):
            c_ = DS[qq]["t"]; u_ = DS[ix[ten_mot(v["khuyet"]["huong_toi"])]]["t"] - c_; u_ = u_ / u_.norm().clamp_min(1e-9)
            rr_ = torch.linalg.norm(X_ - c_, dim=-1); core = torch.maximum(core, ((X_ - c_) @ u_) - v["khuyet"]["cos"] * rr_)
    if ong is not None:
        LE_ = SH["le_ong"]; lo_ = X_.min(0).values - LE_; hi_ = X_.max(0).values + LE_
        a_, b_ = ong[:, :3], ong[:, 3:]; mk = ((torch.minimum(a_, b_) <= hi_) & (torch.maximum(a_, b_) >= lo_)).all(1)
        if mk.any():
            a_, b_ = a_[mk], b_[mk]; dmin = torch.full((X_.shape[0],), LE_, dtype=F64, device=DEV)
            for c0 in range(0, a_.shape[0], 64):
                aa, ab = a_[c0:c0 + 64], b_[c0:c0 + 64] - a_[c0:c0 + 64]; pa = X_[:, None, :] - aa[None]
                tt = ((pa * ab[None]).sum(-1) / (ab * ab).sum(-1).clamp_min(1e-18)[None]).clamp(0, 1)
                dmin = torch.minimum(dmin, torch.linalg.norm(pa - tt[..., None] * ab[None], dim=-1).min(1).values)
            co_ = dmin - r_ong; core = co_ if core is None else torch.minimum(core, co_)
        else: co_ = torch.full((X_.shape[0],), LE_ - r_ong, dtype=F64, device=DEV); core = co_ if core is None else torch.minimum(core, co_)
    if core is not None: f = torch.maximum(f, core)
    tru_, mien_, prior = tru_cua(v, truoc)
    gmin = gi.min(0).values; gmax = gi.max(0).values
    for q in tru_:
        if (LOt[q] > gmax).any() or (HIt[q] < gmin).any(): continue
        f = torch.maximum(f, v["cach_bo_phan"] + BV - tra_D_luoi(q, gi))
    for q in mien_:
        if (LOt[q] > gmax).any() or (HIt[q] < gmin).any(): continue
        f = torch.maximum(f, -tra_D_luoi(q, gi))
    for pv in prior: f = torch.maximum(f, v["cach_bo_phan"] + BV - pv["truong"](X_, gi))
    return f
def mc_vung(fn, L_lip):
    """marching cubes theo ô OT³ trên miền; chỉ ô có |f(tâm)| ≤ L·nửa đường chéo + lề."""
    nt = np.ceil((GHI - GLO) / OT).astype(int); hd = math.sqrt(3) * OT * h / 2 * L_lip + SH["le_o"]
    c3 = np.stack(np.meshgrid(*[np.arange(n_) for n_ in nt], indexing="ij"), -1).reshape(-1, 3); gc_ = GLO + c3 * OT + OT // 2
    with torch.no_grad():
        fc = torch.cat([fn(T(gc_[a0:a0 + 4096]) * h, T(gc_[a0:a0 + 4096], torch.long)) for a0 in range(0, len(gc_), 4096)]).cpu().numpy()
    chon = c3[np.abs(fc) <= hd]; trong = int((fc < -hd).sum()); Vs, Fs, nv_ = [], [], 0; tg("vung_chon_o")
    ar = np.arange(OT + 1)
    for c in chon:
        g0 = GLO + c * OT; gi = torch.as_tensor(np.stack(np.meshgrid(ar, ar, ar, indexing="ij"), -1).reshape(-1, 3) + g0, device=DEV)
        with torch.no_grad(): f = fn(gi.to(F64) * h, gi).reshape(OT + 1, OT + 1, OT + 1).cpu().numpy().astype(np.float32)
        tg("vung_truong_gpu")
        if f.min() >= 0 or f.max() <= 0: continue
        v_, fa_, _, _ = marching_cubes(f, 0.0, allow_degenerate=False); Vs.append((v_ + g0) * h); Fs.append(fa_ + nv_); nv_ += len(v_); tg("vung_mc_cpu")
    if not Vs: return None, {"so_o": int(len(chon)), "o_trong": trong}
    m = trimesh.Trimesh(np.concatenate(Vs), np.concatenate(Fs), process=True)
    ps = m.split(only_watertight=False); giu = [p for p in ps if abs(p.volume) >= SH["the_tich_manh_min"]]
    m = trimesh.util.concatenate(giu) if giu else None; tg("vung_ghep_trimesh")
    return m, {"so_o": int(len(chon)), "o_trong": trong, "bo_manh_vun": len(ps) - len(giu)}
def so_khoi(m):
    if m is None: return 0
    return int(sum(p.volume > 0 for p in m.split(only_watertight=False)))
VUNG = []; _gv = np.random.default_rng(ARG.seed * SH["hat_vung"][0] + SH["hat_vung"][1]); t2 = time.time()
for v in VPT:
    day_q = None
    if v.get("quanh"): day_q = float(_gv.uniform(*v["day"])) if isinstance(v["day"], list) else v["day"]
    tg("khac"); ong, r_ong, info = (xay_ong(v, day_q if day_q is not None else 0.0) if v.get("ong_ngoai_vi") else (None, 0.0, None)); tg("vung_do_thi_ong")
    truoc = list(VUNG)
    fn = (lambda X_, gi, v=v, ong=ong, r_ong=r_ong, day_q=day_q, truoc=truoc: truong_vung(v, X_, gi, ong, r_ong, day_q, truoc))
    L_lip = 1.0 + (v["khuyet"]["cos"] if v.get("khuyet") else 0.0) + (SH["lip_vo_them"] if VO["hinh"] != "hop_bo_tron" else 0.0)
    m, inf2 = mc_vung(fn, L_lip)
    VUNG.append({"ten": v["ten"], "truong": fn, "mesh": m, "day_quanh": day_q, "ong": info, "o": inf2, "v": v})
    VB[v["ten"]] = m; VB0[v["ten"]] = None; loai[v["ten"]] = v["ten"]
    print("VUNG", v["ten"], "mat", 0 if m is None else len(m.faces), inf2, info, round(time.time() - t2, 1), flush=True)

def mesh_hop(lui, buoc):
    b = NC - lui; rr = RBO - lui; core = b - rr; V_, F_ = [], []
    for ax in range(3):
        for sg in (-1.0, 1.0):
            o1, o2 = [k for k in range(3) if k != ax]; n1 = max(2, int(2 * b[o1] / buoc) + 1); n2 = max(2, int(2 * b[o2] / buoc) + 1)
            u, w_ = np.meshgrid(np.linspace(-b[o1], b[o1], n1), np.linspace(-b[o2], b[o2], n2), indexing="ij")
            P = np.zeros((n1, n2, 3)); P[..., ax] = sg * b[ax]; P[..., o1] = u; P[..., o2] = w_; P = P.reshape(-1, 3)
            cc = np.clip(P, -core, core); dv_ = P - cc; P = cc + rr * dv_ / np.maximum(np.linalg.norm(dv_, axis=1, keepdims=True), 1e-12)
            idx = np.arange(n1 * n2).reshape(n1, n2); q = np.stack([idx[:-1, :-1], idx[1:, :-1], idx[1:, 1:], idx[:-1, 1:]], -1).reshape(-1, 4)
            f = np.concatenate([q[:, [0, 1, 2]], q[:, [0, 2, 3]]]); f = f if (sg > 0) == (ax != 1) else f[:, ::-1]
            F_.append(f + sum(len(x) for x in V_)); V_.append(P)
    m = trimesh.Trimesh(np.concatenate(V_), np.concatenate(F_), process=True); m.merge_vertices(); trimesh.repair.fix_normals(m); return m
if VO["hinh"] == "hop_bo_tron":
    for lp in LOP_VO:
        mo = mesh_hop(lp["ngoai"], SH["buoc_vo"]); mi = mesh_hop(lp["ngoai"] + lp["day"], SH["buoc_vo"]); mi.invert()
        VB[lp["ten"]] = trimesh.util.concatenate([mo, mi]); VB0[lp["ten"]] = VB[lp["ten"]]
else:
    axv = np.arange(GLO[0], GHI[0] + 1) * h; Xv = torch.stack(torch.meshgrid(T(axv), T(axv), T(axv), indexing="ij"), -1)
    with torch.no_grad(): fv = phi_vo(Xv.reshape(-1, 3)).reshape(Xv.shape[:3]).cpu().numpy().astype(np.float32)
    v_, fa_, _, _ = marching_cubes(np.pad(fv, 1, constant_values=1.0), 0.0); VB[VO_TEN] = trimesh.Trimesh((v_ - 1 + GLO) * h, fa_, process=True); VB0[VO_TEN] = VB[VO_TEN]
def luu(M_, ten):
    np.savez_compressed(f"{OUT}/{ten}", **{f"{k}__v": np.asarray(M_[k].vertices, np.float32) for k in M_ if M_[k] is not None},
                        **{f"{k}__f": np.asarray(M_[k].faces, np.int32) for k in M_ if M_[k] is not None})
T2 = ({}, {}, {}, [])
if ARG.tang2:
    T2 = dung_tang2([(d["ten"], d["loai"], d["khuon"], d["t"].cpu().numpy(), d["R"].cpu().numpy(), d["S"], d["tam"]) for d in DS], VB)
    for k_, m_ in T2[0].items(): VB[k_] = m_; VB0[k_] = m_
tg("khac"); luu(VB, "mesh_vo_bao.npz"); luu(VB0, "mesh_vo_bao_truoc.npz"); tg("luu_npz")

def sd(Pt, m): return igl.signed_distance(np.asarray(Pt, float), np.asarray(m.vertices, float), np.asarray(m.faces, np.int64))[0]
def kc_vo(m):
    with torch.no_grad(): return -phi_vo(T(np.asarray(m.vertices))).cpu().numpy()
def sd_cat(Pt, m, lo_, hi_, le=1.5):
    fm = ((m.vertices >= lo_ - le) & (m.vertices <= hi_ + le)).all(1)[m.faces].any(1)
    if not fm.any(): return np.full(len(Pt), 9.0)
    d = igl.signed_distance(np.asarray(Pt, float), np.asarray(m.vertices, float), np.asarray(m.faces[fm], np.int64))[0]
    return np.where(np.abs(d) < le - 0.5, d, np.abs(d))
def kc_hai(ma, mb):
    lo_ = np.maximum(ma.bounds[0], mb.bounds[0]) - 1.0; hi_ = np.minimum(ma.bounds[1], mb.bounds[1]) + 1.0
    if (lo_ > hi_).any(): return 9.0
    va = ma.vertices[((ma.vertices >= lo_) & (ma.vertices <= hi_)).all(1)]; vb = mb.vertices[((mb.vertices >= lo_) & (mb.vertices <= hi_)).all(1)]
    x = [sd_cat(va, mb, lo_, hi_).min()] if len(va) else []; x += [sd_cat(vb, ma, lo_, hi_).min()] if len(vb) else []
    return float(min(x)) if x else 9.0
do = []
def ghi_do(lo_, a, b, can, src, d):
    do.append({"loai": lo_, "a": a, "b": b, "can": can, "do": None if d is None else round(d, 3), "dat": d is not None and bool(d >= can - SH["tol_do"]),
               "dat_khong_chong": d is not None and bool(d >= -SH["tol_do"]), "nguon": src})
for lo_, a, b, can, src in ([] if ARG.khong_tu_do else mota):
    if lo_ == "giữ chỗ": continue
    if VB.get(a) is None or (b != VO_TEN and VB.get(b) is None): ghi_do(lo_, a, b, can, src, None); continue
    if lo_ == "chứa": d = float(kc_vo(VB[a]).min()) if b == VO_TEN else float((-sd(VB[a].vertices, VB[b])).min())
    else: d = kc_hai(VB[a], VB[b])
    ghi_do(lo_, a, b, can, src, d)
tg("do_hang_ran_igl")
for V in ([] if ARG.khong_tu_do else VUNG):
    v = V["v"]; m = V["mesh"]
    ghi_do("chứa", v["ten"], VO_TEN, v["cach_vo"], v.get("nguon", ""), None if m is None else float(kc_vo(m).min()))
    tru_, mien_, prior = tru_cua(v, [x for x in VUNG if x is not V and VUNG.index(x) < VUNG.index(V)])
    for q in tru_:
        if VB[bien[q]] is None or m is None: continue
        d = kc_hai(m, VB[bien[q]])
        if d < v["cach_bo_phan"] + 1.0: ghi_do("rời", v["ten"], bien[q], v["cach_bo_phan"], v.get("nguon", ""), d)
    for pv in prior:
        if pv["mesh"] is not None and m is not None: ghi_do("rời", v["ten"], pv["ten"], v["cach_bo_phan"], v.get("nguon", ""), kc_hai(m, pv["mesh"]))
def V_hop2(lui):
    a_ = NC - lui - (RBO - lui); rr = RBO - lui
    return 8 * a_.prod() + 8 * (a_[0] * a_[1] + a_[1] * a_[2] + a_[0] * a_[2]) * rr + 2 * math.pi * a_.sum() * rr ** 2 + 4 / 3 * math.pi * rr ** 3
V_TRONG = V_hop2(D_VO) if VO["hinh"] == "hop_bo_tron" else None
vung_do = {}
for V in VUNG:
    v = V["v"]; m = V["mesh"]; e = {"so_mat": 0 if m is None else int(len(m.faces)), "so_khoi": so_khoi(m), "o": V["o"], "ong": V["ong"], "day_quanh": V["day_quanh"]}
    if m is not None:
        e["the_tich"] = float(m.volume); e["kin"] = bool(m.is_watertight)
        if V_TRONG: e["ti_le_the_tich"] = e["the_tich"] / V_TRONG
        if v.get("quanh") and VB.get(ten_mot(v["quanh"])) is not None and not ARG.khong_tu_do:
            mq = VB[ten_mot(v["quanh"])]; bx = ((m.vertices >= mq.bounds[0] - 0.5) & (m.vertices <= mq.bounds[1] + 0.5)).all(1)
            dv = np.full(len(m.vertices), 9.0); dv[bx] = np.abs(sd(m.vertices[bx], mq)) if bx.any() else 9.0; gan_ = (dv < SH["tol_cham"])[m.faces].all(1)
            e["dien_tich_cham_um2"] = float(m.area_faces[gan_].sum())
    e["lien_khoi_can"] = v.get("lien_khoi"); e["ti_le_can"] = v.get("ti_le_the_tich"); vung_do[v["ten"]] = e
tg("do_vung_igl")
mem = {}
def kc_den(a, b_list):
    out = 9e9
    for b in b_list:
        if b == VO_TEN: out = min(out, float(kc_vo(VB[a]).min())); continue
        cand = sorted([q for q in range(NP) if DS[q]["loai"] == b and bien[q] != a and VB[bien[q]] is not None], key=lambda q: np.linalg.norm(ctr[bien[q]] - ctr[a]))
        for q in cand[:SH["gan_so_lang_gieng"]]: out = min(out, float(sd(VB[a].vertices, VB[bien[q]]).min()))
    return out
gan_do = []
for r in ([] if ARG.khong_tu_do else QH.get("gan", []) + QH.get("gan_nho", [])):
    bs = r["b"] if isinstance(r["b"], list) else [r["b"]]
    for k in bien:
        if loai[k] != r["a"] or VB[k] is None: continue
        dk = kc_den(k, bs); gan_do.append({"a": k, "b": bs, "toi_da": r["toi_da"], "do": round(dk, 3), "dat": bool(dk <= r["toi_da"] + SH["tol_do"]), "nguon": r.get("nguon", "")})
huong_do = []
for r in ([] if ARG.khong_tu_do else QH.get("huong", [])):
    for q in range(NP):
        if DS[q]["loai"] != r["a"]: continue
        with torch.no_grad():
            cs = float((DS[q]["R"][:, 2] * phap_vo(DS[q]["t"][None])[0]).sum().abs())
        g_ = math.degrees(math.acos(min(1.0, cs))); huong_do.append({"a": bien[q], "goc_do": round(g_, 2), "dat": bool(g_ <= r["lech_toi_da_do"]), "nguon": r.get("nguon", "")})
tg("do_gan_huong")
tt_luat = {}
def gop(src, dat):
    for m_ in str(src).replace("(", " ").replace(")", " ").replace("+", " ").replace(",", " ").replace("/", " ").split():
        if any(ch.isdigit() for ch in m_) or "." in m_:
            e = tt_luat.setdefault(m_, {"dat": 0, "truot": 0}); e["dat" if dat else "truot"] += 1
for x in do: gop(x["nguon"], x["dat"])
for x in gan_do: gop(x["nguon"], x["dat"])
for x in huong_do: gop(x["nguon"], x["dat"])
for V in VUNG:
    v = V["v"]; e = vung_do[v["ten"]]
    for key, ok in (("nguon_ti_le", e.get("ti_le_the_tich") is not None and v.get("ti_le_the_tich") and v["ti_le_the_tich"][0] <= e["ti_le_the_tich"] <= v["ti_le_the_tich"][1]),
                    ("nguon_lien_khoi", (not v.get("lien_khoi")) or e["so_khoi"] == 1),
                    ("nguon_cham", (v.get("cham_toi_thieu_um2") is None) or e.get("dien_tich_cham_um2", 0.0) >= v["cham_toi_thieu_um2"])):
        if v.get(key): gop(v[key], bool(ok))
json.dump({"seed": ARG.seed, "luat": tt_luat, "ghi_chu": "bộ dựng tự đo trên mesh; chấm cuối bằng bộ kiểm độc lập"}, open(f"{OUT}/trang_thai_luat.json", "w"), ensure_ascii=False, indent=1)
tu_the = {d["ten"]: {"khuon": d["khuon"], "t": d["t"].cpu().numpy().round(7).tolist(), "s": round(d["S"], 7), "R": d["R"].cpu().numpy().round(9).tolist()} for d in DS}
LKM = set(QH.get("lien_khoi_mien", []))
kq = {"ban": "botri_chung_v121/cuc_bo", "ban_v121": V121, "dac_ta": ARG.dac_ta, "phien_ban_dac_ta": DT.get("phien_ban"), "seed": ARG.seed, "cach": ARG.cach, "thu_vien": ARG.thu_vien, "khoi_tao": ARG.khoi_tao,
      "so_luong": {l: SO[l] for l in LON + NHO}, "khong_chieu": ARG.khong_chieu, "luoi": LUOI, "luoi_day_du": ARG.luoi_day_du,
      "trong_so": ARG.trong_so, "h1": h, "beta": BETA, "beta_V": BV, "K": ARG.K, "so_luat": len(do), "dat": sum(x["dat"] for x in do), "dat_khong_chong": sum(x["dat_khong_chong"] for x in do),
      "lop_vo": LOP_VO if VO["hinh"] == "hop_bo_tron" else None,
      "diem_vi_pham_ban_dau": n_vp, "diem_chieu_theo_luot": so_luot, "diem_vi_pham_con_lai": con_vp, "diem_khong_kha_thi_truoc_chieu": khong_kha_thi,
      "giay_bo_tri": round(t_botri, 1), "giay_chieu": round(t_chieu, 1), "giay_tong": round(time.time() - t00, 1),
      "bo_nho_gb": round(torch.cuda.max_memory_allocated() / 2 ** 30, 2), "phat_cung_truoc_chieu": phat_truoc, "bao_cao_bo_tri": bao, "do": do, "gan": gan_do, "huong": huong_do,
      "vung_phu_thuoc": vung_do, "the_tich_truoc_sau": vol, "b0": b0, "doi_mat_um": doi,
      "duoi_luoi": [k for k in bien if 2 * co[k] < SH["duoi_luoi_h"] * h], "mat_bo_phan": [k for k in bien if vol[k][1] <= SH["mat_bo_phan_o"] * h ** 3],
      "con_lai_duoi_50pt": [k for k in bien if vol[k][1] < SH["con_lai_ti_le"] * vol[k][0]],
      "vo_manh": [k for k in bien if b0[k] > 1 and loai[k] not in LKM], "b0_truoc": b0t, "iou_truoc_sau": iou,
      "loai": {**{k: loai[k] for k in bien}, **T2[1]}, "tu_the": tu_the, "cha": T2[2], "tang2": T2[3],
      "ban_v12": {"chieu": ARG.chieu, "bien_mang_he_so": HS_M, "L_M": L_M, "beta": BETA, "beta_M": DEL_M, "tang2": ARG.tang2, "khong_tu_do": ARG.khong_tu_do},
      "kha_thi": khong_kha_thi == 0, "thong_ke_chieu": tom_tat_chieu() if ARG.chieu == "kkt" else None, "con_vp_max_um": con_max, "kiem_chieu": KIEM_CHIEU, "thoi_gian_buoc_giay": TG, "quy_uoc_mesh": "lớp vỏ = hai mặt (ngoài + trong đảo chiều); vùng phụ thuộc = mesh kín"}
tg("khac"); kq["thoi_gian_buoc_giay"] = TG; json.dump(kq, open(f"{OUT}/thong_tin.json", "w"), ensure_ascii=False, indent=1); print("THOI_GIAN", TG, flush=True)
vv = {k: (round(e.get("ti_le_the_tich", -1), 3), e["so_khoi"]) for k, e in vung_do.items()}
print("V12 kha_thi", khong_kha_thi == 0, "khong_kha_thi", khong_kha_thi, "kiem_chieu", KIEM_CHIEU, "thong_ke", json.dumps(kq["thong_ke_chieu"], ensure_ascii=False), "giay_dat", bao.get("giay_dat"), flush=True)
print(f"XONG_BOTRI {ARG.cach}/{ARG.khoi_tao} cuc_bo dat {kq['dat']}/{kq['so_luat']} con_vp {con_vp} kkt {khong_kha_thi} vung {vv} "
      f"gan {sum(x['dat'] for x in gan_do)}/{len(gan_do)} huong {sum(x['dat'] for x in huong_do)}/{len(huong_do)} t {kq['giay_bo_tri']}+{kq['giay_chieu']}s tong {kq['giay_tong']} gb {kq['bo_nho_gb']}", flush=True)
