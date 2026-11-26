# -*- coding: utf-8 -*-
"""giam_doc_quyen.py (22/11/2026, bản 2 sau rà soát) — bằng chứng ĐỘC QUYỀN cho lượt đo giờ (chỉ đọc /proc, /sys, nvidia-smi).
Ba chế độ:
 (1) --yen <cpu> --gpu <g|-> --giay S --cua_so N [--cho_toi_da T] [--nen ra.json]
     Cổng: N cửa sổ S giây LIÊN TIẾP mà GPU g không có tiến trình tính toán và luồng của NGƯỜI DÙNG (uid của ta, ngoài cây ta) trên vùng <cpu>
     dùng ≤ 0,02 lõi và uid khác ngoài danh sách trắng ≤ 0,05 lõi. Tiến trình HỆ THỐNG trong danh sách trắng (uid khác: falco, tailscaled, kworker…) không ghim được ⇒ không chặn cổng, chỉ đo: với --nen, ghi tải nền
     mỗi lõi (bận user+nice+system) của các cửa sổ yên ⇒ p50/p95 theo lõi ⇒ ngưỡng từng việc. Thoát 0 khi đạt; 3 khi quá T giây.
 (2) --pid P --gpu g --loi L --igl I --du_tru D --ra f.jsonl --tom f.json [--nen nen.json] [--chu_ky 2]
     Lấy mẫu tới khi P thoát. Mỗi mẫu: quét /proc lấy cây ta, rồi đọc /proc/stat và tick cây ta SÁT NHAU. Tick ta = u+s+cu+cs của mọi tiến trình
     sống trong cây (cu, cs giữ tick của con đã thoát); con của ta biến mất ⇒ trừ tick đã tính của nó khỏi lớp của cha (tránh tính hai lần).
     Lớp theo affinity: L (lõi việc), I (igl/nhập), B (bao: lõi bộ điều khiển). Lạ theo vùng = Σ bận − Σ tick ta, CỘNG CÓ DẤU, chỉ cắt về 0 ở cuối.
     Anh em SMT của L ∪ I: ta không bao giờ chạy ở đó ⇒ mọi tick bận là lạ. Gán chủ (người dùng / hệ thống) bằng quét luồng mỗi 10 s.
     DOC_QUYEN (bản 3b) ⇔ GPU g không có tiến trình lạ ở mọi mẫu ∧ tính toán người dùng ≤ 0,01 lõi ∧ uid khác ngoài danh sách trắng ≤ 0,03 lõi
                 ∧ lạ/lõi trên L ≤ nền_p95(L) + 0,02 ∧ mọi tiến trình ta có tick đều ghim trong L ∪ I ∪ B. Igl/anh em vượt nền + 0,05: chỉ cờ "vuot_nen".
     Chế độ (1) có --bo_qua_cmd <regex>: bỏ qua tiến trình có dòng lệnh khớp (bước chờ xả CPU trong khi GPU còn chạy việc tdg).
"""
import os, sys, time, json, subprocess, argparse
HZ = os.sysconf("SC_CLK_TCK"); UID = os.getuid()
ap = argparse.ArgumentParser()
ap.add_argument("--yen"); ap.add_argument("--giay", type=float, default=60); ap.add_argument("--cua_so", type=int, default=1); ap.add_argument("--cho_toi_da", type=float, default=1e12)
ap.add_argument("--nen"); ap.add_argument("--pid", type=int); ap.add_argument("--gpu", default="-"); ap.add_argument("--loi", default=""); ap.add_argument("--igl", default="")
ap.add_argument("--bao", default="0,1,16,17"); ap.add_argument("--bo_qua_cmd", default=""); ap.add_argument("--du_tru", default=""); ap.add_argument("--ra"); ap.add_argument("--tom"); ap.add_argument("--chu_ky", type=float, default=2.0)
A = ap.parse_args()

def ds(s):
    r = set()
    for p in (s or "").replace(" ", "").split(","):
        if p: a, _, b = p.partition("-"); r |= set(range(int(a), int(b or a) + 1))
    return r
def anh_em(cs):
    r = set()
    for c in cs:
        try: r |= ds(open(f"/sys/devices/system/cpu/cpu{c}/topology/thread_siblings_list").read().strip())
        except OSError: pass
    return r - set(cs)
def stat_cpu():
    o = {}
    for l in open("/proc/stat"):
        if l.startswith("cpu") and l[3].isdigit():
            p = l.split(); v = [int(x) for x in p[1:9]]; o[int(p[0][3:])] = (v[0] + v[1] + v[2], v[5] + v[6])
    return o
def doc_stat(path):
    try: s = open(path).read()
    except OSError: return None
    i = s.rfind(")"); f = s[i + 2:].split()
    return {"comm": s[s.find("(") + 1:i], "ppid": int(f[1]), "tick": int(f[11]) + int(f[12]), "con": int(f[13]) + int(f[14]), "cpu": int(f[36])}
def tat_ca():
    pp = {}
    for d in os.listdir("/proc"):
        if d.isdigit():
            x = doc_stat(f"/proc/{d}/stat")
            if x: pp[int(d)] = x
    return pp
def cay(root, pp):
    con = {}
    for p, x in pp.items(): con.setdefault(x["ppid"], []).append(p)
    r, st = set(), [root]
    while st:
        p = st.pop()
        if p in r or p not in pp: continue
        r.add(p); st += con.get(p, [])
    return r
UIDC = {}
def uid(p):
    if p not in UIDC:
        try: UIDC[p] = int([l for l in open(f"/proc/{p}/status") if l.startswith("Uid:")][0].split()[1])
        except (OSError, IndexError): UIDC[p] = -1
    return UIDC[p]
def affinity(p):
    try:
        for l in open(f"/proc/{p}/status"):
            if l.startswith("Cpus_allowed_list"): return l.split(":")[1].strip()
    except OSError: return None
import re
HE_THONG_OK = re.compile(r"^(falco|tailscaled|kworker|ksoftirqd|rcu_|migration|cpuhp|idle_inject|kthreadd|kswapd|jbd2|irq/|systemd|containerd|dockerd|calico|multipathd|"
                         r"atopacctd|rsyslogd|polkitd|snapd|journald|udevd|kipmi|watchdog|khugepaged|kcompactd|scsi_|nvidia|nv_queue|irqbalance|cron|sshd|dbus|"
                         r"cadvisor|redis|postgres|runc|kube|bird|confd|node_exporter|tini|gdm|Xorg|thermald|upowerd|accounts-daemon)")
BO_QUA = re.compile(A.bo_qua_cmd) if A.bo_qua_cmd else None
CMDC = {}
def bo_qua(p):
    if BO_QUA is None: return False
    if p not in CMDC:
        try: CMDC[p] = open(f"/proc/{p}/cmdline", "rb").read().replace(b"\0", b" ").decode(errors="replace")
        except OSError: CMDC[p] = ""
    return bool(BO_QUA.search(CMDC[p]))
def luong(vung, bo):
    o = {}
    for d in os.listdir("/proc"):
        if not d.isdigit() or int(d) in bo or bo_qua(int(d)): continue
        try: ts = os.listdir(f"/proc/{d}/task")
        except OSError: continue
        for t in ts:
            x = doc_stat(f"/proc/{d}/task/{t}/stat")
            if x: o[(int(d), int(t))] = (x["comm"], x["tick"], x["cpu"])
    return o
def la_theo_chu(l0, l1, vung):
    r = {"nguoi_dung": 0, "he_thong": 0, "khac": 0}; top = {}
    for k, (comm, tk, cpu) in l1.items():
        d = tk - l0.get(k, (comm, tk, cpu))[1] if k in l0 else 0
        if d <= 0 or cpu not in vung: continue
        c = "nguoi_dung" if uid(k[0]) == UID else ("he_thong" if HE_THONG_OK.match(comm) else "khac"); r[c] += d; top[f"{c}:{comm}[{k[0]}]@{cpu}"] = top.get(f"{c}:{comm}[{k[0]}]@{cpu}", 0) + d
    return r, top
def gpu():
    q = subprocess.run(["nvidia-smi", "--query-gpu=index,uuid,utilization.gpu,memory.used,clocks.sm,temperature.gpu,power.draw", "--format=csv,noheader,nounits"],
                       capture_output=True, text=True).stdout
    g, u2i = {}, {}
    for l in q.strip().splitlines():
        p = [x.strip() for x in l.split(",")]; u2i[p[1]] = p[0]
        def f_(x):
            try: return float(x)
            except ValueError: return None
        g[p[0]] = {"util": f_(p[2]), "mem": f_(p[3]), "sm": f_(p[4]), "nhiet": f_(p[5]), "w": f_(p[6]), "apps": []}
    q = subprocess.run(["nvidia-smi", "--query-compute-apps=gpu_uuid,pid,used_memory", "--format=csv,noheader,nounits"], capture_output=True, text=True).stdout
    for l in q.strip().splitlines():
        p = [x.strip() for x in l.split(",")]
        if len(p) == 3 and p[0] in u2i:
            try: m = float(p[2])
            except ValueError: m = 0.0
            g[u2i[p[0]]]["apps"].append((int(p[1]), m))
    return g

if A.yen is not None:
    D = ds(A.yen); t_bd = time.time(); lien = 0; nen = {c: [] for c in D}; me = os.getpid()
    while True:
        bo = cay(me, tat_ca()); a = stat_cpu(); l0 = luong(D, bo); time.sleep(A.giay); b = stat_cpu(); l1 = luong(D, bo)
        ban = {c: (b[c][0] - a[c][0]) / (A.giay * HZ) for c in D if c in a}
        chu, top = la_theo_chu(l0, l1, D)
        gg = {} if A.gpu == "-" else gpu(); apps = [a for x in A.gpu.split(",") if x != "-" for a in gg.get(x, {}).get("apps", [])]
        ok = not apps and chu["nguoi_dung"] / HZ / A.giay <= 0.02 and chu["khac"] / HZ / A.giay <= 0.05
        lien = lien + 1 if ok else 0
        if ok:
            for c, v in ban.items(): nen[c].append(v)
        else: nen = {c: [] for c in D}
        print(time.strftime("%T"), "YEN" if ok else "BAN", f"lien {lien}/{A.cua_so}", "nguoi_dung_loi", round(chu["nguoi_dung"] / HZ / A.giay, 3), "khac_loi", round(chu["khac"] / HZ / A.giay, 3),
              "he_thong_loi", round(chu["he_thong"] / HZ / A.giay, 3), "max_ban", round(max(ban.values() or [0]), 3), "gpu_apps", apps,
              "top", sorted(((round(v / HZ, 2), k) for k, v in top.items()), reverse=True)[:6], flush=True)
        if lien >= A.cua_so:
            if A.nen:
                q = lambda v, p: sorted(v)[min(len(v) - 1, int(p * len(v)))] if v else None
                json.dump({"giay": A.giay, "cua_so": A.cua_so, "loi": {str(c): {"p50": q(v, 0.5), "p95": q(v, 0.95)} for c, v in nen.items()}}, open(A.nen, "w"), indent=1)
            sys.exit(0)
        if time.time() - t_bd > A.cho_toi_da: sys.exit(3)

L, I, BAO, D = ds(A.loi), ds(A.igl), ds(A.bao), ds(A.du_tru); SIB = anh_em(L | I) - BAO; KHAC = D - L - I - SIB
NEN = {}
try: NEN = {int(c): v["p95"] or 0 for c, v in json.load(open(A.nen))["loi"].items()} if A.nen else {}
except (OSError, ValueError, KeyError): NEN = {}
nen_tb = lambda S: (sum(NEN.get(c, 0.0) for c in S) / len(S)) if S else 0.0
def lop(p, aff):
    a = ds(aff or "")
    if a and a <= L: return "L"
    if a and I and a <= I: return "I"
    if a and a <= BAO: return "B"
    return "?"
f = open(A.ra, "w"); t0 = time.time(); me = os.getpid()
pp = tat_ca(); ta = cay(A.pid, pp); aff = {p: affinity(p) for p in ta}; cls = {p: lop(p, aff[p]) for p in ta}
c_prev = {p: pp[p]["tick"] + pp[p]["con"] for p in ta}; cpu0 = stat_cpu(); PARENT = {p: pp[p]["ppid"] for p in ta}
tong = {"ban_L": 0, "ban_I": 0, "ban_SIB": 0, "ban_KHAC": 0, "irq_LI": 0, "ta_L": 0, "ta_I": 0, "ta_B": 0, "ta_?": 0}
chu_tong = {"nguoi_dung": 0, "he_thong": 0, "khac": 0}; top_tong = {}
gpu_la_max = 0; gpu_la_mib = 0.0; gpu_ta_dinh = 0.0; gpu_tong_dinh = 0.0; gpu_kia_max = 0; util, sm, nhiet, w, tan = [], [], [], [], []
n = 0; l0 = None; aff_ghi = {}
def con_song(p):
    try: os.kill(p, 0); return True
    except OSError: return False
while con_song(A.pid):
    time.sleep(A.chu_ky); n += 1
    pp = tat_ca(); ta1 = cay(A.pid, pp)
    for p in ta1:
        if p not in aff: aff[p] = affinity(p); cls[p] = lop(p, aff[p])
    cpu1 = stat_cpu(); now = {}
    for p in ta1:
        x = doc_stat(f"/proc/{p}/stat")
        if x: now[p] = x["tick"] + x["con"]
    dt = {"L": 0, "I": 0, "B": 0, "?": 0}
    for p, v in now.items(): dt[cls[p]] += v - c_prev.get(p, 0)
    for p in set(c_prev) - set(now):
        par = PARENT.get(p)
        if par in now: dt[cls[par]] -= c_prev[p]
    PARENT = {p: pp[p]["ppid"] for p in ta1 if p in pp}
    c_prev = now
    for p in ta1:
        if now.get(p, 0) > 0: aff_ghi[f"{pp[p]['comm']}:{aff.get(p)}"] = cls[p]
    ban = {c: cpu1[c][0] - cpu0[c][0] for c in D}; irq = {c: cpu1[c][1] - cpu0[c][1] for c in D}; cpu0 = cpu1
    s = lambda S: sum(ban.get(c, 0) for c in S)
    tong["ban_L"] += s(L); tong["ban_I"] += s(I); tong["ban_SIB"] += s(SIB); tong["ban_KHAC"] += s(KHAC); tong["irq_LI"] += sum(irq.get(c, 0) for c in L | I)
    for k in dt: tong[f"ta_{k}"] += dt[k]
    mau = {"t": round(time.time() - t0, 1), "ta": len(ta1), "dt_ta": dt, "ban_L": s(L), "ban_I": s(I), "ban_SIB": s(SIB)}
    if n % max(1, int(round(10 / A.chu_ky))) == 1 or l0 is None:
        l1 = luong(D, ta1 | {me})
        if l0 is not None:
            chu, top = la_theo_chu(l0, l1, D)
            for k in chu: chu_tong[k] += chu[k]
            for k, v in top.items(): top_tong[k] = top_tong.get(k, 0) + v
        l0 = l1
    if A.gpu != "-":
        g = gpu(); meg = g.get(A.gpu, {"apps": []})
        la = [(p, m) for p, m in meg["apps"] if p not in ta1]; ta_m = sum(m for p, m in meg["apps"] if p in ta1)
        gpu_la_max = max(gpu_la_max, len(la)); gpu_la_mib = max(gpu_la_mib, sum(m for _, m in la)); gpu_ta_dinh = max(gpu_ta_dinh, ta_m)
        gpu_tong_dinh = max(gpu_tong_dinh, meg.get("mem") or 0); gpu_kia_max = max([gpu_kia_max] + [len(v["apps"]) for k, v in g.items() if k != A.gpu])
        if meg.get("util") is not None: util.append(meg["util"]); sm.append(meg["sm"] or 0); nhiet.append(meg["nhiet"] or 0); w.append(meg["w"] or 0)
        mau.update({"gpu_la": la, "gpu_ta_mib": ta_m, "gpu_util": meg.get("util"), "sm": meg.get("sm"), "nhiet": meg.get("nhiet")})
    try: tan.append(sum(int(open(f"/sys/devices/system/cpu/cpu{c}/cpufreq/scaling_cur_freq").read()) for c in L) / max(len(L), 1) / 1000)
    except (OSError, ValueError): pass
    f.write(json.dumps(mau) + "\n"); f.flush()
wall = max(time.time() - t0, 1e-9); cs = lambda x: round(x / HZ, 2)
la = {"L": max(0, tong["ban_L"] - tong["ta_L"]), "I": max(0, tong["ban_I"] - tong["ta_I"]), "SIB": tong["ban_SIB"], "KHAC": tong["ban_KHAC"]}
fr = lambda tk, S: round(tk / HZ / (wall * len(S)), 5) if S else None
ti = {"la_L": fr(la["L"], L), "la_I": fr(la["I"], I), "la_SIB": fr(la["SIB"], SIB), "la_KHAC": fr(la["KHAC"], KHAC)}
nguong = {"la_L": nen_tb(L) + 0.02, "la_I": nen_tb(I) + 0.05, "la_SIB": nen_tb(SIB) + 0.05}
khong_ghim = sorted(k for k, v in aff_ghi.items() if v == "?")
tom = {"wall_s": round(wall, 1), "so_mau": n, "loi": A.loi, "igl": A.igl, "anh_em": sorted(SIB), "khac": sorted(KHAC),
       "cpu_s": {**{k: cs(v) for k, v in tong.items()}, **{f"la_{k}": cs(v) for k, v in la.items()}}, "ti_le": ti, "nguong": nguong, "nen_p95_tb": {"L": nen_tb(L), "I": nen_tb(I), "SIB": nen_tb(SIB)},
       "la_theo_chu_loi": {k: round(v / HZ / wall, 4) for k, v in chu_tong.items()}, "ke_la_top": sorted(((cs(v), k) for k, v in top_tong.items()), reverse=True)[:12],
       "gpu": A.gpu, "gpu_la_max": gpu_la_max, "gpu_la_mib_max": gpu_la_mib, "gpu_ta_dinh_mib": gpu_ta_dinh, "gpu_tong_dinh_mib": gpu_tong_dinh,
       "gpu_kia_so_tien_trinh_max": gpu_kia_max, "gpu_util_tb": round(sum(util) / len(util), 1) if util else None, "sm_mhz_tb": round(sum(sm) / len(sm)) if sm else None,
       "nhiet_max": max(nhiet) if nhiet else None, "w_tb": round(sum(w) / len(w), 1) if w else None, "tan_so_L_mhz_tb": round(sum(tan) / len(tan)) if tan else None,
       "affinity": aff_ghi, "khong_ghim": khong_ghim}
dk = {"gpu_khong_la": A.gpu == "-" or gpu_la_max == 0, "nguoi_dung<=0.01_loi": tom["la_theo_chu_loi"]["nguoi_dung"] <= 0.01, "uid_khac_ngoai_ds<=0.03_loi": tom["la_theo_chu_loi"]["khac"] <= 0.03,
      "la_L<=nen+0.02": (ti["la_L"] or 0) <= nguong["la_L"], "ghim": not khong_ghim}
tom["vuot_nen"] = {k: (ti[k] or 0) > nguong[k] for k in ("la_I", "la_SIB")}
tom["dieu_kien"] = dk; tom["ket_luan"] = "DOC_QUYEN" if all(dk.values()) else "CHONG"
json.dump(tom, open(A.tom, "w"), ensure_ascii=False, indent=1); print(tom["ket_luan"], json.dumps(ti))
